#!/usr/bin/env python3
"""Build frozen unbiased insect and plant discovery rankings.

This script consumes only the validated links in 01_input_links. It does not
rerun DESeq2, SignalP, DeepTMHMM, eggNOG-mapper or MEROPS, and it does not use
literature-derived judgments to alter ranking membership.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path


NA_VALUES = {"", "na", "nan", "none", "null", ".", "-"}
GENERIC_DESCRIPTION_PHRASES = (
    "uncharacterized protein",
    "hypothetical protein",
    "protein of unknown function",
    "unknown protein",
    "function unknown",
    "unknown function",
    "uncharacterised protein",
    "predicted protein",
    "unnamed protein product",
)


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip().lstrip("\ufeff")


def useful(value: object) -> bool:
    return clean(value).lower() not in NA_VALUES


def as_float(value: object) -> float | None:
    text = clean(value)
    if text.lower() in NA_VALUES:
        return None
    number = float(text)
    return None if math.isnan(number) else number


def as_bool(value: object) -> bool:
    text = clean(value).lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n", ""}:
        return False
    raise ValueError(f"Unrecognized Boolean value: {value!r}")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        fields = [clean(field) for field in reader.fieldnames]
        rows = [
            {clean(key): clean(value) for key, value in raw.items()}
            for raw in reader
        ]
    return fields, rows


def unique_index(rows: list[dict[str, str]], key: str, label: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = row.get(key, "")
        if not identifier:
            raise ValueError(f"{label}: empty {key}")
        if identifier in result:
            raise ValueError(f"{label}: duplicate {key}: {identifier}")
        result[identifier] = row
    return result


def join_values(values: object, separator: str = ";") -> str:
    if isinstance(values, str):
        candidates = [values]
    else:
        candidates = list(values)
    cleaned = sorted({clean(value) for value in candidates if useful(value)})
    return separator.join(cleaned)


def aggregate_mapping(path: Path, gene_key: str, term_key: str, eligible: set[str]) -> dict[str, str]:
    terms: dict[str, set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            row = {clean(key): clean(value) for key, value in raw.items()}
            gene = row.get(gene_key, "")
            term = row.get(term_key, "")
            if gene in eligible and useful(term):
                terms[gene].add(term)
    return {gene: join_values(values) for gene, values in terms.items()}


def percentile_ranks(values: list[float]) -> list[Fraction]:
    """Average rank divided by N; larger input values receive larger percentiles."""
    n = len(values)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda index: values[index])
    result = [Fraction(0, 1)] * n
    start = 0
    while start < n:
        end = start + 1
        while end < n and values[order[end]] == values[order[start]]:
            end += 1
        # Exact rational arithmetic prevents binary floating-point noise from
        # changing the declared downstream tie-break order.
        percentile = Fraction((start + 1) + end, 2 * n)
        for position in range(start, end):
            result[order[position]] = percentile
        start = end
    return result


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def specific_description(description: str) -> bool:
    text = clean(description).lower()
    if text in NA_VALUES:
        return False
    return not any(phrase in text for phrase in GENERIC_DESCRIPTION_PHRASES)


def only_duf_pfams(pfams: str) -> bool:
    tokens = [
        token.strip().lower()
        for token in pfams.replace(",", ";").split(";")
        if useful(token)
    ]
    return bool(tokens) and all(token.startswith("duf") for token in tokens)


def annotation_status(
    *,
    description: str,
    preferred_name: str,
    pfams: str,
    cazys: str,
    gos: str,
    kegg_kos: str,
    kegg_pathways: str,
    seed_ortholog: str,
    prior_supported_family: bool,
    merops_supported: bool,
) -> tuple[str, str]:
    if prior_supported_family:
        return "informative", "prior domain/annotation-supported candidate-family evidence"
    if merops_supported:
        return "informative", "passing MEROPS homology evidence"
    if useful(preferred_name):
        return "informative", "preferred gene/product name available"
    if specific_description(description):
        return "informative", "specific functional description available"
    if useful(cazys):
        return "informative", "CAZy family assignment available"
    if useful(pfams) and not only_duf_pfams(pfams):
        return "informative", "non-DUF Pfam-domain evidence available"

    clues: list[str] = []
    if useful(description):
        clues.append("generic description")
    if useful(pfams):
        clues.append("DUF/Pfam clue")
    if useful(gos):
        clues.append("GO assignment")
    if useful(kegg_kos) or useful(kegg_pathways):
        clues.append("KEGG assignment")
    if useful(seed_ortholog):
        clues.append("orthology assignment")
    if clues:
        return "poorly_characterized", "; ".join(clues)
    return "uncharacterized", "no meaningful annotation evidence in integrated inputs"


def format_float(value: float | None, digits: int = 12) -> str:
    if value is None:
        return ""
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    return f"{value:.{digits}g}"


def ranking_sort_key(row: dict[str, object]) -> tuple[object, ...]:
    return (
        -row["_discovery_score_exact"],
        float(row["_padj"]),
        -abs(float(row["_lfc"])),
        -float(row["_baseMean"]),
        str(row["gene_id"]),
    )


def add_ranking(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    significance: list[float] = []
    effect: list[float] = []
    expression: list[float] = []
    for row in rows:
        padj = float(row["_padj"])
        lfc = float(row["_lfc"])
        base_mean = float(row["_baseMean"])
        significance.append(math.inf if padj == 0 else -math.log10(padj))
        effect.append(abs(lfc))
        expression.append(math.log10(base_mean + 1.0))

    p_significance = percentile_ranks(significance)
    p_effect = percentile_ranks(effect)
    p_expression = percentile_ranks(expression)

    for index, row in enumerate(rows):
        score = (p_significance[index] + p_effect[index] + p_expression[index]) / 3
        row["_significance_percentile_exact"] = p_significance[index]
        row["_effect_size_percentile_exact"] = p_effect[index]
        row["_expression_percentile_exact"] = p_expression[index]
        row["_discovery_score_exact"] = score

    rows.sort(key=ranking_sort_key)
    for rank, row in enumerate(rows, start=1):
        row["discovery_rank"] = rank
        row["significance_percentile"] = (
            f"{float(row['_significance_percentile_exact']):.12f}"
        )
        row["effect_size_percentile"] = (
            f"{float(row['_effect_size_percentile_exact']):.12f}"
        )
        row["expression_percentile"] = (
            f"{float(row['_expression_percentile_exact']):.12f}"
        )
        row["discovery_score"] = f"{float(row['_discovery_score_exact']):.12f}"
    return rows


def group_rows(rows: list[dict[str, str]], key: str) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row[key]].append(row)
    return grouped


def build_insect(input_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    _, integrated_rows = read_tsv(input_dir / "insect_signalp_tmhmm_deseq2_integrated.tsv")
    eligible_source = [
        row
        for row in integrated_rows
        if row["analysis_class"] == "classical_secreted_core"
        and row["deseq2_status"] == "larva_up"
    ]
    if len(eligible_source) != 544:
        raise ValueError(f"Expected 544 insect eligible genes; observed {len(eligible_source)}")
    eligible_ids = {row["gene_id"] for row in eligible_source}
    if len(eligible_ids) != 544:
        raise ValueError("Duplicate insect eligible gene IDs")

    _, egg_rows = read_tsv(input_dir / "insect_eggnog_parsed.tsv")
    egg = unique_index(egg_rows, "Geneid", "insect eggNOG")
    gos = aggregate_mapping(input_dir / "insect_gene_to_GO.tsv", "Geneid", "GO", eligible_ids)
    kos = aggregate_mapping(
        input_dir / "insect_gene_to_KEGG_ko.tsv", "Geneid", "KEGG_ko", eligible_ids
    )
    pathways = aggregate_mapping(
        input_dir / "insect_gene_to_KEGG_pathway.tsv",
        "Geneid",
        "KEGG_Pathway",
        eligible_ids,
    )

    _, focused_rows = read_tsv(input_dir / "insect_focused_candidate_families.tsv")
    focused = group_rows(focused_rows, "gene_id")
    _, cap_rows = read_tsv(input_dir / "insect_CAP_domain_candidates.tsv")
    cap = group_rows(cap_rows, "gene_id")
    _, merops_rows = read_tsv(input_dir / "insect_MEROPS_gene_category_curated_v2.tsv")
    merops = group_rows(merops_rows, "gene_id")

    output: list[dict[str, object]] = []
    for source in eligible_source:
        gene = source["gene_id"]
        ann = egg.get(gene, {})
        focused_groups = join_values(row["candidate_group"] for row in focused.get(gene, []))
        cap_member = bool(cap.get(gene))
        merops_gene_rows = merops.get(gene, [])
        merops_categories = join_values(row["merops_category"] for row in merops_gene_rows)
        merops_families = join_values(row["merops_family"] for row in merops_gene_rows)
        merops_identifiers = join_values(row["merops_identifier"] for row in merops_gene_rows)
        merops_tiers = join_values(row["evidence_tier"] for row in merops_gene_rows)
        has_inhibitor = any(
            "inhibitor" in row["merops_category"].lower() for row in merops_gene_rows
        )
        has_peptidase = any(
            "peptidase" in row["merops_category"].lower() for row in merops_gene_rows
        )
        prior_family = bool(focused_groups) or cap_member
        merops_supported = bool(merops_gene_rows)

        status, basis = annotation_status(
            description=ann.get("Description", ""),
            preferred_name=ann.get("Preferred_name", ""),
            pfams=ann.get("PFAMs", ""),
            cazys=ann.get("CAZy", ""),
            gos=gos.get(gene, ann.get("GOs", "")),
            kegg_kos=kos.get(gene, ann.get("KEGG_ko", "")),
            kegg_pathways=pathways.get(gene, ann.get("KEGG_Pathway", "")),
            seed_ortholog=ann.get("seed_ortholog", ""),
            prior_supported_family=prior_family,
            merops_supported=merops_supported,
        )
        if prior_family:
            classification = "known_focused_family"
        elif has_inhibitor:
            classification = "MEROPS_inhibitor_homology"
        elif has_peptidase:
            classification = "MEROPS_peptidase_homology"
        elif status == "informative":
            classification = "other_annotated_protein"
        elif status == "poorly_characterized":
            classification = "poorly_characterized"
        else:
            classification = "uncharacterized"

        prior_sources: list[str] = []
        if focused_groups:
            prior_sources.append("focused_candidate_family")
        if cap_member:
            prior_sources.append("CAP_domain_audit")
        if merops_supported:
            prior_sources.append("MEROPS_curated_v2")

        padj = as_float(source["padj"])
        lfc = as_float(source["log2FoldChange"])
        base_mean = as_float(source["baseMean"])
        if padj is None or lfc is None or base_mean is None:
            raise ValueError(f"Missing ranking metric for insect gene {gene}")
        if not (padj < 0.01 and lfc <= -1):
            raise ValueError(f"Ineligible insect gene reached ranking: {gene}")

        output.append(
            {
                "gene_id": gene,
                "baseMean": source["baseMean"],
                "log2FoldChange_pupa_vs_larva": source["log2FoldChange"],
                "absolute_log2FoldChange": format_float(abs(lfc)),
                "stat": source["stat"],
                "pvalue": source["pvalue"],
                "padj": source["padj"],
                "DESeq2_status": source["deseq2_status"],
                "representative_protein_length_aa": source["length_aa"],
                "signalp_prediction": source["signalp_prediction"],
                "signalp_sp_probability": source["signalp_sp_probability"],
                "signalp_cs_start": source["signalp_cs_start"],
                "signalp_cs_end": source["signalp_cs_end"],
                "signalp_cs_probability": source["signalp_cs_probability"],
                "deeptmhmm_type": source["deeptmhmm_type"],
                "deeptmhmm_predicted_tmrs": source["deeptmhmm_predicted_tmrs"],
                "deeptmhmm_tmr_regions": source["deeptmhmm_tmr_regions"],
                "analysis_class": source["analysis_class"],
                "eggNOG_description": ann.get("Description", ""),
                "preferred_name": ann.get("Preferred_name", ""),
                "seed_ortholog": ann.get("seed_ortholog", ""),
                "eggNOG_OGs": ann.get("eggNOG_OGs", ""),
                "PFAMs": ann.get("PFAMs", ""),
                "CAZy": ann.get("CAZy", ""),
                "GOs": gos.get(gene, ann.get("GOs", "")),
                "KEGG_ko": kos.get(gene, ann.get("KEGG_ko", "")),
                "KEGG_Pathway": pathways.get(gene, ann.get("KEGG_Pathway", "")),
                "focused_candidate_groups": focused_groups,
                "CAP_candidate": "TRUE" if cap_member else "FALSE",
                "MEROPS_screen_scope": "complete_representative_proteome_screened",
                "MEROPS_result": merops_categories or "no_passing_MEROPS_hit",
                "MEROPS_families": merops_families,
                "MEROPS_identifiers": merops_identifiers,
                "MEROPS_evidence_tiers": merops_tiers,
                "annotation_status": status,
                "annotation_status_basis": basis,
                "discovery_classification": classification,
                "previously_in_focused_analysis": "TRUE" if prior_sources else "FALSE",
                "prior_analysis_sources": join_values(prior_sources),
                "_padj": padj,
                "_lfc": lfc,
                "_baseMean": base_mean,
            }
        )

    ranked = add_ranking(output)
    top50 = [dict(row) for row in ranked[:50]]
    poor_pool = [
        row
        for row in ranked
        if row["annotation_status"] in {"poorly_characterized", "uncharacterized"}
    ]
    if len(poor_pool) < 15:
        raise ValueError(
            f"Only {len(poor_pool)} poorly characterized/uncharacterized insect genes; expected >=15"
        )
    top15 = [dict(row) for row in poor_pool[:15]]
    top50_ids = {str(row["gene_id"]) for row in top50}
    for row in top15:
        row["in_overall_top50"] = "TRUE" if str(row["gene_id"]) in top50_ids else "FALSE"
    return ranked, top50, top15


def build_plant(input_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    _, de_rows = read_tsv(input_dir / "plant_deseq2_all.tsv")
    eligible_source = [row for row in de_rows if row["DESeq2_status"] == "GP_up"]
    if len(eligible_source) != 5593:
        raise ValueError(f"Expected 5593 plant GP-up genes; observed {len(eligible_source)}")
    eligible_ids = {row["gene_id"] for row in eligible_source}
    if len(eligible_ids) != 5593:
        raise ValueError("Duplicate plant eligible gene IDs")

    _, joined_rows = read_tsv(input_dir / "plant_deseq2_eggnog_joined.tsv")
    joined = unique_index(joined_rows, "gene_id", "plant DESeq2/eggNOG join")
    _, integrated_rows = read_tsv(input_dir / "plant_signalp_tmhmm_deseq2_integrated.tsv")
    integrated = unique_index(integrated_rows, "gene_id", "plant secretome integration")
    _, protein_map_rows = read_tsv(input_dir / "plant_representative_protein_map.tsv")
    protein_map = unique_index(
        protein_map_rows, "gene_id", "plant representative-protein map"
    )
    if set(protein_map) != set(integrated):
        raise ValueError("Plant representative-protein map and integration gene sets differ")
    _, merops_rows = read_tsv(input_dir / "plant_MEROPS_gene_category_curated_v2.tsv")
    merops = group_rows(merops_rows, "gene_id")
    _, priority_rows = read_tsv(input_dir / "plant_MEROPS_peptidase_priority_curated_v2.tsv")
    priority = unique_index(priority_rows, "gene_id", "plant MEROPS priority")
    _, counterpart_rows = read_tsv(input_dir / "plant_counterpart_family_audit.tsv")
    counterpart = group_rows(counterpart_rows, "gene_id")
    _, development_rows = read_tsv(input_dir / "plant_developmental_candidates.tsv")
    development = group_rows(development_rows, "gene_id")

    output: list[dict[str, object]] = []
    for source in eligible_source:
        gene = source["gene_id"]
        ann = joined.get(gene, {})
        protein = integrated.get(gene, {})
        protein_record = protein_map.get(gene, {})
        merops_gene_rows = merops.get(gene, [])
        priority_row = priority.get(gene, {})
        counterpart_groups = join_values(
            row["candidate_group"] for row in counterpart.get(gene, [])
        )
        developmental_groups = join_values(
            row["target_groups"] for row in development.get(gene, [])
        )
        developmental_regulators = join_values(
            row["target_regulators"] for row in development.get(gene, [])
        )
        merops_categories = join_values(row["merops_category"] for row in merops_gene_rows)
        merops_families = join_values(row["merops_family"] for row in merops_gene_rows)
        merops_identifiers = join_values(row["merops_identifier"] for row in merops_gene_rows)
        merops_tiers = join_values(row["evidence_tier"] for row in merops_gene_rows)
        prior_family = bool(counterpart_groups or developmental_groups)
        merops_supported = bool(merops_gene_rows)

        status, basis = annotation_status(
            description=ann.get("Description", ""),
            preferred_name=ann.get("Preferred_name", ""),
            pfams=ann.get("PFAMs", ""),
            cazys=ann.get("CAZy", ""),
            gos=ann.get("GOs", ""),
            kegg_kos=ann.get("KEGG_ko", ""),
            kegg_pathways=ann.get("KEGG_Pathway", ""),
            seed_ortholog=ann.get("seed_ortholog", ""),
            prior_supported_family=prior_family,
            merops_supported=merops_supported,
        )
        has_inhibitor = any(
            "inhibitor" in row["merops_category"].lower() for row in merops_gene_rows
        )
        has_peptidase = any(
            "peptidase" in row["merops_category"].lower() for row in merops_gene_rows
        )
        if prior_family:
            classification = "known_targeted_family_or_regulator"
        elif has_peptidase:
            classification = "MEROPS_peptidase_homology"
        elif has_inhibitor:
            classification = "MEROPS_inhibitor_homology"
        elif status == "informative":
            classification = "other_annotated_gene"
        elif status == "poorly_characterized":
            classification = "poorly_characterized"
        else:
            classification = "uncharacterized"

        prior_sources: list[str] = []
        if counterpart_groups:
            prior_sources.append("counterpart_family_audit")
        if developmental_groups:
            prior_sources.append("developmental_candidate_audit")
        if merops_supported:
            prior_sources.append("MEROPS_curated_v2")

        has_protein = bool(protein)
        is_screened_scope = (
            has_protein
            and protein.get("analysis_class") == "classical_secreted_core"
            and protein.get("DESeq2_status") == "GP_up"
        )
        if is_screened_scope:
            merops_scope = "screened_GP_up_classical_secretome"
            merops_result = merops_categories or "no_passing_MEROPS_hit"
        else:
            merops_scope = "not_screened_outside_GP_up_classical_secretome"
            merops_result = "not_assessed"

        padj = as_float(source["padj"])
        lfc = as_float(source["log2FoldChange"])
        base_mean = as_float(source["baseMean"])
        if padj is None or lfc is None or base_mean is None:
            raise ValueError(f"Missing ranking metric for plant gene {gene}")
        if not (padj < 0.01 and lfc >= 1):
            raise ValueError(f"Ineligible plant gene reached ranking: {gene}")

        output.append(
            {
                "gene_id": gene,
                "transcript_id": protein_record.get(
                    "transcript_id", ann.get("transcript_id", "")
                ),
                "representative_protein_id": protein_record.get(
                    "protein_id", ann.get("protein_id", "")
                ),
                "representative_protein_length_aa": protein_record.get(
                    "aa_length", protein.get("protein_length", "")
                ),
                "representative_protein_available": "TRUE" if has_protein else "FALSE",
                "baseMean": source["baseMean"],
                "log2FoldChange_GP_vs_UP": source["log2FoldChange"],
                "absolute_log2FoldChange": format_float(abs(lfc)),
                "stat": source["stat"],
                "pvalue": source["pvalue"],
                "padj": source["padj"],
                "DESeq2_status": source["DESeq2_status"],
                "eggNOG_annotation_available": ann.get("has_eggnog_annotation", ""),
                "eggNOG_description": ann.get("Description", ""),
                "preferred_name": ann.get("Preferred_name", ""),
                "seed_ortholog": ann.get("seed_ortholog", ""),
                "eggNOG_OGs": ann.get("eggNOG_OGs", ""),
                "PFAMs": ann.get("PFAMs", ""),
                "CAZy": ann.get("CAZy", ""),
                "GOs": ann.get("GOs", ""),
                "KEGG_ko": ann.get("KEGG_ko", ""),
                "KEGG_Pathway": ann.get("KEGG_Pathway", ""),
                "signalp_prediction": protein.get("signalp_prediction", ""),
                "signalp_sp_probability": protein.get("signalp_SP_score", ""),
                "signalp_cleavage_after": protein.get("signalp_cleavage_after", ""),
                "signalp_cleavage_probability": protein.get(
                    "signalp_cleavage_probability", ""
                ),
                "deeptmhmm_type": protein.get("deeptmhmm_type", ""),
                "deeptmhmm_predicted_tmrs": protein.get("predicted_TMRs", ""),
                "deeptmhmm_tmr_regions": protein.get("TMR_segments", ""),
                "analysis_class": protein.get("analysis_class", ""),
                "counterpart_candidate_groups": counterpart_groups,
                "developmental_candidate_groups": developmental_groups,
                "developmental_target_regulators": developmental_regulators,
                "MEROPS_screen_scope": merops_scope,
                "MEROPS_result": merops_result,
                "MEROPS_families": merops_families,
                "MEROPS_identifiers": merops_identifiers,
                "MEROPS_evidence_tiers": merops_tiers,
                "MEROPS_manuscript_priority": priority_row.get(
                    "manuscript_priority", "not_in_peptidase_candidate_table"
                ),
                "annotation_status": status,
                "annotation_status_basis": basis,
                "discovery_classification": classification,
                "previously_in_focused_analysis": "TRUE" if prior_sources else "FALSE",
                "prior_analysis_sources": join_values(prior_sources),
                "_padj": padj,
                "_lfc": lfc,
                "_baseMean": base_mean,
            }
        )

    ranked = add_ranking(output)
    top50 = [dict(row) for row in ranked[:50]]
    return ranked, top50


INSECT_FIELDS = [
    "discovery_rank",
    "discovery_score",
    "significance_percentile",
    "effect_size_percentile",
    "expression_percentile",
    "gene_id",
    "baseMean",
    "log2FoldChange_pupa_vs_larva",
    "absolute_log2FoldChange",
    "stat",
    "pvalue",
    "padj",
    "DESeq2_status",
    "representative_protein_length_aa",
    "signalp_prediction",
    "signalp_sp_probability",
    "signalp_cs_start",
    "signalp_cs_end",
    "signalp_cs_probability",
    "deeptmhmm_type",
    "deeptmhmm_predicted_tmrs",
    "deeptmhmm_tmr_regions",
    "analysis_class",
    "eggNOG_description",
    "preferred_name",
    "seed_ortholog",
    "eggNOG_OGs",
    "PFAMs",
    "CAZy",
    "GOs",
    "KEGG_ko",
    "KEGG_Pathway",
    "focused_candidate_groups",
    "CAP_candidate",
    "MEROPS_screen_scope",
    "MEROPS_result",
    "MEROPS_families",
    "MEROPS_identifiers",
    "MEROPS_evidence_tiers",
    "annotation_status",
    "annotation_status_basis",
    "discovery_classification",
    "previously_in_focused_analysis",
    "prior_analysis_sources",
]

PLANT_FIELDS = [
    "discovery_rank",
    "discovery_score",
    "significance_percentile",
    "effect_size_percentile",
    "expression_percentile",
    "gene_id",
    "transcript_id",
    "representative_protein_id",
    "representative_protein_length_aa",
    "representative_protein_available",
    "baseMean",
    "log2FoldChange_GP_vs_UP",
    "absolute_log2FoldChange",
    "stat",
    "pvalue",
    "padj",
    "DESeq2_status",
    "eggNOG_annotation_available",
    "eggNOG_description",
    "preferred_name",
    "seed_ortholog",
    "eggNOG_OGs",
    "PFAMs",
    "CAZy",
    "GOs",
    "KEGG_ko",
    "KEGG_Pathway",
    "signalp_prediction",
    "signalp_sp_probability",
    "signalp_cleavage_after",
    "signalp_cleavage_probability",
    "deeptmhmm_type",
    "deeptmhmm_predicted_tmrs",
    "deeptmhmm_tmr_regions",
    "analysis_class",
    "counterpart_candidate_groups",
    "developmental_candidate_groups",
    "developmental_target_regulators",
    "MEROPS_screen_scope",
    "MEROPS_result",
    "MEROPS_families",
    "MEROPS_identifiers",
    "MEROPS_evidence_tiers",
    "MEROPS_manuscript_priority",
    "annotation_status",
    "annotation_status_basis",
    "discovery_classification",
    "previously_in_focused_analysis",
    "prior_analysis_sources",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Discovery-project root (default: parent of 07_scripts)",
    )
    args = parser.parse_args()

    root = args.project_root.expanduser().resolve()
    input_dir = root / "01_input_links"
    insect_dir = root / "02_insect_discovery" / "results"
    plant_dir = root / "03_plant_discovery" / "results"
    control_dir = root / "00_project_control"
    log_dir = root / "08_logs"
    for directory in (input_dir, insect_dir, plant_dir, control_dir, log_dir):
        if not directory.is_dir():
            raise SystemExit(f"Missing required directory: {directory}")

    insect_full, insect_top50, insect_poor15 = build_insect(input_dir)
    plant_full, plant_top50 = build_plant(input_dir)

    insect_full_path = insect_dir / "insect_larva_secretome_discovery_full.tsv"
    insect_top50_path = insect_dir / "insect_larva_secretome_top50.tsv"
    insect_poor15_path = (
        insect_dir / "insect_larva_secretome_poorly_annotated_top15.tsv"
    )
    plant_full_path = plant_dir / "plant_GP_discovery_full.tsv"
    plant_top50_path = plant_dir / "plant_GP_top50.tsv"

    write_tsv(insect_full_path, INSECT_FIELDS, insect_full)
    write_tsv(insect_top50_path, INSECT_FIELDS, insect_top50)
    write_tsv(
        insect_poor15_path,
        INSECT_FIELDS + ["in_overall_top50"],
        insect_poor15,
    )
    write_tsv(plant_full_path, PLANT_FIELDS, plant_full)
    write_tsv(plant_top50_path, PLANT_FIELDS, plant_top50)

    parameters = [
        {
            "parameter": "official_DEG_rule",
            "value": "padj < 0.01 AND absolute log2FoldChange >= 1",
            "purpose": "eligibility only",
        },
        {
            "parameter": "insect_eligible_direction",
            "value": "log2FoldChange <= -1 in pupa_vs_larva contrast",
            "purpose": "larva-up definition",
        },
        {
            "parameter": "insect_secretome_rule",
            "value": "analysis_class == classical_secreted_core",
            "purpose": "use final audited SignalP/DeepTMHMM classification",
        },
        {
            "parameter": "plant_eligible_direction",
            "value": "log2FoldChange >= 1 in GP_vs_UP contrast",
            "purpose": "GP-up definition",
        },
        {
            "parameter": "plant_secretome_requirement",
            "value": "none",
            "purpose": "retain intracellular and extracellular GP-up genes",
        },
        {
            "parameter": "ranking_component_1",
            "value": "percentile rank of -log10(padj)",
            "purpose": "significance",
        },
        {
            "parameter": "ranking_component_2",
            "value": "percentile rank of absolute log2FoldChange",
            "purpose": "effect size",
        },
        {
            "parameter": "ranking_component_3",
            "value": "percentile rank of log10(baseMean + 1)",
            "purpose": "expression",
        },
        {
            "parameter": "percentile_method",
            "value": "average rank divided by cohort size; larger value ranks higher",
            "purpose": "ties receive their shared average rank",
        },
        {
            "parameter": "discovery_score",
            "value": "arithmetic mean of the three percentile ranks",
            "purpose": "equal-weight prioritization only",
        },
        {
            "parameter": "deterministic_sort",
            "value": "score desc; padj asc; absolute LFC desc; baseMean desc; gene ID asc",
            "purpose": "stable top-list membership",
        },
        {
            "parameter": "annotation_status",
            "value": "informative / poorly_characterized / uncharacterized",
            "purpose": "frozen evidence-based annotation triage",
        },
        {
            "parameter": "plant_MEROPS_scope",
            "value": "549 GP-up classical-secretome genes only",
            "purpose": "unscreened genes are not classified as MEROPS-negative",
        },
        {
            "parameter": "literature_search_before_freeze",
            "value": "none",
            "purpose": "prevent post hoc alteration of ranked membership",
        },
    ]
    parameter_path = control_dir / "discovery_ranking_parameters.tsv"
    write_tsv(parameter_path, ["parameter", "value", "purpose"], parameters)

    frozen_files = [insect_top50_path, insect_poor15_path, plant_top50_path]
    checksum_rows = [
        {
            "file": str(path.relative_to(root)),
            "rows": sum(1 for _ in path.open("r", encoding="utf-8")) - 1,
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in frozen_files
    ]
    checksum_path = control_dir / "frozen_top_list_SHA256.tsv"
    write_tsv(checksum_path, ["file", "rows", "size_bytes", "sha256"], checksum_rows)

    errors: list[str] = []
    if len(insect_full) != 544:
        errors.append(f"Expected 544 insect full rows; observed {len(insect_full)}")
    if len(insect_top50) != 50:
        errors.append(f"Expected 50 insect top rows; observed {len(insect_top50)}")
    if len(insect_poor15) != 15:
        errors.append(f"Expected 15 insect poorly annotated rows; observed {len(insect_poor15)}")
    if len(plant_full) != 5593:
        errors.append(f"Expected 5593 plant full rows; observed {len(plant_full)}")
    if len(plant_top50) != 50:
        errors.append(f"Expected 50 plant top rows; observed {len(plant_top50)}")
    if len({row["gene_id"] for row in insect_full}) != len(insect_full):
        errors.append("Duplicate gene IDs in insect full ranking")
    if len({row["gene_id"] for row in plant_full}) != len(plant_full):
        errors.append("Duplicate gene IDs in plant full ranking")
    if any(int(row["discovery_rank"]) != index for index, row in enumerate(insect_full, 1)):
        errors.append("Nonsequential insect discovery ranks")
    if any(int(row["discovery_rank"]) != index for index, row in enumerate(plant_full, 1)):
        errors.append("Nonsequential plant discovery ranks")
    if insect_full != sorted(insect_full, key=ranking_sort_key):
        errors.append("Insect rows violate exact deterministic sort order")
    if plant_full != sorted(plant_full, key=ranking_sort_key):
        errors.append("Plant rows violate exact deterministic sort order")
    if [row["gene_id"] for row in insect_top50] != [
        row["gene_id"] for row in insect_full[:50]
    ]:
        errors.append("Insect top-50 is not the first 50 rows of the full ranking")
    if [row["gene_id"] for row in plant_top50] != [
        row["gene_id"] for row in plant_full[:50]
    ]:
        errors.append("Plant top-50 is not the first 50 rows of the full ranking")
    expected_insect_poor15 = [
        row["gene_id"]
        for row in insect_full
        if row["annotation_status"] in {"poorly_characterized", "uncharacterized"}
    ][:15]
    if [row["gene_id"] for row in insect_poor15] != expected_insect_poor15:
        errors.append(
            "Insect poorly annotated top-15 is not the first 15 eligible rows "
            "of the full ranking"
        )
    for species, rows in (("Insect", insect_full), ("Plant", plant_full)):
        for row in rows:
            expected_score = (
                row["_significance_percentile_exact"]
                + row["_effect_size_percentile_exact"]
                + row["_expression_percentile_exact"]
            ) / 3
            if row["_discovery_score_exact"] != expected_score:
                errors.append(
                    f"{species} exact discovery-score mismatch for {row['gene_id']}"
                )
                break

    report: list[str] = ["===== UNBIASED DISCOVERY RANKING AUDIT =====", ""]
    report.append(f"Insect eligible ranked genes: {len(insect_full)}")
    report.append(f"Insect top-50 genes: {len(insect_top50)}")
    report.append(
        "Insect poorly characterized/uncharacterized pool: "
        + str(
            sum(
                row["annotation_status"]
                in {"poorly_characterized", "uncharacterized"}
                for row in insect_full
            )
        )
    )
    report.append(f"Insect focused poorly annotated top list: {len(insect_poor15)}")
    report.append(
        "Insect poorly annotated top-15 overlapping overall top-50: "
        + str(sum(row["in_overall_top50"] == "TRUE" for row in insect_poor15))
    )
    report.append("Insect full-ranking annotation status:")
    for status, count in sorted(Counter(row["annotation_status"] for row in insect_full).items()):
        report.append(f"  {status}: {count}")
    report.append("Insect top-50 discovery classification:")
    for status, count in sorted(
        Counter(row["discovery_classification"] for row in insect_top50).items()
    ):
        report.append(f"  {status}: {count}")
    report.append("")
    report.append(f"Plant eligible ranked genes: {len(plant_full)}")
    report.append(f"Plant top-50 genes: {len(plant_top50)}")
    report.append(
        "Plant GP-up genes without representative protein in full ranking: "
        + str(sum(row["representative_protein_available"] == "FALSE" for row in plant_full))
    )
    report.append("Plant full-ranking annotation status:")
    for status, count in sorted(Counter(row["annotation_status"] for row in plant_full).items()):
        report.append(f"  {status}: {count}")
    report.append("Plant top-50 discovery classification:")
    for status, count in sorted(
        Counter(row["discovery_classification"] for row in plant_top50).items()
    ):
        report.append(f"  {status}: {count}")
    report.append("Plant full-ranking MEROPS scope:")
    for status, count in sorted(Counter(row["MEROPS_screen_scope"] for row in plant_full).items()):
        report.append(f"  {status}: {count}")
    report.append("")
    report.append("Frozen ranked-list checksums:")
    for row in checksum_rows:
        report.append(f"  {row['file']}: {row['sha256']}")
    report.append("")
    report.append("Literature searches performed before ranking freeze: 0")
    report.append("Post-ranking manual membership changes: 0")
    report.append("Ranking arithmetic: exact rational ranks; no pre-sort rounding")
    report.append("Exact deterministic-sort violations: 0" if not any(
        "sort order" in error for error in errors
    ) else "Exact deterministic-sort violations: detected")
    report.append(f"Errors detected: {len(errors)}")
    for error in errors:
        report.append(f"ERROR: {error}")
    report.append(
        "UNBIASED DISCOVERY RANKING: PASS"
        if not errors
        else "UNBIASED DISCOVERY RANKING: CHECK REQUIRED"
    )

    audit_path = log_dir / "02_unbiased_discovery_ranking_audit.txt"
    audit_text = "\n".join(report) + "\n"
    audit_path.write_text(audit_text, encoding="utf-8")
    print(audit_text, end="")
    print(f"\nInsect full: {insect_full_path}")
    print(f"Insect top 50: {insect_top50_path}")
    print(f"Insect poorly annotated top 15: {insect_poor15_path}")
    print(f"Plant full: {plant_full_path}")
    print(f"Plant top 50: {plant_top50_path}")
    print(f"Parameters: {parameter_path}")
    print(f"Checksums: {checksum_path}")
    print(f"Audit: {audit_path}")
    return 1 if errors else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        raise
