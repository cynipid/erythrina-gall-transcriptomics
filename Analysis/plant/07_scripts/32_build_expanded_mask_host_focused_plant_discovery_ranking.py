#!/usr/bin/env python3
"""Build the final expanded-mask host-focused E. variegata GP-up ranking.

The ranking is restricted to the already-defined expanded-mask GP-up DEGs and
uses the frozen equal-weight formula from the preceding host-focused ranking.
Biological annotations and candidate-family overlays are descriptive only and
cannot affect eligibility or rank. No DESeq2 fitting, annotation search,
SignalP/DeepTMHMM run, MEROPS search, or candidate-family search is performed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_DESEQ_ROWS = 24_931
EXPECTED_GP_UP = 2_349
EXPECTED_PREVIOUS_GP_UP = 3_179
EXPECTED_SHARED_GP_UP = 2_343
EXPECTED_PREVIOUS_ONLY = 836
EXPECTED_CURRENT_ONLY = 6
EXPECTED_PREVIOUS_ONLY_MASKED = 813
PADJ_THRESHOLD = 0.01
LFC_THRESHOLD = 1.0
STATUS_COLUMN = "DESeq2_status_host_focused_expanded_mask"
MASK_EXCLUDE_STATUSES = {
    "review_exclude_existing_confirmed_multilevel_sequence",
    "review_exclude_additional_strong_nucleotide_similarity",
}


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "na", "nan", "none", "null", "-", "--"} else text


def number(value: object) -> float | None:
    text = clean(value)
    if not text:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input has no header: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def require_columns(path: Path, fields: list[str], required: set[str]) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")


def unique_index(
    path: Path, fields: list[str], rows: list[dict[str, str]], key: str = "gene_id"
) -> dict[str, dict[str, str]]:
    require_columns(path, fields, {key})
    output: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = clean(row.get(key))
        if not identifier:
            raise ValueError(f"Blank {key} in {path}")
        if identifier in output:
            raise ValueError(f"Duplicate {key} in {path}: {identifier}")
        output[identifier] = row
    return output


def unique_join(values: list[object]) -> str:
    output: list[str] = []
    seen: set[str] = set()
    for value in values:
        for item in re.split(r"[;,]", clean(value)):
            item = item.strip()
            if item and item not in seen:
                seen.add(item)
                output.append(item)
    return ";".join(output)


def aggregate_rows(
    path: Path,
    fields: list[str],
    rows: list[dict[str, str]],
    requested: dict[str, str],
) -> dict[str, dict[str, str]]:
    require_columns(path, fields, {"gene_id"})
    present = {source: target for source, target in requested.items() if source in fields}
    buckets: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: {target: [] for target in set(present.values())}
    )
    for row in rows:
        gene = clean(row.get("gene_id"))
        if not gene:
            continue
        for source, target in present.items():
            buckets[gene][target].append(row.get(source, ""))
    return {
        gene: {target: unique_join(values) for target, values in columns.items()}
        for gene, columns in buckets.items()
    }


def percentile_rank(values: list[float]) -> list[float]:
    """Return average-tie percentile ranks, smallest=1/n and largest=1."""
    size = len(values)
    order = sorted(range(size), key=lambda index: values[index])
    output = [0.0] * size
    start = 0
    while start < size:
        end = start + 1
        while end < size and values[order[end]] == values[order[start]]:
            end += 1
        percentile = (((start + 1) + end) / 2.0) / size
        for position in range(start, end):
            output[order[position]] = percentile
        start = end
    return output


def annotation_status(row: dict[str, str]) -> str:
    description = clean(row.get("Description"))
    preferred = clean(row.get("Preferred_name"))
    pfams = clean(row.get("PFAMs"))
    combined = " ".join([description, preferred, pfams]).strip().lower()
    if not combined or any(
        pattern in combined
        for pattern in (
            "uncharacterized protein",
            "hypothetical protein",
            "unknown protein",
            "protein of unknown function",
            "function unknown",
        )
    ):
        return "uncharacterized"
    duf_tokens = [token.strip() for token in re.split(r"[,;]", pfams) if token.strip()]
    duf_only = bool(duf_tokens) and all(
        re.fullmatch(r"DUF\d+", token, flags=re.IGNORECASE) for token in duf_tokens
    )
    if duf_only or any(
        pattern in combined
        for pattern in (
            "poorly characterized",
            "putative protein",
            "predicted protein",
            "conserved protein",
            "domain of unknown function",
        )
    ):
        return "poorly_characterized"
    return "informative"


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--discovery-root",
        type=Path,
        default=Path.home() / "projects" / "qe_cross_species_discovery",
    )
    parser.add_argument(
        "--plant-root", type=Path, default=Path.home() / "projects" / "qe_plant_denovo"
    )
    parser.add_argument(
        "--merops-root", type=Path, default=Path.home() / "projects" / "qe_merops_screen"
    )
    args = parser.parse_args()

    root = args.discovery_root.resolve()
    plant_root = args.plant_root.resolve()
    merops_root = args.merops_root.resolve()
    results = root / "05_host_focused_DESeq2_sensitivity" / "results"
    logs = root / "05_host_focused_DESeq2_sensitivity" / "logs"
    control = root / "00_project_control"

    deseq_input = results / "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_eggnog_joined.tsv"
    signal_input = plant_root / "09_signalp_tmhmm/results/signalp6_deeptmhmm_plant_deseq2_primary_41330.tsv"
    merops_input = merops_root / "results/plant_GPup_secretome_MEROPS_gene_category_curated_v2.tsv"
    defense_input = root / "05_plant_defense_audit/results/plant_expanded_mask_defense_candidate_gene_summary.tsv"
    developmental_input = plant_root / "09_signalp_tmhmm/results/plant_developmental_cambial_candidate_audit_nonredundant.tsv"
    counterpart_input = plant_root / "09_signalp_tmhmm/results/plant_targeted_counterpart_family_taxonomic_annotation_audit.tsv"
    previous_ranking_input = results / "plant_GP_host_focused_discovery_full.tsv"
    mask_input = root / "04_cross_species_integration/results/plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv"

    full_output = results / "plant_GP_host_focused_expanded_mask_discovery_full.tsv"
    top50_output = results / "plant_GP_host_focused_expanded_mask_top50.tsv"
    poor15_output = results / "plant_GP_host_focused_expanded_mask_poorly_annotated_top15.tsv"
    comparison_output = results / "plant_GP_host_focused_current_vs_expanded_mask_discovery_rank_comparison.tsv"
    parameters_output = control / "32_expanded_mask_host_focused_plant_discovery_parameters.tsv"
    audit_output = logs / "32_expanded_mask_host_focused_plant_discovery_ranking_audit.txt"

    if audit_output.is_file():
        previous_audit = audit_output.read_text(encoding="utf-8").splitlines()
        if previous_audit and previous_audit[-1] == "EXPANDED-MASK HOST-FOCUSED PLANT UNBIASED DISCOVERY RANKING: PASS":
            print("Validated expanded-mask discovery-ranking outputs already exist; no files were changed.")
            return 1

    required_inputs = [
        deseq_input,
        signal_input,
        merops_input,
        defense_input,
        previous_ranking_input,
        mask_input,
    ]
    input_hashes = {str(path): sha256(path) for path in required_inputs}

    deseq_fields, deseq_rows = read_tsv(deseq_input)
    require_columns(
        deseq_input,
        deseq_fields,
        {
            "gene_id", "baseMean", "log2FoldChange", "padj", STATUS_COLUMN,
            "Description", "Preferred_name", "PFAMs",
        },
    )
    deseq = unique_index(deseq_input, deseq_fields, deseq_rows)
    if len(deseq_rows) != EXPECTED_DESEQ_ROWS:
        raise ValueError(f"Expanded-mask DESeq2 rows: {len(deseq_rows)}; expected {EXPECTED_DESEQ_ROWS}")
    expected_statuses = Counter(
        {"GP_up": 2349, "UP_up": 1386, "non_DE": 21022, "NA_padj_not_testable": 174}
    )
    observed_statuses = Counter(clean(row.get(STATUS_COLUMN)) for row in deseq_rows)
    if observed_statuses != expected_statuses:
        raise ValueError(f"Expanded-mask DESeq2 statuses differ: {dict(observed_statuses)}")

    signal_fields, signal_rows = read_tsv(signal_input)
    signal = unique_index(signal_input, signal_fields, signal_rows)
    merops_fields, merops_rows = read_tsv(merops_input)
    merops = aggregate_rows(
        merops_input,
        merops_fields,
        merops_rows,
        {
            "merops_identifier": "MEROPS_identifiers",
            "merops_family": "MEROPS_families",
            "merops_category": "MEROPS_categories",
            "evidence_tier": "MEROPS_evidence_tiers",
            "manuscript_priority": "MEROPS_manuscript_priority",
        },
    )
    defense_fields, defense_rows = read_tsv(defense_input)
    defense = aggregate_rows(
        defense_input,
        defense_fields,
        defense_rows,
        {"defense_groups": "defense_groups", "candidate_subgroups": "defense_candidate_subgroups"},
    )
    if len(defense) != 806:
        raise ValueError(f"Expanded-mask defense overlay has {len(defense)} unique genes; expected 806")

    optional_specs = [
        (
            "developmental",
            developmental_input,
            {
                "target_group": "developmental_target_groups",
                "candidate_group": "developmental_target_groups",
                "target_symbol": "developmental_target_symbols",
                "candidate_name": "developmental_target_symbols",
            },
        ),
        (
            "counterpart",
            counterpart_input,
            {
                "candidate_family": "counterpart_candidate_families",
                "family": "counterpart_candidate_families",
                "candidate_group": "counterpart_candidate_families",
            },
        ),
    ]
    optional: dict[str, dict[str, dict[str, str]]] = {}
    optional_status: list[str] = []
    for label, path, requested in optional_specs:
        if path.is_file():
            fields, rows = read_tsv(path)
            optional[label] = aggregate_rows(path, fields, rows, requested)
            input_hashes[str(path)] = sha256(path)
            optional_status.append(f"  {label}: {path} [used; rows={len(rows)}]")
        else:
            optional[label] = {}
            optional_status.append(f"  {label}: {path} [not found; overlay left blank]")

    mask_fields, mask_rows = read_tsv(mask_input)
    require_columns(mask_input, mask_fields, {"gene_id", "proposed_expanded_mask_status"})
    mask = unique_index(mask_input, mask_fields, mask_rows)
    if len(mask) != 118_082:
        raise ValueError(f"Expanded-mask universe has {len(mask)} genes; expected 118082")
    excluded = {
        gene
        for gene, row in mask.items()
        if clean(row.get("proposed_expanded_mask_status")) in MASK_EXCLUDE_STATUSES
    }
    if len(excluded) != 16_112:
        raise ValueError(f"Expanded exclusion mask has {len(excluded)} genes; expected 16112")
    if set(deseq) & excluded:
        raise ValueError("One or more expanded-mask genes remain in the final DESeq2 table")

    eligible: list[dict[str, object]] = []
    missing_signal: list[str] = []
    for gene, source in deseq.items():
        if clean(source.get(STATUS_COLUMN)) != "GP_up":
            continue
        base_mean = number(source.get("baseMean"))
        lfc = number(source.get("log2FoldChange"))
        padj = number(source.get("padj"))
        if base_mean is None or lfc is None or padj is None:
            raise ValueError(f"GP-up gene lacks finite DESeq2 values: {gene}")
        if not (padj < PADJ_THRESHOLD and lfc >= LFC_THRESHOLD):
            raise ValueError(f"GP-up status disagrees with frozen thresholds: {gene}")

        signal_row = signal.get(gene, {})
        if not signal_row:
            missing_signal.append(gene)
        merops_row = merops.get(gene, {})
        defense_row = defense.get(gene, {})
        developmental_row = optional["developmental"].get(gene, {})
        counterpart_row = optional["counterpart"].get(gene, {})
        row: dict[str, object] = {
            "gene_id": gene,
            "transcript_id": clean(source.get("transcript_id")),
            "protein_id": clean(source.get("protein_id")),
            "aa_length": clean(source.get("aa_length")),
            "baseMean": base_mean,
            "log2FoldChange": lfc,
            "lfcSE": clean(source.get("lfcSE")),
            "stat": clean(source.get("stat")),
            "pvalue": clean(source.get("pvalue")),
            "padj": padj,
            STATUS_COLUMN: "GP_up",
            "Description": clean(source.get("Description")),
            "Preferred_name": clean(source.get("Preferred_name")),
            "seed_ortholog": clean(source.get("seed_ortholog")),
            "eggNOG_OGs": clean(source.get("eggNOG_OGs")),
            "GOs": clean(source.get("GOs")),
            "KEGG_ko": clean(source.get("KEGG_ko")),
            "KEGG_Pathway": clean(source.get("KEGG_Pathway")),
            "CAZy": clean(source.get("CAZy")),
            "PFAMs": clean(source.get("PFAMs")),
            "annotation_status": annotation_status(source),
            "protein_prediction_available": "TRUE" if signal_row else "FALSE",
            "signalp_prediction": clean(signal_row.get("signalp_prediction")),
            "signalp_SP_score": clean(signal_row.get("signalp_SP_score")),
            "signalp_cleavage_after": clean(signal_row.get("signalp_cleavage_after")),
            "deeptmhmm_type": clean(signal_row.get("deeptmhmm_type")),
            "predicted_TMRs": clean(signal_row.get("predicted_TMRs")),
            "has_TMR": clean(signal_row.get("has_TMR")),
            "TMR_segments": clean(signal_row.get("TMR_segments")),
            "analysis_class": clean(signal_row.get("analysis_class")),
            "primary_classical_secretome": clean(signal_row.get("primary_classical_secretome")),
            "MEROPS_overlay_present": "TRUE" if merops_row else "FALSE",
            "MEROPS_identifiers": merops_row.get("MEROPS_identifiers", ""),
            "MEROPS_families": merops_row.get("MEROPS_families", ""),
            "MEROPS_categories": merops_row.get("MEROPS_categories", ""),
            "MEROPS_evidence_tiers": merops_row.get("MEROPS_evidence_tiers", ""),
            "MEROPS_manuscript_priority": merops_row.get("MEROPS_manuscript_priority", ""),
            "defense_groups": defense_row.get("defense_groups", ""),
            "defense_candidate_subgroups": defense_row.get("defense_candidate_subgroups", ""),
            "developmental_target_groups": developmental_row.get("developmental_target_groups", ""),
            "developmental_target_symbols": developmental_row.get("developmental_target_symbols", ""),
            "counterpart_candidate_families": counterpart_row.get("counterpart_candidate_families", ""),
            "metric_neg_log10_padj": -math.log10(max(padj, 1e-300)),
            "metric_abs_log2FoldChange": abs(lfc),
            "metric_log10_baseMean_plus1": math.log10(base_mean + 1.0),
        }
        eligible.append(row)

    if len(eligible) != EXPECTED_GP_UP:
        raise ValueError(f"Eligible expanded-mask GP-up genes: {len(eligible)}; expected {EXPECTED_GP_UP}")
    metric_columns = [
        "metric_neg_log10_padj", "metric_abs_log2FoldChange", "metric_log10_baseMean_plus1"
    ]
    percentile_columns = [
        "percentile_significance", "percentile_effect_size", "percentile_expression"
    ]
    for metric, percentile in zip(metric_columns, percentile_columns):
        ranks = percentile_rank([float(row[metric]) for row in eligible])
        for row, rank in zip(eligible, ranks):
            row[percentile] = rank
    for row in eligible:
        row["discovery_score"] = sum(float(row[field]) for field in percentile_columns) / 3.0
    eligible.sort(
        key=lambda row: (
            -float(row["discovery_score"]), float(row["padj"]),
            -float(row["metric_abs_log2FoldChange"]), -float(row["baseMean"]), str(row["gene_id"]),
        )
    )
    for rank, row in enumerate(eligible, start=1):
        row["discovery_rank"] = rank

    fields = [
        "discovery_rank", "discovery_score", "percentile_significance",
        "percentile_effect_size", "percentile_expression", "gene_id", "transcript_id",
        "protein_id", "aa_length", "baseMean", "log2FoldChange", "lfcSE", "stat",
        "pvalue", "padj", STATUS_COLUMN, "annotation_status", "Description",
        "Preferred_name", "seed_ortholog", "eggNOG_OGs", "PFAMs", "GOs", "KEGG_ko",
        "KEGG_Pathway", "CAZy", "protein_prediction_available", "signalp_prediction",
        "signalp_SP_score", "signalp_cleavage_after", "deeptmhmm_type", "predicted_TMRs",
        "has_TMR", "TMR_segments", "analysis_class", "primary_classical_secretome",
        "MEROPS_overlay_present", "MEROPS_identifiers", "MEROPS_families",
        "MEROPS_categories", "MEROPS_evidence_tiers", "MEROPS_manuscript_priority",
        "defense_groups", "defense_candidate_subgroups", "developmental_target_groups",
        "developmental_target_symbols", "counterpart_candidate_families",
        "metric_neg_log10_padj", "metric_abs_log2FoldChange", "metric_log10_baseMean_plus1",
    ]
    top50 = eligible[:50]
    poor15 = [
        row for row in eligible
        if row["annotation_status"] in {"poorly_characterized", "uncharacterized"}
    ][:15]
    if len(top50) != 50 or len(poor15) != 15:
        raise ValueError(f"Expected 50 top rows and 15 poorly annotated rows; observed {len(top50)} and {len(poor15)}")

    previous_fields, previous_rows = read_tsv(previous_ranking_input)
    require_columns(
        previous_ranking_input,
        previous_fields,
        {"gene_id", "discovery_rank", "discovery_score", "DESeq2_status_host_focused"},
    )
    previous = unique_index(previous_ranking_input, previous_fields, previous_rows)
    if len(previous) != EXPECTED_PREVIOUS_GP_UP:
        raise ValueError(f"Previous ranking has {len(previous)} rows; expected {EXPECTED_PREVIOUS_GP_UP}")
    current = {str(row["gene_id"]): row for row in eligible}
    shared = set(previous) & set(current)
    previous_only = set(previous) - set(current)
    current_only = set(current) - set(previous)
    if (len(shared), len(previous_only), len(current_only)) != (
        EXPECTED_SHARED_GP_UP, EXPECTED_PREVIOUS_ONLY, EXPECTED_CURRENT_ONLY
    ):
        raise ValueError(
            "Previous/current GP-up membership comparison differs: "
            f"shared={len(shared)}, previous_only={len(previous_only)}, current_only={len(current_only)}"
        )
    previous_only_masked = previous_only & excluded
    if len(previous_only_masked) != EXPECTED_PREVIOUS_ONLY_MASKED:
        raise ValueError(
            f"Previous-only GP-up genes directly excluded by expanded mask: {len(previous_only_masked)}; expected 813"
        )
    if len(previous_only - excluded) != 23:
        raise ValueError(
            f"Previous GP-up genes retained but no longer GP-up after refit: {len(previous_only - excluded)}; expected 23"
        )

    previous_top50 = {
        gene for gene, row in previous.items() if number(row.get("discovery_rank")) is not None
        and float(row["discovery_rank"]) <= 50
    }
    current_top50 = {str(row["gene_id"]) for row in top50}
    if len(previous_top50) != 50 or len(current_top50) != 50:
        raise ValueError(
            f"Top-50 membership sizes differ: previous={len(previous_top50)}, expanded={len(current_top50)}"
        )
    shared_top50 = previous_top50 & current_top50

    comparison_rows: list[dict[str, object]] = []
    for gene in sorted(set(previous) | set(current)):
        old = previous.get(gene, {})
        new = current.get(gene, {})
        if gene in shared:
            transition = "retained_GP_up"
        elif gene in excluded:
            transition = "GP_up_to_excluded_by_expanded_mask"
        elif gene in previous_only:
            transition = f"GP_up_to_{clean(deseq.get(gene, {}).get(STATUS_COLUMN)) or 'not_in_final_tested_universe'}"
        else:
            transition = "new_GP_up_after_refit"
        old_rank = number(old.get("discovery_rank"))
        new_rank = number(new.get("discovery_rank"))
        comparison_rows.append({
            "gene_id": gene,
            "comparison_transition": transition,
            "previous_GP_up_member": "TRUE" if gene in previous else "FALSE",
            "expanded_mask_GP_up_member": "TRUE" if gene in current else "FALSE",
            "expanded_nucleotide_mask_member": "TRUE" if gene in excluded else "FALSE",
            "previous_discovery_rank": clean(old.get("discovery_rank")),
            "expanded_mask_discovery_rank": clean(new.get("discovery_rank")),
            "rank_change_positive_means_rise": (old_rank - new_rank) if old_rank is not None and new_rank is not None else "",
            "previous_discovery_score": clean(old.get("discovery_score")),
            "expanded_mask_discovery_score": clean(new.get("discovery_score")),
            "expanded_mask_DESeq2_status": (
                "excluded_by_expanded_nucleotide_mask"
                if gene in excluded else clean(deseq.get(gene, {}).get(STATUS_COLUMN))
            ),
        })
    comparison_rows.sort(
        key=lambda row: (
            0 if row["expanded_mask_GP_up_member"] == "TRUE" else 1,
            number(row["expanded_mask_discovery_rank"]) or math.inf,
            number(row["previous_discovery_rank"]) or math.inf,
            str(row["gene_id"]),
        )
    )

    write_tsv(full_output, eligible, fields)
    write_tsv(top50_output, top50, fields)
    write_tsv(poor15_output, poor15, fields)
    write_tsv(comparison_output, comparison_rows, list(comparison_rows[0]))
    parameter_rows: list[dict[str, object]] = [
        {"parameter": "official_padj_rule", "value": "padj < 0.01"},
        {"parameter": "official_GP_direction_rule", "value": "log2FoldChange >= 1"},
        {"parameter": "rank_scope", "value": "expanded-mask host-focused GP-up genes only"},
        {
            "parameter": "ranking_method",
            "value": "mean of average-tie percentile ranks for -log10(padj), abs(log2FoldChange), and log10(baseMean+1)",
        },
        {"parameter": "ranking_weights", "value": "equal; one third per component"},
        {"parameter": "top_list_size", "value": 50},
        {"parameter": "poor_annotation_review_size", "value": 15},
        {
            "parameter": "MEROPS_scope_caveat",
            "value": "retained overlay came from original GP-up classical-secretome-restricted screen; absence is not negative evidence",
        },
    ]
    write_tsv(parameters_output, parameter_rows, ["parameter", "value"])

    if any(sha256(Path(path)) != digest for path, digest in input_hashes.items()):
        raise ValueError("One or more inputs changed during output generation")
    annotation_counts = Counter(str(row["annotation_status"]) for row in eligible)
    top_annotation_counts = Counter(str(row["annotation_status"]) for row in top50)
    comparison_counts = Counter(str(row["comparison_transition"]) for row in comparison_rows)
    lines = [
        "===== EXPANDED-MASK HOST-FOCUSED PLANT UNBIASED DISCOVERY RANKING =====",
        "Scope: final expanded-mask GP-up ranking; no DESeq2 fitting or DEG reclassification.",
        "No annotation search, candidate-family search, SignalP/DeepTMHMM run or MEROPS search was performed.",
        "Annotation and candidate overlays were descriptive and did not affect the ranking.",
        "",
        f"Expanded-mask DESeq2 genes: {len(deseq_rows)}",
        f"Expanded-mask GP-up genes ranked: {len(eligible)}",
        f"GP-up genes with representative-protein predictions: {len(eligible) - len(missing_signal)}",
        f"GP-up genes without representative-protein predictions: {len(missing_signal)}",
        f"Top-50 rows: {len(top50)}",
        f"Poorly characterized/uncharacterized review rows: {len(poor15)}",
        "",
        "Previous-to-expanded ranking membership:",
        f"  Retained GP-up genes: {len(shared)}",
        f"  Previous GP-up genes absent from final GP-up set: {len(previous_only)}",
        f"  Directly excluded by expanded nucleotide mask: {len(previous_only_masked)}",
        f"  Retained but no longer GP-up after refit: {len(previous_only - excluded)}",
        f"  Newly GP-up after refit: {len(current_only)}",
        f"  Genes retained in both top-50 lists: {len(shared_top50)}",
        f"  Genes leaving the previous top 50: {len(previous_top50 - current_top50)}",
        f"  Genes entering the expanded-mask top 50: {len(current_top50 - previous_top50)}",
        "",
        "Transition inventory:",
    ]
    for label in sorted(comparison_counts):
        lines.append(f"  {label}: {comparison_counts[label]}")
    lines.extend(["", "Annotation-status distribution among all ranked genes:"])
    for label in ("informative", "poorly_characterized", "uncharacterized"):
        lines.append(f"  {label}: {annotation_counts[label]}")
    lines.append("")
    lines.append("Annotation-status distribution in top 50:")
    for label in ("informative", "poorly_characterized", "uncharacterized"):
        lines.append(f"  {label}: {top_annotation_counts[label]}")
    lines.extend(
        [
            "",
            "Ranking formula:",
            "  discovery_score = mean(percentile[-log10(padj)], percentile[abs(log2FC)], percentile[log10(baseMean+1)])",
            "  Average ranks were used for ties; the highest value in each component approaches 1.",
            "",
            "Optional descriptive overlays:",
            *optional_status,
            "",
            "Interpretive boundary:",
            "  The discovery score is a prioritization device, not a significance test.",
            "  The top 50 was selected mechanically before biological or literature interpretation.",
            "  MEROPS absence remains non-negative evidence because the retained plant screen had restricted original scope.",
            "  Genes without representative proteins remained eligible; protein-prediction fields are blank.",
            "  Earlier mixed-provenance and 3,948-mask rankings remain preserved as provenance/sensitivity stages.",
            "",
            "Input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "EXPANDED-MASK HOST-FOCUSED PLANT UNBIASED DISCOVERY RANKING: PASS",
        ]
    )
    logs.mkdir(parents=True, exist_ok=True)
    audit_output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    for path in [full_output, top50_output, poor15_output, comparison_output, parameters_output, audit_output]:
        print(f"Output: {path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("===== EXPANDED-MASK HOST-FOCUSED PLANT UNBIASED DISCOVERY RANKING =====", file=sys.stderr)
        print(f"ERROR: {error}", file=sys.stderr)
        print("EXPANDED-MASK HOST-FOCUSED PLANT UNBIASED DISCOVERY RANKING: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
