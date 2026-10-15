#!/usr/bin/env python3
"""Build the frozen host-focused E. variegata GP-up discovery ranking.

This script does not fit DESeq2, alter DEG calls, rerun sequence-prediction
software, or apply a biological candidate-family filter.  It ranks the
already-defined host-focused GP-up genes with the same equal-weight method
used for the original discovery scan.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


EXPECTED_HOST_ROWS = 25_779
EXPECTED_GP_UP = 3_179
PADJ_THRESHOLD = 0.01
LFC_THRESHOLD = 1.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input has no header: {path}")
        rows = [dict(row) for row in reader]
        return list(reader.fieldnames), rows


def require_columns(path: Path, fields: list[str], required: set[str]) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "na", "nan", "none", "null", "-", "--"}:
        return ""
    return text


def numeric(value: object) -> float | None:
    text = clean(value)
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def unique_join(values: list[object], separator: str = ";") -> str:
    found: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = clean(value)
        if not text:
            continue
        for item in re.split(r"[;,]", text):
            item = item.strip()
            if item and item not in seen:
                seen.add(item)
                found.append(item)
    return separator.join(found)


def index_unique(
    path: Path, fields: list[str], rows: list[dict[str, str]], key: str = "gene_id"
) -> dict[str, dict[str, str]]:
    require_columns(path, fields, {key})
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = clean(row.get(key))
        if not identifier:
            raise ValueError(f"Blank {key} in {path}")
        if identifier in result:
            raise ValueError(f"Duplicate {key} in {path}: {identifier}")
        result[identifier] = row
    return result


def aggregate_rows(
    path: Path,
    fields: list[str],
    rows: list[dict[str, str]],
    requested: dict[str, str],
    key: str = "gene_id",
) -> dict[str, dict[str, str]]:
    """Aggregate selected columns from a possibly one-to-many evidence table."""
    require_columns(path, fields, {key})
    buckets: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: {output: [] for output in requested.values()}
    )
    for row in rows:
        identifier = clean(row.get(key))
        if not identifier:
            continue
        for source, output in requested.items():
            if source in row:
                buckets[identifier][output].append(row.get(source, ""))
    return {
        identifier: {column: unique_join(values) for column, values in columns.items()}
        for identifier, columns in buckets.items()
    }


def percentile_rank(values: list[float]) -> list[float]:
    """Average-tie percentile ranks, smallest=1/n and largest=1."""
    n = len(values)
    order = sorted(range(n), key=lambda index: values[index])
    output = [0.0] * n
    start = 0
    while start < n:
        end = start + 1
        while end < n and values[order[end]] == values[order[start]]:
            end += 1
        average_one_based_rank = ((start + 1) + end) / 2.0
        percentile = average_one_based_rank / n
        for position in range(start, end):
            output[order[position]] = percentile
        start = end
    return output


def annotation_status(row: dict[str, str]) -> str:
    description = clean(row.get("Description"))
    preferred = clean(row.get("Preferred_name"))
    pfams = clean(row.get("PFAMs"))
    combined = " ".join([description, preferred, pfams]).strip().lower()

    uncharacterized_patterns = (
        "uncharacterized protein",
        "hypothetical protein",
        "unknown protein",
        "protein of unknown function",
        "function unknown",
    )
    if not combined or any(pattern in combined for pattern in uncharacterized_patterns):
        return "uncharacterized"

    poorly_patterns = (
        "poorly characterized",
        "putative protein",
        "predicted protein",
        "conserved protein",
        "domain of unknown function",
    )
    duf_only = bool(pfams) and all(
        re.fullmatch(r"DUF\d+", token, flags=re.IGNORECASE)
        for token in re.split(r"[,;]", pfams)
        if token.strip()
    )
    if any(pattern in combined for pattern in poorly_patterns) or duf_only:
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
        "--plant-root",
        type=Path,
        default=Path.home() / "projects" / "qe_plant_denovo",
    )
    parser.add_argument(
        "--merops-root",
        type=Path,
        default=Path.home() / "projects" / "qe_merops_screen",
    )
    args = parser.parse_args()

    discovery_root = args.discovery_root.resolve()
    plant_root = args.plant_root.resolve()
    merops_root = args.merops_root.resolve()

    host_input = discovery_root / (
        "05_host_focused_DESeq2_sensitivity/results/"
        "plant_gene_DESeq2_GP_vs_UP_host_focused_eggnog_joined.tsv"
    )
    signal_input = plant_root / (
        "09_signalp_tmhmm/results/"
        "signalp6_deeptmhmm_plant_deseq2_primary_41330.tsv"
    )
    merops_input = merops_root / (
        "results/plant_GPup_secretome_MEROPS_gene_category_curated_v2.tsv"
    )
    defense_input = discovery_root / (
        "05_plant_defense_audit/results/"
        "plant_host_focused_defense_candidate_gene_summary.tsv"
    )
    developmental_input = plant_root / (
        "09_signalp_tmhmm/results/"
        "plant_developmental_cambial_candidate_audit_nonredundant.tsv"
    )
    counterpart_input = plant_root / (
        "09_signalp_tmhmm/results/"
        "plant_targeted_counterpart_family_taxonomic_annotation_audit.tsv"
    )

    results_dir = discovery_root / "05_host_focused_DESeq2_sensitivity/results"
    control_dir = discovery_root / "00_project_control"
    logs_dir = discovery_root / "05_host_focused_DESeq2_sensitivity/logs"

    full_output = results_dir / "plant_GP_host_focused_discovery_full.tsv"
    top50_output = results_dir / "plant_GP_host_focused_top50.tsv"
    poor15_output = results_dir / "plant_GP_host_focused_poorly_annotated_top15.tsv"
    parameters_output = control_dir / "20_host_focused_plant_discovery_parameters.tsv"
    audit_output = logs_dir / "20_host_focused_plant_discovery_ranking_audit.txt"

    required_inputs = [host_input, signal_input, merops_input]
    for path in required_inputs:
        if not path.is_file():
            raise FileNotFoundError(f"Required input is missing: {path}")

    host_fields, host_rows = read_tsv(host_input)
    require_columns(
        host_input,
        host_fields,
        {
            "gene_id",
            "baseMean",
            "log2FoldChange",
            "padj",
            "DESeq2_status_host_focused",
            "Description",
            "Preferred_name",
            "PFAMs",
        },
    )
    host_index = index_unique(host_input, host_fields, host_rows)
    if len(host_rows) != EXPECTED_HOST_ROWS:
        raise ValueError(
            f"Host-focused input has {len(host_rows)} rows; expected {EXPECTED_HOST_ROWS}"
        )

    signal_fields, signal_rows = read_tsv(signal_input)
    signal_index = index_unique(signal_input, signal_fields, signal_rows)

    merops_fields, merops_rows = read_tsv(merops_input)
    merops_index = aggregate_rows(
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

    optional_inputs: list[tuple[str, Path, dict[str, str]]] = [
        (
            "defense",
            defense_input,
            {
                "defense_groups": "defense_groups",
                "candidate_subgroups": "defense_candidate_subgroups",
            },
        ),
        (
            "developmental",
            developmental_input,
            {
                "target_group": "developmental_target_groups",
                "target_symbol": "developmental_target_symbols",
                "candidate_group": "developmental_target_groups",
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
    optional_indexes: dict[str, dict[str, dict[str, str]]] = {}
    optional_status: list[tuple[str, str, str]] = []
    for label, path, requested in optional_inputs:
        if path.is_file():
            fields, rows = read_tsv(path)
            present = {source: output for source, output in requested.items() if source in fields}
            optional_indexes[label] = aggregate_rows(path, fields, rows, present)
            optional_status.append((label, str(path), f"used; rows={len(rows)}"))
        else:
            optional_indexes[label] = {}
            optional_status.append((label, str(path), "not found; overlay left blank"))

    eligible: list[dict[str, object]] = []
    status_gp = sum(
        clean(row.get("DESeq2_status_host_focused")) == "GP_up" for row in host_rows
    )
    if status_gp != EXPECTED_GP_UP:
        raise ValueError(f"Host-focused GP_up status count is {status_gp}; expected {EXPECTED_GP_UP}")

    direction_mismatches: list[str] = []
    missing_signal: list[str] = []
    for gene_id, source in host_index.items():
        status = clean(source.get("DESeq2_status_host_focused"))
        if status != "GP_up":
            continue
        base_mean = numeric(source.get("baseMean"))
        lfc = numeric(source.get("log2FoldChange"))
        padj = numeric(source.get("padj"))
        if base_mean is None or lfc is None or padj is None:
            direction_mismatches.append(gene_id)
            continue
        if not (padj < PADJ_THRESHOLD and lfc >= LFC_THRESHOLD):
            direction_mismatches.append(gene_id)
            continue

        signal = signal_index.get(gene_id, {})
        if not signal:
            missing_signal.append(gene_id)
        merops = merops_index.get(gene_id, {})
        defense = optional_indexes["defense"].get(gene_id, {})
        developmental = optional_indexes["developmental"].get(gene_id, {})
        counterpart = optional_indexes["counterpart"].get(gene_id, {})

        row: dict[str, object] = {
            "gene_id": gene_id,
            "transcript_id": clean(source.get("transcript_id")),
            "protein_id": clean(source.get("protein_id")),
            "aa_length": clean(source.get("aa_length")),
            "baseMean": base_mean,
            "log2FoldChange": lfc,
            "lfcSE": clean(source.get("lfcSE")),
            "stat": clean(source.get("stat")),
            "pvalue": clean(source.get("pvalue")),
            "padj": padj,
            "DESeq2_status_host_focused": status,
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
            "protein_prediction_available": "TRUE" if signal else "FALSE",
            "signalp_prediction": clean(signal.get("signalp_prediction")),
            "signalp_SP_score": clean(signal.get("signalp_SP_score")),
            "signalp_cleavage_after": clean(signal.get("signalp_cleavage_after")),
            "deeptmhmm_type": clean(signal.get("deeptmhmm_type")),
            "predicted_TMRs": clean(signal.get("predicted_TMRs")),
            "has_TMR": clean(signal.get("has_TMR")),
            "TMR_segments": clean(signal.get("TMR_segments")),
            "analysis_class": clean(signal.get("analysis_class")),
            "primary_classical_secretome": clean(
                signal.get("primary_classical_secretome")
            ),
            "MEROPS_overlay_present": "TRUE" if merops else "FALSE",
            "MEROPS_identifiers": merops.get("MEROPS_identifiers", ""),
            "MEROPS_families": merops.get("MEROPS_families", ""),
            "MEROPS_categories": merops.get("MEROPS_categories", ""),
            "MEROPS_evidence_tiers": merops.get("MEROPS_evidence_tiers", ""),
            "MEROPS_manuscript_priority": merops.get(
                "MEROPS_manuscript_priority", ""
            ),
            "defense_groups": defense.get("defense_groups", ""),
            "defense_candidate_subgroups": defense.get(
                "defense_candidate_subgroups", ""
            ),
            "developmental_target_groups": developmental.get(
                "developmental_target_groups", ""
            ),
            "developmental_target_symbols": developmental.get(
                "developmental_target_symbols", ""
            ),
            "counterpart_candidate_families": counterpart.get(
                "counterpart_candidate_families", ""
            ),
        }
        row["metric_neg_log10_padj"] = -math.log10(max(padj, 1e-300))
        row["metric_abs_log2FoldChange"] = abs(lfc)
        row["metric_log10_baseMean_plus1"] = math.log10(base_mean + 1.0)
        eligible.append(row)

    if direction_mismatches:
        raise ValueError(
            "GP_up status/threshold disagreement for: "
            + ", ".join(direction_mismatches[:20])
        )
    if len(eligible) != EXPECTED_GP_UP:
        raise ValueError(f"Eligible GP-up rows are {len(eligible)}; expected {EXPECTED_GP_UP}")
    metric_columns = [
        "metric_neg_log10_padj",
        "metric_abs_log2FoldChange",
        "metric_log10_baseMean_plus1",
    ]
    percentile_columns = [
        "percentile_significance",
        "percentile_effect_size",
        "percentile_expression",
    ]
    for metric, percentile in zip(metric_columns, percentile_columns):
        ranks = percentile_rank([float(row[metric]) for row in eligible])
        for row, value in zip(eligible, ranks):
            row[percentile] = value

    for row in eligible:
        row["discovery_score"] = sum(float(row[col]) for col in percentile_columns) / 3.0

    eligible.sort(
        key=lambda row: (
            -float(row["discovery_score"]),
            float(row["padj"]),
            -float(row["metric_abs_log2FoldChange"]),
            -float(row["baseMean"]),
            str(row["gene_id"]),
        )
    )
    for rank, row in enumerate(eligible, start=1):
        row["discovery_rank"] = rank

    fields = [
        "discovery_rank",
        "discovery_score",
        "percentile_significance",
        "percentile_effect_size",
        "percentile_expression",
        "gene_id",
        "transcript_id",
        "protein_id",
        "aa_length",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "DESeq2_status_host_focused",
        "annotation_status",
        "Description",
        "Preferred_name",
        "seed_ortholog",
        "eggNOG_OGs",
        "PFAMs",
        "GOs",
        "KEGG_ko",
        "KEGG_Pathway",
        "CAZy",
        "protein_prediction_available",
        "signalp_prediction",
        "signalp_SP_score",
        "signalp_cleavage_after",
        "deeptmhmm_type",
        "predicted_TMRs",
        "has_TMR",
        "TMR_segments",
        "analysis_class",
        "primary_classical_secretome",
        "MEROPS_overlay_present",
        "MEROPS_identifiers",
        "MEROPS_families",
        "MEROPS_categories",
        "MEROPS_evidence_tiers",
        "MEROPS_manuscript_priority",
        "defense_groups",
        "defense_candidate_subgroups",
        "developmental_target_groups",
        "developmental_target_symbols",
        "counterpart_candidate_families",
        "metric_neg_log10_padj",
        "metric_abs_log2FoldChange",
        "metric_log10_baseMean_plus1",
    ]

    top50 = eligible[:50]
    poor_candidates = [
        row
        for row in eligible
        if row["annotation_status"] in {"poorly_characterized", "uncharacterized"}
    ]
    poor15 = poor_candidates[:15]

    write_tsv(full_output, eligible, fields)
    write_tsv(top50_output, top50, fields)
    write_tsv(poor15_output, poor15, fields)

    parameter_rows: list[dict[str, object]] = [
        {"parameter": "official_padj_rule", "value": "padj < 0.01"},
        {"parameter": "official_GP_direction_rule", "value": "log2FoldChange >= 1"},
        {
            "parameter": "ranking_method",
            "value": (
                "mean of average-tie percentile ranks for -log10(padj), "
                "abs(log2FoldChange), and log10(baseMean+1)"
            ),
        },
        {"parameter": "ranking_weights", "value": "equal; one third per component"},
        {"parameter": "rank_scope", "value": "host-focused GP-up genes only"},
        {"parameter": "top_list_size", "value": 50},
        {"parameter": "poor_annotation_review_size", "value": 15},
        {
            "parameter": "MEROPS_scope_caveat",
            "value": (
                "existing overlay restricted to original GP-up classical-secretome screen; "
                "absence is not negative evidence"
            ),
        },
    ]
    write_tsv(parameters_output, parameter_rows, ["parameter", "value"])

    annotation_counts = defaultdict(int)
    for row in eligible:
        annotation_counts[str(row["annotation_status"])] += 1
    top_annotation_counts = defaultdict(int)
    for row in top50:
        top_annotation_counts[str(row["annotation_status"])] += 1

    audit_lines = [
        "===== HOST-FOCUSED PLANT UNBIASED DISCOVERY RANKING =====",
        f"Discovery root: {discovery_root}",
        "",
        "Scope:",
        "  Ranking of the frozen host-focused GP-up DEGs; no DESeq2 fitting or DEG reclassification.",
        "  No SignalP, DeepTMHMM or MEROPS analysis was rerun.",
        "  Annotation and candidate-family evidence were descriptive and did not affect ranking.",
        "",
        "Validated inputs:",
        f"  Host-focused DESeq2 + eggNOG: {host_input}",
        f"    Rows: {len(host_rows)}; SHA-256: {sha256(host_input)}",
        f"  SignalP/DeepTMHMM: {signal_input}",
        f"    Rows: {len(signal_rows)}; SHA-256: {sha256(signal_input)}",
        f"  Existing plant MEROPS overlay: {merops_input}",
        f"    Rows: {len(merops_rows)}; unique genes: {len(merops_index)}; SHA-256: {sha256(merops_input)}",
        "",
        "Optional descriptive overlays:",
    ]
    for label, path, status in optional_status:
        audit_lines.append(f"  {label}: {path} [{status}]")
    audit_lines.extend(
        [
            "",
            f"Host-focused GP-up genes ranked: {len(eligible)}",
            f"GP-up genes with representative-protein predictions: {len(eligible) - len(missing_signal)}",
            f"GP-up genes without representative-protein predictions: {len(missing_signal)}",
            f"Top-50 rows: {len(top50)}",
            f"Poorly characterized/uncharacterized review rows: {len(poor15)}",
            "",
            "Annotation-status distribution among all ranked genes:",
        ]
    )
    for label in ("informative", "poorly_characterized", "uncharacterized"):
        audit_lines.append(f"  {label}: {annotation_counts[label]}")
    audit_lines.append("")
    audit_lines.append("Annotation-status distribution in frozen top 50:")
    for label in ("informative", "poorly_characterized", "uncharacterized"):
        audit_lines.append(f"  {label}: {top_annotation_counts[label]}")
    audit_lines.extend(
        [
            "",
            "Ranking formula:",
            "  discovery_score = mean(percentile[-log10(padj)], percentile[abs(log2FC)], percentile[log10(baseMean+1)])",
            "  Average ranks were used for ties; the highest value in each component approaches 1.",
            "",
            "Interpretive boundary:",
            "  The discovery score is a prioritization device, not a significance test.",
            "  The top 50 was selected mechanically before biological or literature interpretation.",
            "  Existing MEROPS absence is not negative evidence because the plant MEROPS screen had a restricted original scope.",
            "  Genes without a representative protein remained eligible; their SignalP/DeepTMHMM fields are blank.",
            "  The earlier mixed-provenance plant ranking remains unchanged and must not be used as a host-specific top-gene list.",
            "",
            "Errors detected: 0",
            "",
            f"Output: {full_output}",
            f"Output: {top50_output}",
            f"Output: {poor15_output}",
            f"Output: {parameters_output}",
            "HOST-FOCUSED PLANT UNBIASED DISCOVERY RANKING: PASS",
        ]
    )
    logs_dir.mkdir(parents=True, exist_ok=True)
    audit_output.write_text("\n".join(audit_lines) + "\n", encoding="utf-8")
    print("\n".join(audit_lines))
    print(f"Output: {audit_output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("===== HOST-FOCUSED PLANT UNBIASED DISCOVERY RANKING =====", file=sys.stderr)
        print(f"ERROR: {error}", file=sys.stderr)
        print("HOST-FOCUSED PLANT UNBIASED DISCOVERY RANKING: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
