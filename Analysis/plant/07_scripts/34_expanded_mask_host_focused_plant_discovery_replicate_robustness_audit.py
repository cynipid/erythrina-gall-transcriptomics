#!/usr/bin/env python3
"""Audit replicate-level support for final expanded-mask plant discovery lists.

The audit uses DESeq2-normalized counts descriptively. It performs no new
hypothesis test, changes no DEG call, removes no candidate, and does not alter
or recompute the frozen discovery ranking.
"""

from __future__ import annotations

import csv
import hashlib
import math
import statistics
import sys
from collections import Counter
from pathlib import Path


ROOT = Path.home() / "projects" / "qe_cross_species_discovery"
RESULTS = ROOT / "05_host_focused_DESeq2_sensitivity" / "results"
CONTROL = ROOT / "00_project_control"
LOGS = ROOT / "05_host_focused_DESeq2_sensitivity" / "logs"

TOP50 = RESULTS / "plant_GP_host_focused_expanded_mask_top50.tsv"
POOR15 = RESULTS / "plant_GP_host_focused_expanded_mask_poorly_annotated_top15.tsv"
COUNTS = RESULTS / "plant_gene_normalized_counts_host_focused_expanded_mask.tsv"

TOP50_OUTPUT = RESULTS / "plant_GP_host_focused_expanded_mask_top50_replicate_audit.tsv"
POOR15_OUTPUT = RESULTS / "plant_GP_host_focused_expanded_mask_poorly_annotated_top15_replicate_audit.tsv"
MAPPING_OUTPUT = CONTROL / "34_expanded_mask_host_focused_discovery_sample_mapping.tsv"
SUMMARY_OUTPUT = CONTROL / "34_expanded_mask_host_focused_plant_discovery_replicate_robustness_summary.tsv"
AUDIT_OUTPUT = LOGS / "34_expanded_mask_host_focused_plant_discovery_replicate_robustness_audit.txt"

EXPECTED_COUNT_ROWS = 24_931
EXPECTED_TOP50_ROWS = 50
EXPECTED_POOR15_ROWS = 15
PSEUDOCOUNT = 0.5
DOMINANCE_REVIEW_THRESHOLD = 0.80

SAMPLE_MAP = {
    "LTS20_VA01_Sample_S7": "GP1",
    "LTS20_VA01_Sample_S8": "GP2",
    "LTS20_VA01_Sample_S9": "GP3",
    "LTS20_VA01_Sample_S10": "UP1",
    "LTS20_VA01_Sample_S11": "UP2",
    "LTS20_VA01_Sample_S12": "UP3",
}
GP_SAMPLES = list(SAMPLE_MAP)[:3]
UP_SAMPLES = list(SAMPLE_MAP)[3:]


class AuditError(RuntimeError):
    pass


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip().lstrip("\ufeff")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise AuditError(f"Required input is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def require_columns(path: Path, fields: list[str], required: set[str]) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise AuditError(f"Missing columns in {path}: {', '.join(missing)}")


def unique_index(
    path: Path, fields: list[str], rows: list[dict[str, str]], key: str = "gene_id"
) -> dict[str, dict[str, str]]:
    require_columns(path, fields, {key})
    output: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = clean(row.get(key))
        if not identifier:
            raise AuditError(f"Blank {key} in {path}")
        if identifier in output:
            raise AuditError(f"Duplicate {key} in {path}: {identifier}")
        output[identifier] = row
    return output


def number(value: object, label: str) -> float:
    try:
        result = float(clean(value))
    except (TypeError, ValueError) as exc:
        raise AuditError(f"Non-numeric value for {label}: {value!r}") from exc
    if not math.isfinite(result) or result < 0:
        raise AuditError(f"Invalid normalized count for {label}: {value!r}")
    return result


def truth(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def fmt(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def mean(values: list[float]) -> float:
    return statistics.fmean(values)


def coefficient_of_variation(values: list[float]) -> float | None:
    group_mean = mean(values)
    return None if group_mean == 0 else statistics.stdev(values) / group_mean


def dominance(values: list[float], samples: list[str]) -> tuple[str, float]:
    total = sum(values)
    maximum_index = max(range(len(values)), key=lambda index: values[index])
    fraction = values[maximum_index] / total if total > 0 else 0.0
    return SAMPLE_MAP[samples[maximum_index]], fraction


def audit_candidate(
    candidate: dict[str, str], count_row: dict[str, str]
) -> dict[str, object]:
    gene = candidate["gene_id"]
    gp = [number(count_row[sample], f"{gene}/{sample}") for sample in GP_SAMPLES]
    up = [number(count_row[sample], f"{gene}/{sample}") for sample in UP_SAMPLES]
    gp_mean = mean(gp)
    up_mean = mean(up)
    mean_direction = gp_mean > up_mean
    pairwise_support = sum(gp_value > up_value for gp_value in gp for up_value in up)

    leave_one_out_details: list[str] = []
    leave_one_out_support = 0
    leave_one_out_log2_ratios: list[float] = []
    for index, sample in enumerate(GP_SAMPLES):
        reduced_gp_mean = mean(gp[:index] + gp[index + 1 :])
        concordant = reduced_gp_mean > up_mean
        leave_one_out_support += int(concordant)
        leave_one_out_details.append(f"without_{SAMPLE_MAP[sample]}={truth(concordant)}")
        leave_one_out_log2_ratios.append(
            math.log2((reduced_gp_mean + PSEUDOCOUNT) / (up_mean + PSEUDOCOUNT))
        )
    for index, sample in enumerate(UP_SAMPLES):
        reduced_up_mean = mean(up[:index] + up[index + 1 :])
        concordant = gp_mean > reduced_up_mean
        leave_one_out_support += int(concordant)
        leave_one_out_details.append(f"without_{SAMPLE_MAP[sample]}={truth(concordant)}")
        leave_one_out_log2_ratios.append(
            math.log2((gp_mean + PSEUDOCOUNT) / (reduced_up_mean + PSEUDOCOUNT))
        )

    gp_dominant_sample, gp_dominant_fraction = dominance(gp, GP_SAMPLES)
    up_dominant_sample, up_dominant_fraction = dominance(up, UP_SAMPLES)
    gp_nonzero = sum(value > 0 for value in gp)
    up_nonzero = sum(value > 0 for value in up)
    complete_pairwise = pairwise_support == 9
    majority_pairwise = pairwise_support >= 5

    flags: list[str] = []
    if not mean_direction:
        flags.append("mean_direction_not_GP_up")
    if not majority_pairwise:
        flags.append("pairwise_support_below_5_of_9")
    elif not complete_pairwise:
        flags.append("incomplete_pairwise_separation")
    if leave_one_out_support < 6:
        flags.append("leave_one_sample_out_direction_failure")
    if gp_nonzero < 3:
        flags.append("GP_not_detected_in_all_replicates")
    if gp_dominant_fraction >= DOMINANCE_REVIEW_THRESHOLD:
        flags.append("GP_single_replicate_ge80pct_of_group_total")

    result: dict[str, object] = dict(candidate)
    result.update(
        {
            "GP1_normalized_count": gp[0],
            "GP2_normalized_count": gp[1],
            "GP3_normalized_count": gp[2],
            "UP1_normalized_count": up[0],
            "UP2_normalized_count": up[1],
            "UP3_normalized_count": up[2],
            "GP_mean_normalized_count": gp_mean,
            "GP_median_normalized_count": statistics.median(gp),
            "GP_coefficient_of_variation": fmt(coefficient_of_variation(gp)),
            "GP_nonzero_replicates": gp_nonzero,
            "UP_mean_normalized_count": up_mean,
            "UP_median_normalized_count": statistics.median(up),
            "UP_coefficient_of_variation": fmt(coefficient_of_variation(up)),
            "UP_nonzero_replicates": up_nonzero,
            "descriptive_log2_GP_UP_mean_ratio_pseudocount_0.5": math.log2(
                (gp_mean + PSEUDOCOUNT) / (up_mean + PSEUDOCOUNT)
            ),
            "mean_direction_concordant": truth(mean_direction),
            "GP_vs_UP_pairwise_support": pairwise_support,
            "GP_vs_UP_pairwise_total": 9,
            "GP_vs_UP_pairwise_support_fraction": pairwise_support / 9.0,
            "complete_pairwise_separation": truth(complete_pairwise),
            "majority_pairwise_support": truth(majority_pairwise),
            "leave_one_sample_out_concordant": leave_one_out_support,
            "leave_one_sample_out_total": 6,
            "all_leave_one_sample_out_concordant": truth(leave_one_out_support == 6),
            "leave_one_sample_out_minimum_descriptive_log2_mean_ratio": min(
                leave_one_out_log2_ratios
            ),
            "leave_one_sample_out_details": ";".join(leave_one_out_details),
            "largest_GP_replicate_fraction": gp_dominant_fraction,
            "dominant_GP_sample": gp_dominant_sample,
            "largest_UP_replicate_fraction": up_dominant_fraction,
            "dominant_UP_sample": up_dominant_sample,
            "replicate_review_flag": truth(bool(flags)),
            "replicate_review_reasons": ";".join(flags),
        }
    )
    return result


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def summarize(label: str, rows: list[dict[str, object]]) -> dict[str, object]:
    return {
        "candidate_set": label,
        "audited_candidates": len(rows),
        "mean_direction_concordant": sum(row["mean_direction_concordant"] == "TRUE" for row in rows),
        "complete_pairwise_separation": sum(row["complete_pairwise_separation"] == "TRUE" for row in rows),
        "majority_pairwise_support": sum(row["majority_pairwise_support"] == "TRUE" for row in rows),
        "all_leave_one_sample_out_concordant": sum(row["all_leave_one_sample_out_concordant"] == "TRUE" for row in rows),
        "GP_detected_in_all_replicates": sum(int(row["GP_nonzero_replicates"]) == 3 for row in rows),
        "GP_single_replicate_dominance_ge80pct": sum(float(row["largest_GP_replicate_fraction"]) >= DOMINANCE_REVIEW_THRESHOLD for row in rows),
        "candidates_with_any_review_flag": sum(row["replicate_review_flag"] == "TRUE" for row in rows),
    }


def main() -> int:
    required = [TOP50, POOR15, COUNTS]
    hashes_before = {str(path): sha256(path) for path in required}
    top_fields, top_rows = read_tsv(TOP50)
    poor_fields, poor_rows = read_tsv(POOR15)
    count_fields, count_rows = read_tsv(COUNTS)

    required_candidate = {
        "gene_id", "discovery_rank", "discovery_score", "padj", "log2FoldChange",
        "DESeq2_status_host_focused_expanded_mask",
    }
    require_columns(TOP50, top_fields, required_candidate)
    require_columns(POOR15, poor_fields, required_candidate | {"annotation_status"})
    require_columns(COUNTS, count_fields, {"gene_id", *SAMPLE_MAP.keys()})
    top = unique_index(TOP50, top_fields, top_rows)
    poor = unique_index(POOR15, poor_fields, poor_rows)
    counts = unique_index(COUNTS, count_fields, count_rows)

    if len(top_rows) != EXPECTED_TOP50_ROWS:
        raise AuditError(f"Top-50 rows: {len(top_rows)}; expected 50")
    if len(poor_rows) != EXPECTED_POOR15_ROWS:
        raise AuditError(f"Poorly annotated top-15 rows: {len(poor_rows)}; expected 15")
    if len(count_rows) != EXPECTED_COUNT_ROWS:
        raise AuditError(f"Expanded-mask normalized-count rows: {len(count_rows)}; expected 24931")
    if not set(poor).issubset(top):
        raise AuditError(
            f"Poorly annotated top-15 contains {len(set(poor)-set(top))} genes outside the top 50"
        )
    if [int(float(row["discovery_rank"])) for row in top_rows] != list(range(1, 51)):
        raise AuditError("Top-50 discovery ranks are not exactly 1 through 50")
    if any(
        clean(row.get("DESeq2_status_host_focused_expanded_mask")) != "GP_up"
        or float(row["padj"]) >= 0.01
        or float(row["log2FoldChange"]) < 1
        for row in [*top_rows, *poor_rows]
    ):
        raise AuditError("One or more candidate rows disagree with the frozen GP-up definition")
    if any(clean(row.get("annotation_status")) not in {"poorly_characterized", "uncharacterized"} for row in poor_rows):
        raise AuditError("Poorly annotated top-15 contains an informative candidate")

    missing = sorted((set(top) | set(poor)).difference(counts))
    if missing:
        raise AuditError(f"Candidates missing from normalized counts: {', '.join(missing[:20])}")

    audited_top = [audit_candidate(row, counts[row["gene_id"]]) for row in top_rows]
    audited_poor = [audit_candidate(row, counts[row["gene_id"]]) for row in poor_rows]
    added_fields = [field for field in audited_top[0] if field not in top_fields]
    write_tsv(TOP50_OUTPUT, audited_top, top_fields + added_fields)
    write_tsv(POOR15_OUTPUT, audited_poor, poor_fields + added_fields)

    mapping_rows = [
        {
            "source_count_column": source,
            "manuscript_sample_label": manuscript,
            "condition": "GP" if manuscript.startswith("GP") else "UP",
        }
        for source, manuscript in SAMPLE_MAP.items()
    ]
    write_tsv(
        MAPPING_OUTPUT,
        mapping_rows,
        ["source_count_column", "manuscript_sample_label", "condition"],
    )

    summaries = [
        summarize("plant_GP_host_focused_expanded_mask_top50", audited_top),
        summarize("plant_GP_host_focused_expanded_mask_poorly_annotated_top15", audited_poor),
    ]
    write_tsv(SUMMARY_OUTPUT, summaries, list(summaries[0]))

    if any(sha256(Path(path)) != digest for path, digest in hashes_before.items()):
        raise AuditError("One or more frozen input files changed during output generation")

    top_flag_reasons = Counter(
        reason
        for row in audited_top
        for reason in clean(row["replicate_review_reasons"]).split(";")
        if reason
    )
    lines = [
        "===== EXPANDED-MASK HOST-FOCUSED PLANT DISCOVERY REPLICATE ROBUSTNESS AUDIT =====",
        "Scope: descriptive replicate support for the frozen expanded-mask top 50 and poorly annotated top 15.",
        "No new hypothesis test, DEG reclassification, candidate removal or reranking was performed.",
        "",
        "Sample mapping:",
    ]
    lines.extend(f"  {source} -> {manuscript}" for source, manuscript in SAMPLE_MAP.items())
    lines.extend(
        [
            "",
            f"Expanded-mask normalized-count rows: {len(count_rows)}",
            f"Top-50 candidates audited: {len(audited_top)}",
            f"Poorly annotated top-15 candidates audited: {len(audited_poor)}",
            "Candidates missing from normalized counts: 0",
        ]
    )
    for summary in summaries:
        lines.extend(
            [
                "",
                f"{summary['candidate_set']}:",
                f"  Mean direction concordant: {summary['mean_direction_concordant']} / {summary['audited_candidates']}",
                f"  Complete GP-versus-UP pairwise separation: {summary['complete_pairwise_separation']} / {summary['audited_candidates']}",
                f"  Majority pairwise support: {summary['majority_pairwise_support']} / {summary['audited_candidates']}",
                f"  All six leave-one-sample-out comparisons concordant: {summary['all_leave_one_sample_out_concordant']} / {summary['audited_candidates']}",
                f"  Detected in all three GP replicates: {summary['GP_detected_in_all_replicates']} / {summary['audited_candidates']}",
                f"  GP single-replicate dominance >=80%: {summary['GP_single_replicate_dominance_ge80pct']} / {summary['audited_candidates']}",
                f"  Candidates with any review flag: {summary['candidates_with_any_review_flag']} / {summary['audited_candidates']}",
            ]
        )
    lines.extend(["", "Top-50 review-flag reasons:"])
    if top_flag_reasons:
        for reason, count in sorted(top_flag_reasons.items()):
            lines.append(f"  {reason}: {count}")
    else:
        lines.append("  none: 0")
    lines.extend(
        [
            "",
            "Interpretive boundary:",
            "  These metrics are descriptive checks based on DESeq2-normalized counts.",
            "  They are not independent significance tests and do not alter the frozen ranking.",
            "  A review flag requests closer inspection; it is not an automatic exclusion criterion.",
            "  UP single-replicate dominance is reported in gene tables but is not flagged because near-zero UP totals can yield misleading fractions.",
            "  The pseudocount of 0.5 is used only for descriptive log2 mean ratios.",
            "",
            "Frozen-input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "EXPANDED-MASK HOST-FOCUSED PLANT DISCOVERY REPLICATE ROBUSTNESS AUDIT: PASS",
        ]
    )
    AUDIT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT_OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    for path in [TOP50_OUTPUT, POOR15_OUTPUT, MAPPING_OUTPUT, SUMMARY_OUTPUT, AUDIT_OUTPUT]:
        print(f"Output: {path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(
            "===== EXPANDED-MASK HOST-FOCUSED PLANT DISCOVERY REPLICATE ROBUSTNESS AUDIT =====",
            file=sys.stderr,
        )
        print(f"ERROR: {error}", file=sys.stderr)
        print(
            "EXPANDED-MASK HOST-FOCUSED PLANT DISCOVERY REPLICATE ROBUSTNESS AUDIT: CHECK REQUIRED",
            file=sys.stderr,
        )
        raise SystemExit(1)
