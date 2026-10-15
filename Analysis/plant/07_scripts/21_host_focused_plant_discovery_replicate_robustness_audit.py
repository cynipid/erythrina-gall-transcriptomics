#!/usr/bin/env python3
"""Descriptive replicate-robustness audit for frozen host-focused plant lists.

The audit uses DESeq2-normalized counts only. It performs no new hypothesis
test, changes no DEG call, removes no candidate and does not rerank either
frozen discovery list.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from pathlib import Path


EXPECTED_COUNT_ROWS = 25_779
EXPECTED_TOP50_ROWS = 50
EXPECTED_REVIEW_ROWS = 15
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
GP_SOURCE_SAMPLES = list(SAMPLE_MAP)[:3]
UP_SOURCE_SAMPLES = list(SAMPLE_MAP)[3:]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "na", "nan", "none", "null"} else text


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
    path: Path, fields: list[str], rows: list[dict[str, str]], key: str
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


def number(value: object, context: str) -> float:
    text = clean(value)
    if not text:
        raise ValueError(f"Missing numeric value: {context}")
    try:
        result = float(text)
    except ValueError as exc:
        raise ValueError(f"Invalid numeric value for {context}: {text}") from exc
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"Invalid normalized count for {context}: {text}")
    return result


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def dominance(values: list[float], source_samples: list[str]) -> tuple[str, float]:
    total = sum(values)
    maximum_index = max(range(len(values)), key=lambda index: values[index])
    label = SAMPLE_MAP[source_samples[maximum_index]]
    fraction = values[maximum_index] / total if total > 0 else 0.0
    return label, fraction


def bool_text(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def audit_candidate(
    candidate: dict[str, str], count_row: dict[str, str]
) -> dict[str, object]:
    gene_id = candidate["gene_id"]
    gp = [number(count_row[sample], f"{gene_id}/{sample}") for sample in GP_SOURCE_SAMPLES]
    up = [number(count_row[sample], f"{gene_id}/{sample}") for sample in UP_SOURCE_SAMPLES]

    gp_mean = mean(gp)
    up_mean = mean(up)
    mean_direction = gp_mean > up_mean
    log2_ratio = math.log2((gp_mean + PSEUDOCOUNT) / (up_mean + PSEUDOCOUNT))

    pairwise_support = sum(gp_value > up_value for gp_value in gp for up_value in up)
    complete_pairwise = pairwise_support == 9
    majority_pairwise = pairwise_support >= 5

    leave_one_out_labels: list[str] = []
    leave_one_out_support = 0
    for index, source_sample in enumerate(GP_SOURCE_SAMPLES):
        retained_gp = gp[:index] + gp[index + 1 :]
        concordant = mean(retained_gp) > up_mean
        leave_one_out_support += int(concordant)
        leave_one_out_labels.append(f"without_{SAMPLE_MAP[source_sample]}={bool_text(concordant)}")
    for index, source_sample in enumerate(UP_SOURCE_SAMPLES):
        retained_up = up[:index] + up[index + 1 :]
        concordant = gp_mean > mean(retained_up)
        leave_one_out_support += int(concordant)
        leave_one_out_labels.append(f"without_{SAMPLE_MAP[source_sample]}={bool_text(concordant)}")

    gp_dominant_sample, gp_dominant_fraction = dominance(gp, GP_SOURCE_SAMPLES)
    up_dominant_sample, up_dominant_fraction = dominance(up, UP_SOURCE_SAMPLES)

    flags: list[str] = []
    if not mean_direction:
        flags.append("mean_direction_not_GP_up")
    if not majority_pairwise:
        flags.append("pairwise_support_below_5_of_9")
    elif not complete_pairwise:
        flags.append("incomplete_pairwise_separation")
    if leave_one_out_support < 6:
        flags.append("leave_one_sample_out_direction_failure")
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
            "UP_mean_normalized_count": up_mean,
            "descriptive_log2_GP_UP_mean_ratio_pseudocount_0.5": log2_ratio,
            "mean_direction_concordant": bool_text(mean_direction),
            "GP_vs_UP_pairwise_support": pairwise_support,
            "GP_vs_UP_pairwise_total": 9,
            "GP_vs_UP_pairwise_support_fraction": pairwise_support / 9.0,
            "complete_pairwise_separation": bool_text(complete_pairwise),
            "majority_pairwise_support": bool_text(majority_pairwise),
            "leave_one_sample_out_concordant": leave_one_out_support,
            "leave_one_sample_out_total": 6,
            "all_leave_one_sample_out_concordant": bool_text(
                leave_one_out_support == 6
            ),
            "leave_one_sample_out_details": ";".join(leave_one_out_labels),
            "largest_GP_replicate_fraction": gp_dominant_fraction,
            "dominant_GP_sample": gp_dominant_sample,
            "largest_UP_replicate_fraction": up_dominant_fraction,
            "dominant_UP_sample": up_dominant_sample,
            "replicate_review_flag": bool_text(bool(flags)),
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
        "mean_direction_concordant": sum(
            row["mean_direction_concordant"] == "TRUE" for row in rows
        ),
        "complete_pairwise_separation": sum(
            row["complete_pairwise_separation"] == "TRUE" for row in rows
        ),
        "majority_pairwise_support": sum(
            row["majority_pairwise_support"] == "TRUE" for row in rows
        ),
        "all_leave_one_sample_out_concordant": sum(
            row["all_leave_one_sample_out_concordant"] == "TRUE" for row in rows
        ),
        "review_flagged": sum(row["replicate_review_flag"] == "TRUE" for row in rows),
        "GP_dominance_ge80pct": sum(
            float(row["largest_GP_replicate_fraction"]) >= DOMINANCE_REVIEW_THRESHOLD
            for row in rows
        ),
        "UP_dominance_ge80pct": sum(
            float(row["largest_UP_replicate_fraction"]) >= DOMINANCE_REVIEW_THRESHOLD
            and sum(float(row[f"UP{i}_normalized_count"]) for i in range(1, 4)) > 0
            for row in rows
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--discovery-root",
        type=Path,
        default=Path.home() / "projects" / "qe_cross_species_discovery",
    )
    args = parser.parse_args()
    root = args.discovery_root.resolve()

    results_dir = root / "05_host_focused_DESeq2_sensitivity/results"
    control_dir = root / "00_project_control"
    logs_dir = root / "05_host_focused_DESeq2_sensitivity/logs"

    top50_input = results_dir / "plant_GP_host_focused_top50.tsv"
    review15_input = results_dir / "plant_GP_host_focused_poorly_annotated_top15.tsv"
    counts_input = results_dir / "plant_gene_normalized_counts_host_focused.tsv"

    top50_output = results_dir / "plant_GP_host_focused_top50_replicate_audit.tsv"
    review15_output = (
        results_dir
        / "plant_GP_host_focused_poorly_annotated_top15_replicate_audit.tsv"
    )
    mapping_output = control_dir / "21_host_focused_discovery_sample_mapping.tsv"
    summary_output = control_dir / (
        "21_host_focused_plant_discovery_replicate_robustness_summary.tsv"
    )
    audit_output = logs_dir / (
        "21_host_focused_plant_discovery_replicate_robustness_audit.txt"
    )

    top_fields, top_rows = read_tsv(top50_input)
    review_fields, review_rows = read_tsv(review15_input)
    count_fields, count_rows = read_tsv(counts_input)

    require_columns(top50_input, top_fields, {"gene_id", "discovery_rank"})
    require_columns(review15_input, review_fields, {"gene_id", "discovery_rank"})
    require_columns(counts_input, count_fields, {"gene_id", *SAMPLE_MAP.keys()})

    top_index = unique_index(top50_input, top_fields, top_rows, "gene_id")
    review_index = unique_index(review15_input, review_fields, review_rows, "gene_id")
    count_index = unique_index(counts_input, count_fields, count_rows, "gene_id")

    if len(top_rows) != EXPECTED_TOP50_ROWS:
        raise ValueError(f"Top-50 table has {len(top_rows)} rows; expected 50")
    if len(review_rows) != EXPECTED_REVIEW_ROWS:
        raise ValueError(f"Review table has {len(review_rows)} rows; expected 15")
    if len(count_rows) != EXPECTED_COUNT_ROWS:
        raise ValueError(
            f"Host-focused normalized-count table has {len(count_rows)} rows; "
            f"expected {EXPECTED_COUNT_ROWS}"
        )
    if not set(review_index).issubset(top_index):
        outside = sorted(set(review_index).difference(top_index))
        raise ValueError(
            "Poorly annotated top-15 is not a subset of the frozen top 50: "
            + ", ".join(outside)
        )

    missing_top = sorted(set(top_index).difference(count_index))
    missing_review = sorted(set(review_index).difference(count_index))
    if missing_top or missing_review:
        raise ValueError(
            f"Candidates missing from normalized counts: top50={len(missing_top)}, "
            f"review15={len(missing_review)}"
        )

    input_hashes = {
        str(top50_input): sha256(top50_input),
        str(review15_input): sha256(review15_input),
        str(counts_input): sha256(counts_input),
    }

    audited_top = [audit_candidate(row, count_index[row["gene_id"]]) for row in top_rows]
    audited_review = [
        audit_candidate(row, count_index[row["gene_id"]]) for row in review_rows
    ]

    added_fields = [
        "GP1_normalized_count",
        "GP2_normalized_count",
        "GP3_normalized_count",
        "UP1_normalized_count",
        "UP2_normalized_count",
        "UP3_normalized_count",
        "GP_mean_normalized_count",
        "UP_mean_normalized_count",
        "descriptive_log2_GP_UP_mean_ratio_pseudocount_0.5",
        "mean_direction_concordant",
        "GP_vs_UP_pairwise_support",
        "GP_vs_UP_pairwise_total",
        "GP_vs_UP_pairwise_support_fraction",
        "complete_pairwise_separation",
        "majority_pairwise_support",
        "leave_one_sample_out_concordant",
        "leave_one_sample_out_total",
        "all_leave_one_sample_out_concordant",
        "leave_one_sample_out_details",
        "largest_GP_replicate_fraction",
        "dominant_GP_sample",
        "largest_UP_replicate_fraction",
        "dominant_UP_sample",
        "replicate_review_flag",
        "replicate_review_reasons",
    ]
    write_tsv(top50_output, audited_top, top_fields + added_fields)
    write_tsv(review15_output, audited_review, review_fields + added_fields)

    mapping_rows = [
        {
            "source_count_column": source,
            "manuscript_sample_label": manuscript,
            "condition": "GP" if manuscript.startswith("GP") else "UP",
        }
        for source, manuscript in SAMPLE_MAP.items()
    ]
    write_tsv(
        mapping_output,
        mapping_rows,
        ["source_count_column", "manuscript_sample_label", "condition"],
    )

    summaries = [
        summarize("plant_GP_host_focused_top50", audited_top),
        summarize(
            "plant_GP_host_focused_poorly_annotated_top15", audited_review
        ),
    ]
    summary_fields = list(summaries[0])
    write_tsv(summary_output, summaries, summary_fields)

    if any(sha256(Path(path)) != digest for path, digest in input_hashes.items()):
        raise ValueError("One or more frozen input checksums changed during output generation")

    lines = [
        "===== HOST-FOCUSED PLANT DISCOVERY REPLICATE ROBUSTNESS AUDIT =====",
        f"Discovery root: {root}",
        "",
        "Authoritative sample-label mapping:",
    ]
    for source, manuscript in SAMPLE_MAP.items():
        lines.append(f"  {source} -> {manuscript}")
    lines.extend(
        [
            "",
            "Interpretive scope:",
            "  Metrics are descriptive sensitivity checks based on host-focused DESeq2-normalized counts.",
            "  They are not new significance tests and do not alter DEG calls, candidate membership or discovery ranks.",
            "  The descriptive mean ratio uses a pseudocount of 0.5 only to avoid undefined logarithms.",
            "  A GP single-replicate share >=80% is a review flag, not an exclusion criterion.",
            "  UP dominance is reported descriptively but is not flagged because near-zero background counts can give misleading fractions.",
            "",
            f"Host-focused normalized-count rows: {len(count_rows)}",
            f"Frozen top-50 rows: {len(top_rows)}",
            f"Frozen poorly annotated/uncharacterized top-15 rows: {len(review_rows)}",
            "Missing candidate count rows: 0",
        ]
    )
    for summary in summaries:
        lines.extend(
            [
                "",
                f"{summary['candidate_set']}:",
                f"  Audited candidates: {summary['audited_candidates']}",
                f"  Mean direction concordant: {summary['mean_direction_concordant']}",
                f"  Complete pairwise separation: {summary['complete_pairwise_separation']}",
                f"  Majority pairwise support: {summary['majority_pairwise_support']}",
                "  All six leave-one-sample-out comparisons concordant: "
                f"{summary['all_leave_one_sample_out_concordant']}",
                f"  GP single-replicate dominance flags: {summary['GP_dominance_ge80pct']}",
                "  UP single-replicate dominance descriptions (not review flags): "
                f"{summary['UP_dominance_ge80pct']}",
                f"  Candidates with any review flag: {summary['review_flagged']}",
            ]
        )
    lines.extend(
        [
            "",
            "Frozen-input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "",
            f"Output: {top50_output}",
            f"Output: {review15_output}",
            f"Output: {mapping_output}",
            f"Output: {summary_output}",
            f"Output: {audit_output}",
            "HOST-FOCUSED PLANT DISCOVERY REPLICATE ROBUSTNESS AUDIT: PASS",
        ]
    )
    logs_dir.mkdir(parents=True, exist_ok=True)
    audit_output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(
            "===== HOST-FOCUSED PLANT DISCOVERY REPLICATE ROBUSTNESS AUDIT =====",
            file=sys.stderr,
        )
        print(f"ERROR: {error}", file=sys.stderr)
        print(
            "HOST-FOCUSED PLANT DISCOVERY REPLICATE ROBUSTNESS AUDIT: CHECK REQUIRED",
            file=sys.stderr,
        )
        raise SystemExit(1)
