#!/usr/bin/env python3
"""Replicate-level robustness audit for frozen insect and plant discovery lists.

This script adds descriptive replicate-support metrics to frozen candidate
tables. It does not test significance, filter candidates, or alter rankings.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import statistics
import sys
from pathlib import Path


HOME = Path.home()
ROOT = HOME / "projects/qe_cross_species_discovery"

INPUTS = {
    "insect_top50": ROOT / "02_insect_discovery/results/insect_larva_secretome_top50.tsv",
    "insect_top15": ROOT / "02_insect_discovery/results/insect_larva_secretome_poorly_annotated_top15.tsv",
    "plant_top50": ROOT / "03_plant_discovery/results/plant_GP_top50.tsv",
    "insect_counts": HOME / "projects/qe_genome_guided/08_results_primary/normalized_counts.tsv",
    "plant_counts": HOME / "projects/qe_plant_denovo/08_results_primary/plant_gene_normalized_counts_primary.tsv",
    "insect_metadata": HOME / "projects/qe_genome_guided/06_counts/sample_metadata.tsv",
    "plant_metadata": HOME / "projects/qe_plant_denovo/01_input/expression/plant_deseq2_samples.tsv",
}

EXPECTED_SHA256 = {
    "insect_top50": "b9dfa7c4b0f617c296c9b721eb505b4ebd6afed41f64e1ad92b2d16d16ad9b1d",
    "insect_top15": "999e95487f565dc3a38f4b9b5514403106d4b8119a4df11eadf35212437e6b8a",
    "plant_top50": "1d4aff7ecab264a97c5c7944a74919293c6f8514aded0f7b56e7cf8a672aa4bb",
    "insect_counts": "3cb87e132ab43caca58c819536e1be1986772917a9ceef406bb23178aa9f1ad2",
    "plant_counts": "05b0034c5f286bc4d9035e1978ccd177ad1943e6657962fcef87382e45aa5497",
    "insect_metadata": "9aa4da092a346a912b5df8b06c8a8e0e2e0d2335988f15f60aad81bcebd73ef3",
    "plant_metadata": "f9c74c1ef5c8a24fbb82daff221c6e07c739a28df561fa0d9a23d74620592ef0",
}

INSECT_RAW_SAMPLES = ["L1", "L3", "L4", "P3", "P4", "P5"]
INSECT_MANUSCRIPT_LABELS = {
    "L1": "L1",
    "L3": "L2",
    "L4": "L3",
    "P3": "P1",
    "P4": "P2",
    "P5": "P3",
}
INSECT_EXPECTED_GROUPS = {
    "L1": "larva",
    "L3": "larva",
    "L4": "larva",
    "P3": "pupa",
    "P4": "pupa",
    "P5": "pupa",
}

PLANT_RAW_SAMPLES = [
    "LTS20_VA01_Sample_S7",
    "LTS20_VA01_Sample_S8",
    "LTS20_VA01_Sample_S9",
    "LTS20_VA01_Sample_S10",
    "LTS20_VA01_Sample_S11",
    "LTS20_VA01_Sample_S12",
]
PLANT_MANUSCRIPT_LABELS = {
    "LTS20_VA01_Sample_S7": "GP1",
    "LTS20_VA01_Sample_S8": "GP2",
    "LTS20_VA01_Sample_S9": "GP3",
    "LTS20_VA01_Sample_S10": "UP1",
    "LTS20_VA01_Sample_S11": "UP2",
    "LTS20_VA01_Sample_S12": "UP3",
}
PLANT_EXPECTED_GROUPS = {
    "LTS20_VA01_Sample_S7": "SG",
    "LTS20_VA01_Sample_S8": "SG",
    "LTS20_VA01_Sample_S9": "SG",
    "LTS20_VA01_Sample_S10": "S",
    "LTS20_VA01_Sample_S11": "S",
    "LTS20_VA01_Sample_S12": "S",
}

OUTPUTS = {
    "insect_top50": ROOT / "02_insect_discovery/results/insect_larva_secretome_top50_replicate_audit.tsv",
    "insect_top15": ROOT / "02_insect_discovery/results/insect_larva_secretome_poorly_annotated_top15_replicate_audit.tsv",
    "plant_top50": ROOT / "03_plant_discovery/results/plant_GP_top50_replicate_audit.tsv",
    "sample_mapping": ROOT / "00_project_control/discovery_replicate_sample_mapping.tsv",
    "summary": ROOT / "00_project_control/03_discovery_replicate_robustness_summary.tsv",
    "audit": ROOT / "00_project_control/03_discovery_replicate_robustness_audit.txt",
}

PSEUDOCOUNT = 0.5
SINGLE_REPLICATE_SHARE_REVIEW_THRESHOLD = 0.80


class AuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditError(f"Missing TSV header: {path}")
        rows = list(reader)
        return list(reader.fieldnames), rows


def require_unique(rows: list[dict[str, str]], key: str, label: str) -> None:
    values = [row.get(key, "") for row in rows]
    if any(not value for value in values):
        raise AuditError(f"Blank {key} in {label}")
    if len(values) != len(set(values)):
        raise AuditError(f"Duplicate {key} in {label}")


def number(value: str, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise AuditError(f"Non-numeric value for {label}: {value!r}") from exc
    if not math.isfinite(result):
        raise AuditError(f"Non-finite value for {label}: {value!r}")
    return result


def verify_inputs() -> dict[str, str]:
    observed: dict[str, str] = {}
    for name, path in INPUTS.items():
        if not path.is_file():
            raise AuditError(f"Missing required input: {path}")
        observed[name] = sha256(path)
        if observed[name] != EXPECTED_SHA256[name]:
            raise AuditError(
                f"SHA-256 mismatch for {name}: expected {EXPECTED_SHA256[name]}, "
                f"observed {observed[name]}"
            )
    return observed


def verify_metadata() -> None:
    insect_fields, insect_rows = read_tsv(INPUTS["insect_metadata"])
    if insect_fields != ["sample", "stage"]:
        raise AuditError(f"Unexpected insect metadata header: {insect_fields}")
    require_unique(insect_rows, "sample", "insect metadata")
    insect_groups = {row["sample"]: row["stage"] for row in insect_rows}
    if insect_groups != INSECT_EXPECTED_GROUPS:
        raise AuditError(f"Unexpected insect sample-to-stage mapping: {insect_groups}")

    plant_fields, plant_rows = read_tsv(INPUTS["plant_metadata"])
    if plant_fields != ["sample_id", "condition"]:
        raise AuditError(f"Unexpected plant metadata header: {plant_fields}")
    require_unique(plant_rows, "sample_id", "plant metadata")
    plant_groups = {row["sample_id"]: row["condition"] for row in plant_rows}
    if plant_groups != PLANT_EXPECTED_GROUPS:
        raise AuditError(f"Unexpected plant sample-to-condition mapping: {plant_groups}")


def read_count_matrix(
    path: Path,
    id_column: str,
    expected_samples: list[str],
) -> dict[str, dict[str, float]]:
    fields, rows = read_tsv(path)
    expected_header = [id_column, *expected_samples]
    if fields != expected_header:
        raise AuditError(
            f"Unexpected count-matrix header in {path}: expected {expected_header}, observed {fields}"
        )
    require_unique(rows, id_column, str(path))
    result: dict[str, dict[str, float]] = {}
    for row in rows:
        gene = row[id_column]
        values: dict[str, float] = {}
        for sample in expected_samples:
            value = number(row[sample], f"{gene}/{sample}")
            if value < 0:
                raise AuditError(f"Negative normalized count for {gene}/{sample}")
            values[sample] = value
        result[gene] = values
    return result


def fmt(value: float) -> str:
    return format(value, ".12g")


def truth(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def mean(values: list[float]) -> float:
    return statistics.fmean(values)


def group_cv(values: list[float]) -> float | None:
    group_mean = mean(values)
    if group_mean == 0:
        return None
    return statistics.stdev(values) / group_mean


def directional_metrics(
    focal: list[float],
    reference: list[float],
    focal_name: str,
    reference_name: str,
) -> dict[str, str]:
    focal_mean = mean(focal)
    reference_mean = mean(reference)
    focal_median = statistics.median(focal)
    reference_median = statistics.median(reference)
    pairwise_count = sum(a > b for a in focal for b in reference)
    complete_separation = pairwise_count == len(focal) * len(reference)
    mean_concordant = focal_mean > reference_mean

    leave_one_out_results: list[bool] = []
    leave_one_out_ratios: list[float] = []
    for index in range(len(focal)):
        reduced_focal = focal[:index] + focal[index + 1 :]
        reduced_focal_mean = mean(reduced_focal)
        leave_one_out_results.append(reduced_focal_mean > reference_mean)
        leave_one_out_ratios.append(
            math.log2((reduced_focal_mean + PSEUDOCOUNT) / (reference_mean + PSEUDOCOUNT))
        )
    for index in range(len(reference)):
        reduced_reference = reference[:index] + reference[index + 1 :]
        reduced_reference_mean = mean(reduced_reference)
        leave_one_out_results.append(focal_mean > reduced_reference_mean)
        leave_one_out_ratios.append(
            math.log2((focal_mean + PSEUDOCOUNT) / (reduced_reference_mean + PSEUDOCOUNT))
        )

    focal_total = sum(focal)
    reference_total = sum(reference)
    focal_max_share = max(focal) / focal_total if focal_total > 0 else None
    reference_max_share = max(reference) / reference_total if reference_total > 0 else None
    focal_cv = group_cv(focal)
    reference_cv = group_cv(reference)
    focal_nonzero = sum(value > 0 for value in focal)
    reference_nonzero = sum(value > 0 for value in reference)

    if complete_separation:
        support_class = "complete_pairwise_separation"
    elif pairwise_count >= 5:
        support_class = "majority_pairwise_support"
    else:
        support_class = "limited_or_reverse_pairwise_support"

    flags: list[str] = []
    if not mean_concordant:
        flags.append("mean_direction_not_concordant")
    if not all(leave_one_out_results):
        flags.append("leave_one_out_direction_not_fully_concordant")
    if focal_nonzero < len(focal):
        flags.append(f"{focal_name}_not_detected_in_all_replicates")
    if focal_max_share is not None and focal_max_share >= SINGLE_REPLICATE_SHARE_REVIEW_THRESHOLD:
        flags.append(f"{focal_name}_single_replicate_ge80pct_of_group_total")

    def optional(value: float | None) -> str:
        return "" if value is None else fmt(value)

    prefix_f = focal_name
    prefix_r = reference_name
    return {
        f"{prefix_f}_mean_normalized_count": fmt(focal_mean),
        f"{prefix_f}_median_normalized_count": fmt(focal_median),
        f"{prefix_f}_minimum_normalized_count": fmt(min(focal)),
        f"{prefix_f}_maximum_normalized_count": fmt(max(focal)),
        f"{prefix_f}_coefficient_of_variation": optional(focal_cv),
        f"{prefix_f}_nonzero_replicates": str(focal_nonzero),
        f"{prefix_f}_maximum_fraction_of_group_total": optional(focal_max_share),
        f"{prefix_r}_mean_normalized_count": fmt(reference_mean),
        f"{prefix_r}_median_normalized_count": fmt(reference_median),
        f"{prefix_r}_minimum_normalized_count": fmt(min(reference)),
        f"{prefix_r}_maximum_normalized_count": fmt(max(reference)),
        f"{prefix_r}_coefficient_of_variation": optional(reference_cv),
        f"{prefix_r}_nonzero_replicates": str(reference_nonzero),
        f"{prefix_r}_maximum_fraction_of_group_total": optional(reference_max_share),
        f"descriptive_log2_mean_ratio_{prefix_f}_over_{prefix_r}_pseudocount_0.5": fmt(
            math.log2((focal_mean + PSEUDOCOUNT) / (reference_mean + PSEUDOCOUNT))
        ),
        f"pairwise_{prefix_f}_greater_than_{prefix_r}_count_out_of_9": str(pairwise_count),
        f"pairwise_{prefix_f}_greater_than_{prefix_r}_fraction": fmt(pairwise_count / 9),
        "mean_direction_concordant_with_DESeq2": truth(mean_concordant),
        "complete_pairwise_separation": truth(complete_separation),
        "leave_one_sample_out_direction_concordant_count_out_of_6": str(sum(leave_one_out_results)),
        "leave_one_sample_out_all_direction_concordant": truth(all(leave_one_out_results)),
        "leave_one_sample_out_minimum_descriptive_log2_mean_ratio": fmt(min(leave_one_out_ratios)),
        "replicate_support_class": support_class,
        "replicate_review_status": "review" if flags else "no_major_review_flag",
        "replicate_review_flags": ";".join(flags) if flags else "none",
    }


def validate_candidate_table(
    fields: list[str],
    rows: list[dict[str, str]],
    expected_rows: int,
    organism: str,
    table_name: str,
) -> None:
    if "gene_id" not in fields:
        raise AuditError(f"Missing gene_id column in {table_name}")
    require_unique(rows, "gene_id", table_name)
    if len(rows) != expected_rows:
        raise AuditError(f"Unexpected row count in {table_name}: {len(rows)}; expected {expected_rows}")

    previous_score = math.inf
    for row in rows:
        score = number(row["discovery_score"], f"{table_name}/discovery_score")
        if score > previous_score + 1e-12:
            raise AuditError(f"Discovery-score order is not descending in {table_name}")
        previous_score = score
        padj = number(row["padj"], f"{table_name}/padj")
        if not padj < 0.01:
            raise AuditError(f"Non-official DEG in {table_name}: {row['gene_id']}")
        if organism == "insect":
            lfc = number(row["log2FoldChange_pupa_vs_larva"], f"{table_name}/log2FC")
            if lfc > -1:
                raise AuditError(f"Non-larva-up row in {table_name}: {row['gene_id']}")
            if row.get("DESeq2_status") not in {"", "larva_up"}:
                raise AuditError(f"Unexpected DESeq2 status in {table_name}: {row['gene_id']}")
        else:
            lfc = number(row["log2FoldChange_GP_vs_UP"], f"{table_name}/log2FC")
            if lfc < 1:
                raise AuditError(f"Non-GP-up row in {table_name}: {row['gene_id']}")
            if row.get("DESeq2_status") not in {"", "GP_up"}:
                raise AuditError(f"Unexpected DESeq2 status in {table_name}: {row['gene_id']}")


def annotate_table(
    fields: list[str],
    rows: list[dict[str, str]],
    counts: dict[str, dict[str, float]],
    raw_samples: list[str],
    manuscript_labels: dict[str, str],
    focal_samples: list[str],
    reference_samples: list[str],
    focal_name: str,
    reference_name: str,
) -> tuple[list[str], list[dict[str, str]]]:
    missing = [row["gene_id"] for row in rows if row["gene_id"] not in counts]
    if missing:
        raise AuditError(f"Candidate genes missing from normalized counts: {', '.join(missing)}")

    output_rows: list[dict[str, str]] = []
    added_fields = [f"normalized_count_{manuscript_labels[sample]}" for sample in raw_samples]
    metric_fields: list[str] | None = None
    for row in rows:
        gene = row["gene_id"]
        new_row = dict(row)
        for sample in raw_samples:
            new_row[f"normalized_count_{manuscript_labels[sample]}"] = fmt(counts[gene][sample])
        metrics = directional_metrics(
            [counts[gene][sample] for sample in focal_samples],
            [counts[gene][sample] for sample in reference_samples],
            focal_name,
            reference_name,
        )
        if metric_fields is None:
            metric_fields = list(metrics)
        new_row.update(metrics)
        output_rows.append(new_row)

    return [*fields, *added_fields, *(metric_fields or [])], output_rows


def write_tsv_atomic(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def summarize(candidate_set: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    metrics = {
        "audited_candidates": len(rows),
        "mean_direction_concordant": sum(row["mean_direction_concordant_with_DESeq2"] == "TRUE" for row in rows),
        "complete_pairwise_separation": sum(row["complete_pairwise_separation"] == "TRUE" for row in rows),
        "majority_or_complete_pairwise_support": sum(
            row["replicate_support_class"] != "limited_or_reverse_pairwise_support" for row in rows
        ),
        "all_six_leave_one_sample_out_comparisons_concordant": sum(
            row["leave_one_sample_out_all_direction_concordant"] == "TRUE" for row in rows
        ),
        "no_major_review_flag": sum(row["replicate_review_status"] == "no_major_review_flag" for row in rows),
        "review_flagged": sum(row["replicate_review_status"] == "review" for row in rows),
    }
    return [
        {"candidate_set": candidate_set, "metric": metric, "value": str(value)}
        for metric, value in metrics.items()
    ]


def main() -> None:
    observed_hashes = verify_inputs()
    verify_metadata()

    insect_top50_fields, insect_top50_rows = read_tsv(INPUTS["insect_top50"])
    insect_top15_fields, insect_top15_rows = read_tsv(INPUTS["insect_top15"])
    plant_top50_fields, plant_top50_rows = read_tsv(INPUTS["plant_top50"])

    validate_candidate_table(insect_top50_fields, insect_top50_rows, 50, "insect", "insect top 50")
    validate_candidate_table(insect_top15_fields, insect_top15_rows, 15, "insect", "insect poorly annotated top 15")
    validate_candidate_table(plant_top50_fields, plant_top50_rows, 50, "plant", "plant top 50")

    expected_insect_ranks = list(range(1, 51))
    observed_insect_ranks = [int(row["discovery_rank"]) for row in insect_top50_rows]
    if observed_insect_ranks != expected_insect_ranks:
        raise AuditError("Insect top-50 discovery ranks are not exactly 1 through 50")
    observed_plant_ranks = [int(row["discovery_rank"]) for row in plant_top50_rows]
    if observed_plant_ranks != expected_insect_ranks:
        raise AuditError("Plant top-50 discovery ranks are not exactly 1 through 50")

    insect_counts = read_count_matrix(INPUTS["insect_counts"], "Geneid", INSECT_RAW_SAMPLES)
    plant_counts = read_count_matrix(INPUTS["plant_counts"], "gene_id", PLANT_RAW_SAMPLES)

    insect_top50_out_fields, insect_top50_out = annotate_table(
        insect_top50_fields,
        insect_top50_rows,
        insect_counts,
        INSECT_RAW_SAMPLES,
        INSECT_MANUSCRIPT_LABELS,
        ["L1", "L3", "L4"],
        ["P3", "P4", "P5"],
        "larva",
        "pupa",
    )
    insect_top15_out_fields, insect_top15_out = annotate_table(
        insect_top15_fields,
        insect_top15_rows,
        insect_counts,
        INSECT_RAW_SAMPLES,
        INSECT_MANUSCRIPT_LABELS,
        ["L1", "L3", "L4"],
        ["P3", "P4", "P5"],
        "larva",
        "pupa",
    )
    plant_top50_out_fields, plant_top50_out = annotate_table(
        plant_top50_fields,
        plant_top50_rows,
        plant_counts,
        PLANT_RAW_SAMPLES,
        PLANT_MANUSCRIPT_LABELS,
        PLANT_RAW_SAMPLES[:3],
        PLANT_RAW_SAMPLES[3:],
        "GP",
        "UP",
    )

    mapping_rows = []
    for raw in INSECT_RAW_SAMPLES:
        mapping_rows.append(
            {
                "organism": "Q. erythrinae",
                "count_matrix_sample_id": raw,
                "metadata_group": INSECT_EXPECTED_GROUPS[raw],
                "manuscript_sample_label": INSECT_MANUSCRIPT_LABELS[raw],
            }
        )
    for raw in PLANT_RAW_SAMPLES:
        mapping_rows.append(
            {
                "organism": "E. variegata",
                "count_matrix_sample_id": raw,
                "metadata_group": PLANT_EXPECTED_GROUPS[raw],
                "manuscript_sample_label": PLANT_MANUSCRIPT_LABELS[raw],
            }
        )

    summary_rows = []
    summary_rows.extend(summarize("insect_larva_secretome_top50", insect_top50_out))
    summary_rows.extend(summarize("insect_larva_secretome_poorly_annotated_top15", insect_top15_out))
    summary_rows.extend(summarize("plant_GP_top50", plant_top50_out))

    write_tsv_atomic(OUTPUTS["insect_top50"], insect_top50_out_fields, insect_top50_out)
    write_tsv_atomic(OUTPUTS["insect_top15"], insect_top15_out_fields, insect_top15_out)
    write_tsv_atomic(OUTPUTS["plant_top50"], plant_top50_out_fields, plant_top50_out)
    write_tsv_atomic(
        OUTPUTS["sample_mapping"],
        ["organism", "count_matrix_sample_id", "metadata_group", "manuscript_sample_label"],
        mapping_rows,
    )
    write_tsv_atomic(OUTPUTS["summary"], ["candidate_set", "metric", "value"], summary_rows)

    # Confirm that the frozen ranked inputs remain byte-for-byte unchanged.
    for name in ("insect_top50", "insect_top15", "plant_top50"):
        if sha256(INPUTS[name]) != observed_hashes[name]:
            raise AuditError(f"Frozen ranked input changed during audit: {INPUTS[name]}")

    summary_lookup = {
        (row["candidate_set"], row["metric"]): row["value"] for row in summary_rows
    }
    report = [
        "===== DISCOVERY REPLICATE ROBUSTNESS AUDIT =====",
        f"Discovery root: {ROOT}",
        "",
        "Authoritative sample-label mappings:",
        "  Insect count IDs: L1->L1, L3->L2, L4->L3, P3->P1, P4->P2, P5->P3",
        "  Plant count IDs: S7->GP1, S8->GP2, S9->GP3, S10->UP1, S11->UP2, S12->UP3",
        "  Plant metadata terminology: SG=GP and S=UP",
        "",
        "Interpretive scope:",
        "  Replicate metrics are descriptive sensitivity checks based on DESeq2-normalized counts.",
        "  They do not constitute new significance tests, alter the official DEG calls, remove genes, or rerank candidates.",
        "  The descriptive mean-ratio calculation uses a pseudocount of 0.5 only to avoid undefined logarithms.",
        "  A single-replicate share >=80% is flagged for review but is not an exclusion criterion.",
        "",
    ]
    for candidate_set in (
        "insect_larva_secretome_top50",
        "insect_larva_secretome_poorly_annotated_top15",
        "plant_GP_top50",
    ):
        report.extend(
            [
                f"{candidate_set}:",
                f"  Audited candidates: {summary_lookup[(candidate_set, 'audited_candidates')]}",
                f"  Mean direction concordant: {summary_lookup[(candidate_set, 'mean_direction_concordant')]}",
                f"  Complete pairwise separation: {summary_lookup[(candidate_set, 'complete_pairwise_separation')]}",
                f"  Majority or complete pairwise support: {summary_lookup[(candidate_set, 'majority_or_complete_pairwise_support')]}",
                f"  All six leave-one-sample-out comparisons concordant: {summary_lookup[(candidate_set, 'all_six_leave_one_sample_out_comparisons_concordant')]}",
                f"  Review flagged: {summary_lookup[(candidate_set, 'review_flagged')]}",
                "",
            ]
        )
    report.extend(
        [
            "Frozen ranked-input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "",
            f"Output: {OUTPUTS['insect_top50']}",
            f"Output: {OUTPUTS['insect_top15']}",
            f"Output: {OUTPUTS['plant_top50']}",
            f"Output: {OUTPUTS['sample_mapping']}",
            f"Output: {OUTPUTS['summary']}",
            f"Output: {OUTPUTS['audit']}",
            "DISCOVERY REPLICATE ROBUSTNESS AUDIT: PASS",
        ]
    )
    OUTPUTS["audit"].parent.mkdir(parents=True, exist_ok=True)
    audit_tmp = OUTPUTS["audit"].with_name(OUTPUTS["audit"].name + ".tmp")
    audit_tmp.write_text("\n".join(report) + "\n", encoding="utf-8")
    os.replace(audit_tmp, OUTPUTS["audit"])
    print("\n".join(report))


if __name__ == "__main__":
    try:
        main()
    except AuditError as error:
        print("===== DISCOVERY REPLICATE ROBUSTNESS AUDIT =====", file=sys.stderr)
        print(f"ERROR: {error}", file=sys.stderr)
        print("DISCOVERY REPLICATE ROBUSTNESS AUDIT: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
