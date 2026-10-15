#!/usr/bin/env python3
"""Describe confirmed insect-like sequence abundance in GP and UP libraries.

This is an evidence-only audit of DESeq2-normalized counts. It neither
reruns DESeq2 nor changes the confirmed multilevel sequence mask.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import sys
from collections import Counter
from pathlib import Path


HOME = Path.home()
ROOT = HOME / "projects/qe_cross_species_discovery"
INPUTS = {
    "confirmed_mask": ROOT
    / "04_cross_species_integration/results/plant_confirmed_multilevel_insect_like_sequence_mask.tsv",
    "normalized_counts": HOME
    / "projects/qe_plant_denovo/08_results_primary/plant_gene_normalized_counts_primary.tsv",
}
OUTPUTS = {
    "gene_audit": ROOT
    / "03_plant_discovery/results/plant_confirmed_insectlike_sequence_GP_UP_expression_audit.tsv",
    "summary": ROOT
    / "04_cross_species_integration/results/plant_confirmed_insectlike_sequence_GP_UP_expression_summary.tsv",
    "parameters": ROOT / "00_project_control/09_confirmed_insectlike_sequence_UP_audit_parameters.tsv",
    "audit": ROOT / "08_logs/09_confirmed_insectlike_sequence_UP_audit.txt",
}

# These headers are the frozen primary normalized-count columns.
GP = ["LTS20_VA01_Sample_S7", "LTS20_VA01_Sample_S8", "LTS20_VA01_Sample_S9"]
UP = ["LTS20_VA01_Sample_S10", "LTS20_VA01_Sample_S11", "LTS20_VA01_Sample_S12"]
DISPLAY = {
    "LTS20_VA01_Sample_S7": "GP1",
    "LTS20_VA01_Sample_S8": "GP2",
    "LTS20_VA01_Sample_S9": "GP3",
    "LTS20_VA01_Sample_S10": "UP1",
    "LTS20_VA01_Sample_S11": "UP2",
    "LTS20_VA01_Sample_S12": "UP3",
}


class AuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise AuditError(f"Missing header: {path}")
        fields = [item.strip().lstrip("\ufeff") for item in reader.fieldnames]
        rows = []
        for row in reader:
            rows.append({k.strip().lstrip("\ufeff"): (v or "").strip() for k, v in row.items()})
    return fields, rows


def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def fmt(value: float) -> str:
    return format(value, ".10g")


def mean(values):
    return sum(values) / len(values)


def detected(values):
    return sum(value > 0 for value in values)


def detection_class(gp, up):
    gp_any, up_any = any(value > 0 for value in gp), any(value > 0 for value in up)
    if gp_any and not up_any:
        return "detected_in_GP_only"
    if gp_any and up_any:
        return "detected_in_both_GP_and_UP"
    if not gp_any and up_any:
        return "detected_in_UP_only"
    return "not_detected_in_normalized_count_table"


def main():
    for name, path in INPUTS.items():
        if not path.is_file():
            raise AuditError(f"Missing input {name}: {path}")
    input_hashes = {name: sha256(path) for name, path in INPUTS.items()}

    mask_fields, mask_rows = read_tsv(INPUTS["confirmed_mask"])
    required_mask = {
        "gene_id",
        "DESeq2_status",
        "host_reanalysis_mask_status",
    }
    if not required_mask.issubset(mask_fields):
        raise AuditError("Confirmed mask lacks required columns")
    confirmed = {
        row["gene_id"]: row
        for row in mask_rows
        if row["host_reanalysis_mask_status"]
        == "review_exclude_confirmed_multilevel_insect_like_sequence"
    }
    if len(confirmed) != 3948:
        raise AuditError(f"Expected 3,948 confirmed masked genes; found {len(confirmed)}")

    count_fields, count_rows = read_tsv(INPUTS["normalized_counts"])
    if "gene_id" not in count_fields:
        raise AuditError("Normalized-count table lacks gene_id")
    required_counts = set(GP + UP)
    if not required_counts.issubset(count_fields):
        raise AuditError("Expected frozen GP/UP normalized-count columns were not found")
    counts = {}
    all_totals = Counter()
    for row in count_rows:
        gene = row["gene_id"]
        if not gene or gene in counts:
            raise AuditError(f"Invalid/duplicate gene_id in normalized counts: {gene!r}")
        try:
            values = {sample: float(row[sample]) for sample in GP + UP}
        except ValueError as error:
            raise AuditError(f"Non-numeric normalized count for {gene}") from error
        if any(value < 0 or not math.isfinite(value) for value in values.values()):
            raise AuditError(f"Invalid normalized count for {gene}")
        counts[gene] = values
        all_totals.update(values)

    missing = sorted(set(confirmed) - set(counts))
    missing_excluded_low_count = [
        gene for gene in missing if confirmed[gene]["DESeq2_status"] == "excluded_low_count"
    ]
    unexpected_missing = sorted(set(missing) - set(missing_excluded_low_count))
    audit_rows = []
    masked_totals = Counter()
    class_counts = Counter()
    status_counts = Counter()
    for gene in sorted(confirmed):
        status = confirmed[gene]["DESeq2_status"]
        status_counts[status] += 1
        if gene not in counts:
            continue
        values = counts[gene]
        gp_values = [values[sample] for sample in GP]
        up_values = [values[sample] for sample in UP]
        gp_mean, up_mean = mean(gp_values), mean(up_values)
        ratio = math.log2((gp_mean + 0.5) / (up_mean + 0.5))
        cls = detection_class(gp_values, up_values)
        class_counts[cls] += 1
        masked_totals.update(values)
        audit_rows.append(
            {
                "gene_id": gene,
                "DESeq2_status": status,
                "GP1_normalized_count": fmt(gp_values[0]),
                "GP2_normalized_count": fmt(gp_values[1]),
                "GP3_normalized_count": fmt(gp_values[2]),
                "UP1_normalized_count": fmt(up_values[0]),
                "UP2_normalized_count": fmt(up_values[1]),
                "UP3_normalized_count": fmt(up_values[2]),
                "GP_mean_normalized_count": fmt(gp_mean),
                "UP_mean_normalized_count": fmt(up_mean),
                "descriptive_log2_GP_to_UP_mean_ratio_pseudocount_0.5": fmt(ratio),
                "GP_replicates_with_count_gt0": detected(gp_values),
                "UP_replicates_with_count_gt0": detected(up_values),
                "normalized_count_detection_class": cls,
            }
        )
    if unexpected_missing:
        raise AuditError(
            "Confirmed-mask genes unexpectedly missing from normalized-count table: "
            f"{len(unexpected_missing)}; " + ", ".join(unexpected_missing[:20])
        )

    audit_fields = [
        "gene_id",
        "DESeq2_status",
        "GP1_normalized_count",
        "GP2_normalized_count",
        "GP3_normalized_count",
        "UP1_normalized_count",
        "UP2_normalized_count",
        "UP3_normalized_count",
        "GP_mean_normalized_count",
        "UP_mean_normalized_count",
        "descriptive_log2_GP_to_UP_mean_ratio_pseudocount_0.5",
        "GP_replicates_with_count_gt0",
        "UP_replicates_with_count_gt0",
        "normalized_count_detection_class",
    ]
    write_tsv(OUTPUTS["gene_audit"], audit_fields, audit_rows)

    summary = []
    for status in sorted(status_counts):
        summary.append(
            {
                "section": "confirmed_mask_DESeq2_status",
                "metric": status,
                "count_or_value": status_counts[status],
                "detail": "Confirmed multilevel insect-like sequence genes",
            }
        )
    for cls in sorted(class_counts):
        summary.append(
            {
                "section": "confirmed_mask_normalized_count_detection",
                "metric": cls,
                "count_or_value": class_counts[cls],
                "detail": "Detection means normalized count >0 in one or more replicates",
            }
        )
    summary.append(
        {
            "section": "normalized_count_table_coverage",
            "metric": "confirmed_mask_genes_excluded_before_primary_normalized_count_table",
            "count_or_value": len(missing_excluded_low_count),
            "detail": (
                "These genes have frozen DESeq2 status excluded_low_count and are absent "
                "because the primary normalized-count table retains the post-filter gene set"
            ),
        }
    )
    for sample in GP + UP:
        fraction = 100 * masked_totals[sample] / all_totals[sample] if all_totals[sample] else float("nan")
        summary.append(
            {
                "section": "per_sample_normalized_count_burden",
                "metric": DISPLAY[sample],
                "count_or_value": fmt(masked_totals[sample]),
                "detail": (
                    f"Confirmed-mask normalized-count sum; {fmt(fraction)}% of the "
                    "total normalized-count sum in this sample"
                ),
            }
        )
    summary.append(
        {
            "section": "official_UP_up_status",
            "metric": "confirmed_mask_genes_classified_UP_up",
            "count_or_value": status_counts["UP_up"],
            "detail": "Existing frozen DESeq2 status; not a new test",
        }
    )
    write_tsv(OUTPUTS["summary"], ["section", "metric", "count_or_value", "detail"], summary)

    parameters = [
        {
            "parameter": "scope",
            "value": "Descriptive normalized-count audit only; no DESeq2 rerun or filtering",
        },
        {"parameter": "confirmed_mask_genes_expected", "value": "3948"},
        {"parameter": "normalized_count_detection_rule", "value": "count > 0"},
        {
            "parameter": "mean_ratio_note",
            "value": "Pseudocount 0.5 is used only for a finite descriptive log2 mean ratio",
        },
        {"parameter": "GP_samples", "value": "GP1=S7, GP2=S8, GP3=S9"},
        {"parameter": "UP_samples", "value": "UP1=S10, UP2=S11, UP3=S12"},
    ]
    parameters.extend(
        {"parameter": "SHA256_" + name, "value": value}
        for name, value in sorted(input_hashes.items())
    )
    write_tsv(OUTPUTS["parameters"], ["parameter", "value"], parameters)

    for name, path in INPUTS.items():
        if sha256(path) != input_hashes[name]:
            raise AuditError(f"Input changed during audit: {name}")

    report = [
        "===== CONFIRMED INSECT-LIKE SEQUENCE GP/UP EXPRESSION AUDIT =====",
        f"Discovery root: {ROOT}",
        "",
        "Scope:",
        "  Descriptive audit of frozen DESeq2-normalized counts; no new statistical test or DEG call.",
        "",
        f"Confirmed multilevel insect-like sequence genes: {len(confirmed)}",
        f"Confirmed genes found in normalized-count table: {len(audit_rows)}",
        (
            "Confirmed genes absent because of the frozen low-count filter: "
            f"{len(missing_excluded_low_count)}"
        ),
        "",
        "Confirmed genes by existing frozen DESeq2 status:",
    ]
    for status in sorted(status_counts):
        report.append(f"  {status}: {status_counts[status]}")
    report.extend(["", "Normalized-count detection across conditions:"])
    for cls in sorted(class_counts):
        report.append(f"  {cls}: {class_counts[cls]}")
    report.extend(["", "Confirmed-mask normalized-count burden by sample:"])
    for sample in GP + UP:
        fraction = 100 * masked_totals[sample] / all_totals[sample] if all_totals[sample] else float("nan")
        report.append(
            f"  {DISPLAY[sample]}: {fmt(masked_totals[sample])} "
            f"({fmt(fraction)}% of total normalized-count sum)"
        )
    report.extend(
        [
            "",
            "Interpretive boundary:",
            "  Detection in an UP library, if present, is not a taxonomic assignment or proof of contamination.",
            "  The audit asks whether the confirmed sequence component is GP-concentrated; it does not identify its physical RNA source.",
            "  Existing DESeq2 results and the confirmed mask were not changed.",
            "",
            "Input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "",
        ]
    )
    for path in OUTPUTS.values():
        report.append(f"Output: {path}")
    report.append("CONFIRMED INSECT-LIKE SEQUENCE GP/UP EXPRESSION AUDIT: PASS")
    OUTPUTS["audit"].parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUTS["audit"].with_name(OUTPUTS["audit"].name + ".tmp")
    temporary.write_text("\n".join(report) + "\n", encoding="utf-8")
    os.replace(temporary, OUTPUTS["audit"])
    print("\n".join(report))


if __name__ == "__main__":
    try:
        main()
    except AuditError as error:
        print("===== CONFIRMED INSECT-LIKE SEQUENCE GP/UP EXPRESSION AUDIT =====", file=sys.stderr)
        print(f"ERROR: {error}", file=sys.stderr)
        print("CONFIRMED INSECT-LIKE SEQUENCE GP/UP EXPRESSION AUDIT: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
