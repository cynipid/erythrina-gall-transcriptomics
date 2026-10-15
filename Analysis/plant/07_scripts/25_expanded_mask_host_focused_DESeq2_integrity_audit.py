#!/usr/bin/env python3
"""Independent integrity audit of the expanded-mask host-focused DESeq2 refit."""

from __future__ import annotations

import csv
import hashlib
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path


HOME = Path.home()
ROOT = HOME / "projects/qe_cross_species_discovery"
PLANT = HOME / "projects/qe_plant_denovo"
BRANCH = ROOT / "05_host_focused_DESeq2_sensitivity"
RESULTS = BRANCH / "results"
LOGS = BRANCH / "logs"

INPUTS = {
    "expanded_mask": ROOT
    / "04_cross_species_integration/results/plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv",
    "current_host_results": RESULTS
    / "plant_gene_DESeq2_GP_vs_UP_host_focused_all.tsv",
    "expanded_results": RESULTS
    / "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_all.tsv",
    "expanded_normalized_counts": RESULTS
    / "plant_gene_normalized_counts_host_focused_expanded_mask.tsv",
    "filtering_summary": RESULTS
    / "plant_host_focused_expanded_mask_DESeq2_filtering_summary.tsv",
    "mask_exclusion_audit": RESULTS
    / "plant_expanded_nucleotide_mask_exclusion_audit.tsv",
    "current_comparison": RESULTS
    / "plant_host_focused_expanded_mask_vs_current_DESeq2_comparison.tsv",
    "r_input_identity": LOGS
    / "24_host_focused_expanded_mask_DESeq2_input_identity.tsv",
    "r_audit": LOGS / "24_host_focused_expanded_mask_DESeq2_audit.txt",
}
R_ANALYSIS_INPUTS = {
    "plant_gene_txi.rds": PLANT / "07_deseq2/plant_gene_txi.rds",
    "plant_deseq2_samples.tsv": PLANT
    / "01_input/expression/plant_deseq2_samples.tsv",
    "plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv": INPUTS[
        "expanded_mask"
    ],
    "plant_gene_DESeq2_GP_vs_UP_host_focused_all.tsv": INPUTS[
        "current_host_results"
    ],
}
OUTPUTS = {
    "status_contingency": RESULTS
    / "plant_host_focused_current_vs_expanded_mask_status_contingency.tsv",
    "integrity_summary": ROOT
    / "00_project_control/25_expanded_mask_host_focused_DESeq2_integrity_summary.tsv",
    "audit": LOGS / "25_expanded_mask_host_focused_DESeq2_integrity_audit.txt",
}

EXPECTED_CURRENT = {
    "GP_up": 3179,
    "UP_up": 1396,
    "non_DE": 21008,
    "NA_padj_not_testable": 196,
}
EXPECTED_EXPANDED = {
    "GP_up": 2349,
    "UP_up": 1386,
    "non_DE": 21022,
    "NA_padj_not_testable": 174,
}
EXPECTED_REMOVED = {
    "GP_up": 813,
    "UP_up": 0,
    "non_DE": 14,
    "NA_padj_not_testable": 21,
}
EXPECTED_SAMPLES = [
    "LTS20_VA01_Sample_S7",
    "LTS20_VA01_Sample_S8",
    "LTS20_VA01_Sample_S9",
    "LTS20_VA01_Sample_S10",
    "LTS20_VA01_Sample_S11",
    "LTS20_VA01_Sample_S12",
]


class AuditError(RuntimeError):
    pass


def digest(path: Path, algorithm: str = "sha256") -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditError(f"Missing TSV header: {path}")
        fields = [field.strip().lstrip("\ufeff") for field in reader.fieldnames]
        rows = [
            {
                key.strip().lstrip("\ufeff"): (value or "").strip()
                for key, value in row.items()
            }
            for row in reader
        ]
    return fields, rows


def require_columns(path: Path, fields: list[str], required: set[str]) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise AuditError(f"Missing columns in {path}: {', '.join(missing)}")


def index_unique(
    rows: list[dict[str, str]], key: str, label: str
) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        value = row.get(key, "")
        if not value:
            raise AuditError(f"Blank {key} in {label}")
        if value in indexed:
            raise AuditError(f"Duplicate {key} in {label}: {value}")
        indexed[value] = row
    return indexed


def missing_number(value: str) -> bool:
    return value.strip().upper() in {"", "NA", "NAN"}


def classify(padj: str, log2fc: str) -> str:
    if missing_number(padj):
        return "NA_padj_not_testable"
    if missing_number(log2fc):
        raise AuditError("A tested result has missing log2FoldChange")
    p_value = float(padj)
    effect = float(log2fc)
    if not math.isfinite(p_value) or not math.isfinite(effect):
        raise AuditError("A tested result contains a non-finite padj/log2FoldChange")
    if p_value < 0.01 and effect >= 1.0:
        return "GP_up"
    if p_value < 0.01 and effect <= -1.0:
        return "UP_up"
    return "non_DE"


def parse_bool(value: str, label: str) -> bool:
    normalized = value.strip().upper()
    if normalized == "TRUE":
        return True
    if normalized == "FALSE":
        return False
    raise AuditError(f"Invalid Boolean value in {label}: {value!r}")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
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


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def main() -> None:
    for label, path in {**INPUTS, **R_ANALYSIS_INPUTS}.items():
        if not path.is_file():
            raise AuditError(f"Missing required input {label}: {path}")
    hashes_before = {label: digest(path) for label, path in INPUTS.items()}

    mask_fields, mask_rows = read_tsv(INPUTS["expanded_mask"])
    current_fields, current_rows = read_tsv(INPUTS["current_host_results"])
    expanded_fields, expanded_rows = read_tsv(INPUTS["expanded_results"])
    norm_fields, norm_rows = read_tsv(INPUTS["expanded_normalized_counts"])
    summary_fields, summary_rows = read_tsv(INPUTS["filtering_summary"])
    exclusion_fields, exclusion_rows = read_tsv(INPUTS["mask_exclusion_audit"])
    comparison_fields, comparison_rows = read_tsv(INPUTS["current_comparison"])
    identity_fields, identity_rows = read_tsv(INPUTS["r_input_identity"])

    require_columns(
        INPUTS["expanded_mask"],
        mask_fields,
        {"gene_id", "proposed_expanded_mask_status"},
    )
    result_required = {"gene_id", "baseMean", "log2FoldChange", "padj"}
    require_columns(
        INPUTS["current_host_results"],
        current_fields,
        result_required | {"DESeq2_status_host_focused"},
    )
    require_columns(
        INPUTS["expanded_results"],
        expanded_fields,
        result_required | {"DESeq2_status_host_focused_expanded_mask"},
    )
    require_columns(
        INPUTS["expanded_normalized_counts"],
        norm_fields,
        {"gene_id", *EXPECTED_SAMPLES},
    )
    require_columns(
        INPUTS["filtering_summary"], summary_fields, {"metric", "value"}
    )
    require_columns(
        INPUTS["mask_exclusion_audit"],
        exclusion_fields,
        {
            "gene_id",
            "mask_source",
            "present_in_original_tximport",
            "removed_before_low_count_filter",
            "retained_in_expanded_mask_results",
        },
    )
    require_columns(
        INPUTS["current_comparison"],
        comparison_fields,
        {
            "gene_id",
            "DESeq2_status_host_focused",
            "DESeq2_status_host_focused_expanded_mask",
            "additional_strong_nucleotide_masked",
            "comparison_status",
        },
    )
    require_columns(
        INPUTS["r_input_identity"],
        identity_fields,
        {"input", "path", "MD5_before_analysis"},
    )

    mask = index_unique(mask_rows, "gene_id", "expanded mask")
    current = index_unique(current_rows, "gene_id", "current host results")
    expanded = index_unique(expanded_rows, "gene_id", "expanded results")
    normalized = index_unique(norm_rows, "gene_id", "expanded normalized counts")
    exclusion = index_unique(exclusion_rows, "gene_id", "mask exclusion audit")
    comparison = index_unique(comparison_rows, "gene_id", "comparison table")
    identity = index_unique(identity_rows, "input", "R input identity")

    errors: list[str] = []
    if len(mask) != 118_082:
        errors.append(f"Expanded mask universe has {len(mask)} rows; expected 118082")
    existing_mask = {
        gene
        for gene, row in mask.items()
        if row["proposed_expanded_mask_status"]
        == "review_exclude_existing_confirmed_multilevel_sequence"
    }
    additional_mask = {
        gene
        for gene, row in mask.items()
        if row["proposed_expanded_mask_status"]
        == "review_exclude_additional_strong_nucleotide_similarity"
    }
    expanded_mask = existing_mask | additional_mask
    if len(existing_mask) != 3_948:
        errors.append(f"Existing mask count is {len(existing_mask)}; expected 3948")
    if len(additional_mask) != 12_164:
        errors.append(f"Additional mask count is {len(additional_mask)}; expected 12164")
    if len(expanded_mask) != 16_112:
        errors.append(f"Expanded mask count is {len(expanded_mask)}; expected 16112")

    if len(current) != 25_779:
        errors.append(f"Current host result rows: {len(current)}; expected 25779")
    if len(expanded) != 24_931:
        errors.append(f"Expanded result rows: {len(expanded)}; expected 24931")
    if set(expanded) != set(normalized):
        errors.append("Expanded result and normalized-count gene sets differ")
    if set(expanded) & expanded_mask:
        errors.append("Expanded-mask genes remain in expanded DESeq2 results")
    if set(normalized) & expanded_mask:
        errors.append("Expanded-mask genes remain in normalized counts")
    if set(exclusion) != expanded_mask:
        errors.append("Exclusion-audit membership differs from expanded mask")
    if set(comparison) != set(current):
        errors.append("Comparison table membership differs from current host results")

    current_calls: dict[str, str] = {}
    for gene, row in current.items():
        calculated = classify(row["padj"], row["log2FoldChange"])
        current_calls[gene] = calculated
        if calculated != row["DESeq2_status_host_focused"]:
            errors.append(f"Current stored DEG status mismatch: {gene}")
            if len(errors) >= 50:
                break
    expanded_calls: dict[str, str] = {}
    for gene, row in expanded.items():
        calculated = classify(row["padj"], row["log2FoldChange"])
        expanded_calls[gene] = calculated
        if calculated != row["DESeq2_status_host_focused_expanded_mask"]:
            errors.append(f"Expanded stored DEG status mismatch: {gene}")
            if len(errors) >= 50:
                break

    current_counts = Counter(current_calls.values())
    expanded_counts = Counter(expanded_calls.values())
    if dict(current_counts) != EXPECTED_CURRENT:
        errors.append(
            f"Current status counts differ: {dict(current_counts)}; "
            f"expected {EXPECTED_CURRENT}"
        )
    if dict(expanded_counts) != EXPECTED_EXPANDED:
        errors.append(
            f"Expanded status counts differ: {dict(expanded_counts)}; "
            f"expected {EXPECTED_EXPANDED}"
        )

    directly_removed = set(current).difference(expanded)
    expected_directly_removed = additional_mask.intersection(current)
    if directly_removed != expected_directly_removed:
        errors.append(
            "Current-to-expanded gene loss does not exactly equal the additional "
            "mask among tested genes"
        )
    removed_counts = Counter(current_calls[gene] for gene in directly_removed)
    if dict(removed_counts) != {key: value for key, value in EXPECTED_REMOVED.items() if value}:
        errors.append(
            f"Direct-removal status counts differ: {dict(removed_counts)}; "
            f"expected {EXPECTED_REMOVED}"
        )

    for row in exclusion_rows:
        try:
            if not parse_bool(row["present_in_original_tximport"], "exclusion audit"):
                errors.append(f"Masked gene absent from original tximport: {row['gene_id']}")
            if not parse_bool(row["removed_before_low_count_filter"], "exclusion audit"):
                errors.append(f"Masked gene not marked removed: {row['gene_id']}")
            if parse_bool(row["retained_in_expanded_mask_results"], "exclusion audit"):
                errors.append(f"Masked gene marked retained: {row['gene_id']}")
        except AuditError as error:
            errors.append(str(error))
        if len(errors) >= 50:
            break

    # Validate normalized-count columns and numerical values.
    if norm_fields != ["gene_id", *EXPECTED_SAMPLES]:
        errors.append(
            "Normalized-count columns/order differ from the six frozen sample IDs"
        )
    invalid_count_values = 0
    for row in norm_rows:
        for sample in EXPECTED_SAMPLES:
            try:
                value = float(row[sample])
            except ValueError:
                invalid_count_values += 1
                continue
            if not math.isfinite(value) or value < 0:
                invalid_count_values += 1
    if invalid_count_values:
        errors.append(
            f"Normalized-count table contains {invalid_count_values} invalid values"
        )

    summary = index_unique(summary_rows, "metric", "filtering summary")
    expected_summary = {
        "tximport_genes_before_mask": 118_082,
        "existing_confirmed_multilevel_mask_genes": 3_948,
        "additional_strong_nucleotide_mask_genes": 12_164,
        "expanded_mask_genes_total": 16_112,
        "expanded_mask_genes_present_in_tximport": 16_112,
        "genes_after_expanded_mask_before_low_count_filter": 101_970,
        "genes_after_frozen_low_count_filter": 24_931,
        "expanded_mask_host_focused_total_DEGs": 3_735,
        "expanded_mask_host_focused_GP_up_DEGs": 2_349,
        "expanded_mask_host_focused_UP_up_DEGs": 1_386,
        "expanded_mask_host_focused_non_DE": 21_022,
        "expanded_mask_host_focused_NA_padj_not_testable": 174,
    }
    for metric, expected in expected_summary.items():
        if metric not in summary:
            errors.append(f"Filtering summary lacks metric: {metric}")
            continue
        try:
            observed = int(float(summary[metric]["value"]))
        except ValueError:
            errors.append(f"Filtering summary has nonnumeric value for {metric}")
            continue
        if observed != expected:
            errors.append(
                f"Filtering summary mismatch for {metric}: {observed} != {expected}"
            )

    # Confirm that the R-side input identity records still match the source files.
    r_md5_ok = True
    if set(identity) != set(R_ANALYSIS_INPUTS):
        errors.append("R input-identity membership does not match expected inputs")
        r_md5_ok = False
    else:
        for name, path in R_ANALYSIS_INPUTS.items():
            observed_md5 = identity[name]["MD5_before_analysis"].lower()
            current_md5 = digest(path, "md5")
            if observed_md5 != current_md5:
                errors.append(f"R input MD5 mismatch after analysis: {name}")
                r_md5_ok = False

    r_audit_lines = INPUTS["r_audit"].read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()
    if not r_audit_lines or r_audit_lines[-1] != (
        "EXPANDED-MASK HOST-FOCUSED PLANT DESEQ2 SENSITIVITY ANALYSIS: PASS"
    ):
        errors.append("R analysis audit does not end in PASS")

    transition_counts: dict[str, Counter[str]] = defaultdict(Counter)
    changed_retained = 0
    for gene, old_status in current_calls.items():
        if gene in directly_removed:
            new_status = "removed_by_expanded_mask"
        else:
            new_status = expanded_calls[gene]
            if new_status != old_status:
                changed_retained += 1
        transition_counts[old_status][new_status] += 1

        comparison_row = comparison[gene]
        try:
            comparison_mask_flag = parse_bool(
                comparison_row["additional_strong_nucleotide_masked"],
                "comparison table",
            )
        except AuditError as error:
            errors.append(str(error))
            comparison_mask_flag = False
        if comparison_mask_flag != (gene in additional_mask):
            errors.append(f"Comparison mask-membership mismatch: {gene}")
        expected_comparison = (
            "removed_by_expanded_mask"
            if gene in directly_removed
            else (
                "status_unchanged"
                if old_status == new_status
                else "status_changed_after_refit"
            )
        )
        if comparison_row["comparison_status"] != expected_comparison:
            errors.append(f"Comparison-status mismatch: {gene}")
            if len(errors) >= 50:
                break

    contingency_rows: list[dict[str, object]] = []
    status_order = [
        "GP_up",
        "UP_up",
        "non_DE",
        "NA_padj_not_testable",
    ]
    new_order = [*status_order, "removed_by_expanded_mask"]
    for old_status in status_order:
        for new_status in new_order:
            count = transition_counts[old_status][new_status]
            if count:
                contingency_rows.append(
                    {
                        "current_3948_mask_status": old_status,
                        "expanded_16112_mask_status": new_status,
                        "genes": count,
                    }
                )
    write_tsv(
        OUTPUTS["status_contingency"],
        ["current_3948_mask_status", "expanded_16112_mask_status", "genes"],
        contingency_rows,
    )

    stable = all(
        digest(path) == hashes_before[label] for label, path in INPUTS.items()
    )
    if not stable:
        errors.append("One or more audited input files changed during the audit")

    integrity_rows = [
        {"metric": "current_tested_genes", "value": len(current)},
        {"metric": "expanded_tested_genes", "value": len(expanded)},
        {"metric": "directly_removed_tested_genes", "value": len(directly_removed)},
        {"metric": "retained_genes_with_changed_status", "value": changed_retained},
        {"metric": "current_total_DEGs", "value": current_counts["GP_up"] + current_counts["UP_up"]},
        {"metric": "expanded_total_DEGs", "value": expanded_counts["GP_up"] + expanded_counts["UP_up"]},
        {"metric": "current_GP_up", "value": current_counts["GP_up"]},
        {"metric": "expanded_GP_up", "value": expanded_counts["GP_up"]},
        {"metric": "current_UP_up", "value": current_counts["UP_up"]},
        {"metric": "expanded_UP_up", "value": expanded_counts["UP_up"]},
        {"metric": "errors_detected", "value": len(errors)},
    ]
    write_tsv(OUTPUTS["integrity_summary"], ["metric", "value"], integrity_rows)

    report = [
        "===== EXPANDED-MASK HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT =====",
        f"Discovery root: {ROOT}",
        "",
        f"Expanded mask genes: {len(expanded_mask)}",
        f"Current tested genes: {len(current)}",
        f"Expanded-mask tested genes: {len(expanded)}",
        f"Tested genes removed directly by expanded mask: {len(directly_removed)}",
        "",
        "Current 3,948-mask statuses:",
        f"  GP_up: {current_counts['GP_up']}",
        f"  UP_up: {current_counts['UP_up']}",
        f"  non_DE: {current_counts['non_DE']}",
        f"  NA_padj_not_testable: {current_counts['NA_padj_not_testable']}",
        "",
        "Expanded 16,112-mask statuses:",
        f"  GP_up: {expanded_counts['GP_up']}",
        f"  UP_up: {expanded_counts['UP_up']}",
        f"  non_DE: {expanded_counts['non_DE']}",
        f"  NA_padj_not_testable: {expanded_counts['NA_padj_not_testable']}",
        "",
        "Direct removals by previous status:",
        f"  GP_up: {removed_counts['GP_up']}",
        f"  UP_up: {removed_counts['UP_up']}",
        f"  non_DE: {removed_counts['non_DE']}",
        f"  NA_padj_not_testable: {removed_counts['NA_padj_not_testable']}",
        f"Retained genes changing status after refit: {changed_retained}",
        "",
        "Nonzero current-to-expanded transitions:",
    ]
    for row in contingency_rows:
        report.append(
            f"  {row['current_3948_mask_status']} -> "
            f"{row['expanded_16112_mask_status']}: {row['genes']}"
        )
    report.extend(
        [
            "",
            f"Masked genes retained in expanded results: {len(set(expanded) & expanded_mask)}",
            f"Result/count identifier agreement: {'PASS' if set(expanded) == set(normalized) else 'FAIL'}",
            f"R input MD5 verification: {'PASS' if r_md5_ok else 'FAIL'}",
            f"Input SHA-256 verification after output generation: {'PASS' if stable else 'FAIL'}",
            f"Errors detected: {len(errors)}",
        ]
    )
    report.extend(f"ERROR: {error}" for error in errors)
    report.extend(
        [
            "",
            f"Output: {OUTPUTS['status_contingency']}",
            f"Output: {OUTPUTS['integrity_summary']}",
            f"Output: {OUTPUTS['audit']}",
            "EXPANDED-MASK HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT: "
            + ("PASS" if not errors and stable else "CHECK REQUIRED"),
        ]
    )
    write_text(OUTPUTS["audit"], "\n".join(report) + "\n")
    print("\n".join(report))
    if errors or not stable:
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        main()
    except (AuditError, OSError, ValueError, KeyError) as error:
        print(
            "===== EXPANDED-MASK HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT =====\n"
            f"ERROR: {error}\n"
            "EXPANDED-MASK HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT: CHECK REQUIRED",
            file=sys.stderr,
        )
        raise SystemExit(1)
