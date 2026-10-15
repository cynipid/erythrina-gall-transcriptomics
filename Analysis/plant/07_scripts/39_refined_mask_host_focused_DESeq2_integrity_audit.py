#!/usr/bin/env python3
"""Independent read-only audit of the final refined-mask plant DESeq2 refit."""

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
INTEGRATION = ROOT / "04_cross_species_integration/results"
BRANCH = ROOT / "05_host_focused_DESeq2_sensitivity"
RESULTS = BRANCH / "results"
LOGS = BRANCH / "logs"
CONTROL = ROOT / "00_project_control"

INPUTS = {
    "coding_candidates": INTEGRATION
    / "plant_coding_denominator_additional_sequence_review_candidates.tsv",
    "refined_mask": INTEGRATION
    / "plant_refined_coding_denominator_sequence_review_mask.tsv",
    "expanded_results": RESULTS
    / "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_all.tsv",
    "refined_results": RESULTS
    / "plant_gene_DESeq2_GP_vs_UP_host_focused_refined_mask_all.tsv",
    "refined_normalized_counts": RESULTS
    / "plant_gene_normalized_counts_host_focused_refined_mask.tsv",
    "filtering_summary": RESULTS
    / "plant_host_focused_refined_mask_DESeq2_filtering_summary.tsv",
    "mask_exclusion_audit": RESULTS
    / "plant_refined_coding_denominator_mask_exclusion_audit.tsv",
    "expanded_comparison": RESULTS
    / "plant_host_focused_expanded_vs_refined_mask_DESeq2_comparison.tsv",
    "r_input_identity": LOGS
    / "38_host_focused_refined_mask_DESeq2_input_identity.tsv",
    "r_audit": LOGS / "38_host_focused_refined_mask_DESeq2_audit.txt",
}

R_ANALYSIS_INPUTS = {
    "plant_gene_txi.rds": PLANT / "07_deseq2/plant_gene_txi.rds",
    "plant_deseq2_samples.tsv": PLANT
    / "01_input/expression/plant_deseq2_samples.tsv",
    "plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv": INTEGRATION
    / "plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv",
    "plant_coding_denominator_additional_sequence_review_candidates.tsv": INPUTS[
        "coding_candidates"
    ],
    "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_all.tsv": INPUTS[
        "expanded_results"
    ],
}

OUTPUTS = {
    "status_contingency": RESULTS
    / "plant_host_focused_expanded_vs_refined_mask_status_contingency.tsv",
    "integrity_summary": CONTROL
    / "39_refined_mask_host_focused_DESeq2_integrity_summary.tsv",
    "audit": LOGS / "39_refined_mask_host_focused_DESeq2_integrity_audit.txt",
}

EXPECTED_EXPANDED = {
    "GP_up": 2349,
    "UP_up": 1386,
    "non_DE": 21022,
    "NA_padj_not_testable": 174,
}
EXPECTED_REFINED = {
    "GP_up": 2147,
    "UP_up": 1383,
    "non_DE": 21032,
    "NA_padj_not_testable": 172,
}
EXPECTED_DIRECT_REMOVALS = {
    "GP_up": 192,
    "UP_up": 0,
    "non_DE": 3,
    "NA_padj_not_testable": 2,
}
EXPECTED_SAMPLES = [
    "LTS20_VA01_Sample_S7",
    "LTS20_VA01_Sample_S8",
    "LTS20_VA01_Sample_S9",
    "LTS20_VA01_Sample_S10",
    "LTS20_VA01_Sample_S11",
    "LTS20_VA01_Sample_S12",
]
STATUS_ORDER = ["GP_up", "UP_up", "non_DE", "NA_padj_not_testable"]


class AuditError(RuntimeError):
    pass


def digest(path: Path, algorithm: str = "sha256") -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise AuditError(f"Missing required input: {path}")
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


def unique_index(
    rows: list[dict[str, str]], key: str, label: str
) -> dict[str, dict[str, str]]:
    output: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = row.get(key, "")
        if not identifier:
            raise AuditError(f"Blank {key} in {label}")
        if identifier in output:
            raise AuditError(f"Duplicate {key} in {label}: {identifier}")
        output[identifier] = row
    return output


def is_missing(value: str) -> bool:
    return value.strip().upper() in {"", "NA", "NAN", "NONE", "NULL"}


def classify(padj: str, log2fc: str) -> str:
    if is_missing(padj):
        return "NA_padj_not_testable"
    if is_missing(log2fc):
        raise AuditError("A result with valid padj has missing log2FoldChange")
    try:
        adjusted = float(padj)
        effect = float(log2fc)
    except ValueError as error:
        raise AuditError(f"Invalid DESeq2 numeric value: {error}") from error
    if not math.isfinite(adjusted) or not math.isfinite(effect):
        raise AuditError("A result contains non-finite padj/log2FoldChange")
    if adjusted < 0.01 and effect >= 1.0:
        return "GP_up"
    if adjusted < 0.01 and effect <= -1.0:
        return "UP_up"
    return "non_DE"


def parse_bool(value: str, label: str) -> bool:
    normalized = value.strip().upper()
    if normalized == "TRUE":
        return True
    if normalized == "FALSE":
        return False
    raise AuditError(f"Invalid Boolean in {label}: {value!r}")


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
    for path in [*INPUTS.values(), *R_ANALYSIS_INPUTS.values()]:
        if not path.is_file():
            raise AuditError(f"Missing required input: {path}")
    hashes_before = {label: digest(path) for label, path in INPUTS.items()}

    candidate_fields, candidate_rows = read_tsv(INPUTS["coding_candidates"])
    mask_fields, mask_rows = read_tsv(INPUTS["refined_mask"])
    expanded_fields, expanded_rows = read_tsv(INPUTS["expanded_results"])
    refined_fields, refined_rows = read_tsv(INPUTS["refined_results"])
    norm_fields, norm_rows = read_tsv(INPUTS["refined_normalized_counts"])
    summary_fields, summary_rows = read_tsv(INPUTS["filtering_summary"])
    exclusion_fields, exclusion_rows = read_tsv(INPUTS["mask_exclusion_audit"])
    comparison_fields, comparison_rows = read_tsv(INPUTS["expanded_comparison"])
    identity_fields, identity_rows = read_tsv(INPUTS["r_input_identity"])

    require_columns(
        INPUTS["coding_candidates"],
        candidate_fields,
        {"gene_id", "current_expanded_mask_DESeq2_status", "proposed_review_status"},
    )
    require_columns(
        INPUTS["refined_mask"],
        mask_fields,
        {
            "gene_id",
            "proposed_expanded_mask_status",
            "additional_coding_denominator_mask_overlap",
            "refined_mask_status",
        },
    )
    result_fields = {"gene_id", "baseMean", "log2FoldChange", "padj"}
    require_columns(
        INPUTS["expanded_results"],
        expanded_fields,
        result_fields | {"DESeq2_status_host_focused_expanded_mask"},
    )
    require_columns(
        INPUTS["refined_results"],
        refined_fields,
        result_fields | {"DESeq2_status_host_focused_refined_mask"},
    )
    require_columns(
        INPUTS["refined_normalized_counts"],
        norm_fields,
        {"gene_id", *EXPECTED_SAMPLES},
    )
    require_columns(INPUTS["filtering_summary"], summary_fields, {"metric", "value"})
    require_columns(
        INPUTS["mask_exclusion_audit"],
        exclusion_fields,
        {
            "gene_id",
            "mask_source",
            "present_in_original_tximport",
            "removed_before_low_count_filter",
            "retained_in_refined_mask_results",
        },
    )
    require_columns(
        INPUTS["expanded_comparison"],
        comparison_fields,
        {
            "gene_id",
            "DESeq2_status_host_focused_expanded_mask",
            "DESeq2_status_host_focused_refined_mask",
            "additional_coding_denominator_masked",
            "comparison_status",
        },
    )
    require_columns(
        INPUTS["r_input_identity"],
        identity_fields,
        {"input", "path", "MD5_before_analysis"},
    )

    candidates = unique_index(candidate_rows, "gene_id", "coding candidates")
    mask = unique_index(mask_rows, "gene_id", "refined mask")
    expanded = unique_index(expanded_rows, "gene_id", "expanded results")
    refined = unique_index(refined_rows, "gene_id", "refined results")
    normalized = unique_index(norm_rows, "gene_id", "normalized counts")
    exclusion = unique_index(exclusion_rows, "gene_id", "exclusion audit")
    comparison = unique_index(comparison_rows, "gene_id", "comparison table")
    identity = unique_index(identity_rows, "input", "R input identity")

    errors: list[str] = []

    expanded_mask = {
        gene
        for gene, row in mask.items()
        if row["proposed_expanded_mask_status"]
        in {
            "review_exclude_existing_confirmed_multilevel_sequence",
            "review_exclude_additional_strong_nucleotide_similarity",
        }
    }
    coding_mask = {
        gene
        for gene, row in mask.items()
        if row["additional_coding_denominator_mask_overlap"].upper() == "TRUE"
    }
    refined_mask = {
        gene
        for gene, row in mask.items()
        if row["refined_mask_status"]
        in {
            "review_exclude_frozen_expanded_nucleotide_mask",
            "review_exclude_additional_coding_denominator_supported_sequence",
        }
    }

    if len(mask) != 118_082:
        errors.append(f"Refined-mask universe has {len(mask)} genes; expected 118082")
    if len(expanded_mask) != 16_112:
        errors.append(f"Frozen expanded mask has {len(expanded_mask)} genes; expected 16112")
    if len(candidates) != 292 or len(coding_mask) != 292:
        errors.append(
            f"Coding additions differ: candidates={len(candidates)}, mask={len(coding_mask)}; expected 292"
        )
    if set(candidates) != coding_mask:
        errors.append("Coding-candidate membership differs from refined-mask additions")
    if expanded_mask & coding_mask:
        errors.append("Expanded-mask and coding-denominator additions overlap")
    if len(refined_mask) != 16_404 or refined_mask != expanded_mask | coding_mask:
        errors.append("Refined-mask membership is not the exact 16,404-gene union")

    if len(expanded) != 24_931:
        errors.append(f"Expanded result rows: {len(expanded)}; expected 24931")
    if len(refined) != 24_734:
        errors.append(f"Refined result rows: {len(refined)}; expected 24734")
    if set(refined) != set(normalized):
        errors.append("Refined result and normalized-count identifiers differ")
    if set(refined) & refined_mask:
        errors.append("One or more refined-mask genes remain in refined results")
    if set(normalized) & refined_mask:
        errors.append("One or more refined-mask genes remain in normalized counts")
    if set(exclusion) != refined_mask:
        errors.append("Exclusion-audit membership differs from the refined mask")
    if set(comparison) != set(expanded):
        errors.append("Comparison-table membership differs from expanded results")

    expanded_calls: dict[str, str] = {}
    for gene, row in expanded.items():
        calculated = classify(row["padj"], row["log2FoldChange"])
        expanded_calls[gene] = calculated
        if calculated != row["DESeq2_status_host_focused_expanded_mask"]:
            errors.append(f"Expanded stored status mismatch: {gene}")
            if len(errors) >= 50:
                break

    refined_calls: dict[str, str] = {}
    for gene, row in refined.items():
        calculated = classify(row["padj"], row["log2FoldChange"])
        refined_calls[gene] = calculated
        if calculated != row["DESeq2_status_host_focused_refined_mask"]:
            errors.append(f"Refined stored status mismatch: {gene}")
            if len(errors) >= 50:
                break

    expanded_counts = Counter(expanded_calls.values())
    refined_counts = Counter(refined_calls.values())
    if dict(expanded_counts) != EXPECTED_EXPANDED:
        errors.append(
            f"Expanded status counts differ: {dict(expanded_counts)} != {EXPECTED_EXPANDED}"
        )
    if dict(refined_counts) != EXPECTED_REFINED:
        errors.append(
            f"Refined status counts differ: {dict(refined_counts)} != {EXPECTED_REFINED}"
        )

    directly_removed = set(expanded) - set(refined)
    expected_removed = coding_mask & set(expanded)
    if directly_removed != expected_removed:
        errors.append("Expanded-to-refined gene loss is not exactly the 197 tested coding additions")
    removed_counts = Counter(expanded_calls[gene] for gene in directly_removed)
    expected_nonzero = {key: value for key, value in EXPECTED_DIRECT_REMOVALS.items() if value}
    if dict(removed_counts) != expected_nonzero:
        errors.append(
            f"Direct-removal statuses differ: {dict(removed_counts)} != {expected_nonzero}"
        )

    for gene, row in candidates.items():
        expected_status = (
            expanded_calls[gene]
            if gene in expanded_calls
            else "not_in_expanded_mask_host_focused_tested_universe"
        )
        if row["current_expanded_mask_DESeq2_status"] != expected_status:
            errors.append(f"Candidate status disagrees with expanded results: {gene}")
        if row["proposed_review_status"] != "review_add_coding_denominator_supported_sequence":
            errors.append(f"Unexpected candidate review status: {gene}")
        if len(errors) >= 50:
            break

    for gene, row in exclusion.items():
        try:
            if not parse_bool(row["present_in_original_tximport"], "exclusion audit"):
                errors.append(f"Excluded gene absent from original tximport: {gene}")
            if not parse_bool(row["removed_before_low_count_filter"], "exclusion audit"):
                errors.append(f"Excluded gene not marked removed: {gene}")
            if parse_bool(row["retained_in_refined_mask_results"], "exclusion audit"):
                errors.append(f"Excluded gene marked retained: {gene}")
        except AuditError as error:
            errors.append(str(error))
        if len(errors) >= 50:
            break

    if norm_fields != ["gene_id", *EXPECTED_SAMPLES]:
        errors.append("Normalized-count columns/order differ from the six frozen samples")
    invalid_counts = 0
    for row in norm_rows:
        for sample in EXPECTED_SAMPLES:
            try:
                value = float(row[sample])
            except ValueError:
                invalid_counts += 1
                continue
            if not math.isfinite(value) or value < 0:
                invalid_counts += 1
    if invalid_counts:
        errors.append(f"Normalized-count table contains {invalid_counts} invalid values")

    summary = unique_index(summary_rows, "metric", "filtering summary")
    expected_summary = {
        "tximport_genes_before_mask": 118_082,
        "frozen_expanded_mask_genes": 16_112,
        "additional_coding_denominator_mask_genes": 292,
        "refined_mask_genes_total": 16_404,
        "additional_candidates_in_expanded_tested_universe": 197,
        "additional_candidates_not_in_expanded_tested_universe": 95,
        "genes_after_refined_mask_before_low_count_filter": 101_678,
        "genes_after_frozen_low_count_filter": 24_734,
        "refined_mask_host_focused_total_DEGs": 3_530,
        "refined_mask_host_focused_GP_up_DEGs": 2_147,
        "refined_mask_host_focused_UP_up_DEGs": 1_383,
        "refined_mask_host_focused_non_DE": 21_032,
        "refined_mask_host_focused_NA_padj_not_testable": 172,
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
            errors.append(f"Filtering summary mismatch for {metric}: {observed} != {expected}")

    r_md5_ok = True
    if set(identity) != set(R_ANALYSIS_INPUTS):
        errors.append("R input-identity membership differs from the five expected inputs")
        r_md5_ok = False
    else:
        for name, path in R_ANALYSIS_INPUTS.items():
            observed_md5 = identity[name]["MD5_before_analysis"].lower()
            if observed_md5 != digest(path, "md5"):
                errors.append(f"R input MD5 mismatch after analysis: {name}")
                r_md5_ok = False

    r_audit_lines = INPUTS["r_audit"].read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()
    expected_r_pass = "REFINED-MASK HOST-FOCUSED PLANT DESEQ2 FINAL ANALYSIS: PASS"
    if not r_audit_lines or r_audit_lines[-1] != expected_r_pass:
        errors.append("R analysis audit does not end in PASS")

    transitions: dict[str, Counter[str]] = defaultdict(Counter)
    changed_retained = 0
    for gene, old_status in expanded_calls.items():
        if gene in directly_removed:
            new_status = "removed_by_refined_coding_denominator_mask"
        else:
            new_status = refined_calls[gene]
            if new_status != old_status:
                changed_retained += 1
        transitions[old_status][new_status] += 1

        row = comparison[gene]
        try:
            flag = parse_bool(
                row["additional_coding_denominator_masked"], "comparison table"
            )
        except AuditError as error:
            errors.append(str(error))
            flag = False
        if flag != (gene in coding_mask):
            errors.append(f"Comparison mask flag mismatch: {gene}")
        expected_comparison = (
            "removed_by_refined_coding_denominator_mask"
            if gene in directly_removed
            else (
                "status_unchanged"
                if old_status == new_status
                else "status_changed_after_refit"
            )
        )
        if row["comparison_status"] != expected_comparison:
            errors.append(f"Comparison status mismatch: {gene}")
        if len(errors) >= 50:
            break

    contingency_rows: list[dict[str, object]] = []
    new_order = [*STATUS_ORDER, "removed_by_refined_coding_denominator_mask"]
    for old_status in STATUS_ORDER:
        for new_status in new_order:
            count = transitions[old_status][new_status]
            if count:
                contingency_rows.append(
                    {
                        "expanded_16112_mask_status": old_status,
                        "refined_16404_mask_status": new_status,
                        "genes": count,
                    }
                )
    write_tsv(
        OUTPUTS["status_contingency"],
        ["expanded_16112_mask_status", "refined_16404_mask_status", "genes"],
        contingency_rows,
    )

    if any(digest(path) != hashes_before[label] for label, path in INPUTS.items()):
        errors.append("One or more audited inputs changed during the audit")

    integrity_rows = [
        {"metric": "refined_mask_genes", "value": len(refined_mask)},
        {"metric": "expanded_tested_genes", "value": len(expanded)},
        {"metric": "refined_tested_genes", "value": len(refined)},
        {"metric": "directly_removed_tested_genes", "value": len(directly_removed)},
        {"metric": "retained_genes_changing_status", "value": changed_retained},
        {"metric": "expanded_total_DEGs", "value": expanded_counts["GP_up"] + expanded_counts["UP_up"]},
        {"metric": "refined_total_DEGs", "value": refined_counts["GP_up"] + refined_counts["UP_up"]},
        {"metric": "refined_GP_up", "value": refined_counts["GP_up"]},
        {"metric": "refined_UP_up", "value": refined_counts["UP_up"]},
        {"metric": "errors_detected", "value": len(errors)},
    ]
    write_tsv(OUTPUTS["integrity_summary"], ["metric", "value"], integrity_rows)

    report = [
        "===== REFINED-MASK HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT =====",
        f"Discovery root: {ROOT}",
        "",
        f"Refined mask genes: {len(refined_mask)}",
        f"Expanded-mask tested genes: {len(expanded)}",
        f"Refined-mask tested genes: {len(refined)}",
        f"Tested genes removed directly by coding-denominator additions: {len(directly_removed)}",
        "",
        "Expanded 16,112-mask statuses:",
        *[f"  {status}: {expanded_counts[status]}" for status in STATUS_ORDER],
        "",
        "Refined 16,404-mask statuses:",
        *[f"  {status}: {refined_counts[status]}" for status in STATUS_ORDER],
        "",
        "Direct removals by expanded-mask status:",
        *[f"  {status}: {removed_counts[status]}" for status in STATUS_ORDER],
        f"Retained genes changing status after refit: {changed_retained}",
        "",
        "Nonzero expanded-to-refined transitions:",
    ]
    for row in contingency_rows:
        report.append(
            f"  {row['expanded_16112_mask_status']} -> "
            f"{row['refined_16404_mask_status']}: {row['genes']}"
        )
    report.extend(
        [
            "",
            f"Masked genes retained in refined results: {len(set(refined) & refined_mask)}",
            "Result/count identifier agreement: " + ("PASS" if set(refined) == set(normalized) else "FAIL"),
            "R input MD5 verification: " + ("PASS" if r_md5_ok else "FAIL"),
            "Input SHA-256 verification after audit: "
            + ("PASS" if not any(digest(path) != hashes_before[label] for label, path in INPUTS.items()) else "FAIL"),
            f"Errors detected: {len(errors)}",
        ]
    )
    if errors:
        report.extend(["", "Errors:", *[f"  ERROR: {error}" for error in errors]])
        final_status = "CHECK REQUIRED"
    else:
        final_status = "PASS"
    report.extend(
        [
            "",
            "Interpretive boundary:",
            "  This audit verifies computation and membership only; it does not revise the frozen mask.",
            "  No sequence search, DESeq2 fitting, annotation search or reranking was performed.",
            f"REFINED-MASK HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT: {final_status}",
        ]
    )

    report_text = "\n".join(report) + "\n"
    write_text(OUTPUTS["audit"], report_text)
    print(report_text, end="")
    for path in OUTPUTS.values():
        print(f"Output: {path}")
    if errors:
        raise AuditError(f"Integrity audit detected {len(errors)} error(s)")


if __name__ == "__main__":
    try:
        main()
    except AuditError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        print(
            "REFINED-MASK HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT: CHECK REQUIRED",
            file=sys.stderr,
        )
        raise SystemExit(1)

