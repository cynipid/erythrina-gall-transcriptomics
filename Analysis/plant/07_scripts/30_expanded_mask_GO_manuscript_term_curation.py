#!/usr/bin/env python3
"""Build manuscript-ready current-resolved GO tables without retesting ORA.

The complete expanded-mask enrichment tables remain authoritative statistical
records. This script removes only GO rows whose identifiers were unresolved by
the GO.db version used for naming, retaining the original Fisher P-values and
BH-adjusted P-values unchanged in separate interpretation tables.
"""

from __future__ import annotations

import csv
import hashlib
import math
import sys
from pathlib import Path


ROOT = Path.home() / "projects" / "qe_cross_species_discovery"
BRANCH = ROOT / "05_host_focused_DESeq2_sensitivity"
ENRICHMENT = BRANCH / "enrichment"
LOGS = BRANCH / "logs"
CONTROL = ROOT / "00_project_control"

GP_ALL = ENRICHMENT / "host_focused_expanded_mask_GO_GP_up_all_terms_named.tsv"
UP_ALL = ENRICHMENT / "host_focused_expanded_mask_GO_UP_up_all_terms_named.tsv"
GP_SIG_SOURCE = ENRICHMENT / "host_focused_expanded_mask_GO_GP_up_FDR_0.05_named.tsv"
UP_SIG_SOURCE = ENRICHMENT / "host_focused_expanded_mask_GO_UP_up_FDR_0.05_named.tsv"
GO_REVIEW = ENRICHMENT / "host_focused_expanded_mask_unresolved_obsolete_GO_term_review.tsv"

OUTPUTS = {
    "gp_current": ENRICHMENT / "host_focused_expanded_mask_GO_GP_up_FDR_0.05_current_resolved.tsv",
    "up_current": ENRICHMENT / "host_focused_expanded_mask_GO_UP_up_FDR_0.05_current_resolved.tsv",
    "gp_top25": ENRICHMENT / "host_focused_expanded_mask_GO_GP_up_top25_current_resolved.tsv",
    "up_top25": ENRICHMENT / "host_focused_expanded_mask_GO_UP_up_top25_current_resolved.tsv",
    "summary": ENRICHMENT / "host_focused_expanded_mask_GO_manuscript_term_curation_summary.tsv",
    "control": CONTROL / "30_expanded_mask_GO_manuscript_term_curation_summary.tsv",
    "audit": LOGS / "30_expanded_mask_GO_manuscript_term_curation_audit.txt",
}

AUDIT_TITLE = "EXPANDED-MASK GO MANUSCRIPT TERM CURATION"
REQUIRED_COLUMNS = {
    "term_id", "GO_term", "ontology", "foreground_genes_with_term",
    "foreground_gene_set_size", "background_genes_with_term",
    "background_gene_set_size", "fold_enrichment", "odds_ratio",
    "pvalue", "padj", "significant_BH_0.05",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fnum(value: str) -> float:
    if value in {"", "NA", "NaN"}:
        return math.nan
    return float(value)


def row_signature(row: dict[str, str], fields: list[str]) -> tuple[str, ...]:
    return tuple(row.get(field, "") for field in fields)


def is_current_resolved(row: dict[str, str]) -> bool:
    return (
        row.get("GO_term", "") not in {"", "NA", "Unresolved in current GO.db"}
        and row.get("ontology", "") not in {"", "NA", "Unresolved"}
    )


def main() -> int:
    ENRICHMENT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    CONTROL.mkdir(parents=True, exist_ok=True)

    if OUTPUTS["audit"].is_file():
        previous = OUTPUTS["audit"].read_text(encoding="utf-8").splitlines()
        if previous and previous[-1] == f"{AUDIT_TITLE}: PASS":
            print("Validated GO manuscript-curation outputs already exist; no files were changed.")
            return 1

    errors: list[str] = []
    required = [GP_ALL, UP_ALL, GP_SIG_SOURCE, UP_SIG_SOURCE, GO_REVIEW]
    missing = [path for path in required if not path.is_file()]
    if missing:
        errors.extend(f"Missing required file: {path}" for path in missing)

    summary_rows: list[dict[str, object]] = []
    output_rows: dict[str, list[dict[str, str]]] = {}
    source_hashes: list[str] = []

    if not errors:
        source_hashes = [sha256(path) for path in required]
        review_rows = read_tsv(GO_REVIEW)
        unresolved_review_ids = {row.get("term_id", "") for row in review_rows}
        unresolved_review_ids.discard("")
        if len(unresolved_review_ids) != 684:
            errors.append(
                f"Unresolved GO review contains {len(unresolved_review_ids)} IDs; expected 684"
            )

        configurations = [
            ("GP_up", GP_ALL, GP_SIG_SOURCE, "gp_current", "gp_top25", 64, 9, 55),
            ("UP_up", UP_ALL, UP_SIG_SOURCE, "up_current", "up_top25", 34, 6, 28),
        ]

        for direction, all_path, sig_path, current_key, top_key, expected_sig, expected_unresolved_sig, expected_current in configurations:
            all_rows = read_tsv(all_path)
            sig_source_rows = read_tsv(sig_path)
            if len(all_rows) != 5667:
                errors.append(f"{direction} all-term table has {len(all_rows)} rows; expected 5667")
            if all_rows and not REQUIRED_COLUMNS.issubset(all_rows[0]):
                errors.append(
                    f"{direction} all-term table lacks columns: "
                    + ", ".join(sorted(REQUIRED_COLUMNS - set(all_rows[0])))
                )
                continue

            significant_rows = [row for row in all_rows if fnum(row.get("padj", "NA")) < 0.05]
            if len(significant_rows) != expected_sig:
                errors.append(
                    f"{direction} significant rows reconstructed from all terms: "
                    f"{len(significant_rows)}; expected {expected_sig}"
                )

            fields = list(all_rows[0]) if all_rows else []
            reconstructed_signatures = [row_signature(row, fields) for row in significant_rows]
            source_signatures = [row_signature(row, fields) for row in sig_source_rows]
            if reconstructed_signatures != source_signatures:
                errors.append(f"{direction} stored significant table differs from the all-term source")

            unresolved_all_ids = {
                row["term_id"] for row in all_rows if not is_current_resolved(row)
            }
            if unresolved_all_ids != unresolved_review_ids:
                errors.append(
                    f"{direction} unresolved-ID membership differs from the 684-ID review inventory"
                )

            unresolved_significant = [
                row for row in significant_rows if not is_current_resolved(row)
            ]
            current_significant = [
                row for row in significant_rows if is_current_resolved(row)
            ]
            if len(unresolved_significant) != expected_unresolved_sig:
                errors.append(
                    f"{direction} significant unresolved terms: {len(unresolved_significant)}; "
                    f"expected {expected_unresolved_sig}"
                )
            if len(current_significant) != expected_current:
                errors.append(
                    f"{direction} current-resolved significant terms: {len(current_significant)}; "
                    f"expected {expected_current}"
                )
            if len(current_significant) + len(unresolved_significant) != len(significant_rows):
                errors.append(f"{direction} significant-term partition is incomplete")

            # The curation is post-statistical: every retained row must remain
            # byte-equivalent at field level to its source row.
            source_by_id = {row["term_id"]: row for row in significant_rows}
            for row in current_significant:
                if row_signature(row, fields) != row_signature(source_by_id[row["term_id"]], fields):
                    errors.append(f"{direction} statistical field changed for {row['term_id']}")
                    break

            output_rows[current_key] = current_significant
            output_rows[top_key] = current_significant[:25]
            summary_rows.append({
                "direction": direction,
                "tested_GO_categories_full_statistical_table": len(all_rows),
                "significant_GO_categories_full_statistical_table": len(significant_rows),
                "significant_unresolved_GO_categories_excluded_from_interpretation": len(unresolved_significant),
                "significant_current_resolved_GO_categories_for_interpretation": len(current_significant),
                "statistical_recalculation_performed": "FALSE",
                "BH_values_changed": "FALSE",
                "interpretation_policy": "use current-resolved table for manuscript text and figures; retain full table for reproducibility",
            })

        if not errors:
            table_fields = list(read_tsv(GP_ALL)[0])
            write_tsv(OUTPUTS["gp_current"], output_rows["gp_current"], table_fields)
            write_tsv(OUTPUTS["up_current"], output_rows["up_current"], table_fields)
            write_tsv(OUTPUTS["gp_top25"], output_rows["gp_top25"], table_fields)
            write_tsv(OUTPUTS["up_top25"], output_rows["up_top25"], table_fields)
            summary_fields = list(summary_rows[0])
            write_tsv(OUTPUTS["summary"], summary_rows, summary_fields)
            write_tsv(OUTPUTS["control"], summary_rows, summary_fields)

            for key in ("gp_current", "up_current", "gp_top25", "up_top25", "summary", "control"):
                if not OUTPUTS[key].is_file():
                    errors.append(f"Expected curation output missing: {OUTPUTS[key]}")

            if source_hashes != [sha256(path) for path in required]:
                errors.append("One or more source tables changed during curation")

    lines = [
        f"===== {AUDIT_TITLE} =====",
        "Scope: post-statistical curation only; no ORA or BH recalculation.",
        "Full statistical enrichment tables were preserved unchanged.",
    ]
    if summary_rows:
        for row in summary_rows:
            lines.extend([
                "",
                f"{row['direction']}:",
                f"  Tested GO categories in full table: {row['tested_GO_categories_full_statistical_table']}",
                f"  Significant GO categories in full table: {row['significant_GO_categories_full_statistical_table']}",
                "  Significant unresolved categories excluded from interpretation: "
                f"{row['significant_unresolved_GO_categories_excluded_from_interpretation']}",
                "  Significant current-resolved categories retained for interpretation: "
                f"{row['significant_current_resolved_GO_categories_for_interpretation']}",
            ])
    lines.extend([
        "",
        "Policy: report full enrichment accounting for reproducibility and use only current-resolved GO labels in manuscript text and figures.",
        "Original Fisher P-values and BH-adjusted P-values retained unchanged: PASS" if not errors else
        "Original Fisher P-values and BH-adjusted P-values retained unchanged: CHECK REQUIRED",
        f"Errors detected: {len(errors)}",
    ])
    lines.extend(f"ERROR: {error}" for error in errors)
    lines.append(f"{AUDIT_TITLE}: {'PASS' if not errors else 'CHECK REQUIRED'}")
    OUTPUTS["audit"].write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    if not errors:
        for key, path in OUTPUTS.items():
            print(f"Output: {path}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

