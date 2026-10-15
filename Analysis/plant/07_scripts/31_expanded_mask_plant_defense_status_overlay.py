#!/usr/bin/env python3
"""Overlay final expanded-mask DESeq2 status on frozen plant-defense candidates.

Candidate membership and family assignments are imported unchanged from the
validated gene-level defense inventory. No annotation search, family redefinition,
MEROPS search, SignalP/DeepTMHMM run, or differential-expression test is performed.
"""

from __future__ import annotations

import csv
import hashlib
import math
import sys
from collections import Counter
from pathlib import Path


ROOT = Path.home() / "projects" / "qe_cross_species_discovery"
DEFENSE_ROOT = ROOT / "05_plant_defense_audit"
RESULTS = DEFENSE_ROOT / "results"
LOGS = ROOT / "05_host_focused_DESeq2_sensitivity" / "logs"
CONTROL = ROOT / "00_project_control"

FROZEN_CANDIDATES = RESULTS / "plant_host_focused_defense_candidate_gene_summary.tsv"
EXPANDED_DESEQ = (
    ROOT / "05_host_focused_DESeq2_sensitivity" / "results"
    / "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_eggnog_joined.tsv"
)
EXPANDED_MASK = (
    ROOT / "04_cross_species_integration" / "results"
    / "plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv"
)

OUTPUTS = {
    "gene_summary": RESULTS / "plant_expanded_mask_defense_candidate_gene_summary.tsv",
    "status_summary": RESULTS / "plant_expanded_mask_defense_unique_gene_status_summary.tsv",
    "gp_candidates": RESULTS / "plant_expanded_mask_defense_GP_up_candidates.tsv",
    "up_candidates": RESULTS / "plant_expanded_mask_defense_UP_up_candidates.tsv",
    "transitions": RESULTS / "plant_defense_current_to_expanded_mask_status_transitions.tsv",
    "parameters": CONTROL / "31_expanded_mask_plant_defense_status_overlay_parameters.tsv",
    "audit": DEFENSE_ROOT / "31_expanded_mask_plant_defense_status_overlay_audit.txt",
}

AUDIT_TITLE = "EXPANDED-MASK PLANT-DEFENSE STATUS OVERLAY"
STATUS_COLUMN = "DESeq2_status_host_focused_expanded_mask"
MASK_EXCLUDE_STATUSES = {
    "review_exclude_existing_confirmed_multilevel_sequence",
    "review_exclude_additional_strong_nucleotide_similarity",
}
FINAL_STATUS_ORDER = [
    "GP_up", "UP_up", "non_DE", "NA_padj_not_testable",
    "excluded_by_expanded_nucleotide_mask",
    "not_in_expanded_mask_host_focused_tested_universe",
]
DEFENSE_GROUPS = [
    "PR_protein", "jasmonate", "salicylate",
    "peptidase", "peptidase_inhibitor",
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


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


def fnum(value: str | None) -> float:
    if value is None or value.strip() in {"", "NA", "NaN"}:
        return math.nan
    return float(value)


def index_unique(rows: list[dict[str, str]], key: str, label: str, errors: list[str]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        value = row.get(key, "")
        if not value:
            errors.append(f"{label} contains a blank {key}")
            continue
        if value in result:
            errors.append(f"{label} contains duplicate {key}: {value}")
            continue
        result[value] = row
    return result


def split_groups(value: str) -> set[str]:
    return {part.strip() for part in value.split(";") if part.strip()}


def status_rank(status: str) -> int:
    try:
        return FINAL_STATUS_ORDER.index(status)
    except ValueError:
        return len(FINAL_STATUS_ORDER)


def numeric_sort_value(value: str | None, default: float) -> float:
    parsed = fnum(value)
    return default if math.isnan(parsed) else parsed


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    CONTROL.mkdir(parents=True, exist_ok=True)

    if OUTPUTS["audit"].is_file():
        previous = OUTPUTS["audit"].read_text(encoding="utf-8").splitlines()
        if previous and previous[-1] == f"{AUDIT_TITLE}: PASS":
            print("Validated expanded-mask defense-overlay outputs already exist; no files were changed.")
            return 1

    errors: list[str] = []
    required = [FROZEN_CANDIDATES, EXPANDED_DESEQ, EXPANDED_MASK]
    missing = [path for path in required if not path.is_file()]
    if missing:
        errors.extend(f"Missing required file: {path}" for path in missing)

    output_rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    transition_rows: list[dict[str, object]] = []
    input_hashes_before: list[str] = []

    if not errors:
        input_hashes_before = [sha256(path) for path in required]
        frozen_fields, frozen_rows = read_tsv(FROZEN_CANDIDATES)
        deseq_fields, deseq_rows = read_tsv(EXPANDED_DESEQ)
        mask_fields, mask_rows = read_tsv(EXPANDED_MASK)

        required_frozen = {
            "gene_id", "defense_groups", "candidate_subgroups",
            "host_focused_status", "host_baseMean", "host_log2FoldChange",
            "host_padj", "signalp_prediction", "has_TMR",
            "primary_classical_secretome", "MEROPS_overlay_present",
        }
        required_deseq = {
            "gene_id", "baseMean", "log2FoldChange", "padj", STATUS_COLUMN,
        }
        required_mask = {"gene_id", "proposed_expanded_mask_status"}
        if not required_frozen.issubset(frozen_fields):
            errors.append(
                "Frozen candidate table lacks columns: "
                + ", ".join(sorted(required_frozen - set(frozen_fields)))
            )
        if not required_deseq.issubset(deseq_fields):
            errors.append(
                "Expanded DESeq2 table lacks columns: "
                + ", ".join(sorted(required_deseq - set(deseq_fields)))
            )
        if not required_mask.issubset(mask_fields):
            errors.append(
                "Expanded-mask table lacks columns: "
                + ", ".join(sorted(required_mask - set(mask_fields)))
            )

        frozen = index_unique(frozen_rows, "gene_id", "frozen defense candidates", errors)
        deseq = index_unique(deseq_rows, "gene_id", "expanded-mask DESeq2", errors)
        mask = index_unique(mask_rows, "gene_id", "expanded nucleotide mask", errors)

        # The final corrected candidate inventory used by the discovery overlay
        # contains 806 unique genes. Stop if an earlier 817-gene draft is supplied.
        if len(frozen) != 806:
            errors.append(
                f"Frozen defense candidate genes: {len(frozen)}; expected final corrected inventory of 806"
            )
        if len(deseq) != 24931:
            errors.append(f"Expanded-mask DESeq2 genes: {len(deseq)}; expected 24931")
        if len(mask) != 118082:
            errors.append(f"Expanded-mask universe: {len(mask)}; expected 118082")

        expanded_mask_ids = {
            gene for gene, row in mask.items()
            if row.get("proposed_expanded_mask_status") in MASK_EXCLUDE_STATUSES
        }
        if len(expanded_mask_ids) != 16112:
            errors.append(f"Expanded exclusion mask genes: {len(expanded_mask_ids)}; expected 16112")
        if set(deseq) & expanded_mask_ids:
            errors.append("Expanded-mask genes are unexpectedly present in final DESeq2 results")

        expected_deseq_counts = Counter({
            "GP_up": 2349, "UP_up": 1386,
            "non_DE": 21022, "NA_padj_not_testable": 174,
        })
        observed_deseq_counts = Counter(row.get(STATUS_COLUMN, "") for row in deseq_rows)
        if observed_deseq_counts != expected_deseq_counts:
            errors.append(f"Expanded-mask DESeq2 status counts differ: {dict(observed_deseq_counts)}")

        allowed_groups = set(DEFENSE_GROUPS)
        observed_groups: set[str] = set()
        for row in frozen_rows:
            observed_groups.update(split_groups(row.get("defense_groups", "")))
        if observed_groups != allowed_groups:
            errors.append(
                "Frozen defense-group membership differs: "
                f"observed {sorted(observed_groups)}, expected {sorted(allowed_groups)}"
            )

        appended_fields = [
            "candidate_membership_source",
            "candidate_membership_changed",
            "previous_3948_mask_host_status",
            "expanded_mask_host_status",
            "expanded_mask_host_baseMean",
            "expanded_mask_host_log2FoldChange",
            "expanded_mask_host_padj",
            "expanded_nucleotide_mask_membership",
        ]
        output_fields = [*frozen_fields, *[field for field in appended_fields if field not in frozen_fields]]

        for gene, source in frozen.items():
            updated: dict[str, object] = dict(source)
            expanded = deseq.get(gene)
            in_mask = gene in expanded_mask_ids
            if expanded is not None:
                final_status = expanded[STATUS_COLUMN]
                final_base_mean = expanded.get("baseMean", "")
                final_lfc = expanded.get("log2FoldChange", "")
                final_padj = expanded.get("padj", "")
            elif in_mask:
                final_status = "excluded_by_expanded_nucleotide_mask"
                final_base_mean = "NA"
                final_lfc = "NA"
                final_padj = "NA"
            else:
                final_status = "not_in_expanded_mask_host_focused_tested_universe"
                final_base_mean = "NA"
                final_lfc = "NA"
                final_padj = "NA"

            updated.update({
                "candidate_membership_source": "frozen_final_corrected_host_focused_defense_gene_summary",
                "candidate_membership_changed": "FALSE",
                "previous_3948_mask_host_status": source.get("host_focused_status", ""),
                "expanded_mask_host_status": final_status,
                "expanded_mask_host_baseMean": final_base_mean,
                "expanded_mask_host_log2FoldChange": final_lfc,
                "expanded_mask_host_padj": final_padj,
                "expanded_nucleotide_mask_membership": "TRUE" if in_mask else "FALSE",
            })
            output_rows.append(updated)

        output_rows.sort(key=lambda row: (
            status_rank(str(row["expanded_mask_host_status"])),
            numeric_sort_value(str(row["expanded_mask_host_padj"]), 1.0),
            -abs(numeric_sort_value(str(row["expanded_mask_host_log2FoldChange"]), 0.0)),
            str(row["gene_id"]),
        ))

        if {str(row["gene_id"]) for row in output_rows} != set(frozen):
            errors.append("Output candidate membership differs from the frozen 806-gene inventory")
        if any(row["candidate_membership_changed"] != "FALSE" for row in output_rows):
            errors.append("One or more candidates are marked as having changed membership")

        # Verify that every original membership/annotation field is preserved.
        for row in output_rows:
            source = frozen[str(row["gene_id"])]
            for field in frozen_fields:
                if str(row.get(field, "")) != source.get(field, ""):
                    errors.append(f"Frozen field {field} changed for {row['gene_id']}")
                    break
            if errors:
                break

        for group in DEFENSE_GROUPS:
            for status in FINAL_STATUS_ORDER:
                count = sum(
                    group in split_groups(str(row["defense_groups"]))
                    and row["expanded_mask_host_status"] == status
                    for row in output_rows
                )
                status_rows.append({
                    "defense_group": group,
                    "expanded_mask_host_status": status,
                    "unique_candidate_genes": count,
                })

        transitions = Counter(
            (str(row["previous_3948_mask_host_status"]), str(row["expanded_mask_host_status"]))
            for row in output_rows
        )
        for (previous, current), count in sorted(
            transitions.items(), key=lambda item: (item[0][0], status_rank(item[0][1]))
        ):
            transition_rows.append({
                "previous_3948_mask_host_status": previous,
                "expanded_16112_mask_host_status": current,
                "candidate_genes": count,
            })

        if not errors:
            write_tsv(OUTPUTS["gene_summary"], output_rows, output_fields)
            write_tsv(OUTPUTS["status_summary"], status_rows, list(status_rows[0]))
            write_tsv(
                OUTPUTS["gp_candidates"],
                [row for row in output_rows if row["expanded_mask_host_status"] == "GP_up"],
                output_fields,
            )
            write_tsv(
                OUTPUTS["up_candidates"],
                [row for row in output_rows if row["expanded_mask_host_status"] == "UP_up"],
                output_fields,
            )
            write_tsv(OUTPUTS["transitions"], transition_rows, list(transition_rows[0]))

            parameters = [
                {"input_role": "frozen_defense_candidate_membership", "path": str(FROZEN_CANDIDATES), "sha256": input_hashes_before[0]},
                {"input_role": "expanded_mask_DESeq2", "path": str(EXPANDED_DESEQ), "sha256": input_hashes_before[1]},
                {"input_role": "expanded_nucleotide_mask", "path": str(EXPANDED_MASK), "sha256": input_hashes_before[2]},
            ]
            write_tsv(OUTPUTS["parameters"], parameters, list(parameters[0]))

            if input_hashes_before != [sha256(path) for path in required]:
                errors.append("One or more inputs changed during status-overlay generation")

    status_counts = Counter(str(row.get("expanded_mask_host_status", "")) for row in output_rows)
    previous_counts = Counter(str(row.get("previous_3948_mask_host_status", "")) for row in output_rows)
    lines = [
        f"===== {AUDIT_TITLE} =====",
        "Scope: frozen candidate membership with final expanded-mask DESeq2 status overlay.",
        "No candidate-family search, MEROPS search, SignalP/DeepTMHMM run or DESeq2 analysis was performed.",
        "",
        f"Frozen candidate genes: {len(output_rows)}",
        "Candidate membership changed: 0",
        f"Candidates overlapping expanded nucleotide mask: {status_counts.get('excluded_by_expanded_nucleotide_mask', 0)}",
        "",
        "Previous 3,948-mask candidate statuses:",
    ]
    for status in sorted(previous_counts):
        lines.append(f"  {status}: {previous_counts[status]}")
    lines.extend(["", "Final expanded-mask candidate statuses:"])
    for status in FINAL_STATUS_ORDER:
        lines.append(f"  {status}: {status_counts.get(status, 0)}")
    lines.extend([
        "",
        "MEROPS boundary: the retained overlay came from the original GP-up classical-secretome-restricted screen; absence remains non-negative evidence.",
        "Candidate membership and original annotation fields preserved: " + ("PASS" if not errors else "CHECK REQUIRED"),
        f"Errors detected: {len(errors)}",
    ])
    lines.extend(f"ERROR: {error}" for error in errors)
    lines.append(f"{AUDIT_TITLE}: {'PASS' if not errors else 'CHECK REQUIRED'}")
    OUTPUTS["audit"].write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    if not errors:
        for path in OUTPUTS.values():
            print(f"Output: {path}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

