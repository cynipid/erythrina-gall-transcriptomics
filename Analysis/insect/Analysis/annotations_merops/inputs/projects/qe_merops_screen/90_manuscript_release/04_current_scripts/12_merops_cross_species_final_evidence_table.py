#!/usr/bin/env python3
"""Create a manuscript-ready cross-species MEROPS evidence table.

The table juxtaposes, but does NOT pair, the final candidate sets:
  - 20 priority plant GP-up classically secreted peptidase-related candidates
  - 9 larva-up classically secreted insect inhibitor-family candidates

No row represents a demonstrated plant--insect protein interaction.

Run in Linux:
    python3 "$HOME/projects/qe_merops_screen/scripts/12_merops_cross_species_final_evidence_table.py"
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


HOME = Path.home()
BASE = HOME / "projects" / "qe_merops_screen"
RESULTS = BASE / "results"

PLANT_INPUT = RESULTS / "plant_GPup_MEROPS_peptidase_priority_curated_v2.tsv"
INSECT_INPUT = RESULTS / "insect_larva_up_inhibitor_annotation_domain_audit_v2.tsv"
PUPA_FAMILY_INPUT = RESULTS / "insect_pupa_up_secreted_inhibitor_family_summary_curated_v2.tsv"

TABLE_OUT = RESULTS / "MEROPS_cross_species_final_evidence_table_v2.tsv"
SUMMARY_OUT = RESULTS / "MEROPS_cross_species_final_evidence_summary_v2.tsv"
AUDIT_OUT = RESULTS / "MEROPS_cross_species_final_evidence_audit_v2.txt"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def get(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name, "")
        if value not in ("", "-", "NA"):
            return value
    return ""


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def main() -> None:
    for path in (PLANT_INPUT, INSECT_INPUT, PUPA_FAMILY_INPUT):
        if not path.exists():
            raise FileNotFoundError(f"Missing required input: {path}")

    plant_all = read_tsv(PLANT_INPUT)
    insect = read_tsv(INSECT_INPUT)
    pupa_families = read_tsv(PUPA_FAMILY_INPUT)

    required_plant = {"gene_id", "merops_family", "manuscript_priority"}
    required_insect = {"gene_id", "merops_family", "domain_evidence_status"}
    missing_plant = required_plant.difference(plant_all[0] if plant_all else {})
    missing_insect = required_insect.difference(insect[0] if insect else {})
    if missing_plant:
        raise RuntimeError("Plant input missing: " + ", ".join(sorted(missing_plant)))
    if missing_insect:
        raise RuntimeError("Insect input missing: " + ", ".join(sorted(missing_insect)))

    # The input retains all 98 plant homology/review rows.  The manuscript
    # evidence table intentionally includes only the pre-defined Priority 1
    # and Priority 2 subset, not unresolved or non-peptidase review entries.
    plant = [
        row for row in plant_all
        if row["manuscript_priority"].startswith(("priority_1_", "priority_2_"))
    ]

    rows: list[dict[str, object]] = []
    for row in plant:
        priority = row["manuscript_priority"]
        evidence = (
            "Plant seed, named MEROPS reference, and annotation support"
            if priority.startswith("priority_1")
            else "Plant seed and annotation support; unassigned MEROPS-family reference"
        )
        rows.append(
            {
                "organism": "Erythrina variegata",
                "comparison_and_direction": "Galled petiole higher (GP-up)",
                "candidate_role": "Classically secreted peptidase-related candidate",
                "selection_tier": priority,
                "gene_id": row["gene_id"],
                "merops_identifier": get(row, "merops_identifier"),
                "merops_family": row["merops_family"],
                "merops_reference_organism": get(row, "merops_reference_organism"),
                "merops_reference_status": get(row, "merops_reference_status"),
                "family_annotation": get(row, "Description", "eggNOG_Description"),
                "preferred_name": get(row, "Preferred_name", "eggNOG_Preferred_name"),
                "PFAMs": get(row, "PFAMs", "eggNOG_PFAMs"),
                "baseMean": get(row, "baseMean"),
                "log2FoldChange": get(row, "log2FoldChange"),
                "padj": get(row, "padj"),
                "secretion_class": get(row, "analysis_class"),
                "independent_domain_or_annotation_evidence": evidence,
                "interpretation_boundary": (
                    "Candidate protease-related protein; catalytic activity, extracellular "
                    "location in gall tissue, and interaction with an insect protein are untested"
                ),
            }
        )

    for row in insect:
        rows.append(
            {
                "organism": "Quadrastichus erythrinae",
                "comparison_and_direction": "Larva higher (larva-up)",
                "candidate_role": "Classically secreted inhibitor-family candidate",
                "selection_tier": "larva_up_family_and_PFAM_supported",
                "gene_id": row["gene_id"],
                "merops_identifier": get(row, "merops_identifier"),
                "merops_family": row["merops_family"],
                "merops_reference_organism": get(row, "merops_reference_organism"),
                "merops_reference_status": get(row, "merops_reference_status"),
                "family_annotation": get(row, "eggNOG_Description", "Description"),
                "preferred_name": get(row, "eggNOG_Preferred_name", "Preferred_name"),
                "PFAMs": get(row, "eggNOG_PFAMs", "PFAMs"),
                "baseMean": get(row, "baseMean"),
                "log2FoldChange": get(row, "log2FoldChange"),
                "padj": get(row, "padj"),
                "secretion_class": get(row, "analysis_class"),
                "independent_domain_or_annotation_evidence": get(row, "domain_evidence_status"),
                "interpretation_boundary": get(row, "interpretation_boundary"),
            }
        )

    rows.sort(
        key=lambda row: (
            0 if row["organism"] == "Erythrina variegata" else 1,
            -abs(as_float(str(row["log2FoldChange"]))),
            str(row["gene_id"]),
        )
    )

    fields = [
        "organism", "comparison_and_direction", "candidate_role", "selection_tier",
        "gene_id", "merops_identifier", "merops_family",
        "merops_reference_organism", "merops_reference_status", "family_annotation",
        "preferred_name", "PFAMs", "baseMean", "log2FoldChange", "padj",
        "secretion_class", "independent_domain_or_annotation_evidence",
        "interpretation_boundary",
    ]
    write_tsv(TABLE_OUT, rows, fields)

    plant_priority_counts = Counter(row["manuscript_priority"] for row in plant)
    insect_family_counts = Counter(row["merops_family"] for row in insect)
    pupa_family_counts = {
        row["merops_family"]: row["pupa_up_candidate_genes"]
        for row in pupa_families
    }

    summary = [
        {"section": "plant", "metric": "priority_peptidase_related_candidates", "count": len(plant)},
        {"section": "plant", "metric": "priority_1_candidates", "count": plant_priority_counts["priority_1_named_reference_plus_plant_seed_plus_annotation"]},
        {"section": "plant", "metric": "priority_2_candidates", "count": plant_priority_counts["priority_2_unassigned_reference_plus_plant_seed_plus_annotation"]},
        {"section": "insect_larva", "metric": "PFAM_supported_inhibitor_candidates", "count": len(insect)},
        {"section": "insect_larva", "metric": "I01_Kazal_type_candidates", "count": insect_family_counts["I01"]},
        {"section": "insect_larva", "metric": "I04_serpin_candidates", "count": insect_family_counts["I04"]},
        {"section": "insect_larva", "metric": "I19_pacifastin_candidates", "count": insect_family_counts["I19"]},
        {"section": "insect_pupa_context", "metric": "all_secreted_inhibitor_candidates", "count": 68},
        {"section": "insect_pupa_context", "metric": "I43_candidates", "count": pupa_family_counts.get("I43", "")},
        {"section": "insect_pupa_context", "metric": "I63_candidates", "count": pupa_family_counts.get("I63", "")},
        {"section": "insect_pupa_context", "metric": "I08_candidates", "count": pupa_family_counts.get("I08", "")},
    ]
    write_tsv(SUMMARY_OUT, summary, ["section", "metric", "count"])

    errors: list[str] = []
    if len(plant) != 20:
        errors.append("Expected 20 plant priority candidates")
    if len(insect) != 9:
        errors.append("Expected 9 insect larva-up candidates")
    if len({row["gene_id"] for row in plant}) != len(plant):
        errors.append("Plant priority table has duplicate gene IDs")
    if len({row["gene_id"] for row in insect}) != len(insect):
        errors.append("Insect larva table has duplicate gene IDs")
    if set(insect_family_counts) != {"I01", "I04", "I19"}:
        errors.append("Insect larva table has unexpected MEROPS-family membership")
    if pupa_family_counts.get("I43") != "28" or pupa_family_counts.get("I63") != "19" or pupa_family_counts.get("I08") != "10":
        errors.append("Pupa I43/I63/I08 context counts do not match frozen values")

    report = [
        "===== MEROPS CROSS-SPECIES FINAL EVIDENCE TABLE V2 =====",
        f"Plant GP-up priority candidates: {len(plant)}",
        f"Insect larva-up Pfam-supported inhibitor candidates: {len(insect)}",
        f"Pupa-up inhibitor context: I43={pupa_family_counts.get('I43', '')}, "
        f"I63={pupa_family_counts.get('I63', '')}, I08={pupa_family_counts.get('I08', '')}",
        "",
        "Interpretive boundary: this table juxtaposes candidate programmes; it does not pair plant and insect proteins or demonstrate a molecular interaction.",
        "",
        f"Errors detected: {len(errors)}",
        *[f"ERROR: {error}" for error in errors],
        (
            "MEROPS CROSS-SPECIES FINAL EVIDENCE TABLE V2: CHECK REQUIRED"
            if errors else "MEROPS CROSS-SPECIES FINAL EVIDENCE TABLE V2: PASS"
        ),
    ]
    AUDIT_OUT.write_text("\n".join(report) + "\n", encoding="utf-8")

    print("\n".join(report))
    print(f"\nOutput: {TABLE_OUT}")
    print(f"Output: {SUMMARY_OUT}")
    print(f"Output: {AUDIT_OUT}")
    if errors:
        raise SystemExit("MEROPS cross-species final evidence table failed")


if __name__ == "__main__":
    main()
