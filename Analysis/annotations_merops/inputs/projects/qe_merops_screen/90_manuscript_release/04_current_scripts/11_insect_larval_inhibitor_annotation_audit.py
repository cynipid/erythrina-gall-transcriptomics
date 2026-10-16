#!/usr/bin/env python3
"""Final annotation/domain-evidence audit for nine larva-up inhibitor candidates.

This script reuses existing eggNOG annotation. It does not perform a new
genome-wide search, does not infer biochemical activity, and does not identify
plant targets. A missing expected Pfam domain is reported for review, not
treated as a failure, because eggNOG annotation is incomplete for some genes.

Run in Linux:
    python3 "$HOME/projects/qe_merops_screen/11_insect_larval_inhibitor_annotation_audit.py"
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


HOME = Path.home()
MEROPS = HOME / "projects" / "qe_merops_screen"
INSECT = HOME / "projects" / "qe_genome_guided"
RESULTS = MEROPS / "results"

SHORTLIST_FILE = RESULTS / "insect_larva_up_I01_I04_I19_inhibitor_shortlist_v2.tsv"
EGGNOG_FILE = (
    INSECT / "10_functional_annotation" / "eggnog_out" /
    "Qerythrinae.emapper.annotations"
)

DETAIL_OUT = RESULTS / "insect_larva_up_inhibitor_annotation_domain_audit_v2.tsv"
SUMMARY_OUT = RESULTS / "insect_larva_up_inhibitor_annotation_domain_summary_v2.tsv"
AUDIT_OUT = RESULTS / "insect_larva_up_inhibitor_annotation_domain_audit_v2.txt"

EXPECTED_FAMILY_DOMAIN = {
    "I01": ("Kazal", re.compile(r"kazal", re.IGNORECASE)),
    "I04": ("Serpin", re.compile(r"serpin", re.IGNORECASE)),
    "I19": ("Pacifastin", re.compile(r"pacifastin", re.IGNORECASE)),
}
EXPECTED_COUNTS = {"I01": 2, "I04": 4, "I19": 3}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_emapper_annotations(path: Path) -> tuple[dict[str, dict[str, str]], str]:
    """Read eggNOG's commented header safely, preserving rows with extra tabs."""
    header: list[str] | None = None
    annotations: dict[str, dict[str, str]] = {}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line.startswith("#query\t"):
                header = line[1:].split("\t")
                continue
            if not line or line.startswith("#"):
                continue
            if header is None:
                raise RuntimeError("eggNOG annotation file has no #query header")

            values = line.split("\t")
            if len(values) < len(header):
                values.extend([""] * (len(header) - len(values)))
            elif len(values) > len(header):
                # Descriptive text occasionally contains an extra tab. Preserve
                # it in the final field rather than shifting every annotation.
                values = values[: len(header) - 1] + ["\t".join(values[len(header) - 1 :])]

            row = dict(zip(header, values))
            query = row.get("query", "").strip()
            if query:
                annotations[query] = row

    if header is None:
        raise RuntimeError("Could not find a #query header in eggNOG annotation file")
    return annotations, "\t".join(header)


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def main() -> None:
    for path in (SHORTLIST_FILE, EGGNOG_FILE):
        if not path.exists():
            raise FileNotFoundError(f"Missing required input: {path}")

    shortlist = read_tsv(SHORTLIST_FILE)
    if not shortlist:
        raise RuntimeError("Larva-up shortlist is empty")

    required = {"gene_id", "merops_family", "DESeq2_status", "analysis_class"}
    missing = required.difference(shortlist[0])
    if missing:
        raise RuntimeError("Shortlist missing columns: " + ", ".join(sorted(missing)))

    eggnog, eggnog_header = parse_emapper_annotations(EGGNOG_FILE)

    detail: list[dict[str, object]] = []
    for candidate in shortlist:
        gene_id = candidate["gene_id"]
        family = candidate["merops_family"]
        annotation = eggnog.get(gene_id, {})
        pfams = annotation.get("PFAMs", "")
        expected_domain, domain_regex = EXPECTED_FAMILY_DOMAIN.get(
            family, ("not_predefined", re.compile(r"$^"))
        )
        annotation_text = " | ".join(
            [
                pfams,
                annotation.get("Description", ""),
                annotation.get("eggNOG_Description", ""),
                annotation.get("Preferred_name", ""),
            ]
        )
        if not annotation:
            domain_status = "eggNOG_annotation_not_available"
        elif domain_regex.search(pfams):
            domain_status = "expected_domain_detected_in_PFAMs"
        elif domain_regex.search(annotation_text):
            domain_status = "expected_domain_named_outside_PFAMs"
        else:
            domain_status = "expected_domain_not_detected_in_eggNOG"

        row: dict[str, object] = dict(candidate)
        row.update(
            {
                "expected_family_domain": expected_domain,
                "domain_evidence_status": domain_status,
                "eggNOG_annotation_available": "yes" if annotation else "no",
                "eggNOG_query": annotation.get("query", ""),
                "eggNOG_seed_ortholog": annotation.get("seed_ortholog", ""),
                "eggNOG_Description": annotation.get("Description", ""),
                "eggNOG_Preferred_name": annotation.get("Preferred_name", ""),
                "eggNOG_GOs": annotation.get("GOs", ""),
                "eggNOG_KEGG_ko": annotation.get("KEGG_ko", ""),
                "eggNOG_CAZy": annotation.get("CAZy", ""),
                "eggNOG_PFAMs": pfams,
                "eggNOG_OGs": annotation.get("eggNOG_OGs", ""),
                "interpretation_boundary": (
                    "Annotation supports family-level homology only; it does not prove "
                    "inhibitory activity, plant delivery, or a plant target"
                ),
            }
        )
        detail.append(row)

    detail.sort(
        key=lambda row: (
            ["I01", "I04", "I19"].index(str(row["merops_family"])),
            as_float(str(row.get("padj", "nan"))),
            str(row["gene_id"]),
        )
    )

    family_counts = Counter(str(row["merops_family"]) for row in detail)
    status_by_family: dict[str, Counter[str]] = defaultdict(Counter)
    for row in detail:
        status_by_family[str(row["merops_family"])][str(row["domain_evidence_status"])] += 1

    summary: list[dict[str, object]] = []
    for family in ("I01", "I04", "I19"):
        status_counts = status_by_family[family]
        summary.append(
            {
                "merops_family": family,
                "expected_family_domain": EXPECTED_FAMILY_DOMAIN[family][0],
                "larva_up_candidate_genes": family_counts[family],
                "eggNOG_annotation_available": sum(
                    1 for row in detail
                    if row["merops_family"] == family and row["eggNOG_annotation_available"] == "yes"
                ),
                "expected_domain_detected_in_PFAMs": status_counts[
                    "expected_domain_detected_in_PFAMs"
                ],
                "expected_domain_named_outside_PFAMs": status_counts[
                    "expected_domain_named_outside_PFAMs"
                ],
                "expected_domain_not_detected_in_eggNOG": status_counts[
                    "expected_domain_not_detected_in_eggNOG"
                ],
                "eggNOG_annotation_not_available": status_counts[
                    "eggNOG_annotation_not_available"
                ],
            }
        )

    errors: list[str] = []
    if len(detail) != 9:
        errors.append("Expected exactly 9 larva-up I01/I04/I19 candidates")
    if len({str(row["gene_id"]) for row in detail}) != len(detail):
        errors.append("Candidate gene IDs are not unique")
    if family_counts != Counter(EXPECTED_COUNTS):
        errors.append("Family membership differs from frozen counts: I01=2, I04=4, I19=3")
    if any(row["DESeq2_status"] != "larva_up" for row in detail):
        errors.append("Shortlist contains a non-larva-up candidate")
    if any(row["analysis_class"] != "classical_secreted_core" for row in detail):
        errors.append("Shortlist contains a candidate outside classical_secreted_core")

    candidate_fields = list(shortlist[0].keys())
    new_fields = [
        "expected_family_domain", "domain_evidence_status",
        "eggNOG_annotation_available", "eggNOG_query", "eggNOG_seed_ortholog",
        "eggNOG_Description", "eggNOG_Preferred_name", "eggNOG_GOs",
        "eggNOG_KEGG_ko", "eggNOG_CAZy", "eggNOG_PFAMs", "eggNOG_OGs",
        "interpretation_boundary",
    ]
    write_tsv(DETAIL_OUT, detail, candidate_fields + new_fields)
    write_tsv(SUMMARY_OUT, summary, list(summary[0].keys()))

    report = [
        "===== INSECT LARVA-UP INHIBITOR ANNOTATION/DOMAIN AUDIT V2 =====",
        f"Input shortlist: {SHORTLIST_FILE}",
        f"eggNOG annotation input: {EGGNOG_FILE}",
        f"Larva-up candidates audited: {len(detail)}",
        f"eggNOG annotations matched: {sum(row['eggNOG_annotation_available'] == 'yes' for row in detail)}",
        "",
        "Family-level annotation/domain summary:",
    ]
    for row in summary:
        report.append(
            "  {merops_family}: {larva_up_candidate_genes} candidates; "
            "Pfam expected-domain support={expected_domain_detected_in_PFAMs}; "
            "annotation-only support={expected_domain_named_outside_PFAMs}; "
            "not detected={expected_domain_not_detected_in_eggNOG}; "
            "annotation unavailable={eggNOG_annotation_not_available}".format(**row)
        )
    report.extend(
        [
            "",
            f"Errors detected: {len(errors)}",
            *[f"ERROR: {error}" for error in errors],
            (
                "INSECT LARVA-UP INHIBITOR ANNOTATION/DOMAIN AUDIT V2: CHECK REQUIRED"
                if errors
                else "INSECT LARVA-UP INHIBITOR ANNOTATION/DOMAIN AUDIT V2: PASS"
            ),
        ]
    )
    AUDIT_OUT.write_text("\n".join(report) + "\n", encoding="utf-8")

    print("\n".join(report))
    print(f"\nOutput: {DETAIL_OUT}")
    print(f"Output: {SUMMARY_OUT}")
    print(f"Output: {AUDIT_OUT}")

    if errors:
        raise SystemExit("Insect larva-up inhibitor annotation/domain audit failed")


if __name__ == "__main__":
    main()
