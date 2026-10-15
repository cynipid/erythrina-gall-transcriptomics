#!/usr/bin/env python3
"""Attach frozen plant eggNOG annotations to expanded-mask DESeq2 results."""

from __future__ import annotations

import csv
import hashlib
import os
import sys
from pathlib import Path


HOME = Path.home()
ROOT = HOME / "projects/qe_cross_species_discovery"
RESULTS = ROOT / "05_host_focused_DESeq2_sensitivity/results"
LOGS = ROOT / "05_host_focused_DESeq2_sensitivity/logs"

INPUT_RESULTS = RESULTS / "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_all.tsv"
ANNOTATION_DONOR = (
    HOME
    / "projects/qe_plant_denovo/08_results_primary/plant_gene_DESeq2_SG_vs_S_all_primary_eggnog_joined.tsv"
)
OUTPUT_JOINED = RESULTS / "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_eggnog_joined.tsv"
OUTPUT_PARAMETERS = ROOT / "00_project_control/26_expanded_mask_host_focused_annotation_join_parameters.tsv"
OUTPUT_AUDIT = LOGS / "26_expanded_mask_host_focused_annotation_join_audit.txt"

EXPECTED_RESULTS = 24_931
EXPECTED_DONOR = 28_189
DESEQ_COLUMNS = {
    "gene_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "DESeq2_status",
    "DESeq2_status_host_focused",
    "DESeq2_status_host_focused_expanded_mask",
    "eggnog_annotation_match",
}


class AuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditError(f"Missing TSV header: {path}")
        fields = [field.strip().lstrip("\ufeff") for field in reader.fieldnames]
        if len(fields) != len(set(fields)):
            duplicates = sorted({field for field in fields if fields.count(field) > 1})
            raise AuditError(
                f"Duplicate TSV columns in {path}: {', '.join(duplicates)}"
            )
        rows = [
            {
                key.strip().lstrip("\ufeff"): (value or "").strip()
                for key, value in row.items()
            }
            for row in reader
        ]
    return fields, rows


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
    for path in (INPUT_RESULTS, ANNOTATION_DONOR):
        if not path.is_file():
            raise AuditError(f"Missing required input: {path}")
    if OUTPUT_AUDIT.is_file() and OUTPUT_AUDIT.read_text(
        encoding="utf-8", errors="replace"
    ).rstrip().endswith("EXPANDED-MASK HOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN: PASS"):
        raise AuditError(
            f"A completed validated audit already exists: {OUTPUT_AUDIT}. "
            "No outputs were changed."
        )

    hashes = {
        "expanded_mask_DESeq2_results": sha256(INPUT_RESULTS),
        "frozen_eggnog_annotation_donor": sha256(ANNOTATION_DONOR),
    }
    result_fields, result_rows = read_tsv(INPUT_RESULTS)
    donor_fields, donor_rows = read_tsv(ANNOTATION_DONOR)
    if "gene_id" not in result_fields or "gene_id" not in donor_fields:
        raise AuditError("Both inputs must contain gene_id")
    if "DESeq2_status_host_focused_expanded_mask" not in result_fields:
        raise AuditError("Expanded DESeq2 status column is absent")
    if len(result_rows) != EXPECTED_RESULTS:
        raise AuditError(
            f"Expanded result rows: {len(result_rows)}; expected {EXPECTED_RESULTS}"
        )
    if len(donor_rows) != EXPECTED_DONOR:
        raise AuditError(
            f"Frozen annotation-donor rows: {len(donor_rows)}; expected {EXPECTED_DONOR}"
        )

    results = index_unique(result_rows, "gene_id", "expanded DESeq2 results")
    donor = index_unique(donor_rows, "gene_id", "frozen annotation donor")
    missing_donor_rows = sorted(set(results).difference(donor))
    if missing_donor_rows:
        raise AuditError(
            f"Expanded result genes absent from annotation donor: {len(missing_donor_rows)}; "
            + ", ".join(missing_donor_rows[:20])
        )

    annotation_fields = [
        field for field in donor_fields if field not in DESEQ_COLUMNS
    ]
    if not annotation_fields:
        raise AuditError("No annotation columns were identified in the donor")
    if len(annotation_fields) != len(set(annotation_fields)):
        raise AuditError("Annotation field selection contains duplicates")
    output_fields = [
        *result_fields,
        "eggnog_annotation_match",
        *annotation_fields,
    ]
    if len(output_fields) != len(set(output_fields)):
        duplicates = sorted(
            {field for field in output_fields if output_fields.count(field) > 1}
        )
        raise AuditError(
            "Joined output would contain duplicate columns: " + ", ".join(duplicates)
        )

    joined_rows: list[dict[str, object]] = []
    annotation_present = 0
    for result_row in result_rows:
        gene = result_row["gene_id"]
        donor_row = donor[gene]
        joined: dict[str, object] = dict(result_row)
        joined["eggnog_annotation_match"] = "matched"
        for field in annotation_fields:
            joined[field] = donor_row.get(field, "")
        if donor_row.get("has_eggnog_annotation", "").strip().upper() == "TRUE":
            annotation_present += 1
        joined_rows.append(joined)

    write_tsv(OUTPUT_JOINED, output_fields, joined_rows)

    # Re-read the product and verify exact preservation of all DESeq2 fields.
    check_fields, check_rows = read_tsv(OUTPUT_JOINED)
    if check_fields != output_fields or len(check_rows) != EXPECTED_RESULTS:
        raise AuditError("Joined output header or row count failed verification")
    checked = index_unique(check_rows, "gene_id", "joined output")
    for gene, source in results.items():
        output = checked[gene]
        for field in result_fields:
            if output[field] != source[field]:
                raise AuditError(
                    f"DESeq2 field changed during annotation join: {gene}, {field}"
                )

    parameters = [
        {
            "parameter": "join_key",
            "value": "gene_id; one-to-one left join from expanded-mask DESeq2 results",
        },
        {
            "parameter": "annotation_transfer_rule",
            "value": "Only non-DESeq2 donor columns transferred; no DESeq2 value replaced",
        },
        {
            "parameter": "expanded_result_rows",
            "value": len(result_rows),
        },
        {
            "parameter": "annotation_donor_rows",
            "value": len(donor_rows),
        },
        {
            "parameter": "SHA256_expanded_mask_DESeq2_results",
            "value": hashes["expanded_mask_DESeq2_results"],
        },
        {
            "parameter": "SHA256_frozen_eggnog_annotation_donor",
            "value": hashes["frozen_eggnog_annotation_donor"],
        },
    ]
    write_tsv(OUTPUT_PARAMETERS, ["parameter", "value"], parameters)

    stable = (
        sha256(INPUT_RESULTS) == hashes["expanded_mask_DESeq2_results"]
        and sha256(ANNOTATION_DONOR) == hashes["frozen_eggnog_annotation_donor"]
    )
    if not stable:
        raise AuditError("An input changed during the annotation join")

    report = [
        "===== EXPANDED-MASK HOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN =====",
        f"Expanded-mask DESeq2 genes: {len(result_rows)}",
        f"Frozen eggNOG donor rows: {len(donor_rows)}",
        f"Genes with annotation-donor match: {len(joined_rows)}",
        f"Genes without annotation-donor match: {len(missing_donor_rows)}",
        f"Genes marked with eggNOG annotation: {annotation_present}",
        f"Annotation columns transferred: {len(annotation_fields)}",
        "DESeq2 field preservation: PASS",
        "Duplicate output columns: 0",
        "Input SHA-256 verification after output generation: PASS",
        "Errors detected: 0",
        "",
        f"Output: {OUTPUT_JOINED}",
        f"Output: {OUTPUT_PARAMETERS}",
        f"Output: {OUTPUT_AUDIT}",
        "EXPANDED-MASK HOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN: PASS",
    ]
    write_text(OUTPUT_AUDIT, "\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    try:
        main()
    except (AuditError, OSError, ValueError, KeyError) as error:
        print(
            "===== EXPANDED-MASK HOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN =====\n"
            f"ERROR: {error}\n"
            "EXPANDED-MASK HOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN: CHECK REQUIRED",
            file=sys.stderr,
        )
        raise SystemExit(1)
