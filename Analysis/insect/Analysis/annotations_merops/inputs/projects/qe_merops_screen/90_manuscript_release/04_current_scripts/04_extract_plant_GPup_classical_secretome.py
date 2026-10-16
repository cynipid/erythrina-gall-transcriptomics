#!/usr/bin/env python3
"""Extract GP-up classical-secreted E. variegata proteins for MEROPS screening."""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path


if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: 04_extract_plant_GPup_classical_secretome.py "
        "INTEGRATED_PLANT.tsv REPRESENTATIVE_PROTEINS.faa OUTPUT_DIRECTORY"
    )

integrated_path = Path(sys.argv[1])
fasta_path = Path(sys.argv[2])
out_dir = Path(sys.argv[3])
out_dir.mkdir(parents=True, exist_ok=True)

EXPECTED_INTEGRATED_ROWS = 41330
EXPECTED_SELECTED_GENES = 549


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pick_column(fieldnames: list[str], choices: list[str], label: str) -> str:
    for choice in choices:
        if choice in fieldnames:
            return choice
    raise ValueError(f"{label} not found. Expected one of: {', '.join(choices)}")


def fasta_records(path: Path):
    header = None
    sequence: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line
                sequence = []
            elif line:
                sequence.append(line.strip())
    if header is not None:
        yield header, "".join(sequence)


if not integrated_path.is_file():
    raise SystemExit(f"Missing integrated plant table: {integrated_path}")
if not fasta_path.is_file():
    raise SystemExit(f"Missing representative plant FASTA: {fasta_path}")

errors: list[str] = []
with integrated_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if not reader.fieldnames:
        raise SystemExit("Integrated plant table has no header")
    fields = reader.fieldnames
    id_col = pick_column(fields, ["gene_id", "Geneid", "GeneID"], "Gene-ID column")
    status_col = pick_column(
        fields,
        ["DESeq2_status", "deseq2_status", "DE_status", "DESeq2_primary_status"],
        "DESeq2-status column",
    )
    class_col = pick_column(
        fields,
        ["analysis_class", "Analysis_class", "signalp_tmhmm_class"],
        "SignalP/TM-class column",
    )
    integrated_rows = list(reader)

if len(integrated_rows) != EXPECTED_INTEGRATED_ROWS:
    errors.append(
        f"Integrated plant rows={len(integrated_rows)}; expected {EXPECTED_INTEGRATED_ROWS}"
    )
if len({row[id_col] for row in integrated_rows}) != len(integrated_rows):
    errors.append("Integrated plant table contains duplicate gene IDs")

selected_rows = [
    row for row in integrated_rows
    if row[status_col] == "GP_up" and row[class_col] == "classical_secreted_core"
]
selected_ids = {row[id_col] for row in selected_rows}
if len(selected_ids) != EXPECTED_SELECTED_GENES:
    errors.append(
        f"Selected GP-up classical-secreted genes={len(selected_ids)}; expected {EXPECTED_SELECTED_GENES}"
    )

selected_fasta = out_dir / "Evariegata_GP_up_classical_secreted_proteins.faa"
selected_table = out_dir / "Evariegata_GP_up_classical_secreted_proteins.tsv"
audit_file = out_dir / "Evariegata_GP_up_classical_secreted_proteins_audit.txt"

found_ids: set[str] = set()
total_fasta_records = 0
with selected_fasta.open("w") as out_handle:
    for header, sequence in fasta_records(fasta_path):
        total_fasta_records += 1
        gene_id = header[1:].split()[0]
        if gene_id in selected_ids:
            found_ids.add(gene_id)
            out_handle.write(f">{gene_id}\n")
            for start in range(0, len(sequence), 80):
                out_handle.write(sequence[start : start + 80] + "\n")

missing_ids = sorted(selected_ids - found_ids)
if missing_ids:
    errors.append(f"Selected genes absent from representative FASTA={len(missing_ids)}")
if total_fasta_records != EXPECTED_INTEGRATED_ROWS:
    errors.append(
        f"Representative FASTA records={total_fasta_records}; expected {EXPECTED_INTEGRATED_ROWS}"
    )
if len(found_ids) != len(selected_ids):
    errors.append("Selected FASTA ID count does not equal selected integrated ID count")

# Preserve the whole integrated record for the selected proteins, with gene ID
# as the first field for transparent downstream joining.
field_order = [id_col] + [field for field in fields if field != id_col]
with selected_table.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=field_order, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(sorted(selected_rows, key=lambda row: row[id_col]))

with audit_file.open("w") as handle:
    handle.write("===== PLANT GP-UP CLASSICAL-SECRETOME EXTRACTION =====\n")
    handle.write(f"Integrated rows: {len(integrated_rows)}\n")
    handle.write(f"Representative FASTA records: {total_fasta_records}\n")
    handle.write(f"Selected GP-up classical-secreted genes: {len(selected_ids)}\n")
    handle.write(f"Selected FASTA records: {len(found_ids)}\n")
    handle.write(f"Integrated input SHA-256: {sha256(integrated_path)}\n")
    handle.write(f"FASTA input SHA-256: {sha256(fasta_path)}\n")
    handle.write(f"Selected FASTA SHA-256: {sha256(selected_fasta)}\n")
    handle.write(f"Errors detected: {len(errors)}\n")
    for error in errors:
        handle.write(f"ERROR: {error}\n")
    handle.write(
        "PLANT GP-UP CLASSICAL-SECRETOME EXTRACTION: "
        + ("PASS\n" if not errors else "CHECK REQUIRED\n")
    )

print("===== PLANT GP-UP CLASSICAL-SECRETOME EXTRACTION =====")
print(f"Integrated rows: {len(integrated_rows)}")
print(f"Representative FASTA records: {total_fasta_records}")
print(f"Selected GP-up classical-secreted genes: {len(selected_ids)}")
print(f"Selected FASTA records: {len(found_ids)}")
print(f"Errors detected: {len(errors)}")
print(f"Output: {selected_fasta}")
print(f"Output: {selected_table}")
print(f"Output: {audit_file}")
if errors:
    raise SystemExit("PLANT GP-UP CLASSICAL-SECRETOME EXTRACTION: CHECK REQUIRED")
print("PLANT GP-UP CLASSICAL-SECRETOME EXTRACTION: PASS")
