#!/usr/bin/env python3
"""Extract frozen discovery candidates and run bounded cross-system similarity screens.

The screen is descriptive. It does not alter discovery rankings, assign
taxonomic origin, or infer molecular mimicry, horizontal transfer, or
contamination from sequence similarity alone.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


HOME = Path.home()
ROOT = HOME / "projects/qe_cross_species_discovery"

INPUTS = {
    "insect_top15": ROOT / "02_insect_discovery/results/insect_larva_secretome_poorly_annotated_top15.tsv",
    "plant_top50": ROOT / "03_plant_discovery/results/plant_GP_top50.tsv",
    "plant_replicate_audit": ROOT / "03_plant_discovery/results/plant_GP_top50_replicate_audit.tsv",
    "insect_representative_fasta": HOME / "projects/qe_merops_screen/01_input/Qerythrinae_insect_representative_proteins.faa",
    "plant_representative_fasta": HOME / "projects/qe_merops_screen/01_input/Evariegata_plant_representative_proteins.faa",
}

EXPECTED_SHA256 = {
    "insect_top15": "999e95487f565dc3a38f4b9b5514403106d4b8119a4df11eadf35212437e6b8a",
    "plant_top50": "1d4aff7ecab264a97c5c7944a74919293c6f8514aded0f7b56e7cf8a672aa4bb",
    "plant_replicate_audit": "f09f539e87d8e976f9130adb2c956c7862f6af79bacf07fed3028e252527759d",
    "insect_representative_fasta": "fc1b4ac07d4b1c43c97f79f8bd93b39179a65b950cf3ad4e4010b34aef93bdf1",
    "plant_representative_fasta": "cb4683b7a4b1b20dfdfe6720b0bc060a3f46ec808f6bc223fb5f55d7ee8160b4",
}

OUTPUTS = {
    "insect_fasta": ROOT / "02_insect_discovery/results/insect_larva_secretome_poorly_annotated_top15.faa",
    "plant_fasta": ROOT / "03_plant_discovery/results/plant_GP_top50_predicted_proteins.faa",
    "plant_no_protein": ROOT / "03_plant_discovery/results/plant_GP_top50_without_predicted_protein.tsv",
    "plant_flagged_fasta": ROOT / "03_plant_discovery/results/plant_GP_replicate_flagged6_proteins.faa",
    "candidate_inventory": ROOT / "04_cross_species_integration/results/discovery_sequence_provenance_candidate_inventory.tsv",
    "plant_vs_insect_all": ROOT / "04_cross_species_integration/results/plant_GP_top50_proteins_vs_insect_representative_all_hits.tsv",
    "plant_vs_insect_best": ROOT / "04_cross_species_integration/results/plant_GP_top50_proteins_vs_insect_representative_best_hits.tsv",
    "insect_vs_plant_all": ROOT / "04_cross_species_integration/results/insect_larva_top15_vs_plant_representative_all_hits.tsv",
    "insect_vs_plant_best": ROOT / "04_cross_species_integration/results/insect_larva_top15_vs_plant_representative_best_hits.tsv",
    "plant_internal_all": ROOT / "04_cross_species_integration/results/plant_GP_top50_candidate_internal_nonself_hits.tsv",
    "plant_internal_summary": ROOT / "04_cross_species_integration/results/plant_GP_top50_candidate_internal_similarity_summary.tsv",
    "parameters": ROOT / "00_project_control/04_discovery_sequence_provenance_parameters.tsv",
    "audit": ROOT / "08_logs/04_discovery_sequence_provenance_screen_audit.txt",
}

DIAMOND = ["micromamba", "run", "-n", "merops", "diamond"]
EVALUE_LIMIT = 1e-5
MAX_TARGET_SEQS = 25
THREADS = 8

ALIGNMENT_FIELDS = [
    "query_id",
    "subject_id",
    "percent_identity",
    "alignment_length",
    "query_length",
    "subject_length",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bitscore",
]

DERIVED_HIT_FIELDS = [
    "query_span_coverage_percent",
    "subject_span_coverage_percent",
    "cross_system_similarity_class",
    "near_identity_cross_system_review",
]

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AMINO_ACIDS = set("BXZJUO")


class AuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv_atomic(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
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


def read_fasta(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    identifier: str | None = None
    description = ""
    sequence_parts: list[str] = []

    def store() -> None:
        if identifier is None:
            return
        if identifier in records:
            raise AuditError(f"Duplicate FASTA identifier in {path}: {identifier}")
        sequence = "".join(sequence_parts).upper().rstrip("*")
        records[identifier] = {
            "description": description,
            "sequence": sequence,
        }

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                full_header = line[1:]
                pieces = full_header.split(maxsplit=1)
                identifier = pieces[0]
                description = pieces[1] if len(pieces) == 2 else ""
                sequence_parts = []
            else:
                if identifier is None:
                    raise AuditError(f"Sequence before first FASTA header: {path}")
                sequence_parts.append(line)
    store()
    return records


def write_fasta_atomic(
    path: Path,
    identifiers: list[str],
    records: dict[str, dict[str, str]],
    descriptions: dict[str, str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier in identifiers:
            sequence = records[identifier]["sequence"]
            description = descriptions.get(identifier, "")
            header = f">{identifier}"
            if description:
                header += f" {description}"
            handle.write(header + "\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")
    os.replace(temporary, path)


def verify_input_hashes() -> dict[str, str]:
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


def require_unique(rows: list[dict[str, str]], key: str, label: str) -> None:
    values = [row.get(key, "") for row in rows]
    if any(not value for value in values):
        raise AuditError(f"Blank {key} in {label}")
    if len(values) != len(set(values)):
        raise AuditError(f"Duplicate {key} in {label}")


def parse_header_value(description: str, key: str) -> str | None:
    match = re.search(rf"(?:^|\s){re.escape(key)}=([^\s]+)", description)
    return match.group(1) if match else None


def validate_sequence(sequence: str, identifier: str) -> tuple[float, str]:
    if not sequence:
        raise AuditError(f"Empty protein sequence: {identifier}")
    invalid = sorted(set(sequence) - STANDARD_AMINO_ACIDS - AMBIGUOUS_AMINO_ACIDS)
    if invalid:
        raise AuditError(f"Invalid amino-acid characters in {identifier}: {''.join(invalid)}")
    ambiguous_fraction = sum(char in AMBIGUOUS_AMINO_ACIDS for char in sequence) / len(sequence)
    return ambiguous_fraction, "PASS" if ambiguous_fraction <= 0.05 else "REVIEW"


def number(value: str, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise AuditError(f"Non-numeric value for {label}: {value!r}") from exc
    if not math.isfinite(result):
        raise AuditError(f"Non-finite value for {label}: {value!r}")
    return result


def run_command(command: list[str], stdout_path: Path | None = None) -> str:
    if stdout_path is None:
        completed = subprocess.run(
            command,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return completed.stdout.strip()

    with stdout_path.open("w", encoding="utf-8", newline="\n") as output_handle:
        subprocess.run(
            command,
            check=True,
            text=True,
            stdout=output_handle,
            stderr=subprocess.PIPE,
        )
    return ""


def diamond_makedb(input_fasta: Path, database_prefix: Path) -> None:
    run_command(
        [
            *DIAMOND,
            "makedb",
            "--in",
            str(input_fasta),
            "--db",
            str(database_prefix),
            "--threads",
            str(THREADS),
            "--quiet",
        ]
    )


def diamond_blastp(query: Path, database_prefix: Path, output: Path, max_targets: int) -> None:
    run_command(
        [
            *DIAMOND,
            "blastp",
            "--query",
            str(query),
            "--db",
            str(database_prefix),
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qlen",
            "slen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "--evalue",
            str(EVALUE_LIMIT),
            "--max-target-seqs",
            str(max_targets),
            "--sensitive",
            "--threads",
            str(THREADS),
            "--quiet",
        ],
        stdout_path=output,
    )


def similarity_class(
    percent_identity: float,
    query_coverage: float,
    subject_coverage: float,
    evalue: float,
) -> tuple[str, bool]:
    near_identity = (
        percent_identity >= 90.0
        and query_coverage >= 80.0
        and subject_coverage >= 80.0
    )
    if near_identity:
        return "near_identity_cross_system_review", True
    if (
        evalue <= 1e-10
        and percent_identity >= 40.0
        and query_coverage >= 60.0
        and subject_coverage >= 60.0
    ):
        return "substantial_cross_system_homology", False
    if (
        evalue <= EVALUE_LIMIT
        and percent_identity >= 25.0
        and query_coverage >= 40.0
    ):
        return "partial_cross_system_homology", False
    return "weak_or_short_cross_system_match", False


def internal_similarity_label(cross_system_label: str) -> str:
    """Translate shared numerical tiers into within-plant wording."""
    translations = {
        "near_identity_cross_system_review": "near_identity_candidate_pair_review",
        "substantial_cross_system_homology": "substantial_candidate_sequence_homology",
        "partial_cross_system_homology": "partial_candidate_sequence_homology",
        "weak_or_short_cross_system_match": "weak_or_short_candidate_match",
    }
    try:
        return translations[cross_system_label]
    except KeyError as exc:
        raise AuditError(
            f"Unexpected similarity class for internal comparison: {cross_system_label}"
        ) from exc


def parse_diamond_raw(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, values in enumerate(reader, start=1):
            if len(values) != len(ALIGNMENT_FIELDS):
                raise AuditError(
                    f"Unexpected DIAMOND field count at {path}:{line_number}: {len(values)}"
                )
            row = dict(zip(ALIGNMENT_FIELDS, values))
            qlen = int(row["query_length"])
            slen = int(row["subject_length"])
            qstart = int(row["query_start"])
            qend = int(row["query_end"])
            sstart = int(row["subject_start"])
            send = int(row["subject_end"])
            qcov = 100.0 * (abs(qend - qstart) + 1) / qlen
            scov = 100.0 * (abs(send - sstart) + 1) / slen
            classification, near_identity = similarity_class(
                float(row["percent_identity"]), qcov, scov, float(row["evalue"])
            )
            row["query_span_coverage_percent"] = format(qcov, ".10g")
            row["subject_span_coverage_percent"] = format(scov, ".10g")
            row["cross_system_similarity_class"] = classification
            row["near_identity_cross_system_review"] = "TRUE" if near_identity else "FALSE"
            rows.append(row)
    return rows


def best_hits_with_no_hit_rows(
    all_hits: list[dict[str, object]],
    query_ids: list[str],
) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in all_hits:
        grouped[str(row["query_id"])].append(row)

    output: list[dict[str, object]] = []
    for query_id in query_ids:
        hits = grouped.get(query_id, [])
        if not hits:
            blank = {field: "" for field in ALIGNMENT_FIELDS + DERIVED_HIT_FIELDS}
            blank["query_id"] = query_id
            blank["cross_system_similarity_class"] = "no_match_at_evalue_1e-5"
            blank["near_identity_cross_system_review"] = "FALSE"
            output.append(blank)
            continue
        best = sorted(
            hits,
            key=lambda row: (-float(row["bitscore"]), float(row["evalue"])),
        )[0]
        output.append(best)
    return output


def candidate_inventory_row(
    organism: str,
    candidate_set: str,
    source_row: dict[str, str],
    fasta_id: str | None,
    protein_id: str | None,
    sequence: str | None,
    replicate_flag: bool,
) -> dict[str, object]:
    if sequence is None:
        return {
            "organism": organism,
            "candidate_set": candidate_set,
            "discovery_rank": source_row["discovery_rank"],
            "gene_id": source_row["gene_id"],
            "fasta_lookup_id": "",
            "representative_protein_id": protein_id or "",
            "representative_protein_available": "FALSE",
            "sequence_length_aa": "",
            "sequence_sha256": "",
            "ambiguous_amino_acid_fraction": "",
            "sequence_quality_status": "not_applicable_no_predicted_protein",
            "replicate_review_flag": "TRUE" if replicate_flag else "FALSE",
            "annotation_status": source_row.get("annotation_status", ""),
            "discovery_classification": source_row.get("discovery_classification", ""),
        }
    ambiguous_fraction, quality_status = validate_sequence(sequence, fasta_id or "")
    return {
        "organism": organism,
        "candidate_set": candidate_set,
        "discovery_rank": source_row["discovery_rank"],
        "gene_id": source_row["gene_id"],
        "fasta_lookup_id": fasta_id or "",
        "representative_protein_id": protein_id or "",
        "representative_protein_available": "TRUE",
        "sequence_length_aa": str(len(sequence)),
        "sequence_sha256": sequence_sha256(sequence),
        "ambiguous_amino_acid_fraction": format(ambiguous_fraction, ".10g"),
        "sequence_quality_status": quality_status,
        "replicate_review_flag": "TRUE" if replicate_flag else "FALSE",
        "annotation_status": source_row.get("annotation_status", ""),
        "discovery_classification": source_row.get("discovery_classification", ""),
    }


def main() -> None:
    observed_hashes = verify_input_hashes()
    if shutil.which("micromamba") is None:
        raise AuditError("micromamba is not available on PATH")
    diamond_version = run_command([*DIAMOND, "version"])

    insect_fields, insect_rows = read_tsv(INPUTS["insect_top15"])
    plant_fields, plant_rows = read_tsv(INPUTS["plant_top50"])
    audit_fields, audit_rows = read_tsv(INPUTS["plant_replicate_audit"])
    require_unique(insect_rows, "gene_id", "insect top 15")
    require_unique(plant_rows, "gene_id", "plant top 50")
    require_unique(audit_rows, "gene_id", "plant replicate audit")
    if len(insect_rows) != 15:
        raise AuditError(f"Expected 15 insect candidates; found {len(insect_rows)}")
    if len(plant_rows) != 50:
        raise AuditError(f"Expected 50 plant candidates; found {len(plant_rows)}")

    insect_proteome = read_fasta(INPUTS["insect_representative_fasta"])
    plant_proteome = read_fasta(INPUTS["plant_representative_fasta"])
    if len(insect_proteome) != 13455:
        raise AuditError(f"Expected 13,455 insect proteins; found {len(insect_proteome)}")
    if len(plant_proteome) != 41330:
        raise AuditError(f"Expected 41,330 plant proteins; found {len(plant_proteome)}")

    flagged_ids = {
        row["gene_id"]
        for row in audit_rows
        if row.get("replicate_review_status") == "review"
    }
    if len(flagged_ids) != 6:
        raise AuditError(f"Expected six plant replicate-review genes; found {len(flagged_ids)}")

    insect_ids = [row["gene_id"] for row in insect_rows]
    missing_insect = [identifier for identifier in insect_ids if identifier not in insect_proteome]
    if missing_insect:
        raise AuditError(f"Insect proteins missing from FASTA: {', '.join(missing_insect)}")

    plant_protein_rows = [
        row for row in plant_rows if row.get("representative_protein_available") == "TRUE"
    ]
    plant_no_protein_rows = [
        row for row in plant_rows if row.get("representative_protein_available") != "TRUE"
    ]
    if len(plant_protein_rows) != 44 or len(plant_no_protein_rows) != 6:
        raise AuditError(
            f"Unexpected plant protein availability: {len(plant_protein_rows)} with protein, "
            f"{len(plant_no_protein_rows)} without protein"
        )
    plant_ids = [row["gene_id"] for row in plant_protein_rows]
    missing_plant = [identifier for identifier in plant_ids if identifier not in plant_proteome]
    if missing_plant:
        raise AuditError(f"Plant gene IDs missing from FASTA: {', '.join(missing_plant)}")

    inventory_rows: list[dict[str, object]] = []
    insect_descriptions: dict[str, str] = {}
    for row in insect_rows:
        gene_id = row["gene_id"]
        sequence = insect_proteome[gene_id]["sequence"]
        expected_length = int(number(row["representative_protein_length_aa"], gene_id))
        if len(sequence) != expected_length:
            raise AuditError(
                f"Insect length mismatch for {gene_id}: table={expected_length}, FASTA={len(sequence)}"
            )
        insect_descriptions[gene_id] = (
            f"discovery_rank={row['discovery_rank']} annotation_status={row['annotation_status']}"
        )
        inventory_rows.append(
            candidate_inventory_row(
                "Q. erythrinae",
                "larva_up_classical_secretome_poorly_annotated_top15",
                row,
                gene_id,
                gene_id,
                sequence,
                False,
            )
        )

    plant_descriptions: dict[str, str] = {}
    for row in plant_rows:
        gene_id = row["gene_id"]
        protein_id = row.get("representative_protein_id", "") or None
        has_protein = row.get("representative_protein_available") == "TRUE"
        if not has_protein:
            inventory_rows.append(
                candidate_inventory_row(
                    "E. variegata",
                    "GP_up_top50",
                    row,
                    None,
                    protein_id,
                    None,
                    gene_id in flagged_ids,
                )
            )
            continue

        fasta_record = plant_proteome[gene_id]
        sequence = fasta_record["sequence"]
        header_protein_id = parse_header_value(
            fasta_record["description"], "representative_protein"
        )
        header_length = parse_header_value(fasta_record["description"], "length_aa")
        expected_length = int(number(row["representative_protein_length_aa"], gene_id))
        if header_protein_id != protein_id:
            raise AuditError(
                f"Plant header protein-ID mismatch for {gene_id}: table={protein_id}, "
                f"header={header_protein_id}"
            )
        if header_length is None or int(header_length) != expected_length:
            raise AuditError(
                f"Plant header length mismatch for {gene_id}: table={expected_length}, "
                f"header={header_length}"
            )
        if len(sequence) != expected_length:
            raise AuditError(
                f"Plant sequence length mismatch for {gene_id}: table={expected_length}, "
                f"FASTA={len(sequence)}"
            )
        plant_descriptions[gene_id] = (
            f"representative_protein={protein_id} discovery_rank={row['discovery_rank']} "
            f"annotation_status={row['annotation_status']}"
        )
        inventory_rows.append(
            candidate_inventory_row(
                "E. variegata",
                "GP_up_top50",
                row,
                gene_id,
                protein_id,
                sequence,
                gene_id in flagged_ids,
            )
        )

    write_fasta_atomic(
        OUTPUTS["insect_fasta"], insect_ids, insect_proteome, insect_descriptions
    )
    write_fasta_atomic(
        OUTPUTS["plant_fasta"], plant_ids, plant_proteome, plant_descriptions
    )
    flagged_protein_ids = [
        row["gene_id"] for row in plant_protein_rows if row["gene_id"] in flagged_ids
    ]
    if len(flagged_protein_ids) != 6:
        raise AuditError(
            f"Expected all six flagged plant genes to have proteins; found {len(flagged_protein_ids)}"
        )
    write_fasta_atomic(
        OUTPUTS["plant_flagged_fasta"],
        flagged_protein_ids,
        plant_proteome,
        plant_descriptions,
    )

    no_protein_output_rows: list[dict[str, object]] = []
    for row in plant_no_protein_rows:
        out = dict(row)
        out["sequence_provenance_status"] = "no_representative_protein_predicted"
        out["required_followup"] = "transcript_level_sequence_assessment"
        no_protein_output_rows.append(out)
    write_tsv_atomic(
        OUTPUTS["plant_no_protein"],
        [*plant_fields, "sequence_provenance_status", "required_followup"],
        no_protein_output_rows,
    )

    inventory_fields = [
        "organism",
        "candidate_set",
        "discovery_rank",
        "gene_id",
        "fasta_lookup_id",
        "representative_protein_id",
        "representative_protein_available",
        "sequence_length_aa",
        "sequence_sha256",
        "ambiguous_amino_acid_fraction",
        "sequence_quality_status",
        "replicate_review_flag",
        "annotation_status",
        "discovery_classification",
    ]
    write_tsv_atomic(OUTPUTS["candidate_inventory"], inventory_fields, inventory_rows)

    OUTPUTS["audit"].parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="diamond_sequence_provenance_") as temp_name:
        temp_dir = Path(temp_name)
        insect_db = temp_dir / "insect_representative"
        plant_db = temp_dir / "plant_representative"
        plant_candidate_db = temp_dir / "plant_candidate"
        diamond_makedb(INPUTS["insect_representative_fasta"], insect_db)
        diamond_makedb(INPUTS["plant_representative_fasta"], plant_db)
        diamond_makedb(OUTPUTS["plant_fasta"], plant_candidate_db)

        plant_vs_insect_raw = temp_dir / "plant_vs_insect.raw.tsv"
        insect_vs_plant_raw = temp_dir / "insect_vs_plant.raw.tsv"
        plant_internal_raw = temp_dir / "plant_internal.raw.tsv"
        diamond_blastp(
            OUTPUTS["plant_fasta"], insect_db, plant_vs_insect_raw, MAX_TARGET_SEQS
        )
        diamond_blastp(
            OUTPUTS["insect_fasta"], plant_db, insect_vs_plant_raw, MAX_TARGET_SEQS
        )
        diamond_blastp(
            OUTPUTS["plant_fasta"], plant_candidate_db, plant_internal_raw, 100
        )

        plant_vs_insect_hits = parse_diamond_raw(plant_vs_insect_raw)
        insect_vs_plant_hits = parse_diamond_raw(insect_vs_plant_raw)
        plant_internal_hits = [
            row
            for row in parse_diamond_raw(plant_internal_raw)
            if row["query_id"] != row["subject_id"]
        ]
        for row in plant_internal_hits:
            row["internal_similarity_class"] = internal_similarity_label(
                str(row.pop("cross_system_similarity_class"))
            )
            row["near_identity_candidate_pair_review"] = row.pop(
                "near_identity_cross_system_review"
            )

    hit_fields = [*ALIGNMENT_FIELDS, *DERIVED_HIT_FIELDS]
    write_tsv_atomic(
        OUTPUTS["plant_vs_insect_all"], hit_fields, plant_vs_insect_hits
    )
    write_tsv_atomic(
        OUTPUTS["insect_vs_plant_all"], hit_fields, insect_vs_plant_hits
    )
    plant_best = best_hits_with_no_hit_rows(plant_vs_insect_hits, plant_ids)
    insect_best = best_hits_with_no_hit_rows(insect_vs_plant_hits, insect_ids)
    write_tsv_atomic(OUTPUTS["plant_vs_insect_best"], hit_fields, plant_best)
    write_tsv_atomic(OUTPUTS["insect_vs_plant_best"], hit_fields, insect_best)
    internal_hit_fields = [
        *ALIGNMENT_FIELDS,
        "query_span_coverage_percent",
        "subject_span_coverage_percent",
        "internal_similarity_class",
        "near_identity_candidate_pair_review",
    ]
    write_tsv_atomic(
        OUTPUTS["plant_internal_all"], internal_hit_fields, plant_internal_hits
    )

    internal_by_query: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in plant_internal_hits:
        internal_by_query[str(row["query_id"])].append(row)
    internal_summary: list[dict[str, object]] = []
    for row in plant_protein_rows:
        query_id = row["gene_id"]
        hits = internal_by_query.get(query_id, [])
        if not hits:
            internal_summary.append(
                {
                    "discovery_rank": row["discovery_rank"],
                    "gene_id": query_id,
                    "best_nonself_candidate_match": "",
                    "percent_identity": "",
                    "query_span_coverage_percent": "",
                    "subject_span_coverage_percent": "",
                    "evalue": "",
                    "bitscore": "",
                    "internal_similarity_class": "no_nonself_match_at_evalue_1e-5",
                }
            )
            continue
        best = sorted(
            hits,
            key=lambda hit: (-float(hit["bitscore"]), float(hit["evalue"])),
        )[0]
        internal_summary.append(
            {
                "discovery_rank": row["discovery_rank"],
                "gene_id": query_id,
                "best_nonself_candidate_match": best["subject_id"],
                "percent_identity": best["percent_identity"],
                "query_span_coverage_percent": best["query_span_coverage_percent"],
                "subject_span_coverage_percent": best["subject_span_coverage_percent"],
                "evalue": best["evalue"],
                "bitscore": best["bitscore"],
                "internal_similarity_class": best["internal_similarity_class"],
            }
        )
    write_tsv_atomic(
        OUTPUTS["plant_internal_summary"],
        [
            "discovery_rank",
            "gene_id",
            "best_nonself_candidate_match",
            "percent_identity",
            "query_span_coverage_percent",
            "subject_span_coverage_percent",
            "evalue",
            "bitscore",
            "internal_similarity_class",
        ],
        internal_summary,
    )

    parameter_rows = [
        {"parameter": "DIAMOND_version", "value": diamond_version},
        {"parameter": "DIAMOND_mode", "value": "blastp --sensitive"},
        {"parameter": "evalue_limit", "value": str(EVALUE_LIMIT)},
        {"parameter": "maximum_target_sequences", "value": str(MAX_TARGET_SEQS)},
        {"parameter": "threads", "value": str(THREADS)},
        {
            "parameter": "near_identity_cross_system_review",
            "value": "identity>=90%; query span coverage>=80%; subject span coverage>=80%",
        },
        {
            "parameter": "substantial_cross_system_homology",
            "value": "evalue<=1e-10; identity>=40%; query and subject span coverage>=60%",
        },
        {
            "parameter": "partial_cross_system_homology",
            "value": "evalue<=1e-5; identity>=25%; query span coverage>=40%",
        },
        {
            "parameter": "interpretive_boundary",
            "value": "Similarity classes are descriptive review tiers, not taxonomic assignments",
        },
    ]
    write_tsv_atomic(OUTPUTS["parameters"], ["parameter", "value"], parameter_rows)

    # Confirm frozen inputs remain byte-for-byte unchanged after all output generation.
    for name, path in INPUTS.items():
        if sha256(path) != observed_hashes[name]:
            raise AuditError(f"Frozen input changed during analysis: {path}")

    plant_best_classes = Counter(
        str(row["cross_system_similarity_class"]) for row in plant_best
    )
    insect_best_classes = Counter(
        str(row["cross_system_similarity_class"]) for row in insect_best
    )
    internal_pairs = {
        tuple(sorted((str(row["query_id"]), str(row["subject_id"]))))
        for row in plant_internal_hits
    }
    same_component_pair = (
        "TRINITY_DN276_c0_g1",
        "TRINITY_DN276_c0_g3",
    )
    dn276_detected = tuple(sorted(same_component_pair)) in internal_pairs

    report = [
        "===== DISCOVERY SEQUENCE-PROVENANCE SCREEN =====",
        f"Discovery root: {ROOT}",
        "",
        f"Insect poorly annotated proteins extracted: {len(insect_ids)}",
        f"Plant top-50 predicted proteins extracted: {len(plant_ids)}",
        f"Plant top-50 genes without predicted protein: {len(plant_no_protein_rows)}",
        f"Plant replicate-flagged proteins extracted: {len(flagged_protein_ids)}",
        f"Candidate inventory rows: {len(inventory_rows)}",
        "",
        "Plant candidate proteins versus insect representative proteome:",
    ]
    for key in sorted(plant_best_classes):
        report.append(f"  {key}: {plant_best_classes[key]}")
    report.extend(
        [
            "",
            "Insect candidate proteins versus plant representative proteome:",
        ]
    )
    for key in sorted(insect_best_classes):
        report.append(f"  {key}: {insect_best_classes[key]}")
    report.extend(
        [
            "",
            f"Plant candidate internal nonself hit rows: {len(plant_internal_hits)}",
            f"Plant candidate internal unique pairs: {len(internal_pairs)}",
            f"TRINITY_DN276_c0_g1/g3 similarity detected: {'YES' if dn276_detected else 'NO'}",
            "",
            "Interpretive boundary:",
            "  Cross-system similarity is not a taxonomic assignment.",
            "  It does not demonstrate contamination, molecular mimicry, horizontal transfer, or functional equivalence.",
            "  Near-identity matches are review flags only and require independent taxonomic evidence.",
            "  Six plant genes without predicted proteins remain eligible discovery candidates and require transcript-level analysis.",
            "",
            "Frozen-input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "",
        ]
    )
    for key in OUTPUTS:
        report.append(f"Output: {OUTPUTS[key]}")
    report.append("DISCOVERY SEQUENCE-PROVENANCE SCREEN: PASS")

    OUTPUTS["audit"].parent.mkdir(parents=True, exist_ok=True)
    audit_tmp = OUTPUTS["audit"].with_name(OUTPUTS["audit"].name + ".tmp")
    audit_tmp.write_text("\n".join(report) + "\n", encoding="utf-8")
    os.replace(audit_tmp, OUTPUTS["audit"])
    print("\n".join(report))


if __name__ == "__main__":
    try:
        main()
    except (AuditError, subprocess.CalledProcessError) as error:
        print("===== DISCOVERY SEQUENCE-PROVENANCE SCREEN =====", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError):
            print(f"ERROR: command failed: {' '.join(error.cmd)}", file=sys.stderr)
            if error.stderr:
                print(error.stderr.strip(), file=sys.stderr)
        else:
            print(f"ERROR: {error}", file=sys.stderr)
        print("DISCOVERY SEQUENCE-PROVENANCE SCREEN: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
