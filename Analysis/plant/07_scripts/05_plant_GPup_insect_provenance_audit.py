#!/usr/bin/env python3
"""Audit insect-like sequence evidence among protein-bearing plant GP-up DEGs.

This is a taxonomic-provenance quality-control screen. It does not rerun
DESeq2, alter the frozen discovery ranking, remove any gene, or establish the
taxonomic origin of a transcript. A strong match to the Q. erythrinae
representative proteome is a review flag that requires contextual evidence.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


HOME = Path.home()
ROOT = HOME / "projects/qe_cross_species_discovery"
LINKS = ROOT / "01_input_links"

INPUTS = {
    "plant_integration": LINKS / "plant_signalp_tmhmm_deseq2_integrated.tsv",
    "plant_protein_map": LINKS / "plant_representative_protein_map.tsv",
    "plant_top50": ROOT / "03_plant_discovery/results/plant_GP_top50.tsv",
    "prior_top50_best_hits": (
        ROOT
        / "04_cross_species_integration/results/plant_GP_top50_proteins_vs_insect_representative_best_hits.tsv"
    ),
    "plant_representative_fasta": (
        HOME
        / "projects/qe_merops_screen/01_input/Evariegata_plant_representative_proteins.faa"
    ),
    "insect_representative_fasta": (
        HOME
        / "projects/qe_merops_screen/01_input/Qerythrinae_insect_representative_proteins.faa"
    ),
}

EXPECTED_FROZEN_SHA256 = {
    "plant_top50": "1d4aff7ecab264a97c5c7944a74919293c6f8514aded0f7b56e7cf8a672aa4bb",
    "plant_representative_fasta": "cb4683b7a4b1b20dfdfe6720b0bc060a3f46ec808f6bc223fb5f55d7ee8160b4",
    "insect_representative_fasta": "fc1b4ac07d4b1c43c97f79f8bd93b39179a65b950cf3ad4e4010b34aef93bdf1",
}

OUTPUTS = {
    "gpup_query_fasta": (
        ROOT / "03_plant_discovery/results/plant_GPup_predicted_proteins_provenance_audit.faa"
    ),
    "gpup_best": (
        ROOT
        / "04_cross_species_integration/results/plant_GPup_proteins_vs_Qerythrinae_best_hits_provenance.tsv"
    ),
    "near_identity": (
        ROOT
        / "04_cross_species_integration/results/plant_GPup_insect_like_near_identity_review.tsv"
    ),
    "top50_annotated": (
        ROOT / "03_plant_discovery/results/plant_GP_top50_provenance_annotated.tsv"
    ),
    "top50_summary": (
        ROOT / "04_cross_species_integration/results/plant_GP_top50_provenance_summary.tsv"
    ),
    "gpup_summary": (
        ROOT / "04_cross_species_integration/results/plant_GPup_insect_provenance_summary.tsv"
    ),
    "parameters": ROOT / "00_project_control/05_plant_GPup_insect_provenance_parameters.tsv",
    "audit": ROOT / "08_logs/05_plant_GPup_insect_provenance_audit.txt",
}

DIAMOND = ["micromamba", "run", "-n", "merops", "diamond"]
EVALUE_LIMIT = 1e-5
MAX_TARGET_SEQS = 25
THREADS = 8

HIT_FIELDS = [
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
DERIVED_FIELDS = [
    "query_span_coverage_percent",
    "subject_span_coverage_percent",
    "Qerythrinae_similarity_review_class",
    "insect_like_near_identity_review",
]
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = set("BXZJUO")


class AuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean(value: object) -> str:
    return "" if value is None else str(value).strip().lstrip("\ufeff")


def as_float(value: object, label: str) -> float:
    try:
        number = float(clean(value))
    except ValueError as exc:
        raise AuditError(f"Non-numeric {label}: {value!r}") from exc
    if not math.isfinite(number):
        raise AuditError(f"Non-finite {label}: {value!r}")
    return number


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditError(f"Missing TSV header: {path}")
        fields = [clean(field) for field in reader.fieldnames]
        rows = [{clean(key): clean(value) for key, value in row.items()} for row in reader]
    return fields, rows


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


def index_unique(rows: list[dict[str, str]], key: str, label: str) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        value = row.get(key, "")
        if not value:
            raise AuditError(f"Blank {key} in {label}")
        if value in indexed:
            raise AuditError(f"Duplicate {key} in {label}: {value}")
        indexed[value] = row
    return indexed


def parse_header_value(description: str, key: str) -> str | None:
    match = re.search(rf"(?:^|\s){re.escape(key)}=([^\s]+)", description)
    return match.group(1) if match else None


def read_fasta(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    identifier: str | None = None
    description = ""
    parts: list[str] = []

    def store() -> None:
        if identifier is None:
            return
        if identifier in records:
            raise AuditError(f"Duplicate FASTA identifier in {path}: {identifier}")
        sequence = "".join(parts).upper().rstrip("*")
        if not sequence:
            raise AuditError(f"Empty FASTA record in {path}: {identifier}")
        records[identifier] = {"description": description, "sequence": sequence}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                split = line[1:].split(maxsplit=1)
                identifier = split[0]
                description = split[1] if len(split) > 1 else ""
                parts = []
            else:
                if identifier is None:
                    raise AuditError(f"Sequence before a FASTA header in {path}")
                parts.append(line)
    store()
    return records


def validate_protein(sequence: str, gene_id: str) -> None:
    invalid = sorted(set(sequence) - STANDARD_AA - AMBIGUOUS_AA)
    if invalid:
        raise AuditError(f"Invalid amino-acid characters in {gene_id}: {''.join(invalid)}")
    if sum(letter in AMBIGUOUS_AA for letter in sequence) / len(sequence) > 0.05:
        raise AuditError(f"More than 5% ambiguous amino acids in {gene_id}")


def write_fasta(path: Path, rows: list[dict[str, str]], records: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            gene_id = row["gene_id"]
            sequence = records[gene_id]["sequence"]
            handle.write(
                f">{gene_id} representative_protein={row['protein_id']} "
                f"DESeq2_status={row['DESeq2_status']}\n"
            )
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")
    os.replace(temporary, path)


def run(command: list[str], stdout_path: Path | None = None) -> str:
    if stdout_path is None:
        completed = subprocess.run(
            command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        return completed.stdout.strip()
    with stdout_path.open("w", encoding="utf-8", newline="\n") as handle:
        subprocess.run(command, check=True, text=True, stdout=handle, stderr=subprocess.PIPE)
    return ""


def diamond_makedb(input_fasta: Path, prefix: Path) -> None:
    run([
        *DIAMOND, "makedb", "--in", str(input_fasta), "--db", str(prefix),
        "--threads", str(THREADS), "--quiet",
    ])


def diamond_blastp(query: Path, database: Path, output: Path) -> None:
    run([
        *DIAMOND, "blastp", "--query", str(query), "--db", str(database),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "--evalue", str(EVALUE_LIMIT), "--max-target-seqs", str(MAX_TARGET_SEQS),
        "--sensitive", "--threads", str(THREADS), "--quiet",
    ], stdout_path=output)


def classify_hit(identity: float, qcov: float, scov: float, evalue: float) -> tuple[str, str]:
    if identity >= 90.0 and qcov >= 80.0 and scov >= 80.0:
        return "insect_like_near_identity_review", "TRUE"
    if identity >= 40.0 and qcov >= 60.0 and scov >= 60.0 and evalue <= 1e-10:
        return "insect_like_substantial_homology_review", "FALSE"
    if identity >= 25.0 and qcov >= 40.0 and evalue <= EVALUE_LIMIT:
        return "partial_homology_to_Qerythrinae_reference", "FALSE"
    return "weak_or_short_match_to_Qerythrinae_reference", "FALSE"


def parse_hits(path: Path) -> list[dict[str, str]]:
    parsed: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, values in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if len(values) != len(HIT_FIELDS):
                raise AuditError(f"Unexpected DIAMOND field count at {path}:{line_number}")
            hit = dict(zip(HIT_FIELDS, values))
            qlen, slen = int(hit["query_length"]), int(hit["subject_length"])
            qspan = abs(int(hit["query_end"]) - int(hit["query_start"])) + 1
            sspan = abs(int(hit["subject_end"]) - int(hit["subject_start"])) + 1
            qcov, scov = 100 * qspan / qlen, 100 * sspan / slen
            status, near = classify_hit(
                float(hit["percent_identity"]), qcov, scov, float(hit["evalue"])
            )
            hit.update({
                "query_span_coverage_percent": format(qcov, ".10g"),
                "subject_span_coverage_percent": format(scov, ".10g"),
                "Qerythrinae_similarity_review_class": status,
                "insect_like_near_identity_review": near,
            })
            parsed.append(hit)
    return parsed


def choose_best(hits: list[dict[str, str]], query_ids: list[str]) -> dict[str, dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for hit in hits:
        grouped.setdefault(hit["query_id"], []).append(hit)
    result: dict[str, dict[str, str]] = {}
    for query in query_ids:
        candidates = grouped.get(query, [])
        if not candidates:
            blank = {field: "" for field in HIT_FIELDS + DERIVED_FIELDS}
            blank.update({
                "query_id": query,
                "Qerythrinae_similarity_review_class": "no_match_to_Qerythrinae_reference_at_evalue_1e-5",
                "insect_like_near_identity_review": "FALSE",
            })
            result[query] = blank
        else:
            result[query] = sorted(
                candidates, key=lambda row: (-float(row["bitscore"]), float(row["evalue"]))
            )[0]
    return result


def add_context(hit: dict[str, str], row: dict[str, str]) -> dict[str, object]:
    output = {
        "gene_id": row["gene_id"],
        "representative_protein_id": row["protein_id"],
        "representative_protein_length_aa": row["aa_length"],
        "baseMean": row["baseMean"],
        "log2FoldChange_GP_vs_UP": row["log2FoldChange"],
        "padj": row["padj"],
        "DESeq2_status": row["DESeq2_status"],
        "analysis_class": row["analysis_class"],
    }
    output.update(hit)
    return output


def main() -> None:
    for name, path in INPUTS.items():
        if not path.is_file():
            raise AuditError(f"Missing input {name}: {path}")
    observed_hashes = {name: sha256(path) for name, path in INPUTS.items()}
    for name, expected in EXPECTED_FROZEN_SHA256.items():
        if observed_hashes[name] != expected:
            raise AuditError(
                f"Frozen-input checksum mismatch for {name}: expected {expected}; "
                f"observed {observed_hashes[name]}"
            )
    if shutil.which("micromamba") is None:
        raise AuditError("micromamba is not available on PATH")
    diamond_version = run([*DIAMOND, "version"])

    _, integration_rows = read_tsv(INPUTS["plant_integration"])
    _, map_rows = read_tsv(INPUTS["plant_protein_map"])
    top_fields, top_rows = read_tsv(INPUTS["plant_top50"])
    _, prior_best_rows = read_tsv(INPUTS["prior_top50_best_hits"])
    integration = index_unique(integration_rows, "gene_id", "plant integration")
    protein_map = index_unique(map_rows, "gene_id", "plant protein map")
    top50 = index_unique(top_rows, "gene_id", "frozen plant top-50")
    prior_best = index_unique(prior_best_rows, "query_id", "prior top-50 best hits")
    if len(integration) != 41330 or len(protein_map) != 41330:
        raise AuditError("Expected 41,330 plant integrated genes and protein-map rows")
    if set(integration) != set(protein_map):
        raise AuditError("Plant integration and representative-protein map gene sets differ")
    if len(top50) != 50:
        raise AuditError(f"Expected 50 frozen plant top candidates; found {len(top50)}")

    plant_fasta = read_fasta(INPUTS["plant_representative_fasta"])
    insect_fasta = read_fasta(INPUTS["insect_representative_fasta"])
    if len(plant_fasta) != 41330 or len(insect_fasta) != 13455:
        raise AuditError("Representative FASTA record counts differ from validated universes")

    gp_rows: list[dict[str, str]] = []
    for gene_id, integrated in integration.items():
        if integrated.get("DESeq2_status") != "GP_up":
            continue
        mapping = protein_map[gene_id]
        if gene_id not in plant_fasta:
            raise AuditError(f"GP-up gene missing from plant representative FASTA: {gene_id}")
        sequence = plant_fasta[gene_id]["sequence"]
        validate_protein(sequence, gene_id)
        header_protein = parse_header_value(plant_fasta[gene_id]["description"], "representative_protein")
        header_length = parse_header_value(plant_fasta[gene_id]["description"], "length_aa")
        if header_protein != mapping["protein_id"]:
            raise AuditError(f"Representative-protein mismatch for {gene_id}")
        if header_length != mapping["aa_length"] or len(sequence) != int(mapping["aa_length"]):
            raise AuditError(f"Representative-protein length mismatch for {gene_id}")
        gp_rows.append({**integrated, **mapping})
    gp_rows.sort(key=lambda row: row["gene_id"])
    if len(gp_rows) != 4926:
        raise AuditError(f"Expected 4,926 GP-up genes with proteins; found {len(gp_rows)}")
    if len({row["gene_id"] for row in gp_rows}) != len(gp_rows):
        raise AuditError("Duplicate GP-up protein-bearing gene IDs")

    write_fasta(OUTPUTS["gpup_query_fasta"], gp_rows, plant_fasta)
    with tempfile.TemporaryDirectory(prefix="gpup_insect_provenance_") as temporary_name:
        temporary = Path(temporary_name)
        insect_db = temporary / "Qerythrinae_representative"
        raw_hits = temporary / "plant_GPup_vs_Qerythrinae.raw.tsv"
        diamond_makedb(INPUTS["insect_representative_fasta"], insect_db)
        diamond_blastp(OUTPUTS["gpup_query_fasta"], insect_db, raw_hits)
        best_by_gene = choose_best(parse_hits(raw_hits), [row["gene_id"] for row in gp_rows])

    full_rows = [add_context(best_by_gene[row["gene_id"]], row) for row in gp_rows]
    full_fields = [
        "gene_id", "representative_protein_id", "representative_protein_length_aa",
        "baseMean", "log2FoldChange_GP_vs_UP", "padj", "DESeq2_status", "analysis_class",
        *HIT_FIELDS, *DERIVED_FIELDS,
    ]
    write_tsv(OUTPUTS["gpup_best"], full_fields, full_rows)
    near_rows = [
        row for row in full_rows if row["Qerythrinae_similarity_review_class"] == "insect_like_near_identity_review"
    ]
    write_tsv(OUTPUTS["near_identity"], full_fields, near_rows)

    top50_annotated: list[dict[str, object]] = []
    top50_protein_count = 0
    prior_class_equivalent = {
        "insect_like_near_identity_review": "near_identity_cross_system_review",
        "insect_like_substantial_homology_review": "substantial_cross_system_homology",
        "partial_homology_to_Qerythrinae_reference": "partial_cross_system_homology",
        "weak_or_short_match_to_Qerythrinae_reference": "weak_or_short_cross_system_match",
        "no_match_to_Qerythrinae_reference_at_evalue_1e-5": "no_match_at_evalue_1e-5",
    }
    for gene_id, top in sorted(top50.items(), key=lambda item: int(item[1]["discovery_rank"])):
        combined = dict(top)
        if top.get("representative_protein_available") == "TRUE":
            top50_protein_count += 1
            hit = best_by_gene.get(gene_id)
            if hit is None:
                raise AuditError(f"Top-50 protein missing from GP-up provenance scan: {gene_id}")
            combined.update({f"provenance_{key}": value for key, value in hit.items()})
            old = prior_best.get(gene_id)
            if old is None:
                raise AuditError(f"Top-50 protein missing from prior screen: {gene_id}")
            if old.get("subject_id", "") != hit.get("subject_id", ""):
                raise AuditError(
                    f"Top-50 cross-screen inconsistency for {gene_id} (subject_id): "
                    f"prior={old.get('subject_id', '')}; full-GP scan={hit.get('subject_id', '')}"
                )
            expected_prior_class = prior_class_equivalent[
                hit["Qerythrinae_similarity_review_class"]
            ]
            if old.get("cross_system_similarity_class", "") != expected_prior_class:
                raise AuditError(
                    f"Top-50 cross-screen inconsistency for {gene_id} (similarity class): "
                    f"prior={old.get('cross_system_similarity_class', '')}; "
                    f"full-GP scan={expected_prior_class}"
                )
            combined["plant_provenance_interpretation"] = hit["Qerythrinae_similarity_review_class"]
        else:
            combined["plant_provenance_interpretation"] = "no_predicted_protein_transcript_level_review_required"
            for key in HIT_FIELDS + DERIVED_FIELDS:
                combined[f"provenance_{key}"] = ""
        top50_annotated.append(combined)
    if top50_protein_count != 44:
        raise AuditError(f"Expected 44 protein-bearing top-50 genes; found {top50_protein_count}")
    top50_fields = [
        *top_fields,
        "plant_provenance_interpretation",
        *[f"provenance_{key}" for key in HIT_FIELDS + DERIVED_FIELDS],
    ]
    write_tsv(OUTPUTS["top50_annotated"], top50_fields, top50_annotated)

    gp_classes = Counter(str(row["Qerythrinae_similarity_review_class"]) for row in full_rows)
    top_classes = Counter(str(row["plant_provenance_interpretation"]) for row in top50_annotated)
    class_order = [
        "insect_like_near_identity_review",
        "insect_like_substantial_homology_review",
        "partial_homology_to_Qerythrinae_reference",
        "weak_or_short_match_to_Qerythrinae_reference",
        "no_match_to_Qerythrinae_reference_at_evalue_1e-5",
        "no_predicted_protein_transcript_level_review_required",
    ]
    write_tsv(
        OUTPUTS["gpup_summary"],
        ["candidate_set", "provenance_class", "candidate_genes", "percent_of_set"],
        [
            {
                "candidate_set": "all_GP_up_genes_with_representative_protein",
                "provenance_class": label,
                "candidate_genes": str(gp_classes[label]),
                "percent_of_set": format(100 * gp_classes[label] / len(full_rows), ".6g"),
            }
            for label in class_order
            if label in gp_classes
        ],
    )
    write_tsv(
        OUTPUTS["top50_summary"],
        ["candidate_set", "provenance_class", "candidate_genes", "percent_of_set"],
        [
            {
                "candidate_set": "frozen_GP_top50",
                "provenance_class": label,
                "candidate_genes": str(top_classes[label]),
                "percent_of_set": format(100 * top_classes[label] / 50, ".6g"),
            }
            for label in class_order
            if label in top_classes
        ],
    )

    parameter_rows = [
        {"parameter": "DIAMOND_version", "value": diamond_version},
        {"parameter": "DIAMOND_mode", "value": "blastp --sensitive"},
        {"parameter": "evalue_limit", "value": str(EVALUE_LIMIT)},
        {"parameter": "maximum_target_sequences", "value": str(MAX_TARGET_SEQS)},
        {"parameter": "threads", "value": str(THREADS)},
        {"parameter": "GP_up_protein_bearing_genes", "value": "4926"},
        {"parameter": "insect_like_near_identity_review", "value": "identity>=90%; query span coverage>=80%; subject span coverage>=80%"},
        {"parameter": "insect_like_substantial_homology_review", "value": "evalue<=1e-10; identity>=40%; query and subject span coverage>=60%"},
        {"parameter": "interpretive_boundary", "value": "Reference-proteome similarity is a provenance review flag, not definitive taxonomic assignment or evidence of contamination"},
    ]
    write_tsv(OUTPUTS["parameters"], ["parameter", "value"], parameter_rows)

    for name, path in INPUTS.items():
        if sha256(path) != observed_hashes[name]:
            raise AuditError(f"Input changed during analysis: {name}")

    report = [
        "===== PLANT GP-UP INSECT-PROVENANCE AUDIT =====",
        f"Discovery root: {ROOT}",
        "",
        "Scope:",
        "  Protein-bearing GP-up genes were screened against the Q. erythrinae representative proteome.",
        "  The frozen GP top-50 ranking was not changed.",
        "",
        f"Plant representative proteins: {len(plant_fasta)}",
        f"Insect representative proteins: {len(insect_fasta)}",
        f"GP-up genes with representative proteins audited: {len(full_rows)}",
        f"Frozen GP top-50 genes: {len(top50_annotated)}",
        f"Frozen GP top-50 genes with predicted protein: {top50_protein_count}",
        f"Frozen GP top-50 genes requiring transcript-level audit: {50 - top50_protein_count}",
        "",
        "All GP-up protein-bearing genes by Q. erythrinae similarity review class:",
    ]
    for label in class_order:
        if label in gp_classes:
            report.append(f"  {label}: {gp_classes[label]}")
    report.extend(["", "Frozen GP top-50 by provenance interpretation:"])
    for label in class_order:
        if label in top_classes:
            report.append(f"  {label}: {top_classes[label]}")
    report.extend([
        "",
        "Interpretive boundary:",
        "  Strong similarity to the insect reference proteome is a review flag, not a taxonomic assignment.",
        "  It does not independently prove contamination, horizontal transfer, molecular mimicry, or functional equivalence.",
        "  No strong match to the Q. erythrinae reference proteome is not proof of plant origin.",
        "  Six frozen top-50 genes have no predicted representative protein and remain unresolved pending a transcript-level audit.",
        "",
        "Frozen-input SHA-256 verification after output generation: PASS",
        "Errors detected: 0",
        "",
    ])
    for path in OUTPUTS.values():
        report.append(f"Output: {path}")
    report.append("PLANT GP-UP INSECT-PROVENANCE AUDIT: PASS")
    OUTPUTS["audit"].parent.mkdir(parents=True, exist_ok=True)
    temporary_audit = OUTPUTS["audit"].with_name(OUTPUTS["audit"].name + ".tmp")
    temporary_audit.write_text("\n".join(report) + "\n", encoding="utf-8")
    os.replace(temporary_audit, OUTPUTS["audit"])
    print("\n".join(report))


if __name__ == "__main__":
    try:
        main()
    except (AuditError, subprocess.CalledProcessError) as error:
        print("===== PLANT GP-UP INSECT-PROVENANCE AUDIT =====", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError):
            print(f"ERROR: command failed: {' '.join(error.cmd)}", file=sys.stderr)
            if error.stderr:
                print(error.stderr.strip(), file=sys.stderr)
        else:
            print(f"ERROR: {error}", file=sys.stderr)
        print("PLANT GP-UP INSECT-PROVENANCE AUDIT: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
