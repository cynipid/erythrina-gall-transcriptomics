#!/usr/bin/env python3
"""Confirm selected plant Trinity transcripts against the Q. erythrinae CDS set.

The selected set is the union of (i) protein-level near-identity GP-up review
genes and (ii) all genes in the frozen GP top-50. This script is an
independent nucleotide-level provenance audit. It does not filter counts,
rerun DESeq2, or modify discovery rankings.
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
from collections import Counter, defaultdict
from pathlib import Path


HOME = Path.home()
ROOT = HOME / "projects/qe_cross_species_discovery"

INPUTS = {
    "proposed_mask": ROOT / "04_cross_species_integration/results/plant_host_reanalysis_proposed_insect_like_sequence_mask.tsv",
    "top50": ROOT / "03_plant_discovery/results/plant_GP_top50.tsv",
    "top50_status": ROOT / "03_plant_discovery/results/plant_GP_top50_provenance_reanalysis_status.tsv",
    "trinity_assembly": HOME / "projects/qe_plant_busco_embryophyta/input/Trinity_plant.fasta",
    "insect_cds": HOME / "projects/qe_genome_guided/01_reference/complete.cds.fna",
}

EXPECTED_TOP50_SHA256 = "1d4aff7ecab264a97c5c7944a74919293c6f8514aded0f7b56e7cf8a672aa4bb"

OUTPUTS = {
    "selected_transcripts": ROOT / "03_plant_discovery/results/plant_GPup_insectlike_plus_top50_transcripts_for_CDS_audit.fasta",
    "transcript_mappings": ROOT / "04_cross_species_integration/results/plant_selected_transcripts_vs_Qerythrinae_CDS_mappings.tsv",
    "gene_best": ROOT / "04_cross_species_integration/results/plant_GPup_insectlike_plus_top50_CDS_provenance_best_by_gene.tsv",
    "top50_nucleotide": ROOT / "03_plant_discovery/results/plant_GP_top50_nucleotide_provenance_audit.tsv",
    "summary": ROOT / "04_cross_species_integration/results/plant_GPup_insectlike_CDS_provenance_summary.tsv",
    "parameters": ROOT / "00_project_control/07_plant_transcript_insect_CDS_provenance_parameters.tsv",
    "audit": ROOT / "08_logs/07_plant_transcript_insect_CDS_provenance_audit.txt",
}

MINIMAP = ["micromamba", "run", "-n", "qe_provenance", "minimap2"]
THREADS = 8
TRINITY_GENE = re.compile(r"^(TRINITY_DN\d+_c\d+_g\d+)_i\d+$")
PAF_FIELDS = [
    "query_transcript_id", "query_length_nt", "query_start", "query_end",
    "strand", "target_CDS_id", "target_length_nt", "target_start", "target_end",
    "matching_bases", "alignment_block_length", "mapping_quality",
]
DERIVED_FIELDS = [
    "percent_identity", "query_span_coverage_percent", "target_span_coverage_percent",
    "nucleotide_provenance_review_class", "high_identity_Qerythrinae_CDS_review",
]


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


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditError(f"Missing TSV header: {path}")
        fields = [clean(field) for field in reader.fieldnames]
        return fields, [{clean(k): clean(v) for k, v in row.items()} for row in reader]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def index_unique(rows: list[dict[str, str]], key: str, label: str) -> dict[str, dict[str, str]]:
    result = {}
    for row in rows:
        identifier = row.get(key, "")
        if not identifier or identifier in result:
            raise AuditError(f"Invalid or duplicate {key} in {label}: {identifier!r}")
        result[identifier] = row
    return result


def gene_from_transcript(transcript_id: str) -> str:
    match = TRINITY_GENE.fullmatch(transcript_id)
    if match is None:
        raise AuditError(f"Unrecognised Trinity transcript identifier: {transcript_id}")
    return match.group(1)


def extract_selected_transcripts(
    assembly: Path, selected_genes: set[str], output: Path
) -> dict[str, list[str]]:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    gene_transcripts: dict[str, list[str]] = defaultdict(list)
    transcript_id: str | None = None
    transcript_gene: str | None = None
    sequence_parts: list[str] = []

    def store() -> None:
        if transcript_id is None or transcript_gene not in selected_genes:
            return
        sequence = "".join(sequence_parts).upper()
        if not sequence:
            raise AuditError(f"Empty selected transcript: {transcript_id}")
        invalid = sorted(set(sequence) - set("ACGTRYSWKMBDHVN"))
        if invalid:
            raise AuditError(f"Invalid nucleotide characters in {transcript_id}: {''.join(invalid)}")
        out_handle.write(f">{transcript_id} gene_id={transcript_gene}\n")
        for start in range(0, len(sequence), 80):
            out_handle.write(sequence[start:start + 80] + "\n")
        gene_transcripts[transcript_gene].append(transcript_id)

    with temporary.open("w", encoding="utf-8", newline="\n") as out_handle:
        with assembly.open("r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    store()
                    transcript_id = line[1:].split(maxsplit=1)[0]
                    transcript_gene = gene_from_transcript(transcript_id)
                    sequence_parts = []
                else:
                    if transcript_id is None:
                        raise AuditError("Sequence occurred before first Trinity FASTA header")
                    sequence_parts.append(line)
            store()
    os.replace(temporary, output)
    return dict(gene_transcripts)


def run(command: list[str], stdout_path: Path | None = None) -> str:
    if stdout_path is None:
        completed = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return completed.stdout.strip()
    with stdout_path.open("w", encoding="utf-8", newline="\n") as handle:
        subprocess.run(command, check=True, text=True, stdout=handle, stderr=subprocess.PIPE)
    return ""


def classify(identity: float, query_coverage: float, alignment_length: int) -> tuple[str, str]:
    if identity >= 95.0 and query_coverage >= 60.0 and alignment_length >= 300:
        return "high_identity_Qerythrinae_CDS_review", "TRUE"
    if identity >= 90.0 and query_coverage >= 40.0 and alignment_length >= 200:
        return "substantial_Qerythrinae_CDS_similarity_review", "FALSE"
    if identity >= 70.0 and alignment_length >= 150:
        return "partial_Qerythrinae_CDS_similarity", "FALSE"
    return "weak_or_short_CDS_match", "FALSE"


def parse_paf(path: Path) -> list[dict[str, str]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, values in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if len(values) < 12:
                raise AuditError(f"Malformed PAF row at {path}:{line_number}")
            row = dict(zip(PAF_FIELDS, values[:12]))
            qlen, tlen = int(row["query_length_nt"]), int(row["target_length_nt"])
            qspan = int(row["query_end"]) - int(row["query_start"])
            tspan = int(row["target_end"]) - int(row["target_start"])
            matches, length = int(row["matching_bases"]), int(row["alignment_block_length"])
            identity, qcov, tcov = 100 * matches / length, 100 * qspan / qlen, 100 * tspan / tlen
            status, strict = classify(identity, qcov, length)
            row.update({
                "percent_identity": format(identity, ".10g"),
                "query_span_coverage_percent": format(qcov, ".10g"),
                "target_span_coverage_percent": format(tcov, ".10g"),
                "nucleotide_provenance_review_class": status,
                "high_identity_Qerythrinae_CDS_review": strict,
            })
            rows.append(row)
    return rows


def best_per_gene(
    mappings: list[dict[str, str]], gene_transcripts: dict[str, list[str]], selected_genes: set[str]
) -> dict[str, dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mappings:
        grouped[gene_from_transcript(row["query_transcript_id"])].append(row)
    result = {}
    for gene in sorted(selected_genes):
        candidates = grouped.get(gene, [])
        if not candidates:
            result[gene] = {
                "best_query_transcript_id": "", "target_CDS_id": "", "percent_identity": "",
                "query_span_coverage_percent": "", "target_span_coverage_percent": "",
                "alignment_block_length": "", "mapping_quality": "",
                "nucleotide_provenance_review_class": "no_CDS_match_reported_by_minimap2",
                "high_identity_Qerythrinae_CDS_review": "FALSE",
            }
            continue
        best = sorted(
            candidates,
            key=lambda row: (
                row["high_identity_Qerythrinae_CDS_review"] != "TRUE",
                -float(row["query_span_coverage_percent"]),
                -float(row["percent_identity"]),
                -int(row["alignment_block_length"]),
                -int(row["mapping_quality"]),
            ),
        )[0]
        result[gene] = {
            "best_query_transcript_id": best["query_transcript_id"],
            "target_CDS_id": best["target_CDS_id"],
            "percent_identity": best["percent_identity"],
            "query_span_coverage_percent": best["query_span_coverage_percent"],
            "target_span_coverage_percent": best["target_span_coverage_percent"],
            "alignment_block_length": best["alignment_block_length"],
            "mapping_quality": best["mapping_quality"],
            "nucleotide_provenance_review_class": best["nucleotide_provenance_review_class"],
            "high_identity_Qerythrinae_CDS_review": best["high_identity_Qerythrinae_CDS_review"],
        }
    return result


def main() -> None:
    for name, path in INPUTS.items():
        if not path.is_file():
            raise AuditError(f"Missing input {name}: {path}")
    hashes = {name: sha256(path) for name, path in INPUTS.items()}
    if hashes["top50"] != EXPECTED_TOP50_SHA256:
        raise AuditError("Frozen plant top-50 checksum mismatch")
    if shutil.which("micromamba") is None:
        raise AuditError("micromamba is not available on PATH")
    minimap_version = run([*MINIMAP, "--version"])

    _, mask_rows = read_tsv(INPUTS["proposed_mask"])
    top_fields, top_rows = read_tsv(INPUTS["top50"])
    _, top_status_rows = read_tsv(INPUTS["top50_status"])
    mask = index_unique(mask_rows, "gene_id", "proposed host-reanalysis mask")
    top = index_unique(top_rows, "gene_id", "frozen top-50")
    top_status = index_unique(top_status_rows, "gene_id", "top-50 provenance status")
    if len(mask) != 41330 or len(top) != 50 or set(top) != set(top_status):
        raise AuditError("Mask or frozen top-50 table membership does not match expected values")

    strict_gp = {
        gene for gene, row in mask.items()
        if row["DESeq2_status"] == "GP_up"
        and row["Qerythrinae_similarity_review_class"] == "insect_like_near_identity_review"
    }
    selected = strict_gp | set(top)
    if len(strict_gp) != 2534:
        raise AuditError(f"Expected 2,534 strict GP-up protein review genes; observed {len(strict_gp)}")
    if len(selected) != 2543:
        raise AuditError(f"Expected 2,543 union genes; observed {len(selected)}")

    gene_transcripts = extract_selected_transcripts(
        INPUTS["trinity_assembly"], selected, OUTPUTS["selected_transcripts"]
    )
    missing = sorted(selected - set(gene_transcripts))
    if missing:
        raise AuditError("Selected genes missing from Trinity assembly: " + ", ".join(missing))
    top_no_protein = {
        gene for gene, row in top.items() if row.get("representative_protein_available") != "TRUE"
    }
    if len(top_no_protein) != 6:
        raise AuditError(f"Expected six top-50 genes without a predicted protein; observed {len(top_no_protein)}")

    with tempfile.TemporaryDirectory(prefix="plant_transcript_CDS_provenance_") as temporary_name:
        paf = Path(temporary_name) / "selected_transcripts_vs_Qerythrinae_CDS.paf"
        run([
            *MINIMAP, "-x", "asm5", "--secondary=no", "-t", str(THREADS),
            str(INPUTS["insect_cds"]), str(OUTPUTS["selected_transcripts"]),
        ], stdout_path=paf)
        mappings = parse_paf(paf)

    mapping_fields = [*PAF_FIELDS, *DERIVED_FIELDS]
    write_tsv(OUTPUTS["transcript_mappings"], mapping_fields, mappings)
    best = best_per_gene(mappings, gene_transcripts, selected)

    gene_rows = []
    for gene in sorted(selected):
        mask_row = mask.get(gene, {})
        frozen_top_row = top.get(gene, {})
        reasons = []
        if gene in strict_gp:
            reasons.append("GP_up_protein_near_identity_review")
        if gene in top:
            reasons.append("frozen_GP_top50")
        gene_rows.append({
            "gene_id": gene,
            "selection_reason": ";".join(reasons),
            "DESeq2_status": mask_row.get("DESeq2_status", frozen_top_row.get("DESeq2_status", "")),
            "protein_provenance_review_class": mask_row.get(
                "Qerythrinae_similarity_review_class",
                "no_predicted_representative_protein",
            ),
            "Trinity_transcripts_extracted": str(len(gene_transcripts[gene])),
            **best[gene],
        })
    best_fields = [
        "gene_id", "selection_reason", "DESeq2_status", "protein_provenance_review_class",
        "Trinity_transcripts_extracted", "best_query_transcript_id", "target_CDS_id",
        "percent_identity", "query_span_coverage_percent", "target_span_coverage_percent",
        "alignment_block_length", "mapping_quality", "nucleotide_provenance_review_class",
        "high_identity_Qerythrinae_CDS_review",
    ]
    write_tsv(OUTPUTS["gene_best"], best_fields, gene_rows)

    top_nucleotide_rows = []
    for gene, row in sorted(top.items(), key=lambda item: int(item[1]["discovery_rank"])):
        entry = dict(row)
        entry.update({f"nucleotide_{key}": value for key, value in best[gene].items()})
        entry["top50_protein_provenance_status"] = top_status[gene]["full_proteome_provenance_status"]
        top_nucleotide_rows.append(entry)
    top_nucleotide_fields = [
        *top_fields, "top50_protein_provenance_status",
        *[f"nucleotide_{key}" for key in best_fields if key != "gene_id"],
    ]
    write_tsv(OUTPUTS["top50_nucleotide"], top_nucleotide_fields, top_nucleotide_rows)

    strict_rows = [row for row in gene_rows if row["gene_id"] in strict_gp]
    top_rows_gene = [row for row in gene_rows if row["gene_id"] in top]
    summary_sets = {
        "strict_GP_up_protein_near_identity_review": strict_rows,
        "frozen_GP_top50": top_rows_gene,
        "frozen_GP_top50_without_predicted_protein": [row for row in top_rows_gene if row["gene_id"] in top_no_protein],
    }
    labels = [
        "high_identity_Qerythrinae_CDS_review", "substantial_Qerythrinae_CDS_similarity_review",
        "partial_Qerythrinae_CDS_similarity", "weak_or_short_CDS_match",
        "no_CDS_match_reported_by_minimap2",
    ]
    summary_rows = []
    for name, rows in summary_sets.items():
        counts = Counter(row["nucleotide_provenance_review_class"] for row in rows)
        for label in labels:
            summary_rows.append({
                "candidate_set": name, "nucleotide_provenance_review_class": label,
                "candidate_genes": str(counts[label]),
                "percent_of_candidate_set": format(100 * counts[label] / len(rows), ".6g"),
            })
    write_tsv(OUTPUTS["summary"], ["candidate_set", "nucleotide_provenance_review_class", "candidate_genes", "percent_of_candidate_set"], summary_rows)

    parameter_rows = [
        {"parameter": "minimap2_version", "value": minimap_version},
        {"parameter": "preset", "value": "asm5"},
        {"parameter": "secondary_alignments", "value": "disabled"},
        {"parameter": "threads", "value": str(THREADS)},
        {"parameter": "strict_protein_review_genes", "value": "2534"},
        {"parameter": "union_with_frozen_top50_genes", "value": "2543"},
        {"parameter": "high_identity_CDS_review", "value": "identity>=95%; query span coverage>=60%; alignment block length>=300 nt"},
        {"parameter": "interpretive_boundary", "value": "Nucleotide similarity validates sequence concordance with the Q. erythrinae CDS reference but does not by itself identify the physical source of RNA"},
    ]
    for name, checksum in sorted(hashes.items()):
        parameter_rows.append({"parameter": f"SHA256_{name}", "value": checksum})
    write_tsv(OUTPUTS["parameters"], ["parameter", "value"], parameter_rows)

    for name, path in INPUTS.items():
        if sha256(path) != hashes[name]:
            raise AuditError(f"Input changed during analysis: {name}")
    strict_counts = Counter(row["nucleotide_provenance_review_class"] for row in strict_rows)
    top_counts = Counter(row["nucleotide_provenance_review_class"] for row in top_rows_gene)
    no_protein_counts = Counter(row["nucleotide_provenance_review_class"] for row in summary_sets["frozen_GP_top50_without_predicted_protein"])
    report = [
        "===== PLANT TRANSCRIPT–INSECT CDS PROVENANCE AUDIT =====", f"Discovery root: {ROOT}", "",
        f"Trinity assembly: {INPUTS['trinity_assembly']}", f"Insect CDS reference: {INPUTS['insect_cds']}",
        f"Strict protein-level GP-up review genes: {len(strict_rows)}", f"Frozen top-50 genes: {len(top_rows_gene)}",
        f"Selected genes in union: {len(selected)}", f"Selected Trinity transcripts extracted: {sum(len(items) for items in gene_transcripts.values())}", "",
        "Strict protein-level GP-up review genes by nucleotide provenance class:",
    ]
    report.extend(f"  {label}: {strict_counts[label]}" for label in labels)
    report.extend(["", "Frozen top-50 genes by nucleotide provenance class:"])
    report.extend(f"  {label}: {top_counts[label]}" for label in labels)
    report.extend(["", "Top-50 genes without predicted protein by nucleotide provenance class:"])
    report.extend(f"  {label}: {no_protein_counts[label]}" for label in labels)
    report.extend([
        "", "Interpretive boundary:",
        "  High-identity nucleotide mapping is an independent sequence-provenance confirmation, not proof of the physical RNA source.",
        "  It does not establish dissection error, library carryover, or a causal biological mechanism.",
        "  This audit does not change the frozen rankings or any DESeq2 result.",
        "", "Input SHA-256 verification after output generation: PASS", "Errors detected: 0", "",
    ])
    for path in OUTPUTS.values():
        report.append(f"Output: {path}")
    report.append("PLANT TRANSCRIPT–INSECT CDS PROVENANCE AUDIT: PASS")
    OUTPUTS["audit"].parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUTS["audit"].with_name(OUTPUTS["audit"].name + ".tmp")
    temporary.write_text("\n".join(report) + "\n", encoding="utf-8")
    os.replace(temporary, OUTPUTS["audit"])
    print("\n".join(report))


if __name__ == "__main__":
    try:
        main()
    except (AuditError, subprocess.CalledProcessError) as error:
        print("===== PLANT TRANSCRIPT–INSECT CDS PROVENANCE AUDIT =====", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError):
            print(f"ERROR: command failed: {' '.join(error.cmd)}", file=sys.stderr)
            if error.stderr:
                print(error.stderr.strip(), file=sys.stderr)
        else:
            print(f"ERROR: {error}", file=sys.stderr)
        print("PLANT TRANSCRIPT–INSECT CDS PROVENANCE AUDIT: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
