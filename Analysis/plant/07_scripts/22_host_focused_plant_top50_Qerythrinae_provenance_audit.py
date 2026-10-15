#!/usr/bin/env python3
"""Audit the frozen host-focused plant top 50 against Q. erythrinae sequences.

The script reuses the validated genome-wide protein-similarity table and runs
only a nucleotide comparison of transcripts belonging to the frozen top 50.
It performs no filtering, DESeq2 fitting, reranking, or taxonomic assignment.
"""

from __future__ import annotations

import argparse
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


EXPECTED_TOP50 = 50
EXPECTED_REVIEW15 = 15
EXPECTED_PROTEIN_BEARING_TOP50 = 32
EXPECTED_NO_PROTEIN_TOP50 = 18
EXPECTED_PROTEIN_UNIVERSE = 41_330
EXPECTED_MASK_UNIVERSE = 41_330
EXPECTED_CONFIRMED_MASK = 3_948
THREADS = 8
MINIMAP = ["micromamba", "run", "-n", "qe_provenance", "minimap2"]
TRINITY_GENE = re.compile(r"^(TRINITY_DN\d+_c\d+_g\d+)_i\d+$")

PAF_FIELDS = [
    "query_transcript_id",
    "query_length_nt",
    "query_start",
    "query_end",
    "strand",
    "target_CDS_id",
    "target_length_nt",
    "target_start",
    "target_end",
    "matching_bases",
    "alignment_block_length",
    "mapping_quality",
]
DERIVED_FIELDS = [
    "percent_identity",
    "query_span_coverage_percent",
    "target_span_coverage_percent",
    "nucleotide_provenance_review_class",
    "high_identity_Qerythrinae_CDS_review",
]

NUCLEOTIDE_CLASSES = [
    "high_identity_Qerythrinae_CDS_review",
    "substantial_Qerythrinae_CDS_similarity_review",
    "partial_Qerythrinae_CDS_similarity",
    "weak_or_short_CDS_match",
    "no_CDS_match_reported_by_minimap2",
]


class AuditError(RuntimeError):
    pass


def clean(value: object) -> str:
    return "" if value is None else str(value).strip().lstrip("\ufeff")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise AuditError(f"Required input is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditError(f"Missing TSV header: {path}")
        fields = [clean(field) for field in reader.fieldnames]
        rows = [{clean(key): clean(value) for key, value in row.items()} for row in reader]
    return fields, rows


def require_columns(path: Path, fields: list[str], required: set[str]) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise AuditError(f"Missing columns in {path}: {', '.join(missing)}")


def index_unique(
    rows: list[dict[str, str]], key: str, label: str
) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = row.get(key, "")
        if not identifier or identifier in result:
            raise AuditError(f"Invalid or duplicate {key} in {label}: {identifier!r}")
        result[identifier] = row
    return result


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

    with temporary.open("w", encoding="utf-8", newline="\n") as out_handle:

        def store() -> None:
            if transcript_id is None or transcript_gene not in selected_genes:
                return
            sequence = "".join(sequence_parts).upper()
            if not sequence:
                raise AuditError(f"Empty selected transcript: {transcript_id}")
            invalid = sorted(set(sequence).difference("ACGTRYSWKMBDHVN"))
            if invalid:
                raise AuditError(
                    f"Invalid nucleotide characters in {transcript_id}: {''.join(invalid)}"
                )
            out_handle.write(f">{transcript_id} gene_id={transcript_gene}\n")
            for start in range(0, len(sequence), 80):
                out_handle.write(sequence[start : start + 80] + "\n")
            gene_transcripts[transcript_gene].append(transcript_id)

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
        completed = subprocess.run(
            command,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return completed.stdout.strip()
    with stdout_path.open("w", encoding="utf-8", newline="\n") as handle:
        subprocess.run(
            command,
            check=True,
            text=True,
            stdout=handle,
            stderr=subprocess.PIPE,
        )
    return ""


def nucleotide_class(
    identity: float, query_coverage: float, alignment_length: int
) -> tuple[str, str]:
    # Frozen definitions from scripts 07 and 08.
    if identity >= 95.0 and query_coverage >= 60.0 and alignment_length >= 300:
        return "high_identity_Qerythrinae_CDS_review", "TRUE"
    if identity >= 90.0 and query_coverage >= 40.0 and alignment_length >= 200:
        return "substantial_Qerythrinae_CDS_similarity_review", "FALSE"
    if identity >= 70.0 and alignment_length >= 150:
        return "partial_Qerythrinae_CDS_similarity", "FALSE"
    return "weak_or_short_CDS_match", "FALSE"


def parse_paf(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, values in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if len(values) < 12:
                raise AuditError(f"Malformed PAF row at {path}:{line_number}")
            row = dict(zip(PAF_FIELDS, values[:12]))
            query_length = int(row["query_length_nt"])
            target_length = int(row["target_length_nt"])
            query_span = int(row["query_end"]) - int(row["query_start"])
            target_span = int(row["target_end"]) - int(row["target_start"])
            matches = int(row["matching_bases"])
            alignment_length = int(row["alignment_block_length"])
            if query_length <= 0 or target_length <= 0 or alignment_length <= 0:
                raise AuditError(f"Invalid PAF lengths at {path}:{line_number}")
            identity = 100.0 * matches / alignment_length
            query_coverage = 100.0 * query_span / query_length
            target_coverage = 100.0 * target_span / target_length
            label, strict = nucleotide_class(identity, query_coverage, alignment_length)
            row.update(
                {
                    "percent_identity": format(identity, ".10g"),
                    "query_span_coverage_percent": format(query_coverage, ".10g"),
                    "target_span_coverage_percent": format(target_coverage, ".10g"),
                    "nucleotide_provenance_review_class": label,
                    "high_identity_Qerythrinae_CDS_review": strict,
                }
            )
            rows.append(row)
    return rows


def best_per_gene(
    mappings: list[dict[str, str]], selected_genes: set[str]
) -> dict[str, dict[str, str]]:
    class_priority = {
        "high_identity_Qerythrinae_CDS_review": 4,
        "substantial_Qerythrinae_CDS_similarity_review": 3,
        "partial_Qerythrinae_CDS_similarity": 2,
        "weak_or_short_CDS_match": 1,
    }
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mappings:
        grouped[gene_from_transcript(row["query_transcript_id"])].append(row)

    result: dict[str, dict[str, str]] = {}
    for gene in sorted(selected_genes):
        candidates = grouped.get(gene, [])
        if not candidates:
            result[gene] = {
                "best_query_transcript_id": "",
                "target_CDS_id": "",
                "percent_identity": "",
                "query_span_coverage_percent": "",
                "target_span_coverage_percent": "",
                "alignment_block_length": "",
                "mapping_quality": "",
                "nucleotide_provenance_review_class": (
                    "no_CDS_match_reported_by_minimap2"
                ),
                "high_identity_Qerythrinae_CDS_review": "FALSE",
            }
            continue
        best = sorted(
            candidates,
            key=lambda row: (
                -class_priority[row["nucleotide_provenance_review_class"]],
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
            "nucleotide_provenance_review_class": best[
                "nucleotide_provenance_review_class"
            ],
            "high_identity_Qerythrinae_CDS_review": best[
                "high_identity_Qerythrinae_CDS_review"
            ],
        }
    return result


def overall_review_class(protein_class: str, nucleotide_class_name: str) -> str:
    if nucleotide_class_name in {
        "high_identity_Qerythrinae_CDS_review",
        "substantial_Qerythrinae_CDS_similarity_review",
    }:
        return "additional_Qerythrinae_like_sequence_review"
    if protein_class == "insect_like_near_identity_review":
        return "protein_near_identity_without_strong_nucleotide_confirmation_review"
    if (
        protein_class == "insect_like_substantial_homology_review"
        or nucleotide_class_name == "partial_Qerythrinae_CDS_similarity"
    ):
        return "cross_species_homology_review"
    return "no_strong_Qerythrinae_similarity_in_current_screen"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--discovery-root",
        type=Path,
        default=Path.home() / "projects" / "qe_cross_species_discovery",
    )
    parser.add_argument(
        "--trinity-assembly",
        type=Path,
        default=(
            Path.home()
            / "projects/qe_plant_busco_embryophyta/input/Trinity_plant.fasta"
        ),
    )
    parser.add_argument(
        "--insect-cds",
        type=Path,
        default=(
            Path.home()
            / "projects/qe_genome_guided/01_reference/complete.cds.fna"
        ),
    )
    args = parser.parse_args()

    root = args.discovery_root.resolve()
    host_results = root / "05_host_focused_DESeq2_sensitivity/results"
    control_dir = root / "00_project_control"
    logs_dir = root / "05_host_focused_DESeq2_sensitivity/logs"

    inputs = {
        "top50": host_results / "plant_GP_host_focused_top50.tsv",
        "review15": host_results
        / "plant_GP_host_focused_poorly_annotated_top15.tsv",
        "protein_similarity": root
        / "04_cross_species_integration/results/plant_representative_proteome_vs_Qerythrinae_best_hits.tsv",
        "confirmed_mask": root
        / "04_cross_species_integration/results/plant_confirmed_multilevel_insect_like_sequence_mask.tsv",
        "trinity_assembly": args.trinity_assembly.resolve(),
        "insect_cds": args.insect_cds.resolve(),
    }
    for label, path in inputs.items():
        if not path.is_file():
            raise AuditError(f"Missing input {label}: {path}")

    outputs = {
        "selected_transcripts": host_results
        / "plant_GP_host_focused_top50_transcripts_for_Qerythrinae_CDS_audit.fasta",
        "mappings": host_results
        / "plant_GP_host_focused_top50_vs_Qerythrinae_CDS_mappings.tsv",
        "top50_audit": host_results
        / "plant_GP_host_focused_top50_sequence_provenance_audit.tsv",
        "review15_audit": host_results
        / "plant_GP_host_focused_poorly_annotated_top15_sequence_provenance_audit.tsv",
        "summary": host_results
        / "plant_GP_host_focused_top50_sequence_provenance_summary.tsv",
        "parameters": control_dir
        / "22_host_focused_plant_top50_Qerythrinae_provenance_parameters.tsv",
        "audit": logs_dir
        / "22_host_focused_plant_top50_Qerythrinae_provenance_audit.txt",
    }

    hashes = {label: sha256(path) for label, path in inputs.items()}
    if shutil.which("micromamba") is None:
        raise AuditError("micromamba is not available on PATH")
    minimap_version = run([*MINIMAP, "--version"])

    top_fields, top_rows = read_tsv(inputs["top50"])
    review_fields, review_rows = read_tsv(inputs["review15"])
    protein_fields, protein_rows = read_tsv(inputs["protein_similarity"])
    mask_fields, mask_rows = read_tsv(inputs["confirmed_mask"])

    require_columns(
        inputs["top50"],
        top_fields,
        {
            "gene_id",
            "discovery_rank",
            "annotation_status",
            "protein_prediction_available",
        },
    )
    require_columns(inputs["review15"], review_fields, {"gene_id", "discovery_rank"})
    require_columns(
        inputs["protein_similarity"],
        protein_fields,
        {
            "gene_id",
            "representative_protein_id",
            "Qerythrinae_similarity_review_class",
        },
    )
    require_columns(
        inputs["confirmed_mask"],
        mask_fields,
        {"gene_id", "host_reanalysis_mask_status"},
    )

    top = index_unique(top_rows, "gene_id", "host-focused top 50")
    review = index_unique(review_rows, "gene_id", "host-focused review top 15")
    protein = index_unique(protein_rows, "gene_id", "genome-wide protein similarity")
    mask = index_unique(mask_rows, "gene_id", "confirmed multilevel mask")

    if len(top) != EXPECTED_TOP50:
        raise AuditError(f"Top-50 rows: {len(top)}; expected {EXPECTED_TOP50}")
    if len(review) != EXPECTED_REVIEW15 or not set(review).issubset(top):
        raise AuditError("Poorly annotated/uncharacterized top-15 membership is invalid")
    if len(protein) != EXPECTED_PROTEIN_UNIVERSE:
        raise AuditError(
            f"Protein-similarity universe: {len(protein)}; expected {EXPECTED_PROTEIN_UNIVERSE}"
        )
    if len(mask) != EXPECTED_MASK_UNIVERSE:
        raise AuditError(
            f"Full mask-universe rows: {len(mask)}; expected {EXPECTED_MASK_UNIVERSE}"
        )
    confirmed_mask_genes = {
        gene
        for gene, row in mask.items()
        if row["host_reanalysis_mask_status"]
        == "review_exclude_confirmed_multilevel_insect_like_sequence"
    }
    if len(confirmed_mask_genes) != EXPECTED_CONFIRMED_MASK:
        raise AuditError(
            f"Confirmed exclusion flags: {len(confirmed_mask_genes)}; "
            f"expected {EXPECTED_CONFIRMED_MASK}"
        )
    mask_overlap = sorted(set(top).intersection(confirmed_mask_genes))
    if mask_overlap:
        raise AuditError(
            "Host-focused top 50 overlaps the excluded confirmed mask: "
            + ", ".join(mask_overlap)
        )

    protein_bearing = {
        gene
        for gene, row in top.items()
        if row["protein_prediction_available"] == "TRUE"
    }
    no_protein = set(top).difference(protein_bearing)
    if len(protein_bearing) != EXPECTED_PROTEIN_BEARING_TOP50:
        raise AuditError(
            f"Protein-bearing top-50 genes: {len(protein_bearing)}; "
            f"expected {EXPECTED_PROTEIN_BEARING_TOP50}"
        )
    if len(no_protein) != EXPECTED_NO_PROTEIN_TOP50:
        raise AuditError(
            f"Top-50 genes without proteins: {len(no_protein)}; "
            f"expected {EXPECTED_NO_PROTEIN_TOP50}"
        )
    missing_protein_rows = sorted(protein_bearing.difference(protein))
    unexpected_protein_rows = sorted(no_protein.intersection(protein))
    if missing_protein_rows or unexpected_protein_rows:
        raise AuditError(
            "Protein availability disagrees with the genome-wide similarity table: "
            f"missing={len(missing_protein_rows)}, unexpected={len(unexpected_protein_rows)}"
        )

    gene_transcripts = extract_selected_transcripts(
        inputs["trinity_assembly"], set(top), outputs["selected_transcripts"]
    )
    missing_transcripts = sorted(set(top).difference(gene_transcripts))
    if missing_transcripts:
        raise AuditError(
            "Top-50 genes missing from Trinity assembly: "
            + ", ".join(missing_transcripts)
        )

    with tempfile.TemporaryDirectory(prefix="host_top50_Qe_CDS_") as temporary_name:
        paf = Path(temporary_name) / "top50_vs_Qerythrinae_CDS.paf"
        run(
            [
                *MINIMAP,
                "-x",
                "asm5",
                "--secondary=no",
                "-t",
                str(THREADS),
                str(inputs["insect_cds"]),
                str(outputs["selected_transcripts"]),
            ],
            stdout_path=paf,
        )
        mappings = parse_paf(paf)

    write_tsv(outputs["mappings"], [*PAF_FIELDS, *DERIVED_FIELDS], mappings)
    best = best_per_gene(mappings, set(top))

    appended_fields = [
        "protein_similarity_available",
        "protein_representative_id",
        "protein_Qerythrinae_subject_id",
        "protein_percent_identity",
        "protein_query_span_coverage_percent",
        "protein_subject_span_coverage_percent",
        "protein_evalue",
        "protein_Qerythrinae_similarity_review_class",
        "Trinity_transcripts_extracted",
        "nucleotide_best_query_transcript_id",
        "nucleotide_target_CDS_id",
        "nucleotide_percent_identity",
        "nucleotide_query_span_coverage_percent",
        "nucleotide_target_span_coverage_percent",
        "nucleotide_alignment_block_length",
        "nucleotide_mapping_quality",
        "nucleotide_provenance_review_class",
        "nucleotide_high_identity_Qerythrinae_CDS_review",
        "overall_Qerythrinae_sequence_review_class",
        "confirmed_mask_overlap",
    ]

    audited_rows: list[dict[str, object]] = []
    for source in sorted(top_rows, key=lambda row: int(row["discovery_rank"])):
        gene = source["gene_id"]
        protein_row = protein.get(gene, {})
        protein_class = protein_row.get(
            "Qerythrinae_similarity_review_class", "no_predicted_representative_protein"
        )
        nucleotide = best[gene]
        nucleotide_class_name = nucleotide["nucleotide_provenance_review_class"]
        output: dict[str, object] = dict(source)
        output.update(
            {
                "protein_similarity_available": "TRUE" if protein_row else "FALSE",
                "protein_representative_id": protein_row.get(
                    "representative_protein_id", ""
                ),
                "protein_Qerythrinae_subject_id": protein_row.get("subject_id", ""),
                "protein_percent_identity": protein_row.get("percent_identity", ""),
                "protein_query_span_coverage_percent": protein_row.get(
                    "query_span_coverage_percent", ""
                ),
                "protein_subject_span_coverage_percent": protein_row.get(
                    "subject_span_coverage_percent", ""
                ),
                "protein_evalue": protein_row.get("evalue", ""),
                "protein_Qerythrinae_similarity_review_class": protein_class,
                "Trinity_transcripts_extracted": len(gene_transcripts[gene]),
                "nucleotide_best_query_transcript_id": nucleotide[
                    "best_query_transcript_id"
                ],
                "nucleotide_target_CDS_id": nucleotide["target_CDS_id"],
                "nucleotide_percent_identity": nucleotide["percent_identity"],
                "nucleotide_query_span_coverage_percent": nucleotide[
                    "query_span_coverage_percent"
                ],
                "nucleotide_target_span_coverage_percent": nucleotide[
                    "target_span_coverage_percent"
                ],
                "nucleotide_alignment_block_length": nucleotide[
                    "alignment_block_length"
                ],
                "nucleotide_mapping_quality": nucleotide["mapping_quality"],
                "nucleotide_provenance_review_class": nucleotide_class_name,
                "nucleotide_high_identity_Qerythrinae_CDS_review": nucleotide[
                    "high_identity_Qerythrinae_CDS_review"
                ],
                "overall_Qerythrinae_sequence_review_class": overall_review_class(
                    protein_class, nucleotide_class_name
                ),
                "confirmed_mask_overlap": "FALSE",
            }
        )
        audited_rows.append(output)

    audited_index = {str(row["gene_id"]): row for row in audited_rows}
    review_rows_out = [audited_index[row["gene_id"]] for row in review_rows]
    write_tsv(outputs["top50_audit"], top_fields + appended_fields, audited_rows)
    write_tsv(outputs["review15_audit"], top_fields + appended_fields, review_rows_out)

    summary_rows: list[dict[str, object]] = []
    for set_name, rows in [
        ("host_focused_GP_top50", audited_rows),
        ("host_focused_GP_poorly_annotated_top15", review_rows_out),
        (
            "host_focused_GP_top50_without_predicted_protein",
            [row for row in audited_rows if row["protein_prediction_available"] == "FALSE"],
        ),
    ]:
        nucleotide_counts = Counter(
            str(row["nucleotide_provenance_review_class"]) for row in rows
        )
        overall_counts = Counter(
            str(row["overall_Qerythrinae_sequence_review_class"]) for row in rows
        )
        for label in NUCLEOTIDE_CLASSES:
            summary_rows.append(
                {
                    "candidate_set": set_name,
                    "summary_type": "nucleotide_provenance_review_class",
                    "review_class": label,
                    "candidate_genes": nucleotide_counts[label],
                    "percent_of_candidate_set": format(
                        100.0 * nucleotide_counts[label] / len(rows), ".6g"
                    ),
                }
            )
        for label in [
            "additional_Qerythrinae_like_sequence_review",
            "protein_near_identity_without_strong_nucleotide_confirmation_review",
            "cross_species_homology_review",
            "no_strong_Qerythrinae_similarity_in_current_screen",
        ]:
            summary_rows.append(
                {
                    "candidate_set": set_name,
                    "summary_type": "overall_Qerythrinae_sequence_review_class",
                    "review_class": label,
                    "candidate_genes": overall_counts[label],
                    "percent_of_candidate_set": format(
                        100.0 * overall_counts[label] / len(rows), ".6g"
                    ),
                }
            )
    write_tsv(
        outputs["summary"],
        [
            "candidate_set",
            "summary_type",
            "review_class",
            "candidate_genes",
            "percent_of_candidate_set",
        ],
        summary_rows,
    )

    parameters = [
        {"parameter": "minimap2_version", "value": minimap_version},
        {"parameter": "preset", "value": "asm5"},
        {"parameter": "secondary_alignments", "value": "disabled"},
        {"parameter": "threads", "value": str(THREADS)},
        {
            "parameter": "protein_near_identity_reused_definition",
            "value": "identity>=90%; query and subject span coverage>=80%",
        },
        {
            "parameter": "protein_substantial_homology_reused_definition",
            "value": "identity>=40%; query and subject coverage>=60%; evalue<=1e-10",
        },
        {
            "parameter": "nucleotide_high_identity_reused_definition",
            "value": "identity>=95%; query span coverage>=60%; alignment block length>=300 nt",
        },
        {
            "parameter": "nucleotide_substantial_similarity_reused_definition",
            "value": "identity>=90%; query span coverage>=40%; alignment block length>=200 nt",
        },
        {
            "parameter": "interpretive_boundary",
            "value": (
                "Similarity to the Q. erythrinae reference is a sequence-review flag, "
                "not a taxonomic assignment or proof of physical RNA source"
            ),
        },
    ]
    for label, digest in sorted(hashes.items()):
        parameters.append({"parameter": f"SHA256_{label}", "value": digest})
    write_tsv(outputs["parameters"], ["parameter", "value"], parameters)

    for label, path in inputs.items():
        if sha256(path) != hashes[label]:
            raise AuditError(f"Input changed during analysis: {label}")

    nucleotide_counts = Counter(
        str(row["nucleotide_provenance_review_class"]) for row in audited_rows
    )
    review_nucleotide_counts = Counter(
        str(row["nucleotide_provenance_review_class"]) for row in review_rows_out
    )
    overall_counts = Counter(
        str(row["overall_Qerythrinae_sequence_review_class"]) for row in audited_rows
    )
    protein_counts = Counter(
        str(row["protein_Qerythrinae_similarity_review_class"])
        for row in audited_rows
    )

    report = [
        "===== HOST-FOCUSED PLANT TOP-50 Q. ERYTHRINAE PROVENANCE AUDIT =====",
        f"Discovery root: {root}",
        "",
        "Scope:",
        "  Reused the validated genome-wide plant-protein versus Q. erythrinae protein results.",
        "  Ran minimap2 only for transcripts belonging to the frozen host-focused GP top 50.",
        "  No DESeq2 fitting, filtering, reranking or taxonomic assignment was performed.",
        "",
        f"Frozen host-focused top-50 genes: {len(audited_rows)}",
        f"Poorly annotated/uncharacterized top-15 genes: {len(review_rows_out)}",
        f"Top-50 genes with predicted proteins: {len(protein_bearing)}",
        f"Top-50 genes without predicted proteins: {len(no_protein)}",
        f"Selected Trinity transcripts extracted: {sum(len(value) for value in gene_transcripts.values())}",
        "Confirmed-mask overlap: 0",
        "",
        "Protein-level Q. erythrinae similarity classes in the top 50:",
    ]
    protein_labels = [
        "insect_like_near_identity_review",
        "insect_like_substantial_homology_review",
        "partial_homology_to_Qerythrinae_reference",
        "weak_or_short_match_to_Qerythrinae_reference",
        "no_match_to_Qerythrinae_reference_at_evalue_1e-5",
        "no_predicted_representative_protein",
    ]
    report.extend(f"  {label}: {protein_counts[label]}" for label in protein_labels)
    report.extend(["", "Nucleotide-level Q. erythrinae CDS classes in the top 50:"])
    report.extend(
        f"  {label}: {nucleotide_counts[label]}" for label in NUCLEOTIDE_CLASSES
    )
    report.extend(
        [
            "",
            "Nucleotide-level classes in the poorly annotated/uncharacterized top 15:",
        ]
    )
    report.extend(
        f"  {label}: {review_nucleotide_counts[label]}"
        for label in NUCLEOTIDE_CLASSES
    )
    report.extend(["", "Combined sequence-review classes in the top 50:"])
    for label in [
        "additional_Qerythrinae_like_sequence_review",
        "protein_near_identity_without_strong_nucleotide_confirmation_review",
        "cross_species_homology_review",
        "no_strong_Qerythrinae_similarity_in_current_screen",
    ]:
        report.append(f"  {label}: {overall_counts[label]}")
    report.extend(
        [
            "",
            "Interpretive boundary:",
            "  Strong Q. erythrinae similarity is a review flag, not proof of the physical RNA source.",
            "  Lack of a strong Q. erythrinae match is not proof of plant origin.",
            "  Ordinary conserved homology does not establish contamination, convergence or functional equivalence.",
            "  Broader taxonomic evidence is required before uncharacterized candidates are described as plant-derived or novel.",
            "  The frozen ranking and host-focused DESeq2 results remain unchanged.",
            "",
            "Frozen-input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "",
        ]
    )
    for path in outputs.values():
        report.append(f"Output: {path}")
    report.append("HOST-FOCUSED PLANT TOP-50 Q. ERYTHRINAE PROVENANCE AUDIT: PASS")

    outputs["audit"].parent.mkdir(parents=True, exist_ok=True)
    temporary_audit = outputs["audit"].with_name(outputs["audit"].name + ".tmp")
    temporary_audit.write_text("\n".join(report) + "\n", encoding="utf-8")
    os.replace(temporary_audit, outputs["audit"])
    print("\n".join(report))


if __name__ == "__main__":
    try:
        main()
    except (AuditError, subprocess.CalledProcessError) as error:
        print(
            "===== HOST-FOCUSED PLANT TOP-50 Q. ERYTHRINAE PROVENANCE AUDIT =====",
            file=sys.stderr,
        )
        if isinstance(error, subprocess.CalledProcessError):
            print(f"ERROR: command failed: {' '.join(error.cmd)}", file=sys.stderr)
            if error.stderr:
                print(error.stderr.strip(), file=sys.stderr)
        else:
            print(f"ERROR: {error}", file=sys.stderr)
        print(
            "HOST-FOCUSED PLANT TOP-50 Q. ERYTHRINAE PROVENANCE AUDIT: CHECK REQUIRED",
            file=sys.stderr,
        )
        raise SystemExit(1)
