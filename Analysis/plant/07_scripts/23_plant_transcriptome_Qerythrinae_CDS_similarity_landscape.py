#!/usr/bin/env python3
"""Audit the complete Trinity transcriptome against Q. erythrinae CDS.

This read-only provenance analysis applies the nucleotide thresholds already
used in the discovery project to every Trinity transcript.  It aggregates
transcript mappings to genes by prioritising evidence class before coverage,
compares the result with the current 3,948-gene multilevel mask, and writes a
proposed expanded review mask.  It does not alter or rerun DESeq2.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


MINIMAP = ["micromamba", "run", "-n", "qe_provenance", "minimap2"]
THREADS = 8
EXPECTED_TRINITY_GENES = 118_082
EXPECTED_PROTEIN_UNIVERSE = 41_330
EXPECTED_CURRENT_MASK = 3_948
EXPECTED_TOP50 = 50

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
    "gene_id",
    "percent_identity",
    "query_span_coverage_percent",
    "target_span_coverage_percent",
    "nucleotide_provenance_review_class",
    "strong_nucleotide_Qerythrinae_similarity_review",
]
CLASS_PRIORITY = {
    "high_identity_Qerythrinae_CDS_review": 4,
    "substantial_Qerythrinae_CDS_similarity_review": 3,
    "partial_Qerythrinae_CDS_similarity": 2,
    "weak_or_short_CDS_match": 1,
    "no_CDS_match_reported_by_minimap2": 0,
}
STRONG_CLASSES = {
    "high_identity_Qerythrinae_CDS_review",
    "substantial_Qerythrinae_CDS_similarity_review",
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
        rows = [
            {
                key.strip().lstrip("\ufeff"): (value or "").strip()
                for key, value in row.items()
            }
            for row in reader
        ]
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
        value = row.get(key, "")
        if not value:
            raise AuditError(f"Blank {key} in {label}")
        if value in result:
            raise AuditError(f"Duplicate {key} in {label}: {value}")
        result[value] = row
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


def write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def gene_from_transcript(identifier: str) -> str:
    match = TRINITY_GENE.fullmatch(identifier)
    if match is None:
        raise AuditError(f"Unrecognised Trinity transcript identifier: {identifier}")
    return match.group(1)


def scan_trinity_headers(path: Path) -> tuple[set[str], Counter[str]]:
    genes: set[str] = set()
    transcript_counts: Counter[str] = Counter()
    seen_transcripts: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.startswith(">"):
                continue
            transcript = raw[1:].split(maxsplit=1)[0]
            if transcript in seen_transcripts:
                raise AuditError(f"Duplicate Trinity transcript identifier: {transcript}")
            seen_transcripts.add(transcript)
            gene = gene_from_transcript(transcript)
            genes.add(gene)
            transcript_counts[gene] += 1
    if not seen_transcripts:
        raise AuditError(f"No FASTA records found: {path}")
    return genes, transcript_counts


def run(command: list[str], stdout_path: Path | None = None) -> str:
    if stdout_path is None:
        result = subprocess.run(
            command,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return result.stdout.strip()
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
) -> str:
    if identity >= 95.0 and query_coverage >= 60.0 and alignment_length >= 300:
        return "high_identity_Qerythrinae_CDS_review"
    if identity >= 90.0 and query_coverage >= 40.0 and alignment_length >= 200:
        return "substantial_Qerythrinae_CDS_similarity_review"
    if identity >= 70.0 and alignment_length >= 150:
        return "partial_Qerythrinae_CDS_similarity"
    return "weak_or_short_CDS_match"


def parse_paf(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, values in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(values) < 12:
                raise AuditError(f"Malformed PAF row at {path}:{line_number}")
            row: dict[str, object] = dict(zip(PAF_FIELDS, values[:12]))
            query_length = int(str(row["query_length_nt"]))
            target_length = int(str(row["target_length_nt"]))
            query_span = int(str(row["query_end"])) - int(str(row["query_start"]))
            target_span = int(str(row["target_end"])) - int(str(row["target_start"]))
            matches = int(str(row["matching_bases"]))
            alignment_length = int(str(row["alignment_block_length"]))
            if min(query_length, target_length, query_span, target_span, alignment_length) <= 0:
                raise AuditError(f"Invalid PAF lengths at {path}:{line_number}")
            identity = 100.0 * matches / alignment_length
            query_coverage = 100.0 * query_span / query_length
            target_coverage = 100.0 * target_span / target_length
            label = nucleotide_class(identity, query_coverage, alignment_length)
            transcript = str(row["query_transcript_id"])
            row.update(
                {
                    "gene_id": gene_from_transcript(transcript),
                    "percent_identity": format(identity, ".10g"),
                    "query_span_coverage_percent": format(query_coverage, ".10g"),
                    "target_span_coverage_percent": format(target_coverage, ".10g"),
                    "nucleotide_provenance_review_class": label,
                    "strong_nucleotide_Qerythrinae_similarity_review": (
                        "TRUE" if label in STRONG_CLASSES else "FALSE"
                    ),
                }
            )
            rows.append(row)
    return rows


def best_mapping_per_gene(
    mappings: list[dict[str, object]], all_genes: set[str]
) -> dict[str, dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in mappings:
        gene = str(row["gene_id"])
        if gene not in all_genes:
            raise AuditError(f"Mapped gene absent from Trinity universe: {gene}")
        grouped[gene].append(row)

    result: dict[str, dict[str, object]] = {}
    for gene in sorted(all_genes):
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
                "strong_nucleotide_Qerythrinae_similarity_review": "FALSE",
            }
            continue
        best = sorted(
            candidates,
            key=lambda row: (
                -CLASS_PRIORITY[str(row["nucleotide_provenance_review_class"])],
                -float(str(row["query_span_coverage_percent"])),
                -float(str(row["percent_identity"])),
                -int(str(row["alignment_block_length"])),
                -int(str(row["mapping_quality"])),
                str(row["query_transcript_id"]),
                str(row["target_CDS_id"]),
            ),
        )[0]
        result[gene] = {
            "best_query_transcript_id": best["query_transcript_id"],
            "target_CDS_id": best["target_CDS_id"],
            "percent_identity": best["percent_identity"],
            "query_span_coverage_percent": best[
                "query_span_coverage_percent"
            ],
            "target_span_coverage_percent": best[
                "target_span_coverage_percent"
            ],
            "alignment_block_length": best["alignment_block_length"],
            "mapping_quality": best["mapping_quality"],
            "nucleotide_provenance_review_class": best[
                "nucleotide_provenance_review_class"
            ],
            "strong_nucleotide_Qerythrinae_similarity_review": best[
                "strong_nucleotide_Qerythrinae_similarity_review"
            ],
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--discovery-root",
        type=Path,
        default=Path.home() / "projects/qe_cross_species_discovery",
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
    parser.add_argument("--threads", type=int, default=THREADS)
    args = parser.parse_args()
    if args.threads < 1:
        raise AuditError("--threads must be at least 1")

    root = args.discovery_root.resolve()
    integration = root / "04_cross_species_integration/results"
    host_results = root / "05_host_focused_DESeq2_sensitivity/results"
    inputs = {
        "trinity_assembly": args.trinity_assembly.resolve(),
        "insect_CDS": args.insect_cds.resolve(),
        "protein_similarity": integration
        / "plant_representative_proteome_vs_Qerythrinae_best_hits.tsv",
        "current_confirmed_mask": integration
        / "plant_confirmed_multilevel_insect_like_sequence_mask.tsv",
        "current_host_focused_DESeq2": host_results
        / "plant_gene_DESeq2_GP_vs_UP_host_focused_all.tsv",
        "corrected_top50_audit": host_results
        / "plant_GP_host_focused_top50_sequence_provenance_audit.tsv",
    }
    outputs = {
        "all_mappings": integration
        / "plant_transcriptome_vs_Qerythrinae_CDS_all_mappings.tsv",
        "gene_best": integration
        / "plant_transcriptome_vs_Qerythrinae_CDS_best_by_gene.tsv",
        "expanded_mask": integration
        / "plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv",
        "status_summary": integration
        / "plant_expanded_nucleotide_review_summary_by_current_host_DESeq2_status.tsv",
        "top50_recheck": host_results
        / "plant_GP_host_focused_top50_transcriptomewide_provenance_recheck.tsv",
        "parameters": root
        / "00_project_control/23_plant_transcriptome_Qerythrinae_CDS_similarity_parameters.tsv",
        "audit": root
        / "08_logs/23_plant_transcriptome_Qerythrinae_CDS_similarity_landscape_audit.txt",
    }

    for label, path in inputs.items():
        if not path.is_file():
            raise AuditError(f"Missing input {label}: {path}")
    if outputs["audit"].is_file() and outputs["audit"].read_text(
        encoding="utf-8", errors="replace"
    ).rstrip().endswith(
        "PLANT TRANSCRIPTOME–Q. ERYTHRINAE CDS SIMILARITY LANDSCAPE: PASS"
    ):
        raise AuditError(
            f"A completed validated audit already exists: {outputs['audit']}. "
            "No outputs were changed."
        )

    hashes = {label: sha256(path) for label, path in inputs.items()}
    if shutil.which("micromamba") is None:
        raise AuditError("micromamba is not available on PATH")
    minimap_version = run([*MINIMAP, "--version"])

    trinity_genes, transcript_counts = scan_trinity_headers(
        inputs["trinity_assembly"]
    )
    if len(trinity_genes) != EXPECTED_TRINITY_GENES:
        raise AuditError(
            f"Trinity gene universe: {len(trinity_genes)}; "
            f"expected {EXPECTED_TRINITY_GENES}"
        )

    protein_fields, protein_rows = read_tsv(inputs["protein_similarity"])
    mask_fields, mask_rows = read_tsv(inputs["current_confirmed_mask"])
    host_fields, host_rows = read_tsv(inputs["current_host_focused_DESeq2"])
    top_fields, top_rows = read_tsv(inputs["corrected_top50_audit"])
    require_columns(
        inputs["protein_similarity"],
        protein_fields,
        {"gene_id", "Qerythrinae_similarity_review_class"},
    )
    require_columns(
        inputs["current_confirmed_mask"],
        mask_fields,
        {"gene_id", "host_reanalysis_mask_status"},
    )
    require_columns(
        inputs["current_host_focused_DESeq2"],
        host_fields,
        {"gene_id", "DESeq2_status_host_focused"},
    )
    require_columns(
        inputs["corrected_top50_audit"],
        top_fields,
        {
            "gene_id",
            "discovery_rank",
            "nucleotide_best_query_transcript_id",
            "nucleotide_target_CDS_id",
            "nucleotide_provenance_review_class",
        },
    )
    protein = index_unique(protein_rows, "gene_id", "protein similarity table")
    current_mask_table = index_unique(mask_rows, "gene_id", "current mask")
    host = index_unique(host_rows, "gene_id", "current host-focused DESeq2")
    top50 = index_unique(top_rows, "gene_id", "corrected host-focused top 50")
    if len(protein) != EXPECTED_PROTEIN_UNIVERSE:
        raise AuditError(
            f"Protein similarity rows: {len(protein)}; expected {EXPECTED_PROTEIN_UNIVERSE}"
        )
    if len(current_mask_table) != EXPECTED_PROTEIN_UNIVERSE:
        raise AuditError(
            f"Current mask-universe rows: {len(current_mask_table)}; "
            f"expected {EXPECTED_PROTEIN_UNIVERSE}"
        )
    if len(top50) != EXPECTED_TOP50:
        raise AuditError(f"Corrected top-50 rows: {len(top50)}; expected {EXPECTED_TOP50}")
    if not set(protein).issubset(trinity_genes):
        raise AuditError("Protein-similarity genes are not a subset of Trinity genes")
    if not set(host).issubset(trinity_genes):
        raise AuditError("Host-focused DESeq2 genes are not a subset of Trinity genes")

    current_mask = {
        gene
        for gene, row in current_mask_table.items()
        if row["host_reanalysis_mask_status"]
        == "review_exclude_confirmed_multilevel_insect_like_sequence"
    }
    if len(current_mask) != EXPECTED_CURRENT_MASK:
        raise AuditError(
            f"Current confirmed mask: {len(current_mask)}; "
            f"expected {EXPECTED_CURRENT_MASK}"
        )

    with tempfile.TemporaryDirectory(prefix="plant_transcriptome_Qe_CDS_") as temp_name:
        paf = Path(temp_name) / "Trinity_vs_Qerythrinae_CDS.paf"
        run(
            [
                *MINIMAP,
                "-x",
                "asm5",
                "--secondary=no",
                "-t",
                str(args.threads),
                str(inputs["insect_CDS"]),
                str(inputs["trinity_assembly"]),
            ],
            stdout_path=paf,
        )
        mapping_rows = parse_paf(paf)

    write_tsv(
        outputs["all_mappings"],
        [*PAF_FIELDS, *DERIVED_FIELDS],
        mapping_rows,
    )
    best = best_mapping_per_gene(mapping_rows, trinity_genes)
    strong_genes = {
        gene
        for gene, row in best.items()
        if row["strong_nucleotide_Qerythrinae_similarity_review"] == "TRUE"
    }
    if not current_mask.issubset(strong_genes):
        missing = sorted(current_mask.difference(strong_genes))
        raise AuditError(
            "Current confirmed-mask genes not recovered by transcriptome-wide "
            f"nucleotide audit: {len(missing)}; {', '.join(missing[:20])}"
        )

    gene_best_fields = [
        "gene_id",
        "Trinity_transcripts",
        "best_query_transcript_id",
        "target_CDS_id",
        "percent_identity",
        "query_span_coverage_percent",
        "target_span_coverage_percent",
        "alignment_block_length",
        "mapping_quality",
        "nucleotide_provenance_review_class",
        "strong_nucleotide_Qerythrinae_similarity_review",
        "representative_protein_available",
        "protein_Qerythrinae_similarity_review_class",
        "current_confirmed_mask_overlap",
        "current_host_focused_DESeq2_status",
    ]
    gene_best_rows: list[dict[str, object]] = []
    expanded_mask_rows: list[dict[str, object]] = []
    for gene in sorted(trinity_genes):
        evidence = best[gene]
        protein_row = protein.get(gene, {})
        if gene in host:
            current_status = host[gene]["DESeq2_status_host_focused"]
        elif gene in current_mask:
            current_status = "excluded_by_current_3948_gene_mask"
        else:
            current_status = "excluded_by_current_host_focused_low_count_filter"
        gene_best_rows.append(
            {
                "gene_id": gene,
                "Trinity_transcripts": transcript_counts[gene],
                **evidence,
                "representative_protein_available": (
                    "TRUE" if protein_row else "FALSE"
                ),
                "protein_Qerythrinae_similarity_review_class": protein_row.get(
                    "Qerythrinae_similarity_review_class", ""
                ),
                "current_confirmed_mask_overlap": (
                    "TRUE" if gene in current_mask else "FALSE"
                ),
                "current_host_focused_DESeq2_status": current_status,
            }
        )
        if gene in current_mask:
            mask_status = "review_exclude_existing_confirmed_multilevel_sequence"
            rationale = (
                "Already present in the validated 3,948-gene multilevel mask; "
                "strong nucleotide similarity recovered"
            )
        elif gene in strong_genes:
            mask_status = "review_exclude_additional_strong_nucleotide_similarity"
            rationale = (
                "High/substantial nucleotide similarity to Q. erythrinae CDS "
                "under the frozen thresholds"
            )
        else:
            mask_status = "retain_no_strong_nucleotide_Qerythrinae_similarity"
            rationale = (
                "No high/substantial Q. erythrinae CDS similarity under the "
                "frozen thresholds"
            )
        expanded_mask_rows.append(
            {
                "gene_id": gene,
                "nucleotide_provenance_review_class": evidence[
                    "nucleotide_provenance_review_class"
                ],
                "current_confirmed_mask_overlap": (
                    "TRUE" if gene in current_mask else "FALSE"
                ),
                "proposed_expanded_mask_status": mask_status,
                "proposed_expanded_mask_rationale": rationale,
                "current_host_focused_DESeq2_status": current_status,
            }
        )
    write_tsv(outputs["gene_best"], gene_best_fields, gene_best_rows)
    expanded_fields = [
        "gene_id",
        "nucleotide_provenance_review_class",
        "current_confirmed_mask_overlap",
        "proposed_expanded_mask_status",
        "proposed_expanded_mask_rationale",
        "current_host_focused_DESeq2_status",
    ]
    write_tsv(outputs["expanded_mask"], expanded_fields, expanded_mask_rows)

    top50_recheck_rows: list[dict[str, object]] = []
    for source in sorted(top_rows, key=lambda row: int(row["discovery_rank"])):
        gene = source["gene_id"]
        evidence = best[gene]
        for old_field, new_field in (
            (
                "nucleotide_best_query_transcript_id",
                "best_query_transcript_id",
            ),
            ("nucleotide_target_CDS_id", "target_CDS_id"),
            (
                "nucleotide_provenance_review_class",
                "nucleotide_provenance_review_class",
            ),
        ):
            if source[old_field] != str(evidence[new_field]):
                raise AuditError(
                    f"Top-50 transcriptome-wide recheck mismatch for {gene} "
                    f"in {old_field}: top50={source[old_field]!r}, "
                    f"full={evidence[new_field]!r}"
                )
        top50_recheck_rows.append(
            {
                "discovery_rank": source["discovery_rank"],
                "gene_id": gene,
                "nucleotide_best_query_transcript_id": evidence[
                    "best_query_transcript_id"
                ],
                "nucleotide_target_CDS_id": evidence["target_CDS_id"],
                "nucleotide_percent_identity": evidence["percent_identity"],
                "nucleotide_query_span_coverage_percent": evidence[
                    "query_span_coverage_percent"
                ],
                "nucleotide_alignment_block_length": evidence[
                    "alignment_block_length"
                ],
                "nucleotide_provenance_review_class": evidence[
                    "nucleotide_provenance_review_class"
                ],
                "proposed_expanded_mask_overlap": (
                    "TRUE" if gene in strong_genes else "FALSE"
                ),
                "current_confirmed_mask_overlap": (
                    "TRUE" if gene in current_mask else "FALSE"
                ),
            }
        )
    top50_fields = list(top50_recheck_rows[0])
    write_tsv(outputs["top50_recheck"], top50_fields, top50_recheck_rows)

    summary_counter: dict[str, Counter[str]] = defaultdict(Counter)
    summary_totals: Counter[str] = Counter()
    for row in expanded_mask_rows:
        status = str(row["current_host_focused_DESeq2_status"])
        label = str(row["proposed_expanded_mask_status"])
        summary_counter[status][label] += 1
        summary_totals[status] += 1
    summary_rows: list[dict[str, object]] = []
    for status in sorted(summary_totals):
        for label in (
            "review_exclude_existing_confirmed_multilevel_sequence",
            "review_exclude_additional_strong_nucleotide_similarity",
            "retain_no_strong_nucleotide_Qerythrinae_similarity",
        ):
            count = summary_counter[status][label]
            summary_rows.append(
                {
                    "current_host_focused_DESeq2_status": status,
                    "proposed_expanded_mask_status": label,
                    "genes": count,
                    "percent_within_current_status": format(
                        100.0 * count / summary_totals[status], ".8g"
                    ),
                }
            )
    write_tsv(
        outputs["status_summary"],
        [
            "current_host_focused_DESeq2_status",
            "proposed_expanded_mask_status",
            "genes",
            "percent_within_current_status",
        ],
        summary_rows,
    )

    class_counts = Counter(
        str(row["nucleotide_provenance_review_class"])
        for row in best.values()
    )
    additional_strong = strong_genes.difference(current_mask)
    top50_strong = sum(
        best[gene]["strong_nucleotide_Qerythrinae_similarity_review"] == "TRUE"
        for gene in top50
    )
    params: list[dict[str, object]] = [
        {"parameter": "minimap2_version", "value": minimap_version},
        {"parameter": "minimap2_preset", "value": "asm5"},
        {"parameter": "secondary_alignments", "value": "disabled"},
        {"parameter": "threads", "value": args.threads},
        {
            "parameter": "gene_level_best_alignment_rule",
            "value": (
                "evidence class first (high > substantial > partial > weak), "
                "then query coverage, identity, alignment length and mapping quality"
            ),
        },
        {
            "parameter": "high_identity_threshold",
            "value": "identity >=95%; query coverage >=60%; alignment length >=300 nt",
        },
        {
            "parameter": "substantial_similarity_threshold",
            "value": "identity >=90%; query coverage >=40%; alignment length >=200 nt",
        },
        {
            "parameter": "partial_similarity_threshold",
            "value": "identity >=70%; alignment length >=150 nt",
        },
        {
            "parameter": "interpretive_boundary",
            "value": (
                "Proposed mask is a sequence-review resource, not a taxonomic "
                "assignment or proof of the physical RNA source"
            ),
        },
    ]
    params.extend(
        {"parameter": f"SHA256_{label}", "value": value}
        for label, value in sorted(hashes.items())
    )
    write_tsv(outputs["parameters"], ["parameter", "value"], params)

    for label, path in inputs.items():
        if sha256(path) != hashes[label]:
            raise AuditError(f"Input changed during analysis: {label}")

    report = [
        "===== PLANT TRANSCRIPTOME–Q. ERYTHRINAE CDS SIMILARITY LANDSCAPE =====",
        f"Discovery root: {root}",
        "",
        "Scope:",
        "  Complete Trinity transcriptome screened against Q. erythrinae CDS.",
        "  Gene-level evidence class was prioritized before coverage and identity.",
        "  The output is a proposed expanded review mask; no DESeq2 result was altered.",
        "",
        f"Trinity genes audited: {len(trinity_genes)}",
        f"Trinity transcripts audited: {sum(transcript_counts.values())}",
        f"Reported minimap2 mapping rows: {len(mapping_rows)}",
        "",
        "Gene-level nucleotide provenance classes:",
    ]
    for label in (
        "high_identity_Qerythrinae_CDS_review",
        "substantial_Qerythrinae_CDS_similarity_review",
        "partial_Qerythrinae_CDS_similarity",
        "weak_or_short_CDS_match",
        "no_CDS_match_reported_by_minimap2",
    ):
        report.append(f"  {label}: {class_counts[label]}")
    report.extend(
        [
            "",
            f"Current confirmed multilevel mask: {len(current_mask)}",
            f"All strong nucleotide-review genes: {len(strong_genes)}",
            f"Additional strong nucleotide-review genes: {len(additional_strong)}",
            f"Proposed expanded mask total: {len(strong_genes)}",
            f"Frozen host-focused top-50 strong-review genes: {top50_strong} / 50",
            "Top-50 transcriptome-wide identity recheck: PASS",
            "Current-mask recovery within expanded mask: PASS",
            "",
            "Interpretive boundary:",
            "  Strong sequence similarity is a review flag, not proof of physical RNA source.",
            "  The proposed expanded mask has not been applied to DESeq2 by this script.",
            "  Ordinary partial homology is retained and is not treated as an exclusion criterion.",
            "  Final mask adoption should be decided only after reviewing these counts and statuses.",
            "",
            "Input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "",
        ]
    )
    report.extend(f"Output: {path}" for path in outputs.values())
    report.append(
        "PLANT TRANSCRIPTOME–Q. ERYTHRINAE CDS SIMILARITY LANDSCAPE: PASS"
    )
    write_text_atomic(outputs["audit"], "\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    try:
        main()
    except (AuditError, subprocess.CalledProcessError) as error:
        print(
            "===== PLANT TRANSCRIPTOME–Q. ERYTHRINAE CDS SIMILARITY LANDSCAPE =====",
            file=sys.stderr,
        )
        if isinstance(error, subprocess.CalledProcessError):
            print(f"ERROR: command failed: {' '.join(error.cmd)}", file=sys.stderr)
            if error.stderr:
                print(error.stderr.strip(), file=sys.stderr)
        else:
            print(f"ERROR: {error}", file=sys.stderr)
        print(
            "PLANT TRANSCRIPTOME–Q. ERYTHRINAE CDS SIMILARITY LANDSCAPE: CHECK REQUIRED",
            file=sys.stderr,
        )
        raise SystemExit(1)
