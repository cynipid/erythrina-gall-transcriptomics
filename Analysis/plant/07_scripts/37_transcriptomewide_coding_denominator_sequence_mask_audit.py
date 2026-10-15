#!/usr/bin/env python3
"""Audit coding-aware nucleotide coverage across the plant transcriptome.

The completed minimap2 gene-best table measured query coverage against whole
Trinity transcripts. This audit asks how many currently retained genes meet
strong sequence-review criteria when coverage is instead supported by the
Q. erythrinae CDS target or by the predicted protein coding length. It changes
no mask and performs no alignment or DESeq2 analysis.
"""

from __future__ import annotations

import csv
import hashlib
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path.home() / "projects" / "qe_cross_species_discovery"
INTEGRATION = ROOT / "04_cross_species_integration" / "results"
RESULTS = ROOT / "05_host_focused_DESeq2_sensitivity" / "results"
CONTROL = ROOT / "00_project_control"
LOGS = ROOT / "05_host_focused_DESeq2_sensitivity" / "logs"

GENE_BEST = INTEGRATION / "plant_transcriptome_vs_Qerythrinae_CDS_best_by_gene.tsv"
PROTEIN = INTEGRATION / "plant_representative_proteome_vs_Qerythrinae_best_hits.tsv"
CURRENT_MASK = INTEGRATION / "plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv"
EXPANDED_DESEQ = RESULTS / "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_eggnog_joined.tsv"
TOP50 = RESULTS / "plant_GP_host_focused_expanded_mask_top50_final_quality_provenance_audit.tsv"

CANDIDATES_OUTPUT = INTEGRATION / "plant_coding_denominator_additional_sequence_review_candidates.tsv"
SUMMARY_OUTPUT = INTEGRATION / "plant_coding_denominator_sequence_review_summary.tsv"
TOP50_OUTPUT = RESULTS / "plant_GP_host_focused_expanded_mask_top50_coding_denominator_recheck.tsv"
PARAMETERS_OUTPUT = CONTROL / "37_transcriptomewide_coding_denominator_sequence_mask_audit_parameters.tsv"
AUDIT_OUTPUT = LOGS / "37_transcriptomewide_coding_denominator_sequence_mask_audit.txt"

CURRENT_MASK_EXCLUDE = {
    "review_exclude_existing_confirmed_multilevel_sequence",
    "review_exclude_additional_strong_nucleotide_similarity",
}
SUPPORTED_PROTEIN_CLASSES = {
    "insect_like_near_identity_review",
    "insect_like_substantial_homology_review",
}

MIN_NUCLEOTIDE_IDENTITY = 90.0
MIN_ALIGNMENT_LENGTH = 200
MIN_TARGET_CDS_COVERAGE = 80.0
MIN_PREDICTED_CODING_COVERAGE = 80.0


class AuditError(RuntimeError):
    pass


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip().lstrip("\ufeff")


def number(value: object) -> float | None:
    text = clean(value)
    if not text or text.lower() in {"na", "nan", "none", "null"}:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


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
        return list(reader.fieldnames), [dict(row) for row in reader]


def require_columns(path: Path, fields: list[str], required: set[str]) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise AuditError(f"Missing columns in {path}: {', '.join(missing)}")


def unique_index(
    path: Path, fields: list[str], rows: list[dict[str, str]], key: str = "gene_id"
) -> dict[str, dict[str, str]]:
    require_columns(path, fields, {key})
    output: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = clean(row.get(key))
        if not identifier:
            raise AuditError(f"Blank {key} in {path}")
        if identifier in output:
            raise AuditError(f"Duplicate {key} in {path}: {identifier}")
        output[identifier] = row
    return output


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    required = [GENE_BEST, PROTEIN, CURRENT_MASK, EXPANDED_DESEQ, TOP50]
    hashes_before = {str(path): sha256(path) for path in required}
    best_fields, best_rows = read_tsv(GENE_BEST)
    protein_fields, protein_rows = read_tsv(PROTEIN)
    mask_fields, mask_rows = read_tsv(CURRENT_MASK)
    deseq_fields, deseq_rows = read_tsv(EXPANDED_DESEQ)
    top_fields, top_rows = read_tsv(TOP50)

    require_columns(
        GENE_BEST,
        best_fields,
        {
            "gene_id", "Trinity_transcripts", "best_query_transcript_id", "target_CDS_id",
            "percent_identity", "query_span_coverage_percent", "target_span_coverage_percent",
            "alignment_block_length", "nucleotide_provenance_review_class",
            "strong_nucleotide_Qerythrinae_similarity_review", "representative_protein_available",
            "protein_Qerythrinae_similarity_review_class",
        },
    )
    require_columns(
        PROTEIN,
        protein_fields,
        {
            "gene_id", "representative_protein_id", "representative_protein_length_aa",
            "percent_identity", "alignment_length", "query_length", "subject_length",
            "query_span_coverage_percent", "subject_span_coverage_percent",
            "Qerythrinae_similarity_review_class",
        },
    )
    require_columns(
        CURRENT_MASK,
        mask_fields,
        {"gene_id", "proposed_expanded_mask_status"},
    )
    require_columns(
        EXPANDED_DESEQ,
        deseq_fields,
        {"gene_id", "DESeq2_status_host_focused_expanded_mask"},
    )
    require_columns(TOP50, top_fields, {"gene_id", "discovery_rank"})

    best = unique_index(GENE_BEST, best_fields, best_rows)
    protein = unique_index(PROTEIN, protein_fields, protein_rows)
    mask = unique_index(CURRENT_MASK, mask_fields, mask_rows)
    deseq = unique_index(EXPANDED_DESEQ, deseq_fields, deseq_rows)
    top = unique_index(TOP50, top_fields, top_rows)

    if len(best_rows) != 118_082 or len(mask_rows) != 118_082:
        raise AuditError(
            f"Transcriptome universe differs: gene_best={len(best_rows)}, mask={len(mask_rows)}"
        )
    if len(protein_rows) != 41_330 or len(deseq_rows) != 24_931 or len(top_rows) != 50:
        raise AuditError(
            f"Input counts differ: protein={len(protein_rows)}, DESeq2={len(deseq_rows)}, top50={len(top_rows)}"
        )
    if set(best) != set(mask):
        raise AuditError("Gene-best and expanded-mask universes differ")
    current_excluded = {
        gene for gene, row in mask.items()
        if clean(row.get("proposed_expanded_mask_status")) in CURRENT_MASK_EXCLUDE
    }
    if len(current_excluded) != 16_112:
        raise AuditError(f"Current expanded mask has {len(current_excluded)} genes; expected 16112")

    candidates: list[dict[str, object]] = []
    for gene, nucleotide in best.items():
        if gene in current_excluded:
            continue
        identity = number(nucleotide.get("percent_identity"))
        alignment_length = number(nucleotide.get("alignment_block_length"))
        target_coverage = number(nucleotide.get("target_span_coverage_percent"))
        query_coverage = number(nucleotide.get("query_span_coverage_percent"))
        if identity is None or alignment_length is None:
            continue
        common = identity >= MIN_NUCLEOTIDE_IDENTITY and alignment_length >= MIN_ALIGNMENT_LENGTH
        if not common:
            continue

        protein_row = protein.get(gene, {})
        protein_length = number(protein_row.get("representative_protein_length_aa"))
        protein_class = clean(protein_row.get("Qerythrinae_similarity_review_class"))
        coding_length = 3.0 * protein_length if protein_length is not None else None
        coding_coverage = (
            100.0 * alignment_length / coding_length
            if coding_length is not None and coding_length > 0
            else None
        )
        target_supported = target_coverage is not None and target_coverage >= MIN_TARGET_CDS_COVERAGE
        coding_supported = (
            coding_coverage is not None
            and coding_coverage >= MIN_PREDICTED_CODING_COVERAGE
            and protein_class in SUPPORTED_PROTEIN_CLASSES
        )
        if not (target_supported or coding_supported):
            continue

        rules: list[str] = []
        if target_supported:
            rules.append("target_CDS_coverage_ge80pct")
        if coding_supported and protein_class == "insect_like_near_identity_review":
            rules.append("protein_near_identity_plus_predicted_coding_coverage_ge80pct")
        if coding_supported and protein_class == "insect_like_substantial_homology_review":
            rules.append("protein_substantial_homology_plus_predicted_coding_coverage_ge80pct")

        if gene in deseq:
            current_status = clean(deseq[gene].get("DESeq2_status_host_focused_expanded_mask"))
        else:
            current_status = "not_in_expanded_mask_host_focused_tested_universe"
        candidates.append(
            {
                "gene_id": gene,
                "current_expanded_mask_status": mask[gene]["proposed_expanded_mask_status"],
                "current_expanded_mask_DESeq2_status": current_status,
                "current_top50_member": "TRUE" if gene in top else "FALSE",
                "current_top50_rank": top.get(gene, {}).get("discovery_rank", ""),
                "best_query_transcript_id": nucleotide["best_query_transcript_id"],
                "target_CDS_id": nucleotide["target_CDS_id"],
                "nucleotide_percent_identity": identity,
                "nucleotide_alignment_block_length": int(alignment_length),
                "whole_transcript_query_coverage_percent": query_coverage if query_coverage is not None else "",
                "target_CDS_coverage_percent": target_coverage if target_coverage is not None else "",
                "representative_protein_available": "TRUE" if protein_row else "FALSE",
                "representative_protein_id": protein_row.get("representative_protein_id", ""),
                "representative_protein_length_aa": protein_length if protein_length is not None else "",
                "predicted_coding_length_3x_aa": coding_length if coding_length is not None else "",
                "nucleotide_alignment_vs_predicted_coding_length_percent": (
                    coding_coverage if coding_coverage is not None else ""
                ),
                "protein_percent_identity": protein_row.get("percent_identity", ""),
                "protein_query_span_coverage_percent": protein_row.get("query_span_coverage_percent", ""),
                "protein_subject_span_coverage_percent": protein_row.get("subject_span_coverage_percent", ""),
                "protein_Qerythrinae_similarity_review_class": protein_class,
                "target_CDS_coverage_rule": "TRUE" if target_supported else "FALSE",
                "protein_coding_coverage_rule": "TRUE" if coding_supported else "FALSE",
                "coding_denominator_review_rules": ";".join(rules),
                "proposed_review_status": "review_add_coding_denominator_supported_sequence",
            }
        )

    candidates.sort(
        key=lambda row: (
            0 if row["current_top50_member"] == "TRUE" else 1,
            int(row["current_top50_rank"]) if row["current_top50_rank"] else 10**9,
            0 if row["current_expanded_mask_DESeq2_status"] == "GP_up" else 1,
            -float(row["nucleotide_percent_identity"]),
            -int(row["nucleotide_alignment_block_length"]),
            str(row["gene_id"]),
        )
    )
    fields = [
        "gene_id", "current_expanded_mask_status", "current_expanded_mask_DESeq2_status",
        "current_top50_member", "current_top50_rank", "best_query_transcript_id", "target_CDS_id",
        "nucleotide_percent_identity", "nucleotide_alignment_block_length",
        "whole_transcript_query_coverage_percent", "target_CDS_coverage_percent",
        "representative_protein_available", "representative_protein_id",
        "representative_protein_length_aa", "predicted_coding_length_3x_aa",
        "nucleotide_alignment_vs_predicted_coding_length_percent", "protein_percent_identity",
        "protein_query_span_coverage_percent", "protein_subject_span_coverage_percent",
        "protein_Qerythrinae_similarity_review_class", "target_CDS_coverage_rule",
        "protein_coding_coverage_rule", "coding_denominator_review_rules",
        "proposed_review_status",
    ]
    write_tsv(CANDIDATES_OUTPUT, candidates, fields)

    top_recheck = [row for row in candidates if row["current_top50_member"] == "TRUE"]
    top_coding_supported = sum(
        row["protein_coding_coverage_rule"] == "TRUE" for row in top_recheck
    )
    top_target_only = sum(
        row["target_CDS_coverage_rule"] == "TRUE"
        and row["protein_coding_coverage_rule"] == "FALSE"
        for row in top_recheck
    )
    if top_coding_supported < 10:
        raise AuditError(
            "Known top-50 protein/coding-denominator cases were not all recovered: "
            f"protein-coding-supported={top_coding_supported}; expected at least 10"
        )
    write_tsv(TOP50_OUTPUT, top_recheck, fields)

    status_counts = Counter(str(row["current_expanded_mask_DESeq2_status"]) for row in candidates)
    protein_class_counts = Counter(
        str(row["protein_Qerythrinae_similarity_review_class"]) or "no_representative_protein"
        for row in candidates
    )
    rule_counts: Counter[str] = Counter()
    rule_status_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in candidates:
        for rule in str(row["coding_denominator_review_rules"]).split(";"):
            if not rule:
                continue
            rule_counts[rule] += 1
            rule_status_counts[rule][str(row["current_expanded_mask_DESeq2_status"])] += 1

    summary_rows: list[dict[str, object]] = []
    for status, count in sorted(status_counts.items()):
        summary_rows.append(
            {"summary_dimension": "current_expanded_mask_DESeq2_status", "category": status, "candidate_genes": count}
        )
    for protein_class, count in sorted(protein_class_counts.items()):
        summary_rows.append(
            {"summary_dimension": "protein_similarity_class", "category": protein_class, "candidate_genes": count}
        )
    for rule, count in sorted(rule_counts.items()):
        summary_rows.append(
            {"summary_dimension": "coding_denominator_review_rule", "category": rule, "candidate_genes": count}
        )
        for status, status_count in sorted(rule_status_counts[rule].items()):
            summary_rows.append(
                {
                    "summary_dimension": f"rule_by_DESeq2_status:{rule}",
                    "category": status,
                    "candidate_genes": status_count,
                }
            )
    write_tsv(
        SUMMARY_OUTPUT,
        summary_rows,
        ["summary_dimension", "category", "candidate_genes"],
    )

    parameters = [
        {"parameter": "scope", "value": "audit only; no mask or DESeq2 result changed"},
        {"parameter": "minimum_nucleotide_identity_percent", "value": MIN_NUCLEOTIDE_IDENTITY},
        {"parameter": "minimum_alignment_block_length_nt", "value": MIN_ALIGNMENT_LENGTH},
        {"parameter": "minimum_target_CDS_coverage_percent", "value": MIN_TARGET_CDS_COVERAGE},
        {"parameter": "minimum_predicted_coding_coverage_percent", "value": MIN_PREDICTED_CODING_COVERAGE},
        {
            "parameter": "coding_coverage_definition",
            "value": "100 * nucleotide alignment block length / (3 * representative protein length aa)",
        },
        {
            "parameter": "protein_support_for_coding_coverage_rule",
            "value": "protein near-identity or substantial-homology class under frozen DIAMOND definitions",
        },
        {
            "parameter": "rationale",
            "value": "whole-Trinity-transcript coverage can be depressed by untranslated sequence",
        },
    ]
    parameters.extend(
        {"parameter": f"SHA256_{Path(path).name}", "value": digest}
        for path, digest in sorted(hashes_before.items())
    )
    write_tsv(PARAMETERS_OUTPUT, parameters, ["parameter", "value"])

    if any(sha256(Path(path)) != digest for path, digest in hashes_before.items()):
        raise AuditError("One or more frozen inputs changed during output generation")

    lines = [
        "===== TRANSCRIPTOME-WIDE CODING-DENOMINATOR SEQUENCE-MASK AUDIT =====",
        "Scope: denominator sensitivity audit using completed nucleotide/protein alignments.",
        "No alignment, mask application, DESeq2 fitting, candidate exclusion or reranking was performed.",
        "",
        f"Trinity genes reviewed: {len(best_rows)}",
        f"Current expanded mask: {len(current_excluded)}",
        f"Additional coding-denominator-supported review genes: {len(candidates)}",
        f"Current top-50 genes recovered: {len(top_recheck)}",
        f"  Protein/coding-length-supported top-50 genes: {top_coding_supported}",
        f"  Additional target-CDS-only top-50 genes: {top_target_only}",
        "",
        "Additional review genes by current expanded-mask DESeq2 status:",
    ]
    for status, count in sorted(status_counts.items()):
        lines.append(f"  {status}: {count}")
    lines.extend(["", "Additional review genes by protein-similarity class:"])
    for label, count in sorted(protein_class_counts.items()):
        lines.append(f"  {label}: {count}")
    lines.extend(["", "Rule recovery counts (rules can overlap):"])
    for rule, count in sorted(rule_counts.items()):
        lines.append(f"  {rule}: {count}")
    lines.extend(["", "Recovered current top-50 genes:"])
    for row in top_recheck:
        lines.append(
            "  "
            + " | ".join(
                [
                    f"rank {row['current_top50_rank']}",
                    str(row["gene_id"]),
                    f"nt={float(row['nucleotide_percent_identity']):.6g}%/{row['nucleotide_alignment_block_length']}nt",
                    f"whole_tx_qcov={float(row['whole_transcript_query_coverage_percent']):.6g}%",
                    f"target_cov={float(row['target_CDS_coverage_percent']):.6g}%",
                    f"coding_cov={float(row['nucleotide_alignment_vs_predicted_coding_length_percent']):.6g}%",
                    str(row["protein_Qerythrinae_similarity_review_class"]),
                    str(row["coding_denominator_review_rules"]),
                ]
            )
        )
    lines.extend(
        [
            "",
            "Interpretive boundary:",
            "  These genes are proposed for sequence-provenance review, not automatically excluded.",
            "  Coding-length coverage is a denominator sensitivity measure and may exceed 100% because alignment blocks can include gaps or adjacent sequence.",
            "  Target-CDS and coding-length support address long-UTR dilution of whole-transcript query coverage.",
            "  A refined mask and another DESeq2 refit should be considered only after reviewing the genome-wide counts and status concentration.",
            "",
            "Frozen-input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "TRANSCRIPTOME-WIDE CODING-DENOMINATOR SEQUENCE-MASK AUDIT: PASS",
        ]
    )
    AUDIT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT_OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    for path in [CANDIDATES_OUTPUT, SUMMARY_OUTPUT, TOP50_OUTPUT, PARAMETERS_OUTPUT, AUDIT_OUTPUT]:
        print(f"Output: {path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("===== TRANSCRIPTOME-WIDE CODING-DENOMINATOR SEQUENCE-MASK AUDIT =====", file=sys.stderr)
        print(f"ERROR: {error}", file=sys.stderr)
        print("TRANSCRIPTOME-WIDE CODING-DENOMINATOR SEQUENCE-MASK AUDIT: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
