#!/usr/bin/env python3
"""Build the final quality/provenance audit for the frozen plant top 50.

This script reuses the completed transcriptome-wide Q. erythrinae similarity
landscape, expanded mask, and replicate audit. It runs no sequence search,
changes no DEG or candidate membership, and does not filter or rerank genes.
"""

from __future__ import annotations

import csv
import hashlib
import math
import sys
from collections import Counter
from pathlib import Path


ROOT = Path.home() / "projects" / "qe_cross_species_discovery"
RESULTS = ROOT / "05_host_focused_DESeq2_sensitivity" / "results"
INTEGRATION = ROOT / "04_cross_species_integration" / "results"
CONTROL = ROOT / "00_project_control"
LOGS = ROOT / "05_host_focused_DESeq2_sensitivity" / "logs"

TOP50 = RESULTS / "plant_GP_host_focused_expanded_mask_top50_replicate_audit.tsv"
POOR15 = RESULTS / "plant_GP_host_focused_expanded_mask_poorly_annotated_top15_replicate_audit.tsv"
GENE_BEST = INTEGRATION / "plant_transcriptome_vs_Qerythrinae_CDS_best_by_gene.tsv"
EXPANDED_MASK = INTEGRATION / "plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv"

TOP50_OUTPUT = RESULTS / "plant_GP_host_focused_expanded_mask_top50_final_quality_provenance_audit.tsv"
POOR15_OUTPUT = RESULTS / "plant_GP_host_focused_expanded_mask_poorly_annotated_top15_final_quality_provenance_audit.tsv"
SUMMARY_OUTPUT = RESULTS / "plant_GP_host_focused_expanded_mask_top50_final_quality_provenance_summary.tsv"
PARAMETERS_OUTPUT = CONTROL / "35_expanded_mask_host_focused_plant_top50_final_quality_provenance_parameters.tsv"
AUDIT_OUTPUT = LOGS / "35_expanded_mask_host_focused_plant_top50_final_quality_provenance_audit.txt"

EXPECTED_TOP50 = 50
EXPECTED_POOR15 = 15
EXPECTED_TRINITY_GENES = 118_082
EXPECTED_EXPANDED_MASK = 16_112
STRONG_NUCLEOTIDE_CLASSES = {
    "high_identity_Qerythrinae_CDS_review",
    "substantial_Qerythrinae_CDS_similarity_review",
}
MASK_EXCLUDE_STATUSES = {
    "review_exclude_existing_confirmed_multilevel_sequence",
    "review_exclude_additional_strong_nucleotide_similarity",
}


class AuditError(RuntimeError):
    pass


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip().lstrip("\ufeff")


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


def numeric(value: object, label: str) -> float:
    try:
        result = float(clean(value))
    except (TypeError, ValueError) as exc:
        raise AuditError(f"Non-numeric value for {label}: {value!r}") from exc
    if not math.isfinite(result):
        raise AuditError(f"Non-finite value for {label}: {value!r}")
    return result


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def transcript_count_bin(count: int) -> str:
    if count == 1:
        return "one_Trinity_transcript"
    if count <= 4:
        return "two_to_four_Trinity_transcripts"
    return "five_or_more_Trinity_transcripts"


def review_dimensions(
    top: dict[str, str], provenance: dict[str, str]
) -> list[str]:
    flags: list[str] = []
    annotation = clean(top.get("annotation_status"))
    if annotation == "uncharacterized":
        flags.append("functional_annotation_uncharacterized")
    elif annotation == "poorly_characterized":
        flags.append("functional_annotation_poorly_characterized")
    if clean(top.get("protein_prediction_available")) != "TRUE":
        flags.append("representative_protein_unavailable")
    if clean(top.get("replicate_review_flag")) == "TRUE":
        flags.append("replicate_profile_review")
    nucleotide_class = clean(provenance.get("nucleotide_provenance_review_class"))
    if nucleotide_class == "partial_Qerythrinae_CDS_similarity":
        flags.append("partial_nucleotide_cross_species_homology_context")
    protein_class = clean(provenance.get("protein_Qerythrinae_similarity_review_class"))
    if protein_class in {
        "insect_like_substantial_homology_review",
        "partial_homology_to_Qerythrinae_reference",
    }:
        flags.append("protein_cross_species_homology_context")
    return flags


def main() -> int:
    required = [TOP50, POOR15, GENE_BEST, EXPANDED_MASK]
    hashes_before = {str(path): sha256(path) for path in required}
    top_fields, top_rows = read_tsv(TOP50)
    poor_fields, poor_rows = read_tsv(POOR15)
    best_fields, best_rows = read_tsv(GENE_BEST)
    mask_fields, mask_rows = read_tsv(EXPANDED_MASK)

    require_columns(
        TOP50,
        top_fields,
        {
            "gene_id", "discovery_rank", "annotation_status", "protein_prediction_available",
            "primary_classical_secretome", "replicate_review_flag", "replicate_review_reasons",
            "DESeq2_status_host_focused_expanded_mask",
        },
    )
    require_columns(POOR15, poor_fields, {"gene_id", "discovery_rank", "annotation_status"})
    require_columns(
        GENE_BEST,
        best_fields,
        {
            "gene_id", "Trinity_transcripts", "best_query_transcript_id", "target_CDS_id",
            "percent_identity", "query_span_coverage_percent", "target_span_coverage_percent",
            "alignment_block_length", "mapping_quality", "nucleotide_provenance_review_class",
            "strong_nucleotide_Qerythrinae_similarity_review", "representative_protein_available",
            "protein_Qerythrinae_similarity_review_class", "current_confirmed_mask_overlap",
        },
    )
    require_columns(
        EXPANDED_MASK,
        mask_fields,
        {"gene_id", "nucleotide_provenance_review_class", "proposed_expanded_mask_status"},
    )

    top = unique_index(TOP50, top_fields, top_rows)
    poor = unique_index(POOR15, poor_fields, poor_rows)
    best = unique_index(GENE_BEST, best_fields, best_rows)
    mask = unique_index(EXPANDED_MASK, mask_fields, mask_rows)
    if len(top_rows) != EXPECTED_TOP50:
        raise AuditError(f"Top-50 rows: {len(top_rows)}; expected 50")
    if len(poor_rows) != EXPECTED_POOR15 or not set(poor).issubset(top):
        raise AuditError("Poorly annotated top-15 membership is invalid")
    if len(best_rows) != EXPECTED_TRINITY_GENES or len(mask_rows) != EXPECTED_TRINITY_GENES:
        raise AuditError(
            f"Transcriptome-wide universe differs: gene_best={len(best_rows)}, mask={len(mask_rows)}"
        )
    if set(best) != set(mask):
        raise AuditError("Gene-best and expanded-mask universes differ")

    excluded = {
        gene for gene, row in mask.items()
        if clean(row.get("proposed_expanded_mask_status")) in MASK_EXCLUDE_STATUSES
    }
    if len(excluded) != EXPECTED_EXPANDED_MASK:
        raise AuditError(f"Expanded-mask genes: {len(excluded)}; expected 16112")
    if set(top) & excluded:
        raise AuditError("One or more final top-50 genes overlap the expanded mask")
    if [int(numeric(row["discovery_rank"], row["gene_id"])) for row in top_rows] != list(range(1, 51)):
        raise AuditError("Top-50 ranks are not exactly 1 through 50")

    appended_fields = [
        "final_candidate_membership",
        "expanded_mask_overlap",
        "expanded_mask_status",
        "provenance_Trinity_transcripts",
        "provenance_Trinity_transcript_count_bin",
        "provenance_best_query_transcript_id",
        "provenance_Qerythrinae_target_CDS_id",
        "provenance_nucleotide_percent_identity",
        "provenance_nucleotide_query_span_coverage_percent",
        "provenance_nucleotide_target_span_coverage_percent",
        "provenance_nucleotide_alignment_block_length",
        "provenance_nucleotide_mapping_quality",
        "provenance_nucleotide_review_class",
        "provenance_strong_nucleotide_Qerythrinae_similarity_review",
        "provenance_protein_review_class",
        "provenance_current_confirmed_mask_overlap",
        "quality_review_dimension_count",
        "quality_review_dimensions",
    ]
    audited_rows: list[dict[str, object]] = []
    for source in top_rows:
        gene = source["gene_id"]
        provenance = best.get(gene)
        mask_row = mask.get(gene)
        if provenance is None or mask_row is None:
            raise AuditError(f"Top-50 gene missing transcriptome-wide provenance: {gene}")
        nucleotide_class = clean(provenance.get("nucleotide_provenance_review_class"))
        if nucleotide_class in STRONG_NUCLEOTIDE_CLASSES:
            raise AuditError(f"Strong nucleotide-review gene remains in final top 50: {gene}")
        if clean(provenance.get("strong_nucleotide_Qerythrinae_similarity_review")) != "FALSE":
            raise AuditError(f"Strong nucleotide flag remains TRUE for {gene}")
        if clean(mask_row.get("proposed_expanded_mask_status")) != "retain_no_strong_nucleotide_Qerythrinae_similarity":
            raise AuditError(f"Unexpected expanded-mask status for {gene}")
        if clean(provenance.get("current_confirmed_mask_overlap")) != "FALSE":
            raise AuditError(f"Current confirmed-mask overlap remains for {gene}")
        if clean(provenance.get("representative_protein_available")) != clean(source.get("protein_prediction_available")):
            raise AuditError(f"Representative-protein availability disagreement for {gene}")

        transcript_count = int(numeric(provenance.get("Trinity_transcripts"), f"{gene}/Trinity_transcripts"))
        flags = review_dimensions(source, provenance)
        output: dict[str, object] = dict(source)
        output.update(
            {
                "final_candidate_membership": "retained_frozen_expanded_mask_top50",
                "expanded_mask_overlap": "FALSE",
                "expanded_mask_status": mask_row["proposed_expanded_mask_status"],
                "provenance_Trinity_transcripts": transcript_count,
                "provenance_Trinity_transcript_count_bin": transcript_count_bin(transcript_count),
                "provenance_best_query_transcript_id": provenance["best_query_transcript_id"],
                "provenance_Qerythrinae_target_CDS_id": provenance["target_CDS_id"],
                "provenance_nucleotide_percent_identity": provenance["percent_identity"],
                "provenance_nucleotide_query_span_coverage_percent": provenance["query_span_coverage_percent"],
                "provenance_nucleotide_target_span_coverage_percent": provenance["target_span_coverage_percent"],
                "provenance_nucleotide_alignment_block_length": provenance["alignment_block_length"],
                "provenance_nucleotide_mapping_quality": provenance["mapping_quality"],
                "provenance_nucleotide_review_class": nucleotide_class,
                "provenance_strong_nucleotide_Qerythrinae_similarity_review": "FALSE",
                "provenance_protein_review_class": provenance["protein_Qerythrinae_similarity_review_class"],
                "provenance_current_confirmed_mask_overlap": "FALSE",
                "quality_review_dimension_count": len(flags),
                "quality_review_dimensions": ";".join(flags),
            }
        )
        audited_rows.append(output)

    audited = {str(row["gene_id"]): row for row in audited_rows}
    poor_rows_out = [audited[row["gene_id"]] for row in poor_rows]
    output_fields = top_fields + appended_fields
    write_tsv(TOP50_OUTPUT, audited_rows, output_fields)
    write_tsv(POOR15_OUTPUT, poor_rows_out, output_fields)

    set_rows = {
        "expanded_mask_GP_top50": audited_rows,
        "expanded_mask_GP_poorly_annotated_top15": poor_rows_out,
    }
    summary_rows: list[dict[str, object]] = []
    summary_fields = ["candidate_set", "summary_dimension", "category", "candidate_genes", "percent_of_set"]
    for set_name, rows in set_rows.items():
        dimensions = {
            "annotation_status": Counter(clean(row.get("annotation_status")) for row in rows),
            "replicate_review_flag": Counter(clean(row.get("replicate_review_flag")) for row in rows),
            "representative_protein_available": Counter(clean(row.get("protein_prediction_available")) for row in rows),
            "primary_classical_secretome": Counter(clean(row.get("primary_classical_secretome")) for row in rows),
            "nucleotide_provenance_class": Counter(clean(row.get("provenance_nucleotide_review_class")) for row in rows),
            "protein_provenance_class": Counter(clean(row.get("provenance_protein_review_class")) or "no_predicted_representative_protein" for row in rows),
            "Trinity_transcript_count_bin": Counter(clean(row.get("provenance_Trinity_transcript_count_bin")) for row in rows),
        }
        for dimension, counts in dimensions.items():
            for category, count in sorted(counts.items()):
                summary_rows.append(
                    {
                        "candidate_set": set_name,
                        "summary_dimension": dimension,
                        "category": category,
                        "candidate_genes": count,
                        "percent_of_set": format(100.0 * count / len(rows), ".8g"),
                    }
                )
    write_tsv(SUMMARY_OUTPUT, summary_rows, summary_fields)

    parameters = [
        {"parameter": "candidate_membership", "value": "frozen expanded-mask GP-up top 50"},
        {"parameter": "sequence_analysis", "value": "none; reused completed transcriptome-wide nucleotide/protein evidence"},
        {"parameter": "strong_nucleotide_definition", "value": "high identity or substantial similarity under frozen script-23 thresholds"},
        {"parameter": "Trinity_transcript_counts", "value": "descriptive; not treated as proof of fragmentation or alternative splicing"},
        {"parameter": "review_flags", "value": "descriptive; no candidate exclusion or reranking"},
    ]
    parameters.extend(
        {"parameter": f"SHA256_{Path(path).name}", "value": digest}
        for path, digest in sorted(hashes_before.items())
    )
    write_tsv(PARAMETERS_OUTPUT, parameters, ["parameter", "value"])

    if any(sha256(Path(path)) != digest for path, digest in hashes_before.items()):
        raise AuditError("One or more frozen inputs changed during output generation")

    top_annotation = Counter(clean(row.get("annotation_status")) for row in audited_rows)
    poor_annotation = Counter(clean(row.get("annotation_status")) for row in poor_rows_out)
    top_nucleotide = Counter(clean(row.get("provenance_nucleotide_review_class")) for row in audited_rows)
    poor_nucleotide = Counter(clean(row.get("provenance_nucleotide_review_class")) for row in poor_rows_out)
    top_protein = Counter(
        clean(row.get("provenance_protein_review_class")) or "no_predicted_representative_protein"
        for row in audited_rows
    )
    transcript_bins = Counter(clean(row.get("provenance_Trinity_transcript_count_bin")) for row in audited_rows)
    replicate_flags = sum(clean(row.get("replicate_review_flag")) == "TRUE" for row in audited_rows)
    poor_replicate_flags = sum(clean(row.get("replicate_review_flag")) == "TRUE" for row in poor_rows_out)
    complete_pairwise = sum(clean(row.get("complete_pairwise_separation")) == "TRUE" for row in audited_rows)
    leave_one_out = sum(clean(row.get("all_leave_one_sample_out_concordant")) == "TRUE" for row in audited_rows)
    protein_available = sum(clean(row.get("protein_prediction_available")) == "TRUE" for row in audited_rows)
    secretome = sum(clean(row.get("primary_classical_secretome")) == "TRUE" for row in audited_rows)
    defense_overlay = sum(bool(clean(row.get("defense_groups"))) for row in audited_rows)
    merops_overlay = sum(clean(row.get("MEROPS_overlay_present")) == "TRUE" for row in audited_rows)
    any_review = sum(int(row["quality_review_dimension_count"]) > 0 for row in audited_rows)

    if complete_pairwise != 50 or leave_one_out != 50:
        raise AuditError(
            f"Replicate-support accounting differs: complete_pairwise={complete_pairwise}, leave_one_out={leave_one_out}"
        )
    if replicate_flags != 6 or poor_replicate_flags != 1:
        raise AuditError(
            f"Replicate-review flags differ: top50={replicate_flags}, poor15={poor_replicate_flags}"
        )
    if top_annotation != Counter({"uncharacterized": 41, "informative": 9}):
        raise AuditError(f"Top-50 annotation-status distribution differs: {dict(top_annotation)}")
    if poor_annotation != Counter({"uncharacterized": 15}):
        raise AuditError(f"Poorly annotated top-15 distribution differs: {dict(poor_annotation)}")
    if protein_available != 30 or secretome != 11 or defense_overlay != 2 or merops_overlay != 2:
        raise AuditError(
            "Top-50 descriptive overlay counts differ: "
            f"protein={protein_available}, secretome={secretome}, defense={defense_overlay}, MEROPS={merops_overlay}"
        )

    lines = [
        "===== EXPANDED-MASK HOST-FOCUSED PLANT TOP-50 FINAL QUALITY/PROVENANCE AUDIT =====",
        "Scope: frozen top-50 evidence overlay using completed transcriptome-wide provenance and replicate results.",
        "No sequence search, DESeq2 analysis, candidate filtering, exclusion or reranking was performed.",
        "",
        f"Final top-50 genes audited: {len(audited_rows)}",
        f"Poorly annotated/uncharacterized top-15 genes audited: {len(poor_rows_out)}",
        "Expanded-mask overlap: 0",
        "Strong nucleotide Q. erythrinae-similarity flags: 0",
        "Current confirmed multilevel-mask overlap: 0",
        "",
        "Expression/replicate evidence:",
        f"  Complete GP-versus-UP pairwise separation: {complete_pairwise} / 50",
        f"  All leave-one-sample-out directions concordant: {leave_one_out} / 50",
        f"  GP single-replicate dominance review: {replicate_flags} / 50",
        f"  Dominance review within poorly annotated top 15: {poor_replicate_flags} / 15",
        "",
        "Functional/protein evidence in top 50:",
        f"  Informative annotation: {top_annotation['informative']}",
        f"  Poorly characterized: {top_annotation['poorly_characterized']}",
        f"  Uncharacterized: {top_annotation['uncharacterized']}",
        f"  Representative protein available: {protein_available}",
        f"  Representative protein unavailable: {50 - protein_available}",
        f"  Primary classical-secretome members: {secretome}",
        f"  Defense-overlay members: {defense_overlay}",
        f"  Retained MEROPS-overlay members: {merops_overlay}",
        f"  Candidates with at least one descriptive review dimension: {any_review}",
        "",
        "Nucleotide-level Q. erythrinae similarity classes in top 50:",
    ]
    for label, count in sorted(top_nucleotide.items()):
        lines.append(f"  {label}: {count}")
    lines.extend(["", "Nucleotide classes in poorly annotated top 15:"])
    for label, count in sorted(poor_nucleotide.items()):
        lines.append(f"  {label}: {count}")
    lines.extend(["", "Protein-level Q. erythrinae similarity classes in top 50:"])
    for label, count in sorted(top_protein.items()):
        lines.append(f"  {label}: {count}")
    lines.extend(["", "Trinity transcript-count description in top 50:"])
    for label, count in sorted(transcript_bins.items()):
        lines.append(f"  {label}: {count}")
    lines.extend(
        [
            "",
            "Interpretive boundary:",
            "  Absence from the expanded mask is not proof of plant origin; it means no strong Q. erythrinae match under the frozen thresholds.",
            "  Partial or ordinary conserved homology is contextual evidence and is not an exclusion criterion.",
            "  Uncharacterized status is not evidence of novelty, secretion or gall-inducing function.",
            "  Missing representative protein limits protein-level interpretation but does not invalidate gene-level expression.",
            "  Trinity transcript counts are descriptive and do not distinguish fragmentation from biological isoforms.",
            "  Replicate-dominance flags require cautious interpretation but do not change the frozen top-50 membership.",
            "",
            "Frozen-input SHA-256 verification after output generation: PASS",
            "Candidate membership preserved: PASS",
            "Errors detected: 0",
            "EXPANDED-MASK HOST-FOCUSED PLANT TOP-50 FINAL QUALITY/PROVENANCE AUDIT: PASS",
        ]
    )
    AUDIT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT_OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    for path in [TOP50_OUTPUT, POOR15_OUTPUT, SUMMARY_OUTPUT, PARAMETERS_OUTPUT, AUDIT_OUTPUT]:
        print(f"Output: {path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(
            "===== EXPANDED-MASK HOST-FOCUSED PLANT TOP-50 FINAL QUALITY/PROVENANCE AUDIT =====",
            file=sys.stderr,
        )
        print(f"ERROR: {error}", file=sys.stderr)
        print(
            "EXPANDED-MASK HOST-FOCUSED PLANT TOP-50 FINAL QUALITY/PROVENANCE AUDIT: CHECK REQUIRED",
            file=sys.stderr,
        )
        raise SystemExit(1)
