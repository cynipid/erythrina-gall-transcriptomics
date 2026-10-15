#!/usr/bin/env python3
"""Extract final top-50 protein/nucleotide discordance cases for review.

This evidence-only step reuses completed DIAMOND and minimap2 results. It does
not run a similarity search, change the expanded mask, exclude candidates, or
alter the frozen ranking.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path


ROOT = Path.home() / "projects" / "qe_cross_species_discovery"
RESULTS = ROOT / "05_host_focused_DESeq2_sensitivity" / "results"
INTEGRATION = ROOT / "04_cross_species_integration" / "results"
CONTROL = ROOT / "00_project_control"
LOGS = ROOT / "05_host_focused_DESeq2_sensitivity" / "logs"

TOP50 = RESULTS / "plant_GP_host_focused_expanded_mask_top50_final_quality_provenance_audit.tsv"
PROTEIN = INTEGRATION / "plant_representative_proteome_vs_Qerythrinae_best_hits.tsv"
NEAR_IDENTITY_CDS = INTEGRATION / "plant_all_protein_near_identity_CDS_best_by_gene.tsv"

REVIEW_OUTPUT = RESULTS / "plant_GP_host_focused_expanded_mask_top50_protein_nucleotide_discordance_review.tsv"
SUMMARY_OUTPUT = RESULTS / "plant_GP_host_focused_expanded_mask_top50_protein_nucleotide_discordance_summary.tsv"
PARAMETERS_OUTPUT = CONTROL / "36_expanded_mask_top50_protein_nucleotide_discordance_review_parameters.tsv"
AUDIT_OUTPUT = LOGS / "36_expanded_mask_top50_protein_nucleotide_discordance_review_audit.txt"

REVIEW_CLASSES = {
    "insect_like_near_identity_review",
    "insect_like_substantial_homology_review",
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


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    required = [TOP50, PROTEIN, NEAR_IDENTITY_CDS]
    hashes_before = {str(path): sha256(path) for path in required}
    top_fields, top_rows = read_tsv(TOP50)
    protein_fields, protein_rows = read_tsv(PROTEIN)
    near_fields, near_rows = read_tsv(NEAR_IDENTITY_CDS)

    require_columns(
        TOP50,
        top_fields,
        {
            "gene_id", "discovery_rank", "annotation_status", "Description", "Preferred_name",
            "PFAMs", "seed_ortholog", "eggNOG_OGs", "protein_prediction_available",
            "replicate_review_flag", "replicate_review_reasons", "GP1_normalized_count",
            "GP2_normalized_count", "GP3_normalized_count", "UP1_normalized_count",
            "UP2_normalized_count", "UP3_normalized_count", "primary_classical_secretome",
            "defense_groups", "MEROPS_overlay_present", "provenance_Trinity_transcripts",
            "provenance_nucleotide_review_class", "provenance_strong_nucleotide_Qerythrinae_similarity_review",
            "provenance_nucleotide_percent_identity",
            "provenance_nucleotide_query_span_coverage_percent",
            "provenance_nucleotide_target_span_coverage_percent",
            "provenance_nucleotide_alignment_block_length",
            "provenance_Qerythrinae_target_CDS_id",
            "provenance_protein_review_class", "expanded_mask_overlap",
        },
    )
    require_columns(
        PROTEIN,
        protein_fields,
        {
            "gene_id", "representative_protein_id", "representative_protein_length_aa",
            "subject_id", "percent_identity", "alignment_length", "query_length",
            "subject_length", "evalue", "bitscore", "query_span_coverage_percent",
            "subject_span_coverage_percent", "Qerythrinae_similarity_review_class",
        },
    )
    require_columns(
        NEAR_IDENTITY_CDS,
        near_fields,
        {
            "gene_id", "protein_provenance_review_class", "nucleotide_provenance_review_class",
            "multilevel_provenance_status", "confirmed_for_host_reanalysis_exclusion",
        },
    )

    top = unique_index(TOP50, top_fields, top_rows)
    protein = unique_index(PROTEIN, protein_fields, protein_rows)
    near = unique_index(NEAR_IDENTITY_CDS, near_fields, near_rows)
    if len(top_rows) != 50 or len(protein_rows) != 41_330 or len(near_rows) != 4_213:
        raise AuditError(
            f"Input counts differ: top50={len(top_rows)}, protein={len(protein_rows)}, near_identity={len(near_rows)}"
        )
    review_genes = {
        gene for gene, row in top.items()
        if clean(row.get("provenance_protein_review_class")) in REVIEW_CLASSES
    }
    class_counts = Counter(
        clean(top[gene].get("provenance_protein_review_class")) for gene in review_genes
    )
    if class_counts != Counter(
        {"insect_like_near_identity_review": 9, "insect_like_substantial_homology_review": 3}
    ):
        raise AuditError(f"Expected 9 near-identity and 3 substantial-homology genes; observed {dict(class_counts)}")

    fields = [
        "discovery_rank", "gene_id", "protein_nucleotide_review_group", "annotation_status",
        "Description", "Preferred_name", "PFAMs", "seed_ortholog", "eggNOG_OGs",
        "representative_protein_id", "representative_protein_length_aa", "Qerythrinae_subject_id",
        "protein_percent_identity", "protein_alignment_length_aa", "protein_query_length_aa",
        "protein_subject_length_aa", "protein_query_span_coverage_percent",
        "protein_subject_span_coverage_percent", "protein_evalue", "protein_bitscore",
        "protein_Qerythrinae_similarity_review_class", "nucleotide_Qerythrinae_similarity_review_class",
        "nucleotide_Qerythrinae_target_CDS_id", "nucleotide_percent_identity",
        "nucleotide_query_span_coverage_percent", "nucleotide_target_span_coverage_percent",
        "nucleotide_alignment_block_length", "nucleotide_strong_similarity_review",
        "original_multilevel_provenance_status",
        "original_confirmed_for_host_reanalysis_exclusion", "expanded_mask_overlap",
        "Trinity_transcripts", "primary_classical_secretome", "defense_groups",
        "MEROPS_overlay_present", "replicate_review_flag", "replicate_review_reasons",
        "GP1_normalized_count", "GP2_normalized_count", "GP3_normalized_count",
        "UP1_normalized_count", "UP2_normalized_count", "UP3_normalized_count",
    ]
    output_rows: list[dict[str, object]] = []
    for source in sorted(
        (top[gene] for gene in review_genes), key=lambda row: int(row["discovery_rank"])
    ):
        gene = source["gene_id"]
        protein_row = protein.get(gene)
        if protein_row is None:
            raise AuditError(f"Protein-review gene missing from protein table: {gene}")
        protein_class = clean(source.get("provenance_protein_review_class"))
        if clean(protein_row.get("Qerythrinae_similarity_review_class")) != protein_class:
            raise AuditError(f"Protein review-class disagreement for {gene}")
        if clean(source.get("expanded_mask_overlap")) != "FALSE":
            raise AuditError(f"Review gene unexpectedly overlaps expanded mask: {gene}")
        if clean(source.get("provenance_strong_nucleotide_Qerythrinae_similarity_review")) != "FALSE":
            raise AuditError(f"Review gene unexpectedly has strong nucleotide evidence: {gene}")

        if protein_class == "insect_like_near_identity_review":
            near_row = near.get(gene)
            if near_row is None:
                raise AuditError(f"Near-identity gene missing from multilevel audit: {gene}")
            if clean(near_row.get("confirmed_for_host_reanalysis_exclusion")) != "FALSE":
                raise AuditError(f"Near-identity top gene was previously confirmed for exclusion: {gene}")
            if clean(near_row.get("multilevel_provenance_status")) != "protein_near_identity_but_nucleotide_evidence_insufficient_for_confirmed_mask":
                raise AuditError(f"Unexpected multilevel status for {gene}")
            review_group = "protein_near_identity_nucleotide_evidence_insufficient"
            multilevel_status = near_row["multilevel_provenance_status"]
            multilevel_exclusion = near_row["confirmed_for_host_reanalysis_exclusion"]
        else:
            near_row = {}
            review_group = "protein_substantial_homology_context"
            multilevel_status = "not_in_protein_near_identity_multilevel_audit"
            multilevel_exclusion = "NOT_APPLICABLE"

        output_rows.append(
            {
                "discovery_rank": source["discovery_rank"],
                "gene_id": gene,
                "protein_nucleotide_review_group": review_group,
                "annotation_status": source["annotation_status"],
                "Description": source["Description"],
                "Preferred_name": source["Preferred_name"],
                "PFAMs": source["PFAMs"],
                "seed_ortholog": source["seed_ortholog"],
                "eggNOG_OGs": source["eggNOG_OGs"],
                "representative_protein_id": protein_row["representative_protein_id"],
                "representative_protein_length_aa": protein_row["representative_protein_length_aa"],
                "Qerythrinae_subject_id": protein_row["subject_id"],
                "protein_percent_identity": protein_row["percent_identity"],
                "protein_alignment_length_aa": protein_row["alignment_length"],
                "protein_query_length_aa": protein_row["query_length"],
                "protein_subject_length_aa": protein_row["subject_length"],
                "protein_query_span_coverage_percent": protein_row["query_span_coverage_percent"],
                "protein_subject_span_coverage_percent": protein_row["subject_span_coverage_percent"],
                "protein_evalue": protein_row["evalue"],
                "protein_bitscore": protein_row["bitscore"],
                "protein_Qerythrinae_similarity_review_class": protein_class,
                "nucleotide_Qerythrinae_similarity_review_class": source["provenance_nucleotide_review_class"],
                "nucleotide_Qerythrinae_target_CDS_id": source["provenance_Qerythrinae_target_CDS_id"],
                "nucleotide_percent_identity": source["provenance_nucleotide_percent_identity"],
                "nucleotide_query_span_coverage_percent": source["provenance_nucleotide_query_span_coverage_percent"],
                "nucleotide_target_span_coverage_percent": source["provenance_nucleotide_target_span_coverage_percent"],
                "nucleotide_alignment_block_length": source["provenance_nucleotide_alignment_block_length"],
                "nucleotide_strong_similarity_review": "FALSE",
                "original_multilevel_provenance_status": multilevel_status,
                "original_confirmed_for_host_reanalysis_exclusion": multilevel_exclusion,
                "expanded_mask_overlap": "FALSE",
                "Trinity_transcripts": source["provenance_Trinity_transcripts"],
                "primary_classical_secretome": source["primary_classical_secretome"],
                "defense_groups": source["defense_groups"],
                "MEROPS_overlay_present": source["MEROPS_overlay_present"],
                "replicate_review_flag": source["replicate_review_flag"],
                "replicate_review_reasons": source["replicate_review_reasons"],
                "GP1_normalized_count": source["GP1_normalized_count"],
                "GP2_normalized_count": source["GP2_normalized_count"],
                "GP3_normalized_count": source["GP3_normalized_count"],
                "UP1_normalized_count": source["UP1_normalized_count"],
                "UP2_normalized_count": source["UP2_normalized_count"],
                "UP3_normalized_count": source["UP3_normalized_count"],
            }
        )
    write_tsv(REVIEW_OUTPUT, output_rows, fields)

    summaries: list[dict[str, object]] = []
    dimensions = {
        "protein_review_class": Counter(str(row["protein_Qerythrinae_similarity_review_class"]) for row in output_rows),
        "nucleotide_review_class": Counter(str(row["nucleotide_Qerythrinae_similarity_review_class"]) for row in output_rows),
        "annotation_status": Counter(str(row["annotation_status"]) for row in output_rows),
        "replicate_review_flag": Counter(str(row["replicate_review_flag"]) for row in output_rows),
        "primary_classical_secretome": Counter(str(row["primary_classical_secretome"]) for row in output_rows),
    }
    for dimension, counts in dimensions.items():
        for category, count in sorted(counts.items()):
            summaries.append(
                {"summary_dimension": dimension, "category": category, "candidate_genes": count}
            )
    write_tsv(SUMMARY_OUTPUT, summaries, ["summary_dimension", "category", "candidate_genes"])

    parameters = [
        {"parameter": "scope", "value": "12 final top-50 genes with protein near-identity or substantial-homology review labels"},
        {"parameter": "sequence_search", "value": "none; completed DIAMOND and minimap2 evidence reused"},
        {"parameter": "near_identity_definition", "value": "protein identity >=90% and query/subject span coverage >=80%"},
        {"parameter": "mask_policy", "value": "near-identity protein evidence required high/substantial nucleotide concordance for the original confirmed multilevel exclusion"},
        {"parameter": "interpretive_policy", "value": "evidence review only; no filtering, taxonomic assignment or reranking"},
    ]
    parameters.extend(
        {"parameter": f"SHA256_{Path(path).name}", "value": digest}
        for path, digest in sorted(hashes_before.items())
    )
    write_tsv(PARAMETERS_OUTPUT, parameters, ["parameter", "value"])
    if any(sha256(Path(path)) != digest for path, digest in hashes_before.items()):
        raise AuditError("One or more frozen inputs changed during output generation")

    annotation_counts = Counter(str(row["annotation_status"]) for row in output_rows)
    nucleotide_counts = Counter(str(row["nucleotide_Qerythrinae_similarity_review_class"]) for row in output_rows)
    replicate_flags = sum(str(row["replicate_review_flag"]) == "TRUE" for row in output_rows)
    secretome = sum(str(row["primary_classical_secretome"]) == "TRUE" for row in output_rows)
    lines = [
        "===== EXPANDED-MASK TOP-50 PROTEIN/NUCLEOTIDE DISCORDANCE REVIEW =====",
        "Scope: focused extraction of existing evidence; no sequence search, filtering, exclusion or reranking.",
        "",
        f"Protein-review genes extracted: {len(output_rows)}",
        f"Protein near-identity with insufficient nucleotide evidence: {class_counts['insect_like_near_identity_review']}",
        f"Protein substantial-homology context: {class_counts['insect_like_substantial_homology_review']}",
        "Expanded-mask overlap: 0",
        "Strong nucleotide-similarity flags: 0",
        "Previously confirmed multilevel exclusions retained: 0",
        f"Replicate-profile review flags: {replicate_flags}",
        f"Primary classical-secretome members: {secretome}",
        "",
        "Annotation status:",
    ]
    for label, count in sorted(annotation_counts.items()):
        lines.append(f"  {label}: {count}")
    lines.extend(["", "Nucleotide evidence classes:"])
    for label, count in sorted(nucleotide_counts.items()):
        lines.append(f"  {label}: {count}")
    lines.extend(
        [
            "",
            "Focused cases (rank, gene, protein identity/alignment/query length, protein coverage, nucleotide identity/alignment/query coverage, nucleotide class, annotation):",
        ]
    )
    for row in output_rows:
        description = clean(row["Description"]) or clean(row["Preferred_name"]) or "no informative label"
        lines.append(
            "  "
            + " | ".join(
                [
                    str(row["discovery_rank"]),
                    str(row["gene_id"]),
                    f"protein={row['protein_percent_identity']}%/{row['protein_alignment_length_aa']}aa/query{row['protein_query_length_aa']}aa",
                    f"qcov={row['protein_query_span_coverage_percent']}%",
                    f"scov={row['protein_subject_span_coverage_percent']}%",
                    f"nucleotide={row['nucleotide_percent_identity'] or 'none'}%/{row['nucleotide_alignment_block_length'] or 'none'}nt/qcov={row['nucleotide_query_span_coverage_percent'] or 'none'}%",
                    str(row["nucleotide_Qerythrinae_similarity_review_class"]),
                    description,
                ]
            )
        )
    lines.extend(
        [
            "",
            "Interpretive boundary:",
            "  Protein near identity without nucleotide concordance remains discordant evidence, not proof of plant or insect origin.",
            "  Universally conserved proteins and short proteins require special caution when interpreting percentage identity.",
            "  Uncharacterized labels do not support novelty or gall-function claims.",
            "  The 12 genes remain members of the frozen top 50 unless a separate, explicitly justified policy is adopted.",
            "",
            "Frozen-input SHA-256 verification after output generation: PASS",
            "Errors detected: 0",
            "EXPANDED-MASK TOP-50 PROTEIN/NUCLEOTIDE DISCORDANCE REVIEW: PASS",
        ]
    )
    AUDIT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT_OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    for path in [REVIEW_OUTPUT, SUMMARY_OUTPUT, PARAMETERS_OUTPUT, AUDIT_OUTPUT]:
        print(f"Output: {path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("===== EXPANDED-MASK TOP-50 PROTEIN/NUCLEOTIDE DISCORDANCE REVIEW =====", file=sys.stderr)
        print(f"ERROR: {error}", file=sys.stderr)
        print("EXPANDED-MASK TOP-50 PROTEIN/NUCLEOTIDE DISCORDANCE REVIEW: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
