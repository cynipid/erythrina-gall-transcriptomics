#!/usr/bin/env python3
"""Independently audit the final expanded-mask plant discovery ranking.

This script is read-only with respect to all scientific result tables. It
reconstructs GP-up membership, the three percentile components, the equal-
weight discovery score, deterministic ordering, top-list selection, and the
previous-to-expanded membership comparison. It performs no DESeq2 fitting,
reranking, annotation search, or candidate-family modification.
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
LOGS = ROOT / "05_host_focused_DESeq2_sensitivity" / "logs"
CONTROL = ROOT / "00_project_control"

DESEQ = RESULTS / "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_eggnog_joined.tsv"
FULL = RESULTS / "plant_GP_host_focused_expanded_mask_discovery_full.tsv"
TOP50 = RESULTS / "plant_GP_host_focused_expanded_mask_top50.tsv"
POOR15 = RESULTS / "plant_GP_host_focused_expanded_mask_poorly_annotated_top15.tsv"
COMPARISON = RESULTS / "plant_GP_host_focused_current_vs_expanded_mask_discovery_rank_comparison.tsv"
PREVIOUS = RESULTS / "plant_GP_host_focused_discovery_full.tsv"
MASK = (
    ROOT
    / "04_cross_species_integration/results/plant_expanded_nucleotide_supported_Qerythrinae_sequence_review_mask.tsv"
)

SUMMARY = CONTROL / "33_expanded_mask_host_focused_plant_discovery_ranking_integrity_summary.tsv"
AUDIT = LOGS / "33_expanded_mask_host_focused_plant_discovery_ranking_integrity_audit.txt"

STATUS_COLUMN = "DESeq2_status_host_focused_expanded_mask"
MASK_EXCLUDE_STATUSES = {
    "review_exclude_existing_confirmed_multilevel_sequence",
    "review_exclude_additional_strong_nucleotide_similarity",
}
EXPECTED_STATUS_COUNTS = Counter(
    {"GP_up": 2349, "UP_up": 1386, "non_DE": 21022, "NA_padj_not_testable": 174}
)
ABS_TOLERANCE = 1e-12
ORDER_EQUIVALENCE_TOLERANCE = 1e-12


class AuditError(RuntimeError):
    pass


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip().lstrip("\ufeff")


def numeric(value: object, label: str) -> float:
    text = clean(value)
    try:
        result = float(text)
    except (TypeError, ValueError) as exc:
        raise AuditError(f"Non-numeric value for {label}: {text!r}") from exc
    if not math.isfinite(result):
        raise AuditError(f"Non-finite value for {label}: {text!r}")
    return result


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


def independent_percentiles(values_by_gene: dict[str, float]) -> dict[str, float]:
    """Calculate average-tie percentile ranks from frequency blocks."""
    counts = Counter(values_by_gene.values())
    before = 0
    percentile_by_value: dict[float, float] = {}
    total = len(values_by_gene)
    for value in sorted(counts):
        count = counts[value]
        first_rank = before + 1
        last_rank = before + count
        percentile_by_value[value] = ((first_rank + last_rank) / 2.0) / total
        before += count
    return {gene: percentile_by_value[value] for gene, value in values_by_gene.items()}


def close(observed: object, expected: float, label: str) -> None:
    value = numeric(observed, label)
    if not math.isclose(value, expected, rel_tol=1e-12, abs_tol=ABS_TOLERANCE):
        raise AuditError(f"Numeric mismatch for {label}: observed={value}, expected={expected}")


def truth(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def rows_identical(
    observed: list[dict[str, str]], expected: list[dict[str, str]], label: str
) -> None:
    if len(observed) != len(expected):
        raise AuditError(f"{label} row count differs: {len(observed)} versus {len(expected)}")
    for index, (left, right) in enumerate(zip(observed, expected), start=1):
        if left != right:
            differing = sorted(
                field for field in set(left) | set(right) if left.get(field, "") != right.get(field, "")
            )
            raise AuditError(
                f"{label} differs at row {index} ({left.get('gene_id', '')}); fields: {', '.join(differing[:10])}"
            )


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    required_paths = [DESEQ, FULL, TOP50, POOR15, COMPARISON, PREVIOUS, MASK]
    hashes_before = {str(path): sha256(path) for path in required_paths}

    deseq_fields, deseq_rows = read_tsv(DESEQ)
    full_fields, full_rows = read_tsv(FULL)
    top_fields, top_rows = read_tsv(TOP50)
    poor_fields, poor_rows = read_tsv(POOR15)
    comparison_fields, comparison_rows = read_tsv(COMPARISON)
    previous_fields, previous_rows = read_tsv(PREVIOUS)
    mask_fields, mask_rows = read_tsv(MASK)

    require_columns(
        DESEQ,
        deseq_fields,
        {"gene_id", "baseMean", "log2FoldChange", "padj", STATUS_COLUMN},
    )
    require_columns(
        FULL,
        full_fields,
        {
            "gene_id", "discovery_rank", "discovery_score", "percentile_significance",
            "percentile_effect_size", "percentile_expression", "metric_neg_log10_padj",
            "metric_abs_log2FoldChange", "metric_log10_baseMean_plus1", "baseMean",
            "log2FoldChange", "padj", STATUS_COLUMN, "annotation_status",
        },
    )
    require_columns(PREVIOUS, previous_fields, {"gene_id", "discovery_rank", "discovery_score"})
    require_columns(MASK, mask_fields, {"gene_id", "proposed_expanded_mask_status"})
    require_columns(
        COMPARISON,
        comparison_fields,
        {
            "gene_id", "comparison_transition", "previous_GP_up_member",
            "expanded_mask_GP_up_member", "expanded_nucleotide_mask_member",
            "previous_discovery_rank", "expanded_mask_discovery_rank",
            "rank_change_positive_means_rise", "expanded_mask_DESeq2_status",
        },
    )

    deseq = unique_index(DESEQ, deseq_fields, deseq_rows)
    full = unique_index(FULL, full_fields, full_rows)
    previous = unique_index(PREVIOUS, previous_fields, previous_rows)
    mask = unique_index(MASK, mask_fields, mask_rows)
    comparison = unique_index(COMPARISON, comparison_fields, comparison_rows)

    if len(deseq_rows) != 24_931:
        raise AuditError(f"Expanded-mask DESeq2 rows: {len(deseq_rows)}; expected 24931")
    status_counts = Counter(clean(row.get(STATUS_COLUMN)) for row in deseq_rows)
    if status_counts != EXPECTED_STATUS_COUNTS:
        raise AuditError(f"Expanded-mask DESeq2 status counts differ: {dict(status_counts)}")

    eligible = {
        gene: row for gene, row in deseq.items() if clean(row.get(STATUS_COLUMN)) == "GP_up"
    }
    if len(eligible) != 2_349 or len(full_rows) != 2_349:
        raise AuditError(
            f"GP-up/full-ranking counts differ: eligible={len(eligible)}, full={len(full_rows)}; expected 2349"
        )
    if set(eligible) != set(full):
        raise AuditError(
            f"Full-ranking membership differs from GP-up set: missing={len(set(eligible)-set(full))}, extra={len(set(full)-set(eligible))}"
        )

    significance: dict[str, float] = {}
    effect: dict[str, float] = {}
    expression: dict[str, float] = {}
    padj_by_gene: dict[str, float] = {}
    lfc_by_gene: dict[str, float] = {}
    base_by_gene: dict[str, float] = {}
    for gene, row in eligible.items():
        padj = numeric(row.get("padj"), f"{gene}/DESeq2_padj")
        lfc = numeric(row.get("log2FoldChange"), f"{gene}/DESeq2_log2FoldChange")
        base = numeric(row.get("baseMean"), f"{gene}/DESeq2_baseMean")
        if not (padj < 0.01 and lfc >= 1 and base >= 0):
            raise AuditError(f"GP-up threshold disagreement for {gene}")
        padj_by_gene[gene] = padj
        lfc_by_gene[gene] = lfc
        base_by_gene[gene] = base
        significance[gene] = -math.log10(max(padj, 1e-300))
        effect[gene] = abs(lfc)
        expression[gene] = math.log10(base + 1.0)

    p_significance = independent_percentiles(significance)
    p_effect = independent_percentiles(effect)
    p_expression = independent_percentiles(expression)
    scores = {
        gene: (p_significance[gene] + p_effect[gene] + p_expression[gene]) / 3.0
        for gene in eligible
    }
    expected_order = sorted(
        eligible,
        key=lambda gene: (
            -scores[gene], padj_by_gene[gene], -effect[gene], -base_by_gene[gene], gene
        ),
    )
    observed_order = [row["gene_id"] for row in full_rows]
    equivalent_reordered_positions: list[int] = []
    maximum_equivalent_score_difference = 0.0
    for position, (observed, expected) in enumerate(
        zip(observed_order, expected_order), start=1
    ):
        if observed == expected:
            continue
        score_difference = abs(scores[observed] - scores[expected])
        if score_difference > ORDER_EQUIVALENCE_TOLERANCE:
            raise AuditError(
                "Material ranking-order difference at rank "
                f"{position}: observed={observed} score={scores[observed]:.17g}; "
                f"expected={expected} score={scores[expected]:.17g}; "
                f"absolute_difference={score_difference:.17g}"
            )
        equivalent_reordered_positions.append(position)
        maximum_equivalent_score_difference = max(
            maximum_equivalent_score_difference, score_difference
        )

    for gene in expected_order:
        source = eligible[gene]
        ranked = full[gene]
        for field in ("baseMean", "log2FoldChange", "padj"):
            close(ranked.get(field), numeric(source.get(field), f"{gene}/source_{field}"), f"{gene}/{field}")
        close(ranked.get("metric_neg_log10_padj"), significance[gene], f"{gene}/metric_significance")
        close(ranked.get("metric_abs_log2FoldChange"), effect[gene], f"{gene}/metric_effect")
        close(ranked.get("metric_log10_baseMean_plus1"), expression[gene], f"{gene}/metric_expression")
        close(ranked.get("percentile_significance"), p_significance[gene], f"{gene}/percentile_significance")
        close(ranked.get("percentile_effect_size"), p_effect[gene], f"{gene}/percentile_effect")
        close(ranked.get("percentile_expression"), p_expression[gene], f"{gene}/percentile_expression")
        close(ranked.get("discovery_score"), scores[gene], f"{gene}/discovery_score")
        if clean(ranked.get(STATUS_COLUMN)) != "GP_up":
            raise AuditError(f"Non-GP-up status in ranking for {gene}")

    for observed_rank, gene in enumerate(observed_order, start=1):
        if int(numeric(full[gene].get("discovery_rank"), f"{gene}/rank")) != observed_rank:
            raise AuditError(f"Non-sequential discovery rank for {gene}")

    if top_fields != full_fields or poor_fields != full_fields:
        raise AuditError("Top-50 or poorly annotated top-15 header differs from the full ranking")
    if len(top_rows) != 50 or len(poor_rows) != 15:
        raise AuditError(f"Top-list sizes differ: top50={len(top_rows)}, poor15={len(poor_rows)}")
    rows_identical(top_rows, full_rows[:50], "Top-50 selection")
    top50_order_matches_independent_reconstruction = (
        observed_order[:50] == expected_order[:50]
    )
    if not top50_order_matches_independent_reconstruction:
        for position, (observed, expected) in enumerate(
            zip(observed_order[:50], expected_order[:50]), start=1
        ):
            if observed != expected and abs(scores[observed] - scores[expected]) > ORDER_EQUIVALENCE_TOLERANCE:
                raise AuditError(
                    f"Material independently reconstructed top-50 difference at rank {position}"
                )
    expected_poor = [
        row for row in full_rows
        if clean(row.get("annotation_status")) in {"poorly_characterized", "uncharacterized"}
    ][:15]
    rows_identical(poor_rows, expected_poor, "Poorly annotated top-15 selection")

    excluded = {
        gene for gene, row in mask.items()
        if clean(row.get("proposed_expanded_mask_status")) in MASK_EXCLUDE_STATUSES
    }
    if len(mask_rows) != 118_082 or len(excluded) != 16_112:
        raise AuditError(f"Mask accounting differs: universe={len(mask_rows)}, excluded={len(excluded)}")
    if set(full) & excluded:
        raise AuditError("Expanded-mask genes are present in the final ranking")
    if len(previous_rows) != 3_179:
        raise AuditError(f"Previous ranking rows: {len(previous_rows)}; expected 3179")

    old_set = set(previous)
    new_set = set(full)
    shared = old_set & new_set
    old_only = old_set - new_set
    new_only = new_set - old_set
    if (len(shared), len(old_only), len(new_only), len(old_only & excluded)) != (2343, 836, 6, 813):
        raise AuditError(
            "Ranking membership transitions differ: "
            f"shared={len(shared)}, old_only={len(old_only)}, new_only={len(new_only)}, masked={len(old_only & excluded)}"
        )
    union = old_set | new_set
    if set(comparison) != union or len(comparison_rows) != 3_185:
        raise AuditError(
            f"Comparison membership differs: rows={len(comparison_rows)}, expected union={len(union)}"
        )

    transition_counts: Counter[str] = Counter()
    for gene in union:
        row = comparison[gene]
        if gene in shared:
            transition = "retained_GP_up"
        elif gene in excluded:
            transition = "GP_up_to_excluded_by_expanded_mask"
        elif gene in old_only:
            transition = f"GP_up_to_{clean(deseq.get(gene, {}).get(STATUS_COLUMN)) or 'not_in_final_tested_universe'}"
        else:
            transition = "new_GP_up_after_refit"
        transition_counts[transition] += 1
        if clean(row.get("comparison_transition")) != transition:
            raise AuditError(f"Incorrect comparison transition for {gene}")
        if clean(row.get("previous_GP_up_member")) != truth(gene in old_set):
            raise AuditError(f"Incorrect previous membership flag for {gene}")
        if clean(row.get("expanded_mask_GP_up_member")) != truth(gene in new_set):
            raise AuditError(f"Incorrect expanded membership flag for {gene}")
        if clean(row.get("expanded_nucleotide_mask_member")) != truth(gene in excluded):
            raise AuditError(f"Incorrect mask flag for {gene}")

        old_rank = numeric(previous[gene]["discovery_rank"], f"{gene}/previous_rank") if gene in old_set else None
        new_rank = numeric(full[gene]["discovery_rank"], f"{gene}/expanded_rank") if gene in new_set else None
        if clean(row.get("previous_discovery_rank")) != (clean(previous[gene]["discovery_rank"]) if gene in old_set else ""):
            raise AuditError(f"Incorrect previous rank in comparison for {gene}")
        if clean(row.get("expanded_mask_discovery_rank")) != (clean(full[gene]["discovery_rank"]) if gene in new_set else ""):
            raise AuditError(f"Incorrect expanded rank in comparison for {gene}")
        if old_rank is not None and new_rank is not None:
            close(
                row.get("rank_change_positive_means_rise"),
                old_rank - new_rank,
                f"{gene}/rank_change",
            )
        elif clean(row.get("rank_change_positive_means_rise")):
            raise AuditError(f"Unexpected rank-change value for nonshared gene {gene}")

    expected_transitions = Counter(
        {
            "retained_GP_up": 2343,
            "GP_up_to_excluded_by_expanded_mask": 813,
            "GP_up_to_non_DE": 23,
            "new_GP_up_after_refit": 6,
        }
    )
    if transition_counts != expected_transitions:
        raise AuditError(f"Comparison transition counts differ: {dict(transition_counts)}")

    old_top50 = {
        gene for gene, row in previous.items()
        if numeric(row.get("discovery_rank"), f"{gene}/old_rank") <= 50
    }
    new_top50 = {row["gene_id"] for row in top_rows}
    if len(old_top50 & new_top50) != 35:
        raise AuditError(f"Shared top-50 genes: {len(old_top50 & new_top50)}; expected 35")

    if any(sha256(Path(path)) != digest for path, digest in hashes_before.items()):
        raise AuditError("One or more scientific input/result files changed during the audit")

    annotation_counts = Counter(clean(row.get("annotation_status")) for row in full_rows)
    top_annotation_counts = Counter(clean(row.get("annotation_status")) for row in top_rows)
    top_protein = sum(clean(row.get("protein_prediction_available")) == "TRUE" for row in top_rows)
    top_secretome = sum(clean(row.get("primary_classical_secretome")) == "TRUE" for row in top_rows)
    top_defense = sum(bool(clean(row.get("defense_groups"))) for row in top_rows)
    top_merops = sum(clean(row.get("MEROPS_overlay_present")) == "TRUE" for row in top_rows)

    summary_rows: list[dict[str, object]] = [
        {"metric": "expanded_mask_GP_up_genes_ranked", "value": len(full_rows)},
        {
            "metric": "ranking_order_exact_match",
            "value": truth(not equivalent_reordered_positions),
        },
        {"metric": "ranking_order_validation", "value": "PASS"},
        {
            "metric": "numerically_equivalent_reordered_positions",
            "value": len(equivalent_reordered_positions),
        },
        {
            "metric": "maximum_equivalent_score_difference",
            "value": maximum_equivalent_score_difference,
        },
        {"metric": "ranking_metrics_and_scores_match", "value": "PASS"},
        {"metric": "top50_exact_selection", "value": "PASS"},
        {"metric": "poorly_annotated_top15_exact_selection", "value": "PASS"},
        {"metric": "shared_previous_and_expanded_GP_up", "value": len(shared)},
        {"metric": "previous_only_GP_up", "value": len(old_only)},
        {"metric": "new_GP_up_after_refit", "value": len(new_only)},
        {"metric": "shared_top50", "value": len(old_top50 & new_top50)},
        {"metric": "top50_uncharacterized", "value": top_annotation_counts["uncharacterized"]},
        {"metric": "top50_representative_protein_available", "value": top_protein},
        {"metric": "top50_primary_classical_secretome", "value": top_secretome},
        {"metric": "top50_defense_overlay_present", "value": top_defense},
        {"metric": "top50_MEROPS_overlay_present", "value": top_merops},
    ]
    write_tsv(SUMMARY, summary_rows, ["metric", "value"])

    lines = [
        "===== EXPANDED-MASK HOST-FOCUSED PLANT DISCOVERY RANKING INTEGRITY AUDIT =====",
        "Scope: independent read-only reconstruction of ranking membership, scores, order and list selection.",
        "No DESeq2 fitting, reranking, annotation search or candidate-family modification was performed.",
        "",
        f"Expanded-mask DESeq2 rows: {len(deseq_rows)}",
        f"Official expanded-mask GP-up genes: {len(eligible)}",
        f"Full discovery-ranking rows: {len(full_rows)}",
        "Expanded-mask genes retained in ranking: 0",
        "",
        "Independent ranking reconstruction:",
        "  GP-up membership: PASS",
        "  Three percentile components: PASS",
        "  Equal-weight discovery scores: PASS",
        "  Ranking order: PASS",
        f"  Numerically equivalent reordered positions (score difference <= {ORDER_EQUIVALENCE_TOLERANCE:g}): {len(equivalent_reordered_positions)}",
        f"  Maximum score difference within equivalent reordered positions: {maximum_equivalent_score_difference:.17g}",
        "  Material ranking-order reversals: 0",
        "  Sequential ranks: PASS",
        "  DESeq2 numeric-field agreement: PASS",
        "  Top-50 exact selection: PASS",
        "  Independent top-50 order: "
        + ("exact match" if top50_order_matches_independent_reconstruction else "numerically equivalent match"),
        "  Poorly annotated top-15 exact selection: PASS",
        "",
        "Previous-to-expanded comparison:",
        f"  Retained GP-up genes: {len(shared)}",
        f"  Previous-only GP-up genes: {len(old_only)}",
        f"  Direct mask exclusions: {len(old_only & excluded)}",
        f"  GP-up to non-DE after refit: {transition_counts['GP_up_to_non_DE']}",
        f"  Newly GP-up after refit: {len(new_only)}",
        f"  Shared top-50 genes: {len(old_top50 & new_top50)}",
        f"  Top-50 exits: {len(old_top50 - new_top50)}",
        f"  Top-50 entrants: {len(new_top50 - old_top50)}",
        "  Comparison-table reconstruction: PASS",
        "",
        "Descriptive composition of all 2,349 ranked genes:",
        f"  Informative: {annotation_counts['informative']}",
        f"  Poorly characterized: {annotation_counts['poorly_characterized']}",
        f"  Uncharacterized: {annotation_counts['uncharacterized']}",
        "",
        "Descriptive composition of the top 50:",
        f"  Informative: {top_annotation_counts['informative']}",
        f"  Poorly characterized: {top_annotation_counts['poorly_characterized']}",
        f"  Uncharacterized: {top_annotation_counts['uncharacterized']}",
        f"  Representative-protein prediction available: {top_protein}",
        f"  Primary classical-secretome members: {top_secretome}",
        f"  Defense-overlay members: {top_defense}",
        f"  Retained MEROPS-overlay members: {top_merops}",
        "",
        "Interpretive boundary:",
        "  The score is a prioritization device, not a new significance test.",
        "  High placement of uncharacterized genes is not functional evidence.",
        "  Protein-prediction absence may reflect lack of a representative predicted protein and is not evidence of noncoding status.",
        "  MEROPS absence remains non-negative evidence because the retained screen had restricted original scope.",
        "  Replicate-level robustness must be audited separately before individual top genes are prioritized.",
        "",
        "Scientific input/result SHA-256 verification after audit: PASS",
        "Errors detected: 0",
        "EXPANDED-MASK HOST-FOCUSED PLANT DISCOVERY RANKING INTEGRITY AUDIT: PASS",
    ]
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"Output: {SUMMARY}")
    print(f"Output: {AUDIT}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(
            "===== EXPANDED-MASK HOST-FOCUSED PLANT DISCOVERY RANKING INTEGRITY AUDIT =====",
            file=sys.stderr,
        )
        print(f"ERROR: {error}", file=sys.stderr)
        print(
            "EXPANDED-MASK HOST-FOCUSED PLANT DISCOVERY RANKING INTEGRITY AUDIT: CHECK REQUIRED",
            file=sys.stderr,
        )
        raise SystemExit(1)
