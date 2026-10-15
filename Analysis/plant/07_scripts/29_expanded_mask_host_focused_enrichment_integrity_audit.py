#!/usr/bin/env python3
"""Independent integrity audit of expanded-mask plant GO/KEGG enrichment.

The script does not modify DESeq2 or enrichment result tables. It rebuilds the
valid-padj annotation backgrounds and foregrounds, independently recalculates
the one-sided hypergeometric/Fisher probabilities and BH adjustments, and
inventories GO identifiers unresolved by the current GO.db annotation step.
"""

from __future__ import annotations

import csv
import hashlib
import math
import sys
from collections import Counter
from pathlib import Path


ROOT = Path.home() / "projects" / "qe_cross_species_discovery"
BRANCH = ROOT / "05_host_focused_DESeq2_sensitivity"
RESULTS = BRANCH / "results"
ENRICHMENT = BRANCH / "enrichment"
LOGS = BRANCH / "logs"
CONTROL = ROOT / "00_project_control"
PLANT = Path.home() / "projects" / "qe_plant_denovo"

DESEQ = RESULTS / "plant_gene_DESeq2_GP_vs_UP_host_focused_expanded_mask_eggnog_joined.tsv"
GO_MAP = (
    PLANT / "04_enrichment" / "primary_FDR001_LFC1" / "go"
    / "plant_tested_gene_to_GO.tsv"
)
KEGG_MAP = (
    PLANT / "04_enrichment" / "primary_FDR001_LFC1" / "kegg"
    / "plant_tested_gene_to_KEGG_pathway_canonical.tsv"
)

TABLES = {
    ("GO", "GP"): ENRICHMENT / "host_focused_expanded_mask_GO_GP_up_all_terms_named.tsv",
    ("GO", "UP"): ENRICHMENT / "host_focused_expanded_mask_GO_UP_up_all_terms_named.tsv",
    ("KEGG", "GP"): ENRICHMENT / "host_focused_expanded_mask_KEGG_GP_up_all_pathways_named.tsv",
    ("KEGG", "UP"): ENRICHMENT / "host_focused_expanded_mask_KEGG_UP_up_all_pathways_named.tsv",
}

EXPECTED = {
    ("GO", "GP"): {"foreground": 840, "background": 11431, "tested": 5667, "significant": 64},
    ("GO", "UP"): {"foreground": 572, "background": 11431, "tested": 5667, "significant": 34},
    ("KEGG", "GP"): {"foreground": 433, "background": 5988, "tested": 358, "significant": 6},
    ("KEGG", "UP"): {"foreground": 239, "background": 5988, "tested": 358, "significant": 3},
}

# These five IDs were visible among the reported leading categories and are
# obsolete in the current ontology. Any additional unresolved IDs are retained
# in the review output without an assumed historical label.
KNOWN_OBSOLETE = {
    "GO:0055114": ("obsolete oxidation-reduction process", "BP"),
    "GO:0031226": ("obsolete intrinsic component of plasma membrane", "CC"),
    "GO:1901575": ("obsolete organic substance catabolic process", "BP"),
    "GO:0090305": ("obsolete nucleic acid phosphodiester bond hydrolysis", "BP"),
    "GO:0097159": ("obsolete organic cyclic compound binding", "MF"),
}

AUDIT = LOGS / "29_expanded_mask_host_focused_enrichment_integrity_audit.txt"
SUMMARY = CONTROL / "29_expanded_mask_host_focused_enrichment_integrity_summary.tsv"
GO_REVIEW = ENRICHMENT / "host_focused_expanded_mask_unresolved_obsolete_GO_term_review.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value: str | None) -> float:
    if value is None or value.strip() in {"", "NA", "NaN"}:
        return math.nan
    return float(value)


def isclose(observed: float, expected: float, rel: float = 2e-7, abs_: float = 1e-12) -> bool:
    return math.isclose(observed, expected, rel_tol=rel, abs_tol=abs_)


def clean_mapping(rows: list[dict[str, str]]) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for row in rows:
        gene = row.get("gene_id", "").strip()
        term = row.get("term_id", "").strip()
        if gene and gene != "NA" and term and term not in {"NA", "-"}:
            pairs.add((gene, term))
    return pairs


def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_upper_tail(a: int, n: int, k: int, population: int) -> float:
    """P[X >= a] for X~Hypergeom(population, k, n)."""
    lower = max(0, n - (population - k))
    upper = min(n, k)
    if a <= lower:
        return 1.0
    if a > upper:
        return 0.0
    log_denom = log_comb(population, n)
    logs = [log_comb(k, x) + log_comb(population - k, n - x) - log_denom
            for x in range(a, upper + 1)]
    maximum = max(logs)
    if maximum == -math.inf:
        return 0.0
    result = math.exp(maximum) * math.fsum(math.exp(value - maximum) for value in logs)
    return min(1.0, result)


def bh_adjust(pvalues: list[float]) -> list[float]:
    count = len(pvalues)
    order = sorted(range(count), key=lambda idx: (pvalues[idx], idx))
    adjusted = [math.nan] * count
    running = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, pvalues[original_index] * count / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def build_term_sets(pairs: set[tuple[str, str]], valid: set[str]) -> tuple[set[str], dict[str, set[str]]]:
    background = {gene for gene, _ in pairs if gene in valid}
    terms: dict[str, set[str]] = {}
    for gene, term in pairs:
        if gene in background:
            terms.setdefault(term, set()).add(gene)
    terms = {term: genes for term, genes in terms.items() if len(genes) >= 5}
    return background, terms


def main() -> int:
    LOGS.mkdir(parents=True, exist_ok=True)
    CONTROL.mkdir(parents=True, exist_ok=True)
    ENRICHMENT.mkdir(parents=True, exist_ok=True)

    lines = ["===== EXPANDED-MASK HOST-FOCUSED PLANT ENRICHMENT INTEGRITY AUDIT ====="]
    errors: list[str] = []
    warnings: list[str] = []
    summary_rows: list[dict[str, object]] = []

    required = [DESEQ, GO_MAP, KEGG_MAP, *TABLES.values()]
    missing = [path for path in required if not path.is_file()]
    if missing:
        errors.extend(f"Missing required file: {path}" for path in missing)
        deseq_rows: list[dict[str, str]] = []
    else:
        deseq_rows = read_tsv(DESEQ)

    if not missing:
        gene_ids = [row.get("gene_id", "") for row in deseq_rows]
        if len(deseq_rows) != 24931:
            errors.append(f"Expanded-mask DESeq2 rows: {len(deseq_rows)}; expected 24931")
        if len(gene_ids) != len(set(gene_ids)):
            errors.append("Expanded-mask DESeq2 gene IDs are not unique")

        valid: set[str] = set()
        gp: set[str] = set()
        up: set[str] = set()
        reconstructed = Counter()
        stored = Counter(row.get("DESeq2_status_host_focused_expanded_mask", "") for row in deseq_rows)
        for row in deseq_rows:
            gene = row["gene_id"]
            padj = fnum(row.get("padj"))
            lfc = fnum(row.get("log2FoldChange"))
            if math.isnan(padj):
                status = "NA_padj_not_testable"
            elif padj < 0.01 and lfc >= 1:
                status = "GP_up"
                gp.add(gene)
                valid.add(gene)
            elif padj < 0.01 and lfc <= -1:
                status = "UP_up"
                up.add(gene)
                valid.add(gene)
            else:
                status = "non_DE"
                valid.add(gene)
            reconstructed[status] += 1
            if row.get("DESeq2_status_host_focused_expanded_mask", "") != status:
                errors.append(f"Stored status mismatch for {gene}")
                break

        expected_statuses = Counter({
            "GP_up": 2349, "UP_up": 1386,
            "non_DE": 21022, "NA_padj_not_testable": 174,
        })
        if stored != expected_statuses:
            errors.append(f"Stored DESeq2 status counts differ: {dict(stored)}")
        if reconstructed != expected_statuses:
            errors.append(f"Reconstructed DESeq2 status counts differ: {dict(reconstructed)}")
        if len(valid) != 24757:
            errors.append(f"Valid-padj genes: {len(valid)}; expected 24757")

        go_pairs = clean_mapping(read_tsv(GO_MAP))
        if len(go_pairs) != 684983:
            errors.append(f"GO source pairs: {len(go_pairs)}; expected 684983")
        if len({g for g, _ in go_pairs}) != 11440:
            errors.append("GO source gene accounting differs from 11440")
        if len({t for _, t in go_pairs}) != 11253:
            errors.append("GO source term accounting differs from 11253")

        raw_kegg_pairs = clean_mapping(read_tsv(KEGG_MAP))
        if len(raw_kegg_pairs) != 24683:
            errors.append(f"KEGG source pairs: {len(raw_kegg_pairs)}; expected 24683")
        kegg_pairs = {
            (gene, "map00650" if term == "map00072" else
             "map01110" if term == "map01130" else term)
            for gene, term in raw_kegg_pairs
        }
        if len(kegg_pairs) != 24264:
            errors.append(f"Canonical KEGG pairs: {len(kegg_pairs)}; expected 24264")
        if len({t for _, t in kegg_pairs}) != 391:
            errors.append("Canonical KEGG term accounting differs from 391")

        go_background, go_terms = build_term_sets(go_pairs, valid)
        kegg_background, kegg_terms = build_term_sets(kegg_pairs, valid)

        if len(go_background) != 11431 or len(go_terms) != 5667:
            errors.append(
                f"Rebuilt GO background/tested terms: {len(go_background)}/{len(go_terms)}; "
                "expected 11431/5667"
            )
        if len(kegg_background) != 5988 or len(kegg_terms) != 358:
            errors.append(
                f"Rebuilt KEGG background/tested terms: {len(kegg_background)}/{len(kegg_terms)}; "
                "expected 5988/358"
            )

        unresolved: dict[str, dict[str, object]] = {}
        for key, path in TABLES.items():
            kind, direction = key
            expected = EXPECTED[key]
            rows = read_tsv(path)
            mapping_terms = go_terms if kind == "GO" else kegg_terms
            background = go_background if kind == "GO" else kegg_background
            direction_ids = gp if direction == "GP" else up
            foreground = direction_ids & background
            term_ids = [row.get("term_id", "") for row in rows]

            if len(rows) != expected["tested"]:
                errors.append(f"{direction}_{kind} rows: {len(rows)}; expected {expected['tested']}")
            if len(term_ids) != len(set(term_ids)):
                errors.append(f"{direction}_{kind} contains duplicate term IDs")
            if set(term_ids) != set(mapping_terms):
                errors.append(f"{direction}_{kind} term membership differs from the rebuilt mapping")
            if len(foreground) != expected["foreground"]:
                errors.append(
                    f"{direction}_{kind} foreground: {len(foreground)}; expected {expected['foreground']}"
                )
            if len(background) != expected["background"]:
                errors.append(
                    f"{direction}_{kind} background: {len(background)}; expected {expected['background']}"
                )

            pvalues: list[float] = []
            for row in rows:
                term = row["term_id"]
                term_genes = mapping_terms.get(term, set())
                observed_a = int(row["foreground_genes_with_term"])
                expected_a = len(foreground & term_genes)
                observed_k = int(row["background_genes_with_term"])
                expected_k = len(term_genes)
                if observed_a != expected_a or observed_k != expected_k:
                    errors.append(f"{direction}_{kind} count mismatch for {term}")
                    break
                if int(row["foreground_gene_set_size"]) != len(foreground):
                    errors.append(f"{direction}_{kind} foreground-size field mismatch for {term}")
                    break
                if int(row["background_gene_set_size"]) != len(background):
                    errors.append(f"{direction}_{kind} background-size field mismatch for {term}")
                    break

                expected_fold = 0.0 if expected_a == 0 else (
                    (expected_a / len(foreground)) / (expected_k / len(background))
                )
                if not isclose(fnum(row.get("fold_enrichment")), expected_fold, rel=2e-8):
                    errors.append(f"{direction}_{kind} fold-enrichment mismatch for {term}")
                    break

                expected_p = hypergeom_upper_tail(
                    expected_a, len(foreground), expected_k, len(background)
                )
                observed_p = fnum(row.get("pvalue"))
                if not isclose(observed_p, expected_p):
                    errors.append(
                        f"{direction}_{kind} Fisher P mismatch for {term}: "
                        f"observed {observed_p}, expected {expected_p}"
                    )
                    break
                pvalues.append(observed_p)

                if kind == "GO" and row.get("GO_term") == "Unresolved in current GO.db":
                    entry = unresolved.setdefault(term, {
                        "term_id": term,
                        "current_table_label": row.get("GO_term", ""),
                        "historical_or_current_review_name": KNOWN_OBSOLETE.get(term, ("unresolved identifier", "Unresolved"))[0],
                        "ontology_review": KNOWN_OBSOLETE.get(term, ("unresolved identifier", "Unresolved"))[1],
                        "current_review_status": "known_obsolete" if term in KNOWN_OBSOLETE else "requires_manual_GO_review",
                        "significant_in": [],
                        "minimum_padj": 1.0,
                        "interpretation_policy": "exclude from manuscript biological interpretation unless remapped to an active GO term",
                    })
                    if fnum(row.get("padj")) < 0.05:
                        entry["significant_in"].append(f"{direction}_up")
                    entry["minimum_padj"] = min(entry["minimum_padj"], fnum(row.get("padj")))

            if len(pvalues) == len(rows):
                recalculated = bh_adjust(pvalues)
                for row, expected_padj in zip(rows, recalculated):
                    if not isclose(fnum(row.get("padj")), expected_padj, rel=3e-7):
                        errors.append(
                            f"{direction}_{kind} BH-adjusted P mismatch for {row['term_id']}"
                        )
                        break

            significant = sum(fnum(row.get("padj")) < 0.05 for row in rows)
            if significant != expected["significant"]:
                errors.append(
                    f"{direction}_{kind} significant categories: {significant}; "
                    f"expected {expected['significant']}"
                )

            # R orders using full in-memory precision, whereas write.table can
            # make nearly tied values identical in the saved text. Validate
            # monotonic numeric ordering with tolerance; do not demand a new
            # term-ID order for values tied only after serialization.
            for previous, current in zip(rows, rows[1:]):
                previous_padj = fnum(previous.get("padj"))
                current_padj = fnum(current.get("padj"))
                padj_tied = isclose(previous_padj, current_padj, rel=1e-10, abs_=1e-13)
                if current_padj < previous_padj and not padj_tied:
                    errors.append(
                        f"{direction}_{kind} has a material padj ordering reversal: "
                        f"{previous['term_id']} -> {current['term_id']}"
                    )
                    break
                if padj_tied:
                    previous_p = fnum(previous.get("pvalue"))
                    current_p = fnum(current.get("pvalue"))
                    p_tied = isclose(previous_p, current_p, rel=1e-10, abs_=1e-13)
                    if current_p < previous_p and not p_tied:
                        errors.append(
                            f"{direction}_{kind} has a material pvalue ordering reversal "
                            f"within a padj tie: {previous['term_id']} -> {current['term_id']}"
                        )
                        break

            summary_rows.append({
                "analysis": f"{direction}_{kind}",
                "table_path": str(path),
                "tested_categories": len(rows),
                "annotation_eligible_foreground": len(foreground),
                "annotation_background": len(background),
                "significant_categories_BH_lt_0.05": significant,
                "sha256": sha256(path),
                "status": "PASS" if not errors else "CHECK",
            })
            lines.extend([
                "",
                f"{direction}_{kind}:",
                f"  Tested categories: {len(rows)}",
                f"  Annotation-eligible foreground genes: {len(foreground)}",
                f"  Annotation background genes: {len(background)}",
                f"  Significant categories (BH < 0.05): {significant}",
                f"  SHA-256: {sha256(path)}",
            ])

        review_rows: list[dict[str, object]] = []
        for term in sorted(unresolved):
            entry = unresolved[term]
            entry["significant_in"] = ",".join(entry["significant_in"]) or "none"
            review_rows.append(entry)
        write_tsv(
            GO_REVIEW,
            review_rows,
            [
                "term_id", "current_table_label", "historical_or_current_review_name",
                "ontology_review", "current_review_status", "significant_in",
                "minimum_padj", "interpretation_policy",
            ],
        )
        if review_rows:
            known_obsolete_count = sum(
                row["current_review_status"] == "known_obsolete" for row in review_rows
            )
            manual_review_count = len(review_rows) - known_obsolete_count
            significant_gp_unresolved = sum(
                "GP_up" in str(row["significant_in"]).split(",") for row in review_rows
            )
            significant_up_unresolved = sum(
                "UP_up" in str(row["significant_in"]).split(",") for row in review_rows
            )
            significant_any_unresolved = sum(
                row["significant_in"] != "none" for row in review_rows
            )
            warnings.append(
                f"{len(review_rows)} GO identifiers were unresolved by current GO.db; "
                "they were inventoried and must not be interpreted as current GO terms"
            )
            lines.extend([
                "",
                "GO ontology-version review:",
                f"  Unresolved GO identifiers: {len(review_rows)}",
                f"  Confirmed known-obsolete identifiers in frozen review dictionary: {known_obsolete_count}",
                f"  Additional identifiers requiring ontology review: {manual_review_count}",
                f"  Significant unresolved identifiers in GP-up: {significant_gp_unresolved}",
                f"  Significant unresolved identifiers in UP-up: {significant_up_unresolved}",
                f"  Significant unresolved identifiers in either direction: {significant_any_unresolved}",
                "  Policy: retain statistical audit, but do not interpret unresolved terms as current GO categories.",
            ])

        lines = lines[:1] + [
            f"Expanded-mask DESeq2 rows: {len(deseq_rows)}",
            f"Genes with valid padj: {len(valid)}",
            f"Official GP-up genes: {len(gp)}",
            f"Official UP-up genes: {len(up)}",
            f"DESeq2 input SHA-256: {sha256(DESEQ)}",
            "Independent Fisher/BH recalculation: completed",
        ] + lines[1:]

    # Report status is evaluated after all tables so a late error marks every
    # row consistently in the compact summary.
    final_status = "PASS" if not errors else "CHECK"
    for row in summary_rows:
        row["status"] = final_status
    write_tsv(
        SUMMARY,
        summary_rows,
        [
            "analysis", "table_path", "tested_categories",
            "annotation_eligible_foreground", "annotation_background",
            "significant_categories_BH_lt_0.05", "sha256", "status",
        ],
    )

    lines.extend(["", f"Warnings: {len(warnings)}"])
    lines.extend(f"WARNING: {warning}" for warning in warnings)
    lines.append(f"Errors detected: {len(errors)}")
    lines.extend(f"ERROR: {error}" for error in errors)
    lines.append(
        "EXPANDED-MASK HOST-FOCUSED PLANT ENRICHMENT INTEGRITY AUDIT: "
        + ("PASS" if not errors else "CHECK REQUIRED")
    )
    AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\nOutput: {GO_REVIEW}")
    print(f"Output: {SUMMARY}")
    print(f"Output: {AUDIT}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
