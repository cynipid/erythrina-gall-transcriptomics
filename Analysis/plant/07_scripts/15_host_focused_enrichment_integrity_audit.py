#!/usr/bin/env python3
"""Integrity audit for the host-focused GP-versus-UP GO/KEGG enrichment branch.

This script is read-only with respect to scientific result tables.  It verifies
the host-focused DESeq2 foreground sizes and the structure/accounting of the
four GO/KEGG enrichment outputs produced by 14_host_focused_GO_KEGG_enrichment.R.
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
ENRICHMENT = ROOT / "05_host_focused_DESeq2_sensitivity" / "enrichment"
LOGS = ROOT / "05_host_focused_DESeq2_sensitivity" / "logs"
CONTROL = ROOT / "00_project_control"
DESEQ = RESULTS / "plant_gene_DESeq2_GP_vs_UP_host_focused_eggnog_joined.tsv"

# These are the frozen, observed values from the completed host-focused branch.
EXPECTED = {
    ("GO", "GP"): {"foreground": 854, "background": 11439, "significant": 67},
    ("GO", "UP"): {"foreground": 582, "background": 11439, "significant": 30},
    ("KEGG", "GP"): {"foreground": 441, "background": 5993, "significant": 6},
    ("KEGG", "UP"): {"foreground": 245, "background": 5993, "significant": 4},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fnum(value: str) -> float:
    if value in {"", "NA", "NaN"}:
        return math.nan
    return float(value)


def resolve_enrichment_tables() -> dict[tuple[str, str], Path]:
    """Find output tables by contents, not an assumed output filename.

    The enrichment script is allowed to use any informative filename.  A table
    is identified from its standard ORA columns, GO versus KEGG term syntax,
    and its annotation-eligible foreground size (which is unique here).
    """
    required = {
        "term_id", "foreground_genes_with_term", "foreground_gene_set_size",
        "background_genes_with_term", "background_gene_set_size", "pvalue", "padj",
    }
    resolved: dict[tuple[str, str], list[Path]] = {key: [] for key in EXPECTED}
    for root in (ENRICHMENT, RESULTS):
        for path in root.rglob("*.tsv"):
            try:
                rows = read_tsv(path)
            except (UnicodeDecodeError, csv.Error):
                continue
            if not rows or not required.issubset(rows[0]):
                continue
            term_ids = [r.get("term_id", "") for r in rows]
            kind = "GO" if all(t.startswith("GO:") for t in term_ids) else "KEGG" if all(t.startswith("map") for t in term_ids) else None
            if kind is None:
                continue
            fg_sizes = {r.get("foreground_gene_set_size", "") for r in rows}
            bg_sizes = {r.get("background_gene_set_size", "") for r in rows}
            for (candidate_kind, direction), expected in EXPECTED.items():
                if (kind == candidate_kind and fg_sizes == {str(expected["foreground"])}
                        and bg_sizes == {str(expected["background"])}):
                    resolved[(candidate_kind, direction)].append(path)
    result: dict[tuple[str, str], Path] = {}
    for key, paths in resolved.items():
        if len(paths) == 1:
            result[key] = paths[0]
    return result


def main() -> int:
    LOGS.mkdir(parents=True, exist_ok=True)
    CONTROL.mkdir(parents=True, exist_ok=True)
    lines = ["===== HOST-FOCUSED PLANT ENRICHMENT INTEGRITY AUDIT ====="]
    errors: list[str] = []
    records: list[dict[str, object]] = []

    if not DESEQ.is_file():
        errors.append(f"Missing host-focused annotated DESeq2 input: {DESEQ}")
    else:
        deseq = read_tsv(DESEQ)
        if len(deseq) != 25779:
            errors.append(f"Host-focused DESeq2 row count is {len(deseq)}; expected 25779")
        gp = sum(fnum(r.get("padj", "NA")) < 0.01 and fnum(r.get("log2FoldChange", "NA")) >= 1 for r in deseq)
        up = sum(fnum(r.get("padj", "NA")) < 0.01 and fnum(r.get("log2FoldChange", "NA")) <= -1 for r in deseq)
        lines += [
            f"Host-focused DESeq2 rows: {len(deseq)}",
            f"Official GP-up genes: {gp}",
            f"Official UP-up genes: {up}",
            f"Input SHA-256: {sha256(DESEQ)}",
        ]
        if gp != 3179:
            errors.append(f"Official GP-up count is {gp}; expected 3179")
        if up != 1396:
            errors.append(f"Official UP-up count is {up}; expected 1396")

    table_paths = resolve_enrichment_tables()
    required_cols = {
        "term_id", "foreground_genes_with_term", "foreground_gene_set_size",
        "background_genes_with_term", "background_gene_set_size", "pvalue", "padj",
    }
    for (kind, direction), expected in EXPECTED.items():
        table = table_paths.get((kind, direction))
        label = f"{direction}_{kind}"
        if table is None:
            errors.append(
                f"Could not uniquely identify {label} enrichment table by its ORA contents in {ENRICHMENT}"
            )
            records.append({"table": label, "path": "NOT_FOUND", "tested": "NA", "significant": "NA", "status": "FAIL"})
            continue
        rows = read_tsv(table)
        missing = required_cols - set(rows[0]) if rows else required_cols
        if missing:
            errors.append(f"{table.name} lacks required columns: {', '.join(sorted(missing))}")
        term_ids = [r.get("term_id", "") for r in rows]
        if len(term_ids) != len(set(term_ids)):
            errors.append(f"{table.name} has duplicate term IDs")
        fg_sizes = {r.get("foreground_gene_set_size", "") for r in rows}
        bg_sizes = {r.get("background_gene_set_size", "") for r in rows}
        if fg_sizes != {str(expected["foreground"])}:
            errors.append(f"{table.name} foreground size {sorted(fg_sizes)}; expected {expected['foreground']}")
        if bg_sizes != {str(expected["background"])}:
            errors.append(f"{table.name} background size {sorted(bg_sizes)}; expected {expected['background']}")
        for r in rows:
            for key in ("pvalue", "padj"):
                value = fnum(r.get(key, "NA"))
                if not (0 <= value <= 1):
                    errors.append(f"{table.name} has invalid {key}: {r.get(key)}")
                    break
            if fnum(r.get("foreground_genes_with_term", "NA")) > fnum(r.get("foreground_gene_set_size", "NA")):
                errors.append(f"{table.name} contains a foreground count larger than its foreground set")
                break
            if fnum(r.get("background_genes_with_term", "NA")) > fnum(r.get("background_gene_set_size", "NA")):
                errors.append(f"{table.name} contains a background count larger than its background set")
                break
        sig = sum(fnum(r.get("padj", "NA")) < 0.05 for r in rows)
        if sig != expected["significant"]:
            errors.append(f"{table.name} significant-category count is {sig}; expected {expected['significant']}")
        lines += [
            "",
            f"{label}:",
            f"  Table: {table}",
            f"  Tested categories: {len(rows)}",
            f"  Annotation-eligible foreground genes: {next(iter(fg_sizes), 'NA')}",
            f"  Annotation background genes: {next(iter(bg_sizes), 'NA')}",
            f"  Significant categories (BH < 0.05): {sig}",
            f"  SHA-256: {sha256(table)}",
        ]
        records.append({"table": label, "path": str(table), "tested": len(rows), "significant": sig, "status": "PASS"})

    lines += ["", f"Errors detected: {len(errors)}"]
    if errors:
        lines.extend(f"ERROR: {e}" for e in errors)
        lines.append("HOST-FOCUSED PLANT ENRICHMENT INTEGRITY AUDIT: CHECK REQUIRED")
    else:
        lines.append("HOST-FOCUSED PLANT ENRICHMENT INTEGRITY AUDIT: PASS")
    audit = LOGS / "15_host_focused_enrichment_integrity_audit.txt"
    audit.write_text("\n".join(lines) + "\n", encoding="utf-8")
    summary = CONTROL / "15_host_focused_enrichment_integrity_summary.tsv"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["table", "path", "tested", "significant", "status"], delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print("\n".join(lines))
    print(f"\nOutput: {audit}")
    print(f"Output: {summary}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
