#!/usr/bin/env python3
"""Read-only scientific preflight for the cross-species discovery scan.

The script validates the registered input links, frozen checksums, core row
counts, identifier relationships, DEG directions, secretome membership,
MEROPS curation, and previously audited candidate sets. It never modifies any
source analysis file. Its only output is a text audit written under 08_logs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from collections import Counter
from pathlib import Path


EXPECTED_INPUTS = {
    "insect_CAP_domain_candidates.tsv",
    "insect_MEROPS_gene_category_curated_v2.tsv",
    "insect_deseq2_all.tsv",
    "insect_eggnog_parsed.tsv",
    "insect_focused_candidate_families.tsv",
    "insect_gene_to_GO.tsv",
    "insect_gene_to_KEGG_ko.tsv",
    "insect_gene_to_KEGG_pathway.tsv",
    "insect_signalp_tmhmm_deseq2_integrated.tsv",
    "plant_MEROPS_gene_category_curated_v2.tsv",
    "plant_MEROPS_peptidase_priority_curated_v2.tsv",
    "plant_counterpart_family_audit.tsv",
    "plant_deseq2_all.tsv",
    "plant_deseq2_eggnog_joined.tsv",
    "plant_developmental_candidates.tsv",
    "plant_representative_protein_map.tsv",
    "plant_signalp_tmhmm_deseq2_integrated.tsv",
}

NA_VALUES = {"", "na", "nan", "none", "null", "."}


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip().lstrip("\ufeff")


def as_float(value: object) -> float | None:
    text = clean(value)
    if text.lower() in NA_VALUES:
        return None
    number = float(text)
    return None if math.isnan(number) else number


def as_bool(value: object) -> bool:
    text = clean(value).lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n", ""}:
        return False
    raise ValueError(f"Unrecognized Boolean value: {value!r}")


def same_number(left: object, right: object) -> bool:
    a = as_float(left)
    b = as_float(right)
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        fields = [clean(name) for name in reader.fieldnames]
        rows: list[dict[str, str]] = []
        for raw in reader:
            row = {clean(key): clean(value) for key, value in raw.items()}
            rows.append(row)
    return fields, rows


def index_unique(
    rows: list[dict[str, str]], key: str, label: str, errors: list[str]
) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    duplicates: list[str] = []
    missing = 0
    for row in rows:
        identifier = clean(row.get(key))
        if not identifier:
            missing += 1
        elif identifier in result:
            duplicates.append(identifier)
        else:
            result[identifier] = row
    if missing:
        errors.append(f"{label}: {missing} row(s) have an empty {key}")
    if duplicates:
        examples = ", ".join(sorted(set(duplicates))[:5])
        errors.append(
            f"{label}: {len(duplicates)} duplicate {key} occurrence(s); "
            f"examples: {examples}"
        )
    return result


def append_counter(report: list[str], counter: Counter[str], indent: str = "  ") -> None:
    if not counter:
        report.append(f"{indent}[none]: 0")
        return
    for key, count in sorted(counter.items(), key=lambda item: (-item[1], item[0])):
        report.append(f"{indent}{key or '[blank]'}: {count}")


def expect_equal(
    report: list[str], errors: list[str], label: str, observed: object, expected: object
) -> None:
    report.append(f"{label}: {observed}")
    if observed != expected:
        errors.append(f"{label}: expected {expected}, observed {observed}")


def require_columns(
    label: str, fields: list[str], required: set[str], errors: list[str]
) -> None:
    missing = sorted(required - set(fields))
    if missing:
        errors.append(f"{label}: missing columns: {', '.join(missing)}")
    duplicates = [name for name, count in Counter(fields).items() if count > 1]
    if duplicates:
        errors.append(f"{label}: duplicate columns: {', '.join(sorted(duplicates))}")


def official_direction(row: dict[str, str]) -> str:
    padj = as_float(row.get("padj"))
    lfc = as_float(row.get("log2FoldChange"))
    if padj is None or lfc is None or padj >= 0.01 or abs(lfc) < 1:
        return "not_official_DEG"
    return "positive_up" if lfc >= 1 else "negative_up"


def compare_numeric_columns(
    label: str,
    left: dict[str, dict[str, str]],
    right: dict[str, dict[str, str]],
    columns: tuple[str, ...],
    errors: list[str],
) -> int:
    mismatches = 0
    examples: list[str] = []
    for identifier in sorted(left.keys() & right.keys()):
        for column in columns:
            if not same_number(left[identifier].get(column), right[identifier].get(column)):
                mismatches += 1
                if len(examples) < 5:
                    examples.append(f"{identifier}:{column}")
    if mismatches:
        errors.append(
            f"{label}: {mismatches} numeric mismatch(es); examples: "
            + ", ".join(examples)
        )
    return mismatches


def scan_mapping(
    path: Path,
    gene_column: str,
    term_column: str,
    allowed_genes: set[str],
    universe_label: str,
    expected_rows: int,
    report: list[str],
    errors: list[str],
) -> None:
    row_count = 0
    genes: set[str] = set()
    empty_gene = 0
    empty_term = 0
    outside: set[str] = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = [clean(name) for name in (reader.fieldnames or [])]
        if gene_column not in fields or term_column not in fields:
            errors.append(
                f"{path.name}: required mapping columns not found; observed {fields}"
            )
            return
        for raw in reader:
            row = {clean(key): clean(value) for key, value in raw.items()}
            row_count += 1
            gene = row.get(gene_column, "")
            term = row.get(term_column, "")
            if not gene:
                empty_gene += 1
            else:
                genes.add(gene)
                if gene not in allowed_genes:
                    outside.add(gene)
            if not term:
                empty_term += 1
    expect_equal(report, errors, f"{path.name} rows", row_count, expected_rows)
    report.append(f"{path.name} unique genes: {len(genes)}")
    report.append(f"{path.name} genes outside {universe_label}: {len(outside)}")
    if empty_gene:
        errors.append(f"{path.name}: {empty_gene} empty gene identifiers")
    if empty_term:
        errors.append(f"{path.name}: {empty_term} empty term identifiers")
    if outside:
        errors.append(
            f"{path.name}: {len(outside)} gene(s) outside {universe_label}"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Discovery-project root (default: parent of 07_scripts)",
    )
    args = parser.parse_args()

    root = args.project_root.expanduser().resolve()
    input_dir = root / "01_input_links"
    control_dir = root / "00_project_control"
    log_dir = root / "08_logs"
    manifest_path = control_dir / "validated_input_manifest.tsv"
    log_path = log_dir / "01_discovery_input_preflight.txt"

    report: list[str] = ["===== CROSS-SPECIES DISCOVERY INPUT PREFLIGHT =====", ""]
    errors: list[str] = []
    warnings: list[str] = []

    for directory in (input_dir, control_dir, log_dir):
        if not directory.is_dir():
            errors.append(f"Missing required directory: {directory}")
    if not manifest_path.is_file():
        errors.append(f"Missing input manifest: {manifest_path}")

    if errors:
        report.extend(f"ERROR: {error}" for error in errors)
        print("\n".join(report))
        return 1

    # ------------------------------------------------------------------
    # Manifest, symlink, size and checksum validation
    # ------------------------------------------------------------------
    report.append("===== 1. VALIDATED INPUT MANIFEST =====")
    manifest_fields, manifest_rows = read_tsv(manifest_path)
    require_columns(
        "validated input manifest",
        manifest_fields,
        {"link_name", "source_path", "size_bytes", "sha256"},
        errors,
    )
    manifest_index = index_unique(
        manifest_rows, "link_name", "validated input manifest", errors
    )
    expect_equal(report, errors, "Manifest records", len(manifest_rows), 17)
    observed_names = set(manifest_index)
    missing_names = sorted(EXPECTED_INPUTS - observed_names)
    extra_names = sorted(observed_names - EXPECTED_INPUTS)
    report.append(f"Missing expected manifest entries: {len(missing_names)}")
    report.append(f"Unexpected manifest entries: {len(extra_names)}")
    if missing_names:
        errors.append("Missing manifest entries: " + ", ".join(missing_names))
    if extra_names:
        errors.append("Unexpected manifest entries: " + ", ".join(extra_names))

    manifest_failures = 0
    for name in sorted(EXPECTED_INPUTS & observed_names):
        row = manifest_index[name]
        link = input_dir / name
        source = Path(row["source_path"]).expanduser()
        if not link.is_symlink():
            errors.append(f"Input is not a symbolic link: {link}")
            manifest_failures += 1
            continue
        if not source.is_file():
            errors.append(f"Manifest source is missing: {source}")
            manifest_failures += 1
            continue
        if not link.exists():
            errors.append(f"Broken input link: {link}")
            manifest_failures += 1
            continue
        if link.resolve() != source.resolve():
            errors.append(f"Link target mismatch for {name}")
            manifest_failures += 1
        try:
            recorded_size = int(row["size_bytes"])
        except ValueError:
            recorded_size = -1
        if source.stat().st_size != recorded_size:
            errors.append(f"Size mismatch for {name}")
            manifest_failures += 1
        if sha256_file(source) != row["sha256"].lower():
            errors.append(f"SHA-256 mismatch for {name}")
            manifest_failures += 1
    report.append(f"Manifest/link/checksum failures: {manifest_failures}")
    report.append("")

    def load(name: str) -> tuple[list[str], list[dict[str, str]]]:
        return read_tsv(input_dir / name)

    # ------------------------------------------------------------------
    # Insect DESeq2 and secretome
    # ------------------------------------------------------------------
    report.append("===== 2. INSECT DESEQ2 AND CLASSICAL SECRETOME =====")
    insect_de_fields, insect_de_rows = load("insect_deseq2_all.tsv")
    require_columns(
        "insect DESeq2",
        insect_de_fields,
        {"Geneid", "baseMean", "log2FoldChange", "stat", "pvalue", "padj"},
        errors,
    )
    insect_de = index_unique(insect_de_rows, "Geneid", "insect DESeq2", errors)
    expect_equal(report, errors, "Insect DESeq2 rows", len(insect_de_rows), 13919)
    valid_insect_padj = sum(as_float(row["padj"]) is not None for row in insect_de_rows)
    insect_directions = Counter(official_direction(row) for row in insect_de_rows)
    report.append(f"Insect genes with valid padj: {valid_insect_padj}")
    expect_equal(
        report,
        errors,
        "Insect official DEGs",
        insect_directions["positive_up"] + insect_directions["negative_up"],
        4406,
    )
    expect_equal(report, errors, "Insect pupa-up DEGs", insect_directions["positive_up"], 2497)
    expect_equal(report, errors, "Insect larva-up DEGs", insect_directions["negative_up"], 1909)

    insect_int_fields, insect_int_rows = load(
        "insect_signalp_tmhmm_deseq2_integrated.tsv"
    )
    require_columns(
        "insect SignalP/DeepTMHMM/DESeq2 integration",
        insect_int_fields,
        {
            "gene_id",
            "analysis_class",
            "primary_classical_secreted",
            "baseMean",
            "log2FoldChange",
            "stat",
            "pvalue",
            "padj",
            "deseq2_status",
            "official_DEG",
        },
        errors,
    )
    insect_int = index_unique(
        insect_int_rows,
        "gene_id",
        "insect SignalP/DeepTMHMM/DESeq2 integration",
        errors,
    )
    expect_equal(report, errors, "Insect integrated protein rows", len(insect_int_rows), 13455)
    insect_missing_protein = set(insect_de) - set(insect_int)
    insect_extra_protein = set(insect_int) - set(insect_de)
    expect_equal(
        report,
        errors,
        "Insect DESeq2 genes without representative protein",
        len(insect_missing_protein),
        1068,
    )
    expect_equal(
        report,
        errors,
        "Insect representative proteins outside DESeq2 table",
        len(insect_extra_protein),
        604,
    )
    expect_equal(
        report,
        errors,
        "Insect DESeq2/representative-protein intersection",
        len(set(insect_de) & set(insect_int)),
        12851,
    )
    extra_insect_status = Counter(
        insect_int[gene]["deseq2_status"] for gene in insect_extra_protein
    )
    report.append("Insect representative proteins outside DESeq2 by status:")
    append_counter(report, extra_insect_status)
    expect_equal(
        report,
        errors,
        "Insect outside-DESeq2 proteins classified as excluded_low_count",
        extra_insect_status["excluded_low_count"],
        604,
    )
    expect_equal(
        report,
        errors,
        "Insect outside-DESeq2 proteins with any other status",
        sum(
            count
            for status, count in extra_insect_status.items()
            if status != "excluded_low_count"
        ),
        0,
    )
    insect_numeric_mismatches = compare_numeric_columns(
        "Insect DESeq2-to-integrated comparison",
        insect_de,
        insect_int,
        ("baseMean", "log2FoldChange", "stat", "pvalue", "padj"),
        errors,
    )
    report.append(f"Insect DESeq2/integrated numeric mismatches: {insect_numeric_mismatches}")

    insect_classes = Counter(row["analysis_class"] for row in insect_int_rows)
    report.append("Insect final analysis classes:")
    append_counter(report, insect_classes)
    expect_equal(
        report,
        errors,
        "Insect classical-secretome proteins",
        insect_classes["classical_secreted_core"],
        2039,
    )
    insect_primary_mismatch = sum(
        (row["analysis_class"] == "classical_secreted_core")
        != as_bool(row["primary_classical_secreted"])
        for row in insect_int_rows
    )
    expect_equal(
        report,
        errors,
        "Insect analysis-class/primary-secretome mismatches",
        insect_primary_mismatch,
        0,
    )

    insect_status_mismatch = 0
    insect_official_mismatch = 0
    for row in insect_int_rows:
        direction = official_direction(row)
        status = row["deseq2_status"]
        if direction == "positive_up" and status != "pupa_up":
            insect_status_mismatch += 1
        elif direction == "negative_up" and status != "larva_up":
            insect_status_mismatch += 1
        if as_bool(row["official_DEG"]) != (direction != "not_official_DEG"):
            insect_official_mismatch += 1
    expect_equal(
        report,
        errors,
        "Insect official-direction/status mismatches",
        insect_status_mismatch,
        0,
    )
    expect_equal(
        report,
        errors,
        "Insect official-DEG flag mismatches",
        insect_official_mismatch,
        0,
    )

    insect_classical = [
        row for row in insect_int_rows if row["analysis_class"] == "classical_secreted_core"
    ]
    insect_classical_status = Counter(row["deseq2_status"] for row in insect_classical)
    report.append("Insect classical-secretome DESeq2 statuses:")
    append_counter(report, insect_classical_status)
    expect_equal(
        report,
        errors,
        "Eligible larva-up classical-secretome cohort",
        insect_classical_status["larva_up"],
        544,
    )
    expect_equal(
        report,
        errors,
        "Pupa-up classical-secretome context",
        insect_classical_status["pupa_up"],
        597,
    )
    report.append("")

    # ------------------------------------------------------------------
    # Insect annotations and focused candidate families
    # ------------------------------------------------------------------
    report.append("===== 3. INSECT ANNOTATIONS AND PRIOR FOCUSED SETS =====")
    insect_egg_fields, insect_egg_rows = load("insect_eggnog_parsed.tsv")
    require_columns(
        "insect eggNOG",
        insect_egg_fields,
        {
            "Geneid",
            "seed_ortholog",
            "Description",
            "Preferred_name",
            "GOs",
            "KEGG_ko",
            "KEGG_Pathway",
            "CAZy",
            "PFAMs",
        },
        errors,
    )
    insect_egg = index_unique(insect_egg_rows, "Geneid", "insect eggNOG", errors)
    expect_equal(report, errors, "Insect representative eggNOG rows", len(insect_egg_rows), 12352)
    expect_equal(
        report,
        errors,
        "Insect eggNOG genes outside representative-protein universe",
        len(set(insect_egg) - set(insect_int)),
        0,
    )

    scan_mapping(
        input_dir / "insect_gene_to_GO.tsv",
        "Geneid",
        "GO",
        set(insect_int),
        "representative-protein universe",
        1180785,
        report,
        errors,
    )
    scan_mapping(
        input_dir / "insect_gene_to_KEGG_ko.tsv",
        "Geneid",
        "KEGG_ko",
        set(insect_int),
        "representative-protein universe",
        11086,
        report,
        errors,
    )
    scan_mapping(
        input_dir / "insect_gene_to_KEGG_pathway.tsv",
        "Geneid",
        "KEGG_Pathway",
        set(insect_int),
        "representative-protein universe",
        45900,
        report,
        errors,
    )

    focused_fields, focused_rows = load("insect_focused_candidate_families.tsv")
    require_columns(
        "insect focused candidates",
        focused_fields,
        {"gene_id", "candidate_group", "analysis_class", "deseq2_status"},
        errors,
    )
    focused_index = index_unique(
        focused_rows, "gene_id", "insect focused candidates", errors
    )
    expect_equal(report, errors, "Insect focused-family candidate genes", len(focused_rows), 50)
    focused_groups = Counter(row["candidate_group"] for row in focused_rows)
    report.append("Insect focused candidates by group:")
    append_counter(report, focused_groups)
    if sorted(focused_groups.values()) != [9, 13, 13, 15]:
        errors.append(
            "Insect focused-family distribution: expected group sizes 9, 13, 13, 15; "
            f"observed {sorted(focused_groups.values())}"
        )
    expect_equal(
        report,
        errors,
        "Insect focused candidates outside representative-protein universe",
        len(set(focused_index) - set(insect_int)),
        0,
    )

    cap_fields, cap_rows = load("insect_CAP_domain_candidates.tsv")
    require_columns(
        "insect CAP candidates",
        cap_fields,
        {"gene_id", "CAP_evidence", "PFAM_CAP", "CAP_annotation", "analysis_class"},
        errors,
    )
    cap_index = index_unique(cap_rows, "gene_id", "insect CAP candidates", errors)
    expect_equal(report, errors, "Insect CAP candidate genes", len(cap_rows), 16)
    expect_equal(
        report,
        errors,
        "Insect CAP candidates outside representative-protein universe",
        len(set(cap_index) - set(insect_int)),
        0,
    )
    report.append("")

    # ------------------------------------------------------------------
    # Insect MEROPS
    # ------------------------------------------------------------------
    report.append("===== 4. INSECT MEROPS CURATED V2 =====")
    insect_merops_fields, insect_merops_rows = load(
        "insect_MEROPS_gene_category_curated_v2.tsv"
    )
    require_columns(
        "insect MEROPS",
        insect_merops_fields,
        {
            "gene_id",
            "merops_category",
            "merops_family",
            "evidence_tier",
            "passes_corrected_filter",
            "curated_AB_candidate",
            "DESeq2_status",
            "analysis_class",
            "classical_secreted_core",
        },
        errors,
    )
    expect_equal(report, errors, "Insect corrected-v2 gene-category rows", len(insect_merops_rows), 1870)
    insect_merops_keys = [
        (row["gene_id"], row["merops_category"]) for row in insect_merops_rows
    ]
    expect_equal(
        report,
        errors,
        "Duplicate insect gene-by-MEROPS-category rows",
        len(insect_merops_keys) - len(set(insect_merops_keys)),
        0,
    )
    expect_equal(
        report,
        errors,
        "Insect MEROPS rows failing corrected filter",
        sum(not as_bool(row["passes_corrected_filter"]) for row in insect_merops_rows),
        0,
    )
    insect_merops_categories = Counter(row["merops_category"] for row in insect_merops_rows)
    report.append("Insect corrected-v2 MEROPS rows by category:")
    append_counter(report, insect_merops_categories)
    insect_secreted_inhibitors = [
        row
        for row in insect_merops_rows
        if "inhibitor" in row["merops_category"].lower()
        and row["analysis_class"] == "classical_secreted_core"
        and as_bool(row["curated_AB_candidate"])
    ]
    insect_secreted_inhibitor_genes = {row["gene_id"] for row in insect_secreted_inhibitors}
    expect_equal(
        report,
        errors,
        "Insect curated A+B classical-secreted inhibitor rows",
        len(insect_secreted_inhibitors),
        149,
    )
    expect_equal(
        report,
        errors,
        "Insect curated A+B classical-secreted inhibitor genes",
        len(insect_secreted_inhibitor_genes),
        149,
    )
    insect_inhibitor_tiers = Counter(row["evidence_tier"] for row in insect_secreted_inhibitors)
    report.append("Insect curated secreted inhibitors by evidence tier:")
    append_counter(report, insect_inhibitor_tiers)
    expect_equal(
        report,
        errors,
        "Insect inhibitor Tier A genes",
        sum(count for tier, count in insect_inhibitor_tiers.items() if tier.startswith("A_")),
        64,
    )
    expect_equal(
        report,
        errors,
        "Insect inhibitor Tier B genes",
        sum(count for tier, count in insect_inhibitor_tiers.items() if tier.startswith("B_")),
        85,
    )
    insect_inhibitor_status = Counter(
        row["DESeq2_status"] for row in insect_secreted_inhibitors
    )
    report.append("Insect curated secreted inhibitors by DESeq2 status:")
    append_counter(report, insect_inhibitor_status)
    for status, expected in {
        "pupa_up": 68,
        "larva_up": 23,
        "non_DE": 57,
        "excluded_low_count": 1,
        "NA_padj_not_testable": 0,
    }.items():
        expect_equal(
            report,
            errors,
            f"Insect inhibitor {status}",
            insect_inhibitor_status[status],
            expected,
        )
    expect_equal(
        report,
        errors,
        "Insect MEROPS genes outside representative-protein universe",
        len({row["gene_id"] for row in insect_merops_rows} - set(insect_int)),
        0,
    )
    report.append("")

    # ------------------------------------------------------------------
    # Plant DESeq2, annotation and protein universe
    # ------------------------------------------------------------------
    report.append("===== 5. PLANT DESEQ2, ANNOTATION AND PROTEIN UNIVERSE =====")
    plant_de_fields, plant_de_rows = load("plant_deseq2_all.tsv")
    require_columns(
        "plant DESeq2",
        plant_de_fields,
        {
            "gene_id",
            "baseMean",
            "log2FoldChange",
            "stat",
            "pvalue",
            "padj",
            "DESeq2_status",
        },
        errors,
    )
    plant_de = index_unique(plant_de_rows, "gene_id", "plant DESeq2", errors)
    expect_equal(report, errors, "Plant DESeq2 rows", len(plant_de_rows), 28189)
    valid_plant_padj = sum(as_float(row["padj"]) is not None for row in plant_de_rows)
    expect_equal(report, errors, "Plant genes with valid padj", valid_plant_padj, 27948)
    expect_equal(report, errors, "Plant genes with NA padj", len(plant_de_rows) - valid_plant_padj, 241)
    plant_directions = Counter(official_direction(row) for row in plant_de_rows)
    expect_equal(
        report,
        errors,
        "Plant official DEGs",
        plant_directions["positive_up"] + plant_directions["negative_up"],
        7060,
    )
    expect_equal(report, errors, "Plant GP-up DEGs", plant_directions["positive_up"], 5593)
    expect_equal(report, errors, "Plant UP-up DEGs", plant_directions["negative_up"], 1467)

    plant_direction_status_mismatch = 0
    for row in plant_de_rows:
        direction = official_direction(row)
        status = row["DESeq2_status"]
        if direction == "positive_up" and status != "GP_up":
            plant_direction_status_mismatch += 1
        elif direction == "negative_up" and status != "UP_up":
            plant_direction_status_mismatch += 1
    expect_equal(
        report,
        errors,
        "Plant official-direction/status mismatches",
        plant_direction_status_mismatch,
        0,
    )
    report.append("Verified plant contrast direction: positive log2FoldChange = GP-up")

    plant_join_fields, plant_join_rows = load("plant_deseq2_eggnog_joined.tsv")
    require_columns(
        "plant DESeq2/eggNOG join",
        plant_join_fields,
        {
            "gene_id",
            "baseMean",
            "log2FoldChange",
            "stat",
            "pvalue",
            "padj",
            "DESeq2_status",
            "has_eggnog_annotation",
            "seed_ortholog",
            "Description",
            "Preferred_name",
            "GOs",
            "KEGG_ko",
            "KEGG_Pathway",
            "CAZy",
            "PFAMs",
        },
        errors,
    )
    plant_join = index_unique(
        plant_join_rows, "gene_id", "plant DESeq2/eggNOG join", errors
    )
    expect_equal(report, errors, "Plant DESeq2/eggNOG joined rows", len(plant_join_rows), 28189)
    expect_equal(
        report,
        errors,
        "Plant DESeq2/joined gene-set differences",
        len(set(plant_de) ^ set(plant_join)),
        0,
    )
    plant_numeric_mismatches = compare_numeric_columns(
        "Plant DESeq2-to-eggNOG comparison",
        plant_de,
        plant_join,
        ("baseMean", "log2FoldChange", "stat", "pvalue", "padj"),
        errors,
    )
    report.append(f"Plant DESeq2/joined numeric mismatches: {plant_numeric_mismatches}")
    plant_annotated = {
        row["gene_id"] for row in plant_join_rows if as_bool(row["has_eggnog_annotation"])
    }
    plant_deg_ids = {
        gene for gene, row in plant_de.items() if official_direction(row) != "not_official_DEG"
    }
    expect_equal(report, errors, "Plant eggNOG-annotated tested genes", len(plant_annotated), 18904)
    expect_equal(
        report,
        errors,
        "Plant eggNOG-annotated official DEGs",
        len(plant_annotated & plant_deg_ids),
        2414,
    )
    expect_equal(
        report,
        errors,
        "Plant genes without eggNOG annotation",
        len(plant_join_rows) - len(plant_annotated),
        9285,
    )
    expect_equal(
        report,
        errors,
        "Plant official DEGs without eggNOG annotation",
        len(plant_deg_ids - plant_annotated),
        4646,
    )

    plant_map_fields, plant_map_rows = load("plant_representative_protein_map.tsv")
    require_columns(
        "plant representative-protein map",
        plant_map_fields,
        {"protein_id", "transcript_id", "gene_id", "aa_length"},
        errors,
    )
    plant_map = index_unique(
        plant_map_rows, "gene_id", "plant representative-protein map", errors
    )
    expect_equal(report, errors, "Plant representative-protein rows", len(plant_map_rows), 41330)
    protein_ids = [row["protein_id"] for row in plant_map_rows]
    expect_equal(
        report,
        errors,
        "Duplicate plant representative protein IDs",
        len(protein_ids) - len(set(protein_ids)),
        0,
    )

    plant_int_fields, plant_int_rows = load(
        "plant_signalp_tmhmm_deseq2_integrated.tsv"
    )
    require_columns(
        "plant SignalP/DeepTMHMM/DESeq2 integration",
        plant_int_fields,
        {
            "gene_id",
            "analysis_class",
            "primary_classical_secretome",
            "baseMean",
            "log2FoldChange",
            "stat",
            "pvalue",
            "padj",
            "DESeq2_status",
        },
        errors,
    )
    plant_int = index_unique(
        plant_int_rows,
        "gene_id",
        "plant SignalP/DeepTMHMM/DESeq2 integration",
        errors,
    )
    expect_equal(report, errors, "Plant integrated protein rows", len(plant_int_rows), 41330)
    expect_equal(
        report,
        errors,
        "Plant protein-map/integrated gene-set differences",
        len(set(plant_map) ^ set(plant_int)),
        0,
    )
    plant_classes = Counter(row["analysis_class"] for row in plant_int_rows)
    report.append("Plant final analysis classes:")
    append_counter(report, plant_classes)
    expect_equal(
        report,
        errors,
        "Plant classical-secretome proteins",
        plant_classes["classical_secreted_core"],
        2171,
    )
    plant_primary_mismatch = sum(
        (row["analysis_class"] == "classical_secreted_core")
        != as_bool(row["primary_classical_secretome"])
        for row in plant_int_rows
    )
    expect_equal(
        report,
        errors,
        "Plant analysis-class/primary-secretome mismatches",
        plant_primary_mismatch,
        0,
    )
    plant_status = Counter(row["DESeq2_status"] for row in plant_int_rows)
    report.append("Plant representative-protein DESeq2 statuses:")
    append_counter(report, plant_status)
    for status, expected in {
        "GP_up": 4926,
        "UP_up": 1083,
        "non_DE": 17006,
        "NA_padj_not_testable": 224,
        "excluded_low_count": 18091,
    }.items():
        expect_equal(report, errors, f"Plant protein status {status}", plant_status[status], expected)
    plant_gp_ids = {
        gene for gene, row in plant_de.items() if official_direction(row) == "positive_up"
    }
    plant_gp_with_protein = plant_gp_ids & set(plant_int)
    expect_equal(
        report,
        errors,
        "Plant GP-up genes with representative proteins",
        len(plant_gp_with_protein),
        4926,
    )
    expect_equal(
        report,
        errors,
        "Plant GP-up genes without representative proteins",
        len(plant_gp_ids - set(plant_int)),
        667,
    )
    plant_gp_classical = [
        row
        for row in plant_int_rows
        if row["DESeq2_status"] == "GP_up"
        and row["analysis_class"] == "classical_secreted_core"
    ]
    expect_equal(
        report,
        errors,
        "Plant GP-up classical-secretome genes previously screened by MEROPS",
        len(plant_gp_classical),
        549,
    )
    report.append("")

    # ------------------------------------------------------------------
    # Plant MEROPS
    # ------------------------------------------------------------------
    report.append("===== 6. PLANT MEROPS CURATED V2 AND PRIORITY SET =====")
    plant_merops_fields, plant_merops_rows = load(
        "plant_MEROPS_gene_category_curated_v2.tsv"
    )
    require_columns(
        "plant MEROPS gene-category table",
        plant_merops_fields,
        {
            "gene_id",
            "merops_category",
            "evidence_tier",
            "passes_corrected_filter",
            "curated_AB_candidate",
            "DESeq2_status",
            "analysis_class",
        },
        errors,
    )
    expect_equal(report, errors, "Plant corrected-v2 gene-category rows", len(plant_merops_rows), 138)
    plant_merops_keys = [
        (row["gene_id"], row["merops_category"]) for row in plant_merops_rows
    ]
    expect_equal(
        report,
        errors,
        "Duplicate plant gene-by-MEROPS-category rows",
        len(plant_merops_keys) - len(set(plant_merops_keys)),
        0,
    )
    expect_equal(
        report,
        errors,
        "Plant MEROPS rows failing corrected filter",
        sum(not as_bool(row["passes_corrected_filter"]) for row in plant_merops_rows),
        0,
    )
    plant_merops_categories = Counter(row["merops_category"] for row in plant_merops_rows)
    report.append("Plant corrected-v2 MEROPS rows by category:")
    append_counter(report, plant_merops_categories)
    plant_peptidases = [
        row for row in plant_merops_rows if "peptidase" in row["merops_category"].lower()
    ]
    plant_peptidase_ids = {row["gene_id"] for row in plant_peptidases}
    expect_equal(report, errors, "Plant broad peptidase-homologue rows", len(plant_peptidases), 98)
    expect_equal(report, errors, "Plant broad peptidase-homologue genes", len(plant_peptidase_ids), 98)
    plant_peptidase_tiers = Counter(row["evidence_tier"] for row in plant_peptidases)
    report.append("Plant peptidase homologues by evidence tier:")
    append_counter(report, plant_peptidase_tiers)
    for prefix, expected in {"A_": 82, "B_": 11, "C_": 5}.items():
        expect_equal(
            report,
            errors,
            f"Plant peptidase Tier {prefix[0]} genes",
            sum(count for tier, count in plant_peptidase_tiers.items() if tier.startswith(prefix)),
            expected,
        )
    expect_equal(
        report,
        errors,
        "Plant curated A+B peptidase genes",
        len({row["gene_id"] for row in plant_peptidases if as_bool(row["curated_AB_candidate"])}),
        93,
    )

    priority_fields, priority_rows = load(
        "plant_MEROPS_peptidase_priority_curated_v2.tsv"
    )
    require_columns(
        "plant MEROPS priority table",
        priority_fields,
        {
            "gene_id",
            "evidence_tier",
            "curated_AB_candidate",
            "manuscript_priority",
            "DESeq2_status",
            "analysis_class",
        },
        errors,
    )
    priority_index = index_unique(
        priority_rows, "gene_id", "plant MEROPS priority table", errors
    )
    expect_equal(report, errors, "Plant MEROPS priority-table rows", len(priority_rows), 98)
    expect_equal(
        report,
        errors,
        "Plant peptidase/priority-table gene-set differences",
        len(plant_peptidase_ids ^ set(priority_index)),
        0,
    )
    priority_counts = Counter(row["manuscript_priority"] for row in priority_rows)
    report.append("Plant manuscript-priority distribution:")
    append_counter(report, priority_counts)
    priority_1 = sum(
        count for label, count in priority_counts.items() if label.startswith("priority_1")
    )
    priority_2 = sum(
        count for label, count in priority_counts.items() if label.startswith("priority_2")
    )
    expect_equal(report, errors, "Plant Priority 1 candidates", priority_1, 19)
    expect_equal(report, errors, "Plant Priority 2 candidates", priority_2, 1)
    expect_equal(report, errors, "Plant total manuscript-priority candidates", priority_1 + priority_2, 20)
    priority_gene_ids = {
        row["gene_id"]
        for row in priority_rows
        if row["manuscript_priority"].startswith(("priority_1", "priority_2"))
    }
    expect_equal(
        report,
        errors,
        "Plant priority genes outside GP-up DEG set",
        len(priority_gene_ids - plant_gp_ids),
        0,
    )
    expect_equal(
        report,
        errors,
        "Plant priority genes outside classical secretome",
        len(
            priority_gene_ids
            - {
                row["gene_id"]
                for row in plant_int_rows
                if row["analysis_class"] == "classical_secreted_core"
            }
        ),
        0,
    )
    report.append(
        "Plant MEROPS scope boundary: only the 549 GP-up classical-secretome genes "
        "were screened; other plant genes are unassessed, not MEROPS-negative."
    )
    report.append("")

    # ------------------------------------------------------------------
    # Plant hypothesis-guided audit tables
    # ------------------------------------------------------------------
    report.append("===== 7. PLANT PRIOR HYPOTHESIS-GUIDED AUDITS =====")
    counterpart_fields, counterpart_rows = load("plant_counterpart_family_audit.tsv")
    require_columns(
        "plant counterpart-family audit",
        counterpart_fields,
        {"candidate_group", "gene_id", "family_evidence", "DESeq2_status", "analysis_class"},
        errors,
    )
    expect_equal(report, errors, "Plant counterpart-family audit rows", len(counterpart_rows), 244)
    counterpart_keys = [(row["candidate_group"], row["gene_id"]) for row in counterpart_rows]
    expect_equal(
        report,
        errors,
        "Duplicate plant counterpart group-by-gene rows",
        len(counterpart_keys) - len(set(counterpart_keys)),
        0,
    )
    expect_equal(
        report,
        errors,
        "Plant counterpart-audit genes outside representative-protein universe",
        len({row["gene_id"] for row in counterpart_rows} - set(plant_int)),
        0,
    )
    counterpart_groups = Counter(row["candidate_group"] for row in counterpart_rows)
    report.append("Plant counterpart-audit rows by candidate group:")
    append_counter(report, counterpart_groups)

    development_fields, development_rows = load("plant_developmental_candidates.tsv")
    require_columns(
        "plant developmental candidates",
        development_fields,
        {"gene_id", "DESeq2_status", "target_groups", "target_regulators"},
        errors,
    )
    development_index = index_unique(
        development_rows, "gene_id", "plant developmental candidates", errors
    )
    expect_equal(report, errors, "Plant developmental candidate genes", len(development_rows), 8)
    expect_equal(
        report,
        errors,
        "Plant developmental candidates outside representative-protein universe",
        len(set(development_index) - set(plant_int)),
        0,
    )
    report.append("")

    # ------------------------------------------------------------------
    # Eligibility and stop/go decision
    # ------------------------------------------------------------------
    report.append("===== 8. DISCOVERY ELIGIBILITY FREEZE =====")
    expect_equal(
        report,
        errors,
        "Insect eligible cohort for discovery ranking",
        len(
            {
                row["gene_id"]
                for row in insect_int_rows
                if row["analysis_class"] == "classical_secreted_core"
                and row["deseq2_status"] == "larva_up"
            }
        ),
        544,
    )
    expect_equal(
        report,
        errors,
        "Plant eligible cohort for discovery ranking",
        len(plant_gp_ids),
        5593,
    )
    report.append("Ranking inputs frozen: padj, absolute log2FoldChange and baseMean")
    report.append("No ranking was performed by this preflight script.")
    report.append("")

    report.append(f"Warnings: {len(warnings)}")
    for warning in warnings:
        report.append(f"WARNING: {warning}")
    report.append(f"Errors detected: {len(errors)}")
    for error in errors:
        report.append(f"ERROR: {error}")
    if errors:
        report.append("DISCOVERY INPUT PREFLIGHT: CHECK REQUIRED")
    else:
        report.append("DISCOVERY INPUT PREFLIGHT: PASS")

    text = "\n".join(report) + "\n"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path.write_text(text, encoding="utf-8")
    print(text, end="")
    print(f"Audit: {log_path}")
    return 1 if errors else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        raise
