#!/usr/bin/env python3
"""Reparse and curate the existing insect and plant MEROPS DIAMOND screens.

This script does not rerun DIAMOND and does not overwrite the original MEROPS
outputs. It corrects single-HSP MEROPS-unit coverage, distinguishes named,
unassigned, and explicitly non-functional MEROPS reference records, and writes
v2 candidate tables plus membership comparisons with the original sets.

Run on Linux/WSL:
    python3 ~/projects/qe_merops_screen/scripts/07_merops_reparse_curation_v2.py
"""

from __future__ import annotations

import csv
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


# ---------------------------------------------------------------------------
# 1. Frozen project paths and thresholds
# ---------------------------------------------------------------------------

USER_HOME = Path.home()
BASE = USER_HOME / "projects" / "qe_merops_screen"
RESULTS = BASE / "results"

INSECT_RAW_DIR = BASE / "03_insect"
PLANT_RAW_DIR = BASE / "04_plant"

INSECT_INTEGRATED = (
    USER_HOME
    / "projects/qe_genome_guided/11_signalp_tmhmm/results"
    / "signalp6_deeptmhmm_deseq2_primary_13455.tsv"
)
PLANT_SELECTED = (
    PLANT_RAW_DIR / "Evariegata_GP_up_classical_secreted_proteins.tsv"
)
PLANT_EGGNOG = (
    USER_HOME
    / "projects/qe_plant_denovo/02_eggnog/processed"
    / "E_variegata_longest_protein_per_gene.emapper.annotations"
)

PROVENANCE = BASE / "05_provenance" / "v1_membership_baselines"

OLD_INSECT_CANDIDATES = (
    PROVENANCE / "insect_MEROPS_classical_secreted_inhibitor_candidates.tsv"
)
OLD_PLANT_CANDIDATES = (
    PROVENANCE / "plant_GPup_secretome_MEROPS_peptidase_candidates.tsv"
)

EXPECTED_INSECT_RAW_ROWS = 37204
EXPECTED_INSECT_RAW_HIT_GENES = 1920
EXPECTED_INSECT_INTEGRATED_ROWS = 13455
EXPECTED_PLANT_RAW_ROWS = 3288
EXPECTED_PLANT_RAW_HIT_GENES = 139
EXPECTED_PLANT_SELECTED_GENES = 549

MIN_IDENTITY = 25.0
MIN_SUBJECT_SPAN_COVERAGE = 40.0
MAX_EVALUE = 1e-5

RAW_COLUMNS = [
    "gene_id",
    "merops_reference_id",
    "percent_identity",
    "alignment_length",
    "query_length",
    "subject_length",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bitscore",
    "diamond_query_coverage_percent",
    "merops_title",
]

PEPTIDASE_PREFIXES = set("ACGMNPSTU")
VALID_DESEQ2_STATUSES = {
    "pupa_up",
    "larva_up",
    "GP_up",
    "UP_up",
    "non_DE",
    "NA_padj_not_testable",
    "excluded_low_count",
    "absent_from_quantification",
}


# ---------------------------------------------------------------------------
# 2. General helpers
# ---------------------------------------------------------------------------

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_one_raw(directory: Path, label: str) -> Path:
    candidates = sorted(
        path
        for path in directory.glob("*raw_top25.tsv")
        if path.is_file() and "curated" not in path.name.lower()
    )
    if len(candidates) != 1:
        shown = ", ".join(str(path) for path in candidates) or "none"
        raise SystemExit(
            f"Could not identify exactly one {label} raw_top25 TSV in "
            f"{directory}. Found: {shown}"
        )
    return candidates[0]


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"Missing {label}: {path}")


def pick_column(fieldnames: list[str], choices: Iterable[str], label: str) -> str:
    for choice in choices:
        if choice in fieldnames:
            return choice
    raise ValueError(f"{label} not found; expected one of: {', '.join(choices)}")


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {path}")
        return list(reader.fieldnames), list(reader)


def read_emapper(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(errors="replace").splitlines()
    header_index = next(
        (index for index, line in enumerate(lines) if line.startswith("#query\t")),
        None,
    )
    if header_index is None:
        raise ValueError(f"Could not find #query header in eggNOG file: {path}")

    headers = lines[header_index][1:].split("\t")
    rows: list[dict[str, str]] = []
    for fields in csv.reader(lines[header_index + 1 :], delimiter="\t"):
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) < len(headers):
            fields += [""] * (len(headers) - len(fields))
        elif len(fields) > len(headers):
            fields = fields[: len(headers) - 1] + [
                " ".join(fields[len(headers) - 1 :])
            ]
        rows.append(dict(zip(headers, fields)))
    return rows


# ---------------------------------------------------------------------------
# 3. MEROPS parsing and evidence tiers
# ---------------------------------------------------------------------------

def parse_merops_title(title: str) -> tuple[str, str, str, str, str]:
    identifier_match = re.search(r"\[([^\]]+)\]", title)
    family_match = re.search(r"\]#([^#]+)#", title)
    organism_match = re.search(r"\(([^()]*)\)\s+\[[^\]]+\]#", title)

    identifier = identifier_match.group(1).strip() if identifier_match else "unparsed"
    family = family_match.group(1).strip() if family_match else "unparsed"
    organism = organism_match.group(1).strip() if organism_match else "unparsed"

    prefix = identifier[:1].upper()
    if prefix == "I" or family.upper().startswith("I"):
        category = "inhibitor_homology"
    elif prefix in PEPTIDASE_PREFIXES:
        category = "peptidase_homology"
    else:
        category = "other_or_unparsed"

    lower_title = title.lower()
    code_suffix = identifier.split(".", 1)[1].upper() if "." in identifier else ""

    explicitly_nonfunctional = (
        "non-peptidase homologue" in lower_title
        or "non-peptidase homolog" in lower_title
        or "non-inhibitory homologue" in lower_title
        or "non-inhibitory homolog" in lower_title
        or code_suffix.startswith("UN")
    )
    explicitly_unassigned = (
        "unassigned" in lower_title or code_suffix.startswith("UP")
    )

    if category == "inhibitor_homology":
        if explicitly_nonfunctional:
            status = "non_inhibitory_homologue_reference"
        elif explicitly_unassigned:
            status = "unassigned_inhibitor_family_reference"
        else:
            status = "named_inhibitor_reference"
    elif category == "peptidase_homology":
        if explicitly_nonfunctional:
            status = "non_peptidase_homologue_reference"
        elif explicitly_unassigned:
            status = "unassigned_peptidase_family_reference"
        else:
            status = "named_peptidase_reference"
    else:
        status = "unresolved_reference"

    return identifier, family, category, organism, status


def evidence_tier(reference_status: str) -> str:
    if reference_status.startswith("named_"):
        return "A_named_MEROPS_reference_match"
    if reference_status.startswith("unassigned_"):
        return "B_unassigned_family_homology"
    if reference_status.startswith("non_"):
        return "C_nonfunctional_reference_homology_only"
    return "D_unresolved_reference"


def tier_rank(reference_status: str) -> int:
    return {
        "named_inhibitor_reference": 0,
        "named_peptidase_reference": 0,
        "unassigned_inhibitor_family_reference": 1,
        "unassigned_peptidase_family_reference": 1,
        "non_inhibitory_homologue_reference": 2,
        "non_peptidase_homologue_reference": 2,
        "unresolved_reference": 3,
    }.get(reference_status, 4)


def alignment_sort_key(row: dict[str, object]) -> tuple[float, float, float, float]:
    return (
        float(row["evalue"]),
        -float(row["bitscore"]),
        -float(row["subject_span_coverage_percent"]),
        -float(row["percent_identity"]),
    )


def curated_sort_key(row: dict[str, object]) -> tuple[int, float, float, float, float]:
    return (tier_rank(str(row["merops_reference_status"])),) + alignment_sort_key(row)


def parse_raw_hits(path: Path) -> tuple[list[dict[str, object]], int, int]:
    rows: list[dict[str, object]] = []
    malformed_rows = 0
    invalid_coordinate_rows = 0

    with path.open(newline="") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            # The first 13 fields are fixed. A MEROPS title may contain a tab.
            if len(fields) > len(RAW_COLUMNS):
                fields = fields[: len(RAW_COLUMNS) - 1] + [
                    " ".join(fields[len(RAW_COLUMNS) - 1 :])
                ]
            if len(fields) != len(RAW_COLUMNS):
                malformed_rows += 1
                continue

            row: dict[str, object] = dict(zip(RAW_COLUMNS, fields))
            try:
                for column in (
                    "alignment_length",
                    "query_length",
                    "subject_length",
                    "query_start",
                    "query_end",
                    "subject_start",
                    "subject_end",
                ):
                    row[column] = int(str(row[column]))
                for column in (
                    "percent_identity",
                    "evalue",
                    "bitscore",
                    "diamond_query_coverage_percent",
                ):
                    row[column] = float(str(row[column]))
            except ValueError:
                malformed_rows += 1
                continue

            query_length = int(row["query_length"])
            subject_length = int(row["subject_length"])
            query_start = int(row["query_start"])
            query_end = int(row["query_end"])
            subject_start = int(row["subject_start"])
            subject_end = int(row["subject_end"])

            coordinates_valid = (
                query_length > 0
                and subject_length > 0
                and 1 <= query_start <= query_length
                and 1 <= query_end <= query_length
                and 1 <= subject_start <= subject_length
                and 1 <= subject_end <= subject_length
            )
            if not coordinates_valid:
                invalid_coordinate_rows += 1
                continue

            query_span = abs(query_end - query_start) + 1
            subject_span = abs(subject_end - subject_start) + 1
            row["query_coordinate_span"] = query_span
            row["subject_coordinate_span"] = subject_span
            row["query_span_coverage_percent"] = 100.0 * query_span / query_length
            row["subject_span_coverage_percent"] = (
                100.0 * subject_span / subject_length
            )
            # Retained only to document why the earlier value could exceed 100%.
            row["legacy_alignment_length_subject_coverage_percent"] = (
                100.0 * int(row["alignment_length"]) / subject_length
            )

            identifier, family, category, organism, status = parse_merops_title(
                str(row["merops_title"])
            )
            row["merops_identifier"] = identifier
            row["merops_family"] = family
            row["merops_category"] = category
            row["merops_reference_organism"] = organism
            row["merops_reference_status"] = status
            row["evidence_tier"] = evidence_tier(status)
            row["passes_corrected_filter"] = (
                float(row["evalue"]) <= MAX_EVALUE
                and float(row["percent_identity"]) >= MIN_IDENTITY
                and float(row["subject_span_coverage_percent"])
                >= MIN_SUBJECT_SPAN_COVERAGE
            )
            row["passes_legacy_filter_for_comparison"] = (
                float(row["evalue"]) <= MAX_EVALUE
                and float(row["percent_identity"]) >= MIN_IDENTITY
                and float(row["legacy_alignment_length_subject_coverage_percent"])
                >= MIN_SUBJECT_SPAN_COVERAGE
            )
            rows.append(row)

    return rows, malformed_rows, invalid_coordinate_rows


def collapse_gene_category(
    filtered_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in filtered_rows:
        grouped[(str(row["gene_id"]), str(row["merops_category"]))].append(row)

    collapsed: list[dict[str, object]] = []
    for (_, _), hits in grouped.items():
        best_overall = min(hits, key=alignment_sort_key)
        representative = min(hits, key=curated_sort_key)
        output = dict(representative)
        output["best_overall_merops_reference_id"] = best_overall[
            "merops_reference_id"
        ]
        output["best_overall_merops_identifier"] = best_overall[
            "merops_identifier"
        ]
        output["best_overall_merops_family"] = best_overall["merops_family"]
        output["best_overall_reference_status"] = best_overall[
            "merops_reference_status"
        ]
        output["best_overall_evalue"] = best_overall["evalue"]
        output["best_overall_bitscore"] = best_overall["bitscore"]
        output["passing_hits_total"] = len(hits)
        output["passing_named_reference_hits"] = sum(
            str(hit["merops_reference_status"]).startswith("named_")
            for hit in hits
        )
        output["passing_unassigned_reference_hits"] = sum(
            str(hit["merops_reference_status"]).startswith("unassigned_")
            for hit in hits
        )
        output["passing_nonfunctional_reference_hits"] = sum(
            str(hit["merops_reference_status"]).startswith("non_")
            for hit in hits
        )
        output["curated_AB_candidate"] = tier_rank(
            str(representative["merops_reference_status"])
        ) <= 1
        collapsed.append(output)

    return sorted(
        collapsed,
        key=lambda row: (
            str(row["merops_category"]),
            tier_rank(str(row["merops_reference_status"])),
            alignment_sort_key(row),
            str(row["gene_id"]),
        ),
    )


RAW_DERIVED_COLUMNS = [
    "query_coordinate_span",
    "subject_coordinate_span",
    "query_span_coverage_percent",
    "subject_span_coverage_percent",
    "legacy_alignment_length_subject_coverage_percent",
    "merops_identifier",
    "merops_family",
    "merops_category",
    "merops_reference_organism",
    "merops_reference_status",
    "evidence_tier",
    "passes_corrected_filter",
    "passes_legacy_filter_for_comparison",
]

COLLAPSE_COLUMNS = [
    "best_overall_merops_reference_id",
    "best_overall_merops_identifier",
    "best_overall_merops_family",
    "best_overall_reference_status",
    "best_overall_evalue",
    "best_overall_bitscore",
    "passing_hits_total",
    "passing_named_reference_hits",
    "passing_unassigned_reference_hits",
    "passing_nonfunctional_reference_hits",
    "curated_AB_candidate",
]


# ---------------------------------------------------------------------------
# 4. Metadata readers and membership comparison
# ---------------------------------------------------------------------------

def read_old_ids(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    fields, rows = read_tsv(path)
    id_column = pick_column(fields, ["gene_id", "Geneid", "GeneID"], "Old gene ID")
    return {row[id_column] for row in rows if row[id_column]}


def membership_rows(
    old_ids: set[str],
    broad_rows: list[dict[str, object]],
    curated_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    broad_by_id = {str(row["gene_id"]): row for row in broad_rows}
    curated_ids = {str(row["gene_id"]) for row in curated_rows}
    all_ids = sorted(old_ids | set(broad_by_id))
    output: list[dict[str, object]] = []

    for gene_id in all_ids:
        old_member = gene_id in old_ids
        broad_member = gene_id in broad_by_id
        curated_member = gene_id in curated_ids
        row = broad_by_id.get(gene_id, {})

        if old_member and not broad_member:
            status = "removed_by_corrected_coverage_filter"
        elif old_member and broad_member and not curated_member:
            status = "retained_homology_only_nonfunctional_or_unresolved_reference"
        elif old_member and curated_member:
            status = "retained_curated_AB_candidate"
        elif not old_member and curated_member:
            status = "new_curated_AB_candidate_after_reselection"
        elif not old_member and broad_member:
            status = "new_broad_homology_after_reselection"
        else:
            status = "unclassified"

        output.append(
            {
                "gene_id": gene_id,
                "old_candidate_member": old_member,
                "v2_broad_homology_member": broad_member,
                "v2_curated_AB_candidate": curated_member,
                "membership_status": status,
                "v2_merops_identifier": row.get("merops_identifier", ""),
                "v2_merops_family": row.get("merops_family", ""),
                "v2_reference_status": row.get("merops_reference_status", ""),
                "v2_evidence_tier": row.get("evidence_tier", ""),
                "v2_subject_span_coverage_percent": row.get(
                    "subject_span_coverage_percent", ""
                ),
            }
        )
    return output


def attach_insect_metadata(
    collapsed: list[dict[str, object]],
    integrated_path: Path,
    errors: list[str],
) -> tuple[list[dict[str, object]], int]:
    fields, rows = read_tsv(integrated_path)
    id_column = pick_column(fields, ["gene_id", "Geneid", "GeneID"], "Insect gene ID")
    status_column = pick_column(
        fields,
        ["DESeq2_status", "deseq2_status", "DE_status", "DESeq2_primary_status"],
        "Insect DESeq2 status",
    )
    class_column = pick_column(
        fields,
        ["analysis_class", "Analysis_class", "signalp_tmhmm_class"],
        "Insect analysis class",
    )
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        gene_id = row[id_column]
        if gene_id in index:
            errors.append(f"Duplicate insect integrated gene ID: {gene_id}")
        index[gene_id] = row

    optional = {
        "signalp_prediction": ["SignalP_prediction", "signalp_prediction", "SignalP"],
        "deeptmhmm_type": ["deeptmhmm_type", "DeepTMHMM", "deeptmhmm_prediction"],
        "predicted_TMRs": ["predicted_TMRs", "TMRs", "predicted_tmrs"],
        "baseMean": ["baseMean", "base_mean"],
        "log2FoldChange": ["log2FoldChange", "log2FC"],
        "padj": ["padj", "adjusted_pvalue"],
    }

    output: list[dict[str, object]] = []
    for candidate in collapsed:
        gene_id = str(candidate["gene_id"])
        metadata = index.get(gene_id)
        if metadata is None:
            errors.append(f"Insect candidate absent from integrated table: {gene_id}")
            continue
        row = dict(candidate)
        row["DESeq2_status"] = metadata[status_column]
        row["analysis_class"] = metadata[class_column]
        for output_name, choices in optional.items():
            source = next((choice for choice in choices if choice in fields), None)
            row[output_name] = metadata.get(source, "") if source else ""
        row["classical_secreted_core"] = (
            row["analysis_class"] == "classical_secreted_core"
        )
        output.append(row)

    return output, len(index)


def attach_plant_metadata(
    collapsed: list[dict[str, object]],
    selected_path: Path,
    annotation_path: Path,
    errors: list[str],
) -> tuple[list[dict[str, object]], int, int]:
    fields, rows = read_tsv(selected_path)
    id_column = pick_column(fields, ["gene_id", "Geneid", "GeneID"], "Plant gene ID")
    selected: dict[str, dict[str, str]] = {}
    for row in rows:
        gene_id = row[id_column]
        if gene_id in selected:
            errors.append(f"Duplicate selected plant gene ID: {gene_id}")
        selected[gene_id] = row

    annotation_index: dict[str, dict[str, str]] = {}
    for annotation in read_emapper(annotation_path):
        query = annotation.get("query", annotation.get("#query", "")).split()[0]
        gene_id = re.sub(r"_i[0-9]+\.p[0-9]+$", "", query)
        if gene_id:
            if gene_id in annotation_index:
                errors.append(f"Multiple eggNOG annotations map to plant gene: {gene_id}")
            annotation_index[gene_id] = annotation

    output: list[dict[str, object]] = []
    for candidate in collapsed:
        gene_id = str(candidate["gene_id"])
        metadata = selected.get(gene_id)
        if metadata is None:
            errors.append(f"Plant candidate absent from selected set: {gene_id}")
            continue
        annotation = annotation_index.get(gene_id, {})
        row = dict(candidate)
        row["DESeq2_status"] = metadata.get(
            "DESeq2_status", metadata.get("deseq2_status", "")
        )
        row["analysis_class"] = metadata.get("analysis_class", "")
        for column in ("baseMean", "log2FoldChange", "padj"):
            row[column] = metadata.get(column, "")
        row["seed_ortholog"] = annotation.get("seed_ortholog", "")
        row["seed_taxonomy_code"] = (
            str(row["seed_ortholog"]).split(".", 1)[0]
            if row["seed_ortholog"] not in {"", "-"}
            else "no_seed_ortholog"
        )
        for column in (
            "Description",
            "Preferred_name",
            "PFAMs",
            "CAZy",
            "GOs",
            "KEGG_ko",
            "eggNOG_OGs",
        ):
            row[column] = annotation.get(column, "")
        output.append(row)

    return output, len(selected), len(annotation_index)


# ---------------------------------------------------------------------------
# 5. Run both reparses
# ---------------------------------------------------------------------------

RESULTS.mkdir(parents=True, exist_ok=True)

insect_raw = find_one_raw(INSECT_RAW_DIR, "insect")
plant_raw = find_one_raw(PLANT_RAW_DIR, "plant")
for required_path, label in (
    (INSECT_INTEGRATED, "insect integrated table"),
    (PLANT_SELECTED, "plant selected-secretome table"),
    (PLANT_EGGNOG, "plant representative eggNOG annotation"),
):
    require_file(required_path, label)

errors: list[str] = []
warnings: list[str] = []

insect_raw_rows, insect_malformed, insect_invalid_coordinates = parse_raw_hits(
    insect_raw
)
plant_raw_rows, plant_malformed, plant_invalid_coordinates = parse_raw_hits(
    plant_raw
)

insect_filtered = [
    row for row in insect_raw_rows if bool(row["passes_corrected_filter"])
]
plant_filtered = [
    row for row in plant_raw_rows if bool(row["passes_corrected_filter"])
]
insect_legacy_filtered = [
    row
    for row in insect_raw_rows
    if bool(row["passes_legacy_filter_for_comparison"])
]
plant_legacy_filtered = [
    row
    for row in plant_raw_rows
    if bool(row["passes_legacy_filter_for_comparison"])
]

insect_collapsed = collapse_gene_category(insect_filtered)
plant_collapsed = collapse_gene_category(plant_filtered)

insect_all, insect_integrated_count = attach_insect_metadata(
    insect_collapsed, INSECT_INTEGRATED, errors
)
plant_all, plant_selected_count, annotation_count = attach_plant_metadata(
    plant_collapsed, PLANT_SELECTED, PLANT_EGGNOG, errors
)

insect_broad = [
    row
    for row in insect_all
    if row["merops_category"] == "inhibitor_homology"
    and bool(row["classical_secreted_core"])
]
insect_curated_ab = [
    row for row in insect_broad if bool(row["curated_AB_candidate"])
]

plant_broad = [
    row for row in plant_all if row["merops_category"] == "peptidase_homology"
]
plant_curated_ab = [
    row for row in plant_broad if bool(row["curated_AB_candidate"])
]

old_insect_ids = read_old_ids(OLD_INSECT_CANDIDATES)
old_plant_ids = read_old_ids(OLD_PLANT_CANDIDATES)
if not old_insect_ids:
    warnings.append("Original insect candidate table was absent; membership comparison is incomplete")
if not old_plant_ids:
    warnings.append("Original plant candidate table was absent; membership comparison is incomplete")

insect_changes = membership_rows(old_insect_ids, insect_broad, insect_curated_ab)
plant_changes = membership_rows(old_plant_ids, plant_broad, plant_curated_ab)


# ---------------------------------------------------------------------------
# 6. Structural audit
# ---------------------------------------------------------------------------

def unique_gene_count(rows: list[dict[str, object]]) -> int:
    return len({str(row["gene_id"]) for row in rows})


insect_raw_hit_genes = unique_gene_count(insect_raw_rows)
plant_raw_hit_genes = unique_gene_count(plant_raw_rows)

if len(insect_raw_rows) != EXPECTED_INSECT_RAW_ROWS:
    errors.append(
        f"Insect raw rows={len(insect_raw_rows)}; expected {EXPECTED_INSECT_RAW_ROWS}"
    )
if insect_raw_hit_genes != EXPECTED_INSECT_RAW_HIT_GENES:
    errors.append(
        f"Insect raw-hit genes={insect_raw_hit_genes}; "
        f"expected {EXPECTED_INSECT_RAW_HIT_GENES}"
    )
if insect_integrated_count != EXPECTED_INSECT_INTEGRATED_ROWS:
    errors.append(
        f"Insect integrated genes={insect_integrated_count}; "
        f"expected {EXPECTED_INSECT_INTEGRATED_ROWS}"
    )
if len(plant_raw_rows) != EXPECTED_PLANT_RAW_ROWS:
    errors.append(
        f"Plant raw rows={len(plant_raw_rows)}; expected {EXPECTED_PLANT_RAW_ROWS}"
    )
if plant_raw_hit_genes != EXPECTED_PLANT_RAW_HIT_GENES:
    errors.append(
        f"Plant raw-hit genes={plant_raw_hit_genes}; "
        f"expected {EXPECTED_PLANT_RAW_HIT_GENES}"
    )
if plant_selected_count != EXPECTED_PLANT_SELECTED_GENES:
    errors.append(
        f"Plant selected genes={plant_selected_count}; "
        f"expected {EXPECTED_PLANT_SELECTED_GENES}"
    )
if insect_malformed or plant_malformed:
    errors.append(
        f"Malformed rows: insect={insect_malformed}, plant={plant_malformed}"
    )
if insect_invalid_coordinates or plant_invalid_coordinates:
    errors.append(
        "Invalid-coordinate rows: "
        f"insect={insect_invalid_coordinates}, plant={plant_invalid_coordinates}"
    )
if any(float(row["subject_span_coverage_percent"]) > 100.0 for row in insect_raw_rows + plant_raw_rows):
    errors.append("Corrected subject-span coverage exceeds 100%")
if unique_gene_count(insect_broad) != len(insect_broad):
    errors.append("Insect broad secreted inhibitor set contains duplicate genes")
if unique_gene_count(plant_broad) != len(plant_broad):
    errors.append("Plant broad peptidase set contains duplicate genes")
if not {str(row["gene_id"]) for row in insect_curated_ab}.issubset(
    {str(row["gene_id"]) for row in insect_broad}
):
    errors.append("Insect curated A+B set is not a subset of its broad set")
if not {str(row["gene_id"]) for row in plant_curated_ab}.issubset(
    {str(row["gene_id"]) for row in plant_broad}
):
    errors.append("Plant curated A+B set is not a subset of its broad set")
if any(
    str(row.get("DESeq2_status", "")) not in VALID_DESEQ2_STATUSES
    for row in insect_broad + plant_broad
):
    errors.append("Unexpected or missing DESeq2 status in a broad candidate set")


# ---------------------------------------------------------------------------
# 7. Write new v2 outputs only
# ---------------------------------------------------------------------------

raw_output_columns = RAW_COLUMNS + RAW_DERIVED_COLUMNS
insect_candidate_columns = raw_output_columns + COLLAPSE_COLUMNS + [
    "DESeq2_status",
    "analysis_class",
    "classical_secreted_core",
    "signalp_prediction",
    "deeptmhmm_type",
    "predicted_TMRs",
    "baseMean",
    "log2FoldChange",
    "padj",
]
plant_candidate_columns = raw_output_columns + COLLAPSE_COLUMNS + [
    "DESeq2_status",
    "analysis_class",
    "baseMean",
    "log2FoldChange",
    "padj",
    "seed_taxonomy_code",
    "seed_ortholog",
    "Description",
    "Preferred_name",
    "PFAMs",
    "CAZy",
    "GOs",
    "KEGG_ko",
    "eggNOG_OGs",
]
change_columns = [
    "gene_id",
    "old_candidate_member",
    "v2_broad_homology_member",
    "v2_curated_AB_candidate",
    "membership_status",
    "v2_merops_identifier",
    "v2_merops_family",
    "v2_reference_status",
    "v2_evidence_tier",
    "v2_subject_span_coverage_percent",
]

output_paths = {
    "insect_filtered": RESULTS / "insect_MEROPS_conservative_hits_curated_v2.tsv",
    "insect_all": RESULTS / "insect_MEROPS_gene_category_curated_v2.tsv",
    "insect_broad": RESULTS / "insect_MEROPS_classical_secreted_inhibitor_homologues_curated_v2.tsv",
    "insect_curated": RESULTS / "insect_MEROPS_classical_secreted_inhibitor_candidates_curated_v2.tsv",
    "insect_changes": RESULTS / "insect_MEROPS_candidate_membership_changes_v2.tsv",
    "plant_filtered": RESULTS / "plant_GPup_secretome_MEROPS_conservative_hits_curated_v2.tsv",
    "plant_all": RESULTS / "plant_GPup_secretome_MEROPS_gene_category_curated_v2.tsv",
    "plant_broad": RESULTS / "plant_GPup_secretome_MEROPS_peptidase_homologues_curated_v2.tsv",
    "plant_curated": RESULTS / "plant_GPup_secretome_MEROPS_peptidase_candidates_curated_v2.tsv",
    "plant_changes": RESULTS / "plant_GPup_MEROPS_candidate_membership_changes_v2.tsv",
    "summary": RESULTS / "MEROPS_reparse_curation_v2_summary.tsv",
    "parameters": RESULTS / "MEROPS_reparse_curation_v2_parameters.tsv",
    "audit": RESULTS / "MEROPS_reparse_curation_v2_audit.txt",
}

write_tsv(output_paths["insect_filtered"], insect_filtered, raw_output_columns)
write_tsv(output_paths["insect_all"], insect_all, insect_candidate_columns)
write_tsv(output_paths["insect_broad"], insect_broad, insect_candidate_columns)
write_tsv(output_paths["insect_curated"], insect_curated_ab, insect_candidate_columns)
write_tsv(output_paths["insect_changes"], insect_changes, change_columns)
write_tsv(output_paths["plant_filtered"], plant_filtered, raw_output_columns)
write_tsv(output_paths["plant_all"], plant_all, plant_candidate_columns)
write_tsv(output_paths["plant_broad"], plant_broad, plant_candidate_columns)
write_tsv(output_paths["plant_curated"], plant_curated_ab, plant_candidate_columns)
write_tsv(output_paths["plant_changes"], plant_changes, change_columns)

insect_tiers = Counter(str(row["evidence_tier"]) for row in insect_broad)
plant_tiers = Counter(str(row["evidence_tier"]) for row in plant_broad)
insect_statuses = Counter(str(row["DESeq2_status"]) for row in insect_curated_ab)

summary_rows: list[dict[str, object]] = [
    {"section": "insect", "metric": "raw_rows", "count": len(insect_raw_rows)},
    {"section": "insect", "metric": "legacy_filter_passing_hits", "count": len(insect_legacy_filtered)},
    {"section": "insect", "metric": "corrected_filter_passing_hits", "count": len(insect_filtered)},
    {"section": "insect", "metric": "legacy_passing_hits_removed_by_corrected_coverage", "count": sum(bool(row["passes_legacy_filter_for_comparison"]) and not bool(row["passes_corrected_filter"]) for row in insect_raw_rows)},
    {"section": "insect", "metric": "raw_hits_with_legacy_coverage_over_100", "count": sum(float(row["legacy_alignment_length_subject_coverage_percent"]) > 100.0 for row in insect_raw_rows)},
    {"section": "insect", "metric": "old_secreted_inhibitor_candidates", "count": len(old_insect_ids)},
    {"section": "insect", "metric": "v2_broad_secreted_inhibitor_homologues", "count": len(insect_broad)},
    {"section": "insect", "metric": "v2_curated_AB_secreted_inhibitor_candidates", "count": len(insect_curated_ab)},
    {"section": "plant", "metric": "raw_rows", "count": len(plant_raw_rows)},
    {"section": "plant", "metric": "legacy_filter_passing_hits", "count": len(plant_legacy_filtered)},
    {"section": "plant", "metric": "corrected_filter_passing_hits", "count": len(plant_filtered)},
    {"section": "plant", "metric": "legacy_passing_hits_removed_by_corrected_coverage", "count": sum(bool(row["passes_legacy_filter_for_comparison"]) and not bool(row["passes_corrected_filter"]) for row in plant_raw_rows)},
    {"section": "plant", "metric": "raw_hits_with_legacy_coverage_over_100", "count": sum(float(row["legacy_alignment_length_subject_coverage_percent"]) > 100.0 for row in plant_raw_rows)},
    {"section": "plant", "metric": "old_peptidase_candidates", "count": len(old_plant_ids)},
    {"section": "plant", "metric": "v2_broad_peptidase_homologues", "count": len(plant_broad)},
    {"section": "plant", "metric": "v2_curated_AB_peptidase_candidates", "count": len(plant_curated_ab)},
]
for tier in sorted(set(insect_tiers) | set(plant_tiers)):
    summary_rows.append(
        {"section": "insect", "metric": f"broad_{tier}", "count": insect_tiers[tier]}
    )
    summary_rows.append(
        {"section": "plant", "metric": f"broad_{tier}", "count": plant_tiers[tier]}
    )
for status in sorted(insect_statuses):
    summary_rows.append(
        {"section": "insect", "metric": f"curated_AB_DESeq2_{status}", "count": insect_statuses[status]}
    )
write_tsv(output_paths["summary"], summary_rows, ["section", "metric", "count"])

parameter_rows = [
    {"parameter": "scope", "value": "existing DIAMOND hits only; DIAMOND not rerun"},
    {"parameter": "maximum_evalue", "value": str(MAX_EVALUE)},
    {"parameter": "minimum_percent_identity", "value": str(MIN_IDENTITY)},
    {"parameter": "minimum_subject_span_coverage_percent", "value": str(MIN_SUBJECT_SPAN_COVERAGE)},
    {"parameter": "coverage_definition", "value": "100 * (abs(subject_end-subject_start)+1) / subject_length"},
    {"parameter": "tier_A", "value": "best passing named MEROPS reference match; predicted function, not experimentally validated query activity"},
    {"parameter": "tier_B", "value": "best passing unassigned family reference match; family homology only"},
    {"parameter": "tier_C", "value": "only explicitly non-peptidase/non-inhibitory reference support"},
    {"parameter": "curated_AB_definition", "value": "tier A or B; tier C/D retained only in broad homology tables"},
    {"parameter": "representative_hit_selection", "value": "best evidence tier first, then evalue, bitscore, corrected subject-span coverage, percent identity; best overall score fields retained separately"},
    {"parameter": "insect_raw_sha256", "value": sha256(insect_raw)},
    {"parameter": "plant_raw_sha256", "value": sha256(plant_raw)},
    {"parameter": "insect_integrated_sha256", "value": sha256(INSECT_INTEGRATED)},
    {"parameter": "plant_selected_sha256", "value": sha256(PLANT_SELECTED)},
    {"parameter": "plant_eggnog_sha256", "value": sha256(PLANT_EGGNOG)},
]
write_tsv(output_paths["parameters"], parameter_rows, ["parameter", "value"])

report = [
    "===== MEROPS REPARSE AND CURATION V2 =====",
    f"Insect raw file: {insect_raw}",
    f"Plant raw file: {plant_raw}",
    "",
    f"Insect raw rows: {len(insect_raw_rows)}",
    f"Insect legacy-filter passing hits: {len(insect_legacy_filtered)}",
    f"Insect corrected-filter passing hits: {len(insect_filtered)}",
    f"Insect legacy-passing hits removed by corrected coverage: {sum(bool(row['passes_legacy_filter_for_comparison']) and not bool(row['passes_corrected_filter']) for row in insect_raw_rows)}",
    f"Insect raw hits with legacy coverage over 100%: {sum(float(row['legacy_alignment_length_subject_coverage_percent']) > 100.0 for row in insect_raw_rows)}",
    f"Insect old secreted inhibitor candidates: {len(old_insect_ids)}",
    f"Insect v2 broad secreted inhibitor homologues: {len(insect_broad)}",
    f"Insect v2 curated A+B inhibitor candidates: {len(insect_curated_ab)}",
    "Insect broad evidence tiers:",
]
report.extend(f"  {tier}: {count}" for tier, count in sorted(insect_tiers.items()))
report.extend(
    [
        "Insect curated A+B DESeq2 statuses:",
        *(
            f"  {status}: {count}"
            for status, count in sorted(insect_statuses.items())
        ),
        "",
        f"Plant raw rows: {len(plant_raw_rows)}",
        f"Plant legacy-filter passing hits: {len(plant_legacy_filtered)}",
        f"Plant corrected-filter passing hits: {len(plant_filtered)}",
        f"Plant legacy-passing hits removed by corrected coverage: {sum(bool(row['passes_legacy_filter_for_comparison']) and not bool(row['passes_corrected_filter']) for row in plant_raw_rows)}",
        f"Plant raw hits with legacy coverage over 100%: {sum(float(row['legacy_alignment_length_subject_coverage_percent']) > 100.0 for row in plant_raw_rows)}",
        f"Plant old peptidase candidates: {len(old_plant_ids)}",
        f"Plant v2 broad peptidase homologues: {len(plant_broad)}",
        f"Plant v2 curated A+B peptidase candidates: {len(plant_curated_ab)}",
        "Plant broad evidence tiers:",
    ]
)
report.extend(f"  {tier}: {count}" for tier, count in sorted(plant_tiers.items()))
report.extend(
    [
        "",
        f"Representative plant eggNOG annotations indexed: {annotation_count}",
        f"Warnings: {len(warnings)}",
        *(f"WARNING: {warning}" for warning in warnings),
        f"Errors detected: {len(errors)}",
        *(f"ERROR: {error}" for error in errors),
        "MEROPS REPARSE AND CURATION V2: " + ("PASS" if not errors else "CHECK REQUIRED"),
    ]
)
output_paths["audit"].write_text("\n".join(report) + "\n")

print("\n".join(report))
print(f"\nSummary: {output_paths['summary']}")
print(f"Insect membership changes: {output_paths['insect_changes']}")
print(f"Plant membership changes: {output_paths['plant_changes']}")
print(f"Audit: {output_paths['audit']}")

if errors:
    raise SystemExit(1)
