#!/usr/bin/env python3
"""Reconstruct and audit the final insect SignalP 6/DeepTMHMM classification.

The script reads the frozen 13,455-row objective merge table, applies the
coordinate-based classification used in the analysis, appends the three
manually reviewed SignalP/TMR cases, and writes the final classification table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import Counter
from pathlib import Path


EXPECTED_RECORDS = 13_455
EXPECTED_INPUT_SHA256 = (
    "4ef04acf43fa4152c68e51a39b261cbfcec1f879ee75fd88c8b4850204a6daf2"
)
EXPECTED_OUTPUT_SHA256 = (
    "440116b8580abfc8a19f33536ed9b59a97323cd83ebc33f245321ab618de883d"
)
EXPECTED_CLASSES = {
    "classical_secreted_core": 2_039,
    "membrane_without_signalp": 1_971,
    "no_signalp_no_tmr": 9_037,
    "signalp_tm_discordant_review": 3,
    "signalp_with_downstream_tmr": 405,
}
EXPECTED_INPUT_COLUMNS = [
    "gene_id",
    "length_aa",
    "signalp_prediction",
    "signalp_other_probability",
    "signalp_sp_probability",
    "signalp_cs_start",
    "signalp_cs_end",
    "signalp_cs_probability",
    "deeptmhmm_type",
    "deeptmhmm_predicted_tmrs",
    "deeptmhmm_signal_regions",
    "deeptmhmm_tmr_regions",
    "combined_group",
]
OUTPUT_COLUMNS = EXPECTED_INPUT_COLUMNS + [
    "downstream_tmr_count",
    "analysis_class",
    "primary_classical_secreted",
    "manual_review",
    "review_interpretation",
    "review_note",
]

MANUAL_REVIEWS = {
    "egapxtmp_007701": (
        "likely_membrane_anchor",
        "ABHD2-like alpha/beta hydrolase; N-terminal membrane anchor supported",
    ),
    "egapxtmp_008737": (
        "likely_secreted",
        "Nudel-like serine protease with Trypsin and LDL-receptor domains",
    ),
    "egapxtmp_014139": (
        "likely_secreted",
        "DUF229 protein with extracellular-region and extracellular-space GO annotations",
    ),
}

SEGMENT_PATTERN = re.compile(r"([MB]):(\d+)-(\d+)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    default_base = Path.home() / "projects/qe_genome_guided/11_signalp_tmhmm"
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct the final 13,455-protein insect SignalP 6/DeepTMHMM "
            "classification from the frozen objective merge table."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=default_base / "results/signalp6_deeptmhmm_genomewide_13455.tsv",
        help="Objective SignalP/DeepTMHMM merge table.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            default_base
            / "results/signalp6_deeptmhmm_final_classification_13455.tsv"
        ),
        help="Final classification table to create or verify.",
    )
    return parser.parse_args()


def parse_segments(value: str, gene_id: str, errors: list[str]) -> list[tuple[str, int, int]]:
    if not value:
        return []
    tokens = value.split(";")
    segments: list[tuple[str, int, int]] = []
    for token in tokens:
        match = SEGMENT_PATTERN.fullmatch(token)
        if not match:
            errors.append(f"{gene_id}: malformed TMR segment: {token}")
            continue
        symbol, start_text, end_text = match.groups()
        start, end = int(start_text), int(end_text)
        if start < 1 or end < start:
            errors.append(f"{gene_id}: invalid TMR coordinates: {token}")
        segments.append((symbol, start, end))
    return segments


def classify_row(row: dict[str, str], errors: list[str]) -> dict[str, str]:
    gene_id = row["gene_id"]
    signalp = row["signalp_prediction"]

    try:
        tmr_count = int(row["deeptmhmm_predicted_tmrs"])
    except ValueError:
        errors.append(f"{gene_id}: invalid DeepTMHMM TMR count")
        tmr_count = -1

    segments = parse_segments(row["deeptmhmm_tmr_regions"], gene_id, errors)
    if len(segments) != tmr_count:
        errors.append(
            f"{gene_id}: reported TMR count {tmr_count} differs from "
            f"the {len(segments)} parsed segments"
        )

    if signalp == "SP":
        try:
            cleavage_after = int(row["signalp_cs_start"])
            cleavage_right = int(row["signalp_cs_end"])
        except ValueError:
            errors.append(f"{gene_id}: SP prediction lacks valid cleavage coordinates")
            cleavage_after = cleavage_right = 0
        if cleavage_right != cleavage_after + 1:
            errors.append(f"{gene_id}: nonadjacent SignalP cleavage coordinates")
        # A downstream TMR is one whose first residue occurs after the predicted
        # cleavage boundary. The three N-terminal overlap cases therefore have
        # a downstream count of zero and enter manual review.
        downstream_tmr_count = sum(start > cleavage_after for _, start, _ in segments)
    elif signalp == "OTHER":
        downstream_tmr_count = 0
    else:
        errors.append(f"{gene_id}: unexpected SignalP prediction: {signalp}")
        downstream_tmr_count = 0

    expected_group = (
        f"SignalP_{signalp}__DeepTMHMM_"
        f"{'TMR_present' if tmr_count > 0 else 'no_TMR'}"
    )
    if row["combined_group"] != expected_group:
        errors.append(
            f"{gene_id}: combined group {row['combined_group']} differs from "
            f"{expected_group}"
        )

    if signalp == "SP" and tmr_count == 0:
        analysis_class = "classical_secreted_core"
    elif signalp == "SP" and downstream_tmr_count > 0:
        analysis_class = "signalp_with_downstream_tmr"
    elif signalp == "SP" and tmr_count > 0:
        analysis_class = "signalp_tm_discordant_review"
    elif signalp == "OTHER" and tmr_count > 0:
        analysis_class = "membrane_without_signalp"
    elif signalp == "OTHER" and tmr_count == 0:
        analysis_class = "no_signalp_no_tmr"
    else:
        analysis_class = "unclassified"
        errors.append(f"{gene_id}: could not assign an analysis class")

    manual_review = analysis_class == "signalp_tm_discordant_review"
    if manual_review:
        if gene_id not in MANUAL_REVIEWS:
            errors.append(f"{gene_id}: discordant case lacks a manual interpretation")
            interpretation, note = "", ""
        else:
            interpretation, note = MANUAL_REVIEWS[gene_id]
    else:
        interpretation, note = "", ""

    result = dict(row)
    result.update(
        {
            "downstream_tmr_count": str(downstream_tmr_count),
            "analysis_class": analysis_class,
            "primary_classical_secreted": (
                "yes" if analysis_class == "classical_secreted_core" else "no"
            ),
            "manual_review": "yes" if manual_review else "no",
            "review_interpretation": interpretation,
            "review_note": note,
        }
    )
    return result


def main() -> int:
    args = parse_args()
    input_path = args.input.expanduser().resolve()
    output_path = args.output.expanduser().resolve()
    errors: list[str] = []

    if not input_path.is_file():
        print(f"ERROR: Missing objective merge table: {input_path}")
        return 1
    observed_input_hash = sha256(input_path)
    if observed_input_hash != EXPECTED_INPUT_SHA256:
        errors.append(
            "Objective merge-table SHA-256 differs from the frozen analysis input"
        )

    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_columns = reader.fieldnames or []
        rows = list(reader)

    if input_columns != EXPECTED_INPUT_COLUMNS:
        errors.append(
            "Unexpected input columns or column order: " + ", ".join(input_columns)
        )
    gene_ids = [row.get("gene_id", "") for row in rows]
    if len(rows) != EXPECTED_RECORDS:
        errors.append(f"Input rows={len(rows)}; expected {EXPECTED_RECORDS}")
    if len(set(gene_ids)) != EXPECTED_RECORDS:
        errors.append(
            f"Unique input IDs={len(set(gene_ids))}; expected {EXPECTED_RECORDS}"
        )

    output_rows = [classify_row(row, errors) for row in rows]
    class_counts = Counter(row["analysis_class"] for row in output_rows)
    if dict(class_counts) != EXPECTED_CLASSES:
        errors.append(
            f"Analysis-class distribution differs: {dict(class_counts)}"
        )

    review_ids = {
        row["gene_id"] for row in output_rows if row["manual_review"] == "yes"
    }
    if review_ids != set(MANUAL_REVIEWS):
        errors.append(
            "Manual-review ID set differs: " + ", ".join(sorted(review_ids))
        )
    primary_count = sum(
        row["primary_classical_secreted"] == "yes" for row in output_rows
    )
    if primary_count != EXPECTED_CLASSES["classical_secreted_core"]:
        errors.append(f"Primary secretome={primary_count}; expected 2039")

    print("===== INSECT SIGNALP-DEEPTMHMM FINAL CLASSIFICATION =====")
    print(f"Input rows: {len(rows)}")
    print(f"Unique IDs: {len(set(gene_ids))}")
    print("\nFinal analysis-class distribution:")
    for analysis_class in sorted(EXPECTED_CLASSES):
        count = class_counts.get(analysis_class, 0)
        print(
            f"{analysis_class}: {count} "
            f"({100 * count / EXPECTED_RECORDS:.2f}%)"
        )
    print(f"\nPrimary high-confidence classical insect secretome: {primary_count}")
    print(f"Manual SP/TM review tier: {len(review_ids)}")
    print("\nManual-review proteins:")
    for gene_id in sorted(MANUAL_REVIEWS):
        interpretation, note = MANUAL_REVIEWS[gene_id]
        print(f"{gene_id}\t{interpretation}\t{note}")

    if errors:
        print(f"\nErrors detected: {len(errors)}")
        for error in errors[:100]:
            print("ERROR:", error)
        print("INSECT SIGNALP-DEEPTMHMM FINAL CLASSIFICATION: CHECK REQUIRED")
        return 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".reconstructed.tmp")
    with temporary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    reconstructed_hash = sha256(temporary_path)
    if reconstructed_hash != EXPECTED_OUTPUT_SHA256:
        print("\nErrors detected: 1")
        print(
            "ERROR: Reconstructed output SHA-256 differs from the frozen final table"
        )
        print(f"Reconstructed temporary output: {temporary_path}")
        print(f"Observed SHA-256: {reconstructed_hash}")
        print(f"Expected SHA-256: {EXPECTED_OUTPUT_SHA256}")
        print("INSECT SIGNALP-DEEPTMHMM FINAL CLASSIFICATION: CHECK REQUIRED")
        return 1

    temporary_path.replace(output_path)
    print("\nErrors detected: 0")
    print(f"Output: {output_path}")
    print(f"Output rows: {len(output_rows)}")
    print(f"Input SHA-256: {observed_input_hash}")
    print(f"Output SHA-256: {sha256(output_path)}")
    print("INSECT SIGNALP-DEEPTMHMM FINAL CLASSIFICATION: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
