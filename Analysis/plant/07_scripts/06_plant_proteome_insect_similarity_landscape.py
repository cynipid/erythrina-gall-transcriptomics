#!/usr/bin/env python3
"""Map the plant representative proteome to the Q. erythrinae proteome.

This is a provenance quality-control analysis. It creates review labels and a
proposed conservative reanalysis mask; it does not modify a count matrix,
rerun DESeq2, or change any frozen discovery ranking.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


HOME = Path.home()
ROOT = HOME / "projects/qe_cross_species_discovery"
LINKS = ROOT / "01_input_links"

INPUTS = {
    "plant_integration": LINKS / "plant_signalp_tmhmm_deseq2_integrated.tsv",
    "plant_protein_map": LINKS / "plant_representative_protein_map.tsv",
    "plant_fasta": HOME / "projects/qe_merops_screen/01_input/Evariegata_plant_representative_proteins.faa",
    "insect_fasta": HOME / "projects/qe_merops_screen/01_input/Qerythrinae_insect_representative_proteins.faa",
    "GPup_prior_audit": ROOT / "04_cross_species_integration/results/plant_GPup_proteins_vs_Qerythrinae_best_hits_provenance.tsv",
    "top50_provenance": ROOT / "03_plant_discovery/results/plant_GP_top50_provenance_annotated.tsv",
}

EXPECTED_SHA256 = {
    "plant_fasta": "cb4683b7a4b1b20dfdfe6720b0bc060a3f46ec808f6bc223fb5f55d7ee8160b4",
    "insect_fasta": "fc1b4ac07d4b1c43c97f79f8bd93b39179a65b950cf3ad4e4010b34aef93bdf1",
}

OUTPUTS = {
    "all_best": ROOT / "04_cross_species_integration/results/plant_representative_proteome_vs_Qerythrinae_best_hits.tsv",
    "by_status": ROOT / "04_cross_species_integration/results/plant_representative_proteome_insect_similarity_by_DESeq2_status.tsv",
    "by_class": ROOT / "04_cross_species_integration/results/plant_representative_proteome_insect_similarity_by_analysis_class.tsv",
    "fisher": ROOT / "04_cross_species_integration/results/plant_insect_like_sequence_GPup_enrichment_fisher.tsv",
    "mask": ROOT / "04_cross_species_integration/results/plant_host_reanalysis_proposed_insect_like_sequence_mask.tsv",
    "top50_status": ROOT / "03_plant_discovery/results/plant_GP_top50_provenance_reanalysis_status.tsv",
    "parameters": ROOT / "00_project_control/06_plant_proteome_insect_similarity_parameters.tsv",
    "audit": ROOT / "08_logs/06_plant_proteome_insect_similarity_landscape_audit.txt",
}

DIAMOND = ["micromamba", "run", "-n", "merops", "diamond"]
EVALUE = 1e-5
TARGETS = 25
THREADS = 8
HIT_FIELDS = ["query_id", "subject_id", "percent_identity", "alignment_length", "query_length", "subject_length", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"]
DERIVED = ["query_span_coverage_percent", "subject_span_coverage_percent", "Qerythrinae_similarity_review_class", "insect_like_near_identity_review"]


class AuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean(value: object) -> str:
    return "" if value is None else str(value).strip().lstrip("\ufeff")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AuditError(f"Missing TSV header: {path}")
        fields = [clean(item) for item in reader.fieldnames]
        return fields, [{clean(k): clean(v) for k, v in row.items()} for row in reader]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def index_unique(rows: list[dict[str, str]], key: str, label: str) -> dict[str, dict[str, str]]:
    indexed = {}
    for row in rows:
        value = row.get(key, "")
        if not value or value in indexed:
            raise AuditError(f"Invalid or duplicate {key} in {label}: {value!r}")
        indexed[value] = row
    return indexed


def fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].split(maxsplit=1)[0]
                if not identifier or identifier in ids:
                    raise AuditError(f"Invalid or duplicate FASTA identifier: {identifier!r}")
                ids.add(identifier)
    return ids


def run(command: list[str], stdout_path: Path | None = None) -> str:
    if stdout_path is None:
        completed = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return completed.stdout.strip()
    with stdout_path.open("w", encoding="utf-8", newline="\n") as handle:
        subprocess.run(command, check=True, text=True, stdout=handle, stderr=subprocess.PIPE)
    return ""


def similarity_class(identity: float, qcov: float, scov: float, evalue: float) -> tuple[str, str]:
    if identity >= 90 and qcov >= 80 and scov >= 80:
        return "insect_like_near_identity_review", "TRUE"
    if identity >= 40 and qcov >= 60 and scov >= 60 and evalue <= 1e-10:
        return "insect_like_substantial_homology_review", "FALSE"
    if identity >= 25 and qcov >= 40 and evalue <= EVALUE:
        return "partial_homology_to_Qerythrinae_reference", "FALSE"
    return "weak_or_short_match_to_Qerythrinae_reference", "FALSE"


def parse_hits(path: Path) -> list[dict[str, str]]:
    hits = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for number, values in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(values) != len(HIT_FIELDS):
                raise AuditError(f"Unexpected DIAMOND field count at {path}:{number}")
            hit = dict(zip(HIT_FIELDS, values))
            qlen, slen = int(hit["query_length"]), int(hit["subject_length"])
            qspan = abs(int(hit["query_end"]) - int(hit["query_start"])) + 1
            sspan = abs(int(hit["subject_end"]) - int(hit["subject_start"])) + 1
            qcov, scov = 100 * qspan / qlen, 100 * sspan / slen
            label, near = similarity_class(float(hit["percent_identity"]), qcov, scov, float(hit["evalue"]))
            hit.update({
                "query_span_coverage_percent": format(qcov, ".10g"),
                "subject_span_coverage_percent": format(scov, ".10g"),
                "Qerythrinae_similarity_review_class": label,
                "insect_like_near_identity_review": near,
            })
            hits.append(hit)
    return hits


def best_by_query(hits: list[dict[str, str]], queries: list[str]) -> dict[str, dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for hit in hits:
        grouped.setdefault(hit["query_id"], []).append(hit)
    result = {}
    for query in queries:
        candidates = grouped.get(query, [])
        if candidates:
            result[query] = sorted(
                candidates,
                key=lambda h: (-float(h["bitscore"]), float(h["evalue"])),
            )[0]
        else:
            row = {field: "" for field in HIT_FIELDS + DERIVED}
            row.update({"query_id": query, "Qerythrinae_similarity_review_class": "no_match_to_Qerythrinae_reference_at_evalue_1e-5", "insect_like_near_identity_review": "FALSE"})
            result[query] = row
    return result


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact P using probabilities no greater than observed."""
    population, successes, draws = a + b + c + d, a + c, a + b
    low, high = max(0, draws - (population - successes)), min(draws, successes)
    def logp(x: int) -> float:
        return log_choose(successes, x) + log_choose(population - successes, draws - x) - log_choose(population, draws)
    observed = logp(a)
    accepted = [logp(x) for x in range(low, high + 1) if logp(x) <= observed + 1e-12]
    peak = max(accepted)
    return math.exp(peak) * sum(math.exp(value - peak) for value in accepted)


def main() -> None:
    for name, path in INPUTS.items():
        if not path.is_file():
            raise AuditError(f"Missing input {name}: {path}")
    hashes = {name: sha256(path) for name, path in INPUTS.items()}
    for name, expected in EXPECTED_SHA256.items():
        if hashes[name] != expected:
            raise AuditError(f"Checksum mismatch for {name}")
    if shutil.which("micromamba") is None:
        raise AuditError("micromamba is not available on PATH")
    diamond_version = run([*DIAMOND, "version"])

    _, integrated_rows = read_tsv(INPUTS["plant_integration"])
    _, mapping_rows = read_tsv(INPUTS["plant_protein_map"])
    _, gp_prior_rows = read_tsv(INPUTS["GPup_prior_audit"])
    top_fields, top_rows = read_tsv(INPUTS["top50_provenance"])
    integration = index_unique(integrated_rows, "gene_id", "plant integration")
    mapping = index_unique(mapping_rows, "gene_id", "plant protein map")
    prior_gp = index_unique(gp_prior_rows, "gene_id", "prior GP-up audit")
    top = index_unique(top_rows, "gene_id", "top-50 provenance")
    plant_ids, insect_ids = fasta_ids(INPUTS["plant_fasta"]), fasta_ids(INPUTS["insect_fasta"])
    if len(integration) != 41330 or len(mapping) != 41330 or len(plant_ids) != 41330 or len(insect_ids) != 13455:
        raise AuditError("Representative universe counts differ from validated values")
    if set(integration) != set(mapping) or set(integration) != plant_ids:
        raise AuditError("Plant integrated, mapping, and FASTA identifier sets differ")
    if len(top) != 50:
        raise AuditError("Frozen plant top-50 table does not have 50 rows")

    with tempfile.TemporaryDirectory(prefix="plant_proteome_insect_similarity_") as temporary_name:
        temporary = Path(temporary_name)
        db, raw = temporary / "Qerythrinae_representative", temporary / "plant_vs_insect.raw.tsv"
        run([*DIAMOND, "makedb", "--in", str(INPUTS["insect_fasta"]), "--db", str(db), "--threads", str(THREADS), "--quiet"])
        run([*DIAMOND, "blastp", "--query", str(INPUTS["plant_fasta"]), "--db", str(db), "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "--evalue", str(EVALUE), "--max-target-seqs", str(TARGETS), "--sensitive", "--threads", str(THREADS), "--quiet"], stdout_path=raw)
        best = best_by_query(parse_hits(raw), sorted(plant_ids))

    full_rows = []
    for gene_id in sorted(plant_ids):
        row, protein = integration[gene_id], mapping[gene_id]
        hit = best[gene_id]
        full_rows.append({
            "gene_id": gene_id, "representative_protein_id": protein["protein_id"], "representative_protein_length_aa": protein["aa_length"],
            "baseMean": row["baseMean"], "log2FoldChange": row["log2FoldChange"], "padj": row["padj"], "DESeq2_status": row["DESeq2_status"], "analysis_class": row["analysis_class"],
            **hit,
        })
    full_fields = ["gene_id", "representative_protein_id", "representative_protein_length_aa", "baseMean", "log2FoldChange", "padj", "DESeq2_status", "analysis_class", *HIT_FIELDS, *DERIVED]
    write_tsv(OUTPUTS["all_best"], full_fields, full_rows)

    gp_full = {row["gene_id"]: row for row in full_rows if row["DESeq2_status"] == "GP_up"}
    if len(gp_full) != 4926 or set(gp_full) != set(prior_gp):
        raise AuditError("GP-up membership differs from the preceding validated provenance audit")
    for gene_id, prior in prior_gp.items():
        current = gp_full[gene_id]
        for field in ("subject_id", "Qerythrinae_similarity_review_class"):
            if prior.get(field, "") != current.get(field, ""):
                raise AuditError(f"GP-up audit mismatch for {gene_id} in {field}")

    labels = ["insect_like_near_identity_review", "insect_like_substantial_homology_review", "partial_homology_to_Qerythrinae_reference", "weak_or_short_match_to_Qerythrinae_reference", "no_match_to_Qerythrinae_reference_at_evalue_1e-5"]
    def summarize(rows: list[dict[str, object]], field: str) -> list[dict[str, object]]:
        groups: dict[str, Counter[str]] = {}
        totals: Counter[str] = Counter()
        for row in rows:
            group = str(row[field])
            totals[group] += 1
            groups.setdefault(group, Counter())[str(row["Qerythrinae_similarity_review_class"])] += 1
        output = []
        for group in sorted(groups):
            for label in labels:
                count = groups[group][label]
                output.append({field: group, "provenance_class": label, "candidate_genes": count, "percent_within_group": format(100 * count / totals[group], ".6g")})
        return output
    status_summary = summarize(full_rows, "DESeq2_status")
    class_summary = summarize(full_rows, "analysis_class")
    write_tsv(OUTPUTS["by_status"], ["DESeq2_status", "provenance_class", "candidate_genes", "percent_within_group"], status_summary)
    write_tsv(OUTPUTS["by_class"], ["analysis_class", "provenance_class", "candidate_genes", "percent_within_group"], class_summary)

    tested = [row for row in full_rows if row["DESeq2_status"] in {"GP_up", "UP_up", "non_DE"}]
    gp = [row for row in tested if row["DESeq2_status"] == "GP_up"]
    other = [row for row in tested if row["DESeq2_status"] != "GP_up"]
    a = sum(row["Qerythrinae_similarity_review_class"] == "insect_like_near_identity_review" for row in gp)
    b, c = len(gp) - a, sum(row["Qerythrinae_similarity_review_class"] == "insect_like_near_identity_review" for row in other)
    d = len(other) - c
    odds_ratio = math.inf if b * c == 0 else (a * d) / (b * c)
    write_tsv(OUTPUTS["fisher"], ["comparison", "GPup_near_identity", "GPup_not_near_identity", "other_tested_near_identity", "other_tested_not_near_identity", "odds_ratio", "fisher_two_sided_pvalue"], [{"comparison": "GP_up_vs_UP_up_plus_non_DE_among_representative_protein_genes", "GPup_near_identity": a, "GPup_not_near_identity": b, "other_tested_near_identity": c, "other_tested_not_near_identity": d, "odds_ratio": format(odds_ratio, ".10g"), "fisher_two_sided_pvalue": format(fisher_two_sided(a,b,c,d), ".10g")}])

    mask_rows = []
    for row in full_rows:
        label = row["Qerythrinae_similarity_review_class"]
        mask_rows.append({
            "gene_id": row["gene_id"], "DESeq2_status": row["DESeq2_status"], "Qerythrinae_similarity_review_class": label,
            "proposed_host_reanalysis_mask": "review_exclude_high_confidence_insect_like_sequence" if label == "insect_like_near_identity_review" else "retain_pending_additional_taxonomic_evidence",
            "mask_rationale": "Near identity to the Q. erythrinae representative proteome across >=80% of both sequences" if label == "insect_like_near_identity_review" else "Does not meet the predeclared near-identity review criterion",
        })
    write_tsv(OUTPUTS["mask"], ["gene_id", "DESeq2_status", "Qerythrinae_similarity_review_class", "proposed_host_reanalysis_mask", "mask_rationale"], mask_rows)

    top_rows_out = []
    for gene_id, row in sorted(top.items(), key=lambda item: int(item[1]["discovery_rank"])):
        current = best.get(gene_id)
        output = dict(row)
        if current is None:
            output["full_proteome_provenance_status"] = "no_predicted_protein_transcript_level_review_required"
            output["proposed_host_reanalysis_mask"] = "retain_pending_transcript_level_taxonomic_evidence"
        else:
            output["full_proteome_provenance_status"] = current["Qerythrinae_similarity_review_class"]
            output["proposed_host_reanalysis_mask"] = "review_exclude_high_confidence_insect_like_sequence" if current["Qerythrinae_similarity_review_class"] == "insect_like_near_identity_review" else "retain_pending_additional_taxonomic_evidence"
        top_rows_out.append(output)
    write_tsv(OUTPUTS["top50_status"], [*top_fields, "full_proteome_provenance_status", "proposed_host_reanalysis_mask"], top_rows_out)

    parameters = [
        {"parameter": "DIAMOND_version", "value": diamond_version}, {"parameter": "DIAMOND_mode", "value": "blastp --sensitive"},
        {"parameter": "evalue_limit", "value": str(EVALUE)}, {"parameter": "maximum_target_sequences", "value": str(TARGETS)}, {"parameter": "threads", "value": str(THREADS)},
        {"parameter": "insect_like_near_identity_review", "value": "identity>=90%; query and subject span coverage>=80%"},
        {"parameter": "scope", "value": "all 41,330 plant representative proteins versus the 13,455-protein Q. erythrinae reference"},
    ]
    write_tsv(OUTPUTS["parameters"], ["parameter", "value"], parameters)

    for name, path in INPUTS.items():
        if sha256(path) != hashes[name]:
            raise AuditError(f"Input changed during analysis: {name}")
    all_counts = Counter(str(row["Qerythrinae_similarity_review_class"]) for row in full_rows)
    report = [
        "===== PLANT PROTEOME–INSECT SIMILARITY LANDSCAPE =====", f"Discovery root: {ROOT}", "",
        f"Plant representative proteins screened: {len(full_rows)}", f"Insect representative proteins: {len(insect_ids)}", f"Protein-bearing GP-up genes: {len(gp)}", "",
        "All plant representative proteins by Q. erythrinae similarity review class:",
    ]
    report.extend(f"  {label}: {all_counts[label]}" for label in labels)
    report.extend(["", "GP-up near-identity enrichment among tested representative-protein genes:", f"  GP-up near identity: {a} / {len(gp)}", f"  Other tested near identity: {c} / {len(other)}", f"  Odds ratio: {format(odds_ratio, '.6g')}", f"  Fisher two-sided P: {format(fisher_two_sided(a,b,c,d), '.6g')}", "", "Interpretive boundary:", "  The proposed mask is a conservative review resource, not an applied filtering step.", "  Strong Q. erythrinae similarity does not alone prove taxonomic origin, but broad near identity is incompatible with a convergence claim.", "  DESeq2, enrichment, and the frozen discovery ranks remain unchanged by this script.", "", "Frozen-input SHA-256 verification after output generation: PASS", "Errors detected: 0", ""])
    for path in OUTPUTS.values(): report.append(f"Output: {path}")
    report.append("PLANT PROTEOME–INSECT SIMILARITY LANDSCAPE: PASS")
    OUTPUTS["audit"].parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUTS["audit"].with_name(OUTPUTS["audit"].name + ".tmp")
    temporary.write_text("\n".join(report) + "\n", encoding="utf-8")
    os.replace(temporary, OUTPUTS["audit"])
    print("\n".join(report))


if __name__ == "__main__":
    try:
        main()
    except (AuditError, subprocess.CalledProcessError) as error:
        print("===== PLANT PROTEOME–INSECT SIMILARITY LANDSCAPE =====", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError):
            print(f"ERROR: command failed: {' '.join(error.cmd)}", file=sys.stderr)
            if error.stderr: print(error.stderr.strip(), file=sys.stderr)
        else:
            print(f"ERROR: {error}", file=sys.stderr)
        print("PLANT PROTEOME–INSECT SIMILARITY LANDSCAPE: CHECK REQUIRED", file=sys.stderr)
        raise SystemExit(1)
