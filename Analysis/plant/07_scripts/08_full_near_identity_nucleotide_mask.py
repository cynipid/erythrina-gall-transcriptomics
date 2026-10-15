#!/usr/bin/env python3
"""Build a conservative whole-proteome insect-like sequence review mask.

The script confirms, at nucleotide level, every plant representative protein
that met the predeclared protein near-identity criterion. The resulting mask
is a provenance resource for a later host-focused DESeq2 sensitivity analysis;
it does not alter current counts, DESeq2 output, or frozen rankings.
"""

from __future__ import annotations

import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


HOME = Path.home()
ROOT = HOME / "projects/qe_cross_species_discovery"
INPUTS = {
    "protein_mask": ROOT / "04_cross_species_integration/results/plant_host_reanalysis_proposed_insect_like_sequence_mask.tsv",
    "prior_CDS_audit": ROOT / "04_cross_species_integration/results/plant_GPup_insectlike_plus_top50_CDS_provenance_best_by_gene.tsv",
    "trinity_assembly": HOME / "projects/qe_plant_busco_embryophyta/input/Trinity_plant.fasta",
    "insect_CDS": HOME / "projects/qe_genome_guided/01_reference/complete.cds.fna",
}
OUTPUTS = {
    "selected_transcripts": ROOT / "03_plant_discovery/results/plant_all_protein_near_identity_transcripts_for_CDS_audit.fasta",
    "transcript_mappings": ROOT / "04_cross_species_integration/results/plant_all_protein_near_identity_transcripts_vs_Qerythrinae_CDS.tsv",
    "gene_best": ROOT / "04_cross_species_integration/results/plant_all_protein_near_identity_CDS_best_by_gene.tsv",
    "confirmed_mask": ROOT / "04_cross_species_integration/results/plant_confirmed_multilevel_insect_like_sequence_mask.tsv",
    "status_summary": ROOT / "04_cross_species_integration/results/plant_confirmed_multilevel_insect_like_sequence_summary_by_DESeq2_status.tsv",
    "parameters": ROOT / "00_project_control/08_full_near_identity_nucleotide_mask_parameters.tsv",
    "audit": ROOT / "08_logs/08_full_near_identity_nucleotide_mask_audit.txt",
}
MINIMAP = ["micromamba", "run", "-n", "qe_provenance", "minimap2"]
THREADS = 8
TRINITY_GENE = re.compile(r"^(TRINITY_DN\d+_c\d+_g\d+)_i\d+$")
PAF = ["query_transcript_id", "query_length_nt", "query_start", "query_end", "strand", "target_CDS_id", "target_length_nt", "target_start", "target_end", "matching_bases", "alignment_block_length", "mapping_quality"]
DERIVED = ["percent_identity", "query_span_coverage_percent", "target_span_coverage_percent", "nucleotide_provenance_review_class", "high_identity_Qerythrinae_CDS_review"]


class AuditError(RuntimeError): pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None: raise AuditError(f"Missing TSV header: {path}")
        fields = [item.strip().lstrip("\ufeff") for item in reader.fieldnames]
        return fields, [{k.strip().lstrip("\ufeff"): (v or "").strip() for k, v in row.items()} for row in reader]


def unique(rows, key, label):
    result = {}
    for row in rows:
        value = row.get(key, "")
        if not value or value in result: raise AuditError(f"Invalid/duplicate {key} in {label}: {value!r}")
        result[value] = row
    return result


def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader(); writer.writerows(rows)
    os.replace(temporary, path)


def gene_from_transcript(identifier: str) -> str:
    match = TRINITY_GENE.fullmatch(identifier)
    if match is None: raise AuditError(f"Unrecognised Trinity transcript identifier: {identifier}")
    return match.group(1)


def extract(assembly: Path, selected: set[str], output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    found = defaultdict(list); tid = gene = None; parts = []
    def store(handle):
        if tid is None or gene not in selected: return
        seq = "".join(parts).upper()
        if not seq: raise AuditError(f"Empty selected transcript: {tid}")
        if set(seq) - set("ACGTRYSWKMBDHVN"): raise AuditError(f"Invalid nucleotide character in {tid}")
        handle.write(f">{tid} gene_id={gene}\n")
        for start in range(0, len(seq), 80): handle.write(seq[start:start+80] + "\n")
        found[gene].append(tid)
    with temporary.open("w", encoding="utf-8", newline="\n") as out:
        with assembly.open("r", encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                line = raw.strip()
                if not line: continue
                if line.startswith(">"):
                    store(out); tid = line[1:].split(maxsplit=1)[0]; gene = gene_from_transcript(tid); parts = []
                else:
                    if tid is None: raise AuditError("Sequence occurred before FASTA header")
                    parts.append(line)
            store(out)
    os.replace(temporary, output)
    return dict(found)


def run(command, stdout=None):
    if stdout is None:
        result = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout.strip()
    with stdout.open("w", encoding="utf-8", newline="\n") as handle:
        subprocess.run(command, check=True, text=True, stdout=handle, stderr=subprocess.PIPE)
    return ""


def classification(identity, qcov, length):
    if identity >= 95 and qcov >= 60 and length >= 300: return "high_identity_Qerythrinae_CDS_review", "TRUE"
    if identity >= 90 and qcov >= 40 and length >= 200: return "substantial_Qerythrinae_CDS_similarity_review", "FALSE"
    if identity >= 70 and length >= 150: return "partial_Qerythrinae_CDS_similarity", "FALSE"
    return "weak_or_short_CDS_match", "FALSE"


def parse_paf(path: Path):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for n, values in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(values) < 12: raise AuditError(f"Malformed PAF at {path}:{n}")
            row = dict(zip(PAF, values[:12])); qlen, tlen = int(row["query_length_nt"]), int(row["target_length_nt"])
            qspan, tspan = int(row["query_end"])-int(row["query_start"]), int(row["target_end"])-int(row["target_start"])
            matches, length = int(row["matching_bases"]), int(row["alignment_block_length"])
            identity, qcov, tcov = 100*matches/length, 100*qspan/qlen, 100*tspan/tlen
            label, strict = classification(identity, qcov, length)
            row.update({"percent_identity":format(identity,".10g"), "query_span_coverage_percent":format(qcov,".10g"), "target_span_coverage_percent":format(tcov,".10g"), "nucleotide_provenance_review_class":label, "high_identity_Qerythrinae_CDS_review":strict})
            rows.append(row)
    return rows


def gene_best(mappings, selected):
    by_gene = defaultdict(list)
    for row in mappings: by_gene[gene_from_transcript(row["query_transcript_id"])].append(row)
    result = {}
    for gene in selected:
        candidates = by_gene.get(gene, [])
        if not candidates:
            result[gene] = {"best_query_transcript_id":"", "target_CDS_id":"", "percent_identity":"", "query_span_coverage_percent":"", "target_span_coverage_percent":"", "alignment_block_length":"", "mapping_quality":"", "nucleotide_provenance_review_class":"no_CDS_match_reported_by_minimap2", "high_identity_Qerythrinae_CDS_review":"FALSE"}
            continue
        best = sorted(candidates, key=lambda x:(x["high_identity_Qerythrinae_CDS_review"] != "TRUE", -float(x["query_span_coverage_percent"]), -float(x["percent_identity"]), -int(x["alignment_block_length"]), -int(x["mapping_quality"])))[0]
        result[gene] = {"best_query_transcript_id":best["query_transcript_id"], "target_CDS_id":best["target_CDS_id"], "percent_identity":best["percent_identity"], "query_span_coverage_percent":best["query_span_coverage_percent"], "target_span_coverage_percent":best["target_span_coverage_percent"], "alignment_block_length":best["alignment_block_length"], "mapping_quality":best["mapping_quality"], "nucleotide_provenance_review_class":best["nucleotide_provenance_review_class"], "high_identity_Qerythrinae_CDS_review":best["high_identity_Qerythrinae_CDS_review"]}
    return result


def main():
    for name, path in INPUTS.items():
        if not path.is_file(): raise AuditError(f"Missing input {name}: {path}")
    hashes = {name:sha256(path) for name,path in INPUTS.items()}
    if shutil.which("micromamba") is None: raise AuditError("micromamba is not available on PATH")
    minimap_version = run([*MINIMAP,"--version"])
    _, mask_rows = read_tsv(INPUTS["protein_mask"]); _, prior_rows = read_tsv(INPUTS["prior_CDS_audit"])
    mask, prior = unique(mask_rows,"gene_id","protein mask"), unique(prior_rows,"gene_id","prior CDS audit")
    if len(mask) != 41330: raise AuditError(f"Expected 41,330 protein-mask rows; found {len(mask)}")
    selected = {gene for gene,row in mask.items() if row["Qerythrinae_similarity_review_class"] == "insect_like_near_identity_review"}
    if len(selected) != 4213: raise AuditError(f"Expected 4,213 protein near-identity genes; found {len(selected)}")
    found = extract(INPUTS["trinity_assembly"], selected, OUTPUTS["selected_transcripts"])
    missing = sorted(selected-set(found))
    if missing: raise AuditError("Near-identity genes missing from Trinity assembly: " + ", ".join(missing))
    with tempfile.TemporaryDirectory(prefix="full_near_identity_CDS_") as d:
        paf = Path(d)/"selected.paf"
        run([*MINIMAP,"-x","asm5","--secondary=no","-t",str(THREADS),str(INPUTS["insect_CDS"]),str(OUTPUTS["selected_transcripts"])], stdout=paf)
        mappings = parse_paf(paf)
    write_tsv(OUTPUTS["transcript_mappings"], [*PAF,*DERIVED], mappings)
    best = gene_best(mappings, selected)
    for gene,row in prior.items():
        if gene not in selected: continue
        current = best[gene]
        for field in ("best_query_transcript_id", "target_CDS_id", "nucleotide_provenance_review_class"):
            if row.get(field,"") != current.get(field,""):
                raise AuditError(f"Prior CDS audit mismatch for {gene} in {field}")
    best_rows = []
    for gene in sorted(selected):
        source = mask[gene]; info = best[gene]
        confirmed = info["nucleotide_provenance_review_class"] in {"high_identity_Qerythrinae_CDS_review", "substantial_Qerythrinae_CDS_similarity_review"}
        best_rows.append({"gene_id":gene,"DESeq2_status":source["DESeq2_status"],"protein_provenance_review_class":source["Qerythrinae_similarity_review_class"],"Trinity_transcripts_extracted":len(found[gene]),**info,"multilevel_provenance_status":"confirmed_multilevel_insect_like_sequence" if confirmed else "protein_near_identity_but_nucleotide_evidence_insufficient_for_confirmed_mask","confirmed_for_host_reanalysis_exclusion":"TRUE" if confirmed else "FALSE"})
    best_fields = ["gene_id","DESeq2_status","protein_provenance_review_class","Trinity_transcripts_extracted","best_query_transcript_id","target_CDS_id","percent_identity","query_span_coverage_percent","target_span_coverage_percent","alignment_block_length","mapping_quality","nucleotide_provenance_review_class","high_identity_Qerythrinae_CDS_review","multilevel_provenance_status","confirmed_for_host_reanalysis_exclusion"]
    write_tsv(OUTPUTS["gene_best"], best_fields, best_rows)
    confirmed = {row["gene_id"] for row in best_rows if row["confirmed_for_host_reanalysis_exclusion"] == "TRUE"}
    full_mask = []
    for gene, row in sorted(mask.items()):
        if gene in confirmed:
            status = "review_exclude_confirmed_multilevel_insect_like_sequence"
            reason = "Protein near-identity plus high/substantial nucleotide concordance to the Q. erythrinae CDS reference"
        elif gene in selected:
            status = "retain_provisional_nucleotide_evidence_insufficient"
            reason = "Protein near-identity but no high/substantial nucleotide concordance"
        else:
            status = "retain_not_in_protein_near_identity_set"
            reason = "Does not meet predeclared protein near-identity criterion"
        full_mask.append({"gene_id":gene,"DESeq2_status":row["DESeq2_status"],"protein_provenance_review_class":row["Qerythrinae_similarity_review_class"],"host_reanalysis_mask_status":status,"mask_rationale":reason})
    write_tsv(OUTPUTS["confirmed_mask"],["gene_id","DESeq2_status","protein_provenance_review_class","host_reanalysis_mask_status","mask_rationale"],full_mask)
    labels = ["confirmed_multilevel_insect_like_sequence","protein_near_identity_but_nucleotide_evidence_insufficient_for_confirmed_mask"]
    counter = defaultdict(Counter); totals = Counter()
    for row in best_rows:
        totals[row["DESeq2_status"]] += 1; counter[row["DESeq2_status"]][row["multilevel_provenance_status"]] += 1
    summary = []
    for status in sorted(totals):
        for label in labels:
            count=counter[status][label]
            summary.append({"DESeq2_status":status,"multilevel_provenance_status":label,"candidate_genes":count,"percent_within_status":format(100*count/totals[status],".6g")})
    write_tsv(OUTPUTS["status_summary"],["DESeq2_status","multilevel_provenance_status","candidate_genes","percent_within_status"],summary)
    params=[{"parameter":"minimap2_version","value":minimap_version},{"parameter":"preset","value":"asm5"},{"parameter":"secondary_alignments","value":"disabled"},{"parameter":"threads","value":str(THREADS)},{"parameter":"protein_near_identity_genes","value":"4213"},{"parameter":"confirmed_multilevel_insect_like_sequence","value":"protein near identity plus nucleotide class high_identity or substantial"},{"parameter":"interpretive_boundary","value":"Mask is a conservative sequence-provenance resource; it does not identify the physical RNA source or prove a technical error"}]
    params.extend({"parameter":"SHA256_"+name,"value":value} for name,value in sorted(hashes.items()))
    write_tsv(OUTPUTS["parameters"],["parameter","value"],params)
    for name,path in INPUTS.items():
        if sha256(path) != hashes[name]: raise AuditError(f"Input changed during analysis: {name}")
    classes=Counter(row["nucleotide_provenance_review_class"] for row in best_rows); status_counts=Counter(row["DESeq2_status"] for row in best_rows); confirmed_status=Counter(row["DESeq2_status"] for row in best_rows if row["confirmed_for_host_reanalysis_exclusion"]=="TRUE")
    report=["===== FULL PLANT PROTEIN-NEAR-IDENTITY NUCLEOTIDE MASK =====",f"Discovery root: {ROOT}","",f"Protein near-identity genes audited: {len(best_rows)}",f"Selected Trinity transcripts extracted: {sum(len(v) for v in found.values())}",f"Confirmed multilevel insect-like sequences: {len(confirmed)}","","Nucleotide provenance classes:"]
    for label in ["high_identity_Qerythrinae_CDS_review","substantial_Qerythrinae_CDS_similarity_review","partial_Qerythrinae_CDS_similarity","weak_or_short_CDS_match","no_CDS_match_reported_by_minimap2"]: report.append(f"  {label}: {classes[label]}")
    report.extend(["","Confirmed multilevel insect-like sequences by DESeq2 status:"])
    for status in sorted(status_counts): report.append(f"  {status}: {confirmed_status[status]} / {status_counts[status]}")
    report.extend(["","Interpretive boundary:","  The confirmed mask is a conservative, multi-level sequence-provenance resource for a later host-focused sensitivity analysis.","  It does not prove the physical source of RNA or establish dissection or laboratory error.","  No DESeq2 result or frozen ranking has been changed.","","Input SHA-256 verification after output generation: PASS","Errors detected: 0",""])
    for path in OUTPUTS.values(): report.append(f"Output: {path}")
    report.append("FULL PLANT PROTEIN-NEAR-IDENTITY NUCLEOTIDE MASK: PASS")
    OUTPUTS["audit"].parent.mkdir(parents=True,exist_ok=True); temporary=OUTPUTS["audit"].with_name(OUTPUTS["audit"].name+".tmp"); temporary.write_text("\n".join(report)+"\n",encoding="utf-8"); os.replace(temporary,OUTPUTS["audit"]); print("\n".join(report))


if __name__ == "__main__":
    try: main()
    except (AuditError, subprocess.CalledProcessError) as error:
        print("===== FULL PLANT PROTEIN-NEAR-IDENTITY NUCLEOTIDE MASK =====",file=sys.stderr)
        if isinstance(error,subprocess.CalledProcessError):
            print(f"ERROR: command failed: {' '.join(error.cmd)}",file=sys.stderr)
            if error.stderr: print(error.stderr.strip(),file=sys.stderr)
        else: print(f"ERROR: {error}",file=sys.stderr)
        print("FULL PLANT PROTEIN-NEAR-IDENTITY NUCLEOTIDE MASK: CHECK REQUIRED",file=sys.stderr); raise SystemExit(1)
