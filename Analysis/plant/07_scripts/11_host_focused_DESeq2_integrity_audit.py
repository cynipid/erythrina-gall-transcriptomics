#!/usr/bin/env python3
"""Read-only integrity audit for the host-focused plant DESeq2 sensitivity branch."""
import csv, hashlib, os, sys
from collections import Counter, defaultdict
from pathlib import Path

H=Path.home(); R=H/"projects/qe_cross_species_discovery"; P=H/"projects/qe_plant_denovo"
base=R/"05_host_focused_DESeq2_sensitivity/results"
paths={
 "mask":R/"04_cross_species_integration/results/plant_confirmed_multilevel_insect_like_sequence_mask.tsv",
 "mixed":P/"08_results_primary/plant_gene_DESeq2_SG_vs_S_all_primary.tsv",
 "host":base/"plant_gene_DESeq2_GP_vs_UP_host_focused_all.tsv",
 "norm":base/"plant_gene_normalized_counts_host_focused.tsv",
 "summary":base/"plant_host_focused_DESeq2_filtering_summary.tsv",
 "exclude":base/"plant_confirmed_multilevel_mask_exclusion_audit.tsv"}
out=base/"plant_host_focused_vs_mixed_status_contingency.tsv"
audit=R/"05_host_focused_DESeq2_sensitivity/logs/11_host_focused_DESeq2_integrity_audit.txt"

def read(path):
    with path.open(encoding="utf-8-sig",newline="") as h:
        x=csv.DictReader(h,delimiter="\t")
        if not x.fieldnames: raise RuntimeError("Missing header: "+str(path))
        return list(x)

def index(rows,label):
    result={}
    for x in rows:
        g=(x.get("gene_id") or "").strip()
        if not g or g in result: raise RuntimeError("Invalid/duplicate gene ID in "+label)
        result[g]=x
    return result

def call(x,suffix=""):
    p=x.get("padj"+suffix) or ""; l=x.get("log2FoldChange"+suffix) or ""
    if p.strip().upper() in {"", "NA", "NAN"}: return "NA_padj_not_testable"
    p,l=float(p),float(l)
    return "GP_up" if p<.01 and l>=1 else "UP_up" if p<.01 and l<=-1 else "non_DE"

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()

def main():
    for p in paths.values():
        if not p.is_file(): raise RuntimeError("Missing input: "+str(p))
    hashes={k:sha(v) for k,v in paths.items()}
    mr,xr,hr,nr,sr,er=(read(paths[k]) for k in ("mask","mixed","host","norm","summary","exclude"))
    mask={x["gene_id"] for x in mr if x.get("host_reanalysis_mask_status")=="review_exclude_confirmed_multilevel_insect_like_sequence"}
    mixed,host,norm,excluded=index(xr,"mixed"),index(hr,"host"),index(nr,"normalized"),index(er,"exclusion")
    errors=[]
    if len(mask)!=3948: errors.append("Mask count is not 3948")
    if set(host)!=set(norm): errors.append("Host result and normalized-count identifiers differ")
    if set(host)&mask: errors.append("Masked genes remain in host results")
    if set(norm)&mask: errors.append("Masked genes remain in host normalized counts")
    if set(excluded)!=mask: errors.append("Mask/exclusion-audit identifiers differ")
    if any(x.get("retained_for_host_focused_analysis")!="FALSE" for x in er): errors.append("A masked gene was retained")
    hs=Counter(call(x) for x in hr); ms=Counter(call(x) for x in xr)
    if dict(hs)!={"GP_up":3179,"UP_up":1396,"non_DE":21008,"NA_padj_not_testable":196}: errors.append("Unexpected host status counts: "+str(dict(hs)))
    if (ms["GP_up"],ms["UP_up"])!=(5593,1467): errors.append("Unexpected mixed DEG counts: "+str(dict(ms)))
    summary={x["metric"]:x["value"] for x in sr}
    for k,v in {"host_focused_GP_up_DEGs":3179,"host_focused_UP_up_DEGs":1396,"genes_after_frozen_low_count_filter":25779}.items():
        if int(float(summary.get(k,"-1")))!=v: errors.append("Summary mismatch: "+k)
    c=defaultdict(Counter)
    for g in set(host)&set(mixed): c[call(mixed[g])][call(host[g])]+=1
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["mixed_primary_status","host_focused_status","genes"],delimiter="\t",lineterminator="\n"); w.writeheader()
        for a in sorted(c):
            for b,n in sorted(c[a].items()): w.writerow({"mixed_primary_status":a,"host_focused_status":b,"genes":n})
    stable=all(sha(paths[k])==hashes[k] for k in paths)
    report=["===== HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT =====",f"Masked genes: {len(mask)}",f"Host tested genes: {len(host)}",f"Host GP-up: {hs['GP_up']}",f"Host UP-up: {hs['UP_up']}",f"Mixed baseline GP-up: {ms['GP_up']}",f"Mixed baseline UP-up: {ms['UP_up']}",f"Masked genes retained in host results: {len(set(host)&mask)}",f"Masked genes retained in host normalized counts: {len(set(norm)&mask)}",f"Input SHA-256 verification after output generation: {'PASS' if stable else 'FAIL'}",f"Errors detected: {len(errors)}"]+[f"ERROR: {e}" for e in errors]+[f"Output: {out}",f"Output: {audit}","HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT: "+("PASS" if not errors and stable else "CHECK REQUIRED")]
    audit.parent.mkdir(parents=True,exist_ok=True); audit.write_text("\n".join(report)+"\n",encoding="utf-8"); print("\n".join(report))
    if errors or not stable: raise SystemExit(1)
if __name__=="__main__":
    try: main()
    except Exception as e:
        print("===== HOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT =====\nERROR: "+str(e)+"\nHOST-FOCUSED PLANT DESEQ2 INTEGRITY AUDIT: CHECK REQUIRED",file=sys.stderr); raise SystemExit(1)
