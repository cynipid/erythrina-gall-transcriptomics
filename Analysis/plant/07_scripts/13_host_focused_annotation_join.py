#!/usr/bin/env python3
"""Attach frozen plant eggNOG fields to host-focused DESeq2 results."""
import csv, hashlib, os, sys
from pathlib import Path

H=Path.home()
ROOT=H/"projects/qe_cross_species_discovery"
HOST=ROOT/"05_host_focused_DESeq2_sensitivity/results/plant_gene_DESeq2_GP_vs_UP_host_focused_all.tsv"
DONOR=H/"projects/qe_plant_denovo/08_results_primary/plant_gene_DESeq2_SG_vs_S_all_primary_eggnog_joined.tsv"
OUT=ROOT/"05_host_focused_DESeq2_sensitivity/results/plant_gene_DESeq2_GP_vs_UP_host_focused_eggnog_joined.tsv"
AUDIT=ROOT/"05_host_focused_DESeq2_sensitivity/logs/13_host_focused_annotation_join_audit.txt"

def read(p):
 with p.open(encoding="utf-8-sig",newline="") as f:
  r=csv.DictReader(f,delimiter="\t")
  if not r.fieldnames: raise ValueError("Missing header: "+str(p))
  return r.fieldnames,[{k:(v or "").strip() for k,v in x.items()} for x in r]
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""): h.update(b)
 return h.hexdigest()
def main():
 if not HOST.is_file() or not DONOR.is_file(): raise ValueError("Missing host result or frozen annotation donor")
 before={HOST:sha(HOST),DONOR:sha(DONOR)}
 hf,hr=read(HOST); df,dr=read(DONOR)
 if "gene_id" not in hf or "gene_id" not in df: raise ValueError("Both inputs must have gene_id")
 d={}
 for x in dr:
  g=x["gene_id"]
  if not g or g in d: raise ValueError("Invalid/duplicate donor gene ID: "+g)
  d[g]=x
 # Only annotation columns are transferred. DESeq2 columns remain host-focused.
 blocked={"gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","DESeq2_status","DESeq2_status_host_focused"}
 annotation=[x for x in df if x not in blocked]
 rows=[]; matches=0
 for x in hr:
  y=dict(x); donor=d.get(x["gene_id"])
  y["eggnog_annotation_match"]= "matched" if donor else "unannotated_or_absent_from_frozen_donor"
  if donor:
   matches+=1
   for k in annotation: y[k]=donor.get(k,"")
  else:
   for k in annotation: y[k]=""
  rows.append(y)
 OUT.parent.mkdir(parents=True,exist_ok=True)
 with OUT.open("w",encoding="utf-8",newline="") as f:
  w=csv.DictWriter(f,fieldnames=hf+["eggnog_annotation_match"]+annotation,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
 if sha(HOST)!=before[HOST] or sha(DONOR)!=before[DONOR]: raise ValueError("Input changed during join")
 report=["===== HOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN =====",f"Host-focused DESeq2 genes: {len(hr)}",f"Frozen eggNOG donor rows: {len(dr)}",f"Host-focused genes with annotation donor match: {matches}",f"Host-focused genes without donor match: {len(hr)-matches}","Errors detected: 0",f"Output: {OUT}",f"Output: {AUDIT}","HOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN: PASS"]
 AUDIT.parent.mkdir(parents=True,exist_ok=True);AUDIT.write_text("\n".join(report)+"\n",encoding="utf-8");print("\n".join(report))
if __name__=="__main__":
 try:main()
 except Exception as e:
  print("===== HOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN =====\nERROR: "+str(e)+"\nHOST-FOCUSED PLANT EGGNOG ANNOTATION JOIN: CHECK REQUIRED",file=sys.stderr);raise SystemExit(1)
