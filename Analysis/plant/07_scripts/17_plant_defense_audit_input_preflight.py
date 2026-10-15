#!/usr/bin/env python3
"""Read-only input audit for the host-focused plant-defense candidate screen."""
from pathlib import Path
import csv, hashlib, sys

ROOT=Path.home()/"projects"/"qe_cross_species_discovery"
OUT=ROOT/"05_plant_defense_audit"
FILES={
 "full_gene_level_eggnog":Path.home()/"projects/qe_plant_denovo/02_eggnog/processed/E_variegata_gene_level_eggnog.tsv",
 "plant_signalp_deeptmhmm":Path.home()/"projects/qe_plant_denovo/09_signalp_tmhmm/results/signalp6_deeptmhmm_plant_deseq2_primary_41330.tsv",
 "host_focused_deseq2_eggnog":ROOT/"05_host_focused_DESeq2_sensitivity/results/plant_gene_DESeq2_GP_vs_UP_host_focused_eggnog_joined.tsv",
 "plant_merops_overlay":Path.home()/"projects/qe_merops_screen/results/plant_GPup_secretome_MEROPS_gene_category_curated_v2.tsv",
}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
 OUT.mkdir(exist_ok=True); (OUT/'results').mkdir(exist_ok=True)
 rows=[]; errors=[]; lines=['===== HOST-FOCUSED PLANT-DEFENSE INPUT PREFLIGHT =====']
 for label,p in FILES.items():
  if not p.exists(): errors.append(f'Missing: {p}');continue
  x=read(p); cols=list(x[0]) if x else []
  if 'gene_id' not in cols:errors.append(f'{label} lacks gene_id')
  ids=[r.get('gene_id','') for r in x]
  lines += [f'\n{label}:',f'  Path: {p}',f'  Data rows: {len(x)}',f'  Unique gene IDs: {len(set(ids))}',f'  SHA-256: {sha(p)}']
  rows.append({'resource':label,'path':str(p),'data_rows':len(x),'unique_gene_ids':len(set(ids)),'sha256':sha(p)})
 # frozen expected universes
 expected={'full_gene_level_eggnog':23379,'plant_signalp_deeptmhmm':41330,'host_focused_deseq2_eggnog':25779}
 for k,n in expected.items():
  q=next((r for r in rows if r['resource']==k),None)
  if q and q['data_rows']!=n:errors.append(f'{k} has {q["data_rows"]} rows; expected {n}')
 lines += ['', 'Interpretive scope:', 'Candidate families will be defined from the full annotated gene universe, then overlaid with host-focused DESeq2 status.', 'SignalP/DeepTMHMM covers all representative proteins. MEROPS is an existing overlay restricted to the original GP-up classical-secretome universe.', '',f'Errors detected: {len(errors)}']
 lines += [f'ERROR: {e}' for e in errors]
 lines.append('HOST-FOCUSED PLANT-DEFENSE INPUT PREFLIGHT: '+('PASS' if not errors else 'CHECK REQUIRED'))
 audit=OUT/'17_plant_defense_input_preflight.txt';audit.write_text('\n'.join(lines)+'\n')
 with open(OUT/'plant_defense_input_resource_manifest.tsv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['resource','path','data_rows','unique_gene_ids','sha256'],delimiter='\t');w.writeheader();w.writerows(rows)
 print('\n'.join(lines));print('\nOutput:',audit)
 return 1 if errors else 0
if __name__=='__main__':sys.exit(main())
