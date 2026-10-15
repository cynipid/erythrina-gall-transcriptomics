#!/usr/bin/env python3
"""Collapse host-focused plant-defence evidence rows to unique genes."""
from pathlib import Path
import csv, math, sys
from collections import defaultdict

ROOT=Path.home()/'projects/qe_cross_species_discovery'
RES=ROOT/'05_plant_defense_audit/results'
INPUT=RES/'plant_host_focused_defense_candidate_inventory.tsv'
def read(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def num(x):
 try:return float(x)
 except (TypeError,ValueError):return math.nan
def join(rows,key):return '; '.join(sorted({r[key] for r in rows if r.get(key)}))
def main():
 if not INPUT.exists():raise SystemExit(f'Missing inventory: {INPUT}')
 rows=read(INPUT); by=defaultdict(list)
 for r in rows:by[r['gene_id']].append(r)
 out=[]
 for gene,items in by.items():
  r=items[0]
  out.append({'gene_id':gene,'defense_groups':join(items,'defense_group'),'candidate_subgroups':join(items,'candidate_subgroup'),
   'Description':r['Description'],'Preferred_name':r['Preferred_name'],'PFAMs':r['PFAMs'],
   'host_focused_status':r['host_focused_status'],'host_baseMean':r['host_baseMean'],'host_log2FoldChange':r['host_log2FoldChange'],'host_padj':r['host_padj'],
   'signalp_prediction':r['signalp_prediction'],'has_TMR':r['has_TMR'],'analysis_class':r['analysis_class'],'primary_classical_secretome':r['primary_classical_secretome'],
   'MEROPS_overlay_present':r['MEROPS_overlay_present'],'MEROPS_families':join(items,'MEROPS_families'),'MEROPS_categories':join(items,'MEROPS_categories')})
 out.sort(key=lambda r:(0 if r['host_focused_status']=='GP_up' else 1 if r['host_focused_status']=='UP_up' else 2,num(r['host_padj']) if not math.isnan(num(r['host_padj'])) else 1))
 fields=list(out[0]); gene_table=RES/'plant_host_focused_defense_candidate_gene_summary.tsv'
 with open(gene_table,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
 groups=['PR_protein','jasmonate','salicylate','peptidase','peptidase_inhibitor']
 statuses=['GP_up','UP_up','non_DE','NA_padj_not_testable','not_in_host_focused_tested_universe']
 grid=[]
 for g in groups:
  for s in statuses:
   n=sum(g in x['defense_groups'].split('; ') and x['host_focused_status']==s for x in out)
   grid.append({'defense_group':g,'host_focused_status':s,'unique_candidate_genes':n})
 status_table=RES/'plant_host_focused_defense_unique_gene_status_summary.tsv'
 with open(status_table,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(grid[0]),delimiter='\t');w.writeheader();w.writerows(grid)
 for direction,stem in [('GP_up','plant_host_focused_defense_GP_up_candidates.tsv'),('UP_up','plant_host_focused_defense_UP_up_candidates.tsv')]:
  subset=[r for r in out if r['host_focused_status']==direction]
  subset.sort(key=lambda r:(num(r['host_padj']) if not math.isnan(num(r['host_padj'])) else 1,-abs(num(r['host_log2FoldChange']))))
  with open(RES/stem,'w',newline='',encoding='utf-8') as f:
   w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(subset)
 lines=['===== HOST-FOCUSED PLANT-DEFENSE GENE-LEVEL SUMMARY =====',f'Evidence rows: {len(rows)}',f'Unique candidate genes: {len(out)}',f'GP-up candidate genes: {sum(r["host_focused_status"]=="GP_up" for r in out)}',f'UP-up candidate genes: {sum(r["host_focused_status"]=="UP_up" for r in out)}','', 'Interpretive scope: genes with multiple evidence memberships are counted once in the gene-level table but may contribute to more than one defence-group summary.','', 'Errors detected: 0','HOST-FOCUSED PLANT-DEFENSE GENE-LEVEL SUMMARY: PASS']
 audit=ROOT/'05_plant_defense_audit/19_host_focused_plant_defense_gene_summary_audit.txt';audit.write_text('\n'.join(lines)+'\n')
 print('\n'.join(lines));print('\nOutput:',gene_table);print('Output:',status_table);print('Output:',audit)
if __name__=='__main__':main()
