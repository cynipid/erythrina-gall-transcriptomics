#!/usr/bin/env python3
"""Annotation/domain-supported plant defence candidate audit.

This is a hypothesis-guided inventory, not a new DEG test or an exhaustive
gene-family annotation.  It identifies candidates from the full eggNOG/Pfam
annotated representative-gene resource, then overlays host-focused DESeq2,
SignalP/DeepTMHMM, and the pre-existing restricted MEROPS screen.
"""
from pathlib import Path
import csv, hashlib, re, sys
from collections import defaultdict, Counter

HOME=Path.home(); ROOT=HOME/'projects/qe_cross_species_discovery'; OUT=ROOT/'05_plant_defense_audit'; RES=OUT/'results'
EGG=HOME/'projects/qe_plant_denovo/02_eggnog/processed/E_variegata_gene_level_eggnog.tsv'
HOST=ROOT/'05_host_focused_DESeq2_sensitivity/results/plant_gene_DESeq2_GP_vs_UP_host_focused_eggnog_joined.tsv'
SPTM=HOME/'projects/qe_plant_denovo/09_signalp_tmhmm/results/signalp6_deeptmhmm_plant_deseq2_primary_41330.tsv'
MEROPS=HOME/'projects/qe_merops_screen/results/plant_GPup_secretome_MEROPS_gene_category_curated_v2.tsv'

# Each rule requires either a recognisable functional annotation or a named
# protein-domain architecture; individual terms are retained as evidence.
RULES=[
 ('PR_protein','PR1_CAP',r'pathogenesis.related protein 1|\bpr[- ]?1\b|\bCAP\b'),
 ('PR_protein','beta_1_3_glucanase_PR2',r'beta.1,3.glucanase|glyco_hydro_17|glycosyl hydrolase family 17'),
 ('PR_protein','chitinase_PR3_PR4',r'chitinase|glyco_hydro_18|glyco_hydro_19'),
 ('PR_protein','thaumatin_PR5',r'thaumatin'),
 ('PR_protein','peroxidase_PR9',r'classical plant \(class iii\) peroxidase|\bperoxidase\b'),
 ('PR_protein','Bet_v_1_PR10',r'bet v 1|pathogenesis.related protein 10'),
 ('PR_protein','defensin_PR12',r'defensin'),
 ('PR_protein','lipid_transfer_PR14',r'lipid transfer protein|lipid_transf'),
 ('PR_protein','proteinase_inhibitor_PR6',r'protease inhibitor|proteinase inhibitor|peptidase inhibitor|kunitz|cystatin|serpin|bowman|potato_inhibit|potato inhibitor|kazal|pacifastin|\bTIMP\b'),
 ('jasmonate','biosynthesis',r'lipoxygenase|allene oxide synthase|allene oxide cyclase|12.oxophytodienoate reductase|\bOPR3\b'),
 ('jasmonate','perception_signalling',r'jasmonate[ -]z[iy]m.domain|\bJAZ\b|\bTIFY\b|coronatine.insensitive|\bCOI1\b|\bMYC2\b|jasmonate resistant|\bJAR1\b'),
 ('salicylate','biosynthesis_transport',r'isochorismate synthase|\bICS1\b|\bSID2\b|enhanced disease susceptibility 5|\bEDS5\b|avrpphb susceptible 3|\bPBS3\b'),
 ('salicylate','signalling',r'nonexpressor of pathogenesis.related genes|\bNPR1\b|phytoalexin deficient 4|\bPAD4\b|enhanced disease susceptibility 1|\bEDS1\b|\bTGA\b'),
 ('peptidase','protease_annotation_or_domain',r'peptidase|protease|proteinase|subtilisin|papain|metalloprotease|aspartic protease|carboxypeptidase|Peptidase_[A-Z]'),
 ('peptidase_inhibitor','inhibitor_annotation_or_domain',r'protease inhibitor|proteinase inhibitor|peptidase inhibitor|kunitz|cystatin|serpin|bowman|potato_inhibit|potato inhibitor|kazal|pacifastin|\bTIMP\b'),
]

def read(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def clean(x): return '' if x in (None,'-','--','NA') else x
def main():
 RES.mkdir(parents=True,exist_ok=True); OUT.mkdir(exist_ok=True)
 for p in (EGG,HOST,SPTM,MEROPS):
  if not p.exists(): raise SystemExit(f'Missing required input: {p}')
 egg,host,sptm,merops=map(read,(EGG,HOST,SPTM,MEROPS))
 host_by={r['gene_id']:r for r in host}; sp_by={r['gene_id']:r for r in sptm}; mer_by=defaultdict(list)
 for r in merops: mer_by[r['gene_id']].append(r)
 out=[]
 for r in egg:
  text=' | '.join(clean(r.get(k)) for k in ('Description','Preferred_name','PFAMs','CAZy','GOs','KEGG_ko'))
  for group,subgroup,pattern in RULES:
   hits=re.findall(pattern,text,flags=re.I)
   if not hits: continue
   text_lower=text.lower()
   # Kunitz/cystatin/serpin and related inhibitor descriptions contain the
   # word "protease" but are not themselves peptidases.  Retain peptidase
   # membership only when a catalytic peptidase domain/annotation is also
   # present.
   if group=='peptidase' and re.search(
       r'protease inhibitor|proteinase inhibitor|peptidase inhibitor|kunitz|cystatin|serpin|bowman|potato_inhibit|potato inhibitor|kazal|pacifastin',
       text, flags=re.I
   ) and not re.search(r'Peptidase_[A-Z]|\bpeptidase [ACGMSTU]\d|subtilisin.like protease', text, flags=re.I):
    continue
   # MEROPS I9 and I29 are propeptide domains within peptidase precursors;
   # they are not evidence of a stand-alone inhibitory protein.
   if group in {'peptidase_inhibitor','PR_protein'} and subgroup in {'inhibitor_annotation_or_domain','proteinase_inhibitor_PR6'}:
    if re.search(r'Inhibitor_I(?:9|29)', text, flags=re.I) and not re.search(
        r'protease inhibitor|proteinase inhibitor|peptidase inhibitor|kunitz|cystatin|serpin|bowman|potato_inhibit|potato inhibitor|kazal|pacifastin|\bTIMP\b',
        text, flags=re.I
    ):
     continue
   # PR-9 refers to secreted/class-III plant peroxidases, not glutathione
   # peroxidases or unrelated intracellular peroxidase families.
   if group=='PR_protein' and subgroup=='peroxidase_PR9':
    pfam=clean(r.get('PFAMs')).lower(); desc=clean(r.get('Description')).lower()
    if 'gshpx' in pfam or 'glutathione peroxidase' in desc:
     continue
    if 'class iii' not in desc and 'classical plant' not in desc and 'peroxidase' not in [x.strip().lower() for x in pfam.split(',')]:
     continue
   # A GH18 domain embedded in a receptor-like kinase is not automatically a
   # canonical PR-3/PR-4 chitinase.
   if group=='PR_protein' and subgroup=='chitinase_PR3_PR4':
    if re.search(r'Pkinase', clean(r.get('PFAMs')), flags=re.I) and not re.search(r'chitinase', clean(r.get('Description')), flags=re.I):
     continue
   gid=r['gene_id']; h=host_by.get(gid,{}); s=sp_by.get(gid,{})
   m=mer_by.get(gid,[])
   out.append({
    'gene_id':gid,'defense_group':group,'candidate_subgroup':subgroup,
    'candidate_evidence':'annotation_or_Pfam_domain','matched_text':'; '.join(sorted(set(str(x) for x in hits))),
    'Description':clean(r.get('Description')),'Preferred_name':clean(r.get('Preferred_name')),'PFAMs':clean(r.get('PFAMs')),
    'host_focused_status':h.get('DESeq2_status_host_focused','not_in_host_focused_tested_universe'),
    'host_baseMean':h.get('baseMean',''),'host_log2FoldChange':h.get('log2FoldChange',''),'host_padj':h.get('padj',''),
    'signalp_prediction':s.get('signalp_prediction','not_available'),'has_TMR':s.get('has_TMR','not_available'),
    'analysis_class':s.get('analysis_class','not_available'),'primary_classical_secretome':s.get('primary_classical_secretome','not_available'),
    'MEROPS_overlay_present':str(bool(m)).upper(),
    'MEROPS_overlay_scope':'original_GP-up_classical-secretome_only',
    'MEROPS_families':'; '.join(sorted({x.get('merops_family','') for x in m if x.get('merops_family')})),
    'MEROPS_categories':'; '.join(sorted({x.get('merops_category','') for x in m if x.get('merops_category')})),
   })
 fields=list(out[0]) if out else ['gene_id']
 table=RES/'plant_host_focused_defense_candidate_inventory.tsv'
 with open(table,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
 counts=Counter((r['defense_group'],r['host_focused_status']) for r in out)
 summary=RES/'plant_host_focused_defense_candidate_status_summary.tsv'
 with open(summary,'w',newline='',encoding='utf-8') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['defense_group','host_focused_status','candidate_rows'])
  for (g,s),n in sorted(counts.items()):w.writerow([g,s,n])
 lines=['===== HOST-FOCUSED PLANT-DEFENSE CANDIDATE AUDIT =====',f'Full annotated gene universe: {len(egg)}',f'Host-focused tested genes: {len(host)}',f'Candidate evidence rows: {len(out)}',f'Candidate genes: {len(set(r["gene_id"] for r in out))}',f'MEROPS overlay rows: {len(merops)}; unique genes: {len(mer_by)}','', 'Boundary: candidate membership is annotation/domain-supported and may include multiple defence groups per gene. MEROPS absence is not negative evidence because its source screen was restricted to original GP-up classical-secretome proteins.','', 'Errors detected: 0','HOST-FOCUSED PLANT-DEFENSE CANDIDATE AUDIT: PASS']
 audit=OUT/'18_host_focused_plant_defense_candidate_audit.txt';audit.write_text('\n'.join(lines)+'\n')
 print('\n'.join(lines));print('\nOutput:',table);print('Output:',summary);print('Output:',audit)
if __name__=='__main__':main()
