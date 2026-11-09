"""Freeze reusable, previously audited numerical utilities before updating entries."""
from pathlib import Path
import re, shutil, hashlib, json
ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'code_history' / 'before_two_generation_20261024'
DEST.mkdir(parents=True, exist_ok=True)
paths = list((ROOT/'scripts').glob('[0-9][0-9]_*.R')) + list((ROOT/'R').glob('*.R')) + [ROOT/'config.R', ROOT/'AGENTS.md']
for p in paths:
    to = DEST / p.relative_to(ROOT)
    to.parent.mkdir(parents=True,exist_ok=True)
    if not to.exists(): shutil.copy2(p,to)
groups = {
 'scripts/21_ast_overall_only.R':['ov_bind','ov_design','ov_extract','ov_expand','ov_species_label','ov_export'],
 'R/cross_generation_clinical.R':['fl_cg_make_comparison','fl_cg_describe','fl_cg_fit'],
 'R/microbiome_steatosis_diversity.R':['fl_dv_permutations','fl_dv_fast_setup','fl_dv_fast_null'],
 'R/microbiome_similarity_adjusted.R':['fl_msa_node_permutations','fl_msa_lm_check','fl_msa_kernel'],
}
chunks=['# Frozen audited numerical helpers. Active designs are defined in reanalysis_20261024.R.']
for rel,names in groups.items():
    s=(DEST/rel).read_text(encoding='utf-8-sig')
    starts=list(re.finditer(r'^([A-Za-z][\w.]*)\s*<-\s*function\s*\(',s,re.M))
    for name in names:
        i=next(i for i,m in enumerate(starts) if m[1]==name)
        end=starts[i+1].start() if i+1<len(starts) else len(s)
        chunks.append(s[starts[i].start():end].strip())
(ROOT/'R/reanalysis_numerical_helpers.R').write_text('\n\n'.join(chunks)+'\n',encoding='utf-8')
print('Frozen helpers and original source backup:',DEST)
