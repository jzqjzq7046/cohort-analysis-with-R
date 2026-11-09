"""Independent validation against the paper's 4,332 public sample scores; no plotting."""
from pathlib import Path
import pandas as pd
import numpy as np
root=Path(__file__).resolve().parents[1]
d=root/'literature_reviews/external_rank_redesign_20261029'
a=pd.read_csv(d/'CaseControl_Subject_Data_Table_Jan21.tsv.bz2',sep='\t',index_col=0,low_memory=False)
a.index=[str(i).split('|')[-1].replace('t__','') for i in a.index]
s=pd.read_csv(d/'all_normalized_scores_ZOE_ranks_disease_dataset.tsv',sep='\t',index_col=0)
r=pd.read_excel(root/'literature_reviews/ABC_primary_20261027/Asnicar2026_Supplementary_Tables.xlsx',sheet_name='S5',skiprows=3).set_index('SGB')['HEALTH_ranks']
assert len(r)==661 and r.index.is_unique and set(r.index).issubset(a.index)
x=a.loc[r.index,s.columns].to_numpy(float)/100
assert np.isfinite(x).all() and (x>=0).all() and (x<=1).all()
w=1-2*(r-r.min())/(r.max()-r.min())
rows=[]
for target,y in [('cardio_minusone_to_one',w.to_numpy()@(x>0)),('cardio_minusone_to_one_arcsin',w.to_numpy()@np.arcsin(np.sqrt(x)))]:
 z=s.loc[target].to_numpy(float)
 err=float(np.max(abs(y-z)))
 assert err<1e-10,(target,err)
 rows.append(dict(score=target,N=len(z),max_absolute_error=err,passed=True))
pd.DataFrame(rows).to_csv(d/'public_score_reproduction_audit.csv',index=False)
pd.DataFrame({'SGB':r.index,'HEALTH_ranks':r.values,'fixed_weight':w.values}).to_csv(d/'frozen_health_rank_full_precision.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
