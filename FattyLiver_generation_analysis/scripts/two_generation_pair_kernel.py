"""Node-permuted rectangular MRQAP with additive endpoint nuisance spaces.

No dyad-independent P values. Algebraically identical to residual MRQAP but
avoids materializing a millions-of-dyads by dozens-of-columns design matrix.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
os.environ.setdefault('OMP_NUM_THREADS','4')
import argparse
import time
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser()
    for k in ['prefix']:p.add_argument('--'+k,required=True)
    for k in ['na','nb','ka','kb','permutations','df']:p.add_argument('--'+k,type=int,required=True)
    a=p.parse_args();base=Path(a.prefix)
    def read(suffix,dtype,shape):return np.fromfile(str(base)+suffix,dtype=dtype).reshape(shape,order='F')
    E=read('_E.bin','<f8',(a.na,a.nb));X=read('_X.bin','<f8',(a.na,a.nb))
    Qa=read('_QA.bin','<f8',(a.na,a.ka));Qb=read('_QB.bin','<f8',(a.nb,a.kb))
    rows=read('_ROW.bin','<i4',(a.permutations,a.na))-1;cols=read('_COL.bin','<i4',(a.permutations,a.nb))-1
    assert np.all(np.sort(rows,axis=1)==np.arange(a.na)) and np.all(np.sort(cols,axis=1)==np.arange(a.nb))
    assert np.allclose(Qa.T@Qa,np.eye(a.ka)) and np.allclose(Qb.T@Qb,np.eye(a.kb))
    assert np.isclose(np.sum(X*X),1)
    norm2=np.sum(E*E);rs=E.sum(axis=1)/np.sqrt(a.nb);cs=E.sum(axis=0)/np.sqrt(a.na)
    constant=np.sum(E)**2/(a.na*a.nb);out=np.empty(a.permutations);started=time.monotonic()
    for i,(r,c) in enumerate(zip(rows,cols)):
        ep=E[np.ix_(r,c)]
        numerator=np.sum(X*ep)
        # Row and column nuisance projections share the intercept, counted once.
        projection=np.sum((Qa.T@rs[r])**2)+np.sum((Qb.T@cs[c])**2)-constant
        rss=norm2-projection-numerator*numerator
        if rss<=0:raise ValueError('Non-positive permutation residual variance')
        out[i]=numerator/np.sqrt(rss/a.df)
        if (i+1)%1000==0 or i==0 or i+1==len(out):print(f'Pair MRQAP {i+1}/{len(out)}; {time.monotonic()-started:.1f}s',flush=True)
    assert np.isfinite(out).all()
    np.savetxt(str(base)+'_t.csv',out,delimiter=',',header='t_null',comments='',fmt='%.17g')
if __name__=='__main__':main()
