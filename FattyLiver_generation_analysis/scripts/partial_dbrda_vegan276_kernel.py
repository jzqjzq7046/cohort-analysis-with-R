"""Exact reduced-model partial dbRDA permutations for vegan 2.7-6.

vegan permutes conditions with the reduced response, then re-residualizes the
fixed exposure against those permuted conditions. This differs from keeping
the nuisance projection fixed. Family directions are invariant under the
required within-family permutations and are removed algebraically.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
os.environ.setdefault('OMP_NUM_THREADS','4')
import argparse,time
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--prefix',required=True)
    for k in ['n','k','permutations','df']:p.add_argument('--'+k,type=int,required=True)
    a=p.parse_args();prefix=Path(a.prefix)
    def read(suffix,dtype,shape):return np.fromfile(str(prefix)+suffix,dtype=dtype).reshape(shape,order='F')
    E=read('_E.bin','<f8',(a.n,a.n));Q=read('_Q.bin','<f8',(a.n,a.k));r=read('_X.bin','<f8',(a.n,))
    P=read('_P.bin','<i4',(a.permutations,a.n))-1
    assert np.isfinite(E).all() and np.isfinite(Q).all() and np.isfinite(r).all()
    total=np.trace(E);rnorm=np.sum(r*r);result=np.empty(a.permutations);start=time.monotonic()
    for lo in range(0,a.permutations,40):
        hi=min(lo+40,a.permutations)
        X=np.column_stack([r[np.argsort(v)] for v in P[lo:hi]])
        residual_norm=rnorm-np.sum((Q.T@X)**2,axis=0)
        if np.any(residual_norm<=1e-10):raise ValueError('No permuted exposure residual variation')
        numerator=np.sum(X*(E@X),axis=0)/residual_norm
        result[lo:hi]=numerator/((total-numerator)/a.df)
        if lo==0 or hi%1000==0 or hi==a.permutations:print(f'Partial dbRDA vegan 2.7-6: {hi}/{a.permutations}; {time.monotonic()-start:.1f}s',flush=True)
    assert np.isfinite(result).all()
    np.savetxt(str(prefix)+'_F.csv',result,delimiter=',',header='F_null',comments='',fmt='%.17g')
if __name__=='__main__':main()
