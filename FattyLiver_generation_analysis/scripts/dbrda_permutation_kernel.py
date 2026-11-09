"""Numerically equivalent reduced-model dbRDA permutation statistics.

Called by standalone R entry 08; NumPy OpenBLAS accelerates matrix products.
All files stay in the run's private folder. No network or raw data changes.
The R caller verifies every statistic in the first 10 draws against vegan.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import argparse
import time
from pathlib import Path
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--permutations", type=int, required=True)
    parser.add_argument("--df", type=int, required=True)
    args = parser.parse_args()
    prefix = Path(args.prefix)
    n, k, B = args.n, args.k, args.permutations
    E = np.fromfile(str(prefix) + "_E.bin", dtype="<f8").reshape((n, n), order="F")
    Q = np.fromfile(str(prefix) + "_Q.bin", dtype="<f8").reshape((n, k), order="F")
    P = np.fromfile(str(prefix) + "_P.bin", dtype="<i4").reshape((B, n), order="F") - 1
    if not (np.isfinite(E).all() and np.isfinite(Q).all() and (P >= 0).all() and (P < n).all()):
        raise ValueError("Invalid kernel inputs")
    total = np.trace(E)
    output = np.empty(B)
    start_time = time.monotonic()
    # Ep = E[p,p]. Quadratic forms use inverse-permuted Q instead of moving E.
    # Q[:,0] is the partial steatosis direction. Other columns span nuisance
    # covariates after family centering. Family directions contribute zero,
    # since E has zero within-family row/column sums and p stays in families.
    for start in range(0, B, 20):
        stop = min(start + 20, B)
        block = np.concatenate([Q[np.argsort(p), :] for p in P[start:stop]], axis=1)
        quadratic = np.sum(block * (E @ block), axis=0).reshape((stop - start, k))
        numerator = quadratic[:, 0]
        residual = total - quadratic[:, 1:].sum(axis=1) - numerator
        if (residual <= 0).any():
            raise ValueError("Non-positive residual inertia")
        output[start:stop] = numerator / (residual / args.df)
        if start == 0 or stop % 1000 == 0 or stop == B:
            print(f"08 optimized dbRDA: {stop}/{B}; {time.monotonic()-start_time:.1f}s", flush=True)
    if not np.isfinite(output).all():
        raise ValueError("Non-finite F statistics")
    np.savetxt(str(prefix) + "_F.csv", output, delimiter=",", header="F_null", comments="", fmt="%.17g")
    print(f"NumPy {np.__version__}; all {B} statistics complete", flush=True)

if __name__ == "__main__":
    main()
