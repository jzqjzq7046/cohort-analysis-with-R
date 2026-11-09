"""Freedman--Lane MRQAP studentized statistics for a rectangular dyadic matrix.

The R caller supplies reduced-model residuals, orthonormal nuisance directions,
the residualized unit-length exposure, and role-specific node permutations.
Every node has one permutation identity across all its dyads.  These are not
independent edge shuffles.  Inputs and outputs remain in a private run folder.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import argparse
import time
from pathlib import Path
import numpy as np


def studentized_statistics(E, Q, x_unit, rows, cols, df, chunk_size=20):
    na, nb = E.shape
    if Q.shape[0] != na * nb or x_unit.size != na * nb:
        raise ValueError("Residual and design dimensions disagree")
    if rows.shape[0] != cols.shape[0] or df <= 0:
        raise ValueError("Invalid permutation dimensions or residual degrees of freedom")
    if not all(np.isfinite(z).all() for z in (E, Q, x_unit)):
        raise ValueError("Non-finite model inputs")
    if not (np.all(np.sort(rows, axis=1) == np.arange(na)) and
            np.all(np.sort(cols, axis=1) == np.arange(nb))):
        raise ValueError("Every draw must permute complete nodes within each role")
    B = rows.shape[0]
    statistics = np.empty(B)
    norm2 = np.sum(E * E)
    started = time.monotonic()
    for start in range(0, B, chunk_size):
        stop = min(start + chunk_size, B)
        block = np.empty((na * nb, stop - start), order="F")
        for j, draw in enumerate(range(start, stop)):
            block[:, j] = E[np.ix_(rows[draw], cols[draw])].ravel(order="F")
        numerator = x_unit @ block
        projections = Q.T @ block
        rss = norm2 - np.sum(projections * projections, axis=0) - numerator * numerator
        if np.any(rss <= 0):
            raise ValueError("Non-positive residual sum of squares")
        statistics[start:stop] = numerator / np.sqrt(rss / df)
        if start == 0 or stop % 1000 == 0 or stop == B:
            print(f"07 MRQAP: {stop}/{B}; {time.monotonic()-started:.1f}s", flush=True)
    if not np.isfinite(statistics).all():
        raise ValueError("Non-finite permutation t statistics")
    return statistics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--na", type=int, required=True)
    parser.add_argument("--nb", type=int, required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--permutations", type=int, required=True)
    parser.add_argument("--df", type=int, required=True)
    parser.add_argument("--chunk", type=int, default=20)
    args = parser.parse_args()
    prefix = Path(args.prefix)
    na, nb, k, B = args.na, args.nb, args.k, args.permutations
    def read(suffix, dtype, shape):
        return np.fromfile(str(prefix) + suffix, dtype=dtype).reshape(shape, order="F")
    E = read("_E.bin", "<f8", (na, nb))
    Q = read("_Q.bin", "<f8", (na * nb, k))
    x_unit = read("_X.bin", "<f8", (na * nb,))
    rows = read("_ROW.bin", "<i4", (B, na)) - 1
    cols = read("_COL.bin", "<i4", (B, nb)) - 1
    output = studentized_statistics(E, Q, x_unit, rows, cols, args.df, args.chunk)
    np.savetxt(str(prefix) + "_t.csv", output, delimiter=",", header="t_null",
               comments="", fmt="%.17g")
    print(f"NumPy {np.__version__}; all {B} statistics complete", flush=True)


if __name__ == "__main__":
    main()
