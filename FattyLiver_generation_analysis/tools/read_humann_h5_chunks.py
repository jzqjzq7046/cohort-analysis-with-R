"""Read-only, indexed HDF5 import in bounded feature chunks; no analysis/plots."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0, 'E:/YFS/workspace/.python_packages')
import h5py
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('h5', type=Path)
p.add_argument('manifest', type=Path)
p.add_argument('output', type=Path)
a = p.parse_args()
m = json.loads(a.manifest.read_text(encoding='utf-8-sig'))
rows = np.asarray(m['zero_based_row_indices'], dtype=int)
samples = np.asarray(m['zero_based_sample_indices'], dtype=int)
assert len(rows) and np.all(np.diff(rows) > 0)
assert len(samples) and len(set(samples.tolist())) == len(samples)
a.output.mkdir(parents=True, exist_ok=True)
chunks = []
with h5py.File(a.h5, 'r') as f:
    ds = f[m['dataset']]
    assert tuple(ds.shape) == (m['n_samples'], m['n_source_rows'])
    assert samples.min() >= 0 and samples.max() < ds.shape[0]
    for start in range(0, len(rows), 768):
        indices = rows[start:start + 768]
        data = np.asarray(ds[:, indices], dtype='<f8')[samples, :]
        assert data.shape == (len(samples), len(indices))
        assert np.isfinite(data).all() and (data >= 0).all()
        name = f'chunk_{start // 768 + 1:03d}.bin'
        data.tofile(a.output / name)
        chunks.append(dict(file=name, start_one_based=start + 1,
                           rows=len(indices), samples=len(samples)))
report = dict(dataset=m['dataset'], original_HDF5_read_only=True,
              selected_samples=len(samples), selected_rows=len(rows),
              h5py_version=h5py.__version__, chunks=chunks)
(a.output / 'chunks.json').write_text(json.dumps(report, indent=2))
print(json.dumps({k: v for k, v in report.items() if k != 'chunks'}))
