"""Read an explicitly indexed existing HUMAnN HDF5 dataset; never draws or fits.

R writes the row-index manifest and metadata, and reads the feature x sample
binary result. Data remain local and original HDF5/RDS files remain read-only.
"""
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
indices = np.asarray(m['zero_based_row_indices'], dtype=int)
assert len(indices) > 0 and np.all(np.diff(indices) > 0)
with h5py.File(a.h5, 'r') as f:
    ds = f[m['dataset']]
    assert tuple(ds.shape) == (m['n_samples'], m['n_source_rows'])
    data = np.asarray(ds[:, indices], dtype='<f8')
assert data.shape == (m['n_samples'], len(indices))
assert np.isfinite(data).all() and (data >= 0).all()
a.output.parent.mkdir(parents=True, exist_ok=True)
data.tofile(a.output)
report = {'dataset': m['dataset'], 'h5_shape': [m['n_samples'], m['n_source_rows']],
          'binary_R_dimensions': [len(indices), m['n_samples']], 'dtype': 'float64 little-endian',
          'minimum': float(data.min()), 'maximum': float(data.max()),
          'sample_sum_quantiles': np.quantile(data.sum(axis=1), [0, .25, .5, .75, 1]).tolist(),
          'original_HDF5_read_only': True, 'h5py_version': h5py.__version__}
a.output.with_suffix('.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report))
