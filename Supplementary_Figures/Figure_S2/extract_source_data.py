#!/usr/bin/env python3
"""Read the stored VST assay and export its Pearson sample-correlation matrix.

This does not run DESeq2, VST, PCA, fitting, normalization, or gene filtering.
Dependencies: rdata, numpy, pandas, xarray. Run python extract_source_data.py.
The matrix-only plot.py needs only numpy, pandas and matplotlib.
"""
from pathlib import Path
import gzip, json, hashlib, warnings
import numpy as np
import pandas as pd
import rdata

HERE = Path(__file__).resolve().parent
path = HERE / 'insect_gene_vsd.rds.gz'
expected_hash = '758e139a6eb7eaeecec576c4d3c8c97eb73834826c39e1a3fc3493e1bc8c3e10'
raw = path.read_bytes()
assert hashlib.sha256(raw).hexdigest() == expected_hash
unwrapped = raw
gzip_layers = 0
while unwrapped[:2] == b'\x1f\x8b':
    unwrapped = gzip.decompress(unwrapped)
    gzip_layers += 1
with warnings.catch_warnings(record=True) as captured:
    warnings.simplefilter('always')
    obj = rdata.conversion.convert(rdata.parser.parse_data(unwrapped))
assays = obj.assays.data.listData
assert len(assays) == 1
assay = assays[0]
assert assay.shape == (13919, 6)
genes = list(map(str, assay.coords['dim_0'].values))
original_ids = list(map(str, assay.coords['dim_1'].values))
assert len(set(genes)) == 13919
assert original_ids == ['L1', 'L3', 'L4', 'P3', 'P4', 'P5']
assert list(map(str, obj.colData.rownames)) == original_ids
assert list(map(str, obj.colData.listData['stage'])) == ['larva'] * 3 + ['pupa'] * 3
key = pd.read_csv(HERE / 'sample_key.tsv', sep='\t')
assert key['sample_id'].tolist() == original_ids
display_ids = key['sample_label'].tolist()
assert display_ids == ['L1', 'L2', 'L3', 'P1', 'P2', 'P3']
values = np.asarray(assay.values, dtype=float)
assert np.isfinite(values).all()
matrix = np.corrcoef(values, rowvar=False)
accepted_rounded = np.array([
    [1,.98,.98,.84,.83,.84], [.98,1,.98,.82,.81,.82],
    [.98,.98,1,.82,.82,.82], [.84,.82,.82,1,.97,.98],
    [.83,.81,.82,.97,1,.99], [.84,.82,.82,.98,.99,1],
])
assert np.allclose(np.round(matrix, 2), accepted_rounded, atol=0, rtol=0)
pd.DataFrame(matrix, index=display_ids, columns=display_ids).to_csv(
    HERE / 'Supplementary_Figure_S2_source_data.tsv', sep='\t', float_format='%.17g')
audit = {
    'frozen_rds_sha256': expected_hash, 'decoder': 'rdata ' + rdata.__version__,
    'numpy': np.__version__, 'assay_dimensions': [13919, 6],
    'original_sample_ids': original_ids, 'display_sample_ids': display_ids,
    'sample_order_checked_against': 'RDS assay dimnames, RDS colData, frozen PCA sample key',
    'gzip_layers_unwrapped': gzip_layers,
    'calculation': 'Pearson sample correlations of stored frozen VST assay values',
    'every_rounded_coefficient_matches_retained_artwork': True,
    'model_or_vst_or_pca_recomputed': False,
    'decoder_warnings': [str(w.message) for w in captured],
}
(HERE / 'source_extraction_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps(audit, indent=2))
