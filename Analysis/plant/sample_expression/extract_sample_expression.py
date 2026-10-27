#!/usr/bin/env python3
"""Extract and verify frozen sample expression; no model fitting."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((HERE / 'source_manifest.json').read_text())
    tables = {}
    for name, meta in manifest.items():
        compressed = (HERE / meta['path']).read_bytes()
        assert hashlib.sha256(compressed).hexdigest() == meta['sha256_compressed'], name
        raw = gzip.decompress(compressed)
        assert hashlib.sha256(raw).hexdigest() == meta['sha256_uncompressed'], name
        records = list(csv.DictReader(raw.decode().splitlines(), delimiter='\t'))
        assert len(records) == meta['rows'], name
        tables[name] = {r['gene_id']: r for r in records}
        assert len(tables[name]) == len(records), name
    source = json.loads((HERE / 'sample_expression_source.json').read_text())
    counts = tables[source['sources'][0]]
    stats = tables[source['sources'][1]]
    membership = list(csv.DictReader((HERE / 'candidate_membership.tsv').open(), delimiter='\t'))
    assert len(membership) == 32
    rows = []
    for member in membership:
        gene = member['gene_id']
        c, s = counts[gene], stats[gene]
        lfc, padj = float(s['log2FoldChange']), float(s['padj'])
        status = ('GP-up' if lfc > 0 else 'UP-up') if padj < .01 and abs(lfc) >= 1 else 'Non-DE'
        row = [gene, status] + [member[k].upper() == 'TRUE' for k in source['headers'][2:5]]
        row += [float(c[k]) for k in source['sample_mapping']]
        row += [float(s[k]) for k in ('baseMean', 'log2FoldChange', 'lfcSE', 'padj')]
        rows.append(row)
    expected = {r[0]: r for r in source['rows']}
    for row in rows:
        for actual, reference in zip(row, expected[row[0]]):
            if isinstance(actual, float):
                assert math.isclose(actual, reference, rel_tol=1e-12, abs_tol=1e-12), row[0]
            else:
                assert actual == reference, row[0]
    args.outdir.mkdir(parents=True, exist_ok=True)
    with (args.outdir / 'sample_expression.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(source['headers'])
        writer.writerows(rows)
    print('PASS: 32 clusters, 192 normalized counts and frozen statistics verified; no fitting.')

if __name__ == '__main__':
    main()
