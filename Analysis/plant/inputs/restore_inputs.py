#!/usr/bin/env python3
"""Restore frozen input bytes into a new directory; run no analyses."""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parent


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            result.update(block)
    return result.hexdigest()


def safe_relative(value):
    path = Path(value)
    if path.is_absolute() or '..' in path.parts:
        raise ValueError('Invalid relative path in manifest')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    destination = args.destination.expanduser().resolve()
    if destination.exists():
        parser.error('Destination already exists; choose a new directory.')
    with (ROOT / 'manifest.tsv').open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for row in rows:
        source = ROOT / safe_relative(row['compressed_path'])
        safe_relative(row['project_path'])
        if digest(source) != row['sha256_compressed']:
            raise ValueError('Compressed file checksum mismatch: ' + source.name)
    destination.mkdir(parents=True)
    for row in rows:
        target = destination / safe_relative(row['project_path'])
        target.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(ROOT / row['compressed_path'], 'rb') as source:
            with target.open('xb') as output:
                shutil.copyfileobj(source, output)
        if target.stat().st_size != int(row['original_bytes']):
            raise ValueError('Restored size mismatch: ' + str(target))
        if digest(target) != row['sha256_original']:
            raise ValueError('Restored checksum mismatch: ' + str(target))
    print(f'Restored and verified {len(rows)} files in {destination}')
    print('No analysis scripts were run.')


if __name__ == '__main__':
    main()
