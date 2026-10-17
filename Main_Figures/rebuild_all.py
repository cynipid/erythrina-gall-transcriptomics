#!/usr/bin/env python3
"""Rebuild this package's figures from its frozen inputs."""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--outdir', type=Path, help='Optional new output directory; inputs remain in this package')
args = parser.parse_args()
folders = sorted((p for p in ROOT.glob('Figure_*') if p.is_dir()), key=lambda p:int(p.name.split('_')[-1].lstrip('S')))
for folder in folders:
    target = (args.outdir.resolve()/folder.name) if args.outdir else folder
    target.mkdir(parents=True, exist_ok=True)
    print('Rebuilding '+folder.name, flush=True)
    subprocess.run([sys.executable, str(folder/'plot.py'), '--outdir', str(target)], check=True)
print(str(len(folders))+' figures rebuilt successfully.')
