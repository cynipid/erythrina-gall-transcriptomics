#!/usr/bin/env python3
"""Rebuild the final six main and eight supplementary figures from frozen inputs."""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--group', choices=('all', 'main', 'supplementary'), default='all')
parser.add_argument('--outdir', type=Path, help='Optional new output directory')
args = parser.parse_args()
count = 0
for group, parent, prefix, total in [('main', 'Main_Figures', 'Figure_', 6), ('supplementary', 'Supplementary_Figures', 'Figure_S', 8)]:
    if args.group not in ('all', group):
        continue
    for number in range(1, total+1):
        folder = ROOT / parent / (prefix+str(number))
        target = args.outdir.resolve()/parent/folder.name if args.outdir else folder
        target.mkdir(parents=True, exist_ok=True)
        print('Rebuilding '+parent+'/'+folder.name, flush=True)
        subprocess.run([sys.executable, str(folder/'plot.py'), '--outdir', str(target)], check=True)
        count += 1
print(str(count)+' figures rebuilt successfully.')
