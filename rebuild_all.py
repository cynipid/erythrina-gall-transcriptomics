#!/usr/bin/env python3
"""Rebuild all ten main figures using only files in this package."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
for number in range(1, 11):
    folder = ROOT / 'Main_Figures' / f'Figure_{number}'
    name = 'plot.py' if number in {2, 3, 6, 7} else f'Figure_{number}_plotting_script.py'
    print(f'Rebuilding Figure {number}', flush=True)
    subprocess.run([sys.executable, str(folder/name)], cwd=ROOT, check=True)
print('All ten figures rebuilt successfully.')
