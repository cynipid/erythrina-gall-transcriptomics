Figure 5

Files: Figure_5.pdf (vector PDF), Figure_5.png (full-resolution PNG), plot.py, caption.txt and adjacent frozen source input(s).
Rebuild: python plot.py
Test without overwriting archived outputs: python plot.py --outdir /path/to/test_output
The script resolves frozen inputs relative to itself and runs from any working directory. No differential-expression fitting, PCA, enrichment, secretion prediction, MEROPS search, or other statistical analysis is performed.

Dependencies are documented in the script header and the package requirements.txt.
The script uses the frozen 806-gene defense overlay and regenerates the figure source-summary TSV in the chosen output directory. Ten displayed candidates include TRINITY_DN27036_c0_g1; the obsolete TRINITY_DN3635_c0_g1 substitution is absent. Lowercase panel labels match the approved caption. No candidate or numeric value changed.
