# Frozen plant sample expression

The membership table selects 32 unique clusters: 18 named in the main text, 10 shown in Figure 5c and 19 in Figure 6c (sets overlap). The source JSON retains 192 normalized counts, unshrunk log2 fold changes, lfcSE, adjusted P values and membership flags. This is extraction from the retained final fit, with no model fitting or new candidate selection.

Run from the repository root:

```bash
python Analysis/plant/sample_expression/extract_sample_expression.py --outdir /path/to/new/expression
```

The standard-library helper verifies both compressed and uncompressed input hashes, row counts and all selected values against the frozen source JSON before exporting TSV. `source_manifest.json` identifies inputs; `sample_expression_validation.json` preserves the earlier validation. Sample IDs S7–S9 map to GP1–GP3 and S10–S12 to UP1–UP3. The 806-candidate inventory is retained in S13 and in the Figure 5 overlay. Normalized counts describe sample variation; these six petioles originated from one tree and one sampling occasion.
