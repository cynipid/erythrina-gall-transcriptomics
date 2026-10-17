Supplementary Figure S2 - insect sample correlations

Production PDF/PNG were cosmetically refreshed with the included portable plotting script. The complete correlation matrix
was recovered from the archived insect VST assay using rdata 1.1.0, with named assay
dimensions and colData/sample-key order checked. Every rounded coefficient matches
the accepted artwork. No VST, model, PCA, normalization or gene filter was rerun.
Frozen RDS SHA256: 758e139a6eb7eaeecec576c4d3c8c97eb73834826c39e1a3fc3493e1bc8c3e10

Run python plot.py for a portable rendering of the supplied full-precision matrix.
The numeric values and labels match the original R-generated snapshot; fonts, color
interpolation and legend geometry were refreshed for portable script/output consistency.
Optional: python extract_source_data.py regenerates the source matrix from the
included frozen RDS and sample_key.tsv. It additionally requires rdata, xarray,
numpy and pandas. Extraction warnings about unregistered Bioconductor class
constructors are recorded; the named numeric assay is decoded and validated.
