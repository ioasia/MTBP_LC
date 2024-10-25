## MTBP_LC

# MTBP signature analysis

1) Cell line marker extraction ('High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer_markers_mmc3.xlsx'). Other lists available as well.
2) timsTOF data(log10 transform, median/mean sweep, plot distributions) and log2 mRNA for overlapping samples)
3) Correlation calculation (share functions for symbol harmonization/correlation)
4)	Threshold definition for significant correlations based on a specific number of samples (use known LC markers for biologically driven decision, alternatively distribution-driven)
5)	Filter to cell type lists with at least one protein (marker-marker correlation)
6)	Gene-list based GSVA (ssgsea) for cell type enrichment score (distribution of scores across cases)
(maybe Xcell, BayesPrism , Bisque, CIBERSORTx, DeconRNAseq, dtangle, DWLS, EPIC, MCP- counter, Scaden, and SCDC, BayesDeBulk, ESTIMATE). 
7)	Similar analysis for pathway estimation ([PROGENY](https://www.nature.com/articles/s41467-017-02391-6))
