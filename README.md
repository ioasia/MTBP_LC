## MTBP_LC

# MTBP signature analysis

•	Cell line marker extraction ('High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer_markers_mmc3.xlsx'). Other lists available as well.
•	 timsTOF data(log10 transform, median/mean sweep, plot distributions) and log2 mRNA for overlapping samples)
•	Correlation calculation (share functions for symbol harmonization/correlation)
•	Threshold definition for significant correlations based on a specific number of samples (use known LC markers for biologically driven decision, alternatively distribution-driven)
•	Filter to cell type lists with at least one protein (marker-marker correlation)
•	Gene-list based GSVA (ssgsea) for cell type enrichment score (distribution of scores across cases)
(maybe Xcell, BayesPrism , Bisque, CIBERSORTx, DeconRNAseq, dtangle, DWLS, EPIC, MCP- counter, Scaden, and SCDC, BayesDeBulk, ESTIMATE). 
•	Similar analysis for pathway estimation (PROGENY)
