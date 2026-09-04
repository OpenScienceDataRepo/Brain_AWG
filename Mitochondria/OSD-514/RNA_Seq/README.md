# FOLDER PATHS
```
OSD-514/
└── RNA_Seq/
		├── RESULTS_OSD514/		     
				├── tables/
				└── figs/
		└── GSEA/					
			├── tables/
			└── figs/
└── Proteomics/
		├── TMT_all/
		├── RESULTS_OSD514/		     
				├── tables/
				└── figs/
		└── GSEA/					
			├── tables/
			└── figs/
```
# Extract Raw Counts

OSD-514 did not have a clean raw counts file to run in DESeq2. It had an "unnormalized counts" file but that had non-integer values, meaning it was a convenience matrix of RSEM “expected counts” (decimals).

The 24 files compiled are also not raw counts files. The *.genes.results files from RSEM contain expected counts (fractional). We used them because DESeq2 supports this workflow via tximport. We gave DESeq2 the expected counts plus the average transcript lengths that tximport extracts. DESeq2 then uses these as offsets (length-aware normalization).

Just rounding the “Unnormalized Counts” CSV throws away the information RSEM provides. It biases low counts, meaning small fractional differences can flip to 0 or 1 after rounding, distorting dispersion and p-values. It also breaks length-aware normalization since tximport+DESeq2 can properly correct for gene length and library size.

# QC & Downstream Analysis

## 01. General Pathway
The script for the full pathway can be found here.

The following processes were run on the data:
1. Quality Control
2. DESeq2
3. Visualization (Volcano Plots, Heatmaps, etc)
4. Mito-Specific GO ORA
5. GSEA

## 02. Proteomics
The proteomics processing can be found here.

A similar "README" will be stored in that folder as well.
