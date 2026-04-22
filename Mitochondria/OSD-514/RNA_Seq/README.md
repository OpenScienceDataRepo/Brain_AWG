# FOLDER PATHS
```
OSD-514/
└── RNA_Seq/
		├── RawCounts/
		├── Collapsed_Counts/
				├── counts_DESeq_ready.csv
				├── metadata_from_filenames.csv
				└── OSD514_RSEM_expected_counts.csv
		├── RESULTS_OSD514/		     # General pathway, sex factored out
				├── tables/
				├── figs/
				├── GSEA/
						├── tables/
						└── figs/
    ├── QC_Sex_Cov/		             # Sex-stratified QC, sex as covariate
				├── tables/
				└── figs/
		├── DEG_Sex_Cov/		     # Sex-stratified DEG, sex as covariate
				├── tables/
				└── figs/
    └── GSEA_Sex_Cov/		          # Sex-stratified GSEA, sex as covariate
				├── tables/
				└── figs/
└── Proteomics/
		├── TMT_all/
		└── TMT_COMBINED/
```
# Extract Raw Counts

OSD-514 did not have a clean raw counts file to run in DESeq2. It had an "unnormalized counts" file but that had non-integer values, meaning it was a convenience matrix of RSEM “expected counts” (decimals).

The 24 files compiled are also not raw counts files. The *.genes.results files from RSEM contain expected counts (fractional). We used them because DESeq2 supports this workflow via tximport. We gave DESeq2 the expected counts plus the average transcript lengths that tximport extracts. DESeq2 then uses these as offsets (length-aware normalization).

Just rounding the “Unnormalized Counts” CSV throws away the information RSEM provides. It biases low counts, meaning small fractional differences can flip to 0 or 1 after rounding, distorting dispersion and p-values. It also breaks length-aware normalization since tximport+DESeq2 can properly correct for gene length and library size.

Because of this, the 24 raw count results were downloaded for each sample in each condition (saved into the path OSD-514 > RNA-Seq > RawCounts) and collapsed into one matrix manually using the Collapsed_Counts.r code in this folder.

# QC & Downstream Analysis
From here, we ran two different paths of downstream analysis:
1. **General:** Tests what changes with condition (ignoring sex)
2. **Sex-Stratified:** Tests what changes with condition (after accounting for sex differences)?

## 01. General Pathway
The script for the full pathway can be found here.

The following processes were run on the data:
1. Quality Control
2. DESeq2
3. Visualization (Volcano Plots, Heatmaps, etc)
4. Mito-Specific GO ORA
5. GSEA

We also created a script that would quantify the DEGs and DEPs by sex and condition. This was done to create a comparative figure for the manuscript, and can be found here.

## 02. Sex-Stratified Pathway
The script for the full pathway can be found here.

The following processes were run on the data:
1. Quality Control
2. DESeq2
3. Visualization (Volcano Plots, Heatmaps, etc)
4. Mito-Specific GO ORA
5. GSEA

## 03. Proteomics
The proteomics processing followed the same logic— one general, one sex-stratified— and can be found here.

A similar "README" will be stored in that folder as well.
