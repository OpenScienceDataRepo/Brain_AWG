# PLEASE NOTE THIS FILE IS INCOMPLETE, UPDATES COMING SOON
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
		├── RESULTS_Sex_Strat/		     # Sex-stratified DEG, sex as interaction term
				├── tables/
				└── figs/
    └── GSEA_Sex_Strat/		          # Sex-stratified GSEA, sex as interaction term
				├── tables/
				└── figs/
└── Proteomics/
		├── TMT_all/
		└── TMT_COMBINED/
```
# Data Collection
The following processed proteomic files from OSD-514 were used for analysis:

1. GLDS-514_proteomics_TMTc.tar.gz
2. GLDS-514_proteomics_TMTb.tar.gz
3. GLDS-514_proteomics_TMTa.tar.gz
4. OSD-514_metadata_OSD-514-ISA.zip (contained metadata)

The TMT files were manually added to a folder called "TMT_all" in order to be combined into one matrix using [this script.](Mitochondria/OSD-514/Proteomics/CombineTMT.r)
