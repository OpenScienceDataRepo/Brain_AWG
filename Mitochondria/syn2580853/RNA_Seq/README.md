# FOLDER PATHS
```
syn2580853/
└── RNA_Seq/
    ├── Collapsed_Counts/
    │   ├── Trim_*_counts.txt
    │   └── metadata_from_filenames.csv
    ├── RESULTS_syn2580853/
    │   ├── tables/
    │   └── figs/
    ├── GSEA/
    │   ├── tables/
    │   └── figs/
    ├── TauAnalysis.R
    └── README.md
```

# Data Collection

The following files from Synapse project **syn2580853** (EmoryDrosophilaTau) are used for analysis:

1. **18 RNA-seq count files** — `Trim_*_counts.txt` (e.g. `Trim_DE01D01_counts.txt`)
2. **RNA-seq metadata** — `EmoryDrosophilaTau_RNAseq_metadata.csv` (optional reference)

Download the 18 `Trim_*_counts.txt` files into a local folder called `Collapsed_Counts/` before running the pipeline.

# Sample Naming

| Code | Meaning |
|---|---|
| `DE` | Control (Elav-GAL4 driver only) |
| `DT` | TauWT |
| `01` / `10` / `20` | Day 1, 10, or 20 |
| `D01`–`D03` | Biological replicate |

Example: `Trim_DT20D02_counts.txt` = TauWT, Day 20, replicate 2.

`TauAnalysis.R` parses sample names from filenames and writes `metadata_from_filenames.csv` into `Collapsed_Counts/`.

# QC & Downstream Analysis

The full RNA-seq pathway is in [TauAnalysis.R](TauAnalysis.R).

The following processes are run on the data:

1. Quality Control (library sizes, PCA)
2. DESeq2 (`~ day + genotype`, Tau vs Control)
3. Visualization (volcano plots, heatmaps)
4. Timepoint-specific contrasts (Day01, Day10, Day20)
5. Mitochondrial gene subset
6. enrichR GO ORA
7. GSEA (fgsea, GO Biological Process)

## Key Output Tables

| File | Description |
|---|---|
| `DE_Tau_vs_Control.csv` | Overall Tau vs Control |
| `DE_Tau_vs_Control_shrunk.csv` | apeglm-shrunk LFC |
| `DE_Tau_vs_Control_Day*.csv` | Per-timepoint contrasts |
| `DE_mitochondrial_genes.csv` | Curated mitochondrial gene subset |
| `sample_metadata.csv` | Sample annotations |
| `fgsea_GO_BP_Tau_vs_Control.csv` | GSEA pathway results |

## How to Run

1. Download count files into `Collapsed_Counts/`
2. Update `setwd()` in `TauAnalysis.R` to your local `RNA_Seq/` path
3. Run the script in RStudio

```r
setwd("Your_Working_Directory/RNA_Seq")
source("TauAnalysis.R")
```
