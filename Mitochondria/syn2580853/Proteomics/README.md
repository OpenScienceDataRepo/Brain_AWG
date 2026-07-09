# FOLDER PATHS
```
syn2580853/
└── Proteomics/
    ├── TMT_all/
    │   └── baylor_proteomics_fly_proteinoutput.txt
    ├── RESULTS_syn2580853/
    │   ├── tables/
    │   └── figs/
    ├── TauProteomics.R
    ├── EmoryDrosophilaTau_biospecimen_metadata.csv
    └── README.md
```

# Data Collection

The following files from Synapse project **syn2580853** (EmoryDrosophilaTau) are used for analysis:

1. **Baylor proteomics output** — `baylor_proteomics_fly_proteinoutput.txt` (MaxQuant LFQ intensities, 27 samples)
2. **Biospecimen metadata** — `EmoryDrosophilaTau_biospecimen_metadata.csv` (genotype, age, assay)

Place the Baylor file in a local folder called `TMT_all/`. Place the biospecimen metadata in the `Proteomics/` folder alongside `TauProteomics.R`.

# Genotypes & Timepoints

| Genotype | Construct | Samples |
|---|---|---|
| Control | Elav-GAL4/+ | Days 1, 10, 20 (n=3 each) |
| TauWT | Elav-GAL4/+; UAS-TauWT/+ | Days 1, 10, 20 (n=3 each) |
| TauR406W | Elav-GAL4/+; UAS-TauR406W/+ | Days 1, 10, 20 (n=3 each) |

# QC & Downstream Analysis

The full proteomics pathway is in [TauProteomics.R](TauProteomics.R).

The following processes are run on the data:

1. Quality Control (PCA)
2. limma differential expression (9 contrasts)
3. Visualization (volcano plots, heatmaps)
4. Mitochondrial protein subset
5. enrichR GO ORA
6. GSEA (fgsea, supplementary contrasts)

## Contrasts

| Contrast | Description |
|---|---|
| `DE_TauWT_vs_Control_Day*` | Wild-type Tau vs Control |
| `DE_TauR406W_vs_Control_Day*` | R406W mutant vs Control |
| `DE_TauR406W_vs_TauWT_Day*` | Mutant vs wild-type Tau (supplementary) |

## Key Output Tables

| File | Description |
|---|---|
| `proteomics_metadata.csv` | Sample genotype and day |
| `DE_TauWT_vs_Control_Day*.csv` | TauWT DE results |
| `DE_TauR406W_vs_Control_Day*.csv` | TauR406W DE results |
| `DE_mitochondrial_proteins_*.csv` | Mitochondrial protein subsets |
| `GO_proteomics_*.csv` | enrichR GO results |

## How to Run

1. Download data files into `TMT_all/` and `Proteomics/`
2. Update `setwd()` in `TauProteomics.R` to your local `Proteomics/` path
3. Run the script in RStudio

```r
setwd("Your_Working_Directory/Proteomics")
source("TauProteomics.R")
```

## Proteomics

The RNA-seq processing for this dataset can be found [here](../RNA_Seq).
