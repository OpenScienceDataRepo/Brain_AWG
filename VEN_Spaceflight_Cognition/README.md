# VEN Fatigue Hypothesis: Social Cognition Paradox in Long-Duration Spaceflight

**Paper:** "The Social Cognition Paradox in Long-Duration Spaceflight: A VEN Fatigue Hypothesis for Duration-Dependent Emotion Recognition Decline"  
**Authors:** Esila Keskin (UWE Bristol), Nilufar Ali (Boise State / NASA GeneLab Brain AWG), Margaret Windy McNerney (Stanford)  

---

## Overview

This repository contains all analysis code and data (where distributable) for the paper.

The study derives and tests the **VEN Fatigue Hypothesis**: Von Economo neurons (VENs), optimised for fast social decisions on minimal evidence, are overactivated by gravity-disrupted social cues during spaceflight. This drives a compensatory myelination upregulation that sustains Emotion Recognition Task (ERT) speed within 6-month missions. Beyond this threshold, fatigue accumulates faster than compensation, producing selective late-inflight ERT speed decline.

Seven convergent lines of evidence are tested across five independent datasets spanning cognitive, transcriptomic, proteomic, and human blood-level analyses.

---

## Seven lines of evidence

1. **ERT cognitive paradox** -- ERT speed is stable across 6-month ISS missions (N=24; +0.106 SD) but shows the largest early-to-late inflight decline of any cognitive speed domain in the 340-day NASA Twins Study (N=1; -1.8 SD; within-individual cross-domain t=11.41, df=8)
2. **ISS rodent myelination** -- Myelination genes are specifically upregulated in ISS frontal cortex (GSE239336; +0.381 log2FC; 3.89 SD above genome-wide null; permutation p=0.0001)
3. **Ground analogue control** -- No specific myelination signal in OSD-202 hindlimb unloading + radiation (permutation p=0.164), with directional reversal consistent with a qualitatively distinct molecular response
4. **Human cortical organoids** -- Layer V projection genes massively upregulated in ISS cortical organoids (GSE259421; +1.347 log2FC; 6.17 SD above null; permutation p=0.0008)
5. **Dopaminergic organoid replication** -- Layer V signal replicates in dopaminergic (midbrain) organoids from the same mission (7.12 SD; p=0.0002), confirming the response is not cortical-lineage-specific
6. **Proteomics (PXD069807)** -- VEN-associated structural proteins in iPSC organoid proteomics (Jourdon et al. 2026) follow global ISS trends; Layer V transcription factors fall below mass-spectrometry detection threshold as predicted
7. **SOMA atlas human blood** -- All eight queried VEN panel genes reach p<0.05 in human astronaut blood (NASA Twins Study + Inspiration4); MOG log2FC +5.82 (p=4.59e-27) in year-long ISS, replicating in Inspiration4 PBMC (MBP log2FC +0.60)

---

## Repository structure

```
ven-spaceflight-cognition/

  figures/ -- Publication figures (PNG)
    fig1_ert_paradox.png
    fig2_molecular_dissociation.png
    fig3_organoid_permutation.png
    fig4_domain_specificity.png
    fig5_soma_blood.png
    fig5_proteomics_permutation.png
    dopaminergic_vs_cortical_comparison.png

  data/raw/
    twins_cognitive_heatmap.csv -- Digitised from Garrett-Bakelman et al. 2019, Fig 10B
    dev2024_raw_scores.csv -- From Dev et al. 2024, Table 2
    GSE239336_FCT_GCvsFLT-SAL_DEanalysis.txt  -- ISS rodent DE analysis (NCBI GEO)
    GLDS-202_...csv -- OSD-202 (NASA OSDR) [gitignored, 134 MB]
    GSE259421_all_counts.txt -- Human organoid counts (NCBI GEO) [gitignored, 45 MB]

  results/ -- All analysis outputs (CSV, JSON)

  step1_create_cognitive_csvs.py -- Creates cognitive data CSVs from digitised values
  step2_run_analysis.py -- ISS rodent (GSE239336) and OSD-202 analyses
  step2_v2_run_analysis.py -- Updated molecular analysis pipeline
  step3_make_figures.py -- Generates publication figures
  make_paper_figures.py -- Updated figure generation script
  ven_organoid_analysis.py -- Core cortical organoid VEN panel pipeline
  ven_dopaminergic_analysis.py -- Dopaminergic organoid replication analysis
  soma_analysis.py -- SOMA atlas VEN panel survey (human astronaut blood)
  alysson_proteomics_analysis.py -- PatternLab V export parser for PXD069807 proteomics
  parse_sepr2.py -- MS-NRBF binary parser for .sepr2 PatternLab files
  inspect_sepr2.py -- Diagnostic utility for .sepr2 binary structure
```

---

## Data availability

| Dataset | Source | Size | Status |
|---------|--------|------|--------|
| `twins_cognitive_heatmap.csv` | Digitised from Garrett-Bakelman et al. 2019, Science | 1 KB | Included |
| `dev2024_raw_scores.csv` | Dev et al. 2024, Table 2 | 1 KB | Included |
| `GSE239336_FCT_GCvsFLT-SAL_DEanalysis.txt` | NCBI GEO: GSE239336 | ~1 MB | Included |
| `GLDS-202_rna_seq_...csv` | NASA OSDR: OSD-202 | 134 MB | Gitignored -- download from [NASA OSDR](https://osdr.nasa.gov) |
| `GSE259421_all_counts.txt` | NCBI GEO: GSE259421 | 45 MB | Gitignored -- download from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE259421) |
| `ISS_WT83_proteins.txt` / `Ground_WT83_proteins.txt` | PRIDE: PXD069807 (Jourdon et al. 2026) | varies | Gitignored -- download from [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD069807) |
| SOMA atlas queries | soma.weill.cornell.edu (Overbey et al. Nature 2024) | N/A | Values extracted manually; hardcoded in `soma_analysis.py` |

---

## Reproducing the analysis

### 1. Install dependencies

```bash
pip install numpy pandas scipy matplotlib mygene
```

### 2. Download large data files (gitignored)

- OSD-202: download `GLDS-202_rna_seq_differential_expression_GLbulkRNAseq.csv` from [NASA OSDR](https://osdr.nasa.gov), place in `data/raw/`
- GSE259421: download `GSE259421_all_counts.txt` from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE259421), place in the root directory
- PXD069807: download PatternLab V exports (`ISS_WT83_proteins.txt`, `Ground_WT83_proteins.txt`) from [PRIDE PXD069807](https://www.ebi.ac.uk/pride/archive/projects/PXD069807), place in the root directory

### 3. Run the analysis pipeline

```bash
# Step 1: create cognitive performance CSVs from digitised values
python step1_create_cognitive_csvs.py

# Step 2a: ISS rodent transcriptomics (GSE239336) and OSD-202 ground analogue
python step2_run_analysis.py

# Step 2b: human cortical organoid VEN panel analysis (GSE259421)
python ven_organoid_analysis.py

# Step 2c: dopaminergic organoid replication
python ven_dopaminergic_analysis.py

# Step 2d: proteomics analysis (PXD069807)
python alysson_proteomics_analysis.py

# Step 2e: SOMA atlas human blood survey
python soma_analysis.py

# Step 3: generate publication figures
python make_paper_figures.py
```

Results are saved to `results/` and figures to `figures/`.

---

## VEN gene panel

The 31-gene panel was defined in [Keskin 2026 (arXiv:2604.09229)](https://arxiv.org/abs/2604.09229) prior to accessing any dataset analysed here. All directional predictions were pre-registered via the NASA OSDR Brain AWG forum (April 2026).

| Category | Genes |
|----------|-------|
| Myelination (n=7) | MBP, MOG, PLP1, MAG, CNP, MOBP, ERMN |
| Fast Signalling (n=7) | SCN1A, KCNQ2, ANK3, NEFH, NEFM, NEFL, SNCG |
| Social Circuit (n=6) | OXTR, AVPR1A, HTR2A, DRD1, CHRM1, GABRB2 |
| Layer V Projection (n=5) | FEZF2, BCL11B, TBR1, SATB2, CUX1 |
| Metabolic Support (n=6) | VDAC1, ATP2B2, SLC17A7, SNAP25, SYP, NRXN1 |

NOS1 (direct VEN biochemical marker) tracked individually outside the panel.

Statistical approach: one-sample t-tests against zero and permutation tests (N=10,000 random gene sets of equal size drawn from the full expressed transcriptome; seed 42). Results reported as standard deviations above the genome-wide permutation null alongside two-tailed permutation p-values.

---

## Citation

If you use this code or data, please cite:

> Keskin E, Ali N, McNerney MW. The Social Cognition Paradox in Long-Duration Spaceflight: A VEN Fatigue Hypothesis for Duration-Dependent Emotion Recognition Decline. 2026. Preprint.

Related preprint (Fast Lane Hypothesis SNN model):

> Keskin E. Fast Lane Hypothesis: Von Economo Neurons as a Biological Speed-Accuracy Tradeoff Mechanism. arXiv:2604.09229. 2026.

---

## License

Code: MIT License. See LICENSE file.  
Data: subject to original dataset terms (NCBI GEO / NASA OSDR open data policies; PRIDE proteomics repository terms).
