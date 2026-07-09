# OSD514 TMT proteomics processing workflow

This repository documents and preserves the processing workflow used to generate
the files in `RESULTS_OSD514_ms3fix` from the NASA/OSD514 Drosophila TMT
proteomics raw data.

The workflow starts from Thermo RAW files named like `NASA_Flies_TMT*.raw`,
converts them to `mzML`, processes the multiplexed TMT data with FragPipe, exports
Proteome Discoverer-like per-sample protein abundance tables, and then performs
limma differential expression and GO enrichment analysis in R.

## Repository layout

```text
.
|-- a_OSD-514_protein-expression-profiling_mass-spectrometry_Orbitrap Fusion.txt
|-- fragpipe_dmel_tmt/
|   |-- run_fragpipe_dmel_tmt.sh
|   |-- build_tmt_annotation.py
|   |-- extract_tmt_protein_matrix.py
|   `-- export_tmt_all_like_from_psm.py
|-- TMT_all_from_psm_pdlike_ms3fix/
|-- TMT_expression_matrix_ms3fix.csv
|-- TMT_expression_matrix_ms3fix_batch_corrected.csv
|-- sample_metadata_ms3fix.csv
|-- GEN_PROTEOMICS_SCRIPT_ms3fix.r
`-- RESULTS_OSD514_ms3fix/
```

## Inputs

The workflow uses these main inputs:

- `NASA_Flies_TMT*.raw`: Thermo Orbitrap Fusion raw files for the TMT plexes.
  These raw files are not stored in this repository.
- `mzML/*.mzML`: RAW files converted to mzML, for example with
  ThermoRawFileParser. The FragPipe runner expects the file names to contain
  `TMTA`, `TMTB`, or `TMTC` so files can be assigned to the correct TMT plex.
- `a_OSD-514_protein-expression-profiling_mass-spectrometry_Orbitrap Fusion.txt`:
  assay metadata mapping each sample to its TMT run and reporter channel.
- `dmel_UP000000803_uniprot_target_decoy.fasta`: D. melanogaster target-decoy
  FASTA used by FragPipe.

## Processing overview

```text
NASA_Flies_TMT*.raw
  -> mzML conversion with ThermoRawFileParser
  -> FragPipe TMT10-MS3 workflow
  -> FragPipe per-plex protein.tsv and psm.tsv
  -> PD-like per-sample TMT_all files
  -> protein x sample expression matrix
  -> limma differential expression
  -> PCA, heatmap, volcano plots, fgsea, FlyEnrichr outputs
  -> RESULTS_OSD514_ms3fix/
```

## 1. Convert RAW files to mzML

The raw-to-mzML conversion is required before running FragPipe, but this
conversion step is not currently scripted in this repository. Thermo RAW files
can be converted with ThermoRawFileParser, which supports Thermo `.raw` input
and mzML output on Linux, macOS, and Windows.

Directory conversion:

```bash
mkdir -p mzML

ThermoRawFileParser \
  -d=/path/to/NASA_Flies_raw_files \
  -o=/path/to/mzML \
  -f=1
```

The downstream FragPipe script assumes converted files are in an `mzML`
directory and that each filename identifies the TMT plex, for example `TMTA`,
`TMTB`, or `TMTC`.

## 2. Run FragPipe TMT processing

FragPipe processing is configured in:

```text
fragpipe_dmel_tmt/run_fragpipe_dmel_tmt.sh
```

Important implementation details:

- The runner uses the built-in FragPipe `TMT10-MS3.workflow` template.
- The workflow is patched to use reporter-ion mode:
  - `tmtintegrator.ms1_int=false`
  - `tmtintegrator.log2transformed=false`
  - `tmtintegrator.ref_tag=ctrl_pool`
  - relaxed TMT-Integrator peptide summarization gates
- `build_tmt_annotation.py` reads the assay metadata and writes one annotation
  file per plex:
  - `TMTA_annotation.txt`
  - `TMTB_annotation.txt`
  - `TMTC_annotation.txt`
- A FragPipe manifest is built by assigning every `.mzML` file to `TMTA`,
  `TMTB`, or `TMTC` from the filename.

Before rerunning elsewhere, update the path variables at the top of
`fragpipe_dmel_tmt/run_fragpipe_dmel_tmt.sh`, especially `BASE_DIR`,
`MZML_DIR`, `ASSAY_LABEL_FILE`, `TOOLS_DIR`, `SIF`, and `FALLBACK_DB`.

Run:

```bash
bash fragpipe_dmel_tmt/run_fragpipe_dmel_tmt.sh
```

Expected FragPipe-stage outputs:

```text
fragpipe_dmel_tmt/work/results/
fragpipe_dmel_tmt/logs/fragpipe_tmt_run.log
fragpipe_dmel_tmt/work/protein_abundance_tmt_matrix.tsv
TMT_all_from_psm_pdlike_ms3fix/
```

## 3. Export PD-like per-sample TMT files

The final step of `run_fragpipe_dmel_tmt.sh` calls:

```text
fragpipe_dmel_tmt/export_tmt_all_like_from_psm.py
```

This script converts FragPipe outputs into per-sample tables that mimic the
published `TMT_all`/Proteome Discoverer layout. It reads each plex's
`psm.tsv` and `protein.tsv`, applies PD-like PSM filters, aggregates reporter
ion intensities to the protein level, and writes one text file per biological
sample plus one pool file per plex.

Default PSM filters:

- `Qvalue <= 0.01`
- `Probability >= 0.90`
- `Purity >= 0.50` when a `Purity` column is present

The resulting directory in this repository is:

```text
TMT_all_from_psm_pdlike_ms3fix/
```

It contains 18 biological sample files and 3 pool files:

```text
Earth_F1_TMTa.txt      Earth_M1_TMTa.txt
Earth_F2_TMTb.txt      Earth_M2_TMTb.txt
Earth_F3_TMTc.txt      Earth_M3_TMTc.txt
SF1g_F1_TMTa.txt       SF1g_M1_TMTa.txt
SF1g_F2_TMTb.txt       SF1g_M2_TMTb.txt
SF1g_F3_TMTc.txt       SF1g_M3_TMTc.txt
SFug_F1_TMTa.txt       SFug_M1_TMTa.txt
SFug_F2_TMTb.txt       SFug_M2_TMTb.txt
SFug_F3_TMTc.txt       SFug_M3_TMTc.txt
pool_TMTa.txt          pool_TMTb.txt          pool_TMTc.txt
```

Each sample table includes protein accession, description, coverage, peptide
counts, PSM counts, unique peptide counts, reporter abundance, and reporter
abundance count.

## 4. Build the expression matrix

The downstream R workflow is:

```text
GEN_PROTEOMICS_SCRIPT_ms3fix.r
```

This script reads the assay metadata and all files in
`TMT_all_from_psm_pdlike_ms3fix/`, then builds a protein x sample matrix using
these processing rules:

1. Remove pool files from the analysis matrix.
2. Keep proteins with at least 2 unique peptides in each source table.
3. Keep finite, positive abundance values.
4. Collapse duplicate protein/sample/run rows using the median abundance.
5. Keep proteins observed in all three TMT runs.
6. Apply within-run median scaling by sample.
7. Keep only complete-case proteins across the 18 biological samples.
8. Transform abundance values with `log2(value + 1)`.

The resulting matrix has 1,668 proteins and 18 biological samples.

Outputs:

```text
RESULTS_OSD514_ms3fix/tables/TMT_expression_matrix_ms3fix.csv
RESULTS_OSD514_ms3fix/tables/sample_metadata_ms3fix.csv
TMT_expression_matrix_ms3fix.csv
sample_metadata_ms3fix.csv
```

The repository-level copies are convenience copies of the same processed matrix
and metadata.

## 5. PCA and batch-corrected matrix

The R script generates PCA plots from the raw log2 matrix and a batch-corrected
matrix. Batch correction uses:

```r
limma::removeBatchEffect(
  expr_limma,
  batch = meta_limma$tmt_run,
  design = model.matrix(~ condition + sex, data = meta_limma)
)
```

Outputs:

```text
RESULTS_OSD514_ms3fix/figs/PCA_TMT_condition_batch_ms3fix_raw.png
RESULTS_OSD514_ms3fix/figs/PCA_TMT_condition_batch_ms3fix_batch_corrected.png
RESULTS_OSD514_ms3fix/tables/TMT_expression_matrix_ms3fix_batch_corrected.csv
TMT_expression_matrix_ms3fix_batch_corrected.csv
```

The batch-corrected matrix is used for visualization/export. Differential
expression is modeled with limma design matrices that include batch terms rather
than using the batch-corrected matrix directly.

## 6. limma differential expression

The R script tests three condition-level contrasts:

- `SF1g_vs_Earth`
- `SFug_vs_Earth`
- `SF1g_vs_SFug`

It runs three model specifications:

```text
condition_only       ~0 + condition
condition_batch      ~0 + condition + tmt_run
condition_batch_sex  ~0 + condition + tmt_run + sex
```

For each model and contrast, the script writes:

- the limma design matrix
- the model metadata
- the full limma result table from `topTable`
- significant proteins
- a volcano plot
- a per-model significant-count summary

Model-specific outputs are stored under:

```text
RESULTS_OSD514_ms3fix/tables/model_outputs/
RESULTS_OSD514_ms3fix/figs/model_outputs/
```

Model-output significance uses:

```text
adj.P.Val < 0.05 and abs(logFC) >= 0.5
```

The primary model for top-level results is:

```text
condition_batch_sex
```

Primary-model tables are copied to:

```text
RESULTS_OSD514_ms3fix/tables/Limma_SF1g_vs_Earth_results.csv
RESULTS_OSD514_ms3fix/tables/Limma_SFug_vs_Earth_results.csv
RESULTS_OSD514_ms3fix/tables/Limma_SF1g_vs_SFug_results.csv
RESULTS_OSD514_ms3fix/tables/Limma_SF1g_vs_Earth_significant.csv
RESULTS_OSD514_ms3fix/tables/Limma_SFug_vs_Earth_significant.csv
RESULTS_OSD514_ms3fix/tables/Limma_SF1g_vs_SFug_significant.csv
```

The top-level primary-model significant files use a stricter fold-change cutoff:

```text
adj.P.Val < 0.05 and abs(logFC) >= 1
```

## 7. Sex-stratified limma models

The R script also runs sex-stratified contrasts with a group model:

```text
~0 + group + tmt_run
```

Contrasts:

- `SFug_females_vs_Earth_females`
- `SFug_males_vs_Earth_males`
- `SF1g_females_vs_Earth_females`
- `SF1g_males_vs_Earth_males`

Outputs:

```text
RESULTS_OSD514_ms3fix/tables/model_outputs/sex_stratified_condition_batch/
RESULTS_OSD514_ms3fix/figs/model_outputs/sex_stratified_condition_batch/
RESULTS_OSD514_ms3fix/tables/sex_stratified_significant_counts_summary.csv
RESULTS_OSD514_ms3fix/tables/sex_stratified_significant_counts_direction_summary.csv
```

The sex-stratified output includes both FDR-only significant tables and
FDR-plus-fold-change significant tables.

## 8. Heatmap

The heatmap uses the top 50 proteins from the primary `SF1g_vs_Earth` limma
results, ordered by adjusted p-value. Values are row-scaled z-scores from the
raw log2 expression matrix, with columns ordered by condition, sex, and TMT run.

Output:

```text
RESULTS_OSD514_ms3fix/figs/Heatmap_Top50_ms3fix.png
```

## 9. GSEA and FlyEnrichr enrichment

The R script runs GO enrichment analyses for:

- `SF1g_vs_Earth`
- `SFug_vs_Earth`
- `SF1g_vs_SFug`

Ranked gene lists are built from limma results as:

```text
sign(logFC) * -log10(P.Value)
```

The script then:

1. maps protein IDs through `org.Dm.eg.db` using `UNIPROT`
2. builds GO Biological Process gene sets
3. runs `fgseaMultilevel`
4. saves fgsea result tables and plots
5. selects significant proteins using `adj.P.Val < 0.05` and `abs(logFC) >= 1`
6. maps UniProt IDs to fly gene symbols
7. submits those gene lists to FlyEnrichr GO libraries
8. writes FlyEnrichr tables and biological-process network plots

Outputs:

```text
RESULTS_OSD514_ms3fix/GSEA/tables/fgsea_SF1g_vs_Earth.csv
RESULTS_OSD514_ms3fix/GSEA/tables/fgsea_SFug_vs_Earth.csv
RESULTS_OSD514_ms3fix/GSEA/tables/fgsea_SF1g_vs_SFug.csv
RESULTS_OSD514_ms3fix/GSEA/tables/*GO*.csv
RESULTS_OSD514_ms3fix/GSEA/figs/fgsea_*.png
RESULTS_OSD514_ms3fix/GSEA/figs/enrichment_*.png
RESULTS_OSD514_ms3fix/GSEA/figs/network_*.png
```

FlyEnrichr requires network access when rerunning the enrichment step.

## Reproducing `RESULTS_OSD514_ms3fix`

After the FragPipe stage has produced `TMT_all_from_psm_pdlike_ms3fix/`, run:

```bash
Rscript GEN_PROTEOMICS_SCRIPT_ms3fix.r
```

Before rerunning on a new machine, update the `BASE_DIR` value near the top of
`GEN_PROTEOMICS_SCRIPT_ms3fix.r`. The current script also installs missing CRAN
and Bioconductor packages automatically.

The script writes:

```text
RESULTS_OSD514_ms3fix/tables/
RESULTS_OSD514_ms3fix/figs/
RESULTS_OSD514_ms3fix/GSEA/
```

## Notes on retained output files

Some files in `RESULTS_OSD514_ms3fix/tables/` have `_ms3fix` in the filename
and some do not. The current main R script writes the expression-matrix files
with `_ms3fix` names, while the primary limma result files are written without
the suffix and are copied from the `condition_batch_sex` model output. The
`_ms3fix`-suffixed limma tables are retained analysis artifacts from the same
MS3-fixed workflow.

## Key scripts

- `fragpipe_dmel_tmt/run_fragpipe_dmel_tmt.sh`: end-to-end FragPipe TMT runner.
- `fragpipe_dmel_tmt/build_tmt_annotation.py`: creates per-plex TMT reporter
  annotation files from the assay metadata.
- `fragpipe_dmel_tmt/export_tmt_all_like_from_psm.py`: exports FragPipe PSM
  reporter intensities into PD-like per-sample protein abundance tables.
- `fragpipe_dmel_tmt/extract_tmt_protein_matrix.py`: extracts a sample matrix
  from FragPipe TMT protein report output.
- `GEN_PROTEOMICS_SCRIPT_ms3fix.r`: builds the final expression matrix, runs
  limma models, generates plots, and runs GO enrichment.
