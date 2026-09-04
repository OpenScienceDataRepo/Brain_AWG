# FragPipe D. melanogaster TMT configuration (brain_awg)

This folder contains a reproducible **TMT-focused** FragPipe run for your `mzML` files.

It is designed to avoid the high zero inflation you saw from per-file LFQ by using TMT quantification (`TMT10.workflow` + `TMT-Integrator`).

## Upstream quant coverage fix (important)

If protein coverage is very sparse (e.g., only ~500 proteins), use an **MS3 reporter-ion workflow**:

- Workflow template: `TMT10-MS3.workflow` (not `TMT10.workflow`)
- `tmtintegrator.ms1_int=false` (reporter-ion mode)
- `tmtintegrator.log2transformed=false`
- Relax summarization gates for recovery:
  - `tmtintegrator.min_pep_prob=0.0`
  - `tmtintegrator.max_pep_prob_thres=1.0`
  - `tmtintegrator.min_percent=0.0`

The runner now defaults to `WORKFLOW_TEMPLATE=TMT10-MS3.workflow` and patches these settings automatically.

## What this run does

1. Reuses your D. melanogaster target-decoy FASTA.
2. Builds a TMT manifest grouping files by plex (`TMTA`, `TMTB`, `TMTC`).
3. Auto-generates TMT experiment annotation from assay metadata (`a_OSD-514_protein-expression-profiling_mass-spectrometry_Orbitrap Fusion.txt`).
4. Runs FragPipe headless using built-in `TMT10` workflow template.
5. Extracts a clean protein abundance matrix from the TMT output.
6. Exports per-sample `TMT_all`-like text files using PD-like aggregation from `psm.tsv`.

## Files

- `run_fragpipe_dmel_tmt.sh` – end-to-end TMT runner.
- `build_tmt_annotation.py` – builds per-plex annotation files from assay metadata (fallback: `TMT_all` headers).
- `extract_tmt_protein_matrix.py` – extracts protein matrix from FragPipe TMT report.
- `export_tmt_all_like_from_psm.py` – high-coverage PD-like export from `psm.tsv` (recommended).
- `export_tmt_all_like.py` – legacy export from `tmt-report` table (lower coverage, not recommended for this dataset).

## Expected outputs

- Main FragPipe output directory:
  - `work/results/`
- FragPipe run log:
  - `logs/fragpipe_tmt_run.log`
- Final protein matrix for downstream DE:
  - `work/protein_abundance_tmt_matrix.tsv`
- TMT_all-like per-sample files:
  - `/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/TMT_all_from_psm_pdlike_<RUN_SUBDIR>/`

## How the authors likely processed `TMT_all`

Based on `i_Investigation.txt` and assay metadata, the published `TMT_all` files were produced with:

- Thermo **Proteome Discoverer 2.1**
- **SequestHT** database search
- Reporter Ions Quantifier with TMT10 settings (MS3 multi-notch method)
- No cross-channel normalization/scaling in PD export

This repo reproduces analogous outputs with FragPipe and exports them in the same per-sample text layout.

