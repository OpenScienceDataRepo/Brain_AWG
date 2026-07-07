# syn2580853 — Emory Drosophila Tau Model

Drosophila melanogaster model of Tau-driven neurodegeneration (EmoryTau / AMP-AD). This directory holds analysis code and documentation for the Synapse project **syn2580853**.

## Study Overview

Human Tau (wild-type and R406W FTLD mutant) is expressed in neurons via Elav-GAL4. RNA-seq and LC-MS/MS proteomics were collected at Day 1, Day 10, and Day 20.


| Genotype | Construct                   |
| -------- | --------------------------- |
| Control  | Elav-GAL4/+                 |
| TauWT    | Elav-GAL4/+; UAS-TauWT/+    |
| TauR406W | Elav-GAL4/+; UAS-TauR406W/+ |


## Subdirectories

- [RNA_Seq](RNA_Seq) — transcriptomics pipeline (`TauAnalysis.R`)
- [Proteomics](Proteomics) — proteomics pipeline (`TauProteomics.R`)

Raw data and analysis outputs are stored on Synapse and generated locally when the scripts are run; they are not committed to this repository.