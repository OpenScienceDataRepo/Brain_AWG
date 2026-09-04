#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User-configurable paths
# -----------------------------
BASE_DIR="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/fragpipe_dmel_tmt"
MZML_DIR="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/mzML"
TMT_ALL_DIR="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/TMT_all"
ASSAY_LABEL_FILE="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/a_OSD-514_protein-expression-profiling_mass-spectrometry_Orbitrap Fusion.txt"
TOOLS_DIR="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/fragpipe_tools"
SIF="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/brain_cDNA_discovery/singularity_containers/fragpipe_latest.sif"
THREADS="${THREADS:-24}"
WORKFLOW_TEMPLATE="${WORKFLOW_TEMPLATE:-TMT10-MS3.workflow}"

# Reuse FASTA from LFQ setup if already present.
FALLBACK_DB="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/fragpipe_dmel_lfq/db/dmel_UP000000803_uniprot_target_decoy.fasta"

MANIFEST_DIR="${BASE_DIR}/manifests"
DB_DIR="${BASE_DIR}/db"
WORK_DIR="${BASE_DIR}/work"
LOG_DIR="${BASE_DIR}/logs"
RUN_SUBDIR="${RUN_SUBDIR:-results}"
OUT_DIR="${WORK_DIR}/${RUN_SUBDIR}"
CACHE_BIND_DIR="${WORK_DIR}/fragpipe_cache"
JOBS_BIND_DIR="${WORK_DIR}/fragpipe_jobs"

TARGET_DECOY_FASTA="${DB_DIR}/dmel_UP000000803_uniprot_target_decoy.fasta"
MANIFEST_FILE="${MANIFEST_DIR}/brain_awg_tmt.fp-manifest"
WORKFLOW_FILE="${MANIFEST_DIR}/brain_awg_dmel_tmt10.workflow"

mkdir -p "${MANIFEST_DIR}" "${DB_DIR}" "${WORK_DIR}" "${OUT_DIR}" "${LOG_DIR}" "${CACHE_BIND_DIR}" "${JOBS_BIND_DIR}"

if [[ ! -f "${SIF}" ]]; then
  echo "ERROR: FragPipe container not found: ${SIF}" >&2
  exit 1
fi
if [[ ! -d "${TOOLS_DIR}" ]]; then
  echo "ERROR: FragPipe tools folder not found: ${TOOLS_DIR}" >&2
  exit 1
fi

# 1) FASTA
if [[ -s "${TARGET_DECOY_FASTA}" ]]; then
  echo "[1/6] Using existing FASTA: ${TARGET_DECOY_FASTA}"
elif [[ -s "${FALLBACK_DB}" ]]; then
  echo "[1/6] Copying FASTA from LFQ setup"
  cp -f "${FALLBACK_DB}" "${TARGET_DECOY_FASTA}"
else
  echo "ERROR: target-decoy FASTA not found. Expected one of:" >&2
  echo "  - ${TARGET_DECOY_FASTA}" >&2
  echo "  - ${FALLBACK_DB}" >&2
  exit 1
fi

# 2) Manifest: all mzML files grouped by plex (TMTA/TMTB/TMTC)
echo "[2/6] Generating TMT manifest..."
python3 - << 'PY'
from pathlib import Path
import re

mzml_dir = Path('/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/mzML')
out = Path('/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/fragpipe_dmel_tmt/manifests/brain_awg_tmt.fp-manifest')

rows = []
for p in sorted(mzml_dir.glob('*.mzML')):
    n = p.name.upper()
    if 'TMTA' in n:
        plex = 'TMTA'
    elif 'TMTB' in n:
        plex = 'TMTB'
    elif 'TMTC' in n:
        plex = 'TMTC'
    else:
        continue
    rows.append((str(p), plex))

if not rows:
    raise SystemExit('No TMTA/TMTB/TMTC mzML files found')

with out.open('w') as w:
    for p, plex in rows:
        w.write(f"{p}\t{plex}\n")

print(f"Wrote {len(rows)} manifest rows -> {out}")
PY

# 3) Build per-plex TMT annotation files in mzML folder (FragPipe expects <EXP>_annotation.txt)
echo "[3/6] Generating TMT annotation file..."
# Clear old annotation files to avoid parser ambiguity
find "${MZML_DIR}" -maxdepth 1 -type f -name '*annotation.txt' -delete

python3 "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/fragpipe_dmel_tmt/build_tmt_annotation.py" \
  --tmt-all-dir "${TMT_ALL_DIR}" \
  --assay-file "${ASSAY_LABEL_FILE}" \
  --mzml-dir "${MZML_DIR}"

# Ensure expected per-plex annotation files exist.
for plex in TMTA TMTB TMTC; do
  if [[ ! -s "${MZML_DIR}/${plex}_annotation.txt" ]]; then
    echo "ERROR: missing annotation file ${MZML_DIR}/${plex}_annotation.txt" >&2
    exit 1
  fi

  # Headless FragPipe expects annotation paths under workdir/<EXP>/<EXP>_annotation.txt.
  mkdir -p "${OUT_DIR}/${plex}"
  cp -f "${MZML_DIR}/${plex}_annotation.txt" "${OUT_DIR}/${plex}/${plex}_annotation.txt"
done

extra_ann_count=$(find "${MZML_DIR}" -maxdepth 1 -type f -name '*annotation.txt' | wc -l)
if [[ "${extra_ann_count}" -lt 3 ]]; then
  echo "ERROR: expected at least 3 per-plex annotation files in ${MZML_DIR}, found ${extra_ann_count}" >&2
  find "${MZML_DIR}" -maxdepth 1 -type f -name '*annotation.txt' -print >&2 || true
  exit 1
fi

# 4) Build workflow from FragPipe TMT template and patch database/ref tag
echo "[4/7] Creating workflow from FragPipe ${WORKFLOW_TEMPLATE} template..."
singularity exec "${SIF}" bash -lc \
  "cat /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/workflows/${WORKFLOW_TEMPLATE}" > "${WORKFLOW_FILE}"

python3 - << 'PY'
from pathlib import Path
wf = Path('/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/fragpipe_dmel_tmt/manifests/brain_awg_dmel_tmt10.workflow')
db = '/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/fragpipe_dmel_tmt/db/dmel_UP000000803_uniprot_target_decoy.fasta'
lines = wf.read_text().splitlines()
out = []
set_db = False
set_ref = False
for line in lines:
    if line.startswith('database.db-path='):
        out.append(f'database.db-path={db}')
        set_db = True
    elif line.startswith('tmtintegrator.ref_tag='):
        out.append('tmtintegrator.ref_tag=ctrl_pool')
        set_ref = True
    elif line.startswith('tmtintegrator.ms1_int='):
        out.append('tmtintegrator.ms1_int=false')
    elif line.startswith('tmtintegrator.log2transformed='):
        out.append('tmtintegrator.log2transformed=false')
    elif line.startswith('tmtintegrator.min_percent='):
        out.append('tmtintegrator.min_percent=0.0')
    elif line.startswith('tmtintegrator.max_pep_prob_thres='):
        out.append('tmtintegrator.max_pep_prob_thres=1.0')
    elif line.startswith('tmtintegrator.min_pep_prob='):
        out.append('tmtintegrator.min_pep_prob=0.0')
    else:
        out.append(line)
if not set_db:
    out.insert(0, f'database.db-path={db}')
if not set_ref:
    out.append('tmtintegrator.ref_tag=ctrl_pool')
# Ensure reporter-ion quant mode settings are present even if missing in template.
if not any(l.startswith('tmtintegrator.ms1_int=') for l in out):
    out.append('tmtintegrator.ms1_int=false')
if not any(l.startswith('tmtintegrator.log2transformed=') for l in out):
    out.append('tmtintegrator.log2transformed=false')
if not any(l.startswith('tmtintegrator.min_percent=') for l in out):
    out.append('tmtintegrator.min_percent=0.0')
if not any(l.startswith('tmtintegrator.max_pep_prob_thres=') for l in out):
    out.append('tmtintegrator.max_pep_prob_thres=1.0')
if not any(l.startswith('tmtintegrator.min_pep_prob=') for l in out):
    out.append('tmtintegrator.min_pep_prob=0.0')
wf.write_text('\n'.join(out) + '\n')
print(f'Patched workflow: {wf}')
PY

# 5) Run FragPipe headless
# Keep MZML + annotation path bound and set locale for MSFragger runtime stability.
echo "[5/7] Running FragPipe TMT headless (this can take many hours)..."
set -x
singularity exec \
  --cleanenv \
  --bind "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/" \
  --bind "${CACHE_BIND_DIR}:/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/cache" \
  --bind "${JOBS_BIND_DIR}:/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/jobs" \
  --env "LC_ALL=C.utf8" \
  --env "LANG=C.utf8" \
  "${SIF}" \
  /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/bin/fragpipe --headless \
    --workflow "${WORKFLOW_FILE}" \
    --manifest "${MANIFEST_FILE}" \
    --workdir "${OUT_DIR}" \
    --config-tools-folder "${TOOLS_DIR}" \
    --threads "${THREADS}" \
    |& tee "${LOG_DIR}/fragpipe_tmt_run.log"
set +x

# 6) Extract a clean protein matrix from TMT report
echo "[6/7] Extracting protein matrix from TMT output..."
python3 "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/fragpipe_dmel_tmt/extract_tmt_protein_matrix.py" \
  --workdir "${OUT_DIR}" \
  --out "${WORK_DIR}/protein_abundance_tmt_matrix.tsv"

# 7) Export per-sample files in TMT_all-like format
echo "[7/7] Exporting TMT_all-like per-sample txt files..."
TMT_ALL_EXPORT_DIR="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/TMT_all_from_psm_pdlike_${RUN_SUBDIR}"
python3 "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/brain_awg/proteomics/fragpipe_dmel_tmt/export_tmt_all_like_from_psm.py" \
  --results-dir "${OUT_DIR}" \
  --assay-file "${ASSAY_LABEL_FILE}" \
  --out-dir "${TMT_ALL_EXPORT_DIR}"

echo "Done. Outputs:"
echo "  - FragPipe result dir: ${OUT_DIR}"
echo "  - Run log: ${LOG_DIR}/fragpipe_tmt_run.log"
echo "  - Protein abundance matrix: ${WORK_DIR}/protein_abundance_tmt_matrix.tsv"
echo "  - TMT_all-like files (PD-like PSM aggregation): ${TMT_ALL_EXPORT_DIR}"
