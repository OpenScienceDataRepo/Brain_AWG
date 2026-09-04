"""
alysson_proteomics_analysis.py
Analyses Jourdon et al. (2026) iPSC brain organoid proteomics (PXD069807)
for VEN panel protein abundance changes in WT83 ISS vs Ground.

Data exported from PatternLab V SEPro files as tab-delimited protein tables.
Input files: ISS_WT83_proteins.txt, Ground_WT83_proteins.txt
"""

import os, re, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

os.makedirs("results", exist_ok=True)

VEN_PANEL = {
    "Myelination": ["MBP","MOG","PLP1","MAG","CNP","MOBP","ERMN"],
    "FastSignalling": ["SCN1A","KCNQ2","ANK3","NEFH","NEFM","NEFL","SNCG"],
    "SocialCircuit": ["OXTR","AVPR1A","HTR2A","DRD1","CHRM1","GABRB2"],
    "LayerVProjection": ["FEZF2","BCL11B","TBR1","SATB2","CUX1"],
    "MetabolicSupport": ["VDAC1","ATP2B2","SLC17A7","SNAP25","SYP","NRXN1"],
}
ALL_PANEL_GENES = [g for gs in VEN_PANEL.values() for g in gs]
N_PERM, RNG_SEED = 10_000, 42
rng = np.random.default_rng(RNG_SEED)

def parse_patternlab_txt(filepath):
    """
    Parse a PatternLab V tab-delimited protein export.
    Returns dict: gene_symbol -> SpectrumCount (int).
    Skips header/comment lines (starting with # or non-data content).
    Protein rows: Locus | Group | Length | MolWt | SequenceCount | SpectrumCount | ...
    Gene extracted from GN= tag in Description column.
    """
    gene_counts = {}
    gene_re = re.compile(r'GN=([A-Z][A-Z0-9]{1,15})')

    with open(filepath, encoding='utf-8', errors='replace') as f:
        in_data = False
        for line in f:
            line = line.rstrip('\n')
            # Data header line starts with "#Locus"
            if line.startswith('#Locus'):
                in_data = True
                continue
            if not in_data:
                continue
            # Skip peptide sub-rows (start with tab) and comment lines
            if line.startswith('\t') or line.startswith('#') or not line.strip():
                continue

            parts = line.split('\t')
            if len(parts) < 10:
                continue

            locus = parts[0].strip()
            try:
                spec_count = int(float(parts[5].strip()))
            except (ValueError, IndexError):
                continue

            description = parts[9] if len(parts) > 9 else ''

            # Extract gene from GN= in description
            gn = gene_re.search(description)
            if gn:
                gene = gn.group(1)
            else:
                # Fall back to locus: sp|ACC|GENE_HUMAN
                m = re.search(r'\|([A-Z][A-Z0-9]+)_HUMAN', locus)
                gene = m.group(1) if m else ''

            if gene and spec_count > 0:
                # If gene appears multiple times (isoforms), sum counts
                gene_counts[gene] = gene_counts.get(gene, 0) + spec_count

    return gene_counts


for fname in ["ISS_WT83_proteins.txt", "Ground_WT83_proteins.txt"]:
    if not os.path.exists(fname):
        print(f"Missing: {fname}")
        sys.exit(1)

print("Loading ISS_WT83_proteins.txt ...")
iss_counts = parse_patternlab_txt("ISS_WT83_proteins.txt")
print(f" {len(iss_counts):,} proteins with gene names")

print("Loading Ground_WT83_proteins.txt ...")
ground_counts = parse_patternlab_txt("Ground_WT83_proteins.txt")
print(f" {len(ground_counts):,} proteins with gene names")

PSEUDO = 1
all_genes = set(iss_counts) | set(ground_counts)
log2fc = {}
for gene in all_genes:
    ic = iss_counts.get(gene, 0) + PSEUDO
    gc = ground_counts.get(gene, 0) + PSEUDO
    log2fc[gene] = np.log2(ic / gc)

log2fc_s = pd.Series(log2fc, name="log2FC_spectral")
log2fc_s.to_csv("results/alysson_all_proteins_log2FC.csv", header=True)

print(f"\nTotal proteins with log2FC: {len(log2fc_s):,}")
print(f"Genome-wide: mean={log2fc_s.mean():+.4f}  SD={log2fc_s.std():.4f}")

print("\nVEN panel gene detection (ISS spectral count / Ground spectral count):")
for cat, genes in VEN_PANEL.items():
    present = [g for g in genes if g in iss_counts or g in ground_counts]
    print(f" {cat}: {len(present)}/{len(genes)}")
    for g in genes:
        ic = iss_counts.get(g, 0); gc = ground_counts.get(g, 0)
        fc = log2fc.get(g, float('nan'))
        print(f"    {g:10s}: ISS={ic:4d}  Ground={gc:4d}  log2FC={fc:+.3f}")

background = log2fc_s.values
rows = []

print(f"\n{'-'*65}")
print(f"VEN panel permutation test (N={N_PERM:,}, seed {RNG_SEED})")
print(f"{'Category':20s} {'n':>4s} {'mean log2FC':>12s} {'perm p':>10s} {'SD above null':>14s}")

for category, gene_list in VEN_PANEL.items():
    present = [g for g in gene_list if g in log2fc_s.index]
    if not present:
        print(f"{category:20s} 0 --- --- ---")
        rows.append({"Category":category,"n_panel":len(gene_list),"n_present":0,
                     "mean_log2FC":np.nan,"SD_above_null":np.nan,"perm_p":np.nan,"sig":""})
        continue
    cat_fc = log2fc_s[present].values
    cat_mean = cat_fc.mean()
    t_stat, p_ttest = stats.ttest_1samp(cat_fc, 0) if len(present) > 1 else (np.nan, np.nan)
    null_means = np.array([
        rng.choice(background, size=len(present), replace=False).mean()
        for _ in range(N_PERM)
    ])
    null_sd  = null_means.std()
    sd_above = (cat_mean - null_means.mean()) / null_sd if null_sd > 0 else np.nan
    perm_p = np.mean(np.abs(null_means - null_means.mean()) >= abs(cat_mean - null_means.mean()))
    sig = "***" if perm_p<0.001 else ("**" if perm_p<0.01 else ("*" if perm_p<0.05 else ("+" if perm_p<0.10 else "")))
    rows.append({"Category":category,"n_panel":len(gene_list),"n_present":len(present),
                 "mean_log2FC":round(cat_mean,4),"t_stat":round(t_stat,3) if not np.isnan(t_stat) else np.nan,
                 "p_ttest":round(p_ttest,4) if not np.isnan(p_ttest) else np.nan,
                 "SD_above_null":round(sd_above,3),"perm_p":round(perm_p,4),"sig":sig,
                 "genes_present":", ".join(present)})
    print(f"{category:20s} {len(present):>4d} {cat_mean:>+12.4f} {perm_p:>10.4f}{sig} {sd_above:>14.2f}")

results_df = pd.DataFrame(rows)
results_df.to_csv("results/alysson_VEN_panel_WT83_ISS_vs_Ground.csv", index=False)
print("\nSaved: results/alysson_VEN_panel_WT83_ISS_vs_Ground.csv")

print("\nIndividual VEN marker genes:")
for gene in ALL_PANEL_GENES + ["NOS1"]:
    ic = iss_counts.get(gene, 0); gc = ground_counts.get(gene, 0)
    fc = log2fc.get(gene)
    detected = f"ISS={ic}  Ground={gc}  log2FC={fc:+.3f}" if fc is not None else "not detected"
    print(f"  {gene:10s}: {detected}")

CAT_ORDER  = ["Myelination","FastSignalling","SocialCircuit","LayerVProjection","MetabolicSupport"]
CAT_COLORS = {"Myelination":"#4472C4","FastSignalling":"#70AD47",
              "SocialCircuit":"#ED7D31","LayerVProjection":"#FFC000",
              "MetabolicSupport":"#7030A0"}
CAT_LABELS = ["Myelin-\nation","Fast\nSignal.","Social\nCircuit","Layer V\nProj.","Metabolic\nSupport"]

res_idx = results_df.set_index("Category")
sds = [float(res_idx.loc[c,"SD_above_null"]) if c in res_idx.index else np.nan for c in CAT_ORDER]
pps = [float(res_idx.loc[c,"perm_p"]) if c in res_idx.index else np.nan for c in CAT_ORDER]

fig, ax = plt.subplots(figsize=(9,5))
valid = [v for v in sds if not np.isnan(v)]
bar_vals = [v if not np.isnan(v) else 0 for v in sds]
ax.bar(range(len(CAT_ORDER)), bar_vals, color=[CAT_COLORS[c] for c in CAT_ORDER],
       edgecolor="black", lw=0.8, alpha=0.88)
ax.axhline(0, color="black",   lw=0.9)
ax.axhline( 1.96, color="#888888", lw=0.9, ls="--", alpha=0.65, label="p=0.05 (±1.96 SD)")
ax.axhline(-1.96, color="#888888", lw=0.9, ls="--", alpha=0.65)
if valid:
    ylo, yhi = min(min(valid)-0.5, -2.5), max(max(valid)+1.2, 2.5)
    ax.set_ylim(ylo, yhi)
    for i,(sd,pp) in enumerate(zip(sds,pps)):
        if np.isnan(sd) or np.isnan(pp): continue
        if pp < 0.10:
            sym = "***" if pp<0.001 else ("**" if pp<0.01 else ("*" if pp<0.05 else "+"))
            ax.text(i, sd+(yhi-ylo)*0.025, sym, ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.set_xticks(range(len(CAT_ORDER)))
ax.set_xticklabels(CAT_LABELS, fontsize=9, linespacing=1.3)
ax.set_ylabel("SD above genome-wide permutation null\n(N=10,000 random proteins, seed 42)", fontsize=9)
ax.set_title("VEN Gene Panel Alysson Lab Proteomics (PXD069807)\n"
             "WT83 ISS vs Ground (30 days) · Jourdon et al. 2026\n"
             "Spectral count quantification", fontsize=10, fontweight="bold")
ax.legend(loc="upper right", fontsize=8.5)
ax.yaxis.grid(True, alpha=0.22, lw=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("results/alysson_VEN_panel_permutation.png", dpi=200, bbox_inches="tight")
plt.close()
print("Saved: results/alysson_VEN_panel_permutation.png")
print("COMPLETE Alysson lab proteomics VEN panel analysis")
