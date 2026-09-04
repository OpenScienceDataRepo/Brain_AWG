"""
ven_dopaminergic_analysis.py
Analyses the DOPAMINERGIC (midbrain) organoid samples from GSE259421
(Marotta et al. 2024) using the same VEN gene panel and permutation
methodology as ven_organoid_analysis.py (cortical organoids).

Purpose: test whether VEN panel signals observed in cortical organoids
are cortex-specific (supporting VEN-pathway interpretation) or arise from
general microgravity effects on any neural organoid type.

Dopaminergic organoids do not contain VENs (which are restricted to human
ACC and frontal insula); any VEN panel signal in this dataset would imply
a general neural/microgravity effect rather than VEN-specific biology.

Subjects S3 and S4 (dopaminergic iPSC line) -- ISS 38 days vs Ground.
Sample counts match those of the cortical analysis (S1/S2).

DEPENDENCIES: pip install pandas numpy scipy matplotlib mygene
REQUIRES: GSE259421_all_counts.txt in the same directory
          results/ensembl_to_symbol.csv (from ven_organoid_analysis.py run)
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

os.makedirs("results", exist_ok=True)

#  VEN gene panel (identical to cortical analysis) 
VEN_PANEL = {
    "Myelination": ["MBP","MOG","PLP1","MAG","CNP","MOBP","ERMN"],
    "FastSignalling": ["SCN1A","KCNQ2","ANK3","NEFH","NEFM","NEFL","SNCG"],
    "SocialCircuit": ["OXTR","AVPR1A","HTR2A","DRD1","CHRM1","GABRB2"],
    "LayerVProjection": ["FEZF2","BCL11B","TBR1","SATB2","CUX1"],
    "MetabolicSupport": ["VDAC1","ATP2B2","SLC17A7","SNAP25","SYP","NRXN1"],
}
ALL_PANEL_GENES = [g for gs in VEN_PANEL.values() for g in gs]
INDIVIDUAL_MARKERS = ["NOS1", "TH", "SLC6A3"]  # TH and SLC6A3 = dopaminergic markers

#  Sample metadata (full 42 samples, BAM order) 
BAM_ORDER = [
    "a1","a13","a14","a15","a16","a2","a23","a24","a25","a29",
    "a30","a31","a36","a37","a38","a42","a43","a48","a49","a7",
    "a8","a9","b1","b15","b16","b17","b18","b2","b25","b26",
    "b3","b31","b33","b39","b40","b41","b44","b45","b50","b51",
    "b8","b9"
]

SAMPLE_METADATA = [
    ("S1","Cortical","NoMicroglia","Ground"),
    ("S1","Cortical","NoMicroglia","Ground"),
    ("S1","Cortical","NoMicroglia","Ground"),
    ("S1","Cortical","WithMicroglia","Ground"),
    ("S1","Cortical","WithMicroglia","Ground"),
    ("S1","Cortical","WithMicroglia","Ground"),
    ("S2","Cortical","NoMicroglia","Ground"),
    ("S2","Cortical","NoMicroglia","Ground"),
    ("S2","Cortical","WithMicroglia","Ground"),
    ("S2","Cortical","WithMicroglia","Ground"),
    ("S1","Cortical","NoMicroglia","LEO"),
    ("S1","Cortical","NoMicroglia","LEO"),
    ("S1","Cortical","WithMicroglia","LEO"),
    ("S1","Cortical","WithMicroglia","LEO"),
    ("S1","Cortical","WithMicroglia","LEO"),
    ("S2","Cortical","NoMicroglia","LEO"),
    ("S2","Cortical","NoMicroglia","LEO"),
    ("S2","Cortical","WithMicroglia","LEO"),
    ("S2","Cortical","WithMicroglia","LEO"),
    ("S3","Dopaminergic","NoMicroglia","Ground"),
    ("S3","Dopaminergic","NoMicroglia","Ground"),
    ("S3","Dopaminergic","WithMicroglia","Ground"),
    ("S3","Dopaminergic","WithMicroglia","Ground"),
    ("S3","Dopaminergic","WithMicroglia","Ground"),
    ("S4","Dopaminergic","NoMicroglia","Ground"),
    ("S4","Dopaminergic","NoMicroglia","Ground"),
    ("S4","Dopaminergic","NoMicroglia","Ground"),
    ("S4","Dopaminergic","NoMicroglia","Ground"),
    ("S4","Dopaminergic","WithMicroglia","Ground"),
    ("S4","Dopaminergic","WithMicroglia","Ground"),
    ("S4","Dopaminergic","WithMicroglia","Ground"),
    ("S3","Dopaminergic","NoMicroglia","LEO"),
    ("S3","Dopaminergic","NoMicroglia","LEO"),
    ("S3","Dopaminergic","NoMicroglia","LEO"),
    ("S3","Dopaminergic","WithMicroglia","LEO"),
    ("S3","Dopaminergic","WithMicroglia","LEO"),
    ("S4","Dopaminergic","NoMicroglia","LEO"),
    ("S4","Dopaminergic","NoMicroglia","LEO"),
    ("S4","Dopaminergic","NoMicroglia","LEO"),
    ("S4","Dopaminergic","NoMicroglia","LEO"),
    ("S4","Dopaminergic","WithMicroglia","LEO"),
    ("S4","Dopaminergic","WithMicroglia","LEO"),
]

meta_df = pd.DataFrame(
    SAMPLE_METADATA,
    columns=["Subject","OrgType","Microglia","Condition"],
    index=[b + ".bam" for b in BAM_ORDER]
)

#  Load counts 
for fname in ["GSE259421_all_counts.txt","GSE259421_all_counts.txt.gz"]:
    if os.path.exists(fname):
        counts_file = fname
        break
else:
    sys.exit("ERROR: GSE259421_all_counts.txt not found.")

print(f"\n{'='*60}")
print(f"Loading: {counts_file}")
raw = pd.read_csv(counts_file, sep="\t", comment="#", index_col=0)
bam_cols = [c for c in raw.columns if c.endswith(".bam")]
raw = raw[bam_cols]
raw.columns = [c.split("/")[-1] for c in raw.columns]

if raw.index[0].startswith("ENSG"):
    raw.index = raw.index.str.replace(r"\.\d+$","", regex=True)

#  Load cached gene symbol mapping (must have run ven_organoid_analysis.py first)
MAPPING_CACHE = "results/ensembl_to_symbol.csv"
if not os.path.exists(MAPPING_CACHE):
    print("NOTE: results/ensembl_to_symbol.csv not found.")
    print("Run ven_organoid_analysis.py first to build the symbol mapping cache.")
    try:
        import mygene
        print("Falling back to live MyGene.info query...")
        mg = mygene.MyGeneInfo()
        result = mg.querymany(raw.index.tolist(), scopes="ensembl.gene",
                              fields="symbol", species="human",
                              as_dataframe=True, verbose=False)
        symbol_map = result["symbol"].dropna().to_dict()
    except ImportError:
        sys.exit("ERROR: mygene not installed. Run: pip install mygene")
else:
    cached = pd.read_csv(MAPPING_CACHE).set_index("EnsemblID")["Symbol"].to_dict()
    coverage = len(set(raw.index) & set(cached)) / len(raw.index)
    print(f"Gene symbol cache coverage: {coverage:.1%}")
    symbol_map = cached

raw.index = [symbol_map.get(g, g) for g in raw.index]
raw = raw.groupby(raw.index).sum()

# normalise 
cpm = raw.div(raw.sum(axis=0), axis=1) * 1_000_000
log2cpm = np.log2(cpm + 1)

#   Sample selection helpers  
def select_dopaminergic(microglia_group, condition):
    mask = (
        (meta_df["OrgType"] == "Dopaminergic") &
        (meta_df["Microglia"] == microglia_group) &
        (meta_df["Condition"] == condition)
    )
    return [c for c in meta_df[mask].index if c in log2cpm.columns]

print("\nDopaminergic organoid sample counts:")
for mic in ["NoMicroglia","WithMicroglia"]:
    for cond in ["Ground","LEO"]:
        n = len(select_dopaminergic(mic, cond))
        print(f"  {mic:15s} {cond:6s}: n = {n}")

def compute_log2fc(leo_cols, gnd_cols):
    return log2cpm[leo_cols].mean(axis=1) - log2cpm[gnd_cols].mean(axis=1)

fc_groups = {
    "NoMicroglia": compute_log2fc(
        select_dopaminergic("NoMicroglia","LEO"),
        select_dopaminergic("NoMicroglia","Ground")
    ),
    "WithMicroglia": compute_log2fc(
        select_dopaminergic("WithMicroglia","LEO"),
        select_dopaminergic("WithMicroglia","Ground")
    ),
}
fc_groups["Combined"] = compute_log2fc(
    select_dopaminergic("NoMicroglia","LEO") + select_dopaminergic("WithMicroglia","LEO"),
    select_dopaminergic("NoMicroglia","Ground") + select_dopaminergic("WithMicroglia","Ground")
)

#  VEN panel analysis 
N_PERM  = 10_000
RNG_SEED = 42
rng = np.random.default_rng(RNG_SEED)

def run_panel_analysis(fc_series, panel_dict, n_perm=N_PERM):
    background_fc = fc_series.dropna().values
    rows = []
    for category, gene_list in panel_dict.items():
        present = [g for g in gene_list if g in fc_series.index
                   and not np.isnan(fc_series[g])]
        if not present:
            rows.append({"Category":category,"n_panel":len(gene_list),
                         "n_present":0,"mean_log2FC":np.nan,
                         "t_stat":np.nan,"p_ttest":np.nan,
                         "SD_above_null":np.nan,"perm_p":np.nan,"sig":"",
                         "genes_present":"none found"})
            continue
        category_fc = fc_series[present].values
        category_mean = category_fc.mean()
        t_stat, p_ttest = (stats.ttest_1samp(category_fc, 0)
                           if len(present) > 1 else (np.nan, np.nan))
        null_means = np.array([
            rng.choice(background_fc, size=len(present), replace=False).mean()
            for _ in range(n_perm)
        ])
        null_mean = null_means.mean()
        null_sd = null_means.std()
        sd_above  = (category_mean - null_mean) / null_sd if null_sd > 0 else np.nan
        perm_p = np.mean(np.abs(null_means - null_mean) >= np.abs(category_mean - null_mean))
        if perm_p < 0.001: sig = "***"
        elif perm_p < 0.01: sig = "**"
        elif perm_p < 0.05: sig = "*"
        elif perm_p < 0.10: sig = "+"
        else: sig = ""
        rows.append({
            "Category": category,
            "n_panel": len(gene_list),
            "n_present": len(present),
            "mean_log2FC": round(category_mean, 4),
            "t_stat": round(t_stat, 3) if not np.isnan(t_stat) else np.nan,
            "p_ttest": round(p_ttest, 4) if not np.isnan(p_ttest) else np.nan,
            "SD_above_null": round(sd_above, 3),
            "perm_p": round(perm_p, 4),
            "sig": sig,
            "genes_present": ", ".join(present),
        })
    return pd.DataFrame(rows)

panel_results = {}
for grp, fc in fc_groups.items():
    print(f"\n{'='*60}\nGROUP: {grp}\n{'='*60}")
    df = run_panel_analysis(fc, VEN_PANEL)
    panel_results[grp] = df
    print(df[["Category","n_present","mean_log2FC","p_ttest",
              "SD_above_null","perm_p","sig"]].to_string(index=False))
    df.to_csv(f"results/dopaminergic_VEN_panel_{grp}.csv", index=False)

#  Individual marker genes (including dopaminergic markers TH, SLC6A3) 
print(f"\n{'='*60}\nINDIVIDUAL MARKER GENES (including dopaminergic markers)\n{'='*60}")
individual_genes = ALL_PANEL_GENES + INDIVIDUAL_MARKERS
marker_rows = []
for grp, fc in fc_groups.items():
    for gene in individual_genes:
        val = float(fc[gene]) if gene in fc.index else np.nan
        marker_rows.append({"Group":grp,"Gene":gene,"log2FC":val})
marker_df = pd.DataFrame(marker_rows)
try:
    pivot = marker_df.pivot(index="Gene",columns="Group",values="log2FC").round(3)
    print(pivot.to_string())
    pivot.to_csv("results/dopaminergic_individual_genes_pivot.csv")
except Exception:
    print(marker_df.to_string())
marker_df.to_csv("results/dopaminergic_individual_genes_all.csv", index=False)

print("CROSS-ORGANOID-TYPE COMPARISON: CORTICAL vs DOPAMINERGIC")
print("(Combined group, SD above genome-wide null)")

CORTICAL_CSV = "results/VEN_panel_Combined.csv"
if os.path.exists(CORTICAL_CSV):
    cort = pd.read_csv(CORTICAL_CSV).set_index("Category")
    dopa = panel_results["Combined"].set_index("Category")
    print(f"\n  {'Category':20s} {'Cortical SD':>12s} {'Cort perm-p':>12s} "
          f"{'Dopamin. SD':>12s} {'Dopa perm-p':>12s}")
    print("  " + "-"*72)
    for cat in ["Myelination","FastSignalling","SocialCircuit",
                "LayerVProjection","MetabolicSupport"]:
        c_sd = float(cort.loc[cat,"SD_above_null"]) if cat in cort.index else np.nan
        c_pp = float(cort.loc[cat,"perm_p"]) if cat in cort.index else np.nan
        d_sd = float(dopa.loc[cat,"SD_above_null"])  if cat in dopa.index else np.nan
        d_pp = float(dopa.loc[cat,"perm_p"]) if cat in dopa.index else np.nan
        def fmt(v): return f"{v:+.2f}" if not np.isnan(v) else "N/A"
        def fmtp(v): return f"{v:.4f}"  if not np.isnan(v) else "N/A"
        print(f" {cat:20s} {fmt(c_sd):>12s} {fmtp(c_pp):>12s} "
              f"{fmt(d_sd):>12s} {fmtp(d_pp):>12s}")
    print()
    print(" Interpretation guide:")
    print(" - Signal in Cortical but NOT Dopaminergic -> cortex-specific (supports VEN hypothesis)")
    print(" - Signal in BOTH -> general microgravity effect on any neural organoid")
    print(" - Signal in Dopaminergic but NOT Cortical -> unexpected, warrants investigation")
else:
    print(" NOTE: results/VEN_panel_Combined.csv not found.")
    print(" Run ven_organoid_analysis.py first for cross-organoid comparison.")

def make_comparison_figure():
    if not os.path.exists(CORTICAL_CSV):
        return
    cort = pd.read_csv(CORTICAL_CSV).set_index("Category")
    dopa = panel_results["Combined"].set_index("Category")

    CAT_ORDER = ["Myelination","FastSignalling","SocialCircuit",
                 "LayerVProjection","MetabolicSupport"]
    CAT_LABELS = ["Myelin-\nation","Fast\nSignal.","Social\nCircuit",
                  "Layer V\nProj.","Metabolic\nSupport"]
    CAT_COLORS = {
        "Myelination": "#4472C4",
        "FastSignalling": "#70AD47",
        "SocialCircuit": "#ED7D31",
        "LayerVProjection": "#FFC000",
        "MetabolicSupport": "#7030A0",
    }
    color_list = [CAT_COLORS[c] for c in CAT_ORDER]

    def get_vals(df, col):
        return [float(df.loc[c,col]) if c in df.index else np.nan for c in CAT_ORDER]

    c_sd = get_vals(cort, "SD_above_null")
    d_sd = get_vals(dopa, "SD_above_null")
    c_pp = get_vals(cort, "perm_p")
    d_pp = get_vals(dopa, "perm_p")

    x = np.arange(len(CAT_ORDER))
    w = 0.32

    fig, ax = plt.subplots(figsize=(10, 5.5))

    ax.bar(x - w/2, c_sd, w, color=color_list, alpha=0.92,
           edgecolor="white", lw=0.5, label="Cortical organoids (S1/S2)")
    ax.bar(x + w/2, d_sd, w, color=color_list, alpha=0.50,
           edgecolor="black", lw=0.8, hatch="//",
           label="Dopaminergic organoids (S3/S4)")

    ax.axhline(0,     color="black",   lw=0.9)
    ax.axhline( 1.96, color="#888888", lw=0.9, ls="--", alpha=0.65,
               label="p=0.05 (+-1.96 SD)")
    ax.axhline(-1.96, color="#888888", lw=0.9, ls="--", alpha=0.65)

    def sig_star(p):
        if np.isnan(p): return ""
        if p < 0.001: return "***"
        if p < 0.01: return "**"
        if p < 0.05: return "*"
        if p < 0.10: return "+"
        return ""

    ylo = min([v for v in c_sd + d_sd if not np.isnan(v)]) - 0.4
    yhi = max([v for v in c_sd + d_sd if not np.isnan(v)]) + 1.4
    ax.set_ylim(ylo, yhi)

    for i, (cs, ds, cp, dp) in enumerate(zip(c_sd, d_sd, c_pp, d_pp)):
        star_c = sig_star(cp)
        star_d = sig_star(dp)
        pad = (yhi - ylo) * 0.025
        if star_c and not np.isnan(cs):
            ax.text(i - w/2, cs + pad, star_c,
                    ha="center", va="bottom", fontsize=9, fontweight="bold")
        if star_d and not np.isnan(ds):
            ax.text(i + w/2, ds + pad, star_d,
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(CAT_LABELS, fontsize=9, linespacing=1.3)
    ax.set_ylabel("SD above genome-wide permutation null\n"
                  "(N=10,000 random gene sets, seed 42)", fontsize=9)
    ax.set_title(
        "VEN Gene Panel Specificity: Cortical vs Dopaminergic Organoids\n"
        "GSE259421 (ISS 38 days vs Ground) -- Combined group",
        fontsize=10, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8.5, frameon=True, framealpha=0.92)
    ax.yaxis.grid(True, alpha=0.22, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig("results/dopaminergic_vs_cortical_comparison.png",
                dpi=200, bbox_inches="tight")
    plt.close()
    print("\nSaved: results/dopaminergic_vs_cortical_comparison.png")

make_comparison_figure()

print(" SUMMARY (Combined, Dopaminergic organoids)")
summary = panel_results["Combined"][[
    "Category","n_present","mean_log2FC","t_stat","p_ttest",
    "SD_above_null","perm_p","sig","genes_present"
]].copy()
summary.columns = ["Category","N genes","Mean log2FC","t","p (t-test)",
                   "SD above null","perm p","","Genes detected"]
print(summary.to_string(index=False))
summary.to_csv("results/dopaminergic_SUMMARY_combined.csv", index=False)

metadata = {
    "dataset": "GSE259421",
    "organoid_type": "Dopaminergic (midbrain iPSC-derived)",
    "subjects": "S3, S4",
    "analysis_script": "ven_dopaminergic_analysis.py",
    "permutation_n": N_PERM,
    "random_seed": RNG_SEED,
    "purpose": "Cortex-specificity control, do VEN panel signals arise "
               "in non-VEN-containing organoid type under same spaceflight conditions?"
}
with open("results/dopaminergic_analysis_metadata.json","w") as f:
    json.dump(metadata, f, indent=2)

print(f"\n{'='*60}")
print("ANALYSIS COMPLETE Dopaminergic organoids (GSE259421)")
print(f"{'='*60}")
print("Outputs in results/:")
for f in sorted(os.listdir("results")):
    if "dopaminergic" in f.lower():
        size = os.path.getsize(f"results/{f}")
        print(f" {f} ({size:,} bytes)")
