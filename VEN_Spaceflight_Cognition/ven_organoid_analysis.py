"""
VEN Gene Panel Analysis GSE259421
Human iPSC-derived cortical organoids on the ISS (Marotta et al. 2024)

REPRODUCIBILITY NOTES
* Gene symbol mapping uses MyGene.info (mygene Python package),
  a standard, citable bioinformatics resource (Xin et al. 2016).
* VEN gene panel is pre-defined from Keskin (2026) Fast Lane Hypothesis
  and Keskin (2026) VEN Fatigue Hypothesis, NOT selected post-hoc.
* Permutation test (N=10,000) follows identical methodology to
  Keskin (2026) VEN Fatigue Hypothesis, enabling direct comparison.
* Random seed fixed at 42 for full reproducibility.
* All intermediate files saved for audit.

DEPENDENCIES: pip install pandas numpy scipy matplotlib mygene

USAGE
    Place GSE259421_all_counts.txt (or .txt.gz) in the same directory.
    Run: python ven_organoid_analysis.py
    Results saved to: ./results/

CITATION FOR GENE MAPPING
    Xin J, Mark A, Afrasiabi C, et al. High-performance web services for
    querying gene and variant annotation. Genome Biology. 2016;17:91.
    https://mygene.info
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

# SECTION 1: SAMPLE METADATA
# Sequential mapping: BAM file order → GEO sample order (GSM8115892–GSM8115933)
# Source: GSE259421 GEO page, verified manually against 42 listed samples

BAM_ORDER = [
    "a1","a13","a14","a15","a16","a2","a23","a24","a25","a29",
    "a30","a31","a36","a37","a38","a42","a43","a48","a49","a7",
    "a8","a9","b1","b15","b16","b17","b18","b2","b25","b26",
    "b3","b31","b33","b39","b40","b41","b44","b45","b50","b51",
    "b8","b9"
]

# Columns: Subject | OrgType | Microglia | Condition
# Derived from GEO sample titles GSM8115892–GSM8115933
SAMPLE_METADATA = [
    ("S1","Cortical","NoMicroglia","Ground"),   # GSM8115892
    ("S1","Cortical","NoMicroglia","Ground"),   # GSM8115893
    ("S1","Cortical","NoMicroglia","Ground"),   # GSM8115894
    ("S1","Cortical","WithMicroglia","Ground"), # GSM8115895
    ("S1","Cortical","WithMicroglia","Ground"), # GSM8115896
    ("S1","Cortical","WithMicroglia","Ground"), # GSM8115897
    ("S2","Cortical","NoMicroglia","Ground"),   # GSM8115898
    ("S2","Cortical","NoMicroglia","Ground"),   # GSM8115899
    ("S2","Cortical","WithMicroglia","Ground"), # GSM8115900
    ("S2","Cortical","WithMicroglia","Ground"), # GSM8115901
    ("S1","Cortical","NoMicroglia","LEO"),      # GSM8115902
    ("S1","Cortical","NoMicroglia","LEO"),      # GSM8115903
    ("S1","Cortical","WithMicroglia","LEO"),    # GSM8115904
    ("S1","Cortical","WithMicroglia","LEO"),    # GSM8115905
    ("S1","Cortical","WithMicroglia","LEO"),    # GSM8115906
    ("S2","Cortical","NoMicroglia","LEO"),      # GSM8115907
    ("S2","Cortical","NoMicroglia","LEO"),      # GSM8115908
    ("S2","Cortical","WithMicroglia","LEO"),    # GSM8115909
    ("S2","Cortical","WithMicroglia","LEO"),    # GSM8115910
    ("S3","Dopaminergic","NoMicroglia","Ground"),   # GSM8115911
    ("S3","Dopaminergic","NoMicroglia","Ground"),   # GSM8115912
    ("S3","Dopaminergic","WithMicroglia","Ground"), # GSM8115913
    ("S3","Dopaminergic","WithMicroglia","Ground"), # GSM8115914
    ("S3","Dopaminergic","WithMicroglia","Ground"), # GSM8115915
    ("S4","Dopaminergic","NoMicroglia","Ground"),   # GSM8115916
    ("S4","Dopaminergic","NoMicroglia","Ground"),   # GSM8115917
    ("S4","Dopaminergic","NoMicroglia","Ground"),   # GSM8115918
    ("S4","Dopaminergic","NoMicroglia","Ground"),   # GSM8115919
    ("S4","Dopaminergic","WithMicroglia","Ground"), # GSM8115920
    ("S4","Dopaminergic","WithMicroglia","Ground"), # GSM8115921
    ("S4","Dopaminergic","WithMicroglia","Ground"), # GSM8115922
    ("S3","Dopaminergic","NoMicroglia","LEO"),      # GSM8115923
    ("S3","Dopaminergic","NoMicroglia","LEO"),      # GSM8115924
    ("S3","Dopaminergic","NoMicroglia","LEO"),      # GSM8115925
    ("S3","Dopaminergic","WithMicroglia","LEO"),    # GSM8115926
    ("S3","Dopaminergic","WithMicroglia","LEO"),    # GSM8115927
    ("S4","Dopaminergic","NoMicroglia","LEO"),      # GSM8115928
    ("S4","Dopaminergic","NoMicroglia","LEO"),      # GSM8115929
    ("S4","Dopaminergic","NoMicroglia","LEO"),      # GSM8115930
    ("S4","Dopaminergic","NoMicroglia","LEO"),      # GSM8115931
    ("S4","Dopaminergic","WithMicroglia","LEO"),    # GSM8115932
    ("S4","Dopaminergic","WithMicroglia","LEO"),    # GSM8115933
]

assert len(BAM_ORDER) == len(SAMPLE_METADATA) == 42, "Metadata length mismatch"

meta_df = pd.DataFrame(
    SAMPLE_METADATA,
    columns=["Subject","OrgType","Microglia","Condition"],
    index=[b + ".bam" for b in BAM_ORDER]
)

# SECTION 2: VEN GENE PANEL
# Pre-defined in Keskin (2026) Fast Lane Hypothesis & VEN Fatigue Hypothesis.
# These gene symbols are looked up programmatically  NOT by Ensembl ID.

VEN_PANEL = {
    "Myelination": ["MBP","MOG","PLP1","MAG","CNP","MOBP","ERMN"],
    "FastSignalling": ["SCN1A","KCNQ2","ANK3","NEFH","NEFM","NEFL","SNCG"],
    "SocialCircuit": ["OXTR","AVPR1A","HTR2A","DRD1","CHRM1","GABRB2"],
    "LayerVProjection": ["FEZF2","BCL11B","TBR1","SATB2","CUX1"],
    "MetabolicSupport": ["VDAC1","ATP2B2","SLC17A7","SNAP25","SYP","NRXN1"],
}

# Additional direct VEN biochemical markers to report individually
# NOS1 is not in the standard panel but is a known VEN marker (Stimpson et al. 2011)
INDIVIDUAL_MARKERS = ["NOS1"]

ALL_PANEL_GENES = [g for genes in VEN_PANEL.values() for g in genes]

# SECTION 3: LOAD COUNTS FILE

def find_counts_file():
    for fname in ["GSE259421_all_counts.txt","GSE259421_all_counts.txt.gz"]:
        if os.path.exists(fname):
            return fname
    sys.exit("ERROR: Could not find GSE259421_all_counts.txt — place it in this directory.")

counts_file = find_counts_file()
print(f"Loading: {counts_file}")
raw = pd.read_csv(counts_file, sep="\t", comment="#", index_col=0)

# Keep only count columns (BAM files), drop annotation columns
bam_cols = [c for c in raw.columns if c.endswith(".bam")]
raw = raw[bam_cols]
raw.columns = [c.split("/")[-1] for c in raw.columns]  # strip full path

print(f"Raw matrix: {raw.shape[0]:,} genes × {raw.shape[1]} samples")
print(f"Index type: {raw.index[0]} ({'Ensembl' if raw.index[0].startswith('ENSG') else 'symbol'})")

# Strip Ensembl version suffix if present (ENSG00000xxx.12 → ENSG00000xxx)
if raw.index[0].startswith("ENSG"):
    raw.index = raw.index.str.replace(r"\.\d+$", "", regex=True)
    print("Stripped Ensembl version suffixes.")

# SECTION 4: GENE SYMBOL MAPPING
# Uses MyGene.info via the mygene Python package.
# Maps ALL Ensembl IDs in the matrix not just VEN panel genes.
# Mapping saved to results/ensembl_to_symbol.csv for full transparency.

MAPPING_CACHE = "results/ensembl_to_symbol.csv"

def build_symbol_map(ensembl_ids):
    """Query MyGene.info to map Ensembl IDs → gene symbols."""
    try:
        import mygene
    except ImportError:
        sys.exit("ERROR: mygene not installed. Run: pip install mygene")

    print(f"\nQuerying MyGene.info for {len(ensembl_ids):,} Ensembl IDs...")
    print("(This may take 1–2 minutes for the full transcriptome.)")
    mg = mygene.MyGeneInfo()
    result = mg.querymany(
        ensembl_ids,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        as_dataframe=True,
        verbose=False
    )
    symbol_map = result["symbol"].dropna().to_dict()
    print(f"  Mapped: {len(symbol_map):,} / {len(ensembl_ids):,} Ensembl IDs to gene symbols")

    # Save full mapping for reproducibility audit
    mapping_df = pd.DataFrame(
        list(symbol_map.items()), columns=["EnsemblID","Symbol"]
    )
    mapping_df.to_csv(MAPPING_CACHE, index=False)
    print(f"  Full mapping saved to: {MAPPING_CACHE}")
    return symbol_map

if os.path.exists(MAPPING_CACHE):
    print(f"\nLoading cached gene symbol mapping from {MAPPING_CACHE}")
    cached = pd.read_csv(MAPPING_CACHE).set_index("EnsemblID")["Symbol"].to_dict()
    # Check if cache covers current dataset
    current_ids = set(raw.index)
    cached_ids  = set(cached.keys())
    coverage = len(current_ids & cached_ids) / len(current_ids)
    if coverage > 0.9:
        symbol_map = cached
        print(f"  Cache coverage: {coverage:.1%} — using cached mapping.")
    else:
        print(f"  Cache coverage only {coverage:.1%} — re-querying.")
        symbol_map = build_symbol_map(raw.index.tolist())
else:
    symbol_map = build_symbol_map(raw.index.tolist())

# Apply mapping
raw.index = [symbol_map.get(g, g) for g in raw.index]

# Collapse duplicate symbols by summing read counts
raw = raw.groupby(raw.index).sum()
print(f"Genes after symbol mapping and deduplication: {raw.shape[0]:,}")

# Check VEN panel coverage
found  = [g for g in ALL_PANEL_GENES + INDIVIDUAL_MARKERS if g in raw.index]
missing = [g for g in ALL_PANEL_GENES + INDIVIDUAL_MARKERS if g not in raw.index]
print(f"\nVEN panel genes found in data: {len(found)}/{len(ALL_PANEL_GENES + INDIVIDUAL_MARKERS)}")
if found:   print(f"  Found:   {found}")
if missing: print(f"  Missing: {missing} (low/absent expression or mapping gap)")

# SECTION 5: NORMALISATION (CPM + log2)

cpm = raw.div(raw.sum(axis=0), axis=1) * 1_000_000
log2cpm = np.log2(cpm + 1)

print(f"\nNormalisation: CPM + log2(CPM+1)")
print(f"Library sizes (mean ± SD): "
      f"{raw.sum(axis=0).mean()/1e6:.1f}M ± {raw.sum(axis=0).std()/1e6:.1f}M reads")

# SECTION 6: SAMPLE SELECTION & LOG2FC

def select_samples(microglia_group, condition):
    mask = (
        (meta_df["OrgType"] == "Cortical") &
        (meta_df["Microglia"] == microglia_group) &
        (meta_df["Condition"] == condition)
    )
    return [c for c in meta_df[mask].index if c in log2cpm.columns]

# Report sample counts
print("\nCortical organoid sample counts:")
for mic in ["NoMicroglia", "WithMicroglia"]:
    for cond in ["Ground", "LEO"]:
        n = len(select_samples(mic, cond))
        print(f"  {mic:15s} {cond:6s}: n = {n}")

def compute_log2fc(leo_cols, gnd_cols):
    return log2cpm[leo_cols].mean(axis=1) - log2cpm[gnd_cols].mean(axis=1)

fc_groups = {
    "NoMicroglia": compute_log2fc(
        select_samples("NoMicroglia","LEO"),
        select_samples("NoMicroglia","Ground")
    ),
    "WithMicroglia": compute_log2fc(
        select_samples("WithMicroglia","LEO"),
        select_samples("WithMicroglia","Ground")
    ),
}
fc_groups["Combined"] = compute_log2fc(
    select_samples("NoMicroglia","LEO") + select_samples("WithMicroglia","LEO"),
    select_samples("NoMicroglia","Ground") + select_samples("WithMicroglia","Ground")
)

# Save all log2FC values for full transparency
for grp, fc in fc_groups.items():
    fc.rename("log2FC").to_csv(f"results/all_genes_log2FC_{grp}.csv", header=True)
print("\nAll-gene log2FC vectors saved to results/all_genes_log2FC_*.csv")

# SECTION 7: VEN PANEL ANALYSIS
# One-sample t-test + permutation test (N=10,000)
# Permutation tests the SPECIFICITY of the VEN panel signal:
#   Is the category mean more extreme than expected from a random gene set of the same size drawn from the full transcriptome?

N_PERM  = 10_000
RNG_SEED = 42
rng = np.random.default_rng(RNG_SEED)

def run_panel_analysis(fc_series, panel_dict, n_perm=N_PERM):
    """
    For each gene category in panel_dict:
      1. Identify genes present in the data
      2. One-sample t-test (H0: mean log2FC = 0)
      3. Permutation test: 10,000 random gene sets of equal size
         drawn from the full transcriptome background
    Returns DataFrame of results.
    """
    background_fc = fc_series.dropna().values
    rows = []

    for category, gene_list in panel_dict.items():
        present = [g for g in gene_list if g in fc_series.index
                   and not np.isnan(fc_series[g])]
        n_genes_panel    = len(gene_list)
        n_genes_present  = len(present)

        if n_genes_present == 0:
            rows.append({
                "Category": category, "n_panel": n_genes_panel,
                "n_present": 0, "mean_log2FC": np.nan,
                "t_stat": np.nan, "p_ttest": np.nan,
                "SD_above_null": np.nan, "perm_p": np.nan,
                "sig": "", "genes_present": "none found"
            })
            continue

        category_fc  = fc_series[present].values
        category_mean = category_fc.mean()

        # One-sample t-test
        if n_genes_present > 1:
            t_stat, p_ttest = stats.ttest_1samp(category_fc, popmean=0)
        else:
            t_stat, p_ttest = np.nan, np.nan

        # Permutation test
        null_means = np.array([
            rng.choice(background_fc, size=n_genes_present, replace=False).mean()
            for _ in range(n_perm)
        ])
        null_mean  = null_means.mean()
        null_sd = null_means.std()
        sd_above   = (category_mean - null_mean) / null_sd if null_sd > 0 else np.nan
        perm_p = np.mean(
            np.abs(null_means - null_mean) >= np.abs(category_mean - null_mean)
        )

        # Significance label
        if perm_p < 0.001: sig = "***"
        elif perm_p < 0.01: sig = "**"
        elif perm_p < 0.05: sig = "*"
        elif perm_p < 0.10: sig = "†"
        else:                  sig = ""

        rows.append({
            "Category": category,
            "n_panel": n_genes_panel,
            "n_present": n_genes_present,
            "mean_log2FC": round(category_mean, 4),
            "t_stat": round(t_stat, 3) if not np.isnan(t_stat) else np.nan,
            "p_ttest": round(p_ttest, 4) if not np.isnan(p_ttest) else np.nan,
            "SD_above_null":  round(sd_above, 3),
            "perm_p": round(perm_p, 4),
            "sig": sig,
            "genes_present":  ", ".join(present),
        })

    return pd.DataFrame(rows)

panel_results = {}
for group_name, fc in fc_groups.items():
    print(f"\n{'='*60}")
    print(f"GROUP: {group_name}")
    print(f"{'='*60}")
    df = run_panel_analysis(fc, VEN_PANEL)
    panel_results[group_name] = df
    print(df[["Category","n_present","mean_log2FC","p_ttest",
              "SD_above_null","perm_p","sig"]].to_string(index=False))
    df.to_csv(f"results/VEN_panel_{group_name}.csv", index=False)

print(f"\n{'='*60}")
print("INDIVIDUAL MARKER GENES")
print(f"{'='*60}")
individual_genes = ALL_PANEL_GENES + INDIVIDUAL_MARKERS
marker_rows = []
for group_name, fc in fc_groups.items():
    for gene in individual_genes:
        val = float(fc[gene]) if gene in fc.index else np.nan
        marker_rows.append({"Group": group_name, "Gene": gene, "log2FC": val})
marker_df = pd.DataFrame(marker_rows)

try:
    pivot = marker_df.pivot(index="Gene", columns="Group", values="log2FC").round(3)
    print(pivot.to_string())
    pivot.to_csv("results/individual_genes_pivot.csv")
except Exception as e:
    print(marker_df.to_string())
marker_df.to_csv("results/individual_genes_all.csv", index=False)

-print("PAPER-READY SUMMARY TABLE (Combined group)")
-summary = panel_results["Combined"][[
    "Category","n_present","mean_log2FC","t_stat","p_ttest",
    "SD_above_null","perm_p","sig","genes_present"
]].copy()
summary.columns = [
    "Category","N genes","Mean log2FC","t","p (t-test)",
    "SD above null","perm p","","Genes detected"
]
print(summary.to_string(index=False))
summary.to_csv("results/SUMMARY_TABLE_combined.csv", index=False)

metadata = {
    "dataset": "GSE259421",
    "paper": "Marotta et al. Stem Cells Transl Med 2024",
    "analysis_date": pd.Timestamp.now().strftime("%Y-%m-%d"),
    "n_genes_in_matrix": int(raw.shape[0]),
    "gene_mapping_source": "MyGene.info (mygene Python package, Xin et al. 2016)",
    "normalisation": "CPM + log2(CPM+1)",
    "permutation_n": N_PERM,
    "random_seed": RNG_SEED,
    "VEN_panel_source": "Keskin 2026 Fast Lane Hypothesis + VEN Fatigue Hypothesis",
    "sample_groups_analysed": {
        "NoMicroglia": {
            "LEO_n": len(select_samples("NoMicroglia","LEO")),
            "Ground_n": len(select_samples("NoMicroglia","Ground"))
        },
        "WithMicroglia": {
            "LEO_n": len(select_samples("WithMicroglia","LEO")),
            "Ground_n": len(select_samples("WithMicroglia","Ground"))
        }
    }
}
with open("results/analysis_metadata.json","w") as f:
    json.dump(metadata, f, indent=2)

CATEGORY_COLORS = {
    "Myelination": "#2166AC",
    "FastSignalling": "#4DAC26",
    "SocialCircuit": "#D01C8B",
    "LayerVProjection": "#F4A582",
    "MetabolicSupport": "#878787",
}

GROUPS_TO_PLOT = ["NoMicroglia", "WithMicroglia", "Combined"]
GROUP_LABELS = {
    "NoMicroglia": "Without Microglia",
    "WithMicroglia": "With Microglia",
    "Combined": "Combined",
}

def make_bar_figure(value_key, ylabel, title, filename, reference_lines=None):
    fig, axes = plt.subplots(1, 3, figsize=(17, 6), sharey=False)
    for ax, grp in zip(axes, GROUPS_TO_PLOT):
        df = panel_results[grp]
        categories = df["Category"].tolist()
        values = df[value_key].fillna(0).tolist()
        perm_ps = df["perm_p"].tolist()
        colors = [CATEGORY_COLORS[c] for c in categories]

        bars = ax.bar(range(len(categories)), values,
                      color=colors, edgecolor="black", linewidth=0.9, alpha=0.88)

        # Significance markers
        for i, (bar, pp) in enumerate(zip(bars, perm_ps)):
            if isinstance(pp, float) and pp < 0.10:
                sym = "***" if pp<0.001 else ("**" if pp<0.01 else ("*" if pp<0.05 else "†"))
                offset = max(abs(bar.get_height()), 0.05) * 0.15
                ypos   = bar.get_height() + offset
                ax.text(i, ypos, sym, ha="center", va="bottom",
                        fontsize=13, fontweight="bold", color="black")

        ax.axhline(0, color="black", linewidth=0.9)

        if reference_lines:
            for yval, style in reference_lines:
                ax.axhline(yval, color="grey", linewidth=0.7, linestyle=style)
                ax.axhline(-yval, color="grey", linewidth=0.7, linestyle=style)

        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels(
            [c.replace("Projection","Proj.").replace("Signalling","Sig.") for c in categories],
            rotation=35, ha="right", fontsize=9
        )
        ax.set_title(GROUP_LABELS[grp], fontsize=11, fontweight="bold", pad=8)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Category legend
    patches = [
        plt.Rectangle((0,0),1,1, fc=CATEGORY_COLORS[c], ec="black", lw=0.5)
        for c in CATEGORY_COLORS
    ]
    axes[2].legend(patches, list(CATEGORY_COLORS.keys()),
                   loc="best", fontsize=8, framealpha=0.8, title="VEN category")

    fig.suptitle(title, fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches="tight")
    plt.close()
    print(f" Saved: {filename}")

print("\nGenerating figures...")

make_bar_figure(
    value_key="SD_above_null",
    ylabel="SDs above genome-wide null\n(permutation test, N=10,000)",
    title=("VEN Gene Panel Specificity GSE259421 Human Cortical Organoids (ISS vs Ground)\n"
           "Keskin (2026) | VEN Fatigue Hypothesis\n"
           "Dotted lines: ±1.96 SD (p≈0.05 equivalent)"),
    filename="results/Fig1_permutation_specificity.png",
    reference_lines=[(1.96, ":")]
)

make_bar_figure(
    value_key="mean_log2FC",
    ylabel="Mean log₂FC (ISS vs Ground)",
    title=("VEN Gene Panel Expression Change — GSE259421 Human Cortical Organoids\n"
           "Keskin (2026) | VEN Fatigue Hypothesis\n"
           "* permutation p<0.05  † permutation p<0.10"),
    filename="results/Fig2_log2FC_raw.png",
)

def make_marker_heatmap():
    highlight = ["MBP","MAG","CNP","MOBP","SNAP25","NEFH","NEFL","BCL11B","FEZF2","NOS1"]
    plot_genes = [g for g in highlight if g in marker_df["Gene"].values]
    if not plot_genes:
        print("  No marker genes found — skipping heatmap.")
        return

    fig, ax = plt.subplots(figsize=(7, max(4, len(plot_genes)*0.55)))
    data_rows = []
    for gene in plot_genes:
        row = {}
        for grp in GROUPS_TO_PLOT:
            sub = marker_df[(marker_df.Gene==gene) & (marker_df.Group==grp)]
            row[GROUP_LABELS[grp]] = float(sub["log2FC"].values[0]) if len(sub) else np.nan
        data_rows.append(row)
    heat_df = pd.DataFrame(data_rows, index=plot_genes)

    vmax = max(abs(heat_df.values[~np.isnan(heat_df.values)]).max(), 0.5)
    im = ax.imshow(heat_df.values.astype(float), cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, aspect="auto")
    plt.colorbar(im, ax=ax, label="log₂FC (ISS vs Ground)")
    ax.set_xticks(range(len(heat_df.columns)))
    ax.set_xticklabels(heat_df.columns, rotation=20, ha="right", fontsize=10)
    ax.set_yticks(range(len(plot_genes)))
    ax.set_yticklabels(plot_genes, fontsize=9)
    ax.set_title("Key VEN Marker Gene Expression\nGSE259421 Cortical Organoids (ISS vs Ground)",
                 fontweight="bold", fontsize=10)

    for i, gene in enumerate(plot_genes):
        for j, col in enumerate(heat_df.columns):
            val = heat_df.loc[gene, col]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if abs(val) > vmax*0.55 else "black")

    plt.tight_layout()
    plt.savefig("results/Fig3_marker_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved: results/Fig3_marker_heatmap.png")

make_marker_heatmap()


print("ALL OUTPUT FILES")
for f in sorted(os.listdir("results")):
    size = os.path.getsize(f"results/{f}")
    print(f"  results/{f}  ({size:,} bytes)")

print("ANALYSIS COMPLETE")
print(f"Dataset : GSE259421 (Marotta et al. 2024, PMID 39441987)")
print(f"Organism : Homo sapiens (iPSC-derived cortical organoids)")
print(f"Comparison : ISS (38 days LEO) vs Ground control")
print(f"Gene map : MyGene.info — saved to results/ensembl_to_symbol.csv")
print(f"Panel : {len(ALL_PANEL_GENES)} genes across 5 VEN categories (pre-defined)")
print(f"Permutation: N={N_PERM:,}, seed={RNG_SEED}")
