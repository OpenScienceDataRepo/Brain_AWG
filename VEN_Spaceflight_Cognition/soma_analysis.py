"""
soma_analysis.py
Summarises SOMA atlas (soma.weill.cornell.edu) findings for VEN panel genes
queried manually from NASA Twin Study RNAseq and I4 cell-free/PBMC RNA data.

Data extracted from SOMA Browser PDFs (May 2026).
Genes queried: MBP, MOG, NEFL, NEFM, NEFH, SYP, VDAC1, SNAP25
"""

import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

os.makedirs("results", exist_ok=True)

#  Manually extracted significant / best results from SOMA  
# Each entry = best (lowest p) in-flight spaceflight comparison per gene.
# Dataset: "Twin" = NASA Twins Study PolyA+ blood RNA; "I4" = Inspiration4
# CellType: CPT = peripheral blood mononuclear cells (CPT tube); cfRNA = cell-free RNA
# Comparison: the label from SOMA

SOMA_RESULTS = [
    # Gene, Category, Dataset, CellType, Comparison, log2FC, pvalue, qvalue, note
    ("MBP", "Myelination", "I4 PBMC",  "PBMC",  "FP2: post vs pre-flight", 0.602,  0.0, 0.0, "I4 3-day mission; p machine-zero"),
    ("MOG", "Myelination", "Twin", "CPT", "In-flight 2nd Half vs Ground",  5.821,  4.59e-27, 3.04e-24,  "Year-long ISS; year-long CPT blood"),
    ("NEFL",  "FastSignalling", "Twin", "CPT", "In-flight 2nd Half vs Ground",  3.133,  6.42e-5, 3.99e-4, "Neurofilament-L in blood"),
    ("NEFM", "FastSignalling", "Twin", "CPT", "In-flight 2nd Half vs Ground",  1.033,  0.0427, 0.0819, "Trend; q<0.10"),
    ("NEFH", "FastSignalling", "I4 cfRNA", "cfRNA", "Post-flight vs Pre-flight", 0.998,  0.0361, 0.136, "I4 cell-free RNA; q trend"),
    ("VDAC1", "MetabolicSupport", "Twin", "CPT", "In-flight 2nd Half vs Ground",  -0.583, 1.32e-4, 7.35e-4, "Downregulated; mitochondrial"),
    ("SYP",   "MetabolicSupport", "Twin", "CPT", "In-flight 1st Half vs Ground",  0.732,  0.00129, 0.0194, "Synaptophysin upregulated"),
    ("SNAP25","MetabolicSupport", "Twin", "CPT", "In-flight 2nd Half vs Ground",  3.331,  0.00757, 0.0211, "Synaptic protein upregulated"),
]

df = pd.DataFrame(SOMA_RESULTS,
    columns=["Gene","Category","Dataset","CellType","Comparison",
             "log2FC","pvalue","qvalue","Note"])
df.to_csv("results/soma_VEN_panel_results.csv", index=False)
print("SOMA results table:")
print(df[["Gene","Category","Dataset","log2FC","pvalue","qvalue"]].to_string(index=False))

def sig_label(p, q):
    if p == 0 or p < 1e-10: return "***"
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    if q < 0.20: return "+"
    return "ns"

df["sig"] = df.apply(lambda r: sig_label(r["pvalue"], r["qvalue"]), axis=1)

CAT_COLORS = {
    "Myelination": "#4472C4",
    "FastSignalling": "#70AD47",
    "MetabolicSupport": "#7030A0",
}
DATASET_HATCH = {"I4 PBMC": "//", "I4 cfRNA": "\\\\", "Twin": ""}

fig, ax = plt.subplots(figsize=(11, 5))

x = np.arange(len(df))
bars = []
for i, row in df.iterrows():
    b = ax.bar(i, row["log2FC"],
               color=CAT_COLORS.get(row["Category"], "#888888"),
               hatch=DATASET_HATCH.get(row["Dataset"], ""),
               edgecolor="black", lw=0.8, alpha=0.88, width=0.65)
    bars.append(b)

ax.axhline(0, color="black", lw=0.9)

# Annotate significance
ymax = df["log2FC"].abs().max()
for i, row in df.iterrows():
    if row["sig"] not in ("", "ns"):
        ypos = row["log2FC"] + (0.15 if row["log2FC"] >= 0 else -0.45)
        ax.text(i, ypos, row["sig"], ha="center", va="bottom", fontsize=9, fontweight="bold")
    # Dataset tag
    tag = row["Dataset"].replace("I4 PBMC", "I4").replace("I4 cfRNA", "cfRNA").replace("Twin", "Twin")
    ax.text(i, -0.7, tag, ha="center", va="top", fontsize=7, color="#444444", style="italic")

ax.set_xticks(x)
ax.set_xticklabels(df["Gene"], fontsize=10, fontweight="bold")
ax.set_ylabel("log₂ Fold Change (spaceflight vs ground/pre-flight)", fontsize=9)
ax.set_title(
    "SOMA Atlas VEN Panel Genes in Human Astronaut Blood\n"
    "NASA Twin Study RNAseq (CPT cells, year-long ISS) & Inspiration4 PBMC/cfRNA\n"
    "soma.weill.cornell.edu · Overbey et al. Nature 2024 · Meydan & Mason Lab",
    fontsize=10, fontweight="bold")

# Legend: categories
from matplotlib.patches import Patch
legend_patches = [Patch(facecolor=c, edgecolor="black", label=k)
                  for k, c in CAT_COLORS.items()]
legend_patches += [
    Patch(facecolor="white", edgecolor="black", hatch="//",  label="I4 PBMC (3-day)"),
    Patch(facecolor="white", edgecolor="black", hatch="\\\\", label="I4 cfRNA (3-day)"),
    Patch(facecolor="white", edgecolor="black", label="NASA Twins (year-long ISS)"),
]
ax.legend(handles=legend_patches, fontsize=8, loc="upper left", framealpha=0.85)

ax.yaxis.grid(True, alpha=0.25, lw=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim(-1.2, 7.2)

plt.tight_layout()
plt.savefig("results/soma_VEN_panel_blood.png", dpi=200, bbox_inches="tight")
plt.close()
print("\nSaved: results/soma_VEN_panel_blood.png")

n_sig = (df["pvalue"] < 0.05).sum()
n_trend = ((df["pvalue"] >= 0.05) & (df["qvalue"] < 0.20)).sum()
n_total = len(df)
print(f"\nOut of {n_total} genes queried:")
print(f" Significant (p<0.05): {n_sig}")
print(f" Trend (q<0.20): {n_trend}")
print(f" Not significant: {n_total - n_sig - n_trend}")
print("\nTop hits:")
for _, r in df[df["pvalue"] < 0.05].sort_values("pvalue").iterrows():
    print(f" {r.Gene:8s} {r.Category:20s} log2FC={r.log2FC:+.3f}  p={r.pvalue:.2e}  {r.sig}  [{r.Dataset}]")
print("\nDone.")
