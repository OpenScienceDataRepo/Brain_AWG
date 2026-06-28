"""
make_paper_figures.py
Generates all four publication-ready figures for:
  "The Social Cognition Paradox in Long-Duration Spaceflight:
   A VEN Fatigue Hypothesis for Duration-Dependent Emotion
   Recognition Decline"

OUTPUT
 
  figures/fig1_ert_paradox.png/.pdf -- Figure 1: ERT temporal profiles
  figures/fig2_domain_specificity.png/.pdf -- Figure 2: Domain specificity
  figures/fig3_molecular_dissociation.png/.pdf -- Figure 3: ISS vs OSD-202
  figures/fig4_organoid_permutation.png/.pdf   -- Figure 4: Organoid panel

REQUIRES
 
  data/raw/twins_cognitive_heatmap.csv
  data/raw/dev2024_raw_scores.csv
  results/gse239336_ven_signature.json
  results/osd202_ven_signature.json
  results/VEN_panel_Combined.csv
  results/VEN_panel_NoMicroglia.csv
  results/VEN_panel_WithMicroglia.csv

Run step1_create_cognitive_csvs.py and step2_run_analysis.py first.

Usage
   python make_paper_figures.py
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

os.makedirs("figures", exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth":  0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "figure.dpi":      150,
})

CAT_COLORS = {
    "Myelination": "#4E9AF1",
    "FastSignalling":  "#52C869",
    "SocialCircuit": "#E05252",
    "LayerVProjection": "#F5A623",
    "MetabolicSupport": "#9B59B6",
}
CAT_LABELS = {
    "Myelination": "Myelination",
    "FastSignalling":  "Fast\nSignalling",
    "SocialCircuit": "Social\nCircuit",
    "LayerVProjection": "Layer V\nProjection",
    "MetabolicSupport": "Metabolic\nSupport",
}

def sig_star(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    if p < 0.10: return "†"
    return ""

def add_sig(ax, x, y, label, color="black", fontsize=11, dy=0.08):
    if not label:
        return
    sign = 1 if y >= 0 else -1
    ax.text(x, y + sign * dy, label, ha="center", va="bottom" if y >= 0 else "top",
            fontsize=fontsize, color=color, fontweight="bold")

def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(f"figures/{stem}.{ext}", dpi=150, bbox_inches="tight")
    print(f"  Saved figures/{stem}.pdf/.png")
    plt.close(fig)

def fig1_ert_paradox():
    print("Generating Figure 1: ERT Paradox...")
    twins = pd.read_csv("data/raw/twins_cognitive_heatmap.csv")
    dev = pd.read_csv("data/raw/dev2024_raw_scores.csv")

    ert_tw = twins[(twins["task"] == "ERT") & (twins["metric"] == "speed")]
    phases_tw  = ["pre_flight", "inf_early_1_6", "inf_late_7_12", "post_flight"]
    labels_tw  = ["Pre-\nFlight", "Early\nIn-Flight\n(1-6 mo)", "Late\nIn-Flight\n(7-12 mo)", "Post-\nFlight"]
    vals_tw = [float(ert_tw[ert_tw["phase"] == p]["value"].values[0]) for p in phases_tw]

    ert_dev = dev[(dev["task"] == "ERT") & (dev["metric"] == "speed")].copy()
    pre_mean = float(ert_dev[ert_dev["phase"] == "pre_flight"]["mean"].values[0])
    pre_sd = float(ert_dev[ert_dev["phase"] == "pre_flight"]["sd"].values[0])
    ert_dev["z"] = -(ert_dev["mean"] - pre_mean) / pre_sd
    phases_dev = ["pre_flight", "early_flight", "late_flight", "early_post_flight", "late_post_flight"]
    labels_dev = ["Pre-\nFlight", "Early\nFlight", "Late\nFlight", "Early\nPost", "Late\nPost"]
    vals_dev = [float(ert_dev[ert_dev["phase"] == p]["z"].values[0]) for p in phases_dev]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    fig.subplots_adjust(wspace=0.12)

    # Panel A: 340-day mission
    ax = axes[0]
    phase_colors_tw = ["#888888", "#3A7FC1", "#E05252", "#9B59B6"]
    ax.bar(range(4), vals_tw, color=phase_colors_tw, alpha=0.88,
           edgecolor="black", linewidth=0.7, zorder=2)
    ax.axhline(0, color="black", lw=0.9, zorder=3)
    ax.axhline(-1.0, color="#E05252", ls="--", alpha=0.55, lw=1.4,
               label=">1 SD decline threshold")
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels_tw, fontsize=9)
    ax.set_ylabel("ERT Speed Change\n(SD units, spaceflight minus ground)", fontsize=10)
    ax.set_ylim(-3.2, 2.0)
    ax.set_title("A  340-day Mission (NASA Twins Study, $N=1$)", fontsize=10, fontweight="bold", loc="left")
    ax.annotate("$-1.8$ SD\n(late in-flight)",
                xy=(2, vals_tw[2]),
                xytext=(2.6, vals_tw[2] + 1.2),
                fontsize=9, color="#E05252",
                arrowprops=dict(arrowstyle="->", color="#E05252", lw=1.4))
    ax.legend(fontsize=8, frameon=False)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, alpha=0.3, lw=0.6)

    # Panel B: 6-month missions
    ax = axes[1]
    ax.bar(range(5), vals_dev, color="#52C869", alpha=0.88,
           edgecolor="black", linewidth=0.7, zorder=2)
    ax.axhline(0, color="black", lw=0.9, zorder=3)
    ax.axhline(-1.0, color="#E05252", ls="--", alpha=0.55, lw=1.4)
    ax.set_xticks(range(5))
    ax.set_xticklabels(labels_dev, fontsize=9)
    ax.set_title("B  6-month Missions (Dev et al. 2024, $N=24$)", fontsize=10, fontweight="bold", loc="left")
    ax.text(2, 1.55, "STABLE", fontsize=14, color="#52C869",
            fontweight="bold", ha="center")
    ax.annotate("$+0.106$ SD", xy=(2, vals_dev[2]),
                xytext=(3.2, vals_dev[2] + 0.35),
                fontsize=9, color="#52C869",
                arrowprops=dict(arrowstyle="->", color="#52C869", lw=1.4))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, alpha=0.3, lw=0.6)

    fig.suptitle("Figure 1  The ERT Paradox: Social Cognition Declines Selectively in Long-Duration Spaceflight",
                 fontsize=10.5, fontweight="bold", y=1.01)
    save(fig, "fig1_ert_paradox")



def fig2_domain_specificity():
    print("Generating Figure 2: Domain Specificity...")
    twins = pd.read_csv("data/raw/twins_cognitive_heatmap.csv")

    TASK_DOMAINS = {
        "ERT": "Social", "VOLT": "Memory", "F2B": "Working Mem.",
        "AM": "Reasoning", "LOT": "Spatial", "MRT": "Reasoning",
        "DSST": "Processing", "PVT": "Attention", "BART": "Risk",
        "MP": "Motor",
    }
    DOMAIN_COLOR = {
        "Social": "#E05252", "Memory": "#F5A623", "Working Mem.": "#F5A623",
        "Spatial": "#3A7FC1", "Reasoning": "#52C869", "Processing": "#9B59B6",
        "Attention": "#1ABC9C", "Risk": "#D68910", "Motor": "#717D7E",
    }

    speed = twins[(twins["metric"] == "speed") & twins["task"].isin(TASK_DOMAINS)].copy()
    speed["domain"] = speed["task"].map(TASK_DOMAINS)

    early = speed[speed["phase"] == "inf_early_1_6"].set_index("task")["value"]
    late  = speed[speed["phase"] == "inf_late_7_12"].set_index("task")["value"]

    tasks = sorted(
        [t for t in early.index if t in late.index],
        key=lambda t: late[t] - early[t]
    )
    early_v = [float(early[t]) for t in tasks]
    late_v  = [float(late[t])  for t in tasks]
    domains = [TASK_DOMAINS[t] for t in tasks]
    colors  = [DOMAIN_COLOR.get(d, "#888") for d in domains]

    x = np.arange(len(tasks))
    w = 0.38

    fig, ax = plt.subplots(figsize=(12, 5.5))

    bars_e = ax.bar(x - w/2, early_v, w, color=colors, alpha=0.60,
                    edgecolor="black", linewidth=0.6, label="Early In-Flight (months 1-6)")
    bars_l = ax.bar(x + w/2, late_v,  w, color=colors, alpha=0.92,
                    edgecolor="black", linewidth=0.6, label="Late In-Flight (months 7-12)")

    ax.axhline(0, color="black", lw=0.9)
    ax.axhline(-1.0, color="red", ls="--", alpha=0.35, lw=1.2)

    xtick_labels = [f"{t}\n({TASK_DOMAINS[t]})" for t in tasks]
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels, fontsize=9, rotation=10, ha="right")
    ax.set_ylabel("Performance (SD units, spaceflight minus ground)", fontsize=10)
    ax.set_title("Figure 2  Domain-Specific Late In-Flight Decline in the 340-day Mission\n"
                 "ERT (Social) shows the largest early-to-late change of any cognitive domain",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, alpha=0.3, lw=0.6)

    ert_pos = tasks.index("ERT")
    ax.annotate(
        "ERT (Social)\nEarly: "
        f"{early_v[ert_pos]:+.1f} SD\nLate: {late_v[ert_pos]:+.1f} SD",
        xy=(ert_pos + w/2, late_v[ert_pos]),
        xytext=(ert_pos + w/2 + 3.5, 0.8),
        fontsize=9, color="#E05252", fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#E05252", lw=1.4),
    )
    save(fig, "fig4_domain_specificity")


 # FIGURE 3  Molecular dissociation: ISS frontal cortex vs ground-based analogue
def fig3_molecular_dissociation():
    print("Generating Figure 3: Molecular Dissociation...")
    with open("results/gse239336_ven_signature.json") as f:
        gse = json.load(f)
    with open("results/osd202_ven_signature.json") as f:
        osd = json.load(f)

    cat_keys  = ["myelination", "fast_signalling", "social_circuit", "layer5_proj", "metabolic"]
    cat_names = ["Myelination", "Fast\nSignalling", "Social\nCircuit", "Layer V\nProjection", "Metabolic\nSupport"]

    iss_fc  = [gse[c]["mean_logfc"] for c in cat_keys]
    gnd_fc  = [osd[c]["mean_logfc"] if c in osd else np.nan for c in cat_keys]
    iss_p = [gse[c]["p_value"] for c in cat_keys]
    iss_pp  = [gse[c]["p_permutation"] for c in cat_keys]
    gnd_p = [osd[c].get("p_value", np.nan) for c in cat_keys]
    gnd_pp  = [osd[c].get("p_permutation", np.nan) if c in osd else np.nan for c in cat_keys]
    iss_sd = [gse[c]["sd_above_null"] for c in cat_keys]
    gnd_sd  = [osd[c].get("sd_above_null", np.nan) if c in osd else np.nan for c in cat_keys]

    x = np.arange(len(cat_keys))
    w = 0.38

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # Left: mean log2FC bars
    ax = axes[0]
    b1 = ax.bar(x - w/2, iss_fc, w, color="#3A7FC1", alpha=0.88,
                edgecolor="black", lw=0.7, label="ISS Frontal Cortex (GSE239336)")
    b2 = ax.bar(x + w/2, gnd_fc, w, color="#E05252", alpha=0.88,
                edgecolor="black", lw=0.7, label="Ground Stress (OSD-202)")
    ax.axhline(0, color="black", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(cat_names, fontsize=9)
    ax.set_ylabel("Mean $\\log_2$FC (ISS minus ground)", fontsize=10)
    ax.set_title("A  Mean $\\log_2$FC per VEN category", fontsize=10, fontweight="bold", loc="left")
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, alpha=0.3, lw=0.6)
    for i, (iv, gv, ip, gp) in enumerate(zip(iss_fc, gnd_fc, iss_p, gnd_p)):
        add_sig(ax, i - w/2, iv, sig_star(ip) if not np.isnan(ip) else "", color="#3A7FC1")
        add_sig(ax, i + w/2, gv, sig_star(gp) if not np.isnan(gp) else "", color="#E05252")
    ax.annotate("Myelination:\nISS +0.381*\nGround -0.254*",
                xy=(0, 0), xytext=(1.2, 0.62),
                fontsize=8.5, color="#333",
                arrowprops=dict(arrowstyle="->", color="#555", lw=1.2))

    # Right: SD above null
    ax = axes[1]
    valid = [(n, is_, gs, ip, gp)
             for n, is_, gs, ip, gp in zip(cat_names, iss_sd, gnd_sd, iss_pp, gnd_pp)
             if not np.isnan(gs)]
    if valid:
        ns = [v[0] for v in valid]
        iss_s = [v[1] for v in valid]
        gnd_s = [v[2] for v in valid]
        ipp = [v[3] for v in valid]
        gpp = [v[4] for v in valid]
        xv = np.arange(len(ns))
        ax.bar(xv - w/2, iss_s, w, color="#3A7FC1", alpha=0.88, edgecolor="black", lw=0.7,
               label="ISS")
        ax.bar(xv + w/2, gnd_s, w, color="#E05252", alpha=0.88, edgecolor="black", lw=0.7,
               label="Ground stress")
        ax.axhline(0, color="black", lw=0.9)
        ax.axhline(1.96, color="gray", ls=":", lw=1, alpha=0.7, label="p=0.05 threshold")
        ax.axhline(-1.96, color="gray", ls=":", lw=1, alpha=0.7)
        ax.set_xticks(xv)
        ax.set_xticklabels(ns, fontsize=9)
        ax.set_ylabel("SD above genome-wide permutation null", fontsize=10)
        ax.set_title("B  Pathway specificity (permutation test)", fontsize=10, fontweight="bold", loc="left")
        ax.legend(fontsize=8, frameon=False)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, alpha=0.3, lw=0.6)
        for i, (is_, gp) in enumerate(zip(iss_s, gpp)):
            add_sig(ax, xv[i] - w/2, is_, sig_star(ipp[i]) if not np.isnan(ipp[i]) else "", color="#3A7FC1")
            add_sig(ax, xv[i] + w/2, gnd_s[i], sig_star(gp) if not np.isnan(gp) else "", color="#E05252")
    else:
        # Fallback: just show ISS SD
        ax.bar(x, iss_sd, color="#3A7FC1", alpha=0.88, edgecolor="black", lw=0.7, label="ISS")
        ax.axhline(0, color="black", lw=0.9)
        ax.set_xticks(x); ax.set_xticklabels(cat_names, fontsize=9)
        ax.set_ylabel("SD above genome-wide permutation null", fontsize=10)
        ax.set_title("B  ISS Pathway specificity", fontsize=10, fontweight="bold", loc="left")

    fig.suptitle("Figure 3  Myelination Upregulation Is Spaceflight-Specific in Mouse Frontal Cortex\n"
                 "ISS tissue shows a targeted myelination response absent in ground-based stress",
                 fontsize=10.5, fontweight="bold", y=1.02)
    save(fig, "fig2_molecular_dissociation")


# FIGURE 4  Organoid permutation specificity: VEN panel in GSE259421
def fig4_organoid_permutation():
    print("Generating Figure 4: Organoid Permutation Specificity...")

    # Load all three groups
    dfs = {}
    for group, fname in [
        ("Combined ($n=9$ LEO, $n=9$ Gnd)", "results/VEN_panel_Combined.csv"),
        ("Without Microglia ($n=4$ LEO, $n=5$ Gnd)", "results/VEN_panel_NoMicroglia.csv"),
        ("With Microglia ($n=5$ LEO, $n=5$ Gnd)",   "results/VEN_panel_WithMicroglia.csv"),
    ]:
        if os.path.exists(fname):
            df = pd.read_csv(fname)
            # Normalise column names
            df.columns = [c.strip() for c in df.columns]
            dfs[group] = df

    if not dfs:
        print("  ERROR: No VEN_panel CSV files found in results/. Skipping Figure 4.")
        return

    # Use canonical category order
    cat_order = ["Myelination", "FastSignalling", "SocialCircuit",
                 "LayerVProjection", "MetabolicSupport"]
    cat_nice = ["Myelination", "Fast\nSignalling", "Social\nCircuit",
                 "Layer V\nProjection", "Metabolic\nSupport"]
    cat_col = [CAT_COLORS[c] for c in cat_order]

    group_names = list(dfs.keys())
    n_groups = len(group_names)
    x = np.arange(len(cat_order))
    bar_w = 0.22
    offsets = np.linspace(-(n_groups-1)/2, (n_groups-1)/2, n_groups) * bar_w
    group_hatches = ["", "//", ".."]
    group_alphas  = [0.9, 0.65, 0.65]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.subplots_adjust(wspace=0.28)

    #   Left panel: SD above null for all three groups  
    ax = axes[0]
    for gi, (gname, offset, hatch, alpha) in enumerate(zip(group_names, offsets, group_hatches, group_alphas)):
        df = dfs[gname]
        df_idx = df.set_index("Category")
        sd_vals = []
        pp_vals = []
        for cat in cat_order:
            if cat in df_idx.index:
                row = df_idx.loc[cat]
                sd_vals.append(float(row["SD_above_null"]))
                pp_vals.append(float(row["perm_p"]))
            else:
                sd_vals.append(np.nan)
                pp_vals.append(np.nan)

        bars = ax.bar(x + offset, sd_vals, bar_w,
                      color=cat_col, alpha=alpha,
                      edgecolor="black", lw=0.6, hatch=hatch,
                      label=gname.split("(")[0].strip())
        for xi, (sv, pp) in enumerate(zip(sd_vals, pp_vals)):
            if not np.isnan(sv) and not np.isnan(pp):
                s = sig_star(pp)
                if s:
                    ax.text(xi + offset, sv + (0.18 if sv >= 0 else -0.18),
                            s, ha="center",
                            va="bottom" if sv >= 0 else "top",
                            fontsize=10, fontweight="bold",
                            color="black")

    ax.axhline(0, color="black", lw=0.9)
    ax.axhline(1.96, color="gray", ls=":", lw=1, alpha=0.65, label="p=0.05 (permutation)")
    ax.axhline(-1.96, color="gray", ls=":", lw=1, alpha=0.65)
    ax.set_xticks(x)
    ax.set_xticklabels(cat_nice, fontsize=9)
    ax.set_ylabel("SD above genome-wide permutation null\n($N=10{,}000$ random gene sets, seed 42)", fontsize=9)
    ax.set_title("A  Permutation specificity by group", fontsize=10, fontweight="bold", loc="left")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, alpha=0.3, lw=0.6)

    ax = axes[1]
    df_c  = dfs[group_names[0]].set_index("Category")
    fc_vals = []
    pp_vals = []
    for cat in cat_order:
        if cat in df_c.index:
            row = df_c.loc[cat]
            fc_vals.append(float(row["mean_log2FC"]))
            pp_vals.append(float(row["perm_p"]))
        else:
            fc_vals.append(np.nan)
            pp_vals.append(np.nan)

    bars = ax.bar(x, fc_vals, 0.55, color=cat_col, alpha=0.88,
                  edgecolor="black", lw=0.7)
    ax.axhline(0, color="black", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(cat_nice, fontsize=9)
    ax.set_ylabel("Mean $\\log_2$FC (ISS minus ground)", fontsize=10)
    ax.set_title("B  Mean log$_2$FC (combined group)", fontsize=10, fontweight="bold", loc="left")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, alpha=0.3, lw=0.6)
    for xi, (fv, pp) in enumerate(zip(fc_vals, pp_vals)):
        if not np.isnan(fv) and not np.isnan(pp):
            s = sig_star(pp)
            if s:
                ax.text(xi, fv + (0.06 if fv >= 0 else -0.06), s,
                        ha="center",
                        va="bottom" if fv >= 0 else "top",
                        fontsize=11, fontweight="bold", color="black")

    # Annotations for key findings
    lv_idx = cat_order.index("LayerVProjection")
    my_idx = cat_order.index("Myelination")
    ax.annotate("Layer V:\n+1.347** (6.17 SD)",
                xy=(lv_idx, fc_vals[lv_idx]),
                xytext=(lv_idx + 0.65, fc_vals[lv_idx] - 0.45),
                fontsize=8, color=CAT_COLORS["LayerVProjection"], fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=CAT_COLORS["LayerVProjection"], lw=1.2))
    ax.annotate("Myelination:\n-0.456* (-2.35 SD)\n(predicted direction)",
                xy=(my_idx, fc_vals[my_idx]),
                xytext=(my_idx + 0.6, 0.28),
                fontsize=8, color=CAT_COLORS["Myelination"], fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=CAT_COLORS["Myelination"], lw=1.2))

    fig.suptitle("Figure 4  VEN Gene Panel in Human iPSC-Derived Cortical Organoids (GSE259421, ISS vs Ground)\n"
                 "Layer V Projection genes 6.17 SD above null; Myelination downregulation matches a priori prediction",
                 fontsize=10, fontweight="bold", y=1.02)
    save(fig, "fig3_organoid_permutation")


if __name__ == "__main__":
    print("=" * 65)
    print("  Generating publication figures for VEN Fatigue Hypothesis")
    print("=" * 65 + "\n")

    missing = []
    required = [
        "data/raw/twins_cognitive_heatmap.csv",
        "data/raw/dev2024_raw_scores.csv",
        "results/gse239336_ven_signature.json",
        "results/osd202_ven_signature.json",
        "results/VEN_panel_Combined.csv",
    ]
    for f in required:
        if not os.path.exists(f):
            missing.append(f)
    if missing:
        print("MISSING REQUIRED FILES:")
        for f in missing:
            print(f"  {f}")
        print("\nRun step1_create_cognitive_csvs.py and step2_run_analysis.py first.")
        import sys; sys.exit(1)

    fig1_ert_paradox()
    fig2_domain_specificity()
    fig3_molecular_dissociation()
    fig4_organoid_permutation()

    print("\nAll figures generated in figures/")
    print("Files: fig1_ert_paradox, fig2_domain_specificity,")
    print(" fig3_molecular_dissociation, fig4_organoid_permutation")
    print("Each saved as .pdf (vector) and .png (raster, 150 dpi)")
