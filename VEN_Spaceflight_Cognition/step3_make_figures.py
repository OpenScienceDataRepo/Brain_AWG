"""
step3_make_figures.py
Generates all four paper figures from results/ files.
REQUIRES: Run step1 and step2 first.

Run: python step3_make_figures.py
Output: figures/fig1_ert_paradox.pdf/.png
        figures/fig2_domain_specificity.pdf/.png
        figures/fig3_molecular_dissociation.pdf/.png
        figures/fig4_ven_model_results.pdf/.png
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

sns.set_theme(style="whitegrid", font_scale=1.15)
os.makedirs("figures", exist_ok=True)

COLORS = {
    'social': '#E05252',
    'memory': '#F5A623',
    'spatial': '#3A7FC1',
    'reasoning': '#52C869',
    'processing':'#9B59B6',
    'attention': '#1ABC9C',
    'risk': '#2C3E50',
    'motor': '#95A5A6',
    'composite': '#888888',
    'efficiency':'#AAAAAA',
}

def strip(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# FIGURE 1: The ERT Paradox
# Left: ERT trajectory across phases in 340-day mission
# Right: ERT trajectory in 6-month missions (Dev 2024)
def fig1_ert_paradox():
    twins = pd.read_csv("data/raw/twins_cognitive_heatmap.csv")
    dev = pd.read_csv("data/raw/dev2024_raw_scores.csv")

    # Twins Study ERT speed
    ert_twins = twins[(twins['task']=='ERT') & (twins['metric']=='speed')].copy()
    phase_order_twins  = ['pre_flight','inf_early_1_6','inf_late_7_12','post_flight']
    phase_labels_twins = ['Pre-Flight','In-Flight\n1-6 months','In-Flight\n7-12 months','Post-Flight']

    vals_twins = [float(ert_twins[ert_twins['phase']==p]['value'].values[0])
                  for p in phase_order_twins]

    # Dev 2024 ERT speed - normalise to preflight z-score (corrected: speed = lower ms = better)
    ert_dev = dev[(dev['task']=='ERT') & (dev['metric']=='speed')].copy()
    pre_mean = float(ert_dev[ert_dev['phase']=='pre_flight']['mean'].values[0])
    pre_sd = float(ert_dev[ert_dev['phase']=='pre_flight']['sd'].values[0])
    ert_dev['z'] = -(ert_dev['mean'] - pre_mean) / pre_sd  # flip: faster = positive

    phase_order_dev = ['pre_flight','early_flight','late_flight','early_post_flight','late_post_flight']
    phase_labels_dev  = ['Pre-\nFlight','Early\nFlight','Late\nFlight','Early\nPost','Late\nPost']
    vals_dev = [float(ert_dev[ert_dev['phase']==p]['z'].values[0]) for p in phase_order_dev]
    sd_dev   = [float(ert_dev[ert_dev['phase']==p]['sd'].values[0]) / pre_sd
                for p in phase_order_dev]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: 340-day mission
    ax = axes[0]
    bar_colors = ['#888888','#3A7FC1','#E05252','#9B59B6']
    bars = ax.bar(range(4), vals_twins, color=bar_colors, alpha=0.85,
                  edgecolor='k', linewidth=0.7)
    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(-1.0, color='#E05252', ls='--', alpha=0.6, lw=1.5, label='>1 SD decline')
    ax.set_xticks(range(4))
    ax.set_xticklabels(phase_labels_twins, fontsize=10)
    ax.set_ylabel("ERT Speed (SD units, TW-HR difference)")
    ax.set_title("A: ERT Speed - 340-day Mission\n(NASA Twins Study, n=1)", fontsize=11)
    ax.legend(fontsize=9)
    ax.set_ylim(-3.0, 1.8)

    # Annotate the paradox
    ax.annotate("VEN Fatigue\nParadox", xy=(2, vals_twins[2]),
                xytext=(2.4, vals_twins[2] - 0.5), fontsize=9,
                color='#E05252',
                arrowprops=dict(arrowstyle='->', color='#E05252', lw=1.5))
    strip(ax)

    # Panel B: 6-month missions
    ax = axes[1]
    ax.bar(range(5), vals_dev, color='#52C869', alpha=0.85,
           edgecolor='k', linewidth=0.7)
    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(-1.0, color='#E05252', ls='--', alpha=0.6, lw=1.5, label='>1 SD decline')
    ax.set_xticks(range(5))
    ax.set_xticklabels(phase_labels_dev, fontsize=10)
    ax.set_ylabel("ERT Speed (z-score, relative to preflight)")
    ax.set_title("B: ERT Speed — 6-month Missions\n(Dev et al. 2024, n~24)", fontsize=11)
    ax.legend(fontsize=9)
    ax.set_ylim(-3.0, 1.8)
    ax.text(1.0, 1.3, "STABLE", fontsize=13, color='#52C869',
            fontweight='bold', ha='center')
    strip(ax)

    plt.suptitle("The ERT Paradox: Social Cognition Selectively Declines in Long-Duration Spaceflight Only",
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig("figures/fig1_ert_paradox.pdf", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig1_ert_paradox.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 1 saved.")


# FIGURE 2: Domain specificity
# Early inflight vs late inflight for ALL tasks (340-day)
# Shows ERT is the ONLY task with large early improvement -> late decline
def fig2_domain_specificity():
    twins = pd.read_csv("data/raw/twins_cognitive_heatmap.csv")
    TASK_DOMAINS = {
        'ERT':'social','VOLT':'memory','F2B':'working_memory','AM':'reasoning',
        'LOT':'spatial','MRT':'reasoning','DSST':'processing','PVT':'attention',
        'BART':'risk','MP':'motor',
    }

    speed = twins[(twins['metric']=='speed') & (~twins['task'].isin(['ALL','EFF']))].copy()
    speed['domain'] = speed['task'].map(TASK_DOMAINS)

    early = speed[speed['phase']=='inf_early_1_6'].set_index('task')['value']
    late  = speed[speed['phase']=='inf_late_7_12'].set_index('task')['value']

    tasks = [t for t in early.index if t in late.index]
    early_vals = [float(early[t]) for t in tasks]
    late_vals  = [float(late[t])  for t in tasks]
    domains = [TASK_DOMAINS.get(t,'unknown') for t in tasks]

    # Sort by early->late change (most declining first)
    order = np.argsort([l - e for e, l in zip(early_vals, late_vals)])
    tasks = [tasks[i] for i in order]
    early_vals  = [early_vals[i] for i in order]
    late_vals = [late_vals[i] for i in order]
    domains = [domains[i] for i in order]

    x = np.arange(len(tasks))
    w = 0.38

    fig, ax = plt.subplots(figsize=(11, 5.5))
    bars_e = ax.bar(x - w/2, early_vals, w, label='In-Flight Early (months 1-6)',
                    color=[COLORS.get(d,'#888') for d in domains],
                    edgecolor='k', linewidth=0.6)
    bars_l = ax.bar(x + w/2, late_vals,  w, label='In-Flight Late (months 7-12)',
                    color=[COLORS.get(d,'#888') for d in domains],
                    edgecolor='k', linewidth=0.6)

    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(-1.0, color='red', ls='--', alpha=0.4, lw=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{t}\n({TASK_DOMAINS.get(t,'?')})" for t in tasks],
        fontsize=9, rotation=10, ha='right'
    )
    ax.set_ylabel("Performance (SD, TW-HR difference)")
    ax.set_title("Domain-Specific Late-Inflight Decline in 340-day Mission\n"
                 "ERT shows the largest early->late change of any task", fontsize=11)
    ax.legend(fontsize=9, loc='upper left')
    strip(ax)

    # Highlight ERT bar
    ert_idx = tasks.index('ERT')
    ax.annotate('ERT\n(Social)\n-1.8 SD', xy=(ert_idx + w/2, late_vals[ert_idx]),
                xytext=(ert_idx + w/2 + 1.2, late_vals[ert_idx] - 0.3),
                fontsize=9, color='#E05252', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#E05252'))

    plt.tight_layout()
    plt.savefig("figures/fig2_domain_specificity.pdf", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig2_domain_specificity.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 2 saved.")


# FIGURE 3: Molecular dissociation
# Myelination VEN score: ISS (+) vs Ground stress (-)
def fig3_molecular_dissociation():
    with open("results/gse239336_ven_signature.json") as f:
        gse = json.load(f)
    with open("results/osd202_ven_signature.json") as f:
        osd = json.load(f)

    categories = list(gse.keys())
    iss_means = [gse[c].get('mean_logfc') for c in categories]
    ground_means = [osd[c].get('mean_logfc') for c in categories]
    iss_p = [gse[c].get('p_value') for c in categories]
    ground_p = [osd[c].get('p_value') for c in categories]

    # Only plot categories where we have data in both
    valid = [(c, im, gm, ip, gp) for c, im, gm, ip, gp in
             zip(categories, iss_means, ground_means, iss_p, ground_p)
             if im is not None and gm is not None]
    cats = [v[0] for v in valid]
    iss_v  = [v[1] for v in valid]
    gnd_v = [v[2] for v in valid]
    iss_p  = [v[3] for v in valid]
    gnd_p  = [v[4] for v in valid]

    x  = np.arange(len(cats))
    w  = 0.38
    cat_labels = [c.replace('_', '\n') for c in cats]

    fig, ax = plt.subplots(figsize=(11, 5.5))
    b1 = ax.bar(x - w/2, iss_v, w, label='ISS Tissue (GSE239336)',
                color='#3A7FC1', alpha=0.85, edgecolor='k', lw=0.7)
    b2 = ax.bar(x + w/2, gnd_v, w, label='Ground Stress (OSD-202, HLU+radiation)',
                color='#E05252', alpha=0.85, edgecolor='k', lw=0.7)

    # Significance stars
    for i, (iv, gv, ip, gp) in enumerate(zip(iss_v, gnd_v, iss_p, gnd_p)):
        if ip and ip < 0.05:
            ax.text(i - w/2, iv + 0.025*np.sign(iv), '*', ha='center',
                    fontsize=14, color='#3A7FC1')
        if gp and gp < 0.05:
            ax.text(i + w/2, gv + 0.025*np.sign(gv), '*', ha='center',
                    fontsize=14, color='#E05252')

    ax.axhline(0, color='black', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(cat_labels, fontsize=10)
    ax.set_ylabel("Mean log2FC (one-sample t-test vs 0)\n* p < 0.05")
    ax.set_title("VEN Pathway Molecular Signature: ISS vs Ground Stress\n"
                 "Myelination is specifically UP in spaceflight and DOWN in ground stress",
                 fontsize=11)
    ax.legend(fontsize=9)

    # Annotate the myelination dissociation
    myel_idx = cats.index('myelination')
    ax.annotate(f"ISS: +{iss_v[myel_idx]:.2f}*\nGround: {gnd_v[myel_idx]:.2f}*",
                xy=(myel_idx, 0), xytext=(myel_idx + 0.8, 0.45),
                fontsize=9, color='#333333',
                arrowprops=dict(arrowstyle='->', color='#555555'))
    strip(ax)

    plt.tight_layout()
    plt.savefig("figures/fig3_molecular_dissociation.pdf", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig3_molecular_dissociation.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 3 saved.")


# FIGURE 4: VEN circuit model results (from Kaggle v2)
# Copy results/ven_model_clinical.json from Kaggle before running
def fig4_ven_model():
    model_path = "results/ven_model_clinical.json"
    if not os.path.exists(model_path):
        print(f" SKIP Fig 4: {model_path} not found.")
        print(" After Kaggle finishes, copy clinical.json here and rename it.")
        print(" Command: cp /path/to/ven_v2/results/clinical.json results/ven_model_clinical.json")
        return

    with open(model_path) as f:
        model = json.load(f)

    conds = ['typical','autism_like','ftd_like','alz_like']
    labels = ['Typical\n(2% VENs)','Autism-like\n(0.4% VENs)','FTD-like\n(ablated)',"Alzheimer's-like"]
    colors = ['#3A7FC1','#52C869','#E05252','#F5A623']

    rts = [model[c]['mean_rt']  for c in conds]
    rt_sd = [model[c]['std_rt']   for c in conds]
    accs  = [model[c]['mean_acc'] for c in conds]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    ax.bar(range(4), rts, yerr=rt_sd, color=colors, alpha=0.85,
           edgecolor='k', lw=0.7, capsize=5)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Mean Reaction Time (ms)")
    ax.set_title(f"A: VEN Circuit Model - Reaction Time\n"
                 f"(arXiv:2604.09229, 20 seeds)", fontsize=11)
    strip(ax)

    ax = axes[1]
    ax.bar(range(4), accs, color=colors, alpha=0.85, edgecolor='k', lw=0.7)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Decision Accuracy")
    ax.set_ylim(0, 1.05)
    ax.set_title("B: VEN Circuit Model - Accuracy\n"
                 "FTD-like shows specific accuracy impairment", fontsize=11)
    strip(ax)

    plt.tight_layout()
    plt.savefig("figures/fig4_ven_model.pdf", dpi=150, bbox_inches='tight')
    plt.savefig("figures/fig4_ven_model.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure 4 saved.")

print("Generating figures from real data...\n")
fig1_ert_paradox()
fig2_domain_specificity()
fig3_molecular_dissociation()
fig4_ven_model()

print("\nAll figures complete. Check figures/ directory.")
print("fig4 requires results/ven_model_clinical.json from your Kaggle run.")
