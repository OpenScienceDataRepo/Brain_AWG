"""
step2_run_analysis.py
ISS rodent (GSE239336) + ground analogue (OSD-202) + cognitive analysis

requires step1 to be run first:
  data/raw/twins_cognitive_heatmap.csv
  data/raw/dev2024_raw_scores.csv
  data/raw/GSE239336_FCT_GCvsFLT-SAL_DEanalysis.txt
  data/raw/GLDS-202_rna_seq_differential_expression_GLbulkRNAseq.csv
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats

os.makedirs("results", exist_ok=True)

N_PERM = 10000
RNG_SEED = 42
rng = np.random.default_rng(RNG_SEED)

# BCL11B is the official HGNC symbol for what is sometimes called CTIP2
VEN_GENES = {
    "myelination": ["MBP", "MOG", "PLP1", "MAG", "CNP", "MOBP", "ERMN"],
    "fast_signalling": ["SCN1A", "KCNQ2", "ANK3", "NEFH", "NEFM", "NEFL", "SNCG"],
    "social_circuit": ["OXTR", "AVPR1A", "HTR2A", "DRD1", "CHRM1", "GABRB2"],
    "layer5_proj": ["FEZF2", "BCL11B", "TBR1", "SATB2", "CUX1"],
    "metabolic": ["VDAC1", "ATP2B2", "SLC17A7", "SNAP25", "SYP", "NRXN1"],
}
ALL_VEN_GENES = [g for gs in VEN_GENES.values() for g in gs]

# curated VEN markers from peer-reviewed literature for cross-reference

OSPINA_MARKERS = ["VAT1L", "CHST8", "LYPD1", "SULF2"]

KIM_HODGE_MARKERS = ["GABRQ", "ADRA1A", "VMAT2"]
CURATED_MARKERS = OSPINA_MARKERS + KIM_HODGE_MARKERS

TASK_DOMAINS = {
    "ERT": "social", "VOLT": "memory", "F2B": "working_memory",
    "AM": "reasoning", "LOT": "spatial", "MRT": "reasoning",
    "DSST": "processing", "PVT": "attention", "BART": "risk",
    "MP": "motor", "MPT": "motor", "ALL": "composite", "EFF": "efficiency",
}


def two_tailed_perm(obs_values, background, n_perm=N_PERM, rng=rng):
    """Two-tailed permutation test against genome-wide null."""
    if len(obs_values) == 0:
        return None, None, None, None

    obs_mean = np.mean(obs_values)
    n = len(obs_values)

    null_means = np.array([
        rng.choice(background, size=n, replace=False).mean()
        for _ in range(n_perm)
    ])
    null_mean = null_means.mean()
    null_sd = null_means.std()

    sd_above = (obs_mean - null_mean) / null_sd if null_sd > 0 else np.nan
    p_perm = float(np.mean(
        np.abs(null_means - null_mean) >= np.abs(obs_mean - null_mean)
    ))
    return obs_mean, sd_above, p_perm, null_mean


def sig_label(p):
    if p is None or np.isnan(p):
        return ""
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    if p < 0.10: return "+"
    return ""


def run_panel(gene_dict, fc_series, background, label):
    rows = []
    for cat, genes in gene_dict.items():
        present = [g for g in genes if g in fc_series.index and not np.isnan(fc_series[g])]
        vals = fc_series[present].values.astype(float) if present else np.array([])

        row = {
            "dataset": label,
            "category": cat,
            "n_panel": len(genes),
            "n_present": len(present),
            "genes_present": ", ".join(present),
            "genes_missing": ", ".join([g for g in genes if g not in fc_series.index]),
            "mean_log2FC": float(np.mean(vals)) if len(vals) > 0 else None,
        }

        if len(vals) >= 2:
            t, p_t = stats.ttest_1samp(vals, 0)
            row["t_stat"] = round(float(t), 3)
            row["p_ttest"] = round(float(p_t), 4)
        else:
            row["t_stat"] = None
            row["p_ttest"] = None

        obs_mean, sd_above, p_perm, _ = two_tailed_perm(vals, background)
        row["SD_above_null"] = round(sd_above, 3) if sd_above is not None else None
        row["perm_p"] = round(p_perm, 4) if p_perm is not None else None
        row["sig"] = sig_label(p_perm)
        rows.append(row)
    return pd.DataFrame(rows)


# cognitive analysis
print("loading twins cognitive data (N=1, hypothesis-generating)")

twins = pd.read_csv("data/raw/twins_cognitive_heatmap.csv")
twins["domain"] = twins["task"].map(TASK_DOMAINS)

twins_wide = twins.pivot_table(
    index=["task", "metric", "domain"],
    columns="phase", values="value"
).reset_index()

twins_wide["early_to_late_change"] = (
    twins_wide["inf_late_7_12"] - twins_wide["inf_early_1_6"]
)
twins_wide.to_csv("results/twins_cognitive_wide.csv", index=False)

# ERT speed specificity
speed = twins_wide[
    (twins_wide["metric"] == "speed") &
    (~twins_wide["task"].isin(["ALL", "EFF"]))
].copy()

ert_row = speed[speed["task"] == "ERT"].iloc[0]
ert_change = float(ert_row["early_to_late_change"])
ert_early = float(ert_row["inf_early_1_6"])
ert_late = float(ert_row["inf_late_7_12"])

other_changes = speed[speed["task"] != "ERT"]["early_to_late_change"].values
se_other = np.std(other_changes, ddof=1) / np.sqrt(len(other_changes))
t_cross_domain = (ert_change - np.mean(other_changes)) / se_other

print(f"ERT early->late speed change: {ert_change:+.1f} SD")
print(f"cross-domain t: {t_cross_domain:.2f} SEs below person-level mean")

# 6-month comparison
dev = pd.read_csv("data/raw/dev2024_raw_scores.csv")
ert_dev = dev[(dev["task"] == "ERT") & (dev["metric"] == "speed")].copy()
pre_mean = float(ert_dev[ert_dev["phase"] == "pre_flight"]["mean"].values[0])
pre_sd = float(ert_dev[ert_dev["phase"] == "pre_flight"]["sd"].values[0])
ert_dev["z"] = -(ert_dev["mean"] - pre_mean) / pre_sd
dev_early = float(ert_dev[ert_dev["phase"] == "early_flight"]["z"].values[0])
dev_late = float(ert_dev[ert_dev["phase"] == "late_flight"]["z"].values[0])
dev_change = dev_late - dev_early
duration_effect = ert_change - dev_change

print(f"6-month ERT change (N=24): {dev_change:+.3f} SD")
print(f"duration effect: {duration_effect:+.2f} SD")

ert_comparison = {
    "twins_340day_N1": {
        "early_inflight": ert_early,
        "late_inflight": ert_late,
        "early_to_late_change": ert_change,
        "note": "N=1, hypothesis-generating only"
    },
    "dev2024_6month_N24": {
        "early_inflight": dev_early,
        "late_inflight": dev_late,
        "early_to_late_change": dev_change,
        "note": "N=24, Dev et al. 2024"
    },
    "duration_effect_SD": duration_effect,
}
with open("results/ert_duration_comparison.json", "w") as f:
    json.dump(ert_comparison, f, indent=2)


# GSE239336 ISS rodent
print("loading GSE239336")

gse_path = "data/raw/GSE239336_FCT_GCvsFLT-SAL_DEanalysis.txt"
if not os.path.exists(gse_path):
    raise FileNotFoundError(f"missing: {gse_path}")

gse_raw = pd.read_csv(gse_path, sep="\t", skiprows=6, low_memory=False)
gse_raw.columns = [c.strip() for c in gse_raw.columns]
gse_raw = gse_raw.rename(columns={"Target name": "gene", "Log2": "log2fc", "Pvalue": "pvalue"})
gse_raw = gse_raw[["gene", "log2fc", "pvalue"]].dropna(subset=["gene", "log2fc"])
gse_raw["gene"] = gse_raw["gene"].astype(str).str.strip()
gse_raw["gene_up"] = gse_raw["gene"].str.upper()
gse_raw["log2fc"] = pd.to_numeric(gse_raw["log2fc"], errors="coerce")
gse_raw = gse_raw.dropna(subset=["log2fc"])
print(f"  {len(gse_raw)} genes, genome-wide mean={gse_raw['log2fc'].mean():+.4f}")

gse_fc = gse_raw.set_index("gene_up")["log2fc"]
gse_bg = gse_raw["log2fc"].values

gse_df = run_panel(VEN_GENES, gse_fc, gse_bg, "GSE239336_ISS")
gse_df.to_csv("results/gse239336_ven_panel.csv", index=False)

pd.DataFrame([
    {"gene": g, "category": cat, "log2fc": float(gse_fc[g]), "dataset": "GSE239336_ISS"}
    for cat, genes in VEN_GENES.items()
    for g in genes if g in gse_fc.index
]).to_csv("results/gse239336_ven_genes.csv", index=False)

gse_json = {r["category"]: r for _, r in gse_df.iterrows()}
with open("results/gse239336_ven_signature.json", "w") as f:
    json.dump({k: dict(v) for k, v in gse_json.items()}, f, indent=2)

print("GSE239336 results:")
for _, r in gse_df.iterrows():
    fc = f"{r['mean_log2FC']:+.3f}" if r["mean_log2FC"] is not None else "N/A"
    pp = f"{r['perm_p']:.4f}" if r["perm_p"] is not None else "N/A"
    sd = f"{r['SD_above_null']:+.2f}" if r["SD_above_null"] is not None else "N/A"
    print(f"  {r['category']:20s}  FC={fc}  perm_p={pp}  SD={sd}  {r['sig']}")


# OSD-202 ground analogue
print("loading OSD-202")

osd_path = "data/raw/GLDS-202_rna_seq_differential_expression_GLbulkRNAseq.csv"
if not os.path.exists(osd_path):
    raise FileNotFoundError(f"missing: {osd_path}")

osd_raw = pd.read_csv(osd_path, low_memory=False)
contrast = "(cobalt-57 gamma radiation & Hindlimb Unloaded & 1 month)v(non-irradiated & Normally Loaded Control & 1 month)"
osd_de = osd_raw[["SYMBOL", f"Log2fc_{contrast}", f"P.value_{contrast}"]].copy()
osd_de.columns = ["gene", "log2fc", "pvalue"]
osd_de = osd_de.dropna(subset=["gene", "log2fc"])
osd_de["gene_up"] = osd_de["gene"].str.upper()
osd_de["log2fc"] = pd.to_numeric(osd_de["log2fc"], errors="coerce")
osd_de = osd_de.dropna(subset=["log2fc"])
print(f"  {len(osd_de)} genes, genome-wide mean={osd_de['log2fc'].mean():+.4f}")

osd_fc = osd_de.set_index("gene_up")["log2fc"]
osd_bg = osd_de["log2fc"].values

osd_df = run_panel(VEN_GENES, osd_fc, osd_bg, "OSD202_Ground")
osd_df.to_csv("results/osd202_ven_panel.csv", index=False)

osd_json = {r["category"]: r for _, r in osd_df.iterrows()}
with open("results/osd202_ven_signature.json", "w") as f:
    json.dump({k: dict(v) for k, v in osd_json.items()}, f, indent=2)

print("OSD-202 results:")
for _, r in osd_df.iterrows():
    fc = f"{r['mean_log2FC']:+.3f}" if r["mean_log2FC"] is not None else "N/A"
    pp = f"{r['perm_p']:.4f}" if r["perm_p"] is not None else "N/A"
    sd = f"{r['SD_above_null']:+.2f}" if r["SD_above_null"] is not None else "N/A"
    print(f"  {r['category']:20s}  FC={fc}  perm_p={pp}  SD={sd}  {r['sig']}")


# cross-dataset comparison
comp_rows = []
for cat in VEN_GENES:
    i = gse_json.get(cat, {})
    g = osd_json.get(cat, {})
    comp_rows.append({
        "category": cat,
        "ISS_mean_log2FC": i.get("mean_log2FC"),
        "ISS_perm_p": i.get("perm_p"),
        "ISS_SD_null": i.get("SD_above_null"),
        "ISS_sig": i.get("sig"),
        "Gnd_mean_log2FC": g.get("mean_log2FC"),
        "Gnd_perm_p": g.get("perm_p"),
        "Gnd_SD_null": g.get("SD_above_null"),
        "Gnd_sig": g.get("sig"),
    })
pd.DataFrame(comp_rows).to_csv("results/molecular_comparison.csv", index=False)


# curated VEN gene cross-reference (Windy feedback)
# Ospina-Perez et al. 2019 and Kim/Hodge 2018/2020
# descriptive only, not used to modify the pre-registered panel
print("checking curated VEN markers from literature")

curated_rows = []
for gene in CURATED_MARKERS:
    g_up = gene.upper()
    src = "Ospina-Perez 2019" if gene in OSPINA_MARKERS else "Kim 2018/Hodge 2020"
    row = {
        "gene": gene,
        "source": src,
        "in_VEN_panel": gene in ALL_VEN_GENES,
        "GSE239336_detected": g_up in gse_fc.index,
        "GSE239336_log2FC": float(gse_fc[g_up]) if g_up in gse_fc.index else None,
        "OSD202_detected": g_up in osd_fc.index,
        "OSD202_log2FC": float(osd_fc[g_up]) if g_up in osd_fc.index else None,
    }
    curated_rows.append(row)
    iss = f"{row['GSE239336_log2FC']:+.3f}" if row["GSE239336_log2FC"] is not None else "not detected"
    gnd = f"{row['OSD202_log2FC']:+.3f}" if row["OSD202_log2FC"] is not None else "not detected"
    print(f"  {gene:10s} ({src})  ISS={iss}  Ground={gnd}")

pd.DataFrame(curated_rows).to_csv("results/curated_VEN_geneset_crossref.csv", index=False)


# save summary
myel_i = gse_json.get("myelination", {})
myel_g = osd_json.get("myelination", {})
soc_i = gse_json.get("social_circuit", {})

main_findings = {
    "permutation_type": "two-tailed, N=10000, seed=42, numpy.default_rng",
    "H1_ISS_myelination_FC": myel_i.get("mean_log2FC"),
    "H1_ISS_myelination_t_p": myel_i.get("p_ttest"),
    "H1_ISS_myelination_perm_p": myel_i.get("perm_p"),
    "H1_ISS_myelination_SD_null": myel_i.get("SD_above_null"),
    "H1_Gnd_myelination_perm_p": myel_g.get("perm_p"),
    "H1_Gnd_myelination_SD_null": myel_g.get("SD_above_null"),
    "H1_ISS_social_circuit_perm_p": soc_i.get("perm_p"),
    "H2_ERT_change_340day_N1": ert_change,
    "H2_ERT_change_6month_N24": dev_change,
    "H2_duration_effect_SD": duration_effect,
    "N1_caveat": "340-day: N=1 astronaut, hypothesis-generating; replicated in N=24 Dev 2024",
    "curated_geneset_refs": [
        "Ospina-Perez et al. 2019 Cereb Cortex DOI:10.1093/cercor/bhy286",
        "Kim et al. 2018 Cereb Cortex PMC6075576",
        "Hodge et al. 2020 Nat Commun DOI:10.1038/s41467-020-14952-3",
    ],
}
with open("results/main_findings.json", "w") as f:
    json.dump(main_findings, f, indent=2)

print("done. all results saved to results/")
