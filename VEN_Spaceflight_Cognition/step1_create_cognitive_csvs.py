"""
step1_create_cognitive_csvs.py
Creates structured CSVs from published cognitive performance data.

SOURCE 1: Image_1 heatmap (NASA Twins Study, Fig 10 equivalent)
  Values are TW-HR standardized scores (SD units relative to 15-astronaut baseline)
  Published in: Garrett-Bakelman et al., Science 2019
  Phases: Pre-Flight, In-Flight 1-6 (early, months 1-6),
          In-Flight 7-12 (late, months 7-12), Post-Flight

SOURCE 2: Table_1.docx Table 2 (Dev et al. 2024, Frontiers Physiology)
  Raw mean (SD) scores for 6-month ISS mission astronauts
  Published in: Dev et al. 2024, doi:10.3389/fphys.2024.1451269

Run: python step1_create_cognitive_csvs.py
Output: data/raw/twins_cognitive_heatmap.csv
        data/raw/dev2024_raw_scores.csv
"""

import os
import pandas as pd
import numpy as np
import re

os.makedirs("data/raw", exist_ok=True)

# SOURCE 1: Twins Study heatmap values
# Read directly from Image_1 uploaded  (values printed in cells)
# These are TW-HR score differences in SD units


# Values read directly from Image 1 (NASA Twins Study heatmap figure)
# task, metric, pre_flight, inf_1_6, inf_7_12, post_flight
TWINS_HEATMAP = [
    # SPEED scores
    ("MP",   "speed", -1.5, -0.2,  -0.3,  -2.3),
    ("VOLT", "speed", 0.1,  0.3,   0.9,   0.3),
    ("F2B",  "speed", -0.3,  0.4,  -0.1,  -0.2),
    ("AM", "speed", 0.4,  1.8,   1.2,  -0.4),
    ("LOT",  "speed", -0.5,  0.0,   0.0,  -1.5),
    ("ERT",  "speed", -1.8,  0.9,  -0.9,  -2.3),
    ("MRT",  "speed", 0.3,  1.0,   0.8,  -0.1),
    ("DSST", "speed", -2.0, -0.4,  -0.2,  -0.9),
    ("PVT",  "speed", -1.2,  0.1,  -0.1,  -1.1),
    ("BART", "speed", -0.6,  0.2,  -0.6,  -2.5),
    ("ALL",  "speed", -0.7,  0.4,  -0.1,  -1.1),
    # ACCURACY scores
    ("MP",   "accuracy", -0.9,  0.9,   0.9,   0.0),
    ("VOLT", "accuracy",  0.0, -1.4,  -1.3,  -2.3),
    ("F2B",  "accuracy", -0.7,  0.3,  -0.7,  -1.2),
    ("AM",   "accuracy",  0.9, -0.7,  -2.9,  -3.1),
    ("LOT",  "accuracy", -0.9,  0.2,   0.7,   0.4),
    ("ERT",  "accuracy",  0.0, -0.9,  -0.1,  -0.2),
    ("MRT",  "accuracy", -0.2, -0.4,   0.1,  -0.8),
    ("DSST", "accuracy",  0.2,  0.0,   1.1,   0.8),
    ("PVT",  "accuracy",  0.0, -0.3,  -0.2,  -2.6),
    ("ALL",  "accuracy", -0.2, -0.3,  -0.4,  -1.1),
    ("BART", "accuracy",  0.9,  1.4,   1.1,   0.7),
    ("EFF",  "accuracy", -0.4,  0.0,  -0.3,  -1.1),
]

rows = []
for task, metric, pre, inf16, inf712, post in TWINS_HEATMAP:
    rows.append({"task": task, "metric": metric,
                 "phase": "pre_flight",    "value": pre})
    rows.append({"task": task, "metric": metric,
                 "phase": "inf_early_1_6", "value": inf16})
    rows.append({"task": task, "metric": metric,
                 "phase": "inf_late_7_12", "value": inf712})
    rows.append({"task": task, "metric": metric,
                 "phase": "post_flight",   "value": post})

twins_df = pd.DataFrame(rows)
twins_df.to_csv("data/raw/twins_cognitive_heatmap.csv", index=False)
print(f"Saved: data/raw/twins_cognitive_heatmap.csv ({len(twins_df)} rows)")
print("  Source: Image_1 heatmap — NASA Twins Study (Garrett-Bakelman et al. 2019)")
print("  Values: TW-HR difference scores in SD units")


# SOURCE 2: Dev et al. 2024 - raw scores for 6-month missions
# Parsed from Table_1.docx, Table 2 (already confirmed by reading docx)

def parse_mean_sd(cell_text):
    """Extract mean and SD from '1627.96 (617.83) 712.20 - 3097.40' format."""
    cell_text = cell_text.replace('\n', ' ').strip()
    match = re.match(r'([\d.]+)\s*\(([\d.]+)\)', cell_text)
    if match:
        return float(match.group(1)), float(match.group(2))
    return None, None


# Raw data from Table 2 of Dev et al. 2024 (parsed from docx)
# Speed is in milliseconds (higher = slower)
# Accuracy is proportion 0-1 (higher = better)
DEV2024_RAW = [
    # (task, metric, pre_mean, pre_sd, early_mean, early_sd,
    # late_mean, late_sd, early_post_mean, early_post_sd,
    # late_post_mean, late_post_sd)
    ("VOLT", "speed",
     1627.96, 617.83, 1683.65, 673.4, 1567.78, 460.3,
     1582.56, 549.95, 1621.22, 439.8),
    ("VOLT", "accuracy",
     0.95, 0.06, 0.94, 0.07, 0.94, 0.05,
     0.95, 0.05, 0.92, 0.07),
    ("F2B", "speed",
     552.94, 60.95, 628.04, 97.03, 616.9, 95.39,
     602.06, 83.21, 627.02, 81.66),
    ("F2B", "accuracy",
     0.92, 0.07, 0.94, 0.05, 0.92, 0.07,
     0.91, 0.07, 0.91, 0.07),
    ("AM", "speed",
     2353.32, 877.69, 2264.35, 724.74, 2014.07, 668.01,
     1954.39, 575.82, 1910.12, 417.89),
    ("AM", "accuracy",
     0.81, 0.11, 0.77, 0.15, 0.80, 0.10,
     0.77, 0.15, 0.75, 0.13),
    ("LOT", "speed",
     4691.3, 1183.95, 4871.07, 1218.39, 4980.0, 1124.30,
     4824.74, 1090.94, 4862.98, 997.43),
    ("LOT", "accuracy",
     0.78, 0.08, 0.80, 0.07, 0.79, 0.07,
     0.77, 0.08, 0.78, 0.09),
    ("ERT", "speed",
     2645.3, 886.4, 2621.02, 757.04, 2527.08, 780.23,
     2282.83, 643.33, 2289.49, 573.68),
    ("ERT", "accuracy",
     0.71, 0.11, 0.69, 0.08, 0.71, 0.12,
     0.70, 0.12, 0.75, 0.11),
    ("MRT", "speed",
     8546.15, 2340.55, 9315.28, 2485.21, 8986.65, 2703.81,
     8232.45, 2736.61, 757.36, 2405.76),  # note: late post looks like typo in source
    ("MRT", "accuracy",
     0.79, 0.14, 0.75, 0.11, 0.73, 0.14,
     0.74, 0.15, 0.74, 0.13),
    ("DSST", "speed",
     1209.38, 145.95, 1322.4, 208.59, 1315.22, 204.36,
     1316.29, 236.28, 1255.77, 190.77),
    ("DSST", "accuracy",
     0.98, 0.02, 0.98, 0.02, 0.99, 0.02,
     0.98, 0.02, 0.98, 0.03),
    ("BART", "speed",
     637.83, 429.65, 682.47, 389.64, 629.39, 372.43,
     718.27, 306.3, 631.53, 335.99),
    ("PVT", "speed",
     5.25, 0.24, 5.38, 0.25, 5.30, 0.34,
     5.34, 0.34, 5.34, 0.41),
    ("PVT", "accuracy",
     0.97, 0.04, 0.96, 0.03, 0.97, 0.02,
     0.96, 0.05, 0.96, 0.03),
    ("MPT", "speed",
     1061.18, 153.53, 1083.93, 140.01, 1056.56, 126.8,
     1078.51, 177.39, 996.36, 129.77),
]

phases = ["pre_flight", "early_flight", "late_flight",
          "early_post_flight", "late_post_flight"]

dev_rows = []
for row in DEV2024_RAW:
    task, metric = row[0], row[1]
    values = row[2:]  # mean, sd pairs
    for i, phase in enumerate(phases):
        mean = values[i * 2]
        sd   = values[i * 2 + 1]
        dev_rows.append({
            "task": task, "metric": metric,
            "phase": phase, "mean": mean, "sd": sd,
        })

dev_df = pd.DataFrame(dev_rows)
dev_df.to_csv("data/raw/dev2024_raw_scores.csv", index=False)
print(f"\nSaved: data/raw/dev2024_raw_scores.csv ({len(dev_df)} rows)")
print(" Source: Dev et al. 2024, Table 2 - 6-month ISS missions (n~24)")
print(" Speed in ms (higher = slower), Accuracy 0-1 (higher = better)")

# Compute ERT z-scores relative to preflight for verification
ert_speed = dev_df[(dev_df['task'] == 'ERT') & (dev_df['metric'] == 'speed')].copy()
pre_mean = float(ert_speed[ert_speed['phase'] == 'pre_flight']['mean'].values[0])
pre_sd   = float(ert_speed[ert_speed['phase'] == 'pre_flight']['sd'].values[0])
ert_speed['z_vs_pre'] = (ert_speed['mean'] - pre_mean) / pre_sd
print("\n  ERT Speed z-scores relative to preflight (Dev 2024, 6-month missions):")
for _, row in ert_speed.iterrows():
    print(f" {row['phase']:25s}: mean={row['mean']:.1f}ms  z={row['z_vs_pre']:+.3f}")
print("\n  Note: ERT speed is STABLE in 6-month missions (z-scores near 0)")

print("\nCSV creation complete.")
