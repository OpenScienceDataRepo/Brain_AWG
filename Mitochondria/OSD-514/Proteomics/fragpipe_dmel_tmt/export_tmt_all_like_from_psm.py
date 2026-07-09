#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd


def accession_from_protein(protein: str) -> str:
    if pd.isna(protein):
        return ""
    first = str(protein).split(";")[0]
    parts = first.split("|")
    if len(parts) >= 3:
        return parts[1]
    return first


def plex_to_display(plex: str) -> str:
    return {"TMTA": "TMTa", "TMTB": "TMTb", "TMTC": "TMTc"}[plex]


def plex_to_fraction(plex: str) -> str:
    # Mirror author headers in TMT_all.
    return {"TMTA": "F5", "TMTB": "F5", "TMTC": "F8"}[plex]


def channel_display(ch: str) -> str:
    # Author files use 131 (not 131N) in header text.
    return "131" if str(ch).upper() == "131N" else str(ch).upper()


def load_mapping(assay_file: Path) -> pd.DataFrame:
    m = pd.read_csv(assay_file, sep="\t")
    out = m[["Sample Name", "Parameter Value[Run Number]", "Label"]].copy()
    out.columns = ["sample", "plex_display", "channel"]
    out["plex"] = out["plex_display"].str.upper()
    out["channel"] = out["channel"].astype(str).str.upper()

    # Add bridge pools per plex, as seen in author files.
    pools = []
    for plex in ["TMTA", "TMTB", "TMTC"]:
        pools.append({"sample": "ctrl_pool_1", "plex": plex, "channel": "130C"})
        pools.append({"sample": "ctrl_pool_2", "plex": plex, "channel": "131N"})
    out = pd.concat([out[["sample", "plex", "channel"]], pd.DataFrame(pools)], ignore_index=True)
    return out


def load_protein_metrics(results_dir: Path, plex: str) -> pd.DataFrame:
    fp = results_dir / plex / "protein.tsv"
    cols = [
        "Protein",
        "Protein Description",
        "Coverage",
        "Total Peptides",
        "Total Spectral Count",
        "Unique Peptides",
    ]
    m = pd.read_csv(fp, sep="\t", usecols=cols)
    return m


def aggregate_psm_to_protein(results_dir: Path, plex: str, q_cut: float, p_cut: float, purity_cut: float) -> pd.DataFrame:
    psm = pd.read_csv(results_dir / plex / "psm.tsv", sep="\t")

    # PD-like quality filters.
    filt = (psm["Qvalue"] <= q_cut) & (psm["Probability"] >= p_cut)
    if "Purity" in psm.columns:
        filt &= psm["Purity"].fillna(0) >= purity_cut
    psm = psm[filt].copy()

    reporter_cols = [f"Intensity {plex}_{c}" for c in ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N"]]
    existing_reporters = [c for c in reporter_cols if c in psm.columns]

    keep_cols = ["Protein", "Qvalue", "Probability"] + existing_reporters
    psm = psm[keep_cols].copy()

    # Sum reporter intensities per protein.
    agg_sum = psm.groupby("Protein", as_index=False)[existing_reporters].sum(min_count=1)

    # Per-channel contributing-PSM count (non-null, >0).
    counts = psm[["Protein"] + existing_reporters].copy()
    for c in existing_reporters:
        counts[c] = pd.to_numeric(counts[c], errors="coerce").fillna(0) > 0
    agg_count = counts.groupby("Protein", as_index=False)[existing_reporters].sum()
    agg_count = agg_count.rename(columns={c: c.replace("Intensity ", "Count ") for c in existing_reporters})

    out = agg_sum.merge(agg_count, on="Protein", how="left")
    return out


def write_non_pool(df: pd.DataFrame, out_file: Path, frac: str, ch: str, sample: str):
    ch_disp = channel_display(ch)
    abundance_col = f"Abundance: {frac}: {ch_disp}, Sample, n/a, {sample}"
    count_col = f"Abundances Count: {frac}: {ch_disp}, Sample, n/a, {sample}"

    out = pd.DataFrame(
        {
            "Accession": df["Accession"],
            "Description": df["Protein Description"],
            "Coverage [%]": df["Coverage"],
            "# Peptides": df["Total Peptides"],
            "# PSMs": df["Total Spectral Count"],
            "# Unique Peptides": df["Unique Peptides"],
            abundance_col: df["ABUNDANCE"],
            count_col: df["ABUND_COUNT"],
            "Modifications": "",
        }
    )
    out.to_csv(out_file, sep="\t", index=False)


def write_pool(df: pd.DataFrame, out_file: Path, frac: str):
    h1 = f"Abundance: {frac}: 130C, Sample, n/a, ctrl_pool_1"
    h2 = f"Abundance: {frac}: 131, Sample, n/a, ctrl_pool_2"

    out = pd.DataFrame(
        {
            "Accession": df["Accession"],
            "Description": df["Protein Description"],
            "Coverage [%]": df["Coverage"],
            "# Peptides": df["Total Peptides"],
            "# PSMs": df["Total Spectral Count"],
            "# Unique Peptides": df["Unique Peptides"],
            h1: df["ABUNDANCE_130C"],
            h2: df["ABUNDANCE_131N"],
            "Modifications": "",
        }
    )
    out.to_csv(out_file, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser(description="Export TMT_all-like files using PD-like aggregation from psm.tsv")
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--assay-file", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--qvalue", type=float, default=0.01)
    ap.add_argument("--probability", type=float, default=0.90)
    ap.add_argument("--purity", type=float, default=0.50)
    args = ap.parse_args()

    mapping = load_mapping(args.assay_file)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for plex in ["TMTA", "TMTB", "TMTC"]:
        disp = plex_to_display(plex)
        frac = plex_to_fraction(plex)

        agg = aggregate_psm_to_protein(args.results_dir, plex, args.qvalue, args.probability, args.purity)
        metrics = load_protein_metrics(args.results_dir, plex)
        merged = metrics.merge(agg, on="Protein", how="left")
        merged["Accession"] = merged["Protein"].map(accession_from_protein)

        # Biological sample files.
        smap = mapping[(mapping["plex"] == plex) & (~mapping["sample"].str.startswith("ctrl_pool"))]
        for _, r in smap.iterrows():
            sample = str(r["sample"])
            ch = str(r["channel"]).upper()
            i_col = f"Intensity {plex}_{ch}"
            c_col = f"Count {plex}_{ch}"
            if i_col not in merged.columns:
                continue

            tmp = merged.copy()
            tmp["ABUNDANCE"] = pd.to_numeric(tmp[i_col], errors="coerce")
            tmp["ABUND_COUNT"] = pd.to_numeric(tmp.get(c_col, 0), errors="coerce").fillna(0).astype(int)
            tmp = tmp[tmp["ABUNDANCE"].notna() & (tmp["ABUNDANCE"] > 0)].copy()
            if tmp.empty:
                continue

            out_file = args.out_dir / f"{sample}_{disp}.txt"
            write_non_pool(tmp, out_file, frac, ch, sample)
            written.append(out_file)

        # Pool file (130C + 131N columns).
        c1 = f"Intensity {plex}_130C"
        c2 = f"Intensity {plex}_131N"
        if c1 in merged.columns and c2 in merged.columns:
            tmp = merged.copy()
            tmp["ABUNDANCE_130C"] = pd.to_numeric(tmp[c1], errors="coerce")
            tmp["ABUNDANCE_131N"] = pd.to_numeric(tmp[c2], errors="coerce")
            tmp = tmp[(tmp["ABUNDANCE_130C"].fillna(0) > 0) | (tmp["ABUNDANCE_131N"].fillna(0) > 0)].copy()
            if not tmp.empty:
                out_file = args.out_dir / f"pool_{disp}.txt"
                write_pool(tmp, out_file, frac)
                written.append(out_file)

    print(f"Wrote {len(written)} files to {args.out_dir}")
    for p in sorted(written):
        print(f"  - {p.name}")


if __name__ == "__main__":
    main()
