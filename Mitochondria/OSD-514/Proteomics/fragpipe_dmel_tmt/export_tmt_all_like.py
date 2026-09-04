#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd


def normalize_plex(label: str) -> str:
    s = str(label).strip().upper()
    if s == "TMTA":
        return "TMTA"
    if s == "TMTB":
        return "TMTB"
    if s == "TMTC":
        return "TMTC"
    raise ValueError(f"Unknown plex label: {label}")


def plex_to_display(plex: str) -> str:
    return {"TMTA": "TMTa", "TMTB": "TMTb", "TMTC": "TMTc"}[plex]


def plex_to_fraction(plex: str) -> str:
    # Mirrors author header style seen in TMT_all.
    return {"TMTA": "F5", "TMTB": "F5", "TMTC": "F8"}[plex]


def accession_from_protein(protein: str) -> str:
    if pd.isna(protein):
        return ""
    first = str(protein).split(";")[0]
    parts = first.split("|")
    if len(parts) >= 3:
        return parts[1]
    return first


def load_plex_protein_metrics(results_dir: Path, plex: str) -> pd.DataFrame:
    fp = results_dir / plex / "protein.tsv"
    if not fp.exists():
        raise FileNotFoundError(fp)

    cols = [
        "Protein",
        "Protein Description",
        "Coverage",
        "Total Peptides",
        "Total Spectral Count",
        "Unique Peptides",
    ]
    df = pd.read_csv(fp, sep="\t", usecols=cols)
    return df


def build_sample_mapping(assay_file: Path) -> pd.DataFrame:
    m = pd.read_csv(assay_file, sep="\t")
    out = m[["Sample Name", "Parameter Value[Run Number]", "Label"]].copy()
    out = out.rename(
        columns={
            "Sample Name": "sample",
            "Parameter Value[Run Number]": "plex_display",
            "Label": "channel",
        }
    )
    out["plex"] = out["plex_display"].str.upper().map({"TMTA": "TMTA", "TMTB": "TMTB", "TMTC": "TMTC"})
    out["channel"] = out["channel"].astype(str).str.upper()
    out = out.dropna(subset=["sample", "plex", "channel"]).copy()

    # Add the two bridge pool channels per plex to mimic author outputs.
    pools = []
    for plex in ["TMTA", "TMTB", "TMTC"]:
        pools.append({"sample": "pool_1", "plex": plex, "channel": "130C"})
        pools.append({"sample": "pool_2", "plex": plex, "channel": "131N"})
    out = pd.concat([out, pd.DataFrame(pools)], ignore_index=True)

    return out


def write_non_pool_file(df: pd.DataFrame, out_file: Path, frac: str, channel: str, sample: str):
    abundance_col = f"Abundance: {frac}: {channel}, Sample, n/a, {sample}"
    count_col = f"Abundances Count: {frac}: {channel}, Sample, n/a, {sample}"

    out = pd.DataFrame(
        {
            "Accession": df["Accession"],
            "Description": df["Protein Description"],
            "Coverage [%]": df["Coverage"],
            "# Peptides": df["Total Peptides"],
            "# PSMs": df["Total Spectral Count"],
            "# Unique Peptides": df["Unique Peptides"],
            abundance_col: df["ABUNDANCE"],
            count_col: df["NumberPSM"],
            "Modifications": "",
        }
    )
    out.to_csv(out_file, sep="\t", index=False)


def write_pool_file(df: pd.DataFrame, out_file: Path, frac: str):
    c1 = "130C"
    c2 = "131N"
    h1 = f"Abundance: {frac}: {c1}, Sample, n/a, ctrl_pool_1"
    h2 = f"Abundance: {frac}: {c2}, Sample, n/a, ctrl_pool_2"

    out = pd.DataFrame(
        {
            "Accession": df["Accession"],
            "Description": df["Protein Description"],
            "Coverage [%]": df["Coverage"],
            "# Peptides": df["Total Peptides"],
            "# PSMs": df["Total Spectral Count"],
            "# Unique Peptides": df["Unique Peptides"],
            h1: df.get("ABUNDANCE_130C", pd.NA),
            h2: df.get("ABUNDANCE_131N", pd.NA),
            "Modifications": "",
        }
    )
    out.to_csv(out_file, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser(description="Export FragPipe results into TMT_all-like per-sample TXT tables")
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--assay-file", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    ab_file = args.results_dir / "tmt-report" / "abundance_protein_MD.tsv"
    if not ab_file.exists():
        raise SystemExit(f"Missing abundance table: {ab_file}")

    ab = pd.read_csv(ab_file, sep="\t")
    ab["Accession"] = ab["Protein"].map(accession_from_protein)

    sample_map = build_sample_mapping(args.assay_file)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    written = []

    for plex in ["TMTA", "TMTB", "TMTC"]:
        frac = plex_to_fraction(plex)
        disp = plex_to_display(plex)

        metrics = load_plex_protein_metrics(args.results_dir, plex)
        merged = ab.merge(metrics, on="Protein", how="left", suffixes=("", "_metrics"))

        # Write 6 biological sample files for this plex.
        sub = sample_map[(sample_map["plex"] == plex) & (~sample_map["sample"].str.startswith("pool_"))]
        for _, r in sub.iterrows():
            ch = r["channel"]
            sample = str(r["sample"])
            ch_col = f"{plex}_{ch}"
            if ch_col not in merged.columns:
                continue

            tmp = merged.copy()
            tmp["ABUNDANCE"] = pd.to_numeric(tmp[ch_col], errors="coerce")
            tmp = tmp[tmp["ABUNDANCE"].notna()].copy()
            if tmp.empty:
                continue

            out_file = args.out_dir / f"{sample}_{disp}.txt"
            write_non_pool_file(tmp, out_file, frac, ch, sample)
            written.append(out_file)

        # Write one pool file per plex (two abundance columns).
        ch1 = f"{plex}_130C"
        ch2 = f"{plex}_131N"
        tmp = merged.copy()
        if ch1 in tmp.columns:
            tmp["ABUNDANCE_130C"] = pd.to_numeric(tmp[ch1], errors="coerce")
        if ch2 in tmp.columns:
            tmp["ABUNDANCE_131N"] = pd.to_numeric(tmp[ch2], errors="coerce")
        keep = tmp[[c for c in ["ABUNDANCE_130C", "ABUNDANCE_131N"] if c in tmp.columns]].notna().any(axis=1)
        tmp = tmp[keep].copy()
        if not tmp.empty:
            out_file = args.out_dir / f"pool_{disp}.txt"
            write_pool_file(tmp, out_file, frac)
            written.append(out_file)

    print(f"Wrote {len(written)} files to: {args.out_dir}")
    for p in sorted(written):
        print(f"  - {p.name}")


if __name__ == "__main__":
    main()
