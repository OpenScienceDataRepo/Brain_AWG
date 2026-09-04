#!/usr/bin/env python3
import argparse
import glob
import os
import re
from pathlib import Path
import csv


def plex_from_filename(name: str) -> str:
    up = name.upper()
    if "TMTA" in up:
        return "TMTA"
    if "TMTB" in up:
        return "TMTB"
    if "TMTC" in up:
        return "TMTC"
    raise ValueError(f"Cannot infer plex from filename: {name}")


def replicate_from_plex(plex: str) -> int:
    return {"TMTA": 1, "TMTB": 2, "TMTC": 3}[plex]


def infer_condition(sample_name: str) -> str:
    s = sample_name.lower()
    if "spaceug" in s:
        return "SPACEFLIGHT_MICROGRAVITY"
    if "space1g" in s:
        return "SPACEFLIGHT_1G"
    if "agc" in s:
        return "EARTH"
    if "pool" in s:
        return "BRIDGE_POOL"
    return "UNKNOWN"


def infer_sample_id(sample_name: str, plex: str) -> str:
    # Stable, compact sample ID used by TMT-Integrator.
    return f"{sample_name}_{plex}".replace(" ", "_")


def canonical_plex(raw: str) -> str:
    s = str(raw).strip().upper()
    if s in {"TMTA", "TMTA"}:
        return "TMTA"
    if s in {"TMTB"}:
        return "TMTB"
    if s in {"TMTC"}:
        return "TMTC"
    return s


def parse_assay_labels(assay_file: str):
    plex_to_channel_to_sample = {}
    with open(assay_file, "r", errors="ignore", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_name = (row.get("Sample Name") or row.get("Labeled Extract Name") or "").strip()
            # Some ISA-Tab exports place plex in Run Number and channel in Label.
            plex_raw = (row.get("Parameter Value[Run Number]") or row.get("Run Number") or row.get("Label") or "").strip()
            channel = (row.get("Label") or row.get("Term Accession Number") or "").strip().upper()
            if not sample_name or not plex_raw or not channel:
                continue

            plex = canonical_plex(plex_raw)
            if plex not in {"TMTA", "TMTB", "TMTC"}:
                continue

            sample = infer_sample_id(sample_name, plex)
            plex_map = plex_to_channel_to_sample.setdefault(plex, {})
            plex_map.setdefault(channel, sample)
    return plex_to_channel_to_sample


def parse_tmt_all_headers(tmt_all_dir: str):
    tmt_files = sorted(glob.glob(os.path.join(tmt_all_dir, "*.txt")))
    if not tmt_files:
        raise SystemExit(f"No .txt files found in {tmt_all_dir}")

    # Extract all entries matching: Abundance: F#: CHANNEL, Sample, n/a, SAMPLE_NAME
    patt = re.compile(r"Abundance:\s*[^:]*:\s*([^,]+),\s*Sample,\s*n/a,\s*([^\"\t]+)")

    # FragPipe auto-discovered *annotation.txt format expects 2 whitespace-separated columns:
    #   <TMT_label> <sample_name>
    # where labels are channels like 126, 127N, ... and sample_name can be NA to ignore a channel.
    plex_to_channel_to_sample = {}
    for fp in tmt_files:
        fn = os.path.basename(fp)
        plex = plex_from_filename(fn)
        with open(fp, "r", errors="ignore") as f:
            header = f.readline().strip()

        matches = patt.findall(header)
        if not matches:
            continue

        for channel_raw, sample_raw in matches:
            channel = channel_raw.strip().upper()
            sample_name = sample_raw.strip()
            sample = infer_sample_id(sample_name, plex)

            plex_map = plex_to_channel_to_sample.setdefault(plex, {})
            plex_map.setdefault(channel, sample)

    return plex_to_channel_to_sample


def main() -> None:
    ap = argparse.ArgumentParser(description="Build FragPipe TMT annotation from TMT_all header files")
    ap.add_argument("--tmt-all-dir", required=True)
    ap.add_argument("--assay-file", required=False, help="OSD assay metadata TSV for authoritative labeling")
    ap.add_argument("--out", required=False, help="Optional single annotation output file")
    ap.add_argument("--mzml-dir", required=False, help="If set, write per-plex files: <PLEX>_annotation.txt")
    args = ap.parse_args()

    plex_to_channel_to_sample = {}
    if args.assay_file and os.path.isfile(args.assay_file):
        plex_to_channel_to_sample = parse_assay_labels(args.assay_file)
        print(f"Using assay metadata labels from: {args.assay_file}")

    if not plex_to_channel_to_sample:
        plex_to_channel_to_sample = parse_tmt_all_headers(args.tmt_all_dir)
        print("Using labels parsed from TMT_all headers")

    if not plex_to_channel_to_sample:
        raise SystemExit("No abundance descriptors parsed from TMT_all headers.")

    # TMT10 labels in canonical order; missing labels are set to NA so they are ignored.
    ordered_labels = ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N"]
    wrote = []

    # Optional single output (uses TMTA if available, otherwise first plex).
    if args.out:
        pref = "TMTA" if "TMTA" in plex_to_channel_to_sample else sorted(plex_to_channel_to_sample.keys())[0]
        rows = [(label, plex_to_channel_to_sample[pref].get(label, "NA")) for label in ordered_labels]
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w") as w:
            for r in rows:
                w.write("\t".join(map(str, r)) + "\n")
        wrote.append(str(out))

    # Preferred mode: one file per plex, e.g. TMTA_annotation.txt.
    if args.mzml_dir:
        mzml_dir = Path(args.mzml_dir)
        mzml_dir.mkdir(parents=True, exist_ok=True)
        for plex in sorted(plex_to_channel_to_sample.keys()):
            rows = [(label, plex_to_channel_to_sample[plex].get(label, "NA")) for label in ordered_labels]
            out = mzml_dir / f"{plex}_annotation.txt"
            with out.open("w") as w:
                for r in rows:
                    w.write("\t".join(map(str, r)) + "\n")
            wrote.append(str(out))

    if not wrote:
        raise SystemExit("Please provide --out and/or --mzml-dir")

    print("Wrote annotation files:")
    for p in wrote:
        print(f"  - {p}")


if __name__ == "__main__":
    main()
