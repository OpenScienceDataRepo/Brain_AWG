#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


def find_tmt_protein_table(workdir: Path):
    candidates = []
    patterns = [
        "*abundance_protein*.tsv",     # e.g. tmt-report/abundance_protein_MD.tsv
        "*ratio_protein*.tsv",         # e.g. tmt-report/ratio_protein_MD.tsv
        "*protein*report*.tsv",
        "*tmt*protein*.tsv",
        "*protein*.tsv",
    ]
    for pat in patterns:
        candidates.extend(workdir.rglob(pat))

    # Keep TSV files only and prioritize abundance protein table from tmt-report.
    uniq = []
    seen = set()
    for p in candidates:
        if p.suffix.lower() != ".tsv":
            continue
        s = str(p)
        if s in seen:
            continue
        seen.add(s)
        uniq.append(p)

    def rank(p: Path):
        n = p.name.lower()
        full = str(p).lower()
        return (
            "abundance_protein" not in n,
            "tmt-report" not in full,
            "ratio_" in n,
            "combined_" in n,
            len(full),
            full,
        )

    candidates = sorted(uniq, key=rank)
    return candidates[0] if candidates else None


def load_channel_rename_map(workdir: Path):
    channel_map = {}
    for plex in ("TMTA", "TMTB", "TMTC"):
        ann = workdir / plex / f"{plex}_annotation.txt"
        if not ann.exists():
            continue
        with ann.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = re.split(r"\s+", line)
                if len(parts) < 2:
                    continue
                label = parts[0].strip().upper()
                sample = parts[1].strip()
                if sample.upper() == "NA":
                    continue
                channel_map[f"{plex}_{label}"] = sample
    return channel_map


def main():
    ap = argparse.ArgumentParser(description="Extract a clean protein x sample matrix from FragPipe TMT output")
    ap.add_argument("--workdir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    inp = find_tmt_protein_table(args.workdir)
    if inp is None:
        raise SystemExit(f"ERROR: could not find TMT protein table under {args.workdir}")

    rename_map = load_channel_rename_map(args.workdir)

    with inp.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames:
            raise SystemExit(f"ERROR: no header in {inp}")

        fields = r.fieldnames
        id_col = "Protein" if "Protein" in fields else fields[0]
        gene_col = "Gene" if "Gene" in fields else ("Mapped Genes" if "Mapped Genes" in fields else None)

        # Keep abundance/intensity-like sample columns.
        tmt_channel_re = re.compile(r"^(TMT[A-Z0-9]+)_(126|127N|127C|128N|128C|129N|129C|130N|130C|131N?|132N|132C|133N|133C|134N)$", re.IGNORECASE)
        sample_cols = []
        for c in fields:
            cl = c.lower()
            if "count" in cl:
                continue
            if any(k in cl for k in ["abundance", "intensity", "reporter", "ratio"]):
                sample_cols.append(c)
                continue
            if tmt_channel_re.match(c):
                sample_cols.append(c)

        if not sample_cols:
            raise SystemExit("ERROR: no abundance/intensity-like columns found in TMT table")

        rows = list(r)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as w:
        wr = csv.writer(w, delimiter="\t")
        matrix_cols = [id_col] + ([gene_col] if gene_col else []) + sample_cols
        out_head = [id_col] + ([gene_col] if gene_col else []) + [rename_map.get(c, c) for c in sample_cols]
        wr.writerow(out_head)
        for row in rows:
            wr.writerow([row.get(c, "") for c in matrix_cols])

    print(f"Input table: {inp}")
    print(f"Output matrix: {args.out}")
    print(f"Rows: {len(rows)} | Sample columns: {len(sample_cols)}")
    if rename_map:
        print(f"Applied metadata label remap for {len(rename_map)} channels from annotation files")


if __name__ == "__main__":
    main()
