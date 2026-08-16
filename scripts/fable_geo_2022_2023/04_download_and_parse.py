#!/usr/bin/env python3
"""
Download the open, processed expression matrices for the analysis cohorts and
parse the full per-sample metadata table from each GEO series matrix.

Cohorts (human lung NSCLC, treated with immune checkpoint inhibitors, with
per-patient outcome annotation):
  GSE161537  n=82  NSCLC, immunotherapy   -> log2 CPM CSV   (RECIST / OS / PFS)
  GSE162520  n=92  NSCLC, PD-1/PD-L1      -> log2 CPM CSV   (OS / PFS)
  GSE207422  n bulk=? NSCLC neoadjuvant anti-PD-1+chemo -> log2 TPM (MPR/NMPR)

For each series we save:
  results/fable_geo_2022_2023/data/<GSE>_expr.<ext>        (raw download)
  results/fable_geo_2022_2023/data/<GSE>_sample_meta.csv   (parsed GSM table)
"""
import gzip
import io
import re
import sys
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "results" / "fable_geo_2022_2023"
DATADIR = OUTDIR / "data"
DATADIR.mkdir(parents=True, exist_ok=True)
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

EXPR_FILES = {
    "GSE161537": "GSE161537_nivobio_log2cpm.csv.gz",
    "GSE162520": "GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz",
    "GSE207422": "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
}


def matrix_prefix(acc):
    num = re.match(r"GSE(\d+)", acc).group(1)
    return "GSE" + (num[:-3] + "nnn" if len(num) > 3 else "nnn")


def download(acc, fname):
    url = f"{FTP}/{matrix_prefix(acc)}/{acc}/suppl/{fname}"
    dest = DATADIR / f"{acc}__{fname}"
    if dest.exists():
        print(f"exists {dest.name}", file=sys.stderr)
        return dest
    print(f"downloading {url}", file=sys.stderr)
    urllib.request.urlretrieve(url, dest)
    return dest


def parse_series_matrix_meta(acc):
    """Return a per-sample DataFrame (index=GSM) from the series matrix."""
    url = f"{FTP}/{matrix_prefix(acc)}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    with urllib.request.urlopen(url, timeout=120) as r:
        raw = r.read()
    gsm = None
    fields = {}  # field_name -> list aligned to samples
    char_rows = []  # list of lists (each characteristics_ch line)
    with gzip.open(io.BytesIO(raw), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0].lstrip("!")
            vals = [p.strip().strip('"') for p in parts[1:]]
            if key == "Sample_geo_accession":
                gsm = vals
            elif key == "Sample_title":
                fields["title"] = vals
            elif key == "Sample_source_name_ch1":
                fields["source_name"] = vals
            elif key.startswith("Sample_characteristics_ch"):
                char_rows.append(vals)
    if gsm is None:
        raise RuntimeError(f"no GSM row for {acc}")
    df = pd.DataFrame({"gsm": gsm})
    for k, v in fields.items():
        if len(v) == len(gsm):
            df[k] = v
    # Explode characteristics "key: value" cells into columns
    parsed = {}
    for row in char_rows:
        if len(row) != len(gsm):
            continue
        for i, cell in enumerate(row):
            if ":" in cell:
                k, val = cell.split(":", 1)
                k = k.strip().lower()
                parsed.setdefault(k, [None] * len(gsm))
                parsed[k][i] = val.strip()
    for k, col in parsed.items():
        df[k] = col
    df.set_index("gsm", inplace=True)
    return df


def main():
    for acc, fname in EXPR_FILES.items():
        download(acc, fname)
        meta = parse_series_matrix_meta(acc)
        meta.to_csv(DATADIR / f"{acc}_sample_meta.csv")
        print(f"{acc}: meta shape {meta.shape}; cols={list(meta.columns)}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
