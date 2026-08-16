#!/usr/bin/env python3
"""Fetch and harmonise public ICI cohorts with PFS from GEO.

Cohorts (both public, no authentication required):
  * GSE135222  advanced NSCLC, anti-PD-1/PD-L1, RNA-seq (TPM), PFS in DAYS
  * GSE190265  advanced NSCLC, anti-PD-1, RNA-seq (TPM), PFS in MONTHS

Outputs (small, committed so that the analysis is reproducible offline):
  data/<gse>_clinical.csv      one row per sample: sample_id, time_months, event
  data/<gse>_genes_tpm.csv     TPM restricted to the illustrative gene panel
  data/provenance.json         URLs, sizes, sha256 of every raw file used

Raw downloads land in data/raw/ and are NOT committed (see .gitignore).

Usage:
  python3 01_fetch_geo.py                 # download (or reuse cache) and harmonise
  python3 01_fetch_geo.py --verify-symbols  # re-resolve symbols via Ensembl REST
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import sys
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RAW = DATA / "raw"

# Primary genes of interest for this project, plus CD8A as a known-direction
# immune comparator (not a TACSTD2/CLDN4 claim). None of these is a validated
# ICI biomarker on these cohorts.
PANEL = ["TACSTD2", "CLDN4", "CD8A"]

# Primary-assembly Ensembl gene IDs, resolved via the Ensembl REST xrefs endpoint
# (rest.ensembl.org/xrefs/symbol/homo_sapiens/<symbol>?object_type=gene).
# Re-check with --verify-symbols.
ENSEMBL = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
    "CD8A": "ENSG00000153563",
}

FILES = {
    "GSE135222_series_matrix": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/"
        "GSE135222_series_matrix.txt.gz"
    ),
    "GSE135222_tpm": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/"
        "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
    ),
    "GSE190265_clinical": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/suppl/"
        "GSE190265_samples_info_France3.csv.gz"
    ),
    "GSE190265_tpm": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/suppl/"
        "GSE190265_TPM_France3.csv.gz"
    ),
}

DAYS_PER_MONTH = 30.4375  # 365.25 / 12


def download(key: str, url: str) -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    dest = RAW / url.rsplit("/", 1)[-1]
    if not dest.exists():
        print(f"[fetch] {key}: {url}", file=sys.stderr)
        with urllib.request.urlopen(url, timeout=120) as fh, dest.open("wb") as out:
            out.write(fh.read())
    return dest


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_series_matrix(path: Path) -> pd.DataFrame:
    """Return the !Sample_* header block of a GEO series matrix as a tidy frame."""
    rows: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            key, values = parts[0], parts[1:]
            # characteristics_ch1 repeats; disambiguate by its "field: value" prefix
            if key == "!Sample_characteristics_ch1" and values and ":" in values[0]:
                key = f"char::{values[0].split(':', 1)[0].strip()}"
                values = [v.split(":", 1)[1].strip() if ":" in v else v for v in values]
            suffix = 2
            base = key
            while key in rows:
                key = f"{base}#{suffix}"
                suffix += 1
            rows[key] = values
    return pd.DataFrame(rows)


def clinical_gse135222(path: Path) -> pd.DataFrame:
    meta = parse_series_matrix(path)
    time_col = next(c for c in meta.columns if "pfs.time" in c.lower())
    # The event indicator is the characteristic literally named
    # "progression-free survival (pfs)"; 1 = progression/death, 0 = censored.
    event_col = next(
        c for c in meta.columns if c.startswith("char::") and "progression-free" in c.lower()
    )
    out = pd.DataFrame(
        {
            "sample_id": meta["!Sample_title"].str.replace(" ", "", regex=False),
            "geo_accession": meta["!Sample_geo_accession"],
            "time_days": pd.to_numeric(meta[time_col]),
            "event": pd.to_numeric(meta[event_col]).astype(int),
            "age": pd.to_numeric(meta[next(c for c in meta.columns if c.endswith("age"))]),
            "sex": meta[next(c for c in meta.columns if c.endswith("gender"))],
        }
    )
    out["time_months"] = out["time_days"] / DAYS_PER_MONTH
    return out


def expression_gse135222(path: Path) -> pd.DataFrame:
    """TPM matrix: rows = versioned Ensembl IDs, cols = samples. Return samples x panel."""
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.str.split(".").str[0]
    wanted = {ENSEMBL[g]: g for g in PANEL}
    missing = [g for eid, g in wanted.items() if eid not in expr.index]
    if missing:
        raise SystemExit(f"GSE135222: panel genes absent from matrix: {missing}")
    sub = expr.loc[list(wanted)].rename(index=wanted).T
    sub.index.name = "sample_id"
    return sub[PANEL]


def clinical_gse190265(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=";")
    out = pd.DataFrame(
        {
            "sample_id": df["sample"].astype(str),
            "time_months": pd.to_numeric(df["time_PFS"]),
            "event": pd.to_numeric(df["evtPFS"]).astype(int),
        }
    )
    out["time_days"] = out["time_months"] * DAYS_PER_MONTH
    return out


def expression_gse190265(path: Path) -> pd.DataFrame:
    """TPM matrix: rows = samples, cols = gene symbols; the header lacks the id column."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split(";")
        body = fh.read()
    expr = pd.read_csv(
        io.StringIO(body), sep=";", header=None, names=["sample_id"] + header, index_col=0
    )
    expr.index = expr.index.astype(str)
    # The deposited header has been through R's make.names(): HLA-DRA became HLA.DRA.
    resolved = {}
    for gene in PANEL:
        for candidate in (gene, gene.replace("-", ".")):
            if candidate in expr.columns:
                resolved[gene] = candidate
                break
    missing = [g for g in PANEL if g not in resolved]
    if missing:
        raise SystemExit(f"GSE190265: panel genes absent from matrix: {missing}")
    sub = expr[[resolved[g] for g in PANEL]].copy()
    sub.columns = PANEL
    sub.index.name = "sample_id"
    return sub


def verify_symbols() -> None:
    for symbol, expected in ENSEMBL.items():
        url = (
            "https://rest.ensembl.org/xrefs/symbol/homo_sapiens/"
            f"{symbol}?object_type=gene;content-type=application/json"
        )
        with urllib.request.urlopen(url, timeout=60) as fh:
            ids = [x["id"] for x in json.load(fh) if x["id"].startswith("ENSG")]
        status = "ok" if expected in ids else "MISMATCH"
        print(f"{symbol:8s} expected={expected} ensembl={ids} {status}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify-symbols", action="store_true", help="re-resolve panel via Ensembl")
    args = ap.parse_args()
    if args.verify_symbols:
        verify_symbols()
        return

    DATA.mkdir(parents=True, exist_ok=True)
    paths = {key: download(key, url) for key, url in FILES.items()}

    clin222 = clinical_gse135222(paths["GSE135222_series_matrix"])
    expr222 = expression_gse135222(paths["GSE135222_tpm"])
    shared = clin222["sample_id"].isin(expr222.index)
    if not shared.all():
        raise SystemExit(f"GSE135222 id mismatch: {clin222.loc[~shared, 'sample_id'].tolist()}")

    clin265 = clinical_gse190265(paths["GSE190265_clinical"])
    expr265 = expression_gse190265(paths["GSE190265_tpm"])
    shared = clin265["sample_id"].isin(expr265.index)
    if not shared.all():
        raise SystemExit(f"GSE190265 id mismatch: {clin265.loc[~shared, 'sample_id'].tolist()}")

    clin222.to_csv(DATA / "gse135222_clinical.csv", index=False)
    expr222.round(4).to_csv(DATA / "gse135222_genes_tpm.csv")
    clin265.to_csv(DATA / "gse190265_clinical.csv", index=False)
    expr265.round(4).to_csv(DATA / "gse190265_genes_tpm.csv")

    provenance = {
        "retrieved": date.today().isoformat(),
        "panel": PANEL,
        "panel_ensembl_ids": ENSEMBL,
        "time_unit_note": {
            "GSE135222": "pfs.time is in days; converted with 365.25/12 days per month",
            "GSE190265": "time_PFS is already in months as deposited",
        },
        "files": {
            key: {
                "url": FILES[key],
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for key, path in paths.items()
        },
        "cohorts": {
            "GSE135222": {"n": int(len(clin222)), "events": int(clin222["event"].sum())},
            "GSE190265": {"n": int(len(clin265)), "events": int(clin265["event"].sum())},
        },
    }
    (DATA / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    for name, clin in (("GSE135222", clin222), ("GSE190265", clin265)):
        print(
            f"{name}: n={len(clin)} events={int(clin['event'].sum())} "
            f"censored={int((clin['event'] == 0).sum())} "
            f"max_followup_months={clin['time_months'].max():.1f}"
        )


if __name__ == "__main__":
    main()
