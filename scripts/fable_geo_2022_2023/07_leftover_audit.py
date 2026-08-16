#!/usr/bin/env python3
"""
Finish leftover 2022-2023 GEO lung ICI series that were not deeply verified
in the first pass.

For each leftover / borderline series:
  * fetch series-matrix metadata (if missing)
  * download open processed expression < 2 GB when it might contain
    TACSTD2/CLDN4 + a per-patient ICI outcome
  * record gene presence, sample mapping, and a written verdict

No fabricated statistics: tests are only run when both genes (or a gene)
and a real per-sample outcome label exist.
"""
import gzip
import io
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "fable_geo_2022_2023"
DATA = OUT / "data"
META = OUT / "series_matrix_meta"
DATA.mkdir(parents=True, exist_ok=True)
META.mkdir(parents=True, exist_ok=True)
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# leftover relevant (empty verdict) + borderline series that look like
# lung ICI outcome studies but failed the first-pass text filter
LEFTOVER = [
    "GSE185204", "GSE186446", "GSE248378", "GSE235048", "GSE228419",
    "GSE164146", "GSE212622", "GSE224099", "GSE193719", "GSE189804",
    "GSE217451", "GSE224216", "GSE150255", "GSE194350", "GSE195770",
    "GSE218402", "GSE224246", "GSE238006", "GSE250254", "GSE193049",
    "GSE229353", "GSE178521", "GSE192591", "GSE192790", "GSE197236",
    "GSE198099", "GSE213590", "GSE213902", "GSE223779", "GSE193707",
    # borderline / possibly missed
    "GSE221733", "GSE221322", "GSE189045", "GSE250262",
]

# processed files worth downloading (<2GB, might have TACSTD2/CLDN4)
DOWNLOAD = {
    "GSE248378": ["GSE248378_Durva_Post_FPKMs.txt.gz"],
    "GSE235048": ["GSE235048_24CFlow_TPM.txt.gz"],
    "GSE186446": ["GSE186446_fcount_aggr.txt.gz"],
    "GSE228419": ["GSE228419_normalized_counts.csv.gz"],
    "GSE193049": [],  # RAW.tar only; inspect after listing
    "GSE221733": ["GSE221733_4301_CTA_initial.csv.gz",
                  "GSE221733_4301_CTA_norm.xlsx",
                  "GSE221733_4301_CTA_QC.xlsx"],
    "GSE221322": ["GSE221322_4301_protein_norm.csv.gz"],
}


def prefix(acc):
    num = re.match(r"GSE(\d+)", acc).group(1)
    return "GSE" + (num[:-3] + "nnn" if len(num) > 3 else "nnn")


def fetch(url, dest, timeout=120):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  exists {dest.name}", file=sys.stderr)
        return dest
    print(f"  GET {url}", file=sys.stderr)
    urllib.request.urlretrieve(url, dest)
    return dest


def list_dir(url):
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            html = r.read().decode("utf-8", "replace")
    except Exception as e:  # noqa: BLE001
        return [], str(e)
    names = re.findall(r'href="([^"?/][^"]*)"', html)
    return [n for n in names if n not in ("../",) and not n.endswith("/")], ""


def parse_series_matrix_meta(acc):
    url = f"{FTP}/{prefix(acc)}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            raw = r.read()
    except Exception as e:  # noqa: BLE001
        return None, str(e)
    gsm = None
    fields = {}
    char_rows = []
    extra = {}
    with gzip.open(io.BytesIO(raw), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_"):
                parts = line.rstrip("\n").split("\t", 1)
                k = parts[0].lstrip("!")
                v = parts[1].strip().strip('"') if len(parts) > 1 else ""
                extra.setdefault(k, []).append(v)
                continue
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0].lstrip("!")
            vals = [p.strip().strip('"') for p in parts[1:]]
            if key == "Sample_geo_accession":
                gsm = vals
            elif key in ("Sample_title", "Sample_source_name_ch1",
                         "Sample_description"):
                fields[key] = vals
            elif key.startswith("Sample_characteristics_ch"):
                char_rows.append(vals)
    if gsm is None:
        return None, "no GSM row"
    df = pd.DataFrame({"gsm": gsm})
    for k, v in fields.items():
        if len(v) == len(gsm):
            df[k] = v
    parsed = {}
    for row in char_rows:
        if len(row) != len(gsm):
            continue
        for i, cell in enumerate(row):
            if ":" in cell:
                k, val = cell.split(":", 1)
                parsed.setdefault(k.strip().lower(), [None] * len(gsm))
                parsed[k.strip().lower()][i] = val.strip()
    for k, col in parsed.items():
        df[k] = col
    df.set_index("gsm", inplace=True)
    return (df, extra), ""


def peek_genes(path, n_header=1):
    """Return (n_rows, n_cols, first_cols, gene_hits) for a matrix file."""
    p = Path(path)
    name = p.name.lower()
    if name.endswith(".xlsx"):
        try:
            x = pd.read_excel(p)
        except Exception as e:  # noqa: BLE001
            return {"error": str(e)}
        genes = set()
        for col in x.columns:
            if str(col).upper() in ("TACSTD2", "CLDN4", "TROP2", "CLAUDIN-4"):
                genes.add(str(col))
        if "Gene" in x.columns or "gene" in x.columns:
            gcol = "Gene" if "Gene" in x.columns else "gene"
            genes |= set(x[gcol].astype(str).str.upper()) & {"TACSTD2", "CLDN4"}
        # also search first column values
        first = x.iloc[:, 0].astype(str).str.upper()
        genes |= set(first) & {"TACSTD2", "CLDN4"}
        return {"kind": "xlsx", "shape": list(x.shape),
                "cols": list(x.columns)[:12],
                "genes": sorted(genes),
                "head": x.head(3).to_dict(orient="list")}
    # text
    try:
        if str(p).endswith(".gz"):
            fh = gzip.open(p, "rt", errors="replace")
        else:
            fh = open(p, "rt", errors="replace")
        with fh:
            header = fh.readline()
            # detect sep
            sep = "\t" if header.count("\t") >= header.count(",") else ","
            if header.count(";") > header.count(sep):
                sep = ";"
            cols = [c.strip().strip('"') for c in header.split(sep)]
            genes = []
            n = 0
            for line in fh:
                n += 1
                tok = line.split(sep, 1)[0].strip().strip('"').upper()
                if tok in ("TACSTD2", "CLDN4", "TROP2"):
                    genes.append(tok)
            return {"kind": "text", "n_genes": n, "n_cols": len(cols) - 1,
                    "first_cols": cols[:8], "genes": genes, "sep": sep}
    except Exception as e:  # noqa: BLE001
        return {"error": str(e)}


def main():
    ver = pd.read_csv(OUT / "geo_verified.csv")
    rows = []
    for acc in LEFTOVER:
        print(f"\n==== {acc} ====", file=sys.stderr)
        rec = {"accession": acc}
        vr = ver[ver.accession == acc]
        if not vr.empty:
            r = vr.iloc[0]
            rec.update({
                "title": r["title"], "n_samples": r["n_samples"],
                "pdat": r["pdat"], "gdstype": r["gdstype"],
                "relevant": bool(r["relevant"]),
                "suppl_files": r["suppl_files"],
                "suppl_total_mb": r["suppl_total_mb"],
            })
        # series matrix
        parsed, err = parse_series_matrix_meta(acc)
        if parsed is None:
            rec["meta_error"] = err
            print(f"  meta error: {err}", file=sys.stderr)
        else:
            df, extra = parsed
            df.to_csv(DATA / f"{acc}_sample_meta.csv")
            rec["n_meta"] = len(df)
            rec["meta_cols"] = ";".join(df.columns)
            rec["series_summary"] = " ".join(extra.get("Series_summary", []))[:400]
            rec["series_overall_design"] = " ".join(
                extra.get("Series_overall_design", []))[:400]
            # write json if missing
            jp = META / f"{acc}.json"
            if not jp.exists():
                chars = {c: sorted(df[c].dropna().astype(str).unique())[:25]
                         for c in df.columns if c not in
                         ("Sample_title", "Sample_source_name_ch1",
                          "Sample_description", "title", "source_name")}
                jp.write_text(json.dumps({
                    "accession": acc, "n_samples": len(df),
                    "characteristics": chars,
                    "title": rec.get("title", ""),
                }, indent=2, ensure_ascii=False))
            print(f"  meta n={len(df)} cols={list(df.columns)}", file=sys.stderr)
        time.sleep(0.15)

        # downloads
        gene_info = {}
        for fname in DOWNLOAD.get(acc, []):
            url = f"{FTP}/{prefix(acc)}/{acc}/suppl/{fname}"
            dest = DATA / f"{acc}__{fname}"
            try:
                fetch(url, dest)
            except Exception as e:  # noqa: BLE001
                gene_info[fname] = {"error": str(e)}
                continue
            gene_info[fname] = peek_genes(dest)
            print(f"  peek {fname}: {gene_info[fname]}", file=sys.stderr)
        rec["gene_info"] = json.dumps(gene_info, ensure_ascii=False)
        rows.append(rec)

    pd.DataFrame(rows).to_csv(OUT / "leftover_audit.csv", index=False)
    print(f"\nWrote leftover_audit.csv ({len(rows)} series)", file=sys.stderr)


if __name__ == "__main__":
    main()
