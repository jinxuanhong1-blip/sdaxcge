#!/usr/bin/env python3
"""
Step 8: for every non-embedded series that exposes a candidate processed
gene-level table, download it (git-ignored), locate TACSTD2 / CLDN4 (by gene
symbol or Ensembl gene id, mouse or human) and record the marker row.

Two table shapes are handled:
  * differential-expression result tables -> capture logFC / p / padj columns
  * count / TPM / FPKM matrices           -> capture per-sample values

Writes:
  notes/geo_sweep/analysis_supp.json          (per-series marker rows + meta)
  results/geo_sweep/marker_rows_supp.tsv      (flat marker row dump)
"""
import gzip
import io
import json
import re
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

NOTES = Path("notes/geo_sweep")
RESULTS = Path("results/geo_sweep")
SUPP = RESULTS / "supp"
SUPP.mkdir(parents=True, exist_ok=True)

MAX_BYTES = 300 * 1024 ** 2

SYM = {
    "TACSTD2": {"TACSTD2", "TROP2", "TROP-2", "GA733-1", "M1S1", "EGP-1", "EGP1"},
    "CLDN4": {"CLDN4", "CLAUDIN-4", "CLAUDIN4", "CPETR1", "CPE-R", "CPER"},
}
# Verified via Ensembl REST lookup/symbol (2026-08-16).
# Mouse Cldn4 is ENSMUSG00000047501 (NOT ENSMUSG00000024959, which is Bad).
ENS = {
    "TACSTD2": {"ENSG00000184292", "ENSMUSG00000051397"},
    "CLDN4": {"ENSG00000189143", "ENSMUSG00000047501"},
}
ENS_TX = {
    "TACSTD2": {"ENST00000371225", "ENSMUST00000058178"},
    "CLDN4": {"ENST00000340958", "ENST00000466411", "ENST00000435050",
              "ENST00000476494", "ENST00000431918", "ENSMUST00000051401"},
}
ENTREZ = {
    "TACSTD2": {"4070", "56753"},   # human, mouse
    "CLDN4": {"1364", "12740"},
}
DE_COLS = re.compile(
    r"log2?[\s_]?fold|log2fc|logfc|fold[\s_]?change|p[\s_.]?val|pvalue|"
    r"padj|adj[\s_.]?p|fdr|q[\s_.]?value|baseMean|stat|t\.value", re.IGNORECASE)


def http_download(url, dest, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-sweep/1.0"})
            with urllib.request.urlopen(req, timeout=240) as r, open(dest, "wb") as f:
                f.write(r.read())
            return dest.stat().st_size
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"download {url}: {last}")


def read_table(path):
    name = path.name.lower()
    if name.endswith(".xlsx"):
        return pd.read_excel(path, engine="openpyxl")
    raw = path.read_bytes()
    if name.endswith(".gz"):
        raw = gzip.decompress(raw)
    text = raw.decode("utf-8", "replace")
    sample = text[:5000]
    if sample.count("\t") >= sample.count(",") and "\t" in sample:
        sep = "\t"
    elif ";" in sample and sample.count(";") > sample.count(","):
        sep = ";"
    else:
        sep = ","
    df = pd.read_csv(io.StringIO(text), sep=sep, low_memory=False)
    # gene ids sometimes land in the index (header has one fewer field than data)
    if not isinstance(df.index, pd.RangeIndex):
        df = df.reset_index()
    return df


ID_COL_RE = re.compile(r"gene[\s_.]?id|entrez|geneid|ensembl|symbol|gene[\s_.]?name|^gene$|^id$",
                       re.IGNORECASE)


def _row_dict(df, i):
    out = {}
    for c in df.columns:
        v = df.at[i, c]
        if pd.isna(v):
            out[c] = None
        elif isinstance(v, (int, float, np.integer, np.floating)):
            out[c] = float(v)
        else:
            out[c] = str(v)
    return out


def find_marker_rows(df):
    """Return dict gene -> list of matched row dicts (id col + full row).

    Matches gene symbol, Ensembl gene id (version-stripped), or Entrez id.
    Identifier columns are any non-numeric column, plus the first column and any
    column named like a gene id (covers integer Entrez id columns).
    """
    id_cols = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
    id_cols.append(df.columns[0])
    for c in df.columns:
        if ID_COL_RE.search(str(c)):
            id_cols.append(c)
    id_cols = list(dict.fromkeys(id_cols))  # unique, keep order

    hits = {}
    for gene in SYM:
        syms, enss, entz = SYM[gene], ENS[gene], ENTREZ[gene]
        txs = ENS_TX[gene]
        matched_idx = set()
        for c in id_cols:
            up = df[c].astype(str).str.upper().str.strip().str.strip('"')
            base = up.str.replace(r"\.\d+$", "", regex=True)
            # entrez columns may be "497097.0" if float-coerced
            ent = base.str.replace(r"\.0$", "", regex=True)
            mask = (up.isin(syms) | base.isin(enss) | base.isin(txs)
                    | ent.isin(entz))
            for i in df.index[mask]:
                matched_idx.add(i)
        if matched_idx:
            hits[gene] = [_row_dict(df, i) for i in sorted(matched_idx)]
    return hits


def find_markers_any_orientation(df):
    """Try normal orientation; if nothing, try transposed (genes as columns)."""
    hits = find_marker_rows(df)
    if hits:
        return hits, False
    # transposed: genes are column headers, samples are rows
    up_cols = {str(c).upper().strip().strip('"'): c for c in df.columns}
    all_syms = SYM["TACSTD2"] | SYM["CLDN4"]
    if any(s in up_cols for s in all_syms):
        first = df.columns[0]
        dft = df.set_index(first).T
        dft.index = dft.index.astype(str)
        dft = dft.reset_index().rename(columns={"index": "gene"})
        for c in dft.columns:
            if c != "gene":
                dft[c] = pd.to_numeric(dft[c], errors="coerce")
        return find_marker_rows(dft), True
    return {}, False


def numeric_sample_cols(df):
    return [c for c in df.columns
            if pd.api.types.is_numeric_dtype(df[c]) and not DE_COLS.search(str(c))]


def main():
    listing = json.load(open(NOTES / "supp_listing.json"))
    results = []
    flat = []
    for o in listing:
        acc = o["accession"]
        cands = sorted(o["candidates"], key=lambda f: (f["size_bytes"] or 1e18))
        if not cands:
            results.append({"accession": acc, "status": "no_candidate_file",
                            "is_ici": o["is_ici"], "is_kd_ko": o["is_kd_ko"]})
            continue
        series_res = {"accession": acc, "is_ici": o["is_ici"],
                      "is_kd_ko": o["is_kd_ko"], "files_tried": [],
                      "markers": {}}
        found_any = False
        for f in cands:
            if (f["size_bytes"] or 0) > MAX_BYTES:
                series_res["files_tried"].append(
                    {"name": f["name"], "status": "too_large"})
                continue
            dest = SUPP / f["name"]
            try:
                if not dest.exists():
                    http_download(f["url"], dest)
                df = read_table(dest)
            except Exception as e:  # noqa: BLE001
                series_res["files_tried"].append(
                    {"name": f["name"], "status": f"read_error: {e}"})
                continue
            hits, transposed = find_markers_any_orientation(df)
            is_de = bool([c for c in df.columns if DE_COLS.search(str(c))]) and not transposed
            scols = numeric_sample_cols(df) if not transposed else []
            ft = {"name": f["name"], "shape": list(df.shape),
                  "is_de_table": is_de, "transposed": transposed,
                  "n_numeric_sample_cols": len(scols),
                  "markers_found": sorted(hits.keys())}
            series_res["files_tried"].append(ft)
            for gene, rows in hits.items():
                if gene in series_res["markers"]:
                    continue
                rec = {"file": f["name"], "is_de_table": is_de,
                       "transposed": transposed,
                       "n_rows_matched": len(rows), "rows": rows[:5]}
                if transposed and not scols:
                    scols = [k for k, v in rows[0].items()
                             if k != "gene" and isinstance(v, (int, float))]
                if not is_de and scols:
                    vals = []
                    for row in rows:
                        for c in scols:
                            v = row.get(c)
                            if isinstance(v, (int, float)) and v is not None and not pd.isna(v):
                                vals.append(float(v))
                    if vals:
                        arr = np.array(vals, dtype=float)
                        rec["per_sample_summary"] = {
                            "n_samples": len(scols), "n_values": int(arr.size),
                            "mean": float(np.mean(arr)),
                            "median": float(np.median(arr)),
                            "min": float(np.min(arr)), "max": float(np.max(arr)),
                        }
                        rec["sample_columns"] = scols[:60]
                        rec["per_sample_values"] = {
                            c: (float(rows[0][c]) if isinstance(rows[0].get(c), (int, float))
                                and rows[0].get(c) is not None else None)
                            for c in scols}
                series_res["markers"][gene] = rec
                found_any = True
                flat.append({"accession": acc, "gene": gene, "file": f["name"],
                             "is_de_table": is_de,
                             "n_rows_matched": len(rows),
                             "first_row": json.dumps(rows[0])[:800]})
            if len(series_res["markers"]) >= 2:
                break
        series_res["status"] = "markers_found" if found_any else "markers_absent"
        results.append(series_res)
        print(f"{acc:12s} ici={o['is_ici']} kd={o['is_kd_ko']} "
              f"-> {series_res['status']} "
              f"({','.join(series_res['markers'].keys()) or '-'})")

    (NOTES / "analysis_supp.json").write_text(json.dumps(results, indent=2))
    if flat:
        pd.DataFrame(flat).to_csv(RESULTS / "marker_rows_supp.tsv",
                                  sep="\t", index=False)
    n_found = sum(1 for r in results if r.get("status") == "markers_found")
    print(f"\nMarkers found in {n_found} series. "
          f"Wrote analysis_supp.json, marker_rows_supp.tsv")


if __name__ == "__main__":
    main()
