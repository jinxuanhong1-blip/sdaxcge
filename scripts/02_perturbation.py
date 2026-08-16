"""Claim A10 / P4: does experimental NKX2-1 loss raise the ELF3 module and TACSTD2/CLDN4?

Public processed matrices only. Sample groups are taken from GEO metadata or from
column names in the author-provided matrix — never inferred from the genes themselves.
"""

from __future__ import annotations

import os
import sys
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import requests
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

GEO = os.path.join(C.RAW, "geo")


def geo_sample_table(gse: str) -> pd.DataFrame:
    """GEO SOFT-lite sample characteristics via eutils esummary (GDS) is incomplete;
    use the series family miniml if present, else the NCBI gds ftp matrix header.
    Fallback: eutils esearch GSM + esummary.
    """
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
           f"?db=gds&retmax=1&term={gse}[ACCN]+AND+gse[ETYP]")
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    uid = ET.fromstring(r.text).findtext(".//Id")
    if not uid:
        return pd.DataFrame()
    # Fetch the series accession page samples via eSearch of GSMs
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
           f"?db=gds&retmax=400&term={gse}[ACCN]+AND+gsm[ETYP]")
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    ids = [e.text for e in ET.fromstring(r.text).findall(".//Id")]
    if not ids:
        return pd.DataFrame()
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
           f"?db=gds&id={','.join(ids)}")
    r = requests.get(url, timeout=90)
    r.raise_for_status()
    rows = []
    for d in ET.fromstring(r.text).findall("DocSum"):
        g = {i.get("Name"): (i.text or "") for i in d.findall("Item")}
        rows.append({"accession": g.get("Accession", ""),
                     "title": g.get("title", ""),
                     "summary": g.get("summary", "")})
    return pd.DataFrame(rows)


def pick_symbol_column(df: pd.DataFrame) -> str:
    for c in df.columns:
        cl = str(c).lower()
        if cl in {"gene_name", "gene", "symbol", "gene_symbol", "gene_short_name",
                  "hgnc_symbol", "mgi_symbol"}:
            return c
    # cufflinks
    if "gene_short_name" in df.columns:
        return "gene_short_name"
    return df.columns[0]


def collapse_symbol_matrix(df: pd.DataFrame, symbol_col: str, value_cols) -> pd.DataFrame:
    sub = df[[symbol_col] + list(value_cols)].copy()
    sub[symbol_col] = sub[symbol_col].astype(str).str.split(",").str[0].str.strip()
    sub = sub[sub[symbol_col].ne("") & sub[symbol_col].ne("nan")]
    for c in value_cols:
        sub[c] = pd.to_numeric(sub[c], errors="coerce")
    g = sub.groupby(symbol_col, sort=False).mean(numeric_only=True)
    return g


def logcpm_from_counts(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.divide(lib, axis=1) * 1e6 + 1.0)


def contrast(mat_log: pd.DataFrame, lo_cols, hi_cols, dataset, comparison, genes):
    """lo_cols = NKX-low / knockdown; hi_cols = control / NKX-high.

    Claim P4: genes in MODULE+TARGET should be HIGHER in lo than hi.
    """
    rows = []
    for g in genes:
        if g not in mat_log.index:
            rows.append({"dataset": dataset, "comparison": comparison, "gene": g,
                         "n_nkx_low": len(lo_cols), "n_nkx_high": len(hi_cols),
                         "detected": False, "log2fc_low_minus_high": np.nan,
                         "cliffs_delta": np.nan, "p_value": np.nan})
            continue
        a = mat_log.loc[g, lo_cols].to_numpy(dtype=float)
        b = mat_log.loc[g, hi_cols].to_numpy(dtype=float)
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if len(a) < 2 or len(b) < 2:
            p = np.nan
            d = np.nan
        else:
            try:
                _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            except ValueError:
                p = np.nan
            d = C.cliffs_delta(a, b)
        rows.append({"dataset": dataset, "comparison": comparison, "gene": g,
                     "n_nkx_low": len(a), "n_nkx_high": len(b), "detected": True,
                     "mean_low": float(np.nanmean(a)), "mean_high": float(np.nanmean(b)),
                     "log2fc_low_minus_high": float(np.nanmean(a) - np.nanmean(b)),
                     "cliffs_delta": d, "p_value": p})
    df = pd.DataFrame(rows)
    if len(df):
        df["fdr"] = C.bh_fdr(df["p_value"].to_numpy())
    return df


def gse129340():
    """H441 and H209 TTF-1/NKX2-1 siRNA. Cufflinks FPKM tracking files.

    These files are per-condition (already contrasted by the authors into one
    tracking table with q-values). We re-read FPKM columns if present; otherwise
    we take the author log2(fold_change) for the claim genes only and mark it
    as author-computed (not re-derived).
    """
    out = []
    notes = []
    for label, fname, genes in [
        ("GSE129340_H441_TTF1KD", "GSE129340_H441_TTF1KD_genes.fpkm_tracking.gz", C.CLAIM_GENES_H),
        ("GSE129340_H209_TTF1KD", "GSE129340_H209_TTF1KD_genes.fpkm_tracking.gz", C.CLAIM_GENES_H),
    ]:
        path = os.path.join(GEO, fname)
        if not os.path.exists(path):
            notes.append(f"MISSING {fname}")
            continue
        df = pd.read_csv(path, sep="\t")
        notes.append(f"{label} columns={list(df.columns)[:20]} n={len(df)}")
        # cufflinks tracking: gene_short_name, FPKM, plus optional sample FPKMs
        fpkm_cols = [c for c in df.columns if "FPKM" in c.upper() or c.endswith("_FPKM")]
        sym = "gene_short_name" if "gene_short_name" in df.columns else pick_symbol_column(df)
        # Prefer per-replicate FPKM columns
        rep = [c for c in df.columns if c not in
               {"tracking_id", "class_code", "nearest_ref_id", "gene_id",
                "gene_short_name", "tss_id", "locus", "length", "coverage"}
               and pd.api.types.is_numeric_dtype(df[c])]
        if len(rep) >= 4:
            mat = collapse_symbol_matrix(df, sym, rep)
            mat = np.log2(mat + 1.0)
            # Guess KD vs control from column names
            lo = [c for c in mat.columns if any(k in c.lower() for k in
                  ("kd", "sirna", "si_ttf", "sittf", "ttf1kd", "nkx", "knock"))]
            hi = [c for c in mat.columns if any(k in c.lower() for k in
                  ("ctrl", "control", "scr", "ntc", "sictrl", "mock"))]
            if not lo or not hi:
                # fall through to author fold-change
                notes.append(f"{label} could not split replicates lo={lo} hi={hi} cols={list(mat.columns)}")
            else:
                out.append(contrast(mat, lo, hi, label, "TTF1/NKX2-1 KD vs control", genes))
                continue
        # Author-provided fold change
        fc_col = None
        for c in df.columns:
            if c.lower() in {"log2(fold_change)", "log2_fold_change", "log2fc",
                             "log2(fold change)"} or "fold" in c.lower():
                fc_col = c
                break
        if fc_col is None:
            notes.append(f"{label} no usable FPKM replicates or fold-change column")
            continue
        sub = df[[sym, fc_col]].copy()
        sub[sym] = sub[sym].astype(str)
        sub[fc_col] = pd.to_numeric(sub[fc_col], errors="coerce")
        g = sub.groupby(sym)[fc_col].mean()
        rows = []
        for gene in genes:
            rows.append({"dataset": label, "comparison": "author log2FC (KD / control)",
                         "gene": gene, "detected": gene in g.index,
                         "log2fc_low_minus_high": float(g[gene]) if gene in g.index else np.nan,
                         "note": "author-computed fold-change; sign is KD minus control if they used KD/ctrl"})
        out.append(pd.DataFrame(rows))
    return (pd.concat(out, ignore_index=True) if out else pd.DataFrame()), notes


def gse229541():
    path = os.path.join(GEO, "GSE229541_NKX2-1_RNAseq_TPM_Final.csv.gz")
    notes = []
    if not os.path.exists(path):
        return pd.DataFrame(), ["MISSING GSE229541"]
    df = pd.read_csv(path)
    notes.append(f"GSE229541 shape={df.shape} cols={list(df.columns)[:25]}")
    sym = pick_symbol_column(df)
    num = [c for c in df.columns if c != sym and pd.api.types.is_numeric_dtype(df[c])]
    mat = collapse_symbol_matrix(df, sym, num)
    mat = np.log2(mat + 1.0)
    # Column names typically encode cell line + NKX2-1 high/low / sg / OE
    cols = list(mat.columns)
    notes.append("GSE229541 numeric columns: " + ",".join(map(str, cols)))
    # Group by tokens
    low_tok = ("low", "kd", "ko", "sg", "shrna", "crispra_off", "loss", "deplete")
    high_tok = ("high", "oe", "over", "wt", "ctrl", "control", "parental", "vec")
    # We will emit every pairwise contrast we can name from the columns, plus
    # a generic low-vs-high if tokens match.
    lo = [c for c in cols if any(t in str(c).lower() for t in low_tok)
          and not any(t in str(c).lower() for t in ("high", "oe"))]
    hi = [c for c in cols if any(t in str(c).lower() for t in ("high", "oe", "wt", "ctrl", "control", "parental"))]
    frames = []
    if lo and hi:
        frames.append(contrast(mat, lo, hi, "GSE229541",
                               f"NKX2-1-low/KD ({len(lo)}) vs high/ctrl ({len(hi)})",
                               C.CLAIM_GENES_H))
    # Also dump per-column mean of claim genes so the report can be honest if grouping failed
    have = [g for g in C.CLAIM_GENES_H if g in mat.index]
    dump = mat.loc[have].T.reset_index().rename(columns={"index": "sample"})
    dump.insert(0, "dataset", "GSE229541")
    dump.to_csv(os.path.join(C.TABLES, "gse229541_claim_genes_log2tpm1p.csv"), index=False)
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), notes


def gse129583():
    path = os.path.join(GEO, "GSE129583_processeddata_matrix_Nkx2-1_RNAseq.xlsx")
    notes = []
    if not os.path.exists(path):
        return pd.DataFrame(), ["MISSING GSE129583"]
    xl = pd.ExcelFile(path)
    notes.append(f"GSE129583 sheets={xl.sheet_names}")
    frames = []
    for sheet in xl.sheet_names:
        df = xl.parse(sheet)
        notes.append(f"sheet {sheet} shape={df.shape} cols={list(df.columns)[:20]}")
        # first column gene
        df = df.rename(columns={df.columns[0]: "gene"})
        df["gene"] = df["gene"].astype(str)
        # mouse symbols often already title-case or uppercase
        num = [c for c in df.columns[1:] if pd.api.types.is_numeric_dtype(df[c])]
        if len(num) < 2:
            continue
        mat = df.set_index("gene")[num]
        # map to mouse claim genes with case-insensitive match
        idx_map = {str(i).lower(): i for i in mat.index}
        want = {}
        for g in C.CLAIM_GENES_M + list(C.MARKERS_M):
            if g.lower() in idx_map:
                want[g] = idx_map[g.lower()]
            elif g.upper() in mat.index:
                want[g] = g.upper()
        mat2 = mat.loc[list(want.values())]
        mat2.index = list(want.keys())
        mat2 = np.log2(mat2.clip(lower=0) + 1.0)
        lo = [c for c in mat2.columns if any(t in str(c).lower() for t in
              ("ko", "cko", "mut", "nkx", "deleted", "cnull", "null", "cre"))]
        hi = [c for c in mat2.columns if any(t in str(c).lower() for t in
              ("wt", "ctrl", "control", "fl", "flox", "het"))]
        # avoid putting the same column in both
        lo = [c for c in lo if c not in hi]
        notes.append(f"sheet {sheet} lo={lo} hi={hi}")
        if lo and hi:
            frames.append(contrast(mat2, lo, hi, f"GSE129583:{sheet}",
                                   "Nkx2-1 loss vs control",
                                   C.CLAIM_GENES_M + list(C.MARKERS_M)))
        dump = mat2.T.reset_index().rename(columns={"index": "sample"})
        dump.insert(0, "dataset", f"GSE129583:{sheet}")
        safe = sheet.replace(" ", "_")[:40]
        dump.to_csv(os.path.join(C.TABLES, f"gse129583_{safe}_claim_genes.csv"), index=False)
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), notes


def mouse_luad_counts(fname, dataset, genes):
    path = os.path.join(GEO, fname)
    notes = []
    if not os.path.exists(path):
        return pd.DataFrame(), [f"MISSING {fname}"]
    df = pd.read_csv(path, sep="\t")
    notes.append(f"{dataset} shape={df.shape} cols={list(df.columns)[:15]}")
    sym = pick_symbol_column(df)
    num = [c for c in df.columns if c != sym and pd.api.types.is_numeric_dtype(df[c])]
    mat = collapse_symbol_matrix(df, sym, num)
    # case-insensitive gene match
    idx = {str(i).lower(): i for i in mat.index}
    keep = {g: idx[g.lower()] for g in genes if g.lower() in idx}
    mat = mat.loc[list(keep.values())]
    mat.index = list(keep.keys())
    log = logcpm_from_counts(mat)
    cols = list(log.columns)
    notes.append(f"{dataset} samples={cols}")
    lo = [c for c in cols if any(t in str(c).lower() for t in
          ("neg", "low", "ko", "null", "loss", "nkx2-1-", "nkx21-"))]
    hi = [c for c in cols if any(t in str(c).lower() for t in
          ("pos", "high", "wt", "ctrl", "nkx2-1+", "nkx21+"))]
    frames = []
    if lo and hi:
        frames.append(contrast(log, lo, hi, dataset, "Nkx2-1-low/neg vs high/pos", genes))
    dump = log.T.reset_index().rename(columns={"index": "sample"})
    dump.insert(0, "dataset", dataset)
    dump.to_csv(os.path.join(C.TABLES, f"{dataset.lower()}_claim_genes_logcpm.csv"), index=False)
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), notes


def main():
    all_notes = []
    frames = []
    for fn, notes in (
        gse129340(),
        gse229541(),
        gse129583(),
        mouse_luad_counts("GSE115899_counts_13438R.txt.gz", "GSE115899",
                          C.CLAIM_GENES_M + list(C.MARKERS_M)),
        mouse_luad_counts("GSE145152_counts_14489R.txt.gz", "GSE145152",
                          C.CLAIM_GENES_M + list(C.MARKERS_M)),
        mouse_luad_counts("GSE188435_17699R_17992R_counts.txt.gz", "GSE188435",
                          C.CLAIM_GENES_M + list(C.MARKERS_M)),
    ):
        if len(fn):
            frames.append(fn)
        all_notes.extend(notes)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if len(out):
        C.write_table(out, "perturbation_contrasts.csv")
    with open(os.path.join(C.LOGS, "perturbation_notes.txt"), "w") as fh:
        fh.write("\n".join(all_notes) + "\n")
    print("\n".join(all_notes))
    print("done", flush=True)


if __name__ == "__main__":
    main()
