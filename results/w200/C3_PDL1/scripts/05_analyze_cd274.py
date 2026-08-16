#!/usr/bin/env python3
"""Quantify CD274/Cd274 after CLDN4 loss or TROP2-ADC / SN-38 in public matrices.

Every contrast is labeled with what it actually is. Association and lookalikes
are not treated as C3.

Outputs under omics/cd274_contrasts/
"""

from __future__ import annotations

import gzip
import io
import json
import os
import sys
import urllib.request
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "omics", "cd274_contrasts")
UA = "C3-PDL1-catalog/1.0"
os.makedirs(DATA, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

GENES_H = {"CD274", "CLDN4", "TACSTD2", "PDCD1LG2", "CD8A", "STAT1", "IRF1"}
GENES_M = {"Cd274", "Cldn4", "Tacstd2", "Pdcd1lg2", "Cd8a", "Stat1", "Irf1"}


def ftp(acc: str, rel: str) -> str:
    n = int(acc[3:])
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{n // 1000}nnn/{acc}/{rel}"


def download(url: str, dest: str) -> str:
    if os.path.exists(dest) and os.path.getsize(dest) > 100:
        return dest
    print(f"GET {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=180) as fh, open(dest, "wb") as out:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    return dest


def log2p(x):
    x = np.asarray(x, dtype=float)
    return np.log2(x + 1.0)


def summarize(values_a, values_b, name_a, name_b):
    a = np.asarray(values_a, dtype=float)
    b = np.asarray(values_b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return None
    # Welch t on the supplied scale (already log-ish if we pass log2p)
    t, p = stats.ttest_ind(b, a, equal_var=False)
    lfc = float(np.mean(b) - np.mean(a))
    return {
        "n_ctrl": int(len(a)),
        "n_treat": int(len(b)),
        "mean_ctrl": float(np.mean(a)),
        "mean_treat": float(np.mean(b)),
        "log2FC_treat_minus_ctrl": lfc,
        "welch_p": float(p) if np.isfinite(p) else None,
        "direction": "UP" if lfc > 0 else ("DOWN" if lfc < 0 else "FLAT"),
        "ctrl_name": name_a,
        "treat_name": name_b,
        "ctrl_values": [float(x) for x in a],
        "treat_values": [float(x) for x in b],
    }


def pick_gene_row(df, candidates):
    """Find a gene row. Tries index and a few annotation columns."""
    idx = df.index.astype(str)
    for c in candidates:
        hit = idx[idx.str.upper() == c.upper()]
        if len(hit):
            return hit[0]
    # symbol in a column
    for col in df.columns[:4]:
        s = df[col].astype(str).str.upper()
        for c in candidates:
            hit = df.index[s == c.upper()]
            if len(hit):
                return hit[0]
    # gene_id|symbol
    for c in candidates:
        hit = idx[idx.str.upper().str.contains(rf"(^|[|; ]){c.upper()}($|[|; ])")]
        if len(hit):
            return hit[0]
    return None


def write_contrast(rows, path):
    cols = [
        "accession", "model", "contrast", "qualifies_as", "gene", "scale",
        "n_ctrl", "n_treat", "mean_ctrl", "mean_treat",
        "log2FC_treat_minus_ctrl", "welch_p", "direction", "note",
    ]
    pd.DataFrame(rows)[cols].to_csv(path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# GSE207704  MCF7 / T47D CLDN4-/- vs parental
# ---------------------------------------------------------------------------
def gse207704():
    acc = "GSE207704"
    dest = os.path.join(DATA, "GSE207704_CLDN4_RNAseq.txt.gz")
    download(ftp(acc, "suppl/GSE207704_CLDN4_RNAseq.txt.gz"), dest)
    df = pd.read_csv(dest, sep="\t", index_col=0)
    # inspect columns
    print("GSE207704 cols", list(df.columns), "shape", df.shape)
    # expected: T47D, T47D, T47D CLDN4-/-, T47D CLDN4-/-, MCF7, MCF7, MCF7 CLDN4-/-, MCF7 CLDN4-/-
    # or gene-level with sample names
    rows = []
    # try to map columns
    cols = list(df.columns)
    print(df.head(2).to_string())
    # gene column may be the index (ensembl or symbol)
    gene_col = None
    if df.index.name and "gene" in str(df.index.name).lower():
        pass
    # add a symbol column guess
    first = df.index.astype(str)
    looks_ensembl = first.str.startswith("ENSG").mean() > 0.5
    if looks_ensembl:
        # last column sometimes gene name; else keep ensembl and match later
        pass
    groups = {
        "T47D_WT": [c for c in cols if "T47D" in c and "CLDN4" not in c and "-/-" not in c and "KO" not in c],
        "T47D_KO": [c for c in cols if "T47D" in c and ("CLDN4" in c or "-/-" in c or "KO" in c)],
        "MCF7_WT": [c for c in cols if "MCF7" in c and "CLDN4" not in c and "-/-" not in c and "KO" not in c],
        "MCF7_KO": [c for c in cols if "MCF7" in c and ("CLDN4" in c or "-/-" in c or "KO" in c)],
    }
    print("groups", {k: v for k, v in groups.items()})
    for gene in ["CD274", "CLDN4", "TACSTD2", "PDCD1LG2", "STAT1", "IRF1"]:
        rid = pick_gene_row(df, [gene])
        if rid is None:
            # try scanning a symbol column
            for col in df.columns:
                if df[col].dtype == object:
                    m = df[col].astype(str).str.upper() == gene
                    if m.any():
                        rid = df.index[m][0]
                        break
        if rid is None:
            print("  missing", gene)
            continue
        series = pd.to_numeric(df.loc[rid, cols], errors="coerce")
        for line, wt, ko in (("T47D", "T47D_WT", "T47D_KO"), ("MCF7", "MCF7_WT", "MCF7_KO")):
            if not groups[wt] or not groups[ko]:
                continue
            a = log2p(series[groups[wt]].values)
            b = log2p(series[groups[ko]].values)
            s = summarize(a, b, "parental", "CLDN4-/-")
            if not s:
                continue
            rows.append({
                "accession": acc, "model": f"{line} breast cancer cell line",
                "contrast": "CLDN4-/- vs parental",
                "qualifies_as": "ARM_A_CLDN4_genetic_loss",
                "gene": gene, "scale": "log2(count+1)",
                **{k: s[k] for k in ("n_ctrl", "n_treat", "mean_ctrl", "mean_treat",
                                     "log2FC_treat_minus_ctrl", "welch_p", "direction")},
                "note": f"n=2 vs 2; values {s['ctrl_values']} vs {s['treat_values']}",
            })
    return rows, {"columns": cols, "groups": groups, "shape": list(df.shape)}


# ---------------------------------------------------------------------------
# GSE22493  SKOV-3 CLDN4 siRNA / overexpression  (two-color log2)
# ---------------------------------------------------------------------------
def gse22493():
    acc = "GSE22493"
    dest = os.path.join(DATA, "GSE22493_series_matrix.txt.gz")
    download(ftp(acc, "matrix/GSE22493_series_matrix.txt.gz"), dest)
    # parse series matrix: sample titles + table
    sample_titles = []
    table_lines = []
    in_table = False
    with gzip.open(dest, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                sample_titles = [x.strip().strip('"') for x in line.split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line)
    df = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", index_col=0)
    print("GSE22493 samples", sample_titles, "shape", df.shape)
    # VALUE is log2(KD/OE). Gene annotation is not in the matrix — need platform.
    plat = os.path.join(DATA, "GPL10555.annot.gz")
    download("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10555/annot/GPL10555.annot.gz", plat)
    # annot may 404; try SOFT
    genes = {}
    try:
        with gzip.open(plat, "rt", errors="replace") as fh:
            header = None
            for line in fh:
                if line.startswith("#") or line.startswith("!"):
                    continue
                if header is None:
                    header = line.rstrip("\n").split("\t")
                    continue
                parts = line.rstrip("\n").split("\t")
                rec = dict(zip(header, parts))
                genes[rec.get("ID", parts[0])] = rec.get("Gene symbol") or rec.get("Gene symbol", "")
    except Exception as e:
        print("annot fail", e)
        # SOFT fallback
        soft = os.path.join(DATA, "GPL10555_family.soft.gz")
        try:
            download("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10555/soft/GPL10555_family.soft.gz", soft)
            with gzip.open(soft, "rt", errors="replace") as fh:
                pid = None
                for line in fh:
                    if line.startswith("^PLATFORM") or line.startswith("!platform_table_begin"):
                        continue
                    if line.startswith("ID\t") or line.startswith("ID "):
                        header = line.rstrip("\n").split("\t")
                        continue
                    if line.startswith("!platform_table_end"):
                        break
                    if line.startswith("#") or line.startswith("!"):
                        continue
                    if "header" in dir() and not line.startswith("^"):
                        parts = line.rstrip("\n").split("\t")
                        rec = dict(zip(header, parts + [""] * 10))
                        sid = rec.get("ID") or parts[0]
                        sym = rec.get("GENE_SYMBOL") or rec.get("Gene Symbol") or rec.get("ORF") or ""
                        genes[sid] = sym
        except Exception as e2:
            print("soft fail", e2)
    # map
    df["symbol"] = [genes.get(str(i), "") for i in df.index]
    print("mapped symbols", (df["symbol"] != "").sum(), "of", len(df))
    rows = []
    for gene in ["CD274", "CLDN4", "TACSTD2", "PDCD1LG2", "STAT1", "IRF1", "PD-L1", "B7-H1"]:
        sub = df[df["symbol"].str.upper() == gene.upper()]
        if sub.empty:
            # also search GB_ACC via contains in symbol field
            continue
        # median across probes, then mean across arrays
        num = sub.drop(columns=["symbol"]).apply(pd.to_numeric, errors="coerce")
        per_array = num.median(axis=0)
        vals = per_array.values.astype(float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        rows.append({
            "accession": acc,
            "model": "SKOV-3-IP-Luc ovarian; CLDN4 siRNA vs CLDN4 overexpression",
            "contrast": "deposited log2(KD / overexpression)",
            "qualifies_as": "ARM_A_CLDN4_genetic_loss_vs_OE_not_WT",
            "gene": gene,
            "scale": "log2(KD/OE) deposited VALUE",
            "n_ctrl": 0,
            "n_treat": int(len(vals)),
            "mean_ctrl": None,
            "mean_treat": float(np.mean(vals)),
            "log2FC_treat_minus_ctrl": float(np.mean(vals)),
            "welch_p": float(stats.ttest_1samp(vals, 0.0).pvalue) if len(vals) > 1 else None,
            "direction": "UP" if np.mean(vals) > 0 else "DOWN",
            "note": f"n_probes={len(sub)} per-array medians={list(map(float, vals))}; control is overexpression not WT",
        })
    if not any(r["gene"] == "CD274" for r in rows):
        rows.append({
            "accession": acc, "model": "SKOV-3-IP-Luc",
            "contrast": "CD274 probe search",
            "qualifies_as": "ARM_A_CLDN4_genetic_loss_vs_OE_not_WT",
            "gene": "CD274", "scale": "NA",
            "n_ctrl": 0, "n_treat": 0, "mean_ctrl": None, "mean_treat": None,
            "log2FC_treat_minus_ctrl": None, "welch_p": None, "direction": "NOT_ON_ARRAY_OR_UNANNOTATED",
            "note": f"mapped { (df['symbol']!='').sum() }/{len(df)} probes; CD274/PD-L1/B7-H1 not found in GPL10555 symbols",
        })
    return rows, {"n_mapped": int((df["symbol"] != "").sum()), "n_probes": int(len(df))}


# ---------------------------------------------------------------------------
# GSE22421  C-CPE treated vs untreated SKOV-3 (same platform)
# ---------------------------------------------------------------------------
def gse22421():
    acc = "GSE22421"
    dest = os.path.join(DATA, "GSE22421_series_matrix.txt.gz")
    download(ftp(acc, "matrix/GSE22421_series_matrix.txt.gz"), dest)
    table_lines = []
    in_table = False
    with gzip.open(dest, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line)
    df = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", index_col=0)
    # reuse GPL map if present
    genes = {}
    for cand in [os.path.join(DATA, "GPL10555.annot.gz"), os.path.join(DATA, "GPL10555_family.soft.gz")]:
        if not os.path.exists(cand):
            continue
        with gzip.open(cand, "rt", errors="replace") as fh:
            header = None
            for line in fh:
                if line.startswith("!") or line.startswith("#") or line.startswith("^"):
                    if line.startswith("!platform_table_begin"):
                        header = None
                    if line.startswith("!platform_table_end"):
                        break
                    continue
                if header is None:
                    header = line.rstrip("\n").split("\t")
                    continue
                parts = line.rstrip("\n").split("\t")
                rec = dict(zip(header, parts + [""] * 8))
                genes[rec.get("ID") or parts[0]] = rec.get("GENE_SYMBOL") or rec.get("Gene symbol") or rec.get("Gene Symbol") or ""
        if genes:
            break
    df["symbol"] = [genes.get(str(i), "") for i in df.index]
    rows = []
    for gene in ["CD274", "CLDN4", "TACSTD2"]:
        sub = df[df["symbol"].str.upper() == gene]
        if sub.empty:
            continue
        num = sub.drop(columns=["symbol"]).apply(pd.to_numeric, errors="coerce")
        vals = num.median(axis=0).values.astype(float)
        vals = vals[np.isfinite(vals)]
        rows.append({
            "accession": acc,
            "model": "SKOV-3 C-CPE (CLDN3/4 binder) vs untreated",
            "contrast": "deposited log2(C-CPE / untreated) if two-color, else as deposited",
            "qualifies_as": "NEAR_MISS_ligand_not_genetic_KD",
            "gene": gene, "scale": "deposited VALUE",
            "n_ctrl": 0, "n_treat": int(len(vals)),
            "mean_ctrl": None, "mean_treat": float(np.mean(vals)),
            "log2FC_treat_minus_ctrl": float(np.mean(vals)),
            "welch_p": float(stats.ttest_1samp(vals, 0.0).pvalue) if len(vals) > 1 else None,
            "direction": "UP" if np.mean(vals) > 0 else "DOWN",
            "note": f"not genetic CLDN4 loss; per-array {list(map(float, vals))}",
        })
    if not rows:
        rows.append({
            "accession": acc, "model": "SKOV-3 C-CPE",
            "contrast": "CD274 probe search",
            "qualifies_as": "NEAR_MISS_ligand_not_genetic_KD",
            "gene": "CD274", "scale": "NA",
            "n_ctrl": 0, "n_treat": 0, "mean_ctrl": None, "mean_treat": None,
            "log2FC_treat_minus_ctrl": None, "welch_p": None,
            "direction": "NOT_ON_ARRAY_OR_UNANNOTATED",
            "note": "same GPL10555 limitation as GSE22493",
        })
    return rows, {"n_mapped": int((df["symbol"] != "").sum())}


# ---------------------------------------------------------------------------
# GSE50927  mouse lung Cldn4 KO — deposited DE tables, not raw counts
# ---------------------------------------------------------------------------
def gse50927():
    acc = "GSE50927"
    dest = os.path.join(DATA, "GSE50927_Cldn4lungWTvsKOgenes.csv.gz")
    download(ftp(acc, "suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz"), dest)
    df = pd.read_csv(dest)
    print("GSE50927 cols", list(df.columns), df.shape)
    print(df.head(2).to_string())
    rows = []
    # try to find Cd274
    textcols = [c for c in df.columns if df[c].dtype == object]
    hit = None
    for c in textcols:
        m = df[c].astype(str).str.upper().isin(["CD274", "PDCD1LG1", "B7H1", "B7-H1"])
        if m.any():
            hit = df[m]
            break
    # also gene symbols like Cd274
    if hit is None:
        for c in textcols:
            m = df[c].astype(str).str.lower().isin(["cd274", "pdcd1lg1"])
            if m.any():
                hit = df[m]
                break
    note = "deposited DE table WT vs Cldn4 KO lung (also VILI-injured; not a cancer model)"
    if hit is None or hit.empty:
        rows.append({
            "accession": acc, "model": "mouse lung WT vs Cldn4 KO ± VILI",
            "contrast": "deposited WTvsKO gene table",
            "qualifies_as": "ARM_A_CLDN4_genetic_loss_NOT_CANCER",
            "gene": "Cd274", "scale": "NA",
            "n_ctrl": 0, "n_treat": 0, "mean_ctrl": None, "mean_treat": None,
            "log2FC_treat_minus_ctrl": None, "welch_p": None,
            "direction": "NOT_IN_DE_TABLE",
            "note": note + f"; columns={list(df.columns)}; Cd274 not listed (table may be filtered)",
        })
    else:
        rec = hit.iloc[0].to_dict()
        rows.append({
            "accession": acc, "model": "mouse lung WT vs Cldn4 KO ± VILI",
            "contrast": "deposited WTvsKO gene table",
            "qualifies_as": "ARM_A_CLDN4_genetic_loss_NOT_CANCER",
            "gene": "Cd274", "scale": "as_deposited",
            "n_ctrl": None, "n_treat": None,
            "mean_ctrl": None, "mean_treat": None,
            "log2FC_treat_minus_ctrl": rec.get("log2FoldChange") or rec.get("logFC") or rec.get("log2FC"),
            "welch_p": rec.get("pvalue") or rec.get("P.Value") or rec.get("pval"),
            "direction": "SEE_NOTE",
            "note": note + " | " + json.dumps({k: str(v)[:80] for k, v in rec.items()}),
        })
    return rows, {"columns": list(df.columns), "n": int(len(df))}


# ---------------------------------------------------------------------------
# GSE274940  EpH4 WT vs complete Cldn-null
# ---------------------------------------------------------------------------
def gse274940():
    acc = "GSE274940"
    dest = os.path.join(DATA, "GSE274940_raw_counts.csv.gz")
    download(ftp(acc, "suppl/GSE274940_raw_counts.csv.gz"), dest)
    df = pd.read_csv(dest, index_col=0)
    print("GSE274940 cols", list(df.columns), df.shape)
    cols = list(df.columns)
    wt = [c for c in cols if "WT" in c or "wt" in c]
    ko = [c for c in cols if "null" in c.lower() or "KO" in c or "Cldn" in c and "WT" not in c]
    if not wt or not ko:
        # fallback by name patterns from sample table
        wt = [c for c in cols if "WT" in str(c)]
        ko = [c for c in cols if "null" in str(c).lower()]
    print("wt", wt, "ko", ko)
    rows = []
    for gene in ["Cd274", "Cldn4", "Tacstd2", "Cd8a", "Stat1"]:
        rid = pick_gene_row(df, [gene, gene.upper()])
        if rid is None:
            print("  missing", gene)
            continue
        series = pd.to_numeric(df.loc[rid], errors="coerce")
        a = log2p(series[wt].values)
        b = log2p(series[ko].values)
        s = summarize(a, b, "WT", "Cldn-null")
        if not s:
            continue
        rows.append({
            "accession": acc,
            "model": "EpH4 mouse mammary epithelium, complete Cldn-null vs WT",
            "contrast": "Cldn-null vs WT",
            "qualifies_as": "NEAR_MISS_all_claudins_not_CLDN4",
            "gene": gene, "scale": "log2(count+1)",
            **{k: s[k] for k in ("n_ctrl", "n_treat", "mean_ctrl", "mean_treat",
                                 "log2FC_treat_minus_ctrl", "welch_p", "direction")},
            "note": f"not CLDN4-specific; {s['ctrl_values']} vs {s['treat_values']}",
        })
    return rows, {"wt": wt, "ko": ko}


# ---------------------------------------------------------------------------
# GSE312098  CX-1 IMMU132 (sacituzumab govitecan) 2-day
# ---------------------------------------------------------------------------
def gse312098():
    acc = "GSE312098"
    dest = os.path.join(DATA, "GSE312098_gene_fpkm.txt.gz")
    download(ftp(acc, "suppl/GSE312098_gene_fpkm.txt.gz"), dest)
    df = pd.read_csv(dest, sep="\t", index_col=0)
    print("GSE312098 cols", list(df.columns), df.shape)
    cols = list(df.columns)
    ctrl = [c for c in cols if "Control" in c or "control" in c]
    sg = [c for c in cols if "IMMU132" in c or "IMMU-132" in c]
    # exclude combination
    sg = [c for c in sg if "Combo" not in c and "Combination" not in c and "GSK" not in c]
    print("ctrl", ctrl, "sg", sg)
    rows = []
    for gene in ["CD274", "CLDN4", "TACSTD2", "PDCD1LG2", "STAT1", "IRF1"]:
        rid = pick_gene_row(df, [gene])
        if rid is None:
            # FPKM tables often have gene_name column
            if "gene_name" in df.columns:
                m = df["gene_name"].astype(str).str.upper() == gene
                if m.any():
                    rid = df.index[m][0]
            if rid is None:
                print("  missing", gene)
                continue
        numeric_cols = [c for c in cols if c not in ("gene_name", "gene_id", "gene_biotype")]
        series = pd.to_numeric(df.loc[rid, numeric_cols], errors="coerce")
        a = log2p(series[ctrl].values)
        b = log2p(series[sg].values)
        s = summarize(a, b, "Control", "IMMU132")
        if not s:
            continue
        rows.append({
            "accession": acc,
            "model": "CX-1 colorectal cell line, 2-day treatment",
            "contrast": "IMMU132 (sacituzumab govitecan) vs Control",
            "qualifies_as": "ARM_B_TROP2_ADC",
            "gene": gene, "scale": "log2(FPKM+1)",
            **{k: s[k] for k in ("n_ctrl", "n_treat", "mean_ctrl", "mean_treat",
                                 "log2FC_treat_minus_ctrl", "welch_p", "direction")},
            "note": f"{s['ctrl_values']} vs {s['treat_values']}",
        })
    return rows, {"ctrl": ctrl, "sg": sg, "all_cols": cols}


# ---------------------------------------------------------------------------
# GSE311016  CRC PDX IMMU132 29-day (paired by PDX)
# ---------------------------------------------------------------------------
def gse311016():
    acc = "GSE311016"
    dest = os.path.join(DATA, "GSE311016_gene_fpkm.txt.gz")
    download(ftp(acc, "suppl/GSE311016_gene_fpkm.txt.gz"), dest)
    df = pd.read_csv(dest, sep="\t", index_col=0)
    print("GSE311016 cols", list(df.columns), df.shape)
    cols = list(df.columns)
    ctrl = [c for c in cols if "Control" in c or "control" in c]
    sg = [c for c in cols if "IMMU132" in c]
    rows = []
    for gene in ["CD274", "CLDN4", "TACSTD2", "PDCD1LG2", "STAT1", "IRF1"]:
        rid = pick_gene_row(df, [gene])
        if rid is None:
            if "gene_name" in df.columns:
                m = df["gene_name"].astype(str).str.upper() == gene
                if m.any():
                    rid = df.index[m][0]
            if rid is None:
                print("  missing", gene)
                continue
        numeric_cols = [c for c in cols if c not in ("gene_name", "gene_id", "gene_biotype")]
        series = pd.to_numeric(df.loc[rid, numeric_cols], errors="coerce")
        a = log2p(series[ctrl].values)
        b = log2p(series[sg].values)
        s = summarize(a, b, "Control", "IMMU132")
        if not s:
            continue
        # paired by PDX id if names contain PDX
        note = f"unpaired Welch; {s['ctrl_values']} vs {s['treat_values']}"
        rows.append({
            "accession": acc,
            "model": "CRC PDX (5 models), 29-day",
            "contrast": "IMMU132 (sacituzumab govitecan) vs Control",
            "qualifies_as": "ARM_B_TROP2_ADC",
            "gene": gene, "scale": "log2(FPKM+1)",
            **{k: s[k] for k in ("n_ctrl", "n_treat", "mean_ctrl", "mean_treat",
                                 "log2FC_treat_minus_ctrl", "welch_p", "direction")},
            "note": note,
        })
    return rows, {"ctrl": ctrl, "sg": sg}


# ---------------------------------------------------------------------------
# GSE222450  SN-38 ± anti-PD1 mouse HNSCC
# ---------------------------------------------------------------------------
def gse222450():
    acc = "GSE222450"
    dest = os.path.join(DATA, "GSE222450_read_count.txt.gz")
    download(ftp(acc, "suppl/GSE222450_read_count.txt.gz"), dest)
    df = pd.read_csv(dest, sep="\t", index_col=0)
    print("GSE222450 cols", list(df.columns), df.shape)
    cols = list(df.columns)
    veh = [c for c in cols if "Vehicle" in c or "vehicle" in c or "Ctrl" in c or "control" in c.lower()]
    sn = [c for c in cols if "SN" in c and "PD" not in c and "anti" not in c.lower()]
    # fallback from known GSM order if names are GSM
    print("veh", veh, "sn", sn)
    rows = []
    for gene in ["Cd274", "Cldn4", "Tacstd2", "Cd8a", "Stat1", "Irf1"]:
        rid = pick_gene_row(df, [gene, gene.upper()])
        if rid is None:
            print("  missing", gene)
            continue
        series = pd.to_numeric(df.loc[rid], errors="coerce")
        if not veh or not sn:
            continue
        a = log2p(series[veh].values)
        b = log2p(series[sn].values)
        s = summarize(a, b, "Vehicle", "SN-38")
        if not s:
            continue
        rows.append({
            "accession": acc,
            "model": "mouse HNSCC, SN-38 vs vehicle (not a TROP2 ADC)",
            "contrast": "SN-38 vs Vehicle",
            "qualifies_as": "MECH_SN38_payload_not_TROP2ADC",
            "gene": gene, "scale": "log2(count+1)",
            **{k: s[k] for k in ("n_ctrl", "n_treat", "mean_ctrl", "mean_treat",
                                 "log2FC_treat_minus_ctrl", "welch_p", "direction")},
            "note": f"{s['ctrl_values']} vs {s['treat_values']}; paper claims SN-38 lowers PD-L1 via FoxO3a",
        })
    return rows, {"cols": cols, "veh": veh, "sn": sn}


# ---------------------------------------------------------------------------
# GSE334497  4T1 Trop2 KO vs WT tumors
# ---------------------------------------------------------------------------
def gse334497():
    acc = "GSE334497"
    dest = os.path.join(DATA, "GSE334497_normalized_counts.csv.gz")
    download(ftp(acc, "suppl/GSE334497_normalized_counts.csv.gz"), dest)
    df = pd.read_csv(dest, index_col=0)
    print("GSE334497 cols", list(df.columns), df.shape)
    cols = list(df.columns)
    wt = [c for c in cols if "WT" in c or "wt" in c]
    ko = [c for c in cols if "KO" in c or "ko" in c]
    rows = []
    for gene in ["Cd274", "Cldn4", "Cldn7", "Tacstd2", "Cd8a", "Stat1"]:
        rid = pick_gene_row(df, [gene, gene.upper()])
        if rid is None:
            print("  missing", gene)
            continue
        series = pd.to_numeric(df.loc[rid], errors="coerce")
        a = log2p(series[wt].values)
        b = log2p(series[ko].values)
        s = summarize(a, b, "WT", "Trop2 KO")
        if not s:
            continue
        rows.append({
            "accession": acc,
            "model": "4T1 mouse mammary tumor, Trop2 KO vs WT",
            "contrast": "Trop2 KO vs WT",
            "qualifies_as": "NEAR_MISS_TROP2_loss_not_ADC_not_CLDN4",
            "gene": gene, "scale": "log2(norm+1)",
            **{k: s[k] for k in ("n_ctrl", "n_treat", "mean_ctrl", "mean_treat",
                                 "log2FC_treat_minus_ctrl", "welch_p", "direction")},
            "note": f"{s['ctrl_values']} vs {s['treat_values']}; paper is TROP2/claudin-7 barrier, not PD-L1 induction",
        })
    return rows, {"wt": wt, "ko": ko}


# ---------------------------------------------------------------------------
# GSE245459  TACSTD2 shRNA SKOV-3
# ---------------------------------------------------------------------------
def gse245459():
    acc = "GSE245459"
    dest = os.path.join(DATA, "GSE245459_fpkm.anno.txt.gz")
    download(ftp(acc, "suppl/GSE245459_fpkm.anno.txt.gz"), dest)
    df = pd.read_csv(dest, sep="\t")
    print("GSE245459 cols", list(df.columns)[:20], df.shape)
    # identify gene and sample columns
    gene_col = None
    for c in df.columns:
        if str(c).lower() in ("gene_name", "gene", "symbol", "gene_symbol"):
            gene_col = c
            break
    if gene_col is None:
        # first column maybe gene
        gene_col = df.columns[0]
    df = df.set_index(gene_col)
    cols = [c for c in df.columns if c.lower() not in ("gene_id", "gene_biotype", "description", "length")]
    shnc = [c for c in cols if "shNC" in c or "shnc" in c or "NC" in c and "DDP" not in c and "CDDP" not in c]
    sh = [c for c in cols if (c.startswith("sh") or "shTAC" in c or "sh1" in c or "sh2" in c or "sh3" in c)
          and "NC" not in c and "DDP" not in c and "CDDP" not in c]
    print("shnc", shnc, "sh", sh)
    rows = []
    for gene in ["CD274", "CLDN4", "TACSTD2"]:
        rid = pick_gene_row(df, [gene])
        if rid is None:
            print("  missing", gene)
            continue
        if not shnc or not sh:
            continue
        series = pd.to_numeric(df.loc[rid, shnc + sh], errors="coerce")
        a = log2p(series[shnc].values)
        b = log2p(series[sh].values)
        s = summarize(a, b, "shNC", "TACSTD2 sh")
        if not s:
            continue
        rows.append({
            "accession": acc,
            "model": "SKOV-3 TACSTD2 shRNA ± cisplatin",
            "contrast": "TACSTD2 sh vs shNC (no cisplatin)",
            "qualifies_as": "NEAR_MISS_TROP2_KD_not_ADC",
            "gene": gene, "scale": "log2(FPKM+1)",
            **{k: s[k] for k in ("n_ctrl", "n_treat", "mean_ctrl", "mean_treat",
                                 "log2FC_treat_minus_ctrl", "welch_p", "direction")},
            "note": f"{s['ctrl_values']} vs {s['treat_values']}",
        })
    return rows, {"cols": cols[:30], "shnc": shnc, "sh": sh}


def main():
    all_rows = []
    meta = {}
    jobs = [
        ("GSE207704", gse207704),
        ("GSE22493", gse22493),
        ("GSE22421", gse22421),
        ("GSE50927", gse50927),
        ("GSE274940", gse274940),
        ("GSE312098", gse312098),
        ("GSE311016", gse311016),
        ("GSE222450", gse222450),
        ("GSE334497", gse334497),
        ("GSE245459", gse245459),
    ]
    for name, fn in jobs:
        print(f"\n######## {name} ########", flush=True)
        try:
            rows, info = fn()
            all_rows.extend(rows)
            meta[name] = {"ok": True, "info": info, "n_rows": len(rows)}
        except Exception as e:
            print("FAILED", name, e, flush=True)
            meta[name] = {"ok": False, "error": repr(e)}
            all_rows.append({
                "accession": name, "model": "", "contrast": "",
                "qualifies_as": "DOWNLOAD_OR_PARSE_FAIL",
                "gene": "CD274", "scale": "NA",
                "n_ctrl": None, "n_treat": None, "mean_ctrl": None, "mean_treat": None,
                "log2FC_treat_minus_ctrl": None, "welch_p": None,
                "direction": "FAIL", "note": repr(e),
            })
    write_contrast(all_rows, os.path.join(OUT, "contrasts.tsv"))
    with open(os.path.join(OUT, "run_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
    print("\nWrote", os.path.join(OUT, "contrasts.tsv"), "n=", len(all_rows))


if __name__ == "__main__":
    main()
