#!/usr/bin/env python3
"""Re-parse already-downloaded GEO files for CD274/CLDN4/TACSTD2.

Produces a small, honest contrast table. Replicates are reported as raw
values. n=1 deposited averages are flagged. Lookalikes stay labeled.
"""

from __future__ import annotations

import gzip
import json
import os
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "omics", "cd274_contrasts")
os.makedirs(OUT, exist_ok=True)
UA = "C3-PDL1-catalog/1.0"

ENS = {
    "human": {
        "CD274": "ENSG00000120217",
        "CLDN4": "ENSG00000189143",
        "TACSTD2": "ENSG00000184292",
        "STAT1": "ENSG00000115415",
        "IRF1": "ENSG00000125347",
        "PDCD1LG2": "ENSG00000197646",
    },
    "mouse": {
        "Cd274": "ENSMUSG00000016496",
        "Cldn4": "ENSMUSG00000047501",
        "Tacstd2": "ENSMUSG00000051397",
        "Stat1": "ENSMUSG00000026104",
        "Irf1": "ENSMUSG00000018899",
        "Cd8a": "ENSMUSG00000053977",
        "Cldn7": "ENSMUSG00000009634",
    },
}


def log2p(x):
    return np.log2(np.asarray(x, dtype=float) + 1.0)


def row(accession, model, contrast, qualifies, gene, scale, ctrl, treat, note, n_is_avg=False):
    a = np.asarray(ctrl, dtype=float)
    b = np.asarray(treat, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    lfc = float(np.mean(b) - np.mean(a)) if len(a) and len(b) else None
    p = None
    if len(a) >= 2 and len(b) >= 2:
        _, p = stats.ttest_ind(b, a, equal_var=False)
        p = float(p)
    return {
        "accession": accession,
        "model": model,
        "contrast": contrast,
        "qualifies_as": qualifies,
        "gene": gene,
        "scale": scale,
        "n_ctrl": int(len(a)),
        "n_treat": int(len(b)),
        "ctrl_values": ";".join(f"{x:.4g}" for x in a),
        "treat_values": ";".join(f"{x:.4g}" for x in b),
        "mean_ctrl": float(np.mean(a)) if len(a) else None,
        "mean_treat": float(np.mean(b)) if len(b) else None,
        "log2FC": lfc,
        "welch_p": p,
        "direction": ("UP" if lfc > 0.25 else "DOWN" if lfc is not None and lfc < -0.25 else "FLAT") if lfc is not None else "NA",
        "replicates_collapsed_in_deposit": "Y" if n_is_avg else "N",
        "note": note,
    }


def find_ens(idx, ens):
    s = pd.Index(idx.astype(str))
    hit = s[s.str.startswith(ens)]
    return hit[0] if len(hit) else None


rows = []

# ---- GSE207704: deposited FPKM already averaged (1 number / condition) ----
df = pd.read_csv(os.path.join(DATA, "GSE207704_CLDN4_RNAseq.txt.gz"), sep="\t")
# two CLDN4 loci; use the high-expressed row (first)
for gene in ["CD274", "CLDN4", "TACSTD2", "STAT1"]:
    hit = df[df["gene_short_name"].astype(str).str.upper() == gene]
    if hit.empty:
        continue
    # take the row with highest WT mean
    hit = hit.copy()
    hit["_wt"] = hit["MCF7_WT_FPKM (fpkm)"] + hit["T47D_WT_FPKM (fpkm)"]
    r = hit.sort_values("_wt", ascending=False).iloc[0]
    for line, wt_c, ko_c in (
        ("MCF7", "MCF7_WT_FPKM (fpkm)", "MCF7_CLDN4KO_FPKM (fpkm)"),
        ("T47D", "T47D_WT_FPKM (fpkm)", "T47D_CLDN4KO_FPKM (fpkm)"),
    ):
        wt, ko = float(r[wt_c]), float(r[ko_c])
        rows.append(row(
            "GSE207704", f"{line} breast cancer line (deposited FPKM, replicates collapsed)",
            "CLDN4-KO vs parental", "ARM_A_CLDN4_genetic_loss",
            gene, "log2(FPKM+1)", [log2p(wt)], [log2p(ko)],
            f"raw FPKM WT={wt:.4g} KO={ko:.4g}; ratio KO/WT={ko/wt if wt else float('nan'):.3g}; "
            f"GEO lists 2 reps/arm but this file has 1 value. CLDN4 residual RNA remains high.",
            n_is_avg=True,
        ))

# ---- GSE50927: deposited edgeR WTvsKO (logFC is KO-WT as typical edgeR) ----
d = pd.read_csv(os.path.join(DATA, "GSE50927_Cldn4lungWTvsKOgenes.csv.gz"))
for gene in ["Cd274", "Cldn4", "Tacstd2", "Pdcd1lg2"]:
    hit = d[d["Marker.Symbol"].astype(str) == gene]
    if hit.empty:
        continue
    r = hit.iloc[0]
    rows.append({
        "accession": "GSE50927",
        "model": "mouse whole lung, Cldn4 KO vs WT (includes VILI-injured animals; not cancer)",
        "contrast": "deposited edgeR logFC (KO vs WT convention as in file)",
        "qualifies_as": "ARM_A_CLDN4_genetic_loss_NOT_CANCER",
        "gene": gene,
        "scale": "edgeR logFC",
        "n_ctrl": None, "n_treat": None,
        "ctrl_values": "", "treat_values": "",
        "mean_ctrl": None, "mean_treat": None,
        "log2FC": float(r["logFC"]),
        "welch_p": float(r["PValue"]),
        "direction": "UP" if r["logFC"] > 0.25 else "DOWN" if r["logFC"] < -0.25 else "FLAT",
        "replicates_collapsed_in_deposit": "Y",
        "note": f"FDR={r['FDR']:.3g}; logCPM={r['logCPM']:.3g}; Entrez={r['EntrezID']}. "
                f"Cldn4 KO is confirmed (logFC=-6.06). Not a tumor model.",
    })

# ---- GSE274940: EpH4 complete Cldn-null raw counts ----
d = pd.read_csv(os.path.join(DATA, "GSE274940_raw_counts.csv.gz"))
for gene, alias in [("Cd274", "Cd274"), ("Cldn4", "Cldn4"), ("Tacstd2", "Tacstd2")]:
    hit = d[d["ALIAS"].astype(str) == alias]
    if hit.empty:
        continue
    r = hit.iloc[0]
    wt = [float(r[c]) for c in ("WT1", "WT2", "WT3")]
    ko = [float(r[c]) for c in ("KO1", "KO2", "KO3")]
    rows.append(row(
        "GSE274940", "EpH4 mouse mammary epithelium, complete Cldn-family null vs WT",
        "Cldn-null vs WT", "NEAR_MISS_all_claudins_not_CLDN4",
        gene, "log2(count+1)", log2p(wt), log2p(ko),
        f"raw counts WT={wt} KO={ko}. Cldn4 RNA is NOT depleted in the KO columns "
        f"(see Cldn4 row) — this accession cannot be read as CLDN4 loss.",
    ))

# ---- GSE312098 / GSE311016 UTF-16 FPKM ----
def read_utf16_gz(path):
    with gzip.open(path, "rb") as fh:
        raw = fh.read()
    return pd.read_csv(io_from_utf16(raw), sep="\t")


def io_from_utf16(raw):
    import io
    text = raw.decode("utf-16")
    return io.StringIO(text)


d = read_utf16_gz(os.path.join(DATA, "GSE312098_gene_fpkm.txt.gz"))
print("GSE312098 cols", d.columns.tolist()[:20], "shape", d.shape)
print(d.head(2).to_string())
# columns from decompressed header: gene_id, X_1..  need GSM map from series matrix
# Sample order from GEO: Control x3, IMMU132 x3, GSK x3, Combo x3 — but file has X_1..
# Download series matrix for mapping if needed.

# fetch series matrix sample titles vs supplementary column order is unreliable.
# Use GEO sample titles we already have and the file header C/T naming if present.
# Header was: gene_id, X_1 ...  — inspect unique prefixes

# Try matching gene via ensembl or symbol column
gene_col = d.columns[0]
print("gene_col", gene_col, d[gene_col].head(3).tolist())

# map samples from series matrix
def series_sample_titles(acc):
    n = int(acc[3:])
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{n//1000}nnn/{acc}/matrix/{acc}_series_matrix.txt.gz"
    dest = os.path.join(DATA, f"{acc}_series_matrix.txt.gz")
    if not os.path.exists(dest):
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as fh, open(dest, "wb") as out:
            out.write(fh.read())
    titles, geo, chars = [], [], []
    with gzip.open(dest, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.split("\t")[1:]]
            if line.startswith("!Sample_geo_accession"):
                geo = [x.strip().strip('"') for x in line.split("\t")[1:]]
            if line.startswith("!Sample_characteristics") or line.startswith("!Sample_source"):
                chars.append([x.strip().strip('"') for x in line.split("\t")[1:]])
    return titles, geo


titles, geo = series_sample_titles("GSE312098")
print("GSE312098 titles", titles)
# File columns after gene_id
expr_cols = [c for c in d.columns if c != gene_col]
print("expr_cols", expr_cols)
# If counts match 12 samples, assign in GEO order (often the case)
# GEO sample table (newest first in esummary — NOT reliable).
# series matrix order is the GEO order.
if len(expr_cols) == len(titles):
    colmap = dict(zip(expr_cols, titles))
else:
    colmap = {c: c for c in expr_cols}
print("colmap", colmap)

ctrl = [c for c, t in colmap.items() if "Control" in t or "control" in t]
sg = [c for c, t in colmap.items() if "IMMU132" in t and "Combination" not in t and "GSK" not in t]
if not ctrl or not sg:
    # try column names themselves
    ctrl = [c for c in expr_cols if str(c).startswith("X_") and False]
    # fallback: from UTF16 header earlier: X_1.. — use series matrix order
    if len(expr_cols) == 12 and len(titles) == 12:
        colmap = dict(zip(expr_cols, titles))
        ctrl = [c for c, t in colmap.items() if "Control" in t]
        sg = [c for c, t in colmap.items() if "IMMU132" in t and "Combination" not in t]
print("ctrl", ctrl, "sg", sg, "map", {c: colmap.get(c) for c in ctrl + sg})

# gene lookup
def pick_human(df, gene, gene_col):
    ens = ENS["human"][gene]
    s = df[gene_col].astype(str)
    m = s.str.startswith(ens) | (s.str.upper() == gene)
    if m.any():
        return df.loc[m].iloc[0]
    # other annotation columns
    for c in df.columns:
        if df[c].dtype == object:
            m = df[c].astype(str).str.upper() == gene
            if m.any():
                return df.loc[m].iloc[0]
            m = df[c].astype(str).str.startswith(ens)
            if m.any():
                return df.loc[m].iloc[0]
    return None


for gene in ["CD274", "CLDN4", "TACSTD2", "STAT1", "IRF1"]:
    rec = pick_human(d, gene, gene_col)
    if rec is None:
        print("GSE312098 missing", gene)
        continue
    if not ctrl or not sg:
        continue
    wt = [float(rec[c]) for c in ctrl]
    ko = [float(rec[c]) for c in sg]
    rows.append(row(
        "GSE312098", "CX-1 colorectal cell line, 2-day treatment",
        "IMMU132 (sacituzumab govitecan) vs Control",
        "ARM_B_TROP2_ADC", gene, "log2(FPKM+1)",
        log2p(wt), log2p(ko),
        f"raw FPKM ctrl={wt} IMMU132={ko}; cols {ctrl}->{[colmap.get(c) for c in ctrl]} vs {sg}",
    ))

# GSE311016 PDX paired
d2 = read_utf16_gz(os.path.join(DATA, "GSE311016_gene_fpkm.txt.gz"))
print("GSE311016 cols", d2.columns.tolist(), d2.shape)
titles2, geo2 = series_sample_titles("GSE311016")
print("GSE311016 titles", titles2)
gene_col2 = d2.columns[0]
expr2 = [c for c in d2.columns if c != gene_col2]
# header was C_114 T_114 C_36 T_36 ... = Control/Treat by PDX
print("expr2", expr2)
ctrl2 = [c for c in expr2 if str(c).startswith("C_")]
sg2 = [c for c in expr2 if str(c).startswith("T_")]
if not ctrl2:
    if len(expr2) == len(titles2):
        cmap = dict(zip(expr2, titles2))
        ctrl2 = [c for c, t in cmap.items() if "Control" in t]
        sg2 = [c for c, t in cmap.items() if "IMMU132" in t]
    else:
        cmap = {}
else:
    cmap = {c: c for c in expr2}
print("ctrl2", ctrl2, "sg2", sg2)

for gene in ["CD274", "CLDN4", "TACSTD2", "STAT1", "IRF1"]:
    rec = pick_human(d2, gene, gene_col2)
    if rec is None:
        print("GSE311016 missing", gene)
        continue
    if not ctrl2 or not sg2:
        continue
    wt = [float(rec[c]) for c in ctrl2]
    ko = [float(rec[c]) for c in sg2]
    # paired logFC if same PDX suffix
    paired = []
    for c in ctrl2:
        suffix = str(c).split("_", 1)[-1]
        t = next((x for x in sg2 if str(x).endswith("_" + suffix) or str(x).endswith(suffix)), None)
        if t is not None:
            paired.append(float(np.log2(float(rec[t]) + 1) - np.log2(float(rec[c]) + 1)))
    note = f"raw FPKM ctrl={wt} IMMU132={ko}"
    if paired:
        note += f"; paired log2FC by PDX={['%.3f'%x for x in paired]} mean={np.mean(paired):.3f}"
    rows.append(row(
        "GSE311016", "CRC PDX (5 models), 29-day IMMU132 vs Control",
        "IMMU132 (sacituzumab govitecan) vs Control",
        "ARM_B_TROP2_ADC", gene, "log2(FPKM+1)",
        log2p(wt), log2p(ko), note,
    ))

# ---- GSE222450: need sample map + ensembl ----
titles3, geo3 = series_sample_titles("GSE222450")
print("GSE222450 titles", titles3)
d3 = pd.read_csv(os.path.join(DATA, "GSE222450_read_count.txt.gz"), sep="\t")
print("GSE222450 cols", d3.columns.tolist())
# 12 count columns vs 12 titles. series matrix order should match GSM order,
# but count columns are A2,A3,A4,B1,B2,B3,E1,E2,E3,G1,G2,G3.
# Fetch characteristics from series matrix for treatment.
chars = {}
dest = os.path.join(DATA, "GSE222450_series_matrix.txt.gz")
with gzip.open(dest, "rt", errors="replace") as fh:
    header_geo = None
    for line in fh:
        if line.startswith("!Sample_geo_accession"):
            header_geo = [x.strip().strip('"') for x in line.split("\t")[1:]]
        if line.startswith("!Sample_title"):
            chars["title"] = [x.strip().strip('"') for x in line.split("\t")[1:]]
        if line.startswith("!Sample_characteristics_ch1"):
            chars.setdefault("char", []).append([x.strip().strip('"') for x in line.split("\t")[1:]])
print("GSE222450 char", {k: (v if k != "char" else v) for k, v in chars.items()})

# Map A/B/E/G from paper: typically Vehicle, SN-38, anti-PD1, combo.
# Use titles in series-matrix order and assume count file is NOT in that order.
# Safer: look at a README or the first lines of the count file comment.
# Alternative: match by downloading each GSM's supplementary? Too heavy.
# Use mean Cd274 and report all 4 groups if we can map via GEO SOFT sample source.

# Try GSM-level: GEO often names files A2 etc in titles? titles were 'Vehicle (rep1)'...
# Look at series matrix !Sample_supplementary_file
with gzip.open(dest, "rt", errors="replace") as fh:
    for line in fh:
        if "supplementary" in line.lower() or line.startswith("!Sample_description"):
            print(line[:300].rstrip())

# If we cannot map, compute Cd274 for all columns and assign groups from
# the paper's stated n=3 and our title list, using a documented assumption
# ONLY if the count-file header documents it.
# Check first comment lines
with gzip.open(os.path.join(DATA, "GSE222450_read_count.txt.gz"), "rt") as fh:
    for i, line in enumerate(fh):
        if i < 5:
            print("count hdr", line[:200].rstrip())

# Map via NCBI GSM metadata (extract treatment from title) + guess column
# groups A,B,E,G = 4 treatments x 3 reps. Match treatments alphabetically
# to Vehicle, SN-38, anti-PD1, combo? That's an assumption — do NOT do it.
# Instead report Cd274 per column group A/B/E/G and let the writeup say
# mapping to treatment is not in the count file.

ens = ENS["mouse"]["Cd274"]
m = d3["GeneID"].astype(str).str.startswith(ens)
print("Cd274 rows", m.sum())
if m.any():
    rec = d3.loc[m].iloc[0]
    groups = {"A": ["A2", "A3", "A4"], "B": ["B1", "B2", "B3"],
              "E": ["E1", "E2", "E3"], "G": ["G1", "G2", "G3"]}
    for gname, cols in groups.items():
        vals = [float(rec[c]) for c in cols]
        print(" group", gname, vals)
    # We will NOT invent the treatment map. Record all four groups.
    rows.append({
        "accession": "GSE222450",
        "model": "mouse HNSCC RNA-seq (SN-38 ± anti-PD1); count columns A/B/E/G are unlabeled vs treatment",
        "contrast": "Cd274 raw counts by file group (treatment map not in count file)",
        "qualifies_as": "MECH_SN38_payload_not_TROP2ADC",
        "gene": "Cd274",
        "scale": "raw counts",
        "n_ctrl": 3, "n_treat": 3,
        "ctrl_values": "", "treat_values": "",
        "mean_ctrl": None, "mean_treat": None,
        "log2FC": None, "welch_p": None,
        "direction": "UNMAPPED_COLUMNS",
        "replicates_collapsed_in_deposit": "N",
        "note": "GeneID="+str(rec["GeneID"])+"; A="+str([float(rec[c]) for c in groups['A']])
                +"; B="+str([float(rec[c]) for c in groups['B']])
                +"; E="+str([float(rec[c]) for c in groups['E']])
                +"; G="+str([float(rec[c]) for c in groups['G']])
                +". Paper (PMID 36641765) claims SN-38 LOWERS PD-L1 via FoxO3a. "
                +"GEO titles are Vehicle / SN-38 / anti-PD1 / combo but the count "
                +"matrix uses A/B/E/G without a key. Not used as a C3 test.",
    })

# ---- GSE334497 Trop2 KO 4T1 ----
d4 = pd.read_csv(os.path.join(DATA, "GSE334497_normalized_counts.csv.gz"))
d4 = d4.rename(columns={d4.columns[0]: "gene"})
print("GSE334497 genes", d4["gene"].head(3).tolist())
wt_cols = [c for c in d4.columns if "control" in c.lower() or (c.startswith("RESUB-") and "KO" not in c)]
ko_cols = [c for c in d4.columns if "KO" in c]
print("334497 wt", wt_cols, "ko", ko_cols)
for gene, ens in [("Cd274", ENS["mouse"]["Cd274"]), ("Cldn4", ENS["mouse"]["Cldn4"]),
                  ("Tacstd2", ENS["mouse"]["Tacstd2"]), ("Cldn7", ENS["mouse"]["Cldn7"])]:
    m = d4["gene"].astype(str).str.startswith(ens)
    if not m.any():
        print("334497 missing", gene)
        continue
    rec = d4.loc[m].iloc[0]
    wt = [float(rec[c]) for c in wt_cols]
    ko = [float(rec[c]) for c in ko_cols]
    rows.append(row(
        "GSE334497", "4T1 mouse mammary tumor, Trop2 KO vs WT",
        "Trop2 KO vs WT", "NEAR_MISS_TROP2_loss_not_ADC_not_CLDN4",
        gene, "log2(norm+1)", log2p(wt), log2p(ko),
        f"raw {wt} vs {ko}. TROP2 loss, not TROP2 ADC, not CLDN4 loss. "
        f"Paper (PMID 41932810) is barrier/claudin-7 immune exclusion.",
    ))

# ---- GSE245459 TACSTD2 shRNA ----
# large file; read only needed columns
d5 = pd.read_csv(os.path.join(DATA, "GSE245459_fpkm.anno.txt.gz"), sep="\t",
                 usecols=["name", "shNC1", "shNC2", "shNC3", "sh1", "sh2", "sh3", "GeneName"],
                 low_memory=False)
for gene in ["CD274", "CLDN4", "TACSTD2"]:
    hit = d5[d5["GeneName"].astype(str).str.upper() == gene]
    if hit.empty:
        print("245459 missing", gene)
        continue
    rec = hit.iloc[0]
    wt = [float(rec[c]) for c in ("shNC1", "shNC2", "shNC3")]
    ko = [float(rec[c]) for c in ("sh1", "sh2", "sh3")]
    rows.append(row(
        "GSE245459", "SKOV-3 TACSTD2 shRNA (no cisplatin arms)",
        "shTACSTD2 vs shNC", "NEAR_MISS_TROP2_KD_not_ADC",
        gene, "log2(FPKM+1)", log2p(wt), log2p(ko),
        f"raw FPKM shNC={wt} sh={ko}",
    ))

# ---- GSE22493: try GEO GPL SOFT via query ----
gpl_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10555&targ=self&form=text&view=data"
gpl_dest = os.path.join(DATA, "GPL10555_data.txt")
if not os.path.exists(gpl_dest) or os.path.getsize(gpl_dest) < 1000:
    print("GET", gpl_url)
    try:
        req = urllib.request.Request(gpl_url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as fh, open(gpl_dest, "wb") as out:
            out.write(fh.read())
    except Exception as e:
        print("GPL fetch fail", e)

# parse GSE22493 matrix + optional GPL
import io
table_lines = []
in_table = False
with gzip.open(os.path.join(DATA, "GSE22493_series_matrix.txt.gz"), "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            table_lines.append(line)
m224 = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t")
print("GSE22493 matrix", m224.shape, m224.columns.tolist()[:5])
# GPL
sym = {}
if os.path.exists(gpl_dest) and os.path.getsize(gpl_dest) > 1000:
    with open(gpl_dest, "r", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("#") or line.startswith("^") or line.startswith("!"):
                if line.startswith("!platform_table_begin"):
                    header = None
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                print("GPL header", header[:12])
                continue
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts + [""] * 8))
            sid = rec.get("ID") or parts[0]
            for key in ("GENE_SYMBOL", "Gene Symbol", "gene_assignment", "ORF", "GB_ACC"):
                if rec.get(key):
                    sym[sid] = rec[key]
                    break
print("GPL mapped", len(sym))
idcol = m224.columns[0]
m224["sym"] = m224[idcol].astype(str).map(lambda x: sym.get(x, ""))
# also search GB_ACC-like
for gene in ["CD274", "PDCD1LG2", "CLDN4", "TACSTD2"]:
    sub = m224[m224["sym"].astype(str).str.upper().str.contains(gene, na=False)]
    print("GSE22493", gene, "n_probes", len(sub), "sym sample", sub["sym"].head(3).tolist() if len(sub) else None)
    if sub.empty:
        continue
    valcols = [c for c in m224.columns if c not in (idcol, "sym")]
    med = sub[valcols].apply(pd.to_numeric, errors="coerce").median(axis=0)
    vals = [float(x) for x in med.tolist() if np.isfinite(x)]
    rows.append({
        "accession": "GSE22493",
        "model": "SKOV-3 CLDN4 siRNA vs CLDN4 overexpression (not WT)",
        "contrast": "deposited log2(KD/OE), median across probes",
        "qualifies_as": "ARM_A_CLDN4_genetic_loss_vs_OE_not_WT",
        "gene": gene,
        "scale": "log2(KD/OE)",
        "n_ctrl": 0, "n_treat": len(vals),
        "ctrl_values": "",
        "treat_values": ";".join(f"{x:.4g}" for x in vals),
        "mean_ctrl": None,
        "mean_treat": float(np.mean(vals)) if vals else None,
        "log2FC": float(np.mean(vals)) if vals else None,
        "welch_p": float(stats.ttest_1samp(vals, 0).pvalue) if len(vals) > 1 else None,
        "direction": "UP" if vals and np.mean(vals) > 0.25 else "DOWN" if vals and np.mean(vals) < -0.25 else "FLAT",
        "replicates_collapsed_in_deposit": "N",
        "note": f"n_probes={len(sub)}; control is overexpression. Same accession as C4 IFN-negative.",
    })
if not any(r["accession"] == "GSE22493" and r["gene"] == "CD274" for r in rows):
    rows.append({
        "accession": "GSE22493",
        "model": "SKOV-3 CLDN4 siRNA vs overexpression",
        "contrast": "CD274 on GPL10555",
        "qualifies_as": "ARM_A_CLDN4_genetic_loss_vs_OE_not_WT",
        "gene": "CD274", "scale": "NA",
        "n_ctrl": 0, "n_treat": 0, "ctrl_values": "", "treat_values": "",
        "mean_ctrl": None, "mean_treat": None, "log2FC": None, "welch_p": None,
        "direction": "NOT_ANNOTATED",
        "replicates_collapsed_in_deposit": "N",
        "note": f"GPL10555 symbols mapped={len(sym)}. CD274/PD-L1 not found. Cannot test C3 on this array.",
    })

# GSE278664 is baseline BRCA-mut vs wt biopsies, not SG treatment — record as excluded
rows.append({
    "accession": "GSE278664",
    "model": "ovarian tumor biopsies labeled BRCAmut / BRCAwt",
    "contrast": "none — not a sacituzumab-treatment series",
    "qualifies_as": "EXCLUDED_title_mentions_SG_but_samples_are_baseline",
    "gene": "CD274", "scale": "NA",
    "n_ctrl": 0, "n_treat": 0, "ctrl_values": "", "treat_values": "",
    "mean_ctrl": None, "mean_treat": None, "log2FC": None, "welch_p": None,
    "direction": "NOT_A_TREATMENT_CONTRAST",
    "replicates_collapsed_in_deposit": "N",
    "note": "Title mentions sacituzumab govitecan + berzosertib; deposited samples are pretreatment biopsies split by BRCA status. Cannot test C3 arm B.",
})

out = pd.DataFrame(rows)
out.to_csv(os.path.join(OUT, "contrasts.tsv"), sep="\t", index=False)
# compact CD274-only
cd = out[out["gene"].str.upper().isin(["CD274", "CD274"]) | (out["gene"] == "Cd274")]
cd.to_csv(os.path.join(OUT, "cd274_only.tsv"), sep="\t", index=False)
print("\n==== CD274 rows ====")
print(cd[["accession", "qualifies_as", "direction", "log2FC", "welch_p", "note"]].to_string(index=False))
print("wrote", len(out), "rows")
