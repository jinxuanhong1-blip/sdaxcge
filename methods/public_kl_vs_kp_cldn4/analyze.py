#!/usr/bin/env python3
"""Public mouse combinatorial contrasts: KL-resistant vs KP/K-sensitive.

ADDITIVE. Cldn4-only. No private 8-KL matrices. No Tacstd2 x Cldn4 dual-high.
One contrast per accession (not a mega-merge). Honest n = mice or tumors.

Accessions
  GSE6135   KL vs KP and KL vs K bulk microarray (Cldn4 + T/NK + IFN/MHC)
  GSE154989 K/KP Smart-seq2 epithelium (no T/NK; K vs KP Cldn4; KP-only IFN/MHC)
  GSE179502 KT;Lkb1 XTR neoplastic epithelium (NonRestored ≈ KL vs Restored)
  GSE267321 LKR13 K / KK / KLK non-malignant (T/NK by genotype; Cldn4 empty)
"""

from __future__ import annotations

import gzip
import json
import re
import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gene_sets import (  # noqa: E402
    CLDN4,
    EPI_MARKERS,
    HOST_LUNG,
    IFN,
    IFN_MHC,
    MHC,
    NK_MARKERS,
    T_MARKERS,
    T_NK,
)

CACHE = Path("/tmp/geo")
OUT = HERE / "tables"
OUT.mkdir(parents=True, exist_ok=True)

GEO = {
    "GSE6135_matrix": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6135/matrix/"
        "GSE6135-GPL8321_series_matrix.txt.gz",
        CACHE / "GSE6135" / "GSE6135-GPL8321_series_matrix.txt.gz",
    ),
    "GPL8321": (
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8321/annot/GPL8321.annot.gz",
        CACHE / "GPL8321" / "GPL8321.annot.gz",
    ),
    "GSE154989_h5": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/"
        "GSE154989_mmLungPlate_fQC_dSp_normTPM.h5",
        CACHE / "GSE154989" / "GSE154989_mmLungPlate_fQC_dSp_normTPM.h5",
    ),
    "GSE154989_smp": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/"
        "GSE154989_mmLungPlate_fQC_smpTable.csv.gz",
        CACHE / "GSE154989" / "GSE154989_mmLungPlate_fQC_smpTable.csv.gz",
    ),
    "GSE154989_gene": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/"
        "GSE154989_mmLungPlate_fQC_geneTable.csv.gz",
        CACHE / "GSE154989" / "GSE154989_mmLungPlate_fQC_geneTable.csv.gz",
    ),
    "GSE179502_mtx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179502/suppl/"
        "GSE179502_XTR_sorted_scRNAseq_matrix.mtx.gz",
        CACHE / "GSE179502" / "GSE179502_XTR_sorted_scRNAseq_matrix.mtx.gz",
    ),
    "GSE179502_feat": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179502/suppl/"
        "GSE179502_XTR_sorted_scRNAseq_features.tsv.gz",
        CACHE / "GSE179502" / "GSE179502_XTR_sorted_scRNAseq_features.tsv.gz",
    ),
    "GSE179502_bc": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179502/suppl/"
        "GSE179502_XTR_sorted_scRNAseq_barcodes.tsv.gz",
        CACHE / "GSE179502" / "GSE179502_XTR_sorted_scRNAseq_barcodes.tsv.gz",
    ),
    "GSE267321_csv": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267321/suppl/"
        "GSE267321_Normalized_expression_matrix_02122026.csv.gz",
        CACHE / "GSE267321" / "GSE267321_Normalized_expression_matrix_02122026.csv.gz",
    ),
}


def download(key: str) -> Path:
    url, dest = GEO[key]
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    import urllib.request

    print(f"downloading {key} -> {dest}")
    urllib.request.urlretrieve(url, dest)
    return dest


def _fmt(x) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def welch_mwu(a, b) -> dict:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_a": float(np.mean(a)) if len(a) else np.nan,
        "mean_b": float(np.mean(b)) if len(b) else np.nan,
        "delta_a_minus_b": np.nan,
        "welch_p": np.nan,
        "mwu_p": np.nan,
        "rank_biserial": np.nan,
    }
    if len(a) and len(b):
        out["delta_a_minus_b"] = out["mean_a"] - out["mean_b"]
    if len(a) >= 2 and len(b) >= 2:
        _, out["welch_p"] = stats.ttest_ind(a, b, equal_var=False)
        u, out["mwu_p"] = stats.mannwhitneyu(a, b, alternative="two-sided")
        # rank-biserial: 1 - 2U/(n1 n2) with U for a>b convention via (mean ranks)
        n1, n2 = len(a), len(b)
        # Positive = group a ranks higher than group b.
        out["rank_biserial"] = float((2.0 * u) / (n1 * n2) - 1.0)
    return out


def spearman(x, y) -> tuple[float, float]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p)


def present(index, genes) -> list[str]:
    idx = set(map(str, index))
    return [g for g in genes if g in idx]


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=1)
    sd = df.std(axis=1, ddof=0).replace(0, np.nan)
    return df.sub(mu, axis=0).div(sd, axis=0)


# ---------------------------------------------------------------------------
# GSE6135
# ---------------------------------------------------------------------------

def parse_gse6135_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    titles = geo = chars = None
    table_lines = []
    in_table = False
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                geo = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and chars is None:
                chars = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                table_lines.append(line)
    expr = pd.read_csv(StringIO("".join(table_lines)), sep="\t")
    expr = expr.set_index("ID_REF")
    expr.columns = [c.strip('"') for c in expr.columns]
    meta = pd.DataFrame({"geo": geo, "title": titles, "characteristics": chars})
    return expr, meta


def load_gpl8321(path: Path) -> pd.DataFrame:
    rows = []
    header = None
    go = False
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                go = True
                header = next(f).rstrip().split("\t")
                continue
            if line.startswith("!platform_table_end"):
                break
            if go:
                rows.append(line.rstrip("\n").split("\t"))
    ann = pd.DataFrame(rows, columns=header)
    return ann[["ID", "Gene symbol"]].rename(columns={"ID": "probe", "Gene symbol": "symbol"})


def classify_gse6135(title: str, characteristics: str) -> dict:
    blob = f"{title} {characteristics}"
    low = blob.lower()
    is_met = bool(re.search(r"\bmet\b", low))
    hist = "unknown"
    if re.search(r"ad-sq|adenosquamous|mixed", low):
        hist = "Ad-sq"
    elif re.search(r"\bsq\b|squamous", low):
        hist = "Sq"
    elif re.search(r"\bad\b|adenocarcinoma", low):
        hist = "Ad"
    lkb = "none"
    m = re.search(r"lkb1?\s*(L/[L+\-])", blob, flags=re.IGNORECASE)
    if m:
        tok = m.group(1).upper()
        lkb = {"L/L": "L/L", "L/-": "L/-", "L/+": "L/+"}.get(tok, tok)
    elif "lkb" in low:
        lkb = "Lkb1_unparsed"
    if "p53" in low:
        geno = "KP"
        arm = "KP"
    elif "p16" in low or "ink4" in low:
        geno = "K_Ink4a"
        arm = "K_Ink4a"
    elif "lkb" in low:
        if lkb in {"L/+"}:
            geno = "KL_het"
            arm = "exclude_het"
        else:
            geno = "KL"
            arm = "KL"
    elif "kras" in low or "k-ras" in low:
        geno = "K"
        arm = "K"
    else:
        geno = "other"
        arm = "other"
    if is_met:
        arm = "exclude_met" if geno == "KL" else arm
    head = title.split("||")[0].strip().replace(" ", "")
    m2 = re.match(r"([0-9]+)", head)
    mouse_num = m2.group(1) if m2 else head
    return {
        "lkb_allele": lkb,
        "is_met": is_met,
        "histology": hist,
        "genotype": geno,
        "primary_arm": arm,
        "mouse_num": mouse_num,
        "mouse_id": f"{geno}_{mouse_num}",
    }


def collapse_symbol_matrix(expr: pd.DataFrame, ann: pd.DataFrame, wanted: list[str]) -> pd.DataFrame:
    probe_to_syms: dict[str, list[str]] = {}
    for _, r in ann.iterrows():
        raw = str(r["symbol"]) if pd.notna(r["symbol"]) else ""
        if not raw or raw == "nan":
            continue
        probe_to_syms[r["probe"]] = [s.strip() for s in raw.split("///") if s.strip()]
    gene_rows = {g: [] for g in wanted}
    for probe, syms in probe_to_syms.items():
        if probe not in expr.index:
            continue
        for s in syms:
            if s in gene_rows:
                gene_rows[s].append(probe)
    out = {}
    for g, probes in gene_rows.items():
        if probes:
            out[g] = expr.loc[probes].mean(axis=0)
    return pd.DataFrame(out).T


def run_gse6135() -> dict:
    expr, meta = parse_gse6135_matrix(download("GSE6135_matrix"))
    ann = load_gpl8321(download("GPL8321"))
    wanted = CLDN4 + T_NK + IFN_MHC
    gene_expr = collapse_symbol_matrix(expr, ann, wanted)
    parsed = [classify_gse6135(t, c) for t, c in zip(meta["title"], meta["characteristics"])]
    meta = pd.concat([meta, pd.DataFrame(parsed)], axis=1)
    meta = meta.set_index("geo")
    z = zscore_rows(gene_expr)
    tnk_present = present(gene_expr.index, T_NK)
    ifn_present = present(gene_expr.index, IFN)
    mhc_present = present(gene_expr.index, MHC)
    ifnmhc_present = present(gene_expr.index, IFN_MHC)
    tumor = meta.copy()
    tumor["cldn4"] = gene_expr.loc["Cldn4"] if "Cldn4" in gene_expr.index else np.nan
    tumor["tnk_score"] = z.loc[tnk_present].mean(axis=0) if tnk_present else np.nan
    tumor["ifn_score"] = z.loc[ifn_present].mean(axis=0) if ifn_present else np.nan
    tumor["mhc_score"] = z.loc[mhc_present].mean(axis=0) if mhc_present else np.nan
    tumor["ifn_mhc_score"] = z.loc[ifnmhc_present].mean(axis=0) if ifnmhc_present else np.nan
    tumor = tumor.reset_index().rename(columns={"geo": "sample_id"})
    tumor["accession"] = "GSE6135"
    tumor["unit"] = "tumor"
    tumor["in_primary_contrast"] = tumor["primary_arm"].isin(["KL", "KP", "K"])

    g = tumor.groupby("mouse_id", sort=False)
    mouse = pd.DataFrame(
        {
            "accession": "GSE6135",
            "unit": "mouse",
            "mouse_id": g.size().index,
            "genotype": g["genotype"].first().values,
            "primary_arm": g["primary_arm"].first().values,
            "lkb_allele": g["lkb_allele"].first().values,
            "histology": g["histology"].apply(lambda s: ",".join(sorted(set(s)))).values,
            "n_tumors": g.size().values,
            "sample_ids": g["sample_id"].apply(lambda s: ",".join(s)).values,
            "cldn4": g["cldn4"].mean().values,
            "tnk_score": g["tnk_score"].mean().values,
            "ifn_score": g["ifn_score"].mean().values,
            "mhc_score": g["mhc_score"].mean().values,
            "ifn_mhc_score": g["ifn_mhc_score"].mean().values,
            "in_primary_contrast": g["in_primary_contrast"].first().values,
        }
    )

    # Mouse-level primary: drop mice whose only tumors are met/het.
    # For KL mouse 592, primary_arm of met row is exclude_met but T1/T2 are KL.
    # Recompute mouse primary_arm from non-met tumors.
    mouse_fix = []
    for mid, sub in tumor.groupby("mouse_id"):
        prim = sub.loc[~sub["is_met"]]
        use = prim if len(prim) else sub
        mouse_fix.append(
            {
                "mouse_id": mid,
                "n_primary_tumors": int((~sub["is_met"]).sum()),
                "primary_arm_mouse": use["primary_arm"].iloc[0]
                if use["primary_arm"].nunique() == 1
                else use.loc[use["primary_arm"] == "KL", "primary_arm"].iloc[0]
                if (use["primary_arm"] == "KL").any()
                else use["primary_arm"].iloc[0],
                "cldn4_primary": float(prim["cldn4"].mean()) if len(prim) else np.nan,
                "tnk_primary": float(prim["tnk_score"].mean()) if len(prim) else np.nan,
                "ifn_primary": float(prim["ifn_score"].mean()) if len(prim) else np.nan,
                "mhc_primary": float(prim["mhc_score"].mean()) if len(prim) else np.nan,
                "ifn_mhc_primary": float(prim["ifn_mhc_score"].mean()) if len(prim) else np.nan,
            }
        )
    mouse = mouse.merge(pd.DataFrame(mouse_fix), on="mouse_id")

    contrasts = []
    for unit_name, df, arm_col, val_map in [
        (
            "tumor",
            tumor.loc[tumor["in_primary_contrast"]],
            "primary_arm",
            {
                "cldn4": "cldn4",
                "tnk_score": "tnk_score",
                "ifn_score": "ifn_score",
                "mhc_score": "mhc_score",
                "ifn_mhc_score": "ifn_mhc_score",
            },
        ),
        (
            "mouse",
            mouse.loc[mouse["primary_arm_mouse"].isin(["KL", "KP", "K"])],
            "primary_arm_mouse",
            {
                "cldn4": "cldn4_primary",
                "tnk_score": "tnk_primary",
                "ifn_score": "ifn_primary",
                "mhc_score": "mhc_primary",
                "ifn_mhc_score": "ifn_mhc_primary",
            },
        ),
    ]:
        for metric, col in val_map.items():
            for a, b, cname in [("KL", "KP", "KL_minus_KP"), ("KL", "K", "KL_minus_K")]:
                aa = df.loc[df[arm_col] == a, col]
                bb = df.loc[df[arm_col] == b, col]
                st = welch_mwu(aa, bb)
                contrasts.append(
                    {
                        "accession": "GSE6135",
                        "contrast": cname,
                        "unit": unit_name,
                        "metric": metric,
                        "n_KL": st["n_a"],
                        "n_sensitive": st["n_b"],
                        "mean_KL": st["mean_a"],
                        "mean_sensitive": st["mean_b"],
                        "delta_KL_minus_sensitive": st["delta_a_minus_b"],
                        "welch_p": st["welch_p"],
                        "mwu_p": st["mwu_p"],
                        "rank_biserial": st["rank_biserial"],
                        "status": "ok",
                        "note": (
                            "KL = Lkb1 L/L or L/- primary tumors; excluded L/+ het and metastasis. "
                            "T/NK and IFN/MHC are mean z-scores of present genes across the 25 arrays."
                        ),
                    }
                )

    tumor_out = tumor[
        [
            "accession",
            "unit",
            "sample_id",
            "mouse_id",
            "genotype",
            "primary_arm",
            "lkb_allele",
            "is_met",
            "histology",
            "title",
            "cldn4",
            "tnk_score",
            "ifn_score",
            "mhc_score",
            "ifn_mhc_score",
            "in_primary_contrast",
        ]
    ]
    mouse_out = mouse[
        [
            "accession",
            "unit",
            "mouse_id",
            "genotype",
            "primary_arm_mouse",
            "lkb_allele",
            "histology",
            "n_tumors",
            "n_primary_tumors",
            "sample_ids",
            "cldn4_primary",
            "tnk_primary",
            "ifn_primary",
            "mhc_primary",
            "ifn_mhc_primary",
        ]
    ].rename(
        columns={
            "primary_arm_mouse": "primary_arm",
            "cldn4_primary": "cldn4",
            "tnk_primary": "tnk_score",
            "ifn_primary": "ifn_score",
            "mhc_primary": "mhc_score",
            "ifn_mhc_primary": "ifn_mhc_score",
        }
    )
    tumor_out.to_csv(OUT / "GSE6135_tumors.tsv", sep="\t", index=False)
    mouse_out.to_csv(OUT / "GSE6135.tsv", sep="\t", index=False)
    pd.DataFrame(contrasts).to_csv(OUT / "GSE6135_contrasts.tsv", sep="\t", index=False)
    units = mouse_out
    return {
        "accession": "GSE6135",
        "present": {
            "Cldn4": present(gene_expr.index, CLDN4),
            "T_NK": tnk_present,
            "IFN": ifn_present,
            "MHC": mhc_present,
        },
        "n_tumors": int(len(tumor)),
        "n_mice": int(len(mouse)),
        "n_KL_mice": int((mouse["primary_arm_mouse"] == "KL").sum()),
        "n_KP_mice": int((mouse["primary_arm_mouse"] == "KP").sum()),
        "n_K_mice": int((mouse["primary_arm_mouse"] == "K").sum()),
        "contrasts": contrasts,
        "units": units,
    }


# ---------------------------------------------------------------------------
# GSE154989
# ---------------------------------------------------------------------------

def _read_csv_gz(path: Path) -> pd.DataFrame:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return pd.read_csv(fh)
    return pd.read_csv(path)


def run_gse154989() -> dict:
    import h5py

    gene_df = _read_csv_gz(download("GSE154989_gene"))
    smp = _read_csv_gz(download("GSE154989_smp"))
    symbols = gene_df["geneSymbol"].astype(str).tolist()
    n_cells = len(smp)
    needed = set(CLDN4 + T_NK + IFN_MHC + ["Epcam", "Ptprc"])
    want_idx = {i + 1: g for i, g in enumerate(symbols) if g in needed}
    vecs = {g: np.zeros(n_cells, dtype=np.float64) for g in needed if g in set(symbols)}
    h5_path = download("GSE154989_h5")
    with h5py.File(h5_path, "r") as f:
        ii = f["i"][0].astype(np.int32)
        jj = f["j"][0].astype(np.int32)
        vv = f["v"][0].astype(np.float64)
    for k in range(len(ii)):
        g = want_idx.get(int(ii[k]))
        if g is None:
            continue
        vecs[g][int(jj[k]) - 1] = vv[k]

    smp = smp.copy()
    smp["genotype"] = smp["mouseID"].astype(str).str.split("_", n=1).str[0]
    # T_early mouseIDs start with T
    smp.loc[smp["typeID"].astype(str).str.startswith("01_T"), "genotype"] = "T"
    smp["mouse"] = smp["mouseID"].astype(str).str.replace(r"_T\d+$", "", regex=True)
    smp["cldn4"] = vecs.get("Cldn4", np.zeros(n_cells))
    ifn_p = [g for g in IFN if g in vecs]
    mhc_p = [g for g in MHC if g in vecs]
    tnk_p = [g for g in T_NK if g in vecs]
    smp["ifn_score"] = np.mean([np.log1p(vecs[g]) for g in ifn_p], axis=0) if ifn_p else np.nan
    smp["mhc_score"] = np.mean([np.log1p(vecs[g]) for g in mhc_p], axis=0) if mhc_p else np.nan
    smp["ifn_mhc_score"] = (
        np.mean([np.log1p(vecs[g]) for g in ifn_p + mhc_p], axis=0) if (ifn_p + mhc_p) else np.nan
    )
    # Leak audit only — not a T/NK fraction.
    leak = np.zeros(n_cells, dtype=bool)
    for g in ["Cd3d", "Cd3e", "Nkg7", "Ptprc"]:
        if g in vecs:
            leak |= vecs[g] > 0
    smp["tnk_leak"] = leak

    rows = []
    for mouse, sub in smp.groupby("mouse"):
        geno = sub["genotype"].iloc[0]
        rows.append(
            {
                "accession": "GSE154989",
                "unit": "mouse",
                "mouse_id": mouse,
                "genotype": geno,
                "primary_arm": geno if geno in {"K", "KP"} else "exclude",
                "n_cells": int(len(sub)),
                "n_tumors_deposited": int(sub["mouseID"].nunique()),
                "cldn4": float(sub["cldn4"].mean()),
                "cldn4_pctpos": float((sub["cldn4"] > 0).mean()),
                "ifn_score": float(sub["ifn_score"].mean()),
                "mhc_score": float(sub["mhc_score"].mean()),
                "ifn_mhc_score": float(sub["ifn_mhc_score"].mean()),
                "tnk_score": np.nan,
                "tnk_fraction": np.nan,
                "n_tnk_leak_cells": int(sub["tnk_leak"].sum()),
            }
        )
    units = pd.DataFrame(rows)
    units["in_primary"] = (units["primary_arm"].isin(["K", "KP"])) & (units["n_cells"] >= 20)
    units.to_csv(OUT / "GSE154989.tsv", sep="\t", index=False)

    contrasts = []
    # No KL arm.
    for metric in ["cldn4", "cldn4_pctpos", "ifn_score", "mhc_score", "ifn_mhc_score"]:
        for cname, note in [
            (
                "KL_minus_KP",
                "no-go: GSE154989 is K/KP epithelium only; there is no KL / Stk11-loss arm.",
            ),
            (
                "KL_minus_K",
                "no-go: GSE154989 is K/KP epithelium only; there is no KL / Stk11-loss arm.",
            ),
        ]:
            contrasts.append(
                {
                    "accession": "GSE154989",
                    "contrast": cname,
                    "unit": "mouse",
                    "metric": metric,
                    "n_KL": 0,
                    "n_sensitive": int(units.loc[units["in_primary"] & (units["genotype"] == ("KP" if "KP" in cname else "K"))].shape[0]),
                    "mean_KL": np.nan,
                    "mean_sensitive": np.nan,
                    "delta_KL_minus_sensitive": np.nan,
                    "welch_p": np.nan,
                    "mwu_p": np.nan,
                    "rank_biserial": np.nan,
                    "status": "no-go",
                    "note": note,
                }
            )
    # Allowed: K vs KP Cldn4 (and IFN/MHC on the same epithelial cells).
    prim = units.loc[units["in_primary"]]
    for metric in ["cldn4", "cldn4_pctpos", "ifn_score", "mhc_score", "ifn_mhc_score"]:
        aa = prim.loc[prim["genotype"] == "KP", metric]
        bb = prim.loc[prim["genotype"] == "K", metric]
        st = welch_mwu(aa, bb)
        contrasts.append(
            {
                "accession": "GSE154989",
                "contrast": "KP_minus_K",
                "unit": "mouse",
                "metric": metric,
                "n_KL": st["n_a"],
                "n_sensitive": st["n_b"],
                "mean_KL": st["mean_a"],
                "mean_sensitive": st["mean_b"],
                "delta_KL_minus_sensitive": st["delta_a_minus_b"],
                "welch_p": st["welch_p"],
                "mwu_p": st["mwu_p"],
                "rank_biserial": st["rank_biserial"],
                "status": "ok_not_KL",
                "note": (
                    "Allowed K vs KP epithelial contrast (no KL). "
                    "n_KL column here is n_KP. T/NK is a design no-go (FACS CD45−)."
                ),
            }
        )
    # T/NK no-go rows
    for cname in ["KL_minus_KP", "KL_minus_K", "KP_minus_K"]:
        contrasts.append(
            {
                "accession": "GSE154989",
                "contrast": cname,
                "unit": "mouse",
                "metric": "tnk_fraction",
                "n_KL": 0,
                "n_sensitive": 0,
                "mean_KL": np.nan,
                "mean_sensitive": np.nan,
                "delta_KL_minus_sensitive": np.nan,
                "welch_p": np.nan,
                "mwu_p": np.nan,
                "rank_biserial": np.nan,
                "status": "no-go",
                "note": "T/NK design no-go: FACS tdTomato+/CD45−/CD11b−/TER119−/CD31−. Do not invent T/NK from leak Cd3d/Nkg7.",
            }
        )
    kp = prim.loc[prim["genotype"] == "KP"]
    r_ifn, p_ifn = spearman(kp["cldn4"], kp["ifn_score"])
    r_mhc, p_mhc = spearman(kp["cldn4"], kp["mhc_score"])
    contrasts.append(
        {
            "accession": "GSE154989",
            "contrast": "KP_only_Cldn4_vs_IFN",
            "unit": "mouse",
            "metric": "spearman_cldn4_ifn",
            "n_KL": 0,
            "n_sensitive": int(len(kp)),
            "mean_KL": np.nan,
            "mean_sensitive": np.nan,
            "delta_KL_minus_sensitive": r_ifn,
            "welch_p": np.nan,
            "mwu_p": p_ifn,
            "rank_biserial": np.nan,
            "status": "ok_not_KL",
            "note": f"Allowed within-KP Cldn4 vs IFN (Spearman ρ in delta column). n={len(kp)} KP mice ≥20 cells.",
        }
    )
    contrasts.append(
        {
            "accession": "GSE154989",
            "contrast": "KP_only_Cldn4_vs_MHC",
            "unit": "mouse",
            "metric": "spearman_cldn4_mhc",
            "n_KL": 0,
            "n_sensitive": int(len(kp)),
            "mean_KL": np.nan,
            "mean_sensitive": np.nan,
            "delta_KL_minus_sensitive": r_mhc,
            "welch_p": np.nan,
            "mwu_p": p_mhc,
            "rank_biserial": np.nan,
            "status": "ok_not_KL",
            "note": f"Allowed within-KP Cldn4 vs MHC (Spearman ρ in delta column). n={len(kp)} KP mice ≥20 cells.",
        }
    )
    pd.DataFrame(contrasts).to_csv(OUT / "GSE154989_contrasts.tsv", sep="\t", index=False)
    return {
        "accession": "GSE154989",
        "n_cells": n_cells,
        "n_mice": int(len(units)),
        "n_K_mice": int(((units["genotype"] == "K") & units["in_primary"]).sum()),
        "n_KP_mice": int(((units["genotype"] == "KP") & units["in_primary"]).sum()),
        "tnk_leak_cells": int(smp["tnk_leak"].sum()),
        "contrasts": contrasts,
        "units": units,
        "present": {
            "Cldn4": ["Cldn4"] if "Cldn4" in vecs else [],
            "T_NK": tnk_p,
            "IFN": ifn_p,
            "MHC": mhc_p,
        },
    }


# ---------------------------------------------------------------------------
# GSE179502
# ---------------------------------------------------------------------------

GSE179502_META = {
    "CM0875": {
        "cohort": "NonRestored",
        "treatment": "Tamoxifen",
        "restorable": False,
        "arm": "KL",
    },
    "CM0879": {
        "cohort": "Restored",
        "treatment": "Tamoxifen",
        "restorable": True,
        "arm": "Restored",
    },
    "CM0884": {
        "cohort": "Restored",
        "treatment": "Tamoxifen",
        "restorable": True,
        "arm": "Restored",
    },
    "ZR1932": {
        "cohort": "NonRestored",
        "treatment": "Vehicle",
        "restorable": False,
        "arm": "KL",
    },
    "ZR1966": {
        "cohort": "NonRestored",
        "treatment": "Vehicle",
        "restorable": True,
        "arm": "KL",
    },
    "ZR1969": {
        "cohort": "Restored",
        "treatment": "Tamoxifen",
        "restorable": True,
        "arm": "Restored",
    },
}


def _open_maybe_gz(path: Path):
    # GEO sometimes ships barcodes as uncompressed text with a .gz suffix.
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def run_gse179502() -> dict:
    from scipy.io import mmread

    feat_path = download("GSE179502_feat")
    bc_path = download("GSE179502_bc")
    mtx_path = download("GSE179502_mtx")
    with gzip.open(feat_path, "rt") as f:
        feats = [line.rstrip("\n").split("\t") for line in f]
    symbols = [r[1] for r in feats]
    with _open_maybe_gz(bc_path) as f:
        barcodes = [line.strip() for line in f if line.strip()]
    mice = [b.split("_")[0] for b in barcodes]
    print(f"GSE179502 loading mtx ({len(symbols)} x {len(barcodes)})")
    mat = mmread(str(mtx_path)).tocsr()
    if mat.shape[0] != len(symbols):
        mat = mat.T.tocsr()
    n_umi = np.asarray(mat.sum(axis=0)).ravel()
    n_genes = np.asarray((mat > 0).sum(axis=0)).ravel()
    mito_idx = [i for i, s in enumerate(symbols) if s.lower().startswith("mt-")]
    mito = np.asarray(mat[mito_idx, :].sum(axis=0)).ravel() / np.maximum(n_umi, 1)
    qc = (n_umi >= 500) & (n_genes >= 200) & (mito <= 0.20)

    needed = list(dict.fromkeys(CLDN4 + T_NK + IFN_MHC + ["Stk11", "Epcam"]))
    sym_to_i = {s: i for i, s in enumerate(symbols)}
    present_genes = [g for g in needed if g in sym_to_i]
    # log1p CP10k for needed genes
    scale = 1e4 / np.maximum(n_umi, 1)
    gene_log = {}
    for g in present_genes:
        raw = np.asarray(mat[sym_to_i[g], :].todense()).ravel()
        gene_log[g] = np.log1p(raw * scale)

    rows = []
    for mouse, meta in GSE179502_META.items():
        idx = np.array([m == mouse for m in mice])
        q = idx & qc
        n_qc = int(q.sum())
        cldn = gene_log.get("Cldn4", np.zeros(len(barcodes)))
        ifn_p = [g for g in IFN if g in gene_log]
        mhc_p = [g for g in MHC if g in gene_log]
        rows.append(
            {
                "accession": "GSE179502",
                "unit": "mouse",
                "mouse_id": mouse,
                "genotype": "KT;Lkb1XTR",
                "primary_arm": meta["arm"],
                "cohort": meta["cohort"],
                "treatment": meta["treatment"],
                "restorable": meta["restorable"],
                "n_barcodes": int(idx.sum()),
                "n_cells": n_qc,
                "cldn4": float(cldn[q].mean()) if n_qc else np.nan,
                "cldn4_pctpos": float((cldn[q] > 0).mean()) if n_qc else np.nan,
                "stk11": float(gene_log["Stk11"][q].mean()) if n_qc and "Stk11" in gene_log else np.nan,
                "ifn_score": float(np.mean([gene_log[g][q] for g in ifn_p], axis=0).mean()) if ifn_p and n_qc else np.nan,
                "mhc_score": float(np.mean([gene_log[g][q] for g in mhc_p], axis=0).mean()) if mhc_p and n_qc else np.nan,
                "ifn_mhc_score": float(np.mean([gene_log[g][q] for g in ifn_p + mhc_p], axis=0).mean())
                if (ifn_p + mhc_p) and n_qc
                else np.nan,
                "tnk_score": np.nan,
                "tnk_fraction": np.nan,
            }
        )
    units = pd.DataFrame(rows)
    units.to_csv(OUT / "GSE179502.tsv", sep="\t", index=False)

    contrasts = []
    kl = units.loc[units["primary_arm"] == "KL"]
    rest = units.loc[units["primary_arm"] == "Restored"]
    for metric in ["cldn4", "cldn4_pctpos", "stk11", "ifn_score", "mhc_score", "ifn_mhc_score"]:
        st = welch_mwu(kl[metric], rest[metric])
        contrasts.append(
            {
                "accession": "GSE179502",
                "contrast": "KL_minus_Restored",
                "unit": "mouse",
                "metric": metric,
                "n_KL": st["n_a"],
                "n_sensitive": st["n_b"],
                "mean_KL": st["mean_a"],
                "mean_sensitive": st["mean_b"],
                "delta_KL_minus_sensitive": st["delta_a_minus_b"],
                "welch_p": st["welch_p"],
                "mwu_p": st["mwu_p"],
                "rank_biserial": st["rank_biserial"],
                "status": "ok_proxy",
                "note": (
                    "NonRestored ≈ KL (Lkb1 off) vs Restored (Lkb1 on). "
                    "Restored is K-like (Kras; Lkb1 restored), not KP (no Trp53 arm). "
                    "n=3 vs 3 mice; MWU cannot beat p=0.1. Sorted neoplastic epithelium."
                ),
            }
        )
        contrasts.append(
            {
                "accession": "GSE179502",
                "contrast": "KL_minus_KP",
                "unit": "mouse",
                "metric": metric,
                "n_KL": 3,
                "n_sensitive": 0,
                "mean_KL": st["mean_a"],
                "mean_sensitive": np.nan,
                "delta_KL_minus_sensitive": np.nan,
                "welch_p": np.nan,
                "mwu_p": np.nan,
                "rank_biserial": np.nan,
                "status": "no-go",
                "note": "no-go for KL vs KP: series is KT;Lkb1 XTR only (no p53 / KP arm).",
            }
        )
        contrasts.append(
            {
                "accession": "GSE179502",
                "contrast": "KL_minus_K",
                "unit": "mouse",
                "metric": metric,
                "n_KL": st["n_a"],
                "n_sensitive": st["n_b"],
                "mean_KL": st["mean_a"],
                "mean_sensitive": st["mean_b"],
                "delta_KL_minus_sensitive": st["delta_a_minus_b"],
                "welch_p": st["welch_p"],
                "mwu_p": st["mwu_p"],
                "rank_biserial": st["rank_biserial"],
                "status": "ok_proxy",
                "note": "Proxy KL vs K: NonRestored minus Restored (Lkb1-on ≈ K). Same numbers as KL_minus_Restored.",
            }
        )
    for cname in ["KL_minus_KP", "KL_minus_K", "KL_minus_Restored"]:
        contrasts.append(
            {
                "accession": "GSE179502",
                "contrast": cname,
                "unit": "mouse",
                "metric": "tnk_fraction",
                "n_KL": 3 if cname != "KL_minus_KP" else 3,
                "n_sensitive": 0 if cname == "KL_minus_KP" else 3,
                "mean_KL": np.nan,
                "mean_sensitive": np.nan,
                "delta_KL_minus_sensitive": np.nan,
                "welch_p": np.nan,
                "mwu_p": np.nan,
                "rank_biserial": np.nan,
                "status": "no-go",
                "note": "T/NK design no-go: FACS-sorted neoplastic epithelium. Do not invent T/NK.",
            }
        )
    pd.DataFrame(contrasts).to_csv(OUT / "GSE179502_contrasts.tsv", sep="\t", index=False)
    return {
        "accession": "GSE179502",
        "n_mice": 6,
        "n_qc": int(qc.sum()),
        "contrasts": contrasts,
        "units": units,
        "present": {
            "Cldn4": ["Cldn4"] if "Cldn4" in gene_log else [],
            "T_NK": [g for g in T_NK if g in gene_log],
            "IFN": [g for g in IFN if g in gene_log],
            "MHC": [g for g in MHC if g in gene_log],
        },
    }


# ---------------------------------------------------------------------------
# GSE267321
# ---------------------------------------------------------------------------

def run_gse267321() -> dict:
    path = download("GSE267321_csv")
    print("GSE267321 streaming needed genes from CSV")
    needed = list(dict.fromkeys(CLDN4 + T_NK + IFN_MHC + EPI_MARKERS + HOST_LUNG + T_MARKERS + NK_MARKERS))
    # header = cells
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split(",")
        cells = header
        gene_rows = {}
        for line in f:
            gene, *vals = line.rstrip("\n").split(",")
            if gene in needed:
                gene_rows[gene] = np.array([float(x) if x else 0.0 for x in vals], dtype=np.float64)
    n = len(cells)
    meta = []
    for c in cells:
        # K_BARCODE.1_LKR13.K.1
        parts = c.split("_")
        geno = parts[0]
        sample = parts[-1].replace(".", "-")
        meta.append({"cell": c, "genotype": geno, "sample": sample})
    meta = pd.DataFrame(meta)
    expr = {g: gene_rows[g] if g in gene_rows else np.zeros(n) for g in needed}

    t_mask = np.zeros(n, dtype=bool)
    for g in T_MARKERS:
        t_mask |= expr[g] > 0
    nk_mask = np.zeros(n, dtype=bool)
    for g in NK_MARKERS:
        nk_mask |= expr[g] > 0
    tnk_mask = t_mask | nk_mask
    epi_mask = (expr.get("Epcam", np.zeros(n)) > 0) | (
        (expr.get("Cdh1", np.zeros(n)) > 0)
        & ((expr.get("Krt8", np.zeros(n)) > 0) | (expr.get("Krt18", np.zeros(n)) > 0) | (expr.get("Krt19", np.zeros(n)) > 0))
    )
    host_mask = np.zeros(n, dtype=bool)
    for g in HOST_LUNG:
        host_mask |= expr[g] > 0
    host_epi = host_mask & epi_mask
    marker_epi = epi_mask & ~host_epi
    cldn = expr["Cldn4"]
    ifn_p = [g for g in IFN if g in gene_rows]
    mhc_p = [g for g in MHC if g in gene_rows]

    rows = []
    for sample, idx in meta.groupby("sample").groups.items():
        idx = np.array(list(idx))
        geno = meta.loc[idx[0], "genotype"]
        n_cells = int(len(idx))
        n_tnk = int(tnk_mask[idx].sum())
        n_epi = int(marker_epi[idx].sum())
        epi_i = idx[marker_epi[idx]]
        rows.append(
            {
                "accession": "GSE267321",
                "unit": "tumor",
                "mouse_id": sample,
                "genotype": geno,
                "primary_arm": "KL" if geno == "KLK" else geno,
                "stk11": "KO" if geno == "KLK" else "WT",
                "keap1": "KO" if geno in {"KK", "KLK"} else "WT",
                "n_cells": n_cells,
                "n_tnk": n_tnk,
                "tnk_fraction": n_tnk / n_cells if n_cells else np.nan,
                "n_epithelial": n_epi,
                "n_host_epithelial": int(host_epi[idx].sum()),
                "cldn4": float(cldn[idx].mean()),
                "cldn4_pctpos": float((cldn[idx] > 0).mean()) * 100.0,
                "n_cldn4_pos": int((cldn[idx] > 0).sum()),
                "cldn4_epithelial": float(cldn[epi_i].mean()) if n_epi else np.nan,
                "ifn_score": float(np.mean([expr[g][epi_i] for g in ifn_p], axis=0).mean())
                if ifn_p and n_epi
                else np.nan,
                "mhc_score": float(np.mean([expr[g][epi_i] for g in mhc_p], axis=0).mean())
                if mhc_p and n_epi
                else np.nan,
                "ifn_mhc_score": float(np.mean([expr[g][epi_i] for g in ifn_p + mhc_p], axis=0).mean())
                if (ifn_p + mhc_p) and n_epi
                else np.nan,
                "tnk_score": np.nan,
            }
        )
    units = pd.DataFrame(rows)
    units.to_csv(OUT / "GSE267321.tsv", sep="\t", index=False)

    contrasts = []
    klk = units.loc[units["genotype"] == "KLK"]
    k = units.loc[units["genotype"] == "K"]
    kk = units.loc[units["genotype"] == "KK"]
    n_cldn4_pos_all = int((cldn > 0).sum())
    cldn4_empty = n_cldn4_pos_all <= 10
    for metric, cldn4_metric in [
        ("tnk_fraction", False),
        ("cldn4", True),
        ("cldn4_pctpos", True),
        ("cldn4_epithelial", True),
        ("ifn_score", False),
        ("mhc_score", False),
        ("ifn_mhc_score", False),
    ]:
        for a_df, b_df, cname, extra in [
            (klk, k, "KL_minus_K", "KLK (STK11+KEAP1) minus K. No STK11-only KL library."),
            (klk, k, "KLK_minus_K", "Same as KL_minus_K; KLK is the only Stk11-loss arm."),
            (klk, kk, "KLK_minus_KK", "KLK minus KEAP1-only KK."),
        ]:
            st = welch_mwu(a_df[metric], b_df[metric])
            status = "ok"
            note = extra + f" Honest n=2 vs 2 tumors; MWU cannot beat p=1/3."
            if cname.startswith("KL_minus") and cname == "KL_minus_K":
                pass
            if cldn4_metric and cldn4_empty:
                status = "no-go"
                note = (
                    f"Cldn4 essentially empty: {n_cldn4_pos_all}/7956 cells >0. "
                    "Non-malignant digest; do not treat Cldn4 as a genotype test."
                )
            if metric in {"ifn_score", "mhc_score", "ifn_mhc_score"}:
                status = "underpowered"
                note += " IFN/MHC scored on leftover marker epithelium (not author-malignant); n=2 vs 2."
            contrasts.append(
                {
                    "accession": "GSE267321",
                    "contrast": cname,
                    "unit": "tumor",
                    "metric": metric,
                    "n_KL": st["n_a"],
                    "n_sensitive": st["n_b"],
                    "mean_KL": st["mean_a"],
                    "mean_sensitive": st["mean_b"],
                    "delta_KL_minus_sensitive": st["delta_a_minus_b"],
                    "welch_p": st["welch_p"],
                    "mwu_p": st["mwu_p"],
                    "rank_biserial": st["rank_biserial"],
                    "status": status,
                    "note": note,
                }
            )
        contrasts.append(
            {
                "accession": "GSE267321",
                "contrast": "KL_minus_KP",
                "unit": "tumor",
                "metric": metric,
                "n_KL": int(len(klk)),
                "n_sensitive": 0,
                "mean_KL": float(klk[metric].mean()) if len(klk) else np.nan,
                "mean_sensitive": np.nan,
                "delta_KL_minus_sensitive": np.nan,
                "welch_p": np.nan,
                "mwu_p": np.nan,
                "rank_biserial": np.nan,
                "status": "no-go",
                "note": "no-go for KL vs KP: genotypes are K / KK / KLK only. No p53 / KP arm.",
            }
        )
    pd.DataFrame(contrasts).to_csv(OUT / "GSE267321_contrasts.tsv", sep="\t", index=False)
    return {
        "accession": "GSE267321",
        "n_cells": n,
        "n_cldn4_pos": n_cldn4_pos_all,
        "n_tnk": int(tnk_mask.sum()),
        "n_epi": int(marker_epi.sum()),
        "contrasts": contrasts,
        "units": units,
        "present": {
            "Cldn4": ["Cldn4"] if "Cldn4" in gene_rows else [],
            "T_NK": [g for g in T_NK if g in gene_rows],
            "IFN": ifn_p,
            "MHC": mhc_p,
        },
    }


def write_summary(results: list[dict]) -> pd.DataFrame:
    rows = []
    for res in results:
        for c in res["contrasts"]:
            if c["contrast"] not in {"KL_minus_KP", "KL_minus_K", "KL_minus_Restored", "KP_minus_K"}:
                continue
            if c["metric"] not in {
                "cldn4",
                "cldn4_pctpos",
                "tnk_score",
                "tnk_fraction",
                "ifn_score",
                "mhc_score",
                "ifn_mhc_score",
            }:
                continue
            rows.append(c)
    summary = pd.DataFrame(rows)
    # Prefer the honest biological unit per series.
    keep = []
    for acc, unit in {
        "GSE6135": "mouse",
        "GSE154989": "mouse",
        "GSE179502": "mouse",
        "GSE267321": "tumor",
    }.items():
        keep.append(summary[(summary["accession"] == acc) & (summary["unit"] == unit)])
    summary = pd.concat(keep, ignore_index=True)
    order = ["accession", "contrast", "metric", "status", "n_KL", "n_sensitive"]
    summary = summary.sort_values(["accession", "contrast", "metric"])
    summary.to_csv(OUT / "contrast_summary.tsv", sep="\t", index=False)

    genes = []
    for res in results:
        for set_name, lst in res.get("present", {}).items():
            genes.append(
                {
                    "accession": res["accession"],
                    "set": set_name,
                    "n_present": len(lst),
                    "present": ",".join(lst),
                }
            )
    pd.DataFrame(genes).to_csv(OUT / "genes_used.tsv", sep="\t", index=False)

    honest = pd.DataFrame(
        [
            {
                "accession": "GSE6135",
                "contrastable_KL_vs_KP": "yes",
                "contrastable_KL_vs_K": "yes",
                "honest_unit": "mouse (primary tumors collapsed)",
                "n_KL": results[0]["n_KL_mice"],
                "n_KP": results[0]["n_KP_mice"],
                "n_K": results[0]["n_K_mice"],
                "cldn4": "yes (array)",
                "tnk": "yes (bulk score)",
                "ifn_mhc": "yes (bulk score)",
                "note": "Only public KL bulk with Cldn4. Excluded L/+ and metastasis from primary KL.",
            },
            {
                "accession": "GSE154989",
                "contrastable_KL_vs_KP": "no-go",
                "contrastable_KL_vs_K": "no-go",
                "honest_unit": "biological mouse (strip _T#; ≥20 cells)",
                "n_KL": 0,
                "n_KP": results[1]["n_KP_mice"],
                "n_K": results[1]["n_K_mice"],
                "cldn4": "yes (epithelium)",
                "tnk": "no-go (CD45− FACS)",
                "ifn_mhc": "yes within KP / K vs KP",
                "note": "No KL arm. Allowed: K vs KP Cldn4; KP-only Cldn4 vs IFN/MHC.",
            },
            {
                "accession": "GSE179502",
                "contrastable_KL_vs_KP": "no-go",
                "contrastable_KL_vs_K": "yes (Restored proxy)",
                "honest_unit": "mouse",
                "n_KL": 3,
                "n_KP": 0,
                "n_K": 3,
                "cldn4": "yes (neoplastic epithelium)",
                "tnk": "no-go (sorted neoplastic)",
                "ifn_mhc": "yes (same epithelium)",
                "note": "NonRestored ≈ KL vs Restored ≈ K. n=3 vs 3.",
            },
            {
                "accession": "GSE267321",
                "contrastable_KL_vs_KP": "no-go",
                "contrastable_KL_vs_K": "yes (T/NK only)",
                "honest_unit": "tumor (2 per genotype)",
                "n_KL": 2,
                "n_KP": 0,
                "n_K": 2,
                "cldn4": "empty",
                "tnk": "yes (fraction)",
                "ifn_mhc": "leftover epi only; underpowered",
                "note": "KLK not KL; KEAP1 co-loss. Cldn4 floor. n=2 vs 2.",
            },
        ]
    )
    honest.to_csv(OUT / "honest_n.tsv", sep="\t", index=False)
    return summary


def main() -> None:
    print("=== GSE6135 ===")
    r6135 = run_gse6135()
    print(
        f"GSE6135 mice KL/KP/K = {r6135['n_KL_mice']}/{r6135['n_KP_mice']}/{r6135['n_K_mice']}"
    )
    print("=== GSE154989 ===")
    r154 = run_gse154989()
    print(f"GSE154989 K/KP mice ≥20 = {r154['n_K_mice']}/{r154['n_KP_mice']}")
    print("=== GSE179502 ===")
    r179 = run_gse179502()
    print(f"GSE179502 QC cells = {r179['n_qc']}")
    print("=== GSE267321 ===")
    r267 = run_gse267321()
    print(f"GSE267321 Cldn4+ cells = {r267['n_cldn4_pos']} / {r267['n_cells']}")
    results = [r6135, r154, r179, r267]
    summary = write_summary(results)
    with open(OUT / "summary.json", "w") as fh:
        json.dump(
            {
                "GSE6135": {k: v for k, v in r6135.items() if k not in {"contrasts", "units"}},
                "GSE154989": {k: v for k, v in r154.items() if k not in {"contrasts", "units"}},
                "GSE179502": {k: v for k, v in r179.items() if k not in {"contrasts", "units"}},
                "GSE267321": {k: v for k, v in r267.items() if k not in {"contrasts", "units"}},
            },
            fh,
            indent=2,
            default=str,
        )
    print("wrote", OUT)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
