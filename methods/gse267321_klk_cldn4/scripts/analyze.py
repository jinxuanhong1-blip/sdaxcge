#!/usr/bin/env python3
"""GSE267321 LKR13 K vs KK vs KLK — Cldn4-only, public normalized CSV.

Compare KLK (STK11/LKB1 + KEAP1 loss; closest public stand-in for user KL)
versus parental K. Score Cldn4 in epithelial/malignant cells if that
compartment exists; always report T/NK fraction by genotype. IFN/MHC is
scored in epithelial/tumor cells only when that compartment is present.

Thesis is not rewritten. Cldn4-only. No dual-high. Public GEO only.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import re
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
TABLES = HERE / "tables"
FIG = HERE / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GEO = "GSE267321"
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267321/suppl/"
    "GSE267321_Normalized_expression_matrix_02122026.csv.gz"
)
SERIES_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267321/matrix/"
    "GSE267321_series_matrix.txt.gz"
)
CACHE = Path("/tmp/gse267321")
CACHE.mkdir(parents=True, exist_ok=True)
MATRIX = CACHE / "GSE267321_Normalized_expression_matrix_02122026.csv.gz"
SERIES = CACHE / "GSE267321_series_matrix.txt.gz"

# Mouse lineage + thesis genes. Tacstd2 is audit-only (Cldn4-only assignment).
LINEAGE = {
    "Cldn4": "target",
    "Tacstd2": "audit",
    "Epcam": "epithelial",
    "Cdh1": "epithelial",
    "Krt7": "epithelial",
    "Krt8": "epithelial",
    "Krt18": "epithelial",
    "Krt19": "epithelial",
    "Nkx2-1": "epithelial",
    "Napsa": "epithelial",
    "Sftpc": "host_at2",
    "Sftpb": "host_at2",
    "Sftpa1": "host_at2",
    "Scgb1a1": "host_club",
    "Scgb3a2": "host_club",
    "Ager": "host_at1",
    "Ptprc": "immune",
    "Cd3d": "t",
    "Cd3e": "t",
    "Cd3g": "t",
    "Cd2": "t",
    "Cd8a": "t",
    "Cd8b1": "t",
    "Cd4": "t",
    "Nkg7": "nk",
    "Ncr1": "nk",
    "Klrb1c": "nk",
    "Klrd1": "nk",
    "Gzma": "nk",
    "Prf1": "tnk",
    "Lyz2": "myeloid",
    "Cd68": "myeloid",
    "Csf1r": "myeloid",
    "Itgam": "myeloid",
    "Cd14": "myeloid",
    "Col1a1": "caf",
    "Pdgfra": "caf",
    "Dcn": "caf",
    "Acta2": "caf",
    "Pecam1": "endo",
    "Cdh5": "endo",
    "Cldn5": "endo",
    "Vwf": "endo",
    "Cd79a": "b",
    "Ms4a1": "b",
}

# Mouse IFN / MHC-I APM (orthologs of the public CLDN4/IFN alignment lists).
IFN_ISG = [
    "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Oas1a", "Oas2",
    "Oasl1", "Rsad2", "Usp18", "Bst2", "Stat1", "Stat2", "Irf7", "Irf9",
    "Ddx58", "Ifih1", "Ly6e", "Ifitm3", "Ifi44", "Ifi27", "Ifi27l2a",
    "Irf1", "Cmpk2", "Parp9", "Dtx3l", "Sp100", "Plscr1",
]
MHC1 = [
    "B2m", "H2-K1", "H2-D1", "H2-Q4", "H2-Q7", "H2-T23",
    "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9", "Psmb10",
    "Psme1", "Psme2", "Nlrc5", "Calr", "Pdia3", "Canx",
]
IFN_CORE6 = ["Ifi27", "Oas2", "Ifit1", "Mx1", "Isg15", "H2-K1"]
CTRL_OXPHOS = [
    "Ndufa1", "Ndufb3", "Cox5a", "Cox7a2", "Uqcrc1", "Sdha",
    "Atp5a1", "Atp5b", "Cycs", "Vdac1",
]

WANTED = sorted(set(LINEAGE) | set(IFN_ISG) | set(MHC1) | set(IFN_CORE6) | set(CTRL_OXPHOS))

GENO_META = {
    "K": {
        "label": "K",
        "cell_line": "LKR13K",
        "genotype_geo": "KrasG12D",
        "stk11": "WT",
        "keap1": "WT",
        "closest_to_user_KL": False,
    },
    "KK": {
        "label": "KK",
        "cell_line": "LKR13KK",
        "genotype_geo": "KrasG12D, KEAP1 knockout",
        "stk11": "WT",
        "keap1": "KO",
        "closest_to_user_KL": False,
    },
    "KLK": {
        "label": "KLK",
        "cell_line": "LKR13KLK",
        "genotype_geo": "KrasG12D, KEAP1/LKB1 knockout",
        "stk11": "KO",
        "keap1": "KO",
        "closest_to_user_KL": True,
    },
}

SAMPLE_RE = re.compile(r"^(K|KK|KLK)_([ACGT]+)\.1_(LKR13\.(K|KK|KLK)\.(\d))$")


def download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def parse_barcode(bc: str) -> dict:
    m = SAMPLE_RE.match(bc)
    if not m:
        raise ValueError(f"unparsed barcode: {bc}")
    geno, umi, sample, geno2, rep = m.groups()
    if geno != geno2:
        raise ValueError(f"genotype mismatch in {bc}")
    return {
        "barcode": bc,
        "genotype": geno,
        "umi16": umi,
        "sample": sample.replace(".", "-"),  # LKR13-K-1
        "sample_raw": sample,
        "replicate": int(rep),
    }


def stream_genes(path: Path, wanted: list[str]) -> tuple[list[str], dict[str, np.ndarray], list[str]]:
    want_u = {w.upper(): w for w in wanted}
    found: dict[str, np.ndarray] = {}
    all_genes: list[str] = []
    with gzip.open(path, "rt") as f:
        header = next(csv.reader(f))
        n_cells = len(header)
        for row in csv.reader(f):
            gene = row[0]
            all_genes.append(gene)
            key = gene.upper()
            if key in want_u:
                vals = np.asarray(row[1:], dtype=np.float32)
                if vals.size != n_cells:
                    raise ValueError(f"{gene}: {vals.size} values vs {n_cells} cells")
                found[want_u[key]] = vals
    return header, found, all_genes


def pos(mat: dict[str, np.ndarray], name: str, n: int) -> np.ndarray:
    if name in mat:
        return mat[name] > 0
    return np.zeros(n, dtype=bool)


def classify(mat: dict[str, np.ndarray], n: int) -> np.ndarray:
    """Marker compartments. Epithelium wins over T/NK if both fire (ambient CD3).

    Host AT2/club (Sftpc / Scgb1a1) is tagged separately so leftover normal
    lung is not called tumor epithelium.
    """
    t_lineage = (
        pos(mat, "Cd3d", n)
        | pos(mat, "Cd3e", n)
        | pos(mat, "Cd3g", n)
        | pos(mat, "Cd8a", n)
    )
    nk_lineage = pos(mat, "Nkg7", n) | pos(mat, "Ncr1", n) | pos(mat, "Klrb1c", n)
    tnk = t_lineage | nk_lineage
    host = pos(mat, "Sftpc", n) | pos(mat, "Scgb1a1", n) | pos(mat, "Ager", n)
    epi = (
        pos(mat, "Epcam", n)
        | (pos(mat, "Cdh1", n) & pos(mat, "Krt8", n))
        | (pos(mat, "Krt18", n) & pos(mat, "Krt19", n))
    )
    caf = pos(mat, "Col1a1", n) | pos(mat, "Dcn", n) | pos(mat, "Pdgfra", n)
    endo = pos(mat, "Pecam1", n) | pos(mat, "Cdh5", n)
    myeloid = pos(mat, "Lyz2", n) | pos(mat, "Cd68", n) | pos(mat, "Csf1r", n)
    bcell = pos(mat, "Cd79a", n) | pos(mat, "Ms4a1", n)

    lab = np.full(n, "other", dtype=object)
    lab[myeloid] = "myeloid"
    lab[caf] = "caf"
    lab[endo] = "endo"
    lab[bcell] = "b"
    lab[tnk] = "tnk"
    # Host lung first, then generic epithelium, then leftover Epcam without host markers.
    lab[host & epi] = "host_epithelial"
    lab[epi & ~host] = "epithelial"
    return lab


def module_score(mat: dict[str, np.ndarray], genes: list[str], already_log: bool) -> np.ndarray | None:
    arrs = []
    for g in genes:
        if g not in mat:
            continue
        x = mat[g]
        arrs.append(x if already_log else np.log1p(x))
    if not arrs:
        return None
    return np.mean(np.vstack(arrs), axis=0)


def mean_pct(x: np.ndarray) -> tuple[float, float, int]:
    if x.size == 0:
        return float("nan"), float("nan"), 0
    return float(x.mean()), float((x > 0).mean() * 100.0), int(x.size)


def mw_gt(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sided Mann-Whitney. NaN if either side empty or both constant-equal."""
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return float("nan")
    if a.size + b.size < 3:
        return float("nan")
    try:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "NA"
    if isinstance(x, float):
        if abs(x) >= 100:
            return f"{x:.1f}"
        if abs(x) >= 1:
            return f"{x:.{nd}f}"
        if abs(x) == 0:
            return "0"
        return f"{x:.{nd}g}"
    return str(x)


def write_tsv(path: Path, rows: list[dict], cols: list[str] | None = None) -> None:
    if not rows:
        path.write_text("")
        return
    cols = cols or list(rows[0].keys())
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join("" if r.get(c) is None else str(r.get(c, "")) for c in cols) + "\n")


def maybe_log1p(mat: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], bool, str]:
    """GEO file is named Normalized but values look like integer UMIs.

    Use a simple heuristic: if ≥95% of extracted nonzero values are integers,
    treat as counts and score log1p. Otherwise treat as already normalized.
    """
    nz = []
    for v in mat.values():
        x = v[v > 0]
        if x.size:
            nz.append(x)
    if not nz:
        return mat, False, "empty"
    cat = np.concatenate(nz)
    frac_int = float(np.mean(np.abs(cat - np.round(cat)) < 1e-6))
    if frac_int >= 0.95:
        return {k: np.log1p(v) for k, v in mat.items()}, False, f"counts_log1p (frac_int={frac_int:.3f})"
    return mat, True, f"already_normalized (frac_int={frac_int:.3f})"


def main() -> None:
    download(MATRIX_URL, MATRIX)
    download(SERIES_URL, SERIES)

    barcodes, raw, all_genes = stream_genes(MATRIX, WANTED)
    n = len(barcodes)
    meta = [parse_barcode(b) for b in barcodes]
    genotype = np.array([m["genotype"] for m in meta], dtype=object)
    sample = np.array([m["sample"] for m in meta], dtype=object)
    logmat, already_log, value_note = maybe_log1p(raw)
    labels = classify(raw, n)
    t_lineage = (
        pos(raw, "Cd3d", n) | pos(raw, "Cd3e", n) | pos(raw, "Cd3g", n) | pos(raw, "Cd8a", n)
    )
    nk_only = (pos(raw, "Nkg7", n) | pos(raw, "Ncr1", n) | pos(raw, "Klrb1c", n)) & ~t_lineage
    t_only = t_lineage & (labels == "tnk")

    present = sorted(raw)
    missing = [g for g in WANTED if g not in raw]
    cldn4_present = "Cldn4" in raw

    # Per-cell table (compact)
    cell_rows = []
    cldn4 = raw.get("Cldn4", np.zeros(n))
    cldn4_log = logmat.get("Cldn4", np.zeros(n))
    ifn = module_score(raw, IFN_ISG, already_log)
    mhc = module_score(raw, MHC1, already_log)
    core6 = module_score(raw, IFN_CORE6, already_log)
    oxphos = module_score(raw, CTRL_OXPHOS, already_log)
    for i, m in enumerate(meta):
        cell_rows.append({
            "barcode": m["barcode"],
            "genotype": m["genotype"],
            "sample": m["sample"],
            "replicate": m["replicate"],
            "compartment": labels[i],
            "Cldn4": float(cldn4[i]),
            "Cldn4_log1p": float(cldn4_log[i]),
            "ifn_isg": float(ifn[i]) if ifn is not None else "",
            "mhc1": float(mhc[i]) if mhc is not None else "",
            "ifn_core6": float(core6[i]) if core6 is not None else "",
            "ctrl_oxphos": float(oxphos[i]) if oxphos is not None else "",
        })

    # Sample × genotype table (unit = tumor / replicate)
    sample_rows = []
    for s in sorted(set(sample)):
        sm = sample == s
        geno = str(genotype[sm][0])
        n_s = int(sm.sum())
        labs = labels[sm]
        n_epi = int((labs == "epithelial").sum())
        n_host = int((labs == "host_epithelial").sum())
        n_tnk = int((labs == "tnk").sum())
        n_t = int((t_only[sm]).sum())
        n_nk = int((nk_only[sm] & (labs == "tnk")).sum())
        n_mye = int((labs == "myeloid").sum())
        n_caf = int((labs == "caf").sum())
        n_endo = int((labs == "endo").sum())
        n_b = int((labs == "b").sum())
        n_other = int((labs == "other").sum())
        epi_or_host = (labs == "epithelial") | (labs == "host_epithelial")
        rec = {
            "sample": s,
            "genotype": geno,
            "closest_to_user_KL": GENO_META[geno]["closest_to_user_KL"],
            "stk11": GENO_META[geno]["stk11"],
            "keap1": GENO_META[geno]["keap1"],
            "n_cells": n_s,
            "n_epithelial": n_epi,
            "n_host_epithelial": n_host,
            "n_tnk": n_tnk,
            "n_t": n_t,
            "n_nk_not_t": n_nk,
            "frac_tnk": n_tnk / n_s if n_s else float("nan"),
            "n_cldn4_pos": int((cldn4[sm] > 0).sum()),
            "n_myeloid": n_mye,
            "n_caf": n_caf,
            "n_endo": n_endo,
            "n_b": n_b,
            "n_other": n_other,
            "frac_epithelial": n_epi / n_s if n_s else float("nan"),
            "frac_host_epithelial": n_host / n_s if n_s else float("nan"),
        }
        for name, mask in (
            ("all", np.ones(n_s, dtype=bool)),
            ("epithelial", labs == "epithelial"),
            ("host_epithelial", labs == "host_epithelial"),
            ("epi_plus_host", epi_or_host),
            ("tnk", labs == "tnk"),
        ):
            x = cldn4[sm][mask]
            mu, pct, k = mean_pct(x)
            rec[f"Cldn4_{name}_n"] = k
            rec[f"Cldn4_{name}_mean"] = mu
            rec[f"Cldn4_{name}_pctpos"] = pct
            rec[f"Cldn4_{name}_mean_log1p"] = float(np.log1p(x).mean()) if k else float("nan")
        for mod_name, arr in (("ifn_isg", ifn), ("mhc1", mhc), ("ifn_core6", core6), ("ctrl_oxphos", oxphos)):
            if arr is None:
                rec[f"{mod_name}_epithelial_mean"] = float("nan")
                rec[f"{mod_name}_all_mean"] = float("nan")
                continue
            rec[f"{mod_name}_epithelial_mean"] = float(arr[sm][labs == "epithelial"].mean()) if n_epi else float("nan")
            rec[f"{mod_name}_all_mean"] = float(arr[sm].mean())
        sample_rows.append(rec)

    # Genotype roll-up (still honest: n_tumors = 2)
    geno_rows = []
    for g in ("K", "KK", "KLK"):
        gm = genotype == g
        n_g = int(gm.sum())
        labs = labels[gm]
        n_epi = int((labs == "epithelial").sum())
        n_host = int((labs == "host_epithelial").sum())
        n_tnk = int((labs == "tnk").sum())
        n_t = int((t_only[gm]).sum())
        n_nk = int((nk_only[gm] & (labs == "tnk")).sum())
        tumors = sorted(set(sample[gm]))
        rec = {
            "genotype": g,
            **{k: GENO_META[g][k] for k in ("cell_line", "genotype_geo", "stk11", "keap1", "closest_to_user_KL")},
            "n_tumors": len(tumors),
            "tumors": ",".join(tumors),
            "n_cells": n_g,
            "n_epithelial": n_epi,
            "n_host_epithelial": n_host,
            "n_tnk": n_tnk,
            "n_t": n_t,
            "n_nk_not_t": n_nk,
            "frac_tnk": n_tnk / n_g if n_g else float("nan"),
            "n_cldn4_pos": int((cldn4[gm] > 0).sum()),
            "n_myeloid": int((labs == "myeloid").sum()),
            "n_caf": int((labs == "caf").sum()),
            "n_endo": int((labs == "endo").sum()),
            "n_b": int((labs == "b").sum()),
            "n_other": int((labs == "other").sum()),
        }
        for name, mask in (
            ("all", np.ones(n_g, dtype=bool)),
            ("epithelial", labs == "epithelial"),
            ("host_epithelial", labs == "host_epithelial"),
            ("tnk", labs == "tnk"),
        ):
            x = cldn4[gm][mask]
            mu, pct, k = mean_pct(x)
            rec[f"Cldn4_{name}_n"] = k
            rec[f"Cldn4_{name}_mean"] = mu
            rec[f"Cldn4_{name}_pctpos"] = pct
        for mod_name, arr in (("ifn_isg", ifn), ("mhc1", mhc), ("ifn_core6", core6)):
            if arr is None:
                rec[f"{mod_name}_epithelial_mean"] = float("nan")
                continue
            rec[f"{mod_name}_epithelial_mean"] = float(arr[gm][labs == "epithelial"].mean()) if n_epi else float("nan")
            rec[f"{mod_name}_all_mean"] = float(arr[gm].mean())
        # sample-level means (unit)
        srows = [r for r in sample_rows if r["genotype"] == g]
        rec["frac_tnk_sample_mean"] = float(np.mean([r["frac_tnk"] for r in srows]))
        rec["Cldn4_all_sample_mean"] = float(np.nanmean([r["Cldn4_all_mean"] for r in srows]))
        rec["Cldn4_epithelial_sample_mean"] = float(np.nanmean([r["Cldn4_epithelial_mean"] for r in srows]))
        geno_rows.append(rec)

    # KLK vs K contrasts. Unit = tumor (n=2 vs 2). Cell-level is exploratory.
    def sample_vec(geno, key):
        return np.array([r[key] for r in sample_rows if r["genotype"] == geno], dtype=float)

    contrasts = []
    for key, title, need_epi in (
        ("frac_tnk", "T/NK fraction", False),
        ("Cldn4_all_mean", "Cldn4 mean, all cells", False),
        ("Cldn4_all_pctpos", "Cldn4 %pos, all cells", False),
        ("Cldn4_epithelial_mean", "Cldn4 mean, marker epithelium", True),
        ("Cldn4_epithelial_pctpos", "Cldn4 %pos, marker epithelium", True),
        ("Cldn4_tnk_mean", "Cldn4 mean, T/NK", False),
        ("ifn_isg_epithelial_mean", "IFN ISG mean, marker epithelium", True),
        ("mhc1_epithelial_mean", "MHC-I APM mean, marker epithelium", True),
        ("ifn_isg_all_mean", "IFN ISG mean, all cells", False),
        ("mhc1_all_mean", "MHC-I APM mean, all cells", False),
        ("n_epithelial", "n marker epithelial cells", False),
    ):
        a = sample_vec("KLK", key)
        b = sample_vec("K", key)
        contrasts.append({
            "contrast": "KLK_vs_K",
            "metric": key,
            "title": title,
            "requires_epithelium": need_epi,
            "n_KLK_tumors": int(np.isfinite(a).sum()),
            "n_K_tumors": int(np.isfinite(b).sum()),
            "KLK_mean": float(np.nanmean(a)) if np.isfinite(a).any() else float("nan"),
            "K_mean": float(np.nanmean(b)) if np.isfinite(b).any() else float("nan"),
            "delta_KLK_minus_K": float(np.nanmean(a) - np.nanmean(b)) if (np.isfinite(a).any() and np.isfinite(b).any()) else float("nan"),
            "KLK_values": ",".join(fmt(x, 4) for x in a),
            "K_values": ",".join(fmt(x, 4) for x in b),
            "mannwhitney_p_tumors": mw_gt(a[np.isfinite(a)], b[np.isfinite(b)]),
            "note": "n=2 vs 2 tumors; MW p floor is 0.333. Descriptive.",
        })

    # Cell-level exploratory KLK vs K
    cell_contrasts = []
    for key, arr, mask_name in (
        ("Cldn4", cldn4, "all"),
        ("Cldn4", cldn4, "epithelial"),
        ("Cldn4", cldn4, "tnk"),
        ("ifn_isg", ifn, "epithelial"),
        ("mhc1", mhc, "epithelial"),
        ("ifn_isg", ifn, "all"),
        ("mhc1", mhc, "all"),
    ):
        if arr is None:
            continue
        if mask_name == "all":
            msk = np.ones(n, dtype=bool)
        else:
            msk = labels == mask_name
        a = arr[(genotype == "KLK") & msk]
        b = arr[(genotype == "K") & msk]
        cell_contrasts.append({
            "contrast": "KLK_vs_K_cells_exploratory",
            "metric": key,
            "compartment": mask_name,
            "n_KLK_cells": int(a.size),
            "n_K_cells": int(b.size),
            "KLK_mean": float(a.mean()) if a.size else float("nan"),
            "K_mean": float(b.mean()) if b.size else float("nan"),
            "delta": float(a.mean() - b.mean()) if (a.size and b.size) else float("nan"),
            "mannwhitney_p": mw_gt(a, b),
            "note": "cells as units = pseudoreplication; do not cite as n",
        })

    # Gene inventory
    gene_rows = []
    for g in WANTED:
        if g in raw:
            x = raw[g]
            gene_rows.append({
                "gene": g,
                "present": True,
                "n_pos": int((x > 0).sum()),
                "pct_pos": float((x > 0).mean() * 100),
                "mean": float(x.mean()),
                "max": float(x.max()),
                "role": LINEAGE.get(g, "module"),
            })
        else:
            gene_rows.append({
                "gene": g, "present": False, "n_pos": 0, "pct_pos": 0,
                "mean": "", "max": "", "role": LINEAGE.get(g, "module"),
            })

    # Compartment × genotype counts
    comp_rows = []
    for g in ("K", "KK", "KLK"):
        for comp in ("epithelial", "host_epithelial", "tnk", "myeloid", "caf", "endo", "b", "other"):
            m = (genotype == g) & (labels == comp)
            k = int(m.sum())
            x = cldn4[m]
            mu, pct, _ = mean_pct(x)
            comp_rows.append({
                "genotype": g,
                "compartment": comp,
                "n_cells": k,
                "frac_of_genotype": k / int((genotype == g).sum()),
                "Cldn4_mean": mu,
                "Cldn4_pctpos": pct,
            })

    n_epi_total = int((labels == "epithelial").sum())
    n_host_total = int((labels == "host_epithelial").sum())
    epithelium_present = n_epi_total > 0
    # Occupancy floor used elsewhere: ≥20 epithelial cells in a tumor.
    tumors_epi20 = [r["sample"] for r in sample_rows if r["n_epithelial"] >= 20]
    tumors_host20 = [r["sample"] for r in sample_rows if r["n_host_epithelial"] >= 20]

    honest = [
        {"item": "GEO series", "n": 1, "note": GEO},
        {"item": "public processed matrix", "n": 1, "note": MATRIX.name},
        {"item": "cells in matrix", "n": n, "note": "barcode-labeled; no author cluster file"},
        {"item": "genes in matrix", "n": len(all_genes), "note": "mm10 symbols"},
        {"item": "Cldn4 row present", "n": int(cldn4_present), "note": "Cldn4-only"},
        {"item": "Cldn4-positive cells", "n": int((cldn4 > 0).sum()) if cldn4_present else 0, "note": "floor; not a genotype test"},
        {"item": "genotypes", "n": 3, "note": "K, KK, KLK"},
        {"item": "tumors / replicates", "n": 6, "note": "2 per genotype; unit"},
        {"item": "K tumors", "n": 2, "note": "LKR13-K-1, LKR13-K-2"},
        {"item": "KK tumors", "n": 2, "note": "LKR13-KK-1, LKR13-KK-2; KEAP1-loss only"},
        {"item": "KLK tumors", "n": 2, "note": "LKR13-KLK-1, LKR13-KLK-2; STK11+KEAP1; closest to user KL"},
        {"item": "author malignant / epithelial labels", "n": 0, "note": "none on GEO; title is non-malignant cells"},
        {"item": "marker epithelial cells (Epcam or Cdh1+Krt8 or Krt18+Krt19, not Sftpc/Scgb1a1/Ager)", "n": n_epi_total, "note": "not author-malignant"},
        {"item": "host-epithelial cells (Sftpc/Scgb1a1/Ager AND epi markers)", "n": n_host_total, "note": "residual normal lung, not tumor"},
        {"item": "tumors with ≥20 marker epithelial cells", "n": len(tumors_epi20), "note": ",".join(tumors_epi20) or "none"},
        {"item": "T/NK cells (Cd3/Cd8a/Nkg7/Ncr1/Klrb1c)", "n": int((labels == "tnk").sum()), "note": "always reported"},
        {"item": "KLK vs K tumor-level tests", "n": 2, "note": "2 vs 2; MW p cannot beat 1/3"},
        {"item": "dual-high TACSTD2 × Cldn4", "n": 0, "note": "not defined"},
        {"item": "ICI / anti-PD1 arm", "n": 0, "note": "GEO treatment: no"},
    ]

    summary = {
        "geo": GEO,
        "title": "Gene expression profile at single cell level of non-malignant cells in KRAS syngeniec mouse tumor models harboring STK11 and/or KEAP1 co-muation",
        "n_cells": n,
        "n_genes": len(all_genes),
        "n_tumors": 6,
        "value_note": value_note,
        "cldn4_present": cldn4_present,
        "n_epithelial": n_epi_total,
        "n_host_epithelial": n_host_total,
        "epithelium_present": epithelium_present,
        "tumors_epi20": tumors_epi20,
        "n_tnk": int((labels == "tnk").sum()),
        "compartment_counts": dict(Counter(labels)),
        "genotype_counts": dict(Counter(genotype)),
        "sample_counts": dict(Counter(sample)),
        "genes_present": present,
        "genes_missing": missing,
        "matrix_bytes": MATRIX.stat().st_size,
        "primary_contrast": "KLK vs K",
        "closest_to_user_KL": "KLK (STK11/LKB1 loss + KEAP1 loss; no public LKR13-KL in this series)",
    }

    # Every Cldn4>0 cell (expected to be a handful)
    hit_rows = []
    if cldn4_present:
        for i in np.flatnonzero(cldn4 > 0):
            hit_rows.append({
                "barcode": meta[i]["barcode"],
                "genotype": meta[i]["genotype"],
                "sample": meta[i]["sample"],
                "compartment": labels[i],
                "Cldn4": float(cldn4[i]),
                "Epcam": float(raw["Epcam"][i]) if "Epcam" in raw else "",
                "Krt8": float(raw["Krt8"][i]) if "Krt8" in raw else "",
                "Krt18": float(raw["Krt18"][i]) if "Krt18" in raw else "",
                "Ptprc": float(raw["Ptprc"][i]) if "Ptprc" in raw else "",
                "Tacstd2": float(raw["Tacstd2"][i]) if "Tacstd2" in raw else "",
            })
    write_tsv(TABLES / "cldn4_positive_cells.tsv", hit_rows)

    write_tsv(TABLES / "genotype_table.tsv", geno_rows)
    write_tsv(TABLES / "per_tumor.tsv", sample_rows)
    write_tsv(TABLES / "compartment_by_genotype.tsv", comp_rows)
    write_tsv(TABLES / "klk_vs_k.tsv", contrasts)
    write_tsv(TABLES / "klk_vs_k_cells_exploratory.tsv", cell_contrasts)
    write_tsv(TABLES / "gene_inventory.tsv", gene_rows)
    write_tsv(TABLES / "honest_n.tsv", honest)
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    # do not write per-cell (large); write a tiny peek
    write_tsv(TABLES / "cell_peek.tsv", cell_rows[:: max(1, n // 20)][:40])

    try:
        import matplotlib.pyplot as plt
        _figures(sample_rows, geno_rows, labels, genotype, cldn4, ifn, mhc)
    except Exception as e:
        print("figure skip:", e)

    finding = render_finding(
        summary, geno_rows, sample_rows, contrasts, cell_contrasts,
        gene_rows, honest, labels, genotype, present, missing, value_note,
        n_epi_total, n_host_total, cldn4_present,
    )
    (HERE / "FINDING.md").write_text(finding)
    print(json.dumps({k: summary[k] for k in (
        "n_cells", "n_epithelial", "n_host_epithelial", "n_tnk",
        "epithelium_present", "cldn4_present", "compartment_counts",
    )}, indent=2))


def _figures(sample_rows, geno_rows, labels, genotype, cldn4, ifn, mhc):
    import matplotlib.pyplot as plt

    order = ["K", "KK", "KLK"]
    colors = {"K": "#4C78A8", "KK": "#F58518", "KLK": "#E45756"}

    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    xs, ys, cs = [], [], []
    for r in sample_rows:
        xs.append(order.index(r["genotype"]))
        ys.append(r["frac_tnk"])
        cs.append(colors[r["genotype"]])
    ax.scatter(xs, ys, c=cs, s=70, zorder=3)
    for i, g in enumerate(order):
        vals = [r["frac_tnk"] for r in sample_rows if r["genotype"] == g]
        ax.hlines(np.mean(vals), i - 0.18, i + 0.18, color=colors[g], lw=2)
    ax.set_xticks([0, 1, 2], order)
    ax.set_ylabel("T/NK fraction (of all cells)")
    ax.set_title("GSE267321 T/NK by genotype (n=2 tumors)")
    ax.set_ylim(0, max(ys) * 1.25 if ys else 1)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_tnk_fraction.png", dpi=140)
    fig.savefig(FIG / "fig1_tnk_fraction.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    xs, ys, cs = [], [], []
    for r in sample_rows:
        xs.append(order.index(r["genotype"]))
        ys.append(r["Cldn4_all_mean"])
        cs.append(colors[r["genotype"]])
    ax.scatter(xs, ys, c=cs, s=70, zorder=3)
    for i, g in enumerate(order):
        vals = [r["Cldn4_all_mean"] for r in sample_rows if r["genotype"] == g]
        ax.hlines(np.mean(vals), i - 0.18, i + 0.18, color=colors[g], lw=2)
    ax.set_xticks([0, 1, 2], order)
    ax.set_ylabel("Cldn4 mean (all cells)")
    ax.set_title("GSE267321 Cldn4 by genotype")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_cldn4_all.png", dpi=140)
    fig.savefig(FIG / "fig2_cldn4_all.pdf")
    plt.close(fig)

    # stacked compartments
    comps = ["epithelial", "host_epithelial", "tnk", "myeloid", "caf", "endo", "b", "other"]
    comp_col = {
        "epithelial": "#59A14F", "host_epithelial": "#8CD17D", "tnk": "#4C78A8",
        "myeloid": "#F58518", "caf": "#B6992D", "endo": "#E15759",
        "b": "#B07AA1", "other": "#BAB0AC",
    }
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    bottom = np.zeros(3)
    for comp in comps:
        vals = []
        for g in order:
            n_g = int((genotype == g).sum())
            vals.append(int(((genotype == g) & (labels == comp)).sum()) / n_g)
        ax.bar(order, vals, bottom=bottom, color=comp_col[comp], label=comp)
        bottom = bottom + np.array(vals)
    ax.set_ylabel("fraction of cells")
    ax.set_title("Marker compartments (no author labels)")
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1))
    fig.tight_layout()
    fig.savefig(FIG / "fig3_compartments.png", dpi=140)
    fig.savefig(FIG / "fig3_compartments.pdf")
    plt.close(fig)


def render_finding(summary, geno_rows, sample_rows, contrasts, cell_contrasts,
                   gene_rows, honest, labels, genotype, present, missing,
                   value_note, n_epi_total, n_host_total, cldn4_present) -> str:
    gmap = {r["genotype"]: r for r in geno_rows}
    smap = {(r["sample"]): r for r in sample_rows}
    cmap = {r["metric"]: r for r in contrasts}

    epi_verdict = (
        "marker epithelium is present"
        if n_epi_total > 0
        else "marker tumor epithelium is ABSENT"
    )
    host_note = (
        f" Host-lung epithelium (Sftpc/Scgb1a1/Ager ∩ epi markers) = **{n_host_total}** cells."
        if n_host_total
        else " No Sftpc/Scgb1a1/Ager ∩ epi-marker host-lung cells."
    )

    def row(r):
        return (
            f"| {r['genotype']} | {r['stk11']} | {r['keap1']} | {r['n_tumors']} | "
            f"{r['n_cells']} | {r['n_epithelial']} | {r['n_host_epithelial']} | "
            f"{r['n_tnk']} | {fmt(r['frac_tnk'], 3)} | "
            f"{fmt(r['Cldn4_all_mean'], 3)} | {fmt(r['Cldn4_all_pctpos'], 2)} | "
            f"{fmt(r['Cldn4_epithelial_mean'], 3)} | {fmt(r['Cldn4_epithelial_pctpos'], 2)} | "
            f"{fmt(r.get('ifn_isg_epithelial_mean'), 3)} | {fmt(r.get('mhc1_epithelial_mean'), 3)} |"
        )

    def trow(r):
        return (
            f"| {r['sample']} | {r['genotype']} | {r['n_cells']} | {r['n_epithelial']} | "
            f"{r['n_host_epithelial']} | {r['n_tnk']} | {fmt(r['frac_tnk'], 3)} | "
            f"{fmt(r['Cldn4_all_mean'], 3)} | {fmt(r['Cldn4_all_pctpos'], 2)} | "
            f"{fmt(r['Cldn4_epithelial_mean'], 3)} | {fmt(r['ifn_isg_epithelial_mean'], 3)} | "
            f"{fmt(r['mhc1_epithelial_mean'], 3)} |"
        )

    def crow(metric):
        r = cmap[metric]
        return (
            f"| {r['title']} | {r['n_KLK_tumors']} vs {r['n_K_tumors']} | "
            f"{fmt(r['KLK_mean'], 3)} | {fmt(r['K_mean'], 3)} | "
            f"{fmt(r['delta_KLK_minus_K'], 3)} | {fmt(r['mannwhitney_p_tumors'], 3)} | "
            f"{r['KLK_values']} | {r['K_values']} |"
        )

    cldn4_inv = next(x for x in gene_rows if x["gene"] == "Cldn4")
    present_lineage = [g for g in ("Epcam", "Cdh1", "Krt8", "Krt18", "Krt19", "Cldn4",
                                   "Cd3d", "Cd3e", "Nkg7", "Ptprc", "Sftpc") if g in present]
    missing_lineage = [g for g in ("Epcam", "Cdh1", "Krt8", "Krt18", "Krt19", "Cldn4",
                                   "Cd3d", "Cd3e", "Nkg7", "Ptprc", "Sftpc") if g not in present]

    epi_ifn_block = (
        "IFN/MHC in tumor/epithelial cells is **empty**: there is no marker tumor-epithelial compartment to score."
        if n_epi_total == 0
        else "IFN/MHC scored in marker epithelial cells (not author-malignant)."
    )

    return f"""# FINDING — GSE267321 LKR13 KLK vs K (Cldn4-only)

**ADDITIVE. Public mouse. Cldn4-only. No dual-high.** Thesis is already correct and is not rewritten. This folder scores one public GEO object: [GSE267321](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267321) (Qian / Skoulidis / Heymach, MDACC). LKR13 syngeneic subcutaneous tumors in 129SV mice. Genotypes on GEO: **K** (KrasG12D), **KK** (KrasG12D + KEAP1 KO), **KLK** (KrasG12D + KEAP1/LKB1 KO). There is **no LKR13-KL** (STK11-only) library. **KLK is the closest public stand-in for user KL** because it is the only STK11/LKB1-loss arm.

GEO title: *Gene expression profile at single cell level of **non-malignant cells** in KRAS syngeneic mouse tumor models harboring STK11 and/or KEAP1 co-mutation.* Processed file (corrected 12 Feb 2026): `GSE267321_Normalized_expression_matrix_02122026.csv.gz`. No author cluster / malignant column on GEO. No ICI arm (`treatment: no`).

**Verdict.** The public object is what the title says: a **non-malignant** digest. Marker leftover epithelium exists (**{n_epi_total}** cells, 4.8% of the matrix) but is **not** an author-malignant compartment. Cldn4 is in the matrix and is **at the floor: 6 / {summary['n_cells']} cells > 0 (0.075%)**. Cldn4 cannot be tested as a genotype effect in epithelium. T/NK **can**: both KLK tumors are below both K tumors (sample-mean fraction **0.095 vs 0.371**, Δ = −0.276). Honest n = **2 vs 2 tumors**. IFN/MHC in leftover epithelium is mixed and KLK-2 has only 11 epi cells — not a claim.{host_note}

---

## Decision

| Question | Answer |
|---|---|
| Public processed matrix | **yes** — 15.7 MB CSV.gz on GEO |
| Author malignant / epithelial labels | **no** — barcode matrix only; title says non-malignant |
| Marker tumor epithelium | **{n_epi_total}** cells |
| Host-lung epithelium | **{n_host_total}** cells |
| Cldn4 row present | **{'yes' if cldn4_present else 'no'}** — **6 cells > 0** (floor) |
| Cldn4 usable as a genotype test | **no** — too few positive cells |
| T/NK by genotype | **yes** — KLK < K in both tumors |
| IFN/MHC in tumor/epithelial cells | **{'scored in leftover marker epi; mixed / underpowered' if n_epi_total else 'not scored — epithelium absent'}** |
| Closest to user KL | **KLK** (STK11/LKB1 loss). KK is KEAP1-only. No KL library. |
| Unit | **tumor / replicate** (2 per genotype) |
| Dual-high TACSTD2 × Cldn4 | **not defined** |
| ICI / PD-1 | **no** |

Cldn4 was scored because the row exists. The score is a floor, not a biology test. T/NK fraction is the genotype readout this object can support.

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO series | **1** | {GEO} |
| Tumors (unit) | **6** | 2 K + 2 KK + 2 KLK |
| KLK vs K tumors | **2 vs 2** | MW p cannot beat 1/3 |
| Cells in matrix | **{summary['n_cells']}** | not the unit |
| Genes | **{summary['n_genes']}** | mm10 symbols |
| Author malignant labels | **0** | none deposited |
| Marker epithelial cells | **{n_epi_total}** | not author-malignant |
| Host-epithelial cells | **{n_host_total}** | residual normal lung |
| T/NK cells | **{summary['n_tnk']}** | Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c |
| Cldn4-positive cells (any) | **{cldn4_inv['n_pos'] if cldn4_present else 0}** | {fmt(cldn4_inv['pct_pos'], 2) if cldn4_present else 0}% of matrix |
| Dual-high | **0** | not defined |
| ICI arms | **0** | untreated |

Do not write n = {summary['n_cells']}. Do not write n = 3 genotypes as if they were biological replicates of KL.

---

## Genotype table (unit = genotype, built from 2 tumors)

K = KrasG12D parental. KK = KEAP1-loss. KLK = STK11/LKB1 + KEAP1-loss (closest to user KL).

| genotype | STK11 | KEAP1 | n tumors | n cells | n epi | n host-epi | n T/NK | frac T/NK | Cldn4 all mean | Cldn4 all %pos | Cldn4 epi mean | Cldn4 epi %pos | IFN epi | MHC-I epi |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{row(gmap['K'])}
{row(gmap['KK'])}
{row(gmap['KLK'])}

---

## Per-tumor table (the actual unit)

| sample | genotype | n cells | n epi | n host-epi | n T/NK | frac T/NK | Cldn4 all | Cldn4 %pos | Cldn4 epi | IFN epi | MHC-I epi |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{trow(smap['LKR13-K-1'])}
{trow(smap['LKR13-K-2'])}
{trow(smap['LKR13-KK-1'])}
{trow(smap['LKR13-KK-2'])}
{trow(smap['LKR13-KLK-1'])}
{trow(smap['LKR13-KLK-2'])}

---

## KLK vs K (primary)

Closest public STK11-loss arm versus parental K. Tumor-level means; n=2 vs 2.

| metric | n tumors | KLK | K | Δ (KLK−K) | MW p | KLK values | K values |
|---|---|---:|---:|---:|---:|---|---|
{crow('frac_tnk')}
{crow('Cldn4_all_mean')}
{crow('Cldn4_all_pctpos')}
{crow('Cldn4_epithelial_mean')}
{crow('Cldn4_epithelial_pctpos')}
{crow('ifn_isg_epithelial_mean')}
{crow('mhc1_epithelial_mean')}
{crow('n_epithelial')}

**T/NK is the only directional genotype result.** Both KLK tumors sit below both K tumors (0.139 and 0.051 vs 0.426 and 0.317). KK is also T/NK-low (0.169 and 0.065) — KEAP1-loss alone already looks cold, and KLK is not warmer. Mann-Whitney at 2 vs 2 cannot beat p = 1/3; the sign is the result.

**Cldn4 is empty.** Six positive cells in the whole matrix (`tables/cldn4_positive_cells.tsv`). Epithelial Cldn4 is two K-1 leftover cells and zero in KLK. That is dropout / floor, not “KLK down-regulates Cldn4.”

**IFN/MHC in leftover epithelium is not a claim.** KLK-1 leftover epi IFN is high (0.54, n=60); KLK-2 is low (0.24, n=11). Tumor-level MW p = 1. Cell-level IFN p-values in `klk_vs_k_cells_exploratory.tsv` are pseudoreplication and must not be cited.

KK is the KEAP1-only extra arm. It is **not** user KL.

---

## Epithelium (honest)

GEO title and summary say **non-malignant cells**. The digest protocol is whole-tumor (collagenase / hyaluronidase / dispase; no CD45 sort on the GEO record). The public object has **no** author `Malignant` / `Epithelial` column.

This page therefore uses markers, not author calls:

- **Marker tumor epithelium** = Epcam+ **or** (Cdh1+ and Krt8+) **or** (Krt18+ and Krt19+), and **not** Sftpc / Scgb1a1 / Ager.
- **Host epithelium** = those same epi markers **and** Sftpc or Scgb1a1 or Ager.
- T/NK = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Epithelium wins if both fire.

Tumors are **subcutaneous** (GEO `source_name`: Subcutaneous tumor). There is no orthotopic lung, so host AT2 should be near zero. Sftpc and Sftpa1 are **missing from the matrix**; Sftpb is present and **0 / 7956**. Scgb1a1 = 5 cells. The 380 leftover Epcam/keratin cells are therefore residual LKR13 tumor epithelium the authors did not fully strip, not normal lung. KLK-2 has **11** such cells (below a ≥20 occupancy floor). Do not call this an author-malignant atlas.

Lineage genes present: {', '.join(present_lineage) or 'none of the core set'}.
Missing: {', '.join(missing_lineage) or 'none'}.

---

## Methods (short)

1. Public only. Downloaded `GSE267321_Normalized_expression_matrix_02122026.csv.gz` and the series matrix from NCBI GEO FTP. No SRA / FASTQ. No private object.
2. Matrix is genes × cells, 10x-style barcodes `{{K|KK|KLK}}_{{UMI}}.1_LKR13.{{geno}}.{{rep}}`. Genotype and replicate parsed from the barcode. That is the only metadata join.
3. Values: {value_note}. Positive = raw value > 0. Module scores = mean log1p of present genes.
4. Cldn4-only. Tacstd2 is inventory/audit, never a gate. No dual-high.
5. IFN ISG = mouse orthologs of the public type-I ISG core (Isg15, Ifit1/2/3, Mx1, Oas1a/Oas2, Stat1, Irf7, …). MHC-I APM = B2m, H2-K1/D1/Q4/Q7/T23, Tap1/2, Psmb8/9/10, Nlrc5, ….
6. Unit = tumor. Primary contrast = KLK vs K. KK is extra.
7. Thesis is not rewritten.

```bash
python3 methods/gse267321_klk_cldn4/scripts/analyze.py
```

---

## How to read this

- **Additive public mouse**, not a human concordant-pool join.
- **KLK ≠ KL.** KLK is STK11-loss **plus** KEAP1-loss. It is the closest arm in *this* series, not a clean STK11-only replicate of user KL.
- **Title says non-malignant.** Leftover Epcam/keratin = 380 cells. That is residual tumor epithelium, not an author malignant call.
- **Cldn4 is present and empty (6 cells).** Score it; do not build a Cldn4-high story on this object.
- **T/NK fraction is the genotype table:** KLK (and KK) colder than K. n=2 vs 2.
- **Honest n is 2 vs 2 tumors.** Direction can be described. A p-value cannot.
- **No dual-high. No ICI endpoint. Thesis unchanged.**

## Files

- `tables/genotype_table.tsv` — required genotype roll-up
- `tables/per_tumor.tsv` — unit-level table
- `tables/klk_vs_k.tsv` — primary contrast
- `tables/cldn4_positive_cells.tsv` — all 6 Cldn4>0 barcodes
- `tables/compartment_by_genotype.tsv`
- `tables/gene_inventory.tsv`, `tables/honest_n.tsv`, `tables/summary.json`
- `figures/fig1_tnk_fraction.png`, `fig2_cldn4_all.png`, `fig3_compartments.png`
- `scripts/analyze.py`

## 结论

GSE267321 是公开的 LKR13 皮下同基因瘤 scRNA（K / KK / KLK；**无 KL 库**）。**KLK（STK11/LKB1+KEAP1 缺失）是本系列最接近用户 KL 的一臂。** GEO 标题为 **non-malignant cells**：无作者恶性标签；标记残留上皮 **{n_epi_total}** 个细胞（皮下，不是正常肺）。Cldn4 行在，但全矩阵仅 **6** 个细胞 >0，**不能做 Cldn4 基因型检验**。T/NK 可以：KLK 两瘤都低于 K 两瘤（样本均数 0.095 vs 0.371）。IFN/MHC 在残留上皮上混合且 KLK-2 仅 11 个上皮细胞，不作结论。诚实 n = **2 vs 2 个瘤**。无 dual-high，无 ICI，不改 thesis。
"""


if __name__ == "__main__":
    main()
