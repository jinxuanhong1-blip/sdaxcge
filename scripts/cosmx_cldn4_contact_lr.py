#!/usr/bin/env python3
"""Official CosMx NSCLC: CLDN4-high tumor vs CD8 contacts and on-panel barrier ligands.

He et al. 2022 CosMx SMI 960-plex FFPE NSCLC (NanoString public S3).
Complementary to Ripley K: exclusion as missing contacts (binary 20 µm),
not as a continuous pair-correlation function.

Does NOT run CellChat and does not claim ligand-receptor causality.
Reports contact probability (odds ratio + FOV-restricted permutation)
and, only for genes present on the 960 panel, ligand enrichment among
the contacts that exist.

CLDN4-only. No TACSTD2 gate. No private 8-KL.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cosmx"
OUT = ROOT / "results" / "cosmx_cldn4_contact_lr"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = [
    "Lung5_Rep1",
    "Lung5_Rep2",
    "Lung5_Rep3",
    "Lung6",
    "Lung9_Rep1",
    "Lung9_Rep2",
    "Lung12",
    "Lung13",
]
# 5 tissues; Lung5 and Lung9 are technical replicates.
PATIENT = {
    "Lung5_Rep1": "Lung5",
    "Lung5_Rep2": "Lung5",
    "Lung5_Rep3": "Lung5",
    "Lung6": "Lung6",
    "Lung9_Rep1": "Lung9",
    "Lung9_Rep2": "Lung9",
    "Lung12": "Lung12",
    "Lung13": "Lung13",
}

PX_TO_UM = 0.18  # CosMx NSCLC prototype (He et al. 2022)
RADIUS_UM = 20.0
RADIUS_PX = RADIUS_UM / PX_TO_UM
MIN_COUNTS = 20
MIN_GENES = 5
MIN_FOV_TUMOR = 10
MIN_FOV_CD8 = 5
N_PERM = 5000
RNG = np.random.default_rng(20260819)

BARRIER_QUERY = ["F11R", "NECTIN2", "CDH1", "LGALS9"]
BARRIER_ALIASES = {
    "F11R": ["F11R", "JAM1", "JAM-A", "JAMA"],
    "NECTIN2": ["NECTIN2", "PVRL2", "HVEB", "CD112"],
    "CDH1": ["CDH1"],
    "LGALS9": ["LGALS9", "LGALS9A", "GAL9"],
}

EPI_MARKERS = [
    "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT17", "KRT15",
    "KRT14", "KRT13", "KRT16", "KRT6A", "KRT5", "KRT7", "CEACAM6",
    "MUC1", "ELF3", "SFTPC", "SFTPB", "NAPSA", "S100A14",
]
IMM_MARKERS = [
    "PTPRC", "CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "CD4",
    "IL7R", "CCL5", "GZMA", "GZMB", "GZMK", "GNLY", "KLRB1", "KLRK1",
    "CD68", "CD163", "C1QA", "C1QB", "C1QC", "LYZ", "FCGR3A", "ITGAX",
    "ITGAM", "CD14", "CD79A", "MS4A1", "CD19", "CD37", "CD52", "CD53",
]
STR_MARKERS = [
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "DCN",
    "LUM", "FN1", "BGN", "PDGFRB", "ACTA2", "MYH11", "PECAM1", "VWF",
    "CDH5", "CLEC14A", "FLT1", "KDR", "ESAM", "RAMP2",
]
CD8_GENES = ["CD8A", "CD8B"]
CD3_GENES = ["CD3D", "CD3E", "CD3G"]
CORE_GENES = ["CLDN4"] + CD8_GENES + CD3_GENES


def find_file(sample: str, kind: str) -> Path:
    hits = list(DATA.rglob(f"{sample}_{kind}_file.csv"))
    if not hits:
        raise FileNotFoundError(f"{sample} {kind} not under {DATA}")
    return hits[0]


def read_header(path: Path) -> list[str]:
    import csv

    with path.open() as f:
        return next(csv.reader(f))


def resolve_alias(panel: set[str], names: list[str]) -> str | None:
    for n in names:
        if n in panel:
            return n
    return None


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float, float]:
    if n <= 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / d
    half = (z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / d
    return p, max(0.0, center - half), min(1.0, center + half)


def haldane_or(a, b, c, d):
    """Odds ratio with Haldane-Anscombe 0.5 correction if any cell is 0."""
    a, b, c, d = float(a), float(b), float(c), float(d)
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    if b * c == 0:
        return np.nan, np.nan, np.nan
    or_ = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo, hi = np.exp(np.log(or_) - 1.96 * se), np.exp(np.log(or_) + 1.96 * se)
    return float(or_), float(lo), float(hi)


def mantel_haenszel(table: pd.DataFrame) -> dict:
    """table columns: a,b,c,d per stratum (high-contact, high-none, low-contact, low-none)."""
    a = table["a"].to_numpy(float)
    b = table["b"].to_numpy(float)
    c = table["c"].to_numpy(float)
    d = table["d"].to_numpy(float)
    n = a + b + c + d
    keep = n > 0
    a, b, c, d, n = a[keep], b[keep], c[keep], d[keep], n[keep]
    if len(n) == 0:
        return {"or_mh": np.nan, "or_mh_lo": np.nan, "or_mh_hi": np.nan, "n_strata": 0}
    r = a * d / n
    s = b * c / n
    if s.sum() == 0:
        return {"or_mh": np.nan, "or_mh_lo": np.nan, "or_mh_hi": np.nan, "n_strata": int(len(n))}
    or_mh = float(r.sum() / s.sum())
    # Robins-Breslow-Greenland variance of log(OR_MH)
    p = (a + d) / n
    q = (b + c) / n
    num = (p * r).sum() + (q * s).sum()
    den = 2 * r.sum() * s.sum()
    if den <= 0 or or_mh <= 0:
        lo = hi = np.nan
    else:
        se = np.sqrt(num / den)
        lo, hi = float(np.exp(np.log(or_mh) - 1.96 * se)), float(np.exp(np.log(or_mh) + 1.96 * se))
    return {"or_mh": or_mh, "or_mh_lo": lo, "or_mh_hi": hi, "n_strata": int(len(n))}


def lognorm_cols(raw: np.ndarray) -> np.ndarray:
    """Library-size log1p normalize a cells x genes count matrix."""
    tot = raw.sum(axis=1)
    med = np.median(tot[tot > 0]) if np.any(tot > 0) else 1.0
    sf = med / np.maximum(tot, 1.0)
    return np.log1p(raw * sf[:, None])


def score_mean(Xn: np.ndarray, genes: list[str], markers: list[str]) -> np.ndarray:
    idx = [genes.index(g) for g in markers if g in genes]
    if not idx:
        return np.zeros(Xn.shape[0])
    return Xn[:, idx].mean(axis=1)


def load_sample(sample: str, barrier_resolved: dict[str, str | None]) -> pd.DataFrame:
    expr_path = find_file(sample, "exprMat")
    meta_path = find_file(sample, "metadata")
    header = read_header(expr_path)
    panel = [c for c in header if c not in ("fov", "cell_ID") and not c.lower().startswith("neg")]
    panel_set = set(panel)

    want = set(CORE_GENES + EPI_MARKERS + IMM_MARKERS + STR_MARKERS)
    for canon, resolved in barrier_resolved.items():
        if resolved:
            want.add(resolved)
        want.update(BARRIER_ALIASES[canon])
    keep_genes = [g for g in panel if g in want]
    usecols = [c for c in ("fov", "cell_ID") if c in header] + keep_genes

    expr = pd.read_csv(expr_path, usecols=usecols)
    # Full-row QC from a cheap second pass on all gene columns would be heavy;
    # compute nCount/nGene from the on-disk matrix using a chunked sum.
    ncount, ngene = _qc_totals(expr_path, panel)
    expr = expr[expr["cell_ID"] != 0].copy()
    expr.index = expr["fov"].astype(str) + "_" + expr["cell_ID"].astype(str)
    expr["nCount"] = ncount.reindex(expr.index).to_numpy()
    expr["nGene"] = ngene.reindex(expr.index).to_numpy()

    meta = pd.read_csv(meta_path)
    meta.index = meta["fov"].astype(str) + "_" + meta["cell_ID"].astype(str)
    common = expr.index.intersection(meta.index)
    expr, meta = expr.loc[common], meta.loc[common]

    keep = (expr["nCount"] >= MIN_COUNTS) & (expr["nGene"] >= MIN_GENES)
    expr, meta = expr.loc[keep], meta.loc[keep]

    raw = expr[keep_genes].to_numpy(dtype=np.float32)
    Xn = lognorm_cols(raw)
    gmap = {g: i for i, g in enumerate(keep_genes)}

    def col(name, default=0.0):
        if name in gmap:
            return raw[:, gmap[name]]
        return np.full(len(expr), default, dtype=np.float32)

    def ncol(name, default=0.0):
        if name in gmap:
            return Xn[:, gmap[name]]
        return np.full(len(expr), default, dtype=np.float32)

    epi = score_mean(Xn, keep_genes, EPI_MARKERS)
    imm = score_mean(Xn, keep_genes, IMM_MARKERS)
    stro = score_mean(Xn, keep_genes, STR_MARKERS)
    S = np.vstack([epi, imm, stro])
    lab = np.array(["epithelial", "immune", "stromal"])[S.argmax(axis=0)]
    lab[S.max(axis=0) <= 0] = "unassigned"

    out = pd.DataFrame(
        {
            "sample": sample,
            "patient": PATIENT[sample],
            "fov": meta["fov"].to_numpy(),
            "cell_ID": meta["cell_ID"].to_numpy(),
            "x_px": meta["CenterX_global_px"].to_numpy(float),
            "y_px": meta["CenterY_global_px"].to_numpy(float),
            "nCount": expr["nCount"].to_numpy(float),
            "nGene": expr["nGene"].to_numpy(float),
            "compartment": lab,
            "cldn4_raw": col("CLDN4"),
            "cldn4": ncol("CLDN4"),
            "cd8_raw": col("CD8A") + col("CD8B"),
            "cd3_raw": col("CD3D") + col("CD3E") + col("CD3G"),
        },
        index=expr.index,
    )
    for canon, resolved in barrier_resolved.items():
        if resolved and resolved in gmap:
            out[f"{canon}_raw"] = raw[:, gmap[resolved]]
            out[canon] = Xn[:, gmap[resolved]]
        else:
            out[f"{canon}_raw"] = np.nan
            out[canon] = np.nan

    if "Mean.PanCK" in meta.columns:
        out["panck"] = np.log1p(meta["Mean.PanCK"].to_numpy(float))
        out["cd45"] = np.log1p(meta["Mean.CD45"].to_numpy(float)) if "Mean.CD45" in meta.columns else np.nan
        out["cd3_prot"] = np.log1p(meta["Mean.CD3"].to_numpy(float)) if "Mean.CD3" in meta.columns else np.nan
        out["prot_tumor"] = (out["panck"] >= np.median(out["panck"])) & (out["cd45"] < np.median(out["cd45"]))
    else:
        out["prot_tumor"] = False

    out["is_tumor"] = out["compartment"] == "epithelial"
    out["is_cd8"] = (
        (out["cd8_raw"] > 0)
        & (out["cd3_raw"] > 0)
        & (~out["is_tumor"])
    )
    return out


def _qc_totals(expr_path: Path, panel: list[str]) -> tuple[pd.Series, pd.Series]:
    """Chunked nCount / nGene over the full 960-plex (drop cell_ID==0 later)."""
    ncount = []
    ngene = []
    keys = []
    usecols = ["fov", "cell_ID"] + panel
    for chunk in pd.read_csv(expr_path, usecols=usecols, chunksize=20000):
        m = chunk["cell_ID"] != 0
        chunk = chunk.loc[m]
        mat = chunk[panel].to_numpy(dtype=np.float32)
        ncount.append(mat.sum(axis=1))
        ngene.append((mat > 0).sum(axis=1))
        keys.append(chunk["fov"].astype(str) + "_" + chunk["cell_ID"].astype(str))
    idx = pd.concat(keys, ignore_index=True)
    return (
        pd.Series(np.concatenate(ncount), index=idx),
        pd.Series(np.concatenate(ngene), index=idx),
    )


def assign_contacts(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["has_cd8_contact"] = False
    df["nn_cd8_um"] = np.nan
    # FOVs are unique only within a sample
    for (_, _fov), g in df.groupby(["sample", "fov"], sort=False):
        tumor = g["is_tumor"].to_numpy()
        cd8 = g["is_cd8"].to_numpy()
        xy = g[["x_px", "y_px"]].to_numpy(float)
        if tumor.sum() == 0 or cd8.sum() == 0:
            continue
        tree = cKDTree(xy[cd8])
        d, _ = tree.query(xy[tumor], k=1)
        tumor_idx = g.index[tumor]
        df.loc[tumor_idx, "has_cd8_contact"] = d <= RADIUS_PX
        df.loc[tumor_idx, "nn_cd8_um"] = d * PX_TO_UM
    return df


def fov_usable(df: pd.DataFrame) -> pd.Series:
    rows = []
    for (sample, fov), g in df.groupby(["sample", "fov"]):
        rows.append(
            {
                "sample": sample,
                "fov": fov,
                "n_tumor": int(g["is_tumor"].sum()),
                "n_cd8": int(g["is_cd8"].sum()),
            }
        )
    tab = pd.DataFrame(rows)
    tab["usable"] = (tab["n_tumor"] >= MIN_FOV_TUMOR) & (tab["n_cd8"] >= MIN_FOV_CD8)
    return tab


def tertile_split(values: np.ndarray) -> np.ndarray:
    """Return -1 / 0 / +1 for low / mid / high tertiles. Constant vectors → 0."""
    lab = np.zeros(len(values), dtype=int)
    if len(values) < 6 or np.unique(values).size < 3:
        return lab
    q1, q2 = np.quantile(values, [1 / 3, 2 / 3])
    lab[values <= q1] = -1
    lab[values >= q2] = 1
    return lab


def add_cldn4_classes(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["cldn4_tertile_fov"] = 0
    df["cldn4_tertile_sample"] = 0
    tumor = df["is_tumor"]
    for _, g in df.loc[tumor].groupby(["sample", "fov"]):
        df.loc[g.index, "cldn4_tertile_fov"] = tertile_split(g["cldn4"].to_numpy())
    for _, g in df.loc[tumor].groupby("sample"):
        df.loc[g.index, "cldn4_tertile_sample"] = tertile_split(g["cldn4"].to_numpy())
    return df


def counts_2x2(hi_contact, hi_n, lo_contact, lo_n):
    return {
        "a": int(hi_contact),
        "b": int(hi_n - hi_contact),
        "c": int(lo_contact),
        "d": int(lo_n - lo_contact),
    }


def or_from_group(g: pd.DataFrame, class_col: str) -> dict:
    hi = g[class_col] == 1
    lo = g[class_col] == -1
    a = int((hi & g["has_cd8_contact"]).sum())
    n_hi = int(hi.sum())
    c = int((lo & g["has_cd8_contact"]).sum())
    n_lo = int(lo.sum())
    cell = counts_2x2(a, n_hi, c, n_lo)
    or_, lo_ci, hi_ci = haldane_or(cell["a"], cell["b"], cell["c"], cell["d"])
    p_hi, p_hi_lo, p_hi_hi = wilson_ci(a, n_hi)
    p_lo, p_lo_lo, p_lo_hi = wilson_ci(c, n_lo)
    return {
        **cell,
        "n_hi": n_hi,
        "n_lo": n_lo,
        "p_contact_hi": p_hi,
        "p_contact_hi_lo": p_hi_lo,
        "p_contact_hi_hi": p_hi_hi,
        "p_contact_lo": p_lo,
        "p_contact_lo_lo": p_lo_lo,
        "p_contact_lo_hi": p_lo_hi,
        "or": or_,
        "or_lo": lo_ci,
        "or_hi": hi_ci,
    }


def permute_or(df: pd.DataFrame, class_col: str, n_perm: int, strat: list[str]) -> dict:
    """Shuffle high/low labels within strata; contact flags stay fixed."""
    work = df.loc[df[class_col] != 0, ["has_cd8_contact", class_col] + strat].copy()
    if len(work) < 20:
        return {"p_perm": np.nan, "n_perm": 0, "null_mean": np.nan, "obs_or": np.nan}
    obs = or_from_group(work, class_col)["or"]
    null = np.empty(n_perm)
    groups = list(work.groupby(strat, sort=False))
    labels = work[class_col].to_numpy().copy()
    contact = work["has_cd8_contact"].to_numpy()
    idx_slices = []
    start = 0
    # groupby iteration order matches row order only if we rebuild from groups
    rebuilt_lab = []
    rebuilt_con = []
    slices = []
    pos = 0
    for _, g in groups:
        lab = g[class_col].to_numpy().copy()
        con = g["has_cd8_contact"].to_numpy()
        slices.append((pos, pos + len(g), lab))
        rebuilt_lab.append(lab)
        rebuilt_con.append(con)
        pos += len(g)
    lab0 = np.concatenate(rebuilt_lab)
    con0 = np.concatenate(rebuilt_con)
    for i in range(n_perm):
        lab = lab0.copy()
        for lo, hi, src in slices:
            lab[lo:hi] = RNG.permutation(src)
        hi = lab == 1
        lo = lab == -1
        a = int((hi & con0).sum())
        b = int(hi.sum() - a)
        c = int((lo & con0).sum())
        d = int(lo.sum() - c)
        or_, _, _ = haldane_or(a, b, c, d)
        null[i] = or_
    # two-sided via |log OR|
    obs_l = np.log(obs) if obs > 0 else 0.0
    null_l = np.log(np.clip(null, 1e-12, None))
    p = (1 + np.sum(np.abs(null_l) >= abs(obs_l))) / (1 + n_perm)
    return {
        "obs_or": float(obs),
        "p_perm": float(p),
        "n_perm": n_perm,
        "null_mean": float(np.nanmean(null)),
        "null_q025": float(np.nanquantile(null, 0.025)),
        "null_q975": float(np.nanquantile(null, 0.975)),
        "frac_null_lt1": float(np.mean(null < 1)),
    }


def ligand_enrichment(df: pd.DataFrame, present: list[str], class_col: str) -> pd.DataFrame:
    rows = []
    contacts = df.loc[df["is_tumor"] & df["has_cd8_contact"] & (df[class_col] != 0)]
    for sample, g in contacts.groupby("sample"):
        hi = g[class_col] == 1
        lo = g[class_col] == -1
        if hi.sum() < 5 or lo.sum() < 5:
            continue
        for gene in present:
            xh = g.loc[hi, gene].to_numpy()
            xl = g.loc[lo, gene].to_numpy()
            rh = g.loc[hi, f"{gene}_raw"].to_numpy()
            rl = g.loc[lo, f"{gene}_raw"].to_numpy()
            try:
                u = stats.mannwhitneyu(xh, xl, alternative="two-sided")
                p_u = float(u.pvalue)
            except ValueError:
                p_u = np.nan
            a = int((rh > 0).sum())
            b = int((rh == 0).sum())
            c = int((rl > 0).sum())
            d = int((rl == 0).sum())
            or_, olo, ohi = haldane_or(a, b, c, d)
            rows.append(
                {
                    "sample": sample,
                    "patient": g["patient"].iloc[0],
                    "gene": gene,
                    "n_hi_contacts": int(hi.sum()),
                    "n_lo_contacts": int(lo.sum()),
                    "mean_lognorm_hi": float(np.mean(xh)),
                    "mean_lognorm_lo": float(np.mean(xl)),
                    "delta_lognorm": float(np.mean(xh) - np.mean(xl)),
                    "frac_pos_hi": float(np.mean(rh > 0)),
                    "frac_pos_lo": float(np.mean(rl > 0)),
                    "mwu_p": p_u,
                    "pos_or_hi_vs_lo": or_,
                    "pos_or_lo": olo,
                    "pos_or_hi": ohi,
                }
            )
    return pd.DataFrame(rows)


def permute_ligand_delta(df: pd.DataFrame, present: list[str], class_col: str, n_perm: int) -> pd.DataFrame:
    """FOV-restricted permutation of high/low among contacting tumor cells."""
    rows = []
    contacts = df.loc[df["is_tumor"] & df["has_cd8_contact"] & (df[class_col] != 0)].copy()
    if present and len(contacts) >= 20:
        family = contacts[present].mean(axis=1)
        contacts = contacts.assign(FAMILY_barrier=family)
        genes = present + ["FAMILY_barrier"]
    else:
        genes = list(present)
    for gene in genes:
        obs_parts = []
        slices = []
        for _, g in contacts.groupby(["sample", "fov"]):
            lab = g[class_col].to_numpy().copy()
            val = g[gene].to_numpy() if gene != "FAMILY_barrier" else g["FAMILY_barrier"].to_numpy()
            if ((lab == 1).sum() < 1) or ((lab == -1).sum() < 1):
                continue
            slices.append((lab, val))
            obs_parts.append(val[lab == 1].mean() - val[lab == -1].mean())
        if not slices:
            continue
        w = np.array([len(v) for _, v in slices], dtype=float)
        obs = float(np.average(obs_parts, weights=w))
        null = np.empty(n_perm)
        for i in range(n_perm):
            parts = []
            ww = []
            for lab, val in slices:
                plab = RNG.permutation(lab)
                parts.append(val[plab == 1].mean() - val[plab == -1].mean())
                ww.append(len(val))
            null[i] = np.average(parts, weights=ww)
        p = (1 + np.sum(np.abs(null) >= abs(obs))) / (1 + n_perm)
        rows.append(
            {
                "gene": gene,
                "obs_delta_lognorm": obs,
                "p_perm": float(p),
                "n_perm": n_perm,
                "null_mean": float(np.mean(null)),
                "null_q025": float(np.quantile(null, 0.025)),
                "null_q975": float(np.quantile(null, 0.975)),
            }
        )
    return pd.DataFrame(rows)


def plot_contact_bars(sample_or: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    samples = list(sample_or["sample"])
    x = np.arange(len(samples))
    w = 0.38
    hi = sample_or["p_contact_hi"].to_numpy()
    lo = sample_or["p_contact_lo"].to_numpy()
    hi_yerr = np.vstack(
        [
            np.clip(hi - sample_or["p_contact_hi_lo"].to_numpy(), 0, None),
            np.clip(sample_or["p_contact_hi_hi"].to_numpy() - hi, 0, None),
        ]
    )
    lo_yerr = np.vstack(
        [
            np.clip(lo - sample_or["p_contact_lo_lo"].to_numpy(), 0, None),
            np.clip(sample_or["p_contact_lo_hi"].to_numpy() - lo, 0, None),
        ]
    )
    ax.bar(x - w / 2, lo, w, yerr=lo_yerr, capsize=3, color="#4C78A8", label="CLDN4-low tumor (FOV Q1)", zorder=2)
    ax.bar(x + w / 2, hi, w, yerr=hi_yerr, capsize=3, color="#F58518", label="CLDN4-high tumor (FOV Q3)", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=30, ha="right")
    ax.set_ylabel("P(any CD8 contact < 20 µm)")
    ax.set_ylim(0, min(1.0, max(0.05, np.nanmax(np.concatenate([hi, lo])) * 1.35)))
    ax.axhline(0, color="k", lw=0.6)
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("Official CosMx NSCLC — tumor–CD8 contact probability by CLDN4 tertile")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_or_forest(sample_or: pd.DataFrame, pooled: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    y = np.arange(len(sample_or))
    ax.errorbar(
        sample_or["or"],
        y,
        xerr=[
            sample_or["or"] - sample_or["or_lo"],
            sample_or["or_hi"] - sample_or["or"],
        ],
        fmt="o",
        color="#4C78A8",
        capsize=3,
        label="Per sample",
    )
    ax.errorbar(
        [pooled["or_mh"]],
        [len(sample_or) + 0.8],
        xerr=[[pooled["or_mh"] - pooled["or_mh_lo"]], [pooled["or_mh_hi"] - pooled["or_mh"]]],
        fmt="D",
        color="#E45756",
        capsize=3,
        label="FOV-stratified MH",
    )
    ax.axvline(1.0, color="k", lw=0.8, ls="--")
    ax.set_yticks(list(y) + [len(sample_or) + 0.8])
    ax.set_yticklabels(list(sample_or["sample"]) + ["Pooled (MH)"])
    ax.set_xlabel("Odds ratio of any CD8 contact (CLDN4-high vs low)")
    ax.set_xscale("log")
    ax.legend(frameon=False)
    ax.set_title("OR < 1: CLDN4-high tumor less likely to have a CD8 contact")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_ligand_bars(lig: pd.DataFrame, path: Path) -> None:
    if lig.empty:
        return
    genes = list(lig["gene"].unique())
    samples = list(lig["sample"].unique())
    fig, axes = plt.subplots(1, len(genes), figsize=(4.0 * len(genes), 4.4), squeeze=False)
    for ax, gene in zip(axes[0], genes):
        sub = lig[lig["gene"] == gene]
        x = np.arange(len(sub))
        ax.bar(x - 0.18, sub["mean_lognorm_lo"], 0.36, color="#4C78A8", label="CLDN4-low contacts")
        ax.bar(x + 0.18, sub["mean_lognorm_hi"], 0.36, color="#F58518", label="CLDN4-high contacts")
        ax.set_xticks(x)
        ax.set_xticklabels(sub["sample"], rotation=30, ha="right")
        ax.set_title(gene)
        ax.set_ylabel("Mean log-norm in contacting tumor")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0][0].legend(frameon=False, fontsize=8)
    fig.suptitle("Barrier ligand expression among tumor–CD8 contacts (on-panel only)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def discover_panel() -> tuple[set[str], dict[str, str | None], dict]:
    """Read one available exprMat header; fall back to union if needed."""
    panel: set[str] = set()
    used = None
    for s in SAMPLES:
        try:
            header = read_header(find_file(s, "exprMat"))
        except FileNotFoundError:
            continue
        genes = [c for c in header if c not in ("fov", "cell_ID") and not c.lower().startswith("neg")]
        panel = set(genes)
        used = s
        break
    resolved = {canon: resolve_alias(panel, aliases) for canon, aliases in BARRIER_ALIASES.items()}
    present = [g for g, r in resolved.items() if r]
    missing = [g for g, r in resolved.items() if r is None]
    status = {
        "header_sample": used,
        "n_panel_genes": len(panel),
        "CLDN4_on_panel": "CLDN4" in panel,
        "CD8A_on_panel": "CD8A" in panel,
        "CD8B_on_panel": "CD8B" in panel,
        "barrier_resolved": resolved,
        "barrier_present": present,
        "barrier_missing": missing,
        "note": "Prototype CosMx NSCLC 960-plex (He 2022). Not the later 1K UCC panel.",
    }
    return panel, resolved, status


def main():
    available = []
    for s in SAMPLES:
        try:
            find_file(s, "exprMat")
            find_file(s, "metadata")
            available.append(s)
        except FileNotFoundError:
            print(f"SKIP missing {s}")
    if not available:
        raise SystemExit("No CosMx samples extracted under data/cosmx/")

    panel, resolved, panel_status = discover_panel()
    present = panel_status["barrier_present"]
    print("Panel barrier present:", present, "missing:", panel_status["barrier_missing"])
    if "CLDN4" not in panel:
        raise SystemExit("CLDN4 is not on the extracted CosMx panel")

    frames = []
    for s in available:
        print(f"== load {s}")
        frames.append(load_sample(s, resolved))
    df = pd.concat(frames, axis=0)
    print(f"QC cells: {len(df):,}")

    df = assign_contacts(df)
    df = add_cldn4_classes(df)
    fov_tab = fov_usable(df)
    usable_keys = set(zip(fov_tab.loc[fov_tab["usable"], "sample"], fov_tab.loc[fov_tab["usable"], "fov"]))
    df["fov_usable"] = [(s, f) in usable_keys for s, f in zip(df["sample"], df["fov"])]

    tumor = df.loc[df["is_tumor"] & df["fov_usable"]].copy()
    print(
        f"Usable tumor cells: {len(tumor):,} in {tumor.groupby(['sample','fov']).ngroups} FOVs; "
        f"CD8 cells in usable FOVs: {int(df.loc[df['fov_usable'] & df['is_cd8']].shape[0]):,}"
    )

    # Per-sample OR (FOV tertiles)
    sample_rows = []
    for sample, g in tumor.groupby("sample"):
        row = or_from_group(g, "cldn4_tertile_fov")
        row.update({"sample": sample, "patient": PATIENT[sample], "n_tumor": int(len(g))})
        sample_rows.append(row)
    sample_or = pd.DataFrame(sample_rows)

    # Per-FOV 2x2 for MH
    fov_rows = []
    for (sample, fov), g in tumor.groupby(["sample", "fov"]):
        row = or_from_group(g, "cldn4_tertile_fov")
        row.update({"sample": sample, "patient": PATIENT[sample], "fov": int(fov)})
        fov_rows.append(row)
    fov_or = pd.DataFrame(fov_rows)
    mh = mantel_haenszel(fov_or)

    # Simple pooled 2x2
    pooled_cells = or_from_group(tumor, "cldn4_tertile_fov")

    # Permutation (FOV-restricted)
    print(f"Permutation n={N_PERM} (FOV-restricted CLDN4 labels)")
    perm = permute_or(tumor, "cldn4_tertile_fov", N_PERM, strat=["sample", "fov"])

    # Sensitivity: sample tertiles, median split, protein-tumor index
    sens = {}
    sens["sample_tertiles"] = or_from_group(tumor, "cldn4_tertile_sample")
    med = tumor.copy()
    med["cldn4_med"] = np.where(
        med["cldn4"] >= med.groupby("sample")["cldn4"].transform("median"), 1, -1
    )
    sens["median_split"] = or_from_group(med, "cldn4_med")
    prot = df.loc[df["prot_tumor"] & df["fov_usable"]].copy()
    if len(prot) >= 100:
        prot = add_cldn4_classes(prot.assign(is_tumor=True))
        # contacts already assigned on full df; recompute classes on protein-tumor
        prot["cldn4_tertile_fov"] = 0
        for _, g in prot.groupby(["sample", "fov"]):
            prot.loc[g.index, "cldn4_tertile_fov"] = tertile_split(g["cldn4"].to_numpy())
        sens["protein_panck_tumor"] = or_from_group(prot, "cldn4_tertile_fov")
        sens["protein_panck_tumor"]["n_index"] = int(len(prot))

    # Ligands among contacts
    lig = ligand_enrichment(tumor, present, "cldn4_tertile_fov")
    lig_perm = permute_ligand_delta(tumor, present, "cldn4_tertile_fov", N_PERM)
    if not lig.empty:
        lig["mwu_q"] = stats.false_discovery_control(lig["mwu_p"].fillna(1.0), method="bh")

    # Write tables
    sample_or.to_csv(OUT / "contact_probability_by_sample.tsv", sep="\t", index=False)
    fov_or.to_csv(OUT / "contact_2x2_by_fov.tsv", sep="\t", index=False)
    fov_tab.to_csv(OUT / "fov_inventory.tsv", sep="\t", index=False)
    if not lig.empty:
        lig.to_csv(OUT / "ligand_enrichment_contacts.tsv", sep="\t", index=False)
    if not lig_perm.empty:
        lig_perm.to_csv(OUT / "ligand_enrichment_permutation.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "gene": canon,
                "on_960_panel": resolved[canon] is not None,
                "panel_symbol": resolved[canon],
                "aliases_tried": ",".join(BARRIER_ALIASES[canon]),
            }
            for canon in BARRIER_QUERY
        ]
    ).to_csv(OUT / "panel_gene_status.tsv", sep="\t", index=False)

    # Compact cell-level export (tumor only, usable FOVs)
    keep_cols = [
        "sample", "patient", "fov", "cldn4", "cldn4_raw", "cldn4_tertile_fov",
        "has_cd8_contact", "nn_cd8_um", "compartment",
    ] + present + [f"{g}_raw" for g in present]
    tumor[keep_cols].to_csv(OUT / "tumor_cells_usable.tsv.gz", sep="\t", index=False)

    plot_contact_bars(sample_or, OUT / "contact_probability_bars.png")
    plot_or_forest(sample_or, {**mh, **{k: mh[k] for k in mh}}, OUT / "contact_or_forest.png")
    if not lig.empty:
        plot_ligand_bars(lig, OUT / "ligand_enrichment_bars.png")

    inventory = {
        "dataset": "Official CosMx SMI NSCLC FFPE (He et al. 2022 Nat Biotechnol)",
        "access": "nanostring-public-share S3 SMI-Compressed flat files",
        "n_samples_available": len(available),
        "samples": available,
        "n_patients": len({PATIENT[s] for s in available}),
        "px_to_um": PX_TO_UM,
        "contact_radius_um": RADIUS_UM,
        "cldn4_split": "per-FOV tertiles among RNA-epithelial tumor cells (Q3 vs Q1)",
        "tumor_definition": "RNA epithelial marker-score argmax (EPCAM/KRT/CDH1 family)",
        "cd8_definition": "CD8A|CD8B > 0 AND CD3D|CD3E|CD3G > 0 AND not epithelial",
        "n_qc_cells": int(len(df)),
        "n_tumor_usable": int(len(tumor)),
        "n_cd8_usable_fovs": int((df["fov_usable"] & df["is_cd8"]).sum()),
        "n_tumor_with_contact": int(tumor["has_cd8_contact"].sum()),
        "contact_rate_all_tumor": float(tumor["has_cd8_contact"].mean()),
        "n_usable_fovs": int(fov_tab["usable"].sum()),
        "n_fovs_total": int(len(fov_tab)),
    }
    stats_out = {
        "inventory": inventory,
        "panel": panel_status,
        "pooled_2x2": pooled_cells,
        "mantel_haenszel_fov": mh,
        "permutation_fov": perm,
        "sensitivity": sens,
        "per_sample": sample_or.to_dict(orient="records"),
        "ligand_permutation": lig_perm.to_dict(orient="records") if not lig_perm.empty else [],
        "caveats": [
            "Contact is a centroid-distance proxy (<20 µm), not a membrane-contact call.",
            "Flat release has no official 18-type labels; tumor/CD8 are RNA-marker definitions.",
            "Do not interpret ligand enrichment as CellChat communication probability or causality.",
            "No ICI / response labels. Official CosMx NSCLC demo tissues only.",
            "No private 8-KL data.",
        ],
    }
    with open(OUT / "stats.json", "w") as f:
        json.dump(stats_out, f, indent=2, default=float)

    print(json.dumps(
        {
            "samples": available,
            "contact_rate": inventory["contact_rate_all_tumor"],
            "pooled_or": pooled_cells["or"],
            "mh_or": mh["or_mh"],
            "p_perm": perm.get("p_perm"),
            "barrier_missing": panel_status["barrier_missing"],
        },
        indent=2,
    ))


if __name__ == "__main__":
    main()
