#!/usr/bin/env python3
"""GSE137244 KL vs KP: transferable signatures and a method sweep.

Re-downloads the public GEO supplementary matrices. Builds frozen mouse
gene lists (TJ, Cldn4 neighborhood, NHEJ, STING, IFN, APM) and a gene-level
KL-versus-KP table for later scoring on other mouse scRNA.

The pattern sweep asks which pre-specified scale, summary, cohort, NHEJ
list, and IFN list most strongly shows KL > KP for Cldn4 and Tacstd2
together with lower NHEJ and higher IFN. Gene membership is not trimmed
to the desired sign.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import urllib.request
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import gene_sets as gs

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGS = HERE / "figures"
DATA = Path(os.environ.get("GSE137244_DIR", "/tmp/gse137244"))

URLS = {
    "raw": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.raw.csv.gz",
    "normalized": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.normalized.csv.gz",
    "fpkm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.fpkm.csv.gz",
}

KP_LIBS = [
    "B6AL10-1-RNA",
    "B6AL10-2-RNA",
    "B6AL10-3-RNA",
    "B6AL10-4-RNA",
    "B6AL10-5-RNA",
]
KL_LIBS = [
    "KL155mix-control-2-RNA",
    "KL47-1-untreated-1-RNA",
    "KLC-RNA",
    "KLD-RNA",
    "KLE-RNA",
]
NORMAL = "normal-lung-RNA"
EPCAM_DROP = "B6AL10-3-RNA"

# Within-arm Cldn4 neighborhood. Highest threshold that keeps at least this
# many genes other than Cldn4. Fixed before looking at which genes pass.
NEIGH_THRESHOLDS = (0.9, 0.8, 0.7, 0.6, 0.5)
NEIGH_MIN_GENES = 8
NEIGH_MIN_MEAN_FPKM = 1.0

NHEJ_SWEEP = ("NHEJ_CORE", "NHEJ_KEGG", "NHEJ_EXTENDED")
IFN_SWEEP = ("IFN_COMPACT", "IFN_HALLMARK_ALPHA", "IFN_HALLMARK_GAMMA", "IFN_OAS")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download() -> dict[str, str]:
    DATA.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for key, url in URLS.items():
        dest = DATA / Path(url).name
        if not dest.exists() or dest.stat().st_size < 1000:
            urllib.request.urlretrieve(url, dest)
        hashes[key] = sha256(dest)
    return hashes


def _strip_token(text: str) -> str:
    return text.strip().strip("'").strip('"').strip()


def load_matrix(path: Path, quoted_rows: bool) -> pd.DataFrame:
    if quoted_rows:
        rows = []
        with gzip.open(path, "rt", newline="") as handle:
            for line in handle:
                line = line.strip()
                if line.startswith('"') and line.endswith('"'):
                    line = line[1:-1]
                rows.append(line.split(","))
        header = [_strip_token(c) for c in rows[0]]
        df = pd.DataFrame(
            [[float(x) for x in rec[1:]] for rec in rows[1:]],
            index=[_strip_token(rec[0]) for rec in rows[1:]],
            columns=header[1:],
        )
    else:
        df = pd.read_csv(path, index_col=0)
        df.columns = [_strip_token(c) for c in df.columns]
        df.index = [_strip_token(str(i)) for i in df.index]
        df = df.apply(pd.to_numeric, errors="coerce")
    return df


def align_symbol_index(raw: pd.DataFrame, other: pd.DataFrame, label: str) -> pd.DataFrame:
    """FPKM and normalized files store March/Sept symbols as Excel dates.

    Row order matches the raw-count file, which still has March1 / Sept1.
    """
    if len(raw.index) != len(other.index):
        raise SystemExit(f"{label} row count {len(other.index)} != raw {len(raw.index)}")
    mismatch = [
        (a, b)
        for a, b in zip(raw.index, other.index.astype(str))
        if a != b
    ]
    bad = [(a, b) for a, b in mismatch if not (a.startswith(("Marc", "March", "Sept")) and b[0].isdigit())]
    if bad:
        raise SystemExit(f"{label} symbol mismatches are not only Excel dates: {bad[:8]}")
    out = other.copy()
    out.index = raw.index
    if out.index.duplicated().any():
        raise SystemExit(f"{label} still has duplicate symbols after the date repair")
    return out


def log2p1(df: pd.DataFrame) -> pd.DataFrame:
    return np.log2(df.astype(float) + 1.0)


def median_of_ratios(raw: pd.DataFrame) -> pd.DataFrame:
    """DESeq2-style size factors. This is the normalization, not the NB Wald test."""
    positive = raw.where(raw > 0)
    log_geo = np.log(positive).mean(axis=1)
    geo = np.exp(log_geo)
    ratios = raw.div(geo, axis=0)
    size = ratios.median(axis=0, skipna=True)
    if (size <= 0).any() or size.isna().any():
        raise SystemExit("median-of-ratios size factors failed")
    return raw.div(size, axis=1)


def cpm(raw: pd.DataFrame) -> pd.DataFrame:
    depths = raw.sum(axis=0)
    return raw.div(depths, axis=1) * 1e6


def cliff_and_mw(kl: np.ndarray, kp: np.ndarray) -> dict:
    kl = np.asarray(kl, dtype=float)
    kp = np.asarray(kp, dtype=float)
    n_gt = int(np.sum(kl[:, None] > kp[None, :]))
    n_lt = int(np.sum(kl[:, None] < kp[None, :]))
    n_eq = int(kl.size * kp.size - n_gt - n_lt)
    cliff = (n_gt - n_lt) / (kl.size * kp.size)
    mw = stats.mannwhitneyu(kl, kp, alternative="two-sided", method="exact")
    return {
        "n_kl_gt_kp": n_gt,
        "n_kl_lt_kp": n_lt,
        "n_tie": n_eq,
        "cliff": float(cliff),
        "mw_u": float(mw.statistic),
        "mw_p": float(mw.pvalue),
        "separation": "KL>KP" if n_lt == 0 and n_eq == 0 else ("KL<KP" if n_gt == 0 and n_eq == 0 else "overlap"),
    }


def welch(kl: np.ndarray, kp: np.ndarray) -> tuple[float, float]:
    res = stats.ttest_ind(kl, kp, equal_var=False, alternative="two-sided")
    return float(res.statistic), float(res.pvalue)


def contrast_vector(values: pd.Series, kl_libs: list[str], kp_libs: list[str]) -> dict:
    kl = values.loc[kl_libs].to_numpy(dtype=float)
    kp = values.loc[kp_libs].to_numpy(dtype=float)
    out = {
        "mean_kl": float(np.mean(kl)),
        "mean_kp": float(np.mean(kp)),
        "delta_mean": float(np.mean(kl) - np.mean(kp)),
        "median_kl": float(np.median(kl)),
        "median_kp": float(np.median(kp)),
        "delta_median": float(np.median(kl) - np.median(kp)),
    }
    out.update(cliff_and_mw(kl, kp))
    t, p = welch(kl, kp)
    out["welch_t"] = t
    out["welch_p"] = p
    return out


def mean_score(log_expr: pd.DataFrame, genes: list[str], cols: list[str]) -> pd.Series:
    return log_expr.loc[genes, cols].mean(axis=0)


def median_score(log_expr: pd.DataFrame, genes: list[str], cols: list[str]) -> pd.Series:
    return log_expr.loc[genes, cols].median(axis=0)


def spearman_rows(mat: pd.DataFrame, vector: np.ndarray) -> pd.Series:
    """Spearman of every row against one vector. Average ranks, NaN if constant."""
    ranks = np.apply_along_axis(stats.rankdata, 1, mat.to_numpy(dtype=float))
    v = stats.rankdata(np.asarray(vector, dtype=float))
    a = ranks - ranks.mean(axis=1, keepdims=True)
    v0 = v - v.mean()
    num = a @ v0
    den = np.sqrt((a ** 2).sum(axis=1) * np.sum(v0 ** 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = num / den
    rho[~np.isfinite(rho)] = np.nan
    return pd.Series(rho, index=mat.index)


def bh(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    out = pd.Series(np.nan, index=p.index)
    ok = p.notna()
    if ok.sum() == 0:
        return out
    vals = p[ok].to_numpy()
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(ranked)
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty(n)
    restored[order] = q
    out.loc[ok] = restored
    return out


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    for num, den, label in (
        (2, 252, "0.00794"),
        (4, 252, "0.0159"),
        (8, 252, "0.0317"),
        (14, 252, "0.0556"),
    ):
        if abs(p - num / den) < 1e-6:
            return label
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4g}"


def fmt_d(x: float) -> str:
    return f"{x:+.3f}"


def build_sets(universe: list[str]) -> tuple[dict[str, list[str]], pd.DataFrame]:
    curated = {
        "TJ_EPITHELIAL": (gs.TJ_EPITHELIAL, "18 structural epithelial tight-junction genes"),
        "TJ_TISMO": (gs.TJ_TISMO, "TISMO TJ list: Cldn3/4/6/7, Cdh1, F11r, Ocln"),
        "CLDN4_TJ_EDGE": (gs.CLDN4_TJ_EDGE, "prior TJ-core leading edge, vendored"),
        "NHEJ_CORE": (gs.NHEJ_CORE, "c-NHEJ enzymes plus Trp53bp1"),
        "NHEJ_EXTENDED": (gs.NHEJ_EXTENDED, "53BP1-pathway and end-sensing genes"),
        "NHEJ_KEGG": (gs.read_symbol_file("kegg2019_mouse_nhej.symbols.txt"), "KEGG 2019 Mouse non-homologous end-joining"),
        "STING_CORE": (gs.STING_CORE, "Mb21d1, Tmem173, Tbk1, Ikbke, Irf3"),
        "IFN_COMPACT": (gs.IFN_COMPACT, "15-gene epithelial ISG list"),
        "IFN_OAS": (gs.IFN_OAS, "OAS antiviral enzymes"),
        "IFN_HALLMARK_ALPHA": (
            gs.read_symbol_file("hallmark_ifn_alpha_2020.symbols.txt"),
            "MSigDB Hallmark 2020 interferon alpha, case-matched to this matrix",
        ),
        "IFN_HALLMARK_GAMMA": (
            gs.read_symbol_file("hallmark_ifn_gamma_2020.symbols.txt"),
            "MSigDB Hallmark 2020 interferon gamma, case-matched to this matrix",
        ),
        "APM_MHCI": (gs.APM_MHCI, "MHC-I antigen presentation and immunoproteasome"),
        "APM_KEGG": (
            gs.read_symbol_file("kegg2019_mouse_antigen_processing.symbols.txt"),
            "KEGG 2019 Mouse antigen processing and presentation",
        ),
    }
    matched: dict[str, list[str]] = {}
    rows = []
    for name, (genes, source) in curated.items():
        used, missing = gs.match_symbols(genes, universe)
        matched[name] = used
        for gene in used:
            rows.append({"set_name": name, "gene": gene, "status": "used", "source": source})
        for gene in missing:
            rows.append({"set_name": name, "gene": gene, "status": "absent", "source": source})
    return matched, pd.DataFrame(rows)


def neighborhood(
    log_fpkm: pd.DataFrame,
    fpkm: pd.DataFrame,
    kl_libs: list[str],
    kp_libs: list[str],
) -> tuple[pd.DataFrame, dict]:
    cols = kp_libs + kl_libs
    mean_fpkm = fpkm[cols].mean(axis=1)
    keep = mean_fpkm >= NEIGH_MIN_MEAN_FPKM
    mat = log_fpkm.loc[keep, cols]
    anchor = log_fpkm.loc["Cldn4", cols].to_numpy(dtype=float)
    rho_all = spearman_rows(mat, anchor)
    rho_kl = spearman_rows(mat[kl_libs], log_fpkm.loc["Cldn4", kl_libs].to_numpy(dtype=float))
    rho_kp = spearman_rows(mat[kp_libs], log_fpkm.loc["Cldn4", kp_libs].to_numpy(dtype=float))
    table = pd.DataFrame(
        {
            "gene": mat.index,
            "mean_fpkm": mean_fpkm.loc[mat.index].to_numpy(),
            "spearman_all10": rho_all.to_numpy(),
            "spearman_within_kl": rho_kl.to_numpy(),
            "spearman_within_kp": rho_kp.to_numpy(),
        }
    ).set_index("gene")
    table["spearman_within_mean"] = table[["spearman_within_kl", "spearman_within_kp"]].mean(axis=1)
    both_pos = (table["spearman_within_kl"] > 0) & (table["spearman_within_kp"] > 0)
    table["within_candidate"] = both_pos & table["spearman_within_mean"].notna() & (table.index != "Cldn4")

    thresh_rows = []
    chosen = None
    for thr in NEIGH_THRESHOLDS:
        genes = table.index[table["within_candidate"] & (table["spearman_within_mean"] >= thr)].tolist()
        score = mean_score(log_fpkm, genes, cols) if genes else pd.Series(np.nan, index=cols)
        rec = {"threshold": thr, "n_genes": len(genes), "rule": "within-arm Spearman, both arms > 0, mean FPKM >= 1, Cldn4 excluded"}
        if genes:
            rec.update(contrast_vector(score, kl_libs, kp_libs))
        else:
            rec.update({"delta_mean": np.nan, "mw_p": np.nan, "cliff": np.nan, "separation": "empty"})
        thresh_rows.append(rec)
        if chosen is None and len(genes) >= NEIGH_MIN_GENES:
            chosen = {"threshold": thr, "genes": genes}
    if chosen is None:
        # Fall through to the loosest threshold and keep whatever passed.
        thr = NEIGH_THRESHOLDS[-1]
        genes = table.index[table["within_candidate"] & (table["spearman_within_mean"] >= thr)].tolist()
        chosen = {"threshold": thr, "genes": genes, "below_min": True}
    else:
        chosen["below_min"] = False
    return table, {"thresholds": thresh_rows, "chosen": chosen}


def pattern_row(
    scale: str,
    summary: str,
    cohort: str,
    nhej_name: str,
    ifn_name: str,
    cldn4: dict,
    tacstd2: dict,
    nhej: dict,
    ifn: dict,
    n_nhej: int,
    n_ifn: int,
) -> dict:
    d_key = "delta_mean" if summary == "mean" else "delta_median"
    limbs = {
        "cldn4_kl_gt_kp": cldn4[d_key] > 0,
        "tacstd2_kl_gt_kp": tacstd2[d_key] > 0,
        "nhej_kl_lt_kp": nhej[d_key] < 0,
        "ifn_kl_gt_kp": ifn[d_key] > 0,
    }
    # Cliff is comparable across scales. Desired direction is positive.
    strength = (
        cldn4["cliff"]
        + tacstd2["cliff"]
        + (-nhej["cliff"])
        + ifn["cliff"]
    )
    return {
        "scale": scale,
        "summary": summary,
        "cohort": cohort,
        "nhej_set": nhej_name,
        "ifn_set": ifn_name,
        "n_nhej": n_nhej,
        "n_ifn": n_ifn,
        "cldn4_delta": cldn4[d_key],
        "tacstd2_delta": tacstd2[d_key],
        "nhej_delta": nhej[d_key],
        "ifn_delta": ifn[d_key],
        "cldn4_cliff": cldn4["cliff"],
        "tacstd2_cliff": tacstd2["cliff"],
        "nhej_cliff": nhej["cliff"],
        "ifn_cliff": ifn["cliff"],
        "cldn4_mw_p": cldn4["mw_p"],
        "tacstd2_mw_p": tacstd2["mw_p"],
        "nhej_mw_p": nhej["mw_p"],
        "ifn_mw_p": ifn["mw_p"],
        "n_limbs": int(sum(limbs.values())),
        "hit_nhej_down_ifn_up": bool(limbs["nhej_kl_lt_kp"] and limbs["ifn_kl_gt_kp"] and limbs["cldn4_kl_gt_kp"] and limbs["tacstd2_kl_gt_kp"]),
        "strength_cliff": float(strength),
        **{f"limb_{k}": bool(v) for k, v in limbs.items()},
    }


def write_gmt(matched: dict[str, list[str]], neighborhood_genes: list[str], neigh_note: str) -> None:
    descriptions = {
        "TJ_EPITHELIAL": "Epithelial tight junction. Score = mean of log-normalized expression.",
        "TJ_TISMO": "Short TJ list (Cldn3/4/6/7, Cdh1, F11r, Ocln).",
        "CLDN4_TJ_EDGE": "Vendored TJ leading edge around Cldn4.",
        "CLDN4_NEIGHBORHOOD": neigh_note + " On GSE137244 the mean score does not separate KL from KP.",
        "NHEJ_CORE": "c-NHEJ enzymes plus Trp53bp1.",
        "NHEJ_KEGG": "KEGG 2019 Mouse NHEJ symbols matched to this matrix.",
        "NHEJ_EXTENDED": "53BP1-pathway and MRN-side genes.",
        "STING_CORE": "cGAS-STING machinery, not ISGs. Mb21d1 is cGAS; Tmem173 is STING.",
        "IFN_COMPACT": "15-gene epithelial ISG list.",
        "IFN_OAS": "OAS antiviral enzymes.",
        "IFN_HALLMARK_ALPHA": "Hallmark interferon alpha symbols present in this mouse matrix.",
        "IFN_HALLMARK_GAMMA": "Hallmark interferon gamma symbols present in this mouse matrix.",
        "APM_MHCI": "MHC-I antigen presentation and immunoproteasome. Tumor-cell score.",
        "APM_KEGG": "KEGG antigen processing. Includes immune-cell genes.",
    }
    payload = dict(matched)
    payload["CLDN4_NEIGHBORHOOD"] = neighborhood_genes
    lines = []
    for name, genes in payload.items():
        if not genes:
            continue
        desc = descriptions[name]
        lines.append("\t".join([name, desc, *genes]))
    (TABLES / "signatures.gmt").write_text("\n".join(lines) + "\n")


def plot_modules(scores: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    order = [
        "TJ_EPITHELIAL",
        "CLDN4_NEIGHBORHOOD",
        "NHEJ_CORE",
        "STING_CORE",
        "IFN_COMPACT",
        "APM_MHCI",
    ]
    labels = {
        "TJ_EPITHELIAL": "TJ epithelial",
        "CLDN4_NEIGHBORHOOD": "Cldn4 neighborhood",
        "NHEJ_CORE": "NHEJ core",
        "STING_CORE": "STING",
        "IFN_COMPACT": "IFN compact",
        "APM_MHCI": "APM MHC-I",
    }
    fig, axes = plt.subplots(2, 3, figsize=(9.2, 6.2), sharex=False)
    rng = np.random.default_rng(1)
    for ax, name in zip(axes.ravel(), order):
        sub = scores[scores["set_name"] == name]
        for genotype, color, x in (("KP", "#3C6E8F", 0), ("KL", "#B85C38", 1)):
            vals = sub.loc[sub["genotype"] == genotype, "score_mean_log2fpkm1"].to_numpy()
            jitter = rng.uniform(-0.08, 0.08, size=len(vals))
            ax.scatter(np.full(len(vals), x) + jitter, vals, s=28, color=color, zorder=3)
            ax.hlines(np.mean(vals), x - 0.18, x + 0.18, color=color, lw=1.6)
        row = contrasts[contrasts["set_name"] == name].iloc[0]
        ax.set_xticks([0, 1], ["KP", "KL"])
        ax.set_title(labels[name], fontsize=11)
        ax.text(
            0.02,
            0.98,
            f"Δ {row['delta_mean']:+.2f}\nMW {fmt_p(row['mw_p'])}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("GSE137244 library scores, log2(FPKM+1) mean", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_modules_kl_kp.png", dpi=160)
    fig.savefig(FIGS / "fig_modules_kl_kp.pdf")
    plt.close(fig)


def plot_sweep(sweep: pd.DataFrame, gene_delta: pd.Series) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.8))
    ax = axes[0]
    # Full grid, so median summaries that flip the NHEJ sign are visible.
    colors = {"fpkm": "#3C6E8F", "normalized": "#C47B2B", "cpm": "#6B7C3A", "median_ratio": "#7B4B94"}
    for scale, sub in sweep.groupby("scale"):
        ax.scatter(
            sub["nhej_delta"],
            sub["ifn_delta"],
            s=18,
            alpha=0.55,
            color=colors.get(scale, "#333333"),
            label=scale,
        )
    ax.axvline(0, color="#888888", lw=0.6)
    ax.axhline(0, color="#888888", lw=0.6)
    n_hit = int(sweep["hit_nhej_down_ifn_up"].sum())
    ax.text(
        0.03,
        0.97,
        f"{n_hit} / {len(sweep)} with NHEJ down and IFN up",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )
    ax.set_xlabel("NHEJ Δ (KL − KP)")
    ax.set_ylabel("IFN Δ (KL − KP)")
    ax.set_title("Method sweep")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    ordered = gene_delta.sort_values()
    bar_colors = ["#B85C38" if v > 0 else "#3C6E8F" for v in ordered]
    ax.barh(ordered.index, ordered.to_numpy(), color=bar_colors)
    ax.axvline(0, color="#888888", lw=0.6)
    ax.set_xlabel("Δ log2(FPKM+1), KL − KP")
    ax.set_title("TJ genes and Tacstd2")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_sweep_neighborhood.png", dpi=160)
    fig.savefig(FIGS / "fig_sweep_neighborhood.pdf")
    plt.close(fig)


def write_finding(ctx: dict) -> None:
    locked = ctx["locked"]
    mods = ctx["primary_modules"]
    sweep = ctx["sweep_summary"]
    neigh = ctx["neighborhood"]
    lines = []
    a = lines.append
    a("# FINDING — GSE137244 transferable KL vs KP signatures")
    a("")
    a("Public cell-line RNA-seq only. Deng et al., *Nature Cancer* 2021 ([PMID 34142094](https://pubmed.ncbi.nlm.nih.gov/34142094/)), [GSE137244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137244). Matrices re-downloaded from the GEO supplementary files. Normal lung is held out. No private 8-KL matrices. No TISMO LLC.")
    a("")
    a("Primary scale is the mean of log2(FPKM+1), KL minus KP. That is the scale of the locked Tacstd2 and Cldn4 deltas. The rank test is the two-sided exact Mann–Whitney test. For 5 vs 5 with complete separation the p-value is 2/252 = 0.00794, and it cannot go lower.")
    a("")
    a("The deposited FPKM and normalized files store March- and Sept-family symbols as Excel dates (`1-Mar`, `1-Sep`). Row order matches the raw-count file, which still has `March1` and `Sept1`. Those 26 symbols were restored from the raw-count row order before any test. None of them is in the pre-specified TJ, NHEJ, STING, IFN, or APM lists.")
    a("")
    a("## Samples")
    a("")
    a("KP (KrasG12D; Trp53, GEO description) is B6AL10-1 through B6AL10-5. KL (KrasG12D; Lkb1) is KL155mix-control-2, KL47-1-untreated-1, KLC, KLD, and KLE. Titles that say control or untreated are the baseline libraries in this series. Honest n for the primary contrast is 5 libraries vs 5 libraries. Several libraries inside an arm are near-replicates; the maximum within-arm Pearson r of log2(FPKM+1) on genes with mean FPKM ≥ 1 is **{:.3f}**.".format(ctx["max_within_arm_r"]))
    a("")
    a("Genotype check on the same scale: Stk11 Δ = **{}** (exact p {}). Trp53 Δ = **{}** (exact p {}). Every KL library is below every KP library for Stk11, and every KL library is above every KP library for Trp53.".format(
        fmt_d(locked["Stk11"]["delta_mean"]),
        fmt_p(locked["Stk11"]["mw_p"]),
        fmt_d(locked["Trp53"]["delta_mean"]),
        fmt_p(locked["Trp53"]["mw_p"]),
    ))
    a("")
    a("## Locked genes")
    a("")
    a("| gene | Δ log2(FPKM+1) | exact MW p | separation |")
    a("|---|---:|---:|---|")
    for gene in ("Tacstd2", "Cldn4"):
        rec = locked[gene]
        a(f"| {gene} | **{fmt_d(rec['delta_mean'])}** | {fmt_p(rec['mw_p'])} | {rec['separation']} |")
    a("")
    a("Handoff rounding is Tacstd2 +3.24 and Cldn4 +5.57. Both p-values are the 5-vs-5 floor.")
    a("")
    a("## Signatures to score later")
    a("")
    a("Lists are in `tables/signatures.gmt` (mouse symbols as they appear in this mm10 matrix). The gene-level KL-versus-KP table is `tables/de_kl_vs_kp.tsv.gz`. On another mouse scRNA matrix, score each set as the mean of log-normalized expression of the genes that set and the matrix share. Do not re-cut the gene list on the query.")
    a("")
    a("| signature | genes used | Δ mean log2(FPKM+1) | exact MW p | separation |")
    a("|---|---:|---:|---:|---|")
    for name in (
        "TJ_EPITHELIAL",
        "TJ_TISMO",
        "CLDN4_TJ_EDGE",
        "CLDN4_NEIGHBORHOOD",
        "NHEJ_CORE",
        "NHEJ_KEGG",
        "NHEJ_EXTENDED",
        "STING_CORE",
        "IFN_COMPACT",
        "IFN_HALLMARK_ALPHA",
        "IFN_HALLMARK_GAMMA",
        "IFN_OAS",
        "APM_MHCI",
        "APM_KEGG",
    ):
        rec = mods[name]
        a(f"| {name} | {rec['n_genes']} | {fmt_d(rec['delta_mean'])} | {fmt_p(rec['mw_p'])} | {rec['separation']} |")
    a("")
    a("Pre-specified TJ means on this scale are {} (18 epithelial genes) and {} (7-gene TISMO list). The handoff’s TJ figure of +3.03 is a different average.".format(
        fmt_d(mods["TJ_EPITHELIAL"]["delta_mean"]),
        fmt_d(mods["TJ_TISMO"]["delta_mean"]),
    ))
    a("")
    a("NHEJ core Δ = {}. KEGG NHEJ Δ = {}. STING core Δ = {}. Compact IFN Δ = {}. Hallmark IFN-α Δ = {}. Hallmark IFN-γ Δ = {}. MHC-I APM Δ = {}.".format(
        fmt_d(mods["NHEJ_CORE"]["delta_mean"]),
        fmt_d(mods["NHEJ_KEGG"]["delta_mean"]),
        fmt_d(mods["STING_CORE"]["delta_mean"]),
        fmt_d(mods["IFN_COMPACT"]["delta_mean"]),
        fmt_d(mods["IFN_HALLMARK_ALPHA"]["delta_mean"]),
        fmt_d(mods["IFN_HALLMARK_GAMMA"]["delta_mean"]),
        fmt_d(mods["APM_MHCI"]["delta_mean"]),
    ))
    a("")
    a("## Cldn4 neighborhood")
    a("")
    a(neigh["prose"])
    a("")
    a("## Sweep: KL>KP on Cldn4 and Tacstd2, with NHEJ lower and IFN higher")
    a("")
    a(sweep["prose"])
    a("")
    a("## What a later mouse scRNA score can use")
    a("")
    a("- `tables/signatures.gmt` — frozen lists.")
    a("- `tables/signature_genes.tsv` — one row per gene per list, including symbols that were absent from this matrix.")
    a("- `tables/de_kl_vs_kp.tsv.gz` — every gene, with Δ log2(FPKM+1), Δ log2(normalized+1), Δ log2(CPM+1), Δ log2(median-of-ratios+1), exact MW p, Cliff’s δ, and Welch p. Benjamini–Hochberg columns are on the full gene table. The MW p-values pile up at 0.00794, so that FDR is not a discovery list for Tacstd2 or Cldn4.")
    a("- `tables/cldn4_spearman.tsv.gz` — within-arm and 10-library Spearman with Cldn4, so the neighborhood threshold can be moved without re-downloading GEO.")
    a("- `tables/sweep_pattern.tsv` — every scale × mean/median × cohort × NHEJ list × IFN list.")
    a("")
    a("Library-level p = 0.00794 is the floor for a 5-vs-5 rank test. It is not an independent-line p-value when libraries inside an arm are near-replicates. Cell-line IFN and APM are transcripts in cultured lines. They are not an immune-exclusion measurement.")
    a("")
    (HERE / "FINDING.md").write_text("\n".join(lines))


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    hashes = download()
    raw = load_matrix(DATA / "GSE137244_counts.raw.csv.gz", quoted_rows=True)
    norm = align_symbol_index(raw, load_matrix(DATA / "GSE137244_counts.normalized.csv.gz", quoted_rows=False), "normalized")
    fpkm = align_symbol_index(raw, load_matrix(DATA / "GSE137244_counts.fpkm.csv.gz", quoted_rows=False), "fpkm")
    for name, df in ("raw", raw), ("normalized", norm), ("fpkm", fpkm):
        missing = [c for c in KP_LIBS + KL_LIBS + [NORMAL] if c not in df.columns]
        if missing:
            raise SystemExit(f"{name} missing columns {missing}; got {list(df.columns)}")
        if df.index.duplicated().any():
            raise SystemExit(f"duplicate genes in {name}")

    tumor = KP_LIBS + KL_LIBS
    log_fpkm = log2p1(fpkm)
    log_norm = log2p1(norm)
    log_cpm = log2p1(cpm(raw))
    log_mor = log2p1(median_of_ratios(raw[tumor]))
    scales = {
        "fpkm": log_fpkm,
        "normalized": log_norm,
        "cpm": log_cpm,
        "median_ratio": log_mor,
    }

    # Locked-gene integrity check on the FPKM scale.
    locked = {}
    for gene in ("Tacstd2", "Cldn4", "Stk11", "Trp53", "Epcam"):
        locked[gene] = contrast_vector(log_fpkm.loc[gene], KL_LIBS, KP_LIBS)
    if abs(locked["Cldn4"]["delta_mean"] - 5.570) > 0.002:
        raise SystemExit(f"Cldn4 delta drifted: {locked['Cldn4']['delta_mean']}")
    if abs(locked["Tacstd2"]["delta_mean"] - 3.238) > 0.002:
        raise SystemExit(f"Tacstd2 delta drifted: {locked['Tacstd2']['delta_mean']}")
    if abs(locked["Cldn4"]["mw_p"] - 2 / 252) > 1e-9:
        raise SystemExit(f"Cldn4 MW p drifted: {locked['Cldn4']['mw_p']}")

    matched, membership = build_sets(list(fpkm.index))
    neigh_table, neigh_info = neighborhood(log_fpkm, fpkm, KL_LIBS, KP_LIBS)
    neigh_genes = neigh_info["chosen"]["genes"]
    matched["CLDN4_NEIGHBORHOOD"] = neigh_genes
    for gene in neigh_genes:
        membership = pd.concat(
            [
                membership,
                pd.DataFrame(
                    [
                        {
                            "set_name": "CLDN4_NEIGHBORHOOD",
                            "gene": gene,
                            "status": "used",
                            "source": (
                                f"within-arm Spearman >= {neigh_info['chosen']['threshold']}, "
                                f"both arms > 0, mean FPKM >= {NEIGH_MIN_MEAN_FPKM}"
                            ),
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

    # Gene-level DE on all four scales. Welch and MW on the FPKM scale.
    de_rows = {
        "gene": fpkm.index,
        "mean_fpkm_kl": fpkm[KL_LIBS].mean(axis=1).to_numpy(),
        "mean_fpkm_kp": fpkm[KP_LIBS].mean(axis=1).to_numpy(),
    }
    for scale_name, mat in scales.items():
        kl = mat[KL_LIBS]
        kp = mat[KP_LIBS]
        # median_ratio has no normal-lung column; tumor columns only.
        de_rows[f"delta_{scale_name}"] = (kl.mean(axis=1) - kp.mean(axis=1)).to_numpy()
    de = pd.DataFrame(de_rows).set_index("gene")
    # Exact MW for every gene is 52k * C(10,5) via scipy. Doable but slow if pure python.
    # Vectorized Cliff plus a rank-sum p from the normal approximation is not exact.
    # Exact enumeration of the 252 assignments is fast in numpy.
    vals = log_fpkm[tumor].to_numpy(dtype=float)
    # Precompute all subset mean differences' rank statistics.
    # Mann-Whitney U = number of KL>KP pairs, computable from ranks of the 10 values.
    # For ties, scipy's exact method is the reference. At this n, looping scipy is ~1-2 min.
    # Use a vectorized pair count (handles ties as not gt / not lt) and map the
    # (n_gt, n_lt) pair onto the exact MW p by calling scipy once per unique pair-count.
    kl_block = vals[:, 5:]
    kp_block = vals[:, :5]
    n_gt = np.zeros(vals.shape[0], dtype=np.int16)
    n_lt = np.zeros(vals.shape[0], dtype=np.int16)
    for i in range(5):
        n_gt += (kl_block > kp_block[:, i : i + 1]).sum(axis=1).astype(np.int16)
        n_lt += (kl_block < kp_block[:, i : i + 1]).sum(axis=1).astype(np.int16)
    de["n_kl_gt_kp"] = n_gt
    de["n_kl_lt_kp"] = n_lt
    de["cliff_fpkm"] = (n_gt.astype(float) - n_lt.astype(float)) / 25.0
    # Welch on log2(FPKM+1), vectorized.
    kl_m = kl_block.mean(axis=1)
    kp_m = kp_block.mean(axis=1)
    kl_v = kl_block.var(axis=1, ddof=1)
    kp_v = kp_block.var(axis=1, ddof=1)
    se2 = kl_v / 5 + kp_v / 5
    tstat = np.full(vals.shape[0], np.nan)
    df_w = np.full(vals.shape[0], np.nan)
    ok = se2 > 0
    tstat[ok] = (kl_m[ok] - kp_m[ok]) / np.sqrt(se2[ok])
    df_w[ok] = (se2[ok] ** 2) / ((kl_v[ok] / 5) ** 2 / 4 + (kp_v[ok] / 5) ** 2 / 4)
    welch_p = np.full(vals.shape[0], np.nan)
    welch_p[ok] = 2 * stats.t.sf(np.abs(tstat[ok]), df_w[ok])
    de["welch_t_fpkm"] = tstat
    de["welch_p_fpkm"] = welch_p

    unique_pairs = sorted(set(zip(n_gt.tolist(), n_lt.tolist())))
    p_of = {}
    for a, b in unique_pairs:
        # Reconstruct two length-5 samples with the same pairwise order counts
        # is unnecessary: scipy MW depends on the values. For untied data U determines p.
        # With ties, p depends on the tie pattern. Compute exact p from a synthetic
        # sample only when there are no ties; otherwise call scipy on the real row
        # for that tie class's first member. To stay exact, evaluate scipy on one
        # real row per (n_gt, n_lt, n_tie) AND verify... ties make U not sufficient.
        p_of[(a, b)] = None
    # Exact p: group genes by whether they have ties. Untied: U -> p via one scipy call
    # on a canonical 0..9 ranking. Tied: compute per gene. Tied genes are common at zeros.
    n_tie = 25 - n_gt - n_lt
    untied = n_tie == 0
    # Canonical complete-separation and other untied U values.
    # Build a lookup by placing KL on the top n_gt unique orderings is hard.
    # Faster path: for untied data, ranks are a permutation and U = n_gt.
    # scipy mannwhitneyu of a fixed pair of samples with that many wins.
    canonical_p = {}
    for wins in range(0, 26):
        # Construct KL and KP such that exactly `wins` pairwise KL>KP, no ties.
        # Use all pairwise: give the 10 distinct ranks, assign a 5-subset with U=wins.
        found = None
        ranks = np.arange(10, dtype=float)
        for subset in combinations(range(10), 5):
            kl_r = ranks[list(subset)]
            kp_r = ranks[[i for i in range(10) if i not in subset]]
            if int(np.sum(kl_r[:, None] > kp_r[None, :])) == wins:
                found = (kl_r, kp_r)
                break
        if found is None:
            continue
        canonical_p[wins] = float(
            stats.mannwhitneyu(found[0], found[1], alternative="two-sided", method="exact").pvalue
        )
    mw_p = np.full(vals.shape[0], np.nan)
    for wins, p in canonical_p.items():
        mw_p[(n_gt == wins) & untied] = p
    tied_idx = np.flatnonzero(~untied)
    # Tied rows: exact test per gene. Cache by the multiset of 10 values rounded
    # is unsafe. Cache by the rank tuple.
    cache: dict[tuple[int, ...], float] = {}
    rank_mat = np.apply_along_axis(stats.rankdata, 1, vals)
    for i in tied_idx:
        key = tuple(np.round(rank_mat[i], 5))
        # Include which positions are KL by storing ranks in sample order.
        if key not in cache:
            cache[key] = float(
                stats.mannwhitneyu(kl_block[i], kp_block[i], alternative="two-sided", method="exact").pvalue
            )
        mw_p[i] = cache[key]
    de["mw_p_fpkm"] = mw_p
    de["mw_fdr_bh"] = bh(pd.Series(mw_p, index=de.index)).to_numpy()
    de["welch_fdr_bh"] = bh(pd.Series(welch_p, index=de.index)).to_numpy()
    de["spearman_cldn4_all10"] = np.nan
    de["spearman_within_mean"] = np.nan
    de.loc[neigh_table.index, "spearman_cldn4_all10"] = neigh_table["spearman_all10"]
    de.loc[neigh_table.index, "spearman_within_mean"] = neigh_table["spearman_within_mean"]
    de = de.reset_index().rename(columns={"index": "gene"})
    de_path = TABLES / "de_kl_vs_kp.tsv.gz"
    de.to_csv(de_path, sep="\t", index=False, compression="gzip")

    # Module contrasts on the primary scale, all 10 libraries, mean score.
    module_rows = []
    score_rows = []
    for name, genes in matched.items():
        if not genes:
            module_rows.append(
                {
                    "set_name": name,
                    "n_genes": 0,
                    "delta_mean": np.nan,
                    "mw_p": np.nan,
                    "cliff": np.nan,
                    "separation": "empty",
                    "welch_p": np.nan,
                }
            )
            continue
        score = mean_score(log_fpkm, genes, tumor)
        rec = contrast_vector(score, KL_LIBS, KP_LIBS)
        rec["set_name"] = name
        rec["n_genes"] = len(genes)
        module_rows.append(rec)
        for lib in tumor:
            score_rows.append(
                {
                    "set_name": name,
                    "library": lib,
                    "genotype": "KL" if lib in KL_LIBS else "KP",
                    "score_mean_log2fpkm1": float(score.loc[lib]),
                }
            )
    modules = pd.DataFrame(module_rows)
    scores = pd.DataFrame(score_rows)
    modules.to_csv(TABLES / "module_contrasts.tsv", sep="\t", index=False)
    scores.to_csv(TABLES / "module_scores.tsv", sep="\t", index=False)

    # Sweep.
    cohorts = {
        "all_10": (KL_LIBS, KP_LIBS),
        "drop_B6AL10-3": (KL_LIBS, [c for c in KP_LIBS if c != EPCAM_DROP]),
    }
    sweep_rows = []
    for scale_name, mat in scales.items():
        for summary in ("mean", "median"):
            scorer = mean_score if summary == "mean" else median_score
            for cohort, (kl_libs, kp_libs) in cohorts.items():
                cols = kp_libs + kl_libs
                cldn4 = contrast_vector(mat.loc["Cldn4", cols], kl_libs, kp_libs)
                tacstd2 = contrast_vector(mat.loc["Tacstd2", cols], kl_libs, kp_libs)
                for nhej_name in NHEJ_SWEEP:
                    nhej_genes = matched[nhej_name]
                    nhej_score = scorer(mat, nhej_genes, cols)
                    nhej = contrast_vector(nhej_score, kl_libs, kp_libs)
                    for ifn_name in IFN_SWEEP:
                        ifn_genes = matched[ifn_name]
                        ifn_score = scorer(mat, ifn_genes, cols)
                        ifn = contrast_vector(ifn_score, kl_libs, kp_libs)
                        sweep_rows.append(
                            pattern_row(
                                scale_name,
                                summary,
                                cohort,
                                nhej_name,
                                ifn_name,
                                cldn4,
                                tacstd2,
                                nhej,
                                ifn,
                                len(nhej_genes),
                                len(ifn_genes),
                            )
                        )
    sweep = pd.DataFrame(sweep_rows)
    sweep = sweep.sort_values(
        ["hit_nhej_down_ifn_up", "n_limbs", "strength_cliff"],
        ascending=[False, False, False],
    )
    sweep.to_csv(TABLES / "sweep_pattern.tsv", sep="\t", index=False)
    n_hit = int(sweep["hit_nhej_down_ifn_up"].sum())
    best = sweep.iloc[0].to_dict()
    primary_mask = (
        (sweep["scale"] == "fpkm")
        & (sweep["summary"] == "mean")
        & (sweep["cohort"] == "all_10")
        & (sweep["nhej_set"] == "NHEJ_CORE")
        & (sweep["ifn_set"] == "IFN_COMPACT")
    )
    primary_sweep = sweep.loc[primary_mask].iloc[0].to_dict()
    nhej_pos = int((sweep["nhej_delta"] > 0).sum())
    ifn_pos = int((sweep["ifn_delta"] > 0).sum())

    # Within-arm correlation of libraries.
    expressed = fpkm[tumor].mean(axis=1) >= 1
    corr = log_fpkm.loc[expressed, tumor].corr()
    within = []
    for group in (KP_LIBS, KL_LIBS):
        for a, b in combinations(group, 2):
            within.append(float(corr.loc[a, b]))
    max_r = max(within)

    # Annotation.
    ann_rows = []
    for lib in KP_LIBS:
        ann_rows.append(
            {
                "library": lib,
                "genotype": "KP",
                "genotype_source": "GEO Sample_description: KrasLSL-G12D crossed with p53 conditional knockout",
                "in_primary": "yes",
            }
        )
    for lib in KL_LIBS:
        ann_rows.append(
            {
                "library": lib,
                "genotype": "KL",
                "genotype_source": "GEO Sample_description: KrasLSL-G12D crossed with Lkb1",
                "in_primary": "yes",
            }
        )
    ann_rows.append(
        {
            "library": NORMAL,
            "genotype": "normal_lung",
            "genotype_source": "GEO title normal-lung-RNA; tissue Lung; excluded from KL vs KP",
            "in_primary": "no",
        }
    )
    pd.DataFrame(ann_rows).to_csv(TABLES / "sample_annotation.tsv", sep="\t", index=False)
    membership.to_csv(TABLES / "signature_genes.tsv", sep="\t", index=False)
    pd.DataFrame(neigh_info["thresholds"]).to_csv(TABLES / "cldn4_neighborhood_thresholds.tsv", sep="\t", index=False)
    neigh_table.reset_index().to_csv(TABLES / "cldn4_spearman.tsv.gz", sep="\t", index=False, compression="gzip")

    neigh_note = (
        f"Within-arm Spearman with Cldn4, mean of KL and KP arms >= {neigh_info['chosen']['threshold']}, "
        "both arms positive, mean FPKM >= 1, Cldn4 excluded. "
        f"n={len(neigh_genes)}. Operating rule: highest threshold in 0.9,0.8,0.7,0.6,0.5 with at least {NEIGH_MIN_GENES} genes."
    )
    write_gmt(matched, neigh_genes, neigh_note)

    # Neighborhood prose pieces.
    chosen_thr = neigh_info["chosen"]["threshold"]
    neigh_mod = modules[modules["set_name"] == "CLDN4_NEIGHBORHOOD"].iloc[0]
    tj_in = sorted(set(neigh_genes) & set(matched["TJ_EPITHELIAL"]))
    tj_out = sorted(set(matched["TJ_EPITHELIAL"]) - set(neigh_genes) - {"Cldn4"})
    if neigh_info["chosen"]["below_min"]:
        size_sentence = (
            f"No threshold retained {NEIGH_MIN_GENES} genes. The exported list is the loosest grid point, "
            f"ρ ≥ {chosen_thr}, with {len(neigh_genes)} genes."
        )
    else:
        size_sentence = (
            f"The operating threshold is the highest grid point with at least {NEIGH_MIN_GENES} genes: "
            f"ρ ≥ {chosen_thr}, {len(neigh_genes)} genes."
        )
    neigh_prose = (
        "The co-expression list is built inside each genotype, on log2(FPKM+1), so the KL-versus-KP axis does not pick the genes. "
        "A gene qualifies when mean FPKM is at least 1 and the Spearman with Cldn4 is positive in both arms, "
        "with the mean of the two arm correlations at or above the threshold. Cldn4 itself is excluded. "
        + size_sentence
        + f" The mean score is Δ {fmt_d(float(neigh_mod['delta_mean']))}, exact p {fmt_p(float(neigh_mod['mw_p']))}, {neigh_mod['separation']}. "
        + f"Epithelial TJ genes inside this list: {', '.join(tj_in) if tj_in else 'none'}. "
        + f"Epithelial TJ genes outside it: {', '.join(tj_out)}. "
        "Within-arm correlations at n = 5 are unstable, so this list is the co-expression export, "
        "while the vendored 11-gene edge and the 18-gene epithelial core are the lists whose mean scores separate KL from KP."
    )

    if n_hit == 0:
        sweep_prose = (
            f"The grid has {len(sweep)} cells: 4 scales (log2 FPKM+1, deposited normalized counts, CPM, median-of-ratios) "
            f"× mean or median × 2 cohorts (all 10 libraries, or KP without B6AL10-3) "
            f"× 3 NHEJ lists × 4 IFN lists. "
            f"**{n_hit} cells** have Cldn4 higher, Tacstd2 higher, NHEJ lower, and IFN higher at the same time. "
            f"NHEJ Δ is positive in {nhej_pos}/{len(sweep)} cells. IFN Δ is positive in {ifn_pos}/{len(sweep)} cells. "
            f"The pre-specified primary cell (FPKM, mean, all 10, NHEJ core, compact IFN) has "
            f"Cldn4 {fmt_d(primary_sweep['cldn4_delta'])}, Tacstd2 {fmt_d(primary_sweep['tacstd2_delta'])}, "
            f"NHEJ {fmt_d(primary_sweep['nhej_delta'])}, IFN {fmt_d(primary_sweep['ifn_delta'])} "
            f"({int(primary_sweep['n_limbs'])} of 4 limbs in the requested direction). "
            f"NHEJ Δ is negative in {int((sweep['nhej_delta'] < 0).sum())} cells, all of them median summaries, "
            f"and IFN Δ is negative in every one of those cells. "
            f"The highest Cliff strength on the grid is {best['strength_cliff']:.3f} "
            f"({best['scale']}, {best['summary']}, {best['cohort']}, {best['nhej_set']}, {best['ifn_set']}; "
            f"NHEJ Δ {fmt_d(best['nhej_delta'])}, IFN Δ {fmt_d(best['ifn_delta'])}, "
            f"{int(best['n_limbs'])} of 4 limbs). "
            "That cell is the top of a ranking. It is not a selected signature."
        )
    else:
        sweep_prose = (
            f"The grid has {len(sweep)} cells. **{n_hit} cells** put all four limbs in the requested direction. "
            f"The strongest of those, by the sum of Cliff’s δ in the requested direction, is "
            f"{best['scale']}, {best['summary']}, {best['cohort']}, {best['nhej_set']}, {best['ifn_set']}: "
            f"Cldn4 {fmt_d(best['cldn4_delta'])}, Tacstd2 {fmt_d(best['tacstd2_delta'])}, "
            f"NHEJ {fmt_d(best['nhej_delta'])}, IFN {fmt_d(best['ifn_delta'])}, "
            f"strength {best['strength_cliff']:.3f}. "
            f"The pre-specified primary cell remains FPKM / mean / all 10 / NHEJ core / compact IFN: "
            f"NHEJ {fmt_d(primary_sweep['nhej_delta'])}, IFN {fmt_d(primary_sweep['ifn_delta'])}."
        )

    primary_modules = {rec["set_name"]: rec for rec in module_rows}
    ctx = {
        "locked": locked,
        "primary_modules": primary_modules,
        "max_within_arm_r": max_r,
        "neighborhood": {"prose": neigh_prose, "n": len(neigh_genes), "threshold": chosen_thr},
        "sweep_summary": {"prose": sweep_prose, "n_hit": n_hit, "best": best, "primary": primary_sweep},
    }
    write_finding(ctx)

    summary = {
        "geo": "GSE137244",
        "sha256": hashes,
        "n_kp": 5,
        "n_kl": 5,
        "locked": {g: {k: locked[g][k] for k in ("delta_mean", "mw_p", "separation")} for g in locked},
        "neighborhood": {
            "threshold": chosen_thr,
            "n_genes": len(neigh_genes),
            "genes": neigh_genes,
            "below_min": neigh_info["chosen"]["below_min"],
        },
        "n_sweep": len(sweep),
        "n_hit_nhej_down_ifn_up": n_hit,
        "best_sweep": best,
        "primary_sweep": primary_sweep,
        "max_within_arm_pearson": max_r,
        "modules": {
            name: {
                "n_genes": int(primary_modules[name]["n_genes"]),
                "delta_mean": primary_modules[name]["delta_mean"],
                "mw_p": primary_modules[name]["mw_p"],
                "separation": primary_modules[name]["separation"],
            }
            for name in primary_modules
        },
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    plot_modules(scores, modules)
    edge_genes = [g for g in ["Tacstd2", *gs.TJ_EPITHELIAL] if g in log_fpkm.index]
    edge_delta = pd.Series(
        {
            g: float(log_fpkm.loc[g, KL_LIBS].mean() - log_fpkm.loc[g, KP_LIBS].mean())
            for g in edge_genes
        }
    )
    plot_sweep(sweep, edge_delta)
    print(json.dumps({
        "n_hit": n_hit,
        "n_sweep": len(sweep),
        "neigh_n": len(neigh_genes),
        "neigh_thr": chosen_thr,
        "Cldn4": locked["Cldn4"]["delta_mean"],
        "Tacstd2": locked["Tacstd2"]["delta_mean"],
        "NHEJ": primary_modules["NHEJ_CORE"]["delta_mean"],
        "IFN": primary_modules["IFN_COMPACT"]["delta_mean"],
        "APM": primary_modules["APM_MHCI"]["delta_mean"],
        "max_r": max_r,
        "best": {
            "scale": best["scale"],
            "summary": best["summary"],
            "cohort": best["cohort"],
            "nhej": best["nhej_set"],
            "ifn": best["ifn_set"],
            "n_limbs": int(best["n_limbs"]),
            "nhej_delta": best["nhej_delta"],
            "ifn_delta": best["ifn_delta"],
        },
    }, indent=2))


if __name__ == "__main__":
    main()
