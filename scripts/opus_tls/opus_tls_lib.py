"""Shared definitions and statistics helpers for the TACSTD2/CLDN4 vs TLS analysis.

Everything here is deliberately dependency-light (numpy / pandas / scipy /
statsmodels / lifelines) so the pipeline runs on a bare python install.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = os.environ.get("OPUS_TLS_DATA", "/tmp/opus_tls_data")
REPO_DIR = os.environ.get(
    "OPUS_TLS_REPO", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
RESULTS_DIR = os.path.join(REPO_DIR, "results", "opus_tls")
PROC_DIR = os.path.join(DATA_DIR, "processed")

for _d in (RESULTS_DIR, PROC_DIR, os.path.join(RESULTS_DIR, "tables"), os.path.join(RESULTS_DIR, "figures")):
    os.makedirs(_d, exist_ok=True)

RNG_SEED = 20260816

# ---------------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------------
# TLS / B-lineage signatures. Sources are given in notes/opus_tls/METHODS.md.

SIGNATURES: dict[str, list[str]] = {
    # Cabrita et al. Nature 2020 melanoma TLS signature.
    "TLS_Cabrita": ["CCL19", "CCL21", "CXCL13", "CCL17", "CCR7", "CXCR5", "SELL", "LAMP3"],
    # Coppola/Messina 12-chemokine TLS signature.
    "TLS_12chemokine": [
        "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18",
        "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13",
    ],
    # Meylan et al. Immunity 2022 / Petitprez et al. Nature 2020 "TLS imprint"
    # (B-cell + follicular helper + follicular dendritic cell + HEV programme).
    "TLS_imprint": [
        "CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS",
    ],
    # B-lineage (MCP-counter B lineage core + canonical B markers).
    "B_cell": [
        "MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "BLK", "CD22", "CR2",
        "FCRL2", "FCRL5", "PAX5", "TCL1A", "TNFRSF13C", "VPREB3", "CD37", "IRF8",
    ],
    # Plasma-cell / immunoglobulin programme.
    "Plasma_cell": [
        "MZB1", "JCHAIN", "DERL3", "TNFRSF17", "POU2AF1", "XBP1", "PRDM1",
        "SDC1", "IGHG1", "IGHG3", "IGKC", "IGHA1", "SSR4", "FKBP11",
    ],
    # Germinal-centre reaction (dark/light zone).
    "GC_reaction": ["AICDA", "BCL6", "RGS13", "MEF2B", "LMO2", "MKI67", "STMN1", "CD38"],
    # T follicular helper.
    "Tfh": ["CXCL13", "CD200", "BTLA", "ICOS", "PDCD1", "CD40LG", "MAF", "IL21", "TOX2"],
    # High endothelial venules / lymphoid stroma.
    "HEV_stroma": ["CHST4", "GLYCAM1", "MADCAM1", "ACKR1", "SELP", "LTB", "LTA", "TNFSF14"],
    # Positive control immune programmes (should also fall with tumour purity).
    "T_cell_CD8": ["CD8A", "CD8B", "GZMK", "GZMA", "PRF1", "CD2", "CD3D", "CD3E"],
    "IFNg_Ayers": ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"],
    "Myeloid": ["CD68", "CD163", "CSF1R", "ITGAM", "AIF1", "MRC1", "FCGR3A"],
    # Tumour-cell content proxy (pan-epithelial).
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"],
}

TLS_SIGNATURES = [
    "TLS_Cabrita", "TLS_12chemokine", "TLS_imprint", "B_cell",
    "Plasma_cell", "GC_reaction", "Tfh", "HEV_stroma",
]
CONTROL_SIGNATURES = ["T_cell_CD8", "IFNg_Ayers", "Myeloid", "Epithelial"]

# Genes of interest: primary hypothesis genes + epithelial/ADC-target comparators.
FOCUS_GENES = ["TACSTD2", "CLDN4"]
COMPARATOR_GENES = [
    "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "CLDN3", "CLDN7",
    "MUC1", "ERBB2", "FOLR1", "CEACAM5", "MET", "NECTIN4", "CD274",
]


@dataclass
class Cohort:
    """A harmonised expression cohort."""

    name: str
    expr: pd.DataFrame  # genes x samples, log2 scale
    pheno: pd.DataFrame  # samples x variables
    kind: str  # "atlas" or "ici"
    platform: str
    note: str = ""

    def __post_init__(self) -> None:
        shared = [s for s in self.expr.columns if s in self.pheno.index]
        self.expr = self.expr[shared]
        self.pheno = self.pheno.loc[shared]

    @property
    def n(self) -> int:
        return self.expr.shape[1]


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def signature_scores(expr: pd.DataFrame, signatures: dict[str, list[str]] | None = None,
                     min_genes: int = 3) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Mean-of-z signature scores. Returns (samples x signature, genes used)."""
    signatures = signatures or SIGNATURES
    z = zscore_rows(expr)
    out, used = {}, {}
    for name, genes in signatures.items():
        present = [g for g in genes if g in z.index]
        if len(present) < min_genes:
            continue
        out[name] = z.loc[present].mean(axis=0)
        used[name] = present
    return pd.DataFrame(out), used


def rank_signature_scores(expr: pd.DataFrame, signatures: dict[str, list[str]] | None = None,
                          min_genes: int = 3) -> pd.DataFrame:
    """Sensitivity scoring: mean within-sample percentile rank of signature genes."""
    signatures = signatures or SIGNATURES
    ranks = expr.rank(axis=0, pct=True)
    out = {}
    for name, genes in signatures.items():
        present = [g for g in genes if g in ranks.index]
        if len(present) < min_genes:
            continue
        out[name] = ranks.loc[present].mean(axis=0)
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    q[ok] = out
    return q


def fisher_z(r: float) -> float:
    r = float(np.clip(r, -0.999999, 0.999999))
    return np.arctanh(r)


def spearman_ci(x, y, n_boot: int = 2000, seed: int = RNG_SEED) -> tuple[float, float, float, int]:
    """Spearman rho with p-value and bootstrap 95% CI."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = x.size
    if n < 8:
        return np.nan, np.nan, np.nan, n
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(x[idx]).size < 3 or np.unique(y[idx]).size < 3:
            boots[i] = np.nan
            continue
        boots[i] = stats.spearmanr(x[idx], y[idx])[0]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return rho, p, (lo, hi), n


def partial_spearman(x, y, covars: pd.DataFrame | np.ndarray) -> tuple[float, float, int]:
    """Spearman partial correlation: rank-transform, then correlate OLS residuals.

    p-value from the t-statistic of the residual Pearson correlation with
    n - 2 - k degrees of freedom (k = number of covariates).
    """
    x = pd.Series(np.asarray(x, float))
    y = pd.Series(np.asarray(y, float))
    C = pd.DataFrame(covars).reset_index(drop=True)
    df = pd.concat([x.rename("x").reset_index(drop=True), y.rename("y").reset_index(drop=True), C], axis=1)
    df = df.dropna()
    if df.shape[0] < 10 + C.shape[1]:
        return np.nan, np.nan, df.shape[0]
    n = df.shape[0]
    xr = stats.rankdata(df["x"].values)
    yr = stats.rankdata(df["y"].values)
    Cr = np.column_stack([stats.rankdata(df[c].values) for c in C.columns])
    Cr = np.column_stack([np.ones(n), Cr])
    bx, *_ = np.linalg.lstsq(Cr, xr, rcond=None)
    by, *_ = np.linalg.lstsq(Cr, yr, rcond=None)
    rx = xr - Cr @ bx
    ry = yr - Cr @ by
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan, np.nan, n
    r = float(np.corrcoef(rx, ry)[0, 1])
    k = Cr.shape[1] - 1
    dof = n - 2 - k
    if dof <= 0:
        return r, np.nan, n
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), dof)
    return r, float(p), n


def mann_whitney_effect(a, b) -> dict:
    """Mann-Whitney U with rank-biserial correlation and Hodges-Lehmann shift."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return {"n1": a.size, "n2": b.size, "U": np.nan, "p": np.nan, "rbc": np.nan, "hl_shift": np.nan}
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rbc = 2 * U / (a.size * b.size) - 1  # >0 means group a larger
    hl = float(np.median(np.subtract.outer(a, b)))
    return {"n1": int(a.size), "n2": int(b.size), "U": float(U), "p": float(p),
            "rbc": float(rbc), "hl_shift": hl,
            "median1": float(np.median(a)), "median2": float(np.median(b))}


def meta_fisher_z(rhos: list[float], ns: list[int]) -> dict:
    """Fixed- and random-effects (DerSimonian-Laird) meta-analysis of correlations."""
    rhos = np.asarray(rhos, float)
    ns = np.asarray(ns, float)
    ok = np.isfinite(rhos) & (ns > 4)
    rhos, ns = rhos[ok], ns[ok]
    if rhos.size == 0:
        return {}
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    # Variance of Fisher z for Spearman (Bonett & Wright 2000): (1 + rho^2/2)/(n-3)
    v = (1 + rhos ** 2 / 2) / (ns - 3)
    w = 1 / v
    z_fe = float(np.sum(w * z) / np.sum(w))
    se_fe = float(np.sqrt(1 / np.sum(w)))
    Q = float(np.sum(w * (z - z_fe) ** 2))
    dof = rhos.size - 1
    C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - dof) / C) if C > 0 else 0.0
    wr = 1 / (v + tau2)
    z_re = float(np.sum(wr * z) / np.sum(wr))
    se_re = float(np.sqrt(1 / np.sum(wr)))
    I2 = max(0.0, (Q - dof) / Q * 100) if Q > 0 else 0.0
    return {
        "k": int(rhos.size),
        "n_total": int(ns.sum()),
        "rho_fe": float(np.tanh(z_fe)),
        "p_fe": float(2 * stats.norm.sf(abs(z_fe / se_fe))),
        "ci_fe": (float(np.tanh(z_fe - 1.96 * se_fe)), float(np.tanh(z_fe + 1.96 * se_fe))),
        "rho_re": float(np.tanh(z_re)),
        "p_re": float(2 * stats.norm.sf(abs(z_re / se_re))),
        "ci_re": (float(np.tanh(z_re - 1.96 * se_re)), float(np.tanh(z_re + 1.96 * se_re))),
        "Q": Q, "p_Q": float(stats.chi2.sf(Q, dof)) if dof > 0 else np.nan,
        "I2": I2, "tau2": tau2,
    }


def genomewide_rank(target_rho: float, all_rhos: pd.Series) -> dict:
    """Where does the target gene sit in the transcriptome-wide null of
    correlations with the same outcome? Two-sided empirical p."""
    vals = all_rhos.dropna().values
    if vals.size == 0 or not np.isfinite(target_rho):
        return {}
    frac_below = float(np.mean(vals <= target_rho))
    frac_extreme = float(np.mean(np.abs(vals) >= abs(target_rho)))
    return {
        "percentile": 100 * frac_below,
        "rank_from_bottom": int(np.sum(vals <= target_rho)),
        "n_genes": int(vals.size),
        "p_empirical_two_sided": max(frac_extreme, 1 / vals.size),
    }


def ols_table(y: pd.Series, X: pd.DataFrame, label: str = "") -> pd.DataFrame:
    """Multivariable OLS on rank-normalised y with HC3 robust SEs."""
    import statsmodels.api as sm

    df = pd.concat([y.rename("_y"), X], axis=1).dropna()
    if df.shape[0] < 20:
        return pd.DataFrame()
    yv = df["_y"].values
    Xv = df.drop(columns="_y")
    Xv = sm.add_constant(Xv.astype(float), has_constant="add")
    model = sm.OLS(yv, Xv).fit(cov_type="HC3")
    out = pd.DataFrame({
        "term": model.params.index,
        "beta": model.params.values,
        "se": model.bse.values,
        "t": model.tvalues.values,
        "p": model.pvalues.values,
        "ci_low": model.conf_int()[0].values,
        "ci_high": model.conf_int()[1].values,
    })
    out["n"] = df.shape[0]
    out["r2"] = model.rsquared
    out["model"] = label
    return out


def inverse_normal(x: pd.Series) -> pd.Series:
    """Rank-based inverse normal transform (Blom), robust to skew/outliers."""
    r = x.rank(method="average")
    n = r.notna().sum()
    return pd.Series(stats.norm.ppf((r - 0.375) / (n + 0.25)), index=x.index)


def tertile_quartile_split(v: pd.Series, q: int = 4) -> pd.Series:
    try:
        return pd.qcut(v, q, labels=[f"Q{i+1}" for i in range(q)])
    except ValueError:
        return pd.Series(index=v.index, dtype="object")


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"
