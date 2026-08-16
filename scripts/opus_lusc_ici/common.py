"""Shared helpers: GEO parsing, expression normalisation and small-sample statistics."""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from config import RAW_DIR, SCRIPTS_DIR

# --------------------------------------------------------------------------
# GEO series-matrix parsing
# --------------------------------------------------------------------------


def read_series_matrix(gse: str) -> pd.DataFrame:
    """Return one row per GSM with title and every ``characteristics`` field."""
    path = RAW_DIR / f"{gse}_series_matrix.txt.gz"
    with gzip.open(path, "rt", encoding="utf8", errors="replace") as fh:
        lines = fh.read().split("\n")

    titles, gsms, chars = None, None, {}
    for line in lines:
        if line.startswith("!Sample_title"):
            titles = re.findall(r'"([^"]*)"', line)
        elif line.startswith("!Sample_geo_accession"):
            gsms = re.findall(r'"([^"]*)"', line)
        elif line.startswith("!Sample_characteristics_ch1"):
            values = re.findall(r'"([^"]*)"', line)
            key = values[0].split(":", 1)[0].strip() if ":" in values[0] else "unknown"
            parsed = [v.split(":", 1)[1].strip() if ":" in v else np.nan for v in values]
            # Repeated keys (GEO allows this) get a numeric suffix.
            base, n = key, 1
            while key in chars:
                n += 1
                key = f"{base}_{n}"
            chars[key] = parsed

    frame = pd.DataFrame({"gsm": gsms, "title": titles})
    for key, values in chars.items():
        frame[key] = values
    return frame


# --------------------------------------------------------------------------
# Expression normalisation
# --------------------------------------------------------------------------


def load_probemap() -> pd.Series:
    """Ensembl gene id (version stripped) -> HGNC symbol, from the Xena GDC probemap."""
    path = RAW_DIR / "gencode.v36.annotation.gtf.gene.probemap"
    frame = pd.read_csv(path, sep="\t", usecols=["id", "gene"])
    frame["ens"] = frame["id"].str.split(".").str[0]
    return frame.drop_duplicates("ens").set_index("ens")["gene"]


def collapse_to_symbols(expr: pd.DataFrame, mapping: pd.Series | None = None) -> pd.DataFrame:
    """Map row ids to gene symbols and keep, per symbol, the most expressed row.

    ``expr`` is genes x samples on a log scale (or any scale where the row mean
    orders rows sensibly).
    """
    if mapping is not None:
        index = expr.index.to_series().str.split(".").str[0].map(mapping)
        expr = expr.loc[index.notna()]
        expr.index = index.dropna().values
    expr = expr[~expr.index.isna()]
    if expr.index.has_duplicates:
        order = expr.mean(axis=1).groupby(level=0).rank(ascending=False, method="first")
        expr = expr.loc[order.values == 1]
    return expr.sort_index()


def counts_to_logcpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Library-size normalise raw counts and return log2(CPM + 1)."""
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def to_log2p1(values: pd.DataFrame) -> pd.DataFrame:
    """log2(x + 1) for linear TPM/FPKM input."""
    return np.log2(values.clip(lower=0) + 1.0)


def quantile_normalise_samples(expr: pd.DataFrame) -> pd.DataFrame:
    """Quantile-normalise columns (samples) onto a common distribution."""
    ranks = expr.rank(axis=0, method="average")
    reference = np.sort(expr.to_numpy(), axis=0).mean(axis=1)
    positions = (ranks - 1).to_numpy()
    lower = np.floor(positions).astype(int)
    upper = np.ceil(positions).astype(int)
    frac = positions - lower
    n = len(reference)
    lower = np.clip(lower, 0, n - 1)
    upper = np.clip(upper, 0, n - 1)
    out = reference[lower] * (1 - frac) + reference[upper] * frac
    return pd.DataFrame(out, index=expr.index, columns=expr.columns)


def sample_percentile_ranks(expr: pd.DataFrame) -> pd.DataFrame:
    """Within-sample percentile rank of every gene (platform-independent scale)."""
    return expr.rank(axis=0, pct=True, method="average")


# --------------------------------------------------------------------------
# Signature scoring
# --------------------------------------------------------------------------


def load_signatures() -> dict:
    raw = json.loads((SCRIPTS_DIR / "signatures.json").read_text())
    return {k: v for k, v in raw.items() if not k.startswith("_")}


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def signature_score(expr: pd.DataFrame, genes) -> pd.Series:
    """Mean within-cohort z-score across the genes of a signature.

    Genes absent from the matrix are dropped; the number retained is reported by
    :func:`signature_coverage`.
    """
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns)
    return zscore_rows(expr.loc[present]).mean(axis=0)


def signature_coverage(expr: pd.DataFrame, genes) -> tuple[int, int]:
    present = [g for g in genes if g in expr.index]
    return len(present), len(genes)


def ssgsea(expr: pd.DataFrame, genes, alpha: float = 0.25) -> pd.Series:
    """Single-sample GSEA enrichment score (Barbie et al., Nature 2009).

    Implemented directly so the pipeline has no R dependency.
    """
    present = [g for g in genes if g in expr.index]
    if len(present) < 3:
        return pd.Series(np.nan, index=expr.columns)

    ranks = expr.rank(axis=0, method="average")
    n_genes = ranks.shape[0]
    scores = {}
    gene_set = set(present)
    for sample in expr.columns:
        order = ranks[sample].sort_values(ascending=False)
        in_set = order.index.isin(gene_set)
        weights = np.where(in_set, order.to_numpy() ** alpha, 0.0)
        total = weights.sum()
        if total == 0:
            scores[sample] = np.nan
            continue
        cdf_in = np.cumsum(weights) / total
        cdf_out = np.cumsum(~in_set) / max(n_genes - in_set.sum(), 1)
        scores[sample] = float(np.sum(cdf_in - cdf_out))
    return pd.Series(scores)


# --------------------------------------------------------------------------
# Small-sample statistics
# --------------------------------------------------------------------------


def mann_whitney(x, y):
    """Two-sided Mann-Whitney U with rank-biserial effect size and AUC."""
    from scipy import stats

    x = np.asarray(pd.Series(x).dropna(), dtype=float)
    y = np.asarray(pd.Series(y).dropna(), dtype=float)
    if len(x) < 2 or len(y) < 2:
        return dict(n_x=len(x), n_y=len(y), u=np.nan, p=np.nan, auc=np.nan, rbc=np.nan)
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    auc = u / (len(x) * len(y))
    return dict(n_x=len(x), n_y=len(y), u=float(u), p=float(p),
                auc=float(auc), rbc=float(2 * auc - 1))


def hedges_g(x, y):
    x = np.asarray(pd.Series(x).dropna(), dtype=float)
    y = np.asarray(pd.Series(y).dropna(), dtype=float)
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan, np.nan
    sp2 = ((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2)
    if sp2 <= 0:
        return np.nan, np.nan
    d = (x.mean() - y.mean()) / np.sqrt(sp2)
    j = 1 - 3 / (4 * (nx + ny) - 9)
    g = j * d
    se = np.sqrt((nx + ny) / (nx * ny) + g**2 / (2 * (nx + ny - 2)))
    return float(g), float(se)


def permutation_pvalue(x, y, statistic, n_perm: int = 20000, seed: int = 0):
    """Exact-style two-sided permutation p-value for a two-group statistic."""
    rng = np.random.default_rng(seed)
    x = np.asarray(pd.Series(x).dropna(), dtype=float)
    y = np.asarray(pd.Series(y).dropna(), dtype=float)
    if len(x) < 2 or len(y) < 2:
        return np.nan
    observed = statistic(x, y)
    pooled = np.concatenate([x, y])
    nx = len(x)
    count = 0
    for _ in range(n_perm):
        rng.shuffle(pooled)
        if abs(statistic(pooled[:nx], pooled[nx:])) >= abs(observed) - 1e-12:
            count += 1
    return (count + 1) / (n_perm + 1)


def spearman(x, y):
    from scipy import stats

    frame = pd.DataFrame({"x": pd.Series(x, dtype=float), "y": pd.Series(y, dtype=float)}).dropna()
    if len(frame) < 4:
        return dict(n=len(frame), rho=np.nan, p=np.nan)
    rho, p = stats.spearmanr(frame["x"], frame["y"])
    return dict(n=int(len(frame)), rho=float(rho), p=float(p))


def partial_spearman(x, y, covariates):
    """Spearman correlation of x and y after linear removal of covariates on ranks."""
    from scipy import stats

    frame = pd.DataFrame({"x": pd.Series(x, dtype=float), "y": pd.Series(y, dtype=float)})
    covariates = pd.DataFrame(covariates)
    frame = frame.join(covariates, how="inner").dropna()
    if len(frame) < 6:
        return dict(n=len(frame), rho=np.nan, p=np.nan)

    ranked = frame.rank()
    design = np.column_stack([np.ones(len(ranked)), ranked[covariates.columns].to_numpy()])

    def residual(col):
        beta, *_ = np.linalg.lstsq(design, ranked[col].to_numpy(), rcond=None)
        return ranked[col].to_numpy() - design @ beta

    rx, ry = residual("x"), residual("y")
    if np.std(rx) == 0 or np.std(ry) == 0:
        return dict(n=len(frame), rho=np.nan, p=np.nan)
    rho, _ = stats.pearsonr(rx, ry)
    dof = len(frame) - 2 - covariates.shape[1]
    if dof <= 0:
        return dict(n=len(frame), rho=float(rho), p=np.nan)
    t = rho * np.sqrt(dof / max(1 - rho**2, 1e-12))
    p = 2 * stats.t.sf(abs(t), dof)
    return dict(n=int(len(frame)), rho=float(rho), p=float(p))


def benjamini_hochberg(pvalues) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = ~np.isnan(p)
    if mask.sum() == 0:
        return out
    vals = p[mask]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(vals)
    adjusted = ranked * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(n)
    result[order] = np.clip(adjusted, 0, 1)
    out[mask] = result
    return out


def random_effects_meta(effects, variances):
    """DerSimonian-Laird random-effects pooling."""
    effects = np.asarray(effects, dtype=float)
    variances = np.asarray(variances, dtype=float)
    ok = ~(np.isnan(effects) | np.isnan(variances) | (variances <= 0))
    effects, variances = effects[ok], variances[ok]
    k = len(effects)
    if k == 0:
        return dict(k=0, estimate=np.nan, se=np.nan, ci_low=np.nan, ci_high=np.nan,
                    p=np.nan, tau2=np.nan, q=np.nan, i2=np.nan)
    if k == 1:
        from scipy import stats

        se = float(np.sqrt(variances[0]))
        est = float(effects[0])
        return dict(k=1, estimate=est, se=se, ci_low=est - 1.96 * se,
                    ci_high=est + 1.96 * se, p=float(2 * stats.norm.sf(abs(est / se))),
                    tau2=0.0, q=0.0, i2=0.0)

    from scipy import stats

    w = 1 / variances
    fixed = np.sum(w * effects) / np.sum(w)
    q = float(np.sum(w * (effects - fixed) ** 2))
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    w_star = 1 / (variances + tau2)
    est = float(np.sum(w_star * effects) / np.sum(w_star))
    se = float(np.sqrt(1 / np.sum(w_star)))
    i2 = float(max(0.0, (q - (k - 1)) / q) * 100) if q > 0 else 0.0
    return dict(k=int(k), estimate=est, se=se, ci_low=est - 1.96 * se,
                ci_high=est + 1.96 * se, p=float(2 * stats.norm.sf(abs(est / se))),
                tau2=float(tau2), q=q, i2=i2)


def write_table(frame: pd.DataFrame, path: Path, **kw) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=kw.pop("index", False), **kw)
    print(f"  wrote {path.relative_to(path.parents[2])}  ({len(frame)} rows)")
