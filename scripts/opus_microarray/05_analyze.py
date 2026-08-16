#!/usr/bin/env python3
"""Step 5 - TACSTD2 / CLDN4 (and control genes) vs ICI response.

For each curated series with a usable response label we:

  * collapse probes → gene (median of log2 intensities)
  * code two binary endpoints
        OR  = CR/PR (or R / response) vs PD/NR
        DCB = CR/PR/SD vs PD   (when RECIST is available)
  * restrict to a pre-specified subset (lung-only, pre-treatment, ...)
  * test each gene with two-sided Mann–Whitney U, report AUC, cliffs-delta,
    and a 2000-shuffle permutation p-value
  * fit a univariate logistic regression (Firth-like fallback: skip if
    either class has n<3 or complete separation)
  * write boxplots

A series whose platform lacks TACSTD2/CLDN4 still runs the control-gene
analysis so the pipeline is shown to work on the same samples.

Outputs
  results/opus_microarray/association_results.tsv
  results/opus_microarray/sample_level.tsv
  results/opus_microarray/figures/*.png
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from catalog import CATALOG  # noqa: E402
from common import CONTROL_GENES, RESULTS_DIR, TARGET_GENES, ensure_dirs  # noqa: E402
from geo_io import maybe_log2  # noqa: E402

RNG = np.random.default_rng(20260816)
N_PERM = 2000


def collapse_to_genes(expr: pd.DataFrame, pmap: pd.DataFrame) -> pd.DataFrame:
    """Median of probes per gene. Index = gene, columns = samples."""
    logged, _ = maybe_log2(expr)
    pieces = {}
    if pmap is not None and len(pmap):
        for gene, sub in pmap.groupby("gene"):
            pids = [p for p in sub.probe_id.astype(str) if p in logged.index]
            if not pids:
                continue
            pieces[gene] = logged.loc[pids].median(axis=0)
    # NanoString: probe ID == gene symbol
    for g in TARGET_GENES + CONTROL_GENES:
        if g in logged.index and g not in pieces:
            pieces[g] = logged.loc[g]
    if not pieces:
        return pd.DataFrame(columns=logged.columns)
    return pd.DataFrame(pieces).T


def code_response(pheno: pd.DataFrame, rec: dict, which: str) -> pd.Series:
    labels = rec.get(f"{which}_labels") or {}
    fields = rec.get("response_fields") or []
    out = pd.Series(index=pheno.index, dtype="float")
    for fld in fields:
        if fld not in pheno.columns:
            # try a prefix match (GEO keys can be truncated / numbered)
            matches = [c for c in pheno.columns if c.startswith(fld) or fld in c]
            if not matches:
                continue
            fld = matches[0]
        mapped = pheno[fld].astype(str).str.strip().map(lambda v: labels.get(v, np.nan))
        out = out.where(~mapped.notna(), mapped)
    return out


def apply_filters(pheno: pd.DataFrame, rec: dict) -> pd.Series:
    keep = pd.Series(True, index=pheno.index)
    if rec.get("subset_lung"):
        fld = rec.get("lung_field", "source_name")
        pat = rec.get("lung_pattern", r"LUNG")
        if fld in pheno.columns:
            keep &= pheno[fld].astype(str).str.contains(pat, case=False, regex=True)
    for col, pat in (rec.get("require") or {}).items():
        if col in pheno.columns:
            keep &= pheno[col].astype(str).str.contains(pat, regex=True)
    return keep


def mannwhitney(x: np.ndarray, y: np.ndarray) -> dict:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    n1, n0 = int(x.size), int(y.size)
    if n1 < 2 or n0 < 2:
        return {"n_pos": n1, "n_neg": n0, "u": np.nan, "p_mw": np.nan, "auc": np.nan, "delta": np.nan}
    res = stats.mannwhitneyu(x, y, alternative="two-sided", method="auto")
    u = float(res.statistic)
    auc = u / (n1 * n0)
    # Cliff's delta = 2*AUC - 1
    delta = 2 * auc - 1
    return {"n_pos": n1, "n_neg": n0, "u": u, "p_mw": float(res.pvalue), "auc": auc, "delta": delta}


def perm_p(x: np.ndarray, y: np.ndarray, obs_u: float, n: int = N_PERM) -> float:
    if not np.isfinite(obs_u):
        return float("nan")
    pool = np.concatenate([x, y])
    n1 = x.size
    more = 0
    for _ in range(n):
        RNG.shuffle(pool)
        u = stats.mannwhitneyu(pool[:n1], pool[n1:], alternative="two-sided", method="asymptotic").statistic
        if abs(u - (n1 * (pool.size - n1) / 2)) >= abs(obs_u - (n1 * (pool.size - n1) / 2)):
            more += 1
    return (more + 1) / (n + 1)


def cliffs_and_logit(x: np.ndarray, y: np.ndarray) -> dict:
    """Univariate logistic P(response=1 | gene). Returns empty dict on failure."""
    try:
        import statsmodels.api as sm
    except Exception:  # noqa: BLE001
        return {}
    yy = np.concatenate([np.ones(x.size), np.zeros(y.size)])
    xx = np.concatenate([x, y])
    mask = np.isfinite(xx) & np.isfinite(yy)
    xx, yy = xx[mask], yy[mask]
    if yy.sum() < 3 or (1 - yy).sum() < 3:
        return {}
    try:
        model = sm.Logit(yy, sm.add_constant(xx)).fit(disp=False, maxiter=100)
        return {
            "logit_coef": float(model.params[1]),
            "logit_p": float(model.pvalues[1]),
            "logit_or_per_unit": float(math.exp(model.params[1])),
        }
    except Exception:  # noqa: BLE001
        return {}


def boxplot(df: pd.DataFrame, gene: str, endpoint: str, title: str, dest: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    groups = []
    labels = []
    for lab, name in ((0, "non-resp"), (1, "resp")):
        v = df.loc[df[endpoint] == lab, gene].dropna().to_numpy()
        groups.append(v)
        labels.append(f"{name}\nn={v.size}")
    ax.boxplot(groups, tick_labels=labels, widths=0.55, patch_artist=True,
               boxprops=dict(facecolor="#dbeafe", edgecolor="#1e3a5f"),
               medianprops=dict(color="#b45309", linewidth=2))
    for i, v in enumerate(groups, start=1):
        jitter = RNG.normal(0, 0.04, size=v.size)
        ax.scatter(np.full(v.size, i) + jitter, v, s=18, c="#1e3a5f", alpha=0.7, zorder=3)
    ax.set_ylabel(f"{gene} (log2 intensity)")
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def main() -> int:
    ensure_dirs()
    ext = RESULTS_DIR / "extracted"
    fig_dir = RESULTS_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    assoc_rows = []
    sample_rows = []

    for rec in CATALOG:
        gse = rec["gse"]
        if rec["decision"] not in {"analyze", "document_absent"}:
            continue
        pheno_p = ext / f"{gse}_pheno.tsv"
        expr_p = ext / f"{gse}_expr.tsv"
        pmap_p = ext / f"{gse}_probe_map.tsv"
        if not pheno_p.exists() or not expr_p.exists():
            print(f"[{gse}] missing extracted tables, skip")
            continue
        pheno = pd.read_csv(pheno_p, sep="\t", index_col=0)
        expr = pd.read_csv(expr_p, sep="\t", index_col=0)
        pmap = pd.read_csv(pmap_p, sep="\t") if pmap_p.exists() and pmap_p.stat().st_size > 1 else pd.DataFrame()
        genes = collapse_to_genes(expr, pmap)
        keep = apply_filters(pheno, rec)
        pheno_f = pheno.loc[keep]
        # align expression columns to pheno index (GSM)
        common = [c for c in pheno_f.index if c in genes.columns]
        if not common:
            # some matrices use title as column; try geo_accession column
            if "geo_accession" in pheno_f.columns:
                remap = dict(zip(pheno_f["geo_accession"], pheno_f.index))
                genes = genes.rename(columns=remap)
                common = [c for c in pheno_f.index if c in genes.columns]
        pheno_f = pheno_f.loc[common]
        genes = genes[common]

        or_y = code_response(pheno_f, rec, "or")
        dcb_y = code_response(pheno_f, rec, "dcb")
        print(f"[{gse}] n_filt={len(pheno_f)}  OR pos/neg={int((or_y==1).sum())}/{(or_y==0).sum()}  "
              f"genes={list(genes.index)}")

        for endpoint, y in (("OR", or_y), ("DCB", dcb_y)):
            if y.notna().sum() < 4 or y.nunique(dropna=True) < 2:
                continue
            # skip duplicate DCB when it equals OR (binary R/NR cohorts)
            if endpoint == "DCB" and y.equals(or_y):
                continue
            for gene in list(genes.index):
                vec = genes.loc[gene]
                x = vec[y == 1].to_numpy(dtype=float)
                z = vec[y == 0].to_numpy(dtype=float)
                mw = mannwhitney(x, z)
                pp = perm_p(x[np.isfinite(x)], z[np.isfinite(z)], mw["u"])
                lg = cliffs_and_logit(x[np.isfinite(x)], z[np.isfinite(z)])
                mean_pos = float(np.nanmean(x)) if x.size else np.nan
                mean_neg = float(np.nanmean(z)) if z.size else np.nan
                # Epithelial genes (TACSTD2/CLDN4) sit near the floor in sorted
                # lymphocytes / whole blood; flag so we do not over-interpret.
                floor = int(
                    gene in TARGET_GENES
                    and np.isfinite(mean_pos)
                    and np.isfinite(mean_neg)
                    and max(mean_pos, mean_neg) < 4.5
                )
                row = {
                    "gse": gse,
                    "tissue": rec["tissue"],
                    "decision": rec["decision"],
                    "endpoint": endpoint,
                    "gene": gene,
                    "is_target": int(gene in TARGET_GENES),
                    "n_filt": int(len(pheno_f)),
                    "mean_pos": mean_pos,
                    "mean_neg": mean_neg,
                    "delta_mean": (mean_pos - mean_neg) if np.isfinite(mean_pos) and np.isfinite(mean_neg) else np.nan,
                    **mw,
                    "p_perm": pp,
                    **lg,
                    "near_detection_floor": floor,
                    "reason": rec["reason"][:160],
                }
                assoc_rows.append(row)
                if gene in TARGET_GENES or gene in {"CD274", "CD8A"}:
                    boxplot(
                        pd.DataFrame({gene: vec, endpoint: y}),
                        gene,
                        endpoint,
                        f"{gse} {gene} vs {endpoint}",
                        fig_dir / f"{gse}_{gene}_{endpoint}.png",
                    )

        # sample-level dump for verification
        sl = pheno_f.copy()
        sl["gse"] = gse
        for gene in genes.index:
            sl[f"expr_{gene}"] = genes.loc[gene]
        sl["OR"] = or_y
        sl["DCB"] = dcb_y
        sample_rows.append(sl.reset_index().rename(columns={"index": "sample_id"}))

    assoc = pd.DataFrame(assoc_rows)
    assoc.to_csv(RESULTS_DIR / "association_results.tsv", sep="\t", index=False)
    if sample_rows:
        pd.concat(sample_rows, ignore_index=True, sort=False).to_csv(
            RESULTS_DIR / "sample_level.tsv", sep="\t", index=False
        )
    print(f"[write] association_results.tsv ({len(assoc)} tests)")
    if len(assoc):
        show = assoc[assoc.is_target == 1][["gse", "endpoint", "gene", "n_pos", "n_neg", "auc", "p_mw", "p_perm"]]
        print(show.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
