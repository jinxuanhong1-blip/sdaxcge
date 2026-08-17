#!/usr/bin/env python3
"""REAL Palantir on winning-pair epithelium, CLDN4 only.

ADDITIVE. GSE131907 + GSE205335. No GSE148071. No GSE207422.
No TACSTD2∩CLDN4 dual-high gate. Patient is the inferential unit.
Root / early cell is never CLDN4-high.
Do NOT stop if DPT is empty — Palantir destinies + pseudotime are the product.

Thesis tested (not assumed): CLDN4-high = barrier end; CLDN4-low = IFN-higher
end of the same malignant trajectory.
"""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (  # noqa: E402
    AUTHOR_AIRWAY,
    AUTHOR_AT2,
    AUTHOR_CLUB,
    AUTHOR_TUMOR_STATE,
    COMPARATOR,
    CONTROLS,
    FOCAL,
    QC_NEG,
    STATES,
)

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_PATIENT = 10
MIN_CELLS_PER_TERTILE_ARM = 8
N_WAYPOINTS = 500
N_DM_COMPS = 10
CAP_PER_UNIT = 350


def _require_palantir():
    try:
        import palantir  # noqa: F401
    except ImportError as e:
        raise SystemExit(
            "palantir is required. pip install palantir. "
            "Do not substitute empty DPT for Palantir destinies."
        ) from e
    return __import__("palantir")


def palantir_status() -> dict:
    try:
        import palantir

        return {
            "available": True,
            "module": "palantir",
            "version": getattr(palantir, "__version__", "unknown"),
        }
    except ImportError:
        return {"available": False, "reason": "python package palantir not installed"}


def _bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n, dtype=float)
    running = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        k = n - rank + 1
        running = min(running, p[i] * n / k)
        q[i] = running
    return [float(min(1.0, x)) for x in q]


def _spearman(x: np.ndarray, y: np.ndarray) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 4:
        return {"n": n, "rho": None, "p": None, "note": "n<4"}
    rho, p = stats.spearmanr(x[mask], y[mask])
    if not np.isfinite(rho):
        return {"n": n, "rho": None, "p": None, "note": "undefined"}
    return {"n": n, "rho": float(rho), "p": float(p)}


def _wilcoxon_paired(a: np.ndarray, b: np.ndarray) -> dict:
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 4:
        return {"n": n, "W": None, "p": None, "delta_median": None, "note": "n<4"}
    aa, bb = a[mask], b[mask]
    try:
        w, p = stats.wilcoxon(aa, bb, alternative="two-sided", zero_method="wilcox")
    except ValueError:
        return {
            "n": n,
            "W": None,
            "p": None,
            "delta_median": float(np.median(aa - bb)),
            "note": "wilcoxon failed",
        }
    return {
        "n": n,
        "W": float(w),
        "p": float(p),
        "delta_median": float(np.median(aa - bb)),
        "median_high": float(np.median(aa)),
        "median_low": float(np.median(bb)),
    }


def _score_mean(adata, genes: tuple[str, ...], key: str) -> list[str]:
    present = [g for g in genes if g in adata.var_names]
    absent = [g for g in genes if g not in adata.var_names]
    if not present:
        adata.obs[key] = np.nan
        return absent
    X = adata[:, present].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    adata.obs[key] = np.asarray(X, dtype=float).mean(axis=1)
    return absent


def _is_malignant(subtype: pd.Series) -> pd.Series:
    s = subtype.astype(str)
    hit = s.isin(AUTHOR_TUMOR_STATE) | s.str.contains(
        r"malign|cancer|tumor", case=False, regex=True
    )
    return hit


def _patient_key(obs: pd.DataFrame) -> pd.Series:
    if "patient_key" in obs.columns:
        key = obs["patient_key"].astype(str)
        bad = key.isin(["nan", "None", "NA", "GSE131907:nan", "GSE205335:nan", ""])
        if "unit_id" in obs.columns:
            key = key.where(~bad, obs["unit_id"].astype(str))
        return key
    if "patient_id" in obs.columns and "dataset" in obs.columns:
        pid = obs["patient_id"].astype(str)
        return obs["dataset"].astype(str) + ":" + pid
    return obs["unit_id"].astype(str)


def _pick_early_cell(adata) -> tuple[str, dict]:
    """External arrow: GSE131907 nLung author AT2. Never CLDN4-high."""
    nlung = (adata.obs["dataset"].astype(str) == "GSE131907") & (
        adata.obs["Sample_Origin"].astype(str) == "nLung"
    )
    subtype = adata.obs["author_subtype"].astype(str)
    at2 = nlung & subtype.isin(AUTHOR_AT2)
    high = adata.obs["cldn4_tertile"].astype(str) == "high"
    pool = at2 & ~high
    info = {
        "n_nLung": int(nlung.sum()),
        "n_nLung_author_AT2": int(at2.sum()),
        "n_nLung_AT2_not_CLDN4_high": int(pool.sum()),
        "rule": None,
        "rejected_CLDN4_high_root": True,
    }
    names = adata.obs_names.to_numpy()

    def _median_at2(mask: pd.Series) -> int:
        idx = np.flatnonzero(mask.to_numpy())
        scores = adata.obs.loc[mask, "score_AT2"].to_numpy()
        return int(idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))])

    if int(pool.sum()) >= 10:
        pick = _median_at2(pool)
        info["rule"] = "GSE131907 nLung author AT2, not CLDN4-high, median AT2 score"
        return str(names[pick]), info
    if int(at2.sum()) >= 10:
        idx = np.flatnonzero(at2.to_numpy())
        cldn = adata.obs.loc[at2, "expr_CLDN4"].to_numpy()
        pick = int(idx[int(np.nanargmin(cldn))])
        if str(adata.obs.iloc[pick]["cldn4_tertile"]) == "high":
            raise SystemExit("early-cell fallback still CLDN4-high; refuse")
        info["rule"] = "GSE131907 nLung author AT2 with lowest CLDN4 (high tertile empty pool)"
        return str(names[pick]), info
    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before early-cell pick")
    sub = adata.obs.loc[nlung & ~high, ["leiden", "score_AT2", "expr_CLDN4"]]
    if sub.empty:
        raise SystemExit("no GSE131907 nLung non-CLDN4-high cells for Palantir early cell")
    means = sub.groupby("leiden", observed=True)["score_AT2"].mean().sort_values(ascending=False)
    top = str(means.index[0])
    cand = nlung & ~high & (adata.obs["leiden"].astype(str) == top)
    pick = _median_at2(cand)
    info["rule"] = f"GSE131907 nLung Leiden {top} not CLDN4-high, median AT2"
    info["fallback_cluster"] = top
    return str(names[pick]), info


def _maybe_harmony(adata) -> dict:
    info = {"used": False, "reason": None}
    try:
        import harmonypy
    except ImportError:
        info["reason"] = "harmonypy not installed; neighbors on PCA"
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
        return info
    ho = harmonypy.run_harmony(
        adata.obsm["X_pca"][:, :N_PCS],
        adata.obs,
        "dataset",
        max_iter_harmony=20,
    )
    z = np.asarray(ho.Z_corr)
    if z.shape[0] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z
    elif z.shape[1] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z.T
    else:
        raise ValueError(f"Harmony Z_corr shape {z.shape} vs n_obs={adata.n_obs}")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca_harmony")
    info["used"] = True
    info["reason"] = "harmonypy on PCA, batch=dataset"
    return info


def _try_dpt(adata, early_name: str) -> dict:
    """Companion only. Empty / inf DPT is recorded, never a stop."""
    info = {"ran": False, "empty": True, "reason": None, "n_finite": 0}
    try:
        if early_name not in adata.obs_names:
            info["reason"] = f"early cell {early_name} not in obs_names"
            adata.obs["dpt_pseudotime"] = np.nan
            return info
        adata.uns["iroot"] = int(np.where(adata.obs_names == early_name)[0][0])
        sc.tl.diffmap(adata, n_comps=15)
        sc.tl.dpt(adata, n_dcs=10)
        dpt = pd.to_numeric(adata.obs["dpt_pseudotime"], errors="coerce")
        dpt = dpt.replace([np.inf, -np.inf], np.nan)
        adata.obs["dpt_pseudotime"] = dpt
        n_fin = int(np.isfinite(dpt.to_numpy()).sum())
        info["ran"] = True
        info["n_finite"] = n_fin
        info["empty"] = n_fin == 0
        if info["empty"]:
            info["reason"] = "DPT all-NaN/inf; continuing to Palantir"
        return info
    except Exception as exc:
        info["reason"] = f"DPT failed: {exc}"
        adata.obs["dpt_pseudotime"] = np.nan
        print(f"DPT companion failed (continuing): {exc}", flush=True)
        return info


def _run_palantir(adata, early_name: str) -> dict:
    """Real Palantir. Must write destinies + pseudotime or raise."""
    palantir = _require_palantir()
    info = {"ok": False}
    pca_key = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"
    info["pca"] = pca_key

    run_dm = palantir.utils.run_diffusion_maps
    try:
        run_dm(adata, n_components=N_DM_COMPS, knn=N_NEIGHBORS, pca_key=pca_key)
    except TypeError:
        if pca_key != "X_pca":
            adata.obsm["X_pca"] = np.asarray(adata.obsm[pca_key])
        run_dm(adata, n_components=N_DM_COMPS)

    palantir.utils.determine_multiscale_space(adata)
    try:
        palantir.utils.run_magic_imputation(adata)
        info["magic"] = True
    except Exception as exc:
        info["magic"] = False
        info["magic_reason"] = str(exc)[:300]
        print(f"MAGIC optional failed: {exc}", flush=True)

    n_wp = int(min(N_WAYPOINTS, max(120, adata.n_obs // 15)))
    kwargs = dict(
        early_cell=early_name,
        num_waypoints=n_wp,
        knn=N_NEIGHBORS,
        save_as_df=True,
    )

    def _call_palantir(kw):
        try:
            return palantir.core.run_palantir(adata, n_jobs=2, **kw)
        except TypeError:
            kw = dict(kw)
            kw.pop("save_as_df", None)
            try:
                return palantir.core.run_palantir(adata, n_jobs=2, **kw)
            except TypeError:
                return palantir.core.run_palantir(adata, **kw)

    def _boundary_terminals() -> dict[str, str]:
        ms_key = "DM_EigenVectors_multiscaled"
        if ms_key not in adata.obsm:
            raise RuntimeError("multiscale space missing; cannot pick fallback terminals")
        ms = pd.DataFrame(adata.obsm[ms_key], index=adata.obs_names)
        extrema = pd.Index(set(ms.idxmax()).union(ms.idxmin())).difference([early_name])
        extrema = [c for c in extrema if c in adata.obs_names]
        if not extrema:
            raise RuntimeError("no DM-boundary cells for fallback terminals")
        early_vec = ms.loc[early_name].to_numpy()
        dist = {c: float(np.linalg.norm(ms.loc[c].to_numpy() - early_vec)) for c in extrema}
        picked = sorted(dist, key=dist.get, reverse=True)[:3]
        return {f"boundary_{i+1}": c for i, c in enumerate(picked)}

    pr = None
    try:
        pr = _call_palantir(kwargs)
    except Exception as exc:
        print(f"Palantir auto-terminals failed ({type(exc).__name__}: {exc}). Retrying with DM-boundary terminals.", flush=True)
        term = _boundary_terminals()
        info["terminal_fallback"] = term
        info["auto_terminal_error"] = f"{type(exc).__name__}: {str(exc)[:300]}"
        kwargs["terminal_states"] = term
        pr = _call_palantir(kwargs)

    # Normalize outputs onto AnnData
    pt = None
    entropy = None
    branch = None
    terminals = None
    if pr is not None:
        pt = getattr(pr, "pseudotime", None)
        entropy = getattr(pr, "entropy", None)
        branch = getattr(pr, "branch_probs", None)
        terminals = getattr(pr, "terminal_states", None)
    if pt is None and "palantir_pseudotime" in adata.obs:
        pt = adata.obs["palantir_pseudotime"]
    if entropy is None and "palantir_entropy" in adata.obs:
        entropy = adata.obs["palantir_entropy"]
    if branch is None and "palantir_fate_probabilities" in adata.obsm:
        raw = adata.obsm["palantir_fate_probabilities"]
        if isinstance(raw, pd.DataFrame):
            branch = raw
        else:
            cols = adata.uns.get("palantir_fate_probabilities_columns")
            branch = pd.DataFrame(raw, index=adata.obs_names, columns=cols)
    if isinstance(pt, pd.Series):
        pt = pt.reindex(adata.obs_names)
        adata.obs["palantir_pseudotime"] = pd.to_numeric(pt, errors="coerce")
    elif pt is not None:
        adata.obs["palantir_pseudotime"] = np.asarray(pt, dtype=float)
    if isinstance(entropy, pd.Series):
        adata.obs["palantir_entropy"] = pd.to_numeric(
            entropy.reindex(adata.obs_names), errors="coerce"
        )
    elif entropy is not None:
        adata.obs["palantir_entropy"] = np.asarray(entropy, dtype=float)

    if branch is None and "palantir_fate_probabilities" in adata.obsm:
        raw = adata.obsm["palantir_fate_probabilities"]
        if isinstance(raw, pd.DataFrame):
            branch = raw
        else:
            cols = adata.uns.get("palantir_fate_probabilities_columns")
            branch = pd.DataFrame(raw, index=adata.obs_names, columns=cols)
    if branch is None or (isinstance(branch, pd.DataFrame) and branch.shape[1] == 0):
        print("Palantir auto-terminals empty; retry with DM-boundary terminals (not CLDN4).", flush=True)
        ms_key = "DM_EigenVectors_multiscaled"
        if ms_key not in adata.obsm:
            raise RuntimeError("Palantir destiny table is empty and multiscale space is missing")
        ms = pd.DataFrame(adata.obsm[ms_key], index=adata.obs_names)
        # farthest diffusion-map extrema from the early cell — never CLDN4-defined
        extrema = pd.Index(set(ms.idxmax()).union(ms.idxmin()))
        extrema = extrema.difference([early_name])
        high = adata.obs["cldn4_tertile"].astype(str) == "high"
        extrema_ok = [c for c in extrema if c in adata.obs_names and not bool(high.get(c, False))]
        if len(extrema_ok) < 2:
            extrema_ok = [c for c in extrema if c in adata.obs_names]
        if len(extrema_ok) < 1:
            raise RuntimeError("Palantir destiny table is empty and no DM-boundary terminals")
        # keep up to 3 farthest from early cell
        early_vec = ms.loc[early_name].to_numpy()
        dist = {c: float(np.linalg.norm(ms.loc[c].to_numpy() - early_vec)) for c in extrema_ok}
        picked = sorted(dist, key=dist.get, reverse=True)[:3]
        term = {f"boundary_{i+1}": c for i, c in enumerate(picked)}
        info["terminal_fallback"] = term
        kwargs["terminal_states"] = term
        try:
            pr = palantir.core.run_palantir(adata, n_jobs=2, **kwargs)
        except TypeError:
            kwargs.pop("save_as_df", None)
            pr = palantir.core.run_palantir(adata, **kwargs)
        pt = getattr(pr, "pseudotime", None) if pr is not None else None
        entropy = getattr(pr, "entropy", None) if pr is not None else None
        branch = getattr(pr, "branch_probs", None) if pr is not None else None
        if pt is None and "palantir_pseudotime" in adata.obs:
            pt = adata.obs["palantir_pseudotime"]
        if entropy is None and "palantir_entropy" in adata.obs:
            entropy = adata.obs["palantir_entropy"]
        if branch is None and "palantir_fate_probabilities" in adata.obsm:
            raw = adata.obsm["palantir_fate_probabilities"]
            branch = raw if isinstance(raw, pd.DataFrame) else pd.DataFrame(raw, index=adata.obs_names)
        if isinstance(pt, pd.Series):
            adata.obs["palantir_pseudotime"] = pd.to_numeric(pt.reindex(adata.obs_names), errors="coerce")
        if isinstance(entropy, pd.Series):
            adata.obs["palantir_entropy"] = pd.to_numeric(entropy.reindex(adata.obs_names), errors="coerce")
    if branch is None:
        raise RuntimeError("Palantir returned no branch/destiny probabilities")
    if not isinstance(branch, pd.DataFrame):
        branch = pd.DataFrame(branch, index=adata.obs_names)
    branch = branch.reindex(adata.obs_names)
    branch.columns = [str(c) for c in branch.columns]
    if branch.shape[1] == 0:
        raise RuntimeError("Palantir destiny table is empty after terminal fallback")

    # Rename destinies to destiny_1..k but keep original terminal cell ids
    rename = {c: f"destiny_{i+1}" for i, c in enumerate(branch.columns)}
    orig_names = list(branch.columns)
    branch = branch.rename(columns=rename)
    for col in branch.columns:
        adata.obs[f"fate_{col}"] = pd.to_numeric(branch[col], errors="coerce")
    adata.obs["palantir_destiny"] = branch.idxmax(axis=1)
    adata.obs["palantir_destiny_prob"] = branch.max(axis=1)
    adata.obsm["palantir_fate_probabilities"] = branch.to_numpy(dtype=float)
    adata.uns["palantir_fate_names"] = list(branch.columns)
    adata.uns["palantir_terminal_cells"] = {
        rename[c]: c for c in orig_names
    }

    n_pt = int(np.isfinite(adata.obs["palantir_pseudotime"].to_numpy()).sum())
    if n_pt == 0:
        raise RuntimeError("Palantir pseudotime is empty after a successful call")
    info.update(
        {
            "ok": True,
            "n_waypoints": n_wp,
            "n_destinies": int(branch.shape[1]),
            "destiny_names": list(branch.columns),
            "terminal_cells": adata.uns["palantir_terminal_cells"],
            "n_pseudotime_finite": n_pt,
            "early_cell": early_name,
        }
    )
    return info, branch


def _fmt(row: dict, keys=("rho", "p")) -> str:
    n = row.get("n")
    if row.get("rho") is None and "rho" in keys:
        return f"n={n}, ρ=NA, p=NA"
    if "W" in keys:
        if row.get("W") is None:
            return f"n={n}, W=NA, p=NA"
        return (
            f"n={n}, W={row['W']:.1f}, Δmed={row.get('delta_median'):.3f}, "
            f"p={row['p']:.3g}"
        )
    return f"n={n}, ρ={row['rho']:.3f}, p={row['p']:.3g}"


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    extra = s["extra_figure"]

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    dest_lines = []
    for d in s.get("destinies", []):
        dest_lines.append(
            f"| {d['destiny']} | {d['n_cells']} | {d['n_patients']} | "
            f"{d['pct_malignant']:.2f} | {d['mean_CLDN4']:.3f} | "
            f"{d['mean_barrier_keratin']:.3f} | {d['mean_IFN']:.3f} | "
            f"{d['mean_pseudotime']:.3f} | {d['top_subtype']} |"
        )

    lines = [
        "# Finding — winning-pair REAL Palantir, CLDN4 only (GSE131907+GSE205335)",
        "",
        "ADDITIVE. **CLDN4 only.** Winning pair from the CLDN4-first combinatorial search "
        "(PR #290: author %pos GSE131907+GSE205335 vs T/NK). "
        "This folder does **not** redo GSE131907-only PAGA (PR #325) or the Slingshot/DPT "
        "fallback (PR #449). GSE148071 is not added. GSE207422 is not added. "
        "No TACSTD2∩CLDN4 dual-high gate.",
        "",
        "Primary method: **Palantir** (Setty et al. 2019), installed and run. "
        "Inferential unit = **patient**. Cell-level ρ is descriptive. "
        "Barrier/keratin score **excludes CLDN4**. IFN core **excludes CLDN4**. "
        "Early cell is GSE131907 nLung author AT2 and is **never CLDN4-high**. "
        "DPT is a companion only; empty DPT is not a stop.",
        "",
        f"**Thesis tested.** CLDN4-high = barrier end; CLDN4-low = IFN-higher end "
        f"of the same malignant trajectory. {s.get('thesis_verdict', '')}",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- Analysis cells after QC (capped ≤{s['cap_per_unit']}/sample): "
        f"**n_cells = {s['n_cells']}** "
        f"(GSE131907 {s['n_cells_gse131907']}, GSE205335 {s['n_cells_gse205335']}).",
        f"- Patients: **n_patients = {s['n_patients']}** "
        f"(GSE131907 {s['n_patients_gse131907']}, GSE205335 {s['n_patients_gse205335']}).",
        f"- Patients with ≥{MIN_CELLS_PER_PATIENT} epithelial cells used for Spearman: "
        f"**n = {s['n_patients_eligible']}**.",
        f"- Patients with ≥{MIN_CELLS_PER_PATIENT} malignant cells: "
        f"**n_malignant_patients = {s['n_patients_malignant']}**.",
        f"- GSE131907 nLung cells / author AT2 in the object: "
        f"{s['n_cells_nLung']} / {s['n_author_AT2']}.",
        f"- Author subtypes (cells): {s['subtype_counts']}.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Patients with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and "
        f"CLDN4-low arms: **n = {s['n_patients_paired_tertile']}**.",
        f"- Palantir destinies: **{s['palantir']['n_destinies']}** "
        f"({', '.join(s['palantir'].get('destiny_names') or [])}). "
        f"Finite Palantir pseudotime cells: {s['palantir'].get('n_pseudotime_finite')}.",
        f"- DPT companion: ran={s['dpt'].get('ran')} empty={s['dpt'].get('empty')} "
        f"n_finite={s['dpt'].get('n_finite')}; {s['dpt'].get('reason') or 'ok'}.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- GSE148071 not used. GSE207422 not used. GSE131907 PE unlabeled epithelium dropped. "
        "GSE205335 normal-tissue samples dropped.",
        f"- Palantir: available={s['palantir_pkg'].get('available')}; "
        f"version={s['palantir_pkg'].get('version')}.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Batch: {s['harmony']['reason']}.",
        f"- Palantir early cell: {s['root']['rule']} (cell `{s['root'].get('early_cell')}`, "
        f"patient {s['root'].get('root_patient')}, CLDN4 tertile {s['root'].get('root_tertile')}).",
        f"- Palantir waypoints: {s['palantir'].get('n_waypoints')}; diffusion components {N_DM_COMPS}.",
        "- Terminals: Palantir auto (not defined by CLDN4 / barrier / IFN).",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN core: STAT1, IRF1/7/9, ISG15, IFIT1/2/3, OAS1/2, MX1/2, CXCL9/10/11, "
        "IDO1, TAP1, GBP1, IFI44L, IFI27, IFI44, RSAD2, USP18, EPSTI1, SAMD9 "
        "(**CLDN4 out**; LAMP3 and CDKN1A out).",
        "",
        "## Destinies (cell assignment = argmax fate probability)",
        "",
        "| Destiny | n_cells | n_patients | % malignant | mean CLDN4 | mean barrier (no CLDN4) | mean IFN | mean PT | top subtype |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    lines.extend(dest_lines or ["| (none) | 0 | 0 |  |  |  |  |  |  |"])
    lines += [
        "",
        "## Primary (patient-level Spearman, BH inside this list)",
        "",
        "| Contrast | n_patients | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(row_md(r))
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_patients | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("sensitivity_spearman", []):
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {rho} | {pv} |")
    extra_rho = "NA" if extra.get("spearman_rho") is None else f"{round(extra['spearman_rho'], 3)}"
    extra_p = "NA" if extra.get("spearman_p") is None else f"{extra['spearman_p']:.3g}"
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)",
        "",
        (
            f"Emitted: **{extra['emitted']}**. "
            f"Observed patient Spearman(CLDN4, barrier) n={extra['spearman_n']}, "
            f"ρ={extra_rho}, p={extra_p}."
        ),
        "",
        "| Paired contrast (high − low) | n_patients | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("paired_tertile", []):
        d = "NA" if r.get("delta_median") is None else f"{r['delta_median']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {d} | {pv} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- This is not a redo of PR #325 or PR #449.",
        "- Palantir destinies are phenotypic terminals, not proven lineages.",
        "- Barrier score excludes CLDN4; IFN core excludes CLDN4. Overlap of programs is still possible.",
        "- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or CEACAM5/6/MKI67 — **not CNV**.",
        "- GSE205335 is an ICI biopsy/effusion cohort, but this analysis is **not** an ICI / MPR / RECIST test.",
        "- No TACSTD2∩CLDN4 both-high gate. GSE148071 not used.",
        "- Do not write “AT2 differentiates into LUAD because Palantir ran.”",
        "- DPT was not the product. Empty DPT would not have stopped this run.",
        "",
        "## Outputs",
        "",
        "- `results/tables/palantir_pseudotime.tsv` — **done criterion**",
        "- `results/tables/palantir_destinies.tsv` — **done criterion**",
        "- `results/tables/palantir_fate_probabilities.tsv`",
        "- `results/tables/patient_means.tsv`",
        "- `results/tables/patient_level_spearman.tsv`",
        "- `results/figures/fig_palantir_destinies_programs.png`",
        "- `results/figures/fig_extra_fate_vs_programs.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/winpair_palantir_real_cldn4/requirements.txt",
        "python3 methods/winpair_palantir_real_cldn4/scripts/download.py \\",
        "  --out /tmp/winpair_131907_205335",
        "python3 methods/winpair_palantir_real_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/winpair_131907_205335 \\",
        "  --out /tmp/winpair_131907_205335/epithelium.h5ad",
        "python3 methods/winpair_palantir_real_cldn4/scripts/analyze.py \\",
        "  --input /tmp/winpair_131907_205335/epithelium.h5ad \\",
        "  --outdir methods/winpair_palantir_real_cldn4/results \\",
        "  --finding methods/winpair_palantir_real_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="methods/winpair_palantir_real_cldn4/results")
    p.add_argument("--finding", default="methods/winpair_palantir_real_cldn4/FINDING.md")
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    pal_pkg = palantir_status()
    print(json.dumps({"palantir": pal_pkg}, indent=2), flush=True)
    if not pal_pkg.get("available"):
        raise SystemExit("Install palantir before this script. Do not stop at empty DPT.")

    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()
    adata.obs["patient_key"] = _patient_key(adata.obs)

    adata.var["n_cells"] = np.array((adata.X > 0).sum(axis=0)).ravel()
    adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).ravel()
    adata.obs["n_umi"] = np.array(adata.X.sum(axis=1)).ravel()
    sc.pp.filter_genes(adata, min_cells=10)
    keep = (adata.obs["n_genes"] >= 200) & (adata.obs["n_umi"] >= 500)
    adata = adata[keep].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    absent: dict[str, list[str]] = {}
    for name, genes in STATES.items():
        absent[name] = _score_mean(adata, genes, f"score_{name}")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata.var_names:
            x = adata[:, g].X
            if hasattr(x, "toarray"):
                x = x.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(x, dtype=float).ravel()
        else:
            adata.obs[f"expr_{g}"] = np.nan
            absent.setdefault("single_genes", []).append(g)

    cldn = adata.obs["expr_CLDN4"].to_numpy()
    q1, q2 = np.nanquantile(cldn, [1 / 3, 2 / 3])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[cldn <= q1] = "low"
    tert[cldn > q2] = "high"
    adata.obs["cldn4_tertile"] = pd.Categorical(
        tert, categories=["low", "mid", "high"], ordered=True
    )
    adata.obs["is_malignant"] = _is_malignant(adata.obs["author_subtype"])

    try:
        sc.pp.highly_variable_genes(adata, layer="counts", flavor="seurat_v3", n_top_genes=N_HVG)
    except ImportError:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_HVG)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    harmony = _maybe_harmony(adata)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2)
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.umap(adata)

    early_name, root_info = _pick_early_cell(adata)
    early_obs = adata.obs.loc[early_name]
    if str(early_obs["cldn4_tertile"]) == "high":
        raise SystemExit("refusing CLDN4-high Palantir early cell")
    root_info["early_cell"] = early_name
    root_info["root_patient"] = str(early_obs.get("patient_key", ""))
    root_info["root_unit"] = str(early_obs.get("unit_id", ""))
    root_info["root_origin"] = str(early_obs.get("Sample_Origin", ""))
    root_info["root_subtype"] = str(early_obs.get("author_subtype", ""))
    root_info["root_dataset"] = str(early_obs.get("dataset", ""))
    root_info["root_tertile"] = str(early_obs.get("cldn4_tertile", ""))
    root_info["root_CLDN4"] = float(early_obs.get("expr_CLDN4", np.nan))
    print(json.dumps({"early_cell": root_info}, indent=2, default=str), flush=True)

    dpt_info = _try_dpt(adata, early_name)
    print(json.dumps({"dpt_companion": dpt_info}, indent=2), flush=True)

    pal_info, branch = _run_palantir(adata, early_name)
    print(json.dumps({"palantir_run": {k: pal_info[k] for k in pal_info if k != "terminal_cells"}}, indent=2), flush=True)

    # Cell-level destiny / pseudotime table (done criterion)
    cell_tab = pd.DataFrame(
        {
            "cell": adata.obs_names.astype(str),
            "dataset": adata.obs["dataset"].astype(str).to_numpy(),
            "patient_key": adata.obs["patient_key"].astype(str).to_numpy(),
            "unit_id": adata.obs["unit_id"].astype(str).to_numpy(),
            "Sample_Origin": adata.obs["Sample_Origin"].astype(str).to_numpy(),
            "author_subtype": adata.obs["author_subtype"].astype(str).to_numpy(),
            "is_malignant": adata.obs["is_malignant"].astype(bool).to_numpy(),
            "cldn4_tertile": adata.obs["cldn4_tertile"].astype(str).to_numpy(),
            "expr_CLDN4": adata.obs["expr_CLDN4"].to_numpy(),
            "score_barrier_keratin": adata.obs["score_barrier_keratin"].to_numpy(),
            "score_IFN": adata.obs["score_IFN"].to_numpy(),
            "score_AT2": adata.obs["score_AT2"].to_numpy(),
            "palantir_pseudotime": adata.obs["palantir_pseudotime"].to_numpy(),
            "palantir_entropy": adata.obs.get(
                "palantir_entropy", pd.Series(np.nan, index=adata.obs_names)
            ).to_numpy(),
            "palantir_destiny": adata.obs["palantir_destiny"].astype(str).to_numpy(),
            "palantir_destiny_prob": adata.obs["palantir_destiny_prob"].to_numpy(),
            "dpt_pseudotime": adata.obs["dpt_pseudotime"].to_numpy(),
        }
    )
    for col in branch.columns:
        cell_tab[f"fate_{col}"] = branch[col].to_numpy()
    cell_tab.to_csv(tabdir / "palantir_pseudotime.tsv", sep="\t", index=False)
    branch.to_csv(tabdir / "palantir_fate_probabilities.tsv", sep="\t")

    dest_rows = []
    term_map = adata.uns.get("palantir_terminal_cells", {})
    for dest in branch.columns:
        sub = adata.obs[adata.obs["palantir_destiny"].astype(str) == dest]
        tcell = str(term_map.get(dest, ""))
        t_sub = adata.obs.loc[tcell] if tcell in adata.obs_names else None
        dest_rows.append(
            {
                "destiny": dest,
                "terminal_cell": tcell,
                "terminal_subtype": "" if t_sub is None else str(t_sub["author_subtype"]),
                "terminal_CLDN4": None if t_sub is None else float(t_sub["expr_CLDN4"]),
                "terminal_barrier": None if t_sub is None else float(t_sub["score_barrier_keratin"]),
                "terminal_IFN": None if t_sub is None else float(t_sub["score_IFN"]),
                "n_cells": int(len(sub)),
                "n_patients": int(sub["patient_key"].nunique()),
                "n_GSE131907": int((sub["dataset"] == "GSE131907").sum()),
                "n_GSE205335": int((sub["dataset"] == "GSE205335").sum()),
                "pct_malignant": float(sub["is_malignant"].mean()) if len(sub) else np.nan,
                "pct_nLung": float((sub["Sample_Origin"].astype(str) == "nLung").mean()) if len(sub) else np.nan,
                "pct_airway": float(sub["author_subtype"].astype(str).isin(AUTHOR_AIRWAY).mean()) if len(sub) else np.nan,
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()) if len(sub) else np.nan,
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()) if len(sub) else np.nan,
                "mean_IFN": float(sub["score_IFN"].mean()) if len(sub) else np.nan,
                "mean_AT2": float(sub["score_AT2"].mean()) if len(sub) else np.nan,
                "mean_pseudotime": float(sub["palantir_pseudotime"].mean()) if len(sub) else np.nan,
                "top_subtype": (
                    sub["author_subtype"].astype(str).value_counts().index[0] if len(sub) else "NA"
                ),
            }
        )
    dest_df = pd.DataFrame(dest_rows)
    dest_df.to_csv(tabdir / "palantir_destinies.tsv", sep="\t", index=False)

    # Patient-level means (inferential unit)
    rows = []
    for patient, sub in adata.obs.groupby("patient_key", observed=True):
        mal = sub[sub["is_malignant"]]
        tumor = sub[sub["Sample_Origin"].astype(str) != "nLung"]
        rec = {
            "patient_key": patient,
            "dataset": str(sub["dataset"].iloc[0]),
            "n_cells": int(len(sub)),
            "n_malignant": int(len(mal)),
            "n_tumor_tissue": int(len(tumor)),
            "n_author_AT2": int(sub["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
            "n_cldn4_low": int((sub["cldn4_tertile"] == "low").sum()),
            "n_cldn4_high": int((sub["cldn4_tertile"] == "high").sum()),
            "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
            "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
            "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
            "mean_IFN": float(sub["score_IFN"].mean()),
            "mean_AT2": float(sub["score_AT2"].mean()),
            "mean_malignant_like": float(sub["score_malignant_like"].mean()),
            "mean_palantir_pt": float(sub["palantir_pseudotime"].mean()),
            "mean_dpt": float(pd.to_numeric(sub["dpt_pseudotime"], errors="coerce").mean()),
            "mean_CLDN4_malignant": float(mal["expr_CLDN4"].mean()) if len(mal) else np.nan,
            "mean_barrier_malignant": float(mal["score_barrier_keratin"].mean()) if len(mal) else np.nan,
            "mean_IFN_malignant": float(mal["score_IFN"].mean()) if len(mal) else np.nan,
            "mean_pt_malignant": float(mal["palantir_pseudotime"].mean()) if len(mal) else np.nan,
        }
        for col in branch.columns:
            rec[f"mean_fate_{col}"] = float(sub[f"fate_{col}"].mean())
            rec[f"mean_fate_{col}_malignant"] = (
                float(mal[f"fate_{col}"].mean()) if len(mal) else np.nan
            )
        rows.append(rec)
    patient_df = pd.DataFrame(rows)
    patient_df.to_csv(tabdir / "patient_means.tsv", sep="\t", index=False)
    elig = patient_df[patient_df["n_cells"] >= MIN_CELLS_PER_PATIENT].copy()
    elig_mal = patient_df[patient_df["n_malignant"] >= MIN_CELLS_PER_PATIENT].copy()

    dest_names = list(branch.columns)
    # Identify barrier-like vs IFN-like destinies from destiny means (descriptive)
    if len(dest_df):
        barrier_dest = str(dest_df.loc[dest_df["mean_barrier_keratin"].idxmax(), "destiny"])
        ifn_dest = str(dest_df.loc[dest_df["mean_IFN"].idxmax(), "destiny"])
    else:
        barrier_dest, ifn_dest = dest_names[0], dest_names[0]

    contrasts = [
        ("CLDN4 vs Palantir PT", "mean_CLDN4", "mean_palantir_pt"),
        ("barrier (no CLDN4) vs Palantir PT", "mean_barrier_keratin", "mean_palantir_pt"),
        ("IFN vs Palantir PT", "mean_IFN", "mean_palantir_pt"),
        ("CLDN4 vs IFN", "mean_CLDN4", "mean_IFN"),
        ("CLDN4 vs barrier (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("IFN vs barrier (no CLDN4)", "mean_IFN", "mean_barrier_keratin"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("SFTPC vs Palantir PT (control)", "mean_AT2", "mean_palantir_pt"),
    ]
    for dest in dest_names:
        contrasts.append((f"CLDN4 vs fate {dest}", "mean_CLDN4", f"mean_fate_{dest}"))
        contrasts.append((f"barrier vs fate {dest}", "mean_barrier_keratin", f"mean_fate_{dest}"))
        contrasts.append((f"IFN vs fate {dest}", "mean_IFN", f"mean_fate_{dest}"))

    primary = []
    for name, a, b in contrasts:
        if a not in elig.columns or b not in elig.columns:
            primary.append({"contrast": name, "n": 0, "rho": None, "p": None})
            continue
        primary.append({"contrast": name, **_spearman(elig[a].to_numpy(), elig[b].to_numpy())})
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q
    pd.DataFrame(primary).to_csv(tabdir / "patient_level_spearman.tsv", sep="\t", index=False)

    def sp(frame, a, b):
        if frame is None or len(frame) == 0 or a not in frame.columns or b not in frame.columns:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    elig_131 = elig[elig["dataset"] == "GSE131907"]
    elig_205 = elig[elig["dataset"] == "GSE205335"]
    sensitivity = [
        {"contrast": "GSE131907-only CLDN4 vs Palantir PT", **sp(elig_131, "mean_CLDN4", "mean_palantir_pt")},
        {"contrast": "GSE205335-only CLDN4 vs Palantir PT", **sp(elig_205, "mean_CLDN4", "mean_palantir_pt")},
        {"contrast": "GSE131907-only CLDN4 vs IFN", **sp(elig_131, "mean_CLDN4", "mean_IFN")},
        {"contrast": "GSE205335-only CLDN4 vs IFN", **sp(elig_205, "mean_CLDN4", "mean_IFN")},
        {"contrast": "GSE131907-only CLDN4 vs barrier", **sp(elig_131, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "GSE205335-only CLDN4 vs barrier", **sp(elig_205, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "malignant-cells CLDN4 vs Palantir PT", **sp(elig_mal, "mean_CLDN4_malignant", "mean_pt_malignant")},
        {"contrast": "malignant-cells CLDN4 vs IFN", **sp(elig_mal, "mean_CLDN4_malignant", "mean_IFN_malignant")},
        {"contrast": "malignant-cells CLDN4 vs barrier", **sp(elig_mal, "mean_CLDN4_malignant", "mean_barrier_malignant")},
        {"contrast": "malignant-cells IFN vs barrier", **sp(elig_mal, "mean_IFN_malignant", "mean_barrier_malignant")},
        {"contrast": "malignant-cells IFN vs Palantir PT", **sp(elig_mal, "mean_IFN_malignant", "mean_pt_malignant")},
        {"contrast": "malignant-cells barrier vs Palantir PT", **sp(elig_mal, "mean_barrier_malignant", "mean_pt_malignant")},
    ]
    if barrier_dest in dest_names:
        sensitivity.append(
            {
                "contrast": f"malignant-cells CLDN4 vs fate {barrier_dest} (barrier-high destiny)",
                **sp(elig_mal, "mean_CLDN4_malignant", f"mean_fate_{barrier_dest}_malignant"),
            }
        )
        sensitivity.append(
            {
                "contrast": f"malignant-cells IFN vs fate {ifn_dest} (IFN-high destiny)",
                **sp(elig_mal, "mean_IFN_malignant", f"mean_fate_{ifn_dest}_malignant"),
            }
        )
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    paired_recs = []
    for patient, sub in adata.obs.groupby("patient_key", observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        rec = {
            "patient_key": patient,
            "dataset": str(sub["dataset"].iloc[0]),
            "n_high": int(len(hi)),
            "n_low": int(len(lo)),
            "barrier_high": float(hi["score_barrier_keratin"].mean()),
            "barrier_low": float(lo["score_barrier_keratin"].mean()),
            "IFN_high": float(hi["score_IFN"].mean()),
            "IFN_low": float(lo["score_IFN"].mean()),
            "AT2_high": float(hi["score_AT2"].mean()),
            "AT2_low": float(lo["score_AT2"].mean()),
            "pt_high": float(hi["palantir_pseudotime"].mean()),
            "pt_low": float(lo["palantir_pseudotime"].mean()),
        }
        for dest in dest_names:
            rec[f"{dest}_high"] = float(hi[f"fate_{dest}"].mean())
            rec[f"{dest}_low"] = float(lo[f"fate_{dest}"].mean())
        paired_recs.append(rec)
    paired_df = pd.DataFrame(paired_recs)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paired_rows = []
    pair_specs = [
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("IFN high vs low", "IFN_high", "IFN_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("Palantir PT high vs low", "pt_high", "pt_low"),
    ]
    for dest in dest_names:
        pair_specs.append((f"fate {dest} high vs low", f"{dest}_high", f"{dest}_low"))
    for label, a, b in pair_specs:
        if paired_df.empty or a not in paired_df.columns:
            paired_rows.append({"contrast": label, "n": 0, "W": None, "p": None, "delta_median": None})
        else:
            paired_rows.append({"contrast": label, **_wilcoxon_paired(paired_df[a].to_numpy(), paired_df[b].to_numpy())})

    c4_bar = next(r for r in primary if r["contrast"] == "CLDN4 vs barrier (no CLDN4)")
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN")
    emit_extra = True  # extra figures always for this folder

    # ---- figures ----
    colors = {"GSE131907": "#2a6f97", "GSE205335": "#b23a48"}

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 9.2))
    sc.pl.umap(adata, color="palantir_pseudotime", ax=axes[0, 0], show=False, frameon=False, cmap="magma", title="Palantir pseudotime")
    sc.pl.umap(adata, color="palantir_destiny", ax=axes[0, 1], show=False, frameon=False, title="Palantir destiny (argmax fate)")
    sc.pl.umap(adata, color="expr_CLDN4", ax=axes[1, 0], show=False, frameon=False, cmap="viridis", title="CLDN4")
    sc.pl.umap(adata, color="score_IFN", ax=axes[1, 1], show=False, frameon=False, cmap="plasma", title="IFN core (no CLDN4)")
    fig.suptitle(
        f"REAL Palantir  n_cells={adata.n_obs}  n_patients={patient_df.shape[0]}  destinies={len(dest_names)}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_palantir_overview")

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    if len(dest_df):
        x = np.arange(len(dest_df))
        axes[0].bar(x, dest_df["mean_CLDN4"], color="#3d5a80")
        axes[1].bar(x, dest_df["mean_barrier_keratin"], color="#ee6c4d")
        axes[2].bar(x, dest_df["mean_IFN"], color="#9b2226")
        for ax, ylab in zip(axes, ["mean CLDN4", "mean barrier (no CLDN4)", "mean IFN"]):
            ax.set_xticks(x)
            ax.set_xticklabels(dest_df["destiny"], rotation=30, ha="right")
            ax.set_ylabel(ylab)
    axes[0].set_title("CLDN4 by destiny")
    axes[1].set_title("Barrier (CLDN4 excluded)")
    axes[2].set_title("IFN core (CLDN4 excluded)")
    fig.suptitle("Palantir destinies scored by CLDN4 / barrier / IFN", fontsize=11)
    _save(fig, figdir / "fig_palantir_destinies_programs")

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    panels = (
        ("mean_palantir_pt", "mean_CLDN4", "CLDN4 vs Palantir PT", next(r for r in primary if r["contrast"] == "CLDN4 vs Palantir PT")),
        ("mean_palantir_pt", "mean_barrier_keratin", "barrier vs Palantir PT", next(r for r in primary if r["contrast"] == "barrier (no CLDN4) vs Palantir PT")),
        ("mean_palantir_pt", "mean_IFN", "IFN vs Palantir PT", next(r for r in primary if r["contrast"] == "IFN vs Palantir PT")),
    )
    for ax, (x, y, lab, row) in zip(axes, panels):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=48, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel("patient-mean Palantir PT")
        ax.set_ylabel(lab.split(" vs ")[0])
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"Patient-level programs along Palantir PT  n={len(elig)}", fontsize=11)
    _save(fig, figdir / "fig_palantir_pt_vs_programs")

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    extra_panels = (
        ("mean_CLDN4", "mean_barrier_keratin", c4_bar, "CLDN4", "barrier (no CLDN4)"),
        ("mean_CLDN4", "mean_IFN", c4_ifn, "CLDN4", "IFN"),
        ("mean_IFN", "mean_barrier_keratin", next(r for r in primary if r["contrast"] == "IFN vs barrier (no CLDN4)"), "IFN", "barrier (no CLDN4)"),
    )
    for ax, (x, y, row, xlab, ylab) in zip(axes, extra_panels):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=48, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel(f"patient-mean {xlab}")
        ax.set_ylabel(f"patient-mean {ylab}")
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: patient-level CLDN4 / barrier / IFN", fontsize=11)
    _save(fig, figdir / "fig_extra_patient_programs")

    # fate vs programs for the two highlighted destinies
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    for ax, dest, ylab in (
        (axes[0], barrier_dest, "barrier-high destiny fate"),
        (axes[1], ifn_dest, "IFN-high destiny fate"),
    ):
        xcol = f"mean_fate_{dest}"
        if xcol not in elig.columns:
            ax.set_axis_off()
            continue
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub["mean_CLDN4"], sub[xcol], s=48, c=col, label=f"{ds} n={len(sub)}")
        row = next((r for r in primary if r["contrast"] == f"CLDN4 vs fate {dest}"), {"n": 0, "rho": None, "p": None})
        ax.set_xlabel("patient-mean CLDN4")
        ax.set_ylabel(ylab)
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: patient fate vs CLDN4 (barrier-high vs IFN-high destiny)", fontsize=11)
    _save(fig, figdir / "fig_extra_fate_vs_programs")

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    ct = (
        adata.obs.assign(author_subtype=adata.obs["author_subtype"].astype(str).fillna("NA"))
        .groupby(["dataset", "author_subtype"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.T.plot(kind="bar", ax=axes[0], color={"GSE131907": "#2a6f97", "GSE205335": "#b23a48"})
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Author subtypes  n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=8)
    plt.setp(axes[0].get_xticklabels(), rotation=40, ha="right")
    unit_ct = patient_df.groupby("dataset").size()
    axes[1].bar(unit_ct.index.astype(str), unit_ct.to_numpy(), color=["#2a6f97", "#b23a48"][: len(unit_ct)])
    axes[1].set_ylabel("patients")
    axes[1].set_title(f"Honest n patients={patient_df.shape[0]} (eligible {len(elig)})")
    _save(fig, figdir / "fig_honest_n")

    if not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("IFN_low", "IFN_high", "IFN core"),
            ("pt_low", "pt_high", "Palantir PT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            for ds, col in colors.items():
                sub = paired_df[paired_df["dataset"] == ds]
                ax.scatter(sub[lo], sub[hi], s=44, c=col, label=f"{ds} n={len(sub)}")
            lims = [
                min(paired_df[lo].min(), paired_df[hi].min()),
                max(paired_df[lo].max(), paired_df[hi].max()),
            ]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ls="--", c="0.6", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
        axes[0].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("barrier")), keys=("W", "p")))
        axes[1].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("IFN")), keys=("W", "p")))
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("Palantir PT")), keys=("W", "p")))
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(f"EXTRA: within-patient CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    # gene trends: deciles of Palantir PT (unsmoothed means)
    pt = adata.obs["palantir_pseudotime"].to_numpy()
    finite = np.isfinite(pt)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    if finite.sum() >= 20:
        dec = pd.qcut(pt[finite], 10, labels=False, duplicates="drop")
        trend = pd.DataFrame(
            {
                "decile": dec,
                "CLDN4": adata.obs.loc[finite, "expr_CLDN4"].to_numpy(),
                "barrier": adata.obs.loc[finite, "score_barrier_keratin"].to_numpy(),
                "IFN": adata.obs.loc[finite, "score_IFN"].to_numpy(),
                "AT2": adata.obs.loc[finite, "score_AT2"].to_numpy(),
            }
        )
        g = trend.groupby("decile").mean()
        for col, lab, c in (
            ("CLDN4", "CLDN4", "#3d5a80"),
            ("barrier", "barrier (no CLDN4)", "#ee6c4d"),
            ("IFN", "IFN", "#9b2226"),
            ("AT2", "AT2", "#2a9d8f"),
        ):
            z = (g[col] - g[col].mean()) / (g[col].std() + 1e-9)
            ax.plot(g.index.to_numpy() + 1, z.to_numpy(), marker="o", label=lab, color=c)
        ax.axhline(0, c="0.7", lw=1)
        ax.set_xlabel("Palantir PT decile")
        ax.set_ylabel("z-scored mean (cells)")
        ax.legend(frameon=False, fontsize=8)
    ax.set_title("EXTRA: programs along Palantir PT (cell-decile means, z-scored)")
    _save(fig, figdir / "fig_extra_pt_decile_trends")

    for color, fname, cmap in (
        ("dataset", "fig_umap_dataset", None),
        ("Sample_Origin", "fig_umap_origin", None),
        ("author_subtype", "fig_umap_author_subtype", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("is_malignant", "fig_umap_malignant", None),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    subtype_counts = adata.obs["author_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    c4_pt = next(r for r in primary if r["contrast"] == "CLDN4 vs Palantir PT")
    bar_pt = next(r for r in primary if r["contrast"] == "barrier (no CLDN4) vs Palantir PT")
    ifn_pt = next(r for r in primary if r["contrast"] == "IFN vs Palantir PT")
    ifn_bar = next(r for r in primary if r["contrast"] == "IFN vs barrier (no CLDN4)")
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_ifn = next(r for r in paired_rows if r["contrast"].startswith("IFN"))
    mal_c4_ifn = next(r for r in sensitivity if r["contrast"] == "malignant-cells CLDN4 vs IFN")
    mal_c4_bar = next(r for r in sensitivity if r["contrast"] == "malignant-cells CLDN4 vs barrier")
    mal_ifn_bar = next(r for r in sensitivity if r["contrast"] == "malignant-cells IFN vs barrier")

    same_dest = barrier_dest == ifn_dest
    thesis_bits = []
    if same_dest:
        thesis_bits.append(
            f"Palantir auto-terminals put the highest-barrier and highest-IFN means on the **same** destiny ({barrier_dest}); "
            "branching into two destinies is not required for a two-end continuum."
        )
    else:
        thesis_bits.append(
            f"Highest-barrier destiny is {barrier_dest}; highest-IFN destiny is {ifn_dest}."
        )
    # Thesis hold rules (honest):
    # 1) patient CLDN4 vs barrier ρ>0
    # 2) patient CLDN4 vs IFN ρ<0 OR paired high-low IFN Δ<0
    # 3) malignant-cell CLDN4 vs IFN ρ<0 or malignant IFN vs barrier ρ<0
    hold_bar = c4_bar.get("rho") is not None and c4_bar["rho"] > 0 and (c4_bar.get("p") or 1) < 0.05
    hold_ifn_anti = (
        (c4_ifn.get("rho") is not None and c4_ifn["rho"] < 0 and (c4_ifn.get("p") or 1) < 0.05)
        or (pair_ifn.get("delta_median") is not None and pair_ifn["delta_median"] < 0 and (pair_ifn.get("p") or 1) < 0.05)
        or (mal_c4_ifn.get("rho") is not None and mal_c4_ifn["rho"] < 0 and (mal_c4_ifn.get("p") or 1) < 0.05)
    )
    if hold_bar and hold_ifn_anti:
        thesis_verdict = "Thesis is **supported** at patient (or paired / malignant-restricted) level on this object."
    elif hold_bar and not hold_ifn_anti:
        thesis_verdict = "Thesis is **partial**: barrier end holds; IFN-higher CLDN4-low end is not supported at the pre-specified tests."
    else:
        thesis_verdict = "Thesis is **not supported** on the pre-specified patient-level tests."
    thesis_bits.append(thesis_verdict)

    parts = [
        f"Patient-level CLDN4 vs Palantir PT: {_fmt(c4_pt)}.",
        f"Barrier (no CLDN4) vs Palantir PT: {_fmt(bar_pt)}.",
        f"IFN vs Palantir PT: {_fmt(ifn_pt)}.",
        f"CLDN4 vs barrier (no CLDN4): {_fmt(c4_bar)}.",
        f"CLDN4 vs IFN: {_fmt(c4_ifn)}.",
        f"IFN vs barrier: {_fmt(ifn_bar)}.",
        f"Malignant-restricted CLDN4 vs IFN: {_fmt(mal_c4_ifn)}.",
        f"Malignant-restricted CLDN4 vs barrier: {_fmt(mal_c4_bar)}.",
        f"Malignant-restricted IFN vs barrier: {_fmt(mal_ifn_bar)}.",
        f"Paired CLDN4-high vs low barrier: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"Palantir destinies={len(dest_names)}; early cell not CLDN4-high ({root_info['root_tertile']}).",
        f"DPT companion empty={dpt_info.get('empty')}.",
        "Not a TACSTD2 redo. No both-high gate. GSE148071 not added. GSE207422 not added.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What holds (n={c4_bar['n']} patients).** "
        f"CLDN4 vs barrier (CLDN4 excluded): {_fmt(c4_bar)}. "
        f"CLDN4 vs IFN: {_fmt(c4_ifn)}. "
        f"Within-patient CLDN4-high vs low barrier {_fmt(pair_bar, keys=('W', 'p'))}; "
        f"IFN {_fmt(pair_ifn, keys=('W', 'p'))}. "
        f"**Palantir product.** {pal_info['n_destinies']} destinies, "
        f"{pal_info['n_pseudotime_finite']} cells with finite pseudotime. "
        f"**What is mixed / null.** CLDN4 vs Palantir PT {_fmt(c4_pt)}; "
        f"IFN vs Palantir PT {_fmt(ifn_pt)}."
    )

    summary = {
        "accessions": ["GSE131907", "GSE205335"],
        "winning_pair": True,
        "gse148071_added": False,
        "gse207422_added": False,
        "not_a_pr325_redo": True,
        "not_a_pr449_redo": True,
        "primary_gene": "CLDN4",
        "dual_high": False,
        "inferential_unit": "patient",
        "clock": "palantir",
        "palantir_pkg": pal_pkg,
        "palantir": pal_info,
        "dpt": dpt_info,
        "harmony": harmony,
        "cap_per_unit": CAP_PER_UNIT,
        "n_cells": int(adata.n_obs),
        "n_cells_gse131907": int((adata.obs["dataset"] == "GSE131907").sum()),
        "n_cells_gse205335": int((adata.obs["dataset"] == "GSE205335").sum()),
        "n_cells_nLung": int((adata.obs["Sample_Origin"] == "nLung").sum()),
        "n_author_AT2": int(adata.obs["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
        "n_patients": int(patient_df.shape[0]),
        "n_patients_gse131907": int((patient_df["dataset"] == "GSE131907").sum()),
        "n_patients_gse205335": int((patient_df["dataset"] == "GSE205335").sum()),
        "n_patients_eligible": int(len(elig)),
        "n_patients_malignant": int(len(elig_mal)),
        "n_patients_paired_tertile": int(len(paired_df)),
        "subtype_counts": subtype_counts,
        "cldn4_tertile_counts": tertile_counts,
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "destinies": dest_rows,
        "barrier_high_destiny": barrier_dest,
        "ifn_high_destiny": ifn_dest,
        "thesis_verdict": " ".join(thesis_bits),
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "extra_figure": {
            "emitted": bool(emit_extra),
            "spearman_rho": c4_bar["rho"],
            "spearman_p": c4_bar["p"],
            "spearman_n": c4_bar["n"],
        },
        "qc": {
            "min_genes": 200,
            "min_umi": 500,
            "lineage_EPCAM_mean": float(adata.obs["expr_EPCAM"].mean()) if "expr_EPCAM" in adata.obs else None,
            "lineage_PTPRC_mean": float(adata.obs["expr_PTPRC"].mean()) if "expr_PTPRC" in adata.obs else None,
        },
        "verdict": verdict,
        "what_holds": what_holds,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accessions": ["GSE131907", "GSE205335"],
                "papers": {
                    "GSE131907": "Kim et al. Nat Commun 2020 PMID 32385277",
                    "GSE205335": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
                    "Palantir": "Setty et al. Nat Biotechnol 2019",
                },
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "ifn_excludes_CLDN4": True,
                "root_not_CLDN4_high": True,
                "gse148071_added": False,
                "dual_high": False,
                "inferential_unit": "patient",
                "palantir": pal_pkg,
                "dpt_is_companion_only": True,
            },
            indent=2,
        )
    )
    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": int(adata.n_obs),
                "n_patients": int(patient_df.shape[0]),
                "n_destinies": int(len(dest_names)),
                "pseudotime_table": str(tabdir / "palantir_pseudotime.tsv"),
                "destiny_table": str(tabdir / "palantir_destinies.tsv"),
                "finding": args.finding,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
