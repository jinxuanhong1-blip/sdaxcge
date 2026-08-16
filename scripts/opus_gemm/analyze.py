#!/usr/bin/env python3
"""Tacstd2 / Cldn4 versus ICI treatment, genotype and immune signatures.

All p-values are computed from the per-dataset log2 matrices produced by
build_matrices.py. Contrasts with n<2 in either arm are reported as
descriptive (no p-value). Multiple-testing correction is Benjamini-Hochberg
within a pre-specified family (ici / genotype / resistance / mechanism).

Usage:
    python3 scripts/opus_gemm/analyze.py \
        --matrices results/opus_gemm/data/matrices \
        --outdir results/opus_gemm
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from datasets import DATASETS, Dataset  # noqa: E402
from genes import FOCUS, SIGNATURES, load_ids  # noqa: E402


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    """Hedges' g (small-sample corrected Cohen's d), a minus b."""
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return float("nan")
    s1, s2 = float(np.var(a, ddof=1)), float(np.var(b, ddof=1))
    pooled = math.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    if pooled == 0:
        return 0.0
    d = (float(np.mean(a)) - float(np.mean(b))) / pooled
    j = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0) if (n1 + n2) > 2 else 1.0
    return d * j


def find_gene(expr: pd.DataFrame, symbol: str, gene_id_type: str, ids: dict[str, str]) -> str | None:
    if gene_id_type == "ensembl":
        eid = ids.get(symbol)
        if eid and eid in expr.index:
            return eid
        # version-stripped already; try prefix
        hits = [i for i in expr.index if str(i).startswith(eid or "___never___")]
        return hits[0] if hits else None
    upper = {str(i).upper(): i for i in expr.index}
    return upper.get(symbol.upper())


def signature_scores(expr: pd.DataFrame, gene_id_type: str, ids: dict[str, str]) -> pd.DataFrame:
    z = (expr - expr.mean(axis=1).values[:, None]) / expr.std(axis=1, ddof=1).replace(0, np.nan).values[:, None]
    out = {}
    for name, members in SIGNATURES.items():
        rows = [find_gene(expr, g, gene_id_type, ids) for g in members]
        rows = [r for r in rows if r is not None]
        if len(rows) < 3:
            continue
        out[name] = z.loc[rows].mean(axis=0)
    return pd.DataFrame(out)


def contrast_test(values: pd.Series, design: pd.DataFrame, contrast, paired_by: str | None) -> dict:
    test_cols = design.loc[design["group"] == contrast.test, "column"].tolist()
    ref_cols = design.loc[design["group"] == contrast.ref, "column"].tolist()
    a = values.reindex(test_cols).dropna()
    b = values.reindex(ref_cols).dropna()
    row = {
        "n_test": int(a.size),
        "n_ref": int(b.size),
        "mean_test": float(a.mean()) if a.size else float("nan"),
        "mean_ref": float(b.mean()) if b.size else float("nan"),
        "sd_test": float(a.std(ddof=1)) if a.size > 1 else float("nan"),
        "sd_ref": float(b.std(ddof=1)) if b.size > 1 else float("nan"),
        "delta_log2": float(a.mean() - b.mean()) if a.size and b.size else float("nan"),
        "hedges_g": hedges_g(a.to_numpy(), b.to_numpy()),
        "test_name": "",
        "statistic": float("nan"),
        "p_value": float("nan"),
        "note": "",
    }
    if a.size < 2 or b.size < 2:
        row["note"] = "n<2 in at least one arm; descriptive only"
        return row

    if paired_by:
        keys_a = {re.search(paired_by, c).group(0) if re.search(paired_by, c) else c: c for c in a.index}
        # use the first capturing group(s) when present
        def pkey(c: str) -> str | None:
            m = re.search(paired_by, c)
            if not m:
                return None
            return m.group(1) if m.lastindex else m.group(0)

        pairs = []
        used_ref = set()
        for c in a.index:
            k = pkey(c)
            if k is None:
                continue
            matches = [r for r in b.index if pkey(r) == k and r not in used_ref]
            if matches:
                pairs.append((c, matches[0]))
                used_ref.add(matches[0])
        if len(pairs) >= 2 and len(pairs) == a.size == b.size:
            av = values.loc[[p[0] for p in pairs]].to_numpy()
            bv = values.loc[[p[1] for p in pairs]].to_numpy()
            t, p = stats.ttest_rel(av, bv)
            row.update(test_name="paired_t", statistic=float(t), p_value=float(p),
                       n_test=len(pairs), n_ref=len(pairs),
                       delta_log2=float(np.mean(av - bv)),
                       hedges_g=hedges_g(av, bv))
            return row
        row["note"] = "pairing incomplete; fell back to Welch t"

    t, p = stats.ttest_ind(a.to_numpy(), b.to_numpy(), equal_var=False)
    row.update(test_name="welch_t", statistic=float(t), p_value=float(p))
    if a.size < 3 or b.size < 3:
        row["note"] = (row["note"] + "; n<3 per arm, underpowered").strip("; ")
    return row


def bh(pvals: list[float]) -> list[float]:
    mask = [p == p and p is not None for p in pvals]  # not NaN
    out = [float("nan")] * len(pvals)
    if sum(mask) == 0:
        return out
    _, adj, _, _ = multipletests([pvals[i] for i, m in enumerate(mask) if m], method="fdr_bh")
    j = 0
    for i, m in enumerate(mask):
        if m:
            out[i] = float(adj[j])
            j += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrices", default="results/opus_gemm/data/matrices")
    ap.add_argument("--ids", default="results/opus_gemm/gene_ids.tsv")
    ap.add_argument("--outdir", default="results/opus_gemm")
    args = ap.parse_args()

    ids = load_ids(args.ids)
    os.makedirs(args.outdir, exist_ok=True)

    detection_rows, value_rows, contrast_rows, corr_rows = [], [], [], []

    for ds in DATASETS:
        expr_path = os.path.join(args.matrices, f"{ds.key}.expr.tsv.gz")
        design_path = os.path.join(args.matrices, f"{ds.key}.design.tsv")
        if not os.path.exists(expr_path):
            print(f"[skip] {ds.key}: no matrix", file=sys.stderr)
            continue
        expr = pd.read_csv(expr_path, sep="\t", index_col=0)
        design = pd.read_csv(design_path, sep="\t")
        design = design[design["column"].isin(expr.columns)]

        gene_rows = {}
        for sym in FOCUS:
            gid = find_gene(expr, sym, ds.gene_id_type, ids)
            detection_rows.append({
                "dataset": ds.key, "accession": ds.acc, "symbol": sym,
                "resolved_id": gid or "", "detected": "yes" if gid else "no",
                "median_log2": f"{float(expr.loc[gid].median()):.4f}" if gid else "",
            })
            if gid:
                gene_rows[sym] = expr.loc[gid]

        sig = signature_scores(expr, ds.gene_id_type, ids)

        for _, s in design.iterrows():
            rec = {
                "dataset": ds.key, "accession": ds.acc, "column": s["column"],
                "gsm": s.get("gsm", ""), "group": s["group"],
                "model_family": ds.model_family, "setting": ds.setting,
            }
            for sym, ser in gene_rows.items():
                rec[sym] = f"{float(ser[s['column']]):.5f}"
            for name in sig.columns:
                rec[name] = f"{float(sig.loc[s['column'], name]):.5f}"
            value_rows.append(rec)

        for contrast in ds.contrasts:
            for sym, ser in gene_rows.items():
                stats_row = contrast_test(ser, design, contrast, ds.paired_by)
                contrast_rows.append({
                    "dataset": ds.key, "accession": ds.acc, "model_family": ds.model_family,
                    "setting": ds.setting, "symbol": sym, "contrast": contrast.name,
                    "family": contrast.family, "test_group": contrast.test, "ref_group": contrast.ref,
                    **{k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in stats_row.items()},
                    "contrast_note": contrast.note,
                })

        # Spearman: focus genes vs signatures (in-vivo / sorted only — cultured
        # lines have no immune infiltrate)
        if ds.setting != "cultured_cells" and not sig.empty:
            for sym, ser in gene_rows.items():
                x = ser.reindex(design["column"]).to_numpy(dtype=float)
                for name in sig.columns:
                    y = sig.reindex(design["column"])[name].to_numpy(dtype=float)
                    mask = np.isfinite(x) & np.isfinite(y)
                    if mask.sum() < 5:
                        continue
                    rho, p = stats.spearmanr(x[mask], y[mask])
                    corr_rows.append({
                        "dataset": ds.key, "accession": ds.acc, "model_family": ds.model_family,
                        "symbol": sym, "signature": name, "n": int(mask.sum()),
                        "spearman_rho": f"{float(rho):.5f}", "p_value": f"{float(p):.6g}",
                    })

        print(f"[ok] {ds.key}  focus={list(gene_rows)}  signatures={list(sig.columns)}")

    # BH within family for focus-gene contrasts
    by_fam: dict[str, list[int]] = {}
    for i, r in enumerate(contrast_rows):
        by_fam.setdefault(r["family"], []).append(i)
    for fam, idxs in by_fam.items():
        pvals = [float(contrast_rows[i]["p_value"]) if contrast_rows[i]["p_value"] != "" else float("nan")
                 for i in idxs]
        adj = bh(pvals)
        for i, q in zip(idxs, adj):
            contrast_rows[i]["fdr_bh_within_family"] = f"{q:.6g}" if q == q else ""

    if corr_rows:
        pvals = [float(r["p_value"]) for r in corr_rows]
        adj = bh(pvals)
        for r, q in zip(corr_rows, adj):
            r["fdr_bh"] = f"{q:.6g}" if q == q else ""

    def write_tsv(path: str, rows: list[dict]) -> None:
        if not rows:
            return
        fields = list(rows[0])
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {path}  ({len(rows)} rows)")

    write_tsv(os.path.join(args.outdir, "gene_detection.tsv"), detection_rows)
    write_tsv(os.path.join(args.outdir, "sample_values.tsv"), value_rows)
    write_tsv(os.path.join(args.outdir, "focus_contrasts.tsv"), contrast_rows)
    write_tsv(os.path.join(args.outdir, "signature_correlations.tsv"), corr_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
