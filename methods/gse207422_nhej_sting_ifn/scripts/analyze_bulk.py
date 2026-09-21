#!/usr/bin/env python3
"""Baseline bulk log2TPM on GSE207422: NHEJ / STING / IFN vs TACSTD2 and CLDN4.

Secondary to the single-cell split. These are pre-treatment biopsies from the
same series, not the 12 post-surgery scRNA tumors. EPCAM is a purity control.
Pathologic response is not tested here.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import MODULES, PRIMARY_MODULES, STING_ALIASES, requested_symbols

HERE = Path(__file__).resolve().parents[1]


def resolve(present: dict[str, np.ndarray], name: str) -> str | None:
    options = STING_ALIASES.get(name, (name,))
    hits = [g for g in options if g in present]
    if not hits:
        return None
    return max(hits, key=lambda g: float(np.nansum(present[g])))


def spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(len(x))
    if n < 5 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p), n


def partial_spearman(x, y, z) -> tuple[float, float]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]
    if len(x) < 6 or np.unique(z).size < 2:
        return np.nan, np.nan

    def resid(a, cov):
        A = np.column_stack([np.ones(len(cov)), stats.rankdata(cov)])
        coef, *_ = np.linalg.lstsq(A, stats.rankdata(a), rcond=None)
        return stats.rankdata(a) - A @ coef

    rx, ry = resid(x, z), resid(y, z)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan, np.nan
    rho, p = stats.pearsonr(rx, ry)
    return float(rho), float(p)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()
    expr_path = args.workdir / "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz"
    meta_path = args.workdir / "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx"
    tabdir = args.outdir / "tables"
    figdir = args.outdir / "figures"
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    if not expr_path.exists():
        print("bulk expression missing; skip", flush=True)
        return

    want = set(requested_symbols()) | {"EPCAM"}
    with gzip.open(expr_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = [s.strip('"') for s in header[1:]]
        rows: dict[str, np.ndarray] = {}
        ngenes = 0
        for line in fh:
            ngenes += 1
            parts = line.rstrip("\n").split("\t")
            gene = parts[0].strip('"')
            if gene in want:
                rows[gene] = np.array([float(x) for x in parts[1:]])
    print(f"bulk genes scanned {ngenes}; kept {sorted(rows)}", flush=True)

    if "TACSTD2" not in rows or "CLDN4" not in rows:
        pd.DataFrame([{
            "status": "gate_stop",
            "TACSTD2_present": "TACSTD2" in rows,
            "CLDN4_present": "CLDN4" in rows,
            "n_genes": ngenes,
        }]).to_csv(tabdir / "bulk_gate.tsv", sep="\t", index=False)
        print("bulk gate stop", flush=True)
        return

    meta = pd.read_excel(meta_path)
    meta = meta.dropna(subset=["Sample"]).copy()
    meta["Sample"] = meta["Sample"].astype(str).str.strip()
    meta = meta[~meta["Sample"].str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    meta["timing"] = np.where(
        meta["Resource"].astype(str).str.contains("Pre", case=False, na=False),
        "pre",
        "post",
    )
    meta = meta.drop_duplicates("Sample").set_index("Sample")

    common = [s for s in samples if s in meta.index]
    idx = [samples.index(s) for s in common]
    df = meta.loc[common, ["Patient", "Resource", "timing", "Pathology", "Pathologic Response"]].copy()
    df["EPCAM"] = rows["EPCAM"][idx] if "EPCAM" in rows else np.nan
    df["CLDN4"] = rows["CLDN4"][idx]
    df["TACSTD2"] = rows["TACSTD2"][idx]

    used = {}
    for name, glist in MODULES.items():
        cols = []
        symbols = []
        for g in glist:
            sym = resolve(rows, g)
            if sym is None:
                continue
            symbols.append(sym)
            cols.append(rows[sym][idx])
        used[name] = symbols
        df[name] = np.mean(np.vstack(cols), axis=0) if cols else np.nan

    # Primary bulk cohort: pre-treatment biopsies that are in the matrix.
    pre = df[df["timing"] == "pre"].copy()
    df.to_csv(tabdir / "bulk_per_sample.tsv", sep="\t")

    rows_out = []
    for splitter in ("CLDN4", "TACSTD2"):
        for name in list(MODULES):
            rho, p, n = spearman(pre[splitter], pre[name])
            pr, pp = partial_spearman(pre[splitter], pre[name], pre["EPCAM"])
            rows_out.append({
                "cohort": "pre_biopsy",
                "splitter": splitter,
                "module": name,
                "n": n,
                "spearman_rho": rho,
                "spearman_p": p,
                "partial_epcam_rho": pr,
                "partial_epcam_p": pp,
                "n_genes_used": len(used[name]),
                "genes_used": ",".join(used[name]),
            })
        rho, p, n = spearman(pre["CLDN4"], pre["TACSTD2"])
        rows_out.append({
            "cohort": "pre_biopsy",
            "splitter": splitter,
            "module": "the_other_focal" if splitter == "CLDN4" else "CLDN4_repeat",
            "n": n,
            "spearman_rho": rho,
            "spearman_p": p,
            "partial_epcam_rho": partial_spearman(pre["CLDN4"], pre["TACSTD2"], pre["EPCAM"])[0],
            "partial_epcam_p": partial_spearman(pre["CLDN4"], pre["TACSTD2"], pre["EPCAM"])[1],
            "n_genes_used": 1,
            "genes_used": "TACSTD2,CLDN4",
        })
    out = pd.DataFrame(rows_out)
    # Drop the duplicate CLDN4-vs-TACSTD2 row produced by the TACSTD2 loop.
    out = out[out["module"] != "CLDN4_repeat"]
    out.to_csv(tabdir / "bulk_tests.tsv", sep="\t", index=False)

    gate = {
        "status": "run",
        "n_genes_in_matrix": ngenes,
        "n_samples_in_matrix": len(samples),
        "n_samples_with_meta": len(common),
        "n_pre": int((df["timing"] == "pre").sum()),
        "n_post": int((df["timing"] == "post").sum()),
        "TACSTD2_present": True,
        "CLDN4_present": True,
        "note": "Pre-treatment biopsies. Not the scRNA post-surgery cohort. Response is not tested.",
        "module_genes": used,
    }
    (tabdir / "bulk_summary.json").write_text(json.dumps(gate, indent=2))
    print(out[out["module"].isin(PRIMARY_MODULES)].to_string(index=False), flush=True)
    print(json.dumps({k: gate[k] for k in ("n_pre", "n_post", "n_samples_in_matrix")}, indent=2), flush=True)

    if len(pre) >= 5:
        fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.3))
        titles = {"nhej": "NHEJ", "sting": "STING", "ifn_effector": "IFN effectors"}
        for ax, name in zip(axes, PRIMARY_MODULES):
            ax.scatter(pre["CLDN4"], pre[name], s=22, c="#3d405b", zorder=3)
            ax.set_xlabel("CLDN4 log2TPM")
            ax.set_ylabel(titles[name])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        fig.suptitle(f"Baseline bulk, pre-treatment n={len(pre)}", y=1.02, fontsize=11)
        fig.tight_layout()
        fig.savefig(figdir / "fig_bulk_cldn4_modules.png", dpi=160, bbox_inches="tight")
        fig.savefig(figdir / "fig_bulk_cldn4_modules.pdf", bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
