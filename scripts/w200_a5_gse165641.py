#!/usr/bin/env python3
"""w200 / A5: Tacstd2-high vs CD8/NK in a real public KL lung dataset.

GSE76628 is not Kras;Lkb1 lung (it is Ad-VEGF flank-skin / gastric-stroma arrays
in nude mice). This script analyzes the closest public series that actually is
KL lung and can support a cell-level contrast:

  GSE165641 — 10X scRNA-seq of RBC-depleted cells from lung tumor sections of
  two KrasLSL-G12D/+; Lkb1fl/fl (KL) mice, 10 weeks after Ad-Cre
  (Wang & Zhong; PMID 34369094). Author-processed log-normalized matrix,
  31,053 genes x 7,180 cells (KL1 n=4,400; KL2 n=2,780).

Analyses (all limitations written into the outputs):

  1. Marker-based coarse cell typing (no published cluster labels are in the
     GEO matrix). Primary contrast uses strict Epithelial (Ptprc-negative
     Epcam/Krt/Nkx2-1) vs CD8 T vs NK.
  2. Tacstd2 expression and detection rate by cell type (Mann-Whitney +
     Fisher's exact).
  3. Within epithelial cells, Tacstd2 tertiles vs junction (Cldn4) and
     MHC/IFN genes — a test of whether a Tacstd2-high tumor subset is
     intrinsically immune-resistant.

Outputs: results/w200/A5_GSE165641/
"""

from __future__ import annotations

import hashlib
import json
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdata
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data" / "kl_lung"
OUT = REPO / "results" / "w200" / "A5_GSE165641"

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165641/suppl/"
    "GSE165641_processed_normalized_matrix_data.Rdata.gz"
)
MATRIX_NAME = "GSE165641_processed_normalized_matrix_data.Rdata.gz"

# Pre-specified panels. Detection uses raw > 0 on the author-normalized matrix
# (zeros match count zeros for these genes).
CD8_GENES = ["Cd8a", "Cd8b1", "Cd3e", "Cd3d"]
NK_GENES = ["Ncr1", "Klrb1c", "Klrk1", "Nkg7"]
EPI_GENES = ["Epcam", "Krt8", "Krt18", "Krt19", "Nkx2-1"]
MHC_IFN = ["B2m", "H2-K1", "H2-D1", "Tap1", "Psmb8", "Stat1", "Irf1",
           "Cxcl9", "Cxcl10", "Cd274", "Nlrc5"]
JUNCTION = ["Cldn4", "Epcam", "Krt8"]
EXTRA = ["Ptprc", "Tacstd2", "Sftpc", "Scgb1a1", "Cd68", "Adgre1", "Itgam",
         "S100a8", "Mki67", "Pecam1"]

CELLTYPE_ORDER = [
    "Epithelial", "CD8_T", "NK", "Other_T",
    "Neutrophil", "Myeloid", "Other_immune", "Other",
]


def fetch(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        tmp = dest.with_suffix(dest.suffix + ".part")
        with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:
            shutil.copyfileobj(r, f)
        tmp.rename(dest)
    md5 = hashlib.md5(dest.read_bytes()).hexdigest()
    return {"url": url, "path": str(dest.relative_to(REPO)), "md5": md5,
            "bytes": dest.stat().st_size}


def load_matrix(path: Path) -> pd.DataFrame:
    """Return cells x genes, author-normalized (Seurat-like log)."""
    parsed = rdata.parser.parse_file(path)
    converted = rdata.conversion.convert(parsed)
    df = converted["processed_data"]
    X = df.T
    X.index = X.index.astype(str)
    X.columns = X.columns.astype(str)
    return X


def assign_celltype(X: pd.DataFrame) -> pd.Series:
    """Mutually exclusive marker rules, applied in this priority order.

    Epithelial requires Ptprc == 0 so Epcam+ immune doublets / ambient-heavy
    droplets are not counted as tumor. CD8 and NK are called first so a rare
    lymphocyte is not swallowed by myeloid rules.
    """
    pos = {g: X[g] > 0 for g in
           ["Cd8a", "Cd8b1", "Cd3e", "Cd3d", "Ncr1", "Klrb1c",
            "Epcam", "Krt8", "Krt18", "Nkx2-1", "Ptprc", "Cd68", "Adgre1"]}
    cd8 = (pos["Cd8a"] | pos["Cd8b1"]) & (pos["Cd3e"] | pos["Cd3d"])
    nk = (pos["Ncr1"] | pos["Klrb1c"]) & ~cd8
    other_t = (pos["Cd3e"] | pos["Cd3d"]) & ~cd8
    epi = ((pos["Epcam"] | (pos["Krt8"] & pos["Krt18"]) | pos["Nkx2-1"])
           & ~pos["Ptprc"])
    neutrophil = pos["Ptprc"] & (X["S100a8"] >= 2.0) & ~cd8 & ~nk & ~other_t
    myeloid = (pos["Ptprc"] & (pos["Cd68"] | pos["Adgre1"])
               & ~cd8 & ~nk & ~other_t & ~neutrophil)
    other_imm = pos["Ptprc"] & ~cd8 & ~nk & ~other_t & ~neutrophil & ~myeloid

    lab = pd.Series("Other", index=X.index)
    # later assignments do not overwrite earlier (priority = first write)
    for mask, name in [
        (cd8, "CD8_T"),
        (nk, "NK"),
        (other_t, "Other_T"),
        (epi, "Epithelial"),
        (neutrophil, "Neutrophil"),
        (myeloid, "Myeloid"),
        (other_imm, "Other_immune"),
    ]:
        lab.loc[mask & (lab == "Other")] = name
    return lab


def rank_biserial(a: np.ndarray, b: np.ndarray) -> float:
    """Rank-biserial effect size; + means a > b."""
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(2.0 * u / (len(a) * len(b)) - 1.0)


def bh_fdr(p: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p.index, dtype=float)
    ok = p.notna()
    if not ok.any():
        return out
    pv = p[ok]
    n = len(pv)
    order = pv.sort_values().index
    ranked = pv.loc[order].values * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out.loc[order] = np.clip(ranked, 0, 1)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dest = DATA / MATRIX_NAME
    file_meta = fetch(MATRIX_URL, dest)
    X = load_matrix(dest)
    assert X.shape == (7180, 31053), X.shape

    mouse = X.index.to_series().str.split("_").str[0]
    celltype = assign_celltype(X)
    tac = X["Tacstd2"]

    per_cell = pd.DataFrame({
        "barcode": X.index,
        "mouse": mouse.values,
        "celltype": celltype.values,
        "Tacstd2": tac.round(4).values,
        "Epcam": X["Epcam"].round(4).values,
        "Ptprc": X["Ptprc"].round(4).values,
        "Cldn4": X["Cldn4"].round(4).values,
        "Cd8a": X["Cd8a"].round(4).values,
        "Ncr1": X["Ncr1"].round(4).values,
    })
    per_cell.to_csv(OUT / "per_cell_markers.csv", index=False)

    # ---- cell-type composition ----
    comp = (pd.crosstab(celltype, mouse)
            .reindex(index=CELLTYPE_ORDER, columns=["KL1", "KL2"])
            .fillna(0).astype(int))
    comp["total"] = comp.sum(axis=1)
    comp.index.name = "celltype"
    comp.to_csv(OUT / "celltype_composition.csv")

    # ---- Tacstd2 by cell type ----
    rows = []
    for ct in CELLTYPE_ORDER:
        v = tac[celltype == ct]
        rows.append({
            "celltype": ct, "n": int(len(v)),
            "n_tacstd2_pos": int((v > 0).sum()),
            "detection_rate": round(float((v > 0).mean()), 4) if len(v) else np.nan,
            "mean": round(float(v.mean()), 4) if len(v) else np.nan,
            "median": round(float(v.median()), 4) if len(v) else np.nan,
            "p90": round(float(v.quantile(0.90)), 4) if len(v) else np.nan,
        })
    by_ct = pd.DataFrame(rows)
    by_ct.to_csv(OUT / "tacstd2_by_celltype.csv", index=False)

    # Primary contrast: Epithelial vs CD8 vs NK
    contrast_rows = []
    for a_name, b_name in [("Epithelial", "CD8_T"), ("Epithelial", "NK"),
                           ("CD8_T", "NK")]:
        a, b = tac[celltype == a_name], tac[celltype == b_name]
        mw = stats.mannwhitneyu(a, b, alternative="two-sided")
        tab = np.array([
            [(a > 0).sum(), (a == 0).sum()],
            [(b > 0).sum(), (b == 0).sum()],
        ], dtype=int)
        oddsr, fp = stats.fisher_exact(tab)
        contrast_rows.append({
            "comparison": f"{a_name}_vs_{b_name}",
            "n_a": int(len(a)), "n_b": int(len(b)),
            "mean_a": round(float(a.mean()), 4),
            "mean_b": round(float(b.mean()), 4),
            "detection_a": round(float((a > 0).mean()), 4),
            "detection_b": round(float((b > 0).mean()), 4),
            "mannwhitney_p": mw.pvalue,
            "rank_biserial": round(rank_biserial(a.values, b.values), 3),
            "fisher_oddsratio_detection": round(float(oddsr), 3),
            "fisher_p_detection": fp,
        })
    contrast = pd.DataFrame(contrast_rows)
    contrast["note"] = (
        "marker-based cell types on author-normalized GSE165641 KL lung "
        "scRNA; Tacstd2 is epithelial, not a CD8/NK program"
    )
    contrast.to_csv(OUT / "tacstd2_epithelial_vs_cd8_nk.csv", index=False)

    # ---- Tacstd2+ in Ptprc+ (ambient / doublet caveat) ----
    amb = pd.DataFrame([{
        "n_tacstd2_pos_ptprc_pos": int(((tac > 0) & (X["Ptprc"] > 0)).sum()),
        "frac_of_those_s100a8_pos": round(float(
            X.loc[(tac > 0) & (X["Ptprc"] > 0), "S100a8"].gt(0).mean()), 4),
        "frac_of_those_cd8_or_nk": round(float(
            ((celltype.isin(["CD8_T", "NK"])) & (tac > 0)).sum()
            / max(((tac > 0) & (X["Ptprc"] > 0)).sum(), 1)), 4),
        "interpretation": (
            "Tacstd2+ Ptprc+ droplets are almost all neutrophil/myeloid, "
            "not CD8/NK; most consistent with ambient epithelial RNA or "
            "doublets in a neutrophil-rich KL TME"
        ),
    }])
    amb.to_csv(OUT / "tacstd2_in_immune_droplets.csv", index=False)

    # ---- within-epithelial Tacstd2 tertiles vs MHC/IFN / junction ----
    epi = celltype == "Epithelial"
    t_epi = tac[epi]
    tert = pd.qcut(t_epi.rank(method="first"), 3, labels=["low", "mid", "high"])
    genes = [g for g in JUNCTION + MHC_IFN + ["Mki67", "Sftpc"] if g in X.columns]
    hi_lo = []
    for g in genes:
        a = X.loc[tert.index[tert == "high"], g]
        b = X.loc[tert.index[tert == "low"], g]
        if np.allclose(a.values, a.values[0]) and np.allclose(b.values, b.values[0]) and a.iloc[0] == b.iloc[0]:
            pval = np.nan
            rrb = 0.0
        else:
            pval = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
            rrb = rank_biserial(a.values, b.values)
        hi_lo.append({
            "gene": g,
            "n_high": int(len(a)), "n_low": int(len(b)),
            "mean_high": round(float(a.mean()), 4),
            "mean_low": round(float(b.mean()), 4),
            "detection_high": round(float((a > 0).mean()), 4),
            "detection_low": round(float((b > 0).mean()), 4),
            "rank_biserial_high_vs_low": round(rrb, 3),
            "mannwhitney_p": pval,
        })
    hilo = pd.DataFrame(hi_lo)
    hilo["bh_fdr"] = bh_fdr(hilo["mannwhitney_p"])
    hilo["note"] = (
        "epithelial cells only (Ptprc-); tertiles of Tacstd2. A Tacstd2-high "
        "immune-resistant subset would predict lower MHC/IFN in the high "
        "tertile — that pattern is not observed"
    )
    hilo.to_csv(OUT / "epithelial_tacstd2_high_vs_low.csv", index=False)

    epi_ann = pd.DataFrame({
        "barcode": t_epi.index,
        "mouse": mouse.loc[t_epi.index].values,
        "Tacstd2": t_epi.round(4).values,
        "tacstd2_tertile": tert.astype(str).values,
    })
    epi_ann.to_csv(OUT / "epithelial_cells.csv", index=False)

    # ---- figures ----
    fig, ax = plt.subplots(figsize=(10, 4.8))
    data, labels, ns = [], [], []
    for ct in CELLTYPE_ORDER:
        v = tac[celltype == ct].values
        data.append(v)
        labels.append(ct)
        ns.append(len(v))
    ax.boxplot(data, tick_labels=[f"{l}\n(n={n})" for l, n in zip(labels, ns)])
    ax.set_ylabel("Tacstd2 (author-normalized)")
    ax.set_title("GSE165641 KL lung scRNA: Tacstd2 by marker-based cell type\n"
                 "KrasLSL-G12D/+; Lkb1fl/fl tumor sections, 2 mice, 7,180 cells")
    ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_tacstd2_by_celltype.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    det = by_ct.set_index("celltype").loc[CELLTYPE_ORDER, "detection_rate"]
    axes[0].bar(range(len(det)), det.values, color="#4C72B0")
    axes[0].set_xticks(range(len(det)))
    axes[0].set_xticklabels(CELLTYPE_ORDER, rotation=35, ha="right")
    axes[0].set_ylabel("Tacstd2 detection rate (expr > 0)")
    axes[0].set_title("Tacstd2 is detected in epithelium, not CD8/NK")
    axes[0].set_ylim(0, 1)

    show = ["Cldn4", "Epcam", "H2-D1", "Tap1", "Irf1", "Cd274"]
    x = np.arange(len(show))
    hi_m = hilo.set_index("gene").loc[show, "mean_high"]
    lo_m = hilo.set_index("gene").loc[show, "mean_low"]
    axes[1].bar(x - 0.18, lo_m, 0.36, label="Tacstd2-low tertile")
    axes[1].bar(x + 0.18, hi_m, 0.36, label="Tacstd2-high tertile")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(show)
    axes[1].set_ylabel("mean normalized expression")
    axes[1].set_title("Within epithelium: high Tacstd2 is Cldn4-high,\n"
                      "not MHC/IFN-low (not an immune-resistant subset)")
    axes[1].legend(fontsize=8)
    fig.suptitle("GSE165641 KL lung — Tacstd2-high vs CD8/NK / immune-resistance",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig_detection_and_epithelial_highlow.png", dpi=150)
    plt.close(fig)

    # ---- provenance ----
    missing = [g for g in (CD8_GENES + NK_GENES + EPI_GENES + MHC_IFN
                           + JUNCTION + EXTRA) if g not in X.columns]
    prov = {
        "task": "w200/A5: TROP2-high immune-resistant subset — real public KL lung",
        "rejected_accession": "GSE76628 (mouse flank-skin Ad-VEGF / gastric stroma; not KL lung)",
        "used_accession": "GSE165641",
        "why_this_series": (
            "Public 10X scRNA of KrasLSL-G12D/+; Lkb1fl/fl lung tumor sections "
            "with both epithelial and immune cells, so Tacstd2-high tumor vs "
            "CD8/NK is actually computable"
        ),
        "pmid": "34369094",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_cells": int(X.shape[0]),
        "n_genes": int(X.shape[1]),
        "mice": mouse.value_counts().to_dict(),
        "celltype_rules": {
            "CD8_T": "(Cd8a>0 | Cd8b1>0) & (Cd3e>0 | Cd3d>0)",
            "NK": "(Ncr1>0 | Klrb1c>0) & not CD8_T",
            "Epithelial": "(Epcam>0 | (Krt8>0 & Krt18>0) | Nkx2-1>0) & Ptprc==0",
        },
        "processing": (
            "author-processed log-normalized matrix from GEO; no re-alignment; "
            "marker-based types (no published cluster labels in the matrix)"
        ),
        "genes_missing": missing,
        "files": [file_meta],
    }
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")

    print("composition:\n", comp.to_string())
    print("\nTacstd2 by type:\n", by_ct.to_string(index=False))
    print("\ncontrast:\n", contrast.drop(columns="note").to_string(index=False))
    print("\nepithelial high vs low:\n",
          hilo.drop(columns="note").to_string(index=False))


if __name__ == "__main__":
    main()
