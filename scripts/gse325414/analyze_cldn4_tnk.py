#!/usr/bin/env python3
"""GSE325414 donor-level malignant CLDN4 vs T/NK.

Jimenez et al., Clin Cancer Res (GEO GSE325414). BD Rhapsody WTA of early-stage
NSCLC in a pulsed-electric-field treat-and-resect study. Author labels are in
the deposited cell metadata. GEO overall design:

  INXBX = index diagnostic biopsy
  RESRT = resection, treatment (ablated) area
  RESRU = resection, untreated tumor
  RESRI = resection, indeterminate

Primary unit is the donor on RESRU (untreated resected tumor). RESRT is
reported separately because ablation recruits lymphocytes. The patient is the
unit. Cell-level p-values are not computed.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "GSE325414"
OUT = ROOT / "results" / "gse325414"

GENES = ["CLDN4", "TACSTD2", "EPCAM", "KRT8", "KRT18", "KRT19"]
MIN_MALIGNANT = 50
MIN_TNK = 30
PRIMARY = "RESRU"


def read_features(path: Path) -> list[str]:
    genes = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    return genes


def stream_counts(mtx_path: Path, gene_rows: dict[str, int], n_cells: int) -> dict[str, np.ndarray]:
    """One pass over a genes x cells MTX. Keep selected-gene counts and nUMI."""
    counts = {g: np.zeros(n_cells, dtype=np.float32) for g in gene_rows}
    n_umi = np.zeros(n_cells, dtype=np.float64)
    row_to_gene = {r: g for g, r in gene_rows.items()}
    with gzip.open(mtx_path, "rt") as fh:
        first = fh.readline()
        if not first.startswith("%%MatrixMarket"):
            raise ValueError(f"not MTX: {first[:40]!r}")
        for line in fh:
            if line.startswith("%"):
                continue
            n_rows, n_cols, nnz = map(int, line.split())
            break
        else:
            raise ValueError("MTX header missing")
        if n_cols != n_cells:
            raise ValueError(f"MTX cells {n_cols} != barcodes {n_cells}")
        seen = 0
        for line in fh:
            if not line or line.startswith("%"):
                continue
            r_s, c_s, v_s = line.split()
            r = int(r_s)
            c = int(c_s) - 1
            v = float(v_s)
            n_umi[c] += v
            g = row_to_gene.get(r)
            if g is not None:
                counts[g][c] = v
            seen += 1
            if seen % 20_000_000 == 0:
                print(f"  mtx entries {seen:,}/{nnz:,}", flush=True)
        if seen != nnz:
            print(f"warning: read {seen} entries, header nnz {nnz}", flush=True)
    counts["nUMI"] = n_umi
    return counts


def log_cp10k(raw: np.ndarray, n_umi: np.ndarray) -> np.ndarray:
    scale = np.divide(1e4, n_umi, out=np.zeros_like(n_umi), where=n_umi > 0)
    return np.log1p(raw * scale).astype(np.float32)


def spearman(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 4:
        return {"n": int(x.size), "rho": None, "p": None}
    rho, p = stats.spearmanr(x, y)
    return {"n": int(x.size), "rho": float(rho), "p": float(p)}


def summarize(df: pd.DataFrame, compartment: str) -> pd.DataFrame:
    rows = []
    sub = df[df["compartment"] == compartment]
    for donor, g in sub.groupby("donor"):
        n = len(g)
        mal = g["is_tumor_epi"].to_numpy()
        tnk = g["is_tnk"].to_numpy()
        cd8 = g["is_cd8"].to_numpy()
        n_mal = int(mal.sum())
        n_tnk = int(tnk.sum())
        n_cd8 = int(cd8.sum())
        cld = g.loc[mal, "CLDN4_log"]
        rows.append(
            {
                "donor": donor,
                "compartment": compartment,
                "histology": g["histology"].iloc[0],
                "stage": g["stage"].iloc[0],
                "sex": g["sex"].iloc[0],
                "n_cells": n,
                "n_tumor_epi": n_mal,
                "n_tnk": n_tnk,
                "n_cd8": n_cd8,
                "n_nk": int(g["is_nk"].sum()),
                "tnk_frac": n_tnk / n if n else np.nan,
                "cd8_frac": n_cd8 / n if n else np.nan,
                "cldn4_mean": float(cld.mean()) if n_mal else np.nan,
                "cldn4_pct": float((g.loc[mal, "CLDN4_raw"] > 0).mean()) if n_mal else np.nan,
                "tacstd2_mean": float(g.loc[mal, "TACSTD2_log"].mean()) if n_mal else np.nan,
                "epcam_mean": float(g.loc[mal, "EPCAM_log"].mean()) if n_mal else np.nan,
                "passes_gate": n_mal >= MIN_MALIGNANT and n_tnk >= MIN_TNK,
            }
        )
    return pd.DataFrame(rows)


def leave_one_out(df: pd.DataFrame, x: str, y: str) -> pd.DataFrame:
    rows = []
    if len(df) < 5:
        return pd.DataFrame(rows)
    full = spearman(df[x], df[y])
    for donor in df["donor"]:
        sub = df[df["donor"] != donor]
        sp = spearman(sub[x], sub[y])
        rows.append(
            {
                "dropped_donor": donor,
                "x": x,
                "y": y,
                "n": sp["n"],
                "rho": sp["rho"],
                "p": sp["p"],
                "full_rho": full["rho"],
                "delta_rho": None if sp["rho"] is None or full["rho"] is None else sp["rho"] - full["rho"],
            }
        )
    return pd.DataFrame(rows)


def scatter(ax, df: pd.DataFrame, title: str) -> None:
    colors = {
        "Adenocarcinoma": "#1b4f72",
        "Mucinous adenocarcinoma": "#148f77",
        "Squamous cell carcinoma": "#b9770e",
    }
    for hist, g in df.groupby("histology"):
        ax.scatter(
            g["tnk_frac"],
            g["cldn4_mean"],
            s=46,
            c=colors.get(hist, "#555555"),
            label=f"{hist} (n={len(g)})",
            zorder=3,
        )
    sp = spearman(df["cldn4_mean"], df["tnk_frac"])
    rho_s = "NA" if sp["rho"] is None else f"{sp['rho']:.2f}"
    p_s = "NA" if sp["p"] is None else f"{sp['p']:.3g}"
    ax.set_title(f"{title}\nρ={rho_s}, p={p_s}, n={sp['n']}")
    ax.set_xlabel("T/NK fraction")
    ax.set_ylabel("Malignant CLDN4 mean log1p(CP10k)")
    ax.legend(frameon=False, fontsize=7, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = read_features(DATA / "GSE325414_features.tsv.gz")
    missing = [g for g in GENES if g not in genes]
    if missing:
        raise SystemExit(f"genes absent from WTA matrix: {missing}")
    gene_rows = {g: genes.index(g) + 1 for g in GENES}
    print("gene rows", gene_rows, flush=True)

    with gzip.open(DATA / "GSE325414_barcodes.tsv.gz", "rt") as fh:
        barcodes = [line.strip().split("\t")[0] for line in fh]
    print(f"barcodes {len(barcodes):,}", flush=True)
    print("streaming matrix", flush=True)
    counts = stream_counts(DATA / "GSE325414_matrix.mtx.gz", gene_rows, len(barcodes))

    meta = pd.read_csv(
        DATA / "GSE325414_metadata_individual_cells.txt.gz", sep="\t", low_memory=False
    )
    clin = pd.read_csv(Path(__file__).resolve().parent / "donor_clinical.tsv", sep="\t")
    if not meta["cell"].isin(barcodes).all():
        raise SystemExit("metadata cells are not a subset of barcodes")
    # barcodes and metadata are the same set; align by barcode order
    meta = meta.set_index("cell").loc[barcodes].reset_index()
    meta = meta.merge(clin, on="donor", how="left")
    if meta["histology"].isna().any():
        raise SystemExit("donor missing clinical table")

    level2 = meta["sub.pop.level2"].astype(str)
    level3 = meta["sub.pop.level3"].astype(str)
    n_umi = counts["nUMI"]
    frame = pd.DataFrame(
        {
            "cell": barcodes,
            "donor": meta["donor"].to_numpy(),
            "compartment": meta["sample.type.3"].astype(str).to_numpy(),
            "histology": meta["histology"].to_numpy(),
            "stage": meta["stage"].to_numpy(),
            "sex": meta["sex"].to_numpy(),
            "level2": level2.to_numpy(),
            "level3": level3.to_numpy(),
            "nUMI": n_umi,
            "is_tumor_epi": level2.eq("EpithelialCells_TumorCells").to_numpy(),
            "is_tnk": level2.isin(["Tcell", "NKcell"]).to_numpy(),
            "is_cd8": level3.eq("Tcell_CD8").to_numpy(),
            "is_nk": level2.eq("NKcell").to_numpy(),
        }
    )
    for g in GENES:
        frame[f"{g}_raw"] = counts[g]
        frame[f"{g}_log"] = log_cp10k(counts[g], n_umi)

    # sanity: tumor epithelium should carry epithelial genes
    epi = frame["is_tumor_epi"].to_numpy()
    other = ~epi
    sanity = {
        g: {
            "tumor_epi_mean_log": float(frame.loc[epi, f"{g}_log"].mean()),
            "other_mean_log": float(frame.loc[other, f"{g}_log"].mean()),
            "tumor_epi_pct": float((frame.loc[epi, f"{g}_raw"] > 0).mean()),
        }
        for g in ["CLDN4", "EPCAM", "KRT8", "TACSTD2"]
    }
    print("sanity", json.dumps(sanity, indent=2), flush=True)

    pieces = []
    for comp in ["RESRU", "RESRT", "INXBX"]:
        pieces.append(summarize(frame, comp))
    per = pd.concat(pieces, ignore_index=True)
    per.to_csv(OUT / "per_donor_compartment.tsv", sep="\t", index=False)

    tests = []
    for comp in ["RESRU", "RESRT", "INXBX"]:
        gated = per[(per["compartment"] == comp) & (per["passes_gate"])]
        luad = gated[gated["histology"].str.contains("denocarcinoma", case=False)]
        squamous = gated[gated["histology"].str.contains("Squamous", case=False)]
        for label, block in [
            ("all_histology", gated),
            ("adenocarcinoma", luad),
            ("squamous", squamous),
        ]:
            for xname, yname in [
                ("cldn4_mean", "tnk_frac"),
                ("cldn4_pct", "tnk_frac"),
                ("cldn4_mean", "cd8_frac"),
                ("tacstd2_mean", "tnk_frac"),
            ]:
                sp = spearman(block[xname], block[yname])
                tests.append(
                    {
                        "compartment": comp,
                        "subset": label,
                        "x": xname,
                        "y": yname,
                        "gate": f"tumor_epi>={MIN_MALIGNANT}, T/NK>={MIN_TNK}",
                        **sp,
                        "n_donors_present": int((per["compartment"] == comp).sum()),
                        "n_fail_gate": int(((per["compartment"] == comp) & ~per["passes_gate"]).sum()),
                    }
                )
    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(OUT / "stats.tsv", sep="\t", index=False)
    print(tests_df.to_string(index=False), flush=True)

    prim = per[(per["compartment"] == PRIMARY) & (per["passes_gate"])].copy()
    luad = prim[prim["histology"].str.contains("denocarcinoma", case=False)]
    loo = pd.concat(
        [
            leave_one_out(prim, "cldn4_mean", "tnk_frac").assign(subset="all_histology"),
            leave_one_out(luad, "cldn4_mean", "tnk_frac").assign(subset="adenocarcinoma"),
            leave_one_out(luad, "cldn4_pct", "tnk_frac").assign(subset="adenocarcinoma"),
        ],
        ignore_index=True,
    )
    loo.to_csv(OUT / "leave_one_out_resru.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.6), sharey=True)
    scatter(axes[0], prim, "RESRU, all histology")
    scatter(axes[1], luad, "RESRU, adenocarcinoma")
    fig.suptitle("GSE325414 untreated resected tumor: malignant CLDN4 vs T/NK", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "fig_resru_cldn4_vs_tnk.png", dpi=160)
    plt.close(fig)

    summary = {
        "series": "GSE325414",
        "primary_compartment": PRIMARY,
        "primary_definition": "resection, untreated tumor",
        "malignant": "author sub.pop.level2 EpithelialCells_TumorCells",
        "tnk": "author sub.pop.level2 Tcell or NKcell",
        "expression": "log1p(CP10k) within malignant cells",
        "gate": {"min_malignant": MIN_MALIGNANT, "min_tnk": MIN_TNK},
        "n_cells": int(len(frame)),
        "n_donors": int(frame["donor"].nunique()),
        "sanity_epithelial_markers": sanity,
        "primary_spearman_cldn4_mean_vs_tnk": tests_df[
            (tests_df.compartment == PRIMARY)
            & (tests_df.subset == "all_histology")
            & (tests_df.x == "cldn4_mean")
            & (tests_df.y == "tnk_frac")
        ].iloc[0].to_dict(),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
