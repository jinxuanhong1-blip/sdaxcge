#!/usr/bin/env python3
"""LUAD-only stability sweep for GSE325414 RESRU CLDN4 vs T/NK.

The patient is the unit. Squamous donors are excluded. This does not change
the pre-specified all-histology test and is not merged into concordant-4.

Stability is fixed before ranking:

  stable inverse = n >= 8, rho < 0, and every single-donor deletion
  still has rho < 0 and p < 0.05.

A result that is nominal only while T13 is included is called fragile to T13.
The grid is a sensitivity analysis. The smallest p in the grid is not a
discovery p-value.
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

GENES = [
    "CLDN4",
    "TACSTD2",
    "EPCAM",
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "CDH1",
    "PTPRC",
    "SFTPA1",
    "SFTPA2",
    "SFTPB",
    "SFTPC",
    "AGER",
    "SCGB1A1",
    "SCGB3A1",
    "TPPP3",
    "FOXJ1",
]
EPI_GENES = ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1"]
NORMAL_GENES = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3", "FOXJ1"]
MIN_MAL_GRID = [1, 20, 30, 50, 100, 200, 500]
MIN_TNK_GRID = [1, 10, 20, 30, 50, 100]
# Expression cutoffs for percent-positive, fixed before ranking.
LOG_CUTS = [0.0, 1.0]


def read_features(path: Path) -> list[str]:
    genes = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    return genes


def stream_counts(mtx_path: Path, gene_rows: dict[str, int], n_cells: int) -> dict[str, np.ndarray]:
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
            _n_rows, n_cols, nnz = map(int, line.split())
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
            c = int(c_s) - 1
            v = float(v_s)
            n_umi[c] += v
            g = row_to_gene.get(int(r_s))
            if g is not None:
                counts[g][c] = v
            seen += 1
            if seen % 20_000_000 == 0:
                print(f"  mtx entries {seen:,}/{nnz:,}", flush=True)
    counts["nUMI"] = n_umi
    return counts


def log_cp10k(raw: np.ndarray, n_umi: np.ndarray) -> np.ndarray:
    scale = np.divide(1e4, n_umi, out=np.zeros_like(n_umi), where=n_umi > 0)
    return np.log1p(raw * scale).astype(np.float32)


def module(frame: pd.DataFrame, genes: list[str]) -> np.ndarray:
    mats = [frame[f"{g}_log"].to_numpy() for g in genes if f"{g}_log" in frame.columns]
    if not mats:
        return np.zeros(len(frame), dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0).astype(np.float32)


def spearman(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 4:
        return {"n": int(x.size), "rho": None, "p": None}
    rho, p = stats.spearmanr(x, y)
    return {"n": int(x.size), "rho": float(rho), "p": float(p)}


def load_resru_cells() -> pd.DataFrame:
    genes = read_features(DATA / "GSE325414_features.tsv.gz")
    missing = [g for g in GENES if g not in genes]
    if missing:
        raise SystemExit(f"genes absent: {missing}")
    gene_rows = {g: genes.index(g) + 1 for g in GENES}
    with gzip.open(DATA / "GSE325414_barcodes.tsv.gz", "rt") as fh:
        barcodes = [line.strip().split("\t")[0] for line in fh]
    print(f"streaming {len(barcodes):,} cells", flush=True)
    counts = stream_counts(DATA / "GSE325414_matrix.mtx.gz", gene_rows, len(barcodes))
    meta = pd.read_csv(
        DATA / "GSE325414_metadata_individual_cells.txt.gz", sep="\t", low_memory=False
    )
    clin = pd.read_csv(Path(__file__).resolve().parent / "donor_clinical.tsv", sep="\t")
    meta = meta.set_index("cell").loc[barcodes].reset_index()
    meta = meta.merge(clin, on="donor", how="left")
    keep = meta["sample.type.3"].astype(str).eq("RESRU") & meta["histology"].str.contains(
        "denocarcinoma", case=False
    )
    frame = pd.DataFrame(
        {
            "donor": meta.loc[keep, "donor"].to_numpy(),
            "histology": meta.loc[keep, "histology"].to_numpy(),
            "author_tumor": meta.loc[keep, "sub.pop.level2"].astype(str).eq("EpithelialCells_TumorCells").to_numpy(),
            "is_tnk": meta.loc[keep, "sub.pop.level2"].astype(str).isin(["Tcell", "NKcell"]).to_numpy(),
            "is_cd8": meta.loc[keep, "sub.pop.level3"].astype(str).eq("Tcell_CD8").to_numpy(),
            "nUMI": counts["nUMI"][keep.to_numpy()],
        }
    )
    for g in GENES:
        raw = counts[g][keep.to_numpy()]
        frame[f"{g}_raw"] = raw
        frame[f"{g}_log"] = log_cp10k(raw, frame["nUMI"].to_numpy())
    frame["epi_score"] = module(frame, EPI_GENES)
    frame["normal_score"] = module(frame, NORMAL_GENES)
    print(
        f"RESRU adenocarcinoma-family cells {len(frame):,} donors {frame['donor'].nunique()}",
        flush=True,
    )
    return frame


def malignant_masks(frame: pd.DataFrame) -> dict[str, np.ndarray]:
    author = frame["author_tumor"].to_numpy()
    normal = frame["normal_score"].to_numpy()
    epi = frame["epi_score"].to_numpy()
    ptprc = frame["PTPRC_log"].to_numpy()
    # Within each donor, drop the top quartile of normal-lung score inside author tumor.
    drop_high = author.copy()
    for donor, idx in frame.groupby("donor").groups.items():
        idx = np.asarray(list(idx))
        local = idx[author[idx]]
        if local.size < 20:
            continue
        cut = float(np.quantile(normal[local], 0.75))
        drop_high[local[normal[local] > cut]] = False
    return {
        "author_tumor": author,
        "author_drop_normal_q75": drop_high,
        "author_and_normal_lt_0.25": author & (normal < 0.25),
        "marker_epi": (epi >= 1.0) & (ptprc < 0.3) & (normal < 0.3),
        "marker_and_author": (epi >= 1.0) & (ptprc < 0.3) & (normal < 0.3) & author,
    }


def per_donor(frame: pd.DataFrame, mal: np.ndarray) -> pd.DataFrame:
    rows = []
    work = frame.copy()
    work["mal"] = mal
    for donor, g in work.groupby("donor"):
        n = len(g)
        m = g["mal"].to_numpy()
        n_mal = int(m.sum())
        n_tnk = int(g["is_tnk"].sum())
        n_cd8 = int(g["is_cd8"].sum())
        logs = g.loc[m, "CLDN4_log"]
        raw = g.loc[m, "CLDN4_raw"]
        row = {
            "donor": donor,
            "histology": g["histology"].iloc[0],
            "strict_luad": g["histology"].iloc[0] == "Adenocarcinoma",
            "n_cells": n,
            "n_mal": n_mal,
            "n_tnk": n_tnk,
            "n_cd8": n_cd8,
            "tnk_frac": n_tnk / n if n else np.nan,
            "cd8_frac": n_cd8 / n if n else np.nan,
            "cldn4_mean": float(logs.mean()) if n_mal else np.nan,
            "cldn4_pct0": float((raw > 0).mean()) if n_mal else np.nan,
            "cldn4_pct1": float((logs > 1.0).mean()) if n_mal else np.nan,
        }
        rows.append(row)
    return pd.DataFrame(rows)


def eval_block(block: pd.DataFrame, x: str, y: str) -> dict:
    sp = spearman(block[x], block[y])
    out = {
        "n": sp["n"],
        "rho": sp["rho"],
        "p": sp["p"],
        "worst_loo_donor": None,
        "worst_loo_rho": None,
        "worst_loo_p": None,
        "t13_in_set": bool((block["donor"] == "T13").any()) if len(block) else False,
        "rho_without_t13": None,
        "p_without_t13": None,
        "stable_inverse": False,
        "fragile_to_t13": False,
    }
    if sp["n"] is None or sp["n"] < 5 or sp["rho"] is None:
        return out
    loos = []
    for donor in block["donor"]:
        sub = block[block["donor"] != donor]
        one = spearman(sub[x], sub[y])
        loos.append((donor, one["rho"], one["p"]))
        if donor == "T13":
            out["rho_without_t13"] = one["rho"]
            out["p_without_t13"] = one["p"]
    # Worst = largest p (least significant). Ties break toward the less negative rho.
    loos_ok = [t for t in loos if t[1] is not None]
    if loos_ok:
        worst = max(loos_ok, key=lambda t: (t[2] if t[2] is not None else 1.0, -(t[1] or 0)))
        out["worst_loo_donor"] = worst[0]
        out["worst_loo_rho"] = worst[1]
        out["worst_loo_p"] = worst[2]
        all_sig = all(r is not None and r < 0 and p is not None and p < 0.05 for _, r, p in loos_ok)
        out["stable_inverse"] = bool(sp["n"] >= 8 and sp["rho"] < 0 and sp["p"] < 0.05 and all_sig and len(loos_ok) == len(block))
    if out["t13_in_set"] and sp["p"] is not None and sp["p"] < 0.05 and sp["rho"] < 0:
        p_drop = out["p_without_t13"]
        out["fragile_to_t13"] = p_drop is None or p_drop >= 0.05 or (out["rho_without_t13"] is not None and out["rho_without_t13"] >= 0)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame = load_resru_cells()
    masks = malignant_masks(frame)
    coverage = {name: int(mask.sum()) for name, mask in masks.items()}
    print("malignant cell counts", coverage, flush=True)

    donor_tables = {name: per_donor(frame, mask) for name, mask in masks.items()}
    # Audit table at the original gate for each definition.
    audit_rows = []
    grid_rows = []
    for definition, table in donor_tables.items():
        for histology, hist_block in [
            ("luad_incl_mucinous", table),
            ("luad_strict", table[table["strict_luad"]]),
        ]:
            for min_mal in MIN_MAL_GRID:
                for min_tnk in MIN_TNK_GRID:
                    block = hist_block[(hist_block["n_mal"] >= min_mal) & (hist_block["n_tnk"] >= min_tnk)]
                    for x, y in [
                        ("cldn4_mean", "tnk_frac"),
                        ("cldn4_pct0", "tnk_frac"),
                        ("cldn4_pct1", "tnk_frac"),
                        ("cldn4_mean", "cd8_frac"),
                    ]:
                        ev = eval_block(block, x, y)
                        grid_rows.append(
                            {
                                "definition": definition,
                                "histology": histology,
                                "min_mal": min_mal,
                                "min_tnk": min_tnk,
                                "x": x,
                                "y": y,
                                **ev,
                            }
                        )
                        if min_mal == 50 and min_tnk == 30 and y == "tnk_frac":
                            audit_rows.append(grid_rows[-1])

    grid = pd.DataFrame(grid_rows)
    audit = pd.DataFrame(audit_rows)
    grid.to_csv(OUT / "luad_sweep_grid.tsv", sep="\t", index=False)
    audit.to_csv(OUT / "luad_sweep_original_gate.tsv", sep="\t", index=False)

    target = grid[(grid["histology"] == "luad_incl_mucinous") & (grid["y"] == "tnk_frac")].copy()
    nominal = target[(target["p"].notna()) & (target["p"] < 0.05) & (target["rho"] < 0)]
    stable = target[target["stable_inverse"] == True]  # noqa: E712
    fragile = nominal[nominal["fragile_to_t13"] == True]  # noqa: E712
    survives_t13 = nominal[nominal["fragile_to_t13"] == False]  # noqa: E712

    def strongest(df: pd.DataFrame) -> dict | None:
        if df.empty:
            return None
        hit = df.sort_values(["rho", "p"], ascending=[True, True]).iloc[0]
        return hit.to_dict()

    summary = {
        "compartment": "RESRU",
        "histology_primary": "adenocarcinoma including mucinous",
        "stability_rule": "n>=8, rho<0, p<0.05, and every leave-one-out rho<0 and p<0.05",
        "malignant_cell_counts": coverage,
        "n_grid_tnk": int(len(target)),
        "n_nominal_inverse_p_lt_0.05": int(len(nominal)),
        "n_stable_inverse": int(len(stable)),
        "n_nominal_fragile_to_t13": int(len(fragile)),
        "n_nominal_still_p_lt_0.05_without_t13": int(len(survives_t13)),
        "strongest_stable": strongest(stable),
        "strongest_nominal": strongest(nominal),
        "strongest_nominal_that_survives_t13": strongest(survives_t13),
    }
    (OUT / "luad_sweep_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    # Figure: original-gate LUAD points with T13 marked, and gate path with/without T13.
    author = donor_tables["author_tumor"]
    base = author[(author["n_mal"] >= 50) & (author["n_tnk"] >= 30)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6))
    ax = axes[0]
    colors = {"Adenocarcinoma": "#1b4f72", "Mucinous adenocarcinoma": "#148f77"}
    for hist, g in base.groupby("histology"):
        ax.scatter(g["tnk_frac"], g["cldn4_mean"], s=50, c=colors.get(hist, "#555"), label=hist, zorder=3)
    t13 = base[base["donor"] == "T13"]
    if len(t13):
        ax.scatter(t13["tnk_frac"], t13["cldn4_mean"], s=120, facecolors="none", edgecolors="#c0392b", linewidths=1.6, zorder=4, label="T13")
        ax.annotate("T13", (t13["tnk_frac"].iloc[0], t13["cldn4_mean"].iloc[0]), textcoords="offset points", xytext=(6, -8), color="#c0392b", fontsize=8)
    sp = spearman(base["cldn4_mean"], base["tnk_frac"])
    ax.set_title(f"Author tumor, gate 50/30\nρ={sp['rho']:.2f}, p={sp['p']:.3g}, n={sp['n']}")
    ax.set_xlabel("T/NK fraction")
    ax.set_ylabel("CLDN4 mean log1p(CP10k)")
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    path_rows = []
    for min_mal in MIN_MAL_GRID:
        for label, block0 in [("with T13", author), ("without T13", author[author["donor"] != "T13"])]:
            block = block0[(block0["n_mal"] >= min_mal) & (block0["n_tnk"] >= 30)]
            one = spearman(block["cldn4_mean"], block["tnk_frac"])
            path_rows.append({"min_mal": min_mal, "label": label, **one})
    path = pd.DataFrame(path_rows)
    path.to_csv(OUT / "luad_gate_path_mean.tsv", sep="\t", index=False)
    for label, g in path.groupby("label"):
        ax.plot(g["min_mal"], g["rho"], marker="o", label=label)
    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xlabel("Minimum malignant cells (T/NK ≥ 30)")
    ax.set_ylabel("Spearman ρ, CLDN4 mean vs T/NK")
    ax.set_title("Same author call, gate path")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.suptitle("GSE325414 RESRU LUAD: CLDN4 mean vs T/NK is fragile to T13", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "fig_luad_t13_sensitivity.png", dpi=160)
    plt.close(fig)

    print(json.dumps({k: summary[k] for k in summary if k not in ("strongest_stable", "strongest_nominal", "strongest_nominal_that_survives_t13")}, indent=2))
    print("--- strongest stable ---")
    print(summary["strongest_stable"])
    print("--- strongest nominal ---")
    print(summary["strongest_nominal"])
    print("--- strongest nominal surviving T13 ---")
    print(summary["strongest_nominal_that_survives_t13"])
    print("original gate audit:")
    cols = ["definition", "histology", "x", "n", "rho", "p", "p_without_t13", "fragile_to_t13", "stable_inverse", "worst_loo_donor", "worst_loo_p"]
    print(audit[cols].to_string(index=False))


if __name__ == "__main__":
    main()
