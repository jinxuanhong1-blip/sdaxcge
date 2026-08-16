#!/usr/bin/env python3
"""Visium leftover: GSE263196 (SCLC) and GSE273378 (stage I LUAD)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    B_GENES,
    BROAD_EPI,
    EPI,
    T_GENES,
    dump_json,
    hex_neighbors,
    log_cp10k,
    partial_spearman,
    present,
    spearman,
    wilcoxon_signed,
)

OUT = ROOT / "results" / "leftover_spatial"
MIN_GENES = 200
MIN_NEI = 3


def read_positions(path):
    df = pd.read_csv(path, header=None)
    if df.shape[1] < 6:
        raise ValueError(f"unexpected positions shape {df.shape} {path}")
    # drop header row if present
    if str(df.iloc[0, 0]).lower() in {"barcode", "barcodes"}:
        df = df.iloc[1:].reset_index(drop=True)
    df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][: df.shape[1]]
    df["barcode"] = df["barcode"].astype(str)
    df["in_tissue"] = pd.to_numeric(df["in_tissue"], errors="coerce").fillna(0).astype(int)
    df["array_row"] = pd.to_numeric(df["array_row"], errors="coerce")
    df["array_col"] = pd.to_numeric(df["array_col"], errors="coerce")
    return df


def load_section(prefix_barcodes, prefix_features, prefix_mtx, prefix_pos):
    barcodes = pd.read_csv(prefix_barcodes, header=None)[0].astype(str)
    feat = pd.read_csv(prefix_features, sep="\t", header=None)
    genes = feat[1].astype(str) if feat.shape[1] > 1 else feat[0].astype(str)
    mtx = io.mmread(prefix_mtx).tocsr()
    if mtx.shape[0] == len(barcodes) and mtx.shape[1] == len(genes):
        pass
    elif mtx.shape[1] == len(barcodes) and mtx.shape[0] == len(genes):
        mtx = mtx.T.tocsr()
    else:
        raise ValueError(f"shape mismatch mtx={mtx.shape} bc={len(barcodes)} genes={len(genes)}")
    pos = read_positions(prefix_pos)
    pos = pos.set_index("barcode")
    return barcodes, genes, mtx, pos


def gene_index(genes):
    idx = {}
    for i, g in enumerate(genes):
        idx.setdefault(str(g).upper(), i)
    return idx


def score_from_log(logx, gidx, names):
    ii = [gidx[g] for g in names if g in gidx]
    if not ii:
        return np.full(logx.shape[0], np.nan)
    return logx[:, ii].mean(axis=1)


def analyze_section(name, barcodes, genes, mtx, pos):
    gidx = gene_index(genes)
    n_genes = np.asarray((mtx > 0).sum(axis=1)).ravel()
    in_t = pos.reindex(barcodes)["in_tissue"].fillna(0).astype(int).to_numpy()
    keep = (in_t == 1) & (n_genes >= MIN_GENES)
    if keep.sum() < 50:
        return None
    mtx_k = mtx[keep]
    bc_k = barcodes[keep].to_numpy()
    pos_k = pos.reindex(bc_k)
    logx = log_cp10k(mtx_k)
    cl4 = score_from_log(logx, gidx, ["CLDN4"])
    tac = score_from_log(logx, gidx, ["TACSTD2"])
    epi = score_from_log(logx, gidx, present(EPI, gidx))
    broad = score_from_log(logx, gidx, present(BROAD_EPI, gidx))
    tsc = score_from_log(logx, gidx, present(T_GENES, gidx))
    bsc = score_from_log(logx, gidx, present(B_GENES, gidx))
    tb = np.nanmean(np.vstack([tsc, bsc]), axis=0)

    same = {
        "section": name,
        "n_spots": int(keep.sum()),
        "genes_CLDN4": "CLDN4" in gidx,
        "genes_TACSTD2": "TACSTD2" in gidx,
        "same_CLDN4_TB": spearman(cl4, tb),
        "same_TACSTD2_TB": spearman(tac, tb),
        "same_epi_TB": spearman(epi, tb),
        "same_CLDN4_T": spearman(cl4, tsc),
        "same_CLDN4_B": spearman(cl4, bsc),
        "same_TACSTD2_T": spearman(tac, tsc),
        "same_TACSTD2_B": spearman(tac, bsc),
        "same_epi_T": spearman(epi, tsc),
        "same_broad_TB": spearman(broad, tb),
    }

    coord = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(pos_k["array_row"], pos_k["array_col"])) if pd.notna(r) and pd.notna(c)}
    rings = []
    for k in (1, 2, 3):
        nei_tb = np.full(len(bc_k), np.nan)
        nnei = np.zeros(len(bc_k), dtype=int)
        for (r, c), i in coord.items():
            acc = []
            if k == 1:
                nbs = hex_neighbors(r, c)
            else:
                # graph distance exactly k
                seen = {(r, c)}
                frontier = {(r, c)}
                for _ in range(k):
                    nxt = set()
                    for rr, cc in frontier:
                        for nb in hex_neighbors(rr, cc):
                            if nb not in seen:
                                nxt.add(nb)
                                seen.add(nb)
                    frontier = nxt
                nbs = list(frontier)
            for nb in nbs:
                j = coord.get(nb)
                if j is not None:
                    acc.append(tb[j])
            nnei[i] = len(acc)
            if len(acc) >= MIN_NEI:
                nei_tb[i] = float(np.nanmean(acc))
        ok = np.isfinite(nei_tb)
        rings.append(
            {
                "section": name,
                "ring": k,
                "n": int(ok.sum()),
                "CLDN4_vs_neiTB": spearman(cl4[ok], nei_tb[ok]),
                "TACSTD2_vs_neiTB": spearman(tac[ok], nei_tb[ok]),
                "epi_vs_neiTB": spearman(epi[ok], nei_tb[ok]),
                "partial_CLDN4_vs_neiTB_ctrl_broad": partial_spearman(cl4[ok], nei_tb[ok], broad[ok]),
                "partial_TACSTD2_vs_neiTB_ctrl_broad": partial_spearman(tac[ok], nei_tb[ok], broad[ok]),
            }
        )
    return {"same": same, "rings": rings}


def run_gse263196():
    d = ROOT / "data" / "GSE263196"
    sections = [
        ("SCLC3", "GSM8187469"),
        ("SCLC4", "GSM8187470"),
        ("SCLC8", "GSM8187471"),
        ("SCLC9", "GSM8187472"),
        ("SCLC12", "GSM8187473"),
    ]
    same_rows, ring_rows = [], []
    for name, gsm in sections:
        print("GSE263196", name, flush=True)
        bc, genes, mtx, pos = load_section(
            d / f"{gsm}_{name}_barcodes.tsv.gz",
            d / f"{gsm}_{name}_features.tsv.gz",
            d / f"{gsm}_{name}_matrix.mtx.gz",
            d / f"{gsm}_{name}_tissue_positions_list.csv.gz",
        )
        res = analyze_section(name, bc, genes, mtx, pos)
        if res is None:
            continue
        same_rows.append(res["same"])
        ring_rows.extend(res["rings"])
    return same_rows, ring_rows


def run_gse273378():
    d = ROOT / "data" / "GSE273378"
    files = sorted(d.glob("*_matrix.mtx.gz"))
    same_rows, ring_rows = [], []
    for mtx_p in files:
        stem = mtx_p.name.replace("_matrix.mtx.gz", "")
        print("GSE273378", stem, flush=True)
        bc, genes, mtx, pos = load_section(
            d / f"{stem}_barcodes.tsv.gz",
            d / f"{stem}_features.tsv.gz",
            mtx_p,
            d / f"{stem}_tissue_positions_list.csv.gz",
        )
        res = analyze_section(stem.split("_", 1)[-1] if "_" in stem else stem, bc, genes, mtx, pos)
        if res is None:
            continue
        # pathology overlap if present
        path_p = d / f"{stem}_pathology.csv.gz"
        if path_p.exists():
            path = pd.read_csv(path_p)
            res["same"]["pathology_labels"] = int(path.iloc[:, 1].notna().sum()) if path.shape[1] > 1 else 0
        same_rows.append(res["same"])
        ring_rows.extend(res["rings"])
    return same_rows, ring_rows


def flatten_same(rows, series):
    out = []
    for r in rows:
        rec = {"series": series, "section": r["section"], "n_spots": r["n_spots"]}
        for k, v in r.items():
            if isinstance(v, dict) and "rho" in v:
                rec[f"{k}_n"] = v["n"]
                rec[f"{k}_rho"] = v["rho"]
                rec[f"{k}_p"] = v["p"]
        out.append(rec)
    return pd.DataFrame(out)


def flatten_rings(rows, series):
    out = []
    for r in rows:
        rec = {"series": series, "section": r["section"], "ring": r["ring"], "n": r["n"]}
        for k, v in r.items():
            if isinstance(v, dict) and "rho" in v:
                rec[f"{k}_n"] = v["n"]
                rec[f"{k}_rho"] = v["rho"]
                rec[f"{k}_p"] = v["p"]
        out.append(rec)
    return pd.DataFrame(out)


def summarize(same_df, ring_df, series):
    s = same_df[same_df.series == series]
    r = ring_df[ring_df.series == series]
    summary = {"series": series, "n_sections": int(s.section.nunique())}
    for col in [
        "same_CLDN4_TB_rho",
        "same_TACSTD2_TB_rho",
        "same_epi_TB_rho",
        "same_CLDN4_T_rho",
        "same_TACSTD2_T_rho",
    ]:
        if col in s:
            w = wilcoxon_signed(s[col].to_numpy())
            summary[col.replace("_rho", "")] = w
    for k in (1, 2, 3):
        sub = r[r.ring == k]
        for col in [
            "CLDN4_vs_neiTB_rho",
            "TACSTD2_vs_neiTB_rho",
            "partial_CLDN4_vs_neiTB_ctrl_broad_rho",
            "partial_TACSTD2_vs_neiTB_ctrl_broad_rho",
        ]:
            if col in sub:
                summary[f"ring{k}_{col.replace('_rho','')}"] = wilcoxon_signed(sub[col].to_numpy())
    return summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    s1, r1 = run_gse263196()
    s2, r2 = run_gse273378()
    same = pd.concat(
        [flatten_same(s1, "GSE263196"), flatten_same(s2, "GSE273378")],
        ignore_index=True,
    )
    rings = pd.concat(
        [flatten_rings(r1, "GSE263196"), flatten_rings(r2, "GSE273378")],
        ignore_index=True,
    )
    same.to_csv(OUT / "tables" / "visium_samespot.csv", index=False)
    rings.to_csv(OUT / "tables" / "visium_rings.csv", index=False)
    summary = {
        "GSE263196": summarize(same, rings, "GSE263196"),
        "GSE273378": summarize(same, rings, "GSE273378"),
    }
    dump_json(OUT / "tables" / "visium_summary.json", summary)
    print(json_preview(summary))


def json_preview(summary):
    import json

    return json.dumps(summary, indent=2, default=str)


if __name__ == "__main__":
    main()
