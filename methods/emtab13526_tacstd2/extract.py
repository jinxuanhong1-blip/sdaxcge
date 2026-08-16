#!/usr/bin/env python3
"""Stream E-MTAB-13526 unfiltered 10x mtx files to a lineage/target gene panel.

Deposited matrices are Cell Ranger 3.1.0 *raw* feature-barcode matrices
(33,538 genes × 6,794,880 whitelist barcodes). Author QC (De Zuani et al.
Nat Commun 2024 methods): UMI 400–100,000, genes 180–6,000, mito < 20%.
We apply the same UMI/gene/mito gates. Barcode sequences are not required;
column index + sample id is enough for patient-level summaries.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

MIN_UMI = 400
MAX_UMI = 100_000
MIN_GENES = 180
MAX_GENES = 6_000
MAX_MITO = 0.20

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T": ["CD3D", "CD3E", "CD2"],
    "NK": ["NKG7", "GNLY", "FGFBP2"],
    "B": ["CD79A", "MS4A1"],
    "myeloid": ["LYZ", "CD68", "CD14"],
    "fibroblast": ["COL1A1", "DCN"],
    "endothelial": ["VWF", "PECAM1"],
}
NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
TARGETS = ["TACSTD2", "CLDN4"]
PANEL = sorted({g for genes in LINEAGES.values() for g in genes} | set(NORMAL_LUNG) | set(TARGETS))

CD235A_TUMOR = [
    "P4_T2",
    "P4_T3",
    "P8_T2",
    "P15_T2",
    "P16_T2",
    "P17_T2",
    "P17_T3",
    "P18_T2",
    "P19_T2",
    "P20_T2",
    "P21_T1",
    "P21_T2",
    "P22_T1",
    "P23_T1",
    "P24_T1",
]


def load_features(path: Path) -> tuple[dict[int, str], list[int], dict[str, int]]:
    symbols: list[str] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            symbols.append(parts[1] if len(parts) > 1 else parts[0])
    symbol_to_idx = {}
    for i, s in enumerate(symbols):
        symbol_to_idx.setdefault(s, i)
    panel_map = {}
    missing = []
    for g in PANEL:
        if g in symbol_to_idx:
            panel_map[symbol_to_idx[g]] = g
        else:
            missing.append(g)
    mt_idx = [i for i, s in enumerate(symbols) if s.startswith("MT-")]
    return panel_map, mt_idx, {
        "n_features": len(symbols),
        "present": sorted(panel_map.values()),
        "missing": missing,
        "n_mt_genes": len(mt_idx),
        "index": {g: symbol_to_idx[g] for g in panel_map.values()},
    }


def stream_mtx(
    mtx_path: Path,
    panel_map: dict[int, str],
    mt_idx: list[int],
) -> tuple[np.ndarray, dict[str, np.ndarray], dict]:
    mt_set = set(mt_idx)
    panel_coo: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(mtx_path, "rt") as fh:
        header = fh.readline()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"not Matrix Market: {mtx_path}")
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_genes, n_bc, n_nz = (int(x) for x in line.split())
        total = np.zeros(n_bc, dtype=np.uint32)
        n_feat = np.zeros(n_bc, dtype=np.uint16)
        mito = np.zeros(n_bc, dtype=np.uint32)
        for raw in fh:
            a, b, c = raw.split()
            gi = int(a) - 1
            bi = int(b) - 1
            cv = int(c)
            total[bi] += cv
            n_feat[bi] += 1
            if gi in mt_set:
                mito[bi] += cv
            name = panel_map.get(gi)
            if name is not None:
                panel_coo[name].append((bi, cv))

    mito_frac = np.divide(mito, total, out=np.zeros(n_bc, dtype=np.float32), where=total > 0)
    keep = np.where(
        (total >= MIN_UMI)
        & (total <= MAX_UMI)
        & (n_feat >= MIN_GENES)
        & (n_feat <= MAX_GENES)
        & (mito_frac <= MAX_MITO)
    )[0]
    pos = {int(i): k for k, i in enumerate(keep.tolist())}
    n_keep = len(keep)
    expr = {g: np.zeros(n_keep, dtype=np.float32) for g in PANEL}
    for g, pairs in panel_coo.items():
        arr = expr[g]
        for bi, cv in pairs:
            j = pos.get(bi)
            if j is not None:
                arr[j] = cv
    meta = {
        "barcode_col": keep.astype(np.int32),
        "total": total[keep].astype(np.int32),
        "n_genes": n_feat[keep].astype(np.int32),
        "mito_frac": mito_frac[keep],
    }
    qc = {
        "mtx": mtx_path.name,
        "n_genes_matrix": n_genes,
        "n_barcodes_whitelist": n_bc,
        "n_nonzero": n_nz,
        "n_pass_qc": int(n_keep),
        "n_umi_ge1": int((total >= 1).sum()),
        "n_umi_ge400": int((total >= MIN_UMI).sum()),
    }
    return expr, meta, qc


def patient_from_sample(sample: str) -> str:
    return sample.split("_")[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("data/emtab13526/raw"))
    ap.add_argument("--out", type=Path, default=Path("data/emtab13526/extracted"))
    args = ap.parse_args()
    raw: Path = args.raw
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)

    feat_path = raw / "features.tsv.gz"
    if not feat_path.exists():
        feat_path = raw / "P4_T2-features.tsv.gz"
    panel_map, mt_idx, gene_info = load_features(feat_path)
    (out / "gene_index.json").write_text(json.dumps(gene_info, indent=2) + "\n")
    print("genes present:", gene_info["present"], "missing:", gene_info["missing"], flush=True)

    sample_table = Path("data/emtab13526/sample_metadata.tsv")
    sdrf_map: dict[str, dict] = {}
    if sample_table.exists():
        sdf = pd.read_csv(sample_table, sep="\t")
        sdrf_map = {r["sample"]: r.to_dict() for _, r in sdf.iterrows()}

    all_expr: dict[str, list[np.ndarray]] = {g: [] for g in PANEL}
    rows: list[dict] = []
    qc_all = []
    for sample in CD235A_TUMOR:
        mtx = raw / f"{sample}-matrix.mtx.gz"
        if not mtx.exists():
            print(f"SKIP missing {mtx}", flush=True)
            continue
        print(f"extract {sample}", flush=True)
        expr, meta, qc = stream_mtx(mtx, panel_map, mt_idx)
        qc["sample"] = sample
        qc_all.append(qc)
        n = len(meta["barcode_col"])
        info = sdrf_map.get(sample, {})
        for i in range(n):
            rows.append(
                {
                    "sample": sample,
                    "patient": patient_from_sample(sample),
                    "tissue": "Tumor",
                    "facs": info.get("facs", "CD235a-"),
                    "disease": info.get("disease", ""),
                    "barcode_col": int(meta["barcode_col"][i]),
                    "total": int(meta["total"][i]),
                    "n_genes": int(meta["n_genes"][i]),
                    "mito_frac": float(meta["mito_frac"][i]),
                }
            )
        for g in PANEL:
            all_expr[g].append(expr[g])
        print(f"  pass QC {n} cells", flush=True)

    if not rows:
        raise SystemExit("no cells extracted")

    meta_df = pd.DataFrame(rows)
    packed = {g: np.concatenate(all_expr[g]) for g in PANEL if all_expr[g]}
    assert all(len(v) == len(meta_df) for v in packed.values())
    np.savez_compressed(out / "gene_panel.npz", **packed)
    meta_df.to_csv(out / "cell_metadata.tsv", sep="\t", index=False)
    (out / "sample_qc.json").write_text(json.dumps(qc_all, indent=2) + "\n")
    summary = {
        "n_cells": int(len(meta_df)),
        "n_patients": int(meta_df["patient"].nunique()),
        "n_lanes": int(meta_df["sample"].nunique()),
        "qc_gates": {
            "min_umi": MIN_UMI,
            "max_umi": MAX_UMI,
            "min_genes": MIN_GENES,
            "max_genes": MAX_GENES,
            "max_mito": MAX_MITO,
        },
        "per_sample": qc_all,
        "genes": gene_info,
    }
    (out / "extract_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("n_cells", "n_patients", "n_lanes")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
