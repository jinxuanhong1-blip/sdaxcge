#!/usr/bin/env python3
"""Score Cldn4 vs T / exclusion on GSE261890 public spatial libraries."""
from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    DATA,
    EXCLUSION_T,
    OUT,
    T_SCORE_GENES,
    dump_json,
    gene_index,
    hex_neighbors,
    high_low_cut,
    load_panel,
    mtx_shape,
    read_features,
    spearman,
    stream_mtx_panel,
)

SPATIAL_DIR = DATA / "GSE261890"
EXTRACT = SPATIAL_DIR / "extracted"

TARS = [
    dict(
        gse="GSE261890",
        gsm="GSM8153570",
        sample="Sen.Res.Veh",
        arm="vehicle",
        tar="GSM8153570_Sen.Res.Veh.Spatial.tar.gz",
        note="Public library title Sen.Res.Veh; overall-design lists sen+veh / sen+ICI / res+veh / res+ICI but only 2 GSM deposited",
    ),
    dict(
        gse="GSE261890",
        gsm="GSM8153571",
        sample="Sen.Res.ICI",
        arm="ICI",
        tar="GSM8153571_Sen.Res.ICI.Spatial.tar.gz",
        note="Public library title Sen.Res.ICI; n=1 section",
    ),
]

# BMKMANU S1000 deposits many bin levels. L7 is the primary public cell-like bin
# (Veh ~3–8k barcodes). L1 is subcellular and too large for this leftover.
PRIMARY_LEVEL = "L7_heAuto"
SENS_LEVEL = "L5_heAuto"


def _is_level_matrix(name: str, level: str) -> bool:
    n = name.replace("\\", "/")
    if f"/{level}/" not in n and not n.endswith(f"/{level}"):
        return False
    base = Path(n).name.lower()
    return base in {
        "matrix.mtx.gz",
        "matrix.mtx",
        "features.tsv.gz",
        "features.tsv",
        "barcodes.tsv.gz",
        "barcodes.tsv",
        "barcodes_pos.tsv.gz",
        "barcodes_pos.tsv",
    }


def extract_level(tar_path: Path, dest: Path, level: str) -> dict[str, Path]:
    dest.mkdir(parents=True, exist_ok=True)
    out = {}
    with tarfile.open(tar_path, "r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile() or not _is_level_matrix(m.name, level):
                continue
            target = dest / f"{level}_{Path(m.name).name}"
            if not (target.exists() and target.stat().st_size > 0):
                src = tf.extractfile(m)
                if src is None:
                    continue
                target.write_bytes(src.read())
            key = Path(m.name).name.split(".")[0]  # matrix / features / barcodes / barcodes_pos
            if "barcodes_pos" in Path(m.name).name:
                key = "barcodes_pos"
            elif Path(m.name).name.startswith("barcodes"):
                key = "barcodes"
            elif Path(m.name).name.startswith("features"):
                key = "features"
            elif "matrix" in Path(m.name).name:
                key = "matrix"
            out[key] = target
    return out


def neighbor_mean(x, y, values, k=6):
    """Mean of k nearest neighbors (exclude self). Grid-hash for 10k bins."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    values = np.asarray(values, float)
    n = x.size
    out = np.full(n, np.nan)
    if n < k + 2:
        return out
    # cell size ~ median nearest spacing
    # hash into bins of width ~ 3 * median delta
    order = np.argsort(x)
    dx = np.diff(x[order])
    dx = dx[dx > 0]
    cell = float(np.median(dx)) * 3 if dx.size else 10.0
    cell = max(cell, 1e-6)
    buckets = {}
    ix = np.floor(x / cell).astype(int)
    iy = np.floor(y / cell).astype(int)
    for i, (a, b) in enumerate(zip(ix, iy)):
        buckets.setdefault((a, b), []).append(i)
    for i in range(n):
        cand = []
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                cand.extend(buckets.get((ix[i] + da, iy[i] + db), []))
        if len(cand) <= 1:
            continue
        cand = np.asarray(cand, int)
        cand = cand[cand != i]
        if cand.size == 0:
            continue
        d2 = (x[cand] - x[i]) ** 2 + (y[cand] - y[i]) ** 2
        take = cand[np.argsort(d2)[:k]]
        out[i] = float(np.nanmean(values[take]))
    return out


def score_mtx_section(name, mtx, features, barcodes, pos_path, wanted, level="L7"):
    genes = read_features(features)
    ng, nc, nz = mtx_shape(mtx)
    rec = {
        "sample": name,
        "fmt": "mtx",
        "mtx_shape": f"{ng}x{nc}",
        "n_genes_features": len(genes),
    }
    if ng != len(genes):
        rec["status"] = "features_mismatch"
        rec["detail"] = f"mtx {ng} vs features {len(genes)}"
        return rec, None
    idx = gene_index(genes, wanted)
    rec["genes_found"] = ",".join(sorted(idx))
    rec["Cldn4_in_matrix"] = "Cldn4" in idx
    mat, umi, n_genes, n_use, n_drop, keep = stream_mtx_panel(mtx, idx, nc, min_umi=50)
    rec["n_spots_raw"] = nc
    rec["n_spots_used"] = n_use
    rec["n_spots_dropped"] = n_drop
    rec["median_umi"] = float(np.median(umi)) if n_use else np.nan
    rec["status"] = "ok"
    logmat = {k: np.log1p(v) for k, v in mat.items()}
    cldn = logmat.get("Cldn4")
    tscore = None
    t_arrs = [logmat[g] for g in T_SCORE_GENES if g in logmat]
    if t_arrs:
        tscore = np.mean(np.vstack(t_arrs), axis=0)
    epi_arrs = [logmat[g] for g in ["Epcam", "Cdh1", "Krt8", "Krt18"] if g in logmat]
    epi = np.mean(np.vstack(epi_arrs), axis=0) if epi_arrs else None
    rec["level"] = level
    rec["Cldn4_mean"] = float(cldn.mean()) if cldn is not None else np.nan
    rec["Cldn4_pct"] = float((mat["Cldn4"] > 0).mean() * 100) if "Cldn4" in mat else np.nan
    rec["Tscore_mean"] = float(tscore.mean()) if tscore is not None else np.nan
    rec["Epi_mean"] = float(epi.mean()) if epi is not None else np.nan
    if cldn is not None and tscore is not None:
        rec["same_Cldn4_T"] = spearman(cldn, tscore)
        rec["exclusion_Cldn4_minus_T"] = float((cldn - tscore).mean())
        for how in ("median", "quartile"):
            _, meta = high_low_cut(cldn, tscore, how=how)
            rec[f"hiCldn4_loT_{how}"] = meta
    xs = ys = None
    if barcodes is not None and Path(barcodes).exists() and pos_path is not None and Path(pos_path).exists():
        bc = pd.read_csv(barcodes, header=None)[0].astype(str).to_numpy()
        pos = pd.read_csv(pos_path, sep="\t", header=None)
        pos.columns = ["barcode", "x", "y"] + [f"c{i}" for i in range(max(0, pos.shape[1] - 3))]
        pos["barcode"] = pos["barcode"].astype(str)
        lutx = dict(zip(pos["barcode"], pos["x"].astype(float)))
        luty = dict(zip(pos["barcode"], pos["y"].astype(float)))
        if len(bc) == nc:
            bc_k = bc[keep]
            xs = np.array([lutx.get(b, np.nan) for b in bc_k], float)
            ys = np.array([luty.get(b, np.nan) for b in bc_k], float)
            rec["n_with_xy"] = int(np.isfinite(xs).sum())
            if tscore is not None and cldn is not None and rec["n_with_xy"] >= 50:
                nei = neighbor_mean(xs, ys, tscore, k=6)
                rec["nei_Cldn4_vs_T"] = spearman(cldn, nei)
    spots = pd.DataFrame(
        {
            "sample": name,
            "level": level,
            "Cldn4": cldn if cldn is not None else np.nan,
            "Tscore": tscore if tscore is not None else np.nan,
            "Epi": epi if epi is not None else np.nan,
            "umi": umi,
            "x": xs if xs is not None else np.nan,
            "y": ys if ys is not None else np.nan,
        }
    )
    return rec, spots


def score_one(cfg, wanted, level=PRIMARY_LEVEL):
    tar_path = SPATIAL_DIR / cfg["tar"]
    rec = {k: cfg.get(k) for k in ["gse", "gsm", "sample", "arm", "note"]}
    rec["level"] = level
    if not tar_path.exists():
        rec["status"] = "missing_tar"
        return rec, None
    dest = EXTRACT / cfg["sample"]
    files = extract_level(tar_path, dest, level)
    rec["extracted"] = {k: v.name for k, v in files.items()}
    if "matrix" not in files or "features" not in files:
        rec["status"] = "no_matrix_in_tar"
        rec["detail"] = f"level {level} missing MTX/features"
        return rec, None
    scored, spots = score_mtx_section(
        cfg["sample"],
        files["matrix"],
        files["features"],
        files.get("barcodes"),
        files.get("barcodes_pos"),
        wanted,
        level=level,
    )
    rec.update(scored)
    return rec, spots


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wanted, _ = load_panel()
    rows = []
    frames = []
    for cfg in TARS:
        print("SPATIAL", cfg["sample"], flush=True)
        rec, spots = score_one(cfg, wanted)
        rows.append(rec)
        if spots is not None:
            frames.append(spots)
        print(" ", rec.get("status"), rec.get("fmt"), rec.get("Cldn4_mean"), rec.get("same_Cldn4_T"), flush=True)
    pd.DataFrame([{k: v for k, v in r.items() if not isinstance(v, (list, dict))} for r in rows]).to_csv(
        OUT / "spatial_sample_inventory.tsv", sep="\t", index=False
    )
    dump_json(OUT / "spatial_sample_records.json", rows)
    if frames:
        pd.concat(frames, ignore_index=True).to_csv(OUT / "spatial_spot_scores.tsv.gz", sep="\t", index=False)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
