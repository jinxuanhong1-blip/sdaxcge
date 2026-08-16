#!/usr/bin/env python3
"""Table-first TROP2/CLDN4 multiplex IF / IMC spatial stats.

表优先的 TROP2/CLDN4 多重 IF / IMC 空间统计。

Runs only if a cell table exists (or --synthetic-smoke for parser QA).
Bessede mIF images and cell tables are not public. This script will not
download, reconstruct, or fabricate them.

仅在细胞表存在时运行（或 --synthetic-smoke 做解析器自检）。
Bessede mIF 图像与细胞表不公开。本脚本不下载、不重建、不编造。

See ../playbook.md.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, KDTree

HERE = Path(__file__).resolve().parent
MIF_ROOT = HERE.parent
DEFAULT_DATA = MIF_ROOT / "data"
DEFAULT_OUT = MIF_ROOT / "out"

SKIP_EN = """# SKIPPED: no cell table

Bessede et al. Clin Cancer Res 2024;30:779-785 multiplex IF images and
cell-level tables are **not public** (ethics / consent; request A. Italiano
+ ethics committee). That panel is PanCK / TROP2 / CD8 / PD-L1 / DAPI —
**not CLDN4**.

No public TROP2+CLDN4 multiplex IF or IMC cell table was found to ship here
(checked 2026-08-16). Jackson-Fischer 2020 IMC and QuPath LuCa-7color are
not substitutes.

Drop a user table at `methods/mif/data/cells.csv` (or `.parquet`) and re-run.
Do not invent pixels or cells.
"""

SKIP_ZH = """# 已跳过：没有细胞表

Bessede 等 Clin Cancer Res 2024 多重免疫荧光图像与单细胞表**不公开**
（伦理 / 知情同意；需联系 A. Italiano 并经伦理委员会批准）。该面板为
PanCK / TROP2 / CD8 / PD-L1 / DAPI，**没有 CLDN4**。

未发现可随仓库分发的公开 TROP2+CLDN4 多重 IF / IMC 细胞表（检索日期
2026-08-16）。Jackson-Fischer 2020 IMC 与 QuPath LuCa-7color 不能替代。

将使用者提供的表放到 `methods/mif/data/cells.csv`（或 `.parquet`）后重跑。
不要编造像素或细胞。
"""


def _norm(name: str) -> str:
    return (
        str(name)
        .strip()
        .lower()
        .replace("µ", "u")
        .replace("μ", "u")
        .replace("-", "_")
        .replace(" ", "_")
        .replace(":", "_")
        .replace("/", "_")
    )


def _has_any(n: str, tokens: tuple[str, ...]) -> bool:
    return any(t in n for t in tokens)


def _is_summary(n: str) -> bool:
    if _has_any(n, ("std", "stdev", "max", "min", "range", "haralick")):
        return False
    return _has_any(n, ("mean", "median", "avg")) or not _has_any(
        n, ("area", "perimeter", "eccentricity", "circularity", "solidity")
    )


def map_columns(columns: list[str]) -> dict[str, str]:
    """Map raw headers → canonical names. First match wins."""
    mapping: dict[str, str] = {}
    used_raw: set[str] = set()

    def take(canonical: str, raw: str) -> None:
        if canonical not in mapping and raw not in used_raw:
            mapping[canonical] = raw
            used_raw.add(raw)

    norms = {c: _norm(c) for c in columns}

    exact = {
        "cell_id": ("cell_id", "cellid", "object_id", "objectid", "cell", "name"),
        "image_id": ("image_id", "roi_id", "slide_id", "filename", "image", "roi", "parent"),
        "sample_id": ("sample_id", "patient_id", "case_id", "donor"),
        "x": ("x", "centroid_x", "x_centroid", "x_um", "cell__centroid_x_um", "centroid_x_um"),
        "y": ("y", "centroid_y", "y_centroid", "y_um", "cell__centroid_y_um", "centroid_y_um"),
        "phenotype": ("phenotype", "classification", "class", "celltype", "cell_type"),
        "area_um2": ("area_um2", "area", "cell__area", "cell_area"),
        "nucleus_area": ("nucleus_area", "nucleus__area"),
    }
    # Prefer earlier aliases (Object ID over Name; image_id over Parent).
    # 靠前的别名优先（Object ID 优于 Name；image_id 优于 Parent）。
    for canon, aliases in exact.items():
        for alias in aliases:
            for raw, n in norms.items():
                if n == alias:
                    take(canon, raw)

    # QuPath-style "Cell: Centroid X µm"
    for raw, n in norms.items():
        if "centroid" in n and (n.endswith("_x") or n.endswith("x") or "_x_" in n):
            take("x", raw)
        if "centroid" in n and (n.endswith("_y") or n.endswith("y") or "_y_" in n):
            take("y", raw)

    def marker_slot(n: str, marker_tokens: tuple[str, ...]) -> str | None:
        if not _has_any(n, marker_tokens):
            return None
        if not _is_summary(n) and _has_any(n, ("area", "perimeter")):
            return None
        compartment = "cell"
        if _has_any(n, ("nucleus", "nuclear")):
            compartment = "nucleus"
        elif _has_any(n, ("membrane", "membr")):
            compartment = "membrane"
        elif _has_any(n, ("cytoplasm", "cyto")):
            compartment = "cytoplasm"
        return compartment

    marker_specs = {
        "trop2": ("trop2", "tacstd2"),
        "cldn4": ("cldn4", "claudin4", "claudin_4", "claudin"),
        "panck": ("panck", "pan_ck", "cytokeratin", "ae1", "pancytokeratin"),
        "cd8": ("cd8",),
        "pd_l1": ("pd_l1", "pdl1", "cd274"),
        "dapi": ("dapi", "ir193", "ir191", "dna1", "dna2"),
    }
    # Prefer specific claudin-4 over generic "claudin" if both exist — handled by token order.
    for raw, n in norms.items():
        if "claudin" in n and not _has_any(n, ("cldn4", "claudin4", "claudin_4")):
            if "claudin7" in n or "cldn7" in n or "claudin1" in n or "cldn1" in n:
                continue
        for canon, tokens in marker_specs.items():
            slot = marker_slot(n, tokens)
            if slot is None:
                continue
            if canon == "cldn4" and "claudin" in tokens:
                if not _has_any(n, ("cldn4", "claudin4", "claudin_4", "claudin")):
                    continue
                if _has_any(n, ("cldn7", "claudin7", "cldn1", "claudin1")):
                    continue
            take(canon if slot == "cell" else f"{canon}_{slot}", raw)
            if slot == "cell":
                take(canon, raw)

    return mapping


def load_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".parquet", ".pq"}:
        df = pd.read_parquet(path)
    elif path.suffix.lower() in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)
    mapping = map_columns(list(df.columns))
    rename = {raw: canon for canon, raw in mapping.items()}
    out = df.rename(columns=rename).copy()
    # Keep unmapped originals for debugging; drop exact duplicate canonicals.
    out.attrs["column_map"] = mapping
    return out


def parse_phenotype_tokens(series: pd.Series) -> pd.DataFrame:
    s = series.astype(str).str.lower()
    return pd.DataFrame(
        {
            "is_tumor": s.str.contains(r"panck|pan_ck|cytokeratin|tumor|epithelial|ae1"),
            "is_trop2": s.str.contains(r"trop2|tacstd2"),
            "is_cldn4": s.str.contains(r"cldn4|claudin.?4"),
            "is_cd8": s.str.contains(r"cd8"),
        },
        index=series.index,
    )


def threshold_positive(df: pd.DataFrame, col: str, q: float) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    x = pd.to_numeric(df[col], errors="coerce")
    flags = pd.Series(False, index=df.index)
    if "image_id" in df.columns:
        for _, idx in df.groupby("image_id", dropna=False).groups.items():
            vals = x.loc[idx]
            thr = vals.quantile(q)
            flags.loc[idx] = vals >= thr
    else:
        flags = x >= x.quantile(q)
    return flags.fillna(False)


def assign_phenotypes(df: pd.DataFrame, q: float) -> pd.DataFrame:
    out = df.copy()
    if "phenotype" in out.columns:
        tokens = parse_phenotype_tokens(out["phenotype"])
        for c in tokens.columns:
            out[c] = tokens[c]
    else:
        for c in ("is_tumor", "is_trop2", "is_cldn4", "is_cd8"):
            out[c] = False

    if not out["is_tumor"].any() and "panck" in out.columns:
        out["is_tumor"] = threshold_positive(out, "panck", q)
    if not out["is_trop2"].any():
        src = next((c for c in ("trop2_membrane", "trop2", "trop2_cytoplasm") if c in out.columns), None)
        if src:
            out["is_trop2"] = threshold_positive(out, src, q)
    if not out["is_cldn4"].any():
        src = next((c for c in ("cldn4_membrane", "cldn4", "cldn4_cytoplasm") if c in out.columns), None)
        if src:
            out["is_cldn4"] = threshold_positive(out, src, q)
    if not out["is_cd8"].any() and "cd8" in out.columns:
        out["is_cd8"] = threshold_positive(out, "cd8", q)

    # Tumor gate: if PanCK exists, TROP2/CLDN4 positivity for co-expression is among tumor.
    out["is_dp"] = out["is_tumor"] & out["is_trop2"] & out["is_cldn4"]
    if not out["is_tumor"].any():
        out["is_dp"] = out["is_trop2"] & out["is_cldn4"]
    return out


def _xy(df: pd.DataFrame, pixel_size_um: float | None) -> np.ndarray | None:
    if "x" not in df.columns or "y" not in df.columns:
        return None
    xy = df[["x", "y"]].apply(pd.to_numeric, errors="coerce")
    if xy.isna().any().any():
        return None
    arr = xy.to_numpy(dtype=float)
    if pixel_size_um is not None:
        arr = arr * float(pixel_size_um)
    return arr


def nn_distances(src_xy: np.ndarray, tgt_xy: np.ndarray) -> np.ndarray:
    if len(src_xy) == 0 or len(tgt_xy) == 0:
        return np.array([])
    tree = KDTree(tgt_xy)
    d, _ = tree.query(src_xy, k=1)
    return np.asarray(d, dtype=float)


def g_function(src_xy: np.ndarray, tgt_xy: np.ndarray, radii: list[float]) -> dict[str, float]:
    if len(src_xy) == 0 or len(tgt_xy) == 0:
        return {f"G_r{int(r)}": float("nan") for r in radii}
    d = nn_distances(src_xy, tgt_xy)
    return {f"G_r{int(r)}": float(np.mean(d <= r)) for r in radii}


def mixing_score(src_xy: np.ndarray, cd8_xy: np.ndarray, other_xy: np.ndarray, k: int) -> float:
    """Fraction of k-NN (among all other cells) that are CD8, averaged over src."""
    if len(src_xy) == 0:
        return float("nan")
    catalog = []
    labels = []
    if len(cd8_xy):
        catalog.append(cd8_xy)
        labels.append(np.ones(len(cd8_xy), dtype=int))
    if len(other_xy):
        catalog.append(other_xy)
        labels.append(np.zeros(len(other_xy), dtype=int))
    if not catalog:
        return float("nan")
    pts = np.vstack(catalog)
    lab = np.concatenate(labels)
    kk = min(k, len(pts))
    if kk < 1:
        return float("nan")
    tree = KDTree(pts)
    _, idx = tree.query(src_xy, k=kk)
    idx = np.atleast_2d(idx)
    return float(lab[idx].mean())


def ripley_k_translation(xy: np.ndarray, radii: list[float]) -> dict[str, float]:
    """Ripley K with simple translation edge correction on the axis-aligned bbox.

    Publication-grade K should use spatstat. This is a documented approximation.
    发表级 K 请用 spatstat。此处为有文档说明的近似。
    """
    out = {f"K_r{int(r)}": float("nan") for r in radii}
    n = len(xy)
    if n < 10:
        return out
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    area = max((xmax - xmin) * (ymax - ymin), 1e-9)
    tree = KDTree(xy)
    for r in radii:
        counts = tree.query_ball_point(xy, r, return_length=True)
        # exclude self
        mean_count = float(np.mean(np.asarray(counts) - 1))
        out[f"K_r{int(r)}"] = area * mean_count / max(n - 1, 1)
    return out


def image_area_mm2(xy: np.ndarray) -> float:
    if len(xy) < 3:
        return float("nan")
    try:
        hull = ConvexHull(xy)
        # ConvexHull.volume is area in 2D
        return float(hull.volume) / 1e6
    except Exception:
        xmin, ymin = xy.min(axis=0)
        xmax, ymax = xy.max(axis=0)
        return float(max((xmax - xmin) * (ymax - ymin), 0.0)) / 1e6


def summarize_image(df: pd.DataFrame, radii: list[float], k: int, n_perm: int, rng: np.random.Generator) -> dict:
    rec: dict = {
        "n_cells": int(len(df)),
        "n_tumor": int(df["is_tumor"].sum()) if "is_tumor" in df else 0,
        "n_trop2": int(df["is_trop2"].sum()) if "is_trop2" in df else 0,
        "n_cldn4": int(df["is_cldn4"].sum()) if "is_cldn4" in df else 0,
        "n_dp": int(df["is_dp"].sum()) if "is_dp" in df else 0,
        "n_cd8": int(df["is_cd8"].sum()) if "is_cd8" in df else 0,
    }
    tumor = df["is_tumor"] if "is_tumor" in df else pd.Series(False, index=df.index)
    if tumor.any():
        rec["frac_tumor_trop2"] = float(df.loc[tumor, "is_trop2"].mean())
        rec["frac_tumor_cldn4"] = float(df.loc[tumor, "is_cldn4"].mean())
        rec["frac_tumor_dp"] = float(df.loc[tumor, "is_dp"].mean())
    else:
        rec["frac_tumor_trop2"] = rec["frac_tumor_cldn4"] = rec["frac_tumor_dp"] = float("nan")

    if "trop2_nucleus" in df.columns and tumor.any():
        rec["median_trop2_nucleus_tumor"] = float(
            pd.to_numeric(df.loc[tumor, "trop2_nucleus"], errors="coerce").median()
        )
    if "trop2_membrane" in df.columns and tumor.any():
        rec["median_trop2_membrane_tumor"] = float(
            pd.to_numeric(df.loc[tumor, "trop2_membrane"], errors="coerce").median()
        )
    if "cldn4_membrane" in df.columns and "cldn4_cytoplasm" in df.columns and tumor.any():
        mem = pd.to_numeric(df.loc[tumor, "cldn4_membrane"], errors="coerce")
        cyto = pd.to_numeric(df.loc[tumor, "cldn4_cytoplasm"], errors="coerce")
        rec["median_cldn4_membrane_ratio_tumor"] = float((mem / cyto.clip(lower=1e-9)).median())

    xy = _xy(df, None)
    if xy is None:
        rec["spatial"] = "skipped_no_xy"
        return rec

    area = image_area_mm2(xy)
    rec["area_mm2_convexhull"] = area
    if area and not math.isnan(area) and area > 0:
        rec["density_tumor_per_mm2"] = rec["n_tumor"] / area
        rec["density_cd8_per_mm2"] = rec["n_cd8"] / area
        rec["density_dp_per_mm2"] = rec["n_dp"] / area

    cd8_xy = xy[df["is_cd8"].to_numpy()] if rec["n_cd8"] else np.empty((0, 2))
    dp_xy = xy[df["is_dp"].to_numpy()] if rec["n_dp"] else np.empty((0, 2))
    tumor_xy = xy[tumor.to_numpy()] if rec["n_tumor"] else np.empty((0, 2))
    non_cd8_xy = xy[~df["is_cd8"].to_numpy()]

    d_cd8_dp = nn_distances(cd8_xy, dp_xy)
    d_cd8_tumor = nn_distances(cd8_xy, tumor_xy)
    rec["median_nn_cd8_to_dp_um"] = float(np.median(d_cd8_dp)) if len(d_cd8_dp) else float("nan")
    rec["median_nn_cd8_to_tumor_um"] = float(np.median(d_cd8_tumor)) if len(d_cd8_tumor) else float("nan")
    rec.update(g_function(cd8_xy, dp_xy, radii))
    rec["mixing_dp_kNN_frac_cd8"] = mixing_score(dp_xy, cd8_xy, non_cd8_xy, k)
    rec.update({f"dp_{kk}": vv for kk, vv in ripley_k_translation(dp_xy, radii).items()})

    # Label permutation null for G at 20 µm (CD8 labels shuffled among non-tumor if possible).
    if rec["n_cd8"] and rec["n_dp"] and n_perm > 0:
        pool = np.flatnonzero((~tumor).to_numpy()) if (~tumor).sum() >= rec["n_cd8"] else np.arange(len(df))
        g_null = []
        n_cd8 = rec["n_cd8"]
        for _ in range(n_perm):
            pick = rng.choice(pool, size=n_cd8, replace=False)
            g_null.append(g_function(xy[pick], dp_xy, [20.0])["G_r20"])
        rec["G_r20_null_mean"] = float(np.mean(g_null))
        rec["G_r20_null_p95"] = float(np.quantile(g_null, 0.95))
        rec["G_r20_obs_gt_null95"] = bool(rec.get("G_r20", 0) > rec["G_r20_null_p95"])
    return rec


def write_skip(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    text = SKIP_EN + "\n---\n\n" + SKIP_ZH
    (out_dir / "SKIPPED_NO_TABLE.md").write_text(text, encoding="utf-8")
    print(text)


def synthetic_smoke_frame(n: int = 80, seed: int = 0) -> pd.DataFrame:
    """Random points for parser/KD-tree QA. Not biology. 随机点，不是生物学。"""
    rng = np.random.default_rng(seed)
    x = rng.uniform(0, 400, n)
    y = rng.uniform(0, 400, n)
    panck = rng.normal(2.0, 0.8, n)
    trop2 = rng.normal(1.0, 0.6, n)
    cldn4 = rng.normal(1.0, 0.6, n)
    cd8 = rng.normal(0.2, 0.3, n)
    return pd.DataFrame(
        {
            "cell_id": [f"s{i}" for i in range(n)],
            "image_id": "synthetic_smoke",
            "x": x,
            "y": y,
            "panck": panck,
            "trop2": trop2,
            "cldn4": cldn4,
            "cd8": cd8,
            "area_um2": rng.uniform(40, 120, n),
        }
    )


def find_default_table() -> Path | None:
    for name in ("cells.csv", "cells.tsv", "cells.parquet", "cells.txt"):
        p = DEFAULT_DATA / name
        if p.is_file() and p.stat().st_size > 0:
            # header-only schema copy does not live here; require >1 line for csv
            if p.suffix.lower() in {".csv", ".tsv", ".txt"}:
                nlines = sum(1 for _ in p.open(encoding="utf-8", errors="ignore"))
                if nlines < 2:
                    continue
            return p
    return None


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--table", type=Path, default=None, help="Cell table CSV/TSV/parquet")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--q", type=float, default=0.75, help="Per-image positive quantile if no phenotype column")
    p.add_argument("--pixel-size-um", type=float, default=None, help="Multiply x,y if they are pixels")
    p.add_argument("--k", type=int, default=10, help="Neighbors for mixing score")
    p.add_argument("--n-perm", type=int, default=50, help="Permutations for G-function null")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument(
        "--synthetic-smoke",
        action="store_true",
        help="Run parser/KD-tree on random points. NOT a biological result.",
    )
    args = p.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    radii = [10.0, 20.0, 30.0, 50.0]

    if args.synthetic_smoke:
        print(
            "SYNTHETIC SMOKE ONLY — random geometry, not TROP2/CLDN4 biology, not Bessede.\n"
            "仅合成自检 — 随机几何，不是 TROP2/CLDN4 生物学，不是 Bessede。"
        )
        df = synthetic_smoke_frame(seed=args.seed)
        out_dir = args.out / "synthetic_smoke"
        source = "synthetic_smoke"
    else:
        table = args.table or find_default_table()
        if table is None:
            write_skip(args.out)
            return 0
        df = load_table(table)
        out_dir = args.out
        source = str(table)

    column_map = dict(getattr(df, "attrs", {}).get("column_map", {}))

    if args.pixel_size_um is not None and "x" in df.columns:
        df["x"] = pd.to_numeric(df["x"], errors="coerce") * args.pixel_size_um
        df["y"] = pd.to_numeric(df["y"], errors="coerce") * args.pixel_size_um

    if "cell_id" not in df.columns:
        df["cell_id"] = [f"row{i}" for i in range(len(df))]
    if "image_id" not in df.columns:
        df["image_id"] = "image0"

    df = assign_phenotypes(df, q=args.q)

    rows = []
    for image_id, sub in df.groupby("image_id", dropna=False):
        rec = summarize_image(sub, radii=radii, k=args.k, n_perm=args.n_perm, rng=rng)
        rec["image_id"] = image_id
        rows.append(rec)

    out_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows)
    summary_path = out_dir / "per_image_spatial_summary.csv"
    summary.to_csv(summary_path, index=False)

    meta = {
        "source": source,
        "n_rows": int(len(df)),
        "n_images": int(df["image_id"].nunique()),
        "column_map": column_map,
        "positive_quantile_if_thresholded": args.q,
        "radii_um": radii,
        "mixing_k": args.k,
        "n_perm": args.n_perm,
        "synthetic_smoke": bool(args.synthetic_smoke),
        "note_en": (
            "Methods-only spatial summaries from a cell table. "
            "Not a Bessede re-analysis. No public TROP2/CLDN4 images were used."
        ),
        "note_zh": "仅方法：由细胞表得到的空间摘要。不是 Bessede 再分析。未使用公开 TROP2/CLDN4 图像。",
    }
    (out_dir / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Wrote {summary_path}")
    print(f"Images: {meta['n_images']}  cells: {meta['n_rows']}")
    if args.synthetic_smoke:
        print("Remember: synthetic_smoke is QA, not a result. 合成自检不是结果。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
