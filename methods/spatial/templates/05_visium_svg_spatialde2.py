"""05 · Spatially variable genes (SVG). 空间可变基因。

Primary: SpatialDE2 (variance-component test + optional tissue segmentation).
Fallback: squidpy Moran's I autocorrelation (fast, dependency-light).
主：SpatialDE2（方差分量检验 + 可选组织分区）。备选：squidpy Moran's I。

Tools: SpatialDE2 (PMBio/SpatialDE), squidpy, scanpy. No results are fabricated;
the test statistics come directly from the tools.
工具：SpatialDE2、squidpy、scanpy。不伪造结果；统计量直接来自工具。

Usage:
    python 05_visium_svg_spatialde2.py --input <sample_norm.h5ad> --sample <id>
        [--method spatialde2|moran]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc

from _utils import load_config, resolve, ensure_dir


def run_spatialde2(adata, cfg, out: Path, sample: str) -> bool:
    try:
        import SpatialDE
    except Exception as exc:
        print(f"[05_svg] SpatialDE2 unavailable ({exc}); use --method moran instead.")
        return False

    # SpatialDE2 operates on counts + coordinates. 使用 counts 与坐标。
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    svg_full, _ = SpatialDE.test(adata, omnibus=True)
    svg_full = svg_full.sort_values("padj" if "padj" in svg_full else "pval")
    svg_full.to_csv(out / f"{sample}_spatialde2_svg.csv", index=False)

    alpha = cfg["svg"]["fdr_alpha"]
    padj_col = "padj" if "padj" in svg_full else "pval"
    n_sig = int((svg_full[padj_col] < alpha).sum())
    print(f"[05_svg] SpatialDE2: {n_sig} genes with {padj_col} < {alpha}")

    # Optional tissue segmentation into expression-based regions.
    try:
        seg = SpatialDE.tissue_segmentation(adata)
        adata.obs["spatialde2_region"] = np.asarray(seg.segments).astype(str)
    except Exception as exc:
        print(f"[05_svg] tissue_segmentation skipped: {exc}")
    return True


def run_moran(adata, cfg, out: Path, sample: str) -> None:
    import squidpy as sq

    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(
            adata, coord_type=cfg["niche"]["coord_type"], n_neighs=cfg["niche"]["n_neighs"]
        )
    if adata.var.get("highly_variable") is not None:
        genes = adata.var_names[adata.var["highly_variable"]].tolist()
    else:
        genes = None
    sq.gr.spatial_autocorr(adata, mode="moran", genes=genes, seed=cfg["project"]["random_seed"])
    res = adata.uns["moranI"].sort_values("pval_norm")
    res.to_csv(out / f"{sample}_moranI_svg.csv")
    n_sig = int((res["pval_norm_fdr_bh"] < cfg["svg"]["fdr_alpha"]).sum())
    print(f"[05_svg] Moran's I: {n_sig} genes with FDR < {cfg['svg']['fdr_alpha']}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--method", default=None, choices=["spatialde2", "moran"])
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    out = ensure_dir(resolve(cfg["paths"]["output_dir"]) / "visium" / args.sample)
    adata = sc.read_h5ad(args.input)

    method = args.method or cfg["svg"]["tool"]
    if method == "spatialde2":
        ok = run_spatialde2(adata, cfg, out, args.sample)
        if not ok:
            run_moran(adata, cfg, out, args.sample)
    else:
        run_moran(adata, cfg, out, args.sample)

    adata.write(out / f"{args.sample}_svg.h5ad")


if __name__ == "__main__":
    main()
