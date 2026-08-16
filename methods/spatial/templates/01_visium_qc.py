"""01 · Visium QC. Visium 质控。

Load 10x Visium (spaceranger 'outs' or .h5ad), compute per-spot QC metrics, and
apply the thresholds in config.yaml. Writes a filtered .h5ad plus QC plots.
读取 10x Visium 数据，计算每个 spot 的质控指标，按 config.yaml 阈值过滤，输出
过滤后的 .h5ad 与质控图。

Tools: scanpy >= 1.10, squidpy >= 1.7, anndata.  No results are fabricated.
工具：scanpy、squidpy、anndata。不伪造任何结果。

Usage:
    python 01_visium_qc.py --input <outs_dir_or_h5ad> --sample <id> [--config config.yaml]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

from _utils import load_config, resolve, ensure_dir, assert_panel


def read_visium(input_path: Path):
    """Read either a spaceranger 'outs' directory or a pre-made .h5ad."""
    if input_path.is_dir():
        adata = sc.read_visium(input_path)
    else:
        adata = sc.read_h5ad(input_path)
    adata.var_names_make_unique()
    return adata


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="spaceranger outs dir or .h5ad")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    q = cfg["qc"]["visium"]
    out = ensure_dir(resolve(cfg["paths"]["output_dir"]) / "visium" / args.sample)

    adata = read_visium(Path(args.input))
    adata.obs["sample"] = args.sample

    # Whole-transcriptome data must carry TACSTD2/CLDN4; fail early if it does not.
    assert_panel(adata.var_names, cfg, platform="visium")

    # QC metrics. 质控指标。
    adata.var["mito"] = adata.var_names.str.startswith(q["mito_prefix"])
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mito"], inplace=True, percent_top=None, log1p=False
    )

    # Report distributions BEFORE filtering so thresholds are auditable.
    # 过滤前先记录分布，使阈值可审计。
    qc_summary = adata.obs[
        ["total_counts", "n_genes_by_counts", "pct_counts_mito"]
    ].describe()
    qc_summary.to_csv(out / "qc_metrics_summary.csv")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        for ax, col in zip(
            axes, ["total_counts", "n_genes_by_counts", "pct_counts_mito"]
        ):
            adata.obs[col].hist(bins=60, ax=ax)
            ax.set_title(col)
        fig.tight_layout()
        fig.savefig(out / "qc_distributions.png", dpi=150)
        plt.close(fig)
    except Exception as exc:  # plotting is optional, never block QC on it
        print(f"[warn] skipped QC plot: {exc}")

    n0 = adata.n_obs
    keep = (
        (adata.obs["total_counts"] >= q["min_counts_per_spot"])
        & (adata.obs["n_genes_by_counts"] >= q["min_genes_per_spot"])
        & (adata.obs["pct_counts_mito"] <= q["max_pct_mito"])
    )
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=q["min_spots_per_gene"])

    # Keep raw counts available for downstream models (cell2location needs them).
    # 保留原始 counts 供下游模型使用（cell2location 需要）。
    adata.layers["counts"] = adata.X.copy()

    adata.write(out / f"{args.sample}_qc.h5ad")
    print(
        f"[01_qc] {args.sample}: kept {adata.n_obs}/{n0} spots, "
        f"{adata.n_vars} genes -> {out / (args.sample + '_qc.h5ad')}"
    )


if __name__ == "__main__":
    main()
