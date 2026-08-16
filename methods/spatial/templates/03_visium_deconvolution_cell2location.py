"""03 · Visium deconvolution with cell2location. Visium 去卷积（cell2location）。

Two stages, matching the E-MTAB-13530 study design (Visium deconvolved with a
matched scRNA-seq reference via cell2location):
两阶段，与 E-MTAB-13530 研究一致（用匹配的 scRNA-seq 参考经 cell2location 去卷积）：

  A) estimate reference cell-type signatures from an annotated scRNA-seq .h5ad
     从带注释的 scRNA-seq 参考估计细胞类型特征
  B) map cell-type abundances onto Visium spots (raw counts required)
     将细胞类型丰度映射到 Visium spot（需原始 counts）

Tools: cell2location 0.1.4/0.1.5 (2024-2025), scvi-tools, scanpy. GPU strongly
recommended. No results are fabricated; if inputs are absent the script exits.
工具：cell2location、scvi-tools、scanpy。强烈建议使用 GPU。缺少输入则退出，不伪造。

Usage:
    # Stage A (once per reference):
    python 03_visium_deconvolution_cell2location.py ref \
        --reference <annotated_scrna.h5ad> --label-key cell_type
    # Stage B (per Visium sample):
    python 03_visium_deconvolution_cell2location.py map \
        --input <sample_qc.h5ad> --sample <id>
"""
from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

from _utils import load_config, resolve, ensure_dir


def stage_ref(args, cfg) -> None:
    import cell2location
    from cell2location.models import RegressionModel

    out = ensure_dir(resolve(cfg["paths"]["output_dir"]) / "cell2location" / "ref")
    adata_ref = sc.read_h5ad(args.reference)

    # cell2location trains on raw counts. 使用原始 counts 训练。
    from cell2location.utils.filtering import filter_genes

    selected = filter_genes(
        adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    adata_ref = adata_ref[:, selected].copy()

    RegressionModel.setup_anndata(
        adata=adata_ref,
        batch_key=args.batch_key if args.batch_key else None,
        labels_key=args.label_key,
    )
    mod = RegressionModel(adata_ref)
    mod.train(max_epochs=cfg["deconvolution"]["visium"]["cell2location"]["max_epochs_ref"])
    adata_ref = mod.export_posterior(adata_ref)
    mod.save(str(out), overwrite=True)
    adata_ref.write(out / "reference_signatures.h5ad")
    print(f"[03_ref] reference signatures -> {out}")


def stage_map(args, cfg) -> None:
    import cell2location

    c2l = cfg["deconvolution"]["visium"]["cell2location"]
    ref_dir = resolve(cfg["paths"]["output_dir"]) / "cell2location" / "ref"
    adata_ref = sc.read_h5ad(ref_dir / "reference_signatures.h5ad")

    # Per-cell-type reference expression from the fitted regression model.
    if "means_per_cluster_mu_fg" in adata_ref.varm:
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{c}" for c in adata_ref.uns["mod"]["factor_names"]]
        ].copy()
    else:
        inf_aver = adata_ref.var[
            [f"means_per_cluster_mu_fg_{c}" for c in adata_ref.uns["mod"]["factor_names"]]
        ].copy()
    inf_aver.columns = adata_ref.uns["mod"]["factor_names"]

    adata_vis = sc.read_h5ad(args.input)
    adata_vis.X = adata_vis.layers["counts"].copy()  # raw counts required
    shared = [g for g in inf_aver.index if g in adata_vis.var_names]
    adata_vis = adata_vis[:, shared].copy()
    inf_aver = inf_aver.loc[shared]

    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=c2l["n_cells_per_location"],
        detection_alpha=c2l["detection_alpha"],
    )
    mod.train(max_epochs=c2l["max_epochs_map"], batch_size=None, train_size=1)
    adata_vis = mod.export_posterior(adata_vis)

    # q05 posterior of abundances is the recommended point estimate.
    abund = adata_vis.obsm["q05_cell_abundance_w_sf"]
    abund.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in abund.columns]
    adata_vis.obsm["cell_abundance"] = abund
    # Dominant cell type per spot -> obs['cell_type'] for the niche template.
    adata_vis.obs["cell_type"] = abund.idxmax(axis=1).astype("category")

    out = ensure_dir(resolve(cfg["paths"]["output_dir"]) / "visium" / args.sample)
    adata_vis.write(out / f"{args.sample}_c2l.h5ad")
    print(f"[03_map] {args.sample}: abundances -> {out / (args.sample + '_c2l.h5ad')}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_ref = sub.add_parser("ref", help="estimate reference signatures")
    p_ref.add_argument("--reference", required=True)
    p_ref.add_argument("--label-key", required=True)
    p_ref.add_argument("--batch-key", default=None)
    p_ref.add_argument("--config", default=None)

    p_map = sub.add_parser("map", help="map abundances onto a Visium sample")
    p_map.add_argument("--input", required=True)
    p_map.add_argument("--sample", required=True)
    p_map.add_argument("--config", default=None)

    args = ap.parse_args()
    cfg = load_config(args.config)
    if args.cmd == "ref":
        stage_ref(args, cfg)
    else:
        stage_map(args, cfg)


if __name__ == "__main__":
    main()
