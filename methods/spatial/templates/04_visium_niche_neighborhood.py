"""04 · Niche / neighbourhood of TACSTD2-CLDN4 tumour vs T/B/TLS spots.
04 · TACSTD2-CLDN4 肿瘤 spot 相对 T/B/TLS spot 的生态位/邻域分析。

This is the central question of the playbook. Steps:
本手册的核心问题。步骤：

  1. Label each spot as one of {tumor_TACSTD2_CLDN4, T_cell, B_cell, TLS, other}
     using marker-set scores (from template 02) and TACSTD2/CLDN4 positivity.
     用标记评分与 TACSTD2/CLDN4 阳性，将每个 spot 标注为上述类别之一。
  2. Build the Visium hex-grid spatial graph (squidpy). 构建 Visium 六边形空间图。
  3. Neighbourhood enrichment: are tumour spots preferentially adjacent to TLS/B/T,
     or excluded from them? 邻域富集：肿瘤 spot 是否偏好与 TLS/B/T 相邻或被排斥？
  4. Co-occurrence vs distance, and a compositional niche clustering.
     共现-距离曲线，以及基于邻域组成的生态位聚类。

Tools: squidpy >= 1.7, scanpy. Works on marker-score labels or on
deconvolution argmax (obs['cell_type']) if present. No results are fabricated.
工具：squidpy、scanpy。可用标记评分标签，或（若存在）去卷积 argmax。绝不伪造结果。

Usage:
    python 04_visium_niche_neighborhood.py --input <sample_norm_or_c2l.h5ad> --sample <id>
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import squidpy as sq

from _utils import load_config, resolve, ensure_dir, assert_panel, genes_in


def label_spots(adata, cfg) -> None:
    """Assign a categorical niche label to every spot. 为每个 spot 赋类别标签。

    Priority: a spot is 'tumor_TACSTD2_CLDN4' only if BOTH TACSTD2 and CLDN4 are
    detected AND the tumour score is the winning lineage score. Otherwise it is
    assigned to whichever of T/B/TLS score is highest above zero, else 'other'.
    优先级：仅当 TACSTD2 与 CLDN4 均检出且肿瘤评分为最高谱系评分时标为肿瘤；
    否则取 T/B/TLS 中最高（且>0）者，均不满足则为 other。
    """
    import scipy.sparse as sp

    def expr(gene):
        X = adata[:, gene].X
        return np.asarray(X.todense()).ravel() if sp.issparse(X) else np.asarray(X).ravel()

    lineage_scores = {
        "tumor_TACSTD2_CLDN4": "score_tumor_epithelial",
        "T_cell": "score_t_cell",
        "B_cell": "score_b_cell",
        "TLS": "score_tls",
    }
    have = {k: v for k, v in lineage_scores.items() if v in adata.obs}
    S = adata.obs[list(have.values())].to_numpy()
    winners = np.array(list(have.keys()))[S.argmax(axis=1)]
    winners[S.max(axis=1) <= 0] = "other"

    # Enforce TACSTD2 AND CLDN4 co-positivity for the tumour label.
    tacstd2 = expr("TACSTD2") > 0
    cldn4 = expr("CLDN4") > 0
    is_tumor_call = winners == "tumor_TACSTD2_CLDN4"
    downgrade = is_tumor_call & ~(tacstd2 & cldn4)
    winners[downgrade] = "other"

    adata.obs["niche_label"] = winners
    adata.obs["niche_label"] = adata.obs["niche_label"].astype("category")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--config", default=None)
    ap.add_argument(
        "--cluster-key", default=None,
        help="obs key to use as categories; default from config.niche.cluster_key "
             "if present, else the marker-derived 'niche_label'",
    )
    args = ap.parse_args()

    cfg = load_config(args.config)
    nz = cfg["niche"]
    out = ensure_dir(resolve(cfg["paths"]["output_dir"]) / "visium" / args.sample)

    adata = sc.read_h5ad(args.input)
    assert_panel(adata.var_names, cfg, platform="visium")

    # Ensure marker scores exist (recompute if input skipped template 02).
    for label, genes in cfg["markers"].items():
        key = f"score_{label}"
        if key not in adata.obs:
            present = genes_in(adata.var_names, genes)
            if present:
                sc.tl.score_genes(adata, present, score_name=key)

    # Decide the categorical key. 选择用于分析的类别键。
    if args.cluster_key and args.cluster_key in adata.obs:
        cluster_key = args.cluster_key
    elif nz["cluster_key"] in adata.obs:
        cluster_key = nz["cluster_key"]
    else:
        label_spots(adata, cfg)
        cluster_key = "niche_label"
    print(f"[04_niche] using categorical key: {cluster_key}")

    # 2) Spatial graph. 空间图。
    sq.gr.spatial_neighbors(
        adata, coord_type=nz["coord_type"], n_neighs=nz["n_neighs"]
    )

    # 3) Neighbourhood enrichment (permutation z-scores). 邻域富集（置换 z 分数）。
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, seed=cfg["project"]["random_seed"])
    zscore = adata.uns[f"{cluster_key}_nhood_enrichment"]["zscore"]
    cats = list(adata.obs[cluster_key].cat.categories)
    import pandas as pd

    pd.DataFrame(zscore, index=cats, columns=cats).to_csv(
        out / f"{args.sample}_nhood_enrichment_zscore.csv"
    )

    # 4a) Co-occurrence probability vs distance. 共现概率随距离变化。
    try:
        sq.gr.co_occurrence(adata, cluster_key=cluster_key)
    except Exception as exc:
        print(f"[warn] co_occurrence skipped: {exc}")

    # 4b) Compositional niches: cluster spots by neighbour label fractions.
    #     组成型生态位：按邻居标签比例对 spot 聚类。
    conn = adata.obsp["spatial_connectivities"]
    onehot = pd.get_dummies(adata.obs[cluster_key]).to_numpy(dtype=float)
    neigh_comp = conn.dot(onehot)
    rowsum = neigh_comp.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    neigh_comp = neigh_comp / rowsum
    from sklearn.cluster import KMeans

    km = KMeans(n_clusters=nz["n_niches"], random_state=cfg["project"]["random_seed"], n_init=10)
    adata.obs["niche"] = km.fit_predict(neigh_comp).astype(str)
    adata.obs["niche"] = adata.obs["niche"].astype("category")

    # Which niches are enriched for tumour-adjacent TLS/B/T? Save the composition.
    comp_df = pd.DataFrame(neigh_comp, columns=cats, index=adata.obs_names)
    comp_df["niche"] = adata.obs["niche"].values
    comp_df.groupby("niche").mean().to_csv(out / f"{args.sample}_niche_composition.csv")

    adata.write(out / f"{args.sample}_niche.h5ad")
    print(
        f"[04_niche] {args.sample}: labels={adata.obs[cluster_key].value_counts().to_dict()}\n"
        f"           enrichment + composition CSVs -> {out}"
    )


if __name__ == "__main__":
    main()
