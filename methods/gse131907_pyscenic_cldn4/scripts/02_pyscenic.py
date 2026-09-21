#!/usr/bin/env python3
"""GRNBoost2 + cisTarget + AUCell on GSE131907 malignant cells.

pySCENIC 0.12.1 touches numpy.object, which NumPy 2 removed. The alias is
restored before import. The GRN math is unchanged.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

# Must run before pyscenic imports transform.py.
try:
    _ = np.object  # noqa: SLF001
except AttributeError:
    np.object = object

import pandas as pd
import dask.dataframe as dd
from arboreto.algo import grnboost2
import arboreto.core as arboreto_core
from ctxcore.genesig import GeneSignature, Regulon
from ctxcore.rnkdb import FeatherRankingDatabase
from dask.distributed import Client, LocalCluster
from pyscenic.aucell import aucell
from pyscenic.prune import df2regulons, prune2df
from pyscenic import utils as pyscenic_utils
from pyscenic.utils import modules_from_adjacencies

PREP = Path("/tmp/gse131907_pyscenic/prepared")
CIS = Path("/tmp/gse131907_pyscenic/cistarget")
GMT = Path("/tmp/gse131907_pyscenic/genesets/h.all.symbols.gmt")
ROOT = Path(__file__).resolve().parents[3]
PROG = ROOT / "methods/gse131907_pyscenic_cldn4/resources/program_sets.json"

ADJ = PREP / "adjacencies.tsv.gz"
REG_TSV = PREP / "regulons.tsv"
AUC_NPY = PREP / "auc.npy"
AUC_COLS = PREP / "auc_columns.txt"
SEED = 131907
NES_THRESHOLD = 3.0
AUC_THRESHOLD = 0.05
DBS = [
    (
        CIS / "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "hg38_500bp_v10_clust",
    ),
    (
        CIS / "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "hg38_10kbp_v10_clust",
    ),
]
MOTIF_TBL = CIS / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"


def _patch_arboreto_dask() -> None:
    """arboreto 0.1.6 always builds a meta dataframe, even when it is unused.

    dask 2026 from_delayed rejects an empty list. The GRNBoost2 regressor is
    unchanged; only the empty meta frame is skipped.
    """
    original = dd.from_delayed

    def from_delayed_compat(dfs, meta=None, **kwargs):
        if isinstance(dfs, (list, tuple)) and len(dfs) == 0:
            empty = meta if isinstance(meta, pd.DataFrame) else pd.DataFrame()
            return dd.from_pandas(empty, npartitions=1)
        return original(dfs, meta=meta, **kwargs)

    arboreto_core.from_delayed = from_delayed_compat


def _patch_pyscenic_pandas3() -> None:
    """pandas 3 drops the grouping column inside DataFrameGroupBy.apply.

    modules4top_factors needs the target column after selecting the top
    regulators of each gene. The replacement does that selection without apply.
    """
    def modules4top_factors(adjacencies, n, context=frozenset()):
        parts = []
        weight = pyscenic_utils.COLUMN_NAME_WEIGHT
        target = pyscenic_utils.COLUMN_NAME_TARGET
        tf_col = pyscenic_utils.COLUMN_NAME_TF
        for _, grp in adjacencies.groupby(by=target, sort=False):
            parts.append(grp.nlargest(n, weight))
        df = pd.concat(parts, ignore_index=True) if parts else adjacencies.iloc[0:0]
        for tf_name, df_grp in df.groupby(by=tf_col, sort=False):
            if len(df_grp) == 0:
                continue
            yield pyscenic_utils.Regulon(
                name=tf_name,
                context=frozenset(["top{}perTarget".format(n)]).union(context),
                transcription_factor=tf_name,
                gene2weight=list(zip(df_grp[target].values, df_grp[weight].values)),
                gene2occurrence=[],
            )

    pyscenic_utils.modules4top_factors = modules4top_factors


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def holdout_cldn4(reg: Regulon) -> Regulon | None:
    pairs = [(g, float(w)) for g, w in reg.gene2weight.items() if g != "CLDN4"]
    if len(pairs) < 10:
        return None
    return reg.copy(name=f"{reg.name}|CLDN4out", gene2weight=pairs)


def load_programs(universe: set[str]) -> list[GeneSignature]:
    raw = json.loads(PROG.read_text())
    sets: dict[str, list[str]] = {}
    for key, genes in raw.items():
        if key == "note":
            continue
        sets[key] = [g for g in genes if g in universe and g != "CLDN4"]
    if GMT.exists():
        for line in GMT.read_text().splitlines():
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            if parts[0] in {"HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"}:
                sets[parts[0]] = [g for g in parts[2:] if g in universe and g != "CLDN4"]
        alpha = set(sets.get("HALLMARK_INTERFERON_ALPHA_RESPONSE", []))
        gamma = set(sets.get("HALLMARK_INTERFERON_GAMMA_RESPONSE", []))
        sets["IFN_UNION"] = sorted(alpha | gamma)
    sigs = []
    used = {}
    for name, genes in sets.items():
        genes = sorted(set(genes))
        used[name] = genes
        if len(genes) < 8:
            log(f"skip program {name}: only {len(genes)} genes in universe")
            continue
        sigs.append(GeneSignature(name=f"PROGRAM::{name}", gene2weight=genes))
    (PREP / "program_genes_used.json").write_text(json.dumps(used, indent=2) + "\n")
    return sigs


def regulons_to_frame(regs: list[Regulon]) -> pd.DataFrame:
    rows = []
    for reg in regs:
        genes = list(reg.genes)
        weights = list(reg.weights)
        rows.append({
            "name": reg.name,
            "tf": reg.transcription_factor,
            "nes": float(getattr(reg, "nes", np.nan)),
            "score": float(reg.score),
            "n_targets": len(genes),
            "contains_CLDN4": "CLDN4" in set(genes),
            "context": ";".join(sorted(reg.context)) if reg.context else "",
            "targets": ";".join(genes),
            "weights": ";".join(f"{w:.6g}" for w in weights),
        })
    return pd.DataFrame(rows)


def main() -> None:
    _patch_arboreto_dask()
    _patch_pyscenic_pandas3()
    genes = (PREP / "genes.txt").read_text().splitlines()
    grn_genes = (PREP / "grn_genes.txt").read_text().splitlines()
    grn_cells = np.load(PREP / "grn_cells.npy")
    obs = pd.read_csv(PREP / "obs_malignant.tsv.gz", sep="\t")
    gene_index = {g: i for i, g in enumerate(genes)}
    grn_idx = np.array([gene_index[g] for g in grn_genes], dtype=int)
    log(f"mmap counts; grn cells={grn_cells.size} genes={len(grn_genes)} universe={len(genes)}")
    X = np.load(PREP / "X_counts.npy", mmap_mode="r")
    sub = np.array(X[np.ix_(grn_cells, grn_idx)], dtype=np.float64)
    del X
    lib = obs["total_umi"].to_numpy(dtype=np.float64)[grn_cells]
    lib[lib <= 0] = np.nan
    logX = np.log1p(sub / lib[:, None] * 1e4)
    tf_table = pd.read_csv(PREP / "gene_table.tsv.gz", sep="\t")
    tfs = tf_table.loc[tf_table["is_tf_regulator"] & tf_table["gene"].isin(grn_genes), "gene"].tolist()
    log(f"TFs for GRNBoost2: {len(tfs)}")
    ex = pd.DataFrame(logX, columns=grn_genes)
    # Drop columns with no variance so GBM does not see constants.
    nunique = (ex.var(axis=0) > 0).to_numpy()
    ex = ex.loc[:, nunique]
    tfs = [t for t in tfs if t in ex.columns]
    log(f"GRN matrix {ex.shape} after dropping constant genes; TFs={len(tfs)}")

    if ADJ.exists() and ADJ.stat().st_size > 0:
        log(f"load cached adjacencies {ADJ}")
        adj = pd.read_csv(ADJ, sep="\t")
    else:
        log("GRNBoost2 starting")
        cluster = LocalCluster(
            n_workers=3,
            threads_per_worker=1,
            memory_limit="1.6GiB",
            processes=True,
            silence_logs=40,
        )
        client = Client(cluster)
        try:
            t0 = time.time()
            adj = grnboost2(
                ex,
                tf_names=tfs,
                client_or_address=client,
                seed=SEED,
                verbose=True,
            )
            log(f"GRNBoost2 done in {(time.time() - t0) / 60:.1f} min, links={len(adj)}")
        finally:
            client.close()
            cluster.close()
        adj.to_csv(ADJ, sep="\t", index=False)
        log(f"wrote {ADJ}")

    need_prune = not (REG_TSV.exists() and REG_TSV.stat().st_size > 0)
    if need_prune:
        log("modules_from_adjacencies")
        modules = modules_from_adjacencies(adj, ex, rho_mask_dropouts=False)
        log(f"modules={len(modules)}")
        dbs = [FeatherRankingDatabase(str(path), name=name) for path, name in DBS]
        log("cisTarget prune2df NES>=%s rank_threshold=1500" % NES_THRESHOLD)
        del ex
        t1 = time.time()
        df = prune2df(
            dbs,
            modules,
            str(MOTIF_TBL),
            rank_threshold=1500,
            auc_threshold=AUC_THRESHOLD,
            nes_threshold=NES_THRESHOLD,
            client_or_address="custom_multiprocessing",
            num_workers=2,
            module_chunksize=40,
        )
        log(f"prune2df rows={len(df)} in {(time.time() - t1) / 60:.1f} min")
        if df.empty:
            raise SystemExit("cisTarget returned no enriched regulons at NES 3.0")
        regs = df2regulons(df, save_columns=["NES"])
        log(f"regulons={len(regs)}")
        frame = regulons_to_frame(regs)
        frame.to_csv(REG_TSV, sep="\t", index=False)
        # Keep the enriched-feature table for audit.
        df.to_csv(PREP / "cistarget_features.tsv.gz", sep="\t", index=False)
    else:
        log("skip prune; reload regulons.tsv")
        frame = pd.read_csv(REG_TSV, sep="\t")
        regs = []
        for rec in frame.itertuples(index=False):
            genes_r = str(rec.targets).split(";") if isinstance(rec.targets, str) else []
            weights = [float(x) for x in str(rec.weights).split(";")] if isinstance(rec.weights, str) else []
            pairs = list(zip(genes_r, weights))
            if not pairs:
                continue
            regs.append(Regulon(
                name=rec.name,
                gene2weight=pairs,
                gene2occurrence=[],
                transcription_factor=rec.tf,
                score=float(rec.score) if pd.notna(rec.score) else 0.0,
                nes=float(rec.nes) if pd.notna(rec.nes) else 0.0,
            ))

    if AUC_NPY.exists() and AUC_COLS.exists():
        log("auc already cached")
        return

    scored: list = []
    for reg in regs:
        scored.append(reg)
        held = holdout_cldn4(reg)
        if held is not None:
            scored.append(held)
    programs = load_programs(set(genes))
    scored.extend(programs)
    log(f"AUCell signatures={len(scored)} (regulons + CLDN4-held-out + programs)")

    X = np.load(PREP / "X_counts.npy", mmap_mode="r")
    n_cells = X.shape[0]
    cols = None
    chunks = []
    t2 = time.time()
    step = 2000
    for start in range(0, n_cells, step):
        stop = min(start + step, n_cells)
        block = pd.DataFrame(np.array(X[start:stop], dtype=np.float32), columns=genes)
        part = aucell(
            block,
            scored,
            auc_threshold=AUC_THRESHOLD,
            noweights=False,
            normalize=False,
            seed=SEED,
            num_workers=1,
        )
        # Align columns on the first chunk.
        if cols is None:
            cols = list(part.columns)
        part = part.reindex(columns=cols)
        chunks.append(part.to_numpy(dtype=np.float32))
        log(f"  aucell cells {stop}/{n_cells}")
    auc = np.vstack(chunks)
    np.save(AUC_NPY, auc)
    AUC_COLS.write_text("\n".join(map(str, cols)) + "\n")
    meta = {
        "n_cells": int(n_cells),
        "n_signatures": len(cols),
        "nes_threshold": NES_THRESHOLD,
        "auc_threshold": AUC_THRESHOLD,
        "rank_threshold": 1500,
        "n_regulons": int(frame.shape[0]),
        "n_tfs_grn": len(tfs),
        "grn_matrix": "log1p(CP10k) on the stratified malignant subsample; library size from all genes",
        "aucell_matrix": "raw counts of genes detected in >=1% of QC malignant cells; within-cell ranks match log1p(CP10k)",
        "databases": [name for _, name in DBS],
        "motif_annotation": MOTIF_TBL.name,
        "seed": SEED,
        "minutes_aucell": (time.time() - t2) / 60,
        "pyscenic": "0.12.1",
        "numpy_object_shim": True,
        "cistarget_run": True,
        "scenicplus_run": False,
        "scenicplus_reason": "GSE131907 has no matched scATAC on GEO",
    }
    (PREP / "pyscenic_run.json").write_text(json.dumps(meta, indent=2) + "\n")
    log(f"AUCell wrote {auc.shape} in {meta['minutes_aucell']:.1f} min")


if __name__ == "__main__":
    main()
