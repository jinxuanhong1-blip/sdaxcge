#!/usr/bin/env python3
"""scirpy 0.25: 10x VDJ -> clone_id + within-sample expansion.

Reads a folder of per-sample `filtered_contig_annotations.csv` (or a single
CSV). Writes clone-level and cell-level tables. Does not download data.

Example:
  python3 01_scirpy_clonotypes.py \\
      --contigs-glob 'data/public/gse243013/tcr/*/filtered_contig_annotations.csv' \\
      --sample-from-parent \\
      --out workdir/scirpy
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import anndata as ad
import pandas as pd
import scirpy as ir


def _sample_id_from_path(path: Path, sample_from_parent: bool) -> str:
    if sample_from_parent:
        return path.parent.name
    stem = path.name
    for suffix in (
        "_filtered_contig_annotations.csv",
        "filtered_contig_annotations.csv",
        ".csv",
    ):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)] or path.stem
            break
    return stem or path.stem


def load_airr(paths: list[Path], sample_from_parent: bool) -> ad.AnnData:
    parts = []
    for path in paths:
        airr = ir.io.read_10x_vdj(path)
        sid = _sample_id_from_path(path, sample_from_parent)
        airr.obs["sample_id"] = sid
        # 10x barcodes repeat across captures; make them unique for concat
        airr.obs_names = pd.Index([f"{sid}_{b}" for b in airr.obs_names], name="cell_id")
        parts.append(airr)
    if len(parts) == 1:
        return parts[0]
    return ad.concat(parts, index_unique=None)


def define_clones(adata: ad.AnnData) -> ad.AnnData:
    ir.pp.index_chains(adata)
    ir.tl.chain_qc(adata)
    pairing = adata.obs["chain_pairing"]
    rtype = adata.obs["receptor_type"]
    keep = rtype.isin(["TCR"]) & ~pairing.isin(
        ["multichain", "orphan VDJ", "orphan VJ", "no IR"]
    )
    adata = adata[keep].copy()
    ir.pp.ir_dist(adata, metric="identity", sequence="nt")
    ir.tl.define_clonotypes(
        adata,
        receptor_arms="all",
        dual_ir="primary_only",
        key_added="clone_id",
    )
    ir.tl.clonal_expansion(
        adata,
        target_col="clone_id",
        expanded_in="sample_id",
        breakpoints=(1, 2),
        key_added="clonal_expansion",
    )
    return adata


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contigs-glob", required=True, help="Glob of 10x contig CSVs")
    p.add_argument(
        "--sample-from-parent",
        action="store_true",
        help="Use the parent directory name as sample_id",
    )
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    paths = sorted(Path(p) for p in glob.glob(args.contigs_glob))
    if not paths:
        raise SystemExit(f"no files matched {args.contigs_glob}")

    adata = load_airr(paths, args.sample_from_parent)
    adata = define_clones(adata)

    # Single-modality AnnData: columns are unprefixed.
    # MuData users must keep the airr: prefix (see playbook §3.3).
    airr_tbl = ir.get.airr(adata, ["junction_aa", "v_call", "j_call", "umi_count"], ["VJ_1", "VDJ_1"])
    cells = adata.obs[["sample_id", "clone_id", "clonal_expansion", "receptor_type", "chain_pairing"]].copy()
    cells = cells.join(airr_tbl)
    cells["clone_size"] = cells.groupby(["sample_id", "clone_id"], observed=True)["clone_id"].transform(
        "size"
    )
    cells["expanded"] = cells["clone_size"] >= 2

    clones = (
        cells.groupby(["sample_id", "clone_id"], observed=True)
        .agg(
            clone_size=("clone_id", "size"),
            VJ_1_junction_aa=("VJ_1_junction_aa", "first"),
            VDJ_1_junction_aa=("VDJ_1_junction_aa", "first"),
            VJ_1_v_call=("VJ_1_v_call", "first"),
            VDJ_1_v_call=("VDJ_1_v_call", "first"),
        )
        .reset_index()
    )

    args.out.mkdir(parents=True, exist_ok=True)
    cells.to_csv(args.out / "cells.tsv", sep="\t")
    clones.to_csv(args.out / "clone_table.tsv", sep="\t", index=False)
    qc = (
        adata.obs.groupby(["sample_id", "chain_pairing"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    qc.to_csv(args.out / "qc_chain_pairing.tsv", sep="\t", index=False)
    print(f"wrote {args.out}  cells={cells.shape[0]} clones={clones.shape[0]}")


if __name__ == "__main__":
    main()
