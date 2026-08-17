#!/usr/bin/env python3
"""Stream GSE131907 raw UMI → epithelium / malignant AnnData (Kim 2020).

Keeps author Epithelial cells and author Malignant cells on
nLung / tLung / tL/B / mLN / mBrain. PE unlabeled epithelium is dropped.
Immune and stroma are not extracted. nLung AT2 is kept for the Slingshot root.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

KEEP_ORIGINS = ("tLung", "nLung", "tL/B", "mLN", "mBrain")
EPI_TYPES = {"Epithelial cells", "Malignant cells"}
AUTHOR_AT2 = {"AT2"}
CAP_PER_SAMPLE = 400
SEED = 1


def _read_annotation(path: Path) -> pd.DataFrame:
    ann = pd.read_csv(path, sep="\t", compression="gzip")
    if "Index" not in ann.columns:
        raise SystemExit(f"annotation missing Index: {list(ann.columns)}")
    return ann


def stream_umi_subset(
    umi_path: Path,
    keep_ids: list[str],
    gene_cap: int | None = None,
) -> tuple[list[str], list[str], sparse.csr_matrix]:
    keep_set = set(keep_ids)
    with gzip.open(umi_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep_set]
        if not col_idx:
            raise SystemExit("no requested cell IDs in UMI header")
        ordered_cells = [cell_ids[i] for i in col_idx]
        missing = keep_set.difference(cell_ids)
        print(
            f"UMI header cells={len(cell_ids)} keep={len(col_idx)} "
            f"missing_in_matrix={len(missing)}",
            flush=True,
        )
        genes: list[str] = []
        data: list[np.ndarray] = []
        indices: list[np.ndarray] = []
        indptr = [0]
        nnz_total = 0
        for gi, line in enumerate(f, start=1):
            raw = line.rstrip("\n")
            if not raw:
                continue
            tab0 = raw.find("\t")
            gene = raw[:tab0]
            rest = raw[tab0 + 1 :]
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != len(cell_ids):
                raise SystemExit(
                    f"row {gi} {gene}: expected {len(cell_ids)} values, got {vals.size}"
                )
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                data.append(sub[nz].astype(np.float32, copy=False))
                indices.append(nz.astype(np.int32, copy=False))
                nnz_total += int(nz.size)
            indptr.append(nnz_total)
            genes.append(gene)
            if gi % 2000 == 0:
                print(f"  streamed {gi} genes, nnz={nnz_total}", flush=True)
            if gene_cap is not None and gi >= gene_cap:
                break
    if data:
        data_a = np.concatenate(data)
        indices_a = np.concatenate(indices)
    else:
        data_a = np.array([], dtype=np.float32)
        indices_a = np.array([], dtype=np.int32)
    mat = sparse.csr_matrix(
        (data_a, indices_a, np.asarray(indptr, dtype=np.int64)),
        shape=(len(genes), len(ordered_cells)),
        dtype=np.float32,
    )
    return ordered_cells, genes, mat


def cap_barcodes(
    frame: pd.DataFrame,
    cap: int,
    seed: int,
    protect: pd.Series,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    keep_idx: list[int] = []
    for _, sub in frame.groupby("Sample", observed=True):
        prot = sub.index[protect.loc[sub.index].to_numpy()]
        prot = pd.Index(prot)
        rest = sub.index.difference(prot)
        n_rest = max(0, cap - len(prot))
        if len(rest) > n_rest:
            chosen = rng.choice(rest.to_numpy(), size=n_rest, replace=False)
            picked = prot.append(pd.Index(chosen))
        else:
            picked = prot.append(rest)
        keep_idx.extend(picked.tolist())
    return frame.loc[keep_idx].copy()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", required=True, help="download directory")
    p.add_argument("--out", required=True, help="output .h5ad")
    p.add_argument("--gene-cap", type=int, default=None, help="debug only")
    p.add_argument("--cap-per-sample", type=int, default=CAP_PER_SAMPLE)
    args = p.parse_args()
    data = Path(args.data)
    ann_path = data / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = data / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    if not ann_path.is_file() or not umi_path.is_file():
        raise SystemExit("run scripts/download.py first")

    ann = _read_annotation(ann_path)
    cell_type = ann["Cell_type"].astype(str)
    subtype = (
        ann["Cell_subtype"].astype(str)
        if "Cell_subtype" in ann.columns
        else pd.Series("NA", index=ann.index)
    )
    keep_type = cell_type.isin(EPI_TYPES) | subtype.isin(
        {"tS1", "tS2", "tS3", "Malignant cells"}
    )
    epi = ann.loc[keep_type & ann["Sample_Origin"].isin(KEEP_ORIGINS)].copy()
    catalog = {
        "n_ann_rows": int(len(ann)),
        "n_epithelium_malignant_catalog": int(len(epi)),
        "by_origin": epi["Sample_Origin"].value_counts().to_dict(),
        "by_cell_type": epi["Cell_type"].astype(str).value_counts().to_dict(),
        "by_subtype": subtype.loc[epi.index].fillna("NA").value_counts().to_dict(),
        "n_samples_catalog": int(epi["Sample"].nunique()) if "Sample" in epi.columns else None,
        "dropped_origins": sorted(
            set(ann.loc[keep_type, "Sample_Origin"].astype(str).unique())
            - set(KEEP_ORIGINS)
        ),
    }
    print(json.dumps({"catalog": catalog}, indent=2), flush=True)

    protect = epi["Cell_subtype"].astype(str).isin(AUTHOR_AT2) & (
        epi["Sample_Origin"].astype(str) == "nLung"
    )
    if args.cap_per_sample and args.cap_per_sample > 0:
        epi_use = cap_barcodes(epi, args.cap_per_sample, SEED, protect)
    else:
        epi_use = epi
    inventory = {
        **catalog,
        "cap_per_sample": int(args.cap_per_sample),
        "n_cells_after_cap": int(len(epi_use)),
        "n_protected_nLung_AT2": int(protect.sum()),
        "n_protected_kept": int(
            epi_use["Cell_subtype"].astype(str).isin(AUTHOR_AT2).sum()
        ),
        "by_origin_after_cap": epi_use["Sample_Origin"].value_counts().to_dict(),
        "n_samples_after_cap": int(epi_use["Sample"].nunique()),
    }
    print(json.dumps({"extract": inventory}, indent=2), flush=True)
    keep_ids = epi_use["Index"].astype(str).tolist()

    cells, genes, mat = stream_umi_subset(umi_path, keep_ids, gene_cap=args.gene_cap)
    X = mat.T.tocsr()
    obs = epi_use.set_index("Index").reindex(cells)
    if obs["Cell_type"].isna().any():
        n_miss = int(obs["Cell_type"].isna().sum())
        print(f"warning: {n_miss} extracted cells missing annotation", file=sys.stderr)
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))

    try:
        import anndata as ad
    except ImportError as e:
        raise SystemExit("anndata is required to write h5ad") from e

    def _sanitize(obj):
        if isinstance(obj, dict):
            return {str(k).replace("/", "_"): _sanitize(v) for k, v in obj.items()}
        return obj

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs["dataset"] = "GSE131907"
    adata.obs["histology"] = "LUAD"
    adata.obs["author_subtype"] = adata.obs.get(
        "Cell_subtype", pd.Series("NA", index=adata.obs.index)
    ).astype(str)
    adata.layers["counts"] = adata.X.copy()
    # h5ad uns keys cannot contain '/'; keep a sanitized copy for analyze.py
    adata.uns["extract_inventory"] = _sanitize(inventory)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    inventory.update(
        {
            "n_cells_written": int(adata.n_obs),
            "n_genes_written": int(adata.n_vars),
            "nnz": int(adata.X.nnz) if sparse.issparse(adata.X) else int(adata.n_obs * adata.n_vars),
            "out": str(out),
        }
    )
    inv_path = out.with_suffix(".inventory.json")
    inv_path.write_text(json.dumps(inventory, indent=2))
    adata.write_h5ad(out)
    print(f"wrote {out}  cells={adata.n_obs} genes={adata.n_vars}", flush=True)


if __name__ == "__main__":
    main()
