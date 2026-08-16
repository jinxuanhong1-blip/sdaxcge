#!/usr/bin/env python3
"""Stream GSE131907 raw UMI text → epithelial nLung+tLung AnnData.

ADDITIVE public slice. Does not load the 2.9 GB log2TPM matrix.
Keeps author labels. LUAD atlas only (Kim et al. 2020).
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


KEEP_ORIGINS = ("tLung", "nLung")
EPI_TYPE = "Epithelial cells"


def _read_annotation(path: Path) -> pd.DataFrame:
    ann = pd.read_csv(path, sep="\t", compression="gzip")
    if "Index" not in ann.columns:
        raise SystemExit(f"annotation missing Index: {list(ann.columns)}")
    return ann


def _parse_series_matrix(path: Path) -> pd.DataFrame:
    """Sample-level GEO characteristics (patient / stage). Best-effort."""
    rows: dict[str, list[str]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0].lstrip("!")
            rows[key] = [p.strip().strip('"') for p in parts[1:]]
    if "Sample_title" not in rows:
        return pd.DataFrame()
    titles = rows["Sample_title"]
    n = len(titles)
    out = pd.DataFrame({"Sample": titles})
    chars = rows.get("Sample_characteristics_ch1", [])
    # GEO repeats characteristics as multiple rows with the same key; some
    # series matrices flatten them. Handle both "tag: value" tokens.
    extra: dict[str, list[str]] = {}
    for key, vals in rows.items():
        if key in {"Sample_title", "Sample_geo_accession", "Sample_source_name_ch1"}:
            if len(vals) == n:
                extra[key] = vals
    # characteristics may be one row per tag with n columns, or n*k mixed.
    # Kim series uses several !Sample_characteristics_ch1 lines in SOFT;
    # the series matrix often concatenates. Parse "key: value" if aligned.
    if chars and len(chars) == n:
        parsed = []
        for cell in chars:
            d = {}
            for tok in str(cell).split(";"):
                if ":" in tok:
                    k, v = tok.split(":", 1)
                    d[k.strip()] = v.strip()
            parsed.append(d)
        char_df = pd.DataFrame(parsed)
        out = pd.concat([out, char_df], axis=1)
    for k, vals in extra.items():
        out[k] = vals
    return out


def stream_umi_subset(
    umi_path: Path,
    keep_ids: list[str],
    gene_cap: int | None = None,
) -> tuple[list[str], list[str], sparse.csr_matrix]:
    keep_set = set(keep_ids)
    with gzip.open(umi_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[0] not in {"Index", "index", "GENE", "gene"}:
            # first field is the gene-name column label
            pass
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
            # Fast path: split once; gene is field 0.
            # Use numpy on the keep columns only.
            raw = line.rstrip("\n")
            if not raw:
                continue
            tab0 = raw.find("\t")
            gene = raw[:tab0]
            rest = raw[tab0 + 1 :]
            # np.fromstring on the full row then take col_idx is simpler and
            # still OK at 30k genes × 208k cells if we stream one row.
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


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", required=True, help="download directory")
    p.add_argument("--out", required=True, help="output .h5ad")
    p.add_argument("--gene-cap", type=int, default=None, help="debug only")
    args = p.parse_args()
    data = Path(args.data)
    ann_path = data / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = data / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    sm_path = data / "GSE131907_series_matrix.txt.gz"
    if not ann_path.is_file() or not umi_path.is_file():
        raise SystemExit("run scripts/scrna_paga/download.sh first")

    ann = _read_annotation(ann_path)
    epi = ann.loc[
        (ann["Cell_type"] == EPI_TYPE) & (ann["Sample_Origin"].isin(KEEP_ORIGINS))
    ].copy()
    inventory = {
        "n_ann_rows": int(len(ann)),
        "n_epithelial_nLung_tLung": int(len(epi)),
        "by_origin": epi["Sample_Origin"].value_counts().to_dict(),
        "by_subtype": epi.get("Cell_subtype", pd.Series(dtype=str))
        .fillna("NA")
        .value_counts()
        .to_dict(),
        "n_samples": int(epi["Sample"].nunique()) if "Sample" in epi.columns else None,
    }
    print(json.dumps(inventory, indent=2), flush=True)
    keep_ids = epi["Index"].astype(str).tolist()

    cells, genes, mat = stream_umi_subset(umi_path, keep_ids, gene_cap=args.gene_cap)
    # cells × genes for AnnData
    X = mat.T.tocsr()
    obs = epi.set_index("Index").reindex(cells)
    if obs["Cell_type"].isna().any():
        n_miss = int(obs["Cell_type"].isna().sum())
        print(f"warning: {n_miss} extracted cells missing annotation", file=sys.stderr)
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))

    try:
        import anndata as ad
    except ImportError as e:
        raise SystemExit("anndata is required to write h5ad") from e

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs["dataset"] = "GSE131907"
    adata.obs["histology"] = "LUAD"
    adata.layers["counts"] = adata.X.copy()

    if sm_path.is_file():
        meta = _parse_series_matrix(sm_path)
        if not meta.empty and "Sample" in adata.obs.columns:
            # titles in this series are sample IDs (e.g. LUNG_T18)
            sample_map = {}
            for _, row in meta.iterrows():
                sample_map[str(row.get("Sample", ""))] = {
                    k: row[k] for k in meta.columns if k != "Sample"
                }
            # attach only if we can match
            matched = adata.obs["Sample"].astype(str).isin(sample_map)
            inventory["series_matrix_sample_match"] = int(matched.sum())

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    inv_path = out.with_suffix(".inventory.json")
    inventory.update(
        {
            "n_cells_written": int(adata.n_obs),
            "n_genes_written": int(adata.n_vars),
            "nnz": int(adata.X.nnz) if sparse.issparse(adata.X) else int(adata.n_obs * adata.n_vars),
            "out": str(out),
        }
    )
    inv_path.write_text(json.dumps(inventory, indent=2))
    print(f"wrote {out}  cells={adata.n_obs} genes={adata.n_vars}", flush=True)


if __name__ == "__main__":
    main()
