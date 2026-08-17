#!/usr/bin/env python3
"""Extract epithelium for GSE205335 and GSE207422 separately.

Platforms differ. Do not Harmony-merge. GSE148071 is not added.
GSE205335: author Epithelial cells on non-normal tissues.
GSE207422: marker epithelial (Hu canonical argmax); leftover AT2/club/ciliated
kept for the DPT root. Author CopyKAT is not public.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (  # noqa: E402
    A3_NORMAL_LUNG,
    LINEAGES_207422,
    PAPER_GROUP_207422,
)

NORMAL_TISSUE_PREFIX = ("Normal ",)
CAP_PER_UNIT = 400
SEED = 1


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def cap_barcodes(
    frame: pd.DataFrame,
    unit_col: str,
    cap: int,
    seed: int,
    protect: pd.Series | None = None,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    keep_idx: list[int] = []
    for _, sub in frame.groupby(unit_col, observed=True):
        if protect is not None:
            prot = pd.Index(sub.index[protect.loc[sub.index].to_numpy()])
        else:
            prot = pd.Index(sub.index[:0])
        # Leftover is protected for the DPT root, but one leftover-heavy
        # patient must not exceed the cap (P07 leftover would otherwise dominate).
        max_prot = min(len(prot), max(50, cap // 4))
        if len(prot) > max_prot:
            prot = pd.Index(rng.choice(prot.to_numpy(), size=max_prot, replace=False))
        rest = sub.index.difference(prot)
        n_rest = max(0, cap - len(prot))
        if len(rest) > n_rest:
            chosen = rng.choice(rest.to_numpy(), size=n_rest, replace=False)
            picked = prot.append(pd.Index(chosen))
        else:
            picked = prot.append(rest)
        keep_idx.extend(picked.tolist())
    return frame.loc[keep_idx].copy()


def load_rds_matrix(path: Path):
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if str(path).endswith(".gz"):
            matrix_path = Path(tmp) / "matrix.rds"
            print(f"decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("read GSE205335 RDS", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    return matrix, genes, barcodes


def _sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_bool_dtype(out[col]) or str(out[col].dtype) == "boolean":
            out[col] = out[col].map({True: "True", False: "False"}).astype(str)
        elif out[col].dtype == object or str(out[col].dtype) == "category":
            out[col] = out[col].astype(str).replace({"nan": "NA", "None": "NA", "<NA>": "NA"})
    return out


def extract_gse205335(data: Path, cap: int):
    import anndata as ad

    ident_path = data / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = data / "gse205335" / "GSE205335_family.soft.gz"
    rds_path = data / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ident = pd.read_csv(ident_path, sep="\t")
    meta = parse_geo_soft(soft_path)
    cells = ident.merge(
        meta[
            [
                "orig.ident",
                "gsm",
                "patient",
                "tissue",
                "recist",
                "cancer_subtype",
                "tumor_stage",
                "platform",
            ]
        ],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise SystemExit(f"GSE205335 identity samples missing GEO metadata: {missing}")
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith(NORMAL_TISSUE_PREFIX)
    cells["is_epi"] = cells["lineage.total"].eq("Epithelial cells")
    epi = cells.loc[cells["is_epi"] & ~cells["is_normal_tissue"]].copy()
    leftover = ~epi["lineage.sub"].astype(str).eq("Malignant cells")
    epi = epi.set_index("barcode", drop=False)
    leftover.index = epi.index
    catalog = {
        "n_ident": int(len(ident)),
        "n_epithelial_non_normal": int(len(epi)),
        "n_leftover_nonmalignant": int(leftover.sum()),
        "by_lineage_sub": epi["lineage.sub"].value_counts().to_dict(),
        "by_tissue": epi["tissue"].value_counts().to_dict(),
        "by_recist": epi.groupby("patient")["recist"].first().value_counts().to_dict(),
        "n_patients": int(epi["patient"].nunique()),
        "n_normal_epi_dropped": int((cells["is_epi"] & cells["is_normal_tissue"]).sum()),
        "n_normal_only_patients_dropped": int(
            cells.loc[cells["is_normal_tissue"], "patient"].nunique()
            - cells.loc[~cells["is_normal_tissue"] & cells["is_epi"], "patient"].nunique()
        ),
    }
    print(json.dumps({"GSE205335_catalog": catalog}, indent=2), flush=True)
    epi = cap_barcodes(epi, "patient", cap, SEED, leftover)
    keep = set(epi["barcode"].astype(str))
    matrix, genes, barcodes = load_rds_matrix(rds_path)
    col_idx = np.array([i for i, b in enumerate(barcodes) if b in keep], dtype=np.int64)
    print(f"GSE205335 subset columns {len(col_idx)} / {len(barcodes)}", flush=True)
    sub = matrix[:, col_idx].tocsr()
    del matrix
    ordered = [str(barcodes[i]) for i in col_idx]
    obs = epi.set_index("barcode").reindex(ordered)
    adata = ad.AnnData(
        X=sub.T.tocsr(),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata.obs["dataset"] = "GSE205335"
    adata.obs["histology"] = adata.obs["cancer_subtype"].astype(str)
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["unit_id"] = adata.obs["patient"].astype(str)
    adata.obs["author_lineage"] = adata.obs["lineage.total"].astype(str)
    adata.obs["author_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["is_malignant"] = adata.obs["lineage.sub"].astype(str).eq("Malignant cells")
    adata.obs["is_leftover"] = ~adata.obs["is_malignant"]
    adata.obs["timing"] = "palliative_ici"
    adata.obs["benefit_raw"] = adata.obs["recist"].astype(str)
    adata.obs["platform"] = adata.obs["platform"].astype(str)
    adata.layers["counts"] = adata.X.copy()
    adata.obs = _sanitize_obs(adata.obs)
    print(f"GSE205335 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata, catalog


def _score_lineage(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES_207422)
    scores = np.vstack([_score_lineage(expr, LINEAGES_207422[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def stream_207422_markers(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, np.ndarray]:
    keep = set()
    for genes in LINEAGES_207422.values():
        keep.update(genes)
    keep.update(A3_NORMAL_LUNG)
    keep.update(["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "SFTPC", "SFTPA1", "SFTPB"])
    print(f"[GSE207422] marker pass {path}", flush=True)
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        n_genes_cell = np.zeros(n, dtype=np.int32)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, rest = line.split("\t", 1)
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"column mismatch for {gene}: {vals.size} != {n}")
            totals += vals
            n_genes_cell += vals > 0
            if gene in keep:
                store[gene] = vals
            if n_genes % 4000 == 0:
                print(f"  marker pass {n_genes} genes, kept {len(store)}", flush=True)
    print(f"[GSE207422] cells={n} genes={n_genes} markers={len(store)}", flush=True)
    return cells, store, totals, n_genes_cell


def stream_207422_subset(path: Path, keep_ids: list[str]):
    keep_set = set(keep_ids)
    print(f"[GSE207422] full-gene pass keep={len(keep_set)}", flush=True)
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep_set]
        if not col_idx:
            raise SystemExit("no requested GSE207422 cell IDs in UMI header")
        ordered_cells = [cell_ids[i] for i in col_idx]
        genes: list[str] = []
        data: list[np.ndarray] = []
        indices: list[np.ndarray] = []
        indptr = [0]
        nnz_total = 0
        for gi, line in enumerate(fh, start=1):
            raw = line.rstrip("\n")
            if not raw:
                continue
            tab0 = raw.find("\t")
            gene = raw[:tab0]
            rest = raw[tab0 + 1 :]
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != len(cell_ids):
                raise SystemExit(f"row {gi} {gene}: expected {len(cell_ids)} values, got {vals.size}")
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                data.append(sub[nz].astype(np.float32, copy=False))
                indices.append(nz.astype(np.int32, copy=False))
                nnz_total += int(nz.size)
            indptr.append(nnz_total)
            genes.append(gene)
            if gi % 4000 == 0:
                print(f"  streamed {gi} genes, nnz={nnz_total}", flush=True)
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


def load_207422_sample_meta(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path)
    raw = raw.dropna(subset=["Sample"]).copy()
    raw = raw[~raw["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    raw["Sample"] = raw["Sample"].astype(str)
    raw["paper_group"] = raw["Sample"].map(PAPER_GROUP_207422)
    raw["path_response"] = raw["Pathologic Response"].replace({"pCR": "MPR"})
    raw["timing"] = np.where(
        raw["Resource"].astype(str).str.contains("Pre", case=False, na=False),
        "pre",
        "post",
    )
    return raw


def extract_gse207422(data: Path, cap: int, post_only: bool = True):
    import anndata as ad

    umi_path = data / "gse207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_path = data / "gse207422" / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    meta = load_207422_sample_meta(meta_path)
    cells, expr, totals, n_genes_cell = stream_207422_markers(umi_path)
    n = len(cells)
    sample = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    lineage = assign_lineage(expr, n)
    epi_mask = lineage == "epithelial"
    normal_umi = np.zeros(n, dtype=np.float64)
    for g in A3_NORMAL_LUNG:
        if g in expr:
            normal_umi += expr[g]
    leftover = epi_mask & (normal_umi > 0)
    malig = epi_mask & (normal_umi == 0)
    timing = pd.Series(sample).map(meta.set_index("Sample")["timing"]).to_numpy()
    paper = pd.Series(sample).map(PAPER_GROUP_207422).fillna("NA").to_numpy()
    if post_only:
        keep_time = timing == "post"
    else:
        keep_time = np.ones(n, dtype=bool)
    keep_epi = epi_mask & keep_time
    catalog = {
        "n_cells_matrix": int(n),
        "n_epithelial": int(epi_mask.sum()),
        "n_epithelial_post": int((epi_mask & (timing == "post")).sum()),
        "n_epithelial_pre": int((epi_mask & (timing == "pre")).sum()),
        "n_a3_malignant_like": int(malig.sum()),
        "n_leftover": int(leftover.sum()),
        "n_leftover_post": int((leftover & (timing == "post")).sum()),
        "by_sample_epithelial": {
            str(s): int(((sample == s) & epi_mask).sum()) for s in sorted(set(sample))
        },
        "post_only": post_only,
    }
    print(json.dumps({"GSE207422_catalog": catalog}, indent=2), flush=True)

    frame = pd.DataFrame(
        {
            "barcode": cells[keep_epi],
            "sample": sample[keep_epi],
            "lineage": lineage[keep_epi],
            "is_malignant": malig[keep_epi],
            "is_leftover": leftover[keep_epi],
            "timing": timing[keep_epi],
            "paper_group": paper[keep_epi],
            "n_umi": totals[keep_epi],
            "n_genes": n_genes_cell[keep_epi],
        }
    )
    frame["patient_id"] = frame["sample"].str.replace("BD_immune", "P", regex=False)
    frame["unit_id"] = frame["patient_id"]
    frame = frame.set_index("barcode", drop=False)
    protect = frame["is_leftover"].astype(bool)
    frame = cap_barcodes(frame, "unit_id", cap, SEED + 2, protect)
    keep_ids = frame["barcode"].astype(str).tolist()
    ordered, genes, mat = stream_207422_subset(umi_path, keep_ids)
    X = mat.T.tocsr()
    obs = frame.set_index("barcode").reindex(ordered)
    # re-attach sample-level fields
    meta_idx = meta.set_index("Sample")
    obs["histology"] = obs["sample"].map(meta_idx["Pathology"]).astype(str)
    obs["recist"] = obs["sample"].map(meta_idx["RECIST"]).astype(str)
    obs["path_response"] = obs["sample"].map(meta_idx["path_response"]).astype(str)
    obs["dataset"] = "GSE207422"
    obs["author_lineage"] = "marker_epithelial"
    obs["author_subtype"] = np.where(obs["is_malignant"].astype(str).eq("True"), "A3_malignant_like", "leftover_epi")
    obs["benefit_raw"] = obs["paper_group"].astype(str)
    obs["platform"] = "10x_GEX_Hu2023"
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.layers["counts"] = adata.X.copy()
    adata.obs = _sanitize_obs(adata.obs)
    print(f"GSE207422 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata, catalog


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/ici_pair_205335_207422"))
    p.add_argument("--out", type=Path, default=Path("/tmp/ici_pair_205335_207422/extracted"))
    p.add_argument("--cap-per-unit", type=int, default=CAP_PER_UNIT)
    p.add_argument("--include-207422-pre", action="store_true")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    inv = {"cap_per_unit": args.cap_per_unit, "harmony_joint": False, "gse148071": False}

    cache_a = args.out / "GSE205335_epithelium.h5ad"
    if cache_a.is_file():
        print(f"reuse {cache_a}", flush=True)
        inv["GSE205335"] = {"cached": True, "path": str(cache_a)}
    else:
        a, cat_a = extract_gse205335(args.data, args.cap_per_unit)
        a.write_h5ad(cache_a)
        inv["GSE205335"] = {"cached": False, "path": str(cache_a), "n_cells": int(a.n_obs), "catalog": cat_a}
        del a

    cache_b = args.out / "GSE207422_epithelium.h5ad"
    if cache_b.is_file():
        print(f"reuse {cache_b}", flush=True)
        inv["GSE207422"] = {"cached": True, "path": str(cache_b)}
    else:
        b, cat_b = extract_gse207422(args.data, args.cap_per_unit, post_only=not args.include_207422_pre)
        b.write_h5ad(cache_b)
        inv["GSE207422"] = {"cached": False, "path": str(cache_b), "n_cells": int(b.n_obs), "catalog": cat_b}
        del b

    (args.out / "extract_inventory.json").write_text(json.dumps(inv, indent=2, default=str))
    print(json.dumps(inv, indent=2, default=str), flush=True)


if __name__ == "__main__":
    main()
