#!/usr/bin/env python3
"""Prepare GSE131907 author-malignant counts for pySCENIC.

Two streaming passes over the public UMI matrix. Pass 1 keeps only
per-gene detection and the library size. Pass 2 writes QC malignant
cells into a float32 memmap, genes detected in at least 1% of those cells.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
GEO = Path("/tmp/gse131907_pyscenic/geo")
PREP = Path("/tmp/gse131907_pyscenic/prepared")
CISTARGET = Path("/tmp/gse131907_pyscenic/cistarget")

MATRIX = GEO / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
ANN = GEO / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
SERIES = GEO / "GSE131907_series_matrix.txt.gz"
TF_LIST = CISTARGET / "allTFs_hg38.txt"
GMT = Path("/tmp/gse131907_pyscenic/genesets/h.all.symbols.gmt")

MALIGNANT = {"Malignant cells", "tS1", "tS2", "tS3"}
MIN_UMI = 200
MIN_DET = 0.01
SEED = 131907
N_GRN_CELLS = 3000
N_HVG = 2000
TF_MIN_DET = 0.05


def gene_name(raw: bytes) -> str:
    return raw.decode("ascii", errors="replace").split(".")[0]


def parse_series(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            if key not in fields:
                fields[key] = values
    n = len(fields.get("title", []))
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def program_symbols() -> set[str]:
    genes = set()
    raw = json.loads((ROOT / "methods/gse131907_pyscenic_cldn4/resources/program_sets.json").read_text())
    for key, vals in raw.items():
        if key == "note":
            continue
        genes.update(vals)
    if GMT.exists():
        for line in GMT.read_text().splitlines():
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] in {"HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"}:
                genes.update(parts[2:])
    return genes


def stream_genes(cells_n: int, mal_idx: np.ndarray):
    with gzip.open(MATRIX, "rb") as fh:
        fh.readline()
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = gene_name(raw[:tab])
            arr = np.fromstring(raw[tab + 1:].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != cells_n:
                raise SystemExit(f"{gene}: {arr.size} values, header {cells_n}")
            yield n_genes, gene, arr[mal_idx]


def main() -> None:
    PREP.mkdir(parents=True, exist_ok=True)
    ann = pd.read_csv(ANN, sep="\t", dtype=str)
    ann["author_malignant"] = ann["Cell_subtype"].isin(MALIGNANT)
    print(f"[ann] cells={len(ann)} author_malignant={int(ann.author_malignant.sum())}", flush=True)

    meta = parse_series(SERIES).rename(columns={"title": "sample"})
    keep = [c for c in ["sample", "geo_accession", "patient_id", "tumor_stage",
                        "tissue_origin_abbrevation", "source_name_ch1",
                        "lung_cancer_subtype"] if c in meta.columns]
    meta = meta[keep].drop_duplicates("sample")
    print(f"[series] samples={len(meta)} patients={meta['patient_id'].nunique()}", flush=True)

    with gzip.open(MATRIX, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
    cells = header.split("\t")[1:]
    print(f"[matrix] cells_in_header={len(cells)}", flush=True)
    ann_i = ann.set_index("Index")
    missing = [c for c in cells if c not in ann_i.index]
    if missing:
        raise SystemExit(f"{len(missing)} barcodes missing from annotation, e.g. {missing[:3]}")
    obs = ann_i.reindex(cells).reset_index().rename(columns={"Index": "barcode", "Sample": "sample"})
    obs = obs.merge(meta, on="sample", how="left")
    mal_mask = obs["author_malignant"].fillna(False).to_numpy(dtype=bool)
    mal_idx = np.flatnonzero(mal_mask)
    print(f"[mal] columns={mal_idx.size}", flush=True)

    total = np.zeros(mal_idx.size, dtype=np.float64)
    # Scalars only. A boolean mask per gene would be ~1 GB and a count vector
    # per gene would exceed RAM. Author-malignant cells in this matrix already
    # sit at UMI>=200, so detection on all author-malignant cells matches QC.
    nnz = {}
    order = []
    n_genes = 0
    for n_genes, gene, sub in stream_genes(len(cells), mal_idx):
        total += sub.astype(np.float64)
        hit_n = int((sub > 0).sum())
        nnz[gene] = nnz.get(gene, 0) + hit_n
        if gene not in order:
            order.append(gene)
        if n_genes % 2000 == 0:
            print(f"  pass1 scanned={n_genes}", flush=True)
    print(f"[pass1] genes={n_genes} symbols={len(order)}", flush=True)

    qc = total >= MIN_UMI
    n_qc = int(qc.sum())
    print(f"[qc] UMI>={MIN_UMI}: {n_qc}/{mal_idx.size}", flush=True)
    if n_qc != mal_idx.size:
        print("[qc] warning: some author-malignant cells failed UMI; detection threshold uses all author-malignant cells",
              flush=True)
    obs_m = obs.iloc[mal_idx].reset_index(drop=True)
    obs_m = obs_m.loc[qc].reset_index(drop=True)
    obs_m["total_umi"] = total[qc]

    keep_genes = []
    n_mal = mal_idx.size
    for gene in order:
        det = nnz[gene] / n_mal if n_mal else 0.0
        if det >= MIN_DET or gene == "CLDN4":
            keep_genes.append(gene)
    print(f"[genes] universe={len(keep_genes)} (detection>={MIN_DET})", flush=True)
    del nnz

    gene_col = {g: i for i, g in enumerate(keep_genes)}
    x_path = PREP / "X_counts.npy"
    if x_path.exists():
        x_path.unlink()
    X = np.lib.format.open_memmap(x_path, mode="w+", dtype=np.float32, shape=(n_qc, len(keep_genes)))
    X[:] = 0
    seen = 0
    for n_scanned, gene, sub in stream_genes(len(cells), mal_idx):
        j = gene_col.get(gene)
        if j is None:
            continue
        X[:, j] += sub[qc]
        seen += 1
        if n_scanned % 2000 == 0:
            print(f"  pass2 scanned={n_scanned} written_rows={seen}", flush=True)
    X.flush()
    print(f"[pass2] memmap {X.shape}", flush=True)

    if "CLDN4" not in gene_col:
        raise SystemExit("CLDN4 absent after the detection filter")
    cldn4 = np.array(X[:, gene_col["CLDN4"]], dtype=np.float64)
    obs_m["cldn4_count"] = cldn4
    lib = obs_m["total_umi"].to_numpy(dtype=np.float64)
    lib_safe = lib.copy()
    lib_safe[lib_safe <= 0] = np.nan
    obs_m["cldn4_log1p_cp10k"] = np.log1p(cldn4 / lib_safe * 1e4)

    # Variance of log1p(CP10k) in row chunks so the full dense copy is never built.
    sum_x = np.zeros(len(keep_genes), dtype=np.float64)
    sum_x2 = np.zeros(len(keep_genes), dtype=np.float64)
    step = 2000
    for start in range(0, n_qc, step):
        stop = min(start + step, n_qc)
        block = np.log1p(np.asarray(X[start:stop], dtype=np.float64) / lib_safe[start:stop, None] * 1e4)
        sum_x += np.nansum(block, axis=0)
        sum_x2 += np.nansum(block * block, axis=0)
        del block
    mean = sum_x / n_qc
    var = np.maximum(sum_x2 / n_qc - mean * mean, 0)
    det = (np.asarray(X > 0).sum(axis=0) / n_qc)

    tfs = {ln.strip() for ln in TF_LIST.read_text().splitlines() if ln.strip()}
    hvg = []
    for j in np.argsort(-var):
        g = keep_genes[j]
        if g == "CLDN4" or not np.isfinite(var[j]) or var[j] <= 0:
            continue
        hvg.append(g)
        if len(hvg) >= N_HVG:
            break
    tf_keep = sorted(g for g, d in zip(keep_genes, det) if g in tfs and d >= TF_MIN_DET)
    extra = program_symbols()
    grn_genes = []
    seen_g = set()
    for g in hvg + tf_keep + sorted(extra) + ["CLDN4"]:
        if g in gene_col and g not in seen_g:
            grn_genes.append(g)
            seen_g.add(g)
    print(f"[grn genes] hvg={len(hvg)} tfs={len(tf_keep)} grn={len(grn_genes)}", flush=True)

    rng = np.random.default_rng(SEED)
    counts = obs_m.groupby("sample").size()
    eligible = counts[counts >= 20]
    weights = eligible.to_numpy(dtype=float)
    take = np.maximum(np.floor(weights / weights.sum() * N_GRN_CELLS).astype(int), 8)
    take = np.minimum(take, eligible.to_numpy())
    short = int(N_GRN_CELLS - int(take.sum()))
    order_s = np.argsort(-eligible.to_numpy())
    while short > 0:
        progressed = False
        for i in order_s:
            if take[i] < int(eligible.iloc[i]):
                take[i] += 1
                short -= 1
                progressed = True
                if short == 0:
                    break
        if not progressed:
            break
    chosen = []
    for sample, n_take in zip(eligible.index, take):
        idx = np.flatnonzero(obs_m["sample"].to_numpy() == sample)
        chosen.append(np.sort(rng.choice(idx, size=int(n_take), replace=False)))
    grn_cells = np.sort(np.concatenate(chosen))
    print(f"[grn cells] {grn_cells.size} from {len(eligible)} samples", flush=True)

    (PREP / "genes.txt").write_text("\n".join(keep_genes) + "\n")
    (PREP / "grn_genes.txt").write_text("\n".join(grn_genes) + "\n")
    np.save(PREP / "grn_cells.npy", grn_cells.astype(np.int32))
    obs_m.to_csv(PREP / "obs_malignant.tsv.gz", sep="\t", index=False)
    pd.DataFrame({
        "gene": keep_genes,
        "detection_qc": det,
        "var_log1p": var,
        "is_hvg": [g in set(hvg) for g in keep_genes],
        "is_tf_regulator": [g in set(tf_keep) for g in keep_genes],
        "in_grn": [g in seen_g for g in keep_genes],
    }).to_csv(PREP / "gene_table.tsv.gz", sep="\t", index=False)
    X.flush()
    del X

    summary = {
        "n_cells_matrix": len(cells),
        "n_genes_matrix": n_genes,
        "n_author_malignant": int(mal_mask.sum()),
        "n_malignant_umi_ge_200": n_qc,
        "n_genes_universe": len(keep_genes),
        "min_detection": MIN_DET,
        "n_grn_cells": int(grn_cells.size),
        "n_grn_genes": len(grn_genes),
        "n_hvg": len(hvg),
        "n_tf_regulators": len(tf_keep),
        "n_samples_ge20": int(eligible.size),
        "n_patients_in_qc": int(obs_m["patient_id"].nunique()),
        "seed": SEED,
        "malignant_definition": "Cell_subtype in {Malignant cells, tS1, tS2, tS3}",
    }
    (PREP / "prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
