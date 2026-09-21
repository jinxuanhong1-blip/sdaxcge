#!/usr/bin/env python3
"""Extract QC-pass malignant cells from the concordant-4 public matrices.

Malignant definition matches the locked concordant-4 analysis:
  GSE123902, GSE189357: (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0
  GSE131907: author Cell_subtype == "Malignant cells" on tumor-bearing samples
  GSE205335: author lineage.sub == "Malignant cells", normal tissue dropped

Counts saved are the gene universe only (TFs, program genes, CLDN4, controls).
Per-cell n_count / n_feat / n_mt are from the full transcriptome so CP10k
is not computed on the subset.
"""
from __future__ import annotations

import gzip
import json
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
RAW = Path("/tmp/concordant4_raw")
OUT = Path("/tmp/cldn4_vko")
QC_MIN_GENES = 200
QC_MIN_UMI = 500
QC_MAX_MT = 20.0
ELIG_131907 = {"tLung", "tL/B", "mLN", "PE", "mBrain"}


def load_universe() -> list[str]:
    programs = json.loads((ROOT / "data" / "programs.json").read_text())
    tfs = [
        g.strip().upper()
        for g in (ROOT / "data" / "human_tfs.txt").read_text().splitlines()
        if g.strip()
    ]
    genes = set(tfs)
    genes.add("CLDN4")
    genes.add("PVRL2")
    for key, val in programs.items():
        if isinstance(val, list):
            genes.update(g.upper() for g in val)
    for g in ["EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]:
        genes.add(g)
    return sorted(genes)


def save_cohort(dataset: str, genes: list[str], chunks: list[tuple[np.ndarray, np.ndarray, np.ndarray]]) -> None:
    """chunks: (counts cells x genes float32, patient str array, n_count float64)."""
    OUT.mkdir(parents=True, exist_ok=True)
    if not chunks:
        raise RuntimeError(f"no cells for {dataset}")
    mats, pats, ncounts = zip(*chunks)
    X = sparse.vstack([sparse.csr_matrix(m) for m in mats], format="csr")
    patient = np.concatenate(pats)
    n_count = np.concatenate(ncounts).astype(np.float64)
    sparse.save_npz(OUT / f"{dataset}.npz", X)
    np.savez(
        OUT / f"{dataset}.meta.npz",
        genes=np.array(genes),
        patient=patient,
        n_count=n_count,
        dataset=np.array(dataset),
    )
    print(
        f"saved {dataset}: cells={X.shape[0]} genes={X.shape[1]} "
        f"nnz={X.nnz} patients={len(set(patient.tolist()))}",
        flush=True,
    )


def qc_mask(n_count: np.ndarray, n_feat: np.ndarray, n_mt: np.ndarray) -> np.ndarray:
    mt_pct = np.divide(100.0 * n_mt, np.maximum(n_count, 1.0))
    return (n_feat >= QC_MIN_GENES) & (n_count >= QC_MIN_UMI) & (mt_pct < QC_MAX_MT)


def align_genes(symbols: list[str], universe: list[str]) -> dict[str, int]:
    """Map upper symbol -> column in the universe. First occurrence wins."""
    want = set(universe)
    out = {}
    for i, s in enumerate(symbols):
        u = s.upper()
        if u in want and u not in out:
            out[u] = i
    return out


def matrix_from_full(
    full: sparse.spmatrix,
    symbols: list[str],
    universe: list[str],
    n_count: np.ndarray,
    n_feat: np.ndarray,
    n_mt: np.ndarray,
    malignant: np.ndarray,
    patients: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    keep = malignant & qc_mask(n_count, n_feat, n_mt)
    if keep.sum() == 0:
        return None
    mapping = align_genes(symbols, universe)
    cols = []
    for g in universe:
        src = mapping.get(g)
        if src is None:
            cols.append(np.zeros(int(keep.sum()), dtype=np.float32))
        else:
            cols.append(np.asarray(full[src, keep].todense()).ravel().astype(np.float32))
    X = np.column_stack(cols)
    return X, patients[keep].astype(str), n_count[keep].astype(np.float64)


def extract_gse123902(universe: list[str], rows: list[dict]) -> None:
    locked = pd.read_csv(ROOT / "data" / "GSE123902_marker_units.tsv", sep="\t")
    tumor = locked[locked["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient")
    tumor = tumor[(tumor["n_malignant"] >= 20) & (tumor["n_tnk"] >= 20)]
    tar_path = RAW / "GSE123902" / "GSE123902_RAW.tar"
    extract_dir = Path("/tmp/gse123902_files")
    extract_dir.mkdir(exist_ok=True)
    chunks = []
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers()}
        for rec in tumor.itertuples(index=False):
            member = members.get(rec.file)
            if member is None:
                rows.append(_row("GSE123902", rec.patient, "donor", rec.tissue, 0, 0, 0, "missing_file"))
                continue
            dest = extract_dir / rec.file
            if not dest.exists():
                tf.extract(member, path=extract_dir)
                # tar extract keeps the name at the top
                extracted = extract_dir / rec.file
                if not extracted.exists():
                    # extract put it in extract_dir/name
                    pass
            df = pd.read_csv(extract_dir / rec.file, index_col=0)
            symbols = [str(c) for c in df.columns]
            # cells x genes -> genes x cells sparse
            mat = sparse.csr_matrix(df.to_numpy(dtype=np.float32).T)
            n_count = np.asarray(mat.sum(axis=1)).ravel() if False else np.asarray(mat.sum(axis=0)).ravel()
            n_feat = np.asarray((mat > 0).sum(axis=0)).ravel()
            mt = np.array([s.upper().startswith("MT-") for s in symbols])
            n_mt = np.asarray(mat[mt].sum(axis=0)).ravel() if mt.any() else np.zeros(mat.shape[1])
            sym_u = [s.upper() for s in symbols]
            index = {g: i for i, g in enumerate(sym_u) if g not in {}}  # first wins below
            first = {}
            for i, g in enumerate(sym_u):
                if g not in first:
                    first[g] = i

            def gene_row(name: str) -> np.ndarray:
                i = first.get(name)
                if i is None:
                    return np.zeros(mat.shape[1], dtype=np.float32)
                return np.asarray(mat[i].todense()).ravel()

            mal = (
                ((gene_row("EPCAM") > 0) | (gene_row("KRT8") > 0) | (gene_row("KRT18") > 0) | (gene_row("KRT19") > 0))
                & (gene_row("PTPRC") == 0)
            )
            tnk = (
                (gene_row("CD3D") > 0)
                | (gene_row("CD3E") > 0)
                | (gene_row("CD8A") > 0)
                | (gene_row("NKG7") > 0)
                | (gene_row("GNLY") > 0)
                | (gene_row("KLRD1") > 0)
            ) & ~mal
            patients = np.array([rec.patient] * mat.shape[1])
            packed = matrix_from_full(mat, symbols, universe, n_count, n_feat, n_mt, mal, patients)
            n_qc = 0 if packed is None else packed[0].shape[0]
            cldn = gene_row("CLDN4")
            mal_qc = mal & qc_mask(n_count, n_feat, n_mt)
            pct = float(100.0 * np.mean(cldn[mal_qc] > 0)) if mal_qc.any() else float("nan")
            rows.append(
                _row(
                    "GSE123902",
                    rec.patient,
                    "donor",
                    rec.tissue,
                    int(mat.shape[1]),
                    int(mal.sum()),
                    n_qc,
                    "ok" if n_qc else "no_qc_malignant",
                    n_tnk=int(tnk.sum()),
                    cldn4_pct=pct,
                    histology="LUAD",
                )
            )
            if packed is not None:
                chunks.append(packed)
            print(f"  GSE123902 {rec.patient} mal={int(mal.sum())} qc={n_qc}", flush=True)
            del df, mat
    save_cohort("GSE123902", universe, chunks)


def _row(dataset, unit_id, unit_type, tissue, n_cells, n_mal, n_qc, status, n_tnk=0, cldn4_pct=float("nan"), histology=""):
    return {
        "dataset": dataset,
        "unit_id": unit_id,
        "unit_type": unit_type,
        "tissue": tissue,
        "histology": histology,
        "n_cells": n_cells,
        "n_malignant": n_mal,
        "n_malignant_qc": n_qc,
        "n_tnk": n_tnk,
        "cldn4_pct_qc": cldn4_pct,
        "status": status,
    }


def extract_gse189357(universe: list[str], rows: list[dict]) -> None:
    locked = pd.read_csv(ROOT / "data" / "GSE189357_marker_units.tsv", sep="\t")
    tar_path = RAW / "GSE189357" / "GSE189357_RAW.tar"
    extract_dir = Path("/tmp/gse189357_files")
    extract_dir.mkdir(exist_ok=True)
    with tarfile.open(tar_path) as tf:
        names = tf.getnames()
        # extract once
        if not any(extract_dir.iterdir()):
            tf.extractall(extract_dir)
    chunks = []
    for rec in locked.itertuples(index=False):
        pat = rec.patient
        mtx = list(extract_dir.glob(f"*_{pat}_matrix.mtx.gz"))
        feat = list(extract_dir.glob(f"*_{pat}_features.tsv.gz"))
        bc = list(extract_dir.glob(f"*_{pat}_barcodes.tsv.gz"))
        if not mtx or not feat or not bc:
            rows.append(_row("GSE189357", pat, "patient", "TUMOR", 0, 0, 0, "missing_10x"))
            continue
        mat = spio.mmread(gzip.open(mtx[0], "rb")).tocsr()
        symbols = []
        with gzip.open(feat[0], "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                symbols.append(parts[1] if len(parts) > 1 else parts[0])
        if mat.shape[0] != len(symbols):
            raise RuntimeError(f"{pat} feature/matrix mismatch {mat.shape} vs {len(symbols)}")
        # collapse duplicate symbols by summing into the first row index
        upper = [s.upper() for s in symbols]
        n_count = np.asarray(mat.sum(axis=0)).ravel()
        n_feat = np.asarray((mat > 0).sum(axis=0)).ravel()
        mt = np.array([s.startswith("MT-") for s in upper])
        n_mt = np.asarray(mat[mt].sum(axis=0)).ravel() if mt.any() else np.zeros(mat.shape[1])
        first: dict[str, int] = {}
        for i, g in enumerate(upper):
            if g not in first:
                first[g] = i

        def gene_row(name: str) -> np.ndarray:
            i = first.get(name)
            if i is None:
                return np.zeros(mat.shape[1], dtype=np.float32)
            return np.asarray(mat[i].todense()).ravel().astype(np.float32)

        mal = (
            ((gene_row("EPCAM") > 0) | (gene_row("KRT8") > 0) | (gene_row("KRT18") > 0) | (gene_row("KRT19") > 0))
            & (gene_row("PTPRC") == 0)
        )
        tnk = (
            (gene_row("CD3D") > 0)
            | (gene_row("CD3E") > 0)
            | (gene_row("CD8A") > 0)
            | (gene_row("NKG7") > 0)
            | (gene_row("GNLY") > 0)
            | (gene_row("KLRD1") > 0)
        ) & ~mal
        patients = np.array([pat] * mat.shape[1])
        # rebuild a matrix with collapsed symbols only for universe genes inside matrix_from_full
        packed = matrix_from_full(mat, symbols, universe, n_count, n_feat, n_mt, mal, patients)
        n_qc = 0 if packed is None else int(packed[0].shape[0])
        mal_qc = mal & qc_mask(n_count, n_feat, n_mt)
        cldn = gene_row("CLDN4")
        pct = float(100.0 * np.mean(cldn[mal_qc] > 0)) if mal_qc.any() else float("nan")
        rows.append(
            _row(
                "GSE189357",
                pat,
                "patient",
                "TUMOR",
                int(mat.shape[1]),
                int(mal.sum()),
                n_qc,
                "ok" if n_qc else "no_qc_malignant",
                n_tnk=int(tnk.sum()),
                cldn4_pct=pct,
                histology="NSCLC",
            )
        )
        if packed is not None:
            chunks.append(packed)
        print(f"  GSE189357 {pat} mal={int(mal.sum())} qc={n_qc}", flush=True)
        del mat
    save_cohort("GSE189357", universe, chunks)


def extract_gse131907(universe: list[str], rows: list[dict]) -> None:
    ann_path = RAW / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = RAW / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    cells = []
    with gzip.open(ann_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            cells.append(
                {
                    "index": p[idx["Index"]],
                    "sample": p[idx["Sample"]],
                    "origin": p[idx["Sample_Origin"]],
                    "malignant": p[idx["Cell_subtype"]] == "Malignant cells",
                    "tnk": p[idx["Cell_type"]] in {"T lymphocytes", "NK cells"},
                }
            )
    by_sample: dict[str, list[dict]] = defaultdict(list)
    for rec in cells:
        by_sample[rec["sample"]].append(rec)
    eligible = []
    for sample, recs in sorted(by_sample.items()):
        origin = recs[0]["origin"]
        n_mal = sum(r["malignant"] for r in recs)
        n_tnk = sum(r["tnk"] for r in recs)
        ok = origin in ELIG_131907 and n_mal >= 20
        if ok:
            eligible.append(sample)
        rows.append(
            _row(
                "GSE131907",
                sample,
                "sample",
                origin,
                len(recs),
                n_mal,
                0,
                "pending" if ok else "not_eligible_origin_or_n",
                n_tnk=n_tnk,
                histology="LUAD",
            )
        )
    elig = set(eligible)
    mal_ids = [r["index"] for r in cells if r["malignant"] and r["sample"] in elig]
    sample_of = {r["index"]: r["sample"] for r in cells}
    print(f"  GSE131907 eligible samples={len(elig)} malignant cells={len(mal_ids)}", flush=True)

    universe_set = set(universe)
    gene_to_row = {g: i for i, g in enumerate(universe)}
    # COO builders
    coo_r: list[int] = []
    coo_c: list[int] = []
    coo_v: list[float] = []
    n_count = np.zeros(len(mal_ids), dtype=np.float64)
    n_feat = np.zeros(len(mal_ids), dtype=np.int32)
    n_mt = np.zeros(len(mal_ids), dtype=np.float64)

    with gzip.open(mat_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        col_of = {b: i for i, b in enumerate(barcodes)}
        keep_cols = [col_of[b] for b in mal_ids if b in col_of]
        if len(keep_cols) != len(mal_ids):
            missing = len(mal_ids) - len(keep_cols)
            raise RuntimeError(f"GSE131907 missing {missing} malignant barcodes in the matrix header")
        # local cell index for each matrix column
        targets = np.array(sorted(keep_cols), dtype=np.int32)
        col_to_local = {int(c): i for i, c in enumerate(targets)}
        # the mal_ids order must match local indices used below
        # rebuild mal_ids in target-column order
        bc_of_col = barcodes
        ordered_ids = [bc_of_col[int(c)] for c in targets]
        patients = np.array([sample_of[b] for b in ordered_ids])

        n_genes_seen = 0
        seen_genes: set[str] = set()
        for line in fh:
            line = line.rstrip("\n")
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab].upper()
            if gene in seen_genes:
                n_genes_seen += 1
                continue
            seen_genes.add(gene)
            store_row = gene_to_row.get(gene) if gene in universe_set else None
            is_mt = gene.startswith("MT-")
            # walk fields; values start at field 0
            field = 0
            i = tab + 1
            n = len(line)
            ti = 0
            n_t = len(targets)
            while ti < n_t and i < n:
                want = int(targets[ti])
                while field < want:
                    j = line.find("\t", i)
                    if j < 0:
                        i = n
                        break
                    i = j + 1
                    field += 1
                if field != want or i >= n:
                    break
                j = line.find("\t", i)
                if j < 0:
                    j = n
                token = line[i:j]
                if token and token not in {"0", "0.0"}:
                    val = float(token)
                    if val != 0.0:
                        local = col_to_local[want]
                        n_count[local] += val
                        n_feat[local] += 1
                        if is_mt:
                            n_mt[local] += val
                        if store_row is not None:
                            coo_r.append(store_row)
                            coo_c.append(local)
                            coo_v.append(val)
                i = j + 1
                field += 1
                ti += 1
            n_genes_seen += 1
            if n_genes_seen % 4000 == 0:
                print(f"    streamed {n_genes_seen} genes nnz={len(coo_v)}", flush=True)

    X = sparse.coo_matrix(
        (np.asarray(coo_v, dtype=np.float32), (np.asarray(coo_r), np.asarray(coo_c))),
        shape=(len(universe), len(targets)),
    ).tocsr()
    # cells x genes
    Xt = X.T.tocsr()
    keep = qc_mask(n_count, n_feat, n_mt)
    print(f"  GSE131907 qc {int(keep.sum())} / {len(targets)}", flush=True)
    # CLDN4 pct per sample on QC cells
    cldn_i = universe.index("CLDN4")
    cldn = np.asarray(Xt.getcol(cldn_i).todense()).ravel()
    by = defaultdict(lambda: [0, 0])
    for local, ok in enumerate(keep):
        if not ok:
            continue
        s = patients[local]
        by[s][1] += 1
        if cldn[local] > 0:
            by[s][0] += 1
    for row in rows:
        if row["dataset"] != "GSE131907" or row["status"] != "pending":
            continue
        pos, n = by.get(row["unit_id"], [0, 0])
        row["n_malignant_qc"] = n
        row["cldn4_pct_qc"] = (100.0 * pos / n) if n else float("nan")
        row["status"] = "ok" if n else "no_qc_malignant"
    save_cohort(
        "GSE131907",
        universe,
        [(Xt[keep].toarray().astype(np.float32), patients[keep], n_count[keep])],
    )


def extract_gse205335(universe: list[str], rows: list[dict]) -> None:
    import subprocess

    ident = pd.read_csv(RAW / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    soft_text = gzip.open(RAW / "GSE205335" / "GSE205335_family.soft.gz", "rt").read().splitlines()
    soft = []
    cur: dict[str, str] = {}
    for line in soft_text:
        if line.startswith("^SAMPLE"):
            if cur.get("title"):
                soft.append(cur)
            cur = {}
        elif line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1"):
            val = line.split("=", 1)[1].strip()
            if ": " in val:
                k, v = val.split(": ", 1)
                cur[k] = v
    if cur.get("title"):
        soft.append(cur)
    code_map = {}
    for s in soft:
        title = s.get("title", "")
        code = title.split(" ", 1)[1] if " " in title else ""
        code = code.upper().replace("-", "_")
        code_map[code] = s
    locked = pd.read_csv(ROOT / "data" / "GSE205335_patients.tsv", sep="\t")
    keep_pats = set(locked["patient"].astype(str))
    hist = dict(zip(locked["patient"].astype(str), locked["cancer_subtype"].astype(str)))

    def norm_code(x: str) -> str:
        x = str(x).upper().replace("-", "_")
        if x.endswith("_3P") or x.endswith("_5P"):
            x = x[:-3]
        return x

    ident["code"] = ident["orig.ident"].map(norm_code)
    ident["patient"] = ident["code"].map(lambda c: code_map.get(c, {}).get("patient", ""))
    ident["tissue"] = ident["code"].map(lambda c: code_map.get(c, {}).get("tissue", ""))
    ident = ident[ident["patient"].isin(keep_pats)].copy()
    ident["is_normal"] = ident["tissue"].str.startswith("Normal", na=False)
    tumor = ident[~ident["is_normal"]].copy()
    mal = tumor[tumor["lineage.sub"] == "Malignant cells"].copy()
    print(
        f"  GSE205335 ident matched {len(ident)} tumor {len(tumor)} malignant {len(mal)}",
        flush=True,
    )
    if len(mal) < 1000:
        raise RuntimeError("GSE205335 malignant barcode map failed")

    rds = Path("/tmp/gse205335_matrix.rds")
    if not rds.exists():
        raise RuntimeError("missing inflated RDS /tmp/gse205335_matrix.rds")
    work = Path("/tmp/gse205335_subset")
    work.mkdir(exist_ok=True)
    bc_file = work / "barcodes.tsv"
    mal[["barcode", "patient"]].to_csv(bc_file, sep="\t", index=False)
    gene_file = work / "genes.txt"
    gene_file.write_text("\n".join(universe) + "\n")
    subprocess.check_call(
        ["Rscript", str(ROOT / "scripts" / "subset_gse205335.R"), str(rds), str(bc_file), str(gene_file), str(work)]
    )
    qc = pd.read_csv(work / "cell_qc.tsv", sep="\t")
    genes = [g.strip() for g in (work / "genes.txt").read_text().splitlines() if g.strip()]
    mtx = spio.mmread(work / "matrix.mtx").tocsr()  # genes x cells
    if mtx.shape[1] != len(qc) or mtx.shape[0] != len(genes):
        raise RuntimeError(f"GSE205335 mtx shape {mtx.shape} vs genes {len(genes)} cells {len(qc)}")
    n_count = qc["n_count"].to_numpy(dtype=np.float64)
    n_feat = qc["n_feat"].to_numpy(dtype=np.float64)
    n_mt = qc["n_mt"].to_numpy(dtype=np.float64)
    keep = qc_mask(n_count, n_feat, n_mt)
    # place genes into universe columns
    present = {g: i for i, g in enumerate(genes)}
    # build cells x universe on the QC subset only, per patient chunks to limit RAM
    patients = qc["patient"].astype(str).to_numpy()
    Xt = mtx.T.tocsr()
    chunks = []
    for pat in sorted(keep_pats):
        sel = (patients == pat) & keep
        n_mal_pat = int((patients == pat).sum())
        n_qc = int(sel.sum())
        tissue = ",".join(sorted(set(tumor.loc[tumor["patient"] == pat, "tissue"].astype(str))))
        n_tnk = int(((tumor["patient"] == pat) & (tumor["lineage.total"] == "T/NK cells")).sum())
        n_cells = int((tumor["patient"] == pat).sum())
        if n_qc == 0:
            rows.append(
                _row(
                    "GSE205335",
                    pat,
                    "patient",
                    tissue,
                    n_cells,
                    n_mal_pat,
                    0,
                    "no_qc_malignant",
                    n_tnk=n_tnk,
                    histology=hist.get(pat, ""),
                )
            )
            continue
        block = np.zeros((n_qc, len(universe)), dtype=np.float32)
        sub = Xt[sel]
        for j, g in enumerate(universe):
            src = present.get(g)
            if src is not None:
                block[:, j] = np.asarray(sub.getcol(src).todense()).ravel()
        cldn_j = universe.index("CLDN4")
        pct = float(100.0 * np.mean(block[:, cldn_j] > 0))
        rows.append(
            _row(
                "GSE205335",
                pat,
                "patient",
                tissue,
                n_cells,
                n_mal_pat,
                n_qc,
                "ok",
                n_tnk=n_tnk,
                cldn4_pct=pct,
                histology=hist.get(pat, ""),
            )
        )
        chunks.append((block, np.array([pat] * n_qc), n_count[sel]))
        print(f"  GSE205335 {pat} mal={n_mal_pat} qc={n_qc} cldn4%={pct:.1f}", flush=True)
    save_cohort("GSE205335", universe, chunks)


def main() -> None:
    universe = load_universe()
    print(f"universe genes {len(universe)}", flush=True)
    rows: list[dict] = []
    extract_gse123902(universe, rows)
    extract_gse189357(universe, rows)
    extract_gse205335(universe, rows)
    extract_gse131907(universe, rows)
    inv = pd.DataFrame(rows)
    out_tab = ROOT / "results" / "tables"
    out_tab.mkdir(parents=True, exist_ok=True)
    inv.to_csv(out_tab / "unit_inventory.tsv", sep="\t", index=False)
    print(inv.groupby("dataset")[["n_malignant", "n_malignant_qc"]].sum())
    print("inventory", out_tab / "unit_inventory.tsv")


if __name__ == "__main__":
    main()
