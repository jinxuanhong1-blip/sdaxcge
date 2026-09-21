#!/usr/bin/env python3
"""Extract locked-unit malignant counts for a fixed gene panel.

Compartments match the locked concordant-4 Seurat run. Counts are checked
against data/locked_units.tsv before anything is written.
"""

from __future__ import annotations

import gzip
import os
import re
import subprocess
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
LOCKED_PATH = ROOT / "data" / "locked_units.tsv"

PANEL = [
    "ELF3", "TACSTD2", "CLDN4", "CLDN7",
    "CDH1", "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT5",
    "GRHL2", "GRHL1", "GRHL3", "CLDN1", "CLDN3", "CDH3", "MUC1",
    "KLF5", "TJP1", "F11R", "SPDEF", "OCLN", "OVOL2",
]
MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
GENES = list(dict.fromkeys(PANEL + MARKERS))

GSE123902_FILES = {
    "LX255B": "GSM3516668_MSK_LX255B_METASTASIS_dense.csv.gz",
    "LX653": "GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz",
    "LX661": "GSM3516663_MSK_LX661_PRIMARY_TUMOUR_dense.csv.gz",
    "LX666": "GSM3516664_MSK_LX666_METASTASIS_dense.csv.gz",
    "LX675": "GSM3516665_MSK_LX675_PRIMARY_TUMOUR_dense.csv.gz",
    "LX676": "GSM3516667_MSK_LX676_PRIMARY_TUMOUR_dense.csv.gz",
    "LX679": "GSM3516669_MSK_LX679_PRIMARY_TUMOUR_dense.csv.gz",
    "LX680": "GSM3516670_MSK_LX680_PRIMARY_TUMOUR_dense.csv.gz",
    "LX681": "GSM3516671_MSK_LX681_METASTASIS_dense.csv.gz",
    "LX682": "GSM3516672_MSK_LX682_PRIMARY_TUMOUR_dense.csv.gz",
    "LX684": "GSM3516674_MSK_LX684_PRIMARY_TUMOUR_dense.csv.gz",
    "LX699": "GSM3516677_MSK_LX699_METASTASIS_dense.csv.gz",
    "LX701": "GSM3516678_MSK_LX701_METASTASIS_dense.csv.gz",
}


def say(msg: str) -> None:
    print(msg, flush=True)


def load_locked(dataset: str) -> pd.DataFrame:
    df = pd.read_csv(LOCKED_PATH, sep="\t", dtype={"unit_id": str})
    df = df[df.dataset == dataset].copy()
    for col in ["n_cells", "n_malignant", "n_tnk"]:
        df[col] = df[col].astype(int)
    return df


def check_units(dataset: str, got: pd.DataFrame) -> None:
    exp = load_locked(dataset).set_index("unit_id")
    g = got.set_index("unit_id")
    missing = [u for u in exp.index if u not in g.index]
    if missing:
        raise SystemExit(f"{dataset} missing units {missing}")
    g = g.loc[exp.index]
    for col in ["n_cells", "n_malignant", "n_tnk"]:
        if not np.array_equal(g[col].to_numpy(), exp[col].to_numpy()):
            bad = g.index[g[col].to_numpy() != exp[col].to_numpy()]
            raise SystemExit(f"{dataset} {col} mismatch {list(bad)}")
    say(f"locked counts ok {dataset} n={len(exp)}")


def save_panel(
    out: Path,
    dataset: str,
    genes: list[str],
    unit_id: np.ndarray,
    counts: np.ndarray,
    lib: np.ndarray,
    measured_units: list[str],
    measured: np.ndarray,
) -> None:
    out.mkdir(parents=True, exist_ok=True)
    counts = np.nan_to_num(counts, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    np.savez_compressed(
        out / f"{dataset}.npz",
        genes=np.array(genes, dtype=str),
        unit_id=np.asarray(unit_id, dtype=str),
        counts=counts,
        lib=np.asarray(lib, dtype=np.float32),
        measured_units=np.array(measured_units, dtype=str),
        measured=np.asarray(measured, dtype=bool),
    )
    say(f"wrote {dataset} cells={counts.shape[0]} genes={counts.shape[1]}")


def align_genes(present: list[str], mat: np.ndarray) -> np.ndarray:
    """mat columns follow `present`. Return columns in GENES order, 0 if absent."""
    idx = {g: i for i, g in enumerate(present)}
    out = np.zeros((mat.shape[0], len(GENES)), dtype=np.float32)
    for j, g in enumerate(GENES):
        i = idx.get(g)
        if i is not None:
            out[:, j] = mat[:, i]
    return out


def marker_masks(get, n: int):
    def g(name):
        v = get(name)
        if v is None:
            return np.zeros(n, dtype=np.float32)
        return np.asarray(v, dtype=np.float32)

    mal = ((g("EPCAM") > 0) | (g("KRT8") > 0) | (g("KRT18") > 0) | (g("KRT19") > 0)) & (g("PTPRC") == 0)
    tnk = (
        (g("CD3D") > 0) | (g("CD3E") > 0) | (g("CD8A") > 0) | (g("NKG7") > 0) | (g("GNLY") > 0) | (g("KLRD1") > 0)
    ) & ~mal
    return mal, tnk


def _gene_matrix(df: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    cols = {c.upper(): c for c in df.columns}
    present = [g for g in genes if g in cols]
    mat = np.column_stack([df[cols[g]].to_numpy(dtype=np.float32) for g in present]) if present else np.zeros((len(df), 0), np.float32)
    return mat, present


def extract_gse123902(geo: Path, out: Path) -> None:
    dataset = "GSE123902"
    locked = load_locked(dataset)
    ids = locked.unit_id.tolist()
    say(f"== {dataset} ==")
    unit_rows = []
    blocks_u, blocks_c, blocks_lib = [], [], []
    measured_rows = []
    with tarfile.open(geo / "gse123902" / "GSE123902_RAW.tar") as tar:
        for unit in ids:
            say(f"  {unit}")
            raw = tar.extractfile(GSE123902_FILES[unit])
            with gzip.GzipFile(fileobj=raw) as gz:
                df = pd.read_csv(gz, index_col=0)
            df.columns = [str(c).upper() for c in df.columns]
            df = df.loc[:, ~df.columns.duplicated()]
            n = len(df)

            def getter(name, _df=df):
                if name not in _df.columns:
                    return None
                return _df[name].to_numpy(dtype=np.float32)

            mal, tnk = marker_masks(getter, n)
            unit_rows.append(dict(
                unit_id=unit, n_cells=int(n), n_malignant=int(mal.sum()), n_tnk=int(tnk.sum()),
            ))
            mal_df = df.loc[mal]
            mat, present = _gene_matrix(mal_df, GENES)
            measured_rows.append([g in set(present) for g in GENES])
            blocks_u.append(np.full(mal_df.shape[0], unit, dtype=object))
            blocks_c.append(align_genes(present, mat))
            blocks_lib.append(df.to_numpy(dtype=np.float32).sum(axis=1)[mal])
            del df, mal_df
    check_units(dataset, pd.DataFrame(unit_rows))
    save_panel(
        out, dataset, GENES, np.concatenate(blocks_u), np.vstack(blocks_c), np.concatenate(blocks_lib),
        ids, np.asarray(measured_rows, dtype=bool),
    )


def extract_gse189357(geo: Path, out: Path) -> None:
    dataset = "GSE189357"
    locked = load_locked(dataset)
    ids = locked.unit_id.tolist()
    raw_dir = geo / "gse189357" / "raw"
    if not any(raw_dir.glob("*_matrix.mtx.gz")):
        raw_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(geo / "gse189357" / "GSE189357_RAW.tar") as tar:
            tar.extractall(raw_dir)
    say(f"== {dataset} ==")
    unit_rows = []
    blocks_u, blocks_c, blocks_lib = [], [], []
    measured_rows = []
    for unit in ids:
        say(f"  {unit}")
        mtx = next(raw_dir.glob(f"*_{unit}_matrix.mtx.gz"))
        feat = next(raw_dir.glob(f"*_{unit}_features.tsv.gz"))
        symbols = []
        with gzip.open(feat, "rt") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                symbols.append(parts[1].upper() if len(parts) > 1 else "")
        with gzip.open(mtx, "rb") as handle:
            mat = mmread(handle).tocsr()
        if mat.shape[0] != len(symbols):
            raise SystemExit(f"{unit} features {len(symbols)} != matrix {mat.shape[0]}")
        first = {}
        for i, g in enumerate(symbols):
            if g and g not in first:
                first[g] = i
        n = mat.shape[1]

        def getter(name, _first=first, _mat=mat):
            i = _first.get(name)
            if i is None:
                return None
            return np.asarray(_mat.getrow(i).toarray()).ravel()

        mal, tnk = marker_masks(getter, n)
        unit_rows.append(dict(unit_id=unit, n_cells=int(n), n_malignant=int(mal.sum()), n_tnk=int(tnk.sum())))
        mal_idx = np.flatnonzero(mal)
        present = [g for g in GENES if g in first]
        measured_rows.append([g in first for g in GENES])
        raw = np.vstack([np.asarray(mat.getrow(first[g]).toarray()).ravel()[mal_idx] for g in present]).T if present else np.zeros((mal_idx.size, 0), np.float32)
        lib = np.asarray(mat.sum(axis=0)).ravel()[mal_idx]
        blocks_u.append(np.full(mal_idx.size, unit, dtype=object))
        blocks_c.append(align_genes(present, raw))
        blocks_lib.append(lib)
        del mat
    check_units(dataset, pd.DataFrame(unit_rows))
    save_panel(
        out, dataset, GENES, np.concatenate(blocks_u), np.vstack(blocks_c), np.concatenate(blocks_lib),
        ids, np.asarray(measured_rows, dtype=bool),
    )


def extract_gse131907(geo: Path, out: Path) -> None:
    dataset = "GSE131907"
    locked = load_locked(dataset)
    ids = locked.unit_id.tolist()
    say(f"== {dataset} ==")
    ann = pd.read_csv(geo / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str)
    ann = ann[ann["Sample"].isin(ids)]
    n_cells = ann.groupby("Sample").size().reindex(ids).astype(int)
    n_mal = ann["Cell_subtype"].eq("Malignant cells").groupby(ann["Sample"]).sum().reindex(ids).astype(int)
    n_tnk = ann["Cell_type"].isin(["T lymphocytes", "NK cells"]).groupby(ann["Sample"]).sum().reindex(ids).astype(int)
    got = pd.DataFrame({"unit_id": ids, "n_cells": n_cells.to_numpy(), "n_malignant": n_mal.to_numpy(), "n_tnk": n_tnk.to_numpy()})
    check_units(dataset, got)

    mat_path = geo / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    with gzip.open(mat_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    barcodes = header[1:]
    meta = ann.set_index("Index").reindex(barcodes)
    keep = meta["Sample"].isin(ids).to_numpy() & meta["Cell_subtype"].eq("Malignant cells").to_numpy()
    keep = np.where(meta["Sample"].notna().to_numpy(), keep, False)
    keep_idx = np.flatnonzero(keep)
    unit_of = meta["Sample"].to_numpy()[keep_idx]
    order = np.argsort(unit_of, kind="mergesort")
    unit_of = unit_of[order]
    keep_idx = keep_idx[order]
    if keep_idx.size != int(n_mal.sum()):
        raise SystemExit(f"GSE131907 malignant cells {keep_idx.size} != annotation {int(n_mal.sum())}")
    say(f"  streaming malignant cells={keep_idx.size} genes-of-interest={len(GENES)}")
    stored = {}
    lib = np.zeros(keep_idx.size, dtype=np.float64)
    seen = set()
    with gzip.open(mat_path, "rt") as handle:
        handle.readline()
        n_genes = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.upper()
            if not gene or gene in seen:
                continue
            seen.add(gene)
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != len(barcodes):
                raise SystemExit(f"{gene}: {arr.size} != {len(barcodes)}")
            kept = arr[keep_idx]
            lib += kept
            if gene in GENES:
                stored[gene] = kept
            n_genes += 1
            if n_genes % 5000 == 0:
                say(f"  genes={n_genes} stored={len(stored)}")
    say(f"  stream done genes={n_genes} stored={len(stored)}")
    present = [g for g in GENES if g in stored]
    raw = np.column_stack([stored[g] for g in present]) if present else np.zeros((unit_of.size, 0), np.float32)
    # A gene is measured in a unit if it occurred in the shared matrix (one feature list).
    measured = np.tile(np.array([g in stored for g in GENES], dtype=bool), (len(ids), 1))
    save_panel(out, dataset, GENES, unit_of, align_genes(present, raw), lib, ids, measured)


def build_gse205335_map(geo: Path) -> Path:
    ident = pd.read_csv(geo / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t", dtype=str)
    samples = []
    cur = {}
    with gzip.open(geo / "gse205335" / "GSE205335_family.soft.gz", "rt") as handle:
        lines = handle.read().splitlines()
    for line in lines:
        if line.startswith("^SAMPLE"):
            if cur.get("title"):
                samples.append(cur)
            cur = {}
        elif line.startswith("!Sample_title"):
            cur["title"] = line.split(" = ", 1)[1]
        elif line.startswith("!Sample_characteristics_ch1"):
            val = line.split(" = ", 1)[1]
            if ": " in val:
                k, v = val.split(": ", 1)
                cur[k] = v
    if cur.get("title"):
        samples.append(cur)
    rows = []
    for s in samples:
        title = s.get("title", "")
        code = title.split(" ", 1)[1] if " " in title else ""
        rows.append({
            "patient": s.get("patient"),
            "tissue": s.get("tissue", ""),
            "code": code.upper().replace("-", "_"),
        })
    soft = pd.DataFrame(rows)
    if soft["code"].duplicated().any():
        raise SystemExit("GSE205335 SOFT codes are not unique")

    def norm_code(x: str) -> str:
        x = str(x).upper().replace("-", "_")
        return re.sub(r"_[35]P$", "", x)

    ident = ident.copy()
    ident["orig_code"] = ident["orig.ident"].map(norm_code)
    m = ident.merge(soft, left_on="orig_code", right_on="code", how="left")
    if m["patient"].isna().any():
        raise SystemExit(f"GSE205335 unmapped cells {int(m['patient'].isna().sum())}")
    keep = set(load_locked("GSE205335").unit_id)
    m = m[m.patient.isin(keep)].copy()
    m["is_normal"] = m["tissue"].fillna("").str.startswith("Normal")
    m["is_mal"] = m["lineage.sub"].eq("Malignant cells") & ~m["is_normal"]
    m["is_tnk"] = m["lineage.total"].eq("T/NK cells") & ~m["is_normal"]
    out = geo / "gse205335" / "map.tsv"
    m[["barcode", "patient", "is_mal", "is_tnk", "is_normal"]].rename(columns={"patient": "unit_id"}).to_csv(
        out, sep="\t", index=False
    )
    say(f"wrote map {out} rows={len(m)}")
    return out


def peel_rds(gz_path: Path, rds_path: Path) -> None:
    if rds_path.exists() and rds_path.stat().st_size > 0:
        return
    say("peeling double-gzip RDS")
    with gzip.open(gz_path, "rb") as outer:
        magic = outer.read(2)
        outer.seek(0)
        if magic == b"\x1f\x8b":
            # gzip.open().mode is an int; pass mode so the inner layer is not inferred from it.
            with gzip.GzipFile(fileobj=outer, mode="rb") as inner, rds_path.open("wb") as dest:
                while True:
                    chunk = inner.read(16 * 1024 * 1024)
                    if not chunk:
                        break
                    dest.write(chunk)
        else:
            with rds_path.open("wb") as dest:
                while True:
                    chunk = outer.read(16 * 1024 * 1024)
                    if not chunk:
                        break
                    dest.write(chunk)


def extract_gse205335(geo: Path, out: Path) -> None:
    dataset = "GSE205335"
    say(f"== {dataset} ==")
    build_gse205335_map(geo)
    gz = geo / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    rds = geo / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds"
    peel_rds(gz, rds)
    genes = geo / "gse205335" / "panel_genes.txt"
    genes.write_text("\n".join(GENES) + "\n")
    tsv_gz = out / "GSE205335_panel.tsv.gz"
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update({
        "RDS": str(rds),
        "MAP": str(geo / "gse205335" / "map.tsv"),
        "OUT": str(tsv_gz),
        "GENES": str(genes),
    })
    subprocess.check_call(["Rscript", str(ROOT / "scripts" / "extract_gse205335.R")], env=env)
    say("reading panel table")
    df = pd.read_csv(tsv_gz, sep="\t", dtype={"unit_id": str})
    present = [c for c in df.columns if c not in ("unit_id", "lib")]
    raw = df[present].to_numpy(dtype=np.float32) if present else np.zeros((len(df), 0), np.float32)
    got = df.groupby("unit_id").size()
    locked = load_locked(dataset).set_index("unit_id")
    if not np.array_equal(got.reindex(locked.index).to_numpy(), locked.n_malignant.to_numpy()):
        bad = locked.index[got.reindex(locked.index).to_numpy() != locked.n_malignant.to_numpy()]
        raise SystemExit(f"GSE205335 malignant count mismatch {list(bad)}")
    say(f"locked malignant counts ok {dataset}")
    measured = np.tile(np.array([g in set(present) for g in GENES], dtype=bool), (len(locked), 1))
    save_panel(
        out, dataset, GENES, df.unit_id.to_numpy(), align_genes(present, raw),
        df["lib"].to_numpy(dtype=np.float32), list(locked.index), measured,
    )


def main() -> None:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--geo", type=Path, default=Path("/tmp/geo_c4"))
    p.add_argument("--out", type=Path, default=Path("/tmp/geo_c4/panel"))
    p.add_argument("--datasets", nargs="*", default=["GSE123902", "GSE189357", "GSE131907", "GSE205335"])
    args = p.parse_args()
    runners = {
        "GSE123902": extract_gse123902,
        "GSE189357": extract_gse189357,
        "GSE131907": extract_gse131907,
        "GSE205335": extract_gse205335,
    }
    for ds in args.datasets:
        runners[ds](args.geo, args.out)


if __name__ == "__main__":
    main()
