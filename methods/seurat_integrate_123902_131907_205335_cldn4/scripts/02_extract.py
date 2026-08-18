#!/usr/bin/env python3
"""Extract-only. Stream GSE131907 UMI TSV and GSE123902 dense CSVs to MTX.

Primary analysis is Seurat + Harmony in 03_seurat_integrate.R.
No Python-only primary. No dual-high. Not GSE148071 / GSE189357.
"""
from __future__ import annotations

import csv
import gzip
import random
import tarfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GEO = Path("/tmp/triple_geo")
EXT = ROOT / "extract"
DATA = ROOT / "data"

EPI_MARKERS = ("EPCAM", "KRT8", "KRT18", "KRT19")
TNK_MARKERS = ("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1")
CAP_MAL = 80
CAP_TNK = 80
SEED = 1


def write_mtx(path: Path, rows: np.ndarray, cols: np.ndarray, vals: np.ndarray,
              n_row: int, n_col: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("%%MatrixMarket matrix coordinate real general\n")
        fh.write(f"{n_row} {n_col} {len(vals)}\n")
        for r, c, v in zip(rows, cols, vals):
            fh.write(f"{int(r) + 1} {int(c) + 1} {float(v)}\n")


def extract_gse131907() -> None:
    keep_path = EXT / "keep_GSE131907.tsv"
    umi = GEO / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    out = EXT / "GSE131907"
    keep_ids: list[str] = []
    with keep_path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        bcol = header.index("barcode")
        for line in fh:
            keep_ids.append(line.rstrip("\n").split("\t")[bcol])
    keep_set = set(keep_ids)
    print(f"[131907] keep barcodes={len(keep_ids)} unique={len(keep_set)}", flush=True)

    with gzip.open(umi, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep_set]
        if not col_idx:
            raise SystemExit("no requested GSE131907 cell IDs in UMI header")
        ordered = [cell_ids[i] for i in col_idx]
        print(f"[131907] header cells={len(cell_ids)} keep={len(col_idx)}", flush=True)
        genes: list[str] = []
        all_rows: list[np.ndarray] = []
        all_cols: list[np.ndarray] = []
        all_vals: list[np.ndarray] = []
        nnz_total = 0
        for gi, line in enumerate(f, start=0):
            raw = line.rstrip("\n")
            if not raw:
                continue
            tab0 = raw.find("\t")
            gene = raw[:tab0]
            vals = np.fromstring(raw[tab0 + 1 :], sep="\t", dtype=np.float32)
            if vals.size != len(cell_ids):
                raise SystemExit(f"row {gi} {gene}: {vals.size} != {len(cell_ids)}")
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                all_rows.append(np.full(nz.size, gi, dtype=np.int32))
                all_cols.append(nz.astype(np.int32, copy=False))
                all_vals.append(sub[nz].astype(np.float32, copy=False))
                nnz_total += int(nz.size)
            genes.append(gene)
            if (gi + 1) % 2000 == 0:
                print(f"  streamed {gi + 1} genes nnz={nnz_total}", flush=True)

    if all_vals:
        rows = np.concatenate(all_rows)
        cols = np.concatenate(all_cols)
        vals = np.concatenate(all_vals)
    else:
        rows = np.array([], dtype=np.int32)
        cols = np.array([], dtype=np.int32)
        vals = np.array([], dtype=np.float32)
    write_mtx(out / "matrix.mtx", rows, cols, vals, len(genes), len(ordered))
    (out / "genes.tsv").write_text("\n".join(genes) + "\n")
    (out / "barcodes.tsv").write_text("\n".join(ordered) + "\n")
    print(f"[131907] wrote {out} genes={len(genes)} cells={len(ordered)} nnz={len(vals)}", flush=True)


def _load_locked_123902() -> list[dict]:
    path = DATA / "GSE123902_marker_units.tsv"
    rows = []
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for rec in rdr:
            rec["eligible"] = str(rec["eligible"]).lower() == "true"
            rows.append(rec)
    tumor = [r for r in rows if r["tissue"] in ("PRIMARY", "METASTASIS") and r["eligible"]]
    # PRIMARY preferred over METASTASIS
    tumor.sort(key=lambda r: (r["patient"], 0 if r["tissue"] == "PRIMARY" else 1))
    seen = set()
    keep = []
    for r in tumor:
        if r["patient"] in seen:
            continue
        seen.add(r["patient"])
        keep.append(r)
    return keep


def _read_dense_csv(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
        genes = [g.strip().strip('"') for g in header[1:]]
        barcodes: list[str] = []
        chunks: list[np.ndarray] = []
        for line in fh:
            raw = line.rstrip("\n")
            if not raw:
                continue
            comma = raw.find(",")
            barcodes.append(raw[:comma].strip().strip('"'))
            vals = np.fromstring(raw[comma + 1 :], sep=",", dtype=np.float32)
            if vals.size != len(genes):
                raise SystemExit(f"{path.name}: row {len(barcodes)} width {vals.size} != {len(genes)}")
            chunks.append(vals)
    mat = np.vstack(chunks) if chunks else np.zeros((0, len(genes)), dtype=np.float32)
    return barcodes, genes, mat


def extract_gse123902() -> None:
    csv_dir = GEO / "GSE123902" / "csv"
    if not any(csv_dir.glob("GSM*_dense.csv.gz")):
        tar = GEO / "GSE123902" / "GSE123902_RAW.tar"
        csv_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar) as tf:
            tf.extractall(csv_dir)
    units = _load_locked_123902()
    print(f"[123902] locked tumor donors={len(units)}", flush=True)
    rng = random.Random(SEED + 10)

    gene_union: list[str] = []
    gene_index: dict[str, int] = {}
    cell_barcodes: list[str] = []
    cell_meta: list[dict] = []
    # COO accumulators (gene, cell, value)
    coo_r: list[np.ndarray] = []
    coo_c: list[np.ndarray] = []
    coo_v: list[np.ndarray] = []

    for rec in units:
        fp = csv_dir / rec["file"]
        if not fp.exists():
            matches = list(csv_dir.glob(f"{rec['gsm']}_*_dense.csv.gz"))
            if not matches:
                raise SystemExit(f"missing {rec['file']}")
            fp = matches[0]
        print(f"  {fp.name}", flush=True)
        barcodes, genes, mat = _read_dense_csv(fp)
        gmap = {g.upper(): j for j, g in enumerate(genes)}

        def col_or_zero(name: str) -> np.ndarray:
            j = gmap.get(name)
            if j is None:
                return np.zeros(mat.shape[0], dtype=np.float32)
            return mat[:, j]

        epi = np.zeros(mat.shape[0], dtype=bool)
        for g in EPI_MARKERS:
            epi |= col_or_zero(g) > 0
        ptprc = col_or_zero("PTPRC")
        mal = epi & (ptprc == 0)
        tnk = np.zeros(mat.shape[0], dtype=bool)
        for g in TNK_MARKERS:
            tnk |= col_or_zero(g) > 0
        tnk = tnk & (~mal)

        mal_idx = np.flatnonzero(mal).tolist()
        tnk_idx = np.flatnonzero(tnk).tolist()
        if len(mal_idx) > CAP_MAL:
            mal_idx = rng.sample(mal_idx, CAP_MAL)
        if len(tnk_idx) > CAP_TNK:
            tnk_idx = rng.sample(tnk_idx, CAP_TNK)
        keep_idx = mal_idx + tnk_idx
        if not keep_idx:
            print(f"  skip {rec['patient']}: no marker mal/TNK", flush=True)
            continue
        keep_idx = sorted(keep_idx)
        sub = mat[keep_idx, :]
        patient = rec["patient"]
        for local_i, src_i in enumerate(keep_idx):
            bc = f"{rec['gsm']}_{barcodes[src_i]}"
            cell_barcodes.append(bc)
            cell_meta.append({
                "barcode": bc,
                "dataset": "GSE123902",
                "unit_id": patient,
                "patient_id": patient,
                "compartment": "malignant" if src_i in set(mal_idx) else "TNK",
                "gsm": rec["gsm"],
                "tissue": rec["tissue"],
                "file": rec["file"],
            })
            col = len(cell_barcodes) - 1
            nz = np.flatnonzero(sub[local_i])
            if not nz.size:
                continue
            mapped = []
            for j in nz:
                g = genes[int(j)]
                if g not in gene_index:
                    gene_index[g] = len(gene_union)
                    gene_union.append(g)
                mapped.append(gene_index[g])
            coo_r.append(np.asarray(mapped, dtype=np.int32))
            coo_c.append(np.full(len(mapped), col, dtype=np.int32))
            coo_v.append(sub[local_i, nz].astype(np.float32, copy=False))
        del mat, sub
    if not cell_barcodes:
        raise SystemExit("GSE123902 extract produced zero cells")
    rows = np.concatenate(coo_r) if coo_r else np.array([], dtype=np.int32)
    cols = np.concatenate(coo_c) if coo_c else np.array([], dtype=np.int32)
    vals = np.concatenate(coo_v) if coo_v else np.array([], dtype=np.float32)
    out = EXT / "GSE123902"
    write_mtx(out / "matrix.mtx", rows, cols, vals, len(gene_union), len(cell_barcodes))
    (out / "genes.tsv").write_text("\n".join(gene_union) + "\n")
    (out / "barcodes.tsv").write_text("\n".join(cell_barcodes) + "\n")
    keep_path = EXT / "keep_GSE123902.tsv"
    with keep_path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["barcode", "dataset", "unit_id", "patient_id",
                        "compartment", "gsm", "tissue", "file"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(cell_meta)
    print(
        f"[123902] wrote {out} genes={len(gene_union)} cells={len(cell_barcodes)} "
        f"nnz={len(vals)} patients={len({m['patient_id'] for m in cell_meta})}",
        flush=True,
    )


def main() -> int:
    EXT.mkdir(parents=True, exist_ok=True)
    extract_gse131907()
    extract_gse123902()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
