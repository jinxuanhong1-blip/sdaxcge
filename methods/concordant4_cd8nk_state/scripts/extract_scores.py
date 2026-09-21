#!/usr/bin/env python3
"""Patient-level CD8 / NK gene sums for the concordant-4 state test.

Exposure is not computed here. Malignant CLDN4 %pos is recomputed only as a
gate QC against the locked patient table. Scores are mean log1p(CP10k),
accumulated as sums so analyze_state.py can form modules and cell-weighted pools.
"""
from __future__ import annotations

import argparse
import gzip
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio

ROOT = Path(__file__).resolve().parents[1]

PANEL = [
    "CLDN4", "EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC",
    "CD3D", "CD3E", "CD8A", "CD8B", "CD4",
    "KLRD1", "NCR1", "KLRF1", "NCAM1", "GNLY", "NKG7",
    "GZMB", "PRF1", "IFNG", "GZMA", "GZMH",
    "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CTLA4", "ENTPD1", "CXCL13",
    "IL7R", "TCF7", "CCR7", "SELL", "LEF1",
    "CX3CR1", "FGFBP2",
]
PANEL_SET = set(PANEL)

CD8_131907 = {
    "Cytotoxic CD8+ T",
    "Exhausted CD8+ T",
    "Naive CD8+ T",
    "CD8 low T",
}

GSE189357 = {
    "TD1": "GSM5699777",
    "TD2": "GSM5699778",
    "TD3": "GSM5699779",
    "TD4": "GSM5699780",
    "TD5": "GSM5699781",
    "TD6": "GSM5699782",
    "TD7": "GSM5699783",
    "TD8": "GSM5699784",
    "TD9": "GSM5699785",
}


class Sums:
    def __init__(self) -> None:
        self.rows: list[dict] = []
        self.qc: list[dict] = []
        self.notes: list[str] = []

    def note(self, msg: str) -> None:
        print(msg, flush=True)
        self.notes.append(msg)

    def add_slice(
        self,
        dataset: str,
        unit: str,
        compartment: str,
        slice_name: str,
        logv: dict[str, np.ndarray],
        raw: dict[str, np.ndarray],
        libs: np.ndarray,
    ) -> None:
        n = int(libs.shape[0])
        if n == 0:
            return
        base = {
            "dataset": dataset,
            "unit_id": unit,
            "compartment": compartment,
            "slice": slice_name,
            "n_cells": n,
            "median_lib": float(np.median(libs)),
            "mean_lib": float(np.mean(libs)),
        }
        for gene, lv in logv.items():
            rv = raw[gene]
            pos = rv > 0
            n_pos = int(pos.sum())
            self.rows.append(
                {
                    **base,
                    "gene": gene,
                    "sum_log1p": float(lv.sum()),
                    "n_pos": n_pos,
                    "sum_log1p_pos": float(lv[pos].sum()) if n_pos else 0.0,
                }
            )

    def add_malignant(
        self,
        dataset: str,
        unit: str,
        n_cells_matrix: int,
        cldn4: np.ndarray,
    ) -> None:
        n = int(cldn4.shape[0])
        n_pos = int((cldn4 > 0).sum()) if n else 0
        self.qc.append(
            {
                "dataset": dataset,
                "unit_id": unit,
                "n_cells_matrix": int(n_cells_matrix),
                "n_malignant": n,
                "n_cldn4_pos": n_pos,
                "mal_CLDN4_pct_recomputed": (100.0 * n_pos / n) if n else np.nan,
                "mal_CLDN4_mean_log1p_raw": float(np.log1p(cldn4).mean()) if n else np.nan,
            }
        )


def _gene(raw: dict[str, np.ndarray], name: str, n: int) -> np.ndarray:
    if name not in raw:
        return np.zeros(n, dtype=np.float64)
    return raw[name]


def marker_masks(raw: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Locked malignant gate, then CD8 and NK that do not use the CosMx effector genes as the gate."""
    epi = (
        (_gene(raw, "EPCAM", n) > 0)
        | (_gene(raw, "KRT8", n) > 0)
        | (_gene(raw, "KRT18", n) > 0)
        | (_gene(raw, "KRT19", n) > 0)
    )
    mal = epi & (_gene(raw, "PTPRC", n) == 0)
    cd8a = _gene(raw, "CD8A", n)
    cd8b = _gene(raw, "CD8B", n)
    cd4 = _gene(raw, "CD4", n)
    cd3 = (_gene(raw, "CD3D", n) > 0) | (_gene(raw, "CD3E", n) > 0)
    cd8_signal = (cd8a > 0) | (cd8b > 0)
    cd4_dominant = (cd4 > cd8a) & (cd4 > cd8b)
    cd8 = (~mal) & cd8_signal & cd3 & (~cd4_dominant)
    nk_id = (
        (_gene(raw, "KLRD1", n) > 0)
        | (_gene(raw, "NCR1", n) > 0)
        | (_gene(raw, "KLRF1", n) > 0)
        | (_gene(raw, "GNLY", n) > 0)
        | ((_gene(raw, "NCAM1", n) > 0) & (_gene(raw, "PTPRC", n) > 0))
    )
    nk = (
        (~mal)
        & (~cd8)
        & (_gene(raw, "CD3D", n) == 0)
        & (_gene(raw, "CD3E", n) == 0)
        & (cd8a == 0)
        & (cd8b == 0)
        & nk_id
    )
    return mal, cd8, nk


def consume_compartment(
    sums: Sums,
    dataset: str,
    unit: str,
    compartment: str,
    raw: dict[str, np.ndarray],
    libs: np.ndarray,
    author: np.ndarray | None,
) -> None:
    ok = libs > 0
    if not np.any(ok):
        return
    libs = libs[ok]
    raw_ok = {g: v[ok] for g, v in raw.items()}
    author_ok = None if author is None else author[ok]
    logv = {g: np.log1p(raw_ok[g] / libs * 10000.0) for g in raw_ok}

    def take(mask: np.ndarray) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], np.ndarray]:
        return (
            {g: lv[mask] for g, lv in logv.items()},
            {g: rv[mask] for g, rv in raw_ok.items()},
            libs[mask],
        )

    sums.add_slice(dataset, unit, compartment, "ALL", *take(np.ones(libs.shape[0], dtype=bool)))
    umi = libs >= 500
    if np.any(umi):
        sums.add_slice(dataset, unit, compartment, "UMI500", *take(umi))
    if compartment == "cd8" and "GNLY" in raw_ok:
        gnly = raw_ok["GNLY"] > 0
        if np.any(gnly):
            sums.add_slice(dataset, unit, compartment, "GNLY_POS", *take(gnly))
    if author_ok is not None:
        for sl in pd.unique(author_ok):
            if sl in {"ALL", "UMI500", "GNLY_POS"}:
                continue
            mask = author_ok == sl
            sums.add_slice(dataset, unit, compartment, str(sl), *take(mask))


def score_marker_matrix(
    sums: Sums,
    dataset: str,
    unit: str,
    raw: dict[str, np.ndarray],
    libs: np.ndarray,
) -> None:
    n = int(libs.shape[0])
    mal, cd8, nk = marker_masks(raw, n)
    ok = libs > 0
    sums.add_malignant(dataset, unit, n, _gene(raw, "CLDN4", n)[mal & ok] if np.any(mal & ok) else np.zeros(0))
    present = sorted(set(PANEL) & set(raw))
    sums.note(
        f"  {dataset} {unit}: cells={n} mal={int(mal.sum())} cd8={int(cd8.sum())} nk={int(nk.sum())} genes={len(present)}"
    )
    for comp, mask in (("cd8", cd8), ("nk", nk)):
        if not np.any(mask):
            continue
        consume_compartment(
            sums,
            dataset,
            unit,
            comp,
            {g: raw[g][mask] for g in raw},
            libs[mask],
            None,
        )


def score_labeled_frame(sums: Sums, df: pd.DataFrame, dataset: str) -> None:
    meta = {"dataset", "unit_id", "compartment", "author_subset", "libsize"}
    genes = [c for c in df.columns if c not in meta]
    missing_score = sorted({"GZMB", "PRF1", "NKG7", "IFNG", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX"} - set(genes))
    sums.note(f"{dataset} genes in export {len(genes)}; missing score genes: {missing_score or 'none'}")
    for unit, sub_u in df.groupby("unit_id", sort=False):
        n_matrix = int(len(sub_u))
        mal = sub_u[sub_u["compartment"] == "malignant"]
        if len(mal) and "CLDN4" in mal.columns:
            cldn = mal["CLDN4"].to_numpy(dtype=np.float64)
            libs = mal["libsize"].to_numpy(dtype=np.float64)
            cldn = cldn[libs > 0]
        else:
            cldn = np.zeros(0)
        sums.add_malignant(dataset, str(unit), n_matrix, cldn)
        for comp in ("cd8", "nk"):
            ss = sub_u[sub_u["compartment"] == comp]
            if ss.empty:
                sums.note(f"  {dataset} {unit}: {comp}=0")
                continue
            libs = ss["libsize"].to_numpy(dtype=np.float64)
            raw = {g: ss[g].to_numpy(dtype=np.float64) for g in genes}
            author = ss["author_subset"].astype(str).to_numpy()
            sums.note(f"  {dataset} {unit}: {comp}={len(ss)}")
            consume_compartment(sums, dataset, str(unit), comp, raw, libs, author)


def extract_gse205335(sums: Sums, cells_path: Path) -> None:
    sums.note(f"read {cells_path}")
    df = pd.read_csv(cells_path, sep="\t")
    score_labeled_frame(sums, df, "GSE205335")


def extract_gse123902(sums: Sums, raw_dir: Path) -> None:
    units = pd.read_csv(ROOT / "data" / "GSE123902_marker_units.tsv", sep="\t")
    locked = pd.read_csv(ROOT / "data" / "locked_patient_units.tsv", sep="\t")
    keep = set(locked.loc[locked["dataset"] == "GSE123902", "unit_id"])
    tumor = units[units["tissue"].isin(["PRIMARY", "METASTASIS"]) & units["unit_id"].isin(keep)]
    if tumor["unit_id"].duplicated().any():
        raise SystemExit("GSE123902 has more than one tumor file per locked donor")
    tar_path = raw_dir / "GSE123902" / "GSE123902_RAW.tar"
    csv_dir = raw_dir / "GSE123902" / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    if not list(csv_dir.glob("*.csv.gz")):
        sums.note(f"untar {tar_path}")
        with tarfile.open(tar_path) as tf:
            tf.extractall(csv_dir)
    for rec in tumor.itertuples(index=False):
        path = csv_dir / rec.file
        sums.note(f"GSE123902 {rec.unit_id} {path.name}")
        df = pd.read_csv(path, index_col=0)
        df.columns = [str(c).upper() for c in df.columns]
        if df.columns.duplicated().any():
            df = df.T.groupby(level=0).sum().T
        libs = df.sum(axis=1).to_numpy(dtype=np.float64)
        raw = {g: df[g].to_numpy(dtype=np.float64) for g in PANEL if g in df.columns}
        del df
        score_marker_matrix(sums, "GSE123902", str(rec.unit_id), raw, libs)


def _extract_tar_member(tf: tarfile.TarFile, name: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    src = tf.extractfile(name)
    if src is None:
        raise SystemExit(f"missing {name} in tar")
    with dest.open("wb") as out:
        out.write(src.read())


def extract_gse189357(sums: Sums, raw_dir: Path) -> None:
    tar_path = raw_dir / "GSE189357" / "GSE189357_RAW.tar"
    ex = raw_dir / "GSE189357" / "raw"
    ex.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tf:
        names = set(tf.getnames())
        for unit, gsm in GSE189357.items():
            prefix = f"{gsm}_{unit}"
            for suffix in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
                member = f"{prefix}_{suffix}"
                if member not in names:
                    alt = [n for n in names if n.endswith(f"{unit}_{suffix}")]
                    if not alt:
                        raise SystemExit(f"missing {member}")
                    member = alt[0]
                _extract_tar_member(tf, member, ex / f"{prefix}_{suffix}")
            sums.note(f"GSE189357 {unit}")
            feat = pd.read_csv(ex / f"{prefix}_features.tsv.gz", sep="\t", header=None)
            symbols = feat.iloc[:, 1].astype(str).str.upper()
            X = spio.mmread(ex / f"{prefix}_matrix.mtx.gz")
            if hasattr(X, "tocsr"):
                X = X.tocsr()
            else:
                X = spio.mmread(ex / f"{prefix}_matrix.mtx.gz").tocsr()
            if X.shape[0] != len(symbols):
                X = X.T.tocsr()
            if X.shape[0] != len(symbols):
                raise SystemExit(f"{unit} matrix {X.shape} vs features {len(symbols)}")
            libs = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
            wanted = defaultdict(list)
            for i, sym in enumerate(symbols):
                if sym in PANEL_SET:
                    wanted[sym].append(i)
            raw = {g: np.asarray(X[rows].sum(axis=0)).ravel().astype(np.float64) for g, rows in wanted.items()}
            del X
            score_marker_matrix(sums, "GSE189357", unit, raw, libs)


def extract_gse131907(sums: Sums, raw_dir: Path) -> None:
    locked = pd.read_csv(ROOT / "data" / "locked_patient_units.tsv", sep="\t")
    keep_samples = set(locked.loc[locked["dataset"] == "GSE131907", "unit_id"])
    ann = raw_dir / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat = raw_dir / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    by_index: dict[str, tuple[str, str, str]] = {}
    with gzip.open(ann, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            sample = parts[idx["Sample"]]
            if sample not in keep_samples:
                continue
            subtype = parts[idx["Cell_subtype"]]
            cell_type = parts[idx["Cell_type"]]
            if subtype == "Malignant cells":
                comp, sl = "malignant", "Malignant"
            elif subtype in CD8_131907:
                comp, sl = "cd8", subtype
            elif cell_type == "NK cells" and subtype not in CD8_131907:
                comp, sl = "nk", subtype if subtype not in {"", "NA"} else "NK"
            else:
                continue
            by_index[parts[idx["Index"]]] = (sample, comp, sl)
    sums.note(f"GSE131907 annotated keep {len(by_index)}")

    with gzip.open(mat, "rt") as handle:
        barcodes = handle.readline().rstrip("\n").split("\t")[1:]
        keep_cols: list[int] = []
        meta: list[tuple[str, str, str]] = []
        for i, bc in enumerate(barcodes):
            rec = by_index.get(bc)
            if rec is not None:
                keep_cols.append(i)
                meta.append(rec)
        n_keep = len(keep_cols)
        sums.note(f"GSE131907 matrix keep {n_keep} / {len(barcodes)}")
        if n_keep < 1000:
            raise SystemExit("GSE131907 keep set too small; barcode key mismatch")
        keep_cols_arr = np.asarray(keep_cols, dtype=np.int32)
        lib = np.zeros(n_keep, dtype=np.float64)
        panel: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in handle:
            n_genes += 1
            parts = line.rstrip("\n").split("\t")
            gene = parts[0].upper()
            vals = np.fromiter(
                (float(parts[c + 1]) if parts[c + 1] else 0.0 for c in keep_cols_arr),
                dtype=np.float64,
                count=n_keep,
            )
            lib += vals
            if gene in PANEL_SET:
                if gene in panel:
                    panel[gene] += vals
                else:
                    panel[gene] = vals
            if n_genes % 4000 == 0:
                sums.note(f"  GSE131907 scanned {n_genes} genes ({gene})")
        sums.note(f"GSE131907 scanned {n_genes} genes; panel hit {sorted(panel)}")

    meta_df = pd.DataFrame(meta, columns=["unit_id", "compartment", "author_subset"])
    # group
    for (unit, comp), sub in meta_df.groupby(["unit_id", "compartment"], sort=False):
        ix = sub.index.to_numpy()
        if comp == "malignant":
            cldn = panel["CLDN4"][ix] if "CLDN4" in panel else np.zeros(len(ix))
            libs = lib[ix]
            cldn = cldn[libs > 0]
            # n_cells_matrix is not the full sample; QC uses n_malignant only
            sums.add_malignant("GSE131907", str(unit), int(len(ix)), cldn)
            continue
        raw = {g: arr[ix] for g, arr in panel.items()}
        author = sub["author_subset"].astype(str).to_numpy()
        sums.note(f"  GSE131907 {unit}: {comp}={len(ix)}")
        consume_compartment(sums, "GSE131907", str(unit), str(comp), raw, lib[ix], author)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("/tmp/concordant4_raw"))
    ap.add_argument("--cells205", type=Path, default=Path("/tmp/concordant4_state/gse205335_cells.tsv.gz"))
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "tables")
    ap.add_argument("--dataset", default="all", choices=["all", "GSE123902", "GSE131907", "GSE189357", "GSE205335"])
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    sums = Sums()
    which = args.dataset
    if which in {"all", "GSE123902"}:
        extract_gse123902(sums, args.raw)
    if which in {"all", "GSE189357"}:
        extract_gse189357(sums, args.raw)
    if which in {"all", "GSE131907"}:
        extract_gse131907(sums, args.raw)
    if which in {"all", "GSE205335"}:
        if not args.cells205.exists():
            raise SystemExit(f"missing {args.cells205}; run export_gse205335.R first")
        extract_gse205335(sums, args.cells205)

    tag = "" if which == "all" else f".{which}"
    sums_path = args.outdir / f"compartment_gene_sums{tag}.tsv"
    qc_path = args.outdir / f"malignant_qc{tag}.tsv"
    log_path = args.outdir / f"extract_log{tag}.txt"
    pd.DataFrame(sums.rows).to_csv(sums_path, sep="\t", index=False)
    pd.DataFrame(sums.qc).to_csv(qc_path, sep="\t", index=False)
    log_path.write_text("\n".join(sums.notes) + "\n")
    sums.note(f"wrote {sums_path} rows={len(sums.rows)} qc={len(sums.qc)}")


if __name__ == "__main__":
    main()
