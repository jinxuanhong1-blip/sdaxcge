#!/usr/bin/env python3
"""Full-unit T / NK / myeloid counts for the locked concordant-4 (n=65).

Author lineages: GSE131907 (sample) and GSE205335 (patient, Normal* dropped).
Marker lineages: GSE123902 (donor tumor) and GSE189357 (patient).
The malignant CLDN4 score and quartile are copied from the locked patient table.
They are not recomputed.
"""
from __future__ import annotations

import gzip
import io as pyio
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
GEO = Path("/tmp/geo_c4")
OUT = ROOT / "results" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

MARKER_MAL = ("EPCAM", "KRT8", "KRT18", "KRT19")
MARKER_T = ("CD3D", "CD3E", "CD3G", "CD4", "CD8A")
MARKER_NK = ("NKG7", "GNLY", "KLRD1", "NCAM1", "NCR1")
MARKER_MYE = ("LYZ", "CD68", "CD14", "FCGR3A", "CD163", "C1QA", "MARCO", "AIF1")
MARKER_B = ("MS4A1", "CD79A", "CD19")


def _pos(mat_row: np.ndarray) -> np.ndarray:
    return np.asarray(mat_row).ravel() > 0


def classify_marker(gene_index: dict[str, int], get_row) -> dict[str, np.ndarray]:
    """Mutually exclusive classes. Locked malignant gate is unchanged.

    T is claimed before NK so cytotoxic T cells that co-express NKG7 stay T.
    Locked T/NK genes are a subset of T ∪ NK, so n_T + n_NK recovers that gate
    when CD8A/CD3/CD4 cover the T side and NKG7/GNLY/KLRD1 cover the NK side.
    """

    def any_of(genes: tuple[str, ...]) -> np.ndarray:
        acc = None
        for g in genes:
            if g not in gene_index:
                continue
            hit = _pos(get_row(gene_index[g]))
            acc = hit if acc is None else (acc | hit)
        if acc is None:
            raise RuntimeError(f"none of {genes} found")
        return acc

    n = int(len(get_row(next(iter(gene_index.values())))))
    if "PTPRC" in gene_index:
        ptprc = _pos(get_row(gene_index["PTPRC"]))
    else:
        ptprc = np.zeros(n, dtype=bool)
    mal = any_of(MARKER_MAL) & ~ptprc
    t = ~mal & any_of(MARKER_T)
    nk = ~mal & ~t & any_of(MARKER_NK)
    mye = ~mal & ~t & ~nk & any_of(MARKER_MYE)
    b = ~mal & ~t & ~nk & ~mye & any_of(MARKER_B)
    locked_tnk_genes = ("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1")
    locked = ~mal & any_of(locked_tnk_genes)
    return {
        "malignant": mal,
        "T": t,
        "NK": nk,
        "myeloid": mye,
        "B": b,
        "locked_tnk": locked,
        "n": n,
    }


def _row_from_series(df: pd.DataFrame, gene: str) -> np.ndarray:
    # df columns are cells, index is uppercase genes. Duplicate symbols summed.
    if gene not in df.index:
        return np.zeros(df.shape[1], dtype=np.float64)
    sub = df.loc[gene]
    if isinstance(sub, pd.DataFrame):
        return sub.to_numpy(dtype=np.float64).sum(axis=0)
    return sub.to_numpy(dtype=np.float64)


def counts_from_flags(flags: dict[str, np.ndarray]) -> dict[str, int]:
    n = int(flags["n"] if not isinstance(flags["n"], np.ndarray) else flags["n"])
    n_t = int(flags["T"].sum())
    n_nk = int(flags["NK"].sum())
    n_mye = int(flags["myeloid"].sum())
    n_b = int(flags["B"].sum())
    n_mal = int(flags["malignant"].sum())
    n_locked = int(flags["locked_tnk"].sum())
    return {
        "n_cells": n,
        "n_T": n_t,
        "n_NK": n_nk,
        "n_myeloid": n_mye,
        "n_B": n_b,
        "n_malignant": n_mal,
        "n_tnk_locked_gate": n_locked,
        "n_tnk_split": n_t + n_nk,
        "n_rest": n - n_t - n_nk - n_mye,
        "label_source": "marker",
    }


def load_123902(files: list[str]) -> dict[str, dict[str, int]]:
    tar_path = GEO / "GSE123902_RAW.tar"
    want = set(files)
    out: dict[str, dict[str, int]] = {}
    with tarfile.open(tar_path) as tar:
        for member in tar.getmembers():
            if member.name not in want:
                continue
            raw = gzip.GzipFile(fileobj=tar.extractfile(member)).read()
            df = pd.read_csv(pyio.BytesIO(raw), index_col=0)
            # cells x genes -> genes x cells
            df = df.T
            df.index = df.index.astype(str).str.upper()
            df = df.groupby(level=0).sum()
            gene_index = {g: i for i, g in enumerate(df.index)}
            cache: dict[int, np.ndarray] = {}

            def get_row(i: int, _df=df, _cache=cache) -> np.ndarray:
                if i not in _cache:
                    _cache[i] = _df.iloc[i].to_numpy(dtype=np.float64)
                return _cache[i]

            flags = classify_marker(gene_index, get_row)
            flags["n"] = df.shape[1]
            out[member.name] = counts_from_flags(flags)
            print(f"  GSE123902 {member.name} cells={df.shape[1]}", flush=True)
    missing = want - set(out)
    if missing:
        raise SystemExit(f"GSE123902 missing from tar: {sorted(missing)}")
    return out


def _read_10x_from_tar(tar: tarfile.TarFile, prefix: str) -> sparse.csc_matrix:
    def member(suffix: str):
        name = f"{prefix}{suffix}"
        m = tar.extractfile(name)
        if m is None:
            raise SystemExit(f"missing {name}")
        return gzip.GzipFile(fileobj=m)

    mtx = spio.mmread(member("_matrix.mtx.gz")).tocsc()
    features = pd.read_csv(member("_features.tsv.gz"), sep="\t", header=None)
    # 10x v3: id, symbol, type. v2 may be two columns.
    symbol_col = 1 if features.shape[1] > 1 else 0
    symbols = features.iloc[:, symbol_col].astype(str).str.upper()
    return mtx, symbols


def load_189357(patients: list[str]) -> dict[str, dict[str, int]]:
    tar_path = GEO / "GSE189357_RAW.tar"
    out: dict[str, dict[str, int]] = {}
    with tarfile.open(tar_path) as tar:
        names = tar.getnames()
        for pat in patients:
            hits = [n for n in names if n.endswith(f"_{pat}_matrix.mtx.gz")]
            if len(hits) != 1:
                raise SystemExit(f"GSE189357 matrix count for {pat}: {hits}")
            prefix = hits[0][: -len("_matrix.mtx.gz")]
            mtx, symbols = _read_10x_from_tar(tar, prefix)
            # sum duplicate symbols into a compact gene x cell matrix of needed genes
            need = set(MARKER_MAL + MARKER_T + MARKER_NK + MARKER_MYE + MARKER_B + ("PTPRC",))
            groups: dict[str, list[int]] = defaultdict(list)
            for i, g in enumerate(symbols.tolist()):
                if g in need:
                    groups[g].append(i)
            rows = []
            names_g = []
            for g, idxs in groups.items():
                if len(idxs) == 1:
                    rows.append(mtx[idxs[0]])
                else:
                    rows.append(sparse.csr_matrix(mtx[idxs].sum(axis=0)))
                names_g.append(g)
            if not rows:
                raise SystemExit(f"no marker genes in {pat}")
            sub = sparse.vstack(rows).tocsr()
            gene_index = {g: i for i, g in enumerate(names_g)}

            def get_row(i: int, _sub=sub) -> np.ndarray:
                return np.asarray(_sub[i].todense()).ravel()

            flags = classify_marker(gene_index, get_row)
            flags["n"] = mtx.shape[1]
            out[pat] = counts_from_flags(flags)
            print(f"  GSE189357 {pat} cells={mtx.shape[1]}", flush=True)
    return out


def load_131907(samples: list[str]) -> dict[str, dict[str, int]]:
    path = GEO / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    acc: dict[str, dict[str, int]] = {
        s: defaultdict(int) for s in samples
    }
    want = set(samples)
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            sample = p[idx["Sample"]]
            if sample not in want:
                continue
            ctype = p[idx["Cell_type"]]
            sub = p[idx["Cell_subtype"]]
            rec = acc[sample]
            rec["n_cells"] += 1
            if sub == "Malignant cells":
                rec["n_malignant"] += 1
            if ctype == "T lymphocytes":
                rec["n_T"] += 1
            elif ctype == "NK cells":
                rec["n_NK"] += 1
            elif ctype == "Myeloid cells":
                rec["n_myeloid"] += 1
            elif ctype == "B lymphocytes":
                rec["n_B"] += 1
            elif ctype == "MAST cells":
                rec["n_mast"] += 1
    out = {}
    for s, rec in acc.items():
        n_t = int(rec["n_T"])
        n_nk = int(rec["n_NK"])
        n_mye = int(rec["n_myeloid"])
        n = int(rec["n_cells"])
        out[s] = {
            "n_cells": n,
            "n_T": n_t,
            "n_NK": n_nk,
            "n_myeloid": n_mye,
            "n_B": int(rec["n_B"]),
            "n_mast": int(rec["n_mast"]),
            "n_malignant": int(rec["n_malignant"]),
            "n_tnk_locked_gate": n_t + n_nk,
            "n_tnk_split": n_t + n_nk,
            "n_rest": n - n_t - n_nk - n_mye,
            "label_source": "author_Cell_type",
        }
    return out


def load_205335(patients: list[str]) -> dict[str, dict[str, int]]:
    gsm = pd.read_csv(DATA / "GSE205335_gsm_map.tsv", sep="\t")
    normal = set(gsm.loc[gsm["tissue"].str.startswith("Normal"), "orig.ident"])
    patient_of = dict(zip(gsm["orig.ident"], gsm["patient"]))
    want = set(patients)
    acc: dict[str, dict[str, int]] = {p: defaultdict(int) for p in patients}
    path = GEO / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            orig = p[idx["orig.ident"]]
            if orig in normal or orig not in patient_of:
                continue
            pat = patient_of[orig]
            if pat not in want:
                continue
            total = p[idx["lineage.total"]]
            sub = p[idx["lineage.sub"]]
            rec = acc[pat]
            rec["n_cells"] += 1
            if sub == "Malignant cells":
                rec["n_malignant"] += 1
            if sub == "NK cells":
                rec["n_NK"] += 1
            elif total == "T/NK cells":
                rec["n_T"] += 1
            elif total == "Myeloid cells":
                rec["n_myeloid"] += 1
            elif total == "B/Plasma cells":
                rec["n_B"] += 1
            elif total == "Mast cells":
                rec["n_mast"] += 1
    out = {}
    for pat, rec in acc.items():
        n_t = int(rec["n_T"])
        n_nk = int(rec["n_NK"])
        n_mye = int(rec["n_myeloid"])
        n = int(rec["n_cells"])
        out[pat] = {
            "n_cells": n,
            "n_T": n_t,
            "n_NK": n_nk,
            "n_myeloid": n_mye,
            "n_B": int(rec["n_B"]),
            "n_mast": int(rec["n_mast"]),
            "n_malignant": int(rec["n_malignant"]),
            "n_tnk_locked_gate": n_t + n_nk,
            "n_tnk_split": n_t + n_nk,
            "n_rest": n - n_t - n_nk - n_mye,
            "label_source": "author_lineage",
        }
    return out


def main() -> None:
    locked = pd.read_csv(DATA / "locked_patient_units.tsv", sep="\t")
    if len(locked) != 65:
        raise SystemExit(f"expected 65 locked units, got {len(locked)}")

    print("author GSE131907", flush=True)
    c131 = load_131907(locked.loc[locked.dataset == "GSE131907", "unit_id"].tolist())
    print("author GSE205335", flush=True)
    c205 = load_205335(locked.loc[locked.dataset == "GSE205335", "unit_id"].tolist())

    marker = pd.read_csv(GEO / "GSE123902_marker_units.tsv", sep="\t")
    marker = marker[marker.tissue.isin(["PRIMARY", "METASTASIS"])]
    marker = marker.sort_values(["patient", "tissue"]).drop_duplicates("patient")
    file_of = dict(zip(marker.patient, marker.file))
    donors = locked.loc[locked.dataset == "GSE123902", "unit_id"].tolist()
    missing_files = [d for d in donors if d not in file_of]
    if missing_files:
        raise SystemExit(f"no tumor file for {missing_files}")
    print("marker GSE123902", flush=True)
    by_file = load_123902([file_of[d] for d in donors])
    c123 = {d: by_file[file_of[d]] for d in donors}

    print("marker GSE189357", flush=True)
    pats = locked.loc[locked.dataset == "GSE189357", "unit_id"].tolist()
    c189 = load_189357(pats)

    bundles = {
        "GSE123902": c123,
        "GSE131907": c131,
        "GSE205335": c205,
        "GSE189357": c189,
    }
    rows = []
    for rec in locked.itertuples(index=False):
        got = bundles[rec.dataset][rec.unit_id]
        rows.append(
            {
                "dataset": rec.dataset,
                "unit_id": rec.unit_id,
                "unit_type": rec.unit_type,
                "tissue": rec.tissue,
                "mal_CLDN4_pct": rec.mal_CLDN4_pct,
                "cldn4_quartile": rec.cldn4_quartile,
                "locked_n_cells": rec.n_cells,
                "locked_n_malignant": rec.n_malignant,
                "locked_n_tnk": rec.n_tnk,
                "locked_frac_tnk": rec.frac_tnk,
                **got,
            }
        )
    tab = pd.DataFrame(rows)
    tab["frac_T"] = tab["n_T"] / tab["n_cells"]
    tab["frac_NK"] = tab["n_NK"] / tab["n_cells"]
    tab["frac_myeloid"] = tab["n_myeloid"] / tab["n_cells"]
    tab["frac_rest"] = tab["n_rest"] / tab["n_cells"]
    tab["frac_tnk_split"] = (tab["n_T"] + tab["n_NK"]) / tab["n_cells"]
    tab["focus_n"] = tab["n_T"] + tab["n_NK"] + tab["n_myeloid"]

    # Author cohorts must reproduce the locked T/NK count exactly.
    author = tab[tab.dataset.isin(["GSE131907", "GSE205335"])]
    bad = author[author.n_tnk_locked_gate != author.locked_n_tnk]
    if len(bad):
        print(bad[["dataset", "unit_id", "n_tnk_locked_gate", "locked_n_tnk", "n_cells", "locked_n_cells"]])
        raise SystemExit("author T/NK counts do not match the locked table")
    bad_n = author[author.n_cells != author.locked_n_cells]
    if len(bad_n):
        print(bad_n[["dataset", "unit_id", "n_cells", "locked_n_cells"]])
        raise SystemExit("author cell totals do not match the locked table")

    path = OUT / "composition_counts.tsv"
    tab.to_csv(path, sep="\t", index=False)
    print(f"wrote {path} n={len(tab)}", flush=True)
    for ds, sub in tab.groupby("dataset"):
        d_tnk = int((sub.n_tnk_split - sub.locked_n_tnk).abs().sum())
        print(
            f"  {ds} units={len(sub)} L1|split-locked T/NK|={d_tnk} "
            f"median frac T/NK/mye="
            f"{sub.frac_T.median():.3f}/{sub.frac_NK.median():.3f}/{sub.frac_myeloid.median():.3f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
