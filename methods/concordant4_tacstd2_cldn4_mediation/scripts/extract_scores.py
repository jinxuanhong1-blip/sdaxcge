#!/usr/bin/env python3
"""Malignant TACSTD2 and CLDN4 scores on the locked concordant-4 units.

Cohorts, and only these: GSE123902, GSE131907, GSE205335, GSE189357.
Gates match the locked patient / donor / sample tables (n = 65).
GEO matrices are read from /tmp/geo_dl and are not written into the repo.
"""
from __future__ import annotations

import gzip
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
GEO = Path("/tmp/geo_dl")
OUT = HERE / "results"

MARKER_EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
MARKER_TNK = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
FOCUS = ["TACSTD2", "CLDN4"]
KERATIN = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT6A"]
PANEL = list(dict.fromkeys(FOCUS + KERATIN + ["PTPRC"] + MARKER_TNK))
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}

LOCKED_IDS = {
    "GSE123902": [
        "LX255B", "LX653", "LX661", "LX666", "LX675", "LX676", "LX679",
        "LX680", "LX681", "LX682", "LX684", "LX699", "LX701",
    ],
    "GSE131907": [
        "BRONCHO_11", "BRONCHO_58", "EBUS_06", "EBUS_10", "EBUS_12", "EBUS_13",
        "EBUS_15", "EBUS_19", "EBUS_28", "EBUS_49", "EBUS_51", "NS_02", "NS_03",
        "NS_04", "NS_06", "NS_07", "NS_12", "NS_13", "NS_16", "NS_17", "NS_19",
    ],
    "GSE205335": [
        "P0031", "P1006", "P1015", "P1016", "P1017", "P1018", "P1025", "P1027",
        "P1030", "P1037", "P1056", "P1062", "P1063", "P1072", "P1076", "P1079",
        "P1084", "P1089", "P1090", "P1115", "P1119", "P4001",
    ],
    "GSE189357": [f"TD{i}" for i in range(1, 10)],
}


def _upper_gene(name: str) -> str:
    return str(name).upper().split(".")[0].strip()


def _as_int_umi(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=np.float64)
    rounded = np.rint(arr)
    if arr.size and np.nanmax(np.abs(arr - rounded)) > 1e-3:
        raise ValueError("expression is not integer UMI")
    return rounded


def _hist(arr: np.ndarray) -> str:
    if arr.size == 0:
        return ""
    vals, cts = np.unique(arr.astype(np.int64), return_counts=True)
    return ",".join(f"{int(v)}:{int(c)}" for v, c in zip(vals, cts))


def _pct_mask(mask: np.ndarray) -> float:
    if mask.size == 0:
        return float("nan")
    return float(np.mean(mask))


def score_gene(x: np.ndarray, lib: np.ndarray) -> dict:
    x = _as_int_umi(x)
    lib = np.maximum(lib.astype(np.float64), 1.0)
    cp = x / lib * 1e4
    pos = x > 0
    pb = float(x.sum() / lib.sum() * 1e4) if lib.size else float("nan")
    out = {
        "pct_gt0": _pct_mask(x > 0),
        "pct_ge2": _pct_mask(x >= 2),
        "pct_ge3": _pct_mask(x >= 3),
        "pct_ge5": _pct_mask(x >= 5),
        "pct_ge10": _pct_mask(x >= 10),
        "pct_cp_ge1": _pct_mask(cp >= 1),
        "pct_cp_ge5": _pct_mask(cp >= 5),
        "pct_cp_ge10": _pct_mask(cp >= 10),
        "mean_log1p_raw": float(np.log1p(x).mean()) if x.size else float("nan"),
        "mean_log1p_cp10k": float(np.log1p(cp).mean()) if x.size else float("nan"),
        "mean_log1p_pos_raw": float(np.log1p(x[pos]).mean()) if np.any(pos) else float("nan"),
        "p90_log1p_raw": float(np.quantile(np.log1p(x), 0.9)) if x.size else float("nan"),
        "p90_log1p_cp10k": float(np.quantile(np.log1p(cp), 0.9)) if x.size else float("nan"),
        "pb_cp10k": pb,
        "pb_log1p_cp10k": float(np.log1p(pb)) if np.isfinite(pb) else float("nan"),
        "hist": _hist(x),
    }
    return out


def score_unit(genes: dict[str, np.ndarray], lib: np.ndarray, mask: np.ndarray) -> dict:
    L = lib[mask]
    row: dict = {"n_malignant": int(mask.sum())}
    xt = _as_int_umi(genes["TACSTD2"][mask])
    xc = _as_int_umi(genes["CLDN4"][mask])
    n_mal = int(xt.size)
    if n_mal == 0:
        both = t_only = c_only = float("nan")
    else:
        both = float(np.mean((xt > 0) & (xc > 0)))
        t_only = float(np.mean((xt > 0) & (xc == 0)))
        c_only = float(np.mean((xc > 0) & (xt == 0)))
    row["pct_both_gt0"] = both
    row["pct_tacstd2_only_gt0"] = t_only
    row["pct_cldn4_only_gt0"] = c_only
    for gene in FOCUS:
        if gene not in genes:
            raise SystemExit(f"missing {gene}")
        scored = score_gene(genes[gene][mask], L)
        for key, val in scored.items():
            row[f"{gene.lower()}_{key}"] = val
    for gene in KERATIN:
        key = gene.lower()
        if gene not in genes:
            row[f"{key}_mean_log1p_raw"] = float("nan")
            row[f"{key}_mean_log1p_cp10k"] = float("nan")
            continue
        g = _as_int_umi(genes[gene][mask])
        lib_m = np.maximum(L.astype(np.float64), 1.0)
        row[f"{key}_mean_log1p_raw"] = float(np.log1p(g).mean()) if g.size else float("nan")
        row[f"{key}_mean_log1p_cp10k"] = (
            float(np.log1p(g / lib_m * 1e4).mean()) if g.size else float("nan")
        )
    return row


def marker_masks(genes: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    epi = np.zeros(n, dtype=bool)
    for g in MARKER_EPI:
        if g in genes:
            epi |= genes[g] > 0
    ptprc = genes["PTPRC"] if "PTPRC" in genes else np.zeros(n)
    mal = epi & (ptprc == 0)
    tnk = np.zeros(n, dtype=bool)
    for g in MARKER_TNK:
        if g in genes:
            tnk |= genes[g] > 0
    tnk = tnk & ~mal
    return mal, tnk


def _collapse_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [_upper_gene(c) for c in df.columns]
    if df.columns.duplicated().any():
        df = df.T.groupby(level=0).sum().T
    return df


def _finish(rec: dict, cohort: str, unit_id: str, n_cells: int, n_tnk: int, source: str) -> dict:
    rec.update(
        {
            "cohort": cohort,
            "unit_id": unit_id,
            "n_cells": int(n_cells),
            "n_tnk": int(n_tnk),
            "frac_tnk": float(n_tnk / n_cells) if n_cells else float("nan"),
            "source": source,
        }
    )
    return rec


def gse123902() -> list[dict]:
    tar_path = GEO / "GSE123902" / "GSE123902_RAW.tar"
    wanted = set(LOCKED_IDS["GSE123902"])
    rows = []
    with tarfile.open(tar_path) as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".csv.gz")]
        for m in members:
            name = Path(m.name).name
            if "NORMAL" in name:
                continue
            if "PRIMARY" not in name and "METASTASIS" not in name:
                continue
            donor = name.split("_")[2]
            if donor not in wanted:
                continue
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            df = pd.read_csv(raw, index_col=0)
            lib = df.sum(axis=1).to_numpy(dtype=np.float64)
            df = _collapse_columns(df)
            genes = {g: df[g].to_numpy(dtype=np.float64) for g in PANEL if g in df.columns}
            missing = [g for g in FOCUS + ["PTPRC"] + MARKER_EPI + MARKER_TNK if g not in genes]
            if missing:
                raise SystemExit(f"{name} missing {missing}")
            mal, tnk = marker_masks(genes, len(df))
            rec = score_unit(genes, lib, mal)
            rows.append(_finish(rec, "GSE123902", donor, len(df), int(tnk.sum()), name))
            print(
                f"  {donor}: cells={len(df)} mal={rec['n_malignant']} tnk={int(tnk.sum())} "
                f"cldn4%={rec['cldn4_pct_gt0']:.4f} tacstd2%={rec['tacstd2_pct_gt0']:.4f}",
                flush=True,
            )
    return rows


def _tar_members(tf: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    return {Path(m.name).name: m for m in tf.getmembers()}


def _read_gz_lines(tf: tarfile.TarFile, member: tarfile.TarInfo):
    with gzip.GzipFile(fileobj=tf.extractfile(member)) as handle:
        for line in handle:
            yield line.decode()


def gse189357() -> list[dict]:
    tar_path = GEO / "GSE189357" / "GSE189357_RAW.tar"
    rows = []
    with tarfile.open(tar_path) as tf:
        members = _tar_members(tf)
        for i in range(1, 10):
            sample = f"TD{i}"
            feat_name = next(n for n in members if n.endswith(f"_{sample}_features.tsv.gz"))
            bar_name = next(n for n in members if n.endswith(f"_{sample}_barcodes.tsv.gz"))
            mtx_name = next(n for n in members if n.endswith(f"_{sample}_matrix.mtx.gz"))
            genes = []
            for line in _read_gz_lines(tf, members[feat_name]):
                parts = line.rstrip("\n").split("\t")
                sym = parts[1] if len(parts) > 1 else parts[0]
                genes.append(_upper_gene(sym))
            n_cells = sum(1 for _ in _read_gz_lines(tf, members[bar_name]))
            gene_to_cols: dict[str, list[int]] = {}
            for gi, g in enumerate(genes):
                if g in PANEL:
                    gene_to_cols.setdefault(g, []).append(gi)
            counts = {g: np.zeros(n_cells, dtype=np.float64) for g in gene_to_cols}
            lib = np.zeros(n_cells, dtype=np.float64)
            row_to_gene = {gi: g for g, gis in gene_to_cols.items() for gi in gis}
            with gzip.GzipFile(fileobj=tf.extractfile(members[mtx_name])) as handle:
                dims = None
                for line in handle:
                    if line.startswith(b"%"):
                        continue
                    dims = line.decode().split()
                    break
                if dims is None:
                    raise SystemExit(f"{sample} mtx missing dims")
                n_rows, n_cols = int(dims[0]), int(dims[1])
                if n_rows != len(genes) or n_cols != n_cells:
                    raise SystemExit(
                        f"{sample} mtx shape {n_rows}x{n_cols} != genes {len(genes)} x cells {n_cells}"
                    )
                for line in handle:
                    a, b, v = line.split()
                    gi = int(a) - 1
                    ci = int(b) - 1
                    val = float(v)
                    lib[ci] += val
                    g = row_to_gene.get(gi)
                    if g is not None:
                        counts[g][ci] += val
            missing = [g for g in FOCUS + ["PTPRC"] + MARKER_EPI + MARKER_TNK if g not in counts]
            if missing:
                raise SystemExit(f"{sample} missing {missing}")
            mal, tnk = marker_masks(counts, n_cells)
            rec = score_unit(counts, lib, mal)
            rows.append(_finish(rec, "GSE189357", sample, n_cells, int(tnk.sum()), mtx_name))
            print(
                f"  {sample}: cells={n_cells} mal={rec['n_malignant']} tnk={int(tnk.sum())} "
                f"cldn4%={rec['cldn4_pct_gt0']:.4f} tacstd2%={rec['tacstd2_pct_gt0']:.4f}",
                flush=True,
            )
    return rows


def gse131907() -> list[dict]:
    ann = pd.read_csv(GEO / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    matrix = GEO / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    print(f"streaming {matrix.name}", flush=True)
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        n_cells = len(cells)
        lib = np.zeros(n_cells, dtype=np.float64)
        n_rows = 0
        panel = set(PANEL)
        for line in handle:
            tab = line.find("\t")
            gene = _upper_gene(line[:tab])
            arr = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise SystemExit(f"{gene}: {arr.size} values, expected {n_cells}")
            lib += arr
            n_rows += 1
            if gene in panel:
                if gene in found:
                    found[gene] = found[gene] + arr
                else:
                    found[gene] = arr.astype(np.float64, copy=False)
            if n_rows % 4000 == 0:
                print(f"  streamed {n_rows} genes", flush=True)
    print(f"scanned {n_rows} genes; panel found {sorted(found)}", flush=True)
    missing = [g for g in FOCUS if g not in found]
    if missing:
        raise SystemExit(f"GSE131907 missing {missing}")
    ann = ann.set_index("Index").reindex(cells)
    if ann["Sample"].isna().any():
        raise SystemExit("GSE131907 barcodes do not match annotation Index")
    mal = ann["Cell_subtype"].eq("Malignant cells").to_numpy()
    tnk = ann["Cell_type"].isin(["T lymphocytes", "NK cells"]).to_numpy()
    sample = ann["Sample"].astype(str).to_numpy()
    origin = ann["Sample_Origin"].astype(str).to_numpy()
    wanted = set(LOCKED_IDS["GSE131907"])
    rows = []
    for sid in wanted:
        sm = sample == sid
        if not np.any(sm):
            raise SystemExit(f"GSE131907 missing sample {sid}")
        ori = str(origin[sm][0])
        if ori not in TUMOR_ORIGINS:
            raise SystemExit(f"{sid} origin {ori} is not a locked tumor-bearing site")
        mal_m = sm & mal
        rec = score_unit(found, lib, mal_m)
        n = int(sm.sum())
        n_tnk = int((sm & tnk).sum())
        rows.append(_finish(rec, "GSE131907", sid, n, n_tnk, matrix.name))
        print(
            f"  {sid} {ori}: cells={n} mal={rec['n_malignant']} tnk={n_tnk} "
            f"cldn4%={rec['cldn4_pct_gt0']:.4f} tacstd2%={rec['tacstd2_pct_gt0']:.4f}",
            flush=True,
        )
    return rows


def gse205335() -> list[dict]:
    import rdata
    from scipy import sparse

    meta = pd.read_csv(DATA / "GSE205335_gsm_metadata.csv")
    tumor = meta[~meta["tissue"].astype(str).str.startswith("Normal")].copy()
    keep_ident = set(tumor["orig.ident"].astype(str))
    ident = pd.read_csv(GEO / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    ident["barcode"] = ident["barcode"].astype(str)
    ident["orig.ident"] = ident["orig.ident"].astype(str)
    ident = ident[ident["orig.ident"].isin(keep_ident)].copy()
    ident = ident.merge(tumor[["orig.ident", "patient"]], on="orig.ident", how="left")
    if ident["patient"].isna().any():
        raise SystemExit("GSE205335 identity rows missing patient")
    wanted = set(LOCKED_IDS["GSE205335"])
    ident = ident[ident["patient"].astype(str).isin(wanted)].copy()

    rds_gz = GEO / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    rds_path = GEO / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds"
    if not rds_path.exists() or rds_path.stat().st_size < 1000:
        print("decompress GSE205335 RDS", flush=True)
        with gzip.open(rds_gz, "rb") as src, rds_path.open("wb") as dst:
            while True:
                chunk = src.read(16 * 1024 * 1024)
                if not chunk:
                    break
                dst.write(chunk)
    print("read GSE205335 RDS", flush=True)
    obj = rdata.read_rds(str(rds_path))
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise SystemExit(f"RDS missing dgCMatrix fields: {sorted(vars(obj))[:30]}")
    genes = np.array([_upper_gene(g) for g in np.asarray(obj.Dimnames[0], dtype=str)])
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    lib = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    wanted_genes = {}
    for gene in PANEL:
        hits = np.flatnonzero(genes == gene)
        if hits.size == 0:
            continue
        acc = np.zeros(matrix.shape[1], dtype=np.float64)
        for h in hits:
            acc += np.asarray(matrix.getrow(int(h)).todense()).ravel()
        wanted_genes[gene] = acc
    del obj, matrix
    missing = [g for g in FOCUS if g not in wanted_genes]
    if missing:
        raise SystemExit(f"GSE205335 RDS missing {missing}")
    bc_index = pd.Index(barcodes.astype(str))
    missing_bc = int((~ident["barcode"].isin(bc_index)).sum())
    if missing_bc:
        raise SystemExit(f"GSE205335 identity barcodes missing from RDS: {missing_bc}")
    ident = ident.copy()
    ident["pos"] = bc_index.get_indexer(ident["barcode"])
    rows = []
    for patient, sub in ident.groupby("patient"):
        patient = str(patient)
        if patient not in wanted:
            continue
        idx = sub["pos"].to_numpy()
        mal_local = sub["lineage.sub"].eq("Malignant cells").to_numpy()
        tnk_local = sub["lineage.total"].eq("T/NK cells").to_numpy()
        mask = np.zeros(len(barcodes), dtype=bool)
        mask[idx[mal_local]] = True
        rec = score_unit(wanted_genes, lib, mask)
        n = int(len(sub))
        n_tnk = int(tnk_local.sum())
        rows.append(_finish(rec, "GSE205335", patient, n, n_tnk, rds_gz.name))
        print(
            f"  {patient}: cells={n} mal={rec['n_malignant']} tnk={n_tnk} "
            f"cldn4%={rec['cldn4_pct_gt0']:.4f} tacstd2%={rec['tacstd2_pct_gt0']:.4f}",
            flush=True,
        )
    return rows


def _compare_locked(units: pd.DataFrame) -> pd.DataFrame:
    ref = pd.read_csv(DATA / "locked" / "elf3_patient_scores.tsv", sep="\t")
    ref["unit_id"] = ref["unit_id"].astype(str)
    rows = []
    for cohort, ids in LOCKED_IDS.items():
        mine = units[units["cohort"] == cohort].set_index("unit_id")
        theirs = ref[ref["dataset"] == cohort].set_index("unit_id")
        missing = sorted(set(ids) - set(mine.index))
        extra = sorted(set(mine.index) - set(ids))
        both = [i for i in ids if i in mine.index and i in theirs.index]
        pct_m = mine.loc[both, "cldn4_pct_gt0"].to_numpy(float) * 100.0
        pct_t = theirs.loc[both, "pct_CLDN4"].to_numpy(float)
        mean_m = mine.loc[both, "cldn4_mean_log1p_raw"].to_numpy(float)
        mean_t = theirs.loc[both, "mean_CLDN4"].to_numpy(float)
        tac_m = mine.loc[both, "tacstd2_pct_gt0"].to_numpy(float) * 100.0
        tac_t = theirs.loc[both, "pct_TACSTD2"].to_numpy(float)
        frac_m = mine.loc[both, "frac_tnk"].to_numpy(float)
        frac_t = theirs.loc[both, "frac_tnk"].to_numpy(float)
        nmal_m = mine.loc[both, "n_malignant"].to_numpy(float)
        nmal_t = theirs.loc[both, "n_malignant"].to_numpy(float)
        rows.append(
            {
                "cohort": cohort,
                "n": len(both),
                "missing": ",".join(missing),
                "extra": ",".join(extra),
                "max_abs_cldn4_pct": float(np.max(np.abs(pct_m - pct_t))) if both else float("nan"),
                "max_abs_tacstd2_pct": float(np.max(np.abs(tac_m - tac_t))) if both else float("nan"),
                "max_abs_cldn4_mean_vs_log1p_raw": float(np.max(np.abs(mean_m - mean_t))) if both else float("nan"),
                "max_abs_frac_tnk": float(np.max(np.abs(frac_m - frac_t))) if both else float("nan"),
                "max_abs_n_mal": float(np.max(np.abs(nmal_m - nmal_t))) if both else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    print("GSE123902", flush=True)
    rows.extend(gse123902())
    print("GSE189357", flush=True)
    rows.extend(gse189357())
    print("GSE205335", flush=True)
    rows.extend(gse205335())
    print("GSE131907", flush=True)
    rows.extend(gse131907())
    units = pd.DataFrame(rows)
    if len(units) != 65:
        raise SystemExit(f"expected 65 units, got {len(units)}")
    units.to_csv(OUT / "unit_scores.tsv", sep="\t", index=False)
    comp = _compare_locked(units)
    comp.to_csv(OUT / "locked_reproduction.tsv", sep="\t", index=False)
    print(comp.to_string(index=False), flush=True)
    bad = comp[
        (comp["max_abs_cldn4_pct"] > 1e-6)
        | (comp["max_abs_tacstd2_pct"] > 1e-6)
        | (comp["max_abs_frac_tnk"] > 1e-8)
        | (comp["max_abs_n_mal"] > 0)
        | (comp["missing"] != "")
    ]
    prov = {
        "n_units": int(len(units)),
        "by_cohort": units.groupby("cohort").size().astype(int).to_dict(),
        "locked_match_ok": bool(bad.empty),
        "sources": {
            "GSE123902": "GEO GSE123902_RAW.tar dense UMI CSV; marker malignant",
            "GSE131907": "GEO raw UMI matrix + author malignant / T/NK labels",
            "GSE205335": "GEO dgCMatrix RDS + author malignant / T/NK labels",
            "GSE189357": "GEO 10x MTX; marker malignant",
        },
    }
    (OUT / "extract_provenance.json").write_text(json.dumps(prov, indent=2))
    if not bad.empty:
        raise SystemExit("locked reproduction failed:\n" + bad.to_string(index=False))
    print(f"wrote {OUT / 'unit_scores.tsv'} n={len(units)}", flush=True)


if __name__ == "__main__":
    main()
