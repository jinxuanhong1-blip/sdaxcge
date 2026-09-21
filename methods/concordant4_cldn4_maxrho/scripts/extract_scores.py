#!/usr/bin/env python3
"""Cell-level CLDN4 / keratin scores for the locked concordant-4 units.

Cohorts (not expanded): GSE123902, GSE131907, GSE205335, GSE189357.
Malignant and T/NK gates match the locked concordant-4 tables.

Scores written per unit:
  percent-positive at raw-UMI cuts, mean log1p(UMI), mean log1p(CP10K),
  keratin means on both scales, and the malignant CLDN4 UMI histogram
  (for cohort-quantile cuts in the sweep).
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
SCORE_GENES = ["CLDN4", "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT6A"]
PANEL = list(dict.fromkeys(SCORE_GENES + ["PTPRC"] + MARKER_TNK))
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}


def _upper_gene(name: str) -> str:
    return str(name).upper().split(".")[0].strip()


def _as_int_umi(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=np.float64)
    rounded = np.rint(arr)
    if arr.size and np.nanmax(np.abs(arr - rounded)) > 1e-3:
        raise ValueError("expression is not integer UMI")
    return rounded


def _hist(arr: np.ndarray) -> dict[str, int]:
    if arr.size == 0:
        return {}
    vals, cts = np.unique(arr.astype(np.int64), return_counts=True)
    return {str(int(v)): int(c) for v, c in zip(vals, cts)}


def _mean_log1p(x: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    return float(np.log1p(x).mean())


def _mean_log1p_cp10k(x: np.ndarray, lib: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    lib = np.maximum(lib.astype(np.float64), 1.0)
    return float(np.log1p(x / lib * 1e4).mean())


def _mean_log1p_pos(x: np.ndarray) -> float:
    pos = x[x > 0]
    if pos.size == 0:
        return float("nan")
    return float(np.log1p(pos).mean())


def _p90_log1p(x: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    return float(np.quantile(np.log1p(x), 0.9))


def _p90_log1p_cp10k(x: np.ndarray, lib: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    lib = np.maximum(lib.astype(np.float64), 1.0)
    return float(np.quantile(np.log1p(x / lib * 1e4), 0.9))


def _pct(x: np.ndarray, cut: float, ge: bool) -> float:
    if x.size == 0:
        return float("nan")
    if ge:
        return float(np.mean(x >= cut))
    return float(np.mean(x > cut))


def score_malignant(genes: dict[str, np.ndarray], lib: np.ndarray, mask: np.ndarray) -> dict:
    """genes: name -> full-cell vector aligned with lib and mask."""
    x = _as_int_umi(genes["CLDN4"][mask])
    L = lib[mask]
    row = {
        "n_malignant": int(mask.sum()),
        "pct_gt0": _pct(x, 0, ge=False),
        "pct_ge2": _pct(x, 2, ge=True),
        "pct_ge3": _pct(x, 3, ge=True),
        "pct_ge5": _pct(x, 5, ge=True),
        "pct_ge10": _pct(x, 10, ge=True),
        "mean_log1p_raw": _mean_log1p(x),
        "mean_log1p_cp10k": _mean_log1p_cp10k(x, L),
        "mean_log1p_pos_raw": _mean_log1p_pos(x),
        "mean_log1p_pos_cp10k": _mean_log1p_pos_cp10k(x, L),
        "p90_log1p_raw": _p90_log1p(x),
        "p90_log1p_cp10k": _p90_log1p_cp10k(x, L),
        "cldn4_hist": json.dumps(_hist(x), separators=(",", ":")),
    }
    for gene in SCORE_GENES:
        if gene == "CLDN4":
            continue
        if gene not in genes:
            row[f"{gene}_raw"] = float("nan")
            row[f"{gene}_cp10k"] = float("nan")
            continue
        g = _as_int_umi(genes[gene][mask])
        row[f"{gene}_raw"] = _mean_log1p(g)
        row[f"{gene}_cp10k"] = _mean_log1p_cp10k(g, L)
    return row


def _mean_log1p_pos_cp10k(x: np.ndarray, lib: np.ndarray) -> float:
    pos = x > 0
    if not np.any(pos):
        return float("nan")
    return _mean_log1p_cp10k(x[pos], lib[pos])


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


def gse123902() -> list[dict]:
    tar_path = GEO / "GSE123902" / "GSE123902_RAW.tar"
    rows = []
    with tarfile.open(tar_path) as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".csv.gz")]
        for m in members:
            name = Path(m.name).name
            parts = name.split("_")
            donor = parts[2]
            if "NORMAL" in name:
                tissue = "NORMAL"
            elif "METASTASIS" in name:
                tissue = "METASTASIS"
            elif "PRIMARY" in name:
                tissue = "PRIMARY"
            else:
                tissue = "OTHER"
            if tissue not in {"PRIMARY", "METASTASIS"}:
                continue
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            df = pd.read_csv(raw, index_col=0)
            lib = df.sum(axis=1).to_numpy(dtype=np.float64)
            df = _collapse_columns(df)
            missing = [g for g in PANEL if g not in df.columns]
            if "CLDN4" in missing:
                raise SystemExit(f"{name} missing CLDN4; also missing {missing}")
            genes = {g: df[g].to_numpy(dtype=np.float64) for g in PANEL if g in df.columns}
            mal, tnk = marker_masks(genes, len(df))
            rec = score_malignant(genes, lib, mal)
            rec.update(
                {
                    "cohort": "GSE123902",
                    "unit_id": donor,
                    "histology": "LUAD",
                    "tissue": tissue,
                    "n_cells": int(len(df)),
                    "n_tnk": int(tnk.sum()),
                    "frac_tnk": float(tnk.mean()) if len(df) else float("nan"),
                    "file": name,
                }
            )
            rows.append(rec)
            print(
                f"  {name}: cells={rec['n_cells']} mal={rec['n_malignant']} tnk={rec['n_tnk']}",
                flush=True,
            )
    frame = pd.DataFrame(rows)
    frame = frame.sort_values(["unit_id", "tissue"]).drop_duplicates("unit_id", keep="first")
    frame = frame[(frame["n_malignant"] >= 20) & (frame["n_tnk"] >= 20)].copy()
    print(f"GSE123902 locked-style donors n={len(frame)}", flush=True)
    return frame.to_dict(orient="records")


def _parse_189357_histology() -> dict[str, str]:
    out = {}
    for path in sorted((DATA / "geo_meta").glob("GSM5699*.txt")):
        patient = None
        hist = None
        for line in path.read_text().splitlines():
            if line.startswith("!Sample_title = "):
                patient = line.split("=", 1)[1].strip().split()[0]
            if "histolgical type:" in line or "histological type:" in line:
                hist = line.split(":", 1)[1].strip()
        if not patient or not hist:
            raise SystemExit(f"histology parse failed for {path.name}")
        out[patient] = hist
    if set(out) != {f"TD{i}" for i in range(1, 10)}:
        raise SystemExit(f"unexpected GSE189357 histology keys: {sorted(out)}")
    return out


def _tar_members(tf: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    return {Path(m.name).name: m for m in tf.getmembers()}


def _read_gz_lines(tf: tarfile.TarFile, member: tarfile.TarInfo):
    with gzip.GzipFile(fileobj=tf.extractfile(member)) as handle:
        for line in handle:
            yield line.decode()


def gse189357() -> list[dict]:
    hist_map = _parse_189357_histology()
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
            if "CLDN4" not in counts:
                raise SystemExit(f"{sample} missing CLDN4")
            mal, tnk = marker_masks(counts, n_cells)
            rec = score_malignant(counts, lib, mal)
            rec.update(
                {
                    "cohort": "GSE189357",
                    "unit_id": sample,
                    "histology": hist_map[sample],
                    "tissue": "TUMOR",
                    "n_cells": n_cells,
                    "n_tnk": int(tnk.sum()),
                    "frac_tnk": float(tnk.sum() / n_cells),
                    "file": mtx_name,
                }
            )
            rows.append(rec)
            print(
                f"  {sample} {rec['histology']}: cells={n_cells} mal={rec['n_malignant']} tnk={rec['n_tnk']}",
                flush=True,
            )
    frame = pd.DataFrame(rows)
    frame = frame[(frame["n_malignant"] >= 20) & (frame["n_tnk"] >= 20)].copy()
    return frame.to_dict(orient="records")


def gse131907() -> list[dict]:
    ann = pd.read_csv(
        GEO / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t"
    )
    matrix = GEO / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    print(f"streaming {matrix.name}", flush=True)
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        n_cells = len(cells)
        lib = np.zeros(n_cells, dtype=np.float64)
        n_rows = 0
        for line in handle:
            tab = line.find("\t")
            gene = _upper_gene(line[:tab])
            arr = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise SystemExit(f"{gene}: {arr.size} values, expected {n_cells}")
            lib += arr
            n_rows += 1
            if gene in PANEL:
                if gene in found:
                    found[gene] = found[gene] + arr
                else:
                    found[gene] = arr.astype(np.float64, copy=False)
            if n_rows % 2000 == 0:
                print(f"  streamed {n_rows} genes", flush=True)
    print(f"scanned {n_rows} genes; panel found {sorted(found)}", flush=True)
    if "CLDN4" not in found:
        raise SystemExit("GSE131907 matrix has no CLDN4")
    ann = ann.set_index("Index").reindex(cells)
    if ann["Sample"].isna().any():
        raise SystemExit("GSE131907 barcodes do not match annotation Index")
    mal = ann["Cell_subtype"].eq("Malignant cells").to_numpy()
    tnk = ann["Cell_type"].isin(["T lymphocytes", "NK cells"]).to_numpy()
    sample = ann["Sample"].astype(str).to_numpy()
    origin = ann["Sample_Origin"].astype(str).to_numpy()
    rows = []
    for sid in pd.unique(sample):
        sm = sample == sid
        ori = str(origin[sm][0])
        if ori not in TUMOR_ORIGINS:
            continue
        mal_m = sm & mal
        if int(mal_m.sum()) < 20:
            continue
        rec = score_malignant(found, lib, mal_m)
        n = int(sm.sum())
        n_tnk = int((sm & tnk).sum())
        rec.update(
            {
                "cohort": "GSE131907",
                "unit_id": str(sid),
                "histology": "LUAD",
                "tissue": ori,
                "n_cells": n,
                "n_tnk": n_tnk,
                "frac_tnk": n_tnk / n,
                "file": matrix.name,
            }
        )
        rows.append(rec)
    print(f"GSE131907 tumor samples n_mal>=20: {len(rows)}", flush=True)
    return rows


def gse205335() -> list[dict]:
    import rdata

    meta = pd.read_csv(DATA / "GSE205335_gsm_metadata.csv")
    tumor = meta[~meta["tissue"].astype(str).str.startswith("Normal")].copy()
    subtype = tumor.groupby("patient")["cancer_subtype"].nunique()
    if int(subtype.max()) != 1:
        raise SystemExit("a GSE205335 patient has conflicting subtypes")
    hist = tumor.groupby("patient")["cancer_subtype"].first().to_dict()
    tissue = tumor.groupby("patient")["tissue"].apply(lambda s: ",".join(sorted(set(s)))).to_dict()
    keep_ident = set(tumor["orig.ident"].astype(str))
    ident = pd.read_csv(
        GEO / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t"
    )
    ident["barcode"] = ident["barcode"].astype(str)
    ident["orig.ident"] = ident["orig.ident"].astype(str)
    ident = ident[ident["orig.ident"].isin(keep_ident)].copy()
    ident = ident.merge(tumor[["orig.ident", "patient"]], on="orig.ident", how="left")
    if ident["patient"].isna().any():
        raise SystemExit("GSE205335 identity rows missing patient")

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
        raise SystemExit(f"RDS missing dgCMatrix fields: {sorted(vars(obj))[:20]}")
    from scipy import sparse

    genes = np.array([_upper_gene(g) for g in np.asarray(obj.Dimnames[0], dtype=str)])
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    lib = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    wanted = {}
    for gene in PANEL:
        hits = np.flatnonzero(genes == gene)
        if hits.size == 0:
            continue
        acc = np.zeros(matrix.shape[1], dtype=np.float64)
        for h in hits:
            acc += np.asarray(matrix.getrow(int(h)).todense()).ravel()
        wanted[gene] = acc
    del obj, matrix
    if "CLDN4" not in wanted:
        raise SystemExit("GSE205335 RDS has no CLDN4")
    bc_index = pd.Index(barcodes.astype(str))
    missing = int((~ident["barcode"].isin(bc_index)).sum())
    if missing:
        raise SystemExit(f"GSE205335 identity barcodes missing from RDS: {missing}")
    pos = bc_index.get_indexer(ident["barcode"])
    ident = ident.copy()
    ident["pos"] = pos
    rows = []
    for patient, sub in ident.groupby("patient"):
        idx = sub["pos"].to_numpy()
        mal_local = sub["lineage.sub"].eq("Malignant cells").to_numpy()
        tnk_local = sub["lineage.total"].eq("T/NK cells").to_numpy()
        mask = np.zeros(len(barcodes), dtype=bool)
        mask[idx[mal_local]] = True
        if int(mask.sum()) < 1:
            continue
        rec = score_malignant(wanted, lib, mask)
        n = int(len(sub))
        n_tnk = int(tnk_local.sum())
        rec.update(
            {
                "cohort": "GSE205335",
                "unit_id": str(patient),
                "histology": str(hist[patient]),
                "tissue": tissue[patient],
                "n_cells": n,
                "n_tnk": n_tnk,
                "frac_tnk": n_tnk / n if n else float("nan"),
                "file": rds_gz.name,
            }
        )
        rows.append(rec)
        print(
            f"  {patient} {rec['histology']}: cells={n} mal={rec['n_malignant']} tnk={rec['n_tnk']}",
            flush=True,
        )
    return rows


def _locked_compare(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    specs = [
        ("GSE123902", "GSE123902_marker_units.tsv", "patient", "mal_CLDN4_pct", "mal_CLDN4_mean", False),
        ("GSE131907", "GSE131907_samples.tsv", "sample", "mal_CLDN4_pct", "mal_CLDN4_mean", True),
        ("GSE205335", "GSE205335_patients.tsv", "patient", "mal_CLDN4_pct_pos", "mal_CLDN4_mean", True),
        ("GSE189357", "GSE189357_marker_units.tsv", "patient", "mal_CLDN4_pct", "mal_CLDN4_mean", False),
    ]
    for cohort, fname, idcol, pctcol, meancol, pct_is_percent in specs:
        locked = pd.read_csv(DATA / "locked" / fname, sep="\t")
        if cohort == "GSE123902":
            locked = locked[locked["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
            locked = locked.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
            locked = locked[locked["eligible"].astype(str).str.lower() == "true"]
        elif cohort == "GSE131907":
            locked = locked[locked["origin"].isin(TUMOR_ORIGINS) & (locked["n_malignant"] >= 20)]
        elif cohort == "GSE189357":
            locked = locked[locked["eligible"].astype(str).str.lower() == "true"]
        mine = units[units["cohort"] == cohort].set_index("unit_id")
        locked = locked.set_index(idcol)
        ids = sorted(set(mine.index) & set(locked.index.astype(str)))
        only_mine = sorted(set(mine.index) - set(locked.index.astype(str)))
        only_locked = sorted(set(locked.index.astype(str)) - set(mine.index))
        pct_locked = locked.loc[ids, pctcol].to_numpy(float)
        if pct_is_percent:
            pct_locked = pct_locked / 100.0
        pct_mine = mine.loc[ids, "pct_gt0"].to_numpy(float)
        mean_locked = locked.loc[ids, meancol].to_numpy(float)
        raw = mine.loc[ids, "mean_log1p_raw"].to_numpy(float)
        cp = mine.loc[ids, "mean_log1p_cp10k"].to_numpy(float)
        frac_l = locked.loc[ids, "frac_tnk"].to_numpy(float)
        frac_m = mine.loc[ids, "frac_tnk"].to_numpy(float)
        nmal_l = locked.loc[ids, "n_malignant"].to_numpy(float)
        nmal_m = mine.loc[ids, "n_malignant"].to_numpy(float)
        rows.append(
            {
                "cohort": cohort,
                "n_overlap": len(ids),
                "only_mine": ",".join(only_mine),
                "only_locked": ",".join(only_locked),
                "max_abs_pct": float(np.max(np.abs(pct_mine - pct_locked))) if ids else float("nan"),
                "max_abs_frac_tnk": float(np.max(np.abs(frac_m - frac_l))) if ids else float("nan"),
                "max_abs_n_mal": float(np.max(np.abs(nmal_m - nmal_l))) if ids else float("nan"),
                "max_abs_mean_vs_raw": float(np.max(np.abs(raw - mean_locked))) if ids else float("nan"),
                "max_abs_mean_vs_cp10k": float(np.max(np.abs(cp - mean_locked))) if ids else float("nan"),
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
    print("GSE131907", flush=True)
    rows.extend(gse131907())
    print("GSE205335", flush=True)
    rows.extend(gse205335())
    units = pd.DataFrame(rows)
    units.to_csv(OUT / "unit_scores.tsv", sep="\t", index=False)
    comp = _locked_compare(units)
    comp.to_csv(OUT / "locked_reproduction.tsv", sep="\t", index=False)
    print(comp.to_string(index=False), flush=True)
    prov = {
        "n_units": int(len(units)),
        "by_cohort": units.groupby("cohort").size().astype(int).to_dict(),
        "sources": {
            "GSE123902": "GEO GSE123902_RAW.tar dense UMI CSV",
            "GSE131907": "GEO GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz + cell annotation",
            "GSE205335": "GEO GSE205335_Lung_IO_UMI_matrix.rds.gz + CellIdentity",
            "GSE189357": "GEO GSE189357_RAW.tar 10x MTX",
        },
    }
    (OUT / "extract_provenance.json").write_text(json.dumps(prov, indent=2))
    print(f"wrote {OUT / 'unit_scores.tsv'} n={len(units)}", flush=True)


if __name__ == "__main__":
    main()
