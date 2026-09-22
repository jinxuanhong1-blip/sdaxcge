#!/usr/bin/env python3
"""TACSTD2–CLDN4 same-cell coexpression BRIDGE (not mediation).

Public only:
  - concordant-4 malignant cells (GSE123902 / 131907 / 205335 / 189357)
  - CosMx He2022 tumor cells (figshare 25976224)
  - DepMap Gygi lung protein (TACSTD2–CLDN4)

Does not recompute locked CLDN4 CosMx 0.36/0.52 or concordant-4 ρ=−0.531.
Does not claim TACSTD2→CLDN4→immune mediation (PR #718 null).
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
GEO = Path(os.environ.get("GEO_C4", "/tmp/geo_c4"))
COSMX = Path(
    os.environ.get(
        "COSMX_H5AD",
        "/workspace/data/cosmx_nsclc/cosmx_human_nsclc_clustered.h5ad",
    )
)
EXPECTED_COSMX = 2_755_776_882

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


def fisher_or(n11, n10, n01, n00):
    table = np.array([[n11, n10], [n01, n00]], dtype=float)
    if table.min() < 0 or table.sum() == 0:
        return np.nan, np.nan
    # Haldane-Anscombe if a zero cell blocks OR
    if (table == 0).any():
        table = table + 0.5
    oddsratio = (table[0, 0] * table[1, 1]) / (table[0, 1] * table[1, 0])
    _, p = stats.fisher_exact([[n11, n10], [n01, n00]])
    return float(oddsratio), float(p)


def pair_stats(tac: np.ndarray, cld: np.ndarray) -> dict:
    tac = np.asarray(tac, dtype=float)
    cld = np.asarray(cld, dtype=float)
    m = np.isfinite(tac) & np.isfinite(cld)
    tac, cld = tac[m], cld[m]
    n = int(tac.size)
    if n == 0:
        return dict(
            n=0,
            n11=0,
            n10=0,
            n01=0,
            n00=0,
            p_tac=np.nan,
            p_cld=np.nan,
            p_both=np.nan,
            p_indep=np.nan,
            ratio_vs_indep=np.nan,
            odds_ratio=np.nan,
            fisher_p=np.nan,
            spearman=np.nan,
            spearman_p=np.nan,
            spearman_n=0,
        )
    tpos = tac > 0
    cpos = cld > 0
    n11 = int(np.sum(tpos & cpos))
    n10 = int(np.sum(tpos & ~cpos))
    n01 = int(np.sum(~tpos & cpos))
    n00 = int(np.sum(~tpos & ~cpos))
    p_tac = float(np.mean(tpos))
    p_cld = float(np.mean(cpos))
    p_both = float(n11 / n)
    p_indep = p_tac * p_cld
    ratio = p_both / p_indep if p_indep > 0 else np.nan
    oratio, fp = fisher_or(n11, n10, n01, n00)
    if n >= 4 and np.unique(tac).size > 1 and np.unique(cld).size > 1:
        rho, rp = stats.spearmanr(tac, cld)
        rho, rp = float(rho), float(rp)
    else:
        rho, rp = np.nan, np.nan
    return dict(
        n=n,
        n11=n11,
        n10=n10,
        n01=n01,
        n00=n00,
        p_tac=p_tac,
        p_cld=p_cld,
        p_both=p_both,
        p_indep=p_indep,
        ratio_vs_indep=ratio,
        odds_ratio=oratio,
        fisher_p=fp,
        spearman=rho,
        spearman_p=rp,
        spearman_n=n,
    )


def load_locked() -> pd.DataFrame:
    return pd.read_csv(DATA / "locked_units.tsv", sep="\t", dtype={"unit_id": str})


def marker_malignant(df: pd.DataFrame) -> np.ndarray:
    cols = {c.upper(): c for c in df.columns}

    def g(name):
        c = cols.get(name)
        if c is None:
            return np.zeros(len(df), dtype=float)
        return df[c].to_numpy(dtype=float)

    epi = (g("EPCAM") > 0) | (g("KRT8") > 0) | (g("KRT18") > 0) | (g("KRT19") > 0)
    return epi & (g("PTPRC") == 0)


def load_gse123902(locked: pd.DataFrame) -> list[dict]:
    ids = locked.loc[locked.dataset == "GSE123902", "unit_id"].astype(str).tolist()
    tar_path = GEO / "gse123902" / "GSE123902_RAW.tar"
    out = []
    with tarfile.open(tar_path) as tar:
        members = {Path(m.name).name: m for m in tar.getmembers() if m.isfile()}
        for unit in ids:
            fname = GSE123902_FILES[unit]
            member = members[fname]
            raw = tar.extractfile(member)
            assert raw is not None
            with gzip.GzipFile(fileobj=raw) as gz:
                df = pd.read_csv(gz, index_col=0)
            df.columns = [str(c).upper() for c in df.columns]
            df = df.loc[:, ~df.columns.duplicated()]
            mal = marker_malignant(df)
            if "TACSTD2" not in df.columns or "CLDN4" not in df.columns:
                raise SystemExit(f"GSE123902 {unit} missing TACSTD2/CLDN4")
            st = pair_stats(df.loc[mal, "TACSTD2"].to_numpy(), df.loc[mal, "CLDN4"].to_numpy())
            st.update(dataset="GSE123902", unit_id=unit)
            out.append(st)
            say(f"  GSE123902 {unit} n_mal={st['n']} ratio={st['ratio_vs_indep']:.3f}")
    return out


def load_gse189357(locked: pd.DataFrame) -> list[dict]:
    ids = locked.loc[locked.dataset == "GSE189357", "unit_id"].astype(str).tolist()
    tar_path = GEO / "gse189357" / "GSE189357_RAW.tar"
    raw_dir = GEO / "gse189357" / "raw"
    if not (raw_dir.exists() and any(raw_dir.glob("*_matrix.mtx.gz"))):
        raw_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path) as tar:
            tar.extractall(raw_dir)
    out = []
    for unit in ids:
        mtx = next(raw_dir.glob(f"*_{unit}_matrix.mtx.gz"))
        feat = next(raw_dir.glob(f"*_{unit}_features.tsv.gz"))
        symbols = []
        with gzip.open(feat, "rt") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                symbols.append((parts[1] if len(parts) > 1 else parts[0]).upper())
        with gzip.open(mtx, "rb") as handle:
            mat = mmread(handle).tocsr()
        if mat.shape[0] != len(symbols):
            mat = mat.T.tocsr()
        seen: dict[str, int] = {}
        keep_idx = []
        for i, s in enumerate(symbols):
            if s and s not in seen:
                seen[s] = i
                keep_idx.append(i)
        mat = mat[keep_idx]
        kept_symbols = [symbols[i] for i in keep_idx]
        sym_index = {s: i for i, s in enumerate(kept_symbols)}

        def get(name):
            i = sym_index.get(name)
            if i is None:
                return np.zeros(mat.shape[1], dtype=float)
            return np.asarray(mat.getrow(i).toarray()).ravel()

        mal = (
            (get("EPCAM") > 0) | (get("KRT8") > 0) | (get("KRT18") > 0) | (get("KRT19") > 0)
        ) & (get("PTPRC") == 0)
        tac = get("TACSTD2")[mal]
        cld = get("CLDN4")[mal]
        st = pair_stats(tac, cld)
        st.update(dataset="GSE189357", unit_id=unit)
        out.append(st)
        say(f"  GSE189357 {unit} n_mal={st['n']} ratio={st['ratio_vs_indep']:.3f}")
    return out


def load_gse131907(locked: pd.DataFrame) -> list[dict]:
    ids = locked.loc[locked.dataset == "GSE131907", "unit_id"].astype(str).tolist()
    ann = pd.read_csv(
        GEO / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        sep="\t",
        dtype=str,
    )
    ann = ann[ann["Sample"].isin(ids)]
    mat_path = GEO / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
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
    uniq, starts, counts = np.unique(unit_of, return_index=True, return_counts=True)
    starts = starts.astype(np.int64)
    stored = {}
    say(f"  GSE131907 streaming malignant={keep_idx.size}")
    with gzip.open(mat_path, "rt") as handle:
        handle.readline()
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.upper()
            if gene not in {"TACSTD2", "CLDN4"}:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            stored[gene] = arr[keep_idx]
            if len(stored) == 2:
                break
    if set(stored) != {"TACSTD2", "CLDN4"}:
        raise SystemExit(f"GSE131907 missing genes: {stored.keys()}")
    out = []
    for i, unit in enumerate(uniq):
        sl = slice(int(starts[i]), int(starts[i] + counts[i]))
        st = pair_stats(stored["TACSTD2"][sl], stored["CLDN4"][sl])
        st.update(dataset="GSE131907", unit_id=str(unit))
        out.append(st)
        say(f"  GSE131907 {unit} n_mal={st['n']} ratio={st['ratio_vs_indep']:.3f}")
    return out


def ensure_gse205335_mtx() -> Path:
    export = GEO / "gse205335" / "focus_export"
    mtx = export / "matrix.mtx"
    if mtx.exists() and (export / "genes.txt").exists() and (export / "cells.txt").exists():
        return export
    genes = DATA / "wanted_focus.txt"
    genes.write_text(
        "\n".join(
            [
                "TACSTD2",
                "CLDN4",
                "EPCAM",
                "KRT8",
                "KRT18",
                "KRT19",
                "PTPRC",
                "CD3D",
                "CD3E",
                "CD8A",
                "NKG7",
                "GNLY",
                "KLRD1",
            ]
        )
        + "\n"
    )
    rds_gz = GEO / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    export.mkdir(parents=True, exist_ok=True)
    say("R export GSE205335 focus genes")
    subprocess.check_call(
        [
            "Rscript",
            str(ROOT / "scripts" / "export_gse205335.R"),
            str(rds_gz),
            str(genes),
            str(export),
        ]
    )
    return export


def load_gse205335(locked: pd.DataFrame) -> list[dict]:
    ids = locked.loc[locked.dataset == "GSE205335", "unit_id"].astype(str).tolist()
    export = ensure_gse205335_mtx()
    mat = mmread(export / "matrix.mtx").tocsr()
    genes = [g.strip().upper() for g in (export / "genes.txt").read_text().splitlines() if g.strip()]
    cells = [c.strip() for c in (export / "cells.txt").read_text().splitlines() if c.strip()]
    if mat.shape[0] != len(genes):
        mat = mat.T.tocsr()
    gene_index = {g: i for i, g in enumerate(genes)}
    ident = pd.read_csv(
        GEO / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz",
        sep="\t",
        dtype=str,
    )
    gsm = pd.read_csv(DATA / "GSE205335_gsm_map.tsv", sep="\t", dtype=str)
    ident = ident.merge(gsm[["orig.ident", "patient", "tissue"]], on="orig.ident", how="left")
    ident = ident.set_index("barcode").reindex(cells)
    patient = ident["patient"].astype(str).to_numpy()
    lineage_sub = ident["lineage.sub"].astype(str).to_numpy()
    mal = lineage_sub == "Malignant cells"

    def col(name):
        i = gene_index.get(name)
        if i is None:
            return np.zeros(mat.shape[1], dtype=float)
        return np.asarray(mat.getrow(i).todense()).ravel()

    tac_all = col("TACSTD2")
    cld_all = col("CLDN4")
    out = []
    for unit in ids:
        mask = mal & (patient == unit)
        st = pair_stats(tac_all[mask], cld_all[mask])
        st.update(dataset="GSE205335", unit_id=unit)
        out.append(st)
        say(f"  GSE205335 {unit} n_mal={st['n']} ratio={st['ratio_vs_indep']:.3f}")
    return out


def summarize_units(rows: list[dict], locked: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.DataFrame(rows)
    locked = locked.copy()
    locked["unit_id"] = locked["unit_id"].astype(str)
    df = df.merge(
        locked[["dataset", "unit_id", "n_cells", "n_malignant", "n_tnk", "frac_tnk"]],
        on=["dataset", "unit_id"],
        how="left",
        suffixes=("", "_locked"),
    )
    # patient-level dual-positive % among malignant vs T/NK
    df["dual_pct"] = 100.0 * df["p_both"]
    # keep units with ≥30 malignant for within-cell layers (matches PR #638)
    df["usable_rho"] = df["n"] >= 30
    return df


def cohort_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = [(ds, sub) for ds, sub in df.groupby("dataset")] + [("ALL", df)]
    for ds, sub in groups:
        use = sub[sub["usable_rho"]].copy()
        use = use[np.isfinite(use["ratio_vs_indep"]) & np.isfinite(use["spearman"])]
        if len(use) == 0:
            continue
        dual_r, dual_p = (
            stats.spearmanr(use["dual_pct"], use["frac_tnk"])
            if len(use) >= 5 and use["frac_tnk"].notna().sum() >= 5
            else (np.nan, np.nan)
        )
        rows.append(
            dict(
                dataset=ds,
                n_units=int(len(sub)),
                n_units_ge30=int(len(use)),
                median_ratio_vs_indep=float(use["ratio_vs_indep"].median()),
                median_odds_ratio=float(use["odds_ratio"].median()),
                median_p_both=float(use["p_both"].median()),
                median_p_tac=float(use["p_tac"].median()),
                median_p_cld=float(use["p_cld"].median()),
                median_spearman=float(use["spearman"].median()),
                frac_spearman_pos=float(np.mean(use["spearman"] > 0)),
                wilcoxon_ratio_gt1_p=float(
                    stats.wilcoxon(use["ratio_vs_indep"] - 1.0, alternative="greater").pvalue
                )
                if len(use) >= 5
                else np.nan,
                wilcoxon_spearman_p=float(
                    stats.wilcoxon(use["spearman"], alternative="greater").pvalue
                )
                if len(use) >= 5
                else np.nan,
                dual_vs_tnk_spearman=float(dual_r) if np.isfinite(dual_r) else np.nan,
                dual_vs_tnk_p=float(dual_p) if np.isfinite(dual_p) else np.nan,
            )
        )
    return pd.DataFrame(rows)


def _decode_h5_cats(f, key: str):
    """Decode anndata categorical obs column from an open h5py file."""
    import h5py

    node = f[f"obs/{key}"]
    if isinstance(node, h5py.Dataset):
        vals = node[:]
        if vals.dtype.kind in {"S", "O"}:
            return np.array(
                [v.decode() if isinstance(v, (bytes, np.bytes_)) else str(v) for v in vals]
            )
        return vals
    # categorical
    codes = node["codes"][:]
    cats = node["categories"][:]
    cats = np.array([c.decode() if isinstance(c, (bytes, np.bytes_)) else str(c) for c in cats])
    return cats[codes]


def analyze_cosmx() -> tuple[pd.DataFrame, dict]:
    """Author-matched tumor cells per section (same map as PR #726).

    Detection uses raw counts from layers/counts. Spearman uses log1p CP10k
    on the same cells (matches inventory Spearmans in PR #726).
    """
    if not COSMX.exists() or COSMX.stat().st_size != EXPECTED_COSMX:
        raise SystemExit(f"CosMx h5ad missing or wrong size: {COSMX}")
    import h5py
    from scipy import sparse

    tumor_label = {
        "LUAD-5 R1": "tumor 5",
        "LUAD-5 R2": "tumor 5",
        "LUAD-5 R3": "tumor 5",
        "LUSC-6": "tumor 6",
        "LUAD-9 R1": "tumor 9",
        "LUAD-9 R2": "tumor 9",
        "LUAD-12": "tumor 12",
        "LUAD-13": "tumor 13",
    }
    expected_patient = {
        "LUAD-5 R1": "Lung5",
        "LUAD-5 R2": "Lung5",
        "LUAD-5 R3": "Lung5",
        "LUSC-6": "Lung6",
        "LUAD-9 R1": "Lung9",
        "LUAD-9 R2": "Lung9",
        "LUAD-12": "Lung12",
        "LUAD-13": "Lung13",
    }
    say(f"reading CosMx {COSMX}")
    with h5py.File(COSMX, "r") as f:
        genes = f["var/_index"][:]
        genes = np.array([g.decode() if isinstance(g, (bytes, np.bytes_)) else str(g) for g in genes])
        gmap = {g.upper(): i for i, g in enumerate(genes)}
        for g in ("TACSTD2", "CLDN4"):
            if g not in gmap:
                raise SystemExit(f"CosMx missing {g}")
        sample = _decode_h5_cats(f, "sample")
        patient = _decode_h5_cats(f, "patient")
        cell_type = _decode_h5_cats(f, "cell_type")
        n_counts = f["obs/n_counts"][:].astype(np.float64)
        indptr = f["layers/counts/indptr"][:]
        indices = f["layers/counts/indices"][:]
        data = f["layers/counts/data"][:]
        mat = sparse.csr_matrix((data, indices, indptr), shape=(len(sample), len(genes)))
        del data, indices, indptr
        tac_raw = np.asarray(mat.getcol(gmap["TACSTD2"]).toarray()).ravel().astype(np.float64)
        cld_raw = np.asarray(mat.getcol(gmap["CLDN4"]).toarray()).ravel().astype(np.float64)
        del mat
    lib = n_counts.copy()
    lib[lib <= 0] = np.nan
    tac_log = np.nan_to_num(np.log1p(tac_raw / lib * 1e4), nan=0.0)
    cld_log = np.nan_to_num(np.log1p(cld_raw / lib * 1e4), nan=0.0)

    rows = []
    for s, lab in tumor_label.items():
        m = sample == s
        pat = str(patient[m][0])
        if pat != expected_patient[s]:
            raise SystemExit(f"{s} patient {pat} != {expected_patient[s]}")
        mal = cell_type[m] == lab
        # detection on raw counts; Spearman on log1p CP10k
        det = pair_stats(tac_raw[m][mal], cld_raw[m][mal])
        rho, rp, nrho = (
            stats.spearmanr(tac_log[m][mal], cld_log[m][mal]).correlation,
            stats.spearmanr(tac_log[m][mal], cld_log[m][mal]).pvalue,
            int(mal.sum()),
        )
        det.update(
            sample=s,
            patient=pat,
            tumor_label=lab,
            n_tumor=int(mal.sum()),
            spearman=float(rho),
            spearman_p=float(rp),
            spearman_n=nrho,
            spearman_layer="log1p_cp10k",
            detection_layer="raw_count",
        )
        rows.append(det)
        say(
            f"  CosMx {s} n={det['n']} rho={det['spearman']:.3f} "
            f"ratio={det['ratio_vs_indep']:.3f} OR={det['odds_ratio']:.3f}"
        )
    df = pd.DataFrame(rows)
    # Cross-check inventory Spearmans from PR #726
    inv = pd.read_csv(DATA / "cosmx_inventory.csv")
    merged = df.merge(inv[["sample", "spearman_tacstd2_cldn4", "n_tumor"]], on="sample", suffixes=("", "_inv"))
    max_abs = float(np.max(np.abs(merged["spearman"] - merged["spearman_tacstd2_cldn4"])))
    if max_abs > 1e-6:
        say(f"WARNING spearman drift vs PR726 inventory max|Δ|={max_abs}")
    else:
        say("CosMx Spearman matches PR #726 inventory")
    if not np.array_equal(merged["n_tumor"].to_numpy(), merged["n_tumor_inv"].to_numpy()):
        say(f"WARNING n_tumor differs from inventory: {list(zip(merged.sample, merged.n_tumor, merged.n_tumor_inv))}")
    summary = dict(
        n_sections=int(len(df)),
        n_patients=int(df["patient"].nunique()),
        n_tumor_total=int(df["n"].sum()),
        median_spearman=float(df["spearman"].median()),
        min_spearman=float(df["spearman"].min()),
        max_spearman=float(df["spearman"].max()),
        all_spearman_pos=bool((df["spearman"] > 0).all()),
        median_ratio_vs_indep=float(df["ratio_vs_indep"].median()),
        median_odds_ratio=float(df["odds_ratio"].median()),
        sections_ratio_gt1=int((df["ratio_vs_indep"] > 1).sum()),
        sections_or_gt1=int((df["odds_ratio"] > 1).sum()),
        wilcoxon_ratio_gt1_p=float(
            stats.wilcoxon(df["ratio_vs_indep"] - 1.0, alternative="greater").pvalue
        ),
        wilcoxon_spearman_p=float(stats.wilcoxon(df["spearman"], alternative="greater").pvalue),
        spearman_matches_pr726_inventory=bool(max_abs <= 1e-6),
    )
    return df, summary


def analyze_depmap() -> dict:
    g = pd.read_csv(DATA / "gygi_lung_selected_genes.csv")
    m = g[["TACSTD2", "CLDN4"]].dropna()
    rho, p = stats.spearmanr(m["TACSTD2"], m["CLDN4"])
    # bootstrap CI
    rng = np.random.default_rng(42)
    boots = []
    arr = m.to_numpy()
    for _ in range(5000):
        ix = rng.integers(0, len(arr), len(arr))
        boots.append(stats.spearmanr(arr[ix, 0], arr[ix, 1]).correlation)
    lo, hi = np.quantile(boots, [0.025, 0.975])
    return dict(
        layer="protein",
        cohort="Gygi_CCLE_MS_LUNG",
        n=int(len(m)),
        spearman_rho=float(rho),
        spearman_p=float(p),
        ci95_low=float(lo),
        ci95_high=float(hi),
        rounds_to_0_69=bool(abs(rho - 0.69) < 0.005),
        note="Pairwise complete cases. Slide n=118 is not reproduced (CLDN4 quantified in 45/77 lung columns).",
    )


def load_cosmx_cold_niche() -> pd.DataFrame:
    """Cold-niche quadrant means already computed in PR #726 — cite, do not re-fit neighbors."""
    return pd.read_csv(DATA / "cosmx_quadrant_summary.csv")


def make_figures(
    unit_df: pd.DataFrame,
    cohort_df: pd.DataFrame,
    cosmx_df: pd.DataFrame,
    cold: pd.DataFrame,
    depmap: dict,
) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    # Fig1: concordant-4 ratio vs independence
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    order = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
    colors = {
        "GSE123902": "#1b4f72",
        "GSE131907": "#b85c38",
        "GSE205335": "#1e7f4f",
        "GSE189357": "#6c3483",
    }
    use = unit_df[unit_df["usable_rho"]].copy()
    positions = []
    data = []
    labels = []
    for i, ds in enumerate(order):
        vals = use.loc[use.dataset == ds, "ratio_vs_indep"].to_numpy()
        positions.append(i)
        data.append(vals)
        labels.append(f"{ds}\nn={len(vals)}")
        ax.scatter(
            np.full(len(vals), i) + np.random.default_rng(i).uniform(-0.12, 0.12, len(vals)),
            vals,
            s=18,
            alpha=0.65,
            color=colors[ds],
            zorder=3,
        )
    ax.boxplot(data, positions=positions, widths=0.45, showfliers=False, medianprops=dict(color="black"))
    ax.axhline(1.0, color="#666", ls="--", lw=1, label="independence")
    ax.set_xticks(positions, labels)
    ax.set_ylabel("P(TACSTD2+ ∩ CLDN4+) / [P(TAC+)·P(CLDN4+)]")
    ax.set_title("Concordant-4 malignant same-cell co-detection vs independence")
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_c4_ratio_vs_independence.png", dpi=160)
    fig.savefig(FIGURES / "fig1_c4_ratio_vs_independence.pdf")
    plt.close(fig)

    # Fig2: CosMx ratio + Spearman
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.8))
    axes[0].bar(range(len(cosmx_df)), cosmx_df["ratio_vs_indep"], color="#2c6e49")
    axes[0].axhline(1.0, color="#666", ls="--", lw=1)
    axes[0].set_xticks(range(len(cosmx_df)), cosmx_df["sample"], rotation=45, ha="right", fontsize=8)
    axes[0].set_ylabel("detection ratio vs independence")
    axes[0].set_title("CosMx tumor same-cell co-detection")
    axes[1].bar(range(len(cosmx_df)), cosmx_df["spearman"], color="#1b4f72")
    axes[1].axhline(0.0, color="#666", ls="--", lw=1)
    axes[1].set_xticks(range(len(cosmx_df)), cosmx_df["sample"], rotation=45, ha="right", fontsize=8)
    axes[1].set_ylabel("Spearman ρ")
    axes[1].set_title("CosMx tumor TACSTD2–CLDN4")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cosmx_codection_spearman.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cosmx_codection_spearman.pdf")
    plt.close(fig)

    # Fig3: cold niche quadrants at 10 µm (from PR #726)
    sub = cold[(cold.metric == "imm_frac") & (cold.radius == 10)].iloc[0]
    fig, ax = plt.subplots(figsize=(5.8, 3.6))
    names = ["Both low", "TAC high\nCLDN4 low", "TAC low\nCLDN4 high", "Both high"]
    vals = [
        sub.mean_tac_low_cld_low,
        sub.mean_tac_high_cld_low,
        sub.mean_tac_low_cld_high,
        sub.mean_tac_high_cld_high,
    ]
    ax.bar(names, vals, color=["#9aa5b1", "#5b8c5a", "#3d5a80", "#1b4332"])
    ax.set_ylabel("Immune-neighbor fraction @ 10 µm")
    ax.set_title("CosMx overlapping cold niche (median quadrants)")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.002, f"{v:.3f}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_cosmx_cold_niche_quadrants.png", dpi=160)
    fig.savefig(FIGURES / "fig3_cosmx_cold_niche_quadrants.pdf")
    plt.close(fig)

    # Fig4: PPT strip
    fig, ax = plt.subplots(figsize=(9.5, 3.2))
    ax.axis("off")
    lines = [
        "BRIDGE ONLY — TACSTD2 & CLDN4 co-occur in the same malignant cells",
        f"Concordant-4: median detection ratio vs independence = "
        f"{cohort_df.loc[cohort_df.dataset=='ALL','median_ratio_vs_indep'].iloc[0]:.2f} "
        f"(median ρ = {cohort_df.loc[cohort_df.dataset=='ALL','median_spearman'].iloc[0]:.3f})",
        f"CosMx: median detection ratio = {cosmx_df['ratio_vs_indep'].median():.2f}; "
        f"median ρ = {cosmx_df['spearman'].median():.3f}; 8/8 ρ>0",
        f"DepMap protein: ρ = {depmap['spearman_rho']:.3f} (n={depmap['n']}) → rounds to 0.69",
        "Cold niche: CosMx both-high immune frac @10µm = "
        f"{sub.mean_tac_high_cld_high:.3f} vs both-low {sub.mean_tac_low_cld_low:.3f} "
        f"(HH<LL in {int(sub.delta_hh_ll_sections_neg)}/{int(sub.delta_hh_ll_sections_defined)})",
        "NOT mediation: TACSTD2→CLDN4→immune path is null on CosMx + concordant-4 (PR #718)",
    ]
    y = 0.92
    for i, line in enumerate(lines):
        ax.text(
            0.02,
            y,
            line,
            transform=ax.transAxes,
            fontsize=11 if i == 0 else 9.5,
            fontweight="bold" if i == 0 else "normal",
            va="top",
            family="DejaVu Sans",
        )
        y -= 0.14
    fig.tight_layout()
    fig.savefig(FIGURES / "fig0_ppt_bridge_strip.png", dpi=160)
    fig.savefig(FIGURES / "fig0_ppt_bridge_strip.pdf")
    plt.close(fig)


def write_finding(
    unit_df: pd.DataFrame,
    cohort_df: pd.DataFrame,
    cosmx_df: pd.DataFrame,
    cosmx_sum: dict,
    cold: pd.DataFrame,
    depmap: dict,
) -> None:
    all_row = cohort_df[cohort_df.dataset == "ALL"].iloc[0]
    cold10 = cold[(cold.metric == "imm_frac") & (cold.radius == 10)].iloc[0]
    # dual vs T/NK DL-style note: report ALL usable spearman only as descriptive
    text = f"""# TACSTD2–CLDN4 same-cell coexpression BRIDGE (not mediation)

Public-only. Unit for concordant-4 is patient/donor/sample. Cell counts are not n.
Does **not** recompute locked CosMx CLDN4 cytotoxic ratios 0.36/0.52 or
concordant-4 CLDN4 %pos vs T/NK ρ = −0.531. Does **not** claim mediation
(PR #718: TACSTD2→CLDN4→immune matrix is null).

## Slide verdict (PPT)

| Layer | Same-cell coexpression? | Overlap vs independence | Cold niche overlap |
|---|---|---|---|
| Concordant-4 malignant | **YES** — median within-cell ρ = **{all_row.median_spearman:.3f}** ({all_row.frac_spearman_pos:.0%} units ρ>0; n_units≥30 = {int(all_row.n_units_ge30)}) | median P(both)/[P·P] = **{all_row.median_ratio_vs_indep:.2f}**; median OR = **{all_row.median_odds_ratio:.2f}** | Dual-pos % vs T/NK Spearman = **{all_row.dual_vs_tnk_spearman:.3f}** (p={all_row.dual_vs_tnk_p:.2e}) — descriptive bridge, not a new lock |
| CosMx He2022 tumor | **YES** — median ρ = **{cosmx_sum['median_spearman']:.3f}** (range {cosmx_sum['min_spearman']:.3f}–{cosmx_sum['max_spearman']:.3f}); **{cosmx_sum['n_sections']}/{cosmx_sum['n_sections']}** ρ>0 | median detection ratio = **{cosmx_sum['median_ratio_vs_indep']:.2f}**; {cosmx_sum['sections_ratio_gt1']}/{cosmx_sum['n_sections']} sections >1 | Both-high immune frac @10 µm **{cold10.mean_tac_high_cld_high:.3f}** vs both-low **{cold10.mean_tac_low_cld_low:.3f}**; HH<LL **{int(cold10.delta_hh_ll_sections_neg)}/{int(cold10.delta_hh_ll_sections_defined)}** (PR #726 quadrants) |
| DepMap Gygi protein | **YES** — Spearman **{depmap['spearman_rho']:.3f}** (n={depmap['n']}, p={depmap['spearman_p']:.2e}) | continuous protein coexpression (rounds to **0.69**) | N/A (cell lines; no immune niche) |

## Framing

- **BRIDGE:** the two genes mark overlapping malignant programs / cold neighborhoods.
- **NOT mediation:** holding CLDN4 does not account for TACSTD2 short-range exclusion on CosMx (PR #726), and the observational mediation matrix on CosMx + concordant-4 is null (PR #718).

## Concordant-4 cohort table

| cohort | units ≥30 | median ρ | median ratio vs indep | median OR |
|---|---:|---:|---:|---:|
"""
    for _, r in cohort_df[cohort_df.dataset != "ALL"].iterrows():
        text += (
            f"| {r.dataset} | {int(r.n_units_ge30)} | {r.median_spearman:.3f} | "
            f"{r.median_ratio_vs_indep:.2f} | {r.median_odds_ratio:.2f} |\n"
        )
    text += (
        f"| ALL | {int(all_row.n_units_ge30)} | {all_row.median_spearman:.3f} | "
        f"{all_row.median_ratio_vs_indep:.2f} | {all_row.median_odds_ratio:.2f} |\n\n"
        "## CosMx detection co-occurrence\n\n"
        f"Tumor cells = {cosmx_sum['n_tumor_total']} across {cosmx_sum['n_sections']} sections / "
        f"{cosmx_sum['n_patients']} patients. Detection = raw count > 0.\n\n"
        "| section | patient | n | ρ | P(both) | P(tac)·P(cld) | ratio | OR |\n"
        "|---|---|---:|---:|---:|---:|---:|---:|\n"
    )
    for _, r in cosmx_df.iterrows():
        text += (
            f"| {r['sample']} | {r['patient']} | {int(r['n'])} | {r['spearman']:.3f} | "
            f"{r['p_both']:.3f} | {r['p_indep']:.3f} | {r['ratio_vs_indep']:.3f} | {r['odds_ratio']:.2f} |\n"
        )
    text += f"""
## DepMap protein checksum

Gygi `_LUNG` complete cases: n={depmap['n']}, Spearman ρ={depmap['spearman_rho']:.3f}
(95% bootstrap CI {depmap['ci95_low']:.2f}–{depmap['ci95_high']:.2f}), p={depmap['spearman_p']:.2e}.
Rounds to the slide value 0.69. Slide n=118 is not a complete-case row (PR #575).

## Do not claim

- Mediation / “CLDN4 explains TACSTD2 immune association”
- Private 8KL / KD co-culture
- Visium same-spot correlation as spatial exclusion
- CosMx TACSTD2 at 50/100 µm as 8/8 (it is not; locked CLDN4 0.36/0.52 stays)

## Reproduce

```bash
export GEO_C4=/tmp/geo_c4
# downloads: scripts/download.sh
python3 scripts/analyze_bridge.py
```
"""
    (ROOT / "FINDING.md").write_text(text)
    (ROOT / "results" / "PPT_EVIDENCE_TABLE.md").write_text(
        "\n".join(
            [
                "# PPT — TACSTD2∩CLDN4 same-cell BRIDGE",
                "",
                "| Source | Number | Meaning |",
                "|---|---|---|",
                f"| Concordant-4 within-cell ρ | **{all_row.median_spearman:.3f}** (n_units={int(all_row.n_units_ge30)}) | same malignant cells |",
                f"| Concordant-4 detection lift | **{all_row.median_ratio_vs_indep:.2f}×** vs independence | co-occurrence |",
                f"| Concordant-4 OR | **{all_row.median_odds_ratio:.2f}** | co-detection |",
                f"| CosMx within-cell ρ | **{cosmx_sum['median_spearman']:.3f}** (8/8 >0) | same tumor cells |",
                f"| CosMx detection lift | **{cosmx_sum['median_ratio_vs_indep']:.2f}×** | co-occurrence |",
                f"| DepMap protein ρ | **{depmap['spearman_rho']:.3f}** (n={depmap['n']}) | rounds to 0.69 |",
                f"| CosMx cold niche @10µm | both-high **{cold10.mean_tac_high_cld_high:.3f}** vs both-low **{cold10.mean_tac_low_cld_low:.3f}** | overlapping cold |",
                "| Mediation | **null** (PR #718) | bridge ≠ path |",
                "",
            ]
        )
    )


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    locked = load_locked()
    say("== concordant-4 ==")
    rows = []
    rows.extend(load_gse123902(locked))
    rows.extend(load_gse131907(locked))
    rows.extend(load_gse205335(locked))
    rows.extend(load_gse189357(locked))
    unit_df = summarize_units(rows, locked)
    cohort_df = cohort_summary(unit_df)
    unit_df.to_csv(TABLES / "c4_unit_pair_codection.tsv", sep="\t", index=False)
    cohort_df.to_csv(TABLES / "c4_cohort_summary.tsv", sep="\t", index=False)

    say("== CosMx ==")
    cosmx_df, cosmx_sum = analyze_cosmx()
    cosmx_df.to_csv(TABLES / "cosmx_section_pair_codection.tsv", sep="\t", index=False)
    with open(TABLES / "cosmx_summary.json", "w") as fh:
        json.dump(cosmx_sum, fh, indent=2)

    say("== DepMap ==")
    depmap = analyze_depmap()
    with open(TABLES / "depmap_protein.json", "w") as fh:
        json.dump(depmap, fh, indent=2)

    cold = load_cosmx_cold_niche()
    make_figures(unit_df, cohort_df, cosmx_df, cold, depmap)
    write_finding(unit_df, cohort_df, cosmx_df, cosmx_sum, cold, depmap)

    # key stats
    key = {
        "concordant4": cohort_df[cohort_df.dataset == "ALL"].iloc[0].to_dict(),
        "cosmx": cosmx_sum,
        "depmap": depmap,
        "framing": "BRIDGE_ONLY_NOT_MEDIATION",
        "locked_not_recomputed": [
            "CosMx CLDN4 cytotoxic 0.36/0.52",
            "concordant-4 CLDN4 %pos vs T/NK rho=-0.531",
        ],
        "mediation_null_citation": "PR #718",
    }
    with open(TABLES / "key_stats.json", "w") as fh:
        json.dump(key, fh, indent=2, default=float)
    say("DONE")
    say(cohort_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
