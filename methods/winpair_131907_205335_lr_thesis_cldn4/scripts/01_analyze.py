#!/usr/bin/env python3
"""Thesis-aligned CLDN4-only ligand families on GSE131907 + GSE205335.

Both pre-specified arms:
  1) barrier / inhibitory pairs — expected MORE from CLDN4-high malignant → T/NK
  2) IFN / T-recruit / MHC-I pairs — expected MORE from CLDN4-low (KD-like) malignant → T/NK

Patient is the unit. Honest n. No dual-high. No GSE148071.
PR #320 T/NK ρ is given and is not re-audited.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]

TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
GIVEN_Q4Q1 = {"n_compared": 23, "n_q1": 12, "n_q4": 11, "r": -0.705, "pr": 320}
EXPR_PROP = 0.10
MIN_MALIG_BIN = 10
MIN_TNK = 20
STATE_GENES = {"CLDN4", "TACSTD2"}
ALIASES = {"PVRL2": "NECTIN2"}


def log(msg: str) -> None:
    print(msg, flush=True)


def fmt_p(value) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_delta(value) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:+.3f}"


def parse_series_matrix(path: Path) -> pd.DataFrame:
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
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


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


def stream_umi(path: Path, keep: set[str]) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, int]:
    log(f"[stream] {path} keep={sorted(keep)}")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            gene = ALIASES.get(gene, gene)
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"column mismatch for {gene}: {vals.size} != {n}")
            totals += vals
            if gene in keep and gene not in store:
                store[gene] = vals
            if n_genes % 5000 == 0:
                log(f"[stream] genes_seen={n_genes} kept={len(store)}")
    log(f"[stream] cells={n} genes_in_file={n_genes} genes_kept={len(store)}")
    return cells, store, totals, n_genes


def load_rds_genes(path: Path, wanted: set[str]) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, list[str]]:
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            log(f"decompress {path.name}")
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        log("read RDS")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    log(f"build CSC {tuple(obj.Dim)}")
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {ALIASES.get(g, g): i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
    log(f"extracted {len(extracted)} / {len(wanted)} genes")
    return extracted, library_umi, barcodes, genes.tolist()


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        cp = np.where(total > 0, umi / total * 1e4, 0.0)
    return np.log1p(cp).astype(np.float32)


def group_gene_stats(logx: dict[str, np.ndarray], mask: np.ndarray) -> tuple[dict[str, float], dict[str, float], int]:
    n = int(mask.sum())
    means, fracs = {}, {}
    if n == 0:
        return means, fracs, 0
    for g, arr in logx.items():
        v = arr[mask]
        means[g] = float(np.mean(v))
        fracs[g] = float(np.mean(v > 0))
    return means, fracs, n


def partner_from_stats(units: list[str], means: dict[str, float], fracs: dict[str, float]) -> tuple[float, float]:
    m, f = [], []
    for g in units:
        if g not in means:
            return np.nan, np.nan
        m.append(means[g])
        f.append(fracs[g])
    return float(np.min(m)), float(np.min(f))


def score_pairs(pairs: pd.DataFrame, logx: dict[str, np.ndarray], sender: np.ndarray, receiver: np.ndarray) -> pd.DataFrame:
    s_mean, s_frac, n_s = group_gene_stats(logx, sender)
    r_mean, r_frac, n_r = group_gene_stats(logx, receiver)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac = partner_from_stats(lig_u, s_mean, s_frac)
        rec_m, rec_f = partner_from_stats(rec_u, r_mean, r_frac)
        if not np.isfinite(l_mean) or not np.isfinite(rec_m):
            rows.append(
                {
                    "pair": rec.pair,
                    "ligand": rec.ligand,
                    "receptor": rec.receptor,
                    "family": rec.family,
                    "axis": rec.axis,
                    "thesis_direction": rec.thesis_direction,
                    "n_sender": n_s,
                    "n_receiver": n_r,
                    "ligand_mean": np.nan,
                    "receptor_mean": np.nan,
                    "ligand_frac": np.nan,
                    "receptor_frac": np.nan,
                    "cpdb_mean_score": np.nan,
                    "pass_expr_prop": False,
                    "gene_missing": True,
                }
            )
            continue
        pass_prop = (l_frac >= EXPR_PROP) and (rec_f >= EXPR_PROP)
        rows.append(
            {
                "pair": rec.pair,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "family": rec.family,
                "axis": rec.axis,
                "thesis_direction": rec.thesis_direction,
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": rec_m,
                "ligand_frac": l_frac,
                "receptor_frac": rec_f,
                "cpdb_mean_score": 0.5 * (l_mean + rec_m),
                "pass_expr_prop": bool(pass_prop),
                "gene_missing": False,
            }
        )
    return pd.DataFrame(rows)


def wilcoxon_safe(a, b, alternative: str = "two-sided") -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3 or np.allclose(a, b):
        return np.nan
    try:
        return float(stats.wilcoxon(a, b, zero_method="wilcox", alternative=alternative).pvalue)
    except ValueError:
        return np.nan


def load_gse131907(workdir: Path, keep: set[str], locked_samples: set[str]) -> dict:
    ann = pd.read_csv(workdir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str)
    meta = parse_series_matrix(workdir / "GSE131907_series_matrix.txt.gz")
    sample_meta = meta.rename(columns={"title": "Sample"})
    keep_cols = [c for c in ["Sample", "geo_accession", "patient_id", "tumor_stage"] if c in sample_meta.columns]
    sample_meta = sample_meta[keep_cols].drop_duplicates("Sample")
    cells, expr, totals, n_genes = stream_umi(workdir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", keep)
    per = ann.set_index("Index").reindex(cells)
    if per["Sample"].isna().any():
        raise SystemExit("GSE131907 matrix cell IDs do not align with annotation Index")
    sample = per["Sample"].to_numpy()
    origin = per["Sample_Origin"].to_numpy()
    ctype = per["Cell_type"].fillna("").to_numpy()
    csub = per["Cell_subtype"].fillna("").to_numpy()
    patient_map = sample_meta.set_index("Sample")["patient_id"].to_dict()
    patient = np.array([patient_map.get(s, s) for s in sample], dtype=object)
    is_tumor = np.isin(origin, list(TUMOR_ORIGINS))
    is_locked_sample = np.isin(sample, list(locked_samples))
    is_malig = is_tumor & is_locked_sample & np.isin(csub, list(MALIGNANT_SUBTYPES))
    is_tnk = is_locked_sample & np.isin(ctype, ["T lymphocytes", "NK cells"])
    logx = {g: log1p_cp10k(arr, totals) for g, arr in expr.items()}
    if "CLDN4" not in logx:
        raise SystemExit("CLDN4 missing from GSE131907")
    return {
        "cohort": "GSE131907",
        "cells": cells,
        "logx": logx,
        "sample": sample,
        "patient": patient,
        "is_malig": is_malig,
        "is_tnk": is_tnk,
        "n_genes": n_genes,
        "n_cells_matrix": int(len(cells)),
        "n_locked_samples": int(len(locked_samples)),
        "n_locked_patients_geo": int(pd.Series([patient_map.get(s, s) for s in locked_samples]).nunique()),
        "genes_kept": sorted(logx),
    }


def load_gse205335(identities: Path, soft: Path, matrix: Path, keep: set[str], locked_patients: set[str]) -> dict:
    ident = pd.read_csv(identities, sep="\t")
    if ident["barcode"].duplicated().any():
        raise ValueError("GSE205335 Cell-identity barcodes are not unique")
    metadata = parse_geo_soft(soft)
    extracted, library_umi, barcodes, all_genes = load_rds_genes(matrix, keep)
    indexed = ident.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise ValueError(f"GSE205335 matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra")
    cells = indexed.loc[barcodes].reset_index()
    cells = cells.merge(
        metadata[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise ValueError("GSE205335 identity-table samples did not match GEO metadata")
    is_malig = cells["lineage.sub"].eq("Malignant cells").to_numpy() & cells["patient"].isin(locked_patients).to_numpy()
    is_tnk = cells["lineage.total"].eq("T/NK cells").to_numpy() & cells["patient"].isin(locked_patients).to_numpy()
    logx = {g: log1p_cp10k(extracted[g], library_umi) for g in extracted}
    if "CLDN4" not in logx:
        raise SystemExit("CLDN4 missing from GSE205335")
    return {
        "cohort": "GSE205335",
        "cells": barcodes,
        "logx": logx,
        "sample": cells["orig.ident"].to_numpy(),
        "patient": cells["patient"].to_numpy(),
        "is_malig": is_malig,
        "is_tnk": is_tnk,
        "n_genes": int(len(all_genes)),
        "n_cells_matrix": int(len(barcodes)),
        "n_locked_patients": int(len(locked_patients)),
        "genes_kept": sorted(logx),
    }


def split_cldn4(obj: dict) -> None:
    cldn4 = obj["logx"]["CLDN4"]
    mal = obj["is_malig"]
    thr = float(np.median(cldn4[mal])) if mal.any() else float("nan")
    obj["cldn4"] = cldn4
    obj["cldn4_threshold"] = thr
    obj["cld_hi"] = mal & (cldn4 >= thr)
    obj["cld_lo"] = mal & (cldn4 < thr)
    log(
        f"[{obj['cohort']}] CLDN4 median={thr:.4f} "
        f"mal={int(mal.sum())} high={int(obj['cld_hi'].sum())} low={int(obj['cld_lo'].sum())} "
        f"tnk={int(obj['is_tnk'].sum())}"
    )


def patient_n_table(obj: dict) -> pd.DataFrame:
    rows = []
    for pid in pd.unique(obj["patient"]):
        m = obj["patient"] == pid
        if not (m & (obj["is_malig"] | obj["is_tnk"])).any():
            continue
        n_hi = int((m & obj["cld_hi"]).sum())
        n_lo = int((m & obj["cld_lo"]).sum())
        n_tnk = int((m & obj["is_tnk"]).sum())
        n_mal = int((m & obj["is_malig"]).sum())
        samples = sorted(set(obj["sample"][m & (obj["is_malig"] | obj["is_tnk"])].astype(str)))
        rows.append(
            {
                "cohort": obj["cohort"],
                "patient_id": pid,
                "n_samples": len(samples),
                "samples": ",".join(samples),
                "n_malig": n_mal,
                "n_cldn4_high": n_hi,
                "n_cldn4_low": n_lo,
                "n_T_NK": n_tnk,
                "mean_cldn4_malig": float(np.mean(obj["cldn4"][m & obj["is_malig"]])) if n_mal else np.nan,
                "eligible_paired": n_hi >= MIN_MALIG_BIN and n_lo >= MIN_MALIG_BIN and n_tnk >= MIN_TNK and n_mal > 0,
            }
        )
    return pd.DataFrame(rows)


def score_patients(obj: dict, pairs: pd.DataFrame, eligible: set[str]) -> pd.DataFrame:
    chunks = []
    for pid in sorted(eligible):
        m = obj["patient"] == pid
        h, l, t = m & obj["cld_hi"], m & obj["cld_lo"], m & obj["is_tnk"]
        hdf = score_pairs(pairs, obj["logx"], h, t)
        ldf = score_pairs(pairs, obj["logx"], l, t)
        key = ["pair", "ligand", "receptor", "family", "axis", "thesis_direction"]
        d = hdf.merge(ldf, on=key, suffixes=("_high", "_low"))
        d["delta_high_minus_low"] = d["cpdb_mean_score_high"] - d["cpdb_mean_score_low"]
        d["pass_either"] = d["pass_expr_prop_high"] | d["pass_expr_prop_low"]
        d["pass_both"] = d["pass_expr_prop_high"] & d["pass_expr_prop_low"]
        d["gene_missing"] = d["gene_missing_high"] | d["gene_missing_low"]
        d["cohort"] = obj["cohort"]
        d["patient_id"] = pid
        d["unit_id"] = f"{obj['cohort']}:{pid}"
        d["n_high"] = int(h.sum())
        d["n_low"] = int(l.sum())
        d["n_tnk"] = int(t.sum())
        chunks.append(d)
        log(f"  scored {obj['cohort']} {pid} high={int(h.sum())} low={int(l.sum())} tnk={int(t.sum())}")
    return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()


def onesided_alt(thesis_direction: str) -> str:
    return "greater" if thesis_direction == "high>low" else "less"


def agrees_thesis(delta: float, thesis_direction: str) -> bool:
    if not np.isfinite(delta) or delta == 0:
        return False
    return (delta > 0) if thesis_direction == "high>low" else (delta < 0)


def rank_family(raw: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rec in pairs.itertuples(index=False):
        sub = raw[raw["pair"] == rec.pair].copy()
        finite = sub[np.isfinite(sub["delta_high_minus_low"])]
        n = int(finite["unit_id"].nunique())
        n_131 = int(finite.loc[finite["cohort"] == "GSE131907", "unit_id"].nunique())
        n_205 = int(finite.loc[finite["cohort"] == "GSE205335", "unit_id"].nunique())
        n_nonzero = int((finite["delta_high_minus_low"].abs() > 0).sum())
        n_pass = int(finite.loc[finite["pass_either"], "unit_id"].nunique())
        p_two = wilcoxon_safe(finite["cpdb_mean_score_high"], finite["cpdb_mean_score_low"], "two-sided")
        p_one = wilcoxon_safe(
            finite["cpdb_mean_score_high"],
            finite["cpdb_mean_score_low"],
            onesided_alt(rec.thesis_direction),
        )
        med = float(np.median(finite["delta_high_minus_low"])) if n else np.nan
        rows.append(
            {
                "family": rec.family,
                "pair": rec.pair,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "axis": rec.axis,
                "thesis_direction": rec.thesis_direction,
                "n": n,
                "n_GSE131907": n_131,
                "n_GSE205335": n_205,
                "n_nonzero_delta": n_nonzero,
                "n_pass_expr_prop": n_pass,
                "delta": med,
                "mean_delta": float(np.mean(finite["delta_high_minus_low"])) if n else np.nan,
                "mean_score_high": float(finite["cpdb_mean_score_high"].mean()) if n else np.nan,
                "mean_score_low": float(finite["cpdb_mean_score_low"].mean()) if n else np.nan,
                "median_ligand_frac_high": float(finite["ligand_frac_high"].median()) if n else np.nan,
                "median_ligand_frac_low": float(finite["ligand_frac_low"].median()) if n else np.nan,
                "median_receptor_frac": float(finite["receptor_frac_high"].median()) if n else np.nan,
                "p": p_two,
                "p_onesided_thesis": p_one,
                "agrees_thesis": agrees_thesis(med, rec.thesis_direction) if n else False,
                "sparse_sender": bool(n and float(finite["ligand_frac_high"].median()) < EXPR_PROP and float(finite["ligand_frac_low"].median()) < EXPR_PROP),
                "gene_missing": bool(sub["gene_missing"].all()) if len(sub) else True,
                "note": rec.note,
            }
        )
    out = pd.DataFrame(rows)
    if len(out) and out["p"].notna().any():
        mask = out["p"].notna()
        out.loc[mask, "padj_family"] = np.nan
        for fam, idx in out.groupby("family").groups.items():
            submask = out.index.isin(idx) & mask
            if submask.any():
                out.loc[submask, "padj_family"] = multipletests(out.loc[submask, "p"], method="fdr_bh")[1]
    else:
        out["padj_family"] = np.nan
    return out


def family_composite(raw: pd.DataFrame, family: str) -> dict:
    sub = raw[raw["family"] == family].copy()
    if sub.empty:
        return {"family": family, "n": 0, "delta": np.nan, "p": np.nan}
    wide = (
        sub.groupby("unit_id")
        .agg(
            score_high=("cpdb_mean_score_high", "mean"),
            score_low=("cpdb_mean_score_low", "mean"),
            cohort=("cohort", "first"),
        )
        .reset_index()
    )
    wide["delta"] = wide["score_high"] - wide["score_low"]
    thesis = "high>low" if family == "barrier_inhibitory" else "low>high"
    p_two = wilcoxon_safe(wide["score_high"], wide["score_low"], "two-sided")
    p_one = wilcoxon_safe(wide["score_high"], wide["score_low"], onesided_alt(thesis))
    med = float(np.median(wide["delta"]))
    return {
        "family": family,
        "thesis_direction": thesis,
        "n": int(len(wide)),
        "n_GSE131907": int((wide["cohort"] == "GSE131907").sum()),
        "n_GSE205335": int((wide["cohort"] == "GSE205335").sum()),
        "delta": med,
        "mean_delta": float(np.mean(wide["delta"])),
        "p": p_two,
        "p_onesided_thesis": p_one,
        "agrees_thesis": agrees_thesis(med, thesis),
        "composite": "mean of pre-specified pair scores per patient",
    }


def md_table(df: pd.DataFrame, cols: list[str], fmt: dict[str, str] | None = None) -> list[str]:
    if df is None or df.empty:
        return ["_(empty)_", ""]
    fmt = fmt or {}
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for rec in df.itertuples(index=False):
        row = []
        for c in cols:
            v = getattr(rec, c)
            if c in fmt and pd.notna(v):
                row.append(fmt[c].format(v))
            elif isinstance(v, float) and pd.notna(v):
                row.append(f"{v:.3g}")
            else:
                row.append(str(v))
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    return lines


def plot_n_cells(ntab: pd.DataFrame, path: Path) -> None:
    show = ntab.sort_values(["cohort", "patient_id"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10.5, 4.0))
    x = np.arange(len(show))
    ax.bar(x - 0.2, show["n_malig"], width=0.4, label="malignant", color="#4C72B0")
    ax.bar(x + 0.2, show["n_T_NK"], width=0.4, label="T/NK", color="#DD8452")
    ax.set_xticks(x)
    labels = [f"{c}:{p}" for c, p in zip(show["cohort"], show["patient_id"])]
    ax.set_xticklabels(labels, rotation=80, ha="right", fontsize=6)
    ax.set_ylabel("cells")
    ax.legend(frameon=False)
    ax.set_title("Locked slice cells (author malignant vs same-patient T/NK)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_family_bars(fam: pd.DataFrame, path: Path, title: str, expected_sign: int) -> None:
    show = fam.copy()
    fig, ax = plt.subplots(figsize=(8.6, max(3.2, 0.42 * len(show) + 1.4)))
    y = np.arange(len(show))
    colors = []
    for rec in show.itertuples(index=False):
        if rec.agrees_thesis:
            colors.append("#b2182b" if expected_sign > 0 else "#2166ac")
        elif np.isfinite(rec.delta) and rec.delta != 0:
            colors.append("#999999")
        else:
            colors.append("#dddddd")
    ax.barh(y, show["delta"].fillna(0), color=colors, edgecolor="none")
    ax.axvline(0, color="k", lw=0.8)
    labels = []
    for rec in show.itertuples(index=False):
        labels.append(f"{rec.pair}  n={int(rec.n)}  p={fmt_p(rec.p)}")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("median patient Δ (CLDN4-high − low)")
    ax.set_title(title)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_paired_boxes(raw: pd.DataFrame, pairs: list[str], path: Path, title: str) -> None:
    n = len(pairs)
    fig, axes = plt.subplots(2, int(np.ceil(n / 2)), figsize=(3.1 * int(np.ceil(n / 2)), 6.2), sharey=False)
    axes = np.atleast_1d(axes).ravel()
    for i, pair in enumerate(pairs):
        ax = axes[i]
        sub = raw[raw["pair"] == pair]
        if sub.empty:
            ax.set_title(pair, fontsize=8)
            ax.axis("off")
            continue
        vals = [sub["cpdb_mean_score_low"].to_numpy(), sub["cpdb_mean_score_high"].to_numpy()]
        ax.boxplot(vals, tick_labels=["low", "high"], widths=0.55, showfliers=False)
        for rec in sub.itertuples(index=False):
            ax.plot([1, 2], [rec.cpdb_mean_score_low, rec.cpdb_mean_score_high], color="0.7", lw=0.5, alpha=0.7)
        ax.set_title(f"{pair}\nn={sub['unit_id'].nunique()}", fontsize=8)
        ax.set_ylabel("CPDB-style score" if i % int(np.ceil(n / 2)) == 0 else "")
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_cohort_split(raw: pd.DataFrame, pairs: pd.DataFrame, path: Path) -> None:
    rows = []
    for rec in pairs.itertuples(index=False):
        sub = raw[raw["pair"] == rec.pair]
        for cohort, csub in sub.groupby("cohort"):
            if csub.empty:
                continue
            rows.append(
                {
                    "pair": rec.pair,
                    "family": rec.family,
                    "cohort": cohort,
                    "delta": float(np.median(csub["delta_high_minus_low"])),
                    "n": int(csub["unit_id"].nunique()),
                }
            )
    d = pd.DataFrame(rows)
    if d.empty:
        return
    order = list(pairs["pair"])
    fig, axes = plt.subplots(1, 2, figsize=(11.2, max(3.6, 0.32 * len(order) + 1.2)), sharex=True)
    for ax, fam, title in (
        (axes[0], "barrier_inhibitory", "Barrier / inhibitory (thesis: Δ > 0)"),
        (axes[1], "ifn_recruit", "IFN / T-recruit / MHC-I (thesis: Δ < 0)"),
    ):
        fam_pairs = [p for p in order if p in set(d.loc[d["family"] == fam, "pair"])]
        y = np.arange(len(fam_pairs))
        for i, cohort in enumerate(["GSE131907", "GSE205335"]):
            vals = []
            for pair in fam_pairs:
                hit = d[(d["pair"] == pair) & (d["cohort"] == cohort)]
                vals.append(float(hit["delta"].iloc[0]) if len(hit) else np.nan)
            ax.barh(y + (i - 0.5) * 0.36, vals, height=0.34, label=cohort, alpha=0.85)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(fam_pairs, fontsize=7)
        ax.set_xlabel("median Δ (high − low)")
        ax.set_title(title, fontsize=9)
        ax.legend(frameon=False, fontsize=7)
        ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_two_arms(barrier: pd.DataFrame, recruit: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.4), sharex=True)
    for ax, fam, title, exp in (
        (axes[0], barrier, "Barrier / inhibitory\nthesis: more from CLDN4-high", +1),
        (axes[1], recruit, "IFN / T-recruit / MHC-I\nthesis: more from CLDN4-low", -1),
    ):
        y = np.arange(len(fam))
        colors = ["#b2182b" if (exp > 0 and a) or (exp < 0 and a) else "#7f7f7f" for a in fam["agrees_thesis"]]
        ax.barh(y, fam["delta"].fillna(0), color=colors, edgecolor="none")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{p} (n={int(n)})" for p, n in zip(fam["pair"], fam["n"])], fontsize=7)
        ax.set_xlabel("median Δ (high − low)")
        ax.set_title(title, fontsize=10)
        ax.invert_yaxis()
    fig.suptitle("Both pre-specified families (patient-paired Wilcoxon)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(
    outdir: Path,
    summary: dict,
    ntab: pd.DataFrame,
    barrier: pd.DataFrame,
    recruit: pd.DataFrame,
    composites: pd.DataFrame,
) -> None:
    n_paired = int(summary["n_patients_paired"])
    n_131 = int(summary["n_patients_paired_GSE131907"])
    n_205 = int(summary["n_patients_paired_GSE205335"])
    b_hit = int(barrier["agrees_thesis"].sum())
    r_hit = int(recruit["agrees_thesis"].sum())
    b_sig = int(((barrier["p"] < 0.05) & barrier["agrees_thesis"]).sum())
    r_sig = int(((recruit["p"] < 0.05) & recruit["agrees_thesis"]).sum())
    b_comp = composites[composites["family"] == "barrier_inhibitory"].iloc[0]
    r_comp = composites[composites["family"] == "ifn_recruit"].iloc[0]

    verdict = (
        f"On {n_paired} paired patients ({n_131} GSE131907 + {n_205} GSE205335), "
        f"the barrier/inhibitory family is **{b_hit}/{len(barrier)} pairs in the thesis direction** "
        f"(high > low; composite Δ={fmt_delta(b_comp['delta'])}, p={fmt_p(b_comp['p'])}; "
        f"{b_sig} pairs p<0.05 and Δ>0). "
        f"The IFN/T-recruit/MHC-I family is **{r_hit}/{len(recruit)} pairs in the thesis direction** "
        f"(low > high; composite Δ={fmt_delta(r_comp['delta'])}, p={fmt_p(r_comp['p'])}; "
        f"{r_sig} pairs p<0.05 and Δ<0). "
        "Barrier-up-in-high is a primary arm, not a recruit leftover. "
        "Scores are CellPhoneDB-style co-expression, not secretion or contact."
    )

    lines = [
        "# FINDING — winning-pair thesis LR (CLDN4-only, both families)",
        "",
        f"**Verdict:** {verdict}",
        "",
        "ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No GSE148071. "
        "Patient is the unit. p-values are descriptive.",
        "",
        "The between-patient abundance cut on this same slice is **given** and is not re-audited: "
        f"PR #{GIVEN_Q4Q1['pr']} Q4 vs Q1 author %pos vs T/NK "
        f"**n={GIVEN_Q4Q1['n_compared']} ({GIVEN_Q4Q1['n_q1']}/{GIVEN_Q4Q1['n_q4']}) r={GIVEN_Q4Q1['r']:+.3f}**. "
        "This folder tests a different contrast: outgoing CellPhoneDB-style scores from "
        "CLDN4-high vs CLDN4-low malignant cells to **same-patient** T/NK, "
        "on two **pre-specified** pair families.",
        "",
        "---",
        "",
        "## Thesis (pre-specified; both arms)",
        "",
        "1. **CLDN4-high** malignant cells send **more barrier / inhibitory** pairs to T/NK: "
        "F11R–ITGAL/ITGB2, NECTIN2–TIGIT, CDH1–ITGAE, LGALS9–CD45/CD44.",
        "2. **CLDN4-low (KD-like)** malignant cells send **more IFN / T-recruit / MHC-I** pairs: "
        "CXCL9/10/11–CXCR3, CCL5–CCR5, HLA-A/B/C–CD8.",
        "",
        "Both families are scored and tabulated. Barrier-up-in-high hits are not filed as "
        "\"recruit-up, skip.\"",
        "",
        "## Honest n",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        f"| GSE131907 cells in UMI matrix | {summary['n_cells_GSE131907']:,} | author barcodes |",
        f"| GSE205335 cells in UMI matrix | {summary['n_cells_GSE205335']:,} | author barcodes |",
        f"| GSE131907 locked samples (PR #320) | {summary['n_locked_samples_GSE131907']} | tumor-origin, n_malignant ≥20 |",
        f"| GSE131907 unique patients in locked samples | {summary['n_locked_patients_GSE131907']} | GEO patient_id; extract is sample-level |",
        f"| GSE205335 locked patients (PR #320) | {summary['n_locked_patients_GSE205335']} | tumor extract |",
        f"| Author malignant cells (locked slice) | {summary['n_malig']:,} | {summary['n_malig_GSE131907']:,} + {summary['n_malig_GSE205335']:,} |",
        f"| CLDN4-high / low malignant | {summary['n_cldn4_high']:,} / {summary['n_cldn4_low']:,} | cohort-specific global medians |",
        f"| T/NK cells (locked slice) | {summary['n_tnk']:,} | {summary['n_tnk_GSE131907']:,} + {summary['n_tnk_GSE205335']:,} |",
        f"| **Paired patients (Wilcoxon unit)** | **{n_paired}** | ≥10 high, ≥10 low, ≥20 T/NK |",
        f"|  … GSE131907 / GSE205335 | {n_131} / {n_205} | GSE131907 cells pooled by patient_id |",
        f"| Pre-specified pairs | {summary['n_pairs']} | 8 barrier/inhibitory + 10 IFN/recruit/MHC-I |",
        "",
        f"GSE131907 CLDN4 median log1p(CP10k) among locked malignant cells = {summary['cldn4_threshold_GSE131907']:.3f}. "
        f"GSE205335 median = {summary['cldn4_threshold_GSE205335']:.3f}. "
        "Thresholds are not shared across platforms.",
        "",
        "Per-patient counts: `results/n_cells_patients.tsv`.",
        "",
        "### Locked patients that failed the paired LR gate",
        "",
    ]
    dropped = ntab[~ntab["eligible_paired"]].copy()
    if dropped.empty:
        lines += ["None. Every locked patient had both CLDN4 bins and ≥20 T/NK.", ""]
    else:
        lines += [
            "| cohort | patient | n_high | n_low | n_T/NK | reason |",
            "|---|---|---:|---:|---:|---|",
        ]
        for rec in dropped.itertuples():
            reasons = []
            if rec.n_cldn4_high < MIN_MALIG_BIN:
                reasons.append("high<10")
            if rec.n_cldn4_low < MIN_MALIG_BIN:
                reasons.append("low<10")
            if rec.n_T_NK < MIN_TNK:
                reasons.append("T/NK<20")
            lines.append(
                f"| {rec.cohort} | {rec.patient_id} | {int(rec.n_cldn4_high)} | "
                f"{int(rec.n_cldn4_low)} | {int(rec.n_T_NK)} | {','.join(reasons) or 'other'} |"
            )
        lines.append("")

    lines += [
        "## Method",
        "",
        "Documented CellPhoneDB-style score on log1p(CP10k): partner expression = **min of subunit means**; "
        "pair score = **mean of the two partner means** (Efremova et al. 2020; Garcia-Alonso et al. 2022). "
        "Every pre-specified pair is scored on **every eligible patient** (no expr_prop drop from the primary n). "
        "`n_pass_expr_prop` is the number of patients in whom both partners are in ≥10% of cells in high or low; "
        "it is a sparsity flag, not a license to skip the pair. "
        "The test is a **paired Wilcoxon** of high vs low **per patient**. "
        "Δ = median(high − low). Positive Δ = stronger from CLDN4-high. "
        "Family FDR is Benjamini–Hochberg **within family**. "
        "Family composite = mean of that family's pair scores per patient, then the same Wilcoxon. "
        "This is **not** a CellChat communication probability and does **not** observe secretion or spatial contact. "
        "GSE131907 `PVRL2` is aliased to `NECTIN2`.",
        "",
        "## Primary table 1 — barrier / inhibitory (thesis: Δ > 0)",
        "",
        f"Family composite: n={int(b_comp['n'])}, Δ={fmt_delta(b_comp['delta'])}, p={fmt_p(b_comp['p'])} "
        f"(one-sided thesis p={fmt_p(b_comp['p_onesided_thesis'])}).",
        "",
        "Pairs in the thesis direction (Δ > 0) are **hits for this arm**, not a reason to skip to recruit.",
        "",
    ]
    lines += md_table(
        barrier,
        ["pair", "n", "n_GSE131907", "n_GSE205335", "delta", "p", "padj_family", "agrees_thesis", "n_pass_expr_prop", "sparse_sender"],
        {"delta": "{:+.3f}", "p": "{:.3g}", "padj_family": "{:.3g}"},
    )
    lines += [
        "Full table: `results/barrier_inhibitory_high.tsv`.",
        "",
        "## Primary table 2 — IFN / T-recruit / MHC-I (thesis: Δ < 0)",
        "",
        f"Family composite: n={int(r_comp['n'])}, Δ={fmt_delta(r_comp['delta'])}, p={fmt_p(r_comp['p'])} "
        f"(one-sided thesis p={fmt_p(r_comp['p_onesided_thesis'])}).",
        "",
        "A positive Δ here is a **miss for this arm** (higher from CLDN4-high, opposite the KD-like prediction). "
        "Sparse CXCR3 ligands keep their honest n; they are not dropped.",
        "",
    ]
    lines += md_table(
        recruit,
        ["pair", "n", "n_GSE131907", "n_GSE205335", "delta", "p", "padj_family", "agrees_thesis", "n_pass_expr_prop", "sparse_sender"],
        {"delta": "{:+.3f}", "p": "{:.3g}", "padj_family": "{:.3g}"},
    )
    lines += [
        "Full table: `results/ifn_recruit_low.tsv`.",
        "",
        "## Family composites",
        "",
    ]
    lines += md_table(
        composites,
        ["family", "thesis_direction", "n", "n_GSE131907", "n_GSE205335", "delta", "p", "p_onesided_thesis", "agrees_thesis"],
        {"delta": "{:+.3f}", "p": "{:.3g}", "p_onesided_thesis": "{:.3g}"},
    )
    lines += [
        "## What this is not",
        "",
        "- Not a re-audit of PR #320 Q4 vs Q1 n=23 r=−0.705.",
        "- Not dual-high TACSTD2×CLDN4.",
        "- Not GSE148071.",
        "- Not inferCNV/CopyKAT recomputed malignant IDs.",
        "- Not CellChat / LIANA permutation p-values. Those are not patient tests.",
        "- Not spatial proximity or protein secretion.",
        "- Cells are not the sample size. GSE131907 PR #320 extract is sample-level; this Wilcoxon uses unique patients.",
        "",
        "## Extra figures",
        "",
        "- `figures/n_cells_by_patient.png`",
        "- `figures/barrier_inhibitory_high.png`",
        "- `figures/ifn_recruit_low.png`",
        "- `figures/thesis_two_arms.png`",
        "- `figures/barrier_paired_boxes.png`",
        "- `figures/recruit_paired_boxes.png`",
        "- `figures/family_by_cohort.png`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/winpair_131907_205335_lr_thesis_cldn4/scripts/00_download.py",
        "python3 methods/winpair_131907_205335_lr_thesis_cldn4/scripts/01_analyze.py",
        "```",
        "",
    ]
    text = "\n".join(lines)
    (HERE / "FINDING.md").write_text(text)
    (outdir / "FINDING.md").write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse131907", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--gse205335", type=Path, default=Path("/tmp/gse205335"))
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    figdir = HERE / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "figures").mkdir(parents=True, exist_ok=True)

    pairs = pd.read_csv(HERE / "config" / "thesis_pairs.tsv", sep="\t")
    keep = set(STATE_GENES)
    for rec in pairs.itertuples(index=False):
        keep.update(str(rec.ligand).split("+"))
        keep.update(str(rec.receptor).split("+"))
    keep.update(ALIASES.keys())

    s131 = pd.read_csv(HERE / "data" / "GSE131907_samples.tsv", sep="\t")
    locked_samples = set(
        s131.loc[s131["origin"].isin(TUMOR_ORIGINS) & (s131["n_malignant"] >= 20), "sample"].astype(str)
    )
    p205 = pd.read_csv(HERE / "data" / "GSE205335_patients.tsv", sep="\t")
    locked_patients = set(p205["patient"].astype(str))
    log(f"[slice] GSE131907 locked samples={len(locked_samples)} GSE205335 locked patients={len(locked_patients)}")

    obj131 = load_gse131907(args.gse131907, keep, locked_samples)
    obj205 = load_gse205335(
        args.gse205335 / "GSE205335_Lung_IO_CellIdentity.txt.gz",
        args.gse205335 / "GSE205335_family.soft.gz",
        args.gse205335 / "GSE205335_Lung_IO_UMI_matrix.rds.gz",
        keep,
        locked_patients,
    )
    split_cldn4(obj131)
    split_cldn4(obj205)

    ntab = pd.concat([patient_n_table(obj131), patient_n_table(obj205)], ignore_index=True)
    ntab.to_csv(args.outdir / "n_cells_patients.tsv", sep="\t", index=False)
    eligible = ntab.loc[ntab["eligible_paired"], ["cohort", "patient_id"]]
    elig_131 = set(eligible.loc[eligible["cohort"] == "GSE131907", "patient_id"])
    elig_205 = set(eligible.loc[eligible["cohort"] == "GSE205335", "patient_id"])
    n_paired = int(len(elig_131) + len(elig_205))
    log(f"[n] paired patients={n_paired} (GSE131907={len(elig_131)} GSE205335={len(elig_205)})")

    raw = pd.concat(
        [score_patients(obj131, pairs, elig_131), score_patients(obj205, pairs, elig_205)],
        ignore_index=True,
    )
    raw.to_csv(args.outdir / "patient_pair_scores.tsv", sep="\t", index=False)

    ranked = rank_family(raw, pairs)
    barrier = ranked[ranked["family"] == "barrier_inhibitory"].copy()
    recruit = ranked[ranked["family"] == "ifn_recruit"].copy()
    barrier.to_csv(args.outdir / "barrier_inhibitory_high.tsv", sep="\t", index=False)
    recruit.to_csv(args.outdir / "ifn_recruit_low.tsv", sep="\t", index=False)
    ranked.to_csv(args.outdir / "both_families.tsv", sep="\t", index=False)

    composites = pd.DataFrame(
        [family_composite(raw, "barrier_inhibitory"), family_composite(raw, "ifn_recruit")]
    )
    composites.to_csv(args.outdir / "family_composites.tsv", sep="\t", index=False)

    plot_n_cells(ntab, figdir / "n_cells_by_patient.png")
    plot_family_bars(
        barrier,
        figdir / "barrier_inhibitory_high.png",
        "Barrier / inhibitory  (thesis: more from CLDN4-high)",
        +1,
    )
    plot_family_bars(
        recruit,
        figdir / "ifn_recruit_low.png",
        "IFN / T-recruit / MHC-I  (thesis: more from CLDN4-low)",
        -1,
    )
    plot_two_arms(barrier, recruit, figdir / "thesis_two_arms.png")
    plot_paired_boxes(
        raw,
        ["F11R–ITGAL+ITGB2", "NECTIN2–TIGIT", "CDH1–ITGAE+ITGB7", "LGALS9–PTPRC"],
        figdir / "barrier_paired_boxes.png",
        "Barrier / inhibitory pairs — patient-paired high vs low",
    )
    plot_paired_boxes(
        raw,
        ["CXCL9–CXCR3", "CXCL10–CXCR3", "CCL5–CCR5", "HLA-A–CD8A"],
        figdir / "recruit_paired_boxes.png",
        "IFN / T-recruit / MHC-I pairs — patient-paired high vs low",
    )
    plot_cohort_split(raw, pairs, figdir / "family_by_cohort.png")
    for name in (
        "n_cells_by_patient.png",
        "barrier_inhibitory_high.png",
        "ifn_recruit_low.png",
        "thesis_two_arms.png",
        "barrier_paired_boxes.png",
        "recruit_paired_boxes.png",
        "family_by_cohort.png",
    ):
        src = figdir / name
        if src.exists():
            shutil.copy2(src, args.outdir / "figures" / name)

    summary = {
        "datasets": ["GSE131907", "GSE205335"],
        "additive": True,
        "marker": "CLDN4",
        "no_dual_high": True,
        "no_gse148071": True,
        "given_pr320_q4q1": GIVEN_Q4Q1,
        "n_cells_GSE131907": obj131["n_cells_matrix"],
        "n_cells_GSE205335": obj205["n_cells_matrix"],
        "n_locked_samples_GSE131907": int(len(locked_samples)),
        "n_locked_patients_GSE131907": int(obj131["n_locked_patients_geo"]),
        "n_locked_patients_GSE205335": int(len(locked_patients)),
        "n_malig": int(obj131["is_malig"].sum() + obj205["is_malig"].sum()),
        "n_malig_GSE131907": int(obj131["is_malig"].sum()),
        "n_malig_GSE205335": int(obj205["is_malig"].sum()),
        "n_cldn4_high": int(obj131["cld_hi"].sum() + obj205["cld_hi"].sum()),
        "n_cldn4_low": int(obj131["cld_lo"].sum() + obj205["cld_lo"].sum()),
        "n_tnk": int(obj131["is_tnk"].sum() + obj205["is_tnk"].sum()),
        "n_tnk_GSE131907": int(obj131["is_tnk"].sum()),
        "n_tnk_GSE205335": int(obj205["is_tnk"].sum()),
        "cldn4_threshold_GSE131907": obj131["cldn4_threshold"],
        "cldn4_threshold_GSE205335": obj205["cldn4_threshold"],
        "n_patients_paired": n_paired,
        "n_patients_paired_GSE131907": int(len(elig_131)),
        "n_patients_paired_GSE205335": int(len(elig_205)),
        "n_pairs": int(len(pairs)),
        "genes_GSE131907": obj131["genes_kept"],
        "genes_GSE205335": obj205["genes_kept"],
        "method": "cellphonedb_mean_of_means_on_log1p_cp10k",
        "unit_of_inference": "patient; GSE131907 cells pooled by patient_id from locked samples",
        "primary_tables": ["results/barrier_inhibitory_high.tsv", "results/ifn_recruit_low.tsv"],
        "barrier_n_agree": int(barrier["agrees_thesis"].sum()),
        "recruit_n_agree": int(recruit["agrees_thesis"].sum()),
        "family_composites": composites.to_dict(orient="records"),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(outdir=args.outdir, summary=summary, ntab=ntab, barrier=barrier, recruit=recruit, composites=composites)
    (HERE / "README.md").write_text(
        "# Winning-pair thesis LR (CLDN4-only)\n\n"
        "GSE131907 + GSE205335. Both pre-specified families. "
        "See `FINDING.md`, `results/barrier_inhibitory_high.tsv`, `results/ifn_recruit_low.tsv`.\n\n"
        "```bash\n"
        "python3 methods/winpair_131907_205335_lr_thesis_cldn4/scripts/00_download.py\n"
        "python3 methods/winpair_131907_205335_lr_thesis_cldn4/scripts/01_analyze.py\n"
        "```\n"
    )
    (args.outdir / "README.md").write_text(
        "# Thesis LR tables\n\n"
        "Primary: `barrier_inhibitory_high.tsv` and `ifn_recruit_low.tsv` (n / Δ / p).\n"
        "See `../FINDING.md`.\n"
    )
    log(json.dumps({k: summary[k] for k in ("n_patients_paired", "barrier_n_agree", "recruit_n_agree")}, indent=2))
    log(f"[done] {args.outdir}")


if __name__ == "__main__":
    main()
