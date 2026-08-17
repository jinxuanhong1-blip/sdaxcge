#!/usr/bin/env python3
"""Merged GSE131907+GSE205335 CLDN4-only LIANA / CellPhoneDB.

Outgoing CLDN4-high malignant → same-patient T/NK.
Patient-level paired Wilcoxon. Honest n.
No dual-high. No GSE207422.
PR #320 Q4 vs Q1 n=23 r=−0.705 is given and is not re-audited.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import traceback
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]

TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES = {"T lymphocytes", "NK cells"}
GIVEN_Q4Q1 = {"n_compared": 23, "n_q1": 12, "n_q4": 11, "r": -0.705}


def log(msg: str) -> None:
    print(msg, flush=True)


def fmt_p(value: float) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_delta(value: float) -> str:
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
    log(f"[stream] {path} keep={len(keep)} symbols")
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
    name_to_row = {g: i for i, g in enumerate(genes)}
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


def score_pairs(
    pairs: pd.DataFrame,
    logx: dict[str, np.ndarray],
    sender: np.ndarray,
    receiver: np.ndarray,
    expr_prop: float,
) -> pd.DataFrame:
    s_mean, s_frac, n_s = group_gene_stats(logx, sender)
    r_mean, r_frac, n_r = group_gene_stats(logx, receiver)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac = partner_from_stats(lig_u, s_mean, s_frac)
        rec_m, rec_f = partner_from_stats(rec_u, r_mean, r_frac)
        if not np.isfinite(l_mean) or not np.isfinite(rec_m):
            continue
        pass_prop = (l_frac >= expr_prop) and (rec_f >= expr_prop)
        rows.append(
            {
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "pathway": rec.pathway,
                "pair_origin": rec.pair_origin,
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": rec_m,
                "ligand_frac": l_frac,
                "receptor_frac": rec_f,
                "cpdb_mean_score": 0.5 * (l_mean + rec_m),
                "product_score": l_mean * rec_m,
                "pass_expr_prop": pass_prop,
            }
        )
    return pd.DataFrame(rows)


def wilcoxon_safe(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 3 or np.allclose(a, b):
        return np.nan
    try:
        return float(stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


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


def run_liana(adata, groupby: str, out_csv: Path, n_perms: int) -> str:
    try:
        import liana as li
    except Exception as exc:
        return f"LIANA_IMPORT_FAILED: {exc}"
    try:
        li.mt.cellphonedb(
            adata,
            groupby=groupby,
            resource_name="cellphonedb",
            expr_prop=0.10,
            n_perms=n_perms,
            use_raw=False,
            verbose=True,
            key_added="liana_res",
        )
        res = adata.uns["liana_res"].copy()
        res.to_csv(out_csv, index=False)
        return f"LIANA_OK n_edges={len(res)} file={out_csv.name}"
    except Exception as exc:
        return f"LIANA_RUN_FAILED: {exc}\n{traceback.format_exc()}"


def try_liana_cohort(
    cohort: str,
    logx: dict[str, np.ndarray],
    cld_hi: np.ndarray,
    cld_lo: np.ndarray,
    is_t: np.ndarray,
    is_nk: np.ndarray,
    cells: np.ndarray,
    patient: np.ndarray,
    params: dict,
    outdir: Path,
) -> str:
    try:
        import anndata as ad
    except Exception as exc:
        return f"LIANA_IMPORT_FAILED: anndata {exc}"
    mask = cld_hi | cld_lo | is_t | is_nk
    if int(mask.sum()) < 50:
        return "LIANA_SKIP_TOO_FEW_CELLS"
    genes = sorted(logx)
    X = np.vstack([logx[g][mask] for g in genes]).T
    obs = pd.DataFrame({"patient_id": patient[mask]}, index=pd.Index(cells[mask], name="cell"))
    grp = np.array(["other"] * int(mask.sum()), dtype=object)
    grp[cld_hi[mask]] = "Malig_CLDN4high"
    grp[cld_lo[mask]] = "Malig_CLDN4low"
    grp[is_t[mask]] = "T"
    grp[is_nk[mask]] = "NK"
    obs["cc_group"] = grp
    rng = np.random.default_rng(params["random_seed"])
    keep_idx = []
    for g, idx in obs.groupby("cc_group", observed=True).indices.items():
        if g == "other":
            continue
        if len(idx) > params["liana_max_cells_per_group"]:
            idx = rng.choice(idx, size=params["liana_max_cells_per_group"], replace=False)
        keep_idx.append(np.asarray(idx))
    keep_idx = np.sort(np.concatenate(keep_idx)) if keep_idx else np.array([], dtype=int)
    adata = ad.AnnData(X=X[keep_idx], obs=obs.iloc[keep_idx].copy(), var=pd.DataFrame(index=genes))
    adata.obs["cc_group"] = pd.Categorical(adata.obs["cc_group"])
    adata.uns["log1p"] = {"base": None}
    log(f"[liana {cohort}] groups {adata.obs['cc_group'].value_counts().to_dict()}")
    note = run_liana(adata, "cc_group", outdir / f"liana_cellphonedb_{cohort}.csv", params["liana_n_perms"])
    if note.startswith("LIANA_OK"):
        res = pd.read_csv(outdir / f"liana_cellphonedb_{cohort}.csv")
        outg = res[res["source"].eq("Malig_CLDN4high") & res["target"].isin(["T", "NK"])].copy()
        outg.to_csv(outdir / f"liana_cldn4high_to_tnk_{cohort}.csv", index=False)
    return f"{cohort}: {note}"


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
    is_t = is_locked_sample & (ctype == "T lymphocytes")
    is_nk = is_locked_sample & (ctype == "NK cells")
    is_tnk = is_t | is_nk
    logx = {g: log1p_cp10k(arr, totals) for g, arr in expr.items()}
    if "CLDN4" not in logx:
        raise SystemExit("CLDN4 missing from GSE131907")
    return {
        "cohort": "GSE131907",
        "cells": cells,
        "logx": logx,
        "sample": sample,
        "patient": patient,
        "origin": origin,
        "is_malig": is_malig,
        "is_t": is_t,
        "is_nk": is_nk,
        "is_tnk": is_tnk,
        "n_genes": n_genes,
        "n_cells_matrix": int(len(cells)),
        "sample_meta": sample_meta,
        "locked_samples": sorted(locked_samples),
        "n_locked_samples": int(len(locked_samples)),
        "n_locked_patients_geo": int(pd.Series([patient_map.get(s, s) for s in locked_samples]).nunique()),
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
    is_t = is_tnk.copy()
    is_nk = np.zeros(len(cells), dtype=bool)
    logx = {g: log1p_cp10k(extracted[g], library_umi) for g in extracted}
    if "CLDN4" not in logx:
        raise SystemExit("CLDN4 missing from GSE205335")
    return {
        "cohort": "GSE205335",
        "cells": barcodes,
        "logx": logx,
        "sample": cells["orig.ident"].to_numpy(),
        "patient": cells["patient"].to_numpy(),
        "origin": cells["tissue"].fillna("").to_numpy(),
        "is_malig": is_malig,
        "is_t": is_t,
        "is_nk": is_nk,
        "is_tnk": is_tnk,
        "n_genes": int(len(all_genes)),
        "n_cells_matrix": int(len(barcodes)),
        "sample_meta": metadata,
        "locked_patients": sorted(locked_patients),
        "n_locked_patients": int(len(locked_patients)),
        "histology": cells["cancer_subtype"].to_numpy(),
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


def patient_n_table(obj: dict, min_hi: int, min_lo: int, min_tnk: int) -> pd.DataFrame:
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
                "eligible_paired": n_hi >= min_hi and n_lo >= min_lo and n_tnk >= min_tnk and n_mal > 0,
            }
        )
    return pd.DataFrame(rows)


def score_patients(obj: dict, pairs: pd.DataFrame, expr_prop: float, eligible: set[str]) -> pd.DataFrame:
    chunks = []
    for pid in sorted(eligible):
        m = obj["patient"] == pid
        h, l, t = m & obj["cld_hi"], m & obj["cld_lo"], m & obj["is_tnk"]
        hdf = score_pairs(pairs, obj["logx"], h, t, expr_prop)
        ldf = score_pairs(pairs, obj["logx"], l, t, expr_prop)
        key = ["ligand", "receptor", "pathway"]
        d = hdf[key + ["cpdb_mean_score", "pass_expr_prop", "ligand_frac", "receptor_frac", "n_sender", "n_receiver"]].merge(
            ldf[key + ["cpdb_mean_score", "pass_expr_prop", "ligand_frac", "receptor_frac", "n_sender", "n_receiver"]],
            on=key,
            suffixes=("_high", "_low"),
        )
        d["delta_high_minus_low"] = d["cpdb_mean_score_high"] - d["cpdb_mean_score_low"]
        d["pass_either"] = d["pass_expr_prop_high"] | d["pass_expr_prop_low"]
        d["pass_both"] = d["pass_expr_prop_high"] & d["pass_expr_prop_low"]
        d["cohort"] = obj["cohort"]
        d["patient_id"] = pid
        d["unit_id"] = f"{obj['cohort']}:{pid}"
        d["n_high"] = int(h.sum())
        d["n_low"] = int(l.sum())
        d["n_tnk"] = int(t.sum())
        chunks.append(d)
        log(f"  scored {obj['cohort']} {pid} high={int(h.sum())} low={int(l.sum())} tnk={int(t.sum())}")
    return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()


def rank_pairs(raw: pd.DataFrame) -> pd.DataFrame:
    rank_rows = []
    if raw.empty:
        return pd.DataFrame()
    for (lig, recp, path), sub in raw.groupby(["ligand", "receptor", "pathway"], observed=True):
        keep_sub = sub[sub["pass_either"]].copy()
        if len(keep_sub) < 3:
            continue
        if path not in {"T_recruit", "IFN", "MHC_I"} and keep_sub["pass_both"].sum() < 3:
            continue
        p = wilcoxon_safe(keep_sub["cpdb_mean_score_high"], keep_sub["cpdb_mean_score_low"])
        n_131 = int(keep_sub.loc[keep_sub["cohort"] == "GSE131907", "unit_id"].nunique())
        n_205 = int(keep_sub.loc[keep_sub["cohort"] == "GSE205335", "unit_id"].nunique())
        rank_rows.append(
            {
                "ligand": lig,
                "receptor": recp,
                "pathway": path,
                "n_patients": int(keep_sub["unit_id"].nunique()),
                "n_GSE131907": n_131,
                "n_GSE205335": n_205,
                "median_delta": float(np.median(keep_sub["delta_high_minus_low"])),
                "mean_delta": float(np.mean(keep_sub["delta_high_minus_low"])),
                "mean_score_high": float(keep_sub["cpdb_mean_score_high"].mean()),
                "mean_score_low": float(keep_sub["cpdb_mean_score_low"].mean()),
                "frac_pass_high": float(keep_sub["pass_expr_prop_high"].mean()),
                "frac_pass_low": float(keep_sub["pass_expr_prop_low"].mean()),
                "pval": p,
            }
        )
    ranks = pd.DataFrame(rank_rows)
    if len(ranks) and ranks["pval"].notna().any():
        mask = ranks["pval"].notna()
        ranks.loc[mask, "padj"] = multipletests(ranks.loc[mask, "pval"], method="fdr_bh")[1]
    else:
        ranks["padj"] = np.nan
    if len(ranks):
        ranks = ranks.sort_values(["pathway", "median_delta", "pval"], na_position="last")
    return ranks


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


def plot_focus(foc: pd.DataFrame, n_pat: int, path: Path, title: str) -> None:
    if foc.empty:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.axis("off")
        ax.text(0.5, 0.5, "No focus pairs passed the paired gate", ha="center")
        fig.savefig(path, dpi=160)
        plt.close(fig)
        return
    foc_plot = foc.sort_values("median_delta")
    colors = {"T_recruit": "#4C72B0", "IFN": "#55A868", "MHC_I": "#C44E52"}
    fig, ax = plt.subplots(figsize=(8.2, max(3.2, 0.28 * len(foc_plot) + 1.2)))
    y = np.arange(len(foc_plot))
    ax.barh(y, foc_plot["median_delta"], color=[colors.get(p, "#999") for p in foc_plot["pathway"]], edgecolor="none")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{a}–{b} ({p}, n={n})" for a, b, p, n in zip(foc_plot["ligand"], foc_plot["receptor"], foc_plot["pathway"], foc_plot["n_patients"])],
        fontsize=7,
    )
    ax.set_xlabel("median patient Δ score (CLDN4-high − low)")
    ax.set_title(f"{title}  paired n={n_pat}")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_top_hits(ranks: pd.DataFrame, path: Path, n_pat: int) -> None:
    if ranks.empty:
        return
    show = ranks.assign(_abs=ranks["median_delta"].abs()).sort_values(["pval", "_abs"], ascending=[True, False]).head(20)
    fig, ax = plt.subplots(figsize=(8.4, max(3.4, 0.32 * len(show) + 1.2)))
    y = np.arange(len(show))
    colors = np.where(show["median_delta"] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, show["median_delta"], color=colors, alpha=0.85)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{a}–{b} ({p}, n={n})" for a, b, p, n in zip(show["ligand"], show["receptor"], show["pathway"], show["n_patients"])],
        fontsize=7,
    )
    ax.set_xlabel("median patient Δ (CLDN4-high − low)")
    ax.set_title(f"Top paired hits by Wilcoxon p  n={n_pat}")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_cohort_split(raw: pd.DataFrame, focus_pairs: pd.DataFrame, path: Path) -> None:
    if raw.empty or focus_pairs.empty:
        return
    keys = focus_pairs[["ligand", "receptor", "pathway"]].drop_duplicates()
    merged = raw.merge(keys, on=["ligand", "receptor", "pathway"])
    merged = merged[merged["pass_either"]]
    if merged.empty:
        return
    rows = []
    for (lig, recp, pathw, cohort), sub in merged.groupby(["ligand", "receptor", "pathway", "cohort"], observed=True):
        rows.append(
            {
                "pair": f"{lig}–{recp}",
                "pathway": pathw,
                "cohort": cohort,
                "median_delta": float(np.median(sub["delta_high_minus_low"])),
                "n": int(sub["unit_id"].nunique()),
            }
        )
    d = pd.DataFrame(rows)
    pairs = list(d["pair"].drop_duplicates())
    fig, ax = plt.subplots(figsize=(8.6, max(3.4, 0.32 * len(pairs) + 1.4)))
    y = np.arange(len(pairs))
    for i, cohort in enumerate(["GSE131907", "GSE205335"]):
        vals = []
        for pair in pairs:
            hit = d[(d["pair"] == pair) & (d["cohort"] == cohort)]
            vals.append(float(hit["median_delta"].iloc[0]) if len(hit) else np.nan)
        ax.barh(y + (i - 0.5) * 0.36, vals, height=0.34, label=cohort, alpha=0.85)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(pairs, fontsize=7)
    ax.set_xlabel("median Δ (high − low)")
    ax.set_title("Focus pairs by cohort (descriptive; not a second test)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(
    outdir: Path,
    summary: dict,
    ntab: pd.DataFrame,
    ranks: pd.DataFrame,
    focus: pd.DataFrame,
    top: pd.DataFrame,
) -> None:
    n_paired = int(summary["n_patients_paired"])
    n_131 = int(summary["n_patients_paired_GSE131907"])
    n_205 = int(summary["n_patients_paired_GSE205335"])
    n_neg = int(summary.get("n_focus_delta_neg", 0))
    n_focus = int(summary.get("n_focus_tested", 0))
    n_sig = int(summary.get("n_focus_fdr05", 0))
    n_sig_any = int(summary.get("n_focus_fdr05_any", 0))

    if n_paired == 0:
        verdict = (
            "No patient met the paired gate (≥10 CLDN4-high malignant, "
            "≥10 CLDN4-low malignant, ≥20 T/NK). The LR table is empty."
        )
    elif n_sig == 0 and n_focus:
        verdict = (
            f"On {n_paired} paired patients ({n_131} GSE131907 + {n_205} GSE205335), "
            f"CLDN4-high vs CLDN4-low malignant → same-patient T/NK focus pairs do "
            f"**not** support a coordinated T-recruit / MHC-I drop "
            f"({n_neg}/{n_focus} median Δ < 0; **{n_sig} FDR < 0.05** weaker from high). "
            "The LR table is a ranked co-expression list, not a causal claim."
        )
    else:
        verdict = (
            f"Paired n = {n_paired} patients ({n_131} + {n_205}). "
            f"Focus pairs: {n_neg}/{n_focus} median Δ < 0; {n_sig} FDR < 0.05 weaker from high "
            f"({n_sig_any} FDR < 0.05 in either direction). "
            "Read the table; do not upgrade this to a recruitment mechanism."
        )

    foc_show = focus.sort_values(["pathway", "median_delta"]) if len(focus) else focus
    top_show = top.head(15) if len(top) else top

    lines = [
        "# FINDING — merged GSE131907 + GSE205335 LIANA/LR from CLDN4-high malignant to T/NK",
        "",
        f"**Verdict:** {verdict}",
        "",
        "ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No GSE207422. "
        "Patient is the unit. p-values are descriptive.",
        "",
        "The between-patient abundance cut on this same slice is **given** and is not re-audited: "
        f"PR #320 Q4 vs Q1 author %pos vs T/NK **n={GIVEN_Q4Q1['n_compared']} "
        f"({GIVEN_Q4Q1['n_q1']}/{GIVEN_Q4Q1['n_q4']}) r={GIVEN_Q4Q1['r']:+.3f}**. "
        "This folder tests a different contrast: outgoing CellPhoneDB-style scores from "
        "CLDN4-high vs CLDN4-low malignant cells to **same-patient** T/NK.",
        "",
        "---",
        "",
        "## Honest n",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        f"| GSE131907 cells in UMI matrix | {summary['n_cells_GSE131907']:,} | author barcodes |",
        f"| GSE205335 cells in UMI matrix | {summary['n_cells_GSE205335']:,} | author barcodes |",
        f"| GSE131907 locked samples (PR #320) | {summary['n_locked_samples_GSE131907']} | tumor-origin, n_malignant ≥20 |",
        f"| GSE131907 unique patients in locked samples | {summary['n_locked_patients_GSE131907']} | GEO patient_id; extract is sample-level |",
        f"| GSE205335 locked patients (PR #320) | {summary['n_locked_patients_GSE205335']} | ≥20 malignant and ≥20 T/NK |",
        f"| Author malignant cells (locked slice) | {summary['n_malig']:,} | {summary['n_malig_GSE131907']:,} + {summary['n_malig_GSE205335']:,} |",
        f"| CLDN4-high / low malignant | {summary['n_cldn4_high']:,} / {summary['n_cldn4_low']:,} | cohort-specific global medians |",
        f"| T/NK cells (locked slice) | {summary['n_tnk']:,} | {summary['n_tnk_GSE131907']:,} + {summary['n_tnk_GSE205335']:,} |",
        f"| **Paired patients (Wilcoxon unit)** | **{n_paired}** | ≥10 high, ≥10 low, ≥20 T/NK |",
        f"|  … GSE131907 / GSE205335 | {n_131} / {n_205} | GSE131907 cells pooled by patient_id |",
        f"| Pairs scored / ranked | {summary['n_pairs_scored']} / {summary['n_pairs_ranked']} | CellPhoneDB v5 + overlay |",
        f"| LIANA | {summary['liana_status']} | secondary; not the patient test |",
        f"| CellChat | {summary['cellchat_status']} | not run |",
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
            if rec.n_cldn4_high < summary["min_malig_per_state"]:
                reasons.append("high<10")
            if rec.n_cldn4_low < summary["min_malig_per_state"]:
                reasons.append("low<10")
            if rec.n_T_NK < summary["min_tnk"]:
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
        "`pass_expr_prop` requires both partners in ≥10% of cells in their group. "
        "The test is a **paired Wilcoxon** of high vs low **per patient**. "
        "FDR is Benjamini–Hochberg within the outgoing contrast. "
        "This is **not** a CellChat communication probability and does **not** observe secretion or spatial contact.",
        "",
        f"LIANA `mt.cellphonedb` was {'run' if str(summary['liana_status']).startswith('ran') else 'not run / did not import'}. "
        "If present, LIANA p-values are within-object specificity, not patient-level tests.",
        "",
        "## Primary LR table — outgoing CLDN4-high malignant → T/NK",
        "",
        "Median Δ = high − low. Negative = weaker from CLDN4-high. "
        f"Paired n = **{n_paired}** ({n_131} GSE131907 + {n_205} GSE205335).",
        "",
        "### MHC-I and T-recruit focus pairs",
        "",
    ]
    lines += md_table(
        foc_show,
        ["pathway", "ligand", "receptor", "n_patients", "n_GSE131907", "n_GSE205335", "median_delta", "pval", "padj"],
        {"median_delta": "{:+.3f}", "pval": "{:.3g}", "padj": "{:.3g}"},
    )
    lines += [
        f"Focus pairs with median Δ < 0: **{n_neg}/{n_focus}**. "
        f"FDR < 0.05 and Δ < 0: **{n_sig}/{n_focus}**. "
        f"FDR < 0.05 in either direction: **{n_sig_any}/{n_focus}**.",
        "",
        "### Top paired hits (all pathways, by Wilcoxon p)",
        "",
    ]
    lines += md_table(
        top_show,
        ["pathway", "ligand", "receptor", "n_patients", "n_GSE131907", "n_GSE205335", "median_delta", "pval", "padj"],
        {"median_delta": "{:+.3f}", "pval": "{:.3g}", "padj": "{:.3g}"},
    )
    lines += [
        "Full ranked table: `results/lr_table.tsv`. Per-patient scores: `results/patient_outgoing.tsv`.",
        "",
        "## What this is not",
        "",
        "- Not a re-audit of PR #320 Q4 vs Q1 n=23 r=−0.705.",
        "- Not dual-high TACSTD2×CLDN4.",
        "- Not GSE207422.",
        "- Not inferCNV/CopyKAT recomputed malignant IDs.",
        "- Not CellChat. LIANA permutation p-values (if present) are not patient tests.",
        "- Not spatial proximity or protein secretion.",
        "- Cells are not the sample size. GSE131907 PR #320 extract is sample-level; this Wilcoxon uses unique patients.",
        "",
        "## Extra figures",
        "",
        "- `figures/n_cells_by_patient.png`",
        "- `figures/cldn4_outgoing_focus.png`",
        "- `figures/top_hits_delta.png`",
        "- `figures/focus_by_cohort.png`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/merge_131907_205335_liana_cldn4/scripts/00_download.py",
        "python3 methods/merge_131907_205335_liana_cldn4/scripts/01_analyze.py",
        "```",
        "",
    ]
    (HERE / "FINDING.md").write_text("\n".join(lines))
    (outdir / "FINDING.md").write_text("\n".join(lines))


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

    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    P = cfg["params"]
    pairs = pd.read_csv(HERE / "resources" / "cellphonedb_v5_lr_pairs.tsv", sep="\t")
    lr_genes = set(Path(HERE / "resources" / "lr_genes.txt").read_text().split())
    keep = lr_genes | set(cfg["state_genes"])

    s131 = pd.read_csv(HERE / "data" / "GSE131907_samples.tsv", sep="\t")
    locked_samples = set(
        s131.loc[
            s131["origin"].isin(TUMOR_ORIGINS) & (s131["n_malignant"] >= 20),
            "sample",
        ].astype(str)
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

    ntab = pd.concat(
        [
            patient_n_table(obj131, P["min_malig_per_state"], P["min_malig_per_state"], P["min_tnk"]),
            patient_n_table(obj205, P["min_malig_per_state"], P["min_malig_per_state"], P["min_tnk"]),
        ],
        ignore_index=True,
    )
    ntab.to_csv(args.outdir / "n_cells_patients.tsv", sep="\t", index=False)
    eligible = ntab.loc[ntab["eligible_paired"], ["cohort", "patient_id"]]
    elig_131 = set(eligible.loc[eligible["cohort"] == "GSE131907", "patient_id"])
    elig_205 = set(eligible.loc[eligible["cohort"] == "GSE205335", "patient_id"])
    n_paired = int(len(elig_131) + len(elig_205))
    log(f"[n] paired patients={n_paired} (GSE131907={len(elig_131)} GSE205335={len(elig_205)})")

    raw = pd.concat(
        [
            score_patients(obj131, pairs, P["expr_prop"], elig_131),
            score_patients(obj205, pairs, P["expr_prop"], elig_205),
        ],
        ignore_index=True,
    )
    raw.to_csv(args.outdir / "patient_outgoing.tsv", sep="\t", index=False)

    ranks = rank_pairs(raw)
    ranks.to_csv(args.outdir / "lr_table.tsv", sep="\t", index=False)
    ranks.to_csv(args.outdir / "lr_table_high_vs_low.tsv", sep="\t", index=False)

    focus = ranks[ranks["pathway"].isin(["T_recruit", "MHC_I"])].copy() if len(ranks) else ranks
    focus.to_csv(args.outdir / "lr_table_focus_outgoing.tsv", sep="\t", index=False)
    n_focus = int(len(focus))
    n_neg = int((focus["median_delta"] < 0).sum()) if n_focus else 0
    n_sig = int(((focus["padj"] < 0.05) & (focus["median_delta"] < 0)).sum()) if n_focus and focus["padj"].notna().any() else 0
    n_sig_any = int((focus["padj"] < 0.05).sum()) if n_focus and focus["padj"].notna().any() else 0
    top = ranks.assign(_abs=ranks["median_delta"].abs()).sort_values(["pval", "_abs"], ascending=[True, False]).drop(columns="_abs") if len(ranks) else ranks
    top.head(50).to_csv(args.outdir / "lr_table_top_hits.tsv", sep="\t", index=False)

    # Pooled descriptive (not the test): CLDN4-high malignant → T/NK stacked within cohort then listed
    pooled_rows = []
    for obj, elig in ((obj131, elig_131), (obj205, elig_205)):
        if not elig:
            continue
        m = np.isin(obj["patient"], list(elig))
        pooled = score_pairs(pairs, obj["logx"], obj["cld_hi"] & m, obj["is_tnk"] & m, P["expr_prop"])
        pooled["cohort"] = obj["cohort"]
        pooled_rows.append(pooled)
    pooled = pd.concat(pooled_rows, ignore_index=True) if pooled_rows else pd.DataFrame()
    if len(pooled):
        pooled.to_csv(args.outdir / "lr_table_pooled_descriptive.tsv", sep="\t", index=False)

    liana_notes = []
    try:
        liana_notes.append(
            try_liana_cohort(
                "GSE131907",
                obj131["logx"],
                obj131["cld_hi"],
                obj131["cld_lo"],
                obj131["is_t"],
                obj131["is_nk"],
                obj131["cells"],
                obj131["patient"],
                P,
                args.outdir,
            )
        )
    except Exception as exc:
        liana_notes.append(f"GSE131907: LIANA_WRAPPER_FAILED: {exc}")
    try:
        liana_notes.append(
            try_liana_cohort(
                "GSE205335",
                obj205["logx"],
                obj205["cld_hi"],
                obj205["cld_lo"],
                obj205["is_t"],
                obj205["is_nk"],
                obj205["cells"],
                obj205["patient"],
                P,
                args.outdir,
            )
        )
    except Exception as exc:
        liana_notes.append(f"GSE205335: LIANA_WRAPPER_FAILED: {exc}")
    if any("LIANA_OK" in n for n in liana_notes):
        liana_status = "ran_cellphonedb_method"
    elif any("LIANA_IMPORT_FAILED" in n for n in liana_notes):
        liana_status = "LIANA_IMPORT_FAILED"
    else:
        liana_status = liana_notes[0].split(":")[0] if liana_notes else "not_run"

    plot_n_cells(ntab, figdir / "n_cells_by_patient.png")
    plot_focus(focus, n_paired, figdir / "cldn4_outgoing_focus.png", "Outgoing malignant → T/NK focus")
    plot_top_hits(ranks, figdir / "top_hits_delta.png", n_paired)
    plot_cohort_split(raw, focus, figdir / "focus_by_cohort.png")
    for name in ("n_cells_by_patient.png", "cldn4_outgoing_focus.png", "top_hits_delta.png", "focus_by_cohort.png"):
        src = figdir / name
        if src.exists():
            shutil.copy2(src, args.outdir / "figures" / name)

    n_table = pd.DataFrame(
        [
            {"item": "cells_GSE131907", "n": obj131["n_cells_matrix"], "note": "author barcodes"},
            {"item": "cells_GSE205335", "n": obj205["n_cells_matrix"], "note": "author barcodes"},
            {"item": "locked_samples_GSE131907", "n": len(locked_samples), "note": "PR #320 tumor-origin n_mal>=20"},
            {"item": "locked_patients_GSE131907", "n": obj131["n_locked_patients_geo"], "note": "unique patient_id in 21 samples"},
            {"item": "locked_patients_GSE205335", "n": len(locked_patients), "note": "PR #320 extract"},
            {"item": "malignant_cells", "n": int(obj131["is_malig"].sum() + obj205["is_malig"].sum()), "note": "locked slice"},
            {"item": "cldn4_high", "n": int(obj131["cld_hi"].sum() + obj205["cld_hi"].sum()), "note": "cohort medians"},
            {"item": "cldn4_low", "n": int(obj131["cld_lo"].sum() + obj205["cld_lo"].sum()), "note": "below cohort median"},
            {"item": "T_NK", "n": int(obj131["is_tnk"].sum() + obj205["is_tnk"].sum()), "note": "locked slice"},
            {"item": "paired_patients", "n": n_paired, "note": "unit of Wilcoxon"},
            {"item": "paired_GSE131907", "n": int(len(elig_131)), "note": "patient_id pooled"},
            {"item": "paired_GSE205335", "n": int(len(elig_205)), "note": "patient"},
            {"item": "pairs_ranked", "n": int(len(ranks)), "note": "outgoing high vs low"},
        ]
    )
    n_table.to_csv(args.outdir / "n_table.tsv", sep="\t", index=False)

    summary = {
        "datasets": ["GSE131907", "GSE205335"],
        "additive": True,
        "marker": "CLDN4",
        "no_dual_high": True,
        "no_gse207422": True,
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
        "min_malig_per_state": P["min_malig_per_state"],
        "min_tnk": P["min_tnk"],
        "n_pairs_scored": int(raw[["ligand", "receptor"]].drop_duplicates().shape[0]) if len(raw) else 0,
        "n_pairs_ranked": int(len(ranks)),
        "n_focus_tested": n_focus,
        "n_focus_delta_neg": n_neg,
        "n_focus_fdr05": n_sig,
        "n_focus_fdr05_any": n_sig_any,
        "liana_status": liana_status,
        "liana_notes": liana_notes,
        "cellchat_status": "not_run_R_unavailable",
        "method": "cellphonedb_mean_of_means_on_log1p_cp10k",
        "unit_of_inference": "patient; GSE131907 cells pooled by patient_id from locked samples",
        "paired_patients": ntab.loc[ntab["eligible_paired"], ["cohort", "patient_id", "n_cldn4_high", "n_cldn4_low", "n_T_NK"]].to_dict(orient="records"),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    write_finding(args.outdir, summary, ntab, ranks, focus, top)
    (args.outdir / "README.md").write_text(
        "# Merged GSE131907+GSE205335 CLDN4-high malignant → T/NK LR\n\n"
        "See `../FINDING.md` and `lr_table.tsv`.\n"
    )
    log(json.dumps({k: summary[k] for k in ("n_patients_paired", "n_patients_paired_GSE131907", "n_patients_paired_GSE205335", "n_pairs_ranked", "liana_status")}, indent=2))
    log(f"[done] {args.outdir}")


if __name__ == "__main__":
    main()
