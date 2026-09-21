#!/usr/bin/env python3
"""GSE131907 MPE vs primary: CLDN4 with TACSTD2, ELF3, EPCAM, and T/NK.

HRA006761 is controlled access and is not in this run.
Author Cell_subtype 'Malignant cells' is absent from pleural effusion and from
primary tLung (those epithelial states are tS1/tS2/tS3). MPE uses the author
epithelial label. That distinction is kept in every table.
"""
from __future__ import annotations

import argparse
import gzip
import json
from itertools import permutations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, rankdata, spearmanr

PARTNERS = ["TACSTD2", "ELF3", "EPCAM"]
AXIS = ["CLDN4", "TACSTD2", "ELF3", "EPCAM"]
MIN_EPI = 20
SEED = 131907


def parse_series(path: Path) -> pd.DataFrame:
    title = geo = patient = stage = origin = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0][len("!Sample_") :]
            vals = [v.strip('"') for v in parts[1:]]
            if key == "title":
                title = vals
            elif key == "geo_accession":
                geo = vals
            elif key == "characteristics_ch1":
                if vals and vals[0].startswith("patient id:"):
                    patient = [v.split(": ", 1)[1] for v in vals]
                elif vals and vals[0].startswith("tumor stage:"):
                    stage = [v.split(": ", 1)[1] for v in vals]
                elif vals and vals[0].startswith("tissue origin abbrevation:"):
                    origin = [v.split(": ", 1)[1] for v in vals]
    if any(v is None for v in (title, geo, patient, stage, origin)):
        raise SystemExit("series matrix missing sample fields")
    return pd.DataFrame(
        {
            "Sample": title,
            "geo_accession": geo,
            "patient_id": patient,
            "tumor_stage": stage,
            "Sample_Origin_geo": origin,
        }
    )


def site_group(origin: str) -> str:
    return {
        "PE": "MPE",
        "tLung": "primary",
        "tL/B": "advanced_biopsy",
        "mLN": "advanced_biopsy",
        "mBrain": "brain_met",
        "nLung": "normal_lung",
        "nLN": "normal_ln",
    }.get(origin, "other")


def lognorm(counts: np.ndarray, libsize: np.ndarray) -> np.ndarray:
    scale = np.divide(10000.0, libsize, out=np.zeros_like(libsize), where=libsize > 0)
    return np.log1p(counts.astype(np.float64) * scale)


def rank_corr(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(x)
    ry = rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt(np.dot(rx, rx) * np.dot(ry, ry))
    if den == 0:
        return float("nan")
    return float(np.dot(rx, ry) / den)


def spearman_test(x: np.ndarray, y: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    if n < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return {"n": n, "rho": float("nan"), "p": float("nan"), "method": "undefined"}
    rho = rank_corr(x, y)
    # Exact two-sided permutation for the sample sizes in this MPE comparison.
    if n <= 7:
        obs = abs(rho)
        count = 0
        total = 0
        for perm in permutations(y.tolist()):
            total += 1
            r = rank_corr(x, np.asarray(perm, dtype=float))
            if abs(r) + 1e-9 >= obs:
                count += 1
        return {"n": n, "rho": rho, "p": count / total, "method": "exact_permutation"}
    res = spearmanr(x, y)
    return {
        "n": n,
        "rho": float(res.statistic),
        "p": float(res.pvalue),
        "method": "spearmanr_asymptotic",
    }


def fisher_codetect(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a) > 0
    b = np.asarray(b) > 0
    if a.size == 0:
        return {
            "n": 0,
            "n11": 0,
            "n10": 0,
            "n01": 0,
            "n00": 0,
            "pct_partner_given_cldn4": float("nan"),
            "pct_partner_given_cldn4neg": float("nan"),
            "odds_ratio": float("nan"),
            "fisher_p": float("nan"),
        }
    n11 = int(np.sum(a & b))
    n10 = int(np.sum(a & ~b))
    n01 = int(np.sum(~a & b))
    n00 = int(np.sum(~a & ~b))
    table = np.array([[n11, n10], [n01, n00]])
    oddsratio, p = fisher_exact(table, alternative="two-sided")
    p_b_given_a = n11 / (n11 + n10) if (n11 + n10) else float("nan")
    p_b_given_nota = n01 / (n01 + n00) if (n01 + n00) else float("nan")
    return {
        "n": int(a.size),
        "n11": n11,
        "n10": n10,
        "n01": n01,
        "n00": n00,
        "pct_partner_given_cldn4": 100.0 * p_b_given_a,
        "pct_partner_given_cldn4neg": 100.0 * p_b_given_nota,
        "odds_ratio": float(oddsratio),
        "fisher_p": float(p),
    }


def mean_or_nan(x: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    return float(np.mean(x))


def pct_pos(x: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    return float(100.0 * np.mean(x > 0))


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, float_format="%.6g")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--subset", type=Path, default=Path("/tmp/gse131907/mpe_axis"))
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("methods/gse131907_mpe_axis/results"),
    )
    args = ap.parse_args()
    tab = args.outdir / "tables"
    fig = args.outdir / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)

    blob = np.load(args.subset / "gene_counts.npz", allow_pickle=False)
    genes = [str(g) for g in blob["genes"].tolist()]
    barcodes = [str(b) for b in blob["barcodes"].tolist()]
    counts = blob["counts"].astype(np.float64)
    libsize = blob["libsize"].astype(np.float64)
    gidx = {g: i for i, g in enumerate(genes)}
    logx = np.vstack([lognorm(counts[gidx[g]], libsize) for g in genes])

    ann = pd.read_csv(
        args.datadir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    series = parse_series(args.datadir / "GSE131907_series_matrix.txt.gz")
    if ann["Index"].duplicated().any():
        raise SystemExit("duplicate annotation Index")
    ann = ann.set_index("Index", drop=False)
    missing = [b for b in barcodes if b not in ann.index]
    if missing:
        raise SystemExit(f"{len(missing)} matrix barcodes missing from annotation")
    if len(barcodes) != len(ann):
        raise SystemExit(f"cell count mismatch matrix {len(barcodes)} ann {len(ann)}")
    ann = ann.loc[barcodes].reset_index(drop=True)
    meta = series.set_index("Sample")
    ann["patient_id"] = ann["Sample"].map(meta["patient_id"])
    ann["tumor_stage"] = ann["Sample"].map(meta["tumor_stage"])
    if ann["patient_id"].isna().any():
        raise SystemExit("patient_id join failed")
    if not np.array_equal(meta.loc[ann["Sample"], "Sample_Origin_geo"].to_numpy(), ann["Sample_Origin"].to_numpy()):
        raise SystemExit("Sample_Origin disagrees with GEO series matrix")

    origin = ann["Sample_Origin"].to_numpy()
    ctype = ann["Cell_type"].to_numpy()
    csub = ann["Cell_subtype"].to_numpy()
    sample = ann["Sample"].to_numpy()
    patient = ann["patient_id"].to_numpy()

    pe_epi = (origin == "PE") & (ctype == "Epithelial cells")
    primary_ts = (origin == "tLung") & np.isin(csub, ["tS1", "tS2", "tS3"])
    primary_epi = (origin == "tLung") & (ctype == "Epithelial cells")
    author_mal = csub == "Malignant cells"
    at2 = (origin == "nLung") & (csub == "AT2")
    adv_epi = np.isin(origin, ["tL/B", "mLN"]) & (ctype == "Epithelial cells")
    is_tnk = np.isin(ctype, ["T lymphocytes", "NK cells"])

    # Annotation facts already counted from the public cell table.
    checks = {
        "n_cells": int(len(ann)),
        "n_pe_epithelial": int(pe_epi.sum()),
        "n_author_malignant_in_pe": int((author_mal & (origin == "PE")).sum()),
        "n_author_malignant_in_tlung": int((author_mal & (origin == "tLung")).sum()),
        "n_primary_ts": int(primary_ts.sum()),
        "n_author_malignant": int(author_mal.sum()),
        "n_at2": int(at2.sum()),
    }
    expected = {
        "n_cells": 208506,
        "n_pe_epithelial": 396,
        "n_author_malignant_in_pe": 0,
        "n_author_malignant_in_tlung": 0,
        "n_primary_ts": 6352,
        "n_author_malignant": 24784,
        "n_at2": 2020,
    }
    if checks != expected:
        raise SystemExit(f"annotation checks failed: {checks} vs {expected}")

    umi = {g: counts[gidx[g]] for g in genes}
    logg = {g: logx[gidx[g]] for g in genes}

    # Author subtype is NA on every PE epithelial cell. Marker split:
    # carcinoma-like = EPCAM+ and WT1- and CALB2-; mesothelial-like = EPCAM-
    # and (WT1+ or CALB2+). EPCAM is part of the carcinoma-like gate, so an
    # EPCAM co-detection odds ratio inside that gate is not interpretable.
    mpe_carcinoma = pe_epi & (umi["EPCAM"] > 0) & (umi["WT1"] == 0) & (umi["CALB2"] == 0)
    mpe_mesothelial = pe_epi & (umi["EPCAM"] == 0) & ((umi["WT1"] > 0) | (umi["CALB2"] > 0))
    # Same marker gate on author primary tumor states, so an MPE-vs-primary
    # mean is not created by requiring EPCAM only on the MPE side.
    primary_matched = primary_ts & (umi["EPCAM"] > 0) & (umi["WT1"] == 0) & (umi["CALB2"] == 0)
    if int(mpe_carcinoma.sum()) != 259 or int(mpe_mesothelial.sum()) != 76:
        raise SystemExit(
            f"unexpected MPE marker split: carcinoma {int(mpe_carcinoma.sum())} mesothelial {int(mpe_mesothelial.sum())}"
        )

    compartments = {
        "MPE_epithelial": pe_epi,
        "MPE_epithelial_EPCAM_pos": pe_epi & (umi["EPCAM"] > 0),
        "MPE_carcinoma_like": mpe_carcinoma,
        "MPE_mesothelial_like": mpe_mesothelial,
        "primary_tS": primary_ts,
        "primary_epithelial_all": primary_epi,
        "normal_AT2": at2,
        "author_malignant_advanced": author_mal,
    }

    # --- sample table ---
    rows = []
    for samp, idx in pd.Series(np.arange(len(ann))).groupby(sample, sort=False).groups.items():
        idx = np.asarray(list(idx))
        ori = origin[idx][0]
        for comp_name, mask_all in (
            ("MPE_epithelial", pe_epi),
            ("MPE_carcinoma_like", mpe_carcinoma),
            ("MPE_mesothelial_like", mpe_mesothelial),
            ("primary_tS", primary_ts),
            ("primary_tS_matched_gate", primary_matched),
            ("epithelial_same_gate", ctype == "Epithelial cells"),
        ):
            if comp_name.startswith("MPE_") and ori != "PE":
                continue
            if comp_name in ("primary_tS", "primary_tS_matched_gate") and ori != "tLung":
                continue
            if comp_name == "epithelial_same_gate" and ori not in ("PE", "tLung", "tL/B", "mLN", "nLung"):
                continue
            comp = idx[mask_all[idx]]
            rec = {
                "sample": samp,
                "patient_id": patient[idx][0],
                "origin": ori,
                "site_group": site_group(ori),
                "tumor_stage": ann["tumor_stage"].to_numpy()[idx][0],
                "compartment": comp_name,
                "n_cells": int(idx.size),
                "n_tnk": int(is_tnk[idx].sum()),
                "frac_tnk": float(is_tnk[idx].mean()),
                "n_epithelial": int((ctype[idx] == "Epithelial cells").sum()),
                "n_author_malignant": int(author_mal[idx].sum()),
                "n_comp": int(comp.size),
            }
            for g in AXIS + ["KRT19", "MSLN", "CALB2", "WT1"]:
                rec[f"{g}_pct"] = pct_pos(umi[g][comp])
                rec[f"{g}_mean_log"] = mean_or_nan(logg[g][comp])
            rows.append(rec)
    sample_df = pd.DataFrame(rows)
    write_tsv(sample_df, tab / "sample_compartment.tsv")

    # --- cell co-detection and spearman ---
    code_rows = []
    spear_rows = []
    for comp_name, mask in compartments.items():
        for partner in PARTNERS:
            rec = fisher_codetect(umi["CLDN4"][mask], umi[partner][mask])
            rec.update({"compartment": comp_name, "gene": partner})
            code_rows.append(rec)
            both = mask & (umi["CLDN4"] > 0) & (umi[partner] > 0)
            sp_all = spearman_test(logg["CLDN4"][mask], logg[partner][mask])
            sp_both = spearman_test(logg["CLDN4"][both], logg[partner][both])
            spear_rows.append(
                {
                    "compartment": comp_name,
                    "gene": partner,
                    "n_cells": int(mask.sum()),
                    "spearman_log_all_rho": sp_all["rho"],
                    "spearman_log_all_p": sp_all["p"],
                    "spearman_log_all_method": sp_all["method"],
                    "n_double_pos": int(both.sum()),
                    "spearman_log_double_pos_rho": sp_both["rho"],
                    "spearman_log_double_pos_p": sp_both["p"],
                    "spearman_log_double_pos_method": sp_both["method"],
                }
            )
    code_df = pd.DataFrame(code_rows)
    spear_df = pd.DataFrame(spear_rows)
    write_tsv(code_df, tab / "cell_codetection.tsv")
    write_tsv(spear_df, tab / "cell_spearman.tsv")

    # --- sample pseudobulk spearman and MPE vs primary ---
    def eligible(comp: str, site: str) -> pd.DataFrame:
        d = sample_df[(sample_df.compartment == comp) & (sample_df.site_group == site) & (sample_df.n_comp >= MIN_EPI)]
        return d.copy()

    pseudo_rows = []
    for comp, site in (
        ("MPE_epithelial", "MPE"),
        ("MPE_carcinoma_like", "MPE"),
        ("primary_tS", "primary"),
        ("epithelial_same_gate", "MPE"),
        ("epithelial_same_gate", "primary"),
    ):
        d = eligible(comp, site)
        for partner in PARTNERS:
            sp = spearman_test(d[f"CLDN4_mean_log"].to_numpy(), d[f"{partner}_mean_log"].to_numpy())
            pseudo_rows.append(
                {
                    "compartment": comp,
                    "site_group": site,
                    "gene": partner,
                    "min_comp_cells": MIN_EPI,
                    "n_samples": sp["n"],
                    "rho": sp["rho"],
                    "p": sp["p"],
                    "method": sp["method"],
                    "x": "CLDN4_mean_log",
                    "y": f"{partner}_mean_log",
                }
            )
    pseudo_df = pd.DataFrame(pseudo_rows)
    write_tsv(pseudo_df, tab / "sample_pseudobulk_spearman.tsv")

    mw_rows = []
    contrasts = [
        ("MPE_epithelial", "primary_tS", "MPE_epithelial_vs_primary_tS"),
        ("MPE_carcinoma_like", "primary_tS", "MPE_carcinoma_like_vs_primary_tS"),
        ("MPE_carcinoma_like", "primary_tS_matched_gate", "MPE_carcinoma_like_vs_primary_matched_gate"),
        ("epithelial_same_gate", "epithelial_same_gate", "same_gate_epithelial_MPE_vs_primary"),
    ]
    for comp_a, comp_b, name in contrasts:
        a = eligible(comp_a, "MPE")
        b = eligible(comp_b, "primary")
        for g in AXIS:
            for metric in (f"{g}_mean_log", f"{g}_pct"):
                xa = a[metric].to_numpy()
                xb = b[metric].to_numpy()
                if len(xa) >= 3 and len(xb) >= 3 and np.isfinite(xa).all() and np.isfinite(xb).all():
                    stat = mannwhitneyu(xa, xb, alternative="two-sided", method="exact")
                    p = float(stat.pvalue)
                    u = float(stat.statistic)
                else:
                    p, u = float("nan"), float("nan")
                mw_rows.append(
                    {
                        "contrast": name,
                        "metric": metric,
                        "n_mpe": int(len(xa)),
                        "n_primary": int(len(xb)),
                        "median_mpe": float(np.median(xa)) if len(xa) else float("nan"),
                        "median_primary": float(np.median(xb)) if len(xb) else float("nan"),
                        "delta_median_mpe_minus_primary": float(np.median(xa) - np.median(xb))
                        if len(xa) and len(xb)
                        else float("nan"),
                        "mannwhitney_U": u,
                        "p": p,
                    }
                )
    mw_df = pd.DataFrame(mw_rows)
    write_tsv(mw_df, tab / "mpe_vs_primary.tsv")

    # Immune inverse inside MPE. Primary pre-specified set: n_comp >= 20.
    immune_rows = []
    mpe_all = sample_df[sample_df.compartment == "MPE_epithelial"].copy()
    mpe_carc = sample_df[sample_df.compartment == "MPE_carcinoma_like"].copy()
    for label, d in (
        ("epithelial_all_MPE_samples", mpe_all),
        ("epithelial_n_ge_20", mpe_all[mpe_all.n_comp >= MIN_EPI]),
        ("carcinoma_like_n_comp_gt_0", mpe_carc[mpe_carc.n_comp > 0]),
        ("carcinoma_like_n_ge_20", mpe_carc[mpe_carc.n_comp >= MIN_EPI]),
    ):
        sp_pct = spearman_test(d["CLDN4_pct"].to_numpy(), d["frac_tnk"].to_numpy())
        sp_mean = spearman_test(d["CLDN4_mean_log"].to_numpy(), d["frac_tnk"].to_numpy())
        immune_rows.append(
            {
                "set": label,
                "n_samples": int(len(d)),
                "rho_cldn4_pct_vs_frac_tnk": sp_pct["rho"],
                "p_cldn4_pct_vs_frac_tnk": sp_pct["p"],
                "method_pct": sp_pct["method"],
                "rho_cldn4_mean_vs_frac_tnk": sp_mean["rho"],
                "p_cldn4_mean_vs_frac_tnk": sp_mean["p"],
                "method_mean": sp_mean["method"],
                "median_frac_tnk": float(d["frac_tnk"].median()) if len(d) else float("nan"),
                "min_frac_tnk": float(d["frac_tnk"].min()) if len(d) else float("nan"),
                "max_frac_tnk": float(d["frac_tnk"].max()) if len(d) else float("nan"),
                "median_cldn4_pct": float(d["CLDN4_pct"].median()) if len(d) else float("nan"),
            }
        )
    pri_ge = sample_df[(sample_df.compartment == "primary_tS") & (sample_df.n_comp >= MIN_EPI)]
    sp_pct = spearman_test(pri_ge["CLDN4_pct"].to_numpy(), pri_ge["frac_tnk"].to_numpy())
    sp_mean = spearman_test(pri_ge["CLDN4_mean_log"].to_numpy(), pri_ge["frac_tnk"].to_numpy())
    immune_rows.append(
        {
            "set": "primary_tS_n_ge_20",
            "n_samples": int(len(pri_ge)),
            "rho_cldn4_pct_vs_frac_tnk": sp_pct["rho"],
            "p_cldn4_pct_vs_frac_tnk": sp_pct["p"],
            "method_pct": sp_pct["method"],
            "rho_cldn4_mean_vs_frac_tnk": sp_mean["rho"],
            "p_cldn4_mean_vs_frac_tnk": sp_mean["p"],
            "method_mean": sp_mean["method"],
            "median_frac_tnk": float(pri_ge["frac_tnk"].median()) if len(pri_ge) else float("nan"),
            "min_frac_tnk": float(pri_ge["frac_tnk"].min()) if len(pri_ge) else float("nan"),
            "max_frac_tnk": float(pri_ge["frac_tnk"].max()) if len(pri_ge) else float("nan"),
            "median_cldn4_pct": float(pri_ge["CLDN4_pct"].median()) if len(pri_ge) else float("nan"),
        }
    )
    immune_df = pd.DataFrame(immune_rows)
    write_tsv(immune_df, tab / "mpe_cldn4_vs_tnk.tsv")
    write_tsv(
        mpe_all[
            [
                "sample",
                "patient_id",
                "tumor_stage",
                "n_cells",
                "n_comp",
                "n_tnk",
                "frac_tnk",
                "n_author_malignant",
                "CLDN4_pct",
                "CLDN4_mean_log",
                "TACSTD2_pct",
                "TACSTD2_mean_log",
                "ELF3_pct",
                "ELF3_mean_log",
                "EPCAM_pct",
                "EPCAM_mean_log",
                "KRT19_pct",
                "MSLN_pct",
            ]
        ].sort_values("sample"),
        tab / "mpe_sample_detail.tsv",
    )
    detail_cols = [
        "sample",
        "patient_id",
        "n_cells",
        "n_comp",
        "n_tnk",
        "frac_tnk",
        "CLDN4_pct",
        "CLDN4_mean_log",
        "TACSTD2_pct",
        "TACSTD2_mean_log",
        "ELF3_pct",
        "ELF3_mean_log",
        "EPCAM_pct",
        "EPCAM_mean_log",
        "KRT19_pct",
        "MSLN_pct",
        "CALB2_pct",
        "WT1_pct",
    ]
    write_tsv(mpe_carc[detail_cols].sort_values("sample"), tab / "mpe_carcinoma_like_detail.tsv")
    write_tsv(
        sample_df.loc[sample_df.compartment == "MPE_mesothelial_like", detail_cols].sort_values("sample"),
        tab / "mpe_mesothelial_like_detail.tsv",
    )

    # Paired MPE vs matched advanced biopsy, same epithelial gate, n patients with both.
    epi_gate = sample_df[sample_df.compartment == "epithelial_same_gate"]
    paired_rows = []
    pe_samples = epi_gate[epi_gate.origin == "PE"]
    for _, pe_row in pe_samples.iterrows():
        match = epi_gate[(epi_gate.patient_id == pe_row.patient_id) & (epi_gate.site_group == "advanced_biopsy")]
        for _, bx in match.iterrows():
            rec = {
                "patient_id": pe_row.patient_id,
                "mpe_sample": pe_row["sample"],
                "biopsy_sample": bx["sample"],
                "biopsy_origin": bx["origin"],
                "n_mpe_epi": int(pe_row.n_comp),
                "n_biopsy_epi": int(bx.n_comp),
            }
            for g in AXIS:
                rec[f"{g}_mean_mpe"] = float(pe_row[f"{g}_mean_log"])
                rec[f"{g}_mean_biopsy"] = float(bx[f"{g}_mean_log"])
                rec[f"{g}_delta_mpe_minus_biopsy"] = float(pe_row[f"{g}_mean_log"] - bx[f"{g}_mean_log"])
            paired_rows.append(rec)
    paired_df = pd.DataFrame(paired_rows)
    write_tsv(paired_df, tab / "paired_mpe_vs_biopsy.tsv")

    # Identity summaries: cell-level medians of log expression.
    id_rows = []
    for comp_name, mask in {
        "MPE_epithelial": pe_epi,
        "MPE_carcinoma_like": mpe_carcinoma,
        "MPE_mesothelial_like": mpe_mesothelial,
        "primary_tS": primary_ts,
        "normal_AT2": at2,
        "author_malignant_advanced": author_mal,
        "advanced_biopsy_epithelial": adv_epi,
    }.items():
        rec = {"compartment": comp_name, "n_cells": int(mask.sum())}
        for g in AXIS + ["KRT19", "MSLN", "CALB2", "WT1", "KRT7"]:
            rec[f"{g}_pct"] = pct_pos(umi[g][mask])
            rec[f"{g}_mean_log"] = mean_or_nan(logg[g][mask])
            rec[f"{g}_median_log"] = float(np.median(logg[g][mask])) if mask.any() else float("nan")
        id_rows.append(rec)
    id_df = pd.DataFrame(id_rows)
    write_tsv(id_df, tab / "identity_cell_summary.tsv")

    # --- figures ---
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 160,
            "savefig.dpi": 160,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    c_mpe = "#c45c26"
    c_pri = "#2c5f8a"

    # 1. CLDN4 % vs T/NK. MPE x-axis is the carcinoma-like gate.
    # EFFUSION_13 has zero carcinoma-like cells and is omitted.
    fig1, ax = plt.subplots(figsize=(6.4, 4.6))
    pri = sample_df[(sample_df.compartment == "primary_tS")]
    mpe = sample_df[(sample_df.compartment == "MPE_carcinoma_like") & (sample_df.n_comp > 0)]
    ax.scatter(
        pri["CLDN4_pct"],
        pri["frac_tnk"],
        s=np.clip(pri["n_comp"] / 8, 28, 140),
        c=c_pri,
        alpha=0.75,
        label="Primary tS1/tS2/tS3",
        zorder=2,
    )
    ax.scatter(
        mpe["CLDN4_pct"],
        mpe["frac_tnk"],
        s=np.clip(mpe["n_comp"] * 1.4, 40, 180),
        c=c_mpe,
        alpha=0.95,
        label="MPE carcinoma-like",
        zorder=3,
        edgecolors="black",
        linewidths=0.4,
    )
    for _, r in mpe.iterrows():
        ax.annotate(
            f"{r['sample'].replace('EFFUSION_', 'E')} (n={int(r['n_comp'])})",
            (r["CLDN4_pct"], r["frac_tnk"]),
            textcoords="offset points",
            xytext=(5, 4),
            fontsize=8,
            color=c_mpe,
        )
    ax.set_xlabel("CLDN4+ % in compartment")
    ax.set_ylabel("T/NK fraction of the sample")
    xmin = float(np.nanmin([pri["CLDN4_pct"].min(), mpe["CLDN4_pct"].min()]))
    ax.set_xlim(min(-3, xmin - 5), 105)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, loc="lower left")
    ax.set_title("GSE131907: CLDN4+ % vs T/NK")
    fig1.tight_layout()
    fig1.savefig(fig / "fig_cldn4_vs_tnk.png")
    fig1.savefig(fig / "fig_cldn4_vs_tnk.pdf")
    plt.close(fig1)

    # 2. Co-detection bars. EPCAM is the MPE carcinoma-like gate, so it is omitted there.
    fig2, axes = plt.subplots(1, 2, figsize=(8.2, 4.2), sharey=True)
    panels = (
        (axes[0], "MPE_carcinoma_like", "MPE carcinoma-like", ["TACSTD2", "ELF3"]),
        (axes[1], "primary_tS", "Primary tS1/tS2/tS3", PARTNERS),
    )
    for ax, comp, title, partners in panels:
        sub = code_df[code_df.compartment == comp].set_index("gene").loc[partners]
        x = np.arange(len(partners))
        w = 0.36
        ax.bar(x - w / 2, sub["pct_partner_given_cldn4"], width=w, color="#1b4f72", label="CLDN4+")
        ax.bar(x + w / 2, sub["pct_partner_given_cldn4neg"], width=w, color="#d5d8dc", label="CLDN4−")
        ax.set_xticks(x, partners)
        ax.set_title(title)
        ax.set_ylabel("% cells detected" if ax is axes[0] else "")
        ax.set_ylim(0, 100)
    axes[1].legend(frameon=False, loc="upper right")
    fig2.suptitle("Partner detection given CLDN4 status", y=1.02)
    fig2.tight_layout()
    fig2.savefig(fig / "fig_codetection.png", bbox_inches="tight")
    fig2.savefig(fig / "fig_codetection.pdf", bbox_inches="tight")
    plt.close(fig2)

    # 3. Sample mean log, MPE vs primary
    fig3, ax = plt.subplots(figsize=(7.2, 4.4))
    mpe_e = eligible("MPE_carcinoma_like", "MPE")
    pri_e = eligible("primary_tS_matched_gate", "primary")
    rng = np.random.default_rng(SEED)
    for i, g in enumerate(AXIS):
        y1 = mpe_e[f"{g}_mean_log"].to_numpy()
        y2 = pri_e[f"{g}_mean_log"].to_numpy()
        ax.scatter(i - 0.12 + rng.uniform(-0.03, 0.03, size=y1.size), y1, c=c_mpe, s=36, zorder=3)
        ax.scatter(i + 0.12 + rng.uniform(-0.03, 0.03, size=y2.size), y2, c=c_pri, s=28, zorder=2, alpha=0.85)
        ax.plot([i - 0.22, i - 0.02], [np.median(y1), np.median(y1)], color=c_mpe, lw=2)
        ax.plot([i + 0.02, i + 0.22], [np.median(y2), np.median(y2)], color=c_pri, lw=2)
    ax.scatter([], [], c=c_mpe, s=36, label=f"MPE carcinoma-like, n={len(mpe_e)}")
    ax.scatter([], [], c=c_pri, s=28, label=f"Primary tS, same gate, n={len(pri_e)}")
    ax.set_xticks(range(len(AXIS)), AXIS)
    ax.set_ylabel("Sample mean log1p(CP10k)")
    ax.legend(frameon=False, loc="best")
    ax.set_title("Same gate (EPCAM+ WT1− CALB2−), samples with ≥20 cells")
    fig3.tight_layout()
    fig3.savefig(fig / "fig_mpe_vs_primary_means.png")
    fig3.savefig(fig / "fig_mpe_vs_primary_means.pdf")
    plt.close(fig3)

    # 4. Identity distributions for CLDN4
    fig4, ax = plt.subplots(figsize=(6.8, 4.4))
    order = [
        ("normal AT2", at2, "#7f8c8d"),
        ("primary tS", primary_ts, c_pri),
        ("MPE carcinoma-like", mpe_carcinoma, c_mpe),
        ("MPE mesothelial-like", mpe_mesothelial, "#6c3483"),
    ]
    data = [logg["CLDN4"][m] for _, m, _ in order]
    parts = ax.violinplot(data, showmedians=True, showextrema=False)
    for body, (_, _, col) in zip(parts["bodies"], order):
        body.set_facecolor(col)
        body.set_alpha(0.8)
        body.set_edgecolor("black")
        body.set_linewidth(0.4)
    parts["cmedians"].set_color("black")
    ax.set_xticks(range(1, 5), [lab for lab, _, _ in order])
    ax.set_ylabel("CLDN4 log1p(CP10k)")
    ax.set_title("CLDN4 in annotated compartments")
    fig4.tight_layout()
    fig4.savefig(fig / "fig_cldn4_identity.png")
    fig4.savefig(fig / "fig_cldn4_identity.pdf")
    plt.close(fig4)

    summary = {
        "dataset": "GSE131907",
        "citation": "Kim et al. Nat Commun 2020",
        "primary_dataset_attempt": {
            "accession": "HRA006761",
            "repository": "GSA-Human",
            "bioproject": "PRJCA023797",
            "accessibility": "controlled",
            "dac": "HDAC002197",
            "analyzed": False,
            "reason": "GSA-Human page states controlled access / request data. No public matrix.",
        },
        "normalization": "log1p(UMI / full_library_size * 10000)",
        "min_compartment_cells_for_sample_tests": MIN_EPI,
        "annotation_checks": checks,
        "n_mpe_samples": int((sample_df.compartment == "MPE_epithelial").sum()),
        "n_mpe_samples_ge_20_epi": int(((sample_df.compartment == "MPE_epithelial") & (sample_df.n_comp >= MIN_EPI)).sum()),
        "n_primary_ts_samples_ge_20": int(((sample_df.compartment == "primary_tS") & (sample_df.n_comp >= MIN_EPI)).sum()),
        "author_malignant_cells_in_mpe": 0,
        "mpe_carcinoma_like_gate": "PE epithelial AND EPCAM>0 AND WT1==0 AND CALB2==0",
        "mpe_mesothelial_like_gate": "PE epithelial AND EPCAM==0 AND (WT1>0 OR CALB2>0)",
        "n_mpe_carcinoma_like_cells": int(mpe_carcinoma.sum()),
        "n_mpe_mesothelial_like_cells": int(mpe_mesothelial.sum()),
        "n_mpe_carcinoma_like_samples_ge_20": int(((sample_df.compartment == "MPE_carcinoma_like") & (sample_df.n_comp >= MIN_EPI)).sum()),
        "note": "MPE epithelial Cell_subtype is NA. They are not the author Malignant cells label used in the locked GSE131907 n=21 analysis. Carcinoma-like vs mesothelial-like is a marker gate on that unlabeled epithelial set.",
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("codetection MPE / primary")
    print(
        code_df[code_df.compartment.isin(["MPE_epithelial", "MPE_carcinoma_like", "primary_tS"])].to_string(
            index=False
        )
    )
    print("immune")
    print(immune_df.to_string(index=False))
    print("mw CLDN4/TACSTD2/ELF3/EPCAM mean")
    print(mw_df[mw_df.metric.str.endswith("_mean_log")].to_string(index=False))


if __name__ == "__main__":
    main()
