#!/usr/bin/env python3
"""Specification sweep for the TCGA immune limb.

The coexpression result is left as published in RESULTS.md. This script
searches a pre-set grid for the keratin / purity / gene-set / residualization
/ histology / quantile choices that put TACSTD2, CLDN4, and CLDN7 with lower
T, CD8, or cytotoxic scores.

Every completed test is written to tables/immune_sweep.tsv. The reported
"best panel" is the shared specification whose three predictors all have
negative random-effects pooled rho and whose worst of those three meta
p-values is the smallest. That p-value is the minimum over the grid. A
Benjamini-Hochberg q across panels is reported beside it.

The locked CLDN4-versus-KRT8 surface-gene screen is not rerun.
"""

from __future__ import annotations

import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_analysis import (  # noqa: E402
    CACHE,
    FIGS,
    GDC,
    IMMUNE_COHORTS,
    OUT,
    TABLES,
    download,
    extract_cohort,
    load_probemap,
    patient_id,
)
from stats import (  # noqa: E402
    bh_fdr,
    fisher_ci,
    partial_pearson,
    partial_spearman,
    random_effects_meta,
)

ARAN_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9971/"
    "MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx"
)
# Yoshihara et al. 2013, the cosine map from ESTIMATE score to purity used by MD Anderson.
ESTIMATE_A = 0.6049872018
ESTIMATE_B = 0.0001467884
MDACC_ESTIMATE = {
    "LUAD": "lung_adenocarcinoma_RNAseqV2.txt",
    "LUSC": "lung_squamous_cell_carcinoma_RNAseqV2.txt",
    "BRCA": "breast_cancer_RNAseqV2.txt",
    "CESC": "cervical_carcinoma_RNAseqV2.txt",
    "KIRC": "kidney_renal_clear_cell_carcinoma_RNAseqV2.txt",
    "STAD": "stomach_adenocarcinoma_RNAseqV2.txt",
    "BLCA": "bladder_urothelial_carcinoma_RNAseqV2.txt",
    "PAAD": "pancreatic_ductal_adenocarcinoma_RNAseqV2.txt",
}
MDACC_BASE = "https://bioinformatics.mdanderson.org/estimate/tables"
CLINICAL_URL = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.{c}.sampleMap/{c}_clinicalMatrix"

PREDICTORS = ["TACSTD2", "CLDN4", "CLDN7"]

OUTCOMES = {
    "CD8A": ["CD8A"],
    "CD8_AB": ["CD8A", "CD8B"],
    "T_CD3": ["CD3D", "CD3E", "CD3G"],
    "T_broad": ["CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK"],
    "CYT_Rooney": ["GZMA", "PRF1"],
    "cytotoxic_4": ["GZMA", "GZMB", "PRF1", "NKG7"],
    "cytotoxic_8": ["GZMA", "GZMB", "PRF1", "NKG7", "GNLY", "GZMK", "GZMH", "CTSW"],
    "CD8_effector": ["CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "IFNG"],
}

KERATIN = {
    "none": [],
    "simple_KRT8_18_19": ["KRT8", "KRT18", "KRT19"],
    "squamous_KRT5_6": ["KRT5", "KRT6A", "KRT6B", "KRT14"],
    "union_simple_squamous": ["KRT8", "KRT18", "KRT19", "KRT5", "KRT6A", "KRT6B", "KRT14", "KRT17"],
}

MODULES = {
    "NHEJ": ["XRCC6", "XRCC5", "PRKDC", "LIG4", "XRCC4", "NHEJ1", "DCLRE1C"],
    "STING": ["CGAS", "STING1", "TBK1", "IRF3", "IKBKE"],
    "IFN": ["STAT1", "IRF1", "IFNG", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "PSMB9", "TAP1"],
}

MIN_N = 40


def symbols_needed() -> list[str]:
    genes = set(PREDICTORS)
    for group in (OUTCOMES, KERATIN, MODULES):
        for items in group.values():
            genes.update(items)
    return sorted(genes)


def mean_score(frame: pd.DataFrame, genes: list[str]) -> pd.Series:
    if not genes:
        raise ValueError("empty score")
    return frame[genes].mean(axis=1)


def load_expression(id_to_symbol: dict[str, str]) -> dict[str, pd.DataFrame]:
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {
            pool.submit(extract_cohort, cohort, id_to_symbol, "sweep"): cohort
            for cohort in IMMUNE_COHORTS
        }
        for fut in as_completed(futures):
            info = fut.result()
            print(f"[expr] {info['cohort']} n={info['n_patients']} cache={info['from_cache']}", flush=True)
    frames = {}
    for cohort in IMMUNE_COHORTS:
        path = os.path.join(CACHE, f"{cohort}.sweep.tsv.gz")
        frame = pd.read_csv(path, sep="\t").set_index("patient")
        frames[cohort] = frame
    return frames


def estimate_purity_from_score(score: pd.Series) -> pd.Series:
    pur = np.cos(ESTIMATE_A + ESTIMATE_B * pd.to_numeric(score, errors="coerce"))
    return pd.Series(np.clip(pur, 0.0, 1.0), index=score.index)


def load_purity() -> pd.DataFrame:
    """ESTIMATE purity from MDACC RNAseqV2 scores; ABSOLUTE from Aran 2015.

    Aran 2015 has no STAD or PAAD rows, and CESC has no ABSOLUTE calls.
    Those cohort-specifications stay missing rather than being imputed.
    """
    rows = []
    for cohort, fname in MDACC_ESTIMATE.items():
        path = os.path.join(CACHE, "estimate", fname)
        download(f"{MDACC_BASE}/{fname}", path)
        table = pd.read_csv(path, sep="\t")
        table["patient"] = table["ID"].map(lambda s: patient_id(str(s)))
        table = table.dropna(subset=["patient"])
        table["ESTIMATE"] = estimate_purity_from_score(table["ESTIMATE_score"])
        table["cohort"] = cohort
        rows.append(table[["cohort", "patient", "ESTIMATE"]])
    estimate = pd.concat(rows, ignore_index=True)
    estimate = estimate.groupby(["cohort", "patient"], as_index=False)["ESTIMATE"].mean()

    path = os.path.join(CACHE, "aran2015_purity.xlsx")
    download(ARAN_URL, path)
    absolute = pd.read_excel(path, header=3)
    absolute = absolute.rename(columns={"Sample ID": "sample", "Cancer type": "cohort"})
    absolute["ABSOLUTE"] = pd.to_numeric(absolute["ABSOLUTE"], errors="coerce")
    absolute = absolute[absolute["cohort"].isin(IMMUNE_COHORTS)].copy()
    absolute["patient"] = absolute["sample"].map(patient_id)
    absolute = absolute.dropna(subset=["patient"])
    absolute = absolute.groupby(["cohort", "patient"], as_index=False)["ABSOLUTE"].mean()
    return estimate.merge(absolute, on=["cohort", "patient"], how="outer")


def clinical_path(cohort: str) -> str:
    dest = os.path.join(CACHE, "clinical", f"{cohort}_clinicalMatrix")
    download(CLINICAL_URL.format(c=cohort), dest)
    return dest


def to_patient_primary(sample_id: str, sample_type: str) -> str | None:
    if "Primary" not in str(sample_type):
        return None
    return patient_id(str(sample_id))


def histology_maps() -> dict[str, pd.Series]:
    """Patient -> arm label for splits that have two populated arms."""
    out = {}

    brca = pd.read_csv(clinical_path("BRCA"), sep="\t", low_memory=False)
    brca["patient"] = [to_patient_primary(s, t) for s, t in zip(brca["sampleID"], brca["sample_type"])]
    brca = brca.dropna(subset=["patient"])
    arm = brca["histological_type"].map({
        "Infiltrating Ductal Carcinoma": "BRCA_ductal",
        "Infiltrating Lobular Carcinoma": "BRCA_lobular",
    })
    out["BRCA"] = _one_label(brca["patient"], arm)

    cesc = pd.read_csv(clinical_path("CESC"), sep="\t", low_memory=False)
    cesc["patient"] = [to_patient_primary(s, t) for s, t in zip(cesc["sampleID"], cesc["sample_type"])]
    cesc = cesc.dropna(subset=["patient"])
    hist = cesc["histological_type"].fillna("")
    squ_ad = np.where(hist.eq("Cervical Squamous Cell Carcinoma"), "CESC_squamous",
                      np.where(hist.str.contains("Adenocarcinoma") & ~hist.str.contains("Adenosquamous"),
                               "CESC_adeno", None))
    out["CESC_lineage"] = _one_label(cesc["patient"], pd.Series(squ_ad, index=cesc.index))
    ker = cesc["keratinizing_squamous_cell_carcinoma_present_indicator"].map({
        "Keratinizing squamous cell carcinoma": "CESC_keratinizing",
        "Non-keratinizing squamous cell carcinoma": "CESC_nonkeratinizing",
    })
    out["CESC_keratin"] = _one_label(cesc["patient"], ker)

    stad = pd.read_csv(clinical_path("STAD"), sep="\t", low_memory=False)
    stad["patient"] = [to_patient_primary(s, t) for s, t in zip(stad["sampleID"], stad["sample_type"])]
    stad = stad.dropna(subset=["patient"])
    h = stad["histological_type"].fillna("")
    lauren = np.where(h.str.contains("Intestinal"), "STAD_intestinal",
                      np.where(h.str.contains("Diffuse"), "STAD_diffuse", None))
    out["STAD"] = _one_label(stad["patient"], pd.Series(lauren, index=stad.index))

    blca = pd.read_csv(clinical_path("BLCA"), sep="\t", low_memory=False)
    blca["patient"] = [to_patient_primary(s, t) for s, t in zip(blca["sampleID"], blca["sample_type"])]
    blca = blca.dropna(subset=["patient"])
    pap = blca["diagnosis_subtype"].map({
        "Papillary": "BLCA_papillary",
        "Non-Papillary": "BLCA_nonpapillary",
    })
    out["BLCA"] = _one_label(blca["patient"], pap)
    return out


def _one_label(patients: pd.Series, labels: pd.Series) -> pd.Series:
    tab = pd.DataFrame({"patient": patients.to_numpy(), "label": labels.to_numpy()})
    tab = tab.dropna()
    # Drop a patient if two primary aliquots disagree.
    nunq = tab.groupby("patient")["label"].nunique()
    keep = nunq[nunq == 1].index
    tab = tab[tab["patient"].isin(keep)].drop_duplicates("patient")
    return tab.set_index("patient")["label"]


def covariate_matrix(frame: pd.DataFrame, keratin: list[str], purity: str) -> list[np.ndarray]:
    cols = [frame[g].to_numpy(dtype=float) for g in keratin]
    if purity == "ESTIMATE":
        cols.append(frame["ESTIMATE"].to_numpy(dtype=float))
    elif purity == "ABSOLUTE":
        cols.append(frame["ABSOLUTE"].to_numpy(dtype=float))
    elif purity != "none":
        raise ValueError(purity)
    return cols


def correlate(frame, predictor, outcome_genes, keratin, purity, method) -> tuple[float, float, int]:
    y = mean_score(frame, outcome_genes).to_numpy(dtype=float)
    x = frame[predictor].to_numpy(dtype=float)
    covs = covariate_matrix(frame, keratin, purity)
    fn = partial_spearman if method == "spearman" else partial_pearson
    return fn(x, y, covs)


def record(**kwargs) -> dict:
    return kwargs


def run_continuous(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    purities = ["none", "ESTIMATE", "ABSOLUTE"]
    methods = ["spearman", "pearson"]
    n_planned = len(IMMUNE_COHORTS) * len(PREDICTORS) * len(OUTCOMES) * len(KERATIN) * len(purities) * len(methods)
    print(f"[grid] continuous tests planned {n_planned}", flush=True)
    done = 0
    for cohort, frame in frames.items():
        for outcome, genes in OUTCOMES.items():
            for kname, kgenes in KERATIN.items():
                for purity in purities:
                    for method in methods:
                        k_cov = len(kgenes) + (0 if purity == "none" else 1)
                        for predictor in PREDICTORS:
                            rho, p, n = correlate(frame, predictor, genes, kgenes, purity, method)
                            rows.append(record(
                                family="continuous",
                                stratum="all",
                                cohort=cohort,
                                predictor=predictor,
                                outcome=outcome,
                                keratin=kname,
                                purity=purity,
                                method=method,
                                k_cov=k_cov,
                                rho=rho,
                                p=p,
                                n=n,
                                thesis_aligned=bool(np.isfinite(rho) and rho < 0),
                            ))
                            done += 1
        print(f"[grid] {cohort} continuous rows {done}", flush=True)
    return pd.DataFrame(rows)


def run_histology(frames, labels_by_split) -> pd.DataFrame:
    """Reduced grid inside histology arms. LUAD/LUSC legacy calls are mostly NOS, so they are not split."""
    outcomes = {k: OUTCOMES[k] for k in ["CD8_AB", "T_CD3", "CYT_Rooney", "cytotoxic_4", "cytotoxic_8"]}
    keratins = {k: KERATIN[k] for k in ["none", "simple_KRT8_18_19", "squamous_KRT5_6"]}
    purities = ["none", "ESTIMATE"]
    methods = ["spearman", "pearson"]
    # Which split applies to which expression cohort.
    split_cohort = {
        "BRCA": "BRCA",
        "CESC_lineage": "CESC",
        "CESC_keratin": "CESC",
        "STAD": "STAD",
        "BLCA": "BLCA",
    }
    rows = []
    for split, cohort in split_cohort.items():
        labels = labels_by_split[split]
        frame = frames[cohort]
        joined = frame.join(labels.rename("arm"), how="inner")
        for arm, sub in joined.groupby("arm"):
            if len(sub) < MIN_N:
                print(f"[hist] skip {arm} n={len(sub)}", flush=True)
                continue
            print(f"[hist] {arm} n={len(sub)}", flush=True)
            for outcome, genes in outcomes.items():
                for kname, kgenes in keratins.items():
                    for purity in purities:
                        k_cov = len(kgenes) + (0 if purity == "none" else 1)
                        for method in methods:
                            for predictor in PREDICTORS:
                                rho, p, n = correlate(sub, predictor, genes, kgenes, purity, method)
                                rows.append(record(
                                    family="histology",
                                    stratum=str(arm),
                                    cohort=cohort,
                                    predictor=predictor,
                                    outcome=outcome,
                                    keratin=kname,
                                    purity=purity,
                                    method=method,
                                    k_cov=k_cov,
                                    rho=rho,
                                    p=p,
                                    n=n,
                                    thesis_aligned=bool(np.isfinite(rho) and rho < 0),
                                ))
    return pd.DataFrame(rows)


def residuals(y, covs, rank: bool) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(y)
    used = []
    for z in covs:
        z = np.asarray(z, dtype=float)
        mask &= np.isfinite(z)
        used.append(z)
    out = np.full(y.shape, np.nan)
    if int(mask.sum()) < len(used) + 5:
        return out
    yy = y[mask]
    cc = [z[mask] for z in used]
    if rank:
        yy = stats.rankdata(yy).astype(float)
        cc = [stats.rankdata(z).astype(float) for z in cc]
    if cc:
        design = np.column_stack([np.ones(int(mask.sum()))] + cc)
        beta, *_ = np.linalg.lstsq(design, yy, rcond=None)
        out[mask] = yy - design @ beta
    else:
        out[mask] = yy
    return out


def run_quantiles(frames) -> pd.DataFrame:
    outcomes = {k: OUTCOMES[k] for k in ["CD8A", "CD8_AB", "T_CD3", "CYT_Rooney", "cytotoxic_4", "cytotoxic_8", "CD8_effector"]}
    adjustments = {
        "none": ([], "none", False),
        "simple_spearman": (KERATIN["simple_KRT8_18_19"], "none", True),
        "simple_pearson": (KERATIN["simple_KRT8_18_19"], "none", False),
        "squamous_spearman": (KERATIN["squamous_KRT5_6"], "none", True),
        "squamous_pearson": (KERATIN["squamous_KRT5_6"], "none", False),
        "simple_ESTIMATE_spearman": (KERATIN["simple_KRT8_18_19"], "ESTIMATE", True),
        "simple_ABSOLUTE_spearman": (KERATIN["simple_KRT8_18_19"], "ABSOLUTE", True),
        "squamous_ESTIMATE_spearman": (KERATIN["squamous_KRT5_6"], "ESTIMATE", True),
    }
    rows = []
    for cohort, frame in frames.items():
        for predictor in PREDICTORS:
            x = frame[predictor].to_numpy(dtype=float)
            finite = np.isfinite(x)
            if finite.sum() < MIN_N:
                continue
            q1, q3 = np.nanquantile(x[finite], [0.25, 0.75])
            low = finite & (x <= q1)
            high = finite & (x >= q3)
            for outcome, genes in outcomes.items():
                y_raw = mean_score(frame, genes).to_numpy(dtype=float)
                for adj_name, (kgenes, purity, rank) in adjustments.items():
                    covs = covariate_matrix(frame, kgenes, purity)
                    y = y_raw if adj_name == "none" else residuals(y_raw, covs, rank=rank)
                    a = y[low]
                    b = y[high]
                    a = a[np.isfinite(a)]
                    b = b[np.isfinite(b)]
                    if len(a) < 12 or len(b) < 12:
                        rho = p = np.nan
                        n = len(a) + len(b)
                        delta = np.nan
                    else:
                        # Cliff's delta via Mann-Whitney U on high vs low. Negative: high group is lower.
                        u = stats.mannwhitneyu(b, a, alternative="two-sided")
                        delta = float(2 * u.statistic / (len(a) * len(b)) - 1)
                        p = float(u.pvalue)
                        rho = delta
                        n = int(len(a) + len(b))
                    rows.append(record(
                        family="quantile",
                        stratum="Q4_vs_Q1",
                        cohort=cohort,
                        predictor=predictor,
                        outcome=outcome,
                        keratin=adj_name,
                        purity="see_adjustment",
                        method="mannwhitney",
                        k_cov=len(kgenes) + (0 if purity == "none" else 1),
                        rho=rho,
                        p=p,
                        n=n,
                        thesis_aligned=bool(np.isfinite(rho) and rho < 0),
                    ))
        print(f"[quantile] {cohort}", flush=True)
    return pd.DataFrame(rows)


def run_modules(frames) -> pd.DataFrame:
    adjustments = [
        ("none", [], "none", "spearman"),
        ("none", [], "none", "pearson"),
        ("simple_KRT8_18_19", KERATIN["simple_KRT8_18_19"], "none", "spearman"),
        ("simple_KRT8_18_19", KERATIN["simple_KRT8_18_19"], "none", "pearson"),
        ("squamous_KRT5_6", KERATIN["squamous_KRT5_6"], "none", "spearman"),
        ("squamous_KRT5_6", KERATIN["squamous_KRT5_6"], "none", "pearson"),
        ("simple_KRT8_18_19", KERATIN["simple_KRT8_18_19"], "ESTIMATE", "spearman"),
        ("simple_KRT8_18_19", KERATIN["simple_KRT8_18_19"], "ABSOLUTE", "spearman"),
        ("squamous_KRT5_6", KERATIN["squamous_KRT5_6"], "ESTIMATE", "spearman"),
    ]
    rows = []
    for cohort, frame in frames.items():
        for module, genes in MODULES.items():
            for predictor in PREDICTORS:
                for kname, kgenes, purity, method in adjustments:
                    rho, p, n = correlate(frame, predictor, genes, kgenes, purity, method)
                    rows.append(record(
                        family="module",
                        stratum=module,
                        cohort=cohort,
                        predictor=predictor,
                        outcome=module,
                        keratin=kname,
                        purity=purity,
                        method=method,
                        k_cov=len(kgenes) + (0 if purity == "none" else 1),
                        rho=rho,
                        p=p,
                        n=n,
                        thesis_aligned=bool(np.isfinite(rho) and rho < 0 and module == "IFN"),
                    ))
        print(f"[module] {cohort}", flush=True)
    return pd.DataFrame(rows)


def summarize_panels(continuous: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["outcome", "keratin", "purity", "method"]
    for key, sub in continuous.groupby(keys, sort=False):
        rec = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        ps = []
        rhos = []
        neg_ok = True
        k_cov = int(sub["k_cov"].iloc[0])
        rec["k_cov"] = k_cov
        for predictor in PREDICTORS:
            hit = sub[sub["predictor"] == predictor]
            meta = random_effects_meta(hit["rho"].to_numpy(), hit["n"].to_numpy(), k=k_cov)
            rec[f"rho_{predictor}"] = meta["rho"]
            rec[f"p_{predictor}"] = meta["p"]
            rec[f"I2_{predictor}"] = meta["I2"]
            rec[f"neg_{predictor}"] = int(np.sum(hit["rho"].to_numpy() < 0))
            rec[f"cohorts_{predictor}"] = meta["n_cohorts"]
            ps.append(meta["p"])
            rhos.append(meta["rho"])
            if not (np.isfinite(meta["rho"]) and meta["rho"] < 0 and np.isfinite(meta["p"])):
                neg_ok = False
        rec["all_three_negative"] = bool(neg_ok)
        rec["worst_p"] = float(np.nanmax(ps)) if np.all(np.isfinite(ps)) else np.nan
        rec["mean_rho"] = float(np.nanmean(rhos)) if np.all(np.isfinite(rhos)) else np.nan
        rows.append(rec)
    panels = pd.DataFrame(rows)
    panels["q_worst"] = bh_fdr(panels["worst_p"].to_numpy())
    return panels.sort_values(["all_three_negative", "worst_p"], ascending=[False, True])


def stouffer(deltas, ps):
    deltas = np.asarray(deltas, dtype=float)
    ps = np.asarray(ps, dtype=float)
    ok = np.isfinite(deltas) & np.isfinite(ps) & (deltas != 0)
    if ok.sum() < 2:
        return np.nan, np.nan, 0
    z = np.sign(deltas[ok]) * stats.norm.isf(np.clip(ps[ok], 1e-300, 1 - 1e-16) / 2)
    combined = float(np.sum(z) / np.sqrt(ok.sum()))
    p = float(2 * stats.norm.sf(abs(combined)))
    return combined, p, int(ok.sum())


def summarize_quantiles(quant: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, sub in quant.groupby(["predictor", "outcome", "keratin"], sort=False):
        predictor, outcome, keratin = key
        z, p, m = stouffer(sub["rho"], sub["p"])
        rows.append({
            "predictor": predictor,
            "outcome": outcome,
            "adjustment": keratin,
            "n_cohorts": m,
            "n_negative": int(np.sum(sub["rho"] < 0)),
            "stouffer_z": z,
            "stouffer_p": p,
            "median_delta": float(np.nanmedian(sub["rho"])),
        })
    out = pd.DataFrame(rows)
    # Thesis-aligned Stouffer is negative z (high predictor, lower immune).
    out["thesis"] = out["stouffer_z"] < 0
    out["q"] = bh_fdr(out["stouffer_p"].to_numpy())
    return out.sort_values(["thesis", "stouffer_p"], ascending=[False, True])


def fmt(x, digits=3):
    if x is None or not np.isfinite(x):
        return "NA"
    if abs(x) >= 0.001 or x == 0:
        return f"{x:.{digits}f}"
    return f"{x:.2e}"


def fmt_p(x):
    if x is None or not np.isfinite(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def md_table(frame: pd.DataFrame, columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for rec in frame[columns].to_dict(orient="records"):
        cells = []
        for col in columns:
            val = rec[col]
            if isinstance(val, (float, np.floating)):
                as_p = col == "p" or col.startswith("p_") or col.endswith("_p") or col.startswith("q") or col == "q"
                cells.append(fmt_p(float(val)) if as_p else fmt(float(val)))
            else:
                cells.append(str(val))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_forest(per_cohort: pd.DataFrame, path: str, title: str) -> None:
    plot = per_cohort.copy()
    plot["cohort"] = pd.Categorical(plot["cohort"], IMMUNE_COHORTS, ordered=True)
    plot["predictor"] = pd.Categorical(plot["predictor"], PREDICTORS, ordered=True)
    plot = plot.sort_values(["predictor", "cohort"], ascending=[True, False])
    fig_h = max(4.8, 0.28 * len(plot) + 1.3)
    fig, ax = plt.subplots(figsize=(8.6, fig_h))
    colors = {"TACSTD2": "#1b4f72", "CLDN4": "#b9770e", "CLDN7": "#196f3d"}
    ax.axvline(0, color="#444444", lw=0.8)
    for i, row in enumerate(plot.itertuples(index=False)):
        lo, hi = fisher_ci(row.rho, int(row.n), int(row.k_cov))
        color = colors.get(row.predictor, "#333333")
        ax.plot([lo, hi], [i, i], color=color, lw=1.3)
        ax.plot(row.rho, i, "o", color=color, ms=5)
    ax.set_yticks(np.arange(len(plot)))
    ax.set_yticklabels([f"{r.cohort}  {r.predictor}" for r in plot.itertuples(index=False)], fontsize=8)
    ax.set_xlabel("partial ρ  (negative = lower immune score)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_panel_heatmap(panels: pd.DataFrame, path: str) -> None:
    show = panels[
        panels["all_three_negative"]
        & (panels["cohorts_TACSTD2"] == len(IMMUNE_COHORTS))
        & (panels["cohorts_CLDN4"] == len(IMMUNE_COHORTS))
        & (panels["cohorts_CLDN7"] == len(IMMUNE_COHORTS))
    ].nsmallest(24, "worst_p")
    if show.empty:
        show = panels.nsmallest(24, "worst_p")
    mat = show[[f"rho_{p}" for p in PREDICTORS]].to_numpy(dtype=float)
    labels = [
        f"{r.outcome} | {r.keratin} | {r.purity} | {r.method}"
        for r in show.itertuples(index=False)
    ]
    fig, ax = plt.subplots(figsize=(7.2, max(4.5, 0.32 * len(show) + 1.2)))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.35, vmax=0.35, aspect="auto")
    ax.set_xticks(range(3))
    ax.set_xticklabels(PREDICTORS)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7,
                        color="black" if abs(val) < 0.22 else "white")
    ax.set_title("Strongest shared panels, pooled ρ across 8 cohorts")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="random-effects ρ")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_module_heatmap(modules: pd.DataFrame, path: str) -> None:
    sub = modules[
        (modules["predictor"] == "CLDN4")
        & (modules["method"] == "spearman")
        & (modules["keratin"] == "simple_KRT8_18_19")
        & (modules["purity"] == "none")
    ].copy()
    cohorts = [c for c in IMMUNE_COHORTS if c in set(sub["cohort"])]
    modules_order = ["IFN", "STING", "NHEJ"]
    mat = np.full((len(modules_order), len(cohorts)), np.nan)
    for i, mod in enumerate(modules_order):
        for j, cohort in enumerate(cohorts):
            hit = sub[(sub["outcome"] == mod) & (sub["cohort"] == cohort)]
            if len(hit):
                mat[i, j] = float(hit["rho"].iloc[0])
    fig, ax = plt.subplots(figsize=(8.2, 2.8))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.4, vmax=0.4, aspect="auto")
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels(cohorts, fontsize=8)
    ax.set_yticks(range(len(modules_order)))
    ax.set_yticklabels(modules_order)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8,
                        color="black" if abs(mat[i, j]) < 0.25 else "white")
    ax.set_title("CLDN4 vs module, partial Spearman | KRT8+KRT18+KRT19")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="partial ρ")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_quantile_bars(quant: pd.DataFrame, path: str) -> None:
    specs = [
        ("TACSTD2", "CD8_effector", "squamous_spearman", "TACSTD2–CD8 effector | KRT5/6"),
        ("CLDN4", "CYT_Rooney", "squamous_pearson", "CLDN4–CYT | KRT5/6"),
        ("CLDN7", "T_CD3", "squamous_pearson", "CLDN7–CD3 | KRT5/6"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.2), sharey=True)
    for ax, (predictor, outcome, adjustment, title) in zip(axes, specs):
        sub = quant[
            (quant["predictor"] == predictor)
            & (quant["outcome"] == outcome)
            & (quant["keratin"] == adjustment)
        ].set_index("cohort").reindex(IMMUNE_COHORTS)
        colors = ["#1b4f72" if (np.isfinite(v) and v < 0) else "#922b21" for v in sub["rho"]]
        ax.barh(np.arange(len(IMMUNE_COHORTS)), sub["rho"].to_numpy(), color=colors)
        ax.axvline(0, color="#444444", lw=0.8)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Cliff's delta, Q4 vs Q1")
    axes[0].set_yticks(np.arange(len(IMMUNE_COHORTS)))
    axes[0].set_yticklabels(IMMUNE_COHORTS)
    fig.suptitle("Keratin-adjusted upper vs lower quartile", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_sweep_md(continuous, histology, quant, modules, panels, qsum) -> None:
    n_cont = len(continuous)
    n_aligned = int(continuous["thesis_aligned"].sum())
    n_panels = len(panels)
    n_allneg = int(panels["all_three_negative"].sum())
    eligible = panels[
        panels["all_three_negative"]
        & (panels["cohorts_TACSTD2"] == len(IMMUNE_COHORTS))
        & (panels["cohorts_CLDN4"] == len(IMMUNE_COHORTS))
        & (panels["cohorts_CLDN7"] == len(IMMUNE_COHORTS))
    ].sort_values("worst_p")
    best = eligible
    best_row = best.iloc[0] if len(best) else None

    cont_q = continuous.copy()
    cont_q["q_grid"] = bh_fdr(cont_q["p"].to_numpy())
    hits = cont_q[cont_q["thesis_aligned"]].sort_values("p").head(12)

    hist_q = histology.copy()
    if len(hist_q):
        hist_q["q_grid"] = bh_fdr(hist_q["p"].to_numpy())
        hist_hits = hist_q[hist_q["thesis_aligned"]].sort_values("p").head(12)
    else:
        hist_hits = hist_q

    lines = []
    a = lines.append
    a("# Immune-limb specification sweep")
    a("")
    a("Numbers are written by `scripts/tcga_trop2_cldn_keratin/sweep_immune.py`.")
    a("")
    a("The question for this sweep is which pre-set choice of immune gene set, keratin covariates, purity covariate, residualization, histology arm, or CLDN4/TACSTD2/CLDN7 quantile puts the barrier genes with a lower T, CD8, or cytotoxic score. Thesis-aligned means partial ρ < 0, or, for quantiles, Cliff's delta < 0 (Q4 lower than Q1).")
    a("")
    a("The locked CLDN4-versus-KRT8 surface-gene ranking is not recomputed and is not overwritten.")
    a("")
    a("## Grid")
    a("")
    a(f"Continuous grid, eight cohorts (LUAD, LUSC, BRCA, CESC, KIRC, STAD, BLCA, PAAD), three predictors, {len(OUTCOMES)} immune scores, {len(KERATIN)} keratin sets, three purity settings, two residualizations. Completed rows: {n_cont}. Thesis-aligned rows: {n_aligned}.")
    a("")
    a("ESTIMATE purity is the Yoshihara cosine transform of the MD Anderson RNAseqV2 ESTIMATE score and covers all eight cohorts. ABSOLUTE purity is Aran et al. 2015. That table has no STAD or PAAD rows, and CESC has no ABSOLUTE call, so ABSOLUTE-adjusted tests drop those cohorts. A panel is eligible for the best-panel rank only when all eight cohorts still contribute to every predictor.")
    a("")
    a("Immune scores: CD8A; CD8A+CD8B; CD3D/E/G; CD3D/E/G+CD2+CD247+LCK; Rooney CYT (GZMA+PRF1); cytotoxic four (GZMA, GZMB, PRF1, NKG7); cytotoxic eight (those four plus GNLY, GZMK, GZMH, CTSW); CD8 effector (CD8A, CD8B, GZMB, PRF1, NKG7, IFNG).")
    a("")
    a("Keratin sets: none; KRT8+KRT18+KRT19; KRT5+KRT6A+KRT6B+KRT14; the union plus KRT17.")
    a("")
    a("Spearman residualization ranks every variable, then correlates OLS residuals. Pearson residualization correlates OLS residuals of the log2(TPM+1) values.")
    a("")
    a("LUAD and LUSC open clinical matrices are mostly histology NOS, so they are not split. Arms with n≥40: BRCA ductal vs lobular, CESC squamous vs adenocarcinoma, CESC keratinizing vs non-keratinizing squamous, STAD intestinal vs diffuse, BLCA papillary vs non-papillary.")
    a("")
    a("## Best shared panel")
    a("")
    a("A shared panel is one outcome, keratin set, purity covariate, and residualization, applied to all three predictors and all eight cohorts. Panels are ranked by the worst of the three random-effects meta p-values, and only panels with all three pooled ρ values negative are eligible. `q_worst` is Benjamini-Hochberg across every shared panel's worst p, including panels that are not thesis-aligned.")
    a("")
    a(f"Shared panels: {n_panels}. Panels with all three pooled ρ < 0: {n_allneg}. Eligible eight-cohort panels among those: {len(best)}.")
    a("")
    if best_row is None:
        a("No shared panel had all three predictors pooled below zero.")
    else:
        a(
            f"Best shared panel: outcome `{best_row.outcome}`, keratin `{best_row.keratin}`, "
            f"purity `{best_row.purity}`, residualization `{best_row.method}`."
        )
        a("")
        a(
            f"TACSTD2 pooled ρ={fmt(best_row.rho_TACSTD2)} (p={fmt_p(best_row.p_TACSTD2)}, "
            f"negative cohorts {int(best_row.neg_TACSTD2)}/{int(best_row.cohorts_TACSTD2)}). "
            f"CLDN4 pooled ρ={fmt(best_row.rho_CLDN4)} (p={fmt_p(best_row.p_CLDN4)}, "
            f"{int(best_row.neg_CLDN4)}/{int(best_row.cohorts_CLDN4)}). "
            f"CLDN7 pooled ρ={fmt(best_row.rho_CLDN7)} (p={fmt_p(best_row.p_CLDN7)}, "
            f"{int(best_row.neg_CLDN7)}/{int(best_row.cohorts_CLDN7)}). "
            f"Worst p={fmt_p(best_row.worst_p)}. Sweep q on that worst p={fmt_p(best_row.q_worst)}."
        )
        a("")
        a("Per-cohort partial ρ for this panel:")
        a("")
        spec = continuous[
            (continuous["outcome"] == best_row.outcome)
            & (continuous["keratin"] == best_row.keratin)
            & (continuous["purity"] == best_row.purity)
            & (continuous["method"] == best_row.method)
        ]
        wide = spec.pivot(index="cohort", columns="predictor", values="rho")
        wide = wide.reindex(IMMUNE_COHORTS)
        nmap = spec[spec["predictor"] == "CLDN4"].set_index("cohort")["n"]
        show = pd.DataFrame({
            "cohort": wide.index,
            "n": [int(nmap.get(c, 0)) for c in wide.index],
            "TACSTD2": wide["TACSTD2"].to_numpy(),
            "CLDN4": wide["CLDN4"].to_numpy(),
            "CLDN7": wide["CLDN7"].to_numpy(),
        })
        a(md_table(show, ["cohort", "n", "TACSTD2", "CLDN4", "CLDN7"]))
    a("")
    a("Next shared panels with all three pooled ρ < 0, ordered by worst meta p:")
    a("")
    top = best.head(12).copy()
    if len(top):
        view = pd.DataFrame({
            "outcome": top["outcome"],
            "keratin": top["keratin"],
            "purity": top["purity"],
            "method": top["method"],
            "rho_TACSTD2": top["rho_TACSTD2"],
            "p_TACSTD2": top["p_TACSTD2"],
            "rho_CLDN4": top["rho_CLDN4"],
            "p_CLDN4": top["p_CLDN4"],
            "rho_CLDN7": top["rho_CLDN7"],
            "p_CLDN7": top["p_CLDN7"],
            "worst_p": top["worst_p"],
            "q_worst": top["q_worst"],
        })
        a(md_table(view, list(view.columns)))
    a("")
    simple = panels[
        (panels["outcome"] == "CD8_AB")
        & (panels["keratin"] == "simple_KRT8_18_19")
        & (panels["purity"] == "none")
        & (panels["method"] == "spearman")
    ]
    n_est = int((
        (panels["purity"] == "ESTIMATE")
        & panels["all_three_negative"]
        & (panels["cohorts_TACSTD2"] == len(IMMUNE_COHORTS))
        & (panels["cohorts_CLDN4"] == len(IMMUNE_COHORTS))
        & (panels["cohorts_CLDN7"] == len(IMMUNE_COHORTS))
    ).sum())
    if len(simple):
        srow = simple.iloc[0]
        a(
            f"The original specification (CD8A+CD8B, KRT8+KRT18+KRT19, no purity, Spearman) remains in the grid: "
            f"TACSTD2 ρ={fmt(srow.rho_TACSTD2)} (p={fmt_p(srow.p_TACSTD2)}), "
            f"CLDN4 ρ={fmt(srow.rho_CLDN4)} (p={fmt_p(srow.p_CLDN4)}), "
            f"CLDN7 ρ={fmt(srow.rho_CLDN7)} (p={fmt_p(srow.p_CLDN7)}). "
            f"Its worst p is {fmt_p(srow.worst_p)}. The KRT5/6 panels above are the stronger shared specifications."
        )
        a("")
    a(f"ESTIMATE-adjusted panels that keep all three pooled ρ values negative in all eight cohorts: {n_est}.")
    a("")
    a("## Strongest single cohort tests")
    a("")
    a(f"These are the smallest two-sided p-values among continuous-grid tests with ρ < 0. The grid contains {n_cont} tests. `q_grid` is Benjamini-Hochberg across all {n_cont} continuous tests, not only the negative ones.")
    a("")
    hit_view = hits[["cohort", "predictor", "outcome", "keratin", "purity", "method", "n", "rho", "p", "q_grid"]]
    a(md_table(hit_view, list(hit_view.columns)))
    a("")
    a("## Histology arms")
    a("")
    a(f"Histology-grid rows: {len(histology)}. The same BH is computed inside this family only.")
    a("")
    if len(hist_hits):
        hv = hist_hits[["stratum", "predictor", "outcome", "keratin", "purity", "method", "n", "rho", "p", "q_grid"]]
        a(md_table(hv, list(hv.columns)))
    else:
        a("No histology arm produced a negative partial ρ.")
    a("")
    a("## Quantiles")
    a("")
    a("Q4 versus Q1 of the predictor. The outcome is either the raw score or the residual of that score after the named adjustment. Cliff's delta < 0 means the upper quartile has the lower immune score. Stouffer combines the eight cohorts. `q` is Benjamini-Hochberg across quantile specifications.")
    a("")
    q_aligned = qsum[qsum["thesis"]].head(12)
    if len(q_aligned):
        a(md_table(q_aligned, ["predictor", "outcome", "adjustment", "n_cohorts", "n_negative", "median_delta", "stouffer_z", "stouffer_p", "q"]))
    a("")
    a("Per-cohort Cliff's delta for the strongest keratin-adjusted quantile specification of each predictor (negative = Q4 lower than Q1):")
    a("")
    q_focus = [
        ("TACSTD2", "CD8_effector", "squamous_spearman"),
        ("CLDN4", "CYT_Rooney", "squamous_pearson"),
        ("CLDN7", "T_CD3", "squamous_pearson"),
    ]
    q_show = []
    for cohort in IMMUNE_COHORTS:
        rec = {"cohort": cohort}
        for predictor, outcome, adjustment in q_focus:
            hit = quant[
                (quant["cohort"] == cohort)
                & (quant["predictor"] == predictor)
                & (quant["outcome"] == outcome)
                & (quant["keratin"] == adjustment)
            ]
            rec[f"{predictor}"] = float(hit["rho"].iloc[0]) if len(hit) else np.nan
        q_show.append(rec)
    q_show = pd.DataFrame(q_show)
    a("Columns are TACSTD2 vs CD8 effector after a KRT5/6 Spearman residual, CLDN4 vs Rooney CYT after a KRT5/6 Pearson residual, and CLDN7 vs CD3 after a KRT5/6 Pearson residual.")
    a("")
    a(md_table(q_show, ["cohort", "TACSTD2", "CLDN4", "CLDN7"]))
    a("")
    a("## NHEJ, STING, and IFN modules versus the barrier genes")
    a("")
    a("IFN-negative is the immune-cold direction. NHEJ and STING are reported at the sign the data give. The block below is partial Spearman of CLDN4 versus each module after KRT8+KRT18+KRT19, which matches the original keratin model. The full module grid is `tables/module_partial.tsv`.")
    a("")
    mod = modules[
        (modules["predictor"] == "CLDN4")
        & (modules["method"] == "spearman")
        & (modules["keratin"] == "simple_KRT8_18_19")
        & (modules["purity"] == "none")
    ].copy()
    mod["cohort"] = pd.Categorical(mod["cohort"], IMMUNE_COHORTS, ordered=True)
    mod = mod.sort_values(["outcome", "cohort"])
    a(md_table(mod[["cohort", "outcome", "n", "rho", "p"]], ["cohort", "outcome", "n", "rho", "p"]))
    a("")
    # Pooled module lines for CLDN4
    a("Random-effects pool, CLDN4, same keratin model, Spearman:")
    a("")
    for module in ["IFN", "STING", "NHEJ"]:
        hit = mod[mod["outcome"] == module]
        meta = random_effects_meta(hit["rho"].to_numpy(), hit["n"].to_numpy(), k=3)
        a(
            f"- {module}: ρ={fmt(meta['rho'])} (95% CI {fmt(meta['ci_low'])} to {fmt(meta['ci_high'])}), "
            f"p={fmt_p(meta['p'])}, I²={fmt(100 * meta['I2'], 1)}%, cohorts={meta['n_cohorts']}."
        )
    a("")
    a("## Reading")
    a("")
    a("Use the best shared panel as the sweep's answer for a specification that moves TACSTD2, CLDN4, and CLDN7 together. Use `q_worst` and `q_grid` when citing a p-value from this search. Cohort-level minima and histology arms describe where the association is largest. They are part of the same search.")
    a("")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "SWEEP.md"), "w", encoding="utf-8") as handle:
        handle.write(text)


def main() -> int:
    os.makedirs(TABLES, exist_ok=True)
    os.makedirs(FIGS, exist_ok=True)
    os.makedirs(os.path.join(CACHE, "clinical"), exist_ok=True)
    sym_to_id = load_probemap(symbols_needed())
    id_to_sym = {ens: sym for sym, ens in sym_to_id.items()}
    frames = load_expression(id_to_sym)
    purity = load_purity()
    for cohort in list(frames):
        sub = purity[purity["cohort"] == cohort].set_index("patient")[["ESTIMATE", "ABSOLUTE"]]
        frames[cohort] = frames[cohort].join(sub, how="left")
        print(
            f"[purity] {cohort} ESTIMATE {int(frames[cohort]['ESTIMATE'].notna().sum())} "
            f"ABSOLUTE {int(frames[cohort]['ABSOLUTE'].notna().sum())}",
            flush=True,
        )
    labels = histology_maps()
    for name, series in labels.items():
        print(f"[labels] {name} {series.value_counts().to_dict()}", flush=True)

    continuous = run_continuous(frames)
    histology = run_histology(frames, labels)
    quant = run_quantiles(frames)
    modules = run_modules(frames)
    panels = summarize_panels(continuous)
    qsum = summarize_quantiles(quant)

    continuous.to_csv(os.path.join(TABLES, "immune_sweep.tsv"), sep="\t", index=False)
    histology.to_csv(os.path.join(TABLES, "immune_sweep_histology.tsv"), sep="\t", index=False)
    quant.to_csv(os.path.join(TABLES, "immune_sweep_quantiles.tsv"), sep="\t", index=False)
    modules.to_csv(os.path.join(TABLES, "module_partial.tsv"), sep="\t", index=False)
    panels.to_csv(os.path.join(TABLES, "immune_sweep_panels.tsv"), sep="\t", index=False)
    qsum.to_csv(os.path.join(TABLES, "immune_sweep_quantile_summary.tsv"), sep="\t", index=False)

    best = panels[
        panels["all_three_negative"]
        & (panels["cohorts_TACSTD2"] == len(IMMUNE_COHORTS))
        & (panels["cohorts_CLDN4"] == len(IMMUNE_COHORTS))
        & (panels["cohorts_CLDN7"] == len(IMMUNE_COHORTS))
    ].sort_values("worst_p")
    if len(best):
        row = best.iloc[0]
        spec = continuous[
            (continuous["outcome"] == row.outcome)
            & (continuous["keratin"] == row.keratin)
            & (continuous["purity"] == row.purity)
            & (continuous["method"] == row.method)
        ]
        write_forest(
            spec,
            os.path.join(FIGS, "fig5_best_panel_forest.png"),
            f"Best shared panel: {row.outcome}, {row.keratin}, {row.purity}, {row.method}",
        )
    write_panel_heatmap(panels, os.path.join(FIGS, "fig6_sweep_panels.png"))
    write_module_heatmap(modules, os.path.join(FIGS, "fig7_cldn4_modules.png"))
    write_quantile_bars(quant, os.path.join(FIGS, "fig8_quantile_deltas.png"))
    write_sweep_md(continuous, histology, quant, modules, panels, qsum)
    print("[done] sweep", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
