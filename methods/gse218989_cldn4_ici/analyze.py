#!/usr/bin/env python3
"""GSE218989 SMC-KAIST NSCLC ICI bulk TPM: CLDN4 + TACSTD2 (patient unit).

Public inputs only:
  - GEO series GSE218989 TPM matrix + series-matrix characteristics
  - Nature Communications Supplementary Data 8 (author-correction XLSX)
    for OS / Death / PFS time. Histology is not in either source.

Does not invent clinical labels. LUAD vs LUSC is not split.
"""
from __future__ import annotations

import gzip
import json
import math
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
CACHE = Path("/tmp/gse218989")
GMT = HERE / "SI_geneset.gmt"

TPM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/suppl/"
    "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/matrix/"
    "GSE218989_series_matrix.txt.gz"
)
# Author correction (2025): updated Supplementary Data 8
CLIN_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41467-025-58068-y/MediaObjects/"
    "41467_2025_58068_MOESM3_ESM.xlsx"
)

AYERS6 = ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG"]
MHCI = ["HLA-A", "HLA-B", "HLA-C"]
TJ = ["CLDN3", "CLDN4", "CLDN7", "OCLN", "TJP1", "F11R", "MARVELD2", "CRB3", "PARD3", "CGN"]


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    print(f"download {url} -> {dest}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as f:
        f.write(r.read())
    return dest


def load_tpm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    df.index = df.index.astype(str)
    df = df[~df.index.duplicated(keep="first")]
    return df.astype(float)


def load_geo_meta(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        lines = f.read().splitlines()
    titles = None
    cols: dict[str, list[str]] = {}
    for line in lines:
        if line.startswith("!Sample_title"):
            titles = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            cols["gsm"] = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip('"') for x in line.split("\t")[1:]]
            key = vals[0].split(":", 1)[0].strip()
            cols[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    meta = pd.DataFrame({"patient": titles, **cols})
    return meta


def load_gmt(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            out[p[0]] = p[2:]
    return out


def mean_z(logx: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in logx.index]
    if not present:
        return pd.Series(np.nan, index=logx.columns), present
    z = logx.loc[present].T.apply(lambda s: (s - s.mean()) / (s.std(ddof=0) or 1.0), axis=0)
    return z.mean(axis=1), present


def ssgsea_matrix(logx: pd.DataFrame, genes: list[str], alpha: float = 0.25) -> pd.Series:
    """Barbie/GSVA-style ssGSEA (sum of running ES); not the ESTIMATE R package."""
    present = [g for g in genes if g in logx.index]
    if not present:
        return pd.Series(np.nan, index=logx.columns)
    X = logx.to_numpy()
    gene_idx = {g: i for i, g in enumerate(logx.index)}
    mask = np.zeros(X.shape[0], dtype=bool)
    for g in present:
        mask[gene_idx[g]] = True
    n = X.shape[0]
    nh = int(mask.sum())
    scores = np.empty(X.shape[1], dtype=float)
    for j in range(X.shape[1]):
        x = X[:, j]
        # average ranks, 1..n
        order = np.argsort(x, kind="mergesort")
        ranks = np.empty(n, dtype=float)
        ranks[order] = np.arange(1, n + 1, dtype=float)
        # high expression first
        ord_desc = np.argsort(-ranks, kind="mergesort")
        in_set = mask[ord_desc]
        ranked = ranks[ord_desc] ** alpha
        hits = ranked * in_set
        hit_sum = hits.sum()
        if hit_sum == 0:
            scores[j] = np.nan
            continue
        p_hit = np.cumsum(hits) / hit_sum
        p_miss = np.cumsum(~in_set) / (n - nh)
        scores[j] = float((p_hit - p_miss).sum())
    s = pd.Series(scores, index=logx.columns)
    # sample-wise max-abs normalization (GSVA default for ssGSEA)
    m = s.abs().max()
    if m and np.isfinite(m) and m != 0:
        s = s / m
    return s


def spearman(a: pd.Series, b: pd.Series) -> dict:
    m = pd.concat([a, b], axis=1).dropna()
    if len(m) < 5:
        return {"n": int(len(m)), "rho": np.nan, "p": np.nan}
    r, p = stats.spearmanr(m.iloc[:, 0], m.iloc[:, 1])
    return {"n": int(len(m)), "rho": float(r), "p": float(p)}


def mwu_auc(score: pd.Series, y: pd.Series) -> dict:
    """y=1 responder. AUC = P(score_R > score_NR) with ties 0.5."""
    m = pd.concat([score, y], axis=1).dropna()
    m.columns = ["s", "y"]
    a = m.loc[m.y == 1, "s"]
    b = m.loc[m.y == 0, "s"]
    if len(a) < 2 or len(b) < 2:
        return {"n_pos": int(len(a)), "n_neg": int(len(b)), "median_pos": np.nan,
                "median_neg": np.nan, "U": np.nan, "p": np.nan, "auc": np.nan}
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    auc = float(U / (len(a) * len(b)))
    return {
        "n_pos": int(len(a)),
        "n_neg": int(len(b)),
        "median_pos": float(a.median()),
        "median_neg": float(b.median()),
        "U": float(U),
        "p": float(p),
        "auc": auc,
    }


def cox_hr(score: pd.Series, time: pd.Series, event: pd.Series, label: str) -> dict:
    d = pd.DataFrame({"score": score, "time": time, "event": event}).dropna()
    d = d[d["time"] > 0]
    out = {"n": int(len(d)), "n_event": int(d["event"].sum()), "hr": np.nan,
           "hr_lo": np.nan, "hr_hi": np.nan, "p": np.nan, "label": label}
    if len(d) < 10 or d["event"].sum() < 5 or d["score"].std() == 0:
        return out
    z = (d["score"] - d["score"].mean()) / d["score"].std(ddof=0)
    dd = pd.DataFrame({"z": z, "time": d["time"], "event": d["event"]})
    cph = CoxPHFitter()
    cph.fit(dd, duration_col="time", event_col="event")
    row = cph.summary.loc["z"]
    out.update({
        "hr": float(row["exp(coef)"]),
        "hr_lo": float(row["exp(coef) lower 95%"]),
        "hr_hi": float(row["exp(coef) upper 95%"]),
        "p": float(row["p"]),
    })
    return out


def logrank_median(score: pd.Series, time: pd.Series, event: pd.Series) -> dict:
    d = pd.DataFrame({"score": score, "time": time, "event": event}).dropna()
    d = d[d["time"] > 0]
    med = float(d["score"].median())
    hi = d["score"] >= med
    res = logrank_test(d.loc[hi, "time"], d.loc[~hi, "time"],
                       d.loc[hi, "event"], d.loc[~hi, "event"])
    return {
        "n_hi": int(hi.sum()),
        "n_lo": int((~hi).sum()),
        "median_cut": med,
        "p": float(res.p_value),
        "test_stat": float(res.test_statistic),
    }


def fmt_p(p: float) -> str:
    if p is None or (isinstance(p, float) and (math.isnan(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def fmt_num(x: float, nd: int = 3) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x))):
        return "NA"
    return f"{x:.{nd}f}"


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    tpm_path = download(TPM_URL, CACHE / "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz")
    mtx_path = download(MATRIX_URL, CACHE / "GSE218989_series_matrix.txt.gz")
    clin_path = download(CLIN_URL, CACHE / "supp_data8_correction.xlsx")

    tpm = load_tpm(tpm_path)
    geo = load_geo_meta(mtx_path).set_index("patient")
    clin = pd.read_excel(clin_path).set_index("PatientID")
    clin.index = clin.index.astype(str)

    # patient unit: one column per patient
    patients = list(tpm.columns)
    assert len(patients) == len(set(patients))
    assert set(patients) == set(geo.index)

    logx = np.log2(tpm + 1.0)

    gene_presence = {
        "n_genes_matrix": int(tpm.shape[0]),
        "n_patients": int(tpm.shape[1]),
        "CLDN4": "CLDN4" in tpm.index,
        "TACSTD2": "TACSTD2" in tpm.index,
        "CD8A": "CD8A" in tpm.index,
        "IFNG": "IFNG" in tpm.index,
    }

    gmt = load_gmt(GMT)
    stromal_genes = gmt.get("StromalSignature", [])
    immune_genes = gmt.get("ImmuneSignature", [])

    cldn4 = logx.loc["CLDN4"]
    tacstd2 = logx.loc["TACSTD2"]
    cd8a = logx.loc["CD8A"]
    ifng = logx.loc["IFNG"]
    ayers, ayers_used = mean_z(logx, AYERS6)
    mhci, mhci_used = mean_z(logx, MHCI)
    tj, tj_used = mean_z(logx, [g for g in TJ if g != "CLDN4"])
    immune_z, immune_used = mean_z(logx, immune_genes)
    stromal_z, stromal_used = mean_z(logx, stromal_genes)
    immune_ss = ssgsea_matrix(logx, immune_genes)
    stromal_ss = ssgsea_matrix(logx, stromal_genes)

    geo_r = geo.loc[patients, "treatment outcome"].map(
        {"Responder": 1, "Non-responder": 0}
    ).astype(int)
    auth_r = clin.loc[patients, "Responder"].astype(int)
    os_days = clin.loc[patients, "Overall survival (days)"].astype(float)
    death = clin.loc[patients, "Death"].astype(int)
    pfs_days = clin.loc[patients, "Progression-free survival (days)"].astype(float)
    # Reconstructed PFS event: not deposited. Progression if PFS < OS; death counts.
    pfs_event_recon = ((pfs_days < os_days) | (death == 1)).astype(int)

    per = pd.DataFrame({
        "patient": patients,
        "gsm": geo.loc[patients, "gsm"].values,
        "geo_response": geo.loc[patients, "treatment outcome"].values,
        "geo_responder": geo_r.values,
        "supp8_responder": auth_r.values,
        "response_source_disagree": (geo_r.values != auth_r.values),
        "os_days": os_days.values,
        "death": death.values,
        "pfs_days": pfs_days.values,
        "pfs_event_reconstructed": pfs_event_recon.values,
        "histology_public": "not_deposited",
        "CLDN4_log2tpm1": cldn4.values,
        "TACSTD2_log2tpm1": tacstd2.values,
        "CD8A_log2tpm1": cd8a.values,
        "IFNG_log2tpm1": ifng.values,
        "IFNG_Ayers6_meanz": ayers.values,
        "MHCI_meanz": mhci.values,
        "TJ_noCLDN4_meanz": tj.values,
        "ESTIMATE_Immune_meanz": immune_z.values,
        "ESTIMATE_Stromal_meanz": stromal_z.values,
        "ESTIMATE_Immune_ssgsea": immune_ss.values,
        "ESTIMATE_Stromal_ssgsea": stromal_ss.values,
    }).set_index("patient")
    per.to_csv(TABLES / "per_patient.tsv", sep="\t")

    tests: list[dict] = []

    def add(row: dict) -> None:
        tests.append(row)

    # Immune-axis Spearman (pre-specified)
    axes = [
        ("CD8A", cd8a),
        ("IFNG", ifng),
        ("IFNG_Ayers6", ayers),
        ("MHCI", mhci),
        ("ESTIMATE_Immune_meanz", immune_z),
        ("ESTIMATE_Stromal_meanz", stromal_z),
        ("ESTIMATE_Immune_ssgsea", immune_ss),
        ("TJ_noCLDN4", tj),
    ]
    for gene_name, gene_s in [("CLDN4", cldn4), ("TACSTD2", tacstd2)]:
        for ax_name, ax_s in axes:
            sp = spearman(gene_s, ax_s)
            add({"family": "spearman_immune", "feature": gene_name, "axis": ax_name, **sp})
    add({"family": "spearman_immune", "feature": "CLDN4", "axis": "TACSTD2",
         **spearman(cldn4, tacstd2)})

    # Response: GEO primary; Supp8 sensitivity
    for label_name, y in [("GEO_response", geo_r), ("Supp8_response", auth_r)]:
        for feat, s in [("CLDN4", cldn4), ("TACSTD2", tacstd2), ("CD8A", cd8a),
                        ("IFNG", ifng), ("IFNG_Ayers6", ayers), ("MHCI", mhci),
                        ("ESTIMATE_Immune_meanz", immune_z)]:
            mw = mwu_auc(s, y)
            add({"family": "response_mwu", "label": label_name, "feature": feat, **mw})

    # Survival
    for feat, s in [("CLDN4", cldn4), ("TACSTD2", tacstd2), ("CD8A", cd8a),
                    ("IFNG_Ayers6", ayers), ("ESTIMATE_Immune_meanz", immune_z)]:
        add({"family": "cox_os", "feature": feat, **cox_hr(s, os_days, death, "OS")})
        add({"family": "logrank_os_median", "feature": feat, **logrank_median(s, os_days, death)})
        add({"family": "cox_pfs_reconstructed", "feature": feat,
             **cox_hr(s, pfs_days, pfs_event_recon, "PFS_reconstructed")})
        add({"family": "spearman_pfs_time", "feature": feat, **spearman(s, pfs_days)})
        add({"family": "spearman_os_time", "feature": feat, **spearman(s, os_days)})

    stats_df = pd.DataFrame(tests)
    # BH within immune Spearman for the two genes x core axes
    core_axes = {"CD8A", "IFNG", "IFNG_Ayers6", "MHCI", "ESTIMATE_Immune_meanz"}
    mask = (
        (stats_df["family"] == "spearman_immune")
        & stats_df["feature"].isin(["CLDN4", "TACSTD2"])
        & stats_df["axis"].isin(core_axes)
    )
    pvals = stats_df.loc[mask, "p"].astype(float).to_numpy()
    if len(pvals):
        order = np.argsort(pvals)
        q = np.empty_like(pvals)
        n = len(pvals)
        prev = 1.0
        for rank, idx in enumerate(order[::-1], start=0):
            i = order[::-1][rank]
            adj = pvals[i] * n / (n - rank)
            prev = min(prev, adj)
            q[i] = prev
        stats_df.loc[mask, "q_bh"] = q
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    def pick(family, **kw) -> pd.Series:
        m = stats_df["family"] == family
        for k, v in kw.items():
            m = m & (stats_df[k] == v)
        rows = stats_df.loc[m]
        if rows.empty:
            return pd.Series(dtype=float)
        return rows.iloc[0]

    c4_cd8 = pick("spearman_immune", feature="CLDN4", axis="CD8A")
    t2_cd8 = pick("spearman_immune", feature="TACSTD2", axis="CD8A")
    c4_ifn = pick("spearman_immune", feature="CLDN4", axis="IFNG_Ayers6")
    t2_ifn = pick("spearman_immune", feature="TACSTD2", axis="IFNG_Ayers6")
    c4_mhci = pick("spearman_immune", feature="CLDN4", axis="MHCI")
    t2_mhci = pick("spearman_immune", feature="TACSTD2", axis="MHCI")
    c4_imm = pick("spearman_immune", feature="CLDN4", axis="ESTIMATE_Immune_meanz")
    t2_imm = pick("spearman_immune", feature="TACSTD2", axis="ESTIMATE_Immune_meanz")
    c4_tj = pick("spearman_immune", feature="CLDN4", axis="TJ_noCLDN4")
    c4_t2 = pick("spearman_immune", feature="CLDN4", axis="TACSTD2")
    c4_geo = pick("response_mwu", label="GEO_response", feature="CLDN4")
    t2_geo = pick("response_mwu", label="GEO_response", feature="TACSTD2")
    cd8_geo = pick("response_mwu", label="GEO_response", feature="CD8A")
    ifn_geo = pick("response_mwu", label="GEO_response", feature="IFNG_Ayers6")
    c4_os = pick("cox_os", feature="CLDN4")
    t2_os = pick("cox_os", feature="TACSTD2")
    c4_pfs = pick("cox_pfs_reconstructed", feature="CLDN4")
    t2_pfs = pick("cox_pfs_reconstructed", feature="TACSTD2")
    c4_lr = pick("logrank_os_median", feature="CLDN4")
    t2_lr = pick("logrank_os_median", feature="TACSTD2")

    one = pd.DataFrame([{
        "dataset": "GSE218989_SMC_KAIST",
        "unit": "patient",
        "n": int(len(patients)),
        "n_genes": int(tpm.shape[0]),
        "CLDN4_on_matrix": True,
        "TACSTD2_on_matrix": True,
        "n_GEO_R": int((geo_r == 1).sum()),
        "n_GEO_NR": int((geo_r == 0).sum()),
        "n_Supp8_R": int((auth_r == 1).sum()),
        "n_Supp8_NR": int((auth_r == 0).sum()),
        "n_response_label_disagree": int((geo_r != auth_r).sum()),
        "histology_public": "not_deposited",
        "PFS_event_public": "not_deposited",
        "OS_event_public": "SuppData8_Death",
        "CLDN4_vs_GEO_R_median": c4_geo.get("median_pos", np.nan),
        "CLDN4_vs_GEO_NR_median": c4_geo.get("median_neg", np.nan),
        "CLDN4_vs_GEO_MWU_p": c4_geo.get("p", np.nan),
        "CLDN4_vs_GEO_AUC": c4_geo.get("auc", np.nan),
        "TACSTD2_vs_GEO_MWU_p": t2_geo.get("p", np.nan),
        "TACSTD2_vs_GEO_AUC": t2_geo.get("auc", np.nan),
        "CLDN4_vs_CD8A_rho": c4_cd8.get("rho", np.nan),
        "CLDN4_vs_CD8A_p": c4_cd8.get("p", np.nan),
        "TACSTD2_vs_CD8A_rho": t2_cd8.get("rho", np.nan),
        "TACSTD2_vs_CD8A_p": t2_cd8.get("p", np.nan),
        "CLDN4_vs_Ayers6_rho": c4_ifn.get("rho", np.nan),
        "CLDN4_vs_Ayers6_p": c4_ifn.get("p", np.nan),
        "TACSTD2_vs_Ayers6_rho": t2_ifn.get("rho", np.nan),
        "TACSTD2_vs_Ayers6_p": t2_ifn.get("p", np.nan),
        "CLDN4_vs_MHCI_rho": c4_mhci.get("rho", np.nan),
        "CLDN4_vs_MHCI_p": c4_mhci.get("p", np.nan),
        "TACSTD2_vs_MHCI_rho": t2_mhci.get("rho", np.nan),
        "TACSTD2_vs_MHCI_p": t2_mhci.get("p", np.nan),
        "CLDN4_vs_Immune_rho": c4_imm.get("rho", np.nan),
        "CLDN4_vs_Immune_p": c4_imm.get("p", np.nan),
        "TACSTD2_vs_Immune_rho": t2_imm.get("rho", np.nan),
        "TACSTD2_vs_Immune_p": t2_imm.get("p", np.nan),
        "CLDN4_vs_TJ_rho": c4_tj.get("rho", np.nan),
        "CLDN4_vs_TJ_p": c4_tj.get("p", np.nan),
        "CLDN4_vs_TACSTD2_rho": c4_t2.get("rho", np.nan),
        "CLDN4_vs_TACSTD2_p": c4_t2.get("p", np.nan),
        "CLDN4_OS_HR_perSD": c4_os.get("hr", np.nan),
        "CLDN4_OS_HR_lo": c4_os.get("hr_lo", np.nan),
        "CLDN4_OS_HR_hi": c4_os.get("hr_hi", np.nan),
        "CLDN4_OS_Cox_p": c4_os.get("p", np.nan),
        "CLDN4_OS_n_event": c4_os.get("n_event", np.nan),
        "TACSTD2_OS_HR_perSD": t2_os.get("hr", np.nan),
        "TACSTD2_OS_HR_lo": t2_os.get("hr_lo", np.nan),
        "TACSTD2_OS_HR_hi": t2_os.get("hr_hi", np.nan),
        "TACSTD2_OS_Cox_p": t2_os.get("p", np.nan),
        "TACSTD2_OS_n_event": t2_os.get("n_event", np.nan),
        "CLDN4_PFS_recon_HR_perSD": c4_pfs.get("hr", np.nan),
        "CLDN4_PFS_recon_p": c4_pfs.get("p", np.nan),
        "TACSTD2_PFS_recon_HR_perSD": t2_pfs.get("hr", np.nan),
        "TACSTD2_PFS_recon_p": t2_pfs.get("p", np.nan),
        "CD8A_vs_GEO_MWU_p": cd8_geo.get("p", np.nan),
        "Ayers6_vs_GEO_MWU_p": ifn_geo.get("p", np.nan),
    }])
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    inventory = pd.DataFrame([
        {"field": "GEO_treatment_outcome", "present": True,
         "n": 355, "notes": "Responder 168 / Non-responder 187; series matrix"},
        {"field": "GEO_treatment", "present": True, "n": 355,
         "notes": "PD-1/PD-L1 inhibitor for all 355"},
        {"field": "GEO_histology", "present": False, "n": 0,
         "notes": "not a series-matrix characteristic; LUAD vs LUSC not split"},
        {"field": "GEO_PFS", "present": False, "n": 0, "notes": "not on GEO"},
        {"field": "GEO_OS", "present": False, "n": 0, "notes": "not on GEO"},
        {"field": "SuppData8_Responder", "present": True, "n": 355,
         "notes": "497-row published table; 355 overlap TPM; 3 GEO disagreements (author correction)"},
        {"field": "SuppData8_OS_days_Death", "present": True, "n": 355,
         "notes": "Overall survival (days) + Death"},
        {"field": "SuppData8_PFS_days", "present": True, "n": 355,
         "notes": "PFS time only; no PFS event column deposited"},
        {"field": "SuppData8_histology", "present": False, "n": 0,
         "notes": "Clinical_table_v230613 has no histology column"},
    ])
    inventory.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    gene_cov = pd.DataFrame([
        {"set": "Ayers6_IFNG", "n_listed": len(AYERS6), "n_present": len(ayers_used),
         "genes": ",".join(ayers_used)},
        {"set": "MHCI", "n_listed": len(MHCI), "n_present": len(mhci_used),
         "genes": ",".join(mhci_used)},
        {"set": "TJ_noCLDN4", "n_listed": len([g for g in TJ if g != "CLDN4"]),
         "n_present": len(tj_used), "genes": ",".join(tj_used)},
        {"set": "ESTIMATE_Immune", "n_listed": len(immune_genes),
         "n_present": len(immune_used), "genes": ""},
        {"set": "ESTIMATE_Stromal", "n_listed": len(stromal_genes),
         "n_present": len(stromal_used), "genes": ""},
    ])
    gene_cov.to_csv(TABLES / "signature_coverage.tsv", sep="\t", index=False)

    summary = {
        "accession": "GSE218989",
        "matrix": "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz",
        "n_patients": int(len(patients)),
        "n_genes": int(tpm.shape[0]),
        "CLDN4_present": True,
        "TACSTD2_present": True,
        "histology": "not_deposited",
        "geo_R_NR": [int((geo_r == 1).sum()), int((geo_r == 0).sum())],
        "disagree_ids": per.index[per["response_source_disagree"]].tolist(),
        "signature_coverage": gene_cov.to_dict(orient="records"),
        "one_row": one.iloc[0].to_dict(),
        "gene_presence": gene_presence,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    # ---------- figures ----------
    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 160,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    # Fig 1: response boxplots
    fig, axes = plt.subplots(1, 4, figsize=(11.2, 3.6))
    pairs = [
        (cldn4, "CLDN4 log2(TPM+1)"),
        (tacstd2, "TACSTD2 log2(TPM+1)"),
        (cd8a, "CD8A log2(TPM+1)"),
        (ayers, "IFN-γ Ayers6 mean-z"),
    ]
    rng = np.random.default_rng(1)
    for ax, (s, title) in zip(axes, pairs):
        mw = mwu_auc(s, geo_r)
        for i, (lab, yv) in enumerate([("NR", 0), ("R", 1)]):
            v = s[geo_r == yv]
            ax.boxplot(v, positions=[i], widths=0.55, showfliers=False,
                       medianprops={"color": "black"})
            ax.scatter(np.full(len(v), i) + rng.uniform(-0.12, 0.12, len(v)),
                       v, s=8, alpha=0.35, c="#4C78A8" if yv == 0 else "#F58518")
        ax.set_xticks([0, 1], [f"NR\nn={mw['n_neg']}", f"R\nn={mw['n_pos']}"])
        ax.set_title(title, fontsize=9)
        ax.set_ylabel("score")
        ax.text(0.5, 0.98, f"MWU p={fmt_p(mw['p'])}\nAUC={fmt_num(mw['auc'])}",
                transform=ax.transAxes, ha="center", va="top", fontsize=8)
    fig.suptitle("GSE218989 GEO ICI response (patient unit, n=355)", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_response_boxplots.png", bbox_inches="tight")
    plt.close(fig)

    # Fig 2: CLDN4 vs immune
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.4))
    scat = [
        (cd8a, "CD8A log2(TPM+1)", c4_cd8),
        (ayers, "IFN-γ Ayers6 mean-z", c4_ifn),
        (mhci, "MHC-I (HLA-A/B/C) mean-z", c4_mhci),
        (immune_z, "ESTIMATE Immune mean-z", c4_imm),
    ]
    for ax, (y, ylab, sp) in zip(axes.ravel(), scat):
        ax.scatter(cldn4, y, s=10, alpha=0.4, c="#4C78A8", linewidths=0)
        ax.set_xlabel("CLDN4 log2(TPM+1)")
        ax.set_ylabel(ylab)
        ax.set_title(f"ρ={fmt_num(sp.get('rho', np.nan))}  p={fmt_p(sp.get('p', np.nan))}  n={int(sp.get('n', 0))}",
                     fontsize=9)
    fig.suptitle("CLDN4 vs immune axes — GSE218989 n=355", y=1.01)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_immune_scatter.png", bbox_inches="tight")
    plt.close(fig)

    # Fig 3: TACSTD2 vs immune
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.4))
    scat = [
        (cd8a, "CD8A log2(TPM+1)", t2_cd8),
        (ayers, "IFN-γ Ayers6 mean-z", t2_ifn),
        (mhci, "MHC-I (HLA-A/B/C) mean-z", t2_mhci),
        (immune_z, "ESTIMATE Immune mean-z", t2_imm),
    ]
    for ax, (y, ylab, sp) in zip(axes.ravel(), scat):
        ax.scatter(tacstd2, y, s=10, alpha=0.4, c="#54A24B", linewidths=0)
        ax.set_xlabel("TACSTD2 log2(TPM+1)")
        ax.set_ylabel(ylab)
        ax.set_title(f"ρ={fmt_num(sp.get('rho', np.nan))}  p={fmt_p(sp.get('p', np.nan))}  n={int(sp.get('n', 0))}",
                     fontsize=9)
    fig.suptitle("TACSTD2 (TROP2) vs immune axes — GSE218989 n=355", y=1.01)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_tacstd2_immune_scatter.png", bbox_inches="tight")
    plt.close(fig)

    # Fig 4: KM OS
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    for ax, s, name, lr in [(axes[0], cldn4, "CLDN4", c4_lr),
                            (axes[1], tacstd2, "TACSTD2", t2_lr)]:
        d = pd.DataFrame({"s": s, "t": os_days, "e": death}).dropna()
        med = d["s"].median()
        km = KaplanMeierFitter()
        for lab, msk, c in [(f"low (<median) n={(d.s < med).sum()}", d.s < med, "#4C78A8"),
                            (f"high (≥median) n={(d.s >= med).sum()}", d.s >= med, "#E45756")]:
            km.fit(d.loc[msk, "t"], d.loc[msk, "e"], label=lab)
            km.plot_survival_function(ax=ax, ci_show=False, color=c)
        ax.set_xlabel("OS (days)")
        ax.set_ylabel("Survival")
        ax.set_title(f"{name} median split\nlog-rank p={fmt_p(lr.get('p', np.nan))}", fontsize=10)
        ax.legend(fontsize=7, loc="upper right")
    fig.suptitle("OS from Supplementary Data 8 Death (n=355)", y=1.03)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_km_os.png", bbox_inches="tight")
    plt.close(fig)

    # Fig 5: correlation heatmap
    mat = pd.DataFrame({
        "CLDN4": cldn4, "TACSTD2": tacstd2, "CD8A": cd8a, "IFNG": ifng,
        "Ayers6": ayers, "MHC-I": mhci, "Immune_z": immune_z,
        "Stromal_z": stromal_z, "TJ_noCLDN4": tj,
    })
    cols = list(mat.columns)
    R = np.zeros((len(cols), len(cols)))
    P = np.ones((len(cols), len(cols)))
    for i, a in enumerate(cols):
        for j, b in enumerate(cols):
            r, p = stats.spearmanr(mat[a], mat[b])
            R[i, j] = r
            P[i, j] = p
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    im = ax.imshow(R, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(cols)), cols, rotation=45, ha="right")
    ax.set_yticks(range(len(cols)), cols)
    for i in range(len(cols)):
        for j in range(len(cols)):
            ax.text(j, i, f"{R[i, j]:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, label="Spearman ρ")
    ax.set_title("GSE218989 n=355 patient-level Spearman")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig5_correlation_heatmap.png", bbox_inches="tight")
    plt.close(fig)

    # Fig 6: forest of immune Spearman
    rows = []
    for feat, color in [("CLDN4", "#4C78A8"), ("TACSTD2", "#54A24B")]:
        for axn in ["CD8A", "IFNG", "IFNG_Ayers6", "MHCI", "ESTIMATE_Immune_meanz",
                    "ESTIMATE_Stromal_meanz", "TJ_noCLDN4"]:
            r = pick("spearman_immune", feature=feat, axis=axn)
            rows.append((f"{feat} vs {axn}", r.get("rho", np.nan), r.get("p", np.nan), color))
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    y = np.arange(len(rows))
    ax.axvline(0, color="0.5", lw=1)
    for i, (lab, rho, p, c) in enumerate(rows):
        ax.plot(rho, i, "o", color=c, ms=7)
        ax.plot([0, rho], [i, i], color=c, lw=2)
    ax.set_yticks(y, [r[0] for r in rows], fontsize=8)
    ax.set_xlabel("Spearman ρ (n=355)")
    ax.set_title("Immune / TJ axes")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURES / "fig6_spearman_forest.png", bbox_inches="tight")
    plt.close(fig)

    # Fig 7: CLDN4 vs TACSTD2 colored by response
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for yv, lab, c in [(0, "GEO NR", "#4C78A8"), (1, "GEO R", "#F58518")]:
        m = geo_r == yv
        ax.scatter(cldn4[m], tacstd2[m], s=12, alpha=0.45, c=c, label=lab, linewidths=0)
    ax.set_xlabel("CLDN4 log2(TPM+1)")
    ax.set_ylabel("TACSTD2 log2(TPM+1)")
    ax.legend(frameon=False)
    ax.set_title(f"CLDN4 vs TACSTD2  ρ={fmt_num(c4_t2.get('rho', np.nan))} p={fmt_p(c4_t2.get('p', np.nan))}")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig7_cldn4_vs_tacstd2.png", bbox_inches="tight")
    plt.close(fig)

    print("wrote", TABLES)
    print("wrote", FIGURES)
    print(one.T.to_string())


if __name__ == "__main__":
    main()
