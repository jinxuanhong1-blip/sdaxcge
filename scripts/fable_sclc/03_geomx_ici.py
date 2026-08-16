"""GeoMx DSP spatial cohorts of ICI-treated ES-SCLC: TACSTD2 vs outcome.

Cohorts (both first-line chemo-immunotherapy, pre-treatment FFPE ROIs):
  - GSE261345: CANTABRICO trial (durvalumab + platinum/etoposide), 26 patients
  - GSE261348: IMfirst trial (atezolizumab + platinum/etoposide), 32 patients

CLDN4 is NOT in the GeoMx CTA panel; TACSTD2 is. Analyses:
  - ROI-level Spearman correlation of TACSTD2 with immune genes/signatures
  - Patient-level (median across ROIs): responders (CR/PR) vs
    non-responders (SD/PD) Mann-Whitney U
  - Cox proportional hazards for PFS and OS, per-cohort z-scored TACSTD2,
    pooled with cohort covariate; Kaplan-Meier by cohort-wise median split

Clinical fields are parsed from the GEO series matrices (per-GSM
characteristics include best RECIST response and dosing/progression/
follow-up dates).
"""

import gzip
import re

import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from common import (DATA_DIR, TABLES_DIR, FIGURES_DIR, TEFF_GENES,
                    SUBTYPE_TFS, ensure_dirs, bh_fdr, zscore)

COHORTS = {
    "CANTABRICO_GSE261345": {
        "counts": "GSE261345_CANTABRICO_DSP_normalizedcounts.xlsx",
        "series": "GSE261345_series_matrix.txt.gz",
    },
    "IMfirst_GSE261348": {
        "counts": "GSE261348_IMfirst_DSP_normalizedcounts.xlsx",
        "series": "GSE261348_series_matrix.txt.gz",
    },
}

IMMUNE_GENES = ["CD8A", "GZMA", "GZMB", "PRF1", "IFNG", "CXCL9", "CD274",
                "STAT1", "CD68", "CD163"]
RESPONDER = {"Complete response", "Partial response"}
NONRESPONDER = {"Stable disease", "Progressive disease"}


def parse_series_matrix(path):
    """Return per-GSM clinical dataframe keyed by SegmentDisplayName.

    Characteristics lines in GEO series matrices are positional and shift
    when a field is absent for a sample (e.g. no progression date), so each
    'label: value' pair is parsed per sample rather than per line.
    """
    char_lines, meta = [], {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key, vals = parts[0], [p.strip('"') for p in parts[1:]]
            if key == "!Sample_title":
                meta["title"] = vals
            elif key == "!Sample_geo_accession":
                meta["gsm"] = vals
            elif key == "!Sample_description":
                meta["segment"] = vals
            elif key == "!Sample_characteristics_ch1":
                char_lines.append(vals)
    records = []
    for i in range(len(meta["gsm"])):
        rec = {}
        for vals in char_lines:
            v = vals[i]
            if ": " in v:
                label, value = v.split(": ", 1)
                rec[label.lower().replace(" ", "_")] = value
        records.append(rec)
    df = pd.concat([pd.DataFrame(meta), pd.DataFrame(records)], axis=1)
    for col in [c for c in df.columns if c.startswith("date_")]:
        df[col] = pd.to_datetime(df[col], format="%m/%d/%Y", errors="coerce")
    return df


def load_counts(path):
    xl = pd.ExcelFile(path)
    tcm = xl.parse("TargetCountMatrix").set_index("TargetName")
    return np.log2(tcm + 1).T  # segments x genes, log2 normalized counts


def build_patient_table(clin):
    """One row per patient with response and PFS/OS times (months)."""
    p = clin.drop_duplicates("patient_id").set_index("patient_id").copy()
    day = pd.Timedelta(days=1)
    t0 = p["date_of_first_dose_of_treatment"]
    prog = p["date_of_disease_progression_or_death"]
    fup = p["date_of_last_follow-up"]
    pfs_event = (p["disease_progression_or_death"] == "Yes").astype(int)
    pfs_end = prog.where(pfs_event == 1, fup)
    os_event = (p["death"] == "Yes").astype(int)
    out = pd.DataFrame({
        "response": p["best_recist_response_to_treatment"],
        "pfs_months": (pfs_end - t0) / day / 30.44,
        "pfs_event": pfs_event,
        "os_months": (fup - t0) / day / 30.44,
        "os_event": os_event,
    })
    return out


def main():
    ensure_dirs()
    roi_frames, patient_frames, corr_rows = [], [], []

    for cohort, files in COHORTS.items():
        clin = parse_series_matrix(DATA_DIR / files["series"])
        expr = load_counts(DATA_DIR / files["counts"])
        clin = clin[clin["segment"].isin(expr.index)].copy()
        expr = expr.loc[clin["segment"]]
        expr.index = clin["gsm"].values
        clin = clin.set_index("gsm")
        print(f"[{cohort}] {expr.shape[0]} ROIs, "
              f"{clin['patient_id'].nunique()} patients")

        # ROI-level immune correlations of TACSTD2
        teff = zscore(expr[[g for g in TEFF_GENES if g in expr.columns]]).mean(axis=1)
        for cov_name, cov in ([("Teff_score", teff)] +
                              [(g, expr[g]) for g in IMMUNE_GENES
                               if g in expr.columns]):
            rho, p = stats.spearmanr(expr["TACSTD2"], cov)
            corr_rows.append({"cohort": cohort, "gene": "TACSTD2",
                              "covariate": cov_name, "n_roi": len(expr),
                              "spearman_rho": rho, "spearman_p": p})

        roi = pd.DataFrame({
            "cohort": cohort,
            "patient_id": clin["patient_id"],
            "TACSTD2_log2": expr["TACSTD2"],
            "Teff_score": teff,
            **{f"{g}_log2": expr[g] for g in SUBTYPE_TFS + ["YAP1"]},
        })
        roi_frames.append(roi)

        # patient level: median TACSTD2 over ROIs + outcomes
        med = roi.groupby("patient_id")[["TACSTD2_log2", "Teff_score"]].median()
        med["n_roi"] = roi.groupby("patient_id").size()
        pt = build_patient_table(clin).join(med)
        pt["cohort"] = cohort
        pt["TACSTD2_z"] = zscore(pt["TACSTD2_log2"])
        patient_frames.append(pt)

    corr = pd.DataFrame(corr_rows)
    corr["spearman_p_BH"] = bh_fdr(corr["spearman_p"])
    corr.to_csv(TABLES_DIR / "geomx_roi_immune_correlations.csv", index=False)

    pts = pd.concat(patient_frames)
    pts.index.name = "patient_id"
    pts.to_csv(TABLES_DIR / "geomx_patient_level.csv")

    # response comparison (pooled and per cohort)
    resp_rows = []
    for scope, d in [("pooled", pts)] + list(pts.groupby("cohort")):
        r = d[d["response"].isin(RESPONDER)]["TACSTD2_log2"].dropna()
        nr = d[d["response"].isin(NONRESPONDER)]["TACSTD2_log2"].dropna()
        if len(r) >= 3 and len(nr) >= 3:
            u, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
            resp_rows.append({"scope": scope, "n_resp": len(r),
                              "n_nonresp": len(nr),
                              "median_resp": r.median(),
                              "median_nonresp": nr.median(),
                              "mwu_U": u, "mwu_p": p})
    resp = pd.DataFrame(resp_rows)
    resp["mwu_p_BH"] = bh_fdr(resp["mwu_p"])
    resp.to_csv(TABLES_DIR / "geomx_response_stats.csv", index=False)
    print(resp)

    # Cox PH, pooled with cohort covariate
    cox_rows = []
    surv = pts.dropna(subset=["TACSTD2_z"]).copy()
    surv["cohort_flag"] = (surv["cohort"] == "IMfirst_GSE261348").astype(int)
    for endpoint in ["pfs", "os"]:
        d = surv.dropna(subset=[f"{endpoint}_months"]).copy()
        d = d[d[f"{endpoint}_months"] > 0]
        cph = CoxPHFitter()
        cph.fit(d[[f"{endpoint}_months", f"{endpoint}_event",
                   "TACSTD2_z", "cohort_flag"]],
                duration_col=f"{endpoint}_months",
                event_col=f"{endpoint}_event")
        s = cph.summary.loc["TACSTD2_z"]
        cox_rows.append({
            "endpoint": endpoint.upper(), "n": len(d),
            "n_events": int(d[f"{endpoint}_event"].sum()),
            "HR_per_SD_TACSTD2": np.exp(s["coef"]),
            "HR_CI95_low": np.exp(s["coef lower 95%"]),
            "HR_CI95_high": np.exp(s["coef upper 95%"]),
            "cox_p": s["p"],
        })
    cox = pd.DataFrame(cox_rows)
    cox.to_csv(TABLES_DIR / "geomx_cox_tacstd2.csv", index=False)
    print(cox)

    # figures: response boxplot + KM by cohort-wise median split (PFS)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    d = pts[pts["response"].isin(RESPONDER | NONRESPONDER)].copy()
    d["resp_group"] = np.where(d["response"].isin(RESPONDER),
                               "CR/PR", "SD/PD")
    sns.boxplot(data=d, x="resp_group", y="TACSTD2_log2", hue="cohort",
                showfliers=False, ax=axes[0])
    sns.stripplot(data=d, x="resp_group", y="TACSTD2_log2", hue="cohort",
                  dodge=True, palette='dark:k', size=3, legend=False, ax=axes[0])
    axes[0].set_title("TACSTD2 (patient median, GeoMx) vs best response")

    km_d = surv.dropna(subset=["pfs_months"]).copy()
    km_d = km_d[km_d["pfs_months"] > 0]
    km_d["group"] = km_d.groupby("cohort")["TACSTD2_log2"].transform(
        lambda x: np.where(x > x.median(), "TACSTD2 high", "TACSTD2 low"))
    kmf = KaplanMeierFitter()
    for grp, dd in km_d.groupby("group"):
        kmf.fit(dd["pfs_months"], dd["pfs_event"], label=f"{grp} (n={len(dd)})")
        kmf.plot_survival_function(ax=axes[1], ci_show=False)
    lr = logrank_test(
        km_d.loc[km_d["group"] == "TACSTD2 high", "pfs_months"],
        km_d.loc[km_d["group"] == "TACSTD2 low", "pfs_months"],
        km_d.loc[km_d["group"] == "TACSTD2 high", "pfs_event"],
        km_d.loc[km_d["group"] == "TACSTD2 low", "pfs_event"])
    axes[1].set_title(f"PFS by TACSTD2 median split (logrank p={lr.p_value:.3f})")
    axes[1].set_xlabel("months")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "geomx_tacstd2_response_pfs.png", dpi=170)
    plt.close(fig)

    with open(TABLES_DIR / "geomx_km_logrank.txt", "w") as fh:
        fh.write(f"PFS logrank TACSTD2 high vs low (cohort-wise median split): "
                 f"statistic={lr.test_statistic:.3f}, p={lr.p_value:.4f}\n")
        fh.write(f"n_high={int((km_d['group'] == 'TACSTD2 high').sum())}, "
                 f"n_low={int((km_d['group'] == 'TACSTD2 low').sum())}\n")

    pd.concat(roi_frames).to_csv(TABLES_DIR / "geomx_roi_level.csv", index=False)
    print("[done] GeoMx tables/figures written")


if __name__ == "__main__":
    main()
