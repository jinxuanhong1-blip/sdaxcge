"""
IMvigor210 (atezolizumab / anti-PD-L1, metastatic urothelial carcinoma, n=348)
Test: is TACSTD2 (TROP2) expression associated with WORSE benefit on PD-L1 blockade?

This is the closest *genuinely open* atezolizumab cohort with per-sample TACSTD2
expression + clinical outcomes. It is a PROXY for the controlled OAK/POPLAR data
(EGAS00001005013) used by Bessede et al. 2024 (Clin Cancer Res). Caveats (README):
  - single-arm (everyone gets atezo) -> tests the ASSOCIATION on-treatment, NOT the
    treatment-by-biomarker interaction (predictive-vs-prognostic) that Bessede made.
  - urothelial carcinoma, not NSCLC.

Source: IMvigor210CoreBiologies (Mariathasan et al., Nature 2018). Authoritative
CountDataSet `cds` (data/cds.RData) from mirror github.com/SiYangming/
IMvigor210CoreBiologies. Raw DESeq counts (Entrez rownames) + DESeq sizeFactors +
clinical, extracted to CSV in scripts/00_extract_imvigor210_cds.R.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

RAW = "results/hunt_oak/data/raw"
DERIVED = "results/hunt_oak/data/derived"
OUT = "results/hunt_oak/outputs"

counts = pd.read_csv(f"{RAW}/imvigor210_counts_selected.csv")
pData = pd.read_csv(f"{RAW}/imvigor210_pdata_full.csv")

df = pData.merge(counts, on="sample", how="inner")

# sanity: housekeeping ACTB must dwarf TACSTD2 (guards against a broken extraction)
assert df["ACTB"].median() > df["TACSTD2"].median() * 10, "count extraction sanity failed"

# DESeq-normalize (counts / sizeFactor) then log2(x+1)
for g in ["TACSTD2", "CD274", "EPCAM", "PTPRC", "CD8A"]:
    df[f"log2_{g}"] = np.log2(df[g] / df["sizeFactor"] + 1.0)

df["responder"] = df["binaryResponse"].map({"CR/PR": 1, "SD/PD": 0})
df["os_time"] = pd.to_numeric(df["os"], errors="coerce")
df["os_event"] = pd.to_numeric(df["censOS"], errors="coerce")  # 1 = death, 0 = censored

med = df["log2_TACSTD2"].median()
df["TACSTD2_group"] = np.where(df["log2_TACSTD2"] >= med, "high", "low")

keep = ["sample", "TACSTD2", "log2_TACSTD2", "log2_CD274", "log2_EPCAM",
        "log2_CD8A", "TACSTD2_group", "binaryResponse", "responder",
        "os_time", "os_event", "IC Level", "TC Level", "TCGA Subtype",
        "Immune phenotype"]
df[keep].to_csv(f"{DERIVED}/imvigor210_tacstd2_persample.csv", index=False)

res = {"dataset": "IMvigor210 (atezolizumab, mUC, single-arm, n=348)",
       "gene": "TACSTD2 (TROP2)", "n_total": int(len(df))}

# --- (A) response ----------------------------------------------------------
r = df.dropna(subset=["responder"])
hi = r.loc[r.responder == 1, "log2_TACSTD2"]
lo = r.loc[r.responder == 0, "log2_TACSTD2"]
u, p_mw = stats.mannwhitneyu(hi, lo, alternative="two-sided")
res["response_evaluable_n"] = int(len(r))
res["response_TACSTD2_median_responders"] = round(float(hi.median()), 4)
res["response_TACSTD2_median_nonresponders"] = round(float(lo.median()), 4)
res["response_mannwhitney_p"] = round(float(p_mw), 5)

ct = pd.crosstab(r["TACSTD2_group"], r["responder"]).reindex(index=["high", "low"]).fillna(0)
for c in (0, 1):
    if c not in ct.columns:
        ct[c] = 0
odds, p_fish = stats.fisher_exact(ct.loc[["high", "low"], [1, 0]].values)
res["dcb_rate_TACSTD2_high"] = round(float(r.loc[r.TACSTD2_group == "high", "responder"].mean()), 4)
res["dcb_rate_TACSTD2_low"] = round(float(r.loc[r.TACSTD2_group == "low", "responder"].mean()), 4)
res["dcb_fisher_or"] = round(float(odds), 4)
res["dcb_fisher_p"] = round(float(p_fish), 5)

# --- (B) overall survival --------------------------------------------------
s = df.dropna(subset=["os_time", "os_event"]).copy()
res["os_evaluable_n"] = int(len(s))
res["os_events"] = int(s["os_event"].sum())

s["z_TACSTD2"] = (s.log2_TACSTD2 - s.log2_TACSTD2.mean()) / s.log2_TACSTD2.std()
cph = CoxPHFitter().fit(s[["os_time", "os_event", "z_TACSTD2"]], "os_time", "os_event")
row = cph.summary.loc["z_TACSTD2"]
res["os_cox_perSD_HR"] = round(float(row["exp(coef)"]), 4)
res["os_cox_perSD_HR_95CI"] = [round(float(row["exp(coef) lower 95%"]), 4),
                               round(float(row["exp(coef) upper 95%"]), 4)]
res["os_cox_perSD_p"] = round(float(row["p"]), 5)

hi_m, lo_m = s[s.TACSTD2_group == "high"], s[s.TACSTD2_group == "low"]
lr = logrank_test(hi_m.os_time, lo_m.os_time, hi_m.os_event, lo_m.os_event)
res["os_logrank_p_medianSplit"] = round(float(lr.p_value), 5)
kmf = KaplanMeierFitter()
for name, g in (("high", hi_m), ("low", lo_m)):
    kmf.fit(g.os_time, g.os_event)
    res[f"os_median_TACSTD2_{name}"] = round(float(kmf.median_survival_time_), 3)

s["is_high"] = (s.TACSTD2_group == "high").astype(int)
cph2 = CoxPHFitter().fit(s[["os_time", "os_event", "is_high"]], "os_time", "os_event")
row2 = cph2.summary.loc["is_high"]
res["os_cox_highVSlow_HR"] = round(float(row2["exp(coef)"]), 4)
res["os_cox_highVSlow_HR_95CI"] = [round(float(row2["exp(coef) lower 95%"]), 4),
                                   round(float(row2["exp(coef) upper 95%"]), 4)]
res["os_cox_highVSlow_p"] = round(float(row2["p"]), 5)

rho, prho = stats.spearmanr(df.log2_TACSTD2, df.log2_CD274)
res["spearman_TACSTD2_vs_CD274_rho"] = round(float(rho), 4)
res["spearman_TACSTD2_vs_CD274_p"] = round(float(prho), 5)
rho2, prho2 = stats.spearmanr(df.log2_TACSTD2, df.log2_CD8A)
res["spearman_TACSTD2_vs_CD8A_rho"] = round(float(rho2), 4)
res["spearman_TACSTD2_vs_CD8A_p"] = round(float(prho2), 5)

with open(f"{OUT}/imvigor210_results.json", "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
