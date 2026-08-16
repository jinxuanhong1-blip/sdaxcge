"""
GSE135222 (NSCLC, anti-PD-1/PD-L1, n=27) - secondary OPEN check in lung cancer.
Test: is TACSTD2 (TROP2) expression associated with shorter PFS on ICI?

Open TPM matrix + per-sample PFS event/time from GEO. This is NSCLC + PD-(L)1
blockade (the disease/therapy class Bessede studied) but small, single-arm-style,
mostly anti-PD-1, and NOT the OAK/POPLAR data. Hypothesis-generating only.

Source: GEO GSE135222 (Jung et al., Nat Commun 2019). TPM matrix
GSE135222_GEO_RNA-seq_omicslab_exp.tsv; PFS from series matrix characteristics.
"""
import json
import re
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

RAW = "results/hunt_oak/data/raw"
DERIVED = "results/hunt_oak/data/derived"
OUT = "results/hunt_oak/outputs"

ENSG = {"TACSTD2": "ENSG00000184292", "CD274": "ENSG00000120217",
        "CD8A": "ENSG00000153563", "EPCAM": "ENSG00000119888",
        "PTPRC": "ENSG00000081237", "ACTB": "ENSG00000075624"}

exp = pd.read_csv(f"{RAW}/GSE135222_exp.tsv", sep="\t", index_col=0)
exp.index = exp.index.str.split(".").str[0]  # strip Ensembl version
gene_tpm = {g: exp.loc[e] for g, e in ENSG.items() if e in exp.index}
tpm = pd.DataFrame(gene_tpm)  # samples x genes
tpm.index.name = "exp_sample"

# --- parse PFS event/time + sample titles from series matrix ---------------
titles, pfs_event, pfs_time = None, None, None
with open(f"{RAW}/GSE135222_series_matrix.txt") as fh:
    for line in fh:
        parts = [p.strip('"') for p in line.rstrip("\n").split("\t")]
        tag = parts[0]
        vals = parts[1:]
        if tag == "!Sample_title":
            titles = vals
        elif tag == "!Sample_characteristics_ch1":
            if any(v.startswith("progression-free survival (pfs):") for v in vals):
                pfs_event = [int(v.split(":")[1]) for v in vals]
            elif any(v.startswith("pfs.time:") for v in vals):
                pfs_time = [float(v.split(":")[1]) for v in vals]

clin = pd.DataFrame({"title": titles, "pfs_event": pfs_event, "pfs_time": pfs_time})
clin["exp_sample"] = clin["title"].str.replace(" ", "", regex=False)  # "NSCLC 990"->"NSCLC990"

df = clin.merge(tpm.reset_index(), on="exp_sample", how="inner")
assert len(df) == len(clin), "sample-title <-> expression-column mapping incomplete"

for g in ENSG:
    if g in df.columns:
        df[f"log2_{g}"] = np.log2(df[g] + 1.0)

# sanity: ACTB TPM must dwarf TACSTD2
assert df["ACTB"].median() > df["TACSTD2"].median(), "TPM extraction sanity failed"

med = df["log2_TACSTD2"].median()
df["TACSTD2_group"] = np.where(df["log2_TACSTD2"] >= med, "high", "low")

df.to_csv(f"{DERIVED}/gse135222_tacstd2_persample.csv", index=False)

res = {"dataset": "GSE135222 (NSCLC, anti-PD-1/PD-L1, single-arm, n=%d)" % len(df),
       "gene": "TACSTD2 (TROP2)", "n_total": int(len(df)),
       "pfs_events": int(df["pfs_event"].sum()),
       "note": "pfs_event coded 1=progression event, 0=censored (per GEO characteristics)"}

# PFS Cox, per-SD standardized continuous TACSTD2
d = df.copy()
d["z"] = (d.log2_TACSTD2 - d.log2_TACSTD2.mean()) / d.log2_TACSTD2.std()
cph = CoxPHFitter().fit(d[["pfs_time", "pfs_event", "z"]], "pfs_time", "pfs_event")
row = cph.summary.loc["z"]
res["pfs_cox_perSD_HR"] = round(float(row["exp(coef)"]), 4)
res["pfs_cox_perSD_HR_95CI"] = [round(float(row["exp(coef) lower 95%"]), 4),
                                round(float(row["exp(coef) upper 95%"]), 4)]
res["pfs_cox_perSD_p"] = round(float(row["p"]), 5)

hi, lo = d[d.TACSTD2_group == "high"], d[d.TACSTD2_group == "low"]
lr = logrank_test(hi.pfs_time, lo.pfs_time, hi.pfs_event, lo.pfs_event)
res["pfs_logrank_p_medianSplit"] = round(float(lr.p_value), 5)
kmf = KaplanMeierFitter()
for name, g in (("high", hi), ("low", lo)):
    kmf.fit(g.pfs_time, g.pfs_event)
    m = kmf.median_survival_time_
    res[f"pfs_median_days_TACSTD2_{name}"] = None if np.isinf(m) else round(float(m), 1)

for other in ["CD8A", "CD274"]:
    if f"log2_{other}" in df.columns:
        rho, p = stats.spearmanr(df.log2_TACSTD2, df[f"log2_{other}"])
        res[f"spearman_TACSTD2_vs_{other}_rho"] = round(float(rho), 4)
        res[f"spearman_TACSTD2_vs_{other}_p"] = round(float(p), 5)

with open(f"{OUT}/gse135222_results.json", "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
