"""
Hunt: TACSTD2 (TROP2) vs immune signature, and purity-partial Spearman,
in openly downloadable NSCLC RNA datasets.

Datasets actually used (open, processed, << 2 GB):
  1. GSE248378  - neoadjuvant DURVALUMAB (anti-PD-L1) +/- radiation, NSCLC,
                  bulk RNA-seq FPKM of resected tumours (post-treatment). OPEN (GEO).
  2. luad_oncosg_2020 (OncoSG, Nat Genet 2020) - LUAD; cBioPortal ships per-sample
                  tumour PURITY and IMSIG immune-cell signatures + RNA-seq. OPEN.

The specific anti-PD-L1 cohort behind the reported TACSTD2/TROP2-immune result
(Bessede et al., Clin Cancer Res 2024) is POPLAR+OAK atezolizumab RNA-seq, which is
EGA CONTROLLED-ACCESS (EGAS00001005013) -- NOT an open <2GB download. The PACIFIC
durvalumab trial did not collect mandatory tissue and its data are only available via
AstraZeneca's request portal. SUBMARINE (WJOG11518L) durvalumab transcriptome has no
public accession. See README for full provenance.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

import estimate_ssgsea as est

DATA = "data"
OUT = "."


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def partial_spearman(x, y, z):
    """Rank-based partial Spearman of x,y controlling for z (one covariate)."""
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]
    n = m.sum()
    rx = stats.rankdata(x); ry = stats.rankdata(y); rz = stats.rankdata(z)
    rxy = np.corrcoef(rx, ry)[0, 1]
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]
    denom = np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    rp = (rxy - rxz * ryz) / denom
    k = 1  # covariates
    df = n - 2 - k
    t = rp * np.sqrt(df / (1 - rp ** 2))
    pval = 2 * stats.t.sf(abs(t), df)
    return float(rp), float(pval), int(n)


results = {"datasets": {}}

# ============================================================================
# 1. GSE248378 - DURVALUMAB (anti-PD-L1) NSCLC bulk RNA-seq (FPKM)
# ============================================================================
print("=== GSE248378 (durvalumab, anti-PD-L1) ===")
fpkm = pd.read_csv(f"{DATA}/GSE248378_Durva_Post_FPKMs.txt.gz", sep="\t", index_col=0)
fpkm.index = fpkm.index.astype(str)
gsets = est.load_gmt(f"{DATA}/ESTIMATE_SI_geneset.gmt")
common = est.load_common_genes(f"{DATA}/ESTIMATE_common_genes.txt")
scores, ng = est.estimate_scores(fpkm, gsets, common, platform="illumina")
scores.to_csv(f"{OUT}/GSE248378_estimate_scores.csv")

if "TACSTD2" not in fpkm.index:
    raise SystemExit("TACSTD2 not found in GSE248378 FPKM matrix")
tac = fpkm.loc["TACSTD2"]
if isinstance(tac, pd.DataFrame):
    tac = tac.mean(axis=0)
df = scores.copy()
df["TACSTD2_FPKM"] = tac.reindex(df.index).astype(float)
df.to_csv(f"{OUT}/GSE248378_merged_table.csv")

r_imm, p_imm, n_imm = spearman(df["TACSTD2_FPKM"], df["ImmuneScore"])
r_est, p_est, _ = spearman(df["TACSTD2_FPKM"], df["ESTIMATEScore"])
r_pur, p_pur, _ = spearman(df["TACSTD2_FPKM"], df["TumorPurity_ESTIMATE_affyproxy"])
r_ip, p_ip, _ = spearman(df["ImmuneScore"], df["TumorPurity_ESTIMATE_affyproxy"])
pr_imm, pp_imm, n_p = partial_spearman(
    df["TACSTD2_FPKM"], df["ImmuneScore"], df["TumorPurity_ESTIMATE_affyproxy"])

print(f"  n = {n_imm}")
print(f"  TACSTD2 vs ImmuneScore   rho = {r_imm:+.3f}  p = {p_imm:.2e}")
print(f"  TACSTD2 vs ESTIMATEScore rho = {r_est:+.3f}  p = {p_est:.2e}")
print(f"  purity-partial (ctrl ESTIMATE-affy purity) rho = {pr_imm:+.3f}  p = {pp_imm:.2e}")

results["datasets"]["GSE248378_durvalumab"] = {
    "drug": "durvalumab (anti-PD-L1)", "context": "neoadjuvant +/- radiation, resected NSCLC (post-treatment)",
    "source": "GEO GSE248378 (open); file GSE248378_Durva_Post_FPKMs.txt.gz",
    "n_samples": n_imm, "n_genes_estimate_common": ng,
    "immune_metric": "ESTIMATE ImmuneScore (ssGSEA)",
    "purity_metric": "ESTIMATE ESTIMATEScore->Affy purity formula (PROXY; Affy-calibrated, off-label for RNA-seq)",
    "spearman_TACSTD2_vs_ImmuneScore": {"rho": r_imm, "p": p_imm},
    "spearman_TACSTD2_vs_ESTIMATEScore": {"rho": r_est, "p": p_est},
    "spearman_TACSTD2_vs_purityProxy": {"rho": r_pur, "p": p_pur},
    "spearman_ImmuneScore_vs_purityProxy": {"rho": r_ip, "p": p_ip},
    "purity_partial_spearman_TACSTD2_vs_ImmuneScore": {"rho": pr_imm, "p": pp_imm, "n": n_p},
    "caveats": ["ESTIMATE tumour-purity formula is Affymetrix-calibrated; on RNA-seq it is only a rank proxy.",
                "Small n; post-treatment resected tumours, not the reported anti-PD-L1 discovery cohort."],
}

# ============================================================================
# 2. OncoSG luad_oncosg_2020 - LUAD (real purity + IMSIG immune signatures)
# ============================================================================
print("\n=== OncoSG luad_oncosg_2020 (LUAD) ===")
clin_lines = [l.rstrip("\n") for l in open(f"{DATA}/oncosg_data_clinical_sample.txt")
              if not l.startswith("#")]
hdr = clin_lines[0].split("\t")
clin = pd.DataFrame([l.split("\t") for l in clin_lines[1:]], columns=hdr)
clin = clin.set_index("SAMPLE_ID")
tac_z = pd.read_csv(f"{DATA}/oncosg_TACSTD2_zscore_all.tsv", sep="\t",
                    header=None, names=["SAMPLE_ID", "TACSTD2_z"]).set_index("SAMPLE_ID")

imsig_cols = ["IMSIG_B_CELLS", "IMSIG_T_CELLS", "IMSIG_NK_CELLS", "IMSIG_MACROPHAGES",
              "IMSIG_MONOCYTES", "IMSIG_NEUTROPHILS", "IMSIG_PLASMA_CELLS", "IMSIG_INTERFERON"]
o = clin.join(tac_z, how="inner")
for c in imsig_cols + ["PURITY"]:
    o[c] = pd.to_numeric(o[c], errors="coerce")
o["TACSTD2_z"] = pd.to_numeric(o["TACSTD2_z"], errors="coerce")
# composite immune infiltration = mean of z-scored IMSIG immune-cell signatures
z = o[imsig_cols].apply(lambda s: (s - s.mean()) / s.std(ddof=0))
o["IMSIG_immune_composite"] = z.mean(axis=1)
o = o.dropna(subset=["TACSTD2_z", "PURITY", "IMSIG_immune_composite"])
o.to_csv(f"{OUT}/OncoSG_merged_table.csv")

r_imm, p_imm, n_imm = spearman(o["TACSTD2_z"], o["IMSIG_immune_composite"])
r_pur, p_pur, _ = spearman(o["TACSTD2_z"], o["PURITY"])
r_ip, p_ip, _ = spearman(o["IMSIG_immune_composite"], o["PURITY"])
pr_imm, pp_imm, n_p = partial_spearman(o["TACSTD2_z"], o["IMSIG_immune_composite"], o["PURITY"])
print(f"  n = {n_imm}")
print(f"  TACSTD2 vs IMSIG immune composite  rho = {r_imm:+.3f}  p = {p_imm:.2e}")
print(f"  TACSTD2 vs PURITY                  rho = {r_pur:+.3f}  p = {p_pur:.2e}")
print(f"  purity-partial (ctrl real PURITY)  rho = {pr_imm:+.3f}  p = {pp_imm:.2e}")

per_sig = {}
for c in imsig_cols:
    rr, pp, nn = spearman(o["TACSTD2_z"], o[c])
    prr, ppp, _ = partial_spearman(o["TACSTD2_z"], o[c], o["PURITY"])
    per_sig[c] = {"spearman_rho": rr, "spearman_p": pp,
                  "purity_partial_rho": prr, "purity_partial_p": ppp}
    print(f"    {c:18s} rho={rr:+.3f} (p={pp:.1e})  partial={prr:+.3f} (p={ppp:.1e})")

results["datasets"]["OncoSG_luad_oncosg_2020"] = {
    "drug": "none (treatment-naive LUAD resection cohort, East-Asian)",
    "source": "cBioPortal study luad_oncosg_2020 (Chen et al., Nat Genet 2020); open",
    "n_samples": n_imm,
    "immune_metric": "IMSIG immune-cell composite (mean of z-scored IMSIG immune signatures), from cBioPortal clinical",
    "purity_metric": "OncoSG per-sample tumour PURITY (real, from cBioPortal clinical) -- independent of immune score",
    "spearman_TACSTD2_vs_immuneComposite": {"rho": r_imm, "p": p_imm},
    "spearman_TACSTD2_vs_PURITY": {"rho": r_pur, "p": p_pur},
    "spearman_immuneComposite_vs_PURITY": {"rho": r_ip, "p": p_ip},
    "purity_partial_spearman_TACSTD2_vs_immuneComposite": {"rho": pr_imm, "p": pp_imm, "n": n_p},
    "per_signature": per_sig,
}

with open(f"{OUT}/results_summary.json", "w") as fh:
    json.dump(results, fh, indent=2)
print("\nWrote results_summary.json")
