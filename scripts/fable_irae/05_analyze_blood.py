#!/usr/bin/env python3
"""GSE319496 - pretreatment whole-blood bulk RNA-seq, irAE Yes vs No.

mRCC patients on nivolumab + ipilimumab, stratified by whether they later developed
an immune-related adverse event (irAE). This is a systemic/blood compartment
verification of epithelial-marker measurability and a clean case/control DE.

X = raw gene counts (22,986 genes x 51 samples). Sample columns are Sample01..51;
irAE status comes from the series-matrix characteristics (aligned by sample title).

Outputs:
  results/fable_irae/tables/blood_GSE319496_goi.tsv
  results/fable_irae/tables/blood_GSE319496_de_summary.tsv
  results/fable_irae/figures/blood_GSE319496_CLDN4.png
"""
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] / "results" / "fable_irae"
RAW = ROOT / "data" / "raw"
TAB = ROOT / "tables"; TAB.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "figures"; FIG.mkdir(parents=True, exist_ok=True)

GOI = ["TACSTD2", "CLDN4"]


def load_phenotype():
    txt = gzip.decompress((RAW / "GSE319496_series_matrix.txt.gz").read_bytes()).decode(errors="replace")
    titles = accs = irae = None
    for line in txt.splitlines():
        cells = [x.strip('"') for x in line.split("\t")]
        if line.startswith("!Sample_title"):
            titles = cells[1:]
        elif line.startswith("!Sample_geo_accession"):
            accs = cells[1:]
        elif line.startswith("!Sample_characteristics_ch1") and "irae status" in line:
            irae = [c.split(":", 1)[1].strip() for c in cells[1:]]
    # sample id from title "... Sample01"
    samp = [t.split(",")[-1].strip() for t in titles]
    pheno = pd.DataFrame({"gsm": accs, "sample": samp, "irae": irae})
    return pheno


def main():
    counts = pd.read_csv(RAW / "GSE319496_GEO_raw_counts_SampleID.csv.gz", index_col=0)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    pheno = load_phenotype()
    pheno = pheno[pheno["sample"].isin(counts.columns)].reset_index(drop=True)
    print("phenotype irAE counts:", pheno["irae"].value_counts().to_dict())

    # CPM + log2
    lib = counts.sum(axis=0)
    cpm = counts.divide(lib, axis=1) * 1e6
    logcpm = np.log2(cpm + 1.0)

    # measurability of GOI
    goi_rows = []
    present = {g: (g in counts.index) for g in GOI}
    for g in GOI:
        if present[g]:
            raw_g = counts.loc[g]
            goi_rows.append(dict(gene=g, in_matrix=True,
                                 n_samples_detected=int((raw_g > 0).sum()),
                                 pct_detected=float((raw_g > 0).mean()),
                                 mean_cpm=float(cpm.loc[g].mean()),
                                 median_cpm=float(cpm.loc[g].median())))
        else:
            goi_rows.append(dict(gene=g, in_matrix=False, n_samples_detected=0,
                                 pct_detected=0.0, mean_cpm=np.nan, median_cpm=np.nan))
    pd.DataFrame(goi_rows).to_csv(TAB / "blood_GSE319496_goi.tsv", sep="\t", index=False)
    print("[write] blood GOI measurability")
    print(pd.DataFrame(goi_rows).to_string(index=False))

    # DE irAE Yes vs No for measurable GOI (+ genome-wide context for CLDN4)
    yes = pheno.loc[pheno["irae"] == "Yes", "sample"]
    no = pheno.loc[pheno["irae"] == "No", "sample"]
    de_rows = []
    for g in GOI:
        if not present[g]:
            de_rows.append(dict(gene=g, measurable=False))
            continue
        a = logcpm.loc[g, yes].values
        b = logcpm.loc[g, no].values
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        l2fc = cpm.loc[g, yes].mean()
        l2fc = np.log2((cpm.loc[g, yes].mean() + 1e-3) / (cpm.loc[g, no].mean() + 1e-3))
        de_rows.append(dict(gene=g, measurable=True, n_yes=len(a), n_no=len(b),
                            mean_cpm_yes=cpm.loc[g, yes].mean(), mean_cpm_no=cpm.loc[g, no].mean(),
                            log2FC_yes_over_no=l2fc, mwu_U=u, p_value=p))
    de = pd.DataFrame(de_rows)
    de.to_csv(TAB / "blood_GSE319496_de_summary.tsv", sep="\t", index=False)
    print("[write] blood DE summary")
    print(de.to_string(index=False))

    # figure for CLDN4 if measurable
    if present["CLDN4"]:
        fig, ax = plt.subplots(figsize=(4, 4))
        data = [logcpm.loc["CLDN4", no].values, logcpm.loc["CLDN4", yes].values]
        ax.boxplot(data, tick_labels=["irAE: No", "irAE: Yes"], showfliers=False)
        for i, d in enumerate(data, 1):
            ax.scatter(np.random.normal(i, 0.06, len(d)), d, s=20, alpha=0.7,
                       color="#c0392b" if i == 2 else "#2c7fb8")
        p = de.loc[de.gene == "CLDN4", "p_value"].values[0]
        ax.set_title(f"CLDN4 whole blood (GSE319496)\nMWU p={p:.3g}")
        ax.set_ylabel("log2(CPM+1)")
        fig.tight_layout()
        fig.savefig(FIG / "blood_GSE319496_CLDN4.png", dpi=140)
        print("[write] blood figure")


if __name__ == "__main__":
    main()
