#!/usr/bin/env python3
"""GSE283829 leftover: TACSTD2 (TROP2) and CLDN4 vs ICI response in NSCLC.

Cohort: 27 ICI-treated lung-cancer RNA-seq samples deposited with
Lindberg et al., J Thorac Oncol 2025 (PMID 39743139; GEO GSE283829).
The series matrix has no expression table; counts come from the public
supplementary raw-count matrix.

GEO labels the RECIST-like field "disease stage". Values are CR / SD / PD
(no PR). That field is used as best response, not TNM stage.

Primary contrast (pre-specified): CR vs PD on log2(CPM+1).
Sensitivity: CR vs SD+PD; CR+SD vs PD; Kruskal-Wallis across CR/SD/PD;
Spearman vs ordinal response (CR=2, SD=1, PD=0).
Exploratory: PD1–PD-L1 PLA high vs low (the authors' grouping).
Context genes (not primary): CD274, PDCD1, and the paper's PLA-high
non-responder mediators EOMES, HAVCR1, JAML, FCRL1.

Usage: python3 scripts/gse283829_tacstd2_cldn4_response.py
Outputs: results/w200/GSE283829/
"""

from __future__ import annotations

import gzip
import hashlib
import platform
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "GSE283829"
OUT_DIR = ROOT / "results" / "w200" / "GSE283829"

SERIES_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc=GSE283829&targ=self&form=text&view=brief"
)
SAMPLES_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc=GSE283829&targ=gsm&form=text&view=brief"
)
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283829/suppl/"
    "GSE283829_raw_express_matrix_all_samples.txt.gz"
)
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
    "hgnc_complete_set.tsv"
)

PRIMARY = {"TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143"}
CONTEXT = {
    "CD274": "ENSG00000120217",
    "PDCD1": "ENSG00000188389",
    "EOMES": "ENSG00000163508",
    "HAVCR1": "ENSG00000138722",
    "JAML": "ENSG00000160593",
    "FCRL1": "ENSG00000163534",
}
ALL_GENES = {**PRIMARY, **CONTEXT}

N_BOOT = 2000
SEED = 20260816
ORDINAL = {"PD": 0, "SD": 1, "CR": 2}


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        urllib.request.urlretrieve(url, dest)
    return dest


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_soft_samples(path: Path) -> pd.DataFrame:
    text = path.read_text(errors="replace")
    blocks = re.split(r"\n(?=\^SAMPLE = )", text)
    rows = []
    for block in blocks:
        m = re.search(r"\^SAMPLE = (GSM\d+)", block)
        if not m:
            continue
        gsm = m.group(1)
        title = re.search(r"!Sample_title = (.+)", block)
        chars = {}
        for ln in re.findall(r"!Sample_characteristics_ch1 = (.+)", block):
            if ": " in ln:
                k, v = ln.split(": ", 1)
                chars[k.strip()] = v.strip()
        rows.append(
            {
                "gsm": gsm,
                "title": title.group(1).strip() if title else "",
                "sex": chars.get("Sex"),
                "response": chars.get("disease stage"),
                "tumor_type": chars.get("tumor type"),
                "pla": chars.get("pla"),
                "batch": chars.get("batch"),
            }
        )
    df = pd.DataFrame(rows)
    if len(df) != 27:
        raise ValueError(f"expected 27 GSM records, got {len(df)}")
    if set(df["response"]) != {"CR", "SD", "PD"}:
        raise ValueError(f"unexpected response labels: {sorted(df.response.unique())}")
    df["sample"] = "S" + df["title"].str.replace("-", "_")
    df["response_ordinal"] = df["response"].map(ORDINAL)
    df["is_CR"] = (df["response"] == "CR").astype(int)
    df["is_PD"] = (df["response"] == "PD").astype(int)
    df["is_DCR"] = df["response"].isin(["CR", "SD"]).astype(int)
    return df


def load_counts(path: Path) -> pd.DataFrame:
    counts = pd.read_csv(path, sep="\t", index_col=0)
    counts.index = [str(i).strip('"') for i in counts.index]
    counts.columns = [str(c).strip('"') for c in counts.columns]
    if counts.shape[1] != 27:
        raise ValueError(f"expected 27 count columns, got {counts.shape[1]}")
    return counts.astype(float)


def cpm_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0)
    return counts.divide(lib, axis=1) * 1e6


def load_hgnc(path: Path) -> pd.DataFrame:
    hgnc = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    hgnc = hgnc.dropna(subset=["ensembl_gene_id", "symbol"])
    return hgnc


def auc_mwu(values: np.ndarray, labels: np.ndarray) -> float:
    """AUC with labels==1 as the positive class (higher value -> more positive)."""
    pos, neg = values[labels == 1], values[labels == 0]
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(u / (len(pos) * len(neg)))


def bootstrap_auc_ci(values, labels, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    pos_idx, neg_idx = np.where(labels == 1)[0], np.where(labels == 0)[0]
    aucs = []
    for _ in range(n_boot):
        idx = np.concatenate(
            [
                rng.choice(pos_idx, len(pos_idx), replace=True),
                rng.choice(neg_idx, len(neg_idx), replace=True),
            ]
        )
        aucs.append(auc_mwu(values[idx], labels[idx]))
    return tuple(np.percentile(aucs, [2.5, 97.5]))


def bh_adjust(p_values: list[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1.0)
    return out


def contrast_stats(df: pd.DataFrame, col: str, label_col: str, pos_name: str, neg_name: str) -> dict:
    labels = df[label_col].to_numpy()
    v = df[col].to_numpy(float)
    pos, neg = v[labels == 1], v[labels == 0]
    mwu = stats.mannwhitneyu(pos, neg, alternative="two-sided", method="exact")
    auc = auc_mwu(v, labels)
    lo, hi = bootstrap_auc_ci(v, labels)
    return {
        "feature": col.replace("_log2CPM1", ""),
        "contrast": f"{pos_name}_vs_{neg_name}",
        "n_pos": int(len(pos)),
        "n_neg": int(len(neg)),
        "median_pos": float(np.median(pos)),
        "median_neg": float(np.median(neg)),
        "delta_median_pos_minus_neg": float(np.median(pos) - np.median(neg)),
        "mean_pos": float(np.mean(pos)),
        "mean_neg": float(np.mean(neg)),
        "mannwhitney_U": float(mwu.statistic),
        "p_exact_two_sided": float(mwu.pvalue),
        "AUC_high_predicts_pos": auc,
        "AUC_95CI_lo": float(lo),
        "AUC_95CI_hi": float(hi),
    }


def genome_wide_rank(cpm: pd.DataFrame, labels: np.ndarray, protein_ids: set[str]) -> pd.DataFrame:
    keep = [g for g in cpm.index if g in protein_ids and float(cpm.loc[g].median()) >= 1.0]
    mat = np.log2(cpm.loc[keep].to_numpy() + 1.0)
    pos = labels == 1
    neg = labels == 0
    rows = []
    for i, gid in enumerate(keep):
        a, b = mat[i, pos], mat[i, neg]
        mwu = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append(
            {
                "ensembl_gene_id": gid,
                "median_pos": float(np.median(a)),
                "median_neg": float(np.median(b)),
                "delta_median": float(np.median(a) - np.median(b)),
                "p": float(mwu.pvalue),
                "U": float(mwu.statistic),
            }
        )
    out = pd.DataFrame(rows)
    out["rank_by_p"] = out["p"].rank(method="min").astype(int)
    out["n_tested"] = len(out)
    return out.sort_values("p")


def jitter(n, rng, scale=0.08):
    return rng.uniform(-scale, scale, n)


def box_strip(ax, groups, labels, colors, rng):
    ax.boxplot(groups, tick_labels=labels, showfliers=False, widths=0.55)
    for i, (g, c) in enumerate(zip(groups, colors), start=1):
        ax.scatter(
            np.full(len(g), i) + jitter(len(g), rng),
            g,
            s=28,
            alpha=0.85,
            color=c,
            zorder=3,
            edgecolors="k",
            linewidths=0.3,
        )


def fmt(x, d=3):
    if pd.isna(x):
        return "NA"
    x = float(x)
    if x == 0:
        return "0"
    if abs(x) < 0.001:
        return f"{x:.2e}"
    return f"{x:.{d}f}"


def write_readme(df, stats_df, gw_crpd, gw_cr, spearman_rows, kw_rows, corr):
    def row(feature, contrast):
        hit = stats_df[(stats_df.feature == feature) & (stats_df.contrast == contrast)]
        if hit.empty:
            raise KeyError((feature, contrast))
        return hit.iloc[0]

    n = len(df)
    n_cr = int((df.response == "CR").sum())
    n_sd = int((df.response == "SD").sum())
    n_pd = int((df.response == "PD").sum())
    n_pla_h = int((df.pla == "high").sum())
    n_pla_l = int((df.pla == "low").sum())
    n_ac = int((df.tumor_type == "AC").sum())
    n_sq = int((df.tumor_type == "SqCC").sum())
    n_oth = int((df.tumor_type == "other").sum())

    tac_crpd = row("TACSTD2", "CR_vs_PD")
    cld_crpd = row("CLDN4", "CR_vs_PD")
    z_crpd = row("combined_z", "CR_vs_PD")
    tac_cr = row("TACSTD2", "CR_vs_nonCR")
    cld_cr = row("CLDN4", "CR_vs_nonCR")
    tac_dcr = row("TACSTD2", "DCR_vs_PD")
    cld_dcr = row("CLDN4", "DCR_vs_PD")
    tac_pla = row("TACSTD2", "PLA_high_vs_low")
    cld_pla = row("CLDN4", "PLA_high_vs_low")

    def gw_line(table, gene, ensg):
        hit = table[table.ensembl_gene_id == ensg]
        if hit.empty:
            return f"{gene} not in the filtered protein-coding set"
        r = hit.iloc[0]
        return (
            f"{gene} rank {int(r.rank_by_p)} / {int(r.n_tested)} "
            f"(p={fmt(r.p)}, Δmedian={fmt(r.delta_median)})"
        )

    tac_sp = next(x for x in spearman_rows if x["feature"] == "TACSTD2")
    cld_sp = next(x for x in spearman_rows if x["feature"] == "CLDN4")
    tac_kw = next(x for x in kw_rows if x["feature"] == "TACSTD2")
    cld_kw = next(x for x in kw_rows if x["feature"] == "CLDN4")

    # Direction honesty
    def dir_note(acc):
        d = acc.delta_median_pos_minus_neg
        if d > 0:
            return "higher in the positive class"
        if d < 0:
            return "lower in the positive class"
        return "identical medians"

    readme = f"""# GSE283829 leftover: TACSTD2 / CLDN4 vs ICI response (NSCLC)

**Bottom line: no persuasive association.** In this public leftover ICI
RNA-seq cohort (n = {n}; CR {n_cr} / SD {n_sd} / PD {n_pd}), neither
TACSTD2 nor CLDN4 tumor expression separates complete responders from
progressors. Point estimates are small, confidence intervals include
chance, and both genes sit in the middle of the genome-wide rank — not
near a response-associated tail.

This is a leftover GEO series: the series matrix has **no expression
table**. Counts were taken from the public supplementary raw-count
matrix. The analysis is open and was run from those files.

## Cohort

- GEO [GSE283829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283829),
  Lindberg et al., *J Thorac Oncol* (2025), PMID 39743139.
- 27 ICI-treated lung-cancer RNA-seq samples (Illumina NovaSeq 6000,
  GPL24676). The parent study profiled PD1–PD-L1 interactions by in situ
  PLA in a larger biopsy series; this deposit is the RNA-seq subset used
  to look for resistance programs in PLA-high non-responders.
- GEO field `disease stage` is **not TNM**. Values are CR / SD / PD and
  are treated as best response. There is **no PR** category in the
  deposit (either none were sequenced or PR was not labeled).
- Histology: AC {n_ac}, SqCC {n_sq}, other {n_oth}. PLA: high {n_pla_h},
  low {n_pla_l}. Two library batches (batch 2 / batch 3).
- Genes: TACSTD2 = ENSG00000184292, CLDN4 = ENSG00000189143. Expression
  is log2(CPM+1) from the deposited raw counts. No gene-length / TPM
  conversion is possible from this file alone.

## Pre-specified tests

Primary: two-sided exact Mann–Whitney U on log2(CPM+1), **CR vs PD**.
BH q is across the two primary genes only. AUC treats CR as the positive
class (higher expression predicting CR). Sensitivity contrasts and the
PLA split are reported as such and were not folded into the two-gene FDR.

| Feature | Contrast | n | Median pos | Median neg | Δ median | Exact MWU p | BH q (2 genes) | AUC (95% CI) |
|---|---|---|---:|---:|---:|---:|---:|---|
| TACSTD2 | CR vs PD | {int(tac_crpd.n_pos)} vs {int(tac_crpd.n_neg)} | {fmt(tac_crpd.median_pos)} | {fmt(tac_crpd.median_neg)} | {fmt(tac_crpd.delta_median_pos_minus_neg)} | {fmt(tac_crpd.p_exact_two_sided)} | {fmt(tac_crpd.q_BH_2genes)} | {fmt(tac_crpd.AUC_high_predicts_pos)} ({fmt(tac_crpd.AUC_95CI_lo)}–{fmt(tac_crpd.AUC_95CI_hi)}) |
| CLDN4 | CR vs PD | {int(cld_crpd.n_pos)} vs {int(cld_crpd.n_neg)} | {fmt(cld_crpd.median_pos)} | {fmt(cld_crpd.median_neg)} | {fmt(cld_crpd.delta_median_pos_minus_neg)} | {fmt(cld_crpd.p_exact_two_sided)} | {fmt(cld_crpd.q_BH_2genes)} | {fmt(cld_crpd.AUC_high_predicts_pos)} ({fmt(cld_crpd.AUC_95CI_lo)}–{fmt(cld_crpd.AUC_95CI_hi)}) |
| combined z (exploratory) | CR vs PD | {int(z_crpd.n_pos)} vs {int(z_crpd.n_neg)} | {fmt(z_crpd.median_pos)} | {fmt(z_crpd.median_neg)} | {fmt(z_crpd.delta_median_pos_minus_neg)} | {fmt(z_crpd.p_exact_two_sided)} | — | {fmt(z_crpd.AUC_high_predicts_pos)} ({fmt(z_crpd.AUC_95CI_lo)}–{fmt(z_crpd.AUC_95CI_hi)}) |
| TACSTD2 | CR vs SD+PD | {int(tac_cr.n_pos)} vs {int(tac_cr.n_neg)} | {fmt(tac_cr.median_pos)} | {fmt(tac_cr.median_neg)} | {fmt(tac_cr.delta_median_pos_minus_neg)} | {fmt(tac_cr.p_exact_two_sided)} | — | {fmt(tac_cr.AUC_high_predicts_pos)} ({fmt(tac_cr.AUC_95CI_lo)}–{fmt(tac_cr.AUC_95CI_hi)}) |
| CLDN4 | CR vs SD+PD | {int(cld_cr.n_pos)} vs {int(cld_cr.n_neg)} | {fmt(cld_cr.median_pos)} | {fmt(cld_cr.median_neg)} | {fmt(cld_cr.delta_median_pos_minus_neg)} | {fmt(cld_cr.p_exact_two_sided)} | — | {fmt(cld_cr.AUC_high_predicts_pos)} ({fmt(cld_cr.AUC_95CI_lo)}–{fmt(cld_cr.AUC_95CI_hi)}) |
| TACSTD2 | CR+SD vs PD | {int(tac_dcr.n_pos)} vs {int(tac_dcr.n_neg)} | {fmt(tac_dcr.median_pos)} | {fmt(tac_dcr.median_neg)} | {fmt(tac_dcr.delta_median_pos_minus_neg)} | {fmt(tac_dcr.p_exact_two_sided)} | — | {fmt(tac_dcr.AUC_high_predicts_pos)} ({fmt(tac_dcr.AUC_95CI_lo)}–{fmt(tac_dcr.AUC_95CI_hi)}) |
| CLDN4 | CR+SD vs PD | {int(cld_dcr.n_pos)} vs {int(cld_dcr.n_neg)} | {fmt(cld_dcr.median_pos)} | {fmt(cld_dcr.median_neg)} | {fmt(cld_dcr.delta_median_pos_minus_neg)} | {fmt(cld_dcr.p_exact_two_sided)} | — | {fmt(cld_dcr.AUC_high_predicts_pos)} ({fmt(cld_dcr.AUC_95CI_lo)}–{fmt(cld_dcr.AUC_95CI_hi)}) |

Kruskal–Wallis across CR/SD/PD: TACSTD2 p={fmt(tac_kw['p'])}, CLDN4
p={fmt(cld_kw['p'])}. Spearman vs ordinal response (CR=2, SD=1, PD=0):
TACSTD2 ρ={fmt(tac_sp['spearman_r'])} (p={fmt(tac_sp['p'])}), CLDN4
ρ={fmt(cld_sp['spearman_r'])} (p={fmt(cld_sp['p'])}).

TACSTD2 vs CLDN4 Spearman ρ={fmt(corr.statistic)} (p={fmt(corr.pvalue)}).
The two genes are correlated, so they are not independent tests of a
tight-junction / TROP2 axis.

Direction on the primary contrast: TACSTD2 is {dir_note(tac_crpd)};
CLDN4 is {dir_note(cld_crpd)}. Neither is a large, consistent shift.

## PLA split (exploratory; authors' grouping)

The deposit is grouped by PD1–PD-L1 PLA, not by a TACSTD2/CLDN4
hypothesis. PLA-high vs PLA-low:

| Feature | Median high | Median low | Exact MWU p | AUC (high as +) |
|---|---:|---:|---:|---|
| TACSTD2 | {fmt(tac_pla.median_pos)} | {fmt(tac_pla.median_neg)} | {fmt(tac_pla.p_exact_two_sided)} | {fmt(tac_pla.AUC_high_predicts_pos)} ({fmt(tac_pla.AUC_95CI_lo)}–{fmt(tac_pla.AUC_95CI_hi)}) |
| CLDN4 | {fmt(cld_pla.median_pos)} | {fmt(cld_pla.median_neg)} | {fmt(cld_pla.p_exact_two_sided)} | {fmt(cld_pla.AUC_high_predicts_pos)} ({fmt(cld_pla.AUC_95CI_lo)}–{fmt(cld_pla.AUC_95CI_hi)}) |

This is not a response test. It is included so a leftover scan does not
confuse the authors' PLA grouping with RECIST.

## Genome-wide calibration (protein-coding, median CPM ≥ 1)

CR vs PD: {gw_line(gw_crpd, 'TACSTD2', PRIMARY['TACSTD2'])};
{gw_line(gw_crpd, 'CLDN4', PRIMARY['CLDN4'])}.

CR vs SD+PD: {gw_line(gw_cr, 'TACSTD2', PRIMARY['TACSTD2'])};
{gw_line(gw_cr, 'CLDN4', PRIMARY['CLDN4'])}.

A leftover claim that these two genes mark ICI response would require
them to stand out. They do not.

## Honest caveats

- **n = 27, 7 CRs.** Only a large effect (AUC ≳ 0.85) would be reliably
  detected. A null here does not prove no association; it also does not
  support one. AUC intervals all include 0.5 on the primary contrast.
- No PR labels, no PFS/OS in GEO, no PD-L1 IHC percent, no treatment
  line, no EGFR/ALK, no purity. Single-arm ICI series: prognostic vs
  predictive cannot be separated.
- Two library batches. No batch-adjusted primary model was fit; that
  would be underpowered and easy to overfit. Library sizes and
  per-sample values are in `sample_data.csv`.
- Bulk diagnostic-biopsy RNA. TACSTD2/CLDN4 are epithelial; composition
  and purity can move both genes without a cell-intrinsic ICI effect.
- Context genes (CD274, PDCD1, EOMES, HAVCR1, JAML, FCRL1) are in
  `context_stats.csv` as a sanity check against the paper, not as a
  new biomarker hunt.
- No cutpoint search, no multivariable classifier, no in-sample ROC
  optimization.

## Files

- `sample_data.csv` — per-sample GEO metadata, library size, CPM and
  log2(CPM+1) for the primary and context genes, combined z.
- `stats.csv` — Mann–Whitney / AUC for all reported contrasts.
- `kruskal_spearman.csv` — three-level and ordinal tests.
- `gene_correlation.csv` — TACSTD2 vs CLDN4 Spearman.
- `genomewide_rank_CR_vs_PD.csv`, `genomewide_rank_CR_vs_nonCR.csv` —
  ranks of the two genes plus the rest of the filtered set.
- `context_stats.csv` — context-gene CR vs PD only.
- `boxplots_response.png`, `boxplots_cr_vs_pd.png`,
  `scatter_tacstd2_cldn4.png`.
- `analysis_manifest.txt` — checksums, versions, seed.
- Code: `scripts/gse283829_tacstd2_cldn4_response.py` (downloads if
  needed; seed {SEED}).

## Sources

1. GEO [GSE283829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283829).
2. Lindberg et al. *J Thorac Oncol* (2025), PMID 39743139.
"""
    (OUT_DIR / "README.md").write_text(readme)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    series_path = fetch(SERIES_URL, DATA_DIR / "GSE283829_series.txt")
    samples_path = fetch(SAMPLES_URL, DATA_DIR / "GSE283829_samples.txt")
    counts_path = fetch(COUNTS_URL, DATA_DIR / "GSE283829_raw_express_matrix_all_samples.txt.gz")
    hgnc_path = fetch(HGNC_URL, DATA_DIR / "hgnc_complete_set.tsv")

    clin = parse_soft_samples(samples_path)
    counts = load_counts(counts_path)
    missing = [s for s in clin["sample"] if s not in counts.columns]
    if missing:
        # titles are 105691_026; matrix is S105691_026
        clin["sample"] = "S" + clin["title"]
        missing = [s for s in clin["sample"] if s not in counts.columns]
        if missing:
            raise ValueError(f"sample IDs not in count matrix: {missing}")
    clin = clin.set_index("sample").loc[counts.columns].reset_index()
    if len(clin) != 27:
        raise ValueError("failed to align 27 samples")

    cpm = cpm_matrix(counts)
    lib = counts.sum(axis=0)
    clin["library_size"] = clin["sample"].map(lib)

    for name, ensg in ALL_GENES.items():
        if ensg not in counts.index:
            raise ValueError(f"{name}/{ensg} missing from count matrix")
        clin[f"{name}_raw"] = clin["sample"].map(counts.loc[ensg])
        clin[f"{name}_CPM"] = clin["sample"].map(cpm.loc[ensg])
        clin[f"{name}_log2CPM1"] = np.log2(clin[f"{name}_CPM"] + 1)

    z_tac = stats.zscore(clin["TACSTD2_log2CPM1"], ddof=1)
    z_cld = stats.zscore(clin["CLDN4_log2CPM1"], ddof=1)
    clin["combined_z"] = (z_tac + z_cld) / 2.0

    clin = clin.sort_values(["response_ordinal", "sample"], ascending=[False, True])
    clin.to_csv(OUT_DIR / "sample_data.csv", index=False)

    features = [(n, f"{n}_log2CPM1") for n in PRIMARY] + [("combined_z", "combined_z")]
    contrasts = [
        ("is_CR", "CR", "PD", clin["response"].isin(["CR", "PD"])),
        ("is_CR", "CR", "nonCR", pd.Series(True, index=clin.index)),
        ("is_DCR", "DCR", "PD", pd.Series(True, index=clin.index)),
        ("pla_high", "PLA_high", "low", pd.Series(True, index=clin.index)),
    ]
    clin = clin.copy()
    clin["pla_high"] = (clin["pla"] == "high").astype(int)

    rows = []
    for label_col, pos_name, neg_name, mask in contrasts:
        sub = clin.loc[mask].copy()
        for feat_name, col in features:
            rec = contrast_stats(sub, col, label_col, pos_name, neg_name)
            rec["feature"] = feat_name
            rec["family"] = "exploratory" if feat_name == "combined_z" or pos_name.startswith("PLA") else "primary_or_sensitivity"
            rows.append(rec)
    stats_df = pd.DataFrame(rows)

    primary_mask = (stats_df.contrast == "CR_vs_PD") & (stats_df.feature.isin(PRIMARY))
    stats_df.loc[primary_mask, "q_BH_2genes"] = bh_adjust(
        stats_df.loc[primary_mask, "p_exact_two_sided"].tolist()
    )
    stats_df.to_csv(OUT_DIR / "stats.csv", index=False)

    kw_rows = []
    spearman_rows = []
    for feat_name, col in features:
        groups = [clin.loc[clin.response == r, col].to_numpy() for r in ("CR", "SD", "PD")]
        kw = stats.kruskal(*groups)
        kw_rows.append({"feature": feat_name, "test": "kruskal_CR_SD_PD", "p": float(kw.pvalue), "stat": float(kw.statistic)})
        r, p = stats.spearmanr(clin[col], clin["response_ordinal"])
        spearman_rows.append({"feature": feat_name, "test": "spearman_ordinal_CR2_SD1_PD0", "spearman_r": float(r), "p": float(p)})
    pd.DataFrame(kw_rows + [{**x, "stat": x["spearman_r"]} for x in spearman_rows]).to_csv(
        OUT_DIR / "kruskal_spearman.csv", index=False
    )

    corr = stats.spearmanr(clin["TACSTD2_log2CPM1"], clin["CLDN4_log2CPM1"])
    pd.DataFrame(
        [{"comparison": "TACSTD2_vs_CLDN4", "spearman_rho": corr.statistic, "p_two_sided": corr.pvalue}]
    ).to_csv(OUT_DIR / "gene_correlation.csv", index=False)

    ctx_rows = []
    sub = clin[clin.response.isin(["CR", "PD"])]
    for name in CONTEXT:
        ctx_rows.append(contrast_stats(sub, f"{name}_log2CPM1", "is_CR", "CR", "PD"))
    pd.DataFrame(ctx_rows).to_csv(OUT_DIR / "context_stats.csv", index=False)

    hgnc = load_hgnc(hgnc_path)
    protein_ids = set(hgnc.loc[hgnc["locus_group"] == "protein-coding gene", "ensembl_gene_id"])
    # strip version if any
    protein_ids = {re.sub(r"\.\d+$", "", x) for x in protein_ids if isinstance(x, str)}

    cpm_aligned = cpm[clin["sample"].tolist()]
    gw_crpd = genome_wide_rank(
        cpm_aligned.loc[:, clin.response.isin(["CR", "PD"])],
        clin.loc[clin.response.isin(["CR", "PD"]), "is_CR"].to_numpy(),
        protein_ids,
    )
    gw_cr = genome_wide_rank(cpm_aligned, clin["is_CR"].to_numpy(), protein_ids)
    gw_crpd.to_csv(OUT_DIR / "genomewide_rank_CR_vs_PD.csv", index=False)
    gw_cr.to_csv(OUT_DIR / "genomewide_rank_CR_vs_nonCR.csv", index=False)

    # ---- plots ----
    rng = np.random.default_rng(SEED)
    colors_r = {"CR": "#2ca02c", "SD": "#ff7f0e", "PD": "#d62728"}
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2))
    for ax, (feat_name, col) in zip(axes, features):
        groups = [clin.loc[clin.response == r, col].to_numpy() for r in ("CR", "SD", "PD")]
        box_strip(
            ax,
            groups,
            [f"{r}\n(n={(clin.response == r).sum()})" for r in ("CR", "SD", "PD")],
            [colors_r[r] for r in ("CR", "SD", "PD")],
            rng,
        )
        kw = next(x for x in kw_rows if x["feature"] == feat_name)
        ax.set_title(f"{feat_name}\nKruskal p={kw['p']:.3f}")
        ax.set_ylabel("log2(CPM+1)" if col != "combined_z" else "mean z-score")
    fig.suptitle("GSE283829 leftover (NSCLC ICI): TACSTD2 / CLDN4 vs RECIST-like response")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "boxplots_response.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2))
    sub = clin[clin.response.isin(["CR", "PD"])]
    for ax, (feat_name, col) in zip(axes, features):
        rec = stats_df[(stats_df.feature == feat_name) & (stats_df.contrast == "CR_vs_PD")].iloc[0]
        groups = [sub.loc[sub.response == r, col].to_numpy() for r in ("PD", "CR")]
        box_strip(ax, groups, [f"PD\n(n={int(rec.n_neg)})", f"CR\n(n={int(rec.n_pos)})"],
                  ["#d62728", "#2ca02c"], rng)
        ax.set_title(f"{feat_name}\nMWU p={rec.p_exact_two_sided:.3f}, AUC={rec.AUC_high_predicts_pos:.2f}")
        ax.set_ylabel("log2(CPM+1)" if col != "combined_z" else "mean z-score")
    fig.suptitle("GSE283829 leftover: primary contrast CR vs PD")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "boxplots_cr_vs_pd.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    for r, c in colors_r.items():
        s = clin[clin.response == r]
        ax.scatter(s["TACSTD2_log2CPM1"], s["CLDN4_log2CPM1"], s=42, color=c, label=f"{r} (n={len(s)})",
                   edgecolors="k", linewidths=0.3, alpha=0.9)
    ax.set_xlabel("TACSTD2 log2(CPM+1)")
    ax.set_ylabel("CLDN4 log2(CPM+1)")
    ax.legend(frameon=False)
    ax.set_title(f"TACSTD2 vs CLDN4 (Spearman ρ={corr.statistic:.2f})")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scatter_tacstd2_cldn4.png", dpi=200)
    plt.close(fig)

    write_readme(clin, stats_df, gw_crpd, gw_cr, spearman_rows, kw_rows, corr)

    import scipy

    manifest = [
        f"python={platform.python_version()}",
        f"pandas={pd.__version__}",
        f"numpy={np.__version__}",
        f"scipy={scipy.__version__}",
        f"matplotlib={matplotlib.__version__}",
        f"random_seed={SEED}",
        f"bootstrap_replicates={N_BOOT}",
        f"sha256  {sha256(counts_path)}  data/GSE283829/{counts_path.name}",
        f"sha256  {sha256(samples_path)}  data/GSE283829/{samples_path.name}",
        f"sha256  {sha256(series_path)}  data/GSE283829/{series_path.name}",
        f"sha256  {sha256(hgnc_path)}  data/GSE283829/{hgnc_path.name}",
    ]
    (OUT_DIR / "analysis_manifest.txt").write_text("\n".join(manifest) + "\n")

    print(stats_df.to_string(index=False))
    print()
    print(pd.DataFrame(kw_rows).to_string(index=False))
    print(pd.DataFrame(spearman_rows).to_string(index=False))
    for table, name in ((gw_crpd, "CR_vs_PD"), (gw_cr, "CR_vs_nonCR")):
        for gene, ensg in PRIMARY.items():
            hit = table[table.ensembl_gene_id == ensg]
            if hit.empty:
                print(f"{name} {gene}: not in filtered set")
            else:
                r = hit.iloc[0]
                print(f"{name} {gene}: rank {int(r.rank_by_p)}/{int(r.n_tested)} p={r.p:.4g}")


if __name__ == "__main__":
    main()
