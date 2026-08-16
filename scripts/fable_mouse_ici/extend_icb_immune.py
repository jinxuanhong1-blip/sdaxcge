#!/usr/bin/env python3
"""Hunt follow-up: extra paired ICB sets + TROP2-high/immune-low + combined ICB stats.

Adds three new processed <2 GB tumour RNA series with a control-vs-ICB arm
(GSE114601, GSE157880, GSE309199) plus the GSE169196 total-tumour IgG vs
A2V+aPD1 arm. Extracts Tacstd2/Cldn4 and a small cytotoxic/IFN immune panel,
tests Tacstd2 vs immune-score correlation (TROP2-high / immune-low), and
combines independent ICB-vs-control Tacstd2 contrasts with a sign test and
Stouffer Z (honest: one primary contrast per dataset; no fabricated p-values).
"""
import csv
import gzip
import io
import json
import tarfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/tmp/fable_ici_data")
RESULTS = Path(__file__).resolve().parents[2] / "results" / "fable_mouse_ici"
NOTES = Path(__file__).resolve().parents[2] / "notes" / "fable_mouse_ici"
FIG = RESULTS / "figures"
FIG.mkdir(parents=True, exist_ok=True)

ENSEMBL = {
    "Tacstd2": "ENSMUSG00000051397",
    "Cldn4": "ENSMUSG00000041378",
    "Cd8a": "ENSMUSG00000053044",
    "Cd3e": "ENSMUSG00000032093",
    "Cd3d": "ENSMUSG00000032094",
    "Gzmb": "ENSMUSG00000015437",
    "Prf1": "ENSMUSG00000037202",
    "Ifng": "ENSMUSG00000055170",
    "Cd274": "ENSMUSG00000016496",
    "Pdcd1": "ENSMUSG00000026285",
    "Cxcl9": "ENSMUSG00000029417",
    "Cxcl10": "ENSMUSG00000034855",
    "Nkg7": "ENSMUSG00000036363",
}
IMMUNE = ["Cd8a", "Cd3e", "Cd3d", "Gzmb", "Prf1", "Ifng", "Cd274",
          "Pdcd1", "Cxcl9", "Cxcl10", "Nkg7"]
TARGETS = ["Tacstd2", "Cldn4"]
ALIASES = {g.lower(): g for g in list(ENSEMBL)}
ALIASES["trop2"] = "Tacstd2"
ALIASES["pd-l1"] = "Cd274"
ALIASES["pdl1"] = "Cd274"
ALIASES["pd1"] = "Pdcd1"

LONG = []  # dataset, gene, sample, group, is_control, value_log2, value_type, is_icb_arm


def add(ds, gene, sample, group, is_control, val, vtype, is_icb_arm=False):
    LONG.append(dict(dataset=ds, gene=gene, sample=str(sample), group=group,
                     is_control=bool(is_control), value_log2=float(val),
                     value_type=vtype, is_icb_arm=bool(is_icb_arm)))


def match_sym(name):
    if name is None:
        return None
    return ALIASES.get(str(name).strip().lower())


def log2p1(x):
    return np.log2(np.asarray(x, float) + 1.0)


def cpm_log(counts):
    counts = np.asarray(counts, float)
    lib = counts.sum(axis=0)
    lib[lib == 0] = 1.0
    return np.log2(counts / lib * 1e6 + 1.0)


# ---------------------------------------------------------------------------
# GSE114601 — GEMM NSCLC, DESeq-normalised counts, n=2 / arm
# ---------------------------------------------------------------------------
def load_gse114601():
    ds = "GSE114601"
    grp = {
        "s1795": ("anti-PD1", False, True),
        "s2521": ("anti-PD1", False, True),
        "s2596": ("Vehicle", True, False),
        "s2617": ("Vehicle", True, False),
        "s2503": ("JQ1", False, False),
        "s2625": ("JQ1", False, False),
        "s2467": ("anti-PD1+JQ1", False, True),
        "s2669": ("anti-PD1+JQ1", False, True),
    }
    with gzip.open(DATA / "GSE114601_counts.normalized.csv.gz", "rt") as f:
        df = pd.read_csv(f, index_col=0)
    want = {g: None for g in TARGETS + IMMUNE}
    for idx in df.index:
        g = match_sym(idx)
        if g in want:
            want[g] = df.loc[idx]
    for g, ser in want.items():
        if ser is None:
            continue
        for col, val in ser.items():
            gname, isc, icb = grp[col]
            add(ds, g, col, gname, isc, log2p1(val), "log2(norm+1)", icb)


# ---------------------------------------------------------------------------
# GSE157880 — KP lung, IgG vs anti-PD-1 ± 0/4/8 Gy (already abundance-like)
# ---------------------------------------------------------------------------
def load_gse157880():
    ds = "GSE157880"
    # title map from series matrix
    colgrp = {
        "0-1_S13": ("IgG 0Gy", True, False),
        "0-2_S14": ("IgG 0Gy", True, False),
        "0-3_S15": ("IgG 0Gy", True, False),
        "4_2_S26": ("IgG 4Gy", False, False),
        "4_3_S27": ("IgG 4Gy", False, False),
        "8_1_S31": ("IgG 8Gy", False, False),
        "8_2_S32": ("IgG 8Gy", False, False),
        "8_3_S33": ("IgG 8Gy", False, False),
        "0-4_S16": ("PD-1 0Gy", False, True),
        "0-5_S17": ("PD-1 0Gy", False, True),
        "4_4_S28": ("PD-1 4Gy", False, True),
        "4_5_S29": ("PD-1 4Gy", False, True),
        "4_6_S30": ("PD-1 4Gy", False, True),
        "8_4_S34": ("PD-1 8Gy", False, True),
        "8_5_S35": ("PD-1 8Gy", False, True),
        "8_6_S36": ("PD-1 8Gy", False, True),
    }
    with gzip.open(DATA / "GSE157880_Bulk048.txt.gz", "rt") as f:
        df = pd.read_csv(f, sep="\t")
    # gene symbol column
    sc = "Gene Symbol" if "Gene Symbol" in df.columns else "gene_name"
    df = df.set_index(sc)
    for idx in df.index:
        g = match_sym(idx)
        if g not in TARGETS + IMMUNE:
            continue
        for col, (gname, isc, icb) in colgrp.items():
            if col not in df.columns:
                continue
            add(ds, g, col, gname, isc, log2p1(df.loc[idx, col]),
                "log2(abundance+1)", icb)


# ---------------------------------------------------------------------------
# GSE309199 — RPM SCLC tumours, Ctrl / aPD1 / ent / combo (n=3)
# ---------------------------------------------------------------------------
def load_gse309199():
    ds = "GSE309199"
    # column name in matrix -> (group, is_control, is_icb)
    colgrp = {
        "1_296701_S2": ("Ctrl", True, False),
        "2_299538_S3": ("Ctrl", True, False),
        "3_300334_S4": ("Ctrl", True, False),
        "4_295947_S5": ("aPD1", False, True),
        "5_295948_S1": ("aPD1", False, True),
        "6_299534_S6": ("aPD1", False, True),
        "7_295936_S7": ("entinostat", False, False),
        "8_295932_S8": ("entinostat", False, False),
        "9_295951_S9": ("entinostat", False, False),
        "10_295933_S10": ("aPD1+entinostat", False, True),
        "11_296689_S11": ("aPD1+entinostat", False, True),
        "12_298090_S12": ("aPD1+entinostat", False, True),
    }
    with gzip.open(DATA / "GSE309199_raw_counts.tsv.gz", "rt") as f:
        df = pd.read_csv(f, sep="\t", index_col=0)
    # CPM from all genes
    log = pd.DataFrame(cpm_log(df.values), index=df.index, columns=df.columns)
    for idx in log.index:
        g = match_sym(idx)
        if g not in TARGETS + IMMUNE:
            continue
        for col, (gname, isc, icb) in colgrp.items():
            add(ds, g, col, gname, isc, log.loc[idx, col], "log2(CPM+1)", icb)


# ---------------------------------------------------------------------------
# GSE169196 — KPM lung, total viable cells only (IgG / A2V / A2V+aPD1)
# ---------------------------------------------------------------------------
def load_gse169196():
    ds = "GSE169196"
    totals = {
        "GSM5182757": ("IgG", True, False),
        "GSM5182761": ("IgG", True, False),
        "GSM5182765": ("IgG", True, False),
        "GSM5182768": ("A2V", False, False),
        "GSM5182772": ("A2V", False, False),
        "GSM5182776": ("A2V", False, False),
        "GSM5182783": ("A2V+aPD1", False, True),
        "GSM5182787": ("A2V+aPD1", False, True),
        "GSM5182791": ("A2V+aPD1", False, True),
    }
    ens2g = {v: k for k, v in ENSEMBL.items()}
    with tarfile.open(DATA / "GSE169196_RAW.tar") as t:
        for m in t.getmembers():
            gsm = m.name.split("_")[0]
            if gsm not in totals:
                continue
            raw = gzip.decompress(t.extractfile(m).read()).decode()
            counts = {}
            lib = 0.0
            for line in raw.splitlines():
                p = line.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                ens = p[0].split(".")[0]
                try:
                    v = float(p[1])
                except ValueError:
                    continue
                lib += v
                if ens in ens2g:
                    counts[ens2g[ens]] = v
            for g, v in counts.items():
                cpm = v / lib * 1e6 if lib else 0.0
                gname, isc, icb = totals[gsm]
                add(ds, g, gsm, gname, isc, np.log2(cpm + 1.0),
                    "log2(CPM+1)", icb)


# ---------------------------------------------------------------------------
# Immune panel from already-cached bulk matrices (GSE239485, E-MTAB, GSE197260)
# ---------------------------------------------------------------------------
def load_gse239485_panel():
    ds = "GSE239485"
    grp = {"C": ("Control vehicle", True, False),
           "D": ("PolyIC+anti-PD1", False, True),
           "T": ("PolyIC+anti-PD1+anti-C5aR1", False, True)}
    import openpyxl
    wb = openpyxl.load_workbook(DATA / "GSE239485_Processed_data.xlsx",
                                read_only=True, data_only=True)
    ws = wb["DataNorm"]
    rows = ws.iter_rows(values_only=True)
    hdr = next(rows)
    gi = hdr.index("gene_name")
    scols = list(range(11, len(hdr)))
    for row in rows:
        g = match_sym(row[gi])
        if g not in TARGETS + IMMUNE:
            continue
        for ci in scols:
            col = hdr[ci]
            if col is None:
                continue
            pref = str(col).split("_")[0]
            if pref not in grp:
                continue
            gname, isc, icb = grp[pref]
            add(ds, g, str(col), gname, isc, float(row[ci]),
                "log2_normalised", icb)
    wb.close()


def load_emtab_panel():
    ds = "E-MTAB-13704"
    r2grp = {}
    with open(DATA / "E-MTAB-13704.sdrf.txt") as f:
        rd = list(csv.reader(f, delimiter="\t"))
    h = rd[0]
    scan_i = h.index("Scan Name")
    fac_i = [i for i, c in enumerate(h) if c.startswith("Factor Value")][0]
    for row in rd[1:]:
        r2grp[row[scan_i].split("_")[0]] = row[fac_i].strip()
    icb_grps = {"aPD-L1", "ATRi/aPD-L1", "Cisplatin/aPD-L1 /aCTLA4",
                "VEGFRi/aPD-L1"}
    with open(DATA / "E-MTAB-13704_GEMMS_raw_counts.csv") as f:
        rd = csv.reader(f)
        cols = next(rd)[1:]
        mat = []
        genes = []
        for row in rd:
            ens = row[0].split(".")[0]
            g = {v: k for k, v in ENSEMBL.items()}.get(ens)
            if g:
                genes.append(g)
                mat.append([float(x) for x in row[1:]])
            else:
                # still need libsize — accumulate all
                pass
    # libsize from full file
    with open(DATA / "E-MTAB-13704_GEMMS_raw_counts.csv") as f:
        rd = csv.reader(f)
        next(rd)
        lib = np.zeros(len(cols))
        keep = {}
        ens2g = {v: k for k, v in ENSEMBL.items()}
        for row in rd:
            ens = row[0].split(".")[0]
            vals = np.array([float(x) for x in row[1:]])
            lib += vals
            if ens in ens2g:
                keep[ens2g[ens]] = vals
    for g, vals in keep.items():
        log = np.log2(vals / lib * 1e6 + 1.0)
        for j, c in enumerate(cols):
            grp = r2grp.get(c, "?")
            add(ds, g, c, grp, grp == "vehicle", log[j], "log2(CPM+1)",
                grp in icb_grps)


def load_gse197260_panel():
    ds = "GSE197260"
    labels = {
        "vehicle_d3": ("vehicle", True, False),
        "gef_d3": ("gefitinib d3", False, False),
        "gef_d14": ("gefitinib d14", False, False),
        "gef_vehicle_d21": ("gefitinib+vehicle", True, False),
        "gef_dc101_d21": ("gefitinib+anti-VEGFR2", False, False),
        "gef_4h2_d21": ("gefitinib+anti-PD-1", False, True),
        "gef_comb_d21": ("gefitinib+anti-PD-1+anti-VEGFR2", False, True),
    }
    with gzip.open(DATA / "GSE197260_TPM.txt.gz", "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            g = match_sym(p[0])
            if g not in TARGETS + IMMUNE:
                continue
            for ci in range(1, len(hdr)):
                col = hdr[ci]
                gname, isc, icb = labels[col]
                add(ds, g, col, gname, isc, log2p1(float(p[ci])),
                    "log2(TPM+1)", icb)


def welch_onesided_up(treat, ctrl):
    """Two-sided Welch + one-sided p for treat > ctrl."""
    t = stats.ttest_ind(treat, ctrl, equal_var=False)
    # one-sided: P(T > t_obs) if we want treat>ctrl, using the t statistic
    # scipy ttest_ind statistic is (mean_treat - mean_ctrl) / se when we pass treat, ctrl
    p_two = float(t.pvalue)
    # one-sided from t.sf
    df = t.df
    p_up = float(stats.t.sf(t.statistic, df))  # P(T >= t_obs)
    return float(t.statistic), p_two, p_up


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    q = np.full_like(p, np.nan)
    idx = np.where(~np.isnan(p))[0]
    if len(idx) == 0:
        return q
    pv = p[idx]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(pv)
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    q[idx] = out
    return q


def contrasts(df):
    """All treat-vs-control contrasts with n>=2 both sides, for Tacstd2/Cldn4."""
    rows = []
    for ds, dsub in df.groupby("dataset"):
        for gene, gsub in dsub.groupby("gene"):
            if gene not in TARGETS:
                continue
            ctrls = gsub[gsub.is_control]
            if ctrls.empty:
                continue
            # if multiple control labels, keep the one with most samples
            ctrl_lab = ctrls.groupby("group").size().idxmax()
            cvals = ctrls.loc[ctrls.group == ctrl_lab, "value_log2"].to_numpy()
            for grp, arm in gsub.groupby("group"):
                if grp == ctrl_lab:
                    continue
                tvals = arm["value_log2"].to_numpy()
                rec = dict(dataset=ds, gene=gene, comparison=f"{grp} vs {ctrl_lab}",
                           n_control=len(cvals), n_treat=len(tvals),
                           mean_control=round(float(np.mean(cvals)), 4),
                           mean_treat=round(float(np.mean(tvals)), 4),
                           log2FC=round(float(np.mean(tvals) - np.mean(cvals)), 4),
                           is_icb_arm=bool(arm.is_icb_arm.iloc[0]),
                           value_type=arm.value_type.iloc[0])
                if len(cvals) >= 2 and len(tvals) >= 2:
                    tstat, p2, pup = welch_onesided_up(tvals, cvals)
                    rec["welch_t"] = round(tstat, 4)
                    rec["welch_p_two"] = p2
                    rec["welch_p_up"] = pup
                    try:
                        rec["mwu_p"] = float(stats.mannwhitneyu(
                            tvals, cvals, alternative="two-sided").pvalue)
                    except ValueError:
                        rec["mwu_p"] = np.nan
                else:
                    rec["welch_t"] = rec["welch_p_two"] = rec["welch_p_up"] = rec["mwu_p"] = np.nan
                rows.append(rec)
    return pd.DataFrame(rows)


# One primary ICB-vs-control contrast per independent tumour RNA dataset
PRIMARY = [
    ("GSE239485", "PolyIC+anti-PD1 vs Control vehicle"),
    ("GSE297630", "anti-PD-1 vs Control"),
    ("E-MTAB-13704", "aPD-L1 vs vehicle"),
    ("GSE114601", "anti-PD1 vs Vehicle"),
    ("GSE157880", "PD-1 0Gy vs IgG 0Gy"),
    ("GSE309199", "aPD1 vs Ctrl"),
    ("GSE169196", "A2V+aPD1 vs IgG"),
]


def combine_primary(new_stats, old_stats):
    """Merge old Tacstd2/Cldn4 stats (GSE297630 already computed) with new."""
    # old_stats has log2FC_treat_minus_control, welch_p
    frames = []
    if new_stats is not None and len(new_stats):
        frames.append(new_stats)
    if old_stats is not None and len(old_stats):
        o = old_stats.rename(columns={
            "log2FC_treat_minus_control": "log2FC",
            "welch_p": "welch_p_two",
            "welch_t_stat": "welch_t",
        })
        # recompute one-sided from t if present
        if "welch_t" in o.columns:
            # df approx n1+n2-2
            df_approx = o["n_control"] + o["n_treat"] - 2
            o["welch_p_up"] = [
                float(stats.t.sf(t, d)) if pd.notna(t) else np.nan
                for t, d in zip(o["welch_t"], df_approx)
            ]
        o["is_icb_arm"] = True
        frames.append(o[["dataset", "gene", "comparison", "n_control", "n_treat",
                         "mean_control", "mean_treat", "log2FC", "welch_t",
                         "welch_p_two", "welch_p_up", "mwu_p"]])
    alls = pd.concat(frames, ignore_index=True)
    # drop exact duplicate dataset+gene+comparison keeping first (new preferred)
    alls = alls.drop_duplicates(["dataset", "gene", "comparison"], keep="first")
    return alls


def stouffer(p_up):
    p = np.asarray(p_up, float)
    p = np.clip(p, 1e-15, 1 - 1e-15)
    z = stats.norm.isf(p)  # one-sided z
    zc = z.sum() / np.sqrt(len(z))
    return float(zc), float(stats.norm.sf(zc))


def immune_analysis(df):
    """Per-dataset Spearman Tacstd2 vs immune score; median-split immune test."""
    rows = []
    split_rows = []
    bulk = df[df.dataset.isin([
        "GSE239485", "E-MTAB-13704", "GSE114601", "GSE157880",
        "GSE309199", "GSE169196", "GSE197260",
    ])]
    for ds, dsub in bulk.groupby("dataset"):
        wide = dsub.pivot_table(index=["sample", "group", "is_control"],
                                columns="gene", values="value_log2",
                                aggfunc="mean")
        if "Tacstd2" not in wide.columns:
            continue
        present = [g for g in IMMUNE if g in wide.columns]
        if len(present) < 3:
            continue
        score = wide[present].mean(axis=1)
        tac = wide["Tacstd2"]
        # all samples
        rho, p = stats.spearmanr(tac, score)
        rows.append(dict(dataset=ds, subset="all_samples", n=len(tac),
                         n_immune_genes=len(present),
                         immune_genes=",".join(present),
                         spearman_rho=float(rho), spearman_p=float(p)))
        # control-only (baseline TROP2-high / immune-low)
        ctrl_idx = [i for i in wide.index if i[2] is True or i[2] == True]
        if len(ctrl_idx) >= 4:
            rho_c, p_c = stats.spearmanr(tac.loc[ctrl_idx], score.loc[ctrl_idx])
            rows.append(dict(dataset=ds, subset="control_only", n=len(ctrl_idx),
                             n_immune_genes=len(present),
                             immune_genes=",".join(present),
                             spearman_rho=float(rho_c), spearman_p=float(p_c)))
        # median-split Tacstd2-high vs low on ALL samples
        med = float(tac.median())
        hi = score[tac >= med]
        lo = score[tac < med]
        if len(hi) >= 2 and len(lo) >= 2:
            t = stats.ttest_ind(hi, lo, equal_var=False)
            split_rows.append(dict(
                dataset=ds, n_high=len(hi), n_low=len(lo),
                tacstd2_median=round(med, 4),
                immune_mean_Tacstd2_high=round(float(hi.mean()), 4),
                immune_mean_Tacstd2_low=round(float(lo.mean()), 4),
                immune_delta_high_minus_low=round(float(hi.mean() - lo.mean()), 4),
                welch_p=float(t.pvalue),
                trop2_high_immune_low=bool(hi.mean() < lo.mean()),
            ))
    return pd.DataFrame(rows), pd.DataFrame(split_rows)


def plot_primary(prim):
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    y = np.arange(len(prim))
    colors = ["#c0392b" if (pd.notna(p) and p < 0.05) else "#7f8c8d"
              for p in prim.welch_p_two]
    ax.axvline(0, color="k", lw=0.8)
    ax.scatter(prim.log2FC, y, c=colors, s=55, zorder=3)
    for yi, (_, r) in enumerate(prim.iterrows()):
        lab = f"p={r.welch_p_two:.3g}" if pd.notna(r.welch_p_two) else "n<2"
        ax.text(r.log2FC, yi + 0.18, lab, fontsize=7, ha="center")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.dataset} | {r.comparison}" for _, r in prim.iterrows()],
                       fontsize=8)
    ax.set_xlabel("Tacstd2 log2FC (ICB − control)")
    ax.set_title("Primary ICB-vs-control Tacstd2 contrasts (red = two-sided Welch p<0.05)")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "primary_icb_tacstd2.png", dpi=140)
    plt.close(fig)


def plot_immune(corr, split):
    if corr.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    sub = corr[corr.subset == "all_samples"].copy()
    y = np.arange(len(sub))
    colors = ["#1a5276" if r < 0 else "#b9770e" for r in sub.spearman_rho]
    ax.axvline(0, color="k", lw=0.8)
    ax.scatter(sub.spearman_rho, y, c=colors, s=50, zorder=3)
    for yi, (_, r) in enumerate(sub.iterrows()):
        ax.text(r.spearman_rho, yi + 0.18, f"p={r.spearman_p:.3g} n={r.n}",
                fontsize=7, ha="center")
    ax.set_yticks(y)
    ax.set_yticklabels(sub.dataset, fontsize=9)
    ax.set_xlabel("Spearman ρ  Tacstd2 vs immune score")
    ax.set_title("TROP2-high / immune-low?  (negative ρ supports the subset)")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "tacstd2_vs_immune.png", dpi=140)
    plt.close(fig)


def main():
    print("loading new + panel ...")
    load_gse114601()
    load_gse157880()
    load_gse309199()
    load_gse169196()
    load_gse239485_panel()
    load_emtab_panel()
    load_gse197260_panel()
    df = pd.DataFrame(LONG)
    df.to_csv(RESULTS / "extended_per_sample.csv", index=False)
    print("extended rows", len(df), df.groupby(["dataset", "gene"]).size().head(20))

    st = contrasts(df)
    old = pd.read_csv(RESULTS / "stats_results.csv")
    # keep GSE297630 (array, not re-extracted here) and GSE241978/GSE330658 as non-primary
    merged = combine_primary(st, old)
    merged.to_csv(RESULTS / "extended_stats.csv", index=False)

    # primary Tacstd2 table
    prim_rows = []
    for ds, comp in PRIMARY:
        hit = merged[(merged.dataset == ds) & (merged.gene == "Tacstd2")
                     & (merged.comparison == comp)]
        if hit.empty:
            prim_rows.append(dict(dataset=ds, comparison=comp, note="missing"))
        else:
            prim_rows.append(hit.iloc[0].to_dict())
    prim = pd.DataFrame(prim_rows)
    prim.to_csv(RESULTS / "primary_icb_tacstd2.csv", index=False)

    usable = prim.dropna(subset=["log2FC"])
    n_up = int((usable.log2FC > 0).sum())
    n = int(len(usable))
    sign_p = float(stats.binomtest(n_up, n, 0.5, alternative="greater").pvalue)
    # Stouffer on one-sided up p (only rows with a real p)
    with_p = usable.dropna(subset=["welch_p_up"])
    if len(with_p):
        zc, pz = stouffer(with_p.welch_p_up.to_numpy())
    else:
        zc, pz = np.nan, np.nan
    combo = dict(
        n_primary_contrasts=n,
        n_tacstd2_up=n_up,
        sign_test_p_greater=sign_p,
        stouffer_Z=zc,
        stouffer_p_up=pz,
        datasets=", ".join(usable.dataset.tolist()),
        note=("One primary ICB-vs-control contrast per independent tumour RNA "
              "dataset. GSE239485 ICB arm is Poly I:C + anti-PD-1 (Poly I:C "
              "confounder). GSE169196 has no aPD1-monotherapy total-tumour arm "
              "(A2V+aPD1 vs IgG). TISMO 49/64 p=5.8e-5 is an external prior; "
              "the TISMO expression matrix was not independently re-downloadable "
              "from the current site (SPA + github scripts only)."),
    )
    (RESULTS / "combined_icb_tacstd2.json").write_text(json.dumps(combo, indent=2))
    print("COMBINED", combo)

    corr, split = immune_analysis(df)
    corr.to_csv(RESULTS / "tacstd2_immune_correlation.csv", index=False)
    split.to_csv(RESULTS / "trop2_high_immune_split.csv", index=False)
    print("=== immune corr ===")
    print(corr.to_string(index=False))
    print("=== median split ===")
    print(split.to_string(index=False))
    print("=== primary Tacstd2 ===")
    print(prim.to_string(index=False))

    plot_primary(usable)
    plot_immune(corr, split)

    # hunt triage note
    hunt = json.loads((NOTES / "hunt_candidates.json").read_text())
    skip = {
        "GSE190264": "in-vitro LLC1 chemo/MEK, no ICB antibody",
        "GSE114300": "sorted CD4 T cells, not tumour epithelium",
        "GSE277610": "sorted CD8 T cells",
        "GSE309192": "control vs tumour tissue, no ICB",
        "GSE184000": "whole-lung irAE, not tumour",
        "GSE315010": "RAS(ON) inhibitors, no ICB antibody arm in the mouse matrix",
        "GSE303940": "CMT167R PKC inhibitor, n=1, no ICB",
        "GSE194166": "scRNA/TCR, 238 MB, n=1/arm",
        "GSE129298": "scRNA n=1/arm (companion of already-included GSE129297)",
        "GSE157883": "parent SuperSeries of GSE157880 (used the processed Bulk048)",
    }
    (NOTES / "hunt_triage.md").write_text(
        "# Extra accession hunt (2026-08-16)\n\n"
        "Queried NCBI GEO for mouse lung + anti-PD-1/PD-L1 RNA series. "
        "Processed files <2 GB only.\n\n"
        "## Included in this extension\n"
        "- GSE114601 GEMM NSCLC Vehicle vs anti-PD1 (n=2)\n"
        "- GSE157880 KP lung IgG vs PD-1 ± RT (0 Gy pair is clean ICB)\n"
        "- GSE309199 RPM SCLC Ctrl vs aPD1 (n=3)\n"
        "- GSE169196 KPM total-tumour IgG vs A2V+aPD1 (no aPD1-mono arm)\n\n"
        "## Explicitly skipped (honest)\n"
        + "\n".join(f"- {k}: {v}" for k, v in skip.items())
        + "\n"
    )


if __name__ == "__main__":
    main()
