#!/usr/bin/env python3
"""Micronuclei / cytosolic DNA / cGAS GEO: claudin and tight-junction panel.

Public matrices only. One contrast per accession (no mega-merge).
Counts stay in /tmp/geo_mn. This script writes tables under results/.

Direct paper link (summary names cGAS/STING or cytosolic DNA AND claudin
or tight junctions): GSE279172. GSE296527 reuses that high-density WT
library and adds IFNL1 KO.

Other accessions are landmark micronucleus / cytosolic-DNA / cGAS series.
Claudin and tight-junction genes are tested there; those papers do not
claim a claudin mechanism. n=1 matrices are descriptive folds only.
"""

from __future__ import annotations

import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "mn_cydna_cgas_tj"
CACHE = Path("/tmp/geo_mn")
OUT.mkdir(parents=True, exist_ok=True)

HUMAN_PANEL = {
    "claudin": [
        "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN6", "CLDN7",
        "CLDN8", "CLDN9", "CLDN10", "CLDN11", "CLDN12", "CLDN14",
        "CLDN15", "CLDN16", "CLDN18", "CLDN23",
    ],
    "surface": ["TACSTD2"],
    "tight_junction": [
        "TJP1", "TJP2", "TJP3", "OCLN", "CDH1", "F11R", "CGN", "MARVELD2",
    ],
    "sensor": ["CGAS", "MB21D1", "STING1", "TMEM173", "TBK1", "IRF3"],
    "ifnl": ["IFNL1", "IFNL2", "IFNL3", "IFNLR1", "IL10RB"],
    "isg_control": ["ISG15", "IFIT1", "MX1", "OAS1", "CXCL10", "CCL5", "CXCL8"],
    "autophagy_aa": ["MAP1LC3B", "SQSTM1", "BECN1", "SLC1A5", "SLC7A5"],
}

MOUSE_PANEL = {
    "claudin": [
        "Cldn1", "Cldn2", "Cldn3", "Cldn4", "Cldn5", "Cldn6", "Cldn7",
        "Cldn8", "Cldn9", "Cldn10", "Cldn11", "Cldn12", "Cldn14",
        "Cldn15", "Cldn16", "Cldn18", "Cldn23",
    ],
    "surface": ["Tacstd2"],
    "tight_junction": [
        "Tjp1", "Tjp2", "Tjp3", "Ocln", "Cdh1", "F11r", "Cgn", "Marveld2",
    ],
    "sensor": ["Cgas", "Mb21d1", "Sting1", "Tmem173", "Tbk1", "Irf3"],
    "ifnl": ["Ifnl2", "Ifnl3", "Ifnlr1", "Il10rb"],
    "isg_control": ["Isg15", "Ifit1", "Mx1", "Oas1", "Cxcl10", "Ccl5"],
    "autophagy_aa": ["Map1lc3b", "Sqstm1", "Becn1", "Slc1a5", "Slc7a5"],
}


def panel_lookup(organism: str) -> dict[str, str]:
    src = HUMAN_PANEL if organism == "human" else MOUSE_PANEL
    out = {}
    for cls, genes in src.items():
        for g in genes:
            out[g] = cls
    return out


def bh(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, float)
    out = np.full(arr.shape, np.nan)
    ok = np.isfinite(arr)
    if ok.sum() == 0:
        return out.tolist()
    _, q, _, _ = multipletests(arr[ok], method="fdr_bh")
    out[ok] = q
    return out.tolist()


def direction(log2fc: float, q: float, low: bool) -> str:
    if low:
        return "LOW"
    if not np.isfinite(log2fc):
        return "NA"
    if np.isfinite(q) and q < 0.05 and log2fc >= 0.5:
        return "UP"
    if np.isfinite(q) and q < 0.05 and log2fc <= -0.5:
        return "DOWN"
    if np.isfinite(q) and q < 0.05:
        return "SIG_SMALL_FC"
    return "NS"


def read_table(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        first = fh.readline()
    sep = "," if first.count(",") > first.count("\t") else "\t"
    df = pd.read_csv(path, sep=sep)
    df.columns = [str(c).replace("\ufeff", "") for c in df.columns]
    return df


def collapse_symbols(df: pd.DataFrame, symbol_col: str, count_cols: list[str]) -> pd.DataFrame:
    sym = df[symbol_col].astype(str).str.replace("'", "", regex=False).str.strip()
    # drop empty / nan symbols
    keep = sym.notna() & ~sym.str.lower().isin(["", "nan", "none", "-"])
    sub = df.loc[keep, count_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    sub.index = sym[keep]
    # expected counts can be fractional; sum duplicates then round
    summed = sub.groupby(level=0).sum()
    return summed


def deseq_contrast(counts: pd.DataFrame, case: list[str], control: list[str]) -> pd.DataFrame:
    """counts: genes x samples, non-negative. Returns gene-indexed results."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    cols = list(control) + list(case)
    mat = counts[cols].copy()
    mat = mat.loc[mat.sum(axis=1) >= 10]
    # PyDESeq2 wants samples x genes, integers
    rounded = mat.T.round().astype(int)
    rounded.index = rounded.index.astype(str)
    meta = pd.DataFrame(
        {"condition": ["control"] * len(control) + ["case"] * len(case)},
        index=rounded.index,
    )
    dds = DeseqDataSet(
        counts=rounded,
        metadata=meta,
        design="~condition",
        quiet=True,
        refit_cooks=True,
    )
    dds.deseq2()
    statres = DeseqStats(
        dds,
        contrast=["condition", "case", "control"],
        quiet=True,
        cooks_filter=False,
        independent_filter=False,
    )
    statres.summary()
    res = statres.results_df.copy()
    res.index = res.index.astype(str)
    return res


def rows_from_deseq(
    res: pd.DataFrame,
    counts: pd.DataFrame,
    case: list[str],
    control: list[str],
    organism: str,
    meta: dict,
) -> list[dict]:
    lookup = panel_lookup(organism)
    # size-factor-free means of raw (or expected) counts for reporting
    recs = []
    pvals = []
    genes = []
    for gene, cls in lookup.items():
        if gene not in counts.index and gene not in res.index:
            recs.append(_blank(meta, gene, cls, "ABSENT"))
            continue
        mean_case = float(counts.loc[gene, case].mean()) if gene in counts.index else np.nan
        mean_ctrl = float(counts.loc[gene, control].mean()) if gene in counts.index else np.nan
        low = (np.nan_to_num(mean_case) + np.nan_to_num(mean_ctrl)) < 1.0
        if gene not in res.index:
            log2fc = math.log2((mean_case + 0.5) / (mean_ctrl + 0.5)) if np.isfinite(mean_case) else np.nan
            p = np.nan
            qg = np.nan
            statv = np.nan
            base = np.nan
        else:
            log2fc = float(res.loc[gene, "log2FoldChange"])
            p = float(res.loc[gene, "pvalue"]) if np.isfinite(res.loc[gene, "pvalue"]) else np.nan
            qg = float(res.loc[gene, "padj"]) if "padj" in res.columns and np.isfinite(res.loc[gene, "padj"]) else np.nan
            statv = float(res.loc[gene, "stat"]) if "stat" in res.columns and np.isfinite(res.loc[gene, "stat"]) else np.nan
            base = float(res.loc[gene, "baseMean"]) if "baseMean" in res.columns else np.nan
        recs.append(
            {
                **meta,
                "gene": gene,
                "gene_class": cls,
                "mean_case": mean_case,
                "mean_control": mean_ctrl,
                "baseMean": base,
                "log2fc_case_vs_control": log2fc,
                "stat": statv,
                "pvalue": p,
                "padj_genome": qg,
                "low_expression": low,
                "_tested": True,
            }
        )
        pvals.append(p)
    q_panel = bh(pvals)
    qi = 0
    for rec in recs:
        if not rec.get("_tested"):
            continue
        q = q_panel[qi]
        qi += 1
        rec["padj_panel"] = q
        rec["direction"] = direction(rec["log2fc_case_vs_control"], q, rec["low_expression"])
        rec.pop("_tested", None)
    assert qi == len(pvals)
    return recs


def _blank(meta, gene, cls, direction_s) -> dict:
    return {
        **meta,
        "gene": gene,
        "gene_class": cls,
        "mean_case": np.nan,
        "mean_control": np.nan,
        "baseMean": np.nan,
        "log2fc_case_vs_control": np.nan,
        "stat": np.nan,
        "pvalue": np.nan,
        "padj_genome": np.nan,
        "padj_panel": np.nan,
        "low_expression": True,
        "direction": direction_s,
    }


def run_deseq(counts, case, control, organism, meta) -> list[dict]:
    try:
        res = deseq_contrast(counts, case, control)
        meta = {**meta, "method": "pydeseq2_wald"}
        return rows_from_deseq(res, counts, case, control, organism, meta)
    except Exception as exc:  # noqa: BLE001
        print(f"DESeq failed {meta.get('accession')} {meta.get('contrast')}: {exc}")
        return run_logcpm_welch(counts, case, control, organism, {**meta, "method": f"logcpm_welch_fallback:{type(exc).__name__}"})


def run_logcpm_welch(counts, case, control, organism, meta) -> list[dict]:
    lookup = panel_lookup(organism)
    lib = counts.sum(axis=0).replace(0, np.nan)
    logcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)
    recs = []
    pvals = []
    for gene, cls in lookup.items():
        if gene not in counts.index:
            recs.append(_blank(meta, gene, cls, "ABSENT"))
            continue
        a = logcpm.loc[gene, case].astype(float).to_numpy()
        b = logcpm.loc[gene, control].astype(float).to_numpy()
        mean_case = float(counts.loc[gene, case].mean())
        mean_ctrl = float(counts.loc[gene, control].mean())
        low = (mean_case + mean_ctrl) < 1.0
        log2fc = math.log2((mean_case + 0.5) / (mean_ctrl + 0.5))
        if len(a) >= 2 and len(b) >= 2 and np.nanstd(a) + np.nanstd(b) > 0:
            statv, p = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
            statv, p = float(statv), float(p)
        else:
            statv, p = np.nan, np.nan
        recs.append(
            {
                **meta,
                "gene": gene,
                "gene_class": cls,
                "mean_case": mean_case,
                "mean_control": mean_ctrl,
                "baseMean": np.nan,
                "log2fc_case_vs_control": log2fc,
                "stat": statv,
                "pvalue": p,
                "padj_genome": np.nan,
                "low_expression": low,
                "_tested": True,
            }
        )
        pvals.append(p)
    q_panel = bh(pvals)
    qi = 0
    for rec in recs:
        if not rec.get("_tested"):
            continue
        q = q_panel[qi]
        qi += 1
        rec["padj_panel"] = q
        rec["direction"] = direction(rec["log2fc_case_vs_control"], q, rec["low_expression"])
        rec.pop("_tested", None)
    return recs


def run_mwu_logcpm(counts, case, control, organism, meta) -> list[dict]:
    """Cell-level test. Not a biological-replicate DE."""
    lookup = panel_lookup(organism)
    lib = counts.sum(axis=0).replace(0, np.nan)
    logcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)
    recs = []
    pvals = []
    for gene, cls in lookup.items():
        if gene not in counts.index:
            recs.append(_blank(meta, gene, cls, "ABSENT"))
            continue
        a = counts.loc[gene, case].astype(float).to_numpy()
        b = counts.loc[gene, control].astype(float).to_numpy()
        mean_case = float(np.mean(a))
        mean_ctrl = float(np.mean(b))
        low = (mean_case + mean_ctrl) < 0.05  # single-cell counts are sparse
        log2fc = math.log2((mean_case + 0.5) / (mean_ctrl + 0.5))
        xa = logcpm.loc[gene, case].astype(float).to_numpy()
        xb = logcpm.loc[gene, control].astype(float).to_numpy()
        if np.all(a == 0) and np.all(b == 0):
            statv, p, rbc = np.nan, np.nan, np.nan
        else:
            u, p = stats.mannwhitneyu(xa, xb, alternative="two-sided")
            statv = float(u)
            p = float(p)
            rbc = float((2.0 * u) / (len(xa) * len(xb)) - 1.0)
        recs.append(
            {
                **meta,
                "gene": gene,
                "gene_class": cls,
                "mean_case": mean_case,
                "mean_control": mean_ctrl,
                "baseMean": np.nan,
                "log2fc_case_vs_control": log2fc,
                "stat": statv,
                "pvalue": p,
                "padj_genome": np.nan,
                "low_expression": low,
                "rank_biserial": rbc if gene in counts.index else np.nan,
                "detect_case": float(np.mean(a > 0)),
                "detect_control": float(np.mean(b > 0)),
                "_tested": True,
            }
        )
        pvals.append(p)
    q_panel = bh(pvals)
    qi = 0
    for rec in recs:
        if not rec.get("_tested"):
            continue
        q = q_panel[qi]
        qi += 1
        rec["padj_panel"] = q
        rec["direction"] = direction(rec["log2fc_case_vs_control"], q, rec["low_expression"])
        rec.pop("_tested", None)
    return recs


def descriptive_n1(counts, case, control, organism, meta) -> list[dict]:
    lookup = panel_lookup(organism)
    recs = []
    for gene, cls in lookup.items():
        if gene not in counts.index:
            recs.append(_blank(meta, gene, cls, "ABSENT"))
            continue
        mean_case = float(counts.loc[gene, case].mean())
        mean_ctrl = float(counts.loc[gene, control].mean())
        low = (mean_case + mean_ctrl) < 1.0
        log2fc = math.log2((mean_case + 0.5) / (mean_ctrl + 0.5))
        if low:
            direc = "LOW"
        elif log2fc >= 0.5:
            direc = "UP_DESCRIPTIVE"
        elif log2fc <= -0.5:
            direc = "DOWN_DESCRIPTIVE"
        else:
            direc = "FLAT_DESCRIPTIVE"
        recs.append(
            {
                **meta,
                "gene": gene,
                "gene_class": cls,
                "mean_case": mean_case,
                "mean_control": mean_ctrl,
                "baseMean": np.nan,
                "log2fc_case_vs_control": log2fc,
                "stat": np.nan,
                "pvalue": np.nan,
                "padj_genome": np.nan,
                "padj_panel": np.nan,
                "low_expression": low,
                "direction": direc,
            }
        )
    return recs


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_gse279172() -> tuple[pd.DataFrame, dict[str, list[str]]]:
    df = read_table(CACHE / "GSE279172_gene_count.txt.gz")
    # series matrix !Sample_description is "sample 1".."sample 18" in this order
    titles = [
        "WT_high_1", "WT_high_2", "WT_high_3",
        "WT_low_1", "WT_low_2", "WT_low_3",
        "IFNL23KO_high_1", "IFNL23KO_high_2", "IFNL23KO_high_3",
        "IFNL23KO_low_1", "IFNL23KO_low_2", "IFNL23KO_low_3",
        "IFNLRKO_high_1", "IFNLRKO_high_2", "IFNLRKO_high_3",
        "IFNLRKO_low_1", "IFNLRKO_low_2", "IFNLRKO_low_3",
    ]
    cols = [f"sample {i}" for i in range(1, 19)]
    counts = collapse_symbols(df, "gene_name", cols)
    counts.columns = titles
    groups = {
        "WT_high": [c for c in titles if c.startswith("WT_high")],
        "WT_low": [c for c in titles if c.startswith("WT_low")],
        "IFNL23KO_high": [c for c in titles if c.startswith("IFNL23KO_high")],
        "IFNL23KO_low": [c for c in titles if c.startswith("IFNL23KO_low")],
        "IFNLRKO_high": [c for c in titles if c.startswith("IFNLRKO_high")],
        "IFNLRKO_low": [c for c in titles if c.startswith("IFNLRKO_low")],
    }
    return counts, groups


def load_gse296527() -> pd.DataFrame:
    df = read_table(CACHE / "GSE296527_gene_count.txt.gz")
    cols = [c for c in df.columns if c != "gene_id"]
    # gene_name is not in this file; it is a subset. Join symbols from GSE279172
    # for shared Ensembl ids, and keep IFNL1-only columns.
    ref = read_table(CACHE / "GSE279172_gene_count.txt.gz")[["gene_id", "gene_name"]]
    df = df.merge(ref.drop_duplicates("gene_id"), on="gene_id", how="left")
    df["gene_name"] = df["gene_name"].fillna(df["gene_id"])
    return collapse_symbols(df, "gene_name", [c for c in cols if c != "gene_name"])


def load_symbol_matrix(path: Path, symbol_col: str, drop_cols: set[str]) -> pd.DataFrame:
    df = read_table(path)
    cols = [c for c in df.columns if c not in drop_cols and c != symbol_col]
    return collapse_symbols(df, symbol_col, cols)


def load_gse100771() -> tuple[pd.DataFrame, pd.Series]:
    counts = read_table(CACHE / "GSE100771_counts.csv.gz")
    counts = counts.rename(columns={counts.columns[0]: "gene"})
    ph = read_table(CACHE / "GSE100771_Pheno.csv.gz")
    # description holds the MR library id used as the count column
    ph["mr"] = ph["description"].astype(str)
    ph["mn"] = ph["characteristics_ch1.2"].astype(str).str.contains("no micronuclei")
    ph["group"] = np.where(ph["mn"], "MN_neg", "MN_pos")
    # characteristics text "no micronuclei" is the negative set.
    # The positive set says "micronuclei" without "no".
    # Fix: "no micronuclei" contains the substring micronuclei, so check no first.
    gene_cols = [c for c in counts.columns if c != "gene"]
    mat = collapse_symbols(counts, "gene", gene_cols)
    group = ph.set_index("mr")["group"]
    return mat, group


def main() -> None:
    all_rows: list[dict] = []

    print("GSE279172")
    c279, g279 = load_gse279172()
    contrasts_279 = [
        ("WT_high_vs_WT_low", g279["WT_high"], g279["WT_low"],
         "Density contrast. High density is the state where basal IFNL depends on cGAS-STING mtDNA in the paper."),
        ("IFNL23KO_vs_WT_high", g279["IFNL23KO_high"], g279["WT_high"],
         "Ligand KO vs WT at high density, where basal IFNL2/3 is on."),
        ("IFNLRKO_vs_WT_high", g279["IFNLRKO_high"], g279["WT_high"],
         "Receptor KO vs WT at high density."),
        ("IFNL23KO_vs_WT_low", g279["IFNL23KO_low"], g279["WT_low"],
         "Same ligand KO at low density, where basal IFNL is off."),
        ("IFNLRKO_vs_WT_low", g279["IFNLRKO_low"], g279["WT_low"],
         "Same receptor KO at low density."),
    ]
    for name, case, ctrl, note in contrasts_279:
        meta = dict(
            accession="GSE279172", contrast=name, organism="human",
            n_case=len(case), n_control=len(ctrl), unit="bulk_library",
            note=note,
        )
        all_rows.extend(run_deseq(c279, case, ctrl, "human", meta))

    print("GSE296527 IFNL1 only")
    c296 = load_gse296527()
    # WT_high_* in this file match GSE279172 high-density WT (shared libraries).
    meta = dict(
        accession="GSE296527", contrast="IFNL1KO_vs_WT",
        organism="human", n_case=3, n_control=3, unit="bulk_library",
        note="IFNL1 KO is new. WT_high columns are the same three libraries as GSE279172 WT high, not an independent cohort.",
    )
    all_rows.extend(run_deseq(
        c296, ["dL1_high_1", "dL1_high_2", "dL1_high_3"],
        ["WT_high_1", "WT_high_2", "WT_high_3"], "human", meta,
    ))

    print("GSE315942")
    c315 = load_symbol_matrix(
        CACHE / "GSE315942_gene_count.txt.gz",
        "gene_name",
        {"gene_id", "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_length", "gene_biotype", "gene_description", "tf_family"},
    )
    for line, prefix, label in [("OE33", "O", "OE33"), ("SKGT4", "S", "SK-GT-4")]:
        specs = [
            (f"{line}_Cas9_reversine_vs_DMSO",
             [f"{prefix}_Cas9_Rev_{i}" for i in (1, 2, 3)],
             [f"{prefix}_Cas9_Ctr_{i}" for i in (1, 2, 3)],
             f"{label} Cas9: reversine (MPS1 inhibitor, micronuclei/CIN) vs DMSO."),
            (f"{line}_cGASKO_reversine_vs_DMSO",
             [f"{prefix}_KO_Rev_{i}" for i in (1, 2, 3)],
             [f"{prefix}_KO_Ctr_{i}" for i in (1, 2, 3)],
             f"{label} cGAS KO: reversine vs DMSO. Tests whether the reversine shift needs cGAS."),
            (f"{line}_cGASKO_vs_Cas9_DMSO",
             [f"{prefix}_KO_Ctr_{i}" for i in (1, 2, 3)],
             [f"{prefix}_Cas9_Ctr_{i}" for i in (1, 2, 3)],
             f"{label} cGAS KO vs Cas9 in DMSO."),
            (f"{line}_cGASKO_vs_Cas9_reversine",
             [f"{prefix}_KO_Rev_{i}" for i in (1, 2, 3)],
             [f"{prefix}_Cas9_Rev_{i}" for i in (1, 2, 3)],
             f"{label} cGAS KO vs Cas9 under reversine."),
        ]
        for name, case, ctrl, note in specs:
            meta = dict(
                accession="GSE315942", contrast=name, organism="human",
                n_case=3, n_control=3, unit="bulk_library", note=note,
            )
            all_rows.extend(run_deseq(c315, case, ctrl, "human", meta))

    print("GSE316127 clone pseudobulk")
    c316 = load_symbol_matrix(
        CACHE / "GSE316127_gene_count.txt.gz",
        "gene_name",
        {"gene_id", "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_length", "gene_biotype", "gene_description", "tf_family"},
    )
    # File suffixes are _1 _2 _4 (no _3). Sum = one pseudobulk per clone.
    clones_p53 = ["P53c2", "P53c4", "P53c9"]
    suff = ["1", "2", "4"]

    def pseudobulk(prefixes: list[str]) -> pd.DataFrame:
        cols = []
        for p in prefixes:
            parts = [f"{p}_{s}" for s in suff]
            missing = [x for x in parts if x not in c316.columns]
            if missing:
                raise SystemExit(f"missing columns {missing}")
            cols.append(c316[parts].sum(axis=1).rename(p))
        return pd.concat(cols, axis=1)

    pb_dn = pseudobulk([f"{c}_dn" for c in clones_p53])
    pb_p53 = pseudobulk(clones_p53)
    pb = pd.concat([pb_dn, pb_p53], axis=1)
    meta = dict(
        accession="GSE316127",
        contrast="CP-A_p53KO_dnMCAK_vs_p53KO",
        organism="human",
        n_case=3,
        n_control=3,
        unit="clone_pseudobulk",
        note="CP-A Barrett cells. dnMCAK raises chromosome missegregation (CIN) on a p53-KO background. Three clones; replicates _1/_2/_4 summed. Cas9 and parental columns are one line each and are not in this test. Same paper as the cGAS-micronucleus EAC study; this matrix is not a cGAS knockout.",
    )
    all_rows.extend(run_deseq(
        pb,
        [f"{c}_dn" for c in clones_p53],
        clones_p53,
        "human",
        meta,
    ))

    print("GSE237615")
    id8 = read_table(CACHE / "GSE237615_ID8.txt.gz")
    id8ko = read_table(CACHE / "GSE237615_ID8ko.txt.gz")
    kpc = read_table(CACHE / "GSE237615_KPC.txt.gz")
    kpcko = read_table(CACHE / "GSE237615_KPCko.txt.gz")
    id8 = id8.rename(columns={"Gene Symbol": "symbol"})
    id8ko = id8ko.rename(columns={"Gene Symbol": "symbol"})
    kpc = kpc.rename(columns={"Gene Symbol": "symbol"})
    kpcko = kpcko.rename(columns={"Gene Symbol": "symbol"})
    id8_counts = collapse_symbols(id8, "symbol", [c for c in id8.columns if "Read Count" in c])
    id8ko_counts = collapse_symbols(id8ko, "symbol", [c for c in id8ko.columns if "Read Count" in c])
    kpc_counts = collapse_symbols(kpc, "symbol", [c for c in kpc.columns if "Read Count" in c])
    kpcko_counts = collapse_symbols(kpcko, "symbol", [c for c in kpcko.columns if "Read Count" in c])
    id8_counts.columns = [f"ID8_ctrl_{i}" for i in range(1, id8_counts.shape[1] + 1)]
    id8ko_counts.columns = [f"ID8_ko_{i}" for i in range(1, id8ko_counts.shape[1] + 1)]
    kpc_counts.columns = [f"KPC_ctrl_{i}" for i in range(1, kpc_counts.shape[1] + 1)]
    kpcko_counts.columns = [f"KPC_ko_{i}" for i in range(1, kpcko_counts.shape[1] + 1)]
    both_id8 = id8ko_counts.join(id8_counts, how="outer").fillna(0)
    both_kpc = kpcko_counts.join(kpc_counts, how="outer").fillna(0)
    all_rows.extend(run_deseq(
        both_id8, list(id8ko_counts.columns), list(id8_counts.columns), "mouse",
        dict(accession="GSE237615", contrast="ID8_53BP1KO_vs_ctrl", organism="mouse",
             n_case=id8ko_counts.shape[1], n_control=id8_counts.shape[1],
             unit="tumor_library",
             note="ID8 ovarian tumors. 53BP1 KO increases micronuclei and cytosolic dsDNA and is cGAS-STING dependent in the paper. Bulk tumor RNA, not sorted epithelium. Deposited values are fractional expected counts, rounded to integers."),
    ))
    all_rows.extend(run_deseq(
        both_kpc, list(kpcko_counts.columns), list(kpc_counts.columns), "mouse",
        dict(accession="GSE237615", contrast="KPC_53BP1KO_vs_ctrl", organism="mouse",
             n_case=kpcko_counts.shape[1], n_control=kpc_counts.shape[1],
             unit="tumor_library",
             note="KPC pancreatic tumors. Same 53BP1-KO micronucleus design. Bulk tumor RNA. Expected counts rounded to integers."),
    ))

    print("GSE313807")
    c313 = load_symbol_matrix(
        CACHE / "GSE313807_gene_count.txt.gz",
        "gene_name",
        {"gene_id", "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_length", "gene_biotype", "gene_description", "tf_family", "AA_1", "AA_2", "Plas_1", "Plas_2"},
    )
    # SRM order from series matrix library names
    srm = {
        "SRM1": "MSI_MSI_MN", "SRM2": "MSI_CIN_MN", "SRM3": "MSI_MSI_free", "SRM4": "MSI_CIN_free",
        "SRM5": "CIN_MSI_MN", "SRM6": "CIN_CIN_MN", "SRM7": "CIN_MSI_free", "SRM8": "CIN_CIN_free",
        "SRM9": "MSI_MSI_MN", "SRM10": "MSI_CIN_MN", "SRM11": "MSI_MSI_free", "SRM12": "MSI_CIN_free",
        "SRM13": "CIN_MSI_MN", "SRM14": "CIN_CIN_MN", "SRM15": "CIN_MSI_free", "SRM16": "CIN_CIN_free",
    }
    # keep only SRM columns
    c313 = c313[[c for c in c313.columns if c in srm]]
    mn_cols = [c for c, lab in srm.items() if lab.endswith("_MN") and c in c313.columns]
    free_cols = [c for c, lab in srm.items() if lab.endswith("_free") and c in c313.columns]
    msi_cols = [c for c, lab in srm.items() if lab.startswith("MSI_") and c in c313.columns]
    cin_cols = [c for c, lab in srm.items() if lab.startswith("CIN_") and c in c313.columns]
    all_rows.extend(run_deseq(
        c313, mn_cols, free_cols, "mouse",
        dict(accession="GSE313807", contrast="micronuclei_vs_free_cyDNA", organism="mouse",
             n_case=len(mn_cols), n_control=len(free_cols), unit="tumor_library",
             note="Orthotopic MC38 tumors after BMDC transfer of micronuclei vs free cytosolic DNA. Pools MSI and CIN tumors and both DNA sources. n=2 libraries per cell of the full design; this pool is n=8 vs 8. AA_* and Plas_* columns are not in the 16-library sample list and are excluded. Bulk tumor, immune confound."),
    ))
    all_rows.extend(run_deseq(
        c313, msi_cols, cin_cols, "mouse",
        dict(accession="GSE313807", contrast="MSI_tumor_vs_CIN_tumor", organism="mouse",
             n_case=len(msi_cols), n_control=len(cin_cols), unit="tumor_library",
             note="MSI vs CIN MC38 tumors, stimuli pooled. Bulk tumor."),
    ))

    print("GSE99028")
    c990 = read_table(CACHE / "GSE99028_read_count.txt.gz")
    c990 = collapse_symbols(c990, "refGene", [c for c in c990.columns if c != "refGene"])
    all_rows.extend(run_deseq(
        c990, ["shcGas-etoposide-rep1", "shcGas-etoposide-rep2"],
        ["SC-etoposide-rep1", "SC-etoposide-rep2"], "human",
        dict(accession="GSE99028", contrast="shcGAS_vs_shNTC_etoposide", organism="human",
             n_case=2, n_control=2, unit="bulk_library",
             note="IMR90 senescence (Dou 2017 cytosolic chromatin). cGAS knockdown vs control under etoposide. n=2. Underpowered."),
    ))

    print("GSE100102")
    df102 = read_table(CACHE / "GSE100102_Genes_expression_GEO.txt.gz")
    count_cols = [c for c in df102.columns if c.startswith("counts.")]
    mat102 = collapse_symbols(df102, "GeneName", count_cols)
    mat102.columns = [c.replace("counts.", "") for c in mat102.columns]
    for tp in ("tp1", "tp3", "tp4"):
        case = [f"cgas_{tp}.{i}" for i in (1, 2, 3, 4)]
        ctrl = [f"wt_{tp}.{i}" for i in (1, 2, 3, 4)]
        all_rows.extend(run_deseq(
            mat102, case, ctrl, "mouse",
            dict(accession="GSE100102", contrast=f"cGAS_KO_vs_WT_{tp}", organism="mouse",
                 n_case=4, n_control=4, unit="bulk_library",
                 note="Primary MEFs, 20% oxygen passage (Glück 2017 cytosolic chromatin fragments). tp1/tp3/tp4 are passages 11/12/14. Deposited 'counts' are fractional (expected-count scale, sum ~1e7) and were rounded. Not a claudin paper."),
        ))

    print("GSE100771 cell-level")
    mat771, group = load_gse100771()
    # keep columns that are in both
    common = [c for c in mat771.columns if c in set(group.index)]
    mat771 = mat771[common]
    case = [c for c in common if group.loc[c] == "MN_pos"]
    ctrl = [c for c in common if group.loc[c] == "MN_neg"]
    all_rows.extend(run_mwu_logcpm(
        mat771, case, ctrl, "mouse",
        dict(accession="GSE100771", contrast="MN_pos_vs_MN_neg_MEF", organism="mouse",
             n_case=len(case), n_control=len(ctrl), unit="single_cell_one_dish",
             method="mannwhitney_log2cpm",
             note="Mackenzie 2017. Laser-captured MEFs from one irradiated dish, MN+ vs MN-. Cell-level test, not biological-replicate DE. Claudins are not an expected MEF program."),
    ))

    print("GSE98183 descriptive n=1")
    pe = read_table(CACHE / "GSE98183_Counts.geneSymbols.75bpPE.csv.gz")
    # first column may be unnamed gene symbol
    symcol = pe.columns[0]
    pe_counts = collapse_symbols(pe, symcol, [c for c in pe.columns if c != symcol])
    # columns MKH (CIN-high dnMCAK), cont (CIN-medium control)
    all_rows.extend(descriptive_n1(
        pe_counts, ["MKH"], ["cont"], "human",
        dict(accession="GSE98183", contrast="MKH_CINhigh_vs_cont_n1", organism="human",
             n_case=1, n_control=1, unit="one_library",
             method="descriptive_log2_ratio_n1",
             note="Bakhoum 2018 MDA-MB-231 75 bp PE. MKH = dominant-negative MCAK, CIN-high, more micronuclei. cont = control, CIN-medium. One library per genotype. No p-value. 101 bp SE is a second protocol, not a biological replicate, and is not averaged in."),
    ))

    res = pd.DataFrame(all_rows)
    # stable column order
    front = [
        "accession", "contrast", "organism", "unit", "n_case", "n_control",
        "method", "gene", "gene_class", "mean_case", "mean_control", "baseMean",
        "log2fc_case_vs_control", "stat", "pvalue", "padj_genome", "padj_panel",
        "low_expression", "direction", "detect_case", "detect_control",
        "rank_biserial", "note",
    ]
    for col in front:
        if col not in res.columns:
            res[col] = np.nan
    res = res[front]
    res.to_csv(OUT / "de_panel.tsv", sep="\t", index=False)

    key_genes = {
        "CLDN2", "CLDN3", "CLDN4", "CLDN7", "TACSTD2", "TJP1", "OCLN", "CDH1",
        "CGAS", "STING1", "IFNL1", "IFNL2", "IFNL3", "IFNLR1", "ISG15", "MX1",
        "CXCL10", "Cldn2", "Cldn3", "Cldn4", "Cldn7", "Tacstd2", "Tjp1", "Ocln",
        "Cdh1", "Cgas", "Mb21d1", "Sting1", "Isg15", "Mx1", "Cxcl10",
    }
    key = res[res["gene"].isin(key_genes)].copy()
    key.to_csv(OUT / "de_key_genes.tsv", sep="\t", index=False)

    # compact call table: non-LOW claudin/TJ/surface only
    focus = res[res["gene_class"].isin(["claudin", "tight_junction", "surface"])].copy()
    focus.to_csv(OUT / "de_claudin_tj.tsv", sep="\t", index=False)
    print("wrote", OUT, "rows", len(res))
    # print key directions that are not NS/LOW/ABSENT
    show = key[~key["direction"].isin(["NS", "LOW", "ABSENT", "FLAT_DESCRIPTIVE"])]
    cols = ["accession", "contrast", "gene", "log2fc_case_vs_control", "pvalue", "padj_panel", "direction", "n_case", "n_control"]
    print(show[cols].to_string(index=False))


if __name__ == "__main__":
    main()
