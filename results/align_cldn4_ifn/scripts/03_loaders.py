#!/usr/bin/env python3
"""Dataset loaders.

Each loader returns a Contrast with either
  * a per-sample log2 expression matrix (`logmat`) + group labels, or
  * a precomputed per-gene log2FC series (`lfc`) when GEO only deposited
    collapsed / already-differential data.

Row index is always a gene SYMBOL for the relevant species.
"""
import gzip
import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

DATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
sys_path_marker = None


@dataclass
class Contrast:
    key: str
    gse: str
    arm: str                 # "CLDN4" or "TACSTD2"
    species: str             # "human" / "mouse"
    label: str
    perturb_gene: str        # gene expected to go DOWN
    logmat: pd.DataFrame = None
    grp_test: list = field(default_factory=list)
    grp_ref: list = field(default_factory=list)
    lfc: pd.Series = None    # supplied directly when no per-sample data
    p: pd.Series = None
    caveats: str = ""
    n_test: int = 0
    n_ref: int = 0


def _collapse(df, how="max_mean"):
    """Collapse duplicate symbols: keep the row with the highest mean."""
    if not df.index.has_duplicates:
        return df
    order = df.mean(axis=1).sort_values(ascending=False).index
    return df.loc[order][~df.loc[order].index.duplicated()].sort_index()


# --------------------------------------------------------------- GSE207704
def gse207704():
    """Human T47D / MCF7 CLDN4-/- vs WT. GEO deposit is cuffdiff FPKM already
    collapsed to one column per group, so replicate-level stats are impossible
    from the processed file."""
    df = pd.read_csv(os.path.join(DATA, "GSE207704_CLDN4_RNAseq.txt.gz"),
                     sep="\t", low_memory=False)
    df = df.rename(columns={"gene_short_name": "symbol"})
    cols = {c: c.split(" ")[0] for c in df.columns if c.endswith("(fpkm)")}
    df = df.rename(columns=cols)
    v = df.groupby("symbol")[list(cols.values())].max()
    out = []
    for line in ("T47D", "MCF7"):
        ko, wt = f"{line}_CLDN4KO_FPKM", f"{line}_WT_FPKM"
        sub = v[[ko, wt]]
        keep = sub.max(axis=1) >= 1.0          # expressed in at least one arm
        sub = sub[keep]
        lfc = np.log2((sub[ko] + 0.1) / (sub[wt] + 0.1))
        out.append(Contrast(
            key=f"GSE207704_{line}", gse="GSE207704", arm="CLDN4",
            species="human", label=f"{line} CLDN4-/- vs WT (breast)",
            perturb_gene="CLDN4", lfc=lfc, n_test=2, n_ref=2,
            caveats=("GEO deposit is group-collapsed cuffdiff FPKM (one value "
                     "per genotype); the 2 biological replicates per arm are "
                     "not separable, so no replicate-level test is possible.")))
    return out


# ---------------------------------------------------------------- GSE22493
_ALIAS = {"G1P2": "ISG15", "G1P3": "IFI6"}


def _gse22493_symbol(orf, desc):
    s = (orf or "").strip().upper()
    if not s:
        s = (desc or "").split("--")[0].strip().upper()
    return _ALIAS.get(s, s)


def gse22493():
    """SKOV-3 two-colour array. Primary numbers are the GEO-deposited VALUE
    (author-processed log2 knockdown/control; series matrix). ScanArray
    intensities are a sensitivity check only and are written separately."""
    plat_path = os.path.join(DATA, "GPL10555_platform.tsv.gz")
    mat_path = os.path.join(DATA, "GSE22493_series_matrix.txt.gz")
    plat = pd.read_csv(plat_path, sep="\t", dtype=str)
    rows, hdr = [], None
    import gzip as _gz
    with _gz.open(mat_path, "rt") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                hdr = next(fh).rstrip("\n").replace('"', "").split("\t")
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if hdr is not None:
                rows.append(line.rstrip("\n").replace('"', "").split("\t"))
    mat = pd.DataFrame(rows, columns=hdr)
    for c in mat.columns[1:]:
        mat[c] = pd.to_numeric(mat[c], errors="coerce")
    plat["ID"] = plat["ID"].astype(str)
    mat["ID_REF"] = mat["ID_REF"].astype(str)
    m = mat.merge(plat, left_on="ID_REF", right_on="ID", how="left")
    m["symbol"] = [_gse22493_symbol(o, d) for o, d in
                   zip(m["ORF"].fillna(""), m["DESCRIPTION"].fillna(""))]
    m = m[m["symbol"].str.len() > 0]
    gsms = ["GSM558700", "GSM558701", "GSM558702"]
    gene = m.groupby("symbol")[gsms].median()
    return [Contrast(
        key="GSE22493_SKOV3", gse="GSE22493", arm="CLDN4", species="human",
        label="SKOV-3 CLDN4-KD vs CLDN4-OE control (ovarian, deposited VALUE)",
        perturb_gene="CLDN4", logmat=gene,
        grp_test=[], grp_ref=[], lfc=gene.mean(axis=1), n_test=3, n_ref=3,
        caveats=("Two-colour Operon v3 array, 2006. VALUE is author-processed "
                 "log2(KD/control). Control is CLDN4-OVEREXPRESSING cells, not "
                 "parental/scramble. CLDN4 probe is missing on GSM558700; "
                 "GSM558701 deposited vs raw ScanArray disagree in sign. "
                 "Ovarian, not lung."))]


# ---------------------------------------------------------------- GSE50927
def gse50927():
    """Mouse whole lung, Cldn4 KO vs WT (edgeR tables deposited by authors).
    Sign convention verified against Cldn4 itself (logFC = -6.06)."""
    out = []
    spec = [("Cldn4lungWTvsKOgenes", "Cldn4-KO vs WT lung, no VILI"),
            ("VILIwtkohiGenes", "Cldn4-KO(high-injury) vs WT, after VILI"),
            ("VILIwtkoloGenes", "Cldn4-KO(low-injury) vs WT, after VILI")]
    for fn, lab in spec:
        d = pd.read_csv(os.path.join(DATA, f"GSE50927_{fn}.csv.gz"))
        d = d.dropna(subset=["Marker.Symbol"])
        lfc = d.groupby("Marker.Symbol")["logFC"].mean()
        p = d.groupby("Marker.Symbol")["PValue"].min()
        out.append(Contrast(
            key=f"GSE50927_{fn.replace('Genes','').replace('genes','')}",
            gse="GSE50927", arm="CLDN4", species="mouse", label=lab,
            perturb_gene="Cldn4", lfc=lfc, p=p, n_test=1, n_ref=1,
            caveats=("Whole mouse lung, n=1 per group in GEO (5 GSMs total); "
                     "author edgeR p-values come from a fixed-dispersion "
                     "no-replicate fit and are not trustworthy. Bulk tissue, "
                     "so any IFN change also reflects immune-cell content.")))
    return out


# --------------------------------------------------------------- GSE274940
def gse274940():
    """Mouse EpH4 mammary epithelium, complete Cldn-null vs WT (raw counts)."""
    d = pd.read_csv(os.path.join(DATA, "GSE274940_raw_counts.csv.gz"))
    d["symbol"] = d["ALIAS"].fillna(d["ENSEMBL"])
    samples = ["WT1", "WT2", "WT3", "KO1", "KO2", "KO3"]
    cnt = d.groupby("symbol")[samples].sum()
    cnt = cnt[cnt.sum(axis=1) > 0]
    return [Contrast(
        key="GSE274940_EpH4", gse="GSE274940", arm="CLDN4", species="mouse",
        label="EpH4 Cldn-null vs WT (mouse mammary epithelium)",
        perturb_gene="Cldn4", logmat=None, grp_test=["KO1", "KO2", "KO3"],
        grp_ref=["WT1", "WT2", "WT3"], n_test=3, n_ref=3,
        lfc=None, p=None,
        caveats=("NOT a CLDN4-specific perturbation: this line is null for the "
                 "entire claudin family (all 27 members), so it cannot isolate "
                 "CLDN4. Non-transformed mouse mammary epithelium."))], cnt


# --------------------------------------------------------------- GSE334497
def gse334497(sym_map):
    """Mouse 4T1 Trop2-KO vs WT tumours, immunocompetent BALB/c, n=5 v 5."""
    d = pd.read_csv(os.path.join(DATA, "GSE334497_normalized_counts.csv.gz"),
                    index_col=0)
    ko = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
    wt = ["RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R", "control170"]
    assert set(ko + wt) == set(d.columns), sorted(d.columns)
    d.index = [sym_map.get(i, i) for i in d.index]
    d = d.groupby(level=0).sum()
    return Contrast(
        key="GSE334497_4T1", gse="GSE334497", arm="TACSTD2", species="mouse",
        label="4T1 Trop2-KO vs WT tumours in vivo (BALB/c, immunocompetent)",
        perturb_gene="Tacstd2", grp_test=ko, grp_ref=wt, n_test=5, n_ref=5,
        caveats=("Bulk tumour from immunocompetent mice: an IFN increase can "
                 "come from infiltrating immune cells rather than from the "
                 "tumour cell itself. This is the source dataset of the "
                 "TROP2/claudin immune-exclusion paper, so it is not an "
                 "independent test of that paper's claim.")), d


# --------------------------------------------------------------- GSE289287
def gse289287():
    """Human T-47D Trop-2 KO vs WT xenografts (author DESeq2 table)."""
    d = pd.read_csv(os.path.join(DATA, "GSE289287_Trop2KO_tumors_vs_WT.tsv.gz"),
                    sep="\t", low_memory=False)
    d = d.dropna(subset=["Feature_name"])
    wt = ["OV.2808.RNA_normCounts", "OV.2810.RNA_normCounts", "OV.2812.RNA_normCounts"]
    ko = ["OV.2807.RNA_normCounts", "OV.2815.RNA_normCounts",
          "OV.2817.RNA_normCounts", "OV.2818.RNA_normCounts"]
    agg = {c: "sum" for c in wt + ko}
    agg.update({"log2FoldChange": "mean", "pvalue": "min", "baseMean": "sum"})
    g = d.groupby("Feature_name").agg(agg)
    g = g[g["baseMean"] > 0]
    c = Contrast(
        key="GSE289287_T47Dxeno", gse="GSE289287", arm="TACSTD2",
        species="human", label="T-47D Trop-2-KO vs WT xenografts (NRG mice)",
        perturb_gene="TACSTD2", grp_test=ko, grp_ref=wt, n_test=4, n_ref=3,
        lfc=g["log2FoldChange"], p=g["pvalue"],
        caveats=("Xenografts in NRG mice (no T/B/NK cells), so this reads the "
                 "tumour-cell-intrinsic IFN program only. GEO holds no Trop-2 "
                 "KO in-vitro arm (the in-vitro samples are WT and DSG2-KO)."))
    counts = g[wt + ko]
    return c, counts


# --------------------------------------------------------------- GSE245459
def gse245459():
    """Human SKOV3 shTACSTD2 vs shNC, +/- cisplatin (author FPKM)."""
    d = pd.read_csv(os.path.join(DATA, "GSE245459_fpkm.anno.txt.gz"), sep="\t",
                    low_memory=False)
    d = d.dropna(subset=["GeneName"])
    sam = ["shNC1", "shNC2", "shNC3", "shNCDDP1", "shNCDDP2", "shNCDDP3",
           "sh1", "sh2", "sh3", "shDDP1", "shDDP2", "shDDP3"]
    for c in sam:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    fp = d.groupby("GeneName")[sam].max()
    out = []
    for key, test, ref, lab in [
        ("GSE245459_SKOV3", ["sh1", "sh2", "sh3"], ["shNC1", "shNC2", "shNC3"],
         "SKOV3 shTACSTD2 vs shNC (ovarian, untreated)"),
        ("GSE245459_SKOV3_DDP", ["shDDP1", "shDDP2", "shDDP3"],
         ["shNCDDP1", "shNCDDP2", "shNCDDP3"],
         "SKOV3 shTACSTD2 vs shNC, both + cisplatin")]:
        sub = fp[test + ref]
        sub = sub[sub.max(axis=1) >= 0.1]
        out.append((Contrast(
            key=key, gse="GSE245459", arm="TACSTD2", species="human",
            label=lab, perturb_gene="TACSTD2", grp_test=test, grp_ref=ref,
            n_test=3, n_ref=3,
            caveats=("Author-supplied FPKM only (no raw counts); n=3 v 3.")),
            sub))
    return out


def mouse_symbol_map():
    d = pd.read_csv(os.path.join(DATA, "GSE274940_raw_counts.csv.gz"),
                    usecols=["ENSEMBL", "ALIAS"])
    d = d.dropna()
    return dict(zip(d["ENSEMBL"], d["ALIAS"]))
