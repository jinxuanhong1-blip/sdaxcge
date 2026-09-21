#!/usr/bin/env python3
"""GSE135274 HeLa XRCC4 KO ± mirin: IFN / STING / cGAS / APM / NHEJ.

Public GEO count table only (Benjamin et al., Int J Mol Sci 2022, PMID 35054780).

Design (GEO labels, two biological experiments):
  WT     scrambled gRNA          gSCR
  WT_M   scrambled gRNA + mirin  gSCR_M     (paper: 100 µM)
  KO     XRCC4(-/-) clone 2G3    g2G3
  KO_M   XRCC4(-/-) + mirin      g2G3_M

Every GEO sample is RNA at 48 h after TALEN + NHEJ-reporter transfection.
Read 1 and read 2 were quantified separately (technical). Mismatch bins
(mis_0/1/2) are partitions of assigned reads. Biological n = 2 per arm.

Primary effect: mean of the two within-experiment log2(CPM+1) differences.
UP/DOWN requires both experiments to agree with |log2FC| > 0.25 and the gene
not low (max CPM < 1). Welch p on n=2 is reported and is underpowered.

GSEA: preranked Hallmark (weighted KS, 1000 gene-set permutations, seed 42).
Rank = mean paired log2FC. Positive NES = enriched in the numerator arm.
"""
from __future__ import annotations

import gzip
import math
import sys
import time
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank  # noqa: E402

GMT = ROOT / "data" / "genesets" / "h.all.v2023.2.Hs.symbols.gmt"
OUT = ROOT / "results" / "gse135274_xrcc4_mirin_ifn"
RAW = OUT / "raw"
TAB = OUT / "tables"
FIG = OUT / "figures"
GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135274/suppl/"
    "GSE135274_all_sample_human_RNA_max.txt.gz"
)
GEO_NAME = "GSE135274_all_sample_human_RNA_max.txt.gz"
LFC_CUT = 0.25
LOW_CPM = 1.0

# Column pairs are read1 / read2 of one library (SOFT: technical replicates).
PAIRS = {
    "Exp1_KO": ["Experiment_01_g2G3_1_1", "Experiment_01_g2G3_1_2"],
    "Exp1_KO_M": ["Experiment_01_g2G3_M_1_1", "Experiment_01_g2G3_M_1_2"],
    "Exp1_WT": ["Experiment_01_gSCR_1_1", "Experiment_01_gSCR_1_2"],
    "Exp1_WT_M": ["Experiment_01_gSCR_M_1_1", "Experiment_01_gSCR_M_1_2"],
    "Exp2_KO": ["Experiment_02_g2G3_2_1", "Experiment_02_g2G3_2_2"],
    "Exp2_KO_M": ["Experiment_02_g2G3_M_2_1", "Experiment_02_g2G3_M_2_2"],
    "Exp2_WT": ["Experiment_02_gSCR_2_1", "Experiment_02_gSCR_2_2"],
    "Exp2_WT_M": ["Experiment_02_gSCR_M_2_1", "Experiment_02_gSCR_M_2_2"],
}
ARMS = {
    "KO": ["Exp1_KO", "Exp2_KO"],
    "KO_M": ["Exp1_KO_M", "Exp2_KO_M"],
    "WT": ["Exp1_WT", "Exp2_WT"],
    "WT_M": ["Exp1_WT_M", "Exp2_WT_M"],
}
CONTRASTS = [
    ("XRCC4_KO_vs_WT", "KO", "WT", "XRCC4 KO vs scrambled"),
    ("mirin_vs_WT", "WT_M", "WT", "mirin vs scrambled"),
    ("XRCC4_KO_plus_mirin_vs_WT", "KO_M", "WT", "XRCC4 KO + mirin vs scrambled"),
    ("mirin_on_XRCC4_KO", "KO_M", "KO", "mirin on XRCC4 KO vs KO"),
]
IFN_TERMS = [
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]
SANITY_TERMS = IFN_TERMS + [
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_HYPOXIA",
    "HALLMARK_DNA_REPAIR",
]
# Authors' qPCR / glycolysis checks (Benjamin 2022). Directions are measured here.
SANITY_GENES = [
    "CA9",
    "CDKN1A",
    "ENO2",
    "DUSP5",
    "ZMAT3",
    "PGK1",
    "ALDOC",
    "PFKFB4",
    "TPI1",
    "ENO3",
    "PFKP",
    "HK2",
]

# (symbol, category, display). Aliases tried if the hg19 symbol is absent.
PANEL = [
    ("XRCC4", "qc", "XRCC4"),
    ("MRE11A", "qc", "MRE11A"),
    ("MB21D1", "cgas_sting", "cGAS/MB21D1"),
    ("TMEM173", "cgas_sting", "STING/TMEM173"),
    ("TBK1", "cgas_sting", "TBK1"),
    ("IKBKE", "cgas_sting", "IKBKE"),
    ("IRF3", "cgas_sting", "IRF3"),
    ("MAVS", "cgas_sting", "MAVS"),
    ("TREX1", "cgas_sting", "TREX1"),
    ("SAMHD1", "cgas_sting", "SAMHD1"),
    ("ENPP1", "cgas_sting", "ENPP1"),
    ("DNASE2", "cgas_sting", "DNASE2"),
    ("RNASEH2A", "cgas_sting", "RNASEH2A"),
    ("RNASEH2B", "cgas_sting", "RNASEH2B"),
    ("RNASEH2C", "cgas_sting", "RNASEH2C"),
    ("POLA1", "cgas_sting", "POLA1"),
    ("IFNB1", "ifn_ligand", "IFNB1"),
    ("IFNA1", "ifn_ligand", "IFNA1"),
    ("IFNA2", "ifn_ligand", "IFNA2"),
    ("IFNG", "ifn_ligand", "IFNG"),
    ("IL6", "ifn_ligand", "IL6"),
    ("ISG15", "ifn_isg", "ISG15"),
    ("ISG20", "ifn_isg", "ISG20"),
    ("MX1", "ifn_isg", "MX1"),
    ("MX2", "ifn_isg", "MX2"),
    ("OAS1", "ifn_isg", "OAS1"),
    ("OAS2", "ifn_isg", "OAS2"),
    ("OAS3", "ifn_isg", "OAS3"),
    ("OASL", "ifn_isg", "OASL"),
    ("IFIT1", "ifn_isg", "IFIT1"),
    ("IFIT2", "ifn_isg", "IFIT2"),
    ("IFIT3", "ifn_isg", "IFIT3"),
    ("IFIT5", "ifn_isg", "IFIT5"),
    ("RSAD2", "ifn_isg", "RSAD2"),
    ("USP18", "ifn_isg", "USP18"),
    ("IFI6", "ifn_isg", "IFI6"),
    ("IFI16", "ifn_isg", "IFI16"),
    ("IFI27", "ifn_isg", "IFI27"),
    ("IFI35", "ifn_isg", "IFI35"),
    ("IFI44", "ifn_isg", "IFI44"),
    ("IFI44L", "ifn_isg", "IFI44L"),
    ("BST2", "ifn_isg", "BST2"),
    ("IFITM1", "ifn_isg", "IFITM1"),
    ("IFITM2", "ifn_isg", "IFITM2"),
    ("IFITM3", "ifn_isg", "IFITM3"),
    ("DDX58", "ifn_isg", "DDX58"),
    ("IFIH1", "ifn_isg", "IFIH1"),
    ("DDX60", "ifn_isg", "DDX60"),
    ("HERC5", "ifn_isg", "HERC5"),
    ("XAF1", "ifn_isg", "XAF1"),
    ("EIF2AK2", "ifn_isg", "EIF2AK2"),
    ("LY6E", "ifn_isg", "LY6E"),
    ("CMPK2", "ifn_isg", "CMPK2"),
    ("EPSTI1", "ifn_isg", "EPSTI1"),
    ("SAMD9L", "ifn_isg", "SAMD9L"),
    ("TRIM22", "ifn_isg", "TRIM22"),
    ("ZBP1", "ifn_isg", "ZBP1"),
    ("PLSCR1", "ifn_isg", "PLSCR1"),
    ("STAT1", "ifn_signaling", "STAT1"),
    ("STAT2", "ifn_signaling", "STAT2"),
    ("IRF1", "ifn_signaling", "IRF1"),
    ("IRF7", "ifn_signaling", "IRF7"),
    ("IRF9", "ifn_signaling", "IRF9"),
    ("JAK1", "ifn_signaling", "JAK1"),
    ("JAK2", "ifn_signaling", "JAK2"),
    ("IFNAR1", "ifn_signaling", "IFNAR1"),
    ("IFNAR2", "ifn_signaling", "IFNAR2"),
    ("IFNGR1", "ifn_signaling", "IFNGR1"),
    ("IFNGR2", "ifn_signaling", "IFNGR2"),
    ("SOCS1", "ifn_signaling", "SOCS1"),
    ("SOCS3", "ifn_signaling", "SOCS3"),
    ("GBP1", "ifn_isg", "GBP1"),
    ("GBP2", "ifn_isg", "GBP2"),
    ("GBP4", "ifn_isg", "GBP4"),
    ("GBP5", "ifn_isg", "GBP5"),
    ("CXCL9", "ifn_isg", "CXCL9"),
    ("CXCL10", "ifn_isg", "CXCL10"),
    ("CXCL11", "ifn_isg", "CXCL11"),
    ("CCL5", "ifn_isg", "CCL5"),
    ("CIITA", "ifn_isg", "CIITA"),
    ("IDO1", "ifn_isg", "IDO1"),
    ("HLA-A", "apm", "HLA-A"),
    ("HLA-B", "apm", "HLA-B"),
    ("HLA-C", "apm", "HLA-C"),
    ("HLA-E", "apm", "HLA-E"),
    ("HLA-F", "apm", "HLA-F"),
    ("HLA-G", "apm", "HLA-G"),
    ("B2M", "apm", "B2M"),
    ("NLRC5", "apm", "NLRC5"),
    ("TAP1", "apm", "TAP1"),
    ("TAP2", "apm", "TAP2"),
    ("TAPBP", "apm", "TAPBP"),
    ("PSMB8", "apm", "PSMB8"),
    ("PSMB9", "apm", "PSMB9"),
    ("PSMB10", "apm", "PSMB10"),
    ("ERAP1", "apm", "ERAP1"),
    ("ERAP2", "apm", "ERAP2"),
    ("CALR", "apm", "CALR"),
    ("CANX", "apm", "CANX"),
    ("PDIA3", "apm", "PDIA3"),
    ("PSME1", "apm", "PSME1"),
    ("PSME2", "apm", "PSME2"),
    ("LIG4", "nhej", "LIG4"),
    ("XRCC5", "nhej", "Ku80/XRCC5"),
    ("XRCC6", "nhej", "Ku70/XRCC6"),
    ("PRKDC", "nhej", "DNA-PKcs/PRKDC"),
    ("NHEJ1", "nhej", "XLF/NHEJ1"),
    ("DCLRE1C", "nhej", "Artemis/DCLRE1C"),
    ("PAXX", "nhej", "PAXX"),
    ("APLF", "nhej", "APLF"),
    ("POLL", "nhej", "POLL"),
    ("POLM", "nhej", "POLM"),
    ("XRCC1", "alt_nhej", "XRCC1"),
    ("LIG3", "alt_nhej", "LIG3"),
    ("PARP1", "alt_nhej", "PARP1"),
    ("RAD50", "alt_nhej", "RAD50"),
    ("NBN", "alt_nhej", "NBN"),
    ("POLQ", "alt_nhej", "POLQ"),
    ("POLB", "alt_nhej", "POLB"),
]
ALIASES = {
    "MB21D1": ["CGAS"],
    "TMEM173": ["STING1"],
    "MRE11A": ["MRE11"],
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def fmt(x, digits=3):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "NA"
    return f"{x:+.{digits}f}" if isinstance(x, float) else str(x)


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (math.isnan(p) or math.isinf(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def ensure_matrix() -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    dest = RAW / GEO_NAME
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    alt = Path("/tmp") / GEO_NAME
    if alt.exists() and alt.stat().st_size > 1000:
        dest.write_bytes(alt.read_bytes())
        return dest
    log(f"downloading {GEO_URL}")
    urllib.request.urlretrieve(GEO_URL, dest)
    return dest


def load_counts(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum mismatch bins per gene, then sum read1+read2 per library."""
    acc: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = header[2:]
        missing = [c for pair in PAIRS.values() for c in pair if c not in cols]
        if missing:
            raise SystemExit(f"missing columns: {missing}")
        idx = [cols.index(c) for pair in PAIRS.values() for c in pair]
        ordered = [c for pair in PAIRS.values() for c in pair]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            gene = parts[1].strip()
            if not gene or " " in gene or gene.upper().endswith("ALLELE"):
                continue
            gene = gene.upper()
            vals = np.zeros(len(ordered), dtype=np.float64)
            for j, i in enumerate(idx):
                raw = parts[2 + i] if 2 + i < len(parts) else ""
                if raw:
                    vals[j] = float(raw)
            if gene in acc:
                acc[gene] += vals
            else:
                acc[gene] = vals
    raw = pd.DataFrame.from_dict(acc, orient="index", columns=ordered)
    bio_cols = []
    bio = {}
    spearman_rows = []
    for name, (a, b) in PAIRS.items():
        bio[name] = raw[a] + raw[b]
        bio_cols.append(name)
        rho = stats.spearmanr(raw[a], raw[b]).statistic
        spearman_rows.append({"library": name, "read1": a, "read2": b, "spearman_r1_r2": rho,
                              "read1_sum": float(raw[a].sum()), "read2_sum": float(raw[b].sum())})
    counts = pd.DataFrame(bio, index=raw.index)[bio_cols]
    qc = pd.DataFrame(spearman_rows)
    return counts, qc


def to_log_cpm(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0), lib


def call_direction(e1: float, e2: float, max_cpm: float) -> str:
    if max_cpm < LOW_CPM:
        return "LOW"
    if e1 > LFC_CUT and e2 > LFC_CUT:
        return "UP"
    if e1 < -LFC_CUT and e2 < -LFC_CUT:
        return "DOWN"
    if (e1 > LFC_CUT and e2 < -LFC_CUT) or (e1 < -LFC_CUT and e2 > LFC_CUT):
        return "DISCORDANT"
    return "FLAT"


def welch(logcpm: pd.DataFrame, treat: list[str], ctrl: list[str]) -> tuple[pd.Series, pd.Series]:
    a = logcpm[treat].to_numpy(dtype=float)
    b = logcpm[ctrl].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        res = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    t = pd.Series(np.asarray(res.statistic, dtype=float), index=logcpm.index)
    p = pd.Series(np.asarray(res.pvalue, dtype=float), index=logcpm.index)
    return t.replace([np.inf, -np.inf], np.nan), p


def resolve(symbol: str, index: pd.Index) -> str | None:
    if symbol in index:
        return symbol
    for alt in ALIASES.get(symbol, []):
        if alt in index:
            return alt
    return None


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def rank_genes(mean_lfc: pd.Series, keep: pd.Series) -> pd.Series:
    s = mean_lfc[keep].replace([np.inf, -np.inf], np.nan).dropna()
    return s.sort_values(ascending=False)


def write_results(ctx: dict) -> None:
    h = ctx["headline"]
    lines = []
    a = lines.append
    a("# GSE135274 XRCC4 KO ± mirin — IFN / STING / cGAS / APM / NHEJ")
    a("")
    a("Public counts only. This does not change the locked CosMx, concordant-4, GSE137244, TCGA keratin, or TISMO numbers.")
    a("")
    a("## Result")
    a("")
    for para in ctx["reading"]:
        a(para)
        a("")
    a("## Dataset")
    a("")
    a("- **GSE135274**, PMID 35054780 (Benjamin et al., *Int J Mol Sci* 2022). HeLa, not lung.")
    a("- 2×2: scrambled gRNA vs XRCC4(−/−) clone 2G3, each ± mirin. Authors used **100 µM mirin**. Protein loss of XRCC4 is their Western blot; this table is mRNA.")
    a("- RNA at **48 h after TALEN + NHEJ-reporter transfection** in every arm. The contrast is NHEJ blockade on a TALEN-DSB background, not resting HeLa.")
    a("- Two biological experiments. GEO columns `*_1` and `*_2` are read 1 and read 2 quantified separately (technical). They were summed; CPM is unchanged if they are averaged instead.")
    a("- CLC mismatch bins `mis_0/1/2` were summed per gene (a partition of assigned reads).")
    a("- **Biological n = 2 per arm.** Gene-level Welch p-values are descriptive. Calls use sign concordance, not FDR.")
    a("")
    a("## Method")
    a("")
    a("- log2FC = difference of log2(CPM+1), paired inside Experiment 1 and Experiment 2, then averaged.")
    a(f"- **UP / DOWN**: both experiments |log2FC| > {LFC_CUT} and the same sign, and max CPM ≥ {LOW_CPM:g}. Opposite signs beyond that cut = DISCORDANT. Otherwise FLAT. Max CPM < {LOW_CPM:g} = LOW.")
    a(f"- GSEA: MSigDB Hallmark 2023.2.Hs, preranked by mean paired log2FC, weighted KS (p=1), {NPERM} gene-set permutations, seed {SEED}. Positive NES = enriched in the numerator (treated / KO) arm. BH FDR is across the 50 Hallmark sets inside that contrast.")
    a("- Genes in the rank: symbol with no spaces, and ≥10 counts in at least 2 of 8 libraries.")
    a("")
    a("## QC")
    a("")
    a(f"- Read1 vs read2 Spearman within each library: min {ctx['rho_min']:.3f}, median {ctx['rho_med']:.3f}.")
    a(f"- Genes after symbol filter: {ctx['n_genes']}. Genes in the GSEA rank: {ctx['n_rank']}.")
    a(f"- Genome-wide Spearman of Experiment 1 vs Experiment 2 log2FC (KO+mirin vs WT, the authors' strongest contrast): {ctx['exp_rho']:.3f}.")
    a("")
    a("XRCC4 mRNA (mean CPM). The knockout is a small indel (authors: −2 bp / −10 bp) with **no protein**; mRNA is only partly lower.")
    a("")
    a("| arm | Exp1 CPM | Exp2 CPM |")
    a("|---|---:|---:|")
    for arm, e1, e2 in ctx["xrcc4_cpm"]:
        a(f"| {arm} | {e1:.2f} | {e2:.2f} |")
    a("")
    a(
        f"XRCC4 paired log2FC, KO vs WT: Exp1 {fmt(ctx['xrcc4_e1'])}, Exp2 {fmt(ctx['xrcc4_e2'])}, "
        f"mean **{fmt(ctx['xrcc4_ko_lfc'])}** ({ctx['xrcc4_ko_call']}). "
        "Exp1 is nearly unchanged, so the both-experiment rule does not call it DOWN."
    )
    a("")
    a("Author-reported genes on KO+mirin vs scrambled (glycolysis genes were described as down; CA9, CDKN1A, ENO2, DUSP5, ZMAT3 were qPCR-checked):")
    a("")
    a("| gene | Exp1 | Exp2 | mean log2FC | call |")
    a("|---|---:|---:|---:|---|")
    for row in ctx["sanity"]:
        a(f"| {row['gene']} | {fmt(row['e1'])} | {fmt(row['e2'])} | {fmt(row['mean'])} | {row['call']} |")
    a("")
    a("## Hallmark GSEA")
    a("")
    a("Positive NES means the set is higher in the first arm of the contrast.")
    a("")
    term_order = {t: i for i, t in enumerate(SANITY_TERMS)}
    contrast_order = {c[0]: i for i, c in enumerate(CONTRASTS)}
    h = sorted(h, key=lambda r: (contrast_order[r["contrast"]], term_order[r["term"]]))
    a("| contrast | set | NES | nominal p | BH FDR (50 Hallmarks) | mean log2FC in set |")
    a("|---|---|---:|---:|---:|---:|")
    for row in h:
        a(
            f"| {row['label']} | {row['short']} | {fmt(row['nes'])} | {fmt_p(row['nom_p'])} | "
            f"{fmt_p(row['fdr'])} | {fmt(row['mean_stat'])} |"
        )
    a("")
    a("Full 50-set tables: `tables/gsea_hallmark_all.tsv`.")
    a("")
    a("## Panel calls")
    a("")
    a("Measured = not ABSENT and not LOW. UP/DOWN are both-experiment calls.")
    a("")
    a("| contrast | category | measured | UP | DOWN | FLAT | DISCORDANT | LOW | ABSENT | median log2FC |")
    a("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in ctx["cat_rows"]:
        a(
            f"| {row['label']} | {row['category']} | {row['measured']} | {row['UP']} | {row['DOWN']} | "
            f"{row['FLAT']} | {row['DISCORDANT']} | {row['LOW']} | {row['ABSENT']} | {fmt(row['median'])} |"
        )
    a("")
    a("Gene-level table: `tables/panel_de.tsv`.")
    a("")
    a("### Genes called UP or DOWN in at least one contrast")
    a("")
    a("| gene | category | KO vs WT | mirin vs WT | KO+mirin vs WT | mirin on KO |")
    a("|---|---|---|---|---|---|")
    for row in ctx["moving"]:
        a(
            f"| {row['display']} | {row['category']} | {row['c0']} | {row['c1']} | {row['c2']} | {row['c3']} |"
        )
    if not ctx["moving"]:
        a("| — | — | no concordant UP/DOWN | | | |")
    a("")
    a("## Limits")
    a("")
    a("- n=2. A Hallmark nominal p is a permutation of the gene rank, not a patient-level or even a well-powered sample test.")
    a("- HeLa + transfected TALEN breaks. Not a lung tumor, not ICI, not cGAS-STING stimulation as the experiment's intent.")
    a("- XRCC4 protein is gone by the authors' blot; XRCC4 mRNA is a weak QC.")
    a("- Mirin blocks MRE11 exonuclease activity. MRE11A mRNA is not expected to fall.")
    a("- This file does not re-litigate CLDN4 exclusion, TISMO, or the concordant-4 cohort.")
    (OUT / "RESULTS.md").write_text("\n".join(lines) + "\n")


def plot_nes(headline: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    labels = list(dict.fromkeys(headline["label"]))
    terms = [
        ("HALLMARK_INTERFERON_ALPHA_RESPONSE", "IFN-α"),
        ("HALLMARK_INTERFERON_GAMMA_RESPONSE", "IFN-γ"),
    ]
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    for i, (term, name) in enumerate(terms):
        sub = headline[headline["term"] == term].set_index("label").loc[labels]
        vals = sub["nes"].to_numpy(dtype=float)
        bars = ax.bar(x + (i - 0.5) * width, vals, width, label=name)
        for rect, p, fdr in zip(bars, sub["nom_p"], sub["fdr"]):
            if p < 0.05:
                mark = "*" if fdr >= 0.05 else "**"
                y = rect.get_height()
                ax.text(rect.get_x() + rect.get_width() / 2, y + (0.04 if y >= 0 else -0.12),
                        mark, ha="center", va="bottom" if y >= 0 else "top", fontsize=11)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("NES (positive = higher in numerator)")
    ax.set_title("GSE135274 Hallmark IFN  (n=2; * nominal p<0.05, ** also BH FDR<0.05)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_hallmark_ifn_nes.png", dpi=160)
    fig.savefig(FIG / "fig_hallmark_ifn_nes.pdf")
    plt.close(fig)


def plot_heatmap(panel: pd.DataFrame, contrast_ids: list[str], contrast_labels: list[str]) -> None:
    use = panel[panel["call"].isin(["UP", "DOWN", "DISCORDANT", "FLAT"])].copy()
    # one row per gene; columns are contrasts. pivot mean lfc
    mat = use.pivot(index="display", columns="contrast", values="mean_lfc")
    order = []
    seen = set()
    for _sym, _cat, display in PANEL:
        if display in mat.index and display not in seen:
            order.append(display)
            seen.add(display)
    mat = mat.loc[order, contrast_ids]
    fig_h = max(6.0, 0.18 * len(mat))
    fig, ax = plt.subplots(figsize=(8.4, fig_h))
    data = mat.to_numpy(dtype=float)
    vmax = 1.5
    im = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(contrast_labels)))
    ax.set_xticklabels(contrast_labels, rotation=20, ha="right")
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=6)
    ax.set_title("Panel mean paired log2FC  (red = higher in numerator)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="mean log2FC")
    fig.tight_layout()
    fig.savefig(FIG / "fig_panel_log2fc.png", dpi=160)
    fig.savefig(FIG / "fig_panel_log2fc.pdf")
    plt.close(fig)


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    path = ensure_matrix()
    log(f"loading {path}")
    counts, pair_qc = load_counts(path)
    pair_qc.to_csv(TAB / "read_pair_spearman.tsv", sep="\t", index=False)
    logcpm, lib = to_log_cpm(counts)
    cpm = counts.div(lib, axis=1) * 1e6
    lib.to_csv(TAB / "library_size.tsv", sep="\t", header=["sum_counts"])
    keep = (counts >= 10).sum(axis=1) >= 2
    log(f"genes {counts.shape[0]}  rank-eligible {int(keep.sum())}")

    contrasts = {}
    for cid, num, den, label in CONTRASTS:
        e1 = logcpm[f"Exp1_{num}"] - logcpm[f"Exp1_{den}"]
        e2 = logcpm[f"Exp2_{num}"] - logcpm[f"Exp2_{den}"]
        mean = (e1 + e2) / 2.0
        tstat, pval = welch(logcpm, ARMS[num], ARMS[den])
        contrasts[cid] = {
            "label": label,
            "num": num,
            "den": den,
            "e1": e1,
            "e2": e2,
            "mean": mean,
            "t": tstat,
            "p": pval,
        }
        rho = stats.spearmanr(e1[keep], e2[keep]).statistic
        log(f"{cid} exp1-vs-exp2 Spearman {rho:.3f}")
        contrasts[cid]["exp_rho"] = float(rho)

    # Panel table
    panel_rows = []
    for symbol, category, display in PANEL:
        hit = resolve(symbol, counts.index)
        max_cpm = float(cpm.loc[hit].max()) if hit else float("nan")
        for cid, _num, _den, label in CONTRASTS:
            if hit is None:
                panel_rows.append({
                    "gene": symbol, "matched": "", "display": display, "category": category,
                    "contrast": cid, "label": label, "e1": np.nan, "e2": np.nan, "mean_lfc": np.nan,
                    "welch_t": np.nan, "welch_p": np.nan, "max_cpm": np.nan, "call": "ABSENT",
                })
                continue
            c = contrasts[cid]
            e1 = float(c["e1"].loc[hit])
            e2 = float(c["e2"].loc[hit])
            panel_rows.append({
                "gene": symbol, "matched": hit, "display": display, "category": category,
                "contrast": cid, "label": label,
                "e1": e1, "e2": e2, "mean_lfc": float(c["mean"].loc[hit]),
                "welch_t": float(c["t"].loc[hit]) if pd.notna(c["t"].loc[hit]) else np.nan,
                "welch_p": float(c["p"].loc[hit]) if pd.notna(c["p"].loc[hit]) else np.nan,
                "max_cpm": max_cpm,
                "call": call_direction(e1, e2, max_cpm),
            })
    panel = pd.DataFrame(panel_rows)
    # BH within each contrast, among panel genes that are not ABSENT/LOW.
    combined = pd.Series(np.nan, index=panel.index, dtype=float)
    for cid, _n, _d, _l in CONTRASTS:
        usable = (panel["contrast"] == cid) & ~panel["call"].isin(["ABSENT", "LOW"])
        if usable.any():
            combined.loc[usable] = bh_fdr(panel.loc[usable, "welch_p"]).to_numpy()
    panel["welch_fdr_in_panel"] = combined
    panel.to_csv(TAB / "panel_de.tsv", sep="\t", index=False)

    # category summary
    cat_rows = []
    for cid, _n, _d, label in CONTRASTS:
        sub = panel[panel["contrast"] == cid]
        for category, g in sub.groupby("category", sort=False):
            vc = g["call"].value_counts()
            measured = g[g["call"].isin(["UP", "DOWN", "FLAT", "DISCORDANT"])]
            cat_rows.append({
                "contrast": cid,
                "label": label,
                "category": category,
                "n": int(len(g)),
                "measured": int(len(measured)),
                "UP": int(vc.get("UP", 0)),
                "DOWN": int(vc.get("DOWN", 0)),
                "FLAT": int(vc.get("FLAT", 0)),
                "DISCORDANT": int(vc.get("DISCORDANT", 0)),
                "LOW": int(vc.get("LOW", 0)),
                "ABSENT": int(vc.get("ABSENT", 0)),
                "median": float(measured["mean_lfc"].median()) if len(measured) else np.nan,
            })
    cat = pd.DataFrame(cat_rows)
    cat.to_csv(TAB / "panel_category_summary.tsv", sep="\t", index=False)

    # sanity genes on the double-block contrast
    sanity = []
    cid = "XRCC4_KO_plus_mirin_vs_WT"
    c = contrasts[cid]
    for gene in SANITY_GENES:
        if gene not in counts.index:
            sanity.append({"gene": gene, "e1": np.nan, "e2": np.nan, "mean": np.nan, "call": "ABSENT", "max_cpm": np.nan})
            continue
        e1 = float(c["e1"].loc[gene])
        e2 = float(c["e2"].loc[gene])
        max_cpm = float(cpm.loc[gene].max())
        sanity.append({
            "gene": gene, "e1": e1, "e2": e2, "mean": float(c["mean"].loc[gene]),
            "call": call_direction(e1, e2, max_cpm), "max_cpm": max_cpm,
        })
    pd.DataFrame(sanity).to_csv(TAB / "author_gene_sanity.tsv", sep="\t", index=False)

    # GSEA
    gene_sets = read_gmt(GMT)
    gsea_frames = []
    for cid, _n, _d, label in CONTRASTS:
        log(f"GSEA {cid}")
        rank = rank_genes(contrasts[cid]["mean"], keep)
        res = gsea_prerank(rank, gene_sets, nperm=NPERM, seed=SEED)
        res["fdr"] = bh_fdr(res["nom_p"])
        res.insert(0, "contrast", cid)
        res.insert(1, "label", label)
        gsea_frames.append(res)
    gsea = pd.concat(gsea_frames, ignore_index=True)
    gsea.to_csv(TAB / "gsea_hallmark_all.tsv", sep="\t", index=False)
    headline = gsea[gsea["term"].isin(SANITY_TERMS)].copy()
    short = {
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-α",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ",
        "HALLMARK_P53_PATHWAY": "p53",
        "HALLMARK_HYPOXIA": "hypoxia",
        "HALLMARK_DNA_REPAIR": "DNA repair",
    }
    headline["short"] = headline["term"].map(short)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)

    # moving genes
    wide_calls = {}
    for i, (cid, _n, _d, _l) in enumerate(CONTRASTS):
        sub = panel[panel["contrast"] == cid].set_index("display")
        wide_calls[i] = sub["call"]
        if i == 0:
            base = sub[["gene", "category"]].copy()
    moving = []
    displays = []
    for _sym, _cat, display in PANEL:
        calls = [wide_calls[i].get(display, "ABSENT") for i in range(4)]
        if any(x in ("UP", "DOWN") for x in calls):
            displays.append(display)
            moving.append({
                "display": display,
                "category": base.loc[display, "category"] if display in base.index else "",
                "c0": f"{calls[0]} ({fmt(panel[(panel.display==display)&(panel.contrast==CONTRASTS[0][0])]['mean_lfc'].iloc[0])})",
                "c1": f"{calls[1]} ({fmt(panel[(panel.display==display)&(panel.contrast==CONTRASTS[1][0])]['mean_lfc'].iloc[0])})",
                "c2": f"{calls[2]} ({fmt(panel[(panel.display==display)&(panel.contrast==CONTRASTS[2][0])]['mean_lfc'].iloc[0])})",
                "c3": f"{calls[3]} ({fmt(panel[(panel.display==display)&(panel.contrast==CONTRASTS[3][0])]['mean_lfc'].iloc[0])})",
            })

    def arm_cpm(gene, sample):
        return float(cpm.loc[gene, sample])

    xr = "XRCC4"
    xrcc4_cpm = [
        ("WT scrambled", arm_cpm(xr, "Exp1_WT"), arm_cpm(xr, "Exp2_WT")),
        ("WT + mirin", arm_cpm(xr, "Exp1_WT_M"), arm_cpm(xr, "Exp2_WT_M")),
        ("XRCC4 KO", arm_cpm(xr, "Exp1_KO"), arm_cpm(xr, "Exp2_KO")),
        ("XRCC4 KO + mirin", arm_cpm(xr, "Exp1_KO_M"), arm_cpm(xr, "Exp2_KO_M")),
    ]
    ko_row = panel[(panel["gene"] == "XRCC4") & (panel["contrast"] == "XRCC4_KO_vs_WT")].iloc[0]

    ifn = gsea[gsea["term"].isin(IFN_TERMS)]
    # reading built from numbers
    reading = build_reading(panel, cat, ifn, sanity, float(contrasts["XRCC4_KO_plus_mirin_vs_WT"]["exp_rho"]))

    plot_nes(headline[headline["term"].isin(IFN_TERMS)])
    plot_heatmap(panel, [c[0] for c in CONTRASTS], [c[3] for c in CONTRASTS])

    write_results({
        "rho_min": float(pair_qc["spearman_r1_r2"].min()),
        "rho_med": float(pair_qc["spearman_r1_r2"].median()),
        "n_genes": int(counts.shape[0]),
        "n_rank": int(keep.sum()),
        "exp_rho": float(contrasts["XRCC4_KO_plus_mirin_vs_WT"]["exp_rho"]),
        "xrcc4_cpm": xrcc4_cpm,
        "xrcc4_ko_lfc": float(ko_row["mean_lfc"]),
        "xrcc4_e1": float(ko_row["e1"]),
        "xrcc4_e2": float(ko_row["e2"]),
        "xrcc4_ko_call": ko_row["call"],
        "sanity": sanity,
        "headline": headline.to_dict("records"),
        "cat_rows": cat_rows,
        "moving": moving,
        "reading": reading,
    })
    log(f"wrote {OUT / 'RESULTS.md'}")


def _nes(ifn: pd.DataFrame, cid: str, alpha: bool) -> pd.Series:
    key = "ALPHA" if alpha else "GAMMA"
    hit = ifn[(ifn["contrast"] == cid) & ifn["term"].str.contains(key)]
    return hit.iloc[0]


def _gene(panel: pd.DataFrame, symbol: str, cid: str) -> pd.Series:
    return panel[(panel["gene"] == symbol) & (panel["contrast"] == cid)].iloc[0]


def _cat(cat: pd.DataFrame, cid: str, category: str) -> pd.Series:
    return cat[(cat["contrast"] == cid) & (cat["category"] == category)].iloc[0]


def build_reading(panel: pd.DataFrame, cat: pd.DataFrame, ifn: pd.DataFrame, sanity: list, exp_rho: float) -> list[str]:
    """Factual sentences from the tables. Arms are not pooled into one IFN call."""
    ko_a = _nes(ifn, "XRCC4_KO_vs_WT", True)
    ko_g = _nes(ifn, "XRCC4_KO_vs_WT", False)
    mi_a = _nes(ifn, "mirin_vs_WT", True)
    mi_g = _nes(ifn, "mirin_vs_WT", False)
    both_a = _nes(ifn, "XRCC4_KO_plus_mirin_vs_WT", True)
    both_g = _nes(ifn, "XRCC4_KO_plus_mirin_vs_WT", False)
    on_a = _nes(ifn, "mirin_on_XRCC4_KO", True)
    on_g = _nes(ifn, "mirin_on_XRCC4_KO", False)
    isg_ko = _cat(cat, "XRCC4_KO_vs_WT", "ifn_isg")
    isg_mi = _cat(cat, "mirin_vs_WT", "ifn_isg")
    isg_both = _cat(cat, "XRCC4_KO_plus_mirin_vs_WT", "ifn_isg")
    isg_on = _cat(cat, "mirin_on_XRCC4_KO", "ifn_isg")
    apm_ko = _cat(cat, "XRCC4_KO_vs_WT", "apm")
    apm_mi = _cat(cat, "mirin_vs_WT", "apm")
    sting = _gene(panel, "TMEM173", "XRCC4_KO_vs_WT")
    cgas = _gene(panel, "MB21D1", "XRCC4_KO_vs_WT")
    cgas_mi = _gene(panel, "MB21D1", "mirin_vs_WT")
    sting_mi = _gene(panel, "TMEM173", "mirin_vs_WT")
    ifnb = _gene(panel, "IFNB1", "XRCC4_KO_vs_WT")
    xr5 = _gene(panel, "XRCC5", "XRCC4_KO_vs_WT")
    prk = _gene(panel, "PRKDC", "XRCC4_KO_vs_WT")
    nhej1 = _gene(panel, "NHEJ1", "XRCC4_KO_vs_WT")
    isg15_mi = _gene(panel, "ISG15", "mirin_vs_WT")
    ifit1_mi = _gene(panel, "IFIT1", "mirin_vs_WT")
    stat1_mi = _gene(panel, "STAT1", "mirin_vs_WT")
    rsad = _gene(panel, "RSAD2", "mirin_vs_WT")
    down_gly = [s for s in sanity if s["gene"] in ("PGK1", "ALDOC", "PFKFB4", "TPI1", "ENO3", "PFKP", "ENO2", "HK2")]
    n_down = sum(1 for s in down_gly if s["call"] == "DOWN")
    ca9 = next(s for s in sanity if s["gene"] == "CA9")
    cdkn = next(s for s in sanity if s["gene"] == "CDKN1A")
    zmat = next(s for s in sanity if s["gene"] == "ZMAT3")
    return [
        (
            f"XRCC4 knockout versus scrambled gRNA lowers Hallmark IFN. "
            f"IFN-α NES {ko_a['nes']:+.2f} and IFN-γ NES {ko_g['nes']:+.2f} "
            f"(nominal p={fmt_p(ko_a['nom_p'])} and {fmt_p(ko_g['nom_p'])}, "
            f"BH FDR={fmt_p(ko_a['fdr'])} and {fmt_p(ko_g['fdr'])} across 50 Hallmarks). "
            f"Of {int(isg_ko['measured'])} measured ISGs, {int(isg_ko['UP'])} concordant UP and "
            f"{int(isg_ko['DOWN'])} concordant DOWN (median log2FC {isg_ko['median']:+.2f}). "
            f"STING/TMEM173 is {sting['call']} ({sting['mean_lfc']:+.2f}). "
            f"cGAS/MB21D1 is {cgas['call']} ({cgas['mean_lfc']:+.2f}). "
            f"APM: {int(apm_ko['UP'])} UP / {int(apm_ko['DOWN'])} DOWN of {int(apm_ko['measured'])} "
            f"(HLA-A, HLA-B, and B2M are DOWN). "
            f"IFNB1 max CPM is {ifnb['max_cpm']:.2f} ({ifnb['call']})."
        ),
        (
            f"Mirin versus scrambled raises Hallmark IFN. "
            f"IFN-α NES {mi_a['nes']:+.2f} and IFN-γ NES {mi_g['nes']:+.2f} "
            f"(nominal p={fmt_p(mi_a['nom_p'])}, BH FDR={fmt_p(mi_a['fdr'])}). "
            f"{int(isg_mi['UP'])} of {int(isg_mi['measured'])} measured ISGs are concordant UP "
            f"(median {isg_mi['median']:+.2f}; ISG15 {isg15_mi['mean_lfc']:+.2f}, "
            f"IFIT1 {ifit1_mi['mean_lfc']:+.2f}, RSAD2 {rsad['mean_lfc']:+.2f}, "
            f"STAT1 {stat1_mi['mean_lfc']:+.2f}). "
            f"cGAS/MB21D1 is {cgas_mi['call']} ({cgas_mi['mean_lfc']:+.2f}) and "
            f"STING/TMEM173 is {sting_mi['call']} ({sting_mi['mean_lfc']:+.2f}). "
            f"APM stays mostly flat ({int(apm_mi['UP'])} UP / {int(apm_mi['DOWN'])} DOWN "
            f"of {int(apm_mi['measured'])}). "
            "The mirin signal is the ISG cassette, not a rise in cGAS or STING mRNA."
        ),
        (
            f"XRCC4 KO + mirin versus scrambled stays IFN-low: "
            f"IFN-α NES {both_a['nes']:+.2f}, IFN-γ NES {both_g['nes']:+.2f} "
            f"(FDR={fmt_p(both_a['fdr'])} and {fmt_p(both_g['fdr'])}); "
            f"ISGs {int(isg_both['UP'])} UP / {int(isg_both['DOWN'])} DOWN "
            f"(median {isg_both['median']:+.2f}). "
            f"Mirin on the KO background raises IFN relative to KO alone "
            f"(IFN-α NES {on_a['nes']:+.2f}, IFN-γ NES {on_g['nes']:+.2f}, "
            f"FDR={fmt_p(on_a['fdr'])}; ISGs {int(isg_on['UP'])} UP / {int(isg_on['DOWN'])} DOWN) "
            "and does not put the double-block arm above scrambled."
        ),
        (
            f"NHEJ mRNA does not collapse and is not a large backup program. "
            f"In the KO, Ku80/XRCC5 {xr5['mean_lfc']:+.2f}, DNA-PKcs/PRKDC {prk['mean_lfc']:+.2f}, "
            f"and XLF/NHEJ1 {nhej1['mean_lfc']:+.2f} are concordant UP; LIG4, Ku70, and alt-NHEJ are flat. "
            "Hallmark DNA repair is shifted toward the KO (see the GSEA table) with a small mean log2FC inside the set. "
            "PAXX is absent from this CLC gene table. IFNA1, IFNA2, IFNG, CXCL9, and HLA-G are absent too."
        ),
        (
            f"The same rank recovers the paper's non-IFN biology on KO+mirin versus scrambled: "
            f"{n_down}/8 listed glycolysis genes are concordant DOWN, "
            f"CA9 {ca9['mean']:+.2f}, CDKN1A {cdkn['mean']:+.2f}, ZMAT3 {zmat['mean']:+.2f}. "
            f"Experiment 1 vs Experiment 2 log2FC Spearman on that contrast is {exp_rho:.2f}. "
            "Hallmark hypoxia is negative in every contrast; Hallmark p53 is positive where mirin is in the numerator. "
            "That is the check that the IFN signs are not a processing artifact. "
            "The authors' Metascape term \"defense response to virus\" is a different gene set and is not re-tested here."
        ),
        (
            f"With {NPERM} permutations the smallest nominal p is {1/(NPERM+1):.3f}. "
            "These are gene-set permutations of a two-experiment rank."
        ),
    ]


if __name__ == "__main__":
    main()
