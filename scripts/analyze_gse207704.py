#!/usr/bin/env python3
"""IFN / MHC-I/APM / tight-junction reanalysis of GSE207704 CLDN4 knockout.

Primary question: in the public breast-cancer CLDN4 CRISPR knockout
(T47D and MCF7; Kage et al., Breast Cancer Res 2023, PMID 37059993),
do interferon-alpha, interferon-gamma, MHC-I antigen-presentation, or
tight-junction genes move relative to wild type?

The deposited GEO table collapses the two replicates and omits several
canonical HLA/APM genes, so counts are rebuilt from the eight SRA
libraries. Calls below were fixed before looking at those counts.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

import gseapy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import hypergeom, mannwhitneyu, spearmanr

SAMPLES = [
    ("SRR20029125", "GSM6310640", "T47D", "WT", "T47D_WT_1"),
    ("SRR20029124", "GSM6310641", "T47D", "WT", "T47D_WT_2"),
    ("SRR20029123", "GSM6310642", "T47D", "KO", "T47D_KO_1"),
    ("SRR20029122", "GSM6310643", "T47D", "KO", "T47D_KO_2"),
    ("SRR20029121", "GSM6310644", "MCF7", "WT", "MCF7_WT_1"),
    ("SRR20029120", "GSM6310645", "MCF7", "WT", "MCF7_WT_2"),
    ("SRR20029119", "GSM6310646", "MCF7", "KO", "MCF7_KO_1"),
    ("SRR20029118", "GSM6310647", "MCF7", "KO", "MCF7_KO_2"),
]

# Classical MHC-I, beta-2-microglobulin, peptide transport, tapasin,
# immunoproteasome, ER aminopeptidases, loading chaperones, and the
# MHC-I transactivator. CLDN4 is not in this list.
APM_GENES = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M",
    "TAP1", "TAP2", "TAPBP", "TAPBPL",
    "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2",
    "ERAP1", "ERAP2", "CALR", "CANX", "PDIA3", "NLRC5",
]

# Canonical tight-junction components. CLDN4 is the knockout control
# and is scored separately so it cannot create a TJ "down" call by itself.
TJ_GENES = [
    "CLDN1", "CLDN2", "CLDN3", "CLDN5", "CLDN6", "CLDN7", "CLDN8",
    "CLDN9", "CLDN10", "CLDN11", "CLDN12", "CLDN14", "CLDN15",
    "CLDN16", "CLDN17", "CLDN18", "CLDN19", "CLDN20", "CLDN23",
    "OCLN", "TJP1", "TJP2", "TJP3", "CGN", "CGNL1",
    "MARVELD2", "MARVELD3", "F11R", "JAM2", "JAM3", "MPDZ",
]

EPITHELIAL_CONTEXT = ["EPCAM", "TACSTD2", "CDH1", "KRT8", "KRT18", "KRT19", "CLDN4"]

IFN_ALPHA = "Interferon Alpha Response"
IFN_GAMMA = "Interferon Gamma Response"
HALLMARK_FOCUS = [
    IFN_ALPHA,
    IFN_GAMMA,
    "Inflammatory Response",
    "IL-6/JAK/STAT3 Signaling",
    "TNF-alpha Signaling via NF-kB",
    "Allograft Rejection",
    "Complement",
    "Apical Junction",
    "Cholesterol Homeostasis",
    "Fatty Acid Metabolism",
    "Bile Acid Metabolism",
    "Estrogen Response Early",
    "Estrogen Response Late",
]

DESIGN_COMBINED = (
    "~ C(cell_line, Treatment(reference='MCF7')) + "
    "C(condition, Treatment(reference='WT'))"
)
DESIGN_LINE = "~ C(condition, Treatment(reference='WT'))"
KO_COEFF = "C(condition, Treatment(reference='WT'))[T.KO]"

# Prespecified effect and significance gates.
MEAN_LFC_GATE = 0.25
FDR_GATE = 0.05
DISCORD_MEDIAN_GATE = 0.25
GSEA_PERM = 1000
GSEA_SEED = 207704


def classify(mean_lfc, nes, fdr, median_t47d, median_mcf7) -> str:
    """Direction call. Fixed before the GSE207704 counts were inspected."""
    discordant = (
        median_t47d * median_mcf7 < 0
        and abs(median_t47d) >= DISCORD_MEDIAN_GATE
        and abs(median_mcf7) >= DISCORD_MEDIAN_GATE
    )
    if discordant:
        return "discordant"
    same_up = median_t47d > 0 and median_mcf7 > 0
    same_down = median_t47d < 0 and median_mcf7 < 0
    strong = abs(mean_lfc) >= MEAN_LFC_GATE
    sig = fdr < FDR_GATE
    if strong and sig and nes > 0 and mean_lfc > 0 and same_up:
        return "up"
    if strong and sig and nes < 0 and mean_lfc < 0 and same_down:
        return "down"
    if strong and sig and nes > 0 and mean_lfc > 0 and not same_up:
        return "line-dependent"
    if strong and sig and nes < 0 and mean_lfc < 0 and not same_down:
        return "line-dependent"
    if strong and sig and nes * mean_lfc < 0:
        return "unresolved"
    if sig and not strong:
        return "weak"
    if strong and not sig:
        return "suggestive"
    if (not strong) and (not sig):
        return "null"
    return "unresolved"


def _self_check() -> None:
    assert classify(0.05, 0.2, 0.4, 0.02, 0.01) == "null"
    assert classify(0.8, 1.5, 0.01, 0.4, 0.3) == "up"
    assert classify(-0.8, -1.5, 0.01, -0.4, -0.3) == "down"
    assert classify(0.8, 1.5, 0.01, 0.4, -0.4) == "discordant"
    assert classify(-0.8, -1.5, 0.01, 0.4, -0.4) == "discordant"
    assert classify(0.05, 2.0, 0.01, 0.02, 0.01) == "weak"
    assert classify(0.6, 1.2, 0.2, 0.3, 0.2) == "suggestive"
    assert classify(0.6, 1.2, 0.01, 0.4, -0.05) == "line-dependent"
    assert classify(0.6, -1.2, 0.01, 0.4, 0.3) == "unresolved"


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            genes = []
            seen = set()
            for gene in parts[2:]:
                if gene and gene not in seen:
                    seen.add(gene)
                    genes.append(gene)
            sets[parts[0]] = genes
    return sets


def parse_tx2gene(fasta: Path) -> pd.DataFrame:
    rows = []
    opener = gzip.open if str(fasta).endswith(".gz") else open
    with opener(fasta, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            tokens = line[1:].split()
            tid = tokens[0]
            symbol = ensg = None
            for token in tokens:
                if token.startswith("gene_symbol:"):
                    symbol = token.split(":", 1)[1]
                elif token.startswith("gene:"):
                    ensg = token.split(":", 1)[1].split(".")[0]
            if symbol and ensg and symbol not in {"-", ""}:
                rows.append((tid, ensg, symbol))
    out = pd.DataFrame(rows, columns=["transcript", "ensg", "symbol"])
    if out.empty:
        raise RuntimeError(f"no gene_symbol headers in {fasta}")
    return out


def load_kallisto(quant_dir: Path, tx2gene: pd.DataFrame):
    tx_map = tx2gene.set_index("transcript")
    count_frames = []
    tpm_frames = []
    qc_rows = []
    for srr, gsm, line, cond, name in SAMPLES:
        abund = quant_dir / srr / "abundance.tsv"
        info = quant_dir / srr / "run_info.json"
        if not abund.exists():
            raise FileNotFoundError(abund)
        df = pd.read_csv(abund, sep="\t", usecols=["target_id", "est_counts", "tpm"])
        df = df.join(tx_map, on="target_id")
        missing = df["ensg"].isna().mean()
        gene = df.dropna(subset=["ensg"]).groupby("ensg", as_index=True)[["est_counts", "tpm"]].sum()
        count_frames.append(gene["est_counts"].rename(name))
        tpm_frames.append(gene["tpm"].rename(name))
        meta = json.loads(info.read_text()) if info.exists() else {}
        qc_rows.append(
            {
                "sample": name,
                "srr": srr,
                "gsm": gsm,
                "cell_line": line,
                "condition": cond,
                "n_processed": meta.get("n_processed"),
                "n_pseudoaligned": meta.get("n_pseudoaligned"),
                "p_pseudoaligned": meta.get("p_pseudoaligned"),
                "p_unique": meta.get("p_unique"),
                "fraction_tx_without_symbol": round(float(missing), 4),
            }
        )
    counts = pd.concat(count_frames, axis=1).fillna(0.0)
    tpm = pd.concat(tpm_frames, axis=1).fillna(0.0)
    symbol = tx2gene.drop_duplicates("ensg").set_index("ensg")["symbol"]
    counts = counts.join(symbol)
    tpm = tpm.join(symbol)
    # If two Ensembl genes share a symbol, keep the gene with more reads.
    collisions = []
    keep_idx = []
    for sym, block in counts.groupby("symbol"):
        if len(block) == 1:
            keep_idx.append(block.index[0])
            continue
        totals = block.drop(columns=["symbol"]).sum(axis=1)
        winner = totals.idxmax()
        keep_idx.append(winner)
        for ensg in block.index:
            if ensg != winner:
                collisions.append({"symbol": sym, "dropped_ensg": ensg, "kept_ensg": winner})
    counts = counts.loc[keep_idx].set_index("symbol")
    tpm = tpm.loc[keep_idx].set_index("symbol")
    counts = counts.groupby(level=0).sum()
    tpm = tpm.groupby(level=0).sum()
    # Round estimated counts to integers for DESeq2.
    rounded = np.floor(counts.to_numpy() + 0.5).astype(int)
    counts = pd.DataFrame(rounded, index=counts.index, columns=counts.columns)
    qc = pd.DataFrame(qc_rows)
    return counts, tpm, qc, pd.DataFrame(collisions)


def metadata() -> pd.DataFrame:
    rows = []
    for _srr, _gsm, line, cond, name in SAMPLES:
        rows.append({"sample": name, "cell_line": line, "condition": cond})
    return pd.DataFrame(rows).set_index("sample")


def run_deseq(counts: pd.DataFrame, meta: pd.DataFrame, design: str) -> pd.DataFrame:
    # Samples x genes, as pydeseq2 expects.
    use = counts.loc[:, meta.index]
    keep = use.sum(axis=1) >= 10
    mat = use.loc[keep].T
    dds = DeseqDataSet(
        counts=mat,
        metadata=meta,
        design=design,
        refit_cooks=False,
        n_cpus=4,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=["condition", "KO", "WT"],
        cooks_filter=False,
        quiet=True,
        n_cpus=4,
    )
    stats.summary()
    unshrunk = stats.results_df["log2FoldChange"].copy()
    wald = stats.results_df["stat"].copy()
    pvalue = stats.results_df["pvalue"].copy()
    padj = stats.results_df["padj"].copy()
    basemean = stats.results_df["baseMean"].copy()
    lfcse = stats.results_df["lfcSE"].copy()
    if KO_COEFF not in stats.LFC.columns:
        raise RuntimeError(f"missing KO coefficient in {list(stats.LFC.columns)}")
    stats.lfc_shrink(KO_COEFF)
    shrunk = stats.LFC @ stats.contrast_vector / np.log(2)
    out = pd.DataFrame(
        {
            "baseMean": basemean,
            "log2FoldChange": unshrunk,
            "log2FoldChange_shrunk": shrunk,
            "lfcSE": lfcse,
            "stat": wald,
            "pvalue": pvalue,
            "padj": padj,
        }
    )
    out.index.name = "gene"
    return out.sort_values("stat", ascending=False)


def gsea_table(ranked: pd.Series, gene_sets: dict[str, list[str]], min_size: int) -> pd.DataFrame:
    rnk = ranked.dropna().astype(float)
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    frame = pd.DataFrame({"gene": rnk.index, "score": rnk.values})
    pre = gseapy.prerank(
        rnk=frame,
        gene_sets=gene_sets,
        permutation_num=GSEA_PERM,
        min_size=min_size,
        max_size=500,
        outdir=None,
        no_plot=True,
        seed=GSEA_SEED,
        threads=4,
        method="multilevel",
        verbose=False,
    )
    res = pre.res2d.copy()
    for col in ["ES", "NES", "NOM p-val", "FDR q-val", "Tag %", "Gene %"]:
        if col in res.columns and col not in {"Tag %", "Gene %"}:
            res[col] = pd.to_numeric(res[col], errors="coerce")
    res = res.rename(
        columns={
            "Term": "geneset",
            "ES": "ES",
            "NES": "NES",
            "NOM p-val": "pval",
            "FDR q-val": "fdr",
            "Tag %": "tag_pct",
            "Gene %": "gene_pct",
            "Lead_genes": "lead_genes",
        }
    )
    res["n_genes_ranked"] = len(rnk)
    res = res.sort_values("NES", ascending=False).reset_index(drop=True)
    res.insert(0, "nes_rank", np.arange(1, len(res) + 1))
    return res


def overlap_table(ranked: pd.Series, gene_sets: dict[str, list[str]], tail_frac: float = 0.10) -> pd.DataFrame:
    rnk = ranked.dropna().astype(float).sort_values(ascending=False)
    universe = list(rnk.index)
    n_univ = len(universe)
    n_tail = max(1, int(round(tail_frac * n_univ)))
    tail = set(universe[:n_tail])
    rows = []
    for name, genes in gene_sets.items():
        inset = [g for g in genes if g in set(universe)]
        k_set = len(inset)
        if k_set == 0:
            continue
        k_hit = sum(g in tail for g in inset)
        # P(X >= k_hit)
        p_up = float(hypergeom.sf(k_hit - 1, n_univ, k_set, n_tail))
        p_down_tail = set(universe[-n_tail:])
        k_down = sum(g in p_down_tail for g in inset)
        p_down = float(hypergeom.sf(k_down - 1, n_univ, k_set, n_tail))
        rows.append(
            {
                "geneset": name,
                "set_in_universe": k_set,
                "set_size_gmt": len(genes),
                "tail_n": n_tail,
                "universe_n": n_univ,
                "overlap_up_tail": k_hit,
                "hypergeom_p_up": p_up,
                "overlap_down_tail": k_down,
                "hypergeom_p_down": p_down,
            }
        )
    out = pd.DataFrame(rows)
    out["fdr_up"] = bh(out["hypergeom_p_up"].to_numpy())
    out["fdr_down"] = bh(out["hypergeom_p_down"].to_numpy())
    out = out.sort_values(["hypergeom_p_up", "overlap_up_tail"], ascending=[True, False])
    out.insert(0, "overlap_up_rank", np.arange(1, len(out) + 1))
    return out.reset_index(drop=True)


def bh(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        adj[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0, 1)
    return out


def set_effect(de: pd.DataFrame, genes: list[str]) -> dict:
    present = [g for g in genes if g in de.index]
    measured = de.loc[present] if present else de.iloc[0:0]
    finite = measured[np.isfinite(measured["log2FoldChange"])]
    if finite.empty:
        return {
            "n_panel": len(genes),
            "n_tested": 0,
            "mean_lfc": np.nan,
            "median_lfc": np.nan,
            "mean_lfc_shrunk": np.nan,
            "frac_lfc_pos": np.nan,
            "n_lfc_gt_0.25": 0,
            "n_lfc_lt_-0.25": 0,
            "mw_p": np.nan,
            "mw_median_minus_background": np.nan,
        }
    background = de.loc[~de.index.isin(finite.index), "log2FoldChange"].dropna()
    mw = mannwhitneyu(finite["log2FoldChange"], background, alternative="two-sided")
    return {
        "n_panel": len(genes),
        "n_tested": int(len(finite)),
        "mean_lfc": float(finite["log2FoldChange"].mean()),
        "median_lfc": float(finite["log2FoldChange"].median()),
        "mean_lfc_shrunk": float(finite["log2FoldChange_shrunk"].mean()),
        "frac_lfc_pos": float((finite["log2FoldChange"] > 0).mean()),
        "n_lfc_gt_0.25": int((finite["log2FoldChange"] >= MEAN_LFC_GATE).sum()),
        "n_lfc_lt_-0.25": int((finite["log2FoldChange"] <= -MEAN_LFC_GATE).sum()),
        "mw_p": float(mw.pvalue),
        "mw_median_minus_background": float(
            finite["log2FoldChange"].median() - background.median()
        ),
    }


def gene_rows(des: dict[str, pd.DataFrame], genes: list[str], panel: str) -> pd.DataFrame:
    rows = []
    for gene in genes:
        row = {"panel": panel, "gene": gene}
        for key, de in des.items():
            if gene not in de.index:
                row[f"{key}_tested"] = False
                continue
            row[f"{key}_tested"] = True
            rec = de.loc[gene]
            row[f"{key}_baseMean"] = rec["baseMean"]
            row[f"{key}_log2FC"] = rec["log2FoldChange"]
            row[f"{key}_log2FC_shrunk"] = rec["log2FoldChange_shrunk"]
            row[f"{key}_stat"] = rec["stat"]
            row[f"{key}_padj"] = rec["padj"]
        rows.append(row)
    return pd.DataFrame(rows)


def load_author_fpkm(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    totals: dict[str, list[float]] = {}
    with opener(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cols = {name: i for i, name in enumerate(header)}
        need = [
            "gene_short_name",
            "MCF7_CLDN4KO_FPKM (fpkm)",
            "MCF7_WT_FPKM (fpkm)",
            "T47D_CLDN4KO_FPKM (fpkm)",
            "T47D_WT_FPKM (fpkm)",
        ]
        idx = [cols[n] for n in need]
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            name = parts[idx[0]].strip().strip('"')
            if not name or name == "-" or "," in name:
                continue
            vals = []
            for j in idx[1:]:
                raw = parts[j]
                vals.append(float(raw) if raw not in {"", "NA", "nan"} else 0.0)
            if name not in totals:
                totals[name] = [0.0, 0.0, 0.0, 0.0]
            for i, v in enumerate(vals):
                totals[name][i] += v
    df = pd.DataFrame.from_dict(
        totals,
        orient="index",
        columns=["MCF7_KO", "MCF7_WT", "T47D_KO", "T47D_WT"],
    )
    pc = 0.1
    df["MCF7_log2FC"] = np.log2((df["MCF7_KO"] + pc) / (df["MCF7_WT"] + pc))
    df["T47D_log2FC"] = np.log2((df["T47D_KO"] + pc) / (df["T47D_WT"] + pc))
    df["mean_log2FC"] = (df["MCF7_log2FC"] + df["T47D_log2FC"]) / 2
    return df


def fmt_p(value) -> str:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "NA"
    value = float(value)
    if value == 0:
        return "<1e-300"
    if value < 1e-3:
        return f"{value:.2e}"
    return f"{value:.3f}"


def fmt_n(value, digits=2) -> str:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "NA"
    return f"{float(value):.{digits}f}"


def de_counts(de: pd.DataFrame) -> dict:
    finite = de[np.isfinite(de["pvalue"])]
    return {
        "n_tested": int(len(de)),
        "n_p05": int((finite["pvalue"] < 0.05).sum()),
        "n_p05_up": int(((finite["pvalue"] < 0.05) & (finite["log2FoldChange"] > 0)).sum()),
        "n_p05_down": int(((finite["pvalue"] < 0.05) & (finite["log2FoldChange"] < 0)).sum()),
        "n_fdr10": int((de["padj"] < 0.10).sum()) if de["padj"].notna().any() else 0,
        "n_fdr05": int((de["padj"] < 0.05).sum()) if de["padj"].notna().any() else 0,
        "n_fdr05_up": int(((de["padj"] < 0.05) & (de["log2FoldChange"] > 0)).sum()),
        "n_fdr05_down": int(((de["padj"] < 0.05) & (de["log2FoldChange"] < 0)).sum()),
    }


def plot_nes(gsea: pd.DataFrame, path: Path) -> None:
    df = gsea.sort_values("NES")
    colors = []
    for name in df["geneset"]:
        if name in {IFN_ALPHA, IFN_GAMMA}:
            colors.append("#b2182b")
        elif name in {
            "Inflammatory Response",
            "IL-6/JAK/STAT3 Signaling",
            "TNF-alpha Signaling via NF-kB",
            "Allograft Rejection",
            "Complement",
        }:
            colors.append("#ef8a62")
        elif name in {"Cholesterol Homeostasis", "Fatty Acid Metabolism", "Bile Acid Metabolism"}:
            colors.append("#2166ac")
        elif name == "Apical Junction":
            colors.append("#4dac26")
        else:
            colors.append("#bdbdbd")
    fig, ax = plt.subplots(figsize=(8.2, 10.5))
    ax.barh(df["geneset"], df["NES"], color=colors, height=0.78)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("NES, CLDN4 KO versus WT (cell line adjusted)")
    ax.set_title("MSigDB Hallmark enrichment in GSE207704")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_concordance(gene_table: pd.DataFrame, path: Path) -> None:
    df = gene_table[gene_table["panel"].isin(["IFN_alpha", "IFN_gamma", "MHC-I/APM"])].copy()
    df = df[df["T47D_tested"] & df["MCF7_tested"]]
    priority = {"IFN_alpha": 0, "IFN_gamma": 1, "MHC-I/APM": 2}
    df["_priority"] = df["panel"].map(priority)
    df = df.sort_values("_priority").drop_duplicates("gene", keep="last")
    fig, ax = plt.subplots(figsize=(6.4, 6.2))
    colors = {"IFN_alpha": "#b2182b", "IFN_gamma": "#ef8a62", "MHC-I/APM": "#2166ac"}
    for panel, block in df.groupby("panel"):
        ax.scatter(
            block["MCF7_log2FC"],
            block["T47D_log2FC"],
            s=28,
            c=colors[panel],
            label=panel.replace("_", "-"),
            alpha=0.85,
            linewidths=0,
        )
    lim = np.nanmax(np.abs(df[["MCF7_log2FC", "T47D_log2FC"]].to_numpy()))
    lim = max(1.0, float(lim) * 1.05)
    ax.plot([-lim, lim], [-lim, lim], color="#888888", lw=0.7)
    ax.axhline(0, color="#cccccc", lw=0.5)
    ax.axvline(0, color="#cccccc", lw=0.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("MCF7 log2FC (KO / WT)")
    ax.set_ylabel("T47D log2FC (KO / WT)")
    ax.set_title("IFN and MHC-I/APM genes")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def hallmark_line(gsea: pd.DataFrame, name: str) -> pd.Series:
    hit = gsea[gsea["geneset"] == name]
    if hit.empty:
        raise KeyError(name)
    return hit.iloc[0]


def relation_phrase(label: str) -> str:
    return {
        "up": "up in the knockout (same direction as the private CLDN4-loss IFN/APM increase)",
        "down": "down in the knockout (opposite the private CLDN4-loss IFN/APM increase)",
        "null": "null",
        "weak": "a weak coordinated shift (FDR < 0.05, mean |log2FC| < 0.25)",
        "suggestive": "a suggestive mean shift that does not pass the Hallmark FDR gate",
        "discordant": "discordant between T47D and MCF7",
        "line-dependent": "driven by one cell line",
        "unresolved": "unresolved, because the mean shift and the enrichment sign disagree",
    }[label]


def write_results(path: Path, ctx: dict) -> None:
    qc = ctx["qc"]
    comb = ctx["de_counts"]["combined"]
    a = ctx["sets"][IFN_ALPHA]
    g = ctx["sets"][IFN_GAMMA]
    apm = ctx["sets"]["MHC-I/APM"]
    tj = ctx["sets"]["Tight junction (excluding CLDN4)"]
    cldn4 = ctx["cldn4"]
    lines = []
    lines.append("# GSE207704 CLDN4 knockout: IFN, MHC-I/APM, and tight junctions")
    lines.append("")
    lines.append("## Result")
    lines.append("")
    lines.append(ctx["lead"])
    lines.append("")
    lines.append(
        "This is a breast-cancer cell-line knockout (T47D and MCF7), not a lung-cancer knockdown. "
        "It does not re-estimate the private CLDN4-knockdown IFN/APM result."
    )
    lines.append("")
    lines.append("## Question")
    lines.append("")
    lines.append(
        "A private CLDN4 knockdown was described as upregulating IFN/APM "
        "(347 genes up, 736 down; 49/52 IFN/APM genes shared with an SKB264 signature, r = 0.52). "
        "No public lung-cancer CLDN4 knockdown or knockout transcriptome is deposited. "
        "The nearest public CLDN4 loss-of-function RNA-seq series is GSE207704 "
        "(Kage et al., Breast Cancer Research 2023, PMID 37059993): CRISPR knockout of CLDN4 "
        "in T47D and MCF7, two replicates each, compared with wild type."
    )
    lines.append("")
    lines.append("## Data")
    lines.append("")
    lines.append(
        "Eight single-end NovaSeq libraries, SRA SRR20029118–SRR20029125 "
        "(BioProject PRJNA856719). SRA runinfo reports a mean read length of 51 bp "
        "and about 44–47 million reads per library. "
        "The paper text describes index-trimmed 100 bp reads; the public runs are 51 bp."
    )
    lines.append("")
    lines.append(
        "Reads were extracted with sratoolkit 3.4.1 `fastq-dump` in chunks of at most 5 million "
        "spots (1 million when a chunk aborted), then concatenated. Full-file `fasterq-dump` and "
        "`fastq-dump` runs segfaulted in this environment; the merged files contain the same "
        "spot counts as the SRA runinfo. "
        "Reads were quantified with kallisto 0.52.0 against Ensembl release 116 GRCh38 cDNA "
        "(`Homo_sapiens.GRCh38.cdna.all.fa.gz`), k = 31, single-end fragment prior "
        "`-l 200 -s 30`. Transcript estimated counts were summed to gene symbols. "
        "A low transcript-level unique-mapping rate is expected for 51 bp reads across isoforms; "
        "the differential-expression input is the gene-level sum of kallisto estimated counts. "
        "When two Ensembl genes shared a symbol, the gene with more total reads was kept "
        f"({ctx['n_symbol_collisions']} symbols had a dropped gene id)."
    )
    lines.append("")
    lines.append("Pseudoalignment rate:")
    lines.append("")
    lines.append("| Sample | Genotype | Reads | Pseudoaligned |")
    lines.append("|---|---|---:|---:|")
    for _, row in qc.iterrows():
        lines.append(
            f"| {row['sample']} | {row['cell_line']} {row['condition']} | "
            f"{int(row['n_processed']):,} | {float(row['p_pseudoaligned']):.1f}% |"
        )
    lines.append("")
    lines.append(
        "The GEO supplementary file `GSE207704_CLDN4_RNAseq.txt.gz` is a Cufflinks-style "
        "FPKM matrix with four columns (the two replicates are already averaged) and it does not "
        "contain HLA-A, HLA-B, B2M, TAP1, TAP2, or MX1. It is used only as a sensitivity check "
        "on the genes it does contain."
    )
    lines.append("")
    lines.append("## Statistics")
    lines.append("")
    lines.append(
        "Counts were analyzed with PyDESeq2 0.5.4. The primary model is "
        "`~ cell_line + condition` on all eight samples, log2 fold change = KO / WT, "
        "wild type as the condition reference. T47D and MCF7 were also fit separately "
        "(`~ condition`, two versus two). Cook's filtering was off, matching DESeq2 when "
        "a group has fewer than three replicates. Wald statistics used maximum-likelihood "
        "log2 fold changes. apeglm shrinkage was applied only to the reported shrunk log2 fold change."
    )
    lines.append("")
    lines.append(
        "Hallmark enrichment is preranked GSEA (gseapy 1.3.1, fgsea multilevel, "
        f"{GSEA_PERM} permutations, seed {GSEA_SEED}) on the Wald statistic. "
        "The gene sets are Enrichr `MSigDB_Hallmark_2020`, the Hallmark collection of "
        "Liberzon et al. (Cell Systems 2015), with the large sets capped at 200 genes "
        "in that export. FDR is across the 50 Hallmark sets. "
        "A second rank is the hypergeometric overlap of each set with the top 10% of genes "
        "by the same Wald statistic."
    )
    lines.append("")
    lines.append(
        "A set is called **up** or **down** only when three things agree: "
        f"the combined-model mean |log2FC| is at least {MEAN_LFC_GATE}, "
        f"the Hallmark (or panel) GSEA FDR is below {FDR_GATE} with a matching NES sign, "
        "and the median log2FC is in that same direction in both cell lines. "
        "Opposite medians of at least 0.25 log2FC are **discordant**. "
        "A significant GSEA with a smaller mean is **weak**. "
        "A mean past the gate without GSEA FDR support is **suggestive**. "
        "Otherwise the call is **null**. These gates were written into the script before the counts were scored."
    )
    lines.append("")
    lines.append("## Knockout check")
    lines.append("")
    lines.append(
        f"CLDN4 combined log2FC = {fmt_n(cldn4['combined_log2FC'])} "
        f"(shrunk {fmt_n(cldn4['combined_shrunk'])}, Wald p = {fmt_p(cldn4['combined_p'])}, "
        f"FDR = {fmt_p(cldn4['combined_fdr'])}). "
        f"T47D log2FC = {fmt_n(cldn4['T47D_log2FC'])}; "
        f"MCF7 log2FC = {fmt_n(cldn4['MCF7_log2FC'])}. "
        f"Author FPKM log2FC (pseudocount 0.1) was {fmt_n(cldn4['author_T47D'])} in T47D and "
        f"{fmt_n(cldn4['author_MCF7'])} in MCF7. "
        + ctx["cldn4_sentence"]
    )
    lines.append("")
    lines.append("## Global differential expression")
    lines.append("")
    lines.append(
        f"Combined model, {comb['n_tested']:,} genes with total count at least 10. "
        f"Uncorrected p < 0.05: {comb['n_p05_up']:,} up and {comb['n_p05_down']:,} down. "
        f"FDR < 0.05: {comb['n_fdr05_up']:,} up and {comb['n_fdr05_down']:,} down. "
        f"FDR < 0.10: {comb['n_fdr10']:,} genes. "
        "The combined model has four knockout and four wild-type libraries, so single-gene "
        "FDRs are not empty. The per-line models have two replicates and much less power."
    )
    lines.append("")
    lines.append(
        f"Spearman correlation of unshrunk log2FC, T47D versus MCF7, "
        f"genes tested in both: r = {fmt_n(ctx['spearman_all'][0], 3)}, "
        f"p = {fmt_p(ctx['spearman_all'][1])}, n = {ctx['spearman_all_n']:,}. "
        f"Among genes with combined FDR < 0.05 that were tested in both lines, "
        f"r = {fmt_n(ctx['spearman_sig'][0], 3)} "
        f"(p {fmt_p(ctx['spearman_sig'][1])}, n = {ctx['spearman_sig_n']:,}, "
        f"{fmt_n(100 * ctx['spearman_sig_same_sign'], 1)}% same sign). "
        "The genome-wide correlation is near zero because most genes do not move. "
        "The genes the combined model calls significant do move in the same direction in both lines."
    )
    lines.append("")
    lines.append("## Hallmark rank")
    lines.append("")
    lines.append(
        "NES rank 1 is the Hallmark most enriched among genes that rise in the knockout. "
        f"Interferon Alpha Response is NES rank {int(a['nes_rank'])} of {ctx['n_hallmark']} "
        f"(NES {fmt_n(a['NES'])}, FDR {fmt_p(a['fdr'])}; "
        f"overlap-up rank {int(a['overlap_up_rank'])} of {ctx['n_hallmark']}; "
        f"{int(a['overlap_up_tail'])}/{int(a['set_in_universe'])} genes in the top 10%; "
        f"{int(a['overlap_down_tail'])}/{int(a['set_in_universe'])} in the bottom 10%; "
        f"bottom-decile FDR {fmt_p(a['fdr_down'])}). "
        f"Interferon Gamma Response is NES rank {int(g['nes_rank'])} of {ctx['n_hallmark']} "
        f"(NES {fmt_n(g['NES'])}, FDR {fmt_p(g['fdr'])}; "
        f"overlap-up rank {int(g['overlap_up_rank'])} of {ctx['n_hallmark']}; "
        f"{int(g['overlap_up_tail'])}/{int(g['set_in_universe'])} genes in the top 10%; "
        f"{int(g['overlap_down_tail'])}/{int(g['set_in_universe'])} in the bottom 10%; "
        f"bottom-decile FDR {fmt_p(g['fdr_down'])})."
    )
    lines.append("")
    lines.append("| NES rank | Hallmark | NES | FDR | Overlap-up rank | Top-decile overlap | Call |")
    lines.append("|---:|---|---:|---:|---:|---:|---|")
    for row in ctx["hallmark_rows"]:
        lines.append(
            f"| {int(row['nes_rank'])} | {row['geneset']} | {fmt_n(row['NES'])} | "
            f"{fmt_p(row['fdr'])} | {int(row['overlap_up_rank'])} | "
            f"{int(row['overlap_up_tail'])}/{int(row['set_in_universe'])} | {row['call']} |"
        )
    lines.append("")
    lines.append("Full 50-set tables are in `results/gse207704/`.")
    lines.append("")
    lines.append("## Interferon and MHC-I/APM")
    lines.append("")
    lines.append("| Set | Call | Mean log2FC | Median T47D | Median MCF7 | NES | FDR | Genes log2FC > 0 |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|")
    for key in [IFN_ALPHA, IFN_GAMMA, "MHC-I/APM"]:
        row = ctx["sets"][key]
        lines.append(
            f"| {key} | {row['call']} | {fmt_n(row['mean_lfc'])} | {fmt_n(row['median_T47D'])} | "
            f"{fmt_n(row['median_MCF7'])} | {fmt_n(row['NES'])} | {fmt_p(row['fdr'])} | "
            f"{fmt_n(100 * row['frac_lfc_pos'], 0)}% of {int(row['n_tested'])} |"
        )
    lines.append("")
    lines.append(
        f"Interferon-alpha is {relation_phrase(a['call'])}, and the enrichment is "
        f"{'negative' if a['NES'] < 0 else 'positive'} "
        f"(mean log2FC {fmt_n(a['mean_lfc'])}). "
        f"Interferon-gamma is {relation_phrase(g['call'])}, and the enrichment is "
        f"{'negative' if g['NES'] < 0 else 'positive'} "
        f"(mean log2FC {fmt_n(g['mean_lfc'])}). "
        f"The MHC-I/APM panel is {relation_phrase(apm['call'])} "
        f"(mean log2FC {fmt_n(apm['mean_lfc'])}). "
        "Negative NES means the set is enriched among genes that fall in the knockout."
    )
    lines.append("")
    lines.append(
        f"MHC-I/APM genes on the panel: {apm['n_tested']} of {apm['n_panel']} were tested "
        f"(total count ≥ 10 in the combined matrix). "
        f"Mean shrunk log2FC = {fmt_n(apm['mean_lfc_shrunk'])}. "
        f"Mann–Whitney p versus other genes = {fmt_p(apm['mw_p'])} "
        f"(median minus background {fmt_n(apm['mw_median_minus_background'])}). "
        f"Panel genes absent from the tested matrix: {ctx['apm_missing'] or 'none'}."
    )
    lines.append("")
    lines.append(
        "Per-gene log2 fold changes for the interferon Hallmarks, the MHC-I/APM panel, "
        "and the tight-junction panel are in `results/gse207704/genes_ifn_apm_tj.tsv`."
    )
    lines.append("")
    if ctx["leading"]:
        lines.append("Largest absolute combined Wald statistics inside the two interferon Hallmarks:")
        lines.append("")
        lines.append("| Gene | Set membership | log2FC | Shrunk log2FC | FDR |")
        lines.append("|---|---|---:|---:|---:|")
        for row in ctx["leading"]:
            lines.append(
                f"| {row['gene']} | {row['which']} | {fmt_n(row['log2FC'])} | "
                f"{fmt_n(row['shrunk'])} | {fmt_p(row['padj'])} |"
            )
        lines.append("")
    lines.append("## Tight junctions")
    lines.append("")
    lines.append(
        f"The tight-junction panel excludes CLDN4. Call: **{tj['call']}**. "
        f"Unshrunk mean log2FC = {fmt_n(tj['mean_lfc'])}; "
        f"apeglm-shrunk mean = {fmt_n(tj['mean_lfc_shrunk'])} "
        f"(median T47D {fmt_n(tj['median_T47D'])}, median MCF7 {fmt_n(tj['median_MCF7'])}). "
        "The shrunk mean is closer to zero because several low-count claudins have large unshrunk fold changes. "
        f"NES = {fmt_n(tj['NES'])}, panel GSEA FDR = {fmt_p(tj['fdr'])} "
        "(this FDR is within the custom panels, not the Hallmark family). "
        f"{tj['n_tested']} of {tj['n_panel']} panel genes were tested. "
        f"Absent: {ctx['tj_missing'] or 'none'}."
    )
    lines.append("")
    lines.append(
        f"Hallmark Apical Junction, which is broader than the claudin/TJ list, "
        f"is NES rank {int(ctx['sets']['Apical Junction']['nes_rank'])} "
        f"(NES {fmt_n(ctx['sets']['Apical Junction']['NES'])}, "
        f"FDR {fmt_p(ctx['sets']['Apical Junction']['fdr'])}, "
        f"call {ctx['sets']['Apical Junction']['call']})."
    )
    lines.append("")
    lines.append("Epithelial context genes (not part of the TJ call):")
    lines.append("")
    lines.append("| Gene | Combined log2FC | T47D log2FC | MCF7 log2FC |")
    lines.append("|---|---:|---:|---:|")
    for row in ctx["epithelial"]:
        lines.append(
            f"| {row['gene']} | {fmt_n(row['combined'])} | {fmt_n(row['T47D'])} | {fmt_n(row['MCF7'])} |"
        )
    lines.append("")
    lines.append("## Author FPKM sensitivity")
    lines.append("")
    lines.append(
        f"On genes present in both the kallisto model and the deposited FPKM table, "
        f"Spearman correlation of the two-line mean log2FC is "
        f"r = {fmt_n(ctx['fpkm_spearman'][0], 3)} (p {fmt_p(ctx['fpkm_spearman'][1])}, "
        f"n = {ctx['fpkm_spearman_n']:,}). "
        f"Author-matrix GSEA on that mean FPKM log2FC places Interferon Alpha Response at "
        f"NES rank {ctx['fpkm_alpha_rank']} (NES {fmt_n(ctx['fpkm_alpha_nes'])}) and "
        f"Interferon Gamma Response at NES rank {ctx['fpkm_gamma_rank']} "
        f"(NES {fmt_n(ctx['fpkm_gamma_nes'])})."
    )
    lines.append("")
    lines.append(
        f"Of the interferon-alpha Hallmark, {ctx['fpkm_alpha_missing_n']} genes are missing "
        f"from the deposited table. Of interferon-gamma, {ctx['fpkm_gamma_missing_n']} are missing. "
        f"Of the MHC-I/APM panel, {ctx['fpkm_apm_missing_n']} are missing"
        + (f" ({ctx['fpkm_apm_missing']})" if ctx["fpkm_apm_missing"] else "")
        + ". A conclusion drawn only from that table would not have measured those genes."
    )
    lines.append("")
    lines.append("## Limits")
    lines.append("")
    lines.append(
        "Two biological replicates per genotype. The cultures have no immune cells, so this is "
        "tumor-cell-intrinsic transcription only. The perturbation is a coding-sequence knockout "
        "in two ER-positive breast lines, not a knockdown in lung cancer. "
        "Ensembl 116 cDNA does not include every noncoding transcript. "
        "Single-end 51 bp reads limit isoform assignment; gene-level sums are the endpoint. "
        "The private knockdown counts and the SKB264 correlation were not recomputed here."
    )
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `scripts/quantify_gse207704.sh` — kallisto index and quantification")
    lines.append("- `scripts/analyze_gse207704.py` — DESeq2, GSEA, tables, this document")
    lines.append("- `resources/MSigDB_Hallmark_2020.gmt` — Hallmark sets used for the ranks")
    lines.append("- `results/gse207704/` — count matrix, DE tables, GSEA tables, figures")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def build_set_record(name, genes, de_c, de_t, de_m, gsea_c, gsea_custom, overlap) -> dict:
    eff = set_effect(de_c, genes)
    eff_t = set_effect(de_t, genes)
    eff_m = set_effect(de_m, genes)
    if name in set(gsea_c["geneset"]):
        grow = hallmark_line(gsea_c, name)
        nes, fdr, nes_rank = float(grow["NES"]), float(grow["fdr"]), int(grow["nes_rank"])
    else:
        grow = hallmark_line(gsea_custom, name)
        nes, fdr, nes_rank = float(grow["NES"]), float(grow["fdr"]), int(grow["nes_rank"])
    ov = overlap[overlap["geneset"] == name]
    if ov.empty:
        ov_rank, ov_hit, ov_set, ov_down, fdr_down = -1, 0, 0, 0, float("nan")
    else:
        ov_rank = int(ov.iloc[0]["overlap_up_rank"])
        ov_hit = int(ov.iloc[0]["overlap_up_tail"])
        ov_set = int(ov.iloc[0]["set_in_universe"])
        ov_down = int(ov.iloc[0]["overlap_down_tail"])
        fdr_down = float(ov.iloc[0]["fdr_down"])
    call = classify(eff["mean_lfc"], nes, fdr, eff_t["median_lfc"], eff_m["median_lfc"])
    return {
        "geneset": name,
        "call": call,
        "mean_lfc": eff["mean_lfc"],
        "median_lfc": eff["median_lfc"],
        "mean_lfc_shrunk": eff["mean_lfc_shrunk"],
        "median_T47D": eff_t["median_lfc"],
        "median_MCF7": eff_m["median_lfc"],
        "mean_T47D": eff_t["mean_lfc"],
        "mean_MCF7": eff_m["mean_lfc"],
        "frac_lfc_pos": eff["frac_lfc_pos"],
        "n_tested": eff["n_tested"],
        "n_panel": eff["n_panel"],
        "NES": nes,
        "fdr": fdr,
        "nes_rank": nes_rank,
        "overlap_up_rank": ov_rank,
        "overlap_up_tail": ov_hit,
        "overlap_down_tail": ov_down,
        "fdr_down": fdr_down,
        "set_in_universe": ov_set,
        "mw_p": eff["mw_p"],
        "mw_median_minus_background": eff["mw_median_minus_background"],
        "n_up_0.25": eff["n_lfc_gt_0.25"],
        "n_down_0.25": eff["n_lfc_lt_-0.25"],
    }


def lead_sentence(sets: dict) -> str:
    a = sets[IFN_ALPHA]
    g = sets[IFN_GAMMA]
    p = sets["MHC-I/APM"]
    tj = sets["Tight junction (excluding CLDN4)"]
    if {a["call"], g["call"], p["call"]} == {"up"}:
        opener = "CLDN4 knockout raises interferon and MHC-I/APM genes in both breast lines."
    else:
        opener = "CLDN4 knockout does not upregulate interferon or MHC-I/APM genes."
    return (
        f"In GSE207704, {opener} "
        f"Interferon Alpha Response is NES rank {int(a['nes_rank'])} of 50 "
        f"(NES {a['NES']:.2f}, mean log2FC {a['mean_lfc']:.2f}, call {a['call']}). "
        f"Interferon Gamma Response is NES rank {int(g['nes_rank'])} of 50 "
        f"(NES {g['NES']:.2f}, mean log2FC {g['mean_lfc']:.2f}, call {g['call']}). "
        f"The MHC-I/APM panel is {p['call']} (mean log2FC {p['mean_lfc']:.2f}). "
        f"The tight-junction panel, with CLDN4 removed, is {tj['call']} "
        f"(unshrunk mean log2FC {tj['mean_lfc']:.2f}, shrunk mean {tj['mean_lfc_shrunk']:.2f}). "
        "Negative NES means the gene set is enriched among genes that fall after knockout."
    )


def main() -> None:
    _self_check()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quant-dir", type=Path, default=Path("/tmp/gse207704/quant"))
    parser.add_argument("--fasta", type=Path, default=Path("/tmp/gse207704/ref/cdna.fa.gz"))
    parser.add_argument(
        "--hallmark",
        type=Path,
        default=Path("resources/MSigDB_Hallmark_2020.gmt"),
    )
    parser.add_argument(
        "--author-fpkm",
        type=Path,
        default=Path("/tmp/gse207704/GSE207704_CLDN4_RNAseq.txt.gz"),
    )
    parser.add_argument("--outdir", type=Path, default=Path("results/gse207704"))
    parser.add_argument("--results-md", type=Path, default=Path("RESULTS.md"))
    args = parser.parse_args()

    outdir = args.outdir
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    print("loading transcriptome map", flush=True)
    tx2gene = parse_tx2gene(args.fasta)
    print(f"transcripts with symbols: {len(tx2gene):,}", flush=True)
    counts, tpm, qc, collisions = load_kallisto(args.quant_dir, tx2gene)
    print(f"genes: {counts.shape[0]:,}", flush=True)
    counts.to_csv(outdir / "gene_counts.tsv.gz", sep="\t")
    tpm.to_csv(outdir / "gene_tpm.tsv.gz", sep="\t")
    qc.to_csv(outdir / "sample_qc.tsv", sep="\t", index=False)
    collisions.to_csv(outdir / "symbol_collisions.tsv", sep="\t", index=False)

    meta = metadata()
    print("DESeq2 combined", flush=True)
    de_c = run_deseq(counts, meta, DESIGN_COMBINED)
    print("DESeq2 T47D", flush=True)
    de_t = run_deseq(counts, meta.loc[meta["cell_line"] == "T47D"], DESIGN_LINE)
    print("DESeq2 MCF7", flush=True)
    de_m = run_deseq(counts, meta.loc[meta["cell_line"] == "MCF7"], DESIGN_LINE)
    de_c.to_csv(outdir / "de_combined.tsv.gz", sep="\t")
    de_t.to_csv(outdir / "de_T47D.tsv.gz", sep="\t")
    de_m.to_csv(outdir / "de_MCF7.tsv.gz", sep="\t")

    hallmark = read_gmt(args.hallmark)
    custom = {
        "MHC-I/APM": APM_GENES,
        "Tight junction (excluding CLDN4)": TJ_GENES,
    }
    print("GSEA hallmark combined", flush=True)
    gsea_c = gsea_table(de_c["stat"], hallmark, min_size=15)
    print("GSEA hallmark T47D", flush=True)
    gsea_t = gsea_table(de_t["stat"], hallmark, min_size=15)
    print("GSEA hallmark MCF7", flush=True)
    gsea_m = gsea_table(de_m["stat"], hallmark, min_size=15)
    print("GSEA custom panels", flush=True)
    gsea_custom = gsea_table(de_c["stat"], custom, min_size=8)
    gsea_c.to_csv(outdir / "hallmark_gsea_combined.tsv", sep="\t", index=False)
    gsea_t.to_csv(outdir / "hallmark_gsea_T47D.tsv", sep="\t", index=False)
    gsea_m.to_csv(outdir / "hallmark_gsea_MCF7.tsv", sep="\t", index=False)
    gsea_custom.to_csv(outdir / "custom_gsea_combined.tsv", sep="\t", index=False)

    print("overlap", flush=True)
    overlap_h = overlap_table(de_c["stat"], hallmark)
    overlap_h.insert(1, "family", "hallmark")
    overlap_c = overlap_table(de_c["stat"], custom)
    overlap_c.insert(1, "family", "custom")
    overlap = pd.concat([overlap_h, overlap_c], ignore_index=True)
    overlap.to_csv(outdir / "overlap_top_decile.tsv", sep="\t", index=False)

    set_names = {
        IFN_ALPHA: hallmark[IFN_ALPHA],
        IFN_GAMMA: hallmark[IFN_GAMMA],
        "MHC-I/APM": APM_GENES,
        "Tight junction (excluding CLDN4)": TJ_GENES,
        "Apical Junction": hallmark["Apical Junction"],
        "Cholesterol Homeostasis": hallmark["Cholesterol Homeostasis"],
        "Inflammatory Response": hallmark["Inflammatory Response"],
    }
    sets = {}
    for name, genes in set_names.items():
        sets[name] = build_set_record(name, genes, de_c, de_t, de_m, gsea_c, gsea_custom, overlap)
    pd.DataFrame(sets.values()).to_csv(outdir / "set_summary.tsv", sep="\t", index=False)

    gene_table = pd.concat(
        [
            gene_rows({"combined": de_c, "T47D": de_t, "MCF7": de_m}, hallmark[IFN_ALPHA], "IFN_alpha"),
            gene_rows({"combined": de_c, "T47D": de_t, "MCF7": de_m}, hallmark[IFN_GAMMA], "IFN_gamma"),
            gene_rows({"combined": de_c, "T47D": de_t, "MCF7": de_m}, APM_GENES, "MHC-I/APM"),
            gene_rows({"combined": de_c, "T47D": de_t, "MCF7": de_m}, TJ_GENES, "tight_junction"),
            gene_rows({"combined": de_c, "T47D": de_t, "MCF7": de_m}, EPITHELIAL_CONTEXT, "epithelial"),
        ],
        ignore_index=True,
    )
    gene_table.to_csv(outdir / "genes_ifn_apm_tj.tsv", sep="\t", index=False)

    # Leading genes by |Wald stat| inside either IFN hallmark.
    ifn_genes = sorted(set(hallmark[IFN_ALPHA]) | set(hallmark[IFN_GAMMA]))
    ifn_de = de_c.loc[de_c.index.isin(ifn_genes)].copy()
    ifn_de["abs_stat"] = ifn_de["stat"].abs()
    ifn_de = ifn_de.sort_values("abs_stat", ascending=False).head(15)
    leading = []
    for gene, rec in ifn_de.iterrows():
        which = []
        if gene in hallmark[IFN_ALPHA]:
            which.append("alpha")
        if gene in hallmark[IFN_GAMMA]:
            which.append("gamma")
        leading.append(
            {
                "gene": gene,
                "which": "+".join(which),
                "log2FC": float(rec["log2FoldChange"]),
                "shrunk": float(rec["log2FoldChange_shrunk"]),
                "padj": float(rec["padj"]) if pd.notna(rec["padj"]) else np.nan,
            }
        )

    shared = de_t.index.intersection(de_m.index)
    rho_all = spearmanr(de_t.loc[shared, "log2FoldChange"], de_m.loc[shared, "log2FoldChange"])
    sig = de_c.index[de_c["padj"] < 0.05].intersection(shared)
    if len(sig) >= 5:
        rho_sig = spearmanr(de_t.loc[sig, "log2FoldChange"], de_m.loc[sig, "log2FoldChange"])
        same_sign = float(
            ((de_t.loc[sig, "log2FoldChange"] > 0) == (de_m.loc[sig, "log2FoldChange"] > 0)).mean()
        )
        rho_sig_stat = (float(rho_sig.statistic), float(rho_sig.pvalue))
        rho_sig_n = int(len(sig))
    else:
        rho_sig_stat = (float("nan"), float("nan"))
        same_sign = float("nan")
        rho_sig_n = int(len(sig))

    author = load_author_fpkm(args.author_fpkm)
    both = author.index.intersection(de_c.index)
    # Author mean log2FC vs kallisto combined unshrunk log2FC.
    rho_fpkm = spearmanr(author.loc[both, "mean_log2FC"], de_c.loc[both, "log2FoldChange"])
    expressed = author[author[["MCF7_KO", "MCF7_WT", "T47D_KO", "T47D_WT"]].max(axis=1) >= 1]
    author_rank = expressed["mean_log2FC"]
    print("GSEA author FPKM", flush=True)
    gsea_fpkm = gsea_table(author_rank, hallmark, min_size=15)
    gsea_fpkm.to_csv(outdir / "hallmark_gsea_author_fpkm.tsv", sep="\t", index=False)

    def missing(genes, index):
        return [g for g in genes if g not in index]

    apm_missing = missing(APM_GENES, de_c.index)
    tj_missing = missing(TJ_GENES, de_c.index)

    def lfc_of(de, gene):
        if gene not in de.index:
            return np.nan
        return float(de.loc[gene, "log2FoldChange"])

    cldn4 = {
        "combined_log2FC": lfc_of(de_c, "CLDN4"),
        "combined_shrunk": float(de_c.loc["CLDN4", "log2FoldChange_shrunk"]) if "CLDN4" in de_c.index else np.nan,
        "combined_p": float(de_c.loc["CLDN4", "pvalue"]) if "CLDN4" in de_c.index else np.nan,
        "combined_fdr": float(de_c.loc["CLDN4", "padj"]) if "CLDN4" in de_c.index else np.nan,
        "T47D_log2FC": lfc_of(de_t, "CLDN4"),
        "MCF7_log2FC": lfc_of(de_m, "CLDN4"),
        "author_T47D": float(author.loc["CLDN4", "T47D_log2FC"]) if "CLDN4" in author.index else np.nan,
        "author_MCF7": float(author.loc["CLDN4", "MCF7_log2FC"]) if "CLDN4" in author.index else np.nan,
    }
    if cldn4["T47D_log2FC"] < 0 and cldn4["MCF7_log2FC"] < 0:
        cldn4_sentence = "CLDN4 RNA is lower in the knockout in both lines, so the sample labels and the quantification agree with the genotype."
    else:
        cldn4_sentence = "CLDN4 RNA is not lower in both knockout lines. Treat every downstream call as failed QC."

    epithelial = []
    for gene in EPITHELIAL_CONTEXT:
        epithelial.append(
            {
                "gene": gene,
                "combined": lfc_of(de_c, gene),
                "T47D": lfc_of(de_t, gene),
                "MCF7": lfc_of(de_m, gene),
            }
        )

    # Hallmark table: IFN sets, focus sets, and the five most positive and five most negative NES.
    show = set(HALLMARK_FOCUS)
    gsea_sorted = gsea_c.sort_values("NES", ascending=False)
    show.update(gsea_sorted.head(5)["geneset"])
    show.update(gsea_sorted.tail(5)["geneset"])
    hallmark_rows = []
    for _, row in gsea_c.sort_values("nes_rank").iterrows():
        if row["geneset"] not in show:
            continue
        genes = hallmark[row["geneset"]]
        rec = build_set_record(row["geneset"], genes, de_c, de_t, de_m, gsea_c, gsea_custom, overlap)
        hallmark_rows.append(rec)
    # Keep NES rank order.
    hallmark_rows = sorted(hallmark_rows, key=lambda r: r["nes_rank"])

    alpha_fp = hallmark_line(gsea_fpkm, IFN_ALPHA)
    gamma_fp = hallmark_line(gsea_fpkm, IFN_GAMMA)

    plot_nes(gsea_c, figdir / "hallmark_nes.png")
    plot_concordance(gene_table, figdir / "ifn_line_concordance.png")

    ctx = {
        "qc": qc,
        "n_symbol_collisions": int(collisions["symbol"].nunique()) if len(collisions) else 0,
        "de_counts": {
            "combined": de_counts(de_c),
            "T47D": de_counts(de_t),
            "MCF7": de_counts(de_m),
        },
        "sets": sets,
        "cldn4": cldn4,
        "cldn4_sentence": cldn4_sentence,
        "spearman_all": (float(rho_all.statistic), float(rho_all.pvalue)),
        "spearman_all_n": int(len(shared)),
        "spearman_sig": rho_sig_stat,
        "spearman_sig_n": rho_sig_n,
        "spearman_sig_same_sign": same_sign,
        "fpkm_spearman": (float(rho_fpkm.statistic), float(rho_fpkm.pvalue)),
        "fpkm_spearman_n": int(len(both)),
        "fpkm_alpha_rank": int(alpha_fp["nes_rank"]),
        "fpkm_alpha_nes": float(alpha_fp["NES"]),
        "fpkm_gamma_rank": int(gamma_fp["nes_rank"]),
        "fpkm_gamma_nes": float(gamma_fp["NES"]),
        "fpkm_alpha_missing_n": len(missing(hallmark[IFN_ALPHA], author.index)),
        "fpkm_gamma_missing_n": len(missing(hallmark[IFN_GAMMA], author.index)),
        "fpkm_apm_missing_n": len(missing(APM_GENES, author.index)),
        "fpkm_apm_missing": ", ".join(missing(APM_GENES, author.index)),
        "apm_missing": ", ".join(apm_missing),
        "tj_missing": ", ".join(tj_missing),
        "n_hallmark": int(len(gsea_c)),
        "hallmark_rows": hallmark_rows,
        "leading": leading,
        "epithelial": epithelial,
        "lead": lead_sentence(sets),
    }
    # JSON-safe copy of the summary.
    summary = {
        "lead": ctx["lead"],
        "cldn4": cldn4,
        "sets": sets,
        "de_counts": ctx["de_counts"],
        "spearman_lines": {"r": ctx["spearman_all"][0], "p": ctx["spearman_all"][1], "n": ctx["spearman_all_n"]},
        "spearman_vs_author_fpkm": {
            "r": ctx["fpkm_spearman"][0],
            "p": ctx["fpkm_spearman"][1],
            "n": ctx["fpkm_spearman_n"],
        },
        "kallisto": "0.52.0",
        "ensembl": "116",
        "pydeseq2": "0.5.4",
        "gseapy": gseapy.__version__,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    write_results(args.results_md, ctx)
    print(ctx["lead"])
    print(f"wrote {args.results_md}")


if __name__ == "__main__":
    main()
