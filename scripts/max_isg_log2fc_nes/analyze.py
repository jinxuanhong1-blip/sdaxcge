#!/usr/bin/env python3
"""Largest Hallmark-ISG log2 fold change and NES, labeled by cell line.

Recomputes DESeq2 and pre-ranked GSEA from public GEO count matrices.
Does not reuse previously published effect sizes.

Contrasts
---------
GSE294709 HCT116 Ku80 (AID degron and Ku86 flox) and the HEK293 Ku70
arm of the same series, kept as a different cell line.
GSE285698 HCT116 DNA-PKcs (PRKDC) knockout, normoxia and CoCl2.
GSE84986 MCF-7 TP53BP1 knockout, untreated and 5 Gy.

ISG set
-------
Union of MSigDB Hallmark Interferon Alpha Response and Interferon Gamma
Response (Enrichr library MSigDB_Hallmark_2020). The largest log2FC is
the maximum DESeq2 log2FoldChange in that union among genes with
baseMean >= 10. NES is pre-ranked GSEA on the Wald statistic.
"""

from __future__ import annotations

import json
import tarfile
import gzip
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "max_isg_log2fc_nes"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
CACHE = Path("/tmp/max_isg")

ALIASES = {
    "DDX58": ["DDX58", "RIGI"],
    "IFIH1": ["IFIH1", "MDA5"],
    "CGAS": ["CGAS", "MB21D1"],
    "STING1": ["STING1", "TMEM173"],
    "WARS1": ["WARS1", "WARS"],
}

# GEO sample descriptions for GSE285698 (library id -> genotype, condition, rep).
HCT_LIB = {
    "A1_S50": ("WT", "normoxia", "rep1"),
    "B1_S56": ("WT", "normoxia", "rep2"),
    "C1_S62": ("WT", "normoxia", "rep3"),
    "A2_S51": ("KO", "normoxia", "rep1"),
    "B2_S57": ("KO", "normoxia", "rep2"),
    "C2_S63": ("KO", "normoxia", "rep3"),
    "A4_S53": ("WT", "hypoxia", "rep1"),
    "B4_S59": ("WT", "hypoxia", "rep2"),
    "C4_S65": ("WT", "hypoxia", "rep3"),
    "A5_S54": ("KO", "hypoxia", "rep1"),
    "B5_S60": ("KO", "hypoxia", "rep2"),
    "C5_S66": ("KO", "hypoxia", "rep3"),
}

MCF7_TITLE = {
    "GSM2255508": "WT_untreated_1",
    "GSM2255520": "WT_untreated_2",
    "GSM2255532": "WT_untreated_3",
    "GSM2255509": "WT_IR_1",
    "GSM2255521": "WT_IR_2",
    "GSM2255533": "WT_IR_3",
    "GSM2255514": "KO1_untreated_1",
    "GSM2255526": "KO1_untreated_2",
    "GSM2255538": "KO1_untreated_3",
    "GSM2255515": "KO1_IR_1",
    "GSM2255527": "KO1_IR_2",
    "GSM2255539": "KO1_IR_3",
    "GSM2255517": "KO2_untreated_1",
    "GSM2255529": "KO2_untreated_2",
    "GSM2255541": "KO2_untreated_3",
    "GSM2255518": "KO2_IR_1",
    "GSM2255530": "KO2_IR_2",
    "GSM2255542": "KO2_IR_3",
}


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"download {dest.name}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "max-isg/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    return dest


def ensure_inputs() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    _download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294709/suppl/GSE294709_RNAseq-counts-on-gene-AMPSEQ.txt.gz",
        CACHE / "ampseq.txt.gz",
    )
    _download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294709/suppl/GSE294709_est_counts_genes_kallisto_Genome_Center.txt.gz",
        CACHE / "kallisto.txt.gz",
    )
    _download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285698/suppl/GSE285698_raw_counts.txt.gz",
        CACHE / "dnapk_counts.txt.gz",
    )
    _download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84986/suppl/GSE84986_RAW.tar",
        CACHE / "GSE84986_RAW.tar",
    )
    _download(
        "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
        CACHE / "hgnc_complete_set.txt",
    )


def collapse_counts(df: pd.DataFrame, symbol_col: str) -> pd.DataFrame:
    d = df.copy()
    d[symbol_col] = d[symbol_col].astype(str)
    d = d[d[symbol_col].notna() & ~d[symbol_col].isin(["", "nan", "None"])]
    num = d.drop(columns=[symbol_col]).apply(pd.to_numeric, errors="coerce").fillna(0)
    num = num.clip(lower=0).round().astype(int)
    num.index = d[symbol_col].values
    num = num.groupby(level=0).sum()
    return num.loc[num.sum(axis=1) > 0]


def load_amp() -> pd.DataFrame:
    raw = pd.read_csv(CACHE / "ampseq.txt.gz", sep="\t")
    return collapse_counts(raw, "gene_names")


def load_hek() -> pd.DataFrame:
    raw = pd.read_csv(CACHE / "kallisto.txt.gz", sep="\t")
    return collapse_counts(raw, "gene_name")


def load_dnapk() -> pd.DataFrame:
    raw = pd.read_csv(CACHE / "dnapk_counts.txt.gz", sep="\t")
    gene_col = raw.columns[0]
    symbols = []
    for raw_id in raw[gene_col].astype(str):
        symbols.append(raw_id.split("|")[-1] if "|" in raw_id else raw_id)
    raw = raw.drop(columns=[gene_col])
    renamed = {}
    for col in raw.columns:
        parts = str(col).split("_")
        key = "_".join(parts[:2]) if len(parts) >= 2 else str(col)
        if key not in HCT_LIB:
            raise KeyError(col)
        geno, cond, rep = HCT_LIB[key]
        renamed[col] = f"{geno}_{cond}_{rep}"
    raw = raw.rename(columns=renamed)
    raw.insert(0, "symbol", symbols)
    return collapse_counts(raw, "symbol")


def load_hgnc() -> dict[str, str]:
    hgnc = pd.read_csv(
        CACHE / "hgnc_complete_set.txt",
        sep="\t",
        usecols=["symbol", "ensembl_gene_id"],
        dtype=str,
    )
    hgnc = hgnc.dropna(subset=["ensembl_gene_id", "symbol"]).drop_duplicates("ensembl_gene_id")
    return dict(zip(hgnc["ensembl_gene_id"], hgnc["symbol"]))


def load_mcf7(ens_map: dict[str, str]) -> pd.DataFrame:
    columns: dict[str, dict[str, float]] = {}
    with tarfile.open(CACHE / "GSE84986_RAW.tar") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            gsm = Path(member.name).name.split("_", 1)[0]
            if gsm not in MCF7_TITLE:
                continue
            text = gzip.decompress(tar.extractfile(member).read()).decode("utf-8", "replace")
            vals: dict[str, float] = {}
            for line in text.splitlines():
                if not line or line.startswith("#"):
                    continue
                eid, _, rest = line.partition("\t")
                eid = eid.split(".")[0]
                try:
                    vals[eid] = float(rest.strip() or 0)
                except ValueError:
                    continue
            columns[MCF7_TITLE[gsm]] = vals
    symbols = {}
    for sample, vals in columns.items():
        for eid, count in vals.items():
            sym = ens_map.get(eid)
            if not sym:
                continue
            symbols.setdefault(sym, {}).setdefault(sample, 0.0)
            symbols[sym][sample] += count
    mat = pd.DataFrame(symbols).T.fillna(0.0)
    mat = mat.clip(lower=0).round().astype(int)
    mat = mat.loc[mat.sum(axis=1) > 0]
    # Stable column order.
    order = [c for c in MCF7_TITLE.values() if c in mat.columns]
    return mat[order]


def load_hallmark() -> dict[str, list[str]]:
    import gseapy as gp

    lib = gp.get_library("MSigDB_Hallmark_2020")
    alpha = list(lib["Interferon Alpha Response"])
    gamma = list(lib["Interferon Gamma Response"])
    return {"IFN_ALPHA": alpha, "IFN_GAMMA": gamma}


def index_lookup(index: pd.Index) -> dict[str, str]:
    return {str(i).upper(): str(i) for i in index}


def resolve_symbol(symbol: str, lookup: dict[str, str]) -> str | None:
    cands = ALIASES.get(symbol.upper(), [symbol])
    if symbol not in cands:
        cands = [symbol, *cands]
    for cand in cands:
        hit = lookup.get(cand.upper())
        if hit:
            return hit
    return None


def resolve_set(genes: list[str], lookup: dict[str, str]) -> list[str]:
    found = []
    seen = set()
    for gene in genes:
        hit = resolve_symbol(gene, lookup)
        if hit and hit not in seen:
            seen.add(hit)
            found.append(hit)
    return found


def run_deseq(counts_gs: pd.DataFrame, meta: pd.DataFrame, design: str, contrast) -> pd.DataFrame:
    counts = counts_gs.T
    keep = counts.columns[counts.sum(axis=0) >= 10]
    counts = counts[keep]
    meta = meta.loc[counts.index].copy()
    dds = DeseqDataSet(
        counts=counts,
        metadata=meta,
        design=design,
        refit_cooks=False,
        quiet=True,
        n_cpus=1,
        size_factors_fit_type="poscounts",
    )
    try:
        dds.deseq2()
    except Exception:
        dds = DeseqDataSet(
            counts=counts,
            metadata=meta,
            design=design,
            refit_cooks=False,
            quiet=True,
            n_cpus=1,
            size_factors_fit_type="poscounts",
            fit_type="mean",
        )
        dds.deseq2()
    stat = DeseqStats(dds, contrast=contrast, cooks_filter=False, quiet=True, n_cpus=1)
    stat.summary()
    res = stat.results_df.copy()
    res.index.name = "gene"
    return res


def paired_meta(ctrl: list[str], ko: list[str], clones: dict[str, str]) -> pd.DataFrame:
    rows = []
    for sample in ctrl:
        rows.append({"sample": sample, "condition": "ctrl", "clone": clones[sample]})
    for sample in ko:
        rows.append({"sample": sample, "condition": "ko", "clone": clones[sample]})
    meta = pd.DataFrame(rows).set_index("sample")
    meta["condition"] = pd.Categorical(meta["condition"], categories=["ctrl", "ko"])
    return meta


def unpaired_meta(ctrl: list[str], ko: list[str]) -> pd.DataFrame:
    rows = [{"sample": s, "condition": "ctrl"} for s in ctrl] + [
        {"sample": s, "condition": "ko"} for s in ko
    ]
    meta = pd.DataFrame(rows).set_index("sample")
    meta["condition"] = pd.Categorical(meta["condition"], categories=["ctrl", "ko"])
    return meta


def prerank(res: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    import gseapy as gp

    rank = res["stat"].dropna().astype(float)
    lfc = res["log2FoldChange"].reindex(rank.index).astype(float).fillna(0.0)
    rank = rank + 1e-4 * lfc
    rank = rank[~rank.index.duplicated(keep="first")].sort_values(ascending=False)
    lookup = index_lookup(rank.index)
    use = {name: resolve_set(genes, lookup) for name, genes in gene_sets.items()}
    rnk = rank.rename("score").reset_index()
    rnk.columns = ["gene", "score"]
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=use,
        outdir=None,
        min_size=8,
        max_size=500,
        permutation_num=1000,
        weight=1.0,
        seed=123,
        threads=2,
        no_plot=True,
        verbose=False,
    )
    tab = pre.res2d.copy()
    tab["n_in_set"] = tab["Term"].map(lambda t: len(use.get(t, [])))
    return tab


def gene_record(res: pd.DataFrame, symbol: str) -> dict:
    lookup = index_lookup(res.index)
    hit = resolve_symbol(symbol, lookup)
    if hit is None:
        return {
            "gene": symbol,
            "matched_symbol": "",
            "baseMean": np.nan,
            "log2FoldChange": np.nan,
            "lfcSE": np.nan,
            "stat": np.nan,
            "pvalue": np.nan,
            "padj": np.nan,
        }
    row = res.loc[hit]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    return {
        "gene": symbol,
        "matched_symbol": hit,
        "baseMean": float(row["baseMean"]) if pd.notna(row["baseMean"]) else np.nan,
        "log2FoldChange": float(row["log2FoldChange"]) if pd.notna(row["log2FoldChange"]) else np.nan,
        "lfcSE": float(row["lfcSE"]) if pd.notna(row["lfcSE"]) else np.nan,
        "stat": float(row["stat"]) if pd.notna(row["stat"]) else np.nan,
        "pvalue": float(row["pvalue"]) if pd.notna(row["pvalue"]) else np.nan,
        "padj": float(row["padj"]) if pd.notna(row["padj"]) else np.nan,
    }


def isg_table(res: pd.DataFrame, hallmark: dict[str, list[str]], contrast: dict) -> pd.DataFrame:
    lookup = index_lookup(res.index)
    membership: dict[str, set[str]] = {}
    requested: dict[str, set[str]] = {}
    for set_name, genes in hallmark.items():
        for gene in genes:
            hit = resolve_symbol(gene, lookup)
            if not hit:
                continue
            membership.setdefault(hit, set()).add(set_name)
            requested.setdefault(hit, set()).add(gene)
    rows = []
    for symbol, sets in membership.items():
        row = res.loc[symbol]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        rows.append({
            "contrast_id": contrast["contrast_id"],
            "cell_line": contrast["cell_line"],
            "accession": contrast["accession"],
            "perturbation": contrast["perturbation"],
            "condition": contrast["condition"],
            "role": contrast["role"],
            "gene": symbol,
            "requested_symbols": ";".join(sorted(requested[symbol])),
            "in_ifn_alpha": "IFN_ALPHA" in sets,
            "in_ifn_gamma": "IFN_GAMMA" in sets,
            "baseMean": float(row["baseMean"]) if pd.notna(row["baseMean"]) else np.nan,
            "log2FoldChange": float(row["log2FoldChange"]) if pd.notna(row["log2FoldChange"]) else np.nan,
            "lfcSE": float(row["lfcSE"]) if pd.notna(row["lfcSE"]) else np.nan,
            "stat": float(row["stat"]) if pd.notna(row["stat"]) else np.nan,
            "pvalue": float(row["pvalue"]) if pd.notna(row["pvalue"]) else np.nan,
            "padj": float(row["padj"]) if pd.notna(row["padj"]) else np.nan,
        })
    return pd.DataFrame(rows)


def largest_isg(tab: pd.DataFrame) -> dict:
    eligible = tab[tab["baseMean"] >= 10].dropna(subset=["log2FoldChange"])
    if eligible.empty:
        return {}
    top = eligible.sort_values("log2FoldChange", ascending=False).iloc[0]
    sig = eligible[eligible["padj"] < 0.05]
    top_sig = sig.sort_values("log2FoldChange", ascending=False).iloc[0] if len(sig) else None
    out = {
        "n_isg_baseMean_ge_10": int(len(eligible)),
        "n_isg_padj_lt_0.05_up": int((sig["log2FoldChange"] > 0).sum()) if len(sig) else 0,
        "max_gene": top["gene"],
        "max_log2FC": float(top["log2FoldChange"]),
        "max_lfcSE": float(top["lfcSE"]) if pd.notna(top["lfcSE"]) else np.nan,
        "max_baseMean": float(top["baseMean"]),
        "max_padj": float(top["padj"]) if pd.notna(top["padj"]) else np.nan,
        "max_stat": float(top["stat"]) if pd.notna(top["stat"]) else np.nan,
        "max_in_ifn_alpha": bool(top["in_ifn_alpha"]),
        "max_in_ifn_gamma": bool(top["in_ifn_gamma"]),
    }
    if top_sig is not None:
        out.update({
            "max_sig_gene": top_sig["gene"],
            "max_sig_log2FC": float(top_sig["log2FoldChange"]),
            "max_sig_padj": float(top_sig["padj"]) if pd.notna(top_sig["padj"]) else np.nan,
            "max_sig_baseMean": float(top_sig["baseMean"]),
        })
    else:
        out.update({
            "max_sig_gene": "",
            "max_sig_log2FC": np.nan,
            "max_sig_padj": np.nan,
            "max_sig_baseMean": np.nan,
        })
    return out


def gsea_records(gtab: pd.DataFrame, contrast: dict) -> list[dict]:
    colmap = {c.lower(): c for c in gtab.columns}

    def col(*names: str) -> str:
        for name in names:
            if name in gtab.columns:
                return name
            if name.lower() in colmap:
                return colmap[name.lower()]
        raise KeyError(names)

    term_c = col("Term")
    nes_c = col("NES")
    p_c = col("NOM p-val", "pval")
    fdr_c = col("FDR q-val", "fdr")
    es_c = col("ES")
    lead_c = col("Lead_genes")
    rows = []
    for _, row in gtab.iterrows():
        rows.append({
            "contrast_id": contrast["contrast_id"],
            "cell_line": contrast["cell_line"],
            "accession": contrast["accession"],
            "perturbation": contrast["perturbation"],
            "condition": contrast["condition"],
            "role": contrast["role"],
            "geneset": row[term_c],
            "ES": float(row[es_c]),
            "NES": float(row[nes_c]),
            "nom_p": float(row[p_c]),
            "fdr_q": float(row[fdr_c]),
            "n_in_set": int(row["n_in_set"]),
            "lead_genes": row[lead_c],
        })
    return rows


def build_jobs(amp: pd.DataFrame, hek: pd.DataFrame, dnapk: pd.DataFrame, mcf7: pd.DataFrame) -> list[dict]:
    jobs = []

    def add_paired(spec, counts, ko, ctrl, clones):
        spec = dict(spec)
        spec["counts"] = counts[ctrl + ko]
        spec["meta"] = paired_meta(ctrl, ko, clones)
        spec["design"] = "~C(clone) + C(condition, Treatment('ctrl'))"
        spec["contrast"] = ["condition", "ko", "ctrl"]
        spec["n_ko"] = len(ko)
        spec["n_ctrl"] = len(ctrl)
        spec["paired"] = True
        jobs.append(spec)

    def add_unpaired(spec, counts, ko, ctrl):
        spec = dict(spec)
        spec["counts"] = counts[ctrl + ko]
        spec["meta"] = unpaired_meta(ctrl, ko)
        spec["design"] = "~C(condition, Treatment('ctrl'))"
        spec["contrast"] = ["condition", "ko", "ctrl"]
        spec["n_ko"] = len(ko)
        spec["n_ctrl"] = len(ctrl)
        spec["paired"] = False
        jobs.append(spec)

    clones = {
        "Ku_1_CDKi_D4": "c1", "Ku_7_CDKi_D4": "c7", "Ku_9_CDKi_D4": "c9",
        "Ku_1_CDKiIAA_D4": "c1", "Ku_7_CDKiIAA_D4": "c7", "Ku_9_CDKiIAA_D4": "c9",
    }
    add_paired(
        dict(
            contrast_id="HCT116_Ku80AID_CDKi_D4",
            cell_line="HCT116",
            accession="GSE294709",
            perturbation="XRCC5 (Ku80-AID)",
            condition="IAA+Dox+CDKi vs CDKi, day 4, clones 1/7/9",
            role="primary",
            note="Proliferation-matched Ku80 degradation. Same three AID clones.",
        ),
        amp,
        ["Ku_1_CDKiIAA_D4", "Ku_7_CDKiIAA_D4", "Ku_9_CDKiIAA_D4"],
        ["Ku_1_CDKi_D4", "Ku_7_CDKi_D4", "Ku_9_CDKi_D4"],
        clones,
    )
    clones = {
        "Ku_1_DMSO": "c1", "Ku_7_DMSO": "c7", "Ku_9_DMSO": "c9",
        "Ku_1_IAA_D2": "c1", "Ku_7_IAA_D2": "c7", "Ku_9_IAA_D2": "c9",
    }
    add_paired(
        dict(
            contrast_id="HCT116_Ku80AID_IAA_D2",
            cell_line="HCT116",
            accession="GSE294709",
            perturbation="XRCC5 (Ku80-AID)",
            condition="IAA+Dox day 2 vs DMSO, clones 1/7/9",
            role="timecourse",
            note="Acute degron without the CDK4/6 inhibitor. Same AID system, not an independent line.",
        ),
        amp,
        ["Ku_1_IAA_D2", "Ku_7_IAA_D2", "Ku_9_IAA_D2"],
        ["Ku_1_DMSO", "Ku_7_DMSO", "Ku_9_DMSO"],
        clones,
    )
    clones = {
        "Ku_1_DMSO_24h": "c1", "Ku_9_DMSO_24h": "c9",
        "Ku_1_IAA_24h": "c1", "Ku_9_IAA_24h": "c9",
    }
    add_paired(
        dict(
            contrast_id="HCT116_Ku80AID_IAA_24h",
            cell_line="HCT116",
            accession="GSE294709",
            perturbation="XRCC5 (Ku80-AID)",
            condition="IAA 24h vs DMSO 24h, clones 1 and 9",
            role="timecourse",
            note="Clone 7 has no 24h pair in the deposited matrix.",
        ),
        amp,
        ["Ku_1_IAA_24h", "Ku_9_IAA_24h"],
        ["Ku_1_DMSO_24h", "Ku_9_DMSO_24h"],
        clones,
    )
    clones = {
        "Ku_1_CDKi_48h": "c1", "Ku_9_CDKi_48h": "c9",
        "Ku_1_CDKiIAA_48h": "c1", "Ku_9_CDKiIAA_48h": "c9",
    }
    add_paired(
        dict(
            contrast_id="HCT116_Ku80AID_CDKi_48h",
            cell_line="HCT116",
            accession="GSE294709",
            perturbation="XRCC5 (Ku80-AID)",
            condition="CDKi+IAA 48h vs CDKi 48h, clones 1 and 9",
            role="timecourse",
            note="Shorter proliferation-matched degron. Clone 7 is absent.",
        ),
        amp,
        ["Ku_1_CDKiIAA_48h", "Ku_9_CDKiIAA_48h"],
        ["Ku_1_CDKi_48h", "Ku_9_CDKi_48h"],
        clones,
    )
    clones = {
        "MAVS1_DMSO": "c1", "MAVS13_DMSO": "c13", "MAVS15_DMSO": "c15",
        "MAVS1_IAA": "c1", "MAVS13_IAA": "c13", "MAVS15_IAA": "c15",
    }
    add_paired(
        dict(
            contrast_id="HCT116_Ku80AID_MAVSKO",
            cell_line="HCT116",
            accession="GSE294709",
            perturbation="XRCC5 (Ku80-AID), MAVS knockout",
            condition="IAA vs DMSO on MAVS-KO Ku80-AID clones",
            role="epistasis",
            note="Ku80 still degraded. Not a wild-type MAVS contrast.",
        ),
        amp,
        ["MAVS1_IAA", "MAVS13_IAA", "MAVS15_IAA"],
        ["MAVS1_DMSO", "MAVS13_DMSO", "MAVS15_DMSO"],
        clones,
    )
    add_unpaired(
        dict(
            contrast_id="HCT116_Ku80flox_Cre_4OHT",
            cell_line="HCT116",
            accession="GSE294709",
            perturbation="XRCC5 (Ku86flox CreERT2)",
            condition="4OHT day 4 vs EtOH, CreERT2",
            role="unadjusted",
            note="Excision contrast before subtracting the no-Cre 4OHT arm.",
        ),
        amp,
        ["Cre_D4", "Cre_D4_2"],
        ["Cre_D0", "Cre_D0_2"],
    )
    add_unpaired(
        dict(
            contrast_id="HCT116_Ku80flox_noCre_4OHT",
            cell_line="HCT116",
            accession="GSE294709",
            perturbation="none (4OHT, no Cre)",
            condition="4OHT day 4 vs EtOH, no Cre",
            role="negative_control",
            note="Same 4OHT exposure without Ku80 excision.",
        ),
        amp,
        ["Ctrl_D4", "Ctrl_D4_2"],
        ["Ctrl_D0", "Ctrl_D0_2"],
    )

    # Interaction is a custom design.
    cols = ["Cre_D0", "Cre_D0_2", "Cre_D4", "Cre_D4_2", "Ctrl_D0", "Ctrl_D0_2", "Ctrl_D4", "Ctrl_D4_2"]
    meta = pd.DataFrame({
        "genotype": ["Cre", "Cre", "Cre", "Cre", "noCre", "noCre", "noCre", "noCre"],
        "treatment": ["EtOH", "EtOH", "fourOHT", "fourOHT", "EtOH", "EtOH", "fourOHT", "fourOHT"],
    }, index=cols)
    meta["genotype"] = pd.Categorical(meta["genotype"], categories=["noCre", "Cre"])
    meta["treatment"] = pd.Categorical(meta["treatment"], categories=["EtOH", "fourOHT"])
    jobs.append(dict(
        contrast_id="HCT116_Ku80flox_interaction",
        cell_line="HCT116",
        accession="GSE294709",
        perturbation="XRCC5 (Ku86flox CreERT2 minus no-Cre)",
        condition="4OHT effect in Cre beyond 4OHT without Cre",
        role="primary",
        note="Interaction coefficient. This is the Ku-excision effect after the 4OHT-only arm is subtracted.",
        counts=amp[cols],
        meta=meta,
        design="~ C(genotype, Treatment('noCre')) * C(treatment, Treatment('EtOH'))",
        contrast="interaction",
        n_ko=2,
        n_ctrl=2,
        paired=False,
    ))

    clones = {
        "SB_Dox": "SB", "Sa11_Dox": "Sa11", "TI_Dox": "TI",
        "SB_D8": "SB", "Sa11_D8": "Sa11", "TI_D8": "TI",
    }
    add_paired(
        dict(
            contrast_id="HEK293_Ku70_dox_withdrawal",
            cell_line="HEK293",
            accession="GSE294709",
            perturbation="XRCC6 (Ku70, dox-off)",
            condition="8 days without dox vs dox rescue, clones Sa11/SB/TI",
            role="other_cell_line",
            note="Same GEO series, different cell line and Ku subunit. Not an HCT116 result.",
        ),
        hek,
        ["SB_D8", "Sa11_D8", "TI_D8"],
        ["SB_Dox", "Sa11_Dox", "TI_Dox"],
        clones,
    )

    add_unpaired(
        dict(
            contrast_id="HCT116_PRKDC_KO_normoxia",
            cell_line="HCT116",
            accession="GSE285698",
            perturbation="PRKDC (DNA-PKcs KO)",
            condition="normoxia, H2O vehicle, KO vs WT, n=3",
            role="primary",
            note="Unstressed DNA-PKcs knockout. Library ids A1/B1/C1 WT and A2/B2/C2 KO.",
        ),
        dnapk,
        ["KO_normoxia_rep1", "KO_normoxia_rep2", "KO_normoxia_rep3"],
        ["WT_normoxia_rep1", "WT_normoxia_rep2", "WT_normoxia_rep3"],
    )
    add_unpaired(
        dict(
            contrast_id="HCT116_PRKDC_KO_hypoxia",
            cell_line="HCT116",
            accession="GSE285698",
            perturbation="PRKDC (DNA-PKcs KO)",
            condition="CoCl2 200 uM 12h, KO vs WT, n=3",
            role="secondary",
            note="Same knockout under CoCl2. Not the unstressed contrast.",
        ),
        dnapk,
        ["KO_hypoxia_rep1", "KO_hypoxia_rep2", "KO_hypoxia_rep3"],
        ["WT_hypoxia_rep1", "WT_hypoxia_rep2", "WT_hypoxia_rep3"],
    )

    wt_u = ["WT_untreated_1", "WT_untreated_2", "WT_untreated_3"]
    ko_u = [
        "KO1_untreated_1", "KO1_untreated_2", "KO1_untreated_3",
        "KO2_untreated_1", "KO2_untreated_2", "KO2_untreated_3",
    ]
    add_unpaired(
        dict(
            contrast_id="MCF7_53BP1_KO_untreated",
            cell_line="MCF-7",
            accession="GSE84986",
            perturbation="TP53BP1 (two null clones)",
            condition="untreated, clones 1 and 2 pooled vs WT",
            role="primary",
            note="Six KO samples are two clones of one parental line, not six independent knockouts.",
        ),
        mcf7,
        ko_u,
        wt_u,
    )
    add_unpaired(
        dict(
            contrast_id="MCF7_53BP1_KO_clone1_untreated",
            cell_line="MCF-7",
            accession="GSE84986",
            perturbation="TP53BP1 clone 1",
            condition="untreated clone 1 vs WT, n=3",
            role="clone",
            note="Clone-level contrast.",
        ),
        mcf7,
        ["KO1_untreated_1", "KO1_untreated_2", "KO1_untreated_3"],
        wt_u,
    )
    add_unpaired(
        dict(
            contrast_id="MCF7_53BP1_KO_clone2_untreated",
            cell_line="MCF-7",
            accession="GSE84986",
            perturbation="TP53BP1 clone 2",
            condition="untreated clone 2 vs WT, n=3",
            role="clone",
            note="Clone-level contrast. Shares the WT arm with clone 1.",
        ),
        mcf7,
        ["KO2_untreated_1", "KO2_untreated_2", "KO2_untreated_3"],
        wt_u,
    )
    add_unpaired(
        dict(
            contrast_id="MCF7_53BP1_KO_IR",
            cell_line="MCF-7",
            accession="GSE84986",
            perturbation="TP53BP1 (two null clones)",
            condition="5 Gy, 4 h, clones 1 and 2 pooled vs WT",
            role="secondary",
            note="Irradiated arm. Not the baseline contrast.",
        ),
        mcf7,
        [
            "KO1_IR_1", "KO1_IR_2", "KO1_IR_3",
            "KO2_IR_1", "KO2_IR_2", "KO2_IR_3",
        ],
        ["WT_IR_1", "WT_IR_2", "WT_IR_3"],
    )
    return jobs


def run_job(job: dict) -> pd.DataFrame:
    if job["contrast"] == "interaction":
        counts = job["counts"].T
        keep = counts.columns[counts.sum(axis=0) >= 10]
        counts = counts[keep]
        meta = job["meta"].loc[counts.index].copy()
        dds = DeseqDataSet(
            counts=counts,
            metadata=meta,
            design=job["design"],
            refit_cooks=False,
            quiet=True,
            n_cpus=1,
            size_factors_fit_type="poscounts",
        )
        dds.deseq2()
        colnames = list(dds.obsm["design_matrix"].columns)
        inter = [c for c in colnames if ":" in str(c)]
        if len(inter) != 1:
            raise RuntimeError(f"expected one interaction column, got {colnames}")
        contrast = np.zeros(len(colnames))
        contrast[colnames.index(inter[0])] = 1.0
        stat = DeseqStats(dds, contrast=contrast, cooks_filter=False, quiet=True, n_cpus=1)
        stat.summary()
        res = stat.results_df.copy()
        res.index.name = "gene"
        return res
    return run_deseq(job["counts"], job["meta"], job["design"], job["contrast"])


def plot(maxima: pd.DataFrame) -> None:
    show = maxima.copy()
    show["label"] = show["cell_line"] + " | " + show["contrast_id"]
    show = show.iloc[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 7.4), sharey=True)
    colors = {"HCT116": "#1f4e79", "MCF-7": "#b85c38", "HEK293": "#5c6b73"}
    c = [colors.get(x, "#333") for x in show["cell_line"]]
    axes[0].barh(show["label"], show["max_sig_log2FC"], color=c)
    for y, (_, row) in enumerate(show.iterrows()):
        if pd.notna(row["max_sig_log2FC"]):
            axes[0].text(
                row["max_sig_log2FC"] + 0.08,
                y,
                str(row["max_sig_gene"]),
                va="center",
                fontsize=7,
                color="#222",
            )
    axes[0].axvline(0, color="black", lw=0.6)
    axes[0].set_xlabel("Largest Hallmark ISG log2FC with padj < 0.05\n(baseMean ≥ 10)")
    axes[1].barh(show["label"], show["max_NES"], color=c)
    for y, (_, row) in enumerate(show.iterrows()):
        if pd.notna(row["max_NES"]):
            axes[1].text(
                row["max_NES"] + (0.06 if row["max_NES"] >= 0 else -0.06),
                y,
                str(row["max_NES_set"]).replace("IFN_", ""),
                va="center",
                ha="left" if row["max_NES"] >= 0 else "right",
                fontsize=7,
                color="#222",
            )
    axes[1].axvline(0, color="black", lw=0.6)
    axes[1].set_xlabel("Larger of Hallmark IFN-α and IFN-γ NES")
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("Each bar is one cell line and one contrast.", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "max_isg_log2fc_nes.png", dpi=160)
    fig.savefig(FIGS / "max_isg_log2fc_nes.pdf")
    plt.close(fig)


def write_top_isg(isg_all: pd.DataFrame, n: int = 10) -> None:
    parts = []
    for cid, sub in isg_all.groupby("contrast_id", sort=False):
        sig = sub[(sub["baseMean"] >= 10) & (sub["padj"] < 0.05)].dropna(subset=["log2FoldChange"])
        sig = sig.sort_values("log2FoldChange", ascending=False).head(n)
        sig = sig.copy()
        sig.insert(0, "rank", range(1, len(sig) + 1))
        parts.append(sig)
    top = pd.concat(parts, ignore_index=True)
    cols = [
        "contrast_id", "cell_line", "accession", "role", "rank", "gene",
        "log2FoldChange", "lfcSE", "baseMean", "stat", "padj",
        "in_ifn_alpha", "in_ifn_gamma",
    ]
    top[cols].to_csv(TABLES / "top_isg_padj05.tsv", sep="\t", index=False)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    ensure_inputs()
    print("loading matrices", flush=True)
    amp = load_amp()
    hek = load_hek()
    dnapk = load_dnapk()
    ens_map = load_hgnc()
    mcf7 = load_mcf7(ens_map)
    print({
        "amp": amp.shape, "hek": hek.shape, "dnapk": dnapk.shape, "mcf7": mcf7.shape,
    }, flush=True)
    hallmark = load_hallmark()
    (OUT / "genesets.json").write_text(json.dumps(hallmark, indent=1))
    print({k: len(v) for k, v in hallmark.items()}, flush=True)

    jobs = build_jobs(amp, hek, dnapk, mcf7)
    isg_parts = []
    gsea_rows = []
    maxima_rows = []
    target_rows = []
    targets = ["XRCC5", "XRCC6", "PRKDC", "TP53BP1", "MAVS", "ISG15", "IFIT1", "IFI44L", "RSAD2", "CXCL10"]

    for job in jobs:
        print(f"\n=== {job['contrast_id']} ===", flush=True)
        res = run_job(job)
        isg = isg_table(res, hallmark, job)
        isg_parts.append(isg)
        top = largest_isg(isg)
        print("  GSEA", flush=True)
        gtab = prerank(res, hallmark)
        grows = gsea_records(gtab, job)
        gsea_rows.extend(grows)
        nes_map = {r["geneset"]: r for r in grows}
        alpha = nes_map.get("IFN_ALPHA", {})
        gamma = nes_map.get("IFN_GAMMA", {})
        candidates = [r for r in (alpha, gamma) if r]
        best = max(candidates, key=lambda r: r["NES"]) if candidates else {}
        for symbol in targets:
            rec = gene_record(res, symbol)
            rec.update({
                "contrast_id": job["contrast_id"],
                "cell_line": job["cell_line"],
                "accession": job["accession"],
            })
            target_rows.append(rec)
        xrcc5 = gene_record(res, "XRCC5")
        prkdc = gene_record(res, "PRKDC")
        tp53bp1 = gene_record(res, "TP53BP1")
        xrcc6 = gene_record(res, "XRCC6")
        row = {
            "contrast_id": job["contrast_id"],
            "cell_line": job["cell_line"],
            "accession": job["accession"],
            "perturbation": job["perturbation"],
            "condition": job["condition"],
            "role": job["role"],
            "n_ko": job["n_ko"],
            "n_ctrl": job["n_ctrl"],
            "paired": job["paired"],
            "note": job["note"],
            "XRCC5_log2FC": xrcc5["log2FoldChange"],
            "XRCC5_padj": xrcc5["padj"],
            "XRCC6_log2FC": xrcc6["log2FoldChange"],
            "PRKDC_log2FC": prkdc["log2FoldChange"],
            "PRKDC_padj": prkdc["padj"],
            "TP53BP1_log2FC": tp53bp1["log2FoldChange"],
            "TP53BP1_padj": tp53bp1["padj"],
            "IFN_ALPHA_NES": alpha.get("NES", np.nan),
            "IFN_ALPHA_FDR": alpha.get("fdr_q", np.nan),
            "IFN_ALPHA_nom_p": alpha.get("nom_p", np.nan),
            "IFN_GAMMA_NES": gamma.get("NES", np.nan),
            "IFN_GAMMA_FDR": gamma.get("fdr_q", np.nan),
            "IFN_GAMMA_nom_p": gamma.get("nom_p", np.nan),
            "max_NES_set": best.get("geneset", ""),
            "max_NES": best.get("NES", np.nan),
            "max_NES_FDR": best.get("fdr_q", np.nan),
            **top,
        }
        maxima_rows.append(row)
        print(
            f"  max ISG {top.get('max_gene')} {top.get('max_log2FC')}  "
            f"NES {best.get('geneset')} {best.get('NES')}",
            flush=True,
        )

    isg_all = pd.concat(isg_parts, ignore_index=True)
    maxima = pd.DataFrame(maxima_rows)
    gsea = pd.DataFrame(gsea_rows)
    targets_df = pd.DataFrame(target_rows)
    isg_all.to_csv(TABLES / "isg_log2fc.tsv.gz", sep="\t", index=False, compression="gzip")
    maxima.to_csv(TABLES / "contrast_maxima.tsv", sep="\t", index=False)
    gsea.to_csv(TABLES / "gsea_nes.tsv", sep="\t", index=False)
    targets_df.to_csv(TABLES / "target_and_marker_genes.tsv", sep="\t", index=False)
    write_top_isg(isg_all)
    plot(maxima)
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
