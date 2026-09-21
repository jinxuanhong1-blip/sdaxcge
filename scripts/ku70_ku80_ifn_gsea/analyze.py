#!/usr/bin/env python3
"""Ku70/Ku80 (XRCC6/XRCC5) KO/KD RNA-seq meta-analysis: IFN, STING, APM.

Pre-ranked GSEA on DESeq2 Wald statistics. Every accession below was returned
by NCBI GEO or the ENCODE portal (queried 2026-09-21). Nothing here is a
private matrix.

Primary question
----------------
When Ku70 (XRCC6) or Ku80 (XRCC5) is knocked out or knocked down, do
interferon, STING, and antigen-presentation (APM) gene sets move, and in
which direction?

Ranking metric: DESeq2 Wald statistic (KO/KD vs matched control). Positive
NES = the set is enriched among genes up after Ku loss.

Ku/DNA-PKcs genes (XRCC5, XRCC6, PRKDC) are removed from every gene set
before testing so an on-target drop cannot manufacture a negative NES.
"""

from __future__ import annotations

import json
import math
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import norm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parents[2]
DATA = Path("/tmp/ku_data")
OUT = ROOT / "results" / "ku70_ku80_ifn_gsea"
GENESET_DIR = DATA / "genesets"

# Perturbed subunits. Stripped from every tested set.
PERTURBED = {"XRCC5", "XRCC6", "PRKDC"}

# Symbol aliases used by individual matrices (GSE294709 uses RIGI, not DDX58).
ALIASES = {
    "DDX58": ["RIGI", "DDX58"],
    "IFIH1": ["IFIH1", "MDA5"],
    "CGAS": ["CGAS", "MB21D1"],
    "STING1": ["STING1", "TMEM173"],
    "MRE11": ["MRE11", "MRE11A"],
}

FOCUS = [
    "XRCC5", "XRCC6", "PRKDC", "LIG4",
    "ISG15", "IFIT1", "IFIT2", "MX1", "OAS1", "OAS2", "RSAD2", "IFI27",
    "STAT1", "STAT2", "IRF1", "IRF3", "IRF7", "IFNB1", "IFNL1",
    "CGAS", "STING1", "TBK1", "DDX58", "IFIH1", "MAVS",
    "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5",
    "HLA-A", "HLA-B", "HLA-C",
]

MOUSE_MHC_CLASSICAL = ["H2-K1", "H2-D1", "H2-L", "H2-T23"]


def _read_enrichr(path: Path) -> dict[str, list[str]]:
    out = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            name = parts[0]
            genes = [g for g in parts[2:] if g]
            out[name] = genes
    return out


def load_human_sets() -> dict[str, list[str]]:
    hallmark = _read_enrichr(GENESET_DIR / "hallmark.txt")
    reactome = _read_enrichr(GENESET_DIR / "reactome.txt")
    kegg = _read_enrichr(GENESET_DIR / "kegg.txt")
    # Classical MHC-I antigen-presentation machinery. The KEGG set also contains
    # MHC-II, CD4/CD8 and heat-shock genes, which can enrich without TAP/HLA.
    apm_mhci = [
        "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
        "B2M", "TAP1", "TAP2", "TAPBP", "TAPBPL",
        "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2",
        "NLRC5", "ERAP1", "ERAP2",
    ]
    sets = {
        "IFN_ALPHA": hallmark["Interferon Alpha Response"],
        "IFN_GAMMA": hallmark["Interferon Gamma Response"],
        "STING": reactome[
            "STING Mediated Induction Of Host Immune Responses R-HSA-1834941"
        ],
        "APM_MHCI": apm_mhci,
        "APM_KEGG": kegg["Antigen processing and presentation"],
        "RIGI_MDA5": reactome[
            "DDX58/IFIH1-mediated Induction Of Interferon-Alpha/Beta R-HSA-168928"
        ],
    }
    cleaned = {}
    for name, genes in sets.items():
        keep = []
        seen = set()
        for g in genes:
            g = g.strip().upper()
            if not g or g in PERTURBED or g in seen:
                continue
            seen.add(g)
            keep.append(g)
        cleaned[name] = keep
    return cleaned


def fetch_mouse_orthologs(symbols: list[str]) -> dict[str, list[str]]:
    """Human symbol -> mouse symbols.

    Uses mygene homologene. A human gene with exactly one mouse homolog
    contributes that symbol. HLA genes (many H2 paralogs) are not expanded;
    classical mouse MHC-I genes are added once to the APM set by the caller.
    """
    cache = OUT / "mouse_orthologs.json"
    if cache.exists():
        return json.loads(cache.read_text())
    mapping: dict[str, list[str]] = {}
    chunk = 40
    for i in range(0, len(symbols), chunk):
        batch = symbols[i : i + chunk]
        resp = requests.post(
            "https://mygene.info/v3/query",
            data={
                "q": ",".join(batch),
                "scopes": "symbol",
                "fields": "symbol,homologene",
                "species": "human",
                "size": 1,
            },
            timeout=90,
        )
        resp.raise_for_status()
        hits = resp.json()
        entrez_for: dict[str, list[str]] = {}
        all_entrez = []
        for hit in hits:
            if not isinstance(hit, dict) or hit.get("notfound"):
                continue
            q = str(hit.get("query"))
            genes = (hit.get("homologene") or {}).get("genes") or []
            mice = []
            for g in genes:
                if isinstance(g, list) and len(g) >= 2 and g[0] == 10090:
                    mice.append(str(g[1]))
            entrez_for[q] = mice
            all_entrez.extend(mice)
        id_to_sym = {}
        if all_entrez:
            r2 = requests.post(
                "https://mygene.info/v3/gene",
                data={
                    "ids": ",".join(sorted(set(all_entrez))),
                    "fields": "symbol",
                    "species": "mouse",
                },
                timeout=90,
            )
            r2.raise_for_status()
            for hit in r2.json():
                if isinstance(hit, dict) and hit.get("symbol"):
                    id_to_sym[str(hit.get("query"))] = hit["symbol"]
        for q, ids in entrez_for.items():
            syms = [id_to_sym[i] for i in ids if i in id_to_sym]
            # de-duplicate, preserve order
            dedup = list(dict.fromkeys(syms))
            if len(dedup) == 1:
                mapping[q.upper()] = dedup
            else:
                simple = q[0].upper() + q[1:].lower()
                if simple in dedup:
                    mapping[q.upper()] = [simple]
                else:
                    mapping[q.upper()] = []
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(mapping, indent=1))
    return mapping


def to_mouse_sets(human_sets: dict[str, list[str]], orth: dict[str, list[str]]) -> dict[str, list[str]]:
    out = {}
    for name, genes in human_sets.items():
        mice = []
        seen = set()
        for g in genes:
            for m in orth.get(g, []):
                if m not in seen:
                    seen.add(m)
                    mice.append(m)
            # title-case fallback is applied later against the matrix
            simple = g[0].upper() + g[1:].lower()
            if simple not in seen:
                # keep as a candidate; resolve_against_index will drop misses
                seen.add(simple)
                mice.append(simple)
        if name in ("APM_KEGG", "APM_MHCI"):
            for m in MOUSE_MHC_CLASSICAL:
                if m not in seen:
                    mice.append(m)
        out[name] = mice
    return out


def resolve_genes(genes: list[str], index: pd.Index) -> list[str]:
    """Map requested symbols onto the matrix index, applying aliases."""
    lookup = {str(i).upper(): str(i) for i in index}
    found = []
    seen = set()
    for g in genes:
        cands = ALIASES.get(g.upper(), [g])
        # also try the raw symbol
        if g not in cands:
            cands = [g, *cands]
        hit = None
        for c in cands:
            if c.upper() in lookup:
                hit = lookup[c.upper()]
                break
        if hit and hit not in seen:
            seen.add(hit)
            found.append(hit)
    return found


def collapse_symbols(df: pd.DataFrame, symbol_col: str) -> pd.DataFrame:
    """Genes x samples, integer counts, duplicate symbols summed."""
    d = df.copy()
    d[symbol_col] = d[symbol_col].astype(str)
    d = d[d[symbol_col].notna() & (d[symbol_col] != "") & (d[symbol_col] != "nan")]
    num = d.drop(columns=[symbol_col]).apply(pd.to_numeric, errors="coerce").fillna(0)
    num = num.clip(lower=0).round().astype(int)
    num.index = d[symbol_col].values
    num = num.groupby(level=0).sum()
    num = num.loc[num.sum(axis=1) > 0]
    return num


def load_ampseq() -> pd.DataFrame:
    raw = pd.read_csv(DATA / "gse294709" / "counts.txt.gz", sep="\t")
    return collapse_symbols(raw, "gene_names")


def load_kallisto() -> pd.DataFrame:
    raw = pd.read_csv(DATA / "gse294709" / "kallisto_counts.txt.gz", sep="\t")
    return collapse_symbols(raw, "gene_name")


def load_hek() -> pd.DataFrame:
    raw = pd.read_excel(DATA / "gse180581" / "counts.xlsx")
    hgnc = pd.read_csv(
        DATA / "hgnc_complete_set.txt",
        sep="\t",
        usecols=["symbol", "ensembl_gene_id"],
        dtype=str,
    )
    hgnc = hgnc.dropna(subset=["ensembl_gene_id", "symbol"])
    hgnc = hgnc.drop_duplicates("ensembl_gene_id")
    ens_to_sym = dict(zip(hgnc["ensembl_gene_id"], hgnc["symbol"]))
    raw["symbol"] = raw["name"].map(ens_to_sym)
    raw = raw.dropna(subset=["symbol"])
    cols = [c for c in raw.columns if c not in ("name", "symbol")]
    return collapse_symbols(raw[["symbol", *cols]], "symbol")


def load_treg() -> pd.DataFrame:
    raw = pd.read_csv(DATA / "gse247031" / "counts.csv.gz")
    raw.columns = [c.replace(" Read Count", "").replace("\ufeff", "") for c in raw.columns]
    # first column is gene symbol
    sym_col = raw.columns[0]
    raw = raw.rename(columns={sym_col: "symbol"})
    raw["symbol"] = raw["symbol"].astype(str).str.strip("'").str.strip()
    return collapse_symbols(raw, "symbol")


def load_hgnc_maps() -> tuple[dict[str, str], dict[str, str]]:
    """Ensembl gene id (no version) and Entrez id -> HGNC symbol.

    ENCODE GRCh38 V29 gene quantifications use versioned Ensembl ids, with a
    small number of numeric ids and ERCC spike-ins mixed in.
    """
    hgnc = pd.read_csv(
        DATA / "hgnc_complete_set.txt",
        sep="\t",
        usecols=["symbol", "entrez_id", "ensembl_gene_id"],
        dtype=str,
    )
    ens = hgnc.dropna(subset=["ensembl_gene_id", "symbol"]).drop_duplicates("ensembl_gene_id")
    ent = hgnc.dropna(subset=["entrez_id", "symbol"]).drop_duplicates("entrez_id")
    return (
        dict(zip(ens["ensembl_gene_id"], ens["symbol"])),
        dict(zip(ent["entrez_id"], ent["symbol"])),
    )


def load_encode_pair(ko_prefix: str, ctrl_prefix: str, ens_map: dict[str, str], ent_map: dict[str, str]) -> pd.DataFrame:
    frames = []
    for label, prefix in [("ko1", f"{ko_prefix}_1"), ("ko2", f"{ko_prefix}_2"),
                          ("ctrl1", f"{ctrl_prefix}_1"), ("ctrl2", f"{ctrl_prefix}_2")]:
        path = DATA / "encode_counts" / f"{prefix}.tsv"
        df = pd.read_csv(path, sep="\t", usecols=["gene_id", "expected_count"])
        gid = df["gene_id"].astype(str)
        base = gid.str.replace(r"\.\d+$", "", regex=True)
        sym = base.map(ens_map)
        numeric = gid.str.fullmatch(r"\d+")
        sym = sym.where(~numeric, gid.map(ent_map))
        df = df.assign(symbol=sym).dropna(subset=["symbol"])
        s = df.groupby("symbol")["expected_count"].sum()
        s.name = label
        frames.append(s)
    mat = pd.concat(frames, axis=1).fillna(0)
    mat = mat.clip(lower=0).round().astype(int)
    mat = mat.loc[mat.sum(axis=1) > 0]
    if mat.shape[0] < 1000:
        raise RuntimeError(f"{ko_prefix}: only {mat.shape[0]} genes mapped; ID mapping failed")
    return mat


def run_deseq(counts_gs: pd.DataFrame, meta: pd.DataFrame, paired: bool) -> pd.DataFrame:
    """counts_gs: genes x samples. Returns gene-indexed DE table, KO vs ctrl."""
    counts = counts_gs.T  # samples x genes
    # drop genes with < 10 total counts
    keep = counts.columns[counts.sum(axis=0) >= 10]
    counts = counts[keep]
    meta = meta.loc[counts.index].copy()
    meta["condition"] = pd.Categorical(meta["condition"], categories=["ctrl", "ko"])
    if paired:
        design = "~C(clone) + C(condition, Treatment('ctrl'))"
    else:
        design = "~C(condition, Treatment('ctrl'))"
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
    stat = DeseqStats(
        dds,
        contrast=["condition", "ko", "ctrl"],
        cooks_filter=False,
        quiet=True,
        n_cpus=1,
    )
    stat.summary()
    res = stat.results_df.copy()
    res.index.name = "gene"
    return res


def run_flox_interaction(amp: pd.DataFrame) -> pd.DataFrame:
    """Ku-specific effect of 4OHT: (Cre 4OHT − Cre EtOH) − (noCre 4OHT − noCre EtOH)."""
    cols = ["Cre_D0", "Cre_D0_2", "Cre_D4", "Cre_D4_2", "Ctrl_D0", "Ctrl_D0_2", "Ctrl_D4", "Ctrl_D4_2"]
    counts = amp[cols].T
    keep = counts.columns[counts.sum(axis=0) >= 10]
    counts = counts[keep]
    meta = pd.DataFrame({
        "genotype": ["Cre", "Cre", "Cre", "Cre", "noCre", "noCre", "noCre", "noCre"],
        "treatment": ["EtOH", "EtOH", "fourOHT", "fourOHT", "EtOH", "EtOH", "fourOHT", "fourOHT"],
    }, index=cols)
    meta["genotype"] = pd.Categorical(meta["genotype"], categories=["noCre", "Cre"])
    meta["treatment"] = pd.Categorical(meta["treatment"], categories=["EtOH", "fourOHT"])
    design = "~ C(genotype, Treatment('noCre')) * C(treatment, Treatment('EtOH'))"
    dds = DeseqDataSet(
        counts=counts,
        metadata=meta,
        design=design,
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


def prerank(res: pd.DataFrame, gene_sets: dict[str, list[str]], organism_index: pd.Index) -> pd.DataFrame:
    import gseapy as gp

    rank = res["stat"].dropna().astype(float)
    # Break exact Wald-stat ties with a negligible log2FC term so gene order is
    # deterministic. The tie-breaker is far smaller than any distinct stat.
    lfc = res["log2FoldChange"].reindex(rank.index).astype(float).fillna(0.0)
    rank = rank + 1e-4 * lfc
    rank = rank[~rank.index.duplicated(keep="first")]
    rank = rank.sort_values(ascending=False)
    use_sets = {}
    for name, genes in gene_sets.items():
        resolved = resolve_genes(genes, rank.index)
        # also allow genes already in matrix symbols (mouse)
        if len(resolved) < 10:
            # try direct membership for mouse symbols already in `genes`
            extra = [g for g in genes if g in rank.index and g not in resolved]
            resolved = list(dict.fromkeys(resolved + extra))
        use_sets[name] = resolved
    # gseapy wants a DataFrame with gene and score, or a Series
    rnk = rank.rename("score").reset_index()
    rnk.columns = ["gene", "score"]
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=use_sets,
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
    tab["n_in_set_resolved"] = tab["Term"].map(lambda t: len(use_sets.get(t, [])))
    return tab


def stouffer(p_up: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    p = np.clip(np.asarray(p_up, dtype=float), 1e-12, 1 - 1e-12)
    w = np.asarray(weights, dtype=float)
    z = norm.isf(p)
    Z = float(np.sum(w * z) / math.sqrt(np.sum(w ** 2)))
    return Z, float(norm.sf(Z))


def gene_row(res: pd.DataFrame, gene: str) -> dict:
    resolved = resolve_genes([gene], res.index)
    if not resolved:
        # mouse title-case
        simple = gene[0].upper() + gene[1:].lower()
        if simple in res.index:
            resolved = [simple]
        elif gene in res.index:
            resolved = [gene]
    if not resolved:
        return {"gene": gene, "log2FoldChange": np.nan, "stat": np.nan, "pvalue": np.nan, "padj": np.nan}
    g = resolved[0]
    r = res.loc[g]
    if isinstance(r, pd.DataFrame):
        r = r.iloc[0]
    return {
        "gene": gene,
        "matched_symbol": g,
        "log2FoldChange": float(r["log2FoldChange"]),
        "stat": float(r["stat"]) if pd.notna(r["stat"]) else np.nan,
        "pvalue": float(r["pvalue"]) if pd.notna(r["pvalue"]) else np.nan,
        "padj": float(r["padj"]) if pd.notna(r["padj"]) else np.nan,
    }


def build_contrasts(amp, kal, hek, treg, ens_map, ent_map):
    """Return list of dicts with counts (genes x samples) and metadata."""
    specs = []

    def add(spec, counts, samples_ko, samples_ctrl, clones_ko=None, clones_ctrl=None):
        cols = list(samples_ctrl) + list(samples_ko)
        missing = [c for c in cols if c not in counts.columns]
        if missing:
            raise KeyError(f"{spec['contrast_id']} missing columns {missing}")
        sub = counts[cols]
        meta_rows = []
        for s in samples_ctrl:
            meta_rows.append({"sample": s, "condition": "ctrl", "clone": (clones_ctrl or {}).get(s, "c")})
        for s in samples_ko:
            meta_rows.append({"sample": s, "condition": "ko", "clone": (clones_ko or {}).get(s, "c")})
        meta = pd.DataFrame(meta_rows).set_index("sample")
        spec = dict(spec)
        spec["counts"] = sub
        spec["meta"] = meta
        spec["n_ko"] = len(samples_ko)
        spec["n_ctrl"] = len(samples_ctrl)
        specs.append(spec)
        return spec

    # --- GSE294709 HCT116 Ku80-AID, proliferation-matched day 4 (PRIMARY) ---
    clones = {
        "Ku_1_CDKi_D4": "c1", "Ku_7_CDKi_D4": "c7", "Ku_9_CDKi_D4": "c9",
        "Ku_1_CDKiIAA_D4": "c1", "Ku_7_CDKiIAA_D4": "c7", "Ku_9_CDKiIAA_D4": "c9",
    }
    add(
        dict(
            contrast_id="GSE294709_HCT116_Ku80AID_CDKi_D4",
            accession="GSE294709",
            organism="human",
            system="HCT116 Ku80-AID",
            perturbation="XRCC5",
            modality="AID protein degron + CDK4/6 inhibitor",
            tier="primary",
            role="Ku loss",
            note="IAA+Dox+CDKi vs CDKi alone, day 4, clones 1/7/9. Proliferation-matched Ku80 degradation (Zhu et al., Nature 2025).",
        ),
        amp,
        ["Ku_1_CDKiIAA_D4", "Ku_7_CDKiIAA_D4", "Ku_9_CDKiIAA_D4"],
        ["Ku_1_CDKi_D4", "Ku_7_CDKi_D4", "Ku_9_CDKi_D4"],
        clones, clones,
    )
    # acute 48h-equivalent: IAA day 2 vs DMSO (the unlabeled DMSO columns match the day-2 DMSO samples)
    clones = {
        "Ku_1_DMSO": "c1", "Ku_7_DMSO": "c7", "Ku_9_DMSO": "c9",
        "Ku_1_IAA_D2": "c1", "Ku_7_IAA_D2": "c7", "Ku_9_IAA_D2": "c9",
    }
    add(
        dict(
            contrast_id="GSE294709_HCT116_Ku80AID_IAA_D2",
            accession="GSE294709",
            organism="human",
            system="HCT116 Ku80-AID",
            perturbation="XRCC5",
            modality="AID protein degron",
            tier="timecourse",
            role="Ku loss",
            note="IAA+Dox day 2 vs DMSO, clones 1/7/9. Same AID system as the CDKi day-4 contrast; not an independent meta unit.",
        ),
        amp,
        ["Ku_1_IAA_D2", "Ku_7_IAA_D2", "Ku_9_IAA_D2"],
        ["Ku_1_DMSO", "Ku_7_DMSO", "Ku_9_DMSO"],
        clones, clones,
    )
    clones = {
        "Ku_1_DMSO_24h": "c1", "Ku_9_DMSO_24h": "c9",
        "Ku_1_IAA_24h": "c1", "Ku_9_IAA_24h": "c9",
    }
    add(
        dict(
            contrast_id="GSE294709_HCT116_Ku80AID_IAA_24h",
            accession="GSE294709",
            organism="human",
            system="HCT116 Ku80-AID",
            perturbation="XRCC5",
            modality="AID protein degron",
            tier="timecourse",
            role="Ku loss",
            note="IAA 24h vs DMSO 24h, clones 1 and 9 only (clone 7 has no 24h pair).",
        ),
        amp,
        ["Ku_1_IAA_24h", "Ku_9_IAA_24h"],
        ["Ku_1_DMSO_24h", "Ku_9_DMSO_24h"],
        clones, clones,
    )
    clones = {
        "Ku_1_CDKi_48h": "c1", "Ku_9_CDKi_48h": "c9",
        "Ku_1_CDKiIAA_48h": "c1", "Ku_9_CDKiIAA_48h": "c9",
    }
    add(
        dict(
            contrast_id="GSE294709_HCT116_Ku80AID_CDKi_48h",
            accession="GSE294709",
            organism="human",
            system="HCT116 Ku80-AID",
            perturbation="XRCC5",
            modality="AID protein degron + CDK4/6 inhibitor",
            tier="timecourse",
            role="Ku loss",
            note="CDKi+IAA 48h vs CDKi 48h, clones 1 and 9.",
        ),
        amp,
        ["Ku_1_CDKiIAA_48h", "Ku_9_CDKiIAA_48h"],
        ["Ku_1_CDKi_48h", "Ku_9_CDKi_48h"],
        clones, clones,
    )
    # MAVS-KO background: Ku still degraded, IFN should collapse if MAVS-dependent
    clones = {
        "MAVS1_DMSO": "c1", "MAVS13_DMSO": "c13", "MAVS15_DMSO": "c15",
        "MAVS1_IAA": "c1", "MAVS13_IAA": "c13", "MAVS15_IAA": "c15",
    }
    add(
        dict(
            contrast_id="GSE294709_HCT116_Ku80AID_MAVSKO",
            accession="GSE294709",
            organism="human",
            system="HCT116 Ku80-AID MAVS-KO",
            perturbation="XRCC5",
            modality="AID degron on MAVS knockout",
            tier="epistasis",
            role="Ku loss, MAVS absent",
            note="IAA vs DMSO in three MAVS-KO clones on the Ku80-AID background. Mechanistic control, not a primary Ku-loss unit.",
        ),
        amp,
        ["MAVS1_IAA", "MAVS13_IAA", "MAVS15_IAA"],
        ["MAVS1_DMSO", "MAVS13_DMSO", "MAVS15_DMSO"],
        clones, clones,
    )
    # flox Cre: Cre_D4 is 4OHT (XRCC5 mRNA collapses); Cre_D0 is EtOH. Confirmed by CPM.
    add(
        dict(
            contrast_id="GSE294709_HCT116_Ku80flox_Cre_4OHT",
            accession="GSE294709",
            organism="human",
            system="HCT116 Ku86flox CreERT2",
            perturbation="XRCC5",
            modality="4OHT-Cre excision",
            tier="unadjusted",
            role="Ku loss",
            note="Ku86flox/-;CreERT2, 4OHT day 4 (Cre_D4) vs EtOH (Cre_D0). Not the primary flox unit: the matched no-Cre 4OHT arm also moves some ISGs, so the Ku-specific estimate is the interaction contrast.",
        ),
        amp,
        ["Cre_D4", "Cre_D4_2"],
        ["Cre_D0", "Cre_D0_2"],
        None, None,
    )
    add(
        dict(
            contrast_id="GSE294709_HCT116_Ku80flox_noCre_4OHT",
            accession="GSE294709",
            organism="human",
            system="HCT116 Ku86flox no Cre",
            perturbation="none (4OHT control)",
            modality="4OHT without Cre",
            tier="negative_control",
            role="4OHT only",
            note="Same 4OHT vs EtOH in the floxed line that lacks Cre. XRCC5 should not drop. Not a Ku-loss contrast.",
        ),
        amp,
        ["Ctrl_D4", "Ctrl_D4_2"],
        ["Ctrl_D0", "Ctrl_D0_2"],
        None, None,
    )
    # HEK293 Ku70-KO with dox-inducible Ku70 rescue. D8 = 8 days without dox = depleted.
    clones = {
        "SB_Dox": "SB", "Sa11_Dox": "Sa11", "TI_Dox": "TI",
        "SB_D8": "SB", "Sa11_D8": "Sa11", "TI_D8": "TI",
    }
    add(
        dict(
            contrast_id="GSE294709_HEK293_Ku70_dox_withdrawal",
            accession="GSE294709",
            organism="human",
            system="HEK293 Ku70-KO + doxycycline-inducible Ku70",
            perturbation="XRCC6",
            modality="dox withdrawal of ectopic Ku70",
            tier="primary",
            role="Ku loss",
            note="Three clones (Sa11, SB, TI). No-dox day 8 (D8) vs dox (rescue). Series design: without dox the cells are Ku70-depleted.",
        ),
        kal,
        ["SB_D8", "Sa11_D8", "TI_D8"],
        ["SB_Dox", "Sa11_Dox", "TI_Dox"],
        clones, clones,
    )

    # --- GSE180581 HEK293T siRNA on heterozygous KO backgrounds ---
    def hek_add(cid, ko_cols, ctrl_cols, pert, tier, role, note, modality):
        add(
            dict(
                contrast_id=cid,
                accession="GSE180581",
                organism="human",
                system="HEK293T",
                perturbation=pert,
                modality=modality,
                tier=tier,
                role=role,
                note=note,
            ),
            hek,
            ko_cols,
            ctrl_cols,
            None, None,
        )

    hek_add(
        "GSE180581_HEK293T_siKu70",
        ["A07", "A08", "A09"],
        ["A04", "A05", "A06"],
        "XRCC6",
        "primary",
        "Ku loss",
        "Ku70+/- cells, siKu70 vs siControl, n=3. Matched heterozygous background (Anisenko et al., Data in Brief 2021).",
        "siRNA on heterozygous KO",
    )
    hek_add(
        "GSE180581_HEK293T_siKu80",
        ["A13", "A14", "A15"],
        ["A10", "A11", "A12"],
        "XRCC5",
        "primary",
        "Ku loss",
        "Ku80+/- cells, siKu80 vs siControl, n=3.",
        "siRNA on heterozygous KO",
    )
    hek_add(
        "GSE180581_HEK293T_siDNAPKcs",
        ["A19", "A20", "A21"],
        ["A16", "A17", "A18"],
        "PRKDC",
        "comparator",
        "DNA-PKcs loss, not Ku",
        "DNA-PKcs+/- cells, siDNA-PKcs vs siControl, n=3. Same complex, not a Ku gene. Comparator only.",
        "siRNA on heterozygous KO",
    )
    hek_add(
        "GSE180581_HEK293T_Ku70het_vs_WT",
        ["A04", "A05", "A06"],
        ["A01", "A02", "A03"],
        "XRCC6",
        "secondary",
        "heterozygous only",
        "Ku70+/- siControl vs WT siControl. Monoallelic, no further knockdown.",
        "heterozygous deletion",
    )
    hek_add(
        "GSE180581_HEK293T_Ku80het_vs_WT",
        ["A10", "A11", "A12"],
        ["A01", "A02", "A03"],
        "XRCC5",
        "secondary",
        "heterozygous only",
        "Ku80+/- siControl vs WT siControl. Monoallelic, no further knockdown. Shares the WT arm with the Ku70 het contrast.",
        "heterozygous deletion",
    )

    # --- GSE247031 mouse Treg Xrcc6 cKO ---
    add(
        dict(
            contrast_id="GSE247031_Treg_Xrcc6_cKO",
            accession="GSE247031",
            organism="mouse",
            system="mouse lung-tumor Treg, Xrcc6 flox Foxp3-Cre",
            perturbation="XRCC6",
            modality="Treg-specific genetic knockout",
            tier="primary",
            role="Ku loss",
            note="Tumor-infiltrating Tregs, Xrcc6 cKO vs Foxp3-Cre WT, n=3, LLC lung model.",
        ),
        treg,
        ["cKO_1", "cKO_2", "cKO_3"],
        ["WT_1", "WT_2", "WT_3"],
        None, None,
    )

    # --- ENCODE shRNA / CRISPR vs matched non-targeting controls ---
    encode = [
        ("ENCODE_K562_shXRCC5", "k562_shXRCC5", "k562_shNT", "XRCC5", "K562", "shRNA", "GSE88488", "ENCSR715XZS vs ENCSR815CVQ"),
        ("ENCODE_K562_shXRCC6", "k562_shXRCC6", "k562_shNT", "XRCC6", "K562", "shRNA", "GSE88126", "ENCSR232CPD vs ENCSR815CVQ"),
        ("ENCODE_HepG2_shXRCC5", "hepg2_shXRCC5", "hepg2_shNT", "XRCC5", "HepG2", "shRNA", "GSE80921", "ENCSR732IYM vs ENCSR491FOC"),
        ("ENCODE_HepG2_shXRCC6", "hepg2_shXRCC6", "hepg2_shNT", "XRCC6", "HepG2", "shRNA", "GSE80894", "ENCSR500WHE vs ENCSR491FOC"),
        ("ENCODE_K562_crXRCC5", "k562_crXRCC5", "k562_crNT", "XRCC5", "K562", "CRISPR", "GSE176960", "ENCSR276GMG vs ENCSR292PXV"),
        ("ENCODE_K562_crXRCC6", "k562_crXRCC6", "k562_crNT2", "XRCC6", "K562", "CRISPR", "GSE177139", "ENCSR516EPT vs ENCSR404OHQ"),
        ("ENCODE_HepG2_crXRCC6", "hepg2_crXRCC6", "hepg2_crNT", "XRCC6", "HepG2", "CRISPR", "GSE177121", "ENCSR312VLS vs ENCSR964HKT"),
    ]
    for cid, ko_p, ctrl_p, pert, line, modality, gse, note in encode:
        mat = load_encode_pair(ko_p, ctrl_p, ens_map, ent_map)
        add(
            dict(
                contrast_id=cid,
                accession=gse,
                organism="human",
                system=f"{line} ENCODE {modality}",
                perturbation=pert,
                modality=f"ENCODE {modality}",
                tier="encode",
                role="Ku loss",
                note=note + ". Gene quantifications GRCh38 V29 expected_count, rounded to integers. n=2 vs 2. K562 shXRCC5 and shXRCC6 share one non-targeting control.",
            ),
            mat,
            ["ko1", "ko2"],
            ["ctrl1", "ctrl2"],
            None, None,
        )
    return specs


def on_target_pass(row: dict, modality: str, perturbation: str) -> tuple[bool, str]:
    """mRNA on-target rule. AID degrons are protein-level and are exempt."""
    if "AID" in modality or "degron" in modality:
        return True, "protein degron; mRNA drop not required"
    if perturbation not in ("XRCC5", "XRCC6"):
        return True, "not a Ku-gene contrast"
    lfc = row.get("log2FoldChange", np.nan)
    padj = row.get("padj", np.nan)
    if pd.isna(lfc):
        return False, "target gene absent from DE table"
    if lfc < -0.4 or (lfc < 0 and pd.notna(padj) and padj < 0.05):
        return True, "mRNA down"
    return False, f"mRNA not convincingly down (log2FC={lfc:.3f})"


def gsea_p_up(nes: float, nom_p: float) -> float:
    """Convert a same-tail nominal p into a p-value for the UP direction."""
    p = float(nom_p) if nom_p and float(nom_p) > 0 else 1.0 / 1001.0
    p = min(max(p, 1.0 / 1001.0), 1.0)
    if nes >= 0:
        return p
    return 1.0 - p


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "de").mkdir(exist_ok=True)
    print("loading matrices", flush=True)
    human_sets = load_human_sets()
    (OUT / "genesets_human.json").write_text(json.dumps(human_sets, indent=1))
    print({k: len(v) for k, v in human_sets.items()})

    all_genes = sorted({g for gs in human_sets.values() for g in gs} | set(FOCUS))
    print("fetching mouse orthologs", flush=True)
    orth = fetch_mouse_orthologs(all_genes)
    mouse_sets = to_mouse_sets(human_sets, orth)
    (OUT / "genesets_mouse.json").write_text(json.dumps(mouse_sets, indent=1))

    amp = load_ampseq()
    kal = load_kallisto()
    hek = load_hek()
    treg = load_treg()
    ens_map, ent_map = load_hgnc_maps()
    print("shapes", amp.shape, kal.shape, hek.shape, treg.shape, flush=True)

    specs = build_contrasts(amp, kal, hek, treg, ens_map, ent_map)
    specs.append(dict(
        contrast_id="GSE294709_HCT116_Ku80flox_interaction",
        accession="GSE294709",
        organism="human",
        system="HCT116 Ku86flox CreERT2 vs no-Cre",
        perturbation="XRCC5",
        modality="4OHT-Cre excision, interaction vs no-Cre 4OHT",
        tier="primary",
        role="Ku loss",
        note="Interaction: extra log2FC of 4OHT in CreERT2 cells beyond 4OHT in the no-Cre floxed line. This is the Ku-excision effect after subtracting the 4OHT-only effect.",
        n_ko=2,
        n_ctrl=2,
        counts=None,
        meta=None,
        precomputed=True,
    ))
    gsea_rows = []
    focus_rows = []
    inventory = []

    for spec in specs:
        cid = spec["contrast_id"]
        print(f"\n=== {cid} ===", flush=True)
        de_path = OUT / "de" / f"{cid}.tsv.gz"
        if spec.get("precomputed"):
            paired = False
            if de_path.exists():
                res = pd.read_csv(de_path, sep="\t", index_col=0)
                print("  loaded cached interaction DE", flush=True)
            else:
                res = run_flox_interaction(amp)
                res.to_csv(de_path, sep="\t", compression="gzip")
        else:
            paired = spec["meta"]["clone"].nunique() > 1 and spec["meta"].groupby("clone")["condition"].nunique().min() == 2
            if spec["tier"] == "encode" or spec["meta"]["clone"].nunique() == 1:
                paired = False
            if de_path.exists():
                res = pd.read_csv(de_path, sep="\t", index_col=0)
                print("  loaded cached DE", flush=True)
            else:
                res = run_deseq(spec["counts"], spec["meta"], paired=paired)
                res.to_csv(de_path, sep="\t", compression="gzip")

        target_gene = spec["perturbation"] if spec["perturbation"] in ("XRCC5", "XRCC6", "PRKDC") else "XRCC5"
        tgt = gene_row(res, target_gene if target_gene != "PRKDC" else "PRKDC")
        # always also record both Ku genes
        xrcc5 = gene_row(res, "XRCC5")
        xrcc6 = gene_row(res, "XRCC6")
        passed, qc_note = on_target_pass(
            xrcc5 if spec["perturbation"] == "XRCC5" else xrcc6 if spec["perturbation"] == "XRCC6" else tgt,
            spec["modality"],
            spec["perturbation"],
        )
        in_primary_meta = spec["tier"] == "primary" and passed and spec["role"] == "Ku loss"
        # ENCODE joins a secondary meta only if on-target passes
        in_encode_meta = spec["tier"] == "encode" and passed
        print(
            f"  paired={paired} XRCC5 lfc={xrcc5['log2FoldChange']} XRCC6 lfc={xrcc6['log2FoldChange']} qc={passed} {qc_note}",
            flush=True,
        )
        inventory.append({
            "contrast_id": cid,
            "accession": spec["accession"],
            "organism": spec["organism"],
            "system": spec["system"],
            "perturbation": spec["perturbation"],
            "modality": spec["modality"],
            "tier": spec["tier"],
            "role": spec["role"],
            "n_ko": spec["n_ko"],
            "n_ctrl": spec["n_ctrl"],
            "paired_design": paired,
            "n_genes_tested": int(res["stat"].notna().sum()),
            "XRCC5_log2FC": xrcc5["log2FoldChange"],
            "XRCC5_padj": xrcc5["padj"],
            "XRCC6_log2FC": xrcc6["log2FoldChange"],
            "XRCC6_padj": xrcc6["padj"],
            "on_target_pass": passed,
            "qc_note": qc_note,
            "in_primary_meta": in_primary_meta,
            "in_encode_meta": in_encode_meta,
            "note": spec["note"],
        })

        sets = human_sets if spec["organism"] == "human" else mouse_sets
        print("  GSEA", flush=True)
        gtab = prerank(res, sets, res.index)
        # normalize column names across gseapy versions
        colmap = {c.lower(): c for c in gtab.columns}
        def col(*names):
            for n in names:
                if n in gtab.columns:
                    return n
                if n.lower() in colmap:
                    return colmap[n.lower()]
            raise KeyError(names)
        term_c = col("Term")
        nes_c = col("NES")
        p_c = col("NOM p-val", "NOM p-val", "pval")
        fdr_c = col("FDR q-val", "FDR q-val", "fdr")
        es_c = col("ES")
        lead_c = col("Lead_genes") if "Lead_genes" in gtab.columns or "lead_genes" in colmap else None
        for _, row in gtab.iterrows():
            nes = float(row[nes_c])
            nom = float(row[p_c])
            gsea_rows.append({
                "contrast_id": cid,
                "accession": spec["accession"],
                "organism": spec["organism"],
                "tier": spec["tier"],
                "role": spec["role"],
                "perturbation": spec["perturbation"],
                "on_target_pass": passed,
                "in_primary_meta": in_primary_meta,
                "in_encode_meta": in_encode_meta,
                "n_ko": spec["n_ko"],
                "n_ctrl": spec["n_ctrl"],
                "geneset": row[term_c],
                "ES": float(row[es_c]),
                "NES": nes,
                "nom_p": nom,
                "fdr_q": float(row[fdr_c]),
                "p_up": gsea_p_up(nes, nom),
                "n_in_set_resolved": int(row["n_in_set_resolved"]) if "n_in_set_resolved" in gtab.columns else np.nan,
                "lead_genes": row[lead_c] if lead_c else "",
            })

        for gene in FOCUS:
            fr = gene_row(res, gene)
            fr.update({
                "contrast_id": cid,
                "organism": spec["organism"],
                "tier": spec["tier"],
                "in_primary_meta": in_primary_meta,
            })
            focus_rows.append(fr)

    inv = pd.DataFrame(inventory)
    gsea = pd.DataFrame(gsea_rows)
    focus = pd.DataFrame(focus_rows)
    inv.to_csv(OUT / "contrast_inventory.tsv", sep="\t", index=False)
    gsea.to_csv(OUT / "gsea_prerank.tsv", sep="\t", index=False)
    focus.to_csv(OUT / "focus_genes.tsv", sep="\t", index=False)

    # Meta-analysis
    meta_rows = []
    for label, mask in [
        ("primary_Ku_loss", gsea["in_primary_meta"]),
        ("primary_plus_encode_ontarget", gsea["in_primary_meta"] | gsea["in_encode_meta"]),
        ("encode_ontarget_only", gsea["in_encode_meta"]),
        ("timecourse_AID_only", (gsea["tier"] == "timecourse") & (gsea["role"] == "Ku loss")),
    ]:
        sub = gsea.loc[mask]
        for gs, gsub in sub.groupby("geneset"):
            if gsub.empty:
                continue
            Z, p = stouffer(gsub["p_up"].values, np.sqrt(gsub["n_ko"] + gsub["n_ctrl"]))
            # also the DOWN direction
            p_down = 1 - gsub["p_up"].clip(1e-12, 1 - 1e-12)
            Zd, pdn = stouffer(p_down.values, np.sqrt(gsub["n_ko"] + gsub["n_ctrl"]))
            meta_rows.append({
                "meta": label,
                "geneset": gs,
                "n_contrasts": int(len(gsub)),
                "n_NES_positive": int((gsub["NES"] > 0).sum()),
                "n_NES_negative": int((gsub["NES"] < 0).sum()),
                "n_FDR_lt_0.25_up": int(((gsub["NES"] > 0) & (gsub["fdr_q"] < 0.25)).sum()),
                "median_NES": float(gsub["NES"].median()),
                "min_NES": float(gsub["NES"].min()),
                "max_NES": float(gsub["NES"].max()),
                "stouffer_Z_up": Z,
                "stouffer_p_up": p,
                "stouffer_Z_down": Zd,
                "stouffer_p_down": pdn,
                "contrasts": ";".join(gsub["contrast_id"].tolist()),
            })
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(OUT / "meta_stouffer.tsv", sep="\t", index=False)

    make_plots(inv, gsea)
    print("\nWrote", OUT)
    print(meta.to_string(index=False))
    print(inv.loc[inv.tier == "primary", ["contrast_id", "XRCC5_log2FC", "XRCC6_log2FC", "on_target_pass"]].to_string(index=False))


def make_plots(inv: pd.DataFrame, gsea: pd.DataFrame):
    # Heatmap of NES for primary + encode on-target + epistasis + DNA-PKcs + noCre
    show_ids = inv.loc[
        inv["tier"].isin(["primary", "encode", "epistasis", "comparator", "negative_control", "timecourse"])
        & ~((inv["tier"] == "encode") & (~inv["on_target_pass"]))
    , "contrast_id"].tolist()
    # keep a readable subset: primary, epistasis, comparator, negative control, and encode that passed
    sets_order = ["IFN_ALPHA", "IFN_GAMMA", "STING", "APM_MHCI", "APM_KEGG", "RIGI_MDA5"]
    sub = gsea[gsea["contrast_id"].isin(show_ids)].copy()
    # order contrasts
    order = [c for c in inv["contrast_id"] if c in set(show_ids)]
    mat = (
        sub.pivot(index="contrast_id", columns="geneset", values="NES")
        .reindex(index=order, columns=sets_order)
    )
    fdr = (
        sub.pivot(index="contrast_id", columns="geneset", values="fdr_q")
        .reindex(index=order, columns=sets_order)
    )
    fig_h = max(6, 0.38 * len(mat) + 1.5)
    fig, ax = plt.subplots(figsize=(9.4, fig_h))
    data = mat.values.astype(float)
    vmax = np.nanmax(np.abs(data)) if np.isfinite(data).any() else 1
    vmax = max(float(vmax), 1.0)
    im = ax.imshow(data, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(sets_order)))
    ax.set_xticklabels(["IFN-α", "IFN-γ", "STING", "MHC-I APM", "APM KEGG", "RIG-I/MDA5"], fontsize=9)
    ax.set_yticks(range(len(mat.index)))
    labels = []
    inv_i = inv.set_index("contrast_id")
    for cid in mat.index:
        tier = inv_i.loc[cid, "tier"]
        labels.append(f"{cid}  [{tier}]")
    ax.set_yticklabels(labels, fontsize=7)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            q = fdr.values[i, j]
            if not np.isfinite(val):
                continue
            star = ""
            if np.isfinite(q) and q < 0.05:
                star = "*"
            elif np.isfinite(q) and q < 0.25:
                star = "·"
            ax.text(j, i, f"{val:.2f}{star}", ha="center", va="center", fontsize=6.5,
                    color="black" if abs(val) < 0.65 * vmax else "white")
    ax.set_title("Pre-ranked GSEA NES after Ku70/Ku80 loss\n* FDR<0.05   · FDR<0.25   positive = up in KO/KD")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="NES")
    fig.tight_layout()
    fig.savefig(OUT / "fig_nes_heatmap.png", dpi=160)
    fig.savefig(OUT / "fig_nes_heatmap.pdf")
    plt.close()

    # On-target bars for primary + encode
    ont = inv[inv["tier"].isin(["primary", "encode", "epistasis", "comparator", "negative_control"])].copy()
    fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(ont) + 1)))
    y = np.arange(len(ont))
    ax.barh(y - 0.15, ont["XRCC5_log2FC"], height=0.3, label="XRCC5 (Ku80)", color="#1f4e79")
    ax.barh(y + 0.15, ont["XRCC6_log2FC"], height=0.3, label="XRCC6 (Ku70)", color="#c45911")
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(ont["contrast_id"], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("DESeq2 log2 fold-change (KO/KD vs control)")
    ax.set_title("On-target check: Ku70 / Ku80 mRNA")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_ontarget_log2fc.png", dpi=160)
    fig.savefig(OUT / "fig_ontarget_log2fc.pdf")
    plt.close()

    # Forest of IFN-alpha NES for primary meta members
    prim = gsea[(gsea["in_primary_meta"]) & (gsea["geneset"] == "IFN_ALPHA")].copy()
    if len(prim):
        fig, ax = plt.subplots(figsize=(7.2, 0.45 * len(prim) + 1.4))
        y = np.arange(len(prim))
        colors = ["#b2182b" if v > 0 else "#2166ac" for v in prim["NES"]]
        ax.scatter(prim["NES"], y, c=colors, s=40, zorder=3)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(prim["contrast_id"], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("NES, Hallmark interferon-alpha response")
        ax.set_title("Primary Ku-loss contrasts")
        fig.tight_layout()
        fig.savefig(OUT / "fig_ifn_alpha_primary_nes.png", dpi=160)
        fig.savefig(OUT / "fig_ifn_alpha_primary_nes.pdf")
        plt.close()


if __name__ == "__main__":
    main()
