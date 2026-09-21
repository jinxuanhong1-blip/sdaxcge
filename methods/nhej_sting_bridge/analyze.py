#!/usr/bin/env python3
"""NHEJ-component loss versus IFN / STING / antigen-presentation transcripts.

Two public RNA-seq series, reanalyzed from the deposited counts:

* GSE180581 — HEK293T monoallelic KO plus siRNA of Ku70 (XRCC6), Ku80 (XRCC5),
  or DNA-PKcs (PRKDC), each versus siControl in the same heterozygous line
  (n = 3). Anisenko et al., Data in Brief 2021 / Biochimie 2022.
* GSE135274 — HeLa XRCC4(-/-) versus scrambled-gRNA wild type, without mirin
  (n = 2 experiments). Read-1 and read-2 columns are technical and are summed.
  Every sample was also transfected with a TALEN NHEJ reporter. Benjamin &
  Schiller, Int J Mol Sci 2022.

Question: does loss of an NHEJ component shift IFN, STING-axis, or MHC-I APM
transcripts upward? That is the middle step of CLDN4 → NHEJ → IFN. These
series do not perturb CLDN4.

Counts are tested with PyDESeq2 (Wald, BH padj). Panel direction uses the
shrunken-free MLE log2 fold change among genes that pass the count filter
(≥10 counts in ≥ the smaller group size). A second, labeled log2FC includes
genes below that filter so a silent floor is visible. Hallmark IFN-α and IFN-γ
are tested by pre-ranked GSEA on the Wald statistic.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from gene_sets import APM, CORE8, FOCUS_DISPLAY, IFN_CORE, NHEJ, PANELS, STING

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = ROOT / "results" / "nhej_sting_bridge"
GMT = HERE / "gene_sets" / "hallmark_ifn.gmt"

GSE180_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180581/suppl/"
    "GSE180581_all_samples_counts.xlsx"
)
GSE135_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135274/suppl/"
    "GSE135274_all_sample_human_RNA_max.txt.gz"
)
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
    "hgnc_complete_set.txt"
)

# Series-matrix sample order A1..A21 (xlsx columns A01..A21).
GSE180_SAMPLES = [
    ("A01", "parental", "siControl"),
    ("A02", "parental", "siControl"),
    ("A03", "parental", "siControl"),
    ("A04", "Ku70_het", "siControl"),
    ("A05", "Ku70_het", "siControl"),
    ("A06", "Ku70_het", "siControl"),
    ("A07", "Ku70_het", "siKu70"),
    ("A08", "Ku70_het", "siKu70"),
    ("A09", "Ku70_het", "siKu70"),
    ("A10", "Ku80_het", "siControl"),
    ("A11", "Ku80_het", "siControl"),
    ("A12", "Ku80_het", "siControl"),
    ("A13", "Ku80_het", "siKu80"),
    ("A14", "Ku80_het", "siKu80"),
    ("A15", "Ku80_het", "siKu80"),
    ("A16", "DNAPKcs_het", "siControl"),
    ("A17", "DNAPKcs_het", "siControl"),
    ("A18", "DNAPKcs_het", "siControl"),
    ("A19", "DNAPKcs_het", "siDNAPKcs"),
    ("A20", "DNAPKcs_het", "siDNAPKcs"),
    ("A21", "DNAPKcs_het", "siDNAPKcs"),
]

# Primary contrasts: depleted versus matched siControl in the same line.
GSE180_CONTRASTS = [
    ("GSE180581_Ku70", "XRCC6", "Ku70_het", "siKu70", "siControl"),
    ("GSE180581_Ku80", "XRCC5", "Ku80_het", "siKu80", "siControl"),
    ("GSE180581_DNAPKcs", "PRKDC", "DNAPKcs_het", "siDNAPKcs", "siControl"),
]


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    print(f"downloading {url}")
    urllib.request.urlretrieve(url, dest)


def load_hgnc(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    return df


def alias_map(hgnc: pd.DataFrame) -> dict[str, str]:
    """Map a unique previous/alias symbol onto the current symbol.

    A token that is itself a current symbol is left alone. Tokens that point
    at two different current symbols are dropped.
    """
    current = set(hgnc["symbol"].dropna())
    mapping: dict[str, str] = {}
    conflict: set[str] = set()
    for sym, prev, alias in hgnc[["symbol", "prev_symbol", "alias_symbol"]].itertuples(index=False):
        if not isinstance(sym, str):
            continue
        blob = " ".join(x for x in (prev, alias) if isinstance(x, str))
        tokens = [t.strip() for t in blob.replace("|", ",").replace(";", ",").split(",") if t.strip()]
        for tok in tokens:
            if tok == sym or tok in current:
                continue
            if tok in mapping and mapping[tok] != sym:
                conflict.add(tok)
            else:
                mapping[tok] = sym
    for tok in conflict:
        mapping.pop(tok, None)
    # Explicit sensor renames even if HGNC alias parsing misses them.
    mapping.setdefault("TMEM173", "STING1")
    mapping.setdefault("MB21D1", "CGAS")
    mapping.setdefault("C9orf142", "PAXX")
    return mapping


def collapse_symbols(counts: pd.DataFrame, rename: dict[str, str]) -> pd.DataFrame:
    """counts: genes x samples, index = raw symbol."""
    idx = [rename.get(g, g) for g in counts.index]
    out = counts.copy()
    out.index = idx
    out = out.groupby(level=0).sum()
    out = out.loc[out.sum(axis=1) > 0]
    return out


def ensembl_to_symbol(hgnc: pd.DataFrame) -> dict[str, str]:
    sub = hgnc.dropna(subset=["ensembl_gene_id", "symbol"])
    # one symbol per Ensembl id; first approved row wins
    return dict(zip(sub["ensembl_gene_id"], sub["symbol"]))


def load_gse180(path: Path, symbol_of: dict[str, str]) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=0)
    raw = raw.rename(columns={raw.columns[0]: "ensembl"})
    raw["ensembl"] = raw["ensembl"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    raw["symbol"] = raw["ensembl"].map(symbol_of)
    mapped = raw.dropna(subset=["symbol"]).copy()
    sample_cols = [c for c, _, _ in GSE180_SAMPLES]
    # xlsx uses A01..A21
    missing = [c for c in sample_cols if c not in mapped.columns]
    if missing:
        raise SystemExit(f"GSE180581 missing columns: {missing}")
    for c in sample_cols:
        mapped[c] = pd.to_numeric(mapped[c], errors="coerce").fillna(0).astype(int)
    mat = mapped.groupby("symbol")[sample_cols].sum()
    n_ens = raw["ensembl"].nunique()
    n_mapped = raw["symbol"].notna().sum()
    print(f"GSE180581 genes {n_ens}, rows with HGNC symbol {n_mapped}, collapsed symbols {mat.shape[0]}")
    return mat


def load_gse135(path: Path, rename: dict[str, str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip", dtype=str)
    if "geneName" not in df.columns:
        raise SystemExit(f"unexpected GSE135274 header: {df.columns.tolist()[:5]}")
    count_cols = [c for c in df.columns if c not in ("matchType", "geneName")]
    for c in count_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
    df["symbol"] = df["geneName"].map(lambda g: rename.get(g, g) if isinstance(g, str) else g)
    mat = df.groupby("symbol")[count_cols].sum()
    # Sum technical read-1 / read-2 into one biological sample.
    groups = {
        "exp1_XRCC4": ["Experiment_01_g2G3_1_1", "Experiment_01_g2G3_1_2"],
        "exp1_XRCC4_mirin": ["Experiment_01_g2G3_M_1_1", "Experiment_01_g2G3_M_1_2"],
        "exp1_WT": ["Experiment_01_gSCR_1_1", "Experiment_01_gSCR_1_2"],
        "exp1_WT_mirin": ["Experiment_01_gSCR_M_1_1", "Experiment_01_gSCR_M_1_2"],
        "exp2_XRCC4": ["Experiment_02_g2G3_2_1", "Experiment_02_g2G3_2_2"],
        "exp2_XRCC4_mirin": ["Experiment_02_g2G3_M_2_1", "Experiment_02_g2G3_M_2_2"],
        "exp2_WT": ["Experiment_02_gSCR_2_1", "Experiment_02_gSCR_2_2"],
        "exp2_WT_mirin": ["Experiment_02_gSCR_M_2_1", "Experiment_02_gSCR_M_2_2"],
    }
    missing = [c for cols in groups.values() for c in cols if c not in mat.columns]
    if missing:
        raise SystemExit(f"GSE135274 missing columns: {missing}")
    bio = pd.DataFrame({name: mat[cols].sum(axis=1) for name, cols in groups.items()})
    bio = bio.round().astype(int)
    print(f"GSE135274 collapsed symbols {bio.shape[0]}")
    print("library sizes:", {c: int(bio[c].sum()) for c in bio.columns})
    return bio


def load_hallmark(path: Path) -> dict[str, list[str]]:
    sets = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def median_ratio_sf(mat: pd.DataFrame) -> pd.Series:
    """DESeq2 median-of-ratios size factors. mat is genes x samples."""
    vals = mat.to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        logv = np.where(vals > 0, np.log(vals), np.nan)
    n_finite = np.sum(np.isfinite(logv), axis=1)
    geo = np.full(logv.shape[0], np.nan)
    ok = n_finite > 0
    geo[ok] = np.exp(np.nanmean(logv[ok], axis=1))
    usable = np.isfinite(geo) & (geo > 0)
    if usable.sum() < 50:
        totals = mat.sum(axis=0)
        return totals / totals.mean()
    ratio = vals[usable, :] / geo[usable][:, None]
    sf = np.nanmedian(ratio, axis=0)
    sf = np.where(sf <= 0, np.nan, sf)
    if np.isnan(sf).any():
        totals = mat.sum(axis=0).to_numpy(dtype=float)
        sf = totals / np.mean(totals)
    # rescale to geometric mean 1
    sf = sf / np.exp(np.mean(np.log(sf)))
    return pd.Series(sf, index=mat.columns)


def simple_log2fc(mat: pd.DataFrame, treat: list[str], ctrl: list[str], sf: pd.Series) -> pd.Series:
    norm = mat.div(sf, axis=1)
    a = norm[treat].mean(axis=1)
    b = norm[ctrl].mean(axis=1)
    return np.log2(a + 1.0) - np.log2(b + 1.0)


def run_deseq(mat: pd.DataFrame, treat: list[str], ctrl: list[str]) -> pd.DataFrame:
    """Return gene-level DESeq2 results. mat is genes x all samples."""
    sub_cols = treat + ctrl
    sub = mat[sub_cols]
    min_n = min(len(treat), len(ctrl))
    keep = (sub >= 10).sum(axis=1) >= min_n
    sub = sub.loc[keep]
    # drop genes that are all zero in the contrast (should already be gone)
    sub = sub.loc[sub.sum(axis=1) > 0]
    meta = pd.DataFrame(
        {"condition": ["depleted"] * len(treat) + ["control"] * len(ctrl)},
        index=sub_cols,
    )
    counts = sub.T.astype(int)
    last_err = None
    res = None
    for fit_type in ("parametric", "mean"):
        try:
            dds = DeseqDataSet(
                counts=counts,
                metadata=meta,
                design="~condition",
                refit_cooks=False,
                fit_type=fit_type,
                quiet=True,
                n_cpus=1,
            )
            dds.deseq2()
            st = DeseqStats(
                dds,
                contrast=["condition", "depleted", "control"],
                cooks_filter=False,
                independent_filter=True,
                quiet=True,
                n_cpus=1,
            )
            st.summary()
            res = st.results_df.copy()
            res["fit_type"] = fit_type
            print(f"  DESeq2 fit_type={fit_type} genes={res.shape[0]}")
            break
        except Exception as exc:  # noqa: BLE001 — fallback to mean dispersion
            last_err = exc
            print(f"  DESeq2 {fit_type} failed: {type(exc).__name__}: {exc}")
    if res is None:
        raise RuntimeError(f"DESeq2 failed: {last_err}")
    res.index.name = "symbol"
    return res


def panel_stats(lfc: pd.Series, symbols: list[str], background: pd.Series) -> dict:
    """One-sided tests that the panel log2FC is shifted upward."""
    present = []
    missing = []
    for g in dict.fromkeys(symbols):
        if g in lfc.index and np.isfinite(lfc.at[g]):
            present.append(float(lfc.at[g]))
        else:
            missing.append(g)
    arr = np.array(present, dtype=float)
    n = int(arr.size)
    n_up = int((arr > 0).sum()) if n else 0
    out = {
        "n_panel": len(list(dict.fromkeys(symbols))),
        "n_tested": n,
        "n_missing": len(missing),
        "n_up": n_up,
        "n_down": int((arr < 0).sum()) if n else 0,
        "n_zero": int((arr == 0).sum()) if n else 0,
        "median_log2FC": None if n == 0 else round(float(np.median(arr)), 4),
        "mean_log2FC": None if n == 0 else round(float(np.mean(arr)), 4),
        "sign_test_p_up": None,
        "wilcoxon_p_up": None,
        "wilcoxon_p_down": None,
        "mwu_p_panel_gt_bg": None,
        "missing_symbols": ",".join(missing),
    }
    if n >= 5:
        out["sign_test_p_up"] = float(stats.binomtest(n_up, n, 0.5, alternative="greater").pvalue)
        # Wilcoxon needs non-identical values
        if np.nanstd(arr) > 0 and not np.allclose(arr, 0):
            try:
                out["wilcoxon_p_up"] = float(stats.wilcoxon(arr, alternative="greater", zero_method="wilcox").pvalue)
                out["wilcoxon_p_down"] = float(stats.wilcoxon(arr, alternative="less", zero_method="wilcox").pvalue)
            except ValueError:
                out["wilcoxon_p_up"] = None
                out["wilcoxon_p_down"] = None
        bg = background.replace([np.inf, -np.inf], np.nan).dropna()
        bg = bg[~bg.index.isin(list(dict.fromkeys(symbols)))]
        if len(bg) > 20 and n >= 5:
            out["mwu_p_panel_gt_bg"] = float(stats.mannwhitneyu(arr, bg.to_numpy(), alternative="greater").pvalue)
    return out


def prerank_gsea(stat: pd.Series, geneset: list[str], nperm: int = 2000, seed: int = 1) -> dict | None:
    s = stat.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    if s.empty:
        return None
    genes = s.index.to_numpy()
    weights = s.to_numpy(dtype=float)
    gset = set(geneset)
    hit = np.array([g in gset for g in genes])
    nh = int(hit.sum())
    n = int(len(genes))
    if nh < 8 or nh >= n:
        return {"n_genes": nh, "n_universe": n, "ES": None, "NES": None, "nominal_p": None, "note": "set too small or empty in universe"}
    abs_w = np.abs(weights)
    sum_hit = float((abs_w * hit).sum())
    if sum_hit <= 0:
        return None

    def es_of(hit_mask: np.ndarray) -> tuple[float, int]:
        hw = abs_w * hit_mask
        sh = hw.sum()
        step_hit = hw / sh
        step_miss = (~hit_mask).astype(float) / (n - nh)
        running = np.cumsum(np.where(hit_mask, step_hit, -step_miss))
        pos = float(running.max())
        neg = float(running.min())
        if abs(pos) >= abs(neg):
            return pos, int(np.argmax(running))
        return neg, int(np.argmin(running))

    es, peak = es_of(hit)
    rng = np.random.default_rng(seed)
    null = np.empty(nperm)
    for i in range(nperm):
        perm = rng.permutation(hit)
        null[i], _ = es_of(perm)
    if es >= 0:
        same = null[null > 0]
        nes = float(es / same.mean()) if len(same) else float("nan")
        p = float((np.sum(null >= es) + 1) / (nperm + 1))
    else:
        same = null[null < 0]
        nes = float(es / np.abs(same.mean())) if len(same) else float("nan")
        p = float((np.sum(null <= es) + 1) / (nperm + 1))
    # leading edge: members on the side of the peak
    if es >= 0:
        edge = [g for g, h in zip(genes[: peak + 1], hit[: peak + 1]) if h]
    else:
        edge = [g for g, h in zip(genes[peak:], hit[peak:]) if h]
    return {
        "n_genes": nh,
        "n_universe": n,
        "ES": round(es, 4),
        "NES": None if not np.isfinite(nes) else round(nes, 4),
        "nominal_p": p,
        "leading_edge_n": len(edge),
        "leading_edge": ",".join(edge[:25]),
    }


def replicate_cors(mat: pd.DataFrame, cols: list[str]) -> float | None:
    if len(cols) < 2:
        return None
    sf = median_ratio_sf(mat[cols])
    logc = np.log2(mat[cols].div(sf, axis=1) + 1.0)
    # mean pairwise Pearson on genes with some counts
    keep = mat[cols].sum(axis=1) >= 20
    x = logc.loc[keep]
    rhos = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            rhos.append(float(np.corrcoef(x.iloc[:, i], x.iloc[:, j])[0, 1]))
    return float(np.mean(rhos)) if rhos else None


def analyze_contrast(
    name: str,
    target: str,
    mat: pd.DataFrame,
    treat: list[str],
    ctrl: list[str],
    hallmark: dict[str, list[str]],
    note: str,
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    print(f"\n=== {name}  {treat} vs {ctrl}")
    res = run_deseq(mat, treat, ctrl)
    sf = median_ratio_sf(mat[treat + ctrl])
    simple = simple_log2fc(mat, treat, ctrl, sf)
    res = res.join(simple.rename("simple_log2FC_pseudocount"), how="left")
    # mean normalized counts so a floor is visible
    norm = mat[treat + ctrl].div(sf, axis=1)
    res = res.join(norm[treat].mean(axis=1).rename("mean_norm_depleted"), how="left")
    res = res.join(norm[ctrl].mean(axis=1).rename("mean_norm_control"), how="left")
    lfc = res["log2FoldChange"]
    bg = lfc

    panels = {}
    for pname, genes in PANELS.items():
        panels[pname] = panel_stats(lfc, genes, bg)
    # Hallmark on Wald stat
    gsea = {}
    stat = res["stat"] if "stat" in res.columns else lfc
    for hname, genes in hallmark.items():
        gsea[hname] = prerank_gsea(stat, genes, nperm=2000, seed=1)
        # also directional median on the overlap that was tested
        panels[hname] = panel_stats(lfc, genes, bg)

    # on target
    def gene_row(symbol: str) -> dict:
        if symbol not in res.index and symbol not in simple.index:
            return {"symbol": symbol, "status": "absent"}
        row = {"symbol": symbol, "status": "tested" if symbol in res.index else "below_count_filter"}
        if symbol in res.index:
            for col in ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "simple_log2FC_pseudocount", "mean_norm_depleted", "mean_norm_control"):
                if col in res.columns:
                    val = res.at[symbol, col]
                    row[col] = None if pd.isna(val) else float(val)
        else:
            row["simple_log2FC_pseudocount"] = None if symbol not in simple.index else float(simple.at[symbol])
            if symbol in norm.index:
                row["mean_norm_depleted"] = float(norm.loc[symbol, treat].mean())
                row["mean_norm_control"] = float(norm.loc[symbol, ctrl].mean())
        return row

    on_target = gene_row(target)
    nhej_rows = [gene_row(g) for g in NHEJ]
    focus_rows = [gene_row(g) for g in FOCUS_DISPLAY]

    summary = {
        "contrast": name,
        "target": target,
        "n_depleted": len(treat),
        "n_control": len(ctrl),
        "n_genes_tested": int(res.shape[0]),
        "n_padj_lt_0.05": int((res["padj"] < 0.05).sum()) if "padj" in res.columns else None,
        "n_up_padj": int(((res["padj"] < 0.05) & (res["log2FoldChange"] > 0)).sum()) if "padj" in res.columns else None,
        "n_down_padj": int(((res["padj"] < 0.05) & (res["log2FoldChange"] < 0)).sum()) if "padj" in res.columns else None,
        "replicate_pearson_depleted": replicate_cors(mat, treat),
        "replicate_pearson_control": replicate_cors(mat, ctrl),
        "on_target": on_target,
        "panels": panels,
        "gsea": gsea,
        "note": note,
        "fit_type": str(res["fit_type"].iloc[0]) if "fit_type" in res.columns else None,
    }
    focus = pd.DataFrame(focus_rows)
    focus.insert(0, "contrast", name)
    nhej = pd.DataFrame(nhej_rows)
    nhej.insert(0, "contrast", name)
    summary["_nhej_table"] = nhej
    return res, summary, focus


def write_hela_concordance(mat: pd.DataFrame) -> None:
    """Per-experiment log2FC for the XRCC4-null vs scramble arm (no mirin).

    Library-size factors, pseudocount 1. This does not use the DESeq2 fit,
    so a gene that moves in only one experiment is visible.
    """
    cols = {
        "exp1_ko": "exp1_XRCC4",
        "exp1_wt": "exp1_WT",
        "exp2_ko": "exp2_XRCC4",
        "exp2_wt": "exp2_WT",
    }
    sub = mat[list(cols.values())].astype(float)
    sf = sub.sum(axis=0)
    sf = sf / sf.mean()
    norm = sub.div(sf, axis=1)
    rows = []
    for gene in FOCUS_DISPLAY:
        if gene not in norm.index:
            rows.append({"symbol": gene, "status": "absent"})
            continue
        def lfc(ko, wt, g=gene):
            return float(np.log2((norm.at[g, ko] + 1.0) / (norm.at[g, wt] + 1.0)))
        e1 = lfc(cols["exp1_ko"], cols["exp1_wt"])
        e2 = lfc(cols["exp2_ko"], cols["exp2_wt"])
        rows.append({
            "symbol": gene,
            "status": "present",
            "exp1_log2FC": round(e1, 4),
            "exp2_log2FC": round(e2, 4),
            "same_sign": int(np.sign(e1) == np.sign(e2) and e1 != 0),
            "exp1_ko": int(sub.at[gene, cols["exp1_ko"]]),
            "exp1_wt": int(sub.at[gene, cols["exp1_wt"]]),
            "exp2_ko": int(sub.at[gene, cols["exp2_ko"]]),
            "exp2_wt": int(sub.at[gene, cols["exp2_wt"]]),
        })
    pd.DataFrame(rows).to_csv(OUT / "gse135274_experiment_concordance.tsv", sep="\t", index=False)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    _download(GSE180_URL, DATA / "GSE180581_all_samples_counts.xlsx")
    _download(GSE135_URL, DATA / "GSE135274_all_sample_human_RNA_max.txt.gz")
    _download(HGNC_URL, DATA / "hgnc_complete_set.txt")

    hgnc = load_hgnc(DATA / "hgnc_complete_set.txt")
    symbol_of = ensembl_to_symbol(hgnc)
    rename = alias_map(hgnc)
    hallmark = load_hallmark(GMT)

    mat180 = load_gse180(DATA / "GSE180581_all_samples_counts.xlsx", symbol_of)
    mat135 = load_gse135(DATA / "GSE135274_all_sample_human_RNA_max.txt.gz", rename)

    # library sizes for the HEK matrix
    print("GSE180581 library sizes:", {c: int(mat180[c].sum()) for c, _, _ in GSE180_SAMPLES})

    sample_rows = []
    for col, line, sirna in GSE180_SAMPLES:
        sample_rows.append({
            "series": "GSE180581",
            "sample": col,
            "line": line,
            "sirna": sirna,
            "library_size": int(mat180[col].sum()),
        })
    for col in mat135.columns:
        if "XRCC4_mirin" in col:
            geno, drug = "XRCC4_KO", "mirin"
        elif col.endswith("XRCC4"):
            geno, drug = "XRCC4_KO", "none"
        elif "WT_mirin" in col:
            geno, drug = "WT_scramble", "mirin"
        else:
            geno, drug = "WT_scramble", "none"
        sample_rows.append({
            "series": "GSE135274",
            "sample": col,
            "line": geno,
            "sirna": drug,
            "library_size": int(mat135[col].sum()),
        })
    pd.DataFrame(sample_rows).to_csv(OUT / "sample_sheet.tsv", sep="\t", index=False)

    summaries = []
    focus_all = []
    nhej_all = []

    for name, target, line, treat_s, ctrl_s in GSE180_CONTRASTS:
        treat = [c for c, ln, s in GSE180_SAMPLES if ln == line and s == treat_s]
        ctrl = [c for c, ln, s in GSE180_SAMPLES if ln == line and s == ctrl_s]
        note = (
            f"HEK293T {line} {treat_s} vs {ctrl_s}, 72 h, 50 nM siRNA, n=3. "
            "Matched heterozygous background. Parental 293T siControl is not the control."
        )
        res, summary, focus = analyze_contrast(name, target, mat180, treat, ctrl, hallmark, note)
        res.to_csv(OUT / f"{name}_de.tsv.gz", sep="\t", compression="gzip")
        summaries.append(summary)
        focus_all.append(focus)
        nhej_all.append(summary.pop("_nhej_table"))

    he_la = [
        (
            "GSE135274_XRCC4",
            "XRCC4",
            ["exp1_XRCC4", "exp2_XRCC4"],
            ["exp1_WT", "exp2_WT"],
            "HeLa XRCC4(-/-) vs scrambled gRNA, no mirin. n=2 experiments. "
            "Read1+read2 summed. All samples received TALEN + NHEJ-reporter transfection.",
        ),
        (
            "GSE135274_XRCC4_within_mirin",
            "XRCC4",
            ["exp1_XRCC4_mirin", "exp2_XRCC4_mirin"],
            ["exp1_WT_mirin", "exp2_WT_mirin"],
            "Secondary: XRCC4 KO vs scramble, both arms treated with mirin (A-NHEJ block). n=2.",
        ),
        (
            "GSE135274_double_block_vs_WT",
            "XRCC4",
            ["exp1_XRCC4_mirin", "exp2_XRCC4_mirin"],
            ["exp1_WT", "exp2_WT"],
            "Secondary: XRCC4 KO + mirin vs untreated scramble. This is the paper's double-block arm, not pure XRCC4 loss. n=2.",
        ),
    ]
    for name, target, treat, ctrl, note in he_la:
        res, summary, focus = analyze_contrast(name, target, mat135, treat, ctrl, hallmark, note)
        res.to_csv(OUT / f"{name}_de.tsv.gz", sep="\t", compression="gzip")
        summaries.append(summary)
        focus_all.append(focus)
        nhej_all.append(summary.pop("_nhej_table"))

    # tables
    panel_rows = []
    gsea_rows = []
    on_rows = []
    for s in summaries:
        ot = s["on_target"]
        on_rows.append({
            "contrast": s["contrast"],
            "target": s["target"],
            "status": ot.get("status"),
            "log2FoldChange": ot.get("log2FoldChange"),
            "pvalue": ot.get("pvalue"),
            "padj": ot.get("padj"),
            "simple_log2FC": ot.get("simple_log2FC_pseudocount"),
            "mean_norm_depleted": ot.get("mean_norm_depleted"),
            "mean_norm_control": ot.get("mean_norm_control"),
            "n_genes_tested": s["n_genes_tested"],
            "n_padj_lt_0.05": s["n_padj_lt_0.05"],
            "replicate_pearson_depleted": s["replicate_pearson_depleted"],
            "replicate_pearson_control": s["replicate_pearson_control"],
            "fit_type": s["fit_type"],
        })
        for pname, p in s["panels"].items():
            panel_rows.append({
                "contrast": s["contrast"],
                "panel": pname,
                "n_panel": p["n_panel"],
                "n_tested": p["n_tested"],
                "n_missing": p["n_missing"],
                "n_up": p["n_up"],
                "n_down": p["n_down"],
                "median_log2FC": p["median_log2FC"],
                "mean_log2FC": p["mean_log2FC"],
                "sign_test_p_up": p["sign_test_p_up"],
                "wilcoxon_p_up": p["wilcoxon_p_up"],
                "wilcoxon_p_down": p["wilcoxon_p_down"],
                "mwu_p_panel_gt_bg": p["mwu_p_panel_gt_bg"],
                "missing_symbols": p["missing_symbols"],
            })
        for gname, g in s["gsea"].items():
            if g is None:
                continue
            gsea_rows.append({"contrast": s["contrast"], "geneset": gname, **g})

    pd.DataFrame(on_rows).to_csv(OUT / "on_target.tsv", sep="\t", index=False)
    pd.DataFrame(panel_rows).to_csv(OUT / "panel_summary.tsv", sep="\t", index=False)
    gsea_df = pd.DataFrame(gsea_rows)
    if len(gsea_df) and "nominal_p" in gsea_df.columns:
        raw = gsea_df["nominal_p"].to_numpy(dtype=float)
        m = len(gsea_df)
        bh = np.empty(m)
        order = np.argsort(raw, kind="mergesort")
        ranked = raw[order]
        adj = ranked * m / np.arange(1, m + 1)
        for i in range(m - 2, -1, -1):
            adj[i] = min(adj[i], adj[i + 1])
        bh[order] = np.minimum(adj, 1.0)
        gsea_df["bh_q"] = bh
    gsea_df.to_csv(OUT / "gsea_prerank.tsv", sep="\t", index=False)
    write_hela_concordance(mat135)
    pd.concat(focus_all, ignore_index=True).to_csv(OUT / "focus_genes.tsv", sep="\t", index=False)
    pd.concat(nhej_all, ignore_index=True).to_csv(OUT / "nhej_genes.tsv", sep="\t", index=False)

    # JSON without huge missing lists duplicated — keep them, file is small
    slim = []
    for s in summaries:
        sc = dict(s)
        slim.append(sc)
    (OUT / "summary.json").write_text(json.dumps(slim, indent=2, default=str))

    plot_results(pd.DataFrame(panel_rows), pd.DataFrame(on_rows), pd.concat(focus_all, ignore_index=True))
    print("\nWrote", OUT)


def plot_results(panels: pd.DataFrame, on_target: pd.DataFrame, focus: pd.DataFrame) -> None:
    primary = [
        "GSE180581_Ku70",
        "GSE180581_Ku80",
        "GSE180581_DNAPKcs",
        "GSE135274_XRCC4",
    ]
    labels = {
        "GSE180581_Ku70": "Ku70",
        "GSE180581_Ku80": "Ku80",
        "GSE180581_DNAPKcs": "DNA-PKcs",
        "GSE135274_XRCC4": "XRCC4",
    }
    colors = {
        "GSE180581_Ku70": "#2F5D7C",
        "GSE180581_Ku80": "#5B8FA8",
        "GSE180581_DNAPKcs": "#C47B3A",
        "GSE135274_XRCC4": "#8C3A3A",
    }
    show_panels = ["STING", "IFN_CORE", "APM", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"]
    # gmt names
    # actual keys from file: INTERFERON_ALPHA_RESPONSE
    have = set(panels["panel"])
    show_panels = [p for p in [
        "STING", "IFN_CORE", "CORE8", "APM",
        "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE",
    ] if p in have]
    short = {
        "STING": "STING axis",
        "IFN_CORE": "IFN core",
        "CORE8": "Core-8",
        "APM": "APM",
        "INTERFERON_ALPHA_RESPONSE": "Hallmark IFNα",
        "INTERFERON_GAMMA_RESPONSE": "Hallmark IFNγ",
    }

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), gridspec_kw={"width_ratios": [1.45, 1]})

    ax = axes[0]
    x = np.arange(len(show_panels))
    width = 0.18
    for i, contrast in enumerate(primary):
        sub = panels[panels.contrast == contrast].set_index("panel")
        vals = [sub.at[p, "median_log2FC"] if p in sub.index else np.nan for p in show_panels]
        ax.bar(x + (i - 1.5) * width, vals, width=width, color=colors[contrast], label=labels[contrast], zorder=3)
    ax.axhline(0, color="#333", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([short[p] for p in show_panels], rotation=20, ha="right")
    ax.set_ylabel("Median log2 fold change\n(genes passing the count filter)")
    ax.set_title("IFN / STING / APM after NHEJ-component loss")
    ax.legend(frameon=False, ncol=2, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    # on-target log2FC for the perturbed gene
    ot = on_target.set_index("contrast").loc[primary]
    ypos = np.arange(len(primary))
    vals = ot["log2FoldChange"].to_numpy(dtype=float)
    ax.barh(ypos, vals, color=[colors[c] for c in primary], zorder=3)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{labels[c]}\n({ot.at[c, 'target']})" for c in primary])
    ax.set_xlabel("On-target log2 fold change")
    ax.set_title("Perturbation check")
    ax.set_xlim(min(-2.8, float(np.nanmin(vals)) - 0.15), 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for y, v, c in zip(ypos, vals, primary):
        padj = ot.at[c, "padj"]
        tag = f"{v:.2f}"
        if pd.notna(padj):
            tag += f"  padj {padj:.2g}"
        ax.text(0.06, y, tag, va="center", ha="left", fontsize=7.5, color="#222")

    fig.tight_layout()
    fig.savefig(OUT / "fig_nhej_ifn_panels.png", dpi=160)
    fig.savefig(OUT / "fig_nhej_ifn_panels.pdf")
    plt.close(fig)

    # focus heatmap of DESeq2 log2FC
    genes = [g for g in FOCUS_DISPLAY if g in set(focus["symbol"])]
    heat = (
        focus[focus.contrast.isin(primary)]
        .pivot(index="symbol", columns="contrast", values="log2FoldChange")
        .reindex(index=genes, columns=primary)
    )
    fig, ax = plt.subplots(figsize=(6.4, 8.2))
    data = heat.to_numpy(dtype=float)
    finite = data[np.isfinite(data)]
    lim = float(np.nanmax(np.abs(finite))) if finite.size else 1
    lim = max(1.0, min(lim, 3.0))
    im = ax.imshow(data, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(primary)))
    ax.set_xticklabels([labels[c] for c in primary], rotation=0)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            txt = "·" if not np.isfinite(val) else f"{val:.2f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.5, color="#111" if (not np.isfinite(val) or abs(val) < lim * 0.55) else "white")
    ax.set_title("log2FC, depleted vs matched control")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="log2FC")
    fig.tight_layout()
    fig.savefig(OUT / "fig_focus_log2fc.png", dpi=160)
    fig.savefig(OUT / "fig_focus_log2fc.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
