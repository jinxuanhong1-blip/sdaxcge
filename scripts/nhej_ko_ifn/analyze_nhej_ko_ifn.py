#!/usr/bin/env python3
"""Score interferon signatures in GEO RNA-seq of LIG4, XRCC4, PRKDC, and TP53BP1 knockouts.

Primary estimand, fixed before looking at effect sizes: in a cultured
human or mouse cancer cell line, the difference (KO minus matched control)
in the mean log2 expression of a prespecified cell-intrinsic type-I ISG
core. Hallmark interferon-alpha and interferon-gamma sets are secondary.
Tumor bulk, non-malignant lines, and drug-treated arms are scored and
labeled, and are not the primary estimand.

Matrices are processed GEO files (counts or FPKM). Raw FASTQs are not used.
"""

from __future__ import annotations

import gzip
import io
import json
import math
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "results" / "nhej_ko_ifn"
CACHE = OUT / "_cache"
TABLES = OUT / "tables"
HALLMARK_PATH = Path(__file__).resolve().parent / "hallmark_ifn.tsv"

# Cell-intrinsic type-I ISG core. Chosen as antiviral ISGs that epithelial
# and cancer lines actually express, not immune-cell markers.
HUMAN_ISG = [
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
    "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "IFI44", "IFI44L",
    "IFI6", "IFI27", "BST2", "USP18", "IRF7", "STAT1", "STAT2",
    "DDX58", "IFIH1", "XAF1", "HERC5", "CMPK2", "IFITM1",
    "EIF2AK2", "ISG20", "IFNB1",
]
# Mouse stand-ins: Oas1a for OAS1, Oasl2 for OASL, Herc6 for HERC5,
# Ifi27l2a for IFI27. IFI44L and IFI6 have no usable mouse ortholog.
MOUSE_ISG = [
    "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2",
    "Oas1a", "Oas2", "Oas3", "Oasl2", "Rsad2", "Ifi44",
    "Ifi27l2a", "Bst2", "Usp18", "Irf7", "Stat1", "Stat2",
    "Ddx58", "Ifih1", "Xaf1", "Herc6", "Cmpk2", "Ifitm1",
    "Eif2ak2", "Isg20", "Ifnb1",
]
MOUSE_ALIAS = {
    "IFI44L": None,
    "IFI6": None,
    "OAS1": "Oas1a",
    "OASL": "Oasl2",
    "HERC5": "Herc6",
    "IFI27": "Ifi27l2a",
}

# GSE84986 sample titles, from GEO esummary (stable GSM accessions).
MCF7_TITLE = {
    "GSM2255508": "Wild-type untreated replicate 1",
    "GSM2255520": "Wild-type untreated replicate 2",
    "GSM2255532": "Wild-type untreated replicate 3",
    "GSM2255509": "Wild-type IR replicate 1",
    "GSM2255521": "Wild-type IR replicate 2",
    "GSM2255533": "Wild-type IR replicate 3",
    "GSM2255510": "Wild-type Nutlin-3 replicate 1",
    "GSM2255522": "Wild-type Nutlin-3 replicate 2",
    "GSM2255534": "Wild-type Nutlin-3 replicate 3",
    "GSM2255514": "53BP1Δ-1 untreated replicate 1",
    "GSM2255526": "53BP1Δ-1 untreated replicate 2",
    "GSM2255538": "53BP1Δ-1 untreated replicate 3",
    "GSM2255515": "53BP1Δ-1 IR replicate 1",
    "GSM2255527": "53BP1Δ-1 IR replicate 2",
    "GSM2255539": "53BP1Δ-1 IR replicate 3",
    "GSM2255516": "53BP1Δ-1 Nutlin-3 replicate 1",
    "GSM2255528": "53BP1Δ-1 Nutlin-3 replicate 2",
    "GSM2255540": "53BP1Δ-1 Nutlin-3 replicate 3",
    "GSM2255517": "53BP1Δ-2 untreated replicate 1",
    "GSM2255529": "53BP1Δ-2 untreated replicate 2",
    "GSM2255541": "53BP1Δ-2 untreated replicate 3",
    "GSM2255518": "53BP1Δ-2 IR replicate 1",
    "GSM2255530": "53BP1Δ-2 IR replicate 2",
    "GSM2255542": "53BP1Δ-2 IR replicate 3",
    "GSM2255519": "53BP1Δ-2 Nutlin-3 replicate 1",
    "GSM2255531": "53BP1Δ-2 Nutlin-3 replicate 2",
    "GSM2255543": "53BP1Δ-2 Nutlin-3 replicate 3",
    "GSM2255511": "p53 untreated replicate 1",
    "GSM2255523": "p53 untreated replicate 2",
    "GSM2255535": "p53 untreated replicate 3",
    "GSM2255512": "p53 IR replicate 1",
    "GSM2255524": "p53 IR replicate 2",
    "GSM2255536": "p53 IR replicate 3",
    "GSM2255513": "p53 Nutlin-3 replicate 1",
    "GSM2255525": "p53 Nutlin-3 replicate 2",
    "GSM2255537": "p53 Nutlin-3 replicate 3",
}

# GEO Sample_description library ids for GSE285698.
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

DOWNLOADS = {
    "GSE135274_all_sample_human_RNA_max.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135274/suppl/GSE135274_all_sample_human_RNA_max.txt.gz",
    "GSE154443_gene.fpkm.matrix.xlsx":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154443/suppl/GSE154443_gene.fpkm.matrix.xlsx",
    "GSE285698_raw_counts.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285698/suppl/GSE285698_raw_counts.txt.gz",
    "GSE145148_count_matrix.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145148/suppl/GSE145148_count_matrix.csv.gz",
    "GSE237615_Raw_Counts_ID8.xls.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237615/suppl/GSE237615_Raw_Counts_ID8.xls.txt.gz",
    "GSE237615_Raw_Counts_ID8_53BP1_KO.xls.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237615/suppl/GSE237615_Raw_Counts_ID8_53BP1_KO.xls.txt.gz",
    "GSE237615_Raw_Counts_KPC.xls.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237615/suppl/GSE237615_Raw_Counts_KPC.xls.txt.gz",
    "GSE237615_Raw_Counts_KPC_53BP1_KO.xls.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237615/suppl/GSE237615_Raw_Counts_KPC_53BP1_KO.xls.txt.gz",
    "GSE280049_53BP1_KO_RNA-seq.xlsx":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE280nnn/GSE280049/suppl/GSE280049_53BP1_KO_RNA-seq.xlsx",
    "GSE84986_RAW.tar":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84986/suppl/GSE84986_RAW.tar",
}


def download(name: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / name
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = DOWNLOADS[name]
    print(f"download {name}")
    req = urllib.request.Request(url, headers={"User-Agent": "nhej-ko-ifn/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    return dest


def ensembl_map(symbols: list[str], species: str) -> dict[str, str]:
    """Map gene symbols to Ensembl gene ids (no version)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    cache = CACHE / f"ensembl_{species}.json"
    cached = json.loads(cache.read_text()) if cache.exists() else {}
    missing = [s for s in symbols if s not in cached]
    if missing:
        url = f"https://rest.ensembl.org/lookup/symbol/{species}"
        body = json.dumps({"symbols": missing}).encode()
        req = urllib.request.Request(
            url, data=body,
            headers={"Content-Type": "application/json", "Accept": "application/json",
                     "User-Agent": "nhej-ko-ifn/1.0"},
        )
        data = json.loads(urllib.request.urlopen(req, timeout=120).read().decode())
        for sym in missing:
            rec = data.get(sym)
            cached[sym] = rec.get("id") if isinstance(rec, dict) else None
        cache.write_text(json.dumps(cached, indent=1))
    return {s: cached.get(s) for s in symbols}


def collapse_symbols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.index = df.index.astype(str)
    df = df.groupby(level=0).sum(numeric_only=True)
    return df


def counts_to_log2(counts: pd.DataFrame) -> pd.DataFrame:
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts = collapse_symbols(counts)
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def abundance_to_log2(mat: pd.DataFrame) -> pd.DataFrame:
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mat = collapse_symbols(mat)
    return np.log2(mat + 1.0)


def load_hallmark() -> dict[str, list[str]]:
    df = pd.read_csv(HALLMARK_PATH, sep="\t")
    out = {}
    for sig, sub in df.groupby("signature"):
        out[sig] = sub["symbol"].tolist()
    return out


def to_mouse_symbol(sym: str) -> str | None:
    if sym in MOUSE_ALIAS:
        return MOUSE_ALIAS[sym]
    if not sym.isalnum():
        return None
    return sym[0].upper() + sym[1:].lower()


def mouse_hallmark(human_genes: list[str]) -> list[str]:
    out = []
    seen = set()
    for g in human_genes:
        m = to_mouse_symbol(g)
        if m and m not in seen:
            seen.add(m)
            out.append(m)
    return out


def welch(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Return delta (a-b), SE, two-sided Welch p. SE/p are nan if n<2."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    delta = float(a.mean() - b.mean())
    if len(a) < 2 or len(b) < 2:
        return delta, float("nan"), float("nan")
    se = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    if se == 0:
        p = 1.0 if delta == 0 else 0.0
        return delta, 0.0, p
    res = stats.ttest_ind(a, b, equal_var=False)
    return delta, float(se), float(res.pvalue)


def score_genes(logmat: pd.DataFrame, ctrl: list[str], ko: list[str], genes: list[str]) -> dict:
    present = [g for g in genes if g in logmat.index]
    cols = [c for c in ctrl + ko if c in logmat.columns]
    missing_samples = [c for c in ctrl + ko if c not in logmat.columns]
    if missing_samples:
        raise KeyError(f"missing samples {missing_samples}")
    sub = logmat.loc[present, cols]
    if sub.empty:
        detected = []
    else:
        detected = sub.index[sub.max(axis=1) >= 1.0].tolist()
    use = sub.loc[detected] if detected else sub.iloc[0:0]
    gene_lfc = {}
    if len(use):
        gene_lfc = (use[ko].mean(axis=1) - use[ctrl].mean(axis=1)).to_dict()
        scores = use.mean(axis=0)
        delta, se, p = welch(scores[ko].to_numpy(), scores[ctrl].to_numpy())
        sample_scores = scores.to_dict()
    else:
        delta, se, p = float("nan"), float("nan"), float("nan")
        sample_scores = {}
    n_up = sum(1 for v in gene_lfc.values() if v > 0)
    return {
        "n_genes_in_set": len(genes),
        "n_genes_present": len(present),
        "n_genes_detected": len(detected),
        "delta": delta,
        "se": se,
        "p": p,
        "frac_up": (n_up / len(gene_lfc)) if gene_lfc else float("nan"),
        "median_gene_lfc": float(np.median(list(gene_lfc.values()))) if gene_lfc else float("nan"),
        "gene_lfc": gene_lfc,
        "sample_scores": sample_scores,
    }


def target_stats(logmat: pd.DataFrame, ctrl: list[str], ko: list[str], symbol: str) -> dict:
    if symbol not in logmat.index:
        return {"target_ctrl_mean": float("nan"), "target_ko_mean": float("nan"),
                "target_log2_delta": float("nan"), "target_qc": "absent"}
    ctrl_mean = float(logmat.loc[symbol, ctrl].mean())
    ko_mean = float(logmat.loc[symbol, ko].mean())
    delta = ko_mean - ctrl_mean
    # log2(x+1) < 1 means the underlying CPM or FPKM is under 1.
    if ctrl_mean < 1.0 and ko_mean < 1.0:
        qc = "near_floor_both"
    elif delta <= -0.5:
        qc = "rna_down"
    elif delta < -0.25:
        qc = "rna_modestly_down"
    else:
        qc = "rna_not_reduced"
    return {"target_ctrl_mean": ctrl_mean, "target_ko_mean": ko_mean,
            "target_log2_delta": delta, "target_qc": qc}


# ---------------------------------------------------------------------------
# Loaders. Each returns a gene-symbol by sample log2 matrix.
# ---------------------------------------------------------------------------

def load_hela() -> pd.DataFrame:
    path = download("GSE135274_all_sample_human_RNA_max.txt.gz")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # Rows are repeated under mismatch classes (mis_0, mis_1, mis_2).
    # Keep exact matches only so mismatch alignments are not added on top.
    df = df[df["matchType"] == "human_RNA.mis_0"].copy()
    df["geneName"] = df["geneName"].astype(str)
    # Two columns per GEO sample are technical splits. Summing them does not
    # change CPM relative to averaging, because library size scales too.
    pairs = {
        "Exp01_XRCC4KO": ["Experiment_01_g2G3_1_1", "Experiment_01_g2G3_1_2"],
        "Exp01_XRCC4KO_mirin": ["Experiment_01_g2G3_M_1_1", "Experiment_01_g2G3_M_1_2"],
        "Exp01_CTRL": ["Experiment_01_gSCR_1_1", "Experiment_01_gSCR_1_2"],
        "Exp01_CTRL_mirin": ["Experiment_01_gSCR_M_1_1", "Experiment_01_gSCR_M_1_2"],
        "Exp02_XRCC4KO": ["Experiment_02_g2G3_2_1", "Experiment_02_g2G3_2_2"],
        "Exp02_XRCC4KO_mirin": ["Experiment_02_g2G3_M_2_1", "Experiment_02_g2G3_M_2_2"],
        "Exp02_CTRL": ["Experiment_02_gSCR_2_1", "Experiment_02_gSCR_2_2"],
        "Exp02_CTRL_mirin": ["Experiment_02_gSCR_M_2_1", "Experiment_02_gSCR_M_2_2"],
    }
    out = {}
    for name, cols in pairs.items():
        out[name] = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    mat = pd.DataFrame(out)
    mat.index = df["geneName"].values
    return counts_to_log2(mat)


def load_ch12() -> pd.DataFrame:
    path = download("GSE154443_gene.fpkm.matrix.xlsx")
    df = pd.read_excel(path, sheet_name=0)
    df = df.rename(columns={df.columns[0]: "ensembl", df.columns[1]: "CH12F3_WT", df.columns[2]: "CH12F3_Lig4KO"})
    symbols = [
        "Lig4", "Xrcc4", "Prkdc", "Trp53bp1", *MOUSE_ISG,
    ]
    # Also map hallmark later; caller passes genes. Here map a broad mouse set.
    hallmark = load_hallmark()
    extra = mouse_hallmark(hallmark["HALLMARK_IFN_ALPHA"] + hallmark["HALLMARK_IFN_GAMMA"])
    want = list(dict.fromkeys(symbols + extra))
    mapping = ensembl_map(want, "mus_musculus")
    id_to_sym = {}
    for sym, eid in mapping.items():
        if eid and eid not in id_to_sym:
            id_to_sym[eid] = sym
    df["symbol"] = df["ensembl"].map(id_to_sym)
    df = df.dropna(subset=["symbol"])
    mat = df.groupby("symbol")[["CH12F3_WT", "CH12F3_Lig4KO"]].sum()
    return abundance_to_log2(mat)


def load_hct116() -> pd.DataFrame:
    path = download("GSE285698_raw_counts.txt.gz")
    df = pd.read_csv(path, sep="\t")
    gene_col = df.columns[0]
    symbols = []
    for raw in df[gene_col].astype(str):
        symbols.append(raw.split("|")[-1] if "|" in raw else raw)
    df = df.drop(columns=[gene_col])
    df.index = symbols
    renamed = {}
    for col in df.columns:
        key = col
        # A1_S50_L007 -> A1_S50
        parts = col.split("_")
        if len(parts) >= 2:
            key = "_".join(parts[:2])
        if key not in HCT_LIB:
            raise KeyError(col)
        geno, cond, rep = HCT_LIB[key]
        renamed[col] = f"{geno}_{cond}_{rep}"
    df = df.rename(columns=renamed)
    return counts_to_log2(df)


def load_mcf10a() -> pd.DataFrame:
    path = download("GSE145148_count_matrix.csv.gz")
    df = pd.read_csv(path)
    # Cufflinks-style matrix: gene_short_name plus sample columns.
    meta_cols = [c for c in df.columns if c in ("tracking_id", "gene_id", "gene_short_name", "locus") or str(c).startswith("Unnamed")]
    samples = [c for c in df.columns if c not in meta_cols and c != ""]
    df = df[["gene_short_name", *samples]].copy()
    df["gene_short_name"] = df["gene_short_name"].astype(str)
    df = df.groupby("gene_short_name")[samples].sum()
    return counts_to_log2(df)


def _read_bgi_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    sym_col = "Gene Symbol" if "Gene Symbol" in df.columns else df.columns[1]
    df[sym_col] = (
        df[sym_col].astype(str).str.replace("'", "", regex=False).str.strip()
    )
    count_cols = [c for c in df.columns if c not in ("Gene ID", "Gene Symbol", "Type") and c != sym_col]
    # If sym_col was columns[1] and not named Gene Symbol, drop non-counts.
    count_cols = [c for c in df.columns if "Read Count" in c or c.endswith("Count")]
    if not count_cols:
        count_cols = [c for c in df.columns if c not in ("Gene ID", sym_col, "Type")]
    out = df[[sym_col, *count_cols]].groupby(sym_col)[count_cols].sum()
    return out


def load_tumors() -> pd.DataFrame:
    pieces = []
    files = {
        "GSE237615_Raw_Counts_ID8.xls.txt.gz": None,
        "GSE237615_Raw_Counts_ID8_53BP1_KO.xls.txt.gz": None,
        "GSE237615_Raw_Counts_KPC.xls.txt.gz": None,
        "GSE237615_Raw_Counts_KPC_53BP1_KO.xls.txt.gz": None,
    }
    frames = []
    for name in files:
        frames.append(_read_bgi_counts(download(name)))
    mat = frames[0]
    for fr in frames[1:]:
        mat = mat.join(fr, how="outer")
    mat = mat.fillna(0.0)
    # Friendly sample names.
    rename = {}
    for c in mat.columns:
        cl = c.replace(" Read Count", "")
        rename[c] = cl
    mat = mat.rename(columns=rename)
    return counts_to_log2(mat)


def load_rpe1() -> pd.DataFrame:
    path = download("GSE280049_53BP1_KO_RNA-seq.xlsx")
    df = pd.read_excel(path, sheet_name=0)
    df["SYMBOL"] = df["SYMBOL"].astype(str)
    cols = ["RPE-NC1", "RPE-NC2", "RPE-KO1", "RPE-KO2"]
    mat = df.groupby("SYMBOL")[cols].sum()
    return counts_to_log2(mat)


def load_mcf7() -> pd.DataFrame:
    path = download("GSE84986_RAW.tar")
    # Symbols needed: ISG core, hallmarks, targets.
    hallmark = load_hallmark()
    symbols = list(dict.fromkeys(
        HUMAN_ISG + hallmark["HALLMARK_IFN_ALPHA"] + hallmark["HALLMARK_IFN_GAMMA"]
        + ["TP53BP1", "LIG4", "XRCC4", "PRKDC", "TP53"]
    ))
    mapping = ensembl_map(symbols, "homo_sapiens")
    id_to_sym = {}
    for sym, eid in mapping.items():
        if eid and eid not in id_to_sym:
            id_to_sym[eid] = sym
    wanted = set(id_to_sym)
    columns = {}
    with tarfile.open(path) as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            base = Path(member.name).name
            gsm = base.split("_", 1)[0]
            if gsm not in MCF7_TITLE:
                continue
            raw = tar.extractfile(member).read()
            text = gzip.decompress(raw).decode("utf-8", "replace")
            vals = {}
            for line in text.splitlines():
                if not line or line.startswith("#"):
                    continue
                eid, _, rest = line.partition("\t")
                eid = eid.split(".")[0]
                if eid not in wanted:
                    continue
                try:
                    vals[eid] = vals.get(eid, 0.0) + float(rest.strip() or 0)
                except ValueError:
                    continue
            # Library size must use ALL genes, not the signature subset.
            lib = 0.0
            for line in text.splitlines():
                if not line or line.startswith("#"):
                    continue
                _eid, _, rest = line.partition("\t")
                try:
                    lib += float(rest.strip() or 0)
                except ValueError:
                    continue
            columns[gsm] = (vals, lib)
    # Build CPM for mapped symbols only, but library size is full.
    syms = sorted(set(id_to_sym.values()))
    data = {sym: {} for sym in syms}
    for gsm, (vals, lib) in columns.items():
        title = MCF7_TITLE[gsm]
        for eid, sym in id_to_sym.items():
            cpm = (vals.get(eid, 0.0) / lib * 1e6) if lib else 0.0
            data[sym][title] = data[sym].get(title, 0.0) + cpm
    mat = pd.DataFrame(data).T
    mat = collapse_symbols(mat)
    return np.log2(mat + 1.0)


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return ""
    return f"{x:.{nd}g}"


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    hallmark = load_hallmark()
    print("loading matrices")
    hela = load_hela()
    ch12 = load_ch12()
    hct = load_hct116()
    mcf10a = load_mcf10a()
    tumors = load_tumors()
    rpe = load_rpe1()
    mcf7 = load_mcf7()
    print("loaded", {k: v.shape for k, v in {
        "hela": hela, "ch12": ch12, "hct": hct, "mcf10a": mcf10a,
        "tumors": tumors, "rpe": rpe, "mcf7": mcf7,
    }.items()})

    def genes_for(organism: str, signature: str) -> list[str]:
        if signature == "ISG_CORE":
            return HUMAN_ISG if organism == "human" else MOUSE_ISG
        human = hallmark[signature]
        if organism == "human":
            return human
        return mouse_hallmark(human)

    # contrast definitions
    # primary_eligible: cultured cancer line, genetic KO, no drug / no hypoxia
    contrasts = [
        dict(id="GSE135274_HeLa_XRCC4_baseline", accession="GSE135274", gene="XRCC4",
             organism="human", cell_line="HeLa", disease="cervical adenocarcinoma",
             compartment="cultured_cell_line", condition="no mirin",
             eligible="yes", matrix=hela, scale="log2(CPM+1)",
             ctrl=["Exp01_CTRL", "Exp02_CTRL"], ko=["Exp01_XRCC4KO", "Exp02_XRCC4KO"],
             target="XRCC4",
             note="Two biological experiments. Technical column pairs were summed. Only exact-match rows (human_RNA.mis_0) were kept. g2G3 is the XRCC4(-/-) genotype in the GEO sample records."),
        dict(id="GSE135274_HeLa_XRCC4_mirin", accession="GSE135274", gene="XRCC4",
             organism="human", cell_line="HeLa", disease="cervical adenocarcinoma",
             compartment="cultured_cell_line", condition="mirin (MRE11 inhibitor)",
             eligible="no_drug", matrix=hela, scale="log2(CPM+1)",
             ctrl=["Exp01_CTRL_mirin", "Exp02_CTRL_mirin"],
             ko=["Exp01_XRCC4KO_mirin", "Exp02_XRCC4KO_mirin"],
             target="XRCC4",
             note="Same KO on a mirin background. Not the unstressed primary contrast."),
        dict(id="GSE154443_CH12F3_Lig4", accession="GSE154443", gene="LIG4",
             organism="mouse", cell_line="CH12F3", disease="mouse B-cell lymphoma",
             compartment="cultured_cell_line", condition="untreated",
             eligible="yes_unreplicated", matrix=ch12, scale="log2(FPKM+1)",
             ctrl=["CH12F3_WT"], ko=["CH12F3_Lig4KO"],
             target="Lig4",
             note="Deposited matrix is one FPKM column per genotype, so there is no within-study SE. Directional only."),
        dict(id="GSE285698_HCT116_PRKDC_normoxia", accession="GSE285698", gene="PRKDC",
             organism="human", cell_line="HCT116", disease="colorectal carcinoma",
             compartment="cultured_cell_line", condition="normoxia (vehicle)",
             eligible="yes", matrix=hct, scale="log2(CPM+1)",
             ctrl=["WT_normoxia_rep1", "WT_normoxia_rep2", "WT_normoxia_rep3"],
             ko=["KO_normoxia_rep1", "KO_normoxia_rep2", "KO_normoxia_rep3"],
             target="PRKDC",
             note="DNA-PK knockout is PRKDC. Hypoxia (CoCl2) is a separate contrast."),
        dict(id="GSE285698_HCT116_PRKDC_hypoxia", accession="GSE285698", gene="PRKDC",
             organism="human", cell_line="HCT116", disease="colorectal carcinoma",
             compartment="cultured_cell_line", condition="CoCl2 hypoxia",
             eligible="no_drug", matrix=hct, scale="log2(CPM+1)",
             ctrl=["WT_hypoxia_rep1", "WT_hypoxia_rep2", "WT_hypoxia_rep3"],
             ko=["KO_hypoxia_rep1", "KO_hypoxia_rep2", "KO_hypoxia_rep3"],
             target="PRKDC",
             note="KO versus WT under CoCl2. Not the unstressed primary contrast."),
        dict(id="GSE84986_MCF7_53BP1_untreated", accession="GSE84986", gene="TP53BP1",
             organism="human", cell_line="MCF-7", disease="breast carcinoma",
             compartment="cultured_cell_line", condition="untreated",
             eligible="yes", matrix=mcf7, scale="log2(CPM+1)",
             ctrl=["Wild-type untreated replicate 1", "Wild-type untreated replicate 2", "Wild-type untreated replicate 3"],
             ko=[
                 "53BP1Δ-1 untreated replicate 1", "53BP1Δ-1 untreated replicate 2", "53BP1Δ-1 untreated replicate 3",
                 "53BP1Δ-2 untreated replicate 1", "53BP1Δ-2 untreated replicate 2", "53BP1Δ-2 untreated replicate 3",
             ],
             target="TP53BP1",
             note="Two independent 53BP1-null clones pooled (n=6) versus WT (n=3). Clones share a parental line, so the Welch n is not six independent knockouts."),
        dict(id="GSE84986_MCF7_53BP1_clone1_untreated", accession="GSE84986", gene="TP53BP1",
             organism="human", cell_line="MCF-7", disease="breast carcinoma",
             compartment="cultured_cell_line", condition="untreated, clone 1",
             eligible="sensitivity", matrix=mcf7, scale="log2(CPM+1)",
             ctrl=["Wild-type untreated replicate 1", "Wild-type untreated replicate 2", "Wild-type untreated replicate 3"],
             ko=["53BP1Δ-1 untreated replicate 1", "53BP1Δ-1 untreated replicate 2", "53BP1Δ-1 untreated replicate 3"],
             target="TP53BP1",
             note="Clone-level sensitivity."),
        dict(id="GSE84986_MCF7_53BP1_clone2_untreated", accession="GSE84986", gene="TP53BP1",
             organism="human", cell_line="MCF-7", disease="breast carcinoma",
             compartment="cultured_cell_line", condition="untreated, clone 2",
             eligible="sensitivity", matrix=mcf7, scale="log2(CPM+1)",
             ctrl=["Wild-type untreated replicate 1", "Wild-type untreated replicate 2", "Wild-type untreated replicate 3"],
             ko=["53BP1Δ-2 untreated replicate 1", "53BP1Δ-2 untreated replicate 2", "53BP1Δ-2 untreated replicate 3"],
             target="TP53BP1",
             note="Clone-level sensitivity."),
        dict(id="GSE84986_MCF7_53BP1_IR", accession="GSE84986", gene="TP53BP1",
             organism="human", cell_line="MCF-7", disease="breast carcinoma",
             compartment="cultured_cell_line", condition="5 Gy IR, 4 h",
             eligible="no_baseline", matrix=mcf7, scale="log2(CPM+1)",
             ctrl=["Wild-type IR replicate 1", "Wild-type IR replicate 2", "Wild-type IR replicate 3"],
             ko=[
                 "53BP1Δ-1 IR replicate 1", "53BP1Δ-1 IR replicate 2", "53BP1Δ-1 IR replicate 3",
                 "53BP1Δ-2 IR replicate 1", "53BP1Δ-2 IR replicate 2", "53BP1Δ-2 IR replicate 3",
             ],
             target="TP53BP1",
             note="Irradiated arm. Not the baseline primary contrast."),
        dict(id="GSE145148_MCF10A_XRCC4_vs_WT", accession="GSE145148", gene="XRCC4",
             organism="human", cell_line="MCF10A", disease="non-tumorigenic mammary epithelium",
             compartment="cultured_cell_line", condition="unstressed",
             eligible="no_not_cancer", matrix=mcf10a, scale="log2(CPM+1)",
             ctrl=["MCF10A_WT_S1", "MCF10A_WT_S2"],
             ko=["MCF10A_XRCC4KO_S1", "MCF10A_XRCC4KO_S2"],
             target="XRCC4",
             note="MCF10A is not a cancer line. Scored because it is the clearest public XRCC4-KO RNA-seq and is the biological context for cGAS activation after NHEJ loss."),
        dict(id="GSE145148_MCF10A_XRCC4_p53DKO_vs_p53KO", accession="GSE145148", gene="XRCC4",
             organism="human", cell_line="MCF10A", disease="non-tumorigenic mammary epithelium",
             compartment="cultured_cell_line", condition="TP53-null background",
             eligible="no_not_cancer", matrix=mcf10a, scale="log2(CPM+1)",
             ctrl=["MCF10A_p53KO_S1", "MCF10A_p53KO_S2"],
             ko=["MCF10A_XRCC4p53DKO_S1", "MCF10A_XRCC4p53DKO_S2"],
             target="XRCC4",
             note="XRCC4/TP53 double KO versus TP53 KO, so the contrast is XRCC4 loss on a p53-null background."),
        dict(id="GSE280049_RPE1_53BP1", accession="GSE280049", gene="TP53BP1",
             organism="human", cell_line="hTERT-RPE1", disease="non-malignant retinal pigment epithelium",
             compartment="cultured_cell_line", condition="unstressed",
             eligible="no_not_cancer", matrix=rpe, scale="log2(CPM+1)",
             ctrl=["RPE-NC1", "RPE-NC2"], ko=["RPE-KO1", "RPE-KO2"],
             target="TP53BP1",
             note="CRISPR sg53BP1 in hTERT-RPE1. Not a cancer line. Deposited values are integer counts."),
        dict(id="GSE237615_ID8_tumor_53BP1", accession="GSE237615", gene="TP53BP1",
             organism="mouse", cell_line="ID8", disease="ovarian cancer allograft",
             compartment="allograft_tumor", condition="untreated tumor",
             eligible="no_tumor_bulk", matrix=tumors, scale="log2(CPM+1)",
             ctrl=["ID8_Ctrl1", "ID8_Ctrl2", "ID8_Ctrl3"],
             ko=["ID8_53BP1_KO1", "ID8_53BP1_KO2", "ID8_53BP1_KO3"],
             target="Trp53bp1",
             note="53BP1 was knocked out in the ID8 cancer line, but RNA-seq is bulk tumor. Interferon genes can come from infiltrating immune cells."),
        dict(id="GSE237615_KPC_tumor_53BP1", accession="GSE237615", gene="TP53BP1",
             organism="mouse", cell_line="KPC", disease="pancreatic cancer allograft",
             compartment="allograft_tumor", condition="untreated tumor",
             eligible="no_tumor_bulk", matrix=tumors, scale="log2(CPM+1)",
             ctrl=["Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4", "Ctrl_5"],
             ko=["53BP1_KO_1", "53BP1_KO_2", "53BP1_KO_3", "53BP1_KO_4", "53BP1_KO_5"],
             target="Trp53bp1",
             note="Bulk KPC allograft tumors. Same immune-infiltrate caveat as ID8. Control and KO count files were joined on gene symbol."),
    ]

    signatures = ["ISG_CORE", "HALLMARK_IFN_ALPHA", "HALLMARK_IFN_GAMMA"]
    meta_rows = []
    sig_rows = []
    gene_rows = []
    sample_rows = []

    for c in contrasts:
        logmat = c["matrix"]
        tgt = target_stats(logmat, c["ctrl"], c["ko"], c["target"])
        base = {
            "contrast_id": c["id"],
            "accession": c["accession"],
            "gene_ko": c["gene"],
            "organism": c["organism"],
            "cell_line": c["cell_line"],
            "disease": c["disease"],
            "compartment": c["compartment"],
            "condition": c["condition"],
            "eligible": c["eligible"],
            "n_ctrl": len(c["ctrl"]),
            "n_ko": len(c["ko"]),
            "scale": c["scale"],
            "target_symbol": c["target"],
            "target_ctrl_mean": tgt["target_ctrl_mean"],
            "target_ko_mean": tgt["target_ko_mean"],
            "target_log2_delta": tgt["target_log2_delta"],
            "target_qc": tgt["target_qc"],
            "note": c["note"],
        }
        wide = dict(base)
        for sig in signatures:
            genes = genes_for(c["organism"], sig)
            sc = score_genes(logmat, c["ctrl"], c["ko"], genes)
            sig_rows.append({
                **base,
                "signature": sig,
                "n_genes_in_set": sc["n_genes_in_set"],
                "n_genes_present": sc["n_genes_present"],
                "n_genes_detected": sc["n_genes_detected"],
                "delta_log2": sc["delta"],
                "se": sc["se"],
                "welch_p": sc["p"],
                "frac_genes_up": sc["frac_up"],
                "median_gene_log2fc": sc["median_gene_lfc"],
            })
            prefix = {"ISG_CORE": "isg", "HALLMARK_IFN_ALPHA": "ifna", "HALLMARK_IFN_GAMMA": "ifng"}[sig]
            wide[f"{prefix}_n"] = sc["n_genes_detected"]
            wide[f"{prefix}_delta"] = sc["delta"]
            wide[f"{prefix}_se"] = sc["se"]
            wide[f"{prefix}_p"] = sc["p"]
            wide[f"{prefix}_frac_up"] = sc["frac_up"]
            wide[f"{prefix}_median_lfc"] = sc["median_gene_lfc"]
            if sig == "ISG_CORE":
                for g, lfc in sc["gene_lfc"].items():
                    gene_rows.append({
                        "contrast_id": c["id"], "accession": c["accession"], "gene_ko": c["gene"],
                        "signature_gene": g, "log2fc": lfc,
                    })
                for sample, val in sc["sample_scores"].items():
                    arm = "KO" if sample in c["ko"] else "CTRL"
                    sample_rows.append({
                        "contrast_id": c["id"], "sample": sample, "arm": arm, "isg_score": val,
                    })
        meta_rows.append(wide)
        print(
            f"{c['id']}: ISG {wide['isg_delta']:.3f} p={wide['isg_p']} "
            f"target {tgt['target_log2_delta']:.3f} ({tgt['target_qc']}) n={wide['isg_n']}"
        )

    meta = pd.DataFrame(meta_rows)
    sig = pd.DataFrame(sig_rows)
    genes = pd.DataFrame(gene_rows)
    samples = pd.DataFrame(sample_rows)
    meta.to_csv(TABLES / "meta_table.tsv", sep="\t", index=False)
    sig.to_csv(TABLES / "signature_scores.tsv", sep="\t", index=False)
    genes.to_csv(TABLES / "isg_core_gene_log2fc.tsv", sep="\t", index=False)
    samples.to_csv(TABLES / "isg_core_sample_scores.tsv", sep="\t", index=False)

    # One primary row per gene among eligible cultured-line contrasts.
    # CH12 is unreplicated, so it stays in the table but is not pooled.
    primary = meta[meta["eligible"].isin(["yes", "yes_unreplicated"])].copy()
    primary.to_csv(TABLES / "primary_cancer_line_contrasts.tsv", sep="\t", index=False)
    print("wrote", TABLES)


if __name__ == "__main__":
    main()
