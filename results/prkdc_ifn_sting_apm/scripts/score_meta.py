#!/usr/bin/env python3
"""Score IFN, STING, and antigen-presentation signatures on public GEO
DNA-PKcs (PRKDC) inhibitor or knockdown expression profiles, then
random-effects meta-analyse the epithelial-cancer studies.

Expression matrices are read from $PRKDC_DATA (default /tmp/geo/data).
They are the GEO supplementary files listed in scripts/download.sh.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import sqlite3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path(os.environ.get("PRKDC_DATA", "/tmp/geo/data"))
ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "tables"
FIG = ROOT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# Inducible cGAS-STING transcriptional readout. Machinery mRNA (CGAS, STING1,
# TBK1, IRF3) is reported gene-wise and is not the STING score: those genes
# are mostly regulated post-translationally.
STING_RESPONSE = [
    "IFNB1", "IFNL1", "IFNL2", "IFNL3", "CCL4", "CCL5", "CXCL9", "CXCL10",
    "CXCL11", "IL6", "TNF", "TNFAIP3", "NFKBIA", "IRF7", "ZBP1", "IFI16",
]
APM = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M",
    "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2",
    "ERAP1", "ERAP2", "NLRC5", "CALR", "CANX", "PDIA3",
]
# Removed from the IFN scores so MHC/APM genes are not counted twice.
APM_OVERLAP = set(APM) | {
    "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1",
    "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "CIITA", "CD74",
}
KEY_GENES = [
    "PRKDC", "IFNB1", "CXCL10", "CCL5", "STAT1", "IRF7", "ISG15", "MX1",
    "OAS1", "IFI27", "STING1", "CGAS", "TBK1", "IRF3", "HLA-A", "HLA-B",
    "B2M", "TAP1", "PSMB8", "PSMB9", "NLRC5", "CD274",
]
ALIASES = {
    "TMEM173": "STING1",
    "MB21D1": "CGAS",
    "C6orf150": "CGAS",
    "C6ORF150": "CGAS",
}


def load_hallmark():
    gmt = DATA / "h.all.Hs.symbols.gmt"
    sets = {}
    for line in gmt.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        sets[parts[0]] = parts[2:]
    ifna = [g for g in sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"] if g not in APM_OVERLAP]
    ifng = [g for g in sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"] if g not in APM_OVERLAP]
    union = list(dict.fromkeys(ifna + ifng))
    mhci = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
    return {
        "IFN_alpha": ifna, "IFN_gamma": ifng, "IFN": union,
        "STING": STING_RESPONSE, "APM": APM, "MHCI": mhci,
    }


def entrez_to_symbol():
    path = DATA / "Homo_sapiens.gene_info.gz"
    out = {}
    ens = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id = header.index("GeneID")
        i_sym = header.index("Symbol")
        i_xref = header.index("dbXrefs")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            out[p[i_id]] = p[i_sym]
            for x in p[i_xref].split("|"):
                if x.startswith("Ensembl:ENSG"):
                    ens[x.split(":", 1)[1]] = p[i_sym]
    return out, ens


def hugene_probe_map(entrez_sym):
    db = DATA / "hugene20sttranscriptcluster.db/inst/extdata/hugene20sttranscriptcluster.sqlite"
    con = sqlite3.connect(db)
    # Unambiguous transcript clusters only.
    rows = con.execute(
        "SELECT probe_id, gene_id FROM probes WHERE is_multiple = 0 AND gene_id IS NOT NULL"
    ).fetchall()
    mp = {}
    for probe, gid in rows:
        sym = entrez_sym.get(str(gid))
        if sym:
            mp[str(probe)] = sym
    return mp


def illumina_probe_map():
    path = DATA / "GPL10558.annot.gz"
    mp = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                break
        if header is None:
            raise RuntimeError("GPL10558 header missing")
        i_id = header.index("ID")
        i_sym = header.index("Gene symbol")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= i_sym:
                continue
            sym = p[i_sym].strip()
            if sym and sym not in {"NA", "---"}:
                mp[p[i_id]] = sym.split("///")[0].strip()
    return mp


def canon_index(df: pd.DataFrame) -> pd.DataFrame:
    idx = pd.Index([ALIASES.get(str(i), str(i)) for i in df.index])
    df = df.copy()
    df.index = idx
    df = df.groupby(level=0).sum(min_count=1)
    return df


def library_cpm(df: pd.DataFrame) -> pd.DataFrame:
    lib = df.sum(axis=0).replace(0, np.nan)
    return df.div(lib, axis=1) * 1e6


def to_log(df: pd.DataFrame, kind: str) -> pd.DataFrame:
    if kind == "counts":
        x = library_cpm(df)
        return np.log2(x + 1.0)
    if kind == "norm":
        return np.log2(df.clip(lower=0) + 1.0)
    if kind in {"log", "log_signed"}:
        return df
    raise ValueError(kind)


def detected(log_df, linear_df, cols, kind):
    if kind == "log_signed":
        return pd.Series(True, index=log_df.index)
    lin = linear_df[cols]
    # Present in at least half the samples of either arm is applied per contrast.
    return lin


def score_contrast(log_df, linear_df, treat, ctrl, genes, kind):
    treat = [c for c in treat if c in log_df.columns]
    ctrl = [c for c in ctrl if c in log_df.columns]
    genes_present = [g for g in genes if g in log_df.index]
    if len(treat) < 2 or len(ctrl) < 2 or len(genes_present) < 3:
        return None
    sub_lin = linear_df.loc[genes_present, treat + ctrl]
    if kind == "log_signed":
        keep = genes_present
    else:
        def arm_detected(cols):
            n_need = max(1, int(math.ceil(len(cols) / 2)))
            return (sub_lin[cols] > 0).sum(axis=1) >= n_need
        mask = arm_detected(treat) | arm_detected(ctrl)
        keep = [g for g in genes_present if bool(mask.loc[g])]
    if len(keep) < 3:
        return None
    sub = log_df.loc[keep, treat + ctrl].astype(float)
    st = sub[treat].mean(axis=0)
    sc = sub[ctrl].mean(axis=0)
    lfc = float(st.mean() - sc.mean())
    n1, n2 = len(st), len(sc)
    v1 = float(st.var(ddof=1))
    v2 = float(sc.var(ddof=1))
    se = math.sqrt(v1 / n1 + v2 / n2)
    if not math.isfinite(se) or se < 1e-4:
        se = 1e-4
    t_res = stats.ttest_ind(st.values, sc.values, equal_var=False)
    p = float(t_res.pvalue) if t_res.pvalue == t_res.pvalue else float("nan")
    sp2_num = (n1 - 1) * v1 + (n2 - 1) * v2
    sp2_den = n1 + n2 - 2
    sp = math.sqrt(sp2_num / sp2_den) if sp2_den > 0 and sp2_num > 0 else 0.0
    d = 0.0 if sp == 0 else (float(st.mean() - sc.mean()) / sp)
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    g = J * d
    vd = (n1 + n2) / (n1 * n2) + (d ** 2) / (2 * (n1 + n2))
    vg = (J ** 2) * vd
    gene_lfc_s = (sub[treat].mean(axis=1) - sub[ctrl].mean(axis=1)).astype(float)
    return {
        "lfc": lfc,
        "se": se,
        "p": p,
        "g": g,
        "g_var": vg,
        "lfc_w": float(gene_lfc_s.clip(-2, 2).mean()),
        "lfc_med": float(gene_lfc_s.median()),
        "n_treat": n1,
        "n_ctrl": n2,
        "n_genes": len(keep),
        "n_genes_in_set": len(genes),
        "n_genes_present": len(genes_present),
        "coverage": len(keep) / len(genes),
        "gene_lfc": gene_lfc_s,
        "frac_genes_up": float((gene_lfc_s > 0).mean()),
    }


def _global_median(log_df, linear_df, treat, ctrl, kind):
    treat = [c for c in treat if c in log_df.columns]
    ctrl = [c for c in ctrl if c in log_df.columns]
    if len(treat) < 2 or len(ctrl) < 2:
        return float("nan")
    if kind == "log_signed":
        use = log_df
    else:
        lin = linear_df[treat + ctrl]
        n_need_t = max(1, int(math.ceil(len(treat) / 2)))
        n_need_c = max(1, int(math.ceil(len(ctrl) / 2)))
        mask = ((lin[treat] > 0).sum(axis=1) >= n_need_t) | ((lin[ctrl] > 0).sum(axis=1) >= n_need_c)
        use = log_df.loc[mask]
    if use.empty:
        return float("nan")
    lfc = use[treat].mean(axis=1) - use[ctrl].mean(axis=1)
    return float(np.nanmedian(lfc.to_numpy(dtype=float)))


def gene_lfc(log_df, linear_df, treat, ctrl, gene, kind):
    if gene not in log_df.index:
        return float("nan")
    st = log_df.loc[gene, [c for c in treat if c in log_df.columns]]
    sc = log_df.loc[gene, [c for c in ctrl if c in log_df.columns]]
    st = pd.to_numeric(st, errors="coerce")
    sc = pd.to_numeric(sc, errors="coerce")
    if getattr(st, "ndim", 1) == 0 or len(st) == 0 or len(sc) == 0:
        return float("nan")
    return float(np.nanmean(st.to_numpy(dtype=float)) - np.nanmean(sc.to_numpy(dtype=float)))


def dl_random(effects, variances):
    effects = np.asarray(effects, dtype=float)
    variances = np.asarray(variances, dtype=float)
    k = len(effects)
    if k == 0:
        return None
    if k == 1:
        se = math.sqrt(float(variances[0]))
        mu = float(effects[0])
        z = mu / se if se > 0 else float("nan")
        p = float(2 * stats.norm.sf(abs(z))) if se > 0 else float("nan")
        return {
            "k": 1, "mu": mu, "se": se, "lo": mu - 1.96 * se, "hi": mu + 1.96 * se,
            "p": p, "I2": 0.0, "tau2": 0.0, "Q": 0.0, "n_pos": int(mu > 0),
        }
    w = 1.0 / variances
    mu_fe = float(np.sum(w * effects) / np.sum(w))
    Q = float(np.sum(w * (effects - mu_fe) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - (k - 1)) / c) if c > 0 else 0.0
    w_re = 1.0 / (variances + tau2)
    mu = float(np.sum(w_re * effects) / np.sum(w_re))
    se = math.sqrt(1.0 / float(np.sum(w_re)))
    z = mu / se if se > 0 else float("nan")
    p = float(2 * stats.norm.sf(abs(z))) if se > 0 else float("nan")
    I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    return {
        "k": k, "mu": mu, "se": se, "lo": mu - 1.96 * se, "hi": mu + 1.96 * se,
        "p": p, "I2": I2, "tau2": tau2, "Q": Q, "n_pos": int(np.sum(effects > 0)),
    }


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(np.where(np.isnan(p), 1.0, p))
    q = np.empty(m, dtype=float)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        raw_rank = m - rank + 1
        if np.isnan(p[i]):
            q[i] = float("nan")
            continue
        val = min(prev, p[i] * m / raw_rank)
        prev = val
        q[i] = val
    return q


# ---------------------------------------------------------------------------
# Loaders. Each returns (log_df, linear_df, kind) indexed by gene symbol.
# ---------------------------------------------------------------------------

def _prep(df, kind):
    df = canon_index(df.apply(pd.to_numeric, errors="coerce").fillna(0.0))
    if kind == "counts":
        linear = library_cpm(df)
        log_df = np.log2(linear + 1.0)
        return log_df, linear, kind
    if kind == "norm":
        linear = df.clip(lower=0)
        return np.log2(linear + 1.0), linear, kind
    if kind == "log_signed":
        return df, df, kind
    raise ValueError(kind)


def load_gse273409():
    df = pd.read_csv(DATA / "GSE273409_salmon_gene_counts_normalized.csv.gz")
    df = df.drop(columns=["gene_id"]).set_index("gene_name")
    return _prep(df, "norm")


def load_gse242255(ens_map):
    df = pd.read_csv(DATA / "GSE242255_AREK_DNAPK_DESeq2_normalised_counts.csv.gz")
    df["symbol"] = df["Ensembl_ID"].str.replace(r"\.\d+$", "", regex=True).map(ens_map)
    df = df.dropna(subset=["symbol"]).drop(columns=["Ensembl_ID"]).set_index("symbol")
    return _prep(df, "norm")


def load_featurecounts(path, ens_map, versioned=False):
    df = pd.read_csv(path, sep="\t", comment="#")
    # featureCounts files have annotation columns
    ann = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
    cols = [c for c in df.columns if c not in ann]
    sym = df["Geneid"].astype(str)
    if versioned:
        sym = sym.str.replace(r"\.\d+$", "", regex=True)
    df = df[cols]
    df.index = sym.map(lambda x: ens_map.get(x, x if not str(x).startswith("ENSG") else None))
    df = df[df.index.notna()]
    df.columns = [c.split("Aligned")[0].rstrip("_") for c in df.columns]
    # also strip path prefixes
    df.columns = [Path(c).name.replace(".bam", "") for c in df.columns]
    return _prep(df, "counts")


def load_gse116765(ens_map):
    log_df, linear, kind = load_featurecounts(
        DATA / "GSE116765_DNAPK_AllTreats_countTable.txt.gz", ens_map, versioned=False
    )
    # column cleanup happened inside; reload names more carefully
    return log_df, linear, kind


def load_gse116765_22rv1(ens_map):
    return load_featurecounts(
        DATA / "GSE116765_22Rv1_DNAPKi_CountTable.txt.gz", ens_map, versioned=False
    )


def load_gse167956():
    df = pd.read_excel(DATA / "GSE167956_RAW_COUNTS-gene_expressions_in_24_samples.xlsx")
    df = df.drop(columns=["FEATURE_ID"])
    df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
    df = df.set_index("GENE_SYMBOL")
    return _prep(df, "counts")


def load_gse180581(ens_map):
    df = pd.read_excel(DATA / "GSE180581_all_samples_counts.xlsx")
    df["symbol"] = df["name"].map(ens_map)
    df = df.dropna(subset=["symbol"]).drop(columns=["name"]).set_index("symbol")
    return _prep(df, "counts")


def load_gse268415():
    df = pd.read_csv(DATA / "GSE268415_gene_count.txt.gz", sep="\t")
    df = df.set_index("gene_name")
    keep = [c for c in df.columns if c.startswith("hos_")]
    return _prep(df[keep], "counts")


def load_gse285698(ens_map):
    df = pd.read_csv(DATA / "GSE285698_raw_counts.txt.gz", sep="\t")
    # index like ENSG...|SYMBOL
    raw = df.iloc[:, 0].astype(str)
    if df.columns[0] != "gene" and not str(df.columns[0]).startswith("ENSG"):
        # first column may be the index name that is blank -> read as index
        pass
    if df.index.name is None and not str(df.index[0]).startswith("ENSG"):
        df = pd.read_csv(DATA / "GSE285698_raw_counts.txt.gz", sep="\t", index_col=0)
    else:
        df = df.set_index(df.columns[0])
    symbols = []
    for ix in df.index.astype(str):
        if "|" in ix:
            symbols.append(ix.split("|", 1)[1])
        else:
            symbols.append(ens_map.get(ix.split(".")[0], None))
    df.index = symbols
    df = df[df.index.notna()]
    df.columns = [c.split("_L")[0] for c in df.columns]
    return _prep(df, "counts")


def load_gse129436(ens_map):
    df = pd.read_csv(DATA / "GSE129436_norm_matrix.csv.gz", index_col=0)
    df.index = [ens_map.get(str(i).split(".")[0]) for i in df.index]
    df = df[pd.notna(df.index)]
    return _prep(df, "log_signed")


def load_gse315862():
    df = pd.read_excel(DATA / "GSE315862_processed_counts.xlsx")
    df = df.rename(columns={df.columns[0]: "symbol"}).set_index("symbol")
    return _prep(df, "counts")


def load_gse319513(which):
    fn = {
        "24": "GSE319513_KU_24HVSCON_24H_Gene_differential_expression.xlsx",
        "48": "GSE319513_KU_48HVSCON_48H_Gene_differential_expression.xlsx",
    }[which]
    df = pd.read_excel(DATA / fn)
    fpkm = [c for c in df.columns if str(c).startswith("FPKM.")]
    out = df.set_index("gene_name")[fpkm]
    return _prep(out, "norm")


def load_illumina(path, probe_map, skip_comment=False):
    # Some GEO matrices prefix the table with # comments and a quoted note.
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        lines = []
        started = False
        for line in fh:
            if not started:
                if line.startswith("ID_REF") or line.startswith('"ID_REF"'):
                    started = True
                    lines.append(line)
                continue
            lines.append(line)
    from io import StringIO
    df = pd.read_csv(StringIO("".join(lines)), sep="\t")
    df = df.dropna(subset=["ID_REF"])
    symbols = df["ID_REF"].map(probe_map)
    sig_cols = [c for c in df.columns if c != "ID_REF" and "Detection" not in c and not str(c).startswith("Detection")]
    # paired detection columns follow each signal in the raw GEO layout, but
    # pandas names them uniquely. Match p-values by order.
    p_cols = [c for c in df.columns if "Detection" in str(c) or str(c).startswith("Detection")]
    mat = df[sig_cols].apply(pd.to_numeric, errors="coerce")
    mat.index = symbols
    if len(p_cols) == len(sig_cols):
        pmat = df[p_cols].apply(pd.to_numeric, errors="coerce")
        pmat.index = symbols
        # zero-out probes undetected (p >= 0.05) so they do not enter the mean
        undetected = pmat.ge(0.05).to_numpy()
        vals = np.array(mat.to_numpy(dtype=float), copy=True)
        vals[undetected] = np.nan
        mat = pd.DataFrame(vals, index=mat.index, columns=sig_cols)
    mat = mat[mat.index.notna()]
    # mean across probes of the same gene, ignoring undetected probes
    mat = mat.groupby(level=0).mean()
    mat = mat.fillna(0.0)
    return _prep(mat, "norm")


def load_gse63480(probe_map):
    path_url_local = DATA / "GSE63480_series_matrix.txt.gz"
    if not path_url_local.exists():
        raise FileNotFoundError(path_url_local)
    with gzip.open(path_url_local, "rt") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
        df = pd.read_csv(fh, sep="\t")
    df = df.dropna(subset=["ID_REF"])
    df = df[df["ID_REF"] != "!series_matrix_table_end"]
    idref = df["ID_REF"].astype(str).str.replace('"', "", regex=False)
    df = df.drop(columns=["ID_REF"])
    df.columns = [c.replace('"', "") for c in df.columns]
    df.index = idref.map(probe_map)
    df = df[df.index.notna()]
    df = df.apply(pd.to_numeric, errors="coerce")
    # GEO matrix for this array is RMA-like log2 intensity (median well below 20).
    # Confirm and treat as log. If the median is high, log it.
    med = float(np.nanmedian(df.values))
    df = df.groupby(level=0).mean()
    if med > 20:
        return _prep(df, "norm")
    return _prep(df, "log")


def screen_table():
    rows = [
        ("GSE273409", "include", "epithelial cancer", "SCLC lines DMS114 H146 H196 H446 H720 H1341, NU7441 vs DMSO, 5 d, n=3"),
        ("GSE242255", "include", "epithelial cancer", "CWR22Rv1-AR-EK siDNA-PKcs and NU7441/NU5455/AZD7648, n=3"),
        ("GSE116765", "include", "epithelial cancer", "C4-2 NU7441 vs control RNA-seq; CC-115 dual; CC-223 TORKi specificity"),
        ("GSE116765", "include", "epithelial cancer", "22Rv1 CC-115 vs control, dual DNA-PK/mTOR, n=2 vs 3"),
        ("GSE167956", "include", "epithelial cancer", "SUM159 shDNA-PKcs vs shNTC, normoxia primary, hypoxia companion"),
        ("GSE319513", "include", "epithelial cancer", "8505C anaplastic thyroid, KU-57788 vs CON, 24 h primary, 48 h companion. GEO characteristics text is copy-pasted onto both arms; column names KU vs CON are used"),
        ("GSE93075", "include", "epithelial cancer", "HCT116 siPRKDC vs siNeg, Illumina HT-12, n=3"),
        ("GSE93074", "include", "epithelial cancer", "HCT116 NU7441 40 µM 48 h vs vehicle, same paper as GSE93075"),
        ("GSE63480", "include", "epithelial cancer", "C4-2 DNA-PKcs siRNA and NU7441, HuGene 2.0 ST, n=2"),
        ("GSE285698", "include", "epithelial cancer", "HCT116 DNA-PK KO vs WT, normoxia primary, CoCl2 companion"),
        ("GSE85202", "include", "non-epithelial cancer", "Ewing TC32 and SK-N-MC NU7441 40 µM 48 h"),
        ("GSE93075", "include", "non-epithelial cancer", "TC32 Ewing siPRKDC"),
        ("GSE268415", "include", "non-epithelial cancer", "HOS osteosarcoma siPRKDC and AZD7648"),
        ("GSE315862", "include", "non-epithelial cancer", "LPS853 and T778 liposarcoma peposertib 0.5 µM vs NTC, no doxorubicin"),
        ("GSE129436", "include", "non-epithelial cancer", "U937 histiocytic lymphoma, NU7441 vs DMSO, lipofectamine only, 8 h. DNA-transfection arms scored separately"),
        ("GSE180581", "include", "epithelial non-cancer", "HEK293T DNA-PKcs heterozygous + siDNA-PKcs vs siControl"),
        ("GSE287819", "not scored", "epithelial cancer", "C4-2 AZD7648 4 h SLAM-seq, n=2. GEO deposits RAW only, no gene-count matrix"),
        ("GSE141686", "exclude", "not epithelial/cancer", "MC3T3-E1 osteoblasts + NU7441"),
        ("GSE225461", "exclude", "not epithelial/cancer", "primary naive B cells + DNA-PKcs inhibitor"),
        ("GSE183507", "exclude", "wrong drug", "AT2 organoids treated with p38 inhibitors VX-702 and PH-797804, not a DNA-PKcs inhibitor"),
        ("GSE294709", "exclude", "not PRKDC", "Ku70/Ku80 depletion, not DNA-PKcs"),
        ("GSE302992", "exclude", "not PRKDC", "PTK7 knockout in SUM159"),
        ("GSE285638", "exclude", "not PRKDC", "TIMELESS degron plus ATR inhibitor in HCT116"),
        ("GSE303774", "exclude", "not PRKDC", "CHK2 knockout"),
        ("GSE305143", "exclude", "not inhibitor/KD", "Licochalcone A reported to induce PRKDC"),
        ("GSE306599", "exclude", "not PRKDC", "CRTC2 auxin-degron RNA-seq"),
        ("GSE338690", "exclude", "not selective / not epithelial", "3-hydroxyflavone in AML lines"),
        ("GSE243758", "exclude", "not RNA-seq expression", "CRISPR screen, not expression RNA-seq"),
        ("GSE270287", "exclude", "not RNA-seq expression", "CUT&Tag"),
        ("GSE106609", "exclude", "not PRKDC", "PRMT1 knockdown in SK-OV-3; DNA-PK is mentioned as a pathway"),
        ("GSE180443", "exclude", "not PRKDC", "lnc-LEMGC overexpression"),
        ("GSE240510", "exclude", "not a cell-intrinsic PRKDC contrast", "mouse bone-marrow scRNA after AsiDNA"),
        ("GSE265857", "exclude", "not PRKDC", "Prrx1 knockout MEFs"),
    ]
    return pd.DataFrame(rows, columns=["accession", "decision", "lineage", "note"])


def main():
    print("loading annotations")
    signatures = load_hallmark()
    for k, v in signatures.items():
        print(f"  {k}: {len(v)} genes")
    entrez_sym, ens_map = entrez_to_symbol()
    print(f"  ensembl symbols {len(ens_map)}")
    hugene = hugene_probe_map(entrez_sym)
    print(f"  unambiguous hugene probes {len(hugene)}")
    ilmn = illumina_probe_map()
    print(f"  illumina probes with symbols {len(ilmn)}")

    print("loading matrices")
    mats = {}
    mats["GSE273409"] = load_gse273409()
    mats["GSE242255"] = load_gse242255(ens_map)
    mats["GSE116765"] = load_featurecounts(DATA / "GSE116765_DNAPK_AllTreats_countTable.txt.gz", ens_map)
    mats["GSE116765_22"] = load_featurecounts(DATA / "GSE116765_22Rv1_DNAPKi_CountTable.txt.gz", ens_map)
    mats["GSE167956"] = load_gse167956()
    mats["GSE180581"] = load_gse180581(ens_map)
    mats["GSE268415"] = load_gse268415()
    mats["GSE285698"] = load_gse285698(ens_map)
    mats["GSE129436"] = load_gse129436(ens_map)
    mats["GSE315862"] = load_gse315862()
    mats["GSE319513_24"] = load_gse319513("24")
    mats["GSE319513_48"] = load_gse319513("48")
    mats["GSE93074"] = load_illumina(DATA / "GSE93074_non-normalized-nu7441.txt.gz", ilmn)
    mats["GSE93075"] = load_illumina(DATA / "GSE93075_non-normalized-siPRKDCTC32.txt.gz", ilmn)
    mats["GSE85202"] = load_illumina(DATA / "GSE85202_non-normalized.txt.gz", ilmn, skip_comment=True)
    # series matrix may be absent until downloaded
    sm = DATA / "GSE63480_series_matrix.txt.gz"
    if not sm.exists():
        import urllib.request
        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63480/matrix/GSE63480_series_matrix.txt.gz"
        print("downloading", sm.name)
        urllib.request.urlretrieve(url, sm)
    mats["GSE63480"] = load_gse63480(hugene)

    for name, (log_df, linear, kind) in mats.items():
        pr = "PRKDC" in log_df.index
        print(f"  {name:16} genes={log_df.shape[0]:6} samples={log_df.shape[1]:3} kind={kind:10} PRKDC={pr}")

    # Contrast catalogue. Flags select meta membership before QC.
    # study_group pools cell lines that come from one deposit before the
    # across-study meta, so six SCLC lines are one study.
    C = []

    def add(mid, gse, study, cell, lineage, modality, label, treat, ctrl, **flags):
        C.append({
            "id": mid, "matrix": gse, "gse": gse.split("_")[0] if gse.startswith("GSE319513") or gse.startswith("GSE116765") else gse,
            "study": study, "cell": cell, "lineage": lineage, "modality": modality,
            "label": label, "treat": treat, "ctrl": ctrl, **flags,
        })

    lines = ["DMS114", "H146", "H196", "H446", "H720", "H1341"]
    for line in lines:
        add(
            f"GSE273409_{line}_NU7441", "GSE273409", "GSE273409", line,
            "epithelial_cancer", "selective_inhibitor",
            f"{line} NU7441 vs DMSO",
            [f"{line}_NU_{i}" for i in (1, 2, 3)],
            [f"{line}_control_{i}" for i in (1, 2, 3)],
            epi_combined=True, epi_drug=True, epi_genetic=False, nonepi_combined=False,
        )
    add("GSE242255_KD", "GSE242255", "GSE242255", "CWR22Rv1-AR-EK",
        "epithelial_cancer", "genetic", "siDNA-PKcs vs siScr",
        ["KD_1", "KD_2", "KD_3"], ["Control_1", "Control_2", "Control_3"],
        epi_combined=True, epi_drug=False, epi_genetic=True, nonepi_combined=False)
    add("GSE242255_AZD7648", "GSE242255", "GSE242255", "CWR22Rv1-AR-EK",
        "epithelial_cancer", "selective_inhibitor", "AZD7648 vs DMSO",
        ["AZD7648_1", "AZD7648_2", "AZD7648_3"], ["Control_1", "Control_2", "Control_3"],
        epi_combined=False, epi_drug=True, epi_genetic=False, nonepi_combined=False)
    add("GSE242255_NU7441", "GSE242255", "GSE242255", "CWR22Rv1-AR-EK",
        "epithelial_cancer", "selective_inhibitor", "NU7441 vs DMSO",
        ["NU7441_1", "NU7441_2", "NU7441_3"], ["Control_1", "Control_2", "Control_3"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    add("GSE242255_NU5455", "GSE242255", "GSE242255", "CWR22Rv1-AR-EK",
        "epithelial_cancer", "selective_inhibitor", "NU5455 vs DMSO",
        ["NU5455_1", "NU5455_2", "NU5455_3"], ["Control_1", "Control_2", "Control_3"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    add("GSE116765_C42_NU7441", "GSE116765", "GSE116765", "C4-2",
        "epithelial_cancer", "selective_inhibitor", "NU7441 vs control",
        ["Nu_1", "Nu_2", "Nu_3"], ["Ctrl_1", "Ctrl_2", "Ctrl_3"],
        epi_combined=True, epi_drug=True, epi_genetic=False, nonepi_combined=False)
    add("GSE116765_C42_CC115", "GSE116765", "GSE116765", "C4-2",
        "epithelial_cancer", "dual_inhibitor", "CC-115 vs control",
        ["CC-115_1", "CC-115_2", "CC-115_3"], ["Ctrl_1", "Ctrl_2", "Ctrl_3"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    add("GSE116765_C42_CC223", "GSE116765", "GSE116765", "C4-2",
        "epithelial_cancer", "specificity_control", "CC-223 TORKi vs control",
        ["CC-223_1", "CC-223_2", "CC-223_3"], ["Ctrl_1", "Ctrl_2", "Ctrl_3"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    add("GSE116765_22Rv1_CC115", "GSE116765_22", "GSE116765_22Rv1", "22Rv1",
        "epithelial_cancer", "dual_inhibitor", "CC-115 vs control",
        ["22Rv1_CC-115_1", "22Rv1_CC-115_2"],
        ["22Rv1_Control_1", "22Rv1_Control_2", "22Rv1_Control_3"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    # The count matrix spells the hairpin "shDNS-PKcs" (deposit typo). GEO sample
    # titles use shDNA-PKcs. Same three replicates.
    add("GSE167956_normoxia", "GSE167956", "GSE167956", "SUM159",
        "epithelial_cancer", "genetic", "shDNA-PKcs vs shNTC, 20% O2",
        [f"SUM159-shDNS-PKcs-20pct-{i}" for i in (1, 2, 3)],
        [f"SUM159-shNTC-20pct-{i}" for i in (1, 2, 3)],
        epi_combined=True, epi_drug=False, epi_genetic=True, nonepi_combined=False)
    add("GSE167956_hypoxia", "GSE167956", "GSE167956", "SUM159",
        "epithelial_cancer", "genetic", "shDNA-PKcs vs shNTC, 1% O2",
        [f"SUM159-shDNS-PKcs-1pct-{i}" for i in (1, 2, 3)],
        [f"SUM159-shNTC-1pct-{i}" for i in (1, 2, 3)],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    add("GSE319513_24h", "GSE319513_24", "GSE319513", "8505C",
        "epithelial_cancer", "selective_inhibitor", "KU-57788 vs CON, 24 h",
        [f"FPKM.KU_24H{i}" for i in (1, 2, 3)],
        [f"FPKM.CON_24H{i}" for i in (1, 2, 3)],
        epi_combined=True, epi_drug=True, epi_genetic=False, nonepi_combined=False)
    add("GSE319513_48h", "GSE319513_48", "GSE319513", "8505C",
        "epithelial_cancer", "selective_inhibitor", "KU-57788 vs CON, 48 h",
        [f"FPKM.KU_48H{i}" for i in (1, 2, 3)],
        [f"FPKM.CON_48H{i}" for i in (1, 2, 3)],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    add("GSE93075_HCT116_si", "GSE93075", "GSE93075", "HCT116",
        "epithelial_cancer", "genetic", "siPRKDC vs siNeg",
        [f"HCT116_siPRKDC_r{i}" for i in (1, 2, 3)],
        [f"HCT116_siNeg_r{i}" for i in (1, 2, 3)],
        epi_combined=True, epi_drug=False, epi_genetic=True, nonepi_combined=False)
    add("GSE93074_HCT116_NU", "GSE93074", "GSE93076_drug", "HCT116",
        "epithelial_cancer", "selective_inhibitor", "NU7441 vs vehicle",
        [f"HCT116_T{i}" for i in (1, 2, 3)],
        [f"HCT116_V{i}" for i in (1, 2, 3)],
        epi_combined=False, epi_drug=True, epi_genetic=False, nonepi_combined=False)
    add("GSE63480_KD", "GSE63480", "GSE63480", "C4-2",
        "epithelial_cancer", "genetic", "DNA-PKcs siRNA vs control",
        ["GSM1550591", "GSM1550592"], ["GSM1550583", "GSM1550584"],
        epi_combined=True, epi_drug=False, epi_genetic=True, nonepi_combined=False)
    add("GSE63480_NU", "GSE63480", "GSE63480_drug", "C4-2",
        "epithelial_cancer", "selective_inhibitor", "NU7441 vs control (array)",
        ["GSM1550587", "GSM1550588"], ["GSM1550583", "GSM1550584"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    add("GSE285698_normoxia", "GSE285698", "GSE285698", "HCT116",
        "epithelial_cancer", "genetic", "DNA-PK KO vs WT, normoxia",
        ["A2_S51", "B2_S57", "C2_S63"], ["A1_S50", "B1_S56", "C1_S62"],
        epi_combined=True, epi_drug=False, epi_genetic=True, nonepi_combined=False)
    add("GSE285698_hypoxia", "GSE285698", "GSE285698_hyp", "HCT116",
        "epithelial_cancer", "genetic", "DNA-PK KO vs WT, CoCl2",
        ["A5_S54", "B5_S60", "C5_S66"], ["A4_S53", "B4_S59", "C4_S65"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    add("GSE85202_TC32", "GSE85202", "GSE85202", "TC32",
        "non_epithelial_cancer", "selective_inhibitor", "NU7441 vs vehicle",
        [f"TC32_T{i}" for i in (1, 2, 3)], [f"TC32_V{i}" for i in (1, 2, 3)],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False, nonepi_drug=True)
    add("GSE85202_SKNMC", "GSE85202", "GSE85202", "SK-N-MC",
        "non_epithelial_cancer", "selective_inhibitor", "NU7441 vs vehicle",
        [f"SKNMC_T{i}" for i in (1, 2, 3)], [f"SKNMC_V{i}" for i in (1, 2, 3)],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=True, nonepi_drug=True)
    add("GSE93075_TC32", "GSE93075", "GSE93076_TC32", "TC32",
        "non_epithelial_cancer", "genetic", "siPRKDC vs siNeg",
        [f"TC32_siPRKDC_r{i}" for i in (1, 2, 3)], [f"TC32_siNeg_r{i}" for i in (1, 2, 3)],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=True, nonepi_genetic=True)
    add("GSE268415_si", "GSE268415", "GSE268415", "HOS",
        "non_epithelial_cancer", "genetic", "siPRKDC vs siNC",
        ["hos_si1", "hos_si2", "hos_si3"], ["hos_siNC1", "hos_siNC2", "hos_siNC3"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=True, nonepi_genetic=True)
    add("GSE268415_AZD", "GSE268415", "GSE268415_drug", "HOS",
        "non_epithelial_cancer", "selective_inhibitor", "AZD7648 vs DMSO",
        ["hos_zdk641501", "hos_zdk641502", "hos_zdk641503"],
        ["hos_zdk64NC1", "hos_zdk64NC2", "hos_zdk64NC3"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False, nonepi_drug=True)
    add("GSE315862_LPS853", "GSE315862", "GSE315862", "LPS853",
        "non_epithelial_cancer", "selective_inhibitor", "peposertib vs NTC",
        [f"{i}_LPS8530_5_uM_Peposertib_0nM_Doxorubicin_1d_Rep_{r}" for i, r in ((7, 1), (8, 2), (9, 3))],
        [f"{i}_LPS853NTC_0nM_Doxorubicin_1d_Rep_{r}" for i, r in ((1, 1), (2, 2), (3, 3))],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=True, nonepi_drug=True)
    add("GSE315862_T778", "GSE315862", "GSE315862", "T778",
        "non_epithelial_cancer", "selective_inhibitor", "peposertib vs NTC",
        [f"{i}_T7780_5_uM_Peposertib_0nM_Doxorubicin_1d_Rep_{r}" for i, r in ((25, 1), (26, 2), (27, 3))],
        [f"{i}_T778NTC_0nM_Doxorubicin_1d_Rep_{r}" for i, r in ((19, 1), (20, 2), (21, 3))],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=True, nonepi_drug=True)
    add("GSE129436_U937", "GSE129436", "GSE129436", "U937",
        "non_epithelial_cancer", "selective_inhibitor", "NU7441 vs DMSO, no DNA, 8 h",
        ["WT_lipo_2uM_8A", "WT_lipo_2uM_8B", "WT_lipo_2uM_8C"],
        ["WT_lipo_DMSO_8A", "WT_lipo_DMSO_8B", "WT_lipo_DMSO_8C"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=True, nonepi_drug=True)
    add("GSE129436_U937_DNA", "GSE129436", "GSE129436_DNA", "U937",
        "non_epithelial_cancer", "selective_inhibitor", "NU7441 vs DMSO on transfected DNA, 8 h",
        ["WT_DNA_2uM_8A", "WT_DNA_2uM_8B", "WT_DNA_2uM_8C"],
        ["WT_DNA_DMSO_8A", "WT_DNA_DMSO_8B", "WT_DNA_DMSO_8C"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False)
    add("GSE180581_HEK", "GSE180581", "GSE180581", "HEK293T",
        "epithelial_noncancer", "genetic", "siDNA-PKcs vs siControl in DNA-PKcs+/-",
        ["A19", "A20", "A21"], ["A16", "A17", "A18"],
        epi_combined=False, epi_drug=False, epi_genetic=False, nonepi_combined=False, hek=True)

    sig_names = ["IFN", "IFN_alpha", "IFN_gamma", "STING", "APM", "MHCI"]
    rows = []
    gene_rows = []
    missing_cols = []
    for spec in C:
        log_df, linear, kind = mats[spec["matrix"]]
        # column names may have been cleaned; check
        have = set(log_df.columns)
        if any(c not in have for c in spec["treat"] + spec["ctrl"]):
            missing_cols.append((spec["id"], [c for c in spec["treat"] + spec["ctrl"] if c not in have], list(log_df.columns)[:8]))
            print("MISSING", spec["id"], [c for c in spec["treat"] + spec["ctrl"] if c not in have])
            continue
        prkdc = gene_lfc(log_df, linear, spec["treat"], spec["ctrl"], "PRKDC", kind)
        ok = True
        qc = "inhibitor_mRNA_not_required"
        if spec["modality"] == "genetic":
            if prkdc != prkdc:
                ok = False
                qc = "PRKDC_absent"
            elif prkdc >= -0.3:
                ok = False
                qc = "PRKDC_not_down"
            elif prkdc >= -0.5:
                ok = True
                qc = "PRKDC_weak_kd"
            else:
                qc = "PRKDC_down"
        rec = {
            "id": spec["id"], "gse": spec["gse"], "study": spec["study"], "cell": spec["cell"],
            "lineage": spec["lineage"], "modality": spec["modality"], "label": spec["label"],
            "prkdc_lfc": prkdc, "qc": qc, "perturbation_ok": ok,
        }
        for flag in ("epi_combined", "epi_drug", "epi_genetic", "nonepi_combined", "nonepi_drug", "nonepi_genetic", "hek"):
            rec[flag] = bool(spec.get(flag, False)) and ok
        for sn in sig_names:
            sc = score_contrast(log_df, linear, spec["treat"], spec["ctrl"], signatures[sn], kind)
            if sc is None:
                rec[f"{sn}_lfc"] = float("nan")
                rec[f"{sn}_se"] = float("nan")
                rec[f"{sn}_p"] = float("nan")
                rec[f"{sn}_g"] = float("nan")
                rec[f"{sn}_gvar"] = float("nan")
                rec[f"{sn}_ng"] = 0
                rec[f"{sn}_cov"] = 0
                rec[f"{sn}_frac_up"] = float("nan")
            else:
                rec[f"{sn}_lfc"] = sc["lfc"]
                rec[f"{sn}_se"] = sc["se"]
                rec[f"{sn}_p"] = sc["p"]
                rec[f"{sn}_g"] = sc["g"]
                rec[f"{sn}_gvar"] = sc["g_var"]
                rec[f"{sn}_ng"] = sc["n_genes"]
                rec[f"{sn}_cov"] = sc["coverage"]
                rec[f"{sn}_frac_up"] = sc["frac_genes_up"]
                rec[f"{sn}_lfc_w"] = sc["lfc_w"]
                rec[f"{sn}_lfc_med"] = sc["lfc_med"]
                rec["n_treat"] = sc["n_treat"]
                rec["n_ctrl"] = sc["n_ctrl"]
        # Genome-wide median log2 change, detected genes only, as a toxicity/shift check.
        rec["global_median_lfc"] = _global_median(log_df, linear, spec["treat"], spec["ctrl"], kind)
        for g in KEY_GENES:
            gene_rows.append({
                "id": spec["id"], "gse": spec["gse"], "cell": spec["cell"], "gene": g,
                "lfc": gene_lfc(log_df, linear, spec["treat"], spec["ctrl"], g, kind),
            })
        rows.append(rec)
        print(f"{spec['id']:28} PRKDC {prkdc:+6.2f} {qc:28} IFN {rec.get('IFN_lfc', float('nan')):+6.3f} STING {rec.get('STING_lfc', float('nan')):+6.3f} APM {rec.get('APM_lfc', float('nan')):+6.3f}")

    contr = pd.DataFrame(rows)
    contr.to_csv(TAB / "contrasts.tsv", sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(TAB / "gene_lfc.tsv", sep="\t", index=False)
    screen_table().to_csv(TAB / "screen.tsv", sep="\t", index=False)

    def study_pool(df, sig):
        """Pool cell lines inside a study, then return one row per study."""
        out = []
        for study, sub in df.groupby("study", sort=False):
            eff = sub[f"{sig}_lfc"].to_numpy()
            var = (sub[f"{sig}_se"] ** 2).to_numpy()
            if np.any(~np.isfinite(eff)) or np.any(~np.isfinite(var)) or np.any(var <= 0):
                m = np.isfinite(eff) & np.isfinite(var) & (var > 0)
                eff, var = eff[m], var[m]
                sub = sub.loc[m]
            if len(eff) == 0:
                continue
            pooled = dl_random(eff, var)
            out.append({
                "study": study,
                "cells": ",".join(sub["cell"].astype(str)),
                "n_lines": int(len(sub)),
                "modality": ",".join(sorted(set(sub["modality"]))),
                "lfc": pooled["mu"], "se": pooled["se"], "var": pooled["se"] ** 2,
                "line_lfc": ";".join(f"{c}:{v:.3f}" for c, v in zip(sub["cell"], sub[f"{sig}_lfc"])),
            })
        return pd.DataFrame(out)

    meta_rows = []
    loo_rows = []
    study_frames = {}
    analyses = [
        ("epithelial_combined", "epi_combined", "Epithelial cancer, one contrast per study (genetic preferred when both exist)"),
        ("epithelial_drug", "epi_drug", "Epithelial cancer, selective DNA-PKcs inhibitor"),
        ("epithelial_genetic", "epi_genetic", "Epithelial cancer, PRKDC knockdown or knockout"),
        ("nonepithelial_combined", "nonepi_combined", "Non-epithelial cancer, one contrast per study"),
        ("nonepithelial_drug", "nonepi_drug", "Non-epithelial cancer, selective inhibitor"),
        ("nonepithelial_genetic", "nonepi_genetic", "Non-epithelial cancer, PRKDC knockdown"),
        ("hek293t", "hek", "HEK293T epithelial non-cancer"),
    ]
    primary_sigs = ["IFN", "STING", "APM", "MHCI"]
    for sig in ["IFN", "STING", "APM", "MHCI"]:
        if f"{sig}_lfc_w" in contr.columns:
            contr[f"{sig}_w_lfc"] = contr[f"{sig}_lfc_w"]
            contr[f"{sig}_w_se"] = contr[f"{sig}_se"]
    for analysis, flag, desc in analyses:
        sub = contr[contr[flag]].copy()
        for sig in primary_sigs + ["IFN_alpha", "IFN_gamma", "IFN_w", "STING_w", "APM_w", "MHCI_w"]:
            if sub.empty:
                continue
            pooled = study_pool(sub, sig)
            study_frames[(analysis, sig)] = pooled
            if pooled.empty:
                continue
            res = dl_random(pooled["lfc"], pooled["var"])
            meta_rows.append({
                "analysis": analysis, "description": desc, "signature": sig,
                "level": "study", **{k: res[k] for k in ("k", "mu", "se", "lo", "hi", "p", "I2", "tau2", "Q", "n_pos")},
            })
            # cell-line level, for the SCLC-heavy combined set
            eff = sub[f"{sig}_lfc"].to_numpy()
            var = (sub[f"{sig}_se"] ** 2).to_numpy()
            m = np.isfinite(eff) & np.isfinite(var) & (var > 0)
            res_c = dl_random(eff[m], var[m])
            meta_rows.append({
                "analysis": analysis, "description": desc, "signature": sig,
                "level": "cell_line", **{k: res_c[k] for k in ("k", "mu", "se", "lo", "hi", "p", "I2", "tau2", "Q", "n_pos")},
            })
            if analysis in {"epithelial_combined", "epithelial_drug", "epithelial_genetic"} and sig in primary_sigs and len(pooled) > 2:
                for drop in pooled["study"]:
                    keep = pooled[pooled["study"] != drop]
                    r = dl_random(keep["lfc"], keep["var"])
                    loo_rows.append({
                        "analysis": analysis, "signature": sig, "dropped": drop,
                        "k": r["k"], "mu": r["mu"], "se": r["se"], "p": r["p"], "I2": r["I2"],
                        "lo": r["lo"], "hi": r["hi"], "n_pos": r["n_pos"],
                    })

    meta = pd.DataFrame(meta_rows)
    # BH within the three primary signatures of the study-level epithelial combined meta
    mask = (meta["analysis"] == "epithelial_combined") & (meta["level"] == "study") & (meta["signature"].isin(["IFN", "STING", "APM"]))
    if mask.any():
        meta.loc[mask, "q_bh"] = bh(meta.loc[mask, "p"].to_numpy())
    meta.to_csv(TAB / "meta.tsv", sep="\t", index=False)
    pd.DataFrame(loo_rows).to_csv(TAB / "leave_one_study_out.tsv", sep="\t", index=False)

    # study-level estimates used in the headline forest
    pieces = []
    for sig in primary_sigs:
        fr = study_frames.get(("epithelial_combined", sig))
        if fr is None or fr.empty:
            continue
        t = fr.copy()
        t["signature"] = sig
        pieces.append(t)
    if pieces:
        pd.concat(pieces, ignore_index=True).to_csv(TAB / "epithelial_study_estimates.tsv", sep="\t", index=False)

    plot_forest(contr, meta, study_frames, FIG / "01_forest_epithelial_study.png")
    plot_cell_forest(contr, FIG / "02_forest_epithelial_cell_lines.png")
    plot_heatmap(contr, FIG / "03_signature_heatmap.png")
    plot_prkdc(contr, FIG / "04_prkdc_qc.png")

    summary = {
        "n_contrasts_scored": int(len(contr)),
        "signature_sizes": {k: len(v) for k, v in signatures.items()},
        "meta": meta.to_dict(orient="records"),
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2))
    print("wrote", TAB)
    if missing_cols:
        print("MISSING COLS", len(missing_cols))


def plot_forest(contr, meta, study_frames, path):
    sigs = ["IFN", "STING", "APM"]
    titles = {
        "IFN": "IFN\nHallmark α+γ, APM genes removed",
        "STING": "STING\ntranscriptional readout",
        "APM": "APM\nMHC-I machinery",
    }
    base = study_frames.get(("epithelial_combined", "IFN"))
    if base is None or base.empty:
        return
    studies = list(base["study"])
    nice = {
        "GSE273409": "GSE273409 SCLC, NU7441 (6 lines)",
        "GSE242255": "GSE242255 CWR22Rv1, siPRKDC",
        "GSE116765": "GSE116765 C4-2, NU7441",
        "GSE167956": "GSE167956 SUM159, shPRKDC",
        "GSE319513": "GSE319513 8505C, KU-57788",
        "GSE93075": "GSE93075 HCT116, siPRKDC",
        "GSE63480": "GSE63480 C4-2, siPRKDC",
        "GSE285698": "GSE285698 HCT116, DNA-PK KO",
    }
    labels = [nice.get(s, s) for s in studies] + ["Random-effects mean"]
    y = np.arange(len(studies), 0, -1)
    y_sum = 0
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 5.6), sharey=False)
    for ax, sig in zip(axes, sigs):
        fr = study_frames[("epithelial_combined", sig)].set_index("study").loc[studies]
        ax.axvline(0, color="#888888", lw=0.8)
        ax.errorbar(
            fr["lfc"], y, xerr=1.96 * fr["se"],
            fmt="o", color="#1f4e79", ecolor="#1f4e79", elinewidth=1, capsize=2, ms=5,
        )
        mrow = meta[(meta.analysis == "epithelial_combined") & (meta.level == "study") & (meta.signature == sig)].iloc[0]
        ax.plot([mrow["lo"], mrow["hi"]], [y_sum, y_sum], color="#b85c38", lw=2, zorder=2)
        ax.plot(mrow["mu"], y_sum, marker="D", color="#b85c38", ms=7, zorder=3)
        ax.set_yticks(list(y) + [y_sum])
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylim(-0.7, len(studies) + 0.7)
        ax.set_xlim(-1.8, 1.8)
        ax.set_xlabel("log2 change (95% CI)")
        ax.set_title(titles[sig], fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("Epithelial cancer: one estimate per study", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_cell_forest(contr, path):
    sub = contr[contr["epi_combined"]].copy()
    if sub.empty:
        return
    sigs = ["IFN", "STING", "APM"]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 6.4), sharey=True)
    labels = [f"{r.cell} ({r.gse})" for r in sub.itertuples()]
    y = np.arange(len(sub))[::-1]
    for ax, sig in zip(axes, sigs):
        ax.axvline(0, color="#888888", lw=0.8)
        ax.errorbar(
            sub[f"{sig}_lfc"], y, xerr=1.96 * sub[f"{sig}_se"],
            fmt="o", color="#1f4e79", ecolor="#8aa0b8", elinewidth=1, capsize=2, ms=4,
        )
        ax.set_title(sig, fontsize=10)
        ax.set_xlabel("log2 change")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels, fontsize=7)
    fig.suptitle("Epithelial cancer cell-line contrasts in the combined meta", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_heatmap(contr, path):
    show = contr.copy()
    show["row"] = show["cell"] + "  " + show["label"]
    mat = show.set_index("row")[["IFN_lfc", "STING_lfc", "APM_lfc"]]
    fig_h = max(6, 0.28 * len(mat))
    fig, ax = plt.subplots(figsize=(6.4, fig_h))
    vals = mat.to_numpy(dtype=float)
    lim = np.nanmax(np.abs(vals))
    lim = max(0.4, min(1.5, lim))
    im = ax.imshow(vals, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["IFN", "STING", "APM"])
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=7)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6, color="black")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="log2 change")
    ax.set_title("All scored DNA-PKcs inhibitor / PRKDC KD contrasts")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_prkdc(contr, path):
    show = contr.copy()
    show["row"] = show["cell"] + " | " + show["modality"]
    fig, ax = plt.subplots(figsize=(8, max(4, 0.28 * len(show))))
    y = np.arange(len(show))[::-1]
    colors = ["#1f4e79" if m == "genetic" else "#b85c38" if "inhibitor" in m else "#777777" for m in show["modality"]]
    ax.barh(y, show["prkdc_lfc"], color=colors)
    ax.axvline(0, color="black", lw=0.6)
    ax.axvline(-0.5, color="#888888", ls="--", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(show["row"] + " (" + show["gse"] + ")", fontsize=7)
    ax.set_xlabel("PRKDC log2 change (KD/KO should be negative; inhibitors need not be)")
    ax.set_title("Perturbation check")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
