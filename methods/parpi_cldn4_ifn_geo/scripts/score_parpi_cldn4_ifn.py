#!/usr/bin/env python3
"""Score public PARP-inhibitor RNA-seq for CLDN4, tight junctions, and IFN.

Pre-specified genes only. No genome-wide discovery. Counts are library-size
normalized to CPM before log2. Author-normalized matrices are not re-CPM'd.
Mann-Whitney exact p is reported only when both groups have n >= 3.
n < 3 is a descriptive log2FC with p left blank.
"""

from __future__ import annotations

import gzip
import io
import json
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = ROOT / "cache"
OUT = ROOT / "results" / "tables"
CACHE.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

CORE_IFN = [
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
    "OASL", "IFI44", "IFI44L", "RSAD2", "STAT1", "STAT2", "IRF7", "CXCL9",
    "CXCL10", "CXCL11", "CCL5", "CD274", "DDX58", "IFIH1", "ISG20", "BST2",
    "USP18", "IFNB1", "IFNG",
]
CORE_TJ = [
    "CLDN1", "CLDN3", "CLDN7", "CLDN18", "OCLN", "TJP1", "TJP2", "TJP3",
    "F11R", "CGN", "CGNL1", "MARVELD2", "JAM2", "JAM3",
]
REPAIR = [
    "TP53BP1", "XRCC1", "PRKDC", "LIG4", "XRCC4", "NHEJ1", "PARP1", "BRCA1",
    "BRCA2", "RAD51", "RAD51C", "BRIP1",
]
STING = ["CGAS", "STING1", "TMEM173", "TBK1", "IRF3", "MAVS"]
CALLOUT = ["CLDN4", "TACSTD2", "CXCL10", "CCL5", "CD274", "TP53BP1", "XRCC1", "STAT1", "ISG15", "IFNB1"]

MOUSE = {
    "ISG15": "Isg15", "IFIT1": "Ifit1", "IFIT2": "Ifit2", "IFIT3": "Ifit3",
    "MX1": "Mx1", "MX2": "Mx2", "OAS1": "Oas1a", "OAS2": "Oas2", "OAS3": "Oas3",
    "OASL": "Oasl1", "IFI44": "Ifi44", "IFI44L": "Ifi44l", "RSAD2": "Rsad2",
    "STAT1": "Stat1", "STAT2": "Stat2", "IRF7": "Irf7", "CXCL9": "Cxcl9",
    "CXCL10": "Cxcl10", "CXCL11": "Cxcl11", "CCL5": "Ccl5", "CD274": "Cd274",
    "DDX58": "Ddx58", "IFIH1": "Ifih1", "ISG20": "Isg20", "BST2": "Bst2",
    "USP18": "Usp18", "IFNB1": "Ifnb1", "IFNG": "Ifng",
    "CLDN1": "Cldn1", "CLDN3": "Cldn3", "CLDN4": "Cldn4", "CLDN7": "Cldn7",
    "CLDN18": "Cldn18", "OCLN": "Ocln", "TJP1": "Tjp1", "TJP2": "Tjp2",
    "TJP3": "Tjp3", "F11R": "F11r", "CGN": "Cgn", "CGNL1": "Cgnl1",
    "MARVELD2": "Marveld2", "JAM2": "Jam2", "JAM3": "Jam3",
    "TP53BP1": "Trp53bp1", "XRCC1": "Xrcc1", "PRKDC": "Prkdc", "LIG4": "Lig4",
    "XRCC4": "Xrcc4", "NHEJ1": "Nhej1", "PARP1": "Parp1", "BRCA1": "Brca1",
    "BRCA2": "Brca2", "RAD51": "Rad51", "RAD51C": "Rad51c", "BRIP1": "Brip1",
    "CGAS": "Cgas", "STING1": "Sting1", "TMEM173": "Tmem173", "TBK1": "Tbk1",
    "IRF3": "Irf3", "MAVS": "Mavs", "TACSTD2": "Tacstd2",
}


def download(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print("GET", url)
    urllib.request.urlretrieve(url, dest)
    return dest


def read_lines(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def ensembl_map(symbols: list[str]) -> dict[str, str]:
    cache = CACHE / "ensembl_human_symbols.json"
    if cache.exists():
        raw = json.loads(cache.read_text())
    else:
        raw = {}
    missing = [s for s in symbols if s not in raw]
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    for i in range(0, len(missing), 80):
        chunk = missing[i : i + 80]
        req = urllib.request.Request(
            url,
            data=json.dumps({"symbols": chunk}).encode(),
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=90) as resp:
            got = json.loads(resp.read().decode())
        for sym in chunk:
            rec = got.get(sym)
            raw[sym] = None if not rec else rec.get("id")
    cache.write_text(json.dumps(raw))
    return {s: raw[s] for s in symbols if raw.get(s)}


def gzip_text(path: Path):
    return gzip.open(path, "rt")


def load_symbol_matrix(path: Path) -> pd.DataFrame:
    """Read a gene-by-sample table. Some GEO files omit the gene-column header."""
    compression = "gzip" if str(path).endswith(".gz") else None
    peek = pd.read_csv(path, sep="\t", compression=compression, nrows=5)
    first = str(peek.columns[0])
    # A real gene-column header is Gene, Geneid, ID, or similar — not a sample name.
    named_gene = first.lower() in {"gene", "geneid", "gene_id", "id", "symbol", "gene_name"} or first.startswith("Gene")
    if named_gene:
        df = pd.read_csv(path, sep="\t", compression=compression)
        gene = df.columns[0]
        df = df.rename(columns={gene: "gene"})
    else:
        df = pd.read_csv(path, sep="\t", compression=compression, header=None, skiprows=1)
        header = pd.read_csv(path, sep="\t", compression=compression, nrows=0)
        df.columns = ["gene"] + list(header.columns)
    df["gene"] = df["gene"].astype(str).str.split(".").str[0]
    value_cols = [c for c in df.columns if c != "gene"]
    num = df[value_cols].apply(pd.to_numeric, errors="coerce")
    num.insert(0, "gene", df["gene"].values)
    return num.groupby("gene", as_index=True).sum(numeric_only=True)


def load_gse233820(path: Path) -> pd.DataFrame:
    """Header is sample names only; the gene id is an extra leading field."""
    with gzip_text(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = []
        genes = []
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[0].split(".")[0])
            rows.append([float(x) if x not in ("", "NA") else np.nan for x in parts[1:]])
    if len(rows[0]) != len(header):
        raise SystemExit(f"GSE233820 shape mismatch {len(rows[0])} vs {len(header)}")
    df = pd.DataFrame(rows, index=genes, columns=header)
    return df.groupby(level=0).sum()


def cpm(df: pd.DataFrame) -> pd.DataFrame:
    lib = df.sum(axis=0).replace(0, np.nan)
    return df.div(lib, axis=1) * 1e6


def log2p(df: pd.DataFrame, pseudo: float = 1.0) -> pd.DataFrame:
    return np.log2(df.astype(float) + pseudo)


def welch(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    if np.nanstd(a) == 0 and np.nanstd(b) == 0:
        return 1.0 if np.allclose(a.mean(), b.mean()) else np.nan
    res = stats.ttest_ind(a, b, equal_var=False, alternative="two-sided")
    return float(res.pvalue)


def mwu(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan
    if np.all(a == a[0]) and np.all(b == b[0]) and a[0] == b[0]:
        return 1.0
    res = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
    return float(res.pvalue)


def score_contrast(
    log_df: pd.DataFrame,
    treat: list[str],
    ctrl: list[str],
    families: dict[str, list[str]],
    contrast_id: str,
    meta: dict,
) -> tuple[list[dict], list[dict]]:
    missing = [c for c in treat + ctrl if c not in log_df.columns]
    if missing:
        raise SystemExit(f"{contrast_id} missing columns {missing}")
    sub = log_df[treat + ctrl]
    gene_rows = []
    # callout + every family member that is present
    wanted = []
    for genes in families.values():
        wanted.extend(genes)
    wanted = list(dict.fromkeys(wanted))
    present = [g for g in wanted if g in sub.index]
    for gene in present:
        tv = sub.loc[gene, treat].to_numpy(dtype=float)
        cv = sub.loc[gene, ctrl].to_numpy(dtype=float)
        gene_rows.append(
            {
                "contrast_id": contrast_id,
                "gene": gene,
                "n_treat": int(np.isfinite(tv).sum()),
                "n_ctrl": int(np.isfinite(cv).sum()),
                "mean_treat": float(np.nanmean(tv)),
                "mean_ctrl": float(np.nanmean(cv)),
                "log2fc": float(np.nanmean(tv) - np.nanmean(cv)),
                "mw_p": mwu(tv, cv),
                "welch_p": welch(tv, cv),
                **meta,
            }
        )
    fam_rows = []
    for fam, genes in families.items():
        use = []
        for g in genes:
            if g not in sub.index:
                continue
            vals = sub.loc[g, treat + ctrl].to_numpy(dtype=float)
            if np.nanmax(vals) > 0 or np.nanmin(vals) < 0:
                # prelogged matrices can be entirely positive; keep expressed genes
                if np.nanmax(np.abs(vals)) > 0:
                    use.append(g)
        if not use:
            fam_rows.append(
                {
                    "contrast_id": contrast_id,
                    "family": fam,
                    "n_genes": 0,
                    "n_genes_up": 0,
                    "median_gene_log2fc": np.nan,
                    "mean_gene_log2fc": np.nan,
                    "mw_p": np.nan,
                    "welch_p": np.nan,
                    **meta,
                }
            )
            continue
        block = sub.loc[use]
        gene_lfc = block[treat].mean(axis=1) - block[ctrl].mean(axis=1)
        treat_score = block[treat].mean(axis=0).to_numpy(dtype=float)
        ctrl_score = block[ctrl].mean(axis=0).to_numpy(dtype=float)
        fam_rows.append(
            {
                "contrast_id": contrast_id,
                "family": fam,
                "n_genes": len(use),
                "n_genes_up": int((gene_lfc > 0).sum()),
                "median_gene_log2fc": float(np.nanmedian(gene_lfc)),
                "mean_gene_log2fc": float(np.nanmean(gene_lfc)),
                "mw_p": mwu(treat_score, ctrl_score),
                "welch_p": welch(treat_score, ctrl_score),
                **meta,
            }
        )
    return gene_rows, fam_rows


def families_for(index: pd.Index, species: str) -> dict[str, list[str]]:
    ifn = read_lines(DATA / "family_IFN.txt")
    tj = read_lines(DATA / "family_TJ_no_CLDN4.txt")
    if species == "mouse":
        def conv(genes: list[str]) -> list[str]:
            out = []
            for g in genes:
                if g in MOUSE:
                    out.append(MOUSE[g])
            # also keep symbols already mouse-style if present
            return out

        mapped = {
            "CLDN4": conv(["CLDN4"]),
            "TACSTD2": conv(["TACSTD2"]),
            "core_IFN": conv(CORE_IFN),
            "core_TJ": conv(CORE_TJ),
            "repair_NHEJ_HR": conv(REPAIR),
            "STING_axis": conv(STING),
            "recruit": conv(["CXCL9", "CXCL10", "CXCL11", "CCL5", "CD274"]),
        }
    else:
        mapped = {
            "CLDN4": ["CLDN4"],
            "TACSTD2": ["TACSTD2"],
            "core_IFN": CORE_IFN,
            "core_TJ": CORE_TJ,
            "repair_NHEJ_HR": REPAIR,
            "STING_axis": STING,
            "recruit": ["CXCL9", "CXCL10", "CXCL11", "CCL5", "CD274"],
            "thesis_IFN": ifn,
            "thesis_TJ_no_CLDN4": tj,
        }
    # drop genes absent, but keep the name list; scorer filters
    upper = {str(i).upper(): i for i in index}
    resolved = {}
    for fam, genes in mapped.items():
        hit = []
        for g in genes:
            if g in index:
                hit.append(g)
            elif g.upper() in upper:
                hit.append(upper[g.upper()])
        resolved[fam] = hit
    return resolved


def translate_ensembl(df: pd.DataFrame, symbols_needed: list[str]) -> pd.DataFrame:
    mapping = ensembl_map(symbols_needed)
    id_to_sym = {}
    for sym, ensg in mapping.items():
        id_to_sym.setdefault(ensg, sym)
    keep = df.loc[df.index.isin(id_to_sym)].copy()
    keep.index = [id_to_sym[i] for i in keep.index]
    keep = keep.groupby(level=0).sum()
    return keep


def human_symbols_needed() -> list[str]:
    ifn = read_lines(DATA / "family_IFN.txt")
    tj = read_lines(DATA / "family_TJ_no_CLDN4.txt")
    return sorted(set(CORE_IFN + CORE_TJ + REPAIR + STING + CALLOUT + ifn + tj + ["CXCL9", "CXCL11"]))


def prepare_counts_log(df_counts: pd.DataFrame) -> pd.DataFrame:
    return log2p(cpm(df_counts))


def main() -> None:
    gene_rows: list[dict] = []
    fam_rows: list[dict] = []
    notes = []

    # ----- GSE237361 featureCounts, symbols, raw counts -----
    p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237361/suppl/GSE237361_featureCounts_for_DESeq2_annotated.txt.gz",
        CACHE / "GSE237361.txt.gz",
    )
    m237 = load_symbol_matrix(p)
    log237 = prepare_counts_log(m237)
    fam237 = families_for(log237.index, "human")
    specs = [
        (
            "GSE237361_UWB_olaparib_96h",
            [c for c in log237.columns if c.startswith("E") and "UWB_olaparib" in c],
            [c for c in log237.columns if c.startswith("E") and "UWB_DMSO" in c],
            "acute PARPi",
            "UWB1.289 BRCA1-mutant, olaparib 3.5 uM vs DMSO, 96h, n=4. Counts to CPM.",
        ),
        (
            "GSE237361_OVCAR3_olaparib_96h",
            [c for c in log237.columns if "OVC_olaparib" in c],
            [c for c in log237.columns if "OVC_DMSO" in c],
            "acute PARPi",
            "OVCAR3 BRCA1-WT, olaparib 7.5 uM vs DMSO, 96h, n=4. Same line as Yamamoto CLDN4-high model. Counts to CPM.",
        ),
    ]
    for cid, treat, ctrl, kind, note in specs:
        g, f = score_contrast(log237, treat, ctrl, fam237, cid, {"series": "GSE237361", "kind": kind, "species": "human", "note": note})
        gene_rows.extend(g)
        fam_rows.extend(f)
        notes.append({"contrast_id": cid, "kind": kind, "note": note, "n_treat": len(treat), "n_ctrl": len(ctrl)})

    # baseline CLDN4
    base_rows = []
    for label, cols in {
        "UWB_DMSO": [c for c in log237.columns if "UWB_DMSO" in c],
        "OVCAR3_DMSO": [c for c in log237.columns if "OVC_DMSO" in c],
        "UWB_olaparib": [c for c in log237.columns if "UWB_olaparib" in c],
        "OVCAR3_olaparib": [c for c in log237.columns if "OVC_olaparib" in c],
    }.items():
        base_rows.append(
            {
                "series": "GSE237361",
                "group": label,
                "gene": "CLDN4",
                "mean_log2cpm1": float(log237.loc["CLDN4", cols].mean()) if "CLDN4" in log237.index else np.nan,
                "n": len(cols),
            }
        )

    # ----- GSE243208 long counts, ENSG -----
    p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243208/suppl/GSE243208_Olaparib_UWB1289_gene_counts.tsv.gz",
        CACHE / "GSE243208.tsv.gz",
    )
    long = pd.read_csv(p, sep="\t", compression="gzip")
    long["Geneid"] = long["Geneid"].astype(str).str.split(".").str[0]
    wide = long.pivot_table(index="Geneid", columns="sample", values="count", aggfunc="sum")
    sym = translate_ensembl(wide, human_symbols_needed())
    # full library size must use all genes, then subset symbols
    lib = wide.sum(axis=0)
    sym_cpm = sym.div(lib, axis=1) * 1e6
    log243 = log2p(sym_cpm)
    colmap = {
        "MP1237AA4P4": "ola_rep1",
        "MP1237AA5P5": "ola_rep2",
        "MP1237AA6P6": "ola_rep3",
        "MP1237AA28P28": "dmso_rep1",
        "MP1237AA29P29": "dmso_rep2",
        "MP1237AA30P30": "dmso_rep3",
    }
    log243 = log243.rename(columns=colmap)
    fam243 = families_for(log243.index, "human")
    note = (
        "UWB1.289, 24h, n=3. GSM characteristics say Olaparib 4 uM; series overall design says 10 uM. "
        "GSM also says genotype WT; UWB1.289 is the BRCA1-null line, so that field is not used as BRCA1 status. "
        "CPM uses full-gene library size; only pre-specified symbols were ID-mapped."
    )
    g, f = score_contrast(
        log243,
        ["ola_rep1", "ola_rep2", "ola_rep3"],
        ["dmso_rep1", "dmso_rep2", "dmso_rep3"],
        fam243,
        "GSE243208_UWB_olaparib_24h",
        {"series": "GSE243208", "kind": "acute PARPi", "species": "human", "note": note},
    )
    gene_rows.extend(g)
    fam_rows.extend(f)
    notes.append({"contrast_id": "GSE243208_UWB_olaparib_24h", "kind": "acute PARPi", "note": note, "n_treat": 3, "n_ctrl": 3})

    # ----- GSE233820 abundances -----
    p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233820/suppl/GSE233820_Raw_Abundances.txt.gz",
        CACHE / "GSE233820.txt.gz",
    )
    raw233 = load_gse233820(p)
    integer_frac = float(np.mean(np.isclose(raw233.to_numpy(), np.round(raw233.to_numpy()))))
    col_sums = raw233.sum(axis=0)
    # Raw abundances with near-integer values are treated as counts.
    if integer_frac > 0.9:
        log233_all = prepare_counts_log(raw233)
        scale_note = f"integer fraction {integer_frac:.3f}; treated as counts and converted to CPM. Column sums {int(col_sums.min())}-{int(col_sums.max())}."
    else:
        log233_all = log2p(raw233)
        scale_note = f"integer fraction {integer_frac:.3f}; treated as already-scaled abundances, log2(x+1), no CPM."
    sym233 = translate_ensembl(raw233, human_symbols_needed())
    if integer_frac > 0.9:
        lib = raw233.sum(axis=0)
        log233 = log2p(sym233.div(lib, axis=1) * 1e6)
    else:
        log233 = log2p(sym233)
    fam233 = families_for(log233.index, "human")
    for cid, treat_key, ctrl_key, kind in [
        ("GSE233820_SBC5_PARPi_0Gy", "PARPi_0Gy", "DMSO_0Gy", "acute PARPi"),
        ("GSE233820_SBC5_PARPi_on_4Gy", "PARPi_4Gy", "DMSO_4Gy", "PARPi plus radiation"),
    ]:
        treat = [c for c in log233.columns if treat_key in c]
        ctrl = [c for c in log233.columns if ctrl_key in c]
        note = (
            f"SBC5 small-cell lung cancer, olaparib 1 uM, 3 days, n=3. {scale_note} "
            "Not LUAD and not KL. GEO series title says chemokine translation is depressed; "
            "PubMed 40038278 title says chemokine mRNA is stabilized."
        )
        g, f = score_contrast(log233, treat, ctrl, fam233, cid, {"series": "GSE233820", "kind": kind, "species": "human", "note": note})
        gene_rows.extend(g)
        fam_rows.extend(f)
        notes.append({"contrast_id": cid, "kind": kind, "note": note, "n_treat": len(treat), "n_ctrl": len(ctrl)})

    # ----- GSE285827 per-sample counts -----
    tar_path = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285827/suppl/GSE285827_RAW.tar",
        CACHE / "GSE285827_RAW.tar",
    )
    frames = []
    with tarfile.open(tar_path) as tar:
        for mem in tar.getmembers():
            if not mem.name.endswith(".gz"):
                continue
            raw = tar.extractfile(mem).read()
            one = pd.read_csv(io.BytesIO(raw), sep="\t", compression="gzip")
            gene_col = one.columns[0]
            sample_col = one.columns[1]
            one = one.rename(columns={gene_col: "gene", sample_col: sample_col})
            one["gene"] = one["gene"].astype(str)
            one = one.groupby("gene", as_index=True)[sample_col].sum()
            frames.append(one)
    m285 = pd.concat(frames, axis=1).fillna(0.0)
    log285 = prepare_counts_log(m285)
    fam285 = families_for(log285.index, "human")
    for line in ("ovcar3", "caov3"):
        for drug, drug_name in (("talap", "talazoparib"), ("velip", "veliparib")):
            treat = [c for c in log285.columns if c.startswith(f"{line}_{drug}_")]
            ctrl = [c for c in log285.columns if c.startswith(f"{line}_dmso_")]
            cid = f"GSE285827_{line}_{drug_name}"
            note = (
                f"{line.upper()} {drug_name} vs DMSO, three experiment dates as n=3. Counts to CPM. "
                "OVCAR3 is the Yamamoto CLDN4-high BRCA-WT line; CAOV3 is a second HGSOC line. "
                "Paper: PARP inhibitors differentially regulate immune responses by genetic background."
            )
            g, f = score_contrast(log285, treat, ctrl, fam285, cid, {"series": "GSE285827", "kind": "acute PARPi", "species": "human", "note": note})
            gene_rows.extend(g)
            fam_rows.extend(f)
            notes.append({"contrast_id": cid, "kind": "acute PARPi", "note": note, "n_treat": len(treat), "n_ctrl": len(ctrl)})

    # ----- GSE120500 mouse AmpliSeq, columns keyed by SRA spot counts -----
    xls_gz = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120500/suppl/GSE120500_Mouse_Transcriptome_control_and_parpi.xls.gz",
        CACHE / "GSE120500.xls.gz",
    )
    xls_path = CACHE / "GSE120500.xls"
    if not xls_path.exists():
        with gzip.open(xls_gz, "rb") as src, open(xls_path, "wb") as dst:
            dst.write(src.read())
    amp = pd.read_excel(xls_path, engine="xlrd")
    amp = amp.rename(columns={amp.columns[0]: "gene"}).groupby("gene", as_index=True).sum(numeric_only=True)
    spots = pd.read_csv(DATA / "gse120500_sra_spots.tsv", sep="\t")
    col_sums = amp.sum(axis=0)
    assign = {}
    # Within Control-3/4/5 the on-panel sums differ by ~0.2%, so replicate
    # identity is not unique. Group identity is unique: L01-L06 match the six
    # vehicle libraries and L07-L12 match the six olaparib libraries.
    for col, total in col_sums.items():
        rel = (spots["spots"] - float(total)).abs() / spots["spots"]
        j = int(rel.to_numpy().argmin())
        if float(rel.iloc[j]) > 0.02:
            raise SystemExit(f"No SRA spot match for {col}")
        group = "Control" if spots.iloc[j]["title"].startswith("Control") else "PARPi"
        other = spots.loc[~spots["title"].str.startswith(group), "spots"]
        other_rel = float(((other - float(total)).abs() / other).min())
        if other_rel < 0.03:
            raise SystemExit(f"Group call ambiguous for {col}: other_rel={other_rel}")
        assign[col] = f"{group}_{col}"
    amp = amp.rename(columns=assign)
    log120 = prepare_counts_log(amp)
    fam120 = families_for(log120.index, "mouse")
    note = (
        "Ding 2018 Cell Reports. Bulk Brca1-deficient ovarian tumors, olaparib vs vehicle, n=6. "
        "Ion AmpliSeq panel of 4604 genes, not whole transcriptome. "
        "Each count column was assigned to vehicle or olaparib by nearest SRA spot count (relative error < 2% to that group, > 3% to the other group). "
        "Control-3/4/5 are too close to name individual replicates; the six-vs-six group split is the contrast used. "
        "IFN signal can come from infiltrating immune cells. Mouse symbols."
    )
    g, f = score_contrast(
        log120,
        [c for c in log120.columns if c.startswith("PARPi")],
        [c for c in log120.columns if c.startswith("Control")],
        fam120,
        "GSE120500_Brca1def_tumor_olaparib",
        {"series": "GSE120500", "kind": "acute PARPi, bulk tumor", "species": "mouse", "note": note},
    )
    gene_rows.extend(g)
    fam_rows.extend(f)
    notes.append({"contrast_id": "GSE120500_Brca1def_tumor_olaparib", "kind": "acute PARPi, bulk tumor", "note": note, "n_treat": 6, "n_ctrl": 6})

    # ----- GSE120792 BRCA1 MUT vs WT, prelogged -----
    p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120792/suppl/GSE120792_BRCA1_all.txt.gz",
        CACHE / "GSE120792.txt.gz",
    )
    m792 = load_symbol_matrix(p)
    if "ACTB" not in m792.index:
        raise SystemExit("GSE120792 missing ACTB")
    actb = float(m792.loc["ACTB"].mean())
    if not (8.0 <= actb <= 12.0):
        raise SystemExit(f"GSE120792 ACTB mean {actb} is not in the expected log2 band")
    log792 = m792  # already log-scaled; do not log again
    fam792 = families_for(log792.index, "human")
    note = (
        "BRCA1-mutant vs BRCA1-WT ovarian cancer cells, n=3, no olaparib column in this matrix. "
        f"Values left as deposited (ACTB mean {actb:.2f}, consistent with log2 expression). "
        "Difference of means is the contrast. Genotype arm of the bridge, not a drug arm. PMID 34289354."
    )
    g, f = score_contrast(
        log792,
        ["BRCA1_MUT1", "BRCA1_MUT2", "BRCA1_MUT3"],
        ["BRCA1_WT1", "BRCA1_WT2", "BRCA1_WT3"],
        fam792,
        "GSE120792_BRCA1mut_vs_WT",
        {"series": "GSE120792", "kind": "DNA-repair genotype, not drug", "species": "human", "note": note},
    )
    gene_rows.extend(g)
    fam_rows.extend(f)
    notes.append({"contrast_id": "GSE120792_BRCA1mut_vs_WT", "kind": "DNA-repair genotype, not drug", "note": note, "n_treat": 3, "n_ctrl": 3})

    # ----- GSE153867 FPKM, resistant vs parent -----
    p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153867/suppl/GSE153867_fpkm.txt.gz",
        CACHE / "GSE153867.txt.gz",
    )
    m867 = load_symbol_matrix(p)
    log867 = log2p(m867)
    fam867 = families_for(log867.index, "human")
    resist = [f"C-{i}" for i in range(1, 9)]
    parent = [f"O-{i}" for i in range(1, 9)]
    note = (
        "A2780 olaparib-resistant (columns C-1..C-8, GEO 'treated with Olaparib') vs parental untreated "
        "(columns O-1..O-8), n=8. FPKM to log2(FPKM+1). Acquired resistance, not an acute isogenic dose. "
        "C13* shCEBPB columns CS/R are a different experiment and are not in this contrast."
    )
    g, f = score_contrast(
        log867, resist, parent, fam867, "GSE153867_A2780_resistant_vs_parent",
        {"series": "GSE153867", "kind": "acquired PARPi resistance", "species": "human", "note": note},
    )
    gene_rows.extend(g)
    fam_rows.extend(f)
    notes.append({"contrast_id": "GSE153867_A2780_resistant_vs_parent", "kind": "acquired PARPi resistance", "note": note, "n_treat": 8, "n_ctrl": 8})

    # ----- GSE235980 counts, n=2 resistance -----
    p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235980/suppl/GSE235980_CountReads.txt.gz",
        CACHE / "GSE235980_counts.txt.gz",
    )
    m980 = load_symbol_matrix(p)
    log980 = prepare_counts_log(m980)
    fam980 = families_for(log980.index, "human")
    for cid, treat, ctrl, note in [
        (
            "GSE235980_Olres_vs_parent_BRCA1def",
            ["AN14028316", "AN14028318"],
            ["AN14028320", "AN14028322"],
            "UWB1.289 Olres vs parental, drug-free culture, n=2. Counts to CPM. Descriptive only.",
        ),
        (
            "GSE235980_Olres_vs_parent_BRCA1prof",
            ["AN14028324", "AN14028326"],
            ["AN14028328", "AN14028330"],
            "UWB1.289+BRCA1 Olres vs parental, drug-free culture, n=2. Counts to CPM. Descriptive only.",
        ),
    ]:
        g, f = score_contrast(log980, treat, ctrl, fam980, cid, {"series": "GSE235980", "kind": "acquired PARPi resistance", "species": "human", "note": note})
        gene_rows.extend(g)
        fam_rows.extend(f)
        notes.append({"contrast_id": cid, "kind": "acquired PARPi resistance", "note": note, "n_treat": 2, "n_ctrl": 2})

    # author DESeq2 sign check
    author_rows = []
    for label, url, name in [
        ("BRCA1-deficient Olres vs parent (author DESeq2)", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235980/suppl/GSE235980_DeSeq2_Deficient.xlsx", "def.xlsx"),
        ("BRCA1-proficient Olres vs parent (author DESeq2)", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235980/suppl/GSE235980_DeSeq2_Proficient.xlsx", "prof.xlsx"),
    ]:
        path = download(url, CACHE / f"GSE235980_{name}")
        de = pd.read_excel(path)
        de["gene"] = de["gene"].astype(str)
        watch = set(CALLOUT + CORE_IFN + CORE_TJ + REPAIR + STING)
        hit = de[de["gene"].isin(watch)][["gene", "log2FoldChange", "pvalue", "padj", "baseMean"]]
        for rec in hit.to_dict(orient="records"):
            rec["author_contrast"] = label
            author_rows.append(rec)

    # ----- GSE239639 TPM, olaparib vs untreated, n=2 per line -----
    p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239639/suppl/GSE239639_Table_Expression_Level_TPM.txt.gz",
        CACHE / "GSE239639.txt.gz",
    )
    m639 = pd.read_csv(p, sep="\t", compression="gzip")
    m639["ID"] = m639["ID"].astype(str).str.split(".").str[0]
    m639 = m639.groupby("ID").sum(numeric_only=True)
    sym639 = translate_ensembl(m639, human_symbols_needed())
    log639 = log2p(sym639)  # TPM, not counts
    fam639 = families_for(log639.index, "human")
    for line in ("DMR", "HT-1080", "S018"):
        treat = [c for c in log639.columns if c.startswith(f"{line}_OL_")]
        ctrl = [c for c in log639.columns if c.startswith(f"{line}_NT_")]
        cid = f"GSE239639_{line}_olaparib_vs_NT"
        note = (
            f"{line} sarcoma, olaparib (OL) vs not treated (NT), n=2. TPM to log2(TPM+1). "
            "Combo and trabectedin columns are not in this contrast. Descriptive only. Not epithelial lung."
        )
        g, f = score_contrast(log639, treat, ctrl, fam639, cid, {"series": "GSE239639", "kind": "acute PARPi, n=2", "species": "human", "note": note})
        gene_rows.extend(g)
        fam_rows.extend(f)
        notes.append({"contrast_id": cid, "kind": "acute PARPi, n=2", "note": note, "n_treat": 2, "n_ctrl": 2})

    # ----- GSE298546 HeLa n=1 -----
    p = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298546/suppl/GSE298546_count.txt.gz",
        CACHE / "GSE298546.txt.gz",
    )
    m546 = load_symbol_matrix(p)
    log546 = prepare_counts_log(m546)
    fam546 = families_for(log546.index, "human")
    cols = list(log546.columns)
    ola = [c for c in cols if "Olaparib" in c]
    dmso = [c for c in cols if "DMSO" in c]
    note = "HeLa, olaparib vs DMSO, 72h, n=1. Counts to CPM. Descriptive only. PMID 40316566 discusses senescence-associated inflammation in HR-proficient cells."
    g, f = score_contrast(log546, ola, dmso, fam546, "GSE298546_HeLa_olaparib_72h", {"series": "GSE298546", "kind": "acute PARPi, n=1", "species": "human", "note": note})
    gene_rows.extend(g)
    fam_rows.extend(f)
    notes.append({"contrast_id": "GSE298546_HeLa_olaparib_72h", "kind": "acute PARPi, n=1", "note": note, "n_treat": 1, "n_ctrl": 1})

    genes = pd.DataFrame(gene_rows)
    fams = pd.DataFrame(fam_rows)
    genes.to_csv(OUT / "gene_contrasts.tsv", sep="\t", index=False, float_format="%.6g")
    fams.to_csv(OUT / "family_contrasts.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(notes).to_csv(OUT / "contrast_notes.tsv", sep="\t", index=False)
    pd.DataFrame(base_rows).to_csv(OUT / "GSE237361_CLDN4_baseline.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(author_rows).to_csv(OUT / "GSE235980_author_DESeq2_keygenes.tsv", sep="\t", index=False, float_format="%.6g")

    # wide key table
    key_fams = ["CLDN4", "TACSTD2", "core_IFN", "core_TJ", "thesis_IFN", "thesis_TJ_no_CLDN4", "repair_NHEJ_HR", "STING_axis", "recruit"]
    key = fams[fams["family"].isin(key_fams)].pivot(index="contrast_id", columns="family", values="mean_gene_log2fc")
    key_p = fams[fams["family"].isin(["CLDN4", "core_IFN", "core_TJ", "recruit"])].pivot(index="contrast_id", columns="family", values="mw_p")
    key_p = key_p.add_suffix("_mw_p")
    wide = key.join(key_p)
    wide.to_csv(OUT / "family_log2fc_wide.tsv", sep="\t", float_format="%.6g")
    print(wide.to_string())
    print("wrote", OUT)


if __name__ == "__main__":
    main()
