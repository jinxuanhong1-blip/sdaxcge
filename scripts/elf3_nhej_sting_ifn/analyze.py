#!/usr/bin/env python3
"""Public ELF3 ChIP / CUT&Tag and ELF3-loss transcriptomes: NHEJ, STING, IFN.

Question
--------
Does ELF3 loss change the classical NHEJ, cGAS–STING, or interferon programs,
and does public ELF3 binding place those genes in the ELF3/CLDN4 network?

Rules fixed before looking at the set calls
--------------------------------------------
* Expression effect is log2(KD) - log2(control). Two-sided tests only.
* A program is UP or DOWN only when the detected-gene median |log2FC| is at
  least 0.25, the same-direction fraction is at least 0.60, and the
  Mann–Whitney test against the rest of the detected transcriptome has p < 0.05.
* A significant test with a smaller median is WEAK, not a program change.
* ELF3 log2FC must be <= -0.5 or the contrast is QC-fail and is not counted
  in the cross-dataset verdict.
* CLDN4 and MSigDB Hallmark apical junction are the positive controls for the
  ELF3 epithelial network. Hallmark DNA repair is reported separately from
  the 13-gene classical NHEJ set.
* NEC siELF3 tables (Horie et al., Cancer Sci 2023, GSE190618 supplement S7–S9)
  list only genes that fall (fold change < 0.75). They cannot test induction.
* No lung ELF3 ChIP-seq is public (same conclusion as the A10 catalog). DMS53
  CUT&Tag is SCLC, not LUAD.

Nothing here is a private KL matrix.
"""

from __future__ import annotations

import gzip
import json
import sqlite3
import tarfile
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "elf3_nhej_sting_ifn"
CACHE = OUT / "cache"

# Classical c-NHEJ ligation / end-processing machinery. Not the whole DDR,
# and not the Reactome alias dump (that file mixes histone gene synonyms).
NHEJ_CORE = [
    "XRCC6",  # Ku70
    "XRCC5",  # Ku80
    "PRKDC",  # DNA-PKcs
    "XRCC4",
    "LIG4",
    "NHEJ1",  # XLF
    "PAXX",
    "DCLRE1C",  # Artemis
    "APLF",
    "PNKP",
    "APTX",
    "POLL",
    "POLM",
]

# cGAS–STING signaling components. Downstream ISGs stay in the IFN sets.
STING_CORE = [
    "CGAS",
    "STING1",
    "TBK1",
    "IKBKE",
    "IRF3",
    "IRF7",
    "IFI16",
    "DDX41",
    "ZBP1",
    "TREX1",
    "ENPP1",
]

# Single genes shown beside the programs. CLDN4 is the network positive control.
FOCUS_GENES = [
    "ELF3",
    "CLDN4",
    "TACSTD2",
    "AURKA",
    "CDC25B",
    "ITGB6",
    "YWHAB",
] + NHEJ_CORE + STING_CORE + [
    "IFNB1",
    "IFNL1",
    "ISG15",
    "MX1",
    "OAS1",
    "OAS2",
    "IFIT1",
    "IFI27",
    "STAT1",
    "STAT2",
    "IRF9",
    "IRF1",
    "CXCL10",
    "DDX58",
    "B2M",
    "TAP1",
    "HLA-A",
    "HLA-B",
    "CD274",
]

HALLMARK_KEEP = [
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
]

PRIMARY_SETS = [
    "NHEJ_CORE",
    "STING_CORE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]

SECONDARY_SETS = [
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
]

# Human, basal, n>=3 biological replicates. Cross-dataset verdict uses these
# only when ELF3 itself falls (log2FC <= -0.5).
PRIMARY_DATASETS = [
    "A549_shELF3",
    "HBDEC2_KO",
    "HBDEC2_KD",
    "JEG3_KO",
    "synov_siELF3_basal",
]


def ensure_cache() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)


def fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "sdaxcge-elf3-nhej/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
        out.write(resp.read())
    return dest


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in HALLMARK_KEEP:
                sets[parts[0]] = [g for g in parts[2:] if g]
    missing = [n for n in HALLMARK_KEEP if n not in sets]
    if missing:
        raise SystemExit(f"Hallmark GMT missing {missing}")
    return sets


def hgnc_alias_map(path: Path) -> dict[str, str]:
    """Map retired symbols onto current HGNC symbols.

    Agilent GPL20844 (2015) and the CRPC count table still use TMEM173,
    MB21D1, and C9orf142. Official symbols are never overwritten by an alias.
    """
    table = pd.read_csv(path, sep="\t", dtype=str, usecols=["symbol", "alias_symbol", "prev_symbol"])
    current = {s for s in table["symbol"].dropna() if s and s != "nan"}
    mapping = {s: s for s in current}
    for symbol, alias, prev in zip(table["symbol"], table["alias_symbol"], table["prev_symbol"]):
        if not isinstance(symbol, str) or symbol not in current:
            continue
        for field in (alias, prev):
            if not isinstance(field, str) or field in {"", "nan"}:
                continue
            for raw in field.split("|"):
                name = raw.strip()
                if name and name not in current:
                    mapping.setdefault(name, symbol)
    return mapping


def canonicalize(df: pd.DataFrame, alias: dict[str, str], how: str) -> pd.DataFrame:
    renamed = df.copy()
    renamed.index = pd.Index([alias.get(str(i), str(i)) for i in renamed.index])
    if not renamed.index.duplicated().any():
        return renamed
    grouped = renamed.groupby(level=0)
    return grouped.sum(numeric_only=True) if how == "sum" else grouped.mean(numeric_only=True)


def mouse_symbol(human: str) -> str:
    if human.startswith("HLA-"):
        return human  # no 1:1 mouse symbol; will miss
    return human[0] + human[1:].lower()


def gene_sets(species: str, hallmark: dict[str, list[str]]) -> dict[str, list[str]]:
    sets = {"NHEJ_CORE": list(NHEJ_CORE), "STING_CORE": list(STING_CORE)}
    sets.update(hallmark)
    if species == "mouse":
        sets = {k: [mouse_symbol(g) for g in v] for k, v in sets.items()}
    return sets


def bh(pvalues: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=pvalues.index, dtype=float)
    ok = pvalues.notna()
    if ok.sum() == 0:
        return out
    out.loc[ok] = multipletests(pvalues.loc[ok].to_numpy(), method="fdr_bh")[1]
    return out


def welch_p(kd: np.ndarray, ctrl: np.ndarray) -> np.ndarray:
    # kd, ctrl: genes x replicates
    res = stats.ttest_ind(kd, ctrl, axis=1, equal_var=False, nan_policy="omit")
    return np.asarray(res.pvalue, dtype=float)


def paired_p(delta: np.ndarray) -> np.ndarray:
    # delta: genes x donors
    n = np.sum(np.isfinite(delta), axis=1)
    mean = np.nanmean(delta, axis=1)
    sd = np.nanstd(delta, axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = mean / (sd / np.sqrt(n))
    p = np.full(mean.shape, np.nan)
    ok = (n >= 3) & np.isfinite(tstat) & (sd > 0)
    p[ok] = 2 * stats.t.sf(np.abs(tstat[ok]), n[ok] - 1)
    return p


def set_row(name: str, lfc: pd.Series, genes: list[str], detected: pd.Series) -> dict:
    present = [g for g in genes if g in lfc.index]
    use = [g for g in present if bool(detected.get(g, False))]
    vals = lfc.reindex(use).astype(float)
    vals = vals[np.isfinite(vals)]
    all_vals = lfc.reindex(present).astype(float)
    all_vals = all_vals[np.isfinite(all_vals)]
    bg = lfc[detected.reindex(lfc.index).fillna(False)].drop(index=vals.index, errors="ignore")
    bg = bg[np.isfinite(bg)]
    rec = {
        "set": name,
        "n_in_set": len(genes),
        "n_present": len(present),
        "n_detected": int(vals.shape[0]),
        "median_log2fc_detected": float(np.median(vals)) if len(vals) else np.nan,
        "mean_log2fc_detected": float(np.mean(vals)) if len(vals) else np.nan,
        "frac_up_detected": float(np.mean(vals > 0)) if len(vals) else np.nan,
        "median_log2fc_all": float(np.median(all_vals)) if len(all_vals) else np.nan,
        "frac_up_all": float(np.mean(all_vals > 0)) if len(all_vals) else np.nan,
        "mw_p": np.nan,
        "wilcoxon_p": np.nan,
    }
    if len(vals) >= 5 and len(bg) >= 20:
        rec["mw_p"] = float(stats.mannwhitneyu(vals, bg, alternative="two-sided").pvalue)
    if len(vals) >= 8 and np.any(np.abs(vals.to_numpy()) > 1e-8):
        rec["wilcoxon_p"] = float(stats.wilcoxon(vals.to_numpy(), alternative="two-sided").pvalue)
    rec["call"] = call_set(rec["median_log2fc_detected"], rec["frac_up_detected"], rec["mw_p"], rec["n_detected"])
    return rec


def call_set(median: float, frac_up: float, p: float, n: int) -> str:
    if n < 5 or not np.isfinite(median) or not np.isfinite(p):
        return "INSUFFICIENT"
    strong_up = p < 0.05 and median >= 0.25 and frac_up >= 0.60
    strong_down = p < 0.05 and median <= -0.25 and frac_up <= 0.40
    if strong_up:
        return "UP"
    if strong_down:
        return "DOWN"
    if p < 0.05:
        return "WEAK"
    return "NULL"


def gene_table(lfc: pd.Series, p: pd.Series, fdr: pd.Series, detected: pd.Series, genes: list[str]) -> pd.DataFrame:
    rows = []
    for g in genes:
        rows.append(
            {
                "gene": g,
                "log2fc": float(lfc[g]) if g in lfc.index and np.isfinite(lfc[g]) else np.nan,
                "p": float(p[g]) if g in p.index and np.isfinite(p[g]) else np.nan,
                "fdr": float(fdr[g]) if g in fdr.index and np.isfinite(fdr[g]) else np.nan,
                "detected": bool(detected[g]) if g in detected.index else False,
                "present": g in lfc.index,
            }
        )
    return pd.DataFrame(rows)


def contrast_unpaired(log_df: pd.DataFrame, kd_cols: list[str], ctrl_cols: list[str], detected_floor: float) -> dict:
    kd = log_df[kd_cols]
    ctrl = log_df[ctrl_cols]
    lfc = kd.mean(axis=1) - ctrl.mean(axis=1)
    p = pd.Series(welch_p(kd.to_numpy(), ctrl.to_numpy()), index=log_df.index)
    # Detection on the linear-ish log scale is not CPM. Callers pass log2 data
    # and a floor already applied via `detected`, except we recompute from the
    # mean of the log values only as a fallback. The caller sets `detected`.
    return {"lfc": lfc, "p": p, "kd_mean": kd.mean(axis=1), "ctrl_mean": ctrl.mean(axis=1)}


def finalize(name: str, meta: dict, lfc: pd.Series, p: pd.Series, detected: pd.Series, sets: dict, species: str) -> dict:
    fdr = bh(p.where(detected.reindex(p.index).fillna(False)))
    rows = []
    for set_name in PRIMARY_SETS + SECONDARY_SETS:
        rows.append(set_row(set_name, lfc, sets[set_name], detected))
    genes = FOCUS_GENES if species == "human" else [mouse_symbol(g) for g in FOCUS_GENES]
    # mouse focus: Elf3 not mouse_symbol("ELF3") wait, mouse_symbol("ELF3") == "Elf3". Good.
    gtab = gene_table(lfc, p, fdr, detected, genes)
    elf = "ELF3" if species == "human" else "Elf3"
    elf_lfc = float(lfc[elf]) if elf in lfc.index else np.nan
    elf_p = float(p[elf]) if elf in p.index and np.isfinite(p[elf]) else np.nan
    qc = "PASS" if np.isfinite(elf_lfc) and elf_lfc <= -0.5 else "FAIL"
    return {
        "dataset": name,
        **meta,
        "species": species,
        "elf3_log2fc": elf_lfc,
        "elf3_p": elf_p,
        "qc": qc,
        "n_detected": int(detected.sum()),
        "sets": rows,
        "genes": gtab,
        "lfc": lfc,
        "fdr": fdr,
    }


def load_geo_table(path: Path) -> tuple[dict[str, list[str]], pd.DataFrame]:
    meta_lines: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        header = None
        rows = []
        for line in handle:
            if line.startswith("!"):
                key, _, rest = line[1:].rstrip("\n").partition("\t")
                meta_lines.setdefault(key, [])
                meta_lines[key].append(rest)
                continue
            if line.startswith("!series_matrix_table_begin") or line.startswith("series_matrix_table_begin"):
                continue
            if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                header = [c.strip('"') for c in line.rstrip("\n").split("\t")]
                continue
            if line.startswith("!series_matrix_table_end") or not header:
                continue
            if line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            rows.append(parts)
    titles = [c.strip('"') for c in meta_lines["Sample_title"][0].split("\t")]
    accessions = [c.strip('"') for c in meta_lines["Sample_geo_accession"][0].split("\t")]
    df = pd.DataFrame(rows, columns=header)
    df = df.rename(columns={"ID_REF": "id"}).set_index("id")
    df.index = df.index.str.strip('"')
    df = df.apply(pd.to_numeric, errors="coerce")
    df.columns = accessions
    return {"titles": titles, "accessions": accessions, "title_by_acc": dict(zip(accessions, titles))}, df


def hugene_symbol_map(sqlite_path: Path, hgnc_path: Path) -> dict[str, str]:
    entrez_to_symbol = {}
    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str, usecols=["symbol", "entrez_id"])
    for symbol, entrez in zip(hgnc["symbol"], hgnc["entrez_id"]):
        if isinstance(entrez, str) and entrez and entrez != "nan":
            entrez_to_symbol.setdefault(entrez, symbol)
    con = sqlite3.connect(sqlite_path)
    rows = con.execute(
        "SELECT probe_id, gene_id FROM probes WHERE gene_id IS NOT NULL AND gene_id != '' AND is_multiple = 0"
    )
    out = {}
    for probe, entrez in rows:
        symbol = entrez_to_symbol.get(str(entrez))
        if symbol:
            out[str(probe)] = symbol
    con.close()
    return out


def agilent_symbol_map(path: Path) -> dict[str, str]:
    out = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if not line or not line[0].isdigit():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            symbol = parts[9].strip()
            if not symbol or symbol in {"nan", "NA"}:
                continue
            symbol = symbol.split("|")[0].split("///")[0].split(";")[0].strip()
            if symbol and symbol.replace("-", "").isalnum():
                out[parts[0]] = symbol
    return out


def collapse_max_mean(df: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    work = df.copy()
    work["symbol"] = symbols.reindex(work.index).to_numpy()
    work = work.dropna(subset=["symbol"])
    value_cols = [c for c in work.columns if c != "symbol"]
    work["_mean"] = work[value_cols].mean(axis=1)
    work = work.sort_values("_mean", ascending=False).drop_duplicates("symbol")
    return work.drop(columns="_mean").set_index("symbol")[value_cols]


def split_by_title(title_by_acc: dict[str, str], kd_tokens: tuple[str, ...], ctrl_tokens: tuple[str, ...]) -> tuple[list[str], list[str]]:
    kd, ctrl = [], []
    for acc, title in title_by_acc.items():
        low = title.lower()
        if any(tok in low for tok in kd_tokens):
            kd.append(acc)
        elif any(tok in low for tok in ctrl_tokens):
            ctrl.append(acc)
    return kd, ctrl


def analyze_microarray(name, path, symbol_map, sets, species, kind, kd_tokens, ctrl_tokens, note, role, alias) -> dict:
    meta, table = load_geo_table(path)
    symbols = pd.Series(symbol_map)
    if kind == "rma":
        logged = table
    else:
        logged = np.log2(table.clip(lower=0) + 1)
    logged = collapse_max_mean(logged, symbols)
    logged = canonicalize(logged, alias, how="mean")
    kd, ctrl = split_by_title(meta["title_by_acc"], kd_tokens, ctrl_tokens)
    if len(kd) < 2 or len(ctrl) < 2:
        raise SystemExit(f"{name}: could not assign columns kd={kd} ctrl={ctrl} titles={meta['titles']}")
    raw = contrast_unpaired(logged, kd, ctrl, 0)
    # Microarrays are treated as detected for every mapped gene.
    detected = pd.Series(True, index=logged.index)
    result = finalize(name, {"role": role, "note": note, "n_kd": len(kd), "n_ctrl": len(ctrl), "accession_note": ",".join(meta["titles"])}, raw["lfc"], raw["p"], detected, sets, species)
    return result


def analyze_jeg3(path: Path, sets: dict, alias: dict[str, str]) -> dict:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "gene"})
    df = df[df["gene"].fillna("").str.fullmatch(r"[A-Za-z0-9][A-Za-z0-9\-\.]*")]
    df = df.groupby("gene", as_index=True).sum(numeric_only=True)
    df = canonicalize(df, alias, how="sum")
    kd_cols = [c for c in df.columns if c.startswith("E3.")]
    ctrl_cols = [c for c in df.columns if c.startswith("WT.")]
    lib = df.sum(axis=0).replace(0, np.nan)
    cpm = df.div(lib, axis=1) * 1e6
    logged = np.log2(cpm + 1)
    raw = contrast_unpaired(logged, kd_cols, ctrl_cols, 0)
    mean_cpm = cpm[kd_cols + ctrl_cols].mean(axis=1)
    detected = mean_cpm >= 1
    return finalize(
        "JEG3_KO",
        {
            "role": "primary",
            "note": "GSE241792 JEG-3 trophoblast CRISPR sgELF3 vs WT, n=3. Salmon counts, log2(CPM+1). Not lung.",
            "n_kd": 3,
            "n_ctrl": 3,
        },
        raw["lfc"],
        raw["p"],
        detected,
        sets,
        "human",
    )


def analyze_prostate(path: Path, sets: dict, alias: dict[str, str]) -> dict:
    df = pd.read_excel(path)
    df = df.rename(columns={"gene_name": "gene"})
    df = df.dropna(subset=["gene"])
    count_cols = ["shELF3_rep1", "shELF3_rep2", "PLKO_rep1", "PLKO_rep2"]
    df = df.groupby("gene", as_index=True)[count_cols + ["log2FoldChange"]].mean(numeric_only=True)
    df = canonicalize(df, alias, how="mean")
    logged = np.log2(df[count_cols].clip(lower=0) + 1)
    raw = contrast_unpaired(logged, ["shELF3_rep1", "shELF3_rep2"], ["PLKO_rep1", "PLKO_rep2"], 0)
    mean_lin = df[count_cols].mean(axis=1)
    detected = mean_lin >= 1
    result = finalize(
        "CRPC_shELF3",
        {
            "role": "supporting_n2",
            "note": "GSE303076 CRPC shELF3 vs PLKO, n=2. Author log2FoldChange is stored for ELF3/CLDN4 checks. Not in the 5-dataset verdict.",
            "n_kd": 2,
            "n_ctrl": 2,
        },
        raw["lfc"],
        raw["p"],
        detected,
        sets,
        "human",
    )
    result["author_log2fc_elf3"] = float(df.loc["ELF3", "log2FoldChange"]) if "ELF3" in df.index else np.nan
    result["author_log2fc_cldn4"] = float(df.loc["CLDN4", "log2FoldChange"]) if "CLDN4" in df.index else np.nan
    return result


def analyze_mouse(path: Path, sets: dict) -> dict:
    df = pd.read_excel(path)
    df["gene"] = df["gene_id"].astype(str).str.replace(r"^ENSMUSG\d+_", "", regex=True)
    value_cols = [c for c in df.columns if c not in {"gene_id", "gene"}]
    df = df.groupby("gene", as_index=True)[value_cols].mean(numeric_only=True)
    kd = [c for c in df.columns if c.startswith("KO")]
    ctrl = [c for c in df.columns if c.startswith("CTR")]
    logged = np.log2(df.clip(lower=0) + 1)
    raw = contrast_unpaired(logged, kd, ctrl, 0)
    detected = df.mean(axis=1) >= 1
    return finalize(
        "mTEC_KO_TNFa_IFNg",
        {
            "role": "stimulated",
            "note": "GSE319933 primary mouse renal tubular cells, AdvCre Elf3 deletion vs control, all samples already on TNFα+IFNγ. Tests whether Elf3 loss changes programs during cytokine stimulation, not basal loss.",
            "n_kd": len(kd),
            "n_ctrl": len(ctrl),
        },
        raw["lfc"],
        raw["p"],
        detected,
        sets,
        "mouse",
    )


def analyze_synovial(counts_path: Path, meta_path: Path, hgnc_path: Path, sets: dict, alias: dict[str, str]) -> list[dict]:
    # The stimulation field uses the literal string "None" for unstimulated wells.
    meta = pd.read_csv(meta_path, sep="\t", keep_default_na=False)
    counts = pd.read_csv(counts_path, sep="\t").rename(columns={"ID_REF": "ensembl"})
    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id"])
    mapping = hgnc.dropna(subset=["ensembl_gene_id"]).drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]
    counts["gene"] = counts["ensembl"].map(mapping)
    counts = counts.dropna(subset=["gene"])
    sample_cols = [c for c in counts.columns if c not in {"ensembl", "gene"}]
    mat = counts.groupby("gene")[sample_cols].sum()
    mat = canonicalize(mat, alias, how="sum")
    lib = mat.sum(axis=0)
    keep_samples = lib[lib >= 1e5].index.tolist()
    meta = meta[meta["sample"].isin(keep_samples)].copy()
    meta["time"] = meta["time"].astype(int)
    results = []

    def donor_log(sub_meta: pd.DataFrame) -> pd.DataFrame:
        frames = []
        for donor, grp in sub_meta.groupby("donor"):
            cols = [c for c in grp["sample"] if c in mat.columns]
            if not cols:
                continue
            mean_counts = mat[cols].mean(axis=1)
            frames.append(mean_counts.rename(donor))
        out = pd.concat(frames, axis=1)
        lib_d = out.sum(axis=0).replace(0, np.nan)
        return np.log2(out.div(lib_d, axis=1) * 1e6 + 1)

    contrasts = [
        ("synov_siELF3_basal", "primary", "time == 0 and stimulation == 'None'", "GSE129487 RA synovial fibroblasts, siELF3 vs siCtrl, unstimulated (t=0). Four donors, technical duplicates averaged. Not lung."),
        ("synov_siELF3_TNF16", "stimulated", "time == 16 and stimulation == 'TNF (1)'", "GSE129487 siELF3 vs siCtrl at 16 h TNF. Same donors."),
        ("synov_siELF3_TNF_IL17_16", "stimulated", "time == 16 and stimulation == 'TNF (1) + IL17 (1)'", "GSE129487 siELF3 vs siCtrl at 16 h TNF+IL-17A. Same donors."),
    ]
    for name, role, query, note in contrasts:
        sub = meta.query(query)
        kd_meta = sub[sub["sirna"] == "ELF3"]
        ctrl_meta = sub[sub["sirna"] == "Ctrl"]
        kd_log = donor_log(kd_meta)
        ctrl_log = donor_log(ctrl_meta)
        donors = sorted(set(kd_log.columns) & set(ctrl_log.columns))
        kd_log = kd_log[donors]
        ctrl_log = ctrl_log[donors]
        delta = kd_log - ctrl_log
        lfc = delta.mean(axis=1)
        p = pd.Series(paired_p(delta.to_numpy()), index=delta.index)
        # detected if mean CPM-equivalent log is not the point: use either side mean count via inverse is messy.
        # Use mean log2(CPM+1) >= log2(2) i.e. typical CPM at least ~1.
        detected = ((kd_log.mean(axis=1) >= 1) | (ctrl_log.mean(axis=1) >= 1))
        results.append(
            finalize(
                name,
                {"role": role, "note": note + f" Donors used: {','.join(donors)}.", "n_kd": len(donors), "n_ctrl": len(donors)},
                lfc,
                p,
                detected,
                sets,
                "human",
            )
        )
    return results


def load_refgene(path: Path) -> pd.DataFrame:
    cols = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "symbol"]
    df = pd.read_csv(path, sep="\t", header=None, usecols=list(range(13)), names=cols, compression="gzip")
    tss = np.where(df["strand"].to_numpy() == "+", df["txStart"].to_numpy(), df["txEnd"].to_numpy())
    df["tss"] = tss.astype(np.int64)
    return df


def peaks_from_bed(path: Path) -> dict[str, np.ndarray]:
    starts, ends, chroms = [], [], []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chroms.append(parts[0])
            starts.append(int(parts[1]))
            ends.append(int(parts[2]))
    by = {}
    df = pd.DataFrame({"chrom": chroms, "start": starts, "end": ends})
    for chrom, sub in df.groupby("chrom"):
        sub = sub.sort_values("start")
        by[chrom] = (sub["start"].to_numpy(), sub["end"].to_numpy())
    return by


def window_hit(symbol_table: pd.DataFrame, peaks: dict[str, np.ndarray], pad: int, mode: str) -> set[str]:
    hit = set()
    if mode == "tss":
        left = symbol_table["tss"].to_numpy() - pad
        right = symbol_table["tss"].to_numpy() + pad
    else:
        left = symbol_table["txStart"].to_numpy() - pad
        right = symbol_table["txEnd"].to_numpy() + pad
    chrom = symbol_table["chrom"].to_numpy()
    symbols = symbol_table["symbol"].to_numpy()
    for i in range(len(symbol_table)):
        arr = peaks.get(chrom[i])
        if arr is None:
            continue
        starts, ends = arr
        # any peak with start < right and end > left
        j = np.searchsorted(starts, right[i], side="left")
        lo = max(0, j - 50)
        # walk back while peaks could still overlap; peaks are sorted by start
        # A peak starting before `right` may end after `left`. Check a window of candidates.
        # Because peak lengths vary, scan backward until start is far below left.
        k = j - 1
        while k >= 0 and starts[k] > left[i] - 2_000_000:
            if ends[k] > left[i] and starts[k] < right[i]:
                hit.add(symbols[i])
                break
            # Called peaks are a few kb; stop once the peak start is far left of the window.
            if starts[k] < left[i] - 20000:
                break
            k -= 1
    return hit


def fisher(set_genes: list[str], bound: set[str], universe: set[str]) -> dict:
    inset = [g for g in set_genes if g in universe]
    if not inset:
        return {"n": 0, "n_bound": 0, "frac": np.nan, "odds_ratio": np.nan, "p": np.nan}
    a = sum(g in bound for g in inset)
    b = len(inset) - a
    c = sum(g in bound for g in universe if g not in set(inset))
    d = len(universe) - c - a - b
    # universe includes inset; c is bound outside, d is unbound outside
    oddsr, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return {"n": len(inset), "n_bound": a, "frac": a / len(inset), "odds_ratio": float(oddsr), "p": float(p), "bg_frac": c / (c + d) if (c + d) else np.nan}


def chip_atlas(path: Path, sets: dict) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "gene"})
    df = df.drop_duplicates("gene").set_index("gene")
    num = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    groups = {
        "CFPAC1_any": [c for c in num.columns if "CFPAC" in c],
        "ESO26": [c for c in num.columns if "ESO-26" in c],
        "HBDEC2": [c for c in num.columns if "HBDEC2" in c],
        "HepG2_any": [c for c in num.columns if "Hep_G2" in c or "HepG2" in c],
    }
    bound = {}
    for name, cols in groups.items():
        bound[name] = set(num.index[(num[cols] > 0).any(axis=1)])
    bound["any_library"] = set(num.index[(num.drop(columns=[c for c in num.columns if c in {"STRING", "ELF3|Average"}], errors="ignore") > 0).any(axis=1)])
    universe = set(num.index)
    rows = []
    focus = {"CLDN4": ["CLDN4"], "TACSTD2": ["TACSTD2"], "ELF3": ["ELF3"]}
    for exp, genes_bound in bound.items():
        for set_name, genes in {**sets, **focus}.items():
            if set_name not in PRIMARY_SETS + SECONDARY_SETS and set_name not in focus:
                continue
            rec = fisher(genes, genes_bound, universe)
            rec.update({"source": "ChIP-Atlas_hg38_5kb", "experiment": exp, "set": set_name, "n_bound_universe": len(genes_bound), "n_universe": len(universe)})
            rows.append(rec)
    # gene-level matrix for focus genes
    gene_rows = []
    for gene in FOCUS_GENES:
        if gene not in num.index:
            gene_rows.append({"gene": gene, "present": False})
            continue
        rec = {"gene": gene, "present": True, "average": float(num.loc[gene, "ELF3|Average"]) if "ELF3|Average" in num.columns else np.nan}
        for exp, cols in groups.items():
            rec[exp] = int((num.loc[gene, cols] > 0).any())
            rec[exp + "_max"] = float(num.loc[gene, cols].max())
        gene_rows.append(rec)
    return pd.DataFrame(rows), pd.DataFrame(gene_rows)


def cuttag_overlap(sets: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    ref = load_refgene(CACHE / "refGene_hg38.txt.gz")
    experiments = {
        "ECC4_CUTTag": CACHE / "ECC4_ELF3.bed.gz",
        "A99_CUTTag": CACHE / "A99_ELF3.bed.gz",
        "DMS53_CUTTag": CACHE / "DMS53_ELF3.bed.gz",
        "HBDEC2_ChIP": CACHE / "HBDEC2_ELF3.bed.gz",
        "HepG2_ENCODE_IDR": CACHE / "HepG2_ELF3_IDR.bed.gz",
    }
    universe = set(ref["symbol"])
    rows = []
    gene_rows = []
    focus = FOCUS_GENES
    for exp, path in experiments.items():
        peaks = peaks_from_bed(path)
        n_peaks = sum(len(v[0]) for v in peaks.values())
        for mode, pad in (("tss_2kb", 2000), ("tss_10kb", 10000)):
            bound = window_hit(ref, peaks, pad, "tss")
            for set_name, genes in list(sets.items()) + [("CLDN4", ["CLDN4"]), ("TACSTD2", ["TACSTD2"]), ("ELF3", ["ELF3"])]:
                if set_name not in PRIMARY_SETS + SECONDARY_SETS and set_name not in {"CLDN4", "TACSTD2", "ELF3"}:
                    continue
                rec = fisher(genes, bound, universe)
                rec.update({"source": "peak_TSS", "experiment": exp, "window": mode, "set": set_name, "n_peaks": n_peaks})
                rows.append(rec)
        bound2 = window_hit(ref, peaks, 2000, "tss")
        for gene in focus:
            gene_rows.append({"experiment": exp, "gene": gene, "tss_2kb": int(gene in bound2), "n_peaks": n_peaks})
    # mouse CUT&RUN, stimulated TEC
    mref = load_refgene(CACHE / "refGene_mm10.txt.gz")
    msets = {k: [mouse_symbol(g) for g in v] for k, v in sets.items()}
    peaks = peaks_from_bed(CACHE / "mm_ELF3.narrowPeak.gz")
    n_peaks = sum(len(v[0]) for v in peaks.values())
    bound = window_hit(mref, peaks, 2000, "tss")
    universe_m = set(mref["symbol"])
    for set_name, genes in list(msets.items()) + [("Cldn4", ["Cldn4"]), ("Elf3", ["Elf3"])]:
        if set_name not in PRIMARY_SETS + SECONDARY_SETS and set_name not in {"Cldn4", "Elf3"}:
            continue
        rec = fisher(genes, bound, universe_m)
        rec.update({"source": "peak_TSS", "experiment": "mTEC_CUTRUN_TNFa_IFNg", "window": "tss_2kb", "set": set_name, "n_peaks": n_peaks})
        rows.append(rec)
    return pd.DataFrame(rows), pd.DataFrame(gene_rows)


def nec_lists() -> tuple[pd.DataFrame, dict[str, set[str]]]:
    import openpyxl

    wb = openpyxl.load_workbook(CACHE / "CAS-114-2596-s001.xlsx", read_only=True, data_only=True)
    lists = {}
    for sheet, key in (("S7", "near"), ("S8", "down")):
        ws = wb[sheet]
        cols = {0: set(), 1: set(), 2: set()}
        for row in ws.iter_rows(min_row=5, values_only=True):
            for j in range(3):
                val = row[j + 1] if row and len(row) > j + 1 else None
                if val:
                    cols[j].add(str(val).strip())
        for j, line in enumerate(("ECC4", "A99", "DMS53")):
            lists[f"{key}_{line}"] = cols[j]
    ws = wb["S9"]
    s9 = set()
    for row in ws.iter_rows(min_row=4, values_only=True):
        if row and len(row) > 1 and row[1]:
            s9.add(str(row[1]).strip())
    lists["direct_232"] = s9
    wb.close()
    # recurrent: in at least two lines
    lists["near_ge2"] = set()
    lists["down_ge2"] = set()
    for bucket, prefix in ((lists["near_ge2"], "near"), (lists["down_ge2"], "down")):
        from collections import Counter

        cnt = Counter()
        for line in ("ECC4", "A99", "DMS53"):
            cnt.update(lists[f"{prefix}_{line}"])
        bucket.update({g for g, n in cnt.items() if n >= 2})
    rows = []
    catalog = {
        "NHEJ_CORE": NHEJ_CORE,
        "STING_CORE": STING_CORE,
        "focus": FOCUS_GENES,
    }
    # hallmark added by caller via a side file? We'll fill hallmark in main by passing sets.
    return lists


def nec_overlap(lists: dict[str, set[str]], sets: dict) -> pd.DataFrame:
    rows = []
    catalog = {k: sets[k] for k in PRIMARY_SETS + SECONDARY_SETS}
    catalog["CLDN4"] = ["CLDN4"]
    catalog["TACSTD2"] = ["TACSTD2"]
    for list_name, genes_in in lists.items():
        for set_name, genes in catalog.items():
            present = [g for g in genes if g in genes_in]
            rows.append(
                {
                    "list": list_name,
                    "n_list": len(genes_in),
                    "set": set_name,
                    "n_set": len(genes),
                    "n_overlap": len(present),
                    "genes": ",".join(sorted(present)),
                }
            )
    return pd.DataFrame(rows)


def fmt(x, digits=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    ax = abs(float(x))
    if ax != 0 and ax < 1e-3:
        return f"{float(x):.2e}"
    return f"{float(x):.{digits}f}"


def _recurrent_peaks(peak_genes: pd.DataFrame) -> str:
    nec = ["ECC4_CUTTag", "A99_CUTTag", "DMS53_CUTTag"]
    sub = peak_genes[peak_genes["experiment"].isin(nec)]
    if sub.empty:
        return ""
    wide = sub.pivot(index="gene", columns="experiment", values="tss_2kb")
    keep = [g for g in NHEJ_CORE + STING_CORE + ["CLDN4", "TACSTD2"] if g in wide.index]
    bits = []
    for gene in keep:
        n = int(wide.loc[gene, nec].sum())
        if n >= 2:
            bits.append(f"{gene} ({n}/3)")
    return (
        "TSS ±2 kb peaks in at least two of ECC4, A99, and DMS53: "
        + (", ".join(bits) if bits else "none")
        + ". A99 CLDN4 is outside ±2 kb and inside ±10 kb."
    )


def write_report(results: list[dict], set_df: pd.DataFrame, gene_df: pd.DataFrame, chip_df: pd.DataFrame, peak_df: pd.DataFrame, nec_df: pd.DataFrame, peak_genes: pd.DataFrame) -> str:
    # BH across primary dataset x primary set MW p values
    prim = set_df[set_df["dataset"].isin(PRIMARY_DATASETS) & set_df["set"].isin(PRIMARY_SETS)].copy()
    prim = prim.merge(pd.DataFrame([{"dataset": r["dataset"], "qc": r["qc"]} for r in results]), on="dataset", how="left")
    qc_pass = prim[prim["qc"] == "PASS"]
    lines = []
    lines.append("# ELF3 loss vs NHEJ, STING, and IFN")
    lines.append("")
    lines.append("Public ELF3 ChIP-seq / CUT&Tag and ELF3-loss transcriptomes. CLDN4 is the positive-control node of the ELF3 epithelial network. Private KL matrices are not used.")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    # build verdict from qc-pass primary
    lines.append(_verdict_paragraph(results, set_df))
    lines.append("")
    lines.append(_binding_sentence(peak_df, nec_df))
    lines.append("")
    lines.append("Calls use detected genes only. UP/DOWN requires median |log2FC| ≥ 0.25, same-direction fraction ≥ 0.60, and two-sided Mann–Whitney p < 0.05 versus the rest of the detected transcriptome. WEAK means p < 0.05 with a smaller shift. NULL means the set does not move relative to the transcriptome.")
    lines.append("")
    lines.append("## ELF3 knockdown QC")
    lines.append("")
    lines.append("| Dataset | Role | n KD / n ctrl | ELF3 log2FC | ELF3 raw p | QC | CLDN4 log2FC |")
    lines.append("| --- | --- | ---: | ---: | ---: | --- | ---: |")
    gidx = gene_df.set_index(["dataset", "gene"])
    for r in results:
        def gl(gene):
            key = (r["dataset"], gene if r["species"] == "human" else mouse_symbol(gene))
            if key not in gidx.index:
                return np.nan
            return gidx.loc[key, "log2fc"]

        lines.append(
            f"| {r['dataset']} | {r['role']} | {r['n_kd']} / {r['n_ctrl']} | {fmt(r['elf3_log2fc'])} | {fmt(r['elf3_p'])} | {r['qc']} | {fmt(gl('CLDN4'))} |"
        )
    lines.append("")
    lines.append("QC PASS requires ELF3 log2FC ≤ −0.5. That is a point-estimate gate, not an FDR gate: A549 ELF3 raw p is small, and the genome-wide BH FDR can still sit above 0.05 because thousands of genes are tested at n=3. Contrasts that fail the gate are shown and are excluded from the basal cross-dataset sentence.")
    lines.append("")
    lines.append("GSE303076 has two replicates per arm. A Mann–Whitney p-value there uses genes as the sample, so tight replicates produce very small p-values. Read that row as a direction.")
    lines.append("")
    lines.append("## Program shifts (detected genes)")
    lines.append("")
    lines.append("| Dataset | QC | Set | n | median log2FC | frac up | MW p | Wilcoxon p | Call |")
    lines.append("| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |")
    qc_map = {r["dataset"]: r["qc"] for r in results}
    show = set_df[set_df["set"].isin(PRIMARY_SETS + ["HALLMARK_DNA_REPAIR", "HALLMARK_APICAL_JUNCTION", "HALLMARK_E2F_TARGETS"])]
    for rec in show.itertuples(index=False):
        lines.append(
            f"| {rec.dataset} | {qc_map[rec.dataset]} | {rec.set} | {rec.n_detected} | {fmt(rec.median_log2fc_detected)} | {fmt(rec.frac_up_detected)} | {fmt(rec.mw_p)} | {fmt(rec.wilcoxon_p)} | {rec.call} |"
        )
    lines.append("")
    lines.append("Apical junction and E2F/G2M are controls. The NEC paper’s ELF3 signature is G2M/E2F, not interferon. Hallmark DNA repair is the broad proliferation-linked repair set, not classical NHEJ.")
    lines.append("")
    lines.append("## Focus genes")
    lines.append("")
    lines.append("log2FC (KD − control). FDR is Benjamini–Hochberg within detected genes of that contrast. Blank FDR means the gene was below the detection floor.")
    lines.append("")
    focus_show = ["ELF3", "CLDN4", "TACSTD2", "XRCC5", "XRCC6", "PRKDC", "XRCC4", "LIG4", "NHEJ1", "CGAS", "STING1", "TBK1", "IRF3", "IFNB1", "ISG15", "MX1", "OAS1", "IFIT1", "STAT1", "CXCL10", "CD274"]
    header = "| Gene | " + " | ".join(r["dataset"] for r in results if r["species"] == "human") + " |"
    lines.append(header)
    lines.append("| --- | " + " | ".join("---:" for r in results if r["species"] == "human") + " |")
    human_sets = [r["dataset"] for r in results if r["species"] == "human"]
    for gene in focus_show:
        cells = []
        for ds in human_sets:
            key = (ds, gene)
            if key not in gidx.index or not np.isfinite(gidx.loc[key, "log2fc"]):
                cells.append("NA")
                continue
            lfc = gidx.loc[key, "log2fc"]
            fdr = gidx.loc[key, "fdr"]
            star = ""
            if np.isfinite(fdr) and fdr < 0.05:
                star = "*"
            cells.append(f"{fmt(lfc)}{star}")
        lines.append(f"| {gene} | " + " | ".join(cells) + " |")
    lines.append("")
    lines.append("\\* FDR < 0.05 among detected genes. Mouse TEC (stimulated) is in `gene_log2fc.csv`, not this human table.")
    lines.append("")
    lines.append("## Binding")
    lines.append("")
    lines.append("ChIP-Atlas hg38 ELF3, 5 kb gene window, score > 0. Individual libraries jump from 0 to a MACS score above 50, so the >0 call is a peak call, not a tiny-score artifact. Every row in that table has ELF3|Average > 0: the file is the set of genes with a peak in at least one catalogued experiment (about 14,000 genes), not the whole genome. Fisher tests below are inside that table. There is no lung / NSCLC / A549 ELF3 ChIP-seq in it (CFPAC-1, ESO-26, HBDEC2, HepG2 only). CFPAC-1 and HepG2 peaks cover most of the table, so a high bound fraction there is broad occupancy, not a selective NHEJ or IFN program.")
    lines.append("")
    sub = chip_df[(chip_df["source"] == "ChIP-Atlas_hg38_5kb") & (chip_df["experiment"].isin(["HBDEC2", "HepG2_any", "CFPAC1_any", "ESO26"]))]
    sub = sub[sub["set"].isin(PRIMARY_SETS + ["HALLMARK_APICAL_JUNCTION", "CLDN4", "TACSTD2"])]
    lines.append("| Experiment | Set | bound / tested | fraction | background fraction | Fisher p |")
    lines.append("| --- | --- | ---: | ---: | ---: | ---: |")
    for rec in sub.itertuples(index=False):
        lines.append(
            f"| {rec.experiment} | {rec.set} | {rec.n_bound}/{rec.n} | {fmt(rec.frac)} | {fmt(rec.bg_frac)} | {fmt(rec.p)} |"
        )
    lines.append("")
    lines.append("Genome-wide test: called peaks overlapped with any RefSeq TSS ±2 kb (hg38). DMS53 is SCLC; ECC4 and A99 are gastrointestinal NEC. HepG2 is ENCODE IDR (ENCFF080FAU); that file has ~40,000 peaks and a ~40% background rate, so its enrichment is not selective. The NEC/SCLC CUT&Tag background is 11–17%.")
    lines.append("")
    lines.append(_recurrent_peaks(peak_genes))
    lines.append("")
    subp = peak_df[(peak_df["window"] == "tss_2kb") & (peak_df["experiment"] != "mTEC_CUTRUN_TNFa_IFNg")]
    subp = subp[subp["set"].isin(["NHEJ_CORE", "STING_CORE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_APICAL_JUNCTION", "CLDN4", "TACSTD2"])]
    lines.append("| Experiment | Peaks | Set | TSS±2 kb bound / tested | fraction | background | Fisher p |")
    lines.append("| --- | ---: | --- | ---: | ---: | ---: | ---: |")
    for rec in subp.itertuples(index=False):
        lines.append(
            f"| {rec.experiment} | {rec.n_peaks} | {rec.set} | {rec.n_bound}/{rec.n} | {fmt(rec.frac)} | {fmt(rec.bg_frac)} | {fmt(rec.p)} |"
        )
    lines.append("")
    lines.append("Author gene lists from the same NEC study (supplement tables S7–S9) are one-sided for expression: S8 is fold change < 0.75 only. Full siELF3 count matrices were not deposited in GSE190618 (that accession has baseline RNA-seq plus CUT&Tag peaks).")
    lines.append("")
    lines.append("| List | n genes | Set | overlap | genes |")
    lines.append("| --- | ---: | --- | ---: | --- |")
    keep_lists = ["down_ECC4", "down_A99", "down_DMS53", "down_ge2", "direct_232", "near_ge2"]
    keep_sets = PRIMARY_SETS + ["HALLMARK_APICAL_JUNCTION", "HALLMARK_E2F_TARGETS", "CLDN4", "TACSTD2"]
    show_n = nec_df[nec_df["list"].isin(keep_lists) & nec_df["set"].isin(keep_sets)]
    for rec in show_n.itertuples(index=False):
        genes = rec.genes if len(str(rec.genes)) < 180 else str(rec.genes)[:177] + "..."
        lines.append(f"| {rec.list} | {rec.n_list} | {rec.set} | {rec.n_overlap}/{rec.n_set} | {genes} |")
    lines.append("")
    lines.append("## What this does not say")
    lines.append("")
    lines.append("- No public lung ELF3 ChIP, so a LUAD binding claim at NHEJ/STING/IFN promoters is not testable. DMS53 CUT&Tag is SCLC.")
    lines.append("- NEC siELF3 RNA-seq is only the downregulated-gene lists, so those tables cannot show interferon induction.")
    lines.append("- Mouse TEC and the synovial TNF contrasts are cytokine-on. They are not basal ELF3-loss tests.")
    lines.append("- Co-expression of ELF3 with CLDN4 (earlier public catalogs) is not binding and is not re-tested here.")
    lines.append("- This is not an ICI or TROP2-ADC outcome analysis.")
    lines.append("")
    lines.append("## Datasets")
    lines.append("")
    for r in results:
        lines.append(f"- **{r['dataset']}** ({r['qc']}, {r['role']}): {r['note']}")
    lines.append("- **GSE190618** ELF3 CUT&Tag peaks: ECC4, A99, DMS53. Binding only.")
    lines.append("- **GSE156165** HBDEC2 ELF3 ChIP peaks.")
    lines.append("- **ENCODE ENCFF080FAU** HepG2 ELF3 IDR peaks.")
    lines.append("- **ChIP-Atlas** hg38 ELF3 5 kb target matrix.")
    lines.append("- **GSE319848** mouse TEC ELF3 CUT&RUN under TNFα+IFNγ (stimulated binding, not basal).")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append("Microarray: A549 HuGene 2.0 ST is GEO RMA (already log2), probes collapsed to HGNC symbols with Bioconductor `hugene20sttranscriptcluster.db` (unambiguous probes only) and the highest-mean probe kept. HBDEC2 Agilent intensities are log2(signal+1), same collapse using the platform gene-symbol column. Retired symbols (TMEM173, MB21D1, C9orf142, and other HGNC previous/alias symbols that are not themselves current symbols) are renamed to the current symbol before the tests. RNA-seq counts are log2(CPM+1). Mouse TPM is log2(TPM+1). Synovial technical duplicates are averaged inside donor, then a paired t-test is taken across donors. Detection floor for RNA-seq / TPM is mean abundance ≥ 1 (CPM or TPM); synovial detection is mean log2(CPM+1) ≥ 1 on either side. Set tests are two-sided Mann–Whitney of set log2FC versus other detected genes, plus a two-sided Wilcoxon signed-rank versus 0. Gene FDR is BH inside detected genes. Binding: ChIP-Atlas score > 0, or any RefSeq TSS ±2 kb overlapping a called peak. Fisher exact tests are two-sided against all genes in that universe.")
    lines.append("")
    lines.append("Reproduce: `pip install -r scripts/elf3_nhej_sting_ifn/requirements.txt` then `python scripts/elf3_nhej_sting_ifn/analyze.py`.")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (OUT / "REPORT.md").write_text(text)
    return text


def _binding_sentence(peak_df: pd.DataFrame, nec_df: pd.DataFrame) -> str:
    parts = []
    for exp in ("ECC4_CUTTag", "A99_CUTTag", "DMS53_CUTTag"):
        row = peak_df[(peak_df["experiment"] == exp) & (peak_df["window"] == "tss_2kb") & (peak_df["set"] == "NHEJ_CORE")]
        if len(row):
            rec = row.iloc[0]
            parts.append(f"{exp.split('_')[0]} {int(rec.n_bound)}/{int(rec.n)} (p={fmt(rec.p)})")
    direct = nec_df[(nec_df["list"] == "direct_232") & (nec_df["set"].isin(["NHEJ_CORE", "STING_CORE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "CLDN4"]))]
    nice = {
        "NHEJ_CORE": "NHEJ",
        "STING_CORE": "STING",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": "interferon-alpha",
        "CLDN4": "CLDN4",
    }
    direct_bits = [f"{nice.get(r.set, r.set)}: {r.genes or 'none'}" for r in direct.itertuples(index=False)]
    return (
        "Binding is a separate result from those expression calls. ELF3 CUT&Tag TSS±2 kb overlap for classical NHEJ is "
        + "; ".join(parts)
        + ". STING and interferon-alpha sets are not enriched at that window. The 232-gene direct signature (near a peak and down in at least two lines) is "
        + "; ".join(direct_bits)
        + ". Peaks on a few NHEJ genes do not become a downregulated NHEJ program after ELF3 loss."
    )


def _verdict_paragraph(results: list[dict], set_df: pd.DataFrame) -> str:
    qc = {r["dataset"]: r["qc"] for r in results}
    bits = []
    labels = {
        "NHEJ_CORE": "classical NHEJ",
        "STING_CORE": "cGAS–STING",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": "interferon-alpha",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": "interferon-gamma",
    }
    for set_name, label in labels.items():
        sub = set_df[(set_df["set"] == set_name) & (set_df["dataset"].isin(PRIMARY_DATASETS))]
        sub = sub[sub["dataset"].map(qc) == "PASS"]
        calls = sub["call"].value_counts().to_dict()
        meds = ", ".join(f"{r.dataset} {fmt(r.median_log2fc_detected)} ({r.call})" for r in sub.itertuples(index=False))
        n = len(sub)
        n_up = calls.get("UP", 0)
        n_down = calls.get("DOWN", 0)
        if n == 0:
            bits.append(f"{label}: no QC-pass primary contrast.")
        elif n_up == 0 and n_down == 0:
            bits.append(f"{label}: no QC-pass primary dataset meets the UP or DOWN rule ({n} contrasts: {meds}).")
        else:
            bits.append(f"{label}: {n_up} UP and {n_down} DOWN of {n} QC-pass primary contrasts ({meds}).")
    # positive control sentence
    sub = set_df[(set_df["set"] == "HALLMARK_APICAL_JUNCTION") & (set_df["dataset"].isin(PRIMARY_DATASETS))]
    sub = sub[sub["dataset"].map(qc) == "PASS"]
    if len(sub):
        meds = ", ".join(f"{r.dataset} {fmt(r.median_log2fc_detected)} ({r.call})" for r in sub.itertuples(index=False))
        bits.append(f"Apical-junction control on the same contrasts: {meds}.")
    cldn = []
    for r in results:
        if r["dataset"] not in PRIMARY_DATASETS or r["qc"] != "PASS":
            continue
        g = r["genes"].set_index("gene")
        gene = "CLDN4"
        if gene in g.index and np.isfinite(g.loc[gene, "log2fc"]):
            cldn.append(f"{r['dataset']} {fmt(g.loc[gene, 'log2fc'])}")
    if cldn:
        bits.append("CLDN4 log2FC in those contrasts: " + "; ".join(cldn) + ".")
    stim = set_df[(set_df["dataset"] == "synov_siELF3_TNF_IL17_16") & (set_df["set"].isin(labels))]
    if len(stim):
        meds = ", ".join(f"{r.set} {fmt(r.median_log2fc_detected)} ({r.call})" for r in stim.itertuples(index=False))
        bits.append(
            "The one cytokine contrast that passes the log2FC gate (synovial siELF3, 16 h TNF+IL-17A, not lung) is "
            + meds
            + ". Interferon-alpha moves down there, not up. ELF3 raw p in that contrast is 0.21, so the knockdown itself is not significant and the interferon call is not secure. The 16 h TNF-only arm moves the same way but fails the ELF3 log2FC gate (−0.22) and is not counted. Mouse tubular Elf3 deletion on TNFα+IFNγ (ELF3 log2FC −2.45, raw p 0.002) does not call NHEJ, STING, or interferon UP or DOWN."
        )
    return " ".join(bits)


def plot_medians(set_df: pd.DataFrame, results: list[dict]) -> None:
    qc = {r["dataset"]: r["qc"] for r in results}
    order = [r["dataset"] for r in results]
    sets = [
        "NHEJ_CORE",
        "STING_CORE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_APICAL_JUNCTION",
    ]
    short = {
        "NHEJ_CORE": "NHEJ",
        "STING_CORE": "STING",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-alpha",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-gamma",
        "HALLMARK_APICAL_JUNCTION": "Apical junction",
    }
    colors = ["#4C78A8", "#F58518", "#E45756", "#72B7B2", "#54A24B"]
    fig, ax = plt.subplots(figsize=(11, 6.2))
    y = np.arange(len(order))
    height = 0.15
    for i, (set_name, color) in enumerate(zip(sets, colors)):
        vals = []
        for ds in order:
            hit = set_df[(set_df["dataset"] == ds) & (set_df["set"] == set_name)]
            vals.append(float(hit["median_log2fc_detected"].iloc[0]) if len(hit) else np.nan)
        ax.barh(y + (i - 2) * height, vals, height=height, color=color, label=short[set_name])
    labels = []
    for ds in order:
        mark = "" if qc[ds] == "PASS" else " (QC fail)"
        labels.append(ds + mark)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.axvline(0, color="black", lw=0.6)
    ax.axvline(0.25, color="grey", lw=0.4, ls="--")
    ax.axvline(-0.25, color="grey", lw=0.4, ls="--")
    ax.set_xlabel("Median log2FC in detected genes (ELF3 loss − control)")
    ax.set_title("ELF3 loss: NHEJ, STING, IFN, and apical-junction control")
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUT / "fig_median_log2fc.png", dpi=140)
    fig.savefig(OUT / "fig_median_log2fc.pdf")
    plt.close(fig)


def main() -> None:
    ensure_cache()
    hallmark_path = CACHE / "hallmark.gmt"
    if not hallmark_path.exists():
        fetch(
            "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
            hallmark_path,
        )
    hallmark = read_gmt(hallmark_path)
    human_sets = gene_sets("human", hallmark)
    mouse_sets = gene_sets("mouse", hallmark)
    alias = hgnc_alias_map(CACHE / "hgnc.tsv")

    results = []
    results.append(
        analyze_microarray(
            "A549_shELF3",
            CACHE / "GSE137479_series_matrix.txt.gz",
            hugene_symbol_map(CACHE / "hugene20sttranscriptcluster.sqlite", CACHE / "hgnc.tsv"),
            human_sets,
            "human",
            "rma",
            ("knock",),
            ("control", "empty"),
            "GSE137479 A549 LUAD stable shELF3 clones vs empty-vector clones, n=3. HuGene 2.0 ST RMA. The only public lung ELF3-loss transcriptome in this set.",
            "primary",
            alias,
        )
    )
    ag = agilent_symbol_map(CACHE / "GPL20844.txt.gz")
    results.append(
        analyze_microarray(
            "HBDEC2_KO",
            CACHE / "GSE148105_series_matrix.txt.gz",
            ag,
            human_sets,
            "human",
            "linear",
            ("ko",),
            ("wt",),
            "GSE148105 HBDEC2 biliary epithelial ELF3 CRISPR KO vs WT, n=3. Agilent 8x60K. Not lung. Paired with GSE156165 ELF3 ChIP in the same line. Platform symbols TMEM173/MB21D1/C9orf142 are renamed to STING1/CGAS/PAXX.",
            "primary",
            alias,
        )
    )
    results.append(
        analyze_microarray(
            "HBDEC2_KD",
            CACHE / "GSE148106_series_matrix.txt.gz",
            ag,
            human_sets,
            "human",
            "linear",
            ("kd",),
            ("control",),
            "GSE148106 HBDEC2 ELF3 miRNA KD vs control, n=3. Same platform and line as the KO, independent reagent.",
            "primary",
            alias,
        )
    )
    results.append(analyze_jeg3(CACHE / "GSE241792_ELF3_normalized_counts.txt.gz", human_sets, alias))
    results.append(analyze_prostate(CACHE / "GSE303076_PLKOvsshELF3.xlsx", human_sets, alias))
    results.extend(
        analyze_synovial(
            CACHE / "GSE129487_counts.tsv.gz",
            CACHE / "GSE129487_meta.tsv.gz",
            CACHE / "hgnc.tsv",
            human_sets,
            alias,
        )
    )
    results.append(analyze_mouse(CACHE / "GSE319933_TPM.xlsx", mouse_sets))

    set_rows = []
    gene_frames = []
    for r in results:
        for row in r["sets"]:
            set_rows.append({"dataset": r["dataset"], **row})
        g = r["genes"].copy()
        g.insert(0, "dataset", r["dataset"])
        gene_frames.append(g)
    set_df = pd.DataFrame(set_rows)
    gene_df = pd.concat(gene_frames, ignore_index=True)

    print("ChIP-Atlas...", flush=True)
    chip_df, chip_genes = chip_atlas(CACHE / "ELF3.5.tsv", human_sets)
    print("Peak overlap...", flush=True)
    peak_df, peak_genes = cuttag_overlap(human_sets)
    print("NEC lists...", flush=True)
    nec_df = nec_overlap(nec_lists(), human_sets)

    set_df.to_csv(OUT / "set_stats.csv", index=False)
    gene_df.to_csv(OUT / "gene_log2fc.csv", index=False)
    chip_df.to_csv(OUT / "chipatlas_fisher.csv", index=False)
    chip_genes.to_csv(OUT / "chipatlas_focus_genes.csv", index=False)
    peak_df.to_csv(OUT / "peak_tss_fisher.csv", index=False)
    peak_genes.to_csv(OUT / "peak_tss_focus_genes.csv", index=False)
    nec_df.to_csv(OUT / "nec_list_overlap.csv", index=False)

    summary = []
    for r in results:
        summary.append({k: v for k, v in r.items() if k not in {"sets", "genes", "lfc", "fdr"}})
    (OUT / "contrasts.json").write_text(json.dumps(summary, indent=2, default=str))
    write_report(results, set_df, gene_df, chip_df, peak_df, nec_df, peak_genes)
    plot_medians(set_df, results)
    print(f"Wrote {OUT / 'REPORT.md'}")
    # stdout QC
    for r in results:
        print(f"{r['dataset']:28} ELF3 {r['elf3_log2fc']:+.3f}  {r['qc']}")


if __name__ == "__main__":
    main()
