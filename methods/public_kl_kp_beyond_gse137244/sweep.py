#!/usr/bin/env python3
"""Histology, probe, normalization, and TJ-module sweep for public KL vs KP.

Primary test stays two-sided exact Mann-Whitney, delta = mean(KL) - mean(KP).
Welch t is reported beside it and is not substituted for the rank test.
Specifications that compare squamous KL tumors with adenocarcinoma KP tumors
are stored and marked fair_genotype=False.

Does not read or merge private 8KL matrices.
"""

from __future__ import annotations

import gzip
import json
import math
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze as base  # noqa: E402

CACHE = base.CACHE
TABLES = base.TABLES
UA = base.UA

TJ_GENES = ["Cldn1", "Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cldn18", "Ocln", "Tjp1", "Tjp2", "Tjp3", "Marveld2", "Marveld3", "Cgn", "F11r", "Cdh1"]
MODULES = {
    "TJ7": ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"],
    "TJ_CORE": ["Ocln", "Tjp1", "Tjp2", "Tjp3", "F11r", "Cgn", "Marveld2", "Cldn3", "Cldn4", "Cldn7"],
    "TJ_CORE_NO_CLDN4": ["Ocln", "Tjp1", "Tjp2", "Tjp3", "F11r", "Cgn", "Marveld2", "Cldn3", "Cldn7"],
    "CLDN_STRAND": ["Cldn3", "Cldn4", "Cldn7", "Ocln"],
    "PLAQUE": ["Tjp1", "Tjp2", "Tjp3", "Cgn", "Marveld2"],
}
NHEJ = ["Xrcc4", "Xrcc5", "Xrcc6", "Lig4", "Prkdc", "Nhej1", "Dclre1c", "Paxx"]
STING = ["Tmem173", "Sting1", "Mb21d1", "Cgas", "Tbk1", "Irf3"]
IFN = ["Stat1", "Stat2", "Irf1", "Irf7", "Irf9", "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Oasl2", "Rsad2", "Ifih1", "Ddx58", "Ifnb1"]
EPI = ["Epcam", "Krt8"]
HIST = ["Krt5", "Krt14", "Trp63", "Nkx2-1", "Sftpc", "Sftpb"]
SINGLE = ["Tacstd2", "Cldn4"] + TJ_GENES
ALL_SYMBOLS = sorted(set(SINGLE + EPI + HIST + NHEJ + STING + IFN + sum(MODULES.values(), [])))


def mwu(a, b) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 1 or len(b) < 1:
        return float("nan")
    if np.allclose(a, a[0]) and np.allclose(b, b[0]) and float(a[0]) == float(b[0]):
        return 1.0
    return float(mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue)


def welch(a, b) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    if np.std(a) == 0 and np.std(b) == 0:
        return 1.0 if float(a[0]) == float(b[0]) else float("nan")
    return float(ttest_ind(a, b, equal_var=False, alternative="two-sided").pvalue)


def min_p(n: int, m: int) -> float:
    if n < 1 or m < 1:
        return float("nan")
    return 2.0 / math.comb(n + m, n)


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out.tolist()
    idx = np.where(ok)[0]
    pv = p[idx]
    order = np.argsort(pv)
    m = len(pv)
    q = np.empty(m)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        rank = i + 1
        prev = min(prev, float(pv[order[i]]) * m / rank)
        q[i] = prev
    adj = np.empty(m)
    adj[order] = np.clip(q, 0, 1)
    out[idx] = adj
    return out.tolist()


def ensembl_map(symbols: list[str]) -> dict[str, str]:
    dest = CACHE / "ensembl_symbol_ids.json"
    cached: dict[str, str] = {}
    if dest.exists():
        cached = json.loads(dest.read_text())
    missing = [s for s in symbols if s not in cached]
    if missing:
        url = "https://rest.ensembl.org/lookup/symbol/mus_musculus"
        req = urllib.request.Request(
            url,
            data=json.dumps({"symbols": missing}).encode(),
            headers={"Content-Type": "application/json", "Accept": "application/json", "User-Agent": UA},
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            payload = json.loads(resp.read().decode())
        for sym in missing:
            rec = payload.get(sym) or {}
            cached[sym] = rec.get("id") or ""
        dest.write_text(json.dumps(cached, indent=2) + "\n")
    return {k: v for k, v in cached.items() if v}


def homer_symbols(path: Path, symbols: set[str]) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        ann_i = header.index("Annotation/Divergence")
        samples = [h.replace(" FPKM", "").replace("Aligned.out.sam", "") for h in header[ann_i + 1 :]]
        buckets: dict[str, list[list[float]]] = {}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= ann_i:
                continue
            symbol = parts[ann_i].split("|", 1)[0]
            if symbol not in symbols:
                continue
            vals = [float(x) if x not in ("", "NA") else math.nan for x in parts[ann_i + 1 :]]
            buckets.setdefault(symbol, []).append(vals)
    mat = {sym: np.nanmean(np.vstack(rows), axis=0) for sym, rows in buckets.items()}
    out = pd.DataFrame(mat, index=samples).T
    return np.log2(out.astype(float) + 1.0)


def read_symbol_table(df: pd.DataFrame, symbols: list[str], log: bool) -> pd.DataFrame:
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean(numeric_only=True)
    idx = {str(i): i for i in df.index}
    rows = {}
    for sym in symbols:
        key = idx.get(sym)
        if key is None:
            continue
        rows[sym] = pd.to_numeric(df.loc[key], errors="coerce")
    out = pd.DataFrame(rows).T
    out = out.groupby(level=0).mean()
    if log:
        out = np.log2(out.astype(float) + 1.0)
    return out


def gpl_probes(symbols: list[str]) -> dict[str, list[str]]:
    path = base.fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8321/annot/GPL8321.annot.gz",
        CACHE / "GPL8321.annot.gz",
    )
    want = {s.lower(): s for s in symbols}
    found: dict[str, list[str]] = {s: [] for s in symbols}
    with gzip.open(path, "rt", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith(("!", "#", "^")):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            first = rec.get("Gene symbol", "").split("///")[0].strip().lower()
            if first in want:
                found[want[first]].append(rec["ID"])
    return found


def gse6135_probe_matrix(probe_map: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = base.fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6135/matrix/GSE6135-GPL8321_series_matrix.txt.gz",
        CACHE / "GSE6135-GPL8321_series_matrix.txt.gz",
    )
    probes = base._series_matrix_table(path)
    rows = {}
    for gene, ids in probe_map.items():
        for pid in ids:
            if pid in probes.index:
                rows[f"{gene}|{pid}"] = probes.loc[pid]
    return pd.DataFrame(rows).T, probes


def arms_164758(columns: list[str]) -> dict[str, list[str]]:
    kl, kp = [], []
    for col in columns:
        if col.startswith("KL") and not col.startswith("KPL") and "H457" not in col:
            kl.append(col)
        elif col.startswith("KP") and not col.startswith("KPL") and not col.startswith("KPa"):
            kp.append(col)
    return {"KL": kl, "KP": kp}


def z_module(log_mat: pd.DataFrame, genes: list[str], cols: list[str]) -> pd.Series | None:
    present = [g for g in genes if g in log_mat.index]
    if len(present) < 3:
        return None
    block = log_mat.loc[present, cols].astype(float)
    z = block.sub(block.mean(axis=1), axis=0).div(block.std(axis=1).replace(0, np.nan), axis=0)
    if z.notna().sum().sum() == 0:
        return None
    return z.mean(axis=0)


def residual(y: pd.Series, cov: pd.DataFrame) -> pd.Series:
    common = y.index.intersection(cov.index)
    yy = y.loc[common].astype(float)
    xx = cov.loc[common].astype(float)
    keep = [c for c in xx.columns if xx[c].notna().all() and xx[c].std() > 0]
    if not keep or yy.isna().any():
        return pd.Series(dtype=float)
    x = np.column_stack([np.ones(len(common)), xx[keep].to_numpy()])
    beta, *_ = np.linalg.lstsq(x, yy.to_numpy(), rcond=None)
    fitted = x @ beta
    return pd.Series(yy.to_numpy() - fitted, index=common)


class Specs:
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def add(
        self,
        accession: str,
        cohort: str,
        family: str,
        name: str,
        kl: np.ndarray,
        kp: np.ndarray,
        *,
        fair: bool,
        scale: str,
        note: str,
    ) -> None:
        kl = np.asarray(kl, float)
        kp = np.asarray(kp, float)
        kl = kl[np.isfinite(kl)]
        kp = kp[np.isfinite(kp)]
        if len(kl) < 2 or len(kp) < 2:
            return
        delta = float(np.mean(kl) - np.mean(kp))
        self.rows.append(
            {
                "accession": accession,
                "cohort": cohort,
                "family": family,
                "name": name,
                "fair_genotype": bool(fair),
                "scale": scale,
                "n_kl": int(len(kl)),
                "n_kp": int(len(kp)),
                "mean_kl": float(np.mean(kl)),
                "mean_kp": float(np.mean(kp)),
                "delta_mean": delta,
                "kl_gt_kp": bool(delta > 0),
                "mwu_p": mwu(kl, kp),
                "welch_p": welch(kl, kp),
                "mwu_floor": min_p(len(kl), len(kp)),
                "enough_n": bool(min_p(len(kl), len(kp)) < 0.05),
                "note": note,
            }
        )


def gene_and_modules(specs: Specs, acc: str, cohort: str, log_mat: pd.DataFrame, kl: list[str], kp: list[str], scale: str, fair: bool, note: str) -> None:
    if len(kl) < 2 or len(kp) < 2:
        return
    cols = kl + kp
    for gene in list(dict.fromkeys(["Tacstd2", "Cldn4"] + TJ_GENES)):
        if gene not in log_mat.index:
            continue
        specs.add(acc, cohort, "gene", gene, log_mat.loc[gene, kl], log_mat.loc[gene, kp], fair=fair, scale=scale, note=note)
    for mod, genes in MODULES.items():
        score = z_module(log_mat, genes, cols)
        if score is None:
            continue
        present = [g for g in genes if g in log_mat.index]
        specs.add(
            acc,
            cohort,
            "module",
            mod,
            score[kl],
            score[kp],
            fair=fair,
            scale="mean gene z within contrast samples",
            note=note + f" Genes used: {','.join(present)}.",
        )


def epithelial_norm(specs: Specs, acc: str, cohort: str, log_mat: pd.DataFrame, kl: list[str], kp: list[str], fair: bool, note: str) -> None:
    if len(kl) < 2 or len(kp) < 2:
        return
    if "Cldn4" not in log_mat.index or "Epcam" not in log_mat.index:
        return
    cols = kl + kp
    ratio = log_mat.loc["Cldn4", cols] - log_mat.loc["Epcam", cols]
    specs.add(
        acc,
        cohort,
        "norm",
        "Cldn4_minus_Epcam",
        ratio[kl],
        ratio[kp],
        fair=fair,
        scale="Cldn4 minus Epcam on the analysis scale",
        note=note + " Epithelial ratio, not raw Cldn4.",
    )
    if "Krt8" in log_mat.index:
        cov = log_mat.loc[["Epcam", "Krt8"], cols].T
        cov.columns = ["Epcam", "Krt8"]
        resid = residual(log_mat.loc["Cldn4", cols], cov)
        if len(resid):
            specs.add(
                acc,
                cohort,
                "norm",
                "Cldn4_residual_Epcam_Krt8",
                resid[kl],
                resid[kp],
                fair=fair,
                scale="OLS residual of Cldn4 on Epcam and Krt8",
                note=note + " Residual removes linear epithelial-marker content. Fit uses KL and KP samples together.",
            )
    core = z_module(log_mat, MODULES["TJ_CORE"], cols)
    if core is not None and "Epcam" in log_mat.index:
        ep = log_mat.loc["Epcam", cols].astype(float)
        epz = (ep - ep.mean()) / (ep.std() if ep.std() else np.nan)
        adj = core - epz
        specs.add(
            acc,
            cohort,
            "norm",
            "TJ_CORE_minus_Epcam_z",
            adj[kl],
            adj[kp],
            fair=fair,
            scale="TJ_CORE z minus Epcam z",
            note=note + " Module after subtracting the Epcam z-score.",
        )


def pathway_modules(specs: Specs, acc: str, cohort: str, log_mat: pd.DataFrame, kl: list[str], kp: list[str], note: str) -> None:
    if len(kl) < 2 or len(kp) < 2:
        return
    cols = kl + kp
    # Resolve STING aliases to one score without double-counting.
    sting_genes = []
    for a, b in (("Tmem173", "Sting1"), ("Mb21d1", "Cgas")):
        if a in log_mat.index:
            sting_genes.append(a)
        elif b in log_mat.index:
            sting_genes.append(b)
    for extra in ("Tbk1", "Irf3"):
        if extra in log_mat.index:
            sting_genes.append(extra)
    named = {"NHEJ": [g for g in NHEJ if g in log_mat.index], "STING": sting_genes, "IFN": [g for g in IFN if g in log_mat.index]}
    for name, genes in named.items():
        score = z_module(log_mat, genes, cols)
        if score is None:
            continue
        specs.add(
            acc,
            cohort,
            "pathway",
            name,
            score[kl],
            score[kp],
            fair=True,
            scale="mean gene z within contrast samples",
            note=note + f" Genes used: {','.join(genes)}.",
        )
    for gene in list(dict.fromkeys(NHEJ + STING + IFN)):
        if gene not in log_mat.index:
            continue
        specs.add(acc, cohort, "pathway_gene", gene, log_mat.loc[gene, kl], log_mat.loc[gene, kp], fair=True, scale="analysis scale", note=note)


def histology_proxy(log_mat: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    sq = [g for g in ("Krt5", "Krt14", "Trp63") if g in log_mat.index]
    ad = [g for g in ("Nkx2-1", "Sftpc", "Sftpb") if g in log_mat.index]
    block = log_mat.loc[sq + ad, cols].astype(float)
    z = block.sub(block.mean(axis=1), axis=0).div(block.std(axis=1).replace(0, np.nan), axis=0)
    out = pd.DataFrame(index=cols)
    out["squamous_z"] = z.loc[sq].mean(axis=0) if sq else np.nan
    out["adeno_z"] = z.loc[ad].mean(axis=0) if ad else np.nan
    out["call"] = np.where(out["adeno_z"] > out["squamous_z"], "adeno_like", "squamous_like")
    out["sq_genes"] = ",".join(sq)
    out["ad_genes"] = ",".join(ad)
    return out


def mouse_mean(probe_or_gene: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    labels = []
    for mouse, grp in meta.groupby("mouse"):
        pieces.append(probe_or_gene[grp["gsm"].tolist()].mean(axis=1))
        labels.append(mouse)
    out = pd.concat(pieces, axis=1)
    out.columns = labels
    return out


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    specs = Specs()

    # ----- GSE164758 primary tumors -----
    homer = base.fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164758/suppl/GSE164758_primary_tumors_fpkm.txt.gz",
        CACHE / "GSE164758_primary_tumors_fpkm.txt.gz",
    )
    m164 = homer_symbols(homer, set(ALL_SYMBOLS))
    a164 = arms_164758(list(m164.columns))
    note164 = (
        "Eichner primary tumors. No histology field in GEO (genotype, tissue, treatment, guide RNA only). "
        "Bulk FPKM. Unannotated KL-H457 and KPa columns excluded."
    )
    gene_and_modules(specs, "GSE164758", "all_tumors", m164, a164["KL"], a164["KP"], "log2(FPKM+1)", True, note164)
    epithelial_norm(specs, "GSE164758", "all_tumors", m164, a164["KL"], a164["KP"], True, note164)
    pathway_modules(specs, "GSE164758", "all_tumors", m164, a164["KL"], a164["KP"], note164)
    calls = histology_proxy(m164, a164["KL"] + a164["KP"])
    calls["arm"] = ["KL" if i in a164["KL"] else "KP" for i in calls.index]
    calls.to_csv(TABLES / "gse164758_histology_proxy.tsv", sep="\t")
    for call, fair in (("adeno_like", True), ("squamous_like", False)):
        kl = [s for s in a164["KL"] if calls.loc[s, "call"] == call]
        kp = [s for s in a164["KP"] if calls.loc[s, "call"] == call] if fair else a164["KP"]
        # Matched call vs matched call is the fair filter. Squamous-like KL vs all KP is not.
        if not fair:
            kp = a164["KP"]
        tag = "matched_marker_histology" if fair else "squamous_like_KL_vs_all_KP"
        hnote = note164 + f" Marker call {call}. sq={calls['sq_genes'].iloc[0]} ad={calls['ad_genes'].iloc[0]}."
        if fair:
            hnote += " Both arms restricted to the same marker call."
        else:
            hnote += " UNFAIR if KP tumors are adeno-like: squamous program vs genotype."
        gene_and_modules(specs, "GSE164758", tag, m164, kl, kp, "log2(FPKM+1)", fair, hnote)
        epithelial_norm(specs, "GSE164758", tag, m164, kl, kp, fair, hnote)

    # ----- RNA-seq symbol matrices -----
    def load_137244() -> pd.DataFrame:
        path = base.fetch(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.fpkm.csv.gz",
            CACHE / "GSE137244_counts.fpkm.csv.gz",
        )
        df = pd.read_csv(path, index_col=0)
        return read_symbol_table(df, ALL_SYMBOLS, log=True)

    def load_137396() -> pd.DataFrame:
        path = base.fetch(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137396/suppl/GSE137396_Raw_genetable_GEMMnodule.txt.gz",
            CACHE / "GSE137396_Raw_genetable_GEMMnodule.txt.gz",
        )
        df = pd.read_csv(path, sep="\t", index_col=0)
        return read_symbol_table(df, ALL_SYMBOLS, log=True)

    def load_244452() -> pd.DataFrame:
        path = base.fetch(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244452/suppl/GSE244452_KPvsKL_deg_all.txt.gz",
            CACHE / "GSE244452_KPvsKL_deg_all.txt.gz",
        )
        df = pd.read_csv(path, sep="\t").set_index("gene_name")
        cols = ["KP1", "KP2", "KP3", "KL1", "KL2", "KL3"]
        return read_symbol_table(df[cols], ALL_SYMBOLS, log=True)

    def load_274352() -> pd.DataFrame:
        path = base.fetch(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274352/suppl/GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz",
            CACHE / "GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz",
        )
        df = pd.read_csv(path, sep="\t")
        df = df.set_index("external_gene_name")
        empty = [c for c in df.columns if "empty" in str(c)]
        return read_symbol_table(df[empty], ALL_SYMBOLS, log=True)

    def load_274351() -> pd.DataFrame:
        path = base.fetch(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274351/suppl/GSE274351_expredata_TPM_gene.txt.gz",
            CACHE / "GSE274351_expredata_TPM_gene.txt.gz",
        )
        df = pd.read_csv(path, sep="\t", index_col=0)
        ids = ensembl_map(ALL_SYMBOLS)
        rows = {}
        for sym, eid in ids.items():
            if eid in df.index:
                rows[sym] = pd.to_numeric(df.loc[eid], errors="coerce")
        out = pd.DataFrame(rows).T
        return np.log2(out.astype(float) + 1.0)

    m244 = load_244452()
    note244 = (
        "Syngeneic bulk tumors. Deposited table is 7449 genes with per-sample counts, not a full transcriptome. "
        "n=3 vs 3; exact Mann-Whitney floor is 0.10. No histology field."
    )
    gene_and_modules(specs, "GSE244452", "all_tumors", m244, ["KL1", "KL2", "KL3"], ["KP1", "KP2", "KP3"], "log2(count+1)", True, note244)
    epithelial_norm(specs, "GSE244452", "all_tumors", m244, ["KL1", "KL2", "KL3"], ["KP1", "KP2", "KP3"], True, note244)

    m396 = load_137396()
    kl396 = [c for c in m396.columns if str(c).startswith("KL")]
    kp396 = [c for c in m396.columns if "Trp53" in str(c) or str(c).startswith("KP")]
    note396 = "Deng GEMM nodules. KL labeled NA, KP labeled vehicle. log2 of the deposited abundance."
    gene_and_modules(specs, "GSE137396", "all_nodules", m396, kl396, kp396, "log2(abundance+1)", True, note396)
    epithelial_norm(specs, "GSE137396", "all_nodules", m396, kl396, kp396, True, note396)
    pathway_modules(specs, "GSE137396", "all_nodules", m396, kl396, kp396, note396)

    m244_lines = load_137244()
    kl244 = [c for c in m244_lines.columns if str(c).startswith("KL")]
    kp244 = [c for c in m244_lines.columns if str(c).startswith("B6AL10")]
    note_lines = "Locked GSE137244 cell lines. Included so TJ modules use the same definition as the other series."
    gene_and_modules(specs, "GSE137244", "cell_lines", m244_lines, kl244, kp244, "log2(FPKM+1)", True, note_lines)
    pathway_modules(specs, "GSE137244", "cell_lines", m244_lines, kl244, kp244, note_lines)

    m352 = load_274352()
    kl352 = [c for c in m352.columns if str(c).startswith("KL") and "empty" in str(c)]
    kp352 = [c for c in m352.columns if str(c).startswith("KP") and "empty" in str(c)]
    note352 = "Empty-vector cell lines. n=3 vs 3; Mann-Whitney floor 0.10. Not used for pathway claims."
    gene_and_modules(specs, "GSE274352", "empty_vector", m352, kl352, kp352, "log2(normalized_count+1)", True, note352)
    epithelial_norm(specs, "GSE274352", "empty_vector", m352, kl352, kp352, True, note352)

    m351 = load_274351()
    kl351 = [c for c in m351.columns if str(c).startswith("KL")]
    kp351 = [c for c in m351.columns if str(c).startswith("KP")]
    note351 = "LCM early adenomas. Column labels KL1-5 and KP1-5 from the TPM file. GEO SOFT lists KL n=4 and K n=5."
    gene_and_modules(specs, "GSE274351", "lcm_adenoma", m351, kl351, kp351, "log2(TPM+1)", True, note351)
    epithelial_norm(specs, "GSE274351", "lcm_adenoma", m351, kl351, kp351, True, note351)
    pathway_modules(specs, "GSE274351", "lcm_adenoma", m351, kl351, kp351, note351)

    # ----- GSE6135 probes and histology -----
    probe_map = gpl_probes(ALL_SYMBOLS)
    (TABLES / "gse6135_probes.tsv").write_text(
        "gene\tn_probes\tprobes\n" + "\n".join(f"{g}\t{len(ids)}\t{','.join(ids)}" for g, ids in sorted(probe_map.items()) if ids) + "\n"
    )
    pmat, _ = gse6135_probe_matrix(probe_map)
    # gene = mean of its probes
    gene_rows = {}
    for gene, ids in probe_map.items():
        keys = [f"{gene}|{pid}" for pid in ids if f"{gene}|{pid}" in pmat.index]
        if keys:
            gene_rows[gene] = pmat.loc[keys].mean(axis=0)
    gmat = pd.DataFrame(gene_rows).T
    meta = pd.DataFrame(base.GSE6135_TUMORS, columns=["gsm", "mouse", "arm", "histology", "primary"])
    meta = meta[meta["primary"] & meta["gsm"].isin(gmat.columns)]
    cohorts = {
        "mouse_all_histology": (meta, True, "L/L and L/- primaries, all histologies. Mouse mean. L/+ and metastasis excluded."),
        "mouse_adeno_only": (meta[meta["histology"].eq("Ad")], True, "Pathologist adenocarcinoma only. Matched to KP, which are all Ad."),
        "mouse_adeno_plus_mixed": (
            meta[meta["histology"].isin(["Ad", "Ad-sq"])],
            True,
            "Adenocarcinoma plus adenosquamous. Pure squamous dropped. KP remain Ad.",
        ),
        "mouse_sq_mixed_vs_ad_KP": (
            meta[(meta["arm"].eq("KP")) | (meta["histology"].isin(["Sq", "Ad-sq"]))],
            False,
            "UNFAIR histology: squamous or mixed KL versus adenocarcinoma KP.",
        ),
    }
    for cohort, (sub, fair, note) in cohorts.items():
        if sub["arm"].nunique() < 2:
            continue
        mm = mouse_mean(gmat, sub)
        kl = [c for c in mm.columns if str(c).startswith("KL_")]
        kp = [c for c in mm.columns if str(c).startswith("KP_")]
        full_note = "GSE6135 GPL8321 deposited array values, not re-logged. " + note
        gene_and_modules(specs, "GSE6135", cohort, mm, kl, kp, "deposited_array_value", fair, full_note)
        epithelial_norm(specs, "GSE6135", cohort, mm, kl, kp, fair, full_note)
        if fair and cohort in ("mouse_all_histology", "mouse_adeno_only"):
            pathway_modules(specs, "GSE6135", cohort, mm, kl, kp, full_note)
        # probes: mouse means, same cohort
        pm = mouse_mean(pmat, sub)
        for key in pm.index:
            gene, pid = key.split("|", 1)
            if gene not in set(["Tacstd2", "Cldn4"] + TJ_GENES):
                continue
            nprob = len(probe_map.get(gene, []))
            specs.add(
                "GSE6135",
                cohort,
                "probe",
                f"{gene}:{pid}",
                pm.loc[key, kl],
                pm.loc[key, kp],
                fair=fair,
                scale="deposited_array_value",
                note=full_note + f" Single probe. {gene} has {nprob} exact-symbol probe(s) on GPL8321.",
            )
    # tumor-level (pseudoreplication) for the fair all-histology set
    tumor_meta = meta.copy()
    tumor_meta["mouse"] = tumor_meta["gsm"]
    tm = mouse_mean(gmat, tumor_meta)
    kl_t = [c for c in tm.columns if c in set(meta.loc[meta["arm"].eq("KL"), "gsm"])]
    kp_t = [c for c in tm.columns if c in set(meta.loc[meta["arm"].eq("KP"), "gsm"])]
    gene_and_modules(
        specs,
        "GSE6135",
        "tumor_level_all_histology",
        tm,
        kl_t,
        kp_t,
        "deposited_array_value",
        False,
        "Tumor-level, not mouse-level. Two primaries from the same mouse are not independent. Shown so the unit choice is visible.",
    )

    table = pd.DataFrame(specs.rows)
    fdrs = []
    for _, grp in table.groupby(["accession", "cohort", "family"], sort=False):
        fdrs.append(pd.Series(bh(grp["mwu_p"].tolist()), index=grp.index))
    table["mwu_fdr_bh"] = pd.concat(fdrs).sort_index()
    table.to_csv(TABLES / "sweep_specs.tsv", sep="\t", index=False)

    focus_names = {"Cldn4", "Tacstd2", "Cldn4_minus_Epcam", "Cldn4_residual_Epcam_Krt8", "TJ_CORE_minus_Epcam_z"} | set(MODULES)
    focus = table[table["name"].isin(focus_names) | table["name"].astype(str).str.startswith("Cldn4:") | table["family"].isin(["module", "norm", "pathway"])].copy()
    focus.to_csv(TABLES / "sweep_focus.tsv", sep="\t", index=False)

    hits = table[(table["delta_mean"] > 0) & ((table["mwu_p"] < 0.05) | (table["welch_p"] < 0.05))].copy()
    hits = hits.sort_values(["fair_genotype", "mwu_p", "welch_p"], ascending=[False, True, True])
    hits.to_csv(TABLES / "sweep_hits_kl_gt_kp.tsv", sep="\t", index=False)

    print("specs", len(table), "hits", len(hits))
    show = hits[hits["name"].isin(focus_names) | hits["family"].isin(["module", "norm", "pathway", "pathway_gene"]) | hits["name"].astype(str).str.startswith("Cldn")]
    cols = ["accession", "cohort", "name", "fair_genotype", "n_kl", "n_kp", "delta_mean", "mwu_p", "welch_p", "mwu_fdr_bh", "enough_n"]
    print(show[cols].head(80).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
