#!/usr/bin/env python3
"""Public mouse KL vs KP beyond the locked GSE137244 cell-line result.

KL = Kras-mutant, Lkb1/Stk11-null, Trp53 intact.
KP = Kras-mutant, Trp53-null, Lkb1 intact.
KPL = Kras-mutant, Trp53-null and Lkb1-null (secondary only; not called KL).

Expression scale matches the locked GSE137244 numbers:
log2(abundance + 1), delta = mean(KL) - mean(KP), two-sided exact Mann-Whitney.
TJ7 is the mean of Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, Ocln on that scale.
The locked write-up quoted TJ +3.03; this seven-gene mean on the same FPKM is
reported separately and is not relabeled as that +3.03.

Does not read or merge private 8KL matrices.
"""

from __future__ import annotations

import gzip
import json
import math
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables" / "kl_kp"
FIG = ROOT / "figures" / "kl_kp"
CACHE = ROOT / "cache" / "kl_kp"
UA = "public-kl-kp-beyond-gse137244/1.0"

TJ_GENES = ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
FOCUS = ["Tacstd2", "Cldn4"] + TJ_GENES

# GSE6135 mouse tumors. L/+ heterozygotes and the metastasis are excluded.
# Mouse id collapses multiple primaries from the same animal.
GSE6135_TUMORS = [
    # gsm, mouse, arm, histology, primary_cohort
    ("GSM176566", "K_268", "K", "Ad", False),
    ("GSM176567", "K_287", "K", "Ad", False),
    ("GSM176568", "K_405", "K", "Ad", False),
    ("GSM176569", "K_484", "K", "Ad", False),
    ("GSM176571", "K_498", "K", "Ad", False),
    ("GSM176573", "KP_186", "KP", "Ad", True),
    ("GSM176574", "KP_196", "KP", "Ad", True),
    ("GSM176576", "KP_197", "KP", "Ad", True),
    ("GSM176578", "KP_498", "KP", "Ad", True),
    ("GSM176580", "KP_500", "KP", "Ad", True),
    ("GSM176585", "KL_861", "KL", "Ad", True),
    ("GSM176587", "KL_113", "KL", "Sq", True),
    ("GSM176591", "KL_452", "KL", "Ad-sq", True),
    ("GSM176593", "KL_452", "KL", "Ad-sq", True),
    ("GSM176595", "KL_459", "KL", "Ad", True),
    ("GSM176596", "KL_540", "KL", "Sq", True),
    ("GSM176597", "KL_540", "KL", "Sq", True),
    ("GSM176598", "KL_547", "KL", "Ad", True),
    ("GSM176599", "KL_592", "KL", "Ad", True),
    ("GSM176600", "KL_592", "KL", "Ad", True),
]

# Affymetrix Mouse Genome 430 2.0 (GPL8321) probes, exact gene symbol.
GSE6135_PROBES = {
    "Tacstd2": ["1423323_at"],
    "Cldn4": ["1418283_at"],
    "Cldn3": ["1426332_a_at", "1434651_a_at", "1451701_x_at", "1460569_x_at"],
    "Cldn6": ["1417845_at"],
    "Cldn7": ["1448393_at"],
    "Cdh1": ["1448261_at"],
    "F11r": ["1424595_at", "1436374_x_at"],
    "Ocln": ["1448873_at"],
}


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as resp:
        dest.write_bytes(resp.read())
    return dest


def mwu_exact(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 1 or len(b) < 1:
        return float("nan")
    if np.allclose(a, a[0]) and np.allclose(b, b[0]) and a[0] == b[0]:
        return 1.0
    return float(mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue)


def summarize(kl: np.ndarray, kp: np.ndarray) -> dict:
    kl = np.asarray(kl, dtype=float)
    kp = np.asarray(kp, dtype=float)
    kl = kl[np.isfinite(kl)]
    kp = kp[np.isfinite(kp)]
    return {
        "n_kl": int(len(kl)),
        "n_kp": int(len(kp)),
        "mean_kl": float(np.mean(kl)) if len(kl) else float("nan"),
        "mean_kp": float(np.mean(kp)) if len(kp) else float("nan"),
        "median_kl": float(np.median(kl)) if len(kl) else float("nan"),
        "median_kp": float(np.median(kp)) if len(kp) else float("nan"),
        "delta_mean": float(np.mean(kl) - np.mean(kp)) if len(kl) and len(kp) else float("nan"),
        "mwu_p": mwu_exact(kl, kp),
        "kl_values": ";".join(f"{v:.6g}" for v in kl),
        "kp_values": ";".join(f"{v:.6g}" for v in kp),
    }


def _unique_genes(log_mat: pd.DataFrame) -> pd.DataFrame:
    if log_mat.index.duplicated().any():
        return log_mat.groupby(level=0).mean()
    return log_mat


def rows_for_matrix(
    accession: str,
    log_mat: pd.DataFrame,
    kl_cols: list[str],
    kp_cols: list[str],
    *,
    contrast: str,
    unit: str,
    setting: str,
    scale: str,
    role: str,
    note: str,
    genes: list[str] | None = None,
) -> list[dict]:
    """log_mat is genes x samples, already on the comparison scale."""
    log_mat = _unique_genes(log_mat)
    genes = list(dict.fromkeys(genes or (FOCUS + ["TJ7"])))
    present_tj = [g for g in TJ_GENES if g in log_mat.index]
    out = []
    for gene in genes:
        if gene == "TJ7":
            if len(present_tj) < 5:
                continue
            a = log_mat.loc[present_tj, kl_cols].mean(axis=0).to_numpy()
            b = log_mat.loc[present_tj, kp_cols].mean(axis=0).to_numpy()
            gene_note = note + f" TJ7 mean of {','.join(present_tj)}."
        else:
            if gene not in log_mat.index:
                out.append(
                    {
                        "accession": accession,
                        "contrast": contrast,
                        "gene": gene,
                        "unit": unit,
                        "setting": setting,
                        "scale": scale,
                        "role": role,
                        "n_kl": len(kl_cols),
                        "n_kp": len(kp_cols),
                        "mean_kl": np.nan,
                        "mean_kp": np.nan,
                        "median_kl": np.nan,
                        "median_kp": np.nan,
                        "delta_mean": np.nan,
                        "mwu_p": np.nan,
                        "kl_values": "",
                        "kp_values": "",
                        "note": note + " Gene absent from the deposited matrix.",
                    }
                )
                continue
            a = log_mat.loc[gene, kl_cols].to_numpy(dtype=float)
            b = log_mat.loc[gene, kp_cols].to_numpy(dtype=float)
            gene_note = note
        st = summarize(a, b)
        out.append(
            {
                "accession": accession,
                "contrast": contrast,
                "gene": gene,
                "unit": unit,
                "setting": setting,
                "scale": scale,
                "role": role,
                **st,
                "note": gene_note,
            }
        )
    return out


def sample_frame(accession: str, arm_of: dict[str, str], log_mat: pd.DataFrame) -> pd.DataFrame:
    log_mat = _unique_genes(log_mat)
    present = [g for g in FOCUS if g in log_mat.index]
    tj = [g for g in TJ_GENES if g in log_mat.index]
    recs = []
    for col, arm in arm_of.items():
        rec = {"accession": accession, "sample": col, "arm": arm}
        for g in present:
            rec[g] = float(log_mat.loc[g, col])
        if len(tj) >= 5:
            rec["TJ7"] = float(log_mat.loc[tj, col].mean())
        recs.append(rec)
    return pd.DataFrame(recs)


def log2p1(df: pd.DataFrame) -> pd.DataFrame:
    return np.log2(df.astype(float) + 1.0)


def gse137244() -> tuple[list[dict], pd.DataFrame]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.fpkm.csv.gz",
        CACHE / "GSE137244_counts.fpkm.csv.gz",
    )
    df = pd.read_csv(path, index_col=0)
    kl = [c for c in df.columns if str(c).startswith("KL")]
    kp = [c for c in df.columns if str(c).startswith("B6AL10")]
    log_mat = log2p1(df.loc[df.index.intersection(FOCUS)])
    note = (
        "Locked reference, reproduced. KL = KL155mix, KL47-1, KLC, KLD, KLE. "
        "KP = B6AL10-1..5. Normal lung excluded. Not a new dataset."
    )
    rows = rows_for_matrix(
        "GSE137244",
        log_mat,
        kl,
        kp,
        contrast="KL_minus_KP",
        unit="cell_line",
        setting="cultured_cells",
        scale="log2(FPKM+1)",
        role="locked_reference",
        note=note,
    )
    arms = {c: "KL" for c in kl} | {c: "KP" for c in kp}
    return rows, sample_frame("GSE137244", arms, log_mat)


def gse137396() -> tuple[list[dict], pd.DataFrame]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137396/suppl/GSE137396_Raw_genetable_GEMMnodule.txt.gz",
        CACHE / "GSE137396_Raw_genetable_GEMMnodule.txt.gz",
    )
    df = pd.read_csv(path, sep="\t", index_col=0)
    kl = [c for c in df.columns if str(c).startswith("KL")]
    kp = [c for c in df.columns if "Kras-Trp53" in str(c) or str(c).startswith("KP")]
    log_mat = log2p1(df.loc[df.index.intersection(FOCUS)])
    note = (
        "Same Deng et al. paper as GSE137244, but GEMM lung nodules rather than cell lines. "
        "KL samples are labeled NA; KP samples are labeled vehicle. "
        "Previously Welch on a log matrix in an ICI screen (Tacstd2 delta about +1.08)."
    )
    rows = rows_for_matrix(
        "GSE137396",
        log_mat,
        kl,
        kp,
        contrast="KL_minus_KP",
        unit="nodule",
        setting="in_vivo_nodule",
        scale="log2(deposited_abundance+1)",
        role="rescore",
        note=note,
    )
    arms = {c: ("KL" if c in kl else "KP") for c in kl + kp}
    return rows, sample_frame("GSE137396", arms, log_mat)


def gse274351() -> tuple[list[dict], pd.DataFrame]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274351/suppl/GSE274351_expredata_TPM_gene.txt.gz",
        CACHE / "GSE274351_expredata_TPM_gene.txt.gz",
    )
    df = pd.read_csv(path, sep="\t", index_col=0)
    # Supplementary column prefixes. GEO SOFT lists K n=5 and KL n=4, while this
    # file has K1-4 and KL1-5. Stk11 is near the floor in all five KL columns.
    ens = {
        "Tacstd2": "ENSMUSG00000051397",
        "Cldn4": "ENSMUSG00000047501",
        "Cldn3": "ENSMUSG00000070473",
        "Cldn6": "ENSMUSG00000023906",
        "Cldn7": "ENSMUSG00000018569",
        "Cdh1": "ENSMUSG00000000303",
        "F11r": "ENSMUSG00000038235",
        "Ocln": "ENSMUSG00000021638",
        "Stk11": "ENSMUSG00000003068",
    }
    sub = df.loc[[ens[g] for g in ens]].copy()
    sub.index = list(ens.keys())
    kl = [c for c in sub.columns if str(c).startswith("KL")]
    kp = [c for c in sub.columns if str(c).startswith("KP")]
    log_mat = log2p1(sub.loc[list(dict.fromkeys(FOCUS))])
    stk = sub.loc["Stk11", kl].astype(float)
    note = (
        "LCM early adenomas (PNAS 2024). Groups follow supplementary column labels "
        f"(KL n={len(kl)}, KP n={len(kp)}). GEO SOFT instead lists K n=5 and KL n=4. "
        f"Stk11 TPM in the five KL-labeled columns: {', '.join(f'{v:.3g}' for v in stk)}. "
        "Previously scored with Welch; direction was KL < KP."
    )
    rows = rows_for_matrix(
        "GSE274351",
        log_mat,
        kl,
        kp,
        contrast="KL_minus_KP",
        unit="adenoma",
        setting="lcm_early_adenoma",
        scale="log2(TPM+1)",
        role="rescore",
        note=note,
    )
    arms = {c: "KL" for c in kl} | {c: "KP" for c in kp}
    return rows, sample_frame("GSE274351", arms, log_mat)


def gse274352() -> tuple[list[dict], pd.DataFrame]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274352/suppl/GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz",
        CACHE / "GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz",
    )
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.set_index("external_gene_name")
    df = df.groupby(level=0).mean()
    empty = {
        "KL": [c for c in df.columns if str(c).startswith("KL") and "empty" in str(c)],
        "KP": [c for c in df.columns if str(c).startswith("KP") and "empty" in str(c)],
    }
    log_mat = log2p1(df.loc[df.index.intersection(FOCUS)])
    note = (
        "Untreated empty-vector cell lines from the same study as GSE274351. "
        "IFN-beta and STING arms are not in this contrast. n=3 vs 3 cannot reach "
        "two-sided exact Mann-Whitney p < 0.05 (floor 0.10)."
    )
    rows = rows_for_matrix(
        "GSE274352",
        log_mat,
        empty["KL"],
        empty["KP"],
        contrast="KL_minus_KP",
        unit="cell_line",
        setting="cultured_cells_empty_vector",
        scale="log2(normalized_count+1)",
        role="new",
        note=note,
    )
    arms = {c: "KL" for c in empty["KL"]} | {c: "KP" for c in empty["KP"]}
    return rows, sample_frame("GSE274352", arms, log_mat)


def gse244452() -> tuple[list[dict], pd.DataFrame]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244452/suppl/GSE244452_KPvsKL_deg_all.txt.gz",
        CACHE / "GSE244452_KPvsKL_deg_all.txt.gz",
    )
    df = pd.read_csv(path, sep="\t")
    wide = df.set_index("gene_name")[["KP1", "KP2", "KP3", "KL1", "KL2", "KL3"]]
    wide = wide.groupby(level=0).mean()
    log_mat = log2p1(wide)
    present_tj = [g for g in TJ_GENES if g in log_mat.index]
    missing = [g for g in TJ_GENES if g not in log_mat.index]
    note = (
        "Syngeneic KL vs KP bulk tumors, n=3 vs 3. Deposited table has per-sample "
        f"values for {wide.shape[0]} genes, not a full transcriptome. "
        f"TJ genes missing from the deposit: {', '.join(missing) or 'none'}. "
        "Two-sided exact Mann-Whitney floor at n=3 vs 3 is 0.10."
    )
    rows = rows_for_matrix(
        "GSE244452",
        log_mat,
        ["KL1", "KL2", "KL3"],
        ["KP1", "KP2", "KP3"],
        contrast="KL_minus_KP",
        unit="tumor",
        setting="syngeneic_bulk_tumor",
        scale="log2(deposited_count+1)",
        role="new",
        note=note,
        genes=FOCUS + (["TJ7"] if len(present_tj) >= 5 else []),
    )
    arms = {"KL1": "KL", "KL2": "KL", "KL3": "KL", "KP1": "KP", "KP2": "KP", "KP3": "KP"}
    return rows, sample_frame("GSE244452", arms, log_mat)


def _homer_symbol_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        # First field is the Homer command; FPKM names start after Annotation/Divergence.
        try:
            ann_i = header.index("Annotation/Divergence")
        except ValueError:
            ann_i = 8
        sample_names = []
        for name in header[ann_i + 1 :]:
            sample_names.append(name.replace(" FPKM", "").replace("Aligned.out.sam", ""))
        buckets: dict[str, list[list[float]]] = {}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= ann_i:
                continue
            symbol = parts[ann_i].split("|", 1)[0]
            if symbol not in FOCUS:
                continue
            vals = [float(x) if x not in ("", "NA") else math.nan for x in parts[ann_i + 1 :]]
            buckets.setdefault(symbol, []).append(vals)
    cols = sample_names
    mat = {}
    for symbol, rows in buckets.items():
        arr = np.vstack(rows)
        mat[symbol] = np.nanmean(arr, axis=0)
    out = pd.DataFrame(mat, index=cols).T
    out.columns = cols
    return out


def gse164758() -> tuple[list[dict], pd.DataFrame]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164758/suppl/GSE164758_primary_tumors_fpkm.txt.gz",
        CACHE / "GSE164758_primary_tumors_fpkm.txt.gz",
    )
    mat = _homer_symbol_matrix(path)
    # Only libraries that match the 41 GEO samples (treatment: none).
    # Extra KL-H457* and KPa1a* columns in the FPKM file are not in the SOFT sample list.
    kl, kp, kpl = [], [], []
    for col in mat.columns:
        if col.startswith("KL") and not col.startswith("KPL") and "H457" not in col:
            kl.append(col)
        elif col.startswith("KP") and not col.startswith("KPL") and not col.startswith("KPa"):
            kp.append(col)
        elif col.startswith("KPL"):
            kpl.append(col)
    log_mat = log2p1(mat)
    note = (
        "Eichner et al. primary NSCLC tumors, treatment none. "
        f"Annotated GEO samples used: KL n={len(kl)}, KP n={len(kp)}. "
        "Unannotated extra FPKM columns (KL-H457*, KPa1a*) were excluded. "
        "Bulk tumor RNA, so stroma is in the measurement."
    )
    rows = rows_for_matrix(
        "GSE164758",
        log_mat,
        kl,
        kp,
        contrast="KL_minus_KP",
        unit="tumor",
        setting="primary_bulk_tumor",
        scale="log2(FPKM+1)",
        role="new",
        note=note,
    )
    kpl_note = (
        "Secondary contrast on the same matrix. KPL is Kras/p53/Lkb1, not KL. "
        f"KPL n={len(kpl)} vs KP n={len(kp)}."
    )
    rows += rows_for_matrix(
        "GSE164758",
        log_mat,
        kpl,
        kp,
        contrast="KPL_minus_KP",
        unit="tumor",
        setting="primary_bulk_tumor",
        scale="log2(FPKM+1)",
        role="secondary_kpl",
        note=kpl_note,
    )
    arms = {c: "KL" for c in kl} | {c: "KP" for c in kp} | {c: "KPL" for c in kpl}
    return rows, sample_frame("GSE164758", arms, log_mat)


def _series_matrix_table(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as handle:
        lines = []
        started = False
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                started = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if started:
                lines.append(line)
    from io import StringIO

    df = pd.read_csv(StringIO("".join(lines)), sep="\t", index_col=0)
    df.index = [str(i).strip('"') for i in df.index]
    df.columns = [str(c).strip('"') for c in df.columns]
    return df.apply(pd.to_numeric, errors="coerce")


def gse6135() -> tuple[list[dict], pd.DataFrame]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6135/matrix/GSE6135-GPL8321_series_matrix.txt.gz",
        CACHE / "GSE6135-GPL8321_series_matrix.txt.gz",
    )
    probes = _series_matrix_table(path)
    gene_rows = {}
    for gene, ids in GSE6135_PROBES.items():
        present = [i for i in ids if i in probes.index]
        if not present:
            continue
        gene_rows[gene] = probes.loc[present].mean(axis=0)
    log_mat = pd.DataFrame(gene_rows).T  # already log-scale array values
    meta = pd.DataFrame(
        GSE6135_TUMORS, columns=["gsm", "mouse", "arm", "histology", "primary"]
    )
    meta = meta[meta["gsm"].isin(log_mat.columns)]

    def mouse_means(sub: pd.DataFrame) -> pd.DataFrame:
        pieces = []
        for mouse, grp in sub.groupby("mouse"):
            vals = log_mat[grp["gsm"].tolist()].mean(axis=1)
            vals.name = mouse
            pieces.append(vals)
        return pd.concat(pieces, axis=1)

    rows: list[dict] = []
    samples = []
    cohorts = {
        "mouse_all_histology": meta[meta["primary"]],
        "mouse_adeno_only": meta[meta["primary"] & meta["histology"].eq("Ad")],
    }
    for cohort_name, sub in cohorts.items():
        mm = mouse_means(sub)
        kl_m = [m for m in mm.columns if m.startswith("KL_")]
        kp_m = [m for m in mm.columns if m.startswith("KP_")]
        note = (
            "Ji et al. 2007 Mouse Genome 430 2.0. Values are deposited array intensities "
            "(already on a log-like scale; not log-transformed again). "
            "Unit is the mouse: multiple primaries averaged. "
            "L/+ heterozygotes and the 592 metastasis excluded. "
        )
        if cohort_name == "mouse_all_histology":
            note += "KL includes adenocarcinoma, adenosquamous, and squamous primaries (L/L or L/-)."
            role = "new_tacstd2_same_cohort_as_prior_cldn4"
        else:
            note += "Sensitivity: adenocarcinoma primaries only."
            role = "sensitivity_adeno_only"
        rows += rows_for_matrix(
            "GSE6135",
            mm,
            kl_m,
            kp_m,
            contrast="KL_minus_KP",
            unit="mouse",
            setting=cohort_name,
            scale="deposited_array_value",
            role=role,
            note=note,
        )
        arms = {m: ("KL" if m.startswith("KL_") else "KP") for m in kl_m + kp_m}
        sf = sample_frame("GSE6135", arms, mm)
        sf["cohort"] = cohort_name
        samples.append(sf)
    return rows, pd.concat(samples, ignore_index=True)


def context_markers() -> pd.DataFrame:
    """Epithelial and immune context for the two new bulk-tumor matrices."""
    rows = []
    # GSE244452 deposited per-sample counts.
    path = CACHE / "GSE244452_KPvsKL_deg_all.txt.gz"
    df = pd.read_csv(path, sep="\t").set_index("gene_name")
    cols = ["KP1", "KP2", "KP3", "KL1", "KL2", "KL3"]
    for gene in ["Epcam", "Krt8", "Sftpc", "Ptprc", "Cd8a", "Cldn4", "Tacstd2"]:
        if gene not in df.index:
            continue
        vals = df.loc[gene, cols].astype(float)
        rows.append(
            {
                "accession": "GSE244452",
                "gene": gene,
                "scale": "deposited_count",
                "mean_kl": float(vals[["KL1", "KL2", "KL3"]].mean()),
                "mean_kp": float(vals[["KP1", "KP2", "KP3"]].mean()),
                "values": ";".join(f"{c}={vals[c]:.4g}" for c in cols),
            }
        )
    return pd.DataFrame(rows)


def inventory_rows() -> list[dict]:
    """Series inspected and not given a KL vs KP expression test."""
    skipped = [
        (
            "GSE165640",
            "no_processed_matrix",
            "Bulk lungs from KL, KL5, KL9, KP, KP5, KP9. Pure KP is n=2. "
            "Supplementary files are genotype-versus-knockout compare tables, not a KL vs KP matrix. Not scored.",
        ),
        (
            "GSE277929",
            "no_kp_arm",
            "KL GEMM tumors under vehicle, trametinib, entinostat, or the combination. No KP samples.",
        ),
        (
            "GSE193895",
            "no_kp_arm",
            "Kras, Kras/Lkb1, Kras/Keap1, Kras/Keap1/Lkb1. No Trp53-null KP arm.",
        ),
        (
            "GSE21581",
            "no_kp_arm",
            "Carretero microarray of Kras vs Kras/Lkb1 primary tumors and metastases. No KP arm.",
        ),
        (
            "GSE69552",
            "no_kp_arm",
            "KL cell-of-origin histotypes. No KP arm.",
        ),
        (
            "GSE180963",
            "no_kp_arm",
            "scRNA of Kras vs KL, two samples. Not a KL vs KP contrast, and not merged with private 8KL.",
        ),
        (
            "GSE165641",
            "no_kp_arm",
            "KL scRNA. No KP arm. Not merged with private 8KL.",
        ),
        (
            "GSE154977",
            "no_kl_arm",
            "KP scRNA time points. No KL arm. Not merged with private 8KL.",
        ),
        (
            "KP-Tracer (Yang Cell 2022)",
            "not_kl",
            "scRNA of KP vs KPL vs KPA. KPL is Lkb1 loss on a KP background, not KL vs KP. Not scored here.",
        ),
        (
            "organoids",
            "none_found",
            "GEO query for mouse lung organoid plus Lkb1/Stk11 returned patient-derived organoid series "
            "(GSE233468/GSE233665), not Kras/Lkb1 vs Kras/p53 organoids. Dost GSE150425 is Kras-only. "
            "GSE227719 is KP organoids only. No public KL vs KP organoid matrix with usable n.",
        ),
    ]
    rows = []
    for acc, status, note in skipped:
        rows.append(
            {
                "accession": acc,
                "contrast": "KL_minus_KP",
                "gene": "",
                "unit": "",
                "setting": "",
                "scale": "",
                "role": "inventory_not_scored",
                "n_kl": np.nan,
                "n_kp": np.nan,
                "mean_kl": np.nan,
                "mean_kp": np.nan,
                "median_kl": np.nan,
                "median_kp": np.nan,
                "delta_mean": np.nan,
                "mwu_p": np.nan,
                "kl_values": "",
                "kp_values": "",
                "note": f"{status}. {note}",
            }
        )
    return rows


def plot_deltas(contrasts: pd.DataFrame, path: Path) -> None:
    keep_roles = {
        "locked_reference",
        "rescore",
        "new",
        "new_tacstd2_same_cohort_as_prior_cldn4",
        "sensitivity_adeno_only",
        "secondary_kpl",
    }
    sub = contrasts[
        contrasts["role"].isin(keep_roles) & contrasts["gene"].isin(["Tacstd2", "Cldn4"])
    ].copy()
    if sub.empty:
        return
    sub["label"] = sub.apply(
        lambda r: f"{r['accession']} {r['setting']}\n{r['contrast']} n={int(r['n_kl'])} vs {int(r['n_kp'])}",
        axis=1,
    )
    # stable order
    order = list(dict.fromkeys(sub["label"]))
    fig, ax = plt.subplots(figsize=(8.2, 0.48 * len(order) + 1.4))
    ymap = {lab: i for i, lab in enumerate(order[::-1])}
    for gene, color, dx in (("Tacstd2", "#C44E52", -0.12), ("Cldn4", "#4C78A8", 0.12)):
        part = sub[sub["gene"] == gene]
        ys = [ymap[lab] + dx for lab in part["label"]]
        ax.scatter(part["delta_mean"], ys, c=color, s=36, label=gene, zorder=3)
    ax.axvline(0, color="0.45", lw=0.8)
    ax.set_yticks(list(ymap.values()))
    ax.set_yticklabels(list(ymap.keys()), fontsize=8)
    ax.set_xlabel("mean(test arm) − mean(KP) on the dataset scale")
    ax.set_title("Public mouse KL vs KP: Tacstd2 and Cldn4")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / ".gitignore").write_text("*\n!.gitignore\n")

    all_rows: list[dict] = []
    sample_frames = []
    for fn in (gse137244, gse137396, gse164758, gse244452, gse274352, gse274351, gse6135):
        print("running", fn.__name__)
        rows, samples = fn()
        all_rows.extend(rows)
        sample_frames.append(samples)
    all_rows.extend(inventory_rows())

    contrasts = pd.DataFrame(all_rows)
    contrasts.to_csv(TABLES / "contrasts.tsv", sep="\t", index=False)
    samples = pd.concat(sample_frames, ignore_index=True)
    samples.to_csv(TABLES / "sample_values.tsv", sep="\t", index=False)

    scored = contrasts[contrasts["role"] != "inventory_not_scored"].copy()
    focus = scored[scored["gene"].isin(["Tacstd2", "Cldn4", "TJ7"])]
    focus.to_csv(TABLES / "focus_contrasts.tsv", sep="\t", index=False)
    context_markers().to_csv(TABLES / "gse244452_context_markers.tsv", sep="\t", index=False)

    plot_deltas(contrasts, FIG / "kl_minus_kp_delta.png")

    summary = {
        "test": "two-sided exact Mann-Whitney; delta = mean(test arm) - mean(KP)",
        "tj7": TJ_GENES,
        "locked_gse137244_expected": {
            "Tacstd2_delta": 3.24,
            "Cldn4_delta": 5.57,
            "mwu_p": 0.00794,
        },
        "focus": json.loads(focus.to_json(orient="records")),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    show = focus[focus["gene"].isin(["Tacstd2", "Cldn4"])][
        ["accession", "contrast", "setting", "role", "n_kl", "n_kp", "delta_mean", "mwu_p"]
    ]
    print(show.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
