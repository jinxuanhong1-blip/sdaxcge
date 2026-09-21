#!/usr/bin/env python3
"""Cldn4 versus ICB response and versus NHEJ/IFN in Kras-mutant mouse lung.

Public matrices only. TISMO's locked Tacstd2 49/64 pairing is not recomputed.
LLC is KrasG12C/Nras lung carcinoma, not KL. CMT-167 is KrasG12V and has no
ICB rows in the TISMO export. MLE12 is SV40 large-T, not Kras.

Thesis signs, fixed before the sweep:
  response / acquired resistance: Cldn4 higher in the resistant group
  IFN / chemokine / STING transcripts: Cldn4 negatively associated
  c-NHEJ transcripts: Cldn4 positively associated
    (lower NHEJ -> cGAS/STING -> IFN; high Cldn4 tracks the IFN-low state)

The sweep varies module membership, score, correlation, Cldn4 split,
expression scale, and sample filter. The headline is the strongest
thesis-aligned statistic among unconfounded tests. Every tested row is
written out, with BH q inside its question.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from gene_sets import ENSEMBL, IFN_MODULES, MODULE_THESIS_SIGN, MODULES, NHEJ_MODULES

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "kras_lung_icb_cldn4_sweep"
DATA = Path(os.environ.get("KRAS_LUNG_ICB_DATA", "/tmp/kras_icb"))

MIN_ARM = 3
MIN_CORR_N = 8
MIN_MODULE_FRACTION = 0.7


def bh(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg q values. Non-finite p are left as nan and ignored."""
    m_idx = [i for i, p in enumerate(pvalues) if p is not None and math.isfinite(p)]
    q = [float("nan")] * len(pvalues)
    m = len(m_idx)
    if m == 0:
        return q
    order = sorted(m_idx, key=lambda i: pvalues[i])
    running = 1.0
    adj = {}
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * m / rank)
        adj[i] = running
    for i, val in adj.items():
        q[i] = val
    return q


def finite(x: float) -> bool:
    return x is not None and isinstance(x, (int, float)) and math.isfinite(x)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 3 or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return float("nan"), float("nan")
    rho, p = stats.spearmanr(x, y)
    rho = float(rho)
    p = float(p)
    # A perfect correlation makes the t approximation explode to p = 0.
    if (not math.isfinite(p)) or p == 0.0:
        p = permutation_spearman_p(x, y, rho)
    return rho, p


def permutation_spearman_p(x: np.ndarray, y: np.ndarray, rho: float) -> float:
    n = len(x)
    if n <= 7:
        from itertools import permutations

        count = 0
        total = 0
        y = np.asarray(y, dtype=float)
        for order in permutations(range(n)):
            total += 1
            r, _ = stats.spearmanr(x, y[list(order)])
            if abs(float(r)) + 1e-12 >= abs(rho):
                count += 1
        return count / total
    rng = np.random.default_rng(1)
    nperm = 20000
    count = 1
    y = np.asarray(y, dtype=float)
    for _ in range(nperm):
        r, _ = stats.spearmanr(x, rng.permutation(y))
        if abs(float(r)) + 1e-12 >= abs(rho):
            count += 1
    return count / (nperm + 1)


def safe_pearson(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 3 or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return float("nan"), float("nan")
    r, p = stats.pearsonr(x, y)
    return float(r), float(p)


def safe_mwu(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 1 or len(b) < 1:
        return float("nan")
    if np.all(a == a[0]) and np.all(b == b[0]) and a[0] == b[0]:
        return 1.0
    try:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided", method="auto").pvalue)
    except ValueError:
        return float("nan")


def safe_welch(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2 or (np.nanstd(a) == 0 and np.nanstd(b) == 0):
        return float("nan")
    return float(stats.ttest_ind(a, b, equal_var=False, alternative="two-sided").pvalue)


def median_diff(high: np.ndarray, low: np.ndarray) -> float:
    if len(high) == 0 or len(low) == 0:
        return float("nan")
    return float(np.median(high) - np.median(low))


# ---------------------------------------------------------------------------
# Cohort definitions. Roles: resistant groups are the thesis-high Cldn4 arm.
# ---------------------------------------------------------------------------

@dataclass
class SampleAnn:
    name: str
    group: str
    # resistant: thesis-high Cldn4 arm for the response question
    # reference: the comparison arm (responder, parental, or untreated control)
    # other: kept for module correlations, not in the response contrast
    role: str
    ordinal: float | None = None


@dataclass
class Cohort:
    cohort_id: str
    accession: str
    model: str
    kras: str
    is_kl: bool
    setting: str
    response_family: str  # acquired_resistance, tismo_label, tolerant_proxy, none
    confounded: bool
    confound_note: str
    samples: list[SampleAnn]
    genes: list[str]
    # transform name -> genes x samples matrix, columns follow samples
    matrices: dict[str, np.ndarray] = field(default_factory=dict)
    gene_index: dict[str, int] = field(default_factory=dict)

    def assert_not_kl_mislabel(self) -> None:
        if self.is_kl:
            raise RuntimeError(f"{self.cohort_id} was marked KL; this sweep does not include a KL ICB matrix")
        if self.model in {"KL", "KPL", "Lkb1"}:
            raise RuntimeError(f"{self.cohort_id} model label {self.model} is not allowed")


def ann(names_groups: list[tuple[str, str, str]], ordinals: dict[str, float] | None = None) -> list[SampleAnn]:
    ordinals = ordinals or {}
    return [
        SampleAnn(name, group, role, ordinals.get(group))
        for name, group, role in names_groups
    ]


def gse155972_samples(columns: list[str]) -> list[SampleAnn]:
    out = []
    for name in columns:
        low = name.lower()
        treated = "_tx_" in low
        setdb1 = "setdb1" in low
        if treated and setdb1:
            group, role = "setdb1_ko_icb", "reference"  # TISMO Responders
        elif treated and not setdb1:
            group, role = "control_icb", "resistant"  # TISMO Non-responders
        elif setdb1:
            group, role = "setdb1_ko_untreated", "other"
        else:
            group, role = "control_untreated", "other"
        out.append(SampleAnn(name, group, role))
    return out


def gse246922_kp(columns: list[str]) -> list[SampleAnn]:
    out = []
    for name in columns:
        if name.startswith("ResKPlate"):
            group, role, ordinal = "relapse_late", "resistant", 3.0
        elif name.startswith("ResResKP"):
            group, role, ordinal = "relapse_2", "resistant", 2.0
        elif name.startswith("ResKP"):
            group, role, ordinal = "relapse_1", "resistant", 1.0
        elif name.startswith("KPy"):
            group, role, ordinal = "chronic_ifng", "other", None
        elif name.startswith("KP_"):
            group, role, ordinal = "parental", "reference", 0.0
        else:
            raise RuntimeError(f"unmapped KP column {name}")
        out.append(SampleAnn(name, group, role, ordinal))
    return out


def gse246922_llc(columns: list[str]) -> list[SampleAnn]:
    out = []
    for name in columns:
        if name.startswith("ResResLLC"):
            group, role, ordinal = "relapsed", "resistant", 1.0
        elif name.startswith("LLC1y"):
            group, role, ordinal = "chronic_ifng", "other", None
        elif name.startswith("LLC1_"):
            group, role, ordinal = "parental", "reference", 0.0
        else:
            raise RuntimeError(f"unmapped LLC1 column {name}")
        out.append(SampleAnn(name, group, role, ordinal))
    return out


def gse114601_samples(columns: list[str]) -> list[SampleAnn]:
    groups = {
        "s1795": ("anti_pd1", "other"),
        "s2521": ("anti_pd1", "other"),
        "s2596": ("vehicle", "other"),
        "s2617": ("vehicle", "other"),
        "s2503": ("jq1", "other"),
        "s2625": ("jq1", "other"),
        "s2467": ("anti_pd1_jq1", "other"),
        "s2669": ("anti_pd1_jq1", "other"),
    }
    out = []
    for name in columns:
        if name not in groups:
            raise RuntimeError(f"unmapped GSE114601 column {name}")
        group, role = groups[name]
        out.append(SampleAnn(name, group, role))
    return out


GSE157880_GROUP = {
    "0-1_S13": "igg_0gy",
    "0-2_S14": "igg_0gy",
    "0-3_S15": "igg_0gy",
    "0-4_S16": "pd1_0gy",
    "0-5_S17": "pd1_0gy",
    "4_2_S26": "igg_4gy",
    "4_3_S27": "igg_4gy",
    "4_4_S28": "pd1_4gy",
    "4_5_S29": "pd1_4gy",
    "4_6_S30": "pd1_4gy",
    "8_1_S31": "igg_8gy",
    "8_2_S32": "igg_8gy",
    "8_3_S33": "igg_8gy",
    "8_4_S34": "pd1_8gy",
    "8_5_S35": "pd1_8gy",
    "8_6_S36": "pd1_8gy",
}


def gse274960_samples(columns: list[str]) -> list[SampleAnn]:
    out = []
    for name in columns:
        if name.startswith("LL2_IgG"):
            group = "igg"
        elif name.startswith("LL2_E_PD"):
            group = "anti_pd1_entrectinib"
        elif name.startswith("LL2_E_"):
            group = "entrectinib"
        elif name.startswith("LL2_PD_"):
            group = "anti_pd1"
        elif name.startswith("LL2_CNTRL"):
            group = "shctrl"
        elif name.startswith("LL2_shNTRK"):
            group = "shntrk1"
        else:
            raise RuntimeError(f"unmapped GSE274960 column {name}")
        out.append(SampleAnn(name, group, "other"))
    return out


def gse262305_samples(columns: list[str]) -> list[SampleAnn]:
    out = []
    for name in columns:
        if name.startswith("Blank"):
            group = "isotype"
        elif name.startswith("PDL1"):
            group = "anti_pdl1"
        elif name.startswith("BTZ"):
            group = "btz_anti_pdl1"
        else:
            raise RuntimeError(f"unmapped GSE262305 column {name}")
        out.append(SampleAnn(name, group, "other"))
    return out


def gse239485_samples(columns: list[str]) -> list[SampleAnn]:
    out = []
    for name in columns:
        if name.startswith("C_"):
            group = "vehicle"
        elif name.startswith("D_"):
            group = "polyic_anti_pd1"
        elif name.startswith("T_"):
            group = "polyic_anti_pd1_c5ar1"
        else:
            raise RuntimeError(f"unmapped GSE239485 column {name}")
        out.append(SampleAnn(name, group, "other"))
    return out


def gse297630_samples(columns: list[str]) -> list[SampleAnn]:
    out = []
    for name in columns:
        if name.startswith("C-"):
            group, role = "untreated_control", "reference"
        elif name.startswith("P-"):
            group, role = "anti_pd1_tolerant", "resistant"
        else:
            raise RuntimeError(f"unmapped GSE297630 column {name}")
        out.append(SampleAnn(name, group, role))
    return out


def gse260596_samples(columns: list[str]) -> list[SampleAnn]:
    out = []
    for name in columns:
        if name.startswith("aKLH_PD1"):
            group = "anti_pd1_control"
        elif name.startswith("aLair1_PD1"):
            group = "anti_pd1_lair1"
        else:
            raise RuntimeError(f"unmapped GSE260596 column {name}")
        out.append(SampleAnn(name, group, "other"))
    return out


# ---------------------------------------------------------------------------
# Matrix loading
# ---------------------------------------------------------------------------

URLS = {
    "GSE155972_Invivo_LLC_counts_GRCm38.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155972/suppl/GSE155972_Invivo_LLC_counts_GRCm38.txt.gz",
    "GSE114601_counts.normalized.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114601/suppl/GSE114601_counts.normalized.csv.gz",
    "GSE157880_Bulk048.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157880/suppl/GSE157880_Bulk048.txt.gz",
    "GSE246922_KP_RNA_counts_vst.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246922/suppl/GSE246922_KP_RNA_counts_vst.csv.gz",
    "GSE246922_LLC1_RNA_counts_vst.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246922/suppl/GSE246922_LLC1_RNA_counts_vst.csv.gz",
    "GSE239485_Processed_data.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239485/suppl/GSE239485_Processed_data.xlsx",
    "GSE297630_processed_data.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297630/suppl/GSE297630_processed_data.xlsx",
    "GSE260596_AllSamples_Genes_ReadCounts.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260596/suppl/GSE260596_AllSamples_Genes_ReadCounts.txt.gz",
    "GSE274960_gene_count_matrix.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274960/suppl/GSE274960_gene_count_matrix.csv.gz",
    "GSE262305_gene_expression.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262305/suppl/GSE262305_gene_expression.txt.gz",
}

ID_TO_SYMBOL = {}
for _sym, _ens in ENSEMBL.items():
    ID_TO_SYMBOL[_ens] = _sym
    ID_TO_SYMBOL[_sym] = _sym
# Clariom exports and older builds use these symbols.
for _alias, _sym in {
    "H2AFX": "H2ax",
    "TMEM173": "Sting1",
    "MB21D1": "Cgas",
    "1110057K04RIK": "Paxx",
}.items():
    ID_TO_SYMBOL[_alias] = _sym
UPPER_TO_SYMBOL = {key.upper(): val for key, val in ID_TO_SYMBOL.items()}


def ensure_file(name: str) -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    dest = DATA / name
    if not dest.exists() or dest.stat().st_size < 1000:
        urllib.request.urlretrieve(URLS[name], dest)
    return dest


def canonical_gene(token: str) -> str | None:
    token = token.strip().strip('"')
    if not token:
        return None
    if token in ID_TO_SYMBOL:
        return ID_TO_SYMBOL[token]
    bare = token.split(".")[0]
    if bare in ID_TO_SYMBOL:
        return ID_TO_SYMBOL[bare]
    return UPPER_TO_SYMBOL.get(token.upper()) or UPPER_TO_SYMBOL.get(bare.upper())


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def sniff_delimiter(path: Path) -> str:
    with open_text(path) as handle:
        line = handle.readline()
    if line.count("\t") >= line.count(","):
        return "\t"
    return ","


def accumulate(path: Path, gene_col: int, sample_cols: list[int], delimiter: str) -> tuple[list[str], dict[str, np.ndarray]]:
    """Sum duplicate gene rows. Returns sample names and symbol -> values."""
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        samples = [header[i] for i in sample_cols]
        acc: dict[str, np.ndarray] = {}
        for row in reader:
            if gene_col >= len(row):
                continue
            symbol = canonical_gene(row[gene_col])
            if symbol is None:
                # Still keep every gene for within-sample ranks, under its raw id
                # only if it looks like a real gene id. Non-focus genes are stored
                # when the caller asks for a full rank matrix separately.
                continue
            vals = np.empty(len(sample_cols), dtype=np.float64)
            for j, col in enumerate(sample_cols):
                try:
                    vals[j] = float(row[col]) if row[col] not in {"", "NA", "NaN"} else np.nan
                except ValueError:
                    vals[j] = np.nan
            if symbol in acc:
                acc[symbol] = np.nansum(np.vstack([acc[symbol], vals]), axis=0)
            else:
                acc[symbol] = vals
    return samples, acc


def full_matrix_from_delimited(
    path: Path,
    gene_col: int,
    sample_cols: list[int],
    delimiter: str,
) -> tuple[list[str], np.ndarray]:
    """All rows, first occurrence of a canonical symbol wins; other rows kept for ranks."""
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        samples = [header[i] for i in sample_cols]
        rows = []
        names = []
        seen = set()
        for row in reader:
            if max(sample_cols) >= len(row):
                continue
            vals = []
            ok = True
            for col in sample_cols:
                try:
                    vals.append(float(row[col]) if row[col] not in {"", "NA", "NaN"} else np.nan)
                except ValueError:
                    ok = False
                    break
            if not ok:
                continue
            symbol = canonical_gene(row[gene_col]) if gene_col < len(row) else None
            label = symbol if symbol else f"row{len(names)}"
            if symbol and symbol in seen:
                continue
            if symbol:
                seen.add(symbol)
            names.append(label)
            rows.append(vals)
    mat = np.asarray(rows, dtype=np.float64)
    return samples, names, mat


def read_excel_matrix(path: Path, sheet: str, header_row: int, gene_col: int, sample_cols: list[int]):
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    header = None
    rows = []
    names = []
    seen = set()
    for i, row in enumerate(ws.iter_rows(values_only=True)):
        if i < header_row:
            continue
        if i == header_row:
            header = list(row)
            continue
        if row[gene_col] is None:
            continue
        vals = []
        ok = True
        for col in sample_cols:
            cell = row[col] if col < len(row) else None
            if cell is None or cell == "":
                vals.append(np.nan)
            else:
                try:
                    vals.append(float(cell))
                except (TypeError, ValueError):
                    ok = False
                    break
        if not ok:
            continue
        symbol = canonical_gene(str(row[gene_col]))
        label = symbol if symbol else f"row{len(names)}"
        if symbol and symbol in seen:
            continue
        if symbol:
            seen.add(symbol)
        names.append(label)
        rows.append(vals)
    samples = [str(header[i]) for i in sample_cols]
    return samples, names, np.asarray(rows, dtype=np.float64)


def transforms_for(kind: str, mat: np.ndarray) -> dict[str, np.ndarray]:
    if kind == "count":
        lib = np.nansum(mat, axis=0)
        lib[lib <= 0] = np.nan
        cpm = mat / lib * 1e6
        return {
            "log2_count": np.log2(mat + 1.0),
            "log2_cpm": np.log2(cpm + 1.0),
        }
    if kind == "positive":
        return {"log2_value": np.log2(np.clip(mat, a_min=0, a_max=None) + 1.0)}
    if kind == "identity":
        return {"author_scale": mat.copy()}
    raise RuntimeError(kind)


def matrix_to_cohort(
    cohort_id: str,
    accession: str,
    model: str,
    kras: str,
    setting: str,
    response_family: str,
    confounded: bool,
    confound_note: str,
    sample_names: list[str],
    genes: list[str],
    raw: np.ndarray,
    kind: str,
    annotator,
) -> Cohort:
    samples = annotator(sample_names)
    if len(samples) != len(sample_names):
        raise RuntimeError(cohort_id)
    gene_index = {g: i for i, g in enumerate(genes)}
    if "Cldn4" not in gene_index:
        raise RuntimeError(f"{cohort_id} is missing Cldn4")
    cohort = Cohort(
        cohort_id=cohort_id,
        accession=accession,
        model=model,
        kras=kras,
        is_kl=False,
        setting=setting,
        response_family=response_family,
        confounded=confounded,
        confound_note=confound_note,
        samples=samples,
        genes=genes,
        matrices=transforms_for(kind, raw),
        gene_index=gene_index,
    )
    cohort.assert_not_kl_mislabel()
    return cohort


def load_cohorts() -> list[Cohort]:
    cohorts = []

    path = ensure_file("GSE155972_Invivo_LLC_counts_GRCm38.txt.gz")
    samples, genes, mat = _delimited(path, 0, sniff_delimiter(path))
    cohorts.append(
        matrix_to_cohort(
            "GSE155972_LLC",
            "GSE155972",
            "LLC",
            "KrasG12C and NrasQ61; not Stk11/Lkb1",
            "in vivo subcutaneous LLC, anti-PD-1 + anti-CTLA-4",
            "tismo_label",
            True,
            "TISMO Responder is Setdb1-KO and Non-responder is control sgRNA. The label is the genotype.",
            samples,
            genes,
            mat,
            "count",
            gse155972_samples,
        )
    )

    path = ensure_file("GSE246922_KP_RNA_counts_vst.csv.gz")
    samples, genes, mat = _delimited(path, 0, ",")
    cohorts.append(
        matrix_to_cohort(
            "GSE246922_KP",
            "GSE246922",
            "KP",
            "Kras/Trp53 lung cancer cells; not Lkb1",
            "CD45-negative cells from parental KP and from tumors that relapsed after immunotherapy",
            "acquired_resistance",
            False,
            "",
            samples,
            genes,
            mat,
            "identity",
            gse246922_kp,
        )
    )

    path = ensure_file("GSE246922_LLC1_RNA_counts_vst.csv.gz")
    samples, genes, mat = _delimited(path, 0, ",")
    cohorts.append(
        matrix_to_cohort(
            "GSE246922_LLC1",
            "GSE246922",
            "LLC1",
            "KrasG12C and Nras; not Stk11/Lkb1",
            "CD45-negative LLC1 cells, parental versus relapsed after immunotherapy",
            "acquired_resistance",
            False,
            "",
            samples,
            genes,
            mat,
            "identity",
            gse246922_llc,
        )
    )

    path = ensure_file("GSE297630_processed_data.xlsx")
    samples, genes, mat = read_excel_matrix(path, "Expression", 4, 12, list(range(19, 25)))
    # shorten sample names to C-1 .. P-3
    short = []
    for name in samples:
        short.append("C-" + name.split("C-")[1][0] if name.startswith("C-") else "P-" + name.split("P-")[1][0])
    cohorts.append(
        matrix_to_cohort(
            "GSE297630_LLC",
            "GSE297630",
            "LLC",
            "KrasG12C and Nras; not Stk11/Lkb1",
            "subcutaneous LLC tumors; series title calls the anti-PD-1 arm tolerant cells",
            "tolerant_proxy",
            False,
            "Resistant arm is anti-PD-1 treated tumor, not a per-mouse shrinkage label. Reference arm was not treated.",
            short,
            genes,
            mat,
            "identity",
            gse297630_samples,
        )
    )

    path = ensure_file("GSE114601_counts.normalized.csv.gz")
    samples, genes, mat = _delimited(path, 0, ",")
    cohorts.append(
        matrix_to_cohort(
            "GSE114601_KP",
            "GSE114601",
            "KP",
            "Kras/Trp53 GEMM lung tumor; not Lkb1",
            "KP lung nodules; anti-PD-1, JQ1, combination, vehicle. No per-mouse R/NR label.",
            "none",
            False,
            "",
            samples,
            genes,
            mat,
            "positive",
            gse114601_samples,
        )
    )

    path = ensure_file("GSE157880_Bulk048.txt.gz")
    samples, genes, mat = _delimited(path, 4, sniff_delimiter(path), sample_start=8)
    cohorts.append(
        matrix_to_cohort(
            "GSE157880_HKP1",
            "GSE157880",
            "HKP1",
            "KrasG12D/Trp53-null lung line; not Lkb1",
            "HKP1-bearing lungs, anti-PD-1 or IgG, with or without radiation. No per-mouse R/NR label.",
            "none",
            False,
            "",
            samples,
            genes,
            mat,
            "positive",
            lambda cols: [SampleAnn(c, GSE157880_GROUP[c], "other") for c in cols],
        )
    )

    path = ensure_file("GSE274960_gene_count_matrix.csv.gz")
    samples, genes, mat = _delimited(path, 0, ",")
    cohorts.append(
        matrix_to_cohort(
            "GSE274960_LL2",
            "GSE274960",
            "LL/2",
            "LL/2 is Lewis lung, KrasG12C; not Stk11/Lkb1",
            "LL/2 tumors, IgG, anti-PD-1, entrectinib, or the combination. shRNA arms are excluded.",
            "none",
            False,
            "",
            samples,
            genes,
            mat,
            "count",
            gse274960_samples,
        )
    )
    # shCTRL / shNTRK1 is a second experiment with its own baseline. Keep the
    # IgG / entrectinib / anti-PD-1 arms that belong to the ICB study.
    drop_sample_groups(cohorts[-1], {"shctrl", "shntrk1"})

    path = ensure_file("GSE262305_gene_expression.txt.gz")
    # FPKM columns only
    delim = sniff_delimiter(path)
    with open_text(path) as handle:
        header = next(csv.reader(handle, delimiter=delim))
    fpkm_cols = [i for i, h in enumerate(header) if h.endswith("_FPKM")]
    samples, genes, mat = _delimited(path, 1, delim, sample_cols=fpkm_cols)
    cohorts.append(
        matrix_to_cohort(
            "GSE262305_LLC",
            "GSE262305",
            "LLC",
            "KrasG12C and Nras; not Stk11/Lkb1",
            "subcutaneous LLC, isotype, anti-PD-L1, or bortezomib plus anti-PD-L1. No R/NR label.",
            "none",
            False,
            "",
            samples,
            genes,
            mat,
            "positive",
            gse262305_samples,
        )
    )

    path = ensure_file("GSE239485_Processed_data.xlsx")
    samples, genes, mat = read_excel_matrix(path, "DataNorm", 0, 2, list(range(11, 35)))
    cohorts.append(
        matrix_to_cohort(
            "GSE239485_LLC",
            "GSE239485",
            "LLC",
            "KrasG12C and Nras; not Stk11/Lkb1",
            "LLC tumors. Treated arms are Poly I:C plus anti-PD-1, with or without anti-C5aR1. No monotherapy and no R/NR label.",
            "none",
            False,
            "",
            samples,
            genes,
            mat,
            "identity",
            gse239485_samples,
        )
    )

    path = ensure_file("GSE260596_AllSamples_Genes_ReadCounts.txt.gz")
    samples, genes, mat = _delimited(path, 0, sniff_delimiter(path))
    cohorts.append(
        matrix_to_cohort(
            "GSE260596_KP",
            "GSE260596",
            "KP",
            "refractory KP mouse lung tumors (series summary); not Lkb1",
            "All profiled tumors are on anti-PD-1, with control or anti-LAIR1. No untreated arm and no per-mouse R/NR label.",
            "none",
            False,
            "",
            samples,
            genes,
            mat,
            "count",
            gse260596_samples,
        )
    )
    return cohorts


def _delimited(path: Path, gene_col: int, delimiter: str, sample_start: int | None = 1, sample_cols: list[int] | None = None):
    with open_text(path) as handle:
        header = next(csv.reader(handle, delimiter=delimiter))
    if sample_cols is None:
        sample_cols = list(range(sample_start, len(header)))
    samples, names, mat = full_matrix_from_delimited(path, gene_col, sample_cols, delimiter)
    return samples, names, mat


# ---------------------------------------------------------------------------
# Scores and sweep
# ---------------------------------------------------------------------------

def module_indices(cohort: Cohort, module: str) -> list[int]:
    wanted = MODULES[module]
    found = [cohort.gene_index[g] for g in wanted if g in cohort.gene_index]
    if len(found) < max(2, math.ceil(MIN_MODULE_FRACTION * len(wanted))):
        return []
    return found


def score_module(mat: np.ndarray, gene_rows: list[int], cols: np.ndarray, method: str) -> np.ndarray:
    sub = mat[:, cols]
    if method == "pct":
        # Mean within-sample percentile of module genes. Rank is invariant to
        # a per-sample monotone rescaling, so count and CPM share this score.
        n_genes = sub.shape[0]
        out = np.empty(len(cols), dtype=np.float64)
        mod = np.asarray(gene_rows)
        for j in range(sub.shape[1]):
            column = sub[:, j]
            ranks = stats.rankdata(np.nan_to_num(column, nan=-1e30), method="average")
            out[j] = float(np.mean(ranks[mod] / n_genes))
        return out
    block = sub[gene_rows, :]
    if method == "mean":
        return np.nanmean(block, axis=0)
    if method == "zmean":
        mu = np.nanmean(block, axis=1, keepdims=True)
        sd = np.nanstd(block, axis=1, keepdims=True)
        sd[sd == 0] = np.nan
        z = (block - mu) / sd
        return np.nanmean(z, axis=0)
    raise RuntimeError(method)


def response_masks(cohort: Cohort, resistant_groups: set[str] | None = None) -> tuple[np.ndarray, np.ndarray]:
    resistant = []
    reference = []
    for i, sample in enumerate(cohort.samples):
        if sample.role == "resistant" and (resistant_groups is None or sample.group in resistant_groups):
            resistant.append(i)
        elif sample.role == "reference" and resistant_groups is None:
            reference.append(i)
        elif resistant_groups is not None and sample.role == "reference":
            reference.append(i)
    return np.asarray(resistant, dtype=int), np.asarray(reference, dtype=int)


def filter_changes(before: np.ndarray, after: np.ndarray) -> bool:
    if len(before) != len(after):
        return True
    return not np.array_equal(np.sort(before), np.sort(after))


def apply_filter(cohort: Cohort, cols: np.ndarray, cldn4: np.ndarray, filt: str, transform: str) -> np.ndarray:
    if filt == "all":
        return cols
    if filt == "drop_cldn4_floor":
        if transform not in {"log2_count", "log2_cpm", "log2_value"}:
            return cols
        keep = cldn4[cols] > 0
        if keep.sum() == len(cols) or keep.sum() < 4:
            return cols
        return cols[keep]
    if filt == "drop_max_cldn4":
        if len(cols) < 5:
            return cols
        local = cldn4[cols]
        if np.nanmax(local) == np.nanmin(local):
            return cols
        drop_local = int(np.nanargmax(local))
        return np.concatenate([cols[:drop_local], cols[drop_local + 1 :]])
    raise RuntimeError(filt)


def extreme_split(values: np.ndarray, how: str) -> tuple[np.ndarray, np.ndarray]:
    """Return boolean masks high, low inside the provided vector."""
    n = len(values)
    order = np.argsort(values, kind="mergesort")
    if how == "median":
        mid = n // 2
        low_idx = order[:mid]
        high_idx = order[mid:]
        # If odd, median sample stays out of both tails when n is odd? Keeping
        # the upper half including the middle is the usual median split.
    elif how == "tertile":
        a = n // 3
        low_idx = order[:a]
        high_idx = order[-a:]
    elif how == "quartile":
        a = n // 4
        low_idx = order[:a]
        high_idx = order[-a:]
    else:
        raise RuntimeError(how)
    high = np.zeros(n, dtype=bool)
    low = np.zeros(n, dtype=bool)
    high[high_idx] = True
    low[low_idx] = True
    # Overlap can happen with heavy ties. Drop overlap.
    both = high & low
    high[both] = False
    low[both] = False
    return high, low


def record(**kwargs) -> dict:
    return kwargs


def sweep_cohort(cohort: Cohort) -> list[dict]:
    rows = []
    cldn4_row = cohort.gene_index["Cldn4"]
    base_cols = np.arange(len(cohort.samples))

    # Response contrasts
    if cohort.response_family != "none":
        group_sets = [("pooled", None)]
        resistant_groups = sorted({s.group for s in cohort.samples if s.role == "resistant"})
        if len(resistant_groups) > 1:
            for group in resistant_groups:
                group_sets.append((group, {group}))
        for transform, mat in cohort.matrices.items():
            cldn4 = mat[cldn4_row]
            for subset_name, groups in group_sets:
                res_idx, ref_idx = response_masks(cohort, groups)
                if len(res_idx) == 0 or len(ref_idx) == 0:
                    continue
                for filt in ("all", "drop_cldn4_floor", "drop_max_cldn4"):
                    res_f = apply_filter(cohort, res_idx, cldn4, filt, transform)
                    ref_f = apply_filter(cohort, ref_idx, cldn4, filt, transform)
                    # apply_filter on each arm separately can drop the max inside that arm
                    if filt == "drop_max_cldn4":
                        both = np.concatenate([res_idx, ref_idx])
                        kept = apply_filter(cohort, both, cldn4, filt, transform)
                        res_f = np.intersect1d(res_idx, kept)
                        ref_f = np.intersect1d(ref_idx, kept)
                    if filt != "all" and not filter_changes(
                        np.concatenate([res_idx, ref_idx]), np.concatenate([res_f, ref_f])
                    ):
                        continue
                    high = cldn4[res_f]
                    low = cldn4[ref_f]
                    effect = median_diff(high, low)
                    eligible = len(res_f) >= MIN_ARM and len(ref_f) >= MIN_ARM
                    for method, p in (
                        ("mwu", safe_mwu(high, low)),
                        ("welch", safe_welch(high, low)),
                    ):
                        rows.append(
                            _response_row(
                                cohort, subset_name, transform, filt, method, effect, p, len(res_f), len(ref_f), eligible
                            )
                        )
                    # binary spearman, resistant = 1
                    y = np.concatenate([np.ones(len(res_f)), np.zeros(len(ref_f))])
                    x = np.concatenate([high, low])
                    rho, p = safe_spearman(x, y)
                    rows.append(
                        _response_row(
                            cohort,
                            subset_name,
                            transform,
                            filt,
                            "spearman_binary",
                            rho,
                            p,
                            len(res_f),
                            len(ref_f),
                            eligible and len(x) >= MIN_CORR_N,
                            effect_name="spearman_rho",
                        )
                    )
                    for split in ("median", "tertile", "quartile"):
                        if len(x) < 8:
                            continue
                        hi_m, lo_m = extreme_split(x, split)
                        # Fisher: high Cldn4 vs resistant
                        a = int(np.sum(hi_m & (y == 1)))
                        b = int(np.sum(hi_m & (y == 0)))
                        c = int(np.sum(lo_m & (y == 1)))
                        d = int(np.sum(lo_m & (y == 0)))
                        if min(a + b, c + d, a + c, b + d) == 0 or (a + b) < MIN_ARM or (c + d) < MIN_ARM:
                            continue
                        oddsr, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
                        # Thesis: high Cldn4 enriched in resistant, OR > 1
                        rows.append(
                            _response_row(
                                cohort,
                                subset_name,
                                transform,
                                filt,
                                f"fisher_{split}",
                                float(oddsr),
                                float(p),
                                len(res_f),
                                len(ref_f),
                                eligible,
                                effect_name="odds_ratio_high_cldn4_in_resistant",
                                aligned_override=(float(oddsr) > 1),
                            )
                        )
        # ordinal, one row per transform, samples with an ordinal
        for transform, mat in cohort.matrices.items():
            cldn4 = mat[cldn4_row]
            idx = [i for i, s in enumerate(cohort.samples) if s.ordinal is not None]
            if len(idx) < MIN_CORR_N:
                continue
            cols = np.asarray(idx)
            x = cldn4[cols]
            y = np.asarray([cohort.samples[i].ordinal for i in idx], dtype=float)
            rho, p = safe_spearman(x, y)
            rows.append(
                _response_row(
                    cohort,
                    "ordinal_stage",
                    transform,
                    "all",
                    "spearman_ordinal",
                    rho,
                    p,
                    int(np.sum(y > 0)),
                    int(np.sum(y == 0)),
                    True,
                    effect_name="spearman_rho",
                )
            )

    # Module correlations on every cohort
    for transform, mat in cohort.matrices.items():
        cldn4 = mat[cldn4_row]
        for filt in ("all", "drop_cldn4_floor", "drop_max_cldn4"):
            cols = apply_filter(cohort, base_cols, cldn4, filt, transform)
            if filt != "all" and not filter_changes(base_cols, cols):
                continue
            if len(cols) < 6:
                continue
            x = cldn4[cols]
            for module in MODULES:
                gene_rows = module_indices(cohort, module)
                if not gene_rows:
                    continue
                sign = MODULE_THESIS_SIGN[module]
                question = "cldn4_ifn" if module in IFN_MODULES else "cldn4_nhej"
                for score_name in ("zmean", "mean", "pct"):
                    # pct does not depend on the expression transform (monotone
                    # within sample). Emit it once, on the first transform.
                    if score_name == "pct" and transform != next(iter(cohort.matrices)):
                        continue
                    y = score_module(mat, gene_rows, cols, score_name)
                    for assoc, (effect, p) in (
                        ("spearman", safe_spearman(x, y)),
                        ("pearson", safe_pearson(x, y)),
                    ):
                        rows.append(
                            _module_row(
                                cohort,
                                question,
                                module,
                                score_name,
                                assoc,
                                transform if score_name != "pct" else "rank_within_sample",
                                filt,
                                effect,
                                p,
                                len(cols),
                                len(gene_rows),
                                sign,
                                "continuous",
                            )
                        )
                    for split in ("median", "tertile", "quartile"):
                        hi_m, lo_m = extreme_split(x, split)
                        if hi_m.sum() < MIN_ARM or lo_m.sum() < MIN_ARM:
                            continue
                        effect = median_diff(y[hi_m], y[lo_m])
                        p = safe_mwu(y[hi_m], y[lo_m])
                        aligned = finite(effect) and (effect * sign > 0)
                        rows.append(
                            _module_row(
                                cohort,
                                question,
                                module,
                                score_name,
                                f"mwu_{split}",
                                transform if score_name != "pct" else "rank_within_sample",
                                filt,
                                effect,
                                p,
                                len(cols),
                                len(gene_rows),
                                sign,
                                split,
                                aligned_override=aligned,
                                n_high=int(hi_m.sum()),
                                n_low=int(lo_m.sum()),
                            )
                        )
    return rows


def _response_row(
    cohort: Cohort,
    subset: str,
    transform: str,
    filt: str,
    method: str,
    effect: float,
    p: float,
    n_resistant: int,
    n_reference: int,
    eligible: bool,
    effect_name: str = "median_resistant_minus_reference",
    aligned_override: bool | None = None,
) -> dict:
    if aligned_override is None:
        aligned = finite(effect) and effect > 0
    else:
        aligned = aligned_override
    headline_ok = (
        eligible
        and aligned
        and filt == "all"
        and (not cohort.confounded)
        and cohort.response_family in {"acquired_resistance", "tolerant_proxy"}
        and finite(p)
    )
    return {
        "question": "cldn4_response",
        "cohort": cohort.cohort_id,
        "accession": cohort.accession,
        "model": cohort.model,
        "kras": cohort.kras,
        "response_family": cohort.response_family,
        "confounded": cohort.confounded,
        "subset": subset,
        "module": "",
        "score": "",
        "association": method,
        "transform": transform,
        "sample_filter": filt,
        "effect_name": effect_name,
        "effect": effect,
        "p": p,
        "n": n_resistant + n_reference,
        "n_resistant": n_resistant,
        "n_reference": n_reference,
        "n_genes": "",
        "thesis_aligned": aligned,
        "eligible_headline": headline_ok,
        "pre_specified": subset == "pooled" and filt == "all" and method == "mwu" and transform in {"author_scale", "log2_cpm", "log2_value"},
    }


def _module_row(
    cohort,
    question,
    module,
    score_name,
    assoc,
    transform,
    filt,
    effect,
    p,
    n,
    n_genes,
    sign,
    split,
    aligned_override=None,
    n_high="",
    n_low="",
) -> dict:
    if aligned_override is None:
        aligned = finite(effect) and effect * sign > 0
    else:
        aligned = aligned_override
    # Headline may be a threshold split or an alternate module, but not a
    # sample-dropping filter. The response-label confound does not apply here.
    split_ok = True
    if assoc.startswith("mwu_"):
        split_ok = isinstance(n_high, int) and isinstance(n_low, int) and n_high >= MIN_ARM and n_low >= MIN_ARM
    headline_ok = (
        aligned
        and n >= MIN_CORR_N
        and finite(p)
        and filt == "all"
        and split_ok
        and score_name in {"zmean", "mean", "pct"}
        and (
            assoc in {"spearman", "pearson"}
            or assoc.startswith("mwu_")
        )
    )
    return {
        "question": question,
        "cohort": cohort.cohort_id,
        "accession": cohort.accession,
        "model": cohort.model,
        "kras": cohort.kras,
        "response_family": cohort.response_family,
        "confounded": False,
        "subset": split,
        "module": module,
        "score": score_name,
        "association": assoc,
        "transform": transform,
        "sample_filter": filt,
        "effect_name": "rho" if assoc in {"spearman", "pearson"} else "median_high_minus_low",
        "effect": effect,
        "p": p,
        "n": n,
        "n_resistant": n_high,
        "n_reference": n_low,
        "n_genes": n_genes,
        "thesis_aligned": aligned,
        "eligible_headline": headline_ok,
        "pre_specified": (
            filt == "all"
            and module in {"nhej_core", "ifn_chemokine"}
            and score_name == "zmean"
            and assoc == "spearman"
            and transform in {"author_scale", "log2_cpm", "log2_value"}
        ),
    }


def add_q(rows: list[dict]) -> None:
    for question in {r["question"] for r in rows}:
        idx = [i for i, r in enumerate(rows) if r["question"] == question]
        qvals = bh([rows[i]["p"] for i in idx])
        for i, q in zip(idx, qvals):
            rows[i]["q_bh"] = q


def practical_effect(row: dict) -> bool:
    """Drop gaps that are numerically tiny on their own scale."""
    effect = row["effect"]
    if not finite(effect):
        return False
    name = row["effect_name"]
    magnitude = abs(effect)
    if name in {"rho", "spearman_rho"}:
        return magnitude >= 0.30
    if "odds" in name:
        return effect >= 1.5
    if name == "median_high_minus_low":
        if row.get("score") == "pct":
            return magnitude >= 0.02
        return magnitude >= 0.20
    if name == "median_resistant_minus_reference":
        return magnitude >= 0.25
    return magnitude > 0


def pick_winner(rows: list[dict], question: str) -> dict | None:
    pool = [
        r
        for r in rows
        if r["question"] == question and r["eligible_headline"] and r["thesis_aligned"] and practical_effect(r)
    ]
    if not pool:
        return None
    pool.sort(key=lambda r: (r["p"], -abs(r["effect"]) if finite(r["effect"]) else 0))
    return pool[0]


def pre_specified(rows: list[dict]) -> list[dict]:
    return [r for r in rows if r.get("pre_specified")]


def fmt_p(p: float) -> str:
    if not finite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_e(x: float) -> str:
    if not finite(x):
        return "NA"
    return f"{x:.3f}"


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        for row in rows:
            out = {}
            for key in fields:
                val = row.get(key, "")
                if isinstance(val, float):
                    out[key] = f"{val:.6g}" if finite(val) else ""
                elif isinstance(val, bool):
                    out[key] = "1" if val else "0"
                else:
                    out[key] = val
            writer.writerow(out)


def drop_sample_groups(cohort: Cohort, groups: set[str]) -> None:
    keep = [i for i, sample in enumerate(cohort.samples) if sample.group not in groups]
    if len(keep) == len(cohort.samples):
        return
    cohort.samples = [cohort.samples[i] for i in keep]
    for key, mat in list(cohort.matrices.items()):
        cohort.matrices[key] = mat[:, keep]


def sample_table(cohorts: list[Cohort]) -> list[dict]:
    rows = []
    for cohort in cohorts:
        # Prefer CPM for counts, otherwise the only scale
        transform = "log2_cpm" if "log2_cpm" in cohort.matrices else next(iter(cohort.matrices))
        mat = cohort.matrices[transform]
        cldn4 = mat[cohort.gene_index["Cldn4"]]
        tac = cohort.gene_index.get("Tacstd2")
        nhej = module_indices(cohort, "nhej_core")
        ifn = module_indices(cohort, "ifn_chemokine")
        cols = np.arange(len(cohort.samples))
        nhej_s = score_module(mat, nhej, cols, "zmean") if nhej else np.full(len(cols), np.nan)
        ifn_s = score_module(mat, ifn, cols, "zmean") if ifn else np.full(len(cols), np.nan)
        for i, sample in enumerate(cohort.samples):
            rows.append(
                {
                    "cohort": cohort.cohort_id,
                    "accession": cohort.accession,
                    "model": cohort.model,
                    "kras": cohort.kras,
                    "is_kl": 0,
                    "sample": sample.name,
                    "group": sample.group,
                    "role": sample.role,
                    "ordinal": sample.ordinal if sample.ordinal is not None else "",
                    "transform": transform,
                    "cldn4": cldn4[i],
                    "tacstd2": mat[tac, i] if tac is not None else "",
                    "nhej_core_zmean": nhej_s[i],
                    "ifn_chemokine_zmean": ifn_s[i],
                }
            )
    return rows


def hkp1_stratum_note(cohorts: list[Cohort]) -> str:
    cohort = next((c for c in cohorts if c.cohort_id == "GSE157880_HKP1"), None)
    if cohort is None:
        return ""
    mat = cohort.matrices["log2_value"]
    cldn4 = mat[cohort.gene_index["Cldn4"]]
    extended = module_indices(cohort, "nhej_extended")
    core = module_indices(cohort, "nhej_core")
    cols = np.arange(len(cohort.samples))
    ext = score_module(mat, extended, cols, "pct")
    cor = score_module(mat, core, cols, "pct")
    rho_core, p_core = safe_spearman(cldn4, cor)

    def subset(pred):
        idx = np.array([i for i, sample in enumerate(cohort.samples) if pred(sample)])
        rho, p = safe_spearman(cldn4[idx], ext[idx])
        return len(idx), rho, p

    n_igg, rho_igg, p_igg = subset(lambda s: s.group.startswith("igg"))
    n_pd, rho_pd, p_pd = subset(lambda s: s.group.startswith("pd1"))
    n0, rho0, p0 = subset(lambda s: s.group.endswith("0gy"))
    nrt, rho_rt, p_rt = subset(lambda s: not s.group.endswith("0gy"))
    return (
        "The HKP1 BH rows are the extended NHEJ neighborhood (53BP1, ATM, PARP1, and the other genes outside the c-NHEJ core), "
        f"scored as a within-sample percentile. c-NHEJ core on the same 16 samples is Spearman ρ={fmt_e(rho_core)}, p={fmt_p(p_core)}. "
        f"The extended-set correlation is in the IgG lungs (ρ={fmt_e(rho_igg)}, n={n_igg}, p={fmt_p(p_igg)}) and in the anti-PD-1 lungs "
        f"(ρ={fmt_e(rho_pd)}, n={n_pd}, p={fmt_p(p_pd)}). 0 Gy alone is n={n0}, ρ={fmt_e(rho0)}, p={fmt_p(p0)}. "
        f"Radiated lungs are n={nrt}, ρ={fmt_e(rho_rt)}, p={fmt_p(p_rt)}."
    )


def write_results_md(cohorts: list[Cohort], rows: list[dict], winners: dict[str, dict | None]) -> str:
    def line_for(row: dict | None) -> str:
        if row is None:
            return "No thesis-aligned statistic met the reporting rule."
        span = ""
        if finite(row.get("cldn4_min", float("nan"))) and finite(row.get("cldn4_max", float("nan"))):
            width = row["cldn4_max"] - row["cldn4_min"]
            span = f" Cldn4 span on {row.get('cldn4_scale', 'the analysis scale')} is {fmt_e(row['cldn4_min'])} to {fmt_e(row['cldn4_max'])} (width {fmt_e(width)})."
            if width < 1:
                span += " That width is a near-floor gene, so the correlation is carried by a very small expression range."
        if finite(row.get("score_width", float("nan"))):
            span += f" The module score itself spans {fmt_e(row['score_width'])}."
        return (
            f"{row['cohort']} ({row['model']}), {row['association']} / {row['transform']}"
            f"{(' / ' + row['module'] + ' ' + row['score']) if row['module'] else ''}"
            f"{(' / subset ' + row['subset']) if row['subset'] not in {'', 'continuous', 'pooled'} else ''}"
            f", filter {row['sample_filter']}: effect {fmt_e(row['effect'])} ({row['effect_name']}), "
            f"p {fmt_p(row['p'])}, BH q {fmt_p(row['q_bh'])}, n {row['n']}."
            f"{span}"
        )

    n_by_q = {}
    aligned_by_q = {}
    for row in rows:
        n_by_q[row["question"]] = n_by_q.get(row["question"], 0) + 1
        if row["thesis_aligned"]:
            aligned_by_q[row["question"]] = aligned_by_q.get(row["question"], 0) + 1

    prespec = [r for r in rows if r["question"] == "cldn4_response" and r["pre_specified"] and r["association"] == "mwu" and r["subset"] == "pooled" and r["sample_filter"] == "all"]
    # de-duplicate transforms: keep one primary scale per cohort
    seen = set()
    prespec_unique = []
    for row in prespec:
        if row["cohort"] in seen:
            continue
        # prefer author_scale, else log2_cpm, else log2_value
        if row["transform"] not in {"author_scale", "log2_cpm", "log2_value"}:
            continue
        seen.add(row["cohort"])
        prespec_unique.append(row)

    lines = []
    lines.append("# Kras-mutant mouse lung ICB: Cldn4 vs response and vs NHEJ/IFN")
    lines.append("")
    lines.append("Public matrices only. TISMO Tacstd2 49/64 is not recomputed. No series here is a KL (Kras/Lkb1) ICB experiment. LLC and LL/2 are KrasG12C/Nras Lewis lung, not KL. CMT-167 (KrasG12V) has no ICB rows in TISMO. MLE12 is not Kras.")
    lines.append("")
    lines.append("## Strongest thesis-aligned statistics")
    lines.append("")
    lines.append("A headline row has the pre-set sign, no sample dropped, and a non-trivial effect (|ρ|≥0.30, percentile-score gap ≥0.02, other score gaps ≥0.20, resistant-minus-reference median gap ≥0.25, odds ratio ≥1.5). Response headlines also exclude the Setdb1-confounded TISMO label and need n≥3 per arm. Module headlines need n≥8. They may use any panel in that question, any score, Spearman, Pearson, or a median/tertile/quartile split. BH q is computed inside each question across every row of the sweep, aligned or not.")
    lines.append("")
    lines.append(f"- **Cldn4 vs response.** {line_for(winners.get('cldn4_response'))}")
    lines.append(f"- **Cldn4 vs IFN.** {line_for(winners.get('cldn4_ifn'))}")
    lines.append(f"- **Cldn4 vs NHEJ.** {line_for(winners.get('cldn4_nhej'))}")
    lines.append("")
    lines.append("Sweep size: " + ", ".join(f"{q} {n_by_q[q]} tests ({aligned_by_q.get(q, 0)} thesis-aligned)" for q in sorted(n_by_q)) + ".")
    lines.append("")
    surviving = [
        r
        for r in rows
        if finite(r.get("q_bh", float("nan"))) and r["q_bh"] < 0.05 and r["sample_filter"] == "all"
    ]
    surviving.sort(key=lambda r: (r["q_bh"], r["p"]))
    lines.append("## Rows with BH q < 0.05")
    lines.append("")
    if not surviving:
        lines.append("No row in the sweep has BH q < 0.05.")
    else:
        lines.append("These are multiplicity-adjusted hits on the full sample set. The sign column says whether the hit matches the pre-set thesis direction.")
        lines.append("")
        lines.append("| Question | Cohort | Test | Effect | p | q | n | Thesis sign |")
        lines.append("|---|---|---|---:|---:|---:|---:|---|")
        for row in surviving:
            label = row["association"]
            if row["module"]:
                label = f"{row['module']} {row['score']} {row['association']}"
            sign = "matches" if row["thesis_aligned"] else "opposite"
            lines.append(
                f"| {row['question']} | {row['cohort']} | {label} | {fmt_e(row['effect'])} | {fmt_p(row['p'])} | {fmt_p(row['q_bh'])} | {row['n']} | {sign} |"
            )
    note = hkp1_stratum_note(cohorts)
    if note:
        lines.append("")
        lines.append(note)
    lines.append("")
    lines.append("## Pre-specified module correlations")
    lines.append("")
    lines.append("Spearman of Cldn4 with the z-mean of c-NHEJ core or the IFN/chemokine panel, all samples, one expression scale per cohort. Thesis sign is positive for NHEJ and negative for IFN.")
    lines.append("")
    lines.append("| Cohort | Module | n | ρ | p | q | Sign |")
    lines.append("|---|---|---:|---:|---:|---:|---|")
    seen_mod = set()
    for row in rows:
        if not (
            row["pre_specified"]
            and row["association"] == "spearman"
            and row["score"] == "zmean"
            and row["sample_filter"] == "all"
            and row["module"] in {"nhej_core", "ifn_chemokine"}
            and row["transform"] in {"author_scale", "log2_cpm", "log2_value"}
        ):
            continue
        key = (row["cohort"], row["module"])
        if key in seen_mod:
            continue
        seen_mod.add(key)
        sign = "matches" if row["thesis_aligned"] else "opposite"
        lines.append(
            f"| {row['cohort']} | {row['module']} | {row['n']} | {fmt_e(row['effect'])} | {fmt_p(row['p'])} | {fmt_p(row['q_bh'])} | {sign} |"
        )
    lines.append("")
    lines.append("## Pre-specified response contrasts")
    lines.append("")
    lines.append("These are Mann-Whitney tests on the pooled resistant-versus-reference contrast, all samples, one expression scale. They are reported whether or not they win the sweep.")
    lines.append("")
    lines.append("| Cohort | Model | What the arms are | n resistant vs reference | median difference | p | BH q | sign |")
    lines.append("|---|---|---|---:|---:|---:|---:|---|")
    for row in prespec_unique:
        sign = "resistant higher" if row["thesis_aligned"] else "resistant lower or flat"
        note = "Setdb1 genotype, not an independent response call" if row["confounded"] else row["response_family"].replace("_", " ")
        lines.append(
            f"| {row['cohort']} | {row['model']} | {note} | {row['n_resistant']} vs {row['n_reference']} | {fmt_e(row['effect'])} | {fmt_p(row['p'])} | {fmt_p(row['q_bh'])} | {sign} |"
        )
    lines.append("")
    n_high = sum(1 for row in prespec_unique if row["thesis_aligned"])
    n_low = len(prespec_unique) - n_high
    lines.append(
        f"Pre-specified Mann-Whitney sign: Cldn4 is higher in the resistant arm in {n_high}/{len(prespec_unique)} and lower in {n_low}/{len(prespec_unique)}."
    )
    lines.append("")
    lines.append("GSE155972 is in that table and is not eligible for the response headline. TISMO calls Setdb1-KO ICB mice Responders and control-sgRNA ICB mice Non-responders.")
    lines.append("")
    lines.append("## Cohorts")
    lines.append("")
    lines.append("| Cohort | Model | Kras | Response label | Setting |")
    lines.append("|---|---|---|---|---|")
    for cohort in cohorts:
        family = "none (modules and treatment context only)" if cohort.response_family == "none" else cohort.response_family
        if cohort.confounded:
            family += "; confounded"
        lines.append(f"| {cohort.cohort_id} | {cohort.model} | {cohort.kras} | {family} | {cohort.setting} |")
    lines.append("")
    lines.append("## Screened and not scored")
    lines.append("")
    lines.append("- TISMO CMT-167: KrasG12V lung carcinoma, no ICB expression rows.")
    lines.append("- TISMO MLE12: lung adenocarcinoma driven by SV40 large T, not Kras, and no ICB rows.")
    lines.append("- TISMO KPB25L: Kras/p53 mammary, not lung.")
    lines.append("- GSE193895: Kras/Keap1/Lkb1 GEMM tumors, deposited matrix has no ICB arm.")
    lines.append("- GSE277929: KL tumors treated with entinostat and trametinib, not ICB.")
    lines.append("- GSE137244: locked KL-versus-KP cell-line result, not an ICB experiment, not rerun.")
    lines.append("")
    lines.append("## Scales")
    lines.append("")
    lines.append("Count matrices are tested as log2(count+1) and log2(CPM+1). FPKM and normalized counts use log2(value+1). Author VST, RMA, and the GSE239485 processed matrix are used as deposited. Within-sample percentile scores do not depend on a monotone per-sample rescaling, so they are emitted once per cohort.")
    lines.append("")
    lines.append("Chronic IFNG arms in GSE246922 stay out of the resistant group. They are in the module correlations.")
    lines.append("")
    return "\n".join(lines) + "\n"


def plot(cohorts: list[Cohort], rows: list[dict], winners: dict[str, dict | None], sample_rows: list[dict]) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    by_cohort = {c.cohort_id: c for c in cohorts}

    # Response strip plots for the three pre-specified cohorts
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.6), constrained_layout=True)
    specs = [
        ("GSE246922_KP", "KP CD45− relapse", ["parental", "relapse_1", "relapse_2", "relapse_late"]),
        ("GSE246922_LLC1", "LLC1 CD45− relapse", ["parental", "relapsed"]),
        ("GSE155972_LLC", "LLC ICB (Setdb1 label)", ["control_icb", "setdb1_ko_icb"]),
    ]
    for ax, (cid, title, groups) in zip(axes, specs):
        sub = [r for r in sample_rows if r["cohort"] == cid and r["group"] in groups]
        # map group to x
        xpos = {g: i for i, g in enumerate(groups)}
        rng = np.random.default_rng(0)
        for row in sub:
            x = xpos[row["group"]] + float(rng.uniform(-0.12, 0.12))
            color = "#b2182b" if row["role"] == "resistant" else "#2166ac"
            ax.scatter(x, row["cldn4"], s=22, color=color, zorder=3)
        for group, i in xpos.items():
            vals = [float(row["cldn4"]) for row in sub if row["group"] == group]
            if vals:
                ax.hlines(float(np.median(vals)), i - 0.28, i + 0.28, color="black", linewidth=1.4, zorder=2)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([g.replace("_", "\n") for g in groups], fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("Cldn4" if ax is axes[0] else "")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if cid == "GSE155972_LLC":
            hi = [row for row in sub if float(row["cldn4"]) > 1.2]
            ax.set_ylim(-0.05, 1.15)
            if hi:
                ax.text(
                    0.02,
                    0.98,
                    f"one Setdb1 ICB mouse at {float(hi[0]['cldn4']):.2f}, off scale",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=7,
                )
    fig.suptitle("Cldn4 on pre-specified Kras-lung contrasts", fontsize=12)
    fig.savefig(fig_dir / "response_contrasts.png", dpi=160)
    plt.close(fig)

    # Winner module scatters
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), constrained_layout=True)
    for ax, question, ylab in (
        (axes[0], "cldn4_ifn", "IFN score"),
        (axes[1], "cldn4_nhej", "NHEJ score"),
    ):
        winner = winners.get(question)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if winner is None:
            ax.set_title("No aligned result")
            continue
        cohort = by_cohort[winner["cohort"]]
        transform = winner["transform"]
        if transform == "rank_within_sample":
            transform = next(iter(cohort.matrices))
        mat = cohort.matrices[transform]
        cldn4 = mat[cohort.gene_index["Cldn4"]]
        cols = np.arange(len(cohort.samples))
        cols = apply_filter(cohort, cols, cldn4, winner["sample_filter"], transform)
        y = score_module(mat, module_indices(cohort, winner["module"]), cols, winner["score"])
        ax.scatter(cldn4[cols], y, s=22, color="#4d4d4d")
        ax.set_xlabel("Cldn4")
        ax.set_ylabel(ylab)
        span = ""
        if finite(winner.get("cldn4_min", float("nan"))):
            width = winner["cldn4_max"] - winner["cldn4_min"]
            if width < 1:
                span = "\nCldn4 near the floor"
        score_width = float(np.nanmax(y) - np.nanmin(y))
        ax.set_title(
            f"{winner['cohort']}\n{fmt_e(winner['effect'])}, p {fmt_p(winner['p'])}, q {fmt_p(winner['q_bh'])}"
            f"\nscore width {score_width:.3f}{span}",
            fontsize=9,
        )
    fig.suptitle("Strongest thesis-aligned module correlations", fontsize=12)
    fig.savefig(fig_dir / "module_winners.png", dpi=160)
    plt.close(fig)


def inventory_rows() -> list[dict]:
    return [
        {"source": "TISMO", "model": "LLC", "kras": "KrasG12C/Nras", "icb": "yes", "scored": "yes", "note": "Only TISMO lung line with ICB rows. Response label equals Setdb1 genotype."},
        {"source": "TISMO", "model": "CMT-167", "kras": "KrasG12V", "icb": "no", "scored": "no", "note": "No ICB expression rows."},
        {"source": "TISMO", "model": "MLE12", "kras": "not Kras (SV40 large T)", "icb": "no", "scored": "no", "note": "Not a Kras lung line."},
        {"source": "TISMO", "model": "KPB25L", "kras": "Kras/p53", "icb": "yes", "scored": "no", "note": "Mammary, not lung."},
        {"source": "GSE246922", "model": "KP and LLC1", "kras": "Kras/Trp53; LLC1 KrasG12C", "icb": "acquired resistance", "scored": "yes", "note": "CD45- cells. Chronic IFNG is not the resistant arm."},
        {"source": "GSE297630", "model": "LLC", "kras": "KrasG12C/Nras", "icb": "anti-PD-1 tolerant vs untreated", "scored": "yes", "note": "Proxy, not per-mouse response."},
        {"source": "GSE114601", "model": "KP GEMM", "kras": "Kras/Trp53", "icb": "anti-PD-1 vs vehicle", "scored": "modules", "note": "n=2 per monotherapy arm. No R/NR."},
        {"source": "GSE157880", "model": "HKP1", "kras": "KrasG12D/Trp53", "icb": "anti-PD-1 ± radiation", "scored": "modules", "note": "No R/NR. Monotherapy 0 Gy is n=2 vs 3."},
        {"source": "GSE274960", "model": "LL/2", "kras": "KrasG12C", "icb": "anti-PD-1 ± entrectinib", "scored": "modules", "note": "No R/NR."},
        {"source": "GSE262305", "model": "LLC", "kras": "KrasG12C/Nras", "icb": "anti-PD-L1 ± bortezomib", "scored": "modules", "note": "Cldn4 FPKM is near the floor. No R/NR."},
        {"source": "GSE239485", "model": "LLC", "kras": "KrasG12C/Nras", "icb": "Poly I:C + anti-PD-1", "scored": "modules", "note": "No anti-PD-1 monotherapy and no R/NR."},
        {"source": "GSE260596", "model": "KP", "kras": "KP lung, series summary", "icb": "anti-PD-1 ± anti-LAIR1", "scored": "modules", "note": "Every sample is on anti-PD-1. n=7."},
        {"source": "GSE193895", "model": "K / KL / KK / KKL", "kras": "KrasG12D GEMM", "icb": "no", "scored": "no", "note": "Deposited counts have no ICB arm."},
        {"source": "GSE277929", "model": "KL", "kras": "Kras/Lkb1", "icb": "no", "scored": "no", "note": "Entinostat and trametinib, not ICB."},
        {"source": "GSE137244", "model": "KL vs KP lines", "kras": "Kras ± Lkb1", "icb": "no", "scored": "no", "note": "Locked cell-line contrast. Not rerun."},
    ]


def main() -> None:
    cohorts = load_cohorts()
    rows: list[dict] = []
    for cohort in cohorts:
        rows.extend(sweep_cohort(cohort))
    add_q(rows)
    winners = {q: pick_winner(rows, q) for q in ("cldn4_response", "cldn4_ifn", "cldn4_nhej")}
    by_id = {c.cohort_id: c for c in cohorts}
    for winner in winners.values():
        if not winner:
            continue
        cohort = by_id[winner["cohort"]]
        scale = winner["transform"]
        if scale == "rank_within_sample":
            scale = "log2_cpm" if "log2_cpm" in cohort.matrices else next(iter(cohort.matrices))
        values = cohort.matrices[scale][cohort.gene_index["Cldn4"]]
        winner["cldn4_min"] = float(np.nanmin(values))
        winner["cldn4_max"] = float(np.nanmax(values))
        winner["cldn4_scale"] = scale
        if winner.get("module"):
            gene_rows = module_indices(cohort, winner["module"])
            cols = np.arange(len(cohort.samples))
            scored = score_module(cohort.matrices[scale], gene_rows, cols, winner["score"])
            winner["score_width"] = float(np.nanmax(scored) - np.nanmin(scored))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    samples = sample_table(cohorts)
    write_tsv(OUT / "tables" / "sweep.tsv", rows)
    write_tsv(OUT / "tables" / "samples.tsv", samples)
    write_tsv(OUT / "tables" / "inventory.tsv", inventory_rows())
    prespec_rows = [r for r in rows if r["pre_specified"]]
    write_tsv(OUT / "tables" / "prespecified.tsv", prespec_rows)
    text = write_results_md(cohorts, rows, winners)
    (OUT / "RESULTS.md").write_text(text)
    (ROOT / "RESULTS.md").write_text(text)
    summary = {
        "n_tests": len(rows),
        "winners": winners,
        "n_cohorts": len(cohorts),
        "models": sorted({c.model for c in cohorts}),
    }
    # winners contain numpy types? they are plain dicts of python scalars
    def _conv(obj):
        if isinstance(obj, float) and not math.isfinite(obj):
            return None
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        raise TypeError(type(obj))

    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=_conv))
    plot(cohorts, rows, winners, samples)
    print(text)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
