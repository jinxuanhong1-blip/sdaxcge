#!/usr/bin/env python3
"""Public GEO inventory: STING agonist or STING knockout RNA-seq in lung cancer lines.

Question: does CLDN4 change after STING-pathway activation, or after STING knockout?

Open contrasts computed here
- GSE166209 Calu-3 diABZI vs DMSO (author DESeq2; direct STING agonist)
- GSE166209 mouse lung diABZI vs PBS (author DESeq2; not a cancer line)
- GSE271679 NCI-H596 CRISPR STING KO vs WT (PyDESeq2 on featureCounts)
- GSE244945 SCLC STING KO, with or without NOTCH activation (Welch on author log2 RPKM)
- GSE288796 NSCLC lines, herring-testis DNA vs untreated (cGAS ligand, n=2)
- GSE252340 NCI-H1944 MPS1 inhibitor BAY1217389 (adjacent; not a STING agonist)

GSE134129 is in vivo LLC NanoString after cGAMP. The panel has no Cldn4, so no CLDN4 test.

Private 8-KL matrices are not used. Locked KL-vs-KP, CosMx, and concordant-4 results are not reopened.
"""

from __future__ import annotations

import gzip
import json
import math
import tarfile
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "geo_cache"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

PANEL_HUMAN = {
    "CLDN4": "ENSG00000189143",
    "CLDN1": "ENSG00000163347",
    "CLDN3": "ENSG00000165215",
    "CLDN7": "ENSG00000181885",
    "TACSTD2": "ENSG00000184292",
    "TJP1": "ENSG00000104067",
    "OCLN": "ENSG00000197822",
    "KRT8": "ENSG00000170421",
    "EPCAM": "ENSG00000119888",
    "STING1": "ENSG00000184584",
    "CGAS": "ENSG00000164430",
    "CXCL10": "ENSG00000169245",
    "IFNB1": "ENSG00000171855",
    "IFNL1": "ENSG00000182393",
    "IFIT1": "ENSG00000185745",
    "IFIT2": "ENSG00000119922",
    "ISG15": "ENSG00000187608",
    "MX1": "ENSG00000157601",
    "STAT1": "ENSG00000115415",
    "CCL5": "ENSG00000161570",
    "OAS1": "ENSG00000089127",
}


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"download {url}")
    urllib.request.urlretrieve(url, dest)
    return dest


def geo_suppl(accession: str, filename: str) -> Path:
    nn = accession[:-3] + "nnn"
    url = f"{FTP}/{nn}/{accession}/suppl/{filename}"
    return fetch(url, CACHE / filename)


def welch_log(a: np.ndarray, b: np.ndarray) -> dict:
    """Difference of means on an already log2 scale (a minus b)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    delta = float(a.mean() - b.mean())
    if np.allclose(a, b):
        p = 1.0
        se = 0.0
    else:
        res = stats.ttest_ind(a, b, equal_var=False)
        p = float(res.pvalue)
        se = float(np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)))
    return {
        "log2fc": delta,
        "lfcse": se,
        "pvalue": p,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_a": float(a.mean()),
        "mean_b": float(b.mean()),
    }


def ols_interaction(y: np.ndarray, dox: np.ndarray, ko: np.ndarray) -> dict:
    """y ~ 1 + dox + ko + dox:ko. Returns the interaction term (extra dox effect when STING is intact is -interaction if ko is coded 1).

    Coding: dox=1 treated, ko=1 STING knockout.
    Interaction = (dox effect in KO) - (dox effect in WT)? 
    beta_int = E[y|dox=1,ko=1] - E[y|dox=0,ko=1] - E[y|dox=1,ko=0] + E[y|dox=0,ko=0]
             = dox_effect_KO - dox_effect_WT.
    STING-dependent part of the dox effect is -beta_int (how much larger the dox effect is in WT).
    """
    y = np.asarray(y, dtype=float)
    dox = np.asarray(dox, dtype=float)
    ko = np.asarray(ko, dtype=float)
    X = np.column_stack([np.ones(len(y)), dox, ko, dox * ko])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    df = len(y) - X.shape[1]
    sigma2 = float((resid ** 2).sum() / df)
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = float(math.sqrt(cov[3, 3]))
    tstat = float(beta[3] / se) if se > 0 else 0.0
    p = float(2 * stats.t.sf(abs(tstat), df))
    return {
        "interaction_ko_minus_wt": float(beta[3]),
        "sting_dependent_dox_effect": float(-beta[3]),
        "se": se,
        "pvalue": p,
        "df": int(df),
    }


def load_h596_counts() -> tuple[pd.DataFrame, list[str]]:
    tar_path = geo_suppl("GSE271679", "GSE271679_RAW.tar")
    out = CACHE / "h596"
    out.mkdir(parents=True, exist_ok=True)
    wanted = []
    with tarfile.open(tar_path) as tar:
        for member in tar.getmembers():
            if "H596" in member.name and member.name.endswith(".gz"):
                dest = out / Path(member.name).name
                wanted.append(dest.name)
                if not dest.exists():
                    src = tar.extractfile(member)
                    dest.write_bytes(src.read())
    files = sorted(out.glob("*H596*fcCounts.txt.gz"))
    if len(files) != 8:
        raise RuntimeError(f"expected 8 H596 count files, found {len(files)}")
    genes = None
    columns = {}
    for path in files:
        key = ("WT" if "_WT_" in path.name else "KO") + path.name.split("rep")[1][0]
        ids = []
        counts = []
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if not line.startswith("ENSG"):
                    continue
                parts = line.rstrip("\n").split("\t")
                ids.append(parts[0])
                counts.append(int(parts[6]))
        if genes is None:
            genes = ids
        elif genes != ids:
            raise RuntimeError("H596 featureCounts gene order differs across replicates")
        columns[key] = counts
    order = [f"WT{i}" for i in range(1, 5)] + [f"KO{i}" for i in range(1, 5)]
    mat = pd.DataFrame({s: columns[s] for s in order}, index=genes)
    return mat, order


def run_deseq(counts_gene_by_sample: pd.DataFrame, groups: list[str], ref: str, alt: str) -> tuple[pd.DataFrame, pd.Series]:
    """counts: genes x samples. groups aligned to columns. Contrast is alt vs ref."""
    mat = counts_gene_by_sample.loc[counts_gene_by_sample.sum(axis=1) > 0].astype(int)
    counts = mat.T
    meta = pd.DataFrame({"condition": groups}, index=counts.index)
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition", quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", alt, ref], quiet=True)
    stat.summary()
    size = pd.Series(np.asarray(dds.obs["size_factors"]), index=counts.index)
    return stat.results_df, size


def h596_analysis() -> tuple[list[dict], pd.DataFrame]:
    mat, order = load_h596_counts()
    groups = ["WT"] * 4 + ["KO"] * 4
    res, size = run_deseq(mat, groups, ref="WT", alt="KO")
    norm = mat.div(size, axis=1)
    rows = []
    count_rows = []
    for symbol, ens in PANEL_HUMAN.items():
        if ens not in res.index:
            continue
        r = res.loc[ens]
        rec = {
            "accession": "GSE271679",
            "system": "NCI-H596 STING KO vs WT, unstimulated",
            "class": "STING knockout",
            "gene": symbol,
            "log2fc": float(r["log2FoldChange"]),
            "lfcse": float(r["lfcSE"]),
            "pvalue": float(r["pvalue"]) if pd.notna(r["pvalue"]) else None,
            "padj": float(r["padj"]) if pd.notna(r["padj"]) else None,
            "baseMean": float(r["baseMean"]),
            "n_treat": 4,
            "n_ctrl": 4,
            "method": "PyDESeq2 Wald, BH padj; KO versus WT",
        }
        rows.append(rec)
        count_rows.append(
            {
                "gene": symbol,
                "ensembl": ens,
                **{s: float(norm.loc[ens, s]) for s in order},
            }
        )
    n_sig = int((res["padj"] < 0.05).sum())
    n_tested = int(res["padj"].notna().sum())
    rows.append(
        {
            "accession": "GSE271679",
            "system": "NCI-H596 genome-wide",
            "class": "STING knockout",
            "gene": "_n_padj_lt_0.05",
            "log2fc": None,
            "lfcse": None,
            "pvalue": None,
            "padj": None,
            "baseMean": None,
            "n_treat": 4,
            "n_ctrl": 4,
            "method": f"{n_sig} genes with padj<0.05 out of {n_tested} tested",
        }
    )
    return rows, pd.DataFrame(count_rows)


def sclc_analysis() -> list[dict]:
    path = geo_suppl("GSE244945", "GSE244945_Human_SCLC_RSEM_RPKM_log2.csv.gz")
    df = pd.read_csv(path)
    df = df.set_index(df.columns[0])
    groups = {
        "H82_N1_off": [
            "4_H82_NOTCH1_no_dox_R1_001.fastq.gz",
            "5_H82_NOTCH1_no_dox",
            "6_H82_NOTCH1_no_dox",
        ],
        "H82_N1_on": ["7_H82_NOTCH1_dox", "8_H82_NOTCH1_dox", "9_H82_NOTCH1_dox"],
        "COR_DMSO": ["10_CORL88_DMSO", "11_CORL88_DMSO", "12_CORL88_DMSO"],
        "COR_TAS": ["13_CORL88_TAS1440", "14_CORL88_TAS1440", "15_CORL88_TAS1440"],
        "COR_KO_DMSO": [
            "40_CORL88_STING_KO_DMSO",
            "41_CORL88_STING_KO_DMSO",
            "42_CORL88_STING_KO_DMSO",
        ],
        "COR_KO_TAS": [
            "43_CORL88_STING_KO_TAS1440",
            "44_CORL88_STING_KO_TAS1440",
            "45_CORL88_STING_KO_TAS1440",
        ],
        "H82_SKO_off": [
            "55_H82_NOTCH1_STING_KO_no_dox",
            "56_H82_NOTCH1_STING_KO_no_dox",
            "57_H82_NOTCH1_STING_KO_no_dox",
        ],
        "H82_SKO_on": [
            "58_H82_NOTCH1_STING_KO_dox",
            "59_H82_NOTCH1_STING_KO_dox",
            "60_H82_NOTCH1_STING_KO_dox",
        ],
    }
    for cols in groups.values():
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise KeyError(missing)
    contrasts = [
        ("CORL88 STING KO vs WT, DMSO", "COR_KO_DMSO", "COR_DMSO", "STING knockout"),
        ("CORL88 TAS1440 vs DMSO, STING WT", "COR_TAS", "COR_DMSO", "NOTCH activator, not STING agonist"),
        ("CORL88 TAS1440 vs DMSO, STING KO", "COR_KO_TAS", "COR_KO_DMSO", "NOTCH activator in STING KO"),
        ("H82 STING KO vs STING WT, NOTCH off", "H82_SKO_off", "H82_N1_off", "STING knockout"),
        ("H82 STING KO vs STING WT, NOTCH on", "H82_SKO_on", "H82_N1_on", "STING knockout"),
        ("H82 NOTCH dox vs no dox, STING WT", "H82_N1_on", "H82_N1_off", "NOTCH activation"),
        ("H82 NOTCH dox vs no dox, STING KO", "H82_SKO_on", "H82_SKO_off", "NOTCH activation in STING KO"),
    ]
    genes = ["CLDN4", "CLDN7", "CLDN1", "TACSTD2", "TMEM173", "CXCL10", "IFIT1", "ISG15", "KRT8"]
    rows = []
    for gene in genes:
        if gene not in df.index:
            # matrix uses TMEM173, not STING1
            continue
        series = df.loc[gene]
        if isinstance(series, pd.DataFrame):
            series = series.iloc[0]
        for label, a, b, klass in contrasts:
            stats_ab = welch_log(series[groups[a]].to_numpy(float), series[groups[b]].to_numpy(float))
            rows.append(
                {
                    "accession": "GSE244945",
                    "system": label,
                    "class": klass,
                    "gene": gene,
                    "log2fc": stats_ab["log2fc"],
                    "lfcse": stats_ab["lfcse"],
                    "pvalue": stats_ab["pvalue"],
                    "padj": None,
                    "baseMean": None,
                    "n_treat": 3,
                    "n_ctrl": 3,
                    "method": "Welch t-test on author log2 RPKM; delta is mean(log2 treat) - mean(log2 ctrl). Not genome-wide FDR.",
                    "mean_treat": stats_ab["mean_a"],
                    "mean_ctrl": stats_ab["mean_b"],
                }
            )
        # Factorial interaction on the 12 samples that have both factors.
        for design_name, off_wt, on_wt, off_ko, on_ko in [
            ("H82 NOTCH dox x STING KO", "H82_N1_off", "H82_N1_on", "H82_SKO_off", "H82_SKO_on"),
            ("CORL88 TAS1440 x STING KO", "COR_DMSO", "COR_TAS", "COR_KO_DMSO", "COR_KO_TAS"),
        ]:
            pieces = []
            dox = []
            ko = []
            for key, d, k in (
                (off_wt, 0, 0),
                (on_wt, 1, 0),
                (off_ko, 0, 1),
                (on_ko, 1, 1),
            ):
                vals = series[groups[key]].to_numpy(float)
                pieces.append(vals)
                dox.extend([d] * len(vals))
                ko.extend([k] * len(vals))
            y = np.concatenate(pieces)
            inter = ols_interaction(y, np.array(dox), np.array(ko))
            rows.append(
                {
                    "accession": "GSE244945",
                    "system": design_name + " interaction",
                    "class": "STING-dependent component of NOTCH activation",
                    "gene": gene,
                    "log2fc": inter["sting_dependent_dox_effect"],
                    "lfcse": inter["se"],
                    "pvalue": inter["pvalue"],
                    "padj": None,
                    "baseMean": None,
                    "n_treat": 12,
                    "n_ctrl": 12,
                    "method": "OLS y ~ dox + STINGKO + dox:STINGKO on author log2 RPKM. log2fc column is the STING-dependent dox effect (dox effect in WT minus dox effect in KO).",
                }
            )
    return rows


def htdna_analysis() -> list[dict]:
    path = geo_suppl("GSE288796", "GSE288796_Processed_data_TPM.xlsx")
    book = pd.ExcelFile(path)
    sheet1 = book.parse("Sheet1").set_index("Name")
    sheet2 = book.parse("Sheet2").set_index("Name")
    # Column numbers are the prefixes on the GEO sample titles.
    colmap = {
        1: ("HCC827", "UT"),
        2: ("HCC827", "UT"),
        3: ("HCC827", "DNA"),
        4: ("HCC827", "DNA"),
        5: ("H2228", "UT"),
        6: ("H2228", "UT"),
        7: ("H2228", "DNA"),
        8: ("H2228", "DNA"),
        9: ("H1650", "UT"),
        10: ("H1650", "UT"),
        11: ("H1650", "DNA"),
        12: ("H1650", "DNA"),
        13: ("H596", "UT"),
        14: ("H596", "UT"),
        15: ("H596", "DNA"),
        16: ("H596", "DNA"),
        17: ("H1975", "UT"),
        18: ("H1975", "UT"),
        19: ("H1975", "DNA"),
        20: ("H1975", "DNA"),
    }
    genes = ["CLDN4", "CLDN1", "CLDN7", "TACSTD2", "TJP1", "OCLN", "KRT8", "CXCL10", "IFNB1", "IFNL1", "IFIT1", "ISG15", "CCL5", "TMEM173"]
    rows = []

    def one_line(gene: str, line: str, ut: np.ndarray, dna: np.ndarray, sheet: str) -> dict:
        la = np.log2(dna + 1)
        lb = np.log2(ut + 1)
        st = welch_log(la, lb)
        return {
            "accession": "GSE288796",
            "system": f"{line} HT-DNA vs untreated",
            "class": "cGAS ligand (herring-testis DNA), not a CDN STING agonist",
            "gene": gene,
            "log2fc": st["log2fc"],
            "lfcse": st["lfcse"],
            "pvalue": st["pvalue"],
            "padj": None,
            "baseMean": float(np.mean(np.concatenate([ut, dna]))),
            "n_treat": 2,
            "n_ctrl": 2,
            "method": f"Difference of mean log2(TPM+1). n=2. Sheet {sheet}. Column map is the GEO sample-title prefix. Welch p is descriptive.",
            "tpm_ctrl": ";".join(f"{v:.4f}" for v in ut),
            "tpm_treat": ";".join(f"{v:.4f}" for v in dna),
        }

    for gene in genes:
        if gene not in sheet1.index:
            continue
        row = sheet1.loc[gene]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        buckets: dict[tuple[str, str], list[float]] = {}
        for col, (line, cond) in colmap.items():
            buckets.setdefault((line, cond), []).append(float(row[col]))
        for line in ["HCC827", "H2228", "H1650", "H596", "H1975"]:
            rows.append(
                one_line(
                    gene,
                    line,
                    np.array(buckets[(line, "UT")]),
                    np.array(buckets[(line, "DNA")]),
                    "Sheet1",
                )
            )
        if gene in sheet2.index:
            r2 = sheet2.loc[gene]
            if isinstance(r2, pd.DataFrame):
                r2 = r2.iloc[0]
            # S_1/S_2 untreated, S_3/S_4 HT-DNA, from titles S_1_H358_UT_1 ... S_4_H358_DNA_2
            ut = np.array([float(r2["S_1"]), float(r2["S_2"])])
            dna = np.array([float(r2["S_3"]), float(r2["S_4"])])
            rows.append(one_line(gene, "H358", ut, dna, "Sheet2"))
    return rows


def h1944_analysis() -> list[dict]:
    tar_path = geo_suppl("GSE252340", "GSE252340_RAW.tar")
    out = CACHE / "h1944"
    out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tar:
        for member in tar.getmembers():
            if member.name.endswith(".xlsx"):
                dest = out / Path(member.name).name
                if not dest.exists():
                    src = tar.extractfile(member)
                    dest.write_bytes(src.read())
    frames = []
    for path in sorted(out.glob("GSM*_H1944_*.xlsx")):
        df = pd.read_excel(path, sheet_name=0)
        count_col = [c for c in df.columns if "Read_Count" in str(c)][0]
        parts = path.stem.split("_")
        # GSM7999701_H1944_DMSO_1
        label = "_".join(parts[2:4])
        sub = df[["Gene_Symbol", count_col]].copy()
        sub["Gene_Symbol"] = sub["Gene_Symbol"].astype(str)
        sub = sub.groupby("Gene_Symbol", as_index=False)[count_col].sum()
        sub = sub.rename(columns={count_col: label})
        frames.append(sub)
    merged = frames[0]
    for frame in frames[1:]:
        merged = merged.merge(frame, on="Gene_Symbol", how="outer")
    merged = merged.fillna(0).set_index("Gene_Symbol")
    dmso = [c for c in merged.columns if c.startswith("DMSO")]
    bay = [c for c in merged.columns if c.startswith("BAY")]
    mat = merged[dmso + bay]
    groups = ["DMSO"] * len(dmso) + ["BAY"] * len(bay)
    res, size = run_deseq(mat, groups, ref="DMSO", alt="BAY")
    rows = []
    for gene in ["CLDN4", "CLDN1", "CLDN7", "TACSTD2", "KRT8", "STING1", "CXCL10", "IFNB1", "IFNL1", "IFIT1", "ISG15", "MX1", "STAT1", "CCL5"]:
        if gene not in res.index:
            continue
        r = res.loc[gene]
        rows.append(
            {
                "accession": "GSE252340",
                "system": "NCI-H1944 BAY1217389 vs DMSO, 48 h",
                "class": "adjacent: MPS1 inhibitor, micronucleus cGAS-STING activation, not a STING agonist or STING KO",
                "gene": gene,
                "log2fc": float(r["log2FoldChange"]),
                "lfcse": float(r["lfcSE"]),
                "pvalue": float(r["pvalue"]) if pd.notna(r["pvalue"]) else None,
                "padj": float(r["padj"]) if pd.notna(r["padj"]) else None,
                "baseMean": float(r["baseMean"]),
                "n_treat": 2,
                "n_ctrl": 2,
                "method": "PyDESeq2 on author read counts, first sheet of each GEO sample only. n=2, so the dispersion prior carries weight. Duplicate symbols were summed.",
            }
        )
    return rows


def author_contrasts() -> list[dict]:
    specs = [
        (
            "GSE166209_Compare_Treatment_diABZI_6_vs_DMSO.csv.gz",
            "Calu-3 diABZI vs DMSO, 6 h",
            "direct STING agonist",
            "GSE166209",
            2,
            ["CLDN4", "CLDN1", "CLDN3", "CLDN7", "TACSTD2", "TJP1", "OCLN", "KRT8", "EPCAM", "CXCL10", "IFIT1", "IFIT2", "ISG15", "MX1", "IFNB1", "IFNL1", "TMEM173", "OAS1"],
        ),
        (
            "GSE166209_Compare_Treatment_diABZI_12_vs_DMSO.csv.gz",
            "Calu-3 diABZI vs DMSO, 12 h",
            "direct STING agonist",
            "GSE166209",
            2,
            ["CLDN4", "CLDN1", "CLDN3", "CLDN7", "TACSTD2", "TJP1", "OCLN", "KRT8", "EPCAM", "CXCL10", "IFIT1", "IFIT2", "ISG15", "MX1", "IFNB1", "IFNL1", "TMEM173", "OAS1"],
        ),
        (
            "GSE166209_Compare_condition_diABZI.06_vs_mock.06.csv.gz",
            "Mouse lung diABZI vs PBS, 6 h",
            "same series, not a cancer line",
            "GSE166209",
            5,
            ["Cldn4", "Cldn1", "Cldn3", "Cldn7", "Tacstd2", "Tjp1", "Cxcl10", "Ifit1", "Isg15", "Sting1", "Krt8"],
        ),
        (
            "GSE166209_Compare_condition_diABZI.12_mock.12.csv.gz",
            "Mouse lung diABZI vs PBS, 12 h",
            "same series, not a cancer line",
            "GSE166209",
            5,
            ["Cldn4", "Cldn1", "Cldn3", "Cldn7", "Tacstd2", "Tjp1", "Cxcl10", "Ifit1", "Isg15", "Sting1", "Krt8"],
        ),
    ]
    rows = []
    for filename, system, klass, acc, n, genes in specs:
        path = geo_suppl(acc, filename)
        df = pd.read_csv(path)
        name_col = "name"
        for gene in genes:
            hit = df.loc[df[name_col].astype(str) == gene]
            if hit.empty:
                rows.append(
                    {
                        "accession": acc,
                        "system": system,
                        "class": klass,
                        "gene": gene,
                        "log2fc": None,
                        "lfcse": None,
                        "pvalue": None,
                        "padj": None,
                        "baseMean": None,
                        "n_treat": n,
                        "n_ctrl": n,
                        "method": "author DESeq2 supplementary table; gene absent",
                    }
                )
                continue
            r = hit.iloc[0]
            rows.append(
                {
                    "accession": acc,
                    "system": system,
                    "class": klass,
                    "gene": gene,
                    "log2fc": float(r["log2FoldChange"]),
                    "lfcse": float(r["lfcSE"]),
                    "pvalue": float(r["pvalue"]),
                    "padj": float(r["padj"]) if pd.notna(r["padj"]) else None,
                    "baseMean": float(r["baseMean"]),
                    "n_treat": n,
                    "n_ctrl": n,
                    "method": "author DESeq2 supplementary contrast; not recomputed",
                }
            )
    return rows


def llc_note() -> dict:
    path = geo_suppl("GSE134129", "GSE134129_Normalized_matrix.xlsx")
    df = pd.read_excel(path)
    symbols = set(df["ID_REF"].astype(str).str.lower())
    has_cldn4 = "cldn4" in symbols
    def means(gene: str) -> dict | None:
        hit = df.loc[df["ID_REF"].astype(str).str.lower() == gene.lower()]
        if hit.empty:
            return None
        r = hit.iloc[0]
        ctrl = [c for c in df.columns if str(c).startswith("Control")]
        sting = [c for c in df.columns if str(c).startswith("STING_") or str(c).startswith("STING ") ]
        # columns are STING_1 ... and 'STING KO_1'
        sting = [c for c in df.columns if str(c).startswith("STING_")]
        ko = [c for c in df.columns if "KO" in str(c)]
        def lv(cols):
            vals = r[cols].to_numpy(float)
            return float(np.log2(vals + 1).mean()), vals
        c_m, _ = lv(ctrl)
        s_m, _ = lv(sting)
        k_m, _ = lv(ko)
        return {"control_mean_log2": c_m, "sting_mean_log2": s_m, "ko_mean_log2": k_m, "sting_minus_control": s_m - c_m}
    qc = {g: means(g) for g in ["Cxcl10", "Ifnb1", "Ifit1", "Tmem173"]}
    return {"cldn4_on_panel": has_cldn4, "qc_log2": qc, "n_genes": int(df.shape[0])}


def inventory_rows(llc: dict) -> list[dict]:
    return [
        {
            "accession": "GSE166209",
            "match": "direct STING agonist",
            "system": "Calu-3 lung adenocarcinoma line, diABZI 6 h and 12 h vs DMSO, n=2",
            "open_de": "author DESeq2 CSV",
            "cldn4": "tested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166209",
        },
        {
            "accession": "GSE166209",
            "match": "same series, not a cancer line",
            "system": "K18-hACE2 mouse lung tissue, diABZI vs PBS, 6 h and 12 h, n=5",
            "open_de": "author DESeq2 CSV",
            "cldn4": "tested, reported separately",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166209",
        },
        {
            "accession": "GSE271679",
            "match": "STING knockout, lung cancer line",
            "system": "NCI-H596 WT vs CRISPR STING KO, unstimulated, n=4. Same series also has SCC25, OE21, Detroit 562 (not lung) and Calu-3 WT only.",
            "open_de": "featureCounts in RAW.tar; PyDESeq2 run here for H596 only",
            "cldn4": "tested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271679",
        },
        {
            "accession": "GSE244945",
            "match": "STING knockout, SCLC lines",
            "system": "NCI-H82 and CORL88 STING KO, crossed with NOTCH activation (dox or TAS1440). TAS1440 is a NOTCH activator, not a STING agonist.",
            "open_de": "author log2 RPKM matrix",
            "cldn4": "tested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244945",
        },
        {
            "accession": "GSE288796",
            "match": "cGAS ligand, not a CDN STING agonist",
            "system": "HT-DNA vs untreated in HCC827, H2228, H1650, H596, H1975, H358. n=2.",
            "open_de": "author TPM workbook",
            "cldn4": "tested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288796",
        },
        {
            "accession": "GSE252340",
            "match": "adjacent, not agonist or knockout",
            "system": "NCI-H1944 (KL) BAY1217389 MPS1 inhibitor 48 h vs DMSO, n=2. Authors use it to drive micronuclei and cGAS-STING.",
            "open_de": "per-sample read-count xlsx; PyDESeq2 run here",
            "cldn4": "tested, labeled adjacent",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252340",
        },
        {
            "accession": "GSE134129",
            "match": "in vivo LLC tumors, cGAMP, not cultured-line RNA-seq",
            "system": "NanoString PanCancer panel, 750 genes, cGAMP-treated LLC tumors in WT mice vs tumors from STING KO mice. Cldn4 is not on the panel.",
            "open_de": "normalized matrix; CLDN4 test not possible",
            "cldn4": "absent",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134129",
            "note": json.dumps({"cldn4_on_panel": llc["cldn4_on_panel"], "n_genes": llc["n_genes"]}),
        },
        {
            "accession": "GSE254174",
            "match": "excluded, not lung cancer",
            "system": "U2OS and BJ fibroblasts, cGAMP or DMXAA, including STING KO",
            "open_de": "TPM available; not used",
            "cldn4": "not tested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254174",
        },
        {
            "accession": "GSE269551",
            "match": "excluded, not a STING perturbation",
            "system": "H2122 parental vs SMARCA4 KO. STING is in the paper's mechanism, not the contrast.",
            "open_de": "htseq counts; not used",
            "cldn4": "not tested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269551",
        },
        {
            "accession": "GSE305239",
            "match": "excluded, not lung",
            "system": "H4 glioblastoma, diABZI",
            "open_de": "not downloaded",
            "cldn4": "not tested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305239",
        },
        {
            "accession": "GSE308610",
            "match": "excluded, not lung",
            "system": "THP-1 and RPMI-8226, ADU-S100",
            "open_de": "not downloaded",
            "cldn4": "not tested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308610",
        },
        {
            "accession": "GSE147085",
            "match": "excluded, not lung",
            "system": "FaDu HNSCC STING KO",
            "open_de": "not downloaded",
            "cldn4": "not tested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147085",
        },
        {
            "accession": "GSE137244",
            "match": "excluded, different contrast already reported",
            "system": "Mouse KL vs KP cell lines. Baseline genotype, not STING agonist or STING KO.",
            "open_de": "not reanalyzed",
            "cldn4": "not retested",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137244",
        },
    ]


def write_figure(contrasts: pd.DataFrame) -> None:
    order = [
        ("GSE166209", "Calu-3 diABZI vs DMSO, 6 h", "Calu-3 diABZI 6 h", "#0072B2"),
        ("GSE166209", "Calu-3 diABZI vs DMSO, 12 h", "Calu-3 diABZI 12 h", "#0072B2"),
        ("GSE271679", "NCI-H596 STING KO vs WT, unstimulated", "H596 STING KO", "#D55E00"),
        ("GSE288796", "HCC827 HT-DNA vs untreated", "HCC827 HT-DNA", "#009E73"),
        ("GSE288796", "H2228 HT-DNA vs untreated", "H2228 HT-DNA", "#009E73"),
        ("GSE288796", "H1650 HT-DNA vs untreated", "H1650 HT-DNA", "#009E73"),
        ("GSE288796", "H596 HT-DNA vs untreated", "H596 HT-DNA", "#009E73"),
        ("GSE288796", "H1975 HT-DNA vs untreated", "H1975 HT-DNA", "#009E73"),
        ("GSE288796", "H358 HT-DNA vs untreated", "H358 HT-DNA", "#009E73"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 5.2), sharey=True)
    for ax, gene, title in zip(axes, ["CLDN4", "CXCL10"], ["CLDN4", "CXCL10 (pathway control)"]):
        ys = []
        for i, (acc, system, label, color) in enumerate(order):
            hit = contrasts[
                (contrasts.accession == acc)
                & (contrasts.system == system)
                & (contrasts.gene == gene)
            ]
            if hit.empty:
                continue
            r = hit.iloc[0]
            y = len(order) - 1 - i
            ys.append(y)
            x = float(r["log2fc"])
            se = float(r["lfcse"]) if pd.notna(r["lfcse"]) else 0.0
            ax.errorbar(x, y, xerr=se, fmt="o", color=color, ms=7, lw=1.2, capsize=2)
        ax.axvline(0, color="#444444", lw=0.8)
        ax.set_yticks(list(range(len(order))))
        ax.set_yticklabels([lab for *_, lab, _ in order][::-1])
        ax.set_xlabel("log2 fold change")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("Lung cancer lines: STING agonist, STING knockout, or HT-DNA", fontsize=11)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "cldn4_sting_lung_log2fc.png", dpi=160)
    fig.savefig(FIGS / "cldn4_sting_lung_log2fc.pdf")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    author = author_contrasts()
    h596_rows, h596_counts = h596_analysis()
    sclc = sclc_analysis()
    htdna = htdna_analysis()
    h1944 = h1944_analysis()
    llc = llc_note()
    inv = inventory_rows(llc)

    contrasts = pd.DataFrame(author + h596_rows + sclc + htdna + h1944)
    contrasts.to_csv(TABLES / "contrasts.tsv", sep="\t", index=False)
    h596_counts.to_csv(TABLES / "h596_stingko_normalized_counts.tsv", sep="\t", index=False)
    pd.DataFrame(inv).to_csv(TABLES / "inventory.tsv", sep="\t", index=False)
    (TABLES / "gse134129_panel_qc.json").write_text(json.dumps(llc, indent=2))
    write_figure(contrasts)

    def pick(acc, system, gene):
        hit = contrasts[(contrasts.accession == acc) & (contrasts.system == system) & (contrasts.gene == gene)]
        if hit.empty:
            return None
        r = hit.iloc[0]
        return {
            "log2fc": None if pd.isna(r.log2fc) else round(float(r.log2fc), 4),
            "lfcse": None if pd.isna(r.lfcse) else round(float(r.lfcse), 4),
            "pvalue": None if pd.isna(r.pvalue) else float(r.pvalue),
            "padj": None if pd.isna(r.padj) else float(r.padj),
            "baseMean": None if pd.isna(r.baseMean) else round(float(r.baseMean), 2),
        }

    highlight = {
        "calu3_6h_CLDN4": pick("GSE166209", "Calu-3 diABZI vs DMSO, 6 h", "CLDN4"),
        "calu3_12h_CLDN4": pick("GSE166209", "Calu-3 diABZI vs DMSO, 12 h", "CLDN4"),
        "calu3_6h_CXCL10": pick("GSE166209", "Calu-3 diABZI vs DMSO, 6 h", "CXCL10"),
        "calu3_12h_CXCL10": pick("GSE166209", "Calu-3 diABZI vs DMSO, 12 h", "CXCL10"),
        "calu3_6h_IFIT1": pick("GSE166209", "Calu-3 diABZI vs DMSO, 6 h", "IFIT1"),
        "calu3_12h_IFIT1": pick("GSE166209", "Calu-3 diABZI vs DMSO, 12 h", "IFIT1"),
        "calu3_12h_KRT8": pick("GSE166209", "Calu-3 diABZI vs DMSO, 12 h", "KRT8"),
        "calu3_12h_TACSTD2": pick("GSE166209", "Calu-3 diABZI vs DMSO, 12 h", "TACSTD2"),
        "mouse_12h_Cldn4": pick("GSE166209", "Mouse lung diABZI vs PBS, 12 h", "Cldn4"),
        "mouse_6h_Cldn4": pick("GSE166209", "Mouse lung diABZI vs PBS, 6 h", "Cldn4"),
        "h596_CLDN4": pick("GSE271679", "NCI-H596 STING KO vs WT, unstimulated", "CLDN4"),
        "h596_STING1": pick("GSE271679", "NCI-H596 STING KO vs WT, unstimulated", "STING1"),
        "h596_CXCL10": pick("GSE271679", "NCI-H596 STING KO vs WT, unstimulated", "CXCL10"),
        "h596_IFIT1": pick("GSE271679", "NCI-H596 STING KO vs WT, unstimulated", "IFIT1"),
        "h596_TACSTD2": pick("GSE271679", "NCI-H596 STING KO vs WT, unstimulated", "TACSTD2"),
        "h596_CLDN7": pick("GSE271679", "NCI-H596 STING KO vs WT, unstimulated", "CLDN7"),
        "h596_CLDN1": pick("GSE271679", "NCI-H596 STING KO vs WT, unstimulated", "CLDN1"),
        "h596_KRT8": pick("GSE271679", "NCI-H596 STING KO vs WT, unstimulated", "KRT8"),
        "h1944_CLDN4": pick("GSE252340", "NCI-H1944 BAY1217389 vs DMSO, 48 h", "CLDN4"),
        "h1944_IFIT1": pick("GSE252340", "NCI-H1944 BAY1217389 vs DMSO, 48 h", "IFIT1"),
        "h1944_CXCL10": pick("GSE252340", "NCI-H1944 BAY1217389 vs DMSO, 48 h", "CXCL10"),
        "llc": llc,
    }
    for line in ["HCC827", "H2228", "H1650", "H596", "H1975", "H358"]:
        highlight[f"htdna_{line}_CLDN4"] = pick("GSE288796", f"{line} HT-DNA vs untreated", "CLDN4")
        highlight[f"htdna_{line}_CXCL10"] = pick("GSE288796", f"{line} HT-DNA vs untreated", "CXCL10")
    (TABLES / "highlights.json").write_text(json.dumps(highlight, indent=2))
    print(json.dumps(highlight, indent=2))


if __name__ == "__main__":
    main()
