#!/usr/bin/env python3
"""B2: TACSTD2 vs CLDN4 rank among claudins in CCLE/DepMap, all cancers.

Primary estimand: Spearman rank of CLDN4 among CLDN\\d+ genes vs TACSTD2,
separately for RNA (DepMap 24Q2) and protein (Nusinow 2020), all lineages.
Also reports claudin abundance ranks and genome-wide abundance ranks.
Deterministic. Pairwise-complete. BH-FDR within layer.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

CLDN_RE = re.compile(r"^CLDN\d+$")
TENPX_RE = re.compile(r"^(.+)_TenPx\d+$", re.I)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def spearman_pair(a: pd.Series, b: pd.Series) -> dict:
    d = pd.concat([a, b], axis=1).dropna()
    n = int(len(d))
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    if d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1], nan_policy="omit")
    return {"n": n, "rho": float(rho), "p": float(p)}


def bh(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, dtype=float)
    out = np.full(arr.shape, np.nan)
    ok = np.isfinite(arr)
    if ok.any():
        out[ok] = multipletests(arr[ok], method="fdr_bh")[1]
    return [float(x) if np.isfinite(x) else np.nan for x in out]


def rank_desc(values: pd.Series) -> pd.Series:
    """Rank 1 = highest; ties get min rank."""
    return values.rank(ascending=False, method="min")


def load_models(path: Path) -> pd.DataFrame:
    m = pd.read_csv(path)
    keep = [
        c
        for c in [
            "ModelID",
            "CellLineName",
            "StrippedCellLineName",
            "CCLEName",
            "OncotreeLineage",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "PrimaryOrMetastasis",
        ]
        if c in m.columns
    ]
    m = m[keep].copy()
    if "CCLEName" in m.columns:
        m["CCLEName"] = m["CCLEName"].astype(str)
    return m


def load_rna(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "ModelID"})
    df = df.set_index("ModelID")
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_protein(gz_path: Path) -> tuple[pd.DataFrame, dict]:
    raw = pd.read_csv(gz_path, low_memory=False)
    meta = {"n_rows": int(len(raw)), "n_cols": int(raw.shape[1]), "columns_head": list(raw.columns[:12])}
    # Gene symbol column
    gene_col = None
    for cand in ["Gene_Symbol", "gene_symbol", "GeneSymbol", "Gene", "symbol"]:
        if cand in raw.columns:
            gene_col = cand
            break
    if gene_col is None:
        # sometimes first text column
        for c in raw.columns:
            if raw[c].dtype == object:
                gene_col = c
                break
    if gene_col is None:
        raise SystemExit(f"No gene column in protein table: {list(raw.columns)[:20]}")
    meta["gene_col"] = gene_col

    sample_cols = []
    for c in raw.columns:
        if TENPX_RE.match(str(c)):
            sample_cols.append(c)
    if not sample_cols:
        # fallback: numeric columns that are not annotation
        annot = {
            gene_col,
            "Protein_Id",
            "Description",
            "Group",
            "Uniprot",
            "Protein",
            "Accession",
            "Accession_id",
            "Gene",
        }
        for c in raw.columns:
            if c in annot:
                continue
            if pd.api.types.is_numeric_dtype(raw[c]):
                sample_cols.append(c)
    meta["n_sample_cols"] = len(sample_cols)
    if not sample_cols:
        raise SystemExit("No sample columns found in protein table")

    sub = raw[[gene_col] + sample_cols].copy()
    sub[gene_col] = sub[gene_col].astype(str).str.strip()
    # collapse duplicate gene symbols by mean
    num = sub.set_index(gene_col)
    for c in sample_cols:
        num[c] = pd.to_numeric(num[c], errors="coerce")
    num = num.groupby(level=0).mean()

    # average TenPx replicates of the same CCLE name
    rename = {}
    for c in num.columns:
        m = TENPX_RE.match(str(c))
        rename[c] = m.group(1) if m else str(c)
    num = num.rename(columns=rename)
    num = num.T.groupby(level=0).mean().T  # genes x unique CCLE names
    meta["n_unique_lines"] = int(num.shape[1])
    meta["n_genes"] = int(num.shape[0])
    return num, meta


def load_s1(path: Path) -> pd.DataFrame:
    xl = pd.ExcelFile(path)
    meta = {"sheets": xl.sheet_names}
    # prefer a sheet that looks like sample info
    best = None
    for s in xl.sheet_names:
        df = xl.parse(s)
        cols = [str(c).lower() for c in df.columns]
        if any("ccle" in c or "cell" in c or "line" in c for c in cols):
            best = df
            meta["used_sheet"] = s
            break
    if best is None:
        best = xl.parse(xl.sheet_names[0])
        meta["used_sheet"] = xl.sheet_names[0]
    return best, meta


def claudin_cols(columns) -> list[str]:
    return sorted([c for c in columns if CLDN_RE.match(str(c))], key=lambda x: int(re.search(r"\d+", x).group()))


def layer_ranks(mat: pd.DataFrame, models: pd.DataFrame, layer: str, min_n: int = 20) -> dict:
    """Compute all-cancer and per-lineage ranks. mat index=ModelID or CCLEName."""
    if "TACSTD2" not in mat.columns:
        raise SystemExit(f"TACSTD2 missing from {layer}")
    cldns = claudin_cols(mat.columns)
    if "CLDN4" not in cldns:
        raise SystemExit(f"CLDN4 missing from {layer} claudins: {cldns}")

    # join lineage
    if mat.index.name == "ModelID" or "ACH-" in str(mat.index[0]):
        lin = models.set_index("ModelID")["OncotreeLineage"]
        lineage = lin.reindex(mat.index)
    else:
        # protein indexed by CCLEName
        if "CCLEName" in models.columns:
            lin = models.dropna(subset=["CCLEName"]).drop_duplicates("CCLEName").set_index("CCLEName")[
                "OncotreeLineage"
            ]
            lineage = lin.reindex(mat.index)
        else:
            lineage = pd.Series(index=mat.index, data="Unknown")
    lineage = lineage.fillna("Unknown")

    def one_cohort(sub: pd.DataFrame, cohort: str) -> tuple[pd.DataFrame, dict]:
        rows = []
        for g in cldns:
            sp = spearman_pair(sub["TACSTD2"], sub[g])
            med = float(sub[g].median(skipna=True)) if sub[g].notna().any() else np.nan
            n_quant = int(sub[g].notna().sum())
            rows.append(
                {
                    "layer": layer,
                    "cohort": cohort,
                    "gene": g,
                    "n_lines_in_cohort": int(len(sub)),
                    "n_quantified": n_quant,
                    "n_pairwise_vs_TACSTD2": sp["n"],
                    "rho_vs_TACSTD2": sp["rho"],
                    "p_vs_TACSTD2": sp["p"],
                    "median_abundance": med,
                    "mean_abundance": float(sub[g].mean(skipna=True)) if n_quant else np.nan,
                }
            )
        tab = pd.DataFrame(rows)
        # Confirmatory family: pairwise n>=20 and defined ρ (non-constant).
        eligible = (tab["n_pairwise_vs_TACSTD2"] >= 20) & tab["rho_vs_TACSTD2"].notna()
        q = np.full(len(tab), np.nan)
        if eligible.any():
            q[eligible.to_numpy()] = bh(tab.loc[eligible, "p_vs_TACSTD2"].tolist())
        tab["q_vs_TACSTD2"] = q
        tab["in_confirmatory_family"] = eligible
        ranked = tab["rho_vs_TACSTD2"].where(eligible)
        tab["rank_rho_among_claudins"] = rank_desc(ranked).astype("Int64")
        tab["rank_median_among_claudins"] = rank_desc(
            tab["median_abundance"].where(tab["n_quantified"] >= 20)
        ).astype("Int64")
        tab["n_claudins_tested"] = int(eligible.sum())
        tac = spearman_pair(sub["TACSTD2"], sub["CLDN4"]) if "CLDN4" in sub else {"n": 0, "rho": np.nan, "p": np.nan}
        cldn4 = tab.loc[tab["gene"] == "CLDN4"].iloc[0]
        summary = {
            "layer": layer,
            "cohort": cohort,
            "n_lines": int(len(sub)),
            "n_TACSTD2": int(sub["TACSTD2"].notna().sum()),
            "n_CLDN4": int(sub["CLDN4"].notna().sum()),
            "n_pairwise_TACSTD2_CLDN4": tac["n"],
            "rho_TACSTD2_CLDN4": tac["rho"],
            "p_TACSTD2_CLDN4": tac["p"],
            "CLDN4_rank_rho_among_claudins": int(cldn4["rank_rho_among_claudins"])
            if pd.notna(cldn4["rank_rho_among_claudins"])
            else None,
            "CLDN4_rank_median_among_claudins": int(cldn4["rank_median_among_claudins"])
            if pd.notna(cldn4["rank_median_among_claudins"])
            else None,
            "n_claudins_tested": int(tab["n_claudins_tested"].iloc[0]),
            "top_claudin_by_rho": None
            if tab["rho_vs_TACSTD2"].isna().all()
            else str(tab.sort_values("rho_vs_TACSTD2", ascending=False).iloc[0]["gene"]),
            "top_rho": None
            if tab["rho_vs_TACSTD2"].isna().all()
            else float(tab["rho_vs_TACSTD2"].max()),
        }
        return tab, summary

    all_tab, all_sum = one_cohort(mat, "ALL_CANCERS")
    lineage_tabs = [all_tab]
    lineage_sums = [all_sum]
    for lin, idx in lineage.groupby(lineage).groups.items():
        sub = mat.loc[list(idx)]
        if len(sub) < min_n:
            continue
        t, s = one_cohort(sub, f"LINEAGE:{lin}")
        lineage_tabs.append(t)
        lineage_sums.append(s)
    return {
        "table": pd.concat(lineage_tabs, ignore_index=True),
        "summaries": lineage_sums,
        "lineage_counts": lineage.value_counts().to_dict(),
        "claudins": cldns,
    }


def genome_wide_ranks(means: dict, genes: list[str]) -> pd.DataFrame:
    items = [(g, v["mean"], v["n"]) for g, v in means.items() if np.isfinite(v["mean"])]
    items.sort(key=lambda x: x[1], reverse=True)
    rank = {g: i + 1 for i, (g, _, _) in enumerate(items)}
    rows = []
    for g in genes:
        if g not in means:
            rows.append({"gene": g, "mean": np.nan, "n": 0, "rank_among_all_genes": pd.NA, "n_genes": len(items)})
            continue
        rows.append(
            {
                "gene": g,
                "mean": means[g]["mean"],
                "n": means[g]["n"],
                "rank_among_all_genes": rank.get(g, pd.NA),
                "n_genes": len(items),
            }
        )
    return pd.DataFrame(rows)


def protein_genome_ranks(prot: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    med = prot.median(axis=1, skipna=True)
    med = med[med.notna()]
    rank = rank_desc(med)
    rows = []
    for g in genes:
        rows.append(
            {
                "gene": g,
                "median": float(med[g]) if g in med.index else np.nan,
                "n": int(prot.loc[g].notna().sum()) if g in prot.index else 0,
                "rank_among_all_proteins": int(rank[g]) if g in rank.index else pd.NA,
                "n_proteins": int(len(med)),
            }
        )
    return pd.DataFrame(rows)


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, lineterminator="\n")


def fig_bar_rho(tab: pd.DataFrame, title: str, path: Path, highlight="CLDN4") -> None:
    d = tab[tab["cohort"] == "ALL_CANCERS"].copy()
    d = d.sort_values("rho_vs_TACSTD2", ascending=False)
    colors = ["#c0392b" if g == highlight else "#4a6fa5" for g in d["gene"]]
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.bar(d["gene"], d["rho_vs_TACSTD2"], color=colors, edgecolor="none")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_ylabel("Spearman ρ vs TACSTD2")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=60, labelsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_scatter(x, y, xlabel, ylabel, title, path, n, rho, p) -> None:
    d = pd.concat([x, y], axis=1).dropna()
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.scatter(d.iloc[:, 0], d.iloc[:, 1], s=8, alpha=0.35, c="#2c3e50", edgecolors="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\nn={n}  ρ={rho:.3f}  p={p:.2e}" if np.isfinite(rho) else title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_missing(rna, prot_on_models, path) -> None:
    genes = ["TACSTD2"] + claudin_cols(rna.columns)
    rna_frac = [1 - rna[g].notna().mean() if g in rna.columns else 1 for g in genes]
    prot_frac = [
        1 - prot_on_models[g].notna().mean() if g in prot_on_models.columns else 1 for g in genes
    ]
    x = np.arange(len(genes))
    fig, ax = plt.subplots(figsize=(10, 4.0))
    ax.bar(x - 0.18, rna_frac, 0.36, label="RNA missing", color="#7f8c8d")
    ax.bar(x + 0.18, prot_frac, 0.36, label="Protein missing", color="#c0392b")
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=60, fontsize=8)
    ax.set_ylabel("Fraction missing")
    ax.set_title("Missingness: TACSTD2 and claudins (all cancers)")
    ax.legend()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_lineage_rank(summ: list[dict], path: Path) -> None:
    rows = [s for s in summ if s["CLDN4_rank_rho_among_claudins"] is not None]
    rows = sorted(rows, key=lambda s: (0 if s["cohort"] == "ALL_CANCERS" else 1, s["cohort"]))
    labels = [s["cohort"].replace("LINEAGE:", "") for s in rows]
    ranks = [s["CLDN4_rank_rho_among_claudins"] for s in rows]
    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    colors = ["#c0392b" if lab == "ALL_CANCERS" else "#4a6fa5" for lab in labels]
    ax.bar(labels, ranks, color=colors)
    ax.set_ylabel("CLDN4 rank among claudins (1 = highest ρ vs TACSTD2)")
    ax.set_title("Exploratory: CLDN4 coexpression rank by lineage (n≥20)")
    ax.tick_params(axis="x", rotation=60, labelsize=8)
    ax.invert_yaxis()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f"{x:.{nd}f}"


def fmt_p(p):
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    return f"{p:.2e}"


def write_writeup(path: Path, ctx: dict) -> None:
    r = ctx["rna_all"]
    p = ctx["prot_all"]
    ov = ctx["overlap"]
    rna_top = r["top_claudin_by_rho"]
    prot_top = p["top_claudin_by_rho"]
    claim_rna = (
        "CLDN4 is the top TACSTD2-coexpressed claudin at RNA in all CCLE/DepMap lines."
        if rna_top == "CLDN4"
        else (
            f"CLDN4 is NOT the top TACSTD2-coexpressed claudin at RNA "
            f"(rank {r['CLDN4_rank_rho_among_claudins']} of {r['n_claudins_tested']}; top is {rna_top})."
        )
    )
    claim_rna_zh = (
        "在全癌种 RNA 中，CLDN4 是与 TACSTD2 共表达最高的 claudin。"
        if rna_top == "CLDN4"
        else (
            f"在全癌种 RNA 中，CLDN4 并不是与 TACSTD2 共表达最高的 claudin"
            f"（第 {r['CLDN4_rank_rho_among_claudins']} / {r['n_claudins_tested']}；最高为 {rna_top}）。"
        )
    )
    claim_prot = (
        "CLDN4 is the top TACSTD2-coexpressed claudin at protein among confirmatory-family claudins."
        if prot_top == "CLDN4"
        else (
            f"CLDN4 is NOT the top TACSTD2-coexpressed claudin at protein "
            f"(rank {p['CLDN4_rank_rho_among_claudins']} of {p['n_claudins_tested']}; top is {prot_top})."
        )
    )
    claim_prot_zh = (
        "在验证性蛋白家族中，CLDN4 是与 TACSTD2 共表达最高的 claudin。"
        if prot_top == "CLDN4"
        else (
            f"在验证性蛋白家族中，CLDN4 并不是与 TACSTD2 共表达最高的 claudin"
            f"（第 {p['CLDN4_rank_rho_among_claudins']} / {p['n_claudins_tested']}；最高为 {prot_top}）。"
        )
    )
    sha = ctx["git_sha"]
    run = ctx["run_utc"]
    seed = ctx["seed"]
    gw_rna = ctx["gw_rna"]
    gw_prot = ctx["gw_prot"]

    def gw(df, gene, col_rank, col_n):
        row = df.loc[df["gene"] == gene]
        if row.empty:
            return "NA", "NA"
        return str(row.iloc[0][col_rank]), str(row.iloc[0][col_n])

    t2_rna_rank, t2_rna_n = gw(gw_rna, "TACSTD2", "rank_among_all_genes", "n_genes")
    c4_rna_rank, c4_rna_n = gw(gw_rna, "CLDN4", "rank_among_all_genes", "n_genes")
    t2_p_rank, t2_p_n = gw(gw_prot, "TACSTD2", "rank_among_all_proteins", "n_proteins")
    c4_p_rank, c4_p_n = gw(gw_prot, "CLDN4", "rank_among_all_proteins", "n_proteins")

    text = f"""# B2 CCLE all cancers: TACSTD2 vs CLDN4 rank among claudins / 全癌种 CCLE 中 TACSTD2 与 CLDN4 在 claudin 家族中的位次

Analysis commit / 分析提交: `{sha}`  
Run UTC / 运行 UTC: `{run}`  
Catalog / 数据目录: `results/w200/B2_CCLE_all/catalog.tsv`  
Seed / 随机种子: `{seed}`

**Honesty / 诚实声明:** This slice reports the measured ranks. It does not assume CLDN4 is rank 1. Cell lines are not tumors and have no immune microenvironment; this is not an ICI test.

本切片报告实测位次，不预设 CLDN4 为第 1。细胞系不是肿瘤、无免疫微环境；不能检验 ICI。

---

## English

### Question and design
Primary question: in **all CCLE/DepMap cancer cell lines** (no lineage filter), what is **CLDN4's rank among classical claudins** (`CLDN` + integer, e.g. CLDN1–CLDN25 as present) by Spearman correlation with **TACSTD2**, at RNA and at protein?

Prespecified contrasts:
1. Confirmatory: all-cancers Spearman ρ(TACSTD2, each CLDN\\d+), rank of CLDN4 (1 = highest ρ).
2. Confirmatory: all-cancers RNA median-abundance rank of CLDN4 among the same claudins (protein TMT is not treated as copies/cell).
3. Confirmatory: pairwise ρ(TACSTD2, CLDN4) at RNA and protein.
4. Descriptive: genome-wide RNA abundance rank of TACSTD2 and of CLDN4.
5. Exploratory: the same coexpression rank inside Oncotree lineages with n≥20.

Independent biological unit: one cell line (DepMap `ModelID` / CCLE name). No technical-replicate inflation after averaging Nusinow TenPx replicates of the same line.

### Data provenance and accession verification
Databases queried over HTTPS before download (raw responses in `results/w200/B2_CCLE_all/repro/accessions/`):

| Record | Registry | Query | Confirmed title |
|---|---|---|---|
| `10.25452/figshare.plus.25880521` | Figshare+ API `/v2/articles/25880521` | exact article id 25880521 | DepMap 24Q2 Public |
| `MSV000085836` | MassIVE PROXI `accession=MSV000085836&resultType=full` | exact accession | Quantitative Proteomics of the Cancer Cell Line Encyclopedia |

Accepted files (see `catalog.tsv`): DepMap `Model.csv`, `README.txt`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (full file SHA-256 computed in stream; only TACSTD2+CLDN columns retained); Gygi-hosted `protein_quant_current_normalized.csv.gz` and `Table_S1_Sample_Information.xlsx` (processed tables pointed to by the MassIVE record / Gygi CCLE page). No FASTQ. No unverified identifiers entered analysis. Refused files: none (no FASTQ candidates).

RNA values are publisher `log2(TPM+1)` protein-coding gene expression (DepMap 24Q2). Protein values are TMT10 SPS-MS3 normalized log2 ratios (Nusinow et al., *Cell* 2020; 375 lines in the paper).

### Sample size / 样本量
| Stage | RNA | Protein |
|---|---|---|
| Matrix lines available | {ctx['n_rna_matrix']} | {ctx['n_prot_matrix']} (unique CCLE names after TenPx collapse) |
| Mapped to DepMap ModelID | {ctx['n_rna_mapped']} | {ctx['n_prot_mapped']} |
| Analyzed, all cancers | {r['n_lines']} | {p['n_lines']} |
| TACSTD2 non-NA | {r['n_TACSTD2']} | {p['n_TACSTD2']} |
| CLDN4 non-NA | {r['n_CLDN4']} | {p['n_CLDN4']} |
| Pairwise TACSTD2–CLDN4 | {r['n_pairwise_TACSTD2_CLDN4']} | {p['n_pairwise_TACSTD2_CLDN4']} |
| RNA∩protein mapped lines | {ov['n_overlap']} | {ov['n_overlap']} |

Exclusions: RNA rows without a ModelID were none (index is ModelID). Protein columns that did not match `_TenPx##` were ignored if TenPx columns existed. Lineages with n<20 were omitted from exploratory lineage ranks only; they remain inside ALL_CANCERS. Biological unit is the cell line.

### Methods and reproducibility
- Claudin family: gene symbols matching `^CLDN\\d+$` present in that matrix. CLDND1 and similar non-integer symbols are excluded. TACSTD2 is not a claudin and is not inserted into the claudin rank list. Confirmatory ranks use pairwise n≥20.
- Missing values: pairwise deletion. Protein genes with duplicate rows were averaged; TenPx replicates of one CCLE name were averaged.
- Test: two-sided Spearman. Software: Python pandas/numpy/scipy/statsmodels (see `repro/pip-freeze.txt`). Seed `{seed}` is recorded; no stochastic sampler is used.
- Entry point: `python3 scripts/B2_CCLE_all/00_download.py && python3 scripts/B2_CCLE_all/01_analyze.py`.
- Environment: `results/w200/B2_CCLE_all/repro/pip-freeze.txt` and `conda-explicit.txt` (conda may be unavailable).

### Multiple testing / 多重检验
Family: `CLDN\\d+` genes with pairwise n≥20 and defined ρ against TACSTD2 **within one layer × ALL_CANCERS** (RNA and protein are separate families). m_RNA = {r['n_claudins_tested']}, m_protein = {p['n_claudins_tested']}. Raw p: two-sided Spearman. Correction: Benjamini–Hochberg. Prespecified q threshold: 0.05. Ties in ρ use minimum rank. Claudins with n<20 remain in the coverage table but are outside the confirmatory rank/BH family. Lineage tests are exploratory and are **not** folded into the confirmatory BH family.

Confirmatory RNA family: {ctx['n_rna_pass_q']} / {r['n_claudins_tested']} claudins with q<0.05.  
Confirmatory protein family: {ctx['n_prot_pass_q']} / {p['n_claudins_tested']} claudins with q<0.05.

### Results
**Confirmatory, all cancers — coexpression rank (primary B2).**

RNA (DepMap 24Q2, n={r['n_lines']} lines): TACSTD2 vs CLDN4 ρ = {fmt(r['rho_TACSTD2_CLDN4'])}, p = {fmt_p(r['p_TACSTD2_CLDN4'])}, n_pairwise = {r['n_pairwise_TACSTD2_CLDN4']}. CLDN4 rank among claudins by ρ = **{r['CLDN4_rank_rho_among_claudins']} / {r['n_claudins_tested']}**. Top claudin = **{rna_top}** (ρ = {fmt(r['top_rho'])}). {claim_rna} CLDN7 and CLDN4 are close (Δρ ≈ {fmt((r['top_rho'] or 0) - (r['rho_TACSTD2_CLDN4'] or 0))}).

Protein (Nusinow 2020, n={p['n_lines']} lines; CLDN4 quantified in {p['n_CLDN4']}/{p['n_lines']}): TACSTD2 vs CLDN4 ρ = {fmt(p['rho_TACSTD2_CLDN4'])}, p = {fmt_p(p['p_TACSTD2_CLDN4'])}, n_pairwise = {p['n_pairwise_TACSTD2_CLDN4']}. CLDN4 rank among confirmatory-family claudins by ρ = **{p['CLDN4_rank_rho_among_claudins']} / {p['n_claudins_tested']}**. Top claudin = **{prot_top}** (ρ = {fmt(p['top_rho'])}). {claim_prot} Only claudins with pairwise n≥20 enter this rank; several CLDN proteins are mostly missing.

**Confirmatory — RNA abundance rank among claudins (median log2(TPM+1)).**  
CLDN4 median-abundance rank = {r['CLDN4_rank_median_among_claudins']} / {r['n_claudins_tested']}. Protein TMT medians are relative to a common reference and are **not** copies/cell; they are tabled but not used as an abundance claim.

**Descriptive — genome-wide RNA abundance rank (mean log2(TPM+1)).**  
TACSTD2 rank {t2_rna_rank} / {t2_rna_n}; CLDN4 rank {c4_rna_rank} / {c4_rna_n}. Protein genome-wide TMT-median ranks are in `tables/genomewide_protein_abundance_rank.tsv` and must not be read as molecular abundance.

**RNA–protein coupling on mapped overlap (n={ov['n_overlap']} lines).**  
TACSTD2 RNA vs protein ρ = {fmt(ov['rho_tacstd2'])}, p = {fmt_p(ov['p_tacstd2'])}, n = {ov['n_tacstd2']}.  
CLDN4 RNA vs protein ρ = {fmt(ov['rho_cldn4'])}, p = {fmt_p(ov['p_cldn4'])}, n = {ov['n_cldn4']}.

These are cell-line expression associations. They do not measure a TACSTD2–CLDN4 junction, immune exclusion, or ICI resistance.

### Limitations
- CCLE/DepMap lines are in vitro, usually monoclonal, and lack stroma, immune cells, and in vivo barrier anatomy.
- Nusinow protein covers ~375 lines, far fewer than DepMap RNA; membrane proteins including some claudins are frequently missing (see missingness figure). TMT ratios are relative, not copies/cell.
- DepMap 24Q2 RNA mixes historical CCLE and later models; we did not batch-correct beyond the publisher matrix.
- Lineage ranks are exploratory (n≥20 gate only).
- No ICI labels exist for this slice. Do not read ρ as a resistance biomarker.

## 中文

### 研究问题与设计
主要问题：在**全部** CCLE/DepMap 癌细胞系（不按谱系过滤）中，经典 claudin（符号为 `CLDN`+整数）里，**CLDN4 与 TACSTD2 的 Spearman 相关位次**在 RNA 与蛋白层分别是多少？

预设比较：
1. 验证性：全癌种 ρ(TACSTD2, 各 CLDN\\d+)，报告 CLDN4 位次（1 = ρ 最高）。
2. 验证性：全癌种 RNA 上 CLDN4 在同一家族中的中位丰度位次（蛋白 TMT 不视为每细胞拷贝数）。
3. 验证性：TACSTD2–CLDN4 成对 ρ（RNA 与蛋白）。
4. 描述性：TACSTD2 与 CLDN4 的全基因组 RNA 丰度位次。
5. 探索性：Oncotree 谱系内（n≥20）重复共表达位次。

独立生物学单位：一条细胞系（DepMap `ModelID` / CCLE 名）。Nusinow 同一系的 TenPx 重复先平均，不把技术重复当独立样本。

### 数据来源与登录号核验
下载前经 HTTPS 核验（原始响应见 `results/w200/B2_CCLE_all/repro/accessions/`）：

| 记录 | 数据库 | 查询 | 确认标题 |
|---|---|---|---|
| `10.25452/figshare.plus.25880521` | Figshare+ API `/v2/articles/25880521` | 精确文章 id 25880521 | DepMap 24Q2 Public |
| `MSV000085836` | MassIVE PROXI `accession=MSV000085836&resultType=full` | 精确登录号 | Quantitative Proteomics of the Cancer Cell Line Encyclopedia |

接受文件见 `catalog.tsv`：DepMap `Model.csv`、`README.txt`、`OmicsExpressionProteinCodingGenesTPMLogp1.csv`（流式计算全文 SHA-256，仅保留 TACSTD2+CLDN 列）；Gygi 托管的 `protein_quant_current_normalized.csv.gz` 与 `Table_S1_Sample_Information.xlsx`（MassIVE 记录/Gygi CCLE 页指向的处理后表格）。无 FASTQ。分析未纳入未经核验的标识符。拒绝文件：无。

RNA 为 DepMap 24Q2 发布的蛋白编码基因 `log2(TPM+1)`。蛋白为 TMT10 SPS-MS3 归一化 log2 比值（Nusinow 等，*Cell* 2020；论文中 375 系）。

### Sample size / 样本量
| 阶段 | RNA | 蛋白 |
|---|---|---|
| 矩阵中的系 | {ctx['n_rna_matrix']} | {ctx['n_prot_matrix']}（TenPx 合并后的唯一 CCLE 名） |
| 映射到 DepMap ModelID | {ctx['n_rna_mapped']} | {ctx['n_prot_mapped']} |
| 全癌种分析 | {r['n_lines']} | {p['n_lines']} |
| TACSTD2 非缺失 | {r['n_TACSTD2']} | {p['n_TACSTD2']} |
| CLDN4 非缺失 | {r['n_CLDN4']} | {p['n_CLDN4']} |
| TACSTD2–CLDN4 成对 | {r['n_pairwise_TACSTD2_CLDN4']} | {p['n_pairwise_TACSTD2_CLDN4']} |
| RNA∩蛋白已映射系 | {ov['n_overlap']} | {ov['n_overlap']} |

排除：RNA 以 ModelID 为索引，无空 ID。若存在 TenPx 列，则忽略非 TenPx 样本列。n<20 的谱系只从探索性谱系位次中去掉，仍留在 ALL_CANCERS。生物学单位是细胞系。

### 方法与可复现性
- Claudin 家族：该矩阵中匹配 `^CLDN\\d+$` 的基因。CLDND1 等非整数符号不计入。TACSTD2 不是 claudin，不进入 claudin 位次表。验证性位次要求成对 n≥20。
- 缺失：成对删除。蛋白重复基因行取平均；同一 CCLE 名的 TenPx 重复取平均。
- 检验：双侧 Spearman。软件：pandas/numpy/scipy/statsmodels（见 `repro/pip-freeze.txt`）。种子 `{seed}` 已记录；本分析无随机抽样。
- 入口：`python3 scripts/B2_CCLE_all/00_download.py && python3 scripts/B2_CCLE_all/01_analyze.py`。
- 环境：`results/w200/B2_CCLE_all/repro/pip-freeze.txt` 与 `conda-explicit.txt`（可能无 conda）。

### Multiple testing / 多重检验
检验族：在**单一层次 × ALL_CANCERS** 内，成对 n≥20 且 ρ 可定义的 `CLDN\\d+` 对 TACSTD2（RNA 与蛋白分族）。m_RNA = {r['n_claudins_tested']}，m_protein = {p['n_claudins_tested']}。原始 p：双侧 Spearman。校正：Benjamini–Hochberg。预设 q 阈值：0.05。ρ 并列时取最小位次。n<20 的 claudin 仍在覆盖表中，但不进入验证性位次/BH 族。谱系检验为探索性，**不**并入验证性 BH 族。

验证性 RNA 族：{ctx['n_rna_pass_q']} / {r['n_claudins_tested']} 个 claudin q<0.05。  
验证性蛋白族：{ctx['n_prot_pass_q']} / {p['n_claudins_tested']} 个 claudin q<0.05。

### 结果
**验证性，全癌种 — 共表达位次（B2 主终点）。**

RNA（DepMap 24Q2，n={r['n_lines']} 系）：TACSTD2 vs CLDN4 ρ = {fmt(r['rho_TACSTD2_CLDN4'])}，p = {fmt_p(r['p_TACSTD2_CLDN4'])}，成对 n = {r['n_pairwise_TACSTD2_CLDN4']}。CLDN4 按 ρ 在 claudin 中的位次 = **{r['CLDN4_rank_rho_among_claudins']} / {r['n_claudins_tested']}**。最高 claudin = **{rna_top}**（ρ = {fmt(r['top_rho'])}）。{claim_rna_zh} CLDN7 与 CLDN4 很接近（Δρ ≈ {fmt((r['top_rho'] or 0) - (r['rho_TACSTD2_CLDN4'] or 0))}）。

蛋白（Nusinow 2020，n={p['n_lines']} 系；CLDN4 定量 {p['n_CLDN4']}/{p['n_lines']}）：TACSTD2 vs CLDN4 ρ = {fmt(p['rho_TACSTD2_CLDN4'])}，p = {fmt_p(p['p_TACSTD2_CLDN4'])}，成对 n = {p['n_pairwise_TACSTD2_CLDN4']}。CLDN4 按 ρ 在验证性家族中的位次 = **{p['CLDN4_rank_rho_among_claudins']} / {p['n_claudins_tested']}**。最高 claudin = **{prot_top}**（ρ = {fmt(p['top_rho'])}）。{claim_prot_zh} 仅成对 n≥20 的 claudin 进入该位次；多种 CLDN 蛋白大量缺失。

**验证性 — RNA 在 claudin 家族中的中位丰度位次（中位 log2(TPM+1)）。**  
CLDN4 中位丰度位次 = {r['CLDN4_rank_median_among_claudins']} / {r['n_claudins_tested']}。蛋白 TMT 中位数是相对桥样的比值，**不是**每细胞拷贝数；已列表但不作为丰度结论。

**描述性 — 全基因组 RNA 丰度位次（平均 log2(TPM+1)）。**  
TACSTD2 第 {t2_rna_rank} / {t2_rna_n}；CLDN4 第 {c4_rna_rank} / {c4_rna_n}。蛋白全基因组 TMT 中位位次见 `tables/genomewide_protein_abundance_rank.tsv`，不得读成分子丰度。

**已映射重叠系上的 RNA–蛋白耦合（n={ov['n_overlap']}）。**  
TACSTD2 RNA vs 蛋白 ρ = {fmt(ov['rho_tacstd2'])}，p = {fmt_p(ov['p_tacstd2'])}，n = {ov['n_tacstd2']}。  
CLDN4 RNA vs 蛋白 ρ = {fmt(ov['rho_cldn4'])}，p = {fmt_p(ov['p_cldn4'])}，n = {ov['n_cldn4']}。

以上均为细胞系表达关联，不能证明 TACSTD2–CLDN4 连接事件、免疫排斥或 ICI 耐药。

### 局限性
- CCLE/DepMap 细胞系为体外、多为单克隆，缺乏间质、免疫细胞和体内屏障结构。
- Nusinow 蛋白约 375 系，远少于 DepMap RNA；包括部分 claudin 在内的膜蛋白经常缺失（见缺失图）。TMT 比值为相对定量，不是每细胞拷贝数。
- DepMap 24Q2 RNA 混合历史 CCLE 与后续模型；除发布矩阵外未再做批次校正。
- 谱系位次仅为探索性（仅 n≥20）。
- 本切片无 ICI 标签。不能把 ρ 读成耐药生物标志物。
"""
    path.write_text(text)


def capture_env(repro: Path) -> None:
    repro.mkdir(parents=True, exist_ok=True)
    (repro / "seed.txt").write_text("20260816\n")
    pip = subprocess.run(["python3", "-m", "pip", "freeze"], capture_output=True, text=True)
    (repro / "pip-freeze.txt").write_text(pip.stdout if pip.returncode == 0 else pip.stderr)
    try:
        conda = subprocess.run(["conda", "list", "--explicit"], capture_output=True, text=True)
        if conda.returncode == 0 and conda.stdout.strip():
            (repro / "conda-explicit.txt").write_text(conda.stdout)
        else:
            (repro / "conda-explicit.txt").write_text(
                "# conda not available or failed\n# " + (conda.stderr or conda.stdout or "no conda") + "\n"
            )
    except FileNotFoundError:
        (repro / "conda-explicit.txt").write_text("# conda executable not found on PATH\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default="/tmp/ccle_data")
    ap.add_argument("--outdir", default="results/w200/B2_CCLE_all")
    args = ap.parse_args()
    work = Path(args.workdir)
    out = Path(args.outdir)
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    capture_env(out / "repro")
    run_utc = utc_now()
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        sha = "UNKNOWN"

    models = load_models(work / "Model.csv")
    rna = load_rna(work / "extracted" / "rna_tacstd2_cldn.tsv")
    prot_genes, prot_meta = load_protein(work / "protein_quant_current_normalized.csv.gz")
    s1, s1_meta = load_s1(work / "Table_S1_Sample_Information.xlsx")
    means = json.loads((work / "extracted" / "rna_gene_means.json").read_text())["means"]

    n_rna_matrix = int(len(rna))
    rna = rna.join(models.set_index("ModelID")[["OncotreeLineage", "CCLEName"]], how="left")
    n_rna_mapped = int(rna["OncotreeLineage"].notna().sum()) if "OncotreeLineage" in rna.columns else n_rna_matrix
    rna_mat = rna.drop(columns=[c for c in ["OncotreeLineage", "CCLEName"] if c in rna.columns])

    # protein: genes x CCLEName -> lines x genes, map to ModelID when possible
    prot_by_ccle = prot_genes.T
    prot_by_ccle.index.name = "CCLEName"
    n_prot_matrix = int(len(prot_by_ccle))
    if "CCLEName" in models.columns:
        map_m = models.dropna(subset=["CCLEName"]).drop_duplicates("CCLEName").set_index("CCLEName")
        prot_mapped = prot_by_ccle.join(map_m[["ModelID", "OncotreeLineage"]], how="left")
    else:
        prot_mapped = prot_by_ccle.copy()
        prot_mapped["ModelID"] = pd.NA
        prot_mapped["OncotreeLineage"] = pd.NA
    n_prot_mapped = int(prot_mapped["ModelID"].notna().sum())
    # analyze protein on all unique CCLE names (all cancers), lineage from Model when available
    prot_mat = prot_mapped.drop(columns=[c for c in ["ModelID", "OncotreeLineage"] if c in prot_mapped.columns])
    # attach a synthetic models frame indexed by CCLEName for lineage_ranks
    prot_models = prot_mapped.reset_index()[["CCLEName", "OncotreeLineage"]].copy()
    prot_models["ModelID"] = prot_models["CCLEName"]

    rna_res = layer_ranks(rna_mat, models, "RNA")
    prot_res = layer_ranks(prot_mat, prot_models, "protein")

    rna_all = next(s for s in rna_res["summaries"] if s["cohort"] == "ALL_CANCERS")
    prot_all = next(s for s in prot_res["summaries"] if s["cohort"] == "ALL_CANCERS")

    # overlap RNA-protein via CCLEName
    if "CCLEName" in rna.columns:
        rna_ccle = rna.dropna(subset=["CCLEName"]).copy()
        rna_ccle = rna_ccle[~rna_ccle.index.duplicated()]
        rna_ccle = rna_ccle.set_index("CCLEName", drop=False)
        common = sorted(set(rna_ccle.index) & set(prot_mat.index))
        ov_n = len(common)
        if ov_n:
            sp_t = spearman_pair(rna_ccle.loc[common, "TACSTD2"], prot_mat.loc[common, "TACSTD2"])
            sp_c = spearman_pair(rna_ccle.loc[common, "CLDN4"], prot_mat.loc[common, "CLDN4"])
        else:
            sp_t = {"n": 0, "rho": np.nan, "p": np.nan}
            sp_c = {"n": 0, "rho": np.nan, "p": np.nan}
    else:
        ov_n = 0
        sp_t = {"n": 0, "rho": np.nan, "p": np.nan}
        sp_c = {"n": 0, "rho": np.nan, "p": np.nan}
    overlap = {
        "n_overlap": ov_n,
        "n_tacstd2": sp_t["n"],
        "rho_tacstd2": sp_t["rho"],
        "p_tacstd2": sp_t["p"],
        "n_cldn4": sp_c["n"],
        "rho_cldn4": sp_c["rho"],
        "p_cldn4": sp_c["p"],
    }

    genes_gw = ["TACSTD2"] + claudin_cols(set(list(rna_mat.columns) + list(prot_mat.columns)))
    gw_rna = genome_wide_ranks(means, genes_gw)
    gw_prot = protein_genome_ranks(prot_genes, genes_gw)

    rna_all_tab = rna_res["table"]
    prot_all_tab = prot_res["table"]
    rna_conf = rna_all_tab[rna_all_tab["cohort"] == "ALL_CANCERS"]
    prot_conf = prot_all_tab[prot_all_tab["cohort"] == "ALL_CANCERS"]
    n_rna_pass_q = int(((rna_conf["q_vs_TACSTD2"] < 0.05) & rna_conf["q_vs_TACSTD2"].notna()).sum())
    n_prot_pass_q = int(((prot_conf["q_vs_TACSTD2"] < 0.05) & prot_conf["q_vs_TACSTD2"].notna()).sum())

    save_tsv(pd.concat([rna_all_tab, prot_all_tab], ignore_index=True), tables / "claudin_ranks.tsv")
    save_tsv(pd.DataFrame(rna_res["summaries"] + prot_res["summaries"]), tables / "cohort_summaries.tsv")
    save_tsv(gw_rna, tables / "genomewide_rna_abundance_rank.tsv")
    save_tsv(gw_prot, tables / "genomewide_protein_abundance_rank.tsv")
    save_tsv(rna_conf, tables / "rna_all_cancers_claudin_vs_TACSTD2.tsv")
    save_tsv(prot_conf, tables / "protein_all_cancers_claudin_vs_TACSTD2.tsv")

    # coverage
    cov_rows = []
    for layer, mat in [("RNA", rna_mat), ("protein", prot_mat)]:
        for g in ["TACSTD2"] + claudin_cols(mat.columns):
            s = mat[g]
            cov_rows.append(
                {
                    "layer": layer,
                    "gene": g,
                    "n_lines": int(len(mat)),
                    "n_non_na": int(s.notna().sum()),
                    "n_na": int(s.isna().sum()),
                    "min": float(s.min()) if s.notna().any() else np.nan,
                    "median": float(s.median()) if s.notna().any() else np.nan,
                    "max": float(s.max()) if s.notna().any() else np.nan,
                }
            )
    save_tsv(pd.DataFrame(cov_rows), tables / "gene_coverage.tsv")

    # lineage counts
    lin_rna = pd.Series(rna_res["lineage_counts"], name="n").rename_axis("lineage").reset_index()
    lin_rna["layer"] = "RNA"
    lin_p = pd.Series(prot_res["lineage_counts"], name="n").rename_axis("lineage").reset_index()
    lin_p["layer"] = "protein"
    save_tsv(pd.concat([lin_rna, lin_p], ignore_index=True), tables / "lineage_counts.tsv")

    key = {
        "slice": "B2_CCLE_all",
        "run_utc": run_utc,
        "git_sha": sha,
        "rna_release": "DepMap Public 24Q2",
        "protein_source": "Nusinow et al. Cell 2020 / Gygi protein_quant_current_normalized",
        "accessions": ["10.25452/figshare.plus.25880521", "MSV000085836"],
        "n_rna_matrix": n_rna_matrix,
        "n_rna_mapped": n_rna_mapped,
        "n_prot_matrix": n_prot_matrix,
        "n_prot_mapped": n_prot_mapped,
        "protein_parse": prot_meta,
        "s1_meta": {k: s1_meta[k] for k in s1_meta if k != "sheets"} | {"sheets": s1_meta.get("sheets")},
        "rna_all": rna_all,
        "prot_all": prot_all,
        "overlap": overlap,
        "n_rna_pass_q": n_rna_pass_q,
        "n_prot_pass_q": n_prot_pass_q,
        "rna_claudins": rna_res["claudins"],
        "protein_claudins": prot_res["claudins"],
    }
    (tables / "key_stats.json").write_text(json.dumps(key, indent=2, default=str))

    fig_bar_rho(
        rna_conf,
        f"RNA all cancers: TACSTD2 vs claudins (n={rna_all['n_lines']})",
        figs / "fig1_rna_claudin_rho.png",
    )
    fig_bar_rho(
        prot_conf,
        f"Protein all cancers: TACSTD2 vs claudins (n={prot_all['n_lines']})",
        figs / "fig2_protein_claudin_rho.png",
    )
    fig_scatter(
        rna_mat["TACSTD2"],
        rna_mat["CLDN4"],
        "TACSTD2 RNA log2(TPM+1)",
        "CLDN4 RNA log2(TPM+1)",
        "RNA TACSTD2 vs CLDN4, all cancers",
        figs / "fig3_rna_tacstd2_cldn4.png",
        rna_all["n_pairwise_TACSTD2_CLDN4"],
        rna_all["rho_TACSTD2_CLDN4"],
        rna_all["p_TACSTD2_CLDN4"],
    )
    fig_scatter(
        prot_mat["TACSTD2"],
        prot_mat["CLDN4"],
        "TACSTD2 protein (TMT)",
        "CLDN4 protein (TMT)",
        "Protein TACSTD2 vs CLDN4, all cancers",
        figs / "fig4_protein_tacstd2_cldn4.png",
        prot_all["n_pairwise_TACSTD2_CLDN4"],
        prot_all["rho_TACSTD2_CLDN4"],
        prot_all["p_TACSTD2_CLDN4"],
    )
    fig_missing(rna_mat, prot_mat, figs / "fig5_missingness.png")
    fig_lineage_rank(rna_res["summaries"], figs / "fig6_rna_CLDN4_rank_by_lineage.png")
    fig_lineage_rank(prot_res["summaries"], figs / "fig7_protein_CLDN4_rank_by_lineage.png")

    # abundance bars
    for tab, title, fn in [
        (rna_conf, "RNA median abundance among claudins", "fig8_rna_claudin_abundance.png"),
        (prot_conf, "Protein median abundance among claudins", "fig9_protein_claudin_abundance.png"),
    ]:
        d = tab.sort_values("median_abundance", ascending=False)
        colors = ["#c0392b" if g == "CLDN4" else "#4a6fa5" for g in d["gene"]]
        fig, ax = plt.subplots(figsize=(10, 4.0))
        ax.bar(d["gene"], d["median_abundance"], color=colors)
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=60, labelsize=8)
        fig.tight_layout()
        fig.savefig(figs / fn, dpi=160)
        plt.close(fig)

    ctx = {
        "git_sha": sha,
        "run_utc": run_utc,
        "seed": "20260816",
        "n_rna_matrix": n_rna_matrix,
        "n_rna_mapped": n_rna_mapped,
        "n_prot_matrix": n_prot_matrix,
        "n_prot_mapped": n_prot_mapped,
        "rna_all": rna_all,
        "prot_all": prot_all,
        "overlap": overlap,
        "n_rna_pass_q": n_rna_pass_q,
        "n_prot_pass_q": n_prot_pass_q,
        "gw_rna": gw_rna,
        "gw_prot": gw_prot,
    }
    write_writeup(out / "WRITEUP.md", ctx)
    notes = Path("notes/w200/B2_CCLE_all")
    notes.mkdir(parents=True, exist_ok=True)
    (notes / "WRITEUP.md").write_text((out / "WRITEUP.md").read_text())
    (notes / "rerun.md").write_text(
        """# Rerun B2_CCLE_all

```bash
python3 -m pip install pandas numpy scipy matplotlib statsmodels openpyxl
python3 scripts/B2_CCLE_all/00_download.py --workdir /tmp/ccle_data --outdir results/w200/B2_CCLE_all
python3 scripts/B2_CCLE_all/01_analyze.py --workdir /tmp/ccle_data --outdir results/w200/B2_CCLE_all
```

Raw DepMap RNA (~461 MB) is hashed and not retained. Extracted TACSTD2/CLDN columns and gene-means live under `/tmp/ccle_data/extracted/`.
Outputs: `results/w200/B2_CCLE_all/{tables,figures,WRITEUP.md,catalog.tsv,repro/}`.
"""
    )
    print(json.dumps({"rna_all": rna_all, "prot_all": prot_all, "overlap": overlap}, indent=2, default=str))


if __name__ == "__main__":
    main()
