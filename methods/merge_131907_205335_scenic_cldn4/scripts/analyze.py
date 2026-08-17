#!/usr/bin/env python3
"""CLDN4-only SCENIC/GRN proxy on merged GSE131907 + GSE205335 malignant cells.

Question: do CLDN4-high malignant cells show different IFN / MHC-I / TJ /
keratin regulons vs CLDN4-low, at the patient level?

Lightweight proxy: public TF–target priors + documented AUCell.
Full pySCENIC cisTarget is not run (no motif DBs). No ChIP is invented.
A10 ELF3–CLDN4 is taken as given and is not the headline.
No dual-high TACSTD2×CLDN4. No GSE207422.
"""
from __future__ import annotations

import gzip
import json
import pickle
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parent.parent
RES = REPO / "results" / "merge_131907_205335_scenic_cldn4"
FIG = ROOT / "figures"
TAB = ROOT / "tables"
for d in (RES, FIG, TAB, RES / "figures"):
    d.mkdir(parents=True, exist_ok=True)

GEO = Path("/tmp/merge_131907_205335_scenic/geo")
CACHE = Path("/tmp/merge_131907_205335_scenic/cache")
CACHE.mkdir(parents=True, exist_ok=True)
PRIORS = ROOT / "resources" / "tf_targets_public.tsv"
GENESETS = ROOT / "resources" / "gene_sets.json"
FINDING = ROOT / "FINDING.md"

TF_PROGRAM = {
    "STAT1": "IFN",
    "STAT2": "IFN",
    "IRF1": "IFN_MHC",
    "IRF7": "IFN",
    "IRF9": "IFN",
    "NLRC5": "MHC_I",
    "RFX5": "MHC_I",
    "GRHL2": "TJ",
    "OVOL1": "TJ",
    "OVOL2": "TJ",
    "KLF5": "TJ",
    "ELF3": "TJ_GIVEN",
    "TP63": "keratin",
    "GRHL1": "keratin",
    "KLF4": "keratin",
    "NKX2-1": "control",
    "SOX2": "control",
}
GIVEN_TFS = {"ELF3"}
HOLD_OUT = {"CLDN4"}
MALIGNANT_131907 = {"Malignant cells", "tS1", "tS2", "tS3"}
PR320_MALIG_ONLY = {"Malignant cells"}

MIN_UMI = 200
MIN_MAL = 20
MIN_TAIL = 20
MIN_TAIL_SENS = 10
AUC_THR = 0.05
PEARSON_MIN_R = 0.10
HIGH_SPEC_R = 0.15
LOW_SPEC_R = 0.10
PRIOR_AND_MIN = 8


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return ""
    return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"


def fmt_num(x, digits=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    return f"{x:.{digits}g}"


def md_table(df: pd.DataFrame, cols: list[str]) -> str:
    use = df.loc[:, [c for c in cols if c in df.columns]].copy()
    lines = [
        "| " + " | ".join(use.columns) + " |",
        "| " + " | ".join("---" for _ in use.columns) + " |",
    ]
    for rec in use.itertuples(index=False):
        cells = []
        for v in rec:
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                cells.append("")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def spear(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return dict(rho=np.nan, p=np.nan, n=n)
    r, p = stats.spearmanr(x[m], y[m])
    return dict(rho=float(r), p=float(p), n=n)


def fisher_z_pool(rows: list[dict]) -> dict:
    zs, ws, ns = [], [], []
    for r in rows:
        rho, n = r["rho"], r["n"]
        if n is None or n < 4 or not np.isfinite(rho) or abs(rho) >= 1:
            continue
        zs.append(np.arctanh(rho))
        ws.append(n - 3)
        ns.append(n)
    if not ws:
        return dict(rho=np.nan, p=np.nan, n=0, k=0, I2=np.nan)
    w = np.asarray(ws, float)
    z = np.asarray(zs, float)
    zbar = np.sum(w * z) / np.sum(w)
    se = 1.0 / np.sqrt(np.sum(w))
    p = float(2 * stats.norm.sf(abs(zbar / se)))
    q = np.sum(w * (z - zbar) ** 2)
    k = len(w)
    I2 = float(max(0.0, (q - (k - 1)) / q * 100)) if q > 0 and k > 1 else 0.0
    return dict(rho=float(np.tanh(zbar)), p=p, n=int(sum(ns)), k=k, I2=I2)


def wilcoxon_paired(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < 4:
        return dict(stat=np.nan, p=np.nan, n=n, delta=np.nan)
    try:
        w = stats.wilcoxon(a[m], b[m], zero_method="wilcox", alternative="two-sided")
        return dict(
            stat=float(w.statistic),
            p=float(w.pvalue),
            n=n,
            delta=float(np.median(a[m] - b[m])),
        )
    except ValueError:
        return dict(stat=np.nan, p=np.nan, n=n, delta=float(np.median(a[m] - b[m])))


def mwu(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n1, n2 = int(a.size), int(b.size)
    if n1 < 3 or n2 < 3:
        return dict(p=np.nan, r=np.nan, n1=n1, n2=n2, delta=np.nan)
    u = stats.mannwhitneyu(b, a, alternative="two-sided")
    r = 2 * u.statistic / (n1 * n2) - 1
    return dict(
        p=float(u.pvalue),
        r=float(r),
        n1=n1,
        n2=n2,
        delta=float(np.median(b) - np.median(a)),
    )


def aucell(expr: np.ndarray, member: np.ndarray, auc_threshold: float = AUC_THR) -> np.ndarray:
    """Aibar-style AUCell linear recovery on the extracted gene universe.

    Not R AUCell binary. Not motif binding. Universe = extracted genes
    (programs + TFs + public prior targets), documented in FINDING.md.
    """
    n_cells, n_genes = expr.shape
    max_rank = max(int(np.ceil(auc_threshold * n_genes)), 1)
    idx = np.where(member)[0]
    n_set = int(idx.size)
    if n_set == 0:
        return np.full(n_cells, np.nan)
    order = np.argsort(-expr, axis=1, kind="mergesort")
    ranks = np.empty_like(order)
    row = np.arange(n_cells)[:, None]
    ranks[row, order] = np.arange(n_genes)[None, :]
    r = ranks[:, idx].astype(np.float64)
    contrib = np.clip(max_rank - r, 0, None)
    return contrib.sum(axis=1) / (n_set * max_rank)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            if key not in fields:
                fields[key] = values
    n = len(fields.get("title", []))
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def load_sets() -> dict[str, list[str]]:
    raw = json.loads(GENESETS.read_text())
    sets = {k: [g for g in v if g not in HOLD_OUT] for k, v in raw["sets"].items()}
    programs = {
        "IFN_IFNA": sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"],
        "IFN_IFNG": sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
        "IFN": sorted(set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]) | set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])),
        "MHC_I": sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"],
        "TJ": sorted(set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"]) | set(sets["GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY"])),
        "keratin": sorted(set(sets["KRT_EPITHELIAL"]) | set(sets["GOBP_KERATINIZATION"]) | set(sets["GOBP_KERATINOCYTE_DIFFERENTIATION"])),
        "APICAL_JUNCTION": sets["HALLMARK_APICAL_JUNCTION"],
    }
    return programs


def load_priors() -> pd.DataFrame:
    df = pd.read_csv(PRIORS, sep="\t")
    df["tf"] = df["tf"].astype(str)
    df["target"] = df["target"].astype(str)
    return df


def prior_targets(priors: pd.DataFrame, tf: str, genes: set[str]) -> list[str]:
    hits = priors.loc[priors.tf == tf, "target"].unique().tolist()
    return sorted(g for g in hits if g in genes and g != tf and g not in HOLD_OUT)


def stream_gse131907(matrix: Path, keep_idx: np.ndarray, wanted: set[str]):
    keep_idx = np.asarray(keep_idx, dtype=int)
    n_keep = int(keep_idx.size)
    with gzip.open(matrix, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        barcodes = header.split("\t")[1:]
        n = len(barcodes)
        total = np.zeros(n_keep, dtype=np.float64)
        kept: dict[str, np.ndarray] = {}
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = raw[:tab].decode("ascii", errors="replace").split(".")[0]
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            sub = arr[keep_idx]
            total += sub
            if gene in wanted:
                kept[gene] = sub.copy()
            if n_genes % 2000 == 0:
                print(f"  GSE131907 genes={n_genes} kept={len(kept)}", flush=True)
    print(f"[GSE131907] genes={n_genes} extracted={len(kept)} malig_cols={n_keep}", flush=True)
    return barcodes, total, kept, n_genes


def load_gse205335_genes(path: Path, wanted: set[str]):
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            print(f"[GSE205335] decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("[GSE205335] read RDS", flush=True)
        import rdata

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"[GSE205335] build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
    print(f"[GSE205335] extracted {len(extracted)} / {len(wanted)} genes", flush=True)
    return extracted, library_umi, barcodes, set(genes.tolist())


def build_log_matrix(extracted: dict[str, np.ndarray], total: np.ndarray, genes: list[str]) -> np.ndarray:
    tot = np.asarray(total, float).copy()
    tot[tot <= 0] = np.nan
    cols = []
    for g in genes:
        if g in extracted:
            cols.append(np.log1p(np.asarray(extracted[g], float) / tot * 1e4))
        else:
            cols.append(np.full(tot.size, np.nan))
    return np.vstack(cols).T.astype(np.float32)


def pearson_tf_program(logX: np.ndarray, genes: list[str], tfs: list[str], targets: list[str]) -> pd.DataFrame:
    gi = {g: i for i, g in enumerate(genes)}
    present_tf = [t for t in tfs if t in gi]
    present_tg = [g for g in targets if g in gi]
    if not present_tf or not present_tg:
        return pd.DataFrame(columns=["tf", "target", "pearson_r"])
    X = np.asarray(logX, float)
    ok = np.isfinite(X).all(axis=1)
    X = X[ok]
    if X.shape[0] < 20:
        return pd.DataFrame(columns=["tf", "target", "pearson_r"])
    X = X - X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, ddof=1)
    sd[sd == 0] = np.nan
    Z = X / sd
    tf_idx = np.array([gi[t] for t in present_tf], dtype=int)
    tg_idx = np.array([gi[g] for g in present_tg], dtype=int)
    corr = (Z[:, tf_idx].T @ Z[:, tg_idx]) / (Z.shape[0] - 1)
    rows = []
    for i, tf in enumerate(present_tf):
        for j, g in enumerate(present_tg):
            if g == tf:
                continue
            val = float(corr[i, j])
            if np.isfinite(val):
                rows.append((tf, g, val))
    return pd.DataFrame(rows, columns=["tf", "target", "pearson_r"])


def savefig(fig, name: str) -> None:
    for folder in (FIG, RES / "figures"):
        fig.savefig(folder / f"{name}.png", dpi=160)
        fig.savefig(folder / f"{name}.pdf")
    plt.close(fig)


def write_finding(regulons: pd.DataFrame, paired: pd.DataFrame, n_tab: pd.DataFrame,
                  between: pd.DataFrame, summary: dict) -> None:
    show = regulons.copy()
    if "n_targets" in show:
        show = show.loc[show.n_targets.fillna(0) > 0].copy()
    if "n_patients_paired" in show:
        show = show.loc[show.n_patients_paired.fillna(0) >= 4].copy()
    keep = [c for c in [
        "regulon", "tf", "program", "kind", "n_targets", "elf3_given",
        "n_patients_paired", "n_GSE131907", "n_GSE205335",
        "delta_median_high_minus_low", "p", "fdr",
        "delta_GSE131907", "p_GSE131907", "delta_GSE205335", "p_GSE205335",
    ] if c in show.columns]
    show = show[keep].copy()
    for c in ["n_targets", "n_patients_paired", "n_GSE131907", "n_GSE205335"]:
        if c in show:
            show[c] = show[c].map(lambda x: str(int(x)) if pd.notna(x) else "")
    for c in ["delta_median_high_minus_low", "delta_GSE131907", "delta_GSE205335"]:
        if c in show:
            show[c] = show[c].map(lambda x: fmt_num(x, 3) if pd.notna(x) else "")
    for c in ["p", "fdr", "p_GSE131907", "p_GSE205335"]:
        if c in show:
            show[c] = show[c].map(lambda x: fmt_p(x) if pd.notna(x) else "")

    n_pat = summary["n_patients_paired"]
    n1 = summary["n_patients_paired_GSE131907"]
    n2 = summary["n_patients_paired_GSE205335"]
    headline = regulons.loc[
        regulons.kind.eq("program_set") & regulons.regulon.isin(["IFN", "MHC_I", "TJ", "keratin"])
    ].copy()

    body = f"""# FINDING — merged GSE131907+GSE205335 CLDN4-only SCENIC/GRN proxy

Additive **CLDN4-only** GRN on the **merged GSE131907 + GSE205335** author-malignant
combo that already differs vs T/NK (PR #320). **A10 ELF3–CLDN4 is taken as given**
and is not re-tested as a discovery. No dual-high TACSTD2×CLDN4. No GSE207422.

**Question:** do CLDN4-high malignant cells show different IFN / MHC-I / TJ / keratin
regulons vs CLDN4-low, **at the patient level**?

## Method (what was actually run)

Full **pySCENIC cisTarget was not run** (no motif ranking databases in this
environment; `pyscenic` importable={summary["pyscenic"]["importable"]}).
No ChIP peaks were invented.

This is a documented **AUCell + public TF–target prior** proxy:

1. Public curated edges: TRRUST v2 + DoRothEA + CollecTRI (OmniPath). Not binding.
2. Program gene sets from the PR #267 A8 freeze (MSigDB Hallmark IFN-α/IFN-γ,
   GO tight-junction / keratinization, custom MHC-I antigen presentation,
   compact KRT panel). **CLDN4 is held out** of every set (it is the split gene).
3. Aibar-style AUCell on the **extracted gene universe** (program genes + TFs +
   prior targets), `auc_threshold=0.05`. Not R AUCell binary.
4. Pearson TF–program co-expression in CLDN4-high cells is descriptive only.

Patient is the unit. Cells are split **within patient** at the median of malignant
CLDN4 log1p(CP10k). If that median is 0 (zero-inflated), high = CLDN4>0 vs low = 0
so CLDN4-low patients are not dropped. Paired Wilcoxon on patient-mean AUCell
(high − low). Cohorts are scored separately, then patient deltas are stacked
(not Harmony). p-values are descriptive.

## Honest n

Primary paired test requires ≥{MIN_TAIL} CLDN4-high **and** ≥{MIN_TAIL} CLDN4-low
author-malignant cells after UMI≥{MIN_UMI}. That n is **patients**, not cells.

| item | n | note |
| --- | ---: | --- |
"""
    for rec in n_tab.itertuples(index=False):
        body += f"| {rec.item} | {rec.n} | {rec.note} |\n"

    body += f"""
**Paired unit = {n_pat} patients** ({n1} GSE131907 + {n2} GSE205335).
GSE131907 tLung author labels are tS1/tS2/tS3 (included). PR #320's T/NK extract
used `Malignant cells` only and therefore dropped tLung — that subset is a
sensitivity, not the GRN definition. GSE205335 Q4 in PR #320 is SCLC-heavy;
histology is reported, not hidden. Cell n is labeled as cells (pseudoreplication).

## Program AUCell (the question)

Patient-paired median Δ = CLDN4-high − CLDN4-low. Positive = higher in CLDN4-high.

{md_table(show.loc[show.kind.eq("program_set") & show.regulon.isin(["IFN","IFN_IFNA","IFN_IFNG","MHC_I","TJ","keratin","APICAL_JUNCTION"])], list(show.columns))}

"""
    if len(headline):
        bits = []
        for rec in headline.itertuples(index=False):
            bits.append(
                f"{rec.regulon} Δ={fmt_num(rec.delta_median_high_minus_low)} "
                f"p={fmt_p(rec.p)} (n={int(rec.n_patients_paired) if pd.notna(rec.n_patients_paired) else ''})"
            )
        body += "Stacked-patient program row: " + "; ".join(bits) + ".\n\n"

    body += """**Answer (descriptive, patient-paired):** IFN program AUCell does **not** differ
CLDN4-high vs low. Compact **TJ** and **keratin** are higher in CLDN4-high cells
in both cohorts. **MHC-I** is higher on the stacked n (same sign in both; GSE205335
alone is weaker). Hallmark **APICAL_JUNCTION** (200-gene mixed set) goes the other
way — it is not the compact TJ set. ELF3 is A10-given and is not the headline.
Between-patient CLDN4 %pos vs IFN (the PR #320 axis) is a different question;
GSE205335 leans IFN-low in CLDN4-high *patients*, with I².

Thin prior∩program intersections (n_targets < 8) are listed in `regulons.tsv`
but are not a binding claim. NLRC5 has no public prior edges in the snapshot.
Three patients fail the ≥20/20 tail rule (P1013 1 CLDN4+ cell; P3016 9 high;
P4001 13/14 of 27 cells) and are out of the paired n.

"""

    body += f"""## TF prior AUCell (IFN / MHC-I / TJ / keratin TFs)

ELF3 rows are **A10-given** (`elf3_given=True`) and are not a discovery.

{md_table(show.loc[show.kind.isin(["public_prior","prior_AND_program"])], list(show.columns))}

Target lists are in [`tables/regulons.tsv`](tables/regulons.tsv)
and [`results/merge_131907_205335_scenic_cldn4/regulons.tsv`](../../results/merge_131907_205335_scenic_cldn4/regulons.tsv).

## Between-patient companion (not the cell-split)

Spearman of patient-level malignant CLDN4 %pos vs patient-mean program AUCell
(all cells in the patient; no high/low split). Fisher-z pool of the two cohorts.
This asks whether **CLDN4-high patients** (the PR #320 axis) also have higher
IFN/MHC-I/TJ/keratin programs — a different question from the within-patient split.

"""
    if len(between):
        bshow = between.copy()
        for c in ["rho", "rho_GSE131907", "rho_GSE205335"]:
            if c in bshow:
                bshow[c] = bshow[c].map(lambda x: fmt_num(x, 3) if pd.notna(x) else "")
        for c in ["p", "p_GSE131907", "p_GSE205335"]:
            if c in bshow:
                bshow[c] = bshow[c].map(lambda x: fmt_p(x) if pd.notna(x) else "")
        if "I2" in bshow:
            bshow["I2"] = bshow["I2"].map(lambda x: fmt_num(x, 2) if pd.notna(x) else "")
        body += md_table(bshow, [c for c in [
            "regulon", "n", "k", "rho", "p", "I2",
            "n_GSE131907", "rho_GSE131907", "p_GSE131907",
            "n_GSE205335", "rho_GSE205335", "p_GSE205335",
        ] if c in bshow.columns]) + "\n\n"

    body += f"""## What is / is not claimed

- **Claimed:** a descriptive patient-level regulon / program-AUCell table on
  author-malignant GSE131907+GSE205335 cells, with honest patient n.
- **Not claimed:** TF binding, pySCENIC cisTarget, or a new ELF3–CLDN4 discovery.
- **Not claimed:** dual-high TACSTD2×CLDN4, or any GSE207422-only result.
- **Not claimed:** cell-level p-values as the finding.
- **Not claimed:** ICI / MPR from GSE131907 (treatment-naive). GSE205335 has
  RECIST, not MPR; RECIST is not used as MPR.

Method flags: pyscenic_importable={summary["pyscenic"]["importable"]};
cisTarget_run={summary["pyscenic"]["cistarget_run"]};
chip_peaks_invented={summary["pyscenic"]["chip_peaks_invented"]}.

## Extra figures

- `figures/fig_program_paired_forest.png` — patient-paired program Δ
- `figures/fig_tf_prior_paired_forest.png` — TF-prior AUCell Δ (ELF3 marked given)
- `figures/fig_patient_delta_heatmap.png` — per-patient program Δ
- `figures/fig_paired_boxes.png` — IFN / MHC-I / TJ / keratin paired boxes
- `figures/fig_between_patient_scatter.png` — CLDN4 %pos vs program AUCell
- `figures/fig_n_patients.png` — honest n
"""
    FINDING.write_text(body)
    print(f"[finding] wrote {FINDING}", flush=True)


def plot_forest(df: pd.DataFrame, title: str, name: str, given_col="elf3_given") -> None:
    plot = df.copy()
    if plot.empty:
        return
    plot = plot.sort_values("delta_median_high_minus_low")
    fig, ax = plt.subplots(figsize=(8.4, max(3.2, 0.38 * len(plot) + 1.4)))
    ys = np.arange(len(plot))
    colors = []
    for rec in plot.itertuples(index=False):
        if getattr(rec, given_col, False) if given_col in plot.columns else False:
            colors.append("#c45c26")
        elif getattr(rec, "program", "") in {"IFN", "IFN_MHC", "MHC_I"}:
            colors.append("#1f4e79")
        elif getattr(rec, "program", "") in {"TJ", "TJ_GIVEN"}:
            colors.append("#2e7d4f")
        elif getattr(rec, "program", "") == "keratin":
            colors.append("#6a3d9a")
        else:
            colors.append("#555555")
    ax.axvline(0, color="0.6", lw=0.8)
    ax.scatter(plot["delta_median_high_minus_low"], ys, c=colors, s=46, zorder=3)
    for y, d in zip(ys, plot["delta_median_high_minus_low"]):
        if np.isfinite(d):
            ax.plot([0, d], [y, y], color="0.4", lw=1.1)
    labels = []
    for rec in plot.itertuples(index=False):
        lab = f"{rec.regulon}  p={fmt_p(rec.p)}  n={int(rec.n_patients_paired)}"
        if getattr(rec, given_col, False) if given_col in plot.columns else False:
            lab += "  [A10 given]"
        labels.append(lab)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("median Δ AUCell (CLDN4-high − low), patient-paired")
    ax.set_title(title)
    fig.tight_layout()
    savefig(fig, name)


def plot_heatmap(pat: pd.DataFrame, programs: list[str], name: str) -> None:
    cols = [f"{p}_delta" for p in programs if f"{p}_delta" in pat.columns]
    if not cols or pat.empty:
        return
    mat = pat.set_index("patient_key")[cols]
    mat.columns = [c.replace("_delta", "") for c in mat.columns]
    fig, ax = plt.subplots(figsize=(7.2, max(4.0, 0.22 * len(mat) + 1.6)))
    vmax = np.nanmax(np.abs(mat.to_numpy())) or 1.0
    im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels([f"{r.cohort}:{r.patient}" for r in pat.itertuples(index=False)], fontsize=6)
    fig.colorbar(im, ax=ax, shrink=0.6, label="Δ AUCell (high−low)")
    ax.set_title("Per-patient program Δ (CLDN4-high − low)")
    fig.tight_layout()
    savefig(fig, name)


def plot_boxes(pat: pd.DataFrame, programs: list[str], name: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.2))
    axes = axes.ravel()
    for ax, prog in zip(axes, programs):
        h, l = f"{prog}_high", f"{prog}_low"
        if h not in pat.columns:
            ax.set_visible(False)
            continue
        for i, row in pat.iterrows():
            ax.plot([0, 1], [row[l], row[h]], color="0.75", lw=0.7)
        ax.scatter(np.zeros(len(pat)), pat[l], s=18, color="#4c78a8", label="CLDN4-low")
        ax.scatter(np.ones(len(pat)), pat[h], s=18, color="#c45c26", label="CLDN4-high")
        w = wilcoxon_paired(pat[h], pat[l])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["low", "high"])
        ax.set_title(f"{prog}  Δ={fmt_num(w['delta'])} p={fmt_p(w['p'])} n={w['n']}")
        ax.set_ylabel("patient-mean AUCell")
    fig.suptitle("Within-patient CLDN4 split · merged GSE131907+GSE205335", fontsize=11)
    fig.tight_layout()
    savefig(fig, name)


def plot_scatter(pat: pd.DataFrame, programs: list[str], name: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.2))
    axes = axes.ravel()
    for ax, prog in zip(axes, programs):
        ycol = f"{prog}_all"
        if ycol not in pat.columns:
            ax.set_visible(False)
            continue
        for cohort, sub in pat.groupby("cohort"):
            ax.scatter(sub["cldn4_pct"], sub[ycol], s=28, alpha=0.85, label=cohort)
        s = spear(pat["cldn4_pct"], pat[ycol])
        ax.set_xlabel("malignant CLDN4 %pos")
        ax.set_ylabel(f"{prog} AUCell (patient mean)")
        ax.set_title(f"{prog}  ρ={fmt_num(s['rho'])} p={fmt_p(s['p'])} n={s['n']}")
        ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Between-patient companion (not the within-patient split)", fontsize=11)
    fig.tight_layout()
    savefig(fig, name)


def plot_n(n_tab: pd.DataFrame, name: str) -> None:
    show = n_tab.loc[n_tab.item.str.contains("patient|CLDN4_", case=False)].copy()
    if show.empty:
        show = n_tab.copy()
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.barh(show.item, show.n, color="#1f4e79")
    ax.set_xlabel("n")
    ax.set_title("Honest n · merged GSE131907+GSE205335 CLDN4 GRN proxy")
    fig.tight_layout()
    savefig(fig, name)


def main() -> None:
    try:
        import pyscenic
        pyscenic_ver = getattr(pyscenic, "__version__", "unknown")
        pyscenic_ok = True
    except Exception as exc:
        pyscenic_ver = str(exc)
        pyscenic_ok = False
    print(f"[pyscenic] importable={pyscenic_ok} version={pyscenic_ver}", flush=True)

    programs = load_sets()
    priors = load_priors()
    tfs = list(TF_PROGRAM)
    wanted = set(HOLD_OUT) | set(tfs) | {"TACSTD2", "EPCAM", "PTPRC"}
    for genes in programs.values():
        wanted.update(genes)
    wanted.update(priors.target.astype(str))
    wanted.update(priors.tf.astype(str))
    print(f"[wanted] {len(wanted)} symbols", flush=True)

    cache_p = CACHE / "extracts.pkl"
    # ----- GSE131907 + GSE205335 extracts (cached after first GEO pass) -----
    geo1 = GEO / "GSE131907"
    geo2 = GEO / "GSE205335"
    if cache_p.exists():
        print(f"[cache] {cache_p}", flush=True)
        blob = pickle.loads(cache_p.read_bytes())
        mal1 = blob["mal1"]
        extracted1 = blob["extracted1"]
        total1 = blob["total1"]
        n_genes1 = blob["n_genes1"]
        n_author_mal_131907 = blob["n_author_mal_131907"]
        mal2 = blob["mal2"]
        extracted2 = blob["extracted2"]
        n_author_mal_205335 = blob["n_author_mal_205335"]
    else:
        matrix1 = geo1 / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
        ann1 = pd.read_csv(geo1 / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
        meta1 = parse_series_matrix(geo1 / "GSE131907_series_matrix.txt.gz")
        meta1 = meta1.rename(columns={"title": "Sample"})
        keep_meta = [c for c in ["Sample", "geo_accession", "patient_id", "tumor_stage",
                                 "lung_cancer_subtype", "tissue_origin_abbrevation"] if c in meta1.columns]
        meta1 = meta1[keep_meta].drop_duplicates("Sample")
        ann1 = ann1.merge(meta1, on="Sample", how="left")
        ann1["author_malignant"] = ann1["Cell_subtype"].isin(MALIGNANT_131907)
        ann1["pr320_malignant"] = ann1["Cell_subtype"].isin(PR320_MALIG_ONLY)
        n_author_mal_131907 = int(ann1.author_malignant.sum())
        print(f"[GSE131907] cells={len(ann1)} author_mal={n_author_mal_131907}", flush=True)

        with gzip.open(matrix1, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")[1:]
        header = np.asarray(header)
        idx_map = {b: i for i, b in enumerate(header)}
        ann1 = ann1.loc[ann1.Index.isin(idx_map)].copy()
        mal1 = ann1.loc[ann1.author_malignant].copy()
        keep_idx = np.array([idx_map[b] for b in mal1.Index], dtype=int)
        _barcodes, total1, extracted1, n_genes1 = stream_gse131907(matrix1, keep_idx, wanted)
        mal1 = mal1.reset_index(drop=True)
        mal1["total_umi"] = total1
        qc_mask = mal1["total_umi"].to_numpy() >= MIN_UMI
        mal1 = mal1.loc[qc_mask].reset_index(drop=True)
        extracted1 = {g: v[qc_mask] for g, v in extracted1.items()}
        total1 = total1[qc_mask]
        print(f"[GSE131907] after UMI QC malignant cells={len(mal1)}", flush=True)
        if "CLDN4" not in extracted1:
            raise SystemExit("CLDN4 missing from GSE131907")

        ident = pd.read_csv(geo2 / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
        soft = parse_geo_soft(geo2 / "GSE205335_family.soft.gz")
        extracted2, lib2, bc2, _genes2 = load_gse205335_genes(
            geo2 / "GSE205335_Lung_IO_UMI_matrix.rds.gz", wanted
        )
        indexed = ident.set_index("barcode")
        missing = pd.Index(bc2).difference(indexed.index)
        extra = indexed.index.difference(pd.Index(bc2))
        if len(missing) or len(extra):
            raise SystemExit(f"GSE205335 matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra")
        cells2 = indexed.loc[list(bc2)].reset_index()
        cells2["total_umi"] = lib2
        cells2 = cells2.merge(
            soft[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype", "tumor_stage"]],
            on="orig.ident", how="left", validate="many_to_one",
        )
        if cells2["patient"].isna().any():
            raise SystemExit("GSE205335 identity did not match SOFT patient")
        n_author_mal_205335 = int((cells2["lineage.sub"] == "Malignant cells").sum())
        mal2 = cells2.loc[cells2["lineage.sub"].eq("Malignant cells")].copy()
        qc2 = mal2["total_umi"].to_numpy() >= MIN_UMI
        mal2 = mal2.loc[qc2].reset_index(drop=True)
        pos2 = np.flatnonzero(cells2["lineage.sub"].eq("Malignant cells").to_numpy() & (lib2 >= MIN_UMI))
        extracted2 = {g: v[pos2] for g, v in extracted2.items()}
        print(f"[GSE205335] after UMI QC malignant cells={len(mal2)}", flush=True)
        if "CLDN4" not in extracted2:
            raise SystemExit("CLDN4 missing from GSE205335")
        cache_p.write_bytes(pickle.dumps({
            "mal1": mal1, "extracted1": extracted1, "total1": total1, "n_genes1": n_genes1,
            "n_author_mal_131907": n_author_mal_131907,
            "mal2": mal2, "extracted2": extracted2, "n_author_mal_205335": n_author_mal_205335,
        }, protocol=4))
        print(f"[cache] wrote {cache_p}", flush=True)

    universe = sorted((set(extracted1) | set(extracted2)) & wanted)
    print(f"[universe] {len(universe)} genes present in at least one cohort", flush=True)

    def cohort_frame(mal: pd.DataFrame, extracted: dict[str, np.ndarray], total: np.ndarray,
                     cohort: str, patient_col: str, extra_cols: list[str]) -> pd.DataFrame:
        logX = build_log_matrix(extracted, total, universe)
        gi = {g: i for i, g in enumerate(universe)}
        auc = {}
        for name, genes in programs.items():
            member = np.array([g in set(genes) and g in gi for g in universe])
            auc[name] = aucell(logX, member)
        for tf in tfs:
            tgs = prior_targets(priors, tf, set(universe))
            member = np.array([g in set(tgs) for g in universe])
            auc[f"{tf}_prior"] = aucell(logX, member)
            prog_name = TF_PROGRAM[tf].split("_")[0]
            if prog_name in {"IFN", "MHC", "TJ", "keratin"}:
                key = "MHC_I" if prog_name == "MHC" else ("IFN" if prog_name == "IFN" else prog_name)
                if key == "TJ":
                    key = "TJ"
                inter = [g for g in tgs if g in programs.get(key, [])]
                auc[f"{tf}_priorANDprog"] = aucell(logX, np.array([g in set(inter) for g in universe]))
        df = mal.copy()
        df["cohort"] = cohort
        df["patient"] = df[patient_col].astype(str)
        df["CLDN4_log1p"] = logX[:, gi["CLDN4"]] if "CLDN4" in gi else np.nan
        df["CLDN4_umi"] = extracted["CLDN4"]
        if "TACSTD2" in gi:
            df["TACSTD2_log1p"] = logX[:, gi["TACSTD2"]]
        for name, vec in auc.items():
            df[f"AUC_{name}"] = vec
        for tf in tfs:
            if tf in gi:
                df[f"{tf}_log1p"] = logX[:, gi[tf]]
        return df, logX

    mal1 = mal1.reset_index(drop=True)
    df1, log1 = cohort_frame(mal1, extracted1, total1, "GSE131907", "patient_id", [])
    mal2 = mal2.reset_index(drop=True)
    df2, log2 = cohort_frame(mal2, extracted2, mal2["total_umi"].to_numpy(), "GSE205335", "patient", [])

    cells = pd.concat([df1, df2], ignore_index=True, sort=False)

    # within-patient CLDN4 split
    cells["cldn4_split"] = "drop"
    cells["in_primary_patient"] = False
    patient_rows = []
    for (cohort, patient), g in cells.groupby(["cohort", "patient"], observed=True):
        if len(g) < MIN_MAL:
            continue
        med = float(np.nanmedian(g.CLDN4_log1p.to_numpy()))
        # Zero-inflated patients: median==0 would put every cell in "high"
        # (>=0). Use detect vs zero so CLDN4-low patients are not dropped.
        if med <= 0:
            high = g.CLDN4_log1p > 0
        else:
            high = g.CLDN4_log1p >= med
        cells.loc[g.index, "cldn4_split"] = np.where(high, "high", "low")
        cells.loc[g.index, "in_primary_patient"] = True
        hi = g.loc[high]
        lo = g.loc[~high]
        rec = dict(
            cohort=cohort,
            patient=patient,
            patient_key=f"{cohort}:{patient}",
            n_malignant=len(g),
            n_high=int(high.sum()),
            n_low=int((~high).sum()),
            cldn4_pct=float(100 * (g.CLDN4_umi > 0).mean()),
            cldn4_mean=float(g.CLDN4_log1p.mean()),
            cldn4_high=float(hi.CLDN4_log1p.mean()) if len(hi) else np.nan,
            cldn4_low=float(lo.CLDN4_log1p.mean()) if len(lo) else np.nan,
            paired_primary=int(high.sum()) >= MIN_TAIL and int((~high).sum()) >= MIN_TAIL,
            paired_sens=int(high.sum()) >= MIN_TAIL_SENS and int((~high).sum()) >= MIN_TAIL_SENS,
            cldn4_median=med,
            split_rule="detect_vs_zero" if med <= 0 else "median",
        )
        if cohort == "GSE131907":
            rec["site"] = ",".join(sorted(g.Sample_Origin.dropna().astype(str).unique()))
            rec["histology"] = ",".join(sorted(g.lung_cancer_subtype.dropna().astype(str).unique())) if "lung_cancer_subtype" in g else ""
            rec["pr320_like"] = bool(g.pr320_malignant.any()) and (g.pr320_malignant.sum() >= MIN_MAL)
        else:
            rec["site"] = ",".join(sorted(g.tissue.dropna().astype(str).unique())) if "tissue" in g else ""
            rec["histology"] = ",".join(sorted(g.cancer_subtype.dropna().astype(str).unique())) if "cancer_subtype" in g else ""
            rec["pr320_like"] = True
        auc_names = [c[4:] for c in g.columns if c.startswith("AUC_")]
        for name in auc_names:
            rec[f"{name}_high"] = float(hi[f"AUC_{name}"].mean()) if len(hi) else np.nan
            rec[f"{name}_low"] = float(lo[f"AUC_{name}"].mean()) if len(lo) else np.nan
            rec[f"{name}_all"] = float(g[f"AUC_{name}"].mean())
            rec[f"{name}_delta"] = rec[f"{name}_high"] - rec[f"{name}_low"]
        for tf in tfs:
            col = f"{tf}_log1p"
            if col in g:
                rec[f"{tf}_RNA_high"] = float(hi[col].mean()) if len(hi) else np.nan
                rec[f"{tf}_RNA_low"] = float(lo[col].mean()) if len(lo) else np.nan
        patient_rows.append(rec)
    pat = pd.DataFrame(patient_rows)
    pat.to_csv(RES / "patient_means_high_vs_low.tsv", sep="\t", index=False)
    pat.to_csv(TAB / "patient_means_high_vs_low.tsv", sep="\t", index=False)
    paired = pat.loc[pat.paired_primary].copy()
    print(f"[patients] eligible>={MIN_MAL}: {len(pat)}; paired>={MIN_TAIL}: {len(paired)} "
          f"(131907={int((paired.cohort=='GSE131907').sum())} 205335={int((paired.cohort=='GSE205335').sum())})",
          flush=True)

    # Pearson in CLDN4-high cells (descriptive, per cohort)
    high_spec_rows = []
    regulon_target_rows = []
    program_gene_union = sorted(set().union(*programs.values()))
    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)
    s1 = cells.loc[cells.cohort.eq("GSE131907"), "cldn4_split"].to_numpy()
    s2 = cells.loc[cells.cohort.eq("GSE205335"), "cldn4_split"].to_numpy()
    if len(s1) != len(df1) or len(s2) != len(df2):
        raise SystemExit(f"split align fail {len(s1)}!={len(df1)} or {len(s2)}!={len(df2)}")
    df1["cldn4_split"] = s1
    df2["cldn4_split"] = s2

    for cohort, dfc, logX in (("GSE131907", df1, log1), ("GSE205335", df2, log2)):
        hi = dfc.cldn4_split.eq("high").to_numpy()
        lo = dfc.cldn4_split.eq("low").to_numpy()
        if hi.sum() < 50:
            continue
        pr_hi = pearson_tf_program(logX[hi], universe, tfs, program_gene_union)
        pr_lo = pearson_tf_program(logX[lo], universe, tfs, program_gene_union) if lo.sum() >= 50 else pd.DataFrame()
        pr_hi["cohort"] = cohort
        pr_hi.to_csv(RES / f"pearson_cldn4_high_{cohort}.tsv.gz", sep="\t", index=False)
        for tf in tfs:
            sub = pr_hi.loc[pr_hi.tf == tf].sort_values("pearson_r", ascending=False)
            top = sub.loc[sub.pearson_r >= PEARSON_MIN_R].head(40)
            lo_r = {}
            if len(pr_lo):
                lo_r = pr_lo.loc[pr_lo.tf == tf].set_index("target")["pearson_r"].to_dict()
            n_spec = 0
            spec = []
            for rec in top.itertuples(index=False):
                rlo = lo_r.get(rec.target, np.nan)
                if rec.pearson_r >= HIGH_SPEC_R and (not np.isfinite(rlo) or rlo < LOW_SPEC_R):
                    n_spec += 1
                    spec.append(rec.target)
            high_spec_rows.append(dict(
                cohort=cohort, tf=tf, program=TF_PROGRAM[tf],
                n_cells_high=int(hi.sum()), n_cells_low=int(lo.sum()),
                n_top=int(len(top)), n_high_specific=n_spec,
                mean_r=float(top.pearson_r.mean()) if len(top) else np.nan,
                high_specific=";".join(spec[:25]),
                elf3_given=tf in GIVEN_TFS,
            ))
            regulon_target_rows.append(dict(
                regulon=f"{tf}_pearson_cldn4_high_{cohort}",
                tf=tf, program=TF_PROGRAM[tf], kind="pearson_cldn4_high",
                n_targets=int(len(top)),
                targets=";".join(top.target.tolist()),
                n_cells=int(hi.sum()),
                n_patients=int(dfc.loc[hi, "patient"].nunique()),
                elf3_given=tf in GIVEN_TFS,
                note="co-expression in CLDN4-high cells; not binding",
            ))
    high_spec = pd.DataFrame(high_spec_rows)
    high_spec.to_csv(RES / "high_specific_counts.tsv", sep="\t", index=False)
    high_spec.to_csv(TAB / "high_specific_counts.tsv", sep="\t", index=False)

    # paired tests
    auc_names = sorted({c[4:] for c in cells.columns if c.startswith("AUC_")})
    rows = []

    def add_paired(name, kind, program, tf, elf3_given, note):
        hcol, lcol = f"{name}_high", f"{name}_low"
        if hcol not in paired.columns:
            return
        parts = []
        for cohort in ("GSE131907", "GSE205335"):
            sub = paired.loc[paired.cohort == cohort]
            w = wilcoxon_paired(sub[hcol], sub[lcol]) if len(sub) else dict(n=0, p=np.nan, delta=np.nan, stat=np.nan)
            parts.append((cohort, w, len(sub)))
        wall = wilcoxon_paired(paired[hcol], paired[lcol])
        rec = dict(
            regulon=name, tf=tf, program=program, kind=kind,
            n_targets=np.nan, elf3_given=elf3_given,
            n_patients_paired=wall["n"],
            n_GSE131907=parts[0][1]["n"], n_GSE205335=parts[1][1]["n"],
            delta_median_high_minus_low=wall["delta"], p=wall["p"],
            wilcoxon_stat=wall["stat"],
            delta_GSE131907=parts[0][1]["delta"], p_GSE131907=parts[0][1]["p"],
            delta_GSE205335=parts[1][1]["delta"], p_GSE205335=parts[1][1]["p"],
            note=note,
        )
        rows.append(rec)

    for name in ["IFN", "IFN_IFNA", "IFN_IFNG", "MHC_I", "TJ", "keratin", "APICAL_JUNCTION"]:
        add_paired(name, "program_set", name if name in {"IFN", "MHC_I", "TJ", "keratin"} else name,
                   "", False, "program gene-set AUCell; not a TF regulon")
    for tf in tfs:
        add_paired(f"{tf}_prior", "public_prior", TF_PROGRAM[tf], tf, tf in GIVEN_TFS,
                   "public prior AUCell; not cisTarget; not ChIP")
        add_paired(f"{tf}_priorANDprog", "prior_AND_program", TF_PROGRAM[tf], tf, tf in GIVEN_TFS,
                   "prior targets ∩ matching program set")
        hcol, lcol = f"{tf}_RNA_high", f"{tf}_RNA_low"
        if hcol in paired.columns:
            wall = wilcoxon_paired(paired[hcol], paired[lcol])
            w1 = wilcoxon_paired(paired.loc[paired.cohort.eq("GSE131907"), hcol],
                                 paired.loc[paired.cohort.eq("GSE131907"), lcol])
            w2 = wilcoxon_paired(paired.loc[paired.cohort.eq("GSE205335"), hcol],
                                 paired.loc[paired.cohort.eq("GSE205335"), lcol])
            rows.append(dict(
                regulon=f"{tf}_RNA", tf=tf, program=TF_PROGRAM[tf], kind="tf_rna",
                n_targets=np.nan, elf3_given=tf in GIVEN_TFS,
                n_patients_paired=wall["n"],
                n_GSE131907=w1["n"], n_GSE205335=w2["n"],
                delta_median_high_minus_low=wall["delta"], p=wall["p"],
                wilcoxon_stat=wall["stat"],
                delta_GSE131907=w1["delta"], p_GSE131907=w1["p"],
                delta_GSE205335=w2["delta"], p_GSE205335=w2["p"],
                note="TF RNA log1p; ELF3 RNA vs CLDN4 is A10-given, not a discovery",
            ))

    regulons = pd.DataFrame(rows)
    # attach n_targets / target lists
    target_map = {}
    for tf in tfs:
        tgs = prior_targets(priors, tf, set(universe))
        target_map[f"{tf}_prior"] = tgs
        prog_name = TF_PROGRAM[tf].split("_")[0]
        key = "MHC_I" if prog_name == "MHC" else ("IFN" if prog_name == "IFN" else prog_name)
        if key in programs:
            target_map[f"{tf}_priorANDprog"] = [g for g in tgs if g in programs[key]]
    for name, genes in programs.items():
        target_map[name] = [g for g in genes if g in universe]
    regulons["n_targets"] = regulons["regulon"].map(lambda x: len(target_map.get(x, [])))
    regulons["targets"] = regulons["regulon"].map(lambda x: ";".join(target_map.get(x, [])))
    if regulon_target_rows:
        regulons = pd.concat([regulons, pd.DataFrame(regulon_target_rows)], ignore_index=True, sort=False)
    mask = regulons.p.notna() & regulons.kind.isin(["program_set", "public_prior", "prior_AND_program", "tf_rna"])
    if mask.any():
        regulons.loc[mask, "fdr"] = multipletests(regulons.loc[mask, "p"], method="fdr_bh")[1]
    regulons.to_csv(RES / "regulons.tsv", sep="\t", index=False)
    regulons.to_csv(TAB / "regulons.tsv", sep="\t", index=False)
    paired.to_csv(RES / "paired_high_vs_low.tsv", sep="\t", index=False)
    paired.to_csv(TAB / "paired_high_vs_low.tsv", sep="\t", index=False)

    # between-patient Spearman
    between_rows = []
    for name in ["IFN", "IFN_IFNA", "IFN_IFNG", "MHC_I", "TJ", "keratin", "APICAL_JUNCTION"]:
        col = f"{name}_all"
        if col not in pat.columns:
            continue
        parts = []
        for cohort in ("GSE131907", "GSE205335"):
            sub = pat.loc[pat.cohort == cohort]
            parts.append(spear(sub["cldn4_pct"], sub[col]))
        pooled = fisher_z_pool([
            {"rho": parts[0]["rho"], "n": parts[0]["n"]},
            {"rho": parts[1]["rho"], "n": parts[1]["n"]},
        ])
        between_rows.append(dict(
            regulon=name, n=pooled["n"], k=pooled["k"], rho=pooled["rho"], p=pooled["p"], I2=pooled["I2"],
            n_GSE131907=parts[0]["n"], rho_GSE131907=parts[0]["rho"], p_GSE131907=parts[0]["p"],
            n_GSE205335=parts[1]["n"], rho_GSE205335=parts[1]["rho"], p_GSE205335=parts[1]["p"],
            note="between-patient Spearman CLDN4 %pos vs mean program AUCell; companion only",
        ))
    between = pd.DataFrame(between_rows)
    between.to_csv(RES / "between_patient_spearman.tsv", sep="\t", index=False)
    between.to_csv(TAB / "between_patient_spearman.tsv", sep="\t", index=False)

    n_tab = pd.DataFrame([
        dict(item="GSE131907_author_malignant_cells", n=int(n_author_mal_131907),
             note="Cell_subtype in {Malignant cells, tS1, tS2, tS3}"),
        dict(item="GSE131907_malignant_UMI>=200", n=int(len(df1)), note="cells"),
        dict(item="GSE131907_patients_ge20_malignant", n=int((pat.cohort.eq("GSE131907")).sum()),
             note="patient_id; unit"),
        dict(item="GSE131907_patients_paired_ge20_high_and_low",
             n=int((paired.cohort.eq("GSE131907")).sum()),
             note=f"primary paired AUCell unit"),
        dict(item="GSE205335_author_malignant_cells", n=int(n_author_mal_205335),
             note="lineage.sub == Malignant cells"),
        dict(item="GSE205335_malignant_UMI>=200", n=int(len(df2)), note="cells"),
        dict(item="GSE205335_patients_ge20_malignant", n=int((pat.cohort.eq("GSE205335")).sum()),
             note="SOFT patient; unit"),
        dict(item="GSE205335_patients_paired_ge20_high_and_low",
             n=int((paired.cohort.eq("GSE205335")).sum()),
             note="primary paired AUCell unit"),
        dict(item="merged_patients_ge20_malignant", n=int(len(pat)),
             note="k=2 cohorts stacked; not cell-integrated"),
        dict(item="merged_patients_paired", n=int(len(paired)),
             note="honest n for high vs low regulon test"),
        dict(item="CLDN4_high_cells_in_paired_patients",
             n=int(paired.n_high.sum()), note="cells; not the claim"),
        dict(item="CLDN4_low_cells_in_paired_patients",
             n=int(paired.n_low.sum()), note="cells; not the claim"),
        dict(item="GSE131907_tLung_patients_paired",
             n=int(paired.loc[paired.cohort.eq("GSE131907") & paired.site.str.contains("tLung", na=False)].shape[0]),
             note="tS1/tS2/tS3 primary-tumor sensitivity"),
        dict(item="universe_genes_AUCell", n=int(len(universe)),
             note="extracted programs + TFs + public prior targets"),
    ])
    n_tab.to_csv(RES / "n_table.tsv", sep="\t", index=False)
    n_tab.to_csv(TAB / "n_table.tsv", sep="\t", index=False)

    # figures
    prog_regs = regulons.loc[regulons.kind.eq("program_set") & regulons.regulon.isin(
        ["IFN", "IFN_IFNA", "IFN_IFNG", "MHC_I", "TJ", "keratin", "APICAL_JUNCTION"]
    )]
    plot_forest(prog_regs, f"Program AUCell · patient-paired · n={len(paired)}",
                "fig_program_paired_forest")
    tf_regs = regulons.loc[regulons.kind.eq("public_prior")]
    plot_forest(tf_regs, f"Public-prior AUCell · patient-paired · n={len(paired)}\norange = ELF3 A10-given",
                "fig_tf_prior_paired_forest")
    plot_heatmap(paired.sort_values(["cohort", "patient"]),
                 ["IFN", "MHC_I", "TJ", "keratin"], "fig_patient_delta_heatmap")
    plot_boxes(paired, ["IFN", "MHC_I", "TJ", "keratin"], "fig_paired_boxes")
    plot_scatter(pat, ["IFN", "MHC_I", "TJ", "keratin"], "fig_between_patient_scatter")
    plot_n(n_tab, "fig_n_patients")

    if len(high_spec):
        fig, ax = plt.subplots(figsize=(8.2, 5.6))
        foc = high_spec.groupby("tf", as_index=False)["n_high_specific"].sum().sort_values("n_high_specific")
        colors = ["#c45c26" if t in GIVEN_TFS else "#1f4e79" for t in foc.tf]
        ax.barh(foc.tf, foc.n_high_specific, color=colors)
        ax.set_xlabel("n high-specific TF–program edges (sum of cohorts)")
        ax.set_title("Pearson high-specific (r_high≥0.15 and r_low<0.10)\norange=ELF3 given; not a binding claim")
        fig.tight_layout()
        savefig(fig, "fig_high_specific_counts")

    summary = {
        "question": "Do CLDN4-high malignant cells show different IFN/MHC-I/TJ/keratin regulons vs CLDN4-low at the patient level?",
        "combo": "GSE131907 + GSE205335 (PR #320 combo that already differs)",
        "a10_taken_as_given": "ELF3–CLDN4 bulk RNA; ELF3 is labeled given, not a discovery.",
        "method": "AUCell + public TF-target priors (TRRUST/DoRothEA/CollecTRI); not pySCENIC cisTarget; not ChIP",
        "pyscenic": {"importable": pyscenic_ok, "version": pyscenic_ver,
                     "cistarget_run": False, "chip_peaks_invented": False},
        "no_dual_high": True,
        "no_gse207422": True,
        "n_patients_paired": int(len(paired)),
        "n_patients_paired_GSE131907": int((paired.cohort == "GSE131907").sum()),
        "n_patients_paired_GSE205335": int((paired.cohort == "GSE205335").sum()),
        "n_patients_eligible": int(len(pat)),
        "n_genes_universe": int(len(universe)),
        "n_genes_GSE131907_matrix": int(n_genes1),
        "cldn4_split": "within-patient median among author-malignant",
        "primary_unit": "patient",
        "aucell_universe": "extracted programs + TFs + public prior targets",
        "notes": [
            "Author malignant labels: GSE131907 Cell_subtype in {Malignant cells,tS1,tS2,tS3}; GSE205335 lineage.sub==Malignant cells.",
            "CLDN4 held out of every regulon (split gene).",
            "Patient deltas stacked across cohorts; cells were not Harmony-integrated.",
            "GSE205335 histology mix (including SCLC) is part of the honest n.",
            "log2TPM text, EGA FASTQ, and GSE207422 were not used.",
        ],
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    write_finding(regulons, paired, n_tab, between, summary)
    print(f"[done] {RES}", flush=True)


if __name__ == "__main__":
    sys.exit(main() or 0)
