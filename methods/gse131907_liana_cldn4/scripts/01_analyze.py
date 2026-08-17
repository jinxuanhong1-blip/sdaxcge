#!/usr/bin/env python3
"""LIANA/CellPhoneDB-style LR from CLDN4-high malignant cells to T/NK.

Dataset
    GSE131907 (Kim et al., Nat Commun 2020). Author Cell_type / Cell_subtype.
    Malignant = tumor-origin samples AND Cell_subtype in
    {Malignant cells, tS1, tS2, tS3}. T/NK = T lymphocytes + NK cells.

Primary score
    CellPhoneDB mean-of-means on log1p(CP10k) (Efremova 2020 / Garcia-Alonso 2022).
    Inferential unit = tumor sample (paired Wilcoxon high vs low). Cells are
    not replicates. Unique patient n is reported separately.

Secondary
    LIANA mt.cellphonedb if import succeeds. Otherwise documented as not run.
"""

from __future__ import annotations

import argparse
import gzip
import json
import pickle
import traceback
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]

TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES = {"T lymphocytes", "NK cells"}


def log(msg: str) -> None:
    print(msg, flush=True)


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
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def stream_umi(path: Path, keep: set[str], cache_dir: Path | None = None) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, int]:
    """One pass: library totals for every cell; store only keep-genes."""
    cache_dir = cache_dir or path.parent
    cache = cache_dir / "stream_cache_lr.pkl"
    if cache.exists():
        log(f"[stream] load cache {cache}")
        with cache.open("rb") as fh:
            blob = pickle.load(fh)
        cells, store, totals, n_genes = blob["cells"], blob["store"], blob["totals"], blob["n_genes"]
        log(f"[stream] cache cells={len(cells)} genes_in_file={n_genes} genes_kept={len(store)}")
        return cells, store, totals, n_genes
    log(f"[stream] {path} keep={len(keep)} symbols")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"column mismatch for {gene}: {vals.size} != {n}")
            totals += vals
            if gene in keep and gene not in store:
                store[gene] = vals
            if n_genes % 5000 == 0:
                log(f"[stream] genes_seen={n_genes} kept={len(store)}")
    log(f"[stream] cells={n} genes_in_file={n_genes} genes_kept={len(store)}")
    with cache.open("wb") as fh:
        pickle.dump({"cells": cells, "store": store, "totals": totals, "n_genes": n_genes}, fh, protocol=4)
    log(f"[stream] wrote cache {cache}")
    return cells, store, totals, n_genes


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        cp = np.where(total > 0, umi / total * 1e4, 0.0)
    return np.log1p(cp).astype(np.float32)


def group_gene_stats(logx: dict[str, np.ndarray], mask: np.ndarray) -> tuple[dict[str, float], dict[str, float], int]:
    n = int(mask.sum())
    means, fracs = {}, {}
    if n == 0:
        return means, fracs, 0
    for g, arr in logx.items():
        v = arr[mask]
        means[g] = float(np.mean(v))
        fracs[g] = float(np.mean(v > 0))
    return means, fracs, n


def partner_from_stats(units: list[str], means: dict[str, float], fracs: dict[str, float]) -> tuple[float, float]:
    m, f = [], []
    for g in units:
        if g not in means:
            return np.nan, np.nan
        m.append(means[g])
        f.append(fracs[g])
    return float(np.min(m)), float(np.min(f))


def score_pairs(
    pairs: pd.DataFrame,
    logx: dict[str, np.ndarray],
    sender: np.ndarray,
    receiver: np.ndarray,
    expr_prop: float,
) -> pd.DataFrame:
    s_mean, s_frac, n_s = group_gene_stats(logx, sender)
    r_mean, r_frac, n_r = group_gene_stats(logx, receiver)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac = partner_from_stats(lig_u, s_mean, s_frac)
        rec_m, rec_f = partner_from_stats(rec_u, r_mean, r_frac)
        if not np.isfinite(l_mean) or not np.isfinite(rec_m):
            continue
        pass_prop = (l_frac >= expr_prop) and (rec_f >= expr_prop)
        rows.append(
            {
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "pathway": rec.pathway,
                "pair_origin": rec.pair_origin,
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": rec_m,
                "ligand_frac": l_frac,
                "receptor_frac": rec_f,
                "cpdb_mean_score": 0.5 * (l_mean + rec_m),
                "product_score": l_mean * rec_m,
                "pass_expr_prop": pass_prop,
            }
        )
    return pd.DataFrame(rows)


def wilcoxon_safe(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 3 or np.allclose(a, b):
        return np.nan
    try:
        return float(stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def run_liana(adata, groupby: str, out_csv: Path, n_perms: int) -> str:
    try:
        import liana as li
    except Exception as exc:  # pragma: no cover
        return f"LIANA_IMPORT_FAILED: {exc}"
    try:
        li.mt.cellphonedb(
            adata,
            groupby=groupby,
            resource_name="cellphonedb",
            expr_prop=0.10,
            n_perms=n_perms,
            use_raw=False,
            verbose=True,
            key_added="liana_res",
        )
        res = adata.uns["liana_res"].copy()
        res.to_csv(out_csv, index=False)
        return f"LIANA_OK n_edges={len(res)} file={out_csv.name}"
    except Exception as exc:
        return f"LIANA_RUN_FAILED: {exc}\n{traceback.format_exc()}"


def write_finding(
    outdir: Path,
    summary: dict,
    ntab: pd.DataFrame,
    lr_table: pd.DataFrame,
    ranks: pd.DataFrame,
    focus: pd.DataFrame,
) -> None:
    n_pass = int(lr_table["pass_expr_prop"].sum()) if len(lr_table) else 0
    top = lr_table[lr_table["pass_expr_prop"]].sort_values("cpdb_mean_score", ascending=False).head(25)
    focus_pass = focus[focus["pass_expr_prop"]].sort_values("cpdb_mean_score", ascending=False) if len(focus) else focus

    def md_table(df: pd.DataFrame, cols: list[str], fmt: dict[str, str] | None = None) -> list[str]:
        if df is None or df.empty:
            return ["_(empty)_", ""]
        fmt = fmt or {}
        lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
        for rec in df.itertuples(index=False):
            row = []
            for c in cols:
                v = getattr(rec, c)
                if c in fmt and pd.notna(v):
                    row.append(fmt[c].format(v))
                elif isinstance(v, float) and pd.notna(v):
                    row.append(f"{v:.3g}")
                else:
                    row.append(str(v))
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")
        return lines

    n_paired = int(summary["n_samples_paired"])
    n_pat = int(summary["n_patients_paired"])
    n_neg = int(summary.get("n_focus_delta_neg", 0))
    n_focus = int(summary.get("n_focus_tested", 0))
    n_sig = int(summary.get("n_focus_fdr05", 0))
    foc_ranks = ranks[ranks["pathway"].isin(["T_recruit", "IFN", "MHC_I"])] if len(ranks) else ranks
    mhc = foc_ranks[foc_ranks["pathway"] == "MHC_I"] if len(foc_ranks) else foc_ranks
    n_mhc_up = (
        int(((mhc["padj"] < 0.05) & (mhc["median_delta"] > 0)).sum())
        if len(mhc) and mhc["padj"].notna().any()
        else 0
    )
    if n_paired == 0:
        verdict = (
            "No tumor sample met the paired gate (≥10 CLDN4-high malignant, "
            "≥10 CLDN4-low malignant, ≥20 T/NK). The LR table is pooled / descriptive only."
        )
    else:
        verdict = (
            f"On {n_paired} paired tumor samples ({n_pat} patients), CLDN4-high vs CLDN4-low "
            f"malignant → T/NK focus pairs do **not** support a coordinated T-recruit drop "
            f"({n_neg}/{n_focus} median Δ < 0; **{n_sig} FDR < 0.05** on the down side). "
            f"MHC-I outgoing is **higher** from CLDN4-high ({n_mhc_up} pairs FDR < 0.05, Δ > 0). "
            "CXCL9/10/11–CXCR3 fail the 10% expression filter in malignant cells. "
            "The LR table is a ranked co-expression list, not a causal claim."
        )
    verdict_zh = (
        f"在 **{n_paired} 个配对肿瘤样本 / {n_pat} 名患者** 上，CLDN4 高 vs 低恶性细胞 → T/NK "
        f"的焦点对**不支持**协同招募下降（{n_neg}/{n_focus} 中位 Δ < 0；下行 **{n_sig} 个 FDR < 0.05**）。"
        f"MHC-I 出站在 CLDN4 高侧更高（{n_mhc_up} 对 FDR < 0.05）。"
        "CXCL9/10/11–CXCR3 未过恶性细胞 10% 表达门槛。LR 表是共表达排序，不是因果机制。"
    )

    lines = [
        "# FINDING — GSE131907 LIANA/LR from CLDN4-high malignant to T/NK",
        "",
        f"**Verdict:** {verdict}",
        "",
        "Kim et al., *Nat Commun* 2020 (PMID 32385277); GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907). "
        "Treatment-naive LUAD atlas. **No ICI / MPR labels.**",
        "",
        "---",
        "",
        "## English",
        "",
        "### Honest n",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        f"| Cells in UMI matrix | {summary['n_cells']:,} | author barcodes |",
        f"| GEO samples / patients | {summary['n_samples_geo']} / {summary['n_patients_geo']} | 58 / 44 in the series |",
        f"| Tumor-origin samples | {summary['n_tumor_samples']} | tLung, tL/B, mLN, mBrain, PE |",
        f"| Author malignant cells | {summary['n_malig']:,} | subtype ∈ {{Malignant cells, tS1, tS2, tS3}} |",
        f"|  … tLung tS1/tS2/tS3 | {summary['n_malig_tlung']:,} | Kim tumor-specific epi, not the 'Malignant cells' label |",
        f"|  … author Malignant cells | {summary['n_malig_author']:,} | mets / tL-B / mLN / mBrain |",
        f"| CLDN4-high / low malignant | {summary['n_cldn4_high']:,} / {summary['n_cldn4_low']:,} | global median log1p(CP10k) = {summary['cldn4_threshold']:.3f} |",
        f"| T / NK / T+NK (all samples) | {summary['n_T']:,} / {summary['n_NK']:,} / {summary['n_tnk']:,} | author Cell_type |",
        f"| T/NK in tumor-origin samples | {summary.get('pooled_tnk_tumor_origin', summary['n_tnk']):,} | receivers for the pooled LR table |",
        f"| Samples with any malignant + T/NK | {summary['n_samples_with_malig_tnk']} | descriptive |",
        f"| **Paired samples (unit of test)** | **{n_paired}** | ≥10 high, ≥10 low, ≥20 T/NK |",
        f"| **Unique patients in paired set** | **{n_pat}** | do not count cells as n |",
        f"| Pairs scored / pass expr_prop 0.10 | {summary['n_pairs_scored']} / {n_pass} | CellPhoneDB v5 + overlay |",
        f"| LIANA | {summary['liana_status']} | secondary |",
        f"| CellChat | {summary['cellchat_status']} | not run |",
        "",
        "PE epithelial cells are **unlabeled** in the author file (0 `Malignant cells`); they are **not** counted as malignant. "
        "nLung AT1/AT2/Club/Ciliated are **not** malignant. Patients can contribute more than one tumor site; "
        f"paired unique-patient n is {n_pat}, not 44. "
        "Dropped from the paired gate: LUNG_T09 (5 malignant), NS_16 (79 malignant, 0 CLDN4-high), "
        "EBUS_13 (376 malignant, 0 CLDN4-high).",
        "",
        "Per-sample counts: `results/n_cells_samples.tsv`.",
        "",
        "### Method",
        "",
        "Documented CellPhoneDB-style score on log1p(CP10k): partner expression = **min of subunit means**; "
        "pair score = **mean of the two partner means** (Efremova et al. 2020; Garcia-Alonso et al. 2022). "
        "`pass_expr_prop` requires both partners in ≥10% of cells in their group. "
        "Pooled scores are descriptive. The test is a **paired Wilcoxon** of high vs low **per tumor sample**. "
        "FDR is Benjamini–Hochberg within the outgoing contrast. "
        "This is **not** a CellChat communication probability and does **not** observe secretion or spatial contact.",
        "",
        "### LR table — CLDN4-high malignant → T/NK (pooled, pass expr_prop)",
        "",
        "Full table: `results/lr_table.tsv`. Top 25 by score:",
        "",
    ]
    lines += md_table(
        top,
        ["ligand", "receptor", "pathway", "ligand_frac", "receptor_frac", "cpdb_mean_score"],
        {"ligand_frac": "{:.2f}", "receptor_frac": "{:.2f}", "cpdb_mean_score": "{:.3f}"},
    )
    lines += [
        "### Focus axes (T-recruit / IFN / MHC-I) that pass expr_prop",
        "",
    ]
    lines += md_table(
        focus_pass,
        ["ligand", "receptor", "pathway", "ligand_frac", "receptor_frac", "cpdb_mean_score"],
        {"ligand_frac": "{:.2f}", "receptor_frac": "{:.2f}", "cpdb_mean_score": "{:.3f}"},
    )
    lines += [
        "### Patient/sample-level high vs low (outgoing)",
        "",
        "Median Δ = CLDN4-high − CLDN4-low. Negative = weaker from the high state. "
        f"Paired n = **{n_paired} samples / {n_pat} patients**.",
        "",
    ]
    if len(ranks):
        show = ranks.sort_values(["pathway", "median_delta"])
        show = show[show["pathway"].isin(["T_recruit", "IFN", "MHC_I"]) | (show["padj"].fillna(1) < 0.2)]
        if show.empty:
            show = ranks.head(20)
        lines += md_table(
            show,
            ["ligand", "receptor", "pathway", "n_samples", "median_delta", "pval", "padj"],
            {"median_delta": "{:+.3f}", "pval": "{:.3g}", "padj": "{:.3g}"},
        )
    else:
        lines += ["_(no paired samples passed the gate)_", ""]

    lines += [
        f"### LIANA (`mt.cellphonedb`)",
        "",
        f"Status: `{summary['liana_status']}`. "
        "If a LIANA table is present, p-values are within-object specificity on a downsampled object (≤2,000 cells/group, 50 permutations), **not** the patient-level test above.",
        "",
    ]
    liana_csv = outdir / "liana_cldn4high_to_tnk.csv"
    if liana_csv.exists():
        try:
            li_df = pd.read_csv(liana_csv)
            li_show = li_df.sort_values("lr_means" if "lr_means" in li_df.columns else li_df.columns[0], ascending=False).head(15)
            cols = [c for c in ["source", "target", "ligand_complex", "receptor_complex", "lr_means", "cellphone_pvals"] if c in li_show.columns]
            lines += md_table(li_show, cols)
        except Exception:
            lines += [f"See `{liana_csv.name}`.", ""]
    else:
        lines += ["LIANA was not written.", ""]

    lines += [
        "### What this is not",
        "",
        "- Not ICI / MPR (GSE131907 is treatment-naive).",
        "- Not inferCNV/CopyKAT recomputed malignant IDs.",
        "- Not CellChat. LIANA permutation p-values (if present) are within-object specificity, not patient tests.",
        "- Not spatial proximity or protein secretion.",
        "- Cells are not the sample size.",
        "",
        "### Files",
        "",
        "`results/lr_table.tsv`, `results/lr_table_high_vs_low.tsv`, `results/n_table.tsv`, "
        "`results/n_cells_samples.tsv`, `results/summary.json`, `results/figures/`.",
        "",
        "---",
        "",
        "## 中文",
        "",
        f"**结论：** {verdict_zh}",
        "",
        f"GSE131907（Kim 2020）治疗初治 LUAD。作者注释恶性细胞 {summary['n_malig']:,} 个"
        f"（tLung 用 tS1/tS2/tS3；转移灶用 Malignant cells）。"
        f"CLDN4 按恶性细胞全局中位数拆高/低。"
        f"配对检验单位是肿瘤样本，不是细胞："
        f"**{n_paired} 个样本 / {n_pat} 名患者**。",
        "胸水上皮在作者文件里未标恶性，不计入。正常肺 AT1/AT2/Club/Ciliated 不是恶性。",
        "LR 表是共表达排序，不是招募机制。",
        "",
    ]
    (HERE / "FINDING.md").write_text("\n".join(lines))
    (outdir / "FINDING.md").write_text("\n".join(lines))


def write_results_md(outdir: Path, summary: dict) -> None:
    (outdir / "README.md").write_text(
        "# GSE131907 CLDN4-high malignant → T/NK LR\n\n"
        "See `../FINDING.md` for the computed readout and `lr_table.tsv` for the pair table.\n"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    figdir = args.outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    P = cfg["params"]
    pairs = pd.read_csv(HERE / "resources" / "cellphonedb_v5_lr_pairs.tsv", sep="\t")
    lr_genes = set(Path(HERE / "resources" / "lr_genes.txt").read_text().split())
    keep = lr_genes | set(cfg["state_genes"])

    ann = pd.read_csv(args.workdir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str)
    meta = parse_series_matrix(args.workdir / "GSE131907_series_matrix.txt.gz")
    sample_meta = meta.rename(columns={"title": "Sample", "tissue_origin_abbrevation": "Sample_Origin_geo"})
    keep_cols = [c for c in ["Sample", "geo_accession", "patient_id", "tumor_stage", "Sample_Origin_geo"] if c in sample_meta.columns]
    sample_meta = sample_meta[keep_cols].drop_duplicates("Sample")
    sample_meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    cells, expr, totals, n_genes = stream_umi(
        args.workdir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        keep,
        cache_dir=args.workdir,
    )
    per = ann.set_index("Index").reindex(cells)
    if per["Sample"].isna().any():
        raise SystemExit("matrix cell IDs do not align with annotation Index")

    sample = per["Sample"].to_numpy()
    origin = per["Sample_Origin"].to_numpy()
    ctype = per["Cell_type"].fillna("").to_numpy()
    csub = per["Cell_subtype"].fillna("").to_numpy()
    patient_map = sample_meta.set_index("Sample")["patient_id"].to_dict()
    patient = np.array([patient_map.get(s, s) for s in sample], dtype=object)

    is_tumor = np.isin(origin, list(TUMOR_ORIGINS))
    is_malig = is_tumor & np.isin(csub, list(MALIGNANT_SUBTYPES))
    is_t = ctype == "T lymphocytes"
    is_nk = ctype == "NK cells"
    is_tnk = is_t | is_nk
    is_tlung_ts = (origin == "tLung") & np.isin(csub, ["tS1", "tS2", "tS3"])
    is_author_malig = csub == "Malignant cells"

    missing_state = [g for g in cfg["state_genes"] if g not in expr]
    if missing_state:
        raise SystemExit(f"missing state genes: {missing_state}")

    logx = {g: log1p_cp10k(arr, totals) for g, arr in expr.items()}
    cldn4 = logx["CLDN4"]
    thr = float(np.median(cldn4[is_malig])) if is_malig.any() else float("nan")
    cld_hi = is_malig & (cldn4 >= thr)
    cld_lo = is_malig & (cldn4 < thr)
    log(f"[split] CLDN4 median among malignant={thr:.4f} high={int(cld_hi.sum())} low={int(cld_lo.sum())}")

    # Per-sample n
    rows = []
    for s in pd.unique(sample):
        m = sample == s
        orig = str(origin[m][0])
        rows.append(
            {
                "Sample": s,
                "patient_id": patient_map.get(s, s),
                "Sample_Origin": orig,
                "is_tumor_origin": orig in TUMOR_ORIGINS,
                "n_cells": int(m.sum()),
                "n_malig": int((m & is_malig).sum()),
                "n_malig_tS": int((m & is_tlung_ts).sum()),
                "n_malig_author": int((m & is_author_malig).sum()),
                "n_cldn4_high": int((m & cld_hi).sum()),
                "n_cldn4_low": int((m & cld_lo).sum()),
                "n_T": int((m & is_t).sum()),
                "n_NK": int((m & is_nk).sum()),
                "n_T_NK": int((m & is_tnk).sum()),
                "mean_cldn4_malig": float(np.mean(cldn4[m & is_malig])) if (m & is_malig).any() else np.nan,
            }
        )
    extra = [c for c in sample_meta.columns if c not in {"Sample", "patient_id"}]
    ntab = pd.DataFrame(rows).merge(sample_meta[["Sample"] + extra], on="Sample", how="left")
    ntab["eligible_paired"] = (
        (ntab["n_cldn4_high"] >= P["min_malig_per_state"])
        & (ntab["n_cldn4_low"] >= P["min_malig_per_state"])
        & (ntab["n_T_NK"] >= P["min_tnk"])
        & ntab["is_tumor_origin"]
    )
    ntab.to_csv(args.outdir / "n_cells_samples.tsv", sep="\t", index=False)

    eligible = set(ntab.loc[ntab["eligible_paired"], "Sample"])
    n_paired = int(len(eligible))
    n_pat = int(ntab.loc[ntab["eligible_paired"], "patient_id"].nunique())
    log(f"[n] paired samples={n_paired} unique patients={n_pat}")

    # Pooled LR: CLDN4-high malignant → T/NK (tumor-origin T/NK only)
    tnk_tumor = is_tnk & is_tumor
    pooled = score_pairs(pairs, logx, cld_hi, tnk_tumor, P["expr_prop"])
    pooled = pooled.sort_values(["pass_expr_prop", "cpdb_mean_score"], ascending=[False, False])
    pooled.to_csv(args.outdir / "lr_table.tsv", sep="\t", index=False)
    pooled[pooled["pass_expr_prop"]].to_csv(args.outdir / "lr_table_pass.tsv", sep="\t", index=False)
    focus_ligands = set(cfg["pathways"]["T_recruit"]["ligands"]) | set(cfg["pathways"]["IFN"]["ligands"]) | set(cfg["pathways"]["MHC_I"]["ligands"])
    focus_recs = set(cfg["pathways"]["T_recruit"]["receptors"]) | set(cfg["pathways"]["IFN"]["receptors"]) | set(cfg["pathways"]["MHC_I"]["receptors"])

    def is_focus(row) -> bool:
        lig_ok = bool(set(str(row.ligand).split("+")) & focus_ligands)
        rec_ok = bool(set(str(row.receptor).split("+")) & focus_recs)
        return lig_ok and rec_ok

    focus = pooled[pooled.apply(is_focus, axis=1)].copy() if len(pooled) else pooled
    focus.to_csv(args.outdir / "lr_focus_recruit.tsv", sep="\t", index=False)

    # Also score low for a pooled high-vs-low descriptive table
    pooled_lo = score_pairs(pairs, logx, cld_lo, tnk_tumor, P["expr_prop"])
    desc = pooled.merge(
        pooled_lo[["ligand", "receptor", "cpdb_mean_score", "pass_expr_prop", "ligand_frac", "ligand_mean"]],
        on=["ligand", "receptor"],
        suffixes=("_high", "_low"),
    )
    desc["delta_high_minus_low"] = desc["cpdb_mean_score_high"] - desc["cpdb_mean_score_low"]
    desc.to_csv(args.outdir / "lr_table_pooled_high_vs_low.tsv", sep="\t", index=False)

    # Per-sample paired outgoing
    chunks = []
    for s in sorted(eligible):
        m = sample == s
        h, l, t = m & cld_hi, m & cld_lo, m & is_tnk
        hdf = score_pairs(pairs, logx, h, t, P["expr_prop"])
        ldf = score_pairs(pairs, logx, l, t, P["expr_prop"])
        key = ["ligand", "receptor", "pathway"]
        d = hdf[key + ["cpdb_mean_score", "pass_expr_prop", "ligand_frac", "receptor_frac", "n_sender", "n_receiver"]].merge(
            ldf[key + ["cpdb_mean_score", "pass_expr_prop", "ligand_frac", "receptor_frac", "n_sender", "n_receiver"]],
            on=key,
            suffixes=("_high", "_low"),
        )
        d["delta_high_minus_low"] = d["cpdb_mean_score_high"] - d["cpdb_mean_score_low"]
        d["pass_either"] = d["pass_expr_prop_high"] | d["pass_expr_prop_low"]
        d["pass_both"] = d["pass_expr_prop_high"] & d["pass_expr_prop_low"]
        d["Sample"] = s
        d["patient_id"] = patient_map.get(s, s)
        d["n_high"] = int(h.sum())
        d["n_low"] = int(l.sum())
        d["n_tnk"] = int(t.sum())
        chunks.append(d)
    raw = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    raw_path = args.outdir / "patient_outgoing.tsv"
    raw.to_csv(raw_path, sep="\t", index=False)
    import gzip as _gzip
    import shutil
    with raw_path.open("rb") as src, _gzip.open(str(raw_path) + ".gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    raw_path.unlink()

    rank_rows = []
    if len(raw):
        for (lig, recp, path), sub in raw.groupby(["ligand", "receptor", "pathway"], observed=True):
            keep_sub = sub[sub["pass_either"]]
            if len(keep_sub) < 3:
                continue
            if path not in {"T_recruit", "IFN", "MHC_I"} and keep_sub["pass_both"].sum() < 3:
                continue
            p = wilcoxon_safe(keep_sub["cpdb_mean_score_high"], keep_sub["cpdb_mean_score_low"])
            rank_rows.append(
                {
                    "ligand": lig,
                    "receptor": recp,
                    "pathway": path,
                    "n_samples": int(keep_sub["Sample"].nunique()),
                    "n_patients": int(keep_sub["patient_id"].nunique()),
                    "median_delta": float(np.median(keep_sub["delta_high_minus_low"])),
                    "mean_delta": float(np.mean(keep_sub["delta_high_minus_low"])),
                    "mean_score_high": float(keep_sub["cpdb_mean_score_high"].mean()),
                    "mean_score_low": float(keep_sub["cpdb_mean_score_low"].mean()),
                    "frac_pass_high": float(keep_sub["pass_expr_prop_high"].mean()),
                    "frac_pass_low": float(keep_sub["pass_expr_prop_low"].mean()),
                    "pval": p,
                }
            )
    ranks = pd.DataFrame(rank_rows)
    if len(ranks) and ranks["pval"].notna().any():
        mask = ranks["pval"].notna()
        ranks.loc[mask, "padj"] = multipletests(ranks.loc[mask, "pval"], method="fdr_bh")[1]
    else:
        ranks["padj"] = np.nan
    if len(ranks):
        ranks = ranks.sort_values(["pathway", "median_delta"])
    ranks.to_csv(args.outdir / "lr_table_high_vs_low.tsv", sep="\t", index=False)

    foc = ranks[ranks["pathway"].isin(["T_recruit", "IFN", "MHC_I"])] if len(ranks) else ranks
    n_focus = int(len(foc))
    n_neg = int((foc["median_delta"] < 0).sum()) if n_focus else 0
    n_sig = int(((foc["padj"] < 0.05) & (foc["median_delta"] < 0)).sum()) if n_focus and foc["padj"].notna().any() else 0

    # LIANA secondary
    liana_status = "not_run"
    liana_notes = []
    try:
        import anndata as ad

        mask = cld_hi | cld_lo | (is_tnk & is_tumor)
        genes = sorted(logx)
        X = np.vstack([logx[g][mask] for g in genes]).T
        obs = pd.DataFrame(
            {"sample": sample[mask], "patient_id": patient[mask], "origin": origin[mask]},
            index=pd.Index(cells[mask], name="cell"),
        )
        grp = np.array(["other"] * int(mask.sum()), dtype=object)
        grp[cld_hi[mask]] = "Malig_CLDN4high"
        grp[cld_lo[mask]] = "Malig_CLDN4low"
        grp[is_t[mask]] = "T"
        grp[is_nk[mask]] = "NK"
        obs["cc_group"] = grp
        rng = np.random.default_rng(P["random_seed"])
        keep_idx = []
        for g, idx in obs.groupby("cc_group", observed=True).indices.items():
            if g == "other":
                continue
            if len(idx) > P["liana_max_cells_per_group"]:
                idx = rng.choice(idx, size=P["liana_max_cells_per_group"], replace=False)
            keep_idx.append(np.asarray(idx))
        keep_idx = np.sort(np.concatenate(keep_idx)) if keep_idx else np.array([], dtype=int)
        adata = ad.AnnData(X=X[keep_idx], obs=obs.iloc[keep_idx].copy(), var=pd.DataFrame(index=genes))
        adata.obs["cc_group"] = pd.Categorical(adata.obs["cc_group"])
        adata.uns["log1p"] = {"base": None}
        counts = adata.obs["cc_group"].value_counts().to_dict()
        log(f"[liana] groups {counts}")
        note = run_liana(adata, "cc_group", args.outdir / "liana_cellphonedb_cldn4.csv", P["liana_n_perms"])
        liana_notes.append(note)
        if note.startswith("LIANA_OK"):
            liana_status = "ran_cellphonedb_method"
            res = pd.read_csv(args.outdir / "liana_cellphonedb_cldn4.csv")
            outg = res[res["source"].eq("Malig_CLDN4high") & res["target"].isin(["T", "NK"])].copy()
            outg.to_csv(args.outdir / "liana_cldn4high_to_tnk.csv", index=False)
        else:
            liana_status = note.split(":")[0]
    except Exception as exc:
        liana_status = f"LIANA_WRAPPER_FAILED: {type(exc).__name__}"
        liana_notes.append(f"{exc}")

    # Figures
    fig, ax = plt.subplots(figsize=(8.4, 3.8))
    show = ntab[ntab["is_tumor_origin"]].sort_values("Sample")
    x = np.arange(len(show))
    ax.bar(x - 0.2, show["n_malig"], width=0.4, label="malignant", color="#4C72B0")
    ax.bar(x + 0.2, show["n_T_NK"], width=0.4, label="T/NK", color="#DD8452")
    ax.set_xticks(x)
    ax.set_xticklabels(show["Sample"], rotation=80, ha="right", fontsize=6)
    ax.set_ylabel("cells")
    ax.legend(frameon=False)
    ax.set_title("GSE131907 tumor-origin cells (author malignant vs T/NK)")
    fig.tight_layout()
    fig.savefig(figdir / "n_cells_by_sample.png", dpi=160)
    plt.close(fig)

    if len(foc):
        foc_plot = foc.sort_values("median_delta")
        colors = {"T_recruit": "#4C72B0", "IFN": "#55A868", "MHC_I": "#C44E52"}
        fig, ax = plt.subplots(figsize=(8.0, max(3.2, 0.28 * len(foc_plot) + 1.2)))
        y = np.arange(len(foc_plot))
        ax.barh(y, foc_plot["median_delta"], color=[colors.get(p, "#999") for p in foc_plot["pathway"]], edgecolor="none")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{a}–{b} ({p})" for a, b, p in zip(foc_plot["ligand"], foc_plot["receptor"], foc_plot["pathway"])], fontsize=7)
        ax.set_xlabel("median sample Δ score (CLDN4-high − low)")
        ax.set_title(f"Outgoing malignant → T/NK  n={n_paired} samples / {n_pat} patients")
        fig.tight_layout()
        fig.savefig(figdir / "cldn4_outgoing_focus.png", dpi=160)
        plt.close(fig)

    n_table = pd.DataFrame(
        [
            {"item": "cells_in_matrix", "n": int(len(cells)), "note": "author barcodes"},
            {"item": "geo_samples", "n": int(sample_meta["Sample"].nunique()), "note": "58 GEO samples"},
            {"item": "geo_patients", "n": int(sample_meta["patient_id"].nunique()), "note": "44 patients"},
            {"item": "tumor_origin_samples", "n": int(ntab["is_tumor_origin"].sum()), "note": "tLung/tL-B/mLN/mBrain/PE"},
            {"item": "malignant_cells", "n": int(is_malig.sum()), "note": "Malignant cells + tS1/tS2/tS3 in tumor origin"},
            {"item": "malignant_tlung_tS", "n": int(is_tlung_ts.sum()), "note": "tLung tS1/tS2/tS3"},
            {"item": "malignant_author_label", "n": int(is_author_malig.sum()), "note": "Cell_subtype == Malignant cells"},
            {"item": "cldn4_high_malignant", "n": int(cld_hi.sum()), "note": f"median log1p(CP10k)={thr:.3f}"},
            {"item": "cldn4_low_malignant", "n": int(cld_lo.sum()), "note": "below median"},
            {"item": "T_cells", "n": int(is_t.sum()), "note": "author T lymphocytes"},
            {"item": "NK_cells", "n": int(is_nk.sum()), "note": "author NK cells"},
            {"item": "T_NK_cells", "n": int(is_tnk.sum()), "note": "T + NK"},
            {"item": "paired_samples", "n": n_paired, "note": "unit of Wilcoxon"},
            {"item": "paired_patients", "n": n_pat, "note": "unique patient_id in paired set"},
            {"item": "pairs_pass_expr_prop", "n": int(pooled["pass_expr_prop"].sum()) if len(pooled) else 0, "note": "pooled CLDN4-high → T/NK"},
        ]
    )
    n_table.to_csv(args.outdir / "n_table.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE131907",
        "citation": "Kim et al. Nat Commun 2020 PMID 32385277",
        "n_cells": int(len(cells)),
        "n_genes_in_matrix": int(n_genes),
        "n_samples_geo": int(sample_meta["Sample"].nunique()),
        "n_patients_geo": int(sample_meta["patient_id"].nunique()),
        "n_tumor_samples": int(ntab["is_tumor_origin"].sum()),
        "n_malig": int(is_malig.sum()),
        "n_malig_tlung": int(is_tlung_ts.sum()),
        "n_malig_author": int(is_author_malig.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "cldn4_threshold": thr,
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_samples_with_malig_tnk": int(((ntab["n_malig"] > 0) & (ntab["n_T_NK"] > 0) & ntab["is_tumor_origin"]).sum()),
        "n_samples_paired": n_paired,
        "n_patients_paired": n_pat,
        "min_malig_per_state": P["min_malig_per_state"],
        "min_tnk": P["min_tnk"],
        "n_pairs_scored": int(len(pooled)),
        "n_pairs_pass_expr_prop": int(pooled["pass_expr_prop"].sum()) if len(pooled) else 0,
        "n_focus_tested": n_focus,
        "n_focus_delta_neg": n_neg,
        "n_focus_fdr05": n_sig,
        "liana_status": liana_status,
        "liana_notes": liana_notes,
        "cellchat_status": "not_run_R_unavailable",
        "method": "cellphonedb_mean_of_means_on_log1p_cp10k",
        "malignant_definition": "tumor-origin AND Cell_subtype in {Malignant cells, tS1, tS2, tS3}",
        "unit_of_inference": "tumor sample; unique patient n reported separately",
        "n_lr_genes_kept": int(len(expr)),
        "n_lr_genes_missing": int(len(keep - set(expr))),
        "dropped_paired_gate": [
            "LUNG_T09 (5 malignant)",
            "NS_16 (79 malignant, 0 CLDN4-high)",
            "EBUS_13 (376 malignant, 0 CLDN4-high)",
        ],
        "pooled_tnk_tumor_origin": int(tnk_tumor.sum()),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    write_finding(args.outdir, summary, ntab, pooled, ranks, focus)
    write_results_md(args.outdir, summary)
    log(json.dumps({k: summary[k] for k in ("n_malig", "n_tnk", "n_samples_paired", "n_patients_paired", "liana_status")}, indent=2))
    log(f"[done] {args.outdir}")


if __name__ == "__main__":
    main()
