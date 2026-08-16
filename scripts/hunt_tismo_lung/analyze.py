#!/usr/bin/env python3
"""TISMO lung-only Tacstd2/Cldn4 vs all-cancer ICB pairing (claim A4 49/64).

Honest scope
------------
TISMO in-vivo ICB gene CSVs define one comparison group per
(cell_line, study, ICB regimen). Claim A4's 49/64 is the sign of
mean(ICB) - mean(baseline) across those groups for Tacstd2.

Lung ICB in TISMO is LLC / GSE155972 only (two genotype arms).
CMT-167 and MLE12 exist but have no ICB pairing. Non-ICB lung
studies are reported descriptively from the official vivo matrix.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
NOTES = ROOT / "notes" / "hunt_tismo_lung"
RAW = NOTES / "raw"
RESULTS = ROOT / "results" / "hunt_tismo_lung"
ALIYUN = Path("/tmp/tismo_data/aliyun")

LUNG_CANCERS = {"Lung carcinoma", "Lung adenocarcinoma"}
LUNG_LINES = {"LLC", "CMT-167", "MLE12"}
TARGET_GENES = ["Tacstd2", "Cldn4"]
CONTROL_GENES = ["Actb", "Epcam", "Krt8", "Cldn3", "Cldn7", "Cd8a", "Ifng", "Gzmb", "Cd274"]
GROUP_N_RE = re.compile(r"\(n=\d+\)$")


def strip_n(label: str) -> str:
    return GROUP_N_RE.sub("", str(label)).strip()


def load_gene_csv(gene: str) -> pd.DataFrame:
    df = pd.read_csv(RAW / f"{gene}_vivo.csv")
    df = df[df["geneID"].astype(str) == gene].copy()
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df["Baseline"] = pd.to_numeric(df["Baseline"], errors="coerce")
    df["group"] = df["cell_line"].map(strip_n)
    df["cell"] = df["group"].str.split("_").str[0]
    df["is_lung"] = df["cell"].isin(LUNG_LINES)
    return df.dropna(subset=["value", "Baseline"])


def paired_group_table(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    rows = []
    for group, sub in df.groupby("group"):
        base = sub.loc[sub["Baseline"] == 1, "value"]
        icb = sub.loc[sub["Baseline"] == 0, "value"]
        if len(base) == 0 or len(icb) == 0:
            continue
        mb, mi = float(base.mean()), float(icb.mean())
        medb, medi = float(base.median()), float(icb.median())
        # TISMO DESeq2 p is repeated on treated rows
        p_deseq = pd.to_numeric(sub.loc[sub["Baseline"] == 0, "pvalue"], errors="coerce")
        p_deseq = float(p_deseq.dropna().iloc[0]) if p_deseq.notna().any() else np.nan
        gse = ",".join(sorted(sub["GSE_ID"].astype(str).unique()))
        treat = ";".join(sorted(sub["Mouse_treatment"].astype(str).unique()))
        rows.append(
            {
                "gene": gene,
                "group": group,
                "cell_line": sub["cell"].iloc[0],
                "is_lung": bool(sub["is_lung"].iloc[0]),
                "GSE_ID": gse,
                "treatments": treat,
                "n_baseline": int(len(base)),
                "n_icb": int(len(icb)),
                "mean_baseline": mb,
                "mean_icb": mi,
                "delta_mean": mi - mb,
                "median_baseline": medb,
                "median_icb": medi,
                "delta_median": medi - medb,
                "direction_mean": "up" if mi > mb else ("down" if mi < mb else "tie"),
                "direction_median": "up" if medi > medb else ("down" if medi < medb else "tie"),
                "tismo_deseq2_p": p_deseq,
            }
        )
    out = pd.DataFrame(rows).sort_values(["is_lung", "delta_mean"], ascending=[False, False])
    return out.reset_index(drop=True)


def sign_summary(tab: pd.DataFrame, label: str) -> dict:
    ups = int((tab["direction_mean"] == "up").sum())
    downs = int((tab["direction_mean"] == "down").sum())
    ties = int((tab["direction_mean"] == "tie").sum())
    n = ups + downs  # exclude ties for sign test
    deltas = tab["delta_mean"].to_numpy(dtype=float)
    wil = stats.wilcoxon(deltas, alternative="two-sided", zero_method="wilcox") if len(deltas) >= 2 else None
    ttest = stats.ttest_rel(tab["mean_icb"], tab["mean_baseline"]) if len(deltas) >= 2 else None
    binom = stats.binomtest(ups, n, 0.5, alternative="two-sided") if n >= 1 else None
    med_up = int((tab["direction_median"] == "up").sum())
    med_down = int((tab["direction_median"] == "down").sum())
    return {
        "label": label,
        "n_groups": int(len(tab)),
        "n_up_mean": ups,
        "n_down_mean": downs,
        "n_tie_mean": ties,
        "frac_up_mean": ups / len(tab) if len(tab) else np.nan,
        "mean_baseline": float(tab["mean_baseline"].mean()) if len(tab) else np.nan,
        "mean_icb": float(tab["mean_icb"].mean()) if len(tab) else np.nan,
        "mean_delta": float(tab["delta_mean"].mean()) if len(tab) else np.nan,
        "median_delta": float(tab["delta_mean"].median()) if len(tab) else np.nan,
        "wilcoxon_W": float(wil.statistic) if wil else np.nan,
        "wilcoxon_p": float(wil.pvalue) if wil else np.nan,
        "paired_t": float(ttest.statistic) if ttest else np.nan,
        "paired_t_p": float(ttest.pvalue) if ttest else np.nan,
        "binom_n": n,
        "binom_p_two_sided": float(binom.pvalue) if binom else np.nan,
        "n_up_median": med_up,
        "n_down_median": med_down,
        "note": (
            "n=2 groups: signed-rank/t-test are underpowered and are reported only as diagnostics"
            if len(tab) <= 2
            else ""
        ),
    }


def welch_mwu(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2:
        return {"welch_t": np.nan, "welch_p": np.nan, "mwu_u": np.nan, "mwu_p": np.nan}
    t = stats.ttest_ind(b, a, equal_var=False)
    u = stats.mannwhitneyu(b, a, alternative="two-sided")
    return {
        "welch_t": float(t.statistic),
        "welch_p": float(t.pvalue),
        "mwu_u": float(u.statistic),
        "mwu_p": float(u.pvalue),
    }


def sample_level_lung(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    """Per-arm sample-level tests for the only lung ICB study (GSE155972)."""
    lung = df[df["is_lung"]].copy()
    rows = []
    for group, sub in lung.groupby("group"):
        base = sub.loc[sub["Baseline"] == 1, "value"].to_numpy()
        icb = sub.loc[sub["Baseline"] == 0, "value"].to_numpy()
        stats_all = welch_mwu(base, icb)
        # leave-one-out if a treated sample is an extreme outlier (> mean+3sd of the rest)
        loo = {"loo_dropped": "", "loo_n_icb": int(len(icb)), **welch_mwu(base, icb)}
        if len(icb) >= 3:
            for i, srx in enumerate(sub.loc[sub["Baseline"] == 0, "Samples"].astype(str)):
                rest = np.delete(icb, i)
                if icb[i] > rest.mean() + 3 * (rest.std(ddof=1) if len(rest) > 1 else 0) + 1e-9:
                    loo = {"loo_dropped": srx, "loo_n_icb": int(len(rest)), **welch_mwu(base, rest)}
                    break
        rows.append(
            {
                "gene": gene,
                "group": group,
                "n_baseline": int(len(base)),
                "n_icb": int(len(icb)),
                "mean_baseline": float(np.mean(base)),
                "mean_icb": float(np.mean(icb)),
                "delta_mean": float(np.mean(icb) - np.mean(base)),
                **{f"all_{k}": v for k, v in stats_all.items()},
                **loo,
            }
        )
    return pd.DataFrame(rows)


def load_ann_vivo() -> pd.DataFrame | None:
    p = ALIYUN / "TISMO_vivosample_annotations.csv"
    if not p.exists():
        p = RAW / "TISMO_vivosample_annotations.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def load_ann_vitro() -> pd.DataFrame | None:
    p = ALIYUN / "TISMO_vitrosample_annotations.csv"
    if not p.exists():
        p = RAW / "TISMO_vitrosample_annotations.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def extract_from_rds(rds_path: Path, genes: list[str]) -> pd.DataFrame | None:
    if not rds_path.exists():
        return None
    import pyreadr

    obj = list(pyreadr.read_r(str(rds_path)).values())[0]
    present = [g for g in genes if g in obj.index]
    if not present:
        return None
    wide = obj.loc[present]
    long = wide.T.reset_index().rename(columns={"index": "sample"})
    return long


def lung_catalog(ann: pd.DataFrame) -> pd.DataFrame:
    lung = ann[ann["Cancer_type"].isin(LUNG_CANCERS)].copy()
    g = (
        lung.groupby(
            [
                "Study_ID",
                "Cell_Line",
                "Cancer_type",
                "ICB",
                "ICB_study",
                "Mouse_treatment",
                "Cell_genotype",
                "Implantation",
                "Implantation_site",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="n_samples")
    )
    g["icb_paired"] = (g["ICB"].eq(1) | g["ICB_study"].eq("Baseline")) & g["Study_ID"].eq("GSE155972")
    return g.sort_values(["icb_paired", "Study_ID", "Cell_Line"], ascending=[False, True, True])


def descriptive_lung_expr(expr: pd.DataFrame, ann: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """All TISMO lung in-vivo samples (ICB and non-ICB)."""
    lung = ann[ann["Cancer_type"].isin(LUNG_CANCERS)].copy()
    # vivo matrix columns are SRX_ID
    key = "SRX_ID" if set(expr["sample"]).intersection(set(lung["SRX_ID"].astype(str))) else "SampleName"
    merged = lung.merge(expr, left_on=key, right_on="sample", how="inner")
    rows = []
    for (study, line, icb, treat, geno, implant), sub in merged.groupby(
        ["Study_ID", "Cell_Line", "ICB", "Mouse_treatment", "Cell_genotype", "Implantation"],
        dropna=False,
    ):
        rec = {
            "Study_ID": study,
            "Cell_Line": line,
            "ICB": icb,
            "Mouse_treatment": treat,
            "Cell_genotype": geno,
            "Implantation": implant,
            "n": int(len(sub)),
        }
        for g in genes:
            if g in sub.columns:
                rec[f"{g}_mean"] = float(sub[g].mean())
                rec[f"{g}_median"] = float(sub[g].median())
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["Study_ID", "Cell_Line", "ICB"])


def vitro_cytokine(expr: pd.DataFrame, ann: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    lung = ann[ann["Cancer_type"].isin(LUNG_CANCERS) | ann["Cell_Line"].isin(LUNG_LINES)].copy()
    merged = lung.merge(expr, left_on="SampleName", right_on="sample", how="inner")
    rows = []
    for (line, clone, treat, study), sub in merged.groupby(
        ["Cell_Line", "sub_clone", "Cell_treatment", "Study_ID"], dropna=False
    ):
        rec = {
            "Study_ID": study,
            "Cell_Line": line,
            "sub_clone": clone,
            "Cell_treatment": treat,
            "n": int(len(sub)),
        }
        for g in genes:
            if g in sub.columns:
                rec[f"{g}_mean"] = float(sub[g].mean())
        rows.append(rec)
    means = pd.DataFrame(rows)

    # paired cytokine vs no_treatment within line+clone+study
    tests = []
    for (line, clone, study), sub in merged.groupby(["Cell_Line", "sub_clone", "Study_ID"], dropna=False):
        ctrl = sub[sub["Cell_treatment"].isin(["no_treatment", "nostim", "vehicle"])]
        if len(ctrl) < 2:
            continue
        for treat, tsub in sub.groupby("Cell_treatment"):
            if treat in ("no_treatment", "nostim", "vehicle") or len(tsub) < 2:
                continue
            for g in genes:
                if g not in sub.columns:
                    continue
                st = welch_mwu(ctrl[g].to_numpy(), tsub[g].to_numpy())
                tests.append(
                    {
                        "Study_ID": study,
                        "Cell_Line": line,
                        "sub_clone": clone,
                        "treatment": treat,
                        "gene": g,
                        "n_ctrl": int(len(ctrl)),
                        "n_treat": int(len(tsub)),
                        "mean_ctrl": float(ctrl[g].mean()),
                        "mean_treat": float(tsub[g].mean()),
                        "delta": float(tsub[g].mean() - ctrl[g].mean()),
                        **st,
                    }
                )
    return means, pd.DataFrame(tests)


def baseline_rank(expr: pd.DataFrame, ann: pd.DataFrame, gene: str) -> pd.DataFrame:
    """Per cell-line mean of ICB-naive (Baseline==1) in-vivo samples."""
    naive = ann[ann["Baseline"] == 1].copy()
    key = "SRX_ID" if set(expr["sample"]).intersection(set(naive["SRX_ID"].astype(str))) else "SampleName"
    merged = naive.merge(expr[["sample", gene]], left_on=key, right_on="sample", how="inner")
    g = (
        merged.groupby(["Cell_Line", "Cancer_type"], dropna=False)[gene]
        .agg(n="count", mean="mean", median="median")
        .reset_index()
        .sort_values("mean", ascending=False)
    )
    g["rank"] = np.arange(1, len(g) + 1)
    g["is_lung"] = g["Cancer_type"].isin(LUNG_CANCERS) | g["Cell_Line"].isin(LUNG_LINES)
    g["gene"] = gene
    return g


def immune_vs_gene(expr: pd.DataFrame, ann: pd.DataFrame, immune: pd.DataFrame) -> pd.DataFrame:
    lung = ann[ann["Cancer_type"].isin(LUNG_CANCERS)].copy()
    # immune columns are SampleName
    scores = ["CD8 T_mMCPcounter", "T_TIMER", "Cytotoxic_xCell", "ImmuneScore_ESTIMATE"]
    present = [s for s in scores if s in immune.index]
    if not present:
        # pick anything CD8 / ImmuneScore
        present = [i for i in immune.index if re.search(r"CD8 T_mMCP|ImmuneScore_ESTIMATE|Cytotoxic_xCell", i)]
    if not present:
        return pd.DataFrame()
    imm = immune.loc[present].T.reset_index().rename(columns={"index": "SampleName"})
    merged = lung.merge(imm, on="SampleName", how="inner")
    key = "SRX_ID"
    merged = merged.merge(expr, left_on=key, right_on="sample", how="inner")
    rows = []
    for gene in TARGET_GENES:
        if gene not in merged.columns:
            continue
        for score in present:
            if merged[score].nunique() < 3 or merged[gene].nunique() < 3:
                continue
            rho, p = stats.spearmanr(merged[gene], merged[score], nan_policy="omit")
            rows.append(
                {
                    "subset": "all_lung_vivo",
                    "n": int(merged[[gene, score]].dropna().shape[0]),
                    "gene": gene,
                    "score": score,
                    "spearman_rho": float(rho),
                    "spearman_p": float(p),
                }
            )
        icb_study = merged[merged["Study_ID"] == "GSE155972"]
        if len(icb_study) >= 6:
            for score in present:
                rho, p = stats.spearmanr(icb_study[gene], icb_study[score], nan_policy="omit")
                rows.append(
                    {
                        "subset": "GSE155972_only",
                        "n": int(icb_study[[gene, score]].dropna().shape[0]),
                        "gene": gene,
                        "score": score,
                        "spearman_rho": float(rho),
                        "spearman_p": float(p),
                    }
                )
    return pd.DataFrame(rows)


def save(df: pd.DataFrame, name: str) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / name
    df.to_csv(path, index=False)
    print("wrote", path, "rows", len(df))


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    summaries = []
    group_tables = []
    sample_tables = []
    for gene in TARGET_GENES + CONTROL_GENES:
        csv_path = RAW / f"{gene}_vivo.csv"
        if not csv_path.exists():
            print("skip missing", csv_path)
            continue
        df = load_gene_csv(gene)
        tab = paired_group_table(df, gene)
        group_tables.append(tab)
        summaries.append(sign_summary(tab, f"{gene}_all_cancer"))
        lung_tab = tab[tab["is_lung"]]
        summaries.append(sign_summary(lung_tab, f"{gene}_lung_only"))
        if gene in TARGET_GENES:
            sample_tables.append(sample_level_lung(df, gene))

    groups = pd.concat(group_tables, ignore_index=True)
    save(groups, "icb_paired_group_stats.csv")
    save(pd.DataFrame(summaries), "icb_paired_sign_summary.csv")
    if sample_tables:
        save(pd.concat(sample_tables, ignore_index=True), "lung_gse155972_sample_level.csv")

    # compact per-sample ICB values for lung
    lung_samples = []
    for gene in TARGET_GENES + CONTROL_GENES:
        p = RAW / f"{gene}_vivo.csv"
        if not p.exists():
            continue
        df = load_gene_csv(gene)
        sub = df[df["is_lung"]][
            ["Samples", "geneID", "value", "group", "Responder", "Baseline", "GSE_ID", "Mouse_treatment"]
        ]
        lung_samples.append(sub)
    if lung_samples:
        save(pd.concat(lung_samples, ignore_index=True), "lung_icb_per_sample.csv")

    ann_v = load_ann_vivo()
    ann_t = load_ann_vitro()
    genes_needed = TARGET_GENES + CONTROL_GENES

    if ann_v is not None:
        save(lung_catalog(ann_v), "lung_vivo_catalog.csv")
        # copy annotations into notes for reproducibility (small)
        dest = RAW / "TISMO_vivosample_annotations.csv"
        if not dest.exists() and (ALIYUN / "TISMO_vivosample_annotations.csv").exists():
            dest.write_bytes((ALIYUN / "TISMO_vivosample_annotations.csv").read_bytes())
        if ann_t is not None:
            destt = RAW / "TISMO_vitrosample_annotations.csv"
            if not destt.exists() and (ALIYUN / "TISMO_vitrosample_annotations.csv").exists():
                destt.write_bytes((ALIYUN / "TISMO_vitrosample_annotations.csv").read_bytes())

    vivo_expr = extract_from_rds(ALIYUN / "TISMO_expressionvivo_profiles.RDS", genes_needed)
    vitro_expr = extract_from_rds(ALIYUN / "TISMO_expressionvitro_profiles.RDS", genes_needed)
    if vivo_expr is not None:
        save(vivo_expr, "extracted_vivo_gene_matrix.csv")
        if ann_v is not None:
            save(descriptive_lung_expr(vivo_expr, ann_v, TARGET_GENES + ["Actb", "Epcam", "Cd8a"]), "lung_all_studies_descriptive.csv")
            for gene in TARGET_GENES:
                save(baseline_rank(vivo_expr, ann_v, gene), f"baseline_cellline_rank_{gene}.csv")
            try:
                import pyreadr

                imm = list(pyreadr.read_r(str(ALIYUN / "TISMO_immune_infiltration.RDS")).values())[0]
                icorr = immune_vs_gene(vivo_expr, ann_v, imm)
                if len(icorr):
                    save(icorr, "lung_immune_spearman.csv")
            except Exception as e:
                print("immune skip", e)

    if vitro_expr is not None and ann_t is not None:
        save(vitro_expr, "extracted_vitro_gene_matrix.csv")
        means, tests = vitro_cytokine(vitro_expr, ann_t, TARGET_GENES + ["Actb", "Cd274"])
        save(means, "lung_vitro_group_means.csv")
        save(tests, "lung_vitro_cytokine_tests.csv")

    # machine-readable honest verdict
    tac = [s for s in summaries if s["label"] == "Tacstd2_all_cancer"][0]
    tac_l = [s for s in summaries if s["label"] == "Tacstd2_lung_only"][0]
    cld = [s for s in summaries if s["label"] == "Cldn4_all_cancer"][0]
    cld_l = [s for s in summaries if s["label"] == "Cldn4_lung_only"][0]
    verdict = {
        "claim_A4": "Tacstd2 up in 49/64 TISMO ICB groups, p=5.8e-5",
        "all_cancer_Tacstd2": {
            "n_up_over_n": f"{tac['n_up_mean']}/{tac['n_groups']}",
            "wilcoxon_p": tac["wilcoxon_p"],
            "matches_49_of_64": tac["n_up_mean"] == 49 and tac["n_groups"] == 64,
            "wilcoxon_matches_5.8e-5": bool(tac["wilcoxon_p"] and math.isclose(tac["wilcoxon_p"], 5.8e-5, rel_tol=0.05)),
        },
        "lung_only_Tacstd2": {
            "n_up_over_n": f"{tac_l['n_up_mean']}/{tac_l['n_groups']}",
            "wilcoxon_p": tac_l["wilcoxon_p"],
            "cannot_test_49_of_64": True,
            "reason": "TISMO has only 2 lung ICB groups, both LLC GSE155972 genotype arms",
        },
        "all_cancer_Cldn4": {
            "n_up_over_n": f"{cld['n_up_mean']}/{cld['n_groups']}",
            "wilcoxon_p": cld["wilcoxon_p"],
        },
        "lung_only_Cldn4": {
            "n_up_over_n": f"{cld_l['n_up_mean']}/{cld_l['n_groups']}",
            "wilcoxon_p": cld_l["wilcoxon_p"],
        },
        "lung_icb_universe": {
            "cell_lines_with_ICB": ["LLC"],
            "studies_with_ICB": ["GSE155972"],
            "other_lung_lines_in_TISMO": ["CMT-167 (vivo, no ICB)", "MLE12 (vitro only, no ICB)"],
            "implantation_of_ICB_study": "subcutaneous flank, not orthotopic lung",
        },
    }
    (RESULTS / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
