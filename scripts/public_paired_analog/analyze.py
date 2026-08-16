#!/usr/bin/env python3
"""Public paired pre/post ICI (or chemo-IO) lung RNA analog for TACSTD2/CLDN4.

Recomputes extra public analog tables. Does not use or re-score the private
Zhejiang A7 IHC cohort. All inputs are GEO/Nature Communications source data.
"""
from __future__ import annotations

import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

GEO = Path("/tmp/geo")
OUT = Path("results/public_paired_analog/tables")
OUT.mkdir(parents=True, exist_ok=True)

TACSTD2_CLARION = "TC0100014340.hg.1"
CLDN4_CLARION = "TC0700007993.hg.1"


def wilcoxon_safe(delta: np.ndarray) -> float:
    d = np.asarray(delta, dtype=float)
    d = d[np.isfinite(d)]
    if d.size < 3 or np.allclose(d, 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(d, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def mwu(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return float("nan")
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 4:
        return float("nan"), float("nan")
    r, p = stats.spearmanr(a[mask], b[mask])
    return float(r), float(p)


def median(x) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else float("nan")


def parse_series_chars(path: Path) -> pd.DataFrame:
    titles = accs = None
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                key = vals[0].split(":", 1)[0].strip() if vals and ":" in vals[0] else f"char{len(fields)}"
                fields[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
            elif line.startswith("!series_matrix_table_begin"):
                break
    df = pd.DataFrame({"title": titles, "gsm": accs})
    for k, v in fields.items():
        df[k] = v
    return df


def load_fpkm_genes(path: Path, genes: list[str], id_col: str = "gene") -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        want = set(genes)
        rows = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in want:
                rows.append(parts)
                if len(rows) == len(want):
                    break
    df = pd.DataFrame(rows, columns=header)
    df = df.set_index(id_col)
    # drop Entrez if present
    drop = [c for c in df.columns if c.lower().startswith("entrez")]
    df = df.drop(columns=drop)
    return df.astype(float)


def load_gse207422_bulk() -> pd.DataFrame:
    meta = pd.read_excel(GEO / "gse207422_meta.xlsx")
    meta = meta.dropna(subset=["Sample", "Patient"]).copy()
    expr = pd.read_csv(GEO / "gse207422_expr.txt.gz", sep="\t", index_col=0)
    out = meta.copy()
    out["TACSTD2"] = expr.loc["TACSTD2", out["Sample"]].to_numpy()
    out["CLDN4"] = expr.loc["CLDN4", out["Sample"]].to_numpy()
    pr = out["Pathologic Response"].astype(str)
    out["response"] = np.where(pr.str.startswith("NMPR") | pr.eq("NMPR"), "NMPR", "MPR")
    return out


def load_gse166449() -> pd.DataFrame:
    meta = parse_series_chars(GEO / "GSE166449_series_matrix.txt.gz")
    tpm = pd.read_csv(GEO / "GSE166449_Raw_gene_TPM_matrix.txt.gz", sep="\t", index_col=0)
    # titles are Immunotherapy_ResponderN / nonResponderN
    meta["response"] = np.where(meta["title"].str.contains("nonResponder"), "NR", "R")
    # matrix columns may be sample titles
    cols = list(tpm.columns)
    # map by title order if columns match GSM or titles
    if set(meta["title"]).issubset(set(cols)):
        key = "title"
        mat_cols = meta["title"].tolist()
    elif set(meta["gsm"]).issubset(set(cols)):
        key = "gsm"
        mat_cols = meta["gsm"].tolist()
    else:
        # positional
        mat_cols = cols[: len(meta)]
        key = "positional"
    out = meta.copy()
    out["TACSTD2"] = np.log2(tpm.loc["TACSTD2", mat_cols].to_numpy(dtype=float) + 1)
    out["CLDN4"] = np.log2(tpm.loc["CLDN4", mat_cols].to_numpy(dtype=float) + 1)
    out["join"] = key
    return out


def load_gse248249() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = parse_series_chars(GEO / "GSE248249_series_matrix.txt.gz")
    # patient from title
    meta["patient"] = meta["title"].str.extract(r"Patient (\d+)", expand=False)
    meta["timepoint"] = meta["timepoint"].str.replace("-treatment", "", regex=False)
    expr_rows = {}
    with gzip.open(GEO / "GSE248249_series_matrix.txt.gz", "rt", errors="replace") as f:
        in_table = False
        header = None
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                header = next(f).rstrip("\n").split("\t")
                header = [h.strip('"') for h in header]
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            pid = line.split("\t", 1)[0].strip('"')
            if pid in (TACSTD2_CLARION, CLDN4_CLARION):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
                expr_rows[pid] = [float(x) for x in vals[1:]]
    gsms = header[1:]
    expr = pd.DataFrame(expr_rows, index=gsms)
    expr.columns = [{"TC0100014340.hg.1": "TACSTD2", "TC0700007993.hg.1": "CLDN4"}[c] for c in expr.columns]
    meta = meta.set_index("gsm")
    long = meta.join(expr)
    # pairs
    pairs = []
    for pid, g in long.groupby("patient"):
        pre = g[g["timepoint"].str.lower().str.contains("pre")]
        post = g[g["timepoint"].str.lower().str.contains("post")]
        if len(pre) != 1 or len(post) != 1:
            continue
        pre = pre.iloc[0]
        post = post.iloc[0]
        pairs.append(
            {
                "patient": pid,
                "pre_site": pre.get("tumor site"),
                "post_site": post.get("tumor site"),
                "same_site": pre.get("tumor site") == post.get("tumor site"),
                "lung_to_lung": pre.get("tumor site") == "Lung" and post.get("tumor site") == "Lung",
                "TACSTD2_pre": pre["TACSTD2"],
                "TACSTD2_post": post["TACSTD2"],
                "TACSTD2_delta": post["TACSTD2"] - pre["TACSTD2"],
                "CLDN4_pre": pre["CLDN4"],
                "CLDN4_post": post["CLDN4"],
                "CLDN4_delta": post["CLDN4"] - pre["CLDN4"],
            }
        )
    return long.reset_index(), pd.DataFrame(pairs)


def load_altorki() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pre = load_fpkm_genes(GEO / "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz", ["TACSTD2", "CLDN4", "EPCAM"])
    post = load_fpkm_genes(GEO / "GSE248378_Durva_Post_FPKMs.txt.gz", ["TACSTD2", "CLDN4", "EPCAM"])
    pre_meta = parse_series_chars(GEO / "GSE253564_series_matrix.txt.gz")
    post_meta = parse_series_chars(GEO / "GSE248378_series_matrix.txt.gz")

    src = pd.read_excel(GEO / "altorki/41467_source_data.xlsx", sheet_name="Figure 2e")
    src = src.dropna(subset=["Reference number"]).copy()
    src["pid"] = src["Reference number"].str.extract(r"(\d+)", expand=False).astype(int)
    src["mpr"] = np.where(src["Path Response_2Group"].astype(str).eq("Major"), "MPR", "NMPR")
    arm = pd.read_excel(GEO / "altorki/41467_source_data.xlsx", sheet_name="Figure 1b")
    arm = arm.dropna(subset=["Reference number"]).copy()
    arm["pid"] = arm["Reference number"].str.extract(r"(\d+)", expand=False).astype(int)
    mpr_map = dict(zip(src["pid"], src["mpr"]))
    arm_map = dict(zip(arm["pid"], arm["ARM"]))

    def pid_from_pre(name: str) -> int:
        return int(re.search(r"(\d+)", name).group(1))

    def pid_from_post(name: str) -> int:
        m = re.search(r"(\d+)", name)
        return int(m.group(1))

    pre_long = []
    for col in pre.columns:
        pid = pid_from_pre(col)
        pre_long.append(
            {
                "pid": pid,
                "sample": col,
                "time": "pre",
                "TACSTD2": float(pre.loc["TACSTD2", col]),
                "CLDN4": float(pre.loc["CLDN4", col]),
                "EPCAM": float(pre.loc["EPCAM", col]),
                "mpr": mpr_map.get(pid, "unknown"),
                "arm": arm_map.get(pid, "unknown"),
            }
        )
    post_long = []
    for col in post.columns:
        pid = pid_from_post(col)
        post_long.append(
            {
                "pid": pid,
                "sample": col,
                "time": "post",
                "TACSTD2": float(post.loc["TACSTD2", col]),
                "CLDN4": float(post.loc["CLDN4", col]),
                "EPCAM": float(post.loc["EPCAM", col]),
                "mpr": mpr_map.get(pid, "unknown"),
                "arm": arm_map.get(pid, "unknown"),
            }
        )
    pre_df = pd.DataFrame(pre_long)
    post_df = pd.DataFrame(post_long)
    pairs = []
    for pid in sorted(set(pre_df["pid"]) & set(post_df["pid"])):
        a = pre_df[pre_df["pid"] == pid].iloc[0]
        b = post_df[post_df["pid"] == pid].iloc[0]
        pairs.append(
            {
                "pid": pid,
                "pre_sample": a["sample"],
                "post_sample": b["sample"],
                "mpr": a["mpr"],
                "arm": a["arm"],
                "TACSTD2_pre": a["TACSTD2"],
                "TACSTD2_post": b["TACSTD2"],
                "TACSTD2_delta": b["TACSTD2"] - a["TACSTD2"],
                "TACSTD2_log2fc": math.log2((b["TACSTD2"] + 0.1) / (a["TACSTD2"] + 0.1)),
                "CLDN4_pre": a["CLDN4"],
                "CLDN4_post": b["CLDN4"],
                "CLDN4_delta": b["CLDN4"] - a["CLDN4"],
                "CLDN4_log2fc": math.log2((b["CLDN4"] + 0.1) / (a["CLDN4"] + 0.1)),
                "EPCAM_pre": a["EPCAM"],
                "EPCAM_post": b["EPCAM"],
                "EPCAM_log2fc": math.log2((b["EPCAM"] + 0.1) / (a["EPCAM"] + 0.1)),
            }
        )
    return pre_df, post_df, pd.DataFrame(pairs)


def load_gse207422_scrna_sample_means() -> pd.DataFrame:
    """Sample-level mean log1p UMI for TACSTD2/CLDN4; epithelial via EPCAM>0 cells."""
    scmeta = pd.read_excel(GEO / "gse207422_scmeta.xlsx")
    scmeta = scmeta.dropna(subset=["Sample", "Patient"]).copy()
    scmeta["time"] = np.where(scmeta["Resource"].astype(str).str.contains("Pre", case=False), "pre", "post")
    pr = scmeta["Pathologic Response"].astype(str)
    scmeta["response"] = np.where(
        pr.eq("NE"),
        "NE",
        np.where(pr.str.startswith("NMPR") | pr.eq("NMPR"), "NMPR", "MPR"),
    )
    genes = ["TACSTD2", "CLDN4", "EPCAM"]
    # stream UMI matrix: first column gene, remaining cells
    with gzip.open(GEO / "GSE207422_scRNAseq_UMI_matrix.txt.gz", "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cells = header[1:]
        rows = {}
        for line in f:
            g = line.split("\t", 1)[0]
            if g in genes:
                vals = np.array([float(x) for x in line.rstrip("\n").split("\t")[1:]], dtype=np.float32)
                rows[g] = vals
                if len(rows) == len(genes):
                    break
    # sample prefix: BD_immuneNN_...
    sample = pd.Series(cells).str.extract(r"(BD_immune\d+)", expand=False)
    log1p = {g: np.log1p(rows[g]) for g in rows}
    epi = rows["EPCAM"] > 0
    recs = []
    for sid, idx in sample.groupby(sample).groups.items():
        idx = np.array(list(idx))
        rec = {"Sample": sid, "n_cells": int(len(idx)), "n_epcam_pos": int(epi[idx].sum())}
        for g in genes:
            rec[f"{g}_all"] = float(np.mean(log1p[g][idx]))
            rec[f"{g}_epi"] = float(np.mean(log1p[g][idx][epi[idx]])) if epi[idx].any() else float("nan")
        recs.append(rec)
    means = pd.DataFrame(recs)
    return means.merge(scmeta, on="Sample", how="left")


def summarize_pairs(pairs: pd.DataFrame, gene: str, subset_name: str, mask=None) -> dict:
    d = pairs if mask is None else pairs.loc[mask]
    delta = d[f"{gene}_delta"].to_numpy(dtype=float)
    n = int(np.isfinite(delta).sum())
    n_up = int((delta > 0).sum())
    n_down = int((delta < 0).sum())
    return {
        "cohort": "GSE253564_pre__GSE248378_post" if "TACSTD2_log2fc" in d.columns else "GSE248249",
        "subset": subset_name,
        "gene": gene,
        "n_pairs": n,
        "n_up": n_up,
        "n_down": n_down,
        "median_pre": median(d[f"{gene}_pre"]),
        "median_post": median(d[f"{gene}_post"]),
        "median_delta": median(delta),
        "wilcoxon_p": wilcoxon_safe(delta),
    }


def main() -> None:
    inventory = []

    # --- Altorki leftover timepoints: GSE253564 pre + GSE248378 post ---
    pre_df, post_df, alt_pairs = load_altorki()
    pre_df.to_csv(OUT / "altorki_GSE253564_pre_sample.csv", index=False)
    post_df.to_csv(OUT / "altorki_GSE248378_post_sample.csv", index=False)
    alt_pairs.to_csv(OUT / "altorki_paired_pre_post.csv", index=False)

    analog_rows = []
    for gene in ["TACSTD2", "CLDN4"]:
        analog_rows.append(summarize_pairs(alt_pairs, gene, "all_joinable_pairs"))
        analog_rows.append(summarize_pairs(alt_pairs, gene, "NMPR_pairs", alt_pairs["mpr"] == "NMPR"))
        analog_rows.append(summarize_pairs(alt_pairs, gene, "MPR_pairs", alt_pairs["mpr"] == "MPR"))
        analog_rows.append(summarize_pairs(alt_pairs, gene, "Arm1_durva_alone", alt_pairs["arm"].astype(str).str.contains("Durvalumab") & ~alt_pairs["arm"].astype(str).str.contains("SBRT|\\+")))
        analog_rows.append(summarize_pairs(alt_pairs, gene, "Arm2_durva_SBRT", alt_pairs["arm"].astype(str).str.contains("SBRT|\\+")))

    # pretreatment MPR vs NMPR (GSE253564)
    for gene in ["TACSTD2", "CLDN4"]:
        mpr = pre_df.loc[pre_df["mpr"] == "MPR", gene]
        nmpr = pre_df.loc[pre_df["mpr"] == "NMPR", gene]
        analog_rows.append(
            {
                "cohort": "GSE253564_pretreatment",
                "subset": "MPR_vs_NMPR",
                "gene": gene,
                "n_pairs": "",
                "n_mpr": int(mpr.size),
                "n_nmpr": int(nmpr.size),
                "median_mpr": median(mpr),
                "median_nmpr": median(nmpr),
                "mwu_p": mwu(mpr, nmpr),
            }
        )
        mprp = post_df.loc[post_df["mpr"] == "MPR", gene]
        nmprp = post_df.loc[post_df["mpr"] == "NMPR", gene]
        analog_rows.append(
            {
                "cohort": "GSE248378_post_resection",
                "subset": "MPR_vs_NMPR",
                "gene": gene,
                "n_pairs": "",
                "n_mpr": int(mprp.size),
                "n_nmpr": int(nmprp.size),
                "median_mpr": median(mprp),
                "median_nmpr": median(nmprp),
                "mwu_p": mwu(mprp, nmprp),
            }
        )

    inventory.append(
        {
            "resource": "GSE253564 + GSE248378 (Altorki NCT02904954)",
            "cancer": "early NSCLC",
            "treatment": "neoadjuvant durvalumab ± SBRT",
            "paired_pre_post": "yes_after_joining_two_GEO_series",
            "n_pre": int(len(pre_df)),
            "n_post": int(len(post_df)),
            "n_pairs": int(len(alt_pairs)),
            "n_pairs_MPR": int((alt_pairs["mpr"] == "MPR").sum()),
            "n_pairs_NMPR": int((alt_pairs["mpr"] == "NMPR").sum()),
            "note": "Leftover timepoints: pre FPKMs in GSE253564, post FPKMs in GSE248378. Numeric IDs join 20 patients. All 20 joinable pairs are NMPR; MPR posts were not deposited in GSE248378.",
        }
    )

    # --- GSE207422 bulk pretreatment MPR ---
    bulk = load_gse207422_bulk()
    bulk.to_csv(OUT / "gse207422_bulk_pretreatment.csv", index=False)
    for gene in ["TACSTD2", "CLDN4"]:
        mpr = bulk.loc[bulk["response"] == "MPR", gene]
        nmpr = bulk.loc[bulk["response"] == "NMPR", gene]
        analog_rows.append(
            {
                "cohort": "GSE207422_bulk_pretreatment",
                "subset": "MPR_vs_NMPR",
                "gene": gene,
                "n_mpr": int(mpr.size),
                "n_nmpr": int(nmpr.size),
                "median_mpr": median(mpr),
                "median_nmpr": median(nmpr),
                "mwu_p": mwu(mpr, nmpr),
            }
        )
    r, p = spearman(bulk["TACSTD2"], bulk["CLDN4"])
    analog_rows.append(
        {
            "cohort": "GSE207422_bulk_pretreatment",
            "subset": "TACSTD2_vs_CLDN4",
            "gene": "both",
            "n": int(len(bulk)),
            "spearman_rho": r,
            "spearman_p": p,
        }
    )
    inventory.append(
        {
            "resource": "GSE207422 bulk (Hu Genome Med 2023)",
            "cancer": "resectable NSCLC",
            "treatment": "neoadjuvant PD-1 + chemo",
            "paired_pre_post": "no_deposited_pre_only",
            "n_pre": int(len(bulk)),
            "n_post": 0,
            "n_pairs": 0,
            "note": "24 pretreatment biopsies with MPR/NMPR. Not paired.",
        }
    )

    # --- GSE207422 scRNA unpaired ---
    sc = load_gse207422_scrna_sample_means()
    sc.to_csv(OUT / "gse207422_scrna_sample_means.csv", index=False)
    pre = sc[sc["time"] == "pre"]
    post = sc[sc["time"] == "post"]
    post_mpr = post[post["response"] == "MPR"]
    post_nmpr = post[post["response"] == "NMPR"]
    for gene, col in [("TACSTD2", "TACSTD2_epi"), ("CLDN4", "CLDN4_epi")]:
        analog_rows.append(
            {
                "cohort": "GSE207422_scRNA_epithelial_unpaired",
                "subset": "post_vs_pre",
                "gene": gene,
                "n_pre": int(pre[col].notna().sum()),
                "n_post": int(post[col].notna().sum()),
                "median_pre": median(pre[col]),
                "median_post": median(post[col]),
                "mwu_p": mwu(pre[col], post[col]),
            }
        )
        analog_rows.append(
            {
                "cohort": "GSE207422_scRNA_epithelial_unpaired",
                "subset": "post_MPR_vs_NMPR",
                "gene": gene,
                "n_mpr": int(post_mpr[col].notna().sum()),
                "n_nmpr": int(post_nmpr[col].notna().sum()),
                "median_mpr": median(post_mpr[col]),
                "median_nmpr": median(post_nmpr[col]),
                "mwu_p": mwu(post_mpr[col], post_nmpr[col]),
            }
        )
    inventory.append(
        {
            "resource": "GSE207422 scRNA (Hu Genome Med 2023)",
            "cancer": "resectable NSCLC",
            "treatment": "neoadjuvant PD-1 + chemo",
            "paired_pre_post": "no_3pre_12post_different_patients",
            "n_pre": int((sc["time"] == "pre").sum()),
            "n_post": int((sc["time"] == "post").sum()),
            "n_pairs": 0,
            "note": "0 same-patient pairs. Epithelial = EPCAM>0 cells; score = mean log1p UMI.",
        }
    )

    # --- GSE166449 unpaired pembro ---
    luad = load_gse166449()
    luad.to_csv(OUT / "gse166449_pretreatment.csv", index=False)
    for gene in ["TACSTD2", "CLDN4"]:
        r = luad.loc[luad["response"] == "R", gene]
        nr = luad.loc[luad["response"] == "NR", gene]
        analog_rows.append(
            {
                "cohort": "GSE166449_pembro_pretreatment",
                "subset": "RECIST_R_vs_NR",
                "gene": gene,
                "n_R": int(r.size),
                "n_NR": int(nr.size),
                "median_R": median(r),
                "median_NR": median(nr),
                "mwu_p": mwu(r, nr),
            }
        )
    inventory.append(
        {
            "resource": "GSE166449 (Lee Cell 2021 Samsung pembro)",
            "cancer": "advanced LUAD",
            "treatment": "pembrolizumab",
            "paired_pre_post": "no_pretreatment_only",
            "n_pre": int(len(luad)),
            "n_post": 0,
            "n_pairs": 0,
            "note": "22 pretreatment tumors; 7 RECIST responders vs 15 non-responders. Not paired. Not MPR.",
        }
    )

    # --- GSE248249 acquired-resistance pairs ---
    long249, pairs249 = load_gse248249()
    long249.to_csv(OUT / "gse248249_sample_long.csv", index=False)
    pairs249.to_csv(OUT / "gse248249_pairs.csv", index=False)
    for gene in ["TACSTD2", "CLDN4"]:
        analog_rows.append(summarize_pairs(pairs249, gene, "all_pairs"))
        analog_rows.append(summarize_pairs(pairs249, gene, "same_site", pairs249["same_site"]))
        analog_rows.append(summarize_pairs(pairs249, gene, "lung_to_lung", pairs249["lung_to_lung"]))
    inventory.append(
        {
            "resource": "GSE248249 (Memon Cancer Cell 2024)",
            "cancer": "advanced NSCLC",
            "treatment": "PD-(L)1; post = acquired resistance",
            "paired_pre_post": "yes_13_pairs",
            "n_pre": int((long249["timepoint"].str.contains("Pre", case=False)).sum()),
            "n_post": int((long249["timepoint"].str.contains("Post", case=False)).sum()),
            "n_pairs": int(len(pairs249)),
            "n_lung_to_lung": int(pairs249["lung_to_lung"].sum()),
            "note": "Affymetrix Clariom D RMA. 12/13 pairs are not lung-to-lung. No MPR labels (metastatic acquired resistance).",
        }
    )

    # closest non-lung / locked
    inventory.extend(
        [
            {
                "resource": "GSE91061 (Riaz Cell 2017)",
                "cancer": "melanoma",
                "treatment": "nivolumab",
                "paired_pre_post": "yes_pre_on_not_lung",
                "note": "Closest public non-lung paired ICI bulk RNA. Not scored here (wrong tissue).",
            },
            {
                "resource": "GSE291670 (Yang JTM 2025)",
                "cancer": "NSCLC",
                "treatment": "neoadjuvant anlotinib + camrelizumab",
                "paired_pre_post": "no_6_post_only_scRNA",
                "note": "MPR n=3 vs NMPR n=3 post-resection scRNA. Pretreatment n=3 from HRA001033, not this series.",
            },
            {
                "resource": "GSE179994 (Liu Nat Cancer 2022)",
                "cancer": "NSCLC",
                "treatment": "pembro + chemo",
                "paired_pre_post": "partial_Tcell_only",
                "note": "Paired biopsies exist but deposited matrix is T-cell/TCR focused; not a tumor TACSTD2 analog.",
            },
            {
                "resource": "GSE202417",
                "cancer": "NSCLC",
                "treatment": "nivo + bezafibrate",
                "paired_pre_post": "yes_but_PBMC_CD8",
                "note": "Paired blood CD8 arrays, not tumor.",
            },
            {
                "resource": "GSE241934 / HRA007419 (NEOTIDE)",
                "cancer": "EGFR-mutant NSCLC",
                "treatment": "neoadjuvant sintilimab + chemo",
                "paired_pre_post": "bulk_paired_but_GSA_controlled",
                "note": "Paper reports paired bulk RNA; public GEO is post-resection scRNA only. Bulk is HRA007419 (application).",
            },
            {
                "resource": "EGAD00001011302 (CA209-153)",
                "cancer": "advanced NSCLC",
                "treatment": "nivolumab",
                "paired_pre_post": "yes_EGA_DAC",
                "note": "24 pre and 12 on-therapy RNA-seq. Controlled access. Not used.",
            },
            {
                "resource": "HRA006493 (GSA)",
                "cancer": "NSCLC",
                "treatment": "neoadjuvant PD-1 + chemo or anti-VEGFA",
                "paired_pre_post": "yes_GSA_DAC",
                "note": "Paired pre/post scRNA claimed. Controlled access. Not used.",
            },
            {
                "resource": "EGAS00001007753 NEOPREDICT-Lung",
                "cancer": "NSCLC",
                "treatment": "neoadjuvant nivo ± relatlimab",
                "paired_pre_post": "NanoString_17_pairs_EGA",
                "note": "17 paired biopsies, PanCancer Immune/Pathway panels. EGA. Not used.",
            },
            {
                "resource": "Cellular Oncology 2025 CAF/metabolism 13-pair bulk",
                "cancer": "NSCLC",
                "treatment": "neoadjuvant PD-1 + chemo",
                "paired_pre_post": "reported_13_pairs_not_in_GEO",
                "note": "Authors describe 13 paired bulk RNA-seq; no GEO accession found. Not public.",
            },
            {
                "resource": "GSE243013 (Cell 2025)",
                "cancer": "NSCLC",
                "treatment": "neoadjuvant chemo-IO",
                "paired_pre_post": "mostly_post_scRNA",
                "note": "243 post-treatment scRNA/TCR. Matched pre FFPE for a subset is not a downloadable bulk TACSTD2 matrix.",
            },
        ]
    )

    analog = pd.DataFrame(analog_rows)
    analog.to_csv(OUT / "public_analog_TACSTD2_CLDN4.csv", index=False)
    inv = pd.DataFrame(inventory)
    inv.to_csv(OUT / "public_inventory_2023_2026.csv", index=False)

    summary = {
        "A7_private_taken_as_given": {
            "n_pairs": 25,
            "TROP2_Hscore_pre_to_post": "94→121",
            "NMPR_vs_MPR_Hscore": "147 vs 88",
        },
        "altorki_joinable_pairs": int(len(alt_pairs)),
        "altorki_pairs_MPR": int((alt_pairs["mpr"] == "MPR").sum()),
        "altorki_pairs_NMPR": int((alt_pairs["mpr"] == "NMPR").sum()),
        "gse248249_pairs": int(len(pairs249)),
        "gse248249_lung_to_lung": int(pairs249["lung_to_lung"].sum()),
        "gse166449_n": int(len(luad)),
        "gse207422_bulk_n": int(len(bulk)),
        "gse207422_scrna_n": int(len(sc)),
    }
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(analog.to_string(index=False))


if __name__ == "__main__":
    main()
