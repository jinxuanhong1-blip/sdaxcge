#!/usr/bin/env python3
"""C4 analog: GSE245459 SKOV3 shTACSTD2 — CLDN4 and IFN/APM.

Primary contrast is shTACSTD2 vs shNC (no cisplatin). Secondary is the same
knockdown on a cisplatin background. FPKM is not counts; n=3/group. Report
exactly what the numbers say.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "w200" / "C4_GSE245459"
RAW = OUT / "raw" / "GSE245459_fpkm.anno.txt.gz"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import APM, C4_IFN_MHCI, IFN, JUNCTION, RAP1_PI3K_AKT, TARGETS, all_sets

SAMPLE_COLS = [
    "shNC1",
    "shNC2",
    "shNC3",
    "shNCDDP1",
    "shNCDDP2",
    "shNCDDP3",
    "sh1",
    "sh2",
    "sh3",
    "shDDP1",
    "shDDP2",
    "shDDP3",
]

GROUPS = {
    "shNC": ["shNC1", "shNC2", "shNC3"],
    "shTACSTD2": ["sh1", "sh2", "sh3"],
    "shNC_DDP": ["shNCDDP1", "shNCDDP2", "shNCDDP3"],
    "shTACSTD2_DDP": ["shDDP1", "shDDP2", "shDDP3"],
}

# Pre-specified contrasts. First is the C4 analog.
CONTRASTS = [
    ("KD_noDDP", "shTACSTD2", "shNC"),
    ("KD_DDP", "shTACSTD2_DDP", "shNC_DDP"),
    ("DDP_in_NC", "shNC_DDP", "shNC"),
    ("DDP_in_KD", "shTACSTD2_DDP", "shTACSTD2"),
]

PSEUDO = 1.0
MIN_MEAN_FPKM = 1.0


def load_fpkm() -> pd.DataFrame:
    df = pd.read_csv(RAW, sep="\t", compression="gzip", low_memory=False)
    for c in SAMPLE_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["GeneName"] = df["GeneName"].astype(str).str.strip()
    df = df[df["GeneName"].notna() & (df["GeneName"] != "") & (df["GeneName"] != "nan")]
    if "Biotype" in df.columns:
        coding = df[df["Biotype"].astype(str) == "protein_coding"].copy()
        rest = df[df["Biotype"].astype(str) != "protein_coding"].copy()
        df = pd.concat([coding, rest], ignore_index=True)
    df["mean_fpkm"] = df[SAMPLE_COLS].mean(axis=1)
    df = df.sort_values("mean_fpkm", ascending=False)
    df = df.drop_duplicates("GeneName", keep="first")
    return df.set_index("GeneName")


def welch_contrast(logx: pd.DataFrame, treat: list[str], ctrl: list[str]) -> pd.DataFrame:
    a = logx[treat].to_numpy(dtype=float)
    b = logx[ctrl].to_numpy(dtype=float)
    n_a = np.isfinite(a).sum(axis=1)
    n_b = np.isfinite(b).sum(axis=1)
    mean_a = np.nanmean(a, axis=1)
    mean_b = np.nanmean(b, axis=1)
    var_a = np.nanvar(a, axis=1, ddof=1)
    var_b = np.nanvar(b, axis=1, ddof=1)
    log2fc = mean_a - mean_b
    se2 = var_a / n_a + var_b / n_b
    tvals = log2fc / np.sqrt(se2)
    df = se2**2 / ((var_a / n_a) ** 2 / (n_a - 1) + (var_b / n_b) ** 2 / (n_b - 1))
    pvals = 2 * stats.t.sf(np.abs(tvals), df)
    bad = (n_a < 2) | (n_b < 2) | ~np.isfinite(se2) | (se2 == 0)
    tvals[bad] = np.nan
    pvals[bad] = np.nan
    out = pd.DataFrame(
        {
            "gene": logx.index,
            "mean_treat_log2": mean_a,
            "mean_ctrl_log2": mean_b,
            "log2FC": log2fc,
            "t": tvals,
            "p": pvals,
        }
    )
    mask = out["p"].notna()
    q = np.full(len(out), np.nan)
    if mask.sum() > 0:
        q[mask.to_numpy()] = multipletests(out.loc[mask, "p"], method="fdr_bh")[1]
    out["q"] = q
    return out


def geneset_mw(de: pd.DataFrame, genes: list[str], background_ok: pd.Series) -> dict:
    present = [g for g in genes if g in set(de["gene"])]
    set_fc = de.loc[de["gene"].isin(present), "log2FC"].dropna()
    bg_fc = de.loc[background_ok & ~de["gene"].isin(present), "log2FC"].dropna()
    if len(set_fc) < 3 or len(bg_fc) < 20:
        return {
            "n_in_set": len(genes),
            "n_present": len(present),
            "n_tested": int(len(set_fc)),
            "median_log2FC_set": float(set_fc.median()) if len(set_fc) else None,
            "median_log2FC_bg": float(bg_fc.median()) if len(bg_fc) else None,
            "mw_p": None,
            "direction": "untestable",
        }
    u, p = stats.mannwhitneyu(set_fc, bg_fc, alternative="two-sided")
    med_s = float(set_fc.median())
    med_b = float(bg_fc.median())
    if p >= 0.05:
        direction = "null"
    elif med_s > med_b:
        direction = "up"
    else:
        direction = "down"
    return {
        "n_in_set": len(genes),
        "n_present": len(present),
        "n_tested": int(len(set_fc)),
        "median_log2FC_set": med_s,
        "median_log2FC_bg": med_b,
        "delta_median": med_s - med_b,
        "mw_U": float(u),
        "mw_p": float(p),
        "direction": direction,
        "n_up": int((set_fc > 0).sum()),
        "n_down": int((set_fc < 0).sum()),
    }


def zmean_score(logx: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in logx.index]
    sub = logx.loc[present]
    z = (sub.sub(sub.mean(axis=1), axis=0)).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def score_contrast(scores: pd.Series, treat: list[str], ctrl: list[str]) -> dict:
    a = scores[treat].to_numpy(dtype=float)
    b = scores[ctrl].to_numpy(dtype=float)
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return {
        "mean_treat": float(np.mean(a)),
        "mean_ctrl": float(np.mean(b)),
        "delta": float(np.mean(a) - np.mean(b)),
        "t": float(t),
        "p": float(p),
    }


def sign_test(de: pd.DataFrame, genes: list[str]) -> dict:
    sub = de[de["gene"].isin(genes)].dropna(subset=["log2FC"])
    n = len(sub)
    n_up = int((sub["log2FC"] > 0).sum())
    n_down = int((sub["log2FC"] < 0).sum())
    if n == 0:
        return {"n": 0, "n_up": 0, "n_down": 0, "binom_p_two_sided": None}
    p = float(stats.binomtest(n_up, n=n, p=0.5, alternative="two-sided").pvalue)
    return {"n": n, "n_up": n_up, "n_down": n_down, "binom_p_two_sided": p}


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    fpkm = load_fpkm()
    logx = np.log2(fpkm[SAMPLE_COLS] + PSEUDO)
    mean_fpkm = fpkm[SAMPLE_COLS].mean(axis=1)
    expressed = mean_fpkm >= MIN_MEAN_FPKM

    slim = fpkm[SAMPLE_COLS].copy()
    slim.insert(0, "ensembl", fpkm["name"] if "name" in fpkm.columns else "")
    slim.to_csv(TABLES / "fpkm_gene_matrix.tsv.gz", sep="\t", compression="gzip")

    de_tables = {}
    for name, treat_g, ctrl_g in CONTRASTS:
        de = welch_contrast(logx, GROUPS[treat_g], GROUPS[ctrl_g])
        de["mean_fpkm"] = mean_fpkm.reindex(de["gene"]).values
        de["expressed"] = expressed.reindex(de["gene"]).values
        de_tables[name] = de
        de.to_csv(TABLES / f"de_{name}.tsv.gz", sep="\t", index=False, compression="gzip")

    # Marker table for pre-specified genes.
    marker_genes = []
    seen = set()
    for g in TARGETS + C4_IFN_MHCI + APM + IFN + JUNCTION + RAP1_PI3K_AKT:
        if g not in seen:
            marker_genes.append(g)
            seen.add(g)

    rows = []
    for g in marker_genes:
        rec = {"gene": g, "present": g in fpkm.index}
        if g in fpkm.index:
            rec["mean_fpkm"] = float(mean_fpkm.loc[g])
            rec["expressed"] = bool(expressed.loc[g])
            for grp, cols in GROUPS.items():
                rec[f"fpkm_{grp}"] = float(fpkm.loc[g, cols].mean())
                rec[f"log2_{grp}"] = float(logx.loc[g, cols].mean())
            for cname, de in de_tables.items():
                hit = de[de["gene"] == g]
                if len(hit):
                    rec[f"{cname}_log2FC"] = float(hit["log2FC"].iloc[0])
                    rec[f"{cname}_p"] = float(hit["p"].iloc[0])
                    rec[f"{cname}_q"] = float(hit["q"].iloc[0])
        rows.append(rec)
    markers = pd.DataFrame(rows)
    markers.to_csv(TABLES / "markers.tsv", sep="\t", index=False)

    # Gene-set stats per contrast.
    gs_rows = []
    gs_summary = {}
    for cname, de in de_tables.items():
        bg = de["expressed"].fillna(False)
        gs_summary[cname] = {}
        for sname, genes in all_sets().items():
            st = geneset_mw(de, genes, bg)
            st.update({"contrast": cname, "set": sname})
            st.update(sign_test(de, genes))
            gs_rows.append(st)
            gs_summary[cname][sname] = st
    gs_df = pd.DataFrame(gs_rows)
    gs_df.to_csv(TABLES / "geneset_stats.tsv", sep="\t", index=False)

    # Sample-level signature scores (z-mean).
    score_mat = pd.DataFrame(
        {sname: zmean_score(logx, genes) for sname, genes in all_sets().items()}
    )
    score_mat.index.name = "sample"
    meta = []
    for s in score_mat.index:
        if s.startswith("shNCDDP"):
            kd, ddp = "shNC", "DDP"
        elif s.startswith("shNC"):
            kd, ddp = "shNC", "noDDP"
        elif s.startswith("shDDP"):
            kd, ddp = "shTACSTD2", "DDP"
        else:
            kd, ddp = "shTACSTD2", "noDDP"
        meta.append({"sample": s, "kd": kd, "ddp": ddp})
    meta_df = pd.DataFrame(meta).set_index("sample")
    scores_out = score_mat.join(meta_df)
    scores_out.to_csv(TABLES / "signature_scores.tsv", sep="\t")

    score_contrasts = {}
    for cname, treat_g, ctrl_g in CONTRASTS:
        score_contrasts[cname] = {
            sname: score_contrast(score_mat[sname], GROUPS[treat_g], GROUPS[ctrl_g])
            for sname in score_mat.columns
        }

    # Interaction on log2: (KD_DDP - NC_DDP) - (KD - NC), Welch on per-rep diffs
    # is not paired; report mean interaction of group means + bootstrap of
    # unpaired difference-of-differences via sample-level scores only.
    interaction = {}
    for sname in score_mat.columns:
        d_no = score_mat.loc[GROUPS["shTACSTD2"], sname].mean() - score_mat.loc[GROUPS["shNC"], sname].mean()
        d_ddp = (
            score_mat.loc[GROUPS["shTACSTD2_DDP"], sname].mean()
            - score_mat.loc[GROUPS["shNC_DDP"], sname].mean()
        )
        interaction[sname] = {
            "delta_KD_noDDP": float(d_no),
            "delta_KD_DDP": float(d_ddp),
            "interaction_DDP_minus_noDDP": float(d_ddp - d_no),
        }

    # Key-gene compact JSON.
    def gene_pack(g: str) -> dict:
        if g not in fpkm.index:
            return {"present": False}
        pack = {
            "present": True,
            "mean_fpkm": float(mean_fpkm.loc[g]),
            "fpkm_by_group": {k: float(fpkm.loc[g, v].mean()) for k, v in GROUPS.items()},
            "fpkm_replicates": {c: float(fpkm.loc[g, c]) for c in SAMPLE_COLS},
        }
        for cname, de in de_tables.items():
            hit = de[de["gene"] == g].iloc[0]
            pack[cname] = {
                "log2FC": float(hit["log2FC"]),
                "p": float(hit["p"]),
                "q": float(hit["q"]),
            }
        return pack

    key_stats = {
        "dataset": "GSE245459",
        "system": "SKOV3 human ovarian carcinoma, bulk RNA-seq FPKM",
        "perturbation": "shTACSTD2 vs shNC; ± cisplatin 10 µg/ml 48 h (author protocol)",
        "n_per_group": 3,
        "n_genes_after_collapse": int(fpkm.shape[0]),
        "n_protein_coding_approx": int((fpkm.get("Biotype", "") == "protein_coding").sum())
        if "Biotype" in fpkm.columns
        else None,
        "primary_contrast": "KD_noDDP = shTACSTD2 (sh1-3) vs shNC (shNC1-3)",
        "secondary_contrast": "KD_DDP = shDDP vs shNCDDP",
        "method": "log2(FPKM+1); Welch t-test n=3 vs n=3; BH-FDR genome-wide; gene-set Mann-Whitney of log2FC vs expressed background; sample scores = mean of per-gene z across 12 samples",
        "claim_C4": "CLDN4 KD opens IFN/MHC-I (IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A). This slice is a TACSTD2-KD analog, not a CLDN4-KD test.",
        "TACSTD2": gene_pack("TACSTD2"),
        "CLDN4": gene_pack("CLDN4"),
        "C4_panel_genes": {g: gene_pack(g) for g in C4_IFN_MHCI},
        "APM_genes": {g: gene_pack(g) for g in APM},
        "geneset": gs_summary,
        "signature_score_contrasts": score_contrasts,
        "signature_interaction": interaction,
        "knockdown_qc": None,
    }

    # Knockdown QC verdict.
    t2 = key_stats["TACSTD2"]
    if t2.get("present"):
        fc = t2["KD_noDDP"]["log2FC"]
        p = t2["KD_noDDP"]["p"]
        if fc <= -0.5 and p < 0.05:
            qc = "TACSTD2 RNA is down in sh vs shNC — knockdown is visible on FPKM."
        elif fc < 0:
            qc = (
                "TACSTD2 RNA is only modestly lower in sh vs shNC. "
                "Treat knockdown as partial / protein-unverified at the RNA layer."
            )
        else:
            qc = (
                "TACSTD2 RNA is NOT down in sh vs shNC. "
                "Do not treat this series as a confirmed transcriptional knockdown."
            )
        key_stats["knockdown_qc"] = {
            "KD_noDDP_log2FC": fc,
            "KD_noDDP_p": p,
            "KD_DDP_log2FC": t2["KD_DDP"]["log2FC"],
            "KD_DDP_p": t2["KD_DDP"]["p"],
            "verdict": qc,
        }

    # Honest C4 analog verdicts (rule-based, not narrative spin).
    def analog_verdict() -> dict:
        cldn = key_stats["CLDN4"]
        c4set = gs_summary["KD_noDDP"]["C4_IFN_MHCI"]
        apmset = gs_summary["KD_noDDP"]["APM"]
        ifnset = gs_summary["KD_noDDP"]["IFN"]
        cldn_fc = cldn.get("KD_noDDP", {}).get("log2FC")
        cldn_p = cldn.get("KD_noDDP", {}).get("p")
        if cldn_fc is None:
            cldn_call = "untestable"
        elif cldn_p < 0.05 and cldn_fc <= -0.25:
            cldn_call = "CLDN4 down after TACSTD2 KD (supports TACSTD2→CLDN4 in this line)"
        elif cldn_p < 0.05 and cldn_fc >= 0.25:
            cldn_call = "CLDN4 up after TACSTD2 KD (opposite of TACSTD2 maintaining CLDN4)"
        else:
            cldn_call = "CLDN4 not significantly moved (no support for TACSTD2→CLDN4 at RNA)"

        def open_call(st, label):
            if st.get("mw_p") is None:
                return f"{label}: untestable"
            if st["mw_p"] < 0.05 and st["direction"] == "up":
                return f"{label}: UP vs background (direction matches C4 'opening')"
            if st["mw_p"] < 0.05 and st["direction"] == "down":
                return f"{label}: DOWN vs background (opposite of C4 'opening')"
            return f"{label}: null vs background (does not reproduce C4 opening)"

        return {
            "CLDN4": cldn_call,
            "C4_IFN_MHCI": open_call(c4set, "C4 IFN/MHC-I panel"),
            "APM": open_call(apmset, "APM"),
            "IFN": open_call(ifnset, "IFN/ISG"),
            "bottom_line": (
                "This is a TACSTD2-KD analog in ovarian SKOV3, not a CLDN4-KD test "
                "and not lung. Support for C4 requires CLDN4 down and IFN/APM up. "
                "Read the four calls above; do not upgrade a null or opposite result."
            ),
        }

    key_stats["analog_verdict"] = analog_verdict()

    with open(TABLES / "key_stats.json", "w") as f:
        json.dump(key_stats, f, indent=2, allow_nan=True)

    # Figures.
    plt.rcParams.update({"font.size": 9, "figure.dpi": 140})

    def grouped_boxes(ax, gene_or_score, ylabel, title):
        order = ["shNC", "shTACSTD2", "shNC_DDP", "shTACSTD2_DDP"]
        labels = ["shNC", "shTACSTD2", "shNC+DDP", "shTACSTD2+DDP"]
        data = []
        if isinstance(gene_or_score, str) and gene_or_score in logx.index:
            for g in order:
                data.append(logx.loc[gene_or_score, GROUPS[g]].to_numpy())
        else:
            for g in order:
                data.append(gene_or_score.loc[GROUPS[g]].to_numpy())
        ax.boxplot(data, tick_labels=labels, widths=0.55)
        for i, d in enumerate(data, start=1):
            ax.scatter(np.full_like(d, i, dtype=float) + np.random.default_rng(0).uniform(-0.08, 0.08, len(d)), d, s=18, c="k", zorder=3)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=20)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))
    if "TACSTD2" in logx.index:
        grouped_boxes(axes[0], "TACSTD2", "log2(FPKM+1)", "TACSTD2 (knockdown QC)")
    if "CLDN4" in logx.index:
        grouped_boxes(axes[1], "CLDN4", "log2(FPKM+1)", "CLDN4")
    fig.tight_layout()
    fig.savefig(FIGS / "fig1_TACSTD2_CLDN4.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
    grouped_boxes(axes[0], score_mat["C4_IFN_MHCI"], "z-mean", "C4 IFN/MHC-I panel")
    grouped_boxes(axes[1], score_mat["APM"], "z-mean", "APM")
    grouped_boxes(axes[2], score_mat["IFN"], "z-mean", "IFN / ISG")
    fig.tight_layout()
    fig.savefig(FIGS / "fig2_IFN_APM_scores.png")
    plt.close(fig)

    # C4 panel gene heatmap of log2FC.
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    panel = [g for g in C4_IFN_MHCI + ["CLDN4", "TACSTD2"] if g in de_tables["KD_noDDP"]["gene"].values]
    mat = []
    for cname in ["KD_noDDP", "KD_DDP", "DDP_in_NC", "DDP_in_KD"]:
        de = de_tables[cname].set_index("gene")
        mat.append(de.loc[panel, "log2FC"].to_numpy())
    mat = np.vstack(mat)
    vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(panel)))
    ax.set_xticklabels(panel, rotation=45, ha="right")
    ax.set_yticks(range(4))
    ax.set_yticklabels(["KD no DDP", "KD + DDP", "DDP in NC", "DDP in KD"])
    ax.set_title("log2FC (treat − ctrl)")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(FIGS / "fig3_C4_panel_log2FC.png")
    plt.close(fig)

    print("TACSTD2", key_stats["TACSTD2"])
    print("CLDN4", key_stats["CLDN4"])
    print("geneset KD_noDDP", json.dumps(gs_summary["KD_noDDP"], indent=2))
    print("verdict", json.dumps(key_stats["analog_verdict"], indent=2))
    print("wrote", TABLES)
    print("wrote", FIGS)


if __name__ == "__main__":
    main()
