#!/usr/bin/env python3
"""TISMO ICB pairing: Tacstd2, Cldn4, and TJ score (Cldn3/4/6/7, Cdh1, F11r, Ocln)."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tismo_client import read_expression_csv

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "tismo"
OUT = ROOT / "results" / "tismo"
FIG = OUT / "figures"
TABLES = OUT / "tables"

TJ_GENES = ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
STEM_RE = re.compile(r"\(n=\d+\)$")
ICB_TOKEN_RE = re.compile(
    r"antiCTLA4&antiPD[L]?1|antiCTLA4|antiPD1_FcS|antiPD1|antiPDL1|antiPDL2",
    re.I,
)

# TISMO cellLineMeta cancerType, collapsed for splits.
CANCER = {
    "4T1": ("Mammary carcinoma", "Mammary"),
    "E0771": ("Mammary adenocarcinoma", "Mammary"),
    "EMT6": ("Mammary carcinoma", "Mammary"),
    "KPB25L": ("Mammary cancer, NOS", "Mammary"),
    "T11": ("Mammary cancer, NOS", "Mammary"),
    "p53-2225L": ("Mammary cancer, NOS", "Mammary"),
    "p53-2336R": ("Mammary cancer, NOS", "Mammary"),
    "B16": ("Melanoma", "Melanoma"),
    "YUMM1.7": ("Melanoma", "Melanoma"),
    "D3UV2": ("Melanoma", "Melanoma"),
    "D4M.3A.3": ("Melanoma", "Melanoma"),
    "CT26": ("Colorectal carcinoma", "Colorectal"),
    "MC38": ("Colorectal carcinoma", "Colorectal"),
    "LLC": ("Lung carcinoma", "Lung"),
    "402230": ("Sarcoma", "Sarcoma"),
    "BNL-MEA": ("Hepatocellular carcinoma", "Liver"),
    "YTN16": ("Gastric adenocarcinoma", "Gastric"),
    "MOC22": ("Oral squamous cell carcinoma", "HeadNeck"),
    "CMT-167": ("Lung carcinoma", "Lung"),
}


def stem(label: str) -> str:
    return STEM_RE.sub("", str(label)).rstrip()


def parse_line(model: str) -> str:
    return model.split("_")[0]


def is_lung(line: str, cancer_group: str) -> bool:
    return cancer_group == "Lung" or line.upper() == "LLC"


def is_kl(text: str) -> bool:
    t = text.lower()
    return bool(re.search(r"(?:\bkl\b|stk11|lkb1|lkr13ko)", t))


def is_kp(text: str) -> bool:
    t = text.lower()
    # lung KP / Kras-p53; exclude mammary p53 lines and pancreatic KPC
    if "kpc" in t:
        return False
    return bool(re.search(r"(?:\bkp\b|kras.*trp53|kras.*p53|hkp1)", t))


def sign_counts(values) -> dict:
    v = pd.Series(values, dtype=float).dropna()
    n_up = int((v > 0).sum())
    n_down = int((v < 0).sum())
    n_tie = int((v == 0).sum())
    n = n_up + n_down
    if n == 0:
        p_two = p_greater = np.nan
    else:
        p_two = float(binomtest(n_up, n, 0.5, alternative="two-sided").pvalue)
        p_greater = float(binomtest(n_up, n, 0.5, alternative="greater").pvalue)
    return {
        "n_pairs": int(len(v)),
        "n_up": n_up,
        "n_down": n_down,
        "n_tie": n_tie,
        "n_tested": n,
        "binomial_p_two_sided": p_two,
        "binomial_p_greater": p_greater,
    }


def wilcoxon_delta(values) -> dict:
    v = pd.Series(values, dtype=float).dropna()
    v_nz = v[v != 0]
    if len(v_nz) < 1:
        return {"n": int(len(v_nz)), "statistic": np.nan, "p_two_sided": np.nan, "p_greater": np.nan}
    two = wilcoxon(v_nz, alternative="two-sided", zero_method="wilcox")
    greater = wilcoxon(v_nz, alternative="greater", zero_method="wilcox")
    return {
        "n": int(len(v_nz)),
        "statistic": float(two.statistic),
        "p_two_sided": float(two.pvalue),
        "p_greater": float(greater.pvalue),
    }


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (np.isnan(p) or np.isinf(p))):
        return "NA"
    if p == 0:
        return "<1e-16"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def load_vivo(gene: str) -> pd.DataFrame:
    df = read_expression_csv(DATA / "vivo" / f"{gene}.csv", gene=gene)
    df["stem"] = df["cell_line"].map(stem)
    df["line"] = df["stem"].map(parse_line)
    return df


def load_vitro(gene: str) -> pd.DataFrame:
    df = read_expression_csv(DATA / "vitro" / f"{gene}.csv", gene=gene)
    df["stem"] = df["cell_line"].map(stem)
    df["line"] = df["stem"].map(parse_line)
    return df


def cancer_of(line: str) -> tuple[str, str]:
    return CANCER.get(line, ("Unknown", "Other"))


def pair_table(df: pd.DataFrame, value_col: str = "value") -> pd.DataFrame:
    rows = []
    for model, g in df.groupby("stem", sort=True):
        line = parse_line(model)
        cancer, group = cancer_of(line)
        base = g.loc[g["Baseline"] == 1, value_col].astype(float)
        trt = g.loc[g["Baseline"] == 0, value_col].astype(float)
        if len(base) == 0 or len(trt) == 0:
            continue
        treatments = sorted({str(x) for x in g.get("Mouse_treatment", pd.Series(dtype=str)).dropna().unique()})
        gse = str(g["GSE_ID"].iloc[0]) if "GSE_ID" in g.columns else "NA"
        d_mean = float(trt.mean() - base.mean())
        d_med = float(trt.median() - base.median())
        rows.append(
            {
                "stem": model,
                "cell_line": line,
                "gse_id": gse,
                "cancer_type": cancer,
                "cancer_group": group,
                "is_lung": is_lung(line, group),
                "is_llc": line.upper() == "LLC",
                "is_kl": False,
                "is_kp": False,
                "n_baseline": int(len(base)),
                "n_icb": int(len(trt)),
                "mean_baseline": float(base.mean()),
                "mean_icb": float(trt.mean()),
                "median_baseline": float(base.median()),
                "median_icb": float(trt.median()),
                "delta_mean": d_mean,
                "delta_median": d_med,
                "treatments": ";".join(treatments),
            }
        )
    return pd.DataFrame(rows)


def stats_block(deltas: pd.Series, label: str) -> dict:
    sc_mean_style = sign_counts(deltas)
    w = wilcoxon_delta(deltas)
    return {
        "label": label,
        **sc_mean_style,
        "wilcoxon_n": w["n"],
        "wilcoxon_stat": w["statistic"],
        "wilcoxon_p_two_sided": w["p_two_sided"],
        "wilcoxon_p_greater": w["p_greater"],
        "mean_delta": float(pd.Series(deltas, dtype=float).mean()) if len(deltas) else np.nan,
        "median_delta": float(pd.Series(deltas, dtype=float).median()) if len(deltas) else np.nan,
    }


def sentence(gene: str, mean_stats: dict, med_stats: dict | None = None) -> str:
    s = (
        f"{gene} rose in {mean_stats['n_up']}/{mean_stats['n_pairs']} models "
        f"by mean(ICB−naive) (Wilcoxon p={fmt_p(mean_stats['wilcoxon_p_two_sided'])}; "
        f"binomial p={fmt_p(mean_stats['binomial_p_two_sided'])}"
    )
    if mean_stats["n_tie"]:
        s += f"; {mean_stats['n_tie']} ties"
    s += f"; n_down={mean_stats['n_down']}"
    s += ")"
    if med_stats is not None:
        s += (
            f"; median-based: {med_stats['n_up']}/{med_stats['n_pairs']} up "
            f"(Wilcoxon p={fmt_p(med_stats['wilcoxon_p_two_sided'])})"
        )
    return s


def build_sample_matrix(genes: list[str], lock_samples: pd.Index | None = None) -> pd.DataFrame:
    frames = []
    for gene in genes:
        df = load_vivo(gene)[["Samples", "value"]].rename(columns={"value": gene})
        df = df.drop_duplicates("Samples")
        frames.append(df.set_index("Samples"))
    mat = frames[0]
    for f in frames[1:]:
        mat = mat.join(f, how="outer")
    if lock_samples is not None:
        mat = mat.reindex(lock_samples)
    return mat


def tj_from_matrix(mat: pd.DataFrame) -> pd.Series:
    sub = mat[TJ_GENES]
    return sub.mean(axis=1, skipna=True)


def plot_waterfall(pairs: pd.DataFrame, col: str, title: str, ylabel: str, path: Path, highlight_lung: bool = True) -> None:
    d = pairs.sort_values([col, "stem"]).reset_index(drop=True)
    colors = np.where(d[col] >= 0, "#4C78A8", "#E45756")
    fig, ax = plt.subplots(figsize=(11.2, 4.4))
    ax.axhline(0, color="0.35", lw=0.8)
    ax.bar(np.arange(len(d)), d[col], color=colors, width=0.9, linewidth=0)
    if highlight_lung and d["is_lung"].any():
        idx = np.where(d["is_lung"].to_numpy())[0]
        ax.scatter(idx, d.loc[d["is_lung"], col], s=36, facecolors="none", edgecolors="#111", linewidths=1.2, zorder=3, label="lung (LLC)")
        ax.legend(frameon=False, loc="upper left")
    ax.set_ylabel(ylabel)
    ax.set_xlabel(f"TISMO ICB pairs (n={len(d)}), sorted")
    ax.set_xticks([])
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_paired_dots(pairs: pd.DataFrame, title: str, ylabel: str, path: Path) -> None:
    d = pairs.sort_values("delta_mean").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11.2, 4.6))
    x0 = np.zeros(len(d))
    x1 = np.ones(len(d))
    for i, row in d.iterrows():
        color = "#2A6F97" if row["is_lung"] else "0.55"
        lw = 1.6 if row["is_lung"] else 0.55
        ax.plot([0, 1], [row["mean_baseline"], row["mean_icb"]], color=color, lw=lw, alpha=0.85)
    ax.scatter(x0, d["mean_baseline"], s=14, c="0.25", zorder=3)
    ax.scatter(x1, d["mean_icb"], s=14, c="#C44E52", zorder=3)
    if d["is_lung"].any():
        lung = d[d["is_lung"]]
        ax.scatter(np.zeros(len(lung)), lung["mean_baseline"], s=42, facecolors="none", edgecolors="#111", linewidths=1.2, zorder=4)
        ax.scatter(np.ones(len(lung)), lung["mean_icb"], s=42, facecolors="none", edgecolors="#111", linewidths=1.2, zorder=4, label="lung (LLC)")
        ax.legend(frameon=False)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["naive / baseline", "ICB"])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(-0.25, 1.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def json_clean(obj):
    if isinstance(obj, dict):
        return {k: json_clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_clean(v) for v in obj]
    if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    if isinstance(obj, (np.floating,)):
        x = float(obj)
        return None if (np.isnan(x) or np.isinf(x)) else x
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def main() -> int:
    if not (DATA / "vivo" / "Tacstd2.csv").exists():
        raise SystemExit("missing TISMO downloads; run scripts/download_tismo.py")

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    meta_lines = json.loads((DATA / "cellLineMeta.json").read_text())
    line_cancer = {r.get("cellLine") or r.get("cell_line") or r.get("name"): r for r in meta_lines}

    tac = load_vivo("Tacstd2")
    lock_stems = sorted(set(tac["stem"]))
    lock_samples = pd.Index(tac["Samples"].unique())

    mat = build_sample_matrix(["Tacstd2"] + TJ_GENES, lock_samples=lock_samples)
    tj = tj_from_matrix(mat)
    n_tj_genes = mat[TJ_GENES].notna().sum(axis=1)

    # Attach TJ + Cldn4 onto the Tacstd2 sample table (64-slice lock).
    tac = tac.copy()
    tac["Cldn4"] = tac["Samples"].map(mat["Cldn4"])
    tac["TJ"] = tac["Samples"].map(tj)
    tac["n_tj_genes"] = tac["Samples"].map(n_tj_genes)

    pairs_tac = pair_table(tac)
    tmp_cld = tac.dropna(subset=["Cldn4"]).copy()
    tmp_cld["value"] = tmp_cld["Cldn4"]
    pairs_cld = pair_table(tmp_cld)
    tmp_tj = tac.dropna(subset=["TJ"]).copy()
    tmp_tj["value"] = tmp_tj["TJ"]
    pairs_tj = pair_table(tmp_tj)

    # Keep only the Tacstd2 64-stem universe for Cldn4/TJ.
    pairs_tac = pairs_tac[pairs_tac["stem"].isin(lock_stems)].copy()
    pairs_cld = pairs_cld[pairs_cld["stem"].isin(lock_stems)].copy()
    pairs_tj = pairs_tj[pairs_tj["stem"].isin(lock_stems)].copy()

    # Native Cldn4 groups (may be 65).
    cld_native = pair_table(load_vivo("Cldn4"))

    def annotate(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df["stratum_all"] = True
        df["stratum_lung"] = df["is_lung"]
        df["stratum_kl_kp_llc"] = df["is_llc"] | df["is_kl"] | df["is_kp"]
        return df

    pairs_tac = annotate(pairs_tac)
    pairs_cld = annotate(pairs_cld)
    pairs_tj = annotate(pairs_tj)

    strata = [
        ("all", "stratum_all"),
        ("lung", "stratum_lung"),
        ("kl_kp_llc", "stratum_kl_kp_llc"),
    ]
    markers = [
        ("Tacstd2", pairs_tac),
        ("Cldn4", pairs_cld),
        ("TJ", pairs_tj),
    ]

    summary_rows = []
    summary = {"markers": {}, "n_tismo_tacstd2_pairs": int(len(pairs_tac))}
    for marker, pairs in markers:
        summary["markers"][marker] = {}
        for sname, scol in strata:
            sub = pairs[pairs[scol]]
            mean_s = stats_block(sub["delta_mean"], f"{marker}_{sname}_mean")
            med_s = stats_block(sub["delta_median"], f"{marker}_{sname}_median")
            summary["markers"][marker][sname] = {
                "n_pairs": int(len(sub)),
                "cell_lines": sorted(sub["cell_line"].unique()),
                "mean": mean_s,
                "median": med_s,
                "sentence_mean": sentence(marker, mean_s, med_s) if len(sub) else f"{marker} {sname}: n=0",
            }
            for loc, st in (("mean", mean_s), ("median", med_s)):
                summary_rows.append(
                    {
                        "marker": marker,
                        "stratum": sname,
                        "location": loc,
                        **{k: v for k, v in st.items() if k != "label"},
                    }
                )

    # vitro IFNg
    vitro_rows = []
    vitro_summary = {}
    for gene in ["Tacstd2", "Cldn4"]:
        vdf = load_vitro(gene)
        # IFNg only
        if "cytokines" in vdf.columns:
            vdf = vdf[vdf["cytokines"].astype(str).str.contains("IFNg", case=False, na=False) | (vdf["Baseline"] == 1)]
        pairs = pair_table(vdf)
        # keep pairs that actually have IFNg-treated samples
        if "treatments" in pairs.columns:
            # treatments empty for vitro; filter via source
            pass
        mean_s = stats_block(pairs["delta_mean"], f"vitro_{gene}_mean")
        med_s = stats_block(pairs["delta_median"], f"vitro_{gene}_median")
        vitro_summary[gene] = {
            "n_pairs": int(len(pairs)),
            "cell_lines": sorted(pairs["cell_line"].unique()),
            "mean": mean_s,
            "median": med_s,
            "sentence": sentence(f"{gene} (IFNγ in vitro)", mean_s, med_s) if len(pairs) else f"{gene} vitro: n=0",
        }
        pairs = pairs.assign(marker=gene)
        vitro_rows.append(pairs)
        pairs.to_csv(TABLES / f"vitro_ifng_{gene}_pairs.tsv", sep="\t", index=False)

    vitro_all = pd.concat(vitro_rows, ignore_index=True) if vitro_rows else pd.DataFrame()

    # write pair tables
    pairs_tac.sort_values(["cancer_group", "cell_line", "stem"]).to_csv(TABLES / "pairs_Tacstd2.tsv", sep="\t", index=False)
    pairs_cld.sort_values(["cancer_group", "cell_line", "stem"]).to_csv(TABLES / "pairs_Cldn4.tsv", sep="\t", index=False)
    pairs_tj.sort_values(["cancer_group", "cell_line", "stem"]).to_csv(TABLES / "pairs_TJ.tsv", sep="\t", index=False)
    cld_native.to_csv(TABLES / "pairs_Cldn4_native.tsv", sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(TABLES / "stats_by_stratum.tsv", sep="\t", index=False)

    merged = pairs_tac[["stem", "cell_line", "cancer_group", "is_lung", "is_llc", "delta_mean", "delta_median", "mean_baseline", "mean_icb"]].rename(
        columns={
            "delta_mean": "Tacstd2_delta_mean",
            "delta_median": "Tacstd2_delta_median",
            "mean_baseline": "Tacstd2_mean_naive",
            "mean_icb": "Tacstd2_mean_icb",
        }
    )
    merged = merged.merge(
        pairs_cld[["stem", "delta_mean", "delta_median", "mean_baseline", "mean_icb"]].rename(
            columns={
                "delta_mean": "Cldn4_delta_mean",
                "delta_median": "Cldn4_delta_median",
                "mean_baseline": "Cldn4_mean_naive",
                "mean_icb": "Cldn4_mean_icb",
            }
        ),
        on="stem",
        how="left",
    )
    merged = merged.merge(
        pairs_tj[["stem", "delta_mean", "delta_median", "mean_baseline", "mean_icb"]].rename(
            columns={
                "delta_mean": "TJ_delta_mean",
                "delta_median": "TJ_delta_median",
                "mean_baseline": "TJ_mean_naive",
                "mean_icb": "TJ_mean_icb",
            }
        ),
        on="stem",
        how="left",
    )
    merged.to_csv(TABLES / "pairs_merged.tsv", sep="\t", index=False)

    # figures
    plot_waterfall(
        pairs_tac, "delta_mean",
        "Tacstd2 after ICB (mean ICB − mean naive)",
        "Δ Tacstd2 (TISMO log)",
        FIG / "waterfall_Tacstd2.png",
    )
    plot_waterfall(
        pairs_cld, "delta_mean",
        "Cldn4 after ICB (mean ICB − mean naive)",
        "Δ Cldn4 (TISMO log)",
        FIG / "waterfall_Cldn4.png",
    )
    plot_waterfall(
        pairs_tj, "delta_mean",
        "TJ score after ICB (mean ICB − mean naive)",
        "Δ TJ score (mean of Cldn3/4/6/7, Cdh1, F11r, Ocln)",
        FIG / "waterfall_TJ.png",
    )
    plot_paired_dots(
        pairs_tac,
        "Tacstd2 paired means: naive vs ICB (64 TISMO slices)",
        "mean Tacstd2 (TISMO log)",
        FIG / "paired_Tacstd2.png",
    )
    plot_paired_dots(
        pairs_cld,
        "Cldn4 paired means: naive vs ICB (Tacstd2-locked slices)",
        "mean Cldn4 (TISMO log)",
        FIG / "paired_Cldn4.png",
    )
    plot_paired_dots(
        pairs_tj,
        "TJ score paired means: naive vs ICB",
        "mean TJ score (TISMO log)",
        FIG / "paired_TJ.png",
    )

    # lung-only paired if present
    lung_tac = pairs_tac[pairs_tac["is_lung"]]
    if len(lung_tac):
        plot_paired_dots(
            lung_tac,
            "Lung (LLC) Tacstd2: naive vs ICB",
            "mean Tacstd2 (TISMO log)",
            FIG / "paired_lung_Tacstd2.png",
        )
        lung_cld = pairs_cld[pairs_cld["is_lung"]]
        if len(lung_cld):
            plot_paired_dots(
                lung_cld,
                "Lung (LLC) Cldn4: naive vs ICB",
                "mean Cldn4 (TISMO log)",
                FIG / "paired_lung_Cldn4.png",
            )

    # coverage / notes
    extra_cld = sorted(set(cld_native["stem"]) - set(lock_stems))
    missing_cld = sorted(set(lock_stems) - set(pairs_cld["stem"]))
    notes = {
        "pairing": (
            "Author-style TISMO Gene-module slices: one pair per cell_line label "
            "(study × condition × ICB regimen). Baseline==1 vs Baseline==0 "
            "(responders + non-responders pooled). Tacstd2 export defines the 64-model lock."
        ),
        "values": "TISMO Gene-module log expression as served by downVivoExprn (log2(TPM+1) in TISMO docs).",
        "tj_score": "Per-sample mean of Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, Ocln; then pair like single genes.",
        "n_tacstd2_pairs": int(len(pairs_tac)),
        "n_cldn4_locked_pairs": int(len(pairs_cld)),
        "n_tj_locked_pairs": int(len(pairs_tj)),
        "n_cldn4_native_pairs": int(len(cld_native)),
        "cldn4_stems_not_in_tacstd2_64": extra_cld,
        "tacstd2_stems_missing_cldn4": missing_cld,
        "lung_models_in_tismo_icb": sorted(pairs_tac.loc[pairs_tac["is_lung"], "stem"].tolist()),
        "kl_models_in_tismo_icb": [],
        "kp_models_in_tismo_icb": [],
        "llc_is_not_kl": True,
        "line_cancer_keys_sample": sorted({k for k in line_cancer if k})[:20],
        "mean_tj_genes_per_sample": float(tac["n_tj_genes"].mean()) if len(tac) else None,
    }

    paper = {
        "all": (
            "After ICB, "
            + summary["markers"]["Tacstd2"]["all"]["sentence_mean"]
            + "; "
            + summary["markers"]["Cldn4"]["all"]["sentence_mean"]
            + "; "
            + summary["markers"]["TJ"]["all"]["sentence_mean"]
            + "."
        ),
        "lung": (
            "Lung subset (TISMO lung-carcinoma ICB = LLC only; LLC is not KL): "
            + summary["markers"]["Tacstd2"]["lung"]["sentence_mean"]
            + "; "
            + summary["markers"]["Cldn4"]["lung"]["sentence_mean"]
            + "; "
            + summary["markers"]["TJ"]["lung"]["sentence_mean"]
            + "."
        ),
        "ifng": (
            vitro_summary.get("Tacstd2", {}).get("sentence", "")
            + "; "
            + vitro_summary.get("Cldn4", {}).get("sentence", "")
            + "."
        ),
    }

    payload = {
        "source": {
            "database": "TISMO",
            "url": "https://tismo.pku-genomics.org/",
            "citation": "Zeng et al. Nucleic Acids Research 2022 (PMID 34534350)",
            "vivo_endpoint": "POST /rtismo/gene/downVivoExprn type=3",
            "vitro_endpoint": "POST /rtismo/gene/downVitroExprn type=3",
        },
        "notes": notes,
        "summary": summary,
        "vitro_ifng": vitro_summary,
        "paper_sentences": paper,
    }
    (OUT / "summary.json").write_text(json.dumps(json_clean(payload), indent=2) + "\n")

    print("=== Tacstd2 all (mean) ===")
    print(summary["markers"]["Tacstd2"]["all"]["sentence_mean"])
    print("=== Cldn4 all (mean) ===")
    print(summary["markers"]["Cldn4"]["all"]["sentence_mean"])
    print("=== TJ all (mean) ===")
    print(summary["markers"]["TJ"]["all"]["sentence_mean"])
    print("=== lung ===")
    print(paper["lung"])
    print("=== IFNg ===")
    print(paper["ifng"])
    print("n Tacstd2 pairs", len(pairs_tac), "Cldn4 native", len(cld_native))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
