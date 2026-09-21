#!/usr/bin/env python3
"""IFN-limb sensitivity and NHEJ/STING scores for GSE137244.

Does not recompute Tacstd2, Cldn4, or the TJ locks.

IFN catalog (pre-specified, not fit to a delta):
- Mouse Hallmark IFN-alpha and IFN-gamma from MSigDB 2022.1, 2023.1, 2023.2, 2024.1, 2025.1
- Reactome 2024.1 mouse: IFN alpha/beta, IFN gamma, IFN signaling, antiviral ISG,
  ISG15, OAS, DDX58/IFIH1
- COMPACT_ISG from the prior KL vs KP gene list

Thesis-aligned cold = every KL library lower than every KP library on the mean
log2(FPKM+1) score. A panel is a robust hard cold only if that complete split
also holds in every leave-one-library contrast and after the previously reported
KP block collapse (B6AL10-1/2/4/5 averaged, B6AL10-3 kept). Otherwise the call
stays soft.

NHEJ and STING are scored on the same scale. They are not part of the IFN hunt.
"""

from __future__ import annotations

import json
import math
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import analyze as A

CACHE = A.CACHE
TABLES = A.TABLES
FIGS = A.FIGS
GENESETS = A.GENESETS
MSIG = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"

# Symbol renames between this 2019 mm10 matrix and later MSigDB symbols.
# Applied only when the listed symbol is absent and the alias is present.
ALIASES = {
    "Sting1": "Tmem173",
    "Cgas": "Mb21d1",
    "Rigi": "Ddx58",
    "Wars1": "Wars",
    "Tmt1b": "Mettl7b",
    "H2ax": "H2afx",
}

# Enzymatic NHEJ subunits. Histone genes in the Reactome set use symbols this
# matrix does not carry; the core is the repair complex, not a histone average.
CORE_NHEJ = [
    "Xrcc6",
    "Xrcc5",
    "Prkdc",
    "Xrcc4",
    "Lig4",
    "Nhej1",
    "Dclre1c",
    "Paxx",
    "Poll",
    "Polm",
    "Trp53bp1",
]
# cGAS-STING core on the symbols present in this matrix.
CORE_CGAS_STING = ["Mb21d1", "Tmem173", "Tbk1", "Irf3", "Ikbke", "Ifnb1"]

# r >= 0.98 KP block from the library-correlation analysis. Not re-fit here.
KP_BLOCK = [
    "B6AL10-1-RNA",
    "B6AL10-2-RNA",
    "B6AL10-4-RNA",
    "B6AL10-5-RNA",
]
KP_SINGLE = "B6AL10-3-RNA"

HALLMARK_VERSIONS = ("2022.1", "2023.1", "2023.2", "2024.1", "2025.1")
IFN_HALLMARK = (
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
)
REACTOME_IFN = (
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
    "REACTOME_INTERFERON_GAMMA_SIGNALING",
    "REACTOME_INTERFERON_SIGNALING",
    "REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES",
    "REACTOME_ISG15_ANTIVIRAL_MECHANISM",
    "REACTOME_OAS_ANTIVIRAL_RESPONSE",
    "REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA",
)
REACTOME_OTHER = (
    "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES",
    "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
    "REACTOME_HDR_THROUGH_MMEJ_ALT_NHEJ",
)


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return
    urllib.request.urlretrieve(url, dest)


def resolve(genes: list[str], index: set[str]) -> tuple[list[str], list[str], list[str]]:
    """Return (used symbols, missing, alias notes)."""
    used: list[str] = []
    missing: list[str] = []
    notes: list[str] = []
    for gene in genes:
        if gene in index:
            pick = gene
        elif gene in ALIASES and ALIASES[gene] in index:
            pick = ALIASES[gene]
            notes.append(f"{gene}->{pick}")
        else:
            missing.append(gene)
            continue
        if pick not in used:
            used.append(pick)
    return used, missing, notes


def score_vector(log_expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    return A.mean_score(log_expr, genes, A.KP_LIBS + A.KL_LIBS)


def contrast(scores: pd.Series, kl: list[str], kp: list[str]) -> dict:
    stats = A.exact_mw(scores.reindex(kl), scores.reindex(kp))
    kl_max = float(scores.reindex(kl).max())
    kp_min = float(scores.reindex(kp).min())
    stats["complete_kl_lower"] = bool(kl_max < kp_min)
    stats["kl_max"] = kl_max
    stats["kp_min"] = kp_min
    stats["n_kl"] = len(kl)
    stats["n_kp"] = len(kp)
    stats["exact_floor_p"] = 2 / math.comb(len(kl) + len(kp), len(kl))
    return stats


def style_ax(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    A.ensure_fpkm()
    fpkm, _meta = A.load_expression()
    libs = A.KP_LIBS + A.KL_LIBS
    log_expr = np.log2(fpkm[libs] + 1)
    index = set(fpkm.index)

    gmt_dir = CACHE / "msig"
    hallmark_paths = {}
    for ver in HALLMARK_VERSIONS:
        if ver == "2024.1":
            hallmark_paths[ver] = GENESETS / "mh.all.v2024.1.Mm.symbols.gmt"
            continue
        dest = gmt_dir / f"mh.all.v{ver}.Mm.symbols.gmt"
        download(f"{MSIG}/{ver}.Mm/mh.all.v{ver}.Mm.symbols.gmt", dest)
        hallmark_paths[ver] = dest
    reactome_path = gmt_dir / "m2.cp.reactome.v2024.1.Mm.symbols.gmt"
    download(f"{MSIG}/2024.1.Mm/m2.cp.reactome.v2024.1.Mm.symbols.gmt", reactome_path)
    reactome = A.load_gmt(reactome_path)
    compact = A.load_gmt(GENESETS / "epithelial_and_isg.gmt")["COMPACT_ISG"]

    panels: dict[str, dict] = {}
    for ver, path in hallmark_paths.items():
        gmt = A.load_gmt(path)
        tag = ver.replace(".", "_")
        for name in IFN_HALLMARK:
            key = f"{name}__{tag}"
            panels[key] = {"family": "IFN", "source": f"MSigDB mouse hallmark {ver}", "raw": gmt[name]}
    for name in REACTOME_IFN:
        panels[name] = {"family": "IFN", "source": "MSigDB Reactome 2024.1 mouse", "raw": reactome[name]}
    panels["COMPACT_ISG"] = {
        "family": "IFN",
        "source": "prior public KL vs KP compact ISG",
        "raw": compact,
    }
    for name in REACTOME_OTHER:
        family = "STING" if "STING" in name else "NHEJ"
        panels[name] = {"family": family, "source": "MSigDB Reactome 2024.1 mouse", "raw": reactome[name]}
    panels["CORE_NHEJ"] = {"family": "NHEJ", "source": "enzymatic NHEJ subunits", "raw": CORE_NHEJ}
    panels["CORE_CGAS_STING"] = {"family": "STING", "source": "cGAS-STING core, matrix symbols", "raw": CORE_CGAS_STING}

    gene_rows = []
    resolved = {}
    for key, spec in panels.items():
        used, missing, notes = resolve(spec["raw"], index)
        if len(used) < 3:
            raise SystemExit(f"{key} mapped only {len(used)} genes")
        resolved[key] = used
        gene_rows.append(
            {
                "panel": key,
                "family": spec["family"],
                "source": spec["source"],
                "n_listed": len(spec["raw"]),
                "n_used": len(used),
                "n_missing": len(missing),
                "missing": ",".join(missing),
                "aliases": ",".join(notes),
                "genes_used": ",".join(used),
            }
        )
    membership = pd.DataFrame(gene_rows)
    membership.to_csv(TABLES / "ifn_nhej_sting_panels.tsv", sep="\t", index=False)

    # Mean log2 scores, primary 5 vs 5.
    mean_scores = {key: score_vector(log_expr, genes) for key, genes in resolved.items()}
    primary_rows = []
    for key, scores in mean_scores.items():
        stats = contrast(scores, A.KL_LIBS, A.KP_LIBS)
        primary_rows.append(
            {
                "panel": key,
                "family": panels[key]["family"],
                "metric": "mean_log2",
                "n_genes": len(resolved[key]),
                **stats,
                **{A.SHORT[lib]: float(scores.loc[lib]) for lib in libs},
            }
        )

    print("ssGSEA ...", flush=True)
    # min_size=5 so the 5-gene OAS set and the 6–7 gene STING sets are scored.
    # Larger Hallmark and Reactome sets are unaffected.
    import gseapy as gp

    ss = gp.ssgsea(
        data=log_expr[A.KP_LIBS + A.KL_LIBS],
        gene_sets=resolved,
        min_size=5,
        max_size=500,
        sample_norm_method="rank",
        correl_norm_type="rank",
        permutation_num=0,
        weight=0.25,
        threads=4,
        seed=A.SEED,
        no_plot=True,
        outdir=None,
        verbose=False,
    ).res2d.copy()
    ss = ss.rename(columns={"Name": "library", "Term": "panel"})
    missing_terms = [key for key in resolved if key not in set(ss["panel"])]
    if missing_terms:
        raise SystemExit(
            "ssGSEA terms missing: "
            + ", ".join(missing_terms)
            + " | got "
            + ", ".join(map(str, ss["panel"].unique()[:8]))
        )
    ss_rows = []
    for key in resolved:
        sub = ss.loc[ss["panel"] == key].set_index("library")["ES"]
        stats = contrast(sub, A.KL_LIBS, A.KP_LIBS)
        ss_rows.append(
            {
                "panel": key,
                "family": panels[key]["family"],
                "metric": "ssgsea_es",
                "n_genes": len(resolved[key]),
                **stats,
                **{A.SHORT[lib]: float(sub.loc[lib]) for lib in libs},
            }
        )
    primary = pd.DataFrame(primary_rows + ss_rows)
    primary.to_csv(TABLES / "ifn_nhej_sting_primary.tsv", sep="\t", index=False)

    # Leave-one-library on the mean scores. ssGSEA ES does not depend on the other libraries.
    leave_rows = []
    for key, scores in mean_scores.items():
        for drop in libs:
            kl = [x for x in A.KL_LIBS if x != drop]
            kp = [x for x in A.KP_LIBS if x != drop]
            stats = contrast(scores, kl, kp)
            leave_rows.append(
                {
                    "panel": key,
                    "family": panels[key]["family"],
                    "dropped": A.SHORT[drop],
                    "dropped_arm": "KL" if drop in A.KL_LIBS else "KP",
                    **stats,
                }
            )
    leave = pd.DataFrame(leave_rows)
    leave.to_csv(TABLES / "ifn_leaveone.tsv", sep="\t", index=False)

    # KP block collapse: one mean for B6AL10-1/2/4/5, B6AL10-3 kept. 5 KL vs 2 KP units.
    block_rows = []
    for key, scores in mean_scores.items():
        block_mean = float(scores.reindex(KP_BLOCK).mean())
        units_kp = pd.Series({"KP_BLOCK": block_mean, "KP_B6AL10-3": float(scores.loc[KP_SINGLE])})
        units_kl = scores.reindex(A.KL_LIBS)
        stats = A.exact_mw(units_kl.to_numpy(), units_kp.to_numpy())
        kl_max = float(units_kl.max())
        kp_min = float(units_kp.min())
        n_kl_below_block = int((units_kl < block_mean).sum())
        block_rows.append(
            {
                "panel": key,
                "family": panels[key]["family"],
                "kp_block_mean": block_mean,
                "kp_b6al10_3": float(scores.loc[KP_SINGLE]),
                "kl_mean": float(units_kl.mean()),
                "kl_max": kl_max,
                "n_kl_below_block": n_kl_below_block,
                "delta_kl_minus_two_kp_units": stats["delta_mean"],
                "p_5kl_vs_2kp": stats["p"],
                "complete_kl_lower": bool(kl_max < kp_min),
                "exact_floor_p_5v2": 2 / math.comb(7, 5),
                **{A.SHORT[lib]: float(scores.loc[lib]) for lib in A.KL_LIBS},
            }
        )
    block = pd.DataFrame(block_rows)
    block.to_csv(TABLES / "ifn_kp_block.tsv", sep="\t", index=False)

    # Hunt. Hard cold = complete KL-lower on the primary mean. Robust = hard on
    # every leave-one and on the 5-vs-2 block contrast.
    mean_ifn = primary.loc[(primary["family"] == "IFN") & (primary["metric"] == "mean_log2")].copy()
    ss_ifn = primary.loc[(primary["family"] == "IFN") & (primary["metric"] == "ssgsea_es")].copy()

    def robust_mask(panel: str) -> dict:
        sub = leave.loc[leave["panel"] == panel]
        block_row = block.loc[block["panel"] == panel].iloc[0]
        return {
            "leaveone_all_complete_kl_lower": bool(sub["complete_kl_lower"].all()),
            "leaveone_n_complete": int(sub["complete_kl_lower"].sum()),
            "leaveone_n": int(len(sub)),
            "block_complete_kl_lower": bool(block_row["complete_kl_lower"]),
            "block_p": float(block_row["p_5kl_vs_2kp"]),
            "n_kl_below_block": int(block_row["n_kl_below_block"]),
        }

    hunt_rows = []
    for _, row in mean_ifn.iterrows():
        extra = robust_mask(row["panel"])
        ss_row = ss_ifn.loc[ss_ifn["panel"] == row["panel"]].iloc[0]
        hunt_rows.append(
            {
                "panel": row["panel"],
                "n_genes": int(row["n_genes"]),
                "mean_delta": float(row["delta_mean"]),
                "mean_p": float(row["p"]),
                "mean_complete_kl_lower": bool(row["complete_kl_lower"]),
                "mean_cliffs": float(row["cliffs_delta"]),
                "ssgsea_delta": float(ss_row["delta_mean"]),
                "ssgsea_p": float(ss_row["p"]),
                "ssgsea_complete_kl_lower": bool(ss_row["complete_kl_lower"]),
                "robust_hard_cold": bool(
                    row["complete_kl_lower"]
                    and extra["leaveone_all_complete_kl_lower"]
                    and extra["block_complete_kl_lower"]
                ),
                **extra,
            }
        )
    hunt = pd.DataFrame(hunt_rows).sort_values(["mean_complete_kl_lower", "mean_delta"], ascending=[False, True])
    hunt.to_csv(TABLES / "ifn_hunt.tsv", sep="\t", index=False)

    hard = hunt.loc[hunt["mean_complete_kl_lower"]]
    if len(hard):
        strongest = hard.sort_values("mean_delta").iloc[0]
        strongest_kind = "hard_on_primary_mean"
    else:
        strongest = hunt.sort_values("mean_delta").iloc[0]
        strongest_kind = "soft_most_negative_mean"
    ss_hard = hunt.loc[hunt["ssgsea_complete_kl_lower"]]
    n_robust = int(hunt["robust_hard_cold"].sum())
    final_call = "hard_ifn_cold" if n_robust else "soft_ifn"
    summary = {
        "final_call": final_call,
        "n_ifn_panels": int(len(hunt)),
        "n_primary_mean_hard_cold": int(len(hard)),
        "n_ssgsea_hard_cold": int(len(ss_hard)),
        "n_robust_hard_cold": n_robust,
        "strongest_kind": strongest_kind,
        "strongest_panel": strongest["panel"],
        "strongest_mean_delta": float(strongest["mean_delta"]),
        "strongest_mean_p": float(strongest["mean_p"]),
        "strongest_mean_complete_kl_lower": bool(strongest["mean_complete_kl_lower"]),
        "strongest_ssgsea_complete_kl_lower": bool(strongest["ssgsea_complete_kl_lower"]),
        "strongest_leaveone_n_complete": int(strongest["leaveone_n_complete"]),
        "strongest_block_complete_kl_lower": bool(strongest["block_complete_kl_lower"]),
        "ssgsea_hard_panels": ss_hard["panel"].tolist(),
        "rule": (
            "Robust hard IFN-cold requires complete KL-lower mean log2 on all 10 libraries, "
            "in every leave-one-library contrast, and after collapsing the r>=0.98 KP block "
            "to one unit beside B6AL10-3. Otherwise the call is soft."
        ),
    }
    (TABLES / "ifn_final_call.json").write_text(json.dumps(summary, indent=2))

    plot_leaveone(leave, hunt, strongest["panel"])
    plot_nhej_sting(mean_scores)
    print(json.dumps(summary, indent=2))
    print("Done", flush=True)


def plot_leaveone(leave: pd.DataFrame, hunt: pd.DataFrame, strongest: str) -> None:
    show = [
        "HALLMARK_INTERFERON_ALPHA_RESPONSE__2024_1",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE__2024_1",
        "COMPACT_ISG",
        strongest,
    ]
    # unique, keep order
    seen = []
    for name in show:
        if name not in seen and name in set(leave["panel"]):
            seen.append(name)
    fig, axes = plt.subplots(1, len(seen), figsize=(3.3 * len(seen), 4.6), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, panel in zip(axes, seen):
        sub = leave.loc[leave["panel"] == panel].copy()
        # stable library order
        order = [A.SHORT[x] for x in A.KP_LIBS + A.KL_LIBS]
        sub["dropped"] = pd.Categorical(sub["dropped"], order, ordered=True)
        sub = sub.sort_values("dropped")
        colors = ["#4C78A8" if arm == "KP" else "#E45756" for arm in sub["dropped_arm"]]
        y = np.arange(len(sub))
        ax.axvline(0, color="#888888", lw=0.6)
        ax.scatter(sub["delta_mean"], y, c=colors, s=36, zorder=3)
        for i, complete in enumerate(sub["complete_kl_lower"]):
            if complete:
                ax.scatter(
                    [sub["delta_mean"].iloc[i]],
                    [i],
                    s=80,
                    facecolors="none",
                    edgecolors="#222222",
                    linewidths=1.0,
                    zorder=4,
                )
        ax.set_yticks(y)
        ax.set_yticklabels(list(sub["dropped"]), fontsize=8)
        short = panel.replace("HALLMARK_INTERFERON_", "IFN ").replace("_RESPONSE__2024_1", " 2024.1")
        short = short.replace("REACTOME_", "").replace("__", " ")
        row = hunt.loc[hunt["panel"] == panel].iloc[0]
        ax.set_title(f"{short}\nprimary Δ={row['mean_delta']:+.2f}", fontsize=9)
        ax.set_xlabel("Δ when that library is dropped", fontsize=8)
        style_ax(ax)
    axes[0].set_ylabel("Library left out")
    fig.suptitle("Leave-one-library  ·  open ring = that drop makes KL completely lower", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "ifn_leaveone.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGS / "ifn_leaveone.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_nhej_sting(mean_scores: dict[str, pd.Series]) -> None:
    keys = [
        "CORE_NHEJ",
        "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
        "REACTOME_HDR_THROUGH_MMEJ_ALT_NHEJ",
        "CORE_CGAS_STING",
        "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES",
    ]
    labels = {
        "CORE_NHEJ": "NHEJ core",
        "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ": "Reactome NHEJ",
        "REACTOME_HDR_THROUGH_MMEJ_ALT_NHEJ": "Alt-NHEJ",
        "CORE_CGAS_STING": "cGAS-STING core",
        "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES": "Reactome STING",
    }
    fig, axes = plt.subplots(1, len(keys), figsize=(12.5, 4.4), sharey=True)
    ymap = {lib: i for i, lib in enumerate(A.KP_LIBS + A.KL_LIBS)}
    for ax, key in zip(axes, keys):
        scores = mean_scores[key]
        for lib in A.KP_LIBS + A.KL_LIBS:
            arm = "KP" if lib in A.KP_LIBS else "KL"
            ax.scatter(
                float(scores.loc[lib]),
                ymap[lib],
                c="#4C78A8" if arm == "KP" else "#E45756",
                s=28,
                zorder=3,
            )
        stats = contrast(scores, A.KL_LIBS, A.KP_LIBS)
        ax.set_title(f"{labels[key]}\nΔ={stats['delta_mean']:+.2f}  p={stats['p']:.3f}", fontsize=9)
        ax.set_xlabel("mean log2(FPKM+1)", fontsize=8)
        style_ax(ax)
    axes[0].set_yticks(list(ymap.values()))
    axes[0].set_yticklabels([A.SHORT[lib] for lib in ymap], fontsize=8)
    fig.suptitle("NHEJ and STING  ·  KL red, KP blue", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "nhej_sting_per_line.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGS / "nhej_sting_per_line.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
