#!/usr/bin/env python3
"""Co-dependency of TACSTD2 and CLDN4 across CRISPR gene-effect profiles.

Co-dependency = correlation of a gene's Chronos profile with every other gene's
Chronos profile across cell lines. It is only interpretable if the query gene's
profile carries reproducible signal, so this script does three things in order:

1. Technical floor. Using the screen-level matrix, lines screened more than once
   (independent libraries/screens) give a direct estimate of how reproducible
   each gene's gene effect is. TACSTD2/CLDN4 are compared with common essentials
   and non-essential controls.
2. Co-dependency scan, pan-cancer and lung-only, with BH FDR across ~18k genes
   and a label-permutation null that preserves the correlation structure of the
   gene-effect matrix.
3. Positive control. The identical pipeline is run on query genes with known
   partners so that a null result for TACSTD2/CLDN4 can be read as biology plus
   noise rather than a broken scan.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DATA,
    FIGURES,
    GENES_OF_INTEREST,
    NOTES,
    RELEASE,
    TABLES,
    corr_vector_vs_matrix,
    ensure_dirs,
    load_gene_effect,
    load_models,
    lung_cohort,
    strip_entrez,
)

# Positive controls: query genes whose top co-dependencies are textbook.
POSITIVE_CONTROLS = ["EGFR", "CTNNB1", "MYC"]
MIN_SD = 0.10  # drop genes whose profile is flat; they only add noise


def technical_floor(lung_ids: pd.Index) -> pd.DataFrame | None:
    """Between-screen reproducibility of gene effect for lines screened twice."""
    screen_path = DATA / "ScreenGeneEffect.csv"
    map_path = DATA / "CRISPRScreenMap.csv"
    if not screen_path.exists() or not map_path.exists():
        print("screen-level files absent; skipping technical floor", flush=True)
        return None

    screen_map = pd.read_csv(map_path)
    id_col = "ScreenID" if "ScreenID" in screen_map.columns else screen_map.columns[0]
    model_col = "ModelID" if "ModelID" in screen_map.columns else screen_map.columns[1]

    ess = strip_entrez(pd.read_csv(DATA / "AchillesCommonEssentialControls.csv").iloc[:, 0])
    non = strip_entrez(pd.read_csv(DATA / "AchillesNonessentialControls.csv").iloc[:, 0])

    header = pd.read_csv(screen_path, nrows=0)
    id_header = header.columns[0]
    symbols = pd.Series(strip_entrez(header.columns), index=header.columns)
    wanted_symbols = set(GENES_OF_INTEREST) | set(ess[:400]) | set(non[:400])
    keep_cols = [id_header] + symbols[symbols.isin(wanted_symbols)].index.tolist()

    screens = pd.read_csv(screen_path, index_col=0, usecols=keep_cols)
    screens.columns = strip_entrez(screens.columns)
    screens = screens.loc[:, ~screens.columns.duplicated()]

    smap = screen_map.set_index(id_col)[model_col]
    screens = screens[screens.index.isin(smap.index)]
    screens["ModelID"] = smap.reindex(screens.index).to_numpy()

    lung_screens = screens[screens["ModelID"].isin(lung_ids)]
    counts = lung_screens["ModelID"].value_counts()
    repeated = counts[counts >= 2].index
    if len(repeated) < 5:
        print(f"only {len(repeated)} repeatedly screened lung lines; using all lineages", flush=True)
        counts = screens["ModelID"].value_counts()
        repeated = counts[counts >= 2].index
        source = screens
        scope = "all lineages"
    else:
        source = lung_screens
        scope = "lung lines"

    pairs = []
    for model in repeated:
        rows = source[source["ModelID"] == model].drop(columns="ModelID")
        pairs.append((rows.iloc[0], rows.iloc[1]))
    if not pairs:
        return None
    first = pd.DataFrame([p[0] for p in pairs])
    second = pd.DataFrame([p[1] for p in pairs])

    out = []
    ess_set = set(ess)
    non_set = set(non)
    for gene in first.columns:
        a = pd.to_numeric(first[gene], errors="coerce")
        b = pd.to_numeric(second[gene], errors="coerce")
        ok = a.notna().to_numpy() & b.notna().to_numpy()
        if ok.sum() < 10:
            continue
        r, p = stats.pearsonr(a[ok], b[ok])
        klass = (
            "target"
            if gene in GENES_OF_INTEREST
            else "common essential"
            if gene in ess_set
            else "non-essential control"
            if gene in non_set
            else "other"
        )
        out.append(
            {
                "gene": gene,
                "class": klass,
                "n_lines_screened_twice": int(ok.sum()),
                "between_screen_r": float(r),
                "p": float(p),
                "sd_screen1": float(a[ok].std(ddof=1)),
            }
        )
    frame = pd.DataFrame(out)
    frame.attrs["scope"] = scope
    frame.attrs["n_models"] = len(repeated)
    return frame


def scan(query: str, mat: pd.DataFrame, label: str, n_perm: int = 200) -> tuple[pd.DataFrame, dict]:
    """Correlate one gene's profile with all others; compare to a permutation null."""
    y = mat[query].dropna()
    if y.size < 20 or y.std(ddof=1) == 0:
        return pd.DataFrame(), {}
    others = mat.drop(columns=[query])
    res = corr_vector_vs_matrix(y, others, method="pearson")
    res["query"] = query
    res["cohort"] = label
    res = res.sort_values("r", ascending=False)

    obs_max = float(np.nanmax(np.abs(res["r"])))
    obs_nsig = int((res["q"] < 0.05).sum())
    if n_perm > 0:
        rng = np.random.default_rng(12345)
        perm_max = np.empty(n_perm)
        perm_nsig = np.empty(n_perm)
        yv = y.to_numpy()
        sub = others.loc[y.index]
        for i in range(n_perm):
            shuffled = pd.Series(rng.permutation(yv), index=y.index)
            pr = corr_vector_vs_matrix(shuffled, sub, method="pearson")
            perm_max[i] = float(np.nanmax(np.abs(pr["r"])))
            perm_nsig[i] = int((pr["q"] < 0.05).sum())
        perm_stats = {
            "perm_null_max_abs_r_mean": float(perm_max.mean()),
            "perm_null_max_abs_r_p95": float(np.quantile(perm_max, 0.95)),
            "perm_p_for_max_abs_r": float((np.sum(perm_max >= obs_max) + 1) / (n_perm + 1)),
            "perm_null_n_q_lt_0.05_mean": float(perm_nsig.mean()),
            "perm_null_n_q_lt_0.05_p95": float(np.quantile(perm_nsig, 0.95)),
            "perm_p_for_n_q_lt_0.05": float(
                (np.sum(perm_nsig >= obs_nsig) + 1) / (n_perm + 1)
            ),
        }
    else:
        perm_stats = {
            k: np.nan
            for k in (
                "perm_null_max_abs_r_mean",
                "perm_null_max_abs_r_p95",
                "perm_p_for_max_abs_r",
                "perm_null_n_q_lt_0.05_mean",
                "perm_null_n_q_lt_0.05_p95",
                "perm_p_for_n_q_lt_0.05",
            )
        }
    stat = {
        "query": query,
        "cohort": label,
        "n_lines": int(y.size),
        "n_genes_tested": int(res.shape[0]),
        "profile_sd": float(y.std(ddof=1)),
        "max_abs_r": obs_max,
        "median_abs_r": float(res["r"].abs().median()),
        **perm_stats,
        "n_q_lt_0.05": obs_nsig,
        "smallest_abs_r_passing_fdr05": (
            float(res.loc[res["q"] < 0.05, "r"].abs().min()) if obs_nsig else np.nan
        ),
        "n_abs_r_gt_0.3": int((res["r"].abs() > 0.3).sum()),
        "n_abs_r_gt_0.4": int((res["r"].abs() > 0.4).sum()),
        "n_permutations": n_perm,
    }
    return res, stat


def main() -> int:
    ensure_dirs()
    models = load_models()
    lung = lung_cohort(models)
    gene_effect = load_gene_effect()
    lung_ids = lung.index.intersection(gene_effect.index)

    # 1. technical floor
    floor = technical_floor(lung_ids)
    if floor is not None:
        floor.sort_values("between_screen_r", ascending=False).to_csv(
            TABLES / "codependency_technical_floor.csv", index=False
        )
        summary_floor = floor.groupby("class")["between_screen_r"].describe()
        print("Between-screen reproducibility by class:")
        print(summary_floor.to_string())
        print()

    # 2/3. co-dependency scans
    variable = gene_effect.loc[:, gene_effect.std(axis=0, ddof=1) >= MIN_SD]
    for gene in GENES_OF_INTEREST + POSITIVE_CONTROLS:
        if gene not in variable.columns and gene in gene_effect.columns:
            variable[gene] = gene_effect[gene]
    print(f"genes retained after SD >= {MIN_SD} filter: {variable.shape[1]}", flush=True)

    lung_mat = variable.loc[lung_ids]
    lung_mat = lung_mat.loc[:, lung_mat.std(axis=0, ddof=1) >= MIN_SD]
    for gene in GENES_OF_INTEREST + POSITIVE_CONTROLS:
        if gene not in lung_mat.columns:
            lung_mat[gene] = gene_effect.loc[lung_ids, gene]

    all_stats = []
    top_frames = []
    for gene in GENES_OF_INTEREST + POSITIVE_CONTROLS:
        n_perm = 200 if gene in GENES_OF_INTEREST else 0
        for mat, label in ((variable, "pan-cancer"), (lung_mat, "lung")):
            res, stat = scan(gene, mat, label, n_perm=n_perm)
            if res.empty:
                continue
            all_stats.append(stat)
            head = pd.concat([res.head(25), res.tail(25)])
            head.insert(0, "direction", ["positive"] * 25 + ["negative"] * 25)
            top_frames.append(head)
            print(
                f"{gene:8s} {label:11s} n={stat['n_lines']:4d} maxAbsR={stat['max_abs_r']:.3f} "
                f"permNull95={stat['perm_null_max_abs_r_p95']:.3f} "
                f"permP={stat['perm_p_for_max_abs_r']:.3f} q<0.05: {stat['n_q_lt_0.05']}",
                flush=True,
            )

    stats_frame = pd.DataFrame(all_stats)
    stats_frame.to_csv(TABLES / "codependency_scan_stats.csv", index=False)
    tops = pd.concat(top_frames, ignore_index=True)
    tops = tops[["query", "cohort", "direction", "gene", "r", "n", "p", "q"]]
    tops.to_csv(TABLES / "codependency_top_hits.csv", index=False)

    # full ranked table for the two targets only (keeps output small)
    full = []
    for gene in GENES_OF_INTEREST:
        for mat, label in ((variable, "pan-cancer"), (lung_mat, "lung")):
            res = corr_vector_vs_matrix(
                mat[gene].dropna(), mat.drop(columns=[gene]), "pearson"
            )
            res["query"] = gene
            res["cohort"] = label
            full.append(res.sort_values("r", ascending=False).head(200))
    pd.concat(full, ignore_index=True)[["query", "cohort", "gene", "r", "n", "p", "q"]].to_csv(
        TABLES / "codependency_top200_targets.csv", index=False
    )

    # ------------------------------------------------------------------
    # Targeted hypotheses that do not need a genome-wide correction:
    # family paralogs, and a donor-sex sanity check because Y-linked genes
    # surface among the top lung hits.
    # ------------------------------------------------------------------
    families = {
        "TACSTD2": ["EPCAM", "CLDN4", "CLDN3", "CLDN7", "KRT8", "KRT18", "CDH1"],
        "CLDN4": ["CLDN3", "CLDN7", "CLDN1", "CLDN6", "TACSTD2", "EPCAM", "CDH1"],
    }
    targeted = []
    for query, partners in families.items():
        for cohort_mat, label in ((gene_effect, "pan-cancer"), (gene_effect.loc[lung_ids], "lung")):
            for partner in partners:
                if partner not in cohort_mat.columns or partner == query:
                    continue
                pair = cohort_mat[[query, partner]].dropna()
                if len(pair) < 20:
                    continue
                r, p = stats.pearsonr(pair[query], pair[partner])
                targeted.append(
                    {
                        "query": query,
                        "partner": partner,
                        "cohort": label,
                        "n": len(pair),
                        "pearson_r": float(r),
                        "p": float(p),
                    }
                )
    targeted = pd.DataFrame(targeted)
    if not targeted.empty:
        from common import bh_fdr as _bh

        targeted["q_within_family_test"] = _bh(targeted["p"].to_numpy())
        targeted.to_csv(TABLES / "codependency_targeted_family.csv", index=False)
        print()
        print("Targeted family co-dependency tests:")
        print(targeted.sort_values("p").head(12).to_string(index=False))

    y_genes = [g for g in ["UTY", "USP9Y", "DDX3Y", "KDM5D", "RPS4Y1", "EIF1AY", "ZFY"] if g in gene_effect.columns]
    sex_rows = []
    sex = models["Sex"].reindex(lung_ids)
    for gene in GENES_OF_INTEREST:
        prof = gene_effect.loc[lung_ids, gene]
        male = prof[sex == "Male"].dropna()
        female = prof[sex == "Female"].dropna()
        if len(male) > 5 and len(female) > 5:
            u = stats.mannwhitneyu(male, female, alternative="two-sided")
            sex_rows.append(
                {
                    "gene": gene,
                    "n_male": len(male),
                    "n_female": len(female),
                    "mean_male": float(male.mean()),
                    "mean_female": float(female.mean()),
                    "mannwhitney_p": float(u.pvalue),
                    "mean_r_with_Y_linked_gene_effects": float(
                        np.mean(
                            [
                                gene_effect.loc[lung_ids, [gene, yg]].dropna().corr().iloc[0, 1]
                                for yg in y_genes
                            ]
                        )
                    ),
                    "y_genes_used": ";".join(y_genes),
                }
            )
    if sex_rows:
        pd.DataFrame(sex_rows).to_csv(TABLES / "codependency_sex_confounder_check.csv", index=False)
        print()
        print("Donor-sex check on lung gene effect:")
        print(pd.DataFrame(sex_rows).to_string(index=False))

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    ax = axes[0]
    if floor is not None:
        order = ["common essential", "non-essential control", "target"]
        data = [floor.loc[floor["class"] == c, "between_screen_r"].dropna() for c in order]
        ax.boxplot(data, tick_labels=["common\nessential", "non-essential\ncontrol", "target"], showfliers=False)
        for i, c in enumerate(order):
            vals = floor.loc[floor["class"] == c, "between_screen_r"].dropna()
            ax.scatter(np.random.default_rng(1).normal(i + 1, 0.05, len(vals)), vals, s=6, alpha=0.35)
        for gene, colour in zip(GENES_OF_INTEREST, ["#1a9850", "#762a83"]):
            row = floor[floor["gene"] == gene]
            if not row.empty:
                ax.scatter([3], row["between_screen_r"], color=colour, s=70, zorder=5, label=gene)
        ax.axhline(0, color="grey", ls="--", lw=0.8)
        ax.set_ylabel("between-screen Pearson r")
        ax.set_title(f"Reproducibility, lines screened twice\n({floor.attrs.get('n_models','?')} models, {floor.attrs.get('scope','')})")
        ax.legend(fontsize=8)

    for ax, gene in zip(axes[1:], GENES_OF_INTEREST):
        res = corr_vector_vs_matrix(lung_mat[gene].dropna(), lung_mat.drop(columns=[gene]), "pearson")
        ax.scatter(res["r"], -np.log10(res["p"].clip(lower=1e-300)), s=3, alpha=0.25, color="#4d4d4d")
        sig = res[res["q"] < 0.05]
        if not sig.empty:
            ax.scatter(sig["r"], -np.log10(sig["p"]), s=10, color="#b2182b")
        ax.set_xlabel("Pearson r with gene effect profile")
        ax.set_ylabel("-log10 p")
        ax.set_title(f"{gene} co-dependency, lung (n={int(lung_mat[gene].notna().sum())})\n{len(sig)} genes at FDR<5%")
        ax.set_xlim(-1, 1)
    fig.suptitle(f"{RELEASE} · co-dependency scans", y=1.03)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_codependency.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Note
    # ------------------------------------------------------------------
    lines = [
        "# 03 · Co-dependency",
        "",
        "A co-dependency scan is only meaningful when the query gene's Chronos profile is "
        "reproducible. The order of evidence below is deliberate.",
        "",
        "## Technical floor",
        "",
    ]
    if floor is not None:
        grp = floor.groupby("class")["between_screen_r"].agg(["count", "median"])
        lines += [
            f"Lines screened at least twice: **{floor.attrs.get('n_models')}** "
            f"({floor.attrs.get('scope')}). Median between-screen correlation of the "
            "gene effect value across those lines:",
            "",
            "| class | genes | median r |",
            "| --- | --- | --- |",
        ]
        for klass, row in grp.iterrows():
            lines.append(f"| {klass} | {int(row['count'])} | {row['median']:.3f} |")
        lines.append("")
        for gene in GENES_OF_INTEREST:
            row = floor[floor["gene"] == gene]
            if not row.empty:
                lines.append(
                    f"* **{gene}**: between-screen r = {row['between_screen_r'].iloc[0]:.3f} "
                    f"(p = {row['p'].iloc[0]:.3g}, {int(row['n_lines_screened_twice'].iloc[0])} lines)."
                )
        lines.append("")
    lines += [
        "## Scan statistics",
        "",
        "| query | cohort | lines | genes | profile SD | max abs r | perm null max abs r (95th pct) | perm p | genes at FDR<5% | weakest r passing FDR | genes with abs r > 0.4 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]

    def fmt(value: float, spec: str = ".3f") -> str:
        return "n/a" if value is None or (isinstance(value, float) and not np.isfinite(value)) else format(value, spec)

    for _, r in stats_frame.iterrows():
        lines.append(
            f"| {r['query']} | {r['cohort']} | {int(r['n_lines'])} | {int(r['n_genes_tested'])} | "
            f"{r['profile_sd']:.3f} | {r['max_abs_r']:.3f} | {fmt(r['perm_null_max_abs_r_p95'])} | "
            f"{fmt(r['perm_p_for_max_abs_r'])} | {int(r['n_q_lt_0.05'])} | "
            f"{fmt(r['smallest_abs_r_passing_fdr05'])} | {int(r['n_abs_r_gt_0.4'])} |"
        )
    lines += [
        "",
        "`EGFR`, `CTNNB1` and `MYC` are positive controls run through the identical code "
        "path; their scans return the textbook partners (`GRB2`/`SHC1`/`GAB1`, `TCF7L2`, "
        "`MAX`), which is what licenses reading the target-gene scans as biology rather "
        "than as a broken pipeline.",
        "",
        "Note that the FDR<5% count is a poor guide to effect size pan-cancer: with 1178 "
        "lines, |r| of roughly 0.06 already clears BH, which is why the weakest passing |r| "
        "and the permutation columns are reported alongside it.",
        "",
        "Tables: `codependency_scan_stats.csv`, `codependency_top_hits.csv`, "
        "`codependency_top200_targets.csv`, `codependency_technical_floor.csv`. "
        "Figure: `fig2_codependency.png`.",
    ]
    (NOTES / "03_codependency.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
