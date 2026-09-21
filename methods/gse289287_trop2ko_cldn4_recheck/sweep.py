#!/usr/bin/env python3
"""Method / threshold / gene-set sweep for Trop-2 KO collateral pathways.

The family is locked in gene_sets_sweep.json (IFN, APM, NHEJ, STING /
cytosolic DNA). CLDN4 is re-tested under each normalization and stays a
recorded row; it is not the selection criterion.

Primary rank is the author DESeq2 Wald statistic. Other ranks,
normalizations, thresholds, and leave-one-out scores are the sweep.
GSEA is the same prerank engine as gsea_core.py (weighted KS, 1000
gene-set permutations, seed 42). Positive NES = enriched in Trop-2 KO.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import analyze
from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
FAMILY = json.loads((ROOT / "gene_sets_sweep.json").read_text())["sets"]

# HGNC previous symbols that are the names actually used in this table.
# Applied only when the current symbol is absent and the previous symbol is present.
ALIAS = {"STING1": "TMEM173", "CGAS": "MB21D1", "H2AX": "H2AFX"}

SHORT = {
    "Interferon Alpha Response": "IFN-α Hallmark",
    "Interferon Gamma Response": "IFN-γ Hallmark",
    "Interferon Alpha/Beta Signaling R-HSA-909733": "IFN-α/β Reactome",
    "Interferon Gamma Signaling R-HSA-877300": "IFN-γ Reactome",
    "Response To Type I Interferon (GO:0034340)": "IFN-I response GO",
    "Type I Interferon-Mediated Signaling Pathway (GO:0060337)": "IFN-I signaling GO",
    "Antigen processing and presentation": "APM KEGG (I+II)",
    "Antigen Processing And Presentation Of Peptide Antigen Via MHC Class I (GO:0002474)": "APM MHC-I GO",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION": "APM MHC-I custom",
    "STING Mediated Induction Of Host Immune Responses R-HSA-1834941": "STING Reactome",
    "Cytosolic Sensors Of Pathogen-Associated DNA R-HSA-1834949": "DNA sensors Reactome",
    "Regulation Of Innate Immune Responses To Cytosolic DNA R-HSA-3134975": "Cytosolic DNA reg.",
    "IRF3-mediated Induction Of Type I IFN R-HSA-3270619": "IRF3 type-I IFN",
    "Cytosolic DNA-sensing pathway": "DNA sensing KEGG",
    "Non-homologous end-joining": "NHEJ KEGG",
    "Double-Strand Break Repair Via Nonhomologous End Joining (GO:0006303)": "NHEJ GO",
    "Double-Strand Break Repair Via Classical Nonhomologous End Joining (GO:0097680)": "NHEJ classical GO",
    "Nonhomologous End-Joining (NHEJ) R-HSA-5693571": "NHEJ Reactome",
    "HDR Thru MMEJ (alt-NHEJ) R-HSA-5685939": "alt-NHEJ MMEJ",
}
AXIS_COLOR = {"IFN": "#E45756", "APM": "#54A24B", "STING": "#B279A2", "NHEJ": "#4C78A8"}


def symbol_table(df: pd.DataFrame) -> pd.DataFrame:
    pc = df[df["biotype"] == "protein_coding"].copy()
    pc = pc.sort_values("baseMean", ascending=False).drop_duplicates("Feature_name")
    other = df[~df["Feature_name"].isin(pc["Feature_name"])]
    other = other.sort_values("baseMean", ascending=False).drop_duplicates("Feature_name")
    return pd.concat([pc, other], ignore_index=True)


def resolve_genes(genes: list[str], have: set[str]) -> tuple[list[str], list[str]]:
    out, aliased = [], []
    seen = set()
    for g in genes:
        use = g if g in have else ALIAS.get(g)
        if use is None or use not in have or use in seen:
            continue
        if use != g:
            aliased.append(f"{g}->{use}")
        seen.add(use)
        out.append(use)
    return out, aliased


def matrix_from(pc: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    m = pc.set_index("Feature_name")[cols].astype(float)
    m = m.groupby(level=0).mean()
    return m


def median_of_ratios(raw: pd.DataFrame) -> pd.Series:
    """DESeq2 size factors. Genes with a zero in any sample are excluded from the median."""
    x = raw.to_numpy(dtype=float)
    ok = np.all(x > 0, axis=1)
    use = x[ok]
    log_gm = np.log(use).mean(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.log(use) - log_gm[:, None]
    sf = np.exp(np.median(ratio, axis=0))
    sf = sf / np.exp(np.mean(np.log(sf)))
    return pd.Series(sf, index=raw.columns)


def welch_t_matrix(ko: np.ndarray, wt: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n1, n2 = ko.shape[1], wt.shape[1]
    m1, m2 = ko.mean(axis=1), wt.mean(axis=1)
    v1, v2 = ko.var(axis=1, ddof=1), wt.var(axis=1, ddof=1)
    se2 = v1 / n1 + v2 / n2
    t = np.full(ko.shape[0], np.nan)
    p = np.full(ko.shape[0], np.nan)
    ok = se2 > 0
    t[ok] = (m1[ok] - m2[ok]) / np.sqrt(se2[ok])
    df = np.full(ko.shape[0], np.nan)
    df[ok] = (se2[ok] ** 2) / ((v1[ok] / n1) ** 2 / (n1 - 1) + (v2[ok] / n2) ** 2 / (n2 - 1))
    p[ok] = 2 * stats.t.sf(np.abs(t[ok]), df[ok])
    return t, p


def exact_greater(wt: np.ndarray, ko: np.ndarray) -> tuple[float, float, float, int]:
    """Effect KO-WT, two-sided exact p, one-sided (KO>WT) exact p."""
    values = np.concatenate([np.asarray(wt, float), np.asarray(ko, float)])
    n_wt = len(wt)
    obs = float(np.asarray(ko, float).mean() - np.asarray(wt, float).mean())
    n = 0
    ge = 0
    two = 0
    # Enumerate WT index sets. C(7,3)=35.
    from itertools import combinations

    for idx in combinations(range(len(values)), n_wt):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(idx)] = True
        diff = float(values[~mask].mean() - values[mask].mean())
        n += 1
        if diff + 1e-12 >= obs:
            ge += 1
        if abs(diff) + 1e-12 >= abs(obs):
            two += 1
    return obs, two / n, ge / n, n


def score_samples(mat: pd.DataFrame, genes: list[str], wt_cols: list[str], ko_cols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    present = [g for g in genes if g in mat.index]
    sub = np.log2(mat.loc[present].to_numpy(dtype=float) + 1.0)
    # columns follow mat column order; caller passes cols in wt then we index by name
    wt = sub[:, [mat.columns.get_loc(c) for c in wt_cols]].mean(axis=0)
    ko = sub[:, [mat.columns.get_loc(c) for c in ko_cols]].mean(axis=0)
    return wt, ko


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    df = analyze.load_table()
    pc = symbol_table(df)
    have = set(pc["Feature_name"])
    wt_raw_cols, ko_raw_cols = analyze.group_cols(df, "rawCounts")
    wt_norm_cols, ko_norm_cols = analyze.group_cols(df, "normCounts")
    raw_cols = wt_raw_cols + ko_raw_cols
    norm_cols = wt_norm_cols + ko_norm_cols

    raw = matrix_from(pc, raw_cols)
    author = matrix_from(pc, norm_cols)
    lib = raw.sum(axis=0)
    cpm = raw.div(lib, axis=1) * 1e6
    # Size factors from protein-coding genes only.
    pc_symbols = set(pc.loc[pc["biotype"] == "protein_coding", "Feature_name"])
    sf = median_of_ratios(raw.loc[raw.index.isin(pc_symbols)])
    mor = raw.div(sf, axis=1)

    mats = {"author_norm": author, "cpm": cpm, "median_ratio": mor}
    wt_for = {
        "author_norm": wt_norm_cols,
        "cpm": wt_raw_cols,
        "median_ratio": wt_raw_cols,
    }
    ko_for = {
        "author_norm": ko_norm_cols,
        "cpm": ko_raw_cols,
        "median_ratio": ko_raw_cols,
    }

    resolved = {}
    alias_rows = []
    for key, spec in FAMILY.items():
        genes, aliased = resolve_genes(spec["genes"], have)
        resolved[key] = genes
        alias_rows.append(
            {
                "set": spec["name"],
                "axis": spec["axis"],
                "n_listed": len(spec["genes"]),
                "n_in_table": len(genes),
                "n_missing": len(spec["genes"]) - len(genes) - (0 if not aliased else 0),
                "aliased": ";".join(aliased),
            }
        )
    # n_missing recount properly
    for row, (key, spec) in zip(alias_rows, FAMILY.items()):
        row["n_missing"] = len(spec["genes"]) - len(resolved[key])
        # aliased genes are in the table, so they are not missing
        # resolve_genes drops symbols that alias onto an already kept gene; count those as mapped
    pd.DataFrame(alias_rows).to_csv(TABLES / "sweep_set_coverage.tsv", sep="\t", index=False)

    # Ranks. Higher = up in KO.
    stat = pc.drop_duplicates("Feature_name").set_index("Feature_name")["stat"].astype(float)
    lfc = pc.drop_duplicates("Feature_name").set_index("Feature_name")["log2FoldChange"].astype(float)
    ranks = {
        "author_wald": stat.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False),
        "author_log2FC": lfc.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False),
    }
    for norm, mat in mats.items():
        wt_cols, ko_cols = wt_for[norm], ko_for[norm]
        logx = np.log2(mat.to_numpy(dtype=float) + 1.0)
        wt_i = [mat.columns.get_loc(c) for c in wt_cols]
        ko_i = [mat.columns.get_loc(c) for c in ko_cols]
        t, _p = welch_t_matrix(logx[:, ko_i], logx[:, wt_i])
        s = pd.Series(t, index=mat.index).replace([np.inf, -np.inf], np.nan).dropna()
        ranks[f"welch_{norm}"] = s.sort_values(ascending=False)

    gene_sets = {FAMILY[k]["name"]: resolved[k] for k in FAMILY}
    gsea_frames = []
    for rank_name, rank in ranks.items():
        print(f"GSEA {rank_name} n={len(rank)}", flush=True)
        res = gsea_prerank(rank, gene_sets, nperm=NPERM, seed=SEED, min_size=5, max_size=500)
        res["rank"] = rank_name
        res["fdr_within_rank"] = bh_fdr(res["nom_p"])
        gsea_frames.append(res)
    gsea = pd.concat(gsea_frames, ignore_index=True)
    gsea["fdr_sweep"] = bh_fdr(gsea["nom_p"])
    meta = {spec["name"]: spec for spec in FAMILY.values()}
    gsea["axis"] = gsea["term"].map(lambda t: meta[t]["axis"])
    gsea["short"] = gsea["term"].map(SHORT)
    gsea["below_usual_min8"] = gsea["n_set_in_rank"] < 8
    gsea.to_csv(TABLES / "sweep_gsea.tsv", sep="\t", index=False)

    # Sample scores.
    score_rows = []
    for key, spec in FAMILY.items():
        genes = resolved[key]
        if len(genes) < 3:
            continue
        for norm, mat in mats.items():
            wt, ko = score_samples(mat, genes, wt_for[norm], ko_for[norm])
            obs, p_two, p_up, nperm = exact_greater(wt, ko)
            score_rows.append(
                {
                    "set": spec["name"],
                    "short": SHORT[spec["name"]],
                    "axis": spec["axis"],
                    "normalization": norm,
                    "n_genes": len([g for g in genes if g in mat.index]),
                    "score_KO_minus_WT": obs,
                    "exact_p_two": p_two,
                    "exact_p_KO_gt_WT": p_up,
                    "n_permutations": nperm,
                }
            )
    scores = pd.DataFrame(score_rows)
    scores["fdr_two_sweep"] = bh_fdr(scores["exact_p_two"])
    scores["fdr_up_sweep"] = bh_fdr(scores["exact_p_KO_gt_WT"])
    scores.to_csv(TABLES / "sweep_sample_scores.tsv", sep="\t", index=False)

    # Threshold sweep on the author table.
    ann = pc.drop_duplicates("Feature_name").set_index("Feature_name")
    thr_rows = []
    cuts = [("padj", 0.05), ("padj", 0.10), ("padj", 0.20), ("pvalue", 0.05), ("pvalue", 0.10)]
    for key, spec in FAMILY.items():
        genes = [g for g in resolved[key] if g in ann.index]
        sub = ann.loc[genes]
        for col, cut in cuts:
            tested = sub[np.isfinite(sub[col])]
            up = tested[(tested[col] < cut) & (tested["log2FoldChange"] > 0)]
            down = tested[(tested[col] < cut) & (tested["log2FoldChange"] < 0)]
            thr_rows.append(
                {
                    "set": spec["name"],
                    "short": SHORT[spec["name"]],
                    "axis": spec["axis"],
                    "rule": f"{col}<{cut}",
                    "n_genes_in_table": len(genes),
                    "n_with_stat": int(np.isfinite(sub[col]).sum()) if col in sub else 0,
                    "n_up": int(len(up)),
                    "n_down": int(len(down)),
                    "up_genes": ",".join(up.sort_values("log2FoldChange", ascending=False).index[:12]),
                }
            )
        # extra: |log2FC|>0.5 and nominal p<0.05
        hit = sub[(sub["pvalue"] < 0.05) & (sub["log2FoldChange"].abs() > 0.5)]
        thr_rows.append(
            {
                "set": spec["name"],
                "short": SHORT[spec["name"]],
                "axis": spec["axis"],
                "rule": "pvalue<0.05 and |log2FC|>0.5",
                "n_genes_in_table": len(genes),
                "n_with_stat": int(np.isfinite(sub["pvalue"]).sum()),
                "n_up": int((hit["log2FoldChange"] > 0).sum()),
                "n_down": int((hit["log2FoldChange"] < 0).sum()),
                "up_genes": ",".join(hit[hit["log2FoldChange"] > 0].sort_values("log2FoldChange", ascending=False).index[:12]),
            }
        )
    pd.DataFrame(thr_rows).to_csv(TABLES / "sweep_thresholds.tsv", sep="\t", index=False)

    # CLDN4 under each normalization. Not used to pick the pathway panel.
    cldn_rows = []
    row = analyze.one_symbol(df, "CLDN4")
    cldn_rows.append(
        {
            "normalization": "author_DESeq2",
            "log2FC_or_mean_diff": float(row["log2FoldChange"]),
            "pvalue": float(row["pvalue"]),
            "padj_or_note": float(row["padj"]),
            "note": "author Wald p and genome-wide padj",
        }
    )
    for norm, mat in mats.items():
        if "CLDN4" not in mat.index:
            continue
        vals = np.log2(mat.loc["CLDN4"].to_numpy(dtype=float) + 1.0)
        wt = vals[[mat.columns.get_loc(c) for c in wt_for[norm]]]
        ko = vals[[mat.columns.get_loc(c) for c in ko_for[norm]]]
        tt = stats.ttest_ind(ko, wt, equal_var=False)
        obs, p_two, p_up, nperm = exact_greater(wt, ko)
        cldn_rows.append(
            {
                "normalization": norm,
                "log2FC_or_mean_diff": obs,
                "pvalue": float(tt.pvalue),
                "padj_or_note": p_two,
                "note": f"welch p in pvalue column; exact two-sided perm p in padj_or_note; nperm={nperm}; one-sided KO>WT p={p_up:.4f}",
            }
        )
    pd.DataFrame(cldn_rows).to_csv(TABLES / "sweep_cldn4_norms.tsv", sep="\t", index=False)

    # Leave-one-out ISG / STING / NHEJ scores on author_norm, plus GSEA LOO on the Wald rank.
    loo_targets = [
        "Interferon Alpha Response",
        "Interferon Gamma Response",
        "STING Mediated Induction Of Host Immune Responses R-HSA-1834941",
        "Non-homologous end-joining",
    ]
    loo_rows = []
    wald = ranks["author_wald"]
    for name in loo_targets:
        genes = gene_sets[name]
        wt, ko = score_samples(mats["author_norm"], genes, wt_norm_cols, ko_norm_cols)
        full_obs, full_two, full_up, nperm = exact_greater(wt, ko)
        loo_rows.append(
            {
                "set": name,
                "dropped": "(none)",
                "kind": "full_author_norm_score",
                "n_genes": len(genes),
                "effect": full_obs,
                "exact_p_two": full_two,
                "exact_p_KO_gt_WT": full_up,
                "nes": np.nan,
                "nom_p": np.nan,
            }
        )
        for g in genes:
            keep = [x for x in genes if x != g]
            wt_i, ko_i = score_samples(mats["author_norm"], keep, wt_norm_cols, ko_norm_cols)
            obs, p_two, p_up, _ = exact_greater(wt_i, ko_i)
            loo_rows.append(
                {
                    "set": name,
                    "dropped": g,
                    "kind": "loo_gene_author_norm_score",
                    "n_genes": len(keep),
                    "effect": obs,
                    "exact_p_two": p_two,
                    "exact_p_KO_gt_WT": p_up,
                    "nes": np.nan,
                    "nom_p": np.nan,
                }
            )
        # Leave-one-animal-out of the full score. Descriptive effect only.
        logm = np.log2(mats["author_norm"].loc[[g for g in genes if g in mats["author_norm"].index]] + 1.0)
        for col, group in [(c, "WT") for c in wt_norm_cols] + [(c, "KO") for c in ko_norm_cols]:
            wt_cols = [c for c in wt_norm_cols if c != col]
            ko_cols = [c for c in ko_norm_cols if c != col]
            if len(wt_cols) < 2 or len(ko_cols) < 2:
                continue
            obs = float(logm[ko_cols].to_numpy().mean() - logm[wt_cols].to_numpy().mean())
            loo_rows.append(
                {
                    "set": name,
                    "dropped": f"sample_{analyze.animal_of(col)}_{group}",
                    "kind": "loo_sample_author_norm_score",
                    "n_genes": logm.shape[0],
                    "effect": obs,
                    "exact_p_two": np.nan,
                    "exact_p_KO_gt_WT": np.nan,
                    "nes": np.nan,
                    "nom_p": np.nan,
                }
            )
        print(f"GSEA leave-one-out {name}", flush=True)
        base = gsea[(gsea["rank"] == "author_wald") & (gsea["term"] == name)].iloc[0]
        loo_rows.append(
            {
                "set": name,
                "dropped": "(none)",
                "kind": "full_wald_gsea",
                "n_genes": int(base["n_set_in_rank"]),
                "effect": float(base["mean_stat"]),
                "exact_p_two": np.nan,
                "exact_p_KO_gt_WT": np.nan,
                "nes": float(base["nes"]),
                "nom_p": float(base["nom_p"]),
            }
        )
        for g in genes:
            keep = [x for x in genes if x != g]
            res = gsea_prerank(wald, {name: keep}, nperm=NPERM, seed=SEED, min_size=5, max_size=500)
            if res.empty:
                continue
            loo_rows.append(
                {
                    "set": name,
                    "dropped": g,
                    "kind": "loo_gene_wald_gsea",
                    "n_genes": int(res.iloc[0]["n_set_in_rank"]),
                    "effect": float(res.iloc[0]["mean_stat"]),
                    "exact_p_two": np.nan,
                    "exact_p_KO_gt_WT": np.nan,
                    "nes": float(res.iloc[0]["nes"]),
                    "nom_p": float(res.iloc[0]["nom_p"]),
                }
            )
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(TABLES / "sweep_loo.tsv", sep="\t", index=False)

    # Gene-level table for the two tight sets used on the panel.
    gene_rows = []
    for name in [
        "STING Mediated Induction Of Host Immune Responses R-HSA-1834941",
        "Non-homologous end-joining",
        "Interferon Alpha Response",
    ]:
        for g in gene_sets[name]:
            if g not in ann.index:
                gene_rows.append({"set": name, "symbol": g, "present": False})
                continue
            r = ann.loc[g]
            gene_rows.append(
                {
                    "set": name,
                    "symbol": g,
                    "present": True,
                    "log2FC": float(r["log2FoldChange"]),
                    "lfcSE": float(r["lfcSE"]),
                    "pvalue": float(r["pvalue"]),
                    "padj": float(r["padj"]) if np.isfinite(r["padj"]) else np.nan,
                    "stat": float(r["stat"]) if np.isfinite(r["stat"]) else np.nan,
                    "baseMean": float(r["baseMean"]),
                }
            )
    pd.DataFrame(gene_rows).to_csv(TABLES / "sweep_focus_genes.tsv", sep="\t", index=False)

    plot(gsea, scores, loo)
    primary = gsea[gsea["rank"] == "author_wald"].sort_values("nom_p")
    print(primary[["short", "axis", "nes", "nom_p", "fdr_within_rank", "n_set_in_rank"]].to_string(index=False))
    print("--- sample author_norm ---")
    print(
        scores[scores["normalization"] == "author_norm"][
            ["short", "score_KO_minus_WT", "exact_p_two", "exact_p_KO_gt_WT", "fdr_up_sweep"]
        ].to_string(index=False)
    )


def plot(gsea: pd.DataFrame, scores: pd.DataFrame, loo: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    primary = gsea[gsea["rank"] == "author_wald"].copy()
    order = [s for s in SHORT.values() if s in set(primary["short"])]
    primary["short"] = pd.Categorical(primary["short"], order, ordered=True)
    primary = primary.sort_values("short")

    fig = plt.figure(figsize=(11.6, 9.2), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.5, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(primary))
    colors = [AXIS_COLOR[a] for a in primary["axis"]]
    ax.axvline(0, color="#888888", lw=0.8)
    ax.barh(y, primary["nes"], color=colors, height=0.72)
    ax.set_xlim(float(primary["nes"].min()) - 0.35, float(primary["nes"].max()) + 0.7)
    for i, r in enumerate(primary.itertuples()):
        if r.fdr_within_rank < 0.05:
            ax.plot(r.nes + 0.08, i, marker="*", color="#222222", markersize=5)
    ax.set_yticks(y)
    ax.set_yticklabels(primary["short"])
    ax.set_xlabel("NES, author Wald rank (positive = up in Trop-2 KO)")
    ax.set_title("A   Locked family. Three IFN sets pass BH FDR")
    ax.invert_yaxis()

    ax = fig.add_subplot(gs[0, 1])
    rank_order = ["author_wald", "author_log2FC", "welch_author_norm", "welch_cpm", "welch_median_ratio"]
    rank_lab = ["Wald", "log2FC", "Welch\nauthor", "Welch\nCPM", "Welch\nsize factor"]
    heat = gsea.pivot(index="short", columns="rank", values="nes").reindex(index=order, columns=rank_order)
    vmax = float(np.nanmax(np.abs(heat.to_numpy())))
    im = ax.imshow(heat.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(rank_order)))
    ax.set_xticklabels(rank_lab, fontsize=7)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=6.5)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.to_numpy()[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=5.5, color="#111111")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="NES")
    ax.set_title("B   Sweep of five ranks. Sign of the IFN sets holds")

    ax = fig.add_subplot(gs[1, 0])
    sub = loo[(loo["set"] == "Interferon Alpha Response") & (loo["kind"] == "loo_gene_wald_gsea")].copy()
    full = loo[(loo["set"] == "Interferon Alpha Response") & (loo["kind"] == "full_wald_gsea")].iloc[0]
    ax.axvline(full["nes"], color="#E45756", lw=1.1, label=f"full set NES {full['nes']:.2f}")
    ax.hist(sub["nes"], bins=12, color="#F4C2C2", edgecolor="#E45756")
    ax.set_xlabel("NES after dropping one Hallmark IFN-α gene")
    ax.set_ylabel("Genes dropped")
    ax.set_title("C   Leave-one-out. NES stays above 2.17")
    ax.legend(frameon=False, fontsize=7)

    ax = fig.add_subplot(gs[1, 1])
    focus = pd.read_csv(TABLES / "sweep_focus_genes.tsv", sep="\t")
    want = ["TMEM173", "MB21D1", "TBK1", "IRF3", "IFI16", "PRKDC", "POLM", "XRCC4", "XRCC5", "XRCC6", "LIG4", "ISG15", "IFI44L"]
    bits = []
    for sym in want:
        hit = focus[(focus["symbol"] == sym) & (focus["present"] == True)]
        if hit.empty:
            continue
        bits.append(hit.iloc[0])
    block = pd.DataFrame(bits).iloc[::-1]
    yy = np.arange(len(block))
    sig = block["padj"].to_numpy(dtype=float) < 0.05
    ax.axvline(0, color="#888888", lw=0.8)
    ax.errorbar(block["log2FC"], yy, xerr=1.96 * block["lfcSE"], fmt="none", ecolor="#C8C8C8", elinewidth=0.7)
    ax.scatter(block["log2FC"], yy, c=np.where(sig, "#E45756", "#555555"), s=22, zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(block["symbol"])
    ax.set_xlabel("Author log2FC ± 1.96 × lfcSE")
    ax.set_title("D   STING1/cGAS stay flat. Red = padj < 0.05")

    fig.savefig(FIGURES / "fig_sweep_collateral.png", dpi=160)
    fig.savefig(FIGURES / "fig_sweep_collateral.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
