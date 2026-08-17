#!/usr/bin/env python3
"""CellChat-style LR on public GSE131907: malignant CLDN4-high vs low → T/NK.

Author cell labels are used (no lineage reconstruction, no CopyKAT).
CellChat R and LIANA are not run. Probability follows Jin et al. 2021
(10% truncated mean, Hill K_h=0.5, CellChatDB v2 protein pairs).
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import trim_mean

ROOT = Path(__file__).resolve().parents[1]
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
NBOOT = 100
SEED = 1
MIN_GROUP = 25

TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "mBrain")
MALIG_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD3E", "CD8A", "NKG7"]

BARRIER_LIGANDS = {
    "CDH1", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "F11R", "JAM2", "JAM3",
    "CEACAM1", "CEACAM5", "CEACAM6", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4",
    "PVR", "EPCAM", "DSG2", "DSC2", "CADM1",
}
INHIB_LIGANDS = {
    "CD274", "PDCD1LG2", "LGALS9", "HLA-E", "HLA-G", "HLA-F", "TGFB1", "TGFB2",
    "TGFB3", "CD80", "CD86", "CD276", "VSIR", "PVR", "NECTIN2", "CD47", "CDH1",
}
RECRUIT_LIGANDS = {
    "CXCL9", "CXCL10", "CXCL11", "CXCL16", "CCL5", "CCL3", "CCL4", "IL15",
    "IL2", "IL18", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3",
}
ATTACK_LIGANDS = {"IFNG", "TNF", "FASLG", "TNFSF10", "LTA"}


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_lr(db_dir: Path, matrix_genes: set[str]) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    inter = inter.rename(columns={"ligand.symbol": "ligand_symbol", "receptor.symbol": "receptor_symbol"})
    rows = []
    for rec in inter.itertuples(index=False):
        lig = parse_symbols(getattr(rec, "ligand_symbol", rec.ligand))
        recp = parse_symbols(getattr(rec, "receptor_symbol", rec.receptor))
        if not lig:
            lig = parse_symbols(rec.ligand)
        if not recp:
            recp = parse_symbols(rec.receptor)
        if not lig or not recp:
            continue
        if any(g not in matrix_genes for g in lig + recp):
            continue
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": tuple(lig),
                "receptor_genes": tuple(recp),
            }
        )
    return pd.DataFrame(rows)


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


def matrix_gene_list(matrix_path: Path) -> list[str]:
    genes = []
    with gzip.open(matrix_path, "rt") as handle:
        handle.readline()
        for line in handle:
            gene = line.split("\t", 1)[0].split(".")[0]
            genes.append(gene)
    return genes


def stream_matrix(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  stream genes={n_streamed} stored={len(found)}", flush=True)
    print(f"stream done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, n_streamed


def geom_mean_rows(mat: np.ndarray) -> np.ndarray:
    if mat.ndim == 1:
        return mat
    if mat.shape[0] == 1:
        return mat[0]
    out = np.exp(np.mean(np.log(np.clip(mat, 1e-12, None)), axis=0))
    out[np.any(mat <= 0, axis=0)] = 0.0
    return out


def hill_prob(lig: float, rec: float, kh: float = KH) -> float:
    if lig <= 0 or rec <= 0:
        return 0.0
    prod = lig * rec
    return float(prod / (kh + prod))


def group_trim_means(expr: np.ndarray, pos: np.ndarray, labels: np.ndarray, groups: list[str]):
    n_g = expr.shape[0]
    means = np.zeros((n_g, len(groups)), dtype=np.float64)
    props = np.zeros((n_g, len(groups)), dtype=np.float64)
    counts = {}
    for j, g in enumerate(groups):
        idx = np.flatnonzero(labels == g)
        counts[g] = int(idx.size)
        if idx.size == 0:
            continue
        sub = expr[:, idx]
        if idx.size == 1:
            means[:, j] = sub[:, 0]
        else:
            means[:, j] = trim_mean(sub, TRIM, axis=1)
        props[:, j] = pos[:, idx].mean(axis=1)
    return means, props, counts


def complex_value(gene_means, gene_props, gene_index, subunits):
    ix = [gene_index[g] for g in subunits]
    mu = geom_mean_rows(gene_means[ix])
    pr = gene_props[ix].min(axis=0) if len(ix) > 1 else gene_props[ix[0]]
    return mu, pr


def pair_table(lr, gene_means, gene_props, gene_index, groups, counts, src_tgt):
    gpos = {g: i for i, g in enumerate(groups)}
    cache = {}

    def cached(subunits):
        hit = cache.get(subunits)
        if hit is None:
            hit = complex_value(gene_means, gene_props, gene_index, subunits)
            cache[subunits] = hit
        return hit

    rows = []
    for rec in lr.itertuples(index=False):
        lig_mu, lig_pr = cached(rec.ligand_genes)
        rec_mu, rec_pr = cached(rec.receptor_genes)
        for src, tgt in src_tgt:
            if src not in gpos or tgt not in gpos:
                continue
            i, j = gpos[src], gpos[tgt]
            if counts.get(src, 0) < MIN_GROUP or counts.get(tgt, 0) < MIN_GROUP:
                continue
            detected = (lig_pr[i] >= EXPR_PROP) and (rec_pr[j] >= EXPR_PROP)
            prob = hill_prob(float(lig_mu[i]), float(rec_mu[j])) if detected else 0.0
            rows.append(
                {
                    "interaction_name": rec.interaction_name,
                    "pathway_name": rec.pathway_name,
                    "annotation": rec.annotation,
                    "ligand": rec.ligand,
                    "receptor": rec.receptor,
                    "ligand_genes": "|".join(rec.ligand_genes),
                    "receptor_genes": "|".join(rec.receptor_genes),
                    "source": src,
                    "target": tgt,
                    "n_source": counts[src],
                    "n_target": counts[tgt],
                    "ligand_mean": float(lig_mu[i]),
                    "receptor_mean": float(rec_mu[j]),
                    "ligand_prop": float(lig_pr[i]),
                    "receptor_prop": float(rec_pr[j]),
                    "detected": bool(detected),
                    "prob": prob,
                }
            )
    return pd.DataFrame(rows)


def permute_mal_high_low(expr, pos, gene_index, lr, labels, groups, src_tgt, mal_mask, nboot, rng):
    obs_means, obs_props, counts = group_trim_means(expr, pos, labels, groups)
    obs = pair_table(lr, obs_means, obs_props, gene_index, groups, counts, src_tgt)
    obs_prob = obs["prob"].to_numpy()
    ge = np.zeros(len(obs), dtype=np.int32)
    mal_idx = np.flatnonzero(mal_mask)
    mal_labs = labels[mal_idx].copy()
    work = labels.copy()
    for b in range(nboot):
        work[mal_idx] = rng.permutation(mal_labs)
        means, props, counts_b = group_trim_means(expr, pos, work, groups)
        perm = pair_table(lr, means, props, gene_index, groups, counts_b, src_tgt)
        ge += (perm["prob"].to_numpy() >= obs_prob).astype(np.int32)
        if (b + 1) % 25 == 0:
            print(f"    perm {b + 1}/{nboot}", flush=True)
    pval = (ge + 1) / (nboot + 1)
    obs = obs.copy()
    obs["pval"] = pval
    obs["significant"] = (obs["pval"] < 0.05) & (obs["prob"] > 0) & obs["detected"]
    return obs


def ligand_class(genes: str) -> str:
    parts = set(genes.split("|"))
    tags = []
    if parts & BARRIER_LIGANDS:
        tags.append("barrier")
    if parts & INHIB_LIGANDS:
        tags.append("inhibitory")
    if parts & RECRUIT_LIGANDS:
        tags.append("recruit")
    if parts & ATTACK_LIGANDS:
        tags.append("attack")
    return "|".join(tags) if tags else "other"


def contrast_high_low(df: pd.DataFrame, high: str, low: str, partner: str, direction: str) -> pd.DataFrame:
    if direction == "outgoing":
        a = df[(df.source == high) & (df.target == partner)].set_index("interaction_name")
        b = df[(df.source == low) & (df.target == partner)].set_index("interaction_name")
    else:
        a = df[(df.source == partner) & (df.target == high)].set_index("interaction_name")
        b = df[(df.source == partner) & (df.target == low)].set_index("interaction_name")
    common = a.index.intersection(b.index)
    if len(common) == 0:
        return pd.DataFrame()
    out = a.loc[common, ["pathway_name", "annotation", "ligand", "receptor", "ligand_genes", "receptor_genes"]].copy()
    out["direction"] = direction
    out["partner"] = partner
    out["prob_high"] = a.loc[common, "prob"].to_numpy()
    out["prob_low"] = b.loc[common, "prob"].to_numpy()
    out["pval_high"] = a.loc[common, "pval"].to_numpy()
    out["pval_low"] = b.loc[common, "pval"].to_numpy()
    out["sig_high"] = a.loc[common, "significant"].to_numpy()
    out["sig_low"] = b.loc[common, "significant"].to_numpy()
    out["delta_prob"] = out["prob_high"] - out["prob_low"]
    out["ligand_mean_high"] = a.loc[common, "ligand_mean"].to_numpy()
    out["ligand_mean_low"] = b.loc[common, "ligand_mean"].to_numpy()
    out["receptor_mean_partner"] = a.loc[common, "receptor_mean"].to_numpy()
    out["ligand_prop_high"] = a.loc[common, "ligand_prop"].to_numpy()
    out["ligand_prop_low"] = b.loc[common, "ligand_prop"].to_numpy()
    out["receptor_prop_partner"] = a.loc[common, "receptor_prop"].to_numpy()
    out["n_high"] = a.loc[common, "n_source" if direction == "outgoing" else "n_target"].to_numpy()
    out["n_low"] = b.loc[common, "n_source" if direction == "outgoing" else "n_target"].to_numpy()
    out["n_partner"] = a.loc[common, "n_target" if direction == "outgoing" else "n_source"].to_numpy()
    out["ligand_class"] = [ligand_class(x) for x in out["ligand_genes"]]
    out["sig_either"] = out["sig_high"] | out["sig_low"]
    out["sig_diff"] = out["sig_either"] & (out["delta_prob"] != 0)
    return out.reset_index()


def is_malignant(row: pd.Series) -> bool:
    return (row.Cell_type == "Epithelial cells") and (row.Cell_subtype in MALIG_SUBTYPES)


def is_tnk(row: pd.Series) -> bool:
    return row.Cell_type in {"T lymphocytes", "NK cells"}


def make_split(cldn4, origin, mal, tnk, mode: str):
    n = len(cldn4)
    labels = np.array(["drop"] * n, dtype=object)
    if mode == "tumor_tertile":
        keep = np.isin(origin, TUMOR_ORIGINS)
        use_mal = mal & keep
        use_tnk = tnk & keep
        split = "tertile"
        info_samples = "tLung+tL/B+mLN+mBrain"
    elif mode == "tlung_tertile":
        keep = origin == "tLung"
        use_mal = mal & keep
        use_tnk = tnk & keep
        split = "tertile"
        info_samples = "tLung"
    elif mode == "mets_tertile":
        keep = np.isin(origin, ("tL/B", "mLN", "mBrain"))
        use_mal = mal & keep
        use_tnk = tnk & keep
        split = "tertile"
        info_samples = "tL/B+mLN+mBrain (author Malignant cells / tS*)"
    else:
        raise ValueError(mode)

    vals = cldn4[use_mal]
    if vals.size < MIN_GROUP * 2:
        return labels, np.zeros(n, dtype=bool), {"skip": "too_few_malignant"}
    q1, q2 = np.quantile(vals, [1 / 3, 2 / 3])
    labels[use_tnk] = "TNK"
    labels[use_mal & (cldn4 <= q1)] = "Mal_low"
    labels[use_mal & (cldn4 >= q2)] = "Mal_high"
    mal_mask = np.isin(labels, ["Mal_high", "Mal_low"])
    info = {
        "q_low": float(q1),
        "q_high": float(q2),
        "split": split,
        "samples": info_samples,
        "n_malignant_scored": int(use_mal.sum()),
        "n_tnk_scored": int(use_tnk.sum()),
        "n_mal_high": int((labels == "Mal_high").sum()),
        "n_mal_low": int((labels == "Mal_low").sum()),
    }
    return labels, mal_mask, info


def n_table(labels, sample, origin, patient, cldn4) -> pd.DataFrame:
    rows = []
    for g in sorted(set(labels) - {"drop"}):
        idx = labels == g
        rows.append(
            {
                "group": g,
                "n_cells": int(idx.sum()),
                "n_samples": int(pd.Series(sample[idx]).nunique()),
                "n_patients": int(pd.Series(patient[idx]).nunique()),
                "origins": ",".join(sorted(pd.Series(origin[idx]).dropna().unique())),
                "mean_CLDN4_log1p_cp10k": float(np.mean(cldn4[idx])) if idx.any() else None,
            }
        )
    return pd.DataFrame(rows)


def sample_dominance(labels, sample, origin) -> pd.DataFrame:
    rows = []
    for g in sorted(set(labels) - {"drop"}):
        idx = labels == g
        vc = pd.Series(sample[idx]).value_counts()
        if vc.empty:
            continue
        top = vc.index[0]
        rows.append(
            {
                "group": g,
                "top_sample": top,
                "top_n": int(vc.iloc[0]),
                "top_frac": float(vc.iloc[0] / vc.sum()),
                "top_origin": str(pd.Series(origin[idx][np.array(sample[idx]) == top]).iloc[0]) if (np.array(sample[idx]) == top).any() else "",
            }
        )
    return pd.DataFrame(rows)


def plot_counts(n_df: pd.DataFrame, title: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.bar(n_df["group"], n_df["n_cells"], color="#4C72B0")
    ax.set_ylabel("cells")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=20)
    for i, r in n_df.iterrows():
        ax.text(
            i,
            r.n_cells,
            f"n={int(r.n_cells)}\n{int(r.n_samples)} smp / {int(r.n_patients)} pts",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_top_delta(c: pd.DataFrame, path: Path, title: str, k: int = 20) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 6.2))
    if c.empty or not c["sig_diff"].any():
        ax.text(0.5, 0.5, "no significant differential pairs", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(path, dpi=140)
        plt.close(fig)
        return
    s = c[c["sig_diff"]].copy()
    s["absd"] = s["delta_prob"].abs()
    s = s.sort_values("absd", ascending=False).head(k)
    s = s.sort_values("delta_prob")
    colors = ["#C44E52" if d > 0 else "#4C72B0" for d in s["delta_prob"]]
    labels = [f"{a} ({b})" for a, b in zip(s["interaction_name"], s["ligand_class"])]
    ax.barh(range(len(s)), s["delta_prob"], color=colors)
    ax.set_yticks(range(len(s)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("prob(CLDN4-high) − prob(CLDN4-low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def write_ligand_table(outgoing: pd.DataFrame, n_df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Primary deliverable: Mal → T/NK ligands, high vs low CLDN4."""
    if outgoing.empty:
        empty = pd.DataFrame()
        empty.to_csv(path, sep="\t", index=False)
        return empty
    n_high = int(n_df.loc[n_df.group == "Mal_high", "n_cells"].iloc[0]) if (n_df.group == "Mal_high").any() else 0
    n_low = int(n_df.loc[n_df.group == "Mal_low", "n_cells"].iloc[0]) if (n_df.group == "Mal_low").any() else 0
    n_tnk = int(n_df.loc[n_df.group == "TNK", "n_cells"].iloc[0]) if (n_df.group == "TNK").any() else 0
    n_pt_high = int(n_df.loc[n_df.group == "Mal_high", "n_patients"].iloc[0]) if (n_df.group == "Mal_high").any() else 0
    n_pt_low = int(n_df.loc[n_df.group == "Mal_low", "n_patients"].iloc[0]) if (n_df.group == "Mal_low").any() else 0
    n_pt_tnk = int(n_df.loc[n_df.group == "TNK", "n_patients"].iloc[0]) if (n_df.group == "TNK").any() else 0
    tab = outgoing.copy()
    tab["n_mal_high"] = n_high
    tab["n_mal_low"] = n_low
    tab["n_tnk"] = n_tnk
    tab["n_patients_mal_high"] = n_pt_high
    tab["n_patients_mal_low"] = n_pt_low
    tab["n_patients_tnk"] = n_pt_tnk
    tab["source_arm"] = "malignant_CLDN4_high_vs_low"
    tab["target_arm"] = "T_NK"
    tab = tab.sort_values(["sig_diff", "delta_prob"], ascending=[False, False])
    cols = [
        "interaction_name",
        "pathway_name",
        "annotation",
        "ligand",
        "receptor",
        "ligand_genes",
        "receptor_genes",
        "ligand_class",
        "direction",
        "n_mal_high",
        "n_mal_low",
        "n_tnk",
        "n_patients_mal_high",
        "n_patients_mal_low",
        "n_patients_tnk",
        "prob_high",
        "prob_low",
        "delta_prob",
        "pval_high",
        "pval_low",
        "sig_high",
        "sig_low",
        "sig_diff",
        "ligand_mean_high",
        "ligand_mean_low",
        "receptor_mean_partner",
        "ligand_prop_high",
        "ligand_prop_low",
        "receptor_prop_partner",
    ]
    tab = tab[[c for c in cols if c in tab.columns]]
    tab.to_csv(path, sep="\t", index=False)
    return tab


def run_mode(mode, expr, pos, gene_index, lr, cldn4, origin, mal, tnk, sample, patient, rng, out_dir: Path):
    print(f"== {mode} ==", flush=True)
    labels, mal_mask, info = make_split(cldn4, origin, mal, tnk, mode)
    groups = sorted(set(labels) - {"drop"})
    if "Mal_high" not in groups or "TNK" not in groups or mal_mask.sum() < MIN_GROUP * 2:
        print(f"  skip {mode}: {info}", flush=True)
        return None
    n_df = n_table(labels, sample, origin, patient, cldn4)
    n_df.to_csv(out_dir / f"n_cells_{mode}.tsv", sep="\t", index=False)
    sample_dominance(labels, sample, origin).to_csv(out_dir / f"dominance_{mode}.tsv", sep="\t", index=False)
    plot_counts(n_df, f"GSE131907 {mode} cell counts", out_dir / f"fig_n_{mode}.png")
    src_tgt = [("Mal_high", "TNK"), ("Mal_low", "TNK"), ("TNK", "Mal_high"), ("TNK", "Mal_low")]
    used = np.isin(labels, groups)
    pairs = permute_mal_high_low(
        expr[:, used],
        pos[:, used],
        gene_index,
        lr,
        labels[used],
        groups,
        src_tgt,
        mal_mask[used],
        NBOOT,
        rng,
    )
    pairs.to_csv(out_dir / f"lr_pairs_{mode}.tsv", sep="\t", index=False)
    outgoing = contrast_high_low(pairs, "Mal_high", "Mal_low", "TNK", "outgoing")
    incoming = contrast_high_low(pairs, "Mal_high", "Mal_low", "TNK", "incoming")
    if not outgoing.empty:
        outgoing.to_csv(out_dir / f"contrast_{mode}_outgoing.tsv", sep="\t", index=False)
    if not incoming.empty:
        incoming.to_csv(out_dir / f"contrast_{mode}_incoming.tsv", sep="\t", index=False)
    plot_top_delta(
        outgoing,
        out_dir / f"fig_top_outgoing_{mode}.png",
        f"Outgoing Mal→T/NK Δprob (sig only) — {mode}",
    )
    plot_top_delta(
        incoming,
        out_dir / f"fig_top_incoming_{mode}.png",
        f"Incoming T/NK→Mal Δprob (sig only) — {mode}",
    )
    ligand_tab = write_ligand_table(outgoing, n_df, out_dir / f"ligand_table_{mode}.tsv")
    n_sig_out = int(pairs.loc[pairs.source.astype(str).str.startswith("Mal_") & pairs.significant].shape[0])
    n_sig_in = int(pairs.loc[pairs.target.astype(str).str.startswith("Mal_") & pairs.significant].shape[0])
    summary = {
        "mode": mode,
        "info": info,
        "n_cells": n_df.to_dict(orient="records"),
        "n_lr_tested": int(len(pairs)),
        "n_detected": int(pairs["detected"].sum()),
        "n_significant": int(pairs["significant"].sum()),
        "n_sig_outgoing_mal_to_tnk": n_sig_out,
        "n_sig_incoming_tnk_to_mal": n_sig_in,
        "n_sig_diff_outgoing": int(outgoing["sig_diff"].sum()) if not outgoing.empty else 0,
        "n_ligand_table_rows": int(len(ligand_tab)),
        "nboot": NBOOT,
        "expr_prop": EXPR_PROP,
        "trim": TRIM,
        "kh": KH,
    }
    (out_dir / f"summary_{mode}.json").write_text(json.dumps(summary, indent=2))
    print(
        f"  {mode}: tested={len(pairs)} detected={int(pairs.detected.sum())} "
        f"sig={int(pairs.significant.sum())} out_diff={summary['n_sig_diff_outgoing']}",
        flush=True,
    )
    return summary, outgoing


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", type=Path, required=True)
    p.add_argument("--ann", type=Path, required=True)
    p.add_argument("--series", type=Path, required=True)
    p.add_argument("--db", type=Path, default=ROOT / "db")
    p.add_argument("--out", type=Path, default=ROOT / "results")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print("reading gene list", flush=True)
    gene_list = matrix_gene_list(args.matrix)
    matrix_genes = set(gene_list)
    print(f"matrix genes={len(matrix_genes)}", flush=True)
    if "CLDN4" not in matrix_genes:
        raise SystemExit("CLDN4 not in UMI matrix")

    lr_all = load_lr(args.db, matrix_genes)
    print(f"LR pairs with all subunits in matrix: {len(lr_all)}", flush=True)
    wanted = set(EXTRA)
    for rec in lr_all.itertuples(index=False):
        wanted.update(rec.ligand_genes)
        wanted.update(rec.receptor_genes)
    print(f"streaming {len(wanted)} genes", flush=True)
    cell_ids, found, n_umi, n_streamed = stream_matrix(args.matrix, wanted)
    n = len(cell_ids)
    if "CLDN4" not in found:
        raise SystemExit("CLDN4 row not stored from UMI matrix")

    ann = pd.read_csv(args.ann, sep="\t", dtype=str)
    per = ann.set_index("Index").reindex(cell_ids).reset_index()
    if per["Sample"].isna().any():
        raise SystemExit("matrix cell IDs do not align with annotation Index")
    series = parse_series_matrix(args.series)
    sample_meta = series.rename(columns={"title": "Sample"})
    keep_cols = [c for c in ["Sample", "geo_accession", "patient_id", "tumor_stage", "tissue_origin_abbrevation", "source_name_ch1"] if c in sample_meta.columns]
    sample_meta = sample_meta[keep_cols].drop_duplicates("Sample")
    sample_meta.to_csv(args.out / "sample_metadata.tsv", sep="\t", index=False)
    pmap = sample_meta.set_index("Sample")["patient_id"] if "patient_id" in sample_meta.columns else pd.Series(dtype=str)

    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
    cldn4 = log_cp["CLDN4"]
    sample = per["Sample"].to_numpy()
    origin = per["Sample_Origin"].to_numpy()
    patient = np.array([pmap[s] if s in pmap.index else s for s in sample], dtype=object)
    mal = per.apply(is_malignant, axis=1).to_numpy()
    tnk = per.apply(is_tnk, axis=1).to_numpy()

    inventory = {
        "n_barcodes": n,
        "n_samples_matrix": int(pd.Series(sample).nunique()),
        "n_patients_matrix": int(pd.Series(patient).nunique()),
        "n_author_malignant_or_tS": int(mal.sum()),
        "n_author_tnk": int(tnk.sum()),
        "n_mal_by_origin": per.loc[mal].groupby("Sample_Origin").size().to_dict() if mal.any() else {},
        "n_tnk_by_origin": per.loc[tnk].groupby("Sample_Origin").size().to_dict() if tnk.any() else {},
        "n_pe_epithelial_unlabeled": int(
            ((per.Cell_type == "Epithelial cells") & (per.Sample_Origin == "PE") & ~per.Cell_subtype.isin(MALIG_SUBTYPES)).sum()
        ),
        "ici_or_mpr_labels": False,
    }
    (args.out / "inventory.json").write_text(json.dumps(inventory, indent=2))
    print(json.dumps(inventory, indent=2), flush=True)

    lr_genes = sorted({g for rec in lr_all.itertuples(index=False) for g in rec.ligand_genes + rec.receptor_genes})
    gene_index = {g: i for i, g in enumerate(lr_genes)}
    expr = np.vstack([log_cp[g] for g in lr_genes])
    pos = np.vstack([(found[g] > 0).astype(np.float32) for g in lr_genes])

    rng = np.random.default_rng(SEED)
    modes = ["tumor_tertile", "tlung_tertile", "mets_tertile"]
    run_summaries = []
    kept_outgoing = None
    for mode in modes:
        got = run_mode(mode, expr, pos, gene_index, lr_all, cldn4, origin, mal, tnk, sample, patient, rng, args.out)
        if got is not None:
            run_summaries.append(got[0])
            if mode == "tumor_tertile":
                kept_outgoing = got[1]

    kept = next((s for s in run_summaries if s["mode"] == "tumor_tertile"), run_summaries[0] if run_summaries else None)
    if kept is not None and kept_outgoing is not None:
        n_df = pd.read_csv(args.out / f"n_cells_{kept['mode']}.tsv", sep="\t")
        write_ligand_table(kept_outgoing, n_df, args.out / "ligand_table.tsv")

    header = {
        "dataset": "GSE131907",
        "citation": "Kim et al. Nat Commun 2020, PMID 32385277",
        "matrix_cells": n,
        "matrix_genes": n_streamed,
        "lr_pairs_in_matrix": int(len(lr_all)),
        "malignant_rule": (
            "Author Epithelial cells with Cell_subtype in {Malignant cells, tS1, tS2, tS3}. "
            "tS1–tS3 are the primary tLung tumor-epithelial labels; author 'Malignant cells' "
            "are in tL/B, mLN, and mBrain. PE epithelial cells are unlabeled and excluded."
        ),
        "tnk_rule": "Author Cell_type in {T lymphocytes, NK cells}, merged as TNK",
        "primary_samples": "Tumor origins tLung + tL/B + mLN + mBrain. Normal lung/LN and PE excluded from the kept split.",
        "algorithm": (
            "CellChat-like: 10% truncated mean of log1p(CP10k), Hill Kh=0.5, "
            "CellChatDB v2 protein pairs, expr_prop>=0.10, nboot=100. "
            "CLDN4-high/low labels permuted among malignant cells; T/NK fixed."
        ),
        "not_run": "CellChat R; LIANA; CopyKAT; EGA FASTQ; 2.86 GB log2TPM text",
        "kept_split": None if kept is None else kept["mode"],
        "kept_reason": "Primary = tumor-origin author malignant/tS* CLDN4 tertile vs T/NK. tLung-only and mets-only are sensitivities.",
        "inventory": inventory,
        "runs": run_summaries,
        "kept": kept,
    }
    (args.out / "summary.json").write_text(json.dumps(header, indent=2))
    print(
        json.dumps(
            {
                "kept_split": header["kept_split"],
                "runs": [
                    (s["mode"], s["n_significant"], s["n_sig_diff_outgoing"])
                    for s in run_summaries
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
