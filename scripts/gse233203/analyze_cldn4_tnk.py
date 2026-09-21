#!/usr/bin/env python3
"""GSE233203 pleural-fluid LUAD: malignant CLDN4 vs T/NK and ABCP response.

Seven stage-IV adenocarcinoma effusions. Therapeutic response is the GEO label
for atezolizumab + bevacizumab + carboplatin/paclitaxel. This is not solid
tumor and is not merged into concordant-4. n=7 Spearman is reported as-is.
"""
from __future__ import annotations

import gzip
import json
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "GSE233203"
OUT = ROOT / "results" / "gse233203"

SAMPLES = [
    ("GSM7412612_NCCLu_162", "NCCLu_162", "Response"),
    ("GSM7412613_NCCLu_185", "NCCLu_185", "Non-response"),
    ("GSM7412614_NCCLu_327", "NCCLu_327", "Non-response"),
    ("GSM7412615_NCCLu_334", "NCCLu_334", "Non-response"),
    ("GSM7412616_NCCLu_376", "NCCLu_376", "Response"),
    ("GSM7412617_NCCLu_383", "NCCLu_383", "Response"),
    ("GSM7412618_NCCLu_397", "NCCLu_397", "Non-response"),
]
LINEAGE = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "MS4A1", "CD19"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "Fibroblast": ["COL1A1", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
}
NORMAL = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3", "FOXJ1"]
EXTRA = ["CLDN4", "TACSTD2", "CD8A", "PTPRC"]


def wanted() -> list[str]:
    genes = set(EXTRA) | set(NORMAL)
    for vs in LINEAGE.values():
        genes.update(vs)
    return sorted(genes)


def read_10x(prefix: Path):
    mtx = scipy.io.mmread(str(prefix) + "_matrix.mtx.gz").tocsr().astype(np.float32)
    with gzip.open(str(prefix) + "_barcodes.tsv.gz", "rt") as fh:
        barcodes = [line.strip().split("\t")[0] for line in fh]
    genes = []
    with gzip.open(str(prefix) + "_features.tsv.gz", "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    if mtx.shape[0] != len(genes) and mtx.shape[1] == len(genes):
        mtx = mtx.T.tocsr()
    if mtx.shape[0] != len(genes) or mtx.shape[1] != len(barcodes):
        raise ValueError(f"{prefix}: mtx {mtx.shape} genes {len(genes)} barcodes {len(barcodes)}")
    return mtx, np.array(genes), np.array(barcodes)


def gene_index(genes: np.ndarray) -> dict[str, int]:
    idx = {}
    for i, g in enumerate(genes):
        if g not in idx:
            idx[g] = i
    return idx


def extract(mtx, gidx, names, n_umi):
    scale = np.divide(1e4, n_umi, out=np.zeros_like(n_umi), where=n_umi > 0)
    out = {}
    for g in names:
        if g not in gidx:
            continue
        raw = np.asarray(mtx[gidx[g], :].todense()).ravel()
        out[g] = np.log1p(raw * scale).astype(np.float32)
    return out


def module(log_cp, genes, n):
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores, cd3):
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    labels = np.array(names, dtype=object)[best].copy()
    labels[top < 0.12] = "Unassigned"
    t_idx, nk_idx = names.index("T"), names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk = np.isin(labels, ["T", "NK"])
    labels[close & both & tnk & (cd3 > 0.15)] = "T"
    labels[close & both & tnk & (cd3 <= 0.15)] = "NK"
    return labels


def exact_wilcoxon(values, is_response) -> dict:
    values = np.asarray(values, float)
    is_response = np.asarray(is_response, bool)
    ok = np.isfinite(values)
    values, is_response = values[ok], is_response[ok]
    n_r = int(is_response.sum())
    n_n = int((~is_response).sum())
    if n_r < 1 or n_n < 1:
        return {"n_response": n_r, "n_nonresponse": n_n, "p_exact": None}
    obs = float(stats.mannwhitneyu(values[is_response], values[~is_response], alternative="two-sided").statistic)
    expected = n_r * n_n / 2.0
    ext = abs(obs - expected)
    count = total = 0
    for combo in combinations(range(len(values)), n_r):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(combo)] = True
        u = float(stats.mannwhitneyu(values[mask], values[~mask], alternative="two-sided").statistic)
        if abs(u - expected) >= ext - 1e-12:
            count += 1
        total += 1
    return {
        "n_response": n_r,
        "n_nonresponse": n_n,
        "mean_response": float(values[is_response].mean()),
        "mean_nonresponse": float(values[~is_response].mean()),
        "delta_response_minus_nonresponse": float(values[is_response].mean() - values[~is_response].mean()),
        "U": obs,
        "p_exact": count / total,
        "n_perm": total,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tar = DATA / "GSE233203_RAW.tar"
    if tar.exists() and not list(DATA.glob("GSM*_matrix.mtx.gz")):
        import tarfile
        with tarfile.open(tar) as tf:
            tf.extractall(DATA)
    keep = wanted()
    rows = []
    for prefix, sample, response in SAMPLES:
        print("reading", prefix, flush=True)
        mtx, genes, _barcodes = read_10x(DATA / prefix)
        gidx = gene_index(genes)
        n_umi = np.asarray(mtx.sum(axis=0)).ravel()
        n_genes = np.asarray((mtx > 0).sum(axis=0)).ravel()
        mt_idx = [i for i, g in enumerate(genes) if str(g).startswith("MT-")]
        mt = np.asarray(mtx[mt_idx, :].sum(axis=0)).ravel() if mt_idx else np.zeros(mtx.shape[1])
        pct_mt = np.where(n_umi > 0, 100.0 * mt / n_umi, 0.0)
        qc = (n_genes >= 200) & (pct_mt < 20) & (n_genes < 8000)
        mtx = mtx[:, qc]
        n_umi = n_umi[qc]
        n = mtx.shape[1]
        log_cp = extract(mtx, gidx, keep, n_umi)
        scores = {name: module(log_cp, gs, n) for name, gs in LINEAGE.items()}
        cd3 = log_cp.get("CD3E", log_cp.get("CD3D", np.zeros(n)))
        lineage = assign_lineage(scores, cd3)
        normal = module(log_cp, NORMAL, n)
        epi = lineage == "Epithelial"
        if int(epi.sum()) >= 20:
            cut = float(np.quantile(normal[epi], 0.75))
        else:
            cut = 0.3
        malignant = epi & (normal <= cut)
        cld = log_cp.get("CLDN4", np.zeros(n))
        cld_raw_pos = cld > 0
        tnk = np.isin(lineage, ["T", "NK"])
        rows.append(
            {
                "sample": sample,
                "response": response,
                "n_qc": n,
                "n_epithelial": int(epi.sum()),
                "n_malignant": int(malignant.sum()),
                "n_tnk": int(tnk.sum()),
                "tnk_frac": float(tnk.mean()) if n else np.nan,
                "cldn4_mean": float(cld[malignant].mean()) if malignant.any() else np.nan,
                "cldn4_pct": float(cld_raw_pos[malignant].mean()) if malignant.any() else np.nan,
                "tacstd2_mean": float(log_cp.get("TACSTD2", np.zeros(n))[malignant].mean()) if malignant.any() else np.nan,
                "normal_q75_cut": cut,
            }
        )
        print(rows[-1], flush=True)
    per = pd.DataFrame(rows)
    per.to_csv(OUT / "per_sample.tsv", sep="\t", index=False)
    usable = per[per["n_malignant"] >= 20].copy()
    x = usable["cldn4_mean"].to_numpy()
    y = usable["tnk_frac"].to_numpy()
    if len(usable) >= 4:
        rho, p = stats.spearmanr(x, y)
    else:
        rho, p = np.nan, np.nan
    rho_pct, p_pct = (stats.spearmanr(usable["cldn4_pct"], usable["tnk_frac"]) if len(usable) >= 4 else (np.nan, np.nan))
    resp = exact_wilcoxon(usable["cldn4_mean"], usable["response"].eq("Response"))
    resp_pct = exact_wilcoxon(usable["cldn4_pct"], usable["response"].eq("Response"))
    # Leave-one-out on the Spearman.
    loo = []
    if len(usable) >= 5:
        for sample in usable["sample"]:
            sub = usable[usable["sample"] != sample]
            r, pp = stats.spearmanr(sub["cldn4_mean"], sub["tnk_frac"])
            loo.append({"dropped": sample, "n": int(len(sub)), "rho": float(r), "p": float(pp)})
    loo_df = pd.DataFrame(loo)
    if len(loo_df):
        loo_df.to_csv(OUT / "leave_one_out.tsv", sep="\t", index=False)
    summary = {
        "series": "GSE233203",
        "tissue": "pleural fluid, stage IV lung adenocarcinoma",
        "treatment": "atezolizumab + bevacizumab + carboplatin/paclitaxel",
        "n_samples": int(len(per)),
        "n_usable_malignant_ge_20": int(len(usable)),
        "spearman_cldn4_mean_vs_tnk": {"n": int(len(usable)), "rho": None if rho != rho else float(rho), "p": None if p != p else float(p)},
        "spearman_cldn4_pct_vs_tnk": {
            "n": int(len(usable)),
            "rho": None if rho_pct != rho_pct else float(rho_pct),
            "p": None if p_pct != p_pct else float(p_pct),
        },
        "cldn4_mean_response_vs_nonresponse": resp,
        "cldn4_pct_response_vs_nonresponse": resp_pct,
        "not_merged_into_concordant4": True,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    for response, g in usable.groupby("response"):
        color = "#1b4f72" if response == "Response" else "#a04000"
        ax.scatter(g["tnk_frac"], g["cldn4_mean"], s=55, c=color, label=response, zorder=3)
        for _, row in g.iterrows():
            ax.annotate(row["sample"].replace("NCCLu_", ""), (row["tnk_frac"], row["cldn4_mean"]), fontsize=7, alpha=0.8)
    rho_s = "NA" if rho != rho else f"{rho:.2f}"
    p_s = "NA" if p != p else f"{p:.3g}"
    ax.set_title(f"GSE233203 effusion\nCLDN4 vs T/NK ρ={rho_s}, p={p_s}, n={len(usable)}")
    ax.set_xlabel("T/NK fraction")
    ax.set_ylabel("Malignant CLDN4 mean log1p(CP10k)")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_cldn4_vs_tnk.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
