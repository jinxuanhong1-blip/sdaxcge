#!/usr/bin/env python3
"""GSE334497 match test: does 4T1 Trop2 KO drop Cldn4 and raise IFN/immune?

Public author-normalized counts only (GEO supplementary). Contrast is
Trop2 KO minus WT, n=5 vs 5, frozen whole-tumor sections. This is the
Wu et al. JITC 2026 matrix, not a CLDN4 knockdown and not lung.

Match rule, fixed before the call is written:
  Arm A  Cldn4 is lower in KO (log2FC < 0 and one-sided Welch p < 0.05).
  Arm B  Epithelial ISG sample score is higher in KO
         (one-sided exact permutation p < 0.05).
  Joint  both arms. Bulk immune and Hallmark IFN-γ are reported beside
         Arm B so an infiltrate signal is not substituted for an ISG.
"""
from __future__ import annotations

import gzip
import hashlib
import math
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "methods" / "gse334497_cldn4_match"
DATA = OUT / "data"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

COUNTS = DATA / "GSE334497_normalized_counts.csv.gz"
SYMBOLS = DATA / "ensembl_to_symbol.tsv"

# GEO SOFT Sample_description library names (GSE334497_family.soft).
KO = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
WT = ["RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R", "control170"]

# Anchors: official mouse Ensembl gene IDs. Script refuses a bad map.
ANCHORS = {
    "ENSMUSG00000051397": "Tacstd2",
    "ENSMUSG00000047501": "Cldn4",
    "ENSMUSG00000018569": "Cldn7",
    "ENSMUSG00000022512": "Cldn1",
    "ENSMUSG00000029417": "Cxcl9",
}

# Tumor-cell antiviral program a CLDN4 KD would be asked to open.
EPITHELIAL_ISG = [
    "Isg15",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Mx1",
    "Mx2",
    "Oas2",
    "Oas3",
    "Oasl1",
    "Rsad2",
    "Stat1",
    "Irf7",
    "Ifi27",
    "Usp18",
    "Bst2",
    "Ifih1",
]

# Lineage / cytotoxicity / chemokine RNA that bulk sections pick up from infiltrate.
BULK_IMMUNE = [
    "Cd3e",
    "Cd3d",
    "Cd8a",
    "Cd8b1",
    "Cd4",
    "Gzmb",
    "Gzma",
    "Prf1",
    "Ifng",
    "Nkg7",
    "Klrd1",
    "Cxcl9",
    "Cxcl10",
]

# MSigDB Hallmark IFN-γ, mouse symbols (mh.all.v2023.2.Mm).
HALLMARK_IFNG = [
    "Adar", "Apol6", "Arid5b", "Arl4a", "Auts2", "B2m", "Bank1", "Batf2", "Bpgm",
    "Bst2", "Btg1", "C1s1", "Casp1", "Casp3", "Casp4", "Casp7", "Casp8", "Ccl5",
    "Ccl7", "Cd274", "Cd38", "Cd40", "Cd69", "Cd74", "Cd86", "Cdkn1a", "Cfb",
    "Cfh", "Ciita", "Cmklr1", "Cmpk2", "Cmtr1", "Cxcl10", "Cxcl11", "Cxcl9",
    "Ddx60", "Dhx58", "Eif2ak2", "Eif4e3", "Epsti1", "Fas", "Fcgr1", "Fgl2",
    "Fpr1", "Gbp3", "Gch1", "Gpr18", "Gzma", "H2-Aa", "H2-DMa", "Helz2", "Herc6",
    "Hif1a", "Icam1", "Ido1", "Ifi27", "Ifi30", "Ifi35", "Ifi44", "Ifi44l",
    "Ifih1", "Ifit2", "Ifit3", "Ifitm2", "Ifitm3", "Ifnar2", "Il10ra", "Il15",
    "Il15ra", "Il18bp", "Il2rb", "Il4ra", "Il6", "Il7", "Irf1", "Irf2", "Irf4",
    "Irf5", "Irf7", "Irf8", "Irf9", "Isg15", "Isg20", "Isoc1", "Itgb7", "Jak2",
    "Lap3", "Lats2", "Lcp2", "Lgals3bp", "Ly6e", "Lysmd2", "Marchf1", "Mettl7b",
    "Mthfd2", "Mvp", "Mx2", "Myd88", "Nampt", "Ncoa3", "Nfkb1", "Nfkbia", "Nlrc5",
    "Nmi", "Nod1", "Nup93", "Oas2", "Oas3", "Oasl1", "Ogfr", "P2ry14", "Parp12",
    "Parp14", "Pde4b", "Peli1", "Pfkp", "Pim1", "Pla2g4a", "Plscr1", "Pml", "Pnp",
    "Pnpt1", "Psma2", "Psma3", "Psmb10", "Psmb2", "Psmb8", "Psmb9", "Psme1",
    "Psme2", "Ptgs2", "Ptpn1", "Ptpn2", "Ptpn6", "Rapgef6", "Rbck1", "Rigi",
    "Ripk1", "Ripk2", "Rnf213", "Rnf31", "Rsad2", "Rtp4", "Samd9l", "Samhd1",
    "Sectm1a", "Selp", "Serping1", "Slamf7", "Slc25a28", "Socs1", "Socs3", "Sod2",
    "Sp110", "Sppl2a", "Sri", "Sspn", "St3gal5", "St8sia4", "Stat1", "Stat2",
    "Stat3", "Stat4", "Tap1", "Tapbp", "Tdrd7", "Tnfaip2", "Tnfaip3", "Tnfaip6",
    "Tnfsf10", "Tor1b", "Trafd1", "Trim14", "Trim21", "Trim25", "Trim26", "Txnip",
    "Ube2l6", "Upp1", "Usp18", "Vamp5", "Vamp8", "Vcam1", "Wars1", "Xaf1", "Xcl1",
    "Zbp1", "Znfx1",
]

# Infiltrate-leaning members of Hallmark IFN-γ. Not a pure sort.
LEUK_IN_IFNG = [
    "Bank1", "Ccl5", "Ccl7", "Cd274", "Cd38", "Cd40", "Cd69", "Cd74", "Cd86",
    "Ciita", "Cmklr1", "Cxcl9", "Cxcl10", "Cxcl11", "Fas", "Fcgr1", "Fpr1",
    "Gzma", "H2-Aa", "H2-DMa", "Icam1", "Il10ra", "Il15", "Il15ra", "Il2rb",
    "Il7", "Irf4", "Irf8", "Itgb7", "Lcp2", "Nlrc5", "Selp", "Slamf7", "Stat4",
    "Vcam1", "Xcl1",
]

CTRL_OXPHOS = [
    "Cox4i1", "Cox5a", "Cox6a1", "Cox7a2", "Ndufa1", "Ndufb2", "Ndufs2",
    "Sdha", "Sdhb", "Uqcrc1", "Uqcrc2", "Cycs", "Vdac1", "Idh2", "Mdh2", "Fh1",
]

FOCAL = [
    "Tacstd2",
    "Cldn1",
    "Cldn3",
    "Cldn4",
    "Cldn7",
    "Ocln",
    "Tjp1",
    "Epcam",
    "Isg15",
    "Ifit1",
    "Mx1",
    "Oas2",
    "Stat1",
    "Irf7",
    "Ifi27",
    "Cxcl9",
    "Cxcl10",
    "Cd8a",
    "Cd3e",
    "Gzmb",
    "Prf1",
    "Ifng",
    "Nkg7",
    "Cd274",
]

ALPHA = 0.05


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_symbol_map() -> dict[str, str]:
    m: dict[str, str] = {}
    with open(SYMBOLS) as f:
        header = f.readline()
        if not header.startswith("ensembl"):
            raise SystemExit(f"bad symbol header: {header!r}")
        for line in f:
            ens, sym = line.rstrip("\n").split("\t")
            m[ens] = sym
    for ens, sym in ANCHORS.items():
        if m.get(ens) != sym:
            raise SystemExit(f"anchor failed: {ens} -> {m.get(ens)!r}, expected {sym}")
    return m


def welch_full(a: np.ndarray, b: np.ndarray) -> dict:
    diff = float(a.mean() - b.mean())
    n1, n2 = len(a), len(b)
    v1, v2 = float(a.var(ddof=1)), float(b.var(ddof=1))
    se = math.sqrt(v1 / n1 + v2 / n2)
    if se == 0:
        t, df, p_two = 0.0, float(n1 + n2 - 2), 1.0
    else:
        t = diff / se
        df = (v1 / n1 + v2 / n2) ** 2 / (
            (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)
        )
        p_two = float(2 * stats.t.sf(abs(t), df))
    p_down = p_two / 2 if t < 0 else 1 - p_two / 2
    p_up = p_two / 2 if t > 0 else 1 - p_two / 2
    tcrit = float(stats.t.ppf(0.975, df))
    sp = math.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    d = 0.0 if sp == 0 else diff / sp
    try:
        _u, p_mw = stats.mannwhitneyu(a, b, alternative="two-sided")
        p_mw = float(p_mw)
    except ValueError:
        p_mw = 1.0
    return {
        "log2FC": diff,
        "welch_t": float(t),
        "df": float(df),
        "welch_p_two": p_two,
        "welch_p_down": float(p_down),
        "welch_p_up": float(p_up),
        "ci95_lo": diff - tcrit * se,
        "ci95_hi": diff + tcrit * se,
        "cohens_d": float(d),
        "mwu_p_two": p_mw,
        "mean_log2_KO": float(a.mean()),
        "mean_log2_WT": float(b.mean()),
    }


def exact_perm(scores: np.ndarray, n_ko: int, observed: float) -> dict:
    """Exact label permutation of a sample score. 5 vs 5 -> 252 splits."""
    idx = np.arange(len(scores))
    n = n_ge = n_le = n_abs = 0
    for combo in combinations(idx, n_ko):
        mask = np.zeros(len(scores), dtype=bool)
        mask[list(combo)] = True
        delta = float(scores[mask].mean() - scores[~mask].mean())
        n += 1
        if delta >= observed - 1e-12:
            n_ge += 1
        if delta <= observed + 1e-12:
            n_le += 1
        if abs(delta) >= abs(observed) - 1e-12:
            n_abs += 1
    return {
        "n_perm": n,
        "perm_p_up": n_ge / n,
        "perm_p_down": n_le / n,
        "perm_p_two": n_abs / n,
    }


def sample_zscore(log2: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [g for g in genes if g in log2.index]
    if not present:
        raise SystemExit(f"no genes present from {genes[:5]}...")
    sub = log2.loc[present]
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0).fillna(0.0)
    return z.mean(axis=0).to_numpy(), present


def fmt_p(p: float) -> str:
    if p < 0.001:
        return f"{p:.3g}"
    return f"{p:.3f}"


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    if not COUNTS.exists():
        raise SystemExit(f"missing {COUNTS}")

    raw = pd.read_csv(COUNTS, index_col=0)
    missing = [c for c in KO + WT if c not in raw.columns]
    if missing:
        raise SystemExit(f"missing columns {missing}; have {list(raw.columns)}")

    ens2sym = load_symbol_map()
    symbols = pd.Series(raw.index.astype(str).map(ens2sym), index=raw.index)
    expr = raw.copy()
    expr["symbol"] = symbols
    expr = expr.dropna(subset=["symbol"])
    # One row per symbol: keep the Ensembl ID with the higher mean count.
    expr["_mean"] = expr[KO + WT].mean(axis=1)
    expr = expr.sort_values("_mean", ascending=False)
    expr = expr.loc[~expr["symbol"].duplicated()].drop(columns="_mean")
    expr = expr.set_index("symbol")
    mat = expr[KO + WT].astype(float)
    log2 = np.log2(mat + 1.0)

    # Gene-level tests on every mapped gene, then keep a focal table plus claudins.
    gene_rows = []
    for gene in sorted(set(FOCAL) | {g for g in expr.index if str(g).startswith("Cldn")}):
        if gene not in log2.index:
            gene_rows.append({"gene": gene, "present": False})
            continue
        st = welch_full(log2.loc[gene, KO].to_numpy(), log2.loc[gene, WT].to_numpy())
        st["gene"] = gene
        st["present"] = True
        st["mean_norm_KO"] = float(mat.loc[gene, KO].mean())
        st["mean_norm_WT"] = float(mat.loc[gene, WT].mean())
        st["values_KO"] = ",".join(f"{v:.4f}" for v in log2.loc[gene, KO].tolist())
        st["values_WT"] = ",".join(f"{v:.4f}" for v in log2.loc[gene, WT].tolist())
        gene_rows.append(st)
    genes = pd.DataFrame(gene_rows)
    genes.to_csv(TABLES / "focal_genes.tsv", sep="\t", index=False)

    sets = {
        "EPITHELIAL_ISG": (EPITHELIAL_ISG, "up"),
        "BULK_IMMUNE": (BULK_IMMUNE, "up"),
        "HALLMARK_IFNG": (HALLMARK_IFNG, "up"),
        "HALLMARK_IFNG_LEUK": (LEUK_IN_IFNG, "up"),
        "HALLMARK_IFNG_ISG": ([g for g in EPITHELIAL_ISG if g in set(HALLMARK_IFNG)], "up"),
        "CTRL_OXPHOS": (CTRL_OXPHOS, "up"),
    }
    # Remainder of Hallmark IFN-γ after removing infiltrate-leaning and epithelial ISG members.
    used = set(LEUK_IN_IFNG) | set(EPITHELIAL_ISG)
    sets["HALLMARK_IFNG_REST"] = ([g for g in HALLMARK_IFNG if g not in used], "up")

    set_rows = []
    score_frame = {}
    for name, (members, _direction) in sets.items():
        scores, present = sample_zscore(log2, members)
        observed = float(scores[: len(KO)].mean() - scores[len(KO) :].mean())
        # scores are in KO+WT column order because log2 columns are KO+WT.
        perm = exact_perm(scores, len(KO), observed)
        ko_s, wt_s = scores[: len(KO)], scores[len(KO) :]
        d_stats = welch_full(ko_s, wt_s)
        set_rows.append(
            {
                "set": name,
                "n_requested": len(members),
                "n_present": len(present),
                "missing": ",".join(g for g in members if g not in present),
                "score_delta_KO_minus_WT": observed,
                "cohens_d": d_stats["cohens_d"],
                "welch_p_up": d_stats["welch_p_up"],
                "welch_p_two": d_stats["welch_p_two"],
                **perm,
                "genes_present": ",".join(present),
            }
        )
        score_frame[name] = scores
    sets_df = pd.DataFrame(set_rows)
    sets_df.to_csv(TABLES / "set_scores.tsv", sep="\t", index=False)

    score_tbl = pd.DataFrame(score_frame, index=KO + WT)
    score_tbl.insert(0, "group", ["KO"] * len(KO) + ["WT"] * len(WT))
    score_tbl.to_csv(TABLES / "sample_scores.tsv", sep="\t")

    cldn4 = genes.loc[genes["gene"] == "Cldn4"].iloc[0]
    tac = genes.loc[genes["gene"] == "Tacstd2"].iloc[0]
    arm_a = bool(cldn4["present"] and cldn4["log2FC"] < 0 and cldn4["welch_p_down"] < ALPHA)
    isg = sets_df.loc[sets_df["set"] == "EPITHELIAL_ISG"].iloc[0]
    imm = sets_df.loc[sets_df["set"] == "BULK_IMMUNE"].iloc[0]
    ifng = sets_df.loc[sets_df["set"] == "HALLMARK_IFNG"].iloc[0]
    leuk = sets_df.loc[sets_df["set"] == "HALLMARK_IFNG_LEUK"].iloc[0]
    ifng_isg = sets_df.loc[sets_df["set"] == "HALLMARK_IFNG_ISG"].iloc[0]
    rest = sets_df.loc[sets_df["set"] == "HALLMARK_IFNG_REST"].iloc[0]
    ox = sets_df.loc[sets_df["set"] == "CTRL_OXPHOS"].iloc[0]
    arm_b = bool(isg["perm_p_up"] < ALPHA and isg["score_delta_KO_minus_WT"] > 0)
    bulk_up = bool(imm["perm_p_up"] < ALPHA and imm["score_delta_KO_minus_WT"] > 0)
    ifng_up = bool(ifng["perm_p_up"] < ALPHA and ifng["score_delta_KO_minus_WT"] > 0)
    joint = arm_a and arm_b

    call = pd.DataFrame(
        [
            {
                "arm": "A_Cldn4_down",
                "rule": "log2FC<0 and one-sided Welch p<0.05",
                "stat": cldn4["log2FC"],
                "p": cldn4["welch_p_down"],
                "pass": arm_a,
            },
            {
                "arm": "B_epithelial_ISG_up",
                "rule": "sample-score delta>0 and one-sided exact perm p<0.05",
                "stat": isg["score_delta_KO_minus_WT"],
                "p": isg["perm_p_up"],
                "pass": arm_b,
            },
            {
                "arm": "report_bulk_immune_up",
                "rule": "same perm rule; not substituted for arm B",
                "stat": imm["score_delta_KO_minus_WT"],
                "p": imm["perm_p_up"],
                "pass": bulk_up,
            },
            {
                "arm": "report_hallmark_IFNg_up",
                "rule": "same perm rule; not substituted for arm B",
                "stat": ifng["score_delta_KO_minus_WT"],
                "p": ifng["perm_p_up"],
                "pass": ifng_up,
            },
            {
                "arm": "JOINT_Cldn4_and_epithelial_ISG",
                "rule": "A and B",
                "stat": float("nan"),
                "p": float("nan"),
                "pass": joint,
            },
        ]
    )
    call.to_csv(TABLES / "match_call.tsv", sep="\t", index=False)

    plot(genes, score_frame, cldn4, isg, imm, ifng, ox, joint, arm_a, arm_b, bulk_up)
    write_finding(
        cldn4, tac, genes, isg, imm, ifng, leuk, ifng_isg, rest, ox,
        arm_a, arm_b, bulk_up, ifng_up, joint,
    )
    print(call.to_string(index=False))
    print(f"counts md5 {md5(COUNTS)}")


def plot(genes, score_frame, cldn4, isg, imm, ifng, ox, joint, arm_a, arm_b, bulk_up) -> None:
    forest_genes = [
        ("Tacstd2", "QC"),
        ("Cldn1", "barrier"),
        ("Cldn4", "barrier"),
        ("Cldn7", "barrier"),
        ("Ocln", "barrier"),
        ("Isg15", "ISG"),
        ("Ifit1", "ISG"),
        ("Mx1", "ISG"),
        ("Oas2", "ISG"),
        ("Stat1", "ISG"),
        ("Cxcl9", "immune"),
        ("Cd8a", "immune"),
        ("Prf1", "immune"),
        ("Ifng", "immune"),
        ("Gzmb", "immune"),
    ]
    colors = {"QC": "#6b6b6b", "barrier": "#8c6d31", "ISG": "#2a9d8f", "immune": "#c44536"}
    g = genes.set_index("gene")

    fig = plt.figure(figsize=(10.2, 7.2), dpi=160)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.38, wspace=0.32)
    ax0 = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])

    ys = np.arange(len(forest_genes))[::-1]
    for y, (gene, kind) in zip(ys, forest_genes):
        row = g.loc[gene]
        ax0.plot([row["ci95_lo"], row["ci95_hi"]], [y, y], color=colors[kind], lw=1.6)
        ax0.plot(row["log2FC"], y, "o", color=colors[kind], ms=6, zorder=3)
    ax0.axvline(0, color="#888", lw=0.8)
    ax0.set_yticks(ys)
    ax0.set_yticklabels([gname for gname, _ in forest_genes])
    ax0.set_xlabel("log2FC (Trop2 KO − WT), Welch 95% CI")
    ax0.set_title("Gene effects on the public 5 vs 5 matrix")
    handles = [
        plt.Line2D([0], [0], marker="o", color=c, label=lab, lw=0)
        for lab, c in [("Tacstd2 QC", colors["QC"]), ("Barrier", colors["barrier"]),
                       ("Epithelial ISG", colors["ISG"]), ("Bulk immune", colors["immune"])]
    ]
    ax0.legend(handles=handles, frameon=False, loc="lower right", fontsize=8)

    # Cldn4 points. values are stored comma-separated log2.
    ko_vals = np.array([float(x) for x in str(cldn4["values_KO"]).split(",")])
    wt_vals = np.array([float(x) for x in str(cldn4["values_WT"]).split(",")])
    rng = np.random.default_rng(0)
    ax1.scatter(rng.normal(0, 0.04, size=len(wt_vals)), wt_vals, color="#4c78a8", s=36, zorder=3, label="WT")
    ax1.scatter(rng.normal(1, 0.04, size=len(ko_vals)), ko_vals, color="#e45756", s=36, zorder=3, label="Trop2 KO")
    ax1.legend(frameon=False, loc="lower right", fontsize=8)
    ax1.hlines(wt_vals.mean(), -0.2, 0.2, color="#4c78a8", lw=2)
    ax1.hlines(ko_vals.mean(), 0.8, 1.2, color="#e45756", lw=2)
    ax1.set_xticks([0, 1], ["WT", "Trop2 KO"])
    ax1.set_xlim(-0.45, 1.45)
    ax1.set_ylabel("Cldn4, log2(norm + 1)")
    ax1.set_title(
        f"Cldn4  log2FC {cldn4['log2FC']:+.2f}   one-sided p {fmt_p(cldn4['welch_p_down'])}"
    )

    order = [
        ("EPITHELIAL_ISG", "Epithelial\nISG", isg),
        ("HALLMARK_IFNG", "Hallmark\nIFN-γ", ifng),
        ("BULK_IMMUNE", "Bulk\nimmune", imm),
        ("CTRL_OXPHOS", "OXPHOS\ncontrol", ox),
    ]
    all_sc = np.concatenate([score_frame[k] for k, _lab, _row in order])
    y_lo = float(all_sc.min()) - 0.15
    y_hi = float(all_sc.max()) + 0.55
    ax2.set_ylim(y_lo, y_hi)
    for i, (key, label, row) in enumerate(order):
        sc = score_frame[key]
        ax2.scatter(np.full(5, i) + rng.normal(0, 0.04, 5), sc[5:], color="#4c78a8", s=28, zorder=3, label="WT" if i == 0 else None)
        ax2.scatter(np.full(5, i) + rng.normal(0, 0.04, 5), sc[:5], color="#e45756", s=28, zorder=3, label="Trop2 KO" if i == 0 else None)
        ax2.hlines(sc[5:].mean(), i - 0.18, i + 0.18, color="#4c78a8", lw=1.8)
        ax2.hlines(sc[:5].mean(), i - 0.18, i + 0.18, color="#e45756", lw=1.8)
        ax2.text(
            i,
            y_hi - 0.06,
            f"p↑ {fmt_p(row['perm_p_up'])}",
            ha="center",
            va="top",
            fontsize=8,
            color="#222",
        )
    ax2.set_xticks(range(4), [lab for _k, lab, _r in order])
    ax2.set_ylabel("Sample score (mean gene z)")
    ax2.set_title("Exact 252-split permutation, one-sided up in KO")
    ax2.axhline(0, color="#bbb", lw=0.6)
    ax2.legend(frameon=False, loc="lower left", fontsize=8)

    verdict = "JOINT MATCH" if joint else "NO JOINT MATCH"
    detail = (
        f"Cldn4 arm {'pass' if arm_a else 'does not pass'}"
        f"  ·  epithelial ISG {'pass' if arm_b else 'does not pass'}"
        f"  ·  bulk immune {'up' if bulk_up else 'not up'}"
    )
    fig.suptitle(
        f"GSE334497 4T1 Trop2 KO — {verdict}\n{detail}",
        fontsize=12,
        y=1.02,
    )
    fig.savefig(FIGS / "fig_match.png", bbox_inches="tight")
    fig.savefig(FIGS / "fig_match.pdf", bbox_inches="tight")
    plt.close(fig)


def gene_line(genes: pd.DataFrame, name: str, side: str = "down") -> str:
    row = genes.loc[genes["gene"] == name]
    if row.empty or not bool(row.iloc[0]["present"]):
        return f"| {name} | absent | | | | |"
    r = row.iloc[0]
    p_one = r["welch_p_down"] if side == "down" else r["welch_p_up"]
    return (
        f"| {name} | {r['log2FC']:+.3f} | {r['ci95_lo']:+.2f}, {r['ci95_hi']:+.2f} | "
        f"{fmt_p(r['welch_p_two'])} | {fmt_p(p_one)} | {r['cohens_d']:+.2f} |"
    )


def write_finding(cldn4, tac, genes, isg, imm, ifng, leuk, ifng_isg, rest, ox,
                  arm_a, arm_b, bulk_up, ifng_up, joint) -> None:
    def _g(name: str) -> str:
        r = genes.loc[genes["gene"] == name].iloc[0]
        return (
            f"log2FC {r['log2FC']:+.3f}, two-sided *p* = {fmt_p(r['welch_p_two'])}, "
            f"one-sided down *p* = {fmt_p(r['welch_p_down'])}, MWU *p* = {fmt_p(r['mwu_p_two'])}"
        )

    def set_line(label, row, passed):
        return (
            f"| {label} | {int(row['n_present'])}/{int(row['n_requested'])} | "
            f"{row['score_delta_KO_minus_WT']:+.3f} | {row['cohens_d']:+.2f} | "
            f"{fmt_p(row['perm_p_up'])} | {fmt_p(row['perm_p_two'])} | "
            f"{'yes' if passed else 'no'} |"
        )

    claudin_names = [g for g in genes["gene"] if str(g).startswith("Cldn") and bool(
        genes.loc[genes["gene"] == g, "present"].iloc[0]
    )]
    # sort claudins by log2FC
    claudin_names = sorted(
        claudin_names,
        key=lambda g: float(genes.loc[genes["gene"] == g, "log2FC"].iloc[0]),
    )
    claudin_table = "\n".join(gene_line(genes, g, "down") for g in claudin_names)

    text = f"""# FINDING — GSE334497 Trop2 KO vs a CLDN4-KD direction

**Additive public evidence.** Wu *et al.*, *Journal for ImmunoTherapy of Cancer* 2026;14:e012265 ([JITC](https://jitc.bmj.com/content/14/4/e012265)). GEO [GSE334497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334497). CRISPR **Trop2 (Tacstd2)** knockout versus wild-type **4T1** tumors, 5 vs 5, grown 3 weeks in BALB/c, RNA from **frozen whole-tumor sections**. Breast, not lung. This is the depositing paper’s matrix, not an independent cohort, and it is **not a CLDN4 knockdown**.

## Answer

**No.** Trop2 KO does not match a CLDN4 knockdown. *Cldn4* is lower on average, and the 5-vs-5 test does not clear 0.05. The epithelial interferon genes stay flat. Bulk immune RNA does rise.

| Question | Result |
|---|---|
| Did Trop2 KO work? | Yes. *Tacstd2* log2FC **{tac['log2FC']:+.3f}**, one-sided Welch *p* = **{fmt_p(tac['welch_p_down'])}**. |
| Does *Cldn4* drop? | Direction only. log2FC **{cldn4['log2FC']:+.3f}** (95% CI {cldn4['ci95_lo']:+.2f} to {cldn4['ci95_hi']:+.2f}), one-sided Welch *p* = **{fmt_p(cldn4['welch_p_down'])}**, two-sided *p* = {fmt_p(cldn4['welch_p_two'])}, MWU *p* = {fmt_p(cldn4['mwu_p_two'])}. Arm A {'passes' if arm_a else 'does not pass'}. |
| Does epithelial IFN rise? | **No.** Epithelial ISG Δ **{isg['score_delta_KO_minus_WT']:+.3f}**, exact one-sided perm *p* = **{fmt_p(isg['perm_p_up'])}**. Arm B {'passes' if arm_b else 'does not pass'}. |
| Does the full Hallmark IFN-γ score rise? | **Not at the sample level.** Δ **{ifng['score_delta_KO_minus_WT']:+.3f}**, *d* = {ifng['cohens_d']:+.2f}, exact perm *p* = **{fmt_p(ifng['perm_p_up'])}**. |
| What part of IFN-γ does rise? | The infiltrate-leaning slice. Δ **{leuk['score_delta_KO_minus_WT']:+.3f}**, perm *p* = **{fmt_p(leuk['perm_p_up'])}**. The ISG slice inside the same hallmark does not (Δ {ifng_isg['score_delta_KO_minus_WT']:+.3f}, *p* = {fmt_p(ifng_isg['perm_p_up'])}). |
| Does bulk immune RNA rise? | **Yes.** Δ **{imm['score_delta_KO_minus_WT']:+.3f}**, *d* = {imm['cohens_d']:+.2f}, exact perm *p* = **{fmt_p(imm['perm_p_up'])}**. |

Joint rule (both required): *Cldn4* down at one-sided *p* < 0.05, and epithelial ISG up at one-sided exact perm *p* < 0.05. Joint = **{'yes' if joint else 'no'}**.

A prerank gene-set statistic can still put Hallmark IFN-γ at the KO end of the ranking, because *Cxcl9*, MHC-II, and other infiltrate genes sit at that end. That is a statement about gene order. It is not the same as the 10 tumors separating on the mean of all ~184 Hallmark genes. The pre-specified set test here is the exact sample permutation. On that test the full Hallmark score does not pass 0.05, and the epithelial ISG score is null. The score that passes is the bulk immune module, which is the paper’s T-cell infiltration result read out of the same frozen sections. OXPHOS, run with the same permutation, stays null.

## Design

| Item | Choice |
|---|---|
| Matrix | `GSE334497_normalized_counts.csv.gz` (author-normalized counts). No FASTQ. |
| Groups | KO: KO162, KO164, KO165, KO172, RESUB-KO163R. WT: RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R, control170. Library names are GEO `Sample_description`. |
| Transform | log2(normalized count + 1). Ensembl → NCBI symbol (`gene2ensembl` + `Mus_musculus.gene_info`). Duplicate symbols: keep the row with the higher mean count. |
| Gene test | Welch *t* on the 5 vs 5 log2 values. One-sided *p* is the pre-specified arm; two-sided and Mann–Whitney are reported with it. |
| Set test | Per-sample mean of gene-wise z-scores. Exact permutation of the 252 equal-sized label splits. One-sided up = KO higher. |
| Arm B genes | Isg15, Ifit1, Ifit2, Ifit3, Mx1, Mx2, Oas2, Oas3, Oasl1, Rsad2, Stat1, Irf7, Ifi27, Usp18, Bst2, Ifih1. |
| Bulk immune | Cd3e, Cd3d, Cd8a, Cd8b1, Cd4, Gzmb, Gzma, Prf1, Ifng, Nkg7, Klrd1, Cxcl9, Cxcl10. |
| Control | OXPHOS gene score, same permutation, expected null. |
| α | 0.05, pre-specified arms only. No genome-wide FDR claim. |

## Perturbation QC and the Cldn4 arm

| Gene | log2FC | 95% CI | two-sided *p* | one-sided down *p* | *d* |
|---|---:|---|---:|---:|---:|
{gene_line(genes, 'Tacstd2')}
{gene_line(genes, 'Cldn4')}
{gene_line(genes, 'Cldn1')}
{gene_line(genes, 'Cldn7')}
{gene_line(genes, 'Ocln')}
{gene_line(genes, 'Tjp1')}
{gene_line(genes, 'Epcam')}

*Cldn4* is abundant on both sides (mean normalized count KO {cldn4['mean_norm_KO']:.0f}, WT {cldn4['mean_norm_WT']:.0f}). The ten samples overlap, and the Welch interval crosses zero. The paper’s stated RNA mediators are **claudin 1, claudin 7, and occludin**, not claudin 4. *Cldn1* is the large drop ({_g('Cldn1')}). *Cldn7* has almost the same log2FC as *Cldn4* and separates more cleanly on the rank test ({_g('Cldn7')}). *Ocln* is smaller and one-sided only ({_g('Ocln')}). The family table is below; the match arm is *Cldn4* only.

## IFN / immune scores

Positive Δ = higher in Trop2 KO. Permutation *p* is exact (252 splits).

| Set | n present | Δ | *d* | perm *p* up | perm *p* two | passes up |
|---|---:|---:|---:|---:|---:|---|
{set_line('Epithelial ISG (arm B)', isg, arm_b)}
{set_line('Hallmark IFN-γ', ifng, ifng_up)}
{set_line('Hallmark IFN-γ, infiltrate-leaning', leuk, leuk['perm_p_up'] < 0.05 and leuk['score_delta_KO_minus_WT'] > 0)}
{set_line('Hallmark IFN-γ ∩ epithelial ISG', ifng_isg, ifng_isg['perm_p_up'] < 0.05 and ifng_isg['score_delta_KO_minus_WT'] > 0)}
{set_line('Hallmark IFN-γ remainder', rest, rest['perm_p_up'] < 0.05 and rest['score_delta_KO_minus_WT'] > 0)}
{set_line('Bulk immune', imm, bulk_up)}
{set_line('OXPHOS control', ox, ox['perm_p_up'] < 0.05 and ox['score_delta_KO_minus_WT'] > 0)}

Epithelial ISG genes, log2FC (KO − WT):

| Gene | log2FC | 95% CI | two-sided *p* | one-sided up *p* | *d* |
|---|---:|---|---:|---:|---:|
{gene_line(genes, 'Isg15', 'up')}
{gene_line(genes, 'Ifit1', 'up')}
{gene_line(genes, 'Mx1', 'up')}
{gene_line(genes, 'Oas2', 'up')}
{gene_line(genes, 'Stat1', 'up')}
{gene_line(genes, 'Irf7', 'up')}
{gene_line(genes, 'Ifi27', 'up')}

Bulk immune genes:

| Gene | log2FC | 95% CI | two-sided *p* | one-sided up *p* | *d* |
|---|---:|---|---:|---:|---:|
{gene_line(genes, 'Cxcl9', 'up')}
{gene_line(genes, 'Cxcl10', 'up')}
{gene_line(genes, 'Cd8a', 'up')}
{gene_line(genes, 'Cd3e', 'up')}
{gene_line(genes, 'Gzmb', 'up')}
{gene_line(genes, 'Prf1', 'up')}
{gene_line(genes, 'Ifng', 'up')}
{gene_line(genes, 'Nkg7', 'up')}
{gene_line(genes, 'Cd274', 'up')}

*Cxcl9* is the clearest single immune gene. *Cd274* (PD-L1) is in the table because it sits in Hallmark IFN-γ; it is not a pre-specified arm.

## Claudin family (every *Cldn* present)

| Gene | log2FC | 95% CI | two-sided *p* | one-sided down *p* | *d* |
|---|---:|---|---:|---:|---:|
{claudin_table}

## What this is and is not

**Is**

- A pre-specified match test of “Cldn4 down and epithelial IFN up” on the public 4T1 Trop2 KO counts.
- A split of Hallmark IFN-γ. The infiltrate-leaning slice is higher in KO. The epithelial ISG slice is flat. The full 184-gene score falls between them and does not pass 0.05.
- Consistent with Wu *et al.*: Trop2 loss goes with less tight-junction RNA and more T-cell / inflammatory RNA in the bulk tumor.

**Is not**

- A CLDN4 knockdown or knockout.
- A lung model.
- Tumor-cell-intrinsic IFN (no sorted epithelium and no in-vitro 4T1 arm in this accession).
- An independent replication of the JITC paper. These are their tumors.
- A genome-wide discovery list. At n = 5 vs 5, gene-level tests are the pre-specified ones above.

Figure: `figures/fig_match.png`.

Reproduce: `python3 scripts/gse334497_cldn4_match/analyze.py`

## Exploratory sweep (not the pre-specified call)

The pre-specified arms above stay as written: *Cldn4* the gene and the epithelial ISG score do not clear 0.05. A separate grid (`scripts/gse334497_cldn4_match/sweep.py`, write-up `SWEEP.md`) varied the DE method, the gene set, the sample filter, and a *Cldn4* high/low contrast. This is still **Trop2 KO, not a CLDN4 knockdown**. Numbers below are copied from that grid’s all-10-tumor rows.

On all 10 tumors, sample permutation and gene-rank agree for:

| Program | Mean log2FC | Sample perm *p* | Gene-rank *p* |
|---|---:|---:|---:|
| KEGG tight junction down | −0.249 | 0.004 | 5.07×10⁻⁸ |
| Claudin/TJ core down | −0.675 | 0.016 | 8.34×10⁻⁸ |
| Keratinization down | −2.066 | 0.004 | 1.92×10⁻¹⁶ |
| Bulk immune up | +0.648 | 0.008 | 4.99×10⁻⁷ |

Hallmark IFN-γ, IFN-α, and MHC-I/APM are up on the gene-rank test (p = 2.37×10⁻¹², 6.76×10⁻⁶, 3.34×10⁻⁴) and are not up on the sample permutation (p = 0.107, 0.214, 0.175). STING core is a small shift (mean log2FC +0.036, sample p = 0.048, grid FDR about 0.06–0.07); *Cgas* is down and *Irf3* / *Sting1* are up. NHEJ core does not fall.

Figures: `figures/fig_sweep_grid.png`, `figures/fig_sweep_best.png`.
"""
    (OUT / "FINDING.md").write_text(text)
    print("wrote", OUT / "FINDING.md")


if __name__ == "__main__":
    main()
