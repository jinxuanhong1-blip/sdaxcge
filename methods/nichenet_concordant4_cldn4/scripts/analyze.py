#!/usr/bin/env python3
"""Patient-aware NicheNet summary for concordant-4.

Sender split (primary): within each patient, malignant cells in the top vs bottom
quartile of CLDN4 log1p(CP10k), rank ties.method = first.
Receiver geneset: T/NK genes whose patient-level mean tracks malignant CLDN4 %pos,
DerSimonian-Laird across cohorts. NicheNet activity is the sourced nichenetr
predict_ligand_activities (see activity.R).
"""

from __future__ import annotations

import json
import os
import subprocess
from itertools import combinations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata, wilcoxon

HERE = os.environ.get("NICHENET_HERE", "/workspace/methods/nichenet_concordant4_cldn4")
EXT = os.environ.get("NICHENET_EXTRACT", "/tmp/nichenet_work/extract")
PRIOR = os.environ.get("NICHENET_PRIOR", "/tmp/nichenet_prior")
SRC = os.environ.get("NICHENET_SRC", "/tmp/nichenet_src")
TAB = os.path.join(HERE, "results", "tables")
FIG = os.path.join(HERE, "figures")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
ALLOWED = set(COHORTS)

CORE_BARRIER = ["F11R", "NECTIN2", "CDH1", "LGALS9"]
EXT_BARRIER = [
    "PVR", "NECTIN1", "NECTIN3", "CD274", "PDCD1LG2", "CEACAM1",
    "TGFB1", "CD47", "HLA-E", "MIF", "CD24", "LGALS3",
]
CORE_IFN = ["CXCL9", "CXCL10", "CCL5", "IFNG"]
EXT_IFN = ["CXCL11", "CXCL16", "CCL2", "CCL4", "CXCL12", "ICAM1", "TNFSF9", "IL15", "IL18"]
MHC = ["HLA-A", "HLA-B", "HLA-C"]

# Receiver programs. A priori, same lists as the earlier NicheNet gene-set file.
IFN_PROGRAM = [
    "IFNG", "STAT1", "IRF1", "IRF7", "ISG15", "MX1", "IFIT1", "IFIT2", "IFIT3",
    "OAS1", "IFI6", "RSAD2", "CXCL9", "CXCL10", "CXCL11", "GBP1",
]
CYTO_PROGRAM = [
    "GZMB", "GZMA", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7", "IFNG", "FASLG",
    "TNF", "CST7", "FGFBP2", "KLRK1", "GZMM",
]
EXHAUST_PROGRAM = [
    "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CTLA4", "ENTPD1", "LAYN",
    "CXCL13", "CD38", "CD244", "TOX2", "CD160", "BTLA",
]
# Epithelial genes that can leak into a marker T/NK gate. Removed from receiver genesets.
# Lung epithelial / secretory / ciliated genes. Removed from the receiver geneset
# so a doublet or a leaky T/NK gate cannot define the NicheNet response.
EPI_LEAK = {
    "EPCAM", "TACSTD2", "CDH1", "NKX2-1", "FOXA2", "GATA6", "ELF3", "HNF1B", "SOX2",
    "NAPSA", "AGR2", "MUC1", "MUC4", "MUC5B", "MUC16", "CEACAM5", "CEACAM6",
    "WFDC2", "SCGB1A1", "SCGB3A1", "SCGB3A2", "BPIFA1", "BPIFB1",
    "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "SFTA2", "SFTA3",
    "AGER", "AQP4", "AQP1", "FOXJ1", "TMC5", "ARMC3", "EMP2", "RAB25", "GRB7",
    "KRT5", "KRT7", "KRT8", "KRT17", "KRT18", "KRT19", "CLDN3", "CLDN4", "CLDN7",
}
AUTHOR_COHORTS = ["GSE131907", "GSE205335"]

BAR_COLOR = "#1f4e79"
IFN_COLOR = "#b85c38"
MHC_COLOR = "#5c6b73"


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    tmp = np.empty(n)
    tmp[order] = np.clip(q, 0, 1)
    out[np.flatnonzero(ok)] = tmp
    return out


def dl_meta(effects, ses):
    effects = np.asarray(effects, dtype=float)
    ses = np.asarray(ses, dtype=float)
    ok = np.isfinite(effects) & np.isfinite(ses) & (ses > 0)
    effects, ses = effects[ok], ses[ok]
    k = effects.size
    if k < 2:
        return dict(mu=np.nan, se=np.nan, p=np.nan, I2=np.nan, k=int(k), Q=np.nan)
    w = 1.0 / ses ** 2
    mu_fe = np.sum(w * effects) / np.sum(w)
    Q = float(np.sum(w * (effects - mu_fe) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - (k - 1)) / c) if c > 0 else 0.0
    w2 = 1.0 / (ses ** 2 + tau2)
    mu = float(np.sum(w2 * effects) / np.sum(w2))
    se = float(np.sqrt(1.0 / np.sum(w2)))
    z = mu / se if se > 0 else np.nan
    p = float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
    I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    return dict(mu=mu, se=se, p=p, I2=I2, k=int(k), Q=Q)


def fisher_z_meta(rhos, ns):
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    ok = np.isfinite(rhos) & (ns >= 6)
    rhos, ns = rhos[ok], ns[ok]
    if rhos.size == 0:
        return dict(mu=np.nan, se=np.nan, p=np.nan, I2=np.nan, k=0, Q=np.nan, rho=np.nan)
    z = np.arctanh(np.clip(rhos, -0.999, 0.999))
    se = 1.0 / np.sqrt(ns - 3)
    meta = dl_meta(z, se)
    meta["rho"] = float(np.tanh(meta["mu"])) if np.isfinite(meta["mu"]) else np.nan
    return meta


def wilcox_greater_zero(x: np.ndarray):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    nz = x[x != 0]
    if nz.size < 5:
        return dict(n=int(n), n_nonzero=int(nz.size), mean=float(np.mean(x)) if n else np.nan,
                    median=float(np.median(x)) if n else np.nan, p=np.nan, stat=np.nan)
    # R wilcox.test one-sample drops zeros
    st = wilcoxon(nz, alternative="two-sided", zero_method="wilcox", method="auto")
    return dict(n=int(n), n_nonzero=int(nz.size), mean=float(np.mean(x)),
                median=float(np.median(x)), p=float(st.pvalue), stat=float(st.statistic))


def load_units_ligands():
    units, ligs = [], []
    for c in COHORTS:
        up = os.path.join(EXT, f"units_{c}.tsv")
        lp = os.path.join(EXT, f"ligands_{c}.tsv")
        if not os.path.exists(up):
            raise SystemExit(f"missing extract {up}")
        u = pd.read_csv(up, sep="\t")
        if not set(u.cohort).issubset(ALLOWED):
            raise SystemExit("discordant cohort in units")
        units.append(u)
        if os.path.exists(lp) and os.path.getsize(lp) > 10:
            ligs.append(pd.read_csv(lp, sep="\t"))
    return pd.concat(units, ignore_index=True), pd.concat(ligs, ignore_index=True)


def convert_205335_rds():
    rds = os.path.join(EXT, "tnk_GSE205335.rds")
    npz = os.path.join(EXT, "tnk_GSE205335.npz")
    if os.path.exists(npz) or not os.path.exists(rds):
        return
    subprocess.check_call([
        "Rscript", "-e",
        r'''
        obj <- readRDS("%s")
        genes <- obj$genes
        units <- obj$units
        # mean/pct are matrices with dimnames
        write.table(data.frame(gene=genes, obj$mean, check.names=FALSE),
                    "%s", sep="\t", quote=FALSE, row.names=FALSE)
        write.table(data.frame(gene=genes, obj$pct, check.names=FALSE),
                    "%s", sep="\t", quote=FALSE, row.names=FALSE)
        ''' % (
            rds,
            os.path.join(EXT, "tnk_mean_GSE205335.tsv"),
            os.path.join(EXT, "tnk_pct_GSE205335.tsv"),
        )
    ])
    mean = pd.read_csv(os.path.join(EXT, "tnk_mean_GSE205335.tsv"), sep="\t")
    pct = pd.read_csv(os.path.join(EXT, "tnk_pct_GSE205335.tsv"), sep="\t")
    genes = mean["gene"].astype(str).tolist()
    units = [c for c in mean.columns if c != "gene"]
    np.savez_compressed(
        npz,
        genes=np.array(genes, dtype=object),
        units=np.array(units, dtype=object),
        mean=mean[units].to_numpy(dtype=np.float32),
        pct=pct[units].to_numpy(dtype=np.float32),
    )


def load_tnk(cohort: str):
    npz = os.path.join(EXT, f"tnk_{cohort}.npz")
    z = np.load(npz, allow_pickle=True)
    genes = [str(g) for g in z["genes"].tolist()]
    units = [str(u) for u in z["units"].tolist()]
    return genes, units, z["mean"], z["pct"]


def family_of(lig: str) -> str:
    if lig in CORE_BARRIER:
        return "core_barrier"
    if lig in EXT_BARRIER:
        return "ext_barrier"
    if lig in CORE_IFN:
        return "core_ifn_recruit"
    if lig in EXT_IFN:
        return "ext_ifn_recruit"
    if lig in MHC:
        return "mhci"
    return "other"


def patient_family_score(lig: pd.DataFrame, elig: pd.DataFrame, ligands: list[str], col: str) -> pd.Series:
    sub = lig[lig.ligand.isin(ligands)]
    m = sub.pivot_table(index=["cohort", "unit_id"], columns="ligand", values=col, aggfunc="first")
    # align to eligible units
    idx = pd.MultiIndex.from_frame(elig[["cohort", "unit_id"]])
    m = m.reindex(idx)
    return m.mean(axis=1, skipna=True)


def ligand_sender_table(lig: pd.DataFrame, elig: pd.DataFrame, col: str) -> pd.DataFrame:
    elig_key = set(zip(elig.cohort, elig.unit_id))
    d = lig[lig.apply(lambda r: (r.cohort, r.unit_id) in elig_key, axis=1)].copy()
    rows = []
    for ligand, g in d.groupby("ligand"):
        w = wilcox_greater_zero(g[col].to_numpy())
        cohort_mu, cohort_n, cohort_se = [], [], []
        signs = []
        per = {}
        for c, gc in g.groupby("cohort"):
            x = gc[col].to_numpy(dtype=float)
            x = x[np.isfinite(x)]
            per[c] = dict(n=int(x.size), mean=float(np.mean(x)) if x.size else np.nan)
            if x.size >= 3 and np.std(x, ddof=1) > 0:
                cohort_mu.append(float(np.mean(x)))
                cohort_n.append(int(x.size))
                cohort_se.append(float(np.std(x, ddof=1) / np.sqrt(x.size)))
                signs.append(np.sign(np.mean(x)))
        meta = dl_meta(cohort_mu, cohort_se)
        pos = int(np.sum(np.array(signs) > 0))
        neg = int(np.sum(np.array(signs) < 0))
        rows.append({
            "ligand": ligand,
            "family": family_of(ligand),
            "split": col,
            "n_patients": w["n"],
            "n_nonzero": w["n_nonzero"],
            "mean_delta": w["mean"],
            "median_delta": w["median"],
            "p_wilcoxon": w["p"],
            "meta_mu": meta["mu"],
            "meta_p": meta["p"],
            "meta_I2": meta["I2"],
            "meta_k": meta["k"],
            "cohorts_pos": pos,
            "cohorts_neg": neg,
            **{f"mean_{c}": per.get(c, {}).get("mean", np.nan) for c in COHORTS},
            **{f"n_{c}": per.get(c, {}).get("n", 0) for c in COHORTS},
        })
    return pd.DataFrame(rows)


def family_sender_row(lig, elig, ligands, col, name, expect):
    s = patient_family_score(lig, elig, ligands, col)
    # s is indexed by cohort, unit
    vals = s.to_numpy(dtype=float)
    w = wilcox_greater_zero(vals)
    # per cohort
    recs = []
    mus, ses = [], []
    for c in COHORTS:
        try:
            x = s.xs(c, level=0).to_numpy(dtype=float)
        except KeyError:
            x = np.array([])
        x = x[np.isfinite(x)]
        recs.append((c, x.size, float(np.mean(x)) if x.size else np.nan))
        if x.size >= 3 and np.std(x, ddof=1) > 0:
            mus.append(float(np.mean(x)))
            ses.append(float(np.std(x, ddof=1) / np.sqrt(x.size)))
    meta = dl_meta(mus, ses)
    # leave one cohort out
    loco = {}
    for drop in COHORTS:
        keep = s[s.index.get_level_values(0) != drop]
        ww = wilcox_greater_zero(keep.to_numpy(dtype=float))
        loco[drop] = (ww["mean"], ww["p"], ww["n"])
    return {
        "family": name,
        "split": col,
        "expect": expect,
        "n_ligands_requested": len(ligands),
        "n_patients": w["n"],
        "mean_delta": w["mean"],
        "median_delta": w["median"],
        "p_wilcoxon": w["p"],
        "meta_mu": meta["mu"],
        "meta_p": meta["p"],
        "meta_I2": meta["I2"],
        "meta_k": meta["k"],
        "agrees_high_gt_low": bool(np.isfinite(w["mean"]) and w["mean"] > 0),
        "per_cohort": recs,
        "loco": loco,
    }


def build_gene_meta(elig: pd.DataFrame, cohorts=None):
    """Spearman of T/NK mean log vs malignant CLDN4 %pos, per cohort, then DL."""
    frames = []
    use = COHORTS if cohorts is None else list(cohorts)
    for c in use:
        sub = elig[elig.cohort == c]
        if sub.shape[0] < 6:
            continue
        genes, units, mean, pct = load_tnk(c)
        uix = {u: i for i, u in enumerate(units)}
        cols = [uix[u] for u in sub.unit_id if u in uix]
        if len(cols) < 6:
            continue
        y = sub.set_index("unit_id").loc[[units[i] for i in cols], "cldn4_pct"].to_numpy(dtype=float)
        M = np.array(mean[:, cols], dtype=float)
        P = pct[:, cols]
        # Spearman via Pearson of average ranks. Eligible units have finite T/NK means.
        R = rankdata(M, axis=1, method="average")
        R = R - R.mean(axis=1, keepdims=True)
        yr = rankdata(y, method="average")
        yr = yr - yr.mean()
        den = np.sqrt((R ** 2).sum(axis=1) * np.sum(yr ** 2))
        with np.errstate(invalid="ignore", divide="ignore"):
            rhos = (R @ yr) / den
        rhos[~np.isfinite(rhos)] = np.nan
        frames.append(pd.DataFrame({
            "gene": genes, "cohort": c, "rho": rhos, "n": int(len(cols)),
            "pct_median": np.nanmedian(P, axis=1),
        }))
    long = pd.concat(frames, ignore_index=True)
    rows = []
    for gene, g in long.groupby("gene"):
        meta = fisher_z_meta(g.rho.to_numpy(), g.n.to_numpy())
        signs = np.sign(g.rho.to_numpy())
        signs = signs[np.isfinite(g.rho.to_numpy())]
        rows.append({
            "gene": gene,
            "rho": meta["rho"],
            "z": meta["mu"],
            "p": meta["p"],
            "I2": meta["I2"],
            "k": meta["k"],
            "n_pos": int(np.sum(signs > 0)),
            "n_neg": int(np.sum(signs < 0)),
            "epithelial_leak": gene in EPI_LEAK,
        })
    out = pd.DataFrame(rows)
    out["fdr"] = bh(out["p"].to_numpy())
    return out, long


def expressed_fraction(elig, ligand_or_pct_long_builder):
    pass


def background_and_potential(elig, lig, gene_meta_long_pct):
    """gene_meta_long_pct is not used; compute from npz."""
    # T/NK detection: fraction of eligible patients with pct>=0.10, denominator = patients in cohorts where the gene is present
    num = {}
    den = {}
    for c in COHORTS:
        sub = elig[elig.cohort == c]
        genes, units, mean, pct = load_tnk(c)
        uix = {u: i for i, u in enumerate(units)}
        cols = [uix[u] for u in sub.unit_id if u in uix]
        if not cols:
            continue
        P = pct[:, cols]
        for i, g in enumerate(genes):
            den[g] = den.get(g, 0) + len(cols)
            num[g] = num.get(g, 0) + int(np.sum(P[i] >= 0.10))
    bg_rows = []
    for g, d in den.items():
        if d == 0:
            continue
        bg_rows.append({"gene": g, "n_patients_measured": d, "n_detected": num[g], "frac": num[g] / d})
    bg = pd.DataFrame(bg_rows)
    background = set(bg.loc[bg.frac >= 0.10, "gene"])

    # ligands: pct_mal >= 0.10 in >= 10% of eligible patients that have a row
    elig_key = set(zip(elig.cohort, elig.unit_id))
    d = lig[lig.apply(lambda r: (r.cohort, str(r.unit_id)) in elig_key or (r.cohort, r.unit_id) in elig_key, axis=1)]
    # unit ids may be int-like
    lig2 = lig.copy()
    lig2["unit_id"] = lig2["unit_id"].astype(str)
    elig2 = elig.copy()
    elig2["unit_id"] = elig2["unit_id"].astype(str)
    elig_key = set(zip(elig2.cohort, elig2.unit_id))
    d = lig2[lig2.apply(lambda r: (r.cohort, r.unit_id) in elig_key, axis=1)]
    recs = []
    lr = pd.read_csv(os.path.join(HERE, "data", "lr_network_human_21122021.tsv"), sep="\t")
    receptors = lr.groupby("from")["to"].apply(lambda s: sorted(set(s))).to_dict()
    for ligand, g in d.groupby("ligand"):
        frac = float(np.mean(g.pct_mal.to_numpy(dtype=float) >= 0.10))
        n = int(g.shape[0])
        recs_ok = receptors.get(ligand, [])
        rec_frac = 0.0
        best_rec = ""
        for rec in recs_ok:
            if rec not in den:
                continue
            fr = num[rec] / den[rec]
            if fr > rec_frac:
                rec_frac = fr
                best_rec = rec
        expressed_sender = frac >= 0.10 and n >= 5
        expressed_rec = rec_frac >= 0.10
        recs.append({
            "ligand": ligand,
            "family": family_of(ligand),
            "n_patients": n,
            "frac_patients_pct10": frac,
            "best_receptor": best_rec,
            "best_receptor_frac": rec_frac,
            "sender_expressed": expressed_sender,
            "receptor_expressed": expressed_rec,
            "potential": bool(expressed_sender and expressed_rec),
        })
    pot = pd.DataFrame(recs)
    # potential requires the ligand column to exist in the prior; activity.R intersects
    return bg, background, pot


def _sign_ok(g: pd.DataFrame, sign: str) -> pd.Series:
    """Positive (or negative) in every contributing cohort, and in at least 2."""
    if sign == "high":
        return (g.k >= 2) & (g.n_pos == g.k) & (g.rho > 0)
    return (g.k >= 2) & (g.n_neg == g.k) & (g.rho < 0)


def choose_geneset(gene_meta: pd.DataFrame, sign: str) -> tuple[list[str], str]:
    g = gene_meta[~gene_meta.epithelial_leak].copy()
    g = g[_sign_ok(g, sign)]
    direction = "up" if sign == "high" else "down"
    strict = g[g.fdr <= 0.10]
    rule = (
        f"FDR<=0.10, meta rho {direction} with CLDN4, same sign in every cohort that contributed "
        "(k>=2), epithelial genes removed"
    )
    if len(strict) >= 15:
        return strict.sort_values("p").gene.tolist(), rule
    relax = g[g.p <= 0.01]
    rule = (
        "FALLBACK unadjusted p<=0.01, same sign in every contributing cohort, epithelial genes removed "
        "(FDR set had <15 genes)"
    )
    if len(relax) >= 15:
        return relax.sort_values("p").gene.tolist(), rule
    rank = g.sort_values("z", ascending=(sign != "high")).head(100)
    rule = (
        "FALLBACK top 100 |meta z| with the same sign in every contributing cohort "
        "(p<=0.05 set had <15 genes)"
    )
    return rank.gene.tolist(), rule


def exact_label_p(scores: dict[str, float], a: list[str], b: list[str]) -> dict:
    aa = [lig for lig in a if lig in scores and np.isfinite(scores[lig])]
    bb = [lig for lig in b if lig in scores and np.isfinite(scores[lig])]
    if len(aa) < 2 or len(bb) < 2:
        return dict(n_a=len(aa), n_b=len(bb), diff=np.nan, p_two=np.nan, mean_a=np.nan, mean_b=np.nan)
    obs_a = np.mean([scores[x] for x in aa])
    obs_b = np.mean([scores[x] for x in bb])
    diff = obs_a - obs_b
    pool = aa + bb
    n_a = len(aa)
    # exact reassignment of the pool into a group of size n_a
    n_ge = 0
    n_tot = 0
    vals = np.array([scores[x] for x in pool], dtype=float)
    for comb in combinations(range(len(pool)), n_a):
        sel = np.zeros(len(pool), dtype=bool)
        sel[list(comb)] = True
        d = vals[sel].mean() - vals[~sel].mean()
        n_tot += 1
        if abs(d) >= abs(diff) - 1e-15:
            n_ge += 1
    return dict(n_a=len(aa), n_b=len(bb), diff=float(diff), mean_a=float(obs_a), mean_b=float(obs_b),
                p_two=n_ge / n_tot, n_perm=n_tot)


def random_set_p(scores: dict[str, float], ligands: list[str], n_draw=4000, seed=1) -> dict:
    present = [lig for lig in ligands if lig in scores and np.isfinite(scores[lig])]
    pool = np.array([v for k, v in scores.items() if np.isfinite(v)])
    if len(present) < 2 or pool.size < len(present) + 2:
        return dict(mean=np.nan, p_greater=np.nan, n=len(present))
    obs = float(np.mean([scores[x] for x in present]))
    rng = np.random.default_rng(seed)
    k = len(present)
    draws = np.array([rng.choice(pool, size=k, replace=False).mean() for _ in range(n_draw)])
    # one-sided: family mean greater than random ligands
    p = (np.sum(draws >= obs) + 1) / (n_draw + 1)
    return dict(mean=obs, p_greater=float(p), n=k)


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_n(x, d=3, signed=False):
    if x is None or not np.isfinite(x):
        return "NA"
    if signed:
        return f"{x:+.{d}f}"
    return f"{x:.{d}f}"


def main():
    convert_205335_rds()
    units, lig = load_units_ligands()
    units["unit_id"] = units["unit_id"].astype(str)
    lig["unit_id"] = lig["unit_id"].astype(str)
    # Alias collisions (PVRL2 and NECTIN2) can emit two rows. Keep the row with detection.
    lig = lig.sort_values("pct_mal", ascending=False).drop_duplicates(
        ["cohort", "unit_id", "ligand"], keep="first"
    )
    for c in units.cohort.unique():
        if c not in ALLOWED:
            raise SystemExit(c)
    units.to_csv(os.path.join(TAB, "patient_inventory.tsv"), sep="\t", index=False)

    # Quartile arms are a CLDN4 contrast only when at least a quarter of malignant
    # cells express CLDN4. Below that, ties at zero fill the "high" arm.
    # Those units stay in the percent-positive split.
    base = units[units.eligible_q4.astype(bool) & units.cldn4_separated.astype(bool)].copy()
    elig = base[base.cldn4_pct >= 0.25].copy()
    elig_pct = base[(base.cldn4_pct > 0) & (base.cldn4_pct < 1)].copy()
    elig.to_csv(os.path.join(TAB, "eligible_units.tsv"), sep="\t", index=False)

    sender = ligand_sender_table(lig, elig, "delta_q4q1")
    sender_med = ligand_sender_table(lig, elig, "delta_median")
    sender_pct = ligand_sender_table(lig, elig_pct, "delta_pctpos")
    sender_all = pd.concat([sender, sender_med, sender_pct], ignore_index=True)
    sender_all.to_csv(os.path.join(TAB, "ligand_sender_delta.tsv"), sep="\t", index=False)

    families = [
        ("core_barrier", CORE_BARRIER, "high>low"),
        ("ext_barrier", EXT_BARRIER, "high>low"),
        ("core_ifn_recruit", CORE_IFN, "low>high"),
        ("ext_ifn_recruit", EXT_IFN, "low>high"),
        ("mhci", MHC, "low>high"),
    ]
    fam_rows = []
    for col, use in (
        ("delta_q4q1", elig),
        ("delta_median", elig),
        ("delta_pctpos", elig_pct),
    ):
        for name, genes, expect in families:
            fam_rows.append(family_sender_row(lig, use, genes, col, name, expect))
    # serialize loco
    fam_flat = []
    for r in fam_rows:
        base = {k: v for k, v in r.items() if k not in ("per_cohort", "loco")}
        for c, n, mu in r["per_cohort"]:
            base[f"n_{c}"] = n
            base[f"mean_{c}"] = mu
        for c, (mu, p, n) in r["loco"].items():
            base[f"loco_mean_drop_{c}"] = mu
            base[f"loco_p_drop_{c}"] = p
            base[f"loco_n_drop_{c}"] = n
        fam_flat.append(base)
    fam_df = pd.DataFrame(fam_flat)
    fam_df.to_csv(os.path.join(TAB, "family_sender.tsv"), sep="\t", index=False)

    print("building receiver gene meta", flush=True)
    # Activity geneset uses author-annotated T/NK only (GSE131907, GSE205335).
    # Marker-gated cohorts are easier to contaminate with epithelial doublets.
    gene_meta, gene_long = build_gene_meta(elig, AUTHOR_COHORTS)
    gene_meta.to_csv(os.path.join(TAB, "tnk_gene_vs_cldn4.tsv"), sep="\t", index=False)
    gene_meta_all, _ = build_gene_meta(elig)
    gene_meta_all.to_csv(os.path.join(TAB, "tnk_gene_vs_cldn4_all4.tsv"), sep="\t", index=False)
    # gene_long can be large; keep rho only
    gene_long.to_csv(os.path.join(TAB, "tnk_gene_vs_cldn4_by_cohort.tsv.gz"), sep="\t", index=False)

    gs_high, rule_high = choose_geneset(gene_meta, "high")
    gs_low, rule_low = choose_geneset(gene_meta, "low")

    print("background / potential", flush=True)
    bg, background, pot = background_and_potential(elig, lig, None)
    bg.to_csv(os.path.join(TAB, "background_genes.tsv"), sep="\t", index=False)
    pot.to_csv(os.path.join(TAB, "potential_ligands.tsv"), sep="\t", index=False)
    potential = pot.loc[pot.potential, "ligand"].tolist()

    # genesets intersect background
    def _clip(gs):
        return [g for g in gs if g in background]

    gs_high_b = _clip(gs_high)
    gs_low_b = _clip(gs_low)
    programs = {
        "empirical_high": gs_high_b,
        "empirical_low": gs_low_b,
        "a_priori_ifn": _clip(IFN_PROGRAM),
        "a_priori_cytotoxicity": _clip(CYTO_PROGRAM),
        "a_priori_exhaustion": _clip(EXHAUST_PROGRAM),
    }

    act_in = "/tmp/nichenet_work/activity_in"
    act_out = "/tmp/nichenet_work/activity_out"
    os.makedirs(act_in, exist_ok=True)
    os.makedirs(act_out, exist_ok=True)
    with open(os.path.join(act_in, "background_genes.txt"), "w") as f:
        f.write("\n".join(sorted(background)) + "\n")
    with open(os.path.join(act_in, "potential_ligands.txt"), "w") as f:
        f.write("\n".join(potential) + "\n")
    for name, genes in programs.items():
        with open(os.path.join(act_in, f"geneset_{name}.txt"), "w") as f:
            f.write("\n".join(genes) + "\n")
    # continuous z
    cont = gene_meta[(~gene_meta.epithelial_leak) & gene_meta.gene.isin(background) & np.isfinite(gene_meta.z)]
    cont[["gene", "z"]].to_csv(os.path.join(act_in, "continuous_response.tsv"), sep="\t", index=False)

    with open(os.path.join(TAB, "geneset_rules.json"), "w") as f:
        json.dump({
            "empirical_high_rule": rule_high,
            "empirical_low_rule": rule_low,
            "n_high_before_background": len(gs_high),
            "n_low_before_background": len(gs_low),
            "n_high": len(gs_high_b),
            "n_low": len(gs_low_b),
            "n_background": len(background),
            "n_potential": len(potential),
            "programs": {k: len(v) for k, v in programs.items()},
        }, f, indent=2)

    print("running NicheNet activity in R", flush=True)
    subprocess.check_call([
        "Rscript", os.path.join(HERE, "scripts", "activity.R"),
        f"--in={act_in}", f"--out={act_out}", f"--src={SRC}",
        f"--prior={os.path.join(PRIOR, 'ligand_target_matrix_nsga2r_final.rds')}",
    ])

    acts = []
    for name in list(programs) + ["continuous_meta_z"]:
        fp = os.path.join(act_out, f"activity_{name}.tsv")
        if os.path.exists(fp) and os.path.getsize(fp) > 20:
            a = pd.read_csv(fp, sep="\t")
            if "geneset" not in a.columns:
                a["geneset"] = name
            acts.append(a)
    activity = pd.concat(acts, ignore_index=True) if acts else pd.DataFrame()
    activity.to_csv(os.path.join(TAB, "ligand_activity_all.tsv"), sep="\t", index=False)

    # merge sender delta onto activity for prioritization on empirical_high
    send_q = sender.set_index("ligand")
    pri_rows = []
    if not activity.empty and "aupr_corrected" in activity.columns:
        for gs, g in activity.groupby("geneset"):
            if gs == "continuous_meta_z" or "aupr_corrected" not in g.columns:
                continue
            gg = g.copy()
            gg["mean_delta"] = gg["test_ligand"].map(send_q["mean_delta"])
            gg["p_sender"] = gg["test_ligand"].map(send_q["p_wilcoxon"])
            gg["family"] = gg["test_ligand"].map(family_of)
            rec_frac = pot.set_index("ligand")["best_receptor_frac"]
            gg["receptor_frac"] = gg["test_ligand"].map(rec_frac)
            for col, zname in (
                ("aupr_corrected", "z_activity"),
                ("mean_delta", "z_sender"),
                ("receptor_frac", "z_receptor"),
            ):
                x = gg[col].to_numpy(dtype=float)
                mu, sd = np.nanmean(x), np.nanstd(x, ddof=1)
                gg[zname] = (x - mu) / sd if sd and np.isfinite(sd) and sd > 0 else np.nan
            gg["priority"] = gg[["z_activity", "z_sender", "z_receptor"]].sum(axis=1, min_count=1)
            gg["geneset"] = gs
            pri_rows.append(gg)
    priority = pd.concat(pri_rows, ignore_index=True) if pri_rows else pd.DataFrame()
    if not priority.empty:
        priority.to_csv(os.path.join(TAB, "ligand_priority.tsv"), sep="\t", index=False)

    # contrasts
    contrast_rows = []
    if not priority.empty:
        for gs, g in priority.groupby("geneset"):
            score_act = dict(zip(g.test_ligand, g.aupr_corrected))
            score_pri = dict(zip(g.test_ligand, g.priority))
            for metric, scores in (("aupr_corrected", score_act), ("priority", score_pri)):
                lab = exact_label_p(scores, CORE_BARRIER, CORE_IFN)
                contrast_rows.append({"geneset": gs, "metric": metric, "contrast": "core_barrier_minus_core_ifn", **lab})
                for fam_name, fam in (
                    ("core_barrier", CORE_BARRIER),
                    ("core_ifn_recruit", CORE_IFN),
                    ("ext_barrier", EXT_BARRIER),
                    ("mhci", MHC),
                ):
                    rnd = random_set_p(scores, fam, seed=abs(hash(gs + metric + fam_name)) % 10_000)
                    contrast_rows.append({
                        "geneset": gs, "metric": metric, "contrast": f"{fam_name}_vs_random",
                        "n_a": rnd["n"], "n_b": np.nan, "diff": np.nan,
                        "mean_a": rnd["mean"], "mean_b": np.nan, "p_two": rnd["p_greater"], "n_perm": 4000,
                    })
    contrast = pd.DataFrame(contrast_rows)
    contrast.to_csv(os.path.join(TAB, "family_activity_contrast.tsv"), sep="\t", index=False)

    # called-out table
    called = CORE_BARRIER + CORE_IFN + MHC + ["CXCL11", "PVR", "CD274", "TGFB1", "HLA-E"]
    call_rows = []
    for lig_name in called:
        row = {"ligand": lig_name, "family": family_of(lig_name)}
        if lig_name in send_q.index:
            row["mean_delta_q4q1"] = send_q.loc[lig_name, "mean_delta"]
            row["p_sender"] = send_q.loc[lig_name, "p_wilcoxon"]
            row["n_patients"] = send_q.loc[lig_name, "n_patients"]
            row["meta_I2"] = send_q.loc[lig_name, "meta_I2"]
            for c in COHORTS:
                row[f"mean_{c}"] = send_q.loc[lig_name, f"mean_{c}"]
        pot_hit = pot[pot.ligand == lig_name]
        row["potential"] = bool(pot_hit.potential.iloc[0]) if len(pot_hit) else False
        if not activity.empty:
            for gs in ("empirical_high", "empirical_low", "a_priori_ifn", "a_priori_exhaustion", "a_priori_cytotoxicity"):
                sub = activity[(activity.geneset == gs) & (activity.test_ligand == lig_name)]
                if len(sub) and "aupr_corrected" in sub.columns:
                    row[f"aupr_{gs}"] = float(sub.aupr_corrected.iloc[0])
                    row[f"pearson_{gs}"] = float(sub.pearson.iloc[0])
                    row[f"rank_aupr_{gs}"] = int(sub.rank_aupr.iloc[0])
                    row[f"auroc_{gs}"] = float(sub.auroc.iloc[0])
        if not priority.empty:
            sub = priority[(priority.geneset == "empirical_high") & (priority.test_ligand == lig_name)]
            if len(sub):
                row["priority_empirical_high"] = float(sub.priority.iloc[0])
        call_rows.append(row)
    called_df = pd.DataFrame(call_rows)
    called_df.to_csv(os.path.join(TAB, "called_ligands.tsv"), sep="\t", index=False)

    # top ligands by priority and by activity
    if not priority.empty and "empirical_high" in set(priority.geneset):
        top = priority[priority.geneset == "empirical_high"].sort_values("aupr_corrected", ascending=False).head(25)
        top.to_csv(os.path.join(TAB, "top_activity_empirical_high.tsv"), sep="\t", index=False)
        priority[priority.geneset == "empirical_low"].sort_values("aupr_corrected", ascending=False).head(25).to_csv(
            os.path.join(TAB, "top_activity_empirical_low.tsv"), sep="\t", index=False
        )

    write_results(units, elig, fam_df, sender, called_df, contrast, gene_meta, programs, rule_high, rule_low, pot, activity)
    make_figures(fam_df, called_df, sender, priority)
    print("analyze done", flush=True)


def write_results(units, elig, fam_df, sender, called, contrast, gene_meta, programs, rule_high, rule_low, pot, activity):
    def fam(name, split="delta_q4q1"):
        hit = fam_df[(fam_df.family == name) & (fam_df.split == split)]
        return hit.iloc[0] if len(hit) else None

    lines = []
    a = lines.append
    a("# RESULTS — NicheNet ligand activity, concordant-4")
    a("")
    a("CLDN4-high vs CLDN4-low **malignant senders** to **T/NK receivers**.")
    a("Cohorts are only the locked concordant four: GSE123902, GSE131907, GSE205335, GSE189357.")
    a("GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, and GSE207422 are not in this run.")
    a("TACSTD2 is not a gate. There is no dual-high call.")
    a("")
    a("Activity is `predict_ligand_activities` from saeyslab/nichenetr")
    a("`66f90d5eeafef280b2b2f339b3fd70ffec1781dd` (2026-06-01), executed in R:")
    a("`get_single_ligand_importances` → `evaluate_target_prediction` →")
    a("`classification_evaluation_continuous_pred` (ROCR curves, `caTools::trapz` AUPR).")
    a("The prior is NicheNet v2 `ligand_target_matrix_nsga2r_final.rds` (Zenodo 7074291).")
    a("Official rank is **AUPR corrected** (AUPR minus the positive-class prevalence). Pearson and AUROC are reported beside it.")
    a("")
    a("## Honest n")
    a("")
    a("| cohort | units loaded | eligible quartile (n_mal≥40, n_tnk≥20, CLDN4 high>low, ≥25% CLDN4+) |")
    a("|---|---:|---:|")
    for c in COHORTS:
        n_all = int((units.cohort == c).sum())
        n_el = int((elig.cohort == c).sum())
        a(f"| {c} | {n_all} | {n_el} |")
    a(f"| **all** | **{len(units)}** | **{len(elig)}** |")
    a("")
    a("The unit is the locked patient or sample (GSE131907 stays at the sample id used in the concordant-4 inventory, not a further patient collapse).")
    a("Cells are not the independent unit. Primary sender split is **within-unit**: malignant cells ranked by CLDN4 `log1p(CP10k)`, top quartile vs bottom quartile (`rank` ties broken by first occurrence, same rule as the CellChat concordant-4 run).")
    a("A unit enters the quartile test only when CLDN4 itself is higher in the top quartile than the bottom quartile, and at least 25% of malignant cells express CLDN4. Below 25% positive, zero-ties would fill the high arm. Those units are kept for the percent-positive split.")
    a("")
    a("## Sender expression — within-patient Q4 minus Q1")
    a("")
    a("For each eligible unit, the family score is the mean ligand Δ across the pre-specified ligands measured in that unit.")
    a("The test is a Wilcoxon signed-rank test of that patient-level score against 0.")
    a("Cohort means are also combined with DerSimonian–Laird random effects. I² is that meta-analysis, not a cell-pooled test.")
    a("")
    a("| family | expect | n patients | mean Δ | Wilcoxon p | DL mean | DL p | I² | agrees with expect |")
    a("|---|---|---:|---:|---:|---:|---:|---:|---|")
    for name, expect in (
        ("core_barrier", "high>low"),
        ("ext_barrier", "high>low"),
        ("core_ifn_recruit", "low>high"),
        ("ext_ifn_recruit", "low>high"),
        ("mhci", "low>high"),
    ):
        r = fam(name)
        if r is None:
            continue
        agrees = (r.mean_delta > 0) if expect == "high>low" else (r.mean_delta < 0)
        a(
            f"| {name} | {expect} | {int(r.n_patients)} | {fmt_n(r.mean_delta, signed=True)} | {fmt_p(r.p_wilcoxon)} | "
            f"{fmt_n(r.meta_mu, signed=True)} | {fmt_p(r.meta_p)} | {fmt_n(100 * r.meta_I2, 0)}% | {'yes' if agrees else 'no'} |"
        )
    a("")
    a("Core barrier ligands are F11R, NECTIN2, CDH1, LGALS9.")
    a("Core IFN/recruit ligands are CXCL9, CXCL10, CCL5, IFNG.")
    a("MHC-I (HLA-A/B/C) is its own row because the CellChat concordant-4 run found HLA–CD8 communication higher, not lower, from CLDN4-high senders. That opposite direction is what this quartile split shows as well, except GSE189357, where the MHC-I family mean is negative.")
    a("")
    a("Core barrier stays positive in every leave-one-cohort-out (see `family_sender.tsv`). The core IFN/recruit family mean is negative in GSE123902, GSE205335, and GSE189357 and positive in GSE131907, so I² is high and the random-effects p does not match the pooled Wilcoxon.")
    a("CXCL9, CXCL10, and IFNG are detected in ≥10% of malignant cells in fewer than 10% of eligible units, so they are not potential ligands. Their within-patient deltas are near zero. CCL5 is the only core IFN/recruit ligand that passes the expression filter.")
    a("")
    a("### Called ligands")
    a("")
    a("| ligand | family | n | mean Δ Q4−Q1 | p | potential ligand | AUPR-corrected rank on T/NK genes up with CLDN4 | AUPR-corrected rank on T/NK genes down with CLDN4 |")
    a("|---|---|---:|---:|---:|---|---:|---:|")
    if called is not None and len(called):
        for rec in called.itertuples(index=False):
            a(
                f"| {rec.ligand} | {rec.family} | {getattr(rec, 'n_patients', 'NA')} | "
                f"{fmt_n(getattr(rec, 'mean_delta_q4q1', np.nan), signed=True)} | {fmt_p(getattr(rec, 'p_sender', np.nan))} | "
                f"{'yes' if getattr(rec, 'potential', False) else 'no'} | "
                f"{getattr(rec, 'rank_aupr_empirical_high', 'NA')} | {getattr(rec, 'rank_aupr_empirical_low', 'NA')} |"
            )
    a("")
    a("Full per-ligand sender table: `results/tables/ligand_sender_delta.tsv`.")
    a("")
    a("## Receiver geneset (patient-aware)")
    a("")
    a("The empirical geneset is fit on author-annotated T/NK only (GSE131907 and GSE205335). Marker-gated T/NK calls in GSE123902 and GSE189357 are not used to choose the geneset.")
    a("Within each author cohort, Spearman correlation of the T/NK pseudobulk (`mean log1p(CP10k)`) with malignant CLDN4 % positive. The two correlations are Fisher-z combined (DerSimonian–Laird). A gene must keep the same sign in both cohorts.")
    a(f"Genes up with CLDN4: **{len(programs['empirical_high'])}** genes. Rule: {rule_high}.")
    a(f"Genes down with CLDN4: **{len(programs['empirical_low'])}** genes. Rule: {rule_low}.")
    a("Lung epithelial, secretory, and ciliated markers are removed before the cutoff (EPCAM, TACSTD2, CDH1, KRTs, CLDNs, surfactant and secretoglobin genes, ELF3, CEACAM5/6, WFDC2, FOXJ1, TMC5, and the rest of the list in METHODS).")
    a("The all-four correlation table is `tnk_gene_vs_cldn4_all4.tsv`. It is not the activity geneset.")
    # top genes
    g = gene_meta[(~gene_meta.epithelial_leak) & gene_meta.gene.isin(programs["empirical_high"])].sort_values("p").head(12)
    a("")
    a("Highest-confidence T/NK genes **up** with malignant CLDN4 %pos (meta ρ, p):")
    a("")
    if len(g):
        a(", ".join(f"{r.gene} (ρ={r.rho:+.2f}, p={fmt_p(r.p)})" for r in g.itertuples(index=False)))
    g2 = gene_meta[(~gene_meta.epithelial_leak) & gene_meta.gene.isin(programs["empirical_low"])].sort_values("p").head(12)
    a("")
    a("Highest-confidence T/NK genes **down** with malignant CLDN4 %pos:")
    a("")
    if len(g2):
        a(", ".join(f"{r.gene} (ρ={r.rho:+.2f}, p={fmt_p(r.p)})" for r in g2.itertuples(index=False)))
    a("")
    a(f"Background expressed in T/NK (≥10% of eligible units with detection in ≥10% of T/NK cells, gene measured in that unit's cohort): see `background_genes.tsv`.")
    a(f"Potential ligands (sender detected in ≥10% of malignant cells in ≥10% of eligible units, and ≥1 NicheNet-v2 receptor detected in T/NK at the same threshold): **{int(pot.potential.sum()) if len(pot) else 0}**.")
    a("")
    a("## NicheNet activity")
    a("")
    a("Activity asks which malignant-expressed ligands' NicheNet-v2 target profiles match the T/NK geneset.")
    a("It is not a second copy of the sender Δ. A ligand can be higher in CLDN4-high cells and still be a weak predictor of the T/NK program, or the reverse.")
    a("")
    a("Core IFN/recruit contributes only CCL5 as a potential ligand, so a barrier-versus-IFN activity contrast inside that pair of families is not estimated.")
    a("The comparison that is estimated is whether the four core barrier ligands, as a set, have higher AUPR-corrected than a random draw of four potential ligands (one-sided, 4000 draws).")
    a("")
    a("| geneset | n genes | core barrier mean AUPR-corrected | one-sided p vs random ligands |")
    a("|---|---:|---:|---:|")
    if contrast is not None and len(contrast):
        sub = contrast[(contrast.metric == "aupr_corrected") & (contrast.contrast == "core_barrier_vs_random")]
        for rec in sub.itertuples(index=False):
            n_g = len(programs.get(rec.geneset, []))
            a(f"| {rec.geneset} | {n_g} | {fmt_n(rec.mean_a)} | {fmt_p(rec.p_two)} |")
    a("")
    a("Ranks below are among the potential ligands (lower rank = higher AUPR-corrected). Absolute AUPR-corrected values on these genesets are small; the rank is the comparison.")
    a("")
    a("Prioritization on the empirical-high geneset, in the differential-NicheNet style, is the sum of z-scored AUPR-corrected, z-scored within-patient sender Δ (high − low), and z-scored best-receptor detection. Stored in `ligand_priority.tsv`.")
    a("AUPR-corrected remains the activity rank. The sum is a joint score, not a substitute for either arm.")
    a("")
    if not activity.empty and "empirical_high" in set(activity.geneset.astype(str)):
        top = activity[activity.geneset == "empirical_high"].sort_values("aupr_corrected", ascending=False).head(8)
        a("Top potential ligands by AUPR-corrected on T/NK genes **up** with CLDN4:")
        a("")
        a("| rank | ligand | AUPR-corrected | AUROC | Pearson |")
        a("|---:|---|---:|---:|---:|")
        for rec in top.itertuples(index=False):
            a(f"| {int(rec.rank_aupr)} | {rec.test_ligand} | {fmt_n(rec.aupr_corrected)} | {fmt_n(rec.auroc)} | {fmt_n(rec.pearson, signed=True)} |")
        a("")
        top = activity[activity.geneset == "empirical_low"].sort_values("aupr_corrected", ascending=False).head(8)
        if len(top):
            a("Top potential ligands by AUPR-corrected on T/NK genes **down** with CLDN4:")
            a("")
            a("| rank | ligand | AUPR-corrected | AUROC | Pearson |")
            a("|---:|---|---:|---:|---:|")
            for rec in top.itertuples(index=False):
                a(f"| {int(rec.rank_aupr)} | {rec.test_ligand} | {fmt_n(rec.aupr_corrected)} | {fmt_n(rec.auroc)} | {fmt_n(rec.pearson, signed=True)} |")
            a("")
    a("A priori receiver programs (IFN, cytotoxicity, exhaustion) are scored with the same function and the same background. They are not the CLDN4 contrast. Ranks for the called ligands are in `called_ligands.tsv`.")
    a("")
    a("## Other patient-aware splits")
    a("")
    a("Same eligible units. Median split: malignant CLDN4 above vs at-or-below the unit median. Percent-positive split: CLDN4 UMI > 0 vs 0. Family Wilcoxon:")
    a("")
    a("| family | median-split mean Δ | p | %pos-split mean Δ | p |")
    a("|---|---:|---:|---:|---:|")
    for name in ("core_barrier", "core_ifn_recruit", "mhci"):
        rm, rp = fam(name, "delta_median"), fam(name, "delta_pctpos")
        if rm is None:
            continue
        a(f"| {name} | {fmt_n(rm.mean_delta, signed=True)} | {fmt_p(rm.p_wilcoxon)} | {fmt_n(rp.mean_delta, signed=True)} | {fmt_p(rp.p_wilcoxon)} |")
    a("")
    a("Leave-one-cohort-out means for the primary Q4−Q1 family score are columns `loco_mean_drop_*` in `family_sender.tsv`.")
    a("A cohort that flips the family sign when it is the only one left out is visible there. Cohorts are not dropped to manufacture agreement, and no discordant accession is added.")
    a("")
    a("## What this does not claim")
    a("")
    a("- It does not re-estimate the locked concordant-4 CLDN4 %pos vs T/NK fraction correlation.")
    a("- It does not turn a NicheNet target score into a statement that a barrier ligand excludes T cells in space. That evidence is the CosMx contact result, not this table.")
    a("- It does not claim CXCL9/10 are abundant. If they fail the potential-ligand filter, their activity is not scored.")
    a("- Cell-pooled p-values are not reported.")
    a("- The ligand–target matrix is a published prior, not a model fit on these four cohorts.")
    a("")
    a("## Reproduce")
    a("")
    a("```bash")
    a("bash methods/nichenet_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw /tmp/nichenet_prior")
    a("Rscript methods/nichenet_concordant4_cldn4/scripts/extract_gse205335.R --raw=/tmp/concordant4_raw --here=methods/nichenet_concordant4_cldn4 --out=/tmp/nichenet_work/extract")
    a("python3 methods/nichenet_concordant4_cldn4/scripts/extract_rest.py")
    a("python3 methods/nichenet_concordant4_cldn4/scripts/analyze.py")
    a("```")
    a("")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(HERE, "RESULTS.md"), "w") as f:
        f.write(text)
    with open(os.path.join(HERE, "results", "RESULTS.md"), "w") as f:
        f.write(text)


def make_figures(fam_df, called, sender, priority):
    # Figure 1: family mean delta
    want = [
        ("core_barrier", "Barrier\nF11R NECTIN2\nCDH1 LGALS9", BAR_COLOR),
        ("ext_barrier", "Extended\nbarrier", "#5b8fb8"),
        ("core_ifn_recruit", "IFN / recruit\nCXCL9 CXCL10\nCCL5 IFNG", IFN_COLOR),
        ("ext_ifn_recruit", "Extended\nIFN / recruit", "#d4a017"),
        ("mhci", "MHC-I\nHLA-A/B/C", MHC_COLOR),
    ]
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    xs, means, ps, colors, labels = [], [], [], [], []
    for i, (name, lab, col) in enumerate(want):
        hit = fam_df[(fam_df.family == name) & (fam_df.split == "delta_q4q1")]
        if not len(hit):
            continue
        xs.append(i)
        means.append(hit.iloc[0].mean_delta)
        ps.append(hit.iloc[0].p_wilcoxon)
        colors.append(col)
        labels.append(lab)
    ax.axhline(0, color="#888", lw=0.8)
    ax.bar(xs, means, color=colors, width=0.72)
    for x, m, p in zip(xs, means, ps):
        ax.text(x, m + (0.01 if m >= 0 else -0.02), f"p={fmt_p(p)}", ha="center", va="bottom" if m >= 0 else "top", fontsize=8)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Mean within-patient Δ log1p(CP10k)\nCLDN4 Q4 malignant − Q1 malignant")
    ax.set_title("Sender ligands to the same patient's T/NK\npatient is the unit")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_sender_family_delta.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "fig1_sender_family_delta.pdf"))
    plt.close()

    # Figure 2: called ligand sender delta with cohort dots via sender table means
    if called is not None and len(called) and "mean_delta_q4q1" in called.columns:
        show = called[called.ligand.isin(CORE_BARRIER + CORE_IFN + MHC)].copy()
        fig, ax = plt.subplots(figsize=(9.2, 4.6))
        order = [g for g in CORE_BARRIER + CORE_IFN + MHC if g in set(show.ligand)]
        show = show.set_index("ligand").loc[order]
        cols = [BAR_COLOR] * len(CORE_BARRIER) + [IFN_COLOR] * len(CORE_IFN) + [MHC_COLOR] * len(MHC)
        cols = cols[: len(order)]
        ax.axhline(0, color="#888", lw=0.8)
        ax.bar(range(len(order)), show.mean_delta_q4q1, color=cols)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right")
        ax.set_ylabel("Mean Δ (Q4 − Q1 malignant)")
        ax.set_title("Called ligands, within-patient CLDN4 quartile")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fig2_called_ligand_delta.png"), dpi=160)
        fig.savefig(os.path.join(FIG, "fig2_called_ligand_delta.pdf"))
        plt.close()

    # Figure 3: activity vs sender delta for empirical high
    if priority is not None and len(priority) and "empirical_high" in set(priority.geneset):
        g = priority[priority.geneset == "empirical_high"].copy()
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        other = g[g.family == "other"]
        ax.scatter(other.mean_delta, other.aupr_corrected, s=12, c="#c5cdd4", label="other potential", zorder=1)
        for fam, col, lab in (
            ("core_barrier", BAR_COLOR, "core barrier"),
            ("core_ifn_recruit", IFN_COLOR, "core IFN/recruit"),
            ("mhci", MHC_COLOR, "MHC-I"),
        ):
            sub = g[g.family == fam]
            ax.scatter(sub.mean_delta, sub.aupr_corrected, s=46, c=col, label=lab, zorder=2)
            for rec in sub.itertuples(index=False):
                ax.annotate(rec.test_ligand, (rec.mean_delta, rec.aupr_corrected), fontsize=7, xytext=(4, 3), textcoords="offset points")
        ax.axvline(0, color="#888", lw=0.6)
        ax.set_xlabel("Sender Δ, CLDN4-high − CLDN4-low malignant")
        ax.set_ylabel("NicheNet AUPR-corrected\nT/NK genes up with CLDN4")
        ax.legend(frameon=False, fontsize=8)
        ax.set_title("Activity vs sender change")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fig3_activity_vs_sender.png"), dpi=160)
        fig.savefig(os.path.join(FIG, "fig3_activity_vs_sender.pdf"))
        plt.close()


if __name__ == "__main__":
    main()
