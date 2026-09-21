#!/usr/bin/env python3
"""Cross-cohort 4/4 gene lists, differential modules, and Enrichr/gseapy.

Rules are fixed here and are not retuned after looking at the lists.

Patient-effect strict (UP), all four cohorts:
  patient-mean (Q4-Q1) log1p CP10k >= 0.10
  fraction of units with delta > 0 >= 0.60
  pooled-cell AUC >= 0.55
  fraction of CLDN4-high cells with count > 0 >= 0.05
DOWN uses the opposite signs and pct_low.

COMET top decile, per cohort, among genes with the matching detection floor:
  rank by XL-mHG statistic (ascending), then by cell delta in that direction.
  Keep the top 10%. A COMET 4/4 gene is in that decile in every cohort
  and has the matching sign of the patient-mean delta in every cohort.

CLDN4 is the split, so it is reported as a control and left out of marker lists.
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
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import fisher_exact, hypergeom, spearmanr, wilcoxon
import xlmhg.mhg_cython as xlmhg_cython

if not hasattr(np, "float"):
    np.float = np.float64  # type: ignore[attr-defined]
if not hasattr(np, "bool"):
    np.bool = np.bool_  # type: ignore[attr-defined]
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
LIBS = [
    "GO_Biological_Process_2023",
    "KEGG_2021_Human",
    "MSigDB_Hallmark_2020",
    "Reactome_2022",
]
ROOT = Path(__file__).resolve().parents[1]
STRICT_DELTA = 0.10
STRICT_SIGN = 0.60
STRICT_AUC = 0.05
STRICT_PCT = 0.05
DECILE = 0.10
PAIR_P = 0.01
PAIR_FOLD = 1.5
RNG = np.random.default_rng(1)


def bh(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


def is_mt_ribo(gene: str) -> bool:
    g = gene.upper()
    return g.startswith("MT-") or g.startswith(("RPL", "RPS", "MRPL", "MRPS"))


def load_stats(table_dir: Path) -> pd.DataFrame:
    pieces = []
    for c in COHORTS:
        df = pd.read_csv(table_dir / f"{c}_gene_stats.tsv.gz", sep="\t")
        df["gene"] = df["gene"].astype(str).str.upper()
        pieces.append(df)
    wide = None
    for df in pieces:
        c = df["cohort"].iloc[0]
        keep = df.rename(columns={
            "patient_delta_mean": f"{c}__patient_delta",
            "sign_frac_up": f"{c}__sign_up",
            "sign_frac_down": f"{c}__sign_down",
            "cell_delta": f"{c}__cell_delta",
            "auc": f"{c}__auc",
            "pct_high": f"{c}__pct_high",
            "pct_low": f"{c}__pct_low",
            "xlmhg_stat_up": f"{c}__stat_up",
            "xlmhg_stat_down": f"{c}__stat_down",
            "xlmhg_cutoff_up": f"{c}__cut_up",
            "xlmhg_cutoff_down": f"{c}__cut_down",
            "n_units": f"{c}__n_units",
            "n_high": f"{c}__n_high",
            "n_low": f"{c}__n_low",
        })
        cols = ["gene"] + [x for x in keep.columns if x.startswith(c + "__")]
        part = keep[cols].drop_duplicates("gene")
        wide = part if wide is None else wide.merge(part, on="gene", how="inner")
    return wide


def decile_genes(wide: pd.DataFrame, direction: str) -> set[str]:
    sets = []
    for c in COHORTS:
        pct = f"{c}__pct_high" if direction == "up" else f"{c}__pct_low"
        stat = f"{c}__stat_up" if direction == "up" else f"{c}__stat_down"
        effect = f"{c}__cell_delta"
        sub = wide.loc[wide[pct] >= STRICT_PCT, ["gene", stat, effect]].copy()
        ascending_effect = direction == "down"
        sub = sub.sort_values([stat, effect], ascending=[True, ascending_effect], kind="mergesort")
        k = max(1, int(math.ceil(DECILE * len(sub))))
        sets.append(set(sub["gene"].iloc[:k].tolist()))
    hit = set.intersection(*sets) if sets else set()
    sub = wide.loc[wide["gene"].isin(hit)]
    if direction == "up":
        mask = np.ones(len(sub), dtype=bool)
        for c in COHORTS:
            mask &= sub[f"{c}__patient_delta"].to_numpy() > 0
    else:
        mask = np.ones(len(sub), dtype=bool)
        for c in COHORTS:
            mask &= sub[f"{c}__patient_delta"].to_numpy() < 0
    out = set(sub.loc[mask, "gene"].tolist())
    out.discard("CLDN4")
    return out


def passes_strict(row: pd.Series, direction: str) -> bool:
    for c in COHORTS:
        d = row[f"{c}__patient_delta"]
        if direction == "up":
            if d < STRICT_DELTA:
                return False
            if row[f"{c}__sign_up"] < STRICT_SIGN:
                return False
            if row[f"{c}__auc"] < 0.5 + STRICT_AUC:
                return False
            if row[f"{c}__pct_high"] < STRICT_PCT:
                return False
        else:
            if d > -STRICT_DELTA:
                return False
            if row[f"{c}__sign_down"] < STRICT_SIGN:
                return False
            if row[f"{c}__auc"] > 0.5 - STRICT_AUC:
                return False
            if row[f"{c}__pct_low"] < STRICT_PCT:
                return False
    return True


def samesign_mask(wide: pd.DataFrame, direction: str) -> pd.Series:
    if direction == "up":
        m = np.ones(len(wide), dtype=bool)
        for c in COHORTS:
            m &= wide[f"{c}__patient_delta"].to_numpy() > 0
        return pd.Series(m, index=wide.index)
    m = np.ones(len(wide), dtype=bool)
    for c in COHORTS:
        m &= wide[f"{c}__patient_delta"].to_numpy() < 0
    return pd.Series(m, index=wide.index)


def list_frame(wide: pd.DataFrame, genes: list[str], direction: str) -> pd.DataFrame:
    sub = wide.loc[wide["gene"].isin(genes)].copy()
    deltas = np.column_stack([sub[f"{c}__patient_delta"] for c in COHORTS])
    sub["min_patient_delta"] = deltas.min(axis=1)
    sub["max_patient_delta"] = deltas.max(axis=1)
    sub["direction"] = direction
    sub["mt_ribo"] = sub["gene"].map(is_mt_ribo)
    if direction == "up":
        sub = sub.sort_values(["min_patient_delta", "gene"], ascending=[False, True])
    else:
        sub = sub.sort_values(["max_patient_delta", "gene"], ascending=[True, True])
    front = ["gene", "direction", "mt_ribo", "min_patient_delta", "max_patient_delta"]
    rest = [c for c in sub.columns if c not in front]
    return sub[front + rest]


def add_wilcoxon(table_dir: Path, frame: pd.DataFrame, direction: str) -> pd.DataFrame:
    if frame.empty:
        return frame
    out = frame.copy()
    for c in COHORTS:
        npz = np.load(table_dir / f"{c}_patient_delta.npz", allow_pickle=True)
        genes = np.array([str(g).upper() for g in npz["genes"].tolist()])
        pos = {g: i for i, g in enumerate(genes)}
        deltas = npz["deltas"]
        pvals = []
        for g in out["gene"]:
            i = pos.get(g)
            if i is None:
                pvals.append(np.nan)
                continue
            d = deltas[i].astype(float)
            if np.allclose(d, 0):
                pvals.append(1.0)
                continue
            alt = "greater" if direction == "up" else "less"
            pvals.append(float(wilcoxon(d, alternative=alt, zero_method="wilcox").pvalue))
        out[f"{c}__wilcoxon_p"] = pvals
    return out


def write_list(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    sym = path.with_suffix(".symbols.txt")
    sym.write_text("\n".join(df["gene"].tolist()) + ("\n" if len(df) else ""))


def enrichr_or_fisher(genes: list[str], universe: list[str], out_stem: Path, label: str) -> pd.DataFrame:
    genes = [g for g in genes if g != "CLDN4"]
    if len(genes) < 5:
        print(f"  enrich skip {label}: n={len(genes)}", flush=True)
        return pd.DataFrame()
    try:
        import gseapy as gp
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=LIBS,
            organism="human",
            outdir=None,
            cutoff=1.0,
            no_plot=True,
            verbose=False,
        )
        res = enr.results.copy()
        res.insert(0, "input_list", label)
        res.to_csv(out_stem.with_suffix(".tsv"), sep="\t", index=False)
        print(f"  gseapy Enrichr {label} terms {len(res)}", flush=True)
        return res
    except Exception as exc:
        print(f"  gseapy failed ({exc}); Fisher on Enrichr GMT", flush=True)
        return fisher_gmt(genes, universe, out_stem, label)


def fetch_gmt(name: str, cache: Path) -> dict[str, set[str]]:
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / f"{name}.txt"
    if not dest.exists() or dest.stat().st_size < 1000:
        url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={name}"
        print(f"  GET {url}", flush=True)
        urllib.request.urlretrieve(url, dest)
    sets = {}
    for line in dest.read_text(errors="replace").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        term, genes = parts[0], parts[2:]
        if len(genes) == 1 and " " not in genes[0] and "," in genes[0]:
            genes = genes[0].split(",")
        sets[term] = {g.upper() for g in genes if g}
    return sets


def fisher_gmt(genes: list[str], universe: list[str], out_stem: Path, label: str) -> pd.DataFrame:
    gene_set = {g.upper() for g in genes}
    uni = {g.upper() for g in universe}
    rows = []
    cache = ROOT / "data" / "gmt"
    for name in LIBS:
        try:
            lib = fetch_gmt(name, cache)
        except Exception as exc:
            print(f"  GMT {name} failed: {exc}", flush=True)
            continue
        lib_genes = set().union(*lib.values()) if lib else set()
        background = uni & lib_genes if uni else lib_genes
        query = gene_set & background
        if len(query) < 5 or len(background) < 50:
            continue
        for term, members in lib.items():
            m = members & background
            if len(m) < 5:
                continue
            a = len(query & m)
            if a < 2:
                continue
            b = len(query) - a
            c = len(m) - a
            d = len(background) - a - b - c
            if min(a, b, c, d) < 0:
                continue
            oddsr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({
                "input_list": label,
                "Gene_set": name,
                "Term": term,
                "Overlap": f"{a}/{len(m)}",
                "P-value": p,
                "Odds Ratio": oddsr,
                "Genes": ";".join(sorted(query & m)),
                "method": "fisher_tested_universe",
            })
    if not rows:
        return pd.DataFrame()
    res = pd.DataFrame(rows)
    res["Adjusted P-value"] = res.groupby("Gene_set")["P-value"].transform(lambda s: bh(s.to_numpy()))
    res.to_csv(out_stem.with_suffix(".tsv"), sep="\t", index=False)
    return res


def load_arm_log(malig: Path, cohort: str, genes: list[str], n_cap: int = 2000) -> pd.DataFrame | None:
    """Cell x gene log1p CP10k for a stratified Q4/Q1 sample."""
    cdir = malig / cohort
    table_dir = ROOT / "results" / "tables"
    labels = pd.read_csv(table_dir / f"{cohort}_arm_labels.tsv.gz", sep="\t")
    if (cdir / "counts.npz").exists():
        X = sparse.load_npz(cdir / "counts.npz").tocsr()
        all_genes = pd.read_csv(cdir / "genes.tsv", header=None)[0].astype(str).str.upper().to_numpy()
        lib = np.load(cdir / "libsize.npy")
    else:
        from scipy.io import mmread
        X = mmread(cdir / "counts.mtx").tocsr()
        all_genes = pd.read_csv(cdir / "features.tsv", header=None)[0].astype(str).str.upper().to_numpy()
        lib = pd.read_csv(cdir / "libsize.tsv", sep="\t")["libsize"].to_numpy()
    if len(labels) != X.shape[1]:
        raise SystemExit(f"{cohort} label length {len(labels)} != {X.shape[1]}")
    idx = {g: i for i, g in enumerate(all_genes.tolist())}
    use = [g for g in genes if g in idx]
    if len(use) < 4:
        return None
    arm = labels["arm"].to_numpy()
    unit = labels["unit_id"].astype(str).to_numpy()
    keep_idx = []
    for uid in pd.unique(unit):
        for which in ("high", "low"):
            hit = np.flatnonzero((unit == uid) & (arm == which))
            if hit.size == 0:
                continue
            if hit.size > 80:
                hit = RNG.choice(hit, size=80, replace=False)
            keep_idx.extend(hit.tolist())
    keep_idx = np.array(keep_idx, dtype=int)
    if keep_idx.size > n_cap:
        keep_idx = RNG.choice(keep_idx, size=n_cap, replace=False)
    rows = [idx[g] for g in use]
    scale = 1e4 / np.maximum(lib[keep_idx], 1.0)
    sub = X[rows][:, keep_idx].toarray().astype(np.float64)
    sub = np.log1p(sub * scale)
    return pd.DataFrame(sub.T, columns=use)


def modules_for(malig: Path, genes: list[str], direction: str) -> pd.DataFrame:
    genes = [g for g in genes if g != "CLDN4"][:200]
    if len(genes) < 8:
        return pd.DataFrame()
    corrs = []
    used = None
    for c in COHORTS:
        mat = load_arm_log(malig, c, genes)
        if mat is None or mat.shape[0] < 50:
            print(f"  module matrix thin {c}", flush=True)
            return pd.DataFrame()
        if used is None:
            used = mat.columns.tolist()
        mat = mat.reindex(columns=used)
        # Constant genes (no variation in the sample) get correlation 0.
        with np.errstate(invalid="ignore"):
            corr = spearmanr(mat.to_numpy(), axis=0).correlation
        if np.ndim(corr) == 0:
            return pd.DataFrame()
        corr = np.nan_to_num(corr, nan=0.0)
        np.fill_diagonal(corr, 1.0)
        corrs.append(np.clip(corr, -0.999, 0.999))
    z = np.arctanh(np.stack(corrs, axis=0))
    mean_corr = np.tanh(z.mean(axis=0))
    np.fill_diagonal(mean_corr, 1.0)
    dist = np.clip(1.0 - mean_corr, 0, 2)
    np.fill_diagonal(dist, 0)
    dist = (dist + dist.T) / 2
    try:
        Z = linkage(squareform(dist, checks=False), method="average")
    except Exception as exc:
        print(f"  linkage failed: {exc}", flush=True)
        return pd.DataFrame()
    labels = fcluster(Z, t=0.5, criterion="distance")
    rows = []
    for lab in sorted(set(labels.tolist())):
        members = [used[i] for i, v in enumerate(labels.tolist()) if v == lab]
        if len(members) < 4:
            continue
        rows.append({
            "direction": direction,
            "module": f"{direction}_{lab}",
            "n_genes": len(members),
            "genes": ";".join(members),
        })
    return pd.DataFrame(rows)


def comet_pairs(malig: Path, table_dir: Path, genes: list[str]) -> pd.DataFrame:
    genes = [g for g in genes if g != "CLDN4"][:20]
    if len(genes) < 2:
        return pd.DataFrame()
    per = []
    for c in COHORTS:
        labels = pd.read_csv(table_dir / f"{c}_arm_labels.tsv.gz", sep="\t")
        cdir = malig / c
        if (cdir / "counts.npz").exists():
            X = sparse.load_npz(cdir / "counts.npz").tocsr()
            all_genes = pd.read_csv(cdir / "genes.tsv", header=None)[0].astype(str).str.upper().tolist()
            lib = np.load(cdir / "libsize.npy")
        else:
            from scipy.io import mmread
            X = mmread(cdir / "counts.mtx").tocsr()
            all_genes = pd.read_csv(cdir / "features.tsv", header=None)[0].astype(str).str.upper().tolist()
            lib = pd.read_csv(cdir / "libsize.tsv", sep="\t")["libsize"].to_numpy()
        idx = {g: i for i, g in enumerate(all_genes)}
        arm = labels["arm"].to_numpy()
        keep = np.flatnonzero((arm == "high") | (arm == "low"))
        is_high = arm[keep] == "high"
        scale = 1e4 / np.maximum(lib[keep], 1.0)
        binary = {}
        for g in genes:
            if g not in idx:
                binary[g] = None
                continue
            expr = np.asarray(X[idx[g], keep].todense()).ravel() * scale
            expr = np.log1p(expr)
            order = np.lexsort((is_high.astype(np.int8), -expr))
            ix = np.flatnonzero(is_high[order]).astype(np.uint16)
            stat, cutoff = xlmhg_cython.get_xlmhg_stat(
                ix, expr.size, int(ix.size), 1, expr.size, np.longdouble(1e-12)
            )
            if cutoff <= 0:
                binary[g] = np.zeros(expr.size, dtype=bool)
                continue
            thr = float(expr[order[cutoff - 1]])
            binary[g] = expr >= thr if thr > 0 else np.zeros(expr.size, dtype=bool)
        per.append((c, is_high, binary))
    rows = []
    N = None
    for i, a in enumerate(genes):
        for b in genes[i + 1:]:
            rec = {"gene_a": a, "gene_b": b}
            ok = True
            for c, is_high, binary in per:
                ba, bb = binary.get(a), binary.get(b)
                if ba is None or bb is None:
                    ok = False
                    break
                draw = ba & bb
                n_draw = int(draw.sum())
                k = int((draw & is_high).sum())
                N = int(is_high.size)
                K = int(is_high.sum())
                if n_draw == 0 or K == 0:
                    p, fold = 1.0, np.nan
                else:
                    p = float(hypergeom.sf(k - 1, N, K, n_draw))
                    fold = (k / n_draw) / (K / N)
                rec[f"{c}__p"] = p
                rec[f"{c}__fold"] = fold
                rec[f"{c}__k"] = k
                rec[f"{c}__draw"] = n_draw
                if not (p < PAIR_P and fold > PAIR_FOLD and k >= 10):
                    ok = False
            rec["concordant_4of4"] = bool(ok)
            rows.append(rec)
    return pd.DataFrame(rows)


def fmt(x, digits=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    ax = abs(float(x))
    if ax != 0 and ax < 1e-3:
        return f"{float(x):.2e}"
    return f"{float(x):.{digits}f}"


def write_finding(path: Path, control: pd.Series, lists: dict[str, pd.DataFrame], units: pd.DataFrame, enrich_notes: list[str]) -> None:
    def preview(df: pd.DataFrame, n=15) -> str:
        if df is None or df.empty:
            return "_none_\n"
        lines = ["| gene | min patient Δ | max patient Δ | mt/ribo |", "|---|---:|---:|---|"]
        for _, r in df.head(n).iterrows():
            lines.append(
                f"| {r['gene']} | {fmt(r['min_patient_delta'])} | {fmt(r['max_patient_delta'])} | {r['mt_ribo']} |"
            )
        if len(df) > n:
            lines.append(f"| … | {len(df) - n} more in the TSV | | |")
        return "\n".join(lines) + "\n"

    unit_txt = []
    for c in COHORTS:
        sub = units.loc[units["cohort"] == c]
        n_used = int(sub["used"].sum()) if "used" in sub.columns else len(sub)
        unit_txt.append(f"- {c}: {n_used} / {len(sub)} locked units entered Q4 vs Q1")
    lines = [
        "# CLDN4-high malignant markers, concordant-4",
        "",
        "ADDITIVE. CLDN4 only. Not a re-estimate of the locked patient-level T/NK correlation",
        "(ρ = −0.531, n = 65). These lists are genes that move with the CLDN4-high versus",
        "CLDN4-low split inside malignant cells, in every one of GSE123902, GSE131907,",
        "GSE205335, and GSE189357. They are co-markers of that state. They are not a",
        "surface-target ranking and they do not pin CLDN4 over EPCAM, CLDN7, or MUC1.",
        "",
        "The split is within-unit malignant Q4 versus Q1 of log1p(CP10k CLDN4).",
        "Units need ≥40 malignant cells and ≥10 CLDN4-positive malignant cells.",
        "P4001-style thin units and units with essentially no CLDN4 range are out.",
        "Patient / donor / sample remains the unit for the mean Δ. Cell counts below are",
        "the cells that entered the COMET ranking, not n.",
        "",
        "## Units used",
        "",
        *unit_txt,
        "",
        "## CLDN4 positive control",
        "",
        "CLDN4 itself is excluded from the marker lists. In every used unit the Q4 mean",
        "is above the Q1 mean, so the within-unit sign fraction is 1.",
        "",
        "| cohort | units | patient Δ | cell Δ | AUC |",
        "|---|---:|---:|---:|---:|",
    ]
    for c in COHORTS:
        lines.append(
            f"| {c} | {int(control[f'{c}__n_units'])} | {fmt(control[f'{c}__patient_delta'])} | "
            f"{fmt(control[f'{c}__cell_delta'])} | {fmt(control[f'{c}__auc'])} |"
        )
    lines += [
        "",
        "## 4/4 lists",
        "",
        f"- Strict patient-effect UP (Δ≥{STRICT_DELTA}, sign fraction≥{STRICT_SIGN}, AUC≥0.55, detection≥{STRICT_PCT} in all 4): **{len(lists['strict_up'])}** genes.",
        f"- Strict patient-effect DOWN: **{len(lists['strict_down'])}** genes.",
        f"- COMET top-decile UP (XL-mHG decile and patient Δ>0 in all 4): **{len(lists['comet_up'])}** genes.",
        f"- COMET top-decile DOWN: **{len(lists['comet_down'])}** genes.",
        f"- Same-sign patient Δ UP, no effect floor, ranked by the weakest cohort: **{len(lists['samesign_up'])}** genes.",
        f"- Same-sign patient Δ DOWN: **{len(lists['samesign_down'])}** genes.",
        "",
        "Strict UP (strongest minimum Δ first):",
        "",
        preview(lists["strict_up"]),
        "",
        "Strict DOWN:",
        "",
        preview(lists["strict_down"]),
        "",
        "Files: `results/gene_lists/`.",
        "",
        "## Modules and enrichment",
        "",
        "Differential modules are average-linkage clusters of the COMET/strict marker",
        "genes (distance = 1 − mean cross-cohort Spearman, cut at 0.5, modules smaller",
        "than 4 genes dropped). Enrichment is gseapy Enrichr when the API responds,",
        "otherwise a Fisher exact test on the same Enrichr libraries using the genes",
        "tested in all four cohorts as the universe.",
        "",
        *enrich_notes,
        "",
        "Marker-malignant gates in GSE123902 and GSE189357 require EPCAM or a KRT8/18/19",
        "count above zero. A 4/4 gene still has to clear the two author-label cohorts",
        "(GSE131907, GSE205335), so the gate alone cannot put a gene on the list.",
        "",
    ]
    path.write_text("\n".join(lines))


def heatmap(df: pd.DataFrame, path: Path, title: str) -> None:
    if df is None or df.empty:
        return
    show = df.head(25)
    mat = np.column_stack([show[f"{c}__patient_delta"].to_numpy() for c in COHORTS])
    fig, ax = plt.subplots(figsize=(7.2, max(3.5, 0.28 * len(show) + 1.2)))
    lim = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1
    lim = max(float(lim), 0.1)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax.set_xticks(range(4), ["123902", "131907", "205335", "189357"], rotation=0)
    ax.set_yticks(range(len(show)), show["gene"].tolist(), fontsize=8)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="patient mean Δ log1p CP10k")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--malig", default="/tmp/concordant4_malig")
    ap.add_argument("--tables", default=str(ROOT / "results" / "tables"))
    args = ap.parse_args()
    table_dir = Path(args.tables)
    malig = Path(args.malig)
    glist = ROOT / "results" / "gene_lists"
    figdir = ROOT / "results" / "figures"
    wide = load_stats(table_dir)
    wide = wide.loc[wide["gene"] != ""].copy()
    print(f"genes in all 4 cohorts: {len(wide)}", flush=True)
    control = wide.loc[wide["gene"] == "CLDN4"]
    if control.empty:
        raise SystemExit("CLDN4 missing from the intersection")
    control_row = control.iloc[0]
    for c in COHORTS:
        if not (control_row[f"{c}__patient_delta"] > 0):
            raise SystemExit(f"CLDN4 delta not positive in {c}")

    lists = {}
    for direction, key in (("up", "up"), ("down", "down")):
        strict_genes = []
        for _, row in wide.iterrows():
            if row["gene"] == "CLDN4":
                continue
            if passes_strict(row, direction):
                strict_genes.append(row["gene"])
        comet = decile_genes(wide, direction)
        same = wide.loc[samesign_mask(wide, direction) & (wide["gene"] != "CLDN4"), "gene"].tolist()
        lists[f"strict_{key}"] = add_wilcoxon(table_dir, list_frame(wide, strict_genes, direction), direction)
        lists[f"comet_{key}"] = add_wilcoxon(table_dir, list_frame(wide, sorted(comet), direction), direction)
        lists[f"samesign_{key}"] = list_frame(wide, same, direction)
        write_list(lists[f"strict_{key}"], glist / f"patient_effect_4of4_{direction.upper()}_strict.tsv")
        write_list(lists[f"comet_{key}"], glist / f"comet_top_decile_4of4_{direction.upper()}.tsv")
        write_list(lists[f"samesign_{key}"], glist / f"patient_samesign_4of4_{direction.upper()}_ranked.tsv")
        print(
            f"{direction} strict {len(lists[f'strict_{key}'])} "
            f"comet {len(lists[f'comet_{key}'])} samesign {len(lists[f'samesign_{key}'])}",
            flush=True,
        )

    # No-ribo companion of whichever strict list is non-empty, pre-specified sensitivity.
    for direction in ("up", "down"):
        src = lists[f"strict_{direction}"]
        if src.empty:
            continue
        write_list(src.loc[~src["mt_ribo"].astype(bool)], glist / f"patient_effect_4of4_{direction.upper()}_strict_no_mt_ribo.tsv")

    universe = wide["gene"].tolist()
    enrich_notes = []
    for direction in ("up", "down"):
        strict = lists[f"strict_{direction}"]
        comet = lists[f"comet_{direction}"]
        same = lists[f"samesign_{direction}"]
        if len(strict) >= 15:
            query, qname = strict["gene"].tolist(), f"strict_{direction}"
        else:
            query, qname = same["gene"].head(100).tolist(), f"samesign_top100_{direction}"
        res = enrichr_or_fisher(query, universe, glist / f"enrichr_{qname}", qname)
        n_sig = 0 if res.empty or "Adjusted P-value" not in res.columns else int((res["Adjusted P-value"] < 0.05).sum())
        enrich_notes.append(f"- Enrichr input `{qname}` n={len(query)}; terms with adjusted P<0.05: {n_sig}.")
        if len(comet) >= 15 and qname != f"comet_{direction}":
            res2 = enrichr_or_fisher(comet["gene"].tolist(), universe, glist / f"enrichr_comet_{direction}", f"comet_{direction}")
            n2 = 0 if res2.empty or "Adjusted P-value" not in res2.columns else int((res2["Adjusted P-value"] < 0.05).sum())
            enrich_notes.append(f"- Enrichr input `comet_{direction}` n={len(comet)}; adjusted P<0.05: {n2}.")
        # Module source: strict if it has enough genes, else COMET decile, else top same-sign.
        if len(strict) >= 8:
            mod_genes = strict["gene"].tolist()
            mod_from = "strict"
        elif len(comet) >= 8:
            mod_genes = comet["gene"].tolist()
            mod_from = "comet_decile"
        else:
            mod_genes = same["gene"].head(60).tolist()
            mod_from = "samesign_top60"
        print(f"  modules from {mod_from} n={len(mod_genes)}", flush=True)
        mods = modules_for(malig, mod_genes, direction)
        if not mods.empty:
            mods.insert(0, "source", mod_from)
            mods.to_csv(glist / f"modules_{direction}.tsv", sep="\t", index=False)
            enrich_notes.append(
                f"- {direction} modules from `{mod_from}`: " + ", ".join(f"{r.module} (n={r.n_genes})" for r in mods.itertuples())
            )
            for r in mods.itertuples():
                members = r.genes.split(";")
                if len(members) >= 8:
                    enrichr_or_fisher(members, universe, glist / f"enrichr_module_{r.module}", r.module)
        else:
            enrich_notes.append(f"- {direction} modules from `{mod_from}`: no cluster reached 4 genes.")

    pair_source = lists["comet_up"]
    if len(pair_source) < 2:
        pair_source = lists["samesign_up"]
    pairs = comet_pairs(malig, table_dir, pair_source["gene"].tolist())
    if not pairs.empty:
        pairs.to_csv(glist / "comet_pairs_top20.tsv", sep="\t", index=False)
        n4 = int(pairs["concordant_4of4"].sum())
        enrich_notes.append(f"- COMET AND-pairs among the top {min(20, len(pair_source))} UP genes: {n4} / {len(pairs)} pass p<{PAIR_P}, fold>{PAIR_FOLD}, k≥10 in all 4 cohorts.")
    else:
        enrich_notes.append("- COMET pairs: not run (fewer than 2 UP genes).")

    units = pd.concat([pd.read_csv(table_dir / f"{c}_units_used.tsv", sep="\t") for c in COHORTS], ignore_index=True)
    write_finding(ROOT / "FINDING.md", control_row, lists, units, enrich_notes)
    heatmap(lists["strict_up"] if len(lists["strict_up"]) else lists["comet_up"].head(25), figdir / "heatmap_up_patient_delta.png", "CLDN4-high UP markers, patient Δ")
    heatmap(lists["strict_down"] if len(lists["strict_down"]) else lists["comet_down"].head(25), figdir / "heatmap_down_patient_delta.png", "CLDN4-high DOWN markers, patient Δ")
    control.to_csv(glist / "CLDN4_positive_control.tsv", sep="\t", index=False)
    print("wrote FINDING and gene lists", flush=True)


if __name__ == "__main__":
    main()
