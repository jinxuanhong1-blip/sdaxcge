#!/usr/bin/env python3
"""Pairwise and leave-one-out Harmony. Do not force one all-series object.

User A3 direction (pre-specified, not tuned after looking):
  A) malignant-like TACSTD2 NMPR median > MPR median (n>=2 per arm)
  B) TACSTD2 vs T/NK Spearman rho < 0 (n>=4)
KEEP a combo if A or B holds. Extra figures only for KEEP.
GSE205335 is RECIST, not MPR — it can enter T/NK tests only.
"""
from __future__ import annotations

import gzip
import json
import shutil
import sys
import tempfile
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = DATA / "cache"
FIG = ROOT / "figures" / "combos"
TAB = ROOT / "tables"
SCR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCR))

import run_all as ra  # noqa: E402
from gene_sets import TARGETS, core_panel  # noqa: E402

CACHE.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SERIES = ["GSE207422", "GSE241934", "GSE291670", "GSE205335"]


def _leiden_only(Z: np.ndarray, n_neighbors: int = 15, res: float = 0.8) -> np.ndarray:
    import igraph as ig
    import leidenalg
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    nn.fit(Z)
    dist, ind = nn.kneighbors()
    edges, w = [], []
    n = Z.shape[0]
    for i in range(n):
        for d, j in zip(dist[i, 1:], ind[i, 1:]):
            edges.append((i, int(j)))
            w.append(float(1.0 / (d + 1e-6)))
    g = ig.Graph(n=n, edges=edges, directed=False)
    g.es["weight"] = w
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition, weights="weight", resolution_parameter=res, seed=7
    )
    return np.array(part.membership)


def _save_ds(name: str, X: np.ndarray, obs: pd.DataFrame, genes: list[str]) -> None:
    np.savez_compressed(CACHE / f"{name}.npz", X=X, nUMI=obs["nUMI"].to_numpy(np.float32), genes=np.array(genes))
    obs.to_pickle(CACHE / f"{name}_obs.pkl")
    print(f"cached {name} cells={len(obs)} genes={len(genes)}", flush=True)


def _load_ds(name: str) -> tuple[np.ndarray, pd.DataFrame, list[str]]:
    z = np.load(CACHE / f"{name}.npz", allow_pickle=True)
    obs = pd.read_pickle(CACHE / f"{name}_obs.pkl")
    return z["X"], obs, [str(g) for g in z["genes"]]


def cache_three() -> None:
    if (CACHE / "GSE207422.npz").exists() and (CACHE / "GSE241934.npz").exists() and (CACHE / "GSE291670.npz").exists():
        print("three-series cache present", flush=True)
        return
    shared_path = DATA / "GSE207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"

    def tsv_genes(path):
        with gzip.open(path, "rt") as f:
            f.readline()
            return [ln.split("\t", 1)[0] for ln in f]

    g207 = set(tsv_genes(shared_path))
    g241 = set(ra._feat_symbols(DATA / "GSE241934" / "GSE241934_IIT_features.tsv.gz"))
    g291 = set(ra._feat_symbols(next((DATA / "GSE291670" / "raw").glob("*MPR-1_features.tsv.gz"))))
    shared = sorted(g207 & g241 & g291)
    X207, obs207, genes207 = ra.load_gse207422(shared)
    gene_order = []
    seen = set()
    for g in genes207:
        if g in set(shared) and g not in seen:
            seen.add(g)
            gene_order.append(g)
    for g in TARGETS + core_panel():
        if g in set(shared) and g not in seen:
            seen.add(g)
            gene_order.append(g)
    X207 = ra.reindex_X(X207, genes207, gene_order)
    _save_ds("GSE207422", X207, obs207, gene_order)
    del X207

    Xi, oi, gi = ra.load_gse241934("IIT", set(shared), gene_order)
    Xr, orr, gr = ra.load_gse241934("REAL", set(shared), gene_order)
    X241 = np.vstack([ra.reindex_X(Xi, gi, gene_order), ra.reindex_X(Xr, gr, gene_order)])
    obs241 = pd.concat([oi, orr], ignore_index=True)
    _save_ds("GSE241934", X241, obs241, gene_order)
    del X241, Xi, Xr

    X291, o291, g291e = ra.load_gse291670(set(shared), gene_order)
    _save_ds("GSE291670", ra.reindex_X(X291, g291e, gene_order), o291, gene_order)


def _ungzip_until_rds(src: Path, dest: Path) -> None:
    current = src
    tmps = []
    for _ in range(3):
        with open(current, "rb") as fh:
            magic = fh.read(2)
        if magic == b"\x1f\x8b":
            nxt = dest.with_suffix(dest.suffix + f".pass{len(tmps)}")
            with gzip.open(current, "rb") as zin, open(nxt, "wb") as zout:
                shutil.copyfileobj(zin, zout, 16 * 1024 * 1024)
            tmps.append(nxt)
            current = nxt
            continue
        break
    if current != dest:
        shutil.copyfile(current, dest)
    for p in tmps:
        if p.exists() and p != dest:
            p.unlink(missing_ok=True)


def cache_gse205335(gene_order: list[str]) -> dict:
    """Extract shared-panel genes from the public dgCMatrix RDS. Honest skip if parse fails."""
    out = {"included": False, "why": ""}
    dest_npz = CACHE / "GSE205335.npz"
    if dest_npz.exists():
        out["included"] = True
        out["why"] = "cached"
        return out
    rds_gz = DATA / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    if not rds_gz.exists():
        out["why"] = "RDS not downloaded"
        return out
    try:
        import rdata
        from scipy import sparse
    except Exception as e:
        out["why"] = f"rdata/scipy unavailable: {e}"
        return out

    smap = pd.read_csv(TAB / "gse205335_sample_map.tsv", sep="\t")
    smap = smap.set_index("orig.ident")
    ident = pd.read_csv(DATA / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t", dtype=str)
    print("GSE205335 decompress+parse RDS", flush=True)
    try:
        with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
            rds_path = Path(tmp) / "matrix.rds"
            _ungzip_until_rds(rds_gz, rds_path)
            print(f"  rds bytes={rds_path.stat().st_size}", flush=True)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                obj = rdata.read_rds(rds_path)
        if not {"i", "p", "Dim", "Dimnames", "x"}.issubset(vars(obj)):
            out["why"] = f"RDS is not dgCMatrix; keys={list(vars(obj))}"
            return out
        genes = np.asarray(obj.Dimnames[0], dtype=str)
        barcodes = np.asarray(obj.Dimnames[1], dtype=str)
        matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
        print(f"  matrix {matrix.shape}", flush=True)
        gix = {g: i for i, g in enumerate(genes)}
        present = [g for g in gene_order if g in gix]
        if "TACSTD2" not in present or "CLDN4" not in present:
            out["why"] = "TACSTD2/CLDN4 missing from RDS gene names"
            return out
        nUMI = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float32)
        X = np.zeros((matrix.shape[1], len(gene_order)), dtype=np.float32)
        col = {g: j for j, g in enumerate(gene_order)}
        for g in present:
            X[:, col[g]] = np.asarray(matrix.getrow(int(gix[g])).todense()).ravel()
        del matrix
        # align identity
        idf = ident.set_index("barcode")
        # barcodes in RDS may differ slightly
        n = len(barcodes)
        orig = []
        lin = []
        for bc in barcodes:
            if bc in idf.index:
                orig.append(str(idf.loc[bc, "orig.ident"]))
                lin.append(str(idf.loc[bc, "lineage.sub"]))
            else:
                orig.append("")
                lin.append("")
        map_lin = {
            "Malignant cells": "epithelial",
            "Normal epithelial cells": "epithelial",
            "CD4+ T cells": "t",
            "CD8+ T cells": "t",
            "NK cells": "nk",
            "B/Plasma cells": "b",
            "Myeloid cells": "myeloid",
        }
        rows = []
        for i in range(n):
            oid = orig[i]
            if oid in smap.index:
                r = smap.loc[oid]
                patient = str(r["patient"])
                histo = ra.map_histo(r["cancer_subtype"])
                is_tumor = bool(r["is_tumor"])
            else:
                patient = oid or f"unk{i}"
                histo = "unknown"
                is_tumor = True
            al = map_lin.get(lin[i], "")
            rows.append(
                {
                    "barcode": barcodes[i],
                    "sample_id": patient,
                    "patient_id": patient,
                    "dataset": "GSE205335",
                    "cohort": "GSE205335",
                    "mpr": np.nan,  # RECIST is not MPR
                    "pathologic_response": "RECIST_not_MPR",
                    "histology": histo,
                    "is_post": is_tumor,
                    "resource": "palliative ICI (RECIST)",
                    "nUMI": float(nUMI[i]),
                    "author_lineage": al,
                }
            )
        obs = pd.DataFrame(rows)
        _save_ds("GSE205335", X, obs, gene_order)
        out["included"] = True
        out["why"] = f"extracted {n} cells; genes present {len(present)}/{len(gene_order)}"
        out["n_cells"] = int(n)
        out["n_genes_present"] = int(len(present))
        return out
    except Exception as e:
        out["why"] = f"RDS parse failed: {type(e).__name__}: {e}"
        print("GSE205335 FAIL", out["why"], flush=True)
        return out


def combo_names() -> list[tuple[str, tuple[str, ...]]]:
    names = []
    for a, b in combinations(SERIES, 2):
        names.append((f"{a}+{b}", (a, b)))
    for left in SERIES:
        keep = tuple(s for s in SERIES if s != left)
        names.append((f"LOO_minus_{left}", keep))
    return names


def run_one(combo_id: str, members: tuple[str, ...], available: dict) -> dict:
    print(f"\n=== {combo_id} ===", flush=True)
    Xs, obss, genes0 = [], [], None
    for s in members:
        if s not in available:
            return {"combo": combo_id, "members": list(members), "ran": False, "why": f"{s} not available"}
        X, obs, genes = available[s]
        if genes0 is None:
            genes0 = genes
        Xs.append(ra.reindex_X(X, genes, genes0))
        obss.append(obs)
    X = np.vstack(Xs)
    obs = pd.concat(obss, ignore_index=True)
    nUMI = obs["nUMI"].to_numpy(np.float64)
    nUMI = np.where(np.isfinite(nUMI) & (nUMI > 0), nUMI, X.sum(axis=1) + 1)
    L = ra.log_cp10k(X, nUMI)
    del X
    marker_lin = ra.assign_marker_lineage(L, genes0)
    al = obs["author_lineage"].astype(str)
    auth = al.where(al.notna() & (al != "") & (al != "nan"), other=np.nan)
    obs["native_lineage"] = np.where(auth.notna(), auth, marker_lin)

    take = []
    for _, sub in obs.groupby("sample_id"):
        k = min(ra.MAX_PER_SAMPLE, len(sub))
        take.extend(sub.sample(k, random_state=7).index.to_list())
    take = np.array(sorted(take))
    Ls = L[take]
    mu = Ls.mean(axis=0)
    sd = Ls.std(axis=0)
    sd[sd < 1e-6] = 1.0
    Zall = (L - mu) / sd
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors

    pca = PCA(n_components=min(ra.N_PCS, Zall.shape[1] - 1), random_state=7)
    pcs_s = pca.fit_transform(Zall[take])
    Z_h = ra.run_harmony(pcs_s, obs.loc[take, "dataset"].astype(str).tolist())
    cl = _leiden_only(Z_h)
    uv = None
    cl_map = ra.label_clusters(cl, obs.loc[take, "native_lineage"].to_numpy())
    harm_lin_s = np.array([cl_map[int(c)] for c in cl])
    pcs_all = pca.transform(Zall)
    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(pcs_s)
    _, ind = nn.kneighbors(pcs_all)
    obs["harmony_lineage"] = harm_lin_s[ind[:, 0]]
    nl = ra.mean_genes(L, genes0, ra.NORMAL_LUNG)
    is_epi = obs["harmony_lineage"].to_numpy() == "epithelial"
    is_mal = np.zeros(len(obs), dtype=bool)
    if is_epi.sum() >= 20:
        thr = np.quantile(nl[is_epi], 0.75)
        is_mal[is_epi & (nl <= thr)] = True
    else:
        is_mal[is_epi] = True
    obs["is_malignant"] = is_mal
    obs["is_tnk"] = obs["harmony_lineage"].isin(["t", "nk"]).to_numpy()
    per = ra.patient_table(obs, L, genes0)
    # MPR test: only rows with MPR/NMPR (205335 excluded automatically)
    gene = "TACSTD2"
    a = ra.eligible(per, gene)
    nmpr = a.loc[a.mpr == "NMPR", f"{gene}_mean_log1p_cp10k"].to_numpy()
    mpr = a.loc[a.mpr == "MPR", f"{gene}_mean_log1p_cp10k"].to_numpy()
    mpr_st = ra.mwu_report(nmpr, mpr)
    ols = ra.ols_mpr_cohort(a[f"{gene}_mean_log1p_cp10k"], a["mpr"], a["cohort"]) if len(a) else {"n": 0}
    b = per[(per["n_malignant"] >= ra.MIN_MAL) & (per["n_tnk"] >= ra.MIN_TNK) & per[f"{gene}_mean_log1p_cp10k"].notna()]
    # T/NK: include 205335 patients (no MPR required)
    tnk = ra.spearman_report(b[f"{gene}_mean_log1p_cp10k"], b["frac_tnk"])
    part = ra.partial_spearman(b[f"{gene}_mean_log1p_cp10k"], b["frac_tnk"], b["cohort"]) if len(b) else {"n": 0, "rho": np.nan, "p": np.nan}

    dir_mpr = bool(len(nmpr) >= 2 and len(mpr) >= 2 and np.median(nmpr) > np.median(mpr))
    dir_tnk = bool(tnk.get("n", 0) >= 4 and np.isfinite(tnk.get("rho", np.nan)) and tnk["rho"] < 0)
    keep = dir_mpr or dir_tnk

    rec = {
        "combo": combo_id,
        "members": list(members),
        "ran": True,
        "n_cells": int(len(obs)),
        "n_embed": int(len(take)),
        "n_patients": int(per["patient_id"].nunique()),
        "n_mpr": int(mpr_st.get("n_y", 0) or 0),
        "n_nmpr": int(mpr_st.get("n_x", 0) or 0),
        "TACSTD2_median_NMPR": mpr_st.get("median_x"),
        "TACSTD2_median_MPR": mpr_st.get("median_y"),
        "TACSTD2_delta_NMPR_minus_MPR": mpr_st.get("delta_median_x_minus_y"),
        "TACSTD2_mpr_p": mpr_st.get("p"),
        "TACSTD2_ols_beta": ols.get("coef_NMPR_vs_MPR"),
        "TACSTD2_ols_p": ols.get("p"),
        "TACSTD2_tnk_n": tnk.get("n"),
        "TACSTD2_tnk_rho": tnk.get("rho"),
        "TACSTD2_tnk_p": tnk.get("p"),
        "TACSTD2_tnk_partial_rho": part.get("rho"),
        "TACSTD2_tnk_partial_p": part.get("p"),
        "direction_NMPR_gt_MPR": dir_mpr,
        "direction_rho_negative": dir_tnk,
        "KEEP": keep,
        "why_keep": "NMPR>MPR" if dir_mpr else ("rho<0" if dir_tnk else "neither"),
    }
    per.to_csv(TAB / f"combo_{combo_id}_per_sample.tsv", sep="\t", index=False)
    if keep:
        import umap

        uv = umap.UMAP(n_neighbors=15, min_dist=0.3, metric="euclidean", random_state=7).fit_transform(Z_h)
        obs_s = obs.loc[take].copy()
        obs_s["harmony_lineage"] = harm_lin_s
        c4 = ra.eligible(per, "CLDN4")
        stats_d = {
            "TACSTD2_mpr": mpr_st,
            "CLDN4_mpr": ra.mwu_report(
                c4.loc[c4.mpr == "NMPR", "CLDN4_mean_log1p_cp10k"].to_numpy(),
                c4.loc[c4.mpr == "MPR", "CLDN4_mean_log1p_cp10k"].to_numpy(),
            ),
            "TACSTD2_tnk": tnk,
            "CLDN4_tnk": ra.spearman_report(b["CLDN4_mean_log1p_cp10k"], b["frac_tnk"]) if len(b) else {"n": 0, "rho": np.nan, "p": np.nan},
        }
        outdir = FIG / combo_id
        outdir.mkdir(parents=True, exist_ok=True)
        old_fig = ra.FIG
        ra.FIG = outdir
        try:
            ra.make_figures(obs_s.reset_index(drop=True), uv, per, stats_d)
        finally:
            ra.FIG = old_fig
        rec["figures"] = str(outdir)
    print(
        f"  cells={rec['n_cells']} mpr {rec['n_nmpr']} vs {rec['n_mpr']} "
        f"Δ={rec['TACSTD2_delta_NMPR_minus_MPR']} p={rec['TACSTD2_mpr_p']} "
        f"ρ={rec['TACSTD2_tnk_rho']} KEEP={keep}",
        flush=True,
    )
    return rec


def write_finding(rows: list[dict], gse205: dict) -> None:
    keepers = [r for r in rows if r.get("KEEP")]
    lines = [
        "# Pairwise / leave-one-out Harmony — FINDING",
        "",
        "**Public data only. Additive to User A3. The GSE207422 CopyKAT slide was not re-run.**",
        "**Change of plan:** one forced Harmony of every series is not the claim. Each combo was actually integrated.",
        "",
        "## 一句话结论 / TL;DR",
    ]
    if keepers:
        bits = []
        for r in keepers:
            bits.append(
                f"{r['combo']}: "
                + (
                    f"NMPR>MPR Δ={r.get('TACSTD2_delta_NMPR_minus_MPR'):.2f} p={r.get('TACSTD2_mpr_p'):.3g} n={r.get('n_nmpr')} vs {r.get('n_mpr')}"
                    if r.get("direction_NMPR_gt_MPR")
                    else ""
                )
                + (
                    f"{'; ' if r.get('direction_NMPR_gt_MPR') else ''}ρ={r.get('TACSTD2_tnk_rho'):.2f} p={r.get('TACSTD2_tnk_p'):.3g} n={r.get('TACSTD2_tnk_n')}"
                    if r.get("direction_rho_negative")
                    else ""
                )
            )
        lines.append(
            "Combinations that hold the User A3 *direction* (NMPR>MPR and/or ρ<0), not necessarily p<0.05: "
            + " | ".join(bits)
            + ". Extra figures only for these. This does not replace User A3."
        )
    else:
        lines.append(
            "No pairwise or leave-one-out Harmony held the User A3 direction "
            "(TACSTD2 NMPR>MPR or ρ<0 vs T/NK) at the pre-specified floors. "
            "No extra combo figures. This does not replace User A3."
        )
    lines += [
        "",
        "## GSE205335",
        f"- Parse: {gse205.get('why')}",
        "- RECIST is **not** MPR. In any combo that includes it, MPR tests use only neoadjuvant samples; T/NK may use 205335 patients.",
        "",
        "## All combos (honest)",
        "| Combo | cells | patients | NMPR vs MPR n | TACSTD2 Δ (NMPR−MPR) | p | ρ vs T/NK | n | p | KEEP |",
        "|---|---:|---:|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in rows:
        if not r.get("ran"):
            lines.append(f"| {r['combo']} | — | — | not run | — | — | — | — | — | no ({r.get('why','')}) |")
            continue
        lines.append(
            f"| {r['combo']} | {r.get('n_cells',0):,} | {r.get('n_patients','')} | "
            f"{r.get('n_nmpr')} vs {r.get('n_mpr')} | "
            f"{r.get('TACSTD2_delta_NMPR_minus_MPR')} | {r.get('TACSTD2_mpr_p')} | "
            f"{r.get('TACSTD2_tnk_rho')} | {r.get('TACSTD2_tnk_n')} | {r.get('TACSTD2_tnk_p')} | "
            f"{'YES ('+r.get('why_keep','')+')' if r.get('KEEP') else 'no'} |"
        )
    lines += [
        "",
        "## Rule (pre-specified)",
        "- KEEP if TACSTD2 NMPR median > MPR median (n≥2/arm) **or** TACSTD2 vs T/NK ρ<0 (n≥4).",
        "- Extra figures only for KEEP. Cell-level p-values are not used.",
        "- The previous three-series object (429,593 cells) is the LOO_minus_GSE205335 row, re-run here for the same rule.",
        "",
        "## Files",
        "- `tables/combo_screen.json`, `tables/combo_screen.tsv`",
        "- `figures/combos/<KEEP_id>/` only if KEEP",
    ]
    (ROOT / "FINDING_combos.md").write_text("\n".join(lines) + "\n")


def write_snippet(rows: list[dict]) -> None:
    keepers = [r for r in rows if r.get("KEEP")]
    if not keepers:
        txt = """# Paper snippet — pairwise/LOO Harmony (additive)

*Placement: extra note after User A3. Do **not** replace the GSE207422 CopyKAT slide.*

Pairwise and leave-one-out Harmony of public processed neoadjuvant (and, where parseable, GSE205335) scRNA matrices were run instead of one forced four-series object. No combination held the User A3 direction (malignant-like *TACSTD2* NMPR>MPR or ρ<0 versus T/NK) at the pre-specified patient floors. Honest n/ρ/p are in `FINDING_combos.md`. This does not replace User A3.
"""
    else:
        parts = []
        for r in keepers:
            parts.append(
                f"{r['combo']} ({r['n_cells']:,} cells, {r['n_patients']} patients): "
                f"TACSTD2 NMPR vs MPR Δ={r.get('TACSTD2_delta_NMPR_minus_MPR')} "
                f"(n={r.get('n_nmpr')} vs {r.get('n_mpr')}, p={r.get('TACSTD2_mpr_p')}); "
                f"vs T/NK ρ={r.get('TACSTD2_tnk_rho')} (n={r.get('TACSTD2_tnk_n')}, p={r.get('TACSTD2_tnk_p')})"
            )
        txt = (
            "# Paper snippet — pairwise/LOO Harmony (additive)\n\n"
            "*Placement: extra panel after User A3. Do **not** replace the GSE207422 CopyKAT slide.*\n\n"
            "Instead of one forced Harmony of every neoadjuvant series, pairwise and leave-one-out "
            "integrations were run on public processed UMI matrices. Combinations that held the "
            "User A3 *direction* (NMPR>MPR and/or negative TACSTD2–T/NK ρ), not necessarily p<0.05: "
            + "; ".join(parts)
            + ". Extra figures only for these combinations. GSE205335 (RECIST, not MPR) was used only in T/NK tests when the RDS parsed. This does not replace User A3.\n"
        )
    (ROOT / "paper_snippet_combos.md").write_text(txt)


def main() -> int:
    cache_three()
    X0, _, genes0 = _load_ds("GSE207422")
    del X0
    gse205 = cache_gse205335(genes0)
    print("GSE205335", gse205, flush=True)

    available = {}
    for s in SERIES:
        if (CACHE / f"{s}.npz").exists():
            available[s] = _load_ds(s)
            print(f"load {s} {available[s][1].shape[0]}", flush=True)

    rows = []
    for cid, members in combo_names():
        rows.append(run_one(cid, members, available))

    (TAB / "combo_screen.json").write_text(json.dumps({"gse205335": gse205, "combos": rows}, indent=2, default=str))
    pd.DataFrame(rows).to_csv(TAB / "combo_screen.tsv", sep="\t", index=False)
    write_finding(rows, gse205)
    write_snippet(rows)
    print("KEEP", [r["combo"] for r in rows if r.get("KEEP")], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
