#!/usr/bin/env python3
"""B6 radius rework: tumor-spot CLDN4 vs immune neighborhoods at hex rings 1/2/3.

Pre-specified (not tuned to PR67 or the user claim):
  * Visium E-MTAB-13530: tumor sections only (P*_T*).
  * Index spots = epithelial-high AND not immune-rich (immune score < section Q3).
    CLDN4/TACSTD2 are scored only on those index spots.
  * Neighborhoods = Visium hex rings at distance 1, 2, 3 (self excluded).
    Neighborhood immune uses all QC spots in the ring (immune-rich included).
  * GeoMx GSE271689: tumor / PanCK compartment only. No stromal pairing.
    GeoMx has no hex grid, so this is same-AOI co-expression, not rings.

Do not invent results. Do not tune thresholds to enlarge |rho|.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import tarfile
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

TARGETS = ["CLDN4", "TACSTD2"]
# Epithelial score excludes TACSTD2/CLDN4 and close paralogs (CLDN3/7).
EPI_MARKERS = ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "KRT17", "CDH1", "ELF3", "SFN"]
IMMUNE_MARKERS = [
    "PTPRC", "CD3E", "CD3D", "CD2", "TRAC", "IL7R", "CD8A", "CCL5", "NKG7",
    "GZMB", "GZMA", "GZMK", "PRF1", "CD79A", "MS4A1", "LYZ", "C1QA", "C1QB",
    "CD68", "AIF1", "FCER1G", "CD14", "CD52",
]
TEFF_MARKERS = [
    "CD8A", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "IFNG",
    "CXCL9", "CXCL10", "CXCL13", "CD3E", "CD3D", "CD2", "TRAC",
]
TUMOR_SECTIONS = [
    "P10_T1", "P10_T2", "P10_T3", "P10_T4",
    "P11_T1", "P11_T2", "P11_T3", "P11_T4",
    "P15_T1", "P15_T2",
    "P16_T1", "P16_T2",
    "P17_T1", "P17_T2",
    "P19_T1", "P19_T2",
    "P24_T1", "P24_T2",
    "P25_T1", "P25_T2",
]
MIN_COUNTS = 250
MIN_INDEX = 40
MIN_NEIGHBORS = 3
IMMUNE_RICH_Q = 0.75  # pre-specified; not tuned
RINGS = (1, 2, 3)
USER_CLAIM = (
    "CLDN4-high tumor domains spatially avoid immune niches "
    "(expected: significant spatial exclusion, i.e. negative rho of meaningful size)"
)
PR67_MISMATCH = "Visium CLDN4 vs immune neighborhood |median partial rho| <= 0.06"


def visium_hex_distance(r1, c1, r2, c2) -> int:
    """Hex distance on the Visium double-width (array_row, array_col) grid."""
    dr = abs(int(r1) - int(r2))
    dc = abs(int(c1) - int(c2))
    return max(dr, (dc + dr) // 2)


def ring_offsets(k: int) -> list[tuple[int, int]]:
    offs = []
    for dr in range(-k, k + 1):
        for dc in range(-2 * k, 2 * k + 1):
            if visium_hex_distance(0, 0, dr, dc) == k:
                offs.append((dr, dc))
    return offs


def _assert_hex_geometry() -> None:
    assert visium_hex_distance(0, 0, 0, 2) == 1
    assert visium_hex_distance(0, 0, 1, 1) == 1
    assert visium_hex_distance(0, 0, 0, 4) == 2
    assert visium_hex_distance(0, 0, 2, 0) == 2
    assert visium_hex_distance(0, 0, 0, 6) == 3
    assert visium_hex_distance(0, 0, 3, 3) == 3
    assert len(ring_offsets(1)) == 6
    assert len(ring_offsets(2)) == 12
    assert len(ring_offsets(3)) == 18


def spearman(x, y):
    r, p = stats.spearmanr(x, y)
    return float(r), float(p)


def partial_spearman(x, y, covars):
    def _resid(v):
        ranks = stats.rankdata(v).astype(float)
        C = np.column_stack([np.ones(len(ranks))] + [stats.rankdata(c) for c in covars])
        beta, *_ = np.linalg.lstsq(C, ranks, rcond=None)
        return ranks - C @ beta

    rx, ry = _resid(x), _resid(y)
    n = len(x)
    r, _ = stats.pearsonr(rx, ry)
    k = len(covars)
    dof = n - 2 - k
    if dof <= 0 or abs(r) >= 1:
        return float(r), float("nan")
    t = r * np.sqrt(dof / (1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), dof)
    return float(r), float(p)


def stouffer(pvals, signs):
    pvals = np.asarray(pvals, dtype=float)
    signs = np.asarray(signs, dtype=float)
    ok = np.isfinite(pvals) & np.isfinite(signs) & (pvals > 0)
    pvals, signs = pvals[ok], signs[ok]
    if len(pvals) == 0:
        return float("nan"), float("nan")
    z = stats.norm.isf(pvals / 2.0) * np.sign(signs)
    zc = z.sum() / np.sqrt(len(z))
    return float(zc), float(2 * stats.norm.sf(abs(zc)))


def wilcoxon_rhos(rhos):
    rhos = np.asarray([r for r in rhos if np.isfinite(r)])
    if len(rhos) < 5:
        return float("nan"), float("nan")
    try:
        w, p = stats.wilcoxon(rhos)
        return float(w), float(p)
    except ValueError:
        return float("nan"), float("nan")


def collapse_duplicate_genes(X, names):
    names = pd.Index(names)
    if names.is_unique:
        return X, names.tolist()
    order = pd.factorize(names)[0]
    n_unique = int(order.max()) + 1
    M = sp.csr_matrix(
        (np.ones(len(order)), (np.arange(len(order)), order)),
        shape=(len(order), n_unique),
    )
    return sp.csr_matrix(X @ M), list(pd.unique(names))


def score(logX, names_idx, genes):
    cols = [names_idx[g] for g in genes if g in names_idx]
    if not cols:
        return None, []
    sub = np.asarray(logX[:, cols].todense())
    mu, sd = sub.mean(0), sub.std(0)
    sd[sd == 0] = 1.0
    return ((sub - mu) / sd).mean(1), [g for g in genes if g in names_idx]


def load_10x_h5(path: Path):
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        X = sp.csc_matrix(
            (g["data"][:], g["indices"][:], g["indptr"][:]),
            shape=g["shape"][:],
        ).T.tocsr()
        barcodes = [b.decode() if isinstance(b, bytes) else str(b) for b in g["barcodes"][:]]
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in g["features/name"][:]]
    return X, barcodes, names


def load_positions_from_tar(tar_path: Path) -> pd.DataFrame:
    with tarfile.open(tar_path) as t:
        member = None
        for m in t.getmembers():
            if Path(m.name).name in ("tissue_positions.csv", "tissue_positions_list.csv"):
                member = m
                break
        if member is None:
            raise FileNotFoundError(f"no tissue positions in {tar_path}")
        raw = t.extractfile(member).read()
    first = raw.split(b"\n", 1)[0]
    header = 0 if b"barcode" in first else None
    df = pd.read_csv(io.BytesIO(raw), header=header)
    df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][: df.shape[1]]
    return df


def neighborhood_means(rows, cols, values, ring: int) -> np.ndarray:
    """Mean of `values` over Visium hex-ring neighbors (self excluded)."""
    key_to_i = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(rows, cols))}
    offs = ring_offsets(ring)
    out = np.full(len(values), np.nan)
    counts = np.zeros(len(values), dtype=int)
    for i, (r, c) in enumerate(zip(rows, cols)):
        acc = []
        for dr, dc in offs:
            j = key_to_i.get((int(r) + dr, int(c) + dc))
            if j is not None:
                acc.append(values[j])
        if len(acc) >= MIN_NEIGHBORS:
            out[i] = float(np.mean(acc))
            counts[i] = len(acc)
    return out, counts


def k6_neighborhood(pxl_row, pxl_col, values) -> np.ndarray:
    """PR67-like: mean of <=6 nearest pixel neighbors, distance-capped."""
    from scipy.spatial import cKDTree

    xy = np.column_stack([pxl_row, pxl_col]).astype(float)
    tree = cKDTree(xy)
    d, idx = tree.query(xy, k=7)
    ring_d = np.median(d[:, 1])
    out = np.full(len(values), np.nan)
    for i in range(len(values)):
        js = [j for j, dist in zip(idx[i, 1:], d[i, 1:]) if dist <= 2.0 * ring_d]
        if len(js) >= 3:
            out[i] = float(values[js].mean())
    return out


def analyze_visium_section(X, barcodes, names, pos, section: str) -> tuple[list[dict], dict]:
    pos = pos[pos.in_tissue == 1].set_index("barcode")
    keep = [i for i, b in enumerate(barcodes) if b in pos.index]
    X = X[keep]
    bcs = [barcodes[i] for i in keep]
    pos = pos.loc[bcs]
    X, names = collapse_duplicate_genes(X, names)

    tot = np.asarray(X.sum(1)).ravel()
    ok = tot >= MIN_COUNTS
    X, pos, tot = X[ok], pos.iloc[ok], tot[ok]
    n_spots = int(X.shape[0])
    qc = {
        "section": section,
        "patient": section.split("_")[0],
        "n_spots_qc": n_spots,
        "median_counts": float(np.median(tot)) if n_spots else float("nan"),
    }
    if n_spots < 100:
        qc["note"] = "too few spots"
        return [], qc

    Xn = X.multiply(1e4 / tot[:, None]).tocsr()
    logX = Xn.copy()
    logX.data = np.log1p(logX.data)
    names_idx = {g: i for i, g in enumerate(names)}
    epi, epi_used = score(logX, names_idx, EPI_MARKERS)
    imm, imm_used = score(logX, names_idx, IMMUNE_MARKERS)
    teff, teff_used = score(logX, names_idx, TEFF_MARKERS)
    if epi is None or imm is None or teff is None:
        qc["note"] = "missing marker genes"
        return [], qc

    tumor_spot = epi >= np.median(epi)
    immune_rich = imm >= np.quantile(imm, IMMUNE_RICH_Q)
    index = tumor_spot & ~immune_rich
    qc.update(
        n_tumor_spots=int(tumor_spot.sum()),
        n_immune_rich=int(immune_rich.sum()),
        n_index_spots=int(index.sum()),
        n_tumor_and_immune_rich=int((tumor_spot & immune_rich).sum()),
        epi_markers="|".join(epi_used),
        immune_markers="|".join(imm_used),
        teff_markers="|".join(teff_used),
    )
    if int(index.sum()) < MIN_INDEX:
        qc["note"] = "too few index spots after immune-rich exclusion"
        return [], qc

    rows_a = pos["array_row"].to_numpy()
    cols_a = pos["array_col"].to_numpy()
    nbr = {}
    for ring in RINGS:
        nbr[("immune_broad", ring)], _ = neighborhood_means(rows_a, cols_a, imm, ring)
        nbr[("t_effector", ring)], _ = neighborhood_means(rows_a, cols_a, teff, ring)
    # PR67-like sensitivity (not primary): k=6 pixel neighbors, no ring.
    nbr[("immune_broad", "k6")] = k6_neighborhood(pos["pxl_row"], pos["pxl_col"], imm)
    nbr[("t_effector", "k6")] = k6_neighborhood(pos["pxl_row"], pos["pxl_col"], teff)

    rows = []
    for gene in TARGETS:
        if gene not in names_idx:
            continue
        expr = np.asarray(logX[:, names_idx[gene]].todense()).ravel()
        for nb_name, ring in [(a, r) for a in ("immune_broad", "t_effector") for r in list(RINGS) + ["k6"]]:
            nb = nbr[(nb_name, ring)]
            # Primary index: tumor + not immune-rich, finite neighborhood.
            mask = index & np.isfinite(nb)
            # Sensitivity: all tumor spots (immune-rich kept in the CLDN4 score).
            mask_all_tumor = tumor_spot & np.isfinite(nb)
            for mask_name, m in (("index_tumor_not_immune_rich", mask), ("tumor_including_immune_rich", mask_all_tumor)):
                if int(m.sum()) < MIN_INDEX:
                    continue
                e = expr[m]
                nbv = nb[m]
                rho, p = spearman(e, nbv)
                prho, pp = partial_spearman(e, nbv, [epi[m]])
                prho2, pp2 = partial_spearman(e, nbv, [epi[m], imm[m]])
                q1, q2 = np.quantile(nbv, [1 / 3, 2 / 3])
                lo, hi = e[nbv <= q1], e[nbv >= q2]
                try:
                    _, up = stats.mannwhitneyu(hi, lo, alternative="two-sided")
                except ValueError:
                    up = float("nan")
                rows.append(
                    dict(
                        section=section,
                        patient=section.split("_")[0],
                        gene=gene,
                        neighborhood=nb_name,
                        radius=str(ring),
                        index_definition=mask_name,
                        n_index=int(m.sum()),
                        pct_detected=round(float((e > 0).mean()), 4),
                        spearman_rho=float(rho),
                        spearman_p=float(p),
                        partial_rho_epi=float(prho),
                        partial_p_epi=float(pp) if np.isfinite(pp) else float("nan"),
                        partial_rho_epi_imm=float(prho2),
                        partial_p_epi_imm=float(pp2) if np.isfinite(pp2) else float("nan"),
                        hot_minus_cold_median=float(np.median(hi) - np.median(lo)),
                        mannwhitney_p=float(up) if np.isfinite(up) else float("nan"),
                    )
                )
    return rows, qc


def meta_from_rows(res: pd.DataFrame, dataset: str) -> pd.DataFrame:
    meta = []
    keys = ["gene", "neighborhood", "radius", "index_definition"]
    for key, g in res.groupby(keys, dropna=False):
        rec = dict(zip(keys, key))
        rec.update(dataset=dataset, n_sections=int(len(g)), n_patients=int(g["patient"].nunique()))
        for col, out in (
            ("partial_rho_epi", "partial_epi"),
            ("partial_rho_epi_imm", "partial_epi_imm"),
            ("spearman_rho", "spearman"),
        ):
            w, wp = wilcoxon_rhos(g[col])
            zc, zp = stouffer(g[col.replace("rho", "p") if "rho" in col else "spearman_p"], np.sign(g[col]))
            # map p column
            if col == "partial_rho_epi":
                zc, zp = stouffer(g["partial_p_epi"], np.sign(g[col]))
            elif col == "partial_rho_epi_imm":
                zc, zp = stouffer(g["partial_p_epi_imm"], np.sign(g[col]))
            else:
                zc, zp = stouffer(g["spearman_p"], np.sign(g[col]))
            rec[f"median_{out}_rho"] = float(g[col].median())
            rec[f"mean_{out}_rho"] = float(g[col].mean())
            rec[f"n_negative_{out}"] = int((g[col] < 0).sum())
            rec[f"n_positive_{out}"] = int((g[col] > 0).sum())
            rec[f"wilcoxon_p_{out}"] = wp
            rec[f"stouffer_z_{out}"] = zc
            rec[f"stouffer_p_{out}"] = zp
            rec[f"abs_median_{out}_rho"] = abs(float(g[col].median()))
        meta.append(rec)
    return pd.DataFrame(meta)


def run_visium(data: Path, out: Path) -> dict:
    all_rows, qcs = [], []
    missing = []
    for section in TUMOR_SECTIONS:
        h5 = data / "emtab13530" / f"{section}-filtered_feature_bc_matrix.h5"
        tar = data / "emtab13530" / f"{section}-spatial.tar"
        if not h5.exists() or not tar.exists():
            missing.append(section)
            print(f"SKIP {section}: processed h5/spatial not found", flush=True)
            continue
        X, bcs, names = load_10x_h5(h5)
        pos = load_positions_from_tar(tar)
        rows, qc = analyze_visium_section(X, bcs, names, pos, section)
        qcs.append(qc)
        all_rows.extend(rows)
        print(
            f"{section}: spots={qc.get('n_spots_qc')} tumor={qc.get('n_tumor_spots')} "
            f"index={qc.get('n_index_spots')} rows={len(rows)} {qc.get('note','')}",
            flush=True,
        )
    if not all_rows:
        return {"status": "skipped", "reason": "no tumor sections processed", "missing": missing}

    res = pd.DataFrame(all_rows)
    qc_df = pd.DataFrame(qcs)
    meta = meta_from_rows(res, "E-MTAB-13530")
    res.to_csv(out / "visium_per_section.tsv", sep="\t", index=False)
    qc_df.to_csv(out / "visium_section_qc.tsv", sep="\t", index=False)
    meta.to_csv(out / "visium_meta.tsv", sep="\t", index=False)

    # Primary slice for the honest verdict
    prim = meta[
        (meta.gene == "CLDN4")
        & (meta.neighborhood == "immune_broad")
        & (meta.index_definition == "index_tumor_not_immune_rich")
        & (meta.radius.isin(["1", "2", "3"]))
    ].copy()
    return {
        "status": "ran",
        "n_sections_attempted": len(TUMOR_SECTIONS),
        "n_sections_with_rows": int(res["section"].nunique()),
        "missing_sections": missing,
        "primary_rows": prim.to_dict(orient="records"),
    }


def parse_series_meta(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="ignore") as f:
        lines = f.readlines()
    titles = chars = geos = None
    extra = []
    for line in lines:
        if line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            geos = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            extra.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
    df = pd.DataFrame({"geo": geos, "title": titles})
    for row in extra:
        if len(row) != len(geos):
            continue
        keys, vals = [], []
        for item in row:
            if ": " in item:
                k, v = item.split(": ", 1)
            else:
                k, v = "characteristic", item
            keys.append(k)
            vals.append(v)
        key = pd.Series(keys).mode().iloc[0]
        df[key] = vals
    return df


def load_pkc_rts_map(pkc_path: Path, genes: list[str]) -> dict[str, str]:
    obj = json.loads(pkc_path.read_text())
    wanted = {g.upper() for g in genes}
    rts_to_gene = {}
    for t in obj["Targets"]:
        name = str(t.get("DisplayName") or "")
        if name.upper() in wanted:
            for p in t.get("Probes") or []:
                rts = p.get("RTS_ID")
                if rts:
                    rts_to_gene[rts] = name
    return rts_to_gene


def parse_dcc_counts(path: Path, rts_map: dict[str, str]) -> dict[str, float]:
    counts = {g: 0.0 for g in set(rts_map.values())}
    aligned = None
    opener = gzip.open if path.suffix == ".gz" else open
    in_sum = False
    with opener(path, "rt", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith("Aligned,"):
                try:
                    aligned = float(line.split(",", 1)[1])
                except ValueError:
                    pass
            if line == "<Code_Summary>":
                in_sum = True
                continue
            if line == "</Code_Summary>":
                break
            if not in_sum or "," not in line:
                continue
            rts, _, rest = line.partition(",")
            if rts in rts_map:
                try:
                    counts[rts_map[rts]] = float(rest)
                except ValueError:
                    pass
    counts["_aligned"] = aligned if aligned is not None else float("nan")
    return counts


def is_tumor_compartment(label: str) -> bool:
    s = str(label).strip().lower()
    if s in {"na", "nan", "", "ntc", "no template control"}:
        return False
    keys = ("tumor", "tumour", "panck", "pan-ck", "pancytokeratin", "ck+", "ck ", "epithelial")
    # Exact common GeoMx labels
    if s in {"ck", "panck", "pan-ck", "tumor", "tumour", "epithelial", "neoplastic"}:
        return True
    return any(k in s for k in keys) and "cd45" not in s and "stroma" not in s and "leukocyte" not in s


def teff_score_matrix(df: pd.DataFrame, genes: list[str], transform) -> tuple[np.ndarray, list[str]]:
    avail = [g for g in genes if g in df.columns]
    if not avail:
        return None, []
    Z = np.column_stack([transform(df[g].astype(float).to_numpy()) for g in avail])
    Z = (Z - Z.mean(0)) / np.where(Z.std(0) == 0, 1, Z.std(0))
    return Z.mean(1), avail


def geomx_corr_block(expr: pd.Series, score: np.ndarray, cohort: str, gene: str, vs: str, n: int) -> dict:
    rho, p = spearman(expr.to_numpy(float), score)
    return dict(cohort=cohort, gene=gene, vs=vs, n=n, spearman_rho=float(rho), spearman_p=float(p))


def run_geomx_published(xlsx: Path, out: Path) -> dict:
    xl = pd.ExcelFile(xlsx)
    sheets = set(xl.sheet_names)
    info = {"xlsx_sheets": sorted(sheets)}
    # Tumor-compartment sheets from Nat Genet 2025 source data (PR67 mapping).
    wanted = {
        "Yale": "source_data_Figure_6b",
        "Greek": "source_data_Figure_6d",
    }
    # Do not use stromal sheets 6f/6h — tumor compartment only.
    corr_rows = []
    used = {}
    for cohort, sheet in wanted.items():
        if sheet not in sheets:
            info[f"{cohort}_sheet"] = "NOT FOUND"
            continue
        df = pd.read_excel(xl, sheet)
        used[cohort] = {"sheet": sheet, "n_rows": int(len(df)), "n_cols": int(df.shape[1])}
        if cohort == "Yale":
            def transform(v):
                return np.log2(np.asarray(v, float) + 1)
        else:
            # Greek sheets are batch-corrected (can be negative) → rank inverse-normal.
            def transform(v):
                x = np.asarray(v, float)
                r = stats.rankdata(x)
                return stats.norm.ppf((r - 0.5) / len(r))

        imm, imm_genes = teff_score_matrix(df, TEFF_MARKERS, transform)
        broad, broad_genes = teff_score_matrix(df, IMMUNE_MARKERS, transform)
        used[cohort]["teff_genes"] = imm_genes
        used[cohort]["broad_genes"] = broad_genes
        for gene in TARGETS:
            if gene not in df.columns:
                continue
            x = pd.Series(transform(df[gene]))
            if imm is not None:
                corr_rows.append(geomx_corr_block(x, imm, f"{cohort} tumor compartment", gene, "tumor_AOI_Teff_score", len(df)))
            if broad is not None:
                corr_rows.append(geomx_corr_block(x, broad, f"{cohort} tumor compartment", gene, "tumor_AOI_immune_broad_score", len(df)))
            for ig in ["PTPRC", "CD8A", "CD3E", "NKG7", "GZMB"]:
                if ig in df.columns:
                    rho, p = spearman(x, pd.Series(transform(df[ig])))
                    corr_rows.append(dict(
                        cohort=f"{cohort} tumor compartment", gene=gene, vs=f"tumor_AOI_{ig}",
                        n=int(len(df)), spearman_rho=float(rho), spearman_p=float(p),
                    ))
    corr = pd.DataFrame(corr_rows)
    if not corr.empty:
        corr.to_csv(out / "geomx_tumor_compartment_corr.tsv", sep="\t", index=False)
    info["used"] = used
    info["n_corr_rows"] = int(len(corr))
    info["status"] = "ran" if not corr.empty else "skipped"
    if corr.empty:
        info["reason"] = "published tumor-compartment sheets missing or lacked target genes"
    return info, corr


def run_geomx_dcc(data: Path, out: Path) -> dict:
    pkc = data / "geomx_pkc" / "Hs_R_NGS_WTA_v1.0.pkc"
    tar_path = data / "gse271689" / "GSE271689_RAW.tar"
    matrix = data / "gse271689" / "GSE271689_series_matrix.txt.gz"
    if not (pkc.exists() and tar_path.exists() and matrix.exists()):
        return {"status": "skipped", "reason": "DCC/PKC/series matrix not all present"}, pd.DataFrame()

    genes = list(dict.fromkeys(TARGETS + IMMUNE_MARKERS + TEFF_MARKERS))
    rts_map = load_pkc_rts_map(pkc, genes)
    rows = []
    with tarfile.open(tar_path, "r") as tar:
        for m in tar.getmembers():
            if not m.name.endswith(".dcc.gz"):
                continue
            f = tar.extractfile(m)
            tmp = Path("/tmp") / Path(m.name).name
            tmp.write_bytes(f.read())
            counts = parse_dcc_counts(tmp, rts_map)
            gsm = Path(m.name).name.split("_")[0]
            rec = {"geo": gsm, "dcc": Path(m.name).name}
            rec.update(counts)
            rows.append(rec)
            tmp.unlink(missing_ok=True)
    expr = pd.DataFrame(rows)
    meta = parse_series_meta(matrix)
    df = expr.merge(meta, on="geo", how="left")
    if "spotid" in df.columns:
        df = df[~df["spotid"].astype(str).str.contains("No Template", case=False, na=False)]
    cell_col = None
    for c in df.columns:
        if c.lower() in {"cell type", "cell_type", "segment", "aoi", "compartment"}:
            cell_col = c
            break
    labels = sorted(df[cell_col].dropna().astype(str).unique().tolist()) if cell_col else []
    if cell_col:
        tumor = df[df[cell_col].map(is_tumor_compartment)].copy()
    else:
        tumor = df.copy()
    aligned = pd.to_numeric(tumor["_aligned"], errors="coerce")
    for g in [c for c in tumor.columns if c in genes]:
        raw_c = pd.to_numeric(tumor[g], errors="coerce")
        tumor[f"{g}_cpm"] = np.where(aligned > 0, raw_c / aligned * 1e6, np.nan)
    tumor.to_csv(out / "geomx_tumor_aoi_table.tsv", sep="\t", index=False)

    def zmean(frame, markers, suffix="_cpm"):
        cols = [f"{g}{suffix}" for g in markers if f"{g}{suffix}" in frame.columns]
        if not cols:
            return None, []
        Z = frame[cols].apply(pd.to_numeric, errors="coerce")
        Z = (Z - Z.mean()) / Z.std().replace(0, np.nan)
        return Z.mean(axis=1).to_numpy(), cols

    imm, imm_cols = zmean(tumor, TEFF_MARKERS)
    broad, broad_cols = zmean(tumor, IMMUNE_MARKERS)
    corr_rows = []
    for gene in TARGETS:
        col = f"{gene}_cpm"
        if col not in tumor.columns:
            continue
        x = pd.to_numeric(tumor[col], errors="coerce")
        if imm is not None:
            pair = pd.DataFrame({"x": x, "y": imm}).dropna()
            if len(pair) >= 8:
                rho, p = spearman(pair.x, pair.y)
                corr_rows.append(dict(
                    source="GEO_DCC_CPM", cohort="GSE271689 tumor AOI",
                    gene=gene, vs="tumor_AOI_Teff_score", n=int(len(pair)),
                    spearman_rho=float(rho), spearman_p=float(p),
                ))
        if broad is not None:
            pair = pd.DataFrame({"x": x, "y": broad}).dropna()
            if len(pair) >= 8:
                rho, p = spearman(pair.x, pair.y)
                corr_rows.append(dict(
                    source="GEO_DCC_CPM", cohort="GSE271689 tumor AOI",
                    gene=gene, vs="tumor_AOI_immune_broad_score", n=int(len(pair)),
                    spearman_rho=float(rho), spearman_p=float(p),
                ))
    corr = pd.DataFrame(corr_rows)
    if not corr.empty:
        corr.to_csv(out / "geomx_dcc_tumor_aoi_corr.tsv", sep="\t", index=False)
    return {
        "status": "ran" if not corr.empty else "skipped",
        "n_dcc_rows": int(len(expr)),
        "n_tumor_aoi": int(len(tumor)),
        "cell_type_column": cell_col,
        "cell_type_labels": labels,
        "tumor_labels_used": sorted(tumor[cell_col].astype(str).unique().tolist()) if cell_col and len(tumor) else [],
        "teff_cpm_cols": imm_cols,
        "broad_cpm_cols": broad_cols,
        "n_corr_rows": int(len(corr)),
    }, corr


def make_figures(out: Path, visium_res: pd.DataFrame, visium_meta: pd.DataFrame, geomx_pub: pd.DataFrame) -> list[str]:
    figs = []
    figdir = out / "figures"
    figdir.mkdir(exist_ok=True)

    sub = visium_res[
        (visium_res.gene == "CLDN4")
        & (visium_res.neighborhood == "immune_broad")
        & (visium_res.index_definition == "index_tumor_not_immune_rich")
        & (visium_res.radius.isin(["1", "2", "3"]))
    ].copy()
    if not sub.empty:
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        rings = ["1", "2", "3"]
        data = [sub.loc[sub.radius == r, "partial_rho_epi"].to_numpy() for r in rings]
        ax.axhline(0, color="k", lw=0.8)
        ax.boxplot(data, labels=[f"ring {r}" for r in rings], widths=0.55)
        rng = np.random.default_rng(0)
        for i, y in enumerate(data, start=1):
            x = i + (rng.random(len(y)) - 0.5) * 0.18
            ax.scatter(x, y, s=18, c="#1d4e89", alpha=0.75, zorder=3)
        ax.set_ylabel("Per-section partial Spearman ρ\n(CLDN4 vs broad-immune neighborhood)")
        ax.set_title("E-MTAB-13530 tumor spots, immune-rich excluded from CLDN4")
        fig.tight_layout()
        p = figdir / "visium_cldn4_rho_by_ring.png"
        fig.savefig(p, dpi=140)
        fig.savefig(figdir / "visium_cldn4_rho_by_ring.pdf")
        plt.close(fig)
        figs.append(str(p.name))

    if geomx_pub is not None and not geomx_pub.empty:
        g = geomx_pub[geomx_pub.gene == "CLDN4"].copy()
        if not g.empty:
            fig, ax = plt.subplots(figsize=(7.6, 4.2))
            y = np.arange(len(g))
            ax.axvline(0, color="k", lw=0.8)
            ax.scatter(g.spearman_rho, y, s=36, c="#9b2226")
            ax.set_yticks(y)
            ax.set_yticklabels([f"{a} · {b}" for a, b in zip(g.cohort, g.vs)], fontsize=8)
            ax.set_xlabel("Spearman ρ (tumor compartment only)")
            ax.set_title("GSE271689 GeoMx: CLDN4 vs immune scores in tumor AOIs")
            fig.tight_layout()
            p = figdir / "geomx_tumor_cldn4_vs_immune.png"
            fig.savefig(p, dpi=140)
            fig.savefig(figdir / "geomx_tumor_cldn4_vs_immune.pdf")
            plt.close(fig)
            figs.append(str(p.name))
    return figs


def verdict_from_primary(prim: list[dict]) -> dict:
    if not prim:
        return {
            "label": "NOT_TESTABLE",
            "statement": "Primary Visium CLDN4 ring analysis produced no section-level rows.",
        }
    rhos = [r["median_partial_epi_rho"] for r in prim]
    absmax = max(abs(x) for x in rhos)
    # User claim is spatial exclusion of meaningful size. PR67 mismatch was |ρ|≤0.06.
    # We do not treat |ρ|≤0.06 as support. We also do not invent a "match" threshold
    # for a qualitative claim; we report size and direction honestly.
    all_tiny = absmax <= 0.06
    any_neg = any(x < 0 for x in rhos)
    label = "DOES_NOT_SUPPORT_CLAIM"
    if all_tiny:
        statement = (
            f"Reworked Visium tumor-only CLDN4 vs immune-neighborhood partial ρ "
            f"at rings 1/2/3 has |median ρ| ≤ {absmax:.3f} (still ≤ 0.06). "
            "The user spatial-exclusion claim is not supported. "
            "This does not rescue the PR67 mismatch."
        )
    else:
        # Larger than PR67's |ρ|≤0.06, but still judge support by sign + size.
        med = {r["radius"]: r["median_partial_epi_rho"] for r in prim}
        statement = (
            f"Reworked median partial ρ by ring: {med}. "
            "We do not treat this as confirmation of the user claim unless the "
            "effect is a consistent negative exclusion of non-trivial size; "
            "see WRITEUP for the actual numbers."
        )
        if all(x < -0.15 for x in rhos):
            label = "DIRECTIONALLY_CONSISTENT_NOT_A_USER_NUMBER"
        else:
            label = "DOES_NOT_SUPPORT_CLAIM"
    return {
        "label": label,
        "abs_max_median_partial_rho_rings": absmax,
        "median_partial_rho_by_ring": {r["radius"]: r["median_partial_epi_rho"] for r in prim},
        "n_sections_by_ring": {r["radius"]: r["n_sections"] for r in prim},
        "wilcoxon_p_by_ring": {r["radius"]: r["wilcoxon_p_partial_epi"] for r in prim},
        "all_rings_abs_le_0.06": all_tiny,
        "any_ring_negative": any_neg,
        "did_we_tune_thresholds": False,
        "statement": statement,
    }


def main() -> int:
    _assert_hex_geometry()
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/tmp/b6_radius_data")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()
    data = Path(args.data_dir)
    out = Path(args.out_dir) if args.out_dir else Path(__file__).resolve().parent
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)

    visium_info = run_visium(data, out)
    xlsx = data / "gse271689" / "41588_2025_2351_MOESM10_ESM.xlsx"
    if xlsx.exists():
        geomx_pub_info, geomx_pub = run_geomx_published(xlsx, out)
    else:
        geomx_pub_info, geomx_pub = {"status": "skipped", "reason": "MOESM10 xlsx not present"}, pd.DataFrame()
    geomx_dcc_info, _ = run_geomx_dcc(data, out)

    visium_res = pd.read_csv(out / "visium_per_section.tsv", sep="\t") if (out / "visium_per_section.tsv").exists() else pd.DataFrame()
    visium_meta = pd.read_csv(out / "visium_meta.tsv", sep="\t") if (out / "visium_meta.tsv").exists() else pd.DataFrame()
    figs = make_figures(out, visium_res, visium_meta, geomx_pub)

    honest = verdict_from_primary(visium_info.get("primary_rows") or [])
    summary = {
        "task": "B6_radius",
        "claim": USER_CLAIM,
        "pr67_mismatch": PR67_MISMATCH,
        "rework": {
            "visium": "tumor sections + tumor spots only; hex rings 1/2/3; immune-rich spots excluded from CLDN4/TACSTD2 score",
            "geomx": "tumor / PanCK compartment only; same-AOI correlation (no rings, no stroma pairing)",
            "did_we_tune_to_enlarge_rho": False,
        },
        "visium": visium_info,
        "geomx_published_tumor": geomx_pub_info,
        "geomx_dcc_tumor": geomx_dcc_info,
        "figures": figs,
        "honest_verdict": honest,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps({"honest_verdict": honest, "visium_status": visium_info.get("status"), "geomx_pub": geomx_pub_info.get("status"), "geomx_dcc": geomx_dcc_info.get("status")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
