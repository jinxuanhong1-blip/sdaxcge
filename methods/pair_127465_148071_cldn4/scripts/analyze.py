#!/usr/bin/env python3
"""Additive CLDN4-only pairwise combo: GSE127465 + GSE148071.

Patient-level malignant CLDN4 vs T/NK on each public processed matrix.
If the two cohort associations differ, CellChat-style outgoing
(Mal CLDN4-high vs low → T/NK). Not a triple/quad. No dual-high.
Not a cell-level merge (inDrops normalized MTX vs 10x UMI).

CellChat R is not run. Probability is Jin et al. 2021 Hill / mass-action
on truncated means (CellChatDB v2 protein pairs) plus a high/low label
permutation among malignant cells.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import tarfile
import urllib.request
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm, spearmanr

ROOT = Path(__file__).resolve().parents[1]
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
NBOOT = 100
SEED = 1
MIN_GROUP = 25
MIN_EPI = 20
MIN_TNK = 20
MIN_CELLCHAT_ARM = 25

GSE127465 = {
    "acc": "GSE127465",
    "pmid": "30979687",
    "title": "Zilionis et al. Immunity 2019; NSCLC inDrops (human tumor + blood)",
    "base": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/",
    "mtx": "GSE127465_human_counts_normalized_54773x41861.mtx.gz",
    "meta": "GSE127465_human_cell_metadata_54773x25.tsv.gz",
    "genes": "GSE127465_gene_names_human_41861.tsv.gz",
    "bytes_mtx": 528303938,
}
GSE148071 = {
    "acc": "GSE148071",
    "pmid": "33953163",
    "title": "Wu et al. Nat Commun 2021; 42 advanced NSCLC biopsies",
    "tar_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/GSE148071_RAW.tar",
    "tar": "GSE148071_RAW.tar",
    "bytes_tar": 180193280,
}

LINEAGE_MARKERS = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2"],
    "Mast": ["TPSAB1", "CPA3"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}
EXTRA = ["TACSTD2", "CLDN4", "PTPRC", "CD8A", "CD4", "NCAM1", "IFNG", "TNF"]

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

PAT_SPEC_RE = re.compile(r"^Patient(\d+)-specific$")


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
        lig = parse_symbols(getattr(rec, "ligand_symbol", None))
        recp = parse_symbols(getattr(rec, "receptor_symbol", None))
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


def wanted_genes(lr: pd.DataFrame) -> set[str]:
    genes = set(EXTRA)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    for rec in lr.itertuples(index=False):
        genes.update(rec.ligand_genes)
        genes.update(rec.receptor_genes)
    return genes


def download(url: str, dest: Path, expect_bytes: int | None = None) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        if expect_bytes is None or dest.stat().st_size == expect_bytes:
            print(f"exists {dest} ({dest.stat().st_size} B)", flush=True)
            return dest
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, dest)
    print(f"wrote {dest} ({dest.stat().st_size} B)", flush=True)
    return dest


def spearman_safe(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return float("nan"), float("nan"), n
    rho, p = spearmanr(x[m], y[m])
    return float(rho), float(p), n


def fisher_z(rho: float) -> float:
    r = min(0.999999, max(-0.999999, float(rho)))
    return float(np.arctanh(r))


def fisher_z_diff(rho1: float, n1: int, rho2: float, n2: int) -> dict:
    if n1 < 4 or n2 < 4 or not np.isfinite(rho1) or not np.isfinite(rho2):
        return {"z": float("nan"), "p": float("nan"), "delta_z": float("nan")}
    z1, z2 = fisher_z(rho1), fisher_z(rho2)
    se = math.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))
    z = (z1 - z2) / se
    p = float(2.0 * (1.0 - norm.cdf(abs(z))))
    return {"z": float(z), "p": p, "delta_z": float(z1 - z2), "se": se}


def rho_ci(rho: float, n: int) -> tuple[float, float]:
    if n < 4 or not np.isfinite(rho):
        return float("nan"), float("nan")
    se = math.sqrt(1.0 / (n - 3))
    z = fisher_z(rho)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def associations_differ(rho1: float, n1: int, rho2: float, n2: int) -> tuple[bool, str]:
    """True if the two patient-level CLDN4 vs T/NK associations differ."""
    reasons = []
    if not np.isfinite(rho1) or not np.isfinite(rho2):
        return False, "one_or_both_rho_undefined"
    if (rho1 > 0) != (rho2 > 0) and abs(rho1) >= 0.10 and abs(rho2) >= 0.10:
        reasons.append("opposite_sign")
    elif (rho1 > 0) != (rho2 > 0) and (abs(rho1) >= 0.20 or abs(rho2) >= 0.20):
        reasons.append("opposite_sign_one_weak")
    fz = fisher_z_diff(rho1, n1, rho2, n2)
    if np.isfinite(fz["p"]) and fz["p"] < 0.05:
        reasons.append(f"fisher_z_p={fz['p']:.4g}")
    # practical: one |ρ|≥0.30 and the other |ρ|<0.15
    if max(abs(rho1), abs(rho2)) >= 0.30 and min(abs(rho1), abs(rho2)) < 0.15:
        reasons.append("one_moderate_one_near_zero")
    return (len(reasons) > 0), ";".join(reasons) if reasons else "similar"


# ---------------------------------------------------------------------------
# GSE127465
# ---------------------------------------------------------------------------

def load_127465_meta(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if len(df) != 54773:
        raise SystemExit(f"GSE127465 metadata rows {len(df)} != 54773")
    return df


def load_127465_genes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        genes = [ln.strip().split("\t")[0] for ln in handle if ln.strip()]
    if len(genes) != 41861:
        raise SystemExit(f"GSE127465 genes {len(genes)} != 41861")
    return genes


def stream_mtx_columns(mtx_path: Path, wanted_cols_1based: dict[int, str], n_rows: int) -> dict[str, np.ndarray]:
    """Matrix Market coordinate: 54773 cells x 41861 genes, 1-based."""
    out = {name: np.zeros(n_rows, dtype=np.float32) for name in wanted_cols_1based.values()}
    wanted = set(wanted_cols_1based)
    n_hit = 0
    with gzip.open(mtx_path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            # first non-comment is the banner
            banner = line.split()
            if len(banner) == 3 and banner[0].isdigit():
                nrows, ncols, nnz = map(int, banner)
                if nrows != n_rows:
                    raise SystemExit(f"MTX rows {nrows} != {n_rows}")
                break
        for i, line in enumerate(handle, 1):
            parts = line.split()
            if len(parts) < 3:
                continue
            r = int(parts[0])
            c = int(parts[1])
            if c in wanted:
                out[wanted_cols_1based[c]][r - 1] = float(parts[2])
                n_hit += 1
            if i % 10_000_000 == 0:
                print(f"  MTX lines {i:,}", flush=True)
    print(f"  MTX wanted nnz={n_hit}", flush=True)
    return out


def patient_num(patient: str) -> str:
    return str(int(str(patient).lower().lstrip("p")))


def score_127465_patients(meta: pd.DataFrame, expr: dict[str, np.ndarray]) -> pd.DataFrame:
    tumor = meta["Tissue"].astype(str) == "tumor"
    maj = meta["Major cell type"].astype(str)
    pat = meta["Patient"].astype(str)
    cldn4 = expr["CLDN4"]
    rows = []
    for p in sorted(pat.unique(), key=lambda x: int(patient_num(x))):
        own = f"Patient{patient_num(p)}-specific"
        in_p = tumor & (pat == p)
        mal = in_p & (maj == own)
        mal_any_spec = in_p & maj.str.contains("specific", case=False, regex=False)
        tnk = in_p & maj.isin(["tT cells", "tNK cells"])
        t_only = in_p & (maj == "tT cells")
        nk_only = in_p & (maj == "tNK cells")
        n_tumor = int(in_p.sum())
        n_mal = int(mal.sum())
        n_tnk = int(tnk.sum())
        eligible = n_mal >= MIN_EPI and n_tnk >= MIN_TNK
        cldn4_mean = float(np.log1p(cldn4[mal.to_numpy()]).mean()) if n_mal else float("nan")
        cldn4_pct = float((cldn4[mal.to_numpy()] > 0).mean()) if n_mal else float("nan")
        rows.append(
            {
                "cohort": "GSE127465",
                "patient": p,
                "n_tumor": n_tumor,
                "n_malignant": n_mal,
                "n_malignant_any_patient_specific": int(mal_any_spec.sum()),
                "n_T": int(t_only.sum()),
                "n_NK": int(nk_only.sum()),
                "n_TNK": n_tnk,
                "frac_TNK": (n_tnk / n_tumor) if n_tumor else float("nan"),
                "CLDN4_mean_log1p": cldn4_mean,
                "CLDN4_pct_pos": cldn4_pct,
                "eligible": eligible,
                "malignant_def": own,
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# GSE148071
# ---------------------------------------------------------------------------

def extract_148071(tar_path: Path, dest: Path) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    files = sorted(dest.rglob("*_exp.txt.gz"))
    if len(files) >= 42:
        return files
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    files = sorted(dest.rglob("*_exp.txt.gz"))
    if len(files) < 42:
        raise SystemExit(f"GSE148071 expected 42 exp matrices, got {len(files)}")
    return files


def patient_from_filename(path: Path) -> str:
    m = re.search(r"_(P\d+)_", path.name)
    return m.group(1) if m else path.stem


def stream_148071_one(path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        header = [h.strip().strip('"') for h in header if h != ""]
        if header and header[0] in {"", "gene", "Gene", "index", "Index", "GENE"}:
            cell_ids = header[1:]
        else:
            cell_ids = header
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.strip().strip('"').split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                toks = line.rstrip("\n").split("\t")
                arr = np.asarray(toks[1:], dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{path.name} {gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
    return cell_ids, found, n_umi


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores: dict[str, np.ndarray], cd3: np.ndarray) -> np.ndarray:
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    labels = np.array(names, dtype=object)[best].copy()
    t_idx = names.index("T")
    nk_idx = names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both_high = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk_best = np.isin(labels, ["T", "NK"])
    labels[close & both_high & tnk_best & (cd3 > 0.15)] = "T"
    labels[close & both_high & tnk_best & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def load_148071(files: list[Path], wanted: set[str]) -> tuple[pd.DataFrame, dict[str, np.ndarray], np.ndarray]:
    """Return cell table, expr dict (raw UMI for wanted genes), n_umi."""
    parts = []
    expr_blocks: dict[str, list[np.ndarray]] = {g: [] for g in wanted}
    umi_blocks = []
    for path in files:
        pid = patient_from_filename(path)
        print(f"  {path.name} {pid}", flush=True)
        cell_ids, found, n_umi = stream_148071_one(path, wanted)
        n = len(cell_ids)
        scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
        log_cp = {g: np.log1p(found[g] * scale).astype(np.float32) for g in found}
        scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
        cd3 = module_score(log_cp, ["CD3D", "CD3E", "CD3G"], n)
        lineage = assign_lineage(scores, cd3)
        cldn4 = log_cp["CLDN4"] if "CLDN4" in log_cp else np.zeros(n, dtype=np.float32)
        parts.append(
            pd.DataFrame(
                {
                    "patient": pid,
                    "cell": cell_ids,
                    "lineage": lineage,
                    "CLDN4_log1p_cp10k": cldn4,
                    "n_umi": n_umi,
                }
            )
        )
        for g in wanted:
            expr_blocks[g].append(found[g] if g in found else np.zeros(n, dtype=np.float32))
        umi_blocks.append(n_umi.astype(np.float32))
    cells = pd.concat(parts, ignore_index=True)
    expr = {g: np.concatenate(blocks) for g, blocks in expr_blocks.items()}
    n_umi = np.concatenate(umi_blocks)
    return cells, expr, n_umi


def score_148071_patients(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for p, sub in cells.groupby("patient", sort=True):
        mal = sub["lineage"] == "Epithelial"
        tnk = sub["lineage"].isin(["T", "NK"])
        n_mal = int(mal.sum())
        n_tnk = int(tnk.sum())
        n = int(len(sub))
        rows.append(
            {
                "cohort": "GSE148071",
                "patient": p,
                "n_tumor": n,
                "n_malignant": n_mal,
                "n_T": int((sub["lineage"] == "T").sum()),
                "n_NK": int((sub["lineage"] == "NK").sum()),
                "n_TNK": n_tnk,
                "frac_TNK": n_tnk / n if n else float("nan"),
                "CLDN4_mean_log1p": float(sub.loc[mal, "CLDN4_log1p_cp10k"].mean()) if n_mal else float("nan"),
                "CLDN4_pct_pos": float((sub.loc[mal, "CLDN4_log1p_cp10k"] > 0).mean()) if n_mal else float("nan"),
                "eligible": n_mal >= MIN_EPI and n_tnk >= MIN_TNK,
                "malignant_def": "marker_argmax_epithelial_putative",
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# CellChat-style outgoing
# ---------------------------------------------------------------------------

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


def trim_mean_fast(sub: np.ndarray, proportiontocut: float = TRIM) -> np.ndarray:
    n = sub.shape[1]
    if n == 0:
        return np.zeros(sub.shape[0], dtype=np.float64)
    if n == 1:
        return sub[:, 0].astype(np.float64)
    k = int(n * proportiontocut)
    if k == 0:
        return sub.mean(axis=1)
    s = np.sort(sub, axis=1)
    return s[:, k : n - k].mean(axis=1)


def group_trim_means(expr: np.ndarray, pos: np.ndarray, labels: np.ndarray, groups: list[str]):
    means = np.zeros((expr.shape[0], len(groups)), dtype=np.float64)
    props = np.zeros((expr.shape[0], len(groups)), dtype=np.float64)
    counts = {}
    for j, g in enumerate(groups):
        idx = np.flatnonzero(labels == g)
        counts[g] = int(idx.size)
        if idx.size == 0:
            continue
        sub = expr[:, idx]
        means[:, j] = trim_mean_fast(sub)
        props[:, j] = pos[:, idx].mean(axis=1)
    return means, props, counts


def pair_specs(lr: pd.DataFrame, gene_index: dict[str, int]):
    specs = []
    for rec in lr.itertuples(index=False):
        specs.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": rec.ligand_genes,
                "receptor_genes": rec.receptor_genes,
                "lig_ix": [gene_index[g] for g in rec.ligand_genes],
                "rec_ix": [gene_index[g] for g in rec.receptor_genes],
            }
        )
    return specs


def score_pairs(specs, gene_means, gene_props, groups, counts, src_tgt):
    gpos = {g: i for i, g in enumerate(groups)}
    cache = {}

    def cached(ix_key, ix):
        hit = cache.get(ix_key)
        if hit is None:
            mu = geom_mean_rows(gene_means[ix])
            pr = gene_props[ix].min(axis=0) if len(ix) > 1 else gene_props[ix[0]]
            hit = (mu, pr)
            cache[ix_key] = hit
        return hit

    rows = []
    probs = []
    for spec in specs:
        lig_mu, lig_pr = cached(("L", tuple(spec["lig_ix"])), spec["lig_ix"])
        rec_mu, rec_pr = cached(("R", tuple(spec["rec_ix"])), spec["rec_ix"])
        for src, tgt in src_tgt:
            if src not in gpos or tgt not in gpos:
                continue
            i, j = gpos[src], gpos[tgt]
            if counts.get(src, 0) < MIN_GROUP or counts.get(tgt, 0) < MIN_GROUP:
                continue
            detected = (float(lig_pr[i]) >= EXPR_PROP) and (float(rec_pr[j]) >= EXPR_PROP)
            prob = hill_prob(float(lig_mu[i]), float(rec_mu[j])) if detected else 0.0
            rows.append(
                {
                    "interaction_name": spec["interaction_name"],
                    "pathway_name": spec["pathway_name"],
                    "annotation": spec["annotation"],
                    "ligand": spec["ligand"],
                    "receptor": spec["receptor"],
                    "ligand_genes": "|".join(spec["ligand_genes"]),
                    "receptor_genes": "|".join(spec["receptor_genes"]),
                    "source": src,
                    "target": tgt,
                    "n_source": counts[src],
                    "n_target": counts[tgt],
                    "prob": prob,
                    "detected": bool(detected),
                    "ligand_prop": float(lig_pr[i]),
                    "receptor_prop": float(rec_pr[j]),
                }
            )
            probs.append(prob)
    return pd.DataFrame(rows), np.asarray(probs, dtype=np.float64)


def permute_outgoing(expr, pos, specs, labels, groups, src_tgt, mal_mask, nboot, rng):
    obs_means, obs_props, counts = group_trim_means(expr, pos, labels, groups)
    obs, obs_prob = score_pairs(specs, obs_means, obs_props, groups, counts, src_tgt)
    ge = np.zeros(len(obs), dtype=np.int32)
    mal_idx = np.flatnonzero(mal_mask)
    mal_labs = labels[mal_idx].copy()
    work = labels.copy()
    for b in range(nboot):
        work[mal_idx] = rng.permutation(mal_labs)
        means, props, counts_b = group_trim_means(expr, pos, work, groups)
        _, perm_prob = score_pairs(specs, means, props, groups, counts_b, src_tgt)
        ge += (perm_prob >= obs_prob).astype(np.int32)
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
    return "|".join(tags) if tags else "other"


def contrast_outgoing(df: pd.DataFrame) -> pd.DataFrame:
    a = df[(df.source == "Mal_high") & (df.target == "TNK")].set_index("interaction_name")
    b = df[(df.source == "Mal_low") & (df.target == "TNK")].set_index("interaction_name")
    common = a.index.intersection(b.index)
    if len(common) == 0:
        return pd.DataFrame()
    out = a.loc[common, ["pathway_name", "annotation", "ligand", "receptor", "ligand_genes", "receptor_genes"]].copy()
    out["prob_high"] = a.loc[common, "prob"].to_numpy()
    out["prob_low"] = b.loc[common, "prob"].to_numpy()
    out["pval_high"] = a.loc[common, "pval"].to_numpy()
    out["pval_low"] = b.loc[common, "pval"].to_numpy()
    out["delta_prob"] = out["prob_high"] - out["prob_low"]
    out["ligand_class"] = [ligand_class(x) for x in out["ligand_genes"]]
    out["sig_high"] = (out["pval_high"] < 0.05) & (out["prob_high"] > 0)
    out["sig_low"] = (out["pval_low"] < 0.05) & (out["prob_low"] > 0)
    out["sig_diff"] = (out["sig_high"] | out["sig_low"]) & (out["delta_prob"] != 0)
    return out.reset_index().sort_values("delta_prob", ascending=False)


def run_cellchat_outgoing(
    name: str,
    expr_cp: dict[str, np.ndarray],
    keep: np.ndarray,
    mal_mask: np.ndarray,
    tnk_mask: np.ndarray,
    cldn4: np.ndarray,
    lr: pd.DataFrame,
    out_prefix: Path,
) -> pd.DataFrame:
    """keep/mal/tnk are boolean over the full expr arrays. Split malignant by median CLDN4."""
    idx = np.flatnonzero(keep)
    if idx.size == 0:
        return pd.DataFrame()
    genes = sorted({g for rec in lr.itertuples(index=False) for g in rec.ligand_genes + rec.receptor_genes})
    genes = [g for g in genes if g in expr_cp]
    gene_index = {g: i for i, g in enumerate(genes)}
    mat = np.vstack([expr_cp[g][idx] for g in genes]).astype(np.float32)
    pos = (mat > 0).astype(np.float32)
    specs = pair_specs(lr, gene_index)
    labels = np.full(idx.size, "drop", dtype=object)
    mal_l = mal_mask[idx]
    tnk_l = tnk_mask[idx]
    labels[tnk_l] = "TNK"
    vals = cldn4[idx][mal_l]
    if vals.size < MIN_GROUP * 2:
        print(f"  {name}: not enough malignant cells for median split ({vals.size})", flush=True)
        return pd.DataFrame()
    med = float(np.median(vals))
    labels[mal_l & (cldn4[idx] < med)] = "Mal_low"
    labels[mal_l & (cldn4[idx] >= med)] = "Mal_high"
    n_high = int((labels == "Mal_high").sum())
    n_low = int((labels == "Mal_low").sum())
    n_tnk = int((labels == "TNK").sum())
    print(f"  {name} CellChat-style outgoing: Mal_high={n_high} Mal_low={n_low} TNK={n_tnk} cut={med:.4g}", flush=True)
    if min(n_high, n_low, n_tnk) < MIN_GROUP:
        print(f"  {name}: skip CellChat (group < {MIN_GROUP})", flush=True)
        return pd.DataFrame()
    groups = ["Mal_high", "Mal_low", "TNK"]
    src_tgt = [("Mal_high", "TNK"), ("Mal_low", "TNK")]
    mal_split = np.isin(labels, ["Mal_high", "Mal_low"])
    rng = np.random.default_rng(SEED)
    obs = permute_outgoing(mat, pos, specs, labels, groups, src_tgt, mal_split, NBOOT, rng)
    obs.to_csv(out_prefix.parent / f"{out_prefix.name}_pairs.tsv", sep="\t", index=False)
    contr = contrast_outgoing(obs)
    contr.insert(0, "cohort", name)
    contr.to_csv(out_prefix.parent / f"{out_prefix.name}_outgoing.tsv", sep="\t", index=False)
    meta = {
        "cohort": name,
        "n_mal_high": n_high,
        "n_mal_low": n_low,
        "n_tnk": n_tnk,
        "cldn4_median_cut": med,
        "n_pairs_db": int(len(lr)),
        "n_detected_high_or_low": int(((contr.prob_high > 0) | (contr.prob_low > 0)).sum()) if len(contr) else 0,
        "n_sig_diff": int(contr["sig_diff"].sum()) if len(contr) else 0,
    }
    (out_prefix.parent / f"{out_prefix.name}_summary.json").write_text(json.dumps(meta, indent=2))
    return contr


def forest_plot(rows: list[dict], path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 2.6))
    y = np.arange(len(rows))
    for i, r in enumerate(rows):
        lo, hi = r["ci95_lo"], r["ci95_hi"]
        color = "#1f4e79" if r.get("pooled") else "#4a4a4a"
        marker = "D" if r.get("pooled") else "o"
        ax.plot([lo, hi], [y[i], y[i]], color=color, lw=1.6)
        ax.plot(r["rho"], y[i], marker=marker, color=color, ms=7)
        ax.text(1.02, y[i], f"ρ={r['rho']:+.2f}  p={r['p']:.3g}  n={r['n']}", va="center", fontsize=8, family="monospace")
    ax.axvline(0, color="#888", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8)
    ax.set_xlabel("Spearman ρ  (malignant CLDN4 vs T/NK fraction)")
    ax.set_xlim(-1.05, 1.05)
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def scatter_two(p127: pd.DataFrame, p148: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6), sharey=False)
    for ax, df, title in (
        (axes[0], p127, "GSE127465 (author malignant)"),
        (axes[1], p148, "GSE148071 (putative epithelium)"),
    ):
        e = df[df["eligible"]]
        ax.scatter(e["CLDN4_mean_log1p"], e["frac_TNK"], c="#1f4e79", s=36)
        for _, r in e.iterrows():
            ax.annotate(str(r["patient"]), (r["CLDN4_mean_log1p"], r["frac_TNK"]), fontsize=7, alpha=0.8)
        rho, p, n = spearman_safe(e["CLDN4_mean_log1p"], e["frac_TNK"])
        ax.set_title(f"{title}\nn={n}  ρ={rho:+.2f}  p={p:.3g}", fontsize=9)
        ax.set_xlabel("malignant CLDN4 mean log1p")
        ax.set_ylabel("T/NK fraction")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=Path("/tmp/geo_pair"))
    ap.add_argument("--out", type=Path, default=ROOT / "results")
    ap.add_argument("--skip-cellchat", action="store_true")
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    db = ROOT / "db"

    # ---- inventory / gate ----
    d127 = args.cache / "gse127465"
    d148 = args.cache / "gse148071"
    mtx = download(GSE127465["base"] + GSE127465["mtx"], d127 / GSE127465["mtx"], GSE127465["bytes_mtx"])
    meta_p = download(GSE127465["base"] + GSE127465["meta"], d127 / GSE127465["meta"])
    genes_p = download(GSE127465["base"] + GSE127465["genes"], d127 / GSE127465["genes"])
    tar = download(GSE148071["tar_url"], d148 / GSE148071["tar"], GSE148071["bytes_tar"])

    inv = pd.DataFrame(
        [
            {
                "series": "GSE127465",
                "file": GSE127465["mtx"],
                "bytes": mtx.stat().st_size,
                "bytes_gb": mtx.stat().st_size / 1e9,
                "under_2gb": mtx.stat().st_size < 2e9,
                "role": "human normalized MTX (cells x genes)",
            },
            {
                "series": "GSE127465",
                "file": GSE127465["meta"],
                "bytes": meta_p.stat().st_size,
                "bytes_gb": meta_p.stat().st_size / 1e9,
                "under_2gb": True,
                "role": "human cell metadata",
            },
            {
                "series": "GSE148071",
                "file": GSE148071["tar"],
                "bytes": tar.stat().st_size,
                "bytes_gb": tar.stat().st_size / 1e9,
                "under_2gb": tar.stat().st_size < 2e9,
                "role": "42 per-sample UMI TXT",
            },
        ]
    )
    inv.to_csv(out / "geo_file_inventory.tsv", sep="\t", index=False)
    if not bool(inv["under_2gb"].all()):
        raise SystemExit("a processed file is ≥2 GB; stop")

    genes127 = load_127465_genes(genes_p)
    gene_set_127 = set(genes127)
    if "CLDN4" not in gene_set_127:
        raise SystemExit("CLDN4 absent from GSE127465 gene list; stop")

    # Wanted symbols = lineage markers + every CellChatDB protein subunit (filter per matrix later).
    inter_raw = pd.read_csv(db / "interaction_cellchatdb_v2_protein.csv")
    wanted = set(EXTRA)
    for vs in LINEAGE_MARKERS.values():
        wanted.update(vs)
    for col in ("ligand.symbol", "receptor.symbol", "ligand", "receptor"):
        if col in inter_raw.columns:
            for val in inter_raw[col]:
                wanted.update(parse_symbols(val))
    print(f"wanted genes {len(wanted)} (pre-filter)", flush=True)

    # ---- GSE127465 expression ----
    print("GSE127465 metadata + MTX", flush=True)
    meta = load_127465_meta(meta_p)
    col_1based = {}
    for i, g in enumerate(genes127, start=1):
        if g in wanted:
            col_1based[i] = g
    print(f"  streaming {len(col_1based)} gene columns from MTX", flush=True)
    expr127 = stream_mtx_columns(mtx, col_1based, n_rows=len(meta))
    if "CLDN4" not in expr127:
        raise SystemExit("CLDN4 column missing after MTX stream; stop")

    p127 = score_127465_patients(meta, expr127)
    p127.to_csv(out / "gse127465_patients.tsv", sep="\t", index=False)
    e127 = p127[p127["eligible"]]
    rho127, pval127, n127 = spearman_safe(e127["CLDN4_mean_log1p"], e127["frac_TNK"])
    lo127, hi127 = rho_ci(rho127, n127)
    print(f"GSE127465 eligible n={n127} ρ={rho127:+.3f} p={pval127:.4g}", flush=True)

    # ---- GSE148071 ----
    print("GSE148071 matrices", flush=True)
    files148 = extract_148071(tar, d148 / "files")
    # probe genes from first file
    with gzip.open(files148[0], "rt") as handle:
        handle.readline()
        genes148 = []
        for line in handle:
            genes148.append(line.split("\t", 1)[0].strip().strip('"').split(".")[0])
    gene_set_148 = set(genes148)
    if "CLDN4" not in gene_set_148:
        raise SystemExit("CLDN4 absent from GSE148071; stop")
    wanted148 = {g for g in wanted if g in gene_set_148}
    print(f"  GSE148071 genes {len(genes148)}; wanted present {len(wanted148)}", flush=True)
    cells148, expr148_raw, umi148 = load_148071(files148, wanted148)
    lin_tab = (
        cells148.groupby(["patient", "lineage"]).size().rename("n").reset_index()
    )
    lin_tab.to_csv(out / "gse148071_lineage_n.tsv", sep="\t", index=False)

    p148 = score_148071_patients(cells148)
    p148.to_csv(out / "gse148071_patients.tsv", sep="\t", index=False)
    e148 = p148[p148["eligible"]]
    rho148, pval148, n148 = spearman_safe(e148["CLDN4_mean_log1p"], e148["frac_TNK"])
    lo148, hi148 = rho_ci(rho148, n148)
    print(f"GSE148071 eligible n={n148} ρ={rho148:+.3f} p={pval148:.4g}", flush=True)

    fz = fisher_z_diff(rho127, n127, rho148, n148)
    differ, why = associations_differ(rho127, n127, rho148, n148)
    print(f"differ={differ} ({why})  fisher_z={fz}", flush=True)

    patients = pd.concat([p127, p148], ignore_index=True)
    patients.to_csv(out / "patients_cldn4_tnk.tsv", sep="\t", index=False)

    combo = pd.DataFrame(
        [
            {
                "cohort": "GSE127465",
                "assay": "inDrops normalized MTX; author PatientN-specific = malignant",
                "n_deposited": 7,
                "n_eligible": n127,
                "rho": rho127,
                "p": pval127,
                "ci95_lo": lo127,
                "ci95_hi": hi127,
                "endpoint": "T/NK fraction of tumor cells",
            },
            {
                "cohort": "GSE148071",
                "assay": "10x UMI TXT; marker-argmax epithelium = putative malignant",
                "n_deposited": 42,
                "n_eligible": n148,
                "rho": rho148,
                "p": pval148,
                "ci95_lo": lo148,
                "ci95_hi": hi148,
                "endpoint": "T/NK fraction of all cells",
            },
        ]
    )
    combo.to_csv(out / "combo_cldn4_tnk.tsv", sep="\t", index=False)

    forest_plot(
        [
            {
                "label": f"GSE127465  n={n127}/7",
                "rho": rho127,
                "ci95_lo": lo127,
                "ci95_hi": hi127,
                "p": pval127,
                "n": n127,
            },
            {
                "label": f"GSE148071  n={n148}/42",
                "rho": rho148,
                "ci95_lo": lo148,
                "ci95_hi": hi148,
                "p": pval148,
                "n": n148,
            },
        ],
        out / "forest_CLDN4_tnk.png",
        f"CLDN4 vs T/NK  ·  differ={differ} ({why})",
    )
    scatter_two(p127, p148, out / "scatter_patient_cldn4_tnk.png")

    # MW on eligible median split (descriptive)
    mw_rows = []
    for name, e in (("GSE127465", e127), ("GSE148071", e148)):
        if len(e) < 4:
            continue
        med = float(e["CLDN4_mean_log1p"].median())
        hi = e[e["CLDN4_mean_log1p"] >= med]["frac_TNK"]
        lo = e[e["CLDN4_mean_log1p"] < med]["frac_TNK"]
        if len(hi) and len(lo):
            u, p_mw = mannwhitneyu(hi, lo, alternative="two-sided")
            r_rb = 2 * u / (len(hi) * len(lo)) - 1
        else:
            u = p_mw = r_rb = float("nan")
        mw_rows.append(
            {
                "cohort": name,
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "median_TNK_high": float(hi.median()) if len(hi) else float("nan"),
                "median_TNK_low": float(lo.median()) if len(lo) else float("nan"),
                "rank_biserial": float(r_rb),
                "p_mw": float(p_mw),
                "cldn4_median": med,
            }
        )
    mw_df = pd.DataFrame(mw_rows)
    mw_df.to_csv(out / "median_split_tnk.tsv", sep="\t", index=False)

    # ---- CellChat-style if they differ ----
    cellchat_ran = False
    outgoing_paths = []
    if args.skip_cellchat:
        print("skip CellChat by flag", flush=True)
    elif not differ:
        print("cohorts do not differ; CellChat-style outgoing not run", flush=True)
    else:
        cellchat_ran = True
        print("cohorts differ — CellChat-style outgoing per cohort", flush=True)
        lr127 = load_lr(db, set(expr127))
        # 127465: use provided normalized values as expression
        tumor = (meta["Tissue"].astype(str) == "tumor").to_numpy()
        maj = meta["Major cell type"].astype(str)
        pat = meta["Patient"].astype(str)
        own_mal = np.array(
            [maj.iloc[i] == f"Patient{patient_num(pat.iloc[i])}-specific" for i in range(len(meta))]
        )
        tnk127 = maj.isin(["tT cells", "tNK cells"]).to_numpy()
        # eligible patients only
        elig_p = set(e127["patient"])
        elig_cell = tumor & np.array([p in elig_p for p in pat])
        keep127 = elig_cell & (own_mal | tnk127)
        cldn4_127 = np.log1p(expr127["CLDN4"])
        # CellChat on normalized (not log) values
        run_cellchat_outgoing(
            "GSE127465",
            expr127,
            keep127,
            own_mal & elig_cell,
            tnk127 & elig_cell,
            cldn4_127,
            lr127,
            out / "gse127465_cellchat",
        )
        outgoing_paths.append("gse127465_cellchat_outgoing.tsv")

        lr148 = load_lr(db, set(expr148_raw))
        scale = np.where(umi148 > 0, 1e4 / umi148, 0.0)
        expr148_cp = {g: (expr148_raw[g] * scale).astype(np.float32) for g in expr148_raw}
        elig_p148 = set(e148["patient"])
        elig_cell148 = cells148["patient"].isin(elig_p148).to_numpy()
        mal148 = (cells148["lineage"] == "Epithelial").to_numpy()
        tnk148 = cells148["lineage"].isin(["T", "NK"]).to_numpy()
        keep148 = elig_cell148 & (mal148 | tnk148)
        run_cellchat_outgoing(
            "GSE148071",
            expr148_cp,
            keep148,
            mal148 & elig_cell148,
            tnk148 & elig_cell148,
            cells148["CLDN4_log1p_cp10k"].to_numpy(),
            lr148,
            out / "gse148071_cellchat",
        )
        outgoing_paths.append("gse148071_cellchat_outgoing.tsv")

    honest = pd.DataFrame(
        [
            {"item": "GSE127465 patients deposited (human)", "n": 7, "note": "p1–p7; Zilionis 2019"},
            {"item": "GSE127465 tumor cells", "n": int((meta["Tissue"] == "tumor").sum()), "note": "not the test n"},
            {"item": "GSE127465 blood cells", "n": int((meta["Tissue"] == "blood").sum()), "note": "excluded from T/NK and malignant"},
            {
                "item": "GSE127465 author PatientN-specific (tumor, own patient)",
                "n": int(p127["n_malignant"].sum()),
                "note": "primary malignant",
            },
            {"item": "GSE127465 tumor T/NK", "n": int(p127["n_TNK"].sum()), "note": "tT cells + tNK cells"},
            {"item": "GSE127465 eligible patients (≥20 mal + ≥20 T/NK)", "n": n127, "note": "patient-level unit"},
            {"item": "GSE148071 patients deposited", "n": 42, "note": "Wu 2021; one biopsy each"},
            {"item": "GSE148071 cells in tar", "n": int(len(cells148)), "note": "not the test n"},
            {"item": "GSE148071 putative epithelial cells", "n": int((cells148["lineage"] == "Epithelial").sum()), "note": "marker-argmax; no GEO malignant label"},
            {"item": "GSE148071 T/NK cells", "n": int(cells148["lineage"].isin(["T", "NK"]).sum()), "note": "marker-argmax"},
            {"item": "GSE148071 eligible patients (≥20 epi + ≥20 T/NK)", "n": n148, "note": "patient-level unit"},
            {"item": "Pairwise combo patients (sum of eligible)", "n": n127 + n148, "note": "independent patients; not a merged object"},
            {"item": "Dual-high TACSTD2×CLDN4 patients", "n": 0, "note": "not defined"},
            {"item": "Triple/quad series added", "n": 0, "note": "pairwise only"},
            {"item": "CellChat-style outgoing run", "n": int(cellchat_ran), "note": why},
        ]
    )
    honest.to_csv(out / "honest_n.tsv", sep="\t", index=False)

    lineage127 = (
        meta.loc[meta["Tissue"] == "tumor", "Major cell type"].value_counts().rename_axis("lineage").reset_index(name="n")
    )
    lineage127.insert(0, "cohort", "GSE127465")
    lineage127.to_csv(out / "gse127465_tumor_lineage_n.tsv", sep="\t", index=False)

    summary = {
        "additive": True,
        "gene": "CLDN4",
        "dual_high": False,
        "triple_quad": False,
        "cell_level_merge": False,
        "gse127465": {
            "n_deposited": 7,
            "n_eligible": n127,
            "rho": rho127,
            "p": pval127,
            "ci95": [lo127, hi127],
            "malignant": "author PatientN-specific in that patient's tumor",
            "tnk": "tT cells + tNK cells",
            "CLDN4_present": True,
            "processed_bytes": mtx.stat().st_size,
        },
        "gse148071": {
            "n_deposited": 42,
            "n_eligible": n148,
            "rho": rho148,
            "p": pval148,
            "ci95": [lo148, hi148],
            "malignant": "marker-argmax epithelium (putative; no GEO label)",
            "tnk": "marker-argmax T+NK",
            "CLDN4_present": True,
            "processed_bytes": tar.stat().st_size,
        },
        "differ": differ,
        "differ_reason": why,
        "fisher_z": fz,
        "cellchat_outgoing": cellchat_ran,
        "cellchat_files": outgoing_paths,
        "median_split": mw_df.to_dict(orient="records"),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    print("wrote", out)


if __name__ == "__main__":
    main()
