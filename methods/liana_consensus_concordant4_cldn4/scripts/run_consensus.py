#!/usr/bin/env python3
"""LIANA+ / CellPhoneDB / Connectome consensus on concordant-4.

Malignant CLDN4 Q4 vs Q1 senders, receivers = T/NK and myeloid.
Patient / locked sample is the unit. CLDN4 only. No dual-high.
"""

from __future__ import annotations

import math
import os
import tarfile
import traceback
import warnings
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("CONCORDANT4_RAW", "/tmp/concordant4_raw"))
CACHE = Path(os.environ.get("CONCORDANT4_CACHE", "/tmp/concordant4_cache"))

MIN_ARM = 10
MIN_RECV = 20
MIN_MAL_Q4 = 40
MIN_TEST_N = 8
MAX_CELLS = 400
N_PERMS = 1000
SEED = 1337
EXPR_PROP = 0.1

ALIASES = {
    "PVRL2": "NECTIN2",
    "JAM1": "F11R",
    "IL8": "CXCL8",
    "PDL1": "CD274",
    "PDCD1LG1": "CD274",
    "PD1": "PDCD1",
    "CD155": "PVR",
    "CD112": "NECTIN2",
}
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK_MARKERS = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
MYE_MARKERS = ["LYZ", "CD68", "CD14", "FCGR3A", "CSF1R", "AIF1"]
MALIG_SUB = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES = {"T lymphocytes", "NK cells"}
MYE_TYPES = {"Myeloid cells"}

RANK_CLASSES = ["barrier_exclusion", "recruit_effector", "recruit_myeloid"]
CLASS_LABEL = {
    "barrier_exclusion": "Barrier / exclusion",
    "recruit_effector": "Effector recruitment",
    "recruit_myeloid": "Myeloid recruitment",
}
METHOD_COLS = [
    ("cpdb", "d_cpdb", "CellPhoneDB"),
    ("conn", "d_conn", "Connectome"),
    ("liana", "d_liana", "LIANA+"),
]


def log(msg: str) -> None:
    print(msg, flush=True)


def load_edges() -> pd.DataFrame:
    edges = pd.read_csv(HERE / "data" / "edges.tsv", sep="\t")
    edges["ligand"] = edges["ligand"].astype(str)
    edges["receptor"] = edges["receptor"].astype(str)
    return edges


def subunits(token: str) -> list[str]:
    return token.split("_")


def wanted_symbols(edges: pd.DataFrame) -> set[str]:
    genes = set(EPI + TNK_MARKERS + MYE_MARKERS + ["CLDN4", "PTPRC"])
    for col in ("ligand", "receptor"):
        for token in edges[col]:
            genes.update(subunits(token))
    genes.update(ALIASES.keys())
    return genes


def canon_name(gene: str) -> str:
    g = gene.upper()
    return ALIASES.get(g, g)


def collapse_named_rows(named: dict[str, np.ndarray], keep: set[str]) -> dict[str, np.ndarray]:
    """Map raw symbols through aliases. Official symbol wins over its alias."""
    buckets: dict[str, list[tuple[str, np.ndarray]]] = {}
    for gene, arr in named.items():
        g = gene.upper()
        canon = ALIASES.get(g, g)
        if canon not in keep and g not in keep:
            continue
        buckets.setdefault(canon, []).append((g, np.asarray(arr)))
    out: dict[str, np.ndarray] = {}
    for canon, items in buckets.items():
        if canon not in keep:
            continue
        official = [a for g, a in items if g == canon]
        use = official if official else [a for _, a in items]
        out[canon] = use[0] if len(use) == 1 else np.sum(np.vstack(use), axis=0)
    return out


def matrix_from_named(named: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, list[str]]:
    genes = sorted(named)
    if not genes:
        return np.zeros((n, 0), dtype=np.float32), []
    cols = [np.asarray(named[g], dtype=np.float32).reshape(-1) for g in genes]
    mat = np.column_stack(cols)
    if mat.shape[0] != n:
        raise RuntimeError(f"gene matrix rows {mat.shape[0]} != cells {n}")
    return mat, genes


def pos_any(mat: np.ndarray, genes: list[str], markers: list[str]) -> np.ndarray:
    idx = [genes.index(g) for g in markers if g in genes]
    if not idx:
        return np.zeros(mat.shape[0], dtype=bool)
    block = mat[:, idx]
    return np.any(block > 0, axis=1)


def marker_compartments(mat: np.ndarray, genes: list[str]) -> np.ndarray:
    mal = pos_any(mat, genes, EPI) & ~pos_any(mat, genes, ["PTPRC"])
    tnk = pos_any(mat, genes, TNK_MARKERS) & ~mal
    mye = pos_any(mat, genes, MYE_MARKERS) & ~mal & ~tnk
    out = np.full(mat.shape[0], "", dtype=object)
    out[mye] = "MYE"
    out[tnk] = "TNK"
    out[mal] = "MAL"
    return out


def quartile_high_low(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Match the CellChat script: rank ties ordinal, low = rank <= floor(n/4)."""
    n = int(x.size)
    high = np.zeros(n, dtype=bool)
    low = np.zeros(n, dtype=bool)
    if n < 4:
        return high, low
    r = stats.rankdata(np.asarray(x, dtype=float), method="ordinal")
    q1 = math.floor(n * 0.25)
    q4 = math.ceil(n * 0.75)
    if q1 < 1 or q4 > n or q1 >= q4:
        return high, low
    low = r <= q1
    high = r > q4
    return high, low


def lognorm(counts: np.ndarray, lib: np.ndarray) -> np.ndarray:
    scale = np.maximum(lib.astype(np.float64), 1.0)[:, None]
    return np.log1p(counts.astype(np.float64) / scale * 1e4)


def pack_unit(patient: str, mat: np.ndarray, genes: list[str], lib: np.ndarray, comp: np.ndarray) -> dict | None:
    keep = comp != ""
    if not np.any(keep):
        return None
    return {
        "patient": patient,
        "mat": np.asarray(mat[keep], dtype=np.float32),
        "genes": genes,
        "lib": np.asarray(lib[keep], dtype=np.float64),
        "comp": comp[keep],
    }


def read_dense_csv(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    df = pd.read_csv(path)
    genes = [str(c).upper() for c in df.columns[1:]]
    mat = df.iloc[:, 1:].to_numpy(dtype=np.float32)
    lib = mat.sum(axis=1).astype(np.float64)
    return mat, genes, lib


def load_gse123902(edges: pd.DataFrame) -> list[dict]:
    log("==== GSE123902 ====")
    keep = wanted_symbols(edges)
    units = pd.read_csv(HERE / "data" / "GSE123902_marker_units.tsv", sep="\t")
    tumor = units[units["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor["ord"] = np.where(tumor["tissue"].eq("PRIMARY"), 0, 1)
    tumor = tumor.sort_values(["patient", "ord"]).drop_duplicates("patient")
    csv_dir = RAW / "GSE123902" / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    if not list(csv_dir.glob("*_dense.csv.gz")):
        tar = RAW / "GSE123902" / "GSE123902_RAW.tar"
        log(f"  untar {tar.name}")
        with tarfile.open(tar) as tf:
            tf.extractall(csv_dir)
    out = []
    for row in tumor.itertuples(index=False):
        fp = csv_dir / row.file
        if not fp.exists():
            hits = sorted(p for p in csv_dir.glob(f"*{row.patient}*dense.csv.gz") if "NORMAL" not in p.name)
            if not hits:
                log(f"  missing {row.patient}")
                continue
            fp = hits[0]
        log(f"  {row.patient} {fp.name}")
        mat, genes, lib = read_dense_csv(fp)
        named = {g: mat[:, i] for i, g in enumerate(genes)}
        named = collapse_named_rows(named, keep)
        mat, genes = matrix_from_named(named, mat.shape[0])
        comp = marker_compartments(mat, genes)
        packed = pack_unit(str(row.patient), mat, genes, lib, comp)
        if packed:
            packed["cohort"] = "GSE123902"
            out.append(packed)
    return out


def load_gse131907(edges: pd.DataFrame) -> list[dict]:
    log("==== GSE131907 ====")
    keep = wanted_symbols(edges)
    extract = set(keep) | set(ALIASES)
    samples = pd.read_csv(HERE / "data" / "GSE131907_samples.tsv", sep="\t")
    locked = set(samples.loc[samples["n_malignant"] > 0, "sample"].astype(str))
    annot = pd.read_csv(
        RAW / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t"
    )
    path = RAW / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    log("  streaming UMI (wanted genes + library size)")
    import gzip

    with gzip.open(path, "rb") as fh:
        header = fh.readline().rstrip(b"\n").split(b"\t")
        cells = [h.decode() for h in header[1:]]
        n = len(cells)
        lib = np.zeros(n, dtype=np.float64)
        named: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            n_genes += 1
            tab = line.find(b"\t")
            gene = line[:tab].decode().upper()
            arr = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.float64)
            if arr.size != n:
                fixed = np.zeros(n, dtype=np.float64)
                m = min(n, arr.size)
                fixed[:m] = arr[:m]
                arr = fixed
            lib += arr
            if gene in extract:
                named[gene] = arr.astype(np.float32, copy=False)
            if n_genes % 4000 == 0:
                log(f"    streamed {n_genes} genes, kept {len(named)}")
    log(f"  streamed {n_genes} genes, cells {n}, kept {len(named)}")
    named = collapse_named_rows(named, keep)
    mat, genes = matrix_from_named(named, n)
    ann = annot.copy()
    ann["Index"] = ann["Index"].astype(str)
    ann = ann.drop_duplicates("Index").set_index("Index")
    pos = {c: i for i, c in enumerate(cells)}
    common = [c for c in ann.index if c in pos]
    if len(common) < 0.9 * n:
        raise RuntimeError(f"GSE131907 annotation overlap {len(common)}/{n}")
    idx = np.array([pos[c] for c in common])
    sub_ann = ann.loc[common]
    subtype = sub_ann["Cell_subtype"].astype(str)
    ctype = sub_ann["Cell_type"].astype(str)
    sample = sub_ann["Sample"].astype(str)
    comp = np.full(len(common), "", dtype=object)
    mal = subtype.isin(MALIG_SUB).to_numpy()
    tnk = ctype.isin(TNK_TYPES).to_numpy() & ~mal
    mye = ctype.isin(MYE_TYPES).to_numpy() & ~mal & ~tnk
    comp[mye] = "MYE"
    comp[tnk] = "TNK"
    comp[mal] = "MAL"
    mat = mat[idx]
    lib = lib[idx]
    out = []
    for s in sorted(locked):
        which = np.where(sample.to_numpy() == s)[0]
        if which.size == 0:
            log(f"  missing sample {s}")
            continue
        packed = pack_unit(s, mat[which], genes, lib[which], comp[which])
        if packed:
            packed["cohort"] = "GSE131907"
            out.append(packed)
            log(f"  {s} mal {(packed['comp']=='MAL').sum()} tnk {(packed['comp']=='TNK').sum()} mye {(packed['comp']=='MYE').sum()}")
    return out


def load_gse205335(edges: pd.DataFrame) -> list[dict]:
    log("==== GSE205335 ====")
    keep = wanted_symbols(edges)
    export = CACHE / "gse205335"
    mtx = export / "matrix.mtx"
    if not mtx.exists():
        import subprocess

        export.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(
            [
                "Rscript",
                str(HERE / "scripts" / "export_gse205335.R"),
                str(RAW / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"),
                str(HERE / "data" / "genes.txt"),
                str(export),
            ]
        )
    genes = [ln.strip() for ln in (export / "genes.txt").read_text().splitlines() if ln.strip()]
    cells = [ln.strip() for ln in (export / "cells.txt").read_text().splitlines() if ln.strip()]
    log(f"  reading MTX {len(genes)} x {len(cells)}")
    mm = spio.mmread(mtx).tocsr()
    if mm.shape == (len(genes), len(cells)):
        mm = mm.T.tocsr()
    elif mm.shape != (len(cells), len(genes)):
        raise RuntimeError(f"unexpected MTX shape {mm.shape}")
    lib_df = pd.read_csv(export / "libsize.tsv", sep="\t")
    lib_map = dict(zip(lib_df["cell"].astype(str), lib_df["libsize"].astype(float)))
    lib = np.array([lib_map[c] for c in cells], dtype=np.float64)
    named = {g: np.asarray(mm[:, i].todense()).ravel() for i, g in enumerate(genes)}
    named = collapse_named_rows(named, keep)
    mat, genes = matrix_from_named(named, len(cells))
    ident = pd.read_csv(RAW / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    gsm = pd.read_csv(HERE / "data" / "GSE205335_gsm_map.tsv", sep="\t")
    ident = ident.merge(gsm[["orig.ident", "patient", "tissue"]], on="orig.ident", how="left")
    ident["barcode"] = ident["barcode"].astype(str)
    idmap = ident.drop_duplicates("barcode").set_index("barcode")
    common = [i for i, c in enumerate(cells) if c in idmap.index]
    if len(common) < 0.9 * len(cells):
        raise RuntimeError(f"GSE205335 identity overlap {len(common)}/{len(cells)}")
    idx = np.array(common)
    sub = idmap.loc[[cells[i] for i in idx]]
    tissue = sub["tissue"].fillna("").astype(str)
    normal = tissue.str.startswith("Normal").to_numpy()
    lineage_sub = sub["lineage.sub"].astype(str).to_numpy()
    lineage = sub["lineage.total"].astype(str).to_numpy()
    patient = sub["patient"].astype(str).to_numpy()
    comp = np.full(idx.size, "", dtype=object)
    mal = (lineage_sub == "Malignant cells") & ~normal
    tnk = (lineage == "T/NK cells") & ~normal & ~mal
    mye = (lineage == "Myeloid cells") & ~normal & ~mal & ~tnk
    comp[mye] = "MYE"
    comp[tnk] = "TNK"
    comp[mal] = "MAL"
    mat = mat[idx]
    lib = lib[idx]
    locked = pd.read_csv(HERE / "data" / "GSE205335_patients.tsv", sep="\t")
    keep_pt = set(locked.loc[locked["n_malignant"] > 0, "patient"].astype(str))
    out = []
    for pt in sorted(keep_pt):
        which = np.where(patient == pt)[0]
        if which.size == 0:
            log(f"  missing patient {pt}")
            continue
        packed = pack_unit(pt, mat[which], genes, lib[which], comp[which])
        if packed:
            packed["cohort"] = "GSE205335"
            out.append(packed)
    log(f"  units {len(out)}")
    return out


def load_gse189357(edges: pd.DataFrame) -> list[dict]:
    log("==== GSE189357 ====")
    keep = wanted_symbols(edges)
    meta = pd.read_csv(HERE / "data" / "GSE189357_sample_metadata.tsv", sep="\t")
    ex = RAW / "GSE189357" / "raw"
    ex.mkdir(parents=True, exist_ok=True)
    needed = ex / f"{meta.gsm.iloc[0]}_{meta.patient.iloc[0]}_matrix.mtx.gz"
    if not needed.exists():
        tar = RAW / "GSE189357" / "GSE189357_RAW.tar"
        log(f"  untar {tar.name}")
        with tarfile.open(tar) as tf:
            tf.extractall(ex)
    out = []
    for row in meta.itertuples(index=False):
        prefix = ex / f"{row.gsm}_{row.patient}"
        mtx = Path(str(prefix) + "_matrix.mtx.gz")
        feat = Path(str(prefix) + "_features.tsv.gz")
        bar = Path(str(prefix) + "_barcodes.tsv.gz")
        if not mtx.exists():
            hits = list(ex.glob(f"*{row.patient}_matrix.mtx.gz"))
            if not hits:
                log(f"  missing {row.patient}")
                continue
            mtx = hits[0]
            feat = Path(str(mtx).replace("_matrix.mtx.gz", "_features.tsv.gz"))
            bar = Path(str(mtx).replace("_matrix.mtx.gz", "_barcodes.tsv.gz"))
        log(f"  {row.patient}")
        mm = spio.mmread(mtx).tocsr()
        features = pd.read_csv(feat, sep="\t", header=None)
        symbols = features.iloc[:, 1 if features.shape[1] > 1 else 0].astype(str).str.upper()
        if mm.shape[0] != len(symbols):
            if mm.shape[1] == len(symbols):
                mm = mm.T.tocsr()
            else:
                raise RuntimeError(f"{row.patient} mtx {mm.shape} vs features {len(symbols)}")
        n = mm.shape[1]
        lib = np.asarray(mm.sum(axis=0)).ravel().astype(np.float64)
        summed: dict[str, np.ndarray] = {}
        for i, g in enumerate(symbols):
            if g not in keep and g not in ALIASES:
                continue
            arr = np.asarray(mm[i].todense()).ravel()
            if g in summed:
                summed[g] = summed[g] + arr
            else:
                summed[g] = arr
        named = collapse_named_rows(summed, keep)
        mat, genes = matrix_from_named(named, n)
        comp = marker_compartments(mat, genes)
        packed = pack_unit(str(row.patient), mat, genes, lib, comp)
        if packed:
            packed["cohort"] = "GSE189357"
            out.append(packed)
    return out


def qc_against_cellchat(units: list[dict]) -> pd.DataFrame:
    rows = []
    for u in units:
        comp = u["comp"]
        mal = comp == "MAL"
        cldn = gene_vector(u, "CLDN4")
        lib = u["lib"]
        if np.any(mal) and cldn is not None:
            ln = np.log1p(cldn[mal] / np.maximum(lib[mal], 1) * 1e4)
            mean_ln = float(ln.mean())
            pct = float(np.mean(cldn[mal] > 0))
        else:
            mean_ln, pct = np.nan, np.nan
        rows.append(
            {
                "cohort": u["cohort"],
                "patient": u["patient"],
                "n_cells_kept": int(comp.size),
                "n_mal": int(mal.sum()),
                "n_tnk": int((comp == "TNK").sum()),
                "n_mye": int((comp == "MYE").sum()),
                "mal_cldn4_mean": mean_ln,
                "mal_cldn4_pct": pct,
            }
        )
    inv = pd.DataFrame(rows)
    ref = pd.read_csv(HERE / "data" / "reference_cellchat_inventory.tsv", sep="\t")
    m = inv.merge(ref[["cohort", "patient", "n_mal", "n_tnk"]], on=["cohort", "patient"], suffixes=("", "_ref"))
    if m.empty:
        raise RuntimeError("inventory did not match the CellChat reference keys")
    d_mal = (m["n_mal"] - m["n_mal_ref"]).abs()
    d_tnk = (m["n_tnk"] - m["n_tnk_ref"]).abs()
    log(
        f"QC vs CellChat labels: units {len(m)}/{len(ref)} "
        f"max|Δn_mal|={int(d_mal.max())} max|Δn_tnk|={int(d_tnk.max())}"
    )
    bad = m[(d_mal > 0) | (d_tnk > 0)]
    if len(bad):
        log(bad[["cohort", "patient", "n_mal", "n_mal_ref", "n_tnk", "n_tnk_ref"]].head(12).to_string(index=False))
    if int(d_mal.max()) > 0 or int(d_tnk.max()) > 0:
        raise RuntimeError("malignant/T/NK counts do not match the locked CellChat inventory")
    missing = set(zip(ref.cohort, ref.patient.astype(str))) - set(zip(inv.cohort, inv.patient.astype(str)))
    if missing:
        raise RuntimeError(f"missing locked units: {sorted(missing)[:8]}")
    return inv


def gene_vector(unit: dict, gene: str) -> np.ndarray | None:
    if gene not in unit["genes"]:
        return None
    return unit["mat"][:, unit["genes"].index(gene)]


def run_liana_receiver(unit: dict, high: np.ndarray, low: np.ndarray, receiver: str, edges: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame | None:
    comp = unit["comp"]
    groups = {
        "CLDN4_high": np.where(high)[0],
        "CLDN4_low": np.where(low)[0],
        receiver: np.where(comp == receiver)[0],
    }
    if min(len(v) for v in groups.values()) < MIN_ARM:
        return None
    take = []
    labels = []
    for name, idx in groups.items():
        if idx.size > MAX_CELLS:
            idx = rng.choice(idx, size=MAX_CELLS, replace=False)
        take.append(idx)
        labels.extend([name] * idx.size)
    idx = np.concatenate(take)
    counts = unit["mat"][idx]
    lib = unit["lib"][idx]
    expr = lognorm(counts, lib)
    keep_gene = expr.sum(axis=0) > 0
    genes = [g for g, k in zip(unit["genes"], keep_gene) if k]
    expr = expr[:, keep_gene]
    adata = ad.AnnData(expr)
    adata.var_names = genes
    adata.obs_names = [f"c{i}" for i in range(expr.shape[0])]
    adata.obs["label"] = pd.Categorical(labels)
    resource = edges.loc[:, ["ligand", "receptor"]].drop_duplicates()
    pairs = pd.DataFrame({"source": ["CLDN4_high", "CLDN4_low"], "target": [receiver, receiver]})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from liana.method import rank_aggregate

        df = rank_aggregate(
            adata,
            groupby="label",
            resource=resource,
            groupby_pairs=pairs,
            expr_prop=EXPR_PROP,
            min_cells=5,
            n_perms=N_PERMS,
            seed=SEED,
            use_raw=False,
            inplace=False,
            verbose=False,
            return_all_lrs=True,
        )
    return df


def _finite(row: pd.Series | None, col: str) -> float:
    if row is None:
        return np.nan
    v = row[col]
    try:
        v = float(v)
    except (TypeError, ValueError):
        return np.nan
    return v if math.isfinite(v) else np.nan


def _lookup(df: pd.DataFrame, source: str, ligand: str, receptor: str) -> pd.Series | None:
    m = df[
        (df["source"] == source)
        & (df["ligand_complex"] == ligand)
        & (df["receptor_complex"] == receptor)
    ]
    if m.empty:
        return None
    return m.iloc[0]


def score_unit(unit: dict, edges: pd.DataFrame, rng: np.random.Generator) -> list[dict]:
    comp = unit["comp"]
    mal = comp == "MAL"
    n_mal = int(mal.sum())
    cldn = gene_vector(unit, "CLDN4")
    if cldn is None or n_mal < MIN_MAL_Q4:
        return []
    ln = np.log1p(cldn / np.maximum(unit["lib"], 1.0) * 1e4)
    high_m, low_m = quartile_high_low(ln[mal])
    high = np.zeros(comp.size, dtype=bool)
    low = np.zeros(comp.size, dtype=bool)
    high[np.where(mal)[0]] = high_m
    low[np.where(mal)[0]] = low_m
    if high.sum() < MIN_ARM or low.sum() < MIN_ARM:
        return []
    if float(ln[high].mean()) <= float(ln[low].mean()):
        raise RuntimeError(f"CLDN4 quartile inverted in {unit['cohort']} {unit['patient']}")
    rows = []
    for receiver in ("TNK", "MYE"):
        n_recv = int((comp == receiver).sum())
        if n_recv < MIN_RECV:
            continue
        log(
            f"  LIANA {unit['cohort']} {unit['patient']} {receiver} "
            f"high {int(high.sum())} low {int(low.sum())} recv {n_recv}"
        )
        try:
            df = run_liana_receiver(unit, high, low, receiver, edges, rng)
        except Exception:
            log(traceback.format_exc())
            continue
        if df is None or df.empty:
            continue
        pending = []
        for edge in edges.itertuples(index=False):
            hi = _lookup(df, "CLDN4_high", edge.ligand, edge.receptor)
            lo = _lookup(df, "CLDN4_low", edge.ligand, edge.receptor)
            cp_hi, cp_lo = _finite(hi, "lr_means"), _finite(lo, "lr_means")
            cn_hi, cn_lo = _finite(hi, "expr_prod"), _finite(lo, "expr_prod")
            rk_hi, rk_lo = _finite(hi, "magnitude_rank"), _finite(lo, "magnitude_rank")
            sp_hi, sp_lo = _finite(hi, "specificity_rank"), _finite(lo, "specificity_rank")
            p_hi, p_lo = _finite(hi, "cellphone_pvals"), _finite(lo, "cellphone_pvals")

            def delta(a: float, b: float) -> float:
                if not math.isfinite(a) and not math.isfinite(b):
                    return np.nan
                return (0.0 if not math.isfinite(a) else a) - (0.0 if not math.isfinite(b) else b)

            d_rank = (rk_lo - rk_hi) if math.isfinite(rk_hi) and math.isfinite(rk_lo) else np.nan
            d_spec = (sp_lo - sp_hi) if math.isfinite(sp_hi) and math.isfinite(sp_lo) else np.nan
            pending.append(
                {
                    "cohort": unit["cohort"],
                    "patient": unit["patient"],
                    "receiver": receiver,
                    "ligand": edge.ligand,
                    "receptor": edge.receptor,
                    "lr_class": edge.lr_class,
                    "axis": edge.axis,
                    "expect": edge.expect,
                    "n_mal": n_mal,
                    "n_high": int(high.sum()),
                    "n_low": int(low.sum()),
                    "n_recv": n_recv,
                    "cpdb_high": cp_hi,
                    "cpdb_low": cp_lo,
                    "d_cpdb": delta(cp_hi, cp_lo),
                    "conn_high": cn_hi,
                    "conn_low": cn_lo,
                    "d_conn": delta(cn_hi, cn_lo),
                    "liana_rho_high": rk_hi,
                    "liana_rho_low": rk_lo,
                    "d_liana": d_rank,
                    "d_spec": d_spec,
                    "cpdb_p_high": p_hi,
                    "cpdb_p_low": p_lo,
                }
            )
        for arm, pcol, qcol in (
            ("high", "cpdb_p_high", "cpdb_q_high"),
            ("low", "cpdb_p_low", "cpdb_q_low"),
        ):
            pvals = np.array([r[pcol] for r in pending], dtype=float)
            qvals = np.full(pvals.shape, np.nan)
            ok = np.isfinite(pvals)
            if ok.sum() >= 2:
                _, q, _, _ = multipletests(pvals[ok], method="fdr_bh")
                qvals[ok] = q
            elif ok.sum() == 1:
                qvals[ok] = pvals[ok]
            for r, qv in zip(pending, qvals):
                r[qcol] = float(qv) if math.isfinite(qv) else np.nan
        rows.extend(pending)
    return rows


def wilcox_p(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < MIN_TEST_N:
        return np.nan
    if np.all(x == 0) or np.sum(x != 0) < 1:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = stats.wilcoxon(x, zero_method="wilcox", alternative="two-sided", method="auto")
    return float(res.pvalue)


def bh_by(df: pd.DataFrame, pcol: str, qcol: str, keys: list[str]) -> None:
    df[qcol] = np.nan
    for _, sub in df.groupby(keys, sort=False):
        p = sub[pcol].to_numpy(dtype=float)
        ok = np.isfinite(p)
        if ok.sum() == 0:
            continue
        q = np.full(p.shape, np.nan)
        if ok.sum() == 1:
            q[ok] = p[ok]
        else:
            _, qq, _, _ = multipletests(p[ok], method="fdr_bh")
            q[ok] = qq
        df.loc[sub.index, qcol] = q


def summarize_edges(per: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grouped = per.groupby(["receiver", "ligand", "receptor"], sort=False)
    for (receiver, ligand, receptor), sub in grouped:
        rec = {
            "receiver": receiver,
            "ligand": ligand,
            "receptor": receptor,
            "lr_class": sub["lr_class"].iloc[0],
            "axis": sub["axis"].iloc[0],
            "expect": sub["expect"].iloc[0],
        }
        for key, dcol, _ in METHOD_COLS:
            # n is the number of units with a finite delta on ANY method; method-specific n stored too
            x = sub[dcol].to_numpy(dtype=float)
            finite = sub[np.isfinite(sub[dcol])].copy()
            rec[f"n_{key}"] = int(np.isfinite(x).sum())
            rec[f"mean_{dcol}"] = float(np.nanmean(x)) if np.isfinite(x).any() else np.nan
            rec[f"p_{key}"] = wilcox_p(x)
            signs = []
            for cohort, csub in finite.groupby("cohort"):
                if len(csub) < 3:
                    signs.append(f"{cohort}:na")
                    continue
                m = float(csub[dcol].mean())
                signs.append(f"{cohort}:{'+' if m > 0 else '-' if m < 0 else '0'}")
            rec[f"cohort_sign_{key}"] = ";".join(signs)
        any_finite = np.isfinite(sub[["d_cpdb", "d_conn", "d_liana"]].to_numpy(dtype=float)).any(axis=1)
        used = sub.loc[any_finite]
        rec["n"] = int(used[["cohort", "patient"]].drop_duplicates().shape[0])
        for cohort in ("GSE123902", "GSE131907", "GSE205335", "GSE189357"):
            rec[f"n_{cohort}"] = int((used["cohort"] == cohort).sum())
        rows.append(rec)
    out = pd.DataFrame(rows)
    for key, _, _ in METHOD_COLS:
        bh_by(out, f"p_{key}", f"q_{key}", ["receiver"])
    zs = []
    calls = []
    supports = []
    scores = []
    sigs = []
    for rec in out.itertuples(index=False):
        z_parts = []
        sig_high = []
        sig_low = []
        for key, dcol, name in METHOD_COLS:
            p = getattr(rec, f"p_{key}")
            q = getattr(rec, f"q_{key}")
            mean = getattr(rec, f"mean_{dcol}")
            if math.isfinite(p) and math.isfinite(mean) and p > 0:
                z_parts.append(float(np.sign(mean) * stats.norm.isf(min(p, 1) / 2.0)))
            elif math.isfinite(p) and p == 0 and math.isfinite(mean):
                z_parts.append(float(np.sign(mean) * 40))
            if math.isfinite(q) and q < 0.05 and math.isfinite(mean) and mean != 0:
                (sig_high if mean > 0 else sig_low).append(key)
        z = float(np.sum(z_parts) / math.sqrt(len(z_parts))) if z_parts else np.nan
        if len(sig_high) >= 2 and len(sig_low) == 0:
            call = "high>low"
        elif len(sig_low) >= 2 and len(sig_high) == 0:
            call = "low>high"
        else:
            call = "none"
        direction_score = z
        if rec.lr_class == "barrier_exclusion":
            support_score = z
            claim = "barrier_exclusion" if call == "high>low" else "none"
        elif rec.lr_class == "recruit_effector":
            support_score = -z if math.isfinite(z) else np.nan
            if call == "low>high":
                claim = "recruitment_from_low"
            elif call == "high>low":
                claim = "recruitment_from_high"
            else:
                claim = "none"
        elif rec.lr_class == "recruit_myeloid":
            support_score = abs(z) if math.isfinite(z) else np.nan
            if call == "high>low":
                claim = "myeloid_recruitment_from_high"
            elif call == "low>high":
                claim = "myeloid_recruitment_from_low"
            else:
                claim = "none"
        else:
            support_score = z
            claim = "none"
        zs.append(z)
        calls.append(call)
        supports.append(claim)
        scores.append(support_score)
        sigs.append("+".join(sig_high + [f"{s}(low)" for s in sig_low]) if (sig_high or sig_low) else "")
    out["consensus_z"] = zs
    out["consensus_call"] = calls
    out["supports"] = supports
    out["support_score"] = scores
    out["sig_methods"] = sigs
    return out


def summarize_families(per: pd.DataFrame) -> pd.DataFrame:
    rows = []
    classes = ["barrier_exclusion", "recruit_effector", "recruit_myeloid"]
    for receiver in ("TNK", "MYE"):
        for lr_class in classes:
            block = per[(per["receiver"] == receiver) & (per["lr_class"] == lr_class)]
            for key, dcol, _ in METHOD_COLS:
                agg_rows = []
                for (cohort, patient), sub in block.groupby(["cohort", "patient"]):
                    x = sub[dcol].to_numpy(dtype=float)
                    if not np.isfinite(x).any():
                        continue
                    agg_rows.append({"cohort": cohort, "patient": patient, "delta": float(np.nanmean(x))})
                agg = pd.DataFrame(agg_rows)
                n = int(len(agg))
                p = wilcox_p(agg["delta"].to_numpy()) if n else np.nan
                mean = float(agg["delta"].mean()) if n else np.nan
                n_by = agg["cohort"].value_counts().to_dict() if n else {}
                rows.append(
                    {
                        "receiver": receiver,
                        "lr_class": lr_class,
                        "method": key,
                        "n": n,
                        "n_GSE123902": int(n_by.get("GSE123902", 0)),
                        "n_GSE131907": int(n_by.get("GSE131907", 0)),
                        "n_GSE205335": int(n_by.get("GSE205335", 0)),
                        "n_GSE189357": int(n_by.get("GSE189357", 0)),
                        "mean_delta": mean,
                        "p_W": p,
                        "expect": {"barrier_exclusion": "high>low", "recruit_effector": "low>high", "recruit_myeloid": "observed"}[lr_class],
                    }
                )
    out = pd.DataFrame(rows)
    p = out["p_W"].to_numpy(dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() >= 2:
        _, qq, _, _ = multipletests(p[ok], method="fdr_bh")
        q[ok] = qq
    elif ok.sum() == 1:
        q[ok] = p[ok]
    out["q_BH"] = q
    obs = []
    agrees = []
    for rec in out.itertuples(index=False):
        if not math.isfinite(rec.mean_delta):
            o, a = "NA", "thin"
        elif rec.mean_delta > 0:
            o = "high>low"
            a = "yes" if rec.expect == "high>low" else ("observed" if rec.expect == "observed" else "opposite")
        elif rec.mean_delta < 0:
            o = "low>high"
            a = "yes" if rec.expect == "low>high" else ("observed" if rec.expect == "observed" else "opposite")
        else:
            o, a = "tie", "tie"
        obs.append(o)
        agrees.append(a)
    out["observed"] = obs
    out["agrees"] = agrees
    return out


def specificity_summary(per: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (receiver, ligand, receptor), sub in per.groupby(["receiver", "ligand", "receptor"], sort=False):
        def frac(col: str) -> float:
            x = sub[col].to_numpy(dtype=float)
            ok = np.isfinite(x)
            if not ok.any():
                return np.nan
            return float(np.mean(x[ok] < 0.05))

        rows.append(
            {
                "receiver": receiver,
                "ligand": ligand,
                "receptor": receptor,
                "lr_class": sub["lr_class"].iloc[0],
                "axis": sub["axis"].iloc[0],
                "n_p_high": int(np.isfinite(sub["cpdb_p_high"]).sum()),
                "n_p_low": int(np.isfinite(sub["cpdb_p_low"]).sum()),
                "frac_q_high": frac("cpdb_q_high"),
                "frac_q_low": frac("cpdb_q_low"),
            }
        )
    return pd.DataFrame(rows)


def fmt_p(p: float) -> str:
    if not math.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_n(x: float, digits: int = 3) -> str:
    if not math.isfinite(x):
        return "NA"
    return f"{x:+.{digits}f}"


def plot_consensus(summary: pd.DataFrame, path_png: Path, path_pdf: Path) -> None:
    methods = [("q_cpdb", "mean_d_cpdb", "CellPhoneDB"), ("q_conn", "mean_d_conn", "Connectome"), ("q_liana", "mean_d_liana", "LIANA+")]
    receivers = [("TNK", "T/NK"), ("MYE", "Myeloid")]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 11.4), sharex=True)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad("#e6e6e6")
    last_im = None
    for ax, (receiver, title) in zip(axes, receivers):
        blocks = []
        ylabels = []
        seams = []
        sub = summary[summary["receiver"] == receiver]
        for lr_class in RANK_CLASSES:
            block = sub[sub["lr_class"] == lr_class].copy()
            if block.empty:
                continue
            block = block.sort_values("support_score", ascending=False, na_position="last")
            mat = np.full((len(block), 3), np.nan)
            for i, rec in enumerate(block.itertuples(index=False)):
                for j, (qcol, mcol, _) in enumerate(methods):
                    q = getattr(rec, qcol)
                    mean = getattr(rec, mcol)
                    if math.isfinite(q) and math.isfinite(mean) and q > 0:
                        mat[i, j] = np.sign(mean) * min(6.0, -math.log10(q))
                    elif math.isfinite(q) and q == 0 and math.isfinite(mean):
                        mat[i, j] = np.sign(mean) * 6.0
                tag = ""
                if rec.consensus_call == "high>low":
                    tag = " *"
                elif rec.consensus_call == "low>high":
                    tag = " *"
                ylabels.append(f"{rec.axis}{tag}")
            if blocks:
                seams.append(sum(b.shape[0] for b in blocks))
            blocks.append(mat)
        if not blocks:
            ax.set_axis_off()
            continue
        data = np.vstack(blocks)
        masked = np.ma.masked_invalid(data)
        last_im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-6, vmax=6, interpolation="nearest")
        ax.set_yticks(np.arange(data.shape[0]))
        ax.set_yticklabels(ylabels, fontsize=7.5)
        ax.set_xticks(np.arange(3))
        ax.set_xticklabels([m[2] for m in methods], fontsize=9)
        ax.set_title(title, fontsize=12, pad=8)
        for y in seams:
            ax.axhline(y - 0.5, color="#222222", lw=0.8)
        for i in range(data.shape[0]):
            for j in range(3):
                if not np.isfinite(data[i, j]):
                    continue
                if abs(data[i, j]) >= -math.log10(0.05):
                    ax.text(j, i, "·", ha="center", va="center", color="black", fontsize=8)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=axes, fraction=0.03, pad=0.04)
        cbar.set_label("signed −log10 BH q   (red: higher from CLDN4-high)", fontsize=8)
        cbar.ax.tick_params(labelsize=8)
    fig.suptitle(
        "Concordant-4  ·  malignant CLDN4 Q4 vs Q1 outgoing LR",
        fontsize=13,
        y=0.98,
    )
    fig.text(
        0.5,
        0.005,
        "Rows within each class are ranked by support for that class.  * = ≥2 methods BH q<0.05, same sign, no contradicting method.  Gray = not tested (n<8).",
        ha="center",
        fontsize=7.5,
        color="#333333",
    )
    path_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_png, dpi=160, bbox_inches="tight")
    fig.savefig(path_pdf, bbox_inches="tight")
    plt.close(fig)


def write_finding(inv: pd.DataFrame, summary: pd.DataFrame, families: pd.DataFrame, versions: dict) -> None:
    both = inv[(inv["n_mal"] >= MIN_MAL_Q4) & ((inv["n_tnk"] >= MIN_RECV) | (inv["n_mye"] >= MIN_RECV))]
    # Q4 units are those that entered per-edge table; use summary n max as a check, inventory gate is both.

    def fam_line(receiver: str, lr_class: str, method: str) -> str:
        r = families[(families.receiver == receiver) & (families.lr_class == lr_class) & (families.method == method)]
        if r.empty:
            return "| | | | | | | | | | |"
        r = r.iloc[0]
        return (
            f"| {receiver} | {lr_class} | {method} | {int(r.n)} | {int(r.n_GSE123902)} | {int(r.n_GSE131907)} | "
            f"{int(r.n_GSE205335)} | {int(r.n_GSE189357)} | {fmt_n(r.mean_delta)} | {fmt_p(r.p_W)} | {fmt_p(r.q_BH)} | "
            f"{r.observed} | {r.agrees} |"
        )

    hits = summary[summary["supports"] != "none"].sort_values(["supports", "support_score"], ascending=[True, False])

    def hit_line(r: pd.Series) -> str:
        return (
            f"| {r.receiver} | {r.axis} | {r.lr_class} | {r.supports} | {int(r.n)} | "
            f"{fmt_n(r.mean_d_cpdb)} | {fmt_p(r.q_cpdb)} | {fmt_n(r.mean_d_conn)} | {fmt_p(r.q_conn)} | "
            f"{fmt_n(r.mean_d_liana)} | {fmt_p(r.q_liana)} | {fmt_n(r.consensus_z, 2)} | {r.sig_methods} |"
        )

    def top_lines(receiver: str, lr_class: str, k: int = 5) -> list[str]:
        sub = summary[(summary.receiver == receiver) & (summary.lr_class == lr_class)].copy()
        sub = sub.sort_values("support_score", ascending=False, na_position="last").head(k)
        return [hit_line(r) for _, r in sub.iterrows()]

    n_by = inv.groupby("cohort")["patient"].nunique().to_dict()
    lines = [
        "# FINDING — LIANA+ / CellPhoneDB / Connectome consensus, concordant-4",
        "",
        "ADDITIVE. **CLDN4 only. No dual-high.** Concordant four only",
        "(GSE123902 + GSE131907 + GSE205335 + GSE189357).",
        "Do **not** add GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526.",
        "This is an expression ligand–receptor contrast. It is **not** a spatial exclusion test.",
        "",
        f"Engine: Python liana {versions['liana']} `rank_aggregate`",
        "(CellPhoneDB + Connectome + log2FC + NATMI + SingleCellSignalR; RobustRankAggregate).",
        "The three reported scores, from that one call, are:",
        "",
        "- **CellPhoneDB** magnitude `lr_means` (Δ = high − low).",
        "- **Connectome** magnitude `expr_prod` (Δ = high − low).",
        "- **LIANA+** `magnitude_rank` (RobustRankAggregate ρ; Δ = ρ_low − ρ_high, so positive means CLDN4-high ranks stronger).",
        "",
        "Senders = malignant CLDN4 Q4 vs Q1. Receivers = T/NK and myeloid, run separately.",
        "Honest n = patient / locked sample. Primary test = Wilcoxon signed-rank of Δ across units.",
        "BH is within method × receiver across the pre-specified edge panel, and separately",
        "across the family tests, and within each unit for CellPhoneDB permutation p-values.",
        "",
        "Pre-specified classes:",
        "",
        "- **Barrier / exclusion** (expect high > low): junction and inhibitory edges",
        "  (F11R, NECTIN2/PVR–TIGIT/CD96, CDH1, LGALS9, PD-1 ligands, HLA-E–NKG2A,",
        "  CD47–SIRPA, CD24–SIGLEC10, TGFB1, MIF–CD74).",
        "- **Effector recruitment** (expect low > high; the KD-like arm): CXCL9/10/11–CXCR3,",
        "  CXCL16–CXCR6, CCL5–CCR1/CCR5, CX3CL1–CX3CR1.",
        "- **Myeloid recruitment** (no directional thesis): CCL2/CCL7–CCR2, CXCL1/2/8–CXCR1/2,",
        "  CCL3–CCR1, CSF1–CSF1R. Direction is reported, not scored as pass/fail.",
        "- HLA–CD8 and CXCL12–CXCR4 are scored in the table and are not in the support ranks.",
        "",
        "## Honest n",
        "",
        "| gate | n | note |",
        "|---|---:|---|",
        "| Locked four | 65 | 13+21+22+9 |",
        f"| Inventory units | {len(inv)} | "
        + ", ".join(f"{k}={v}" for k, v in n_by.items())
        + " |",
        f"| Q4 eligible (n_mal≥40 and a receiver ≥20) | {len(both)} | not the ligand n |",
        "",
        "GSE123902 / GSE189357: epithelium marker-malignant (EPCAM\\|KRT8\\|KRT18\\|KRT19 > 0 and PTPRC = 0);",
        "T/NK = CD3D\\|CD3E\\|CD8A\\|NKG7\\|GNLY\\|KLRD1 > 0 and not malignant;",
        "myeloid = LYZ\\|CD68\\|CD14\\|FCGR3A\\|CSF1R\\|AIF1 > 0 and not malignant or T/NK.",
        "GSE131907 / GSE205335: author malignant, T/NK, and myeloid labels. Normal tissue in GSE205335 is out.",
        "TACSTD2 is not a gate. PVRL2→NECTIN2, JAM1→F11R, IL8→CXCL8 when the official symbol is absent.",
        "Label QC against the CellChat inventory: n_mal and n_tnk match on every locked unit.",
        "",
        "Each sender or receiver group is capped at 400 cells (seed 1337) inside LIANA.",
        f"CellPhoneDB permutations: n_perms={N_PERMS}. Expression is log1p CP10k using the full-library size.",
        "",
        "## Family tests (pre-specified; BH across this 18-test family)",
        "",
        "Δ is the within-unit mean of detected edge deltas. Positive = higher from CLDN4-high.",
        "",
        "| receiver | class | method | n | 123902 | 131907 | 205335 | 189357 | mean Δ | p_W | q_BH | observed | agrees |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|",
    ]
    for receiver in ("TNK", "MYE"):
        for lr_class in ("barrier_exclusion", "recruit_effector", "recruit_myeloid"):
            for method in ("cpdb", "conn", "liana"):
                lines.append(fam_line(receiver, lr_class, method))
    lines += [
        "",
        "## Consensus hits",
        "",
        "A hit needs ≥2 of {CellPhoneDB, Connectome, LIANA+} with BH q<0.05, the same sign,",
        "and no method significant in the opposite direction. LIANA+ ρ already uses CellPhoneDB",
        "and Connectome magnitudes, so the three columns are not independent votes.",
        "The call is corroboration of the rank aggregate by the two component magnitudes.",
        "",
        "| receiver | edge | class | supports | n | Δ CPDB | q | Δ Connectome | q | Δ LIANA+ ρ | q | z | sig |",
        "|---|---|---|---|---:|---:|---|---:|---|---:|---|---:|---|",
    ]
    if hits.empty:
        lines.append("| | | | none | | | | | | | | | |")
    else:
        for _, r in hits.iterrows():
            lines.append(hit_line(r))
    lines += [
        "",
        "## Ranked within class (top 5 by support score, including non-hits)",
        "",
        "Barrier support score = consensus z (high>low). Effector-recruitment support score = −z",
        "(low>high). Myeloid-recruitment support score = |z|.",
        "",
        "| receiver | edge | class | supports | n | Δ CPDB | q | Δ Connectome | q | Δ LIANA+ ρ | q | z | sig |",
        "|---|---|---|---|---:|---:|---|---:|---|---:|---|---:|---|",
    ]
    for receiver in ("TNK", "MYE"):
        for lr_class in RANK_CLASSES:
            lines.extend(top_lines(receiver, lr_class, 5))
    lines += [
        "",
        "Full rank: `results/tables/consensus_ranked.tsv`.",
        "Per-unit scores: `results/tables/per_patient_edges.tsv`.",
        "CellPhoneDB permutation BH (within unit × sender × receiver): `results/tables/cpdb_specificity_bh.tsv`.",
        "Figure: `results/figures/consensus_barrier_vs_recruitment.png`.",
        "",
        "## What is not claimed",
        "",
        "- This does not measure spatial exclusion, contact, or muzzling.",
        "- TACSTD2 is not used. This is not dual-high.",
        "- GSE148071, GSE127465, and the other non-concordant sets are not added.",
        "- Stouffer z ranks edges. It is not a calibrated meta-analytic p-value:",
        "  the three scores share the same cells, and LIANA+ already aggregates CellPhoneDB and Connectome.",
        "- Myeloid chemokine direction was not given a pass/fail expectation.",
        "- Cell-pooled tests are not the honest n.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "bash methods/liana_consensus_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw",
        "python3 methods/liana_consensus_concordant4_cldn4/scripts/run_consensus.py",
        "```",
        "",
        f"liana {versions['liana']}; anndata {versions['anndata']}; scipy {versions['scipy']}.",
        "GSE205335 RDS is read in R (Matrix) because it is a double-gzipped dgCMatrix.",
        "",
    ]
    (HERE / "FINDING.md").write_text("\n".join(lines))
    log(f"wrote {HERE / 'FINDING.md'}")


def self_check() -> None:
    x = np.array([0, 0, 0, 1, 2, 3, 4, 5], dtype=float)
    high, low = quartile_high_low(x)
    # n=8, q1=2, q4=6; ordinal ranks follow array order for ties
    if high.sum() != 2 or low.sum() != 2:
        raise RuntimeError(f"quartile self-check failed high {high.sum()} low {low.sum()}")
    if not (x[high].mean() > x[low].mean()):
        raise RuntimeError("quartile self-check order")


def main() -> None:
    self_check()
    import anndata
    import liana
    import scipy

    versions = {
        "liana": liana.__version__,
        "anndata": anndata.__version__,
        "scipy": scipy.__version__,
    }
    log("versions " + " ".join(f"{k}={v}" for k, v in versions.items()))
    edges = load_edges()
    units = []
    units += load_gse123902(edges)
    units += load_gse131907(edges)
    units += load_gse205335(edges)
    units += load_gse189357(edges)
    inv = qc_against_cellchat(units)
    rng = np.random.default_rng(SEED)
    per_rows: list[dict] = []
    # stable order
    units = sorted(units, key=lambda u: (u["cohort"], u["patient"]))
    for unit in units:
        per_rows.extend(score_unit(unit, edges, rng))
    if not per_rows:
        raise RuntimeError("no per-patient edge rows")
    per = pd.DataFrame(per_rows)
    summary = summarize_edges(per)
    families = summarize_families(per)
    spec = specificity_summary(per)
    # rank: class, then support score
    summary = summary.sort_values(
        ["receiver", "lr_class", "support_score"],
        ascending=[True, True, False],
        na_position="last",
    )
    tab = HERE / "results" / "tables"
    fig = HERE / "results" / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    inv.to_csv(tab / "patient_inventory.tsv", sep="\t", index=False)
    per.to_csv(tab / "per_patient_edges.tsv", sep="\t", index=False)
    summary.to_csv(tab / "consensus_ranked.tsv", sep="\t", index=False)
    families.to_csv(tab / "family_tests.tsv", sep="\t", index=False)
    spec.to_csv(tab / "cpdb_specificity_bh.tsv", sep="\t", index=False)
    hits = summary[summary["supports"] != "none"]
    hits.to_csv(tab / "consensus_hits.tsv", sep="\t", index=False)
    plot_consensus(
        summary,
        fig / "consensus_barrier_vs_recruitment.png",
        fig / "consensus_barrier_vs_recruitment.pdf",
    )
    (HERE / "results" / "session_versions.txt").write_text(
        "\n".join(f"{k}\t{v}" for k, v in versions.items()) + "\n"
    )
    write_finding(inv, summary, families, versions)
    log(
        f"DONE units {inv.shape[0]} edge-rows {per.shape[0]} "
        f"hits {hits.shape[0]}"
    )


if __name__ == "__main__":
    main()
