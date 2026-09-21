"""Load concordant-4 units for the barrier-family sweep.

CLDN4 only. No dual-high. No GSE148071.
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


def qc_against_cellchat(units: list[dict], require_all: bool = True) -> pd.DataFrame:
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
    if require_all and missing:
        raise RuntimeError(f"missing locked units: {sorted(missing)[:8]}")
    return inv


def gene_vector(unit: dict, gene: str) -> np.ndarray | None:
    if gene not in unit["genes"]:
        return None
    return unit["mat"][:, unit["genes"].index(gene)]

