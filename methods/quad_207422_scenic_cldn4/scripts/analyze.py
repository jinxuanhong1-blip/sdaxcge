#!/usr/bin/env python3
"""QUAD merge CLDN4-only AUCell: patient-level IFN/MHC/TJ regulons.

Cohorts: GSE207422 + GSE131907 + GSE148071 + GSE205335 malignant cells.
A10 ELF3–CLDN4 bulk RNA is taken as given and is not a success criterion.
No TACSTD2 dual-high split. Unit of inference is the patient.
AUCell uses public priors + hardcoded program sets. Not ChIP. Not cisTarget.
"""
from __future__ import annotations

import gzip
import json
import re
import shutil
import tempfile
import traceback
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parent.parent
RES = ROOT / "tables"
FIG = ROOT / "figures"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

DATA = Path("/tmp/quad_scenic")
HOLD_OUT = {"CLDN4", "TACSTD2"}
MIN_CELLS = 20
MIN_PATIENTS = 12
AUC_THR = 0.05
NL_Q = 0.75

IFN_TFS = ["IRF1", "IRF7", "IRF9", "STAT1", "STAT2"]
MHC_TFS = ["NLRC5", "CIITA", "RFX5"]
TJ_TFS = ["ELF3", "GRHL1", "GRHL2", "KLF4", "OVOL1", "OVOL2"]
ALL_TFS = IFN_TFS + MHC_TFS + TJ_TFS
TF_AXIS = {**{t: "IFN" for t in IFN_TFS}, **{t: "MHC" for t in MHC_TFS}, **{t: "TJ" for t in TJ_TFS}}

# Expected direction for CLDN4-high minus CLDN4-low (pre-specified).
# Barrier/TJ up; IFN/MHC down. Controls null. ELF3 is A10-given, not scored.
EXPECTED = {
    "IFN": "down",
    "MHC": "down",
    "TJ": "up",
    "CTRL": "null",
    "A10": "given",
}

LINEAGE = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT7", "CDH1", "MUC1"],
    "T_NK": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "NKG7", "GNLY"],
    "B_Plasma": ["CD79A", "MS4A1", "JCHAIN", "MZB1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA", "FCN1"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
}
NORMAL_LUNG = [
    "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "NAPSA", "AGER",
    "SCGB1A1", "SCGB3A2", "FOXJ1", "CAPS",
]


def log1p_cp10k(counts: np.ndarray, total: np.ndarray) -> np.ndarray:
    tot = np.asarray(total, float)
    tot[tot <= 0] = np.nan
    return np.log1p(np.asarray(counts, float) / tot[:, None] * 1e4)


def aucell(expr: np.ndarray, member: np.ndarray, auc_threshold: float = AUC_THR) -> np.ndarray:
    """Aibar-style AUCell recovery curve. Not the R AUCell binary."""
    n_cells, n_genes = expr.shape
    max_rank = max(int(np.ceil(auc_threshold * n_genes)), 1)
    idx = np.where(member)[0]
    n_set = int(idx.size)
    if n_set == 0:
        return np.full(n_cells, np.nan)
    order = np.argsort(-expr, axis=1, kind="mergesort")
    ranks = np.empty_like(order)
    row = np.arange(n_cells)[:, None]
    ranks[row, order] = np.arange(n_genes)[None, :]
    r = ranks[:, idx].astype(np.float64)
    contrib = np.clip(max_rank - r, 0, None)
    return contrib.sum(axis=1) / (n_set * max_rank)


def stream_genes(matrix: Path, wanted: set[str]) -> tuple[list[str], np.ndarray, dict[str, np.ndarray]]:
    with gzip.open(matrix, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        cells = header.split("\t")[1:]
        n = len(cells)
        total = np.zeros(n, dtype=np.float64)
        kept: dict[str, np.ndarray] = {}
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = raw[:tab].decode("ascii").split(".")[0]
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            total += arr
            if gene in wanted:
                kept[gene] = arr.copy()
            if n_genes % 4000 == 0:
                print(f"  streamed {n_genes} genes, kept {len(kept)}", flush=True)
    print(f"[stream] {matrix.name}: cells={n} genes={n_genes} kept={len(kept)}", flush=True)
    return cells, total, kept


def load_priors() -> pd.DataFrame:
    path = ROOT / "resources" / "tf_targets_public.tsv"
    if not path.exists():
        raise SystemExit("priors missing; run scripts/download_priors.py")
    df = pd.read_csv(path, sep="\t")
    df["tf"] = df["tf"].astype(str)
    df["target"] = df["target"].astype(str)
    return df


def load_program_sets() -> dict[str, list[str]]:
    return json.loads((ROOT / "resources" / "gene_sets.json").read_text())


def build_regulons(priors: pd.DataFrame, programs: dict[str, list[str]], present: set[str]) -> pd.DataFrame:
    rows = []
    for tf in ALL_TFS:
        hits = sorted(
            {
                g
                for g in priors.loc[priors.tf == tf, "target"]
                if g in present and g != tf and g not in HOLD_OUT
            }
        )
        axis = TF_AXIS[tf]
        kind = "public_prior"
        note = "TRRUST+DoRothEA+CollecTRI union; CLDN4/TACSTD2 held out; not ChIP"
        if tf == "ELF3":
            kind = "a10_given_prior"
            note = "A10 ELF3–CLDN4 taken as given; scored but not a success criterion; not ChIP"
            axis = "A10"
        rows.append(
            dict(
                regulon=f"{tf}_prior",
                tf=tf,
                axis=axis,
                kind=kind,
                n_targets=len(hits),
                targets=";".join(hits),
                expected=EXPECTED[axis],
                note=note,
            )
        )
    # axis unions (ELF3 excluded from TJ union — A10 given)
    for axis, tfs in (("IFN", IFN_TFS), ("MHC", MHC_TFS), ("TJ", [t for t in TJ_TFS if t != "ELF3"])):
        members = []
        for tf in tfs:
            rec = next(r for r in rows if r["regulon"] == f"{tf}_prior")
            members.extend(rec["targets"].split(";") if rec["targets"] else [])
        members = sorted({g for g in members if g})
        rows.append(
            dict(
                regulon=f"{axis}_prior_union",
                tf="+".join(tfs),
                axis=axis,
                kind="combinatorial_public_prior",
                n_targets=len(members),
                targets=";".join(members),
                expected=EXPECTED[axis],
                note="union of axis TF priors; CLDN4/TACSTD2 held out; not ChIP",
            )
        )
    prog_meta = {
        "ISG_CORE": ("IFN", "program_ISG", "down", "literature ISG core; not a TF regulon"),
        "MHC1_APM": ("MHC", "program_APM", "down", "MHC-I antigen processing/presentation; not a TF regulon"),
        "TJ_STRUCT": ("TJ", "program_TJ", "up", "structural tight-junction genes; CLDN4 held out"),
        "CTRL_OXPHOS": ("CTRL", "program_control", "null", "OXPHOS control; expected null"),
        "CTRL_RIBO": ("CTRL", "program_control", "null", "ribosome control; expected null"),
    }
    for name, (axis, kind, exp, note) in prog_meta.items():
        genes = [g for g in programs.get(name, []) if g in present and g not in HOLD_OUT]
        # unique preserve order
        seen = set()
        uniq = []
        for g in genes:
            if g not in seen:
                seen.add(g)
                uniq.append(g)
        rows.append(
            dict(
                regulon=name,
                tf="",
                axis=axis,
                kind=kind,
                n_targets=len(uniq),
                targets=";".join(uniq),
                expected=exp,
                note=note,
            )
        )
    return pd.DataFrame(rows)


def wanted_genes(priors: pd.DataFrame, programs: dict) -> set[str]:
    genes = set(HOLD_OUT) | set(ALL_TFS) | set(NORMAL_LUNG)
    for lst in LINEAGE.values():
        genes.update(lst)
    genes.update(priors.target.astype(str))
    genes.update(priors.tf.astype(str))
    for k, v in programs.items():
        if k == "note":
            continue
        genes.update(v)
    return genes


def score_cells(logX: np.ndarray, genes: list[str], regulons: pd.DataFrame) -> pd.DataFrame:
    gi = {g: i for i, g in enumerate(genes)}
    out = {}
    if "CLDN4" in gi:
        out["CLDN4"] = logX[:, gi["CLDN4"]]
    for rec in regulons.itertuples(index=False):
        members = [g for g in str(rec.targets).split(";") if g and g in gi]
        mask = np.zeros(len(genes), dtype=bool)
        for g in members:
            mask[gi[g]] = True
        out[rec.regulon] = aucell(logX, mask)
    return pd.DataFrame(out)


def zscore(x: pd.Series) -> pd.Series:
    s = float(x.std(ddof=1))
    if not np.isfinite(s) or s == 0:
        return pd.Series(np.nan, index=x.index)
    return (x - x.mean()) / s


def wilcox(high: np.ndarray, low: np.ndarray) -> dict:
    high = np.asarray(high, float)
    low = np.asarray(low, float)
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    rec = dict(
        n_high=int(high.size),
        n_low=int(low.size),
        mean_high=float(np.mean(high)) if high.size else np.nan,
        mean_low=float(np.mean(low)) if low.size else np.nan,
        delta=np.nan,
        U=np.nan,
        p=np.nan,
    )
    if high.size and low.size:
        rec["delta"] = rec["mean_high"] - rec["mean_low"]
    if high.size >= 3 and low.size >= 3:
        U, p = stats.mannwhitneyu(high, low, alternative="two-sided")
        rec["U"] = float(U)
        rec["p"] = float(p)
    return rec


def spear(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 6 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return dict(rho=np.nan, p=np.nan, n=n)
    r, p = stats.spearmanr(x[m], y[m])
    return dict(rho=float(r), p=float(p), n=n)


# ---------------------------------------------------------------------------
# Dataset loaders
# ---------------------------------------------------------------------------

def load_gse207422(wanted: set[str]) -> pd.DataFrame:
    meta = pd.read_excel(DATA / "GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    meta = meta.dropna(subset=["Sample", "Patient"]).copy()
    meta["Sample"] = meta["Sample"].astype(str)
    meta["Patient"] = meta["Patient"].astype(str)
    cells, total, kept = stream_genes(
        DATA / "GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz", wanted
    )
    sample = pd.Series(cells).str.replace(r"_\d+$", "", regex=True)
    genes = sorted(kept)
    X = np.vstack([kept[g] for g in genes]).T
    logX = log1p_cp10k(X, total)
    gi = {g: i for i, g in enumerate(genes)}

    def mean_set(names):
        idx = [gi[g] for g in names if g in gi]
        if not idx:
            return np.zeros(logX.shape[0])
        return logX[:, idx].mean(axis=1)

    scores = {name: mean_set(lst) for name, lst in LINEAGE.items()}
    lin = pd.DataFrame(scores)
    lineage = lin.idxmax(axis=1)
    nl = mean_set(NORMAL_LUNG)
    epi = lineage.eq("Epithelial").to_numpy()
    nl_cut = float(np.quantile(nl[epi], NL_Q)) if epi.sum() else np.inf
    mal = epi & (nl <= nl_cut) & (total >= 200)
    df = pd.DataFrame(
        {
            "dataset": "GSE207422",
            "patient": sample.map(dict(zip(meta.Sample, meta.Patient))),
            "sample": sample,
            "malignant": mal,
        }
    )
    df = df.join(pd.DataFrame(logX, columns=genes))
    keep = df.malignant & df.patient.notna()
    print(f"[GSE207422] epi={int(epi.sum())} mal={int(mal.sum())} patients={df.loc[keep,'patient'].nunique()}", flush=True)
    return df.loc[keep].reset_index(drop=True)


def load_gse131907(wanted: set[str]) -> pd.DataFrame:
    ann = pd.read_csv(DATA / "GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    cells, total, kept = stream_genes(
        DATA / "GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", wanted
    )
    genes = sorted(kept)
    X = np.vstack([kept[g] for g in genes]).T
    logX = log1p_cp10k(X, total)
    meta = pd.DataFrame({"Index": cells})
    meta = meta.merge(ann, on="Index", how="left")
    tumor = meta.Sample_Origin.isin({"tLung", "tL/B", "mLN", "mBrain", "PE"})
    mal = tumor & meta.Cell_type.eq("Epithelial cells") & (total >= 200)
    pid = meta.Sample.astype(str).str.extract(r"(\d+)$")[0]
    df = pd.DataFrame(
        {
            "dataset": "GSE131907",
            "patient": pid,
            "sample": meta.Sample.astype(str),
            "malignant": mal.to_numpy(),
        }
    )
    df = df.join(pd.DataFrame(logX, columns=genes))
    keep = df.malignant & df.patient.notna()
    print(f"[GSE131907] tumor-epi={int(mal.sum())} patients={df.loc[keep,'patient'].nunique()}", flush=True)
    return df.loc[keep].reset_index(drop=True)


def load_gse148071(wanted: set[str]) -> pd.DataFrame:
    import h5py
    import scipy.sparse as sp

    meta = pd.read_csv(DATA / "GSE148071/NSCLC_GSE148071_CellMetainfo_table.tsv", sep="\t")
    path = DATA / "GSE148071/NSCLC_GSE148071_expression.h5"
    with h5py.File(path, "r") as f:
        grp = f
        if not all(k in f for k in ("data", "indices", "indptr")):
            for key in f.keys():
                g = f[key]
                if hasattr(g, "keys") and all(k in g for k in ("data", "indices", "indptr")):
                    grp = g
                    break

        def _as_str(arr):
            return np.array(
                [x.decode() if isinstance(x, bytes) else str(x) for x in arr],
                dtype=object,
            )

        def _names(cands):
            for c in cands:
                node = grp[c] if c in grp else (f[c] if c in f else None)
                if node is None:
                    continue
                if hasattr(node, "keys"):
                    for sub in ("name", "id", "gene_names"):
                        if sub in node:
                            return _as_str(node[sub][:])
                    continue
                if hasattr(node, "shape"):
                    return _as_str(node[:])
            return None

        genes_all = _names(["features", "gene_names", "genes", "rownames"])
        barcodes = _names(["barcodes", "cell_names", "colnames"])
        if genes_all is None or barcodes is None:
            raise ValueError(
                f"GSE148071 h5 missing names genes={genes_all is None} "
                f"barcodes={barcodes is None} keys={list(f.keys())}"
            )
        shape = tuple(int(x) for x in grp["shape"][:]) if "shape" in grp else None
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
    n0, n1 = int(shape[0]), int(shape[1])
    if genes_all is not None and barcodes is not None and len(genes_all) == n1 and len(barcodes) == n0:
        n_genes, n_cells = n1, n0
        csc = sp.csc_matrix((data, indices, indptr), shape=(n_cells, n_genes))
        def col(i):
            return np.asarray(csc[:, i].todense()).ravel()
    else:
        n_genes, n_cells = n0, n1
        csc = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))
        def col(i):
            return np.asarray(csc.getrow(i).todense()).ravel()
    present = [g for g in sorted(wanted) if g in set(genes_all)]
    gix = {g: int(np.where(genes_all == g)[0][0]) for g in present}
    expr = {g: col(i) for g, i in gix.items()}
    del csc, data, indices, indptr
    # TISCH values are already log2(TPM/10+1); AUCell ranks them as-is.
    meta = meta.copy()
    meta["Cell"] = meta["Cell"].astype(str)
    bar = pd.Series(barcodes).astype(str)
    # barcodes may or may not match Cell
    if set(bar) & set(meta.Cell):
        order = pd.Index(bar)
        meta = meta.set_index("Cell").reindex(order)
    else:
        meta = meta.reset_index(drop=True)
        if len(meta) != n_cells:
            raise ValueError(f"GSE148071 meta {len(meta)} != cells {n_cells}")
    mal = meta["Celltype (major-lineage)"].astype(str).str.strip().eq("Malignant")
    df = pd.DataFrame(
        {
            "dataset": "GSE148071",
            "patient": meta["Patient"].astype(str).to_numpy(),
            "sample": meta["Sample"].astype(str).to_numpy(),
            "malignant": mal.to_numpy(),
        }
    )
    df = df.join(pd.DataFrame(expr))
    keep = df.malignant & df.patient.notna() & ~df.patient.isin(["nan", "None"])
    print(f"[GSE148071] mal={int(mal.sum())} patients={df.loc[keep,'patient'].nunique()}", flush=True)
    return df.loc[keep].reset_index(drop=True)


def parse_205335_soft(path: Path) -> pd.DataFrame:
    records = []
    current = None
    desc, titles = [], []
    with gzip.open(path, "rt", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = desc[0] if desc else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                desc, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                desc.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    k, v = value.split(": ", 1)
                    current[k] = v
        if current is not None:
            current["description"] = desc[0] if desc else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    md = pd.DataFrame(records)
    read_end = md["platform"].str.extract(r"Single Cell ([35])'")[0]
    md["orig.ident"] = md["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    return md.rename(columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"})


def load_gse205335(wanted: set[str]) -> pd.DataFrame:
    import rdata
    from scipy import sparse

    print("[GSE205335] reading identities + RDS (slow)", flush=True)
    ident = pd.read_csv(DATA / "GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    soft = parse_205335_soft(DATA / "GSE205335/GSE205335_family.soft.gz")
    path = DATA / "GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz"
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = Path(tmp) / path.stem
        with gzip.open(path, "rb") as src, matrix_path.open("wb") as dest:
            shutil.copyfileobj(src, dest, 16 * 1024 * 1024)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    genes_all = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    total = np.asarray(matrix.sum(axis=0)).ravel()
    present = [g for g in sorted(wanted) if g in set(genes_all)]
    expr = {}
    for g in present:
        pos = np.flatnonzero(genes_all == g)
        if len(pos) != 1:
            continue
        expr[g] = np.asarray(matrix.getrow(int(pos[0])).toarray()).ravel()
    del matrix, obj
    genes = list(expr)
    X = np.vstack([expr[g] for g in genes]).T
    logX = log1p_cp10k(X, total)
    cells = pd.DataFrame({"barcode": barcodes, "total_umi": total})
    cells = cells.merge(ident, on="barcode", how="left")
    cells = cells.merge(
        soft[["orig.ident", "patient", "tissue", "cancer_subtype"]],
        on="orig.ident",
        how="left",
    )
    tumor_ok = ~cells.tissue.isin(["Normal LN", "Normal Lung", "Normal Brain"])
    nsclc = cells.cancer_subtype.isin(["ADC", "SQ"])
    mal = cells["lineage.sub"].eq("Malignant cells") & tumor_ok & nsclc & (total >= 200)
    df = pd.DataFrame(
        {
            "dataset": "GSE205335",
            "patient": cells["patient"].astype(str),
            "sample": cells["orig.ident"].astype(str),
            "malignant": mal.to_numpy(),
        }
    )
    df = df.join(pd.DataFrame(logX, columns=genes))
    keep = df.malignant & df.patient.notna() & ~df.patient.isin(["nan", "None"])
    print(f"[GSE205335] mal-NSCLC={int(mal.sum())} patients={df.loc[keep,'patient'].nunique()}", flush=True)
    return df.loc[keep].reset_index(drop=True)


def patient_table(cells: pd.DataFrame, regulons: pd.DataFrame) -> pd.DataFrame:
    score_cols = ["CLDN4"] + [r for r in regulons.regulon if r in cells.columns]
    rows = []
    for (ds, pat), sub in cells.groupby(["dataset", "patient"], observed=True):
        if len(sub) < MIN_CELLS:
            continue
        rec = dict(dataset=ds, patient=str(pat), n_malignant_cells=int(len(sub)))
        for c in score_cols:
            rec[c] = float(sub[c].mean())
        rows.append(rec)
    return pd.DataFrame(rows)


def split_within_cohort(pat: pd.DataFrame) -> pd.DataFrame:
    out = []
    for ds, sub in pat.groupby("dataset"):
        sub = sub.copy()
        med = float(sub.CLDN4.median())
        sub["cldn4_split"] = np.where(sub.CLDN4 >= med, "high", "low")
        sub["cldn4_median"] = med
        # ties at median: if unbalanced, leave as high (>=)
        out.append(sub)
    return pd.concat(out, ignore_index=True)


def run_tests(pat: pd.DataFrame, regulons: pd.DataFrame) -> pd.DataFrame:
    rows = []
    regs = [r for r in regulons.regulon if r in pat.columns]
    # per-cohort
    for ds, sub in pat.groupby("dataset"):
        high = sub[sub.cldn4_split == "high"]
        low = sub[sub.cldn4_split == "low"]
        for reg in regs:
            axis = regulons.loc[regulons.regulon == reg, "axis"].iloc[0]
            expected = regulons.loc[regulons.regulon == reg, "expected"].iloc[0]
            w = wilcox(high[reg].to_numpy(), low[reg].to_numpy())
            s = spear(sub.CLDN4, sub[reg])
            rows.append(
                dict(
                    scope="cohort",
                    dataset=ds,
                    regulon=reg,
                    axis=axis,
                    expected=expected,
                    **w,
                    spearman_rho=s["rho"],
                    spearman_p=s["p"],
                    n_patients=int(len(sub)),
                )
            )
    # pooled on within-cohort z-scores
    z = pat.copy()
    for col in ["CLDN4"] + regs:
        z[col] = z.groupby("dataset")[col].transform(zscore)
    high = z[z.cldn4_split == "high"]
    low = z[z.cldn4_split == "low"]
    for reg in regs:
        axis = regulons.loc[regulons.regulon == reg, "axis"].iloc[0]
        expected = regulons.loc[regulons.regulon == reg, "expected"].iloc[0]
        w = wilcox(high[reg].to_numpy(), low[reg].to_numpy())
        s = spear(z.CLDN4, z[reg])
        rows.append(
            dict(
                scope="pooled_z_within_cohort",
                dataset="QUAD",
                regulon=reg,
                axis=axis,
                expected=expected,
                **w,
                spearman_rho=s["rho"],
                spearman_p=s["p"],
                n_patients=int(z.CLDN4.notna().sum()),
            )
        )
    return pd.DataFrame(rows)


def call_row(row: pd.Series) -> str:
    if row.expected == "given":
        return "A10_GIVEN"
    n = int(row.n_patients) if np.isfinite(row.n_patients) else 0
    if n < MIN_PATIENTS or not np.isfinite(row.p):
        return "UNDERPOWERED"
    delta = row.delta
    p = row.p
    if row.expected == "null":
        return "NULL_OK" if p >= 0.05 else "CTRL_MOVED"
    if p >= 0.05:
        return "NS"
    if row.expected == "down" and delta < 0:
        return "SUPPORT"
    if row.expected == "up" and delta > 0:
        return "SUPPORT"
    return "WRONG_SIGN"


def verdict(tests: pd.DataFrame) -> dict:
    prim = tests[tests.scope == "pooled_z_within_cohort"].copy()
    prim = prim[~prim.regulon.eq("ELF3_prior")]
    # primary regulons: axis unions + programs
    primary = {
        "IFN_prior_union",
        "MHC_prior_union",
        "TJ_prior_union",
        "ISG_CORE",
        "MHC1_APM",
        "TJ_STRUCT",
    }
    sub = prim[prim.regulon.isin(primary)].copy()
    sub["call"] = sub.apply(call_row, axis=1)
    ifn = sub[sub.axis == "IFN"]
    mhc = sub[sub.axis == "MHC"]
    tj = sub[sub.axis == "TJ"]
    ifn_ok = bool((ifn.call == "SUPPORT").any())
    mhc_ok = bool((mhc.call == "SUPPORT").any())
    tj_ok = bool((tj.call == "SUPPORT").any())
    n = int(prim.n_patients.max()) if len(prim) else 0
    if n < MIN_PATIENTS:
        overall = "UNDERPOWERED"
    elif (ifn_ok or mhc_ok) and tj_ok:
        overall = "SUPPORT"
    elif ifn_ok or mhc_ok or tj_ok:
        overall = "PARTIAL"
    else:
        overall = "NOT_SUPPORTED"
    return {
        "overall": overall,
        "ifn_support": ifn_ok,
        "mhc_support": mhc_ok,
        "tj_support": tj_ok,
        "n_patients_pooled": n,
        "primary_calls": sub[["regulon", "axis", "expected", "n_high", "n_low", "delta", "p", "spearman_rho", "spearman_p", "call"]].to_dict("records"),
    }


def write_finding(pat: pd.DataFrame, tests: pd.DataFrame, regulons: pd.DataFrame, verd: dict, counts: dict) -> None:
    tests = tests.copy()
    tests["call"] = tests.apply(call_row, axis=1)
    prim = tests[tests.scope == "pooled_z_within_cohort"]
    n_pat = int(pat.groupby(["dataset", "patient"]).ngroups)
    n_by = pat.groupby("dataset").size().to_dict()
    cells_by = pat.groupby("dataset")["n_malignant_cells"].sum().to_dict()
    high_low = pat.groupby(["dataset", "cldn4_split"]).size().unstack(fill_value=0)

    def fmt_p(p):
        if not np.isfinite(p):
            return "NA"
        return f"{p:.2e}" if p < 0.001 else f"{p:.3f}"

    def fmt_d(x):
        return "NA" if not np.isfinite(x) else f"{x:.3f}"

    lines = []
    lines.append("# FINDING — QUAD CLDN4-only AUCell (IFN / MHC / TJ regulons)")
    lines.append("")
    lines.append(f"**Verdict: `{verd['overall']}`.** Patient-level CLDN4-high vs low on malignant cells from **GSE207422 + GSE131907 + GSE148071 + GSE205335** (207422 included). A10 ELF3–CLDN4 is taken as given and is not a success criterion. No TACSTD2 dual-high split.")
    lines.append("")
    lines.append(f"Pooled n = **{verd['n_patients_pooled']} patients** (within-cohort z-score of patient means; min {MIN_CELLS} malignant cells). IFN support={verd['ifn_support']}; MHC support={verd['mhc_support']}; TJ support={verd['tj_support']}.")
    lines.append("")
    lines.append("Primary unit is the **patient**. Cell counts below are inventory only.")
    lines.append("")
    lines.append("## Pre-specified rules")
    lines.append("")
    lines.append("- Split: within-cohort median of patient-mean malignant CLDN4 (high ≥ median).")
    lines.append("- Pool: z-score patient means within cohort, then Wilcoxon high vs low and Spearman CLDN4 vs regulon.")
    lines.append("- Expected: TJ up in CLDN4-high; IFN/MHC down; OXPHOS/ribo controls null.")
    lines.append("- SUPPORT on an axis: expected direction and p<0.05 on a primary regulon (axis prior-union or program set).")
    lines.append("- Overall SUPPORT: (IFN or MHC) **and** TJ. PARTIAL: one axis. NOT_SUPPORTED: none. UNDERPOWERED: n<12.")
    lines.append("- ELF3_prior is A10-given and excluded from the verdict.")
    lines.append("- AUCell = Aibar recovery curve on public priors / hardcoded programs. Not cisTarget. Not ChIP.")
    lines.append("")
    lines.append("## Honest n")
    lines.append("")
    lines.append("| cohort | patients (≥20 mal. cells) | malignant cells | CLDN4-high | CLDN4-low | malignant definition |")
    lines.append("| --- | ---: | ---: | ---: | ---: | --- |")
    defs = {
        "GSE207422": "marker epithelial AND normal-lung score ≤ p75 (CopyKAT not on GEO)",
        "GSE131907": "author Epithelial cells in tumor/met tissue (tLung/tL/B/mLN/mBrain/PE); no CopyKAT",
        "GSE148071": "TISCH2 major-lineage = Malignant (Wu et al. 2021 via TISCH2)",
        "GSE205335": "author lineage.sub = Malignant cells; ADC+SQ; drop normal tissues",
    }
    for ds in ["GSE207422", "GSE131907", "GSE148071", "GSE205335"]:
        nh = int(high_low.loc[ds, "high"]) if ds in high_low.index and "high" in high_low.columns else 0
        nl = int(high_low.loc[ds, "low"]) if ds in high_low.index and "low" in high_low.columns else 0
        lines.append(
            f"| {ds} | {n_by.get(ds, 0)} | {int(cells_by.get(ds, 0))} | {nh} | {nl} | {defs[ds]} |"
        )
    lines.append(f"| **QUAD pool** | **{n_pat}** | **{int(pat.n_malignant_cells.sum())}** | {int((pat.cldn4_split=='high').sum())} | {int((pat.cldn4_split=='low').sum())} | within-cohort z |")
    lines.append("")
    if counts.get("dropped"):
        lines.append("Dropped (honest): " + "; ".join(counts["dropped"]))
        lines.append("")
    lines.append("## Primary regulon table (pooled, within-cohort z)")
    lines.append("")
    lines.append("| regulon | axis | expected | n_high | n_low | Δ (high−low) | Wilcoxon p | Spearman ρ (p) | call |")
    lines.append("| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |")
    primary = [
        "IFN_prior_union", "ISG_CORE", "IRF1_prior", "STAT1_prior",
        "MHC_prior_union", "MHC1_APM", "NLRC5_prior", "CIITA_prior",
        "TJ_prior_union", "TJ_STRUCT", "GRHL2_prior", "KLF4_prior", "OVOL1_prior",
        "ELF3_prior", "CTRL_OXPHOS", "CTRL_RIBO",
    ]
    shown = set()
    for reg in primary:
        hit = prim[prim.regulon == reg]
        if hit.empty:
            continue
        r = hit.iloc[0]
        shown.add(reg)
        lines.append(
            f"| {reg} | {r.axis} | {r.expected} | {int(r.n_high)} | {int(r.n_low)} | {fmt_d(r.delta)} | {fmt_p(r.p)} | {fmt_d(r.spearman_rho)} ({fmt_p(r.spearman_p)}) | {r.call} |"
        )
    lines.append("")
    lines.append("Full tests: `tables/tests.tsv`. Regulon gene lists: `tables/regulons.tsv`. Patient scores: `tables/patient_scores.tsv`.")
    lines.append("")
    lines.append("## What this is / is not")
    lines.append("")
    lines.append("- **Is** additive CLDN4-only AUCell on a four-cohort malignant merge, patient-level, 207422 included.")
    lines.append("- **Is not** a Harmony/scVI joint embedding (patients are z-scored within cohort; no cell mixing).")
    lines.append("- **Is not** pySCENIC cisTarget or a binding map. Priors are TRRUST/DoRothEA/CollecTRI.")
    lines.append("- **Is not** a re-test of A10 ELF3–CLDN4. ELF3_prior is reported and tagged `A10_GIVEN`.")
    lines.append("- **Is not** a TACSTD2 / dual-high analysis.")
    lines.append("- GSE131907 malignant is tumor-tissue epithelium (author label), not CopyKAT.")
    lines.append("- GSE148071 expression is TISCH2 `log2(TPM/10+1)`; others are `log1p(CP10k)` from UMI. Ranking AUCell is within-matrix; pooling uses within-cohort z.")
    lines.append("- GSE205335 SCLC/NUT and normal tissues were excluded from the primary NSCLC malignant set.")
    lines.append("")
    lines.append("## Regulon inventory")
    lines.append("")
    lines.append("| regulon | axis | kind | n_targets | expected |")
    lines.append("| --- | --- | --- | ---: | --- |")
    for r in regulons.itertuples(index=False):
        lines.append(f"| {r.regulon} | {r.axis} | {r.kind} | {int(r.n_targets)} | {r.expected} |")
    lines.append("")
    (ROOT / "FINDING.md").write_text("\n".join(lines) + "\n")
    print("wrote", ROOT / "FINDING.md")


def make_figures(pat: pd.DataFrame, tests: pd.DataFrame) -> None:
    prim = tests[tests.scope == "pooled_z_within_cohort"].copy()
    prim["call"] = prim.apply(call_row, axis=1)
    keep = [
        "IFN_prior_union", "ISG_CORE", "MHC_prior_union", "MHC1_APM",
        "TJ_prior_union", "TJ_STRUCT", "ELF3_prior", "CTRL_OXPHOS",
    ]
    sub = prim[prim.regulon.isin(keep)].copy()
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 4.2))
    y = np.arange(len(sub))
    colors = []
    for c in sub.call:
        colors.append({"SUPPORT": "#2c7bb6", "WRONG_SIGN": "#d7191c", "NS": "#888888", "A10_GIVEN": "#fdae61", "NULL_OK": "#abdda4", "CTRL_MOVED": "#f46d43"}.get(c, "#bbbbbb"))
    ax.barh(y, sub.delta.fillna(0), color=colors, edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(sub.regulon)
    ax.set_xlabel("CLDN4-high − low (within-cohort z of patient-mean AUCell)")
    ax.set_title(f"QUAD patient-level AUCell  n={int(sub.n_patients.iloc[0])}")
    fig.tight_layout()
    fig.savefig(FIG / "fig_pooled_delta.png", dpi=140)
    fig.savefig(FIG / "fig_pooled_delta.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    for ds, g in pat.groupby("dataset"):
        ax.scatter(g.CLDN4, g.get("ISG_CORE", g.CLDN4 * np.nan), s=22, alpha=0.75, label=f"{ds} n={len(g)}")
    ax.set_xlabel("patient-mean malignant CLDN4")
    ax.set_ylabel("patient-mean ISG_CORE AUCell")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("CLDN4 vs ISG (raw patient means; not z)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_cldn4_vs_isg.png", dpi=140)
    fig.savefig(FIG / "fig_cldn4_vs_isg.pdf")
    plt.close(fig)


def main() -> None:
    priors = load_priors()
    programs = load_program_sets()
    wanted = wanted_genes(priors, programs)
    print(f"[wanted] {len(wanted)} genes", flush=True)

    frames = []
    dropped = []
    cache_dir = DATA / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    loaders = [
        ("GSE207422", load_gse207422),
        ("GSE131907", load_gse131907),
        ("GSE148071", load_gse148071),
        ("GSE205335", load_gse205335),
    ]
    for name, fn in loaders:
        try:
            ck = cache_dir / f"{name}.pkl"
            if ck.exists() and ck.stat().st_size > 0:
                print(f"[cache] {name} {ck}", flush=True)
                frames.append(pd.read_pickle(ck))
                continue
            df = fn(wanted)
            df.to_pickle(ck)
            frames.append(df)
        except Exception as exc:
            traceback.print_exc()
            dropped.append(f"{name} failed: {exc}")
            print(f"[FAIL] {name}: {exc}", flush=True)

    if not frames:
        raise SystemExit("no cohorts loaded")
    # common genes
    gene_cols = None
    for fr in frames:
        cols = {c for c in fr.columns if c not in {"dataset", "patient", "sample", "malignant"}}
        gene_cols = cols if gene_cols is None else gene_cols & cols
    present = set(gene_cols or [])
    print(f"[common genes] {len(present)}", flush=True)
    regulons = build_regulons(priors, programs, present)
    regulons.to_csv(RES / "regulons.tsv", sep="\t", index=False)

    scored = []
    gene_list = sorted(present)
    print(f"[score] {len(gene_list)} common genes; {len(regulons)} regulons", flush=True)
    for fr in frames:
        ds = str(fr.dataset.iloc[0])
        print(f"[score] {ds} cells={len(fr)} ...", flush=True)
        logX = fr.loc[:, gene_list].to_numpy(dtype=np.float32)
        sc = score_cells(logX, gene_list, regulons)
        sc.insert(0, "dataset", fr["dataset"].to_numpy())
        sc.insert(1, "patient", fr["patient"].to_numpy())
        sc.insert(2, "sample", fr["sample"].to_numpy())
        scored.append(sc)
        print(f"[score] {ds} done", flush=True)
    cells = pd.concat(scored, ignore_index=True)
    pat = patient_table(cells, regulons)
    pat = split_within_cohort(pat)
    pat.to_csv(RES / "patient_scores.tsv", sep="\t", index=False)
    tests = run_tests(pat, regulons)
    tests["call"] = tests.apply(call_row, axis=1)
    tests.to_csv(RES / "tests.tsv", sep="\t", index=False)
    verd = verdict(tests)
    counts = {
        "n_patients": int(len(pat)),
        "n_by_dataset": pat.groupby("dataset").size().to_dict(),
        "n_cells_by_dataset": pat.groupby("dataset")["n_malignant_cells"].sum().astype(int).to_dict(),
        "dropped": dropped,
    }
    summary = {
        "question": "Do patient-level CLDN4-high malignant cells differ in IFN/MHC/TJ AUCell on the QUAD merge?",
        "a10_taken_as_given": True,
        "dual_high": False,
        "include_GSE207422": True,
        "unit": "patient (malignant cells, min 20)",
        "aucell": "Aibar recovery curve; priors + program sets; not cisTarget",
        "verdict": verd,
        "counts": counts,
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    make_figures(pat, tests)
    write_finding(pat, tests, regulons, verd, counts)
    print(json.dumps({"verdict": verd["overall"], "n": verd["n_patients_pooled"], "dropped": dropped}, indent=2))


if __name__ == "__main__":
    main()
