#!/usr/bin/env python3
"""GSE205335 malignant CLDN4 Q4 vs same-patient T/NK + CellChat-style ligands.

ADDITIVE. CLDN4 only. Patient is the unit.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parent
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_CELLS = 20
MIN_DETECT_ARM = 3

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
ATTACK_LIGANDS = {"IFNG", "TNF", "FASLG", "TNFSF10", "LTA"}
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD8A", "CD3E", "IFNG", "TNF"]


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_num(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:+.{digits}f}" if value < 0 or value > 0 else f"{value:.{digits}f}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def q4_vs_q1(cldn4, immune) -> dict | None:
    """Within-cohort CLDN4 quartiles; MWU on immune. r_rb < 0 = Q4 immune lower."""
    frame = pd.DataFrame(
        {"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)}
    )
    frame = frame[np.isfinite(frame["c"]) & np.isfinite(frame["i"])].copy()
    n = int(len(frame))
    if n < 6:
        return None
    ranks = frame["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return None
    if qs.nunique() < 4:
        return None
    q1 = frame.loc[qs == "Q1", "i"]
    q4 = frame.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return None
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    return {
        "n": n,
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "median_q1": float(q1.median()),
        "median_q4": float(q4.median()),
        "delta_median": float(q4.median() - q1.median()),
        "mwu_u": float(u),
        "p": float(p),
        "r_rb": float(r_rb),
        "thin": n < 8 or n1 < 3 or n4 < 3,
    }


def assign_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def spearman(x, y) -> tuple[float, float, int]:
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    mask = np.isfinite(xa) & np.isfinite(ya)
    xa, ya = xa[mask], ya[mask]
    n = int(xa.size)
    if n < 3:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(xa, ya)
    return float(rho), float(p), n


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_lr(db_dir: Path, matrix_genes: set[str]) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    inter = inter.rename(
        columns={"ligand.symbol": "ligand_symbol", "receptor.symbol": "receptor_symbol"}
    )
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


def ligand_class(genes) -> str:
    parts = set(genes) if not isinstance(genes, str) else set(genes.split("|"))
    tags = []
    if parts & BARRIER_LIGANDS:
        tags.append("barrier")
    if parts & INHIB_LIGANDS:
        tags.append("inhibitory")
    if parts & RECRUIT_LIGANDS:
        tags.append("recruit")
    if parts & ATTACK_LIGANDS:
        tags.append("attack")
    return "|".join(tags) if tags else "other"


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def load_selected_genes(path: Path, wanted: set[str]) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, list[str]]:
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            print(f"decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("read RDS", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        positions = np.flatnonzero(genes == gene)
        if len(positions) != 1:
            continue
        extracted[gene] = np.asarray(matrix.getrow(int(positions[0])).toarray()).ravel()
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    print(f"extracted {len(extracted)} / {len(wanted)} genes", flush=True)
    return extracted, library_umi, barcodes, genes.tolist()


def trim_mean_1d(values: np.ndarray, proportiontocut: float = TRIM) -> float:
    n = int(values.size)
    if n == 0:
        return 0.0
    if n == 1:
        return float(values[0])
    k = int(n * proportiontocut)
    if k == 0:
        return float(values.mean())
    s = np.sort(values)
    return float(s[k : n - k].mean())


def geom_mean(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return 0.0
    if np.any(arr <= 0):
        return 0.0
    if arr.size == 1:
        return float(arr[0])
    return float(np.exp(np.mean(np.log(arr))))


def hill_prob(lig: float, rec: float, kh: float = KH) -> float:
    if lig <= 0 or rec <= 0:
        return 0.0
    prod = lig * rec
    return float(prod / (kh + prod))


def complex_stats(log_cp: dict[str, np.ndarray], pos: dict[str, np.ndarray], idx: np.ndarray, subunits: tuple[str, ...]):
    means = []
    props = []
    for gene in subunits:
        if gene not in log_cp or idx.size == 0:
            return 0.0, 0.0
        means.append(trim_mean_1d(log_cp[gene][idx]))
        props.append(float(pos[gene][idx].mean()) if idx.size else 0.0)
    return geom_mean(means), float(min(props)) if props else 0.0


def score_patient_pairs(lr: pd.DataFrame, log_cp, pos, mal_idx, tnk_idx) -> list[dict]:
    rows = []
    if mal_idx.size < MIN_CELLS or tnk_idx.size < MIN_CELLS:
        return rows
    cache: dict[tuple[str, tuple[str, ...]], tuple[float, float]] = {}

    def cached(kind: str, subunits: tuple[str, ...], idx: np.ndarray):
        key = (kind, subunits)
        hit = cache.get(key)
        if hit is None:
            hit = complex_stats(log_cp, pos, idx, subunits)
            cache[key] = hit
        return hit

    for rec in lr.itertuples(index=False):
        lig_mal, lig_mal_p = cached("Lmal", rec.ligand_genes, mal_idx)
        rec_tnk, rec_tnk_p = cached("Rtnk", rec.receptor_genes, tnk_idx)
        lig_tnk, lig_tnk_p = cached("Ltnk", rec.ligand_genes, tnk_idx)
        rec_mal, rec_mal_p = cached("Rmal", rec.receptor_genes, mal_idx)
        out_det = lig_mal_p >= EXPR_PROP and rec_tnk_p >= EXPR_PROP
        in_det = lig_tnk_p >= EXPR_PROP and rec_mal_p >= EXPR_PROP
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": "|".join(rec.ligand_genes),
                "receptor_genes": "|".join(rec.receptor_genes),
                "direction": "outgoing",
                "prob": hill_prob(lig_mal, rec_tnk) if out_det else 0.0,
                "detected": bool(out_det),
                "ligand_mean": lig_mal,
                "receptor_mean": rec_tnk,
                "ligand_prop": lig_mal_p,
                "receptor_prop": rec_tnk_p,
            }
        )
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": "|".join(rec.ligand_genes),
                "receptor_genes": "|".join(rec.receptor_genes),
                "direction": "incoming",
                "prob": hill_prob(lig_tnk, rec_mal) if in_det else 0.0,
                "detected": bool(in_det),
                "ligand_mean": lig_tnk,
                "receptor_mean": rec_mal,
                "ligand_prop": lig_tnk_p,
                "receptor_prop": rec_mal_p,
            }
        )
    return rows


def plot_q4q1(q1, q4, ylabel, title, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    bp = ax.boxplot(
        [q1, q4],
        tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"],
        patch_artist=True,
        widths=0.55,
    )
    for patch, color in zip(bp["boxes"], ["#6a8aaa", "#b2182b"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate((q1, q4), start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=16, zorder=3)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_scatter(x, y, hue, xlabel, ylabel, title, path: Path) -> None:
    colors = {"Q1": "#6a8aaa", "Q2": "#bdbdbd", "Q3": "#f4a582", "Q4": "#b2182b"}
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    for q, color in colors.items():
        m = hue == q
        ax.scatter(x[m], y[m], c=color, s=36, label=q, zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.legend(frameon=False, title="CLDN4 quartile")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ligand_table(tbl: pd.DataFrame, path: Path, title: str) -> None:
    if tbl.empty:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.axis("off")
        ax.text(0.5, 0.5, "No differential pairs at the stated gates", ha="center")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    show = tbl.head(18).copy()
    fig, ax = plt.subplots(figsize=(8.4, 0.42 * len(show) + 1.4))
    y = np.arange(len(show))
    colors = np.where(show["delta_median"] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, show["delta_median"], color=colors, alpha=0.85)
    ax.set_yticks(y)
    labels = [
        f"{r.direction[:3]} {r.ligand}–{r.receptor} ({r.ligand_class})"
        for r in show.itertuples()
    ]
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_xlabel("median P(Q4) − median P(Q1)")
    ax.set_title(title, fontsize=9)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def run_q4(patients: pd.DataFrame, out: Path) -> dict:
    patients = patients.copy()
    patients["frac_tnk"] = patients["n_tnk"] / patients["n_cells"]
    patients["q_pct"] = assign_quartiles(patients["mal_CLDN4_pct_pos"])
    patients["q_mean"] = assign_quartiles(patients["mal_CLDN4_mean"])
    rows = []
    for score, cldn_col, q_col in (
        ("pct_pos", "mal_CLDN4_pct_pos", "q_pct"),
        ("mean", "mal_CLDN4_mean", "q_mean"),
    ):
        rho, p_s, n = spearman(patients[cldn_col], patients["frac_tnk"])
        q = q4_vs_q1(patients[cldn_col], patients["frac_tnk"])
        row = {
            "cohort": "GSE205335",
            "malig_def": "author_malig",
            "immune_def": "same_patient_tnk_fraction",
            "score": score,
            "n_patients": n,
            "n_q1": q["n_q1"] if q else None,
            "n_q4": q["n_q4"] if q else None,
            "n_compared": q["n_compared"] if q else None,
            "spearman_rho": rho,
            "spearman_p": p_s,
            "q4q1_r_rb": q["r_rb"] if q else None,
            "q4q1_p": q["p"] if q else None,
            "median_tnk_q1": q["median_q1"] if q else None,
            "median_tnk_q4": q["median_q4"] if q else None,
            "delta_median_tnk": q["delta_median"] if q else None,
            "thin": q["thin"] if q else True,
            "unit": "patient",
        }
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(out / "results" / "q4q1_tnk.tsv", sep="\t", index=False)
    patients.to_csv(out / "results" / "patients_with_quartiles.tsv", sep="\t", index=False)

    q1 = patients.loc[patients["q_pct"] == "Q1", "frac_tnk"].to_numpy()
    q4 = patients.loc[patients["q_pct"] == "Q4", "frac_tnk"].to_numpy()
    prim = table.loc[table["score"] == "pct_pos"].iloc[0]
    plot_q4q1(
        q1,
        q4,
        "Same-patient T/NK fraction",
        f"GSE205335 CLDN4 %pos Q4 vs Q1 T/NK\nr={prim['q4q1_r_rb']:+.3f} p={fmt_p(prim['q4q1_p'])} n=6/6",
        out / "figures" / "q4q1_tnk_pct.png",
    )
    plot_scatter(
        patients["mal_CLDN4_pct_pos"],
        patients["frac_tnk"],
        patients["q_pct"],
        "Malignant CLDN4 % positive",
        "Same-patient T/NK fraction",
        f"GSE205335 n=22  ρ={prim['spearman_rho']:+.3f} p={fmt_p(prim['spearman_p'])}",
        out / "figures" / "scatter_cldn4_tnk.png",
    )
    return {"table": table, "patients": patients, "primary": prim.to_dict()}


def run_cellchat(args, patients: pd.DataFrame, out: Path) -> dict:
    import rdata

    identities = pd.read_csv(args.identities, sep="\t")
    if identities["barcode"].duplicated().any():
        raise ValueError("Cell-identity barcodes are not unique")
    metadata = parse_geo_soft(args.soft)
    metadata.to_csv(out / "results" / "gsm_sample_metadata.csv", index=False)

    # First pass: gene names only, via a cheap RDS read after we know wanted set.
    # Load all CellChat genes + extras; filter to those present after RDS open.
    dummy_genes = set(EXTRA)
    inter = pd.read_csv(args.db / "interaction_cellchatdb_v2_protein.csv")
    for _, rec in inter.iterrows():
        dummy_genes.update(parse_symbols(rec.get("ligand.symbol") or rec["ligand"]))
        dummy_genes.update(parse_symbols(rec.get("receptor.symbol") or rec["receptor"]))

    extracted, library_umi, barcodes, all_genes = load_selected_genes(args.matrix, dummy_genes)
    matrix_genes = set(all_genes)
    lr = load_lr(args.db, matrix_genes)
    print(f"LR pairs with all subunits in matrix: {len(lr)}", flush=True)
    if "CLDN4" not in extracted:
        raise SystemExit("CLDN4 not in UMI matrix")

    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise ValueError(f"Matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra")
    cells = indexed.loc[barcodes].reset_index()
    cells["total_umi"] = library_umi
    cells = cells.merge(
        metadata[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise ValueError("Some identity-table samples did not match GEO metadata")
    cells["compartment"] = np.select(
        [
            cells["lineage.sub"].eq("Malignant cells"),
            cells["lineage.total"].eq("T/NK cells"),
        ],
        ["Malignant", "T/NK"],
        default="Other",
    )
    lib = np.maximum(cells["total_umi"].to_numpy(), 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}

    keep_patients = set(patients["patient"])
    qmap = patients.set_index("patient")["q_pct"].astype(str)
    per_rows = []
    cell_rows = []
    for patient, sub in cells.groupby("patient", observed=True):
        if patient not in keep_patients:
            continue
        idx = sub.index.to_numpy()
        mal = idx[sub["compartment"].to_numpy() == "Malignant"]
        tnk = idx[sub["compartment"].to_numpy() == "T/NK"]
        cell_rows.append(
            {
                "patient": patient,
                "quartile_cldn4_pct": str(qmap.loc[patient]),
                "n_malignant": int(mal.size),
                "n_tnk": int(tnk.size),
                "n_other": int((sub["compartment"] == "Other").sum()),
                "mean_CLDN4_mal": float(log_cp["CLDN4"][mal].mean()) if mal.size else None,
                "pct_CLDN4_mal": float(100 * (extracted["CLDN4"][mal] > 0).mean()) if mal.size else None,
            }
        )
        scored = score_patient_pairs(lr, log_cp, pos, mal, tnk)
        for row in scored:
            row["patient"] = patient
            row["quartile_cldn4_pct"] = str(qmap.loc[patient])
            per_rows.append(row)
        print(f"  scored {patient} mal={mal.size} tnk={tnk.size} q={qmap.loc[patient]}", flush=True)

    per_patient_cells = pd.DataFrame(cell_rows)
    per_patient_cells.to_csv(out / "results" / "per_patient_cells.tsv", sep="\t", index=False)
    long = pd.DataFrame(per_rows)
    long.to_csv(out / "results" / "per_patient_lr.tsv.gz", sep="\t", index=False)

    contrast_rows = []
    for (name, direction), block in long.groupby(["interaction_name", "direction"], observed=True):
        q1 = block[block["quartile_cldn4_pct"] == "Q1"]
        q4 = block[block["quartile_cldn4_pct"] == "Q4"]
        d1 = q1[q1["detected"]]
        d4 = q4[q4["detected"]]
        if len(d1) < MIN_DETECT_ARM or len(d4) < MIN_DETECT_ARM:
            continue
        u, p = stats.mannwhitneyu(d4["prob"].to_numpy(), d1["prob"].to_numpy(), alternative="two-sided")
        n1, n4 = int(len(d1)), int(len(d4))
        r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
        med1 = float(d1["prob"].median())
        med4 = float(d4["prob"].median())
        rec0 = block.iloc[0]
        contrast_rows.append(
            {
                "interaction_name": name,
                "direction": direction,
                "pathway_name": rec0["pathway_name"],
                "annotation": rec0["annotation"],
                "ligand": rec0["ligand"],
                "receptor": rec0["receptor"],
                "ligand_genes": rec0["ligand_genes"],
                "receptor_genes": rec0["receptor_genes"],
                "ligand_class": ligand_class(rec0["ligand_genes"]),
                "n_q1_detected": n1,
                "n_q4_detected": n4,
                "n_compared": n1 + n4,
                "median_prob_q1": med1,
                "median_prob_q4": med4,
                "delta_median": med4 - med1,
                "mwu_u": float(u),
                "p": float(p),
                "r_rb": float(r_rb),
            }
        )
    contrast = pd.DataFrame(contrast_rows)
    if contrast.empty:
        ligand_tbl = contrast
    else:
        contrast = contrast[contrast["delta_median"] != 0].copy()
        contrast = contrast.assign(_absd=contrast["delta_median"].abs()).sort_values(["p", "_absd"], ascending=[True, False]).drop(columns="_absd")
        contrast.to_csv(out / "results" / "lr_q4q1_all.tsv", sep="\t", index=False)
        ligand_tbl = contrast[contrast["p"] < 0.05].copy()
        if ligand_tbl.empty:
            ligand_tbl = contrast.head(20).copy()
            ligand_tbl["note"] = "no pair p<0.05; top 20 by p shown"
        else:
            ligand_tbl["note"] = "patient-level MWU p<0.05; detected in ≥3 Q1 and ≥3 Q4"
        ligand_tbl = ligand_tbl.sort_values(["p", "delta_median"])
    ligand_tbl.to_csv(out / "results" / "ligand_table.tsv", sep="\t", index=False)
    plot_ligand_table(
        ligand_tbl,
        out / "figures" / "fig_extra_ligand_table.png",
        "GSE205335 CellChat-style Mal↔T/NK  Q4 vs Q1 (same-patient; n=6/6)",
    )
    return {
        "n_lr_in_matrix": int(len(lr)),
        "n_patients_scored": int(per_patient_cells.shape[0]),
        "n_pairs_tested_after_detect_gate": int(len(contrast)) if not contrast.empty else 0,
        "n_pairs_p_lt_05": int((contrast["p"] < 0.05).sum()) if not contrast.empty else 0,
        "ligand_table_rows": int(len(ligand_tbl)),
        "cell_counts": per_patient_cells.to_dict(orient="records"),
    }


def write_finding(q4: dict, cellchat: dict | None, out: Path) -> None:
    table = q4["table"]
    patients = q4["patients"]
    prim = table.loc[table["score"] == "pct_pos"].iloc[0]
    sec = table.loc[table["score"] == "mean"].iloc[0]
    q1_pts = patients.loc[patients["q_pct"] == "Q1"].sort_values("mal_CLDN4_pct_pos")
    q4_pts = patients.loc[patients["q_pct"] == "Q4"].sort_values("mal_CLDN4_pct_pos", ascending=False)

    lines = [
        "# FINDING — GSE205335 malignant CLDN4 vs same-patient T/NK Q4 + CellChat-style ligands",
        "",
        "ADDITIVE. **CLDN4 only.** Patient is the unit. Prior TACSTD2 A3 and the",
        "multi-cohort Q4 meta (PR #320) are given and are not re-ranked here.",
        "p-values are descriptive.",
        "",
        "## Honest n",
        "",
        "Locked extract: **22 patients** with ≥20 author-malignant and ≥20 T/NK",
        "cells (PR #279 / #320). Four GEO patients have (near-)zero captured",
        "malignant cells and are already out (3 PR + 1 PD). Q4 vs Q1 uses the",
        "quartile **tails only**: **n=6 vs 6** (n_compared=12), not 22.",
        "MPR/NMPR is unlabeled; RECIST is not used as MPR.",
        "",
        "## Malignant CLDN4 vs same-patient T/NK",
        "",
        "| score | n | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | median T/NK Q1 | median T/NK Q4 | Δ |",
        "|---|---:|---|---|---:|---:|---:|",
        (
            f"| CLDN4 %pos | {int(prim.n_patients)} | "
            f"{prim.spearman_rho:+.3f} ({fmt_p(prim.spearman_p)}) | "
            f"{prim.q4q1_r_rb:+.3f} ({fmt_p(prim.q4q1_p)}; "
            f"{int(prim.n_q1)}/{int(prim.n_q4)}) | "
            f"{prim.median_tnk_q1:.3f} | {prim.median_tnk_q4:.3f} | "
            f"{prim.delta_median_tnk:+.3f} |"
        ),
        (
            f"| CLDN4 mean log1p(CP10k) | {int(sec.n_patients)} | "
            f"{sec.spearman_rho:+.3f} ({fmt_p(sec.spearman_p)}) | "
            f"{sec.q4q1_r_rb:+.3f} ({fmt_p(sec.q4q1_p)}; "
            f"{int(sec.n_q1)}/{int(sec.n_q4)}) | "
            f"{sec.median_tnk_q1:.3f} | {sec.median_tnk_q4:.3f} | "
            f"{sec.delta_median_tnk:+.3f} |"
        ),
        "",
        "Primary cut is **%pos** (same row as PR #320 GSE205335 single:",
        "r=−0.778, p=0.026, 6/6). Mean is the same 22 patients, weaker.",
        "NSCLC-only (ADC+SQ) is not the verdict; it was already thinner in the",
        "given extract.",
        "",
        "### Quartile tails (CLDN4 %pos)",
        "",
        "| tail | patients (histology, RECIST) | n_mal | n_TNK | CLDN4 %pos | T/NK frac |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for rec in q4_pts.itertuples():
        lines.append(
            f"| Q4 | {rec.patient} ({rec.cancer_subtype}, {rec.recist}) | "
            f"{int(rec.n_malignant)} | {int(rec.n_tnk)} | "
            f"{rec.mal_CLDN4_pct_pos:.1f} | {rec.frac_tnk:.3f} |"
        )
    for rec in q1_pts.itertuples():
        lines.append(
            f"| Q1 | {rec.patient} ({rec.cancer_subtype}, {rec.recist}) | "
            f"{int(rec.n_malignant)} | {int(rec.n_tnk)} | "
            f"{rec.mal_CLDN4_pct_pos:.1f} | {rec.frac_tnk:.3f} |"
        )
    lines += [
        "",
        "Q4 mixes SCLC (P1025, P1115, P1016) with ADC. That histology mix is",
        "part of the honest n, not hidden.",
        "",
    ]

    if cellchat is None:
        lines += [
            "## CellChat-style ligands",
            "",
            "Not scored in this write-up (matrix step skipped or failed).",
            "Re-run without `--skip-cellchat` to fill `results/ligand_table.tsv`.",
            "",
        ]
    else:
        lig_path = out / "results" / "ligand_table.tsv"
        lig = pd.read_csv(lig_path, sep="\t") if lig_path.exists() else pd.DataFrame()
        lines += [
            "## CellChat-style ligands (same Q4 vs Q1 patients)",
            "",
            "Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs.",
            "Outgoing = author-malignant → **same-patient** T/NK. Incoming =",
            "T/NK → malignant. Test = Mann–Whitney on per-patient *P*",
            f"(n_Q1/n_Q4 detected ≥{MIN_DETECT_ARM}). CellChat R was not run.",
            "",
            (
                f"Pairs with all subunits in the matrix: {cellchat['n_lr_in_matrix']}. "
                f"After the detect gate: {cellchat['n_pairs_tested_after_detect_gate']} "
                f"direction×pair rows. p<0.05: {cellchat['n_pairs_p_lt_05']}."
            ),
            "",
        ]
        if lig.empty:
            lines.append("No pairs passed the detect gate.")
        else:
            lines += [
                "| direction | pair | class | n_Q1/n_Q4 | median P Q1 | median P Q4 | Δ | r | p |",
                "|---|---|---|---|---:|---:|---:|---:|---|",
            ]
            show = lig.head(25)
            for rec in show.itertuples():
                lines.append(
                    f"| {rec.direction} | {rec.ligand}–{rec.receptor} | {rec.ligand_class} | "
                    f"{int(rec.n_q1_detected)}/{int(rec.n_q4_detected)} | "
                    f"{rec.median_prob_q1:.3f} | {rec.median_prob_q4:.3f} | "
                    f"{rec.delta_median:+.3f} | {rec.r_rb:+.3f} | {fmt_p(rec.p)} |"
                )
            if "note" in lig.columns and len(lig):
                lines += ["", str(lig["note"].iloc[0]), ""]
        lines += [
            "Cell-pooled (stacked Q4 vs Q1) truncated means are **not** the test.",
            "n=6 vs 6 is thin; a single SCLC or LN sample can move a pair.",
            "",
        ]

    lines += [
        "## Files",
        "",
        "- `results/q4q1_tnk.tsv` — Spearman + Q4 vs Q1 T/NK",
        "- `results/patients_with_quartiles.tsv` — 22-patient table + quartile labels",
        "- `results/ligand_table.tsv` — CellChat-style differential pairs",
        "- `figures/q4q1_tnk_pct.png` — Q4 vs Q1 T/NK box",
        "- `figures/fig_extra_ligand_table.png` — extra ligand-table figure",
        "- `METHODS.md` — quartiles, author labels, Hill probability",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--patients", type=Path, default=ROOT / "data" / "GSE205335_patients.tsv")
    p.add_argument("--matrix", type=Path, default=Path("/tmp/gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz"))
    p.add_argument("--identities", type=Path, default=Path("/tmp/gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz"))
    p.add_argument("--soft", type=Path, default=Path("/tmp/gse205335/GSE205335_family.soft.gz"))
    p.add_argument("--db", type=Path, default=ROOT / "db")
    p.add_argument("--out", type=Path, default=ROOT)
    p.add_argument("--skip-cellchat", action="store_true")
    args = p.parse_args()
    (args.out / "results").mkdir(parents=True, exist_ok=True)
    (args.out / "figures").mkdir(parents=True, exist_ok=True)

    patients = pd.read_csv(args.patients, sep="\t")
    q4 = run_q4(patients, args.out)
    print(q4["table"].to_string(index=False), flush=True)

    cellchat = None
    if not args.skip_cellchat:
        cellchat = run_cellchat(args, q4["patients"], args.out)
        print(json.dumps({k: cellchat[k] for k in cellchat if k != "cell_counts"}, indent=2), flush=True)

    write_finding(q4, cellchat, args.out)
    summary = {
        "dataset": "GSE205335",
        "additive": True,
        "marker": "CLDN4",
        "unit": "patient",
        "n_patients": int(len(patients)),
        "q4q1": q4["table"].to_dict(orient="records"),
        "cellchat": None if cellchat is None else {k: v for k, v in cellchat.items() if k != "cell_counts"},
        "algorithm": "pd.qcut on average ranks; MWU rank-biserial; CellChat-like 10% trim mean, Hill Kh=0.5, expr_prop>=0.10, same-patient Mal↔T/NK",
        "not_run": "CellChat R; LIANA; EGA raw; MPR contrast; TACSTD2 split",
    }
    if args.matrix.exists():
        summary["matrix_sha256"] = sha256(args.matrix)
    (args.out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", args.out / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
