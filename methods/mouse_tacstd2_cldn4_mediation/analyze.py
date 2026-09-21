#!/usr/bin/env python3
"""Tacstd2 → Cldn4 → immune at the public-mouse unit, plus GSE137244.

Private 8 KL matrices are not inputs. Public series are not concatenated
with them. Numbers in the tables are computed here.

Primary mouse gate (the public GEMM rescore rule):
  QC: 200 <= nFeature <= 8000, nCount >= 500, mitochondrial fraction < 0.25
  epithelial: Epcam > 0 or Sftpc > 0 or (Krt8 > 0 and Ptprc == 0)
  T/NK: Cd3d > 0 or Cd3e > 0 or Nkg7 > 0 or Ncr1 > 0
  T fraction: n_T/NK / n_QC
  IFN: mean log1p(CP10k) of the 15-gene IFN-only set inside epithelial cells
  Cldn4 % and Tacstd2 %: fraction of epithelial QC cells with count > 0
  Means: mean log1p(CP10k) inside those epithelial cells

Primary mediation, locked before the fit:
  unit = biological mouse
  studies with at least 4 primary mice (WT ATTAC held out; n=2 studies are descriptive)
  exposure = within-study rank z of Tacstd2 %
  mediator = within-study rank z of Cldn4 %
  outcomes = T fraction, and epithelial IFN
  estimator = OLS product of coefficients on the concatenated rank-z values
  uncertainty = stratified bootstrap, and a within-study permutation null
  partial association = Pearson correlation of rank residuals on Cldn4

Raw-value z-scores, gene means instead of percents, equal-study weights,
genotype residuals, and the earlier marker-rule table are sensitivity.
They do not replace the primary rows.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
import tarfile
import bz2
from itertools import permutations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import scipy.sparse as sp
from scipy import stats

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
DATA = Path("/tmp/geo_dl")
SEED = 137244
N_BOOT = 5000
N_PERM = 5000

IFN_GENES = [
    "Stat1", "Stat2", "Irf1", "Irf7", "Irf9", "Isg15", "Ifit1", "Ifit2",
    "Ifit3", "Mx1", "Oasl2", "Rsad2", "Ifih1", "Ddx58", "Ifnb1",
]
TNK_GENES = ["Cd3d", "Cd3e", "Nkg7", "Ncr1"]

# GSE137244 library order in the GEO FPKM. Normal lung is held out of the test.
GSE137244_LIBS = [
    ("B6AL10-1-RNA", "KP"),
    ("B6AL10-2-RNA", "KP"),
    ("B6AL10-3-RNA", "KP"),
    ("B6AL10-4-RNA", "KP"),
    ("B6AL10-5-RNA", "KP"),
    ("KL155mix-control-2-RNA", "KL"),
    ("KL47-1-untreated-1-RNA", "KL"),
    ("KLC-RNA", "KL"),
    ("KLD-RNA", "KL"),
    ("KLE-RNA", "KL"),
    ("normal-lung-RNA", "normal"),
]


def _decode(value) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def _open_text(path: Path):
    raw = path.read_bytes()[:2]
    if raw == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return path.open("rt")


def read_features(path: Path) -> list[str]:
    genes = []
    with _open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    return genes


def read_mtx(matrix: Path, features: Path) -> tuple[sp.csr_matrix, list[str]]:
    genes = read_features(features)
    mat = scipy.io.mmread(matrix)
    X = mat.tocsr() if sp.issparse(mat) else sp.csr_matrix(mat)
    if X.shape[0] != len(genes) and X.shape[1] == len(genes):
        X = X.T.tocsr()
    if X.shape[0] != len(genes):
        raise ValueError(f"gene/matrix mismatch {len(genes)} vs {X.shape} in {matrix}")
    return X, genes


def read_h5(path: Path) -> tuple[sp.csr_matrix, list[str]]:
    import h5py

    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"] if "matrix" in handle else handle[list(handle.keys())[0]]
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = matrix["shape"][:]
        names = [_decode(x) for x in matrix["features"]["name"][:]]
        X = sp.csc_matrix((data, indices, indptr), shape=(int(shape[0]), int(shape[1])))
    if X.shape[0] != len(names) and X.shape[1] == len(names):
        X = X.T
    X = X.tocsr()
    if X.shape[0] != len(names):
        raise ValueError(f"h5 gene/matrix mismatch {len(names)} vs {X.shape} in {path}")
    return X, names


def score_matrix(X: sp.csr_matrix, genes: list[str]) -> dict:
    """Forest gate, plus Tacstd2 % and Tacstd2 mean on the same epithelial cells."""
    index: dict[str, int] = {}
    for i, gene in enumerate(genes):
        index.setdefault(gene, i)
    if "Cldn4" not in index or "Tacstd2" not in index:
        missing = [g for g in ("Cldn4", "Tacstd2") if g not in index]
        return {"error": "missing " + ",".join(missing)}

    n_cells_raw = X.shape[1]
    ncount = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    nfeat = np.diff(X.tocsc().indptr).astype(np.int32)
    mt_rows = [i for i, gene in enumerate(genes) if gene.lower().startswith("mt-")]
    if mt_rows and ncount.sum() > 0:
        mt_count = np.asarray(X[mt_rows].sum(axis=0)).ravel()
        pct_mt = np.divide(mt_count, ncount, out=np.zeros_like(ncount), where=ncount > 0)
    else:
        pct_mt = np.zeros(n_cells_raw)
    qc = (nfeat >= 200) & (nfeat <= 8000) & (ncount >= 500) & (pct_mt < 0.25)
    if int(qc.sum()) < 50:
        return {"error": f"too few QC cells ({int(qc.sum())})", "n_raw": n_cells_raw}

    Xq = X[:, qc]
    ncount_q = ncount[qc]
    n_qc = int(qc.sum())

    def gene_vec(name: str) -> np.ndarray | None:
        i = index.get(name)
        if i is None:
            return None
        return np.asarray(Xq[i].todense()).ravel()

    def positive(names: list[str]) -> np.ndarray:
        mask = np.zeros(n_qc, dtype=bool)
        for name in names:
            vec = gene_vec(name)
            if vec is not None:
                mask |= vec > 0
        return mask

    epi = positive(["Epcam"]) | positive(["Sftpc"]) | (positive(["Krt8"]) & ~positive(["Ptprc"]))
    tnk = positive(TNK_GENES)
    n_epi = int(epi.sum())
    ifn_present = [g for g in IFN_GENES if g in index]
    cldn4_pct = math.nan
    cldn4_mean = math.nan
    tac_pct = math.nan
    tac_mean = math.nan
    ifn = math.nan
    if n_epi >= 1 and ifn_present:
        epi_idx = np.flatnonzero(epi)
        scale = np.maximum(ncount_q[epi_idx], 1.0)

        def mean_log(name: str) -> tuple[float, float]:
            counts = np.asarray(Xq[index[name], epi_idx].todense()).ravel()
            logged = np.log1p(1e4 * counts / scale)
            return float((counts > 0).mean()), float(logged.mean())

        cldn4_pct, cldn4_mean = mean_log("Cldn4")
        tac_pct, tac_mean = mean_log("Tacstd2")
        acc = np.zeros(epi_idx.size, dtype=np.float64)
        for gene in ifn_present:
            counts = np.asarray(Xq[index[gene], epi_idx].todense()).ravel()
            acc += np.log1p(1e4 * counts / scale)
        ifn = float((acc / len(ifn_present)).mean())

    return {
        "n_raw": n_cells_raw,
        "n_qc": n_qc,
        "n_epi": n_epi,
        "n_tnk": int(tnk.sum()),
        "n_ptprc": int(positive(["Ptprc"]).sum()),
        "n_epcam": int(positive(["Epcam"]).sum()),
        "cldn4_pct": cldn4_pct,
        "cldn4_mean": cldn4_mean,
        "tacstd2_pct": tac_pct,
        "tacstd2_mean": tac_mean,
        "ifn": ifn,
        "n_ifn_genes": len(ifn_present),
        "error": "",
    }


def _member_to_bytes(tar: tarfile.TarFile, member: tarfile.TarInfo) -> bytes:
    extracted = tar.extractfile(member)
    if extracted is None:
        raise ValueError(f"cannot read {member.name}")
    return extracted.read()


def score_mtx_pair(matrix: Path, features: Path) -> dict:
    X, genes = read_mtx(matrix, features)
    scored = score_matrix(X, genes)
    del X
    return scored


def score_from_tar_member_tree(blob: bytes) -> dict:
    """Score a gzipped or nested 10x bundle (mtx + features)."""
    raw = blob
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    # GSE295824 members are gzip(bzip2(tar)).
    if raw[:3] == b"BZh":
        raw = bz2.decompress(raw)
    if raw[:5] == b"<?xml" or raw[:1] == b"<":
        raise ValueError("archive member is not a 10x bundle")

    def score_opened(bundle: tarfile.TarFile) -> dict:
        names = bundle.getnames()
        mtx_name = next(n for n in names if n.endswith("matrix.mtx.gz") or n.endswith("matrix.mtx"))
        feat_name = next(
            n for n in names if n.endswith("features.tsv.gz") or n.endswith("features.tsv") or n.endswith("genes.tsv.gz") or n.endswith("genes.tsv")
        )
        dest = DATA / "_tmp_bundle"
        dest.mkdir(parents=True, exist_ok=True)
        for name in (mtx_name, feat_name):
            member = bundle.getmember(name)
            target = dest / Path(name).name
            target.write_bytes(_member_to_bytes(bundle, member))
        try:
            return score_mtx_pair(dest / Path(mtx_name).name, dest / Path(feat_name).name)
        finally:
            for child in dest.glob("*"):
                child.unlink()

    if raw[:262].find(b"ustar") != -1 or raw[257:262] == b"ustar":
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:") as bundle:
            return score_opened(bundle)
    # already-gzipped mtx is not a tar; caller handles flat triplets
    raise ValueError("unrecognized bundle")


def score_h5_member(blob: bytes) -> dict:
    dest = DATA / "_tmp_one.h5"
    dest.write_bytes(blob)
    try:
        X, genes = read_h5(dest)
        return score_matrix(X, genes)
    finally:
        dest.unlink(missing_ok=True)


def add_row(rows: list[dict], accession: str, mouse: str, group: str, primary: bool, note: str, scored: dict) -> None:
    row = {
        "accession": accession,
        "mouse": mouse,
        "group": group,
        "primary_immune": primary,
        "note": note,
    }
    row.update(scored)
    if row.get("n_qc") and row.get("n_tnk") is not None and not row.get("error"):
        row["t_frac"] = row["n_tnk"] / row["n_qc"]
    else:
        row["t_frac"] = math.nan
    rows.append(row)
    print(
        f"  {accession} {mouse}: qc={row.get('n_qc')} epi={row.get('n_epi')} "
        f"tnk={row.get('n_tnk')} tac%={row.get('tacstd2_pct')} cld%={row.get('cldn4_pct')} err={row.get('error')}",
        flush=True,
    )


def score_flat_tar(tar_path: Path, parser) -> list[dict]:
    rows = []
    with tarfile.open(tar_path) as tar:
        members = {m.name: m for m in tar.getmembers() if m.isfile()}
        for spec in parser(members):
            print("scoring", spec["accession"], spec["mouse"], flush=True)
            if spec["kind"] == "h5":
                blob = _member_to_bytes(tar, members[spec["member"]])
                scored = score_h5_member(blob)
            elif spec["kind"] == "bundle":
                blob = _member_to_bytes(tar, members[spec["member"]])
                scored = score_from_tar_member_tree(blob)
            elif spec["kind"] == "triplet":
                dest = DATA / "_tmp_triplet"
                dest.mkdir(parents=True, exist_ok=True)
                paths = {}
                for key in ("matrix", "features"):
                    member = members[spec[key]]
                    target = dest / Path(member.name).name
                    target.write_bytes(_member_to_bytes(tar, member))
                    paths[key] = target
                try:
                    scored = score_mtx_pair(paths["matrix"], paths["features"])
                finally:
                    for child in dest.glob("*"):
                        child.unlink()
            else:
                raise ValueError(spec["kind"])
            add_row(rows, spec["accession"], spec["mouse"], spec["group"], spec["primary"], spec["note"], scored)
    return rows


def gse264739_specs(members: dict) -> list[dict]:
    specs = []
    for mouse, group in [("KP1", "KP"), ("KP2", "KP"), ("KP3", "KP"), ("KPP1", "KPP"), ("KPP2", "KPP"), ("KPP3", "KPP")]:
        gsm = {
            "KP1": "GSM8226902", "KP2": "GSM8226903", "KP3": "GSM8226904",
            "KPP1": "GSM8226905", "KPP2": "GSM8226906", "KPP3": "GSM8226907",
        }[mouse]
        specs.append({
            "kind": "triplet",
            "accession": "GSE264739",
            "mouse": mouse,
            "group": group,
            "primary": True,
            "note": "GEO replicate of autochthonous KP or KPP lung tumor",
            "matrix": f"{gsm}_{mouse}_matrix.mtx.gz",
            "features": f"{gsm}_{mouse}_features.tsv.gz",
        })
    missing = [s["matrix"] for s in specs if s["matrix"] not in members]
    if missing:
        raise SystemExit(f"GSE264739 missing {missing}; have {list(members)[:6]}")
    return specs


def gse295824_specs(members: dict) -> list[dict]:
    specs = []
    for name in sorted(members):
        if not name.endswith(".tar.gz"):
            continue
        stem = Path(name).name.replace(".tar.gz", "")
        group, mouse = stem.split("_", 1)[1].split("_")
        specs.append({
            "kind": "bundle",
            "accession": "GSE295824",
            "mouse": mouse,
            "group": group,
            "primary": True,
            "note": "Sox2 GEMM lung tumor; mouse id is in the archive name",
            "member": name,
        })
    return specs


def gse201247_specs(members: dict) -> list[dict]:
    mapping = [
        ("GSM6056008_1_ATTAC_filtered_feature_bc_matrix.h5", "ATTAC_1", "WT_ATTAC", False),
        ("GSM6056009_2_Kras_filtered_feature_bc_matrix.h5", "Kras_2", "Kras", True),
        ("GSM6056010_3_ATTAC_Kras_filtered_feature_bc_matrix.h5", "ATTAC_Kras_3", "ATTAC_Kras", True),
        ("GSM6056011_4_ATTAC_filtered_feature_bc_matrix.h5", "ATTAC_4", "WT_ATTAC", False),
        ("GSM6056012_5_Kras_filtered_feature_bc_matrix.h5", "Kras_5", "Kras", True),
        ("GSM6056013_6_ATTAC_Kras_filtered_feature_bc_matrix.h5", "ATTAC_Kras_6", "ATTAC_Kras", True),
    ]
    specs = []
    for member, mouse, group, primary in mapping:
        if member not in members:
            raise SystemExit(f"missing {member}")
        specs.append({
            "kind": "h5",
            "accession": "GSE201247",
            "mouse": mouse,
            "group": group,
            "primary": primary,
            "note": "whole-lung digest; WT ATTAC is scored and held out of the primary pool",
            "member": member,
        })
    return specs


def gse_bundle_pair_specs(members: dict, accession: str, pairs: list[tuple[str, str, str]]) -> list[dict]:
    specs = []
    for needle, mouse, group in pairs:
        hits = [name for name in members if needle in name]
        if len(hits) != 1:
            raise SystemExit(f"{accession} expected one member for {needle}, found {hits}")
        specs.append({
            "kind": "bundle",
            "accession": accession,
            "mouse": mouse,
            "group": group,
            "primary": True,
            "note": "mixed lung digest; n=2, descriptive only",
            "member": hits[0],
        })
    return specs


def score_gse266323() -> list[dict]:
    rows = []
    mapping = [
        ("d10_1", "Malat1_CRISPRa"),
        ("d10_2", "Malat1_CRISPRa"),
        ("dTom_1", "Tomato_control"),
        ("dTom_2", "Tomato_control"),
    ]
    for mouse, group in mapping:
        print("scoring GSE266323", mouse, flush=True)
        matrix = DATA / f"GSE266323_{mouse}_matrix.mtx.gz"
        features = DATA / f"GSE266323_{mouse}_features.tsv.gz"
        scored = score_mtx_pair(matrix, features)
        add_row(
            rows, "GSE266323", mouse, group, True,
            "whole-lung digest; d10 is Malat1 CRISPRa and dTom is the tomato control in the published mouse table",
            scored,
        )
    return rows


def score_gse179501() -> list[dict]:
    print("scoring GSE179501", flush=True)
    features = DATA / "GSE179501_features.tsv.gz"
    matrix = DATA / "GSE179501_matrix.mtx.gz"
    barcodes = []
    with _open_text(DATA / "GSE179501_barcodes.tsv.gz") as handle:
        barcodes = [line.split("_", 1)[0].strip() for line in handle if line.strip()]
    X, genes = read_mtx(matrix, features)
    if X.shape[1] != len(barcodes):
        raise SystemExit(f"GSE179501 cells {X.shape[1]} != barcodes {len(barcodes)}")
    groups = {
        "CM2260": "Restored",
        "CM2319": "Restored",
        "CM2324": "Non-Restored",
        "CM2328": "Non-Restored",
    }
    rows = []
    labels = np.array(barcodes)
    for mouse, group in groups.items():
        cols = np.flatnonzero(labels == mouse)
        if cols.size == 0:
            raise SystemExit(f"no cells for {mouse}")
        scored = score_matrix(X[:, cols], genes)
        add_row(
            rows, "GSE179501", mouse, group, True,
            "Lkb1-XTR total viable cells; Restored vs Non-Restored labels are the published mouse table",
            scored,
        )
    del X
    return rows


def score_all() -> list[dict]:
    rows: list[dict] = []
    rows += score_flat_tar(DATA / "GSE264739_RAW.tar", gse264739_specs)
    rows += score_flat_tar(DATA / "GSE295824_RAW.tar", gse295824_specs)
    rows += score_flat_tar(DATA / "GSE201247_RAW.tar", gse201247_specs)
    rows += score_gse266323()
    rows += score_flat_tar(
        DATA / "GSE165641_RAW.tar",
        lambda members: gse_bundle_pair_specs(members, "GSE165641", [("KL1", "KL1", "KL"), ("KL2", "KL2", "KL")]),
    )
    rows += score_flat_tar(
        DATA / "GSE180963_RAW.tar",
        lambda members: gse_bundle_pair_specs(members, "GSE180963", [("_K.tar", "K", "K"), ("_KL.tar", "KL", "KL")]),
    )
    rows += score_gse179501()
    return rows


def validate_gate(rows: list[dict]) -> list[dict]:
    targets = list(csv.DictReader((TABLES / "gate_validation_targets.tsv").open(), delimiter="\t"))
    by_key = {(r["accession"], r["mouse"]): r for r in rows}
    out = []
    worst = 0.0
    for target in targets:
        got = by_key.get((target["accession"], target["mouse"]))
        rec = {"accession": target["accession"], "mouse": target["mouse"], "ok": False}
        if got is None or got.get("error"):
            rec["error"] = "missing" if got is None else got.get("error")
            out.append(rec)
            continue
        pairs = {
            "n_qc": (int(got["n_qc"]), int(target["n_qc"])),
            "n_epi": (int(got["n_epi"]), int(target["n_epi"])),
            "cldn4_pct": (float(got["cldn4_pct"]), float(target["cldn4_pct"])),
            "t_frac": (float(got["t_frac"]), float(target["frac_tnk"])),
            "ifn": (float(got["ifn"]), float(target["ifn"])),
        }
        rec.update({f"{name}_got": a for name, (a, _) in pairs.items()})
        rec.update({f"{name}_target": b for name, (_, b) in pairs.items()})
        ok = (
            pairs["n_qc"][0] == pairs["n_qc"][1]
            and pairs["n_epi"][0] == pairs["n_epi"][1]
            and abs(pairs["cldn4_pct"][0] - pairs["cldn4_pct"][1]) < 1e-8
            and abs(pairs["t_frac"][0] - pairs["t_frac"][1]) < 1e-8
            and abs(pairs["ifn"][0] - pairs["ifn"][1]) < 1e-6
        )
        worst = max(worst, abs(pairs["ifn"][0] - pairs["ifn"][1]), abs(pairs["cldn4_pct"][0] - pairs["cldn4_pct"][1]))
        rec["ok"] = ok
        out.append(rec)
    n_ok = sum(1 for r in out if r["ok"])
    print(f"gate validation {n_ok}/{len(out)} exact, worst abs diff {worst:.3e}", flush=True)
    if n_ok != len(out):
        for rec in out:
            if not rec["ok"]:
                print("  FAIL", rec)
        raise SystemExit("epithelial gate did not reproduce the published rescore; Tacstd2 was not used")
    return out


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def as_float_frame(rows: list[dict]) -> list[dict]:
    numeric = {
        "n_raw", "n_qc", "n_epi", "n_tnk", "n_ptprc", "n_epcam", "cldn4_pct", "cldn4_mean",
        "tacstd2_pct", "tacstd2_mean", "ifn", "n_ifn_genes", "t_frac",
    }
    out = []
    for row in rows:
        rec = dict(row)
        for key in numeric:
            if key in rec and rec[key] != "" and rec[key] is not None:
                try:
                    rec[key] = float(rec[key])
                except (TypeError, ValueError):
                    rec[key] = math.nan
        out.append(rec)
    return out


def primary_mice(rows: list[dict]) -> list[dict]:
    kept = []
    for row in rows:
        if row.get("error"):
            continue
        if not row.get("primary_immune"):
            continue
        if row.get("n_epi", 0) < 20 or row.get("n_qc", 0) < 50:
            continue
        if not math.isfinite(row["tacstd2_pct"]) or not math.isfinite(row["cldn4_pct"]):
            continue
        if not math.isfinite(row["t_frac"]) or not math.isfinite(row["ifn"]):
            continue
        kept.append(row)
    return kept


def group_by(rows: list[dict], key: str) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        grouped.setdefault(str(row[key]), []).append(row)
    return grouped


def rank_z(values: np.ndarray) -> np.ndarray:
    ranks = stats.rankdata(values, method="average")
    sd = float(ranks.std(ddof=1))
    if sd == 0 or not math.isfinite(sd):
        return np.full(values.shape, np.nan)
    return (ranks - ranks.mean()) / sd


def raw_z(values: np.ndarray) -> np.ndarray:
    sd = float(values.std(ddof=1))
    if sd == 0 or not math.isfinite(sd):
        return np.full(values.shape, np.nan)
    return (values - values.mean()) / sd


def attach_z(rows: list[dict], columns: list[str], mode: str) -> list[dict] | None:
    """Return copies with <col>_<mode> z within accession. None if any study has no variance."""
    transform = rank_z if mode == "rank" else raw_z
    grouped = group_by(rows, "accession")
    pieces = []
    for accession, subset in grouped.items():
        if len(subset) < 2:
            return None
        arrays = {col: np.array([r[col] for r in subset], dtype=float) for col in columns}
        zipped = {col: transform(arrays[col]) for col in columns}
        if any(not np.isfinite(z).all() for z in zipped.values()):
            return None
        for i, row in enumerate(subset):
            rec = dict(row)
            for col in columns:
                rec[f"{col}__z"] = float(zipped[col][i])
            pieces.append(rec)
    return pieces


def fit_mediation(x: np.ndarray, m: np.ndarray, y: np.ndarray) -> dict[str, float]:
    if len(x) < 4 or not np.isfinite(x).all() or float(np.dot(x, x)) == 0 or float(np.dot(m, m)) == 0:
        return {k: math.nan for k in ("a", "b", "c", "cp", "ab", "n")}
    a = float(np.dot(x, m) / np.dot(x, x))
    c = float(np.dot(x, y) / np.dot(x, x))
    design = np.column_stack([x, m])
    cp, b = np.linalg.lstsq(design, y, rcond=None)[0]
    return {"a": a, "b": float(b), "c": c, "cp": float(cp), "ab": a * float(b), "n": float(len(x))}


def partial_r(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    def resid(a, b):
        design = np.column_stack([np.ones(len(a)), b])
        coef = np.linalg.lstsq(design, a, rcond=None)[0]
        return a - design @ coef
    xr, yr = resid(x, z), resid(y, z)
    if float(xr.std()) == 0 or float(yr.std()) == 0:
        return math.nan
    return float(stats.pearsonr(xr, yr).statistic)


def concat_z(rows: list[dict], xcol: str, mcol: str, ycol: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.array([r[xcol] for r in rows], dtype=float),
        np.array([r[mcol] for r in rows], dtype=float),
        np.array([r[ycol] for r in rows], dtype=float),
    )


def stratified_indices(groups: dict[str, np.ndarray], rng: np.random.Generator) -> dict[str, np.ndarray]:
    drawn = {}
    for key, idx in groups.items():
        drawn[key] = rng.choice(idx, size=len(idx), replace=True)
    return drawn


def pooled_mediation(rows: list[dict], xcol: str, mcol: str, ycol: str, mode: str, seed: int) -> dict:
    """Mouse-weighted mediation on within-study z-scores. Studies with n < 4 are excluded by the caller."""
    base = attach_z(rows, [xcol, mcol, ycol], mode)
    if base is None:
        return {"error": "a study had no variance"}
    x, m, y = concat_z(base, f"{xcol}__z", f"{mcol}__z", f"{ycol}__z")
    obs = fit_mediation(x, m, y)
    obs["partial_r"] = partial_r(x, y, m)
    obs["pearson_xy"] = float(stats.pearsonr(x, y).statistic)
    obs["pearson_xm"] = float(stats.pearsonr(x, m).statistic)

    by_study = group_by(base, "accession")
    # Equal-study mean of within-study indirect effects.
    study_abs = []
    for subset in by_study.values():
        sx, sm, sy = concat_z(subset, f"{xcol}__z", f"{mcol}__z", f"{ycol}__z")
        study_abs.append(fit_mediation(sx, sm, sy)["ab"])
    obs["ab_equal_study"] = float(np.nanmean(study_abs))
    obs["n_studies"] = float(len(by_study))

    rng = np.random.default_rng(seed)
    study_rows = {k: np.array(v, dtype=object) for k, v in by_study.items()}
    # Use positions within each study list.
    study_pos = {k: np.arange(len(v)) for k, v in by_study.items()}
    boot = []
    for _ in range(N_BOOT):
        drawn_rows = []
        ok = True
        for key, pos in study_pos.items():
            take = rng.choice(pos, size=len(pos), replace=True)
            drawn_rows.extend(study_rows[key][take].tolist())
        zrows = attach_z(drawn_rows, [xcol, mcol, ycol], mode)
        if zrows is None:
            continue
        bx, bm, by = concat_z(zrows, f"{xcol}__z", f"{mcol}__z", f"{ycol}__z")
        boot.append(fit_mediation(bx, bm, by)["ab"])
    boot_arr = np.array([v for v in boot if math.isfinite(v)], dtype=float)
    if boot_arr.size:
        obs["ab_ci_low"] = float(np.quantile(boot_arr, 0.025))
        obs["ab_ci_high"] = float(np.quantile(boot_arr, 0.975))
        obs["n_boot"] = float(boot_arr.size)
    else:
        obs["ab_ci_low"] = math.nan
        obs["ab_ci_high"] = math.nan
        obs["n_boot"] = 0

    perm = []
    for _ in range(N_PERM):
        perm_rows = []
        for key, subset in by_study.items():
            order = rng.permutation(len(subset))
            ys = [subset[i][ycol] for i in order]
            for row, yval in zip(subset, ys):
                rec = dict(row)
                rec[ycol] = yval
                perm_rows.append(rec)
        zrows = attach_z(perm_rows, [xcol, mcol, ycol], mode)
        if zrows is None:
            continue
        px, pm, py = concat_z(zrows, f"{xcol}__z", f"{mcol}__z", f"{ycol}__z")
        perm.append(fit_mediation(px, pm, py)["ab"])
    perm_arr = np.array([v for v in perm if math.isfinite(v)], dtype=float)
    if perm_arr.size and math.isfinite(obs["ab"]):
        obs["ab_perm_p"] = float((1 + np.sum(np.abs(perm_arr) >= abs(obs["ab"]) - 1e-15)) / (1 + perm_arr.size))
    else:
        obs["ab_perm_p"] = math.nan
    obs["n_perm"] = float(perm_arr.size)

    # Permutation p for the total association and the partial association.
    def perm_p(stat_name: str, observed: float) -> float:
        null = []
        for _ in range(N_PERM):
            perm_rows = []
            for subset in by_study.values():
                order = rng.permutation(len(subset))
                ys = [subset[i][ycol] for i in order]
                for row, yval in zip(subset, ys):
                    rec = dict(row)
                    rec[ycol] = yval
                    perm_rows.append(rec)
            zrows = attach_z(perm_rows, [xcol, mcol, ycol], mode)
            if zrows is None:
                continue
            px, pm, py = concat_z(zrows, f"{xcol}__z", f"{mcol}__z", f"{ycol}__z")
            if stat_name == "pearson_xy":
                null.append(float(stats.pearsonr(px, py).statistic))
            else:
                null.append(partial_r(px, py, pm))
        arr = np.array([v for v in null if math.isfinite(v)], dtype=float)
        if not arr.size or not math.isfinite(observed):
            return math.nan
        return float((1 + np.sum(np.abs(arr) >= abs(observed) - 1e-15)) / (1 + arr.size))

    obs["pearson_xy_perm_p"] = perm_p("pearson_xy", obs["pearson_xy"])
    obs["partial_r_perm_p"] = perm_p("partial_r", obs["partial_r"])
    obs["mode"] = mode
    obs["x"] = xcol
    obs["m"] = mcol
    obs["y"] = ycol
    return obs


def spearman_p(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float, str]:
    rho = float(stats.spearmanr(x, y).statistic)
    n = len(x)
    if n <= 7:
        count = 0
        total = 0
        for perm in permutations(range(n)):
            total += 1
            r = float(stats.spearmanr(x, y[list(perm)]).statistic)
            if abs(r) + 1e-12 >= abs(rho):
                count += 1
        return rho, count / total, f"exact/{total}"
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        null[i] = float(stats.spearmanr(x, rng.permutation(y)).statistic)
    p = float((1 + np.sum(np.abs(null) >= abs(rho) - 1e-15)) / (1 + N_PERM))
    return rho, p, f"mc/{N_PERM}"


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    return partial_r(stats.rankdata(x), stats.rankdata(y), stats.rankdata(z))


def partial_spearman_p(x: np.ndarray, y: np.ndarray, z: np.ndarray, rng: np.random.Generator) -> tuple[float, float, str]:
    rho = partial_spearman(x, y, z)
    n = len(x)
    if not math.isfinite(rho):
        return rho, math.nan, "undefined"
    if n <= 7:
        count = 0
        total = 0
        for perm in permutations(range(n)):
            total += 1
            r = partial_spearman(x, y[list(perm)], z)
            if math.isfinite(r) and abs(r) + 1e-12 >= abs(rho):
                count += 1
        return rho, count / total, f"exact/{total}"
    count = 0
    total = 0
    for _ in range(N_PERM):
        r = partial_spearman(x, rng.permutation(y), z)
        if not math.isfinite(r):
            continue
        total += 1
        if abs(r) + 1e-12 >= abs(rho):
            count += 1
    p = float((1 + count) / (1 + total)) if total else math.nan
    return rho, p, f"mc/{total}"


def within_study_table(rows: list[dict], xcol: str, mcol: str, ycol: str, outcome: str) -> list[dict]:
    rng = np.random.default_rng(SEED)
    out = []
    for accession, subset in sorted(group_by(rows, "accession").items()):
        rec = {
            "accession": accession,
            "outcome": outcome,
            "x": xcol,
            "mediator": mcol,
            "n": len(subset),
            "groups": ",".join(sorted({r["group"] for r in subset})),
        }
        if len(subset) < 4:
            rec["spearman"] = ""
            rec["note"] = "n<4; not a correlation test"
            out.append(rec)
            continue
        x = np.array([r[xcol] for r in subset], dtype=float)
        m = np.array([r[mcol] for r in subset], dtype=float)
        y = np.array([r[ycol] for r in subset], dtype=float)
        rho, p, kind = spearman_p(x, y, rng)
        pr, pp, pkind = partial_spearman_p(x, y, m, rng)
        rho_xm, p_xm, _ = spearman_p(x, m, rng)
        fit = fit_mediation(rank_z(x), rank_z(m), rank_z(y))
        rec.update({
            "spearman_xy": rho,
            "spearman_xy_p": p,
            "spearman_xy_p_kind": kind,
            "spearman_xm": rho_xm,
            "spearman_xm_p": p_xm,
            "partial_spearman_xy_m": pr,
            "partial_spearman_p": pp,
            "partial_p_kind": pkind,
            "a": fit["a"],
            "b": fit["b"],
            "c": fit["c"],
            "cp": fit["cp"],
            "ab": fit["ab"],
            "note": "",
        })
        out.append(rec)
    return out


def genotype_residual_rows(rows: list[dict], columns: list[str]) -> list[dict] | None:
    """Residualize columns on group indicators inside each study. Drop a study with too few residual df."""
    kept = []
    for accession, subset in group_by(rows, "accession").items():
        groups = sorted({r["group"] for r in subset})
        df = len(subset) - len(groups)
        if len(groups) < 2 or df < 3:
            continue
        design_cols = groups[1:]
        design = np.column_stack([
            np.ones(len(subset)),
            *[[1.0 if r["group"] == g else 0.0 for r in subset] for g in design_cols],
        ])
        updated = [dict(r) for r in subset]
        for col in columns:
            y = np.array([r[col] for r in subset], dtype=float)
            coef = np.linalg.lstsq(design, y, rcond=None)[0]
            resid = y - design @ coef
            for row, value in zip(updated, resid):
                row[col] = float(value)
        kept.extend(updated)
    if len({r["accession"] for r in kept}) < 2:
        return None
    return kept


def fmt(value, digits=3) -> str:
    if value is None or value == "":
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return "NA"
    if abs(number) != 0 and abs(number) < 0.001:
        return f"{number:.2e}"
    return f"{number:.{digits}f}"


def load_fpkm() -> dict[str, dict[str, float]]:
    path = DATA / "GSE137244_counts.fpkm.csv.gz"
    with gzip.open(path, "rt") as handle:
        header = next(handle).strip().strip('"').split(",")
        samples = [h.strip('"') for h in header[1:]]
        expected = [name for name, _ in GSE137244_LIBS]
        if samples != expected:
            raise SystemExit(f"GSE137244 sample order {samples}")
        table = {name: {} for name in samples}
        for line in handle:
            parts = [x.strip('"') for x in line.strip().split(",")]
            symbol = parts[0]
            for sample, raw in zip(samples, parts[1:]):
                table[sample][symbol] = math.log2(float(raw) + 1.0)
    return table


def gse137244() -> tuple[list[dict], list[dict]]:
    expr = load_fpkm()
    tumor = [(name, geno) for name, geno in GSE137244_LIBS if geno in {"KP", "KL"}]
    libraries = []
    for name, geno in GSE137244_LIBS:
        ifn_vals = [expr[name][g] for g in IFN_GENES if g in expr[name]]
        if len(ifn_vals) != len(IFN_GENES):
            raise SystemExit(f"{name} missing IFN genes")
        libraries.append({
            "library": name,
            "genotype": geno,
            "in_test": geno in {"KP", "KL"},
            "Cldn4": expr[name]["Cldn4"],
            "Tacstd2": expr[name]["Tacstd2"],
            "Epcam": expr[name]["Epcam"],
            "IFN": float(np.mean(ifn_vals)),
        })
    kp = [r for r in libraries if r["genotype"] == "KP"]
    kl = [r for r in libraries if r["genotype"] == "KL"]
    for gene, locked in (("Cldn4", 5.57), ("Tacstd2", 3.24)):
        delta = float(np.mean([r[gene] for r in kl]) - np.mean([r[gene] for r in kp]))
        if abs(delta - locked) > 0.01:
            raise SystemExit(f"GSE137244 {gene} delta {delta} != locked {locked}")
    rng = np.random.default_rng(SEED)
    tests = []

    def pack(label: str, subset: list[dict]) -> dict:
        x = np.array([r["Tacstd2"] for r in subset], dtype=float)
        m = np.array([r["Cldn4"] for r in subset], dtype=float)
        y = np.array([r["IFN"] for r in subset], dtype=float)
        rho, p, kind = spearman_p(x, y, rng)
        pr, pp, pkind = partial_spearman_p(x, y, m, rng)
        rx, px, xkind = spearman_p(x, m, rng)
        fit = fit_mediation(*(rank_z(v) for v in (x, m, y)))
        return {
            "subset": label,
            "n": len(subset),
            "spearman_tac_cldn4": rx,
            "spearman_tac_cldn4_p": px,
            "spearman_tac_cldn4_p_kind": xkind,
            "spearman_tac_ifn": rho,
            "spearman_tac_ifn_p": p,
            "spearman_tac_ifn_p_kind": kind,
            "partial_tac_ifn_cldn4": pr,
            "partial_p": pp,
            "partial_p_kind": pkind,
            "a": fit["a"],
            "b": fit["b"],
            "c": fit["c"],
            "cp": fit["cp"],
            "ab": fit["ab"],
            "t_frac": "not defined; these are cell-line libraries",
        }

    tests.append(pack("10 tumor libraries", [r for r in libraries if r["in_test"]]))
    tests.append(pack("KP only", kp))
    tests.append(pack("KL only", kl))
    # Genotype indicator residual on the 10 libraries.
    ten = [r for r in libraries if r["in_test"]]
    design = np.column_stack([
        np.ones(len(ten)),
        [1.0 if r["genotype"] == "KL" else 0.0 for r in ten],
    ])
    residual = []
    for row in ten:
        residual.append(dict(row))
    for col in ("Tacstd2", "Cldn4", "IFN"):
        y = np.array([r[col] for r in ten], dtype=float)
        coef = np.linalg.lstsq(design, y, rcond=None)[0]
        for rec, value in zip(residual, y - design @ coef):
            rec[col] = float(value)
    tests.append(pack("10 libraries, genotype residual", residual))
    gated = [r for r in ten if r["Epcam"] >= 4.0]
    tests.append(pack("Epcam log2>=4", gated))
    return libraries, tests


def pipeline_b() -> list[dict]:
    """Within-study tests on the earlier public marker-rule mouse table.

    That table is not pooled with the count-gate mice. GSE165641, GSE180963,
    and GSE179501 are omitted here because this script rescores them.
    """
    path = TABLES / "pipeline_b_marker_rule_mouse.tsv"
    raw = list(csv.DictReader(path.open(), delimiter="\t"))
    rows = []
    for rec in raw:
        if rec["dataset"] not in {"GSE136246", "GSE154977", "GSE154989", "GSE179502"}:
            continue
        if rec["epi_score_ok"] != "True":
            continue
        try:
            n_epi = float(rec["n_epi"])
            tac = float(rec["epi_frac_pos_Tacstd2"])
            cld = float(rec["epi_frac_pos_Cldn4"])
            ifn = float(rec["epi_mean_IFN"]) if rec["epi_mean_IFN"] else math.nan
            frac_t = float(rec["frac_T"]) if rec["frac_T"] else math.nan
        except ValueError:
            continue
        if n_epi < 20 or not math.isfinite(tac) or not math.isfinite(cld):
            continue
        genotype = rec["genotype_group"]
        if rec["dataset"] == "GSE154989" and genotype not in {"K", "KP"}:
            continue
        rows.append({
            "accession": rec["dataset"],
            "mouse": rec["mouse_id"],
            "group": genotype,
            "design": rec["design"],
            "scale_mode": rec["scale_mode"],
            "tacstd2_pct": tac,
            "cldn4_pct": cld,
            "t_frac": frac_t,
            "ifn": ifn,
            "primary_immune": rec["design"] == "mixed_lung",
        })
    rng_rows = []
    # T fraction only where the design is an unsorted digest.
    t_rows = [r for r in rows if r["primary_immune"] and math.isfinite(r["t_frac"]) and r["t_frac"] > 0]
    ifn_rows = [r for r in rows if math.isfinite(r["ifn"])]
    out = []
    for subset, ycol, outcome in (
        (t_rows, "t_frac", "T_frac"),
        (ifn_rows, "ifn", "IFN"),
    ):
        table = within_study_table(subset, "tacstd2_pct", "cldn4_pct", ycol, outcome)
        for rec in table:
            rec["pipeline"] = "marker_rule_not_pooled"
        out.extend(table)
    # GSE154989 genotype residual, IFN only, within that one study.
    plate = [r for r in ifn_rows if r["accession"] == "GSE154989"]
    if len(plate) >= 8:
        groups = sorted({r["group"] for r in plate})
        design = np.column_stack([
            np.ones(len(plate)),
            *[[1.0 if r["group"] == g else 0.0 for r in plate] for g in groups[1:]],
        ])
        residual = [dict(r) for r in plate]
        for col in ("tacstd2_pct", "cldn4_pct", "ifn"):
            y = np.array([r[col] for r in plate], dtype=float)
            coef = np.linalg.lstsq(design, y, rcond=None)[0]
            for rec, value in zip(residual, y - design @ coef):
                rec[col] = float(value)
        x = np.array([r["tacstd2_pct"] for r in residual])
        m = np.array([r["cldn4_pct"] for r in residual])
        y = np.array([r["ifn"] for r in residual])
        rng = np.random.default_rng(SEED)
        rho, p, kind = spearman_p(x, y, rng)
        pr, pp, pkind = partial_spearman_p(x, y, m, rng)
        out.append({
            "accession": "GSE154989",
            "outcome": "IFN_genotype_residual",
            "x": "tacstd2_pct",
            "mediator": "cldn4_pct",
            "n": len(residual),
            "groups": ",".join(groups),
            "spearman_xy": rho,
            "spearman_xy_p": p,
            "spearman_xy_p_kind": kind,
            "partial_spearman_xy_m": pr,
            "partial_spearman_p": pp,
            "partial_p_kind": pkind,
            "pipeline": "marker_rule_not_pooled",
            "note": "K vs KP indicator removed before the correlation; plate tumor cells, not a T fraction",
        })
    return out


def forest_plot(within: list[dict], outcome: str, path: Path) -> None:
    rows = [r for r in within if r["outcome"] == outcome and r.get("spearman_xy") not in ("", None)]
    if not rows:
        return
    labels = [f"{r['accession']} (n={r['n']})" for r in rows]
    ypos = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(7.2, 0.55 * len(rows) + 1.6))
    ax.axvline(0, color="#888888", lw=0.8)
    ax.scatter(
        [float(r["spearman_xy"]) for r in rows], ypos + 0.12,
        color="#1f4e79", label="Tacstd2 vs outcome", zorder=3,
    )
    ax.scatter(
        [float(r["partial_spearman_xy_m"]) for r in rows], ypos - 0.12,
        color="#b85c38", label="partial | Cldn4", zorder=3,
    )
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Spearman ρ")
    ax.set_xlim(-1.05, 1.05)
    ax.set_title(f"Within-study Tacstd2 % vs {outcome}")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def path_plot(pooled: list[dict], path: Path) -> None:
    primary = [r for r in pooled if r.get("role") == "primary"]
    if not primary:
        return
    names = {"t_frac": "T fraction", "ifn": "epithelial IFN"}
    labels = [names.get(r["y"], r["y"]) for r in primary]
    fig, ax = plt.subplots(figsize=(7.4, 3.6))
    x = np.arange(len(primary))
    width = 0.18
    series = [("a  Tacstd2→Cldn4", "a"), ("b  Cldn4→outcome | Tacstd2", "b"), ("c  total", "c"), ("c' direct", "cp"), ("a×b indirect", "ab")]
    for i, (label, key) in enumerate(series):
        ax.bar(x + (i - 2) * width, [float(r[key]) for r in primary], width=width, label=label)
    ax.axhline(0, color="#444444", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Standardized coefficient")
    ax.set_title("Primary mouse-unit mediation (rank z)")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def cellline_plot(libraries: list[dict], path: Path) -> None:
    tumor = [r for r in libraries if r["in_test"]]
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))
    colors = {"KP": "#4c78a8", "KL": "#f58518"}
    for ax, ykey, ylabel in (
        (axes[0], "Cldn4", "Cldn4 log2(FPKM+1)"),
        (axes[1], "IFN", "IFN-only mean log2(FPKM+1)"),
    ):
        for geno in ("KP", "KL"):
            subset = [r for r in tumor if r["genotype"] == geno]
            ax.scatter(
                [r["Tacstd2"] for r in subset],
                [r[ykey] for r in subset],
                c=colors[geno], label=geno, s=36,
            )
        ax.set_xlabel("Tacstd2 log2(FPKM+1)")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def path_a_test(rows: list[dict]) -> dict:
    """Within-study rank-z association of Tacstd2 % with Cldn4 %. Same seed as the mediation."""
    by = group_by(rows, "accession")
    xs = []
    ms = []
    for subset in by.values():
        x = np.array([r["tacstd2_pct"] for r in subset], dtype=float)
        m = np.array([r["cldn4_pct"] for r in subset], dtype=float)
        xs.append(rank_z(x))
        ms.append(rank_z(m))
    x = np.concatenate(xs)
    m = np.concatenate(ms)
    a = float(np.dot(x, m) / np.dot(x, x))
    rng = np.random.default_rng(SEED)
    null = []
    subsets = list(by.values())
    for _ in range(N_PERM):
        parts = []
        for subset in subsets:
            mvals = np.array([r["cldn4_pct"] for r in subset], dtype=float)
            parts.append(rank_z(rng.permutation(mvals)))
        mp = np.concatenate(parts)
        null.append(float(np.dot(x, mp) / np.dot(x, x)))
    arr = np.array(null)
    p = float((1 + np.sum(np.abs(arr) >= abs(a) - 1e-15)) / (1 + arr.size))
    return {
        "x": "tacstd2_pct",
        "m": "cldn4_pct",
        "mode": "rank",
        "n": float(len(rows)),
        "n_studies": float(len(by)),
        "a": a,
        "perm_p": p,
        "n_perm": float(arr.size),
    }


def jsonable(value):
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [jsonable(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (np.floating,)):
        number = float(value)
        return None if not math.isfinite(number) else number
    if isinstance(value, (np.integer,)):
        return int(value)
    return value


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rows = score_all()
    rows = as_float_frame(rows)
    for row in rows:
        row["primary_immune"] = bool(row["primary_immune"])
    validation = validate_gate(rows)
    write_tsv(TABLES / "mouse_scores.tsv", rows)
    write_tsv(TABLES / "gate_validation.tsv", validation)

    eligible = primary_mice(rows)
    counts = {k: len(v) for k, v in sorted(group_by(eligible, "accession").items())}
    pool = [r for r in eligible if counts[r["accession"]] >= 4]
    descriptive = [r for r in eligible if counts[r["accession"]] < 4]
    print("primary pool", {k: counts[k] for k in counts if counts[k] >= 4}, "descriptive", len(descriptive), flush=True)

    within = []
    within += within_study_table(pool, "tacstd2_pct", "cldn4_pct", "t_frac", "T_frac")
    within += within_study_table(pool, "tacstd2_pct", "cldn4_pct", "ifn", "IFN")
    # Means are a within-study sensitivity on the same mice.
    within += within_study_table(pool, "tacstd2_mean", "cldn4_mean", "t_frac", "T_frac_mean")
    within += within_study_table(pool, "tacstd2_mean", "cldn4_mean", "ifn", "IFN_mean")
    write_tsv(TABLES / "within_study.tsv", within)

    pooled = []
    specs = [
        ("primary", "rank", "tacstd2_pct", "cldn4_pct", "t_frac"),
        ("primary", "rank", "tacstd2_pct", "cldn4_pct", "ifn"),
        ("sensitivity_raw_z", "raw", "tacstd2_pct", "cldn4_pct", "t_frac"),
        ("sensitivity_raw_z", "raw", "tacstd2_pct", "cldn4_pct", "ifn"),
        ("sensitivity_mean", "rank", "tacstd2_mean", "cldn4_mean", "t_frac"),
        ("sensitivity_mean", "rank", "tacstd2_mean", "cldn4_mean", "ifn"),
    ]
    for role, mode, xcol, mcol, ycol in specs:
        print("pooled", role, mode, ycol, flush=True)
        fit = pooled_mediation(pool, xcol, mcol, ycol, mode, SEED)
        fit["role"] = role
        pooled.append(fit)

    geno_rows = genotype_residual_rows(pool, ["tacstd2_pct", "cldn4_pct", "t_frac", "ifn"])
    if geno_rows is not None:
        for ycol in ("t_frac", "ifn"):
            print("genotype residual", ycol, "n", len(geno_rows), flush=True)
            fit = pooled_mediation(geno_rows, "tacstd2_pct", "cldn4_pct", ycol, "rank", SEED + 1)
            fit["role"] = "sensitivity_genotype_residual"
            pooled.append(fit)
    write_tsv(TABLES / "pooled_mediation.tsv", pooled)

    # Descriptive n=2 signs, not correlations.
    signs = []
    for accession, subset in sorted(group_by(descriptive, "accession").items()):
        ordered = sorted(subset, key=lambda r: (r["group"], r["mouse"]))
        for row in ordered:
            signs.append({
                "accession": accession,
                "mouse": row["mouse"],
                "group": row["group"],
                "n_epi": row["n_epi"],
                "tacstd2_pct": row["tacstd2_pct"],
                "cldn4_pct": row["cldn4_pct"],
                "t_frac": row["t_frac"],
                "ifn": row["ifn"],
                "note": "n<4; listed, not tested",
            })
    write_tsv(TABLES / "descriptive_n_lt4.tsv", signs)

    libraries, cell_tests = gse137244()
    write_tsv(TABLES / "gse137244_libraries.tsv", libraries)
    write_tsv(TABLES / "gse137244_tests.tsv", cell_tests)

    b_rows = pipeline_b()
    write_tsv(TABLES / "pipeline_b_within.tsv", b_rows)

    forest_plot(within, "T_frac", FIGS / "within_T_frac.png")
    forest_plot(within, "IFN", FIGS / "within_IFN.png")
    path_plot(pooled, FIGS / "mediation_paths.png")
    cellline_plot(libraries, FIGS / "gse137244_tacstd2.png")

    path_a = path_a_test(pool)
    write_tsv(TABLES / "path_a.tsv", [path_a])

    summary = {
        "n_scored": len(rows),
        "n_eligible": len(eligible),
        "pool_counts": {k: v for k, v in counts.items() if v >= 4},
        "primary": [r for r in pooled if r.get("role") == "primary"],
        "path_a": path_a,
        "private_8kl_merged": False,
    }
    (TABLES / "summary.json").write_text(json.dumps(jsonable(summary), indent=2) + "\n")
    print(json.dumps(jsonable(summary["primary"]), indent=2))


if __name__ == "__main__":
    main()
