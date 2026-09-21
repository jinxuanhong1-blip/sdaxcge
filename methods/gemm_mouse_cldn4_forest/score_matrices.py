#!/usr/bin/env python3
"""Score public GEMM lung scRNA libraries at the biological-mouse unit.

Primary definitions are locked here and are not chosen by p-value:

- Cldn4 %pos = fraction of epithelial QC cells with Cldn4 count > 0
- epithelial = Epcam > 0 or Sftpc > 0 or (Krt8 > 0 and Ptprc == 0)
- T/NK = Cd3d > 0 or Cd3e > 0 or Nkg7 > 0 or Ncr1 > 0
- IFN = mean log1p(CP10k) of a fixed IFN-only gene set inside epithelial cells
- QC: 200 <= nFeature <= 8000, nCount >= 500, percent.mt < 25
- Mouse unit: libraries that GEO identifies as one biological replicate / mouse.
  Pools, pre/post pairs without a mouse id, and sort fractions of unknown mice
  are scored at library level and flagged out of the mouse forest.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp

EXTRACT = Path("/tmp/geo_dl/extract")
OUT = Path("/workspace/methods/gemm_mouse_cldn4_forest/tables")

IFN_GENES = [
    "Stat1", "Stat2", "Irf1", "Irf7", "Irf9", "Isg15", "Ifit1", "Ifit2",
    "Ifit3", "Mx1", "Oasl2", "Rsad2", "Ifih1", "Ddx58", "Ifnb1",
]
TNK_GENES = ["Cd3d", "Cd3e", "Nkg7", "Ncr1"]


def _decode(x):
    if isinstance(x, bytes):
        return x.decode()
    return str(x)


def read_features(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    genes = []
    with opener(path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    return genes


def read_mtx(matrix: Path, features: Path) -> tuple[sp.csr_matrix, list[str]]:
    genes = read_features(features)
    mat = scipy.io.mmread(matrix)
    if sp.issparse(mat):
        X = mat.tocsr()
    else:
        X = sp.csr_matrix(mat)
    # 10x Market Matrix is genes x cells.
    if X.shape[0] != len(genes) and X.shape[1] == len(genes):
        X = X.T.tocsr()
    if X.shape[0] != len(genes):
        raise ValueError(f"gene/matrix mismatch {len(genes)} vs {X.shape} in {matrix}")
    return X, genes


def read_h5(path: Path) -> tuple[sp.csr_matrix, list[str]]:
    import h5py

    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
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


def _row_index(genes: list[str]) -> dict[str, int]:
    index = {}
    for i, gene in enumerate(genes):
        index.setdefault(gene, i)
    return index


def score_matrix(X: sp.csr_matrix, genes: list[str]) -> dict:
    index = _row_index(genes)
    if "Cldn4" not in index:
        return {"error": "Cldn4 absent"}
    n_cells_raw = X.shape[1]
    ncount = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    nfeat = np.diff(X.tocsc().indptr).astype(np.int32)
    mt_rows = [i for i, g in enumerate(genes) if g.lower().startswith("mt-")]
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

    epcam = positive(["Epcam"])
    sftpc = positive(["Sftpc"])
    krt8 = positive(["Krt8"])
    ptprc = positive(["Ptprc"])
    epi = epcam | sftpc | (krt8 & ~ptprc)
    tnk = positive(TNK_GENES)
    n_epi = int(epi.sum())
    cldn4 = gene_vec("Cldn4")
    assert cldn4 is not None

    ifn_present = [g for g in IFN_GENES if g in index]
    sum_ifn = 0.0
    sum_cldn4_log = 0.0
    n_pos = 0
    if n_epi >= 1 and ifn_present:
        epi_idx = np.flatnonzero(epi)
        scale = ncount_q[epi_idx]
        scale = np.maximum(scale, 1.0)
        acc = np.zeros(epi_idx.size, dtype=np.float64)
        for gene in ifn_present:
            counts = np.asarray(Xq[index[gene], epi_idx].todense()).ravel()
            acc += np.log1p(1e4 * counts / scale)
        sum_ifn = float(acc.mean() * epi_idx.size)  # sum of per-cell means
        # acc is sum of log genes; mean across genes:
        sum_ifn = float((acc / len(ifn_present)).sum())
        cldn4_epi = cldn4[epi_idx]
        n_pos = int((cldn4_epi > 0).sum())
        sum_cldn4_log = float(np.log1p(1e4 * cldn4_epi / scale).sum())

    return {
        "n_raw": n_cells_raw,
        "n_qc": n_qc,
        "n_epi": n_epi,
        "n_cldn4_pos_epi": n_pos,
        "sum_cldn4_log_epi": sum_cldn4_log,
        "sum_ifn_epi": sum_ifn,
        "n_ifn_genes": len(ifn_present),
        "n_tnk": int(tnk.sum()),
        "n_ptprc": int(ptprc.sum()),
        "n_epcam": int(epcam.sum()),
    }


def library_rows() -> list[dict]:
    """Explicit library map. Mouse ids are GEO biological replicates, not inferred from outcome."""
    e = EXTRACT
    rows = []

    def add(accession, library, mouse, group, kind, path, features=None, primary_mouse=True, compartment="gate", note=""):
        rows.append({
            "accession": accession,
            "library": library,
            "mouse": mouse,
            "group": group,
            "kind": kind,
            "path": str(path),
            "features": str(features) if features else "",
            "primary_mouse": primary_mouse,
            "compartment": compartment,
            "note": note,
        })

    # GSE149813: two mice, YFP- and YFP+ epithelial sorts pooled within mouse.
    add("GSE149813", "AM1", "mouse1", "Kras-G12D 7w YFP-", "mtx",
        e / "GSM4513594_AM1_matrix.mtx.gz", e / "GSM4513594_AM1_features.tsv.gz",
        compartment="epi_sorted", note="sorted lung epithelium")
    add("GSE149813", "AM2", "mouse1", "Kras-G12D 7w YFP+", "mtx",
        e / "GSM4513595_AM2_matrix.mtx.gz", e / "GSM4513595_AM2_features.tsv.gz",
        compartment="epi_sorted")
    add("GSE149813", "AM4", "mouse2", "Kras-G12D 7w YFP-", "mtx",
        e / "GSM4513596_AM4_matrix.mtx.gz", e / "GSM4513596_AM4_features.tsv.gz",
        compartment="epi_sorted")
    add("GSE149813", "AM5", "mouse2", "Kras-G12D 7w YFP+", "mtx",
        e / "GSM4513597_AM5_matrix.mtx.gz", e / "GSM4513597_AM5_features.tsv.gz",
        compartment="epi_sorted")

    # GSE188436: condition libraries, no mouse id and no replicate label.
    for gsm, name, group in [
        ("GSM5683062", "18385X1", "KPF1F2_pre"),
        ("GSM5683063", "18836X3", "KP_pre"),
        ("GSM5683064", "18836X4", "KP_post"),
        ("GSM5683065", "18836X5", "KPF1F2_post"),
    ]:
        add("GSE188436", name, name, group, "mtx",
            e / f"{gsm}_{name}_matrix.mtx.gz", e / f"{gsm}_{name}_features.tsv.gz",
            primary_mouse=False, note="GEO sample is a condition, not a named mouse")

    # GSE264739: genotype replicates are biological mice.
    for gsm, name, group in [
        ("GSM8226902", "KP1", "KP"),
        ("GSM8226903", "KP2", "KP"),
        ("GSM8226904", "KP3", "KP"),
        ("GSM8226905", "KPP1", "KPP"),
        ("GSM8226906", "KPP2", "KPP"),
        ("GSM8226907", "KPP3", "KPP"),
    ]:
        add("GSE264739", name, name, group, "mtx",
            e / f"{gsm}_{name}_matrix.mtx.gz", e / f"{gsm}_{name}_features.tsv.gz",
            note="GEO replicate of autochthonous KP or KPP lung tumor")

    # GSE281744: each library pools 2 mice.
    add("GSE281744", "KRAS_21wkON", "pool_ON", "21wk_ON", "h5",
        e / "GSM8627458_KRAS_21wkON_filtered_feature_bc_matrix.h5",
        primary_mouse=False, note="GEO number of mice pooled = 2")
    add("GSE281744", "KRAS_20wkON_1wkOFF", "pool_OFF", "20wkON_1wkOFF", "h5",
        e / "GSM8627459_KRAS_20wkON_1wkOFF_filtered_feature_bc_matrix.h5",
        primary_mouse=False, note="GEO number of mice pooled = 2")

    # GSE281964: pre/post depletion libraries, mouse id not given.
    for gsm, name, group in [
        ("GSM8633704", "KP_pre", "KP_pre"),
        ("GSM8633705", "KP_post", "KP_post"),
        ("GSM8633706", "KPH_pre", "KPH_pre"),
        ("GSM8633707", "KPH_post", "KPH_post"),
    ]:
        add("GSE281964", name, name, group, "mtx",
            e / f"{gsm}_{name}_matrix.mtx.gz", e / f"{gsm}_{name}_features.tsv.gz",
            primary_mouse=False, note="pre/post library without a mouse id")

    # GSE317576: condition libraries, mice not identified.
    add("GSE317576", "PBS", "PBS", "tumor_free_PBS", "h5",
        e / "GSM9474402_sample_filtered_feature_bc_matrix_-_PBS.h5",
        primary_mouse=False, note="condition library; mouse replicates not in GEO")
    add("GSE317576", "PR8", "PR8", "tumor_free_PR8", "h5",
        e / "GSM9474402_sample_filtered_feature_bc_matrix_-_PR8.h5",
        primary_mouse=False, note="condition library")
    add("GSE317576", "PBS-T6", "PBS-T6", "PBS_T6", "h5",
        e / "GSM9474404_filtered_feature_bc_matrix_-_PBS-T6.h5",
        primary_mouse=False, note="condition library")
    add("GSE317576", "PR8-T6", "PR8-T6", "IAV_T6", "h5",
        e / "GSM9474405_filtered_feature_bc_matrix_-_PR8-T6.h5",
        primary_mouse=False, note="condition library")
    add("GSE317576", "PBS-T12", "PBS-T12", "PBS_T12", "h5",
        e / "GSM9474406_filtered_feature_bc_matrix_-_PBS-T12.h5",
        primary_mouse=False, note="condition library")
    add("GSE317576", "PR8-T12", "PR8-T12", "IAV_T12", "h5",
        e / "GSM9474407_filtered_feature_bc_matrix_-_PR8-T12.h5",
        primary_mouse=False, note="condition library")

    # GSE319598: filename prefix is the mouse; tumors from that mouse are pooled.
    # Stromal libraries are mouse-pooled and are not the tumor-cell unit.
    for gsm, name, mouse, group in [
        ("GSM9521009", "2731_T1", "2731", "vehicle"),
        ("GSM9521010", "2731_T2", "2731", "vehicle"),
        ("GSM9521011", "2731_Tmix", "2731", "vehicle"),
        ("GSM9521012", "2733_Tmix", "2733", "vehicle"),
        ("GSM9521013", "2868_Tmix", "2868", "sotorasib"),
        ("GSM9521014", "2873_T2", "2873", "sotorasib"),
        ("GSM9521015", "2873_Tmix", "2873", "sotorasib"),
    ]:
        add("GSE319598", name, mouse, group, "mtx",
            e / f"{gsm}_{name}_matrix.mtx.gz", e / f"{gsm}_{name}_features.tsv.gz",
            compartment="epi_sorted", note="eGFP+ tumor cells; tumors collapsed to filename mouse id")
    for gsm, name, group in [
        ("GSM9521016", "Stromal_control", "vehicle_pool"),
        ("GSM9521017", "Stromal_Treated", "sotorasib_pool"),
    ]:
        add("GSE319598", name, name, group, "mtx",
            e / f"{gsm}_{name}_matrix.mtx.gz", e / f"{gsm}_{name}_features.tsv.gz",
            primary_mouse=False, compartment="stroma_pool",
            note="GEO: stromal cells pooled across mice")

    # GSE338088: control vs knockout replicates are six mice, not paired fractions.
    for gsm, name, group in [
        ("GSM9866384", "T1_NG", "MCT2_control"),
        ("GSM9866385", "T1_Pos", "MCT2_knockout"),
        ("GSM9866386", "T2_NG", "MCT2_control"),
        ("GSM9866387", "T2_Pos", "MCT2_knockout"),
        ("GSM9866388", "T3_NG", "MCT2_control"),
        ("GSM9866389", "T3_Pos", "MCT2_knockout"),
    ]:
        add("GSE338088", name, name, group, "mtx",
            e / f"{gsm}_{name}.matrix.mtx.gz", e / f"{gsm}_{name}.features.tsv.gz",
            note="snRNA biological replicate")

    return rows


def score_all() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lib_path = OUT / "new_library_scores.tsv"
    fields = [
        "accession", "library", "mouse", "group", "primary_mouse", "compartment",
        "note", "n_raw", "n_qc", "n_epi", "n_cldn4_pos_epi", "sum_cldn4_log_epi",
        "sum_ifn_epi", "n_ifn_genes", "n_tnk", "n_ptprc", "n_epcam", "error",
    ]
    lines = ["\t".join(fields)]
    for spec in library_rows():
        print(f"scoring {spec['accession']} {spec['library']}", flush=True)
        path = Path(spec["path"])
        rec = {k: spec.get(k, "") for k in fields}
        rec["error"] = ""
        if not path.exists():
            rec["error"] = "missing file"
            lines.append("\t".join(str(rec[k]) for k in fields))
            print("  MISSING", path)
            continue
        try:
            if spec["kind"] == "h5":
                X, genes = read_h5(path)
            else:
                X, genes = read_mtx(path, Path(spec["features"]))
            scored = score_matrix(X, genes)
            del X
        except Exception as exc:  # noqa: BLE001 — record and continue
            scored = {"error": f"{type(exc).__name__}: {exc}"}
        if scored.get("error"):
            rec["error"] = scored["error"]
        for key, value in scored.items():
            if key in rec:
                rec[key] = value
        lines.append("\t".join("" if rec[k] is None else str(rec[k]) for k in fields))
        lib_path.write_text("\n".join(lines) + "\n")
        print("  ", {k: scored.get(k) for k in ("n_qc", "n_epi", "n_cldn4_pos_epi", "n_tnk", "n_ptprc", "error")}, flush=True)
    print("wrote", lib_path)


def _append_rows(recs: list[dict]) -> None:
    import csv

    path = OUT / "new_library_scores.tsv"
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        existing = list(reader)
    keys = {(r["accession"], r["library"]) for r in existing}
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for rec in existing:
            if rec.get("error"):
                continue
            writer.writerow(rec)
        for rec in recs:
            if (rec["accession"], rec["library"]) in keys and not rec.get("error"):
                continue
            writer.writerow({k: rec.get(k, "") for k in fields})


def score_remaining() -> None:
    """GSE338088 uses dotted 10x names. GSE295824 is gzip(bzip2(tar))."""
    import bz2
    import gzip
    import io
    import tarfile

    recs = []
    fields_spec = [s for s in library_rows() if s["accession"] == "GSE338088"]
    for spec in fields_spec:
        print("scoring", spec["accession"], spec["library"], flush=True)
        path = Path(spec["path"])
        rec = {**spec, "error": ""}
        X, genes = read_mtx(path, Path(spec["features"]))
        scored = score_matrix(X, genes)
        del X
        rec.update(scored)
        recs.append(rec)
        print(" ", {k: scored.get(k) for k in ("n_qc", "n_epi", "n_cldn4_pos_epi", "n_tnk", "error")}, flush=True)

    root = Path("/tmp/geo_dl/gse295824")
    for archive in sorted(root.glob("GSM*.tar.gz")):
        stem = archive.name.replace(".tar.gz", "")
        # GSM8958231_LN2A_M3 -> group LN2A, mouse M3
        parts = stem.split("_", 1)[1].split("_")
        group, mouse = parts[0], parts[1]
        print("scoring GSE295824", mouse, group, flush=True)
        raw = gzip.open(archive, "rb").read()
        payload = bz2.decompress(raw)
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:") as bundle:
            bundle.extractall(path=root / mouse)
        mtx = root / mouse / "matrix.mtx.gz"
        feat = root / mouse / "features.tsv.gz"
        X, genes = read_mtx(mtx, feat)
        scored = score_matrix(X, genes)
        del X
        rec = {
            "accession": "GSE295824",
            "library": mouse,
            "mouse": mouse,
            "group": group,
            "primary_mouse": True,
            "compartment": "gate",
            "note": "Sox2 GEMM lung tumor; mouse id is in the archive name",
            "error": scored.get("error", ""),
        }
        rec.update(scored)
        recs.append(rec)
        print(" ", {k: scored.get(k) for k in ("n_qc", "n_epi", "n_cldn4_pos_epi", "n_tnk", "n_ptprc", "error")}, flush=True)
        for child in (root / mouse).glob("*"):
            child.unlink()
        (root / mouse).rmdir()
    _append_rows(recs)
    print("appended", len(recs))


if __name__ == "__main__":
    import sys

    if "--remaining" in sys.argv:
        score_remaining()
    else:
        score_all()
