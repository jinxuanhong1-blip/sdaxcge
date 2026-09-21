#!/usr/bin/env python3
"""TACSTD2–T/NK attenuation given CLDN4 versus epithelial negative controls.

Concordant-4 only: GSE123902, GSE131907, GSE205335, GSE189357 (n = 65).
The malignant gate and the T/NK definition match the locked CLDN4 % positive
analysis. TACSTD2 and CLDN3 are scores, not gates.

Question
--------
Does the patient-level Spearman association between malignant TACSTD2 and the
T/NK fraction move closer to zero when CLDN4 is partialled out than when
CLDN3, CLDN7, EPCAM, MUC1, or KRT19 is partialled out?

Primary score: malignant percent positive (the locked CLDN4 scale).
Sensitivity: malignant mean log1p.
Primary subset: all four cohorts.
Sensitivity: author-label tumors only (GSE131907 + GSE205335), because the
marker gate uses EPCAM and KRT19.
"""

from __future__ import annotations

import csv
import gzip
import os
import subprocess
import tarfile
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from stats import bh_fdr, dl_meta, partial_spearman, spearman, spearman_p

ROOT = Path(__file__).resolve().parent
GEO = Path(os.environ.get("GEO_DIR", "/tmp/geo_c4"))
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
LOCKED = ROOT / "data" / "locked_pr644_units.tsv"

SCORE_GENES = ["TACSTD2", "CLDN3", "CLDN4", "CLDN7", "EPCAM", "MUC1", "KRT19", "KRT8", "KRT18"]
CONDITIONERS = ["CLDN4", "CLDN3", "CLDN7", "EPCAM", "MUC1", "KRT19"]
CONTROLS = ["CLDN3", "CLDN7", "EPCAM", "MUC1", "KRT19"]
MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
KERATIN = ["KRT8", "KRT18", "KRT19"]
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
N_PERM = 10000
SEED = 1
ELIG_131_ORIGIN = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
MIN_MAL_131 = 20


def say(msg: str) -> None:
    print(msg, flush=True)


def _zeros(n: int) -> np.ndarray:
    return np.zeros(n, dtype=np.float64)


def _gene_cols_from_header(header: list[str], genes: list[str]) -> dict[str, int]:
    first = {}
    for i, name in enumerate(header):
        key = name.strip().upper()
        if key and key not in first:
            first[key] = i
    return {g: first[g] for g in genes if g in first}


def score_from_arrays(counts: dict[str, np.ndarray], malignant: np.ndarray) -> dict:
    idx = np.where(malignant)[0]
    n = int(idx.size)
    out = {"n_malignant_scored": n}
    logs = {}
    for g in SCORE_GENES:
        if n == 0 or g not in counts:
            out[f"mal_{g}_pct"] = float("nan")
            out[f"mal_{g}_mean"] = float("nan")
            continue
        v = counts[g][idx].astype(float)
        out[f"mal_{g}_pct"] = 100.0 * float(np.mean(v > 0))
        logs[g] = np.log1p(v)
        out[f"mal_{g}_mean"] = float(np.mean(logs[g]))
    if n and all(g in logs for g in KERATIN):
        out["mal_KRT_mean"] = float(np.mean((logs["KRT8"] + logs["KRT18"] + logs["KRT19"]) / 3.0))
    else:
        out["mal_KRT_mean"] = float("nan")
    return out


def load_dense_selected(path: Path, genes: list[str]) -> dict[str, np.ndarray]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split(",")
    cols = _gene_cols_from_header(header, genes)
    missing = [g for g in genes if g not in cols]
    usecols = sorted(set(cols.values()))
    arr = np.loadtxt(gzip.open(path, "rt"), delimiter=",", skiprows=1, usecols=usecols)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    pos = {c: i for i, c in enumerate(usecols)}
    n = arr.shape[0]
    out = {}
    for g in genes:
        out[g] = _zeros(n) if g not in cols else arr[:, pos[cols[g]]]
    out["_missing"] = missing
    out["_n"] = n
    return out


def marker_masks(counts: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    n = int(counts["_n"])

    def g(name):
        return counts[name] if name in counts else _zeros(n)

    mal = ((g("EPCAM") > 0) | (g("KRT8") > 0) | (g("KRT18") > 0) | (g("KRT19") > 0)) & (g("PTPRC") == 0)
    tnk = (
        (g("CD3D") > 0) | (g("CD3E") > 0) | (g("CD8A") > 0) | (g("NKG7") > 0) | (g("GNLY") > 0) | (g("KLRD1") > 0)
    ) & (~mal)
    return mal, tnk


def score_gse123902() -> list[dict]:
    say("GSE123902")
    folder = GEO / "gse123902"
    if not any(folder.glob("*.csv.gz")):
        tar = GEO / "GSE123902_RAW.tar"
        if not tar.exists():
            raise SystemExit(f"missing {tar}")
        folder.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar) as tf:
            tf.extractall(folder)
    rows = []
    for path in sorted(folder.glob("*_dense.csv.gz")):
        name = path.name
        if "_NORMAL_" in name:
            continue
        if "_PRIMARY_" in name:
            tissue = "PRIMARY"
        elif "_METASTASIS_" in name:
            tissue = "METASTASIS"
        else:
            continue
        donor = name.split("_")[2]
        counts = load_dense_selected(path, sorted(set(MARKERS + SCORE_GENES)))
        missing_score = [g for g in SCORE_GENES if g in counts.get("_missing", [])]
        if missing_score:
            raise RuntimeError(f"{name} missing score genes {missing_score}")
        mal, tnk = marker_masks(counts)
        rec = {
            "dataset": "GSE123902",
            "unit_id": donor,
            "unit_type": "donor",
            "tissue": tissue,
            "n_cells": int(counts["_n"]),
            "n_malignant": int(mal.sum()),
            "n_tnk": int(tnk.sum()),
            "frac_tnk": float(np.mean(tnk)),
            "malig_def": "marker_malig",
        }
        rec.update(score_from_arrays(counts, mal))
        rows.append(rec)
        say(
            f"  {donor} {tissue} cells={rec['n_cells']} mal={rec['n_malignant']} "
            f"tnk={rec['n_tnk']} CLDN4%={rec['mal_CLDN4_pct']:.2f} TACSTD2%={rec['mal_TACSTD2_pct']:.2f}"
        )
    rows.sort(key=lambda r: (r["unit_id"], r["tissue"]))
    kept = []
    seen = set()
    for rec in rows:
        if rec["unit_id"] in seen:
            continue
        seen.add(rec["unit_id"])
        if rec["n_malignant"] >= 20 and rec["n_tnk"] >= 20:
            kept.append(rec)
    say(f"  kept {len(kept)} donors")
    return kept


def _read_10x_selected(prefix: Path, genes: list[str]) -> dict[str, np.ndarray]:
    symbols = []
    with gzip.open(str(prefix) + "_features.tsv.gz", "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            symbols.append(parts[1].upper() if len(parts) > 1 else parts[0].upper())
    first_row = {}
    for i, sym in enumerate(symbols):
        if sym not in first_row:
            first_row[sym] = i + 1
    wanted = {first_row[g]: g for g in genes if g in first_row}
    missing = [g for g in genes if g not in first_row]
    n_cells = 0
    with gzip.open(str(prefix) + "_barcodes.tsv.gz", "rt") as handle:
        for _ in handle:
            n_cells += 1
    counts = {g: _zeros(n_cells) for g in genes}
    mtx = str(prefix) + "_matrix.mtx.gz"
    with gzip.open(mtx, "rt") as handle:
        n_gene = n_col = None
        for line in handle:
            if line.startswith("%"):
                continue
            if n_gene is None:
                n_gene, n_col, _nnz = map(int, line.split())
                if n_col != n_cells:
                    raise RuntimeError(f"{mtx} columns {n_col} != barcodes {n_cells}")
                continue
            sp = line.find(" ")
            if sp < 0:
                continue
            gene = wanted.get(int(line[:sp]))
            if gene is None:
                continue
            c_s, v_s = line[sp + 1 :].split()
            counts[gene][int(c_s) - 1] = float(v_s)
    counts["_missing"] = missing
    counts["_n"] = n_cells
    return counts


def score_gse189357() -> list[dict]:
    say("GSE189357")
    folder = GEO / "gse189357"
    if not any(folder.glob("*_matrix.mtx.gz")):
        tar = GEO / "GSE189357_RAW.tar"
        if not tar.exists():
            raise SystemExit(f"missing {tar}")
        folder.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar) as tf:
            tf.extractall(folder)
    rows = []
    for pat in [f"TD{i}" for i in range(1, 10)]:
        hits = list(folder.glob(f"*_{pat}_matrix.mtx.gz"))
        if len(hits) != 1:
            raise RuntimeError(f"expected one matrix for {pat}, found {hits}")
        prefix = Path(str(hits[0]).replace("_matrix.mtx.gz", ""))
        counts = _read_10x_selected(prefix, sorted(set(MARKERS + SCORE_GENES)))
        missing_score = [g for g in SCORE_GENES if g in counts.get("_missing", [])]
        if missing_score:
            raise RuntimeError(f"{pat} missing score genes {missing_score}")
        mal, tnk = marker_masks(counts)
        rec = {
            "dataset": "GSE189357",
            "unit_id": pat,
            "unit_type": "patient",
            "tissue": "TUMOR",
            "n_cells": int(counts["_n"]),
            "n_malignant": int(mal.sum()),
            "n_tnk": int(tnk.sum()),
            "frac_tnk": float(np.mean(tnk)),
            "malig_def": "marker_malig",
        }
        rec.update(score_from_arrays(counts, mal))
        rows.append(rec)
        say(
            f"  {pat} cells={rec['n_cells']} mal={rec['n_malignant']} tnk={rec['n_tnk']} "
            f"CLDN4%={rec['mal_CLDN4_pct']:.2f}"
        )
    return rows


def _parse_soft(path: Path) -> dict[str, dict]:
    code_map = {}
    cur = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = {}
            elif cur is None:
                continue
            elif line.startswith("!Sample_title = "):
                cur["title"] = line.split("= ", 1)[1]
            elif line.startswith("!Sample_characteristics_ch1 = "):
                val = line.split("= ", 1)[1]
                if ": " in val:
                    k, v = val.split(": ", 1)
                    cur[k] = v
            elif line.startswith("!Sample_platform_id") and cur.get("title"):
                title = cur["title"]
                code = title.split(" ", 1)[1] if " " in title else ""
                code = code.upper().replace("-", "_")
                if code:
                    code_map[code] = cur
                cur = None
    return code_map


def _norm_code(orig: str) -> str:
    x = orig.upper().replace("-", "_")
    if x.endswith("_3P") or x.endswith("_5P"):
        x = x[:-3]
    return x


def score_gse205335() -> list[dict]:
    say("GSE205335")
    gene_path = GEO / "GSE205335_negctrl_genes.tsv.gz"
    if not gene_path.exists():
        rscript = ROOT / "extract_gse205335.R"
        subprocess.check_call(["Rscript", str(rscript), str(GEO)])
    ident_path = GEO / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = GEO / "GSE205335_family.soft.gz"
    code_map = _parse_soft(soft_path)
    by_patient = defaultdict(lambda: {"n": 0, "mal": 0, "tnk": 0, "tissues": set(), "mal_bc": []})
    n_unmapped = 0
    with gzip.open(ident_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            code = _norm_code(parts[idx["orig.ident"]])
            meta = code_map.get(code)
            if meta is None:
                n_unmapped += 1
                continue
            tissue = meta.get("tissue", "")
            if tissue.startswith("Normal"):
                continue
            pat = meta.get("patient", "")
            rec = by_patient[pat]
            rec["n"] += 1
            rec["tissues"].add(tissue)
            mal = parts[idx["lineage.sub"]] == "Malignant cells"
            tnk = parts[idx["lineage.total"]] == "T/NK cells"
            if mal:
                rec["mal"] += 1
                rec["mal_bc"].append(parts[idx["barcode"]])
            if tnk:
                rec["tnk"] += 1
    say(f"  unmapped identity rows {n_unmapped}; patients with tumor tissue {len(by_patient)}")
    gene_of = {}
    with gzip.open(gene_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        gidx = {name: i for i, name in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            gene_of[parts[0]] = parts
    say(f"  gene table barcodes {len(gene_of)}")
    rows = []
    for pat in sorted(by_patient):
        info = by_patient[pat]
        mal_bc = [b for b in info["mal_bc"] if b in gene_of]
        n_scored = len(mal_bc)
        counts = {g: np.empty(n_scored, dtype=float) for g in SCORE_GENES}
        for i, b in enumerate(mal_bc):
            parts = gene_of[b]
            for g in SCORE_GENES:
                counts[g][i] = float(parts[gidx[g]])
        rec = {
            "dataset": "GSE205335",
            "unit_id": pat,
            "unit_type": "patient",
            "tissue": ",".join(sorted(info["tissues"])),
            "n_cells": info["n"],
            "n_malignant": info["mal"],
            "n_tnk": info["tnk"],
            "frac_tnk": info["tnk"] / info["n"] if info["n"] else float("nan"),
            "malig_def": "author_malig",
        }
        rec.update(score_from_arrays(counts, np.ones(n_scored, dtype=bool)))
        rows.append(rec)
        say(
            f"  {pat} cells={rec['n_cells']} mal={rec['n_malignant']} "
            f"scored={rec['n_malignant_scored']} CLDN4%={rec['mal_CLDN4_pct']:.2f}"
        )
    return rows


def score_gse131907() -> list[dict]:
    say("GSE131907")
    ann_path = GEO / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = GEO / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    cells = []
    with gzip.open(ann_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in handle:
            p = line.rstrip("\n").split("\t")
            cells.append(
                {
                    "index": p[idx["Index"]],
                    "sample": p[idx["Sample"]],
                    "origin": p[idx["Sample_Origin"]],
                    "cell_type": p[idx["Cell_type"]],
                    "cell_subtype": p[idx["Cell_subtype"]],
                }
            )
    by_sample = defaultdict(list)
    for rec in cells:
        by_sample[rec["sample"]].append(rec)
    eligible = []
    for sample, recs in sorted(by_sample.items()):
        origin = recs[0]["origin"]
        n_mal = sum(r["cell_subtype"] == "Malignant cells" for r in recs)
        n_tnk = sum(r["cell_type"] in {"T lymphocytes", "NK cells"} for r in recs)
        if origin in ELIG_131_ORIGIN and n_mal >= MIN_MAL_131:
            eligible.append((sample, origin, len(recs), n_mal, n_tnk))
    say(f"  eligible samples {len(eligible)}")
    mal_index = {}
    for sample, *_ in eligible:
        for r in by_sample[sample]:
            if r["cell_subtype"] == "Malignant cells":
                mal_index[r["index"]] = sample
    say(f"  malignant barcodes {len(mal_index)}")
    targets = set(SCORE_GENES)
    per_sample_vals = {s: {g: [] for g in SCORE_GENES} for s, *_ in eligible}
    with gzip.open(mat_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        keep_cols = []
        order_sample = []
        for i, b in enumerate(barcodes):
            sample = mal_index.get(b)
            if sample is not None:
                keep_cols.append(i)
                order_sample.append(sample)
        say(f"  matrix cells {len(barcodes)} malignant columns matched {len(keep_cols)}")
        found = set()
        n_gene = 0
        for line in handle:
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab].upper()
            n_gene += 1
            if gene not in targets or gene in found:
                if n_gene % 5000 == 0:
                    say(f"  streamed {n_gene} genes, hits {len(found)}")
                continue
            found.add(gene)
            vals = line.rstrip("\n").split("\t")
            buckets = {s: [] for s, *_ in eligible}
            for col, sample in zip(keep_cols, order_sample):
                raw = vals[col + 1]
                fv = 0.0 if raw in {"", "0", "0.0"} else float(raw)
                buckets[sample].append(fv)
            for sample, arr in buckets.items():
                per_sample_vals[sample][gene] = arr
            say(f"  hit {gene} at gene {n_gene}")
        missing = targets - found
        if missing:
            raise RuntimeError(f"GSE131907 missing genes {sorted(missing)}")
    rows = []
    for sample, origin, n_cells, n_mal, n_tnk in eligible:
        counts = {g: np.asarray(per_sample_vals[sample][g], dtype=float) for g in SCORE_GENES}
        n_scored = counts["CLDN4"].size
        rec = {
            "dataset": "GSE131907",
            "unit_id": sample,
            "unit_type": "sample",
            "tissue": origin,
            "n_cells": n_cells,
            "n_malignant": n_mal,
            "n_tnk": n_tnk,
            "frac_tnk": n_tnk / n_cells if n_cells else float("nan"),
            "malig_def": "author_malig",
        }
        rec.update(score_from_arrays(counts, np.ones(n_scored, dtype=bool)))
        rows.append(rec)
        say(f"  {sample} mal={n_mal} scored={rec['n_malignant_scored']} CLDN4%={rec['mal_CLDN4_pct']:.2f}")
    return rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    keys = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = {}
            for k in keys:
                v = row[k]
                if isinstance(v, float):
                    out[k] = f"{v:.8g}"
                else:
                    out[k] = v
            writer.writerow(out)
    say(f"wrote {path}")


def calibrate(units: list[dict]) -> None:
    locked = {}
    with LOCKED.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            locked[(row["dataset"], row["unit_id"])] = row
    keys = [
        "frac_tnk",
        "n_malignant",
        "n_tnk",
        "mal_CLDN4_pct",
        "mal_CLDN7_pct",
        "mal_EPCAM_pct",
        "mal_MUC1_pct",
        "mal_KRT19_pct",
        "mal_KRT8_pct",
        "mal_KRT18_pct",
    ]
    max_diff = {k: 0.0 for k in keys}
    n_match = 0
    missing = []
    for u in units:
        key = (u["dataset"], u["unit_id"])
        ref = locked.get(key)
        if ref is None:
            missing.append(key)
            continue
        n_match += 1
        for k in keys:
            diff = abs(float(u[k]) - float(ref[k]))
            if diff > max_diff[k]:
                max_diff[k] = diff
    extra = sorted(set(locked) - {(u["dataset"], u["unit_id"]) for u in units})
    say(f"calibration matched {n_match} locked units; missing {missing}; extra {extra}")
    for k, v in max_diff.items():
        say(f"  max |Δ| {k} = {v:.6g}")
    if n_match != 65 or missing or extra:
        raise SystemExit("unit set does not match the locked n=65 table")
    if max_diff["mal_CLDN4_pct"] > 0.05 or max_diff["frac_tnk"] > 1e-4:
        raise SystemExit("CLDN4 % or T/NK fraction drifted from the locked table")
    rows = [{"field": k, "max_abs_diff": v, "n_matched": n_match} for k, v in max_diff.items()]
    write_tsv(TAB / "calibration_vs_pr644.tsv", rows)


def _vectors(units: list[dict], kind: str):
    y = np.array([u["frac_tnk"] for u in units], dtype=float)
    x = np.array([u[f"mal_TACSTD2_{kind}"] for u in units], dtype=float)
    cond = {g: np.array([u[f"mal_{g}_{kind}"] for u in units], dtype=float) for g in CONDITIONERS}
    cohorts = np.array([u["dataset"] for u in units])
    return x, y, cond, cohorts


def _cohort_slices(cohorts: np.ndarray, allowed: list[str]):
    slices = []
    for ds in allowed:
        m = np.where(cohorts == ds)[0]
        if m.size >= 5:
            slices.append((ds, m))
    return slices


def cohort_partials(x, y, z, slices, k_cov: int):
    names, rhos, ns, ps = [], [], [], []
    for ds, m in slices:
        rho = partial_spearman(x[m], y[m], [] if z is None else [z[m]])
        n = int(m.size)
        names.append(ds)
        rhos.append(rho)
        ns.append(n)
        ps.append(spearman_p(rho, n, k_cov))
    return names, rhos, ns, ps


def label_swap_p(x, y, z_a, z_b, slices, rng, n_perm: int) -> dict:
    """Exchangeability of two conditioners. Positive delta means A's partial rho is larger."""

    def meta_delta(za, zb):
        ra, na = [], []
        rb, nb = [], []
        for _ds, m in slices:
            ra.append(partial_spearman(x[m], y[m], [za[m]]))
            rb.append(partial_spearman(x[m], y[m], [zb[m]]))
            na.append(int(m.size))
            nb.append(int(m.size))
        ma = dl_meta(ra, na, k_cov=1)
        mb = dl_meta(rb, nb, k_cov=1)
        return ma["rho"] - mb["rho"], ma["rho"], mb["rho"]

    obs, rho_a, rho_b = meta_delta(z_a, z_b)
    count = 0
    za = np.array(z_a, dtype=float, copy=True)
    zb = np.array(z_b, dtype=float, copy=True)
    for _ in range(n_perm):
        a = za.copy()
        b = zb.copy()
        for _ds, m in slices:
            swap = rng.random(m.size) < 0.5
            if not np.any(swap):
                continue
            ii = m[swap]
            a[ii], b[ii] = b[ii], a[ii].copy()
        delta, _, _ = meta_delta(a, b)
        if np.isfinite(delta) and abs(delta) >= abs(obs) - 1e-15:
            count += 1
    return {
        "delta_partial_a_minus_b": obs,
        "meta_partial_a": rho_a,
        "meta_partial_b": rho_b,
        "perm_p_two_sided": (count + 1) / (n_perm + 1),
        "n_perm": n_perm,
    }


def conditioner_shuffle_p(x, y, z, slices, rng, n_perm: int) -> dict:
    """Shuffle the conditioner within cohort. Tests attenuation against a null conditioner."""

    def meta_partial(zz):
        rhos, ns = [], []
        for _ds, m in slices:
            rhos.append(partial_spearman(x[m], y[m], [zz[m]]))
            ns.append(int(m.size))
        return dl_meta(rhos, ns, k_cov=1)["rho"]

    r0_list, n0 = [], []
    for _ds, m in slices:
        r0_list.append(spearman(x[m], y[m]))
        n0.append(int(m.size))
    r0 = dl_meta(r0_list, n0, k_cov=0)["rho"]
    obs_partial = meta_partial(z)
    obs_att = obs_partial - r0
    count = 0
    base = np.asarray(z, dtype=float)
    for _ in range(n_perm):
        zz = base.copy()
        for _ds, m in slices:
            zz[m] = base[m][rng.permutation(m.size)]
        att = meta_partial(zz) - r0
        if np.isfinite(att) and abs(att) >= abs(obs_att) - 1e-15:
            count += 1
    return {
        "meta_unadj": r0,
        "meta_partial": obs_partial,
        "attenuation": obs_att,
        "perm_p_attenuation": (count + 1) / (n_perm + 1),
        "n_perm": n_perm,
    }


def analyze(units: list[dict]) -> None:
    rng = np.random.default_rng(SEED)
    cohort_rows = []
    meta_rows = []
    coli_rows = []
    head_rows = []
    subsets = {
        "all4": COHORTS,
        "author": ["GSE131907", "GSE205335"],
    }
    for kind in ("pct", "mean"):
        x, y, cond, cohorts = _vectors(units, kind)
        for subset, allowed in subsets.items():
            slices = _cohort_slices(cohorts, allowed)
            say(f"stats {kind} {subset} cohorts={[ds for ds, _ in slices]}")
            # unadjusted TACSTD2
            names, rhos, ns, ps = cohort_partials(x, y, None, slices, 0)
            meta = dl_meta(rhos, ns, 0)
            meta_rows.append(
                {
                    "subset": subset,
                    "score": kind,
                    "conditioner": "none",
                    "k_cov": 0,
                    "rho_partial": meta["rho"],
                    "p_partial": meta["p"],
                    "I2_partial": meta["I2"],
                    "ci_lo": meta["ci_lo"],
                    "ci_hi": meta["ci_hi"],
                    "rho_unadj": meta["rho"],
                    "p_unadj": meta["p"],
                    "attenuation": 0.0,
                    "percent_attenuated": 0.0,
                    "k": meta["k"],
                    "N": meta["N"],
                    "perm_p_attenuation": float("nan"),
                }
            )
            for ds, rho, n, p in zip(names, rhos, ns, ps):
                cohort_rows.append(
                    {
                        "subset": subset,
                        "score": kind,
                        "dataset": ds,
                        "conditioner": "none",
                        "n": n,
                        "rho_unadj": rho,
                        "p_unadj": p,
                        "rho_partial": rho,
                        "p_partial": p,
                        "attenuation": 0.0,
                        "percent_attenuated": 0.0,
                    }
                )
            unadj = {ds: rho for ds, rho in zip(names, rhos)}
            unadj_n = {ds: n for ds, n in zip(names, ns)}
            # each conditioner
            for gene in CONDITIONERS:
                say(f"  conditioner {gene}")
                names, rhos, ns, ps = cohort_partials(x, y, cond[gene], slices, 1)
                meta_p = dl_meta(rhos, ns, 1)
                meta_u = dl_meta([unadj[ds] for ds in names], [unadj_n[ds] for ds in names], 0)
                att = meta_p["rho"] - meta_u["rho"]
                if abs(meta_u["rho"]) < 1e-8:
                    pct = float("nan")
                else:
                    pct = 100.0 * (meta_u["rho"] - meta_p["rho"]) / meta_u["rho"]
                shuf = conditioner_shuffle_p(x, y, cond[gene], slices, rng, N_PERM)
                meta_rows.append(
                    {
                        "subset": subset,
                        "score": kind,
                        "conditioner": gene,
                        "k_cov": 1,
                        "rho_partial": meta_p["rho"],
                        "p_partial": meta_p["p"],
                        "I2_partial": meta_p["I2"],
                        "ci_lo": meta_p["ci_lo"],
                        "ci_hi": meta_p["ci_hi"],
                        "rho_unadj": meta_u["rho"],
                        "p_unadj": meta_u["p"],
                        "attenuation": att,
                        "percent_attenuated": pct,
                        "k": meta_p["k"],
                        "N": meta_p["N"],
                        "perm_p_attenuation": shuf["perm_p_attenuation"],
                    }
                )
                for ds, rho, n, p in zip(names, rhos, ns, ps):
                    r0 = unadj[ds]
                    cohort_rows.append(
                        {
                            "subset": subset,
                            "score": kind,
                            "dataset": ds,
                            "conditioner": gene,
                            "n": n,
                            "rho_unadj": r0,
                            "p_unadj": spearman_p(r0, n, 0),
                            "rho_partial": rho,
                            "p_partial": p,
                            "attenuation": rho - r0,
                            "percent_attenuated": float("nan") if abs(r0) < 1e-8 else 100.0 * (r0 - rho) / r0,
                        }
                    )
                # collinearity TACSTD2 vs conditioner
                cr, cn = [], []
                for ds, m in slices:
                    cr.append(spearman(x[m], cond[gene][m]))
                    cn.append(int(m.size))
                cm = dl_meta(cr, cn, 0)
                coli_rows.append(
                    {
                        "subset": subset,
                        "score": kind,
                        "gene": gene,
                        "rho_with_TACSTD2": cm["rho"],
                        "p": cm["p"],
                        "I2": cm["I2"],
                        "ci_lo": cm["ci_lo"],
                        "ci_hi": cm["ci_hi"],
                        "k": cm["k"],
                        "N": cm["N"],
                    }
                )
            for control in CONTROLS:
                say(f"  label-swap CLDN4 vs {control}")
                swap = label_swap_p(x, y, cond["CLDN4"], cond[control], slices, rng, N_PERM)
                # attenuation contrast = partial_CLDN4 - partial_control
                # positive => CLDN4 leaves the less-negative (more attenuated, if unadj < 0) residual
                head_rows.append(
                    {
                        "subset": subset,
                        "score": kind,
                        "index": "CLDN4",
                        "control": control,
                        "meta_partial_CLDN4": swap["meta_partial_a"],
                        "meta_partial_control": swap["meta_partial_b"],
                        "delta_partial_CLDN4_minus_control": swap["delta_partial_a_minus_b"],
                        "perm_p_two_sided": swap["perm_p_two_sided"],
                        "n_perm": swap["n_perm"],
                        "N": int(sum(m.size for _ds, m in slices)),
                    }
                )
    # BH within each subset/score across the five controls
    for subset in subsets:
        for kind in ("pct", "mean"):
            ix = [i for i, r in enumerate(head_rows) if r["subset"] == subset and r["score"] == kind]
            q = bh_fdr([head_rows[i]["perm_p_two_sided"] for i in ix])
            for i, qq in zip(ix, q):
                head_rows[i]["q_bh"] = qq
    write_tsv(TAB / "c4_cohort_attenuation.tsv", cohort_rows)
    write_tsv(TAB / "c4_meta_attenuation.tsv", meta_rows)
    write_tsv(TAB / "c4_collinearity.tsv", coli_rows)
    write_tsv(TAB / "c4_head_to_head.tsv", head_rows)
    plot_c4(meta_rows, cohort_rows)


def plot_c4(meta_rows: list[dict], cohort_rows: list[dict]) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    order = ["none"] + CONDITIONERS
    labels = ["unadjusted"] + CONDITIONERS
    sub = [r for r in meta_rows if r["subset"] == "all4" and r["score"] == "pct"]
    by = {r["conditioner"]: r for r in sub}
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4))
    ax = axes[0]
    ypos = np.arange(len(order))[::-1]
    rhos = [by[g]["rho_partial"] for g in order]
    lo = [by[g]["rho_partial"] - by[g]["ci_lo"] for g in order]
    hi = [by[g]["ci_hi"] - by[g]["rho_partial"] for g in order]
    colors = ["#4d4d4d"] + ["#b2182b" if g == "CLDN4" else "#2166ac" for g in CONDITIONERS]
    ax.errorbar(rhos, ypos, xerr=[lo, hi], fmt="none", ecolor="#666666", elinewidth=1, capsize=2)
    ax.scatter(rhos, ypos, c=colors, s=36, zorder=3)
    ax.axvline(0, color="#888888", lw=0.6)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("DL Spearman ρ  TACSTD2 %pos vs T/NK")
    ax.set_title("Concordant-4 partial ρ")
    ax = axes[1]
    genes = CONDITIONERS
    ypos = np.arange(len(genes))[::-1]
    att = [by[g]["attenuation"] for g in genes]
    colors = ["#b2182b" if g == "CLDN4" else "#2166ac" for g in genes]
    ax.axvline(0, color="#888888", lw=0.6)
    ax.scatter(att, ypos, c=colors, s=36, zorder=3)
    # cohort dots
    for i, g in enumerate(genes):
        xs = [
            r["attenuation"]
            for r in cohort_rows
            if r["subset"] == "all4" and r["score"] == "pct" and r["conditioner"] == g
        ]
        ax.scatter(xs, np.full(len(xs), ypos[i]), s=12, c="#999999", zorder=2)
    ax.set_yticks(ypos)
    ax.set_yticklabels(genes)
    ax.set_xlabel("Change in ρ  (partial − unadjusted)")
    ax.set_title("Concordant-4 shift")
    fig.tight_layout()
    fig.savefig(FIG / "c4_attenuation.png", dpi=160)
    fig.savefig(FIG / "c4_attenuation.pdf")
    plt.close(fig)
    say(f"wrote {FIG / 'c4_attenuation.png'}")


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    units = []
    units.extend(score_gse123902())
    units.extend(score_gse189357())
    units.extend(score_gse205335())
    units.extend(score_gse131907())
    say(f"units {len(units)}")
    counts = defaultdict(int)
    for u in units:
        counts[u["dataset"]] += 1
    say("cohort n " + ", ".join(f"{k}={counts[k]}" for k in COHORTS))
    write_tsv(TAB / "c4_patient_units.tsv", units)
    calibrate(units)
    analyze(units)


if __name__ == "__main__":
    main()
