#!/usr/bin/env python3
"""CLDN4-only AUCell / TF-target proxy on public GSE148071 epithelium.

Question (patient-level): do CLDN4-high putative-malignant cells show
different IFN / MHC / TJ regulon activity vs CLDN4-low?

This is NOT full pySCENIC. cisTarget motif rankings are multi-GB and
were not downloaded. Regulons are public curated TF–target edges
(TRRUST v2 + DoRothEA + CollecTRI) scored with an Aibar-style AUCell
recovery curve. Companion gene-set AUCell (ISG / MHC-I APM / TJ, CLDN4
held out) is reported as pathway activity, not a TF regulon.

CLDN4 is the only splitter. TACSTD2 is never a gate (no dual-high).
GSE207422 and GSE131907 are not used.

Primary n = patients with ≥20 epithelial (putative malignant) cells and
≥8 cells in each within-patient CLDN4 Q4 / Q1 tail. Cell-level tests
are exploratory (pseudoreplication).
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from regulons import (  # noqa: E402
    ALL_TFS,
    CTRL_TFS,
    FOCAL,
    GENESET_MODULES,
    HOLD_OUT,
    IFN_TFS,
    LINEAGE_MARKERS,
    MHC_TFS,
    NORMAL_LUNG,
    TF_PROGRAM,
    TJ_TFS,
    all_always_genes,
)

MIN_UMI = 200
MIN_MAL = 20
MIN_TAIL = 8
AUC_THR = 0.05
MIN_REGULON = 5
NL_Q = 0.75
SEED = 1


def find_data(data_dir: Path) -> Path:
    if data_dir is None:
        cands = [
            ROOT / "data",
            Path("/tmp/gse148071_scenic/data"),
            Path("data/gse148071_scenic"),
        ]
    else:
        cands = [data_dir]
    for p in cands:
        if p and (any(p.rglob("*_exp.txt.gz")) or (p / "GSE148071_RAW.tar").exists()):
            return p
    raise SystemExit("GSE148071 matrices not found. Run scripts/download.py")


def list_exp_files(data_dir: Path) -> list[Path]:
    files = sorted(data_dir.rglob("*_exp.txt.gz"))
    if files:
        return files
    tar_path = data_dir / "GSE148071_RAW.tar"
    dest = data_dir / "GSE148071_files"
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    files = sorted(dest.rglob("*_exp.txt.gz"))
    if not files:
        raise SystemExit(f"no *_exp.txt.gz under {data_dir}")
    return files


def patient_from_filename(path: Path) -> str:
    m = re.search(r"_(P\d+)_", path.name)
    if m:
        return m.group(1)
    stem = path.name.replace(".txt.gz", "")
    parts = stem.split("_")
    return parts[1] if len(parts) >= 2 else stem


def parse_series_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    opener = gzip.open if str(path).endswith(".gz") else open
    titles = accs = None
    chars: list[list[str]] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                break
    n = len(titles or accs or [])
    if n == 0:
        return pd.DataFrame()
    rows = []
    for i in range(n):
        rec = {
            "sample": (titles[i] if titles else f"S{i}"),
            "geo_accession": accs[i] if accs else "",
        }
        for row in chars:
            if i >= len(row):
                continue
            val = row[i]
            if ":" in val:
                k, v = val.split(":", 1)
                rec[k.strip().lower().replace(" ", "_")] = v.strip()
        rows.append(rec)
    return pd.DataFrame(rows)


def load_priors(path: Path) -> dict[str, list[str]]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return {}
    out: dict[str, set[str]] = {}
    for r in df.itertuples(index=False):
        tf = str(r.tf)
        tgt = str(r.target)
        if tf not in ALL_TFS:
            continue
        if tgt in HOLD_OUT or tgt == tf:
            continue
        out.setdefault(tf, set()).add(tgt)
    return {k: sorted(v) for k, v in out.items()}


def aucell(expr: np.ndarray, member: np.ndarray, auc_threshold: float = AUC_THR) -> np.ndarray:
    """Aibar-style AUCell: normalized recovery of set members in the top ranks.

    Ranking is among the genes present in `expr` (this matrix / sample).
    This is not the R AUCell binary and not pySCENIC CLI.
    """
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


def module_mean(log_cp: np.ndarray, gene_index: dict[str, int], genes: list[str]) -> np.ndarray:
    idx = [gene_index[g] for g in genes if g in gene_index]
    if not idx:
        return np.full(log_cp.shape[0], np.nan)
    return log_cp[:, idx].mean(axis=1)


def assign_lineage(log_cp: np.ndarray, gene_index: dict[str, int]) -> np.ndarray:
    names = list(LINEAGE_MARKERS)
    scores = []
    for name in names:
        scores.append(module_mean(log_cp, gene_index, LINEAGE_MARKERS[name]))
    mat = np.vstack(scores)
    best = np.nanargmax(np.nan_to_num(mat, nan=-np.inf), axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    labels = np.array(names, dtype=object)[best].copy()
    labels[~np.isfinite(top) | (top < 0.12)] = "Unassigned"
    return labels


def stream_matrix(path: Path):
    """Load one GSE148071 GSM*_P*_exp.txt.gz as genes x cells (float32)."""
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        header = [h.strip().strip('"') for h in header if h != ""]
        if header and header[0] in {"", "gene", "Gene", "index", "Index", "GENE"}:
            cell_ids = header[1:]
        else:
            cell_ids = header
        n = len(cell_ids)
        genes: list[str] = []
        rows: list[np.ndarray] = []
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
            genes.append(gene)
            rows.append(arr)
    X = np.vstack(rows) if rows else np.zeros((0, n), dtype=np.float32)
    return cell_ids, genes, X


def spear(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return dict(rho=np.nan, p=np.nan, n=n)
    r, p = stats.spearmanr(x[m], y[m])
    return dict(rho=float(r), p=float(p), n=n)


def paired_wilcoxon(high, low):
    high = np.asarray(high, float)
    low = np.asarray(low, float)
    m = np.isfinite(high) & np.isfinite(low)
    h, l = high[m], low[m]
    n = int(m.sum())
    if n < 3:
        return dict(n=n, high_mean=np.nan, low_mean=np.nan, delta_median=np.nan,
                    n_high_gt=0, W=np.nan, p=np.nan)
    d = h - l
    # exact Wilcoxon when n is small and no ties in |d|
    try:
        res = stats.wilcoxon(h, l, zero_method="wilcox", alternative="two-sided", method="exact")
        W, p = float(res.statistic), float(res.pvalue)
    except ValueError:
        res = stats.wilcoxon(h, l, zero_method="pratt", alternative="two-sided", method="auto")
        W, p = float(res.statistic), float(res.pvalue)
    return dict(
        n=n,
        high_mean=float(np.mean(h)),
        low_mean=float(np.mean(l)),
        high_median=float(np.median(h)),
        low_median=float(np.median(l)),
        delta_median=float(np.median(d)),
        n_high_gt=int((d > 0).sum()),
        W=W,
        p=p,
    )


def build_regulons(priors: dict[str, list[str]], genes_present: set[str]) -> dict[str, list[str]]:
    regs: dict[str, list[str]] = {}
    for tf, targets in priors.items():
        keep = [g for g in targets if g in genes_present and g not in HOLD_OUT]
        if len(keep) >= MIN_REGULON:
            regs[f"prior_{tf}"] = keep
    # Combined program unions (unique targets)
    for name, tfs in (("IFN", IFN_TFS), ("MHC", MHC_TFS), ("TJ", TJ_TFS), ("CTRL", CTRL_TFS)):
        union: list[str] = []
        seen = set()
        for tf in tfs:
            for g in priors.get(tf, []):
                if g in genes_present and g not in HOLD_OUT and g not in seen:
                    seen.add(g)
                    union.append(g)
        if len(union) >= MIN_REGULON:
            regs[f"prior_{name}_union"] = union
    # IRF1 also contributes to MHC presentation; extra union
    mhc_plus = list(regs.get("prior_MHC_union", []))
    extra = [g for g in priors.get("IRF1", []) if g in genes_present and g not in HOLD_OUT]
    for g in extra:
        if g not in mhc_plus:
            mhc_plus.append(g)
    if len(mhc_plus) >= MIN_REGULON:
        regs["prior_MHC_IRF1_union"] = mhc_plus
    return regs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=ROOT / "results")
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    figdir = out / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    data_dir = find_data(args.data)
    files = list_exp_files(data_dir)
    print(f"[data] {len(files)} matrices under {data_dir}", flush=True)

    series = parse_series_matrix(data_dir / "GSE148071_series_matrix.txt.gz")
    if series.empty:
        # sometimes extracted next to tar
        sm = list(data_dir.rglob("GSE148071_series_matrix.txt.gz"))
        if sm:
            series = parse_series_matrix(sm[0])
    if not series.empty:
        series.to_csv(out / "sample_metadata.tsv", sep="\t", index=False)

    prior_path = ROOT / "resources" / "tf_targets_public.tsv"
    priors = load_priors(prior_path)
    print(f"[priors] TFs with edges: {sorted(priors)}", flush=True)
    for tf, tgts in sorted(priors.items()):
        print(f"  {tf}: {len(tgts)} targets", flush=True)

    always = set(all_always_genes())
    for tgts in priors.values():
        always.update(tgts)

    cell_rows = []
    lineage_rows = []
    file_audit = []
    genes_seen: set[str] = set()

    for i, path in enumerate(files, 1):
        patient = patient_from_filename(path)
        print(f"[{i}/{len(files)}] {patient} {path.name}", flush=True)
        cell_ids, genes, X = stream_matrix(path)
        genes_seen.update(genes)
        n_umi = X.sum(axis=0)
        keep = n_umi >= MIN_UMI
        n_drop = int((~keep).sum())
        X = X[:, keep]
        cell_ids = [c for c, k in zip(cell_ids, keep) if k]
        n_umi = n_umi[keep]
        if X.shape[1] == 0:
            file_audit.append(dict(patient=patient, file=path.name, n_cells=0, n_qc=0, n_epi=0))
            continue
        tot = n_umi.astype(np.float64)
        tot[tot <= 0] = np.nan
        log_cp = np.log1p((X.astype(np.float64) / tot[None, :]) * 1e4).T  # cells x genes
        gene_index = {g: j for j, g in enumerate(genes)}
        labels = assign_lineage(log_cp, gene_index)
        nl = module_mean(log_cp, gene_index, NORMAL_LUNG)
        cldn4 = log_cp[:, gene_index["CLDN4"]] if "CLDN4" in gene_index else np.full(log_cp.shape[0], np.nan)
        tac = log_cp[:, gene_index["TACSTD2"]] if "TACSTD2" in gene_index else np.full(log_cp.shape[0], np.nan)

        vc = pd.Series(labels).value_counts()
        lineage_rows.append(
            dict(patient=patient, n_qc=int(log_cp.shape[0]), n_dropped_umi=n_drop,
                 **{f"n_{k}": int(vc.get(k, 0)) for k in list(LINEAGE_MARKERS) + ["Unassigned"]})
        )

        epi = labels == "Epithelial"
        n_epi = int(epi.sum())
        file_audit.append(dict(patient=patient, file=path.name, n_cells=len(keep),
                               n_qc=int(keep.sum()) if hasattr(keep, "sum") else int(log_cp.shape[0]),
                               n_epi=n_epi, n_genes=len(genes)))
        if n_epi == 0:
            continue

        # AUCell ranking universe: genes detected in ≥1% of this sample's epithelium,
        # plus always-keep (TFs, prior targets, modules). Caps noise from dropouts.
        epi_X = X[:, epi].T  # cells x genes (raw counts for detection)
        det = (epi_X > 0).mean(axis=0)
        keep_g = (det >= 0.01) | np.array([g in always for g in genes])
        # drop all-zero columns
        keep_g = keep_g & (epi_X.sum(axis=0) > 0)
        expr = log_cp[epi][:, keep_g]
        g_keep = [g for g, k in zip(genes, keep_g) if k]
        g_index = {g: j for j, g in enumerate(g_keep)}

        regs = build_regulons(priors, set(g_keep))
        scores = {}
        for rname, members in regs.items():
            mask = np.array([g in set(members) for g in g_keep])
            scores[rname] = aucell(expr, mask)
        for mname, members in GENESET_MODULES.items():
            mem = [g for g in members if g not in HOLD_OUT]
            mask = np.array([g in set(mem) for g in g_keep])
            scores[mname] = aucell(expr, mask)
            scores[mname.replace("gs_", "mean_")] = module_mean(log_cp[epi], gene_index, mem)

        # TF RNA itself (companion, not the regulon)
        for tf in ALL_TFS:
            scores[f"rna_{tf}"] = (
                log_cp[epi, gene_index[tf]] if tf in gene_index else np.full(n_epi, np.nan)
            )

        epi_ids = [c for c, e in zip(cell_ids, epi) if e]
        epi_nl = nl[epi]
        epi_cldn4 = cldn4[epi]
        epi_tac = tac[epi]
        epi_umi = n_umi[epi]
        # malignant-like: epithelial AND normal-lung score ≤ sample epithelial 75th pct
        nl_cut = np.nanquantile(epi_nl, NL_Q) if np.isfinite(epi_nl).any() else np.inf
        mal_like = epi_nl <= nl_cut

        for j in range(n_epi):
            rec = dict(
                patient=patient,
                barcode=epi_ids[j],
                n_umi=float(epi_umi[j]),
                cldn4=float(epi_cldn4[j]),
                tacstd2=float(epi_tac[j]),
                normal_lung=float(epi_nl[j]),
                mal_like=bool(mal_like[j]),
                n_rank_genes=int(expr.shape[1]),
            )
            for k, v in scores.items():
                rec[k] = float(v[j]) if np.isfinite(v[j]) else np.nan
            cell_rows.append(rec)

    cells = pd.DataFrame(cell_rows)
    lineage = pd.DataFrame(lineage_rows)
    audit = pd.DataFrame(file_audit)
    lineage.to_csv(out / "lineage_counts.tsv", sep="\t", index=False)
    audit.to_csv(out / "file_audit.tsv", sep="\t", index=False)
    cells.to_csv(out / "epithelial_cell_scores.tsv.gz", sep="\t", index=False)
    print(f"[cells] epithelial rows={len(cells)} patients={cells.patient.nunique()}", flush=True)

    # Regulon membership table
    regs_all = build_regulons(priors, genes_seen)
    mem_rows = []
    for rname, members in regs_all.items():
        for g in members:
            mem_rows.append(dict(regulon=rname, target=g, n_targets=len(members)))
    for mname, members in GENESET_MODULES.items():
        mem = [g for g in members if g not in HOLD_OUT]
        for g in mem:
            mem_rows.append(dict(regulon=mname, target=g, n_targets=len(mem)))
    pd.DataFrame(mem_rows).to_csv(out / "regulon_membership.tsv", sep="\t", index=False)
    reg_sizes = {k: len(v) for k, v in regs_all.items()}
    for mname, members in GENESET_MODULES.items():
        reg_sizes[mname] = len([g for g in members if g not in HOLD_OUT])
    print("[regulon sizes]", json.dumps(reg_sizes, indent=2), flush=True)

    score_cols = [
        c for c in cells.columns
        if c.startswith("prior_") or c.startswith("gs_") or c.startswith("mean_")
        or c.startswith("rna_") or c in {"cldn4", "tacstd2", "n_umi", "normal_lung"}
    ]

    def patient_table(df: pd.DataFrame, label: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        rows = []
        paired = []
        elig = []
        for pat, sub in df.groupby("patient"):
            n = len(sub)
            q = sub.cldn4.quantile([0.25, 0.5, 0.75])
            low = sub[sub.cldn4 <= q[0.25]]
            high = sub[sub.cldn4 >= q[0.75]]
            ok = n >= MIN_MAL and len(low) >= MIN_TAIL and len(high) >= MIN_TAIL
            rec = dict(
                split=label,
                patient=pat,
                n_cells=n,
                n_q1=len(low),
                n_q4=len(high),
                eligible=bool(ok),
                cldn4_q1=float(q[0.25]),
                cldn4_q2=float(q[0.50]),
                cldn4_q3=float(q[0.75]),
                cldn4_mean=float(sub.cldn4.mean()),
                tacstd2_mean=float(sub.tacstd2.mean()),
            )
            for col in score_cols:
                rec[f"{col}_mean"] = float(sub[col].mean()) if col in sub else np.nan
                if ok:
                    rec[f"{col}_q4"] = float(high[col].mean()) if col in high else np.nan
                    rec[f"{col}_q1"] = float(low[col].mean()) if col in low else np.nan
                    rec[f"{col}_delta"] = rec[f"{col}_q4"] - rec[f"{col}_q1"]
            rows.append(rec)
            elig.append(dict(split=label, patient=pat, n_cells=n, n_q1=len(low),
                             n_q4=len(high), eligible=bool(ok)))
            if ok:
                prow = dict(split=label, patient=pat, n_cells=n, n_q1=len(low), n_q4=len(high))
                for col in score_cols:
                    prow[f"{col}_q4"] = rec.get(f"{col}_q4", np.nan)
                    prow[f"{col}_q1"] = rec.get(f"{col}_q1", np.nan)
                    prow[f"{col}_delta"] = rec.get(f"{col}_delta", np.nan)
                paired.append(prow)
        return pd.DataFrame(rows), pd.DataFrame(paired), pd.DataFrame(elig)

    per_epi, paired_epi, elig_epi = patient_table(cells, "epithelial")
    mal = cells[cells.mal_like].copy()
    per_mal, paired_mal, elig_mal = patient_table(mal, "mal_like")

    per_patient = pd.concat([per_epi, per_mal], ignore_index=True)
    paired = pd.concat([paired_epi, paired_mal], ignore_index=True)
    eligibility = pd.concat([elig_epi, elig_mal], ignore_index=True)
    per_patient.to_csv(out / "per_patient.tsv", sep="\t", index=False)
    paired.to_csv(out / "paired_high_low.tsv", sep="\t", index=False)
    eligibility.to_csv(out / "eligibility.tsv", sep="\t", index=False)

    # Primary tests: mal_like if ≥6 eligible, else epithelial
    n_mal = int(elig_mal.eligible.sum()) if len(elig_mal) else 0
    n_epi = int(elig_epi.eligible.sum()) if len(elig_epi) else 0
    primary_split = "mal_like" if n_mal >= 6 else "epithelial"
    paired_p = paired[paired.split == primary_split].copy()
    print(f"[n] epithelial eligible={n_epi}  mal_like eligible={n_mal}  primary={primary_split}", flush=True)

    test_cols = [c for c in score_cols if not c.startswith("rna_")]
    tests = []
    for col in test_cols + ["cldn4", "tacstd2", "n_umi"]:
        if f"{col}_q4" not in paired_p.columns:
            continue
        w = paired_wilcoxon(paired_p[f"{col}_q4"], paired_p[f"{col}_q1"])
        tests.append(dict(split=primary_split, feature=col, layer="paired_Q4_vs_Q1", **w))
    # sample-level Spearman of patient-mean CLDN4 vs feature (eligible patients)
    per_p = per_patient[(per_patient.split == primary_split) & (per_patient.eligible)]
    for col in test_cols:
        key = f"{col}_mean"
        if key not in per_p.columns:
            continue
        s = spear(per_p["cldn4_mean"], per_p[key])
        tests.append(dict(split=primary_split, feature=col, layer="sample_spearman_cldn4",
                          n=s["n"], high_mean=np.nan, low_mean=np.nan,
                          delta_median=s["rho"], n_high_gt=np.nan, W=np.nan, p=s["p"],
                          rho=s["rho"]))
    # delta vs delta (TJ union vs IFN union)
    if {"prior_TJ_union_delta", "prior_IFN_union_delta"}.issubset(paired_p.columns):
        s = spear(paired_p["prior_TJ_union_delta"], paired_p["prior_IFN_union_delta"])
        tests.append(dict(split=primary_split, feature="delta_TJ_vs_delta_IFN",
                          layer="paired_delta_spearman", n=s["n"],
                          high_mean=np.nan, low_mean=np.nan, delta_median=s["rho"],
                          n_high_gt=np.nan, W=np.nan, p=s["p"], rho=s["rho"]))
    if {"gs_tj_no_cldn4_delta", "gs_ifn_isg_delta"}.issubset(paired_p.columns):
        s = spear(paired_p["gs_tj_no_cldn4_delta"], paired_p["gs_ifn_isg_delta"])
        tests.append(dict(split=primary_split, feature="delta_gsTJ_vs_delta_gsIFN",
                          layer="paired_delta_spearman", n=s["n"],
                          high_mean=np.nan, low_mean=np.nan, delta_median=s["rho"],
                          n_high_gt=np.nan, W=np.nan, p=s["p"], rho=s["rho"]))

    # cell-level exploratory
    df_p = mal if primary_split == "mal_like" else cells
    for col in test_cols:
        if col not in df_p.columns:
            continue
        s = spear(df_p["cldn4"], df_p[col])
        tests.append(dict(split=primary_split, feature=col, layer="cell_spearman_exploratory",
                          n=s["n"], high_mean=np.nan, low_mean=np.nan,
                          delta_median=s["rho"], n_high_gt=np.nan, W=np.nan, p=s["p"],
                          rho=s["rho"]))

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(out / "tests.tsv", sep="\t", index=False)

    # Compact primary table
    primary_features = [
        "prior_IFN_union", "prior_MHC_union", "prior_MHC_IRF1_union", "prior_TJ_union",
        "prior_STAT1", "prior_IRF1", "prior_IRF7", "prior_NLRC5", "prior_CIITA",
        "prior_GRHL2", "prior_ELF3", "prior_KLF4", "prior_OVOL1", "prior_HIF1A",
        "gs_ifn_isg", "gs_mhc1_apm", "gs_tj_no_cldn4", "gs_oxphos",
        "mean_ifn_isg", "mean_mhc1_apm", "mean_tj_no_cldn4", "mean_oxphos",
        "tacstd2", "n_umi",
    ]
    prim_rows = []
    for feat in primary_features:
        hit = tests_df[(tests_df.feature == feat) & (tests_df.layer == "paired_Q4_vs_Q1")]
        if hit.empty:
            continue
        r = hit.iloc[0]
        prim_rows.append(dict(
            feature=feat,
            n_patients=int(r.n),
            high_mean=r.high_mean,
            low_mean=r.low_mean,
            delta_median=r.delta_median,
            n_high_gt_low=int(r.n_high_gt) if pd.notna(r.n_high_gt) else np.nan,
            W=r.W,
            p=r.p,
            method="AUCell_TF_target_proxy" if feat.startswith("prior_") or feat.startswith("gs_")
            else ("mean_log1p_cp10k" if feat.startswith("mean_") else "companion"),
        ))
    primary = pd.DataFrame(prim_rows)
    primary.to_csv(out / "primary_table.tsv", sep="\t", index=False)

    # Honest n table
    n_table = pd.DataFrame([
        dict(item="patients_deposited", n=42, note="Wu et al. Nat Commun 2021 PMID 33953163"),
        dict(item="exp_matrices", n=len(files), note="GEO GSE148071_RAW.tar GSM*_P*_exp.txt.gz"),
        dict(item="epithelial_cells", n=int(len(cells)), note="marker-argmax; putative malignant"),
        dict(item="mal_like_cells", n=int(cells.mal_like.sum()),
             note=f"epithelial AND normal-lung score ≤ sample Q{NL_Q}"),
        dict(item="epithelial_eligible_patients", n=n_epi,
             note=f"≥{MIN_MAL} epi and ≥{MIN_TAIL} cells per CLDN4 Q1/Q4"),
        dict(item="mal_like_eligible_patients", n=n_mal,
             note=f"≥{MIN_MAL} mal-like and ≥{MIN_TAIL} cells per CLDN4 Q1/Q4"),
        dict(item="primary_n", n=int(len(paired_p)), note=f"split={primary_split}; unit of inference"),
        dict(item="dropped_patients_primary", n=int(42 - len(paired_p)),
             note="below occupancy floor or missing CLDN4 tails"),
    ])
    n_table.to_csv(out / "n_table.tsv", sep="\t", index=False)

    # Figures
    def _paired_plot(feats, title, fname):
        fig, axes = plt.subplots(1, len(feats), figsize=(3.2 * len(feats), 4.2), squeeze=False)
        for ax, feat in zip(axes[0], feats):
            if f"{feat}_q4" not in paired_p.columns:
                ax.set_visible(False)
                continue
            ys = []
            for _, row in paired_p.iterrows():
                ax.plot([0, 1], [row[f"{feat}_q1"], row[f"{feat}_q4"]], color="#888", lw=0.8, alpha=0.7)
                ax.scatter([0, 1], [row[f"{feat}_q1"], row[f"{feat}_q4"]], s=22, color="#1f4e79", zorder=3)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["CLDN4 Q1", "CLDN4 Q4"])
            ax.set_title(feat, fontsize=9)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        fig.suptitle(title, fontsize=11)
        fig.tight_layout()
        fig.savefig(figdir / f"{fname}.png", dpi=140)
        fig.savefig(figdir / f"{fname}.pdf")
        plt.close(fig)

    _paired_plot(
        [c for c in ["prior_IFN_union", "prior_MHC_union", "prior_TJ_union", "gs_oxphos"] if f"{c}_q4" in paired_p.columns],
        f"GSE148071 {primary_split} AUCell  CLDN4 Q4 vs Q1  n={len(paired_p)}",
        "fig_paired_regulons",
    )
    _paired_plot(
        [c for c in ["gs_ifn_isg", "gs_mhc1_apm", "gs_tj_no_cldn4", "gs_oxphos"] if f"{c}_q4" in paired_p.columns],
        f"GSE148071 {primary_split} gene-set AUCell  n={len(paired_p)}",
        "fig_paired_genesets",
    )

    # occupancy
    fig, ax = plt.subplots(figsize=(8, 3.6))
    el = elig_epi.sort_values("n_cells", ascending=False)
    colors = ["#1f4e79" if e else "#c4c4c4" for e in el.eligible]
    ax.bar(el.patient, el.n_cells, color=colors)
    ax.axhline(MIN_MAL, color="#b00", ls="--", lw=0.8, label=f"floor n={MIN_MAL}")
    ax.set_ylabel("epithelial cells")
    ax.set_title("GSE148071 honest n — epithelial occupancy (blue = eligible)")
    ax.tick_params(axis="x", labelrotation=90, labelsize=7)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figdir / "fig_honest_n.png", dpi=140)
    fig.savefig(figdir / "fig_honest_n.pdf")
    plt.close(fig)

    if {"prior_TJ_union_delta", "prior_IFN_union_delta"}.issubset(paired_p.columns):
        fig, ax = plt.subplots(figsize=(4.2, 4.0))
        ax.scatter(paired_p["prior_TJ_union_delta"], paired_p["prior_IFN_union_delta"], s=36, color="#1f4e79")
        for _, r in paired_p.iterrows():
            ax.annotate(r.patient, (r["prior_TJ_union_delta"], r["prior_IFN_union_delta"]),
                        fontsize=7, xytext=(3, 3), textcoords="offset points")
        ax.axhline(0, color="#aaa", lw=0.6)
        ax.axvline(0, color="#aaa", lw=0.6)
        ax.set_xlabel("Δ AUCell prior_TJ_union (Q4−Q1)")
        ax.set_ylabel("Δ AUCell prior_IFN_union (Q4−Q1)")
        ax.set_title(f"same-patient deltas  n={len(paired_p)}")
        fig.tight_layout()
        fig.savefig(figdir / "fig_delta_scatter.png", dpi=140)
        fig.savefig(figdir / "fig_delta_scatter.pdf")
        plt.close(fig)

    # pyscenic availability note
    pyscenic_note = "not_imported"
    try:
        import pyscenic  # type: ignore

        pyscenic_note = f"importable:{getattr(pyscenic, '__version__', 'unknown')}; full grn/ctx/aucell NOT run"
    except Exception as exc:
        pyscenic_note = f"not_available:{type(exc).__name__}"

    summary = {
        "accession": "GSE148071",
        "pmid": "33953163",
        "method": "AUCell + public TF-target priors (TRRUST/DoRothEA/CollecTRI). Not full pySCENIC.",
        "pyscenic": pyscenic_note,
        "splitter": "CLDN4 only (within-patient Q4 vs Q1). No dual-high. TACSTD2 never a gate.",
        "not_used": ["GSE207422", "GSE131907"],
        "primary_split": primary_split,
        "n_patients_deposited": 42,
        "n_epithelial": int(len(cells)),
        "n_mal_like": int(cells.mal_like.sum()),
        "n_eligible_epithelial": n_epi,
        "n_eligible_mal_like": n_mal,
        "n_primary": int(len(paired_p)),
        "min_mal": MIN_MAL,
        "min_tail": MIN_TAIL,
        "auc_threshold": AUC_THR,
        "hold_out": sorted(HOLD_OUT),
        "regulon_sizes": reg_sizes,
        "prior_tfs": {k: len(v) for k, v in priors.items()},
        "eligible_patients": paired_p.patient.tolist(),
        "primary_table": primary.to_dict(orient="records"),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps({k: summary[k] for k in (
        "method", "pyscenic", "primary_split", "n_primary",
        "n_eligible_epithelial", "n_eligible_mal_like", "eligible_patients"
    )}, indent=2))
    print("[done]", out)


if __name__ == "__main__":
    main()
