#!/usr/bin/env python3
"""Threshold sweep of CLDN4 and TACSTD2 vs response in public ICI NSCLC bulk.

Downloads public GEO supplements into data/public_ici_nsclc (gitignored).
Does not read EGA. Does not invent RECIST labels that are absent from the file.
"""

from __future__ import annotations

import gzip
import json
import math
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "public_ici_nsclc"
OUT = ROOT / "results" / "public_ici_nsclc_cldn4"
GENES = ["CLDN4", "TACSTD2"]
KRTS = ["KRT8", "KRT18", "KRT19"]
ENS = {
    "CLDN4": "ENSG00000189143",
    "TACSTD2": "ENSG00000184292",
    "KRT8": "ENSG00000170421",
    "KRT18": "ENSG00000111057",
    "KRT19": "ENSG00000171345",
}
QUANTILES = [0.25, 0.33, 0.50, 0.67, 0.75]


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print("GET", url)
    urllib.request.urlretrieve(url, dest)
    return dest


def series_table(path: Path) -> pd.DataFrame:
    """GEO series matrix characteristics, one row per sample."""
    blocks: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            key, vals = parts[0], parts[1:]
            if key in ("!Sample_title", "!Sample_geo_accession", "!Sample_description"):
                blocks[key] = vals
                continue
            if key.startswith("!Sample_characteristics"):
                # One characteristic per line; key is the text before the first colon.
                name = vals[0].split(":", 1)[0].strip() if vals and ":" in vals[0] else key
                parsed = []
                for v in vals:
                    parsed.append(v.split(":", 1)[1].strip() if ":" in v else v.strip())
                # Duplicate names get a suffix.
                base = name
                i = 2
                while name in blocks:
                    name = f"{base}__{i}"
                    i += 1
                blocks[name] = parsed
    n = len(blocks["!Sample_title"])
    df = pd.DataFrame({k: v for k, v in blocks.items() if len(v) == n})
    df = df.rename(
        columns={
            "!Sample_title": "title",
            "!Sample_geo_accession": "gsm",
            "!Sample_description": "description",
        }
    )
    return df


def woolf(a: int, b: int, c: int, d: int) -> tuple[float, float, float, bool]:
    corrected = min(a, b, c, d) == 0
    if corrected:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    else:
        aa, bb, cc, dd = float(a), float(b), float(c), float(d)
    or_ = (aa * dd) / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    log_or = math.log(or_)
    return or_, math.exp(log_or - 1.96 * se), math.exp(log_or + 1.96 * se), corrected


def contrasts(values: np.ndarray, y: np.ndarray) -> list[dict]:
    rows = []
    mask = np.isfinite(values) & np.isfinite(y)
    v = values[mask]
    yy = y[mask].astype(int)
    if len(v) < 8 or yy.sum() < 2 or (len(yy) - yy.sum()) < 2:
        return rows
    specs = [("quantile", q, np.quantile(v, q)) for q in QUANTILES]
    q1, q3 = np.quantile(v, 0.33), np.quantile(v, 0.67)
    specs.append(("tertile_ends", 0.67, q3))
    for rule, q, cut in specs:
        if rule == "tertile_ends":
            keep = (v >= q3) | (v <= q1)
            if keep.sum() < 8:
                continue
            high = v >= q3
            vv, yk = v[keep], yy[keep]
            high = high[keep]
        else:
            # Ties at the cut go to high. Skip if the cut collapses a group.
            high = v >= cut
            vv, yk = v, yy
        n_high = int(high.sum())
        n_low = int((~high).sum())
        if n_high < 3 or n_low < 3:
            continue
        a = int(((high) & (yk == 1)).sum())
        b = int(((high) & (yk == 0)).sum())
        c = int(((~high) & (yk == 1)).sum())
        d = int(((~high) & (yk == 0)).sum())
        or_, lo, hi, corrected = woolf(a, b, c, d)
        _, p = fisher_exact([[a, b], [c, d]])
        rows.append(
            {
                "rule": "median" if rule == "quantile" and abs(q - 0.5) < 1e-9 else rule,
                "q": q,
                "n": int(len(yk)),
                "n_high": n_high,
                "n_low": n_low,
                "resp_high": a,
                "nonresp_high": b,
                "resp_low": c,
                "nonresp_low": d,
                "OR": or_,
                "OR_lo": lo,
                "OR_hi": hi,
                "haldane": corrected,
                "fisher_p": float(p),
            }
        )
    return rows


def log_expr(series: pd.Series, already_log: bool) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    if already_log:
        return x
    return np.log2(x.clip(lower=0) + 1)


def keratin_residual(expr: pd.DataFrame, gene: str) -> pd.Series | None:
    if gene not in expr.columns or any(k not in expr.columns for k in KRTS):
        return None
    cols = [gene] + KRTS
    sub = expr[cols].apply(pd.to_numeric, errors="coerce").dropna()
    if len(sub) < 8:
        return None
    y = sub[gene].to_numpy(float)
    X = np.column_stack([np.ones(len(sub)), sub[KRTS].to_numpy(float)])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = pd.Series(y - X @ beta, index=sub.index, name=gene)
    return resid


def add_gene_rows(bag: list, cohort: str, endpoint: str, expr: pd.DataFrame, y: pd.Series, note: str) -> None:
    y = y.dropna()
    common = expr.index.intersection(y.index)
    if len(common) < 8:
        raise SystemExit(
            f"{cohort} {endpoint}: only {len(common)} samples joined "
            f"(expr {list(expr.index[:3])} vs label {list(y.index[:3])})"
        )
    y = y.loc[common]
    expr = expr.loc[common]
    for gene in GENES:
        if gene not in expr.columns:
            bag.append(
                {
                    "cohort": cohort,
                    "endpoint": endpoint,
                    "gene": gene,
                    "adjust": "none",
                    "rule": "absent",
                    "note": note,
                }
            )
            continue
        versions = {
            "none": expr[gene],
            "keratin": keratin_residual(expr, gene),
        }
        for adjust, values in versions.items():
            if values is None:
                bag.append(
                    {
                        "cohort": cohort,
                        "endpoint": endpoint,
                        "gene": gene,
                        "adjust": adjust,
                        "rule": "no_keratin",
                        "note": note,
                    }
                )
                continue
            aligned = pd.concat([values.rename("v"), y.rename("y")], axis=1).dropna()
            for row in contrasts(aligned["v"].to_numpy(float), aligned["y"].to_numpy(float)):
                row.update(
                    {
                        "cohort": cohort,
                        "endpoint": endpoint,
                        "gene": gene,
                        "adjust": adjust,
                        "note": note,
                    }
                )
                bag.append(row)


def read_symbols_by_sample(path: Path, already_log: bool) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    df.index = df.index.astype(str)
    keep = [g for g in GENES + KRTS if g in df.index]
    # Exact symbols only. Substring hits such as KRT8P41 are not KRT8.
    sub = df.loc[keep].T
    sub.index = sub.index.astype(str)
    return sub.apply(lambda s: log_expr(s, already_log))


def read_ensembl_by_sample(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    df.index = df.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    rename = {ens: sym for sym, ens in ENS.items() if ens in df.index}
    sub = df.loc[list(rename)].T.rename(columns=rename)
    sub.index = sub.index.astype(str)
    return sub.apply(lambda s: log_expr(s, already_log=False))


def _eu_float(text: str) -> float:
    text = text.strip().replace(",", ".")
    if text == "":
        return np.nan
    return float(text)


def read_sample_by_gene_tpm(path: Path) -> pd.DataFrame:
    """Semicolon TPM. Sample id is either an extra leading field or sits under a blank header.

    France3 rows are one field longer than the gene header. France4 puts the sample id
    under a leading empty header and writes decimals with a comma.
    """
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split(";")
        first = handle.readline()
        if not first:
            raise SystemExit(f"empty TPM file {path}")
        parts = first.rstrip("\n").split(";")
        if len(parts) == len(header) + 1 and header[0] != "":
            genes = header
            def take(fields: list[str]) -> tuple[str, list[str]]:
                return fields[0], fields[1:]
        elif len(parts) == len(header) and header[0] == "":
            genes = header[1:]
            def take(fields: list[str]) -> tuple[str, list[str]]:
                return fields[0], fields[1:]
        else:
            raise SystemExit(f"unexpected width in {path}: {len(parts)} vs {len(header)}")
        ids = []
        rows = []
        sid, vals = take(parts)
        ids.append(sid)
        rows.append([_eu_float(x) for x in vals])
        for line in handle:
            fields = line.rstrip("\n").split(";")
            if len(fields) != len(parts):
                raise SystemExit(f"ragged row in {path}: {len(fields)}")
            sid, vals = take(fields)
            ids.append(sid)
            rows.append([_eu_float(x) for x in vals])
    wide = pd.DataFrame(rows, index=ids, columns=genes)
    keep = [g for g in GENES + KRTS if g in wide.columns]
    return wide[keep].apply(lambda s: log_expr(s, already_log=False))


def pfs_binary(time: pd.Series, event: pd.Series, cut: float) -> pd.Series:
    t = pd.to_numeric(time, errors="coerce")
    e = pd.to_numeric(event, errors="coerce")
    out = pd.Series(np.nan, index=t.index, dtype=float)
    out[(t >= cut)] = 1.0
    out[(t < cut) & (e == 1)] = 0.0
    # Early censor (t < cut, event 0) stays missing.
    return out


def load_gse126044() -> tuple[pd.DataFrame, pd.Series, str]:
    raw = pd.read_csv(DATA / "GSE126044_counts.txt.gz", sep="\t", index_col=0, compression="gzip")
    raw.index = raw.index.astype(str)
    raw = raw.apply(pd.to_numeric, errors="coerce")
    keep = [g for g in GENES + KRTS if g in raw.index]
    # CPM uses the full count library, not the five genes being scored.
    lib = raw.sum(axis=0).replace(0, np.nan)
    counts = raw.loc[keep].T
    expr = np.log2(counts.div(lib, axis=0) * 1e6 + 1)
    expr.index = expr.index.astype(str)
    ph = series_table(DATA / "GSE126044_series_matrix.txt.gz")
    ph["sample"] = ph["title"].str.replace("RNA-seq_", "", regex=False)
    y = ph.set_index("sample")["patient response"].map({"responder": 1.0, "non-responder": 0.0})
    note = "Pre-treatment NSCLC biopsies. Response is the GEO label responder/non-responder, not a RECIST string."
    return expr, y, note


def load_gse166449() -> tuple[pd.DataFrame, pd.Series, str]:
    expr = read_symbols_by_sample(DATA / "GSE166449_Raw_gene_TPM_matrix.txt.gz", already_log=True)
    ph = series_table(DATA / "GSE166449_series_matrix.txt.gz")
    ph["responder"] = ph["title"].str.contains("nonResponder", case=False).map({True: 0.0, False: 1.0})
    # Titles are Immunotherapy_Responder* or Immunotherapy_nonResponder*.
    # nonResponder contains Responder, so test nonResponder first.
    y = ph.set_index("description")["responder"]
    note = "Pre-treatment advanced lung cancer. Response is the GEO title (7 responders, 15 non-responders). Expression values are already on a log scale (CLDN4 max about 3) and were not logged again."
    return expr, y, note


def load_gse135222() -> tuple[pd.DataFrame, pd.Series, str, dict]:
    expr = read_ensembl_by_sample(DATA / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz")
    expr.index = expr.index.str.replace(" ", "", regex=False)
    ph = series_table(DATA / "GSE135222_series_matrix.txt.gz")
    ph["sample"] = ph["title"].str.replace(" ", "", regex=False)
    ph = ph.set_index("sample")
    time = pd.to_numeric(ph["pfs.time"], errors="coerce")
    event = pd.to_numeric(ph["progression-free survival (pfs)"], errors="coerce")
    audit = {
        "mean_time_event1": float(time[event == 1].mean()),
        "mean_time_event0": float(time[event == 0].mean()),
        "n_event1": int((event == 1).sum()),
        "n_event0": int((event == 0).sum()),
        "time_unit": "days, as deposited in pfs.time",
    }
    if not audit["mean_time_event1"] < audit["mean_time_event0"]:
        raise SystemExit("GSE135222 event coding is not shorter-time=1; refusing to binarize")
    y = pfs_binary(time, event, cut=180)
    y.index = expr.index.intersection(y.index).union(y.index)
    note = (
        "Tumor RNA, anti-PD-1/PD-L1. Public GEO has PFS time in days and a 0/1 flag, not RECIST. "
        "1 is treated as the event because those samples have the shorter times. "
        "Benefit = PFS time >= 180 days. Early censoring (time < 180 and flag 0) is excluded. "
        "This is not ORR and not DCR."
    )
    return expr, y, note, audit


def load_gse207422() -> list[tuple[str, pd.DataFrame, pd.Series, str]]:
    expr = read_symbols_by_sample(DATA / "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", already_log=True)
    meta = pd.read_excel(DATA / "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
    meta = meta.dropna(subset=["Sample"]).copy()
    meta["Sample"] = meta["Sample"].astype(str)
    meta = meta[meta["Sample"].isin(expr.index)]
    labeled = meta.dropna(subset=["Pathologic Response"]).copy()
    mpr = labeled.set_index("Sample")["Pathologic Response"].astype(str).str.startswith("MPR").astype(float)
    rec = meta.dropna(subset=["RECIST"]).copy()
    orr = rec.set_index("Sample")["RECIST"].isin(["CR", "PR"]).astype(float)
    base = "Neoadjuvant toripalimab plus platinum chemotherapy. Pre-treatment biopsies only. "
    return [
        (
            "MPR",
            expr,
            mpr,
            base + "MPR includes MPR and MPR (pCR) versus NMPR. Not RECIST ORR.",
        ),
        (
            "ORR",
            expr,
            orr,
            base + "ORR is CR or PR versus SD. The deposited RECIST column has no PD, so DCR is not a contrast.",
        ),
    ]


def load_france(expr_name: str, clinical: pd.DataFrame, id_col: str, time_col: str, event_col: str, cohort_note: str):
    expr = read_sample_by_gene_tpm(DATA / expr_name)
    clinical = clinical.set_index(id_col)
    time = pd.to_numeric(clinical[time_col], errors="coerce")
    event = pd.to_numeric(clinical[event_col], errors="coerce")
    audit = {
        "mean_time_event1": float(time[event == 1].mean()),
        "mean_time_event0": float(time[event == 0].mean()),
        "time_min": float(time.min()),
        "time_max": float(time.max()),
        "n": int(time.notna().sum()),
    }
    if not (audit["mean_time_event1"] < audit["mean_time_event0"]):
        raise SystemExit(f"event coding failed for {expr_name}")
    y = pfs_binary(time, event, cut=6)
    note = cohort_note + (
        f" PFS time ranges {audit['time_min']} to {audit['time_max']}. "
        "Benefit = time >= 6 (deposited units). Event 1 has the shorter times. "
        "Early censoring is excluded. This is not RECIST ORR."
    )
    return expr, y, note, audit


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    ftp = "https://ftp.ncbi.nlm.nih.gov/geo/series"
    files = {
        "GSE126044_counts.txt.gz": f"{ftp}/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
        "GSE126044_series_matrix.txt.gz": f"{ftp}/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz",
        "GSE166449_Raw_gene_TPM_matrix.txt.gz": f"{ftp}/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz",
        "GSE166449_series_matrix.txt.gz": f"{ftp}/GSE166nnn/GSE166449/matrix/GSE166449_series_matrix.txt.gz",
        "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz": f"{ftp}/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
        "GSE135222_series_matrix.txt.gz": f"{ftp}/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz",
        "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz": f"{ftp}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
        "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx": f"{ftp}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
        "GSE190265_TPM_France3.csv.gz": f"{ftp}/GSE190nnn/GSE190265/suppl/GSE190265_TPM_France3.csv.gz",
        "GSE190265_samples_info_France3.csv.gz": f"{ftp}/GSE190nnn/GSE190265/suppl/GSE190265_samples_info_France3.csv.gz",
        "GSE190266_TPM_France4.csv.gz": f"{ftp}/GSE190nnn/GSE190266/suppl/GSE190266_TPM_France4.csv.gz",
        "GSE190266_series_matrix.txt.gz": f"{ftp}/GSE190nnn/GSE190266/matrix/GSE190266_series_matrix.txt.gz",
    }
    for name, url in files.items():
        src = Path("/tmp/geo/dl") / {
            "GSE126044_counts.txt.gz": "GSE126044_exp",
            "GSE166449_Raw_gene_TPM_matrix.txt.gz": "GSE166449_exp",
            "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz": "GSE135222_exp",
            "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz": "GSE207422_exp",
            "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx": "GSE207422_meta",
            "GSE190265_TPM_France3.csv.gz": "GSE190265_exp",
            "GSE190265_samples_info_France3.csv.gz": "GSE190265_info",
            "GSE190266_TPM_France4.csv.gz": "GSE190266_exp",
        }.get(name, name)
        dest = DATA / name
        if src.exists() and not dest.exists():
            dest.write_bytes(src.read_bytes())
        else:
            fetch(url, dest)

    bag: list[dict] = []
    audits = {}

    expr, y, note = load_gse126044()
    add_gene_rows(bag, "GSE126044", "deposited_responder", expr, y, note)

    expr, y, note = load_gse166449()
    add_gene_rows(bag, "GSE166449", "deposited_responder", expr, y, note)

    expr, y, note, audits["GSE135222"] = load_gse135222()
    add_gene_rows(bag, "GSE135222", "PFS_ge_180d", expr, y, note)

    for endpoint, expr, y, note in load_gse207422():
        add_gene_rows(bag, "GSE207422", endpoint, expr, y, note)

    info = pd.read_csv(DATA / "GSE190265_samples_info_France3.csv.gz", sep=";", compression="gzip")
    expr, y, note, audits["GSE190265"] = load_france(
        "GSE190265_TPM_France3.csv.gz",
        info,
        "sample",
        "time_PFS",
        "evtPFS",
        "Dijon France3, anti-PD-1 monotherapy, tumor at diagnosis (GSE190265).",
    )
    add_gene_rows(bag, "GSE190265", "PFS_ge_6", expr, y, note)

    ph = series_table(DATA / "GSE190266_series_matrix.txt.gz")
    # The deposited names are "pfs_time (6 months)" and "pfs_evt (6 months)".
    time_col = [c for c in ph.columns if c.startswith("pfs_time")][0]
    evt_col = [c for c in ph.columns if c.startswith("pfs_evt")][0]
    clin = ph.rename(columns={"title": "sample", time_col: "time_PFS", evt_col: "evtPFS"})
    expr, y, note, audits["GSE190266"] = load_france(
        "GSE190266_TPM_France4.csv.gz",
        clin,
        "sample",
        "time_PFS",
        "evtPFS",
        "Dijon France4, anti-PD-1 monotherapy, lung tumor (GSE190266). Time field name is the GEO string "
        + time_col
        + ".",
    )
    add_gene_rows(bag, "GSE190266", "PFS_ge_6", expr, y, note)

    df = pd.DataFrame(bag)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "sweep.tsv", sep="\t", index=False)

    scored = df[df["rule"].isin(["median", "quantile", "tertile_ends"])].copy()
    stable = scored[(scored["n_high"] >= 5) & (scored["n_low"] >= 5) & (~scored["haldane"].astype(bool))].copy()
    median = scored[(scored["rule"] == "median") & (scored["adjust"] == "none")]

    def _subset(gene: str, adjust: str | None) -> pd.DataFrame:
        out = stable[stable["gene"] == gene]
        if adjust is not None:
            out = out[out["adjust"] == adjust]
        return out

    summary = {
        "n_scored_rows": int(len(scored)),
        "n_stable_rows": int(len(stable)),
        "audits": audits,
        "median_unadjusted": median.to_dict(orient="records"),
        "cldn4_unadjusted_closest_stable_OR_lt_1": _pick(_subset("CLDN4", "none"), "closest"),
        "cldn4_unadjusted_strongest_stable_OR_lt_1": _pick(_subset("CLDN4", "none"), "strongest"),
        "cldn4_any_adjust_closest_stable_OR_lt_1": _pick(_subset("CLDN4", None), "closest"),
        "cldn4_any_adjust_strongest_stable_OR_lt_1": _pick(_subset("CLDN4", None), "strongest"),
        "tacstd2_unadjusted_closest_stable_OR_lt_1": _pick(_subset("TACSTD2", "none"), "closest"),
        "tacstd2_unadjusted_strongest_stable_OR_lt_1": _pick(_subset("TACSTD2", "none"), "strongest"),
        "tacstd2_any_adjust_strongest_stable_OR_lt_1": _pick(_subset("TACSTD2", None), "strongest"),
        "selection_rule": (
            "The gene of interest for the slide comparison is CLDN4. "
            "Stable means n_high>=5, n_low>=5, and no Haldane correction. "
            "Closest minimizes |log OR - log 0.42| among stable rows with OR<1, and ties are all returned. "
            "Strongest minimizes OR among the same rows, and ties are all returned. "
            "Median splits are pre-specified. Other quantiles and tertile ends are a search and were not multiplicity-adjusted. "
            "OR>1 rows are kept in the table; they are not dropped to match the slide."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=_json) + "\n")
    print(json.dumps({
        "cldn4_closest": summary["cldn4_unadjusted_closest_stable_OR_lt_1"],
        "cldn4_strongest": summary["cldn4_unadjusted_strongest_stable_OR_lt_1"],
    }, indent=2, default=_json))
    print("median rows", len(median), "stable", len(stable))


def _pick(frame: pd.DataFrame, how: str) -> list:
    if frame is None or len(frame) == 0:
        return []
    work = frame.copy()
    work = work[work["OR"] < 1]
    if len(work) == 0:
        return []
    if how == "closest":
        work["dist_042"] = (np.log(work["OR"]) - math.log(0.42)).abs()
        work = work.sort_values(["dist_042", "fisher_p"])
        best = float(work["dist_042"].iloc[0])
        tied = work[np.isclose(work["dist_042"].astype(float), best)]
    elif how == "strongest":
        work = work.sort_values(["OR", "fisher_p"])
        best = float(work["OR"].iloc[0])
        tied = work[np.isclose(work["OR"].astype(float), best)]
    else:
        raise SystemExit(how)
    return [_clean(row) for _, row in tied.iterrows()]


def _clean(row):
    if row is None:
        return None
    keep = [
        "cohort",
        "endpoint",
        "gene",
        "adjust",
        "rule",
        "q",
        "n",
        "n_high",
        "n_low",
        "resp_high",
        "nonresp_high",
        "resp_low",
        "nonresp_low",
        "OR",
        "OR_lo",
        "OR_hi",
        "haldane",
        "fisher_p",
        "note",
    ]
    return {k: row.get(k) for k in keep}


def _json(obj):
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    raise TypeError(type(obj))


if __name__ == "__main__":
    main()
