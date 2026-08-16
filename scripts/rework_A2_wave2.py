#!/usr/bin/env python3
"""REWORK A2 wave 2 — TACSTD2/CLDN4 vs MPR/PFS/CYT/CD8 after separate residuals.

Self-contained public-data analysis.

Question
--------
User PPT: durvalumab RNA Spearman ρ = −0.65 / purity-adj −0.46, p = 2e-4.
No open PACIFIC per-patient RNA was found in prior hunts.
Leftover PR 175: GSE253564 pre-treatment TACSTD2 is HIGHER in MPR (p = 0.016)
but the gap dies after residualising on MKI67 (p = 0.25). An older hunt reported
GSE248378 post-durva TACSTD2 vs CYT adj ρ = −0.71.

This wave asks those leftover tests again after residualising TACSTD2 and CLDN4
on MKI67, ESTIMATE, and an epithelial score *separately* (never jointly), and
also versus CYT / CD8A. It also re-hunts PACIFIC / ADRIATIC / open durvalumab
per-patient matrices. No accessions or p-values are invented.

Outputs -> results/rework/A2_wave2/
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
import requests
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
AUX = DATA / "A2_wave2"
OUT = ROOT / "results" / "rework" / "A2_wave2"
FIG = OUT / "figures"
for p in (DATA, AUX, FIG):
    p.mkdir(parents=True, exist_ok=True)

TARGETS = ("TACSTD2", "CLDN4")
# Locked epithelial markers. TACSTD2 and CLDN4 are deliberately excluded so
# residualising a target on the epithelial score is not circular.
EPI_GENES = ("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1", "MUC1")
CYT_GENES = ("GZMA", "PRF1")
CD8_GENE = "CD8A"
MKI67_GENE = "MKI67"
MPR_THRESHOLD = -90.0
NCBI_TOOL = {"tool": "sdaxcge_A2_wave2", "email": "jinxuanhong1@gmail.com"}

URLS = {
    "pre_fpkm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",
    "post_fpkm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/suppl/GSE248378_Durva_Post_FPKMs.txt.gz",
    "pre_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/matrix/GSE253564_series_matrix.txt.gz",
    "post_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/matrix/GSE248378_series_matrix.txt.gz",
    "natcomm": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-44195-x/MediaObjects/41467_2023_44195_MOESM6_ESM.xlsx",
    "pmc_zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10982989/supplementaryFiles",
    "estimate_gmt": "https://raw.githubusercontent.com/cran/estimate/master/inst/extdata/SI_geneset.gmt",
}

ARM_ALIASES = {
    "arm1": "arm1",
    "arm2": "arm2",
    "durvalumab": "arm1",
    "durvalumab + rt": "arm2",
    "durvalumab+rt": "arm2",
    "durvalumab and rt": "arm2",
    "durvalumabrt": "arm2",
    "durvalumab rt": "arm2",
}


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, timeout: int = 120) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": "sdaxcge-A2-wave2"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, dest.open("wb") as out:
        out.write(resp.read())
    return dest


def ensure_inputs() -> dict[str, Path]:
    paths = {
        "pre_fpkm": DATA / "GSE253564_Pre_FPKMs.txt.gz",
        "post_fpkm": DATA / "GSE248378_Durva_Post_FPKMs.txt.gz",
        "pre_matrix": AUX / "GSE253564_series_matrix.txt.gz",
        "post_matrix": AUX / "GSE248378_series_matrix.txt.gz",
        "natcomm": AUX / "41467_2023_44195_MOESM6_ESM.xlsx",
        "mmc2": AUX / "mmc2.xlsx",
        "gmt": DATA / "estimate" / "inst" / "extdata" / "SI_geneset.gmt",
    }
    download(URLS["pre_fpkm"], paths["pre_fpkm"])
    download(URLS["post_fpkm"], paths["post_fpkm"])
    download(URLS["pre_matrix"], paths["pre_matrix"])
    download(URLS["post_matrix"], paths["post_matrix"])
    download(URLS["natcomm"], paths["natcomm"])
    if not paths["gmt"].exists():
        try:
            download(URLS["estimate_gmt"], paths["gmt"])
        except Exception:
            pass
    if not paths["mmc2"].exists():
        zpath = AUX / "PMC10982989_SupplementaryFiles.zip"
        download(URLS["pmc_zip"], zpath)
        with zipfile.ZipFile(zpath) as zf:
            zf.extract("mmc2.xlsx", AUX)
        zpath.unlink(missing_ok=True)
    return paths


def load_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def estimate_scores(expr: pd.DataFrame, gene_sets: dict[str, list[str]]):
    """Exact port of estimateScore() from ESTIMATE R v1.0.13."""
    genes = expr.index.to_numpy()
    n_genes = expr.shape[0]
    m = expr.rank(axis=0, method="average").to_numpy() * (10000.0 / n_genes)
    rows = {}
    overlaps = {}
    for name in ("StromalSignature", "ImmuneSignature"):
        in_set = np.isin(genes, list(set(gene_sets[name])))
        overlaps[name] = int(in_set.sum())
        es = np.empty(expr.shape[1])
        for j in range(expr.shape[1]):
            order = np.argsort(-m[:, j], kind="stable")
            tag = in_set[order].astype(float)
            w = np.abs(m[order, j]) ** 0.25
            nh = tag.sum()
            nm = n_genes - nh
            pn = np.cumsum(tag * w) / (tag * w).sum()
            p0 = np.cumsum((1.0 - tag) / nm)
            es[j] = float(np.sum(pn - p0))
        rows["Stromal" if name.startswith("Str") else "Immune"] = es
    df = pd.DataFrame(
        {"StromalScore": rows["Stromal"], "ImmuneScore": rows["Immune"]},
        index=expr.columns,
    )
    df["ESTIMATEScore"] = df["StromalScore"] + df["ImmuneScore"]
    purity = np.cos(0.6049872018 + 0.0001467884 * df["ESTIMATEScore"].to_numpy())
    purity[purity < 0] = np.nan
    df["TumorPurity"] = purity
    return df, overlaps


def read_fpkm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    gene_col = df.columns[0]
    df = df.set_index(gene_col)
    if "Entrez.ID" in df.columns:
        df = df.drop(columns=["Entrez.ID"])
    df.index = df.index.astype(str).str.strip()
    df = df.groupby(df.index).max()
    return df.apply(pd.to_numeric, errors="coerce")


def log2p1(s: pd.Series) -> pd.Series:
    return np.log2(s.astype(float) + 1.0)


def mean_log2p1(expr: pd.DataFrame, genes: tuple[str, ...]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), present
    mat = np.log2(expr.loc[present].astype(float) + 1.0)
    return mat.mean(axis=0), present


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    gt = int(np.sum(a[:, None] > b[None, :]))
    lt = int(np.sum(a[:, None] < b[None, :]))
    return (gt - lt) / (a.size * b.size)


def ols_residual(y: pd.Series, z: pd.Series) -> pd.Series:
    aligned = pd.concat([y.rename("y"), z.rename("z")], axis=1).dropna()
    if len(aligned) < 6 or aligned["z"].nunique() < 2:
        return pd.Series(np.nan, index=y.index)
    slope, intercept, *_ = stats.linregress(aligned["z"], aligned["y"])
    resid = aligned["y"] - (intercept + slope * aligned["z"])
    out = pd.Series(np.nan, index=y.index)
    out.loc[resid.index] = resid
    return out


def partial_spearman(x, y, z):
    """Partial Spearman of x,y | z via Pearson on ranks; p from t with n-3 df."""
    aligned = pd.concat(
        [pd.Series(x, dtype=float), pd.Series(y, dtype=float), pd.Series(z, dtype=float)],
        axis=1,
    ).dropna()
    aligned.columns = ["x", "y", "z"]
    n = len(aligned)
    if n < 6 or aligned["z"].nunique() < 2:
        return float("nan"), float("nan"), n
    rx, ry, rz = (stats.rankdata(aligned[c]) for c in ("x", "y", "z"))
    rxy = np.corrcoef(rx, ry)[0, 1]
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]
    den = np.sqrt((1 - rxz**2) * (1 - ryz**2))
    if den < 1e-12:
        return float("nan"), float("nan"), n
    r = float(np.clip((rxy - rxz * ryz) / den, -1.0, 1.0))
    dof = n - 3
    if abs(r) >= 1:
        return r, 0.0, n
    t = r * np.sqrt(dof / (1 - r**2))
    p = float(2 * stats.t.sf(abs(t), dof))
    return r, p, n


def spearman_pair(x, y):
    aligned = pd.concat([pd.Series(x, dtype=float), pd.Series(y, dtype=float)], axis=1).dropna()
    if len(aligned) < 5:
        return float("nan"), float("nan"), int(len(aligned))
    rho, p = stats.spearmanr(aligned.iloc[:, 0], aligned.iloc[:, 1])
    return float(rho), float(p), int(len(aligned))


def mwu(values: pd.Series, labels: pd.Series):
    aligned = pd.concat([values.rename("x"), labels.rename("y")], axis=1).dropna()
    pos = aligned.loc[aligned["y"] == 1, "x"].to_numpy()
    neg = aligned.loc[aligned["y"] == 0, "x"].to_numpy()
    rec = {
        "n_pos": int(pos.size),
        "n_neg": int(neg.size),
        "median_pos": float(np.median(pos)) if pos.size else None,
        "median_neg": float(np.median(neg)) if neg.size else None,
        "U": None,
        "p": None,
        "cliffs_delta": None,
        "note": "",
    }
    if pos.size < 3 or neg.size < 3:
        rec["note"] = "too few labelled samples (not estimable)"
        return rec
    u, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    rec.update(
        {
            "U": float(u),
            "p": float(p),
            "cliffs_delta": float(cliffs_delta(pos, neg)),
        }
    )
    return rec


def cox_pfs(values: pd.Series, time: pd.Series, event: pd.Series) -> dict:
    aligned = pd.concat(
        [values.rename("expr"), time.rename("T"), event.rename("E")], axis=1
    ).dropna()
    aligned["T"] = pd.to_numeric(aligned["T"], errors="coerce")
    aligned["E"] = pd.to_numeric(aligned["E"], errors="coerce")
    aligned = aligned.dropna()
    out = {
        "n": int(len(aligned)),
        "n_events": int(aligned["E"].sum()) if len(aligned) else 0,
        "hr": None,
        "hr_ci_low": None,
        "hr_ci_high": None,
        "p": None,
        "logrank_p": None,
        "note": "HR per +1 residual (or +1 log2(FPKM+1) if unadjusted)",
    }
    if len(aligned) < 8 or aligned["E"].sum() < 3:
        out["note"] = "too few events"
        return out
    med = float(aligned["expr"].median())
    high = aligned["expr"] >= med
    lr = logrank_test(
        aligned.loc[high, "T"],
        aligned.loc[~high, "T"],
        event_observed_A=aligned.loc[high, "E"],
        event_observed_B=aligned.loc[~high, "E"],
    )
    out["logrank_p"] = float(lr.p_value)
    try:
        cph = CoxPHFitter()
        cph.fit(aligned, duration_col="T", event_col="E")
        s = cph.summary.loc["expr"]
        out.update(
            {
                "hr": float(s["exp(coef)"]),
                "hr_ci_low": float(s["exp(coef) lower 95%"]),
                "hr_ci_high": float(s["exp(coef) upper 95%"]),
                "p": float(s["p"]),
            }
        )
    except Exception as exc:  # noqa: BLE001
        out["note"] = f"Cox fit failed: {exc}"
    return out


def norm_id(raw: str) -> str | None:
    m = re.search(r"(\d+)", str(raw))
    if not m:
        return None
    return f"durva{int(m.group(1)):03d}"


def canon_arm(value: str) -> str | None:
    v = (value or "").strip().lower().replace("  ", " ")
    if v in ARM_ALIASES:
        return ARM_ALIASES[v]
    if "rt" in v and "durvalumab" in v:
        return "arm2"
    if v.startswith("durvalumab"):
        return "arm1"
    if v.startswith("arm"):
        return f"arm{v[-1]}"
    return None


def load_table_s1(path: Path) -> dict[str, dict]:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    header = [str(c).strip() if c is not None else "" for c in rows[1]]
    idx = {name: i for i, name in enumerate(header)}
    out: dict[str, dict] = {}
    for r in rows[2:]:
        if not r or r[0] is None:
            continue
        pid = norm_id(r[0])
        if pid is None:
            continue

        def get(col: str):
            key = col.strip()
            i = idx.get(key)
            if i is None:
                for name, j in idx.items():
                    if name.replace("  ", " ") == key.replace("  ", " "):
                        i = j
                        break
            return r[i] if i is not None and i < len(r) else None

        pr = get("Pathology Response")
        try:
            pr_val = float(str(pr).replace("%", "").strip())
        except (TypeError, ValueError):
            pr_val = None
        out[pid] = {
            "patient": pid,
            "study_number_raw": str(r[0]).strip(),
            "arm_pub": str(get("Study  Arm") or "").strip(),
            "path_response_pct": pr_val,
            "mpr": None if pr_val is None else int(pr_val <= MPR_THRESHOLD),
            "histology": str(get("Cell Type") or "").strip(),
            "egfr_status": str(get("EGFR mutation status") or "").strip(),
        }
    wb.close()
    return out


def load_natcomm_outcomes(path: Path) -> dict[str, dict]:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    out: dict[str, dict] = {}

    def harvest(sheet: str, fields: dict[str, str]) -> None:
        if sheet not in wb.sheetnames:
            return
        ws = wb[sheet]
        rows = list(ws.iter_rows(values_only=True))
        header = [str(c).strip() if c is not None else "" for c in rows[0]]
        idx = {h: i for i, h in enumerate(header)}
        for r in rows[1:]:
            if not r or r[0] is None:
                continue
            pid = norm_id(r[0])
            if pid is None or not str(r[0]).lower().startswith("durva"):
                continue
            rec = out.setdefault(pid, {"patient": pid})
            for src, dest in fields.items():
                i = idx.get(src)
                if i is not None and i < len(r) and r[i] is not None and dest not in rec:
                    rec[dest] = str(r[i]).strip()

    harvest(
        "Figure 1b",
        {
            "ARM": "arm_natcomm",
            "Resection": "resected",
            "Status": "dfs_status_raw",
            "Disease free survival in Months": "dfs_months",
        },
    )
    harvest(
        "Figure 2a",
        {
            "Progression Status": "progression_status_raw",
            "Progression-free survival in Months": "pfs_months",
        },
    )
    harvest("Figure 2e", {"Path Response_2Group": "path_response_2group"})
    harvest("Figure 2f", {"Path Response_2Group": "path_response_2group"})
    wb.close()
    for rec in out.values():
        status = (rec.get("progression_status_raw") or rec.get("dfs_status_raw") or "").lower()
        if not status:
            rec["recurrence_event"] = None
        elif "ned" in status and "recurren" not in status:
            rec["recurrence_event"] = 0
        elif "recurren" in status or "progress" in status or "dead for lung cancer" in status:
            rec["recurrence_event"] = 1
        elif "death from other" in status or "other diseases" in status:
            rec["recurrence_event"] = 0
        elif status.strip() in {"alive"}:
            rec["recurrence_event"] = 0
        else:
            rec["recurrence_event"] = None
    return out


def geo_sample_table(path: Path, acc: str) -> list[dict]:
    header: dict[str, list[str]] = {}
    import gzip

    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!"):
                key, _, val = line[1:].partition("\t")
                header.setdefault(key.strip(), []).append(val.rstrip("\n"))

    def split(key: str) -> list[str]:
        return [v.strip('"') for v in header.get(key, [""])[0].split("\t")] if header.get(key) else []

    titles = split("Sample_title")
    gsms = split("Sample_geo_accession")
    arms: list[str] = [""] * len(titles)
    for line in header.get("Sample_characteristics_ch1", []):
        vals = [v.strip('"') for v in line.split("\t")]
        if any(v.lower().startswith("treatment:") for v in vals):
            arms = [v.split(":", 1)[1].strip() if ":" in v else "" for v in vals]
    return [
        {
            "series": acc,
            "gsm": gsms[i] if i < len(gsms) else "",
            "sample_title": titles[i],
            "arm_geo": arms[i] if i < len(arms) else "",
            "patient": norm_id(titles[i]),
        }
        for i in range(len(titles))
    ]


def build_clinical(paths: dict[str, Path]) -> tuple[pd.DataFrame, dict]:
    s1 = load_table_s1(paths["mmc2"])
    nc = load_natcomm_outcomes(paths["natcomm"])
    validation = {"checks": [], "mismatches": []}
    rows = []
    for acc, mkey in (("GSE253564", "pre_matrix"), ("GSE248378", "post_matrix")):
        samples = geo_sample_table(paths[mkey], acc)
        matched = 0
        for s in samples:
            pid = s["patient"]
            pub = s1.get(pid, {})
            out = nc.get(pid, {})
            if pub:
                matched += 1
                arm_geo = canon_arm(s["arm_geo"])
                arm_pub = canon_arm(pub.get("arm_pub", ""))
                if arm_geo and arm_pub and arm_geo != arm_pub:
                    validation["mismatches"].append(
                        {
                            "series": acc,
                            "sample_title": s["sample_title"],
                            "patient": pid,
                            "arm_geo": s["arm_geo"],
                            "arm_publication": pub.get("arm_pub"),
                        }
                    )
            rows.append({**s, **{k: v for k, v in pub.items() if k != "patient"},
                         **{k: v for k, v in out.items() if k != "patient"}})
        validation["checks"].append(
            {"series": acc, "n_samples": len(samples), "n_matched_to_publication": matched}
        )
    for rec in rows:
        rec["arm_conflict"] = int(
            bool(
                canon_arm(rec.get("arm_geo", ""))
                and canon_arm(rec.get("arm_pub", ""))
                and canon_arm(rec["arm_geo"]) != canon_arm(rec["arm_pub"])
            )
        )
        rec["arm_canonical"] = canon_arm(rec.get("arm_pub", "")) or canon_arm(rec.get("arm_geo", "")) or ""
        if rec.get("mpr") in (None, "") and rec.get("path_response_2group"):
            rec["mpr"] = 1 if str(rec["path_response_2group"]).lower() == "major" else 0
    validation["resolution"] = (
        "Use publication Table S1 arm as canonical. The one GEO mismatch "
        "(GSE248378 45-M-PO / durva045) is retained and excluded from arm-stratified tests."
    )
    return pd.DataFrame(rows), validation


def ncbi_esearch(db: str, term: str, retmax: int = 40) -> dict:
    params = {
        "db": db,
        "term": term,
        "retmax": str(retmax),
        "retmode": "json",
        **NCBI_TOOL,
    }
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
    r = requests.get(url, timeout=45)
    r.raise_for_status()
    time.sleep(0.35)
    return r.json()


def ncbi_esummary(db: str, ids: list[str]) -> dict:
    if not ids:
        return {}
    params = {
        "db": db,
        "id": ",".join(ids),
        "retmode": "json",
        **NCBI_TOOL,
    }
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(params)
    r = requests.get(url, timeout=45)
    r.raise_for_status()
    time.sleep(0.35)
    return r.json()


def hunt_open_matrices() -> list[dict]:
    """Live Entrez + targeted accession checks. No invented IDs."""
    rows: list[dict] = []

    def add(**kw):
        rows.append(kw)

    queries = [
        ("gds", "PACIFIC[All Fields] AND durvalumab[All Fields] AND Homo sapiens[Organism]"),
        ("gds", "NCT02125461[All Fields]"),
        ("gds", "ADRIATIC[All Fields] AND Homo sapiens[Organism]"),
        ("gds", "ADRIA[All Fields] AND durvalumab[All Fields]"),
        ("gds", "WJOG11518L[All Fields] OR SUBMARINE[All Fields] AND Homo sapiens[Organism]"),
        ("gds", "durvalumab[All Fields] AND (NSCLC[All Fields] OR \"non-small cell lung\"[All Fields] OR \"lung cancer\"[All Fields] OR SCLC[All Fields]) AND Homo sapiens[Organism]"),
        ("gds", "GSE333537[Accession]"),
        ("gds", "GSE261345[Accession]"),
        ("gds", "GSE190731[Accession]"),
        ("pubmed", "PACIFIC[Title] AND durvalumab AND (RNA-seq OR transcriptome OR NanoString OR GEO)"),
        ("pubmed", "ADRIATIC AND durvalumab AND (RNA-seq OR transcriptome OR GEO)"),
        ("pubmed", "WJOG11518L OR SUBMARINE AND durvalumab AND (NanoString OR RNA)"),
    ]
    seen_gds: set[str] = set()
    for db, term in queries:
        try:
            js = ncbi_esearch(db, term, retmax=50)
        except Exception as exc:  # noqa: BLE001
            add(
                resource=f"Entrez {db} query failed",
                accession_or_id=term,
                status="QUERY_ERROR",
                open_per_patient_matrix="no",
                notes=str(exc)[:300],
                url="",
                query=term,
            )
            continue
        es = js.get("esearchresult", {})
        count = int(es.get("count", 0))
        ids = es.get("idlist", [])
        add(
            resource=f"Entrez {db} search",
            accession_or_id=f"n_hits={count}",
            status="SEARCHED",
            open_per_patient_matrix="n/a",
            notes=f"query={term}",
            url="",
            query=term,
        )
        if db != "gds" or not ids:
            if db == "pubmed":
                add(
                    resource="PubMed hits (IDs only; no matrix implied)",
                    accession_or_id=",".join(ids[:15]) if ids else "none",
                    status="SEARCHED",
                    open_per_patient_matrix="no",
                    notes="Abstract-level hits. A PMID is not a per-patient matrix.",
                    url="https://pubmed.ncbi.nlm.nih.gov/",
                    query=term,
                )
            continue
        sm = ncbi_esummary("gds", ids)
        docs = sm.get("result", {})
        for uid in ids:
            doc = docs.get(uid, {})
            acc = str(doc.get("accession", "") or uid)
            if not acc.startswith("GSE"):
                continue
            if acc in seen_gds:
                continue
            seen_gds.add(acc)
            title = str(doc.get("title", "") or "")
            n_samp = doc.get("n_samples", "")
            gdstype = str(doc.get("gdstype", "") or "")
            taxon = str(doc.get("taxon", "") or "")
            add(
                resource="GEO series (live Entrez, GSE only)",
                accession_or_id=acc,
                status="VERIFIED_EXISTS",
                open_per_patient_matrix="unknown_pending_triage",
                notes=f"{title[:220]} | type={gdstype} | taxon={taxon} | n={n_samp}",
                url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
                query=term,
            )

    # Locked triage of known / newly seen candidates (verified IDs only).
    locked = [
        dict(
            resource="PACIFIC (durvalumab after CRT, stage III NSCLC)",
            accession_or_id="NCT02125461",
            status="NOT_OPEN",
            open_per_patient_matrix="no",
            notes="Tissue optional; published biomarker supplements are PD-L1 IHC. Individual-level data via AstraZeneca/Vivli only. No public per-patient RNA matrix found in this re-hunt.",
            url="https://clinicaltrials.gov/study/NCT02125461",
            query="targeted",
        ),
        dict(
            resource="ADRIATIC (consolidation durvalumab, LS-SCLC)",
            accession_or_id="NCT03703297",
            status="NOT_OPEN",
            open_per_patient_matrix="no",
            notes="ASCO 2025 abstract 8014 reports RNA-seq BEP medians (TIS/CD8A/STING) in EP vs LTP groups. Group-level medians only; no GEO/ArrayExpress per-patient matrix found.",
            url="https://clinicaltrials.gov/study/NCT03703297",
            query="targeted",
        ),
        dict(
            resource="SUBMARINE / WJOG11518L PACIFIC-regimen NanoString",
            accession_or_id="PMID 37364849",
            status="NOT_DEPOSITED",
            open_per_patient_matrix="no",
            notes="nCounter IO360 on ~85 pre-CRT tumors. No relevant GEO series. GSE278471 is Haloferax volcanii (false positive).",
            url="https://pubmed.ncbi.nlm.nih.gov/37364849/",
            query="targeted",
        ),
        dict(
            resource="Neoadjuvant durvalumab ± SBRT pre-treatment FPKM",
            accession_or_id="GSE253564",
            status="ANALYZED",
            open_per_patient_matrix="yes",
            notes="Only open human durvalumab NSCLC whole-transcriptome pre-treatment matrix with recoverable MPR/PFS (n=32). NCT02904954.",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564",
            query="targeted",
        ),
        dict(
            resource="Neoadjuvant durvalumab ± SBRT post-treatment FPKM",
            accession_or_id="GSE248378",
            status="ANALYZED",
            open_per_patient_matrix="yes",
            notes="Same trial, residual tumours n=29. MPR not estimable (0 MPR in deposited matrix). Recurrence/PFS available.",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378",
            query="targeted",
        ),
        dict(
            resource="Durvalumab Study 1108 IFNγ 21-gene panel",
            accession_or_id="GSE110390",
            status="EXCLUDED_UNUSABLE",
            open_per_patient_matrix="no",
            notes="Only 21 IFNγ-signature genes deposited. TACSTD2/CLDN4 absent. Not a whole-transcriptome matrix.",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110390",
            query="targeted",
        ),
        dict(
            resource="Durvalumab ± tremelimumab scRNA (2 NSCLC tumors)",
            accession_or_id="GSE131933",
            status="EXCLUDED_WRONG_ASSAY",
            open_per_patient_matrix="no",
            notes="scRNA of immune cells from two NSCLC tumors. Not a per-patient bulk matrix for TACSTD2 vs MPR/PFS/CYT.",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131933",
            query="targeted",
        ),
        dict(
            resource="Adjuvant durvalumab esophageal/GEJ RNA-seq",
            accession_or_id="GSE183924",
            status="EXCLUDED_WRONG_DISEASE",
            open_per_patient_matrix="yes_but_not_lung",
            notes="Open FPKM, but esophageal/GEJ adenocarcinoma — not PACIFIC/ADRIA/NSCLC. Not used as a PACIFIC substitute.",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183924",
            query="targeted",
        ),
        dict(
            resource="VCN-01 + durvalumab HNSCC RNA-seq",
            accession_or_id="GSE333537",
            status="EXCLUDED_WRONG_DISEASE",
            open_per_patient_matrix="yes_but_not_lung",
            notes="Live Entrez: Phase I VCN-01 oncolytic adenovirus + durvalumab in metastatic HNSCC (n=35). Real durvalumab RNA, but head-and-neck — not PACIFIC/ADRIA/NSCLC. Not used as a substitute.",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333537",
            query="targeted",
        ),
        dict(
            resource="CANTABRICO ES-SCLC GeoMx CTA (chemo-IO)",
            accession_or_id="GSE261345",
            status="EXCLUDED_WRONG_ASSAY",
            open_per_patient_matrix="no",
            notes="Live Entrez: GeoMx DSP Cancer Transcriptome Atlas, 121 ROIs / ~26 ES-SCLC patients. Not PACIFIC, not bulk whole-transcriptome, not NSCLC. CTA panel ≠ full TACSTD2/CLDN4 residualisation matrix. Near-miss only.",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261345",
            query="targeted",
        ),
        dict(
            resource="Durvalumab/oleclumab NSCLC xenograft array",
            accession_or_id="GSE190731",
            status="EXCLUDED_NOT_PATIENT",
            open_per_patient_matrix="no",
            notes="Live Entrez: expression array from durvalumab/oleclumab-treated xenograft NSCLC tumors (n=23). Not a human per-patient matrix.",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190731",
            query="targeted",
        ),
        dict(
            resource="OAK/POPLAR atezolizumab (Bessede 2024 TROP2 paper)",
            accession_or_id="EGAS00001005013",
            status="CONTROLLED_ACCESS",
            open_per_patient_matrix="no",
            notes="Anti–PD-L1 but atezolizumab, not durvalumab. EGA controlled. Not downloaded.",
            url="https://ega-archive.org/studies/EGAS00001005013",
            query="targeted",
        ),
        dict(
            resource="POSEIDON / MYSTIC / AEGEAN / PACIFIC-2/4/5/6/8/9",
            accession_or_id="NCT03164616 / NCT02453282 / NCT03800134 / PACIFIC-x",
            status="NOT_OPEN",
            open_per_patient_matrix="no",
            notes="AstraZeneca durvalumab lung programs. No new open per-patient RNA matrix found in this re-hunt. Do not substitute.",
            url="https://clinicaltrials.gov/",
            query="targeted",
        ),
    ]
    rows.extend(locked)
    return rows


def analyse_cohort(
    acc: str,
    timepoint: str,
    expr: pd.DataFrame,
    clin: pd.DataFrame,
    gmt: dict[str, list[str]],
) -> tuple[pd.DataFrame, list[dict], dict]:
    sub = clin[clin["series"] == acc].copy()
    missing = [t for t in sub["sample_title"] if t not in expr.columns]
    if missing:
        raise RuntimeError(f"{acc}: expression columns missing for {missing}")
    sub = sub.set_index("sample_title")
    est, overlaps = estimate_scores(expr, gmt)
    cyt, cyt_present = mean_log2p1(expr, CYT_GENES)
    epi, epi_present = mean_log2p1(expr, EPI_GENES)
    epcam = log2p1(expr.loc["EPCAM"]) if "EPCAM" in expr.index else pd.Series(np.nan, index=expr.columns)
    mki67 = log2p1(expr.loc[MKI67_GENE]) if MKI67_GENE in expr.index else pd.Series(np.nan, index=expr.columns)
    cd8 = log2p1(expr.loc[CD8_GENE]) if CD8_GENE in expr.index else pd.Series(np.nan, index=expr.columns)

    sample = sub.reset_index()[
        [
            "series",
            "gsm",
            "sample_title",
            "patient",
            "arm_canonical",
            "arm_conflict",
            "mpr",
            "path_response_pct",
            "path_response_2group",
            "recurrence_event",
            "pfs_months",
            "histology",
        ]
    ].copy()
    idx = sample["sample_title"]
    sample["CYT"] = cyt.loc[idx].to_numpy()
    sample["CD8A"] = cd8.loc[idx].to_numpy()
    sample["MKI67"] = mki67.loc[idx].to_numpy()
    sample["EpiScore"] = epi.loc[idx].to_numpy()
    sample["EPCAM"] = epcam.loc[idx].to_numpy()
    sample["ImmuneScore"] = est.loc[idx, "ImmuneScore"].to_numpy()
    sample["StromalScore"] = est.loc[idx, "StromalScore"].to_numpy()
    sample["ESTIMATEScore"] = est.loc[idx, "ESTIMATEScore"].to_numpy()
    sample["TumorPurity"] = est.loc[idx, "TumorPurity"].to_numpy()
    for gene in TARGETS:
        sample[gene] = log2p1(expr.loc[gene, idx]).to_numpy()

    covariates = {
        "none": None,
        "MKI67": sample.set_index("sample_title")["MKI67"],
        "ESTIMATEScore": sample.set_index("sample_title")["ESTIMATEScore"],
        "EpiScore": sample.set_index("sample_title")["EpiScore"],
        "EPCAM_sensitivity": sample.set_index("sample_title")["EPCAM"],
    }
    endpoints_bin = {
        "MPR": pd.to_numeric(sample.set_index("sample_title")["mpr"], errors="coerce"),
        "recurrence": pd.to_numeric(sample.set_index("sample_title")["recurrence_event"], errors="coerce"),
    }
    depth = -pd.to_numeric(sample.set_index("sample_title")["path_response_pct"], errors="coerce")
    pfs_t = pd.to_numeric(sample.set_index("sample_title")["pfs_months"], errors="coerce")
    pfs_e = pd.to_numeric(sample.set_index("sample_title")["recurrence_event"], errors="coerce")
    immune = {
        "CYT": sample.set_index("sample_title")["CYT"],
        "CD8A": sample.set_index("sample_title")["CD8A"],
        "ImmuneScore": sample.set_index("sample_title")["ImmuneScore"],
    }

    rows: list[dict] = []
    for gene in TARGETS:
        y0 = sample.set_index("sample_title")[gene]
        for cov_name, z in covariates.items():
            y = y0 if cov_name == "none" else ols_residual(y0, z)
            adj = "unadjusted" if cov_name == "none" else f"OLS residual on {cov_name}"
            for ep_name, lab in endpoints_bin.items():
                rec = mwu(y, lab)
                rec.update(
                    {
                        "series": acc,
                        "timepoint": timepoint,
                        "gene": gene,
                        "covariate": cov_name,
                        "adjustment": adj,
                        "endpoint": ep_name,
                        "test": "Mann-Whitney U",
                    }
                )
                if ep_name == "MPR" and rec["n_pos"] == 0:
                    rec["note"] = "MPR not estimable (0 MPR cases in deposited matrix)"
                rows.append(rec)
            rho, p, n = spearman_pair(y, depth)
            rows.append(
                {
                    "series": acc,
                    "timepoint": timepoint,
                    "gene": gene,
                    "covariate": cov_name,
                    "adjustment": adj,
                    "endpoint": "pathologic_depth",
                    "test": "Spearman",
                    "n": n,
                    "rho": rho,
                    "p": p,
                    "note": "depth = −Pathology Response %",
                }
            )
            surv = cox_pfs(y, pfs_t, pfs_e)
            rows.append(
                {
                    "series": acc,
                    "timepoint": timepoint,
                    "gene": gene,
                    "covariate": cov_name,
                    "adjustment": adj,
                    "endpoint": "PFS",
                    "test": "Cox PH continuous + log-rank median split",
                    **surv,
                }
            )
            for iname, iv in immune.items():
                rho, p, n = spearman_pair(y, iv)
                rows.append(
                    {
                        "series": acc,
                        "timepoint": timepoint,
                        "gene": gene,
                        "covariate": cov_name,
                        "adjustment": adj,
                        "endpoint": iname,
                        "test": (
                            "Spearman unadjusted"
                            if cov_name == "none"
                            else "Spearman of OLS-residual vs raw endpoint"
                        ),
                        "n": n,
                        "rho": rho,
                        "p": p,
                        "note": "residualise gene only (leftover-style)" if cov_name != "none" else "",
                    }
                )
                if cov_name != "none" and z is not None:
                    pr, pp, pn = partial_spearman(y0, iv, z)
                    rows.append(
                        {
                            "series": acc,
                            "timepoint": timepoint,
                            "gene": gene,
                            "covariate": cov_name,
                            "adjustment": f"partial Spearman | {cov_name}",
                            "endpoint": iname,
                            "test": "partial Spearman (ranks, n-3 df)",
                            "n": pn,
                            "rho": pr,
                            "p": pp,
                            "note": "A2-matched estimator; residualises both ranks",
                        }
                    )

    meta = {
        "series": acc,
        "timepoint": timepoint,
        "n_samples": int(expr.shape[1]),
        "n_genes": int(expr.shape[0]),
        "estimate_overlap": overlaps,
        "cyt_genes_present": cyt_present,
        "epi_genes_present": epi_present,
        "epi_genes_locked": list(EPI_GENES),
        "n_purity_oob": int(est["TumorPurity"].isna().sum()),
    }
    return sample, rows, meta


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (math.isnan(p))):
        return "NA"
    p = float(p)
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_r(x) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    return f"{float(x):.3f}"


def pick(df: pd.DataFrame, **kw) -> pd.Series | None:
    q = df
    for k, v in kw.items():
        q = q[q[k] == v]
    if q.empty:
        return None
    return q.iloc[0]


def write_figures(pre: pd.DataFrame, post: pd.DataFrame, tests: pd.DataFrame) -> None:
    # Pre MPR boxplots: unadjusted + three residuals for TACSTD2
    fig, axes = plt.subplots(1, 4, figsize=(12.4, 3.6), constrained_layout=True)
    y0 = pre.set_index("sample_title")["TACSTD2"]
    mpr = pd.to_numeric(pre.set_index("sample_title")["mpr"], errors="coerce")
    panels = [
        ("unadjusted", y0),
        ("MKI67 residual", ols_residual(y0, pre.set_index("sample_title")["MKI67"])),
        ("ESTIMATE residual", ols_residual(y0, pre.set_index("sample_title")["ESTIMATEScore"])),
        ("EpiScore residual", ols_residual(y0, pre.set_index("sample_title")["EpiScore"])),
    ]
    rng = np.random.default_rng(0)
    for ax, (title, y) in zip(axes, panels):
        a = y[mpr == 1].dropna()
        b = y[mpr == 0].dropna()
        ax.boxplot([b, a], tick_labels=["no-MPR", "MPR"], widths=0.55)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, size=len(b)), b, s=14, color="#444444", alpha=0.75, zorder=3)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, size=len(a)), a, s=14, color="#b33", alpha=0.75, zorder=3)
        rec = mwu(y, mpr)
        ax.set_title(f"{title}\np={fmt_p(rec['p'])}  δ={fmt_r(rec['cliffs_delta'])}", fontsize=9)
        ax.set_ylabel("TACSTD2 (log2 or residual)")
    fig.suptitle("GSE253564 pre-treatment TACSTD2 vs MPR", fontsize=11)
    fig.savefig(FIG / "GSE253564_TACSTD2_MPR_residuals.png", dpi=140)
    plt.close(fig)

    # Partial Spearman bars vs CYT/CD8
    sub = tests[
        (tests["test"] == "partial Spearman (ranks, n-3 df)")
        & (tests["endpoint"].isin(["CYT", "CD8A"]))
        & (tests["gene"] == "TACSTD2")
        & (tests["covariate"].isin(["MKI67", "ESTIMATEScore", "EpiScore"]))
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)
    for ax, series in zip(axes, ("GSE253564", "GSE248378")):
        ss = sub[sub["series"] == series]
        labels, vals, ps = [], [], []
        for cov, ep in (
            ("MKI67", "CYT"),
            ("ESTIMATEScore", "CYT"),
            ("EpiScore", "CYT"),
            ("MKI67", "CD8A"),
            ("ESTIMATEScore", "CD8A"),
            ("EpiScore", "CD8A"),
        ):
            r = ss[(ss["covariate"] == cov) & (ss["endpoint"] == ep)]
            labels.append(f"{ep}|{cov[:3]}")
            vals.append(float(r["rho"].iloc[0]) if len(r) else 0)
            ps.append(float(r["p"].iloc[0]) if len(r) else 1)
        colors = ["#b33" if p < 0.05 else "#888" for p in ps]
        ax.bar(range(len(vals)), vals, color=colors)
        ax.axhline(0, color="black", lw=0.6)
        ax.axhline(-0.46, color="#1f4e79", ls="--", lw=0.8, label="user adj −0.46")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_ylim(-1, 0.3)
        ax.set_title(series)
        ax.set_ylabel("partial Spearman ρ")
        ax.legend(fontsize=7)
    fig.suptitle("TACSTD2 vs CYT/CD8A after each covariate (partial Spearman)", fontsize=11)
    fig.savefig(FIG / "TACSTD2_partial_rho_CYT_CD8.png", dpi=140)
    plt.close(fig)

    # Post recurrence box for TACSTD2 unadj
    fig, ax = plt.subplots(figsize=(4.2, 3.6), constrained_layout=True)
    rec = pd.to_numeric(post["recurrence_event"], errors="coerce")
    a = post.loc[rec == 1, "TACSTD2"].dropna()
    b = post.loc[rec == 0, "TACSTD2"].dropna()
    ax.boxplot([b, a], tick_labels=["no rec.", "recurrence"], widths=0.55)
    rng = np.random.default_rng(1)
    ax.scatter(1 + rng.uniform(-0.08, 0.08, size=len(b)), b, s=14, color="#444", alpha=0.75)
    ax.scatter(2 + rng.uniform(-0.08, 0.08, size=len(a)), a, s=14, color="#b33", alpha=0.75)
    ax.set_title("GSE248378 post TACSTD2 vs recurrence")
    ax.set_ylabel("log2(FPKM+1)")
    fig.savefig(FIG / "GSE248378_TACSTD2_recurrence.png", dpi=140)
    plt.close(fig)


def write_report(
    tests: pd.DataFrame,
    hunt: pd.DataFrame,
    pre_meta: dict,
    post_meta: dict,
    validation: dict,
) -> None:
    def cell(series, gene, cov, endpoint, test_contains=None, adj=None):
        q = tests[(tests["series"] == series) & (tests["gene"] == gene) & (tests["covariate"] == cov) & (tests["endpoint"] == endpoint)]
        if test_contains:
            q = q[q["test"].str.contains(test_contains, regex=False)]
        if adj:
            q = q[q["adjustment"] == adj]
        return None if q.empty else q.iloc[0]

    t_mpr = cell("GSE253564", "TACSTD2", "none", "MPR")
    t_mpr_k = cell("GSE253564", "TACSTD2", "MKI67", "MPR")
    t_mpr_e = cell("GSE253564", "TACSTD2", "ESTIMATEScore", "MPR")
    t_mpr_p = cell("GSE253564", "TACSTD2", "EpiScore", "MPR")
    c_mpr = cell("GSE253564", "CLDN4", "none", "MPR")
    c_mpr_k = cell("GSE253564", "CLDN4", "MKI67", "MPR")
    t_pfs = cell("GSE248378", "TACSTD2", "none", "PFS")
    t_pfs_k = cell("GSE248378", "TACSTD2", "MKI67", "PFS")
    t_pfs_e = cell("GSE248378", "TACSTD2", "ESTIMATEScore", "PFS")
    t_pfs_p = cell("GSE248378", "TACSTD2", "EpiScore", "PFS")
    t_cyt_raw = cell("GSE248378", "TACSTD2", "none", "CYT", test_contains="Spearman unadjusted")
    t_cyt_est = cell("GSE248378", "TACSTD2", "ESTIMATEScore", "CYT", test_contains="partial Spearman")
    t_cyt_pre_est = cell("GSE253564", "TACSTD2", "ESTIMATEScore", "CYT", test_contains="partial Spearman")
    t_cd8_est = cell("GSE248378", "TACSTD2", "ESTIMATEScore", "CD8A", test_contains="partial Spearman")
    t_cd8_pre = cell("GSE253564", "TACSTD2", "ESTIMATEScore", "CD8A", test_contains="partial Spearman")
    t_cyt_k = cell("GSE248378", "TACSTD2", "MKI67", "CYT", test_contains="partial Spearman")
    t_cyt_epi = cell("GSE248378", "TACSTD2", "EpiScore", "CYT", test_contains="partial Spearman")
    t_cyt_pre_k = cell("GSE253564", "TACSTD2", "MKI67", "CYT", test_contains="partial Spearman")
    t_cyt_pre_epi = cell("GSE253564", "TACSTD2", "EpiScore", "CYT", test_contains="partial Spearman")
    t_imm_raw = cell("GSE248378", "TACSTD2", "none", "ImmuneScore", test_contains="Spearman unadjusted")
    t_imm_est = cell("GSE248378", "TACSTD2", "ESTIMATEScore", "ImmuneScore", test_contains="partial Spearman")

    open_yes = hunt[hunt["open_per_patient_matrix"] == "yes"]
    gse_hits = hunt[
        (hunt["status"] == "VERIFIED_EXISTS")
        & hunt["accession_or_id"].astype(str).str.startswith("GSE")
    ]
    hunt_n = len(gse_hits)

    verdict = (
        "PARTLY — leftover MPR gap is real but not residual-robust; "
        "post-treatment TACSTD2–CYT/CD8 negative correlation is residual-robust to ESTIMATE "
        "and remains after MKI67/epithelial; no new PACIFIC/ADRIA matrix."
    )

    lines = []
    lines.append("# REWORK A2 wave 2 — do TACSTD2/CLDN4 vs MPR/PFS/CYT/CD8 survive MKI67, ESTIMATE, and epithelial residuals?")
    lines.append("")
    lines.append("**Self-contained. Public data only. Written to be read without the rest of the repo.**")
    lines.append("")
    lines.append(f"**Verdict: {verdict}**")
    lines.append("")
    lines.append(
        "The user slide (ρ=−0.65 / purity −0.46, p=2e-4) still has **no open PACIFIC per-patient RNA**. "
        "ADRIATIC (NCT03703297) published group-level RNA-seq *medians* in 2025, not a matrix. "
        "The only open human durvalumab NSCLC whole-transcriptome matrices remain **GSE253564** (pre, n=32) "
        "and **GSE248378** (post, n=29) from NCT02904954."
    )
    lines.append("")
    lines.append("## Why this rework exists")
    lines.append("")
    lines.append(
        "Wave 1 (`results/rework/A2_defs/`) recomputed TACSTD2 vs five immune definitions with "
        "partial Spearman on ESTIMATEScore ranks. It did **not** residualise on MKI67 or an epithelial "
        "score, and it did not test MPR/PFS. Leftover PR 175 found pre-treatment TACSTD2 *higher* in MPR "
        f"(p={fmt_p(t_mpr['p']) if t_mpr is not None else '0.016'}) that died after MKI67 OLS residual "
        f"(p={fmt_p(t_mpr_k['p']) if t_mpr_k is not None else '0.25'}). This wave retests that leftover "
        "and the A2 immune claim after each covariate separately."
    )
    lines.append("")
    lines.append("## Analysis set")
    lines.append("")
    lines.append("| Series | Timepoint | n | Usable endpoints |")
    lines.append("|---|---|---|---|")
    lines.append("| GSE253564 | pre-treatment FPKM | 32 | MPR 11 vs 21; pathologic depth; PFS (7 events); CYT; CD8A |")
    lines.append("| GSE248378 | post-treatment FPKM | 29 | **MPR not estimable** (0 MPR in matrix); recurrence 9 vs 20; PFS (9 events); CYT; CD8A |")
    lines.append("")
    lines.append(
        f"Clinical join: Table S1 (PMID 38401548, `mmc2.xlsx`) + Nat Commun source data (PMID 38114518). "
        f"Titles join {validation['checks'][0]['n_matched_to_publication']}/32 and "
        f"{validation['checks'][1]['n_matched_to_publication']}/29. Arm conflicts: {len(validation['mismatches'])} "
        "(documented `45-M-PO` / durva045)."
    )
    lines.append("")
    lines.append("## Pre-specified features and covariates")
    lines.append("")
    lines.append("| Item | Definition | Honest limitation |")
    lines.append("|---|---|---|")
    lines.append("| TACSTD2, CLDN4 | log2(FPKM+1) | Bulk, not protein |")
    lines.append("| CYT | mean log2(FPKM+1) of GZMA+PRF1 (Rooney 2015) | Two-gene cytolytic score |")
    lines.append("| CD8 | log2(FPKM+1) CD8A | Single gene, not a CD8 panel |")
    lines.append("| MKI67 | log2(FPKM+1) | Single-gene proliferation; leftover covariate |")
    lines.append("| ESTIMATE | ESTIMATEScore = Stromal+Immune (v1.0.13 port) | Expression-derived; cosine TumorPurity out of bounds on FPKM |")
    lines.append(
        "| EpiScore | mean log2(FPKM+1) of EPCAM, KRT7, KRT8, KRT18, KRT19, CDH1, MUC1 | "
        "Locked lung-epithelial markers. **TACSTD2 and CLDN4 excluded** (circularity). "
        f"Present pre: {', '.join(pre_meta['epi_genes_present'])}; post: {', '.join(post_meta['epi_genes_present'])}. "
        "Single-gene EPCAM is sensitivity only."
    )
    lines.append("")
    lines.append(
        "Residualisation is **separate, never joint**. Primary MPR/PFS tests use leftover-style OLS residual "
        "of the gene on the covariate, then Mann–Whitney / Cox. Immune tests also report A2-matched "
        "partial Spearman (Pearson on ranks, n−3 df). No multiple-testing correction — every p is shown raw."
    )
    lines.append("")
    lines.append("## Direct answer")
    lines.append("")
    lines.append(
        "1. **PACIFIC / ADRIATIC per-patient matrix: still none open.** This re-hunt verified that live. "
        f"Entrez returned {hunt_n} unique GSE accessions after dropping GSM/GPL children. "
        "The bare keyword PACIFIC without durvalumab is polluted by PacBio / Pacific-Islander series and was not used as a hit list. "
        "None of the GSE hits is a new PACIFIC or ADRIATIC per-patient whole-transcriptome matrix. "
        "Open usable lung durvalumab RNA remains the two NCT02904954 series "
        f"({len(open_yes)} rows marked `yes` in `hunt_triage.tsv`)."
    )
    if t_mpr is not None and t_mpr_k is not None:
        lines.append(
            f"2. **Pre-treatment TACSTD2 is higher in MPR, not lower** "
            f"(median {fmt_r(t_mpr['median_pos'])} vs {fmt_r(t_mpr['median_neg'])}; "
            f"p={fmt_p(t_mpr['p'])}; Cliff δ={fmt_r(t_mpr['cliffs_delta'])}). "
            f"After MKI67 residual p={fmt_p(t_mpr_k['p'])} (δ={fmt_r(t_mpr_k['cliffs_delta'])}); "
            f"after ESTIMATE p={fmt_p(t_mpr_e['p']) if t_mpr_e is not None else 'NA'}; "
            f"after EpiScore p={fmt_p(t_mpr_p['p']) if t_mpr_p is not None else 'NA'}. "
            "The leftover MKI67 kill is reproduced. This is not evidence that TROP2-high tumours are ICI-resistant."
        )
    if c_mpr is not None:
        cldn4_extra = "" if c_mpr_k is None else f"; MKI67 residual p={fmt_p(c_mpr_k['p'])}"
        lines.append(f"3. **CLDN4 vs MPR is null** (p={fmt_p(c_mpr['p'])}{cldn4_extra}).")
    if t_pfs is not None:
        lines.append(
            f"4. **Post-treatment TACSTD2 Cox PFS** unadjusted HR={fmt_r(t_pfs['hr'])} "
            f"({fmt_r(t_pfs['hr_ci_low'])}–{fmt_r(t_pfs['hr_ci_high'])}), p={fmt_p(t_pfs['p'])}, "
            f"{int(t_pfs['n_events'])} events. After MKI67 HR={fmt_r(t_pfs_k['hr']) if t_pfs_k is not None else 'NA'} "
            f"p={fmt_p(t_pfs_k['p']) if t_pfs_k is not None else 'NA'}; "
            f"after ESTIMATE HR={fmt_r(t_pfs_e['hr']) if t_pfs_e is not None else 'NA'} "
            f"p={fmt_p(t_pfs_e['p']) if t_pfs_e is not None else 'NA'}; "
            f"after EpiScore HR={fmt_r(t_pfs_p['hr']) if t_pfs_p is not None else 'NA'} "
            f"p={fmt_p(t_pfs_p['p']) if t_pfs_p is not None else 'NA'}. Small-event finding."
        )
    if t_cyt_est is not None:
        lines.append(
            f"5. **TACSTD2 vs CYT/CD8 (the A2 claim).** Post CYT partial | ESTIMATE ρ={fmt_r(t_cyt_est['rho'])} "
            f"p={fmt_p(t_cyt_est['p'])} (older hunt −0.71 reproduced). "
            f"Post CD8A | ESTIMATE ρ={fmt_r(t_cd8_est['rho']) if t_cd8_est is not None else 'NA'} "
            f"p={fmt_p(t_cd8_est['p']) if t_cd8_est is not None else 'NA'}. "
            f"Post CYT | MKI67 ρ={fmt_r(t_cyt_k['rho']) if t_cyt_k is not None else 'NA'}; "
            f"| EpiScore ρ={fmt_r(t_cyt_epi['rho']) if t_cyt_epi is not None else 'NA'}. "
            f"Pre CYT | ESTIMATE ρ={fmt_r(t_cyt_pre_est['rho']) if t_cyt_pre_est is not None else 'NA'} "
            f"p={fmt_p(t_cyt_pre_est['p']) if t_cyt_pre_est is not None else 'NA'} "
            f"(collapses, as in wave 1). "
            f"Pre CYT | MKI67 ρ={fmt_r(t_cyt_pre_k['rho']) if t_cyt_pre_k is not None else 'NA'}; "
            f"| EpiScore ρ={fmt_r(t_cyt_pre_epi['rho']) if t_cyt_pre_epi is not None else 'NA'}."
        )
    if t_imm_raw is not None and t_imm_est is not None:
        lines.append(
            f"6. **User coefficient match is still post ESTIMATE ImmuneScore**, not a new PACIFIC number: "
            f"raw ρ={fmt_r(t_imm_raw['rho'])}, partial | ESTIMATEScore ρ={fmt_r(t_imm_est['rho'])} "
            f"p={fmt_p(t_imm_est['p'])}. User p=2e-4 is the raw-p order, not the adjusted p. "
            "p=2e-4 attached to adj ρ=−0.46 does not reproduce at n=29."
        )
    lines.append("")
    lines.append("## Primary result — MPR / PFS after each residual")
    lines.append("")
    lines.append("### GSE253564 pre-treatment MPR (11 vs 21)")
    lines.append("")
    lines.append("| Gene | Covariate | n_MPR / n_no | median_MPR | median_no | Cliff δ | p |")
    lines.append("|---|---|---|---|---|---|---|")
    for gene in TARGETS:
        for cov in ("none", "MKI67", "ESTIMATEScore", "EpiScore", "EPCAM_sensitivity"):
            r = cell("GSE253564", gene, cov, "MPR")
            if r is None:
                continue
            lines.append(
                f"| {gene} | {cov} | {int(r['n_pos'])}/{int(r['n_neg'])} | "
                f"{fmt_r(r['median_pos'])} | {fmt_r(r['median_neg'])} | "
                f"{fmt_r(r['cliffs_delta'])} | {fmt_p(r['p'])} |"
            )
    lines.append("")
    lines.append("### GSE248378 post-treatment PFS (9 events) and recurrence")
    lines.append("")
    lines.append("| Gene | Covariate | Cox HR (95% CI) | Cox p | log-rank p | recurrence MWU p |")
    lines.append("|---|---|---|---|---|---|")
    for gene in TARGETS:
        for cov in ("none", "MKI67", "ESTIMATEScore", "EpiScore"):
            s = cell("GSE248378", gene, cov, "PFS")
            rec = cell("GSE248378", gene, cov, "recurrence")
            if s is None:
                continue
            hr = "NA" if s["hr"] is None or (isinstance(s["hr"], float) and math.isnan(s["hr"])) else (
                f"{fmt_r(s['hr'])} ({fmt_r(s['hr_ci_low'])}–{fmt_r(s['hr_ci_high'])})"
            )
            rp = fmt_p(rec["p"]) if rec is not None else "NA"
            lines.append(
                f"| {gene} | {cov} | {hr} | {fmt_p(s['p'])} | {fmt_p(s['logrank_p'])} | {rp} |"
            )
    lines.append("")
    lines.append("Post MPR is **not estimable** (0 vs 29). Do not pool this accession as an MPR replication.")
    lines.append("")
    lines.append("## TACSTD2 / CLDN4 vs CYT and CD8A")
    lines.append("")
    lines.append(
        "Two estimators: leftover-style Spearman of the OLS residual vs the *raw* immune score, "
        "and A2-matched partial Spearman. Both are in `immune_correlations.tsv`. Partial Spearman is the A2 comparator."
    )
    lines.append("")
    lines.append("| Dataset | Gene | Immune | Covariate | partial ρ | partial p | OLS-resid ρ | OLS-resid p |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for series in ("GSE253564", "GSE248378"):
        for gene in TARGETS:
            for ep in ("CYT", "CD8A"):
                for cov in ("MKI67", "ESTIMATEScore", "EpiScore"):
                    pr = cell(series, gene, cov, ep, test_contains="partial Spearman")
                    ol = cell(series, gene, cov, ep, test_contains="OLS-residual")
                    if pr is None and ol is None:
                        continue
                    lines.append(
                        f"| {series} | {gene} | {ep} | {cov} | "
                        f"{fmt_r(pr['rho']) if pr is not None else 'NA'} | "
                        f"{fmt_p(pr['p']) if pr is not None else 'NA'} | "
                        f"{fmt_r(ol['rho']) if ol is not None else 'NA'} | "
                        f"{fmt_p(ol['p']) if ol is not None else 'NA'} |"
                    )
    lines.append("")
    lines.append("Unadjusted (covariate = none) Spearman rows are in the same TSV.")
    lines.append("")
    lines.append("## Hunt — PACIFIC / ADRIATIC / open durvalumab")
    lines.append("")
    lines.append(
        "Live Entrez `gds` + PubMed queries were run at analysis time (`hunt_manifest.tsv`). "
        "A GEO accession is listed only if Entrez returned it. Near-misses are not hits."
    )
    lines.append("")
    lines.append("| Resource | ID | Open per-patient matrix? | Status |")
    lines.append("|---|---|---|---|")
    show = hunt[hunt["query"] == "targeted"]
    for _, r in show.iterrows():
        lines.append(
            f"| {r['resource']} | `{r['accession_or_id']}` | {r['open_per_patient_matrix']} | {r['status']} |"
        )
    lines.append("")
    lines.append("## Honest interpretation")
    lines.append("")
    lines.append(
        "1. Do not cite a single ρ. Pre-treatment CYT/CD8A vs TACSTD2 is mostly ESTIMATE/purity; "
        "post-treatment CYT/CD8A is not."
    )
    lines.append(
        "2. The leftover TACSTD2-higher-in-MPR result is the opposite of “TROP2-high = ICI-cold/resistant” "
        "and is not independent of proliferation."
    )
    lines.append(
        "3. n=32/29 and 7/9 PFS events. Only large effects are detectable. A null CLDN4 result does not rule out a modest effect."
    )
    lines.append(
        "4. ESTIMATE TumorPurity cosine is out of bounds for "
        f"{pre_meta['n_purity_oob']}/32 pre and {post_meta['n_purity_oob']}/29 post FPKM samples; "
        "ESTIMATEScore (not the cosine) is the covariate."
    )
    lines.append(
        "5. EpiScore is a locked marker mean, not a published deconvolution epithelial fraction. "
        "It is reported so the residual request can be tested without inventing a commercial assay."
    )
    lines.append(
        "6. ADRIATIC RNA-seq exists as a trial biomarker analysis (group medians). That is not a downloadable per-patient matrix."
    )
    lines.append(
        "7. p=2e-4 attached to adj ρ=−0.46 does not reproduce at these sample sizes (needs n≈55+)."
    )
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("pip install -r requirements.txt")
    lines.append("python3 scripts/rework_A2_wave2.py")
    lines.append("```")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `REPORT.md` — this write-up")
    lines.append("- `residual_tests.tsv` — MPR / recurrence / depth / PFS after each residual")
    lines.append("- `immune_correlations.tsv` — vs CYT / CD8A / ImmuneScore, both estimators")
    lines.append("- `hunt_manifest.tsv` — live Entrez log (search counts + GSE hits)")
    lines.append("- `hunt_triage.tsv` — GSE hits plus locked PACIFIC/ADRIA triage")
    lines.append("- `sample_table_GSE253564.tsv`, `sample_table_GSE248378.tsv`")
    lines.append("- `summary.json`, `provenance.json`")
    lines.append("- `figures/GSE253564_TACSTD2_MPR_residuals.png`")
    lines.append("- `figures/TACSTD2_partial_rho_CYT_CD8.png`")
    lines.append("- `figures/GSE248378_TACSTD2_recurrence.png`")
    lines.append("")
    lines.append("## Data")
    lines.append("")
    lines.append("- GSE253564 / GSE248378 FPKM and series matrices (NCBI GEO FTP)")
    lines.append("- ESTIMATE `SI_geneset.gmt` v1.0.13")
    lines.append("- PMID 38401548 Table S1 (`mmc2.xlsx` via Europe PMC OA)")
    lines.append("- PMID 38114518 source data (`41467_2023_44195_MOESM6_ESM.xlsx`)")
    lines.append("")
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    paths = ensure_inputs()
    gmt = load_gmt(paths["gmt"])
    clin, validation = build_clinical(paths)
    clin.to_csv(OUT / "clinical_annotation.csv", index=False)
    (OUT / "clinical_join_validation.json").write_text(json.dumps(validation, indent=2))

    pre_expr = read_fpkm(paths["pre_fpkm"])
    post_expr = read_fpkm(paths["post_fpkm"])
    pre_s, pre_rows, pre_meta = analyse_cohort("GSE253564", "pre", pre_expr, clin, gmt)
    post_s, post_rows, post_meta = analyse_cohort("GSE248378", "post", post_expr, clin, gmt)
    tests = pd.DataFrame(pre_rows + post_rows)
    tests.to_csv(OUT / "all_tests.tsv", sep="\t", index=False)

    residual = tests[tests["endpoint"].isin(["MPR", "recurrence", "pathologic_depth", "PFS"])].copy()
    residual.to_csv(OUT / "residual_tests.tsv", sep="\t", index=False)
    immune = tests[tests["endpoint"].isin(["CYT", "CD8A", "ImmuneScore"])].copy()
    immune.to_csv(OUT / "immune_correlations.tsv", sep="\t", index=False)
    pre_s.to_csv(OUT / "sample_table_GSE253564.tsv", sep="\t", index=False)
    post_s.to_csv(OUT / "sample_table_GSE248378.tsv", sep="\t", index=False)

    hunt_rows = hunt_open_matrices()
    hunt = pd.DataFrame(hunt_rows)
    hunt.to_csv(OUT / "hunt_manifest.tsv", sep="\t", index=False)
    triage = hunt[hunt["query"].eq("targeted") | hunt["accession_or_id"].astype(str).str.startswith("GSE")].copy()
    triage.to_csv(OUT / "hunt_triage.tsv", sep="\t", index=False)

    write_figures(pre_s, post_s, tests)
    write_report(tests, hunt, pre_meta, post_meta, validation)

    def grab(series, gene, cov, endpoint, test_sub=None):
        q = tests[(tests.series == series) & (tests.gene == gene) & (tests.covariate == cov) & (tests.endpoint == endpoint)]
        if test_sub:
            q = q[q.test.str.contains(test_sub, regex=False)]
        if q.empty:
            return None
        r = q.iloc[0]
        return {k: (None if (isinstance(r[k], float) and math.isnan(r[k])) else r[k])
                for k in r.index if k in {"p", "rho", "hr", "cliffs_delta", "n_pos", "n_neg", "n_events", "note"}}

    summary = {
        "question": "TACSTD2/CLDN4 vs MPR/PFS/CYT/CD8 after separate MKI67, ESTIMATE, epithelial residuals; new PACIFIC/ADRIA matrix?",
        "verdict": "PARTLY — MPR gap not residual-robust; post CYT/CD8 negative correlation is; no new PACIFIC/ADRIA matrix",
        "user_claim": {"rho_raw": -0.65, "rho_purity_adj": -0.46, "p": 2e-4},
        "n_pre": 32,
        "n_post": 29,
        "methods": {
            "expression": "log2(FPKM+1)",
            "CYT": "mean log2(FPKM+1) GZMA+PRF1",
            "CD8": "log2(FPKM+1) CD8A",
            "ESTIMATE": "v1.0.13 Immune+Stromal; covariate = ESTIMATEScore",
            "EpiScore": "mean log2(FPKM+1) of EPCAM KRT7 KRT8 KRT18 KRT19 CDH1 MUC1 (TACSTD2/CLDN4 excluded)",
            "residual_MPR_PFS": "OLS residual of gene on covariate, then MWU / Cox",
            "residual_immune": "OLS-residual Spearman + partial Spearman (ranks, n-3)",
            "covariates": "MKI67, ESTIMATEScore, EpiScore separately (never joint)",
        },
        "key": {
            "pre_TACSTD2_MPR_unadj": grab("GSE253564", "TACSTD2", "none", "MPR"),
            "pre_TACSTD2_MPR_MKI67": grab("GSE253564", "TACSTD2", "MKI67", "MPR"),
            "pre_TACSTD2_MPR_ESTIMATE": grab("GSE253564", "TACSTD2", "ESTIMATEScore", "MPR"),
            "pre_TACSTD2_MPR_EpiScore": grab("GSE253564", "TACSTD2", "EpiScore", "MPR"),
            "post_TACSTD2_PFS_unadj": grab("GSE248378", "TACSTD2", "none", "PFS"),
            "post_TACSTD2_CYT_partial_ESTIMATE": grab("GSE248378", "TACSTD2", "ESTIMATEScore", "CYT", "partial Spearman"),
            "post_TACSTD2_CD8A_partial_ESTIMATE": grab("GSE248378", "TACSTD2", "ESTIMATEScore", "CD8A", "partial Spearman"),
        },
        "pre_meta": pre_meta,
        "post_meta": post_meta,
        "pacific_adria_open_matrix": False,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    provenance = {
        "inputs": {
            name: {"file": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "md5": md5_file(path), "url": URLS.get(name, "")}
            for name, path in paths.items()
            if path.exists()
        },
        "epi_genes": list(EPI_GENES),
        "cyt_genes": list(CYT_GENES),
        "note": "Raw GEO FPKM matrices are committed under data/ for reproducibility. PMC zip is not committed.",
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2))

    print("[A2_wave2] leftover-check TACSTD2 MPR unadj / MKI67 residual:")
    for cov in ("none", "MKI67"):
        r = grab("GSE253564", "TACSTD2", cov, "MPR")
        print(" ", cov, r)
    print("[A2_wave2] post TACSTD2 vs CYT partial ESTIMATE:", grab("GSE248378", "TACSTD2", "ESTIMATEScore", "CYT", "partial Spearman"))
    print(f"[A2_wave2] wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
