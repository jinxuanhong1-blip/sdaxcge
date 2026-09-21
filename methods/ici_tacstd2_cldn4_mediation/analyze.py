#!/usr/bin/env python3
"""TACSTD2 response and TMB/immune associations after CLDN4, in open ICI cohorts.

Question
    In open ICI cohorts that measure both TACSTD2 and CLDN4, does the
    TACSTD2 association with clinical benefit, with TMB, or with immune
    scores shrink once CLDN4 is in the model? Is there a regression-imputation
    indirect effect through CLDN4?

What is not estimated
    OAK and POPLAR linear RNA-seq (EGA EGAS00001005013, datasets
    EGAD00001008390/08391 and the matching counts/CPM/clinical files) are
    controlled by DAC EGAC00001002120. This script does not log into EGA and
    does not compute an OAK or POPLAR odds ratio.

Pre-specified models (no quantile search, cohorts not pooled)
    Expression is log2, then z-scored inside the complete-case sample (sd with
    ddof=1). The primary response contrast is a logistic odds ratio per 1 SD
    of TACSTD2, unadjusted and with z(CLDN4). A median split is reported as a
    2x2 odds ratio; a zero cell stays undefined (no Haldane fill-in).
    Continuous TMB and immune outcomes use Spearman plus a standardized linear
    coefficient. Mediation on a binary outcome is the regression-imputation
    natural indirect effect on the probability scale for a +1 SD shift in
    TACSTD2 (Imai residual imputation, no extra Monte Carlo inside the point
    estimate). The bootstrap resamples patients and re-fits every step.

    Proportion mediated is stored only as a number. The write-up interprets it
    only when the total-effect interval excludes zero.
"""

from __future__ import annotations

import gzip
import json
import math
import urllib.request
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, mannwhitneyu, norm, pearsonr, rankdata, spearmanr, t as student_t

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "data" / "cache"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
N_BOOT = 2000
SEED = 814

TEFF = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]
GENES = ["TACSTD2", "CLDN4", *TEFF]
ENSEMBL = {
    "ENSG00000184292": "TACSTD2",
    "ENSG00000189143": "CLDN4",
    "ENSG00000153563": "CD8A",
    "ENSG00000145649": "GZMA",
    "ENSG00000100453": "GZMB",
    "ENSG00000111537": "IFNG",
    "ENSG00000163508": "EOMES",
    "ENSG00000138755": "CXCL9",
    "ENSG00000169245": "CXCL10",
    "ENSG00000073861": "TBX21",
}

GEO = {
    "GSE126044_counts.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
    "GSE126044_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz",
    "GSE135222_exp.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
    "GSE135222_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz",
    "GSE166449_TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz",
    "GSE166449_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/matrix/GSE166449_series_matrix.txt.gz",
    "GSE190265_TPM.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/suppl/GSE190265_TPM_France3.csv.gz",
    "GSE190265_info.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/suppl/GSE190265_samples_info_France3.csv.gz",
    "GSE190266_TPM.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190266/suppl/GSE190266_TPM_France4.csv.gz",
    "GSE190266_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190266/matrix/GSE190266_series_matrix.txt.gz",
    "GSE207422_log2TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
    "GSE207422_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/matrix/GSE207422_series_matrix.txt.gz",
    "GSE218989_TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/suppl/GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz",
    "GSE218989_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/matrix/GSE218989_series_matrix.txt.gz",
    "cds.RData": "https://raw.githubusercontent.com/SiYangming/IMvigor210CoreBiologies/master/data/cds.RData",
}


def fetch(name: str) -> Path:
    dest = CACHE / name
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    CACHE.mkdir(parents=True, exist_ok=True)
    url = GEO[name]
    print("GET", url)
    urllib.request.urlretrieve(url, dest)
    return dest


def series_samples(path: Path) -> list[dict]:
    buckets: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            buckets.setdefault(key, []).append([v.strip().strip('"') for v in vals])
    titles = buckets["!Sample_title"][0]
    geos = buckets["!Sample_geo_accession"][0]
    chars = buckets.get("!Sample_characteristics_ch1", [])
    descs = buckets.get("!Sample_description", [[]])
    out = []
    for i, title in enumerate(titles):
        row = {
            "title": title,
            "gsm": geos[i],
            "description": descs[0][i] if descs and i < len(descs[0]) else "",
        }
        for crow in chars:
            if i < len(crow) and ":" in crow[i]:
                k, v = crow[i].split(":", 1)
                row[k.strip()] = v.strip()
        out.append(row)
    return out


def _float(token: str) -> float:
    token = token.strip().replace(",", ".")
    if token in {"", "NA", "NaN", "nan", "None"}:
        return math.nan
    return float(token)


def load_genes_by_sample(path: Path, *, sep: str = "\t", ensembl: bool = False) -> tuple[list[str], dict[str, np.ndarray]]:
    found: dict[str, list[np.ndarray]] = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split(sep)
        samples = [h.strip().strip('"') for h in header[1:]]
        n = len(samples)
        for line in handle:
            parts = line.rstrip("\n").split(sep)
            if len(parts) < 2:
                continue
            gene = parts[0].strip().strip('"')
            if ensembl:
                gene = ENSEMBL.get(gene.split(".")[0], "")
            if gene not in GENES:
                continue
            vals = np.array([_float(x) for x in parts[1 : 1 + n]], dtype=float)
            if vals.size != n:
                continue
            found.setdefault(gene, []).append(vals)
    data = {g: np.nanmean(np.vstack(vs), axis=0) for g, vs in found.items()}
    return samples, data


def library_sizes(path: Path) -> np.ndarray:
    with gzip.open(path, "rt", errors="replace") as handle:
        samples = handle.readline().rstrip("\n").split("\t")[1:]
        libs = np.zeros(len(samples), dtype=float)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            vals = parts[1 : 1 + len(samples)]
            libs += np.array([_float(x) if x.strip() not in {"", "NA"} else 0.0 for x in vals], dtype=float)
    return libs


def to_log(values: np.ndarray, *, already_log: bool, cpm_libs: np.ndarray | None = None) -> np.ndarray:
    if already_log:
        return values.astype(float)
    x = values.astype(float)
    if cpm_libs is not None:
        x = x / cpm_libs * 1e6
    return np.log2(x + 1.0)


def teff_score(logged: dict[str, np.ndarray]) -> tuple[np.ndarray, list[str]]:
    rows = []
    present = []
    for gene in TEFF:
        if gene not in logged:
            continue
        row = logged[gene]
        if np.sum(np.isfinite(row)) < 5:
            continue
        sd = np.nanstd(row, ddof=1)
        if not np.isfinite(sd) or sd == 0:
            continue
        present.append(gene)
        rows.append((row - np.nanmean(row)) / sd)
    if len(rows) < 6:
        n = len(next(iter(logged.values())))
        return np.full(n, np.nan), present
    return np.nanmean(np.vstack(rows), axis=0), present


def frame_from_logs(
    samples: list[str],
    logged: dict[str, np.ndarray],
    clinical: dict[str, dict],
) -> pd.DataFrame:
    teff, present = teff_score(logged)
    rows = []
    for i, sid in enumerate(samples):
        clin = clinical.get(sid)
        if clin is None:
            continue
        rec = {
            "sample": sid,
            "TACSTD2": float(logged["TACSTD2"][i]) if "TACSTD2" in logged else math.nan,
            "CLDN4": float(logged["CLDN4"][i]) if "CLDN4" in logged else math.nan,
            "CD8A": float(logged["CD8A"][i]) if "CD8A" in logged else math.nan,
            "Teff": float(teff[i]),
            "teff_genes": ",".join(present),
            **clin,
        }
        rows.append(rec)
    return pd.DataFrame(rows)


def benefit_from_title(title: str) -> int | None:
    t = title.lower().replace(" ", "").replace("_", "")
    if "nonresponder" in t or "non-responder" in t:
        return 0
    if "responder" in t:
        return 1
    return None


def load_gse126044() -> pd.DataFrame:
    path = fetch("GSE126044_counts.txt.gz")
    samples, data = load_genes_by_sample(path)
    libs = library_sizes(path)
    logged = {g: to_log(v, already_log=False, cpm_libs=libs) for g, v in data.items()}
    clin = {}
    for rec in series_samples(fetch("GSE126044_series_matrix.txt.gz")):
        sid = rec["title"].replace("RNA-seq_", "")
        label = rec.get("patient response", "")
        y = 1 if label == "responder" else 0 if label == "non-responder" else None
        clin[sid] = {"benefit": y, "gsm": rec["gsm"]}
    return frame_from_logs(samples, logged, clin)


def load_gse135222() -> pd.DataFrame:
    samples, data = load_genes_by_sample(fetch("GSE135222_exp.tsv.gz"), ensembl=True)
    logged = {g: to_log(v, already_log=False) for g, v in data.items()}
    clin = {}
    for rec in series_samples(fetch("GSE135222_series_matrix.txt.gz")):
        sid = rec["title"].replace(" ", "")
        time = float(rec["pfs.time"])
        event = int(rec["progression-free survival (pfs)"])
        if time >= 180:
            y = 1
        elif event == 1:
            y = 0
        else:
            y = None
        clin[sid] = {"benefit": y, "gsm": rec["gsm"], "pfs_time": time, "pfs_event": event}
    return frame_from_logs(samples, logged, clin)


def load_gse166449() -> pd.DataFrame:
    samples, data = load_genes_by_sample(fetch("GSE166449_TPM.txt.gz"))
    # Deposited CLDN4 max is about 3. Values are already on a log scale.
    logged = {g: to_log(v, already_log=True) for g, v in data.items()}
    clin = {}
    for rec in series_samples(fetch("GSE166449_series_matrix.txt.gz")):
        y = benefit_from_title(rec["title"])
        clin[rec["description"]] = {"benefit": y, "gsm": rec["gsm"]}
    return frame_from_logs(samples, logged, clin)


def load_gse190265() -> pd.DataFrame:
    path = fetch("GSE190265_TPM.csv.gz")
    with gzip.open(path, "rt", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split(";")
        index = {g: i for i, g in enumerate(header) if g in GENES}
        samples = []
        cols = {g: [] for g in index}
        for line in handle:
            parts = line.rstrip("\n").split(";")
            if len(parts) != len(header) + 1:
                raise SystemExit(f"GSE190265 row width {len(parts)} vs header {len(header)}")
            samples.append(parts[0])
            for gene, i in index.items():
                cols[gene].append(_float(parts[i + 1]))
    data = {g: np.array(vs, dtype=float) for g, vs in cols.items()}
    logged = {g: to_log(v, already_log=False) for g, v in data.items()}
    info = pd.read_csv(fetch("GSE190265_info.csv.gz"), sep=";", compression="gzip")
    clin = {}
    for rec in info.to_dict(orient="records"):
        time = float(rec["time_PFS"])
        event = int(rec["evtPFS"])
        if time >= 6:
            y = 1
        elif event == 1:
            y = 0
        else:
            y = None
        clin[str(rec["sample"])] = {"benefit": y, "gsm": "", "pfs_time": time, "pfs_event": event}
    return frame_from_logs(samples, logged, clin)


def load_gse207422() -> pd.DataFrame:
    samples, data = load_genes_by_sample(fetch("GSE207422_log2TPM.txt.gz"))
    logged = {g: to_log(v, already_log=True) for g, v in data.items()}
    clin = {}
    for rec in series_samples(fetch("GSE207422_series_matrix.txt.gz")):
        if rec.get("description") != "Bulk RNAseq":
            continue
        if rec.get("sampling_time") != "Pre-treatment biopsy":
            continue
        path_resp = rec.get("pathologic_response", "")
        if path_resp.startswith("MPR"):
            y = 1
        elif path_resp == "NMPR":
            y = 0
        else:
            continue
        clin[rec["title"]] = {"benefit": y, "gsm": rec["gsm"], "pathologic_response": path_resp}
    return frame_from_logs(samples, logged, clin)


def load_gse218989() -> pd.DataFrame:
    samples, data = load_genes_by_sample(fetch("GSE218989_TPM.txt.gz"))
    logged = {g: to_log(v, already_log=False) for g, v in data.items()}
    clin = {}
    for rec in series_samples(fetch("GSE218989_series_matrix.txt.gz")):
        label = rec.get("treatment outcome", "")
        y = 1 if label == "Responder" else 0 if label == "Non-responder" else None
        clin[rec["title"]] = {
            "benefit": y,
            "gsm": rec["gsm"],
            "treatment": rec.get("treatment", ""),
            "tissue": rec.get("tissue", ""),
        }
    return frame_from_logs(samples, logged, clin)


def gse190266_tacstd2_absent() -> str:
    path = fetch("GSE190266_TPM.csv.gz")
    with gzip.open(path, "rt", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split(";")
    has = "TACSTD2" in header
    last = header[-1]
    if has:
        raise SystemExit("GSE190266 header now contains TACSTD2; the absence note is stale")
    return (
        f"Deposited TPM header has {len(header)} fields and stops at {last}. "
        "TACSTD2 is not in the file, so this cohort is not scored."
    )


def extract_imvigor() -> pd.DataFrame:
    dest = CACHE / "imvigor210_selected.csv"
    rdata = fetch("cds.RData")
    if not dest.exists() or dest.stat().st_mtime < rdata.stat().st_mtime:
        r_code = r"""
load("%s")
counts <- get("counts", envir = attr(cds, "assayData"))
pdat <- attr(attr(cds, "phenoData"), "data")
fdat <- attr(attr(cds, "featureData"), "data")
wanted <- c(TACSTD2=4070, CLDN4=1364, CD8A=925, GZMA=3001, GZMB=3002, IFNG=3458, EOMES=8320, CXCL9=4283, CXCL10=3627, TBX21=30009)
for (nm in names(wanted)) {
  hit <- rownames(fdat)[which(fdat$symbol == nm)]
  if (length(hit) != 1 || hit != as.character(wanted[[nm]])) stop(paste("symbol mismatch", nm, paste(hit, collapse=",")))
}
sel <- t(counts[as.character(wanted), , drop = FALSE])
colnames(sel) <- names(wanted)
keep <- c("binaryResponse", "Immune phenotype", "FMOne mutation burden per MB", "Neoantigen burden per MB", "sizeFactor", "os", "censOS", "Best Confirmed Overall Response")
out <- data.frame(sample = rownames(pdat), sel, pdat[, keep], check.names = FALSE)
write.csv(out, "%s", row.names = FALSE)
""" % (rdata.as_posix(), dest.as_posix())
        script = CACHE / "extract_imvigor.R"
        script.write_text(r_code)
        import subprocess

        subprocess.check_call(["Rscript", str(script)])
    df = pd.read_csv(dest)
    if len(df) != 348:
        raise SystemExit(f"IMvigor210 sample count is {len(df)}, expected 348")
    sf = df["sizeFactor"].astype(float)
    if (sf <= 0).any():
        raise SystemExit("non-positive IMvigor210 sizeFactor")
    for gene in GENES:
        df[gene] = np.log2(df[gene].astype(float) / sf + 1.0)
    logged = {g: df[g].to_numpy(float) for g in TEFF}
    teff, present = teff_score(logged)
    df["Teff"] = teff
    df["teff_genes"] = ",".join(present)
    df["benefit"] = df["binaryResponse"].map({"CR/PR": 1, "SD/PD": 0})
    df["tmb"] = pd.to_numeric(df["FMOne mutation burden per MB"], errors="coerce")
    df["neoantigen"] = pd.to_numeric(df["Neoantigen burden per MB"], errors="coerce")
    df["immune"] = df["Immune phenotype"]
    df["os_time"] = pd.to_numeric(df["os"], errors="coerce")
    df["os_event"] = pd.to_numeric(df["censOS"], errors="coerce")
    both = df[["tmb", "neoantigen"]].dropna()
    rho, p = spearmanr(both["tmb"], both["neoantigen"])
    if not (rho > 0.4):
        raise SystemExit(f"IMvigor TMB vs neoantigen Spearman {rho} failed the column check")
    df.attrs["tmb_neo_rho"] = float(rho)
    df.attrs["tmb_neo_p"] = float(p)
    df.attrs["tmb_neo_n"] = int(len(both))
    return df


def oak_poplar_status() -> dict:
    url = "https://ega-archive.org/datasets/EGAD00001008391"
    req = urllib.request.Request(url, headers={"User-Agent": "sdaxcge-public-ici"})
    with urllib.request.urlopen(req, timeout=40) as resp:
        html = resp.read().decode("utf-8", errors="replace")
        code = resp.status
    dac = "EGAC00001002120" if "EGAC00001002120" in html else ""
    if dac == "" or "Log2 normalized TPM matrix for OAK" not in html:
        raise SystemExit("EGA dataset page did not confirm the OAK log2 TPM matrix and DAC")
    return {
        "study": "EGAS00001005013",
        "checked_url": url,
        "http_status": code,
        "dataset": "EGAD00001008391",
        "dataset_label_on_page": "Log2 normalized TPM matrix for OAK (GO28915)",
        "dac": dac,
        "downloaded": False,
        "note": (
            "The public dataset page names the OAK log2(TPM+1) matrix and DAC "
            "EGAC00001002120. No EGA login was used and the matrix was not downloaded. "
            "POPLAR log2(TPM+1) is EGAD00001008390 on the same study page. "
            "No OAK or POPLAR odds ratio was computed."
        ),
    }


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = np.std(x, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.full_like(x, np.nan)
    return (x - np.mean(x)) / sd


def ols(y: np.ndarray, X: np.ndarray) -> dict | None:
    n, k = X.shape
    if n <= k + 2:
        return None
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - k
    sigma2 = float(np.sum(resid**2) / dof)
    try:
        cov = sigma2 * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return None
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = beta / se
    p = 2 * student_t.sf(np.abs(tstat), dof)
    return {"beta": beta, "se": se, "p": p, "n": n, "dof": dof}


def logit_fit(y: np.ndarray, X: np.ndarray):
    if len(np.unique(y)) < 2 or len(y) <= X.shape[1] + 2:
        return None
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            warnings.filterwarnings("ignore", message=".*[Cc]onverg.*")
            warnings.filterwarnings("ignore", message=".*[Ss]eparat.*")
            res = sm.Logit(y, X).fit(disp=False, maxiter=100)
    except Exception:
        return None
    if not res.mle_retvals.get("converged", False):
        return None
    if np.any(~np.isfinite(res.params)) or np.any(np.abs(res.params) > 20):
        return None
    return res


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 5:
        return math.nan, math.nan, int(mask.sum())
    r, p = spearmanr(x[mask], y[mask])
    return float(r), float(p), int(mask.sum())


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(mask.sum())
    if n < 6:
        return math.nan, math.nan, n
    rx, ry, rz = rankdata(x[mask]), rankdata(y[mask]), rankdata(z[mask])

    def resid(a, b):
        design = np.column_stack([np.ones(len(b)), b])
        coef, *_ = np.linalg.lstsq(design, a, rcond=None)
        return a - design @ coef

    r, p = pearsonr(resid(rx, rz), resid(ry, rz))
    return float(r), float(p), n


def woolf(a: int, b: int, c: int, d: int) -> dict:
    table = np.array([[a, b], [c, d]])
    oddsratio, p = fisher_exact(table)
    if min(a, b, c, d) == 0:
        return {
            "a_resp_high": a,
            "b_nonresp_high": b,
            "c_resp_low": c,
            "d_nonresp_low": d,
            "OR": None,
            "OR_lo": None,
            "OR_hi": None,
            "fisher_p": float(p),
            "zero_cell": True,
        }
    log_or = math.log((a * d) / (b * c))
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return {
        "a_resp_high": a,
        "b_nonresp_high": b,
        "c_resp_low": c,
        "d_nonresp_low": d,
        "OR": float(math.exp(log_or)),
        "OR_lo": float(math.exp(log_or - 1.96 * se)),
        "OR_hi": float(math.exp(log_or + 1.96 * se)),
        "fisher_p": float(p),
        "zero_cell": False,
    }


def median_split_or(y: np.ndarray, x: np.ndarray) -> dict:
    med = float(np.median(x))
    high = x >= med
    a = int(np.sum(high & (y == 1)))
    b = int(np.sum(high & (y == 0)))
    c = int(np.sum(~high & (y == 1)))
    d = int(np.sum(~high & (y == 0)))
    out = woolf(a, b, c, d)
    out["median"] = med
    out["n_high"] = int(high.sum())
    out["n_low"] = int((~high).sum())
    out["n"] = int(len(y))
    return out


def prob_mediation(y: np.ndarray, tac: np.ndarray, cldn: np.ndarray) -> dict | None:
    zt = zscore(tac)
    zm = zscore(cldn)
    if not np.all(np.isfinite(zt)) or not np.all(np.isfinite(zm)):
        return None
    design_m = np.column_stack([np.ones(len(zt)), zt])
    a_fit = ols(zm, design_m)
    if a_fit is None:
        return None
    a0, a1 = a_fit["beta"]
    resid = zm - design_m @ a_fit["beta"]
    unadj = logit_fit(y, np.column_stack([np.ones(len(y)), zt]))
    adj = logit_fit(y, np.column_stack([np.ones(len(y)), zt, zm]))
    if unadj is None or adj is None:
        return None
    b0, bt, bm = (float(v) for v in adj.params)

    def mu(t, m):
        eta = np.clip(b0 + bt * t + bm * m, -50, 50)
        return 1.0 / (1.0 + np.exp(-eta))

    t1 = zt + 1.0
    m1 = zm + a1  # residual imputation: M(t+1) = M(t) + a1
    nie = float(np.mean(mu(t1, m1) - mu(t1, zm)))
    nde = float(np.mean(mu(t1, zm) - mu(zt, zm)))
    total = float(np.mean(mu(t1, m1) - mu(zt, zm)))
    inter = logit_fit(y, np.column_stack([np.ones(len(y)), zt, zm, zt * zm]))
    inter_p = float(inter.pvalues[3]) if inter is not None else math.nan
    r_pear, _ = pearsonr(zt, zm)
    return {
        "n": int(len(y)),
        "n_benefit": int(y.sum()),
        "logOR_unadj": float(unadj.params[1]),
        "logOR_unadj_se": float(unadj.bse[1]),
        "logOR_unadj_p": float(unadj.pvalues[1]),
        "OR_unadj": float(math.exp(unadj.params[1])),
        "OR_unadj_lo": float(math.exp(unadj.conf_int()[1, 0])),
        "OR_unadj_hi": float(math.exp(unadj.conf_int()[1, 1])),
        "logOR_adj": float(adj.params[1]),
        "logOR_adj_p": float(adj.pvalues[1]),
        "OR_adj": float(math.exp(adj.params[1])),
        "OR_adj_lo": float(math.exp(adj.conf_int()[1, 0])),
        "OR_adj_hi": float(math.exp(adj.conf_int()[1, 1])),
        "logOR_CLDN4_given_TACSTD2": float(adj.params[2]),
        "OR_CLDN4": float(math.exp(adj.params[2])),
        "OR_CLDN4_lo": float(math.exp(adj.conf_int()[2, 0])),
        "OR_CLDN4_hi": float(math.exp(adj.conf_int()[2, 1])),
        "OR_CLDN4_p": float(adj.pvalues[2]),
        "a_CLDN4_on_TACSTD2": float(a1),
        "a_p": float(a_fit["p"][1]),
        "pearson_r": float(r_pear),
        "vif": float(1.0 / (1.0 - r_pear**2)) if abs(r_pear) < 0.999 else math.inf,
        "interaction_p": inter_p,
        "NIE": nie,
        "NDE": nde,
        "total_risk_diff": total,
        "diff_logOR": float(unadj.params[1] - adj.params[1]),
    }


def linear_pair(y: np.ndarray, tac: np.ndarray, cldn: np.ndarray) -> dict | None:
    mask = np.isfinite(y) & np.isfinite(tac) & np.isfinite(cldn)
    y, tac, cldn = y[mask], tac[mask], cldn[mask]
    if len(y) < 8:
        return None
    zy, zt, zm = zscore(y), zscore(tac), zscore(cldn)
    unadj = ols(zy, np.column_stack([np.ones(len(zy)), zt]))
    adj = ols(zy, np.column_stack([np.ones(len(zy)), zt, zm]))
    a_fit = ols(zm, np.column_stack([np.ones(len(zm)), zt]))
    if unadj is None or adj is None or a_fit is None:
        return None
    a1 = float(a_fit["beta"][1])
    b1 = float(adj["beta"][2])
    c = float(unadj["beta"][1])
    cp = float(adj["beta"][1])
    rho, rp, n = spearman(tac, y)
    pr, pp, _ = partial_spearman(tac, y, cldn)
    return {
        "n": int(len(y)),
        "spearman": rho,
        "spearman_p": rp,
        "partial_spearman": pr,
        "partial_spearman_p": pp,
        "beta_unadj": c,
        "beta_unadj_p": float(unadj["p"][1]),
        "beta_adj": cp,
        "beta_adj_p": float(adj["p"][1]),
        "beta_CLDN4": b1,
        "beta_CLDN4_p": float(adj["p"][2]),
        "a": a1,
        "indirect": a1 * b1,
        "diff_beta": c - cp,
    }


def boot_indices(n: int, rng: np.random.Generator) -> np.ndarray:
    return rng.integers(0, n, n)


def bootstrap_response(y: np.ndarray, tac: np.ndarray, cldn: np.ndarray, rng: np.random.Generator) -> dict:
    keys = ["logOR_unadj", "logOR_adj", "diff_logOR", "NIE", "NDE", "total_risk_diff", "a_CLDN4_on_TACSTD2"]
    store = {k: [] for k in keys}
    failed = 0
    n = len(y)
    for _ in range(N_BOOT):
        idx = boot_indices(n, rng)
        if y[idx].sum() < 2 or (n - y[idx].sum()) < 2:
            failed += 1
            continue
        fit = prob_mediation(y[idx], tac[idx], cldn[idx])
        if fit is None:
            failed += 1
            continue
        for k in keys:
            store[k].append(fit[k])
    out = {"n_boot": N_BOOT, "n_failed": failed}
    for k, vals in store.items():
        arr = np.array(vals, dtype=float)
        out[k + "_boot_n"] = int(arr.size)
        if arr.size < 100:
            out[k + "_lo"] = math.nan
            out[k + "_hi"] = math.nan
        else:
            lo, hi = np.quantile(arr, [0.025, 0.975])
            out[k + "_lo"] = float(lo)
            out[k + "_hi"] = float(hi)
    return out


def bootstrap_linear(y: np.ndarray, tac: np.ndarray, cldn: np.ndarray, rng: np.random.Generator) -> dict:
    mask = np.isfinite(y) & np.isfinite(tac) & np.isfinite(cldn)
    y, tac, cldn = y[mask], tac[mask], cldn[mask]
    keys = ["beta_unadj", "beta_adj", "diff_beta", "indirect", "partial_spearman"]
    store = {k: [] for k in keys}
    failed = 0
    n = len(y)
    for _ in range(N_BOOT):
        idx = boot_indices(n, rng)
        fit = linear_pair(y[idx], tac[idx], cldn[idx])
        if fit is None:
            failed += 1
            continue
        for k in keys:
            store[k].append(fit[k])
    out = {"n_boot": N_BOOT, "n_failed": failed, "n": n}
    for k, vals in store.items():
        arr = np.array(vals, dtype=float)
        arr = arr[np.isfinite(arr)]
        out[k + "_boot_n"] = int(arr.size)
        if arr.size < 100:
            out[k + "_lo"] = math.nan
            out[k + "_hi"] = math.nan
        else:
            lo, hi = np.quantile(arr, [0.025, 0.975])
            out[k + "_lo"] = float(lo)
            out[k + "_hi"] = float(hi)
    return out


def adjusted_median_or(y: np.ndarray, tac: np.ndarray, cldn: np.ndarray) -> dict | None:
    high = (tac >= np.median(tac)).astype(float)
    zm = zscore(cldn)
    res = logit_fit(y, np.column_stack([np.ones(len(y)), high, zm]))
    if res is None:
        return None
    return {
        "OR": float(math.exp(res.params[1])),
        "OR_lo": float(math.exp(res.conf_int()[1, 0])),
        "OR_hi": float(math.exp(res.conf_int()[1, 1])),
        "p": float(res.pvalues[1]),
    }


def cd8_control(y: np.ndarray, cd8: np.ndarray) -> dict | None:
    mask = np.isfinite(y) & np.isfinite(cd8)
    y, cd8 = y[mask], cd8[mask]
    res = logit_fit(y, np.column_stack([np.ones(len(y)), zscore(cd8)]))
    if res is None:
        return {"n": int(len(y)), "OR": None, "OR_lo": None, "OR_hi": None, "p": None, "unstable": True}
    lo = float(math.exp(res.conf_int()[1, 0]))
    hi = float(math.exp(res.conf_int()[1, 1]))
    unstable = abs(float(res.params[1])) > 5 or (lo > 0 and hi / lo > 100)
    if unstable:
        return {"n": int(len(y)), "OR": None, "OR_lo": None, "OR_hi": None, "p": None, "unstable": True}
    return {
        "n": int(len(y)),
        "OR": float(math.exp(res.params[1])),
        "OR_lo": lo,
        "OR_hi": hi,
        "p": float(res.pvalues[1]),
        "unstable": False,
    }


def analyze_response(df: pd.DataFrame, rng: np.random.Generator) -> dict | None:
    sub = df.dropna(subset=["benefit", "TACSTD2", "CLDN4"]).copy()
    y = sub["benefit"].to_numpy(float)
    tac = sub["TACSTD2"].to_numpy(float)
    cldn = sub["CLDN4"].to_numpy(float)
    if y.sum() < 3 or (len(y) - y.sum()) < 3:
        return None
    point = prob_mediation(y, tac, cldn)
    if point is None:
        point = {"fit": "mle_failed", "n": int(len(y)), "n_benefit": int(y.sum())}
        boot = {}
    else:
        point["fit"] = "mle"
        boot = bootstrap_response(y, tac, cldn, rng)
    med = median_split_or(y, tac)
    med_c = median_split_or(y, cldn)
    adj_med = adjusted_median_or(y, tac, cldn)
    control = cd8_control(y, sub["CD8A"].to_numpy(float))
    return {
        "point": point,
        "bootstrap": boot,
        "median_TACSTD2": med,
        "median_CLDN4": med_c,
        "median_TACSTD2_adj_CLDN4": adj_med,
        "CD8A_perSD": control,
    }


def analyze_continuous(df: pd.DataFrame, outcome: str, rng: np.random.Generator) -> dict | None:
    if outcome not in df.columns:
        return None
    point = linear_pair(df[outcome].to_numpy(float), df["TACSTD2"].to_numpy(float), df["CLDN4"].to_numpy(float))
    if point is None:
        return None
    boot = bootstrap_linear(df[outcome].to_numpy(float), df["TACSTD2"].to_numpy(float), df["CLDN4"].to_numpy(float), rng)
    point["bootstrap"] = boot
    return point


def analyze_inflamed(df: pd.DataFrame, rng: np.random.Generator) -> dict | None:
    if "immune" not in df.columns:
        return None
    sub = df.dropna(subset=["immune", "TACSTD2", "CLDN4"]).copy()
    sub = sub[sub["immune"].isin(["desert", "excluded", "inflamed"])]
    if sub.empty:
        return None
    sub["inflamed"] = (sub["immune"] == "inflamed").astype(float)
    order = {"desert": 0, "excluded": 1, "inflamed": 2}
    sub["ordinal"] = sub["immune"].map(order).astype(float)
    y = sub["inflamed"].to_numpy(float)
    tac = sub["TACSTD2"].to_numpy(float)
    cldn = sub["CLDN4"].to_numpy(float)
    binary = prob_mediation(y, tac, cldn) if y.sum() >= 3 and (len(y) - y.sum()) >= 3 else None
    boot = bootstrap_response(y, tac, cldn, rng) if binary is not None else {}
    ordinal = linear_pair(sub["ordinal"].to_numpy(float), tac, cldn)
    ord_boot = bootstrap_linear(sub["ordinal"].to_numpy(float), tac, cldn, rng) if ordinal is not None else {}
    counts = {k: int((sub["immune"] == k).sum()) for k in ("desert", "excluded", "inflamed")}
    return {
        "n": int(len(sub)),
        "counts": counts,
        "inflamed_vs_not": binary,
        "inflamed_bootstrap": boot,
        "ordinal": ordinal,
        "ordinal_bootstrap": ord_boot,
    }


def imvigor_os(df: pd.DataFrame) -> dict:
    sub = df.dropna(subset=["os_time", "os_event", "TACSTD2", "CLDN4"]).copy()
    sub = sub[sub["os_time"] > 0]
    zt = zscore(sub["TACSTD2"].to_numpy(float))
    zm = zscore(sub["CLDN4"].to_numpy(float))
    csv_path = CACHE / "imvigor_os_input.csv"
    pd.DataFrame(
        {
            "time": sub["os_time"].to_numpy(float),
            "event": sub["os_event"].to_numpy(float),
            "zt": zt,
            "zm": zm,
        }
    ).to_csv(csv_path, index=False)
    r_path = CACHE / "imvigor_os.R"
    out_path = CACHE / "imvigor_os.json"
    r_path.write_text(
        f"""
library(survival)
d <- read.csv("{csv_path.as_posix()}")
fit1 <- coxph(Surv(time, event) ~ zt, data = d)
fit2 <- coxph(Surv(time, event) ~ zt + zm, data = d)
s1 <- summary(fit1)
s2 <- summary(fit2)
set.seed({SEED})
nboot <- {N_BOOT}
diff <- rep(NA_real_, nboot)
n <- nrow(d)
for (i in seq_len(nboot)) {{
  idx <- sample.int(n, n, replace = TRUE)
  b <- d[idx, ]
  if (sum(b$event) < 5) next
  f1 <- tryCatch(coxph(Surv(time, event) ~ zt, data = b), error = function(e) NULL)
  f2 <- tryCatch(coxph(Surv(time, event) ~ zt + zm, data = b), error = function(e) NULL)
  if (is.null(f1) || is.null(f2)) next
  diff[i] <- coef(f1)[["zt"]] - coef(f2)[["zt"]]
}}
diff <- diff[is.finite(diff)]
num <- function(x) formatC(as.numeric(x), digits = 8, format = "fg")
line <- paste0(
  "{{",
  '"n":', num(nrow(d)), ",",
  '"events":', num(sum(d$event)), ",",
  '"HR_unadj":', num(s1$coefficients["zt", "exp(coef)"]), ",",
  '"HR_unadj_lo":', num(s1$conf.int["zt", "lower .95"]), ",",
  '"HR_unadj_hi":', num(s1$conf.int["zt", "upper .95"]), ",",
  '"HR_unadj_p":', num(s1$coefficients["zt", "Pr(>|z|)"]), ",",
  '"HR_adj":', num(s2$coefficients["zt", "exp(coef)"]), ",",
  '"HR_adj_lo":', num(s2$conf.int["zt", "lower .95"]), ",",
  '"HR_adj_hi":', num(s2$conf.int["zt", "upper .95"]), ",",
  '"HR_adj_p":', num(s2$coefficients["zt", "Pr(>|z|)"]), ",",
  '"HR_CLDN4":', num(s2$coefficients["zm", "exp(coef)"]), ",",
  '"HR_CLDN4_p":', num(s2$coefficients["zm", "Pr(>|z|)"]), ",",
  '"diff_logHR":', num(coef(fit1)[["zt"]] - coef(fit2)[["zt"]]), ",",
  '"diff_lo":', num(quantile(diff, 0.025)), ",",
  '"diff_hi":', num(quantile(diff, 0.975)), ",",
  '"diff_boot_n":', length(diff),
  "}}"
)
writeLines(line, "{out_path.as_posix()}")
"""
    )
    import subprocess

    subprocess.check_call(["Rscript", str(r_path)])
    return json.loads(out_path.read_text())


def ci_crosses_zero(lo: float, hi: float) -> bool:
    return (not math.isfinite(lo)) or (not math.isfinite(hi)) or (lo <= 0 <= hi)


def attenuation_call(
    unadj: float,
    adj: float,
    diff_lo: float,
    diff_hi: float,
    unadj_lo: float,
    unadj_hi: float,
    adj_lo: float,
    adj_hi: float,
) -> str:
    total_null = ci_crosses_zero(unadj_lo, unadj_hi)
    change_null = ci_crosses_zero(diff_lo, diff_hi)
    adj_null = ci_crosses_zero(adj_lo, adj_hi)
    shrunk = math.isfinite(unadj) and math.isfinite(adj) and abs(adj) < abs(unadj)
    if total_null and change_null:
        return "no unadjusted association to attenuate"
    if (not total_null) and shrunk and (not change_null) and unadj * adj >= 0:
        if adj_null:
            return "attenuates; adjusted interval includes 0"
        return "partial attenuation; adjusted association remains"
    if shrunk and change_null:
        return "point estimate shrinks; interval for the change includes 0"
    if math.isfinite(adj) and math.isfinite(unadj) and abs(adj) > abs(unadj) and not change_null:
        return "does not attenuate; adjusted coefficient is farther from 0"
    if total_null:
        return "unadjusted interval includes 0"
    return "does not attenuate"


def proportion_mediated(total, total_lo, total_hi, indirect, indirect_lo, indirect_hi):
    if total in (0, None) or indirect in (None,):
        return None
    if not all(math.isfinite(v) for v in (total, total_lo, total_hi, indirect, indirect_lo, indirect_hi)):
        return None
    if ci_crosses_zero(total_lo, total_hi) or ci_crosses_zero(indirect_lo, indirect_hi):
        return None
    if total * indirect <= 0:
        return None
    return indirect / total


def round_or_none(x, nd=4):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return None
    return round(float(x), nd)


def response_row(cohort: str, endpoint: str, block: dict) -> dict:
    p = block["point"]
    b = block.get("bootstrap") or {}
    med = block["median_TACSTD2"]
    adj = block["median_TACSTD2_adj_CLDN4"]
    row = {
        "cohort": cohort,
        "endpoint": endpoint,
        "n": p.get("n"),
        "n_benefit": p.get("n_benefit"),
        "fit": p.get("fit"),
        "pearson_TACSTD2_CLDN4": p.get("pearson_r"),
        "vif": p.get("vif"),
        "a_CLDN4_on_TACSTD2": p.get("a_CLDN4_on_TACSTD2"),
        "a_p": p.get("a_p"),
        "OR_perSD_unadj": p.get("OR_unadj"),
        "OR_perSD_unadj_lo": p.get("OR_unadj_lo"),
        "OR_perSD_unadj_hi": p.get("OR_unadj_hi"),
        "OR_perSD_unadj_p": p.get("logOR_unadj_p"),
        "OR_perSD_adj": p.get("OR_adj"),
        "OR_perSD_adj_lo": p.get("OR_adj_lo"),
        "OR_perSD_adj_hi": p.get("OR_adj_hi"),
        "OR_perSD_adj_p": p.get("logOR_adj_p"),
        "OR_CLDN4_perSD_given_TACSTD2": p.get("OR_CLDN4"),
        "OR_CLDN4_p": p.get("OR_CLDN4_p"),
        "diff_logOR": p.get("diff_logOR"),
        "diff_logOR_lo": b.get("diff_logOR_lo"),
        "diff_logOR_hi": b.get("diff_logOR_hi"),
        "NIE": p.get("NIE"),
        "NIE_lo": b.get("NIE_lo"),
        "NIE_hi": b.get("NIE_hi"),
        "NDE": p.get("NDE"),
        "NDE_lo": b.get("NDE_lo"),
        "NDE_hi": b.get("NDE_hi"),
        "total_risk_diff": p.get("total_risk_diff"),
        "total_risk_diff_lo": b.get("total_risk_diff_lo"),
        "total_risk_diff_hi": b.get("total_risk_diff_hi"),
        "interaction_p": p.get("interaction_p"),
        "boot_failed": b.get("n_failed"),
        "median_OR": med.get("OR") if med else None,
        "median_OR_lo": med.get("OR_lo") if med else None,
        "median_OR_hi": med.get("OR_hi") if med else None,
        "median_fisher_p": med.get("fisher_p") if med else None,
        "median_zero_cell": med.get("zero_cell") if med else None,
        "median_a": med.get("a_resp_high") if med else None,
        "median_b": med.get("b_nonresp_high") if med else None,
        "median_c": med.get("c_resp_low") if med else None,
        "median_d": med.get("d_nonresp_low") if med else None,
        "median_adj_OR": None if adj is None else adj.get("OR"),
        "median_adj_OR_lo": None if adj is None else adj.get("OR_lo"),
        "median_adj_OR_hi": None if adj is None else adj.get("OR_hi"),
        "median_adj_p": None if adj is None else adj.get("p"),
        "CD8A_OR": None if block["CD8A_perSD"] is None else block["CD8A_perSD"]["OR"],
        "CD8A_OR_lo": None if block["CD8A_perSD"] is None else block["CD8A_perSD"]["OR_lo"],
        "CD8A_OR_hi": None if block["CD8A_perSD"] is None else block["CD8A_perSD"]["OR_hi"],
        "CD8A_p": None if block["CD8A_perSD"] is None else block["CD8A_perSD"]["p"],
        "CD8A_unstable": None if block["CD8A_perSD"] is None else block["CD8A_perSD"].get("unstable", False),
    }
    if p.get("fit") == "mle":
        row["attenuation"] = attenuation_call(
            p["logOR_unadj"],
            p["logOR_adj"],
            b.get("diff_logOR_lo", math.nan),
            b.get("diff_logOR_hi", math.nan),
            math.log(p["OR_unadj_lo"]),
            math.log(p["OR_unadj_hi"]),
            math.log(p["OR_adj_lo"]),
            math.log(p["OR_adj_hi"]),
        )
        row["proportion_mediated"] = proportion_mediated(
            p["total_risk_diff"],
            b.get("total_risk_diff_lo", math.nan),
            b.get("total_risk_diff_hi", math.nan),
            p["NIE"],
            b.get("NIE_lo", math.nan),
            b.get("NIE_hi", math.nan),
        )
    else:
        row["attenuation"] = "per-SD logistic MLE did not converge"
        row["proportion_mediated"] = None
    return row


def continuous_row(cohort: str, outcome: str, block: dict) -> dict:
    b = block["bootstrap"]
    # Wald-style interval from the bootstrap of the standardized coefficient.
    row = {
        "cohort": cohort,
        "outcome": outcome,
        "n": block["n"],
        "spearman": block["spearman"],
        "spearman_p": block["spearman_p"],
        "partial_spearman": block["partial_spearman"],
        "partial_spearman_p": block["partial_spearman_p"],
        "beta_unadj": block["beta_unadj"],
        "beta_unadj_lo": b.get("beta_unadj_lo"),
        "beta_unadj_hi": b.get("beta_unadj_hi"),
        "beta_unadj_p": block["beta_unadj_p"],
        "beta_adj": block["beta_adj"],
        "beta_adj_lo": b.get("beta_adj_lo"),
        "beta_adj_hi": b.get("beta_adj_hi"),
        "beta_adj_p": block["beta_adj_p"],
        "indirect": block["indirect"],
        "indirect_lo": b.get("indirect_lo"),
        "indirect_hi": b.get("indirect_hi"),
        "diff_beta": block["diff_beta"],
        "diff_beta_lo": b.get("diff_beta_lo"),
        "diff_beta_hi": b.get("diff_beta_hi"),
        "boot_failed": b.get("n_failed"),
    }
    row["attenuation"] = attenuation_call(
        block["beta_unadj"],
        block["beta_adj"],
        b.get("diff_beta_lo", math.nan),
        b.get("diff_beta_hi", math.nan),
        b.get("beta_unadj_lo", math.nan),
        b.get("beta_unadj_hi", math.nan),
        b.get("beta_adj_lo", math.nan),
        b.get("beta_adj_hi", math.nan),
    )
    row["proportion_mediated"] = proportion_mediated(
        block["beta_unadj"],
        b.get("beta_unadj_lo", math.nan),
        b.get("beta_unadj_hi", math.nan),
        block["indirect"],
        b.get("indirect_lo", math.nan),
        b.get("indirect_hi", math.nan),
    )
    return row


def forest_response(rows: list[dict]) -> None:
    use = [r for r in rows if r.get("fit") == "mle" and r.get("OR_perSD_unadj") is not None]
    if not use:
        return
    FIGS.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.2, 0.55 * len(use) + 1.4))
    for i, r in enumerate(use):
        y0 = i + 0.12
        y1 = i - 0.12
        ax.plot([r["OR_perSD_unadj_lo"], r["OR_perSD_unadj_hi"]], [y0, y0], color="#1f4e79", lw=1.5)
        ax.plot(r["OR_perSD_unadj"], y0, "o", color="#1f4e79", ms=5.5)
        ax.plot([r["OR_perSD_adj_lo"], r["OR_perSD_adj_hi"]], [y1, y1], color="#b86e00", lw=1.5)
        ax.plot(r["OR_perSD_adj"], y1, "o", color="#b86e00", ms=5.5)
    ax.axvline(1, color="#666666", lw=0.8)
    ax.set_yticks(range(len(use)))
    ax.set_yticklabels([f"{r['cohort']}  n={r['n']} ({r['n_benefit']} benefit)" for r in use])
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio per 1 SD TACSTD2  (benefit)")
    ax.set_title("Open ICI: TACSTD2 vs benefit, before and after CLDN4")
    ax.plot([], [], "o", color="#1f4e79", label="TACSTD2 only")
    ax.plot([], [], "o", color="#b86e00", label="TACSTD2 + CLDN4")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18)
    fig.savefig(FIGS / "fig_response_or.png", dpi=160)
    fig.savefig(FIGS / "fig_response_or.pdf")
    plt.close(fig)


def forest_continuous(rows: list[dict]) -> None:
    use = [r for r in rows if r.get("beta_unadj") is not None and r["outcome"] in {"CD8A", "Teff", "log2_TMB", "inflamed_ordinal"}]
    if not use:
        return
    FIGS.mkdir(parents=True, exist_ok=True)
    fig_h = 0.42 * len(use) + 1.3
    fig, ax = plt.subplots(figsize=(8.4, fig_h))
    for i, r in enumerate(use):
        y0 = i + 0.12
        y1 = i - 0.12
        if math.isfinite(r.get("beta_unadj_lo") or math.nan):
            ax.plot([r["beta_unadj_lo"], r["beta_unadj_hi"]], [y0, y0], color="#1f4e79", lw=1.4)
        ax.plot(r["beta_unadj"], y0, "o", color="#1f4e79", ms=5)
        if math.isfinite(r.get("beta_adj_lo") or math.nan):
            ax.plot([r["beta_adj_lo"], r["beta_adj_hi"]], [y1, y1], color="#b86e00", lw=1.4)
        ax.plot(r["beta_adj"], y1, "o", color="#b86e00", ms=5)
    ax.axvline(0, color="#666666", lw=0.8)
    ax.set_yticks(range(len(use)))
    ax.set_yticklabels([f"{r['cohort']}  {r['outcome']}  n={r['n']}" for r in use])
    ax.set_xlabel("Standardized beta, TACSTD2 vs outcome")
    ax.set_title("TACSTD2 vs immune or TMB, before and after CLDN4")
    ax.plot([], [], "o", color="#1f4e79", label="TACSTD2 only")
    ax.plot([], [], "o", color="#b86e00", label="TACSTD2 + CLDN4")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.16)
    fig.savefig(FIGS / "fig_tmb_immune.png", dpi=160)
    fig.savefig(FIGS / "fig_tmb_immune.pdf")
    plt.close(fig)


def dump_tsv(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def main() -> None:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    access = oak_poplar_status()
    absent_190266 = gse190266_tacstd2_absent()

    cohorts = [
        (
            "IMvigor210",
            "CR/PR vs SD/PD",
            "mUC, atezolizumab, single arm. binaryResponse, not durable clinical benefit. Not NSCLC and not a treatment-by-biomarker interaction.",
            extract_imvigor(),
        ),
        (
            "GSE218989",
            "GEO responder vs non-responder",
            "Lung, PD-1/PD-L1 inhibitor, n deposited as responder/non-responder. Pre versus post is not in the sample fields.",
            load_gse218989(),
        ),
        (
            "GSE190265",
            "PFS >= 6 (deposited units)",
            "Dijon anti-PD-1 monotherapy, tumor at diagnosis. Benefit is time >= 6. Early censoring (time < 6, event 0) is excluded. Not RECIST ORR.",
            load_gse190265(),
        ),
        (
            "GSE135222",
            "PFS >= 180 days",
            "NSCLC anti-PD-1/PD-L1. GEO pfs=1 is the progression event because those times are shorter. Benefit is PFS time >= 180 days. Not author-assigned DCB and not RECIST ORR.",
            load_gse135222(),
        ),
        (
            "GSE207422",
            "pathologic MPR vs NMPR",
            "Pretreatment bulk only. Neoadjuvant anti-PD-1 plus chemotherapy, not ICI monotherapy. MPR includes pCR.",
            load_gse207422(),
        ),
        (
            "GSE166449",
            "GEO title responder vs non-responder",
            "Pretreatment advanced lung. Response is the sample title. Expression is already on a log scale (CLDN4 max about 3) and was not logged again.",
            load_gse166449(),
        ),
        (
            "GSE126044",
            "GEO responder vs non-responder",
            "Pretreatment anti-PD-1 NSCLC. Counts converted to log2(CPM+1). FFPE samples were kept.",
            load_gse126044(),
        ),
    ]

    response_rows = []
    continuous_rows = []
    per_rows = []
    scale_checks = []
    for name, endpoint, note, df in cohorts:
        if "TACSTD2" not in df.columns or df["TACSTD2"].notna().sum() == 0 or df["CLDN4"].notna().sum() == 0:
            raise SystemExit(f"{name} is missing TACSTD2 or CLDN4 after the join")
        print(f"{name:12} joined {len(df):4}  labeled {int(df['benefit'].notna().sum())}  TACSTD2 median {df['TACSTD2'].median():.3f}")
        scale_checks.append(
            {
                "cohort": name,
                "n_joined": int(len(df)),
                "TACSTD2_min": float(df["TACSTD2"].min()),
                "TACSTD2_max": float(df["TACSTD2"].max()),
                "CLDN4_min": float(df["CLDN4"].min()),
                "CLDN4_max": float(df["CLDN4"].max()),
                "note": note,
                "endpoint": endpoint,
            }
        )
        # Scale guards. Already-log cohorts stay below ~15. Linear-then-log2 cohorts
        # of these epithelial genes land above that only if the raw matrix was log-scaled twice
        # or left linear. GSE166449 and GSE207422 are the already-log deposits.
        cmax = float(df["CLDN4"].max())
        if name in {"GSE166449", "GSE207422"} and cmax > 20:
            raise SystemExit(f"{name} CLDN4 max {cmax} is too high for an already-log matrix")
        if name not in {"GSE166449", "GSE207422"} and cmax > 20:
            raise SystemExit(f"{name} CLDN4 max {cmax} looks unlogged")
        block = analyze_response(df, rng)
        if block is not None:
            row = response_row(name, endpoint, block)
            row["note"] = note
            response_rows.append(row)
            print(
                f"  OR/SD {row['OR_perSD_unadj']} -> {row['OR_perSD_adj']}  "
                f"median OR {row['median_OR']}  {row['attenuation']}"
            )
        for outcome, label in (
            ("CD8A", "CD8A"),
            ("Teff", "Teff"),
        ):
            cont = analyze_continuous(df, outcome, rng)
            if cont is not None:
                crow = continuous_row(name, label, cont)
                continuous_rows.append(crow)
                print(f"  {label:8} rho {crow['spearman']:.3f} -> partial {crow['partial_spearman']:.3f}  {crow['attenuation']}")
        if "tmb" in df.columns:
            work = df.dropna(subset=["tmb"]).copy()
            work["log2_TMB"] = np.log2(work["tmb"].astype(float) + 1.0)
            cont = analyze_continuous(work, "log2_TMB", rng)
            if cont is not None:
                crow = continuous_row(name, "log2_TMB", cont)
                continuous_rows.append(crow)
                print(f"  log2_TMB rho {crow['spearman']:.3f} -> partial {crow['partial_spearman']:.3f}  {crow['attenuation']}")
            # Median TMB-high as a binary OR, cut = this cohort's median, not a clinical cutoff.
            y = (work["tmb"] >= work["tmb"].median()).to_numpy(float)
            tac = work["TACSTD2"].to_numpy(float)
            cldn = work["CLDN4"].to_numpy(float)
            point = prob_mediation(y, tac, cldn)
            if point is not None:
                boot = bootstrap_response(y, tac, cldn, rng)
                med = median_split_or(y, tac)
                block = {
                    "point": {**point, "fit": "mle"},
                    "bootstrap": boot,
                    "median_TACSTD2": med,
                    "median_CLDN4": median_split_or(y, cldn),
                    "median_TACSTD2_adj_CLDN4": adjusted_median_or(y, tac, cldn),
                    "CD8A_perSD": None,
                }
                row = response_row(name, "TMB >= cohort median", block)
                row["note"] = "Binary TMB uses the within-cohort median. It is not the 10 mut/Mb cutoff."
                response_rows.append(row)
        inflamed = analyze_inflamed(df, rng)
        if inflamed and inflamed["inflamed_vs_not"] is not None:
            point = inflamed["inflamed_vs_not"]
            sub = df.dropna(subset=["immune", "TACSTD2", "CLDN4"])
            sub = sub[sub["immune"].isin(["desert", "excluded", "inflamed"])]
            y = (sub["immune"] == "inflamed").to_numpy(float)
            tac = sub["TACSTD2"].to_numpy(float)
            cldn = sub["CLDN4"].to_numpy(float)
            block = {
                "point": {**point, "fit": "mle"},
                "bootstrap": inflamed["inflamed_bootstrap"],
                "median_TACSTD2": median_split_or(y, tac),
                "median_CLDN4": median_split_or(y, cldn),
                "median_TACSTD2_adj_CLDN4": adjusted_median_or(y, tac, cldn),
                "CD8A_perSD": None,
            }
            row = response_row(name, "immune phenotype inflamed vs not", block)
            row["note"] = "IMvigor210 immune phenotype. Inflamed versus desert or excluded. " + json.dumps(inflamed["counts"])
            response_rows.append(row)
        if inflamed and inflamed["ordinal"] is not None:
            ord_block = dict(inflamed["ordinal"])
            ord_block["bootstrap"] = inflamed["ordinal_bootstrap"]
            continuous_rows.append(continuous_row(name, "inflamed_ordinal", ord_block))

        keep_cols = [
            c
            for c in [
                "sample",
                "gsm",
                "benefit",
                "TACSTD2",
                "CLDN4",
                "CD8A",
                "Teff",
                "tmb",
                "neoantigen",
                "immune",
                "os_time",
                "os_event",
            ]
            if c in df.columns
        ]
        part = df[keep_cols].copy()
        part.insert(0, "cohort", name)
        per_rows.append(part)

    imvigor = cohorts[0][3]
    # Sensitivity: median computed on all 348 tumors, then applied to response-evaluable
    # patients. This is the split used by an earlier TACSTD2-only extraction.
    med_all = float(imvigor["TACSTD2"].median())
    r_all = imvigor.dropna(subset=["benefit"])
    high_all = r_all["TACSTD2"].to_numpy(float) >= med_all
    y_all = r_all["benefit"].to_numpy(float)
    a = int(np.sum(high_all & (y_all == 1)))
    b = int(np.sum(high_all & (y_all == 0)))
    c = int(np.sum(~high_all & (y_all == 1)))
    d = int(np.sum(~high_all & (y_all == 0)))
    median_all = woolf(a, b, c, d)
    median_all["median_on"] = "all 348 tumors, including response-missing"
    median_all["median"] = med_all
    os_res = imvigor_os(imvigor)
    # TMB as a positive control for response in IMvigor.
    resp = imvigor.dropna(subset=["benefit", "tmb"])
    tmb_or = cd8_control(resp["benefit"].to_numpy(float), np.log2(resp["tmb"].to_numpy(float) + 1.0))

    # Response ~ TACSTD2 + log2 TMB, then + CLDN4. One pre-specified covariate model.
    sub = imvigor.dropna(subset=["benefit", "TACSTD2", "CLDN4", "tmb"])
    y = sub["benefit"].to_numpy(float)
    zt = zscore(sub["TACSTD2"].to_numpy(float))
    zm = zscore(sub["CLDN4"].to_numpy(float))
    zb = zscore(np.log2(sub["tmb"].to_numpy(float) + 1.0))
    m1 = logit_fit(y, np.column_stack([np.ones(len(y)), zt, zb]))
    m2 = logit_fit(y, np.column_stack([np.ones(len(y)), zt, zb, zm]))
    tmb_adjusted = None
    if m1 is not None and m2 is not None:
        tmb_adjusted = {
            "n": int(len(y)),
            "n_benefit": int(y.sum()),
            "OR_TACSTD2_given_TMB": float(math.exp(m1.params[1])),
            "OR_TACSTD2_given_TMB_lo": float(math.exp(m1.conf_int()[1, 0])),
            "OR_TACSTD2_given_TMB_hi": float(math.exp(m1.conf_int()[1, 1])),
            "p_TACSTD2_given_TMB": float(m1.pvalues[1]),
            "OR_TMB_given_TACSTD2": float(math.exp(m1.params[2])),
            "p_TMB": float(m1.pvalues[2]),
            "OR_TACSTD2_given_TMB_CLDN4": float(math.exp(m2.params[1])),
            "OR_TACSTD2_given_TMB_CLDN4_lo": float(math.exp(m2.conf_int()[1, 0])),
            "OR_TACSTD2_given_TMB_CLDN4_hi": float(math.exp(m2.conf_int()[1, 1])),
            "p_TACSTD2_given_TMB_CLDN4": float(m2.pvalues[1]),
            "OR_CLDN4_in_that_model": float(math.exp(m2.params[3])),
            "p_CLDN4": float(m2.pvalues[3]),
        }

    payload = {
        "seed": SEED,
        "n_boot": N_BOOT,
        "oak_poplar": access,
        "GSE190266": absent_190266,
        "scale_checks": scale_checks,
        "imvigor_tmb_neoantigen_spearman": {
            "rho": imvigor.attrs.get("tmb_neo_rho"),
            "p": imvigor.attrs.get("tmb_neo_p"),
            "n": imvigor.attrs.get("tmb_neo_n"),
        },
        "imvigor_TMB_perSD_vs_response": tmb_or,
        "imvigor_TACSTD2_response_given_TMB": tmb_adjusted,
        "imvigor_os": os_res,
        "imvigor_median_OR_all_sample_cut": median_all,
        "response": response_rows,
        "continuous": continuous_rows,
    }
    (TABLES / "summary.json").write_text(json.dumps(payload, indent=2, default=float) + "\n")
    dump_tsv(TABLES / "response_or.tsv", response_rows)
    dump_tsv(TABLES / "continuous.tsv", continuous_rows)
    dump_tsv(TABLES / "scale_checks.tsv", scale_checks)
    pd.concat(per_rows, ignore_index=True).to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)
    forest_response([r for r in response_rows if r["endpoint"] != "TMB >= cohort median" and "inflamed" not in r["endpoint"]])
    forest_continuous(continuous_rows)
    print("wrote", TABLES)


if __name__ == "__main__":
    main()
