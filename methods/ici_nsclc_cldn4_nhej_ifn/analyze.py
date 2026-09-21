#!/usr/bin/env python3
"""Open ICI NSCLC bulk: CLDN4 vs response, plus NHEJ and IFN scores.

Patient is the unit. Cohorts are not pooled. POPLAR/OAK linear RNA-seq and
TCCIA linear expression stay closed; those rows are inventory only.

Scores, computed inside each cohort on log2 expression:
  CLDN4          single gene
  NHEJ           mean z of KEGG hsa03450 (Non-homologous end-joining)
  NHEJ core      same, without MRN (MRE11, RAD50), DNTT, FEN1
  IFN Teff       mean z of the POPLAR Teff/IFNg 8 genes
  IFN Hallmark   mean z of MSigDB Hallmark interferon gamma response (200)

Response contrast is Cliff's delta (benefit minus no benefit; positive means
higher in the benefit group) and a two-sided Mann-Whitney test.
"""

from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, pearsonr, rankdata, spearmanr

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "data" / "cache"
GENE_DIR = ROOT / "gene_sets"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"

KEGG_NHEJ = [
    "RAD50",
    "DNTT",
    "FEN1",
    "XRCC6",
    "POLL",
    "POLM",
    "LIG4",
    "MRE11",
    "PRKDC",
    "DCLRE1C",
    "XRCC4",
    "XRCC5",
    "NHEJ1",
]
# Core c-NHEJ. Drops MRN, TdT, and FEN1, which are not NHEJ-specific.
NHEJ_CORE = ["XRCC6", "XRCC5", "PRKDC", "LIG4", "XRCC4", "NHEJ1", "DCLRE1C", "POLL", "POLM"]
TEFF = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]

# GSE135222 is Ensembl. IDs checked against current Ensembl symbols.
# ENSG00000163251 is FZD5 in current Ensembl and is not used for PRKDC.
ENSEMBL = {
    "ENSG00000189143": "CLDN4",
    "ENSG00000153563": "CD8A",
    "ENSG00000145649": "GZMA",
    "ENSG00000100453": "GZMB",
    "ENSG00000111537": "IFNG",
    "ENSG00000163508": "EOMES",
    "ENSG00000138755": "CXCL9",
    "ENSG00000169245": "CXCL10",
    "ENSG00000073861": "TBX21",
    "ENSG00000113522": "RAD50",
    "ENSG00000107447": "DNTT",
    "ENSG00000168496": "FEN1",
    "ENSG00000196419": "XRCC6",
    "ENSG00000166169": "POLL",
    "ENSG00000122678": "POLM",
    "ENSG00000174405": "LIG4",
    "ENSG00000020922": "MRE11",
    "ENSG00000253729": "PRKDC",
    "ENSG00000152457": "DCLRE1C",
    "ENSG00000152422": "XRCC4",
    "ENSG00000079246": "XRCC5",
    "ENSG00000187736": "NHEJ1",
    "ENSG00000148773": "MKI67",
}

# Legacy symbols still used by some GEO matrices.
ALIASES = {
    "MRE11A": "MRE11",
    "WARS": "WARS1",
    "MARCH1": "MARCHF1",
    "RIGI": "DDX58",
}


def expand_ensembl() -> None:
    extra = json.loads((GENE_DIR / "hallmark_ensembl.json").read_text())
    for symbol, ensg in extra.items():
        ENSEMBL.setdefault(ensg, symbol)


def load_hallmark() -> list[str]:
    genes = [
        g.strip()
        for g in (GENE_DIR / "hallmark_ifng.txt").read_text().splitlines()
        if g.strip()
    ]
    if len(genes) != 200:
        raise SystemExit(f"Hallmark IFN-gamma list has {len(genes)} genes, expected 200")
    return genes


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
        row = {"title": title, "gsm": geos[i], "description": descs[0][i] if descs and descs[0] else ""}
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


def load_gene_by_sample(
    path: Path,
    needed: set[str],
    *,
    sep: str = "\t",
    skip_fields: int = 0,
    ensembl: bool = False,
    symbol_strip_underscore: bool = False,
) -> tuple[list[str], dict[str, np.ndarray]]:
    """Return sample ids and gene -> values (raw scale, duplicates averaged)."""
    found: dict[str, list[np.ndarray]] = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split(sep)
        samples = header[1 + skip_fields :]
        n = len(samples)
        for line in handle:
            parts = line.rstrip("\n").split(sep)
            if len(parts) < 2 + skip_fields:
                continue
            gene = parts[0].strip().strip('"')
            if ensembl:
                gene = ENSEMBL.get(gene.split(".")[0], "")
            elif symbol_strip_underscore:
                gene = gene.split("_")[0]
            gene = ALIASES.get(gene, gene)
            if gene not in needed:
                continue
            vals = np.array([_float(x) for x in parts[1 + skip_fields : 1 + skip_fields + n]], dtype=float)
            if vals.size != n:
                continue
            found.setdefault(gene, []).append(vals)
    data = {g: np.nanmean(np.vstack(vs), axis=0) for g, vs in found.items()}
    return samples, data


def load_transposed_semicolon(path: Path, needed: set[str]) -> tuple[list[str], dict[str, np.ndarray]]:
    with gzip.open(path, "rt", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split(";")
        index = {}
        for i, gene in enumerate(header):
            gene = ALIASES.get(gene.strip(), gene.strip())
            if gene in needed and gene not in index:
                index[gene] = i
        samples: list[str] = []
        cols = {g: [] for g in index}
        for line in handle:
            parts = line.rstrip("\n").split(";")
            if not parts or not parts[0]:
                continue
            samples.append(parts[0])
            for gene, i in index.items():
                cols[gene].append(_float(parts[i]) if i < len(parts) else math.nan)
    data = {g: np.array(vs, dtype=float) for g, vs in cols.items()}
    return samples, data


def library_sizes(path: Path) -> tuple[list[str], np.ndarray]:
    with gzip.open(path, "rt", errors="replace") as handle:
        samples = handle.readline().rstrip("\n").split("\t")[1:]
        libs = np.zeros(len(samples), dtype=float)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            vals = parts[1 : 1 + len(samples)]
            libs += np.array([_float(x) if x.strip() not in {"", "NA"} else 0.0 for x in vals], dtype=float)
    return samples, libs


def to_log(values: np.ndarray, *, already_log: bool, cpm_libs: np.ndarray | None) -> np.ndarray:
    if already_log:
        return values.astype(float)
    x = values.astype(float)
    if cpm_libs is not None:
        x = x / cpm_libs * 1e6
    return np.log2(x + 1.0)


def zmean(matrix_genes_by_samples: np.ndarray) -> np.ndarray | None:
    """matrix is genes x samples of log expression. Return mean z across genes."""
    if matrix_genes_by_samples.size == 0:
        return None
    kept = []
    for row in matrix_genes_by_samples:
        if np.sum(np.isfinite(row)) < 3:
            continue
        sd = np.nanstd(row, ddof=1)
        if not np.isfinite(sd) or sd == 0:
            continue
        mu = np.nanmean(row)
        kept.append((row - mu) / sd)
    if not kept:
        return None
    return np.nanmean(np.vstack(kept), axis=0)


def cliffs_delta(benefit: np.ndarray, other: np.ndarray) -> float:
    diff = benefit[:, None] - other[None, :]
    return float((np.sum(diff > 0) - np.sum(diff < 0)) / diff.size)


def bootstrap_delta(benefit: np.ndarray, other: np.ndarray, rng: np.random.Generator, nboot: int = 2000):
    stats = np.empty(nboot, dtype=float)
    nb, no = len(benefit), len(other)
    for i in range(nboot):
        stats[i] = cliffs_delta(benefit[rng.integers(0, nb, nb)], other[rng.integers(0, no, no)])
    lo, hi = np.quantile(stats, [0.025, 0.975])
    return float(lo), float(hi)


def mw_p(benefit: np.ndarray, other: np.ndarray) -> float:
    try:
        res = mannwhitneyu(benefit, other, alternative="two-sided", method="auto")
    except TypeError:
        res = mannwhitneyu(benefit, other, alternative="two-sided")
    return float(res.pvalue)


def spearman(x: np.ndarray, y: np.ndarray):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 5:
        return math.nan, math.nan, int(mask.sum())
    r, p = spearmanr(x[mask], y[mask])
    return float(r), float(p), int(mask.sum())


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray):
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if mask.sum() < 6:
        return math.nan, math.nan
    rx, ry, rz = rankdata(x[mask]), rankdata(y[mask]), rankdata(z[mask])

    def resid(a, b):
        design = np.column_stack([np.ones(len(b)), b])
        coef, *_ = np.linalg.lstsq(design, a, rcond=None)
        return a - design @ coef

    r, p = pearsonr(resid(rx, rz), resid(ry, rz))
    return float(r), float(p)


def benefit_label(text: str) -> int | None:
    t = text.lower().replace(" ", "").replace("_", "")
    if "non-responder" in t or "nonresponder" in t:
        return 0
    if "responder" in t:
        return 1
    return None


def build_expr(samples: list[str], data: dict[str, np.ndarray], *, already_log: bool, libs: np.ndarray | None):
    log = {}
    for gene, vals in data.items():
        use_libs = libs if libs is not None else None
        log[gene] = to_log(vals, already_log=already_log, cpm_libs=use_libs)
    by_sample = []
    for i, sid in enumerate(samples):
        by_sample.append({g: float(v[i]) for g, v in log.items()})
    return by_sample


def align(clinical: list[dict], samples: list[str], by_sample: list[dict], key_fn):
    index = {sid: i for i, sid in enumerate(samples)}
    rows = []
    missed = []
    for rec in clinical:
        key = key_fn(rec)
        if key not in index:
            missed.append(key)
            continue
        expr = by_sample[index[key]]
        rows.append({**rec, "expr_key": key, "expr": expr})
    return rows, missed


def gene_vector(rows: list[dict], gene: str) -> np.ndarray:
    return np.array([r["expr"].get(gene, math.nan) for r in rows], dtype=float)


NHEJ_CORE6 = {"XRCC5", "XRCC6", "PRKDC", "LIG4", "XRCC4", "NHEJ1"}


def score_vector(rows: list[dict], genes: list[str]) -> tuple[np.ndarray | None, list[str]]:
    present = []
    rows_m = []
    min_obs = max(5, int(0.8 * len(rows)))
    for gene in genes:
        vec = gene_vector(rows, gene)
        if np.sum(np.isfinite(vec)) >= min_obs:
            present.append(gene)
            rows_m.append(vec)
    if not rows_m:
        return None, present
    return zmean(np.vstack(rows_m)), present


def score_ok(name: str, present: list[str]) -> bool:
    got = set(present)
    if name == "NHEJ":
        # Ku70/80, DNA-PKcs, Lig4, XRCC4, XLF, plus at least 11/13 KEGG genes.
        return NHEJ_CORE6.issubset(got) and len(got) >= 11
    if name == "NHEJ_core":
        return NHEJ_CORE6.issubset(got) and len(got) >= 7
    if name == "IFN_Teff":
        return len(got) >= 7
    if name == "IFN_Hallmark":
        return len(got) >= 180
    return False


def contrast(rows, values, rng):
    b = np.array([values[i] for i, r in enumerate(rows) if r["benefit"] == 1], dtype=float)
    o = np.array([values[i] for i, r in enumerate(rows) if r["benefit"] == 0], dtype=float)
    b = b[np.isfinite(b)]
    o = o[np.isfinite(o)]
    if len(b) < 3 or len(o) < 3:
        return {
            "n_benefit": int(len(b)),
            "n_other": int(len(o)),
            "delta": math.nan,
            "delta_lo": math.nan,
            "delta_hi": math.nan,
            "p": math.nan,
            "median_benefit": float(np.median(b)) if len(b) else math.nan,
            "median_other": float(np.median(o)) if len(o) else math.nan,
        }
    lo, hi = bootstrap_delta(b, o, rng)
    return {
        "n_benefit": int(len(b)),
        "n_other": int(len(o)),
        "delta": cliffs_delta(b, o),
        "delta_lo": lo,
        "delta_hi": hi,
        "p": mw_p(b, o),
        "median_benefit": float(np.median(b)),
        "median_other": float(np.median(o)),
    }


def prepare_cohorts(hallmark: list[str]):
    needed = set(KEGG_NHEJ + NHEJ_CORE + TEFF + hallmark + ["CLDN4", "MKI67"])
    cohorts = []

    # GSE126044 counts -> log2(CPM+1)
    clin = series_samples(CACHE / "GSE126044_series_matrix.txt.gz")
    for rec in clin:
        rec["benefit"] = benefit_label(rec.get("patient response", ""))
        rec["patient"] = rec["title"]
    samples, libs = library_sizes(CACHE / "GSE126044_counts.txt.gz")
    samples2, data = load_gene_by_sample(CACHE / "GSE126044_counts.txt.gz", needed)
    if samples != samples2:
        raise SystemExit("GSE126044 sample order mismatch")
    by = build_expr(samples, data, already_log=False, libs=libs)
    rows, missed = align(clin, samples, by, lambda r: r["title"].replace("RNA-seq_", ""))
    cohorts.append(dict(
        cohort="GSE126044",
        rows=rows,
        missed=missed,
        role="response",
        n_note="pretreatment anti-PD-1 NSCLC; GEO responder vs non-responder; 5/16 FFPE kept",
        endpoint="GEO ORR responder vs non-responder",
        timepoint="pre-treatment",
        treatment="anti-PD-1",
        scale="log2(CPM+1)",
        benefit_name="responder",
        other_name="non-responder",
    ))

    # GSE135222 TPM, Ensembl. DCB = PFS time >= 180 days.
    clin = series_samples(CACHE / "GSE135222_series_matrix.txt.gz")
    for rec in clin:
        time = float(rec["pfs.time"])
        event = int(rec["progression-free survival (pfs)"])
        rec["pfs_time"] = time
        rec["pfs_event"] = event
        rec["benefit"] = 1 if time >= 180 else 0
        rec["patient"] = rec["title"]
    samples, data = load_gene_by_sample(CACHE / "GSE135222_exp.tsv.gz", needed, ensembl=True)
    by = build_expr(samples, data, already_log=False, libs=None)
    rows, missed = align(clin, samples, by, lambda r: r["title"].replace(" ", ""))
    cohorts.append(dict(
        cohort="GSE135222",
        rows=rows,
        missed=missed,
        role="response",
        n_note="advanced NSCLC anti-PD-1/PD-L1; DCB = PFS time >= 180 days; GEO pfs=1 is the progression event (all 27 evaluable, no censor before 180 days)",
        endpoint="DCB (PFS >= 180 days) vs not",
        timepoint="tumor RNA, timing not split in GEO",
        treatment="anti-PD-1/PD-L1",
        scale="log2(deposited expression + 1)",
        benefit_name="DCB",
        other_name="NDB",
    ))

    # GSE166449 TPM. Response is in the sample title, not a characteristic.
    clin = series_samples(CACHE / "GSE166449_series_matrix.txt.gz")
    for rec in clin:
        rec["benefit"] = benefit_label(rec["title"])
        rec["patient"] = rec["description"]
    samples, data = load_gene_by_sample(CACHE / "GSE166449_TPM.txt.gz", needed)
    by = build_expr(samples, data, already_log=False, libs=None)
    rows, missed = align(clin, samples, by, lambda r: r["description"])
    cohorts.append(dict(
        cohort="GSE166449",
        rows=rows,
        missed=missed,
        role="response",
        n_note="pretreatment advanced lung; response taken from GEO sample title (Responder / nonResponder)",
        endpoint="GEO title responder vs non-responder",
        timepoint="pre-treatment",
        treatment="immunotherapy (GEO title)",
        scale="log2(TPM+1)",
        benefit_name="responder",
        other_name="non-responder",
    ))

    # GSE190266 anti-PD-1 monotherapy. 6-month landmark.
    clin = series_samples(CACHE / "GSE190266_series_matrix.txt.gz")
    kept = []
    for rec in clin:
        time = float(rec["pfs_time (6 months)"])
        event = int(rec["pfs_evt (6 months)"])
        rec["pfs_time"] = time
        rec["pfs_event"] = event
        rec["patient"] = rec["title"]
        if event == 1:
            rec["benefit"] = 0
            kept.append(rec)
        elif event == 0 and time >= 5.99:
            rec["benefit"] = 1
            kept.append(rec)
        # event 0 and time < 6 is an early censor; dropped
    samples, data = load_transposed_semicolon(CACHE / "GSE190266_TPM.csv.gz", needed)
    by = build_expr(samples, data, already_log=False, libs=None)
    rows, missed = align(kept, samples, by, lambda r: r["title"])
    cohorts.append(dict(
        cohort="GSE190266",
        rows=rows,
        missed=missed,
        role="response",
        n_note="Dijon France4, anti-PD-1 monotherapy, tumor at diagnosis. Benefit = no PFS event and time censored at 6 months. No-benefit = PFS event before 6 months. One early censor (time 0.03, event 0) excluded. Deposited TPM is Excel-truncated at 16,384 columns (stops at MTMR14), so NHEJ and Hallmark IFN are not scored",
        endpoint="6-month PFS landmark (no event vs event)",
        timepoint="diagnosis, before anti-PD-1",
        treatment="anti-PD-1 monotherapy",
        scale="log2(TPM+1)",
        benefit_name="no event by 6 mo",
        other_name="event before 6 mo",
    ))

    # GSE207422 pretreatment bulk only. Neoadjuvant anti-PD-1 + chemotherapy.
    clin = series_samples(CACHE / "GSE207422_series_matrix.txt.gz")
    kept = []
    for rec in clin:
        if rec.get("description") != "Bulk RNAseq":
            continue
        if rec.get("sampling_time") != "Pre-treatment biopsy":
            continue
        path_resp = rec.get("pathologic_response", "")
        if path_resp.startswith("MPR"):
            rec["benefit"] = 1
        elif path_resp == "NMPR":
            rec["benefit"] = 0
        else:
            continue
        rec["patient"] = rec.get("patient", rec["title"])
        kept.append(rec)
    samples, data = load_gene_by_sample(CACHE / "GSE207422_bulk_log2TPM.txt.gz", needed)
    # Deposited values are already log2 TPM (filename and range).
    by = build_expr(samples, data, already_log=True, libs=None)
    rows, missed = align(kept, samples, by, lambda r: r["title"])
    cohorts.append(dict(
        cohort="GSE207422",
        rows=rows,
        missed=missed,
        role="response",
        n_note="pretreatment biopsy bulk only (BD Rhapsody single-cell libraries excluded). Neoadjuvant anti-PD-1 plus chemotherapy, not ICI monotherapy. Benefit = pathologic MPR including pCR",
        endpoint="pathologic MPR vs NMPR",
        timepoint="pre-treatment biopsy",
        treatment="anti-PD-1 + chemotherapy",
        scale="deposited log2 TPM",
        benefit_name="MPR",
        other_name="NMPR",
    ))

    # GSE218989 large open ICI lung bulk.
    clin = series_samples(CACHE / "GSE218989_series_matrix.txt.gz")
    for rec in clin:
        rec["benefit"] = benefit_label(rec.get("treatment outcome", ""))
        rec["patient"] = rec["title"]
    samples, data = load_gene_by_sample(CACHE / "GSE218989_TPM.txt.gz", needed)
    by = build_expr(samples, data, already_log=False, libs=None)
    rows, missed = align(clin, samples, by, lambda r: r["title"])
    cohorts.append(dict(
        cohort="GSE218989",
        rows=rows,
        missed=missed,
        role="response",
        n_note="GEO overall design: bulk RNA-seq of immunotherapy-treated lung cancer. Endpoint is the deposited Responder vs Non-responder label. Series title is the parent pan-cancer paper; the sample records say PD-1/PD-L1 inhibitor",
        endpoint="GEO responder vs non-responder",
        timepoint="not stated beyond treatment with PD-1/PD-L1 inhibitor",
        treatment="PD-1/PD-L1 inhibitor",
        scale="log2(TPM+1)",
        benefit_name="responder",
        other_name="non-responder",
    ))

    # Correlation-only: response column is not in the open matrix.
    for cohort, expr_name, skip, scale, note, timepoint, treatment in [
        (
            "GSE253564",
            "GSE253564_FPKM.txt.gz",
            1,
            "log2(FPKM+1)",
            "pretreatment durvalumab +/- SBRT trial RNA. GEO has arm, not MPR or recurrence, so this is not a response test",
            "pre-treatment",
            "durvalumab +/- SBRT (neoadjuvant)",
        ),
        (
            "GSE248378",
            "GSE248378_FPKM.txt.gz",
            0,
            "log2(FPKM+1)",
            "post-resection RNA from the same durvalumab trial. GEO has arm, not recurrence, so this is not a response test and is not a pre/post pair",
            "post-resection",
            "durvalumab +/- SBRT (neoadjuvant)",
        ),
    ]:
        clin = series_samples(CACHE / f"{cohort}_series_matrix.txt.gz")
        for rec in clin:
            rec["benefit"] = None
            rec["patient"] = rec["title"]
        samples, data = load_gene_by_sample(CACHE / expr_name, needed, skip_fields=skip)
        by = build_expr(samples, data, already_log=False, libs=None)
        rows, missed = align(clin, samples, by, lambda r: r["title"])
        cohorts.append(dict(
            cohort=cohort,
            rows=rows,
            missed=missed,
            role="correlation_only",
            n_note=note,
            endpoint="none in GEO",
            timepoint=timepoint,
            treatment=treatment,
            scale=scale,
            benefit_name="",
            other_name="",
        ))
    return cohorts


def summarize(cohort, hallmark, rng):
    rows = [r for r in cohort["rows"] if r.get("benefit") is not None or cohort["role"] != "response"]
    if cohort["role"] == "response":
        rows = [r for r in cohort["rows"] if r["benefit"] in (0, 1)]
    cldn4 = gene_vector(rows, "CLDN4")
    nhej, nhej_genes = score_vector(rows, KEGG_NHEJ)
    core, core_genes = score_vector(rows, NHEJ_CORE)
    teff, teff_genes = score_vector(rows, TEFF)
    hall, hall_genes = score_vector(rows, hallmark)
    if not score_ok("NHEJ", nhej_genes):
        nhej = None
    if not score_ok("NHEJ_core", core_genes):
        core = None
    if not score_ok("IFN_Teff", teff_genes):
        teff = None
    if not score_ok("IFN_Hallmark", hall_genes):
        hall = None
    scores = {
        "CLDN4": cldn4,
        "NHEJ": nhej,
        "NHEJ_core": core,
        "IFN_Teff": teff,
        "IFN_Hallmark": hall,
    }
    out = {
        "cohort": cohort["cohort"],
        "role": cohort["role"],
        "n": len(rows),
        "n_benefit": int(sum(r["benefit"] == 1 for r in rows)) if cohort["role"] == "response" else "",
        "n_other": int(sum(r["benefit"] == 0 for r in rows)) if cohort["role"] == "response" else "",
        "missed_keys": cohort["missed"],
        "endpoint": cohort["endpoint"],
        "timepoint": cohort["timepoint"],
        "treatment": cohort["treatment"],
        "scale": cohort["scale"],
        "note": cohort["n_note"],
        "benefit_name": cohort["benefit_name"],
        "other_name": cohort["other_name"],
        "genes": {
            "NHEJ": nhej_genes,
            "NHEJ_core": core_genes,
            "IFN_Teff": teff_genes,
            "IFN_Hallmark_n": len(hall_genes),
            "CLDN4_present": bool(np.sum(np.isfinite(cldn4)) == len(rows) and len(rows) > 0),
        },
    }
    if "CLDN4" not in rows[0]["expr"] or not np.any(np.isfinite(cldn4)):
        out["error"] = "CLDN4 absent"
        return out, []

    contrasts = {}
    if cohort["role"] == "response":
        for name, vec in scores.items():
            if vec is None:
                contrasts[name] = None
            else:
                contrasts[name] = contrast(rows, vec, rng)
    out["contrast"] = contrasts

    def sp(a, b):
        if a is None or b is None:
            return (math.nan, math.nan, 0)
        return spearman(a, b)

    out["rho"] = {
        "CLDN4_NHEJ": sp(cldn4, nhej),
        "CLDN4_NHEJ_core": sp(cldn4, core),
        "CLDN4_IFN_Teff": sp(cldn4, teff),
        "CLDN4_IFN_Hallmark": sp(cldn4, hall),
        "NHEJ_IFN_Teff": sp(nhej, teff),
        "NHEJ_IFN_Hallmark": sp(nhej, hall),
    }
    mki67 = gene_vector(rows, "MKI67")
    if nhej is not None and teff is not None:
        out["partial_CLDN4_NHEJ_given_Teff"] = partial_spearman(cldn4, nhej, teff)
    else:
        out["partial_CLDN4_NHEJ_given_Teff"] = (math.nan, math.nan)
    if nhej is not None and np.sum(np.isfinite(mki67)) >= 6:
        out["partial_CLDN4_NHEJ_given_MKI67"] = partial_spearman(cldn4, nhej, mki67)
        out["rho_CLDN4_MKI67"] = spearman(cldn4, mki67)
    else:
        out["partial_CLDN4_NHEJ_given_MKI67"] = (math.nan, math.nan)
        out["rho_CLDN4_MKI67"] = (math.nan, math.nan, 0)

    per = []
    for i, rec in enumerate(rows):
        per.append({
            "cohort": cohort["cohort"],
            "patient": rec.get("patient", rec["title"]),
            "gsm": rec["gsm"],
            "benefit": "" if rec["benefit"] is None else rec["benefit"],
            "CLDN4": cldn4[i],
            "NHEJ": "" if nhej is None else nhej[i],
            "NHEJ_core": "" if core is None else core[i],
            "IFN_Teff": "" if teff is None else teff[i],
            "IFN_Hallmark": "" if hall is None else hall[i],
        })
    return out, per


def fmt_p(p):
    if p is None or not isinstance(p, (int, float)) or not math.isfinite(p):
        return ""
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_d(x):
    if x is None or not isinstance(x, (int, float)) or not math.isfinite(x):
        return ""
    return f"{x:.3f}"


def write_tables(results, per_rows):
    TABLES.mkdir(parents=True, exist_ok=True)
    header = [
        "cohort", "role", "n", "n_benefit", "n_other", "endpoint", "timepoint", "treatment", "scale",
        "CLDN4_delta", "CLDN4_p", "NHEJ_delta", "NHEJ_p", "NHEJ_core_delta", "NHEJ_core_p",
        "IFN_Teff_delta", "IFN_Teff_p", "IFN_Hallmark_delta", "IFN_Hallmark_p",
        "rho_CLDN4_NHEJ", "p_CLDN4_NHEJ", "rho_CLDN4_IFN_Teff", "p_CLDN4_IFN_Teff",
        "rho_CLDN4_IFN_Hallmark", "p_CLDN4_IFN_Hallmark",
        "partial_rho_CLDN4_NHEJ_given_Teff", "partial_p_Teff",
        "partial_rho_CLDN4_NHEJ_given_MKI67", "partial_p_MKI67",
        "NHEJ_genes_found", "Teff_genes_found", "Hallmark_genes_found",
        "NHEJ_scored", "Hallmark_scored",
    ]
    lines = ["\t".join(header)]
    for res in results:
        if res.get("error"):
            continue
        c = res.get("contrast") or {}

        def grab(name, field):
            block = c.get(name) if c else None
            if not block:
                return ""
            return fmt_d(block[field]) if field != "p" else fmt_p(block["p"])

        rho = res["rho"]
        pr, pp = res["partial_CLDN4_NHEJ_given_Teff"]
        mr, mp = res.get("partial_CLDN4_NHEJ_given_MKI67", (math.nan, math.nan))
        row = [
            res["cohort"], res["role"], str(res["n"]), str(res["n_benefit"]), str(res["n_other"]),
            res["endpoint"], res["timepoint"], res["treatment"], res["scale"],
            grab("CLDN4", "delta"), grab("CLDN4", "p"),
            grab("NHEJ", "delta"), grab("NHEJ", "p"),
            grab("NHEJ_core", "delta"), grab("NHEJ_core", "p"),
            grab("IFN_Teff", "delta"), grab("IFN_Teff", "p"),
            grab("IFN_Hallmark", "delta"), grab("IFN_Hallmark", "p"),
            fmt_d(rho["CLDN4_NHEJ"][0]), fmt_p(rho["CLDN4_NHEJ"][1]),
            fmt_d(rho["CLDN4_IFN_Teff"][0]), fmt_p(rho["CLDN4_IFN_Teff"][1]),
            fmt_d(rho["CLDN4_IFN_Hallmark"][0]), fmt_p(rho["CLDN4_IFN_Hallmark"][1]),
            fmt_d(pr), fmt_p(pp),
            fmt_d(mr), fmt_p(mp),
            str(len(res["genes"]["NHEJ"])), str(len(res["genes"]["IFN_Teff"])), str(res["genes"]["IFN_Hallmark_n"]),
            "yes" if math.isfinite(rho["CLDN4_NHEJ"][0]) else "no",
            "yes" if math.isfinite(rho["CLDN4_IFN_Hallmark"][0]) else "no",
        ]
        lines.append("\t".join(row))
    (TABLES / "one_row.tsv").write_text("\n".join(lines) + "\n")

    inv = ["cohort\tstatus\twhy"]
    inv.append("POPLAR/OAK RNA-seq\tclosed\tEGAS00001005013 log2(TPM+1), counts, CPM, and clinical (EGAD00001008390/08391/08548/08549/08628-08631) are Genentech DAC EGAC00001002120. No open per-sample CLDN4, NHEJ, or IFN score joined to atezolizumab response")
    inv.append("POPLAR public IFN summary\tpublished summary only\tFehrenbacher Lancet 2016: Teff/IFNg (CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, TBX21) high subgroup, atezolizumab vs docetaxel OS HR 0.43 (95% CI 0.24-0.77). This is a treatment-effect subgroup, not a within-arm CLDN4 test, and it is not an NHEJ result")
    inv.append("Bessede 2024 OAK/POPLAR supplements\topen but not CLDN4\tAACR figshare 10.1158/1078-0432.c.7077754 reports TACSTD2, not CLDN4 or NHEJ, and has no per-sample expression")
    inv.append("TCCIA\tcircRNA deposit only\tZenodo 10.5281/zenodo.7969298 files are TCCIA_Ensemble_circRNAs.tsv.gz and TCCIA_4_methods_circRNAs.tsv.gz. No linear CLDN4 and no NHEJ/IFN score. OAK/POPLAR inside TCCIA use the same controlled EGA matrices")
    inv.append("GSE136961\tCLDN4 absent\tOncomine Immune Response panel, 393 genes. CLDN4 absent. KEGG NHEJ genes XRCC5/XRCC6/PRKDC/LIG4/NHEJ1/DCLRE1C absent. Not scored")
    inv.append("GSE285029\tresponse label absent\t234 ICI lung libraries. GEO characteristics have no response, PFS, or DCB field. Deposited matrix contains negative CLDN4 values, so it is not counts or TPM. Not scored")
    inv.append("GSE182328\tno ICI response label\tGEO label is Akkermansia detectable vs not. That is not ORR, DCB, PFS, or MPR. Not used as a response test")
    inv.append("GSE248249\tnot an open bulk matrix\tSupplementary file is RAW.tar. Not used")
    for res in results:
        inv.append(f"{res['cohort']}\tscored\t{res['note']}")
    (TABLES / "inventory.tsv").write_text("\n".join(inv) + "\n")

    per_header = ["cohort", "patient", "gsm", "benefit", "CLDN4", "NHEJ", "NHEJ_core", "IFN_Teff", "IFN_Hallmark"]
    plines = ["\t".join(per_header)]
    for rec in per_rows:
        plines.append("\t".join(str(rec[h]) for h in per_header))
    (TABLES / "per_sample.tsv").write_text("\n".join(plines) + "\n")

    def _jsonable(obj):
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_jsonable(v) for v in obj]
        if isinstance(obj, float):
            return obj if math.isfinite(obj) else None
        return obj

    (TABLES / "summary.json").write_text(json.dumps(_jsonable(results), indent=2) + "\n")


def forest_plot(results):
    FIGS.mkdir(parents=True, exist_ok=True)
    use = [r for r in results if r["role"] == "response" and not r.get("error")]
    use = sorted(use, key=lambda r: r["n"])
    metrics = [
        ("CLDN4", "CLDN4", "#1f4e79"),
        ("NHEJ", "NHEJ", "#b86e00"),
        ("IFN_Teff", "IFN Teff", "#1e7a46"),
        ("IFN_Hallmark", "IFN Hallmark", "#5b2c6f"),
    ]
    fig_h = 1.15 + 0.55 * len(use)
    fig, ax = plt.subplots(figsize=(8.4, fig_h))
    yticks = []
    ylabels = []
    for i, res in enumerate(use):
        ylabels.append(f"{res['cohort']}  n={res['n']} ({res['n_benefit']}/{res['n_other']})")
        yticks.append(i)
        for j, (name, label, color) in enumerate(metrics):
            block = (res.get("contrast") or {}).get(name)
            if not block or not math.isfinite(block["delta"]):
                continue
            y = i + (j - 1.5) * 0.16
            ax.plot([block["delta_lo"], block["delta_hi"]], [y, y], color=color, lw=1.4, solid_capstyle="round")
            ax.plot(block["delta"], y, "o", color=color, ms=5.5, label=label if i == 0 else None)
    ax.axvline(0, color="#666666", lw=0.8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("Cliff's δ  (higher in benefit →)")
    ax.set_xlim(-1.05, 1.05)
    ax.set_title("Open ICI NSCLC bulk: score versus benefit")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.22)
    fig.savefig(FIGS / "fig_forest_cliffs.png", dpi=160)
    fig.savefig(FIGS / "fig_forest_cliffs.pdf")
    plt.close(fig)


def scatter_large(per_rows):
    rows = [r for r in per_rows if r["cohort"] == "GSE218989"]
    if len(rows) < 20:
        return
    benefit = np.array([int(r["benefit"]) for r in rows])
    cldn4 = np.array([float(r["CLDN4"]) for r in rows])
    nhej = np.array([float(r["NHEJ"]) for r in rows])
    teff = np.array([float(r["IFN_Teff"]) for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0), sharey=False)
    for ax, x, xlab in ((axes[0], nhej, "NHEJ score (KEGG mean z)"), (axes[1], teff, "IFN Teff score (mean z)")):
        for flag, color, label in ((0, "#b03a2e", "non-responder"), (1, "#1f4e79", "responder")):
            m = benefit == flag
            ax.scatter(x[m], cldn4[m], s=12, alpha=0.55, c=color, label=label, linewidths=0)
        ax.set_xlabel(xlab)
        ax.set_ylabel("CLDN4 log2(TPM+1)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[1].legend(frameon=False, loc="best")
    fig.suptitle("GSE218989  n=%d" % len(rows))
    fig.tight_layout()
    fig.savefig(FIGS / "fig_gse218989_cldn4.png", dpi=160)
    fig.savefig(FIGS / "fig_gse218989_cldn4.pdf")
    plt.close(fig)


GEO = {
    "GSE126044_counts.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
    "GSE126044_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz",
    "GSE135222_exp.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
    "GSE135222_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz",
    "GSE166449_TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz",
    "GSE166449_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/matrix/GSE166449_series_matrix.txt.gz",
    "GSE190266_TPM.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190266/suppl/GSE190266_TPM_France4.csv.gz",
    "GSE190266_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190266/matrix/GSE190266_series_matrix.txt.gz",
    "GSE207422_bulk_log2TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
    "GSE207422_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/matrix/GSE207422_series_matrix.txt.gz",
    "GSE218989_TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/suppl/GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz",
    "GSE218989_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/matrix/GSE218989_series_matrix.txt.gz",
    "GSE248378_FPKM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/suppl/GSE248378_Durva_Post_FPKMs.txt.gz",
    "GSE248378_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/matrix/GSE248378_series_matrix.txt.gz",
    "GSE253564_FPKM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",
    "GSE253564_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/matrix/GSE253564_series_matrix.txt.gz",
}


def ensure_cache() -> None:
    import urllib.request

    CACHE.mkdir(parents=True, exist_ok=True)
    for name, url in GEO.items():
        dest = CACHE / name
        if dest.exists() and dest.stat().st_size > 0:
            continue
        print("download", name)
        urllib.request.urlretrieve(url, dest)


def main():
    ensure_cache()
    rng = np.random.default_rng(1)
    expand_ensembl()
    hallmark = load_hallmark()
    cohorts = prepare_cohorts(hallmark)
    results = []
    per_rows = []
    for cohort in cohorts:
        if cohort["missed"]:
            print(f"WARNING {cohort['cohort']} missed {len(cohort['missed'])} keys, e.g. {cohort['missed'][:5]}")
        res, per = summarize(cohort, hallmark, rng)
        results.append(res)
        per_rows.extend(per)
        c = res.get("contrast") or {}
        cld = (c.get("CLDN4") or {})
        nh = (c.get("NHEJ") or {})
        te = (c.get("IFN_Teff") or {})
        print(
            f"{res['cohort']:12} n={res['n']:4} role={res['role']:16} "
            f"CLDN4 δ={cld.get('delta', float('nan'))!s:>8} p={cld.get('p', float('nan'))!s:>10} "
            f"NHEJ δ={nh.get('delta', float('nan'))!s:>8} p={nh.get('p', float('nan'))!s:>10} "
            f"Teff δ={te.get('delta', float('nan'))!s:>8} "
            f"ρ(CLDN4,NHEJ)={res.get('rho', {}).get('CLDN4_NHEJ', ('',))[0]!s:>8} "
            f"ρ(CLDN4,Teff)={res.get('rho', {}).get('CLDN4_IFN_Teff', ('',))[0]!s:>8}"
        )
        if not res["genes"]["CLDN4_present"]:
            print("  CLDN4 missing")
        print(
            f"  genes NHEJ {len(res['genes']['NHEJ'])}/13 core {len(res['genes']['NHEJ_core'])}/9 "
            f"Teff {len(res['genes']['IFN_Teff'])}/8 Hallmark {res['genes']['IFN_Hallmark_n']}/200"
        )
    write_tables(results, per_rows)
    forest_plot(results)
    scatter_large(per_rows)
    print("wrote", TABLES / "one_row.tsv")


if __name__ == "__main__":
    main()
