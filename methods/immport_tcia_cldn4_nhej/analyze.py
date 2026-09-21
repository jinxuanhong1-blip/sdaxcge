#!/usr/bin/env python3
"""ImmPort + TCIA: open NSCLC immunotherapy transcriptomes for CLDN4 vs NHEJ / IFN / APM.

The immunotherapy gate is applied before any correlation. ImmPort shared-study
metadata and TCIA collection pages are queried live. A cohort is scored only
when an open gene-expression matrix exists. Surgical matrices linked from TCIA
are scored and labeled NOT_ICI. They are not an immunotherapy result.

Gene sets are fixed below. A gene enters a score only when it is non-missing
in at least 80% of samples. NA is not filled with zero. RNA-seq values are
linear and are log2(x+1). The Rosetta/Merck matrix is already log2 intensity.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import re
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
CACHE = Path("/tmp/immport_tcia_cldn4_nhej")
TAB = HERE / "tables"
FIG = HERE / "figures"
for d in (TAB, FIG, CACHE):
    d.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0 (research; CLDN4 public-data audit)"}

# Pre-specified programs. IFN has no HLA / TAP / PSMB so it does not collapse into APM.
NHEJ = ["XRCC6", "XRCC5", "PRKDC", "LIG4", "XRCC4", "NHEJ1", "DCLRE1C", "PAXX"]
IFN = [
    "IFNG",
    "STAT1",
    "IRF1",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "IDO1",
    "GBP1",
    "GBP4",
    "GBP5",
    "CXCR6",
    "IL2RG",
]
APM = [
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "TAP1",
    "TAP2",
    "TAPBP",
    "PSMB8",
    "PSMB9",
    "PSMB10",
    "NLRC5",
    "ERAP1",
]
KRT = ["KRT8", "KRT18", "KRT19"]
PROGRAMS = {"NHEJ": NHEJ, "IFN": IFN, "APM": APM}
ALIAS_TO_GENE = {"C9ORF142": "PAXX", "XLF": "NHEJ1"}

MIN_GENE_FRAC = 0.80
MIN_GENES = 4
MIN_SAMPLE_FRAC = 0.50

RNA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103584/suppl/"
    "GSE103584_R01_NSCLC_RNAseq.txt.gz"
)
CLIN_URL = (
    "https://www.cancerimagingarchive.net/wp-content/uploads/"
    "NSCLCR01Radiogenomic_DATA_LABELS_2018-05-22_1500-shifted.csv"
)
ARRAY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58661/matrix/"
    "GSE58661_series_matrix.txt.gz"
)
GPL_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10379/annot/"
    "GPL10379.annot.gz"
)

ICI_TOKENS = (
    "pembro",
    "nivol",
    "atezo",
    "durva",
    "ipilimumab",
    "avelumab",
    "cemiplimab",
    "sintilimab",
    "tislelizumab",
    "toripalimab",
    "camrelizumab",
    "immuno",
    "pd-1",
    "pd1",
    "checkpoint",
    "opdivo",
    "keytruda",
    "tecentriq",
    "imfinzi",
)


def fetch(url: str, dest: Path | None = None, timeout: int = 180) -> bytes:
    if dest is not None and dest.exists() and dest.stat().st_size > 1000:
        return dest.read_bytes()
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    if dest is not None:
        dest.write_bytes(data)
    return data


def fetch_text(url: str, timeout: int = 90) -> tuple[int, str]:
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        return e.code, body


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_rho(r: float) -> str:
    if r != r:
        return "NA"
    return f"{r:+.3f}"


def spearman(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return {"n": n, "rho": float("nan"), "p": float("nan")}
    r, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(r), "p": float(p)}


def partial_spearman(x, y, z) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 5:
        return {"n": n, "rho": float("nan"), "p": float("nan")}
    rx, ry, rz = stats.rankdata(x[m]), stats.rankdata(y[m]), stats.rankdata(z[m])
    Z = np.column_stack([np.ones(n), rz])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    ex, ey = rx - Z @ bx, ry - Z @ by
    if np.std(ex) == 0 or np.std(ey) == 0:
        return {"n": n, "rho": float("nan"), "p": float("nan")}
    r = float(np.corrcoef(ex, ey)[0, 1])
    if abs(r) >= 1:
        return {"n": n, "rho": r, "p": 0.0}
    df = n - 3
    t = r * math.sqrt(df / (1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return {"n": n, "rho": r, "p": p}


def bh(pvals: list[float]) -> list[float]:
    m = len(pvals)
    if m == 0:
        return []
    order = np.argsort([1 if (p != p) else p for p in pvals])
    q = np.ones(m)
    running = 1.0
    for rank_from_end, idx in enumerate(order[::-1]):
        rank = m - rank_from_end
        p = pvals[idx]
        if p != p:
            q[idx] = float("nan")
            continue
        running = min(running, p * m / rank)
        q[idx] = running
    return [float(v) for v in q]


def json_ready(obj):
    if isinstance(obj, dict):
        return {k: json_ready(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_ready(v) for v in obj]
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return None if not math.isfinite(v) else v
    if isinstance(obj, (np.integer,)):
        return int(obj)
    return obj


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for row in rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(
                {
                    k: ("" if isinstance(v, float) and not math.isfinite(v) else v)
                    for k, v in row.items()
                }
            )


def mean_z_score(
    log_genes: dict[str, np.ndarray],
    genes: list[str],
    n_samples: int,
    min_genes: int = MIN_GENES,
) -> dict:
    """log_genes values are already on the analysis scale. NaN = missing."""
    coverage = []
    eligible = []
    for g in genes:
        arr = log_genes.get(g)
        if arr is None:
            coverage.append(
                {"gene": g, "n_finite": 0, "frac": 0.0, "eligible": False, "present": False}
            )
            continue
        n_fin = int(np.isfinite(arr).sum())
        frac = n_fin / n_samples
        ok = frac >= MIN_GENE_FRAC
        coverage.append(
            {
                "gene": g,
                "n_finite": n_fin,
                "frac": frac,
                "eligible": ok,
                "present": True,
            }
        )
        if ok:
            eligible.append(g)
    out = {
        "coverage": coverage,
        "eligible": eligible,
        "score": np.full(n_samples, np.nan),
        "scored": False,
        "reason": "",
    }
    if len(eligible) < min_genes:
        out["reason"] = f"eligible {len(eligible)} < {min_genes}"
        return out
    zrows = []
    for g in eligible:
        arr = log_genes[g].astype(float)
        mu = np.nanmean(arr)
        sd = np.nanstd(arr)
        if not np.isfinite(sd) or sd == 0:
            continue
        z = (arr - mu) / sd
        zrows.append(z)
    if len(zrows) < min_genes:
        out["reason"] = f"nonzero-sd genes {len(zrows)} < {min_genes}"
        out["eligible"] = eligible
        return out
    Z = np.vstack(zrows)
    n_ok = np.isfinite(Z).sum(axis=0)
    need = max(min_genes, math.ceil(MIN_SAMPLE_FRAC * Z.shape[0]))
    with np.errstate(all="ignore"):
        mu = np.nanmean(Z, axis=0)
    mu[n_ok < need] = np.nan
    out["score"] = mu
    out["scored"] = True
    out["eligible"] = eligible
    out["n_z"] = int(Z.shape[0])
    out["min_genes_per_sample"] = int(need)
    return out


def load_rnaseq() -> tuple[list[str], dict[str, np.ndarray]]:
    path = CACHE / "GSE103584_R01_NSCLC_RNAseq.txt.gz"
    fetch(RNA_URL, path)
    want = set(NHEJ + IFN + APM + KRT + ["CLDN4"])
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        found: dict[str, np.ndarray] = {}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            gene = parts[0].strip().upper()
            if gene not in want:
                continue
            vals = []
            for x in parts[1:]:
                if x in ("", "NA", "NaN", "null", "None"):
                    vals.append(np.nan)
                else:
                    vals.append(float(x))
            arr = np.asarray(vals, float)
            # Linear FPKM-like matrix (HLA-A reaches 1e5). Do not treat NA as zero.
            found[gene] = np.log2(arr + 1.0)
    return samples, found


def load_clinical(samples: list[str]) -> dict:
    path = CACHE / "NSCLCR01_clinical.csv"
    fetch(CLIN_URL, path)
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [c.strip() for c in reader.fieldnames]
        rows = list(reader)
    by_id = {r["Case ID"].strip(): r for r in rows}
    blob = path.read_text(errors="replace").lower()
    token_hits = {tok: blob.count(tok) for tok in ICI_TOKENS}
    # "surgery" appears in the column name only.
    matched = []
    missing = []
    for s in samples:
        row = by_id.get(s)
        if row is None:
            missing.append(s)
        else:
            matched.append(row)
    def count_yes(col: str) -> int:
        return sum(1 for r in matched if r.get(col, "").strip().lower() == "yes")

    hist = {}
    for r in matched:
        h = r.get("Histology", "").strip() or "blank"
        hist[h] = hist.get(h, 0) + 1
    surgery_days = 0
    for r in matched:
        v = r.get("Days between CT and surgery", "").strip()
        if v and v.lower() not in {"na", "not collected", ""}:
            surgery_days += 1
    return {
        "n_clinical_rows": len(rows),
        "n_rna": len(samples),
        "n_overlap": len(matched),
        "n_rna_missing_clinical": len(missing),
        "histology": hist,
        "adjuvant_yes": count_yes("Adjuvant Treatment"),
        "chemotherapy_yes": count_yes("Chemotherapy"),
        "radiation_yes": count_yes("Radiation"),
        "surgery_interval_recorded": surgery_days,
        "ici_token_hits": token_hits,
        "ici_token_sum": int(sum(token_hits.values())),
        "histology_by_sample": [
            by_id.get(s, {}).get("Histology", "").strip() if s in by_id else "" for s in samples
        ],
    }


def load_gpl_probes() -> dict[str, list[str]]:
    """Exact single-symbol probes on GPL10379 (same HuRSTA array as GPL15048)."""
    path = CACHE / "GPL10379.annot.gz"
    fetch(GPL_URL, path)
    want = {g.upper() for g in (NHEJ + IFN + APM + KRT + ["CLDN4"])}
    want |= set(ALIAS_TO_GENE)
    probes: dict[str, list[str]] = {g: [] for g in (NHEJ + IFN + APM + KRT + ["CLDN4"])}
    with gzip.open(path, "rt", errors="replace") as f:
        started = False
        for line in f:
            if line.startswith("!platform_table_begin"):
                started = True
                continue
            if not started or line.startswith("ID\t") or line.startswith("!"):
                continue
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            syms = [s.strip() for s in parts[2].split("///") if s.strip()]
            if len(syms) != 1:
                continue
            key = syms[0].upper()
            gene = ALIAS_TO_GENE.get(key, key)
            if gene in probes and key in want:
                probes[gene].append(parts[0])
    return probes


def load_array(probes: dict[str, list[str]]) -> tuple[list[str], dict[str, np.ndarray], dict]:
    path = CACHE / "GSE58661_series_matrix.txt.gz"
    fetch(ARRAY_URL, path)
    wanted = {pid: gene for gene, ids in probes.items() for pid in ids}
    raw: dict[str, np.ndarray] = {}
    samples: list[str] = []
    hist: list[str] = []
    with gzip.open(path, "rt", errors="replace") as f:
        in_table = False
        for line in f:
            if line.startswith("!Sample_characteristics_ch1") and "histology:" in line:
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                hist = [v.split(":", 1)[-1].strip() for v in vals]
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if parts[0] == "ID_REF":
                samples = parts[1:]
                continue
            if parts[0] not in wanted:
                continue
            vals = []
            for x in parts[1:]:
                if x in ("", "NA", "null", "None"):
                    vals.append(np.nan)
                else:
                    vals.append(float(x))
            raw[parts[0]] = np.asarray(vals, float)
    chosen = {}
    probe_meta = {}
    for gene, ids in probes.items():
        present = [i for i in ids if i in raw]
        probe_meta[gene] = {"n_exact_probes": len(ids), "n_in_matrix": len(present), "probe": ""}
        if not present:
            continue
        # Highest median among exact-symbol probes.
        best = max(present, key=lambda i: np.nanmedian(raw[i]))
        chosen[gene] = raw[best]
        probe_meta[gene]["probe"] = best
        med = float(np.nanmedian(raw[best]))
        probe_meta[gene]["median"] = med
    # Rosetta series-matrix values for control and gene probes sit near 2–15 (log2).
    cldn = chosen.get("CLDN4")
    scale = "as_deposited"
    if cldn is not None and np.nanpercentile(cldn, 95) > 30:
        chosen = {g: np.log2(v + 1.0) for g, v in chosen.items()}
        scale = "log2(x+1)"
    return samples, chosen, {"probe_meta": probe_meta, "scale": scale, "histology": hist}


def score_cohort(name: str, samples: list[str], genes: dict[str, np.ndarray], extra: dict) -> dict:
    n = len(samples)
    cldn = genes.get("CLDN4", np.full(n, np.nan))
    programs = {}
    for label, glist in PROGRAMS.items():
        programs[label] = mean_z_score(genes, glist, n, min_genes=MIN_GENES)
    # Epithelial covariate is a 3-gene set. Two detected keratins are enough.
    programs["KRT"] = mean_z_score(genes, KRT, n, min_genes=2)
    tests = []
    for label in ("NHEJ", "IFN", "APM"):
        block = programs[label]
        raw = spearman(cldn, block["score"]) if block["scored"] else {
            "n": 0, "rho": float("nan"), "p": float("nan")
        }
        part = (
            partial_spearman(cldn, block["score"], programs["KRT"]["score"])
            if block["scored"] and programs["KRT"]["scored"]
            else {"n": 0, "rho": float("nan"), "p": float("nan")}
        )
        tests.append({"cohort": name, "program": label, "raw": raw, "partial_krt": part, "block": block})
    return {
        "name": name,
        "n": n,
        "cldn4_n": int(np.isfinite(cldn).sum()),
        "cldn4_median": float(np.nanmedian(cldn)) if np.isfinite(cldn).any() else float("nan"),
        "programs": programs,
        "tests": tests,
        "extra": extra,
        "cldn4": cldn,
        "samples": samples,
    }


def immport_catalog() -> tuple[list[dict], dict]:
    studies: dict[str, dict] = {}

    def add_search(params: dict, query_name: str) -> None:
        q = urllib.parse.urlencode(params)
        url = f"https://www.immport.org/data/query/api/search/study?{q}"
        status, text = fetch_text(url)
        if status != 200:
            return
        data = json.loads(text)
        for hit in data.get("hits", {}).get("hits", []):
            src = hit.get("_source", {})
            acc = src.get("study_accession")
            if not acc:
                continue
            rec = studies.setdefault(
                acc,
                {
                    "study_accession": acc,
                    "brief_title": src.get("brief_title") or "",
                    "condition_or_disease": "|".join(src.get("condition_or_disease") or []),
                    "research_focus": "|".join(src.get("research_focus") or []),
                    "queries": [],
                },
            )
            if query_name not in rec["queries"]:
                rec["queries"].append(query_name)

    add_search(
        {
            "researchFocus": "Oncology",
            "pageSize": 100,
            "sourceFields": "study_accession,brief_title,condition_or_disease,research_focus",
        },
        "researchFocus=Oncology",
    )
    for term in (
        "NSCLC",
        "non-small cell lung",
        "lung cancer",
        "pembrolizumab",
        "nivolumab",
        "atezolizumab",
        "durvalumab",
        "anti-PD-1",
        "checkpoint inhibitor",
    ):
        add_search(
            {
                "term": term,
                "pageSize": 40,
                "sourceFields": "study_accession,brief_title,condition_or_disease,research_focus",
            },
            f"term={term}",
        )
    # File listing is not public without a token.
    status, _ = fetch_text("https://www.immport.org/data/query/api/study/filePath/SDY901")
    notes = {
        "SDY901": "Mouse EGFR lung tumors (PD-1 pathway paper). Not human NSCLC on ICI. Files not public (401).",
        "SDY1079": "Mouse metastatic lung cancer model. Not human ICI RNA.",
        "SDY1298": "Mouse EZH2 lung cancer model. Not human ICI RNA.",
        "SDY1064": "Mouse LUAD metastasis model. Not human ICI RNA.",
        "SDY3403": "CyTOF of PC9 cells after osimertinib. Not ICI RNA-seq.",
        "SDY2310": "Flow cytometry of lung tumors. Not a transcriptome.",
        "SDY1733": "HCC, not NSCLC.",
        "SDY2295": "Glioblastoma, not NSCLC.",
        "SDY3726": "Liver cancer (nivolumab search hit). Not NSCLC.",
        "SDY3359": "Mouse MC38 colorectal. Not NSCLC.",
    }
    rows = []
    for acc in sorted(studies):
        rec = studies[acc]
        title = rec["brief_title"]
        blob = f"{title} {rec['condition_or_disease']} {rec['research_focus']}".lower()
        lungish = any(k in blob for k in ("lung", "nsclc"))
        iciish = any(k in blob for k in ("pd-1", "pd1", "pembro", "nivol", "atezo", "durva", "checkpoint", "immunotherapy"))
        if lungish and "covid" in blob:
            gate = "COVID or SARS-CoV-2 study that mentions lung. Not an NSCLC immunotherapy transcriptome."
        elif acc in notes:
            gate = notes[acc]
        elif lungish:
            gate = "Lung mentioned. Title is not a human NSCLC ICI tumor transcriptome."
        elif iciish:
            gate = "Checkpoint word in a non-lung study. Not scored."
        else:
            gate = "Cataloged. Not NSCLC immunotherapy RNA."
        rec["lungish"] = lungish
        rec["iciish"] = iciish
        rec["gate"] = gate
        rec["score"] = "not_scored"
        rows.append(
            {
                "study_accession": acc,
                "brief_title": title,
                "condition_or_disease": rec["condition_or_disease"],
                "research_focus": rec["research_focus"],
                "queries": ";".join(rec["queries"]),
                "gate": gate,
                "score": "not_scored",
            }
        )
    meta = {"n_studies": len(rows), "filePath_SDY901_http": status}
    return rows, meta


def tcia_catalog() -> list[dict]:
    status, text = fetch_text(
        "https://services.cancerimagingarchive.net/nbia-api/services/v1/getCollectionValues"
    )
    names = []
    if status == 200:
        data = json.loads(text)
        for d in data:
            names.append(d.get("Collection") or d.get("collection") or "")
    curated = {
        "Anti-PD-1_Lung": (
            "Human lung, anti-PD-1, n=46. Modalities PT/CT/SC only. No expression file.",
            "ICI_IMAGING_NO_RNA",
        ),
        "S0819": (
            "Phase III carboplatin/paclitaxel ± cetuximab. Imaging. Cetuximab is not PD-1/PD-L1.",
            "NOT_ICI_IMAGING",
        ),
        "CMB-LCA": (
            "Cancer Moonshot Biobank lung imaging plus clinical. No RNA matrix on the collection page.",
            "IMAGING_NO_RNA",
        ),
        "CPTAC-LUAD": (
            "Resection proteogenomics cohort. External genomics. Not ICI. Not scored.",
            "NOT_ICI",
        ),
        "CPTAC-LSCC": (
            "Resection proteogenomics cohort. External genomics. Not ICI. Not scored.",
            "NOT_ICI",
        ),
        "TCGA-LUAD": (
            "TCGA imaging collection. Treatment-naive surgical TCGA, not an ICI transcriptome.",
            "NOT_ICI",
        ),
        "TCGA-LUSC": (
            "TCGA imaging collection. Not an ICI transcriptome.",
            "NOT_ICI",
        ),
        "NSCLC Radiogenomics": (
            "Open RNA-seq GSE103584, n=130, surgically excised tumors (2008–2014). Scored as NOT_ICI.",
            "SCORED_NOT_ICI",
        ),
        "NSCLC-Radiomics-Genomics": (
            "Open microarray GSE58661 Lung3, n=89, surgery, not ICI. Deposited IDs are merck-EST probes. CLDN4 RefSeq probe is absent, so programs were not scored.",
            "NOT_ICI_CLDN4_ABSENT",
        ),
        "NSCLC-Radiomics": ("Imaging and outcomes. No gene-expression matrix.", "IMAGING_NO_RNA"),
        "NSCLC-Radiomics-Interobserver1": ("Contouring subset. No gene expression.", "IMAGING_NO_RNA"),
        "ACRIN-NSCLC-FDG-PET": ("PET imaging trial. No transcriptome.", "IMAGING_NO_RNA"),
        "4D-Lung": ("4D CT. No transcriptome.", "IMAGING_NO_RNA"),
        "Lung-PET-CT-Dx": ("PET/CT diagnosis. No transcriptome.", "IMAGING_NO_RNA"),
        "LungCT-Diagnosis": ("CT diagnosis. No transcriptome.", "IMAGING_NO_RNA"),
        "Lung-Fused-CT-Pathology": ("CT/pathology fusion. No transcriptome.", "IMAGING_NO_RNA"),
        "QIN LUNG CT": ("Quantitative imaging. No transcriptome.", "IMAGING_NO_RNA"),
        "RIDER Lung CT": ("Repeat CT. No transcriptome.", "IMAGING_NO_RNA"),
        "RIDER Lung PET-CT": ("Repeat PET/CT. No transcriptome.", "IMAGING_NO_RNA"),
        "SPIE-AAPM Lung CT Challenge": ("CT challenge. No transcriptome.", "IMAGING_NO_RNA"),
        "LIDC-IDRI": ("Lung nodule CT. No transcriptome.", "IMAGING_NO_RNA"),
        "NLST": ("Screening CT. No transcriptome.", "IMAGING_NO_RNA"),
    }
    # Confirm the ICI collection page has no expression spreadsheet.
    _, html = fetch_text("https://www.cancerimagingarchive.net/collection/anti-pd-1_lung/")
    hrefs = re.findall(r'href=["\']([^"\']+)["\']', html, flags=re.I)
    expr_hrefs = [
        h
        for h in hrefs
        if re.search(r"\.(csv|tsv|xls|xlsx|txt\.gz)$", h, flags=re.I)
        and not re.search(r"nbia-digest|manifest|\.tcia$", h, flags=re.I)
    ]
    rows = []
    seen = set()
    for name in sorted(set(names) | set(curated)):
        low = name.lower()
        if name not in curated and not any(
            k in low for k in ("lung", "nsclc", "pd-1", "pd1", "ici", "immuno")
        ):
            continue
        note, gate = curated.get(name, ("Lung-named collection was not in the curated gate.", "UNREVIEWED"))
        if name == "Anti-PD-1_Lung":
            note += f" Expression-like hrefs on the collection page: {len(expr_hrefs)}."
        rows.append(
            {
                "collection": name,
                "in_live_api": name in set(names),
                "gate": gate,
                "note": note,
            }
        )
        seen.add(name)
    rows.append(
        {
            "collection": "_api",
            "in_live_api": status == 200,
            "gate": "META",
            "note": f"getCollectionValues HTTP {status}; n_collections={len(names)}; expression_hrefs_anti_pd1={expr_hrefs}",
        }
    )
    return rows


def cohort_rows(result: dict, ici_label: str, verdict: str) -> list[dict]:
    rows = []
    for t in result["tests"]:
        raw, part = t["raw"], t["partial_krt"]
        block = t["block"]
        rows.append(
            {
                "cohort": result["name"],
                "ici": ici_label,
                "verdict": verdict,
                "n_samples": result["n"],
                "cldn4_finite": result["cldn4_n"],
                "program": t["program"],
                "genes_eligible": ",".join(block["eligible"]),
                "n_eligible": len(block["eligible"]),
                "n_in_set": len(PROGRAMS[t["program"]]),
                "score_reason": block.get("reason", ""),
                "n": raw["n"],
                "rho": raw["rho"],
                "p": raw["p"],
                "partial_krt_n": part["n"],
                "partial_krt_rho": part["rho"],
                "partial_krt_p": part["p"],
            }
        )
    return rows


def coverage_rows(cohort: str, result: dict, probe_meta: dict | None = None) -> list[dict]:
    rows = []
    for label, block in result["programs"].items():
        for cov in block["coverage"]:
            meta = (probe_meta or {}).get(cov["gene"], {})
            rows.append(
                {
                    "cohort": cohort,
                    "program": label,
                    "gene": cov["gene"],
                    "present": cov["present"],
                    "n_finite": cov["n_finite"],
                    "n_samples": result["n"],
                    "frac": f"{cov['frac']:.3f}",
                    "eligible_ge_80pct": cov["eligible"],
                    "n_exact_probes": meta.get("n_exact_probes", 1 if cov["present"] else 0),
                    "probe": meta.get("probe", "symbol" if cov["present"] else ""),
                }
            )
    return rows


def plot_scatter(result: dict, path: Path) -> None:
    hist = result["extra"].get("histology_by_sample") or [""] * result["n"]
    colors = []
    for h in hist:
        hl = h.lower()
        if "squamous" in hl:
            colors.append("#d95f02")
        elif "adeno" in hl:
            colors.append("#1b9e77")
        else:
            colors.append("#7570b3")
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.8), constrained_layout=True)
    cldn = result["cldn4"]
    for ax, label in zip(axes, ("NHEJ", "IFN", "APM")):
        y = result["programs"][label]["score"]
        ax.scatter(cldn, y, c=colors, s=18, alpha=0.85, linewidths=0)
        raw = next(t["raw"] for t in result["tests"] if t["program"] == label)
        ax.set_title(f"{label}  ρ={fmt_rho(raw['rho'])}  p={fmt_p(raw['p'])}  n={raw['n']}")
        ax.set_xlabel("CLDN4 log2(x+1)")
        ax.set_ylabel(f"{label} mean z")
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="#1b9e77", markersize=6, label="Adenocarcinoma"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="#d95f02", markersize=6, label="Squamous"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="#7570b3", markersize=6, label="NSCLC NOS"),
    ]
    axes[2].legend(handles=handles, frameon=False, fontsize=8, loc="best")
    fig.suptitle("GSE103584 surgical NSCLC RNA-seq — not immunotherapy")
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_forest(rows: list[dict], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    labels = []
    ys = []
    for i, row in enumerate(rows):
        if row["rho"] != row["rho"]:
            continue
        labels.append(f"{row['cohort']} {row['program']} (n={row['n']})")
        ys.append(row["rho"])
        ax.plot(row["rho"], len(labels) - 1, "o", color="#333333")
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Spearman ρ  CLDN4 vs program")
    ax.set_title("Open TCIA-linked NSCLC expression — surgical, not ICI")
    ax.set_xlim(-1, 1)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> None:
    print("ImmPort catalog", flush=True)
    imm_rows, imm_meta = immport_catalog()
    write_tsv(TAB / "immport_studies.tsv", imm_rows)

    print("TCIA catalog", flush=True)
    tcia_rows = tcia_catalog()
    write_tsv(TAB / "tcia_collections.tsv", tcia_rows)

    print("GSE103584", flush=True)
    samples, rna_genes = load_rnaseq()
    clinical = load_clinical(samples)
    rna = score_cohort(
        "GSE103584",
        samples,
        rna_genes,
        {"histology_by_sample": clinical.pop("histology_by_sample"), "clinical": clinical},
    )

    print("GSE58661", flush=True)
    probes = load_gpl_probes()
    arr_samples, arr_genes, arr_meta = load_array(probes)
    array = score_cohort("GSE58661", arr_samples, arr_genes, {"histology": arr_meta["histology"], "scale": arr_meta["scale"]})

    array_verdict = "NOT_ICI_SURGERY" if array["cldn4_n"] >= 4 else "NOT_ICI_CLDN4_ABSENT"
    primary = []
    for result, ici, verdict in (
        (rna, "no", "NOT_ICI_SURGERY"),
        (array, "no", array_verdict),
    ):
        primary.extend(cohort_rows(result, ici, verdict))
    gene_rows = []
    for gene in NHEJ + IFN + APM + KRT:
        arr = rna_genes.get(gene)
        if arr is None:
            sp = {"n": 0, "rho": float("nan"), "p": float("nan")}
        else:
            sp = spearman(rna["cldn4"], arr)
        program = (
            "NHEJ" if gene in NHEJ else "IFN" if gene in IFN else "APM" if gene in APM else "KRT"
        )
        gene_rows.append(
            {
                "cohort": "GSE103584",
                "program": program,
                "gene": gene,
                "role": "covariate" if program == "KRT" else "descriptive_gene",
                "n": sp["n"],
                "rho": sp["rho"],
                "p": sp["p"],
            }
        )
    gene_idx = [
        i
        for i, row in enumerate(gene_rows)
        if row["program"] != "KRT" and row["p"] == row["p"] and row["n"] >= 100
    ]
    gene_q = bh([gene_rows[i]["p"] for i in gene_idx])
    for row in gene_rows:
        row["q_bh_descriptive"] = ""
    for i, q in zip(gene_idx, gene_q):
        gene_rows[i]["q_bh_descriptive"] = q
    write_tsv(TAB / "gene_spearman.tsv", gene_rows)
    finite_idx = [i for i, row in enumerate(primary) if row["p"] == row["p"]]
    qvals = bh([primary[i]["p"] for i in finite_idx])
    for row in primary:
        row["q_bh_primary"] = ""
    for i, q in zip(finite_idx, qvals):
        primary[i]["q_bh_primary"] = q

    # Immunotherapy gate row: nothing to score.
    gate_row = {
        "cohort": "ImmPort+TCIA ICI",
        "ici": "yes_sought",
        "verdict": "NO_OPEN_TRANSCRIPTOME",
        "n_samples": 0,
        "cldn4_finite": 0,
        "program": "NHEJ|IFN|APM",
        "genes_eligible": "",
        "n_eligible": 0,
        "n_in_set": "",
        "score_reason": "No open human NSCLC immunotherapy expression matrix in ImmPort or TCIA.",
        "n": 0,
        "rho": "",
        "p": "",
        "partial_krt_n": "",
        "partial_krt_rho": "",
        "partial_krt_p": "",
        "q_bh_primary": "",
    }
    write_tsv(TAB / "spearman.tsv", primary)
    write_tsv(TAB / "one_row.tsv", [gate_row] + primary)

    cov = coverage_rows("GSE103584", rna)
    cov += coverage_rows("GSE58661", array, arr_meta["probe_meta"])
    write_tsv(TAB / "gene_coverage.tsv", cov)

    plot_scatter(rna, FIG / "gse103584_cldn4_vs_programs.png")
    plot_forest(primary, FIG / "forest_not_ici.png")

    def pack_test(t):
        return {
            "program": t["program"],
            "eligible": t["block"]["eligible"],
            "reason": t["block"].get("reason", ""),
            "raw": t["raw"],
            "partial_krt": t["partial_krt"],
        }

    summary = {
        "immunotherapy_open_matrices": 0,
        "immport": imm_meta,
        "immport_lung_or_ici_titles": [
            {"id": r["study_accession"], "title": r["brief_title"], "gate": r["gate"]}
            for r in imm_rows
            if any(
                k in (r["brief_title"] + r["condition_or_disease"]).lower()
                for k in ("lung", "nsclc", "pd-1", "pembro", "nivol", "atezo", "durva", "checkpoint")
            )
        ],
        "gse103584": {
            "scale": "log2(x+1) of deposited linear values; NA left missing",
            "clinical": clinical,
            "cldn4_finite": rna["cldn4_n"],
            "cldn4_median_log2": rna["cldn4_median"],
            "krt_eligible": rna["programs"]["KRT"]["eligible"],
            "tests": [pack_test(t) for t in rna["tests"]],
        },
        "gse58661": {
            "scale": arr_meta["scale"],
            "n": array["n"],
            "cldn4_finite": array["cldn4_n"],
            "cldn4_median": array["cldn4_median"],
            "krt_eligible": array["programs"]["KRT"]["eligible"],
            "tests": [pack_test(t) for t in array["tests"]],
        },
        "primary_bh": [
            {"cohort": primary[i]["cohort"], "program": primary[i]["program"], "q": primary[i]["q_bh_primary"]}
            for i in finite_idx
        ],
        "rules": {
            "min_gene_frac": MIN_GENE_FRAC,
            "min_genes": MIN_GENES,
            "min_sample_frac": MIN_SAMPLE_FRAC,
            "na_imputed_as_zero": False,
            "multi_gene_probes": "excluded",
            "ici_response_tested": False,
        },
    }
    summary["gse58661"]["verdict"] = array_verdict
    summary["gse58661"]["cldn4_probes_in_matrix"] = arr_meta["probe_meta"].get("CLDN4", {})
    (TAB / "summary.json").write_text(json.dumps(json_ready(summary), indent=2) + "\n")
    print(json.dumps(json_ready(summary), indent=2))


if __name__ == "__main__":
    main()
