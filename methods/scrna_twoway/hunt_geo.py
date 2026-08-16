#!/usr/bin/env python3
"""Hunt public lung scRNA series that carry BOTH timepoint and response labels.

A series is two-way only if the public record has:
  - a timepoint axis (pre / post / on-treatment / treatment-naive vs treated)
  - a response axis (MPR / NMPR / pCR / RECIST / R vs NR)
  - a processed matrix that includes epithelial / malignant cells (not T-cell-only)

Restricted archives (EGA, dbGaP, GSA HRA without a public processed matrix) are
recorded as not public.
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("results/scrna_twoway")
OUT.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "jinxuanhong1@gmail.com"

QUERIES = [
    (
        "lung_scrna_io",
        '("non-small cell lung" OR NSCLC OR LUAD OR LUSC OR "lung adenocarcinoma" OR "lung squamous" OR "lung cancer") AND (scRNA-seq OR "single-cell RNA" OR "single cell RNA" OR "single-cell transcriptome") AND (immunotherapy OR "PD-1" OR "PD-L1" OR pembrolizumab OR nivolumab OR neoadjuvant OR MPR OR RECIST)',
    ),
    (
        "lung_scrna_mpr",
        '("lung" OR NSCLC) AND (scRNA-seq OR "single-cell") AND (MPR OR "major pathologic" OR "pathologic response" OR pCR)',
    ),
]

# Manually curated near-miss / known series (always scored even if eSearch misses them).
SEED = [
    "GSE207422",
    "GSE291670",
    "GSE205335",
    "GSE243013",
    "GSE146100",
    "GSE179994",
    "GSE253013",
    "GSE248378",
    "GSE171145",
    "GSE148071",
    "GSE131907",
    "GSE127465",
    "GSE123814",
    "GSE176021",
    "GSE200981",
    "GSE154826",
    "GSE162498",
    "GSE189357",
    "GSE217557",
    "GSE221647",
    "GSE235916",
    "GSE271689",
    "GSE283829",
    "GSE285029",
]

TIME_RE = re.compile(
    r"\b(pre[- ]?treatment|post[- ]?treatment|treatment[- ]?naive|baseline|on[- ]?treatment|"
    r"pre[- ]?therapy|post[- ]?therapy|before treatment|after treatment|paired|timepoint|"
    r"pre[- ]?ICI|post[- ]?ICI|TN\b|untreated)\b",
    re.I,
)
RESP_RE = re.compile(
    r"\b(MPR|NMPR|pCR|RECIST|responder|non[- ]?responder|pathologic response|"
    r"major pathologic|partial response|stable disease|progressive disease|"
    r"\bPR\b|\bSD\b|\bPD\b|DCB|NDB)\b",
    re.I,
)
SCRNA_RE = re.compile(r"single[- ]cell|scRNA|10[xX]|BD Rhapsody|Drop[- ]seq", re.I)
LUNG_RE = re.compile(r"lung|NSCLC|LUAD|LUSC|adenocarcinoma", re.I)
T_ONLY_RE = re.compile(r"T cell only|T-cell only|all\.Tcell|scTCR|CD3\+ T", re.I)


def get(url: str, retries: int = 4) -> bytes:
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "scrna-twoway-hunt/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception as e:
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"GET failed {url}: {last}")


def esearch(term: str, retmax: int = 80) -> list[str]:
    q = urllib.parse.urlencode(
        {
            "db": "gds",
            "term": term + " AND gse[ETYP] AND Homo sapiens[ORGN]",
            "retmax": str(retmax),
            "retmode": "json",
            "email": EMAIL,
            "tool": "scrna_twoway",
        }
    )
    data = json.loads(get(f"{EUTILS}/esearch.fcgi?{q}"))
    return data.get("esearchresult", {}).get("idlist", [])


def esummary(ids: list[str]) -> list[dict]:
    out = []
    for i in range(0, len(ids), 20):
        chunk = ids[i : i + 20]
        q = urllib.parse.urlencode(
            {
                "db": "gds",
                "id": ",".join(chunk),
                "retmode": "json",
                "email": EMAIL,
                "tool": "scrna_twoway",
            }
        )
        data = json.loads(get(f"{EUTILS}/esummary.fcgi?{q}"))
        result = data.get("result", {})
        for uid in result.get("uids", []):
            out.append(result[uid])
        time.sleep(0.34)
    return out


def acc_from_summary(rec: dict) -> str | None:
    acc = rec.get("accession") or rec.get("Accession")
    if acc and str(acc).startswith("GSE"):
        return acc
    title = rec.get("title") or ""
    m = re.search(r"(GSE\d+)", title)
    return m.group(1) if m else None


def fetch_geo_soft_brief(acc: str) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief"
    try:
        return get(url).decode("utf-8", "replace")
    except Exception as e:
        return f"FETCH_FAIL {e}"


def classify(acc: str, title: str, summary: str, soft: str) -> dict:
    blob = "\n".join([acc, title, summary, soft])
    has_time = bool(TIME_RE.search(blob))
    has_resp = bool(RESP_RE.search(blob))
    is_scrna = bool(SCRNA_RE.search(blob))
    is_lung = bool(LUNG_RE.search(blob))
    t_only = bool(T_ONLY_RE.search(blob))
    public_matrix = bool(
        re.search(r"Supplementary file|processed data are available|UMI_matrix|rds\.gz|h5ad", blob, re.I)
    )
    # Hard-coded public facts for well-known series (override noisy text matches).
    overrides = {
        "GSE207422": {
            "has_timepoint": True,
            "has_response": True,
            "malignant_matrix": True,
            "t_cell_only": False,
            "public_processed": True,
            "twoway": True,
            "note": "3 pre-biopsy + 12 post-surgery; MPR/NMPR/pCR/NE and RECIST on the same GEO xlsx. Pre patients are not paired to the 12 post patients.",
        },
        "GSE291670": {
            "has_timepoint": False,
            "has_response": True,
            "malignant_matrix": True,
            "t_cell_only": False,
            "public_processed": True,
            "twoway": False,
            "note": "GEO deposit is 6 post-treatment NSCLC (3 MPR + 3 Non-MPR). The 3 treatment-naive samples used in the paper are HRA001033 (GSA), not in this public series.",
        },
        "GSE205335": {
            "has_timepoint": False,
            "has_response": True,
            "malignant_matrix": True,
            "t_cell_only": False,
            "public_processed": True,
            "twoway": False,
            "note": "RECIST on 14 core pre-ICI samples / 11 patients. Post-ICI / neoadjuvant samples were excluded from the response core. Raw in EGA.",
        },
        "GSE243013": {
            "has_timepoint": False,
            "has_response": True,
            "malignant_matrix": False,
            "t_cell_only": False,
            "public_processed": False,
            "twoway": False,
            "note": "234 post-neoadjuvant NSCLC; response exists in the paper. Public GEO record is metadata-oriented; no usable public UMI + response table for malignant cells at hunt time.",
        },
        "GSE146100": {
            "has_timepoint": False,
            "has_response": True,
            "malignant_matrix": True,
            "t_cell_only": False,
            "public_processed": True,
            "twoway": False,
            "note": "One patient, three post-pembrolizumab nodules with different radiographic responses. No pre-treatment scRNA.",
        },
        "GSE179994": {
            "has_timepoint": True,
            "has_response": True,
            "malignant_matrix": False,
            "t_cell_only": True,
            "public_processed": True,
            "twoway": False,
            "note": "Paired pre/post + good/poor response, but deposited matrix is T cells only. No malignant TACSTD2/CLDN4.",
        },
        "GSE253013": {
            "has_timepoint": False,
            "has_response": False,
            "malignant_matrix": True,
            "t_cell_only": False,
            "public_processed": True,
            "twoway": False,
            "note": "Treatment-naive LUAD atlas. No public MPR/RECIST.",
        },
        "GSE229353": {
            "has_timepoint": False,
            "has_response": True,
            "malignant_matrix": False,
            "t_cell_only": False,
            "public_processed": True,
            "twoway": False,
            "note": "scRNA is post-surgery CD45+ immune cells (n=7). Paired pre/post in the paper is FFPE IHC, not scRNA. No malignant compartment.",
        },
        "GSE131907": {
            "has_timepoint": False,
            "has_response": False,
            "malignant_matrix": True,
            "t_cell_only": False,
            "public_processed": True,
            "twoway": False,
            "note": "Early LUAD atlas; not an ICI response series.",
        },
    }
    rec = {
        "accession": acc,
        "title": title[:240],
        "is_lung": is_lung,
        "is_scrna": is_scrna,
        "has_timepoint": has_time,
        "has_response": has_resp,
        "t_cell_only": t_only,
        "public_processed": public_matrix,
        "malignant_matrix": is_scrna and is_lung and not t_only,
        "twoway": False,
        "note": "",
        "soft_chars": len(soft),
    }
    if acc in overrides:
        rec.update(overrides[acc])
    else:
        rec["twoway"] = bool(
            rec["is_lung"]
            and rec["is_scrna"]
            and rec["has_timepoint"]
            and rec["has_response"]
            and rec["malignant_matrix"]
            and rec["public_processed"]
            and not rec["t_cell_only"]
        )
        if rec["twoway"]:
            rec["note"] = "Text-level two-way candidate; inspect GEO files before stacking."
        elif rec["has_timepoint"] and rec["has_response"] and rec["t_cell_only"]:
            rec["note"] = "Has both axes in text but T-cell-restricted deposit."
        elif rec["has_response"] and not rec["has_timepoint"]:
            rec["note"] = rec["note"] or "Response labels without a public pre/post axis."
        elif rec["has_timepoint"] and not rec["has_response"]:
            rec["note"] = rec["note"] or "Timepoints without a public MPR/RECIST axis."
    return rec


def main() -> None:
    found_ids: list[str] = []
    query_hits: dict[str, list[str]] = {}
    for name, term in QUERIES:
        ids = esearch(term)
        query_hits[name] = ids
        found_ids.extend(ids)
        time.sleep(0.34)
    found_ids = list(dict.fromkeys(found_ids))
    summaries = esummary(found_ids) if found_ids else []

    acc_to_sum: dict[str, dict] = {}
    for rec in summaries:
        acc = acc_from_summary(rec)
        if acc:
            acc_to_sum[acc] = rec

    # Resolve seed accessions that eSearch may have missed.
    missing = [a for a in SEED if a not in acc_to_sum]
    if missing:
        # gds search by accession
        extra_ids = []
        for acc in missing:
            extra_ids.extend(esearch(f"{acc}[ACCN]"))
            time.sleep(0.34)
        for rec in esummary(list(dict.fromkeys(extra_ids))):
            acc = acc_from_summary(rec)
            if acc:
                acc_to_sum[acc] = rec

    rows = []
    all_acc = list(dict.fromkeys(SEED + list(acc_to_sum)))
    for acc in all_acc:
        rec = acc_to_sum.get(acc, {})
        title = rec.get("title") or rec.get("entrytitle") or ""
        summary = rec.get("summary") or ""
        # Only pull SOFT for seed + text-level two-axis candidates to keep the hunt cheap.
        pull_soft = acc in SEED or (TIME_RE.search(title + summary) and RESP_RE.search(title + summary))
        soft = fetch_geo_soft_brief(acc) if pull_soft else ""
        time.sleep(0.2)
        row = classify(acc, title, summary, soft)
        row["n_samples"] = rec.get("n_samples")
        row["ftp"] = rec.get("ftplink")
        rows.append(row)

    rows.sort(key=lambda r: (not r["twoway"], not r["has_response"], r["accession"]))
    (OUT / "hunt_series.json").write_text(json.dumps({"queries": query_hits, "series": rows}, indent=2))
    # TSV
    keys = [
        "accession",
        "twoway",
        "has_timepoint",
        "has_response",
        "malignant_matrix",
        "t_cell_only",
        "public_processed",
        "is_lung",
        "is_scrna",
        "n_samples",
        "note",
        "title",
    ]
    lines = ["\t".join(keys)]
    for r in rows:
        lines.append("\t".join("" if r.get(k) is None else str(r.get(k)).replace("\t", " ") for k in keys))
    (OUT / "hunt_series.tsv").write_text("\n".join(lines) + "\n")
    twoway = [r for r in rows if r["twoway"]]
    print(f"scored {len(rows)} series; twoway={len(twoway)} {[r['accession'] for r in twoway]}")


if __name__ == "__main__":
    main()
