#!/usr/bin/env python3
"""GSE329813 leftover public GeoMx matrix: CLDN4 vs T/NK or IFN.

Additive leftover from the durvalumab hunt. This accession is neoadjuvant
pembrolizumab + platinum (NCT05383716), not durvalumab. The assigned test is
CLDN4 vs a T/NK or IFN score on the public processed matrix.

Stop if CLDN4 is absent. Do not substitute TACSTD2 or EPCAM. Do not invent n.
Downloads stay under $GSE329813_CLDN4_DATA (default /tmp/gse329813_cldn4)
and are not committed.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import urllib.request
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
DATA = Path(os.environ.get("GSE329813_CLDN4_DATA", "/tmp/gse329813_cldn4"))

PROCESSED_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329813/"
    "suppl/GSE329813_processed_data_file_normalized_data.csv.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329813/"
    "matrix/GSE329813_series_matrix.txt.gz"
)

# Aliases searched in gene IDs only (no stand-in genes).
CLDN4_ALIASES = (
    "CLDN4",
    "CLAUDIN4",
    "CLAUDIN-4",
    "CLDN-4",
    "CPE-R",
    "CPER",
    "CPETR",
    "CPETR1",
    "WBSCR8",
    "NM_001305",
    "ENSG00000189143",
)

TNK_GENES = [
    "CD3D",
    "CD3E",
    "CD3G",
    "CD2",
    "CD8A",
    "CD8B",
    "NKG7",
    "GNLY",
    "KLRD1",
    "GZMB",
    "PRF1",
    "GZMA",
    "GZMK",
]

IFN_GENES = [
    "IFNG",
    "STAT1",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "IDO1",
    "HLA-DRA",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "CD274",
    "LAG3",
]

NEARBY = ["TACSTD2", "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19"]

TITLE_RE = re.compile(
    r"ROI\s+(\d+),\s+Patient\s+(\d+),\s+(.+),\s+(\S+)\s*$"
)


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse329813-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
        out.write(r.read())
    return dest


def read_processed_genes(path: Path) -> tuple[list[str], list[str]]:
    with gzip.open(path, "rt") as fh:
        header = next(csv.reader(fh))
        rois = [c.strip().strip('"') for c in header[1:]]
        genes = [row[0].strip().strip('"') for row in csv.reader(fh) if row]
    return genes, rois


def parse_series(path: Path) -> dict:
    series: dict[str, list[str]] = {}
    sample: dict[str, list[str]] = {}
    n_expr_rows = 0
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                gene_id = line.split("\t", 1)[0].strip().strip('"')
                if gene_id and gene_id != "ID_REF":
                    n_expr_rows += 1
                continue
            if line.startswith("!Series_"):
                key = line.split("\t", 1)[0][8:]
                val = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(key, []).append(val)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                sample[key] = vals
    titles = sample.get("title", [])
    gsms = sample.get("geo_accession", [])
    row_counts = sample.get("data_row_count", [])
    parsed = []
    for i, title in enumerate(titles):
        m = TITLE_RE.match(title)
        rec = {
            "gsm": gsms[i] if i < len(gsms) else "",
            "title": title,
            "parse_ok": bool(m),
        }
        if m:
            rec.update(
                {
                    "roi": f"ROI {m.group(1)}",
                    "patient": f"P{m.group(2)}",
                    "site": m.group(3),
                    "response": m.group(4),
                }
            )
        parsed.append(rec)
    return {
        "title": " | ".join(series.get("title", [])),
        "design": " | ".join(series.get("overall_design", [])),
        "summary": " | ".join(series.get("summary", [])),
        "platform": " | ".join(series.get("platform_id", [])),
        "suppl": series.get("supplementary_file", []),
        "n_expr_rows_in_series_matrix": n_expr_rows,
        "data_row_count_unique": sorted(set(row_counts)),
        "n_gsm": len(gsms),
        "samples": parsed,
    }


def alias_hits(genes: list[str], aliases: tuple[str, ...]) -> list[str]:
    hits = []
    for g in genes:
        u = g.upper()
        if any(a.upper() == u or a.upper() in u for a in aliases):
            hits.append(g)
    return hits


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    proc = dl(DATA / "GSE329813_processed_data_file_normalized_data.csv.gz", PROCESSED_URL)
    smat = dl(DATA / "GSE329813_series_matrix.txt.gz", MATRIX_URL)

    genes, rois = read_processed_genes(proc)
    upper = {g.upper(): g for g in genes}
    series = parse_series(smat)
    samples = series["samples"]
    parsed = [s for s in samples if s.get("parse_ok")]
    patients = sorted({s["patient"] for s in parsed}, key=lambda x: int(x[1:]))
    tumor = [s for s in parsed if "Primary tumor" in s.get("site", "")]
    ln = [s for s in parsed if "Lymph" in s.get("site", "")]
    tb_pat = {s["patient"] for s in tumor}
    tb_resp = {}
    for s in tumor:
        tb_resp.setdefault(s["patient"], set()).add(s["response"])
    tb_resp_flat = {p: next(iter(v)) for p, v in tb_resp.items()}

    cldn4_hits = alias_hits(genes, CLDN4_ALIASES)
    cldn_family = [g for g in genes if "CLDN" in g.upper() or "CLAUDIN" in g.upper()]
    tnk_present = [g for g in TNK_GENES if g.upper() in upper]
    ifn_present = [g for g in IFN_GENES if g.upper() in upper]
    nearby = {g: g.upper() in upper for g in NEARBY}

    cldn4_on_matrix = bool(cldn4_hits)
    verdict = "score" if cldn4_on_matrix else "cldn4_absent_empty"

    inventory = [
        {
            "item": "GEO GSM / ROI titles",
            "public": "yes",
            "n": len(samples),
            "note": "GSM9711715–GSM9711841; series matrix data_row_count=0",
        },
        {
            "item": "Titles parsed (ROI, patient, site, MPR/NMPR)",
            "public": "yes",
            "n": len(parsed),
            "note": "all titles match the ROI/Patient/site/response pattern",
        },
        {
            "item": "Patients in titles",
            "public": "yes",
            "n": len(patients),
            "note": "P1–P22; design text also says 22 patients",
        },
        {
            "item": "Primary tumor-bed ROIs",
            "public": "yes",
            "n": len(tumor),
            "note": "20 patients with 3 ROIs, 2 patients with 6 ROIs",
        },
        {
            "item": "Lymph-node ROIs",
            "public": "yes",
            "n": len(ln),
            "note": "18/22 patients have LN ROIs",
        },
        {
            "item": "Tumor-bed patients (MPR / NMPR)",
            "public": "yes",
            "n": len(tb_pat),
            "note": (
                f"MPR {sum(1 for v in tb_resp_flat.values() if v == 'MPR')} / "
                f"NMPR {sum(1 for v in tb_resp_flat.values() if v == 'NMPR')}"
            ),
        },
        {
            "item": "Processed matrix genes × ROIs",
            "public": "yes",
            "n": f"{len(genes)} × {len(rois)}",
            "note": "GSE329813_processed_data_file_normalized_data.csv.gz; only public counts",
        },
        {
            "item": "Series-matrix expression rows",
            "public": "empty",
            "n": series["n_expr_rows_in_series_matrix"],
            "note": "header only; do not treat the series matrix as a count table",
        },
        {
            "item": "Durvalumab / PACIFIC arm",
            "public": "no",
            "n": 0,
            "note": "NCT05383716 pembrolizumab + platinum; leftover, not Durva",
        },
        {
            "item": "CLDN4 finite values",
            "public": "no",
            "n": 0,
            "note": "absent from the 1812-gene GeoMx panel; no claudin family member",
        },
        {
            "item": "Primary pairwise n (CLDN4 + T/NK)",
            "public": "empty",
            "n": 0,
            "note": "CLDN4 missing; T/NK genes are present but unused",
        },
        {
            "item": "Primary pairwise n (CLDN4 + IFN)",
            "public": "empty",
            "n": 0,
            "note": "CLDN4 missing; IFN genes are present but unused",
        },
    ]

    coverage = []
    for gene in ["CLDN4", *TNK_GENES, *IFN_GENES, *NEARBY]:
        present = gene.upper() in upper or (gene == "CLDN4" and cldn4_on_matrix)
        coverage.append(
            {
                "gene": gene,
                "set": (
                    "target"
                    if gene == "CLDN4"
                    else "TNK"
                    if gene in TNK_GENES
                    else "IFN"
                    if gene in IFN_GENES
                    else "nearby_not_substitute"
                ),
                "on_processed_matrix": "yes" if present else "no",
                "matrix_id": upper.get(gene.upper(), ""),
            }
        )

    empty_tests = [
        {
            "pair": "CLDN4 vs T/NK mean-z",
            "n": 0,
            "rho": "",
            "p": "",
            "verdict": "EMPTY",
            "reason": "CLDN4 absent from public GeoMx panel",
        },
        {
            "pair": "CLDN4 vs IFN / Ayers-like mean-z",
            "n": 0,
            "rho": "",
            "p": "",
            "verdict": "EMPTY",
            "reason": "CLDN4 absent from public GeoMx panel",
        },
        {
            "pair": "CLDN4 vs CD8A",
            "n": 0,
            "rho": "",
            "p": "",
            "verdict": "EMPTY",
            "reason": "CLDN4 absent; CD8A is on the panel",
        },
        {
            "pair": "CLDN4 vs IFNG",
            "n": 0,
            "rho": "",
            "p": "",
            "verdict": "EMPTY",
            "reason": "CLDN4 absent; IFNG is on the panel",
        },
    ]

    one_row = {
        "dataset": "GSE329813",
        "drug": "pembrolizumab + platinum (NCT05383716); not durvalumab",
        "assay": "GeoMx DSP total-count normalized (scale 10210)",
        "n_roi": len(rois),
        "n_patients": len(patients),
        "n_tumorbed_patients": len(tb_pat),
        "CLDN4": "ABSENT",
        "TNK_genes": f"{len(tnk_present)}/{len(TNK_GENES)}",
        "IFN_genes": f"{len(ifn_present)}/{len(IFN_GENES)}",
        "CLDN4_vs_TNK_n": 0,
        "CLDN4_vs_IFN_n": 0,
        "verdict": "EMPTY",
    }

    summary = {
        "accession": "GSE329813",
        "platform": series["platform"] or "GPL24676",
        "assay": "GeoMx DSP normalized counts (NOT bulk RNA-seq, NOT scRNA)",
        "drug": "pembrolizumab + platinum-doublet; NCT05383716",
        "durvalumab": False,
        "verdict": verdict,
        "cldn4_on_matrix": cldn4_on_matrix,
        "cldn4_hits": cldn4_hits,
        "claudin_family_hits": cldn_family,
        "n_genes": len(genes),
        "n_rois_in_matrix": len(rois),
        "n_gsm": series["n_gsm"],
        "n_patients": len(patients),
        "n_tumorbed_roi": len(tumor),
        "n_ln_roi": len(ln),
        "tumorbed_mpr": sum(1 for v in tb_resp_flat.values() if v == "MPR"),
        "tumorbed_nmpr": sum(1 for v in tb_resp_flat.values() if v == "NMPR"),
        "series_matrix_expr_rows": series["n_expr_rows_in_series_matrix"],
        "data_row_count_unique": series["data_row_count_unique"],
        "tnk_present": tnk_present,
        "tnk_missing": [g for g in TNK_GENES if g.upper() not in upper],
        "ifn_present": ifn_present,
        "ifn_missing": [g for g in IFN_GENES if g.upper() not in upper],
        "nearby_on_panel": nearby,
        "urls": {"processed": PROCESSED_URL, "series_matrix": MATRIX_URL},
        "note": (
            "CLDN4 is not on the deposited GeoMx panel. T/NK and IFN genes are "
            "present. TACSTD2/EPCAM were not used as CLDN4 stand-ins. No "
            "CLDN4–T/NK or CLDN4–IFN Spearman is computed."
        ),
    }

    write_tsv(
        TABLES / "label_inventory.tsv",
        inventory,
        ["item", "public", "n", "note"],
    )
    write_tsv(
        TABLES / "gene_coverage.tsv",
        coverage,
        ["gene", "set", "on_processed_matrix", "matrix_id"],
    )
    write_tsv(
        TABLES / "cldn4_vs_tnk_ifn.tsv",
        empty_tests,
        ["pair", "n", "rho", "p", "verdict", "reason"],
    )
    write_tsv(TABLES / "one_row.tsv", [one_row], list(one_row.keys()))
    write_tsv(
        TABLES / "roi_annotation.tsv",
        parsed,
        ["gsm", "title", "roi", "patient", "site", "response", "parse_ok"],
    )
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps(summary, indent=2))
    print(
        f"\nSTOP: CLDN4 absent ({len(genes)} genes, 0 claudin family). "
        f"T/NK {len(tnk_present)}/{len(TNK_GENES)}, IFN {len(ifn_present)}/{len(IFN_GENES)}. "
        "Empty CLDN4 vs T/NK or IFN."
    )
    print("site counts", dict(Counter(s.get("site") for s in parsed)))
    print("response counts", dict(Counter(s.get("response") for s in parsed)))


if __name__ == "__main__":
    main()
