#!/usr/bin/env python3
"""GSE93157 leftover: is CLDN4 on the public lung ICI NanoString panel?

If CLDN4 is absent, stop. Do not score CLDN4 vs response / CD8 / CD274.
Report honest n from GEO (do not invent p-values).
"""
from __future__ import annotations

import gzip
import json
import re
import urllib.request
from collections import Counter
from pathlib import Path

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93157/"
    "matrix/GSE93157_series_matrix.txt.gz"
)
RAW_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93157/"
    "suppl/GSE93157_raw_data_values.txt.gz"
)
PLATFORM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL19nnn/GPL19965/"
    "soft/GPL19965_family.soft.gz"
)

TARGETS = {
    "CLDN4": (
        "CLDN4",
        "CLAUDIN4",
        "CLAUDIN-4",
        "CPE-R",
        "CPER",
        "CPETR",
        "CPETR1",
        "WBSCR8",
        "NM_001305",
    ),
    "TACSTD2": (
        "TACSTD2",
        "TROP2",
        "TROP-2",
        "EGP1",
        "EGP-1",
        "GA733-1",
        "M1S1",
        "NM_002353",
    ),
}
CONTROLS = ("CD8A", "CD8B", "CD274", "PDCD1", "IFNG", "EPCAM", "CDH1")


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": "gse93157-cldn4-check"})
    with urllib.request.urlopen(req, timeout=120) as resp, dest.open("wb") as out:
        out.write(resp.read())
    return dest


def parse_platform(path: Path) -> list[dict]:
    rows: list[dict] = []
    in_table = False
    header: list[str] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(fh).rstrip("\n").split("\t")
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                parts = line.rstrip("\n").split("\t")
                rec = {header[i]: (parts[i] if i < len(parts) else "") for i in range(len(header))}
                rows.append(rec)
    return rows


def parse_matrix(path: Path) -> dict:
    titles: list[str] = []
    geo_acc: list[str] = []
    sources: list[str] = []
    chars: dict[str, list[str]] = {}
    expr_genes: list[str] = []
    summary = ""
    overall = ""
    platform = ""
    in_table = False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Series_summary"):
                summary += " " + line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Series_overall_design"):
                overall = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Series_platform_id"):
                platform = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                geo_acc = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                sources = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                keys = []
                cleaned = []
                for v in vals:
                    if ":" in v:
                        k, rest = v.split(":", 1)
                        keys.append(k.strip().lower())
                        cleaned.append(rest.strip())
                    else:
                        keys.append("")
                        cleaned.append(v)
                key = Counter(k for k in keys if k).most_common(1)
                if key:
                    chars[key[0][0]] = cleaned
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                gene_id = line.split("\t", 1)[0].strip().strip('"')
                if gene_id in {"ID_REF", ""}:
                    continue
                expr_genes.append(gene_id)
    return {
        "titles": titles,
        "geo_acc": geo_acc,
        "sources": sources,
        "chars": chars,
        "expr_genes": expr_genes,
        "summary": summary.strip(),
        "overall": overall,
        "platform": platform,
    }


def parse_raw_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            first = line.split("\t", 1)[0].strip().strip('"')
            if first in {"", "ID_REF"}:
                continue
            # skip header-ish rows that are not gene symbols
            if first.startswith("FileName") or first.startswith("File"):
                continue
            ids.append(first)
    return ids


def search_rows(rows: list[dict], needles: tuple[str, ...]) -> list[str]:
    hits = []
    for rec in rows:
        blob = " ".join(rec.values()).upper()
        if any(n.upper() in blob for n in needles):
            hits.append(rec.get("ID", ""))
    return sorted({h for h in hits if h})


def search_ids(ids: list[str], needles: tuple[str, ...]) -> list[str]:
    out = []
    for g in ids:
        gu = g.upper()
        if any(n.upper() == gu or n.upper() in gu for n in needles):
            out.append(g)
    return sorted(set(out))


def main() -> None:
    out_dir = Path(__file__).resolve().parent
    cache = out_dir / "cache"
    plat_path = _download(PLATFORM_URL, cache / "GPL19965_family.soft.gz")
    mat_path = _download(MATRIX_URL, cache / "GSE93157_series_matrix.txt.gz")
    raw_path = _download(RAW_URL, cache / "GSE93157_raw_data_values.txt.gz")

    plat = parse_platform(plat_path)
    mat = parse_matrix(mat_path)
    raw_ids = parse_raw_ids(raw_path)
    codeclass = Counter(r.get("Code_Class") or r.get("CodeClass") or "" for r in plat)

    n = len(mat["titles"])
    sources = mat["sources"]
    chars = mat["chars"]
    best = chars.get("best.resp") or chars.get("best.resp.") or []
    geo_resp = chars.get("response") or []
    drug = chars.get("drug") or []
    pfs = chars.get("pfs") or []

    is_lung = [bool(re.search(r"LUNG", s or "", re.I)) for s in sources]
    lung_idx = [i for i, flag in enumerate(is_lung) if flag]
    source_n = dict(Counter(sources))
    lung_source_n = dict(Counter(sources[i] for i in lung_idx))
    lung_best = [best[i] if i < len(best) else "" for i in lung_idx]
    lung_geo_resp = [geo_resp[i] if i < len(geo_resp) else "" for i in lung_idx]
    lung_drug = [drug[i] if i < len(drug) else "" for i in lung_idx]
    lung_pfs = [pfs[i] if i < len(pfs) else "" for i in lung_idx]
    lung_pfs_n = sum(1 for x in lung_pfs if x not in {"", "NA", "na"})

    orr_r = sum(1 for x in lung_best if x in {"CR", "PR"})
    orr_nr = sum(1 for x in lung_best if x in {"SD", "PD"})
    dcb = sum(1 for x in lung_best if x in {"CR", "PR", "SD"})
    ndb = sum(1 for x in lung_best if x == "PD")

    target_hits = {}
    for gene, needles in TARGETS.items():
        plat_hits = search_rows(plat, needles)
        matrix_hits = search_ids(mat["expr_genes"], needles)
        raw_hits = search_ids(raw_ids, needles)
        target_hits[gene] = {
            "on_platform": bool(plat_hits),
            "in_series_matrix": bool(matrix_hits),
            "in_raw_suppl": bool(raw_hits),
            "platform_ids": plat_hits,
            "matrix_ids": matrix_hits,
            "raw_ids": raw_hits,
        }

    controls = {}
    for gene in CONTROLS:
        needles = (gene,)
        controls[gene] = {
            "on_platform": bool(search_rows(plat, needles)),
            "in_series_matrix": gene in mat["expr_genes"],
            "in_raw_suppl": gene in raw_ids,
        }

    claudin_family = sorted(
        {
            r.get("ID", "")
            for r in plat
            if re.search(r"\bCLDN", " ".join(r.values()), re.I)
        }
        | {g for g in mat["expr_genes"] if re.search(r"CLDN", g, re.I)}
        | {g for g in raw_ids if re.search(r"CLDN", g, re.I)}
    )

    present = [
        g
        for g, h in target_hits.items()
        if h["on_platform"] or h["in_series_matrix"] or h["in_raw_suppl"]
    ]
    cldn4_present = any(
        target_hits["CLDN4"][k] for k in ("on_platform", "in_series_matrix", "in_raw_suppl")
    )
    verdict = "score" if cldn4_present else "panel_missing_empty"

    tests = {
        "CLDN4_vs_response": {
            "n": 0 if not cldn4_present else len(lung_idx),
            "status": "empty" if not cldn4_present else "scoreable",
            "reason": "CLDN4 not on GPL19965 / series matrix / raw suppl"
            if not cldn4_present
            else "CLDN4 present",
        },
        "CLDN4_vs_CD8": {
            "n": 0 if not cldn4_present else len(lung_idx),
            "status": "empty" if not cldn4_present else "scoreable",
            "reason": "CLDN4 not measured; CD8A is on the panel but is not a CLDN4 substitute"
            if not cldn4_present
            else "both present",
        },
        "CLDN4_vs_CD274": {
            "n": 0 if not cldn4_present else len(lung_idx),
            "status": "empty" if not cldn4_present else "scoreable",
            "reason": "CLDN4 not measured; CD274 is on the panel but is not a CLDN4 substitute"
            if not cldn4_present
            else "both present",
        },
    }

    report = {
        "accession": "GSE93157",
        "platform": mat["platform"] or "GPL19965",
        "panel": "NanoString nCounter PanCancer Immune Profiling Panel",
        "public_lung_ici": True,
        "verdict": verdict,
        "targets_present": present,
        "targets": target_hits,
        "controls_on_panel": controls,
        "claudin_family_hits": claudin_family,
        "n_geo_arrays": n,
        "n_expr_genes_matrix": len(mat["expr_genes"]),
        "n_raw_feature_ids": len(raw_ids),
        "platform_n_rows": len(plat),
        "platform_codeclass": dict(codeclass),
        "source_n": source_n,
        "n_lung": len(lung_idx),
        "lung_source_n": lung_source_n,
        "lung_best_resp": dict(Counter(lung_best)),
        "lung_geo_response_field": dict(Counter(lung_geo_resp)),
        "lung_drug": dict(Counter(lung_drug)),
        "lung_pfs_n": lung_pfs_n,
        "lung_orr": {"CR_PR": orr_r, "SD_PD": orr_nr},
        "lung_dcb": {"CR_PR_SD": dcb, "PD": ndb},
        "tests": tests,
        "note": (
            "GEO characteristic 'response' is RC_RP_SD on every lung sample, "
            "including PD. Use best.resp (and pfs/pfse), not 'response'. "
            "Do not treat CD8A or CD274 as CLDN4 stand-ins."
        ),
        "citation": {
            "paper": "Prat et al. Cancer Res 2017",
            "pmid": "28487385",
            "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93157",
            "platform": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL19965",
        },
    }

    (out_dir / "panel_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if verdict != "score":
        print(
            "\nEMPTY: CLDN4 is absent from GSE93157. "
            "No CLDN4 vs response / CD8 / CD274 scores."
        )


if __name__ == "__main__":
    main()
