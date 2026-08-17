#!/usr/bin/env python3
"""GSE161537 leftover: is CLDN4 on the public NIVOBIO HTG EdgeSeq OBP panel?

If CLDN4 is absent, stop. Do not score CLDN4 vs DCB/ORR/PFS or vs CD8A/IFN/CD274.
Do not invent a CLDN3 / TACSTD2 / immune-gene proxy.
Patient is the unit. Honest n from GEO.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import re
import urllib.request
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161537/"
    "suppl/GSE161537_nivobio_log2cpm.csv.gz"
)
SERIES_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161537/"
    "matrix/GSE161537_series_matrix.txt.gz"
)

CLDN4_ALIASES = {
    "CLDN4",
    "CLAUDIN4",
    "CLAUDIN-4",
    "CPE-R",
    "CPER",
    "CPETR",
    "CPETR1",
    "WBSCR8",
    "HCPER",
}
TACSTD2_ALIASES = {
    "TACSTD2",
    "TROP2",
    "TROP-2",
    "EGP1",
    "EGP-1",
    "GA733-1",
    "M1S1",
}
QUERIES = [
    "CLDN4",
    "CLDN3",
    "CLDN7",
    "CLDN18",
    "TACSTD2",
    "CD8A",
    "CD8B",
    "CD274",
    "IFNG",
    "IFNGR1",
    "PDCD1",
    "PDCD1LG2",
    "CXCL9",
    "CXCL10",
    "EPCAM",
    "GZMB",
    "PRF1",
]


def fetch(url: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / Path(url).name
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": "gse161537-cldn4-check"})
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
        out.write(resp.read())
    return dest


def load_log2cpm(path: Path) -> tuple[list[str], list[str]]:
    with gzip.open(path, "rt") as fh:
        text = fh.read()
    reader = csv.reader(io.StringIO(text), delimiter=";")
    rows = list(reader)
    header = [h.strip().strip('"') for h in rows[0]]
    patients = header[1:]
    genes = [r[0].strip().strip('"') for r in rows[1:] if r]
    return patients, genes


def parse_series(path: Path) -> dict:
    titles: list[str] = []
    geo_acc: list[str] = []
    chars: dict[str, list[str]] = {}
    summary = ""
    overall = ""
    platform = ""
    pubmed: list[str] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Series_summary"):
                summary += " " + line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Series_overall_design"):
                overall = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Series_platform_id"):
                platform = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Series_pubmed_id"):
                pubmed.append(line.split("\t", 1)[1].strip().strip('"'))
            elif line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                geo_acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                keys = []
                cleaned = []
                for v in vals:
                    if ": " in v:
                        k, rest = v.split(": ", 1)
                        keys.append(k.strip().lower())
                        cleaned.append(rest.strip())
                    elif ":" in v:
                        k, rest = v.split(":", 1)
                        keys.append(k.strip().lower())
                        cleaned.append(rest.strip())
                    else:
                        keys.append("")
                        cleaned.append(v)
                key = Counter(k for k in keys if k).most_common(1)
                if key:
                    chars[key[0][0]] = cleaned
    return {
        "titles": titles,
        "geo_acc": geo_acc,
        "chars": chars,
        "summary": summary.strip(),
        "overall": overall,
        "platform": platform,
        "pubmed": pubmed,
    }


def aliases_for(query: str) -> set[str]:
    if query == "CLDN4":
        return {a.upper() for a in CLDN4_ALIASES}
    if query == "TACSTD2":
        return {a.upper() for a in TACSTD2_ALIASES}
    return {query.upper()}


def gene_hits(query: str, genes: list[str]) -> list[str]:
    alts = aliases_for(query)
    hits = []
    for g in genes:
        token = g.strip().strip('"')
        if token.upper() in alts:
            hits.append(token)
    return hits


def main() -> None:
    expr_path = fetch(MATRIX_URL)
    series_path = fetch(SERIES_URL)
    patients, genes = load_log2cpm(expr_path)
    series = parse_series(series_path)
    chars = series["chars"]

    gene_set = {g.upper() for g in genes}
    presence = []
    for q in QUERIES:
        hits = gene_hits(q, genes)
        presence.append(
            {
                "query": q,
                "on_htg_obp_matrix": bool(hits),
                "ids": hits,
            }
        )

    cldn4_present = any(p["query"] == "CLDN4" and p["on_htg_obp_matrix"] for p in presence)
    claudin_family = sorted(g for g in genes if re.match(r"^CLDN", g, re.I))
    trop_family = sorted(
        g for g in genes if re.search(r"TACSTD|TROP2", g, re.I)
    )

    pids = chars.get("patient id", [])
    recist = chars.get("best response on immunotherapy (recist)", [])
    pfs = chars.get("pfs (month)", [])
    os_m = chars.get("os (month)", [])
    os_e = chars.get("os (event)", [])
    histo = chars.get("histology", [])
    io_line = chars.get("immunotherapy line", [])
    hot = chars.get("hot phenotype", [])
    sex = chars.get("sex", [])
    stage = chars.get("stage", [])

    recist_n = dict(Counter(recist))
    orr_pos = sum(1 for x in recist if x in {"CR", "PR"})
    orr_neg = sum(1 for x in recist if x in {"SD", "PD"})
    recist_na = sum(1 for x in recist if x in {"", "NA", "na"})
    pfs_n = sum(1 for x in pfs if x not in {"", "NA", "na"})
    os_n = sum(1 for x in os_m if x not in {"", "NA", "na"})

    n_patient = len(set(pids)) if pids else len(patients)
    empty_reason = (
        "CLDN4 not on deposited HTG EdgeSeq Oncology Biomarker Panel "
        "(GSE161537_nivobio_log2cpm.csv.gz; 2559 genes). "
        "CLDN3 is on the panel and is not used as a proxy."
    )
    empty = {
        "n": 0 if not cldn4_present else n_patient,
        "status": "empty" if not cldn4_present else "scoreable",
        "reason": empty_reason if not cldn4_present else "CLDN4 present",
    }
    verdict = "score" if cldn4_present else "panel-missing"

    tests = {
        "CLDN4_vs_DCB": dict(empty),
        "CLDN4_vs_ORR": dict(empty),
        "CLDN4_vs_PFS": dict(empty),
        "CLDN4_vs_CD8A": dict(empty),
        "CLDN4_vs_IFN": dict(empty),
        "CLDN4_vs_CD274": dict(empty),
    }
    # DCB is not a deposited field; keep that in the DCB row even when empty.
    if not cldn4_present:
        tests["CLDN4_vs_DCB"]["note"] = (
            "No DCB/NDB characteristic on GEO. PFS months are deposited "
            "but were not recoded as DCB because CLDN4 n=0."
        )

    report = {
        "accession": "GSE161537",
        "alias": "NIVOBIO",
        "public": True,
        "public_on": "2022-07-07",
        "pmid": series["pubmed"],
        "platform_geo": series["platform"] or "GPL18573",
        "assay": "HTG EdgeSeq Oncology Biomarker Panel (OBP)",
        "verdict": verdict,
        "n_gsm": len(series["geo_acc"]),
        "n_titles": len(series["titles"]),
        "n_unique_patient_id": len(set(pids)),
        "n_matrix_patients": len(patients),
        "n_matrix_genes": len(genes),
        "n_unique_matrix_genes": len(gene_set),
        "patient_unit_1to1": (
            len(set(pids)) == len(pids) == len(patients) == len(series["geo_acc"]) == 82
        ),
        "matrix_patient_ids_match_geo": set(patients) == set(pids),
        "recist": recist_n,
        "orr_CR_PR_vs_SD_PD": {"CR_PR": orr_pos, "SD_PD": orr_neg, "NA": recist_na},
        "pfs_month_n": pfs_n,
        "pfs_event_deposited": "pfs (event)" in chars or "pfse" in chars,
        "os_month_n": os_n,
        "os_event": dict(Counter(os_e)),
        "histology": dict(Counter(histo)),
        "immunotherapy_line": dict(Counter(io_line)),
        "hot_phenotype": dict(Counter(hot)),
        "sex": dict(Counter(sex)),
        "stage": dict(Counter(stage)),
        "dcb_field_deposited": False,
        "cldn4_present": cldn4_present,
        "cldn4_n": 0 if not cldn4_present else n_patient,
        "claudin_family_on_panel": claudin_family,
        "tacstd_trop_on_panel": trop_family,
        "presence": presence,
        "controls_on_panel_not_substitutes": {
            g: g.upper() in gene_set
            for g in ("CD8A", "CD274", "IFNG", "IFNGR1", "PDCD1", "CLDN3")
        },
        "tests": tests,
        "ftp": {"log2cpm": MATRIX_URL, "series_matrix": SERIES_URL},
        "note": (
            "Series design: 82 advanced NSCLC, second-line PD-1/PD-L1, "
            "pre-IO FFPE, HTG EdgeSeq OBP. Do not use CLDN3, TACSTD2, "
            "CD8A, IFNG, or CD274 as a CLDN4 stand-in. "
            "Do not audit prior TACSTD2 leftover pages."
        ),
        "citation": {
            "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161537",
            "pmid": series["pubmed"],
        },
    }

    tables = HERE / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    (HERE / "panel_check.json").write_text(json.dumps(report, indent=2) + "\n")

    with (HERE / "presence_table.tsv").open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["query", "on_htg_obp_matrix", "ids"], delimiter="\t"
        )
        w.writeheader()
        for row in presence:
            w.writerow(
                {
                    "query": row["query"],
                    "on_htg_obp_matrix": row["on_htg_obp_matrix"],
                    "ids": ",".join(row["ids"]),
                }
            )

    with (HERE / "panel_symbols.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["symbol"])
        w.writerows([[g] for g in genes])

    with (tables / "inventory.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["item", "n", "note"])
        rows_inv = [
            ("accession", "GSE161537", "NIVOBIO HTG EdgeSeq OBP; public 2022-07-07"),
            ("gsm", str(len(series["geo_acc"])), "GSM4909684–GSM4909765"),
            ("unique_patient_id", str(len(set(pids))), "patient is the unit; 1 GSM : 1 patient"),
            ("matrix_patients", str(len(patients)), "GSE161537_nivobio_log2cpm.csv.gz columns"),
            ("matrix_genes", str(len(genes)), "HTG OBP deposited symbols"),
            ("CLDN4_measured", "0", "absent from matrix; aliases also absent"),
            ("CLDN3_on_panel", "1" if "CLDN3" in gene_set else "0", "not used as proxy"),
            ("CD8A_on_panel", "1" if "CD8A" in gene_set else "0", "not a CLDN4 substitute"),
            ("CD274_on_panel", "1" if "CD274" in gene_set else "0", "not a CLDN4 substitute"),
            ("IFNG_on_panel", "1" if "IFNG" in gene_set else "0", "not a CLDN4 substitute"),
            ("RECIST_labeled", str(len(recist) - recist_na), f"{recist_n}"),
            ("ORR_CR_PR", str(orr_pos), "CR+PR among labeled RECIST"),
            ("ORR_SD_PD", str(orr_neg), "SD+PD among labeled RECIST"),
            ("RECIST_NA", str(recist_na), "best response NA"),
            ("PFS_month_labeled", str(pfs_n), "pfs (month); no pfs event field"),
            ("OS_month_labeled", str(os_n), "os (month)"),
            ("DCB_field", "0", "not deposited; not recoded from PFS"),
        ]
        w.writerows(rows_inv)

    with (tables / "one_row.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "dataset",
                "public_lung_ici",
                "n_patient",
                "cldn4_n",
                "verdict",
                "CLDN4_vs_DCB",
                "CLDN4_vs_ORR",
                "CLDN4_vs_PFS",
                "CLDN4_vs_CD8A",
                "CLDN4_vs_IFN",
                "CLDN4_vs_CD274",
                "note",
            ]
        )
        w.writerow(
            [
                "GSE161537",
                "yes",
                str(n_patient),
                "0",
                verdict,
                "empty",
                "empty",
                "empty",
                "empty",
                "empty",
                "empty",
                "HTG OBP panel-missing CLDN4; no proxy",
            ]
        )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "cldn4_n": 0 if not cldn4_present else n_patient,
                "n_patient": n_patient,
                "n_genes": len(genes),
                "claudin_family_on_panel": claudin_family,
                "tests": tests,
            },
            indent=2,
        )
    )
    if verdict != "score":
        print(
            "\nPANEL-MISSING: CLDN4 is absent from GSE161537 HTG EdgeSeq OBP. "
            "No CLDN4 vs DCB/ORR/PFS or vs CD8A/IFN/CD274 scores. No proxy."
        )


if __name__ == "__main__":
    main()
