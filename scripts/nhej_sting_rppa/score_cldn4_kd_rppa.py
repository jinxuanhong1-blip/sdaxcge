#!/usr/bin/env python3
"""Search PRIDE/MassIVE for the Bitler MCT CLDN4-knockdown RPPA, then score
DNA-PKcs / 53BP1 / XRCC1 from the public supplement when no accession exists.

Yamamoto, Webb, Bitler et al. Mol Cancer Ther 2022 (PMC8988515) measured
OVCAR3 shCtrl vs shCLDN4 by RPPA at the MD Anderson RPPA core. RPPA is not
an MS deposit. This script records the repository search and digitizes the
figshare supplement tables. It does not invent RPPA fold-changes.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "nhej_sting_rppa"
NOTES = ROOT / "notes" / "nhej_sting_rppa"
UA = "nhej-sting-rppa/1.0 (public-supplement-digitization)"

PAPER = {
    "citation": "Yamamoto TM, Webb P, Davis DM, Baumgartner HK, Woodruff ER, Guntupalli SR, Neville M, Behbakht K, Bitler BG. Loss of Claudin-4 Reduces DNA Damage Repair and Increases Sensitivity to PARP Inhibitors. Mol Cancer Ther. 2022;21(4):647-657.",
    "doi": "10.1158/1535-7163.MCT-21-0827",
    "pmcid": "PMC8988515",
    "data_availability": "Data are available upon request from the corresponding author.",
    "rppa": {
        "platform": "MD Anderson RPPA Core (Yiling Lu; NCI CA16672, R50CA221675)",
        "contrast": "OVCAR3 shCtrl (n=3) vs shCLDN4#1 (n=3)",
        "call": "proteins with p<0.05 and FDR<15% (Figure 3A heatmap)",
        "validation": "53BP1 and XRCC1 immunoblot, shCtrl / shCLDN4#1 / shCLDN4#2 (Figure 3B-C)",
    },
}

# Directions below are manuscript statements, not measured folds.
# numeric_fold is null because the supplement has no RPPA matrix.
RPPA_CALLS = [
    {
        "protein": "DNA-PKcs",
        "gene": "PRKDC",
        "repair_role": "NHEJ catalytic kinase",
        "rppa_named_differential": False,
        "kd_protein_direction": "not_reported",
        "kd_logic_score": 0,
        "manuscript_statement": "Not named in the RPPA result, the Figure 3A arrow legend, or any supplement table. Absence is not evidence that the antibody was unchanged or absent from the array.",
    },
    {
        "protein": "53BP1",
        "gene": "TP53BP1",
        "repair_role": "NHEJ promoter; HR antagonist",
        "rppa_named_differential": True,
        "kd_protein_direction": "down",
        "kd_logic_score": -1,
        "manuscript_statement": "Two shCLDN4 constructs significantly reduced 53BP1 protein (Figure 3A-B). Olaparib-induced 53BP1 foci were blunted (Figure 4D-E). Numeric fold is not in the supplement.",
    },
    {
        "protein": "XRCC1",
        "gene": "XRCC1",
        "repair_role": "SSBR/BER scaffold; paper groups it with the DNA-repair RPPA hits",
        "rppa_named_differential": True,
        "kd_protein_direction": "down",
        "kd_logic_score": -1,
        "manuscript_statement": "Two shCLDN4 constructs significantly reduced XRCC1 protein (Figure 3A, 3C). Immunofluorescence after olaparib also showed lower XRCC1 (Figure 4C). Numeric fold is not in the supplement.",
    },
]

# Genes pulled from Table S1 (TCGA RNA, CLDN4-high vs CLDN4-low). Not the RPPA.
PANEL = [
    "PRKDC",
    "TP53BP1",
    "XRCC1",
    "XRCC4",
    "XRCC5",
    "XRCC6",
    "LIG4",
    "NHEJ1",
    "PARP1",
    "H2AX",
    "STING1",
    "CGAS",
    "TBK1",
    "IRF3",
    "CCL5",
    "CXCL10",
    "CLDN4",
    "BRIP1",
    "NUP160",
    "POLR2J",
]

DNA_DRUGS = (
    "olaparib|rucaparib|niraparib|talazoparib|cisplatin|carboplatin|"
    "doxorubicin|teniposide|etoposide|topotecan|gemcitabine|bleomycin|"
    "temozolomide|nu7441|peposertib|nedisertib|m3814|dna-pk"
)

FIGSHARE = {
    "S1": "https://ndownloader.figshare.com/files/39984958",
    "S2": "https://ndownloader.figshare.com/files/39984961",
    "S3": "https://ndownloader.figshare.com/files/39984964",
    "S4": "https://ndownloader.figshare.com/files/39984967",
}

PRIDE_KEYWORDS = [
    "CLDN4",
    "claudin-4",
    "claudin 4",
    "shCLDN4",
    "Bitler",
    "claudin",
    "OVCAR3",
]

OMICSDI_QUERIES = [
    '(source:pride OR source:massive) AND (CLDN4 OR "claudin-4" OR shCLDN4)',
    "(source:pride) AND (CLDN4 OR \"claudin-4\")",
    "(source:massive) AND (CLDN4 OR claudin OR shCLDN4)",
    "shCLDN4",
    'RPPA AND (CLDN4 OR "claudin-4" OR OVCAR3)',
]


def http_json(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()
    return json.loads(raw.decode("utf-8", errors="replace")) if raw else None


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120) as resp:
        dest.write_bytes(resp.read())


def pride_keyword(keyword: str, page_size: int = 100, max_pages: int = 3) -> list[dict]:
    hits: list[dict] = []
    for page in range(max_pages):
        qs = urllib.parse.urlencode({"keyword": keyword, "pageSize": page_size, "page": page})
        url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?{qs}"
        data = http_json(url) or []
        if not isinstance(data, list) or not data:
            break
        hits.extend(data)
        if len(data) < page_size:
            break
        time.sleep(0.05)
    return hits


def ebi_pride_count(query: str) -> int:
    qs = urllib.parse.urlencode({"query": query, "format": "json", "size": 1})
    url = f"https://www.ebi.ac.uk/ebisearch/ws/rest/pride?{qs}"
    data = http_json(url) or {}
    return int(data.get("hitCount") or 0)


def omicsdi(query: str, size: int = 15) -> dict:
    qs = urllib.parse.urlencode({"query": query, "size": size, "start": 0})
    url = f"https://www.omicsdi.org/ws/dataset/search?{qs}"
    data = http_json(url) or {}
    datasets = []
    for ds in data.get("datasets") or []:
        datasets.append(
            {
                "source": ds.get("source"),
                "id": ds.get("id"),
                "title": ds.get("title"),
            }
        )
    return {"query": query, "count": int(data.get("count") or 0), "top": datasets}


def load_table_s1(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, header=1)
    df.columns = [str(c).strip() for c in df.columns]
    df["Gene"] = df["Gene"].astype(str).str.strip()
    return df


def load_table_s3(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, header=3)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "cache"
    cache.mkdir(parents=True, exist_ok=True)

    search_log: dict = {"paper": PAPER, "pride": [], "ebi_pride_phrase": [], "omicsdi": []}

    claudin_rows = []
    for kw in PRIDE_KEYWORDS:
        hits = pride_keyword(kw)
        search_log["pride"].append({"keyword": kw, "n": len(hits)})
        print(f"PRIDE {kw!r}: {len(hits)}", flush=True)
        if kw == "claudin":
            for rec in hits:
                title = rec.get("title") or ""
                acc = rec.get("accession")
                blob = f"{title} {rec.get('projectDescription') or ''}".lower()
                claudin_rows.append(
                    {
                        "accession": acc,
                        "title": title,
                        "mentions_cldn4": ("cldn4" in blob or "claudin-4" in blob or "claudin 4" in blob),
                        "mentions_knockdown": ("knockdown" in blob or "shrna" in blob or "silencing" in blob),
                        "is_rppa": "rppa" in blob or "reverse phase" in blob,
                    }
                )
        time.sleep(0.05)

    for phrase in ['"CLDN4"', '"claudin-4"', '"claudin 4"', "shCLDN4", "Bitler"]:
        n = ebi_pride_count(phrase)
        search_log["ebi_pride_phrase"].append({"query": phrase, "hitCount": n})
        print(f"EBI-PRIDE {phrase}: {n}", flush=True)

    for q in OMICSDI_QUERIES:
        rec = omicsdi(q)
        search_log["omicsdi"].append(rec)
        print(f"OmicsDI {q}: {rec['count']}", flush=True)
        time.sleep(0.15)

    pride_hits = pd.DataFrame(claudin_rows)
    pride_hits.to_csv(OUT / "pride_claudin_keyword_hits.tsv", sep="\t", index=False)

    px_kd = pride_hits[
        pride_hits["mentions_cldn4"] & (pride_hits["mentions_knockdown"] | pride_hits["is_rppa"])
    ] if len(pride_hits) else pride_hits
    omics_px = [
        r for r in search_log["omicsdi"]
        if r["query"].startswith("(source:pride") or r["query"].startswith("(source:massive")
    ]
    accession_found = bool(len(px_kd)) or any(r["count"] > 0 for r in omics_px)
    search_log["cldn4_kd_rppa_accession"] = None
    search_log["accession_found"] = accession_found
    search_log["conclusion"] = (
        "No PRIDE or MassIVE accession for the OVCAR3 shCLDN4 RPPA. "
        "Phrase search of PRIDE for CLDN4 / claudin-4 / shCLDN4 is empty. "
        "OmicsDI restricted to source:pride or source:massive is empty for those terms. "
        "The paper states data are available on request. Supplement tables were digitized instead."
    )
    (OUT / "search_log.json").write_text(json.dumps(search_log, indent=2) + "\n")

    for key, url in FIGSHARE.items():
        download(url, cache / f"table_{key}.xlsx")

    s1 = load_table_s1(cache / "table_S1.xlsx")
    qcol = "q-Value"
    logcol = "Log2 of Ratio (unlogged CLDN4 High/CLDN4 Low)"
    s1[qcol] = pd.to_numeric(s1[qcol], errors="coerce")
    n_fdr = int((s1[qcol] < 0.05).sum())
    panel = s1[s1["Gene"].isin(PANEL)].copy()
    panel["in_paper_fdr_0.05"] = panel[qcol] < 0.05
    panel.to_csv(OUT / "table_s1_nhej_sting_panel.tsv", sep="\t", index=False)

    s3 = load_table_s3(cache / "table_S3.xlsx")
    screen = s3[["Compound", "Class", "Cell Viability in shCTRL", "Cell Viability in shCLDN4"]].copy()
    screen = screen[screen["Compound"].astype(str).str.contains(DNA_DRUGS, case=False, na=False)]
    screen["viability_ratio_kd_over_ctrl"] = (
        pd.to_numeric(screen["Cell Viability in shCLDN4"], errors="coerce")
        / pd.to_numeric(screen["Cell Viability in shCTRL"], errors="coerce")
    )
    screen.to_csv(OUT / "table_s3_dna_damage_screen.tsv", sep="\t", index=False)

    lookup = panel.set_index("Gene")
    rows = []
    for call in RPPA_CALLS:
        gene = call["gene"]
        rec = lookup.loc[gene]
        q = float(rec[qcol])
        log2 = float(rec[logcol])
        if q < 0.05 and log2 < 0:
            tcga_sign = -1  # higher transcript in CLDN4-low
        elif q < 0.05 and log2 > 0:
            tcga_sign = 1
        else:
            tcga_sign = 0
        rows.append(
            {
                **call,
                "kd_numeric_fold": None,
                "kd_numeric_source": "not_in_supplement",
                "tcga_table": "Table S1, TCGA ovarian transcriptome, CLDN4-high vs CLDN4-low; not the RPPA",
                "tcga_log2_high_over_low": log2,
                "tcga_mean_log2_high": float(rec["Mean Log2 Expression in CLDN4 High"]),
                "tcga_mean_log2_low": float(rec["Mean Log2 Expression in CLDN4 Low"]),
                "tcga_pvalue": float(rec["p-Value"]),
                "tcga_qvalue": q,
                "tcga_higher_in": str(rec["Higher expression in"]).strip(),
                "tcga_fdr05": bool(q < 0.05),
                "tcga_sign": tcga_sign,
            }
        )
    score = pd.DataFrame(rows)
    score.to_csv(OUT / "nhej_score.tsv", sep="\t", index=False)

    summary = {
        "accession": None,
        "table_s1_genes": int(s1["Gene"].nunique()),
        "table_s1_fdr05": n_fdr,
        "paper_stated_fdr05": 1582,
        "fdr_count_matches_paper": n_fdr == 1582,
        "scores": rows,
    }
    (OUT / "score_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in summary if k != "scores"}, indent=2))
    print(score[["protein", "gene", "kd_logic_score", "kd_protein_direction", "tcga_log2_high_over_low", "tcga_qvalue", "tcga_sign"]].to_string(index=False))


if __name__ == "__main__":
    main()
