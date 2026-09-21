#!/usr/bin/env python3
"""Record public CLDN4 knockout RNA-seq / proteomics queries.

Queries GEO (gds), SRA, PubMed, Europe PMC, and PRIDE. Writes the raw
counts and a curated call for whether a paired RNA-seq + proteomics
CLDN4 knockout dataset is public.
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

EMAIL = "research@example.com"
NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get_json(url: str) -> dict | list:
    req = urllib.request.Request(url, headers={"User-Agent": "cldn4-kd-match/1.0"})
    with urllib.request.urlopen(req, timeout=90) as resp:
        return json.load(resp)


def esearch(db: str, term: str, retmax: int = 20) -> dict:
    q = urllib.parse.urlencode(
        {
            "db": db,
            "term": term,
            "retmax": retmax,
            "retmode": "json",
            "email": EMAIL,
            "tool": "cldn4_kd_match",
        }
    )
    payload = get_json(NCBI + "esearch.fcgi?" + q)
    res = payload["esearchresult"]
    return {
        "db": db,
        "term": term,
        "count": int(res["count"]),
        "ids": res.get("idlist", []),
    }


def esummary_gds(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    q = urllib.parse.urlencode(
        {"db": "gds", "id": ",".join(ids), "retmode": "json", "email": EMAIL}
    )
    payload = get_json(NCBI + "esummary.fcgi?" + q)
    result = payload["result"]
    rows = []
    for uid in result.get("uids", []):
        rec = result[uid]
        if rec.get("entrytype") != "GSE":
            continue
        rows.append(
            {
                "accession": rec.get("accession"),
                "n_samples": rec.get("n_samples"),
                "title": rec.get("title"),
            }
        )
    return rows


def pride(keyword: str) -> dict:
    q = urllib.parse.urlencode({"keyword": keyword, "pageSize": 50})
    url = "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?" + q
    payload = get_json(url)
    if not isinstance(payload, list):
        return {"keyword": keyword, "count": None, "accessions": [], "raw_type": str(type(payload))}
    return {
        "keyword": keyword,
        "count": len(payload),
        "accessions": [
            {"accession": p.get("accession"), "title": p.get("title")} for p in payload
        ],
    }


def epmc(pmid: str) -> dict:
    q = urllib.parse.urlencode(
        {"query": f"EXT_ID:{pmid} SRC:MED", "format": "json", "resultType": "core"}
    )
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + q
    last = None
    for _ in range(3):
        last = get_json(url)
        hits = (last.get("resultList") or {}).get("result") or []
        if hits:
            hit = hits[0]
            return {
                "pmid": pmid,
                "doi": hit.get("doi"),
                "pmcid": hit.get("pmcid"),
                "hasData": hit.get("hasData"),
                "title": hit.get("title"),
            }
    raise SystemExit(f"EuropePMC returned no hit for {pmid}: {str(last)[:400]}")


def main() -> int:
    out = Path("methods/cldn4_kd_match/tables")
    out.mkdir(parents=True, exist_ok=True)

    queries = [
        ("gds", '(CLDN4 OR Cldn4) AND (knockout OR "knock-out" OR CRISPR OR knockdown OR silencing OR siRNA OR shRNA)'),
        ("gds", 'Cldn4 OR "claudin 4" OR "claudin-4" OR CLDN4'),
        ("gds", "SAA1 AND CLDN4 AND (knockout OR CRISPR OR H1688)"),
        ("gds", "H1688 AND (CLDN4 OR SAA1 OR knockout)"),
        ("sra", "CLDN4 AND (knockout OR CRISPR) AND (lung OR H1688 OR SCLC)"),
        ("sra", "H1688 AND (CLDN4 OR SAA1)"),
        (
            "pubmed",
            '(CLDN4 OR Cldn4 OR "claudin-4") AND (knockout OR CRISPR) AND (proteome OR proteomics OR "mass spectrometry" OR PRIDE OR ProteomeXchange)',
        ),
        (
            "pubmed",
            '(CLDN4 OR Cldn4) AND (knockout OR CRISPR OR knockdown) AND (lung OR H1688 OR SCLC) AND ("RNA-seq" OR RNA-seq OR transcriptome)',
        ),
    ]
    search_rows = []
    gds_titles = {}
    for db, term in queries:
        row = esearch(db, term, retmax=50)
        search_rows.append(row)
        if db == "gds" and row["count"] <= 50:
            gds_titles[term] = esummary_gds(row["ids"])
        print(db, row["count"], term[:80])

    pride_rows = [pride("CLDN4"), pride("claudin-4"), pride("Cldn4")]
    epmc_row = epmc("41016339")

    (out / "search_log.json").write_text(
        json.dumps(
            {
                "searches": search_rows,
                "gds_titles": gds_titles,
                "pride": pride_rows,
                "europepmc_41016339": epmc_row,
            },
            indent=2,
        )
    )

    # Curated call. Accessions below were opened on the GEO record in this run
    # (see gds_titles). They are loss-of-function transcriptomes, not proteomes.
    curated = [
        {
            "accession_or_pmid": "GSE207704",
            "system": "human breast MCF7 and T47D",
            "perturbation": "CRISPR CLDN4-/-",
            "rna_public": "yes",
            "proteomics_public": "no",
            "paired": "no",
            "already_analyzed_pr": "43, 98, 131, 155",
            "note": "RNA-seq FPKM. Not a lung line. No PRIDE companion.",
        },
        {
            "accession_or_pmid": "GSE50927",
            "system": "mouse whole lung, germline Cldn4 KO, ventilator injury",
            "perturbation": "germline knockout",
            "rna_public": "yes",
            "proteomics_public": "no",
            "paired": "no",
            "already_analyzed_pr": "43, 155, 512",
            "note": "RNA-seq. BAL protein leak in the paper is not a mass-spec proteome.",
        },
        {
            "accession_or_pmid": "GSE22493",
            "system": "human ovarian SKOV-3",
            "perturbation": "CLDN4 siRNA vs overexpression arrays",
            "rna_public": "yes_microarray",
            "proteomics_public": "no",
            "paired": "no",
            "already_analyzed_pr": "43, 155",
            "note": "Low-confidence array contrast already catalogued. Not lung.",
        },
        {
            "accession_or_pmid": "PMID:41016339",
            "system": "human SCLC NCI-H1688",
            "perturbation": "CLDN4 knockout, RNA-seq reported",
            "rna_public": "not_found",
            "proteomics_public": "not_described",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "EuropePMC hasData=N, no PMCID. GEO/SRA queries returned no H1688 CLDN4-KO accession. Abstract does not describe proteomics.",
        },
        {
            "accession_or_pmid": "PMID:41697223",
            "system": "mouse germline Cldn4 deletion, abdominal sepsis",
            "perturbation": "germline deletion",
            "rna_public": "not_in_abstract",
            "proteomics_public": "not_in_abstract",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "Gut permeability and flow cytometry. Not a lung-cancer KO transcriptome or proteome.",
        },
        {
            "accession_or_pmid": "PXD066158",
            "system": "germ cell tumor cells",
            "perturbation": "CLDN6-deficient",
            "rna_public": "not_this_record",
            "proteomics_public": "yes_other_gene",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "PRIDE hit for claudin, not a CLDN4 knockout.",
        },
        {
            "accession_or_pmid": "PXD031094",
            "system": "pan-claudin family",
            "perturbation": "interactome",
            "rna_public": "no",
            "proteomics_public": "yes_interactome",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "Shared/specific claudin interactions. Not a CLDN4 KO differential proteome.",
        },
        {
            "accession_or_pmid": "GSE174462",
            "system": "human SCLC NCI-H1688",
            "perturbation": "shFOXM1",
            "rna_public": "yes_other_gene",
            "proteomics_public": "no",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "Only GEO series returned by H1688 AND (CLDN4 OR SAA1 OR knockout). Perturbation is FOXM1, not CLDN4.",
        },
        {
            "accession_or_pmid": "PMID:33006362",
            "system": "A549 and H1299",
            "perturbation": "CRAD knockdown",
            "rna_public": "microarray_other_gene",
            "proteomics_public": "western_blot_only",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "PubMed lung RNA query hit. CLDN4 is a downstream readout of CRAD knockdown, not the perturbed gene.",
        },
        {
            "accession_or_pmid": "PMID:36840413",
            "system": "neuroendocrine carcinoma including SCLC",
            "perturbation": "ELF3 knockdown",
            "rna_public": "yes_other_gene",
            "proteomics_public": "no",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "PubMed lung RNA query hit. Not a CLDN4 knockout.",
        },
        {
            "accession_or_pmid": "PMID:38345099",
            "system": "intestinal epithelium",
            "perturbation": "MUC13 deletion",
            "rna_public": "no",
            "proteomics_public": "not_CLDN4_KO",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "PubMed proteome/knockout query hit. Claudin proteins change after MUC13 loss. Not a CLDN4 knockout proteome.",
        },
        {
            "accession_or_pmid": "PMID:37889067",
            "system": "claudin-low breast cell lines",
            "perturbation": "none (subtype comparison)",
            "rna_public": "CCLE_baseline",
            "proteomics_public": "no",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "PubMed proteome/knockout query hit. Claudin-low subtype, not CLDN4 knockout.",
        },
        {
            "accession_or_pmid": "PMID:19318328",
            "system": "hepatocellular carcinoma FFPE",
            "perturbation": "none",
            "rna_public": "no",
            "proteomics_public": "no",
            "paired": "no",
            "already_analyzed_pr": "",
            "note": "PubMed proteome/knockout query hit. Claudin-4 immunohistochemistry, not a knockout proteome.",
        },
        {
            "accession_or_pmid": "DepMap Expression Public 24Q4",
            "system": "cancer cell lines, including lung",
            "perturbation": "genome-wide CRISPR screen (viability), separate baseline RNA-seq",
            "rna_public": "yes_baseline",
            "proteomics_public": "not_in_24Q4_file_set",
            "paired": "no",
            "already_analyzed_pr": "304, 331, 388",
            "note": "OmicsExpressionProteinCodingGenesTPMLogp1.csv is baseline log2(TPM+1), not expression after CLDN4 knockout. 24Q4 file set has no proteomics matrix.",
        },
    ]
    pd_path = out / "ko_omics_inventory.tsv"
    cols = list(curated[0].keys())
    lines = ["\t".join(cols)]
    for row in curated:
        lines.append("\t".join(str(row[c]).replace("\t", " ") for c in cols))
    pd_path.write_text("\n".join(lines) + "\n")
    print("wrote", pd_path)
    print("europepmc", epmc_row)
    print("pride", [(p["keyword"], p["count"]) for p in pride_rows])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
