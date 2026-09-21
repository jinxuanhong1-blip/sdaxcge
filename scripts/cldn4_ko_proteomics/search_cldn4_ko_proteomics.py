#!/usr/bin/env python3
"""Search ArrayExpress/BioStudies and PRIDE for a CLDN4 knockout or
knockdown proteome in any human epithelial line, then record MHC/IFN
protein scores only when a qualifying protein matrix exists.

Live queries are sent to EBI Search, the BioStudies API, the PRIDE
Archive v2 keyword search, OmicsDI, and PubMed E-utilities. No accessions
are invented. Protein log2 fold-changes are written only from a downloaded
quantification table; this hunt did not find one that is a CLDN4 loss
contrast.
"""

from __future__ import annotations

import csv
import json
import re
import ssl
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "cldn4_ko_proteomics"
SEARCH_DATE = "2026-09-21"
UA = "cldn4-ko-proteomics-hunt/1.0 (public-data screen)"
CTX = ssl.create_default_context()

EXACT_CLDN4 = re.compile(r"\b(CLDN4|Cldn4|claudin-4|claudin 4)\b", re.I)
PERTURB = re.compile(
    r"\b(knock-?out|knock-?down|CRISPR|siRNA|shRNA|silencing|deficient|deletion)\b",
    re.I,
)
PROTEOMICS = re.compile(
    r"proteom|mass spectrom|LC-MS|PRIDE|PXD\d{6}|iTRAQ|TMT\b|label-free",
    re.I,
)

# Same manual panels as the CLDN4/TACSTD2 KD-KO RNA wave (not an MSigDB export).
GENE_SETS = {
    "IFN_ALPHA_TYPE1": [
        "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "MX1", "MX2", "OAS1",
        "OAS2", "OAS3", "OASL", "STAT1", "STAT2", "IRF7", "IRF9", "DDX58",
        "IFIH1", "BST2", "RSAD2", "USP18", "IFI27", "IFI44", "IFI44L", "IFI6",
        "XAF1", "HERC5", "EPSTI1", "SAMD9L", "CMPK2", "LY6E", "IFITM1",
        "IFITM2", "IFITM3", "PLSCR1", "SP100", "SP110", "EIF2AK2", "ZBP1",
        "TRIM22", "HELZ2",
    ],
    "IFN_GAMMA_TYPE2": [
        "GBP1", "GBP2", "GBP3", "GBP4", "GBP5", "CXCL9", "CXCL10", "CXCL11",
        "IDO1", "STAT1", "IRF1", "CIITA", "SOCS1", "SOCS3", "HLA-DRA",
        "HLA-DRB1", "TAP1", "TAP2", "PSMB8", "PSMB9", "PSMB10", "IFNG",
        "IFNGR1", "IFNGR2", "JAK1", "JAK2", "NLRC5", "B2M", "UBD", "VCAM1",
        "ICAM1",
    ],
    "ANTIGEN_PRESENTATION_MHC1": [
        "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "TAP1", "TAP2",
        "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5", "CALR", "CANX", "PDIA3",
        "ERAP1", "ERAP2", "SEC61A1",
    ],
}

# Calls made after reading the public record on SEARCH_DATE. These are not
# hits. They are the records a keyword search returns that still fail the
# design gate (CLDN4/Cldn4 genetic loss, proteome, human epithelial line).
REVIEWED = [
    {
        "accession": "E-GEOD-22493",
        "repository": "ArrayExpress/BioStudies",
        "call": "near_miss_rna",
        "organism_line": "human ovarian SKOV3",
        "why": "CLDN4 siRNA versus overexpression. Microarray RNA (GSE22493). Not a proteome.",
    },
    {
        "accession": "E-GEOD-50927",
        "repository": "ArrayExpress/BioStudies",
        "call": "near_miss_rna",
        "organism_line": "mouse lung",
        "why": "Cldn4 germline knockout RNA-seq (GSE50927). Not human and not protein.",
    },
    {
        "accession": "E-GEOD-60885",
        "repository": "ArrayExpress/BioStudies",
        "call": "not_a_ko_proteome",
        "organism_line": "human trophoblast",
        "why": "DNA methylation study that names claudin-4. Not a CLDN4 knockout proteome.",
    },
    {
        "accession": "E-GEOD-84742",
        "repository": "ArrayExpress/BioStudies",
        "call": "not_a_ko_proteome",
        "organism_line": "mouse colon",
        "why": "Colonic epithelial differentiation RNA-seq. Not a CLDN4 knockout.",
    },
    {
        "accession": "GSE207704",
        "repository": "GEO (not ArrayExpress proteomics, not PRIDE)",
        "call": "near_miss_rna",
        "organism_line": "human breast MCF7 and T47D",
        "why": "CLDN4 CRISPR RNA-seq. No PRIDE or ArrayExpress proteome for this series.",
    },
    {
        "accession": "S-BSST1967",
        "repository": "BioStudies",
        "call": "near_miss_wrong_contrast",
        "organism_line": "human HCC (MHCC97H, PLC/PRF/5, Huh7 and primary cells)",
        "why": (
            "Public files include ExpAll.xlsx described as protein omics (258294 bytes). "
            "The paired paper (PMC12281411, doi:10.1016/j.xcrm.2025.102208) uses MS for "
            "lenvatinib-tolerant versus mock primary cells, a palmitoyl-proteome, and CLDN4 "
            "IP-MS. CLDN4 shRNA is a functional assay, not the proteomic contrast. "
            "ftp.ebi.ac.uk did not serve ExpAll.xlsx in this session (FTP 421 / HTTPS TLS EOF), "
            "so MHC/IFN proteins were not read from that workbook."
        ),
    },
    {
        "accession": "PXD051838",
        "repository": "PRIDE",
        "call": "false_positive_search",
        "organism_line": "mouse brain endothelial cells",
        "why": "Unquoted EBI Search query CLDN4 returned this Foxf2 knockout. The project text names Cldn5, not Cldn4.",
    },
    {
        "accession": "PXD051839",
        "repository": "PRIDE",
        "call": "false_positive_search",
        "organism_line": "mouse brain vessels",
        "why": "Same Foxf2 endothelial knockout series as PXD051838. Text names Cldn5, not Cldn4.",
    },
    {
        "accession": "PXD066158",
        "repository": "PRIDE",
        "call": "near_miss_other_claudin",
        "organism_line": "human germ-cell tumor cells",
        "why": "Proteome of CLDN6-deficient cells. Not CLDN4. Not scored.",
    },
    {
        "accession": "PXD031094",
        "repository": "PRIDE",
        "call": "near_miss_not_ko",
        "organism_line": "MDCK-C7 (dog) plus human reagents",
        "why": "Pan-claudin co-IP interactome, not a CLDN4 knockout abundance proteome.",
    },
    {
        "accession": "PXD005292",
        "repository": "PRIDE",
        "call": "near_miss_not_ko",
        "organism_line": "human breast claudin-low lines",
        "why": "Glycoproteome of the claudin-low subtype. Not a CLDN4 genetic loss contrast.",
    },
    {
        "accession": "PXD024559",
        "repository": "iProX (OmicsDI hit under query CLDN4)",
        "call": "near_miss_not_ko",
        "organism_line": "not a CLDN4 knockout",
        "why": "APEX2 interactors of PEDV M protein. Indexed because a paper mentions CLDN4.",
    },
    {
        "accession": "PMC12603150",
        "repository": "literature; no PXD or GSE in the full text",
        "call": "near_miss_no_proteome_deposit",
        "organism_line": "human high-grade serous ovarian models",
        "why": "CLDN4 knockdown with ISRE reporter, immunoblot, and flow (PMID 41214101). No proteome accession.",
    },
    {
        "accession": "PMID35373300",
        "repository": "literature",
        "call": "near_miss_no_proteome_deposit",
        "organism_line": "human ovarian models",
        "why": "CLDN4 knockdown scored DNA-repair proteins by immunoblot (PMC8988515 abstract). No PRIDE/ArrayExpress proteome.",
    },
    {
        "accession": "PMID38867360",
        "repository": "literature",
        "call": "near_miss_metabolomics",
        "organism_line": "human ovarian models",
        "why": "CLDN4 CRISPRi. The deposited full text describes UHPLC-MS metabolomics, not a proteome (PMC11218812).",
    },
    {
        "accession": "PMID41016339",
        "repository": "literature; no public matrix",
        "call": "near_miss_rna",
        "organism_line": "human SCLC H1688",
        "why": "CLDN4 knockout RNA-seq reported. No ArrayExpress or PRIDE proteome.",
    },
    {
        "accession": "PMID40892111",
        "repository": "literature",
        "call": "near_miss_not_human_epithelial_proteome",
        "organism_line": "mouse 266-6 acinar cells",
        "why": "CLDN4 knockdown RNA-seq in a pancreatitis model. Not a human epithelial proteome.",
    },
    {
        "accession": "PMID30353739",
        "repository": "PubMed title/abstract co-occurrence",
        "call": "not_a_ko_proteome",
        "organism_line": "human intestinal Caco-2",
        "why": "One of two PubMed title/abstract hits for CLDN4 plus knockdown-language plus proteomics-language. Bifidobacterium treatment of Caco-2, not a CLDN4 knockout proteome.",
    },
    {
        "accession": "PMID20511395",
        "repository": "PubMed title/abstract co-occurrence",
        "call": "not_a_ko_proteome",
        "organism_line": "dog MDCK",
        "why": "The other PubMed title/abstract hit. Plasma-membrane proteomics during H-Ras/TGF-beta EMT. Not a human CLDN4 knockout.",
    },
]


def get_json(url: str, timeout: int = 90) -> tuple[dict, dict]:
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout, context=CTX) as resp:
        headers = {k.lower(): v for k, v in resp.headers.items()}
        payload = json.loads(resp.read().decode())
    return headers, payload


def ebisearch(domain: str, query: str, fields: str = "name,description", size: int = 100) -> dict:
    entries = []
    start = 0
    hit_count = None
    while True:
        url = (
            "https://www.ebi.ac.uk/ebisearch/ws/rest/"
            f"{domain}?query={urllib.parse.quote(query)}"
            f"&format=json&size={size}&start={start}&fields={fields}"
        )
        _, js = get_json(url)
        hit_count = js.get("hitCount")
        batch = js.get("entries") or []
        entries.extend(batch)
        start += len(batch)
        if not batch or hit_count is None or start >= int(hit_count) or start > 500:
            break
        time.sleep(0.15)
    return {"domain": domain, "query": query, "hitCount": hit_count, "n_fetched": len(entries), "entries": entries}


def entry_text(entry: dict) -> tuple[str, str]:
    fields = entry.get("fields") or {}
    name = " ".join(fields.get("name") or [])
    desc = " ".join(fields.get("description") or [])
    return name, desc


def classify_blob(title: str, description: str) -> str:
    blob = f"{title} {description}"
    exact = bool(EXACT_CLDN4.search(blob))
    perturb = bool(PERTURB.search(blob))
    proteome = bool(PROTEOMICS.search(blob))
    if exact and perturb and proteome:
        return "exact_cldn4_and_perturbation_and_proteomics_language"
    if exact and perturb:
        return "exact_cldn4_and_perturbation"
    if exact and proteome:
        return "exact_cldn4_and_proteomics_language"
    if exact:
        return "exact_cldn4_mention"
    return "no_exact_cldn4_string"


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def pride_keyword(keyword: str) -> dict:
    url = (
        "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword="
        + urllib.parse.quote(keyword)
        + "&pageSize=100&page=0"
    )
    headers, js = get_json(url)
    if not isinstance(js, list):
        raise RuntimeError(f"unexpected PRIDE payload for {keyword}: {type(js)}")
    return {
        "keyword": keyword,
        "total_records": headers.get("total_records"),
        "n_returned": len(js),
        "accessions": [p.get("accession") for p in js],
        "titles": [p.get("title") for p in js],
    }


def omicsdi_summary(query: str) -> dict:
    sources: dict[str, int] = {}
    start = 0
    total = None
    repos = []
    while True:
        url = (
            "https://www.omicsdi.org/ws/dataset/search?query="
            + urllib.parse.quote(query)
            + f"&start={start}&size=100"
        )
        _, js = get_json(url)
        if total is None:
            total = js.get("count")
            for facet in js.get("facets") or []:
                if facet.get("id") != "repository":
                    continue
                for row in facet.get("facetValues") or []:
                    repos.append({"label": row.get("label"), "count": row.get("count")})
        batch = js.get("datasets") or []
        for dataset in batch:
            source = dataset.get("source") or "unknown"
            sources[source] = sources.get(source, 0) + 1
        start += len(batch)
        if not batch or total is None or start >= int(total) or start > 500:
            break
        time.sleep(0.15)
    return {
        "query": query,
        "count": total,
        "repository_facet": repos,
        "sources_in_fetched_datasets": sources,
        "n_fetched": start,
        "pride_datasets_fetched": sources.get("pride", 0),
    }


def arrayexpress_api_total(query: str) -> dict:
    url = (
        "https://www.ebi.ac.uk/biostudies/api/v1/arrayexpress/search?pageSize=1&page=0&query="
        + urllib.parse.quote(query)
    )
    _, js = get_json(url)
    return {"query": query, "totalHits": js.get("totalHits"), "isTotalHitsExact": js.get("isTotalHitsExact")}


def pubmed_count(name: str, term: str) -> dict:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=20&term="
        + urllib.parse.quote(term)
    )
    _, js = get_json(url)
    result = js["esearchresult"]
    return {
        "name": name,
        "count": result["count"],
        "pmids": ";".join(result.get("idlist") or []),
        "term": term,
    }


def pride_entry_claudin_mentions(accession: str) -> str:
    url = (
        "https://www.ebi.ac.uk/ebisearch/ws/rest/pride/entry/"
        + accession
        + "?format=json&fields=description"
    )
    _, js = get_json(url)
    desc = js["entries"][0]["fields"]["description"][0]
    hits = re.findall(r"Cldn\d+|CLDN\d+|claudin[\s-]*\d+", desc, flags=re.I)
    return ",".join(hits) if hits else "none"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    queries = {
        "pride_quoted_CLDN4": ("pride", '"CLDN4"'),
        "pride_quoted_Cldn4": ("pride", '"Cldn4"'),
        "pride_quoted_claudin4": ("pride", '"claudin-4"'),
        "pride_quoted_phrase_knockdown": ("pride", '"claudin-4 knockdown"'),
        "pride_quoted_phrase_knockout": ("pride", '"claudin-4 knockout"'),
        "pride_unquoted_CLDN4": ("pride", "CLDN4"),
        "pride_CLDN4_AND_ko_language": (
            "pride",
            "CLDN4 AND (knockout OR knockdown OR siRNA OR shRNA OR CRISPR)",
        ),
        "ae_CLDN4": ("biostudies-arrayexpress", "CLDN4"),
        "ae_quoted_CLDN4": ("biostudies-arrayexpress", '"CLDN4"'),
        "ae_CLDN4_AND_proteomics": ("biostudies-arrayexpress", "CLDN4 AND proteomics"),
        "ae_CLDN4_AND_ko": (
            "biostudies-arrayexpress",
            "CLDN4 AND (knockout OR knockdown OR siRNA OR CRISPR)",
        ),
        "bs_quoted_CLDN4_AND_proteomics": ("biostudies", '"CLDN4" AND proteomics'),
        "bs_quoted_CLDN4_AND_proteome": ("biostudies", '"CLDN4" AND proteome'),
        "bs_CLDN4_AND_ko": (
            "biostudies",
            "CLDN4 AND (knockout OR knockdown OR siRNA OR CRISPR)",
        ),
    }

    search_rows = []
    fetched = {}
    for name, (domain, query) in queries.items():
        print(f"EBI {name}", flush=True)
        result = ebisearch(domain, query)
        fetched[name] = result
        search_rows.append(
            {
                "search_date": SEARCH_DATE,
                "source": "ebisearch",
                "name": name,
                "domain": domain,
                "query": query,
                "hitCount": result["hitCount"],
                "n_fetched": result["n_fetched"],
            }
        )
        time.sleep(0.15)

    # Page the broader BioStudies CLDN4 query only to count exact-string classes.
    print("EBI biostudies CLDN4 page", flush=True)
    bs_all = ebisearch("biostudies", "CLDN4")
    search_rows.append(
        {
            "search_date": SEARCH_DATE,
            "source": "ebisearch",
            "name": "bs_CLDN4",
            "domain": "biostudies",
            "query": "CLDN4",
            "hitCount": bs_all["hitCount"],
            "n_fetched": bs_all["n_fetched"],
        }
    )

    screen_rows = []
    for name, result in {**fetched, "bs_CLDN4": bs_all}.items():
        for entry in result["entries"]:
            title, desc = entry_text(entry)
            screen_rows.append(
                {
                    "query_name": name,
                    "domain": result["domain"],
                    "accession": entry.get("id"),
                    "title": title.replace("\t", " "),
                    "class": classify_blob(title, desc),
                    "description_head": desc.replace("\t", " ")[:400],
                }
            )

    print("PRIDE keyword", flush=True)
    for keyword in ["CLDN4", "Cldn4", "claudin-4", "claudin"]:
        pride = pride_keyword(keyword)
        search_rows.append(
            {
                "search_date": SEARCH_DATE,
                "source": "pride_v2_keyword",
                "name": f"keyword:{keyword}",
                "domain": "pride",
                "query": keyword,
                "hitCount": pride["total_records"],
                "n_fetched": pride["n_returned"],
            }
        )
        if keyword == "claudin":
            (OUT / "pride_keyword_claudin_accessions.txt").write_text(
                "\n".join(f"{a}\t{t}" for a, t in zip(pride["accessions"], pride["titles"])) + "\n"
            )
        time.sleep(0.2)

    print("OmicsDI", flush=True)
    omics = omicsdi_summary("CLDN4")
    (OUT / "omicsdi_cldn4_repositories.json").write_text(json.dumps(omics, indent=2) + "\n")
    search_rows.append(
        {
            "search_date": SEARCH_DATE,
            "source": "omicsdi",
            "name": "omicsdi_CLDN4",
            "domain": "omicsdi",
            "query": "CLDN4",
            "hitCount": omics["count"],
            "n_fetched": omics["pride_datasets_fetched"],
        }
    )
    print("ArrayExpress API", flush=True)
    ae_api = arrayexpress_api_total("CLDN4")
    search_rows.append(
        {
            "search_date": SEARCH_DATE,
            "source": "biostudies_arrayexpress_api",
            "name": "arrayexpress_api_CLDN4",
            "domain": "arrayexpress",
            "query": "CLDN4",
            "hitCount": ae_api["totalHits"],
            "n_fetched": ae_api["totalHits"],
        }
    )

    print("PubMed", flush=True)
    pubmed_terms = [
        (
            "title_abstract_ko_and_proteomics",
            '(CLDN4[Title/Abstract] OR "claudin-4"[Title/Abstract] OR "claudin 4"[Title/Abstract] OR Cldn4[Title/Abstract]) AND (knockout[Title/Abstract] OR knockdown[Title/Abstract] OR siRNA[Title/Abstract] OR shRNA[Title/Abstract] OR CRISPR[Title/Abstract]) AND (proteome[Title/Abstract] OR proteomic[Title/Abstract] OR proteomics[Title/Abstract] OR "mass spectrometry"[Title/Abstract])',
        ),
        (
            "all_fields_ko_and_proteomics",
            '(CLDN4 OR "claudin-4" OR "claudin 4" OR Cldn4) AND (knockout OR knockdown OR siRNA OR shRNA OR CRISPR) AND (proteomics OR proteome OR "mass spectrometry" OR "LC-MS/MS")',
        ),
    ]
    pubmed_rows = []
    for name, term in pubmed_terms:
        pubmed_rows.append(pubmed_count(name, term))
        time.sleep(0.4)

    print("Cldn token check", flush=True)
    token_rows = []
    for accession in ["PXD051838", "PXD051839"]:
        token_rows.append(
            {
                "accession": accession,
                "claudin_tokens_in_ebisearch_description": pride_entry_claudin_mentions(accession),
            }
        )

    qualifying = [
        row
        for row in screen_rows
        if row["class"] == "exact_cldn4_and_perturbation_and_proteomics_language"
    ]
    # S-BSST1967 and similar can trip the language flag. Keep them in the
    # qualifying-language table, but the reviewed call decides the score.
    reviewed_ids = {row["accession"] for row in REVIEWED}

    write_tsv(
        OUT / "search_counts.tsv",
        search_rows,
        ["search_date", "source", "name", "domain", "query", "hitCount", "n_fetched"],
    )
    write_tsv(
        OUT / "screened_entries.tsv",
        screen_rows,
        ["query_name", "domain", "accession", "title", "class", "description_head"],
    )
    write_tsv(
        OUT / "language_flags.tsv",
        [row for row in screen_rows if row["class"] != "no_exact_cldn4_string"],
        ["query_name", "domain", "accession", "title", "class", "description_head"],
    )
    write_tsv(
        OUT / "reviewed_calls.tsv",
        REVIEWED,
        ["accession", "repository", "call", "organism_line", "why"],
    )
    write_tsv(
        OUT / "pubmed_counts.tsv",
        pubmed_rows,
        ["name", "count", "pmids", "term"],
    )
    write_tsv(
        OUT / "pride_false_positive_tokens.tsv",
        token_rows,
        ["accession", "claudin_tokens_in_ebisearch_description"],
    )

    score_rows = []
    for set_name, genes in GENE_SETS.items():
        for gene in genes:
            score_rows.append(
                {
                    "gene": gene,
                    "set_name": set_name,
                    "n_qualifying_proteomes": 0,
                    "log2fc_ko_over_control": "NA",
                    "p_value": "NA",
                    "status": "not_scored",
                    "reason": "no CLDN4 or Cldn4 knockout/knockdown proteome of a human epithelial line in ArrayExpress/BioStudies or PRIDE",
                }
            )
    write_tsv(
        OUT / "mhc_ifn_protein_scores.tsv",
        score_rows,
        [
            "gene",
            "set_name",
            "n_qualifying_proteomes",
            "log2fc_ko_over_control",
            "p_value",
            "status",
            "reason",
        ],
    )
    write_tsv(
        OUT / "qualifying_hits.tsv",
        [
            {
                "accession": row["accession"],
                "query_name": row["query_name"],
                "title": row["title"],
                "reviewed": "yes" if row["accession"] in reviewed_ids else "NO_REVIEW_YET",
            }
            for row in qualifying
        ],
        ["accession", "query_name", "title", "reviewed"],
    )

    summary = {
        "search_date": SEARCH_DATE,
        "qualifying_protein_matrices_scored": 0,
        "mhc_ifn_proteins_in_panel": len(score_rows),
        "mhc_ifn_proteins_with_a_log2fc": 0,
        "language_flag_rows": sum(1 for row in screen_rows if row["class"] != "no_exact_cldn4_string"),
        "search_counts": search_rows,
        "pubmed": pubmed_rows,
        "pride_false_positive_tokens": token_rows,
        "omicsdi": omics,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("qualifying_protein_matrices_scored", "mhc_ifn_proteins_in_panel")}, indent=2))


if __name__ == "__main__":
    main()
