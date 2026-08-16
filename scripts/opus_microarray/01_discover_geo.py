#!/usr/bin/env python3
"""Step 1 - discover candidate GEO series.

Goal: enumerate *older, array/panel-based* (i.e. non-RNA-seq) GEO Series that
profile lung cancer patients treated with PD-1/PD-L1 checkpoint inhibitors and
report a clinical response variable.

We query NCBI E-utilities (db=gds) with several complementary term sets, pull
esummary records for every hit, and emit a candidate table. Nothing is filtered
silently: the full hit list and the reason each series was kept or dropped are
written to results/opus_microarray/.

Outputs
  results/opus_microarray/geo_search_hits.tsv        all GSE hits with metadata
  results/opus_microarray/geo_candidates.tsv         array/panel + lung + ICI shortlist
  results/opus_microarray/geo_search_queries.json    exact queries used (provenance)
"""

from __future__ import annotations

import csv
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import CACHE_DIR, RESULTS_DIR, ensure_dirs, eutils, open_maybe_gzip, write_json  # noqa: E402

ICI_TERMS = (
    'nivolumab OR pembrolizumab OR atezolizumab OR durvalumab OR avelumab '
    'OR "PD-1" OR "PD-L1" OR "PD1 blockade" OR "immune checkpoint" '
    'OR ipilimumab OR "checkpoint inhibitor" OR "checkpoint blockade"'
)
LUNG_TERMS = 'lung OR NSCLC OR "non-small cell" OR "non small cell" OR pulmonary OR adenocarcinoma'

QUERIES = {
    # Broad: any human GSE mentioning ICI + lung.
    "broad_ici_lung": f'({ICI_TERMS}) AND ({LUNG_TERMS}) AND "Homo sapiens"[Organism] AND "gse"[Entry Type]',
    # Array-typed only (excludes high-throughput sequencing entries).
    "array_ici_lung": (
        f'({ICI_TERMS}) AND ({LUNG_TERMS}) AND "Homo sapiens"[Organism] AND "gse"[Entry Type] '
        f'AND "Expression profiling by array"[DataSet Type]'
    ),
    # Response/benefit wording, array typed.
    "array_response": (
        '(responder OR nonresponder OR "non-responder" OR RECIST OR "durable clinical benefit" '
        'OR "best response" OR "objective response") AND '
        f'({ICI_TERMS}) AND ({LUNG_TERMS}) AND "Homo sapiens"[Organism] AND "gse"[Entry Type] '
        'AND "Expression profiling by array"[DataSet Type]'
    ),
    # NanoString / targeted panel wording (GSE93157 lives here).
    "panel_ici_lung": (
        f'(nCounter OR NanoString OR "PanCancer" OR HTG OR "EdgeSeq") AND ({ICI_TERMS}) '
        f'AND ({LUNG_TERMS}) AND "Homo sapiens"[Organism] AND "gse"[Entry Type]'
    ),
}

# Series explicitly named in the task or known from the ICI-biomarker literature.
# They are force-included in the hit table even if the queries miss them, so the
# shortlist is not silently dependent on Entrez indexing quirks.
SEED_GSE = [
    "GSE93157",   # Prat 2017, nCounter PanCancer Immune 730, NSCLC/HNSCC/melanoma anti-PD-1
    "GSE119144",  # Hwang 2020, NSCLC anti-PD-1, expression array
    "GSE126044",  # Cho 2020, NSCLC anti-PD-1 (RNA-seq; kept for provenance)
    "GSE135222",  # Jung 2019, NSCLC anti-PD-1 (RNA-seq; kept for provenance)
    "GSE136961",
    "GSE99070",
    "GSE111414",
    "GSE190265",
    "GSE168204",
    "GSE207422",
    "GSE91061",   # melanoma anti-PD-1 (out of tissue scope; provenance only)
    "GSE78220",   # melanoma anti-PD-1 (out of tissue scope; provenance only)
]

ARRAY_TYPE_RE = re.compile(r"expression profiling by array", re.I)
SEQ_TYPE_RE = re.compile(r"high throughput sequencing|expression profiling by high", re.I)
LUNG_RE = re.compile(r"\blung\b|nsclc|non[- ]small cell|pulmonary|lusc|luad", re.I)
ICI_RE = re.compile(
    r"nivolumab|pembrolizumab|atezolizumab|durvalumab|avelumab|ipilimumab|"
    r"pd-?1|pd-?l1|checkpoint|immunotherap",
    re.I,
)
RESP_RE = re.compile(
    r"respon|recist|benefit|progress|refractor|resistan|sensitiv|outcome|survival|pfs|\bos\b", re.I
)


def esearch(name: str, term: str, *, retmax: int = 500) -> list[str]:
    dest = CACHE_DIR / "eutils" / f"esearch_{name}.xml"
    eutils("esearch.fcgi", {"db": "gds", "retmax": str(retmax), "term": term}, dest)
    root = ET.parse(dest).getroot()
    return [e.text for e in root.findall(".//IdList/Id") if e.text]


def esummary(uids: list[str], tag: str) -> list[dict]:
    out: list[dict] = []
    for i in range(0, len(uids), 200):
        chunk = uids[i : i + 200]
        dest = CACHE_DIR / "eutils" / f"esummary_{tag}_{i // 200}.xml"
        eutils("esummary.fcgi", {"db": "gds", "id": ",".join(chunk)}, dest)
        root = ET.parse(dest).getroot()
        for doc in root.findall(".//DocSum"):
            rec: dict[str, str] = {"uid": doc.findtext("Id", "")}
            for item in doc.findall("Item"):
                nm = item.get("Name") or ""
                if nm == "GPL":  # semicolon separated
                    rec["gpl"] = (item.text or "").strip()
                elif nm in {
                    "title",
                    "summary",
                    "gdsType",
                    "Accession",
                    "taxon",
                    "n_samples",
                    "PDAT",
                    "entryType",
                    "GSE",
                    "PubMedIds",
                }:
                    rec[nm] = (item.text or "").strip()
                if nm == "PubMedIds":
                    pmids = [x.text for x in item.findall("int") if x.text]
                    if pmids:
                        rec["PubMedIds"] = ";".join(pmids)
            out.append(rec)
    return out


def search_uids_for_accessions(accessions: list[str]) -> dict[str, str]:
    """Resolve GSE accessions -> gds UIDs so seeds can join the same table."""
    mapping: dict[str, str] = {}
    term = " OR ".join(f"{a}[Accession]" for a in accessions)
    dest = CACHE_DIR / "eutils" / "esearch_seeds.xml"
    eutils("esearch.fcgi", {"db": "gds", "retmax": "200", "term": term}, dest)
    uids = [e.text for e in ET.parse(dest).getroot().findall(".//IdList/Id") if e.text]
    for rec in esummary(uids, "seeds"):
        acc = rec.get("Accession", "")
        if acc.startswith("GSE"):
            mapping[acc] = rec["uid"]
    return mapping


def platform_titles(gpls: list[str]) -> dict[str, str]:
    """Fetch platform titles (cheap, and tells array vs sequencer at a glance)."""
    titles: dict[str, str] = {}
    if not gpls:
        return titles
    term = " OR ".join(f"{g}[Accession]" for g in sorted(set(gpls)))
    dest = CACHE_DIR / "eutils" / "esearch_gpl.xml"
    eutils("esearch.fcgi", {"db": "gds", "retmax": "500", "term": term}, dest)
    uids = [e.text for e in ET.parse(dest).getroot().findall(".//IdList/Id") if e.text]
    for rec in esummary(uids, "gpl"):
        acc = rec.get("Accession", "")
        if acc.startswith("GPL"):
            titles[acc] = rec.get("title", "")
    return titles


def main() -> int:
    ensure_dirs()
    write_json(RESULTS_DIR / "geo_search_queries.json", {"queries": QUERIES, "seed_gse": SEED_GSE})

    hits: dict[str, dict] = {}
    query_hits: dict[str, list[str]] = {}
    for name, term in QUERIES.items():
        uids = esearch(name, term)
        recs = esummary(uids, name)
        accs = []
        for rec in recs:
            acc = rec.get("Accession", "")
            if not acc.startswith("GSE"):
                continue
            accs.append(acc)
            hits.setdefault(acc, rec)
            hits[acc].setdefault("found_by", set())
            hits[acc]["found_by"].add(name)
        query_hits[name] = sorted(set(accs))
        print(f"[query] {name}: {len(uids)} uids -> {len(set(accs))} GSE", flush=True)

    seed_map = search_uids_for_accessions(SEED_GSE)
    for rec in esummary(sorted(seed_map.values()), "seedsum"):
        acc = rec.get("Accession", "")
        if not acc.startswith("GSE"):
            continue
        hits.setdefault(acc, rec)
        hits[acc].setdefault("found_by", set())
        hits[acc]["found_by"].add("seed")
    missing_seeds = [g for g in SEED_GSE if g not in hits]
    if missing_seeds:
        print(f"[warn] seeds not resolvable via Entrez: {missing_seeds}", flush=True)

    gpls = [g for rec in hits.values() for g in (rec.get("gpl", "") or "").split(";") if g]
    gpl_title = platform_titles([f"GPL{g}" if not g.startswith("GPL") else g for g in gpls])

    rows = []
    for acc, rec in sorted(hits.items(), key=lambda kv: int(kv[0][3:])):
        text = " ".join([rec.get("title", ""), rec.get("summary", "")])
        gdst = rec.get("gdsType", "")
        gpl_list = [
            (g if g.startswith("GPL") else f"GPL{g}")
            for g in (rec.get("gpl", "") or "").split(";")
            if g
        ]
        is_array = bool(ARRAY_TYPE_RE.search(gdst)) and not bool(SEQ_TYPE_RE.search(gdst))
        rows.append(
            {
                "gse": acc,
                "n_samples": rec.get("n_samples", ""),
                "pdat": rec.get("PDAT", ""),
                "gds_type": gdst,
                "is_array_type": int(is_array),
                "lung_mention": int(bool(LUNG_RE.search(text))),
                "ici_mention": int(bool(ICI_RE.search(text))),
                "response_mention": int(bool(RESP_RE.search(text))),
                "platforms": ";".join(gpl_list),
                "platform_titles": " | ".join(gpl_title.get(g, "") for g in gpl_list),
                "pubmed": rec.get("PubMedIds", ""),
                "found_by": ",".join(sorted(rec.get("found_by", set()))),
                "title": rec.get("title", "").replace("\t", " "),
                "summary": rec.get("summary", "").replace("\t", " ")[:1500],
            }
        )

    hits_path = RESULTS_DIR / "geo_search_hits.tsv"
    cols = list(rows[0].keys())
    with hits_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"[write] {hits_path} ({len(rows)} rows)")

    cand = [
        r
        for r in rows
        if r["is_array_type"] and r["lung_mention"] and r["ici_mention"] and r["response_mention"]
    ]
    cand_path = RESULTS_DIR / "geo_candidates.tsv"
    with cand_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(cand)
    print(f"[write] {cand_path} ({len(cand)} candidate series)")
    for r in cand:
        print(f"  {r['gse']}  n={r['n_samples']:>4}  {r['pdat']}  {r['title'][:96]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
