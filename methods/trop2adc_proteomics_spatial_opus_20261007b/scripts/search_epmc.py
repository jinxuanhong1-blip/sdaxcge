#!/usr/bin/env python3
"""Mine Europe PMC for TROP2-ADC perturbation datasets (GEO/PXD/etc.).

Strategy: query combinations of (TROP2-ADC drug terms) AND (omics/assay terms),
then for every matching article pull Europe PMC text-mined accession annotations
(GEO, PRIDE/ProteomeXchange, ArrayExpress/BioStudies). Emits a deduplicated table
of candidate accessions with the paper that referenced them.

No large downloads: only JSON metadata.
"""
import json, sys, time, urllib.parse, urllib.request, collections

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"

DRUG_TERMS = [
    "sacituzumab", "IMMU-132", "IMMU132", "Trodelvy",
    "datopotamab", "Dato-DXd", "TROPION",
    "SKB264", "sacituzumab tirumotecan", "sac-TMT", "MK-2870",
    "TROP2 ADC", "TROP2 antibody-drug conjugate", "anti-TROP2 antibody",
]
ASSAY_TERMS = [
    "proteomics", "proteome", "mass spectrometry",
    "spatial transcriptomics", "Visium", "CosMx", "Xenium", "GeoMx",
    "RNA-seq", "transcriptomic", "single-cell",
]


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.load(r)
        except Exception as e:
            if i == tries - 1:
                sys.stderr.write(f"FAIL {url}: {e}\n")
                return None
            time.sleep(2 * (i + 1))
    return None


def search(query, page_size=50):
    q = urllib.parse.quote(query)
    url = f"{EPMC}/search?query={q}&format=json&pageSize={page_size}&resultType=core"
    d = get(url)
    if not d:
        return []
    return d.get("resultList", {}).get("result", [])


def main():
    hits = {}  # key = (src,id) -> paper meta
    for drug in DRUG_TERMS:
        for assay in ASSAY_TERMS:
            query = f'("{drug}") AND ("{assay}")'
            res = search(query)
            for r in res:
                src = r.get("source"); pid = r.get("id")
                if not pid:
                    continue
                key = (src, pid)
                if key not in hits:
                    hits[key] = {
                        "source": src, "id": pid,
                        "pmid": r.get("pmid"), "doi": r.get("doi"),
                        "title": r.get("title", ""),
                        "journal": (r.get("journalInfo", {}) or {}).get("journal", {}).get("title", ""),
                        "year": r.get("pubYear", ""),
                        "matched": set(),
                    }
                hits[key]["matched"].add(f"{drug}+{assay}")
            time.sleep(0.15)
    sys.stderr.write(f"Collected {len(hits)} candidate papers from Europe PMC\n")

    # Pull text-mined accessions for each paper.
    rows = []
    acc_counter = collections.Counter()
    for (src, pid), meta in hits.items():
        url = f"{EPMC}/{src}/{pid}/textMinedTerms?type=Accession%20Numbers&format=json"
        d = get(url)
        accs = []
        if d:
            for sem in d.get("semanticTypeList", {}).get("semanticType", []):
                for t in sem.get("tmSummary", []):
                    accs.append((sem.get("name", ""), t.get("term", "")))
        meta["accessions"] = accs
        for db, term in accs:
            acc_counter[(db, term)] += 1
        rows.append(meta)
        time.sleep(0.1)

    out = {
        "n_papers": len(rows),
        "papers": [
            {**{k: v for k, v in m.items() if k != "matched"},
             "matched": sorted(m["matched"])}
            for m in rows
        ],
    }
    json.dump(out, open("data/epmc_candidates.json", "w"), indent=1)

    # Print accession summary
    print("=== Accessions referenced by candidate papers (db | accession | #papers) ===")
    for (db, term), c in acc_counter.most_common():
        print(f"{db} | {term} | {c}")
    print(f"\nTotal candidate papers: {len(rows)}")


if __name__ == "__main__":
    main()
