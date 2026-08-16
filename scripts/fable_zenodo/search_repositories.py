#!/usr/bin/env python3
"""Search Zenodo, figshare, and OSF for open lung ICI / TACSTD2 / CLDN4 matrices.

Outputs a consolidated JSON of candidate records with access/license metadata to
results/fable_zenodo/search_results.json and a human-readable summary to
notes/fable_zenodo/search_summary.md.

Only the public search APIs are used; nothing is downloaded here.
"""
import json
import time
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "fable_zenodo"
NOTES = ROOT / "notes" / "fable_zenodo"
RESULTS.mkdir(parents=True, exist_ok=True)
NOTES.mkdir(parents=True, exist_ok=True)

UA = "fable-zenodo-search/1.0 (research; contact: cloud-agent)"

QUERIES = [
    "lung immune checkpoint inhibitor single cell",
    "lung cancer ICI single cell RNA",
    "TACSTD2 lung",
    "CLDN4 lung",
    "TACSTD2 expression matrix",
    "CLDN4 expression matrix",
    "NSCLC immunotherapy single-cell matrix",
    "lung adenocarcinoma immune checkpoint scRNA-seq",
]


def http_get_json(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode("utf-8", "replace"))
        except Exception as e:  # noqa
            last = e
            time.sleep(2 ** i)
    print(f"  ! failed {url}: {last}", file=sys.stderr)
    return None


def search_zenodo(query, size=20):
    out = []
    q = urllib.parse.quote(query)
    url = f"https://zenodo.org/api/records?q={q}&size={size}&sort=bestmatch"
    data = http_get_json(url)
    if not data:
        return out
    for hit in data.get("hits", {}).get("hits", []):
        meta = hit.get("metadata", {})
        files = hit.get("files", [])
        access = meta.get("access_right") or hit.get("access", {}).get("record") or "unknown"
        out.append({
            "source": "zenodo",
            "query": query,
            "id": hit.get("id"),
            "title": meta.get("title"),
            "doi": hit.get("doi") or meta.get("doi"),
            "doi_url": hit.get("doi_url") or hit.get("links", {}).get("doi"),
            "access_right": access,
            "license": (meta.get("license") or {}).get("id") if isinstance(meta.get("license"), dict) else meta.get("license"),
            "publication_date": meta.get("publication_date"),
            "resource_type": (meta.get("resource_type") or {}).get("title"),
            "files": [
                {
                    "key": f.get("key"),
                    "size": f.get("size"),
                    "link": (f.get("links") or {}).get("self") or (f.get("links") or {}).get("download"),
                }
                for f in files
            ],
            "num_files": len(files),
            "html": hit.get("links", {}).get("self_html") or hit.get("links", {}).get("html"),
        })
    return out


def search_figshare(query, size=20):
    out = []
    url = "https://api.figshare.com/v2/articles/search"
    body = json.dumps({"search_for": query, "page_size": size}).encode()
    try:
        req = urllib.request.Request(url, data=body, headers={
            "User-Agent": UA, "Content-Type": "application/json", "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            arts = json.loads(r.read().decode("utf-8", "replace"))
    except Exception as e:  # noqa
        print(f"  ! figshare search failed: {e}", file=sys.stderr)
        return out
    for a in arts[:size]:
        detail = http_get_json(f"https://api.figshare.com/v2/articles/{a.get('id')}")
        time.sleep(0.2)
        files = (detail or {}).get("files", []) if detail else []
        out.append({
            "source": "figshare",
            "query": query,
            "id": a.get("id"),
            "title": a.get("title"),
            "doi": a.get("doi") or (detail or {}).get("doi"),
            "doi_url": (detail or {}).get("url"),
            "access_right": "open" if (detail or {}).get("is_public", True) else "restricted",
            "license": ((detail or {}).get("license") or {}).get("name"),
            "publication_date": (detail or {}).get("published_date"),
            "resource_type": (detail or {}).get("defined_type_name"),
            "files": [
                {"key": f.get("name"), "size": f.get("size"), "link": f.get("download_url")}
                for f in files
            ],
            "num_files": len(files),
            "html": a.get("url_public_html") or (detail or {}).get("figshare_url"),
        })
    return out


def search_osf(query, size=20):
    out = []
    q = urllib.parse.quote(query)
    url = f"https://api.osf.io/v2/nodes/?filter[title]={q}&page[size]={size}"
    data = http_get_json(url)
    if not data:
        return out
    for n in data.get("data", []):
        attr = n.get("attributes", {})
        out.append({
            "source": "osf",
            "query": query,
            "id": n.get("id"),
            "title": attr.get("title"),
            "doi": None,
            "doi_url": (n.get("links") or {}).get("html"),
            "access_right": "public" if attr.get("public") else "restricted",
            "license": None,
            "publication_date": attr.get("date_created"),
            "resource_type": attr.get("category"),
            "files": [],
            "num_files": None,
            "html": (n.get("links") or {}).get("html"),
        })
    return out


def main():
    all_records = []
    for q in QUERIES:
        print(f"[zenodo] {q}")
        all_records += search_zenodo(q)
        time.sleep(0.5)
        print(f"[figshare] {q}")
        all_records += search_figshare(q)
        time.sleep(0.5)
        print(f"[osf] {q}")
        all_records += search_osf(q)
        time.sleep(0.5)

    # Deduplicate by (source, id)
    seen = {}
    for r in all_records:
        key = (r["source"], str(r["id"]))
        if key not in seen:
            seen[key] = r
        else:
            # merge query attribution
            prev = seen[key]
            if isinstance(prev.get("query"), str):
                prev["query"] = [prev["query"]]
            if r["query"] not in prev["query"]:
                prev["query"].append(r["query"])
    records = list(seen.values())

    (RESULTS / "search_results.json").write_text(json.dumps(records, indent=2, ensure_ascii=False))
    print(f"\nTotal unique records: {len(records)}")
    print(f"Wrote {RESULTS / 'search_results.json'}")


if __name__ == "__main__":
    main()
