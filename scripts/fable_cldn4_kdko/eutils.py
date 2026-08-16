"""Minimal NCBI E-utilities client (no Entrez Direct needed).

Used to exhaustively search GEO DataSets (gds) and SRA for CLDN4
perturbation studies, and to fetch summaries for verification.
"""
import time
import json
import sys
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# Polite: identify tool + email per NCBI usage policy.
TOOL = "fable_cldn4_kdko"
EMAIL = "cursoragent@cursor.com"


def _get(url, retries=5):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": TOOL})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed after {retries}: {url}\n{last}")


def esearch(db, term, retmax=500):
    q = urllib.parse.urlencode({
        "db": db, "term": term, "retmax": retmax,
        "retmode": "json", "tool": TOOL, "email": EMAIL,
    })
    time.sleep(0.4)
    data = json.loads(_get(f"{BASE}/esearch.fcgi?{q}"))
    return data.get("esearchresult", {})


def esummary(db, ids):
    if not ids:
        return {}
    out = {}
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        q = urllib.parse.urlencode({
            "db": db, "id": ",".join(chunk),
            "retmode": "json", "tool": TOOL, "email": EMAIL,
        })
        time.sleep(0.4)
        data = json.loads(_get(f"{BASE}/esummary.fcgi?{q}"))
        res = data.get("result", {})
        for uid in res.get("uids", []):
            out[uid] = res[uid]
    return out


if __name__ == "__main__":
    # quick CLI: python eutils.py <db> "<term>"
    db, term = sys.argv[1], sys.argv[2]
    r = esearch(db, term)
    print("count:", r.get("count"), "ids:", len(r.get("idlist", [])))
    print(json.dumps(r.get("idlist", []), indent=1))
