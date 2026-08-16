"""Shared NCBI E-utilities helpers (no API key required, rate-limited).

Used by the GEO 2015-2018 human lung ICI/PD-1/PD-L1/CTLA-4 mining pipeline.
All network access is read-only against public NCBI endpoints.
"""
import time
import json
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# Be polite: <= 3 requests/sec without an API key.
_MIN_INTERVAL = 0.4
_last = [0.0]

TOOL = "fable_geo_2015_2018"
EMAIL = "fable-cloud-agent@example.com"


def _throttle():
    dt = time.time() - _last[0]
    if dt < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - dt)
    _last[0] = time.time()


def _get(url, tries=5):
    last_err = None
    for i in range(tries):
        _throttle()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "fable-geo/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last_err = e
            time.sleep(2 ** i)
    raise RuntimeError(f"request failed after {tries} tries: {url}\n{last_err}")


def esearch(db, term, retmax=1000):
    """Return list of UIDs for a term."""
    ids = []
    retstart = 0
    while True:
        params = {
            "db": db,
            "term": term,
            "retmode": "json",
            "retmax": retmax,
            "retstart": retstart,
            "tool": TOOL,
            "email": EMAIL,
        }
        url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)
        data = json.loads(_get(url))
        res = data.get("esearchresult", {})
        batch = res.get("idlist", [])
        ids.extend(batch)
        total = int(res.get("count", "0"))
        retstart += retmax
        if retstart >= total or not batch:
            break
    return ids


def esummary(db, uids):
    """Return dict uid -> summary docsum (JSON esummary)."""
    out = {}
    for i in range(0, len(uids), 300):
        chunk = uids[i:i + 300]
        params = {
            "db": db,
            "id": ",".join(chunk),
            "retmode": "json",
            "tool": TOOL,
            "email": EMAIL,
        }
        url = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(params)
        data = json.loads(_get(url))
        result = data.get("result", {})
        for uid in result.get("uids", []):
            out[uid] = result[uid]
    return out
