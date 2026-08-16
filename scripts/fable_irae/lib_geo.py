"""Small helpers for GEO access via NCBI (no external deps beyond stdlib)."""
import time
import sys
import urllib.request


def fetch_url(url, timeout=90, retries=4, binary=False):
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                data = r.read()
                return data if binary else data.decode("utf-8", errors="replace")
        except Exception as e:  # noqa
            last = e
            wait = 4 * (2 ** attempt)
            print(f"  retry {attempt+1}/{retries} after {wait}s: {e}", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"failed to fetch {url}: {last}")


def geo_text(acc, targ="self", view="brief"):
    """Fetch GEO record as plain text via acc.cgi."""
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
           f"acc={acc}&targ={targ}&form=text&view={view}")
    return fetch_url(url)


def series_suppl_url(acc, fname):
    stub = acc[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/suppl/{fname}"
