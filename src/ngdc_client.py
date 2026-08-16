"""Cached, polite HTTP client for NGDC/CNCB endpoints.

Everything the hunt does over the network goes through here so that (a) reruns are
cheap and reproducible, and (b) we never hammer ngdc.cncb.ac.cn.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterable, Optional, Sequence

import requests

LOG = logging.getLogger("ngdc")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_DIR = os.path.join(ROOT, "data", "cache")
DOWNLOAD_DIR = os.path.join(ROOT, "data", "raw")

USER_AGENT = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/126.0.0.0 Safari/537.36"
)

# ngdc.cncb.ac.cn is in Beijing and occasionally drops connections; be patient.
DEFAULT_TIMEOUT = (30, 300)
MAX_ATTEMPTS = 5

_throttle_lock = threading.Lock()
_last_request_at = [0.0]
MIN_INTERVAL_S = 0.12


def _throttle() -> None:
    with _throttle_lock:
        wait = MIN_INTERVAL_S - (time.time() - _last_request_at[0])
        if wait > 0:
            time.sleep(wait)
        _last_request_at[0] = time.time()


def _cache_path(key: str) -> str:
    h = hashlib.sha1(key.encode("utf-8")).hexdigest()
    return os.path.join(CACHE_DIR, f"{h}.cache")


def _session() -> requests.Session:
    s = requests.Session()
    s.headers.update(
        {
            "User-Agent": USER_AGENT,
            "Accept": "text/html,application/json,application/xhtml+xml,*/*",
            "Accept-Language": "en-US,en;q=0.9",
            "Connection": "keep-alive",
        }
    )
    return s


_thread_local = threading.local()


def session() -> requests.Session:
    s = getattr(_thread_local, "session", None)
    if s is None:
        s = _session()
        _thread_local.session = s
    return s


def fetch(
    url: str,
    *,
    method: str = "GET",
    data: Optional[dict] = None,
    referer: Optional[str] = None,
    use_cache: bool = True,
    allow_error: bool = True,
) -> dict:
    """Fetch a URL, returning {'ok','status','url','text'} and caching the result."""
    key = f"{method} {url} {json.dumps(data, sort_keys=True) if data else ''}"
    path = _cache_path(key)
    if use_cache and os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)

    os.makedirs(CACHE_DIR, exist_ok=True)
    last_err = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            _throttle()
            headers = {"Referer": referer} if referer else None
            resp = session().request(
                method, url, data=data, headers=headers, timeout=DEFAULT_TIMEOUT
            )
            out = {
                "ok": resp.status_code == 200,
                "status": resp.status_code,
                "url": resp.url,
                "text": resp.text,
                "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
            if resp.status_code in (429, 500, 502, 503, 504) and attempt < MAX_ATTEMPTS:
                raise requests.RequestException(f"HTTP {resp.status_code}")
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(out, fh)
            return out
        except Exception as exc:  # noqa: BLE001 - network layer, retry everything
            last_err = exc
            if attempt < MAX_ATTEMPTS:
                time.sleep(min(2 ** attempt, 30))

    LOG.warning("fetch failed after %d attempts: %s (%s)", MAX_ATTEMPTS, url, last_err)
    out = {
        "ok": False,
        "status": 0,
        "url": url,
        "text": "",
        "error": str(last_err),
        "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if not allow_error:
        raise RuntimeError(f"fetch failed: {url}: {last_err}")
    return out


def fetch_many(
    urls: Sequence[str],
    *,
    workers: int = 8,
    referer: Optional[str] = None,
    progress_every: int = 250,
    on_result: Optional[Callable[[str, dict], None]] = None,
) -> dict:
    """Fetch many URLs concurrently. Returns {url: result} unless on_result is given."""
    results: dict = {}
    done = [0]
    lock = threading.Lock()

    def work(u: str):
        r = fetch(u, referer=referer)
        with lock:
            done[0] += 1
            if progress_every and done[0] % progress_every == 0:
                LOG.info("  fetched %d/%d", done[0], len(urls))
            if on_result is not None:
                on_result(u, r)
            else:
                results[u] = r

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(work, urls))
    return results


def download(url: str, dest: str, *, min_bytes: int = 1) -> bool:
    """Stream a file to disk (skipped if already present and non-trivial)."""
    if os.path.exists(dest) and os.path.getsize(dest) >= min_bytes:
        return True
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            _throttle()
            with session().get(url, stream=True, timeout=DEFAULT_TIMEOUT) as resp:
                if resp.status_code != 200:
                    LOG.warning("download HTTP %s: %s", resp.status_code, url)
                    if resp.status_code < 500:
                        return False
                    raise requests.RequestException(f"HTTP {resp.status_code}")
                with open(tmp, "wb") as fh:
                    for chunk in resp.iter_content(1 << 20):
                        if chunk:
                            fh.write(chunk)
            os.replace(tmp, dest)
            return True
        except Exception as exc:  # noqa: BLE001
            LOG.warning("download attempt %d failed (%s): %s", attempt, url, exc)
            if attempt < MAX_ATTEMPTS:
                time.sleep(min(2 ** attempt, 30))
    return False


def setup_logging(level: int = logging.INFO) -> None:
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
