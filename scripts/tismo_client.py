"""Public TISMO API client (no login).

tismo.cistrome.org redirects to tismo.pku-genomics.org. The Gene module uses
/tismo for metadata and /rtismo for expression table downloads.
"""

from __future__ import annotations

import csv
import io
import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

META_BASE = "https://tismo.pku-genomics.org/tismo"
R_BASE = "https://tismo.pku-genomics.org/rtismo"
USER_AGENT = "tismo-icb-recompute/1.0 (public TISMO API)"

ICB_TREATMENTS = [
    "antiCTLA4",
    "antiCTLA4&antiPD1",
    "antiCTLA4&antiPDL1",
    "antiPD1",
    "antiPDL1",
    "antiPDL2",
]


def post(url: str, fields: dict, timeout: int = 180, retries: int = 5) -> bytes:
    body = urlencode(fields).encode()
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = Request(
                url,
                data=body,
                headers={
                    "Content-Type": "application/x-www-form-urlencoded",
                    "User-Agent": USER_AGENT,
                    "Origin": "https://tismo.pku-genomics.org",
                    "Referer": "https://tismo.pku-genomics.org/",
                },
            )
            with urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (HTTPError, URLError, TimeoutError, OSError) as err:
            last = err
            time.sleep(2 ** attempt)
    raise RuntimeError(f"POST {url} failed after {retries} attempts: {last}")


def post_json(url: str, fields: dict | None = None) -> dict:
    payload = json.loads(post(url, fields or {}).decode())
    if payload.get("status") not in (None, 600200):
        raise RuntimeError(f"{url} status={payload.get('status')}: {payload.get('msg')}")
    return payload


def names(endpoint: str, extra: dict | None = None) -> list[str]:
    data = post_json(f"{META_BASE}{endpoint}", extra)["data"]
    return [row["name"] for row in data if row["name"] != "All"]


def get_metadata(kind: str) -> list[dict]:
    return post_json(f"{META_BASE}/metaData/{kind}", {"page": "1", "limit": "100000"})["data"]


def download_vivo_expression(gene: str, treatments: list[str], models: list[str]) -> str:
    raw = post(
        f"{R_BASE}/gene/downVivoExprn",
        {
            "filename": "vivo.csv",
            "type": "3",
            "gene": gene,
            "icbList": json.dumps(treatments),
            "tumorList": json.dumps(models),
        },
    )
    text = raw.decode("utf-8", errors="replace")
    if not text.startswith("Samples"):
        raise RuntimeError(f"unexpected vivo payload for {gene}: {text[:200]}")
    return text


def download_vitro_expression(gene: str, cytokines: list[str], models: list[str]) -> str:
    raw = post(
        f"{R_BASE}/gene/downVitroExprn",
        {
            "filename": "vitro.csv",
            "type": "3",
            "gene": gene,
            "icbList": json.dumps(cytokines),
            "tumorList": json.dumps(models),
            "cytokineList": json.dumps(cytokines),
            "cellList": json.dumps(models),
        },
    )
    text = raw.decode("utf-8", errors="replace")
    if not text.startswith("Samples"):
        raise RuntimeError(f"unexpected vitro payload for {gene}: {text[:200]}")
    return text


def read_expression_csv(text_or_path, gene: str | None = None):
    """Parse TISMO CSV, dropping rare malformed rows."""
    import pandas as pd

    if isinstance(text_or_path, Path):
        raw = text_or_path.read_text(encoding="utf-8", errors="replace")
    else:
        raw = text_or_path
    reader = csv.reader(io.StringIO(raw))
    header = next(reader)
    good = []
    n_bad = 0
    for row in reader:
        if len(row) == len(header):
            good.append(row)
        else:
            n_bad += 1
    df = pd.DataFrame(good, columns=header)
    if gene is not None and "geneID" in df.columns:
        df = df[df["geneID"] == gene]
    for col in ("value", "Baseline", "pvalue", "count", "num"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["value", "Baseline"])
    df.attrs["dropped_malformed_rows"] = n_bad
    return df.reset_index(drop=True)
