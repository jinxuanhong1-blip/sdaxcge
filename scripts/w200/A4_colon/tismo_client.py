"""Minimal public TISMO client (tismo.pku-genomics.org).

tismo.cistrome.org 301-redirects here. The Gene module uses `/tismo` for
metadata and `/rtismo` for expression CSVs. No authentication.
"""

from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

META_BASE = "https://tismo.pku-genomics.org/tismo"
R_BASE = "https://tismo.pku-genomics.org/rtismo"

REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = REPO_ROOT / "results" / "w200" / "A4_colon"
DATA_DIR = RESULTS_DIR / "data"

ICB_TREATMENTS = [
    "antiCTLA4",
    "antiCTLA4&antiPD1",
    "antiCTLA4&antiPDL1",
    "antiPD1",
    "antiPDL1",
    "antiPDL2",
]

_USER_AGENT = "w200-A4-colon/1.0 (+public TISMO API client)"
_COHORT_N_SUFFIX = re.compile(r"\(n=\d+\)$")

# TISMO cellLineMeta labels these four as colorectal carcinoma.
COLON_CANCER_TYPES = {"Colorectal carcinoma"}


def _post(url: str, fields: dict[str, str], timeout: int = 120, retries: int = 5) -> bytes:
    body = urlencode(fields).encode()
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            req = Request(
                url,
                data=body,
                headers={
                    "Content-Type": "application/x-www-form-urlencoded",
                    "User-Agent": _USER_AGENT,
                },
            )
            with urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (HTTPError, URLError, TimeoutError, OSError) as err:
            last_err = err
            time.sleep(2 ** attempt)
    raise RuntimeError(f"POST {url} failed after {retries} attempts: {last_err}")


def _post_json(url: str, fields: dict[str, str]) -> dict:
    payload = json.loads(_post(url, fields).decode())
    status = payload.get("status")
    if status is not None and status != 600200:
        raise RuntimeError(f"{url} returned status={status}: {payload.get('msg')}")
    return payload


def get_vivo_treatments() -> list[str]:
    data = _post_json(f"{META_BASE}/gene/getVivoTreatment", {})["data"]
    return [row["name"] for row in data if row["name"] != "All"]


def get_vivo_cohorts(treatments: list[str]) -> list[str]:
    joined = ",".join(f'"{t}"' for t in treatments)
    data = _post_json(f"{META_BASE}/gene/getVivoCohort", {"treatment": joined})["data"]
    return [row["name"] for row in data if row["name"] != "All"]


def get_metadata(kind: str) -> list[dict]:
    return _post_json(f"{META_BASE}/metaData/{kind}", {"page": "1", "limit": "100000"})["data"]


def download_vivo_expression(gene: str, treatments: list[str], models: list[str]) -> str:
    raw = _post(
        f"{R_BASE}/gene/downVivoExprn",
        {
            "filename": "genetreatment_vivo.csv",
            "type": "csv",
            "gene": gene,
            "icbList": json.dumps(treatments),
            "tumorList": json.dumps(models),
        },
    )
    text = raw.decode("utf-8", errors="replace")
    if text.lstrip().startswith("{") or not text.startswith("Samples,"):
        raise RuntimeError(f"unexpected payload for {gene}: {text[:200]}")
    return text


def cohort_key(label: str) -> str:
    """Strip TISMO's per-gene `(n=N)` suffix so cohorts match across genes."""
    return _COHORT_N_SUFFIX.sub("", str(label)).rstrip()


def load_expression_csv(path: Path, gene: str | None = None):
    import pandas as pd

    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        good, dropped_shape = [], 0
        for row in reader:
            if len(row) == len(header):
                good.append(row)
            else:
                dropped_shape += 1

    df = pd.DataFrame(good, columns=header)
    dropped_gene = 0
    if gene is not None and "geneID" in df.columns:
        keep = df["geneID"] == gene
        dropped_gene = int((~keep).sum())
        df = df[keep]

    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df["Baseline"] = pd.to_numeric(df["Baseline"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df = df.dropna(subset=["value", "Baseline"])
    df["cohort"] = df["cell_line"].map(cohort_key)
    df["model"] = df["cell_line"].map(lambda s: str(s).split("_")[0])
    df.attrs["dropped_malformed_rows"] = dropped_shape
    df.attrs["dropped_wrong_gene_rows"] = dropped_gene
    return df.reset_index(drop=True)


def load_cell_line_cancer_types(path: Path | None = None) -> dict[str, str]:
    p = path or (DATA_DIR / "cellLineMeta.json")
    rows = json.loads(p.read_text())
    return {r["cellLine"]: (r.get("cancerType") or "").strip() for r in rows}


def colon_models(cancer_types: dict[str, str]) -> list[str]:
    return sorted(m for m, ct in cancer_types.items() if ct in COLON_CANCER_TYPES)
