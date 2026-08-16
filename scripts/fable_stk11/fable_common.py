"""Shared utilities for the fable_stk11 slice.

Analysis theme
--------------
Public lung ICI (immune-checkpoint-inhibitor) datasets carrying STK11/KEAP1/KRAS
genotype labels are used to relate the antibody-drug-conjugate (ADC) target genes
TACSTD2 (TROP2) and CLDN4 (Claudin-4) to (a) genotype and (b) ICI response.

All paths are confined to results/fable_stk11/ so the slice never writes outside
its allowed output tree.
"""
from __future__ import annotations

import io
import json
import os
import time
import urllib.request
import urllib.error

# ----------------------------------------------------------------------------
# Paths (all outputs live under results/fable_stk11/)
# ----------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RESULTS = os.path.join(REPO, "results", "fable_stk11")
DATA = os.path.join(RESULTS, "data")
RAW = os.path.join(DATA, "raw")
PROC = os.path.join(DATA, "processed")
FIG = os.path.join(RESULTS, "figures")
TAB = os.path.join(RESULTS, "tables")
NOTES = os.path.join(REPO, "notes", "fable_stk11")

for _d in (RAW, PROC, FIG, TAB, NOTES):
    os.makedirs(_d, exist_ok=True)

# ----------------------------------------------------------------------------
# Genes of interest
# ----------------------------------------------------------------------------
GENES = {
    # genotype / drivers
    "STK11": 6794,
    "KEAP1": 9817,
    "KRAS": 3845,
    "TP53": 7157,
    "EGFR": 1956,
    # ADC targets of interest
    "TACSTD2": 4070,   # TROP2
    "CLDN4": 1364,     # Claudin-4
    # immune context
    "CD274": 29126,    # PD-L1
    "CD8A": 925,
}
ENTREZ_TO_SYMBOL = {v: k for k, v in GENES.items()}

GENOTYPE_GENES = ["STK11", "KEAP1", "KRAS"]
TARGET_GENES = ["TACSTD2", "CLDN4"]
CONTEXT_GENES = ["CD274", "CD8A", "TP53", "EGFR"]

# Ensembl gene ids (for GEO GSE135222, versionless matched on prefix)
ENSEMBL = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
    "CD274": "ENSG00000120217",
    "CD8A": "ENSG00000153563",
    "STK11": "ENSG00000118046",
    "KEAP1": "ENSG00000079999",
    "KRAS": "ENSG00000133703",
}

CBIO = "https://www.cbioportal.org/api"


# ----------------------------------------------------------------------------
# HTTP helpers with retry/backoff
# ----------------------------------------------------------------------------
def _request(url: str, data: bytes | None = None, headers: dict | None = None,
             retries: int = 5, timeout: int = 120) -> bytes:
    hdrs = {"Accept": "application/json"}
    if headers:
        hdrs.update(headers)
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers=hdrs,
                                         method="POST" if data is not None else "GET")
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last = e
            wait = 2 ** i
            print(f"  [retry {i+1}/{retries}] {url} -> {e}; sleeping {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"Request failed after {retries} tries: {url}\n{last}")


def cbio_get(path: str, **params) -> object:
    q = ""
    if params:
        q = "?" + "&".join(f"{k}={v}" for k, v in params.items())
    return json.loads(_request(f"{CBIO}{path}{q}"))


def cbio_post(path: str, payload: object) -> object:
    body = json.dumps(payload).encode()
    return json.loads(_request(f"{CBIO}{path}", data=body,
                               headers={"Content-Type": "application/json"}))


def download(url: str, dest: str, timeout: int = 300) -> str:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  cached: {os.path.basename(dest)}")
        return dest
    print(f"  downloading: {url}")
    raw = _request(url, timeout=timeout)
    with open(dest, "wb") as f:
        f.write(raw)
    return dest


def save_json(obj: object, path: str) -> None:
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=str)


def dir_size_mb(path: str) -> float:
    total = 0
    for root, _dirs, files in os.walk(path):
        for fn in files:
            fp = os.path.join(root, fn)
            if os.path.exists(fp):
                total += os.path.getsize(fp)
    return total / (1024 * 1024)
