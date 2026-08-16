"""Minimal client for the public TISMO web API.

tismo.cistrome.org 301-redirects to tismo.pku-genomics.org. The site is a Vue SPA
whose "Gene" module calls two backends: `/tismo` for metadata/vocabulary lookups
and `/rtismo` for the R-backed expression endpoints. Both accept
form-urlencoded POSTs and require no authentication.
"""

from __future__ import annotations

import csv
import io
import json
import re
import tarfile
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

META_BASE = "https://tismo.pku-genomics.org/tismo"
R_BASE = "https://tismo.pku-genomics.org/rtismo"

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "align_tismo"
DATA_DIR = RESULTS_DIR / "data"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"

# The six ICB regimens exposed by the in vivo Gene module.
ICB_TREATMENTS = [
    "antiCTLA4",
    "antiCTLA4&antiPD1",
    "antiCTLA4&antiPDL1",
    "antiPD1",
    "antiPDL1",
    "antiPDL2",
]

_USER_AGENT = "align-tismo-replication/1.0 (+public TISMO API client)"


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
    """Tumor models (cell lines) that have in vivo samples under `treatments`."""
    joined = ",".join(f'"{t}"' for t in treatments)
    data = _post_json(f"{META_BASE}/gene/getVivoCohort", {"treatment": joined})["data"]
    return [row["name"] for row in data if row["name"] != "All"]


def get_gene_list() -> list[str]:
    data = _post_json(f"{META_BASE}/gene/getGene", {})["data"]
    return [row["name"] for row in data]


def get_metadata(kind: str) -> list[dict]:
    """kind is one of vivoMeta, vitroMeta, cellLineMeta."""
    return _post_json(f"{META_BASE}/metaData/{kind}", {"page": "1", "limit": "100000"})["data"]


def download_vivo_expression(gene: str, treatments: list[str], models: list[str]) -> str:
    """Per-sample in vivo expression for one gene, as returned by the Gene module.

    Columns: Samples, geneID, value, cell_line, Responder, Baseline, GSE_ID,
    Mouse_treatment, pvalue, label, label2, num, count.

    `cell_line` is TISMO's cohort label (cellline_study_condition_regimen(n=N)),
    `Baseline` is 1 for the study's control arm and 0 for the ICB-treated arm,
    and `pvalue` is TISMO's precomputed within-cohort DESeq2 statistic.
    """
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
    if text.lstrip().startswith("{"):
        raise RuntimeError(f"expression request for {gene} returned an error: {text[:200]}")
    if not text.startswith("Samples,"):
        raise RuntimeError(f"unexpected payload for {gene}: {text[:200]}")
    return text


EXPECTED_COLUMNS = [
    "Samples", "geneID", "value", "cell_line", "Responder", "Baseline",
    "GSE_ID", "Mouse_treatment", "pvalue", "label", "label2", "num", "count",
]

_COHORT_N_SUFFIX = re.compile(r"\(n=\d+\)$")


def cohort_key(label: str) -> str:
    """Stable cohort identifier.

    TISMO appends `(n=N)` to each cohort label, but N is recomputed per gene from
    the rows that survive its own filtering, so the same cohort appears as
    `...(n=10)` for one gene and `...(n=11)` for another. Strip it so cohorts can
    be compared across genes.
    """
    return _COHORT_N_SUFFIX.sub("", str(label)).rstrip()


def load_expression_csv(path, gene: str | None = None):
    """Read a downloaded expression CSV, tolerating TISMO's occasional bad rows.

    At least one served row is malformed (an extra field plus a `geneID` from a
    different gene), so rows with the wrong field count or the wrong gene are
    dropped and counted rather than allowed to poison the table.
    """
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
    df.attrs["dropped_malformed_rows"] = dropped_shape
    df.attrs["dropped_wrong_gene_rows"] = dropped_gene
    return df.reset_index(drop=True)


NULL_GENES_DIR = DATA_DIR / "null_genes"
NULL_GENES_ARCHIVE = DATA_DIR / "null_genes.tar.gz"


def pack_null_genes() -> None:
    """Archive the ~500 raw null-panel CSVs.

    Individually they are ~40 MB; the per-gene files repeat the same sample and
    cohort labels, so a single tarball compresses roughly 8x. The archive is what
    gets committed, the loose directory is a working cache.
    """
    if not NULL_GENES_DIR.exists():
        return
    with tarfile.open(NULL_GENES_ARCHIVE, "w:gz") as tar:
        for path in sorted(NULL_GENES_DIR.glob("*.csv")):
            tar.add(path, arcname=f"null_genes/{path.name}")
    size = NULL_GENES_ARCHIVE.stat().st_size
    print(f"  packed null panel -> {NULL_GENES_ARCHIVE.name} ({size / 1e6:.1f} MB)")


def ensure_null_genes() -> Path:
    """Return the null-panel directory, unpacking the archive if needed."""
    if NULL_GENES_DIR.exists() and any(NULL_GENES_DIR.glob("*.csv")):
        return NULL_GENES_DIR
    if NULL_GENES_ARCHIVE.exists():
        with tarfile.open(NULL_GENES_ARCHIVE, "r:gz") as tar:
            tar.extractall(DATA_DIR)
        print(f"  unpacked {NULL_GENES_ARCHIVE.name}")
    return NULL_GENES_DIR


def read_csv_text(text: str):
    import pandas as pd

    return pd.read_csv(io.StringIO(text))
