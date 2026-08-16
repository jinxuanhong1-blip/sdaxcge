"""GEO ICI cohort loaders (GSE126044, GSE135222) and a generic series-matrix parser.

Sample identifiers in a GEO count table and in the series matrix are almost
never in the same order and almost never spelled the same way. These loaders
join on a normalised sample ID, never on column position.
"""

from __future__ import annotations

import gzip
import io
import re
from pathlib import Path
from urllib.request import urlopen, Request

import pandas as pd

__all__ = [
    "download",
    "read_series_matrix",
    "load_gse126044",
    "load_gse135222",
    "GSE126044_COUNTS",
    "GSE126044_MATRIX",
    "GSE135222_EXPR",
    "GSE135222_MATRIX",
]

GSE126044_COUNTS = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/"
    "GSE126044_counts.txt.gz"
)
GSE126044_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/"
    "GSE126044_series_matrix.txt.gz"
)
GSE135222_EXPR = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/"
    "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
)
GSE135222_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/"
    "GSE135222_series_matrix.txt.gz"
)


def download(url: str, dest: str | Path, timeout: int = 120) -> Path:
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = Request(url, headers={"User-Agent": "bulkimmune/1.0"})
    with urlopen(req, timeout=timeout) as resp, open(dest, "wb") as handle:
        handle.write(resp.read())
    return dest


def _open_text(path: str | Path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, encoding="utf-8")


def read_series_matrix(path: str | Path) -> pd.DataFrame:
    """Parse a GEO series matrix into a samples x characteristics table."""
    collected: list[tuple[str, list[str]]] = []
    n_samples: int | None = None
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("!"):
                continue
            key, *rest = line.rstrip("\n").split("\t")
            key = key.lstrip("!")
            values = [v.strip().strip('"') for v in rest]
            if key == "Sample_title":
                n_samples = len(values)
            collected.append((key, values))
    if n_samples is None:
        raise ValueError(f"no Sample_title line in {path}")
    rows: dict[str, list[str]] = {}
    for key, values in collected:
        if len(values) != n_samples:
            continue
        if key in rows:
            n = sum(1 for k in rows if k.startswith(key))
            rows[f"{key}__{n}"] = values
        else:
            rows[key] = values
    frame = pd.DataFrame(rows)
    # Flatten characteristics "key: value" into columns
    extra = {}
    for col in list(frame.columns):
        if "characteristics" not in col.lower():
            continue
        parsed = frame[col].str.split(": ", n=1, expand=True)
        if parsed.shape[1] != 2:
            continue
        keys = parsed[0].dropna().unique()
        if len(keys) != 1:
            continue
        extra[keys[0].strip()] = parsed[1]
    out = pd.concat([frame, pd.DataFrame(extra)], axis=1)
    return out


def _norm_id(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", str(s)).upper()


def load_gse126044(
    counts_path: str | Path,
    matrix_path: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cho et al. Nat Med 2020 -- 16 NSCLC pre-anti-PD-1 biopsies.

    Returns (counts genes x samples, phenotype samples x covariates).
    Phenotype columns: sample_id, geo_accession, response, tissue_preservation
    (fresh vs FFPE -- a real batch covariate), title.
    """
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    pheno = read_series_matrix(matrix_path)
    title_col = "Sample_title" if "Sample_title" in pheno.columns else pheno.columns[0]
    pheno = pheno.copy()
    pheno["sample_id"] = (
        pheno[title_col].str.replace(r"^RNA-seq_", "", regex=True).str.strip()
    )
    pheno["geo_accession"] = pheno.get("Sample_geo_accession", pd.NA)
    # characteristics keys observed in this series
    resp_col = next((c for c in pheno.columns if c.lower().startswith("patient response")), None)
    pres_col = next((c for c in pheno.columns if c.lower() == "sample"), None)
    pheno["response"] = pheno[resp_col].str.lower().str.strip() if resp_col else pd.NA
    pheno["tissue_preservation"] = pheno[pres_col].str.lower().str.strip() if pres_col else pd.NA
    pheno = pheno.set_index("sample_id")
    # Join on normalised IDs -- column order in the count table is NOT the
    # series-matrix order (Dis_07 / Dis_10 / Dis_06 are shuffled).
    cmap = {_norm_id(c): c for c in counts.columns}
    pmap = {_norm_id(i): i for i in pheno.index}
    shared = [k for k in cmap if k in pmap]
    if len(shared) < 8:
        raise ValueError(
            f"GSE126044: only {len(shared)} samples matched between counts and series matrix"
        )
    counts = counts[ [cmap[k] for k in shared] ]
    pheno = pheno.loc[ [pmap[k] for k in shared] ]
    pheno.index = [cmap[k] for k in shared]
    counts.columns = pheno.index
    return counts, pheno[["geo_accession", "response", "tissue_preservation"]]


def load_gse135222(
    expr_path: str | Path,
    matrix_path: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Jung et al. Nat Commun 2019 -- 27 advanced NSCLC anti-PD-1/PD-L1, FPKM + PFS.

    The supplementary matrix is already in FPKM-like units (columns do not sum
    to 1e6). Identifiers are versioned Ensembl IDs.
    """
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    pheno = read_series_matrix(matrix_path)
    title_col = "Sample_title" if "Sample_title" in pheno.columns else pheno.columns[0]
    pheno = pheno.copy()
    pheno["sample_id"] = pheno[title_col].str.replace(" ", "", regex=False).str.strip()
    pheno["geo_accession"] = pheno.get("Sample_geo_accession", pd.NA)
    def _pick(*cands):
        for c in pheno.columns:
            cl = c.lower()
            if any(x in cl for x in cands):
                return c
        return None

    pfs_event = _pick("progression-free survival", "pfs)")
    pfs_time = _pick("pfs.time", "pfs time")
    gender = _pick("gender", "sex")
    age = _pick("age")
    pheno["pfs_event"] = pd.to_numeric(pheno[pfs_event], errors="coerce") if pfs_event else pd.NA
    pheno["pfs_time"] = pd.to_numeric(pheno[pfs_time], errors="coerce") if pfs_time else pd.NA
    pheno["sex"] = pheno[gender].str.lower().str.strip() if gender else pd.NA
    pheno["age"] = pd.to_numeric(pheno[age], errors="coerce") if age else pd.NA
    pheno = pheno.set_index("sample_id")

    cmap = {_norm_id(c): c for c in expr.columns}
    pmap = {_norm_id(i): i for i in pheno.index}
    shared = [k for k in cmap if k in pmap]
    if len(shared) < 10:
        raise ValueError(
            f"GSE135222: only {len(shared)} samples matched between expression and series matrix"
        )
    expr = expr[[cmap[k] for k in shared]]
    pheno = pheno.loc[[pmap[k] for k in shared]]
    pheno.index = [cmap[k] for k in shared]
    expr.columns = pheno.index
    return expr, pheno[["geo_accession", "pfs_event", "pfs_time", "sex", "age"]]
