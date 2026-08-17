"""Download helpers, GEO parsers, and gene-panel extractors."""

from __future__ import annotations

import gzip
import hashlib
import shutil
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = Path("/tmp/multinichenet_cldn4")


def log(msg: str) -> None:
    print(msg, flush=True)


def sha256_head(path: Path, nbytes: int = 1_000_000) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        h.update(fh.read(nbytes))
    return h.hexdigest()[:16]


def download(url: str, dest: Path, min_bytes: int = 1, attempts: int = 4) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        log(f"HAVE {dest} ({dest.stat().st_size} bytes)")
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    wait = 4
    last_err: Exception | None = None
    for i in range(1, attempts + 1):
        try:
            log(f"GET {url} -> {dest} (try {i}/{attempts})")
            import urllib.request

            urllib.request.urlretrieve(url, tmp)
            if tmp.stat().st_size < min_bytes:
                raise RuntimeError(f"too small: {tmp.stat().st_size} < {min_bytes}")
            tmp.replace(dest)
            log(f"OK {dest} {dest.stat().st_size} bytes sha256_1MB={sha256_head(dest)}")
            return dest
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            log(f"fail try {i}: {exc}")
            if tmp.exists():
                tmp.unlink()
            if i < attempts:
                time.sleep(wait)
                wait *= 2
    raise RuntimeError(f"download failed {url}: {last_err}")


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    total = np.maximum(total.astype(np.float64), 1.0)
    return np.log1p(1e4 * umi.astype(np.float64) / total).astype(np.float32)


def stream_dense_umi_gz(
    path: Path, keep: set[str]
) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, int]:
    """One pass over a genes×cells TSV.gz: library totals + keep-gene rows."""
    log(f"[stream] {path} keep={len(keep)}")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"column mismatch for {gene}: {vals.size} != {n}")
            totals += vals
            if gene in keep and gene not in store:
                store[gene] = vals
            if n_genes % 5000 == 0:
                log(f"[stream] genes_seen={n_genes} kept={len(store)}")
    log(f"[stream] cells={n} genes_in_file={n_genes} genes_kept={len(store)}")
    return cells, store, totals, n_genes


def load_gse205335_rds(
    path: Path, keep: set[str]
) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    """Extract keep-genes from the public GSE205335 dgCMatrix RDS.gz."""
    import rdata
    from scipy import sparse

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = Path(tmp) / "GSE205335_Lung_IO_UMI_matrix.rds"
        log(f"decompress {path.name}")
        with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
            shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        log("read RDS")
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    log(f"build CSC {tuple(obj.Dim)}")
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    log("convert CSR for row extract")
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(keep):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel().astype(np.float32)
    del matrix
    log(f"extracted {len(extracted)} / {len(keep)} genes")
    return barcodes, extracted, library_umi


def parse_gse131907_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            if key not in fields:
                fields[key] = values
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def parse_gse205335_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def panel_to_frame(
    barcodes: np.ndarray,
    store: dict[str, np.ndarray],
    totals: np.ndarray,
    meta: pd.DataFrame,
) -> pd.DataFrame:
    """Bind UMI columns + ncount onto a cell metadata frame aligned to barcodes."""
    genes = sorted(store)
    umi = pd.DataFrame({g: store[g] for g in genes}, index=pd.Index(barcodes, name="barcode"))
    umi["ncount"] = totals
    out = meta.join(umi, how="inner")
    return out


def write_cohort_panel(df: pd.DataFrame, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tumor = df[df["is_tumor"]] if "is_tumor" in df.columns else df
    ntab = (
        tumor.groupby("patient", observed=True)
        .agg(
            n_cells=("patient", "size"),
            n_malignant=("is_malignant", "sum"),
            n_tnk=("is_tnk", "sum"),
        )
        .reset_index()
    )
    ntab["cohort"] = df["cohort"].iloc[0] if len(df) else ""
    ntab.to_csv(dest.with_name(dest.stem.replace("_panel", "") + "_patient_n.tsv"), sep="\t", index=False)
    keep = df[df["is_malignant"] | df["is_tnk"]].copy()
    keep.to_parquet(dest)
    log(f"wrote {dest} rows={len(keep)} cols={keep.shape[1]} bytes={dest.stat().st_size}")
