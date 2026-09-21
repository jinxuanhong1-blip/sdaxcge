"""Stream TCGA GDC STAR log2(TPM+1) from the UCSC Xena GDC hub.

Only the requested genes are retained. Primary solid tumor (sample type 01)
is kept, and replicate aliquots are averaged to the patient. Full matrices
are not written to the repository.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import urllib.request

import numpy as np
import pandas as pd

GDC_MIRRORS = [
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download",
    "https://gdc.xenahubs.net/download",
]
PROBEMAP_NAME = "gencode.v36.annotation.gtf.gene.probemap"


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: str, retries: int = 3) -> None:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tcga-trop2-cldn4-mediation/1.0"})
            with urllib.request.urlopen(req, timeout=180) as resp, open(tmp, "wb") as out:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            os.replace(tmp, dest)
            return
        except Exception as exc:  # noqa: BLE001 - retry network errors
            last = exc
            if os.path.exists(tmp):
                os.remove(tmp)
            print(f"[retry {attempt}] {url} :: {exc}", flush=True)
    raise RuntimeError(f"failed to download {url}: {last}")


def load_probemap(cache: str, symbols: list[str]) -> dict[str, str]:
    """Map gene symbols to the longest GENCODE v36 gene on a primary chromosome."""
    path = os.path.join(cache, PROBEMAP_NAME)
    errors = []
    if not (os.path.exists(path) and os.path.getsize(path) > 0):
        for base in GDC_MIRRORS:
            try:
                download(f"{base}/{PROBEMAP_NAME}", path)
                break
            except Exception as exc:  # noqa: BLE001
                errors.append(str(exc))
        if not (os.path.exists(path) and os.path.getsize(path) > 0):
            raise RuntimeError(f"probemap download failed: {errors}")
    table = pd.read_csv(path, sep="\t")
    mapping = {}
    for symbol in symbols:
        sub = table[table["gene"] == symbol].copy()
        if sub.empty:
            raise SystemExit(f"probemap missing {symbol}")
        sub = sub[sub["chrom"].astype(str).str.match(r"^chr([0-9]+|X|Y)$")]
        if sub.empty:
            raise SystemExit(f"no primary-chrom id for {symbol}")
        sub["span"] = sub["chromEnd"] - sub["chromStart"]
        # Longest span wins. Stable sort matches the earlier keratin script when spans tie.
        sub = sub.sort_values("span", ascending=False, kind="mergesort")
        mapping[symbol] = str(sub.iloc[0]["id"])
    if len(set(mapping.values())) != len(mapping):
        raise SystemExit(f"Ensembl id collision: {mapping}")
    return mapping


def patient_id(barcode: str) -> str | None:
    """Return the TCGA patient id for a primary solid tumor barcode, else None."""
    parts = barcode.replace(".", "-").split("-")
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3])


def _open_matrix(cohort: str):
    """Yield a readable binary stream of the STAR TPM gzip and the URL used."""
    errors = []
    for base in GDC_MIRRORS:
        url = f"{base}/TCGA-{cohort}.star_tpm.tsv.gz"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tcga-trop2-cldn4-mediation/1.0"})
            resp = urllib.request.urlopen(req, timeout=180)
            return resp, url
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{url} :: {exc}")
    raise RuntimeError(f"{cohort}: all mirrors failed: {errors}")


def extract_cohort(cohort: str, id_to_gene: dict[str, str], cache: str) -> dict:
    """Cache a patient-level table of the requested genes. Returns provenance."""
    os.makedirs(cache, exist_ok=True)
    cache_path = os.path.join(cache, f"{cohort}.primary01.tsv.gz")
    needed = set(id_to_gene.values())
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
        table = pd.read_csv(cache_path, sep="\t")
        if needed.issubset(table.columns) and "patient" in table.columns:
            return {
                "cohort": cohort,
                "url": "cache",
                "cache": cache_path,
                "sha256": sha256_file(cache_path),
                "n_patients": int(table.shape[0]),
                "n_01_columns": int(table["n_01_aliquots"].sum()) if "n_01_aliquots" in table.columns else None,
                "from_cache": True,
            }
        os.remove(cache_path)

    want = set(id_to_gene)
    resp, url = _open_matrix(cohort)
    print(f"[get] {cohort} {url}", flush=True)
    collected: dict[str, np.ndarray] = {}
    try:
        with gzip.GzipFile(fileobj=resp) as gz:
            header = gz.readline().decode("utf-8").rstrip("\n").split("\t")
            samples = header[1:]
            keep_idx = []
            patients = []
            for i, sample in enumerate(samples):
                pid = patient_id(sample)
                if pid is None:
                    continue
                keep_idx.append(i)
                patients.append(pid)
            n_01 = len(keep_idx)
            if n_01 < 15:
                raise RuntimeError(f"{cohort}: only {n_01} primary-tumor columns")
            keep_idx_arr = np.asarray(keep_idx, dtype=int)
            for raw in gz:
                gid, rest = raw.decode("utf-8").split("\t", 1)
                if gid not in want:
                    continue
                values = np.fromstring(rest, sep="\t", dtype=np.float64)
                if values.size != len(samples):
                    parts = rest.rstrip("\n").split("\t")
                    values = np.array([float(x) if x else np.nan for x in parts], dtype=np.float64)
                collected[id_to_gene[gid]] = values[keep_idx_arr]
    finally:
        resp.close()

    missing = sorted(set(id_to_gene.values()) - set(collected))
    if missing:
        raise RuntimeError(f"{cohort} missing genes: {missing}")
    frame = pd.DataFrame(collected)
    frame.insert(0, "patient", patients)
    collapsed = frame.groupby("patient", as_index=False).mean(numeric_only=True)
    counts = frame.groupby("patient").size().rename("n_01_aliquots")
    collapsed = collapsed.merge(counts, on="patient")
    tmp_path = cache_path + ".tmp"
    collapsed.to_csv(tmp_path, sep="\t", index=False, compression="gzip")
    os.replace(tmp_path, cache_path)
    return {
        "cohort": cohort,
        "url": url,
        "cache": cache_path,
        "sha256": sha256_file(cache_path),
        "n_patients": int(collapsed.shape[0]),
        "n_01_columns": int(n_01),
        "from_cache": False,
    }


def load_cached(cohort: str, cache: str) -> pd.DataFrame:
    path = os.path.join(cache, f"{cohort}.primary01.tsv.gz")
    frame = pd.read_csv(path, sep="\t")
    return frame.set_index("patient")
