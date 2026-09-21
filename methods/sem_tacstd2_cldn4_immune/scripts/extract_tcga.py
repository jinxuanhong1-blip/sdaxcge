#!/usr/bin/env python3
"""Stream TCGA GDC STAR log2(TPM+1) and keep the genes this SEM uses.

Source is the same UCSC Xena GDC hub used by the keratin coexpression
analysis. Primary solid tumor (sample type 01) only. Replicate aliquots
of one patient are averaged. The full matrix is not written to the repo.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CACHE = os.environ.get("TCGA_SEM_CACHE", "/tmp/tcga_sem")
GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PROBEMAP_URL = f"{GDC}/gencode.v36.annotation.gtf.gene.probemap"

# Handoff keratin funnel, plus LUSC as a lung sensitivity cohort.
COHORTS = ["LUAD", "LUSC", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
# Patient counts from the prior Xena extract (results/tcga_trop2_cldn_keratin).
EXPECTED_N = {
    "LUAD": 516,
    "LUSC": 501,
    "BRCA": 1095,
    "CESC": 304,
    "KIRC": 533,
    "STAD": 412,
    "BLCA": 406,
    "PAAD": 178,
}
GENES = [
    "TACSTD2",
    "CLDN4",
    "CD8A",
    "CD8B",
    "CD3D",
    "CD3E",
    "CD3G",
    "NKG7",
    "GZMA",
    "GZMB",
    "PRF1",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT5",
    "KRT6A",
    "KRT6B",
    "KRT14",
]


def download(url: str, dest: str) -> None:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    req = urllib.request.Request(url, headers={"User-Agent": "sem-tacstd2-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp, open(tmp, "wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    os.replace(tmp, dest)


def load_probemap() -> dict[str, str]:
    path = os.path.join(CACHE, "gencode.v36.annotation.gtf.gene.probemap")
    download(PROBEMAP_URL, path)
    pm = pd.read_csv(path, sep="\t")
    mapping = {}
    for symbol in GENES:
        sub = pm[pm["gene"] == symbol].copy()
        if sub.empty:
            raise SystemExit(f"probemap missing {symbol}")
        sub = sub[sub["chrom"].astype(str).str.match(r"^chr([0-9]+|X|Y)$")]
        sub["span"] = sub["chromEnd"] - sub["chromStart"]
        sub = sub.sort_values("span", ascending=False)
        mapping[symbol] = str(sub.iloc[0]["id"])
    if len(set(mapping.values())) != len(mapping):
        raise SystemExit(f"Ensembl id collision: {mapping}")
    return mapping


def patient_id(barcode: str) -> str | None:
    parts = barcode.replace(".", "-").split("-")
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3])


def extract_cohort(cohort: str, id_to_gene: dict[str, str]) -> dict:
    cache_path = os.path.join(CACHE, f"{cohort}.primary01.tsv.gz")
    url = f"{GDC}/TCGA-{cohort}.star_tpm.tsv.gz"
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
        table = pd.read_csv(cache_path, sep="\t")
        n = int(table["patient"].nunique())
        if n != EXPECTED_N[cohort]:
            raise SystemExit(f"{cohort} cache n={n}, expected {EXPECTED_N[cohort]}")
        print(f"[cache] {cohort} n={n}", flush=True)
        return {"cohort": cohort, "url": url, "n_patients": n, "from_cache": True}

    want = {gid.encode("ascii"): symbol for gid, symbol in id_to_gene.items()}
    print(f"[get] {cohort} {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "sem-tacstd2-cldn4/1.0"})
    collected: dict[str, np.ndarray] = {}
    with urllib.request.urlopen(req, timeout=300) as resp:
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
                gid, sep, rest = raw.partition(b"\t")
                if not sep or gid not in want:
                    continue
                values = np.fromstring(rest, sep="\t", dtype=float)
                if values.size != len(samples):
                    parts = rest.decode("utf-8").rstrip("\n").split("\t")
                    values = np.array([float(x) if x else np.nan for x in parts], dtype=float)
                collected[want[gid]] = values[keep_idx_arr]
    missing = [sym for sym in set(id_to_gene.values()) if sym not in collected]
    if missing:
        raise RuntimeError(f"{cohort} missing genes: {missing}")
    frame = pd.DataFrame(collected)
    frame.insert(0, "patient", patients)
    counts = frame.groupby("patient").size().rename("n_01_aliquots")
    collapsed = frame.groupby("patient", as_index=False).mean(numeric_only=True)
    collapsed = collapsed.merge(counts, on="patient")
    if len(collapsed) != EXPECTED_N[cohort]:
        raise SystemExit(f"{cohort} n={len(collapsed)}, expected {EXPECTED_N[cohort]}")
    os.makedirs(CACHE, exist_ok=True)
    collapsed.to_csv(cache_path, sep="\t", index=False, compression="gzip")
    print(f"[ok] {cohort} patients={len(collapsed)} aliquots={n_01}", flush=True)
    return {
        "cohort": cohort,
        "url": url,
        "n_patients": int(len(collapsed)),
        "n_01_columns": int(n_01),
        "from_cache": False,
    }


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    os.makedirs(CACHE, exist_ok=True)
    symbol_to_id = load_probemap()
    id_to_gene = {ens: sym for sym, ens in symbol_to_id.items()}
    provenance = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(extract_cohort, cohort, id_to_gene): cohort for cohort in COHORTS}
        for fut in as_completed(futures):
            provenance.append(fut.result())
    frames = []
    for cohort in COHORTS:
        frame = pd.read_csv(os.path.join(CACHE, f"{cohort}.primary01.tsv.gz"), sep="\t")
        frame.insert(0, "cohort", cohort)
        frames.append(frame)
    out = pd.concat(frames, ignore_index=True)
    dest = os.path.join(ROOT, "data", "tcga_primary01_genes.tsv.gz")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    out.to_csv(dest, sep="\t", index=False, compression="gzip")
    meta = {
        "source": "UCSC Xena GDC hub STAR log2(TPM+1), GENCODE v36 Ensembl IDs",
        "probemap": PROBEMAP_URL,
        "ensembl_ids": symbol_to_id,
        "sample_type": "primary solid tumor, barcode field 4 startswith 01; replicate aliquots averaged",
        "cohorts": provenance,
        "n_rows": int(len(out)),
        "sha256": sha256_file(dest),
        "output": dest,
    }
    with open(os.path.join(ROOT, "data", "tcga_provenance.json"), "w") as handle:
        json.dump(meta, handle, indent=2)
    print(json.dumps({"n_rows": meta["n_rows"], "sha256": meta["sha256"]}, indent=2))


if __name__ == "__main__":
    main()
