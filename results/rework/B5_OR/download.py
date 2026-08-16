#!/usr/bin/env python3
"""Self-contained download + CLDN4 extraction for B5 open-ICI OR rework.

Sources
-------
Lung GEO (processed matrices only; no FASTQ):
  GSE126044, GSE135222, GSE166449, GSE207422
  GSE136961, GSE93157  (targeted panels; expected to lack CLDN4; still fetched)

Zenodo BHK Lab ICB TSVs (CC-BY-4.0), record 10.5281/zenodo.7058399:
  ICB_Mariathasan = IMvigor210 (Mariathasan Nature 2018)
  ICB_Gide, ICB_Liu, ICB_Riaz, ICB_Braun, ICB_Jung (= GSE135222 labels)

Liu raw RNA is dbGaP-controlled; Braun raw WES is EGA-controlled. Both
processed expression+response matrices are openly redistributed on Zenodo.
This script uses those public matrices and records that provenance.

Writes small patient-level tables under processed/. Raw archives stay in raw/.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import sys
import time
import zipfile
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

HERE = Path(__file__).resolve().parent
CLDN4_ENSEMBL = "ENSG00000189143"

ZENODO_FILES = {
    "ICB_Mariathasan.zip": "https://zenodo.org/api/records/7058399/files/ICB_Mariathasan.zip/content",
    "ICB_Gide.zip": "https://zenodo.org/api/records/7058399/files/ICB_Gide.zip/content",
    "ICB_Liu.zip": "https://zenodo.org/api/records/7058399/files/ICB_Liu.zip/content",
    "ICB_Riaz.zip": "https://zenodo.org/api/records/7058399/files/ICB_Riaz.zip/content",
    "ICB_Braun.zip": "https://zenodo.org/api/records/7058399/files/ICB_Braun.zip/content",
    "ICB_Jung.zip": "https://zenodo.org/api/records/7058399/files/ICB_Jung.zip/content",
}

GEO_FILES = {
    "GSE126044_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz",
    "GSE126044_counts.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
    "GSE135222_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz",
    "GSE135222_exp.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
    "GSE166449_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/matrix/GSE166449_series_matrix.txt.gz",
    "GSE166449_TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz",
    "GSE207422_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/matrix/GSE207422_series_matrix.txt.gz",
    "GSE207422_log2TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
    "GSE207422_metadata.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
    "GSE136961_TPM.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136961/suppl/GSE136961_TPM.tsv.gz",
    "GSE93157_raw.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93157/suppl/GSE93157_raw_data_values.txt.gz",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, retries: int = 4) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest.name} ({dest.stat().st_size} bytes)")
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err = None
    for i in range(retries):
        try:
            req = Request(url, headers={"User-Agent": "B5_OR-rework/1.0"})
            with urlopen(req, timeout=180) as r, open(tmp, "wb") as out:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            tmp.replace(dest)
            print(f"OK {dest.name} ({dest.stat().st_size} bytes)")
            return
        except Exception as e:
            last_err = e
            print(f"retry {i + 1} {dest.name}: {e}")
            time.sleep(2 ** (i + 2))
    raise RuntimeError(f"failed {url}: {last_err}")


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    opener = gzip.open if path.suffix == ".gz" or str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if not line.startswith("!Sample_"):
                continue
            key, *rest = line.rstrip("\n").split("\t")
            key = key[len("!Sample_") :]
            vals = [v.strip().strip('"') for v in rest]
            if key in fields:
                n = 0
                while f"{key}_{n}" in fields:
                    n += 1
                key = f"{key}_{n}"
            fields[key] = vals
    if not fields:
        raise ValueError(f"no Sample_ fields in {path}")
    n = max(len(v) for v in fields.values())
    for k, v in fields.items():
        if len(v) < n:
            v.extend([""] * (n - len(v)))
    return pd.DataFrame(fields)


def _is_cldn4(token: str) -> bool:
    t = token.strip().strip('"')
    u = t.upper()
    if u == "CLDN4" or u.startswith("CLDN4|") or u.startswith("CLDN4 "):
        return True
    if t.startswith(CLDN4_ENSEMBL):
        return True
    return False


def extract_gene_row_from_text(handle, sample_ids: list[str] | None = None) -> tuple[str, pd.Series] | None:
    header = handle.readline()
    if not header:
        return None
    parts = header.rstrip("\n").split("\t")
    # first cell may be a gene-id header or empty
    if sample_ids is None:
        if parts[0].strip().strip('"') in {"", "gene", "Gene", "gene_id", "Gene_id", "symbol", "SYMBOL"}:
            sample_ids = [p.strip().strip('"') for p in parts[1:]]
            gene_in_first = True
        else:
            # no gene header: all columns are samples
            sample_ids = [p.strip().strip('"') for p in parts]
            gene_in_first = False
            # this was the first data-less header; continue
    else:
        gene_in_first = True

    if gene_in_first:
        for line in handle:
            tok = line.split("\t", 1)[0]
            if _is_cldn4(tok):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
                gene = vals[0]
                data = pd.to_numeric(vals[1 : 1 + len(sample_ids)], errors="coerce")
                return gene, pd.Series(data, index=sample_ids, name="cldn4_raw")
        return None

    # header-was-samples case: first line already consumed as header; scan body
    for line in handle:
        vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
        if _is_cldn4(vals[0]):
            data = pd.to_numeric(vals[1 : 1 + len(sample_ids)], errors="coerce")
            return vals[0], pd.Series(data, index=sample_ids, name="cldn4_raw")
    return None


def extract_cldn4_from_path(path: Path) -> tuple[str, pd.Series] | None:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return extract_gene_row_from_text(f)
    with open(path, "rt") as f:
        return extract_gene_row_from_text(f)


def extract_cldn4_from_zip(zip_path: Path, inner_name: str) -> tuple[str, pd.Series] | None:
    with zipfile.ZipFile(zip_path) as z:
        with z.open(inner_name) as raw:
            handle = io.TextIOWrapper(raw, encoding="utf-8")
            return extract_gene_row_from_text(handle)


def extract_zip_member(zip_path: Path, member: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as z:
        with z.open(member) as src, open(dest, "wb") as out:
            out.write(src.read())
    return dest


def scan_gene_names(path: Path, limit: int = 200000) -> bool:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        f.readline()
        for i, line in enumerate(f):
            if i > limit:
                break
            if _is_cldn4(line.split("\t", 1)[0]):
                return True
    return False


def icb_table(meta: pd.DataFrame, expr: pd.Series, cohort: str, cancer: str, source: str, assay: str) -> pd.DataFrame:
    meta = meta.copy()
    meta["patientid"] = meta["patientid"].astype(str)
    expr.index = expr.index.astype(str)
    overlap = [p for p in meta["patientid"] if p in expr.index]
    out = meta.loc[meta["patientid"].isin(overlap)].copy()
    out["cldn4_raw"] = out["patientid"].map(expr)
    out["cohort"] = cohort
    out["cancer"] = cancer
    out["source"] = source
    out["assay"] = assay
    keep = [
        "cohort",
        "cancer",
        "source",
        "assay",
        "patientid",
        "cldn4_raw",
        "response",
        "recist",
        "treatment",
        "rna",
        "survival_time_pfs",
        "event_occurred_pfs",
        "survival_unit",
    ]
    for c in keep:
        if c not in out.columns:
            out[c] = pd.NA
    return out[keep]


def build_icb(raw: Path, processed: Path, provenance: list[dict]) -> None:
    specs = [
        (
            "ICB_Mariathasan.zip",
            "ICB_Mariathasan_metadata.tsv",
            "ICB_Mariathasan_expr.tsv",
            "IMvigor210_Mariathasan",
            "urothelial",
            "zenodo7058399",
            "TPM",
        ),
        (
            "ICB_Gide.zip",
            "ICB_Gide_metadata.tsv",
            "ICB_Gide_expr_gene_tpm.tsv",
            "Gide",
            "melanoma",
            "zenodo7058399",
            "TPM",
        ),
        (
            "ICB_Liu.zip",
            "ICB_Liu_metadata.tsv",
            "ICB_Liu_expr.tsv",
            "Liu",
            "melanoma",
            "zenodo7058399",
            "FPKM",
        ),
        (
            "ICB_Riaz.zip",
            "ICB_Riaz_metadata.tsv",
            "ICB_Riaz_expr_gene_tpm.tsv",
            "Riaz",
            "melanoma",
            "zenodo7058399",
            "TPM",
        ),
        (
            "ICB_Braun.zip",
            "ICB_Braun_metadata.tsv",
            "ICB_Braun_expr.tsv",
            "Braun",
            "ccRCC",
            "zenodo7058399",
            "TPM",
        ),
        (
            "ICB_Jung.zip",
            "ICB_Jung_metadata.tsv",
            "ICB_Jung_expr_gene_tpm.tsv",
            "GSE135222_Jung",
            "NSCLC",
            "zenodo7058399+GEO",
            "TPM",
        ),
    ]
    for zip_name, meta_name, expr_name, cohort, cancer, source, assay in specs:
        zpath = raw / zip_name
        meta_dest = raw / f"unz_{cohort}" / meta_name
        extract_zip_member(zpath, meta_name, meta_dest)
        meta = pd.read_csv(meta_dest, sep="\t")
        hit = extract_cldn4_from_zip(zpath, expr_name)
        if hit is None:
            print(f"WARNING: CLDN4 missing in {cohort} {expr_name}")
            continue
        gene, series = hit
        tab = icb_table(meta, series, cohort, cancer, source, assay)
        out = processed / f"{cohort}.csv"
        tab.to_csv(out, index=False)
        print(f"wrote {out.name} n={len(tab)} gene={gene} R/NR={tab['response'].value_counts(dropna=False).to_dict()}")
        provenance.append(
            {
                "cohort": cohort,
                "file": zip_name,
                "inner_expr": expr_name,
                "gene_token": gene,
                "n_rows": int(len(tab)),
                "sha256": sha256_file(zpath),
                "url": ZENODO_FILES[zip_name],
            }
        )


def build_gse126044(raw: Path, processed: Path, provenance: list[dict]) -> None:
    sm = parse_series_matrix(raw / "GSE126044_series_matrix.txt.gz")
    titles = sm["title"].astype(str)
    # titles like RNA-seq_Dis_01
    sample_ids = titles.str.replace(r"^RNA-seq_", "", regex=True)
    resp_col = None
    for c in sm.columns:
        if sm[c].astype(str).str.contains("responder", case=False, na=False).any():
            resp_col = c
            break
    if resp_col is None:
        raise RuntimeError("GSE126044: no responder characteristic")
    resp = (
        sm[resp_col]
        .astype(str)
        .str.replace(r"^patient response:\s*", "", regex=True, flags=re.I)
        .str.lower()
        .str.strip()
    )
    resp_map = resp.map({"responder": "R", "non-responder": "NR", "nonresponder": "NR"})
    hit = extract_cldn4_from_path(raw / "GSE126044_counts.txt.gz")
    if hit is None:
        raise RuntimeError("GSE126044: CLDN4 not found")
    gene, series = hit
    df = pd.DataFrame(
        {
            "cohort": "GSE126044",
            "cancer": "NSCLC",
            "source": "GEO",
            "assay": "counts",
            "patientid": sample_ids.values,
            "cldn4_raw": sample_ids.map(series).values,
            "response": resp_map.values,
            "recist": pd.NA,
            "treatment": "PD-1",
            "rna": "counts",
            "survival_time_pfs": pd.NA,
            "event_occurred_pfs": pd.NA,
            "survival_unit": pd.NA,
            "geo_title": titles.values,
        }
    )
    df = df.dropna(subset=["cldn4_raw"])
    df.to_csv(processed / "GSE126044.csv", index=False)
    print(f"wrote GSE126044 n={len(df)} gene={gene} {df['response'].value_counts(dropna=False).to_dict()}")
    provenance.append(
        {
            "cohort": "GSE126044",
            "file": "GSE126044_counts.txt.gz",
            "gene_token": gene,
            "n_rows": int(len(df)),
            "sha256": sha256_file(raw / "GSE126044_counts.txt.gz"),
            "url": GEO_FILES["GSE126044_counts.txt.gz"],
        }
    )


def build_gse166449(raw: Path, processed: Path, provenance: list[dict]) -> None:
    sm = parse_series_matrix(raw / "GSE166449_series_matrix.txt.gz")
    titles = sm["title"].astype(str)
    desc = sm["description"].astype(str)
    resp = titles.map(
        lambda t: "R" if re.search(r"responder", t, re.I) and not re.search(r"non", t, re.I) else ("NR" if re.search(r"non.?responder", t, re.I) else pd.NA)
    )
    hit = extract_cldn4_from_path(raw / "GSE166449_TPM.txt.gz")
    if hit is None:
        raise RuntimeError("GSE166449: CLDN4 not found")
    gene, series = hit
    df = pd.DataFrame(
        {
            "cohort": "GSE166449",
            "cancer": "NSCLC",
            "source": "GEO",
            "assay": "TPM",
            "patientid": desc.values,
            "cldn4_raw": desc.map(series).values,
            "response": resp.values,
            "recist": pd.NA,
            "treatment": "ICI",
            "rna": "tpm",
            "survival_time_pfs": pd.NA,
            "event_occurred_pfs": pd.NA,
            "survival_unit": pd.NA,
            "geo_title": titles.values,
        }
    )
    df = df.dropna(subset=["cldn4_raw"])
    df.to_csv(processed / "GSE166449.csv", index=False)
    print(f"wrote GSE166449 n={len(df)} gene={gene} {df['response'].value_counts(dropna=False).to_dict()}")
    provenance.append(
        {
            "cohort": "GSE166449",
            "file": "GSE166449_TPM.txt.gz",
            "gene_token": gene,
            "n_rows": int(len(df)),
            "sha256": sha256_file(raw / "GSE166449_TPM.txt.gz"),
            "url": GEO_FILES["GSE166449_TPM.txt.gz"],
        }
    )


def build_gse135222_dcb(raw: Path, processed: Path, provenance: list[dict]) -> None:
    """GEO-only DCB table (PFS >= 180 days). Same patients as ICB_Jung."""
    sm = parse_series_matrix(raw / "GSE135222_series_matrix.txt.gz")
    titles = sm["title"].astype(str)
    # titles like "NSCLC 990"
    sid = titles.str.replace(r"\s+", "", regex=True)
    pfs_event = None
    pfs_time = None
    for c in sm.columns:
        s = sm[c].astype(str)
        if s.str.contains("progression-free survival", case=False, na=False).any() and pfs_event is None:
            pfs_event = pd.to_numeric(s.str.replace(r".*:\s*", "", regex=True), errors="coerce")
        if s.str.contains("pfs.time", case=False, na=False).any():
            pfs_time = pd.to_numeric(s.str.replace(r".*:\s*", "", regex=True), errors="coerce")
    hit = extract_cldn4_from_path(raw / "GSE135222_exp.tsv.gz")
    if hit is None:
        raise RuntimeError("GSE135222: CLDN4 not found")
    gene, series = hit
    # GEO columns NSCLC990 vs titles NSCLC990
    cldn = sid.map(series)
    if cldn.isna().all():
        # try NSCLC_990 / P_990
        cldn = sid.str.replace("NSCLC", "NSCLC", regex=False).map(series)
    dcb = pd.Series(pd.NA, index=sid.index, dtype=object)
    if pfs_time is not None:
        dcb = pd.Series(["R" if t >= 180 else "NR" for t in pfs_time], index=sid.index)
    df = pd.DataFrame(
        {
            "cohort": "GSE135222_DCB180",
            "cancer": "NSCLC",
            "source": "GEO",
            "assay": "TPM",
            "patientid": sid.values,
            "cldn4_raw": cldn.values,
            "response": dcb.values,
            "recist": pd.NA,
            "treatment": "PD-1/PD-L1",
            "rna": "tpm",
            "survival_time_pfs": pfs_time.values if pfs_time is not None else pd.NA,
            "event_occurred_pfs": pfs_event.values if pfs_event is not None else pd.NA,
            "survival_unit": "day",
            "geo_title": titles.values,
        }
    )
    df = df.dropna(subset=["cldn4_raw"])
    df.to_csv(processed / "GSE135222_DCB180.csv", index=False)
    print(f"wrote GSE135222_DCB180 n={len(df)} gene={gene} {df['response'].value_counts(dropna=False).to_dict()}")
    provenance.append(
        {
            "cohort": "GSE135222_DCB180",
            "file": "GSE135222_exp.tsv.gz",
            "gene_token": gene,
            "n_rows": int(len(df)),
            "note": "DCB = PFS>=180 days; same patients as GSE135222_Jung; NOT an independent study",
            "sha256": sha256_file(raw / "GSE135222_exp.tsv.gz"),
            "url": GEO_FILES["GSE135222_exp.tsv.gz"],
        }
    )


def build_gse207422(raw: Path, processed: Path, provenance: list[dict]) -> None:
    meta = pd.read_excel(raw / "GSE207422_metadata.xlsx")
    meta = meta.dropna(subset=["Sample", "Patient"]).copy()
    hit = extract_cldn4_from_path(raw / "GSE207422_log2TPM.txt.gz")
    if hit is None:
        raise RuntimeError("GSE207422: CLDN4 not found")
    gene, series = hit
    meta["cldn4_raw"] = meta["Sample"].astype(str).map(series)
    pr = meta["Pathologic Response"].astype(str)
    resp = pr.map(
        lambda x: "R" if x.upper().startswith("MPR") else ("NR" if x.upper().startswith("NMPR") else pd.NA)
    )
    recist = meta["RECIST"].astype(str).str.upper().replace({"NAN": pd.NA})
    df = pd.DataFrame(
        {
            "cohort": "GSE207422",
            "cancer": "NSCLC",
            "source": "GEO",
            "assay": "log2TPM",
            "patientid": meta["Patient"].astype(str).values,
            "cldn4_raw": meta["cldn4_raw"].values,
            "response": resp.values,
            "recist": recist.values,
            "treatment": "neoadjuvant_ICI_chemo",
            "rna": "log2tpm",
            "survival_time_pfs": pd.NA,
            "event_occurred_pfs": pd.NA,
            "survival_unit": pd.NA,
            "sample_id": meta["Sample"].astype(str).values,
            "pathologic_response": meta["Pathologic Response"].astype(str).values,
        }
    )
    df = df.dropna(subset=["cldn4_raw"])
    df.to_csv(processed / "GSE207422.csv", index=False)
    print(f"wrote GSE207422 n={len(df)} gene={gene} {df['response'].value_counts(dropna=False).to_dict()}")
    provenance.append(
        {
            "cohort": "GSE207422",
            "file": "GSE207422_log2TPM.txt.gz",
            "gene_token": gene,
            "n_rows": int(len(df)),
            "sha256": sha256_file(raw / "GSE207422_log2TPM.txt.gz"),
            "url": GEO_FILES["GSE207422_log2TPM.txt.gz"],
        }
    )


def document_missing_gene(raw: Path, processed: Path, provenance: list[dict]) -> None:
    rows = []
    for acc, fname in [
        ("GSE136961", "GSE136961_TPM.tsv.gz"),
        ("GSE93157", "GSE93157_raw.txt.gz"),
    ]:
        path = raw / fname
        present = scan_gene_names(path)
        rows.append(
            {
                "cohort": acc,
                "cancer": "NSCLC",
                "file": fname,
                "cldn4_present": present,
                "reason": "targeted immune panel; CLDN4 not measured" if not present else "CLDN4 present",
                "usable_for_OR": bool(present),
            }
        )
        provenance.append(
            {
                "cohort": acc,
                "file": fname,
                "cldn4_present": present,
                "sha256": sha256_file(path),
                "url": GEO_FILES[fname],
            }
        )
        print(f"{acc} CLDN4 present={present}")
    pd.DataFrame(rows).to_csv(processed / "unusable_no_CLDN4.csv", index=False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default=str(HERE / "raw"))
    ap.add_argument("--out-dir", default=str(HERE / "processed"))
    args = ap.parse_args()
    raw = Path(args.raw_dir)
    processed = Path(args.out_dir)
    raw.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)

    print("== downloading ==")
    for name, url in {**GEO_FILES, **ZENODO_FILES}.items():
        download(url, raw / name)

    provenance: list[dict] = []
    print("== extracting ==")
    build_icb(raw, processed, provenance)
    build_gse126044(raw, processed, provenance)
    build_gse166449(raw, processed, provenance)
    build_gse135222_dcb(raw, processed, provenance)
    build_gse207422(raw, processed, provenance)
    document_missing_gene(raw, processed, provenance)

    prov_path = processed / "provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"wrote {prov_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
