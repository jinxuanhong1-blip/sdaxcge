#!/usr/bin/env python3
"""Extract panel genes and patient-level summaries for the winning pair.

GSE131907 + GSE205335 only. Patient is the unit (tumor samples pooled).
Sender split is CLDN4-only (global malignant median within each dataset).
"""
from __future__ import annotations

import gzip
import json
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = Path("/tmp/winpair_multinichenet_cldn4")
sys.path.insert(0, str(ROOT / "scripts"))
from gene_sets import (  # noqa: E402
    CYTOTOXICITY,
    EXHAUSTION,
    EXTRA_TNK,
    HU_MALIGNANT_SUB,
    HU_TNK_TOTAL,
    IFN,
    KIM_MALIGNANT_SUBTYPES,
    KIM_TNK_TYPES,
    KIM_TUMOR_ORIGINS,
    LINEAGE,
)

LR = DATA / "lr_network.tsv"
OUT_STATS = CACHE / "patient_gene_stats.parquet"
OUT_PAT = CACHE / "patient_inventory.tsv"
OUT_META = CACHE / "extract_meta.json"

MIN_HIGH = 10
MIN_LOW = 10
MIN_TNK = 20


def log(msg: str) -> None:
    print(msg, flush=True)


def panel_genes() -> set[str]:
    lr = pd.read_csv(LR, sep="\t")
    genes = set(lr["from"].astype(str)) | set(lr["to"].astype(str))
    for block in (CYTOTOXICITY, IFN, EXHAUSTION, EXTRA_TNK, LINEAGE):
        genes.update(block)
    return genes


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    total = np.maximum(total.astype(np.float64), 1.0)
    return np.log1p(1e4 * umi.astype(np.float64) / total).astype(np.float32)


def parse_series_matrix(path: Path) -> pd.DataFrame:
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


def parse_geo_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with gzip.open(path, "rt", errors="replace") as handle:
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
    if "platform" in metadata.columns:
        read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
        metadata["orig.ident_soft"] = (
            metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
        )
    return metadata


def stream_umi(path: Path, keep: set[str]) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, int]:
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


def load_rds_selected(path: Path, wanted: set[str]):
    from scipy import sparse
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.name.endswith(".gz"):
            matrix_path = Path(tmp) / path.name.replace(".gz", "")
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
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel().astype(np.float32)
    log(f"extracted {len(extracted)} / {len(wanted)} genes")
    return extracted, library_umi, barcodes, genes.tolist()


def summarize_patients(
    dataset: str,
    patient: np.ndarray,
    is_malig: np.ndarray,
    is_tnk: np.ndarray,
    cld_hi: np.ndarray,
    cld_lo: np.ndarray,
    totals: np.ndarray,
    store: dict[str, np.ndarray],
    extra_meta: dict[str, np.ndarray] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = sorted(store)
    rows_pat = []
    rows_gene = []
    for pid in pd.unique(patient):
        if pid is None or (isinstance(pid, float) and np.isnan(pid)):
            continue
        m = patient == pid
        hi = m & cld_hi
        lo = m & cld_lo
        tnk = m & is_tnk
        mal = m & is_malig
        rec = {
            "dataset": dataset,
            "patient": str(pid),
            "n_cells": int(m.sum()),
            "n_malignant": int(mal.sum()),
            "n_cldn4_high": int(hi.sum()),
            "n_cldn4_low": int(lo.sum()),
            "n_tnk": int(tnk.sum()),
        }
        if extra_meta:
            for k, arr in extra_meta.items():
                vals = pd.unique(arr[m])
                rec[k] = str(vals[0]) if len(vals) == 1 else "|".join(sorted(map(str, vals))[:6])
        if "CLDN4" in store and mal.any():
            rec["mal_CLDN4_mean"] = float(log1p_cp10k(store["CLDN4"], totals)[mal].mean())
            rec["mal_CLDN4_pct_pos"] = float((store["CLDN4"][mal] > 0).mean())
        rec["tnk_frac"] = float(tnk.sum() / m.sum()) if m.any() else np.nan
        rec["eligible_paired"] = (
            int(hi.sum()) >= MIN_HIGH and int(lo.sum()) >= MIN_LOW and int(tnk.sum()) >= MIN_TNK
        )
        rows_pat.append(rec)
        for compartment, mask in (
            ("mal_high", hi),
            ("mal_low", lo),
            ("tnk", tnk),
            ("malignant", mal),
        ):
            n = int(mask.sum())
            if n == 0:
                continue
            tot = totals[mask]
            for g in genes:
                umi = store[g][mask]
                lx = log1p_cp10k(umi, tot)
                rows_gene.append(
                    {
                        "dataset": dataset,
                        "patient": str(pid),
                        "compartment": compartment,
                        "gene": g,
                        "n_cells": n,
                        "mean_log1p": float(lx.mean()),
                        "frac_pos": float((umi > 0).mean()),
                    }
                )
    return pd.DataFrame(rows_pat), pd.DataFrame(rows_gene)


def extract_gse131907(keep: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    ann = pd.read_csv(
        CACHE / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str
    )
    meta = parse_series_matrix(CACHE / "GSE131907_series_matrix.txt.gz")
    sample_meta = meta.rename(columns={"title": "Sample"})
    keep_cols = [
        c
        for c in ["Sample", "geo_accession", "patient_id", "tumor_stage", "tissue_origin_abbrevation"]
        if c in sample_meta.columns
    ]
    sample_meta = sample_meta[keep_cols].drop_duplicates("Sample")
    cells, store, totals, n_genes = stream_umi(
        CACHE / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", keep
    )
    per = ann.set_index("Index").reindex(cells)
    if per["Sample"].isna().any():
        raise SystemExit("GSE131907 matrix cell IDs do not align with annotation Index")
    sample = per["Sample"].to_numpy()
    origin = per["Sample_Origin"].to_numpy()
    ctype = per["Cell_type"].fillna("").to_numpy()
    csub = per["Cell_subtype"].fillna("").to_numpy()
    patient_map = (
        sample_meta.set_index("Sample")["patient_id"].to_dict()
        if "patient_id" in sample_meta.columns
        else {}
    )
    patient = np.array([patient_map.get(s, s) for s in sample], dtype=object)
    is_tumor = np.isin(origin, list(KIM_TUMOR_ORIGINS))
    is_malig = is_tumor & np.isin(csub, list(KIM_MALIGNANT_SUBTYPES))
    is_tnk = is_tumor & np.isin(ctype, list(KIM_TNK_TYPES))
    if "CLDN4" not in store:
        raise SystemExit("CLDN4 missing from GSE131907 UMI")
    cldn4 = log1p_cp10k(store["CLDN4"], totals)
    thr = float(np.median(cldn4[is_malig])) if is_malig.any() else float("nan")
    cld_hi = is_malig & (cldn4 >= thr)
    cld_lo = is_malig & (cldn4 < thr)
    log(
        f"[GSE131907] cells={len(cells)} mal={int(is_malig.sum())} "
        f"high={int(cld_hi.sum())} low={int(cld_lo.sum())} tnk={int(is_tnk.sum())} thr={thr:.4f}"
    )
    extra = {
        "sample_origin": origin.astype(str),
    }
    pat, gene = summarize_patients(
        "GSE131907", patient, is_malig, is_tnk, cld_hi, cld_lo, totals, store, extra
    )
    info = {
        "n_cells": int(len(cells)),
        "n_genes_in_file": int(n_genes),
        "n_genes_kept": int(len(store)),
        "n_malignant": int(is_malig.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_malig_tlung_ts": int(((origin == "tLung") & np.isin(csub, ["tS1", "tS2", "tS3"])).sum()),
        "n_malig_author_string": int((csub == "Malignant cells").sum()),
        "cldn4_threshold": thr,
        "n_geo_samples": int(sample_meta["Sample"].nunique()) if "Sample" in sample_meta else None,
        "n_geo_patients": int(sample_meta["patient_id"].nunique()) if "patient_id" in sample_meta else None,
    }
    return pat, gene, info


def extract_gse205335(keep: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    ident = pd.read_csv(CACHE / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    log(f"[GSE205335] identity columns={list(ident.columns)}")
    soft = parse_geo_soft(CACHE / "GSE205335_family.soft.gz")
    extracted, library_umi, barcodes, all_genes = load_rds_selected(
        CACHE / "GSE205335_Lung_IO_UMI_matrix.rds.gz", keep
    )
    if "barcode" not in ident.columns:
        # first column is often the barcode
        ident = ident.rename(columns={ident.columns[0]: "barcode"})
    indexed = ident.drop_duplicates("barcode").set_index("barcode")
    cells = indexed.reindex(barcodes).reset_index()
    cells = cells.rename(columns={"index": "barcode"})
    n_miss = int(cells.iloc[:, 1].isna().sum()) if cells.shape[1] > 1 else 0
    log(f"[GSE205335] identity reindex missing={n_miss} / {len(barcodes)}")

    # patient
    if "patient" in cells.columns:
        patient = cells["patient"].astype(str).to_numpy()
    elif "Patient" in cells.columns:
        patient = cells["Patient"].astype(str).to_numpy()
    else:
        # join SOFT
        key = "orig.ident" if "orig.ident" in cells.columns else None
        if key and "orig.ident_soft" in soft.columns:
            cells = cells.merge(
                soft.rename(columns={"orig.ident_soft": "orig.ident", "patient": "patient_soft"}),
                on="orig.ident",
                how="left",
            )
        if "patient" in cells.columns:
            patient = cells["patient"].astype(str).to_numpy()
        elif "patient_soft" in cells.columns:
            patient = cells["patient_soft"].astype(str).to_numpy()
        else:
            raise SystemExit(f"cannot find patient in GSE205335; cols={list(cells.columns)}")

    tissue = None
    for cand in ("tissue", "Tissue", "sample_type"):
        if cand in cells.columns:
            tissue = cells[cand].astype(str).to_numpy()
            break
    if tissue is None and "tissue" in soft.columns and "orig.ident" in cells.columns:
        tmap = soft.set_index("orig.ident_soft")["tissue"].to_dict() if "orig.ident_soft" in soft.columns else {}
        if "orig.ident" in cells.columns:
            tissue = np.array([str(tmap.get(x, "")) for x in cells["orig.ident"]], dtype=object)

    lineage_sub = cells["lineage.sub"].astype(str).to_numpy() if "lineage.sub" in cells.columns else None
    lineage_total = cells["lineage.total"].astype(str).to_numpy() if "lineage.total" in cells.columns else None
    if lineage_sub is None or lineage_total is None:
        raise SystemExit(f"missing lineage columns: {list(cells.columns)}")

    is_tumor = np.ones(len(barcodes), dtype=bool)
    if tissue is not None:
        is_tumor = ~pd.Series(tissue).str.startswith("Normal", na=False).to_numpy()
    is_malig = is_tumor & (lineage_sub == HU_MALIGNANT_SUB)
    is_tnk = is_tumor & (lineage_total == HU_TNK_TOTAL)
    if "CLDN4" not in extracted:
        raise SystemExit("CLDN4 missing from GSE205335 UMI")
    cldn4 = log1p_cp10k(extracted["CLDN4"], library_umi)
    thr = float(np.median(cldn4[is_malig])) if is_malig.any() else float("nan")
    cld_hi = is_malig & (cldn4 >= thr)
    cld_lo = is_malig & (cldn4 < thr)
    log(
        f"[GSE205335] cells={len(barcodes)} mal={int(is_malig.sum())} "
        f"high={int(cld_hi.sum())} low={int(cld_lo.sum())} tnk={int(is_tnk.sum())} thr={thr:.4f}"
    )
    extra = {}
    if tissue is not None:
        extra["tissue"] = tissue.astype(str)
    if "cancer_subtype" in cells.columns:
        extra["cancer_subtype"] = cells["cancer_subtype"].astype(str).to_numpy()
    elif "cancer subtype" in cells.columns:
        extra["cancer_subtype"] = cells["cancer subtype"].astype(str).to_numpy()
    if "recist" in cells.columns:
        extra["recist"] = cells["recist"].astype(str).to_numpy()
    pat, gene = summarize_patients(
        "GSE205335", patient, is_malig, is_tnk, cld_hi, cld_lo, library_umi, extracted, extra
    )
    info = {
        "n_cells": int(len(barcodes)),
        "n_genes_in_file": int(len(all_genes)),
        "n_genes_kept": int(len(extracted)),
        "n_malignant": int(is_malig.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "n_tnk": int(is_tnk.sum()),
        "cldn4_threshold": thr,
        "identity_columns": list(ident.columns),
    }
    return pat, gene, info


def main() -> int:
    keep = panel_genes()
    log(f"panel requested: {len(keep)}")
    p1, g1, i1 = extract_gse131907(keep)
    p2, g2, i2 = extract_gse205335(keep)
    pat = pd.concat([p1, p2], ignore_index=True)
    gene = pd.concat([g1, g2], ignore_index=True)
    pat.to_csv(OUT_PAT, sep="\t", index=False)
    gene.to_parquet(OUT_STATS)
    meta = {
        "GSE131907": i1,
        "GSE205335": i2,
        "n_patient_rows": int(len(pat)),
        "n_gene_stat_rows": int(len(gene)),
        "n_eligible_paired": int(pat["eligible_paired"].sum()),
        "skipped": ["GSE207422"],
    }
    OUT_META.write_text(json.dumps(meta, indent=2) + "\n")
    log(f"wrote {OUT_PAT} n={len(pat)} eligible={int(pat['eligible_paired'].sum())}")
    log(f"wrote {OUT_STATS} rows={len(gene)}")
    print(pat.groupby("dataset")["eligible_paired"].agg(["size", "sum"]).to_string(), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
