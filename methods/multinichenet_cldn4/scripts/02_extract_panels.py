#!/usr/bin/env python3
"""Extract ligand / receptor / T/NK-program panels for the three cohorts.

Keeps only author-malignant + T/NK cells. Runtime cache is /tmp.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (  # noqa: E402
    CYTOTOXICITY,
    EXHAUSTION,
    EXTRA_TNK,
    GSE131907_MALIGNANT_SUBTYPES,
    GSE131907_TNK_TYPES,
    GSE131907_TUMOR_ORIGINS,
    GSE207422_MALIGNANT,
    GSE207422_TNK,
    IFN,
    LINEAGE,
    MPR_LABELS,
    NMPR_LABELS,
)
from lib_io import (  # noqa: E402
    CACHE,
    DATA,
    load_gse205335_rds,
    log,
    parse_gse131907_series_matrix,
    parse_gse205335_soft,
    stream_dense_umi_gz,
    write_cohort_panel,
)

LR = DATA / "lr_network.tsv"


def panel_genes() -> set[str]:
    lr = pd.read_csv(LR, sep="\t")
    genes = set(lr["from"].astype(str)) | set(lr["to"].astype(str))
    for block in (CYTOTOXICITY, IFN, EXHAUSTION, EXTRA_TNK, LINEAGE):
        genes.update(block)
    return genes


def mpr_group(label: str) -> str | None:
    s = str(label).strip()
    if s in MPR_LABELS or s.upper() == "PCR":
        return "MPR"
    if s in NMPR_LABELS:
        return "NMPR"
    return None


def extract_gse131907(keep: set[str]) -> None:
    dest = CACHE / "GSE131907_panel.parquet"
    if dest.exists() and dest.stat().st_size > 1_000_000:
        log(f"HAVE {dest}")
        return
    ann = pd.read_csv(
        CACHE / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str
    )
    meta = parse_gse131907_series_matrix(CACHE / "GSE131907_series_matrix.txt.gz")
    sample_meta = meta.rename(
        columns={"title": "Sample", "tissue_origin_abbrevation": "Sample_Origin_geo"}
    )
    keep_cols = [
        c
        for c in ["Sample", "geo_accession", "patient_id", "tumor_stage", "Sample_Origin_geo"]
        if c in sample_meta.columns
    ]
    sample_meta = sample_meta[keep_cols].drop_duplicates("Sample")
    sample_meta.to_csv(CACHE / "GSE131907_sample_metadata.tsv", sep="\t", index=False)

    cells, store, totals, _ = stream_dense_umi_gz(
        CACHE / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", keep
    )
    per = ann.set_index("Index").reindex(cells)
    if per["Sample"].isna().any():
        raise SystemExit("GSE131907 matrix cell IDs do not align with annotation Index")
    patient_map = sample_meta.set_index("Sample")["patient_id"].to_dict()
    origin = per["Sample_Origin"].to_numpy()
    ctype = per["Cell_type"].fillna("").to_numpy()
    csub = per["Cell_subtype"].fillna("").to_numpy()
    sample = per["Sample"].to_numpy()
    is_tumor = np.isin(origin, list(GSE131907_TUMOR_ORIGINS))
    is_malig = is_tumor & np.isin(csub, list(GSE131907_MALIGNANT_SUBTYPES))
    is_tnk = is_tumor & np.isin(ctype, list(GSE131907_TNK_TYPES))
    meta_df = pd.DataFrame(
        {
            "patient": [patient_map.get(s, s) for s in sample],
            "sample": sample,
            "cohort": "GSE131907",
            "is_tumor": is_tumor,
            "is_malignant": is_malig,
            "is_tnk": is_tnk,
            "celltype": csub,
            "note": origin,
        },
        index=pd.Index(cells, name="barcode"),
    )
    genes = sorted(store)
    umi = pd.DataFrame({g: store[g] for g in genes}, index=meta_df.index)
    umi["ncount"] = totals
    out = meta_df.join(umi, how="inner")
    write_cohort_panel(out, dest)


def extract_gse205335(keep: set[str]) -> None:
    dest = CACHE / "GSE205335_panel.parquet"
    if dest.exists() and dest.stat().st_size > 1_000_000:
        log(f"HAVE {dest}")
        return
    identities = pd.read_csv(CACHE / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    if identities["barcode"].duplicated().any():
        raise SystemExit("GSE205335 Cell-identity barcodes are not unique")
    metadata = parse_gse205335_soft(CACHE / "GSE205335_family.soft.gz")
    metadata.to_csv(CACHE / "GSE205335_gsm_sample_metadata.csv", index=False)

    barcodes, store, totals = load_gse205335_rds(
        CACHE / "GSE205335_Lung_IO_UMI_matrix.rds.gz", keep
    )
    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise SystemExit(
            f"GSE205335 matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra"
        )
    cells = indexed.loc[barcodes].reset_index()
    cells = cells.merge(
        metadata[
            [
                c
                for c in [
                    "orig.ident",
                    "gsm",
                    "patient",
                    "tissue",
                    "recist",
                    "cancer_subtype",
                    "tumor_stage",
                    "platform",
                ]
                if c in metadata.columns
            ]
        ],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise SystemExit("GSE205335 identity samples did not match GEO metadata")
    is_normal = cells["tissue"].astype(str).str.startswith("Normal")
    is_tumor = ~is_normal
    is_malig = is_tumor & cells["lineage.sub"].eq("Malignant cells")
    is_tnk = is_tumor & cells["lineage.total"].eq("T/NK cells")
    meta_df = pd.DataFrame(
        {
            "patient": cells["patient"].to_numpy(),
            "sample": cells["orig.ident"].to_numpy(),
            "cohort": "GSE205335",
            "is_tumor": is_tumor.to_numpy(),
            "is_malignant": is_malig.to_numpy(),
            "is_tnk": is_tnk.to_numpy(),
            "celltype": cells["lineage.sub"].to_numpy(),
            "note": cells["cancer_subtype"].astype(str)
            + "/"
            + cells["recist"].astype(str),
        },
        index=pd.Index(barcodes, name="barcode"),
    )
    genes = sorted(store)
    umi = pd.DataFrame({g: store[g] for g in genes}, index=meta_df.index)
    umi["ncount"] = totals
    out = meta_df.join(umi, how="inner")
    write_cohort_panel(out, dest)


def extract_gse207422(keep: set[str]) -> None:
    dest = CACHE / "GSE207422_panel.parquet"
    if dest.exists() and dest.stat().st_size > 1_000_000:
        log(f"HAVE {dest}")
        return
    ann = pd.read_csv(DATA / "drmref_cell_annotation.tsv.gz", sep="\t")
    meta = pd.read_csv(DATA / "geo_scRNAseq_sample_metadata.tsv", sep="\t")
    meta = meta.rename(columns={"Sample": "orig.ident", "Patient": "patient_meta"})
    meta["mpr_group"] = meta["Pathologic_Response"].map(mpr_group)
    meta["timing"] = np.where(
        meta["Resource"].astype(str).str.contains("Post", case=False), "post", "pre"
    )
    ann = ann.merge(
        meta[["orig.ident", "mpr_group", "timing", "Pathology", "Pathologic_Response"]],
        on="orig.ident",
        how="left",
    )
    ann["mpr_group"] = ann["mpr_group"].fillna(ann["response"].map(mpr_group))
    ann = ann[ann["timing"].eq("post")].set_index("cell_barcode")

    cells, store, totals, _ = stream_dense_umi_gz(
        CACHE / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz", keep
    )
    common = ann.index.intersection(pd.Index(cells))
    log(f"GSE207422 post-tx barcode match {len(common)} / matrix {len(cells)} / annot {len(ann)}")
    pos = pd.Index(cells).get_indexer(common)
    store = {g: arr[pos] for g, arr in store.items()}
    totals = totals[pos]
    per = ann.loc[common]
    is_malig = per["celltype"].eq(GSE207422_MALIGNANT).to_numpy()
    is_tnk = per["celltype"].isin(GSE207422_TNK).to_numpy()
    note = per["mpr_group"].fillna("").astype(str) + "/" + per["Pathology"].astype(str)
    meta_df = pd.DataFrame(
        {
            "patient": per["patient"].to_numpy(),
            "sample": per["orig.ident"].to_numpy(),
            "cohort": "GSE207422",
            "is_tumor": np.ones(len(per), dtype=bool),
            "is_malignant": is_malig,
            "is_tnk": is_tnk,
            "celltype": per["celltype"].to_numpy(),
            "note": note.to_numpy(),
        },
        index=pd.Index(common, name="barcode"),
    )
    genes = sorted(store)
    umi = pd.DataFrame({g: store[g] for g in genes}, index=meta_df.index)
    umi["ncount"] = totals
    out = meta_df.join(umi, how="inner")
    write_cohort_panel(out, dest)


def main() -> int:
    CACHE.mkdir(parents=True, exist_ok=True)
    keep = panel_genes()
    log(f"panel requested: {len(keep)}")
    extract_gse131907(keep)
    extract_gse205335(keep)
    extract_gse207422(keep)
    log("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
