#!/usr/bin/env python3
"""Extract CLDN4 (ENSG00000189143) and clinical labels from PredictIO/ORCESTRA ICB TSVs.

Source: Zenodo 10.5281/zenodo.7199344 (per-study TSV extract of ORCESTRA ICB objects).
CLDN4 is looked up by Ensembl ID, not by symbol, so panel studies that lack the
gene are recorded as missing rather than silently skipped.
"""
from __future__ import annotations

import csv
import io
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ICB = ROOT / "data" / "raw" / "icb"
DER = ROOT / "data" / "derived"
OUT = ROOT / "results" / "claim_B5"
DER.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

CLDN4 = "ENSG00000189143"

# In-scope indications named in the claim. Out-of-scope studies are still
# extracted so the inventory is complete, then excluded from the primary meta.
STUDIES = {
    "Braun": {
        "expr": "ICB_Braun_expr.tsv",
        "genes": "ICB_Braun_expr_genes.tsv",
        "meta": "ICB_Braun_metadata.tsv",
        "indication": "RCC",
        "io_class": "PD-1/PD-L1",
        "citation": "Braun et al., Nat Med 2020 (CheckMate 025/010)",
        "pmid": "32472114",
        "in_scope": True,
    },
    "Miao1": {
        "expr": "ICB_Miao1_expr.tsv",
        "genes": "ICB_Miao1_expr_genes.tsv",
        "meta": "ICB_Miao1_metadata.tsv",
        "indication": "RCC",
        "io_class": "PD-1/PD-L1",
        "citation": "Miao et al., Science 2018",
        "pmid": "29301960",
        "in_scope": True,
    },
    "Shiuan": {
        "expr": "ICB_Shiuan_expr.tsv",
        "genes": "ICB_Shiuan_expr_genes.tsv",
        "meta": "ICB_Shiuan_metadata.tsv",
        "indication": "RCC",
        "io_class": "PD-1/PD-L1",
        "citation": "Shiuan et al., JCI Insight 2020",
        "pmid": "32315290",
        "in_scope": True,
    },
    "Mariathasan": {
        "expr": "ICB_Mariathasan_expr.tsv",
        "genes": "ICB_Mariathasan_expr_genes.tsv",
        "meta": "ICB_Mariathasan_metadata.tsv",
        "indication": "urothelial",
        "io_class": "PD-1/PD-L1",
        "citation": "Mariathasan et al., Nature 2018 (IMvigor210)",
        "pmid": "29443960",
        "in_scope": True,
        "note": "IMvigor210 atezolizumab UC; cancer_type field is biopsy site, not a second primary.",
    },
    "Snyder": {
        "expr": "ICB_Snyder_expr.tsv",
        "genes": "ICB_Snyder_expr_genes.tsv",
        "meta": "ICB_Snyder_metadata.tsv",
        "indication": "urothelial",
        "io_class": "PD-1/PD-L1",
        "citation": "Snyder et al., PLoS Med 2017",
        "pmid": "28542453",
        "in_scope": True,
    },
    "Jung": {
        "expr": "ICB_Jung_expr_gene_tpm.tsv",
        "genes": "ICB_Jung_expr_gene_tpm_genes.tsv",
        "meta": "ICB_Jung_metadata.tsv",
        "indication": "NSCLC",
        "io_class": "PD-1/PD-L1",
        "citation": "Jung et al., Nat Commun 2019",
        "pmid": "31383964",
        "in_scope": True,
    },
    "Fumet1": {
        "expr": "ICB_Fumet1_expr_gene_tpm.tsv",
        "genes": None,
        "meta": "ICB_Fumet1_metadata.tsv",
        "indication": "NSCLC",
        "io_class": "PD-1/PD-L1",
        "citation": "Limagne / Fumet NSCLC anti-PD-1 (ORCESTRA ICB_Fumet1 / ICB_Limagne1)",
        "pmid": "",
        "in_scope": True,
        "alias": "Limagne1",
    },
    "Fumet2": {
        "expr": "ICB_Fumet2_expr_gene_tpm.tsv",
        "genes": None,
        "meta": "ICB_Fumet2_metadata.tsv",
        "indication": "NSCLC",
        "io_class": "PD-1/PD-L1",
        "citation": "Limagne / Fumet NSCLC anti-PD-1/PD-L1 (ORCESTRA ICB_Fumet2 / ICB_Limagne2)",
        "pmid": "",
        "in_scope": True,
        "alias": "Limagne2",
    },
    "Gide": {
        "expr": "ICB_Gide_expr_gene_tpm.tsv",
        "genes": "ICB_Gide_expr_gene_tpm_genes.tsv",
        "meta": "ICB_Gide_metadata.tsv",
        "indication": "melanoma",
        "io_class": "PD-1/PD-L1",
        "citation": "Gide et al., Cancer Cell 2019",
        "pmid": "30753825",
        "in_scope": True,
    },
    "Hugo": {
        "expr": "ICB_Hugo_expr_gene_tpm.tsv",
        "genes": "ICB_Hugo_expr_gene_tpm_genes.tsv",
        "meta": "ICB_Hugo_metadata.tsv",
        "indication": "melanoma",
        "io_class": "PD-1/PD-L1",
        "citation": "Hugo et al., Cell 2016",
        "pmid": "26997480",
        "in_scope": True,
    },
    "Liu": {
        "expr": "ICB_Liu_expr.tsv",
        "genes": "ICB_Liu_expr_genes.tsv",
        "meta": "ICB_Liu_metadata.tsv",
        "indication": "melanoma",
        "io_class": "PD-1/PD-L1",
        "citation": "Liu et al., Nat Med 2019",
        "pmid": "31792460",
        "in_scope": True,
    },
    "Riaz": {
        "expr": "ICB_Riaz_expr_gene_tpm.tsv",
        "genes": "ICB_Riaz_expr_gene_tpm_genes.tsv",
        "meta": "ICB_Riaz_metadata.tsv",
        "indication": "melanoma",
        "io_class": "PD-1/PD-L1",
        "citation": "Riaz et al., Cell 2017",
        "pmid": "29033130",
        "in_scope": True,
    },
    "Puch": {
        "expr": "ICB_Puch_expr.tsv",
        "genes": "ICB_Puch_expr_genes.tsv",
        "meta": "ICB_Puch_metadata.tsv",
        "indication": "melanoma",
        "io_class": "PD-1/PD-L1",
        "citation": "ORCESTRA ICB_Puch (melanoma PD-1/PD-L1)",
        "pmid": "",
        "in_scope": True,
    },
    "Van_Allen": {
        "expr": "ICB_Van_Allen_expr.tsv",
        "genes": "ICB_Van_Allen_expr_genes.tsv",
        "meta": "ICB_Van_Allen_metadata.tsv",
        "indication": "melanoma",
        "io_class": "CTLA-4",
        "citation": "Van Allen et al., Science 2015",
        "pmid": "26359337",
        "in_scope": True,
    },
    "Nathanson": {
        "expr": "ICB_Nathanson_expr.tsv",
        "genes": "ICB_Nathanson_expr_genes.tsv",
        "meta": "ICB_Nathanson_metadata.tsv",
        "indication": "melanoma",
        "io_class": "CTLA-4",
        "citation": "Nathanson et al., Cancer Immunol Res 2017",
        "pmid": "28039162",
        "in_scope": True,
    },
    "Hwang": {
        "expr": "ICB_Hwang_expr.tsv",
        "genes": "ICB_Hwang_expr_genes.tsv",
        "meta": "ICB_Hwang_metadata.tsv",
        "indication": "NSCLC",
        "io_class": "PD-1/PD-L1",
        "citation": "Hwang et al. (immune-gene panel; CLDN4 absent)",
        "pmid": "",
        "in_scope": True,
    },
    "Jerby_Arnon": {
        "expr": "ICB_Jerby_Arnon_expr.tsv",
        "genes": "ICB_Jerby_Arnon_expr_genes.tsv",
        "meta": "ICB_Jerby_Arnon_metadata.tsv",
        "indication": "melanoma",
        "io_class": "PD-1/PD-L1",
        "citation": "Jerby-Arnon et al., Cell 2018 (signature panel; CLDN4 absent)",
        "pmid": "30388455",
        "in_scope": True,
    },
    "Roh": {
        "expr": "ICB_Roh_expr.tsv",
        "genes": "ICB_Roh_expr_genes.tsv",
        "meta": "ICB_Roh_metadata.tsv",
        "indication": "melanoma",
        "io_class": "CTLA-4",
        "citation": "Roh et al., Sci Transl Med 2017 (panel; CLDN4 absent)",
        "pmid": "28448585",
        "in_scope": True,
    },
    "Kim": {
        "expr": "ICB_Kim_expr_gene_tpm.tsv",
        "genes": "ICB_Kim_expr_gene_tpm_genes.tsv",
        "meta": "ICB_Kim_metadata.tsv",
        "indication": "gastric",
        "io_class": "PD-1/PD-L1",
        "citation": "Kim et al., Nat Med 2018",
        "pmid": "30013197",
        "in_scope": False,
    },
    "Padron": {
        "expr": "ICB_Padron_expr.tsv",
        "genes": "ICB_Padron_expr_genes.tsv",
        "meta": "ICB_Padron_metadata.tsv",
        "indication": "pancreas",
        "io_class": "PD-1/PD-L1",
        "citation": "Padron et al. (ORCESTRA ICB_Padron)",
        "pmid": "",
        "in_scope": False,
    },
    "VanDenEnde": {
        "expr": "ICB_VanDenEnde_expr.tsv",
        "genes": "ICB_VanDenEnde_expr_genes.tsv",
        "meta": "ICB_VanDenEnde_metadata.tsv",
        "indication": "esophageal",
        "io_class": "PD-1/PD-L1",
        "citation": "van den Ende et al. (CLDN4 absent from this object's annotation)",
        "pmid": "",
        "in_scope": False,
    },
}


def _parse_header(path: Path) -> list[str]:
    with path.open() as fh:
        rec = next(csv.reader(fh, delimiter="\t", quotechar='"'))
    return rec


def _find_cldn4_row(genes_path: Path | None, expr_path: Path) -> tuple[int | None, str]:
    """Return 0-based data-row index of CLDN4 in the expression matrix, or None."""
    if genes_path is None:
        # First column of the expression file is the Ensembl ID (Fumet objects).
        with expr_path.open() as fh:
            next(fh)  # sample header
            for i, line in enumerate(fh):
                gid = line.split("\t", 1)[0].strip().strip('"').split(".")[0]
                if gid == CLDN4:
                    return i, gid
        return None, ""

    genes = pd.read_csv(genes_path, sep="\t")
    col = None
    for c in ("gene_id_no_ver", "gene_id", "gene_name"):
        if c in genes.columns:
            col = c
            break
    if col is None:
        raise ValueError(f"no gene id column in {genes_path}")
    series = genes[col].astype(str).str.replace(r"\.\d+$", "", regex=True)
    hits = np.flatnonzero(series.to_numpy() == CLDN4)
    if len(hits) == 0 and "gene_name" in genes.columns:
        hits = np.flatnonzero(genes["gene_name"].astype(str).to_numpy() == "CLDN4")
    if len(hits) == 0:
        return None, ""
    return int(hits[0]), str(genes.iloc[int(hits[0])].get("gene_id", CLDN4))


def _looks_like_gene_id(v: str) -> bool:
    s = v.strip().strip('"')
    return s.startswith("ENSG") or s.startswith("ENST") or s.startswith("NM_")


def _read_expr_row(expr_path: Path, row_index: int, n_samples: int) -> np.ndarray:
    with expr_path.open() as fh:
        next(fh)
        for i, line in enumerate(fh):
            if i != row_index:
                continue
            rec = next(csv.reader(io.StringIO(line), delimiter="\t", quotechar='"'))
            if rec and _looks_like_gene_id(rec[0]):
                rec = rec[1:]
            if len(rec) != n_samples:
                # some objects pad an extra index column that is not a gene id
                if len(rec) == n_samples + 1:
                    rec = rec[1:]
            if len(rec) != n_samples:
                raise ValueError(
                    f"{expr_path.name} row {row_index}: got {len(rec)} values, expected {n_samples}"
                )
            out = []
            for v in rec:
                if v in ("NA", "NaN", "", "nan"):
                    out.append(np.nan)
                else:
                    out.append(float(v))
            return np.asarray(out, dtype=np.float64)
    raise IndexError(f"row {row_index} missing in {expr_path}")


def extract_study(name: str, spec: dict) -> tuple[pd.DataFrame, dict]:
    expr_path = ICB / spec["expr"]
    genes_path = ICB / spec["genes"] if spec["genes"] else None
    meta = pd.read_csv(ICB / spec["meta"], sep="\t")
    samples = _parse_header(expr_path)
    if samples and _looks_like_gene_id(samples[0]):
        samples = samples[1:]
    row, gid = _find_cldn4_row(genes_path, expr_path)
    info = {
        "study": name,
        "alias": spec.get("alias", name),
        "indication": spec["indication"],
        "io_class": spec["io_class"],
        "in_scope": spec["in_scope"],
        "citation": spec["citation"],
        "pmid": spec["pmid"],
        "n_meta": int(len(meta)),
        "n_expr": int(len(samples)),
        "cldn4_present": row is not None,
        "cldn4_id": gid,
        "note": spec.get("note", ""),
    }
    if row is None:
        info["n_overlap"] = 0
        info["n_with_response"] = 0
        return pd.DataFrame(), info

    values = _read_expr_row(expr_path, row, len(samples))
    if len(values) != len(samples):
        raise ValueError(f"{name}: expr row length {len(values)} != n samples {len(samples)}")
    expr_df = pd.DataFrame({"patientid": samples, "CLDN4": values})
    meta = meta.copy()
    meta["patientid"] = meta["patientid"].astype(str)
    expr_df["patientid"] = expr_df["patientid"].astype(str)
    merged = expr_df.merge(meta, on="patientid", how="inner")
    keep = [
        "patientid",
        "CLDN4",
        "response",
        "recist",
        "cancer_type",
        "tissueid",
        "treatment",
        "treatmentid",
        "survival_time_os",
        "event_occurred_os",
        "survival_time_pfs",
        "event_occurred_pfs",
    ]
    for c in keep:
        if c not in merged.columns:
            merged[c] = np.nan
    out = merged[keep].copy()
    out.insert(0, "study", name)
    out.insert(1, "indication", spec["indication"])
    out.insert(2, "io_class", spec["io_class"])
    out.insert(3, "in_scope", spec["in_scope"])
    info["n_overlap"] = int(len(out))
    info["n_with_response"] = int(out["response"].isin(["R", "NR"]).sum())
    info["cldn4_median"] = float(np.nanmedian(out["CLDN4"]))
    info["cldn4_iqr"] = float(np.nanpercentile(out["CLDN4"], 75) - np.nanpercentile(out["CLDN4"], 25))
    info["cldn4_min"] = float(np.nanmin(out["CLDN4"]))
    info["cldn4_max"] = float(np.nanmax(out["CLDN4"]))
    return out, info


def main() -> None:
    frames = []
    inventory = []
    for name, spec in STUDIES.items():
        df, info = extract_study(name, spec)
        inventory.append(info)
        if len(df):
            frames.append(df)
        print(
            f"{name:14s} cldn4={info['cldn4_present']} overlap={info['n_overlap']:3d} "
            f"resp={info['n_with_response']:3d}",
            flush=True,
        )

    inv = pd.DataFrame(inventory)
    inv.to_csv(OUT / "cohort_inventory.tsv", sep="\t", index=False)
    if frames:
        all_df = pd.concat(frames, ignore_index=True)
        all_df.to_csv(DER / "cldn4_clinical.tsv", sep="\t", index=False)
        all_df.to_csv(OUT / "cldn4_clinical.tsv", sep="\t", index=False)
        print(f"wrote {len(all_df)} sample rows")
    else:
        print("no samples extracted")


if __name__ == "__main__":
    main()
