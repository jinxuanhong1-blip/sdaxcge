"""Assemble patient-level tumour-segment expression + outcome tables for three ICI cohorts.

Yale (primary)  : GSE271689 DCC files reprocessed here (tumour / PanCK segments),
                  outcomes taken from Source Data Fig. 6b of Aung et al., Nat Genet 2025.
Greek (repl. 1) : Source Data Fig. 6d of the same paper (tumour compartment, not in GEO).
UQ (repl. 2)    : Source Data Fig. 6c (GeoMx CTA panel, GSE221733; no CLDN4 on the panel).

Outputs go to results/hunt_gse271689/data/.
"""

import os

import numpy as np
import pandas as pd
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results", "hunt_gse271689", "data")
os.makedirs(OUT, exist_ok=True)

SRC_XLSX = os.path.join(DATA, "paper", "41588_2025_2351_MOESM10_ESM.xlsx")

YALE_CLINICAL_COLS = [
    "Spot_ID", "No", "Core.Type", "Biopsy_ITx", "Smoking_Status", "PD_L1_TPS_IC_other_biopsy",
    "Mutation_Status", "Stage_at_biopsy", "Stage_at_ITx", "Site_of_biopsy", "Histotype",
    "Metastatic_Site_at_ITx", "Line_ITx", "Agent_ITx", "Concurrent_treatment", "OR_1st_scan",
    "BOR", "response", "Cycles_ITx", "PFS_Days", "PFS_Index", "PFS_2Years_months",
    "PFS_2Years_Index", "PFS_5Years_months", "PFS_5Years_Index", "OS_Days", "OS_Index",
    "OS_5Years_months", "OS_5Years_Index", "OS_2Years_months", "OS_2Years_Index",
    "Previous_surgery", "Previous_RT", "Previous_adjuvant_chemotherapy", "Previous_line",
]


def read_sheet(sheet):
    wb = openpyxl.load_workbook(SRC_XLSX, read_only=True)
    ws = wb[sheet]
    it = ws.iter_rows(values_only=True)
    hdr = list(next(it))
    rows = [list(r) for r in it]
    wb.close()
    hdr = [h if h is not None else "row_index" for h in hdr]
    return pd.DataFrame(rows, columns=hdr)


def yale_expression():
    """Per-patient tumour-segment expression from our own DCC reprocessing."""
    ann = pd.read_csv(os.path.join(DATA, "processed", "aoi_annotation.tsv"), sep="\t")
    expr = pd.read_csv(os.path.join(DATA, "processed", "q3_log2.tsv.gz"), sep="\t", index_col=0)
    ck = ann[(ann["qc_pass"]) & (ann["segment"] == "CK")].copy()
    ck["spot_id"] = ck["spot_id"].astype(str)
    # mean of log2 Q3 expression over the (up to 4) tumour AOIs of a patient
    patient = expr[ck["gsm"]].T.groupby(ck.set_index("gsm")["spot_id"]).mean()
    n_aoi = ck.groupby("spot_id").size().rename("n_tumour_aoi")
    return patient, n_aoi, ck, ann, expr


def main():
    # ---------------- Yale ----------------
    patient, n_aoi, ck, ann, expr = yale_expression()
    fig6b = read_sheet("source_data_Figure_6b")
    clin = fig6b[[c for c in YALE_CLINICAL_COLS if c in fig6b.columns]].copy()
    clin["Spot_ID"] = clin["Spot_ID"].astype(str)
    clin = clin.rename(columns={"Spot_ID": "spot_id"})

    missing = sorted(set(clin["spot_id"]) - set(patient.index))
    if missing:
        print("clinical spot ids without GEO tumour AOIs:", missing)

    yale = clin.merge(patient, left_on="spot_id", right_index=True, how="inner")
    yale = yale.merge(n_aoi, left_on="spot_id", right_index=True, how="left")
    yale.to_csv(os.path.join(OUT, "yale_tumour_patient_level.csv.gz"), index=False)
    print("Yale patients with tumour expression + outcome:", yale.shape[0])

    # patient-level expression for every GEO patient (61), for reference / QC
    patient.assign(n_tumour_aoi=n_aoi).to_csv(os.path.join(OUT, "gse271689_tumour_patient_expression.csv.gz"))

    # CD45 and CD68 segments, same aggregation, for the immune-compartment sensitivity analysis
    for seg in ["CD45", "CD68"]:
        s = ann[(ann["qc_pass"]) & (ann["segment"] == seg)].copy()
        s["spot_id"] = s["spot_id"].astype(str)
        p = expr[s["gsm"]].T.groupby(s.set_index("gsm")["spot_id"]).mean()
        p.to_csv(os.path.join(OUT, f"gse271689_{seg}_patient_expression.csv.gz"))

    # ---------------- Greek ----------------
    fig6d = read_sheet("source_data_Figure_6d")
    fig6d.to_csv(os.path.join(OUT, "greek_tumour_roi_level.csv.gz"), index=False)
    print("Greek tumour ROIs:", fig6d.shape[0])

    # ---------------- UQ ----------------
    fig6c = read_sheet("source_data_Figure_6c")
    fig6c.to_csv(os.path.join(OUT, "uq_tumour_patient_level.csv.gz"), index=False)
    print("UQ tumour samples:", fig6c.shape[0])

    # AOI-level Yale table (tumour) for the clustered-Cox sensitivity analysis
    genes_needed = sorted(set(
        ["TACSTD2", "CLDN4", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10",
         "CXCL11", "CXCL13", "HLA-DRA", "STAT1", "PSMB10", "IDO1", "LAG3", "CD274", "PDCD1",
         "CCL5", "NKG7", "KLRD1", "CMKLR1", "TIGIT", "CD27", "CD276", "CD3D", "CD3E", "CD2",
         "EPCAM", "KRT19", "PTPRC", "MKI67", "VIM", "CDH1"]
    ) & set(expr.index))
    aoi = ck[["gsm", "spot_id", "segment", "plate", "well", "total_counts", "neg_geomean", "loq"]].copy()
    aoi = aoi.merge(expr.loc[genes_needed, ck["gsm"]].T, left_on="gsm", right_index=True)
    aoi.to_csv(os.path.join(OUT, "yale_tumour_aoi_level.csv.gz"), index=False)
    print("Yale tumour AOIs:", aoi.shape[0])


if __name__ == "__main__":
    main()
