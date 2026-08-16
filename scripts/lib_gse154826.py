"""Shared paths, metadata loaders and gene sets for the GSE154826 TACSTD2 analysis."""

import os
import re

import numpy as np
import pandas as pd

ROOT = "/workspace"
DATA = os.path.join(ROOT, "data")
TARS = os.path.join(DATA, "geo", "tars")
PROC = os.path.join(DATA, "processed")
BATCHDIR = os.path.join(PROC, "batches")
META = os.path.join(DATA, "leader_metadata")
RESULTS = os.path.join(ROOT, "results", "hunt_gse154826")
TABLES = os.path.join(RESULTS, "tables")
FIGURES = os.path.join(RESULTS, "figures")

for _d in (PROC, BATCHDIR, TABLES, FIGURES):
    os.makedirs(_d, exist_ok=True)

# LCAM signature, verbatim from Leader et al. scripts/get_LCAM_scores.R
LCAM_HI_SUBTYPES = ["IgG_plasma", "SPP1_mac", "T_activated"]
LCAM_LO_SUBTYPES = ["B", "AM", "cDC2", "AZU1_mac", "cDC1"]
LCAM_BULK_GENES = (
    "IGHG3,IGHG4,IGHG1,MZB1,FAM92B,SPAG4,DERL3,JSRP1,TBCEL,LINC01485,SLC2A5,SPP1,"
    "CCL7,HAMP,CXCL13,ZBED2,GNG4,KRT86,TNFRSF9,LAYN,GZMB,KIR2DL4,GAPT,CTB-133G6.1,"
    "AKR1C3,RSPO3,MCEMP1,RND3,HP,FOLR3,GPD1,GS1-600G8.5,PCOLCE2,CAMP,CCL17,CD1B,"
    "CD1C,CD1E,PLD4,TRPC6,CLEC10A,FCER1A,PLA1A,CFP,CEACAM8,AZU1,PKP2,RETN,LYZ,"
    "ELANE,ATP6AP1L,RP6-159A1.4,THBS1,CFD,P2RY14,DNASE1L3,PTGER3,CST3,BATF3,CPVL,"
    "RGS18,SERPINF2,LGALS2"
).split(",")
# Gene -> subtype blocks, in the same order/lengths as the R script
LCAM_BULK_SUBTYPE = {}
for _sub, _n in zip(LCAM_HI_SUBTYPES + LCAM_LO_SUBTYPES, [10, 4, 8, 2, 10, 10, 10, 9]):
    for _g in LCAM_BULK_GENES[len(LCAM_BULK_SUBTYPE) : len(LCAM_BULK_SUBTYPE) + _n]:
        LCAM_BULK_SUBTYPE[_g] = _sub

# Single-cell (cell-frequency) LCAM definition, from scripts/figure_5abcd_s5a.R
LCAM_HI_CLUSTERS = ["T_activated", "IgG", "MoMac-II"]
LCAM_LO_CLUSTERS = ["B", "AM", "cDC2", "AZU1_mac", "Tcm/naive_II", "cDC1"]
NORM_GROUPS = ["T", "B&plasma", "MNP", "lin_neg"]

MARKERS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7", "ELF3", "CLDN4",
                   "SLPI", "MUC1", "SFTPC", "SFTPB", "SFTPA1", "NAPSA", "SCGB1A1",
                   "SCGB3A2", "KRT5", "KRT17", "AGER", "FOXJ1", "TPPP3"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "RAMP2", "FLT1", "EGFL7", "CCL21"],
    "fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "ACTA2",
                   "MYH11", "FN1", "SPARC"],
    "t_nk": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC2", "IL7R", "CD8A", "CD8B",
             "CD4", "NKG7", "GNLY", "KLRD1", "GZMB", "GZMK"],
    "b_plasma": ["MS4A1", "CD79A", "CD79B", "BANK1", "MZB1", "JCHAIN", "DERL3",
                 "IGHG1", "IGHM", "IGKC"],
    "myeloid": ["LYZ", "CD68", "CD14", "FCGR3A", "AIF1", "C1QA", "C1QB", "MARCO",
                "MRC1", "SPP1", "ITGAX", "FCN1", "S100A8", "S100A9", "CD1C", "CLEC9A",
                "LAMP3"],
    "mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "pdc": ["LILRA4", "IRF7", "CLEC4C", "GZMB"],
    "erythroid": ["HBB", "HBA1", "HBA2"],
    "cycling": ["MKI67", "TOP2A", "STMN1"],
}
# ADC / surface-target genes reported alongside TACSTD2 for context
TARGET_GENES = ["TACSTD2", "PTPRC", "EPCAM", "ERBB2", "EGFR", "MSLN", "CEACAM5",
                "FOLR1", "NECTIN4", "MUC1", "CD274", "CDH1"]

PANEL_GENES = sorted({g for v in MARKERS.values() for g in v} | set(TARGET_GENES))

CD8_ADT = "CD8"
GATE_LINEAGE = "epi_endo_fibro_doublet"


def load_annots_list():
    a = pd.read_csv(os.path.join(META, "annots_list.csv"))
    a["sub_lineage"] = a["sub_lineage"].fillna("")
    # figure_5abcd_s5a.R: MNP clusters are all normalised together (DC included)
    a["norm_group_fig5"] = np.where(a["lineage"] == "MNP", "MNP", a["norm_group"])
    return a


def load_cell_metadata():
    p_gz = os.path.join(META, "cell_metadata.csv.gz")
    p = os.path.join(META, "cell_metadata.csv")
    cm = pd.read_csv(p_gz if os.path.exists(p_gz) else p)
    cm["barcode"] = cm["cell_ID"].str.split("_", n=1).str[1]
    return cm


def load_sample_annots():
    """GEO sample sheet joined with the published Table S1."""
    geo_p = os.path.join(META, "sample_annots.csv")
    if not os.path.exists(geo_p):
        geo_p = os.path.join(DATA, "geo", "sample_annots.csv")
    geo = pd.read_csv(geo_p)
    s1 = pd.read_csv(os.path.join(META, "table_s1_sample_table.csv"))
    s1 = s1[["sample_ID", "prep", "Use.in.Clustering.Model.", "biopsy_site"]].rename(
        columns={"prep": "prep_s1", "Use.in.Clustering.Model.": "use_in_model"}
    )
    return geo.merge(s1, on="sample_ID", how="left")


def batch_tar(batch):
    return os.path.join(TARS, f"GSE154826_amp_batch_ID_{batch}.tar.gz")


def list_batches():
    out = []
    for f in os.listdir(TARS):
        m = re.match(r"GSE154826_amp_batch_ID_(\d+)\.tar\.gz$", f)
        if m:
            out.append(int(m.group(1)))
    return sorted(out)
