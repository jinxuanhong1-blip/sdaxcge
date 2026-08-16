#!/usr/bin/env python3
"""Registry of public KRAS/LKB1 and NSCLC GEMM / GEMM-derived ICI datasets.

Scope (this slice):
  - Kras;Lkb1 (KL / STK11-null) and Kras;Trp53 (KP) NSCLC GEMMs
  - GEMM-derived NSCLC syngeneic lines (344SQ, CMT167, KLA, Egfr-mutant)
  - A real immune-checkpoint (PD-1 / PD-L1) arm, or a genotype / mechanism
    series from a paper that studies ICI in these models
  - Not LLC-only, not SCLC GEMM (RPM / RP) — those were screened and excluded

One entry per analysis unit. Downstream scripts take group rules and contrasts
from here; they do not invent sample labels.
"""
from __future__ import annotations

from dataclasses import dataclass, field

LUAD = "LUAD"
IN_VIVO = "in_vivo_tumor"
SORTED = "sorted_tumor_cells"
CULTURED = "cultured_cells"


@dataclass
class Contrast:
    name: str
    test: str
    ref: str
    family: str          # ici | genotype | resistance | mechanism
    note: str = ""


@dataclass
class Dataset:
    key: str
    acc: str
    file: str
    model: str
    model_family: str    # KP | KL | K | 344SQ | CMT167 | EGFR | KLA
    model_class: str
    setting: str
    value_type: str      # counts | tpm | fpkm | cpm | vst | log2norm
    gene_id_type: str    # symbol | ensembl
    reader: str
    groups: list[tuple[str, str]]
    contrasts: list[Contrast] = field(default_factory=list)
    sep: str = "\t"
    comment: str | None = None
    header_row: int = 0
    gene_col: int | str = 0
    value_col: int | str | None = None
    drop_cols: list[str] = field(default_factory=list)
    member_pattern: str = ""
    col_clean: list[tuple[str, str]] = field(default_factory=list)
    group_source: str = "column"
    gsm_filter: str | None = None
    column_key_regex: str | None = None
    title_key_regex: str | None = None
    resolve_gsm: bool = True
    paired_by: str | None = None
    notes: str = ""


DATASETS: list[Dataset] = [
    # ===================== in-vivo ICI, KRAS/LKB1 ===========================
    Dataset(
        key="GSE182228_KL_aPD1",
        acc="GSE182228",
        file="GSE182228_RAW.tar",
        model="Lkb1(Stk11)-deficient murine LUAD subcutaneous tumours",
        model_family="KL",
        model_class=LUAD,
        setting=IN_VIVO,
        value_type="fpkm",
        gene_id_type="symbol",
        reader="tar_table",
        member_pattern=r"_FPKM\.txt\.gz$",
        gene_col=0,
        value_col=1,
        group_source="characteristic:treatment",
        groups=[
            (r"anti-PD-1\+Palbociclib", "aPD1_palbociclib"),
            (r"^anti-PD-1$", "aPD1"),
            (r"^Palbociclib$", "palbociclib"),
            (r"^Vehicle$", "vehicle"),
        ],
        contrasts=[
            Contrast("aPD1_vs_vehicle", "aPD1", "vehicle", "ici"),
            Contrast("aPD1_palbociclib_vs_palbociclib", "aPD1_palbociclib", "palbociclib", "ici"),
            Contrast("palbociclib_vs_vehicle", "palbociclib", "vehicle", "mechanism"),
            Contrast("aPD1_palbociclib_vs_vehicle", "aPD1_palbociclib", "vehicle", "mechanism"),
        ],
        notes="Primary KRAS/LKB1 ICI series. n=3 per arm. PMID 36871040.",
    ),
    # ===================== in-vivo ICI, KP NSCLC GEMM =======================
    Dataset(
        key="GSE114601_KP_aPD1",
        acc="GSE114601",
        file="GSE114601_counts.raw.csv.gz",
        model="Kras;Trp53 (KP) GEMM autochthonous lung tumour nodules",
        model_family="KP",
        model_class=LUAD,
        setting=IN_VIVO,
        value_type="counts",
        gene_id_type="symbol",
        reader="table",
        sep=",",
        gene_col=0,
        group_source="characteristic:treatment",
        groups=[
            (r"^anti-PD1-JQ1$", "aPD1_JQ1"),
            (r"^anti-PD1$", "aPD1"),
            (r"^JQ1$", "JQ1"),
            (r"^Vehicle$", "vehicle"),
        ],
        contrasts=[
            Contrast("aPD1_vs_vehicle", "aPD1", "vehicle", "ici"),
            Contrast("aPD1_JQ1_vs_JQ1", "aPD1_JQ1", "JQ1", "ici"),
        ],
        notes="n=2 per arm. PMID 30087114.",
    ),
    Dataset(
        key="GSE169194_KPM_total_aPD1",
        acc="GSE169194",
        file="GSE169194_annotated_raw_count.txt.gz",
        model="KrasG12D;Trp53null;Msh2null NSCLC, unfractionated tumour cells",
        model_family="KP",
        model_class=LUAD,
        setting=IN_VIVO,
        value_type="counts",
        gene_id_type="ensembl",
        reader="table",
        gene_col=0,
        group_source="characteristic:treatment",
        gsm_filter=r"total_",
        column_key_regex=r"_(159_\d+)$",
        title_key_regex=r"\[(159_\d+)\]",
        groups=[
            (r"^A2V_aPD1$", "A2V_aPD1"),
            (r"^A2V$", "A2V"),
            (r"^IgG$", "IgG"),
        ],
        contrasts=[
            Contrast("A2V_aPD1_vs_A2V", "A2V_aPD1", "A2V", "ici",
                     "anti-PD-1 added on top of anti-VEGFA/ANGPT2"),
            Contrast("A2V_vs_IgG", "A2V", "IgG", "mechanism"),
        ],
        notes="Unfractionated (DAPI-) tumour suspensions only. PMID 34380768.",
    ),
    Dataset(
        key="GSE246922_KP_ICBrelapse",
        acc="GSE246922",
        file="GSE246922_KP_RNA_counts_vst.csv.gz",
        model="KP murine lung cancer, CD45-negative tumour cells, ICB-relapsed series",
        model_family="KP",
        model_class=LUAD,
        setting=SORTED,
        value_type="vst",
        gene_id_type="ensembl",
        reader="table",
        sep=",",
        gene_col=0,
        group_source="column",
        groups=[
            (r"^ResKPlate", "relapsed_late"),
            (r"^ResResKP", "relapsed_2nd"),
            (r"^ResKP", "relapsed_1st"),
            (r"^KPy", "chronic_IFNg"),
            (r"^KP_", "parental"),
        ],
        contrasts=[
            Contrast("relapsed_1st_vs_parental", "relapsed_1st", "parental", "resistance"),
            Contrast("relapsed_2nd_vs_parental", "relapsed_2nd", "parental", "resistance"),
            Contrast("relapsed_late_vs_parental", "relapsed_late", "parental", "resistance"),
            Contrast("chronic_IFNg_vs_parental", "chronic_IFNg", "parental", "mechanism"),
        ],
        resolve_gsm=False,
        notes=("Column labels only: GEO titles cannot be mapped 1:1 onto ResKP2 / ResResKP / "
               "ResKPlate without guessing the relapse index. PMID 38215748."),
    ),
    Dataset(
        key="GSE260596_344SQ_aLAIR1",
        acc="GSE260596",
        file="GSE260596_AllSamples_Genes_ReadCounts.txt.gz",
        model="344SQ (KrasLA1/+;Trp53R172HdG/+) flank tumours, 129/Sv",
        model_family="344SQ",
        model_class=LUAD,
        setting=IN_VIVO,
        value_type="counts",
        gene_id_type="symbol",
        reader="table",
        gene_col=0,
        group_source="column",
        groups=[(r"^aLair1", "aLAIR1_aPD1"), (r"^aKLH", "control_ab_aPD1")],
        contrasts=[
            Contrast("aLAIR1_aPD1_vs_ctrlAb_aPD1", "aLAIR1_aPD1", "control_ab_aPD1", "mechanism",
                     "both arms receive anti-PD-1; contrast isolates LAIR1 blockade"),
        ],
        column_key_regex=r"_W4_(.+)$",
        title_key_regex=r"mouse (\S+)$",
        notes="All seven tumours are anti-PD-1 treated. PMID 38648067.",
    ),
    Dataset(
        key="GSE197260_EGFR_TKI_aPD1",
        acc="GSE197260",
        file="GSE197260_RNAseqTPM_MM_EGFR-TKI-CD8.txt.gz",
        model="Egfr exon-19-deletion syngeneic lung tumours, C57BL/6J",
        model_family="EGFR",
        model_class=LUAD,
        setting=IN_VIVO,
        value_type="tpm",
        gene_id_type="symbol",
        reader="table",
        gene_col=0,
        group_source="column",
        groups=[
            (r"^gef_comb_d21$", "gefitinib_then_aPD1_aVEGFR2"),
            (r"^gef_4h2_d21$", "gefitinib_then_aPD1"),
            (r"^gef_dc101_d21$", "gefitinib_then_aVEGFR2"),
            (r"^gef_vehicle_d21$", "gefitinib_then_vehicle"),
            (r"^gef_d14$", "gefitinib_d14"),
            (r"^gef_d3$", "gefitinib_d3"),
            (r"^vehicle_d3$", "vehicle"),
        ],
        contrasts=[],
        notes="One tumour per arm — descriptive only, no inferential statistics. PMID 35405020.",
    ),
    # ===================== KRAS/LKB1 genotype (no ICI on the sample) ========
    Dataset(
        key="GSE274351_K_KP_KL_adenoma",
        acc="GSE274351",
        file="GSE274351_expredata_TPM_gene.txt.gz",
        model="Laser-capture microdissected K / KP / KL lung adenomas and normal lung",
        model_family="KL",
        model_class=LUAD,
        setting=IN_VIVO,
        value_type="tpm",
        gene_id_type="ensembl",
        reader="table",
        gene_col=0,
        group_source="column",
        groups=[
            (r"^KP\d", "KP_adenoma"),
            (r"^KL\d", "KL_adenoma"),
            (r"^K\d", "K_adenoma"),
            (r"^NL\d", "normal_lung"),
        ],
        contrasts=[
            Contrast("KL_adenoma_vs_K_adenoma", "KL_adenoma", "K_adenoma", "genotype"),
            Contrast("KP_adenoma_vs_K_adenoma", "KP_adenoma", "K_adenoma", "genotype"),
            Contrast("KL_adenoma_vs_KP_adenoma", "KL_adenoma", "KP_adenoma", "genotype"),
            Contrast("K_adenoma_vs_normal_lung", "K_adenoma", "normal_lung", "genotype"),
            Contrast("KP_adenoma_vs_normal_lung", "KP_adenoma", "normal_lung", "genotype"),
            Contrast("KL_adenoma_vs_normal_lung", "KL_adenoma", "normal_lung", "genotype"),
        ],
        resolve_gsm=False,
        notes=("Column labels only: the file has 4 K and 5 KL columns while GEO lists 5 K and 4 KL "
               "samples. Genotype is taken from the data file. PMID 39186651."),
    ),
    Dataset(
        key="GSE137396_KL_vs_KP_nodules",
        acc="GSE137396",
        file="GSE137396_Raw_genetable_GEMMnodule.txt.gz",
        model="Kras;Lkb1 vs Kras;Trp53 GEMM lung tumour nodules",
        model_family="KL",
        model_class=LUAD,
        setting=IN_VIVO,
        value_type="counts",
        gene_id_type="symbol",
        reader="table",
        gene_col=0,
        group_source="column",
        groups=[(r"^KL", "KL_nodule"), (r"^KP", "KP_nodule")],
        contrasts=[Contrast("KL_vs_KP_nodule", "KL_nodule", "KP_nodule", "genotype")],
        notes="PMID 34142094.",
    ),
    Dataset(
        key="GSE175479_K_KL_AMPK",
        acc="GSE175479",
        file="GSE175479_Raw_gene_count_matrix.txt.gz",
        model="Kras-only vs Kras;Lkb1 vs Kras;AMPKa1/a2 GEMM lung tumours",
        model_family="KL",
        model_class=LUAD,
        setting=IN_VIVO,
        value_type="counts",
        gene_id_type="symbol",
        reader="table",
        gene_col="Gene_symbol",
        group_source="column",
        resolve_gsm=False,
        groups=[(r"^KL-", "KL_tumor"), (r"^KAA-", "KAA_tumor"), (r"^K-", "K_tumor")],
        contrasts=[
            Contrast("KL_vs_K", "KL_tumor", "K_tumor", "genotype"),
            Contrast("KAA_vs_K", "KAA_tumor", "K_tumor", "genotype"),
            Contrast("KL_vs_KAA", "KL_tumor", "KAA_tumor", "genotype"),
        ],
        notes="PMID 34667030. Column labels (K-1..KL-4) used as written in the count matrix.",
    ),
    Dataset(
        key="GSE274352_KP_KL_IFNB",
        acc="GSE274352",
        file="GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz",
        model="KP and KL lung adenocarcinoma lines, inducible IFN-beta",
        model_family="KL",
        model_class=LUAD,
        setting=CULTURED,
        value_type="counts",
        gene_id_type="ensembl",
        reader="table",
        gene_col=0,
        drop_cols=["external_gene_name"],
        group_source="column",
        groups=[(r"IFNBeta$", "IFNbeta"), (r"empty", "empty_vector")],
        contrasts=[Contrast("IFNbeta_vs_empty", "IFNbeta", "empty_vector", "mechanism")],
        paired_by=r"^(K[PL]\d)",
        resolve_gsm=False,
        notes="Paired within cell line (6 lines). Column labels from the processed file. PMID 39186651.",
    ),
    Dataset(
        key="GSE274352_KP_KL_STING",
        acc="GSE274352",
        file="GSE274352_normalizedcounts_genes_STING_vs_emtpy.tsv.gz",
        model="KP and KL lung adenocarcinoma lines, constitutively active STING (V154M)",
        model_family="KL",
        model_class=LUAD,
        setting=CULTURED,
        value_type="counts",
        gene_id_type="ensembl",
        reader="table",
        gene_col=0,
        drop_cols=["external_gene_name"],
        group_source="column",
        groups=[(r"STING", "STING_V154M"), (r"empty", "empty_vector")],
        contrasts=[Contrast("STING_vs_empty", "STING_V154M", "empty_vector", "mechanism")],
        paired_by=r"^(K[PL]\d)",
        resolve_gsm=False,
    ),
    Dataset(
        key="GSE295685_KL_TNG260",
        acc="GSE295685",
        file="GSE295685_TPM.txt.gz",
        model="STK11-null KL cells +/- CoREST inhibitor TNG260",
        model_family="KL",
        model_class=LUAD,
        setting=CULTURED,
        value_type="tpm",
        gene_id_type="ensembl",
        reader="table",
        gene_col=0,
        drop_cols=["gene_type", "gene_name"],
        col_clean=[(r"_TPM$", "")],
        group_source="column",
        groups=[(r"^KLTNG", "TNG260"), (r"^KLCon", "DMSO")],
        contrasts=[Contrast("TNG260_vs_DMSO", "TNG260", "DMSO", "mechanism")],
        resolve_gsm=False,
        notes="TNG260 is reported to sensitise STK11-mutant tumours to anti-PD-1. n=2 per arm. Columns from TPM file.",
    ),
    # ===================== NSCLC syngeneic (CMT167 / EGFR) ==================
    Dataset(
        key="GSE236258_CMT167_trametinib",
        acc="GSE236258",
        file="GSE236258_RAW.tar",
        model="CMT167 and KLA KRAS-mutant LUAD lines +/- trametinib",
        model_family="CMT167",
        model_class=LUAD,
        setting=CULTURED,
        value_type="counts",
        gene_id_type="ensembl",
        reader="tar_table",
        member_pattern=r"\.counts\.txt\.gz$",
        comment="#",
        gene_col=0,
        value_col=-1,
        group_source="title",
        groups=[
            (r"^CTR_", "CMT167_trametinib_resistant"),
            (r"^CT_", "CMT167_trametinib"),
            (r"^CSC_", "CMT167_control"),
            (r"^KT_", "KLA_trametinib"),
            (r"^KSC_", "KLA_control"),
        ],
        contrasts=[
            Contrast("CMT167_trametinib_vs_control", "CMT167_trametinib", "CMT167_control", "mechanism"),
            Contrast("CMT167_trametinib_resistant_vs_control", "CMT167_trametinib_resistant",
                     "CMT167_control", "resistance"),
            Contrast("KLA_trametinib_vs_control", "KLA_trametinib", "KLA_control", "mechanism"),
        ],
        notes=("CT vs CTR labels come from sample names (CTR = trametinib-resistant); "
               "GEO 'treatment' reads trametinib for both. PMID 38643157."),
    ),
    Dataset(
        key="GSE241978_CMT167_AhR",
        acc="GSE241978",
        file="GSE241978_2020-07-21_Sherr_analysis_CMT_KO_vs_Cas9Ctrl.xlsx",
        model="CMT167 cells, AhR knockout vs Cas9 control",
        model_family="CMT167",
        model_class=LUAD,
        setting=CULTURED,
        value_type="log2norm",
        gene_id_type="symbol",
        reader="xlsx",
        header_row=1,
        gene_col="Symbol",
        group_source="column",
        groups=[(r"AHR KO", "AhR_KO"), (r"Control", "cas9_control")],
        contrasts=[Contrast("AhR_KO_vs_control", "AhR_KO", "cas9_control", "mechanism")],
        resolve_gsm=False,
        notes="Normalised (log2) block of the workbook. PMID 39185148.",
    ),
    Dataset(
        key="GSE217405_EGFR_osimertinib",
        acc="GSE217405",
        file="GSE217405_mEGFR_del19.1_del19.2_L860R.1_CPM.xlsx",
        model="Egfr-mutant murine LUAD lines +/- osimertinib",
        model_family="EGFR",
        model_class=LUAD,
        setting=CULTURED,
        value_type="cpm",
        gene_id_type="symbol",
        reader="xlsx",
        header_row=0,
        gene_col=0,
        group_source="column",
        groups=[(r"Osimertinib", "osimertinib"), (r"DMSO", "DMSO")],
        contrasts=[Contrast("osimertinib_vs_DMSO", "osimertinib", "DMSO", "mechanism")],
        paired_by=r"^(mEGFR \S+) .*?(\d+Day)$",
        notes="Paired by cell line and timepoint (3 lines x 4 timepoints). PMID 36657561.",
    ),
    Dataset(
        key="GSE330658_EGFR_invivo",
        acc="GSE330658",
        file="GSE330658_RAW.tar",
        model="Egfr-mutant murine lung tumours, subcutaneous, C57BL/6J",
        model_family="EGFR",
        model_class=LUAD,
        setting=IN_VIVO,
        value_type="tpm",
        gene_id_type="symbol",
        reader="tar_xlsx",
        member_pattern=r"\.xlsx$",
        gene_col="Name",
        value_col="TPM",
        group_source="characteristic:treatment",
        groups=[
            (r"^PTX\+anti-VEGF$", "PTX_aVEGF"),
            (r"^anti-VEGF$", "aVEGF"),
            (r"^PTX$", "PTX"),
            (r"^Control$", "control"),
        ],
        contrasts=[
            Contrast("aVEGF_vs_control", "aVEGF", "control", "mechanism"),
            Contrast("PTX_vs_control", "PTX", "control", "mechanism"),
        ],
        notes=("Series text mentions anti-PD-L1 arms but only control / PTX / anti-VEGF / "
               "PTX+anti-VEGF samples are deposited. PMID 42346215."),
    ),
]

BY_KEY = {d.key: d for d in DATASETS}


def ici_datasets() -> list[Dataset]:
    return [d for d in DATASETS if any(c.family == "ici" for c in d.contrasts)]


if __name__ == "__main__":
    for d in DATASETS:
        cs = ", ".join(f"{c.name}[{c.family}]" for c in d.contrasts) or "descriptive"
        print(f"{d.key:34s} {d.acc:10s} {d.model_family:8s} {cs}")
