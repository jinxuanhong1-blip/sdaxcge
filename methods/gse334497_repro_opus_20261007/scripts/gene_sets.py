"""Gene-set panels for the GSE334497 reanalysis.

The three "requested" panels (APM_core, IFN_core, Immune_core) are the human
lists given in the analysis brief, translated to mouse by namesake symbol
(Ensembl 102 / GRCm38, the submitters' annotation). Where no namesake exists
the substitution is listed in HUMAN_TO_MOUSE_NOTES. Every other panel was fixed
before any statistics were computed; nothing was added or dropped afterwards.
Members missing from the deposited matrix are written to
tables/geneset_membership.csv.
"""

HUMAN_TO_MOUSE_NOTES = {
    "HLA-A/B/C": "H2-K1, H2-D1 (classical MHC-Ia; BALB/c H-2d K/D/L reads map to these "
                 "GRCm38 B6 loci; H2-L has no GRCm38 gene model)",
    "CD8B": "Cd8b1",
    "GBP1": "Gbp2b (MGI symbol of the gene formerly named mouse Gbp1)",
    "OAS1": "Oas1a",
    "MX1": "Mx1",
}

APM_CORE = ["B2m", "H2-K1", "H2-D1", "Tap1", "Tap2", "Psmb8", "Psmb9", "Nlrc5", "Tapbp"]

IFN_CORE = [
    "Stat1", "Irf1", "Cxcl9", "Cxcl10", "Cxcl11",
    "Gbp2b", "Gbp2", "Gbp4", "Irf7", "Isg15", "Mx1", "Oas1a",
    "Ifit1", "Ifit3", "Ido1",
]

IMMUNE_CORE = ["Cd8a", "Cd8b1", "Gzmb", "Prf1", "Nkg7", "Ifng"]

CLAUDIN_FOCAL = ["Cldn1", "Cldn3", "Cldn4", "Cldn7"]

# Paper's own RNA-seq tight-junction read-out (supplementary figure S4D).
PAPER_TJ_TRIO = ["Cldn1", "Cldn7", "Ocln"]

TJ_CLAUDIN_PROGRAM = [
    "Cldn1", "Cldn3", "Cldn4", "Cldn7",
    "Ocln", "Tjp1", "Tjp2", "Tjp3", "F11r", "Marveld2", "Marveld3",
    "Cgn", "Cgnl1", "Crb3", "Pard3", "Pard6b", "Patj", "Llgl2",
]

CD8_T = ["Cd8a", "Cd8b1", "Cd3e", "Cd3d", "Cd3g"]

EFFECTOR = ["Gzma", "Gzmb", "Gzmk", "Prf1", "Nkg7", "Ifng", "Fasl"]

LEUKOCYTE = ["Ptprc", "Coro1a", "Laptm5", "Lcp1", "Cd53", "Itgb2", "Cd48", "Lcp2"]

NK = ["Ncr1", "Klrb1c", "Klrd1", "Klrk1", "Klrc1", "Il2rb"]

MHC_II = ["H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "H2-DMb1", "Cd74", "Ciita"]

CHECKPOINT = ["Pdcd1", "Cd274", "Pdcd1lg2", "Ctla4", "Lag3", "Havcr2", "Tigit"]

PANELS = {
    "APM_core": APM_CORE,
    "IFN_core": IFN_CORE,
    "Immune_core": IMMUNE_CORE,
    "Claudin_1_3_4_7": CLAUDIN_FOCAL,
    "Paper_TJ_trio_Cldn1_Cldn7_Ocln": PAPER_TJ_TRIO,
    "TJ_claudin_program": TJ_CLAUDIN_PROGRAM,
    "CD8_T": CD8_T,
    "Effector": EFFECTOR,
    "Leukocyte": LEUKOCYTE,
    "NK": NK,
    "MHC_II": MHC_II,
    "Checkpoint": CHECKPOINT,
}

# Post hoc tissue-contamination markers, defined after Krt1 / Krtdap / Lypd5 appeared at the
# top of the KO-vs-WT DESeq2 list. Suprabasal epidermis and skeletal muscle (panniculus
# carnosus / chest wall) are not 4T1 products. Diagnostics only, never primary end points.
SKIN_EPIDERMIS = ["Krt1", "Krt10", "Krt77", "Lor", "Flg", "Flg2", "Krtdap", "Dmkn",
                  "Dsg1a", "Dsc3", "Lce1a1", "Sbsn", "Lypd5", "Cdsn", "Spink5"]
SKIN_FLAG_GENES = ["Krt1", "Flg2", "Krt77", "Dsg1a", "Dsc3"]
SKIN_FLAG_MIN_SUM = 50  # summed DESeq2-normalized counts of SKIN_FLAG_GENES
SKELETAL_MUSCLE = ["Acta1", "Ckm", "Myh1", "Myh2", "Myh4", "Myl1", "Mylpf", "Tnnt3",
                   "Tnni2", "Tnnc2", "Actn3", "Des"]

FOCAL_GENES = [
    "Tacstd2", "Epcam",
    "Cldn1", "Cldn2", "Cldn3", "Cldn4", "Cldn5", "Cldn7", "Cldn10", "Cldn12",
    "Cldn15", "Cldn20", "Cldn23",
    "Ocln", "Tjp1", "Tjp2", "Tjp3", "F11r", "Marveld2", "Cgn", "Crb3",
    "Ptprc", "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Gzma", "Gzmb", "Prf1", "Nkg7", "Ifng",
    "Ncr1", "Klrd1", "Pdcd1", "Cd274", "Ctla4", "Lag3", "Havcr2",
    "B2m", "H2-K1", "H2-D1", "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9", "Nlrc5",
    "Stat1", "Irf1", "Irf7", "Cxcl9", "Cxcl10", "Cxcl11", "Gbp2", "Gbp2b", "Gbp4",
    "Isg15", "Mx1", "Oas1a", "Ifit1", "Ifit3", "Ido1",
    "Mki67",
]

# Paper figure 4A / supplementary S4B labels -> candidate MSigDB 2024.1.Hs set
# names (GSEA v4.3.2 remaps mouse identifiers to human orthologs). Where a label
# is ambiguous every candidate is tested and reported.
PAPER_FIG4A = [
    # (label, reported NES, reported nominal p, [candidate set names])
    ("Cell_cell contact zone", 1.61, 0.008, ["GOCC_CELL_CELL_CONTACT_ZONE"]),
    ("Tight junction", 1.49, 0.009, ["KEGG_TIGHT_JUNCTION", "GOCC_TIGHT_JUNCTION"]),
    ("Leukocyte degranulation", -1.51, 0.015, ["GOBP_LEUKOCYTE_DEGRANULATION"]),
    ("Inflammation pathway", -1.55, 0.025, ["BIOCARTA_INFLAM_PATHWAY"]),
    ("Leukocyte mediated cytotoxicity", -1.62, 0.025, ["GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY"]),
    ("Immune response to tumor cell", -1.63, 0.014, ["GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL"]),
    ("Lymphocyte costimulation", -1.64, 0.017, ["GOBP_LYMPHOCYTE_COSTIMULATION"]),
    ("PD1 signaling", -1.64, 0.027, ["REACTOME_PD_1_SIGNALING"]),
    ("NK cell mediated cytotoxicity", -1.65, 0.016, ["KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY"]),
    ("T cell mediated cytotoxicity", -1.66, 0.025, ["GOBP_T_CELL_MEDIATED_CYTOTOXICITY"]),
    ("Th1 cytotoxic module", -1.68, 0.018, ["BOSCO_TH1_CYTOTOXIC_MODULE"]),
    ("Positive regulation of cell killing", -1.72, 0.015, ["GOBP_POSITIVE_REGULATION_OF_CELL_KILLING"]),
]

CONTEXT_SETS = [
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "REACTOME_TIGHT_JUNCTION_INTERACTIONS",
    "GOBP_TIGHT_JUNCTION_ORGANIZATION",
    "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
    "REACTOME_INTERFERON_GAMMA_SIGNALING",
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
]
