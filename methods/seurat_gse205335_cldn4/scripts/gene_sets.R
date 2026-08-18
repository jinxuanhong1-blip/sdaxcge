# Locked gene sets used in prior CLDN4-only IFN/MHC/TJ work (GSE207422 TJ/IFN
# and the pair malignant-family extra). CLDN4 is the splitter and is held out
# of the TJ module. Not a dual-high TACSTD2×CLDN4 list. Not invented here.

# Tight-junction strand / polarity. CLDN4 intentionally absent.
TJ_NO_CLDN4 <- c(
  "CLDN1", "CLDN3", "CLDN7", "OCLN", "TJP1", "TJP2", "F11R", "PARD3",
  "MARVELD2", "CGN", "CRB3", "JAM3"
)

# Compact classical MHC-I / APM panel (CUSTOM_MHC_I_ANTIGEN_PRESENTATION).
MHC1_APM <- c(
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2",
  "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "NLRC5",
  "ERAP1", "ERAP2", "CALR", "PDIA3", "CANX", "IRF1"
)

# Type-I ISG core from the public CLDN4/IFN alignment (align_cldn4_ifn).
ISG_CORE <- c(
  "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
  "IFIT5", "IFITM1", "IFITM2", "IFITM3", "MX1", "MX2", "OAS1", "OAS2",
  "OAS3", "OASL", "RSAD2", "USP18", "BST2", "XAF1", "STAT1", "STAT2",
  "IRF7", "IRF9", "DDX58", "IFIH1", "SAMD9", "SAMD9L", "HERC5", "HERC6",
  "EPSTI1", "CMPK2", "PARP9", "DTX3L", "LY6E", "SP100", "SP110", "PLSCR1"
)

# Specificity control: not everything should move with CLDN4.
CTRL_OXPHOS <- c(
  "NDUFA1", "NDUFB3", "COX5A", "COX7A2", "UQCRC1", "SDHA",
  "ATP5F1A", "ATP5F1B", "CYCS", "VDAC1"
)

MODULES <- list(
  ifn_isg = ISG_CORE,
  mhc1_apm = MHC1_APM,
  tj_no_cldn4 = TJ_NO_CLDN4,
  ctrl_oxphos = CTRL_OXPHOS
)

MODULE_LABEL <- c(
  ifn_isg = "IFN (ISG core)",
  mhc1_apm = "MHC-I/APM",
  tj_no_cldn4 = "TJ (CLDN4 held out)",
  ctrl_oxphos = "OXPHOS control"
)
