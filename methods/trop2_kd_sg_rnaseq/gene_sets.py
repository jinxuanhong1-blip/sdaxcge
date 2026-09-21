"""Pre-specified genes for the TACSTD2 KD / sacituzumab RNA-seq screen.

Human symbols are the row labels. Mouse orthologs are the symbols (or, for
OAS1, Oas1a) used to pull the 4T1 matrix. STING1/CGAS aliases are accepted
when a matrix uses the older name.
"""

# Enzymatic core of classical non-homologous end joining.
NHEJ = ["XRCC6", "XRCC5", "PRKDC", "LIG4", "XRCC4", "NHEJ1", "DCLRE1C", "PAXX"]

# cGAS–STING signaling axis, without IFN target genes.
STING = ["CGAS", "STING1", "TBK1", "IRF3"]

# Compact interferon-stimulated gene set. IRF7 sits here, not in STING.
IFN = [
    "STAT1",
    "STAT2",
    "IRF7",
    "IRF9",
    "ISG15",
    "MX1",
    "OAS1",
    "OAS2",
    "IFIT1",
    "IFIT3",
    "CXCL10",
    "RSAD2",
    "IFI44",
    "BST2",
]

# Asked gene, plus two claudins reported only as epithelial context.
CLDN4 = ["CLDN4"]
CLDN_CONTEXT = ["CLDN3", "CLDN7"]

# Knockdown / payload checks. Not part of the asked panel.
QC = ["TACSTD2", "CDKN1A"]

PANEL = CLDN4 + CLDN_CONTEXT + NHEJ + STING + IFN + QC

SETS = {
    "CLDN4": CLDN4,
    "NHEJ": NHEJ,
    "STING": STING,
    "IFN": IFN,
}

# Mouse symbol used in GSE334497. OAS1 has no single mouse gene named Oas1.
MOUSE_SYMBOL = {
    "CLDN4": "Cldn4",
    "CLDN3": "Cldn3",
    "CLDN7": "Cldn7",
    "XRCC6": "Xrcc6",
    "XRCC5": "Xrcc5",
    "PRKDC": "Prkdc",
    "LIG4": "Lig4",
    "XRCC4": "Xrcc4",
    "NHEJ1": "Nhej1",
    "DCLRE1C": "Dclre1c",
    "PAXX": "Paxx",
    "CGAS": "Cgas",
    "STING1": "Sting1",
    "TBK1": "Tbk1",
    "IRF3": "Irf3",
    "STAT1": "Stat1",
    "STAT2": "Stat2",
    "IRF7": "Irf7",
    "IRF9": "Irf9",
    "ISG15": "Isg15",
    "MX1": "Mx1",
    "OAS1": "Oas1a",
    "OAS2": "Oas2",
    "IFIT1": "Ifit1",
    "IFIT3": "Ifit3",
    "CXCL10": "Cxcl10",
    "RSAD2": "Rsad2",
    "IFI44": "Ifi44",
    "BST2": "Bst2",
    "TACSTD2": "Tacstd2",
    "CDKN1A": "Cdkn1a",
}

MOUSE_ENSEMBL = {
    "Cldn4": "ENSMUSG00000047501",
    "Cldn3": "ENSMUSG00000070473",
    "Cldn7": "ENSMUSG00000018569",
    "Xrcc6": "ENSMUSG00000022471",
    "Xrcc5": "ENSMUSG00000026187",
    "Prkdc": "ENSMUSG00000022672",
    "Lig4": "ENSMUSG00000049717",
    "Xrcc4": "ENSMUSG00000021615",
    "Nhej1": "ENSMUSG00000026162",
    "Dclre1c": "ENSMUSG00000026648",
    "Paxx": "ENSMUSG00000047617",
    "Cgas": "ENSMUSG00000032344",
    "Sting1": "ENSMUSG00000024349",
    "Tbk1": "ENSMUSG00000020115",
    "Irf3": "ENSMUSG00000003184",
    "Stat1": "ENSMUSG00000026104",
    "Stat2": "ENSMUSG00000040033",
    "Irf7": "ENSMUSG00000025498",
    "Irf9": "ENSMUSG00000002325",
    "Isg15": "ENSMUSG00000035692",
    "Mx1": "ENSMUSG00000000386",
    "Oas1a": "ENSMUSG00000052776",
    "Oas2": "ENSMUSG00000032690",
    "Ifit1": "ENSMUSG00000034459",
    "Ifit3": "ENSMUSG00000074896",
    "Cxcl10": "ENSMUSG00000034855",
    "Rsad2": "ENSMUSG00000020641",
    "Ifi44": "ENSMUSG00000028037",
    "Bst2": "ENSMUSG00000046718",
    "Tacstd2": "ENSMUSG00000051397",
    "Cdkn1a": "ENSMUSG00000023067",
}

# Older symbols sometimes used as the gene_name column.
ALIASES = {
    "TMEM173": "STING1",
    "MB21D1": "CGAS",
    "C9ORF142": "PAXX",
    "STING": "STING1",
}
