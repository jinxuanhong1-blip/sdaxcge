"""Mouse gene panels for the Kras-lung ICB Cldn4 sweep.

NHEJ core matches the c-NHEJ machinery used on GSE50927 (Ku70/80, DNA-PKcs,
Artemis, XLF, XRCC4, LIG4, Pol-lambda, Pol-mu, PAXX). Current symbols are
used (Paxx, Cgas, Sting1, H2ax). IFN/chemokine matches that same panel.
"""

from __future__ import annotations

NHEJ_CORE = [
    "Xrcc6",
    "Xrcc5",
    "Prkdc",
    "Dclre1c",
    "Nhej1",
    "Xrcc4",
    "Lig4",
    "Poll",
    "Polm",
    "Paxx",
]

NHEJ_EXTENDED = [
    "Trp53bp1",
    "Rnf8",
    "Rnf168",
    "Rif1",
    "Mad2l2",
    "Atm",
    "H2ax",
    "Mre11a",
    "Rad50",
    "Nbn",
    "Parp1",
    "Pnkp",
    "Aplf",
    "Aptx",
]

IFN_CHEMOKINE = [
    "Ccl5",
    "Cxcl9",
    "Cxcl10",
    "Cxcl11",
    "Ifnb1",
    "Ifng",
    "Isg15",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Rsad2",
    "Stat1",
    "Irf7",
    "Oasl2",
    "Mx1",
    "Gbp4",
]

IFN_CHEMOKINE_ONLY = ["Ccl5", "Cxcl9", "Cxcl10", "Cxcl11"]

IFN_ISG = [
    "Isg15",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Rsad2",
    "Stat1",
    "Irf7",
    "Oasl2",
    "Mx1",
    "Gbp4",
]

STING_PAIR = ["Cgas", "Sting1"]

# Ensembl mouse IDs (release queried 2026-09). Version suffix is stripped
# at load time. Old symbols that are absent from current Ensembl are omitted.
ENSEMBL = {
    "Cldn4": "ENSMUSG00000047501",
    "Tacstd2": "ENSMUSG00000051397",
    "Xrcc6": "ENSMUSG00000022471",
    "Xrcc5": "ENSMUSG00000026187",
    "Prkdc": "ENSMUSG00000022672",
    "Dclre1c": "ENSMUSG00000026648",
    "Nhej1": "ENSMUSG00000026162",
    "Xrcc4": "ENSMUSG00000021615",
    "Lig4": "ENSMUSG00000049717",
    "Poll": "ENSMUSG00000025218",
    "Polm": "ENSMUSG00000020474",
    "Paxx": "ENSMUSG00000047617",
    "Trp53bp1": "ENSMUSG00000043909",
    "Rnf8": "ENSMUSG00000090083",
    "Rnf168": "ENSMUSG00000014074",
    "Rif1": "ENSMUSG00000036202",
    "Mad2l2": "ENSMUSG00000029003",
    "Atm": "ENSMUSG00000034218",
    "H2ax": "ENSMUSG00000049932",
    "Mre11a": "ENSMUSG00000031928",
    "Rad50": "ENSMUSG00000020380",
    "Nbn": "ENSMUSG00000028224",
    "Parp1": "ENSMUSG00000026496",
    "Pnkp": "ENSMUSG00000002963",
    "Aplf": "ENSMUSG00000030051",
    "Aptx": "ENSMUSG00000028411",
    "Cgas": "ENSMUSG00000032344",
    "Sting1": "ENSMUSG00000024349",
    "Tbk1": "ENSMUSG00000020115",
    "Irf3": "ENSMUSG00000003184",
    "Ifnb1": "ENSMUSG00000048806",
    "Zbp1": "ENSMUSG00000027514",
    "Aim2": "ENSMUSG00000037860",
    "Trex1": "ENSMUSG00000049734",
    "Ccl5": "ENSMUSG00000035042",
    "Cxcl9": "ENSMUSG00000029417",
    "Cxcl10": "ENSMUSG00000034855",
    "Cxcl11": "ENSMUSG00000060183",
    "Ifng": "ENSMUSG00000055170",
    "Isg15": "ENSMUSG00000035692",
    "Ifit1": "ENSMUSG00000034459",
    "Ifit2": "ENSMUSG00000045932",
    "Ifit3": "ENSMUSG00000074896",
    "Rsad2": "ENSMUSG00000020641",
    "Stat1": "ENSMUSG00000026104",
    "Irf7": "ENSMUSG00000025498",
    "Oasl2": "ENSMUSG00000029561",
    "Mx1": "ENSMUSG00000000386",
    "Gbp4": "ENSMUSG00000079363",
}

MODULES = {
    "nhej_core": NHEJ_CORE,
    "nhej_extended": NHEJ_EXTENDED,
    "nhej_union": NHEJ_CORE + NHEJ_EXTENDED,
    "ifn_chemokine": IFN_CHEMOKINE,
    "ifn_chemokine_only": IFN_CHEMOKINE_ONLY,
    "ifn_isg": IFN_ISG,
    "sting_pair": STING_PAIR,
}

# Thesis sign for Spearman(Cldn4, module): +1 means high Cldn4 with high module.
MODULE_THESIS_SIGN = {
    "nhej_core": 1,
    "nhej_extended": 1,
    "nhej_union": 1,
    "ifn_chemokine": -1,
    "ifn_chemokine_only": -1,
    "ifn_isg": -1,
    "sting_pair": -1,
}

IFN_MODULES = {"ifn_chemokine", "ifn_chemokine_only", "ifn_isg", "sting_pair"}
NHEJ_MODULES = {"nhej_core", "nhej_extended", "nhej_union"}
