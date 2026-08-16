"""Export compact theme gene sets (from Enrichr/MSigDB via gseapy) to JSON
so single-cell scripts score the same biology used in the TCGA GSEA.
Output: data/theme_genesets.json
"""
import json
import gseapy as gp
import common as C

WANTED = {
    "keratinization": ("GO_Biological_Process_2021",
                       "keratinocyte differentiation (GO:0030216)"),
    "skin_barrier_env": ("GO_Cellular_Component_2021",
                         "cornified envelope (GO:0001533)"),
    "skin_barrier_dev": ("GO_Biological_Process_2021",
                         "epidermis development (GO:0008544)"),
    "tight_junction_cc": ("GO_Cellular_Component_2021",
                          "tight junction (GO:0070160)"),
    "tight_junction_kegg": ("KEGG_2021_Human", "Tight junction"),
    "emt": ("MSigDB_Hallmark_2020", "Epithelial Mesenchymal Transition"),
}

LINEAGE = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"],
    "immune_pan": ["PTPRC"],
    "Tcell": ["CD3D", "CD3E", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "Bcell": ["MS4A1", "CD79A", "CD79B"],
    "plasma": ["MZB1", "IGHG1", "JCHAIN"],
    "myeloid": ["LYZ", "CD68", "CD14", "FCGR3A", "AIF1"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"],
    "mast": ["TPSAB1", "CPA3", "MS4A2"],
}


def main():
    libs = {}
    for _, (lib, _) in WANTED.items():
        if lib not in libs:
            libs[lib] = gp.get_library(name=lib, organism="Human")
    themes = {}
    for name, (lib, term) in WANTED.items():
        genes = libs[lib].get(term)
        if genes is None:
            raise SystemExit(f"term not found: {lib} :: {term}")
        themes[name] = sorted(set(genes))
        print(f"{name}: {len(themes[name])} genes")
    out = {
        "themes": themes,
        "lineage": LINEAGE,
        "user_tj": C.USER_TJ_GENES,
        "user_tfs": C.USER_TFS,
        "trop2": C.TROP2,
    }
    with open(f"{C.DATA}/theme_genesets.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote data/theme_genesets.json")


if __name__ == "__main__":
    main()
