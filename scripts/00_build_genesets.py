#!/usr/bin/env python3
"""Fetch canonical keratin / tight-junction gene sets from Enrichr and save a
local, version-pinned GMT so the downstream analysis is fully reproducible and
network-independent.

Provenance (Enrichr libraries):
  - GO_Biological_Process_2021
  - GO_Cellular_Component_2021
  - KEGG_2021_Human
  - Reactome_2022

The panel of interest for claims A8/A9 is CLDN1, CLDN4, CLDN7, F11R (JAM-A),
PARD3 -- tight-junction / apical-polarity genes.
"""
import os
import json
import gseapy

OUT_DIR = "data/genesets"
os.makedirs(OUT_DIR, exist_ok=True)

# (output_set_name, enrichr_library, exact_term, category)
WANTED = [
    ("GOCC_KERATIN_FILAMENT", "GO_Cellular_Component_2021",
     "keratin filament (GO:0045095)", "keratin"),
    ("REACTOME_KERATINIZATION", "Reactome_2022",
     "Keratinization R-HSA-6805567", "keratin"),
    ("GOCC_CORNIFIED_ENVELOPE", "GO_Cellular_Component_2021",
     "cornified envelope (GO:0001533)", "keratin"),
    ("GOBP_KERATINOCYTE_DIFFERENTIATION", "GO_Biological_Process_2021",
     "keratinocyte differentiation (GO:0030216)", "keratin"),
    ("KEGG_TIGHT_JUNCTION", "KEGG_2021_Human",
     "Tight junction", "tight_junction"),
    ("REACTOME_TIGHT_JUNCTION_INTERACTIONS", "Reactome_2022",
     "Tight Junction Interactions R-HSA-420029", "tight_junction"),
    ("GOCC_BICELLULAR_TIGHT_JUNCTION", "GO_Cellular_Component_2021",
     "bicellular tight junction (GO:0005923)", "tight_junction"),
    ("GOBP_TIGHT_JUNCTION_ORGANIZATION", "GO_Biological_Process_2021",
     "tight junction organization (GO:0120193)", "tight_junction"),
]

PANEL = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]


def main():
    libs = {}
    for _, lib, _, _ in WANTED:
        if lib not in libs:
            libs[lib] = gseapy.get_library(lib)

    gmt_lines = []
    provenance = []
    for name, lib, term, cat in WANTED:
        genes = sorted(set(libs[lib][term]))
        gmt_lines.append("\t".join([name, f"{lib}:{term}"] + genes))
        panel_hits = [g for g in PANEL if g in set(genes)]
        provenance.append({
            "set": name, "category": cat, "library": lib, "term": term,
            "n_genes": len(genes), "panel_members": panel_hits,
        })
        print(f"{name:38s} n={len(genes):4d}  panel_in_set={panel_hits}")

    gmt_path = os.path.join(OUT_DIR, "claim_A8A9_genesets.gmt")
    with open(gmt_path, "w") as fh:
        fh.write("\n".join(gmt_lines) + "\n")
    with open(os.path.join(OUT_DIR, "claim_A8A9_genesets.provenance.json"), "w") as fh:
        json.dump({"panel": PANEL, "sets": provenance}, fh, indent=2)
    print(f"\nWrote {gmt_path} ({len(gmt_lines)} sets)")

    ker = sorted({g for name, lib, term, cat in WANTED if cat == "keratin"
                  for g in libs[lib][term]})
    tj = sorted({g for name, lib, term, cat in WANTED if cat == "tight_junction"
                 for g in libs[lib][term]})
    print(f"Union KERATIN n={len(ker)}; panel in union = {[g for g in PANEL if g in set(ker)]}")
    print(f"Union TIGHT_JUNCTION n={len(tj)}; panel in union = {[g for g in PANEL if g in set(tj)]}")
    with open(os.path.join(OUT_DIR, "claim_A8A9_genesets_union.gmt"), "w") as fh:
        fh.write("\t".join(["KERATIN_UNION", "union_of_keratin_sets"] + ker) + "\n")
        fh.write("\t".join(["TIGHT_JUNCTION_UNION", "union_of_tight_junction_sets"] + tj) + "\n")


if __name__ == "__main__":
    main()
