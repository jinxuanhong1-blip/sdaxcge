"""Tabula Muris FACS Smart-seq2 lung (GSE109774): Tacstd2 / Cldn4 vs Elf3/Grhl1/Klf4/Tfap2a.

Epithelial compartment = cell_ontology_class in
  {epithelial cell of lung, ciliated columnar cell of tracheobronchial tree}.
n is small (~138); this is a tissue-matched mouse check, not a powered DE.

Outputs:
  tables/tabula_muris_lung_coexpr.tsv
  tables/tabula_muris_lung_celltypes.tsv
"""
import os
import tarfile
import numpy as np
import pandas as pd
from scipy import stats
import common as C

TAR = f"{C.DATA}/GSE109774_Lung.tar.gz"
ANN = f"{C.DATA}/tabula_muris_annotations_facs.csv"

MOUSE = {
    "Tacstd2": "TACSTD2", "Cldn1": "CLDN1", "Cldn4": "CLDN4", "Cldn7": "CLDN7",
    "F11r": "F11R", "Pard3": "PARD3",
    "Elf3": "ELF3", "Grhl1": "GRHL1", "Klf4": "KLF4", "Tfap2a": "TFAP2A",
    "Epcam": "EPCAM", "Krt8": "KRT8", "Krt18": "KRT18",
}
EPITHELIAL = {
    "epithelial cell of lung",
    "ciliated columnar cell of tracheobronchial tree",
}


def cell_id_from_member(name):
    # Lung/A1-D043522-3_39_F-1-1.csv -> A1.D043522.3_39_F.1.1
    base = os.path.basename(name).replace(".csv", "")
    return base.replace("-", ".")


def main():
    ann = pd.read_csv(ANN, low_memory=False)
    lung = ann[ann.tissue == "Lung"].copy()
    lung["cell_ontology_class"].value_counts().to_csv(
        f"{C.TABLES}/tabula_muris_lung_celltypes.tsv", sep="\t", header=["n"])
    epi = set(lung.loc[lung.cell_ontology_class.isin(EPITHELIAL), "cell"])
    print(f"lung annotated {len(lung)}; epithelial {len(epi)}")

    want = set(MOUSE)
    series = {}
    with tarfile.open(TAR, "r:gz") as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".csv")]
        for i, m in enumerate(members):
            cid = cell_id_from_member(m.name)
            if cid not in epi:
                continue
            f = tf.extractfile(m)
            df = pd.read_csv(f, header=None, names=["gene", "count"])
            hit = df[df.gene.isin(want)].set_index("gene")["count"]
            series[cid] = hit
            if (i + 1) % 400 == 0:
                print(f"  scanned {i+1}/{len(members)}, kept {len(series)}")
    mat = pd.DataFrame(series).fillna(0)
    print("epithelial matrix", mat.shape, "genes", list(mat.index))
    # log1p of raw Smart-seq2 counts (already relatively deep)
    ln = np.log1p(mat)

    trop = "Tacstd2"
    cldn4 = "Cldn4"
    rows = []
    for g in MOUSE:
        if g not in ln.index:
            rows.append({"mouse_gene": g, "human_ortholog": MOUSE[g], "present": False})
            continue
        r2 = p2 = r4 = p4 = np.nan
        if trop in ln.index and g != trop:
            r2, p2 = stats.spearmanr(ln.loc[g], ln.loc[trop])
        if cldn4 in ln.index and g != cldn4:
            r4, p4 = stats.spearmanr(ln.loc[g], ln.loc[cldn4])
        rows.append({
            "mouse_gene": g, "human_ortholog": MOUSE[g], "present": True,
            "n_cells": ln.shape[1],
            "frac_expr": float((mat.loc[g] > 0).mean()),
            "spearman_r_vs_Tacstd2": r2, "spearman_p_vs_Tacstd2": p2,
            "spearman_r_vs_Cldn4": r4, "spearman_p_vs_Cldn4": p4,
        })
    out = pd.DataFrame(rows)
    out.to_csv(f"{C.TABLES}/tabula_muris_lung_coexpr.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
