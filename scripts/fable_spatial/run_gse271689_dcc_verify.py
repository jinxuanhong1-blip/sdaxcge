"""Verification of GSE271689: rebuild patient-level tumor-compartment (PanCK/CK)
expression from the raw GEO DCC files and compare with the published
patient-level matrix (source data Fig. 6b) used in the survival analysis.

Steps: parse each AOI .dcc (probe barcode -> deduplicated counts), map probes
to genes with the NanoString WTA v1.0 PKC, Q3-normalize each AOI, average CK
AOIs per patient (spotid), and correlate with the published matrix.
"""

import glob
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from common import ensure_results

DATA = "/workspace/data/fable_spatial/gse271689"
OUT = ensure_results()


def load_pkc():
    pkc = json.load(open(os.path.join(DATA, "Hs_R_NGS_WTA_v1.0.pkc")))
    probe2gene = {}
    for t in pkc["Targets"]:
        gene = t["DisplayName"]
        cc = t["CodeClass"]
        for p in t["Probes"]:
            probe2gene[p["RTS_ID"]] = (gene, cc)
    return probe2gene


def parse_dcc(path):
    counts = {}
    in_summary = False
    with gzip.open(path, "rt") as f:
        for ln in f:
            ln = ln.strip()
            if ln == "<Code_Summary>":
                in_summary = True
                continue
            if ln == "</Code_Summary>":
                break
            if in_summary and "," in ln:
                rts, c = ln.split(",")
                counts[rts] = int(c)
    return counts


def main():
    probe2gene = load_pkc()
    ann = pd.read_csv(os.path.join(DATA, "gsm_annotation.csv"), dtype=str)
    ann = ann[(ann.celltype == "CK") & (ann.spotid != "No Template Control")]
    files = {os.path.basename(p).split("_", 1)[1].replace(".dcc.gz", ""): p
             for p in glob.glob(os.path.join(DATA, "dcc", "*.dcc.gz"))}

    per_aoi = {}
    for _, row in ann.iterrows():
        path = files.get(row.dcc)
        if path is None:
            continue
        counts = parse_dcc(path)
        if sum(counts.values()) < 5e4:  # drop failed AOIs
            continue
        gene_counts = {}
        for rts, c in counts.items():
            hit = probe2gene.get(rts)
            if hit is None or not hit[1].startswith("Endogenous"):
                continue
            gene_counts[hit[0]] = gene_counts.get(hit[0], 0) + c
        s = pd.Series(gene_counts)
        q3 = np.quantile(s[s > 0], 0.75)
        per_aoi[(row.spotid, row.dcc)] = s / q3

    df = pd.DataFrame(per_aoi).T
    df.index.names = ["spotid", "aoi"]
    patient = df.groupby(level="spotid").mean()
    print(f"CK AOIs kept: {len(df)}; patients: {len(patient)}")

    pub = pd.read_excel(pd.ExcelFile(os.path.join(DATA, "moesm10.xlsx")),
                        "source_data_Figure_6b").set_index("Spot_ID")
    pub.index = pub.index.astype(str)
    shared_pt = patient.index.intersection(pub.index)
    gene_cols = [c for c in pub.columns[:10130] if c in patient.columns]
    print(f"shared patients: {len(shared_pt)}; shared genes: {len(gene_cols)}")

    # per-patient correlation across genes (log scale)
    rows = []
    for ptid in shared_pt:
        a = np.log1p(patient.loc[ptid, gene_cols].astype(float).to_numpy())
        b = np.log1p(pub.loc[ptid, gene_cols].astype(float).to_numpy())
        r, _ = stats.pearsonr(a, b)
        rows.append(dict(spotid=ptid, pearson_r_log=round(float(r), 4)))
    per_pt = pd.DataFrame(rows)

    # per-gene correlation across patients for the target genes
    tg_rows = []
    for gene in ["TACSTD2", "CLDN4", "CD8A", "CXCL9"]:
        a = patient.loc[shared_pt, gene].astype(float)
        b = pub.loc[shared_pt, gene].astype(float)
        r, p = stats.spearmanr(a, b)
        tg_rows.append(dict(gene=gene, spearman_r=round(float(r), 4), p=float(p),
                            n_patients=len(shared_pt)))
    tg = pd.DataFrame(tg_rows)

    per_pt.to_csv(os.path.join(OUT, "gse271689_dcc_vs_published_per_patient.csv"), index=False)
    tg.to_csv(os.path.join(OUT, "gse271689_dcc_vs_published_targets.csv"), index=False)
    print("\nper-patient log-expression correlation (DCC-derived vs published):")
    print(per_pt.describe().loc[["mean", "min", "max"]].to_string())
    print("\nper-gene across-patient correlation:")
    print(tg.to_string(index=False))


if __name__ == "__main__":
    main()
