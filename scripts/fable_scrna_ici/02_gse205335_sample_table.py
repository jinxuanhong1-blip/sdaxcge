#!/usr/bin/env python3
"""Parse GSE205335 family SOFT into a per-sample table (GSM, patient, tissue,
RECIST, platform) and map each GSM to the orig.ident used in the cell
identity file (e.g. title 'P1006 EBUS_06' + 3' platform -> 'EBUS-06-3P').
Output: results/fable_scrna_ici/gse205335_sample_table.tsv
"""
import gzip
import re
import pandas as pd

SOFT = "/tmp/data_scrna_ici/GSE205335_family.soft.gz"
CELLID = "/tmp/data_scrna_ici/GSE205335_Lung_IO_CellIdentity.txt.gz"
OUT = "/workspace/results/fable_scrna_ici/gse205335_sample_table.tsv"

samples = []
cur = None
with gzip.open(SOFT, "rt") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split("=")[1].strip()}
            samples.append(cur)
        elif cur is not None and line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[1].strip()
        elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
            val = line.split("=", 1)[1].strip()
            key, v = val.split(":", 1)
            cur[key.strip().lower().replace(" ", "_")] = v.strip()

df = pd.DataFrame(samples)

def orig_ident(row):
    token = row["title"].split()[1]          # e.g. EBUS_06, LUNG_N31
    base = token.replace("_", "-")
    if base.startswith("LUNG-"):             # LUNG_T31 -> LUNG-T31
        pass
    suffix = "3P" if "3'" in row["platform"] else "5P"
    return f"{base}-{suffix}"

df["orig.ident"] = df.apply(orig_ident, axis=1)

ci_idents = set(pd.read_csv(CELLID, sep="\t", usecols=["orig.ident"])["orig.ident"])
df["in_cell_identity"] = df["orig.ident"].isin(ci_idents)
unmatched_ci = ci_idents - set(df["orig.ident"])
print("GSMs not matched to cell-identity orig.ident:")
print(df.loc[~df["in_cell_identity"], ["gsm", "title", "orig.ident"]].to_string())
print("orig.idents in cell identity without GSM match:", sorted(unmatched_ci))

df.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT} ({len(df)} samples)")
print(df[["gsm", "title", "orig.ident", "patient", "tissue", "recist",
          "cancer_subtype", "in_cell_identity"]].to_string())
