#!/usr/bin/env python3
"""Fetch public lung datasets into results/fable_stk11/data/processed/.

Sources (all public, no controlled access):
  * cBioPortal REST API
      - luad_tcga_pan_can_atlas_2018  (TCGA LUAD): mutation + RNA-seq  -> genotype<->expression
      - lusc_tcga_pan_can_atlas_2018  (TCGA LUSC): mutation + RNA-seq  -> genotype<->expression (secondary)
      - luad_mskcc_2015 (Rizvi 2015, Science): mutation + ICI response -> genotype<->response
      - nsclc_mskcc_2018 (Hellmann 2018, Cancer Cell): mutation + ICI response
      - nsclc_pd1_msk_2018 (Rizvi 2018, JCO): mutation + ICI response
  * NCBI GEO GSE135222 (Jung 2019, Nat Commun): anti-PD-1/PD-L1 NSCLC RNA-seq + PFS
      -> expression<->response

Only the handful of genes in fable_common.GENES are pulled from cBioPortal, so the
processed footprint stays tiny (well under 2 GB).
"""
from __future__ import annotations

import gzip
import os
import re

import pandas as pd

import fable_common as fc

CBIO_STUDIES = {
    "luad_tcga_pan_can_atlas_2018": {
        "label": "TCGA-LUAD",
        "kind": "tcga",
        "expr_profile": "luad_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "mut_profile": "luad_tcga_pan_can_atlas_2018_mutations",
    },
    "lusc_tcga_pan_can_atlas_2018": {
        "label": "TCGA-LUSC",
        "kind": "tcga",
        "expr_profile": "lusc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "mut_profile": "lusc_tcga_pan_can_atlas_2018_mutations",
    },
    "luad_mskcc_2015": {
        "label": "Rizvi2015",
        "kind": "ici",
        "expr_profile": None,
        "mut_profile": "luad_mskcc_2015_mutations",
    },
    "nsclc_mskcc_2018": {
        "label": "Hellmann2018",
        "kind": "ici",
        "expr_profile": None,
        "mut_profile": "nsclc_mskcc_2018_mutations",
    },
    "nsclc_pd1_msk_2018": {
        "label": "Rizvi2018",
        "kind": "ici",
        "expr_profile": None,
        "mut_profile": "nsclc_pd1_msk_2018_mutations",
    },
}

ALL_ENTREZ = sorted(set(fc.GENES.values()))
MUT_GENES = [fc.GENES[g] for g in (fc.GENOTYPE_GENES + ["TP53", "EGFR"])]


def sequenced_sample_list(study: str) -> str:
    """Return the sample-list id covering samples with mutation data."""
    lists = fc.cbio_get(f"/studies/{study}/sample-lists")
    for want in ("all_cases_with_mutation_data",):
        for sl in lists:
            if sl.get("category") == want:
                return sl["sampleListId"]
    # fallback
    return f"{study}_sequenced"


def fetch_clinical(study: str) -> pd.DataFrame:
    rows = []
    for dtype in ("PATIENT", "SAMPLE"):
        data = fc.cbio_get(f"/studies/{study}/clinical-data",
                           clinicalDataType=dtype, projection="DETAILED", pageSize=100000)
        for d in data:
            rows.append({
                "patientId": d.get("patientId"),
                "sampleId": d.get("sampleId", d.get("patientId")),
                "attr": d["clinicalAttributeId"],
                "value": d["value"],
                "level": dtype,
            })
    df = pd.DataFrame(rows)
    # sample-level attributes
    samp = df[df.level == "SAMPLE"].pivot_table(
        index="sampleId", columns="attr", values="value", aggfunc="first")
    # patient-level attributes: map to samples via a patient->sample table
    smap = (df[df.level == "SAMPLE"][["sampleId", "patientId"]]
            .drop_duplicates().set_index("sampleId")["patientId"])
    pat = df[df.level == "PATIENT"].pivot_table(
        index="patientId", columns="attr", values="value", aggfunc="first")
    out = samp.copy()
    out["patientId"] = smap
    pat_cols = [c for c in pat.columns if c not in out.columns]
    out = out.join(pat[pat_cols], on="patientId")
    out.index.name = "sampleId"
    return out.reset_index()


def fetch_mutations(study: str, profile: str, sample_list: str) -> pd.DataFrame:
    muts = fc.cbio_post(
        f"/molecular-profiles/{profile}/mutations/fetch?projection=DETAILED",
        {"entrezGeneIds": MUT_GENES, "sampleListId": sample_list})
    rows = [{
        "sampleId": m["sampleId"],
        "gene": fc.ENTREZ_TO_SYMBOL.get(m["entrezGeneId"], str(m["entrezGeneId"])),
        "proteinChange": m.get("proteinChange"),
        "mutationType": m.get("mutationType"),
    } for m in muts]
    return pd.DataFrame(rows)


def fetch_expression(study: str, profile: str, sample_list: str) -> pd.DataFrame:
    data = fc.cbio_post(
        f"/molecular-profiles/{profile}/molecular-data/fetch",
        {"entrezGeneIds": ALL_ENTREZ, "sampleListId": sample_list})
    rows = [{
        "sampleId": d["sampleId"],
        "gene": fc.ENTREZ_TO_SYMBOL.get(d["entrezGeneId"], str(d["entrezGeneId"])),
        "value": d["value"],
    } for d in data]
    long = pd.DataFrame(rows)
    if long.empty:
        return long
    wide = long.pivot_table(index="sampleId", columns="gene",
                            values="value", aggfunc="first")
    return wide.reset_index()


# nonsynonymous / functional mutation types kept as "mutant"
NONSYN = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site",
    "Nonstop_Mutation", "Splice_Region",
}


def build_genotype(mut: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    g = pd.DataFrame({"sampleId": samples})
    if not mut.empty:
        mut = mut[mut.mutationType.isin(NONSYN) | mut.mutationType.isna()]
    for gene in fc.GENOTYPE_GENES + ["TP53", "EGFR"]:
        hit = set(mut[mut.gene == gene]["sampleId"]) if not mut.empty else set()
        g[f"{gene}_mut"] = g.sampleId.isin(hit).astype(int)
    return g


def process_cbio() -> None:
    summary = {}
    for study, cfg in CBIO_STUDIES.items():
        print(f"\n=== {cfg['label']} ({study}) ===")
        seq_list = sequenced_sample_list(study)
        print(f"  sequenced sample list: {seq_list}")

        clin = fetch_clinical(study)
        clin.to_csv(os.path.join(fc.PROC, f"{cfg['label']}_clinical.csv"), index=False)

        mut = fetch_mutations(study, cfg["mut_profile"], seq_list)
        mut.to_csv(os.path.join(fc.PROC, f"{cfg['label']}_mutations.csv"), index=False)

        # sequenced samples define the WT/mutant universe
        seq_samples = sorted(set(
            fc.cbio_get(f"/sample-lists/{seq_list}").get("sampleIds", [])))
        if not seq_samples:  # fallback: samples that have any clinical row
            seq_samples = sorted(set(clin.sampleId))
        geno = build_genotype(mut, seq_samples)
        geno.to_csv(os.path.join(fc.PROC, f"{cfg['label']}_genotype.csv"), index=False)

        n_expr = 0
        if cfg["expr_profile"]:
            expr = fetch_expression(study, cfg["expr_profile"], f"{study}_all")
            expr.to_csv(os.path.join(fc.PROC, f"{cfg['label']}_expression.csv"), index=False)
            n_expr = 0 if expr.empty else expr.shape[0]

        summary[cfg["label"]] = {
            "study": study, "kind": cfg["kind"],
            "n_clinical_samples": int(clin.sampleId.nunique()),
            "n_sequenced_samples": len(seq_samples),
            "n_expression_samples": int(n_expr),
            "n_STK11_mut": int(geno.STK11_mut.sum()),
            "n_KEAP1_mut": int(geno.KEAP1_mut.sum()),
            "n_KRAS_mut": int(geno.KRAS_mut.sum()),
        }
        print("  ", summary[cfg["label"]])
    fc.save_json(summary, os.path.join(fc.PROC, "cbio_fetch_summary.json"))


# ---------------------------------------------------------------------------
# GEO GSE135222
# ---------------------------------------------------------------------------
def process_geo() -> None:
    print("\n=== GSE135222 (Jung 2019 anti-PD-1/PD-L1 NSCLC) ===")
    gz = os.path.join(fc.RAW, "GSE135222_exp.tsv.gz")
    fc.download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/"
        "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", gz)

    with gzip.open(gz, "rt") as f:
        expr = pd.read_csv(f, sep="\t")
    expr["ens"] = expr["gene_id"].str.split(".").str[0]
    keep = {v: k for k, v in fc.ENSEMBL.items()}
    sub = expr[expr.ens.isin(keep)].copy()
    sub["gene"] = sub.ens.map(keep)
    sub = sub.drop(columns=["gene_id", "ens"]).set_index("gene").T
    sub.index.name = "sample"
    sub = sub.reset_index()
    sub.to_csv(os.path.join(fc.PROC, "GSE135222_expression.csv"), index=False)

    # sample metadata (PFS event + time) from GEO brief soft record
    soft = os.path.join(fc.RAW, "GSE135222_family.txt")
    if not os.path.exists(soft):
        raw = fc._request(
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222"
            "&targ=gsm&form=text&view=brief")
        with open(soft, "wb") as fh:
            fh.write(raw)
    meta = parse_geo_soft(soft)
    meta.to_csv(os.path.join(fc.PROC, "GSE135222_clinical.csv"), index=False)
    print(f"   expression genes={list(sub.columns[1:])}, samples={sub.shape[0]}")
    print(f"   clinical rows={meta.shape[0]}, events={int(meta.pfs_event.sum())}")


def parse_geo_soft(path: str) -> pd.DataFrame:
    with open(path, "r", errors="ignore") as f:
        txt = f.read()
    blocks = re.split(r"\^SAMPLE = ", txt)
    rows = []
    for b in blocks[1:]:
        title = re.search(r"!Sample_title = (.+)", b)
        pfs = re.search(r"progression-free survival \(pfs\): (\d+)", b, re.I)
        ptime = re.search(r"pfs\.time: (\d+)", b, re.I)
        gender = re.search(r"gender: (\w+)", b, re.I)
        age = re.search(r"age: (\d+)", b, re.I)
        if not (title and pfs and ptime):
            continue
        sample = title.group(1).strip().replace(" ", "")
        rows.append({
            "sample": sample,
            "pfs_event": int(pfs.group(1)),        # 1 = progressed, 0 = censored
            "pfs_time_days": int(ptime.group(1)),
            "gender": gender.group(1) if gender else None,
            "age": int(age.group(1)) if age else None,
        })
    df = pd.DataFrame(rows)
    # durable clinical benefit proxy: PFS >= 6 months (>=180 days)
    df["DCB"] = (df.pfs_time_days >= 180).map({True: "YES", False: "NO"})
    return df


if __name__ == "__main__":
    process_cbio()
    process_geo()
    print(f"\nprocessed dir size = {fc.dir_size_mb(fc.PROC):.1f} MB")
    print(f"raw dir size       = {fc.dir_size_mb(fc.RAW):.1f} MB")
