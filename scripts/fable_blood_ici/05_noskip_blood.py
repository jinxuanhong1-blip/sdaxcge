#!/usr/bin/env python3
"""Noskip catalog + TACSTD2/CLDN4 analyses for human NSCLC blood/PBMC/plasma ICI series.

Does not drop a series if a gene is absent: the catalog records present/absent
per gene. Association tests run only where the gene exists AND a response label
exists. No fabricated p-values.

Reads already-downloaded processed files from $GEO_DIR and $GEO_DIR/noskip.
Writes results/noskip/blood/.
"""
from __future__ import annotations

import gzip
import json
import os
import tarfile
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

GEO = os.environ.get("GEO_DIR", "/tmp/geo")
NS = os.path.join(GEO, "noskip")
OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                   "results", "noskip", "blood"))

CLARIOM = {
    "TACSTD2": "TC0100014340.hg.1",
    "CLDN4": "TC0700007993.hg.1",
    "EPCAM": "TC0200007506.hg.1",
    "PTPRC": "TC0100011105.hg.1",
}


def mw(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    if np.allclose(a, a[0]) and np.allclose(b, b[0]) and np.isclose(a[0], b[0]):
        return np.nan, np.nan
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return float(u), float(p)


def parse_series_meta(path):
    meta = {}
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("!series_matrix_table_begin"):
                break
            if not ln.startswith("!Sample_"):
                continue
            tag = ln.split("\t", 1)[0]
            vals = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            if tag == "!Sample_title":
                meta["title"] = vals
            elif tag == "!Sample_geo_accession":
                meta["gsm"] = vals
            elif tag == "!Sample_characteristics_ch1":
                if vals and ":" in vals[0]:
                    key = vals[0].split(":", 1)[0].strip()
                    meta[key] = [v.split(":", 1)[-1].strip() if ":" in v else v
                                 for v in vals]
    return meta


def extract_series_rows(path, probe_ids):
    """Return {probe: np.array} plus gsm list from a GEO series matrix."""
    want = set(probe_ids) | {"ID_REF"}
    found = {}
    with gzip.open(path, "rt") as fh:
        in_table = False
        for ln in fh:
            if ln.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if ln.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            pid = ln.split("\t", 1)[0].strip('"')
            if pid in want:
                found[pid] = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
    gsm = found.pop("ID_REF")
    out = {p: np.asarray(found[p], float) for p in found}
    return gsm, out


def analyze_gse111414(rows):
    expr = pd.read_csv(os.path.join(NS, "GSE111414_gene_counts.csv.gz"))
    meta = parse_series_meta(os.path.join(NS, "GSE111414_sm.txt.gz"))
    # matrix columns are like 373_N1 matching Patient BS373
    sample_cols = [c for c in expr.columns
                   if c not in ("ENTREZID", "SYMBOL", "GENENAME", "LENGTH")]
    recs = []
    for title, resp, tp, sid in zip(meta["title"], meta["response"],
                                    meta["timepoint"], meta["subject id"]):
        # "Patient BS373 timepoint N1" -> 373_N1
        num = sid.replace("Patient BS", "")
        col = f"{num}_{tp}"
        recs.append({"col": col, "patient": sid, "response": resp,
                     "timepoint": tp, "title": title})
    sm = pd.DataFrame(recs).set_index("col")
    genes = {"TACSTD2": "TACSTD2", "CLDN4": "CLDN4", "EPCAM": "EPCAM"}
    det_rows = []
    for g in genes:
        sub = expr.loc[expr.SYMBOL == g, sample_cols].iloc[0].astype(float)
        det_rows.append({
            "series": "GSE111414", "gene": g,
            "n_samples": int(sub.notna().sum()),
            "n_detected_count_gt0": int((sub > 0).sum()),
            "pct_samples_detected": round(100 * float((sub > 0).mean()), 2),
            "median_raw_count": float(sub.median()),
            "mean_raw_count": float(sub.mean()),
            "max_raw_count": float(sub.max()),
        })
    # baseline N1, one value per patient
    n1 = sm[sm.timepoint == "N1"]
    assoc = []
    for g in ("TACSTD2", "CLDN4"):
        vals = expr.loc[expr.SYMBOL == g, n1.index].iloc[0].astype(float)
        r = vals[n1.response.str.contains("responder;") & ~n1.response.str.contains("non")]
        # response strings: "responder; partial response (PR)" vs "nonresponder; ..."
        is_r = n1.response.str.startswith("responder;")
        rv = vals[is_r]
        nv = vals[~is_r]
        u, p = mw(rv, nv)
        assoc.append({
            "series": "GSE111414", "gene": g, "contrast": "PR_vs_PD_baseline_N1",
            "test": "MannWhitneyU", "n_PR": int(is_r.sum()), "n_PD": int((~is_r).sum()),
            "PR_mean_count": round(float(rv.mean()), 4),
            "PD_mean_count": round(float(nv.mean()), 4),
            "PR_n_detected": int((rv > 0).sum()),
            "PD_n_detected": int((nv > 0).sum()),
            "U": u, "p_value": p,
            "note": "raw counts; 5 PR vs 5 PD patients at timepoint N1",
        })
    rows["detect"].extend(det_rows)
    rows["assoc"].extend(assoc)
    # per-sample table
    out = sm.reset_index()
    for g in ("TACSTD2", "CLDN4", "EPCAM"):
        out[g] = expr.loc[expr.SYMBOL == g, out["col"]].iloc[0].to_numpy(float)
    out.to_csv(os.path.join(OUT, "gse111414_sample_values.csv"), index=False)
    return {
        "tacstd2": "present", "cldn4": "present",
        "n_PR_N1": 5, "n_PD_N1": 5,
        "tacstd2_N1_p": assoc[0]["p_value"],
        "cldn4_N1_p": assoc[1]["p_value"],
    }


def analyze_gse202417(rows):
    gsm, mat = extract_series_rows(os.path.join(NS, "GSE202417_sm.txt.gz"),
                                   list(CLARIOM.values()))
    meta = parse_series_meta(os.path.join(NS, "GSE202417_sm.txt.gz"))
    df = pd.DataFrame({
        "gsm": gsm,
        "title": meta["title"],
        "phenotype": meta["phenotype"],
        "time": meta["time"],
        "treatment": meta["treatment"],
    })
    df["patient"] = df["title"].str.replace(r"-[12]$", "", regex=True)
    for gene, probe in CLARIOM.items():
        df[gene] = mat[probe]
    df.to_csv(os.path.join(OUT, "gse202417_sample_values.csv"), index=False)
    for gene in ("TACSTD2", "CLDN4", "EPCAM", "PTPRC"):
        rows["detect"].append({
            "series": "GSE202417", "gene": gene,
            "n_samples": len(df),
            "median_intensity": round(float(df[gene].median()), 4),
            "mean_intensity": round(float(df[gene].mean()), 4),
            "min_intensity": round(float(df[gene].min()), 4),
            "max_intensity": round(float(df[gene].max()), 4),
        })
    pre = df[df.time == "pre-treatment"]
    for gene in ("TACSTD2", "CLDN4"):
        r = pre.loc[pre.phenotype == "R", gene]
        n = pre.loc[pre.phenotype == "NR", gene]
        u, p = mw(r, n)
        rows["assoc"].append({
            "series": "GSE202417", "gene": gene,
            "contrast": "R_vs_NR_pre_treatment",
            "test": "MannWhitneyU",
            "n_R": int(len(r)), "n_NR": int(len(n)),
            "R_mean": round(float(r.mean()), 4),
            "NR_mean": round(float(n.mean()), 4),
            "U": u, "p_value": p,
            "note": "Clariom S intensity; one pre-treatment sample per patient",
        })
    # paired pre vs post
    pairs = []
    for pid, sub in df.groupby("patient"):
        if set(sub.time) >= {"pre-treatment", "post-treatment"}:
            pairs.append(pid)
    pre_v = df[(df.patient.isin(pairs)) & (df.time == "pre-treatment")].sort_values("patient")
    post_v = df[(df.patient.isin(pairs)) & (df.time == "post-treatment")].sort_values("patient")
    for gene in ("TACSTD2", "CLDN4"):
        w, p = wilcoxon(pre_v[gene].to_numpy(), post_v[gene].to_numpy())
        rows["assoc"].append({
            "series": "GSE202417", "gene": gene,
            "contrast": "pre_vs_post_paired",
            "test": "WilcoxonSignedRank",
            "n_pairs": int(len(pairs)),
            "pre_mean": round(float(pre_v[gene].mean()), 4),
            "post_mean": round(float(post_v[gene].mean()), 4),
            "stat": float(w), "p_value": float(p),
            "note": "14 paired patients (pre untreated vs post nivo+bezafibrate)",
        })
    return {
        "tacstd2": "present", "cldn4": "present",
        "n_R_pre": int((pre.phenotype == "R").sum()),
        "n_NR_pre": int((pre.phenotype == "NR").sum()),
    }


def analyze_gse141479(rows):
    gsm, mat = extract_series_rows(os.path.join(NS, "GSE141479_sm.txt.gz"),
                                   list(CLARIOM.values()))
    meta = parse_series_meta(os.path.join(NS, "GSE141479_sm.txt.gz"))
    df = pd.DataFrame({
        "gsm": gsm, "title": meta["title"],
        "individual": meta["individual"],
        "treatment": meta["treatment"],
    })
    for gene, probe in CLARIOM.items():
        df[gene] = mat[probe]
    df.to_csv(os.path.join(OUT, "gse141479_sample_values.csv"), index=False)
    for gene in ("TACSTD2", "CLDN4", "EPCAM", "PTPRC"):
        rows["detect"].append({
            "series": "GSE141479", "gene": gene,
            "n_samples": len(df),
            "median_intensity": round(float(df[gene].median()), 4),
            "mean_intensity": round(float(df[gene].mean()), 4),
            "min_intensity": round(float(df[gene].min()), 4),
            "max_intensity": round(float(df[gene].max()), 4),
        })
    # paired pre/post where both exist
    both = []
    for pid, sub in df.groupby("individual"):
        if set(sub.treatment) >= {"Pre-treatment", "Post-treatment (nivolumab)"}:
            both.append(pid)
    pre = df[(df.individual.isin(both)) & (df.treatment == "Pre-treatment")].sort_values("individual")
    post = df[(df.individual.isin(both)) & (df.treatment == "Post-treatment (nivolumab)")].sort_values("individual")
    for gene in ("TACSTD2", "CLDN4"):
        w, p = wilcoxon(pre[gene].to_numpy(), post[gene].to_numpy())
        rows["assoc"].append({
            "series": "GSE141479", "gene": gene,
            "contrast": "pre_vs_post_nivolumab_paired",
            "test": "WilcoxonSignedRank",
            "n_pairs": int(len(both)),
            "pre_mean": round(float(pre[gene].mean()), 4),
            "post_mean": round(float(post[gene].mean()), 4),
            "stat": float(w), "p_value": float(p),
            "note": "no public response/survival labels; paired time contrast only",
        })
    return {"tacstd2": "present", "cldn4": "present", "n_paired": len(both)}


def analyze_gse235048(rows):
    expr = pd.read_csv(os.path.join(NS, "GSE235048_24CFlow_TPM.txt.gz"), sep="\t",
                       index_col=0)
    genes_present = {g: (g in expr.index) for g in ("TACSTD2", "CLDN4", "EPCAM", "PTPRC", "CD3D")}
    rec = {"series": "GSE235048"}
    for g, present in genes_present.items():
        if not present:
            rows["detect"].append({
                "series": "GSE235048", "gene": g, "n_samples": expr.shape[1],
                "status": "absent_from_matrix",
            })
            continue
        v = expr.loc[g].astype(float)
        rows["detect"].append({
            "series": "GSE235048", "gene": g, "n_samples": int(v.size),
            "n_detected_tpm_gt0": int((v > 0).sum()),
            "pct_samples_detected": round(100 * float((v > 0).mean()), 2),
            "median_TPM": round(float(v.median()), 4),
            "mean_TPM": round(float(v.mean()), 4),
            "max_TPM": round(float(v.max()), 4),
        })
        rec[g] = v.to_dict()
    # no response labels — do not invent a contrast
    sample = pd.DataFrame({"sample": expr.columns})
    for g in ("TACSTD2", "PTPRC", "CD3D"):
        if g in expr.index:
            sample[g] = expr.loc[g].to_numpy(float)
    sample["CLDN4"] = "ABSENT"
    sample.to_csv(os.path.join(OUT, "gse235048_sample_values.csv"), index=False)
    return {"tacstd2": "present", "cldn4": "absent"}


def analyze_gse216297(rows):
    # 3805-gene TEP panel; both targets absent (Ensembl IDs not in index)
    import pyreadr
    rdata = os.path.join(NS, "GSE216297_TEP_Count_Matrix.RData")
    if not os.path.exists(rdata):
        raw = gzip.decompress(open(os.path.join(NS, "GSE216297_TEP_Count_Matrix.RData.gz"), "rb").read())
        open(rdata, "wb").write(raw)
    mat = pyreadr.read_r(rdata)["TEP_Count_Matrix"]
    idx = set(mat.index.astype(str))
    present = {
        "TACSTD2": "ENSG00000184292" in idx,
        "CLDN4": "ENSG00000189143" in idx,
        "EPCAM": "ENSG00000119888" in idx,
        "PTPRC": "ENSG00000081237" in idx,
    }
    for g, ok in present.items():
        rows["detect"].append({
            "series": "GSE216297", "gene": g,
            "n_genes_on_panel": int(mat.shape[0]),
            "n_samples": int(mat.shape[1]),
            "status": "present" if ok else "absent_from_panel",
            "ensembl": {"TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143",
                        "EPCAM": "ENSG00000119888", "PTPRC": "ENSG00000081237"}[g],
        })
    meta = parse_series_meta(os.path.join(NS, "GSE216297_sm.txt.gz"))
    return {
        "tacstd2": "absent", "cldn4": "absent",
        "n_genes": int(mat.shape[0]), "n_samples": int(mat.shape[1]),
        "n_responder": int(sum(x == "Responder" for x in meta["treatment"])),
        "n_nonresponder": int(sum(x == "nonResponder" for x in meta["treatment"])),
    }


def analyze_gse100860(rows):
    tar_path = os.path.join(NS, "GSE100860_RAW.tar")
    values = {}
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".txt.gz"):
                continue
            fh = tf.extractfile(m)
            raw = gzip.decompress(fh.read())
            rec = {}
            for ln in raw.decode().splitlines()[1:]:
                if not ln.strip():
                    continue
                gene, val = ln.split("\t", 1)
                if gene in ("TACSTD2", "CLDN4", "EPCAM", "PTPRC"):
                    rec[gene] = float(val)
            values[os.path.basename(m.name)] = rec
    df = pd.DataFrame.from_dict(values, orient="index")
    df.index.name = "file"
    df.to_csv(os.path.join(OUT, "gse100860_sample_values.csv"))
    for g in ("TACSTD2", "CLDN4", "EPCAM", "PTPRC"):
        v = df[g].astype(float)
        rows["detect"].append({
            "series": "GSE100860", "gene": g,
            "n_files": int(v.size),
            "n_detected_fpkm_gt0": int((v > 0).sum()),
            "median_FPKM": round(float(v.median()), 4),
            "mean_FPKM": round(float(v.mean()), 4),
            "max_FPKM": round(float(v.max()), 4),
        })
    return {"tacstd2": "present", "cldn4": "present", "n_files": int(len(df))}


def analyze_gse213902(rows):
    """Detection rate in two 10x GEX pools (no per-patient response in GEO)."""
    from scipy.io import mmread

    def one_pool(prefix):
        feat = pd.read_csv(os.path.join(NS, f"{prefix}_features.tsv.gz"),
                           sep="\t", header=None, names=["ensembl", "symbol", "type"])
        mtx = mmread(os.path.join(NS, f"{prefix}_matrix.mtx.gz")).tocsr()
        n_cells = mtx.shape[1]
        out = {}
        for g in ("TACSTD2", "CLDN4", "EPCAM", "PTPRC"):
            hits = feat.index[feat.symbol == g].tolist()
            if not hits:
                out[g] = None
                continue
            # 10x MTX is genes x cells; scipy mmread keeps that orientation
            counts = np.asarray(mtx[hits[0], :].todense()).ravel()
            out[g] = {
                "n_cells": int(n_cells),
                "n_detected": int((counts > 0).sum()),
                "pct_detected": round(100 * float((counts > 0).mean()), 4),
                "mean_umi": round(float(counts.mean()), 6),
            }
        return out

    p1 = one_pool("GSE213902_MQpool1")
    p2 = one_pool("GSE213902_MQpool2")
    for pool, res in (("MQpool1", p1), ("MQpool2", p2)):
        for g, d in res.items():
            rec = {"series": "GSE213902", "gene": g, "pool": pool}
            if d is None:
                rec["status"] = "absent_from_features"
            else:
                rec.update(d)
            rows["detect"].append(rec)
    json.dump({"MQpool1": p1, "MQpool2": p2},
              open(os.path.join(OUT, "gse213902_detectability.json"), "w"), indent=2)
    return {"tacstd2": "present", "cldn4": "present"}


def analyze_gse310370(rows):
    # miRNA array GPL21572 — mRNA genes are not on the platform
    with gzip.open(os.path.join(NS, "GSE310370_sm.txt.gz"), "rt") as fh:
        n_rows = 0
        in_table = False
        for ln in fh:
            if ln.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if ln.startswith("!series_matrix_table_end"):
                break
            if in_table:
                n_rows += 1
    for g in ("TACSTD2", "CLDN4"):
        rows["detect"].append({
            "series": "GSE310370", "gene": g,
            "status": "absent_miRNA_platform_GPL21572",
            "n_probes_on_array": n_rows - 1,
            "n_samples": 4,
            "note": "all 4 plasma samples labelled PD; no mRNA features",
        })
    return {"tacstd2": "absent", "cldn4": "absent"}


def write_catalog(extra):
    """Full blood/PBMC/plasma ICI catalog. Absent genes stay in the table."""
    # pull already-computed primary-slice numbers
    d285 = pd.read_csv(os.path.join(os.path.dirname(OUT), "..", "fable_blood_ici",
                                    "gse285888_gene_detectability.csv"))
    r285 = pd.read_csv(os.path.join(os.path.dirname(OUT), "..", "fable_blood_ici",
                                    "gse285888_response_stats.csv"))
    d305 = pd.read_csv(os.path.join(os.path.dirname(OUT), "..", "fable_blood_ici",
                                    "gse305086_target_detectability.csv"))

    def pct(df, gene):
        return float(df.loc[df.gene == gene, "pct_cells_detected"].iloc[0])

    def pval(df, gene, contrast, metric):
        sub = df[(df.gene == gene) & (df.contrast == contrast) & (df.metric == metric)]
        return float(sub.p_value.iloc[0]) if len(sub) else np.nan

    catalog = [
        {
            "accession": "GSE285888",
            "compartment": "PBMC (scRNA-seq, baseline)",
            "platform": "GPL16791 10x / NovaSeq",
            "n": "222144 cells / 33 pts",
            "TACSTD2": "present (0.0806% cells)",
            "CLDN4": "present (0.1819% cells)",
            "response_or_survival_labels": "yes: CR/DR/PD + irAE severity",
            "analysis": "patient-level MWU vs response/irAE (see results/fable_blood_ici/)",
            "outcome": (f"no association; TACSTD2 CR-vs-PD p="
                        f"{pval(r285,'TACSTD2','CR_vs_PD','detrate'):.2f}; "
                        f"CLDN4 p={pval(r285,'CLDN4','CR_vs_PD','cp10k'):.2f}"),
        },
        {
            "accession": "GSE305086",
            "compartment": "whole blood (PAXgene, array)",
            "platform": "GPL570 HG-U133 Plus 2.0",
            "n": "173 (76 BL + 76 FU + 21 ctrl)",
            "TACSTD2": "present (probes at 20th–42nd percentile = background)",
            "CLDN4": "present (probe-discordant: 26th vs 71st percentile)",
            "response_or_survival_labels": "no public per-patient response/PFS",
            "analysis": "detectability + baseline-vs-control + paired BL-vs-FU",
            "outcome": "near-background; no usable response contrast in GEO",
        },
        {
            "accession": "GSE111414",
            "compartment": "PBMC CD8+ T cells (bulk RNA-seq)",
            "platform": "GPL21290 Illumina HiSeq 3000",
            "n": "20 samples / 10 pts (N1+N2)",
            "TACSTD2": "present (near-zero counts)",
            "CLDN4": "present (near-zero counts)",
            "response_or_survival_labels": "yes: PR vs PD (no survival)",
            "analysis": "MWU on baseline N1 raw counts, 5 PR vs 5 PD",
            "outcome": extra["GSE111414"]["outcome"],
        },
        {
            "accession": "GSE202417",
            "compartment": "PBMC CD8+ T cells (Clariom S)",
            "platform": "GPL23126 Clariom S Human",
            "n": "28 (14 pts × pre/post)",
            "TACSTD2": "present (TC0100014340.hg.1)",
            "CLDN4": "present (TC0700007993.hg.1)",
            "response_or_survival_labels": "yes: R vs NR (no survival)",
            "analysis": "MWU pre-treatment R vs NR; paired Wilcoxon pre vs post",
            "outcome": extra["GSE202417"]["outcome"],
        },
        {
            "accession": "GSE141479",
            "compartment": "blood CD8+ (Clariom S)",
            "platform": "GPL23126 Clariom S Human",
            "n": "74 (41 pre / 33 post nivo)",
            "TACSTD2": "present (TC0100014340.hg.1)",
            "CLDN4": "present (TC0700007993.hg.1)",
            "response_or_survival_labels": "NO response/survival in GEO",
            "analysis": "detectability + paired pre-vs-post only (not dropped)",
            "outcome": extra["GSE141479"]["outcome"],
        },
        {
            "accession": "GSE235048",
            "compartment": "PBMC (bulk RNA-seq TPM)",
            "platform": "GPL11154 Illumina HiSeq 2000",
            "n": "15",
            "TACSTD2": "present (low TPM)",
            "CLDN4": "ABSENT from processed matrix (other CLDNs present)",
            "response_or_survival_labels": "NO response/survival (treatment only)",
            "analysis": "detectability; CLDN4 recorded as absent; no outcome test",
            "outcome": extra["GSE235048"]["outcome"],
        },
        {
            "accession": "GSE216297",
            "compartment": "blood platelets / thrombocytes (TEP RNA-seq)",
            "platform": "GPL20301 Illumina HiSeq 4000",
            "n": "286 (88 Responder / 198 nonResponder)",
            "TACSTD2": "ABSENT (ENSG00000184292 not in 3805-gene TEP panel)",
            "CLDN4": "ABSENT (ENSG00000189143 not in 3805-gene TEP panel)",
            "response_or_survival_labels": "yes: Responder/nonResponder (unused — genes absent)",
            "analysis": "panel membership only; no association computed",
            "outcome": "both genes absent from the processed count matrix",
        },
        {
            "accession": "GSE213902",
            "compartment": "PBMC sorted T cells (10x GEX, 2 hashed pools)",
            "platform": "GPL24676 10x / NovaSeq",
            "n": "2 GEX pools (3903 + 5013 cells); no per-pt GEO labels",
            "TACSTD2": "present on 10x features",
            "CLDN4": "present on 10x features",
            "response_or_survival_labels": "NO per-patient response in GEO (hashed timepoints)",
            "analysis": "per-cell detectability in both GEX pools",
            "outcome": extra["GSE213902"]["outcome"],
        },
        {
            "accession": "GSE100860",
            "compartment": "peripheral blood CD8 T cells (nivo-bound vs unbound)",
            "platform": "GPL16791 Illumina HiSeq 2500",
            "n": "14 FPKM files",
            "TACSTD2": "present (near-zero FPKM)",
            "CLDN4": "present (near-zero FPKM)",
            "response_or_survival_labels": "NO response/survival in GEO",
            "analysis": "detectability across processed FPKM files",
            "outcome": extra["GSE100860"]["outcome"],
        },
        {
            "accession": "GSE310370",
            "compartment": "plasma (OAA-captured miRNA array)",
            "platform": "GPL21572 Affymetrix miRNA 4.0",
            "n": "4 (all labelled PD)",
            "TACSTD2": "ABSENT (miRNA platform, no mRNA probes)",
            "CLDN4": "ABSENT (miRNA platform, no mRNA probes)",
            "response_or_survival_labels": "all PD — no contrast possible",
            "analysis": "catalogued; genes absent",
            "outcome": "both genes absent",
        },
        {
            "accession": "GSE207715",
            "compartment": "plasma / EV miRNA (nivolumab NSCLC)",
            "platform": "GPL21576 / GPL25134 miRNA arrays",
            "n": "282",
            "TACSTD2": "ABSENT (miRNA platform)",
            "CLDN4": "ABSENT (miRNA platform)",
            "response_or_survival_labels": "prognostic miRNA study (not mRNA)",
            "analysis": "catalogued; genes absent — series not dropped",
            "outcome": "both genes absent; processed matrix is miRNA, not mRNA",
        },
        {
            "accession": "GSE306542",
            "compartment": "circulating cells (scRNA, metastatic lung cancer ICI)",
            "platform": "GPL24676 10x",
            "n": "2 GEO samples",
            "TACSTD2": "on standard 10x feature list (not re-quantified; n=2)",
            "CLDN4": "on standard 10x feature list (not re-quantified; n=2)",
            "response_or_survival_labels": "title claims response predictors; no usable n",
            "analysis": "catalogued; n=2 precludes response stats",
            "outcome": "kept in catalog; too few samples for a contrast",
        },
        {
            "accession": "GSE247754",
            "compartment": "systemic / circulating CD8 T cells (scRNA)",
            "platform": "GPL24676 10x",
            "n": "2 GEO samples",
            "TACSTD2": "on standard 10x feature list",
            "CLDN4": "on standard 10x feature list",
            "response_or_survival_labels": "no ICI-response labels in GEO",
            "analysis": "catalogued; n=2",
            "outcome": "kept in catalog; no outcome test",
        },
        {
            "accession": "GSE266035",
            "compartment": "circulating TIL (scRNA)",
            "platform": "GPL24676 10x",
            "n": "9",
            "TACSTD2": "on standard 10x feature list",
            "CLDN4": "on standard 10x feature list",
            "response_or_survival_labels": "kinetics study; no binary response table in GEO",
            "analysis": "catalogued; processed data is MTX inside RAW.tar (~95 MB)",
            "outcome": "kept in catalog; no public response/survival table",
        },
        {
            "accession": "GSE315510",
            "compartment": "T cells (mitochondrial-potential / immunotherapy)",
            "platform": "GPL24676 10x",
            "n": "30",
            "TACSTD2": "on standard 10x feature list",
            "CLDN4": "on standard 10x feature list",
            "response_or_survival_labels": "not an NSCLC-response cohort in GEO metadata",
            "analysis": "catalogued as circulating/T-cell ICI-adjacent; not NSCLC-outcome",
            "outcome": "kept in catalog; no NSCLC ICI response labels",
        },
    ]
    cat = pd.DataFrame(catalog)
    cat.to_csv(os.path.join(OUT, "catalog_blood_ici_series.csv"), index=False)
    return cat


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = {"detect": [], "assoc": []}
    extra = {}

    print("GSE111414", flush=True)
    s = analyze_gse111414(rows)
    extra["GSE111414"] = {
        "outcome": (f"TACSTD2 N1 PR-vs-PD p={s['tacstd2_N1_p']}; "
                    f"CLDN4 p={s['cldn4_N1_p']} (n=5 vs 5, near-zero counts)")
    }

    print("GSE202417", flush=True)
    s = analyze_gse202417(rows)
    extra["GSE202417"] = {"outcome": "see gse202417 rows in association_tests.csv"}

    print("GSE141479", flush=True)
    s = analyze_gse141479(rows)
    extra["GSE141479"] = {"outcome": f"paired pre/post n={s['n_paired']}; no response labels"}

    print("GSE235048", flush=True)
    s = analyze_gse235048(rows)
    extra["GSE235048"] = {"outcome": "TACSTD2 low TPM; CLDN4 ABSENT; no response labels"}

    print("GSE216297", flush=True)
    s = analyze_gse216297(rows)
    extra["GSE216297"] = s

    print("GSE100860", flush=True)
    s = analyze_gse100860(rows)
    extra["GSE100860"] = {"outcome": f"both genes present, near-zero FPKM across {s['n_files']} files"}

    print("GSE213902", flush=True)
    s = analyze_gse213902(rows)
    extra["GSE213902"] = {"outcome": "see gse213902_detectability.json; no per-patient labels"}

    print("GSE310370", flush=True)
    analyze_gse310370(rows)

    det = pd.DataFrame(rows["detect"])
    assoc = pd.DataFrame(rows["assoc"])
    det.to_csv(os.path.join(OUT, "detectability.csv"), index=False)
    assoc.to_csv(os.path.join(OUT, "association_tests.csv"), index=False)

    # fill GSE202417 outcome from the table we just wrote
    a = assoc[assoc.series == "GSE202417"]
    bits = []
    for _, r in a.iterrows():
        bits.append(f"{r.gene} {r.contrast} p={r.p_value}")
    extra["GSE202417"]["outcome"] = "; ".join(bits)

    cat = write_catalog(extra)
    # rewrite catalog now that 202417 outcome is filled
    # (write_catalog already used extra — call again)
    cat = write_catalog(extra)

    summary = {
        "n_series_in_catalog": int(len(cat)),
        "series_with_both_genes_absent": ["GSE216297", "GSE310370", "GSE207715"],
        "series_with_cldn4_absent_tacstd2_present": ["GSE235048"],
        "series_tested_vs_response": ["GSE285888", "GSE111414", "GSE202417"],
        "series_with_response_labels_but_genes_absent": ["GSE216297"],
        "note": "No series was dropped because it is blood. Absent genes are recorded as absent.",
    }
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
    print("\n=== association tests ===")
    print(assoc.to_string(index=False))
    print("\n=== catalog ===")
    print(cat[["accession", "TACSTD2", "CLDN4", "response_or_survival_labels"]].to_string(index=False))


if __name__ == "__main__":
    main()
