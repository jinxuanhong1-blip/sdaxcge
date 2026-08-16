#!/usr/bin/env python3
"""B5_BLCA: harmonise three public urothelial ICI cohorts into one per-sample table.

Cohorts
  IMvigor210 (Mariathasan 2018, atezolizumab)  - counts -> TPM -> log2(TPM+1)
  BACI       (GSE176307, Robertson 2021)       - salmon TPM -> log2(TPM+1)
  Snyder     (Snyder 2017, via PredictIO)      - log2(TPM+0.001) as deposited

Also extracts the PredictIO-harmonised copy of IMvigor210, used solely for the
pipeline-concordance check (S7) and explicitly excluded from pooling.

Outputs (results/w200/B5_BLCA/):
  per_sample_expression.csv  one row per patient, all genes of interest
  cohort_flow.csv            n at every stage: available -> RNA -> evaluable -> analysed
  validation.json            checks of the mirrored IMvigor210 against published values

Nothing here references CLDN4 in a filter, and no response-dependent filtering is
applied anywhere.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import subprocess
import sys

PRIMARY_GENE = "CLDN4"
CLDN4_ENSEMBL = "ENSG00000189143"

CD8_TEFF = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21", "PRF1"]
HOUSEKEEPING = ["ACTB", "GAPDH", "TBP", "RPL13A", "PGK1"]
EXPLORATORY = ["CLDN1", "CLDN2", "CLDN3", "CLDN7", "CLDN18", "TACSTD2",
               "EPCAM", "CDH1", "OCLN", "TJP1", "CD274"]
PANEL = [PRIMARY_GENE] + CD8_TEFF + HOUSEKEEPING + EXPLORATORY

RESPONDER = {"CR", "PR"}
NONRESPONDER = {"SD", "PD"}


def classify(recist):
    """Map a RECIST best-overall-response string to the prespecified ORR endpoint."""
    if recist is None:
        return None
    r = str(recist).strip().upper().rstrip("*")
    if r in RESPONDER:
        return 1
    if r in NONRESPONDER:
        return 0
    return None  # NE / 'Only scanned at baseline' / blank -> non-evaluable


def zscore(vals):
    ok = [v for v in vals if v is not None and math.isfinite(v)]
    if len(ok) < 2:
        return [None] * len(vals)
    m = sum(ok) / len(ok)
    sd = math.sqrt(sum((v - m) ** 2 for v in ok) / (len(ok) - 1))
    if sd == 0:
        return [None] * len(vals)
    return [None if (v is None or not math.isfinite(v)) else (v - m) / sd for v in vals]


# --------------------------------------------------------------------------- #
# IMvigor210
# --------------------------------------------------------------------------- #
def load_imvigor210(interim_dir):
    """Read the tables produced by extract_imvigor210.R."""
    gene_path = os.path.join(interim_dir, "imvigor210_gene_tpm.csv")
    pheno_path = os.path.join(interim_dir, "imvigor210_pheno.csv")
    with open(gene_path) as fh:
        genes = {r["sample_id"]: r for r in csv.DictReader(fh)}
    with open(pheno_path) as fh:
        pheno = {r["sample_id"]: r for r in csv.DictReader(fh)}

    out = []
    for sid, ph in pheno.items():
        g = genes.get(sid, {})
        rec = ph.get("Best Confirmed Overall Response", "")
        row = {
            "cohort": "IMvigor210",
            "sample_id": sid,
            "patient_id": ph.get("ANONPT_ID", sid),
            "expression_scale": "log2(TPM+1)",
            "recist": rec,
            "responder": classify(rec),
            "treatment": "atezolizumab",
            "tissue": ph.get("Tissue", ""),
            "immune_phenotype": ph.get("Immune phenotype", ""),
            "ecog": ph.get("Baseline ECOG Score", ""),
            "received_platinum": ph.get("Received platinum", ""),
            "tmb_per_mb": ph.get("FMOne mutation burden per MB", ""),
            "tcga_subtype": ph.get("TCGA Subtype", ""),
        }
        for gene in PANEL:
            tpm = g.get(f"{gene}__tpm")
            row[gene] = math.log2(float(tpm) + 1.0) if tpm not in (None, "", "NA") else None
        nc = g.get("CLDN4__normcount")
        row["CLDN4_deseq_norm"] = math.log2(float(nc) + 1.0) if nc not in (None, "", "NA") else None
        out.append(row)
    return out


# --------------------------------------------------------------------------- #
# BACI / GSE176307
# --------------------------------------------------------------------------- #
def parse_series_matrix(path):
    """GEO characteristics rows are ragged: parse 'key: value' per sample column."""
    titles, accs = [], []
    chars = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!Sample_title"):
                titles = [f.strip('"') for f in line.split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [f.strip('"') for f in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([f.strip('"') for f in line.split("\t")[1:]])

    n = len(titles)
    samples = []
    for i in range(n):
        attrs = {}
        for row in chars:
            if i < len(row) and row[i] and ":" in row[i]:
                k, v = row[i].split(":", 1)
                attrs[k.strip()] = v.strip()
        title = titles[i]
        baci = title.replace("Patient sample", "").strip()
        samples.append({"baci_id": baci, "gsm": accs[i] if i < len(accs) else "", "attrs": attrs})
    return samples


def load_baci(raw_dir):
    sm = parse_series_matrix(os.path.join(raw_dir, "GSE176307", "GSE176307_series_matrix.txt.gz"))
    clin = {s["baci_id"]: s for s in sm}

    key_path = os.path.join(raw_dir, "GSE176307",
                            "GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz")
    baci_to_rs, rs_to_baci = {}, {}
    with gzip.open(key_path, "rt") as fh:
        for r in csv.DictReader(fh):
            b = (r.get("Sample ID") or "").strip()
            rs = (r.get("Omniseq_RS_ID (RNAseq)") or "").strip()
            if b and rs:
                baci_to_rs[b] = rs
                rs_to_baci[rs] = b

    expr_path = os.path.join(raw_dir, "GSE176307", "GSE176307_salmon_tpm_gene.matrix.tsv.gz")
    want = set(PANEL)
    gene_vals = {}
    with gzip.open(expr_path, "rt") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        cols = [c.strip() for c in hdr[1:]]
        for row in rdr:
            sym = row[0].strip()
            if sym in want:
                vals = {}
                for c, v in zip(cols, row[1:]):
                    try:
                        vals[c] = float(v)
                    except (TypeError, ValueError):
                        vals[c] = None
                # duplicate symbol rows: sum TPM
                if sym in gene_vals:
                    for c in cols:
                        a, b = gene_vals[sym].get(c), vals.get(c)
                        gene_vals[sym][c] = (a or 0.0) + (b or 0.0) if (a is not None or b is not None) else None
                else:
                    gene_vals[sym] = vals

    out = []
    for baci_id, rs in sorted(baci_to_rs.items()):
        c = clin.get(baci_id)
        attrs = c["attrs"] if c else {}
        rec = attrs.get("io.response", "")
        # A few patients have more than one RNA-seq run, recorded as a comma-joined
        # id. Prespecified rule: average technical replicates within donor, never
        # count them as independent n.
        run_ids = [x.strip() for x in str(rs).split(",") if x.strip()]
        row = {
            "cohort": "BACI",
            "sample_id": rs,
            "patient_id": baci_id,
            "expression_scale": "log2(TPM+1)",
            "recist": rec,
            "responder": classify(rec),
            "treatment": attrs.get("io.therapy", ""),
            "tissue": attrs.get("primary tumor location", ""),
            "immune_phenotype": "",
            "ecog": attrs.get("ecog", ""),
            "received_platinum": "",
            "tmb_per_mb": attrs.get("tmb", ""),
            "tcga_subtype": "",
        }
        row["n_rna_runs"] = len(run_ids)
        for gene in PANEL:
            vs = [gene_vals.get(gene, {}).get(rid) for rid in run_ids]
            vs = [v for v in vs if v is not None]
            row[gene] = math.log2(sum(vs) / len(vs) + 1.0) if vs else None
        row["CLDN4_deseq_norm"] = None
        out.append(row)
    return out


# --------------------------------------------------------------------------- #
# PredictIO (Snyder, and the harmonised IMvigor210 copy)
# --------------------------------------------------------------------------- #
def read_predictio_tsv_rows(path):
    """PredictIO TSVs carry an unnamed rowname column, so header[i] aligns to row[i+1]."""
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = [h.strip('"') for h in next(rdr)]
        for row in rdr:
            if not row:
                continue
            name = row[0].strip('"')
            vals = [v.strip('"') for v in row[1:]]
            yield name, dict(zip(hdr, vals))


def load_predictio(raw_dir, study, cohort_label):
    base = os.path.join(raw_dir, "predictio", f"ICB_{study}")
    ens_by_symbol = {}
    for _, rec in read_predictio_tsv_rows(f"{base}_expr_genes.tsv"):
        sym = rec.get("gene_name", "")
        gid = rec.get("gene_id", "")
        if sym in set(PANEL) and gid:
            ens_by_symbol.setdefault(sym, []).append(gid)

    with open(f"{base}_expr.tsv") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        samples = [h.strip('"') for h in next(rdr)]
        wanted_ids = {gid: sym for sym, ids in ens_by_symbol.items() for gid in ids}
        gene_vals = {}
        for row in rdr:
            gid = row[0].strip('"')
            if gid in wanted_ids:
                sym = wanted_ids[gid]
                vals = []
                for v in row[1:]:
                    try:
                        vals.append(float(v))
                    except (TypeError, ValueError):
                        vals.append(None)
                gene_vals.setdefault(sym, []).append(vals)

    meta = {}
    for name, rec in read_predictio_tsv_rows(f"{base}_metadata.tsv"):
        meta[name] = rec

    out = []
    for j, sid in enumerate(samples):
        rec = meta.get(sid, {})
        recist = rec.get("recist", "")
        row = {
            "cohort": cohort_label,
            "sample_id": sid,
            "patient_id": rec.get("patientid", sid),
            "expression_scale": "log2(TPM+0.001)",
            "recist": recist,
            "responder": classify(recist),
            "treatment": rec.get("treatment", ""),
            "tissue": rec.get("cancer_type", ""),
            "immune_phenotype": "",
            "ecog": "",
            "received_platinum": "",
            "tmb_per_mb": rec.get("TMB_perMb", ""),
            "tcga_subtype": "",
        }
        for gene in PANEL:
            mats = gene_vals.get(gene)
            if not mats:
                row[gene] = None
                continue
            # already log scale: average duplicate Ensembl rows rather than summing
            vs = [m[j] for m in mats if m[j] is not None]
            row[gene] = sum(vs) / len(vs) if vs else None
        row["CLDN4_deseq_norm"] = None
        out.append(row)
    return out


# --------------------------------------------------------------------------- #
def add_scores(rows):
    """Per-cohort z-scored control signatures."""
    by_cohort = {}
    for r in rows:
        by_cohort.setdefault(r["cohort"], []).append(r)
    for _, rs in by_cohort.items():
        for label, genes in (("cd8_teff_z", CD8_TEFF), ("housekeeping_z", HOUSEKEEPING)):
            per_gene_z = {}
            for g in genes:
                per_gene_z[g] = zscore([r.get(g) for r in rs])
            for i, r in enumerate(rs):
                zs = [per_gene_z[g][i] for g in genes if per_gene_z[g][i] is not None]
                r[label] = sum(zs) / len(zs) if zs else None
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--interim-dir", default="data/interim")
    ap.add_argument("--out-dir", default="results/w200/B5_BLCA")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.interim_dir, exist_ok=True)

    if not os.path.exists(os.path.join(args.interim_dir, "imvigor210_gene_tpm.csv")):
        print("[run ] extract_imvigor210.R")
        subprocess.run(["Rscript", os.path.join(os.path.dirname(__file__),
                                                "extract_imvigor210.R")], check=True)

    imv = load_imvigor210(args.interim_dir)
    baci = load_baci(args.raw_dir)
    snyder = load_predictio(args.raw_dir, "Snyder", "Snyder")
    imv_harm = load_predictio(args.raw_dir, "Mariathasan", "IMvigor210_PredictIO")

    print(f"[n] IMvigor210={len(imv)} BACI={len(baci)} Snyder={len(snyder)} "
          f"IMvigor210_PredictIO={len(imv_harm)}")

    rows = add_scores(imv + baci + snyder + imv_harm)

    fields = (["cohort", "sample_id", "patient_id", "expression_scale", "recist",
               "responder", "treatment", "tissue", "immune_phenotype", "ecog",
               "received_platinum", "tmb_per_mb", "tcga_subtype"]
              + PANEL + ["CLDN4_deseq_norm", "cd8_teff_z", "housekeeping_z", "n_rna_runs"])
    out_path = os.path.join(args.out_dir, "per_sample_expression.csv")
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})
    print(f"[write] {out_path}")

    # ---- flow table: n at each prespecified stage --------------------------- #
    flow_path = os.path.join(args.out_dir, "cohort_flow.csv")
    with open(flow_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cohort", "n_available", "n_cldn4_present", "n_response_recorded",
                    "n_evaluable_analysed", "n_responders", "n_nonresponders",
                    "orr_percent", "n_excluded_nonevaluable", "expression_scale",
                    "pooled_in_meta"])
        for cohort in ["IMvigor210", "BACI", "Snyder", "IMvigor210_PredictIO"]:
            rs = [r for r in rows if r["cohort"] == cohort]
            n_av = len(rs)
            n_g = sum(1 for r in rs if r.get(PRIMARY_GENE) is not None)
            n_rec = sum(1 for r in rs if str(r.get("recist") or "").strip() not in ("", "NA"))
            an = [r for r in rs if r["responder"] is not None and r.get(PRIMARY_GENE) is not None]
            nr = sum(1 for r in an if r["responder"] == 1)
            nn = len(an) - nr
            w.writerow([cohort, n_av, n_g, n_rec, len(an), nr, nn,
                        round(100.0 * nr / len(an), 2) if an else "",
                        n_g - len(an),
                        rs[0]["expression_scale"] if rs else "",
                        "no (duplicate of IMvigor210)" if cohort == "IMvigor210_PredictIO" else "yes"])
    print(f"[write] {flow_path}")

    # ---- validation of the mirrored IMvigor210 against published values ----- #
    def recist_tally(rs):
        """Count RECIST categories, collapsing every non-evaluable spelling into one
        bucket. The two sources label the same 50 non-evaluable IMvigor210 patients
        differently ('NE' in the original object, 'NA' in the PredictIO harmonisation),
        so comparing raw strings would report a difference that does not exist."""
        tally = {}
        for r in rs:
            k = str(r["recist"] or "").strip().upper()
            if k in ("", "NA", "NE", "N/A"):
                k = "NON_EVALUABLE"
            tally[k] = tally.get(k, 0) + 1
        return tally

    imv_recist = recist_tally(imv)
    harm_recist = recist_tally(imv_harm)
    imv_pheno = {}
    for r in imv:
        k = str(r["immune_phenotype"] or "NA").strip()
        imv_pheno[k] = imv_pheno.get(k, 0) + 1

    expected = {"CR": 25, "PR": 43, "SD": 63, "PD": 167, "NON_EVALUABLE": 50}
    expected_pheno = {"desert": 76, "excluded": 134, "inflamed": 74, "NA": 64}

    validation = {
        "imvigor210_source": "third-party GitHub mirror of the now-404 official IMvigor210CoreBiologies distribution",
        "n_samples_observed": len(imv),
        "n_samples_expected_published": 348,
        "n_samples_match": len(imv) == 348,
        "recist_observed": imv_recist,
        "recist_expected_published": expected,
        "recist_match": {k: imv_recist.get(k, 0) == v for k, v in expected.items()},
        "recist_all_match": all(imv_recist.get(k, 0) == v for k, v in expected.items()),
        "immune_phenotype_observed": imv_pheno,
        "immune_phenotype_expected_published": expected_pheno,
        "immune_phenotype_all_match": all(imv_pheno.get(k, 0) == v for k, v in expected_pheno.items()),
        "independent_cross_check": {
            "description": "RECIST distribution of the same trial as independently harmonised by bhklab/PredictIO from Zenodo 7058399",
            "recist_observed": harm_recist,
            "non_evaluable_label_note": ("the original object labels the 50 non-evaluable patients 'NE' "
                                         "and the PredictIO harmonisation labels them 'NA'; both are "
                                         "collapsed to NON_EVALUABLE before comparison"),
            "agrees_with_mirror": all(harm_recist.get(k, 0) == imv_recist.get(k, 0)
                                      for k in set(harm_recist) | set(imv_recist)),
        },
        "note": ("The mirror is accepted only because sample count, full RECIST distribution and "
                 "immune-phenotype distribution all reproduce the published IMvigor210 values, and "
                 "an independently harmonised copy of the same trial gives an identical RECIST "
                 "distribution."),
    }
    vpath = os.path.join(args.out_dir, "validation.json")
    with open(vpath, "w") as fh:
        json.dump(validation, fh, indent=2)
    print(f"[write] {vpath}")
    print(f"[check] IMvigor210 RECIST reproduces published: {validation['recist_all_match']}")
    print(f"[check] immune phenotype reproduces published: {validation['immune_phenotype_all_match']}")
    print(f"[check] PredictIO copy agrees: {validation['independent_cross_check']['agrees_with_mirror']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
