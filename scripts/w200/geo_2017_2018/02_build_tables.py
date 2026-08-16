#!/usr/bin/env python3
"""Build triage + evidence tables for the 2017-2018 GEO lung ICI leftover slice.

Reads NCBI esummary records and the deposited matrices already downloaded
under results/w200/GEO_2017_2018/downloads/. Does not invent accessions or
outcomes. Writes:

  tables/triage_all_candidates.csv
  tables/gene_presence.csv
  tables/GSE93157_clinical.csv
  tables/GSE93157_lung_recist_counts.csv
  tables/GSE110390_21genes.csv
  tables/analyzable_marker_outcome.csv
"""
from __future__ import annotations

import csv
import gzip
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "w200" / "GEO_2017_2018"
DL = OUT / "downloads"
TABLES = OUT / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

# Manual, evidence-based triage. Every accession is a real NCBI GSE from
# candidates_metadata.json. Verdicts come from series SOFT / deposited files.
TRIAGE = {
    "GSE93157": dict(
        lung="Y", ici_treated="Y", sample_type="bulk tumor FFPE (NSCLC/HNSCC/melanoma)",
        platform="NanoString PanCancer Immune 730 (GPL19965)",
        ici_outcome_on_geo="Y (best.resp RECIST + PFS)",
        tacstd2_cldn4="N (not on 730-gene immune panel)",
        verdict="UNUSABLE — only 2017-2018 lung ICI cohort with deposited outcomes; target genes not measured",
        category="lung-ICI-outcomes-no-markers",
    ),
    "GSE110390": dict(
        lung="Y (97 NSCLC + 62 UC)", ici_treated="Y (durvalumab, NCT01693562)",
        sample_type="tumor needle biopsy RNA-seq",
        platform="Illumina NextSeq 500; only 21 IFNG-axis genes deposited",
        ici_outcome_on_geo="N (no RECIST/PFS/OS on GEO samples)",
        tacstd2_cldn4="N (not among the 21 deposited genes)",
        verdict="UNUSABLE — genuine lung ICI RNA-seq, but GEO holds a 21-gene IFNG subset and no outcomes",
        category="lung-ICI-no-markers-no-outcomes",
    ),
    "GSE91061": dict(
        lung="N (all samples melanoma)", ici_treated="Y (anti-PD-1 ± anti-CTLA-4, BMS-038)",
        sample_type="bulk tumor RNA-seq",
        platform="Illumina; FPKM/raw deposited",
        ici_outcome_on_geo="Y (PD/SD/PRCR)",
        tacstd2_cldn4="not checked (out of scope: not lung)",
        verdict="EXCLUDE — melanoma ICB (Riaz 2017); hit only because the abstract mentions NSCLC",
        category="not-lung",
    ),
    "GSE100860": dict(
        lung="Y (NSCLC patients)", ici_treated="Y (nivolumab, 1 dose)",
        sample_type="peripheral-blood CD8 T cells (bound vs unbound)",
        platform="RNA-seq; series matrix empty (SRA + RAW.tar)",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="not interpretable (T-cell compartment)",
        verdict="EXCLUDE — T-cell binding study, no tumor epithelium, no ICI endpoint",
        category="tcell-no-outcome",
    ),
    "GSE99531": dict(
        lung="Y (NSCLC)", ici_treated="N on profiled samples (paper mentions a separate PD-1 cohort)",
        sample_type="intratumoral CD8 T cells (PD-1 high/int/neg) + healthy TEM",
        platform="RNA-seq; series matrix empty (SRA + RAW.tar)",
        ici_outcome_on_geo="N (no response labels on GSM records)",
        tacstd2_cldn4="not interpretable (T-cell compartment)",
        verdict="EXCLUDE — CD8 subsets, not tumor TACSTD2/CLDN4, no per-sample ICI outcome",
        category="tcell-no-outcome",
    ),
    "GSE90728": dict(
        lung="Y (NSCLC)", ici_treated="N",
        sample_type="CD8 T cells from tumor vs adjacent lung",
        platform="Smart-seq2 bulk of sorted CD8",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="not interpretable (T-cell compartment)",
        verdict="EXCLUDE — untreated CD8 TIL atlas, not ICI-treated tumor epithelium",
        category="tcell-no-ICI",
    ),
    "GSE90729": dict(
        lung="N (HNSCC)", ici_treated="N",
        sample_type="CD8 T cells from HNSCC",
        platform="Smart-seq2",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="n/a",
        verdict="EXCLUDE — HNSCC CD8 companion of GSE90728",
        category="not-lung",
    ),
    "GSE99254": dict(
        lung="Y (NSCLC)", ici_treated="N",
        sample_type="scRNA-seq of T cells (Guo / Zhang lab)",
        platform="Smart-seq2; EGA raw",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="not interpretable (T cells)",
        verdict="EXCLUDE — baseline T-cell landscape, not ICI outcomes, not epithelium",
        category="scrna-tcell",
    ),
    "GSE111360": dict(
        lung="partial (air-liquid / organoid TME; mixed cancers)", ici_treated="N",
        sample_type="scRNA + VDJ organoid cultures",
        platform="10x; RAW.tar",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="not an ICI-outcome test",
        verdict="EXCLUDE — in-vitro organoid TME model",
        category="in-vitro",
    ),
    "GSE124199": dict(
        lung="Y (LUAD)", ici_treated="N",
        sample_type="sorted dendritic-cell subsets",
        platform="RNA-seq of DCs",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="not interpretable (DC compartment)",
        verdict="EXCLUDE — DC Wnt1 study, not ICI-treated patients",
        category="not-ICI",
    ),
    "GSE101929": dict(
        lung="Y (NSCLC)", ici_treated="N",
        sample_type="paired tumor/normal Affymetrix",
        platform="GPL570",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="possible on GPL570 but no ICI endpoint",
        verdict="EXCLUDE — race-comparison microarray, not ICI",
        category="not-ICI",
    ),
    "GSE102286": dict(
        lung="Y (NSCLC)", ici_treated="N",
        sample_type="miRNA array",
        platform="GPL23871",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="N (miRNA platform)",
        verdict="EXCLUDE — miRNA, not ICI, cannot measure TACSTD2/CLDN4 mRNA",
        category="not-ICI",
    ),
    "GSE106420": dict(
        lung="mention only (LN/BM memory T cells vs lung)", ici_treated="N",
        sample_type="memory CD8 from LN/BM",
        platform="RNA-seq",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="n/a",
        verdict="EXCLUDE — healthy-tissue TRM, not lung cancer ICI",
        category="not-cancer-ICI",
    ),
    "GSE108819": dict(
        lung="mention in abstract", ici_treated="N",
        sample_type="in-vitro T cells + breast-tumor DC pulse",
        platform="array",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="n/a",
        verdict="EXCLUDE — breast-tumor DC / T-cell AICD, not lung ICI",
        category="in-vitro",
    ),
    "GSE109010": dict(
        lung="Y (H1299 cell line)", ici_treated="N",
        sample_type="H1299 SMARCA4 / OXPHOS perturbation",
        platform="Affymetrix",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="irrelevant (cell line, no ICI)",
        verdict="EXCLUDE — cell-line metabolism, immunotherapy only in the paper intro",
        category="cell-line",
    ),
    "GSE109020": dict(
        lung="Y (H1299)", ici_treated="N",
        sample_type="SMARCA4 ChIP-seq",
        platform="ChIP-seq",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="N (occupancy, not expression vs ICI)",
        verdict="EXCLUDE — ChIP-seq companion of GSE109010",
        category="cell-line",
    ),
    "GSE113972": dict(
        lung="Y (3 lung biopsies + 4 B-ALL)", ici_treated="N",
        sample_type="primary tumor RNA-seq for TSA discovery",
        platform="HiSeq 2000",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="not an ICI-outcome test",
        verdict="EXCLUDE — antigen discovery, n=3 lung, no ICI treatment/outcome",
        category="not-ICI",
    ),
    "GSE115305": dict(
        lung="Y (NSCLC TIL)", ici_treated="N",
        sample_type="Siglec-9+ vs Siglec-9- CD8 TIL",
        platform="RNA-seq",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="not interpretable (CD8)",
        verdict="EXCLUDE — CD8 Siglec-9 sort, no ICI endpoint",
        category="tcell-no-ICI",
    ),
    "GSE120545": dict(
        lung="N", ici_treated="N",
        sample_type="endothelial cells + anti-HLA/A antibodies",
        platform="array",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="n/a",
        verdict="EXCLUDE — transplant endothelium; false-positive hit on 'PD'",
        category="not-cancer",
    ),
    "GSE121682": dict(
        lung="N (prostate lines)", ici_treated="N (olaparib ± NK/ADCC context)",
        sample_type="22Rv1 / DU145 + olaparib",
        platform="array",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="n/a",
        verdict="EXCLUDE — prostate cell lines",
        category="not-lung",
    ),
    "GSE54781": dict(
        lung="N (AML + healthy tissues)", ici_treated="N",
        sample_type="RT-PCR immunotherapy-target panel",
        platform="custom RT-PCR",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="not this question",
        verdict="EXCLUDE — AML SuperSeries; PDAT 2017 update of older data",
        category="not-lung",
    ),
    "GSE87879": dict(
        lung="Y (H1299)", ici_treated="N (CDK4/6 inhibitor PD 0332991)",
        sample_type="H1299 cell line",
        platform="array",
        ici_outcome_on_geo="N",
        tacstd2_cldn4="irrelevant",
        verdict="EXCLUDE — palbociclib in H1299; 'PD' is the drug name, not PD-1",
        category="cell-line",
    ),
}


def parse_series_chars(path: Path):
    chars = {}
    source = titles = gsms = None
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                gsms = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                source = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                key = vals[0].split(":", 1)[0].strip()
                chars[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    return gsms, titles, source, chars


def extract_gse93157():
    p = DL / "GSE93157" / "GSE93157_series_matrix.txt.gz"
    gsms, titles, source, chars = parse_series_chars(p)
    rows = []
    for i, gsm in enumerate(gsms):
        row = {
            "gsm": gsm,
            "title": titles[i],
            "source": source[i],
            "is_lung": "LUNG" in source[i].upper(),
        }
        for k, v in chars.items():
            row[k] = v[i]
        rows.append(row)
    with (TABLES / "GSE93157_clinical.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    lung = [r for r in rows if r["is_lung"]]
    counts = Counter(r["best.resp"] for r in lung)
    with (TABLES / "GSE93157_lung_recist_counts.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["subset", "n", "CR", "PR", "SD", "PD", "note"])
        w.writerow([
            "lung (non-squamous + squamous)",
            len(lung),
            counts.get("CR", 0),
            counts.get("PR", 0),
            counts.get("SD", 0),
            counts.get("PD", 0),
            "best.resp is the usable RECIST field; GEO 'response' is RC_RP_SD for all 35 lung samples including 14 PD and must not be used",
        ])
        w.writerow([
            "all histologies",
            len(rows),
            Counter(r["best.resp"] for r in rows).get("CR", 0),
            Counter(r["best.resp"] for r in rows).get("PR", 0),
            Counter(r["best.resp"] for r in rows).get("SD", 0),
            Counter(r["best.resp"] for r in rows).get("PD", 0),
            "25 melanoma + 5 HNSCC + 35 lung",
        ])
    return rows


def gene_presence():
    rows = []

    # GSE93157 NanoString
    p = DL / "GSE93157" / "GSE93157_raw_data_values.txt.gz"
    genes = []
    with gzip.open(p, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("ID_REF"):
                continue
            g = line.split("\t", 1)[0].strip().strip('"')
            if g and not g.startswith("#") and g != "ID_REF":
                genes.append(g)
    needles = ["TACSTD2", "TROP2", "CLDN4", "CLAUDIN4", "CLAUDIN-4"]
    hit = [g for g in genes if g.upper().replace("-", "").replace("_", "") in
           {n.replace("-", "") for n in needles} or g.upper() in needles]
    # also substring
    hit += [g for g in genes if any(n in g.upper() for n in ("TACSTD", "TROP2", "CLDN", "CLAUDIN"))]
    hit = sorted(set(hit))
    rows.append({
        "accession": "GSE93157",
        "file": p.name,
        "n_features": len(genes),
        "TACSTD2": "absent",
        "CLDN4": "absent",
        "other_hits": ";".join(hit) if hit else "none",
        "note": "NanoString PanCancer Immune Profiling Panel; epithelial TACSTD2/CLDN4 not on panel",
    })

    # GSE110390 21-gene
    p = DL / "GSE110390" / "GSE110390_CP1108_NSCLC_21gene_expression.txt.gz"
    with gzip.open(p, "rt") as f:
        header = f.readline()
        g21 = [line.split("\t", 1)[0].strip() for line in f if line.strip()]
    n_samples = len(header.strip().split("\t"))
    with (TABLES / "GSE110390_21genes.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "in_deposited_NSCLC_matrix"])
        for g in g21:
            w.writerow([g, "Y"])
        for g in ("TACSTD2", "CLDN4"):
            w.writerow([g, "N"])
    rows.append({
        "accession": "GSE110390",
        "file": p.name,
        "n_features": len(g21),
        "TACSTD2": "absent",
        "CLDN4": "absent",
        "other_hits": "none",
        "note": f"21 IFNG-axis genes x {n_samples} NSCLC samples; full transcriptome not deposited",
    })

    with (TABLES / "gene_presence.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows


def write_triage(meta):
    meta_by = {r["accession"]: r for r in meta}
    out_rows = []
    for acc, rec in sorted(meta_by.items()):
        t = TRIAGE.get(acc, {})
        out_rows.append({
            "accession": acc,
            "pdat": rec.get("pdat"),
            "n_samples": rec.get("n_samples"),
            "gpl": rec.get("gpl"),
            "pubmed": ",".join(rec.get("pubmedids") or []),
            "title": rec.get("title"),
            "lung": t.get("lung", ""),
            "ici_treated": t.get("ici_treated", ""),
            "sample_type": t.get("sample_type", ""),
            "platform": t.get("platform", ""),
            "ici_outcome_on_geo": t.get("ici_outcome_on_geo", ""),
            "tacstd2_cldn4": t.get("tacstd2_cldn4", ""),
            "category": t.get("category", "unreviewed"),
            "verdict": t.get("verdict", ""),
        })
    with (TABLES / "triage_all_candidates.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)

    with (TABLES / "analyzable_marker_outcome.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "accession", "n_lung", "outcome", "TACSTD2_testable", "CLDN4_testable",
            "test_run", "result",
        ])
        w.writerow([
            "NONE",
            0,
            "no 2017-2018 GEO series jointly measures TACSTD2/CLDN4 and a lung ICI endpoint",
            "N",
            "N",
            "N",
            "no association test was performed; a null p-value would be fabricated",
        ])
    return out_rows


def main():
    meta = json.loads((OUT / "candidates_metadata.json").read_text())
    extract_gse93157()
    gene_presence()
    rows = write_triage(meta)
    n_unusable_ici = sum(1 for r in rows if r["category"].startswith("lung-ICI"))
    print(f"triaged {len(rows)} series; lung-ICI leftover cohorts={n_unusable_ici}")
    print("wrote", TABLES)


if __name__ == "__main__":
    main()
