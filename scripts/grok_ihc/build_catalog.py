#!/usr/bin/env python3
"""Build cleaned scoring tables and a source catalog from downloaded files."""
from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "results" / "grok_ihc" / "extracted"
INV = ROOT / "results" / "grok_ihc" / "inventory"
NOTES = ROOT / "notes" / "grok_ihc"


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fields = fieldnames or list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    EXT.mkdir(parents=True, exist_ok=True)
    INV.mkdir(parents=True, exist_ok=True)
    NOTES.mkdir(parents=True, exist_ok=True)

    # Dum 2022 Pathobiology Table 1 — lung/pleura/thymus + lung NET
    # Columns: n on TMA, n analyzable, % negative / weak / moderate / strong
    dum_lung = [
        {
            "source": "Dum 2022 Pathobiology 89:245-258",
            "doi": "10.1159/000522206",
            "pmcid": "PMC9393818",
            "marker": "TROP2",
            "method": "IHC TMA, MSVA-733R; categories negative/weak/moderate/strong",
            "tumor_category": "Adenocarcinoma of the lung",
            "n_on_tma": 196,
            "n_analyzable": 181,
            "pct_negative": 6.1,
            "pct_weak": 4.4,
            "pct_moderate": 33.7,
            "pct_strong": 55.8,
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        },
        {
            "source": "Dum 2022 Pathobiology 89:245-258",
            "doi": "10.1159/000522206",
            "pmcid": "PMC9393818",
            "marker": "TROP2",
            "method": "IHC TMA, MSVA-733R; categories negative/weak/moderate/strong",
            "tumor_category": "Squamous cell carcinoma of the lung",
            "n_on_tma": 80,
            "n_analyzable": 72,
            "pct_negative": 0.0,
            "pct_weak": 12.5,
            "pct_moderate": 6.9,
            "pct_strong": 80.6,
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        },
        {
            "source": "Dum 2022 Pathobiology 89:245-258",
            "doi": "10.1159/000522206",
            "pmcid": "PMC9393818",
            "marker": "TROP2",
            "method": "IHC TMA, MSVA-733R; categories negative/weak/moderate/strong",
            "tumor_category": "Small-cell carcinoma of the lung",
            "n_on_tma": 16,
            "n_analyzable": 11,
            "pct_negative": 27.3,
            "pct_weak": 63.6,
            "pct_moderate": 9.1,
            "pct_strong": 0.0,
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        },
        {
            "source": "Dum 2022 Pathobiology 89:245-258",
            "doi": "10.1159/000522206",
            "pmcid": "PMC9393818",
            "marker": "TROP2",
            "method": "IHC TMA, MSVA-733R; categories negative/weak/moderate/strong",
            "tumor_category": "Lung NET",
            "n_on_tma": 19,
            "n_analyzable": 18,
            "pct_negative": 88.9,
            "pct_weak": 5.6,
            "pct_moderate": 5.6,
            "pct_strong": 0.0,
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        },
        {
            "source": "Dum 2022 Pathobiology 89:245-258",
            "doi": "10.1159/000522206",
            "pmcid": "PMC9393818",
            "marker": "TROP2",
            "method": "IHC TMA, MSVA-733R; categories negative/weak/moderate/strong",
            "tumor_category": "Mesothelioma, epithelioid",
            "n_on_tma": 39,
            "n_analyzable": 31,
            "pct_negative": 83.9,
            "pct_weak": 12.9,
            "pct_moderate": 3.2,
            "pct_strong": 0.0,
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        },
        {
            "source": "Dum 2022 Pathobiology 89:245-258",
            "doi": "10.1159/000522206",
            "pmcid": "PMC9393818",
            "marker": "TROP2",
            "method": "IHC TMA, MSVA-733R; categories negative/weak/moderate/strong",
            "tumor_category": "Mesothelioma, other types",
            "n_on_tma": 76,
            "n_analyzable": 52,
            "pct_negative": 84.6,
            "pct_weak": 13.5,
            "pct_moderate": 1.9,
            "pct_strong": 0.0,
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        },
        {
            "source": "Dum 2022 Pathobiology 89:245-258",
            "doi": "10.1159/000522206",
            "pmcid": "PMC9393818",
            "marker": "TROP2",
            "method": "IHC TMA, MSVA-733R; categories negative/weak/moderate/strong",
            "tumor_category": "Thymoma",
            "n_on_tma": 29,
            "n_analyzable": 23,
            "pct_negative": 39.1,
            "pct_weak": 43.5,
            "pct_moderate": 17.4,
            "pct_strong": 0.0,
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        },
    ]
    write_csv(EXT / "dum2022_trop2_tma" / "dum2022_Table1_lung_pleura_thymus.csv", dum_lung)

    hpa_lung = [
        {
            "source": "Human Protein Atlas pathology.tsv (v23 archive)",
            "url": "https://v23.proteinatlas.org/download/pathology.tsv.zip",
            "gene": "TACSTD2",
            "cancer": "lung cancer",
            "high": 0,
            "medium": 2,
            "low": 3,
            "not_detected": 5,
            "n_annotated": 10,
            "score_type": "ordinal IHC (High/Medium/Low/Not detected)",
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        },
        {
            "source": "Human Protein Atlas pathology.tsv (v23 archive)",
            "url": "https://v23.proteinatlas.org/download/pathology.tsv.zip",
            "gene": "CLDN4",
            "cancer": "lung cancer",
            "high": 0,
            "medium": 9,
            "low": 2,
            "not_detected": 0,
            "n_annotated": 11,
            "score_type": "ordinal IHC (High/Medium/Low/Not detected)",
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        },
    ]
    write_csv(EXT / "hpa" / "HPA_lung_cancer_IHC_counts.csv", hpa_lung)

    hashimoto_summary = [
        {
            "source": "Hashimoto 2025 Sci Rep 15:35427",
            "doi": "10.1038/s41598-025-19362-3",
            "marker": "TROP2",
            "method": "IHC; intensity 0-3 x proportion 1-4; overexpression = product 12",
            "cohort": "advanced NSCLC, nivolumab + ipilimumab",
            "n": 110,
            "n_overexpressed_product12": 46,
            "pct_overexpressed": 41.8,
            "n_intensity3": 64,
            "pct_intensity3": 58.1,
            "n_proportion4": 68,
            "pct_proportion4": 61.8,
            "patient_level_ipd": "no (aggregated Table 2 + Suppl Table A1)",
            "ici_annotated": "yes (Nivo-Ipi)",
        }
    ]
    write_csv(EXT / "hashimoto2025_scirep" / "hashimoto2025_TROP2_IHC_summary.csv", hashimoto_summary)

    pak_summary = [
        {
            "source": "Pak 2012 World J Surg Oncol 10:53",
            "doi": "10.1186/1477-7819-10-53",
            "marker": "TROP2",
            "method": "IHC membrane; total score = proportion(0-4) x intensity(0-3); overexpression >4",
            "histology": "adenocarcinoma",
            "n": 100,
            "n_overexpressed": 23,
            "pct_overexpressed": 23.0,
            "ici_annotated": "no",
            "patient_level_ipd": "no",
        },
        {
            "source": "Pak 2012 World J Surg Oncol 10:53",
            "doi": "10.1186/1477-7819-10-53",
            "marker": "TROP2",
            "method": "IHC membrane; total score = proportion(0-4) x intensity(0-3); overexpression >4",
            "histology": "squamous cell carcinoma",
            "n": 64,
            "n_overexpressed": 41,
            "pct_overexpressed": 64.1,
            "ici_annotated": "no",
            "patient_level_ipd": "no",
        },
    ]
    write_csv(EXT / "pak2012_wjso" / "pak2012_TROP2_overexpression_summary.csv", pak_summary)

    kuo_summary = [
        {
            "source": "Kuo 2025 PLoS ONE 20:e0321555",
            "doi": "10.1371/journal.pone.0321555",
            "marker": "Trop-2",
            "method": "IHC SP295 RPA; membrane H-score 0-300",
            "sample_set": "1 adenocarcinoma",
            "n_ihc": 107,
            "median_hscore": 68,
            "hscore_range": "0-251",
            "pct_any_positive_hscore": 82,
            "ipd_status": "request-only (datarequest@gilead.com)",
            "ici_annotated": "no (PD-L1 22C3 TPS in set 1 only)",
        },
        {
            "source": "Kuo 2025 PLoS ONE 20:e0321555",
            "doi": "10.1371/journal.pone.0321555",
            "marker": "Trop-2",
            "method": "IHC SP295 RPA; membrane H-score 0-300",
            "sample_set": "2 ADC+SCC",
            "n_ihc": 158,
            "median_hscore": "ADC 150; SCC 170",
            "hscore_range": "0-298",
            "pct_any_positive_hscore": 90,
            "ipd_status": "request-only (datarequest@gilead.com)",
            "ici_annotated": "no",
        },
    ]
    write_csv(EXT / "kuo2025_plosone" / "kuo2025_Trop2_Hscore_summary.csv", kuo_summary)

    gyulai_summary = [
        {
            "source": "Gyulai 2023 Pathol Oncol Res 29:1611328",
            "doi": "10.3389/pore.2023.1611328",
            "marker": "CLDN4",
            "method": "IHC H-score 0-300",
            "cohort": "pulmonary ACC n=12 + MEC n=23",
            "finding": "grade 2 vs grade 1 median H-score 60 vs 1; no patient-level table released",
            "patient_level_ipd": "no",
            "ici_annotated": "no",
        }
    ]
    write_csv(EXT / "gyulai2023_cldn_rare_lung" / "gyulai2023_CLDN4_Hscore_summary.csv", gyulai_summary)

    catalog = [
        {
            "id": "bessede2024_mif",
            "citation": "Bessede et al. Clin Cancer Res 2024;30:779-85",
            "doi": "10.1158/1078-0432.CCR-23-2566",
            "marker": "TROP2",
            "assay": "5-plex mIHF (PanCK, TROP2, CD8, PD-L1, DAPI); clone EPR2043",
            "setting": "advanced NSCLC + ICI (BIP n=50 mIHF; OAK/POPLAR RNA)",
            "what_exists": "mIHF images and quantitative scores",
            "public_status": "NOT PUBLIC",
            "access_note": "Paper: immunofluorescence datasets are not publicly available (consent). Access requires request to A. Italiano + ethics approval. RNA: EGA EGAS00001005013 (controlled). Extra clinical: Vivli/Roche.",
            "downloaded": "no (not public)",
            "local_path": "results/grok_ihc/raw/bessede_ccr2024/",
        },
        {
            "id": "bessede2024_figshare_S1S4",
            "citation": "Bessede et al. AACR Figshare supplements S1-S4",
            "doi": "10.1158/1078-0432.25232284 (S1); .25232281 (S2); .25232278 (S3); .25232275 (S4)",
            "marker": "TROP2 (context only)",
            "assay": "clinicopathologic tables (not IHC/mIF scores)",
            "setting": "OAK/POPLAR atezolizumab vs docetaxel; BIP mIHF/plasma",
            "what_exists": "S1 characteristics; S2 multivariate TACSTD2 RNA; S3 BIP mIHF n=50 characteristics; S4 plasma proteomics characteristics",
            "public_status": "PUBLIC (Figshare)",
            "access_note": "Open DOCX on AACR Figshare. S3 is the mIHF cohort clinical table — no per-patient TROP2 mIF scores.",
            "downloaded": "yes",
            "local_path": "results/grok_ihc/raw/bessede_ccr2024/ and extracted/bessede_ccr2024/",
        },
        {
            "id": "zenodo_adc_trop2_hscores",
            "citation": "ADC target profiling in NSCLC (Zenodo 18543127 / 18494664)",
            "doi": "10.5281/zenodo.18543127",
            "marker": "TROP-2",
            "assay": "AI IHC H-score membrane+cytoplasm; TMA images",
            "setting": "resected NSCLC n~1142; OS; no ICI annotation advertised",
            "what_exists": "cohort_data.csv with AI H-scores + TMA spots",
            "public_status": "RESTRICTED",
            "access_note": "Zenodo record public, files require login/access request. Download returned login HTML.",
            "downloaded": "no (restricted)",
            "local_path": "results/grok_ihc/raw/zenodo_adc_trop2/",
        },
        {
            "id": "hashimoto2025_scirep",
            "citation": "Hashimoto et al. Sci Rep 2025;15:35427",
            "doi": "10.1038/s41598-025-19362-3",
            "marker": "TROP2",
            "assay": "IHC intensity 0-3 x proportion 1-4; OE = 12",
            "setting": "advanced NSCLC Nivo-Ipi n=110",
            "what_exists": "aggregated Table 1-3 + Suppl Table A1 (PD-L1/histology x TROP2)",
            "public_status": "PUBLIC",
            "access_note": "CC-BY article + DOCX supplement downloaded.",
            "downloaded": "yes",
            "local_path": "results/grok_ihc/raw/hashimoto2025_scirep/ and extracted/hashimoto2025_scirep/",
        },
        {
            "id": "dum2022_tma",
            "citation": "Dum et al. Pathobiology 2022;89:245-258",
            "doi": "10.1159/000522206",
            "marker": "TROP2",
            "assay": "IHC TMA 18,563 tumors; 4-tier score",
            "setting": "lung ADC/SCC/SCLC/NET + mesothelioma (no ICI)",
            "what_exists": "Table 1 category-level % negative/weak/moderate/strong; raw IPD upon request",
            "public_status": "PUBLIC (paper table); IPD request-only",
            "access_note": "PMC9393818 OA. Category table transcribed from paper. PMC OA tarball 404 at fetch.",
            "downloaded": "yes (table transcribed; full PDF blocked by PMC interstitial)",
            "local_path": "results/grok_ihc/extracted/dum2022_trop2_tma/",
        },
        {
            "id": "kuo2025_plosone",
            "citation": "Kuo et al. PLoS ONE 2025;20:e0321555",
            "doi": "10.1371/journal.pone.0321555",
            "marker": "Trop-2",
            "assay": "IHC SP295 RPA membrane H-score",
            "setting": "NSCLC sample set 1 n=107 IHC; set 2 n=158 IHC",
            "what_exists": "summary H-score stats + S1-S7 figure TIFFs; IPD via Gilead",
            "public_status": "PUBLIC figures/paper; IPD request-only",
            "access_note": "datarequest@gilead.com; no patient-level H-score table in supplements.",
            "downloaded": "yes (PDF + S1-S7 TIF)",
            "local_path": "results/grok_ihc/raw/kuo2025_plosone/",
        },
        {
            "id": "pak2012_wjso",
            "citation": "Pak et al. World J Surg Oncol 2012;10:53",
            "doi": "10.1186/1477-7819-10-53",
            "marker": "TROP2",
            "assay": "IHC total score 0-12; OE >4",
            "setting": "resected NSCLC AdC n=100, SCC n=64; no ICI",
            "what_exists": "aggregated Tables 1-3 in PDF",
            "public_status": "PUBLIC",
            "access_note": "CC article PDF downloaded.",
            "downloaded": "yes",
            "local_path": "results/grok_ihc/raw/pak2012_wjso/",
        },
        {
            "id": "omori2021_jcrco",
            "citation": "Omori et al. J Cancer Res Clin Oncol 2021",
            "doi": "10.1007/s00432-021-03784-3",
            "marker": "TROP2",
            "assay": "IHC score 0-3; high = 3",
            "setting": "lung cancer paired pre/post treatment; 5 ICI cases no TROP2 change",
            "what_exists": "Suppl DOCX tables (aggregated high vs low)",
            "public_status": "PUBLIC supplements",
            "access_note": "MOESM1 docx + MOESM2 pdf downloaded.",
            "downloaded": "yes",
            "local_path": "results/grok_ihc/raw/omori2021_trop2_lung/",
        },
        {
            "id": "inamura2017_oncotarget",
            "citation": "Inamura et al. Oncotarget 2017;8:28725-35",
            "doi": "10.18632/oncotarget.15647",
            "marker": "TROP2",
            "assay": "IHC membrane intensity 0-2; high = I1>=50% or I2>=10%",
            "setting": "ADC 270 / SqCC 201 / HGNET 115; no ICI",
            "what_exists": "aggregated high vs no/low tables in paper",
            "public_status": "PUBLIC paper (CC-BY); publisher 403 at fetch",
            "access_note": "High TROP2: ADC 64% (172/270), SqCC 75% (150/201), HGNET 18% (21/115). No IPD table found.",
            "downloaded": "partial (PMC interstitial HTML only)",
            "local_path": "results/grok_ihc/raw/inamura2017_oncotarget/",
        },
        {
            "id": "hpa_pathology",
            "citation": "Human Protein Atlas pathology.tsv v23",
            "doi": "https://www.proteinatlas.org/about/download",
            "marker": "TACSTD2; CLDN4",
            "assay": "IHC TMA ordinal High/Medium/Low/Not detected",
            "setting": "lung cancer n=10 TACSTD2; n=11 CLDN4; no ICI",
            "what_exists": "cancer-type counts (not H-score, not IPD)",
            "public_status": "PUBLIC",
            "access_note": "v23/v19 pathology.tsv.zip downloaded; current proteinatlas.org/download/pathology.tsv.zip 404.",
            "downloaded": "yes",
            "local_path": "results/grok_ihc/extracted/hpa/",
        },
        {
            "id": "gyulai2023_cldn",
            "citation": "Gyulai et al. Pathol Oncol Res 2023;29:1611328",
            "doi": "10.3389/pore.2023.1611328",
            "marker": "CLDN4 (also CLDN1/2/3/5/7/18)",
            "assay": "IHC H-score 0-300",
            "setting": "rare lung ACC+MEC n=35; no ICI",
            "what_exists": "aggregated H-score comparisons in PDF; no IPD table",
            "public_status": "PUBLIC",
            "access_note": "PDF + HTML Table 1 (clinical only) downloaded.",
            "downloaded": "yes",
            "local_path": "results/grok_ihc/raw/gyulai2023_cldn_rare_lung/",
        },
        {
            "id": "moldvay2016_cldn",
            "citation": "Moldvay / Kovacs et al. Pathol Oncol Res 2017 (online 2016)",
            "doi": "10.1007/s12253-016-0115-0",
            "marker": "CLDN4 (CLDN1-7 panel)",
            "assay": "IHC semiquantitative scores 0-5",
            "setting": "stage I NSCLC n=137 (ADC/L-ADC/SCC); no ICI",
            "what_exists": "aggregated score comparisons; no IPD table",
            "public_status": "PUBLIC",
            "access_note": "PDF downloaded. CLDN4 differs ADC vs L-ADC (p=0.001).",
            "downloaded": "yes",
            "local_path": "results/grok_ihc/raw/moldvay2016_cldn_nsclc/",
        },
        {
            "id": "jung2014_cldn4",
            "citation": "Jung et al. Korean J Thorac Cardiovasc Surg 2014;47:262-8",
            "doi": "10.5090/kjtcs.2014.47.3.262",
            "marker": "CLDN4",
            "assay": "IHC intensity 0-3 + proportion bins; positive = intensity>1 and >50% cells",
            "setting": "resected lung ADC; relapse association; no ICI",
            "what_exists": "aggregated 2x2 tables in paper",
            "public_status": "PUBLIC",
            "access_note": "PMC4157471 fetch returned unrelated interstitial/wrong page. No IPD table indexed on Figshare/Zenodo.",
            "downloaded": "no usable PDF",
            "local_path": "results/grok_ihc/raw/jung2014_cldn4_luad/",
        },
    ]
    write_csv(INV / "source_catalog.csv", catalog)
    (INV / "source_catalog.json").write_text(
        json.dumps(
            {"generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "sources": catalog},
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"catalog n={len(catalog)}")


if __name__ == "__main__":
    main()
