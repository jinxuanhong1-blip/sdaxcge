#!/usr/bin/env python3
"""Download compact GEO processed tables + series matrices for extra mouse-lung ICI RNA."""
from __future__ import annotations

import json
import time
import urllib.request
from pathlib import Path

DATA = Path("notes/a4_extra_mouse_lung_ici/data")
DATA.mkdir(parents=True, exist_ok=True)
MANIFEST = Path("notes/a4_extra_mouse_lung_ici/raw/download_manifest.json")

# Compact processed files only (< ~50 MB expected). No 10x mtx / rds / cloupe.
FILES = [
    # series matrices (metadata + sometimes expression)
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239485/matrix/GSE239485_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297630/matrix/GSE297630_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182228/matrix/GSE182228_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114601/matrix/GSE114601_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114300/matrix/GSE114300_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169194/matrix/GSE169194_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169196/matrix/GSE169196_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246922/matrix/GSE246922_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260596/matrix/GSE260596_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244452/matrix/GSE244452_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157880/matrix/GSE157880_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256071/matrix/GSE256071_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262305/matrix/GSE262305_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233862/matrix/GSE233862_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274960/matrix/GSE274960_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285342/matrix/GSE285342_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317011/matrix/GSE317011_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190264/matrix/GSE190264_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295685/matrix/GSE295685_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298051/matrix/GSE298051_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197260/matrix/GSE197260_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330733/matrix/GSE330733_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330941/matrix/GSE330941_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179500/matrix/GSE179500_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137396/matrix/GSE137396_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/matrix/GSE137244_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274351/matrix/GSE274351_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE338nnn/GSE338923/matrix/GSE338923_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193895/matrix/GSE193895_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236258/matrix/GSE236258_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE209nnn/GSE209766/matrix/GSE209766_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241978/matrix/GSE241978_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277610/matrix/GSE277610_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214613/matrix/GSE214613_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271713/matrix/GSE271713_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266364/matrix/GSE266364_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227534/matrix/GSE227534_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165517/matrix/GSE165517_series_matrix.txt.gz",
    # processed expression
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239485/suppl/GSE239485_Processed_data.xlsx",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297630/suppl/GSE297630_processed_data.xlsx",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330941/suppl/GSE330941_filtered_tablecounts_tpm.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190264/suppl/GSE190264_TPM_LLC1.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE209nnn/GSE209766/suppl/GSE209766_LLC_core_table_gene.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260596/suppl/GSE260596_AllSamples_Genes_ReadCounts.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317011/suppl/GSE317011_counts_anno.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179500/suppl/GSE179500_20210620_XTR_Bulk_RNAseq_DESeq2norm_Counts.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179500/suppl/GSE179500_20210211_XTR_Bulk_RNAseq_DESeq2_results_RestoredvNon-Restored.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137396/suppl/GSE137396_Normalized_logtransformed_medcentred_genetable_GEMMnodule.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.fpkm.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244452/suppl/GSE244452_KPvsKL_deg_all.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114601/suppl/GSE114601_counts.normalized.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169194/suppl/GSE169194_annotated_log2_norm_count.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246922/suppl/GSE246922_KP_RNA_counts_vst.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246922/suppl/GSE246922_LLC1_RNA_counts_vst.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295685/suppl/GSE295685_TPM.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256071/suppl/GSE256071_210602GerA_l2tpm.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262305/suppl/GSE262305_gene_expression.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233862/suppl/GSE233862_Raw_counts.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274960/suppl/GSE274960_gene_count_matrix.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285342/suppl/GSE285342_LLC1_cnt_DDX54.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157880/suppl/GSE157880_Bulk048.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197260/suppl/GSE197260_RNAseqTPM_MM_EGFR-TKI-CD8.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274351/suppl/GSE274351_expredata_TPM_gene.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE338nnn/GSE338923/suppl/GSE338923_Lacun3_STK11_RNAseq_voom_normalized_counts.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330733/suppl/GSE330733_gene.count.annotated.tsv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193895/suppl/GSE193895_Tumor_KrasModel_GEMMs_RNAseq_GenewiseCounts.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277610/suppl/GSE277610_Counts.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214613/suppl/GSE214613_readcount.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271713/suppl/GSE271713_counts.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266364/suppl/GSE266364_All.DEG.Expression.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241978/suppl/GSE241978_2020-07-21_Sherr_analysis_CMT_KO_vs_Cas9Ctrl.xlsx",
]


def download(url: str) -> dict:
    name = url.rstrip("/").split("/")[-1]
    dest = DATA / name
    rec = {"url": url, "file": name, "ok": False, "bytes": 0, "error": None}
    if dest.exists() and dest.stat().st_size > 200:
        rec["ok"] = True
        rec["bytes"] = dest.stat().st_size
        rec["cached"] = True
        return rec
    req = urllib.request.Request(url, headers={"User-Agent": "a4-extra-mouse-lung-ici/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        dest.write_bytes(data)
        rec["ok"] = True
        rec["bytes"] = len(data)
    except Exception as e:
        rec["error"] = str(e)
        dest.write_text(f"ERROR {e}\n{url}\n")
    return rec


def main() -> None:
    rows = []
    for url in FILES:
        rec = download(url)
        rows.append(rec)
        status = "OK" if rec["ok"] else "FAIL"
        print(f"{status}\t{rec['bytes']:>10}\t{rec['file']}\t{rec.get('error') or ''}")
        time.sleep(0.15)
    MANIFEST.write_text(json.dumps(rows, indent=2))
    ok = sum(1 for r in rows if r["ok"])
    print(f"done {ok}/{len(rows)}")


if __name__ == "__main__":
    main()
