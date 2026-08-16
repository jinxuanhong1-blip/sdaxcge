#!/usr/bin/env python3
"""Write an honest leftover catalog for all 53 leftover 2019-2021 lung ICI series."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "noskip" / "GEO_2019_2021"
FABLE = ROOT / "results" / "fable_geo_2019_2021"

inv = json.loads((OUT / "leftover_suppl_inventory.json").read_text())
triage = {r["accession"]: r for r in csv.DictReader(open(FABLE / "tables" / "triage_all_candidates.csv"))}

# Status assigned after download + gene scan + compute. No invented IDs.
STATUS = {
    "GSE190266": "COMPUTED: CLDN4 vs 6-month PFS (n=70). TACSTD2 absent from TPM matrix.",
    "GSE190265": "COMPUTED: both genes vs histology only (no ICI endpoint deposited).",
    "GSE181820": "COMPUTED: both genes vs GEO groups A/B/C (not a deposited ICI response).",
    "GSE146100": "COMPUTED descriptive: scRNA per-nodule means, n=3/1 patient. No inferential test.",
    "GSE145896": "COMPUTED: genes in PD-1+ TIL CD28+/- subsets (not patient ICI endpoint).",
    "GSE131933": "COMPUTED descriptive: tumor-cell means at D0; T1D7 958MB downloaded, genes absent there. Ex-vivo, not clinical survival.",
    "GSE184053": "COMPUTED: TACSTD2 ~0 in Temra; no ICI endpoint. CLDN4 Ensembl not in those count files.",
    "GSE152590": "REVISITED blood: processed TPM downloaded; TACSTD2 and CLDN4 absent.",
    "GSE141479": "REVISITED blood: 57MB series_matrix downloaded; probe IDs, gene symbols absent; no response/PFS.",
    "GSE179994": "REVISITED large: 441MB T-cell RDS downloaded; TACSTD2 string present; no patient-level stats extracted; no response/PFS.",
    "GSE180347": "NO processed matrix (RAW.tar only). Metadata has PD-L1 group. Not computed.",
    "GSE173351": "NO processed matrix (RAW.tar only). Not computed.",
    "GSE176021": "Response status deposited but NO processed gene matrix. Not computed.",
    "GSE176022": "TCR-seq / RAW.tar only. Not a gene-expression leftover compute.",
    "GSE136961": "First-pass: Oncomine panel lacks TACSTD2/CLDN4.",
    "GSE154286": "Chemo (not ICI) + RAW.tar only. Not computed.",
    "GSE144945": "RAW.tar only. Not computed.",
    "GSE120028": "3'SAGE mixed tumors; RAW.tar only. Not computed.",
    "GSE120101": "TCR-seq. Not gene expression.",
    "GSE99995": "IFN-γ / PD-L1 levels, not ICI treatment; series_matrix has no TACSTD2/CLDN4 symbols.",
    "GSE117570": "RAW.tar only; treatment-naive early NSCLC. Not computed.",
    "GSE171650": "Cell-line DE xlsx; not a gene matrix with TACSTD2/CLDN4.",
    "GSE129381": "Mouse/human cell-line counts; TACSTD2 all-zero in human QuantSeq row; not patient ICI.",
    "GSE129968": "H358 cell-line PD-L1 high/low; not patient ICI.",
    "GSE179934": "A549 cell-line counts; not patient ICI.",
    "GSE118933": "IPF fibroblasts, not ICI.",
    "GSE173896": "COPD lung counts; TACSTD2/CLDN4 not in file.",
    "GSE162154": "COPD airway epithelium; not ICI.",
    "GSE168707": "miRNA cisplatin; not ICI / not these genes.",
    "GSE159785": "COVID GeoMx protein; not ICI.",
    "GSE159787": "COVID GeoMx RNA; TACSTD2 present; not ICI.",
    "GSE139327": "HIV BAL/PBMC DE tables; not ICI.",
    "GSE124885": "Pediatric autoimmune BAL; not ICI.",
    "GSE128822": "Treg subset counts; TACSTD2/CLDN4 absent.",
    "GSE115246": "Methylation anti-PD-1; not mRNA of target genes.",
    "GSE119144": "Methylation; not mRNA.",
    "GSE126043": "Methylation beta; not mRNA.",
    "GSE126045": "SuperSeries methylation; not mRNA.",
    "GSE138571": "CRISPR screen; not patient ICI expression.",
    "GSE142620": "In-vitro EMT; RAW.tar.",
    "GSE127465": "scRNA mtx (myeloid atlas); not ICI-outcome cohort.",
    "GSE154826": "CITE-seq batches; not a single processed gene matrix.",
    "GSE148944": "Oral tongue SCC; not lung ICI outcome.",
    "GSE150972": "n=1 neoantigen; RAW.tar.",
    "GSE127825": "RAW.tar only.",
    "GSE141383": "Glioblastoma; not lung.",
    "GSE167235": "RAW.tar only.",
    "GSE129380": "ChIP-seq; not mRNA.",
    "GSE133715": "Cell-line KEAP1; RAW.tar.",
    "GSE135976": "Cell-line; RAW.tar.",
    "GSE136932": "Cell-line pirfenidone; RAW.tar.",
    "GSE136933": "Cell-line; RAW.tar.",
    "GSE136934": "Cell-line; RAW.tar.",
    "GSE153051": "n=2 cell-line; RAW.tar.",
}

rows = []
for acc, d in sorted(inv.items()):
    rows.append({
        "accession": acc,
        "category": d.get("category"),
        "n_samples": d.get("n_samples"),
        "n_suppl": len(d.get("files") or []),
        "n_processed_candidates": len(d.get("processed_candidates") or []),
        "status": STATUS.get(acc, "Leftover inventoried; no additional compute (see inventory)."),
        "title": d.get("title"),
    })

with open(OUT / "leftover_catalog.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"Wrote {len(rows)} leftover catalog rows")
