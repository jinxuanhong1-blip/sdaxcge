#!/usr/bin/env python3
"""Reproduce the CLDN4-loss / IFN-MHC-I public-data hunt.

Uses only Python's standard library. Downloads deposited processed data and
metadata from GEO, then writes compact evidence tables. It intentionally does
not manufacture p-values for datasets that lack replicate-level processed data.
"""

from __future__ import annotations

import csv
import gzip
import io
import math
import statistics
import urllib.request
from pathlib import Path


HERE = Path(__file__).resolve().parent
CACHE = HERE / "_cache"
CACHE.mkdir(exist_ok=True)

URLS = {
    "gse50927_baseline": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
    "gse50927_vili_low": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkoloGenes.csv.gz",
    "gse50927_vili_high": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkohiGenes.csv.gz",
    "gse207704": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz",
    "gse22493_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/matrix/GSE22493_series_matrix.txt.gz",
    "gse22493_soft": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/soft/GSE22493_family.soft.gz",
}

PRIORITY = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]
HUMAN_APM = [
    "HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5", "TAP1", "TAP2", "TAPBP",
    "PSMB8", "PSMB9", "PSMB10", "ERAP1", "ERAP2", "CALR", "CANX", "PDIA3",
]
MOUSE_MAP = {
    "IFI27": "Ifi27l2a",
    "HLA-A": "H2-K1",
    "HLA-B": "H2-D1",
    "HLA-C": "H2-Q7",
}


def fetch(key: str) -> Path:
    path = CACHE / URLS[key].rsplit("/", 1)[-1]
    if not path.exists():
        with urllib.request.urlopen(URLS[key], timeout=120) as response:
            path.write_bytes(response.read())
    return path


def write_tsv(name: str, fieldnames: list[str], rows: list[dict]) -> None:
    with (HERE / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value, digits=4):
    if value in (None, ""):
        return ""
    return f"{float(value):.{digits}g}"


def read_50927(key: str) -> dict[str, dict]:
    with gzip.open(fetch(key), "rt") as handle:
        return {row["Marker.Symbol"].upper(): row for row in csv.DictReader(handle)}


def gse50927_rows() -> tuple[list[dict], list[dict]]:
    contrasts = [
        ("GSE50927_baseline", "gse50927_baseline", "Cldn4 KO / WT; no VILI", "direct, least confounded"),
        ("GSE50927_VILI_low", "gse50927_vili_low", "Cldn4 KO-low / WT; 2h VILI", "direct; same nominal injury"),
        ("GSE50927_VILI_high", "gse50927_vili_high", "Cldn4 KO-high / WT; 2h VILI", "direct; confounded by greater injury"),
    ]
    priority_rows, apm_rows = [], []
    for contrast_id, key, label, caveat in contrasts:
        data = read_50927(key)
        for human_gene in PRIORITY:
            mouse_gene = MOUSE_MAP.get(human_gene, human_gene)
            row = data.get(mouse_gene.upper())
            priority_rows.append({
                "contrast": contrast_id,
                "comparison": label,
                "human_signature_gene": human_gene,
                "measured_gene": mouse_gene,
                "mapping_note": "one-to-many IFI27-family proxy" if human_gene == "IFI27" else ("mouse MHC-I proxy" if human_gene == "HLA-A" else "direct symbol"),
                "log2FC_KO_vs_control": fmt(row["logFC"]) if row else "",
                "PValue": fmt(row["PValue"]) if row else "",
                "FDR": fmt(row["FDR"]) if row else "",
                "direction": "UP" if row and float(row["logFC"]) > 0 else ("DOWN" if row else "NA"),
                "caveat": caveat,
            })
        for human_gene in HUMAN_APM:
            mouse_gene = MOUSE_MAP.get(human_gene, human_gene)
            row = data.get(mouse_gene.upper())
            apm_rows.append({
                "contrast": contrast_id,
                "human_APM_gene": human_gene,
                "measured_gene": mouse_gene,
                "log2FC_KO_vs_control": fmt(row["logFC"]) if row else "",
                "PValue": fmt(row["PValue"]) if row else "",
                "FDR": fmt(row["FDR"]) if row else "",
                "direction": "UP" if row and float(row["logFC"]) > 0 else ("DOWN" if row else "NA"),
            })
    return priority_rows, apm_rows


def gse207704_rows() -> tuple[list[dict], list[dict]]:
    by_gene: dict[str, list[dict]] = {}
    with gzip.open(fetch("gse207704"), "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_gene.setdefault(row["gene_short_name"].upper(), []).append(row)

    priority_rows, apm_rows = [], []
    for cell, ko_col, wt_col in [
        ("MCF7", "MCF7_CLDN4KO_FPKM (fpkm)", "MCF7_WT_FPKM (fpkm)"),
        ("T47D", "T47D_CLDN4KO_FPKM (fpkm)", "T47D_WT_FPKM (fpkm)"),
    ]:
        for gene in PRIORITY:
            records = by_gene.get(gene, [])
            # Median across duplicate transcript/gene records. A 0.1 pseudocount
            # prevents undefined fold changes at very low deposited FPKM.
            lfcs = [
                math.log2((float(r[ko_col]) + 0.1) / (float(r[wt_col]) + 0.1))
                for r in records
            ]
            priority_rows.append({
                "contrast": f"GSE207704_{cell}",
                "comparison": f"{cell} CLDN4 CRISPR KO / WT",
                "human_signature_gene": gene,
                "measured_gene": gene if records else "",
                "mapping_note": "direct symbol" if records else "not present in deposited table",
                "log2FC_KO_vs_control": fmt(statistics.median(lfcs)) if lfcs else "",
                "PValue": "",
                "FDR": "",
                "direction": "UP" if lfcs and statistics.median(lfcs) > 0 else ("DOWN" if lfcs else "NA"),
                "caveat": "deposited FPKM are collapsed across n=2; no valid inferential test",
            })
        for gene in HUMAN_APM:
            records = by_gene.get(gene, [])
            lfcs = [
                math.log2((float(r[ko_col]) + 0.1) / (float(r[wt_col]) + 0.1))
                for r in records
            ]
            apm_rows.append({
                "contrast": f"GSE207704_{cell}",
                "human_APM_gene": gene,
                "measured_gene": gene if records else "",
                "log2FC_KO_vs_control": fmt(statistics.median(lfcs)) if lfcs else "",
                "PValue": "",
                "FDR": "",
                "direction": "UP" if lfcs and statistics.median(lfcs) > 0 else ("DOWN" if lfcs else "NA"),
            })
    return priority_rows, apm_rows


def gse22493_annotations() -> dict[str, str]:
    annotations: dict[str, str] = {}
    in_platform = False
    with gzip.open(fetch("gse22493_soft"), "rt", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_platform, header = True, None
                continue
            if line.startswith("!platform_table_end"):
                in_platform = False
            if not in_platform:
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            # Symbol is usually column 3; old HLA probes only carry a
            # "HLA-A--description" string in column 2.
            symbol = fields[2].strip() if len(fields) > 2 else ""
            if not symbol and len(fields) > 1 and "--" in fields[1]:
                symbol = fields[1].split("--", 1)[0].strip()
            annotations[fields[0]] = symbol.upper()
    return annotations


def gse22493_rows() -> tuple[list[dict], list[dict], list[dict]]:
    annotations = gse22493_annotations()
    values: dict[str, list[float]] = {}
    probe_rows: list[dict] = []
    in_table = False
    with gzip.open(fetch("gse22493_matrix"), "rt") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if not in_table or line.startswith("!") or line.startswith('"ID_REF"'):
                continue
            fields = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            probe, raw = fields[0], fields[1:]
            gene = annotations.get(probe, "")
            if gene not in set(PRIORITY + HUMAN_APM + ["CLDN4"]):
                continue
            numeric = [float(x) for x in raw if x != ""]
            values.setdefault(gene, []).extend(numeric)
            probe_rows.append({
                "probe": probe,
                "gene": gene,
                "replicate_1_log_ratio": raw[0] if len(raw) > 0 else "",
                "replicate_2_log_ratio": raw[1] if len(raw) > 1 else "",
                "replicate_3_log_ratio": raw[2] if len(raw) > 2 else "",
                "direction_definition": "positive = CLDN4 KD / control; deposited normalized log-ratio",
            })

    def summarized(gene: str, family: str) -> dict:
        observed = values.get(gene, [])
        median = statistics.median(observed) if observed else None
        return {
            "contrast": "GSE22493_SKOV3",
            family: gene,
            "measured_gene": gene if observed else "",
            "log2FC_KO_vs_control": fmt(median) if observed else "",
            "PValue": "",
            "FDR": "",
            "direction": "UP" if observed and median > 0 else ("DOWN" if observed else "NA"),
        }

    priority_rows = []
    for gene in PRIORITY:
        row = summarized(gene, "human_signature_gene")
        row.update({
            "comparison": "SKOV-3 CLDN4 lentiviral siRNA / CLDN4-high control",
            "mapping_note": "median of all observed probe-replicate log-ratios" if row["measured_gene"] else "not present/annotated",
            "caveat": "n=3 two-color arrays; inconsistent probes/replicates; no deposited p-values",
        })
        priority_rows.append(row)
    apm_rows = [summarized(gene, "human_APM_gene") for gene in HUMAN_APM]
    return priority_rows, apm_rows, probe_rows


def contrast_summary(priority_rows: list[dict], apm_rows: list[dict]) -> list[dict]:
    contrast_ids = list(dict.fromkeys(row["contrast"] for row in priority_rows))
    out = []
    for contrast in contrast_ids:
        pri = [r for r in priority_rows if r["contrast"] == contrast]
        apm = [r for r in apm_rows if r["contrast"] == contrast]
        measured_pri = [r for r in pri if r["direction"] != "NA"]
        measured_apm = [r for r in apm if r["direction"] != "NA"]
        sig_pri = [r for r in measured_pri if r["FDR"] and float(r["FDR"]) < 0.05 and r["direction"] == "UP"]
        sig_apm = [r for r in measured_apm if r["FDR"] and float(r["FDR"]) < 0.05 and r["direction"] == "UP"]
        up_pri = [r for r in measured_pri if r["direction"] == "UP"]
        up_apm = [r for r in measured_apm if r["direction"] == "UP"]
        if contrast == "GSE50927_baseline":
            verdict = "TOP: partial IFN/APM match; not broad MHC-I"
        elif contrast == "GSE50927_VILI_high":
            verdict = "weak/context-confounded directional match"
        else:
            verdict = "does not reproduce priority IFN/MHC-I pattern"
        out.append({
            "contrast": contrast,
            "priority_measured": len(measured_pri),
            "priority_up": len(up_pri),
            "priority_up_FDR_lt_0.05": len(sig_pri),
            "APM_measured": len(measured_apm),
            "APM_up": len(up_apm),
            "APM_up_FDR_lt_0.05": len(sig_apm),
            "verdict": verdict,
        })
    return out


def main() -> None:
    p509, a509 = gse50927_rows()
    p207, a207 = gse207704_rows()
    p224, a224, probes = gse22493_rows()
    priority_rows = p509 + p207 + p224
    apm_rows = a509 + a207 + a224
    write_tsv(
        "priority_signature.tsv",
        ["contrast", "comparison", "human_signature_gene", "measured_gene", "mapping_note",
         "log2FC_KO_vs_control", "PValue", "FDR", "direction", "caveat"],
        priority_rows,
    )
    write_tsv(
        "apm_signature.tsv",
        ["contrast", "human_APM_gene", "measured_gene", "log2FC_KO_vs_control",
         "PValue", "FDR", "direction"],
        apm_rows,
    )
    write_tsv(
        "contrast_summary.tsv",
        ["contrast", "priority_measured", "priority_up", "priority_up_FDR_lt_0.05",
         "APM_measured", "APM_up", "APM_up_FDR_lt_0.05", "verdict"],
        contrast_summary(priority_rows, apm_rows),
    )
    write_tsv(
        "GSE22493_probe_values.tsv",
        ["probe", "gene", "replicate_1_log_ratio", "replicate_2_log_ratio",
         "replicate_3_log_ratio", "direction_definition"],
        probes,
    )


if __name__ == "__main__":
    main()
