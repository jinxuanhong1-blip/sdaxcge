#!/usr/bin/env python3
"""Inventory GSE76628 for a mouse-level Cldn4 vs T/NK test.

This series is a bulk Affymetrix Ad-VEGF-A164 flank angiogenesis model
(Uhlik et al., Cancer Res 2016). It is not lung Kras/Lkb1 scRNA. The
script writes honest n tables and two inventory figures. It does not
score Cldn4 %pos or T/NK fractions.
"""

from __future__ import annotations

import csv
import gzip
import json
import shutil
import urllib.request
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76628/matrix/"
    "GSE76628_series_matrix.txt.gz"
)
SISTER_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76588/matrix/"
    "GSE76588_series_matrix.txt.gz"
)
ANNOT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1261/annot/"
    "GPL1261.annot.gz"
)
RAW_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76628/suppl/"
    "GSE76628_RAW.tar"
)

PROBE_GENES = [
    "Cldn4",
    "Tacstd2",
    "Cd3d",
    "Cd3e",
    "Nkg7",
    "Klrd1",
    "Epcam",
    "Ncr1",
    "Cd8a",
]


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    tmp.rename(dest)


def parse_series_matrix(path: Path) -> dict[str, list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    meta: dict[str, list[str]] = {}
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                continue
            parts = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            meta[parts[0]] = parts[1:]
    return meta


def parse_probes(annot_path: Path, genes: list[str]) -> dict[str, list[str]]:
    hits: dict[str, list[str]] = {g: [] for g in genes}
    with gzip.open(annot_path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith(("#", "^", "!")):
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            probe = parts[0]
            for g in genes:
                for field in parts[1:8]:
                    toks = [t.strip() for t in field.replace("///", "//").split("//")]
                    if g in toks:
                        hits[g].append(probe)
                        break
    return hits


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def group_from_source(source: str) -> tuple[str, str, str]:
    """Return (day, treatment, group_label)."""
    if "control_day 0" in source:
        return "0", "none", "control_d0"
    day = "NA"
    for token in ("day 5", "day 20", "day 60"):
        if token in source:
            day = token.split()[-1]
    if "DC101" in source:
        treat = "DC101"
    elif "G6" in source:
        treat = "G6"
    elif "No Treatment" in source:
        treat = "NT"
    else:
        treat = "unknown"
    return day, treat, f"AdVEGF_d{day}_{treat}"


def fig_sample_n(groups: list[tuple[str, int]], out: Path) -> None:
    labels = [g for g, _ in groups]
    counts = [n for _, n in groups]
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    ax.bar(range(len(labels)), counts, color="#4C78A8", edgecolor="black", linewidth=0.4)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("n arrays (one per mouse)")
    ax.set_title("GSE76628 public processed matrix — n=78 bulk flank arrays")
    ax.set_ylim(0, max(counts) + 2)
    for i, n in enumerate(counts):
        ax.text(i, n + 0.15, str(n), ha="center", va="bottom", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def fig_nogo_gate(out: Path) -> None:
    rows = [
        ("Public processed matrix present?", "Yes — series matrix, 45,101 probes × 78 samples"),
        ("Cldn4 on platform?", "Yes — 1418283_at (not the failure)"),
        ("T/NK genes on platform?", "Yes — Cd3d/e, Nkg7, Klrd1 (not the failure)"),
        ("Lung Kras/Lkb1 GEMM?", "No — Ad-VEGF-A164 flank skin angiogenesis"),
        ("Epithelial/malignant lung cells?", "No — whole-flank bulk; no tumors"),
        ("Cell-level %pos or T/NK fraction?", "No — Affymetrix bulk MAS5"),
        ("T cells present?", "No — athymic nude (Foxn1nu)"),
        ("Same-study KP-KL deposit?", "No — sister GSE76588 is the same model"),
        ("Mouse-level Cldn4 vs T/NK scored?", "No — stop; table empty"),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(rows) + 0.5)
    ax.axis("off")
    ax.set_title("GSE76628 Cldn4–T/NK gate — NO-GO", loc="left", fontweight="bold")
    for i, (q, a) in enumerate(reversed(rows)):
        y = i + 0.5
        fail = a.startswith("No")
        color = "#F4C2C2" if fail else "#C6E5B3"
        ax.barh(y, 1.0, height=0.86, color=color, edgecolor="#333333", linewidth=0.4)
        ax.text(0.02, y, q, va="center", ha="left", fontsize=8, fontweight="bold")
        ax.text(0.98, y, a, va="center", ha="right", fontsize=8)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    matrix = CACHE / "GSE76628_series_matrix.txt.gz"
    sister = CACHE / "GSE76588_series_matrix.txt.gz"
    annot = CACHE / "GPL1261.annot.gz"
    fetch(MATRIX_URL, matrix)
    fetch(SISTER_MATRIX_URL, sister)
    fetch(ANNOT_URL, annot)

    meta = parse_series_matrix(matrix)
    sister_meta = parse_series_matrix(sister)
    probes = parse_probes(annot, PROBE_GENES)

    titles = meta["!Sample_title"]
    gsm = meta["!Sample_geo_accession"]
    sources = meta["!Sample_source_name_ch1"]
    n = len(gsm)

    sample_rows = []
    group_counter: Counter[str] = Counter()
    for i in range(n):
        day, treat, group = group_from_source(sources[i])
        group_counter[group] += 1
        sample_rows.append(
            {
                "gsm": gsm[i],
                "title": titles[i],
                "source_name": sources[i],
                "day": day,
                "treatment": treat,
                "group": group,
                "organism": "Mus musculus",
                "tissue": "flank skin",
                "model": "Ad-VEGF-A164 tumor-surrogate angiogenesis",
                "strain": "Crl:NU(NCr)-Foxn1nu athymic nude",
                "assay": "bulk Affymetrix Mouse430_2 MAS5",
                "epithelial_malignant_scored": "no",
                "cldn4_pctpos": "",
                "t_fraction": "",
                "nk_fraction": "",
            }
        )

    write_tsv(
        TABLES / "sample_inventory.tsv",
        sample_rows,
        [
            "gsm",
            "title",
            "source_name",
            "day",
            "treatment",
            "group",
            "organism",
            "tissue",
            "model",
            "strain",
            "assay",
            "epithelial_malignant_scored",
            "cldn4_pctpos",
            "t_fraction",
            "nk_fraction",
        ],
    )

    group_order = [
        "control_d0",
        "AdVEGF_d5_NT",
        "AdVEGF_d5_DC101",
        "AdVEGF_d5_G6",
        "AdVEGF_d20_NT",
        "AdVEGF_d20_DC101",
        "AdVEGF_d20_G6",
        "AdVEGF_d60_NT",
        "AdVEGF_d60_DC101",
        "AdVEGF_d60_G6",
    ]
    group_rows = [
        {
            "group": g,
            "n_arrays": group_counter[g],
            "inferential_unit": "one bulk array per mouse",
            "usable_for_cldn4_tnk": "no",
        }
        for g in group_order
    ]
    write_tsv(
        TABLES / "sample_groups.tsv",
        group_rows,
        ["group", "n_arrays", "inferential_unit", "usable_for_cldn4_tnk"],
    )

    write_tsv(
        TABLES / "mouse_units.tsv",
        [],
        [
            "mouse_id",
            "n_epithelial_or_malignant",
            "cldn4_pctpos",
            "cldn4_mean",
            "t_fraction_cd3de",
            "nk_fraction_nkg7_klrd1",
            "ifn_score",
            "mhc_score",
            "tj_score",
            "quartile",
        ],
    )

    write_tsv(
        TABLES / "probe_presence.tsv",
        [
            {
                "gene": g,
                "n_probes": len(probes[g]),
                "probes": ";".join(probes[g]),
                "used_to_score_cells": "no",
            }
            for g in PROBE_GENES
        ],
        ["gene", "n_probes", "probes", "used_to_score_cells"],
    )

    write_tsv(
        TABLES / "public_files.tsv",
        [
            {
                "accession": "GSE76628",
                "file": "GSE76628_series_matrix.txt.gz",
                "role": "processed bulk MAS5 matrix + metadata (used for inventory)",
                "bytes": matrix.stat().st_size,
                "url": MATRIX_URL,
                "downloaded": "yes",
            },
            {
                "accession": "GSE76628",
                "file": "GSE76628_RAW.tar",
                "role": "raw CEL; not a processed matrix; not used",
                "bytes": 305868800,
                "url": RAW_URL,
                "downloaded": "no",
            },
            {
                "accession": "GSE76588",
                "file": "GSE76588_series_matrix.txt.gz",
                "role": "sister processed matrix, same Ad-VEGF model (metadata only)",
                "bytes": sister.stat().st_size,
                "url": SISTER_MATRIX_URL,
                "downloaded": "yes",
            },
            {
                "accession": "GPL1261",
                "file": "GPL1261.annot.gz",
                "role": "probe-to-symbol map; Cldn4/T/NK presence only",
                "bytes": annot.stat().st_size,
                "url": ANNOT_URL,
                "downloaded": "yes",
            },
        ],
        ["accession", "file", "role", "bytes", "url", "downloaded"],
    )

    write_tsv(
        TABLES / "sister_series.tsv",
        [
            {
                "accession": "GSE76628",
                "relation": "this series (part II)",
                "n_samples": n,
                "organism": "Mus musculus",
                "assay": "bulk Affymetrix Mouse430_2",
                "model": "Ad-VEGF-A164 flank/ear angiogenesis, nude mice",
                "kp_kl_lung": "no",
            },
            {
                "accession": "GSE76588",
                "relation": "sister subseries (part I) of SuperSeries GSE76630",
                "n_samples": len(sister_meta.get("!Sample_geo_accession", [])),
                "organism": "Mus musculus",
                "assay": "bulk Affymetrix Mouse430_2",
                "model": "Ad-VEGF-A164 flank angiogenesis, same paper",
                "kp_kl_lung": "no",
            },
            {
                "accession": "GSE76630",
                "relation": "SuperSeries (PMID 27197264)",
                "n_samples": n + len(sister_meta.get("!Sample_geo_accession", [])),
                "organism": "Mus musculus",
                "assay": "bulk Affymetrix Mouse430_2",
                "model": "gastric-cancer stromal signatures via Ad-VEGF mouse model",
                "kp_kl_lung": "no",
            },
        ],
        [
            "accession",
            "relation",
            "n_samples",
            "organism",
            "assay",
            "model",
            "kp_kl_lung",
        ],
    )

    write_tsv(
        TABLES / "nogo_gate.tsv",
        [
            {"criterion": "public_processed_matrix", "result": "yes", "detail": "series matrix 45101 x 78"},
            {"criterion": "Cldn4_on_platform", "result": "yes", "detail": "1418283_at"},
            {"criterion": "T_NK_genes_on_platform", "result": "yes", "detail": "Cd3d/e Nkg7 Klrd1"},
            {"criterion": "lung_Kras_Lkb1_GEMM", "result": "no", "detail": "Ad-VEGF-A164 flank skin"},
            {"criterion": "epithelial_or_malignant_compartment", "result": "no", "detail": "bulk whole-flank; no tumors"},
            {"criterion": "cell_level_pctpos_or_fraction", "result": "no", "detail": "bulk microarray"},
            {"criterion": "T_cells_in_animals", "result": "no", "detail": "athymic nude Foxn1nu"},
            {"criterion": "same_study_KP_KL_deposit", "result": "no", "detail": "GSE76588/GSE76630 same model"},
            {"criterion": "mouse_units_scored", "result": "0", "detail": "stop; empty mouse_units.tsv"},
            {"criterion": "verdict", "result": "NO-GO", "detail": "do not score Cldn4 vs T/NK; not claim-failed"},
        ],
        ["criterion", "result", "detail"],
    )

    fig_sample_n([(g, group_counter[g]) for g in group_order], FIGURES / "fig_sample_n_by_group.png")
    fig_nogo_gate(FIGURES / "fig_nogo_gate.png")

    summary = {
        "accession": "GSE76628",
        "verdict": "NO-GO",
        "n_arrays": n,
        "n_mice_scored": 0,
        "n_epithelial_malignant": 0,
        "cldn4_tnk_scored": False,
        "same_study_kp_kl": False,
        "thesis_audited": False,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
