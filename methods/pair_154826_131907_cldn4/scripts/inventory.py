#!/usr/bin/env python3
"""Honest inventory for GSE154826 + GSE131907 CLDN4 pair gate. No expression matrices opened."""

from __future__ import annotations

import csv
import gzip
import json
import re
import ssl
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
CACHE = Path("/tmp/pair_154826_131907")
AUTHOR = CACHE / "author"
CTX = ssl.create_default_context()


def fetch(url: str, dest: Path | None = None, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "pair-154826-inventory/1.0"})
    with urllib.request.urlopen(req, context=CTX, timeout=timeout) as r:
        data = r.read()
    if dest is not None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
    return data


def head_bytes(url: str, timeout: int = 60) -> int | None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "pair-154826-inventory/1.0"})
    try:
        with urllib.request.urlopen(req, context=CTX, timeout=timeout) as r:
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None


def parse_listing_sizes(html: str) -> list[dict]:
    rows = []
    for m in re.finditer(
        r'href="([^"]+)"[^<]*</a>\s+\S+\s+\S+\s+(\S+)',
        html,
    ):
        name, size = m.group(1), m.group(2)
        if name in ("Parent Directory",) or name.endswith("/"):
            continue
        rows.append({"name": name.split("/")[-1], "listing_size": size})
    return rows


def listing_to_bytes(token: str) -> int | None:
    token = token.strip()
    m = re.fullmatch(r"([0-9.]+)([KMG])", token)
    if not m:
        if token.isdigit():
            return int(token)
        return None
    n, u = float(m.group(1)), m.group(2)
    mul = {"K": 1024, "M": 1024**2, "G": 1024**3}[u]
    return int(n * mul)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    AUTHOR.mkdir(parents=True, exist_ok=True)

    listings = {
        "GSE154826_suppl": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/",
        "GSE154826_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/matrix/",
        "GSE131907_suppl": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/",
    }
    geo_rows = []
    for series_key, url in listings.items():
        html = fetch(url).decode("utf-8", "replace")
        for row in parse_listing_sizes(html):
            file_url = url + row["name"]
            nbytes = head_bytes(file_url)
            geo_rows.append(
                {
                    "series": series_key.split("_")[0],
                    "location": series_key,
                    "file": row["name"],
                    "listing_size": row["listing_size"],
                    "listing_bytes_approx": listing_to_bytes(row["listing_size"]),
                    "content_length": nbytes,
                    "url": file_url,
                }
            )

    seen = {(r["location"], r["file"]) for r in geo_rows}
    extras = [
        (
            "GSE154826",
            "GSE154826_suppl",
            "GSE154826_amp_batch_ID_1.tar.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/GSE154826_amp_batch_ID_1.tar.gz",
        ),
        (
            "GSE154826",
            "GSE154826_matrix",
            "GSE154826-GPL18573_series_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/matrix/GSE154826-GPL18573_series_matrix.txt.gz",
        ),
        (
            "GSE131907",
            "GSE131907_suppl",
            "GSE131907_Lung_Cancer_Feature_Summary.xlsx",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_Feature_Summary.xlsx",
        ),
    ]
    for series, loc, name, file_url in extras:
        if (loc, name) in seen:
            continue
        nbytes = head_bytes(file_url)
        geo_rows.append(
            {
                "series": series,
                "location": loc,
                "file": name,
                "listing_size": str(nbytes) if nbytes else "unknown",
                "listing_bytes_approx": nbytes,
                "content_length": nbytes,
                "url": file_url,
            }
        )

    dropbox = "https://www.dropbox.com/s/vjbide8ro5iwrfh/lung_ldm.rd?dl=1"
    drop_n = head_bytes(dropbox)
    geo_rows.append(
        {
            "series": "GSE154826",
            "location": "author_Dropbox_not_GEO",
            "file": "lung_ldm.rd",
            "listing_size": "unknown" if drop_n is None else str(drop_n),
            "listing_bytes_approx": drop_n,
            "content_length": drop_n,
            "url": dropbox,
        }
    )

    with (TABLES / "geo_file_inventory.tsv").open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "series",
                "location",
                "file",
                "listing_size",
                "listing_bytes_approx",
                "content_length",
                "url",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(geo_rows)

    # --- GSE154826 author metadata ---
    annots_path = AUTHOR / "annots_list.csv"
    if not annots_path.exists():
        fetch(
            "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/annots_list.csv",
            annots_path,
        )
    s1_path = AUTHOR / "table_s1_sample_table.csv"
    if not s1_path.exists():
        fetch(
            "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/table_s1_sample_table.csv",
            s1_path,
        )
    meta_path = AUTHOR / "cell_metadata.csv"
    if not meta_path.exists():
        fetch(
            "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/cell_metadata.csv",
            meta_path,
        )
    de_path = AUTHOR / "immune_vs_ep_de.csv"
    if not de_path.exists():
        fetch(
            "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/immune_vs_ep_de.csv",
            de_path,
        )
    sample_annots = AUTHOR / "GSE154826_sample_annots.csv.gz"
    if not sample_annots.exists():
        fetch(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/GSE154826_sample_annots.csv.gz",
            sample_annots,
        )

    annots = {}
    with annots_path.open() as f:
        for row in csv.DictReader(f):
            annots[str(row["cluster"])] = row

    lineage_n = Counter()
    sub_n = Counter()
    n_cells = 0
    n_samples_meta = set()
    with meta_path.open() as f:
        for row in csv.DictReader(f):
            n_cells += 1
            n_samples_meta.add(row["sample_ID"])
            a = annots.get(str(row["cluster_ID"]), {})
            lin = a.get("lineage") or "NA"
            sub = a.get("sub_lineage") or lin
            lineage_n[lin] += 1
            sub_n[f"{lin}|{sub}"] += 1

    with (TABLES / "gse154826_lineage.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["lineage", "n_cells"])
        for k, v in lineage_n.most_common():
            w.writerow([k, v])

    with (TABLES / "gse154826_sublineage.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["lineage", "sub_lineage", "n_cells"])
        for k, v in sub_n.most_common():
            lin, sub = k.split("|", 1)
            w.writerow([lin, sub, v])

    s1_rows = list(csv.DictReader(s1_path.open()))
    mount_sinai = [
        r
        for r in s1_rows
        if not str(r["patient_ID"]).startswith(("Lambrechts", "zilionis"))
    ]
    patients = sorted({r["patient_ID"] for r in mount_sinai})
    tumor_patients = sorted({r["patient_ID"] for r in mount_sinai if r["tissue"] == "Tumor"})
    digest_patients = sorted(
        {r["patient_ID"] for r in mount_sinai if "digest" in (r.get("prep") or "")}
    )
    beads_tumor = [
        r
        for r in mount_sinai
        if r["tissue"] == "Tumor" and (r.get("prep") or "") == "beads"
    ]
    cite_hint = [
        r
        for r in mount_sinai
        if r.get("HTO") not in (None, "", "NA")
    ]

    with (TABLES / "gse154826_sample_prep.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["prep", "tissue", "n_sample_rows", "n_patients"])
        key_n = Counter((r.get("prep") or "NA", r.get("tissue") or "NA") for r in mount_sinai)
        key_p = {}
        for r in mount_sinai:
            key_p.setdefault((r.get("prep") or "NA", r.get("tissue") or "NA"), set()).add(r["patient_ID"])
        for (prep, tissue), n in sorted(key_n.items()):
            w.writerow([prep, tissue, n, len(key_p[(prep, tissue)])])

    # CLDN4 in author immune vs epithelial DE (compact table, not a matrix)
    cldn4 = None
    with de_path.open() as f:
        for row in csv.DictReader(f):
            gene = row.get("gene") or row.get("Gene") or row.get("") or next(iter(row.values()))
            # file may use first column as gene
            if "CLDN4" in row.values() or gene == "CLDN4":
                cldn4 = row
                break
        if cldn4 is None:
            f.seek(0)
            rdr = csv.reader(f)
            header = next(rdr)
            for rec in rdr:
                if "CLDN4" in rec:
                    cldn4 = dict(zip(header, rec))
                    break

    # series matrices: confirm metadata-only
    sm_dir = CACHE / "geo154826"
    series_n_data_rows = {}
    for sm in sorted(sm_dir.glob("*series_matrix*.gz")):
        n_sample_line = 0
        n_data = 0
        with gzip.open(sm, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_"):
                    n_sample_line += 1
                elif line.startswith("ID_REF") or (line and not line.startswith("!")):
                    n_data += 1
        series_n_data_rows[sm.name] = {"sample_meta_lines": n_sample_line, "non_bang_lines": n_data}

    # GSE131907 annotation only (1.8M) — confirm T/NK + epithelial labels exist; do not open UMI
    ann131 = CACHE / "geo131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    if not ann131.exists():
        fetch(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
            ann131,
        )
    g131_cells = 0
    g131_samples = set()
    g131_origin = Counter()
    g131_celltype = Counter()
    g131_refined = Counter()
    g131_subtype = Counter()
    g131_cols = None
    tlung_malig = 0
    tlung_epi = 0
    tlung_tnk = 0
    mal_by_origin = Counter()
    samples_mal_tnk = set()
    samples_mal = set()
    samples_tnk = set()
    per_sample = {}
    with gzip.open(ann131, "rt") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        g131_cols = rdr.fieldnames
        for row in rdr:
            g131_cells += 1
            samp = row.get("Sample") or ""
            if samp:
                g131_samples.add(samp)
            origin = row.get("Sample_Origin") or "NA"
            ct = row.get("Cell_type") or "NA"
            ref = row.get("Cell_type.refined") or "NA"
            sub = row.get("Cell_subtype") or "NA"
            g131_origin[origin] += 1
            g131_celltype[ct] += 1
            g131_refined[ref] += 1
            g131_subtype[sub] += 1
            if origin == "tLung":
                if sub == "Malignant cells":
                    tlung_malig += 1
                if ct == "Epithelial cells":
                    tlung_epi += 1
                if ct in ("T lymphocytes", "NK cells"):
                    tlung_tnk += 1
            if sub == "Malignant cells":
                mal_by_origin[origin] += 1
                samples_mal.add(samp)
                per_sample.setdefault(samp, {"mal": 0, "tnk": 0})
                per_sample[samp]["mal"] += 1
            if ct in ("T lymphocytes", "NK cells"):
                samples_tnk.add(samp)
                per_sample.setdefault(samp, {"mal": 0, "tnk": 0})
                per_sample[samp]["tnk"] += 1

    with (TABLES / "gse131907_celltype.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Cell_type", "n_cells"])
        for k, v in g131_celltype.most_common():
            w.writerow([k, v])
    with (TABLES / "gse131907_cell_subtype.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Cell_subtype", "n_cells"])
        for k, v in g131_subtype.most_common():
            w.writerow([k, v])
    with (TABLES / "gse131907_sample_origin.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Sample_Origin", "n_cells"])
        for k, v in g131_origin.most_common():
            w.writerow([k, v])
    with (TABLES / "gse131907_malignant_by_origin.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Sample_Origin", "n_malignant"])
        for k, v in mal_by_origin.most_common():
            w.writerow([k, v])
    n_samp_mal_tnk = sum(1 for v in per_sample.values() if v["mal"] > 0 and v["tnk"] > 0)

    mtx_pile = [
        r
        for r in geo_rows
        if r["series"] == "GSE154826" and r["file"].startswith("GSE154826_amp_batch_ID_")
    ]
    mtx_bytes = sum((r["content_length"] or r["listing_bytes_approx"] or 0) for r in mtx_pile)
    g131_eligible = [
        r
        for r in geo_rows
        if r["series"] == "GSE131907"
        and r["file"].endswith((".rds.gz", ".txt.gz"))
        and "matrix" in r["file"]
    ]

    tnk = lineage_n.get("T", 0) + lineage_n.get("NK", 0)
    gate = lineage_n.get("epi_endo_fibro_doublet", 0)
    malignant_labeled = 0  # author table has no Malignant label

    honest = [
        ("public_GSE154826_series", 1, "GEO 2020; Leader/Grout Cancer Cell 2021 PMID 34767762"),
        ("GSE154826_GEO_labeled_processed_TME_matrix", 0, "series matrices are metadata-only; MTX pile unlabeled"),
        ("GSE154826_GEO_series_matrix_expression_rows", 0, "GPL18573 8.6K + GPL24676 6.3K metadata"),
        ("GSE154826_GEO_amp_batch_MTX_tarballs", len(mtx_pile), f"unlabeled 10x feature-barcode MTX; pile ~{mtx_bytes} B"),
        ("GSE154826_author_annotated_barcodes", n_cells, "effiken/Leader_et_al cell_metadata.csv"),
        ("GSE154826_author_samples_in_cell_metadata", len(n_samples_meta), "sample_ID in cell_metadata"),
        ("GSE154826_Mount_Sinai_patients_TableS1", len(patients), "excludes Lambrechts/Zilionis public re-use rows"),
        ("GSE154826_Mount_Sinai_tumor_patients", len(tumor_patients), "tissue=Tumor"),
        ("GSE154826_digest_patients_real_epithelium", len(digest_patients), "prep=digest/dead cell; 695 LUAD, 706 LUSC"),
        ("GSE154826_beads_tumor_sample_rows", len(beads_tumor), "CD45+ beads; epithelium gated out"),
        ("GSE154826_cells_lineage_T", lineage_n.get("T", 0), "author cluster lineage"),
        ("GSE154826_cells_lineage_NK", lineage_n.get("NK", 0), "author cluster lineage"),
        ("GSE154826_cells_lineage_T_plus_NK", tnk, "T/NK present; not the failure"),
        ("GSE154826_cells_epi_endo_fibro_doublet_gate", gate, "not labeled malignant; leftover after CD45 enrichment"),
        ("GSE154826_cells_labeled_malignant", malignant_labeled, "no Malignant field in annots_list.csv"),
        ("GSE154826_usable_malignant_CLDN4_matrix", 0, "stop"),
        ("GSE131907_size_eligible_processed_matrices_lt_2GB", sum(1 for r in g131_eligible if (r["content_length"] or r["listing_bytes_approx"] or 0) < 2 * 1024**3), "raw UMI txt/rds and log2TPM rds; log2TPM txt is 2.9G"),
        ("GSE131907_annotation_cells", g131_cells, "opened annotation only; UMI not downloaded"),
        ("GSE131907_annotation_samples", len(g131_samples), "58 libraries; paper 44 patients"),
        ("GSE131907_epithelial_cells", g131_celltype.get("Epithelial cells", 0), "Cell_type"),
        ("GSE131907_malignant_subtype_cells", g131_subtype.get("Malignant cells", 0), "Cell_subtype; not scored for CLDN4 here"),
        ("GSE131907_T_lymphocytes", g131_celltype.get("T lymphocytes", 0), "Cell_type"),
        ("GSE131907_NK_cells", g131_celltype.get("NK cells", 0), "Cell_type"),
        ("GSE131907_tLung_malignant", tlung_malig, "author Malignant cells are mBrain/tL/B/mLN, not tLung"),
        ("GSE131907_tLung_epithelial", tlung_epi, "primary tumor epithelial barcodes"),
        ("GSE131907_tLung_T_plus_NK", tlung_tnk, "primary tumor T+NK barcodes"),
        ("GSE131907_samples_malignant_and_TNK", n_samp_mal_tnk, "annotation only; would have been pair-eligible if 154826 passed"),
        ("pair_patients_154826_plus_131907", 0, "154826 gate failed; 131907 matrix not opened"),
        ("patient_level_malignant_CLDN4_vs_TNK", 0, "empty"),
        ("CellChat_style_CLDN4high_to_TNK", 0, "not run; rho/Q4 not computed"),
        ("dual_high_TACSTD2_CLDN4", 0, "not defined"),
    ]
    with (TABLES / "honest_n.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["item", "n", "note"])
        w.writerows(honest)

    summary = {
        "verdict": "stop",
        "reason": "GSE154826 has no public labeled processed matrix <2GB usable for patient-level malignant CLDN4 vs T/NK",
        "pair": ["GSE154826", "GSE131907"],
        "cldn4_only": True,
        "dual_high": False,
        "gse154826": {
            "pmid": 34767762,
            "design": "CD45+ immune-enriched CITE-seq / scRNA, Mount Sinai NSCLC",
            "author_cells": n_cells,
            "lineage": dict(lineage_n),
            "mount_sinai_patients": len(patients),
            "tumor_patients": len(tumor_patients),
            "digest_patients": digest_patients,
            "beads_tumor_rows": len(beads_tumor),
            "hashed_sample_rows": len(cite_hint),
            "mtx_tarball_n": len(mtx_pile),
            "mtx_pile_bytes_approx": mtx_bytes,
            "series_matrix": series_n_data_rows,
            "immune_vs_ep_CLDN4": cldn4,
            "dropbox_lung_ldm_rd_bytes": drop_n,
        },
        "gse131907": {
            "annotation_cells": g131_cells,
            "annotation_columns": g131_cols,
            "annotation_n_samples": len(g131_samples),
            "cell_type": dict(g131_celltype),
            "cell_subtype_malignant": g131_subtype.get("Malignant cells", 0),
            "tLung_malignant": tlung_malig,
            "tLung_epithelial": tlung_epi,
            "tLung_T_plus_NK": tlung_tnk,
            "malignant_by_origin": dict(mal_by_origin),
            "samples_malignant_and_TNK": n_samp_mal_tnk,
            "matrix_opened": False,
        },
        "assigned_tests": {
            "pair_rho": None,
            "pair_Q4": None,
            "CellChat_style": "not_run",
        },
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("verdict", "reason")}, indent=2))
    print("cells", n_cells, "lineage", dict(lineage_n))
    print("CLDN4 DE", cldn4)
    print("131907 cells", g131_cells, "cell_type", dict(g131_celltype))
    print("mtx pile bytes", mtx_bytes, "dropbox", drop_n)


if __name__ == "__main__":
    main()
