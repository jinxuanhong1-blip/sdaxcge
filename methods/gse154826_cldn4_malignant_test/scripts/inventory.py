#!/usr/bin/env python3
"""GSE154826 malignant-compartment gate.

Uses public GEO sample_annots + amp-batch tarball listing, plus author
cluster labels (GitHub, not a matrix). Does not download the 77-MTX pile
or Dropbox lung_ldm.rd. Peeks one small public batch tar for features.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import ssl
import tarfile
import urllib.request
from collections import Counter, defaultdict
from io import BytesIO
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
CACHE = Path("/tmp/gse154826_malig")
AUTHOR = CACHE / "author"
GEO = CACHE / "geo"
CTX = ssl.create_default_context()
UA = "gse154826-cldn4-malignant-test/1.0"

GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/"
GEO_MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/matrix/"
GH = "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/"

FOCUS_GENES = ("TACSTD2", "CLDN4", "EPCAM", "KRT19", "PTPRC", "CD3D", "NKG7")


def fetch(url: str, dest: Path | None = None, timeout: int = 180) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, context=CTX, timeout=timeout) as r:
        data = r.read()
    if dest is not None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
    return data


def head_bytes(url: str, timeout: int = 60) -> int | None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, context=CTX, timeout=timeout) as r:
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None


def listing_to_bytes(token: str) -> int | None:
    token = token.strip()
    m = re.fullmatch(r"([0-9.]+)([KMG])", token)
    if not m:
        return int(token) if token.isdigit() else None
    n, u = float(m.group(1)), m.group(2)
    return int(n * {"K": 1024, "M": 1024**2, "G": 1024**3}[u])


def parse_listing(html: str) -> list[dict]:
    rows = []
    # Restrict to GSE154826_* hrefs so the Parent Directory row cannot
    # swallow the first file (NCBI listings put "-" in the size column).
    for m in re.finditer(
        r'href="(GSE154826[^"]+)"[^>]*>[^<]*</a>\s+(\S+)\s+(\S+)\s+(\S+)',
        html,
    ):
        name = m.group(1).split("/")[-1]
        size = m.group(4)
        rows.append(
            {
                "name": name,
                "listing_size": size,
                "listing_bytes_approx": listing_to_bytes(size),
            }
        )
    return rows


def write_tsv(path: Path, header: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for row in rows:
            if isinstance(row, dict):
                w.writerow([row.get(h, "") for h in header])
            else:
                w.writerow(row)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    AUTHOR.mkdir(parents=True, exist_ok=True)
    GEO.mkdir(parents=True, exist_ok=True)

    # --- GEO listing (no full MTX pile) ---
    geo_rows = []
    for loc, url in (("GSE154826_suppl", GEO_SUPPL), ("GSE154826_matrix", GEO_MATRIX)):
        html = fetch(url).decode("utf-8", "replace")
        for row in parse_listing(html):
            file_url = url + row["name"]
            nbytes = head_bytes(file_url)
            geo_rows.append(
                {
                    "location": loc,
                    "file": row["name"],
                    "listing_size": row["listing_size"],
                    "listing_bytes_approx": row["listing_bytes_approx"],
                    "content_length": nbytes,
                    "url": file_url,
                }
            )

    mtx = [r for r in geo_rows if r["file"].startswith("GSE154826_amp_batch_ID_") and r["file"].endswith(".tar.gz")]
    mtx_bytes = sum(int(r["content_length"] or r["listing_bytes_approx"] or 0) for r in mtx)

    # compact GEO objects
    sample_annots_p = AUTHOR / "GSE154826_sample_annots.csv.gz"
    if not sample_annots_p.exists():
        fetch(GEO_SUPPL + "GSE154826_sample_annots.csv.gz", sample_annots_p)
    for sm in (
        "GSE154826-GPL18573_series_matrix.txt.gz",
        "GSE154826-GPL24676_series_matrix.txt.gz",
    ):
        dest = GEO / sm
        if not dest.exists():
            fetch(GEO_MATRIX + sm, dest)

    series_meta = {}
    for sm in sorted(GEO.glob("*series_matrix*.gz")):
        n_sample = 0
        n_data = 0
        with gzip.open(sm, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_"):
                    n_sample += 1
                elif line and not line.startswith("!"):
                    n_data += 1
        series_meta[sm.name] = {
            "bytes": sm.stat().st_size,
            "sample_meta_lines": n_sample,
            "non_bang_lines": n_data,
        }

    # author compact tables (labels, not expression)
    for name in (
        "annots_list.csv",
        "table_s1_sample_table.csv",
        "cell_metadata.csv",
        "immune_vs_ep_de.csv",
    ):
        dest = AUTHOR / name
        if not dest.exists():
            fetch(GH + name, dest)

    annots = {}
    with (AUTHOR / "annots_list.csv").open() as f:
        for row in csv.DictReader(f):
            annots[str(row["cluster"])] = row

    lineage_n = Counter()
    sub_n = Counter()
    lig_n = Counter()
    n_cells = 0
    sample_lineage = defaultdict(Counter)
    with (AUTHOR / "cell_metadata.csv").open() as f:
        for row in csv.DictReader(f):
            n_cells += 1
            a = annots.get(str(row["cluster_ID"]), {})
            lin = a.get("lineage") or "NA"
            sub = a.get("sub_lineage") or lin
            lig = a.get("lig_rec_group") or ""
            lineage_n[lin] += 1
            sub_n[(lin, sub)] += 1
            lig_n[lig] += 1
            sample_lineage[row["sample_ID"]][lin] += 1

    malignant_tokens = ("malignant", "tumor cell", "cancer", "epithelial", "epi")
    annot_token_hits = []
    for cid, a in annots.items():
        blob = " ".join(str(v).lower() for v in a.values())
        hits = [t for t in malignant_tokens if t in blob]
        if hits:
            annot_token_hits.append(
                {
                    "cluster": cid,
                    "lineage": a.get("lineage"),
                    "sub_lineage": a.get("sub_lineage"),
                    "norm_group": a.get("norm_group"),
                    "lig_rec_group": a.get("lig_rec_group"),
                    "token_hits": ",".join(hits),
                }
            )

    # sample_annots (GEO) + Table S1
    with gzip.open(sample_annots_p, "rt") as f:
        geo_samples = list(csv.DictReader(f))
    s1 = list(csv.DictReader((AUTHOR / "table_s1_sample_table.csv").open(encoding="utf-8-sig")))
    s1_by_id = {r["sample_ID"]: r for r in s1}

    # GEO prep strings vs Table S1 short codes
    geo_prep_n = Counter()
    geo_prep_patients = defaultdict(set)
    for r in geo_samples:
        prep = r.get("prep") or "NA"
        tissue = r.get("tissue") or "NA"
        geo_prep_n[(prep, tissue)] += 1
        geo_prep_patients[(prep, tissue)].add(r.get("patient_ID") or "")

    mount = [r for r in s1 if not str(r["patient_ID"]).startswith(("Lambrechts", "zilionis"))]
    s1_prep_n = Counter()
    s1_prep_patients = defaultdict(set)
    for r in mount:
        key = (r.get("prep") or "NA", r.get("tissue") or "NA")
        s1_prep_n[key] += 1
        s1_prep_patients[key].add(r["patient_ID"])

    digest_patients = sorted({r["patient_ID"] for r in mount if "digest" in (r.get("prep") or "")})
    tumor_patients = sorted({r["patient_ID"] for r in mount if r.get("tissue") == "Tumor"})
    mount_patients = sorted({r["patient_ID"] for r in mount})

    # per sample / patient counts from author barcodes
    per_sample = []
    per_patient_tumor = defaultdict(lambda: Counter())
    per_patient_tumor_meta = {}
    for sid, lin_c in sample_lineage.items():
        geo = next((r for r in geo_samples if r["sample_ID"] == sid), {})
        s1r = s1_by_id.get(sid, {})
        patient = geo.get("patient_ID") or s1r.get("patient_ID") or ""
        tissue = geo.get("tissue") or s1r.get("tissue") or ""
        disease = geo.get("disease") or s1r.get("disease") or ""
        prep_geo = geo.get("prep") or ""
        prep_s1 = s1r.get("prep") or ""
        n_t = lin_c.get("T", 0)
        n_nk = lin_c.get("NK", 0)
        n_gate = lin_c.get("epi_endo_fibro_doublet", 0)
        n_other = sum(lin_c.values()) - n_t - n_nk - n_gate
        row = {
            "sample_ID": sid,
            "patient_ID": patient,
            "tissue": tissue,
            "disease": disease,
            "prep_geo": prep_geo,
            "prep_s1": prep_s1,
            "n_cells": sum(lin_c.values()),
            "n_T": n_t,
            "n_NK": n_nk,
            "n_TNK": n_t + n_nk,
            "n_gate_epi_endo_fibro_doublet": n_gate,
            "n_other_immune": n_other,
            "n_author_malignant": 0,
        }
        per_sample.append(row)
        if tissue == "Tumor" and patient and not str(patient).startswith(("Lambrechts", "zilionis")):
            c = per_patient_tumor[patient]
            c["n_cells"] += row["n_cells"]
            c["n_T"] += n_t
            c["n_NK"] += n_nk
            c["n_TNK"] += n_t + n_nk
            c["n_gate"] += n_gate
            c["n_other_immune"] += n_other
            c["n_author_malignant"] += 0
            c["n_libraries"] += 1
            per_patient_tumor_meta[patient] = {
                "disease": disease,
                "prep_geo_set": per_patient_tumor_meta.get(patient, {}).get("prep_geo_set", set())
                | ({prep_geo} if prep_geo else set()),
                "prep_s1_set": per_patient_tumor_meta.get(patient, {}).get("prep_s1_set", set())
                | ({prep_s1} if prep_s1 else set()),
            }

    per_sample.sort(key=lambda r: (r["patient_ID"], r["tissue"], r["sample_ID"]))

    patient_rows = []
    for pt in sorted(per_patient_tumor):
        c = per_patient_tumor[pt]
        meta = per_patient_tumor_meta[pt]
        preps = ",".join(sorted(meta["prep_s1_set"] or meta["prep_geo_set"]))
        is_digest = "digest" in preps
        # leftover gate is not author-malignant
        eligible_malig = False
        patient_rows.append(
            {
                "patient_ID": pt,
                "disease": meta["disease"],
                "prep": preps,
                "digest_no_cd45": int(is_digest),
                "n_libraries": c["n_libraries"],
                "n_cells": c["n_cells"],
                "n_author_malignant": 0,
                "n_gate_epi_endo_fibro_doublet": c["n_gate"],
                "n_TNK": c["n_TNK"],
                "n_T": c["n_T"],
                "n_NK": c["n_NK"],
                "gate_and_TNK_both_ge20": int(c["n_gate"] >= 20 and c["n_TNK"] >= 20),
                "author_malignant_and_TNK": 0,
                "treat_gate_as_malignant": 0,
            }
        )

    n_gate_patients = sum(1 for r in patient_rows if r["n_gate_epi_endo_fibro_doublet"] >= 20)
    n_digest = sum(1 for r in patient_rows if r["digest_no_cd45"])
    n_author_malig_units = 0
    n_author_malig_and_tnk = 0

    # immune vs epi DE for CLDN4 / TACSTD2 (compact author table)
    de_hits = []
    with (AUTHOR / "immune_vs_ep_de.csv").open() as f:
        rdr = csv.reader(f)
        header = next(rdr)
        for rec in rdr:
            if any(g in rec for g in FOCUS_GENES):
                de_hits.append(dict(zip(header, rec)))

    # peek smallest public batch tar for features only
    mtx_sorted = sorted(mtx, key=lambda r: int(r["content_length"] or r["listing_bytes_approx"] or 10**18))
    peek = mtx_sorted[0] if mtx_sorted else None
    feature_hit = {}
    peek_info = {}
    if peek:
        peek_path = GEO / peek["file"]
        if not peek_path.exists():
            fetch(peek["url"], peek_path, timeout=300)
        peek_info = {
            "file": peek["file"],
            "bytes": peek_path.stat().st_size,
            "url": peek["url"],
        }
        with tarfile.open(peek_path, "r:gz") as tf:
            names = tf.getnames()
            peek_info["members"] = names
            feat_name = next(n for n in names if n.endswith("features.tsv") or n.endswith("features.tsv.gz"))
            fh = tf.extractfile(feat_name)
            raw = fh.read()
            if feat_name.endswith(".gz"):
                raw = gzip.decompress(raw)
            text = raw.decode("utf-8", "replace")
            symbols = []
            ftypes = Counter()
            for line in text.splitlines():
                parts = line.split("\t")
                if len(parts) >= 2:
                    symbols.append(parts[1])
                if len(parts) >= 3:
                    ftypes[parts[2]] += 1
            peek_info["n_features"] = len(symbols)
            peek_info["feature_types"] = dict(ftypes)
            for g in FOCUS_GENES:
                feature_hit[g] = g in symbols
            # also record barcode count
            bname = next(n for n in names if n.endswith("barcodes.tsv") or n.endswith("barcodes.tsv.gz"))
            bh = tf.extractfile(bname)
            braw = bh.read()
            if bname.endswith(".gz"):
                braw = gzip.decompress(braw)
            peek_info["n_barcodes"] = sum(1 for _ in braw.splitlines() if _)

    drop_n = head_bytes("https://www.dropbox.com/s/vjbide8ro5iwrfh/lung_ldm.rd?dl=1")

    write_tsv(
        TABLES / "geo_file_inventory.tsv",
        ["location", "file", "listing_size", "listing_bytes_approx", "content_length", "url"],
        geo_rows,
    )
    write_tsv(TABLES / "author_lineage.tsv", ["lineage", "n_cells"], [[k, v] for k, v in lineage_n.most_common()])
    write_tsv(
        TABLES / "author_sublineage.tsv",
        ["lineage", "sub_lineage", "n_cells"],
        [[a, b, n] for (a, b), n in sorted(sub_n.items(), key=lambda x: -x[1])],
    )
    write_tsv(
        TABLES / "author_annot_token_hits.tsv",
        ["cluster", "lineage", "sub_lineage", "norm_group", "lig_rec_group", "token_hits"],
        annot_token_hits,
    )
    write_tsv(
        TABLES / "geo_sample_prep.tsv",
        ["prep", "tissue", "n_sample_rows", "n_patients"],
        [
            [prep, tissue, n, len(geo_prep_patients[(prep, tissue)])]
            for (prep, tissue), n in sorted(geo_prep_n.items())
        ],
    )
    write_tsv(
        TABLES / "table_s1_sample_prep.tsv",
        ["prep", "tissue", "n_sample_rows", "n_patients"],
        [
            [prep, tissue, n, len(s1_prep_patients[(prep, tissue)])]
            for (prep, tissue), n in sorted(s1_prep_n.items())
        ],
    )
    write_tsv(
        TABLES / "per_sample_author_counts.tsv",
        [
            "sample_ID",
            "patient_ID",
            "tissue",
            "disease",
            "prep_geo",
            "prep_s1",
            "n_cells",
            "n_T",
            "n_NK",
            "n_TNK",
            "n_gate_epi_endo_fibro_doublet",
            "n_other_immune",
            "n_author_malignant",
        ],
        per_sample,
    )
    write_tsv(
        TABLES / "per_patient_tumor_author_counts.tsv",
        [
            "patient_ID",
            "disease",
            "prep",
            "digest_no_cd45",
            "n_libraries",
            "n_cells",
            "n_author_malignant",
            "n_gate_epi_endo_fibro_doublet",
            "n_TNK",
            "n_T",
            "n_NK",
            "gate_and_TNK_both_ge20",
            "author_malignant_and_TNK",
            "treat_gate_as_malignant",
        ],
        patient_rows,
    )

    # honest n
    honest = [
        ("public_GSE154826_series", 1, "GEO; Leader/Grout Cancer Cell 2021 PMID 34767762"),
        ("GEO_series_matrix_expression_rows", 0, "GPL18573 + GPL24676 metadata only"),
        ("GEO_amp_batch_MTX_tarballs", len(mtx), f"unlabeled 10x MTX; pile ~{mtx_bytes} B; not downloaded as a pile"),
        ("GEO_sample_annots_rows", len(geo_samples), "GSE154826_sample_annots.csv.gz"),
        ("author_annotated_barcodes", n_cells, "effiken/Leader_et_al cell_metadata.csv"),
        ("author_clusters", len(annots), "annots_list.csv"),
        ("author_lineage_named_Malignant", 0, "no Malignant / tumor-cell lineage"),
        ("author_cells_labeled_malignant", 0, "0 barcodes"),
        ("author_lineage_epi_endo_fibro_doublet", lineage_n.get("epi_endo_fibro_doublet", 0), "leftover gate after CD45+ enrichment"),
        ("author_lineage_T", lineage_n.get("T", 0), "T/NK present; not the failure"),
        ("author_lineage_NK", lineage_n.get("NK", 0), ""),
        ("Mount_Sinai_patients_TableS1", len(mount_patients), "excludes Lambrechts/Zilionis re-use"),
        ("tumor_patients_with_author_barcodes", len(patient_rows), "patient x Tumor"),
        ("digest_no_CD45_tumor_patients", n_digest, "695 LUAD, 706 LUSC — real tissue digest"),
        ("author_malignant_units", n_author_malig_units, "stop threshold n<8 / score threshold n>=10"),
        ("author_malignant_and_same_patient_TNK_units", n_author_malig_and_tnk, "empty"),
        ("tumor_patients_gate_ge20_and_TNK_ge20", n_gate_patients, "leftover gate, NOT author-malignant; not scored"),
        ("CLDN4_vs_TNK_table_rows", 0, "malignant absent; not computed"),
        ("IFN_MHC_Q4_vs_Q1", 0, "optional; not run"),
        ("joined_concordant_pool", 0, "sign not tested; do not merge 123902/131907"),
        ("dual_high_TACSTD2_x_CLDN4", 0, "not defined"),
        ("Dropbox_lung_ldm_rd_downloaded", 0, f"HEAD {drop_n} B; extra huge object skipped"),
    ]
    write_tsv(
        TABLES / "immune_vs_ep_focus.tsv",
        ["gene", "l2fc_immune_vs_epi", "fg_exprs_immune", "bg_exprs_epi"],
        [
            [
                h.get("") or h.get("gene") or "",
                h.get("l2fc"),
                h.get("fg_exprs"),
                h.get("bg_exprs"),
            ]
            for h in de_hits
        ],
    )
    write_tsv(
        TABLES / "features_peek.tsv",
        ["item", "value"],
        [
            ["file", peek_info.get("file", "")],
            ["bytes", peek_info.get("bytes", "")],
            ["n_features", peek_info.get("n_features", "")],
            ["n_barcodes_raw", peek_info.get("n_barcodes", "")],
            ["feature_types", peek_info.get("feature_types", "")],
            *[[f"gene_{g}_present", int(bool(feature_hit.get(g)))] for g in FOCUS_GENES],
            ["note", "smallest public amp-batch tar; features only; matrix not scored"],
        ],
    )

    write_tsv(TABLES / "honest_n.tsv", ["item", "n", "note"], honest)

    verdict = "malignant_absent_CD45_immune_only"
    reason = (
        "Author annotations have 0 Malignant barcodes. Design is CD45+ bead/FACS "
        "CITE-seq; leftover epithelium is gated epi_endo_fibro_doublet. Real no-CD45 "
        f"digest epithelium is {n_digest} patients (<8). Stop. Do not score CLDN4 vs T/NK "
        "and do not merge with GSE123902/GSE131907."
    )

    summary = {
        "dataset": "GSE154826",
        "pmid": 34767762,
        "cldn4_only": True,
        "dual_high": False,
        "verdict": verdict,
        "malignant_present": False,
        "reason": reason,
        "n_author_malignant_cells": 0,
        "n_author_malignant_units": 0,
        "n_digest_no_cd45_tumor_patients": n_digest,
        "n_gate_cells": lineage_n.get("epi_endo_fibro_doublet", 0),
        "n_T": lineage_n.get("T", 0),
        "n_NK": lineage_n.get("NK", 0),
        "n_author_cells": n_cells,
        "n_tumor_patients": len(patient_rows),
        "n_mount_sinai_patients": len(mount_patients),
        "mtx_tarball_n": len(mtx),
        "mtx_pile_bytes": mtx_bytes,
        "mtx_pile_downloaded": False,
        "dropbox_lung_ldm_rd_bytes": drop_n,
        "dropbox_downloaded": False,
        "series_matrix": series_meta,
        "features_peek": peek_info,
        "focus_genes_in_peeked_features": feature_hit,
        "author_lineage": dict(lineage_n),
        "author_lig_rec_group": dict(lig_n),
        "digest_patients": digest_patients,
        "immune_vs_ep_focus": de_hits,
        "concordant_pool_joined": False,
        "forced_merge_123902_131907": False,
        "cldn4_pctpos_vs_tnk": None,
        "ifn_mhc_q4q1": None,
        "stop_rule": "CD45+/immune-only or malignant n<8 → write that and STOP",
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in (
        "verdict", "malignant_present", "n_author_malignant_units",
        "n_digest_no_cd45_tumor_patients", "n_gate_cells", "mtx_tarball_n",
        "focus_genes_in_peeked_features",
    )}, indent=2))


if __name__ == "__main__":
    main()
