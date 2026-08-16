#!/usr/bin/env python3
"""
Step 5: assemble the master catalog (notes/geo_sweep/catalog.tsv).

Merges the outputs of steps 1-4 (search metadata, matrix probe, matrix inspect,
download manifest) into one row per verified extra accession and classifies each
series as ICI-labeled and/or KD/KO, plus which marker terms it mentions.
"""
import json
import re
from pathlib import Path

NOTES = Path("notes/geo_sweep")

KD_KO_RE = re.compile(
    r"knock[\s-]?out|knock[\s-]?down|\bKO\b|\bKD\b|shRNA|siRNA|CRISPR|sgRNA|"
    r"-/-|-\\-|deficient|\bnull\b|floxed|conditional (?:allele|knockout)|"
    r"gene (?:deletion|silencing)|overexpress",
    re.IGNORECASE,
)
ICI_RE = re.compile(
    r"immune checkpoint|checkpoint (?:inhibitor|blockade)|\bICI\b|\bICB\b|"
    r"anti[\s-]?PD[\s-]?L?1|\bPD-1\b|\bPD-L1\b|CTLA-?4|pembrolizumab|nivolumab|"
    r"atezolizumab|durvalumab|ipilimumab|cemiplimab|sintilimab|tislelizumab|"
    r"camrelizumab|immunotherap",
    re.IGNORECASE,
)
MARKERS = ["TACSTD2", "TROP2", "CLDN4", "claudin-4"]


def load(name):
    return json.load(open(NOTES / name))


def main():
    search = load("search_raw.json")["records"]
    probe = {r["accession"]: r for r in load("matrix_probe.json")}
    inspect_list = load("matrix_inspect.json")
    inspect = {}
    for r in inspect_list:
        if "series_title" in r:
            inspect.setdefault(r["accession"], r)  # first (primary) file
    manifest = {m["accession"]: m for m in load("download_manifest.json")}

    rows = []
    for rec in search:
        acc = rec["Accession"]
        pr = probe.get(acc, {})
        ins = inspect.get(acc, {})
        man = manifest.get(acc, {})

        blob = " ".join([
            rec.get("title", "") or "",
            rec.get("summary", "") or "",
            ins.get("series_overall_design", "") or "",
            " ".join(ins.get("sample_titles", []) or []),
            " ".join(ins.get("sample_characteristics", []) or []),
        ])
        is_kd_ko = bool(KD_KO_RE.search(blob))
        matched = rec.get("_matched_terms", [])
        is_ici = ("ICI" in matched) or bool(ICI_RE.search(blob))
        marker_terms = [m for m in MARKERS
                        if re.search(re.escape(m), blob, re.IGNORECASE)
                        or m in matched]

        cats = []
        if is_ici:
            cats.append("ICI")
        if is_kd_ko:
            cats.append("KD/KO")
        if marker_terms and not cats:
            cats.append("marker-only")
        if not cats:
            cats.append("other")

        matrix_files = [f["name"] for f in pr.get("files", [])]
        rows.append({
            "accession": acc,
            "matched_terms": ";".join(matched),
            "category": ";".join(cats),
            "is_ici": "Y" if is_ici else "N",
            "is_kd_ko": "Y" if is_kd_ko else "N",
            "marker_terms": ";".join(marker_terms),
            "organism": rec.get("taxon", ""),
            "gds_type": rec.get("gdsType", ""),
            "n_samples": rec.get("n_samples", ""),
            "platform": (rec.get("GPL") or ""),
            "has_expression_matrix": "Y" if ins.get("has_expression_table") else "N",
            "matrix_n_rows": ins.get("n_data_rows", 0),
            "matrix_total_bytes": pr.get("total_size", 0),
            "under_2gb": "Y" if pr.get("under_2gb") else "N",
            "downloaded": "Y" if man.get("status") == "ok" else "N",
            "matrix_files": ";".join(matrix_files),
            "pubmed_ids": ";".join(
                (rec.get("PubMedIds") or "").split() if isinstance(rec.get("PubMedIds"), str)
                else []),
            "ftp_matrix_url": pr.get("matrix_dir", ""),
            "title": (rec.get("title") or "").replace("\t", " ").strip(),
        })

    rows.sort(key=lambda r: r["accession"])
    cols = ["accession", "matched_terms", "category", "is_ici", "is_kd_ko",
            "marker_terms", "organism", "gds_type", "n_samples", "platform",
            "has_expression_matrix", "matrix_n_rows", "matrix_total_bytes",
            "under_2gb", "downloaded", "matrix_files", "pubmed_ids",
            "ftp_matrix_url", "title"]
    out = NOTES / "catalog.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_ici = sum(1 for r in rows if r["is_ici"] == "Y")
    n_kd = sum(1 for r in rows if r["is_kd_ko"] == "Y")
    n_expr = sum(1 for r in rows if r["has_expression_matrix"] == "Y")
    print(f"Wrote {out} with {len(rows)} accessions")
    print(f"  ICI-labeled: {n_ici}   KD/KO: {n_kd}   "
          f"embedded-expression matrices: {n_expr}")
    print("\nKD/KO sets:")
    for r in rows:
        if r["is_kd_ko"] == "Y":
            print(f"  {r['accession']:12s} expr={r['has_expression_matrix']} "
                  f"{r['title'][:60]}")


if __name__ == "__main__":
    main()
