#!/usr/bin/env python3
"""No-skip inventory of leftover 2019-2021 GEO lung ICI series.

Lists supplementary files (name + size) for every leftover series that was
NOT already analyzed in the first pass. We do not skip for size or tissue
(blood/PBMC). Processed expression matrices are flagged for download.
"""
import csv
import json
import re
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FABLE = ROOT / "results" / "fable_geo_2019_2021"
OUT = ROOT / "results" / "noskip" / "GEO_2019_2021"
OUT.mkdir(parents=True, exist_ok=True)

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

ALREADY = {"GSE126044", "GSE135222", "GSE182328", "GSE111414", "GSE136961"}

# Leftover = lung+ICI categories from first-pass triage (not already analyzed).
# We also keep EXCLUDED_not_ICI / EXCLUDED_not_lung out of the leftover compute
# set, but we still list LUNG_ICI_* leftovers exhaustively.
KEEP_CATS = {
    "LUNG_ICI_other",
    "LUNG_ICI_single_cell",
    "LUNG_ICI_methylation",
    "LUNG_ICI_cellline_invitro",
    "VERIFIED_NO_TARGET_GENES",
}

PROCESSED_RE = re.compile(
    r"(count|tpm|fpkm|rpkm|cpm|expression|expr|matrix|norm|log2|"
    r"normalized|quant|gene.?exp|htseq|featurecount|salmon|kallisto|"
    r"rsem|deseq|edger)",
    re.I,
)
SKIP_NAME_RE = re.compile(
    r"(RAW\.tar|filelist\.txt|\.bam|\.fastq|\.fq\.|\.bai|\.cram|"
    r"methyl|bedgraph|bigwig|bw\.|wig\.|peaks|narrowPeak|broadPeak|"
    r"vcf|maf\.|hic|cool|h5ad|loom|mtx\.|barcodes|features\.tsv|"
    r"filtered_feature|raw_feature|possorted|cloupe|web_summary)",
    re.I,
)


def nnn(gse):
    num = gse[3:]
    return f"GSE{num[:-3]}nnn" if len(num) > 3 else "GSEnnn"


def size_to_bytes(s):
    s = (s or "").strip()
    m = re.match(r"([\d.]+)([KMGT]?)", s)
    if not m:
        return None
    val = float(m.group(1))
    mult = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}[m.group(2)]
    return int(val * mult)


def fetch(url, timeout=60):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:  # noqa: BLE001
            time.sleep(2 ** attempt)
    return None


def list_suppl(gse):
    url = f"{FTP}/{nnn(gse)}/{gse}/suppl/"
    html = fetch(url)
    if html is None:
        return []
    files = []
    for m in re.finditer(
        r'<a href="([^"?/][^"]*)">[^<]+</a>\s+[^\n]*?\s+([\d.]+[KMGT]?|-)\s',
        html,
    ):
        name, size = m.group(1), m.group(2)
        if name.startswith("http") or name in ("Parent Directory",):
            continue
        files.append({
            "name": name,
            "size_str": None if size == "-" else size,
            "size_bytes": size_to_bytes(size) if size != "-" else None,
        })
    if not files:
        for m in re.finditer(r'<a href="([^"?/][^"]*)">', html):
            name = m.group(1)
            if not name.startswith("http"):
                files.append({"name": name, "size_str": None, "size_bytes": None})
    return files


def is_processed_candidate(name):
    if SKIP_NAME_RE.search(name):
        return False
    return bool(PROCESSED_RE.search(name))


def main():
    recs = {r["accession"]: r for r in json.loads((FABLE / "candidates_metadata.json").read_text())}
    triage = list(csv.DictReader(open(FABLE / "tables" / "triage_all_candidates.csv")))
    leftover = [row for row in triage if row["category"] in KEEP_CATS and row["accession"] not in ALREADY]

    inv = {}
    rows = []
    for i, row in enumerate(leftover, 1):
        acc = row["accession"]
        files = list_suppl(acc)
        processed = [f for f in files if is_processed_candidate(f["name"])]
        inv[acc] = {
            "taxon": recs[acc].get("taxon"),
            "title": recs[acc].get("title"),
            "n_samples": recs[acc].get("n_samples"),
            "category": row["category"],
            "reason": row["reason"],
            "files": files,
            "processed_candidates": processed,
        }
        print(f"[{i}/{len(leftover)}] {acc} files={len(files)} processed={len(processed)} "
              f"{[f['name'] for f in processed][:4]}")
        rows.append({
            "accession": acc,
            "category": row["category"],
            "n_samples": recs[acc].get("n_samples"),
            "n_suppl": len(files),
            "n_processed_candidates": len(processed),
            "processed_names": ";".join(f["name"] for f in processed),
            "processed_sizes": ";".join(str(f.get("size_str") or "") for f in processed),
            "title": recs[acc].get("title"),
        })
        time.sleep(0.12)

    (OUT / "leftover_suppl_inventory.json").write_text(json.dumps(inv, indent=2))
    with open(OUT / "leftover_inventory.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["accession"])
        w.writeheader()
        w.writerows(rows)
    print("leftover series:", len(leftover))
    print("with processed candidates:", sum(1 for r in rows if r["n_processed_candidates"] > 0))
    print("Wrote", OUT / "leftover_suppl_inventory.json")


if __name__ == "__main__":
    main()
