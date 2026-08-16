#!/usr/bin/env python3
"""Deep-probe leftover 2020 lung+ICI series.

For each leftover candidate (plus a short list of borderline 2020 series that
the 2019–2021 wave parked as LUNG_ICI_other), download:
  * series-matrix header (characteristics, data table presence)
  * GEO SOFT / FTP file listing
and decide, honestly, whether TACSTD2/CLDN4 vs a per-patient ICI outcome
can be measured from open processed data < 2 GB.

Nothing is analyzed here. This step only writes an audit table.
"""
import csv
import gzip
import io
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results" / "w200" / "GEO_2020"
DL = RES / "downloads"
DL.mkdir(parents=True, exist_ok=True)

# Extra 2020 series the 2019–2021 wave parked as LUNG_ICI_other / borderline.
# Probe them even if the 2020-only search misses a keyword.
FORCE_PROBE = [
    "GSE124885",
    "GSE139327",
    "GSE141479",
    "GSE148944",
    "GSE150972",
    "GSE154286",
    "GSE159785",
    "GSE159787",
    "GSE99995",
    "GSE136961",  # already known unusable; re-document
    "GSE126044",  # already analyzed; confirm, do not re-claim
]

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

OUTCOME_RE = re.compile(
    r"respon|recist|pfs|os |overall survival|progression|dcb|ndb|"
    r"benefit|mpr|pcr|non-?respon|outcome|durval|nivo|pembro|"
    r"atezo|ipilimumab|immunotherap|checkpoint|anti-pd",
    re.I,
)
GENE_RE = re.compile(r"\b(TACSTD2|TROP-?2|CLDN4|claudin-?4)\b", re.I)


def gse_ftp_dir(acc):
    # GSE126044 -> GSE126nnn
    prefix = acc[:6] + "nnn" if len(acc) == 9 else acc[:7] + "nnn"
    if len(acc) == 8:  # GSE99995
        prefix = acc[:5] + "nnn"
    return f"{FTP}/{prefix}/{acc}"


def fetch(url, dest=None, timeout=90):
    req = urllib.request.Request(url, headers={"User-Agent": "w200-GEO-2020/1.0"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
            if dest:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
            return data
        except Exception as e:  # noqa: BLE001
            wait = 2 ** attempt
            print(f"  fetch retry {attempt} {url} ({e})", file=sys.stderr)
            time.sleep(wait)
    return None


def list_ftp(acc):
    base = gse_ftp_dir(acc)
    listing = {}
    for sub in ("matrix", "suppl"):
        url = f"{base}/{sub}/"
        data = fetch(url)
        listing[sub] = []
        if not data:
            continue
        text = data.decode("utf-8", errors="replace")
        for m in re.finditer(r'href="([^"]+)"', text):
            name = m.group(1)
            if name in ("../", "./") or name.endswith("/"):
                continue
            listing[sub].append(name)
    return listing, base


def parse_matrix_header(raw):
    if raw is None:
        return {}
    if raw[:2] == b"\x1f\x8b":
        text = gzip.decompress(raw).decode("utf-8", errors="replace")
    else:
        text = raw.decode("utf-8", errors="replace")
    lines = text.splitlines()
    chars = []
    n_samples = 0
    has_table = False
    table_genes = []
    in_table = False
    header = None
    for line in lines:
        if line.startswith("!Sample_geo_accession"):
            n_samples = max(n_samples, len(line.split("\t")) - 1)
        if line.startswith("!Sample_characteristics") or line.startswith("!Sample_source") or line.startswith("!Sample_title"):
            chars.append(line[:500])
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            has_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            if header is None:
                header = line
                continue
            table_genes.append(line.split("\t", 1)[0].strip('"'))
            if len(table_genes) > 40000:
                break
    joined = "\n".join(chars)
    gene_hits = sorted(set(GENE_RE.findall("\n".join(table_genes))))
    # also scan ID_REF style (may be probes, not symbols)
    return {
        "n_samples_matrix": n_samples,
        "has_embedded_table": has_table,
        "n_table_rows_scanned": len(table_genes),
        "embedded_gene_hits": gene_hits,
        "outcome_like_characteristics": bool(OUTCOME_RE.search(joined)),
        "characteristic_snippets": [c for c in chars if OUTCOME_RE.search(c)][:12],
        "all_characteristic_keys": sorted({
            c.split("\t", 1)[0] for c in chars
        }),
    }


def probe_one(acc, meta):
    ddir = DL / acc
    ddir.mkdir(parents=True, exist_ok=True)
    listing, base = list_ftp(acc)
    matrix_files = listing.get("matrix") or []
    suppl_files = listing.get("suppl") or []
    matrix_name = next((n for n in matrix_files if n.endswith("series_matrix.txt.gz")), None)
    header_info = {}
    if matrix_name:
        dest = ddir / matrix_name
        raw = dest.read_bytes() if dest.exists() else fetch(f"{base}/matrix/{matrix_name}", dest)
        header_info = parse_matrix_header(raw)

    # Scan first 2 MB of each processed-looking suppl file for gene symbols
    # only when the name looks like an expression matrix and is not huge.
    suppl_gene_hits = {}
    for name in suppl_files:
        low = name.lower()
        if not any(k in low for k in ("count", "tpm", "fpkm", "expr", "gene", "matrix", "txt", "csv", "tsv")):
            continue
        if any(k in low for k in ("fastq", "bam", "bed", "bw", "bigwig", "hic", "h5ad", "rds")):
            continue
        url = f"{base}/suppl/{name}"
        # HEAD-ish via range
        req = urllib.request.Request(url, headers={"User-Agent": "w200-GEO-2020/1.0", "Range": "bytes=0-1500000"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                chunk = r.read(1_500_000)
            if chunk[:2] == b"\x1f\x8b":
                try:
                    chunk = gzip.decompress(chunk)
                except Exception:
                    # partial gzip — try incremental
                    try:
                        d = gzip.GzipFile(fileobj=io.BytesIO(chunk))
                        chunk = d.read(800_000)
                    except Exception:
                        chunk = b""
            text = chunk.decode("utf-8", errors="ignore")
            hits = sorted(set(m.group(1).upper() for m in GENE_RE.finditer(text)))
            suppl_gene_hits[name] = hits
        except Exception as e:  # noqa: BLE001
            suppl_gene_hits[name] = [f"ERROR:{e}"]
        time.sleep(0.2)

    verdict, why = decide(acc, meta, header_info, suppl_files, suppl_gene_hits)
    rec = {
        "accession": acc,
        "pdat": meta.get("pdat"),
        "n_samples": meta.get("n_samples"),
        "title": meta.get("title"),
        "matrix_files": ";".join(matrix_files),
        "suppl_files": ";".join(suppl_files),
        "has_embedded_table": header_info.get("has_embedded_table"),
        "n_table_rows_scanned": header_info.get("n_table_rows_scanned"),
        "embedded_gene_hits": ";".join(header_info.get("embedded_gene_hits") or []),
        "outcome_like_characteristics": header_info.get("outcome_like_characteristics"),
        "characteristic_snippets": " | ".join(header_info.get("characteristic_snippets") or [])[:800],
        "suppl_gene_hits": json.dumps(suppl_gene_hits),
        "verdict": verdict,
        "why": why,
    }
    (ddir / "probe.json").write_text(json.dumps({
        **rec,
        "header_info": header_info,
        "listing": listing,
        "ftp_base": base,
    }, indent=2))
    return rec


def decide(acc, meta, header, suppl, suppl_hits):
    title = (meta.get("title") or "") + " " + (meta.get("summary") or "")
    if acc == "GSE126044":
        return "ALREADY_ANALYZED", "Prior 2019-2021 / ici-bulk waves already tested both genes vs anti-PD-1 response"
    if acc == "GSE136961":
        return "NO_TARGET_GENES", "Oncomine Immune Response 395-gene panel; TACSTD2/CLDN4 absent (reconfirmed by prior waves)"
    if re.search(r"SARS-CoV-2|COVID|HIV infection|oral tongue|glioblastoma", title, re.I):
        return "NOT_LUNG_ICI", "Title/summary is not a lung-cancer ICI treatment cohort"
    hits = set(header.get("embedded_gene_hits") or [])
    for v in suppl_hits.values():
        hits.update(x for x in v if not str(x).startswith("ERROR"))
    has_outcome = bool(header.get("outcome_like_characteristics"))
    has_both = ("TACSTD2" in {h.upper() for h in hits} or "TROP2" in {h.upper() for h in hits}) and (
        "CLDN4" in {h.upper() for h in hits}
    )
    has_either = bool(hits)
    if has_both and has_outcome:
        return "ANALYZE", "Open processed matrix appears to contain TACSTD2 and CLDN4 and GEO characteristics look like an ICI outcome"
    if has_either and has_outcome:
        return "MAYBE_ANALYZE", f"Outcome-like metadata plus partial gene hits {sorted(hits)}"
    if has_outcome and not has_either:
        return "OUTCOME_BUT_GENES_UNCONFIRMED", "Outcome-like characteristics; gene presence not confirmed in the first 1.5 MB of processed files"
    if has_either and not has_outcome:
        return "GENES_NO_OUTCOME", f"Genes seen ({sorted(hits)}) but no per-sample ICI outcome in GEO characteristics"
    return "NOT_ANALYZABLE", "No confirmed target genes and no per-sample ICI outcome in the open GEO record"


def main():
    recs = {r["accession"]: r for r in json.loads((RES / "candidates_metadata.json").read_text())}
    triage_path = RES / "tables" / "triage_2020.csv"
    leftover = []
    if triage_path.exists():
        with triage_path.open() as f:
            leftover = [row["accession"] for row in csv.DictReader(f)
                        if row["category"] == "LUNG_ICI_LEFTOVER_CANDIDATE"]
    todo = []
    for acc in leftover + FORCE_PROBE:
        if acc not in todo:
            todo.append(acc)

    # For FORCE_PROBE accessions missing from the 2020 search, fetch esummary by accession.
    missing = [a for a in todo if a not in recs]
    if missing:
        print("Fetching metadata for force-probe accessions not in 2020 search:", missing)
        for acc in missing:
            url = f"{EUTILS}/esearch.fcgi?db=gds&term={acc}[ACCN]+AND+gse[Entry+Type]&retmode=json"
            try:
                with urllib.request.urlopen(url, timeout=30) as r:
                    ids = json.load(r).get("esearchresult", {}).get("idlist", [])
            except Exception as e:  # noqa: BLE001
                print("  esearch fail", acc, e)
                continue
            if not ids:
                recs[acc] = {"accession": acc, "title": "", "summary": "", "pdat": "", "n_samples": ""}
                continue
            url = f"{EUTILS}/esummary.fcgi?db=gds&id={ids[0]}&retmode=json"
            with urllib.request.urlopen(url, timeout=30) as r:
                js = json.load(r)
            uid = js["result"]["uids"][0]
            recs[acc] = js["result"][uid]
            recs[acc]["accession"] = acc
            time.sleep(0.3)

    rows = []
    for acc in todo:
        print(f"== probing {acc} ==")
        rows.append(probe_one(acc, recs.get(acc, {"accession": acc})))
        time.sleep(0.25)

    out = RES / "tables" / "leftover_probe.csv"
    cols = list(rows[0].keys()) if rows else []
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print("Wrote", out)
    for r in rows:
        print(f"  {r['accession']:12s} {r['verdict']:32s} {r['why'][:100]}")


if __name__ == "__main__":
    main()
