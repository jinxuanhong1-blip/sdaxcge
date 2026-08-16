#!/usr/bin/env python3
"""Download the selected GEO series payloads and record a verifiable provenance manifest.

Every download is logged with its URL, HTTP Content-Length, on-disk size and SHA-256, plus the
accession metadata pulled straight from Entrez, so the accessions can be independently verified.
Files above 2 GB and FASTQ/SRA payloads are refused rather than fetched.
"""
from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"
DATA = Path("/tmp/geo_dl/data")
DATA.mkdir(parents=True, exist_ok=True)

FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
MAX_BYTES = 2 * 1024**3
REFUSE = re.compile(r"\.(fastq|fq|sra|bam|cram|bcl)(\.gz|\.bz2)?$", re.I)

# accession -> list of series-level payloads to fetch ("matrix" means the series matrix file(s))
MANIFEST: dict[str, list[str]] = {
    # Tier A: lung tumour, ICI-treated, response/outcome in GEO sample characteristics
    "GSE161537": ["GSE161537_nivobio_log2cpm.csv.gz", "matrix"],
    "GSE309652": ["GSE309652_RAW.tar", "matrix"],
    # Tier A comparator: same platform/publication, NSCLC NOT treated with ICI
    "GSE162520": ["GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz", "matrix"],
    # Tier B: lung tumour, ICI-treated, outcome must be derived from design/publication
    "GSE253564": ["GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz", "matrix"],
    "GSE248378": ["GSE248378_Durva_Post_FPKMs.txt.gz", "matrix"],
    "GSE182328": ["GSE182328_Gene_counts_matrix.txt.gz", "matrix"],
    "GSE248249": ["matrix"],
    # Tier C: non-tumour compartments with explicit response labels (epithelial-gene sanity checks)
    "GSE216297": ["matrix"],
    "GSE111414": ["GSE111414_gene_counts.csv.gz", "matrix"],
}


def stub(acc: str) -> str:
    n = acc[3:]
    return f"GSE{n[:-3]}nnn" if len(n) > 3 else "GSEnnn"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def http_get(url: str, dest: Path) -> tuple[bool, str, int | None]:
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "opus-geo-leftover/1.0"})
            with urllib.request.urlopen(req, timeout=600) as fh:
                declared = fh.headers.get("Content-Length")
                declared_n = int(declared) if declared else None
                if declared_n is not None and declared_n > MAX_BYTES:
                    return False, f"refused: {declared_n} bytes exceeds 2 GB cap", declared_n
                tmp = dest.with_suffix(dest.suffix + ".part")
                with tmp.open("wb") as out:
                    while True:
                        chunk = fh.read(1 << 20)
                        if not chunk:
                            break
                        out.write(chunk)
                        if out.tell() > MAX_BYTES:
                            tmp.unlink(missing_ok=True)
                            return False, "refused: stream exceeded 2 GB cap", declared_n
                tmp.rename(dest)
                return True, "ok", declared_n
        except urllib.error.HTTPError as exc:
            return False, f"http {exc.code}", None
        except Exception as exc:  # noqa: BLE001
            if attempt == 3:
                return False, f"error: {exc}", None
            time.sleep(2 ** attempt)
    return False, "exhausted retries", None


def verify_accession(acc: str) -> dict:
    """Independent accession check straight from Entrez, not from our own cached search results."""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmode=json"
        f"&term={urllib.parse.quote(acc + '[Accession]')}"
    )
    uids = json.loads(urllib.request.urlopen(url, timeout=120).read())["esearchresult"]["idlist"]
    if not uids:
        return {"accession": acc, "verified": False, "error": "no Entrez record"}
    uid = uids[0]
    surl = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&retmode=json"
        f"&id={uid}"
    )
    m = json.loads(urllib.request.urlopen(surl, timeout=120).read())["result"][uid]
    return {
        "accession": acc,
        "verified": m.get("accession") == acc,
        "entrez_accession": m.get("accession"),
        "entrez_uid": uid,
        "title": m.get("title"),
        "gdstype": m.get("gdstype"),
        "taxon": m.get("taxon"),
        "n_samples": m.get("n_samples"),
        "platform": m.get("gpl"),
        "pubmed_ids": m.get("pubmedids"),
        "release_date": m.get("pdat"),
        "bioproject": m.get("bioproject"),
        "ftp": m.get("ftplink"),
        "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
    }


def matrix_urls(acc: str) -> list[str]:
    base = f"{FTP_BASE}/{stub(acc)}/{acc}/matrix/"
    html = urllib.request.urlopen(base, timeout=120).read().decode("utf-8", "replace")
    names = sorted(set(re.findall(r'href="([^"]*series_matrix\.txt\.gz)"', html)))
    return [base + n for n in names]


def main() -> int:
    records = []
    for acc, wanted in MANIFEST.items():
        info = verify_accession(acc)
        time.sleep(0.4)
        acc_dir = DATA / acc
        acc_dir.mkdir(parents=True, exist_ok=True)
        files = []
        urls: list[str] = []
        for w in wanted:
            if w == "matrix":
                urls.extend(matrix_urls(acc))
            else:
                urls.append(f"{FTP_BASE}/{stub(acc)}/{acc}/suppl/{urllib.parse.quote(w)}")
        for url in urls:
            name = urllib.parse.unquote(url.rsplit("/", 1)[1])
            if REFUSE.search(name):
                files.append({"name": name, "url": url, "status": "skipped: raw sequencing payload"})
                continue
            dest = acc_dir / name
            if dest.exists():
                files.append(
                    {
                        "name": name,
                        "url": url,
                        "status": "cached",
                        "bytes": dest.stat().st_size,
                        "sha256": sha256(dest),
                    }
                )
                continue
            ok, msg, declared = http_get(url, dest)
            rec = {"name": name, "url": url, "status": msg, "declared_bytes": declared}
            if ok:
                rec["bytes"] = dest.stat().st_size
                rec["sha256"] = sha256(dest)
            files.append(rec)
            print(f"  {acc} {name}: {msg} ({rec.get('bytes')} bytes)")
        records.append({**info, "files": files})
        print(f"[dl] {acc} verified={info['verified']} n_samples={info.get('n_samples')}")

    (OUT / "download_manifest.json").write_text(json.dumps(records, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
