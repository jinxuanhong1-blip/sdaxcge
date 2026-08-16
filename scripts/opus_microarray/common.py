"""Shared helpers for the opus_microarray slice.

Scope of this slice (do not write outside these trees):
  notes/opus_microarray/, scripts/opus_microarray/, results/opus_microarray/

Downloaded GEO payloads are cached under results/opus_microarray/cache/ (git-ignored
for the bulky files) and every download is recorded with size + sha256 in
results/opus_microarray/download_manifest.tsv so the analysis is auditable.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SLICE = "opus_microarray"
NOTES_DIR = REPO_ROOT / "notes" / SLICE
SCRIPTS_DIR = REPO_ROOT / "scripts" / SLICE
RESULTS_DIR = REPO_ROOT / "results" / SLICE
CACHE_DIR = RESULTS_DIR / "cache"
MANIFEST = RESULTS_DIR / "download_manifest.tsv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
USER_AGENT = "opus-microarray-slice/1.0 (GEO reanalysis; contact: repo maintainer)"

# Genes of interest plus positive controls used to sanity-check each dataset.
TARGET_GENES = ["TACSTD2", "CLDN4"]
CONTROL_GENES = ["CD274", "PDCD1", "CD8A", "GZMB", "IFNG", "CXCL9", "STAT1", "HLA-DRA"]

# Aliases seen in GEO platform annotation for the target genes.
GENE_ALIASES = {
    "TACSTD2": {"TACSTD2", "TROP2", "TROP-2", "GA733-1", "M1S1", "EGP-1", "GP50", "TACSTD-2"},
    "CLDN4": {"CLDN4", "CPE-R", "CPER", "CPETR1", "WBSCR8", "hCPE-R", "CPETR"},
    "CD274": {"CD274", "PD-L1", "PDL1", "B7-H1", "B7H1", "PDCD1L1", "PDCD1LG1"},
    "PDCD1": {"PDCD1", "PD-1", "PD1", "CD279", "SLEB2", "hSLE1"},
    "CD8A": {"CD8A", "CD8", "MAL", "p32", "Leu2"},
    "GZMB": {"GZMB", "CTLA1", "CSPB", "CGL1", "HLP", "SECT", "CGL-1"},
    "IFNG": {"IFNG", "IFG", "IFI", "IFN-gamma", "IFNgamma"},
    "CXCL9": {"CXCL9", "MIG", "SCYB9", "CMK", "crg-10"},
    "STAT1": {"STAT1", "ISGF-3", "STAT91"},
    "HLA-DRA": {"HLA-DRA", "HLA-DRA1", "HLADRA", "MLRW"},
}

DEFAULT_SLEEP = 0.4  # be polite to NCBI (<3 req/s without an API key)


def ensure_dirs() -> None:
    for d in (NOTES_DIR, SCRIPTS_DIR, RESULTS_DIR, CACHE_DIR):
        d.mkdir(parents=True, exist_ok=True)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _record_manifest(url: str, path: Path) -> None:
    ensure_dirs()
    rel = path.relative_to(RESULTS_DIR)
    line = f"{rel}\t{path.stat().st_size}\t{_sha256(path)}\t{url}\n"
    existing = MANIFEST.read_text().splitlines(keepends=True) if MANIFEST.exists() else []
    if not existing:
        existing = ["path\tbytes\tsha256\tsource_url\n"]
    kept = [l for l in existing if not l.startswith(f"{rel}\t")]
    kept.append(line)
    header, rows = kept[0], sorted(set(kept[1:]))
    MANIFEST.write_text(header + "".join(rows))


def fetch(url: str, dest: Path, *, retries: int = 5, force: bool = False) -> Path:
    """Download `url` to `dest` (cached). Retries with exponential backoff."""
    ensure_dirs()
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return dest
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=180) as resp:
                data = resp.read()
            if not data:
                raise RuntimeError("empty response")
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            _record_manifest(url, dest)
            time.sleep(DEFAULT_SLEEP)
            return dest
        except Exception as exc:  # noqa: BLE001 - network layer, want broad retry
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (404, 403):
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to download {url}: {last}")


def eutils(endpoint: str, params: dict[str, str], dest: Path, *, force: bool = False) -> Path:
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    return fetch(url, dest, force=force)


def gse_matrix_url(gse: str, filename: str | None = None) -> str:
    stub = gse[:-3] + "nnn"
    fn = filename or f"{gse}_series_matrix.txt.gz"
    return f"{GEO_FTP}/series/{stub}/{gse}/matrix/{fn}"


def gse_soft_url(gse: str) -> str:
    stub = gse[:-3] + "nnn"
    return f"{GEO_FTP}/series/{stub}/{gse}/soft/{gse}_family.soft.gz"


def gpl_annot_url(gpl: str) -> str:
    stub = gpl[:-3] + "nnn" if len(gpl) > 6 else "GPLnnn"
    return f"{GEO_FTP}/platforms/{stub}/{gpl}/annot/{gpl}.annot.gz"


def gpl_soft_url(gpl: str) -> str:
    stub = gpl[:-3] + "nnn" if len(gpl) > 6 else "GPLnnn"
    return f"{GEO_FTP}/platforms/{stub}/{gpl}/soft/{gpl}_family.soft.gz"


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + "\n")


def norm_gene(token: str) -> str:
    return token.strip().strip('"').upper().replace(" ", "")


def alias_set(gene: str) -> set[str]:
    return {norm_gene(a) for a in GENE_ALIASES.get(gene, {gene})}


def env_flag(name: str, default: bool = False) -> bool:
    v = os.environ.get(name)
    if v is None:
        return default
    return v.lower() in {"1", "true", "yes", "on"}
