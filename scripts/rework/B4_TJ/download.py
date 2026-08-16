#!/usr/bin/env python3
"""Download public GSE126044 / GSE135222 files and documented TJ gene sets."""

from __future__ import annotations

import hashlib
import json
import re
import sys
import urllib.request
from pathlib import Path

from gene_sets import (
    KEGG_REST_URL,
    REACTOME_TJ_GENES,
    REACTOME_TJ_MSIGDB_URL,
    SYMBOL_FALLBACKS,
)

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "rework" / "B4_TJ"
DATA.mkdir(parents=True, exist_ok=True)

URLS = {
    "GSE126044_counts.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/"
        "suppl/GSE126044_counts.txt.gz"
    ),
    "GSE126044_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/"
        "matrix/GSE126044_series_matrix.txt.gz"
    ),
    "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/"
        "suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
    ),
    "GSE135222_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/"
        "matrix/GSE135222_series_matrix.txt.gz"
    ),
    "kegg_hsa04530.txt": KEGG_REST_URL,
}

MSIGDB_JSON = (
    "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp"
    "?geneSetName=REACTOME_TIGHT_JUNCTION_INTERACTIONS&fileType=json"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path, timeout: int = 120) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "B4_TJ_rework/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r, dest.open("wb") as out:
        out.write(r.read())


def parse_kegg_genes(text: str) -> list[str]:
    genes: list[str] = []
    in_gene = False
    for line in text.splitlines():
        if line.startswith("GENE"):
            in_gene = True
            rest = line[4:].strip()
        elif in_gene and (line.startswith(" ") or line.startswith("\t")):
            rest = line.strip()
        elif in_gene:
            break
        else:
            continue
        m = re.match(r"(\d+)\s+([A-Za-z0-9-]+)", rest)
        if m:
            genes.append(m.group(2))
    # unique, keep KEGG order
    seen: set[str] = set()
    out: list[str] = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def mygene_map(symbols: list[str]) -> list[dict]:
    body = json.dumps(
        {
            "q": symbols,
            "scopes": "symbol,alias",
            "fields": "ensembl.gene,symbol,alias,entrezgene",
            "species": "human",
            "size": 1,
        }
    ).encode()
    req = urllib.request.Request(
        "https://mygene.info/v3/query",
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": "B4_TJ_rework/1.0"},
    )
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)


def ensembl_from_hit(hit: dict) -> str | None:
    ens = hit.get("ensembl")
    if isinstance(ens, dict):
        return ens.get("gene")
    if isinstance(ens, list) and ens:
        g = ens[0].get("gene") if isinstance(ens[0], dict) else None
        return g
    return None


def main() -> int:
    manifest: dict = {"files": {}, "urls": URLS}
    for name, url in URLS.items():
        dest = DATA / name
        print(f"GET {url} -> {dest}", flush=True)
        fetch(url, dest)
        manifest["files"][name] = {"url": url, "sha256": sha256(dest), "bytes": dest.stat().st_size}

    msig_path = DATA / "REACTOME_TIGHT_JUNCTION_INTERACTIONS.json"
    print(f"GET {MSIGDB_JSON}", flush=True)
    try:
        fetch(MSIGDB_JSON, msig_path)
        msig = json.loads(msig_path.read_text())
        key = "REACTOME_TIGHT_JUNCTION_INTERACTIONS"
        remote = sorted(msig[key]["geneSymbols"])
        local = sorted(REACTOME_TJ_GENES)
        if remote != local:
            print("WARNING: MSigDB gene list differs from frozen REACTOME_TJ_GENES", file=sys.stderr)
            print(" remote", remote, file=sys.stderr)
            print(" frozen", local, file=sys.stderr)
        manifest["msigdb_reactome_tj"] = {
            "url": REACTOME_TJ_MSIGDB_URL,
            "download": MSIGDB_JSON,
            "sha256": sha256(msig_path),
            "n_genes": len(set(remote)),
            "matches_frozen": remote == local,
        }
    except Exception as exc:
        print(f"MSigDB JSON download failed ({exc}); using frozen list only", file=sys.stderr)
        manifest["msigdb_reactome_tj"] = {"download_failed": str(exc), "used_frozen": True}

    kegg_text = (DATA / "kegg_hsa04530.txt").read_text()
    kegg_genes = parse_kegg_genes(kegg_text)
    (DATA / "kegg_hsa04530_genes.txt").write_text("\n".join(kegg_genes) + "\n")
    manifest["kegg_hsa04530_n_genes"] = len(kegg_genes)

    query_symbols = sorted(
        set(REACTOME_TJ_GENES)
        | set(kegg_genes)
        | {"CLDN4", "CD8A", "CD8B", "TACSTD2", "INADL", "MPP5", "PATJ", "PALS1"}
        | {a for alts in SYMBOL_FALLBACKS.values() for a in alts}
    )
    print(f"mygene.info map for {len(query_symbols)} symbols", flush=True)
    hits = mygene_map(query_symbols)
    (DATA / "mygene_symbol_to_ensembl.json").write_text(json.dumps(hits, indent=2))
    mapping = {}
    for hit in hits:
        q = hit.get("query")
        if not q or hit.get("notfound"):
            continue
        mapping[q] = {
            "symbol": hit.get("symbol"),
            "ensembl": ensembl_from_hit(hit),
            "entrezgene": hit.get("entrezgene"),
        }
    (DATA / "symbol_to_ensembl.json").write_text(json.dumps(mapping, indent=2, sort_keys=True))
    manifest["n_symbol_mapped"] = len(mapping)
    (DATA / "download_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({k: manifest[k] for k in ("kegg_hsa04530_n_genes", "n_symbol_mapped")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
