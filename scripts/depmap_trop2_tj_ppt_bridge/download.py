#!/usr/bin/env python3
"""Download public DepMap 24Q4 RNA + Gygi CCLE protein for the PPT TJ bridge.

Keeps only pre-specified genes, then deletes full matrices.
Does not download private data.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

UA = "sdaxcge-depmap-trop2-tj-ppt-bridge/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"

MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
EXPR_URL = "https://ndownloader.figshare.com/files/51065489"
GMT_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
    "h.all.v2024.1.Hs.symbols.gmt"
)
GYGI_URLS = [
    "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
    "https://gygi.med.harvard.edu/sites/gygi.med.harvard.edu/files/documents/protein_quant_current_normalized.csv.gz",
]
RPPA_AB_URL = "https://data.broadinstitute.org/ccle/CCLE_RPPA_Ab_info_20180123.csv"

RELEASE = "DepMap Public 24Q4"
DOI = "10.25452/figshare.plus.27993248.v1"
HALLMARK_IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE"

# Humanized locked TJ lists (GSE137244 transfer signatures / TISMO).
TJ_EPITHELIAL = [
    "CLDN1",
    "CLDN3",
    "CLDN4",
    "CLDN7",
    "OCLN",
    "MARVELD2",
    "MARVELD3",
    "TJP1",
    "TJP2",
    "TJP3",
    "F11R",
    "JAM2",
    "JAM3",
    "CGN",
    "CGNL1",
    "CRB3",
    "ILDR1",
    "LSR",
]
TJ_TISMO = ["CLDN3", "CLDN4", "CLDN6", "CLDN7", "CDH1", "F11R", "OCLN"]
CLDN4_TJ_EDGE = [
    "CLDN4",
    "CLDN1",
    "CLDN7",
    "CGNL1",
    "MARVELD2",
    "MARVELD3",
    "TJP1",
    "TJP2",
    "ILDR1",
    "CLDN3",
    "OCLN",
]

# Same immune genes as the DepMap TROP2–IFN partial PR (for PPT continuity).
IMMUNE_LIGAND = ["CD274", "PDCD1LG2", "CXCL9", "CXCL10", "CXCL11", "HLA-E", "IDO1"]
MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
ISG = [
    "STAT1",
    "STAT2",
    "IRF1",
    "IRF9",
    "MX1",
    "ISG15",
    "IFIT1",
    "IFIT3",
    "OAS1",
    "OAS2",
    "OAS3",
    "EIF2AK2",
    "IFI35",
    "BST2",
    "SAMHD1",
    "TRIM25",
    "ADAR",
    "TAP1",
    "TAP2",
    "PSMB8",
    "PSMB9",
    "PSMB10",
]
CONTROLS = ["TACSTD2", "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19"]

CORE_GENES = sorted(
    set(CONTROLS)
    | set(TJ_EPITHELIAL)
    | set(TJ_TISMO)
    | set(CLDN4_TJ_EDGE)
    | set(IMMUNE_LIGAND)
    | set(MHC1)
    | set(ISG)
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, min_bytes: int) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    h = hashlib.sha256()
    n = 0
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
            h.update(chunk)
            n += len(chunk)
    if n < min_bytes:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"{url} wrote only {n} bytes")
    tmp.replace(dest)
    print(f"  {dest.name} {n} bytes sha256={h.hexdigest()[:12]}", flush=True)
    return h.hexdigest()


def fetch(urls: list[str], dest: Path, min_bytes: int) -> tuple[str, str]:
    last: Exception | None = None
    for url in urls:
        try:
            digest = download(url, dest, min_bytes)
            return url, digest
        except Exception as exc:
            last = exc
            print(f"WARN {url}: {exc}", flush=True)
            dest.unlink(missing_ok=True)
            time.sleep(2)
    raise SystemExit(f"failed to download {dest.name}: {last}")


def symbol_of(col: str) -> str:
    return col.split(" (")[0].strip()


def parse_gmt_ifng(path: Path) -> list[str]:
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if parts and parts[0] == HALLMARK_IFNG and len(parts) >= 3:
            return [g for g in parts[2:] if g]
    raise SystemExit(f"{HALLMARK_IFNG} not found in {path}")


def extract_model_by_gene(src: Path, wanted: set[str], dest: Path) -> dict:
    with src.open("r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        found: dict[str, int] = {}
        for i, col in enumerate(header[1:], start=1):
            sym = symbol_of(col)
            if sym in wanted and sym not in found:
                found[sym] = i
        keep = [0] + [found[s] for s in sorted(found)]
        missing = sorted(wanted - set(found))
        print(
            f"extract {src.name}: kept {len(found)}/{len(wanted)}; missing={missing}",
            flush=True,
        )
        dest.parent.mkdir(parents=True, exist_ok=True)
        n_rows = 0
        with dest.open("w", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(["ModelID"] + [symbol_of(header[i]) for i in keep[1:]])
            for row in reader:
                if not row:
                    continue
                writer.writerow([row[i] if i < len(row) else "" for i in keep])
                n_rows += 1
    return {
        "n_models": n_rows,
        "n_genes_kept": len(found),
        "n_genes_wanted": len(wanted),
        "genes_kept": sorted(found),
        "missing_genes": missing,
    }


def extract_gygi(src: Path, wanted: set[str], dest: Path) -> dict:
    opener = gzip.open if str(src).endswith(".gz") else open
    with opener(src, "rt", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        if "Gene_Symbol" not in header:
            raise SystemExit(f"Gygi header has no Gene_Symbol: {header[:8]}")
        sym_i = header.index("Gene_Symbol")
        id_i = header.index("Protein_Id") if "Protein_Id" in header else None
        meta = {"Protein_Id", "Gene_Symbol", "Description", "Group_ID", "Uniprot", "Uniprot_Acc"}
        sample_idx = [
            i
            for i, c in enumerate(header)
            if c not in meta and not c.endswith("_Peptides") and not c.startswith("TenPx")
        ]
        quant_idx = [i for i in sample_idx if "_TenPx" in header[i]]
        kept_rows = 0
        genes = []
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("w", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(["Gene_Symbol", "Protein_Id"] + [header[i] for i in quant_idx])
            for row in reader:
                if sym_i >= len(row):
                    continue
                gene = row[sym_i].split(" ")[0].strip()
                if gene not in wanted:
                    continue
                pid = row[id_i] if id_i is not None and id_i < len(row) else ""
                writer.writerow([gene, pid] + [row[i] if i < len(row) else "" for i in quant_idx])
                kept_rows += 1
                genes.append(gene)
    missing = sorted(wanted - set(genes))
    print(
        f"extract Gygi: rows {kept_rows}, genes {len(set(genes))}/{len(wanted)}; missing n={len(missing)}",
        flush=True,
    )
    return {
        "n_quant_columns": len(quant_idx),
        "n_rows_kept": kept_rows,
        "genes_kept": sorted(set(genes)),
        "missing_genes": missing,
        "n_lung_tenpx_columns": sum(1 for i in quant_idx if "_LUNG_" in header[i]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/depmap_trop2_tj_ppt_bridge")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gene_sets = {
        "TJ_EPITHELIAL": TJ_EPITHELIAL,
        "TJ_TISMO": TJ_TISMO,
        "CLDN4_TJ_EDGE": CLDN4_TJ_EDGE,
        "IMMUNE_LIGAND": IMMUNE_LIGAND,
        "MHC1": MHC1,
        "ISG": ISG,
        "CONTROLS": CONTROLS,
    }
    (outdir / "gene_sets.json").write_text(json.dumps(gene_sets, indent=2) + "\n")

    manifest: dict = {"release": RELEASE, "doi": DOI, "files": {}, "core_genes": CORE_GENES}

    model = outdir / "Model.csv"
    if not model.exists() or model.stat().st_size < 100_000:
        url, digest = fetch([MODEL_URL], model, 100_000)
    else:
        url, digest = MODEL_URL, sha256_file(model)
        print(f"OK exists {model}", flush=True)
    manifest["files"]["Model.csv"] = {"url": url, "sha256": digest, "bytes": model.stat().st_size}

    gmt = outdir / "h.all.v2024.1.Hs.symbols.gmt"
    if not gmt.exists() or gmt.stat().st_size < 5_000:
        url, digest = fetch([GMT_URL], gmt, 5_000)
    else:
        url, digest = GMT_URL, sha256_file(gmt)
        print(f"OK exists {gmt}", flush=True)
    ifng = parse_gmt_ifng(gmt)
    (outdir / "hallmark_ifng_genes.txt").write_text("\n".join(ifng) + "\n")
    manifest["files"]["gmt"] = {
        "url": url,
        "sha256": digest,
        "bytes": gmt.stat().st_size,
        "n_hallmark_ifng": len(ifng),
    }

    wanted_expr = set(CORE_GENES) | set(ifng)
    expr_full = outdir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    expr_slim = outdir / "expression_panel.csv"
    if not expr_slim.exists() or expr_slim.stat().st_size < 50_000:
        if not expr_full.exists() or expr_full.stat().st_size < 400_000_000:
            url, digest = fetch([EXPR_URL], expr_full, 400_000_000)
        else:
            url, digest = EXPR_URL, sha256_file(expr_full)
        meta = extract_model_by_gene(expr_full, wanted_expr, expr_slim)
        expr_full.unlink(missing_ok=True)
    else:
        url, digest = EXPR_URL, None
        with expr_slim.open() as fh:
            header = next(csv.reader(fh))
        meta = {
            "n_models": sum(1 for _ in expr_slim.open()) - 1,
            "genes_kept": header[1:],
            "missing_genes": sorted(wanted_expr - set(header[1:])),
            "note": "slim extract already present",
        }
        print(f"OK exists {expr_slim}", flush=True)
    manifest["files"]["expression"] = {"url": url, "sha256": digest, "extract": meta}

    gygi = outdir / "protein_quant_current_normalized.csv.gz"
    gygi_slim = outdir / "gygi_protein_panel.csv"
    if not gygi_slim.exists() or gygi_slim.stat().st_size < 10_000:
        if not gygi.exists() or gygi.stat().st_size < 5_000_000:
            url, digest = fetch(GYGI_URLS, gygi, 5_000_000)
        else:
            url, digest = GYGI_URLS[0], sha256_file(gygi)
        meta = extract_gygi(gygi, wanted_expr | {"TACSTD2", "CLDN4"}, gygi_slim)
        gygi.unlink(missing_ok=True)
    else:
        url, digest = GYGI_URLS[0], None
        meta = {"note": "slim extract already present"}
        print(f"OK exists {gygi_slim}", flush=True)
    manifest["files"]["gygi"] = {"url": url, "sha256": digest, "extract": meta}

    ab = outdir / "CCLE_RPPA_Ab_info_20180123.csv"
    try:
        if not ab.exists() or ab.stat().st_size < 1_000:
            url, digest = fetch([RPPA_AB_URL], ab, 1_000)
        else:
            url, digest = RPPA_AB_URL, sha256_file(ab)
        text = ab.read_text(errors="replace")
        manifest["files"]["rppa_antibody_info"] = {
            "url": url,
            "sha256": digest,
            "bytes": ab.stat().st_size,
            "mentions_CLDN4": "CLDN4" in text or "Claudin-4" in text or "Claudin 4" in text,
            "mentions_TACSTD2_or_TROP2": "TACSTD2" in text or "TROP2" in text,
            "mentions_Claudin7": "Claudin-7" in text or "Claudin 7" in text,
        }
    except SystemExit as exc:
        manifest["files"]["rppa_antibody_info"] = {"error": str(exc)}
        print(f"RPPA antibody info unavailable: {exc}", flush=True)

    (outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("manifest written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
