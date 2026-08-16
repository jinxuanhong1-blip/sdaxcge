#!/usr/bin/env python3
"""Download open LUAD matrices and extract A10 genes.

Sources (all public; raw matrices stay in --cache-dir, not the repo):
  * UCSC Xena GDC hub STAR TPM (TCGA-LUAD, log2(TPM+1), GENCODE v36)
  * PanCanAtlas ABSOLUTE purity (GDC open file)
  * CPTAC pan-cancer freeze v1.2 LUAD RNA + protein (open S3)
  * DepMap Public 24Q4 expression + Model.csv (Figshare+)

Only the extracted gene tables are written under results/w200/A10_NKX21/.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

UA = (
    "sdaxcge-w200-A10-NKX21/1.0 "
    "(reproducible public LUAD extract; "
    "+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

GDC_HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
CPTAC_BASE = (
    "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/"
    "data_freeze_v1.2_reorganized/LUAD"
)
MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
EXPR_URL = "https://ndownloader.figshare.com/files/51065489"

# Primary + exploratory context genes. Entrez used only for DepMap headers.
GENES = {
    "NKX2-1": {"ensembl": "ENSG00000136352", "entrez": "7080"},
    "TACSTD2": {"ensembl": "ENSG00000184292", "entrez": "4070"},
    "CLDN4": {"ensembl": "ENSG00000189143", "entrez": "1364"},
    "SFTPB": {"ensembl": "ENSG00000168878", "entrez": "6439"},
    "NAPSA": {"ensembl": "ENSG00000131400", "entrez": "9476"},
    "KRT5": {"ensembl": "ENSG00000186081", "entrez": "3852"},
    "CLDN7": {"ensembl": "ENSG00000181885", "entrez": "1366"},
}

TCGA_FILES = {
    "TCGA-LUAD.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-LUAD.star_tpm.tsv.gz",
    "gencode.v36.probemap": f"{GDC_HUB}/gencode.v36.annotation.gtf.gene.probemap",
    "tcga_absolute_purity.txt": (
        "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
    ),
}

CPTAC_FILES = [
    "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "LUAD_phenotype.txt",
    "LUAD_meta.txt",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    rec = {"url": url, "path": str(dest)}
    if dest.exists() and dest.stat().st_size > 0:
        rec["status"] = "cached"
        rec["bytes"] = dest.stat().st_size
        rec["sha256"] = sha256_file(dest)
        print(f"[skip] {dest.name} ({rec['bytes']:,} bytes)", flush=True)
        return rec
    print(f"[get ] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    rec["status"] = "downloaded"
    rec["bytes"] = dest.stat().st_size
    rec["sha256"] = sha256_file(dest)
    print(f"[done] {dest.name}  size={rec['bytes']:,}  sha256={rec['sha256'][:12]}", flush=True)
    return rec


def load_probemap(path: Path) -> dict[str, str]:
    """Map Ensembl (with or without version) -> HGNC symbol."""
    ens2sym: dict[str, str] = {}
    with path.open() as f:
        header = f.readline()
        # Xena probemap: id, gene, chrom, chromStart, chromEnd, strand
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            eid, sym = parts[0], parts[1]
            ens2sym[eid] = sym
            ens2sym[eid.split(".")[0]] = sym
    return ens2sym


def extract_tcga_genes(tpm_path: Path, probemap_path: Path, out_path: Path) -> dict:
    ens2sym = load_probemap(probemap_path)
    wanted_ens = {v["ensembl"] for v in GENES.values()}
    wanted_sym = set(GENES)
    found: dict[str, list[str]] = {}
    opener = gzip.open if tpm_path.suffix == ".gz" or tpm_path.name.endswith(".gz") else open
    with opener(tpm_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            eid = parts[0]
            bare = eid.split(".")[0]
            sym = ens2sym.get(eid) or ens2sym.get(bare)
            if bare in wanted_ens or (sym in wanted_sym):
                key = sym if sym in wanted_sym else next(
                    (s for s, meta in GENES.items() if meta["ensembl"] == bare),
                    sym or bare,
                )
                # keep highest-mean row if duplicates
                vals = parts[1:]
                mean = sum(float(x) if x not in ("", "NA") else 0.0 for x in vals) / max(len(vals), 1)
                prev = found.get(key)
                if prev is None or mean > prev[0]:
                    found[key] = [mean, eid, vals]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["sample_id"] + sorted(found))
        for i, sid in enumerate(samples):
            row = [sid]
            for g in sorted(found):
                row.append(found[g][2][i])
            w.writerow(row)

    return {
        "n_samples": len(samples),
        "genes_found": {g: found[g][1] for g in sorted(found)},
        "genes_missing": sorted(set(GENES) - set(found)),
    }


def _read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    # CPTAC freeze files are tab-separated; first column is gene or sample id.
    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    return header, rows


def extract_cptac_gene_matrix(src: Path, out_path: Path, layer: str) -> dict:
    header, rows = _read_table(src)
    # LinkedOmics/CPTAC gene matrices: first header cell is often empty or "gene";
    # remaining headers are sample IDs; first column of each row is the gene symbol.
    samples = header[1:]
    found: dict[str, list[str]] = {}
    aliases = {g.upper(): g for g in GENES}
    aliases["NKX2.1"] = "NKX2-1"
    aliases["TTF1"] = "NKX2-1"
    aliases["TTF-1"] = "NKX2-1"
    aliases["TITF1"] = "NKX2-1"
    ens2sym = {meta["ensembl"]: sym for sym, meta in GENES.items()}
    for row in rows:
        if not row:
            continue
        raw = row[0].strip()
        bare = raw.split(".")[0]
        key = aliases.get(raw.upper()) or ens2sym.get(bare)
        if key is None:
            continue
        vals = row[1:]
        mean = sum(float(x) if x not in ("", "NA", "NaN") else 0.0 for x in vals) / max(len(vals), 1)
        prev = found.get(key)
        if prev is None or mean > prev[0]:
            found[key] = [mean, raw, vals]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["sample_id"] + sorted(found))
        for i, sid in enumerate(samples):
            row = [sid]
            for g in sorted(found):
                row.append(found[g][2][i] if i < len(found[g][2]) else "")
            w.writerow(row)

    return {
        "layer": layer,
        "source": src.name,
        "n_samples": len(samples),
        "genes_found": {g: found[g][1] for g in sorted(found)},
        "genes_missing": sorted(set(GENES) - set(found)),
    }


def match_depmap_column(header: list[str], symbol: str, entrez: str) -> str | None:
    exact = f"{symbol} ({entrez})"
    if exact in header:
        return exact
    cands = [
        h
        for h in header
        if h == symbol or h.startswith(f"{symbol} (") or h.startswith(f"{symbol}(")
    ]
    if len(cands) == 1:
        return cands[0]
    return None


def extract_depmap_genes(expr_path: Path, out_path: Path) -> dict:
    with expr_path.open("r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        cols: dict[str, int] = {}
        missing = []
        for sym, meta in GENES.items():
            col = match_depmap_column(header, sym, meta["entrez"])
            if col is None:
                missing.append(sym)
            else:
                cols[sym] = header.index(col)
        n = 0
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["ModelID"] + sorted(cols))
            for row in reader:
                if not row:
                    continue
                w.writerow([row[0]] + [row[cols[g]] for g in sorted(cols)])
                n += 1
    return {
        "n_models": n,
        "columns": {g: header[i] for g, i in cols.items()},
        "genes_missing": missing,
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/a10_nkx21_data")
    p.add_argument("--out-dir", default="results/w200/A10_NKX21")
    args = p.parse_args()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    manifest: dict = {
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "user_agent": UA,
        "files": [],
        "extracts": {},
    }

    # --- TCGA ---
    for name, url in TCGA_FILES.items():
        rec = fetch(url, cache / name)
        rec["name"] = name
        manifest["files"].append(rec)
        if rec.get("status") == "failed":
            continue

    tpm = cache / "TCGA-LUAD.star_tpm.tsv.gz"
    pmap = cache / "gencode.v36.probemap"
    if tpm.exists() and pmap.exists():
        manifest["extracts"]["tcga_luad"] = extract_tcga_genes(
            tpm, pmap, out / "tcga_luad_genes.csv"
        )
        print("extracted TCGA-LUAD genes:", manifest["extracts"]["tcga_luad"], flush=True)

    purity_src = cache / "tcga_absolute_purity.txt"
    if purity_src.exists():
        dest = out / "tcga_absolute_purity.txt"
        dest.write_bytes(purity_src.read_bytes())

    # --- CPTAC ---
    for name in CPTAC_FILES:
        rec = fetch(f"{CPTAC_BASE}/{name}", cache / name)
        rec["name"] = name
        manifest["files"].append(rec)

    rna = cache / "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
    prot = cache / "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
    if rna.exists():
        manifest["extracts"]["cptac_rna"] = extract_cptac_gene_matrix(
            rna, out / "cptac_luad_rna_genes.csv", "rna"
        )
        print("extracted CPTAC RNA:", manifest["extracts"]["cptac_rna"], flush=True)
    if prot.exists():
        manifest["extracts"]["cptac_protein"] = extract_cptac_gene_matrix(
            prot, out / "cptac_luad_protein_genes.csv", "protein"
        )
        print("extracted CPTAC protein:", manifest["extracts"]["cptac_protein"], flush=True)
    # Phenotype/meta stay in the cache; they are not required for the primary test.

    # --- DepMap ---
    model_path = cache / "Model.csv"
    expr_path = cache / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    rec = fetch(MODEL_URL, model_path)
    rec["name"] = "Model.csv"
    manifest["files"].append(rec)
    rec = fetch(EXPR_URL, expr_path)
    rec["name"] = "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    manifest["files"].append(rec)
    if expr_path.exists():
        manifest["extracts"]["depmap"] = extract_depmap_genes(
            expr_path, out / "depmap24q4_genes.csv"
        )
        print("extracted DepMap genes:", manifest["extracts"]["depmap"], flush=True)

    man_path = out / "download_manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(man_path, flush=True)
    failed = [f for f in manifest["files"] if f.get("status") == "failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
