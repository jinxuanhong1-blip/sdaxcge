#!/usr/bin/env python3
"""Download public lung matrices and extract TACSTD2 + TGF-β genes.

Primary (tumors):
  UCSC Xena GDC hub STAR TPM (log2(TPM+1), GENCODE v36)
    TCGA-LUAD.star_tpm.tsv.gz
    TCGA-LUSC.star_tpm.tsv.gz
    gencode.v36.annotation.gtf.gene.probemap
  GDC open ABSOLUTE purity table
    https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5

Sensitivity (cell lines, optional):
  DepMap Public 24Q4 Figshare+ expression + Model.csv
  Streamed; only requested gene columns are written.

Raw matrices stay in a cache directory and are not committed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

GDC_HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PURITY_URL = "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
EXPR_URL = "https://ndownloader.figshare.com/files/51065489"
UA = "sdaxcge-w200-A11-TGFB/1.0 (public TCGA/DepMap extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"

# Current HGNC symbols (GENCODE v36 / DepMap 24Q4).
# Entrez used only to disambiguate DepMap "SYMBOL (entrez)" headers.
# HALLMARK_TGF_BETA_SIGNALING: MSigDB human set, 54 genes (CC-BY-4.0).
HALLMARK_TGFB = {
    "ACVR1": "90",
    "APC": "324",
    "ARID4B": "51742",
    "BCAR3": "8412",
    "BMP2": "650",
    "BMPR1A": "657",
    "BMPR2": "659",
    "CDH1": "999",
    "CDK9": "1025",
    "CDKN1C": "1028",
    "CTNNB1": "1499",
    "ENG": "2022",
    "FKBP1A": "2280",
    "FNTA": "2339",
    "FURIN": "5045",
    "HDAC1": "3065",
    "HIPK2": "28996",
    "ID1": "3397",
    "ID2": "3398",
    "ID3": "3399",
    "IFNGR2": "3460",
    "JUNB": "3726",
    "KLF10": "7071",
    "LEFTY2": "7044",
    "LTBP2": "4053",
    "MAP3K7": "6885",
    "NCOR2": "9612",
    "NOG": "9241",
    "PMEPA1": "56937",
    "PPM1A": "5494",
    "PPP1CA": "5499",
    "PPP1R15A": "23645",
    "RAB31": "11031",
    "RHOA": "387",
    "SERPINE1": "5054",
    "SKI": "6497",
    "SKIL": "6498",
    "SLC20A1": "6574",
    "SMAD1": "4086",
    "SMAD3": "4088",
    "SMAD6": "4091",
    "SMAD7": "4092",
    "SMURF1": "57154",
    "SMURF2": "64750",
    "SPTBN1": "6711",
    "TGFB1": "7040",
    "TGFBR1": "7046",
    "TGIF1": "7050",
    "THBS1": "7057",
    "TJP1": "7082",
    "TRIM33": "51592",
    "UBE2D3": "7323",
    "WWTR1": "25937",
    "XIAP": "331",
}

# Ligands / receptors / SMADs not all in HALLMARK.
CORE_EXTRA = {
    "TGFB2": "7042",
    "TGFB3": "7043",
    "TGFBR2": "7048",
    "TGFBR3": "7049",
    "SMAD2": "4087",
    "SMAD4": "4089",
}

# Mariathasan 2018 Nature fibroblast TGF-β response (F-TBRS) genes used as a
# published TGF-β *response* signature. Both current and GENCODE-v36 aliases
# are requested so the extract is not empty if the matrix still uses CTGF.
FTBRS = {
    "ACTA2": "59",
    "ACTG2": "72",
    "ADAM12": "8038",
    "ADAM19": "8728",
    "CNN1": "1264",
    "COL4A1": "1282",
    "CCN2": "1490",
    "CTGF": "1490",
    "CTPS1": "1503",
    "RFLNA": "283991",
    "FAM101B": "283991",
    "FSTL3": "10272",
    "HSPB1": "3315",
    "IGFBP3": "3486",
    "IL11": "3589",
    "NT5E": "4907",
    "OLFML2B": "91562",
    "PPP1R13L": "10848",
    "PXDC1": "222166",
    "SH3PXD2A": "9644",
    "TAGLN": "6876",
    "TGFBI": "7045",
    "TNS1": "7145",
    "TPM1": "7168",
}

CONTROLS = {
    "TACSTD2": "4070",
    "CLDN4": "1364",
    "CD8A": "925",
}

GENES = {**HALLMARK_TGFB, **CORE_EXTRA, **FTBRS, **CONTROLS}

TCGA_FILES = {
    "TCGA-LUAD.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-LUAD.star_tpm.tsv.gz",
    "TCGA-LUSC.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-LUSC.star_tpm.tsv.gz",
    "gencode.v36.probemap": f"{GDC_HUB}/gencode.v36.annotation.gtf.gene.probemap",
    "tcga_absolute_purity.txt": PURITY_URL,
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, timeout: int = 600) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest.name} already present ({dest.stat().st_size:,} bytes)", flush=True)
        return
    print(f"[get ] {url} -> {dest}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    print(f"[done] {dest.name}  size={dest.stat().st_size:,}  sha256={sha256_file(dest)[:16]}...", flush=True)


def extract_tcga_genes(tpm_path: Path, probemap_path: Path, out_path: Path) -> dict:
    """Write a compact samples x genes table (log2(TPM+1)) for requested symbols."""
    pm = {}
    with open(probemap_path, "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
        id_i = header.index("id")
        gene_i = header.index("gene")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            pm[parts[id_i]] = parts[gene_i]

    wanted = set(GENES)
    keep_idx: dict[str, int] = {}
    means: dict[str, float] = {}
    with gzip.open(tpm_path, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        for i, row in enumerate(reader):
            ens = row[0]
            sym = pm.get(ens)
            if sym not in wanted:
                continue
            vals = [float(x) if x not in ("", "NA") else float("nan") for x in row[1:]]
            mean = sum(v for v in vals if v == v) / max(1, sum(v == v for v in vals))
            if sym not in keep_idx or mean > means[sym]:
                keep_idx[sym] = i
                means[sym] = mean

    found: dict[str, list[float]] = {}
    with gzip.open(tpm_path, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for i, row in enumerate(reader):
            ens = row[0]
            sym = pm.get(ens)
            if sym in keep_idx and keep_idx[sym] == i:
                found[sym] = [float(x) if x not in ("", "NA") else float("nan") for x in row[1:]]

    missing = sorted(wanted - set(found))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["sample"] + sorted(found))
        for j, sample in enumerate(samples):
            w.writerow([sample] + [found[g][j] for g in sorted(found)])
    return {
        "n_samples": len(samples),
        "genes_found": sorted(found),
        "genes_missing": missing,
        "winning_row_mean": {k: float(v) for k, v in means.items()},
    }


def match_gene_column(header: list[str], symbol: str, entrez: str) -> str | None:
    exact = f"{symbol} ({entrez})"
    if exact in header:
        return exact
    candidates = [
        h
        for h in header
        if h == symbol or h.startswith(f"{symbol} (") or h.startswith(f"{symbol}(")
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        return None
    # Prefer the matching entrez if several aliases exist.
    for h in candidates:
        if f"({entrez})" in h:
            return h
    return candidates[0]


def extract_depmap(expr_path: Path, out_path: Path) -> dict:
    with expr_path.open("r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        cols: dict[str, str] = {}
        missing: list[str] = []
        for sym, entrez in GENES.items():
            col = match_gene_column(header, sym, entrez)
            if col is None:
                missing.append(f"{sym} ({entrez})")
            else:
                cols[sym] = col
        idx = {sym: header.index(col) for sym, col in cols.items()}
        n = 0
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["ModelID"] + list(cols))
            for row in reader:
                if not row:
                    continue
                w.writerow([row[0]] + [row[idx[s]] for s in cols])
                n += 1
    return {"n_models": n, "columns": cols, "genes_missing": missing}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/a11_tgfb_cache")
    p.add_argument("--out-dir", default="results/w200/A11_TGFB")
    p.add_argument("--skip-depmap", action="store_true")
    args = p.parse_args()

    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    manifest: dict = {
        "task": "A11_TGFB",
        "genes": GENES,
        "hallmark_tgfb": sorted(HALLMARK_TGFB),
        "ftbrs": sorted(FTBRS),
        "tcga": {},
        "depmap": None,
    }

    for name, url in TCGA_FILES.items():
        dest = cache / name
        if name == "tcga_absolute_purity.txt" and (out / name).exists() and (out / name).stat().st_size > 0:
            # Prefer a previously copied open GDC table if the API rejects HEAD/GET.
            if not dest.exists():
                dest.write_bytes((out / name).read_bytes())
                print(f"[copy] purity from {out / name}", flush=True)
            manifest["tcga"][name] = {
                "url": url,
                "bytes": dest.stat().st_size,
                "sha256": sha256_file(dest),
                "source": "local_copy_or_cache",
            }
            continue
        try:
            download(url, dest)
            manifest["tcga"][name] = {
                "url": url,
                "bytes": dest.stat().st_size,
                "sha256": sha256_file(dest),
            }
        except Exception as exc:
            if name == "tcga_absolute_purity.txt" and (out / name).exists():
                dest.write_bytes((out / name).read_bytes())
                manifest["tcga"][name] = {
                    "url": url,
                    "error": str(exc),
                    "fallback": "results/w200/A11_TGFB/tcga_absolute_purity.txt",
                    "bytes": dest.stat().st_size,
                    "sha256": sha256_file(dest),
                }
                print(f"[warn] GDC purity download failed ({exc}); used local copy.", flush=True)
            else:
                raise

    pm = cache / "gencode.v36.probemap"
    for cohort in ("LUAD", "LUSC"):
        src = cache / f"TCGA-{cohort}.star_tpm.tsv.gz"
        dest = out / f"tcga_{cohort.lower()}_genes.csv"
        meta = extract_tcga_genes(src, pm, dest)
        manifest["tcga"][f"extract_{cohort}"] = meta
        print(f"[extract] {cohort}: found={len(meta['genes_found'])} missing={meta['genes_missing']}", flush=True)

    purity_src = cache / "tcga_absolute_purity.txt"
    purity_dst = out / "tcga_absolute_purity.txt"
    if purity_src.exists() and (not purity_dst.exists() or purity_dst.stat().st_size != purity_src.stat().st_size):
        purity_dst.write_bytes(purity_src.read_bytes())

    if not args.skip_depmap:
        try:
            model_path = cache / "Model.csv"
            expr_path = cache / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
            download(MODEL_URL, model_path)
            download(EXPR_URL, expr_path, timeout=1800)
            extract_meta = extract_depmap(expr_path, out / "depmap24q4_tgfb_all_models.csv")
            manifest["depmap"] = {
                "release": "DepMap Public 24Q4",
                "citation": "DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. https://doi.org/10.25452/figshare.plus.27993248.v1",
                "extract": extract_meta,
                "Model.csv_sha256": sha256_file(model_path),
            }
        except Exception as exc:
            manifest["depmap"] = {"error": str(exc), "note": "DepMap extract failed; TCGA analysis still proceeds."}
            print(f"[warn] DepMap skipped: {exc}", flush=True)

    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("[done] manifest written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
