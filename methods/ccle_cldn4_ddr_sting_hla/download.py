#!/usr/bin/env python3
"""Download DepMap 24Q4 RNA, MSigDB Hallmark GMT, and MCLP/CCLE RPPA proteins.

γH2AX (H2AX-pS139) is not in the public CCLE RPPA or MCLP antibody lists.
This script records that inventory and pulls the DNA-damage phospho marks
that those panels do measure, plus the RNA genes for DNA-repair, STING, and HLA.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

RELEASE = "DepMap Public 24Q4"
DOI = "10.25452/figshare.plus.27993248.v1"
MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
EXPR_URL = "https://ndownloader.figshare.com/files/51065489"
GMT_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
    "h.all.v2024.1.Hs.symbols.gmt"
)
RPPA_URLS = {
    "CCLE_RPPA_20180123": "https://data.broadinstitute.org/ccle/CCLE_RPPA_20180123.csv",
    "CCLE_RPPA_20181003": "https://data.broadinstitute.org/ccle/CCLE_RPPA_20181003.csv",
}
MCLP_BASE = "https://tcpa.drbioright.org/rppa500mclp"
MCLP_RELEASE = "20221116"
UA = "sdaxcge-ccle-cldn4-ddr-sting-hla/1.0 (public DepMap/MCLP extract)"

# Focused DSB/DDR response. H2AX here is the transcript, not γH2AX protein.
DSB_GENES = [
    "ATM",
    "ATR",
    "CHEK1",
    "CHEK2",
    "H2AX",
    "H2AFX",
    "BRCA1",
    "BRCA2",
    "RAD51",
    "MDC1",
    "TP53BP1",
    "PARP1",
    "PRKDC",
    "XRCC5",
    "XRCC6",
    "RPA1",
    "RPA2",
    "MRE11",
    "RAD50",
    "NBN",
    "RAD17",
]
STING_CORE = ["CGAS", "MB21D1", "STING1", "TMEM173", "TBK1", "IKBKE", "IRF3"]
STING_ISG = ["CXCL10", "CCL5", "ISG15", "IFIT1", "MX1", "OAS1", "IRF7", "IFNB1"]
HLA_GENES = [
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-F",
    "HLA-G",
    "B2M",
    "HLA-DRA",
    "HLA-DRB1",
    "HLA-DPA1",
    "HLA-DPB1",
    "HLA-DQA1",
    "HLA-DQB1",
]
CORE = ["CLDN4", "TACSTD2", "CD274"]
HALLMARKS = ["HALLMARK_DNA_REPAIR", "HALLMARK_INTERFERON_GAMMA_RESPONSE"]

# MCLP protein_id values. validation is recorded from the antibody doc.
MCLP_PROTEINS = [
    "ATMPS1981",
    "RAD17PS645",
    "CHK1PS345",
    "CHK2PT68",
    "ATM",
    "X53BP1",
    "RAD51",
    "RAD50",
    "CGAS",
    "IRF3",
    "IRF1",
    "HLADQA1",
    "HLADRDPDQDX",
    "PDL1",
    "HISTONEH3PS10",
    "CLAUDIN7",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url} -> {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def get_json(url: str):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read().decode())


def symbol_of(col: str) -> str:
    return col.split(" (")[0].strip()


def parse_gmt(path: Path, names: list[str]) -> dict[str, list[str]]:
    found: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if parts and parts[0] in names and len(parts) >= 3:
            found[parts[0]] = [g for g in parts[2:] if g]
    missing = [n for n in names if n not in found]
    if missing:
        raise SystemExit(f"GMT missing {missing}")
    return found


def extract_genes(expr_path: Path, wanted: set[str], out_path: Path) -> dict:
    with expr_path.open("r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
    found: dict[str, str] = {}
    keep_idx = [0]
    for i, col in enumerate(header[1:], start=1):
        sym = symbol_of(col)
        if sym in wanted and sym not in found:
            keep_idx.append(i)
            found[sym] = col
    missing = sorted(wanted - set(found))
    print(
        f"Expression columns kept: {len(found)} / {len(wanted)}; missing={missing}",
        flush=True,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    with expr_path.open("r", newline="") as fh, out_path.open("w", newline="") as out:
        reader = csv.reader(fh)
        writer = csv.writer(out)
        header = next(reader)
        writer.writerow(["ModelID"] + [symbol_of(header[i]) for i in keep_idx[1:]])
        for row in reader:
            writer.writerow([row[i] for i in keep_idx])
            n_rows += 1
    return {
        "n_models": n_rows,
        "n_genes_kept": len(found),
        "n_genes_wanted": len(wanted),
        "missing_genes": missing,
        "found_symbols": sorted(found),
    }


def rppa_header_inventory(path: Path, label: str, url: str) -> dict:
    with path.open("r", newline="") as fh:
        header = next(csv.reader(fh))
        n = sum(1 for _ in fh)
    antibodies = [h for h in header if h]
    hits = [
        h
        for h in antibodies
        if any(k in h.upper().replace(" ", "") for k in ("H2AX", "H2A.X", "GAMMAH2", "PH2AX"))
    ]
    return {
        "source": label,
        "url": url,
        "n_samples": n,
        "n_antibodies": len(antibodies),
        "h2ax_ps139_columns": hits,
        "gamma_h2ax_present": bool(hits),
    }


def mclp_inventory_and_matrix(outdir: Path) -> dict:
    ab = get_json(f"{MCLP_BASE}/CCLE-annotation-antibody")
    antibodies = ab["antibodies"]
    gamma = []
    for a in antibodies:
        blob = f"{a.get('protein_name','')} {a.get('label','')}".upper()
        if any(k in blob for k in ("H2AX", "H2A.X", "GAMMA-H2", "GAMMA H2")):
            gamma.append(a.get("protein_name"))
    histone = [
        a.get("protein_name")
        for a in antibodies
        if "HISTONE" in str(a.get("protein_name", "")).upper() or "H2A" in str(a.get("protein_name", "")).upper()
    ]
    by_label = {}
    for a in antibodies:
        label = str(a.get("label", "")).upper().replace("-", "")
        by_label[label] = {
            "protein_name": a.get("protein_name"),
            "validation_status": a.get("validation_status"),
            "genes": a.get("genes"),
            "catalog_number": a.get("catalog_number"),
            "source": a.get("source"),
        }

    samples_doc = get_json(f"{MCLP_BASE}/data-rppa-rppa-CCLE-all")
    samples = samples_doc["samples"]
    if samples_doc.get("data_release") != MCLP_RELEASE:
        print(f"NOTE MCLP data_release={samples_doc.get('data_release')}", flush=True)

    rows_by_sample = {s: {"sample": s} for s in samples}
    pulled = {}
    for pid in MCLP_PROTEINS:
        doc = get_json(f"{MCLP_BASE}/data-rppa-rppa-CCLE-all-{pid}")
        vals = doc["exprs"]
        if len(vals) != len(samples):
            raise SystemExit(f"{pid}: {len(vals)} values vs {len(samples)} samples")
        meta = by_label.get(pid, {})
        pulled[pid] = {
            "protein_id": pid,
            "protein_name": meta.get("protein_name"),
            "validation_status": meta.get("validation_status"),
            "genes": meta.get("genes"),
            "data_release": doc.get("data_release"),
            "data_norm_level": doc.get("data_norm_level"),
        }
        for s, v in zip(samples, vals):
            rows_by_sample[s][pid] = v
        print(f"MCLP {pid} {meta.get('protein_name')} {meta.get('validation_status')}", flush=True)

    matrix_path = outdir / "mclp_rppa_selected.csv"
    fieldnames = ["sample"] + MCLP_PROTEINS
    with matrix_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for s in samples:
            w.writerow(rows_by_sample[s])

    ann_path = outdir / "mclp_antibody_selected.json"
    ann_path.write_text(json.dumps(pulled, indent=2) + "\n")
    return {
        "source": "MCLP CCLE RPPA",
        "portal": MCLP_BASE,
        "data_release": samples_doc.get("data_release"),
        "data_norm_level": samples_doc.get("data_norm_level"),
        "n_samples": len(samples),
        "n_antibodies_annotated": len(antibodies),
        "gamma_h2ax_present": bool(gamma),
        "gamma_h2ax_name_hits": gamma,
        "histone_antibodies": histone,
        "proteins_pulled": pulled,
        "matrix": matrix_path.name,
    }


def main() -> int:
    outdir = Path(sys.argv[1] if len(sys.argv) > 1 else "data/ccle_cldn4_ddr_sting_hla")
    outdir.mkdir(parents=True, exist_ok=True)

    model_path = outdir / "Model.csv"
    gmt_path = outdir / "h.all.v2024.1.Hs.symbols.gmt"
    expr_full = outdir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    expr_slim = outdir / "expression_panel.csv"

    if not model_path.exists() or model_path.stat().st_size < 100_000:
        download(MODEL_URL, model_path)
    if not gmt_path.exists() or gmt_path.stat().st_size < 5_000:
        download(GMT_URL, gmt_path)

    hallmark = parse_gmt(gmt_path, HALLMARKS)
    for name, genes in hallmark.items():
        (outdir / f"{name}.txt").write_text("\n".join(genes) + "\n")

    wanted = set(CORE) | set(DSB_GENES) | set(STING_CORE) | set(STING_ISG) | set(HLA_GENES)
    for genes in hallmark.values():
        wanted |= set(genes)

    if not expr_slim.exists() or expr_slim.stat().st_size < 50_000:
        if not expr_full.exists() or expr_full.stat().st_size < 400_000_000:
            download(EXPR_URL, expr_full)
        extract_meta = extract_genes(expr_full, wanted, expr_slim)
        try:
            expr_full.unlink()
            print(f"Removed full matrix {expr_full.name}", flush=True)
        except OSError:
            pass
    else:
        with expr_slim.open() as fh:
            slim_header = next(csv.reader(fh))
        found_syms = slim_header[1:]
        extract_meta = {
            "n_models": sum(1 for _ in expr_slim.open()) - 1,
            "n_genes_kept": len(found_syms),
            "n_genes_wanted": len(wanted),
            "missing_genes": sorted(wanted - set(found_syms)),
            "found_symbols": found_syms,
        }
        print(f"OK exists {expr_slim}", flush=True)

    rppa_inv = []
    for label, url in RPPA_URLS.items():
        dest = outdir / f"{label}.csv"
        if not dest.exists() or dest.stat().st_size < 100_000:
            download(url, dest)
        rppa_inv.append(rppa_header_inventory(dest, label, url))
        # Header inventory is enough; drop the 2 MB matrix to keep the cache small.
        dest.unlink()

    mclp = mclp_inventory_and_matrix(outdir)

    manifest = {
        "release": RELEASE,
        "doi": DOI,
        "model_url": MODEL_URL,
        "expression_url": EXPR_URL,
        "gmt_url": GMT_URL,
        "expression_units": "log2(TPM+1)",
        "model_sha256": sha256_file(model_path),
        "expression_slim_sha256": sha256_file(expr_slim),
        "gmt_sha256": sha256_file(gmt_path),
        "hallmark_sizes": {k: len(v) for k, v in hallmark.items()},
        "extract": extract_meta,
        "gamma_h2ax_inventory": rppa_inv + [
            {
                "source": mclp["source"],
                "url": mclp["portal"],
                "data_release": mclp["data_release"],
                "n_samples": mclp["n_samples"],
                "n_antibodies": mclp["n_antibodies_annotated"],
                "gamma_h2ax_present": mclp["gamma_h2ax_present"],
                "histone_antibodies": mclp["histone_antibodies"],
            }
        ],
        "mclp": mclp,
        "note": (
            "H2AX transcript (H2AX/H2AFX) is not γH2AX. "
            "No public CCLE RPPA or MCLP 20221116 antibody is H2AX-pS139."
        ),
    }
    (outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({
        "n_models": extract_meta["n_models"],
        "missing_genes": extract_meta["missing_genes"],
        "gamma_h2ax_present": [x["gamma_h2ax_present"] for x in manifest["gamma_h2ax_inventory"]],
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
