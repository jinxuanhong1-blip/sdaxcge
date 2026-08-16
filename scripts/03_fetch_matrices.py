#!/usr/bin/env python3
"""Stream recount3 gene_sums files and keep only the hunt genes.

One pass per cohort:
  * accumulate per-sample library size (sum of all gene counts)
  * keep the rows listed in scripts/genes.json
  * write a compact TSV: gene symbol x sample raw counts, plus a library-size row

Samples not in the nominated lung-run list are dropped for SRA studies. GTEx and
TCGA keep every column; TCGA tumor/normal is split in the analysis step.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np


TCGA_META = {
    "LUAD": "https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/metadata/AD/LUAD/tcga.tcga.LUAD.MD.gz",
    "LUSC": "https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/metadata/SC/LUSC/tcga.tcga.LUSC.MD.gz",
}


def fetch_tcga_maps(table_dir: Path) -> None:
    """UUID / barcode -> Primary Tumor vs Solid Tissue Normal."""
    table_dir.mkdir(parents=True, exist_ok=True)
    dest = table_dir / "tcga_sample_map.tsv"
    if dest.exists() and dest.stat().st_size > 100:
        return
    rows = ["external_id\ttcga_barcode\tsample_type\tproject\n"]
    for project, url in TCGA_META.items():
        req = urllib.request.Request(url, headers={"User-Agent": "hunt_tf/1.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            with gzip.GzipFile(fileobj=resp) as raw:
                header = raw.readline().decode("ascii", errors="replace").rstrip("\n").split("\t")
                idx = {n: i for i, n in enumerate(header)}
                for line in raw:
                    rec = line.decode("ascii", errors="replace").rstrip("\n").split("\t")
                    ext = rec[idx.get("external_id", 1)]
                    barcode = rec[idx["tcga_barcode"]] if "tcga_barcode" in idx else ""
                    st = rec[idx["cgc_sample_sample_type"]] if "cgc_sample_sample_type" in idx else ""
                    rows.append(f"{ext}\t{barcode}\t{st}\t{project}\n")
    dest.write_text("".join(rows))
    print(f"wrote {dest} n={len(rows)-1}", flush=True)


def stream_extract(url: str, id_to_symbol: dict[str, str], keep_samples: set[str] | None):
    req = urllib.request.Request(url, headers={"User-Agent": "hunt_tf/1.0"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        with gzip.GzipFile(fileobj=resp) as raw:
            header = None
            col_idx = None
            names = None
            lib = None
            kept: dict[str, np.ndarray] = {}
            n_genes = 0
            for line in raw:
                line = line.decode("ascii", errors="replace")
                if line.startswith("#"):
                    continue
                if header is None:
                    header = line.rstrip("\n").split("\t")
                    sample_names = header[1:]
                    if keep_samples:
                        col_idx = [i for i, s in enumerate(sample_names) if s.split(".")[0] in keep_samples or s in keep_samples]
                        if not col_idx:
                            # recount3 SRA columns are usually the SRR / ERR / DRR id
                            col_idx = [i for i, s in enumerate(sample_names) if any(s.startswith(k) or k.startswith(s) for k in list(keep_samples)[:3])]
                        if not col_idx:
                            col_idx = list(range(len(sample_names)))
                        names = [sample_names[i] for i in col_idx]
                    else:
                        col_idx = list(range(len(sample_names)))
                        names = sample_names
                    lib = np.zeros(len(col_idx), dtype=np.float64)
                    continue
                gid, _, rest = line.partition("\t")
                vals = np.fromstring(rest, sep="\t", dtype=np.float64)
                if vals.size != len(header) - 1:
                    continue
                sub = vals[col_idx]
                lib += sub
                n_genes += 1
                symbol = id_to_symbol.get(gid) or id_to_symbol.get(gid.split(".")[0])
                if symbol is not None:
                    kept[symbol] = sub
    return names, lib, kept, n_genes


def write_matrix(path: Path, names, lib, kept) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    symbols = sorted(kept)
    with path.open("w") as fh:
        fh.write("gene\t" + "\t".join(names) + "\n")
        fh.write("__libsize__\t" + "\t".join(f"{x:.0f}" for x in lib) + "\n")
        for sym in symbols:
            fh.write(sym + "\t" + "\t".join(f"{x:.0f}" for x in kept[sym]) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", default="results/hunt_tf/tables/cohorts_to_download.tsv")
    ap.add_argument("--genes", default="scripts/genes.json")
    ap.add_argument("--out-dir", default="results/hunt_tf/matrices")
    ap.add_argument("--cache-dir", default=str(Path.home() / "hunt_data" / "matrices"))
    ap.add_argument("--limit", type=int, default=0, help="optional cap for smoke tests")
    args = ap.parse_args()

    genes = json.loads(Path(args.genes).read_text())
    rows = list(csv.DictReader(open(args.cohorts), delimiter="\t"))
    fetch_tcga_maps(Path(args.out_dir).parent / "tables")
    if args.limit:
        rows = rows[: args.limit]

    out_dir = Path(args.out_dir)
    cache_dir = Path(args.cache_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    manifest = []

    for rec in rows:
        cid = rec["cohort_id"]
        org = rec["organism"]
        dest = cache_dir / f"{cid}.counts.tsv"
        pub = out_dir / f"{cid}.counts.tsv"
        if dest.exists() and dest.stat().st_size > 200:
            print(f"skip existing {cid}", flush=True)
            if not pub.exists():
                pub.write_bytes(dest.read_bytes())
            manifest.append({**rec, "status": "cached", "n_samples_written": "NA"})
            continue

        id_to_symbol = {gid: sym for sym, gid in genes[org].items()}
        id_to_symbol.update({gid.split(".")[0]: sym for sym, gid in genes[org].items()})
        keep = set(rec["selected_runs"].split(",")) if rec["selected_runs"] else None
        print(f"fetch {cid} {rec['url']}", flush=True)
        t0 = time.time()
        try:
            names, lib, kept, n_genes = stream_extract(rec["url"], id_to_symbol, keep)
        except Exception as exc:
            print(f"FAIL {cid}: {exc}", file=sys.stderr, flush=True)
            manifest.append({**rec, "status": f"fail:{exc}", "n_samples_written": 0})
            continue
        write_matrix(dest, names, lib, kept)
        pub.write_bytes(dest.read_bytes())
        print(
            f"  wrote {cid} samples={len(names)} genes_kept={len(kept)}/{len(genes[org])} "
            f"rows_scanned={n_genes} sec={time.time()-t0:.1f}",
            flush=True,
        )
        manifest.append({**rec, "status": "ok", "n_samples_written": len(names), "n_genes_kept": len(kept)})

    Path(out_dir.parent, "tables", "fetch_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
