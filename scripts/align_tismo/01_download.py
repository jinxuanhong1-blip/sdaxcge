#!/usr/bin/env python3
"""Download the TISMO public tables needed to replicate the Tacstd2 ICB claim.

Writes raw, unmodified payloads to results/align_tismo/data/ so every downstream
number is traceable to a file that can be re-fetched.

  --null-genes N   also download N randomly chosen genes (seeded) to build an
                   empirical genome-wide null for the same paired statistic.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import tismo_client as tc

GENE_OF_INTEREST = "Tacstd2"
NULL_SEED = 20240816


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    print(f"  wrote {path.relative_to(tc.REPO_ROOT)} ({len(text):,} bytes)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--null-genes", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()

    tc.DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/5] vocabularies")
    treatments = tc.get_vivo_treatments()
    models = tc.get_vivo_cohorts(treatments)
    genes = tc.get_gene_list()
    _write(
        tc.DATA_DIR / "vocabularies.json",
        json.dumps(
            {
                "icb_treatments": treatments,
                "tumor_models": models,
                "n_genes": len(genes),
                "api_meta_base": tc.META_BASE,
                "api_r_base": tc.R_BASE,
            },
            indent=2,
        ),
    )
    _write(tc.DATA_DIR / "gene_list.txt", "\n".join(genes) + "\n")

    print("[2/5] metadata tables")
    for kind in ("vivoMeta", "vitroMeta", "cellLineMeta"):
        rows = tc.get_metadata(kind)
        _write(tc.DATA_DIR / f"{kind}.json", json.dumps(rows, indent=2))

    print(f"[3/5] {GENE_OF_INTEREST} in vivo expression")
    text = tc.download_vivo_expression(GENE_OF_INTEREST, treatments, models)
    _write(tc.DATA_DIR / f"vivo_expression_{GENE_OF_INTEREST}.csv", text)

    print("[4/5] positive/negative reference genes")
    ref_dir = tc.DATA_DIR / "reference_genes"
    for gene in ("Cd274", "Pdcd1", "Ifng", "Cd8a", "Epcam", "Krt8", "Actb"):
        if gene not in genes:
            print(f"  skip {gene} (absent from TISMO gene list)")
            continue
        _write(ref_dir / f"vivo_expression_{gene}.csv",
               tc.download_vivo_expression(gene, treatments, models))

    if args.null_genes <= 0:
        print("[5/5] null panel skipped (--null-genes 0)")
        return 0

    print(f"[5/5] null panel: {args.null_genes} random genes, {args.workers} workers")
    null_dir = tc.DATA_DIR / "null_genes"
    null_dir.mkdir(parents=True, exist_ok=True)
    pool = [g for g in genes if g != GENE_OF_INTEREST]
    rng = random.Random(NULL_SEED)
    chosen = sorted(rng.sample(pool, min(args.null_genes, len(pool))))
    _write(tc.DATA_DIR / "null_gene_panel.txt", "\n".join(chosen) + "\n")

    todo = [g for g in chosen if not (null_dir / f"{g}.csv").exists()]
    print(f"  {len(chosen) - len(todo)} already cached, fetching {len(todo)}")

    failures: list[str] = []

    def fetch(gene: str) -> tuple[str, str | None]:
        try:
            return gene, tc.download_vivo_expression(gene, treatments, models)
        except Exception as err:  # noqa: BLE001 - a missing gene must not kill the panel
            return gene, f"ERROR:{err}"

    done = 0
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for gene, payload in ex.map(fetch, todo):
            done += 1
            if payload is None or payload.startswith("ERROR:"):
                failures.append(gene)
            else:
                (null_dir / f"{gene}.csv").write_text(payload)
            if done % 25 == 0:
                print(f"  {done}/{len(todo)} ({len(failures)} failed)")

    print(f"  finished: {len(todo) - len(failures)} downloaded, {len(failures)} failed")
    if failures:
        _write(tc.DATA_DIR / "null_gene_failures.txt", "\n".join(failures) + "\n")

    tc.pack_null_genes()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
