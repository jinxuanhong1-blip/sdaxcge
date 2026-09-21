#!/usr/bin/env python3
"""Download TISMO in-vivo ICB tables for Cldn4, Tacstd2, and NHEJ/IFN genes.

One official Gene-module export per symbol (same call as PR #542). Writes a
wide sample x gene matrix. Raw per-gene CSVs stay in a local cache and are
not required for the scored tables.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tismo_client as tc

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = ROOT / "results" / "tismo_cldn4_nhej_ifn" / "tables"
CACHE = ROOT / "data" / "tismo_cldn4_nhej_ifn" / "raw"
META_COLS = [
    "Samples",
    "cell_line",
    "Responder",
    "Baseline",
    "GSE_ID",
    "Mouse_treatment",
    "count",
]

# MSigDB v2026.1 renamed symbols that TISMO still stores under the previous MGI symbol.
ALIASES = {
    "Wars1": "Wars",
    "Rigi": "Ddx58",
}


def genes_to_fetch(gene_sets: dict) -> list[str]:
    wanted = {"Cldn4", "Tacstd2"}
    for rec in gene_sets["sets"].values():
        wanted.update(rec["genes"])
    return sorted(wanted)


def fetch_one(gene: str, treatments: list[str], models: list[str]) -> tuple[str, str | None, str, str]:
    """Return MSigDB symbol, CSV text, error, and the TISMO symbol that responded."""
    tried = [gene]
    if gene in ALIASES:
        tried.append(ALIASES[gene])
    errors = []
    for symbol in tried:
        try:
            text = tc.download_vivo_expression(symbol, treatments, models)
            return gene, text, "", symbol
        except Exception as err:  # noqa: BLE001 — try the previous MGI symbol, then record the miss
            errors.append(f"{symbol}: {err}")
    return gene, None, " | ".join(errors), ""


def matrix_from_frames(frames: list[pd.DataFrame], meta: pd.DataFrame) -> pd.DataFrame:
    wide = pd.concat(frames, axis=1, join="outer")
    wide = wide.reset_index().merge(meta, on="Samples", how="left")
    symbols = sorted(c for c in wide.columns if c not in META_COLS and c != "Samples")
    return wide[META_COLS + symbols]


def metadata_union(parts: list[pd.DataFrame]) -> pd.DataFrame:
    """One metadata row per sample. Later genes fill samples the first gene lacked."""
    meta = pd.concat(parts, ignore_index=True)
    meta = meta.dropna(subset=["cell_line"])
    return meta.drop_duplicates("Samples", keep="first")


def parse_cached(msigdb_symbol: str) -> pd.DataFrame | None:
    source = msigdb_symbol
    path = CACHE / f"{msigdb_symbol}.csv"
    if not path.exists() and msigdb_symbol in ALIASES:
        source = ALIASES[msigdb_symbol]
        path = CACHE / f"{source}.csv"
    if not path.exists():
        return None
    df = tc.read_expression_csv(path, gene=source)
    if df.empty:
        return None
    df.attrs["tismo_symbol"] = source
    return df


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--from-cache",
        action="store_true",
        help="Rebuild the wide matrix from raw CSVs. Fetches only missing alias symbols.",
    )
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    gene_sets = json.loads((HERE / "gene_sets.json").read_text())
    genes = genes_to_fetch(gene_sets)
    print(f"[1] vocabularies; genes={len(genes)}")
    if (OUT / "vocabularies.json").exists() and args.from_cache:
        vocab = json.loads((OUT / "vocabularies.json").read_text())
        treatments = vocab["vivo_treatments"]
        models = vocab["vivo_models"]
    else:
        treatments = tc.names("/gene/getVivoTreatment")
        models = tc.names("/gene/getVivoCohort", {"treatment": ""})
        (OUT / "vocabularies.json").write_text(
            json.dumps({"vivo_treatments": treatments, "vivo_models": models}, indent=2) + "\n"
        )

    if args.from_cache:
        # Pull previous-symbol stand-ins that are not already cached under either name.
        for msigdb_symbol, tismo_symbol in ALIASES.items():
            if (CACHE / f"{msigdb_symbol}.csv").exists() or (CACHE / f"{tismo_symbol}.csv").exists():
                continue
            print(f"  alias fetch {msigdb_symbol} as {tismo_symbol}")
            text = tc.download_vivo_expression(tismo_symbol, treatments, models)
            (CACHE / f"{tismo_symbol}.csv").write_text(text)
        frames = []
        logs = []
        meta_parts = []
        for gene in genes:
            df = parse_cached(gene)
            if df is None:
                logs.append({"gene": gene, "status": "fail", "n_rows": 0, "tismo_symbol": "", "error": "not in cache"})
                continue
            keep = [c for c in META_COLS if c in df.columns]
            meta_parts.append(df[keep].drop_duplicates("Samples"))
            part = df[["Samples", "value"]].drop_duplicates("Samples").rename(columns={"value": gene})
            frames.append(part.set_index("Samples"))
            logs.append(
                {
                    "gene": gene,
                    "status": "ok",
                    "n_rows": int(len(part)),
                    "tismo_symbol": df.attrs.get("tismo_symbol", gene),
                    "n_malformed_dropped": int(df.attrs.get("dropped_malformed_rows", 0)),
                    "error": "",
                }
            )
        wide = matrix_from_frames(frames, metadata_union(meta_parts))
        wide.to_csv(OUT / "expression_wide.tsv.gz", sep="\t", index=False, compression="gzip")
        pd.DataFrame(logs).sort_values("gene").to_csv(OUT / "download_log.tsv", sep="\t", index=False)
        n_ok = sum(1 for row in logs if row["status"] == "ok")
        manifest = {
            "rebuilt_from_cache_unix": time.time(),
            "n_requested": len(genes),
            "n_ok": n_ok,
            "n_fail": len(genes) - n_ok,
            "n_samples": int(wide["Samples"].nunique()),
            "aliases": ALIASES,
            "api": tc.R_BASE + "/gene/downVivoExprn",
        }
        (OUT / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        print(json.dumps(manifest, indent=2))
        return 0

    print(f"[2] download {len(genes)} symbols")
    t0 = time.time()
    logs = []
    frames = []
    meta_parts = []
    done = 0
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = [pool.submit(fetch_one, gene, treatments, models) for gene in genes]
        for fut in as_completed(futures):
            gene, text, err, tismo_symbol = fut.result()
            done += 1
            if text is None:
                logs.append({"gene": gene, "status": "fail", "n_rows": 0, "tismo_symbol": "", "error": err})
                print(f"  FAIL {done}/{len(genes)} {gene}: {err[:160]}")
                continue
            (CACHE / f"{tismo_symbol}.csv").write_text(text)
            df = tc.read_expression_csv(text, gene=tismo_symbol)
            if df.empty:
                logs.append({"gene": gene, "status": "empty", "n_rows": 0, "error": "no rows"})
                print(f"  EMPTY {done}/{len(genes)} {gene}")
                continue
            keep = [c for c in META_COLS if c in df.columns]
            meta_parts.append(df[keep].drop_duplicates("Samples"))
            part = df[["Samples", "value"]].drop_duplicates("Samples").rename(columns={"value": gene})
            frames.append(part.set_index("Samples"))
            logs.append(
                {
                    "gene": gene,
                    "status": "ok",
                    "n_rows": int(len(part)),
                    "tismo_symbol": tismo_symbol,
                    "n_malformed_dropped": int(df.attrs.get("dropped_malformed_rows", 0)),
                    "error": "",
                }
            )
            if done % 25 == 0 or done == len(genes):
                print(f"  {done}/{len(genes)} in {time.time() - t0:.0f}s")

    if not frames or not meta_parts:
        raise SystemExit("no expression tables downloaded")
    wide = matrix_from_frames(frames, metadata_union(meta_parts))
    wide_path = OUT / "expression_wide.tsv.gz"
    wide.to_csv(wide_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(logs).sort_values("gene").to_csv(OUT / "download_log.tsv", sep="\t", index=False)
    n_ok = sum(1 for row in logs if row["status"] == "ok")
    manifest = {
        "downloaded_at_unix": time.time(),
        "n_requested": len(genes),
        "n_ok": n_ok,
        "n_fail": len(genes) - n_ok,
        "n_samples": int(wide["Samples"].nunique()),
        "seconds": round(time.time() - t0, 1),
        "api": tc.R_BASE + "/gene/downVivoExprn",
        "treatments": treatments,
        "n_models": len(models),
    }
    (OUT / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
