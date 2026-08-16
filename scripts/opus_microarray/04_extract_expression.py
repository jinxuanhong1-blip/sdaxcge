#!/usr/bin/env python3
"""Step 4 - download series matrices and extract target / control genes.

For every curated series we:
  1. download the GEO series matrix (cached, sha256 in the manifest)
  2. parse sample-level phenotype
  3. look up probe IDs from platform_probes.tsv (official-symbol matches only)
  4. stream those probe rows out of the matrix
  5. write per-series phenotype + expression tables

Outputs (under results/opus_microarray/)
  catalog_decisions.tsv
  extracted/<GSE>_pheno.tsv
  extracted/<GSE>_expr.tsv          probes x samples (raw matrix values)
  extracted/<GSE>_probe_map.tsv
  extracted/extraction_log.tsv
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from catalog import CATALOG  # noqa: E402
from common import (  # noqa: E402
    CACHE_DIR,
    CONTROL_GENES,
    RESULTS_DIR,
    TARGET_GENES,
    ensure_dirs,
    fetch,
    gse_matrix_url,
)
from geo_io import extract_probe_rows, list_all_probe_ids, parse_series_matrix_header  # noqa: E402

PREFERRED_COLS = {
    "Gene symbol",
    "ILMN_Gene",
    "Symbol",
    "ID",
    "gene_assignment",
    "GB_ACC",
}


def official_probe_map(probes: pd.DataFrame, platform: str) -> pd.DataFrame:
    """Keep high-confidence probe↔gene rows; drop alias-only false hits (MAL, p32)."""
    sub = probes[probes.platform == platform].copy()
    if sub.empty:
        return sub
    keep = []
    for _, r in sub.iterrows():
        col = str(r.matched_column)
        val = str(r.matched_value)
        gene = r.gene
        pid = str(r.probe_id)
        ok = False
        if col in {"Gene symbol", "ILMN_Gene", "Symbol"} and val.split("///")[0].strip() == gene:
            ok = True
        elif col == "ID" and pid == gene:
            ok = True
        elif col == "gene_assignment" and f" // {gene} // " in val:
            ok = True
        elif col == "Alias" and pid == gene:
            # NanoString: ID is the official symbol; Alias lists synonyms.
            ok = True
        if ok:
            keep.append(r)
    return pd.DataFrame(keep)


def download_matrix(gse: str) -> Path:
    dest = CACHE_DIR / "matrix" / f"{gse}_series_matrix.txt.gz"
    return fetch(gse_matrix_url(gse), dest)


def main() -> int:
    ensure_dirs()
    out_dir = RESULTS_DIR / "extracted"
    out_dir.mkdir(parents=True, exist_ok=True)

    probes = pd.read_csv(RESULTS_DIR / "platform_probes.tsv", sep="\t")
    cat = pd.DataFrame(CATALOG)
    cat.to_csv(RESULTS_DIR / "catalog_decisions.tsv", sep="\t", index=False)

    log_rows = []
    for rec in CATALOG:
        gse = rec["gse"]
        print(f"[{gse}] decision={rec['decision']}", flush=True)
        if rec["decision"] == "exclude":
            log_rows.append({"gse": gse, "status": "excluded", "detail": rec["reason"][:180]})
            continue
        try:
            mtx = download_matrix(gse)
        except Exception as exc:  # noqa: BLE001
            log_rows.append({"gse": gse, "status": "download_fail", "detail": str(exc)})
            print(f"  download failed: {exc}")
            continue
        parsed = parse_series_matrix_header(mtx)
        pheno = parsed["pheno"]
        series = parsed["series"]
        pheno.to_csv(out_dir / f"{gse}_pheno.tsv", sep="\t")
        plat = rec.get("platform") or series.get("platform", "")
        pmap = official_probe_map(probes, plat)
        # Always also do an exact ID-level scan of the matrix for the gene symbols
        # themselves (NanoString / custom panels).
        all_ids = None
        extra = []
        if rec["decision"] == "document_absent" or plat in {"GPL19965", "GPL29738"}:
            all_ids = list_all_probe_ids(mtx)
            idset = {i.upper() for i in all_ids}
            for g in TARGET_GENES + CONTROL_GENES:
                extra.append({"gse": gse, "gene": g, "in_matrix_ids": int(g in idset)})
        pd.DataFrame(extra).to_csv(out_dir / f"{gse}_id_scan.tsv", sep="\t", index=False) if extra else None

        genes_wanted = list(TARGET_GENES)
        if rec.get("analyze_controls"):
            genes_wanted += list(CONTROL_GENES)
        pmap = pmap[pmap.gene.isin(genes_wanted)]
        pmap.to_csv(out_dir / f"{gse}_probe_map.tsv", sep="\t", index=False)

        probe_ids = set(pmap.probe_id.astype(str)) if len(pmap) else set()
        # For NanoString, probe ID == gene symbol for genes that are on the panel.
        if all_ids is None and plat in {"GPL19965", "GPL29738"}:
            all_ids = list_all_probe_ids(mtx)
        if all_ids is not None:
            idset = set(all_ids)
            for g in genes_wanted:
                if g in idset:
                    probe_ids.add(g)

        expr = extract_probe_rows(mtx, probe_ids)
        expr.to_csv(out_dir / f"{gse}_expr.tsv", sep="\t")
        log_rows.append(
            {
                "gse": gse,
                "status": "ok",
                "n_samples": int(len(pheno)),
                "n_probes_extracted": int(expr.shape[0]),
                "probes": ";".join(expr.index.astype(str)),
                "platform": plat,
                "title": series.get("title", ""),
                "pubmed": series.get("pubmed", ""),
                "pheno_fields": ";".join(series.get("pheno_fields", [])),
            }
        )
        print(f"  n={len(pheno)}  probes={expr.shape[0]}  fields={series.get('pheno_fields')}")

    pd.DataFrame(log_rows).to_csv(RESULTS_DIR / "extracted" / "extraction_log.tsv", sep="\t", index=False)
    print("[write] extracted/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
