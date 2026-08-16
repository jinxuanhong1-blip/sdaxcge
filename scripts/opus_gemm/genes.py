#!/usr/bin/env python3
"""Genes and gene signatures used by the slice, plus symbol -> Ensembl ID
resolution against the Ensembl REST API (cached to a committed TSV so the
mapping used for the published numbers is auditable).

Run standalone to (re)build the cache:
    python3 scripts/opus_gemm/genes.py --out results/opus_gemm/gene_ids.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
import urllib.request

# Focus genes for this slice.
FOCUS = ["Tacstd2", "Cldn4"]

# Claudin / epithelial-surface context around the focus genes.
CONTEXT = [
    "Cldn1", "Cldn3", "Cldn5", "Cldn6", "Cldn7", "Cldn18", "Epcam", "Cdh1",
    "Krt8", "Krt18", "Krt5", "Muc1", "Sox2", "Nkx2-1", "Vim", "Zeb1", "Snai1",
]

# Signatures scored per sample as the mean z-score of their members.
SIGNATURES: dict[str, list[str]] = {
    "CD8_T_cell": ["Cd8a", "Cd8b1", "Cd3d", "Cd3e", "Cd2", "Gzmk", "Prf1", "Gzmb"],
    "IFNg_response": [
        "Stat1", "Cxcl9", "Cxcl10", "Ido1", "Irf1", "Ciita", "Nlrc5", "Socs1",
        "Gbp2", "Gbp3", "Ifi47", "Igtp", "Iigp1",
    ],
    "IFN_I_response": ["Ifit1", "Ifit3", "Isg15", "Oasl2", "Mx1", "Mx2", "Irf7", "Rsad2", "Usp18"],
    "MHC_I_antigen_presentation": [
        "B2m", "H2-K1", "H2-D1", "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9", "Nlrc5",
    ],
    "MHC_II": ["H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "Ciita"],
    "T_cell_inhibitory": ["Cd274", "Pdcd1lg2", "Pdcd1", "Ctla4", "Lag3", "Havcr2", "Tigit", "Foxp3"],
    "Myeloid_suppressive": ["Arg1", "Mrc1", "Trem2", "Spp1", "Il10", "Vegfa", "Cd68", "Adgre1"],
    "Leukocyte_general": ["Ptprc", "Cd52", "Laptm5", "Coro1a", "Cd53"],
    "Neuroendocrine": ["Ascl1", "Insm1", "Chga", "Syp", "Ncam1", "Calca", "Dll3", "Uchl1"],
    "Non_NE_SCLC": ["Pou2f3", "Neurod1", "Yap1", "Myc", "Rest", "Vim"],
}

EXTRA = ["Ifng", "Ifngr1", "Ifngr2", "Ifnb1", "Tmem173", "Sting1", "Ahr", "Id1", "Lair1",
         "Hdac1", "Kdm1a", "Met", "Egfr", "Kras", "Trp53", "Rb1", "Stk11", "Actb", "Gapdh"]


def all_genes() -> list[str]:
    seen: dict[str, None] = {}
    for g in FOCUS + CONTEXT + EXTRA + [g for v in SIGNATURES.values() for g in v]:
        seen.setdefault(g, None)
    return list(seen)


def _post(url: str, payload: dict, tries: int = 4) -> dict:
    data = json.dumps(payload).encode()
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                url, data=data, headers={"Content-Type": "application/json", "Accept": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=120) as fh:
                return json.loads(fh.read().decode())
        except Exception as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Ensembl POST failed: {last}")


def _get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return json.loads(fh.read().decode())


def resolve(symbols: list[str]) -> tuple[list[dict], str]:
    """Resolve mouse gene symbols to Ensembl gene IDs. Returns rows + release."""
    release = str(_get("https://rest.ensembl.org/info/data?content-type=application/json")["releases"][0])
    rows: list[dict] = []
    for i in range(0, len(symbols), 100):
        chunk = symbols[i : i + 100]
        res = _post("https://rest.ensembl.org/lookup/symbol/mus_musculus", {"symbols": chunk})
        for sym in chunk:
            rec = res.get(sym)
            rows.append(
                {
                    "symbol": sym,
                    "ensembl_gene_id": (rec or {}).get("id", ""),
                    "display_name": (rec or {}).get("display_name", ""),
                    "biotype": (rec or {}).get("biotype", ""),
                    "assembly": (rec or {}).get("assembly_name", ""),
                    "chrom": str((rec or {}).get("seq_region_name", "")),
                    "resolved": "yes" if rec else "no",
                }
            )
        time.sleep(0.3)
    return rows, release


def load_ids(path: str = "results/opus_gemm/gene_ids.tsv") -> dict[str, str]:
    """symbol -> Ensembl gene ID (only resolved entries)."""
    out: dict[str, str] = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        if r.get("ensembl_gene_id"):
            out[r["symbol"]] = r["ensembl_gene_id"]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/opus_gemm/gene_ids.tsv")
    args = ap.parse_args()
    syms = all_genes()
    rows, release = resolve(syms)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        fh.write(f"# Ensembl REST lookup/symbol/mus_musculus, release {release}\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    unresolved = [r["symbol"] for r in rows if r["resolved"] == "no"]
    print(f"resolved {len(rows) - len(unresolved)}/{len(rows)} symbols (Ensembl release {release})")
    if unresolved:
        print("unresolved:", ", ".join(unresolved), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
