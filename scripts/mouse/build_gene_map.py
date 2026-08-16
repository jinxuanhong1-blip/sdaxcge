#!/usr/bin/env python3
"""Fetch verified mouse Ensembl gene IDs for target + immune-module genes.

Used to look up genes in the Ensembl-keyed E-MTAB-13704 count matrix.
Cached to notes/mouse/gene_map.json so we never invent identifiers.
"""
import json
import os
import time
import urllib.request

OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "notes", "mouse", "gene_map.json"))

TARGETS = ["Tacstd2", "Cldn4"]
IMMUNE = [
    "Cd8a", "Cd8b1", "Gzmb", "Gzmk", "Prf1", "Ifng", "Nkg7", "Pdcd1",
    "Cd274", "Cxcl9", "Cxcl10", "Cd3e", "Cd4", "Foxp3", "Ctla4",
    "Havcr2", "Lag3", "Ptprc",
]
UA = {"User-Agent": "Mozilla/5.0 (mouse-ici-catalog)"}


def sym2ens(sym):
    url = f"https://rest.ensembl.org/xrefs/symbol/mus_musculus/{sym}?content-type=application/json"
    req = urllib.request.Request(url, headers=UA)
    for _ in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                data = json.load(r)
            for d in data:
                if d.get("type") == "gene":
                    return d["id"]
            return None
        except Exception:
            time.sleep(2)
    return None


def main():
    m = {}
    for sym in TARGETS + IMMUNE:
        eid = sym2ens(sym)
        m[sym] = eid
        print(sym, "->", eid)
        time.sleep(0.2)
    with open(OUT, "w") as f:
        json.dump({"targets": TARGETS, "immune_module": IMMUNE, "symbol_to_ensembl": m}, f, indent=2)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
