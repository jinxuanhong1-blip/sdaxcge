#!/usr/bin/env python3
"""Refresh public TF–target priors (TRRUST / DoRothEA / CollecTRI). Not ChIP."""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "resources"
FOCUS = ["ELF3", "GRHL1", "KLF4", "TFAP2A", "NKX2-1"]
URLS = {
    "trrust": "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
    "dorothea": "https://omnipathdb.org/interactions?datasets=dorothea&genesymbols=1&fields=sources,references,dorothea_level",
    "collectri": "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&fields=sources,references",
}


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    raw = fetch(URLS["trrust"])
    tr = pd.read_csv(pd.io.common.BytesIO(raw), sep="\t", header=None,
                     names=["tf", "target", "effect", "pmid"])
    for r in tr.itertuples(index=False):
        if r.tf in FOCUS:
            rows.append(dict(tf=r.tf, target=r.target, resource="TRRUST_v2",
                             confidence="", effect=r.effect, reference=str(r.pmid),
                             note="curated literature TF-target; not a ChIP peak"))
    doro = pd.read_csv(pd.io.common.BytesIO(fetch(URLS["dorothea"])), sep="\t")
    for r in doro.itertuples(index=False):
        if r.source_genesymbol in FOCUS:
            eff = "Activation" if r.is_stimulation else ("Repression" if r.is_inhibition else "Unknown")
            rows.append(dict(tf=r.source_genesymbol, target=r.target_genesymbol,
                             resource="DoRothEA", confidence=str(r.dorothea_level),
                             effect=eff, reference=str(r.references)[:300],
                             note="DoRothEA via OmniPath; not a lung ChIP peak"))
    col = pd.read_csv(pd.io.common.BytesIO(fetch(URLS["collectri"])), sep="\t")
    for r in col.itertuples(index=False):
        if r.source_genesymbol in FOCUS:
            eff = "Activation" if r.is_stimulation else ("Repression" if r.is_inhibition else "Unknown")
            rows.append(dict(tf=r.source_genesymbol, target=r.target_genesymbol,
                             resource="CollecTRI", confidence="", effect=eff,
                             reference=str(r.references)[:300],
                             note="CollecTRI via OmniPath; not a lung ChIP peak"))
    df = pd.DataFrame(rows).drop_duplicates(["tf", "target", "resource"])
    df = df.sort_values(["tf", "resource", "target"])
    df.to_csv(OUT / "tf_targets_public.tsv", sep="\t", index=False)
    manifest = {
        "purpose": "Public curated TF-target edges. Not ChIP peaks.",
        "do_not_invent_chip": True,
        "urls": URLS,
        "n_rows": int(len(df)),
        "tacstd2_cldn4_listed_as_targets": bool(df.target.isin(["TACSTD2", "CLDN4"]).any()),
    }
    (OUT / "prior_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(df.groupby(["tf", "resource"]).size().unstack(fill_value=0))
    print("wrote", OUT / "tf_targets_public.tsv")


if __name__ == "__main__":
    main()
