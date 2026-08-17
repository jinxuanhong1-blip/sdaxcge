#!/usr/bin/env python3
"""Public TF–target priors for IFN / MHC-I / TJ / keratin TFs. Not ChIP.

A snapshot is committed under resources/. ELF3–CLDN4 is recorded as A10-given.
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "resources"

# Program TFs. ELF3 is A10-given and is never the headline.
TF_PROGRAM = {
    "STAT1": "IFN",
    "STAT2": "IFN",
    "IRF1": "IFN_MHC",
    "IRF7": "IFN",
    "IRF9": "IFN",
    "NLRC5": "MHC_I",
    "RFX5": "MHC_I",
    "GRHL2": "TJ",
    "OVOL1": "TJ",
    "OVOL2": "TJ",
    "KLF5": "TJ",
    "ELF3": "TJ_GIVEN",
    "TP63": "keratin",
    "GRHL1": "keratin",
    "KLF4": "keratin",
    "NKX2-1": "control",
    "SOX2": "control",
}
FOCUS = list(TF_PROGRAM)
URLS = {
    "trrust": "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
    "dorothea": (
        "https://omnipathdb.org/interactions?datasets=dorothea"
        "&genesymbols=1&fields=sources,references,dorothea_level"
    ),
    "collectri": (
        "https://omnipathdb.org/interactions?datasets=collectri"
        "&genesymbols=1&fields=sources,references"
    ),
}


def fetch(url: str) -> bytes:
    req = urllib.request.Request(
        url, headers={"User-Agent": "sdaxcge-merge-131907-205335-scenic-cldn4"}
    )
    with urllib.request.urlopen(req, timeout=180) as r:
        return r.read()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    raw = fetch(URLS["trrust"])
    tr = pd.read_csv(
        pd.io.common.BytesIO(raw), sep="\t", header=None,
        names=["tf", "target", "effect", "pmid"],
    )
    for r in tr.itertuples(index=False):
        if r.tf in TF_PROGRAM:
            rows.append(dict(
                tf=r.tf, target=r.target, resource="TRRUST_v2",
                confidence="", effect=r.effect, reference=str(r.pmid),
                program=TF_PROGRAM[r.tf],
                note="curated literature TF-target; not a ChIP peak",
            ))
    try:
        doro = pd.read_csv(pd.io.common.BytesIO(fetch(URLS["dorothea"])), sep="\t")
        for r in doro.itertuples(index=False):
            if r.source_genesymbol in TF_PROGRAM:
                eff = "Activation" if r.is_stimulation else (
                    "Repression" if r.is_inhibition else "Unknown"
                )
                rows.append(dict(
                    tf=r.source_genesymbol, target=r.target_genesymbol,
                    resource="DoRothEA", confidence=str(r.dorothea_level),
                    effect=eff, reference=str(r.references)[:300],
                    program=TF_PROGRAM[r.source_genesymbol],
                    note="DoRothEA via OmniPath; not a lung ChIP peak",
                ))
    except Exception as exc:
        print(f"[warn] DoRothEA skipped: {exc}", flush=True)
    try:
        col = pd.read_csv(pd.io.common.BytesIO(fetch(URLS["collectri"])), sep="\t")
        for r in col.itertuples(index=False):
            if r.source_genesymbol in TF_PROGRAM:
                eff = "Activation" if r.is_stimulation else (
                    "Repression" if r.is_inhibition else "Unknown"
                )
                rows.append(dict(
                    tf=r.source_genesymbol, target=r.target_genesymbol,
                    resource="CollecTRI", confidence="", effect=eff,
                    reference=str(r.references)[:300],
                    program=TF_PROGRAM[r.source_genesymbol],
                    note="CollecTRI via OmniPath; not a lung ChIP peak",
                ))
    except Exception as exc:
        print(f"[warn] CollecTRI skipped: {exc}", flush=True)

    df = pd.DataFrame(rows).drop_duplicates(["tf", "target", "resource"])
    df = df.sort_values(["program", "tf", "resource", "target"])
    df.to_csv(OUT / "tf_targets_public.tsv", sep="\t", index=False)
    manifest = {
        "purpose": (
            "Public curated TF-target edges for merged GSE131907+GSE205335 "
            "CLDN4-high GRN proxy. Not ChIP. Not pySCENIC cisTarget."
        ),
        "do_not_invent_chip": True,
        "urls": URLS,
        "n_rows": int(len(df)),
        "tf_program": TF_PROGRAM,
        "cldn4_listed_as_target": bool(df.target.eq("CLDN4").any()),
        "elf3_cldn4_in_prior": bool(((df.tf == "ELF3") & (df.target == "CLDN4")).any()),
        "a10_elf3_cldn4": "given; not a discovery",
    }
    (OUT / "prior_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(df.groupby(["tf", "resource"]).size().unstack(fill_value=0).to_string())
    print("wrote", OUT / "tf_targets_public.tsv", "n=", len(df))


if __name__ == "__main__":
    main()
