#!/usr/bin/env python3
"""Public TF–target priors for IFN / MHC / TJ TFs (TRRUST / DoRothEA / CollecTRI).

Not ChIP. Not cisTarget. Not a binding map. CLDN4 is recorded if listed
as a target and is held out of AUCell regulons in analyze.py.
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "resources"

IFN_TFS = ["STAT1", "STAT2", "IRF1", "IRF3", "IRF7", "IRF9"]
MHC_TFS = ["NLRC5", "CIITA", "RFX5", "RFXANK", "RFXAP"]
TJ_TFS = ["GRHL1", "GRHL2", "ELF3", "KLF4", "OVOL1", "OVOL2"]
CTRL_TFS = ["HIF1A"]  # hypoxia TF; specificity, not a TJ/IFN claim
FOCUS = sorted(set(IFN_TFS + MHC_TFS + TJ_TFS + CTRL_TFS + ["IRF1"]))

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
    req = urllib.request.Request(url, headers={"User-Agent": "gse148071-scenic-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return r.read()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    errors = {}

    try:
        raw = fetch(URLS["trrust"])
        tr = pd.read_csv(
            pd.io.common.BytesIO(raw),
            sep="\t",
            header=None,
            names=["tf", "target", "effect", "pmid"],
        )
        for r in tr.itertuples(index=False):
            if r.tf in FOCUS:
                rows.append(
                    dict(
                        tf=r.tf,
                        target=r.target,
                        resource="TRRUST_v2",
                        confidence="",
                        effect=r.effect,
                        reference=str(r.pmid),
                        note="curated literature TF-target; not a ChIP peak",
                    )
                )
        print(f"TRRUST rows kept: {sum(1 for x in rows if x['resource']=='TRRUST_v2')}", flush=True)
    except Exception as exc:
        errors["trrust"] = repr(exc)
        print(f"TRRUST failed: {exc}", flush=True)

    try:
        doro = pd.read_csv(pd.io.common.BytesIO(fetch(URLS["dorothea"])), sep="\t")
        n0 = len(rows)
        for r in doro.itertuples(index=False):
            if r.source_genesymbol in FOCUS:
                eff = (
                    "Activation"
                    if r.is_stimulation
                    else ("Repression" if r.is_inhibition else "Unknown")
                )
                rows.append(
                    dict(
                        tf=r.source_genesymbol,
                        target=r.target_genesymbol,
                        resource="DoRothEA",
                        confidence=str(r.dorothea_level),
                        effect=eff,
                        reference=str(r.references)[:300],
                        note="DoRothEA via OmniPath; not a lung ChIP peak",
                    )
                )
        print(f"DoRothEA rows kept: {len(rows) - n0}", flush=True)
    except Exception as exc:
        errors["dorothea"] = repr(exc)
        print(f"DoRothEA failed: {exc}", flush=True)

    try:
        col = pd.read_csv(pd.io.common.BytesIO(fetch(URLS["collectri"])), sep="\t")
        n0 = len(rows)
        for r in col.itertuples(index=False):
            if r.source_genesymbol in FOCUS:
                eff = (
                    "Activation"
                    if r.is_stimulation
                    else ("Repression" if r.is_inhibition else "Unknown")
                )
                rows.append(
                    dict(
                        tf=r.source_genesymbol,
                        target=r.target_genesymbol,
                        resource="CollecTRI",
                        confidence="",
                        effect=eff,
                        reference=str(r.references)[:300],
                        note="CollecTRI via OmniPath; not a lung ChIP peak",
                    )
                )
        print(f"CollecTRI rows kept: {len(rows) - n0}", flush=True)
    except Exception as exc:
        errors["collectri"] = repr(exc)
        print(f"CollecTRI failed: {exc}", flush=True)

    df = pd.DataFrame(rows)
    if len(df):
        df = df.drop_duplicates(["tf", "target", "resource"])
        df = df.sort_values(["tf", "resource", "target"])
    else:
        df = pd.DataFrame(
            columns=["tf", "target", "resource", "confidence", "effect", "reference", "note"]
        )
    df.to_csv(OUT / "tf_targets_public.tsv", sep="\t", index=False)
    counts = (
        df.groupby(["tf", "resource"]).size().unstack(fill_value=0).to_dict()
        if len(df)
        else {}
    )
    manifest = {
        "purpose": "Public curated TF-target edges for IFN/MHC/TJ TFs. Not ChIP peaks.",
        "method": "AUCell / TF-target proxy — not full pySCENIC (no cisTarget DB)",
        "do_not_invent_chip": True,
        "tfs": FOCUS,
        "ifn_tfs": IFN_TFS,
        "mhc_tfs": MHC_TFS,
        "tj_tfs": TJ_TFS,
        "urls": URLS,
        "n_rows": int(len(df)),
        "n_unique_tf_target": int(df.drop_duplicates(["tf", "target"]).shape[0]) if len(df) else 0,
        "cldn4_listed_as_target": bool(len(df) and df.target.eq("CLDN4").any()),
        "tacstd2_listed_as_target": bool(len(df) and df.target.eq("TACSTD2").any()),
        "errors": errors,
        "counts_by_tf_resource": counts,
    }
    (OUT / "prior_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ("n_rows", "n_unique_tf_target", "errors")}, indent=2))
    print("wrote", OUT / "tf_targets_public.tsv")


if __name__ == "__main__":
    main()
