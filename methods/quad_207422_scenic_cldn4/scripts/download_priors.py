#!/usr/bin/env python3
"""Download public TF–target priors for IFN / MHC / TJ TFs. Not ChIP."""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "resources"

IFN_TFS = ["IRF1", "IRF7", "IRF9", "STAT1", "STAT2"]
MHC_TFS = ["NLRC5", "CIITA", "RFX5"]
TJ_TFS = ["ELF3", "GRHL1", "GRHL2", "KLF4", "OVOL1", "OVOL2"]
FOCUS = IFN_TFS + MHC_TFS + TJ_TFS

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

AXIS = {tf: "IFN" for tf in IFN_TFS}
AXIS.update({tf: "MHC" for tf in MHC_TFS})
AXIS.update({tf: "TJ" for tf in TJ_TFS})


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=180) as r:
        return r.read()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = Path("/tmp/quad_scenic/priors")
    rows = []

    tr_path = cache / "trrust_rawdata.human.tsv"
    raw = tr_path.read_bytes() if tr_path.exists() else fetch(URLS["trrust"])
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
                    axis=AXIS[r.tf],
                    resource="TRRUST_v2",
                    confidence="",
                    effect=r.effect,
                    reference=str(r.pmid),
                    note="curated literature TF-target; not a ChIP peak",
                )
            )

    doro_path = cache / "dorothea.tsv"
    doro_raw = doro_path.read_bytes() if doro_path.exists() else fetch(URLS["dorothea"])
    doro = pd.read_csv(pd.io.common.BytesIO(doro_raw), sep="\t")
    for r in doro.itertuples(index=False):
        tf = r.source_genesymbol
        if tf not in FOCUS:
            continue
        eff = (
            "Activation"
            if r.is_stimulation
            else ("Repression" if r.is_inhibition else "Unknown")
        )
        rows.append(
            dict(
                tf=tf,
                target=r.target_genesymbol,
                axis=AXIS[tf],
                resource="DoRothEA",
                confidence=str(r.dorothea_level),
                effect=eff,
                reference=str(r.references)[:300],
                note="DoRothEA via OmniPath; not a lung ChIP peak",
            )
        )

    col_path = cache / "collectri.tsv"
    col_raw = col_path.read_bytes() if col_path.exists() else fetch(URLS["collectri"])
    col = pd.read_csv(pd.io.common.BytesIO(col_raw), sep="\t")
    for r in col.itertuples(index=False):
        tf = r.source_genesymbol
        if tf not in FOCUS:
            continue
        eff = (
            "Activation"
            if r.is_stimulation
            else ("Repression" if r.is_inhibition else "Unknown")
        )
        rows.append(
            dict(
                tf=tf,
                target=r.target_genesymbol,
                axis=AXIS[tf],
                resource="CollecTRI",
                confidence="",
                effect=eff,
                reference=str(r.references)[:300],
                note="CollecTRI via OmniPath; not a lung ChIP peak",
            )
        )

    df = pd.DataFrame(rows).drop_duplicates(["tf", "target", "resource"])
    df = df.sort_values(["axis", "tf", "resource", "target"])
    df.to_csv(OUT / "tf_targets_public.tsv", sep="\t", index=False)
    listed = sorted(set(df.loc[df.target.isin(["TACSTD2", "CLDN4"]), "tf"]))
    manifest = {
        "purpose": "Public curated TF-target edges for IFN/MHC/TJ TFs. Not ChIP.",
        "do_not_invent_chip": True,
        "a10_elf3_cldn4_given": True,
        "urls": URLS,
        "n_rows": int(len(df)),
        "n_tfs": int(df.tf.nunique()),
        "tfs": FOCUS,
        "tacstd2_cldn4_listed_as_targets_of": listed,
    }
    (OUT / "prior_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(df.groupby(["axis", "tf"]).size().to_string())
    print("wrote", OUT / "tf_targets_public.tsv", "n=", len(df))


if __name__ == "__main__":
    main()
