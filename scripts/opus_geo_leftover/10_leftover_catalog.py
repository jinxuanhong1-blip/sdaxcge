#!/usr/bin/env python3
"""Disposition table for leftover human lung ICI GEO series.

Records, for every probed lung+ICI expression series, whether it was analysed in this
slice, already finished in a sibling PR, or skipped — and why. Accessions are the live
GEO identifiers from the Entrez/FTP probe, not a hand-written list.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"

EXCLUDED_NAMED = {
    "GSE126044",
    "GSE135222",
    "GSE136961",
    "GSE166449",
    "GSE93157",
    "GSE207422",
    "GSE205335",
}

# Sibling PRs that already finished a TACSTD2/CLDN4 vs-response (or documented-null) analysis
# of that accession. Used only as a disposition label; this slice still independently verified
# the accession against Entrez.
SIBLING = {
    "GSE126044": "PR 4 / 5 / 11 / 36 (named exclusion)",
    "GSE135222": "PR 87 / 4 / 11 / 36 (named exclusion)",
    "GSE136961": "PR 4 / 36 (named exclusion; Oncomine panel lacks genes)",
    "GSE166449": "PR 4 / 11 (named exclusion)",
    "GSE93157": "PR 76 (named exclusion; NanoString panel lacks genes)",
    "GSE207422": "PR 51 / 11 (named exclusion)",
    "GSE205335": "PR 50 (named exclusion; scRNA)",
    "GSE182328": "PR 36 (Akkermansia surrogate, not RECIST)",
    "GSE111414": "PR 36 / 115 (PBMC CD8; genes ~0)",
    "GSE161537": "PR 51 (HTG 2560-gene panel; TACSTD2/CLDN4 absent)",
    "GSE162520": "PR 51 (HTG 2560-gene panel; TACSTD2/CLDN4 absent; ICI-naive comparator)",
    "GSE110390": "PR 76 / 11 (21-gene IFN-γ panel; genes absent)",
    "GSE216297": "PR 115 (platelet RNA; genes filtered from public matrix)",
    "GSE305086": "PR 45 / 115 (whole blood array; TACSTD2 at background)",
    "GSE285888": "PR 45 / 115 (PBMC scRNA)",
    "GSE249262": "PR 115 (CTC)",
    "GSE202417": "PR 115 (CD8 PBMC array)",
    "GSE261345": "PR 35 (GeoMx CTA; TACSTD2 only)",
    "GSE261348": "PR 35 (GeoMx CTA; TACSTD2 only)",
    "GSE233203": "PR 35 (scRNA, n=7)",
    "GSE91061": "PR 76 (melanoma, not lung)",
    "GSE72094": "PR 76 (LUAD survival, not ICI-treated)",
    "GSE81089": "PR 76 (NSCLC survival, not ICI-treated)",
}

ANALYZED_HERE = {
    "GSE253564": "pre-treatment bulk FPKM vs MPR / recurrence / PFS (this slice)",
    "GSE248378": "post-treatment bulk FPKM vs recurrence / PFS; MPR not estimable (0 MPR cases deposited)",
}


def main() -> int:
    probe_path = OUT / "geo_probe.json"
    if not probe_path.exists():
        probe_path = OUT / "geo_probe_slim.json"
    raw = json.loads(probe_path.read_text())
    probe = []
    for r in raw:
        if "characteristics" not in r:
            r = {
                **r,
                "characteristics": [""] * int(r.get("n_characteristics") or 0),
            }
        probe.append(r)
    inv = list(csv.DictReader((OUT / "lung_ici_response_inventory.csv").open()))
    inv_acc = {r["accession"] for r in inv}
    manifest = {r["accession"]: r for r in json.loads((OUT / "download_manifest.json").read_text())}

    rows = []
    for r in probe:
        acc = r["accession"]
        title = r["title"]
        n = r.get("n_characteristics") or max(
            (len(c.split("\t")) - 1 for c in r.get("characteristics") or []), default=0
        )
        if acc in EXCLUDED_NAMED:
            disp, why = "named_exclusion", "on the seven-accession exclusion list"
        elif acc in ANALYZED_HERE:
            disp, why = "analyzed_here", ANALYZED_HERE[acc]
        elif acc in SIBLING:
            disp, why = "already_finished_elsewhere", SIBLING[acc]
        elif acc in inv_acc:
            disp, why = "leftover_not_analyzed", "has an outcome field but not a usable leftover tumour matrix + genes"
        else:
            continue  # not lung+ICI+outcome; keep the table focused

        extra = ""
        if acc == "GSE126045":
            extra = "SuperSeries of excluded GSE126044; not an independent cohort"
            disp, why = "named_exclusion_superseries", extra
        if acc == "GSE309652":
            extra = "NanoString NS_Hs_Metabolism_v1.0 (~768 genes); TACSTD2/CLDN4 absent from RCC files"
            disp, why = "leftover_genes_absent", extra
        if acc == "GSE248249":
            extra = "Clariom D probe IDs; GPL family SOFT 4.8 GB (>2 GB skip); no response field in GEO"
            disp, why = "leftover_skipped", extra
        if acc == "GSE221733":
            extra = "GeoMx CTA spatial; TACSTD2/CLDN4 not on the 68-gene CTA panel"
            disp, why = "leftover_genes_absent", extra

        rec = {
            "accession": acc,
            "disposition": disp,
            "reason": why,
            "n_samples_metadata": n,
            "title": title,
            "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
            "entrez_verified": acc in manifest and manifest[acc].get("verified") is True,
            "pubmed": ";".join(r.get("pubmed") or []),
        }
        rows.append(rec)

    # Guarantee the analysed + named + sibling accessions appear even if probe filtering dropped them.
    have = {r["accession"] for r in rows}
    for acc, why in {**{a: "named exclusion" for a in EXCLUDED_NAMED}, **SIBLING, **ANALYZED_HERE}.items():
        if acc not in have:
            rows.append(
                {
                    "accession": acc,
                    "disposition": (
                        "analyzed_here"
                        if acc in ANALYZED_HERE
                        else "named_exclusion"
                        if acc in EXCLUDED_NAMED
                        else "already_finished_elsewhere"
                    ),
                    "reason": ANALYZED_HERE.get(acc) or SIBLING.get(acc) or why,
                    "n_samples_metadata": "",
                    "title": "",
                    "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
                    "entrez_verified": acc in manifest and manifest[acc].get("verified") is True,
                    "pubmed": "",
                }
            )

    rows.sort(key=lambda r: (r["disposition"], r["accession"]))
    path = OUT / "leftover_disposition.tsv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"[catalog] {len(rows)} rows -> {path.name}")
    from collections import Counter

    print(Counter(r["disposition"] for r in rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
