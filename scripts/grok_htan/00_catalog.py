#!/usr/bin/env python3
"""Catalog HTAN/HCA open lung processed objects on CELLxGENE + note skipped sources."""

from __future__ import annotations

import json
import urllib.request
from datetime import datetime, timezone

from lib_htan import ensure_dirs

CXG = "https://api.cellxgene.cziscience.com/curation/v1/collections"

# Verified CELLxGENE collection IDs (queried live; do not invent).
COLLECTIONS = {
    "62e8f058-9c37-48bc-9200-e767f318a8ec": {
        "label": "HTAN_MSK_SCLC_Chan2021",
        "source": "HTAN",
        "paper": "Chan et al. Cancer Cell 2021; 10.1016/j.ccell.2021.09.008",
    },
    "efd94500-1fdc-4e28-9e9f-a309d0154e21": {
        "label": "HTAN_MSK_Treg_Glasner2023",
        "source": "HTAN",
        "paper": "Glasner et al. Nat Immunol 2023; 10.1038/s41590-023-01504-2",
    },
    "6f6d381a-7701-4781-935c-db10d30de293": {
        "label": "HCA_HLCA_Sikkema2023",
        "source": "HCA",
        "paper": "Sikkema et al. Nat Med 2023; 10.1038/s41591-023-02327-2",
    },
    "5d445965-6f1a-4b68-ba3a-b8f765155d3a": {
        "label": "HCA_Travaglini2020",
        "source": "HCA",
        "paper": "Travaglini et al. Nature 2020; 10.1038/s41586-020-2922-4",
    },
    "4d74781b-8186-4c9a-b659-ff4dc4601d91": {
        "label": "HCA_Madissoon2020_stability",
        "source": "HCA",
        "paper": "Madissoon et al. Genome Biol 2020; 10.1186/s13059-019-1906-x",
    },
    "c1241244-b22d-483d-875b-75699efb9f3c": {
        "label": "HCA_Madissoon2023_spatial",
        "source": "HCA",
        "paper": "Madissoon et al. Nat Genet 2023; 10.1038/s41588-022-01243-4",
    },
    "0bebef1a-4607-4584-9070-dacf89a0d635": {
        "label": "LUAD_histologic_subtypes",
        "source": "open_CXG_lung",
        "paper": "10.1186/s40164-025-00740-6",
    },
    "edb893ee-4066-4128-9aec-5eb2b03f8287": {
        "label": "LuCA_Salcher2022",
        "source": "open_CXG_lung",
        "paper": "Salcher et al. Cancer Cell 2022; 10.1016/j.ccell.2022.10.008",
    },
}

# Manual skip / access notes (verified on official pages).
MANUAL_ROWS = [
    {
        "collection_id": "phs002371",
        "dataset_id": "phs002371",
        "title": "HTAN controlled-access raw sequencing (dbGaP phs002371)",
        "source": "HTAN",
        "assay": "scRNA/spatial raw FASTQ/BAM",
        "tissue": "multi-atlas including lung",
        "cells": "",
        "size_bytes": "",
        "size_gb": "",
        "flag": "SKIP_CONTROLLED_RAW",
        "decision": "skip",
        "url": "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002371.v1.p1",
        "notes": "User instruction: skip dbGaP phs002371 raw. HTAN Level 1/2 sequencing is controlled via this study. Do not download or fabricate.",
        "paper": "HTAN consortium data access",
    },
    {
        "collection_id": "10.5281/zenodo.16546233",
        "dataset_id": "zenodo.19559522",
        "title": "TsankovLab HTAN_Lung processed atlas zip (TP53 LUAD spatial/sc)",
        "source": "HTAN_related",
        "assay": "scRNA + spatial bundle",
        "tissue": "lung adenocarcinoma",
        "cells": "",
        "size_bytes": 18858121175,
        "size_gb": 18.86,
        "flag": "SKIP_GT_2GB",
        "decision": "skip",
        "url": "https://zenodo.org/records/19559522",
        "notes": "Open CC-BY but single zip 18.86GB >2GB cap. H&E zip 30.7GB also skipped.",
        "paper": "Zhao et al. HTAN_Lung; 10.5281/zenodo.16546233",
    },
    {
        "collection_id": "a48f5033-3438-4550-8574-cdff3263fdfd",
        "dataset_id": "HTAN_VUMC_CXG",
        "title": "HTAN VUMC CELLxGENE collection is colorectal polyps, not lung",
        "source": "HTAN",
        "assay": "scRNA",
        "tissue": "colon/rectum (NOT lung)",
        "cells": "",
        "size_bytes": "",
        "size_gb": "",
        "flag": "OUT_OF_SCOPE",
        "decision": "skip",
        "url": "https://cellxgene.cziscience.com/collections/a48f5033-3438-4550-8574-cdff3263fdfd",
        "notes": "Title resembles Sinjab/Kadara lung Cell 2021 but CXG objects are colorectal. Not used.",
        "paper": "VUMC HTAN colorectal",
    },
]


def get_json(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)


def main():
    paths = ensure_dirs()
    rows = list(MANUAL_ROWS)
    for cid, meta in COLLECTIONS.items():
        c = get_json(f"{CXG}/{cid}")
        for d in c.get("datasets") or []:
            assets = d.get("assets") or []
            h5 = next((a for a in assets if str(a.get("filetype", "")).upper() == "H5AD"), None)
            size = h5.get("filesize") if h5 else None
            url = h5.get("url") if h5 else ""
            assays = ";".join(a.get("label", "") for a in (d.get("assay") or []))
            tissues = ";".join(t.get("label", "") for t in (d.get("tissue") or []))
            if size is None:
                flag = "NOSIZE"
                decision = "review"
            elif size >= 2_000_000_000:
                flag = "SKIP_GT_2GB"
                decision = "skip"
            else:
                flag = "OPEN_LT_2GB"
                decision = "eligible"
            rows.append(
                {
                    "collection_id": cid,
                    "dataset_id": d.get("dataset_id"),
                    "title": d.get("title"),
                    "source": meta["source"],
                    "assay": assays,
                    "tissue": tissues,
                    "cells": d.get("cell_count"),
                    "size_bytes": size if size is not None else "",
                    "size_gb": round(size / 1e9, 3) if size else "",
                    "flag": flag,
                    "decision": decision,
                    "url": url,
                    "notes": meta["label"],
                    "paper": meta["paper"],
                }
            )

    import pandas as pd

    df = pd.DataFrame(rows)
    out = paths["notes"] / "catalog.tsv"
    df.to_csv(out, sep="\t", index=False)
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_rows": int(len(df)),
        "n_eligible": int((df["decision"] == "eligible").sum()),
        "n_skip_2gb": int((df["flag"] == "SKIP_GT_2GB").sum()),
        "n_skip_controlled": int((df["flag"] == "SKIP_CONTROLLED_RAW").sum()),
        "catalog": str(out),
    }
    (paths["notes"] / "catalog_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(df["flag"].value_counts().to_string())


if __name__ == "__main__":
    main()
