#!/usr/bin/env python3
"""Score JUMP CLDN4 nuclear/DNA features against controls.

Uses the interpretable (pre-batch-correction) well profiles, which still have
CellProfiler feature names. RxRx3 embeddings are not scored: they have no
DNA or nuclear feature names.

Two contrasts:
  plate_negcon  CLDN4 wells vs documented negative controls on the same plates
  same_well     CLDN4 vs other genes in the same well position
                (CRISPR is always G18; ORF is always B05)

The second contrast is the specificity check. A plate-negcon difference that
is shared by other genes in the same well is a well-position offset.
"""

from __future__ import annotations

import csv
from pathlib import Path

import duckdb
import numpy as np
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "cache" / "jump"
OUT = ROOT / "results" / "tables"

FEATURES = [
    "Nuclei_AreaShape_Area",
    "Nuclei_AreaShape_Compactness",
    "Nuclei_AreaShape_Eccentricity",
    "Nuclei_AreaShape_Solidity",
    "Nuclei_AreaShape_FormFactor",
    "Nuclei_Intensity_MeanIntensity_DNA",
    "Nuclei_Intensity_MedianIntensity_DNA",
    "Nuclei_Intensity_MADIntensity_DNA",
    "Nuclei_Intensity_IntegratedIntensity_DNA",
    "Nuclei_Intensity_StdIntensity_DNA",
    "Nuclei_Texture_Contrast_DNA_5_00_256",
    "Nuclei_Texture_Entropy_DNA_5_00_256",
    "Nuclei_Texture_Variance_DNA_5_00_256",
    "Nuclei_Texture_InfoMeas1_DNA_5_00_256",
    "Nuclei_Granularity_1_DNA",
    "Nuclei_Granularity_3_DNA",
    "Cells_Intensity_MeanIntensity_DNA",
    "Cells_Intensity_IntegratedIntensity_DNA",
    "Cytoplasm_Intensity_MeanIntensity_DNA",
]

PROFILES = {
    "crispr": "https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump-assembled/source_all/workspace/profiles_assembled/CRISPR/v1.0a/profiles_wellpos_cc_var_mad_outlier.parquet",
    "orf": "https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump-assembled/source_all/workspace/profiles_assembled/ORF/v1.0a/profiles_wellpos_cc_var_mad_outlier.parquet",
}
SPECS = {
    "crispr": {
        "source": "source_13",
        "well": "G18",
        "trt": "JCP2022_801379",
        "neg": ("JCP2022_800001", "JCP2022_800002"),
        "plates": (
            "CP-CC9-R1-06",
            "CP-CC9-R2-06",
            "CP-CC9-R3-06",
            "CP-CC9-R4-06",
            "CP-CC9-R5-06",
        ),
    },
    "orf": {
        "source": "source_4",
        "well": "B05",
        "trt": "JCP2022_905637",
        "neg": (
            "JCP2022_915128",
            "JCP2022_915129",
            "JCP2022_915130",
            "JCP2022_915131",
        ),
        "plates": (
            "BR00123512",
            "BR00123513",
            "BR00123514",
            "BR00123515",
            "BR00123516",
        ),
    },
}
ZENODO = {
    "crispr": "https://zenodo.org/api/records/14861664/files/crispr_interpretable_features.parquet/content",
    "orf": "https://zenodo.org/api/records/14861664/files/orf_interpretable_features.parquet/content",
}
SIG_URL = {
    "crispr": "https://zenodo.org/api/records/14861664/files/crispr_interpretable_significance_full.parquet/content",
    "orf": "https://zenodo.org/api/records/14861664/files/orf_interpretable_significance_full.parquet/content",
}


def connect() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute("INSTALL httpfs; LOAD httpfs;")
    return con


def fetch_zenodo(con: duckdb.DuckDBPyConnection, modality: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / f"{modality}_interpretable_features.parquet"
    if not dest.exists() or dest.stat().st_size == 0:
        import urllib.request

        print("download", dest.name)
        urllib.request.urlretrieve(ZENODO[modality], dest)
    return dest


def pull(con: duckdb.DuckDBPyConnection, modality: str) -> tuple[Path, Path]:
    spec = SPECS[modality]
    cols = ", ".join(
        ["Metadata_Source", "Metadata_Plate", "Metadata_Well", "Metadata_JCP2022"]
        + [f'"{c}"' for c in FEATURES]
    )
    same = CACHE / f"{modality}_samewell_panel.parquet"
    plates = CACHE / f"{modality}_plates_panel.parquet"
    if not same.exists():
        print("pull same-well", modality)
        con.execute(
            f"""
            COPY (
              SELECT {cols}
              FROM read_parquet('{PROFILES[modality]}')
              WHERE Metadata_Source = '{spec["source"]}'
                AND Metadata_Well = '{spec["well"]}'
            ) TO '{same}' (FORMAT PARQUET)
            """
        )
    if not plates.exists():
        plist = ", ".join(f"'{p}'" for p in spec["plates"])
        print("pull plates", modality)
        con.execute(
            f"""
            COPY (
              SELECT {cols}
              FROM read_parquet('{PROFILES[modality]}')
              WHERE Metadata_Source = '{spec["source"]}'
                AND Metadata_Plate IN ({plist})
            ) TO '{plates}' (FORMAT PARQUET)
            """
        )
    return same, plates


def mw(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return float("nan")
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def rows_for(df, trt_id: str, control_ids: set[str] | None, features: list[str], contrast: str, modality: str):
    trt = df[df.Metadata_JCP2022 == trt_id]
    if control_ids is None:
        ctrl = df[df.Metadata_JCP2022 != trt_id]
    else:
        ctrl = df[df.Metadata_JCP2022.isin(control_ids)]
    alpha = 0.05 / len(features)
    out = []
    for feature in features:
        x = trt[feature].to_numpy(float)
        y = ctrl[feature].to_numpy(float)
        p = mw(x, y)
        mx = float(np.nanmedian(x))
        my = float(np.nanmedian(y))
        out.append(
            {
                "modality": modality,
                "contrast": contrast,
                "feature": feature,
                "n_cldn4": int(np.isfinite(x).sum()),
                "n_control": int(np.isfinite(y).sum()),
                "median_cldn4": mx,
                "median_control": my,
                "delta": mx - my,
                "mannwhitney_p": p,
                "bonferroni_alpha": alpha,
                "survives_bonferroni": bool(p < alpha),
            }
        )
    return out


def channel_of(name: str) -> str:
    for label in ("DNA", "RNA", "ER", "AGP", "Mito"):
        if name.endswith("_" + label) or f"_{label}_" in name:
            return label
        if label == "Mito" and "mito" in name:
            return label
    return "other"


def significance_floor(con: duckdb.DuckDBPyConnection, modality: str) -> dict:
    import urllib.request
    from collections import Counter

    dest = CACHE / f"{modality}_interpretable_significance_full.parquet"
    if not dest.exists() or dest.stat().st_size == 0:
        print("download", dest.name)
        urllib.request.urlretrieve(SIG_URL[modality], dest)
    row = con.execute(
        f"""
        SELECT * EXCLUDE (Metadata_JCP2022)
        FROM read_parquet('{dest}')
        WHERE Metadata_JCP2022 = '{SPECS[modality]["trt"]}'
        """
    ).fetchdf()
    values = row.iloc[0].astype(float)
    floor = float(values.min())
    tied = values[values <= floor + 1e-9]
    counts = Counter(channel_of(name) for name in tied.index)
    nuclei_dna_intensity = [
        name
        for name in tied.index
        if name.startswith("Nuclei_") and "Intensity" in name and "DNA" in name
    ]
    return {
        "modality": modality,
        "min_feature_significance": floor,
        "n_features": int(len(values)),
        "n_tied_at_floor": int(len(tied)),
        "tied_DNA": counts["DNA"],
        "tied_RNA": counts["RNA"],
        "tied_Mito": counts["Mito"],
        "tied_ER": counts["ER"],
        "tied_AGP": counts["AGP"],
        "tied_other": counts["other"],
        "nuclei_dna_intensity_at_floor": len(nuclei_dna_intensity),
    }


def official(con: duckdb.DuckDBPyConnection, modality: str) -> dict:
    path = fetch_zenodo(con, modality)
    df = con.execute(
        f"""
        SELECT "Corrected p-value" AS q, "Phenotypic activity" AS activity, COUNT(*) AS n
        FROM read_parquet('{path}')
        WHERE "JCP2022 ID" = '{SPECS[modality]["trt"]}'
        GROUP BY 1, 2
        """
    ).fetchdf()
    if len(df) != 1:
        raise SystemExit(f"{modality} official rows: {len(df)}")
    return {
        "modality": modality,
        "jcp": SPECS[modality]["trt"],
        "corrected_p": float(df.q.iloc[0]),
        "phenotypic_activity": float(df.activity.iloc[0]),
        "feature_rows_in_browser_table": int(df.n.iloc[0]),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    con = connect()
    scored = []
    official_rows = []
    floors = []
    for modality in ("crispr", "orf"):
        official_rows.append(official(con, modality))
        floors.append(significance_floor(con, modality))
        same_path, plate_path = pull(con, modality)
        same = con.execute(f"SELECT * FROM read_parquet('{same_path}')").fetchdf()
        plates = con.execute(f"SELECT * FROM read_parquet('{plate_path}')").fetchdf()
        spec = SPECS[modality]
        scored.extend(
            rows_for(same, spec["trt"], None, FEATURES, "same_well", modality)
        )
        scored.extend(
            rows_for(
                plates,
                spec["trt"],
                set(spec["neg"]),
                FEATURES,
                "plate_negcon",
                modality,
            )
        )
        # document the negcon offset on one shape feature and one DNA-intensity feature
        neg = plates[plates.Metadata_JCP2022.isin(spec["neg"])]
        print(
            modality,
            "negcon wells",
            len(neg),
            "same-well others",
            int((same.Metadata_JCP2022 != spec["trt"]).sum()),
        )

    score_path = OUT / "jump_cldn4_nuclear_proxies.tsv"
    with score_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scored[0].keys()))
        writer.writeheader()
        writer.writerows(scored)

    off_path = OUT / "jump_cldn4_official_phenotype.tsv"
    with off_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(official_rows[0].keys()))
        writer.writeheader()
        writer.writerows(official_rows)

    floor_path = OUT / "jump_cldn4_significance_floor.tsv"
    with floor_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(floors[0].keys()))
        writer.writeheader()
        writer.writerows(floors)
    print("wrote", floor_path)

    survivors = [row for row in scored if row["contrast"] == "same_well" and row["survives_bonferroni"]]
    print("same-well Bonferroni survivors", len(survivors))
    for row in survivors:
        print(row["modality"], row["feature"], row["mannwhitney_p"])
    print("wrote", score_path)
    print("wrote", off_path)


if __name__ == "__main__":
    main()
