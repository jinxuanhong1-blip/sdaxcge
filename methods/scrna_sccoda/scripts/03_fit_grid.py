#!/usr/bin/env python3
"""Combinatorial cohort × annotation × contrast composition models.

Primary question: which cohort × annotation recovers T/NK or TLS *down*
in TACSTD2-high. MPR is a separate family. FDR is BH within family;
n is patients, never cells.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lib.dirichlet_multinomial import (  # noqa: E402
    alr_lm,
    fit_dm_glm,
    fraction_tests,
    permutation_p_alr,
    try_sccoda,
)
from lib.fdr import add_fdr, primary_recovery_flag  # noqa: E402

MIN_SAMPLES = 6
MIN_PER_ARM = 2
MIN_CELLS = 50
DROP_RARE_FRAC = 0.005  # drop types with mean fraction < 0.5%


def wide_counts(comp: pd.DataFrame, cohort: str, annotation: str) -> pd.DataFrame:
    sub = comp[(comp["cohort"] == cohort) & (comp["annotation"] == annotation)]
    if sub.empty:
        return pd.DataFrame()
    wide = sub.pivot_table(index="sample_id", columns="celltype", values="n", aggfunc="sum", fill_value=0)
    wide = wide.loc[wide.sum(axis=1) >= MIN_CELLS]
    mean_frac = wide.div(wide.sum(axis=1), axis=0).mean()
    keep = mean_frac[mean_frac >= DROP_RARE_FRAC].index.tolist()
    if len(keep) < 3:
        keep = mean_frac.sort_values(ascending=False).head(4).index.tolist()
    return wide[keep]


def pick_reference(wide: pd.DataFrame, prefer: tuple[str, ...] = ("stromal", "myeloid", "other", "mast")) -> str:
    props = wide.div(wide.sum(axis=1), axis=0)
    mean_p = props.mean()
    sd_p = props.std(ddof=1)
    for name in prefer:
        if name in wide.columns and mean_p[name] >= 0.02:
            return name
    eligible = mean_p[mean_p >= 0.02]
    if eligible.empty:
        eligible = mean_p
    return str(sd_p.loc[eligible.index].idxmin())


def contrast_x(cov: pd.DataFrame, samples: list[str], contrast: str) -> tuple[np.ndarray, pd.Index, dict]:
    c = cov.set_index("sample_id").reindex(samples)
    info = {"contrast": contrast, "n_input": len(samples)}
    if contrast == "tacstd2_high":
        ok = c["eligible_tacstd2"].fillna(0).astype(bool) & c["tacstd2_high"].notna()
        x = c.loc[ok, "tacstd2_high"].astype(float)
        info.update(
            {
                "n_high": int((x == 1).sum()),
                "n_low": int((x == 0).sum()),
                "binary": True,
                "fdr_family": "primary_tacstd2",
            }
        )
        return x.to_numpy(), x.index, info
    if contrast == "tacstd2_z":
        ok = c["eligible_tacstd2"].fillna(0).astype(bool) & c["tacstd2_z"].notna()
        x = c.loc[ok, "tacstd2_z"].astype(float)
        info.update({"n": int(ok.sum()), "binary": False, "fdr_family": "secondary_tacstd2_continuous"})
        return x.to_numpy(), x.index, info
    if contrast == "mpr":
        ok = c["eligible_mpr"].fillna(0).astype(bool) & c["mpr"].isin(["MPR", "NMPR"])
        x = (c.loc[ok, "mpr"] == "MPR").astype(float)
        info.update(
            {
                "n_mpr": int((x == 1).sum()),
                "n_nmpr": int((x == 0).sum()),
                "binary": True,
                "fdr_family": "mpr",
            }
        )
        return x.to_numpy(), x.index, info
    raise ValueError(contrast)


def fit_one(wide: pd.DataFrame, x: np.ndarray, used: pd.Index, info: dict, cohort: str, annotation: str) -> list[dict]:
    wide = wide.loc[used]
    if len(wide) < MIN_SAMPLES:
        return []
    if info.get("binary") and (info.get("n_high", info.get("n_mpr", 0)) < MIN_PER_ARM or info.get("n_low", info.get("n_nmpr", 0)) < MIN_PER_ARM):
        return []
    # drop all-zero columns in this subset
    wide = wide.loc[:, wide.sum(axis=0) > 0]
    if wide.shape[1] < 3:
        return []
    ref = pick_reference(wide)
    counts = wide.to_numpy(dtype=float)
    celltypes = list(wide.columns)
    rows = []
    dm = fit_dm_glm(counts, x, celltypes, reference=ref)
    for i, ct in enumerate(celltypes):
        compartment = ct if ct in {"TNK", "TLS"} else ("TNK" if ct in {"T", "NK"} else ("TLS" if ct in {"B", "plasma"} else ct))
        # DM LRT p-values scale with library size; keep as diagnostic only.
        if info["contrast"] == "tacstd2_high":
            fam = "diagnostic_dm_lrt"
        else:
            fam = "diagnostic_dm_lrt_" + info["fdr_family"]
        rows.append(
            {
                "cohort": cohort,
                "annotation": annotation,
                "contrast": info["contrast"],
                "celltype": ct,
                "compartment": compartment,
                "method": "dm_glm",
                "effect": float(dm.slope[i]),
                "se": float(dm.se_slope[i]) if np.isfinite(dm.se_slope[i]) else np.nan,
                "p_value": float(dm.p_slope[i]) if np.isfinite(dm.p_slope[i]) else np.nan,
                "n": dm.n_samples,
                "n_celltypes": dm.n_celltypes,
                "reference": dm.reference,
                "precision": dm.precision,
                "converged": dm.converged,
                "fdr_family": fam,
                "note": dm.note,
                **{k: v for k, v in info.items() if k not in {"contrast", "fdr_family"}},
            }
        )
    ref_i = celltypes.index(ref)
    for rec in alr_lm(counts, x, celltypes, reference=ref):
        ct = rec["celltype"]
        compartment = ct if ct in {"TNK", "TLS"} else ("TNK" if ct in {"T", "NK"} else ("TLS" if ct in {"B", "plasma"} else ct))
        if info["contrast"] == "tacstd2_high" and compartment in {"TNK", "TLS"} and ct != ref:
            fam = "primary_tacstd2"
            seed = int(hashlib.md5(f"{cohort}|{annotation}|{ct}".encode()).hexdigest()[:8], 16)
            p_perm = permutation_p_alr(counts, x, ref_i, celltypes.index(ct), seed=seed)
            rec["p_value"] = p_perm
            rec["note"] = (rec.get("note") or "") + ";p_patient_perm"
        elif info["contrast"] == "tacstd2_high":
            fam = "secondary_tacstd2_other_types"
            p_perm = np.nan
        else:
            fam = "companion_alr"
            p_perm = np.nan
        rec.update(
            {
                "cohort": cohort,
                "annotation": annotation,
                "contrast": info["contrast"],
                "compartment": compartment,
                "fdr_family": fam,
                "reference": ref,
                "converged": True,
                "p_perm": p_perm,
                "method": "alr_lm" if fam == "companion_alr" else "alr_perm",
            }
        )
        rows.append(rec)
    for rec in fraction_tests(counts, x, celltypes, binary=bool(info.get("binary"))):
        rec.update(
            {
                "cohort": cohort,
                "annotation": annotation,
                "contrast": info["contrast"],
                "compartment": rec["celltype"],
                "fdr_family": "companion_fraction_naive",
                "reference": "",
                "converged": True,
            }
        )
        rows.append(rec)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args()
    res = args.root / "results"
    comp = pd.read_csv(res / "composition_long.tsv", sep="\t")
    cov = pd.read_csv(res / "sample_covariates.tsv", sep="\t")

    pairs = sorted(comp.groupby(["cohort", "annotation"]).size().index.tolist())
    contrasts = ["tacstd2_high", "tacstd2_z", "mpr"]
    all_rows = []
    sccoda_status = {}
    for cohort, annotation in pairs:
        wide = wide_counts(comp, cohort, annotation)
        if wide.empty:
            continue
        cov_c = cov[cov["cohort"] == cohort]
        for contrast in contrasts:
            x, used, info = contrast_x(cov_c, list(wide.index), contrast)
            rows = fit_one(wide, x, used, info, cohort, annotation)
            all_rows.extend(rows)
        # optional scCODA on the primary binary TACSTD2 split only
        x, used, info = contrast_x(cov_c, list(wide.index), "tacstd2_high")
        if len(used) >= MIN_SAMPLES:
            w = wide.loc[used]
            covariates = pd.DataFrame({"tacstd2_high": x}, index=used)
            sccoda_status[f"{cohort}|{annotation}"] = try_sccoda(
                w.to_numpy(), covariates, list(w.columns), "C(tacstd2_high)", pick_reference(w)
            )

    df = pd.DataFrame(all_rows)
    if df.empty:
        raise SystemExit("no models fit")
    df = add_fdr(df)
    # Bonferroni on the primary DM-GLM TNK/TLS tacstd2_high grid
    prim = df[
        (df["method"] == "alr_perm")
        & (df["contrast"] == "tacstd2_high")
        & (df["compartment"].isin(["TNK", "TLS", "T", "NK", "B", "plasma"]))
    ]
    n_prim = int(prim["p_value"].notna().sum())
    df["n_primary_grid"] = n_prim
    df["bonferroni_primary"] = np.nan
    if n_prim:
        mask = df.index.isin(prim.index) & df["p_value"].notna()
        df.loc[mask, "bonferroni_primary"] = np.clip(df.loc[mask, "p_value"] * n_prim, 0, 1)
    df["recovery"] = df.apply(primary_recovery_flag, axis=1)
    # recovery uses BH within primary_tacstd2 family for DM-GLM only
    rec_mask = (df["method"] == "alr_perm") & (df["contrast"] == "tacstd2_high")
    df.loc[~rec_mask, "recovery"] = np.where(df.loc[~rec_mask, "recovery"] == "not_primary", "not_primary", df.loc[~rec_mask, "recovery"])

    df.to_csv(res / "grid_effects.tsv", sep="\t", index=False)

    # recovery table: one row per cohort × annotation × compartment
    rec = df[(df["method"] == "alr_perm") & (df["contrast"] == "tacstd2_high")].copy()
    rec = rec[rec["compartment"].isin(["TNK", "TLS", "T", "NK", "B", "plasma"])]
    rec.to_csv(res / "recovery_tnk_tls.tsv", sep="\t", index=False)

    recovered = rec[rec["recovery"] == "recovered_down"]
    summary = {
        "n_models_rows": int(len(df)),
        "n_primary_grid_tests": n_prim,
        "n_recovered_tnk_or_tls_down": int(len(recovered)),
        "primary_pvalue": "patient-level permutation of the ALR slope (does not scale with n_cells)",
        "recovered": recovered[["cohort", "annotation", "celltype", "compartment", "effect", "p_value", "q_bh", "n"]].to_dict("records"),
        "sccoda": {
            "status": "not_available",
            "n_pairs": len(sccoda_status),
            "reason": next(iter(sccoda_status.values()), {}).get("reason", ""),
        },
        "honest_n_note": "n is patients/samples. Cells are the compositional library size, not the replicate unit.",
        "fdr": {
            "primary_family": "BH over patient-permutation ALR tests labeled primary_tacstd2 (TNK/TLS)",
            "q_cut_recovery": 0.10,
            "bonferroni_column": "bonferroni_primary on the T/NK/B/TLS tacstd2_high ALR-perm grid",
            "dm_lrt": "diagnostic only; not used for FDR (library-size inflated)",
        },
    }
    (res / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps({k: summary[k] for k in summary if k != "recovered"}, indent=2))
    print("recovered rows:", len(recovered))


if __name__ == "__main__":
    main()
