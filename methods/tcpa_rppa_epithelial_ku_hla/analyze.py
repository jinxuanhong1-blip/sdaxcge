#!/usr/bin/env python3
"""TCPA RPPA: epithelial antibodies vs Ku80 / HLA-DQA1 in LUAD and LUSC.

CLDN4, DNA-PKcs, Ku70, STAT1, and HLA class I are tested only if an antibody
with that identity is on the disease level-4 panel. They are not replaced by
STAT3, IRF1, XRCC1, or any other neighbor.

Primary family (per histology, then BH within that family):
  Claudin-7, E-cadherin, EMA  vs  Ku80 and HLA-DQA1.
Sensitivity (not in the BH family):
  within-cohort z-mean of the three epithelial antibodies vs the same two targets.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# TCPA protein_id -> (gene symbol, plain name, role)
ANTIBODY = {
    "CLAUDIN7": ("CLDN7", "Claudin-7", "epithelial"),
    "ECADHERIN": ("CDH1", "E-cadherin", "epithelial"),
    "EMA": ("MUC1", "EMA", "epithelial"),
    "KU80": ("XRCC5", "Ku80", "target"),
    "HLADQA1": ("HLA-DQA1", "HLA-DQA1", "target"),
}

EPITHELIAL = ["CLAUDIN7", "ECADHERIN", "EMA"]
TARGETS = ["KU80", "HLADQA1"]

# Identity search on the TCPA protein_id list. Neighbors are not aliases.
ABSENCE_QUERIES = [
    {
        "requested": "CLDN4",
        "label": "Claudin-4 / CLDN4",
        "exact": ["CLDN4", "CLAUDIN4", "CLAUDIN-4"],
    },
    {
        "requested": "PRKDC",
        "label": "DNA-PKcs / PRKDC",
        "exact": ["PRKDC", "DNAPK", "DNAPKCS", "DNAPKCS", "DNA-PKCS", "DNAPKCSER2056"],
    },
    {
        "requested": "XRCC6",
        "label": "Ku70 / XRCC6",
        "exact": ["KU70", "XRCC6"],
    },
    {
        "requested": "STAT1",
        "label": "STAT1",
        "exact": ["STAT1", "STAT1PY701", "STAT1PS727"],
    },
    {
        "requested": "HLA-I",
        "label": "HLA class I (HLA-A/B/C or B2M)",
        "exact": ["HLAA", "HLAB", "HLAC", "HLA-A", "HLA-B", "HLA-C", "B2M", "HLAABC"],
    },
]

# On the slide, but not used as stand-ins for the absent antigens.
NOT_STANDINS = [
    "STAT3",
    "STAT3PY705",
    "STAT5ALPHA",
    "IRF1",
    "IRF3",
    "IRF3PS396",
    "STING",
    "CGAS",
    "CIITA",
    "NFKBP65PS536",
    "PDL1",
    "XRCC1",
    "PARP1",
    "ATM",
    "ATR",
    "RAD50",
    "RAD51",
    "X53BP1",
    "CHK1",
    "CHK2",
    "RPA32",
    "KAP1",
    "NCADHERIN",
    "PCADHERIN",
    "ALPHACATENIN",
    "P63",
]


def load_ids(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


AUDIT_SUBSTRINGS = (
    "CLAUD",
    "CLDN",
    "KU",
    "XRCC",
    "PRKD",
    "DNAP",
    "STAT",
    "HLA",
    "B2M",
    "MHC",
    "CADHER",
    "EMA",
    "EPCAM",
    "TACSTD",
    "KRT",
    "KERAT",
    "MUC",
)


def audit_ids(ids_by_cohort: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for cohort, ids in ids_by_cohort.items():
        for pid in ids:
            hits = [s for s in AUDIT_SUBSTRINGS if s in pid.upper()]
            if hits:
                rows.append({"cohort": cohort, "protein_id": pid, "matched": ",".join(hits)})
    return pd.DataFrame(rows)


def absence_table(ids_by_cohort: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for q in ABSENCE_QUERIES:
        exact = {x.upper() for x in q["exact"]}
        rec = {
            "requested": q["requested"],
            "label": q["label"],
            "alias_list": ",".join(q["exact"]),
        }
        for cohort, ids in ids_by_cohort.items():
            hits = [i for i in ids if i.upper() in exact]
            rec[f"{cohort}_present"] = bool(hits)
            rec[f"{cohort}_hits"] = ",".join(hits)
        rows.append(rec)
    # positive controls: the antibodies we do test
    for pid, (gene, name, role) in ANTIBODY.items():
        rec = {
            "requested": gene,
            "label": f"{name} ({gene})",
            "alias_list": pid,
        }
        for cohort, ids in ids_by_cohort.items():
            rec[f"{cohort}_present"] = pid in ids
            rec[f"{cohort}_hits"] = pid if pid in ids else ""
        rows.append(rec)
    return pd.DataFrame(rows)


def standin_inventory(ids_by_cohort: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for pid in NOT_STANDINS:
        rec = {"protein_id": pid, "used_as_standin": False}
        for cohort, ids in ids_by_cohort.items():
            rec[f"{cohort}_on_panel"] = pid in set(ids)
        rows.append(rec)
    return pd.DataFrame(rows)


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    n = int(len(d))
    rec = {
        "n_pairwise": n,
        "rho": np.nan,
        "p": np.nan,
        "tested": False,
        "note": "not tested",
    }
    if n == 0:
        rec["note"] = "n=0; Spearman not computed"
        return rec
    if n < 8:
        rec["note"] = f"n={n} < 8; Spearman not computed"
        return rec
    if d["x"].nunique() < 2 or d["y"].nunique() < 2:
        rec["note"] = f"n={n} but no variance"
        return rec
    rho, p = stats.spearmanr(d["x"], d["y"])
    rec.update(
        {
            "rho": float(rho),
            "p": float(p),
            "tested": True,
            "note": "two-sided Spearman, pairwise complete",
        }
    )
    return rec


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return []
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = q
    return out.tolist()


def epithelial_score(wide: pd.DataFrame) -> pd.Series:
    sub = wide[EPITHELIAL]
    z = (sub - sub.mean()) / sub.std(ddof=0)
    # mean of available z; require all three so the score is the same trio
    complete = z.dropna()
    score = pd.Series(np.nan, index=wide.index, name="EPITHELIAL_ZMEAN")
    score.loc[complete.index] = complete.mean(axis=1)
    return score


def concordance(wide_by: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for cohort, wide in wide_by.items():
        for i, a in enumerate(EPITHELIAL):
            for b in EPITHELIAL[i + 1 :]:
                rec = spearman(wide[a], wide[b])
                rec.update(
                    {
                        "cohort": cohort,
                        "protein_a": a,
                        "protein_b": b,
                        "gene_a": ANTIBODY[a][0],
                        "gene_b": ANTIBODY[b][0],
                    }
                )
                rows.append(rec)
    return pd.DataFrame(rows)


def plot(wide_by: dict[str, pd.DataFrame], pairs: pd.DataFrame, path: Path) -> None:
    cohorts = ["LUAD", "LUSC"]
    fig, axes = plt.subplots(4, 3, figsize=(9.2, 11.2))
    row_spec = [
        ("LUAD", "KU80"),
        ("LUSC", "KU80"),
        ("LUAD", "HLADQA1"),
        ("LUSC", "HLADQA1"),
    ]
    for i, (cohort, target) in enumerate(row_spec):
        wide = wide_by[cohort]
        gene_t, name_t, _ = ANTIBODY[target]
        for j, epi in enumerate(EPITHELIAL):
            ax = axes[i, j]
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            gene_e, name_e, _ = ANTIBODY[epi]
            sub = wide[[epi, target]].dropna()
            row = pairs[
                (pairs.cohort == cohort)
                & (pairs.predictor == epi)
                & (pairs.endpoint == target)
                & (pairs.family == "primary")
            ]
            n_target = int(wide[target].notna().sum())
            n_tumors = int(len(wide))
            if sub.empty:
                ax.text(
                    0.5,
                    0.5,
                    f"{name_t} unquantified\n{n_target} / {n_tumors} tumors",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=8,
                )
                ax.set_xticks([])
                ax.set_yticks([])
            else:
                ax.scatter(sub[epi], sub[target], s=12, c="#243447", alpha=0.7, edgecolors="none")
            if len(row) and bool(row.iloc[0]["tested"]):
                rho = float(row.iloc[0]["rho"])
                p = float(row.iloc[0]["p"])
                n = int(row.iloc[0]["n_pairwise"])
                q = row.iloc[0]["q_bh"]
                qtxt = f"  q={q:.3g}" if pd.notna(q) else ""
                ax.set_title(f"{cohort}  {name_e} vs {name_t}\nρ={rho:.3f}  p={p:.3g}{qtxt}  n={n}", fontsize=8)
            else:
                ax.set_title(f"{cohort}  {name_e} vs {name_t}\nnot tested", fontsize=8)
            ax.set_xlabel(f"{name_e} ({gene_e})", fontsize=7)
            ax.set_ylabel(f"{name_t} ({gene_t})", fontsize=7)
            ax.tick_params(labelsize=7)
    fig.suptitle(
        "TCPA RPPA500 level 4  ·  epithelial antibodies vs Ku80 and HLA-DQA1\n"
        "CLDN4, DNA-PKcs, Ku70, STAT1, HLA class I absent; HLA-DQA1 all NA in lung",
        fontsize=11,
    )
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def match_stored(stored: pd.DataFrame, cohort: str, a: str, b: str) -> tuple[float, float]:
    if stored.empty:
        return np.nan, np.nan
    hit = stored[
        (stored.cohort == cohort)
        & (
            ((stored.protein_id == a) & (stored.other_protein_id == b))
            | ((stored.protein_id == b) & (stored.other_protein_id == a))
        )
    ]
    if hit.empty:
        return np.nan, np.nan
    # if both directions exist they should agree; take the first
    return float(hit.iloc[0]["tcpa_stored_rho"]), float(hit.iloc[0]["tcpa_stored_p"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/tcpa_rppa_epithelial_ku_hla")
    ap.add_argument("--outdir", default="methods/tcpa_rppa_epithelial_ku_hla")
    args = ap.parse_args()
    data = Path(args.data)
    res = Path(args.outdir) / "results"
    res.mkdir(parents=True, exist_ok=True)

    ids_by = {c: load_ids(data / f"protein_ids_{c}.txt") for c in ("LUAD", "LUSC")}
    for cohort, ids in ids_by.items():
        (res / f"protein_ids_{cohort}.txt").write_text("\n".join(ids) + "\n")
    presence = absence_table(ids_by)
    presence.to_csv(res / "presence.tsv", sep="\t", index=False)
    standin_inventory(ids_by).to_csv(res / "not_standins.tsv", sep="\t", index=False)
    audit_ids(ids_by).to_csv(res / "panel_name_audit.tsv", sep="\t", index=False)

    long = pd.read_csv(data / "abundances_long.tsv", sep="\t")
    stored = pd.read_csv(data / "tcpa_stored_spearman.tsv", sep="\t")

    n_rows = []
    pair_rows = []
    wide_by = {}
    score_frames = []

    for cohort in ("LUAD", "LUSC"):
        sub = long[long.cohort == cohort]
        wide = sub.pivot_table(index="sample_id", columns="protein_id", values="abundance", aggfunc="mean")
        for pid in EPITHELIAL + TARGETS:
            if pid not in wide.columns:
                wide[pid] = np.nan
        wide = wide[EPITHELIAL + TARGETS]
        n_unique = int(wide.index.nunique())
        n_dup = int(sub.duplicated(["sample_id", "protein_id"]).sum())
        wide_by[cohort] = wide
        score = epithelial_score(wide)
        out = wide.copy()
        out.insert(0, "cohort", cohort)
        out.insert(1, "sample_id", wide.index)
        out["EPITHELIAL_ZMEAN"] = score.to_numpy()
        score_frames.append(out.reset_index(drop=True))

        for pid in EPITHELIAL + TARGETS:
            s = wide[pid]
            n_rows.append(
                {
                    "cohort": cohort,
                    "protein_id": pid,
                    "gene": ANTIBODY[pid][0],
                    "name": ANTIBODY[pid][1],
                    "n_tumors_in_matrix": n_unique,
                    "n_observed": int(s.notna().sum()),
                    "n_na": int(s.isna().sum()),
                    "n_duplicate_sample_rows": n_dup,
                }
            )

        for epi in EPITHELIAL:
            for target in TARGETS:
                rec = spearman(wide[epi], wide[target])
                rho_s, p_s = match_stored(stored, cohort, epi, target)
                if int(wide[target].notna().sum()) == 0:
                    rec["note"] = (
                        "endpoint protein_id is on the TCPA catalog but every "
                        "tumor in this histology is NA; Spearman not computed"
                    )
                rec.update(
                    {
                        "cohort": cohort,
                        "family": "primary",
                        "predictor": epi,
                        "endpoint": target,
                        "predictor_gene": ANTIBODY[epi][0],
                        "endpoint_gene": ANTIBODY[target][0],
                        "tcpa_stored_rho": rho_s,
                        "tcpa_stored_p": p_s,
                    }
                )
                if rec["tested"] and np.isfinite(rho_s):
                    rec["rho_minus_stored"] = float(rec["rho"] - rho_s)
                else:
                    rec["rho_minus_stored"] = np.nan
                pair_rows.append(rec)

        for target in TARGETS:
            rec = spearman(score, wide[target])
            rec.update(
                {
                    "cohort": cohort,
                    "family": "sensitivity_zmean",
                    "predictor": "EPITHELIAL_ZMEAN",
                    "endpoint": target,
                    "predictor_gene": "CLDN7+CDH1+MUC1 z-mean",
                    "endpoint_gene": ANTIBODY[target][0],
                    "tcpa_stored_rho": np.nan,
                    "tcpa_stored_p": np.nan,
                    "rho_minus_stored": np.nan,
                }
            )
            pair_rows.append(rec)

        # explicit untested rows for absent requested antigens
        for q in ABSENCE_QUERIES:
            pair_rows.append(
                {
                    "n_pairwise": 0,
                    "rho": np.nan,
                    "p": np.nan,
                    "tested": False,
                    "note": "antibody absent from TCPA disease L4 panel; not replaced",
                    "cohort": cohort,
                    "family": "absent",
                    "predictor": "epithelial_or_CLDN4",
                    "endpoint": q["requested"],
                    "predictor_gene": "",
                    "endpoint_gene": q["label"],
                    "tcpa_stored_rho": np.nan,
                    "tcpa_stored_p": np.nan,
                    "rho_minus_stored": np.nan,
                }
            )

    pairs = pd.DataFrame(pair_rows)
    pairs["q_bh"] = np.nan
    pairs["bh_family_n"] = np.nan
    for cohort in ("LUAD", "LUSC"):
        mask = (pairs.cohort == cohort) & (pairs.family == "primary") & (pairs.tested)
        idx = pairs.index[mask]
        if len(idx):
            pairs.loc[idx, "q_bh"] = bh(pairs.loc[idx, "p"].tolist())
            pairs.loc[idx, "bh_family_n"] = int(len(idx))

    conc = concordance(wide_by)
    conc.to_csv(res / "epithelial_concordance.tsv", sep="\t", index=False)
    pan_path = data / "pancan_lung_missingness.tsv"
    if pan_path.exists():
        pd.read_csv(pan_path, sep="\t").to_csv(res / "pancan_lung_missingness.tsv", sep="\t", index=False)

    pd.DataFrame(n_rows).to_csv(res / "n_table.tsv", sep="\t", index=False)
    pairs.to_csv(res / "spearman.tsv", sep="\t", index=False)
    pd.concat(score_frames, ignore_index=True).to_csv(res / "sample_scores.tsv", sep="\t", index=False)
    plot(wide_by, pairs, res / "fig_epithelial_vs_ku80_hladqa1.png")

    summary = {
        "source": "TCPA RPPA500 disease-specific level 4 (RBN)",
        "url": "https://tcpa.drbioright.org/rppa500/",
        "cohorts_separate": True,
        "primary_family": "CLDN7, CDH1, MUC1/EMA vs XRCC5/Ku80 and HLA-DQA1",
        "bh": "Benjamini-Hochberg within histology across primary tests that had data (3 Ku80 tests; HLA-DQA1 rows were all NA and were not entered)",
        "cldn4": "absent",
        "dnapkcs": "absent",
        "ku70": "absent",
        "stat1": "absent",
        "hla_class_I": "absent",
        "pairs": pairs[pairs.family != "absent"][
            [
                "cohort",
                "family",
                "predictor",
                "endpoint",
                "n_pairwise",
                "rho",
                "p",
                "q_bh",
                "tested",
                "tcpa_stored_rho",
                "rho_minus_stored",
            ]
        ].to_dict(orient="records"),
    }
    def _jsonable(obj):
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_jsonable(v) for v in obj]
        if isinstance(obj, float) and not np.isfinite(obj):
            return None
        return obj

    (res / "summary.json").write_text(json.dumps(_jsonable(summary), indent=2) + "\n")

    show = pairs[pairs.family != "absent"][
        ["cohort", "family", "predictor", "endpoint", "n_pairwise", "rho", "p", "q_bh", "rho_minus_stored"]
    ]
    with pd.option_context("display.max_rows", 100, "display.width", 160, "display.float_format", lambda v: f"{v:.6g}"):
        print(show.to_string(index=False))
    print("max |rho - stored|", np.nanmax(np.abs(pairs["rho_minus_stored"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
