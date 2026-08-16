# %% [markdown]
# # 07 - Joint readout: SCENIC regulons × LIANA recruitment edges
#
# **Question.** Do the transcription-factor regulons that define the
# TACSTD2/CLDN4-high epithelial state *line up with* weaker epithelium→T
# recruitment edges (CXCL9/10/11–CXCR3, CCL5–CCR5, CXCL16–CXCR6)?
#
# This script does **not** re-run SCENIC or LIANA. It joins their already-
# computed tables and writes a one-page "can we claim reduced recruitment?"
# verdict that is deliberately conservative:
#
# | Evidence | Source | Can support a claim? |
# |---|---|---|
# | Patient-level chemokine pseudobulk high < low | `05` | **yes, primary** |
# | LIANA magnitude_rank worse from high sender | `01` | hypothesis only |
# | High-state regulons anticorrelated with chemokines | `04b` | hypothesis only |
# | Agreement of the three | this script | strengthens the *hypothesis*, still not causation |
#
# Overclaim to avoid: "TF X represses CXCL10, therefore T cells are not
# recruited." That needs motif + perturbation + spatial evidence the
# playbook lists in §8.

# %%
from pathlib import Path
import yaml
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CFG = yaml.safe_load((HERE.parent / "config" / "gene_sets.yaml").read_text())
P = CFG["params"]
chemokines = sorted({g for ax in CFG["t_cell_recruitment"].values()
                     for g in ax["ligands"]})
focus_receptors = {"CXCR3", "CCR5", "CXCR6", "CCR1"}

RES = Path("results")
OUT = RES / "joint"; OUT.mkdir(parents=True, exist_ok=True)


def _read(path):
    p = Path(path)
    return pd.read_csv(p) if p.exists() else None


pb = _read(RES / "pseudobulk" / "chemokine_high_vs_low_pseudobulk.csv")
liana = _read(RES / "liana" / "liana_high_vs_low_patientlevel.csv")
# fall back to the descriptive (not patient-level) table
if liana is None:
    liana = _read(RES / "liana" / "liana_high_vs_low.csv")
reg = _read(RES / "scenic" / "regulons_high_vs_low_patientlevel.csv")
corr = _read(RES / "scenic" / "regulon_chemokine_correlation.csv")

# %% [markdown]
# ## 1. Score each evidence stream independently

# %%
verdict = {"accession_note": "join-only; no new inference"}

if pb is not None and len(pb):
    # Negative delta = lower chemokine in the HIGH state (supports hypothesis).
    sig = pb.dropna(subset=["padj"]) if "padj" in pb.columns else pb
    support = sig[sig["mean_delta_high_minus_low"] < 0]
    verdict["pseudobulk"] = {
        "n_chemokines_tested": int(len(sig)),
        "n_lower_in_high": int(len(support)),
        "n_fdr05_lower_in_high": int(((support.get("padj", 1) < 0.05)).sum())
        if "padj" in support.columns else None,
        "supports_reduced_signal": bool(len(support) >= 2),
        "table": "results/pseudobulk/chemokine_high_vs_low_pseudobulk.csv",
    }
else:
    verdict["pseudobulk"] = {"supports_reduced_signal": None,
                             "reason": "run 05_pseudobulk_guardrail.py first"}

if liana is not None and len(liana):
    df = liana.copy()
    # Accept either the patient-level schema or the descriptive pivot.
    if {"median_rank_hi", "median_rank_lo"}.issubset(df.columns):
        # higher rank number = weaker (LIANA magnitude_rank)
        weaker = df["median_rank_hi"] > df["median_rank_lo"]
        rec_ok = df.get("receptor", df.get("receptor_complex", "")).astype(str)
        rec_ok = rec_ok.apply(lambda x: any(r in x for r in focus_receptors))
        focus = df[rec_ok]
        weaker_f = (focus["median_rank_hi"] > focus["median_rank_lo"]) if len(focus) else weaker
        verdict["liana"] = {
            "n_pairs": int(len(df)),
            "n_focus_pairs": int(len(focus)),
            "frac_weaker_from_high_focus": float(weaker_f.mean()) if len(focus) else None,
            "supports_reduced_signal": bool(len(focus) and weaker_f.mean() >= 0.5),
            "note": "hypothesis only; do not cite as communication proof",
        }
    elif {"high", "low"}.issubset(df.columns):
        weaker = df["high"].fillna(np.inf) > df["low"].fillna(np.inf)
        verdict["liana"] = {
            "n_pairs": int(len(df)),
            "frac_weaker_from_high": float(weaker.mean()),
            "supports_reduced_signal": bool(weaker.mean() >= 0.5),
            "note": "descriptive (not patient-level) — hypothesis only",
        }
    else:
        verdict["liana"] = {"supports_reduced_signal": None,
                            "reason": "unrecognized LIANA table schema"}
else:
    verdict["liana"] = {"supports_reduced_signal": None,
                        "reason": "run 01_liana_lr.py first"}

if corr is not None and len(corr) and reg is not None and len(reg):
    # High-state-up regulons (positive delta) anticorrelated with chemokines.
    up = set(reg.sort_values("delta", ascending=False).head(10)["regulon"])
    sub = corr[corr["regulon"].isin(up) & corr["chemokine"].isin(chemokines)]
    anti = sub[sub["pearson_r"] < -0.15]
    verdict["scenic"] = {
        "n_up_regulon_chemokine_pairs": int(len(sub)),
        "n_anticorrelated_r_lt_-0.15": int(len(anti)),
        "supports_repression_hypothesis": bool(len(anti) >= 2),
        "note": "association only; repression needs motif + perturbation",
    }
    anti.to_csv(OUT / "highstate_regulon_chemokine_anticorrelation.csv", index=False)
else:
    verdict["scenic"] = {"supports_repression_hypothesis": None,
                         "reason": "run 04_pyscenic.sh + 04b first"}

# %% [markdown]
# ## 2. Conservative joint call
# Require the PRIMARY (pseudobulk) to support, plus at least one of LIANA /
# SCENIC as a consistent hypothesis. Never promote the joint call to
# "mechanistic proof".

# %%
primary = verdict["pseudobulk"].get("supports_reduced_signal")
hyp = [verdict["liana"].get("supports_reduced_signal"),
       verdict["scenic"].get("supports_repression_hypothesis")]
if primary is True and any(h is True for h in hyp):
    call = ("CONSISTENT HYPOTHESIS: patient-level chemokine reduction in "
            "TACSTD2/CLDN4-high epithelium, with matching LIANA/SCENIC "
            "direction. Still not proof of recruitment or of TF repression. "
            "Triangulate with spatial / bulk / IHC (playbook §3, §8).")
elif primary is True:
    call = ("PRIMARY ONLY: patient-level chemokine reduction, but LIANA/SCENIC "
            "did not agree or were not run. Report the pseudobulk; do not "
            "over-interpret mechanism.")
elif primary is False:
    call = ("NOT SUPPORTED at the patient-level chemokine test. Do not claim "
            "reduced T-cell recruitment from LIANA/SCENIC/NicheNet alone.")
else:
    call = "INCOMPLETE: run 05 (required) and 01/04b (optional) before judging."

verdict["joint_call"] = call
pd.Series(verdict).to_json(OUT / "joint_verdict.json", indent=2)
print("=== joint verdict ===")
print(call)
print(f"[done] {OUT / 'joint_verdict.json'}")
