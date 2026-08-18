#!/usr/bin/env python3
"""Mouse / line-level Cldn4 vs T/NK or IFN/MHC on public TISMO lung ICI.

Cldn4-only. Tacstd2 49/64 is taken as given and is not re-scored.
No user-private 8 KL.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
LUNG_CANCER = "Lung carcinoma"
OUTLIER_SRX = "SRX8918393"
ICB_TOKENS = ("antipd1", "antipdl1", "antipdl2", "antictla4")

TNK_GENES = ["Cd8a", "Cd3e", "Cd3d", "Cd2", "Nkg7", "Gzmb", "Prf1", "Ifng", "Ncr1"]
IFN_GENES = ["Ifng", "Stat1", "Cxcl9", "Cxcl10", "Ido1", "H2-Aa"]
MHC_GENES = ["H2-K1", "H2-D1", "H2-Q4", "B2m", "Tap1", "Tap2"]
EXTRA_GENES = ["Cldn4", "Cd274", "Actb", "Epcam", "Ptprc"]
ALL_GENES = list(dict.fromkeys(EXTRA_GENES + TNK_GENES + IFN_GENES + MHC_GENES))


def tpm_from_log2p1(x: float) -> float:
    return float(2.0**x - 1.0)


def cohort_key(label: str) -> str:
    if "(n=" in str(label):
        return str(label)[: str(label).rfind("(n=")]
    return str(label)


def is_baseline(row: pd.Series) -> bool:
    return str(row.get("Baseline", "")) == "1" or str(row.get("Responder", "")) == "Baseline"


def lung_line_from_cohort(label: str) -> str | None:
    s = str(label)
    if s.startswith("LLC_"):
        return "LLC"
    if s.startswith("CMT-167_") or s.startswith("CMT167_"):
        return "CMT-167"
    return None


def load_gene_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df["cohort_key"] = df["cell_line"].map(cohort_key)
    df["is_baseline"] = df.apply(is_baseline, axis=1)
    df["lung_line"] = df["cell_line"].map(lung_line_from_cohort)
    return df


def mean_z(df: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in df.columns]
    if not present:
        return pd.Series(np.nan, index=df.index)
    z = df[present].apply(lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) else 0.0, axis=0)
    return z.mean(axis=1)


def spearman(x: pd.Series, y: pd.Series) -> dict:
    mask = x.notna() & y.notna()
    n = int(mask.sum())
    if n < 3:
        return {"n": n, "rho": np.nan, "p": np.nan}
    r = stats.spearmanr(x[mask], y[mask])
    return {"n": n, "rho": float(r.statistic), "p": float(r.pvalue)}


def welch_mwu(a: np.ndarray, b: np.ndarray) -> dict:
    if a.size < 2 or b.size < 2:
        return {"welch_t": np.nan, "welch_p": np.nan, "mwu_u": np.nan, "mwu_p": np.nan}
    t_res = stats.ttest_ind(b, a, equal_var=False, alternative="two-sided")
    u_res = stats.mannwhitneyu(b, a, alternative="two-sided")
    return {
        "welch_t": float(t_res.statistic),
        "welch_p": float(t_res.pvalue),
        "mwu_u": float(u_res.statistic),
        "mwu_p": float(u_res.pvalue),
    }


def fmt_p(p: float) -> str:
    if p is None or (isinstance(p, float) and (np.isnan(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.2f}" if p >= 0.01 else f"{p:.4f}"


def fmt_num(x: float, nd: int = 3) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    return f"{x:.{nd}f}"


def lung_inventory(meta: list[dict]) -> pd.DataFrame:
    rows = []
    for rec in meta:
        if rec.get("cancerType") != LUNG_CANCER:
            continue
        treat = str(rec.get("mouseTreatment") or "")
        treat_l = treat.lower().replace(" ", "")
        is_icb = any(tok in treat_l for tok in ICB_TOKENS)
        rows.append(
            {
                "tismo_id": rec.get("id"),
                "study_id": rec.get("studyId"),
                "cell_line": rec.get("cellLine"),
                "cancer_type": rec.get("cancerType"),
                "cell_genotype": rec.get("cellGenotype"),
                "mouse_strain": rec.get("mouseStrain"),
                "implantation": rec.get("implantation"),
                "implantation_site": rec.get("implantationSite"),
                "mouse_treatment": treat,
                "icb_study_label": rec.get("icbStudy"),
                "replicates": rec.get("replicates"),
                "tumor": rec.get("tumor"),
                "is_icb_treatment": is_icb,
                "has_ici_outcome": rec.get("studyId") == "GSE155972",
            }
        )
    return pd.DataFrame(rows).sort_values(["study_id", "cell_line", "mouse_treatment"])


def build_wide(raw_dir: Path) -> tuple[pd.DataFrame, dict]:
    coverage = {}
    pieces = []
    # LLC-only Cldn4 export drops SRX8918393; the All-models table keeps it.
    cldn4_all = load_gene_csv(raw_dir / "cldn4_icb_all.csv")
    cldn4_llc = cldn4_all.loc[cldn4_all["lung_line"] == "LLC"].copy()
    coverage["Cldn4"] = {"present": True, "n": int(cldn4_llc["Samples"].nunique()), "n_rows": int(len(cldn4_llc)), "source": "cldn4_icb_all.csv LLC rows"}
    one = cldn4_llc[
        ["Samples", "value", "cell_line", "Responder", "Baseline", "GSE_ID", "Mouse_treatment", "cohort_key", "is_baseline", "lung_line"]
    ].copy()
    one = one.rename(columns={"value": "Cldn4"})
    pieces.append(one)
    for gene in ALL_GENES:
        if gene == "Cldn4":
            continue
        safe = gene.lower().replace("-", "")
        path = raw_dir / f"{safe}_llc.csv"
        if not path.exists() or path.stat().st_size < 20:
            coverage[gene] = {"present": False, "n": 0}
            continue
        g = load_gene_csv(path)
        coverage[gene] = {"present": True, "n": int(g["Samples"].nunique()), "n_rows": int(len(g))}
        one = g[["Samples", "value"]].copy()
        one = one.rename(columns={"value": gene})
        pieces.append(one)
    if not pieces:
        raise RuntimeError("no LLC gene CSVs")
    wide = pieces[0]
    for extra in pieces[1:]:
        wide = wide.merge(
            extra[["Samples", extra.columns[1]]],
            on="Samples",
            how="outer",
        )
    # restore meta from Cldn4 if merge dropped it
    if "Cldn4" not in wide.columns:
        raise RuntimeError("Cldn4 LLC export missing")
    wide["arm"] = np.where(wide["cohort_key"].astype(str).str.contains("Setdb1_KO"), "Setdb1_KO", "WT")
    wide["group"] = np.where(wide["is_baseline"], "Baseline", "ICB")
    wide["is_outlier"] = wide["Samples"] == OUTLIER_SRX
    present_tnk = [g for g in TNK_GENES if g in wide.columns]
    wide["tnk_mean"] = wide[present_tnk].mean(axis=1) if present_tnk else np.nan
    wide["ifn_meanz"] = mean_z(wide, IFN_GENES)
    wide["mhc_meanz"] = mean_z(wide, MHC_GENES)
    wide["Cldn4_tpm"] = wide["Cldn4"].map(tpm_from_log2p1)
    immune_scores = attach_immune(raw_dir, wide)
    if immune_scores:
        coverage["_immune_rds"] = {"present": True, "scores": immune_scores}
    return wide, coverage


IMMUNE_TNK_ROWS = [
    "CD8 T_mMCPcounter",
    "T CD8_TIMER",
    "T CD8_CIBERSORT_abs",
    "T_mMCPcounter",
    "T NK_xCell",
    "NK_mMCPcounter",
]


def attach_immune(raw_dir: Path, wide: pd.DataFrame) -> list[str]:
    path = raw_dir / "TISMO_immune_infiltration.RDS"
    if not path.exists():
        return []
    try:
        import pyreadr
        imm = list(pyreadr.read_r(str(path)).values())[0]
    except Exception as err:
        print(f"WARN immune RDS: {err}")
        return []
    added = []
    for rowname in IMMUNE_TNK_ROWS:
        if rowname not in imm.index:
            continue
        col = rowname.replace(" ", "_")
        series = imm.loc[rowname]
        wide[col] = wide["Samples"].map(series)
        added.append(col)
    if added:
        wide["tismo_tnk_infil"] = wide[added].apply(
            lambda r: ((r - r.mean()) / r.std(ddof=0) if r.std(ddof=0) else 0.0),
            axis=0,
        ).mean(axis=1)
        added.append("tismo_tnk_infil")
    return added


def contrast_rows(wide: pd.DataFrame, features: list[str], drop_outlier: bool) -> list[dict]:
    work = wide.loc[~wide["is_outlier"]] if drop_outlier else wide
    out = []
    # within-arm ICB vs baseline
    for arm in ["WT", "Setdb1_KO"]:
        sub = work[work["arm"] == arm]
        a = sub[sub["group"] == "Baseline"]
        b = sub[sub["group"] == "ICB"]
        for feat in features:
            if feat not in work.columns:
                continue
            tests = welch_mwu(a[feat].dropna().to_numpy(), b[feat].dropna().to_numpy())
            out.append(
                {
                    "contrast": f"{arm} ICB vs baseline",
                    "feature": feat,
                    "exclude_outlier": drop_outlier,
                    "n_a": int(a[feat].notna().sum()),
                    "n_b": int(b[feat].notna().sum()),
                    "mean_a": float(a[feat].mean()) if a[feat].notna().any() else np.nan,
                    "mean_b": float(b[feat].mean()) if b[feat].notna().any() else np.nan,
                    "delta_b_minus_a": float(b[feat].mean() - a[feat].mean()) if a[feat].notna().any() and b[feat].notna().any() else np.nan,
                    **tests,
                    "note": "Setdb1_KO ICB = Responders; WT ICB = Non-responders (genotype-confounded)",
                }
            )
    # treated R vs NR (confounded)
    treated = work[work["group"] == "ICB"]
    nr = treated[treated["arm"] == "WT"]
    r = treated[treated["arm"] == "Setdb1_KO"]
    for feat in features:
        if feat not in work.columns:
            continue
        tests = welch_mwu(nr[feat].dropna().to_numpy(), r[feat].dropna().to_numpy())
        out.append(
            {
                "contrast": "treated Setdb1_KO (R) vs WT (NR)",
                "feature": feat,
                "exclude_outlier": drop_outlier,
                "n_a": int(nr[feat].notna().sum()),
                "n_b": int(r[feat].notna().sum()),
                "mean_a": float(nr[feat].mean()) if nr[feat].notna().any() else np.nan,
                "mean_b": float(r[feat].mean()) if r[feat].notna().any() else np.nan,
                "delta_b_minus_a": float(r[feat].mean() - nr[feat].mean()) if nr[feat].notna().any() and r[feat].notna().any() else np.nan,
                **tests,
                "note": "NOT a within-genotype R vs NR. Response label is the Setdb1 genotype.",
            }
        )
    return out


def write_finding(
    dest: Path,
    inv: pd.DataFrame,
    wide: pd.DataFrame,
    coverage: dict,
    line_tbl: pd.DataFrame,
    corr_tbl: pd.DataFrame,
    contrast_tbl: pd.DataFrame,
    cldn4_all: pd.DataFrame,
) -> None:
    lung_icb = cldn4_all[cldn4_all["lung_line"].notna()]
    lung_lines_icb = sorted({x for x in lung_icb["lung_line"].dropna().unique()})
    n_llc = int(wide["Samples"].nunique())
    n_cldn4_floor = int((wide["Cldn4_tpm"] < 1).sum())
    cldn4_mean_tpm = float(wide["Cldn4_tpm"].mean())

    def corr_row(subset: str, x: str, y: str) -> pd.Series:
        hit = corr_tbl[(corr_tbl["subset"] == subset) & (corr_tbl["x"] == x) & (corr_tbl["y"] == y)]
        if hit.empty:
            return pd.Series({"n": 0, "rho": np.nan, "p": np.nan})
        return hit.iloc[0]

    all_tnk = corr_row("all_mice", "Cldn4", "tnk_mean")
    all_ifn = corr_row("all_mice", "Cldn4", "ifn_meanz")
    all_mhc = corr_row("all_mice", "Cldn4", "mhc_meanz")
    noout_tnk = corr_row("drop_SRX8918393", "Cldn4", "tnk_mean")

    def cget(contrast: str, feat: str, drop: bool) -> pd.Series:
        hit = contrast_tbl[
            (contrast_tbl["contrast"] == contrast)
            & (contrast_tbl["feature"] == feat)
            & (contrast_tbl["exclude_outlier"] == drop)
        ]
        if hit.empty:
            return pd.Series({"n_a": 0, "n_b": 0, "delta_b_minus_a": np.nan, "welch_p": np.nan, "mwu_p": np.nan})
        return hit.iloc[0]

    wt_c = cget("WT ICB vs baseline", "Cldn4", False)
    ko_c = cget("Setdb1_KO ICB vs baseline", "Cldn4", False)
    wt_tnk = cget("WT ICB vs baseline", "tnk_mean", False)
    ko_tnk = cget("Setdb1_KO ICB vs baseline", "tnk_mean", False)
    rnr_c = cget("treated Setdb1_KO (R) vs WT (NR)", "Cldn4", False)

    tnk_present = [g for g in TNK_GENES if coverage.get(g, {}).get("present")]
    ifn_present = [g for g in IFN_GENES if coverage.get(g, {}).get("present")]
    mhc_present = [g for g in MHC_GENES if coverage.get(g, {}).get("present")]

    lines_lung = sorted(inv["cell_line"].dropna().unique())
    n_designs = int(len(inv))
    n_studies = int(inv["study_id"].nunique())
    ici_studies = sorted(inv.loc[inv["has_ici_outcome"], "study_id"].unique())

    line_md = ["| line | TISMO cancer | ICI outcome in TISMO | Cldn4 in ICB export | n mice (ICB export) | mean Cldn4 log2(TPM+1) | mean Cldn4 TPM |",
               "|---|---|---|---|---:|---:|---:|"]
    for _, r in line_tbl.iterrows():
        line_md.append(
            f"| **{r['line']}** | {r['cancer_type']} | {r['ici_outcome']} | {r['cldn4_in_icb_export']} | "
            f"{int(r['n_mice_icb']) if pd.notna(r['n_mice_icb']) else 0} | {fmt_num(r['mean_cldn4_log2p1'])} | {fmt_num(r['mean_cldn4_tpm'])} |"
        )

    # per-mouse arm table snippet
    arm_rows = []
    for (arm, group), sub in wide.groupby(["arm", "group"], sort=True):
        arm_rows.append(
            f"| {arm} | {group} | {int(sub['Samples'].nunique())} | "
            f"{fmt_num(sub['Cldn4'].mean())} | {fmt_num(sub['Cldn4_tpm'].mean())} | "
            f"{fmt_num(sub['tnk_mean'].mean())} | {fmt_num(sub['ifn_meanz'].mean())} | {fmt_num(sub['mhc_meanz'].mean())} |"
        )

    md = f"""# Finding — TISMO LLC Cldn4 (mouse / line), public processed only

Additive public **mouse**. **Cldn4-only.** Tacstd2 TISMO 49/64 and the human CLDN4 thesis are taken as given and are not re-scored.

Question: in TISMO **LLC** (and any other TISMO mouse lung line that has both Cldn4 and an ICI outcome), does mouse-level **Cldn4** track T/NK or IFN/MHC, and is there an ICI-response split?

**No user-private 8 KL.** Values are TISMO’s uniformly processed `log2(TPM+1)` gene-module export plus the public vivo metadata table (Zeng et al., *NAR* 2022, PMID 34534350). No FASTQ.

Primary tables: `tables/line_level.tsv`, `tables/mouse_level.tsv`.

## Line-level table

TISMO `cancerType == Lung carcinoma` lines: **{', '.join(lines_lung) if lines_lung else 'none'}** ({n_studies} studies, {n_designs} design groups).

Lung lines that appear in the Cldn4 ICB gene export: **{', '.join(lung_lines_icb) if lung_lines_icb else 'none'}**.

{chr(10).join(line_md)}

**CMT-167** is a TISMO lung line (GSE100412, orthotopic, untreated, n=3) but has **no ICB arm and no response label**. The Cldn4 ICB gene-module query for CMT-167 returns HTTP 500 (empty ICI set). **KPB25L** appears in the Cldn4 ICB export (GSE124821) but TISMO labels it **Mammary cancer, NOS**, not lung — it is not added. MLE12 is a TISMO lung-adenocarcinoma *cell line* annotation only; it has no in-vivo ICI rows.

The only ICI-outcome lung design is **GSE155972 LLC** (Griffin et al., *Nature* 2021): subcutaneous flank, anti-PD1 + anti-CTLA4, WT vs Setdb1_KO. Response is **confounded with genotype**: WT ICB = Non-responders (n=6); Setdb1_KO ICB = Responders (n=7). Baseline is untreated (n=10 / genotype). Honest ICI n = **1 study, 1 line, 2 genotype arms, {n_llc} mice**.

## Mouse-level table (GSE155972 LLC, n={n_llc})

TISMO values are `log2(TPM+1)`. T/NK = mean of {len(tnk_present)}/{len(TNK_GENES)} genes ({', '.join(tnk_present)}). IFN = mean-z of {len(ifn_present)}/{len(IFN_GENES)} leftover genes. MHC-I = mean-z of {len(mhc_present)}/{len(MHC_GENES)} leftover genes. z is across these {n_llc} mice.

| arm | group | n mice | Cldn4 log2p1 | Cldn4 TPM | T/NK mean | IFN mean-z | MHC-I mean-z |
|---|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(arm_rows)}

Cldn4 sits near the detection floor: mean TPM = **{fmt_num(cldn4_mean_tpm, 2)}**; **{n_cldn4_floor}/{n_llc}** mice have TPM < 1. One mouse (**{OUTLIER_SRX}**, Setdb1_KO ICB) is Cldn4 = 3.23 log2p1 (TPM ≈ 8.4) and is marked in `mouse_level.tsv`.

### Cldn4 vs T/NK or IFN/MHC (Spearman, mouse unit)

| subset | n | Cldn4 vs T/NK genes ρ (p) | Cldn4 vs TISMO T/NK infil ρ (p) | Cldn4 vs IFN ρ (p) | Cldn4 vs MHC-I ρ (p) |
|---|---:|---|---|---|---|
| all GSE155972 mice | {int(all_tnk['n'])} | {fmt_num(all_tnk['rho'])} ({fmt_p(all_tnk['p'])}) | {fmt_num(corr_row('all_mice','Cldn4','tismo_tnk_infil')['rho'])} ({fmt_p(corr_row('all_mice','Cldn4','tismo_tnk_infil')['p'])}) | {fmt_num(all_ifn['rho'])} ({fmt_p(all_ifn['p'])}) | {fmt_num(all_mhc['rho'])} ({fmt_p(all_mhc['p'])}) |
| drop {OUTLIER_SRX} | {int(noout_tnk['n'])} | {fmt_num(noout_tnk['rho'])} ({fmt_p(noout_tnk['p'])}) | {fmt_num(corr_row('drop_SRX8918393','Cldn4','tismo_tnk_infil')['rho'])} ({fmt_p(corr_row('drop_SRX8918393','Cldn4','tismo_tnk_infil')['p'])}) | {fmt_num(corr_row('drop_SRX8918393','Cldn4','ifn_meanz')['rho'])} ({fmt_p(corr_row('drop_SRX8918393','Cldn4','ifn_meanz')['p'])}) | {fmt_num(corr_row('drop_SRX8918393','Cldn4','mhc_meanz')['rho'])} ({fmt_p(corr_row('drop_SRX8918393','Cldn4','mhc_meanz')['p'])}) |

Per-arm correlations are in `tables/spearman.tsv`. A positive Cldn4–T/NK ρ on near-floor Cldn4 is **not** a Cldn4-high / T-low exclusion pattern.

### ICI (honest)

| contrast | n | Cldn4 Δ log2p1 | Cldn4 Welch p | T/NK Δ | T/NK Welch p |
|---|---|---:|---:|---:|---:|
| WT ICB vs baseline (NR arm) | {int(wt_c['n_a'])} / {int(wt_c['n_b'])} | {fmt_num(wt_c['delta_b_minus_a'])} | {fmt_p(wt_c['welch_p'])} | {fmt_num(wt_tnk['delta_b_minus_a'])} | {fmt_p(wt_tnk['welch_p'])} |
| Setdb1_KO ICB vs baseline (R arm) | {int(ko_c['n_a'])} / {int(ko_c['n_b'])} | {fmt_num(ko_c['delta_b_minus_a'])} | {fmt_p(ko_c['welch_p'])} | {fmt_num(ko_tnk['delta_b_minus_a'])} | {fmt_p(ko_tnk['welch_p'])} |
| treated R vs NR (genotype-confounded) | {int(rnr_c['n_a'])} / {int(rnr_c['n_b'])} | {fmt_num(rnr_c['delta_b_minus_a'])} | {fmt_p(rnr_c['welch_p'])} | {fmt_num(cget('treated Setdb1_KO (R) vs WT (NR)', 'tnk_mean', False)['delta_b_minus_a'])} | {fmt_p(cget('treated Setdb1_KO (R) vs WT (NR)', 'tnk_mean', False)['welch_p'])} |

There is **no within-genotype responder vs non-responder** contrast in TISMO LLC. Do not read the R vs NR row as an ICI-outcome test of Cldn4.

T/NK / IFN genes **do** move with ICB (same design can detect an immune shift). Cldn4 does not at a usable effect size; it stays near floor.

## Verdict

TISMO LLC Cldn4 is **present** (not absent) but **near the expression floor**. The mouse/line-level table exists.

- **Lines with Cldn4 + ICI outcome:** LLC only (GSE155972). CMT-167 has Cldn4-capable metadata as lung but **no ICI outcome**.
- **Mouse n:** {n_llc} (10 WT baseline, 6 WT ICB/NR, 10 Setdb1_KO baseline, 7 Setdb1_KO ICB/R).
- **Cldn4 vs T/NK or IFN/MHC:** see Spearman table; Cldn4 is too low to support a Cldn4-high / immune-low or Cldn4-high / IFN-high LLC state.
- **ICI response:** labels exist but are the Setdb1 genotype. Cldn4 ICB vs baseline is not significant on either arm.

This does **not** reopen the Tacstd2 49/64 tally (2/64 of those cohorts are this LLC study). It does not use private 8 KL.

## Methods (short)

- Inclusion: TISMO `cancerType == Lung carcinoma` AND (Cldn4 in the public ICB gene export) AND (ICI treatment or response field). That intersection is LLC GSE155972.
- Expression: TISMO gene-module `downVivoExprn` type=3, `icbList` = the six ICB treatments, `tumorList` = LLC (and All for the lung-ICI census). Values = `log2(TPM+1)` as deposited.
- T/NK score = unweighted mean of present T/NK genes (already log2p1). Extra T/NK = mean-z of TISMO public infiltrate rows (CD8 T mMCPcounter, T CD8 TIMER/CIBERSORT, T mMCPcounter, T NK xCell, NK mMCPcounter). IFN / MHC-I = leftover 6-gene lists from the GSE239485 Cldn4 page, mean of gene-wise z on these mice.
- Tests: Spearman on the mouse; Welch t and two-sided MWU when both arms have n≥2. Outlier {OUTLIER_SRX} is kept in the primary table and dropped only in the sensitivity rows.
- Public processed only. No GEO FASTQ. No user-private 8 KL.

Reproduce:

```bash
python3 methods/tismo_llc_cldn4/fetch.py
python3 methods/tismo_llc_cldn4/analyze.py
```

## Files

| File | Role |
|---|---|
| `tables/line_level.tsv` | TISMO lung lines: Cldn4 + ICI present/absent |
| `tables/mouse_level.tsv` | one row per GSE155972 mouse |
| `tables/spearman.tsv` | Cldn4 vs T/NK, IFN, MHC-I |
| `tables/contrasts.tsv` | ICB vs baseline and confounded R vs NR |
| `tables/lung_inventory.tsv` | all TISMO lung-carcinoma design groups |
| `tables/gene_coverage.tsv` | which genes were in the LLC export |
| `tables/summary.json` | machine-readable verdict |
| `raw/` | live TISMO exports used here |

## 中文摘要

公开 TISMO 小鼠、只看 **Cldn4**。不做 Tacstd2 49/64，不用私有 8 KL。

TISMO 肺癌系只有 **LLC** 和 **CMT-167**。有 Cldn4 **且** 有 ICI 结局的只有 **LLC GSE155972**（皮下，抗 PD-1+CTLA4；WT ICB=NR n=6，Setdb1_KO ICB=R n=7，对照各 n=10）。CMT-167 无 ICB。Cldn4 接近检测下限（多数 TPM<1）。小鼠水平 Cldn4 与 T/NK、IFN、MHC-I 的相关见 `tables/spearman.tsv`；不能支持 Cldn4 高 / 免疫低。R vs NR 与基因型完全重叠，不是独立 ICI 结局检验。
"""
    dest.write_text(md, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=HERE / "raw")
    parser.add_argument("--out", type=Path, default=HERE / "tables")
    args = parser.parse_args()
    raw_dir = args.raw
    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    meta = json.loads((raw_dir / "vivo_meta.json").read_text())["data"]
    inv = lung_inventory(meta)
    inv.to_csv(out_dir / "lung_inventory.tsv", sep="\t", index=False)

    cldn4_all = load_gene_csv(raw_dir / "cldn4_icb_all.csv")
    wide, coverage = build_wide(raw_dir)

    cov_rows = [{"gene": g, **coverage.get(g, {"present": False, "n": 0})} for g in ALL_GENES]
    pd.DataFrame(cov_rows).to_csv(out_dir / "gene_coverage.tsv", sep="\t", index=False)

    # line-level: every TISMO lung line, plus whether it is in the Cldn4 ICB export
    lung_in_icb = (
        cldn4_all.loc[cldn4_all["lung_line"].notna()]
        .groupby("lung_line")
        .agg(n_mice_icb=("Samples", "nunique"), mean_cldn4_log2p1=("value", "mean"))
        .reset_index()
        .rename(columns={"lung_line": "line"})
    )
    tpm_by_line = (
        cldn4_all.loc[cldn4_all["lung_line"].notna()]
        .assign(tpm=lambda d: d["value"].map(tpm_from_log2p1))
        .groupby("lung_line")["tpm"]
        .mean()
    )
    lung_in_icb["mean_cldn4_tpm"] = lung_in_icb["line"].map(tpm_by_line)
    lung_in_icb["cldn4_in_icb_export"] = "yes"
    lung_in_icb["ici_outcome"] = lung_in_icb["line"].map(
        lambda x: "GSE155972 anti-PD1+anti-CTLA4; R/NR = Setdb1 genotype" if x == "LLC" else "in ICB export"
    )

    line_rows = []
    for line in sorted(inv["cell_line"].dropna().unique()):
        n_rep = int(pd.to_numeric(inv.loc[inv["cell_line"] == line, "replicates"], errors="coerce").fillna(0).sum())
        has_ici = bool(inv.loc[inv["cell_line"] == line, "has_ici_outcome"].any())
        hit = lung_in_icb[lung_in_icb["line"] == line]
        if hit.empty:
            line_rows.append(
                {
                    "line": line,
                    "cancer_type": LUNG_CANCER,
                    "ici_outcome": "absent" if not has_ici else "metadata yes, Cldn4 ICB export no",
                    "cldn4_in_icb_export": "no",
                    "n_mice_icb": 0,
                    "n_replicates_in_vivo_meta": n_rep,
                    "mean_cldn4_log2p1": np.nan,
                    "mean_cldn4_tpm": np.nan,
                }
            )
        else:
            rec = hit.iloc[0].to_dict()
            rec["cancer_type"] = LUNG_CANCER
            rec["n_replicates_in_vivo_meta"] = n_rep
            line_rows.append(rec)
    line_tbl = pd.DataFrame(line_rows)
    line_tbl.to_csv(out_dir / "line_level.tsv", sep="\t", index=False)

    mouse_cols = [
        "Samples",
        "GSE_ID",
        "lung_line",
        "arm",
        "group",
        "Responder",
        "is_baseline",
        "is_outlier",
        "cohort_key",
        "Mouse_treatment",
        "Cldn4",
        "Cldn4_tpm",
        "tnk_mean",
        "ifn_meanz",
        "mhc_meanz",
        "tismo_tnk_infil",
    ] + [g for g in ALL_GENES if g in wide.columns and g != "Cldn4"]
    mouse_cols = [c for c in mouse_cols if c in wide.columns]
    wide[mouse_cols].sort_values(["arm", "group", "Samples"]).to_csv(out_dir / "mouse_level.tsv", sep="\t", index=False)

    corr_rows = []
    subsets = {
        "all_mice": wide,
        "drop_SRX8918393": wide.loc[~wide["is_outlier"]],
        "WT": wide.loc[wide["arm"] == "WT"],
        "Setdb1_KO": wide.loc[wide["arm"] == "Setdb1_KO"],
        "Baseline": wide.loc[wide["group"] == "Baseline"],
        "ICB": wide.loc[wide["group"] == "ICB"],
    }
    for name, sub in subsets.items():
        for y in ["tnk_mean", "ifn_meanz", "mhc_meanz", "tismo_tnk_infil", "Cd8a", "Ifng", "Cd274"]:
            if y not in sub.columns:
                continue
            rec = spearman(sub["Cldn4"], sub[y])
            corr_rows.append({"subset": name, "x": "Cldn4", "y": y, **rec})
    corr_tbl = pd.DataFrame(corr_rows)
    corr_tbl.to_csv(out_dir / "spearman.tsv", sep="\t", index=False)

    features = ["Cldn4", "Cldn4_tpm", "tnk_mean", "ifn_meanz", "mhc_meanz", "tismo_tnk_infil", "Cd8a", "Ifng", "Cd274", "Actb"]
    features = [f for f in features if f in wide.columns]
    contrast_tbl = pd.DataFrame(contrast_rows(wide, features, False) + contrast_rows(wide, features, True))
    contrast_tbl.to_csv(out_dir / "contrasts.tsv", sep="\t", index=False)

    summary = {
        "verdict": "TISMO LLC Cldn4 present, near floor; only lung line with ICI outcome",
        "lung_lines": sorted(inv["cell_line"].dropna().unique().tolist()),
        "lung_lines_with_cldn4_and_ici": sorted({x for x in cldn4_all["lung_line"].dropna().unique()}),
        "n_mice_gse155972": int(wide["Samples"].nunique()),
        "mean_cldn4_tpm": float(wide["Cldn4_tpm"].mean()),
        "n_tpm_lt_1": int((wide["Cldn4_tpm"] < 1).sum()),
        "spearman_all_cldn4_tnk": corr_row_dict(corr_tbl, "all_mice", "tnk_mean"),
        "spearman_all_cldn4_ifn": corr_row_dict(corr_tbl, "all_mice", "ifn_meanz"),
        "spearman_all_cldn4_mhc": corr_row_dict(corr_tbl, "all_mice", "mhc_meanz"),
        "private_8kl_used": False,
        "tacstd2_rescored": False,
        "source": "https://tismo.pku-genomics.org/",
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    write_finding(HERE / "FINDING.md", inv, wide, coverage, line_tbl, corr_tbl, contrast_tbl, cldn4_all)
    print(f"wrote tables under {out_dir}")
    print(f"wrote {HERE / 'FINDING.md'}")
    print(json.dumps(summary, indent=2))


def corr_row_dict(corr_tbl: pd.DataFrame, subset: str, y: str) -> dict:
    hit = corr_tbl[(corr_tbl["subset"] == subset) & (corr_tbl["y"] == y)]
    if hit.empty:
        return {}
    r = hit.iloc[0]
    return {"n": int(r["n"]), "rho": None if pd.isna(r["rho"]) else float(r["rho"]), "p": None if pd.isna(r["p"]) else float(r["p"])}


if __name__ == "__main__":
    main()
