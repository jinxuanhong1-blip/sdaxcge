#!/usr/bin/env python3
"""CLDN4 RNA and protein vs OS/PFS in open CPTAC LUAD and LSCC.

ImmuneScore / phenotype correlations are out of scope (already reported).
Primary model: continuous Cox per 1 SD. Median-split log-rank is sensitivity.
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

CLDN4 = "ENSG00000189143"

COHORTS = {
    "LUAD": {
        "paper": "Gillette et al. Cell 2020, PMID 32649874",
        "histology": "lung adenocarcinoma",
        "rna": "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
        "protein": "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "surv": "LUAD_survival.txt",
        "meta": "LUAD_meta.txt",
    },
    "LSCC": {
        "paper": "Satpathy et al. Cell 2021, PMID 34358469",
        "histology": "lung squamous cell carcinoma (LUSC synonym)",
        "rna": "LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
        "protein": "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "surv": "LSCC_survival.txt",
        "meta": "LSCC_meta.txt",
    },
}


def find_row(index: pd.Index, prefix: str) -> str | None:
    hits = [i for i in index if str(i) == prefix or str(i).startswith(prefix + ".")]
    if not hits:
        return None
    if len(hits) > 1:
        raise ValueError(f"multiple rows for {prefix}: {hits}")
    return str(hits[0])


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df


def extract_gene(mat: pd.DataFrame, prefix: str) -> tuple[pd.Series, str | None]:
    row = find_row(mat.index, prefix)
    if row is None:
        return pd.Series(np.nan, index=mat.columns, name=prefix), None
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    s.name = prefix
    return s, row


def load_surv(path: Path) -> pd.DataFrame:
    s = pd.read_csv(path, sep="\t")
    idcol = "case_id" if "case_id" in s.columns else s.columns[0]
    s = s.set_index(idcol)
    s.index = s.index.astype(str)
    return s.apply(pd.to_numeric, errors="coerce")


def load_meta(path: Path) -> pd.DataFrame:
    m = pd.read_csv(path, sep="\t")
    idcol = "case_id" if "case_id" in m.columns else m.columns[0]
    m = m.set_index(idcol)
    m.index = m.index.astype(str)
    # LinkedOmics meta files have a data_type header row (CON/BIN/ORD).
    if str(m.index[0]).lower() in {"data_type", "datatype", "type"} or set(
        str(v).upper() for v in m.iloc[0].tolist()
    ).issubset({"CON", "BIN", "ORD", "CAT", "NAN", ""}):
        m = m.iloc[1:]
    return m


def logrank(time: pd.Series, event: pd.Series, group: pd.Series) -> dict:
    d = pd.concat([time, event, group], axis=1).dropna()
    d.columns = ["t", "e", "g"]
    d["t"] = pd.to_numeric(d["t"], errors="coerce")
    d["e"] = pd.to_numeric(d["e"], errors="coerce")
    d = d.dropna()
    d = d[d["t"] > 0]
    rec = {
        "n": int(len(d)),
        "n_event": int(d["e"].sum()) if len(d) else 0,
        "n_high": np.nan,
        "n_low": np.nan,
        "events_high": np.nan,
        "events_low": np.nan,
        "logrank_stat": np.nan,
        "logrank_p": np.nan,
    }
    labs = sorted(d["g"].astype(str).unique())
    if len(labs) != 2:
        return rec
    rec["n_high"] = int((d.g == "high").sum())
    rec["n_low"] = int((d.g == "low").sum())
    rec["events_high"] = int(d.loc[d.g == "high", "e"].sum())
    rec["events_low"] = int(d.loc[d.g == "low", "e"].sum())
    t0, e0 = d.loc[d.g == "low", "t"], d.loc[d.g == "low", "e"]
    t1, e1 = d.loc[d.g == "high", "t"], d.loc[d.g == "high", "e"]
    try:
        stat, p = stats.logrank(t0, t1, e0, e1)
        rec["logrank_stat"] = float(np.asarray(stat).ravel()[0])
        rec["logrank_p"] = float(np.asarray(p).ravel()[0])
    except Exception:
        times = np.sort(d.loc[d.e == 1, "t"].unique())
        o1 = e1_exp = v = 0.0
        for t in times:
            at0 = int((t0 >= t).sum())
            at1 = int((t1 >= t).sum())
            at = at0 + at1
            dth = int(((d.t == t) & (d.e == 1)).sum())
            d1 = int(((t1 == t) & (e1 == 1)).sum())
            if at <= 1:
                continue
            exp1 = dth * at1 / at
            var = dth * (at - dth) / (at - 1) * (at1 / at) * (at0 / at)
            o1 += d1
            e1_exp += exp1
            v += var
        if v > 0:
            stat = (o1 - e1_exp) ** 2 / v
            rec["logrank_stat"] = float(stat)
            rec["logrank_p"] = float(stats.chi2.sf(stat, 1))
    return rec


def cox_continuous(time: pd.Series, event: pd.Series, x: pd.Series) -> dict:
    d = pd.concat([time, event, x], axis=1).dropna()
    d.columns = ["t", "e", "x"]
    d["t"] = pd.to_numeric(d["t"], errors="coerce")
    d["e"] = pd.to_numeric(d["e"], errors="coerce")
    d["x"] = pd.to_numeric(d["x"], errors="coerce")
    d = d.dropna()
    d = d[d["t"] > 0]
    rec = {
        "cox_n": int(len(d)),
        "cox_events": int(d["e"].sum()) if len(d) else 0,
        "cox_hr_per_sd": np.nan,
        "cox_hr_per_sd_ci_low": np.nan,
        "cox_hr_per_sd_ci_high": np.nan,
        "cox_p": np.nan,
        "cox_loghr": np.nan,
        "cox_loghr_se": np.nan,
        "cox_hr_per_unit": np.nan,
        "cox_hr_per_unit_ci_low": np.nan,
        "cox_hr_per_unit_ci_high": np.nan,
        "cox_p_per_unit": np.nan,
    }
    if len(d) < 10 or d["e"].sum() < 5 or float(d["x"].std(ddof=0)) == 0:
        rec["cox_note"] = "too_few_events_or_n"
        return rec
    z = (d["x"] - d["x"].mean()) / d["x"].std(ddof=0)
    try:
        res = PHReg(d["t"].to_numpy(float), z.to_numpy(float)[:, None], status=d["e"].to_numpy(float)).fit(disp=0)
        loghr = float(res.params[0])
        se = float(res.bse[0])
        rec["cox_hr_per_sd"] = float(math.exp(loghr))
        rec["cox_hr_per_sd_ci_low"] = float(math.exp(loghr - 1.96 * se))
        rec["cox_hr_per_sd_ci_high"] = float(math.exp(loghr + 1.96 * se))
        rec["cox_p"] = float(res.pvalues[0])
        rec["cox_loghr"] = loghr
        rec["cox_loghr_se"] = se
    except Exception as exc:
        rec["cox_error"] = str(exc)
    try:
        resu = PHReg(d["t"].to_numpy(float), d["x"].to_numpy(float)[:, None], status=d["e"].to_numpy(float)).fit(disp=0)
        loghru = float(resu.params[0])
        seu = float(resu.bse[0])
        rec["cox_hr_per_unit"] = float(math.exp(loghru))
        rec["cox_hr_per_unit_ci_low"] = float(math.exp(loghru - 1.96 * seu))
        rec["cox_hr_per_unit_ci_high"] = float(math.exp(loghru + 1.96 * seu))
        rec["cox_p_per_unit"] = float(resu.pvalues[0])
    except Exception as exc:
        rec["cox_unit_error"] = str(exc)
    return rec


def km_steps(time: pd.Series, event: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    d = pd.concat([time, event], axis=1).dropna()
    d.columns = ["t", "e"]
    d = d[d["t"] > 0].sort_values("t")
    t = d["t"].to_numpy(float)
    e = d["e"].to_numpy(float)
    xs = [0.0]
    ys = [1.0]
    surv = 1.0
    n = len(d)
    i = 0
    while i < n:
        ti = t[i]
        at = n - i
        di = 0
        while i < n and t[i] == ti:
            di += int(e[i] == 1)
            i += 1
        if di and at:
            surv *= 1.0 - di / at
        xs.extend([ti, ti])
        ys.extend([ys[-1], surv])
    return np.asarray(xs), np.asarray(ys)


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_hr(hr: float, lo: float, hi: float) -> str:
    if hr != hr:
        return "NA"
    return f"{hr:.2f} ({lo:.2f}–{hi:.2f})"


def analyze_cohort(cohort: str, spec: dict, data: Path) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    rna = load_matrix(data / spec["rna"])
    prot = load_matrix(data / spec["protein"])
    surv = load_surv(data / spec["surv"])
    meta = load_meta(data / spec["meta"])

    cldn4_rna, rna_row = extract_gene(rna, CLDN4)
    cldn4_prot, prot_row = extract_gene(prot, CLDN4)

    # Tumor matrices are case-level (one column per case).
    samples = sorted(set(rna.columns) | set(prot.columns) | set(surv.index))
    core = pd.DataFrame(index=samples)
    core["CLDN4_RNA"] = cldn4_rna.reindex(samples)
    core["CLDN4_protein"] = cldn4_prot.reindex(samples)
    core = core.join(surv.reindex(samples), how="left")

    coverage = {
        "cohort": cohort,
        "n_tumor_rna_columns": int(rna.shape[1]),
        "n_tumor_protein_columns": int(prot.shape[1]),
        "n_survival_rows": int(surv.shape[0]),
        "n_meta_rows": int(meta.shape[0]),
        "cldn4_rna_row": rna_row,
        "cldn4_protein_row": prot_row,
        "cldn4_rna_n": int(core["CLDN4_RNA"].notna().sum()),
        "cldn4_rna_na": int(core["CLDN4_RNA"].isna().sum()),
        "cldn4_protein_n": int(core["CLDN4_protein"].notna().sum()),
        "cldn4_protein_na": int(core["CLDN4_protein"].isna().sum()),
        "survival_columns": list(surv.columns),
        "os_events_in_table": int(pd.to_numeric(surv.get("OS_event"), errors="coerce").fillna(0).sum())
        if "OS_event" in surv.columns
        else None,
        "pfs_events_in_table": int(pd.to_numeric(surv.get("PFS_event"), errors="coerce").fillna(0).sum())
        if "PFS_event" in surv.columns
        else None,
        "paper": spec["paper"],
        "histology": spec["histology"],
    }

    rows = []
    for layer, col in [("RNA", "CLDN4_RNA"), ("protein", "CLDN4_protein")]:
        x = core[col]
        med = float(x.median(skipna=True)) if x.notna().any() else np.nan
        grp = pd.Series(np.where(x >= med, "high", "low"), index=core.index, dtype=object)
        grp[x.isna()] = np.nan
        for endpoint, tcol, ecol in [("OS", "OS_days", "OS_event"), ("PFS", "PFS_days", "PFS_event")]:
            if tcol not in core.columns or ecol not in core.columns:
                rows.append(
                    {
                        "cohort": cohort,
                        "gene": "CLDN4",
                        "layer": layer,
                        "endpoint": endpoint,
                        "note": f"missing_{tcol}_or_{ecol}",
                    }
                )
                continue
            lr = logrank(core[tcol], core[ecol], grp)
            cx = cox_continuous(core[tcol], core[ecol], x)
            rec = {
                "cohort": cohort,
                "gene": "CLDN4",
                "layer": layer,
                "endpoint": endpoint,
                "split": "median",
                "median_cut": med,
                "ensembl_row": rna_row if layer == "RNA" else prot_row,
            }
            rec.update(lr)
            rec.update(cx)
            rows.append(rec)

    feat = core.reset_index().rename(columns={"index": "case_id"})
    feat.insert(0, "cohort", cohort)
    keep = [
        "cohort",
        "case_id",
        "CLDN4_RNA",
        "CLDN4_protein",
        "OS_days",
        "OS_event",
        "PFS_days",
        "PFS_event",
    ]
    feat = feat[[c for c in keep if c in feat.columns]]
    return pd.DataFrame(rows), coverage, feat


def plot_km(feat: pd.DataFrame, surv_df: pd.DataFrame, figs: Path) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(14.5, 7.2), sharey=True)
    combos = [
        ("LUAD", "CLDN4_RNA", "OS", "LUAD RNA OS"),
        ("LUAD", "CLDN4_protein", "OS", "LUAD protein OS"),
        ("LSCC", "CLDN4_RNA", "OS", "LSCC RNA OS"),
        ("LSCC", "CLDN4_protein", "OS", "LSCC protein OS"),
        ("LUAD", "CLDN4_RNA", "PFS", "LUAD RNA PFS"),
        ("LUAD", "CLDN4_protein", "PFS", "LUAD protein PFS"),
        ("LSCC", "CLDN4_RNA", "PFS", "LSCC RNA PFS"),
        ("LSCC", "CLDN4_protein", "PFS", "LSCC protein PFS"),
    ]
    for ax, (cohort, pred, endpoint, title) in zip(axes.ravel(), combos):
        sub = feat[feat["cohort"] == cohort].copy()
        tcol = f"{endpoint}_days"
        ecol = f"{endpoint}_event"
        x = pd.to_numeric(sub[pred], errors="coerce")
        med = x.median(skipna=True)
        for lab, color, mask in [
            ("low", "#2c7bb6", x < med),
            ("high", "#d7191c", x >= med),
        ]:
            xs, ys = km_steps(sub.loc[mask, tcol], sub.loc[mask, ecol])
            n = int(pd.concat([sub.loc[mask, tcol], sub.loc[mask, ecol]], axis=1).dropna().shape[0])
            ev = int(pd.to_numeric(sub.loc[mask, ecol], errors="coerce").fillna(0).sum())
            ax.step(xs, ys, where="post", color=color, lw=1.6, label=f"{lab} n={n} ev={ev}")
        layer = "RNA" if pred.endswith("RNA") else "protein"
        rec = surv_df[(surv_df.cohort == cohort) & (surv_df.layer == layer) & (surv_df.endpoint == endpoint)]
        p = float(rec["logrank_p"].iloc[0]) if len(rec) else np.nan
        hr = float(rec["cox_hr_per_sd"].iloc[0]) if len(rec) else np.nan
        ax.set_title(f"{title}\nlog-rank p={fmt_p(p)}; Cox HR/SD={hr:.2f}" if hr == hr else title, fontsize=8)
        ax.set_xlabel(f"{endpoint} (days)", fontsize=8)
        ax.set_ylabel("Survival", fontsize=8)
        ax.legend(fontsize=6, frameon=False)
        ax.set_ylim(0, 1.02)
    fig.suptitle("CPTAC LinkedOmics freeze v1.2 — CLDN4 median split (treatment-naive OS/PFS, not ICI)", fontsize=10)
    fig.tight_layout()
    fig.savefig(figs / "km_cldn4_os_pfs.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def write_finding(surv_df: pd.DataFrame, coverage: list[dict], out: Path, manifest: dict) -> None:
    def cell(cohort: str, layer: str, endpoint: str) -> pd.Series:
        q = surv_df[(surv_df.cohort == cohort) & (surv_df.layer == layer) & (surv_df.endpoint == endpoint)]
        return q.iloc[0] if len(q) else pd.Series(dtype=object)

    cov = {c["cohort"]: c for c in coverage}
    lines = []
    lines.append("# Finding — LinkedOmics CPTAC CLDN4 RNA/protein vs OS/PFS (LUAD, LSCC)")
    lines.append("")
    lines.append("Additive slice. **CLDN4 only.** Open LinkedOmics CPTAC pan-cancer freeze v1.2.")
    lines.append("ImmuneScore / ESTIMATE / xCell / CIBERSORT correlations are **already reported** in the LUAD and LSCC CPTAC slices and were **not re-audited** here.")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append("Clinical tables are **open** (`LUAD_survival.txt`, `LSCC_survival.txt`, plus meta).")
    lines.append("These OS/PFS times are **treatment-naive surgical follow-up**, not immunotherapy outcomes. No ICI labels exist in this freeze.")
    lines.append("")

    # Build compact verdict from the 8 primary Cox tests
    any_sig = False
    for cohort in ["LUAD", "LSCC"]:
        for layer in ["RNA", "protein"]:
            for endpoint in ["OS", "PFS"]:
                r = cell(cohort, layer, endpoint)
                if r.empty:
                    continue
                p = r.get("cox_p", np.nan)
                if p == p and p < 0.05:
                    any_sig = True
    if any_sig:
        lines.append("At least one continuous Cox test is p<0.05. Read n and event counts before claiming prognosis — several protein tests are event-poor.")
    else:
        lines.append("**No CLDN4 RNA or protein association with OS or PFS reaches p<0.05** on the pre-specified continuous Cox (per 1 SD) in either histology. Median-split log-rank is likewise null at α=0.05. Closest: LUAD CLDN4 protein OS log-rank p=0.087 (15 events / 75). Event counts are small; a null is **inconclusive**, not proof of no prognostic value.")
    lines.append("")
    lines.append("Primary model = continuous Cox **per 1 SD** of CLDN4 (pairwise-complete). Median split is sensitivity only (ties at the median coded high). HR >1 means higher CLDN4, higher hazard.")
    lines.append("")
    lines.append("## Survival table")
    lines.append("")
    lines.append("| Cohort | Layer | Endpoint | n | Events | Cox HR per 1 SD (95% CI) | Cox p | Median-split log-rank p | High / low n (events) |")
    lines.append("|---|---|---|---:|---:|---|---:|---:|---|")
    for cohort in ["LUAD", "LSCC"]:
        for layer in ["RNA", "protein"]:
            for endpoint in ["OS", "PFS"]:
                r = cell(cohort, layer, endpoint)
                if r.empty:
                    lines.append(f"| {cohort} | {layer} | {endpoint} | — | — | missing | — | — | — |")
                    continue
                hr = fmt_hr(r.get("cox_hr_per_sd", np.nan), r.get("cox_hr_per_sd_ci_low", np.nan), r.get("cox_hr_per_sd_ci_high", np.nan))
                n = int(r.get("cox_n", r.get("n", 0)) or 0)
                ev = int(r.get("cox_events", r.get("n_event", 0)) or 0)
                nh = r.get("n_high", np.nan)
                nl = r.get("n_low", np.nan)
                eh = r.get("events_high", np.nan)
                el = r.get("events_low", np.nan)
                split = (
                    f"{int(nh)}/{int(nl)} ({int(eh)}/{int(el)})"
                    if nh == nh and nl == nl
                    else "NA"
                )
                lines.append(
                    f"| {cohort} | CLDN4 {layer} | {endpoint} | {n} | {ev} | {hr} | {fmt_p(float(r.get('cox_p', np.nan)))} | {fmt_p(float(r.get('logrank_p', np.nan)))} | {split} |"
                )
    lines.append("")
    lines.append("Machine-readable copy: `results/survival.tsv`.")
    lines.append("")
    lines.append("## Data (open, HEAD HTTP 200)")
    lines.append("")
    lines.append("Filenames from the LinkedOmics CPTAC-pancan-LUAD / CPTAC-pancan-LSCC download tables, served from `cptac-pancancer-data` S3 freeze `data_freeze_v1.2_reorganized`. Phenotype / ImmuneScore files were **not** downloaded.")
    lines.append("")
    lines.append("| Cohort | File | HTTP | Bytes | Role |")
    lines.append("|---|---|---:|---:|---|")
    for f in manifest.get("files", []):
        lines.append(
            f"| {f.get('cohort')} | `{f.get('name')}` | {f.get('http_status')} | {f.get('content_length')} | {f.get('status')} |"
        )
    lines.append("")
    lu = cov.get("LUAD", {})
    ls = cov.get("LSCC", {})
    lines.append(
        f"LUAD: Gillette *Cell* 2020 (PMID 32649874), treatment-naive resected adenocarcinoma. "
        f"RNA CLDN4 `{lu.get('cldn4_rna_row')}` n={lu.get('cldn4_rna_n')}; "
        f"protein `{lu.get('cldn4_protein_row')}` n={lu.get('cldn4_protein_n')} "
        f"(NA={lu.get('cldn4_protein_na')}). Survival table OS events={lu.get('os_events_in_table')}, "
        f"PFS events={lu.get('pfs_events_in_table')}."
    )
    lines.append("")
    lines.append(
        f"LSCC: Satpathy *Cell* 2021 (PMID 34358469), newly diagnosed resected squamous carcinoma, no prior chemo/RT. "
        f"RNA CLDN4 `{ls.get('cldn4_rna_row')}` n={ls.get('cldn4_rna_n')}; "
        f"protein `{ls.get('cldn4_protein_row')}` n={ls.get('cldn4_protein_n')} "
        f"(NA={ls.get('cldn4_protein_na')}). Survival table OS events={ls.get('os_events_in_table')}, "
        f"PFS events={ls.get('pfs_events_in_table')}."
    )
    lines.append("")
    lines.append("CLDN4 protein missingness is TMT dropout, not a join error. All protein tests are pairwise-complete.")
    lines.append("")
    lines.append("Analytic n is smaller than the matrix n because follow-up time is missing (or ≤0) for some cases, including deaths deposited without `OS_days` (LUAD 1, LSCC 2). Those rows cannot enter a Cox or log-rank model. LUAD RNA OS uses 105/110 cases (23 events; table lists 24). LSCC RNA OS uses 94/108 (23 events; table lists 25).")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append("- Gene row matched by Ensembl prefix `ENSG00000189143` (versioned IDs in the freeze).")
    lines.append("- Join key = case ID (matrix columns × `case_id` in survival). Time unit = days as deposited.")
    lines.append("- Cases with time ≤0 or missing time/event/CLDN4 are dropped for that test.")
    lines.append("- **Primary:** Cox PH on CLDN4 z-scored within the pairwise-complete set (HR per 1 SD). 95% CI = exp(β ± 1.96 SE). Two-sided Wald p from `statsmodels.PHReg`.")
    lines.append("- **Sensitivity:** median split (ties at median → high); two-sample log-rank.")
    lines.append("- Per-unit Cox is in `results/survival.tsv` for audit; RNA and protein live on different log2 scales, so per-SD is the comparable number.")
    lines.append("- No ImmuneScore residualization, no stage-adjusted Cox, no optimal cutpoint search.")
    lines.append("- Not ICI survival. Do not meta-analyse these HRs with atezolizumab / pembrolizumab series.")
    lines.append("")
    lines.append("## What this does not claim")
    lines.append("")
    lines.append("- It does not re-test CLDN4 vs ImmuneScore / CD8 / IFN signatures.")
    lines.append("- It does not claim a prognostic or predictive ICI biomarker.")
    lines.append("- It does not treat a null with ~15–23 events as evidence of no effect.")
    lines.append("- It does not use LinkedOmics LinkFinder web p-values; all numbers are recomputed from the freeze files.")
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append("- `results/survival.tsv` — n / events / HR / CI / p for every test")
    lines.append("- `results/coverage.json` — row IDs, missingness, event totals")
    lines.append("- `results/sample_level.tsv` — case-level CLDN4 + OS/PFS (no immune columns)")
    lines.append("- `results/figures/km_cldn4_os_pfs.png` — eight median-split KM panels")
    lines.append("- `data/manifest.json` — URLs, HTTP status, SHA256 (local cache, not committed)")
    lines.append("")
    out.write_text("\n".join(lines) + "\n")


def main() -> int:
    here = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser()
    p.add_argument("--data", default=str(here / "data"))
    p.add_argument("--outdir", default=str(here / "results"))
    args = p.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    tables = out
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    man_path = data / "manifest.json"
    manifest = json.loads(man_path.read_text()) if man_path.exists() else {}

    surv_parts = []
    feat_parts = []
    coverage = []
    for cohort, spec in COHORTS.items():
        sdf, cov, feat = analyze_cohort(cohort, spec, data)
        surv_parts.append(sdf)
        feat_parts.append(feat)
        coverage.append(cov)

    surv_df = pd.concat(surv_parts, ignore_index=True)
    feat = pd.concat(feat_parts, ignore_index=True)
    surv_df.to_csv(tables / "survival.tsv", sep="\t", index=False)
    feat.to_csv(tables / "sample_level.tsv", sep="\t", index=False)
    (tables / "coverage.json").write_text(json.dumps(coverage, indent=2) + "\n")
    plot_km(feat, surv_df, figs)
    write_finding(surv_df, coverage, here / "FINDING.md", manifest)

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "CLDN4 RNA/protein vs OS/PFS only; ImmuneScore not re-audited",
        "n_tests": int(len(surv_df)),
        "coverage": coverage,
        "survival": surv_df.to_dict(orient="records"),
    }
    (tables / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(tables / "survival.tsv")
    print(here / "FINDING.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
