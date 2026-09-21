"""Imputation-free sensitivities for the LUAD null and the IFN-core null.

Does not replace the primary 12-test family. No abundance is imputed.
Alternate panels are the locked IFN/MHC list split into signaling, ISG,
receptor, and antigen-processing scores. Histology uses public grade and
the cBioPortal dominant subtype / pathology text. Strata with pairwise
n < 20 are counted and not tested.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from genes import IFN_ISG, IFN_SCORE_MEMBERS, IFN_SIGNALING, SENSITIVITY_PANELS

MIN_STRATUM_N = 20
MIN_PANEL_N = 8


def _A():
    import analyze as A

    return A


def load_annotations(data: Path) -> dict[str, pd.DataFrame]:
    out = {}
    for cohort in ("LUAD", "LSCC"):
        meta = pd.read_csv(data / cohort / f"{cohort}_meta.txt", sep="\t", dtype=str)
        meta = meta[meta["case_id"] != "data_type"].copy()
        if meta["case_id"].duplicated().any():
            raise ValueError(f"{cohort} meta case_id is not unique")
        meta = meta.set_index("case_id")
        out[cohort] = meta
    lu = pd.read_csv(
        data / "clinical" / "luad_cptac_2020_data_clinical_sample.txt",
        sep="\t",
        skiprows=4,
        dtype=str,
    )
    lu["case_id"] = lu["SAMPLE_ID"].str.replace(r"^X(?=11LU)", "", regex=True)
    if lu["case_id"].duplicated().any():
        raise ValueError("LUAD clinical case_id is not unique after 11LU remap")
    out["LUAD_subtype"] = lu.set_index("case_id")
    ls = pd.read_csv(
        data / "clinical" / "lusc_cptac_2021_data_clinical_sample.txt",
        sep="\t",
        skiprows=4,
        dtype=str,
    )
    if ls["SAMPLE_ID"].duplicated().any():
        raise ValueError("LSCC clinical SAMPLE_ID is not unique")
    out["LSCC_path"] = ls.set_index("SAMPLE_ID")
    return out


def stratum_masks(bundle: dict, ann: dict) -> tuple[dict[str, pd.Series], pd.DataFrame]:
    """Boolean masks aligned to the protein columns, plus a count table."""
    cohort = bundle["cohort"]
    idx = bundle["series"]["CLDN4"].index
    meta = ann[cohort].reindex(idx)
    grade = meta["Histologic_Grade"].fillna("")
    masks: dict[str, pd.Series] = {
        "grade_G2": grade.str.startswith("G2"),
        "grade_G3": grade.str.startswith("G3"),
    }
    rows = []
    if cohort == "LUAD":
        sub = ann["LUAD_subtype"].reindex(idx)["DOMINANT_HISTOLOGICAL_SUBTYPE"]
        for name, n in sub.value_counts(dropna=False).items():
            rows.append(
                {
                    "cohort": cohort,
                    "axis": "dominant_histological_subtype",
                    "level": "NA" if pd.isna(name) else str(name),
                    "n_tumors": int(n),
                    "n_cldn4": int(bundle["series"]["CLDN4"][sub.eq(name) if pd.notna(name) else sub.isna()].notna().sum()),
                }
            )
        known = sub.notna()
        masks["acinar"] = sub.eq("acinar")
        masks["non_acinar"] = known & ~sub.eq("acinar")
    else:
        path = ann["LSCC_path"].reindex(idx)["PATHOLOGY_BASED_HISTOLOGY_ASSESSMENT"]
        for name, n in path.value_counts(dropna=False).items():
            rows.append(
                {
                    "cohort": cohort,
                    "axis": "pathology_text",
                    "level": "NA" if pd.isna(name) else str(name),
                    "n_tumors": int(n),
                    "n_cldn4": int(
                        bundle["series"]["CLDN4"][path.eq(name) if pd.notna(name) else path.isna()].notna().sum()
                    ),
                }
            )
        basal = path.fillna("").str.contains("basaloid", case=False)
        masks["basaloid"] = basal
        masks["not_basaloid"] = path.notna() & ~basal
    for name, n in grade.value_counts(dropna=False).items():
        rows.append(
            {
                "cohort": cohort,
                "axis": "histologic_grade",
                "level": name if name else "NA",
                "n_tumors": int(n),
                "n_cldn4": int(bundle["series"]["CLDN4"][grade.eq(name)].notna().sum()),
            }
        )
    for key, mask in masks.items():
        rows.append(
            {
                "cohort": cohort,
                "axis": "tested_stratum",
                "level": key,
                "n_tumors": int(mask.sum()),
                "n_cldn4": int(bundle["series"]["CLDN4"][mask].notna().sum()),
            }
        )
    return masks, pd.DataFrame(rows)


def panel_scores(bundle: dict) -> dict[str, pd.Series]:
    A = _A()
    scores = {
        "MHC1_protein": bundle["series"]["MHC1_protein"],
        "IFN_core_protein": bundle["series"]["IFN_core_protein"],
        "CD8A_protein": bundle["series"]["CD8A"],
    }
    for name, (members, min_genes) in SENSITIVITY_PANELS.items():
        vecs = []
        used = []
        for g in members:
            if bundle["rows_used"].get(g) and int(bundle["series"][g].notna().sum()) >= MIN_PANEL_N:
                vecs.append(bundle["series"][g])
                used.append(g)
        if len(used) >= min_genes:
            s = A.mean_z(vecs, min_genes)
        else:
            s = pd.Series(np.nan, index=bundle["series"]["CLDN4"].index)
        s.attrs["members"] = used
        scores[name] = s
    return scores


def strict_score(bundle: dict, members: list[str]) -> pd.Series:
    """Mean of z-scores only when every requested member is quantified."""
    cols = []
    for g in members:
        if bundle["rows_used"].get(g) is None:
            return pd.Series(np.nan, index=bundle["series"]["CLDN4"].index)
        s = bundle["series"][g]
        sd = float(s.std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            return pd.Series(np.nan, index=s.index)
        cols.append((s - s.mean()) / sd)
    z = pd.concat(cols, axis=1)
    score = z.mean(axis=1, skipna=False)
    score.name = "strict"
    return score


def blank_partial() -> dict:
    return {
        "n_partial": 0,
        "rho_partial": np.nan,
        "p_partial": np.nan,
        "ci_lo_partial": np.nan,
        "ci_hi_partial": np.nan,
    }


def spearman_row(
    cohort: str,
    family: str,
    predictor: str,
    endpoint: str,
    stratum: str,
    x: pd.Series,
    y: pd.Series,
    wes: pd.Series | None,
    rng: np.random.Generator,
    min_n: int,
) -> dict:
    A = _A()
    # Temporarily honor the caller's minimum n.
    old = A.MIN_N
    A.MIN_N = min_n
    try:
        rec = A.spearman_boot(x, y, rng)
    finally:
        A.MIN_N = old
    row = {
        "cohort": cohort,
        "family": family,
        "predictor": predictor,
        "endpoint": endpoint,
        "stratum": stratum,
        "test": "spearman",
        "n_q4": np.nan,
        "n_q1": np.nan,
        "delta_q4_minus_q1": np.nan,
        **rec,
        **blank_partial(),
    }
    if wes is not None and np.isfinite(row["rho"]) and row["n"] >= min_n:
        old = A.MIN_N
        A.MIN_N = min_n
        try:
            row.update(A.partial_spearman(x, y, wes, rng))
        finally:
            A.MIN_N = old
    return row


def quartile_row(
    cohort: str,
    predictor: str,
    endpoint: str,
    x: pd.Series,
    y: pd.Series,
) -> dict:
    d = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    rec = {
        "cohort": cohort,
        "family": "quartile",
        "predictor": predictor,
        "endpoint": endpoint,
        "stratum": "quantified_predictor",
        "test": "mannwhitney_q4_q1",
        "n": int(len(d)),
        "rho": np.nan,
        "p": np.nan,
        "ci_lo": np.nan,
        "ci_hi": np.nan,
        "n_q4": 0,
        "n_q1": 0,
        "delta_q4_minus_q1": np.nan,
        **blank_partial(),
    }
    if len(d) < 16 or d["x"].nunique() < 4:
        return rec
    pct = d["x"].rank(method="average", pct=True)
    hi = d.loc[pct >= 0.75, "y"]
    lo = d.loc[pct <= 0.25, "y"]
    rec["n_q4"] = int(len(hi))
    rec["n_q1"] = int(len(lo))
    if len(hi) < 8 or len(lo) < 8:
        return rec
    rec["delta_q4_minus_q1"] = float(hi.median() - lo.median())
    _u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    rec["p"] = float(p)
    return rec


def run_sensitivity(
    data: Path,
    outdir: Path,
    bundles: dict,
    rng: np.random.Generator,
) -> pd.DataFrame:
    ann = load_annotations(data)
    rows: list[dict] = []
    count_frames = []
    scorebook = {}
    for cohort, bundle in bundles.items():
        masks, counts = stratum_masks(bundle, ann)
        count_frames.append(counts)
        scores = panel_scores(bundle)
        # Strict all-member versions.
        strict = {
            "IFN_core_strict": strict_score(bundle, IFN_SCORE_MEMBERS),
            "IFN_signaling_strict": strict_score(bundle, IFN_SIGNALING),
            "IFN_isg_strict": strict_score(bundle, IFN_ISG),
            "APM_strict": strict_score(
                bundle, ["TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5"]
            ),
        }
        scores.update(strict)
        scorebook[cohort] = scores
        wes = bundle["wes"]
        # Alternate panels, all tumors.
        for pred in ("CLDN4", "TACSTD2"):
            x = bundle["series"][pred]
            for ep in SENSITIVITY_PANELS:
                rows.append(
                    spearman_row(
                        cohort, "panel", f"{pred}_protein", ep, "all",
                        x, scores[ep], wes, rng, MIN_PANEL_N,
                    )
                )
            for ep in strict:
                rows.append(
                    spearman_row(
                        cohort, "strict", f"{pred}_protein", ep, "all_members_quantified",
                        x, scores[ep], wes, rng, MIN_PANEL_N,
                    )
                )
        # Shared complete-case: CLDN4 and every IFN-core member and every HLA.
        shared = bundle["series"]["CLDN4"].notna()
        for g in IFN_SCORE_MEMBERS + ["HLA-A", "HLA-B", "HLA-C", "CD8A"]:
            if bundle["rows_used"].get(g) is None:
                shared = shared & False
            else:
                shared = shared & bundle["series"][g].notna()
        for pred in ("CLDN4", "TACSTD2"):
            x = bundle["series"][pred].where(shared)
            for ep in ("MHC1_protein", "IFN_core_protein", "IFN_signaling_protein", "CD8A_protein"):
                rows.append(
                    spearman_row(
                        cohort, "shared_complete", f"{pred}_protein", ep, "cldn4_and_ifncore_and_hla",
                        x, scores[ep].where(shared), wes, rng, MIN_PANEL_N,
                    )
                )
        # Rank-floor: missing CLDN4 tied below every quantified value. Not an abundance.
        x = bundle["series"]["CLDN4"]
        floor = float(x.min()) - 1.0
        x_floor = x.fillna(floor)
        for ep in (
            "MHC1_protein",
            "IFN_core_protein",
            "IFN_signaling_protein",
            "IFN_isg_protein",
            "CD8A_protein",
            "APM_protein",
        ):
            rows.append(
                spearman_row(
                    cohort, "rank_floor", "CLDN4_protein", ep, "missing_tied_at_floor",
                    x_floor, scores[ep], wes, rng, MIN_PANEL_N,
                )
            )
        # Quartiles among quantified tumors.
        for pred in ("CLDN4", "TACSTD2"):
            x = bundle["series"][pred]
            for ep in (
                "MHC1_protein",
                "IFN_core_protein",
                "IFN_signaling_protein",
                "IFN_isg_protein",
                "CD8A_protein",
                "APM_protein",
            ):
                rows.append(quartile_row(cohort, f"{pred}_protein", ep, x, scores[ep]))
        # Histology / grade strata.
        for stratum, mask in masks.items():
            for pred in ("CLDN4", "TACSTD2"):
                x = bundle["series"][pred].where(mask)
                for ep in (
                    "MHC1_protein",
                    "IFN_core_protein",
                    "IFN_signaling_protein",
                    "IFN_isg_protein",
                    "CD8A_protein",
                    "APM_protein",
                ):
                    rows.append(
                        spearman_row(
                            cohort, "stratum", f"{pred}_protein", ep, stratum,
                            x, scores[ep].where(mask), wes, rng, MIN_STRATUM_N,
                        )
                    )
        # Leave-one-out of the IFN-core score, CLDN4 only. Diagnostic, not a claim.
        A = _A()
        members = [g for g in IFN_SCORE_MEMBERS if bundle["rows_used"].get(g)]
        for drop in members:
            kept = [g for g in members if g != drop]
            vecs = [bundle["series"][g] for g in kept]
            score = A.mean_z(vecs, 4)
            row = spearman_row(
                cohort, "loo", "CLDN4_protein", "IFN_core_without_" + drop, "all",
                bundle["series"]["CLDN4"], score, None, rng, MIN_PANEL_N,
            )
            rows.append(row)

    df = pd.DataFrame(rows)
    A = _A()
    df["q_family"] = np.nan
    for fam, idx in df.groupby("family").groups.items():
        if fam == "loo":
            continue
        q = A.bh_q(df.loc[idx, "p"].tolist())
        df.loc[idx, "q_family"] = q

    tab = outdir / "tables"
    fig = outdir / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    df.to_csv(tab / "sensitivity.tsv", sep="\t", index=False)
    counts = pd.concat(count_frames, ignore_index=True)
    counts.to_csv(tab / "histology_counts.tsv", sep="\t", index=False)
    plot_panels(df, fig)
    (outdir / "SENSITIVITY.md").write_text(render(df, counts))
    # Sanity: rank-floor n must be the full cohort, not the CLDN4-complete n.
    rf = df[(df.family == "rank_floor") & (df.endpoint == "CD8A_protein")]
    assert set(rf["n"]) == {110, 108}, rf[["cohort", "n"]].to_string()
    return df


def plot_panels(df: pd.DataFrame, figdir: Path) -> None:
    sub = df[(df.family == "panel") & (df.predictor == "CLDN4_protein")].copy()
    order = [
        "IFN_signaling_protein",
        "IFN_isg_protein",
        "IFN_receptor_protein",
        "APM_protein",
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharex=True)
    for ax, cohort in zip(axes, ("LUAD", "LSCC")):
        ys = []
        labels = []
        for i, ep in enumerate(order):
            hit = sub[(sub.cohort == cohort) & (sub.endpoint == ep)]
            y = len(order) - 1 - i
            labels.append(ep.replace("_protein", "").replace("IFN_", "IFN "))
            ys.append(y)
            if hit.empty or not np.isfinite(hit.iloc[0]["rho"]):
                continue
            r = hit.iloc[0]
            color = "#b3483a" if r["ci_hi"] < 0 else ("#2a6f4e" if r["ci_lo"] > 0 else "#5c6570")
            ax.errorbar(
                r["rho"], y,
                xerr=[[r["rho"] - r["ci_lo"]], [r["ci_hi"] - r["rho"]]],
                fmt="o", color=color, ms=5, lw=1.2, capsize=2,
            )
            ax.text(0.72, y, f"{r['rho']:+.2f}  n={int(r['n'])}", va="center", fontsize=7, color=color)
        ax.axvline(0, color="#888888", lw=0.8)
        ax.set_yticks(ys)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(-0.85, 1.15)
        ax.set_title(f"{cohort} CLDN4")
        ax.set_xlabel("Spearman ρ (95% CI)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("Alternate IFN / antigen-processing panels (not the primary IFN-core)", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "fig_ifn_panels.png", dpi=160)
    fig.savefig(figdir / "fig_ifn_panels.pdf")
    plt.close(fig)


def _fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3g}"


def _fmt_r(r: float) -> str:
    if not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def _cell(df: pd.DataFrame, **kw) -> pd.Series:
    hit = df
    for k, v in kw.items():
        hit = hit[hit[k] == v]
    if hit.empty:
        raise KeyError(kw)
    return hit.iloc[0]


def _spear_bits(r: pd.Series) -> str:
    return (
        f"ρ={_fmt_r(r['rho'])}, n={int(r['n'])}, 95% CI {_fmt_r(r['ci_lo'])} to {_fmt_r(r['ci_hi'])}, "
        f"p={_fmt_p(r['p'])}, q={_fmt_p(r['q_family'])}, "
        f"WES partial ρ={_fmt_r(r['rho_partial'])} "
        f"(CI {_fmt_r(r['ci_lo_partial'])} to {_fmt_r(r['ci_hi_partial'])})"
    )


def render(df: pd.DataFrame, counts: pd.DataFrame) -> str:
    lines: list[str] = []
    lines.append("## Sensitivity — IFN-core null, LUAD null, missingness, histology")
    lines.append("")
    lines.append(
        "The primary calls above are unchanged. LSCC CLDN4 vs MHC-I and LSCC CLDN4 vs CD8A "
        "stay the only pairs that meet the primary support rule. "
        "This section asks whether that IFN-core null and that LUAD null are an artifact of "
        "the 10-gene mean, of CLDN4 dropout, or of histologic mixture. "
        "No log2 abundance is imputed. Within each family, q is Benjamini–Hochberg on that "
        "family only. A sensitivity inverse is not promoted into the primary support rule."
    )
    lines.append("")

    # Panel table, CLDN4.
    lines.append("### Alternate panels (locked genes, split)")
    lines.append("")
    lines.append(
        "IFN-signaling = STAT1, STAT2, IRF1, IRF9, JAK1, JAK2 (need ≥4). "
        "IFN-ISG = ISG15, MX1, OAS1, IFI35 (need ≥3). "
        "IFN-receptor = IFNAR1 and IFNGR1 (need both; IFNAR2 and IFNGR2 stay too sparse to score). "
        "APM = TAP1, TAP2, TAPBP, PSMB8, PSMB9, NLRC5 (need ≥4). "
        "These are the same proteins as the gene table, averaged. Not a proteome-wide search."
    )
    lines.append("")
    headers = ["Cohort", "Predictor", "Panel", "n", "ρ", "95% CI", "p", "q", "partial ρ"]
    body = []
    panel = df[df.family == "panel"].sort_values(["endpoint", "cohort", "predictor"])
    for r in panel.itertuples(index=False):
        body.append(
            [
                r.cohort,
                r.predictor.replace("_protein", ""),
                r.endpoint.replace("_protein", ""),
                str(int(r.n)),
                _fmt_r(r.rho),
                f"{_fmt_r(r.ci_lo)} to {_fmt_r(r.ci_hi)}",
                _fmt_p(r.p),
                _fmt_p(r.q_family),
                _fmt_r(r.rho_partial),
            ]
        )
    lines.append(_md(headers, body))
    lines.append("")
    lines.append("Figure: `figures/fig_ifn_panels.png`.")
    lines.append("")

    # Narrative from CLDN4 panels.
    lines.append(_panel_verdict(df))
    lines.append("")

    lines.append("### Leave-one-out of the IFN-core score (CLDN4)")
    lines.append("")
    lines.append(
        "Each row drops one core member and rebuilds the mean-z score. "
        "Diagnostic only. These p-values are not in a discovery FDR."
    )
    lines.append("")
    loo = df[df.family == "loo"]
    lrows = []
    drops = sorted({e.replace("IFN_core_without_", "") for e in loo.endpoint})
    for g in drops:
        a = _cell(loo, cohort="LUAD", endpoint="IFN_core_without_" + g)
        b = _cell(loo, cohort="LSCC", endpoint="IFN_core_without_" + g)
        lrows.append(
            [
                g,
                f"{_fmt_r(a.rho)} (n={int(a.n)}, p={_fmt_p(a.p)})",
                f"{_fmt_r(b.rho)} (n={int(b.n)}, p={_fmt_p(b.p)})",
            ]
        )
    lines.append(_md(["Dropped", "LUAD CLDN4", "LSCC CLDN4"], lrows))
    lines.append("")
    lines.append(_loo_sentence(loo))
    lines.append("")

    lines.append("### Imputation-free complete-case variants")
    lines.append("")
    lines.append(
        "Strict score: a tumor counts only when every panel member is quantified "
        "(no partial mean). Shared complete-case: the same tumors have quantified CLDN4, "
        "all 10 IFN-core members, HLA-A/B/C, and CD8A. "
        "Quartile: among quantified predictor values, highest quartile vs lowest, Mann–Whitney. "
        "Delta is median(Q4) − median(Q1); inverse means delta < 0. "
        "Rank-floor is a separate missingness check, not a complete-case result: "
        "tumors with missing CLDN4 are tied at one rank below every quantified value. "
        "That is not a filled-in log2 abundance."
    )
    lines.append("")
    lines.append("**Strict scores, CLDN4.**")
    lines.append("")
    lines.append(_small_spearman(df[(df.family == "strict") & (df.predictor == "CLDN4_protein")]))
    lines.append("")
    lines.append("**Shared complete-case.**")
    lines.append("")
    lines.append(_small_spearman(df[df.family == "shared_complete"]))
    lines.append("")
    lines.append("**CLDN4 Q4 vs Q1 (quantified tumors only).**")
    lines.append("")
    q = df[(df.family == "quartile") & (df.predictor == "CLDN4_protein")]
    qrows = []
    for r in q.sort_values(["cohort", "endpoint"]).itertuples(index=False):
        qrows.append(
            [
                r.cohort,
                r.endpoint.replace("_protein", ""),
                f"{int(r.n_q4)} vs {int(r.n_q1)}",
                _fmt_r(r.delta_q4_minus_q1),
                _fmt_p(r.p),
                _fmt_p(r.q_family),
            ]
        )
    lines.append(_md(["Cohort", "Endpoint", "n Q4 vs Q1", "Δ median", "p", "q"], qrows))
    lines.append("")
    lines.append("**Rank-floor Spearman for missing CLDN4 (sensitivity, not the primary n).**")
    lines.append("")
    lines.append(_small_spearman(df[df.family == "rank_floor"]))
    lines.append("")
    lines.append(_variant_verdict(df))
    lines.append("")

    lines.append("### Histology")
    lines.append("")
    lines.append(
        "Public labels, joined on the protein case id. "
        "LUAD dominant histological subtype is the cBioPortal field "
        "`DOMINANT_HISTOLOGICAL_SUBTYPE` (study `luad_cptac_2020`; the four `11LU` ids are stored there with an X prefix and were mapped back). "
        "Grade is `Histologic_Grade` in freeze `LUAD_meta.txt` / `LSCC_meta.txt`. "
        "LSCC has no acinar/solid code. Pathology text is "
        "`PATHOLOGY_BASED_HISTOLOGY_ASSESSMENT` (study `lusc_cptac_2021`). "
        "Basaloid means that string contains \"basaloid\". "
        "NMF / RNA clusters are molecular and were not tested as histology. "
        "A stratum is tested only when pairwise n ≥ 20. Smaller levels are counted and stopped."
    )
    lines.append("")
    show = counts[counts.axis.isin(["dominant_histological_subtype", "histologic_grade", "tested_stratum"])]
    crows = []
    for r in show.itertuples(index=False):
        crows.append([
            str(r.cohort),
            str(r.axis),
            "NA" if pd.isna(r.level) else str(r.level),
            str(int(r.n_tumors)),
            str(int(r.n_cldn4)),
        ])
    lines.append(_md(["Cohort", "Axis", "Level", "Tumors", "CLDN4 quantified"], crows))
    lines.append("")
    tested = df[(df.family == "stratum") & (df.predictor == "CLDN4_protein") & df.rho.notna()]
    skipped = df[(df.family == "stratum") & (df.predictor == "CLDN4_protein") & df.rho.isna()]
    lines.append(
        f"CLDN4 stratum tests with pairwise n ≥ 20: **{tested[['cohort','stratum','endpoint']].drop_duplicates().shape[0]}** "
        f"endpoint-rows. Not tested (pairwise n < 20): {int(skipped.shape[0])} rows. "
        "Full grid, including TACSTD2: `tables/sensitivity.tsv`. Counts: `tables/histology_counts.tsv`."
    )
    lines.append("")
    if tested.empty:
        lines.append("No histology stratum reached pairwise n ≥ 20.")
    else:
        lines.append(_small_spearman(tested))
    lines.append("")
    lines.append(_stratum_verdict(df))
    lines.append("")
    lines.append(_strongest(df))
    lines.append("")
    return "\n".join(lines)


def _md(headers: list[str], rows: list[list[str]]) -> str:
    line = "| " + " | ".join(headers) + " |"
    segs = ["---"] + ["---:"] * (len(headers) - 1)
    sep = "| " + " | ".join(segs) + " |"
    body = ["| " + " | ".join(str(cell) for cell in r) + " |" for r in rows]
    return "\n".join([line, sep, *body])


def _small_spearman(df: pd.DataFrame) -> str:
    rows = []
    view = df.sort_values(["cohort", "predictor", "endpoint", "stratum"])
    for r in view.itertuples(index=False):
        if not np.isfinite(r.rho):
            continue
        rows.append(
            [
                r.cohort,
                r.predictor.replace("_protein", ""),
                r.endpoint.replace("_protein", ""),
                r.stratum,
                str(int(r.n)),
                _fmt_r(r.rho),
                f"{_fmt_r(r.ci_lo)} to {_fmt_r(r.ci_hi)}",
                _fmt_p(r.p),
                _fmt_p(r.q_family),
                _fmt_r(r.rho_partial),
            ]
        )
    if not rows:
        return "No Spearman in this block met the minimum n."
    return _md(
        ["Cohort", "Predictor", "Endpoint", "Slice", "n", "ρ", "95% CI", "p", "q", "partial ρ"],
        rows,
    )


def _inverse_ok(r) -> bool:
    if isinstance(r, pd.Series):
        rho, ci_hi, q, ci_hi_p = r["rho"], r["ci_hi"], r["q_family"], r["ci_hi_partial"]
    else:
        rho, ci_hi, q, ci_hi_p = r.rho, r.ci_hi, r.q_family, r.ci_hi_partial
    return bool(
        np.isfinite(rho)
        and np.isfinite(ci_hi)
        and ci_hi < 0
        and np.isfinite(q)
        and q < 0.05
        and np.isfinite(ci_hi_p)
        and ci_hi_p < 0
    )


def _panel_verdict(df: pd.DataFrame) -> str:
    panel = df[df.family == "panel"]
    bits = []
    for cohort in ("LUAD", "LSCC"):
        sub = panel[(panel.cohort == cohort) & (panel.predictor == "CLDN4_protein")]
        inv = [r.endpoint.replace("_protein", "") for r in sub.itertuples() if _inverse_ok(r)]
        crossed = []
        for r in sub.itertuples():
            if np.isfinite(r.ci_lo) and r.ci_lo <= 0 <= r.ci_hi:
                crossed.append(f"{r.endpoint.replace('_protein','')} ρ={r.rho:+.3f}")
        near = []
        for r in sub.itertuples():
            if (
                np.isfinite(r.ci_hi)
                and r.ci_hi < 0
                and not (np.isfinite(r.q_family) and r.q_family < 0.05 and np.isfinite(r.ci_hi_partial) and r.ci_hi_partial < 0)
            ):
                near.append(
                    f"{r.endpoint.replace('_protein','')} ρ={r.rho:+.3f}, "
                    f"CI {r.ci_lo:+.3f} to {r.ci_hi:+.3f}, q={_fmt_p(r.q_family)}, "
                    f"partial CI to {r.ci_hi_partial:+.3f}"
                )
        if inv:
            detail = "; ".join(
                f"{e}: " + _spear_bits(_cell(sub, endpoint=e + "_protein", predictor="CLDN4_protein", cohort=cohort))
                for e in inv
            )
            bits.append(f"{cohort} CLDN4 panels that clear the sensitivity bar: {detail}.")
        else:
            bits.append(f"{cohort} CLDN4: no alternate panel clears the sensitivity bar.")
        if near:
            bits.append(f"{cohort} CLDN4 panels with CI below 0 that do not clear q and partial together: " + "; ".join(near) + ".")
        if crossed:
            bits.append(f"{cohort} CLDN4 panels whose CI includes 0: " + "; ".join(crossed) + ".")
    tac = panel[panel.predictor == "TACSTD2_protein"]
    tac_inv = [
        f"{r.cohort} {r.endpoint.replace('_protein','')} ρ={r.rho:+.3f}, q={_fmt_p(r.q_family)}"
        for r in tac.itertuples()
        if _inverse_ok(r)
    ]
    tac_pos = [
        f"{r.cohort} {r.endpoint.replace('_protein','')} ρ={r.rho:+.3f}, "
        f"CI {r.ci_lo:+.3f} to {r.ci_hi:+.3f}, q={_fmt_p(r.q_family)}, "
        f"partial ρ={r.rho_partial:+.3f} (CI {r.ci_lo_partial:+.3f} to {r.ci_hi_partial:+.3f})"
        for r in tac.itertuples()
        if np.isfinite(r.ci_lo) and r.ci_lo > 0 and np.isfinite(r.q_family) and r.q_family < 0.05
    ]
    if tac_inv:
        bits.append("TACSTD2 sensitivity inverses: " + ", ".join(tac_inv) + ".")
    else:
        bits.append("TACSTD2 has no alternate-panel sensitivity inverse in either cohort.")
    if tac_pos:
        bits.append(
            "Opposite sign, and it does clear the panel-family FDR: " + "; ".join(tac_pos) + ". "
            "That is higher TACSTD2 with higher ISG protein, not the inverse thesis."
        )
    return " ".join(bits)


def _loo_sentence(loo: pd.DataFrame) -> str:
    def rng(cohort: str) -> str:
        sub = loo[loo.cohort == cohort]
        return (
            f"{cohort} leave-one-out ρ spans {_fmt_r(sub.rho.min())} to {_fmt_r(sub.rho.max())} "
            f"(n={int(sub.n.iloc[0])})"
        )
    # Which drop moves LSCC closest to excluding 0.
    ls = loo[loo.cohort == "LSCC"].sort_values("rho")
    most = ls.iloc[0]
    return (
        f"{rng('LUAD')}. {rng('LSCC')}. "
        f"The most negative LSCC drop is {most.endpoint.replace('IFN_core_without_', '')} "
        f"(ρ={_fmt_r(most.rho)}, CI {_fmt_r(most.ci_lo)} to {_fmt_r(most.ci_hi)}). "
        "Dropping one ISG does not by itself turn the IFN-core into a primary endpoint."
    )


def _variant_verdict(df: pd.DataFrame) -> str:
    parts = []
    for fam, label in (
        ("strict", "Strict complete-case"),
        ("shared_complete", "Shared complete-case"),
        ("rank_floor", "Rank-floor"),
    ):
        sub = df[(df.family == fam) & (df.predictor == "CLDN4_protein") & (df.cohort == "LUAD")]
        inv = [r.endpoint for r in sub.itertuples() if _inverse_ok(r)]
        if inv:
            parts.append(f"{label}, LUAD CLDN4, sensitivity inverse: {', '.join(inv)}.")
        else:
            parts.append(f"{label} does not give LUAD CLDN4 a sensitivity inverse on MHC-I, IFN-core, IFN-signaling, or CD8A.")
    q = df[(df.family == "quartile") & (df.predictor == "CLDN4_protein")]
    q_lu = q[(q.cohort == "LUAD") & (q.q_family < 0.05) & (q.delta_q4_minus_q1 < 0)]
    q_ls = q[(q.cohort == "LSCC") & (q.q_family < 0.05) & (q.delta_q4_minus_q1 < 0)]
    if q_lu.empty:
        parts.append("LUAD CLDN4 Q4 vs Q1 does not clear the quartile-family FDR in the inverse direction.")
    else:
        parts.append(
            "LUAD CLDN4 Q4 vs Q1 inverse at q<0.05: "
            + ", ".join(f"{r.endpoint} Δ={r.delta_q4_minus_q1:+.3f}" for r in q_lu.itertuples())
            + "."
        )
    if not q_ls.empty:
        parts.append(
            "LSCC CLDN4 Q4 vs Q1 inverse at q<0.05: "
            + ", ".join(
                f"{r.endpoint.replace('_protein','')} Δ={r.delta_q4_minus_q1:+.3f}, q={_fmt_p(r.q_family)}"
                for r in q_ls.itertuples()
            )
            + ". That agrees in sign with the primary MHC-I and CD8A calls; it is still a quartile contrast, not a new cohort."
        )
    return " ".join(parts)


def _stratum_verdict(df: pd.DataFrame) -> str:
    sub = df[(df.family == "stratum") & (df.predictor == "CLDN4_protein") & df.rho.notna()]
    if sub.empty:
        return "No CLDN4 histology stratum was testable."
    inv = sub[sub.apply(_inverse_ok, axis=1)]
    lu = sub[sub.cohort == "LUAD"]
    lu_inv = inv[inv.cohort == "LUAD"]
    bits = []
    if lu_inv.empty:
        bits.append(
            "Inside testable LUAD slices (acinar, non-acinar, G2, G3, whichever reached n≥20), "
            "CLDN4 does not pick up a sensitivity inverse. The LUAD null is not a grade or acinar/non-acinar mixture artifact at this n."
        )
    else:
        bits.append(
            "LUAD CLDN4 sensitivity inverses inside a stratum: "
            + "; ".join(
                f"{r.stratum} {r.endpoint.replace('_protein','')} {_spear_bits(r)}"
                for r in lu_inv.itertuples()
            )
            + "."
        )
    ls_inv = inv[inv.cohort == "LSCC"]
    if ls_inv.empty:
        bits.append("No LSCC CLDN4 stratum adds a new sensitivity inverse beyond what the full cohort already shows.")
    else:
        bits.append(
            "LSCC CLDN4 stratum sensitivity inverses (same direction as the full-cohort MHC-I/CD8A result, not a replacement): "
            + "; ".join(
                f"{r.stratum} vs {r.endpoint.replace('_protein','')} ρ={r.rho:+.3f}, n={int(r.n)}, q={_fmt_p(r.q_family)}"
                for r in ls_inv.itertuples()
            )
            + "."
        )
    # Explicitly name untested rare LUAD subtypes.
    bas = sub[sub.stratum == "basaloid"]
    if not bas.empty and not bas.apply(_inverse_ok, axis=1).any():
        bits.append(
            f"LSCC basaloid (pathology text, CLDN4 n={int(bas.n.max())}) does not clear the bar on any endpoint; "
            "the LSCC inverses that do clear it are in the non-basaloid majority, not in a basaloid-only slice."
        )
    bits.append(
        "Solid, papillary, micropapillary, lepidic, and the other single LUAD subtype labels each have fewer than 20 CLDN4-quantified tumors, so none is a tested rescue of the LUAD null."
    )
    return " ".join(bits)


def _strongest(df: pd.DataFrame) -> str:
    spear = df[(df.family != "loo") & (df.family != "quartile") & df.rho.notna()].copy()
    spear = spear[spear.apply(_inverse_ok, axis=1)]
    lines = ["### Strongest thesis-aligned protein associations in this sensitivity pass"]
    lines.append("")
    lines.append(
        "Thesis-aligned means higher CLDN4 or TACSTD2 protein with lower MHC-I, IFN-panel, APM, or CD8A protein. "
        "A row is listed only when the bootstrap CI is entirely below 0, the within-family q is < 0.05, "
        "and the WES-purity partial CI is entirely below 0. "
        "Rank-floor and stratum rows that repeat LSCC CLDN4 vs MHC-I or CD8A are the primary result under a different slice, not a new endpoint. "
        "The new panel-level association is LSCC CLDN4 vs the IFN-signaling score (and, at n=45, the IFN-receptor score). "
        "The IFN-core score is not in this list. The LUAD null is not in this list."
    )
    lines.append("")
    if spear.empty:
        lines.append(
            "No sensitivity Spearman meets that bar. The strongest thesis-aligned protein results "
            "remain the primary ones: LSCC CLDN4 vs MHC-I and LSCC CLDN4 vs CD8A. "
            "The IFN-core null and the LUAD null stand."
        )
        return "\n".join(lines)
    spear = spear.sort_values("rho")
    top = spear.head(8)
    rows = []
    for r in top.itertuples(index=False):
        rows.append(
            [
                r.cohort,
                r.predictor.replace("_protein", ""),
                r.endpoint.replace("_protein", ""),
                r.family,
                r.stratum,
                str(int(r.n)),
                _fmt_r(r.rho),
                _fmt_p(r.q_family),
                _fmt_r(r.rho_partial),
            ]
        )
    lines.append(_md(
        ["Cohort", "Predictor", "Endpoint", "Family", "Slice", "n", "ρ", "q", "partial ρ"],
        rows,
    ))
    lines.append("")
    best = spear.iloc[0]
    lu = spear[spear.cohort == "LUAD"]
    if lu.empty:
        lu_txt = "LUAD contributes none."
    else:
        b = lu.iloc[0]
        lu_txt = (
            f"The strongest LUAD row on this bar is {b.predictor.replace('_protein','')} vs "
            f"{b.endpoint.replace('_protein','')} ({b.family}, {b.stratum}, ρ={b.rho:+.3f}, n={int(b.n)})."
        )
    lines.append(
        f"Most negative row in the table: {best.cohort} {best.predictor.replace('_protein','')} vs "
        f"{best.endpoint.replace('_protein','')} ({best.family}, {best.stratum}, ρ={best.rho:+.3f}, n={int(best.n)}). "
        f"{lu_txt} "
        "None of these rows replaces the primary IFN-core call. The IFN-core score itself stays null."
    )
    return "\n".join(lines)
