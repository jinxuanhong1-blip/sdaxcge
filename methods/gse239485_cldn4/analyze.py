#!/usr/bin/env python3
"""GSE239485 LLC bulk: Cldn4 and IFN/MHC on the Tacstd2 +2.09 contrast.

Contrast (given, not re-audited as a discovery): Poly I:C + anti-PD-1 vs vehicle.
Tacstd2 log2FC +2.09 is taken as given and re-measured only to lock the same
samples. TISMO 49/64 is taken as given.

Public processed matrix only: GSE239485_Processed_data.xlsx (DataNorm).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CACHE = Path("/tmp/gse239485")
XLSX = CACHE / "GSE239485_Processed_data.xlsx"
FIG = HERE / "figures"
TAB = HERE / "tables"
FIG.mkdir(exist_ok=True)
TAB.mkdir(exist_ok=True)

# Human leftover lists → mouse symbols on this matrix.
# IFN: leftover IFNG/STAT1/CXCL9/CXCL10/IDO1/HLA-DRA
IFN_GENES = ["Ifng", "Stat1", "Cxcl9", "Cxcl10", "Ido1", "H2-Aa"]
# MHC-I: leftover HLA-A/B/C / B2M / TAP1 / TAP2
MHC_GENES = ["H2-K1", "H2-D1", "H2-Q4", "B2m", "Tap1", "Tap2"]
# Compact Ayers-like (same 6 as leftover IFN; H2-Aa = HLA-DRA)
AYERS6 = IFN_GENES
CONTROLS = ["Cd8a", "Cd274", "Ptprc"]
TARGETS = ["Tacstd2", "Cldn4"]

# Ensembl locks (GRCm38/39).
ENSEMBL = {
    "Tacstd2": "ENSMUSG00000051397",
    "Cldn4": "ENSMUSG00000047501",
    "Ifng": "ENSMUSG00000055170",
    "Stat1": "ENSMUSG00000026104",
    "Cxcl9": "ENSMUSG00000029417",
    "Cxcl10": "ENSMUSG00000034855",
    "Ido1": "ENSMUSG00000031551",
    "H2-Aa": "ENSMUSG00000036594",
    "H2-K1": "ENSMUSG00000061232",
    "H2-D1": "ENSMUSG00000073411",
    "H2-Q4": "ENSMUSG00000035929",
    "B2m": "ENSMUSG00000060802",
    "Tap1": "ENSMUSG00000037321",
    "Tap2": "ENSMUSG00000024339",
    "Cd8a": "ENSMUSG00000053977",
    "Cd274": "ENSMUSG00000016496",
    "Ptprc": "ENSMUSG00000026395",
}


def load_matrix(path: Path) -> pd.DataFrame:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sh = wb[wb.sheetnames[0]]
    rows = [[c for c in r] for r in sh.iter_rows(values_only=True)]
    wb.close()
    header = [str(x) if x is not None else "" for x in rows[0]]
    body = rows[1:]
    df = pd.DataFrame(body, columns=header)
    sample_cols = [c for c in header if c.startswith(("C_", "D_", "T_"))]
    df["ensembl"] = df["gene_id"].astype(str).str.split(".").str[0]
    df["symbol"] = df["gene_name"].astype(str)
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df, sample_cols


def arm(col: str) -> str:
    if col.startswith("C_"):
        return "vehicle"
    if col.startswith("D_"):
        return "polyIC_aPD1"
    if col.startswith("T_"):
        return "triplet"
    return "other"


def pick_gene(df: pd.DataFrame, symbol: str) -> pd.Series | None:
    # Symbol first so a stale Ensembl ID cannot silently swap genes.
    hit = df[df["symbol"].str.lower() == symbol.lower()]
    if len(hit) == 0:
        ens = ENSEMBL.get(symbol)
        if ens is not None:
            hit = df[df["ensembl"] == ens]
    if len(hit) == 0:
        return None
    pc = hit[hit["gene_type"] == "protein_coding"]
    row = (pc if len(pc) else hit).iloc[0]
    ens = ENSEMBL.get(symbol)
    if ens and row["ensembl"] != ens:
        raise SystemExit(f"{symbol}: matrix {row['ensembl']} != locked {ens}")
    return row


def welch_mwu(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_a": int(a.size),
        "n_b": int(b.size),
        "mean_a": float(np.mean(a)) if a.size else None,
        "mean_b": float(np.mean(b)) if b.size else None,
        "median_a": float(np.median(a)) if a.size else None,
        "median_b": float(np.median(b)) if b.size else None,
        "log2fc": None,
        "welch_p": None,
        "mwu_p": None,
    }
    if a.size == 0 or b.size == 0:
        return out
    out["log2fc"] = float(np.mean(b) - np.mean(a))
    if a.size >= 2 and b.size >= 2:
        out["welch_p"] = float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue)
        out["mwu_p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    return out


def mean_z(mat: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    used = [g for g in genes if g in mat.index]
    if not used:
        return pd.Series(np.nan, index=mat.columns), []
    z = mat.loc[used].apply(lambda r: (r - r.mean()) / (r.std(ddof=0) or 1.0), axis=1)
    return z.mean(axis=0), used


def cliff_delta(a, b):
    """Cliff δ = P(b>a) - P(b<a); positive = higher in b."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return None
    gt = sum(1 for x in b for y in a if x > y)
    lt = sum(1 for x in b for y in a if x < y)
    return float((gt - lt) / (a.size * b.size))


def main():
    if not XLSX.exists():
        raise SystemExit(f"missing {XLSX}; download GSE239485_Processed_data.xlsx")

    df, sample_cols = load_matrix(XLSX)
    veh = [c for c in sample_cols if arm(c) == "vehicle"]
    combo = [c for c in sample_cols if arm(c) == "polyIC_aPD1"]
    trip = [c for c in sample_cols if arm(c) == "triplet"]
    assert len(veh) == 8 and len(combo) == 8 and len(trip) == 8, (len(veh), len(combo), len(trip))

    wanted = TARGETS + IFN_GENES + MHC_GENES + CONTROLS
    gene_rows = {}
    coverage = []
    for g in wanted:
        row = pick_gene(df, g)
        if row is None:
            coverage.append({"gene": g, "found": False, "ensembl": ENSEMBL.get(g), "symbol_on_matrix": None})
            continue
        gene_rows[g] = row
        coverage.append(
            {
                "gene": g,
                "found": True,
                "ensembl": row["ensembl"],
                "symbol_on_matrix": row["symbol"],
                "gene_type": row["gene_type"],
            }
        )

    expr = pd.DataFrame({g: gene_rows[g][sample_cols].astype(float) for g in gene_rows})
    expr.index = sample_cols
    expr["arm"] = [arm(c) for c in sample_cols]
    expr["arm_label"] = expr["arm"].map(
        {
            "vehicle": "vehicle",
            "polyIC_aPD1": "Poly I:C + anti-PD-1",
            "triplet": "Poly I:C + anti-PD-1 + C5aR1i",
        }
    )

    ifn_mat = pd.DataFrame({g: expr[g] for g in IFN_GENES if g in expr.columns}).T
    mhc_mat = pd.DataFrame({g: expr[g] for g in MHC_GENES if g in expr.columns}).T
    ifn_score, ifn_used = mean_z(ifn_mat, list(ifn_mat.index))
    mhc_score, mhc_used = mean_z(mhc_mat, list(mhc_mat.index))
    expr["IFN"] = ifn_score
    expr["MHC"] = mhc_score

    contrasts = []

    def add(feature, a_cols, b_cols, design, group_a, group_b, note=""):
        a = expr.loc[a_cols, feature].astype(float)
        b = expr.loc[b_cols, feature].astype(float)
        st = welch_mwu(a, b)
        st.update(
            {
                "feature": feature,
                "design": design,
                "group_a": group_a,
                "group_b": group_b,
                "n": f"{st['n_a']}/{st['n_b']}",
                "cliff_delta": cliff_delta(a, b),
                "note": note,
            }
        )
        contrasts.append(st)

    primary_note = "primary contrast; no aPD-1 monotherapy arm; Poly I:C on every treated sample"
    for feat in list(expr.columns):
        if feat in ("arm", "arm_label"):
            continue
        add(feat, veh, combo, "Poly I:C + anti-PD-1 vs vehicle", "vehicle", "polyIC_aPD1", primary_note)
        add(feat, veh, trip, "Poly I:C + anti-PD-1 + C5aR1i vs vehicle", "vehicle", "triplet", "triple combo; still no aPD-1 monotherapy")

    # Spearman Cldn4 / Tacstd2 vs IFN/MHC/Cd8a on all 24 and on primary 16.
    corrs = []
    for subset_name, cols in (("all_24", sample_cols), ("primary_16", veh + combo)):
        sub = expr.loc[cols]
        for g in ["Cldn4", "Tacstd2"]:
            for axis in ["IFN", "MHC", "Cd8a", "Ifng", "Cd274"]:
                if g not in sub.columns or axis not in sub.columns:
                    continue
                rho, p = stats.spearmanr(sub[g], sub[axis], nan_policy="omit")
                corrs.append(
                    {
                        "subset": subset_name,
                        "n": int(sub[[g, axis]].dropna().shape[0]),
                        "gene": g,
                        "axis": axis,
                        "rho": float(rho),
                        "p": float(p),
                    }
                )

    # Partial Spearman Cldn4 vs IFN/MHC | Cd8a (infiltrate proxy; no ESTIMATE list here).
    partials = []
    for subset_name, cols in (("all_24", sample_cols), ("primary_16", veh + combo)):
        sub = expr.loc[cols].dropna(subset=["Cldn4", "IFN", "MHC", "Cd8a"])
        if len(sub) < 6:
            continue
        for axis in ["IFN", "MHC"]:
            # residualize both on Cd8a
            x = sub["Cldn4"].to_numpy()
            y = sub[axis].to_numpy()
            z = sub["Cd8a"].to_numpy()
            zx = np.vstack([np.ones(len(z)), z]).T
            bx, *_ = np.linalg.lstsq(zx, x, rcond=None)
            by, *_ = np.linalg.lstsq(zx, y, rcond=None)
            rx = x - zx @ bx
            ry = y - zx @ by
            rho, p = stats.spearmanr(rx, ry)
            partials.append(
                {
                    "subset": subset_name,
                    "n": int(len(sub)),
                    "gene": "Cldn4",
                    "axis": axis,
                    "covariate": "Cd8a",
                    "rho": float(rho),
                    "p": float(p),
                }
            )

    contrast_df = pd.DataFrame(contrasts)
    corr_df = pd.DataFrame(corrs)
    part_df = pd.DataFrame(partials)
    cov_df = pd.DataFrame(coverage)

    expr.to_csv(TAB / "per_sample.tsv", sep="\t")
    contrast_df.to_csv(TAB / "contrasts.tsv", sep="\t", index=False)
    corr_df.to_csv(TAB / "spearman.tsv", sep="\t", index=False)
    part_df.to_csv(TAB / "partial_spearman.tsv", sep="\t", index=False)
    cov_df.to_csv(TAB / "gene_coverage.tsv", sep="\t", index=False)

    def row(feature, design="Poly I:C + anti-PD-1 vs vehicle"):
        hit = contrast_df[(contrast_df["feature"] == feature) & (contrast_df["design"] == design)]
        return hit.iloc[0].to_dict() if len(hit) else None

    tac = row("Tacstd2")
    cld = row("Cldn4")
    ifn = row("IFN")
    mhc = row("MHC")
    cd8 = row("Cd8a")
    ifng = row("Ifng")

    one = {
        "dataset": "GSE239485",
        "model": "LLC subcutaneous syngeneic",
        "contrast": "Poly I:C + anti-PD-1 vs vehicle",
        "n_vehicle": 8,
        "n_treated": 8,
        "honest_n": "8/8",
        "Tacstd2_log2FC_given": 2.09,
        "Tacstd2_log2FC": tac["log2fc"] if tac else None,
        "Tacstd2_welch_p": tac["welch_p"] if tac else None,
        "Cldn4_log2FC": cld["log2fc"] if cld else None,
        "Cldn4_welch_p": cld["welch_p"] if cld else None,
        "Cldn4_mwu_p": cld["mwu_p"] if cld else None,
        "IFN_meanz_delta": ifn["log2fc"] if ifn else None,
        "IFN_welch_p": ifn["welch_p"] if ifn else None,
        "IFN_genes": ",".join(ifn_used),
        "IFN_n_genes": len(ifn_used),
        "MHC_meanz_delta": mhc["log2fc"] if mhc else None,
        "MHC_welch_p": mhc["welch_p"] if mhc else None,
        "MHC_genes": ",".join(mhc_used),
        "MHC_n_genes": len(mhc_used),
        "Cd8a_log2FC": cd8["log2fc"] if cd8 else None,
        "Cd8a_welch_p": cd8["welch_p"] if cd8 else None,
        "Ifng_log2FC": ifng["log2fc"] if ifng else None,
        "Ifng_welch_p": ifng["welch_p"] if ifng else None,
        "TISMO_given": "49/64 Tacstd2 up after ICB",
        "no_R_vs_NR": True,
        "no_aPD1_monotherapy": True,
    }
    pd.DataFrame([one]).to_csv(TAB / "one_row.tsv", sep="\t", index=False)
    (TAB / "summary.json").write_text(json.dumps({"one_row": one, "ifn_used": ifn_used, "mhc_used": mhc_used}, indent=2) + "\n")

    # Figures
    def box3(ax, feature, title, ylab):
        groups = [expr.loc[veh, feature], expr.loc[combo, feature], expr.loc[trip, feature]]
        labels = ["vehicle\nn=8", "Poly I:C\n+ aPD-1\nn=8", "triplet\nn=8"]
        colors = ["#9aa0a6", "#1a73e8", "#8ab4f8"]
        bp = ax.boxplot(groups, tick_labels=labels, widths=0.55, showfliers=False, patch_artist=True)
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.45)
        for i, g in enumerate(groups, start=1):
            x = np.random.default_rng(0).normal(i, 0.06, size=len(g))
            ax.scatter(x, g, s=22, c="black", zorder=3)
        ax.set_title(title)
        ax.set_ylabel(ylab)
        ax.axhline(0 if feature in ("IFN", "MHC") else np.nan, color="0.7", lw=0.6)

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4))
    box3(axes[0, 0], "Tacstd2", "Tacstd2 (given +2.09)", "log2 expression")
    box3(axes[0, 1], "Cldn4", "Cldn4, same contrast", "log2 expression")
    box3(axes[1, 0], "IFN", f"IFN mean-z ({len(ifn_used)}/6)", "mean z")
    box3(axes[1, 1], "MHC", f"MHC-I mean-z ({len(mhc_used)}/6)", "mean z")
    fig.suptitle("GSE239485 LLC  ·  vehicle vs Poly I:C + anti-PD-1  ·  n=8/8", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "primary_boxes.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8))
    prim = expr.loc[veh + combo]
    colors = prim["arm"].map({"vehicle": "#9aa0a6", "polyIC_aPD1": "#1a73e8"})
    for ax, axis, title in ((axes[0], "IFN", "Cldn4 vs IFN"), (axes[1], "MHC", "Cldn4 vs MHC-I")):
        ax.scatter(prim["Cldn4"], prim[axis], c=colors, s=36)
        ax.set_xlabel("Cldn4 (log2)")
        ax.set_ylabel(f"{axis} mean-z")
        ax.set_title(title + "  (primary 16)")
        hit = corr_df[(corr_df["subset"] == "primary_16") & (corr_df["gene"] == "Cldn4") & (corr_df["axis"] == axis)].iloc[0]
        ax.text(0.04, 0.96, f"ρ={hit['rho']:+.2f}  p={hit['p']:.3g}\nn=16", transform=ax.transAxes, va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "cldn4_vs_ifn_mhc.png", dpi=160)
    plt.close(fig)

    print(json.dumps(one, indent=2))
    print("IFN used", ifn_used)
    print("MHC used", mhc_used)
    print("wrote", TAB, FIG)


if __name__ == "__main__":
    main()
