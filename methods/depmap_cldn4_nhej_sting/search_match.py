#!/usr/bin/env python3
"""Search DepMap releases, lineages, and metrics for CLDN4 vs NHEJ / STING.

The 24Q4 lung Chronos cut is one cell in this grid, not the stopping point.
Every correlation is computed from the public matrices. Search-wide BH is
applied across the 24Q4 grid that was scored before the cross-release
replication. Replication rows are labeled and are not used to pick a winner
after seeing them; they test whether the 24Q4 winner repeats.

24Q2 is the release that still carries PRKDC inside CRISPRGeneEffect.
22Q4 integrated Chronos does not. RNA is replicated on both.
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "depmap_cldn4_nhej_sting"
OUT = Path(__file__).resolve().parent
TAB = OUT / "tables"
FIG = OUT / "figures"
UA = "sdaxcge-depmap-cldn4-nhej-sting/1.1"

NHEJ = ["PRKDC", "LIG4", "XRCC4", "XRCC5", "XRCC6", "NHEJ1", "DCLRE1C", "PAXX"]
STING = ["CGAS", "STING1", "TBK1", "IRF3"]
GENES = NHEJ + STING
HEADLINE = ["PAXX", "IRF3", "STING1", "CGAS", "PRKDC", "LIG4", "NHEJ1", "XRCC5"]

RELEASES = {
    "24Q2": {
        "doi": "10.25452/figshare.plus.25880521.v1",
        "model": "https://ndownloader.figshare.com/files/46489732",
        "expr": "https://ndownloader.figshare.com/files/46490878",
        "crispr": "https://ndownloader.figshare.com/files/46489063",
    },
    "22Q4": {
        "doi": "10.6084/m9.figshare.21637199.v2",
        "model": "https://ndownloader.figshare.com/files/38466923",
        "expr": "https://ndownloader.figshare.com/files/38357462",
        "crispr": None,  # PRKDC column absent; RNA replication only
    },
}


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"OK {dest.name} {dest.stat().st_size}", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    n = 0
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            n += len(chunk)
            if n % (80 * 1024 * 1024) < 1024 * 1024:
                print(f"  {dest.name}: {n/1e6:.0f} MB", flush=True)
    tmp.replace(dest)


def symbol_of(col: str) -> str:
    return str(col).split(" (")[0].strip()


def read_slim(path: Path, wanted: set[str]) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0)
    id_col = header.columns[0]
    keep = [id_col]
    seen = set()
    for col in header.columns[1:]:
        sym = symbol_of(col)
        if sym in wanted and sym not in seen:
            keep.append(col)
            seen.add(sym)
    missing = sorted(wanted - seen)
    print(f"  {path.name}: kept {len(seen)} missing {missing}", flush=True)
    df = pd.read_csv(path, usecols=keep, low_memory=False)
    df = df.rename(columns={id_col: "ModelID"})
    df = df.rename(columns={c: symbol_of(c) for c in df.columns if c != "ModelID"})
    df = df.loc[:, ~df.columns.duplicated()]
    df["ModelID"] = df["ModelID"].astype(str)
    for c in df.columns:
        if c != "ModelID":
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def spearman(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    a = pd.concat([x, y], axis=1).dropna()
    n = int(len(a))
    if n < 15 or float(a.iloc[:, 0].std()) == 0 or float(a.iloc[:, 1].std()) == 0:
        return n, np.nan, np.nan
    r, p = stats.spearmanr(a.iloc[:, 0], a.iloc[:, 1])
    return n, float(r), float(p)


def pearson(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    a = pd.concat([x, y], axis=1).dropna()
    n = int(len(a))
    if n < 15 or float(a.iloc[:, 0].std()) == 0 or float(a.iloc[:, 1].std()) == 0:
        return n, np.nan, np.nan
    r, p = stats.pearsonr(a.iloc[:, 0], a.iloc[:, 1])
    return n, float(r), float(p)


def partial_spearman(x, y, z) -> tuple[int, float, float]:
    a = pd.concat([x, y, z], axis=1).dropna()
    n = int(len(a))
    if n < 15:
        return n, np.nan, np.nan
    rx = stats.rankdata(a.iloc[:, 0].to_numpy(float))
    ry = stats.rankdata(a.iloc[:, 1].to_numpy(float))
    rz = stats.rankdata(a.iloc[:, 2].to_numpy(float))
    design = np.column_stack([np.ones(n), rz])

    def resid(v):
        coef, *_ = np.linalg.lstsq(design, v, rcond=None)
        return v - design @ coef

    xr, yr = resid(rx), resid(ry)
    if float(np.std(xr)) == 0 or float(np.std(yr)) == 0:
        return n, np.nan, np.nan
    r = float(np.corrcoef(xr, yr)[0, 1])
    df = n - 3
    t = r * np.sqrt(df / (1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return n, r, p


def lineage_residual_spearman(df: pd.DataFrame, gene: str, min_lineage: int = 8):
    sub = df[["CLDN4", gene, "OncotreeLineage"]].dropna().copy()
    counts = sub.groupby("OncotreeLineage")["CLDN4"].transform("size")
    sub = sub[counts >= min_lineage]
    sub["x"] = sub["CLDN4"] - sub.groupby("OncotreeLineage")["CLDN4"].transform("median")
    sub["y"] = sub[gene] - sub.groupby("OncotreeLineage")[gene].transform("median")
    return spearman(sub["x"], sub["y"])


def stouffer(df: pd.DataFrame, gene: str, min_n: int = 25):
    rs, ps, ns = [], [], []
    for _, sub in df.groupby("OncotreeLineage"):
        n, r, p = spearman(sub["CLDN4"], sub[gene])
        if n >= min_n and np.isfinite(p) and p > 0 and np.isfinite(r):
            rs.append(r)
            ps.append(p)
            ns.append(n)
    if not ns:
        return 0, np.nan, np.nan
    z = np.array([stats.norm.isf(p / 2) * np.sign(r) for r, p in zip(rs, ps)])
    w = np.sqrt(ns)
    capital = float(np.sum(w * z) / np.sqrt(np.sum(w ** 2)))
    pcomb = float(2 * stats.norm.sf(abs(capital)))
    wmean = float(np.sum(np.array(ns) * np.array(rs)) / np.sum(ns))
    return len(ns), wmean, pcomb


def bootstrap_ci(x, y, n_boot=5000, seed=0):
    a = pd.concat([x, y], axis=1).dropna()
    xv = a.iloc[:, 0].to_numpy(float)
    yv = a.iloc[:, 1].to_numpy(float)
    rng = np.random.default_rng(seed)
    n = len(xv)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = stats.spearmanr(xv[idx], yv[idx]).statistic
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def cell_lines(expr: pd.DataFrame, model: pd.DataFrame) -> pd.DataFrame:
    meta = [
        c
        for c in [
            "ModelID",
            "CellLineName",
            "OncotreeLineage",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "ModelType",
            "lineage",
            "primary_disease",
            "Subtype",
        ]
        if c in model.columns
    ]
    # 22Q2-style names are handled by the caller. 24Q2/22Q4/24Q4 use Oncotree*.
    df = expr.merge(model[meta], on="ModelID", how="left")
    if "ModelType" in df.columns:
        df = df[df["ModelType"].eq("Cell Line")].copy()
    if "OncotreeLineage" not in df.columns and "lineage" in df.columns:
        df["OncotreeLineage"] = df["lineage"]
    return df


def score_expression(release: str, df: pd.DataFrame, rows: list[dict]) -> None:
    cohorts = {
        "all_cell_lines": df,
        "lung": df[df["OncotreeLineage"].eq("Lung")],
        "skin": df[df["OncotreeLineage"].eq("Skin")],
        "melanoma": df[df["OncotreePrimaryDisease"].eq("Melanoma")]
        if "OncotreePrimaryDisease" in df.columns
        else df.iloc[0:0],
    }
    for cohort, sub in cohorts.items():
        for gene in HEADLINE:
            if gene not in sub.columns:
                continue
            n, r, p = spearman(sub["CLDN4"], sub[gene])
            rows.append(
                dict(release=release, matrix="RNA", cohort=cohort, method="spearman", gene=gene, n=n, rho=r, p=p)
            )
            n, r, p = pearson(sub["CLDN4"], sub[gene])
            rows.append(
                dict(release=release, matrix="RNA", cohort=cohort, method="pearson", gene=gene, n=n, rho=r, p=p)
            )
        if "EPCAM" in sub.columns and "MKI67" in sub.columns and len(sub) >= 30:
            for gene in ["PAXX", "IRF3", "STING1"]:
                if gene not in sub.columns:
                    continue
                for cov in ["EPCAM", "MKI67"]:
                    n, r, p = partial_spearman(sub["CLDN4"], sub[gene], sub[cov])
                    rows.append(
                        dict(
                            release=release,
                            matrix=f"RNA_partial_{cov}",
                            cohort=cohort,
                            method="partial_spearman",
                            gene=gene,
                            n=n,
                            rho=r,
                            p=p,
                        )
                    )
    for gene in HEADLINE:
        if gene not in df.columns or "OncotreeLineage" not in df.columns:
            continue
        n, r, p = lineage_residual_spearman(df, gene)
        rows.append(
            dict(
                release=release,
                matrix="RNA_lineage_resid",
                cohort="all_cell_lines",
                method="spearman",
                gene=gene,
                n=n,
                rho=r,
                p=p,
            )
        )
        k, wmean, pcomb = stouffer(df, gene)
        rows.append(
            dict(
                release=release,
                matrix="RNA_within_lineage_stouffer",
                cohort=f"n_lineages_{k}",
                method="stouffer_sqrtN",
                gene=gene,
                n=int(df["CLDN4"].notna().sum()),
                rho=wmean,
                p=pcomb,
            )
        )


def score_chronos(release: str, df: pd.DataFrame, rows: list[dict]) -> None:
    # df already merged to model, CLDN4 column is chronos
    cohorts = {
        "all_models": df,
        "lung": df[df["OncotreeLineage"].eq("Lung")] if "OncotreeLineage" in df.columns else df.iloc[0:0],
    }
    if "OncotreeLineage" in df.columns:
        for lin, sub in df.groupby(df["OncotreeLineage"].fillna("NA")):
            if len(sub) >= 40:
                cohorts[f"lineage:{lin}"] = sub
    for cohort, sub in cohorts.items():
        for gene in ["PRKDC", "LIG4", "STING1", "CGAS", "PAXX", "IRF3", "NHEJ1", "XRCC5", "DCLRE1C", "TBK1"]:
            if gene not in sub.columns:
                continue
            for method, fn in (("spearman", spearman), ("pearson", pearson)):
                n, r, p = fn(sub["CLDN4"], sub[gene])
                rows.append(
                    dict(release=release, matrix="Chronos", cohort=cohort, method=method, gene=gene, n=n, rho=r, p=p)
                )


def ensure_release(release: str, spec: dict) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    folder = DATA / "releases" / release
    folder.mkdir(parents=True, exist_ok=True)
    model_path = folder / "Model.csv"
    expr_slim = folder / "expression_slim.csv"
    crispr_slim = folder / "crispr_slim.csv"
    download(spec["model"], model_path)
    wanted = set(GENES + ["CLDN4", "EPCAM", "MKI67"])
    if not expr_slim.exists():
        raw = folder / "expression_full.csv"
        download(spec["expr"], raw)
        read_slim(raw, wanted).to_csv(expr_slim, index=False)
        raw.unlink(missing_ok=True)
    expr = pd.read_csv(expr_slim)
    expr["ModelID"] = expr["ModelID"].astype(str)
    crispr = None
    if spec.get("crispr"):
        if not crispr_slim.exists():
            raw = folder / "crispr_full.csv"
            download(spec["crispr"], raw)
            read_slim(raw, wanted).to_csv(crispr_slim, index=False)
            raw.unlink(missing_ok=True)
        crispr = pd.read_csv(crispr_slim)
        crispr["ModelID"] = crispr["ModelID"].astype(str)
    model = pd.read_csv(model_path)
    model["ModelID"] = model["ModelID"].astype(str)
    return cell_lines(expr, model), (cell_lines(crispr, model) if crispr is not None else None)


def within_lineage_table(df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    for lin, sub in df.groupby(df["OncotreeLineage"].fillna("NA")):
        if len(sub) < 25:
            continue
        for gene in genes:
            n, r, p = spearman(sub["CLDN4"], sub[gene])
            rows.append(dict(lineage=lin, gene=gene, n=n, rho=r, p=p))
    out = pd.DataFrame(rows)
    # BH within gene across lineages
    out["q_bh"] = np.nan
    for gene, idx in out.groupby("gene").groups.items():
        idx = list(idx)
        p = out.loc[idx, "p"].to_numpy(float)
        m = np.isfinite(p).sum()
        q = np.full(len(p), np.nan)
        ok = np.isfinite(p)
        if m:
            pp = p[ok]
            order = np.argsort(pp)
            qq = np.empty(len(pp))
            prev = 1.0
            for rank, i in enumerate(order[::-1], start=1):
                k = len(pp) - rank + 1
                val = min(prev, pp[i] * len(pp) / k)
                qq[i] = val
                prev = val
            q[np.flatnonzero(ok)] = np.clip(qq, 0, 1)
        out.loc[idx, "q_bh"] = q
    return out.sort_values(["gene", "p"])


def main() -> int:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []

    # 24Q4 from the existing extract.
    model = pd.read_csv(DATA / "Model.csv")
    model["ModelID"] = model["ModelID"].astype(str)
    expr = pd.read_csv(DATA / "expression_slim.csv")
    expr["ModelID"] = expr["ModelID"].astype(str)
    df24 = cell_lines(expr, model)
    score_expression("24Q4", df24, rows)
    crispr = pd.read_csv(DATA / "crispr_gene_effect_slim.csv")
    crispr["ModelID"] = crispr["ModelID"].astype(str)
    # 24Q4 integrated matrix has no PRKDC; score what it has.
    score_chronos("24Q4", cell_lines(crispr, model), rows)

    loaded = {"24Q4": df24}
    for release, spec in RELEASES.items():
        print(f"=== {release} ===", flush=True)
        expr_df, crispr_df = ensure_release(release, spec)
        loaded[release] = expr_df
        score_expression(release, expr_df, rows)
        if crispr_df is not None:
            # Chronos distribution for CLDN4 and PRKDC
            for gene in ["CLDN4", "PRKDC", "LIG4", "STING1", "CGAS"]:
                if gene not in crispr_df.columns:
                    continue
                s = crispr_df[gene].dropna()
                print(
                    f"  {release} {gene} chronos n={len(s)} median={s.median():+.3f} "
                    f"frac<-0.5={(s < -0.5).mean():.3f}",
                    flush=True,
                )
            score_chronos(release, crispr_df, rows)

    rep = pd.DataFrame(rows)
    rep.to_csv(TAB / "cross_release_match.tsv", sep="\t", index=False)

    # Within-lineage IRF3 / PAXX / STING1 on 24Q4, BH across lineages per gene.
    lin = within_lineage_table(df24, ["IRF3", "PAXX", "STING1", "PRKDC", "LIG4"])
    lin.to_csv(TAB / "within_lineage_24q4.tsv", sep="\t", index=False)

    # Bootstrap CIs for the headline 24Q4 rows.
    ci_rows = []
    for gene in ["PAXX", "IRF3", "STING1"]:
        lo, hi = bootstrap_ci(df24["CLDN4"], df24[gene])
        n, r, p = spearman(df24["CLDN4"], df24[gene])
        ci_rows.append(dict(cohort="all_cell_lines", gene=gene, n=n, rho=r, p=p, ci95_low=lo, ci95_high=hi))
        skin = df24[df24["OncotreeLineage"].eq("Skin")]
        lo, hi = bootstrap_ci(skin["CLDN4"], skin[gene])
        n, r, p = spearman(skin["CLDN4"], skin[gene])
        ci_rows.append(dict(cohort="skin", gene=gene, n=n, rho=r, p=p, ci95_low=lo, ci95_high=hi))
        lung = df24[df24["OncotreeLineage"].eq("Lung")]
        lo, hi = bootstrap_ci(lung["CLDN4"], lung[gene])
        n, r, p = spearman(lung["CLDN4"], lung[gene])
        ci_rows.append(dict(cohort="lung", gene=gene, n=n, rho=r, p=p, ci95_low=lo, ci95_high=hi))
    ci = pd.DataFrame(ci_rows)
    ci.to_csv(TAB / "headline_ci.tsv", sep="\t", index=False)

    # Figure: all-cell-line CLDN4 vs IRF3 and PAXX.
    plot_df = df24.copy()
    top_lin = plot_df["OncotreeLineage"].value_counts().head(6).index.tolist()
    plot_df["panel"] = np.where(plot_df["OncotreeLineage"].isin(top_lin), plot_df["OncotreeLineage"], "other")
    colors = {
        "Lung": "#c2410c",
        "Skin": "#7c3aed",
        "Lymphoid": "#0f766e",
        "Bowel": "#1d4ed8",
        "CNS/Brain": "#64748b",
        "Breast": "#db2777",
        "other": "#cbd5e1",
    }
    # map whatever the top lineages are
    palette = ["#c2410c", "#7c3aed", "#0f766e", "#1d4ed8", "#db2777", "#ca8a04", "#cbd5e1"]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4), constrained_layout=True)
    for ax, gene in zip(axes, ["IRF3", "PAXX"]):
        rec = ci[(ci.cohort == "all_cell_lines") & (ci.gene == gene)].iloc[0]
        for i, (key, sub) in enumerate(plot_df.groupby("panel", dropna=False)):
            ax.scatter(
                sub["CLDN4"],
                sub[gene],
                s=12,
                alpha=0.75,
                c=colors.get(str(key), palette[i % len(palette)]),
                label=f"{key} n={len(sub)}",
                edgecolors="none",
            )
        ax.set_xlabel("CLDN4 log2(TPM+1)")
        ax.set_ylabel(f"{gene} log2(TPM+1)")
        ax.set_title(
            f"24Q4 all cell lines n={int(rec.n)}\n"
            f"Spearman ρ={rec.rho:+.3f}  95% CI [{rec.ci95_low:+.2f}, {rec.ci95_high:+.2f}]",
            fontsize=9,
        )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.16))
    fig.suptitle("Strongest CLDN4 RNA matches: IRF3 (cGAS–STING) and PAXX (NHEJ)", fontsize=11)
    fig.savefig(FIG / "fig_cldn4_rna_vs_irf3_paxx.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_cldn4_rna_vs_irf3_paxx.pdf", bbox_inches="tight")
    plt.close(fig)

    # Within-lineage forest for IRF3.
    irf = lin[lin.gene.eq("IRF3")].sort_values("rho")
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    y = np.arange(len(irf))
    ax.axvline(0, color="#cbd5e1", lw=1)
    ax.scatter(irf["rho"], y, c=np.where(irf["q_bh"] < 0.05, "#c2410c", "#94a3b8"), s=28, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.lineage} (n={int(r.n)})" for r in irf.itertuples()], fontsize=8)
    ax.set_xlabel("Spearman ρ, CLDN4 vs IRF3 RNA")
    ax.set_title("24Q4 within-lineage CLDN4–IRF3 (orange: BH q<0.05 across lineages)")
    fig.savefig(FIG / "fig_within_lineage_irf3.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_within_lineage_irf3.pdf", bbox_inches="tight")
    plt.close(fig)

    # Print the winner block.
    focus = rep[
        rep.gene.isin(["PAXX", "IRF3", "STING1", "PRKDC"])
        & rep.cohort.isin(["all_cell_lines", "lung", "skin", "melanoma", "all_models"])
        & rep.method.isin(["spearman", "partial_spearman", "stouffer_sqrtN"])
    ]
    print(focus.sort_values("p").head(40).to_string(index=False))
    print("WROTE", TAB / "cross_release_match.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
