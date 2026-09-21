#!/usr/bin/env python3
"""TCGA-LUAD CLDN4 vs CD8 fraction from concordant-4 signatures.

Implements PROTOCOL.md. Reads results/reference/cpm_concordant4_*.tsv and the
Xena GDC STAR TPM matrix. Writes results/tables and results/figures.
"""
from __future__ import annotations

import urllib.request
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize, nnls
from sklearn.svm import NuSVR

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
REF = ROOT / "results" / "reference"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
CLASSES = [
    "malignant", "cd8", "cd4", "nk", "b", "myeloid", "endothelial", "fibroblast",
]
HELD = {"CLDN4", "TACSTD2"}
FORCE = [
    "CD8A", "CD8B", "CD3D", "CD3E", "CD4", "FOXP3", "NKG7", "GZMB",
    "MS4A1", "CD79A", "LYZ", "CD68", "KRT8", "KRT18", "KRT19", "KRT7",
    "EPCAM", "PECAM1", "COL1A1", "COL1A2",
]
KERATIN_HOLD = {"KRT8", "KRT18", "KRT19", "KRT7", "EPCAM"}
NU_GRID = (0.25, 0.5, 0.75)


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"GET {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)
    return dest


def select_genes(S: pd.DataFrame, extra_hold: set[str] | None = None) -> list[str]:
    hold = set(HELD)
    if extra_hold:
        hold |= set(extra_hold)
    chosen: list[str] = []
    for k in CLASSES:
        if S[k].isna().all():
            continue
        nextbest = S.drop(columns=k).max(axis=1, skipna=True)
        logfc = np.log2(S[k] + 1.0) - np.log2(nextbest + 1.0)
        ok = (S[k] >= 1.0) & (logfc >= 0.585) & ~S.index.isin(hold) & S[k].notna()
        top = logfc[ok].sort_values(ascending=False).head(120)
        chosen.extend(top.index.tolist())
    for g in FORCE:
        if g in S.index and g not in hold and g not in chosen:
            chosen.append(g)
    # drop genes that are missing in any class the signature still uses
    keep = []
    for g in chosen:
        row = S.loc[g, CLASSES]
        if np.isfinite(row.to_numpy(dtype=float)).all():
            keep.append(g)
    return keep


def assert_reference(S: pd.DataFrame) -> None:
    def winner(gene: str) -> str:
        return str(S.loc[gene, CLASSES].astype(float).idxmax())

    problems = []
    if "CD8A" not in S.index or winner("CD8A") != "cd8":
        problems.append(f"CD8A peaks in {winner('CD8A') if 'CD8A' in S.index else 'MISSING'}")
    bgene = "MS4A1" if "MS4A1" in S.index else "CD79A"
    if bgene not in S.index or winner(bgene) != "b":
        problems.append(f"{bgene} peaks in {winner(bgene) if bgene in S.index else 'MISSING'}")
    egene = "EPCAM" if "EPCAM" in S.index else "KRT19"
    if egene not in S.index or winner(egene) != "malignant":
        problems.append(f"{egene} peaks in {winner(egene) if egene in S.index else 'MISSING'}")
    if problems:
        raise SystemExit("reference QC failed: " + "; ".join(problems))


def nusvr_one(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    best = None
    for nu in NU_GRID:
        model = NuSVR(kernel="linear", nu=nu, C=1.0, tol=1e-3, max_iter=5000)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.fit(X, y)
        w = np.asarray(model.coef_, dtype=float).ravel()
        w = np.clip(w, 0, None)
        if w.sum() <= 0:
            continue
        w = w / w.sum()
        recon = X @ w
        rmse = float(np.sqrt(np.mean((recon - y) ** 2)))
        if best is None or rmse < best[0]:
            if np.std(recon) == 0 or np.std(y) == 0:
                r = 0.0
            else:
                r = float(np.corrcoef(recon, y)[0, 1])
            best = (rmse, w, r)
    if best is None:
        return np.full(X.shape[1], np.nan), np.nan
    return best[1], best[2]


def nnls_one(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    w, _ = nnls(X, y)
    if w.sum() <= 0:
        return np.full(X.shape[1], np.nan), np.nan
    w = w / w.sum()
    recon = X @ w
    if np.std(recon) == 0 or np.std(y) == 0:
        r = 0.0
    else:
        r = float(np.corrcoef(recon, y)[0, 1])
    return w, r


def bayes_one(S_prop: np.ndarray, y: np.ndarray) -> np.ndarray:
    k = S_prop.shape[1]
    y = np.clip(y, 0, None)

    def loss(z: np.ndarray) -> float:
        z = z - z.max()
        theta = np.exp(z)
        theta = theta / theta.sum()
        mu = S_prop @ theta + 1e-12
        return float(-np.dot(y, np.log(mu)))

    res = minimize(loss, np.zeros(k), method="L-BFGS-B")
    z = res.x - res.x.max()
    theta = np.exp(z)
    theta = theta / theta.sum()
    return theta


def run_methods(S: pd.DataFrame, genes: list[str], bulk: pd.DataFrame) -> dict[str, pd.DataFrame]:
    X = S.loc[genes, CLASSES].to_numpy(dtype=float)
    # columns -> proportions for the Bayes MAP; CPM for ν-SVR and NNLS
    colsum = X.sum(axis=0)
    Sprop = X / colsum
    Y = bulk.loc[genes].to_numpy(dtype=float)
    n = Y.shape[1]
    k = len(CLASSES)
    out = {m: np.zeros((n, k)) for m in ("nusvr", "nnls", "bayes")}
    fit_r = {m: np.zeros(n) for m in ("nusvr", "nnls")}
    for i in range(n):
        y = Y[:, i]
        w, r = nusvr_one(X, y)
        out["nusvr"][i] = w
        fit_r["nusvr"][i] = r
        w, r = nnls_one(X, y)
        out["nnls"][i] = w
        fit_r["nnls"][i] = r
        out["bayes"][i] = bayes_one(Sprop, y)
        if (i + 1) % 100 == 0:
            print(f"  deconv {i+1}/{n}", flush=True)
    frames = {}
    for m in out:
        frames[m] = pd.DataFrame(out[m], index=bulk.columns, columns=CLASSES)
    fit = pd.DataFrame(fit_r, index=bulk.columns)
    return {"fractions": frames, "fit": fit}


def spearman(x: np.ndarray, y: np.ndarray) -> dict:
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 10:
        return {"n": n, "rho": np.nan, "p": np.nan, "lo": np.nan, "hi": np.nan, "k": 0}
    rho, p = stats.spearmanr(x[m], y[m])
    z = np.arctanh(np.clip(rho, -0.999, 0.999))
    se = 1 / np.sqrt(n - 3)
    return {
        "n": n, "rho": float(rho), "p": float(p), "k": 0,
        "lo": float(np.tanh(z - 1.96 * se)),
        "hi": float(np.tanh(z + 1.96 * se)),
    }


def partial_spearman(x: np.ndarray, y: np.ndarray, cov: np.ndarray) -> dict:
    if cov.ndim == 1:
        cov = cov[:, None]
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(cov).all(axis=1)
    n = int(m.sum())
    k = cov.shape[1]
    if n < k + 10:
        return {"n": n, "rho": np.nan, "p": np.nan, "lo": np.nan, "hi": np.nan, "k": k}
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    C = np.column_stack([stats.rankdata(cov[m, j]) for j in range(k)])
    A = np.column_stack([np.ones(n), C])

    def resid(v: np.ndarray) -> np.ndarray:
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta

    rx, ry = resid(xr), resid(yr)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan, "lo": np.nan, "hi": np.nan, "k": k}
    rho = float(np.corrcoef(rx, ry)[0, 1])
    df = n - 2 - k
    if abs(rho) >= 1 or df <= 0:
        p = 0.0 if abs(rho) >= 1 else np.nan
    else:
        t = rho * np.sqrt(df / (1 - rho ** 2))
        p = float(2 * stats.t.sf(abs(t), df))
    z = np.arctanh(np.clip(rho, -0.999, 0.999))
    se = 1 / np.sqrt(max(n - 3 - k, 1))
    return {
        "n": n, "rho": rho, "p": p, "k": k,
        "lo": float(np.tanh(z - 1.96 * se)),
        "hi": float(np.tanh(z + 1.96 * se)),
    }


def bh(pvals: list[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(np.where(np.isfinite(p), p, 1.0))
    ranked = np.where(np.isfinite(p[order]), p[order], 1.0)
    q = ranked * m / (np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(m)
    out[order] = q
    return out


def synthetic_check() -> pd.DataFrame:
    rng = np.random.default_rng(472)
    k = 4
    genes = 60
    S = rng.uniform(1, 5, size=(genes, k))
    for i in range(k):
        S[i * 8:(i + 1) * 8, :] = 2
        S[i * 8:(i + 1) * 8, i] = 400
    theta = rng.dirichlet(np.ones(k) * 0.8, size=25)
    Y = S @ theta.T
    Y = np.clip(Y * (1 + rng.normal(0, 0.08, Y.shape)), 0, None)
    rows = []
    est = {"nusvr": [], "nnls": [], "bayes": []}
    Sprop = S / S.sum(axis=0)
    for i in range(Y.shape[1]):
        w, _ = nusvr_one(S, Y[:, i])
        est["nusvr"].append(w)
        w, _ = nnls_one(S, Y[:, i])
        est["nnls"].append(w)
        est["bayes"].append(bayes_one(Sprop, Y[:, i]))
    for m, arr in est.items():
        E = np.vstack(arr)
        rho = float(np.corrcoef(E.ravel(), theta.ravel())[0, 1])
        rows.append({"method": m, "pearson_theta": rho, "n_mixtures": theta.shape[0]})
        if not np.isfinite(rho) or rho < 0.9:
            raise SystemExit(f"synthetic recovery failed for {m}: r={rho}")
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "synthetic_recovery.tsv", sep="\t", index=False)
    print(df.to_string(index=False), flush=True)
    return df


def load_bulk() -> pd.DataFrame:
    path = _download(
        "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz",
        RAW / "TCGA-LUAD.star_tpm.tsv.gz",
    )
    feat = _download(
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5699nnn/GSM5699777/suppl/GSM5699777_TD1_features.tsv.gz",
        RAW / "GSM5699777_TD1_features.tsv.gz",
    )
    features = pd.read_csv(feat, sep="\t", header=None, names=["ensembl", "symbol", "kind"])
    features["ensembl"] = features["ensembl"].astype(str)
    features["symbol"] = features["symbol"].astype(str).str.upper()
    mp = features.drop_duplicates("ensembl").set_index("ensembl")["symbol"]
    print("reading TCGA-LUAD STAR TPM", flush=True)
    expr = pd.read_csv(path, sep="\t", index_col=0)
    ens = pd.Series(expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True), index=expr.index)
    symbol = ens.map(mp)
    keep = symbol.notna() & (symbol != "") & (symbol.str.upper() != "NAN")
    expr = expr.loc[keep]
    symbol = symbol.loc[keep]
    linear = np.clip(np.exp2(expr.to_numpy(dtype=float)) - 1.0, 0, None)
    df = pd.DataFrame(linear, index=symbol.to_numpy(), columns=expr.columns)
    df = df.groupby(level=0).sum()
    # one primary tumor per patient, prefer vial A
    chosen = {}
    for s in df.columns:
        parts = str(s).split("-")
        if len(parts) < 4 or parts[3][:2] != "01":
            continue
        patient = "-".join(parts[:3])
        vial = parts[3]
        prev = chosen.get(patient)
        if prev is None or vial < prev[0]:
            chosen[patient] = (vial, s)
    cols = [v[1] for v in chosen.values()]
    df = df[cols]
    print(f"TCGA-LUAD primary tumors, one per patient: {df.shape[1]}", flush=True)
    return df


def outlier_filter(genes: list[str], S: pd.DataFrame, bulk: pd.DataFrame) -> list[str]:
    keep = []
    dropped = []
    for g in genes:
        if g.startswith(("MT-", "RPL", "RPS", "HBA", "HBB", "HBD")):
            dropped.append(g)
            continue
        if g not in bulk.index:
            dropped.append(g)
            continue
        mx_b = float(np.nanmax(bulk.loc[g].to_numpy(dtype=float)))
        mx_s = float(np.nanmax(S.loc[g, CLASSES].to_numpy(dtype=float)))
        if mx_s <= 0 or mx_b > 50.0 * mx_s:
            dropped.append(g)
            continue
        keep.append(g)
    print(f"signature genes kept {len(keep)} dropped {len(dropped)}", flush=True)
    return keep


def add_row(rows: list, family: str, **kwargs) -> None:
    rows.append({"family": family, **kwargs})


def associations(expr: pd.DataFrame, fractions: dict[str, pd.DataFrame], tag: str) -> list[dict]:
    # expr here is log2(tpm+1) aligned to fraction index
    rows = []
    cldn = expr["CLDN4"].to_numpy()
    cd8a = expr["CD8A"].to_numpy()
    actb = expr["ACTB"].to_numpy()
    krt = expr[["KRT8", "KRT18", "KRT19"]].to_numpy()
    krt18_19 = expr[["KRT18", "KRT19"]].to_numpy()
    for method, frac in fractions.items():
        cd8 = frac["cd8"].to_numpy()
        mal = frac["malignant"].to_numpy()
        cd8nk = (frac["cd8"] + frac["nk"]).to_numpy()
        if method == "nusvr":
            add_row(rows, "biological", tag=tag, method=method, contrast="CLDN4_vs_CD8_fraction",
                    adjustment="none", **spearman(cldn, cd8))
            add_row(rows, "biological", tag=tag, method=method, contrast="CLDN4_vs_CD8_fraction",
                    adjustment="malignant_fraction", **partial_spearman(cldn, cd8, mal))
            add_row(rows, "biological", tag=tag, method=method, contrast="CLDN4_vs_CD8_fraction",
                    adjustment="KRT8_KRT18_KRT19", **partial_spearman(cldn, cd8, krt))
            add_row(rows, "biological", tag=tag, method=method, contrast="CLDN4_vs_CD8_plus_NK",
                    adjustment="none", **spearman(cldn, cd8nk))
        if method in ("nusvr", "bayes", "nnls"):
            add_row(rows, "control", tag=tag, method=method, contrast="CD8A_vs_CD8_fraction",
                    adjustment="none", **spearman(cd8a, cd8))
            add_row(rows, "control", tag=tag, method=method, contrast="ACTB_vs_CD8_fraction",
                    adjustment="none", **spearman(actb, cd8))
            add_row(rows, "control", tag=tag, method=method, contrast="KRT8_vs_malignant_fraction",
                    adjustment="none", **spearman(expr["KRT8"].to_numpy(), mal))
        if method == "bayes":
            add_row(rows, "biological", tag=tag, method=method, contrast="CLDN4_vs_CD8_fraction",
                    adjustment="none", **spearman(cldn, cd8))
            add_row(rows, "biological", tag=tag, method=method, contrast="CLDN4_vs_CD8_fraction",
                    adjustment="malignant_fraction", **partial_spearman(cldn, cd8, mal))
            add_row(rows, "biological", tag=tag, method=method, contrast="CLDN4_vs_CD8_fraction",
                    adjustment="KRT8_KRT18_KRT19", **partial_spearman(cldn, cd8, krt))
        elif method == "nnls":
            add_row(rows, "biological", tag=tag, method=method, contrast="CLDN4_vs_CD8_fraction",
                    adjustment="none", **spearman(cldn, cd8))
        if method == "nusvr":
            add_row(rows, "control", tag=tag, method="nusvr_vs_bayes", contrast="CD8_fraction_concordance",
                    adjustment="none", **spearman(cd8, fractions["bayes"]["cd8"].to_numpy()))
    if tag == "primary":
        add_row(rows, "biological", tag="bulk", method="bulk", contrast="CLDN4_vs_CD8A",
                adjustment="none", **spearman(cldn, cd8a))
        add_row(rows, "biological", tag="bulk", method="bulk", contrast="CLDN4_vs_CD8A",
                adjustment="KRT8_KRT18_KRT19", **partial_spearman(cldn, cd8a, krt))
        add_row(rows, "biological", tag="bulk", method="bulk", contrast="CLDN4_vs_KRT8",
                adjustment="KRT18_KRT19", **partial_spearman(cldn, expr["KRT8"].to_numpy(), krt18_19))
    return rows


def forest_plot(df: pd.DataFrame, path: Path) -> None:
    show = df.loc[df["family"] == "biological"].copy()
    # primary-signature rows plus bulk; drop the keratin-hold and author tags from the main figure
    show = show.loc[show["tag"].isin(["primary", "bulk"])]
    show["label"] = show["method"] + " | " + show["contrast"].str.replace("_", " ") + " | " + show["adjustment"]
    show = show.iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, 7.2))
    y = np.arange(len(show))
    ax.axvline(0, color="#444444", lw=0.8)
    ax.errorbar(
        show["rho"], y,
        xerr=[show["rho"] - show["lo"], show["hi"] - show["rho"]],
        fmt="o", color="#1f4e79", ecolor="#5b7c99", capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(show["label"], fontsize=8)
    ax.set_xlabel("Spearman ρ (Fisher z 95% interval)")
    ax.set_title("TCGA-LUAD CLDN4 vs CD8 — bulk and concordant-4 deconvolution")
    fig.tight_layout()
    fig.subplots_adjust(left=0.48)
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def scatter(x, y, xlabel, ylabel, title, path: Path) -> None:
    m = np.isfinite(x) & np.isfinite(y)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.scatter(x[m], y[m], s=12, alpha=0.55, c="#1f4e79", linewidths=0)
    rho, p = stats.spearmanr(x[m], y[m])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\nρ={rho:.3f}, p={p:.2e}, n={m.sum()}")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    synthetic_check()
    primary = pd.read_csv(REF / "cpm_concordant4_primary.tsv", sep="\t", index_col=0)
    author = pd.read_csv(REF / "cpm_concordant4_author_malignant.tsv", sep="\t", index_col=0)
    assert_reference(primary)
    assert_reference(author)
    bulk_lin = load_bulk()
    need = ["CLDN4", "CD8A", "ACTB", "KRT8", "KRT18", "KRT19"]
    missing = [g for g in need if g not in bulk_lin.index]
    if missing:
        raise SystemExit(f"bulk missing {missing}")
    log = np.log2(bulk_lin.loc[need] + 1.0).T
    ok = np.isfinite(log.to_numpy()).all(axis=1)
    log = log.loc[ok]
    bulk_lin = bulk_lin[log.index]
    print(f"patients with CLDN4/CD8A/ACTB/keratins: {bulk_lin.shape[1]}", flush=True)

    specs = {
        "primary": (primary, select_genes(primary)),
        "keratin_held_out": (primary, select_genes(primary, KERATIN_HOLD)),
        "author_malignant": (author, select_genes(author)),
    }
    all_rows = []
    saved_primary = None
    for tag, (S, genes0) in specs.items():
        genes = outlier_filter(genes0, S, bulk_lin)
        if len(genes) < 50:
            raise SystemExit(f"{tag} signature has only {len(genes)} genes")
        pd.Series(genes, name="gene").to_csv(TAB / f"signature_genes_{tag}.tsv", sep="\t", index=False)
        print(f"running {tag} on {len(genes)} genes", flush=True)
        result = run_methods(S, genes, bulk_lin)
        result["fit"].to_csv(TAB / f"fit_{tag}.tsv", sep="\t")
        print(tag, "median reconstruction r", result["fit"].median().to_dict(), flush=True)
        expr = np.log2(bulk_lin.loc[need] + 1.0).T.loc[result["fractions"]["nusvr"].index]
        all_rows.extend(associations(expr, result["fractions"], tag))
        if tag == "primary":
            saved_primary = (expr, result["fractions"])
            wide = expr.copy()
            wide.index.name = "sample"
            for method, frac in result["fractions"].items():
                ren = frac.add_prefix(method + "_")
                wide = wide.join(ren)
            wide.to_csv(TAB / "luad_sample_level.tsv", sep="\t")
    res = pd.DataFrame(all_rows)
    # BH within biological tests only
    bio = res["family"] == "biological"
    res.loc[bio, "q"] = bh(res.loc[bio, "p"].tolist())
    res.loc[~bio, "q"] = np.nan
    res.to_csv(TAB / "associations.tsv", sep="\t", index=False)
    print(res.loc[res["tag"].isin(["primary", "bulk"]),
                  ["tag", "method", "contrast", "adjustment", "n", "rho", "p", "q"]].to_string(index=False),
          flush=True)
    forest_plot(res, FIG / "cldn4_cd8_forest.png")
    expr, fr = saved_primary
    scatter(
        expr["CLDN4"].to_numpy(), fr["nusvr"]["cd8"].to_numpy(),
        "CLDN4 log2(TPM+1)", "ν-SVR CD8 fraction",
        "TCGA-LUAD primary", FIG / "scatter_cldn4_nusvr_cd8.png",
    )
    scatter(
        expr["CD8A"].to_numpy(), fr["nusvr"]["cd8"].to_numpy(),
        "CD8A log2(TPM+1)", "ν-SVR CD8 fraction",
        "ν-SVR control", FIG / "scatter_cd8a_nusvr_cd8.png",
    )
    scatter(
        expr["CD8A"].to_numpy(), fr["bayes"]["cd8"].to_numpy(),
        "CD8A log2(TPM+1)", "Bayes MAP CD8 fraction",
        "Bayes MAP control", FIG / "scatter_cd8a_bayes_cd8.png",
    )
    scatter(
        expr["CLDN4"].to_numpy(), fr["bayes"]["cd8"].to_numpy(),
        "CLDN4 log2(TPM+1)", "Bayes MAP CD8 fraction",
        "TCGA-LUAD primary", FIG / "scatter_cldn4_bayes_cd8.png",
    )
    print("wrote", TAB, FIG, flush=True)


if __name__ == "__main__":
    main()
