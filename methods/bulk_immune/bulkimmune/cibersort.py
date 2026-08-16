"""CIBERSORT-style support vector regression deconvolution.

This is a reimplementation of the *classic* CIBERSORT core algorithm (Newman et
al., Nature Methods 2015) in relative mode: nu-SVR with a linear kernel over
nu in {0.25, 0.5, 0.75}, best model chosen by RMSE, negative weights clipped and
the remainder rescaled to sum to one, with a permutation p-value for the null
"this mixture contains none of these cell types".

**What this is not.** CIBERSORTx (Newman et al., Nature Biotechnology 2019) adds
B-mode/S-mode batch correction, absolute mode, and signature-matrix inference
from single-cell data. Those are only available from the authors' web server or
licensed Docker image and are *not* reproduced here; ``scripts/R/run_cibersortx.sh``
shows how to drive the official container so the numbers in a paper come from
the tool the paper cites.

**LM22 is not distributed with this repository.** It is available from the
CIBERSORTx portal (https://cibersortx.stanford.edu) after registering and
accepting the Stanford academic licence. Point ``--signature`` at the file you
downloaded. For a licence-free smoke test the pipeline falls back to the TIL10
matrix from quanTIseq (Finotello et al., Genome Medicine 2019, GPL-licensed),
which the fetch script downloads.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd

__all__ = ["load_signature_matrix", "cibersort", "constrained_ls_deconvolve", "CibersortResult"]

NU_GRID = (0.25, 0.5, 0.75)


def load_signature_matrix(path: str) -> pd.DataFrame:
    """Read a tab-separated signature matrix (genes x cell types)."""
    sig = pd.read_csv(path, sep="\t", index_col=0)
    sig = sig.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    if sig.index.has_duplicates:
        sig = sig.groupby(level=0).mean()
    return sig


@dataclass
class CibersortResult:
    fractions: pd.DataFrame       # samples x cell types, rows sum to 1
    correlation: pd.Series        # Pearson r between reconstructed and observed
    rmse: pd.Series
    p_value: pd.Series            # permutation p, NaN if permutations == 0
    n_genes_used: int


def _core(X: np.ndarray, y: np.ndarray, random_state: int = 0) -> tuple[np.ndarray, float, float]:
    """Fit nu-SVR over the nu grid, return (weights, rmse, pearson r)."""
    from sklearn.svm import NuSVR

    best = None
    for nu in NU_GRID:
        model = NuSVR(kernel="linear", nu=nu, C=1.0, tol=1e-3, max_iter=-1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.fit(X, y)
        w = np.asarray(model.coef_).ravel()
        w = np.clip(w, 0, None)
        total = w.sum()
        if total <= 0:
            continue
        w = w / total
        reconstructed = X @ w
        rmse = float(np.sqrt(np.mean((reconstructed - y) ** 2)))
        if np.std(reconstructed) == 0:
            r = 0.0
        else:
            r = float(np.corrcoef(reconstructed, y)[0, 1])
        if best is None or rmse < best[1]:
            best = (w, rmse, r)

    if best is None:
        n = X.shape[1]
        return np.full(n, np.nan), np.nan, np.nan
    return best


def cibersort(
    mixture: pd.DataFrame,
    signature: pd.DataFrame,
    permutations: int = 100,
    random_state: int = 0,
    quantile_normalise_mixture: bool = False,
    verbose: bool = True,
) -> CibersortResult:
    """Run the CIBERSORT core algorithm on every column of ``mixture``.

    Parameters
    ----------
    mixture
        genes x samples in **linear**, per-sample-comparable space (TPM for
        RNA-seq, RMA/MAS5 for arrays). Log-transformed input breaks the linear
        mixing assumption the regression rests on, and the resulting fractions
        are wrong in a way that still looks plausible.
    signature
        genes x cell types signature matrix (LM22 or equivalent).
    permutations
        size of the null distribution for the goodness-of-fit p-value. The
        published default is 100 for exploration and 1,000 for publication.
    quantile_normalise_mixture
        CIBERSORT's own default is TRUE for microarray and FALSE for RNA-seq.
    """
    from .preprocess import quantile_normalise

    if float(np.nanmax(mixture.to_numpy())) < 50:
        raise ValueError(
            "mixture looks log-transformed (max < 50); CIBERSORT requires linear space"
        )

    work = quantile_normalise(mixture) if quantile_normalise_mixture else mixture

    shared = [g for g in signature.index if g in work.index]
    if len(shared) < 50:
        raise ValueError(f"only {len(shared)} signature genes found in the mixture")
    if verbose:
        print(
            f"[cibersort] {len(shared)}/{signature.shape[0]} signature genes matched; "
            f"{signature.shape[1]} cell types; {work.shape[1]} samples"
        )

    X_raw = signature.loc[shared].to_numpy(dtype=float)
    X = (X_raw - X_raw.mean()) / X_raw.std(ddof=1)
    Y = work.loc[shared].to_numpy(dtype=float)

    null_r = None
    if permutations and permutations > 0:
        rng = np.random.default_rng(random_state)
        pool = Y.ravel()
        draws = []
        for _ in range(permutations):
            yr = rng.choice(pool, size=len(shared), replace=False)
            yr = (yr - yr.mean()) / yr.std(ddof=1)
            _, _, r = _core(X, yr)
            draws.append(r)
        null_r = np.sort(np.asarray([d for d in draws if np.isfinite(d)]))
        if verbose:
            print(f"[cibersort] permutation null: n={null_r.size}, max r={null_r.max():.3f}")

    fractions, corrs, rmses, pvals = [], [], [], []
    for j, sample in enumerate(work.columns):
        y = Y[:, j]
        y = (y - y.mean()) / y.std(ddof=1)
        w, rmse, r = _core(X, y)
        fractions.append(w)
        corrs.append(r)
        rmses.append(rmse)
        if null_r is not None and null_r.size:
            pvals.append(float((null_r >= r).sum() + 1) / (null_r.size + 1))
        else:
            pvals.append(np.nan)

    return CibersortResult(
        fractions=pd.DataFrame(fractions, index=work.columns, columns=signature.columns),
        correlation=pd.Series(corrs, index=work.columns, name="correlation"),
        rmse=pd.Series(rmses, index=work.columns, name="rmse"),
        p_value=pd.Series(pvals, index=work.columns, name="p_value"),
        n_genes_used=len(shared),
    )


def constrained_ls_deconvolve(
    mixture: pd.DataFrame,
    signature: pd.DataFrame,
    scaling: pd.Series | None = None,
    simplex: bool = True,
    add_other: bool = True,
) -> pd.DataFrame:
    """Non-negative least squares deconvolution on the simplex (quanTIseq-style).

    Minimise ``||Xw - y||^2`` subject to ``w >= 0`` and, when ``simplex``,
    ``sum(w) <= 1``, so the unexplained mass becomes an explicit
    "Other" compartment (tumour and stroma) instead of being redistributed over
    the immune cell types. That is the structural difference from CIBERSORT
    relative mode, where the fractions sum to one by construction and an
    immune-desert tumour still returns a full set of "fractions".

    ``scaling`` divides each cell type by its mRNA-content factor (quanTIseq's
    TIL10 scaling), converting mRNA proportions into cell proportions.

    This reproduces the estimator *family* used by quanTIseq; it is not the
    quanTIseq pipeline (which also does its own read-level preprocessing). Cite
    and run the original tool if you report quanTIseq numbers.
    """
    from scipy.optimize import minimize, nnls

    shared = [g for g in signature.index if g in mixture.index]
    if len(shared) < 50:
        raise ValueError(f"only {len(shared)} signature genes found in the mixture")
    X = signature.loc[shared].to_numpy(dtype=float)
    Y = mixture.loc[shared].to_numpy(dtype=float)
    n_types = X.shape[1]

    rows = []
    for j in range(Y.shape[1]):
        y = Y[:, j]
        w0, _ = nnls(X, y)
        if simplex and w0.sum() > 1.0:
            result = minimize(
                lambda w: float(np.sum((X @ w - y) ** 2)),
                x0=w0 / max(w0.sum(), 1e-12),
                jac=lambda w: 2.0 * X.T @ (X @ w - y),
                bounds=[(0.0, 1.0)] * n_types,
                constraints=[{"type": "ineq", "fun": lambda w: 1.0 - np.sum(w)}],
                method="SLSQP",
                options={"maxiter": 500, "ftol": 1e-12},
            )
            w0 = np.clip(result.x, 0, None)
        rows.append(w0)

    out = pd.DataFrame(rows, index=mixture.columns, columns=signature.columns)

    if scaling is not None:
        common = [c for c in out.columns if c in scaling.index]
        out[common] = out[common].div(scaling[common], axis=1)

    if add_other:
        out["Other"] = (1.0 - out.sum(axis=1)).clip(lower=0.0)
        total = out.sum(axis=1)
        out = out.div(total.where(total > 0, np.nan), axis=0)
    return out
