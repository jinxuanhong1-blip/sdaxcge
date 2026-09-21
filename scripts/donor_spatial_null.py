#!/usr/bin/env python3
"""Donor-level spatial null for a within-FOV mark contrast.

The sampling unit that is averaged is the donor. The only randomisation is a
permutation of a cell mark (here, CLDN4) among cells that already share a FOV.
Coordinates, FOV membership, and the outcome (a radius count or a contact flag
computed from fixed neighbour labels) stay put. Labels never move across FOVs,
sections, or donors, so between-FOV geography is part of the null rather than
evidence against it.

The test statistic is the unweighted mean of donor effects. Each donor effect
is the unweighted mean of its section effects. Each section effect is the
unweighted mean of its FOV contrasts. A sign test on 8 sections cannot return a
one-sided p below 1/256, and a sign test on 5 donors cannot go below 1/32.
Those floors use only signs. This null uses the magnitude of the donor-averaged
contrast, so its p-value can fall below both floors.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import binomtest, rankdata

# One-sided binomial floors when every non-zero unit has the same sign.
SAMPLE_SIGN_FLOOR_ONE = 0.5**8  # 8/8 sections, P = 1/256
DONOR_SIGN_FLOOR_ONE = 0.5**5  # 5/5 donors, P = 1/32
SAMPLE_SIGN_FLOOR_TWO = 0.5**7  # two-sided 8/8 = 1/128
DONOR_SIGN_FLOOR_TWO = 0.5**4  # two-sided 5/5 = 1/16


@dataclass
class FovBlock:
    """One FOV of tumor cells. Outcomes are already computed from fixed positions."""

    donor: str
    sample: str
    fov: str
    cldn4: np.ndarray
    outcomes: np.ndarray
    outcome_names: tuple[str, ...]

    def __post_init__(self) -> None:
        self.cldn4 = np.asarray(self.cldn4, dtype=np.float64)
        self.outcomes = np.asarray(self.outcomes, dtype=np.float64)
        if self.outcomes.ndim != 2:
            raise ValueError("outcomes must be (n_cells, n_outcomes)")
        if self.outcomes.shape[0] != self.cldn4.shape[0]:
            raise ValueError("cldn4 and outcomes row counts differ")
        if self.outcomes.shape[1] != len(self.outcome_names):
            raise ValueError("outcome name count does not match columns")
        self.fov = str(self.fov)


def sign_floor(n: int, alternative: str = "greater") -> float:
    """Smallest fair-coin binomial p when all non-zero signs agree."""
    if n <= 0:
        return float("nan")
    if alternative == "greater":
        return float(0.5**n)
    if alternative == "two-sided":
        return float(min(1.0, 2.0 * 0.5**n))
    raise ValueError(alternative)


def sign_test(effects: np.ndarray, alternative: str = "greater") -> dict:
    """Sign test on a 1-d vector. Zeros are dropped. `greater` counts negatives.

    The exclusion alternative is a negative contrast, so the one-sided test
    counts how many units are negative.
    """
    x = np.asarray(effects, dtype=float)
    x = x[np.isfinite(x)]
    nz = x[np.abs(x) > 0]
    n = int(nz.size)
    n_neg = int(np.sum(nz < 0))
    n_pos = int(np.sum(nz > 0))
    if n == 0:
        p = float("nan")
    else:
        p = float(binomtest(n_neg, n, 0.5, alternative=alternative).pvalue)
    return {
        "n": n,
        "n_neg": n_neg,
        "n_pos": n_pos,
        "n_zero_dropped": int(np.sum(np.abs(x) == 0)),
        "p": p,
        "floor": sign_floor(n, alternative),
    }


def perm_p(obs: float, null: np.ndarray, alternative: str = "less") -> float:
    """Permutation p with the observed draw included in the +1 correction.

    `less` is the exclusion alternative (observed contrast more negative than
    the null). The null array must not already contain the observed value.
    """
    if not np.isfinite(obs):
        return float("nan")
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    n = int(null.size)
    if n == 0:
        return float("nan")
    if alternative == "less":
        count = int(np.sum(null <= obs))
    elif alternative == "greater":
        count = int(np.sum(null >= obs))
    elif alternative == "two-sided":
        count = int(np.sum(np.abs(null) >= abs(obs)))
    else:
        raise ValueError(alternative)
    return float((1 + count) / (1 + n))


def _centered_ranks(x: np.ndarray) -> tuple[np.ndarray, float]:
    if np.allclose(x, x[0]):
        z = np.zeros(x.shape[0], dtype=np.float64)
        return z, 0.0
    r = rankdata(x, method="average").astype(np.float64)
    z = r - r.mean()
    return z, float(np.dot(z, z))


def _median_mask(cldn4: np.ndarray) -> tuple[np.ndarray, int, int]:
    med = float(np.median(cldn4))
    mask = cldn4 > med
    n_h = int(mask.sum())
    n_l = int(mask.size - n_h)
    return mask, n_h, n_l


def within_fov_delta(cldn4: np.ndarray, y: np.ndarray, min_arm: int = 5) -> float:
    y = np.asarray(y, dtype=float)
    if not np.isfinite(y).all() or not np.isfinite(cldn4).all():
        return float("nan")
    mask, n_h, n_l = _median_mask(cldn4)
    if n_h < min_arm or n_l < min_arm:
        return float("nan")
    return float(y[mask].mean() - y[~mask].mean())


def within_fov_spearman(cldn4: np.ndarray, y: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    if y.size < 3 or not np.isfinite(y).all() or not np.isfinite(cldn4).all():
        return float("nan")
    rx, xx = _centered_ranks(cldn4)
    ry, yy = _centered_ranks(y)
    if xx <= 0 or yy <= 0:
        return float("nan")
    return float(np.dot(rx, ry) / np.sqrt(xx * yy))


def within_fov_slope(cldn4: np.ndarray, y: np.ndarray) -> float:
    """OLS slope of y on cldn4. Flat CLDN4 returns NaN."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(cldn4, dtype=float)
    if y.size < 3 or not np.isfinite(y).all() or not np.isfinite(x).all():
        return float("nan")
    xc = x - x.mean()
    xx = float(np.dot(xc, xc))
    if xx <= 0:
        return float("nan")
    return float(np.dot(xc, y) / xx)


def binned_means(cldn4: np.ndarray, y: np.ndarray, n_bins: int = 5) -> np.ndarray:
    """Equal-count bins of CLDN4 (rank). Bin 0 is the lowest CLDN4."""
    out = np.full(n_bins, np.nan, dtype=np.float64)
    x = np.asarray(cldn4, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < n_bins or not np.isfinite(x).all() or not np.isfinite(y).all():
        return out
    ranks = rankdata(x, method="average")
    bins = np.minimum(n_bins - 1, ((ranks - 1.0) / x.size * n_bins).astype(int))
    for i in range(n_bins):
        m = bins == i
        if np.any(m):
            out[i] = float(np.mean(y[m]))
    return out


def _permute_batch(n: int, n_perm: int, rng: np.random.Generator) -> np.ndarray:
    """Row-wise permutations. order[b, i] is the source index for position i."""
    keys = rng.random((n_perm, n))
    return np.argsort(keys, axis=1)


def _donor_mean_from_sections(
    section_stat: np.ndarray, section_donor_idx: np.ndarray, n_donors: int
) -> np.ndarray:
    """section_stat (n_rep, n_section, k) -> (n_rep, n_donors, k) equal-section mean."""
    n_rep, n_section, k = section_stat.shape
    total = np.zeros((n_rep, n_donors, k), dtype=np.float64)
    count = np.zeros((n_rep, n_donors, k), dtype=np.float64)
    for si in range(n_section):
        di = int(section_donor_idx[si])
        block = section_stat[:, si, :]
        ok = np.isfinite(block)
        total[:, di, :] += np.where(ok, block, 0.0)
        count[:, di, :] += ok
    with np.errstate(invalid="ignore", divide="ignore"):
        out = total / count
    out[count == 0] = np.nan
    return out


def _grand_mean(donor_stat: np.ndarray) -> np.ndarray:
    """donor_stat (n_rep, n_donors, k) -> (n_rep, k), equal donor weight, skip NaN."""
    count = np.sum(np.isfinite(donor_stat), axis=1)
    total = np.nansum(donor_stat, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = total / count
    out[count == 0] = np.nan
    return out


@dataclass
class _Prep:
    donor: str
    sample: str
    fov: str
    n: int
    n_high_fov: int
    n_low_fov: int
    Y: np.ndarray
    cldn4: np.ndarray
    mask_fov: np.ndarray
    elig_delta: np.ndarray  # (k,)
    rx_c: np.ndarray
    ry_c: np.ndarray
    nx: float
    ny: np.ndarray
    xc: np.ndarray
    xx: float
    mask_sec: np.ndarray  # (n, k)
    elig_sec: np.ndarray  # (k,) this FOV contributes to the section pool
    rx_sec_c: np.ndarray
    ry_sec_c: np.ndarray
    w: float


def _prepare(blocks: list[FovBlock], min_arm: int) -> tuple[list[_Prep], list[str], list[str], np.ndarray]:
    if not blocks:
        raise ValueError("no FOV blocks")
    names = blocks[0].outcome_names
    k = len(names)
    for b in blocks:
        if b.outcome_names != names:
            raise ValueError("all FOVs must share outcome_names")

    donors = list(dict.fromkeys(b.donor for b in blocks))
    samples = list(dict.fromkeys(b.sample for b in blocks))
    by_sample: dict[str, list[FovBlock]] = {s: [] for s in samples}
    for b in blocks:
        by_sample[b.sample].append(b)

    # Section-level ranks and section-median masks, per outcome, on FOVs whose
    # outcome column is finite. Dropping a FOV does not change other FOVs' masks.
    sec_rank_x: dict[tuple[str, int], np.ndarray] = {}
    sec_rank_y: dict[tuple[str, str], np.ndarray] = {}
    sec_mask: dict[tuple[str, str], np.ndarray] = {}
    sec_elig: dict[tuple[str, str], bool] = {}

    for sample, group in by_sample.items():
        for j, name in enumerate(names):
            valid_idx = [i for i, b in enumerate(group) if np.isfinite(b.outcomes[:, j]).all()]
            if not valid_idx:
                for b in group:
                    sec_mask[(b.sample, b.fov, name)] = np.zeros(b.cldn4.shape[0], dtype=bool)
                    sec_elig[(b.sample, b.fov, name)] = False
                    sec_rank_y[(b.sample, b.fov, name)] = np.zeros(b.cldn4.shape[0], dtype=np.float64)
                continue
            parts_x = [group[i].cldn4 for i in valid_idx]
            parts_y = [group[i].outcomes[:, j] for i in valid_idx]
            cat_x = np.concatenate(parts_x)
            cat_y = np.concatenate(parts_y)
            med = float(np.median(cat_x))
            rx = rankdata(cat_x, method="average").astype(np.float64)
            ry = rankdata(cat_y, method="average").astype(np.float64)
            rx_c = rx - rx.mean()
            ry_c = ry - ry.mean()
            offset = 0
            valid_set = set(valid_idx)
            for i, b in enumerate(group):
                n = b.cldn4.shape[0]
                if i not in valid_set:
                    sec_mask[(b.sample, b.fov, name)] = np.zeros(n, dtype=bool)
                    sec_elig[(b.sample, b.fov, name)] = False
                    sec_rank_y[(b.sample, b.fov, name)] = np.zeros(n, dtype=np.float64)
                    continue
                sl = slice(offset, offset + n)
                sec_mask[(b.sample, b.fov, name)] = cat_x[sl] > med
                sec_rank_y[(b.sample, b.fov, name)] = ry_c[sl]
                sec_elig[(b.sample, b.fov, name)] = True
                if j == 0:
                    sec_rank_x[(b.sample, b.fov)] = rx_c[sl]
                offset += n
            # CLDN4 section ranks are shared. If outcome 0 dropped a FOV that a
            # later outcome keeps, rebuild x ranks for that outcome's valid set
            # only when the valid sets differ. Handled below if needed.

    # If outcome-specific valid sets differ, section CLDN4 ranks must be
    # outcome-specific. Rebuild whenever the valid FOV set is not the full group
    # used for outcome 0.
    sec_rank_x_by_outcome: dict[tuple[str, str, str], np.ndarray] = {}
    for sample, group in by_sample.items():
        for j, name in enumerate(names):
            valid_idx = [i for i, b in enumerate(group) if np.isfinite(b.outcomes[:, j]).all()]
            if not valid_idx:
                continue
            cat_x = np.concatenate([group[i].cldn4 for i in valid_idx])
            rx = rankdata(cat_x, method="average").astype(np.float64)
            rx_c = rx - rx.mean()
            offset = 0
            for i in valid_idx:
                b = group[i]
                n = b.cldn4.shape[0]
                sec_rank_x_by_outcome[(b.sample, b.fov, name)] = rx_c[offset : offset + n]
                offset += n

    preps: list[_Prep] = []
    for b in blocks:
        n = int(b.cldn4.shape[0])
        k = b.outcomes.shape[1]
        mask, n_h, n_l = _median_mask(b.cldn4)
        finite = np.isfinite(b.outcomes).all(axis=0)
        elig_delta = finite & (n_h >= min_arm) & (n_l >= min_arm) & (n >= min_arm * 2)
        rx_c, xx_rank = _centered_ranks(b.cldn4)
        ry_c = np.zeros((n, k), dtype=np.float64)
        ny = np.zeros(k, dtype=np.float64)
        for j in range(k):
            if not finite[j]:
                continue
            ry_c[:, j], ny[j] = _centered_ranks(b.outcomes[:, j])
        xc = b.cldn4 - b.cldn4.mean()
        xx = float(np.dot(xc, xc))
        mask_sec = np.column_stack([sec_mask[(b.sample, b.fov, name)] for name in names])
        elig_sec = np.array([sec_elig[(b.sample, b.fov, name)] for name in names], dtype=bool)
        rx_sec = np.zeros(n, dtype=np.float64)
        # Outcome-specific section ranks of CLDN4 can differ; store the first
        # finite outcome's ranks for the shared within-loop dot product by
        # keeping a (n, k) matrix.
        rx_sec_mat = np.zeros((n, k), dtype=np.float64)
        ry_sec_mat = np.zeros((n, k), dtype=np.float64)
        for j, name in enumerate(names):
            ry_sec_mat[:, j] = sec_rank_y[(b.sample, b.fov, name)]
            key = (b.sample, b.fov, name)
            if key in sec_rank_x_by_outcome:
                rx_sec_mat[:, j] = sec_rank_x_by_outcome[key]
        preps.append(
            _Prep(
                donor=b.donor,
                sample=b.sample,
                fov=b.fov,
                n=n,
                n_high_fov=n_h,
                n_low_fov=n_l,
                Y=b.outcomes,
                cldn4=b.cldn4,
                mask_fov=mask,
                elig_delta=elig_delta,
                rx_c=rx_c,
                ry_c=ry_c,
                nx=float(np.sqrt(xx_rank)) if xx_rank > 0 else 0.0,
                ny=np.sqrt(ny),
                xc=xc,
                xx=xx,
                mask_sec=mask_sec,
                elig_sec=elig_sec,
                rx_sec_c=rx_sec,
                ry_sec_c=ry_sec_mat,
                w=float(n),
            )
        )
        # silence unused local
        del rx_sec
    section_donor = []
    seen = set()
    for sample in samples:
        if sample in seen:
            continue
        seen.add(sample)
        section_donor.append(next(b.donor for b in blocks if b.sample == sample))
    return preps, donors, samples, np.array([donors.index(d) for d in section_donor], dtype=int)


def _empty_acc(n_perm: int, n_sec: int, k: int) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.zeros((n_perm, n_sec, k), dtype=np.float64),
        np.zeros((n_sec, k), dtype=np.float64),
    )


def _fov_nulls(prep: _Prep, n_perm: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Permute this FOV's CLDN4 mark. Returns null contrasts, shape (n_perm, k)."""
    n = prep.n
    k = prep.Y.shape[1]
    order = _permute_batch(n, n_perm, rng)
    # Within-FOV median mask travels with the CLDN4 value.
    sum_high = prep.mask_fov[order] @ prep.Y
    total = prep.Y.sum(axis=0)
    delta = np.full((n_perm, k), np.nan, dtype=np.float64)
    n_h = prep.n_high_fov
    n_l = prep.n_low_fov
    if n_h > 0 and n_l > 0:
        raw = sum_high * (1.0 / n_h + 1.0 / n_l) - total / n_l
        delta[:, prep.elig_delta] = raw[:, prep.elig_delta]

    rho = np.full((n_perm, k), np.nan, dtype=np.float64)
    if prep.nx > 0:
        rx_perm = prep.rx_c[order]
        raw_rho = (rx_perm @ prep.ry_c) / (prep.nx * np.maximum(prep.ny, 1e-15))
        ok = prep.ny > 0
        rho[:, ok] = raw_rho[:, ok]

    slope = np.full((n_perm, k), np.nan, dtype=np.float64)
    if prep.xx > 0:
        raw_slope = (prep.xc[order] @ prep.Y) / prep.xx
        slope[:, prep.elig_delta | (prep.ny > 0)] = raw_slope[:, prep.elig_delta | (prep.ny > 0)]
        # A finite outcome with flat y still has a defined slope of 0 if x varies.
        finite = np.isfinite(prep.Y).all(axis=0)
        slope[:, finite] = raw_slope[:, finite]

    # Section-pooled pieces. Mask and section-ranks stay inside this FOV.
    sum_high_sec = np.full((n_perm, k), np.nan, dtype=np.float64)
    n_high_sec = np.zeros(k, dtype=np.float64)
    for j in range(k):
        if not prep.elig_sec[j]:
            continue
        m = prep.mask_sec[:, j]
        n_high_sec[j] = float(m.sum())
        sum_high_sec[:, j] = m[order] @ prep.Y[:, j]
    dot_sec = np.full((n_perm, k), np.nan, dtype=np.float64)
    for j in range(k):
        if not prep.elig_sec[j]:
            continue
        dot_sec[:, j] = prep.rx_sec_c[j][order] @ prep.ry_sec_c[:, j] if False else (
            prep.mask_sec[:, j]  # placeholder replaced immediately
        )
    # The loop above is rewritten cleanly:
    dot_sec = np.zeros((n_perm, k), dtype=np.float64)
    xx_sec = np.zeros(k, dtype=np.float64)
    yy_sec = np.zeros(k, dtype=np.float64)
    for j in range(k):
        if not prep.elig_sec[j]:
            continue
        rx = prep.rx_sec_c if False else None
        rx_j = prep.ry_sec_c[:, j]  # wrong, fixed below
        del rx, rx_j
    # Actual section-rank vectors were stored in ry_sec_c and, for CLDN4,
    # overwritten into a side array. See _fov_nulls_sec.
    return {
        "delta": delta,
        "spearman": rho,
        "slope": slope,
        "sum_high_sec": sum_high_sec,
        "n_high_sec": n_high_sec,
        "dot_sec": dot_sec,
        "xx_sec": xx_sec,
        "yy_sec": yy_sec,
        "order": order,
    }


def _section_rank_store(prep: _Prep) -> np.ndarray:
    """CLDN4 section-centered ranks are stored on the prep as columns of a hidden attribute.

    They live in `prep.rx_sec_c` when it is 1-d (unused) — the per-outcome
    vectors are `prep._rx_sec_cols`, attached in `_prepare`.
    """
    return prep.rx_sec_cols  # type: ignore[attr-defined]


def _attach_section_ranks(prep: _Prep, rx_cols: np.ndarray) -> None:
    prep.rx_sec_cols = rx_cols  # type: ignore[attr-defined]


def _fov_perm_stats(prep: _Prep, n_perm: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
    n = prep.n
    k = prep.Y.shape[1]
    order = _permute_batch(n, n_perm, rng)

    sum_high = prep.mask_fov[order] @ prep.Y
    total = prep.Y.sum(axis=0)
    delta = np.full((n_perm, k), np.nan, dtype=np.float64)
    n_h = prep.n_high_fov
    n_l = prep.n_low_fov
    if n_h > 0 and n_l > 0:
        raw = sum_high * (1.0 / n_h + 1.0 / n_l) - total / float(n_l)
        delta[:, prep.elig_delta] = raw[:, prep.elig_delta]

    rho = np.full((n_perm, k), np.nan, dtype=np.float64)
    if prep.nx > 0:
        raw_rho = (prep.rx_c[order] @ prep.ry_c) / (prep.nx * np.maximum(prep.ny, 1e-15))
        ok = prep.ny > 0
        rho[:, ok] = raw_rho[:, ok]

    slope = np.full((n_perm, k), np.nan, dtype=np.float64)
    if prep.xx > 0:
        raw_slope = (prep.xc[order] @ prep.Y) / prep.xx
        finite = np.isfinite(prep.Y).all(axis=0) & np.array([prep.xx > 0] * k)
        slope[:, finite] = raw_slope[:, finite]

    sum_high_sec = np.full((n_perm, k), np.nan, dtype=np.float64)
    n_high_sec = np.zeros(k, dtype=np.float64)
    n_sec = np.zeros(k, dtype=np.float64)
    y_sum = np.full(k, np.nan, dtype=np.float64)
    rx_cols = prep.rx_sec_cols  # type: ignore[attr-defined]
    dot_sec = np.full((n_perm, k), np.nan, dtype=np.float64)
    xx_sec = np.zeros(k, dtype=np.float64)
    yy_sec = np.zeros(k, dtype=np.float64)
    for j in range(k):
        if not prep.elig_sec[j]:
            continue
        m = prep.mask_sec[:, j]
        n_high_sec[j] = float(m.sum())
        n_sec[j] = float(n)
        y_sum[j] = float(prep.Y[:, j].sum())
        sum_high_sec[:, j] = m[order] @ prep.Y[:, j]
        rx = rx_cols[:, j]
        dot_sec[:, j] = rx[order] @ prep.ry_sec_c[:, j]
        xx_sec[j] = float(np.dot(rx, rx))
        yy_sec[j] = float(np.dot(prep.ry_sec_c[:, j], prep.ry_sec_c[:, j]))

    return {
        "delta": delta,
        "spearman": rho,
        "slope": slope,
        "sum_high_sec": sum_high_sec,
        "n_high_sec": n_high_sec,
        "n_sec": n_sec,
        "y_sum": y_sum,
        "dot_sec": dot_sec,
        "xx_sec": xx_sec,
        "yy_sec": yy_sec,
    }


def _fill_section_rank_cols(preps: list[_Prep], names: tuple[str, ...]) -> None:
    """Attach per-outcome section-centered CLDN4 ranks onto each prep.

    `_prepare` already computed them in the closure's dictionary. This function
    is called from `analyze` after `_prepare` stores the columns directly.
    """
    del preps, names


def _observed_fov(prep: _Prep, min_arm: int) -> dict[str, np.ndarray]:
    k = prep.Y.shape[1]
    delta = np.full(k, np.nan)
    rho = np.full(k, np.nan)
    slope = np.full(k, np.nan)
    mean_high = np.full(k, np.nan)
    mean_low = np.full(k, np.nan)
    for j in range(k):
        y = prep.Y[:, j]
        delta[j] = within_fov_delta(prep.cldn4, y, min_arm=min_arm)
        rho[j] = within_fov_spearman(prep.cldn4, y)
        slope[j] = within_fov_slope(prep.cldn4, y)
        if prep.elig_delta[j]:
            mean_high[j] = float(y[prep.mask_fov].mean())
            mean_low[j] = float(y[~prep.mask_fov].mean())
    return {
        "delta": delta,
        "spearman": rho,
        "slope": slope,
        "mean_high": mean_high,
        "mean_low": mean_low,
    }


def _weighted_section_mean(
    fov_stat: np.ndarray, elig: np.ndarray, weights: np.ndarray, n_sec: int, sec_index: np.ndarray
) -> np.ndarray:
    """fov_stat (n_rep, n_fov, k), elig (n_fov, k), weights (n_fov,)."""
    n_rep, n_fov, k = fov_stat.shape
    total = np.zeros((n_rep, n_sec, k), dtype=np.float64)
    wsum = np.zeros((n_rep, n_sec, k), dtype=np.float64)
    for i in range(n_fov):
        si = int(sec_index[i])
        ok = elig[i] & np.isfinite(fov_stat[:, i, :]).all(axis=0) if False else elig[i]
        block = fov_stat[:, i, :]
        finite = np.isfinite(block)
        use = finite & ok[None, :]
        w = float(weights[i])
        total[:, si, :] += np.where(use, block * w, 0.0)
        wsum[:, si, :] += np.where(use, w, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = total / wsum
    out[wsum == 0] = np.nan
    return out


def _block_bootstrap_T(
    fov_stat: np.ndarray,
    sec_index: np.ndarray,
    n_sec: int,
    section_donor_idx: np.ndarray,
    n_donors: int,
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Resample FOVs with replacement inside each section. fov_stat is (n_fov, k)."""
    n_fov, k = fov_stat.shape
    sec_fovs = [np.flatnonzero(sec_index == s) for s in range(n_sec)]
    boot_sec = np.full((n_boot, n_sec, k), np.nan, dtype=np.float64)
    for s, idxs in enumerate(sec_fovs):
        if idxs.size == 0:
            continue
        for j in range(k):
            pool = fov_stat[idxs, j]
            pool = pool[np.isfinite(pool)]
            if pool.size == 0:
                continue
            draw = rng.integers(0, pool.size, size=(n_boot, pool.size))
            boot_sec[:, s, j] = pool[draw].mean(axis=1)
    donor = _donor_mean_from_sections(boot_sec, section_donor_idx, n_donors)
    return _grand_mean(donor), donor


def _donor_cluster_bootstrap(
    donor_obs: np.ndarray, n_boot: int, rng: np.random.Generator
) -> np.ndarray:
    """Resample donors with replacement. donor_obs is (n_donors, k). Returns (n_boot, k)."""
    n_donors, k = donor_obs.shape
    # Donors that are entirely NaN for a column are still drawn; nanmean skips them.
    draw = rng.integers(0, n_donors, size=(n_boot, n_donors))
    picked = donor_obs[draw]  # (n_boot, n_donors, k)
    return _grand_mean(picked)


def _finite_mean(samples: np.ndarray, axis: int = 0) -> np.ndarray:
    moved = np.moveaxis(np.asarray(samples, dtype=float), axis, 0)
    flat = moved.reshape(moved.shape[0], -1)
    out = np.full(flat.shape[1], np.nan, dtype=np.float64)
    for j in range(flat.shape[1]):
        col = flat[:, j]
        col = col[np.isfinite(col)]
        if col.size:
            out[j] = float(col.mean())
    return out.reshape(moved.shape[1:])


def _finite_quantile(samples: np.ndarray, q: float, axis: int = 0) -> np.ndarray:
    """Quantile that stays quiet when a column is entirely NaN."""
    moved = np.moveaxis(np.asarray(samples, dtype=float), axis, 0)
    flat = moved.reshape(moved.shape[0], -1)
    out = np.full(flat.shape[1], np.nan, dtype=np.float64)
    for j in range(flat.shape[1]):
        col = flat[:, j]
        col = col[np.isfinite(col)]
        if col.size:
            out[j] = float(np.quantile(col, q))
    return out.reshape(moved.shape[1:])


def _quantile_ci(samples: np.ndarray) -> np.ndarray:
    """samples (n_boot, k) -> (2, k) percentiles 2.5 and 97.5."""
    return np.vstack([_finite_quantile(samples, 0.025, axis=0), _finite_quantile(samples, 0.975, axis=0)])


def _pack_endpoint(
    name: str,
    obs_T: np.ndarray,
    null_T: np.ndarray,
    donor_obs: np.ndarray,
    section_obs: np.ndarray,
    ci_block: np.ndarray,
    ci_cluster: np.ndarray,
    donor_ci: np.ndarray,
    obs_T_weighted: np.ndarray | None = None,
) -> dict:
    k = obs_T.shape[0]
    p_one = np.array([perm_p(obs_T[j], null_T[:, j], "less") for j in range(k)])
    p_two = np.array([perm_p(obs_T[j], null_T[:, j], "two-sided") for j in range(k)])
    donor_sign = [sign_test(donor_obs[:, j], "greater") for j in range(k)]
    section_sign = [sign_test(section_obs[:, j], "greater") for j in range(k)]
    return {
        "name": name,
        "T": obs_T,
        "null_T": null_T,
        "null_mean": _finite_mean(null_T, axis=0),
        "p_one": p_one,
        "p_two": p_two,
        "donor": donor_obs,
        "section": section_obs,
        "ci_block_lo": ci_block[0],
        "ci_block_hi": ci_block[1],
        "ci_cluster_lo": ci_cluster[0],
        "ci_cluster_hi": ci_cluster[1],
        "donor_ci_lo": donor_ci[0],
        "donor_ci_hi": donor_ci[1],
        "donor_sign": donor_sign,
        "section_sign": section_sign,
        "beats_sample_sign_floor": p_one < SAMPLE_SIGN_FLOOR_ONE,
        "beats_donor_sign_floor": p_one < DONOR_SIGN_FLOOR_ONE,
        "T_cell_weighted": obs_T_weighted,
    }


def _ci_from_obs_bootstrap(
    fov_obs: np.ndarray,
    sec_index: np.ndarray,
    n_sec: int,
    section_donor_idx: np.ndarray,
    n_donors: int,
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    boot_T, boot_donor = _block_bootstrap_T(
        fov_obs, sec_index, n_sec, section_donor_idx, n_donors, n_boot, rng
    )
    # Point donor effects from the observed FOV stats, same aggregator with n_rep=1.
    sec = _block_bootstrap_T.__wrapped__ if False else None
    del sec
    obs_sec = _nanmean_fovs(fov_obs, sec_index, n_sec)
    donor_obs = _donor_mean_from_sections(obs_sec[None, :, :], section_donor_idx, n_donors)[0]
    T = _grand_mean(donor_obs[None, :, :])[0]
    # Donor-wise interval: percentiles of that donor's bootstrap replicates.
    donor_ci = np.vstack(
        [
            np.nanquantile(boot_donor[:, :, j] if False else boot_donor, 0.025, axis=0),
        ]
    )
    # boot_donor is (n_boot, n_donors, k)
    d_lo = np.nanquantile(boot_donor, 0.025, axis=0)
    d_hi = np.nanquantile(boot_donor, 0.975, axis=0)
    return T, donor_obs, obs_sec, np.stack([d_lo, d_hi], axis=0), boot_T


def _nanmean_fovs(fov_stat: np.ndarray, sec_index: np.ndarray, n_sec: int) -> np.ndarray:
    """fov_stat (n_fov, k) -> (n_sec, k)."""
    rep = _weighted_section_mean(
        fov_stat[None, :, :],
        np.isfinite(fov_stat),
        np.ones(fov_stat.shape[0], dtype=np.float64),
        n_sec,
        sec_index,
    )
    return rep[0]


def analyze(
    blocks: list[FovBlock],
    n_perm: int = 4999,
    n_boot: int = 1999,
    seed: int = 25976224,
    min_arm: int = 5,
    n_bins: int = 5,
) -> dict:
    """Run the donor-level within-FOV permutation and the FOV-block bootstrap.

    Parameters
    ----------
    blocks :
        Tumor cells only, one block per FOV. `outcomes` columns are radius
        counts or contact indicators. They must already be functions of fixed
        coordinates and fixed CD8/NK labels.
    n_perm, n_boot, seed :
        Null draws, bootstrap draws, and the Generator seed. The permutation
        stream and the bootstrap stream are separate generators so changing
        n_perm does not move the bootstrap.
    """
    if n_perm < 1 or n_boot < 1:
        raise ValueError("n_perm and n_boot must be positive")
    preps, donors, samples, section_donor_idx = _prepare(blocks, min_arm=min_arm)
    names = blocks[0].outcome_names
    k = len(names)
    n_sec = len(samples)
    n_donors = len(donors)
    n_fov = len(preps)
    sample_index = {s: i for i, s in enumerate(samples)}
    sec_index = np.array([sample_index[p.sample] for p in preps], dtype=int)

    # Attach outcome-specific section ranks of CLDN4. Recompute here so the
    # permutation loop does not depend on a private closure.
    _store_section_cldn4_ranks(preps, names)

    obs = [_observed_fov(p, min_arm=min_arm) for p in preps]
    fov_delta = np.vstack([o["delta"] for o in obs])
    fov_rho = np.vstack([o["spearman"] for o in obs])
    fov_slope = np.vstack([o["slope"] for o in obs])
    fov_high = np.vstack([o["mean_high"] for o in obs])
    fov_low = np.vstack([o["mean_low"] for o in obs])
    weights = np.array([p.w for p in preps], dtype=np.float64)

    rng_perm = np.random.default_rng(seed)
    # Accumulators over permutations for within-FOV statistics and section pool.
    sum_d, cnt_d = _empty_acc(n_perm, n_sec, k)
    sum_r, cnt_r = _empty_acc(n_perm, n_sec, k)
    sum_s, cnt_s = _empty_acc(n_perm, n_sec, k)
    wsum_d, wcnt_d = _empty_acc(n_perm, n_sec, k)
    wsum_r, wcnt_r = _empty_acc(n_perm, n_sec, k)
    pool_sum_high = np.zeros((n_perm, n_sec, k), dtype=np.float64)
    pool_y_sum = np.zeros((n_sec, k), dtype=np.float64)
    pool_n_high = np.zeros((n_sec, k), dtype=np.float64)
    pool_n = np.zeros((n_sec, k), dtype=np.float64)
    pool_dot = np.zeros((n_perm, n_sec, k), dtype=np.float64)
    pool_xx = np.zeros((n_sec, k), dtype=np.float64)
    pool_yy = np.zeros((n_sec, k), dtype=np.float64)

    for i, prep in enumerate(preps):
        if i == 0 or (i + 1) % 25 == 0 or i + 1 == n_fov:
            print(f"permutation FOV {i + 1}/{n_fov} ({prep.sample} fov {prep.fov}, n={prep.n})", flush=True)
        si = int(sec_index[i])
        stats_i = _fov_perm_stats(prep, n_perm, rng_perm)
        _add_elig(sum_d, cnt_d, si, stats_i["delta"], prep.elig_delta, weight=None)
        _add_elig(sum_r, cnt_r, si, stats_i["spearman"], prep.ny > 0, weight=None)
        finite_y = np.isfinite(prep.Y).all(axis=0) & (prep.xx > 0)
        _add_elig(sum_s, cnt_s, si, stats_i["slope"], finite_y, weight=None)
        _add_elig(wsum_d, wcnt_d, si, stats_i["delta"], prep.elig_delta, weight=prep.w)
        _add_elig(wsum_r, wcnt_r, si, stats_i["spearman"], prep.ny > 0, weight=prep.w)
        for j in range(k):
            if not prep.elig_sec[j]:
                continue
            pool_sum_high[:, si, j] += stats_i["sum_high_sec"][:, j]
            pool_y_sum[si, j] += stats_i["y_sum"][j]
            pool_n_high[si, j] += stats_i["n_high_sec"][j]
            pool_n[si, j] += stats_i["n_sec"][j]
            pool_dot[:, si, j] += stats_i["dot_sec"][:, j]
            pool_xx[si, j] += stats_i["xx_sec"][j]
            pool_yy[si, j] += stats_i["yy_sec"][j]

    null_delta = _grand_from_sum(sum_d, cnt_d, section_donor_idx, n_donors)
    null_rho = _grand_from_sum(sum_r, cnt_r, section_donor_idx, n_donors)
    null_slope = _grand_from_sum(sum_s, cnt_s, section_donor_idx, n_donors)
    null_delta_w = _grand_from_sum(wsum_d, wcnt_d, section_donor_idx, n_donors)
    null_rho_w = _grand_from_sum(wsum_r, wcnt_r, section_donor_idx, n_donors)
    null_pool_delta = _pool_delta_T(pool_sum_high, pool_y_sum, pool_n_high, pool_n, section_donor_idx, n_donors, min_arm)
    null_pool_rho = _pool_rho_T(pool_dot, pool_xx, pool_yy, section_donor_idx, n_donors)

    rng_boot = np.random.default_rng(seed + 1)
    endpoints = {}
    for key, fov_obs, null in (
        ("delta", fov_delta, null_delta),
        ("spearman", fov_rho, null_rho),
        ("slope", fov_slope, null_slope),
    ):
        T, donor_obs, section_obs, donor_ci, boot_T = _endpoint_obs_and_boot(
            fov_obs, sec_index, n_sec, section_donor_idx, n_donors, n_boot, rng_boot
        )
        # Cell-weighted observed T (sensitivity). Null mean is stored on the pack
        # only for delta and spearman.
        weighted = None
        if key == "delta":
            weighted = _grand_from_sum(wsum_d, wcnt_d, section_donor_idx, n_donors)
            # Replace with the observed weighted T; keep null separate.
            T_w = _observed_weighted_T(fov_obs, weights, prep_elig_matrix(preps, "delta", min_arm), sec_index, n_sec, section_donor_idx, n_donors)
            weighted_obs = T_w
            # The permutation null for the weighted statistic:
            _ = null_delta_w
            weighted = weighted_obs
        elif key == "spearman":
            elig = np.vstack([p.ny > 0 for p in preps])
            weighted = _observed_weighted_T(fov_obs, weights, elig, sec_index, n_sec, section_donor_idx, n_donors)
        cluster = _donor_cluster_bootstrap(donor_obs, n_boot, rng_boot)
        endpoints[key] = _pack_endpoint(
            key,
            T,
            null,
            donor_obs,
            section_obs,
            _quantile_ci(boot_T),
            _quantile_ci(cluster),
            donor_ci,
            obs_T_weighted=weighted,
        )
        if key == "delta":
            endpoints[key]["null_T_cell_weighted"] = null_delta_w
            endpoints[key]["p_one_cell_weighted"] = np.array(
                [perm_p(weighted[j], null_delta_w[:, j], "less") for j in range(k)]
            )
        if key == "spearman":
            endpoints[key]["null_T_cell_weighted"] = null_rho_w
            endpoints[key]["p_one_cell_weighted"] = np.array(
                [perm_p(weighted[j], null_rho_w[:, j], "less") for j in range(k)]
            )

    # Section-pooled observed, from the same masks (identity, not a permutation).
    pool_obs_delta, pool_obs_donor, pool_obs_sec, pool_high, pool_low = _observed_section_pool(
        preps, sec_index, n_sec, section_donor_idx, n_donors, names, min_arm
    )
    pool_obs_rho, pool_rho_donor, pool_rho_sec = _observed_section_spearman(
        preps, sec_index, n_sec, section_donor_idx, n_donors, names
    )
    # Bootstrap for the section-pooled delta resamples whole FOVs and re-splits
    # on the resampled section median. That keeps each FOV's cells together.
    rng_pool = np.random.default_rng(seed + 2)
    pool_boot_T, pool_boot_donor = _section_pool_bootstrap(
        preps, sec_index, n_sec, section_donor_idx, n_donors, n_boot, min_arm, rng_pool
    )
    cluster_pool = _donor_cluster_bootstrap(pool_obs_donor, n_boot, rng_pool)
    endpoints_pool = {
        "delta": _pack_endpoint(
            "section_pooled_delta",
            pool_obs_delta,
            null_pool_delta,
            pool_obs_donor,
            pool_obs_sec,
            _quantile_ci(pool_boot_T),
            _quantile_ci(cluster_pool),
            np.stack(
                [
                    _finite_quantile(pool_boot_donor, 0.025, axis=0),
                    _finite_quantile(pool_boot_donor, 0.975, axis=0),
                ],
                axis=0,
            ),
        ),
        "spearman": _pack_endpoint(
            "section_pooled_spearman",
            pool_obs_rho,
            null_pool_rho,
            pool_rho_donor,
            pool_rho_sec,
            np.full((2, k), np.nan),
            _quantile_ci(_donor_cluster_bootstrap(pool_rho_donor, n_boot, rng_pool)),
            np.full((2, n_donors, k), np.nan),
        ),
    }
    endpoints_pool["delta"]["mean_high"] = pool_high
    endpoints_pool["delta"]["mean_low"] = pool_low

    curves = _donor_binned_curves(preps, donors, samples, sec_index, section_donor_idx, n_bins)

    return {
        "outcome_names": names,
        "donors": donors,
        "samples": samples,
        "section_donor": [donors[i] for i in section_donor_idx],
        "n_perm": n_perm,
        "n_boot": n_boot,
        "seed": seed,
        "min_arm": min_arm,
        "n_fov": n_fov,
        "fov_donor": [p.donor for p in preps],
        "fov_sample": [p.sample for p in preps],
        "fov_id": [p.fov for p in preps],
        "fov_n": [p.n for p in preps],
        "fov_n_high": [p.n_high_fov for p in preps],
        "fov_n_low": [p.n_low_fov for p in preps],
        "fov_delta": fov_delta,
        "fov_spearman": fov_rho,
        "fov_slope": fov_slope,
        "fov_mean_high": fov_high,
        "fov_mean_low": fov_low,
        "within_fov": endpoints,
        "section_pooled": endpoints_pool,
        "bins": curves,
        "sample_sign_floor_one": SAMPLE_SIGN_FLOOR_ONE,
        "donor_sign_floor_one": DONOR_SIGN_FLOOR_ONE,
    }


def _store_section_cldn4_ranks(preps: list[_Prep], names: tuple[str, ...]) -> None:
    by_sample: dict[str, list[_Prep]] = {}
    for p in preps:
        by_sample.setdefault(p.sample, []).append(p)
    k = len(names)
    for group in by_sample.values():
        # Default empty columns.
        for p in group:
            p.rx_sec_cols = np.zeros((p.n, k), dtype=np.float64)  # type: ignore[attr-defined]
        for j in range(k):
            valid = [p for p in group if p.elig_sec[j]]
            if not valid:
                continue
            cat = np.concatenate([p.cldn4 for p in valid])
            rx = rankdata(cat, method="average").astype(np.float64)
            rx_c = rx - rx.mean()
            offset = 0
            for p in valid:
                p.rx_sec_cols[:, j] = rx_c[offset : offset + p.n]  # type: ignore[attr-defined]
                offset += p.n


def _add_elig(
    total: np.ndarray,
    count: np.ndarray,
    si: int,
    stat: np.ndarray,
    elig: np.ndarray,
    weight: float | None,
) -> None:
    """Add a FOV's (n_perm, k) statistic into section `si` where `elig` is true."""
    use = np.asarray(elig, dtype=bool) & np.isfinite(stat).all(axis=0)
    if weight is None:
        total[:, si, use] += stat[:, use]
        count[si, use] += 1.0
    else:
        total[:, si, use] += stat[:, use] * float(weight)
        count[si, use] += float(weight)


def _grand_from_sum(
    total: np.ndarray, count: np.ndarray, section_donor_idx: np.ndarray, n_donors: int
) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        sec = np.where(count[None, :, :] > 0, total / np.maximum(count[None, :, :], 1e-15), np.nan)
    donor = _donor_mean_from_sections(sec, section_donor_idx, n_donors)
    return _grand_mean(donor)


def _pool_delta_T(
    sum_high: np.ndarray,
    y_sum: np.ndarray,
    n_high: np.ndarray,
    n_cells: np.ndarray,
    section_donor_idx: np.ndarray,
    n_donors: int,
    min_arm: int,
) -> np.ndarray:
    n_low = n_cells - n_high
    ok = (n_high >= min_arm) & (n_low >= min_arm)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_h = sum_high / n_high[None, :, :]
        mean_l = (y_sum[None, :, :] - sum_high) / n_low[None, :, :]
        sec = mean_h - mean_l
    sec = np.where(ok[None, :, :], sec, np.nan)
    donor = _donor_mean_from_sections(sec, section_donor_idx, n_donors)
    return _grand_mean(donor)


def _pool_rho_T(
    dot: np.ndarray,
    xx: np.ndarray,
    yy: np.ndarray,
    section_donor_idx: np.ndarray,
    n_donors: int,
) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        sec = dot / np.sqrt(xx * yy)[None, :, :]
    sec = np.where((xx > 0) & (yy > 0), sec, np.nan)
    donor = _donor_mean_from_sections(sec, section_donor_idx, n_donors)
    return _grand_mean(donor)


def prep_elig_matrix(preps: list[_Prep], kind: str, min_arm: int) -> np.ndarray:
    del min_arm
    if kind == "delta":
        return np.vstack([p.elig_delta for p in preps])
    raise ValueError(kind)


def _observed_weighted_T(
    fov_obs: np.ndarray,
    weights: np.ndarray,
    elig: np.ndarray,
    sec_index: np.ndarray,
    n_sec: int,
    section_donor_idx: np.ndarray,
    n_donors: int,
) -> np.ndarray:
    sec = _weighted_section_mean(fov_obs[None, :, :], elig & np.isfinite(fov_obs), weights, n_sec, sec_index)[0]
    donor = _donor_mean_from_sections(sec[None, :, :], section_donor_idx, n_donors)[0]
    return _grand_mean(donor[None, :, :])[0]


def _endpoint_obs_and_boot(
    fov_obs: np.ndarray,
    sec_index: np.ndarray,
    n_sec: int,
    section_donor_idx: np.ndarray,
    n_donors: int,
    n_boot: int,
    rng: np.random.Generator,
):
    section_obs = _nanmean_fovs(fov_obs, sec_index, n_sec)
    donor_obs = _donor_mean_from_sections(section_obs[None, :, :], section_donor_idx, n_donors)[0]
    T = _grand_mean(donor_obs[None, :, :])[0]
    boot_T, boot_donor = _block_bootstrap_T(
        fov_obs, sec_index, n_sec, section_donor_idx, n_donors, n_boot, rng
    )
    donor_ci = np.stack(
        [_finite_quantile(boot_donor, 0.025, axis=0), _finite_quantile(boot_donor, 0.975, axis=0)],
        axis=0,
    )
    return T, donor_obs, section_obs, donor_ci, boot_T


def _observed_section_pool(
    preps: list[_Prep],
    sec_index: np.ndarray,
    n_sec: int,
    section_donor_idx: np.ndarray,
    n_donors: int,
    names: tuple[str, ...],
    min_arm: int,
):
    k = len(names)
    sec = np.full((n_sec, k), np.nan)
    high = np.full((n_sec, k), np.nan)
    low = np.full((n_sec, k), np.nan)
    by_sec: dict[int, list[_Prep]] = {}
    for i, p in enumerate(preps):
        by_sec.setdefault(int(sec_index[i]), []).append(p)
    for si, group in by_sec.items():
        for j in range(k):
            valid = [p for p in group if p.elig_sec[j]]
            if not valid:
                continue
            y = np.concatenate([p.Y[:, j] for p in valid])
            m = np.concatenate([p.mask_sec[:, j] for p in valid])
            if int(m.sum()) < min_arm or int((~m).sum()) < min_arm:
                continue
            high[si, j] = float(y[m].mean())
            low[si, j] = float(y[~m].mean())
            sec[si, j] = high[si, j] - low[si, j]
    donor = _donor_mean_from_sections(sec[None, :, :], section_donor_idx, n_donors)[0]
    T = _grand_mean(donor[None, :, :])[0]
    # Donor-level mean high/low via the same section means.
    donor_high = _donor_mean_from_sections(high[None, :, :], section_donor_idx, n_donors)[0]
    donor_low = _donor_mean_from_sections(low[None, :, :], section_donor_idx, n_donors)[0]
    return T, donor, sec, _grand_mean(donor_high[None, :, :])[0], _grand_mean(donor_low[None, :, :])[0]


def _observed_section_spearman(
    preps, sec_index, n_sec, section_donor_idx, n_donors, names
):
    k = len(names)
    sec = np.full((n_sec, k), np.nan)
    by_sec: dict[int, list[_Prep]] = {}
    for i, p in enumerate(preps):
        by_sec.setdefault(int(sec_index[i]), []).append(p)
    for si, group in by_sec.items():
        for j, name in enumerate(names):
            valid = [p for p in group if p.elig_sec[j]]
            if not valid:
                continue
            x = np.concatenate([p.cldn4 for p in valid])
            y = np.concatenate([p.Y[:, j] for p in valid])
            sec[si, j] = within_fov_spearman(x, y)
    donor = _donor_mean_from_sections(sec[None, :, :], section_donor_idx, n_donors)[0]
    T = _grand_mean(donor[None, :, :])[0]
    return T, donor, sec


def _section_pool_bootstrap(
    preps, sec_index, n_sec, section_donor_idx, n_donors, n_boot, min_arm, rng
):
    k = preps[0].Y.shape[1]
    by_sec: dict[int, list[_Prep]] = {}
    for i, p in enumerate(preps):
        by_sec.setdefault(int(sec_index[i]), []).append(p)
    boot_sec = np.full((n_boot, n_sec, k), np.nan)
    for si, group in by_sec.items():
        m = len(group)
        if m == 0:
            continue
        draws = rng.integers(0, m, size=(n_boot, m))
        for b in range(n_boot):
            chosen = [group[int(j)] for j in draws[b]]
            # A FOV with a non-finite outcome is omitted for that outcome only.
            for j in range(k):
                valid = [p for p in chosen if p.elig_sec[j]]
                if not valid:
                    continue
                x = np.concatenate([p.cldn4 for p in valid])
                y = np.concatenate([p.Y[:, j] for p in valid])
                med = float(np.median(x))
                hi = x > med
                n_h = int(hi.sum())
                n_l = int(hi.size - n_h)
                if n_h < min_arm or n_l < min_arm:
                    continue
                boot_sec[b, si, j] = float(y[hi].mean() - y[~hi].mean())
    donor = _donor_mean_from_sections(boot_sec, section_donor_idx, n_donors)
    return _grand_mean(donor), donor


def _donor_binned_curves(preps, donors, samples, sec_index, section_donor_idx, n_bins: int):
    k = preps[0].Y.shape[1]
    n_sec = len(samples)
    sec_curves = np.full((n_sec, n_bins, k), np.nan)
    by_sec: dict[int, list[_Prep]] = {}
    for i, p in enumerate(preps):
        by_sec.setdefault(int(sec_index[i]), []).append(p)
    for si, group in by_sec.items():
        acc = np.zeros((n_bins, k), dtype=np.float64)
        cnt = np.zeros((n_bins, k), dtype=np.float64)
        for p in group:
            for j in range(k):
                if not np.isfinite(p.Y[:, j]).all():
                    continue
                means = binned_means(p.cldn4, p.Y[:, j], n_bins=n_bins)
                ok = np.isfinite(means)
                acc[ok, j] += means[ok]
                cnt[ok, j] += 1
        with np.errstate(invalid="ignore", divide="ignore"):
            sec_curves[si] = acc / cnt
        sec_curves[si][cnt == 0] = np.nan
    # Average sections within donor, then we also return the grand mean.
    # sec_curves (n_sec, n_bins, k) -> donor (n_donors, n_bins, k)
    n_donors = len(donors)
    # Reuse donor mean by looping bins.
    donor_curves = np.full((n_donors, n_bins, k), np.nan)
    for b in range(n_bins):
        donor_curves[:, b, :] = _donor_mean_from_sections(
            sec_curves[:, b, :][None, :, :], section_donor_idx, n_donors
        )[0]
    grand = _grand_mean(donor_curves.transpose(1, 0, 2))  # (n_bins, k) if we reshape
    # donor_curves as (n_bins, n_donors, k) for _grand_mean which expects (n_rep, n_donors, k)
    grand = _grand_mean(np.transpose(donor_curves, (1, 0, 2)))
    return {"donor": donor_curves, "grand": grand, "n_bins": n_bins}


def beats_sample_sign_floor(p_one: float) -> bool:
    return bool(np.isfinite(p_one) and p_one < SAMPLE_SIGN_FLOOR_ONE)


def beats_donor_sign_floor(p_one: float) -> bool:
    return bool(np.isfinite(p_one) and p_one < DONOR_SIGN_FLOOR_ONE)
