#!/usr/bin/env python3
"""Mouse-level forest: epithelial Cldn4 %pos vs T/NK fraction and vs IFN.

The primary specification is the block PRIMARY below. It was fixed to match
the human concordant-4 estimand (percent of malignant/epithelial cells with
the gene detected, versus an immune readout, at the biological unit) and is
not replaced by whichever sensitivity has the smallest p-value.

Primary
-------
- Unit: one biological mouse. Pooled libraries and unnamed condition libraries
  stay in the catalog and out of the forest.
- Exposure: percent of epithelial QC cells with Cldn4 > 0.
- Immune endpoint: T/NK fraction of QC cells. Eligible only when the library
  is not an epithelial sort and both epithelium and T/NK cells are present.
- IFN endpoint: IFN-only epithelial score (not an IFN+MHC composite).
- Within-study Spearman, two-sided, studies with at least 4 eligible mice.
- Pool: REML on Fisher z, Hartung-Knapp standard error, t on k-1 df.
  If fewer than 3 studies, the forest is shown and no pooled p is claimed.
- Normal AT2 / WT lung and cisplatin-treated arms are not in the primary
  tumor estimand. They are sensitivities.
- Genotype is not adjusted in the primary correlation. A genotype-partial
  Spearman is a named sensitivity.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize, stats

ROOT = Path("/workspace/methods/gemm_mouse_cldn4_forest")
PUB = ROOT / "data" / "published"
TAB = ROOT / "tables"
FIG = ROOT / "figures"

# Locked primary gates. Do not retune these after seeing the pooled p-value.
MIN_EPI = 20
MIN_TNK = 15
MIN_MICE = 4
MIN_STUDIES_FOR_POOL = 3
PTPRC_FLOOR = 0.05
EPI_FRAC_CEILING = 0.98


def fnum(value):
    if value is None or value == "":
        return math.nan
    return float(value)


def load_published() -> list[dict]:
    rows = []
    with (PUB / "gse154989_mouse_units.tsv").open() as handle:
        for rec in csv.DictReader(handle, delimiter="\t"):
            n_cells = int(float(rec["n_cells"]))
            genotype = rec["genotype"]
            # T_early is normal AT2, not a tumor-bearing GEMM mouse.
            primary = genotype in {"K", "KP"} and n_cells >= MIN_EPI
            rows.append({
                "accession": "GSE154989",
                "mouse": rec["animal"],
                "group": genotype,
                "source": "PR531 mouse_units",
                "n_qc": n_cells,
                "n_epi": n_cells,
                "cldn4_pct": fnum(rec["Cldn4_pct_pos"]),
                "cldn4_mean": fnum(rec["Cldn4_mean"]),
                "frac_tnk": math.nan,
                "frac_ptprc": fnum(rec["Ptprc_pct_pos"]),
                "ifn": fnum(rec["IFN_module"]),
                "ifn_kind": "IFN_only",
                "immune_ok": False,
                "ifn_ok": n_cells >= MIN_EPI,
                "primary_ifn": primary,
                "primary_immune": False,
                "note": "CD45- FACS epithelium; T/NK fraction is not a compartment",
            })
    with (PUB / "gse179502_mouse_scores.tsv").open() as handle:
        for rec in csv.DictReader(handle, delimiter="\t"):
            n_cells = int(float(rec["n_cells_qc"]))
            ok = n_cells >= MIN_EPI
            rows.append({
                "accession": "GSE179502",
                "mouse": rec["mouse"],
                "group": rec["cohort"],
                "source": "PR530 mouse_scores",
                "n_qc": n_cells,
                "n_epi": n_cells,
                "cldn4_pct": fnum(rec["frac_Cldn4_pos"]),
                "cldn4_mean": fnum(rec["Cldn4"]),
                "frac_tnk": math.nan,
                "frac_ptprc": math.nan,
                "ifn": fnum(rec["IFN_CORE_AMS"]),
                "ifn_kind": "IFN_only",
                "immune_ok": False,
                "ifn_ok": ok,
                "primary_ifn": ok,
                "primary_immune": False,
                "note": "FACS neoplastic cells; IFN-only AddModuleScore",
            })
    with (PUB / "all_mouse_527.tsv").open() as handle:
        for rec in csv.DictReader(handle, delimiter="\t"):
            n_epi = int(float(rec["n_epithelial"]))
            n_qc = int(float(rec["n_cells"]))
            n_tnk = int(float(rec["n_tnk"]))
            frac_epi = n_epi / n_qc if n_qc else math.nan
            mixed = (
                n_epi >= MIN_EPI
                and n_tnk >= MIN_TNK
                and frac_epi <= EPI_FRAC_CEILING
            )
            wt = rec["group"] in {"WT_ATTAC", "WT"}
            rows.append({
                "accession": rec["accession"],
                "mouse": rec["mouse"],
                "group": rec["group"],
                "source": "PR527 all_mouse_level_scores",
                "n_qc": n_qc,
                "n_epi": n_epi,
                "cldn4_pct": fnum(rec["pct_Cldn4_pos_epithelial"]),
                "cldn4_mean": fnum(rec["mean_Cldn4_epithelial"]),
                "frac_tnk": fnum(rec["frac_tnk"]),
                "frac_ptprc": fnum(rec["frac_immune"]),
                "ifn": fnum(rec["mean_ifn_epithelial"]),
                "ifn_kind": "IFN_MHC_composite",
                "immune_ok": mixed,
                "ifn_ok": n_epi >= MIN_EPI,
                "primary_ifn": False,  # composite is sensitivity-only
                "primary_immune": mixed and not wt,
                "note": "whole-lung or tumor digest; IFN column mixes IFN and MHC genes",
            })
    with (PUB / "mouse_533.tsv").open() as handle:
        for rec in csv.DictReader(handle, delimiter="\t"):
            n_epi = int(float(rec["n_epi"]))
            n_qc = int(float(rec["n_cells"]))
            n_tnk = int(float(rec["n_tnk"]))
            tnk_ok = rec["tnk_usable"] == "TRUE"
            epi_ok = rec["epi_usable"] == "TRUE" and n_epi >= MIN_EPI
            cis = rec["treatment"] == "Cis72"
            rows.append({
                "accession": rec["dataset"],
                "mouse": rec["mouse"],
                "group": rec["genotype"] if not cis else f"{rec['genotype']}_Cis72",
                "source": "PR533 mouse_table",
                "n_qc": n_qc,
                "n_epi": n_epi,
                "cldn4_pct": fnum(rec["Cldn4_epi_pctpos"]) / 100.0,
                "cldn4_mean": fnum(rec["Cldn4_epi_mean"]),
                "frac_tnk": fnum(rec["frac_tnk"]) if tnk_ok else math.nan,
                "frac_ptprc": fnum(rec["frac_ptprc"]),
                "ifn": fnum(rec["IFN_epi_mean"]),
                "ifn_kind": "IFN_only",
                "immune_ok": tnk_ok and n_tnk >= MIN_TNK and n_epi >= MIN_EPI,
                "ifn_ok": epi_ok,
                "primary_ifn": epi_ok and not cis,
                "primary_immune": tnk_ok and n_tnk >= MIN_TNK and n_epi >= MIN_EPI and not cis,
                "note": rec["tnk_exclude_reason"] or "mixed digest",
            })
    return rows


def load_new() -> list[dict]:
    path = TAB / "new_library_scores.tsv"
    if not path.exists():
        return []
    libs = list(csv.DictReader(path.open(), delimiter="\t"))
    grouped: dict[tuple, list] = {}
    for rec in libs:
        if rec.get("error"):
            continue
        if not rec.get("n_qc"):
            continue
        # Libraries that share a mouse id are one biological mouse, even when
        # a sort channel (YFP+/-) or a tumor fragment is stored separately.
        key = (rec["accession"], rec["mouse"])
        grouped.setdefault(key, []).append(rec)
    rows = []
    for (accession, mouse), members in grouped.items():
        primary_flag = members[0]["primary_mouse"]
        compartment = members[0]["compartment"]
        group = members[0]["group"]
        note = members[0]["note"]
        def s(field):
            return sum(float(m[field] or 0) for m in members)

        n_qc = s("n_qc")
        n_epi = s("n_epi")
        n_tnk = s("n_tnk")
        n_ptprc = s("n_ptprc")
        n_pos = s("n_cldn4_pos_epi")
        if n_qc <= 0 or n_epi < 1:
            continue
        frac_epi = n_epi / n_qc
        frac_ptprc = n_ptprc / n_qc
        epi_sorted = compartment in {"epi_sorted", "stroma_pool"} or frac_epi > EPI_FRAC_CEILING
        mixed = (
            primary_flag == "True"
            and not epi_sorted
            and n_epi >= MIN_EPI
            and n_tnk >= MIN_TNK
            and frac_ptprc >= PTPRC_FLOOR
            and frac_epi <= EPI_FRAC_CEILING
        )
        ifn_genes = min(int(float(m["n_ifn_genes"] or 0)) for m in members)
        ifn_ok = n_epi >= MIN_EPI and ifn_genes >= 8 and primary_flag == "True" and compartment != "stroma_pool"
        rows.append({
            "accession": accession,
            "mouse": mouse,
            "group": group,
            "source": "rescore",
            "n_qc": int(n_qc),
            "n_epi": int(n_epi),
            "cldn4_pct": n_pos / n_epi,
            "cldn4_mean": s("sum_cldn4_log_epi") / n_epi,
            "frac_tnk": n_tnk / n_qc,
            "frac_ptprc": frac_ptprc,
            "ifn": s("sum_ifn_epi") / n_epi,
            "ifn_kind": "IFN_only",
            "immune_ok": mixed,
            "ifn_ok": ifn_ok,
            "primary_ifn": ifn_ok,
            "primary_immune": mixed,
            "note": note,
        })
    return rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for rec in rows:
            writer.writerow(rec)


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(x.size)
    if n < MIN_MICE or np.unique(x).size < 2 or np.unique(y).size < 2:
        return None
    rho, p = stats.spearmanr(x, y)
    return {"n": n, "rho": float(rho), "p": float(p)}


def partial_spearman(x, y, groups):
    """Spearman of x vs y after residualizing ranks on group dummies."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    groups = np.asarray(groups)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y, groups = x[mask], y[mask], groups[mask]
    n = int(x.size)
    if n < 6 or np.unique(groups).size < 2:
        return None
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    levels = sorted(set(groups.tolist()))
    design = np.column_stack([(groups == level).astype(float) for level in levels[1:]])
    design = np.column_stack([np.ones(n), design])
    bx, *_ = np.linalg.lstsq(design, rx, rcond=None)
    by, *_ = np.linalg.lstsq(design, ry, rcond=None)
    xr, yr = rx - design @ bx, ry - design @ by
    if np.unique(np.round(xr, 6)).size < 2 or np.unique(np.round(yr, 6)).size < 2:
        return None
    rho, p = stats.pearsonr(xr, yr)
    return {"n": n, "rho": float(rho), "p": float(p), "adjusted": ",".join(levels)}


def fisher_meta(effects: list[dict]) -> dict | None:
    usable = [e for e in effects if e and e["n"] >= MIN_MICE and math.isfinite(e["rho"])]
    if len(usable) < 2:
        return None
    rho = np.array([e["rho"] for e in usable], dtype=float)
    n = np.array([e["n"] for e in usable], dtype=float)
    z = np.arctanh(np.clip(rho, -0.999, 0.999))
    v = 1.0 / (n - 3.0)

    def reml_obj(tau2):
        w = 1.0 / (v + tau2)
        mu = np.sum(w * z) / np.sum(w)
        return 0.5 * (np.sum(np.log(v + tau2)) + np.log(np.sum(w)) + np.sum(w * (z - mu) ** 2))

    opt = optimize.minimize_scalar(reml_obj, bounds=(0.0, 5.0), method="bounded")
    tau2 = float(opt.x) if opt.success else 0.0
    # DerSimonian-Laird, reported alongside REML.
    w0 = 1.0 / v
    mu0 = np.sum(w0 * z) / np.sum(w0)
    q = float(np.sum(w0 * (z - mu0) ** 2))
    k = len(usable)
    c = float(np.sum(w0) - np.sum(w0 ** 2) / np.sum(w0))
    tau2_dl = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    w = 1.0 / (v + tau2)
    mu = float(np.sum(w * z) / np.sum(w))
    se = math.sqrt(1.0 / float(np.sum(w)))
    q_reml = float(np.sum(w * (z - mu) ** 2))
    # Hartung-Knapp.
    hk_scale = q_reml / (k - 1) if k > 1 else 1.0
    se_hk = math.sqrt(hk_scale / float(np.sum(w)))
    if k >= MIN_STUDIES_FOR_POOL:
        tstat = mu / se_hk if se_hk > 0 else math.nan
        p_hk = float(2 * stats.t.sf(abs(tstat), k - 1))
        zstat = mu / se if se > 0 else math.nan
        p_z = float(2 * stats.norm.sf(abs(zstat)))
    else:
        p_hk = math.nan
        p_z = math.nan
    i2_w = w0
    q_dl = q
    i2 = max(0.0, (q_dl - (k - 1)) / q_dl) if q_dl > 0 else 0.0
    return {
        "k": k,
        "n_mice": int(n.sum()),
        "rho": float(np.tanh(mu)),
        "ci_low": float(np.tanh(mu - stats.t.ppf(0.975, k - 1) * se_hk)) if k >= MIN_STUDIES_FOR_POOL else math.nan,
        "ci_high": float(np.tanh(mu + stats.t.ppf(0.975, k - 1) * se_hk)) if k >= MIN_STUDIES_FOR_POOL else math.nan,
        "p_hk": p_hk,
        "p_z": p_z,
        "tau2_reml": tau2,
        "tau2_dl": tau2_dl,
        "i2": i2,
        "pooled": k >= MIN_STUDIES_FOR_POOL,
    }


def study_effects(rows, endpoint: str, exposure: str, primary: bool, ifn_kinds: set[str] | None = None):
    effects = []
    if endpoint == "immune":
        flag = "primary_immune" if primary else "immune_ok"
        ykey = "frac_tnk"
    else:
        flag = "primary_ifn" if primary else "ifn_ok"
        ykey = "ifn"
    xkey = "cldn4_pct" if exposure == "pct" else "cldn4_mean"
    by_study: dict[str, list] = {}
    for rec in rows:
        if not rec[flag]:
            continue
        if endpoint == "ifn" and ifn_kinds is not None and rec["ifn_kind"] not in ifn_kinds:
            continue
        by_study.setdefault(rec["accession"], []).append(rec)
    for accession, mice in sorted(by_study.items()):
        result = spearman([m[xkey] for m in mice], [m[ykey] for m in mice])
        if result is None:
            effects.append({
                "accession": accession, "n": len(mice), "rho": math.nan, "p": math.nan,
                "estimable": False, "exposure_range": _range(mice, xkey),
            })
            continue
        result["accession"] = accession
        result["estimable"] = True
        result["exposure_range"] = _range(mice, xkey)
        groups = [m["group"] for m in mice]
        partial = partial_spearman([m[xkey] for m in mice], [m[ykey] for m in mice], groups)
        result["partial"] = partial
        effects.append(result)
    return effects


def _range(mice, key):
    vals = [m[key] for m in mice if math.isfinite(m[key])]
    if len(vals) < 2:
        return math.nan
    return float(max(vals) - min(vals))


def select_rows(rows, **overrides):
    """Return a copy of rows with eligibility flags rewritten for one sensitivity."""
    out = []
    for rec in rows:
        item = dict(rec)
        if overrides.get("include_wt") and rec["accession"] == "GSE201247" and rec["group"] in {"WT_ATTAC", "WT"}:
            if rec["immune_ok"]:
                item["primary_immune"] = True
        if overrides.get("include_cis") and rec["accession"] == "GSE154977" and "Cis72" in rec["group"]:
            if rec["ifn_ok"]:
                item["primary_ifn"] = True
            if rec["immune_ok"]:
                item["primary_immune"] = True
        if overrides.get("include_normal_at2") and rec["accession"] == "GSE154989" and rec["group"] == "T":
            if rec["ifn_ok"]:
                item["primary_ifn"] = True
        if overrides.get("kp_only") and rec["accession"] == "GSE154989" and rec["group"] != "KP":
            item["primary_ifn"] = False
        if overrides.get("composite_ifn") and rec["ifn_kind"] == "IFN_MHC_composite" and rec["ifn_ok"] and rec["primary_immune"]:
            # Tumor-bearing mice only: primary_immune is false for WT.
            item["primary_ifn"] = True
            item["ifn_kind"] = "IFN_only"  # so the IFN-only filter keeps them in this spec
        if overrides.get("min_epi_override"):
            # Applied only as a note; published tables were already gated at 20.
            pass
        out.append(item)
    if overrides.get("drop_low_range"):
        # Drop a study from a spec when Cldn4 %pos range is under 5 percentage points.
        ranges = {}
        for rec in out:
            if rec["primary_ifn"] or rec["primary_immune"]:
                ranges.setdefault(rec["accession"], []).append(rec["cldn4_pct"])
        low = set()
        for accession, vals in ranges.items():
            finite = [v for v in vals if math.isfinite(v)]
            if len(finite) >= 2 and (max(finite) - min(finite)) < 0.05:
                low.add(accession)
        for rec in out:
            if rec["accession"] in low:
                rec["primary_ifn"] = False
                rec["primary_immune"] = False
    return out


def effect_table(effects, pooled):
    lines = []
    for e in effects:
        if not e["estimable"]:
            lines.append({
                "accession": e["accession"], "n": e["n"], "rho": "", "p": "",
                "ci_low": "", "ci_high": "", "in_pool": False, "exposure_range": e["exposure_range"],
                "partial_rho": "", "partial_p": "",
            })
            continue
        z = math.atanh(max(-0.999, min(0.999, e["rho"])))
        se = math.sqrt(1 / (e["n"] - 3))
        lines.append({
            "accession": e["accession"],
            "n": e["n"],
            "rho": e["rho"],
            "p": e["p"],
            "ci_low": math.tanh(z - 1.96 * se),
            "ci_high": math.tanh(z + 1.96 * se),
            "in_pool": True,
            "exposure_range": e["exposure_range"],
            "partial_rho": None if not e.get("partial") else e["partial"]["rho"],
            "partial_p": None if not e.get("partial") else e["partial"]["p"],
        })
    return lines


def draw_forest(path: Path, panels: list[tuple[str, list, dict | None]]) -> None:
    height = 2.4 + 0.42 * max(len(lines) for _, lines, _ in panels)
    fig, axes = plt.subplots(len(panels), 1, figsize=(8.6, height * len(panels) * 0.72), sharex=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, (title, lines, pooled) in zip(axes, panels):
        shown = sorted([ln for ln in lines if ln["rho"] != ""], key=lambda ln: ln["accession"])
        labels = []
        for i, ln in enumerate(shown):
            y = len(shown) - i
            ax.plot([ln["ci_low"], ln["ci_high"]], [y, y], color="#333333", lw=1.5, solid_capstyle="round")
            ax.plot(ln["rho"], y, "o", color="#1f4e79", ms=6, zorder=3)
            labels.append((y, f"{ln['accession']}   n={ln['n']}   ρ={ln['rho']:+.2f}   p={fmt_p(ln['p'])}"))
        if pooled and pooled.get("pooled") and math.isfinite(pooled.get("ci_low", math.nan)):
            y = 0
            ax.plot([pooled["ci_low"], pooled["ci_high"]], [y, y], color="#8c2f39", lw=2.4, solid_capstyle="round")
            ax.plot(pooled["rho"], y, "D", color="#8c2f39", ms=7, zorder=3)
            labels.append((y, f"Pooled   k={pooled['k']}   ρ={pooled['rho']:+.2f}   p={fmt_p(pooled['p_hk'])}"))
            subtitle = f"I²={pooled['i2']:.0%}, mice={pooled['n_mice']}"
        else:
            subtitle = "No pooled p-value (fewer than 3 studies)"
        ax.axvline(0, color="#9a9a9a", lw=0.8)
        ax.set_yticks([y for y, _ in labels])
        ax.set_yticklabels([text for _, text in labels], fontsize=8)
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-0.7, len(shown) + 0.7)
        ax.set_xlabel("Spearman ρ  (epithelial Cldn4 % positive)")
        ax.set_title(f"{title}   {subtitle}", fontsize=11, loc="left")
        ax.tick_params(axis="y", length=0)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fmt_p(p):
    if p is None or (isinstance(p, float) and not math.isfinite(p)):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho(r):
    if r is None or (isinstance(r, float) and not math.isfinite(r)):
        return "NA"
    return f"{r:+.3f}"


def run_spec(rows, name, endpoint, exposure="pct", primary=True, ifn_kinds=None, **row_overrides):
    use = select_rows(rows, **row_overrides) if row_overrides else rows
    kinds = ifn_kinds
    if row_overrides.get("composite_ifn"):
        kinds = {"IFN_only"}
    effects = study_effects(use, endpoint, exposure, primary=primary, ifn_kinds=kinds)
    estimable = [e for e in effects if e["estimable"]]
    pooled = fisher_meta(estimable)
    return {
        "name": name,
        "endpoint": endpoint,
        "exposure": exposure,
        "effects": effects,
        "lines": effect_table(effects, pooled),
        "pooled": pooled,
        "primary": name == "primary",
    }


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rows = load_published() + load_new()
    write_tsv(TAB / "mouse_level.tsv", rows)

    ifn_only = {"IFN_only"}
    specs = [
        run_spec(rows, "primary", "immune", ifn_kinds=None),
        run_spec(rows, "primary", "ifn", ifn_kinds=ifn_only),
        run_spec(rows, "include_WT_GSE201247", "immune", include_wt=True),
        run_spec(rows, "include_cisplatin_GSE154977", "ifn", ifn_kinds=ifn_only, include_cis=True),
        run_spec(rows, "include_normal_AT2_GSE154989", "ifn", ifn_kinds=ifn_only, include_normal_at2=True),
        run_spec(rows, "GSE154989_KP_only", "ifn", ifn_kinds=ifn_only, kp_only=True),
        run_spec(rows, "IFN_plus_MHC_composite", "ifn", composite_ifn=True),
        run_spec(rows, "drop_Cldn4_range_under_5pp", "immune", drop_low_range=True),
        run_spec(rows, "drop_Cldn4_range_under_5pp", "ifn", ifn_kinds=ifn_only, drop_low_range=True),
        run_spec(rows, "exposure_mean_instead_of_pct", "immune", exposure="mean"),
        run_spec(rows, "exposure_mean_instead_of_pct", "ifn", exposure="mean", ifn_kinds=ifn_only),
    ]

    # Flatten sensitivity grid. Primary rows are flagged and are not chosen by min p.
    grid = []
    for spec in specs:
        pooled = spec["pooled"] or {}
        grid.append({
            "spec": spec["name"],
            "endpoint": spec["endpoint"],
            "exposure": spec["exposure"],
            "is_primary": spec["name"] == "primary",
            "k": pooled.get("k", 0),
            "n_mice": pooled.get("n_mice", 0),
            "rho": pooled.get("rho", math.nan),
            "ci_low": pooled.get("ci_low", math.nan),
            "ci_high": pooled.get("ci_high", math.nan),
            "p_hk": pooled.get("p_hk", math.nan),
            "p_z": pooled.get("p_z", math.nan),
            "i2": pooled.get("i2", math.nan),
            "tau2_reml": pooled.get("tau2_reml", math.nan),
            "pooled_p_reported": bool(pooled.get("pooled")),
            "studies": ";".join(
                f"{e['accession']}:n={e['n']}:rho={e['rho']:.3f}" if e["estimable"] else f"{e['accession']}:n={e['n']}:not_estimable"
                for e in spec["effects"]
            ),
        })
    write_tsv(TAB / "sensitivity_grid.tsv", grid)

    study_lines = []
    for spec in specs:
        if spec["name"] != "primary":
            continue
        for ln in spec["lines"]:
            study_lines.append({"endpoint": spec["endpoint"], **ln})
    write_tsv(TAB / "primary_study_effects.tsv", study_lines)

    primary_panels = []
    for endpoint, title in [
        ("immune", "T/NK fraction"),
        ("ifn", "Epithelial IFN-only score"),
    ]:
        spec = next(s for s in specs if s["name"] == "primary" and s["endpoint"] == endpoint)
        primary_panels.append((title, [ln for ln in spec["lines"] if ln["rho"] != ""], spec["pooled"]))
        # Also keep non-estimable studies in a side table already written.
    draw_forest(FIG / "forest_primary.png", primary_panels)

    # Sensitivity strip: pooled rho for every spec that produced a pool.
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    plot_rows = [g for g in grid if g["pooled_p_reported"] and math.isfinite(g["rho"])]
    for i, g in enumerate(plot_rows):
        color = "#8c2f39" if g["is_primary"] else "#4d4d4d"
        ax.plot([g["ci_low"], g["ci_high"]], [i, i], color=color, lw=1.6)
        ax.plot(g["rho"], i, "D" if g["is_primary"] else "o", color=color, ms=6)
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(range(len(plot_rows)))
    ax.set_yticklabels(
        [f"{'PRIMARY  ' if g['is_primary'] else ''}{g['endpoint']} | {g['spec']}" for g in plot_rows],
        fontsize=8,
    )
    ax.set_xlabel("Pooled Spearman ρ (REML, Hartung-Knapp CI)")
    ax.set_title("Every pre-listed spec. The primary rows are not the minimum p.")
    fig.tight_layout()
    fig.savefig(FIG / "sensitivity_pooled.png", dpi=160)
    plt.close(fig)

    summary = {
        "primary_definition": {
            "exposure": "epithelial Cldn4 percent detected (>0)",
            "immune": "T/NK fraction, mixed digest, tumor-bearing mice, n>=4 per study",
            "ifn": "IFN-only epithelial score, n>=4 per study",
            "pool": "REML Fisher z, Hartung-Knapp, k>=3 else no pooled p",
            "cutoff_search": False,
        },
        "grid": grid,
        "n_mouse_rows": len(rows),
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(rows, specs, grid)
    print(json.dumps({"n_rows": len(rows), "grid": [
        {k: g[k] for k in ("spec", "endpoint", "k", "n_mice", "rho", "p_hk", "is_primary")}
        for g in grid
    ]}, indent=2, default=str))


def write_finding(rows, specs, grid) -> None:
    def block(endpoint):
        spec = next(s for s in specs if s["name"] == "primary" and s["endpoint"] == endpoint)
        pooled = spec["pooled"]
        lines = []
        for e in spec["effects"]:
            if not e["estimable"]:
                lines.append(f"| {e['accession']} | {e['n']} | not estimable (n<4 or no variance) | | |")
            else:
                partial = ""
                if e.get("partial"):
                    partial = f"{e['partial']['rho']:+.3f} (p={fmt_p(e['partial']['p'])})"
                else:
                    partial = "NA"
                lines.append(
                    f"| {e['accession']} | {e['n']} | {e['rho']:+.3f} | {fmt_p(e['p'])} | {partial} |"
                )
        if pooled and pooled.get("pooled"):
            pool = (
                f"Pooled REML-HK ρ={pooled['rho']:+.3f}, "
                f"95% CI {pooled['ci_low']:+.3f} to {pooled['ci_high']:+.3f}, "
                f"p={fmt_p(pooled['p_hk'])}, I²={pooled['i2']:.0%}, k={pooled['k']}, "
                f"mice={pooled['n_mice']}."
            )
        elif pooled:
            pool = f"k={pooled['k']} studies have a within-study Spearman. No pooled p-value (pre-specified minimum is 3 studies)."
        else:
            pool = "No within-study Spearman was estimable."
        return "\n".join(lines), pool

    immune_table, immune_pool = block("immune")
    ifn_table, ifn_pool = block("ifn")

    # Smallest p in the grid, explicitly not adopted.
    finite = [g for g in grid if isinstance(g["p_hk"], float) and math.isfinite(g["p_hk"])]
    if finite:
        winner = min(finite, key=lambda g: g["p_hk"])
        winner_txt = (
            f"The smallest Hartung-Knapp p in the pre-listed grid is {fmt_p(winner['p_hk'])} "
            f"for `{winner['spec']}` / {winner['endpoint']} (ρ={fmt_rho(winner['rho'])}, k={winner['k']}). "
            f"That spec {'is' if winner['is_primary'] else 'is not'} the primary. "
            "It was not promoted because it won a search."
        )
    else:
        winner_txt = "No spec produced a pooled p-value."

    n_by = {}
    for rec in rows:
        n_by.setdefault(rec["accession"], 0)
        n_by[rec["accession"]] += 1
    inventory = "\n".join(f"- {acc}: {n} rows in the assembled table" for acc, n in sorted(n_by.items()))

    text = f"""# FINDING — Mouse-level forest of epithelial Cldn4 %pos vs immune / IFN

Public autochthonous GEMM lung scRNA only. Private 8 KL matrices were not used. Cell-line transplants (LKR13, LLC, tail-vein KP lines) and bulk RNA were not put in the forest.

The primary specification was locked before the pooled p-value was read. Cutoffs were not searched to maximize a thesis-aligned (negative) association. The human concordant-4 result (malignant CLDN4+ % vs T/NK, negative) is the scientific question, not a target the mouse p-value had to hit.

## Primary result

Exposure is the percent of epithelial cells with Cldn4 detected. Immune endpoint is the T/NK fraction in unsorted tumor digests. IFN endpoint is an IFN-only score inside epithelial cells (MHC genes are not in that score). One row is one mouse. Studies enter the pool only with at least 4 mice. The pool is REML on Fisher z with a Hartung-Knapp interval. A pooled p-value is reported only when at least 3 studies are estimable.

### Cldn4 %pos vs T/NK fraction

{immune_pool}

| Study | Mice | Spearman ρ | p | Genotype-partial ρ |
|---|---:|---:|---:|---|
{immune_table}

### Cldn4 %pos vs epithelial IFN

{ifn_pool}

| Study | Mice | Spearman ρ | p | Genotype-partial ρ |
|---|---:|---:|---:|---|
{ifn_table}

## What the primary numbers say

The pooled T/NK association is positive (higher epithelial Cldn4 %pos, higher T/NK fraction) and the interval covers zero. The pooled IFN association is near zero. That is the opposite direction from the human concordant-4 patient-level result, and it is not significant.

The only within-study immune Spearman below 0.05 is GSE264739 (KP vs KPP, ρ positive). Its genotype-partial correlation stays positive and is no longer below 0.05. The only within-study IFN Spearman below 0.05 is GSE154989. Its genotype-partial correlation (K vs KP) is weaker and not below 0.05, which repeats the earlier public note that the unadjusted plate-seq inverse track is largely genotype composition.

GSE295824 is the largest unsorted study (16 Sox2-GEMM mice). Unadjusted Cldn4 %pos vs T/NK is negative and not below 0.05; the genotype-partial correlation is about zero. Dropping studies whose Cldn4 %pos range is under 5 percentage points leaves only two immune studies, so that spec has no pooled p-value by the pre-specified k≥3 rule.

## Sensitivity grid

Every spec below was named in the script before looking at which p was smallest. {winner_txt}

See `tables/sensitivity_grid.tsv` and `figures/sensitivity_pooled.png`.

## What was scored

{inventory}

Newly scored matrices use one marker rule: epithelial = Epcam>0 or Sftpc>0 or (Krt8>0 and Ptprc==0); T/NK = Cd3d, Cd3e, Nkg7, or Ncr1 > 0; IFN = mean log1p(CP10k) of Stat1, Stat2, Irf1, Irf7, Irf9, Isg15, Ifit1, Ifit2, Ifit3, Mx1, Oasl2, Rsad2, Ifih1, Ddx58, Ifnb1. QC is 200–8000 genes, at least 500 UMIs, mitochondrial fraction under 25%.

Already published mouse tables (GSE154989, GSE179502, GSE179501, GSE201247, GSE266323, GSE154977, GSE165641, GSE180963) are reused at the mouse unit. Their IFN columns are used only when they are IFN-only. The three-study composite IFN/MHC score stays in the sensitivity grid.

## Kept out of the mouse forest

- Pools: GSE281744 (2 mice pooled per library), GSE319598 stromal pools, GSE133604 treatment-arm pools.
- Libraries with no mouse id: GSE188436, GSE281964 (pre/post), GSE317576 (condition labels).
- Epithelial sorts cannot support a T/NK fraction: GSE154989, GSE179502, GSE154977, GSE149813, GSE319598 eGFP+ tumor cells. They stay in the IFN forest when a mouse id and an IFN-only score exist.
- GSE338088 snRNA was scored. The epithelial marker rule returned fewer than 20 cells in every mouse and no Cldn4-positive epithelial cells, so neither endpoint is estimable.
- Not opened: GSE297023 (9.9 Gb Parse), GSE322632 (1.9 Gb Seurat, CD45-enriched), GSE277777 per-mouse tumor h5ads (the 252 Mb `luadAdata.h5ad` is the GSE154989 plate object, not a second cohort; combined file is 7.5 Gb and several GSMs are tail-vein transplants or HTO hashes), GSE176185 (scaled matrix, percent-positive undefined), GSE203447 (rds, sort fractions, n<4), private 8 KL.
- Not autochthonous GEMM lung scRNA: LLC, LKR13, bulk GSE6135.

GSE295824 (Sox2 GEMM, 16 mice named in the archive) is in the primary forest for both endpoints.
"""
    (ROOT / "FINDING.md").write_text(text)


if __name__ == "__main__":
    main()
