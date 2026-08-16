#!/usr/bin/env python3
"""Build WRITEUP.md from committed tables/JSON. No hand-edited numbers."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (pd.isna(p))):
        return "—"
    p = float(p)
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or (isinstance(r, float) and pd.isna(r)):
        return "—"
    return f"{float(r):+.3f}"


def main() -> None:
    out = Path(__file__).resolve().parent
    audit = json.loads((out / "dataset_audit.json").read_text())
    summary = json.loads((out / "summary.json").read_text())
    stats = pd.read_csv(out / "tables" / "per_dataset_spearman.tsv", sep="\t")
    pool = pd.read_csv(out / "tables" / "pooled_spearman.tsv", sep="\t")
    units = pd.read_csv(out / "tables" / "per_unit_metrics.tsv", sep="\t")

    lines = []
    lines.append("# TISCH2 NSCLC pool — epithelial TACSTD2/CLDN4 vs T/NK fraction")
    lines.append("")
    lines.append("**Slice:** `methods/tisch_nsclc_pool/` only. Public TISCH2 h5 + CellMetainfo.")
    lines.append("Numbers below are written from `tables/` and `summary.json`.")
    lines.append("")
    lines.append("## Verdict (honest)")
    lines.append("")

    def pool_row(contrast, method="fisher_z_ivw"):
        hit = pool.loc[(pool["contrast"] == contrast) & (pool["method"] == method)]
        return hit.iloc[0] if len(hit) else None

    t2 = pool_row("TACSTD2_mean_vs_fracTNK")
    c4 = pool_row("CLDN4_mean_vs_fracTNK")
    if t2 is not None:
        lines.append(
            f"- **TACSTD2 mean vs T/NK (Fisher-z pool):** "
            f"n_datasets={int(t2['n_datasets'])}, n_units={int(t2['n_units_sum'])}, "
            f"ρ={fmt_rho(t2['rho'])}, p={fmt_p(t2['p'])}, "
            f"95% CI [{fmt_rho(t2.get('ci95_lo'))}, {fmt_rho(t2.get('ci95_hi'))}], "
            f"I²={t2['I2']:.0f}%."
            if pd.notna(t2.get("rho", float("nan")))
            else f"- **TACSTD2 mean vs T/NK (Fisher-z pool):** {t2.get('note','')}"
        )
    if c4 is not None:
        lines.append(
            f"- **CLDN4 mean vs T/NK (Fisher-z pool):** "
            f"n_datasets={int(c4['n_datasets'])}, n_units={int(c4['n_units_sum'])}, "
            f"ρ={fmt_rho(c4['rho'])}, p={fmt_p(c4['p'])}, "
            f"95% CI [{fmt_rho(c4.get('ci95_lo'))}, {fmt_rho(c4.get('ci95_hi'))}], "
            f"I²={c4['I2']:.0f}%."
            if pd.notna(c4.get("rho", float("nan")))
            else f"- **CLDN4 mean vs T/NK (Fisher-z pool):** {c4.get('note','')}"
        )
    lines.append(
        "- **ICI labels:** TISCH gallery marks GSE151537, GSE146100, and GSE176021_aPD1 as Immunotherapy. "
        "Downloadable CellMetainfo has **no Response/RECIST/MPR column**. "
        "GSE151537 and GSE176021 are T-sorted (no epithelium). GSE146100 is 1 patient / 3 nodules. "
        "This pool is **not** an ICI-response test."
    )
    lines.append("")
    lines.append("## Per-dataset Spearman (eligible tumor-like units)")
    lines.append("")
    lines.append("| Dataset | Unit | Contrast | n | ρ | p | note |")
    lines.append("|---|---|---|---:|---:|---:|---|")
    if len(stats):
        show = stats.loc[stats["contrast"].isin(
            ["TACSTD2_mean_vs_fracTNK", "CLDN4_mean_vs_fracTNK"]
        )].sort_values(["dataset", "contrast"])
        for _, r in show.iterrows():
            lines.append(
                f"| {r['dataset']} | {r['unit']} | {r['contrast']} | "
                f"{int(r['n'])} | {fmt_rho(r['rho'])} | {fmt_p(r['p'])} | {r.get('note','')} |"
            )
    lines.append("")
    lines.append("## ICI / skip audit")
    lines.append("")
    lines.append("| Dataset | Gallery | CellMetainfo ICI cols | Eligible units | In pool | Skip / ICI note |")
    lines.append("|---|---|---|---:|---|---|")
    for d in audit["datasets"]:
        ici = ",".join(d.get("ici_metainfo_columns") or []) or "none"
        lines.append(
            f"| {d['dataset']} | {d.get('ici_gallery')} | {ici} | "
            f"{d.get('n_units_eligible', 0) if d.get('n_units_eligible') is not None else 0} | "
            f"{'yes' if d.get('included_in_pool') else 'no'} | "
            f"{d.get('skip_reason') or d.get('ici_label_honest') or ''} |"
        )
    lines.append("")
    lines.append("## Coverage")
    lines.append("")
    n_elig = int(units["eligible"].sum()) if len(units) and "eligible" in units.columns else 0
    lines.append(f"- Catalog datasets: {audit['n_datasets_catalog']}")
    lines.append(f"- Datasets in Fisher-z pool: {audit['n_included_in_pool']}")
    lines.append(f"- Eligible sample/patient units written: {n_elig}")
    lines.append(f"- Generated: {summary.get('generated_at')}")
    lines.append("")
    lines.append("## Caveats")
    lines.append("")
    lines.append("- TISCH values are MAESTRO `log2(TPM/10+1)`, not raw UMI.")
    lines.append("- GSE131907 has no TISCH Malignant call; tumor-tissue epithelial cells are a proxy.")
    lines.append("- T/NK fraction is composition after dissociation and TISCH annotation, not spatial infiltration.")
    lines.append("- Sample-level tests are the unit. Cell-level p-values are not reported as primary.")
    lines.append("- GSE146100 nodule-level n=3 is listed in `per_unit_metrics.tsv` and excluded from Spearman.")
    lines.append("")
    (out / "WRITEUP.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {out / 'WRITEUP.md'}")


if __name__ == "__main__":
    main()
