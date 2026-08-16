"""Score a prepared expression object with every method the playbook specifies.

A prepared object is a dict with:

    expr_log      genes x samples, log2(x+1) of a library-size-normalised matrix
    expr_linear   the same genes, linear (TPM/CPM/FPKM) -- required for CIBERSORT / CYT
    hgnc_symbols  already collapsed, unique, current symbols

Methods that cannot run (missing optional dependency, missing LM22 file) are
skipped and recorded in ``skipped`` rather than crashing the rest of the run.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from . import estimate as estimate_mod
from . import mcpcounter
from . import signatures as sigmod
from . import ssgsea as ssgsea_mod
from . import tip as tip_mod
from . import xcell as xcell_mod
from .cibersort import cibersort, constrained_ls_deconvolve, load_signature_matrix
from .genes import match_report

__all__ = ["score_all", "TARGETS"]

TARGETS = ("TACSTD2", "CLDN4")


def _read_list(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def score_all(
    expr_log: pd.DataFrame,
    expr_linear: pd.DataFrame,
    resource_dir: str | Path,
    run_xcell: bool = True,
    run_cibersort: bool = True,
    run_tide: bool = True,
    lm22_path: str | Path | None = None,
    permutations: int = 0,
) -> dict[str, Any]:
    resource_dir = Path(resource_dir)
    skipped: dict[str, str] = {}
    blocks: dict[str, pd.DataFrame] = {}
    reports: dict[str, pd.DataFrame] = {}

    # --- curated + hallmark ssGSEA ---------------------------------------
    lib = sigmod.SignatureLibrary.load(
        resource_dir / "curated_signatures.json",
        msigdb_gmt=resource_dir / "h.all.v2025.1.Hs.symbols.gmt"
        if (resource_dir / "h.all.v2025.1.Hs.symbols.gmt").exists()
        else resource_dir / "h.all.v2024.1.Hs.symbols.gmt"
        if (resource_dir / "h.all.v2024.1.Hs.symbols.gmt").exists()
        else None,
        msigdb_keep=[
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_INTERFERON_ALPHA_RESPONSE",
            "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_ALLOGRAFT_REJECTION",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "HALLMARK_TGF_BETA_SIGNALING",
            "HALLMARK_HYPOXIA",
            "HALLMARK_ANGIOGENESIS",
            "HALLMARK_IL6_JAK_STAT3_SIGNALING",
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        ],
    )
    reports["signature_overlap"] = match_report(expr_log.index, lib.sets)
    ssg = ssgsea_mod.ssgsea_scores(expr_log, lib.sets, min_size=2)
    blocks["ssgsea"] = ssg
    try:
        blocks["ssgsea"] = blocks["ssgsea"].assign(CYT=sigmod.cytolytic_activity(expr_linear))
    except Exception as exc:
        skipped["CYT"] = str(exc)

    # --- MCP-counter ------------------------------------------------------
    mcp_path = resource_dir / "mcp_counter_genes.txt"
    if mcp_path.exists():
        markers = mcpcounter.load_markers(mcp_path)
        reports["mcp_overlap"] = match_report(expr_log.index, markers)
        blocks["mcp"] = mcpcounter.mcp_counter(expr_log, markers)
    else:
        skipped["mcp"] = f"missing {mcp_path}"

    # --- ESTIMATE ---------------------------------------------------------
    est_gmt = resource_dir / "estimate_signatures.gmt"
    est_common = resource_dir / "estimate_common_genes.tsv"
    if est_gmt.exists():
        est_sets = estimate_mod.load_signatures(est_gmt)
        common = estimate_mod.load_common_genes(est_common) if est_common.exists() else None
        blocks["estimate"] = estimate_mod.estimate_score(expr_log, est_sets, common_genes=common)
    else:
        skipped["estimate"] = f"missing {est_gmt}"

    # --- xCell ------------------------------------------------------------
    if run_xcell:
        xgmt = resource_dir / "xCell_signatures.gmt"
        xgenes = resource_dir / "xCell_genes.txt"
        xk = resource_dir / "xCell_spillK_rnaseq.tsv"
        xfv = resource_dir / "xCell_fv_rnaseq.tsv"
        if all(p.exists() for p in (xgmt, xgenes, xk, xfv)):
            try:
                xscores = xcell_mod.xcell_analysis(
                    expr_log,
                    xcell_mod.load_signatures(xgmt),
                    _read_list(xgenes),
                    pd.read_csv(xk, sep="\t", index_col=0),
                    pd.read_csv(xfv, sep="\t", index_col=0),
                    cell_types_use=[
                        "B-cells", "CD4+ T-cells", "CD8+ T-cells", "Tregs",
                        "NK cells", "Monocytes", "Macrophages", "Macrophages M1",
                        "Macrophages M2", "DC", "cDC", "pDC", "Neutrophils",
                        "Mast cells", "Eosinophils", "Fibroblasts",
                        "Endothelial cells", "Adipocytes",
                    ],
                )
                blocks["xcell"] = xscores.T
            except Exception as exc:
                skipped["xcell"] = str(exc)
        else:
            skipped["xcell"] = "missing xCell resource files"

    # --- TIDE -------------------------------------------------------------
    if run_tide:
        try:
            from .tide import tide_score

            blocks["tide"] = tide_score(expr_linear, cancer="NSCLC")
        except Exception as exc:
            skipped["tide"] = str(exc)

    # --- TIP --------------------------------------------------------------
    tip_path = resource_dir / "tip_signature_annotation.txt"
    if tip_path.exists():
        ann = tip_mod.load_tip_annotation(tip_path)
        blocks["tip"] = tip_mod.tip_score(expr_log, ann)
    else:
        skipped["tip"] = f"missing {tip_path}"

    # --- CIBERSORT-style / quanTIseq-style --------------------------------
    if run_cibersort:
        sig_path = Path(lm22_path) if lm22_path else resource_dir / "LM22.txt"
        if not sig_path.exists():
            sig_path = resource_dir / "TIL10_signature.txt"
            method_name = "quantiseq_style_TIL10"
        else:
            method_name = "cibersort_LM22"
        if sig_path.exists():
            try:
                sig = load_signature_matrix(sig_path)
                if method_name.startswith("cibersort"):
                    res = cibersort(
                        expr_linear, sig, permutations=permutations, verbose=True
                    )
                    blocks["cibersort"] = res.fractions
                    blocks["cibersort_fit"] = pd.DataFrame(
                        {
                            "correlation": res.correlation,
                            "rmse": res.rmse,
                            "p_value": res.p_value,
                        }
                    )
                else:
                    scaling_path = resource_dir / "TIL10_mRNA_scaling.txt"
                    scaling = None
                    if scaling_path.exists():
                        scaling = pd.read_csv(
                            scaling_path, sep="\t", header=None, index_col=0
                        ).iloc[:, 0]
                    blocks["quantiseq_style"] = constrained_ls_deconvolve(
                        expr_linear, sig, scaling=scaling
                    )
            except Exception as exc:
                skipped[method_name] = str(exc)
        else:
            skipped["deconvolution"] = "neither LM22.txt nor TIL10_signature.txt present"

    # --- targets ----------------------------------------------------------
    present = [g for g in TARGETS if g in expr_log.index]
    targets = expr_log.loc[present].T if present else pd.DataFrame(index=expr_log.columns)
    if "TACSTD2" in expr_log.index and "CLDN4" in expr_log.index:
        targets["TACSTD2_CLDN4_mean"] = expr_log.loc[["TACSTD2", "CLDN4"]].mean(axis=0)

    # flatten score blocks to one samples x scores frame (prefixed)
    flat_parts = []
    for name, frame in blocks.items():
        if name.endswith("_fit"):
            continue
        prefixed = frame.copy()
        prefixed.columns = [f"{name}:{c}" for c in prefixed.columns]
        flat_parts.append(prefixed)
    flat = pd.concat(flat_parts, axis=1) if flat_parts else pd.DataFrame(index=expr_log.columns)

    return {
        "blocks": blocks,
        "scores": flat,
        "targets": targets,
        "skipped": skipped,
        "reports": reports,
        "library": lib,
    }
