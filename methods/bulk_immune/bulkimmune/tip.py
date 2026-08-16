"""TIP -- Tracking Tumor Immunophenotype (Xu et al., Cancer Research 2018).

TIP scores the seven-step cancer-immunity cycle from a published annotation
table (GeneSymbol, Steps, Direction, ImmuneCellType) downloaded from
http://biocc.hrbmu.edu.cn/TIP/download/signature%20annotation.txt.

The original web server computes a per-step ssGSEA of the *positive* genes
minus a per-step ssGSEA of the *negative* genes. That is what is implemented
here. Step 4 (immune-cell recruitment) is also reported per cell type because
that is the biologically interesting resolution for an exclusion analysis.

This is **not** a licensed reimplementation of the TIP web application's
private scoring weights; it is the published gene-set definition scored with
the same ssGSEA statistic the paper describes. Cite Xu 2018 and report the
overlap of each step with your matrix.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .ssgsea import ssgsea

__all__ = ["load_tip_annotation", "tip_sets", "tip_score", "TIP_STEPS"]

TIP_STEPS = {
    1: "Step1_release_of_cancer_antigens",
    2: "Step2_cancer_antigen_presentation",
    3: "Step3_priming_and_activation",
    4: "Step4_trafficking_of_immune_cells",
    5: "Step5_infiltration_of_immune_cells",
    6: "Step6_recognition_of_cancer_cells",
    7: "Step7_killing_of_cancer_cells",
}


def load_tip_annotation(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str)
    required = {"GeneSymbol", "Steps", "Direction"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"TIP annotation missing columns {missing}")
    table["Steps"] = table["Steps"].astype(int)
    table["Direction"] = table["Direction"].str.lower().str.strip()
    return table


def tip_sets(annotation: pd.DataFrame, by_cell_type: bool = True) -> dict[str, list[str]]:
    """Build the ssGSEA gene sets from the TIP annotation table."""
    sets: dict[str, list[str]] = {}
    for step, label in TIP_STEPS.items():
        block = annotation[annotation["Steps"] == step]
        for direction in ("positive", "negative"):
            genes = list(dict.fromkeys(block.loc[block["Direction"] == direction, "GeneSymbol"]))
            if genes:
                sets[f"{label}__{direction}"] = genes
        if by_cell_type and "ImmuneCellType" in block.columns and step == 4:
            for cell, sub in block.groupby("ImmuneCellType"):
                genes = list(dict.fromkeys(sub.loc[sub["Direction"] == "positive", "GeneSymbol"]))
                if genes:
                    sets[f"Step4_{cell.replace(' ', '_')}__positive"] = genes
    return sets


def tip_score(
    expr: pd.DataFrame,
    annotation: pd.DataFrame,
    alpha: float = 0.25,
    by_cell_type: bool = True,
) -> pd.DataFrame:
    """TIP step scores (samples x steps).

    A step score is ``ssGSEA(positive) - ssGSEA(negative)``. Steps with only
    one direction keep that direction's score. Cell-type-resolved Step 4
    scores are appended when ``by_cell_type`` is True.
    """
    sets = tip_sets(annotation, by_cell_type=by_cell_type)
    raw = ssgsea(expr, sets, alpha=alpha, normalize="none", tie_method="average_int")

    columns: dict[str, pd.Series] = {}
    for step, label in TIP_STEPS.items():
        pos_key = f"{label}__positive"
        neg_key = f"{label}__negative"
        pos = raw.loc[pos_key] if pos_key in raw.index else 0.0
        neg = raw.loc[neg_key] if neg_key in raw.index else 0.0
        columns[label] = pos - neg

    if by_cell_type:
        for name in raw.index:
            if name.startswith("Step4_") and name.endswith("__positive"):
                columns[name.replace("__positive", "")] = raw.loc[name]

    out = pd.DataFrame(columns)
    out.attrs["overlap"] = raw.attrs.get("overlap", {})
    return out
