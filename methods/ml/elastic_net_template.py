#!/usr/bin/env python3
"""Leakage-resistant nested-CV elastic-net example for binary ICI outcomes.

This is an educational template, not a clinically validated pipeline. It
deliberately refuses very small or severely imbalanced datasets. Passing the
guards is not evidence that the sample size is adequate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    log_loss,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# Conservative demonstration guards, not a formal sample-size calculation.
MIN_PATIENTS = 40
MIN_PER_CLASS = 15
OUTER_FOLDS = 5
INNER_FOLDS = 4
RANDOM_SEED = 20260816


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run nested-CV elastic-net logistic regression. All learned "
            "preprocessing and tuning occur inside the relevant training fold."
        )
    )
    parser.add_argument("csv", type=Path, help="Patient-level CSV input")
    parser.add_argument(
        "--outcome",
        required=True,
        help="Binary outcome column encoded as 0/1",
    )
    parser.add_argument(
        "--patient-id",
        required=True,
        help="Unique patient identifier column",
    )
    parser.add_argument(
        "--split-group",
        help=(
            "Optional study/site/batch column. If supplied, entire groups are "
            "held out together in both CV loops."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("nested_cv_output"),
        help="Directory for out-of-fold predictions and summary",
    )
    return parser.parse_args()


def validate_data(
    frame: pd.DataFrame,
    outcome_col: str,
    patient_col: str,
    group_col: str | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, list[str]]:
    required = [outcome_col, patient_col] + ([group_col] if group_col else [])
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if frame[patient_col].isna().any() or not frame[patient_col].is_unique:
        raise ValueError(
            "Each row must represent one patient and patient IDs must be unique. "
            "Resolve repeated measures before running this template."
        )

    if len(frame) < MIN_PATIENTS:
        raise ValueError(
            f"REFUSING TO RUN: n={len(frame)} is below the hard demonstration "
            f"minimum of {MIN_PATIENTS}. Do not train a signature on this cohort."
        )

    if frame[outcome_col].isna().any():
        raise ValueError("Outcome labels cannot be missing.")
    observed_labels = set(frame[outcome_col].unique())
    if not observed_labels.issubset({0, 1}) or len(observed_labels) != 2:
        raise ValueError(
            f"Outcome must contain both binary labels 0 and 1; got {observed_labels}."
        )

    y = frame[outcome_col].astype(int).to_numpy()
    counts = np.bincount(y, minlength=2)
    if counts.min() < MIN_PER_CLASS:
        raise ValueError(
            "REFUSING TO RUN: the minority class has "
            f"{counts.min()} patients; at least {MIN_PER_CLASS} are required by "
            "this demonstration guard."
        )

    excluded = {outcome_col, patient_col}
    if group_col:
        excluded.add(group_col)
    feature_cols = [column for column in frame.columns if column not in excluded]
    non_numeric = [
        column
        for column in feature_cols
        if not pd.api.types.is_numeric_dtype(frame[column])
    ]
    if non_numeric:
        raise ValueError(
            "All candidate features must be numeric. Encode only prespecified "
            f"variables before analysis; nonnumeric columns: {non_numeric}"
        )
    if not feature_cols:
        raise ValueError("No numeric candidate features remain.")

    X = frame[feature_cols].to_numpy(dtype=float)
    if np.isinf(X).any():
        raise ValueError("Features contain infinity; resolve this before analysis.")

    groups: np.ndarray | None = None
    if group_col:
        if frame[group_col].isna().any():
            raise ValueError("Split-group labels cannot be missing.")
        groups = frame[group_col].astype(str).to_numpy()
        group_count = np.unique(groups).size
        if group_count < OUTER_FOLDS:
            raise ValueError(
                f"Only {group_count} distinct split groups are available; "
                f"{OUTER_FOLDS} outer folds require at least that many groups."
            )

    return X, y, groups, feature_cols


def make_pipeline() -> Pipeline:
    return Pipeline(
        steps=[
            ("impute", SimpleImputer(strategy="median")),
            ("variance", VarianceThreshold()),
            ("scale", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    penalty="elasticnet",
                    solver="saga",
                    class_weight=None,
                    max_iter=20_000,
                    random_state=RANDOM_SEED,
                ),
            ),
        ]
    )


def make_splitter(folds: int, grouped: bool, seed: int):
    if grouped:
        return StratifiedGroupKFold(
            n_splits=folds,
            shuffle=True,
            random_state=seed,
        )
    return StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)


def run_nested_cv(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, float]]]:
    grouped = groups is not None
    outer = make_splitter(OUTER_FOLDS, grouped, RANDOM_SEED)
    outer_splits = outer.split(X, y, groups) if grouped else outer.split(X, y)

    probabilities = np.full(len(y), np.nan)
    fold_ids = np.full(len(y), -1, dtype=int)
    selected_params: list[dict[str, float]] = []
    parameter_grid = {
        "model__C": np.logspace(-3, 2, 10),
        "model__l1_ratio": [0.05, 0.25, 0.5, 0.75, 0.95],
    }

    for fold, (train_idx, test_idx) in enumerate(outer_splits):
        y_train, y_test = y[train_idx], y[test_idx]
        if np.unique(y_test).size != 2:
            raise ValueError(
                f"Outer fold {fold} lacks one outcome class. Use fewer folds or "
                "redesign the grouped split; do not score this split."
            )

        inner = make_splitter(INNER_FOLDS, grouped, RANDOM_SEED + fold + 1)
        search = GridSearchCV(
            estimator=make_pipeline(),
            param_grid=parameter_grid,
            scoring="neg_log_loss",
            cv=inner,
            n_jobs=-1,
            refit=True,
            error_score="raise",
        )
        if grouped:
            search.fit(X[train_idx], y_train, groups=groups[train_idx])
        else:
            search.fit(X[train_idx], y_train)

        probabilities[test_idx] = search.predict_proba(X[test_idx])[:, 1]
        fold_ids[test_idx] = fold
        selected_params.append(
            {
                "fold": fold,
                "C": float(search.best_params_["model__C"]),
                "l1_ratio": float(search.best_params_["model__l1_ratio"]),
            }
        )

    if np.isnan(probabilities).any() or (fold_ids < 0).any():
        raise RuntimeError("Not every patient received exactly one outer-fold prediction.")
    return probabilities, fold_ids, selected_params


def main() -> None:
    args = parse_args()
    frame = pd.read_csv(args.csv)
    X, y, groups, feature_cols = validate_data(
        frame,
        args.outcome,
        args.patient_id,
        args.split_group,
    )
    probabilities, fold_ids, selected_params = run_nested_cv(X, y, groups)

    ledger_columns = [args.patient_id, args.outcome]
    if args.split_group:
        ledger_columns.append(args.split_group)
    ledger = frame[ledger_columns].copy()
    ledger["outer_fold"] = fold_ids
    ledger["oof_probability"] = probabilities

    summary = {
        "warning": (
            "Exploratory nested-CV results; not external validation and not a "
            "clinically deployable model. Probabilities are not post-hoc calibrated."
        ),
        "n_patients": int(len(y)),
        "n_features": int(len(feature_cols)),
        "class_counts": {
            "0": int((y == 0).sum()),
            "1": int((y == 1).sum()),
        },
        "roc_auc": float(roc_auc_score(y, probabilities)),
        "average_precision": float(average_precision_score(y, probabilities)),
        "brier_score": float(brier_score_loss(y, probabilities)),
        "log_loss": float(log_loss(y, probabilities, labels=[0, 1])),
        "outer_folds": OUTER_FOLDS,
        "inner_folds": INNER_FOLDS,
        "selected_parameters": selected_params,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ledger.to_csv(args.output_dir / "oof_predictions.csv", index=False)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
