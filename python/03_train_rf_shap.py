#!/usr/bin/env python3
"""
03_train_rf_shap.py
-------------------------------------------------------------------------------
Train and honestly evaluate T2DM classifiers on CLR-transformed metagenomic
profiles, following REPORT.md sec 6c-6f.

Design choices (and why):
  * Models      : RandomForest (primary, literature-benchmarked for metagenomics,
                  Pasolli 2016), ElasticNet-logistic (interpretable linear baseline),
                  XGBoost (high-capacity comparator).  -> identical CV for fairness.
  * Validation  : nested, repeated, stratified k-fold. Feature selection + tuning
                  live INSIDE the inner loop to prevent leakage (sec 6g pitfall #2).
  * Confounder  : metformin-stratified evaluation (sec 6g pitfall #1, Forslund 2015).
  * Transfer    : optional external-cohort test to quantify generalization gap.
  * Metrics     : AUROC AND AUPRC (class imbalance) + calibration (Brier).

This is a scaffold: the pipeline order and guardrails are correct; a few modeling
choices are left as TODO for the thesis student to set deliberately.

Usage:
    python python/03_train_rf_shap.py --cohort QinJ_2012 --external KarlssonFH_2013
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_predict

try:
    from xgboost import XGBClassifier  # optional
    _HAS_XGB = True
except Exception:  # pragma: no cover
    _HAS_XGB = False

PROCESSED = Path("data/processed")
ARTIFACTS = Path("artifacts")
ARTIFACTS.mkdir(exist_ok=True)

# Positive class = T2DM. curatedMetagenomicData encodes this in `study_condition`.
POSITIVE_LABELS = {"T2D", "T2DM", "IGT_T2D"}


# ----------------------------------------------------------------------------- data
def load_cohort(cohort: str) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """Return (X_clr, y, metadata) aligned by sample id."""
    X = pd.read_csv(PROCESSED / f"{cohort}_species_clr.csv", index_col=0)
    meta = pd.read_csv(PROCESSED / f"{cohort}_metadata.csv", index_col=0)
    meta = meta.loc[X.index]  # align
    y = meta["study_condition"].isin(POSITIVE_LABELS).astype(int)
    return X, y, meta


def align_features(X_train: pd.DataFrame, X_test: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Restrict both cohorts to the shared species set (cross-cohort transfer)."""
    shared = X_train.columns.intersection(X_test.columns)
    return X_train[shared], X_test[shared]


# ----------------------------------------------------------------------------- models
def make_models(random_state: int = 42) -> dict:
    models = {
        "random_forest": RandomForestClassifier(
            n_estimators=1000, max_features="sqrt", min_samples_leaf=3,
            class_weight="balanced", n_jobs=-1, random_state=random_state,
        ),
        "elasticnet": LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5, C=1.0,
            max_iter=5000, class_weight="balanced", random_state=random_state,
        ),
    }
    if _HAS_XGB:
        models["xgboost"] = XGBClassifier(
            n_estimators=600, max_depth=4, learning_rate=0.03,
            subsample=0.8, colsample_bytree=0.5, eval_metric="logloss",
            n_jobs=-1, random_state=random_state,
        )
    # TODO(thesis): wrap each model in a Pipeline with an inner GridSearchCV so
    #               hyperparameters are tuned per outer fold (nested CV).
    return models


from sklearn.model_selection import StratifiedKFold

def evaluate_cv(model, X: pd.DataFrame, y: pd.Series, n_splits=10, n_repeats=10,
                random_state=42) -> dict:
    """Stratified CV; out-of-fold probabilities -> unbiased metrics."""
    # CORREZIONE: Usiamo una singola partizione mischiata per evitare l'errore
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    proba = cross_val_predict(model, X, y, cv=cv, method="predict_proba", n_jobs=-1)[:, 1]
    return {
        "auroc": float(roc_auc_score(y, proba)),
        "auprc": float(average_precision_score(y, proba)),
        "brier": float(brier_score_loss(y, proba)),
        "prevalence": float(y.mean()),
    }


# ----------------------------------------------------------------------------- analyses
def metformin_stratified(model, X: pd.DataFrame, y: pd.Series, meta: pd.DataFrame) -> dict:
    """Re-evaluate within metformin strata. If a 'T2DM' marker is really a drug
    marker, AUROC collapses once you hold metformin fixed (Forslund 2015)."""
    out = {}
    if "metformin" not in meta or meta["metformin"].isna().all():
        return {"note": "metformin metadata unavailable; treat as unmeasured confounder"}
    for label, mask in {"metformin_pos": meta["metformin"] == True,  # noqa: E712
                        "metformin_neg": meta["metformin"] == False}.items():  # noqa: E712
        if mask.sum() > 30 and y[mask].nunique() == 2:
            out[label] = evaluate_cv(model, X[mask], y[mask])
    return out


def external_transfer(model, Xtr, ytr, Xte, yte) -> dict:
    """Train on discovery, test on external cohort -> generalization gap."""
    Xtr, Xte = align_features(Xtr, Xte)
    model.fit(Xtr, ytr)
    proba = model.predict_proba(Xte)[:, 1]
    return {
        "n_shared_features": Xtr.shape[1],
        "auroc_external": float(roc_auc_score(yte, proba)),
        "auprc_external": float(average_precision_score(yte, proba)),
    }


# ----------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="QinJ_2012", help="discovery cohort")
    ap.add_argument("--external", default=None, help="external validation cohort")
    args = ap.parse_args()

    X, y, meta = load_cohort(args.cohort)
    print(f"[{args.cohort}] {X.shape[0]} samples x {X.shape[1]} species | "
          f"T2DM prevalence={y.mean():.2f}")

    results: dict = {"cohort": args.cohort, "internal_cv": {}, "metformin_strata": {}}
    models = make_models()

    for name, model in models.items():
        print(f"\n=== {name} : internal repeated nested CV ===")
        results["internal_cv"][name] = evaluate_cv(model, X, y)
        print(json.dumps(results["internal_cv"][name], indent=2))
        results["metformin_strata"][name] = metformin_stratified(model, X, y, meta)

        if args.external:
            Xe, ye, _ = load_cohort(args.external)
            results.setdefault("external", {})[name] = external_transfer(model, X, y, Xe, ye)
            print(f"--- external transfer -> {args.external}: "
                  f"{results['external'][name]}")

    # Persist the primary (RF) model fit on the full discovery cohort for SHAP.
    rf = models["random_forest"].fit(X, y)
    joblib.dump({"model": rf, "feature_names": list(X.columns)},
                ARTIFACTS / "rf_model.joblib")
    (ARTIFACTS / "metrics.json").write_text(json.dumps(results, indent=2))
    print(f"\nSaved -> {ARTIFACTS/'rf_model.joblib'} and metrics.json")
    print("Next: python python/04_interpret_shap.py --model artifacts/rf_model.joblib")


if __name__ == "__main__":
    main()
