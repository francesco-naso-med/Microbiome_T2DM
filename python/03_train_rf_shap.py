#!/usr/bin/env python3
"""
03_train_rf_shap.py
-------------------------------------------------------------------------------
Train and evaluate T2DM classifiers on metagenomic profiles (REPORT.md sec 6c-6f).

Validation design
  * Input       : RAW MetaPhlAn relative abundances (<cohort>_species.csv).
                  All compositional preprocessing (prevalence filter, zero
                  replacement, CLR) lives in a Pipeline and is re-fitted inside
                  every training fold -> no leakage from test folds.
  * Internal    : repeated stratified k-fold (default 10 x 10). Each repeat gives
                  a full set of out-of-fold probabilities; metrics are computed
                  per repeat and reported as mean +/- SD across repeats.
  * Tuning      : OFF by default (fixed, literature-standard hyperparameters).
                  With --tune, an inner GridSearchCV runs inside each outer
                  training fold -> genuinely nested CV (slow).
  * External    : the pipeline is fitted on the discovery cohort and applied
                  as-is to the external cohort. Species absent from the external
                  table are treated as not detected (0); the CLR reference is the
                  same retained species set in both cohorts.
  * Labels      : binary T2D vs control from `study_condition`. IGT and samples
                  with missing condition are EXCLUDED (counts printed at runtime).
  * Medication  : metformin / statins derived from `treatment` as THREE-STATE
                  flags (1 / 0 / NaN = not annotated). Unknown is never "untreated".
  * Confounding : no retraining inside tiny strata. The out-of-fold probabilities
                  of the main model (averaged over repeats) are scored separately
                  for T2D subgroups (metformin+, confirmed metformin-, unknown,
                  drug-free) against controls, with stratified bootstrap 95% CI.
                  A model that detects the drug, not the disease, scores much
                  higher on metformin+ cases than on drug-free cases.
  * Metrics     : AUROC (+ bootstrap CI), AUPRC (chance = prevalence), Brier.

Usage:
    python python/03_train_rf_shap.py --cohort QinJ_2012 --external KarlssonFH_2013
    python python/03_train_rf_shap.py --cohort QinJ_2012 --repeats 3 --tune
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import sklearn

# sklearn >= 1.8 deprecates `penalty`: elastic net is implied by 0 < l1_ratio < 1.
_SK = tuple(int(v) for v in sklearn.__version__.split(".")[:2])
_ENET_KW = {} if _SK >= (1, 8) else {"penalty": "elasticnet"}

sys.path.insert(0, str(Path(__file__).resolve().parent))
from clr import CLRTransformer  # noqa: E402

try:
    from xgboost import XGBClassifier
    _HAS_XGB = True
except Exception:  # pragma: no cover
    _HAS_XGB = False

PROCESSED = Path("data/processed")
ARTIFACTS = Path("artifacts")
ARTIFACTS.mkdir(exist_ok=True)

CASE, CONTROL = "T2D", "control"     # anything else (IGT, NaN) is excluded
N_BOOT = 2000
MIN_GROUP = 10                       # smallest case/control group scored


# ----------------------------------------------------------------------------- data
def _drug_flag(treatment: pd.Series, drug: str) -> pd.Series:
    """1 / 0 / NaN. NaN treatment stays NaN: unknown is NOT untreated."""
    t = treatment.astype("string").str.lower()
    return t.str.contains(drug, regex=False).astype("Float64").astype(float)


def load_cohort(cohort: str) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """Return (X_raw_relative_abundance, y, metadata) for T2D vs control only."""
    X = pd.read_csv(PROCESSED / f"{cohort}_species.csv", index_col=0)
    meta = pd.read_csv(PROCESSED / f"{cohort}_metadata.csv", index_col=0).loc[X.index]

    cond = meta["study_condition"]
    keep = cond.isin([CASE, CONTROL])
    dropped = cond[~keep].fillna("<missing>").value_counts().to_dict()
    if dropped:
        print(f"[{cohort}] excluded from T2D-vs-control: {dropped}")
    X, meta = X[keep], meta[keep].copy()

    tr = meta["treatment"] if "treatment" in meta else pd.Series(np.nan, index=meta.index)
    meta["metformin"] = _drug_flag(tr, "metformin")
    meta["statins"] = _drug_flag(tr, "statin")
    y = (meta["study_condition"] == CASE).astype(int)
    return X, y, meta


# ----------------------------------------------------------------------------- models
def make_pipelines(y: pd.Series, tune: bool, random_state: int = 42) -> dict:
    """One Pipeline per model: CLR (+ scaler for the linear model) -> classifier.
    Models use n_jobs=1: parallelism is at the CV level, avoiding oversubscription."""
    pos_weight = float((y == 0).sum() / max((y == 1).sum(), 1))

    specs = {
        "random_forest": (
            [("clr", CLRTransformer()),
             ("clf", RandomForestClassifier(
                 n_estimators=1000, max_features="sqrt", min_samples_leaf=3,
                 class_weight="balanced", n_jobs=1, random_state=random_state))],
            {"clf__max_features": ["sqrt", 0.1], "clf__min_samples_leaf": [1, 3, 5]},
        ),
        "elasticnet": (
            [("clr", CLRTransformer()),
             ("scale", StandardScaler()),
             ("clf", LogisticRegression(
                 **_ENET_KW, solver="saga", l1_ratio=0.5, C=1.0,
                 max_iter=5000, class_weight="balanced", random_state=random_state))],
            {"clf__C": [0.01, 0.1, 1.0], "clf__l1_ratio": [0.1, 0.5, 0.9]},
        ),
    }
    if _HAS_XGB:
        specs["xgboost"] = (
            [("clr", CLRTransformer()),
             ("clf", XGBClassifier(
                 n_estimators=600, max_depth=4, learning_rate=0.03, subsample=0.8,
                 colsample_bytree=0.5, scale_pos_weight=pos_weight,
                 eval_metric="logloss", n_jobs=1, random_state=random_state))],
            {"clf__max_depth": [2, 4], "clf__learning_rate": [0.03, 0.1]},
        )

    out = {}
    for name, (steps, grid) in specs.items():
        pipe = Pipeline(steps)
        if tune:
            pipe = GridSearchCV(
                pipe, grid, scoring="roc_auc", n_jobs=1, refit=True,
                cv=StratifiedKFold(5, shuffle=True, random_state=random_state))
        out[name] = pipe
    return out


# ----------------------------------------------------------------------------- evaluation
def _metrics(y, proba) -> dict:
    return {"auroc": roc_auc_score(y, proba),
            "auprc": average_precision_score(y, proba),
            "brier": brier_score_loss(y, proba)}


def evaluate_cv(model, X: pd.DataFrame, y: pd.Series, n_splits: int = 10,
                n_repeats: int = 10, random_state: int = 42) -> tuple[dict, pd.Series]:
    """Repeated stratified k-fold. cross_val_predict requires a partition, so each
    repeat is one shuffled StratifiedKFold; metrics are aggregated across repeats.
    Also returns the out-of-fold probability averaged over repeats (one per sample)."""
    per_repeat, oof = [], []
    for r in range(n_repeats):
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state + r)
        proba = cross_val_predict(clone(model), X, y, cv=cv,
                                  method="predict_proba", n_jobs=-1)[:, 1]
        per_repeat.append(_metrics(y, proba))
        oof.append(proba)

    df = pd.DataFrame(per_repeat)
    summary = {f"{m}_mean": float(df[m].mean()) for m in df}
    summary.update({f"{m}_sd": float(df[m].std(ddof=1)) if n_repeats > 1 else 0.0 for m in df})
    summary.update({"prevalence": float(y.mean()), "n": int(len(y)),
                    "n_splits": n_splits, "n_repeats": n_repeats})
    return summary, pd.Series(np.mean(oof, axis=0), index=X.index, name="oof_proba")


def bootstrap_auroc(y, p, n_boot: int = N_BOOT, seed: int = 0) -> dict:
    """AUROC with stratified (within-class) percentile bootstrap 95% CI."""
    y, p = np.asarray(y), np.asarray(p)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    rng = np.random.default_rng(seed)
    boots = [roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(neg))],
                           np.r_[p[rng.choice(pos, len(pos))], p[rng.choice(neg, len(neg))]])
             for _ in range(n_boot)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"auroc": float(roc_auc_score(y, p)), "ci95": [float(lo), float(hi)],
            "n_cases": int(len(pos)), "n_controls": int(len(neg))}


def subgroup_analysis(y: pd.Series, proba: pd.Series, meta: pd.DataFrame) -> dict:
    """Score T2D subgroups against controls using the SAME predictions.
    No retraining: tiny strata would give unstable models, not evidence."""
    case, ctrl = y == 1, y == 0
    met, sta = meta["metformin"], meta["statins"]
    groups = {
        "T2D_all_vs_controls":           (case,                          ctrl),
        "T2D_metformin_vs_controls":     (case & (met == 1),             ctrl),
        "T2D_no_metformin_vs_controls":  (case & (met == 0),             ctrl),
        "T2D_metformin_unknown_vs_controls": (case & met.isna(),         ctrl),
        # drug-free on BOTH sides: the cleanest disease-only contrast available
        "T2D_drugfree_vs_controls_drugfree": (case & (met == 0) & (sta == 0),
                                              ctrl & (sta == 0)),
    }
    out = {}
    for name, (cm, km) in groups.items():
        nc, nk = int(cm.sum()), int(km.sum())
        if nc < MIN_GROUP or nk < MIN_GROUP:
            out[name] = {"note": f"not scored (cases={nc}, controls={nk}, min={MIN_GROUP})"}
            continue
        m = cm | km
        out[name] = bootstrap_auroc(y[m], proba[m])
    return out


def _print_subgroups(sub: dict) -> None:
    for name, r in sub.items():
        if "auroc" in r:
            print(f"    {name:<36} AUROC {r['auroc']:.3f} "
                  f"[{r['ci95'][0]:.3f}-{r['ci95'][1]:.3f}]  "
                  f"({r['n_cases']} vs {r['n_controls']})")
        else:
            print(f"    {name:<36} {r['note']}")


def external_transfer(model, Xtr, ytr, Xte, yte, meta_te) -> dict:
    """Fit the whole pipeline on discovery, apply frozen to the external cohort."""
    fitted = clone(model).fit(Xtr, ytr)
    proba = pd.Series(fitted.predict_proba(Xte)[:, 1], index=Xte.index)
    clr = (fitted.best_estimator_ if hasattr(fitted, "best_estimator_") else fitted)["clr"]
    retained = list(clr.get_feature_names_out())
    m = _metrics(yte, proba)
    return {"n_features_model": len(retained),
            "n_features_absent_in_external": int(sum(f not in Xte.columns for f in retained)),
            "auroc_external": bootstrap_auroc(yte, proba),
            "auprc_external": float(m["auprc"]), "brier_external": float(m["brier"]),
            "prevalence_external": float(yte.mean()),
            "subgroups": subgroup_analysis(yte, proba, meta_te)}


# ----------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="QinJ_2012")
    ap.add_argument("--external", default=None)
    ap.add_argument("--folds", type=int, default=10)
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--tune", action="store_true", help="nested CV with inner grid search")
    args = ap.parse_args()

    X, y, meta = load_cohort(args.cohort)
    print(f"[{args.cohort}] {X.shape[0]} samples x {X.shape[1]} species (raw) | "
          f"T2DM prevalence={y.mean():.2f}")
    ext = load_cohort(args.external) if args.external else None
    cv_kw = dict(n_splits=args.folds, n_repeats=args.repeats)

    results: dict = {"cohort": args.cohort, "tuned": args.tune,
                     "labels": {"case": CASE, "control": CONTROL,
                                "n_cases": int(y.sum()), "n_controls": int((y == 0).sum())},
                     "internal_cv": {}, "internal_subgroups": {}}
    pipes = make_pipelines(y, tune=args.tune)

    for name, pipe in pipes.items():
        print(f"\n=== {name}: {args.repeats}x{args.folds} stratified CV"
              f"{' (nested)' if args.tune else ''} ===")
        res, oof = evaluate_cv(pipe, X, y, **cv_kw)
        results["internal_cv"][name] = res
        oof.to_frame().join(y.rename("y")).to_csv(ARTIFACTS / f"oof_{name}.csv")
        print(f"AUROC {res['auroc_mean']:.3f} ± {res['auroc_sd']:.3f} | "
              f"AUPRC {res['auprc_mean']:.3f} ± {res['auprc_sd']:.3f} "
              f"(chance {res['prevalence']:.2f}) | Brier {res['brier_mean']:.3f}")
        sub = subgroup_analysis(y, oof, meta)
        results["internal_subgroups"][name] = sub
        print("  internal subgroups (OOF, bootstrap 95% CI):"); _print_subgroups(sub)

        if ext is not None:
            Xe, ye, me = ext
            r = external_transfer(pipe, X, y, Xe, ye, me)
            results.setdefault("external", {})[name] = r
            a = r["auroc_external"]
            print(f"--- external -> {args.external}: AUROC {a['auroc']:.3f} "
                  f"[{a['ci95'][0]:.3f}-{a['ci95'][1]:.3f}] ({a['n_cases']} vs {a['n_controls']}), "
                  f"{r['n_features_absent_in_external']}/{r['n_features_model']} species absent")
            _print_subgroups(r["subgroups"])

    # Final RF pipeline fitted on the full discovery cohort, for SHAP (04_*).
    rf = clone(pipes["random_forest"]).fit(X, y)
    rf = rf.best_estimator_ if hasattr(rf, "best_estimator_") else rf
    joblib.dump({"pipeline": rf,
                 "feature_names": list(rf["clr"].get_feature_names_out())},
                ARTIFACTS / "rf_model.joblib")
    (ARTIFACTS / "metrics.json").write_text(json.dumps(results, indent=2))
    print(f"\nSaved -> {ARTIFACTS/'rf_model.joblib'} and {ARTIFACTS/'metrics.json'}")
    print("Next: python python/04_interpret_shap.py --cohort", args.cohort)


if __name__ == "__main__":
    main()
