#!/usr/bin/env python3
"""
04_interpret_shap.py
-------------------------------------------------------------------------------
Interpret the trained T2DM Random Forest with SHAP (REPORT.md sec 6e).

Produces:
  figures/shap_beeswarm.png     global importance + direction (which taxa push
                                predicted T2DM risk up vs down)
  figures/shap_bar.png          mean(|SHAP|) ranked feature importance
  figures/shap_waterfall_*.png  per-PATIENT explanation (clinically readable)
  figures/shap_dependence_*.png nonlinearity / interaction (e.g., Akkermansia x metformin)
  artifacts/shap_ranked.csv     ranked taxa with mean|SHAP| and mean signed SHAP

CLINICAL CAVEAT (printed at runtime): SHAP explains the MODEL, not biology. A
high-SHAP taxon is a within-model correlate, not a proven cause. Keep this
separate from the causal grading in REPORT.md sec 2-3.

Usage:
    python python/04_interpret_shap.py --model artifacts/rf_model.joblib \
        --cohort QinJ_2012 --n-patients 3
"""
from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import shap
except Exception as e:  # pragma: no cover
    raise SystemExit("Install plotting deps: pip install shap matplotlib") from e

PROCESSED = Path("data/processed")
FIGS = Path("figures"); FIGS.mkdir(exist_ok=True)
ARTIFACTS = Path("artifacts")


def load_X(cohort: str, feature_names: list[str]) -> pd.DataFrame:
    X = pd.read_csv(PROCESSED / f"{cohort}_species_clr.csv", index_col=0)
    # align to the features the model was trained on
    for f in feature_names:
        if f not in X.columns:
            X[f] = 0.0
    return X[feature_names]


def shorten(names: list[str]) -> list[str]:
    """MetaPhlAn clades are long; keep the species tail for legible plots."""
    return [n.split("|")[-1].replace("s__", "").replace("_", " ") for n in names]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="artifacts/rf_model.joblib")
    ap.add_argument("--cohort", default="QinJ_2012")
    ap.add_argument("--n-patients", type=int, default=3)
    args = ap.parse_args()

    bundle = joblib.load(args.model)
    model, feat = bundle["model"], bundle["feature_names"]
    X = load_X(args.cohort, feat)
    X_disp = X.copy(); X_disp.columns = shorten(feat)

    print("CLINICAL CAVEAT: SHAP explains the model, not biology. "
          "High-SHAP != causal (see REPORT.md sec 2-3).")

    # TreeSHAP: exact and fast for tree ensembles.
    explainer = shap.TreeExplainer(model)
    sv = explainer(X)                       # shap.Explanation
    # binary RF: take the positive-class (T2DM) channel if 3-D
    values = sv.values[..., 1] if sv.values.ndim == 3 else sv.values

    # --- global beeswarm (direction) ---
    plt.figure()
    shap.summary_plot(values, X_disp, show=False, max_display=20)
    plt.tight_layout(); plt.savefig(FIGS / "shap_beeswarm.png", dpi=160); plt.close()

    # --- global bar (mean|SHAP|) ---
    plt.figure()
    shap.summary_plot(values, X_disp, plot_type="bar", show=False, max_display=20)
    plt.tight_layout(); plt.savefig(FIGS / "shap_bar.png", dpi=160); plt.close()

    # --- ranked table: mean|SHAP| and mean signed SHAP (direction) ---
    ranked = (pd.DataFrame({
        "taxon": feat,
        "mean_abs_shap": np.abs(values).mean(axis=0),
        "mean_signed_shap": values.mean(axis=0),  # >0 pushes toward T2DM
    }).sort_values("mean_abs_shap", ascending=False))
    ranked.to_csv(ARTIFACTS / "shap_ranked.csv", index=False)
    print("\nTop 10 model-important taxa (sign>0 => raises predicted T2DM risk):")
    print(ranked.head(10).to_string(index=False))

    # --- per-patient waterfalls (clinically readable explanations) ---
    for i in range(min(args.n_patients, X.shape[0])):
        plt.figure()
        shap.plots.waterfall(sv[i, ..., 1] if sv.values.ndim == 3 else sv[i],
                             max_display=12, show=False)
        plt.tight_layout()
        plt.savefig(FIGS / f"shap_waterfall_patient{i}.png", dpi=160); plt.close()

    # --- dependence plot for the top taxon (nonlinearity / interaction) ---
    top = int(np.argmax(np.abs(values).mean(axis=0)))
    plt.figure()
    shap.dependence_plot(top, values, X_disp, show=False)
    plt.tight_layout()
    plt.savefig(FIGS / f"shap_dependence_{X_disp.columns[top]}.png", dpi=160); plt.close()
    # TODO(thesis): pass interaction_index="metformin" once metformin is a feature
    #               to visualize Akkermansia x metformin confounding (Forslund 2015).

    print(f"\nFigures -> {FIGS}/  |  Ranked taxa -> {ARTIFACTS/'shap_ranked.csv'}")


if __name__ == "__main__":
    main()
