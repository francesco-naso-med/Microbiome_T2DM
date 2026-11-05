"""
clr.py
-------------------------------------------------------------------------------
Leakage-free compositional preprocessing as a scikit-learn transformer.

Everything that LEARNS from data (which species survive the prevalence filter,
how large the zero-replacement value is for each species) is estimated in
`fit` on the TRAINING fold only and then frozen. `transform` applies the frozen
parameters to any sample, including samples from an external cohort.

Steps (REPORT.md sec 6b):
  1. closure          x_i <- x_i / sum(x)         (MetaPhlAn % -> proportions)
  2. prevalence filter (fit): keep species with x > DETECTION in >= PREVALENCE
                       of training samples
  3. re-closure on the retained species
  4. multiplicative zero replacement (Martin-Fernandez et al. 2003):
       zeros -> delta_j,  non-zeros scaled by (1 - sum of deltas in that sample)
     delta_j = DELTA_FRAC * (smallest positive value of species j in training)
  5. CLR:  clr(x)_i = ln x_i - mean_k ln x_k

Why steps 3-5 happen INSIDE transform, after the filter: the CLR reference
(geometric mean) depends on the species set. Computing it on the same retained
set for train, test and external cohorts is what makes the values comparable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


class CLRTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, prevalence: float = 0.10, detection: float = 1e-4,
                 delta_frac: float = 0.65):
        self.prevalence = prevalence
        self.detection = detection
        self.delta_frac = delta_frac

    # ------------------------------------------------------------------ helpers
    @staticmethod
    def _close(X: np.ndarray) -> np.ndarray:
        s = X.sum(axis=1, keepdims=True)
        s[s == 0] = 1.0
        return X / s

    def _align(self, X) -> np.ndarray:
        """Reorder columns to the fitted feature set; absent species = 0."""
        if isinstance(X, pd.DataFrame):
            missing = [c for c in self.features_in_ if c not in X.columns]
            self.n_missing_last_ = len(missing)
            return X.reindex(columns=self.features_in_, fill_value=0.0).to_numpy(float)
        self.n_missing_last_ = 0
        return np.asarray(X, dtype=float)

    # ------------------------------------------------------------------ API
    def fit(self, X, y=None):
        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X)
        self.features_in_ = np.asarray(X.columns)
        P = self._close(X.to_numpy(float))

        keep = (P > self.detection).mean(axis=0) >= self.prevalence
        if keep.sum() < 2:
            raise ValueError("Fewer than 2 species pass the prevalence filter.")
        self.keep_mask_ = keep
        self.feature_names_out_ = self.features_in_[keep]

        Pk = self._close(P[:, keep])
        pos = np.where(Pk > 0, Pk, np.nan)
        min_pos = np.nanmin(pos, axis=0)             # per species, training only
        self.delta_ = self.delta_frac * min_pos
        return self

    def transform(self, X):
        P = self._close(self._align(X))
        Pk = self._close(P[:, self.keep_mask_])

        # samples with none of the retained species: uniform composition
        empty = Pk.sum(axis=1) == 0
        Pk[empty] = 1.0 / Pk.shape[1]

        zero = Pk == 0
        d = np.where(zero, self.delta_[None, :], 0.0)
        Z = np.where(zero, d, Pk * (1.0 - d.sum(axis=1, keepdims=True)))

        L = np.log(Z)
        return L - L.mean(axis=1, keepdims=True)

    def get_feature_names_out(self, input_features=None):
        return self.feature_names_out_
