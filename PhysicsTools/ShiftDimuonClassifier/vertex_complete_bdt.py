# Code style: AI-code. AI-generated full-RECO category-aware BDT benchmark.
"""Fixed-capacity BDTs on all 119 RECO inputs and five explicit categories.

Training takes labels and supplied class/process weights. Inference takes only
the reconstructed matrix and reconstructed vertex category. There is no GEN,
nuisance, validation or test input, and no internal model/cut selection.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer

from features import FEATURE_NAMES as RAW_FEATURE_NAMES, MUON_QUALITY, signed_log


MUON_DETAILS = (
    "nDTHits", "nCSCHits", "nRPCHits", "nGEMHits", "nME0Hits", "nHitsNearEndcap",
    "nHitsBarrel", "nHitsFarEndcap", "nStationsNearEndcap", "nStationsBarrel",
    "nStationsFarEndcap", "hitSpan", "hitDeltaZ", "duplicateGroupSize",
    "duplicateDSACount", "duplicateTraversingCount", "duplicateCosmicCount",
    "nTimingMeasurements", "timingChi2", "timingDeltaChi2", "timingDirectionSign",
)
EVENT_COUNTS = ("nShiftMuon", "nShiftDT", "nShiftCSC", "nShiftRPC", "nShiftGEM")
COMPLETE_FEATURE_NAMES = RAW_FEATURE_NAMES + tuple(
    "mu_" + leg + "_" + key for leg in ("highPt", "lowPt")
    for key in MUON_QUALITY + MUON_DETAILS + ("relativePtErr",)
) + tuple("event_" + key for key in EVENT_COUNTS)
CATEGORY_NAMES = ("Both-Both", "Near-Both", "Both-Far", "Near-Far", "Other")
BDT_CONFIGS = {
    "shallow_d3": {"max_depth": 3, "max_leaf_nodes": 8, "min_samples_leaf": 30,
                   "l2_regularization": 1.0},
    "medium_d5": {"max_depth": 5, "max_leaf_nodes": 32, "min_samples_leaf": 30,
                  "l2_regularization": 1.0},
    "regularized_d5": {"max_depth": 5, "max_leaf_nodes": 16, "min_samples_leaf": 100,
                       "l2_regularization": 10.0},
}
COMMON_BDT_SETTINGS = {
    "loss": "log_loss", "max_iter": 200, "learning_rate": 0.05,
    "max_bins": 255, "early_stopping": False,
}


def _feature_names(names):
    values = np.asarray(names)
    if values.shape != (119,) or tuple(map(str, values)) != COMPLETE_FEATURE_NAMES:
        raise ValueError("The exact ordered 119 reconstructed feature names are required")
    return COMPLETE_FEATURE_NAMES


def _matrix(values):
    raw = np.asarray(values)
    if raw.ndim != 2 or raw.shape[1] != 119 or raw.dtype.kind not in "iuf":
        raise ValueError("X must be a real reconstructed matrix with exactly 119 columns")
    raw = raw.astype(np.float64, copy=False)
    # Keep the existing feature convention: nonfinite fit values are missing
    # observations, filled using positive-weight training rows only.
    return np.where(np.isfinite(raw), raw, np.nan)


def _category_onehot(topology, n_rows):
    values = np.asarray(topology)
    if values.shape != (n_rows,) or values.dtype.kind not in "US":
        raise ValueError("topology must be an aligned vector of reconstructed categories")
    lookup = {name: index for index, name in enumerate(CATEGORY_NAMES)}
    encoded = np.zeros((n_rows, len(CATEGORY_NAMES)), dtype=np.float64)
    for row, category in enumerate(values):
        name = str(category)
        if name not in lookup:
            raise ValueError("Unknown reconstructed category: " + name)
        encoded[row, lookup[name]] = 1.0
    return encoded


class CompleteVertexBDT:
    """Pickleable fixed classifier; no training labels or nuisance arrays stored."""

    def __init__(self, estimator, imputer, feature_names, config_name, seed, training_rows):
        self.estimator = estimator
        self.imputer = imputer
        self.feature_names_in_ = np.asarray(feature_names)
        self.n_features_in_ = len(feature_names)
        self.categories_ = CATEGORY_NAMES
        self.classes_ = np.asarray([0, 1], dtype=np.int64)
        self.config_name = config_name
        self.seed = int(seed)
        self.requires_gen = False
        self.training_metadata_ = {
            "positive_weight_rows": int(training_rows), "reconstructed_features": 119,
            "explicit_category_columns": 5, "estimator_input_columns": 124,
            "feature_transform": "signedlog_then_training_only_median; empty columns filled zero",
            "category_names": list(CATEGORY_NAMES), "config_name": config_name,
            "settings": {**COMMON_BDT_SETTINGS, **BDT_CONFIGS[config_name], "random_state": int(seed)},
            "internal_early_stopping": False, "iterations_completed": int(estimator.n_iter_),
            "validation_or_test_used": False,
            "weight_policy": "Supplied positive weights unchanged; use the same unit-mean convention across fits.",
            "score_interpretation": "Weighted classification score; not a calibrated physical posterior.",
        }

    def predict_proba(self, X, topology, *, feature_names=None):
        if feature_names is not None:
            _feature_names(feature_names)
        raw = _matrix(X)
        categories = _category_onehot(topology, len(raw))
        if not len(raw):
            return np.empty((0, 2), dtype=np.float64)
        values = self.imputer.transform(signed_log(raw))
        return self.estimator.predict_proba(np.column_stack((values, categories)))


def fit_complete_bdt(X119, topology, y, weights, *, feature_names, seed=71, config="shallow_d3"):
    """Fit a predeclared 200-iteration BDT using training rows only.

    ``feature_names`` must match the full 119-input RECO allowlist exactly.
    ``config`` names one of the three declared ``BDT_CONFIGS``; custom overrides
    are not accepted. Categories always use the fixed five-column vocabulary.
    Root-level study code defines physical-event splits and balanced class/
    process weights. Supplied weights are used unchanged, so regularization
    comparisons require a common weight convention (normally unit mean).
    Zero-weight rows are excluded from both preprocessing and estimator fitting.
    """
    names = _feature_names(feature_names)
    if not isinstance(config, str) or config not in BDT_CONFIGS:
        raise ValueError("config must name one of the three predeclared BDT settings")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or not 0 <= seed < 2**32:
        raise ValueError("seed must be an integer in [0, 2**32)")
    raw = _matrix(X119)
    categories = _category_onehot(topology, len(raw))
    labels = np.asarray(y)
    supplied = np.asarray(weights)
    if (labels.shape != (len(raw),) or labels.dtype.kind not in "iu"
            or not np.isin(labels, [0, 1]).all()):
        raise ValueError("Training labels must be aligned integer 0/1; exclude unknown labels")
    if supplied.shape != (len(raw),) or supplied.dtype.kind not in "iuf":
        raise ValueError("Training weights must be an aligned real vector")
    supplied = supplied.astype(np.float64, copy=False)
    if not np.isfinite(supplied).all() or np.any(supplied < 0):
        raise ValueError("Training weights must be finite and nonnegative")
    active = supplied > 0
    if any(not np.any(active & (labels == label)) for label in (0, 1)):
        raise ValueError("Both labels require positive-weight training rows")
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    values = imputer.fit_transform(signed_log(raw[active]))
    estimator = HistGradientBoostingClassifier(
        **COMMON_BDT_SETTINGS, **BDT_CONFIGS[config], random_state=int(seed)
    )
    estimator.fit(np.column_stack((values, categories[active])), labels[active], sample_weight=supplied[active])
    return CompleteVertexBDT(estimator, imputer, names, config, seed, int(active.sum()))
