# Code style: AI-code. AI-generated fixed-RECO-schema neural trainer adapter.
"""Use the stable category trainer with the full 119-input/five-category schema.

The category vocabulary is a declared reconstructed schema, not something
learned from validation. A rare category absent from a training fraction still
has its fixed input/head column. Its predictions are supported by the shared
model; absence from training remains explicit in the returned support counts.

This adapter temporarily replaces the trainer's category-encoding binding and
restores it in ``finally``. Run fits sequentially in a process. Concurrent
threads using this adapter or the original category trainer are unsupported;
separate worker processes may fit independent models safely. No source file is
modified. The original network, loss, preprocessing and checkpoint code run
unchanged, with the exact same train/validation-only argument contract.
"""
from __future__ import annotations

import numpy as np

import neural_vertex_categories as legacy
from vertex_complete_bdt import CATEGORY_NAMES, COMPLETE_FEATURE_NAMES


def fixed_category_arrays(train, val, n_train, n_val):
    """Encode the predeclared RECO schema; never derive vocabulary from a split."""
    train, val = np.asarray(train), np.asarray(val)
    if train.shape != (n_train,) or val.shape != (n_val,):
        raise ValueError("Fixed reconstructed categories must align with train/validation rows")
    for values in (train, val):
        if any(not isinstance(value, (str, np.str_)) for value in values):
            raise ValueError("Fixed reconstructed categories must be strings")
        unknown = sorted(set(values) - set(CATEGORY_NAMES))
        if unknown:
            raise ValueError("Unknown reconstructed category: " + ", ".join(unknown))
    vocabulary = np.asarray(CATEGORY_NAMES)
    return ((train[:, None] == vocabulary).astype(float),
            (val[:, None] == vocabulary).astype(float), list(CATEGORY_NAMES))


def train_category_vertex(
    X_train, y_train, w_train, Z_train, X_val, y_val, w_val, Z_val,
    *, feature_names, train_category, val_category, **options,
):
    """Train the original shared/heads network using fixed 119+5 RECO inputs.

    All optimization arguments retain the original ``train_category_vertex``
    meanings. The study must explicitly declare penalties and settings before
    fitting. ``Z`` and process metadata, when used for training regularization,
    never enter the model's prediction interface. The two permitted
    architectures are ``shared`` and ``heads``.
    """
    names = np.asarray(feature_names)
    if names.shape != (119,) or tuple(map(str, names)) != COMPLETE_FEATURE_NAMES:
        raise ValueError("The exact ordered 119 reconstructed feature names are required")
    if options.get("architecture", "shared") not in ("shared", "heads"):
        raise ValueError("Complete neural architecture must be shared or heads")
    train_x, val_x = np.asarray(X_train), np.asarray(X_val)
    if train_x.ndim != 2 or val_x.ndim != 2 or train_x.shape[1] != 119 or val_x.shape[1] != 119:
        raise ValueError("Complete neural train/validation X must have exactly 119 RECO columns")
    fixed_category_arrays(train_category, val_category, len(train_x), len(val_x))
    original_encoder = legacy._category_arrays
    try:
        legacy._category_arrays = fixed_category_arrays
        model, history = legacy.train_category_vertex(
            X_train, y_train, w_train, Z_train, X_val, y_val, w_val, Z_val,
            feature_names=COMPLETE_FEATURE_NAMES,
            train_category=train_category, val_category=val_category, **options,
        )
    finally:
        legacy._category_arrays = original_encoder
    if tuple(model.categories_) != CATEGORY_NAMES or model.n_features_in_ != 119:
        raise RuntimeError("Complete neural model lost its fixed reconstructed schema")
    model.requires_gen = False
    model.fixed_category_schema_ = CATEGORY_NAMES
    history["fixed_reco_schema"] = {
        "adapter": "vertex_complete_neural.train_category_vertex",
        "features": list(COMPLETE_FEATURE_NAMES), "category_names": list(CATEGORY_NAMES),
        "vocabulary_origin": "Predeclared RECO schema; not fitted to validation or test",
        "features_retained": 119, "explicit_category_columns": 5,
        "missing_training_categories": [str(category) for category in CATEGORY_NAMES
                                        if not np.any(np.asarray(train_category) == category)],
        "original_encoder_restored": legacy._category_arrays is original_encoder,
        "concurrency": "Sequential fits per process; concurrent trainer threads unsupported",
    }
    return model, history


# A descriptive alias for callers that also compare fit_complete_bdt.
train_complete_neural = train_category_vertex
