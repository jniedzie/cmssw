"""Reconstructed-only neural benchmark with a joint DisCo regularizer.

The predictor sees only X. The three columns of Z are reconstructed log1p
mass, signed-log1p vertex z, and log1p vertex radius. A single Euclidean
distance in this three-dimensional space protects their JOINT distribution
separately in each binary training-label class. These are displacement
proxies, not a measurement or guarantee of lifetime independence.
"""

from __future__ import annotations

import math

import numpy as np
import torch
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from torch import nn

from adversarial import (AdversarialClassifier, _classifier, _labels_and_weights,
                         _matrix, _weighted_mean)


MAX_DISTANCE_ROWS = 256
MIN_DISTANCE_ROWS = 4


def _distance_correlation_details(scores, nuisance, weights):
    """Weighted, biased empirical dCor squared and a support/guard reason.

    Double centering uses the same normalized empirical weights in both
    directions. The resulting covariance is measured between scalar-score
    distances and multivariate nuisance distances. Computation is float64 to
    avoid cancellation for nearly constant scores; conversion keeps autograd.
    Missing nuisance observations and zero-weight rows contribute neither
    distances nor centering. No targets are imputed.
    """
    if scores.ndim not in (1, 2) or (scores.ndim == 2 and scores.shape[1] != 1):
        raise ValueError("scores must have one scalar per row")
    if nuisance.ndim != 2 or nuisance.shape[1] < 1:
        raise ValueError("nuisance must be a matrix with at least one column")
    if weights.ndim != 1 or len(weights) != len(scores) or len(nuisance) != len(scores):
        raise ValueError("scores, nuisance and weights must have the same rows")
    if not torch.isfinite(scores).all():
        raise ValueError("scores must be finite")
    if not torch.isfinite(weights).all() or (weights < 0).any():
        raise ValueError("weights must be finite and nonnegative")
    mask = torch.isfinite(nuisance).all(dim=1) & (weights > 0)
    zero = scores.sum() * 0.0
    if int(mask.sum()) < MIN_DISTANCE_ROWS:
        return zero, "fewer than four positive-weight complete nuisance rows"
    score = scores.reshape(-1)[mask].to(torch.float64)
    z = nuisance[mask].to(torch.float64)
    p = weights[mask].to(torch.float64)
    p = p / p.sum()
    # Explicit broadcasting avoids cdist's matrix-multiplication approximation
    # for larger batches, which can give nonzero self-distances.
    a = (score[:, None] - score[None, :]).abs()
    b = torch.linalg.vector_norm(z[:, None, :] - z[None, :, :], dim=2)

    def centered(distance):
        row_mean = distance @ p
        return distance - row_mean[:, None] - row_mean[None, :] + p @ row_mean

    a, b = centered(a), centered(b)
    pair_weights = p[:, None] * p[None, :]
    variance_score = (pair_weights * a.square()).sum()
    variance_nuisance = (pair_weights * b.square()).sum()
    epsilon = torch.finfo(torch.float64).eps
    if float(variance_score.detach()) <= epsilon:
        return zero, "constant or numerically constant scores"
    if float(variance_nuisance.detach()) <= epsilon:
        return zero, "constant or numerically constant joint nuisance"
    covariance = (pair_weights * a * b).sum()
    value = covariance / torch.sqrt(variance_score * variance_nuisance)
    # Rounding can put a nonnegative theoretical covariance just below zero,
    # or perfect correlation just above one. Never clamp an interior value.
    return value.clamp(0.0, 1.0), None


def weighted_distance_correlation_squared(scores, nuisance, weights):
    """Return a differentiable joint dCor-squared scalar for a bounded batch.

    Inputs are torch tensors. This biased finite-sample estimate is a training
    penalty, not an independence test. Undefined/insufficient support returns
    a differentiable zero. Training and diagnostic callers report its reason.
    """
    return _distance_correlation_details(scores, nuisance, weights)[0]


def _support(labels, weights, nuisance, label):
    members = (labels == label) & (weights > 0)
    complete = np.isfinite(nuisance).all(axis=1)
    valid = members & complete
    usable_weights = weights[valid]
    return {
        "rows": int(members.sum()), "valid_rows": int(valid.sum()),
        "missing_rows": int((members & ~complete).sum()),
        "valid_weight_fraction": (float(usable_weights.sum() / weights[members].sum())
                                  if members.any() else None),
        "effective_valid_rows": (float(usable_weights.sum() ** 2 / (usable_weights ** 2).sum())
                                 if len(usable_weights) else 0.0),
    }


def _diagnostic_indices(labels, weights, nuisance, seed):
    """Fixed train/validation-only subsets; quadratic diagnostics stay bounded.

    Uniform row subsampling, followed by the original empirical weights,
    estimates the original weighted distribution without weighting it twice.
    Subset support is explicit; a subsample statistic is never called full-set
    independence. The seed does not consume the classifier-training RNG.
    """
    rng = np.random.default_rng(seed)
    indices = []
    for label in (0, 1):
        valid = np.flatnonzero((labels == label) & (weights > 0) &
                               np.isfinite(nuisance).all(axis=1))
        if len(valid) > MAX_DISTANCE_ROWS:
            valid = np.sort(rng.choice(valid, MAX_DISTANCE_ROWS, replace=False))
        indices.append(torch.from_numpy(valid))
    return indices


def train_disco(X_train, y_train, w_train, Z_train,
                X_val, y_val, w_val, Z_val, *, seed=42, strength=1.0, epochs=150):
    """Return (X-only model, history) for a small deterministic CPU benchmark.

    Train weighted BCE + lambda times the mean available class-conditional
    joint dCor squared. Lambda ramps over the first 20% of epochs; zero lambda
    gives an identical classifier optimization irrespective of nuisance data.
    X imputation/scaling and weighted Z scaling are fitted on training only.
    All inputs are array data; no truth branches or metadata are read. Binary
    labels are the only truth-derived values consumed by this API.
    """
    if not math.isfinite(strength) or strength < 0:
        raise ValueError("strength must be finite and nonnegative")
    if not isinstance(epochs, (int, np.integer)) or epochs < 1:
        raise ValueError("epochs must be a positive integer")
    train_raw = _matrix(X_train, name="X_train")
    val_raw = _matrix(X_val, name="X_val", n_features=train_raw.shape[1])
    if not len(train_raw) or not len(val_raw):
        raise ValueError("training and validation sets must both be nonempty")
    yt, wt = _labels_and_weights(y_train, w_train, len(train_raw), "training")
    yv, wv = _labels_and_weights(y_val, w_val, len(val_raw), "validation")
    if len(np.unique(yt[wt > 0])) != 2:
        raise ValueError("training must have positive-weight examples in both classes")
    zt = _matrix(Z_train, name="Z_train", n_features=3)
    zv = _matrix(Z_val, name="Z_val", n_features=3)
    if len(zt) != len(yt) or len(zv) != len(yv):
        raise ValueError("Z must have one row per classifier example")

    torch.set_num_threads(1)
    torch.manual_seed(int(seed))
    rng = np.random.default_rng(seed)
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    scaler = StandardScaler()
    xt = scaler.fit_transform(imputer.fit_transform(train_raw)).astype(np.float32)
    xv = scaler.transform(imputer.transform(val_raw)).astype(np.float32)
    # A complete-row nuisance scaler preserves the actual joint observations.
    # No validation information, missing-value imputation or gen value enters it.
    valid_train_z = (wt > 0) & np.isfinite(zt).all(axis=1)
    nuisance_scaling = None
    if valid_train_z.any():
        nuisance_scaler = StandardScaler().fit(zt[valid_train_z], sample_weight=wt[valid_train_z])
        nuisance_scaling = {"mean": nuisance_scaler.mean_.tolist(),
                            "scale": nuisance_scaler.scale_.tolist(),
                            "fitted_complete_train_rows": int(valid_train_z.sum())}
        zt = (zt - nuisance_scaler.mean_) / nuisance_scaler.scale_
        zv = (zv - nuisance_scaler.mean_) / nuisance_scaler.scale_

    network = _classifier(xt.shape[1])
    optimizer = torch.optim.Adam(network.parameters(), lr=1e-3)
    bce = nn.BCEWithLogitsLoss(reduction="none")
    tx, vx = torch.from_numpy(xt), torch.from_numpy(xv)
    ty, vy = torch.from_numpy(yt.astype(np.float32)), torch.from_numpy(yv.astype(np.float32))
    tw, vw = torch.from_numpy(wt.astype(np.float32)), torch.from_numpy(wv.astype(np.float32))
    tz, vz = torch.from_numpy(zt), torch.from_numpy(zv)
    diagnostic_indices = {
        "train": _diagnostic_indices(yt, wt, zt, int(seed) + 104729),
        "val": _diagnostic_indices(yv, wv, zv, int(seed) + 130363),
    }
    support = [{"label": label, "train": _support(yt, wt, zt, label),
                "val": _support(yv, wv, zv, label)} for label in (0, 1)]
    records = []
    batch_size = min(MAX_DISTANCE_ROWS, len(yt))
    ramp_epochs = max(1, int(math.ceil(epochs * 0.2)))
    for epoch in range(epochs):
        effective_strength = float(strength * min(1.0, epoch / ramp_epochs))
        network.train()
        order = rng.permutation(len(yt))
        updates_with_penalty = {0: 0, 1: 0}
        for start in range(0, len(yt), batch_size):
            idx = torch.from_numpy(order[start:start + batch_size])
            if float(tw[idx].sum()) <= 0:
                continue
            logits = network(tx[idx]).ravel()
            classification = _weighted_mean(bce(logits, ty[idx]), tw[idx])
            penalties = []
            # At lambda zero nuisance availability cannot affect optimization.
            if effective_strength > 0:
                score = torch.sigmoid(logits)
                for label in (0, 1):
                    members = ty[idx] == label
                    value, reason = _distance_correlation_details(
                        score[members], tz[idx][members], tw[idx][members])
                    if reason is None:
                        penalties.append(value)
                        updates_with_penalty[label] += 1
            penalty = torch.stack(penalties).mean() if penalties else logits.sum() * 0.0
            objective = classification + effective_strength * penalty
            if not torch.isfinite(objective):
                raise RuntimeError("nonfinite DisCo objective")
            optimizer.zero_grad()
            objective.backward()
            if any(parameter.grad is not None and not torch.isfinite(parameter.grad).all()
                   for parameter in network.parameters()):
                raise RuntimeError("nonfinite DisCo classifier gradient")
            optimizer.step()

        network.eval()
        with torch.no_grad():
            train_logits, val_logits = network(tx).ravel(), network(vx).ravel()
            train_score, val_score = torch.sigmoid(train_logits), torch.sigmoid(val_logits)
            train_bce = float(_weighted_mean(bce(train_logits, ty), tw))
            val_bce = float(_weighted_mean(bce(val_logits, vy), vw))
            diagnostics, averages = [], {"train": [], "val": []}
            for label in (0, 1):
                diagnostic = {"label": label, "updates_with_penalty": updates_with_penalty[label]}
                for split, scores, weights, nuisance in (
                    ("train", train_score, tw, tz), ("val", val_score, vw, vz),
                ):
                    idx = diagnostic_indices[split][label]
                    value, reason = _distance_correlation_details(scores[idx], nuisance[idx], weights[idx])
                    used_weights = weights[idx].to(torch.float64)
                    diagnostic[split] = {
                        "sample_rows": int(len(idx)),
                        "effective_sample_rows": (float(used_weights.sum().square() / used_weights.square().sum())
                                                  if len(idx) else 0.0),
                        "dcor_squared": float(value) if reason is None else None,
                        "skip_reason": reason,
                    }
                    if reason is None:
                        averages[split].append(float(value))
                diagnostics.append(diagnostic)
            train_dcor = float(np.mean(averages["train"])) if averages["train"] else None
            val_dcor = float(np.mean(averages["val"])) if averages["val"] else None
            valid_auc = len(np.unique(yv[wv > 0])) == 2
            auc = (float(roc_auc_score(yv, val_score.numpy(), sample_weight=wv))
                   if valid_auc else None)
            records.append({"epoch": epoch + 1, "effective_strength": effective_strength,
                            "train_bce": train_bce, "val_bce": val_bce, "val_auc": auc,
                            "train_dcor_squared": train_dcor, "val_dcor_squared": val_dcor,
                            "decorrelation": diagnostics})

    model = AdversarialClassifier(network, imputer, scaler, train_raw.shape[1])
    history = {
        "method": "weighted_class_conditional_joint_distance_correlation_squared",
        "seed": int(seed), "strength": float(strength), "epochs_requested": int(epochs),
        "architecture": [int(xt.shape[1]), 32, 16, 1], "batch_size": batch_size,
        "nuisance_definition": ["reconstructed log1p mass", "reconstructed signedlog1p vz",
                                "reconstructed log1p vertex radius"],
        "nuisance_scaling": nuisance_scaling, "support": support,
        "diagnostic_max_rows_per_class": MAX_DISTANCE_ROWS,
        "distance_estimator": "biased weighted empirical dCor squared; joint Euclidean nuisance distance",
        "epochs": records,
    }
    return model, history
