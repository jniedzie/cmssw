"""A bounded, single-working-point neighborhood-uniformity AdaBoost study.

This is inspired by Stevens and Williams, https://arxiv.org/abs/1305.7248,
especially their local-efficiency reweighting. It is not the paper's full
uBoost ensemble over many efficiency points, and makes no flatness guarantee.
Only genuine-vertex positives receive a uniformity update. Reconstructed
mass/z/radius ``Z`` defines training neighborhoods; prediction uses ``X`` only.
These observables are displacement proxies, not a measurement of lifetime.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import expit, logsumexp
from sklearn.impute import SimpleImputer
from sklearn.neighbors import NearestNeighbors
from sklearn.tree import DecisionTreeClassifier


def _matrix(values, name, n_features=None):
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 2 or result.shape[1] == 0:
        raise ValueError(f"{name} must be a two-dimensional nonempty-column array")
    if n_features is not None and result.shape[1] != n_features:
        raise ValueError(f"{name} has the wrong number of columns")
    return np.where(np.isfinite(result), result, np.nan)


def _normalize_class_logs(log_weights, labels):
    """Keep the two total training-label weights equal without overflow."""
    result = log_weights.copy()
    for label in (0, 1):
        members = labels == label
        result[members] -= logsumexp(result[members]) + math.log(2.0)
    return result


def _weighted_cut(scores, weights, efficiency):
    """Observed quantile; ties are retained and their achieved rate is reported."""
    order = np.argsort(scores, kind="stable")
    cumulative = np.cumsum(weights[order])
    index = np.searchsorted(cumulative, (1.0 - efficiency) * cumulative[-1])
    return float(scores[order[min(int(index), len(order) - 1)]])


def _local_efficiencies(passed, neighbors, base_weights):
    """Fixed original weights define the target population, not boost weights."""
    weights = base_weights[neighbors]
    return (passed[neighbors] * weights).sum(axis=1) / weights.sum(axis=1)


def _uniformity_log_update(local_efficiency, global_efficiency, strength, rate):
    # Bounded per-round correction is a stabilization choice of this variant.
    # It preserves ordering: underaccepted regions always get the larger update.
    return np.clip((global_efficiency - local_efficiency) * strength * rate, -1.0, 1.0)


class UniformityBDT:
    """Pickle-compatible tree ensemble whose inference needs reconstructed X only.

    ``predict_proba`` supplies a monotonic sigmoid of the raw tree margin
    divided by the sum of absolute tree coefficients. The normalized margin
    lies in [-1,1], avoiding saturated probabilities even in long ensembles.
    Class balancing and uniformity weighting mean it is not a calibrated
    physical posterior probability. Load pickle files only from trusted sources.
    Nuisance values, neighborhood indexes and labels are not kept in this model.
    """

    def __init__(self, trees, coefficients, imputer, n_features):
        self.estimators_ = list(trees)
        self.estimator_weights_ = np.asarray(coefficients, dtype=np.float64)
        self.imputer = imputer
        self.n_features_in_ = int(n_features)
        self.classes_ = np.array([0, 1], dtype=np.int64)

    def decision_function(self, X):
        raw = _matrix(X, "X", self.n_features_in_)
        if not len(raw):
            return np.empty(0, dtype=np.float64)
        features = self.imputer.transform(raw)
        result = np.zeros(len(features), dtype=np.float64)
        for tree, coefficient in zip(self.estimators_, self.estimator_weights_):
            result += coefficient * (2.0 * tree.predict(features) - 1.0)
        return result

    def predict_proba(self, X):
        normalizer = float(np.sum(np.abs(self.estimator_weights_)))
        margin = self.decision_function(X)
        if normalizer > 0:
            margin /= normalizer
        score = expit(2.0 * margin)
        return np.column_stack((1.0 - score, score))


def train_uniform_bdt(X_train, y_train, w_train, Z_train, *, seed=42,
                      strength=1.0, n_estimators=160, target_efficiency=0.9):
    """Return ``(model, history)`` using training rows only.

    ``Z_train`` has three columns: reconstructed log1p mass, signed-log1p
    vertex z, log1p transverse radius, matching ``nuisance_transform``.
    Input labels must be 0/1; exclude unknown labels before calling. All
    genuine vertices, including QCD examples, must be labelled positive.

    Each depth-two tree minimizes its weighted classification criterion.
    Its AdaBoost coefficient is 0.25 * 0.5 * log((1-error)/error); weights
    get the standard exp(-coefficient * signed_label * signed_prediction)
    update. A positive-only neighborhood correction then multiplies weights
    by exp(strength * 0.25 * [achieved_global_efficiency-local_efficiency]).
    Its log is clipped to [-1,1]. Both class totals are reset to 1/2 after
    each round. Thus strength=0 is our equal-class AdaBoost reference, not
    sklearn's unmodified implementation. Original within-class input weights
    define efficiency estimates throughout; changing boost weights do not.

    Positive training nuisance values are standardized with weighted means
    and variances. Neighborhoods contain up to 100 positives (including
    self), with k=max(20,min(100,n_valid//5)), capped at n_valid. Missing
    nuisance rows participate in classification but receive no uniformity
    correction. Fewer than 20 valid positives disables that correction.
    Only one requested operating point is targeted, with no adversary or
    lifetime target and no guarantee outside training-supported regions.
    """
    if not math.isfinite(strength) or strength < 0:
        raise ValueError("strength must be finite and nonnegative")
    if not isinstance(n_estimators, (int, np.integer)) or n_estimators < 1:
        raise ValueError("n_estimators must be a positive integer")
    if not math.isfinite(target_efficiency) or not 0 < target_efficiency < 1:
        raise ValueError("target_efficiency must lie strictly between zero and one")
    raw = _matrix(X_train, "X_train")
    if not len(raw):
        raise ValueError("training set must be nonempty")
    y = np.asarray(y_train)
    w = np.asarray(w_train, dtype=np.float64)
    if y.shape != (len(raw),) or w.shape != (len(raw),):
        raise ValueError("labels and weights must have one entry per row")
    if not np.isin(y, (0, 1)).all():
        raise ValueError("training labels must be binary; exclude unknown labels")
    if not np.isfinite(w).all() or (w < 0).any():
        raise ValueError("training weights must be finite and nonnegative")
    if any(not np.any((y == label) & (w > 0)) for label in (0, 1)):
        raise ValueError("both labels need positive-weight examples")
    z = _matrix(Z_train, "Z_train", 3)
    if len(z) != len(raw):
        raise ValueError("Z_train must have one row per example")

    # Training-only feature fit. A zero-weight row cannot change preprocessing.
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    imputer.fit(raw[w > 0])
    X = imputer.transform(raw)
    y = y.astype(np.int64)
    log_weights = np.full(len(w), -np.inf)
    log_weights[w > 0] = np.log(w[w > 0])
    log_weights = _normalize_class_logs(log_weights, y)
    base_weights = np.exp(log_weights)
    positive = (y == 1) & (w > 0)
    valid = positive & np.isfinite(z).all(axis=1)
    valid_rows = np.flatnonzero(valid)
    n_valid = len(valid_rows)
    metadata = dict(active=False, valid_positive_rows=n_valid,
                    missing_positive_rows=int(np.count_nonzero(positive & ~valid)),
                    n_neighbors=None, center=None, scale=None, constant_columns=None,
                    skip_reason=None, includes_self=True)
    neighbors = None
    if n_valid < 20:
        metadata["skip_reason"] = "fewer than 20 positive-weight rows with all nuisance values"
    else:
        wz = base_weights[valid]
        center = np.average(z[valid], axis=0, weights=wz)
        variance = np.average((z[valid] - center) ** 2, axis=0, weights=wz)
        constant = variance <= np.finfo(np.float64).eps
        scale = np.where(constant, 1.0, np.sqrt(variance))
        standardized = (z[valid] - center) / scale
        k = min(n_valid, max(20, min(100, n_valid // 5)))
        metadata.update(center=center.tolist(), scale=scale.tolist(),
                        constant_columns=np.flatnonzero(constant).tolist(), n_neighbors=k)
        if constant.all():
            metadata["skip_reason"] = "all reconstructed nuisance columns are constant"
        else:
            # Pass training rows explicitly so sklearn includes self.
            neighbors = NearestNeighbors(n_neighbors=k, n_jobs=1).fit(standardized).kneighbors(
                standardized, return_distance=False)
            metadata["active"] = bool(strength > 0)
            if strength == 0:
                metadata["skip_reason"] = "strength is zero; ordinary equal-class AdaBoost reference"

    rate = 0.25
    trees, coefficients, rounds = [], [], []
    score = np.zeros(len(y), dtype=np.float64)
    rng = np.random.default_rng(seed)
    signed_label = 2.0 * y - 1.0
    stopping_reason = "requested number of estimators reached"
    min_leaf = max(2, min(12, int(np.count_nonzero(w > 0)) // 30))
    for iteration in range(n_estimators):
        current_weights = np.exp(log_weights)
        tree = DecisionTreeClassifier(max_depth=2, min_samples_leaf=min_leaf,
                                      random_state=int(rng.integers(0, 2**31 - 1)))
        tree.fit(X, y, sample_weight=current_weights)
        predicted = tree.predict(X)
        error = float(np.sum(current_weights * (predicted != y)))
        if error >= 0.5 - 1e-12:
            stopping_reason = "next tree is not better than chance under current weights"
            break
        bounded_error = np.clip(error, 1e-12, 1.0 - 1e-12)
        coefficient = float(rate * 0.5 * np.log((1.0 - bounded_error) / bounded_error))
        signed_prediction = 2.0 * predicted - 1.0
        trees.append(tree)
        coefficients.append(coefficient)
        score += coefficient * signed_prediction
        log_weights -= coefficient * signed_label * signed_prediction
        cut = _weighted_cut(score[positive], base_weights[positive], target_efficiency)
        passed = score >= cut
        achieved = float(np.average(passed[positive], weights=base_weights[positive]))
        record = dict(iteration=iteration + 1, weighted_error=error,
                      coefficient=coefficient, training_cut=cut,
                      achieved_positive_efficiency=achieved,
                      neighborhood_efficiency_rms=None,
                      neighborhood_efficiency_min=None,
                      neighborhood_efficiency_max=None,
                      uniformity_log_update_min=0.0, uniformity_log_update_max=0.0)
        if neighbors is not None:
            local = _local_efficiencies(passed[valid], neighbors, base_weights[valid])
            # The global target includes missing-nuisance positives. Such rows
            # themselves get no correction; coverage is explicit in history.
            correction = _uniformity_log_update(local, achieved, strength, rate)
            log_weights[valid_rows] += correction
            record.update(neighborhood_efficiency_rms=float(np.sqrt(np.average(
                              (local - achieved) ** 2, weights=base_weights[valid]))),
                          neighborhood_efficiency_min=float(local.min()),
                          neighborhood_efficiency_max=float(local.max()),
                          uniformity_log_update_min=float(correction.min()),
                          uniformity_log_update_max=float(correction.max()))
        log_weights = _normalize_class_logs(log_weights, y)
        updated = np.exp(log_weights)
        record["class_weight_totals"] = [float(updated[y == label].sum()) for label in (0, 1)]
        record["positive_weight_effective_n"] = float(
            updated[y == 1].sum() ** 2 / np.dot(updated[y == 1], updated[y == 1]))
        rounds.append(record)
        if error <= 1e-12:
            stopping_reason = "perfect tree under current training weights"
            break

    if not trees:
        # A neutral score is meaningful for unsupported/no-separation training.
        stopping_reason += "; returning a neutral classifier"
    model = UniformityBDT(trees, coefficients, imputer, X.shape[1])
    history = dict(algorithm="single-working-point neighborhood-uniformity AdaBoost variant",
                   source="https://arxiv.org/abs/1305.7248", seed=int(seed), strength=float(strength),
                   target_efficiency=float(target_efficiency), requested_estimators=int(n_estimators),
                   fitted_estimators=len(trees), learning_rate=rate, max_depth=2,
                   min_samples_leaf=min_leaf, class_totals="0.5 each, renormalized after each update",
                   local_efficiency_weights="fixed original within-class input weights",
                   uniformity_log_update_clip=[-1.0, 1.0], neighborhoods=metadata,
                   training_cut_space="raw sum of weighted signed tree predictions",
                   inference_score="expit(2 * raw margin / sum(abs(tree coefficients))); neutral=0.5",
                   stopping_reason=stopping_reason, rounds=rounds,
                   requires_gen=False, lifetime_independence_validated=False,
                   flatness_guaranteed=False, deployment_ready=False)
    return model, history
