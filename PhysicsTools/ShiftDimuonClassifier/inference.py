"""Score only the reconstructed feature matrix; labels and nuisances are absent."""
import numpy as np
from features import FEATURE_NAMES


class DistilledBDT:
    """Regression trees approximating a teacher; decorrelation needs retesting."""
    def __init__(self, estimator):
        self.estimator = estimator

    def predict_proba(self, X):
        p = np.clip(self.estimator.predict(X), 0, 1)
        return np.column_stack((1 - p, p))


def predict(bundle, X, feature_names=FEATURE_NAMES):
    if tuple(feature_names) != tuple(bundle["feature_names"]) or tuple(feature_names) != FEATURE_NAMES:
        raise ValueError("Feature contract differs from trained model")
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[1] != len(FEATURE_NAMES):
        raise ValueError("Wrong reconstructed feature matrix shape")
    if not len(X):
        return np.empty(0)
    X = np.where(np.isfinite(X), X, np.nan)
    return bundle["model"].predict_proba(X[:, bundle["columns"]])[:, 1]
