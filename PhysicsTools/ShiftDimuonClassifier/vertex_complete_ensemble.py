# Code style: AI-code. AI-generated reconstructed-only equal-weight ensemble.
"""Average fixed full-information fits; no labels or nuisance targets in inference."""
import numpy as np
from vertex_complete_bdt import COMPLETE_FEATURE_NAMES, CATEGORY_NAMES


class CompleteVertexEnsemble:
    def __init__(self, models, explicit_category):
        if not models or len(models) != len(explicit_category):
            raise ValueError('Aligned nonempty ensemble members are required')
        if any(model.n_features_in_ != 119 for model in models):
            raise ValueError('Every member must retain all119 reconstructed features')
        self.models = list(models)
        self.explicit_category = list(explicit_category)
        self.n_features_in_ = 119
        self.feature_names_in_ = np.asarray(COMPLETE_FEATURE_NAMES)
        self.categories_ = CATEGORY_NAMES
        self.classes_ = np.asarray([0, 1])
        self.requires_gen = False

    def predict_proba(self, X, categories):
        x, category = np.asarray(X), np.asarray(categories)
        if x.ndim != 2 or x.shape[1] != 119 or category.shape != (len(x),):
            raise ValueError('Expected aligned RECO119 and reconstructed categories')
        if not set(category).issubset(CATEGORY_NAMES):
            raise ValueError('Unknown reconstructed category')
        if not len(x):
            return np.empty((0, 2))
        predictions = []
        for model, explicit in zip(self.models, self.explicit_category):
            values = model.predict_proba(x, category) if explicit else model.predict_proba(x)
            if values.shape != (len(x), 2) or not np.isfinite(values).all():
                raise ValueError('Invalid member probabilities')
            predictions.append(np.asarray(values, dtype=float))
        result = np.mean(predictions, axis=0)
        if np.any((result < 0) | (result > 1)) or not np.allclose(result.sum(axis=1), 1):
            raise ValueError('Ensemble probabilities do not normalize')
        return result
