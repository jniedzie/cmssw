"""Uniformity-update direction, training boundaries and data-only inference."""

import inspect
import json
import pickle
from pathlib import Path
import sys
import unittest

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.tree import DecisionTreeClassifier

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uniform_bdt import (
    UniformityBDT, _local_efficiencies, _uniformity_log_update, train_uniform_bdt,
)


def synthetic_data(n=120, seed=19):
    rng = np.random.default_rng(seed)
    y = np.repeat([0, 1], n)
    nuisance = rng.uniform(-1, 1, len(y))
    X = np.column_stack((rng.normal(y, 0.9), nuisance + rng.normal(0, 0.2, len(y))))
    Z = np.column_stack((nuisance, nuisance ** 2, np.sin(nuisance)))
    return X, y, np.ones(len(y)), Z


class UniformityMechanism(unittest.TestCase):
    def test_underaccepted_neighborhood_gets_more_attention(self):
        # Two supported neighborhoods have 0% and 100% true-pair acceptance.
        passed = np.array([False, False, True, True])
        neighbors = np.array([[0, 1], [0, 1], [2, 3], [2, 3]])
        weights = np.array([1., 3., 2., 2.])
        local = _local_efficiencies(passed, neighbors, weights)
        update = _uniformity_log_update(local, 0.5, strength=2., rate=0.25)
        self.assertTrue(np.all(update[:2] > 0))
        self.assertTrue(np.all(update[2:] < 0))
        # The global normalization cannot reverse the relative attention shift.
        new_weights = weights * np.exp(update)
        self.assertGreater(new_weights[:2].sum() / new_weights[2:].sum(),
                           weights[:2].sum() / weights[2:].sum())
        np.testing.assert_array_equal(_uniformity_log_update(local, 0.5, 0., .25), 0.)

    def test_fixed_input_weights_define_local_acceptance(self):
        neighbors = np.array([[0, 1], [0, 1]])
        efficiency = _local_efficiencies(np.array([False, True]), neighbors, np.array([1., 3.]))
        np.testing.assert_array_equal(efficiency, np.array([.75, .75]))

    def test_update_changes_fitted_trees_and_preserves_class_totals(self):
        # Process-independent labels; a nuisance-correlated accidental population
        # gives the classifier an easy but nonuniform discrimination shortcut.
        rng = np.random.default_rng(892)
        y = np.repeat([0, 1], 500)
        nuisance = rng.uniform(-1, 1, len(y))
        nuisance[y == 0] = rng.uniform(-1, 0, np.count_nonzero(y == 0))
        X = np.column_stack((rng.normal(y, 1.), nuisance))
        Z = np.column_stack((nuisance, nuisance ** 2, np.sin(nuisance)))
        ordinary, _ = train_uniform_bdt(X, y, np.ones(len(y)), Z, strength=0., n_estimators=60)
        uniform, history = train_uniform_bdt(X, y, np.ones(len(y)), Z, strength=3., n_estimators=60)
        self.assertGreater(np.max(np.abs(ordinary.predict_proba(X) - uniform.predict_proba(X))), 1e-4)
        self.assertTrue(history["neighborhoods"]["active"])
        self.assertTrue(any(r["uniformity_log_update_max"] > 0 for r in history["rounds"]))
        for record in history["rounds"]:
            np.testing.assert_allclose(record["class_weight_totals"], [.5, .5], atol=1e-14)
        # This checks the mechanism, not a guaranteed held-out flatness result.
        self.assertFalse(history["flatness_guaranteed"])
        self.assertFalse(history["lifetime_independence_validated"])
        json.dumps(history, allow_nan=False)


class UniformityBoundaries(unittest.TestCase):
    def test_strength_zero_ignores_nuisance_information(self):
        X, y, w, Z = synthetic_data()
        reference, _ = train_uniform_bdt(X, y, w, Z, strength=0., n_estimators=20)
        changed, _ = train_uniform_bdt(X, y, w, Z[::-1] * 100., strength=0., n_estimators=20)
        missing, _ = train_uniform_bdt(X, y, w, np.full_like(Z, np.nan), strength=0., n_estimators=20)
        np.testing.assert_array_equal(reference.predict_proba(X), changed.predict_proba(X))
        np.testing.assert_array_equal(reference.predict_proba(X), missing.predict_proba(X))

    def test_training_only_imputer_and_nuisance_scale(self):
        X, y, w, Z = synthetic_data()
        X[0, 0] = np.nan
        model, history = train_uniform_bdt(X, y, w, Z, n_estimators=20)
        expected_medians = np.nanmedian(X, axis=0)
        np.testing.assert_allclose(model.imputer.statistics_, expected_medians)
        np.testing.assert_allclose(history["neighborhoods"]["center"], Z[y == 1].mean(axis=0))
        np.testing.assert_allclose(history["neighborhoods"]["scale"], Z[y == 1].std(axis=0))
        model.predict_proba(np.array([[1e10, -1e10], [np.nan, np.inf]]))
        np.testing.assert_allclose(model.imputer.statistics_, expected_medians)

    def test_zero_weight_rows_cannot_change_imputer(self):
        X, y, w, Z = synthetic_data()
        model, _ = train_uniform_bdt(X, y, w, Z, strength=0., n_estimators=5)
        expanded, _ = train_uniform_bdt(np.vstack((X, [1e10, 1e10])), np.append(y, 1),
                                         np.append(w, 0.), np.vstack((Z, [1e10] * 3)),
                                         strength=0., n_estimators=5)
        np.testing.assert_array_equal(model.imputer.statistics_, expanded.imputer.statistics_)

    def test_positive_unit_changes_do_not_change_training_neighborhoods(self):
        X, y, w, Z = synthetic_data()
        model, history = train_uniform_bdt(X, y, w, Z, strength=2., n_estimators=20)
        changed, other = train_uniform_bdt(X, y, w, Z * [10., 1000., 0.25] + [4., -7., 3.],
                                           strength=2., n_estimators=20)
        np.testing.assert_allclose(model.predict_proba(X), changed.predict_proba(X), atol=1e-14)
        self.assertEqual(history["neighborhoods"]["n_neighbors"], other["neighborhoods"]["n_neighbors"])

    def test_missing_targets_are_counted_and_not_imputed(self):
        X, y, w, Z = synthetic_data()
        Z[y == 1, :] = np.nan
        model, history = train_uniform_bdt(X, y, w, Z, strength=3., n_estimators=20)
        ordinary, _ = train_uniform_bdt(X, y, w, Z, strength=0., n_estimators=20)
        self.assertFalse(history["neighborhoods"]["active"])
        self.assertEqual(history["neighborhoods"]["missing_positive_rows"], int((y == 1).sum()))
        self.assertIsNone(history["neighborhoods"]["center"])
        np.testing.assert_array_equal(model.predict_proba(X), ordinary.predict_proba(X))

    def test_unknown_labels_and_invalid_inputs_fail(self):
        X, y, w, Z = synthetic_data()
        for overrides in ({"y_train": np.full_like(y, -1)},
                          {"w_train": -w}, {"w_train": np.zeros_like(w)},
                          {"Z_train": Z[:, :2]}, {"strength": -1.},
                          {"n_estimators": 0}, {"target_efficiency": 1.}):
            with self.subTest(overrides=list(overrides)):
                arguments = dict(X_train=X, y_train=y, w_train=w, Z_train=Z, n_estimators=5)
                arguments.update(overrides)
                with self.assertRaises(ValueError):
                    train_uniform_bdt(**arguments)

    def test_pickle_and_inference_need_only_reconstructed_features(self):
        X, y, w, Z = synthetic_data()
        X[:, 0][::7] = np.nan
        model, _ = train_uniform_bdt(X, y, w, Z, n_estimators=20)
        expected = model.predict_proba(X)
        restored = pickle.loads(pickle.dumps(model))
        np.testing.assert_array_equal(restored.predict_proba(X), expected)
        np.testing.assert_allclose(expected.sum(axis=1), 1.)
        self.assertEqual(list(inspect.signature(restored.predict_proba).parameters), ["X"])
        self.assertEqual(set(restored.__dict__), {"estimators_", "estimator_weights_", "imputer",
                                                "n_features_in_", "classes_"})
        self.assertEqual(restored.predict_proba(np.empty((0, X.shape[1]))).shape, (0, 2))
        with self.assertRaises(ValueError):
            restored.predict_proba(X[:, :1])

    def test_large_margin_scores_preserve_order_without_sigmoid_saturation(self):
        X = np.arange(6., dtype=np.float64).reshape(-1, 1)
        left = DecisionTreeClassifier(max_depth=1, random_state=1).fit(X, np.array([0, 1, 1, 1, 1, 1]))
        right = DecisionTreeClassifier(max_depth=1, random_state=1).fit(X, np.array([0, 0, 0, 0, 0, 1]))
        model = UniformityBDT([left, right], [1e6, 1e6], SimpleImputer().fit(X), 1)
        self.assertGreater(np.max(np.abs(model.decision_function(X))), 1000.)
        probabilities = model.predict_proba(X)[:, 1]
        self.assertTrue(np.all((probabilities > 0) & (probabilities < 1)))
        self.assertLess(probabilities[0], probabilities[1])
        self.assertLess(probabilities[1], probabilities[-1])
        self.assertEqual(len(np.unique(probabilities)), 3)


if __name__ == "__main__":
    unittest.main()
