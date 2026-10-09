"""Focused synthetic checks of the weighted joint DisCo training boundary."""

from pathlib import Path
import pickle
import sys
import unittest

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from adversarial import AdversarialClassifier
from disco import train_disco, weighted_distance_correlation_squared


def small_training_arrays(seed=9):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(104, 6))
    y = (x[:, 0] + 0.5 * x[:, 1] > 0).astype(int)
    weights = rng.uniform(0.25, 2.0, size=len(y))
    z = np.column_stack((x[:, 0], x[:, 1], x[:, 2]))
    # Missing reconstructed inputs are permitted, including an empty column.
    x[:, -1] = np.nan
    return x[:80], y[:80], weights[:80], z[:80], x[80:], y[80:], weights[80:], z[80:]


class WeightedJointDistanceCorrelation(unittest.TestCase):
    def test_penalty_gradients_reach_classifier_scores(self):
        scores = torch.tensor([0.02, 0.12, 0.25, 0.35, 0.8, 0.92],
                              dtype=torch.float64, requires_grad=True)
        z = torch.tensor([[0., 0.], [1., 2.], [3., -1.], [4., 1.], [3., 3.], [5., 4.]],
                         dtype=torch.float64)
        weights = torch.tensor([1., 2., 1., 3., 2., 1.], dtype=torch.float64)
        value = weighted_distance_correlation_squared(scores, z, weights)
        value.backward()
        self.assertGreater(float(value.detach()), 0)
        self.assertTrue(torch.isfinite(scores.grad).all())
        self.assertGreater(float(scores.grad.abs().sum()), 1e-6)

    def test_joint_xor_dependence_is_missed_by_marginals(self):
        z = torch.tensor([[-1., -1.], [-1., 1.], [1., -1.], [1., 1.]],
                         dtype=torch.float64).repeat_interleave(20, dim=0)
        scores = torch.tensor([0., 1., 1., 0.], dtype=torch.float64).repeat_interleave(20)
        weights = torch.ones(len(scores), dtype=torch.float64)
        for column in (0, 1):
            marginal = weighted_distance_correlation_squared(scores, z[:, column:column + 1], weights)
            self.assertLess(float(marginal), 1e-12)
        joint = weighted_distance_correlation_squared(scores, z, weights)
        self.assertGreater(float(joint), 0.2)

    def test_weights_equal_replicated_empirical_observations(self):
        scores = torch.tensor([0.1, 0.3, 0.8, 0.6, 0.2], dtype=torch.float64)
        z = torch.tensor([[0., 2.], [1., -1.], [3., 0.], [1., 2.], [2., 4.]],
                         dtype=torch.float64)
        repeats = torch.tensor([1, 2, 3, 1, 2])
        weighted = weighted_distance_correlation_squared(scores, z, repeats.to(torch.float64))
        replicated_scores = torch.repeat_interleave(scores, repeats)
        replicated_z = torch.repeat_interleave(z, repeats, dim=0)
        replicated = weighted_distance_correlation_squared(
            replicated_scores, replicated_z, torch.ones(len(replicated_scores)))
        self.assertAlmostEqual(float(weighted), float(replicated), places=13)

    def test_missing_and_zero_weight_rows_are_excluded(self):
        scores = torch.tensor([0.1, 0.2, 0.4, 0.5, 0.6, 0.9], dtype=torch.float64)
        z = torch.tensor([[0., 1.], [1., 2.], [2., 4.], [3., -1.], [float("nan"), 0.], [10., 10.]],
                         dtype=torch.float64)
        weights = torch.tensor([1., 2., 1., 2., 3., 0.], dtype=torch.float64)
        complete = weighted_distance_correlation_squared(scores[:4], z[:4], weights[:4])
        filtered = weighted_distance_correlation_squared(scores, z, weights)
        self.assertAlmostEqual(float(complete), float(filtered), places=13)

    def test_constant_and_insufficient_support_have_finite_zero_gradients(self):
        cases = [(torch.ones(6, 3), torch.ones(6)),
                 (torch.arange(18).reshape(6, 3).to(torch.float64), torch.zeros(6)),
                 (torch.full((6, 3), float("nan")), torch.ones(6))]
        for z, weights in cases:
            scores = torch.linspace(0.1, 0.9, 6, dtype=torch.float64, requires_grad=True)
            penalty = weighted_distance_correlation_squared(scores, z, weights)
            penalty.backward()
            self.assertEqual(float(penalty.detach()), 0.)
            self.assertTrue(torch.isfinite(scores.grad).all())
            self.assertEqual(float(scores.grad.abs().sum()), 0.)
        scores = torch.full((6,), 0.5, dtype=torch.float64, requires_grad=True)
        penalty = weighted_distance_correlation_squared(
            scores, torch.arange(18).reshape(6, 3), torch.ones(6))
        penalty.backward()
        self.assertEqual(float(penalty.detach()), 0.)
        self.assertTrue(torch.isfinite(scores.grad).all())

    def test_negative_weights_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            weighted_distance_correlation_squared(torch.arange(4.).requires_grad_(),
                                                   torch.arange(12.).reshape(4, 3),
                                                   torch.tensor([1., -1., 1., 1.]))


class DiscoInferenceContract(unittest.TestCase):
    def test_zero_strength_training_is_independent_of_nuisances(self):
        arrays = small_training_arrays()
        first, _ = train_disco(*arrays, seed=17, strength=0, epochs=4)
        changed = list(arrays)
        changed[3] = np.full_like(changed[3], np.nan)
        changed[7] = np.full_like(changed[7], 1000.)
        second, _ = train_disco(*changed, seed=17, strength=0, epochs=4)
        np.testing.assert_array_equal(first.predict_proba(arrays[4]), second.predict_proba(arrays[4]))

    def test_validation_nuisances_do_not_fit_scaling_or_change_model(self):
        arrays = small_training_arrays()
        first, first_history = train_disco(*arrays, seed=23, strength=0.7, epochs=4)
        changed = list(arrays)
        changed[7] = changed[7] * 1000 + 50
        second, second_history = train_disco(*changed, seed=23, strength=0.7, epochs=4)
        self.assertEqual(first_history["nuisance_scaling"], second_history["nuisance_scaling"])
        np.testing.assert_array_equal(first.predict_proba(arrays[4]), second.predict_proba(arrays[4]))

    def test_regularizer_changes_training_and_records_joint_support(self):
        arrays = small_training_arrays()
        baseline, _ = train_disco(*arrays, seed=31, strength=0, epochs=6)
        protected, history = train_disco(*arrays, seed=31, strength=1, epochs=6)
        self.assertGreater(float(np.max(np.abs(baseline.predict_proba(arrays[4]) -
                                               protected.predict_proba(arrays[4])))), 1e-5)
        self.assertEqual(history["architecture"], [6, 32, 16, 1])
        self.assertEqual(history["nuisance_scaling"]["fitted_complete_train_rows"], 80)
        self.assertEqual(len(history["support"]), 2)
        last = history["epochs"][-1]
        self.assertTrue(np.isfinite(last["val_auc"]))
        self.assertTrue(np.isfinite(last["train_bce"]))
        self.assertTrue(np.isfinite(last["val_bce"]))
        for entry in last["decorrelation"]:
            self.assertGreater(entry["updates_with_penalty"], 0)
            self.assertLessEqual(entry["train"]["sample_rows"], 256)
            self.assertIsNotNone(entry["train"]["dcor_squared"])

    def test_pickle_roundtrip_inference_needs_x_only(self):
        arrays = small_training_arrays()
        model, _ = train_disco(*arrays, seed=2, strength=0.5, epochs=3)
        self.assertIsInstance(model, AdversarialClassifier)
        before = model.predict_proba(arrays[4])
        restored = pickle.loads(pickle.dumps(model))
        np.testing.assert_array_equal(before, restored.predict_proba(arrays[4]))
        self.assertEqual(restored.predict_proba(np.empty((0, 6))).shape, (0, 2))
        self.assertTrue(np.isfinite(before).all())
        np.testing.assert_allclose(before.sum(axis=1), 1., atol=1e-7)

    def test_missing_nuisance_targets_are_reported_without_imputation(self):
        arrays = list(small_training_arrays())
        arrays[3][:, 2] = np.nan
        model, history = train_disco(*arrays, seed=4, strength=1, epochs=3)
        self.assertIsNone(history["nuisance_scaling"])
        for entry in history["support"]:
            self.assertEqual(entry["train"]["valid_rows"], 0)
        for entry in history["epochs"][-1]["decorrelation"]:
            self.assertIsNone(entry["train"]["dcor_squared"])
            self.assertIn("fewer than four", entry["train"]["skip_reason"])
        self.assertTrue(np.isfinite(model.predict_proba(arrays[4])).all())


if __name__ == "__main__":
    unittest.main()
