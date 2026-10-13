# Code style: AI-code. AI-generated tests for neural training boundaries.
"""Synthetic tests; no detector data, truth branches, ROOT or production jobs."""
from pathlib import Path
import pickle
import sys
import unittest

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from neural_vertex import _category_penalty, _groups, _penalty, train_neural_vertex


def arrays(seed=13, features=6):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(180, features))
    y = (x[:, 0] + .4 * x[:, 1] > 0).astype(int)
    w = rng.uniform(.2, 2, len(y))
    z = x[:, :3].copy()
    x[:, -1] = np.nan
    return [x[:120], y[:120], w[:120], z[:120], x[120:], y[120:], w[120:], z[120:]]


class NeuralVertexTraining(unittest.TestCase):
    def test_disabled_categories_cannot_affect_optimization(self):
        data = arrays()
        first, h1 = train_neural_vertex(*data, epochs=3, batch_size=32)
        # Even invalid/misaligned metadata is ignored when the toggle is off.
        second, h2 = train_neural_vertex(*data, epochs=3, batch_size=32,
                                         category_strength=0, train_category=["A"], val_category=["unseen"])
        np.testing.assert_array_equal(first.predict_proba(data[4]), second.predict_proba(data[4]))
        self.assertEqual(h1["best_state_sha256"], h2["best_state_sha256"])

    def test_genuine_category_gradient_and_undefined_single_category(self):
        scores = torch.tensor([.1, .3, .6, .2, .7, .9, .8, .5], dtype=torch.float64, requires_grad=True)
        category = torch.eye(2, dtype=torch.float64)[torch.tensor([0, 0, 0, 0, 1, 1, 1, 1])]
        weights = torch.ones(8, dtype=torch.float64)
        penalty, detail = _category_penalty(scores, category, weights, np.arange(8))
        penalty.backward()
        self.assertGreater(float(penalty.detach()), .1)
        self.assertGreater(float(scores.grad.abs().sum()), 1e-6)
        self.assertTrue(torch.isfinite(scores.grad).all())
        self.assertEqual(detail["populated_categories"], 2)
        value, detail = _category_penalty(scores, category, weights, np.arange(4))
        self.assertEqual(float(value.detach()), 0)
        self.assertIsNone(detail["dcor_squared"])
        self.assertIn("two populated", detail["skip_reason"])

    def test_category_training_encoder_is_train_only_and_not_in_model(self):
        data = arrays()
        categories = np.where(np.arange(180) % 2, "BothBoth", "NearBoth")
        kwargs = dict(epochs=3, batch_size=32, train_category=categories[:120], val_category=categories[120:])
        first, _ = train_neural_vertex(*data, **kwargs)
        second, history = train_neural_vertex(*data, category_strength=.3, **kwargs)
        self.assertGreater(float(np.max(np.abs(first.predict_proba(data[4]) - second.predict_proba(data[4])))), 1e-6)
        self.assertEqual(history["training_categories"], ["BothBoth", "NearBoth"])
        record = history["epochs"][0]
        self.assertGreater(record["val_category_dcor_squared"], 0)
        self.assertLessEqual(record["val_category_diagnostic"]["rows"], int((data[5] == 1).sum()))
        self.assertFalse(hasattr(second, "training_categories"))
        self.assertFalse(hasattr(second, "category_encoder"))
        kwargs["val_category"] = np.full(60, "unseen")
        with self.assertRaisesRegex(ValueError, "unseen validation"):
            train_neural_vertex(*data, category_strength=.3, **kwargs)
        kwargs["train_category"] = np.full(120, "BothBoth")
        kwargs["val_category"] = np.full(60, "BothBoth")
        _, single = train_neural_vertex(*data, category_strength=.3, **kwargs)
        self.assertIsNone(single["epochs"][0]["val_category_dcor_squared"])

    def test_zero_penalty_ignores_nuisances_and_process(self):
        data = arrays()
        first, h1 = train_neural_vertex(*data, strength=0, epochs=5, batch_size=32)
        changed = [x.copy() for x in data]
        changed[3][:] = np.nan
        changed[7] *= 10000
        second, h2 = train_neural_vertex(*changed, strength=0, epochs=5, batch_size=32,
                                         train_process=np.arange(120), val_process=np.arange(60))
        np.testing.assert_array_equal(first.predict_proba(data[4]), second.predict_proba(data[4]))
        self.assertEqual(h1["best_state_sha256"], h2["best_state_sha256"])

    def test_restores_best_epoch_after_validation_deterioration(self):
        data = arrays()
        # Opposite validation labels make learning useful train structure hurt validation.
        data[5] = 1 - data[5]
        model, history = train_neural_vertex(*data, strength=0, epochs=30, patience=4,
                                             learning_rate=.02, dropout=0, batch_size=32)
        objectives = [record["val_objective"] for record in history["epochs"]]
        self.assertEqual(history["best_epoch"], int(np.argmin(objectives)) + 1)
        self.assertLess(history["best_epoch"], history["epochs_completed"])
        self.assertLess(history["epochs_completed"], 30)
        self.assertEqual(history["best_state_sha256"], history["restored_state_sha256"])
        self.assertNotEqual(history["last_state_sha256"], history["restored_state_sha256"])
        score = np.clip(model.predict_proba(data[4])[:, 1].astype(float), 1e-12, 1-1e-12)
        bce = np.average(-data[5] * np.log(score) - (1-data[5]) * np.log1p(-score), weights=data[6])
        self.assertAlmostEqual(bce, history["best_val_objective"], places=6)

    def test_preprocessing_fits_train_only_and_pickle_needs_x_only(self):
        data = arrays()
        model, history = train_neural_vertex(*data, epochs=3, batch_size=32)
        changed = [x.copy() for x in data]
        changed[4][:, :-1] += 1000
        changed[7] *= 10000
        other, other_history = train_neural_vertex(*changed, epochs=3, batch_size=32)
        np.testing.assert_array_equal(model.imputer.statistics_, other.imputer.statistics_)
        np.testing.assert_array_equal(model.scaler.mean_, other.scaler.mean_)
        self.assertEqual(history["nuisance_scaling"], other_history["nuisance_scaling"])
        restored = pickle.loads(pickle.dumps(model))
        np.testing.assert_array_equal(model.predict_proba(data[4]), restored.predict_proba(data[4]))
        self.assertEqual(model.predict_proba(np.empty((0, 6))).shape, (0, 2))
        self.assertFalse(hasattr(model, "train_process"))
        self.assertFalse(hasattr(model, "nuisance_scaling"))

    def test_process_penalties_are_metadata_only_and_record_support(self):
        data = arrays()
        pt = np.where(np.arange(120) % 2, "DY", "Jpsi")
        pv = np.where(np.arange(60) % 2, "DY", "Jpsi")
        model, history = train_neural_vertex(*data, epochs=3, batch_size=32,
                                             penalty_scope="class_process", train_process=pt, val_process=pv)
        groups = history["epochs"][0]["val_groups"]
        self.assertEqual({g["group"] for g in groups},
                         {"class_0", "class_1", "genuine_process_DY", "genuine_process_Jpsi"})
        self.assertTrue(all(g["rows"] >= 4 for g in groups))
        self.assertTrue(np.isfinite(model.predict_proba(data[4])).all())
        # Relabeling process names does not change the partition or training.
        other, _ = train_neural_vertex(*data, epochs=3, batch_size=32,
                                       penalty_scope="class_process", train_process=np.where(pt == "DY", "A", "B"),
                                       val_process=np.where(pv == "DY", "A", "B"))
        np.testing.assert_array_equal(model.predict_proba(data[4]), other.predict_proba(data[4]))

    def test_joint_penalty_gradient_and_class_separation(self):
        z = torch.tensor([[-1., -1., 0], [-1., 1., 0], [1., -1., 0], [1., 1., 0]],
                         dtype=torch.float64).repeat(8, 1)
        score = torch.tensor([.1, .8, .7, .2], dtype=torch.float64).repeat(8).requires_grad_()
        y = np.repeat([0, 1], 16)
        weights = torch.ones(32, dtype=torch.float64)
        groups = _groups(y, weights.numpy(), z.numpy())
        penalty, details = _penalty(score, z, weights, groups)
        penalty.backward()
        self.assertEqual(len(details), 2)
        self.assertGreater(float(penalty.detach()), .2)
        self.assertTrue(torch.isfinite(score.grad).all())
        self.assertGreater(float(score.grad.abs().sum()), 1e-6)
        changed = score.detach().clone()
        changed[16:] = .5
        _, changed_details = _penalty(changed, z, weights, groups)
        self.assertEqual(details[0]["dcor_squared"], changed_details[0]["dcor_squared"])
        self.assertEqual(changed_details[1]["dcor_squared"], 0)

    def test_energy_uses_reconstructed_columns_and_train_scaling(self):
        data = arrays(features=50)
        model, history = train_neural_vertex(*data, epochs=2, batch_size=32, energy_strength=.3)
        self.assertIsNotNone(history["energy_scaling"])
        self.assertEqual(len(history["energy_scaling"]["mean"]), 2)
        self.assertEqual(len(history["epochs"][0]["val_energy_groups"]), 2)
        self.assertTrue(np.isfinite(model.predict_proba(data[4])).all())
        with self.assertRaisesRegex(ValueError, "50-feature"):
            train_neural_vertex(*arrays(), epochs=1, energy_strength=.3)

    def test_bad_weights_and_unknown_training_labels_fail_closed(self):
        for position, value in ((2, -1), (2, np.nan), (6, -1), (1, -1)):
            data = arrays()
            data[position][0] = value
            with self.assertRaises(ValueError):
                train_neural_vertex(*data, epochs=1)
        data = arrays()
        data[2][:] = 0
        with self.assertRaises(ValueError):
            train_neural_vertex(*data, epochs=1)


if __name__ == "__main__":
    unittest.main()
