# Code style: AI-code. AI-generated synthetic optimization-schedule checks.
import inspect
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from features import FEATURE_NAMES
from neural_vertex import NeuralVertexClassifier, train_neural_vertex as original
from neural_vertex_scheduled import train_neural_vertex as scheduled


def arrays():
    rng = np.random.default_rng(12)
    x = rng.normal(size=(120, len(FEATURE_NAMES)))
    y = (x[:, 0] + .2 * x[:, 1] > 0).astype(int)
    w = rng.uniform(.4, 1.7, len(x))
    z = x[:, :3].copy()
    data = (x[:80], y[:80], w[:80], z[:80], x[80:], y[80:], w[80:], z[80:])
    category = np.where(np.arange(120) % 2, 'NearBoth', 'BothBoth')
    process = np.where(np.arange(120) % 3, 'dy', 'jpsi')
    metadata = dict(train_category=category[:80], val_category=category[80:],
                    train_process=process[:80], val_process=process[80:])
    return data, metadata


class ScheduledNeuralTests(unittest.TestCase):
    def test_zero_warmup_bitwise_matches_frozen_trainer_with_all_penalties(self):
        data, metadata = arrays()
        options = dict(epochs=4, batch_size=32, hidden=(8, 4), strength=.3,
                       category_strength=.4, energy_strength=.2, penalty_scope='class_process', **metadata)
        old, old_history = original(*data, **options)
        new, new_history = scheduled(*data, warmup_epochs=0, **options)
        np.testing.assert_array_equal(old.predict_proba(data[4]), new.predict_proba(data[4]))
        self.assertEqual(old_history['best_state_sha256'], new_history['best_state_sha256'])
        for before, after in zip(old_history['epochs'], new_history['epochs']):
            self.assertEqual(before, {key: after[key] for key in before})
        self.assertIsInstance(new, NeuralVertexClassifier)

    def test_ramped_optimization_keeps_full_validation_objective_every_epoch(self):
        data, metadata = arrays()
        _, history = scheduled(*data, epochs=4, warmup_epochs=2, batch_size=32,
            hidden=(8,), strength=.3, energy_strength=.2, category_strength=.4, **metadata)
        self.assertEqual([r['optimization_penalty_fraction'] for r in history['epochs']], [0., .5, 1., 1.])
        for row in history['epochs']:
            expected = (row['val_bce'] + .3 * row['val_dcor_squared']
                        + .2 * row['val_energy_dcor_squared'] + .4 * (row['val_category_dcor_squared'] or 0.))
            self.assertAlmostEqual(row['val_objective'], expected, places=14)
            self.assertEqual(row['optimization_joint_strength'], .3 * row['optimization_penalty_fraction'])
            self.assertGreaterEqual(row['val_score_weighted_std'], 0.)
        self.assertGreater(history['epochs'][0]['val_objective'], history['epochs'][0]['val_bce'])
        self.assertEqual(history['best_val_objective'], min(r['val_objective'] for r in history['epochs']))

    def test_checkpoint_restoration_after_late_validation_deterioration(self):
        data, _ = arrays()
        data = list(data)
        data[5] = 1 - data[5]
        model, history = scheduled(*data, epochs=30, patience=4, warmup_epochs=3,
                                  strength=0., learning_rate=.03, dropout=0., batch_size=32, hidden=(8,))
        self.assertLess(history['best_epoch'], history['epochs_completed'])
        self.assertEqual(history['best_state_sha256'], history['restored_state_sha256'])
        self.assertNotEqual(history['last_state_sha256'], history['restored_state_sha256'])
        self.assertTrue(history['restored_prediction_exact'])
        prediction = np.clip(model.predict_proba(data[4])[:, 1].astype(float), 1e-12, 1-1e-12)
        bce = np.average(-data[5] * np.log(prediction) - (1-data[5]) * np.log1p(-prediction), weights=data[6])
        self.assertAlmostEqual(bce, history['best_val_objective'], places=6)

    def test_preprocessing_stays_train_only_and_api_has_no_test_or_truth(self):
        data, _ = arrays()
        model, history = scheduled(*data, epochs=2, batch_size=32, hidden=(8,))
        changed = [value.copy() for value in data]
        changed[4] += 1000
        changed[7] *= 10000
        other, other_history = scheduled(*changed, epochs=2, batch_size=32, hidden=(8,))
        np.testing.assert_array_equal(model.imputer.statistics_, other.imputer.statistics_)
        np.testing.assert_array_equal(model.scaler.mean_, other.scaler.mean_)
        self.assertEqual(history['nuisance_scaling'], other_history['nuisance_scaling'])
        signature = inspect.signature(scheduled)
        self.assertFalse(any('test' in name.lower() or 'gen' in name.lower() or 'truth' in name.lower()
                             for name in signature.parameters))
        self.assertFalse(any(hasattr(model, name) for name in ('train_category', 'train_process', 'nuisance_scaling')))

    def test_invalid_warmup_rejected(self):
        data, _ = arrays()
        for value in (-1, .5, True):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'warmup_epochs'):
                scheduled(*data, warmup_epochs=value)


if __name__ == '__main__':
    unittest.main()
