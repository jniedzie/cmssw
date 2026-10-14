# Code style: AI-code. AI-generated fixed-category neural adapter tests.
from pathlib import Path
import pickle
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import neural_vertex_categories as legacy
from vertex_complete_neural import fixed_category_arrays, train_category_vertex
from vertex_complete_bdt import CATEGORY_NAMES, COMPLETE_FEATURE_NAMES


def fixture():
    rng = np.random.default_rng(17)
    x = rng.normal(size=(96, 119))
    x[:, -1] = np.nan
    y = np.tile([0, 1], 48)
    w = np.ones(96)
    z = rng.normal(size=(96, 3))
    data = [x[:64], y[:64], w[:64], z[:64], x[64:], y[64:], w[64:], z[64:]]
    options = dict(feature_names=COMPLETE_FEATURE_NAMES,
                   train_category=np.full(64, "Both-Both"),
                   val_category=np.tile(["Both-Both", "Other"], 16),
                   strength=0., category_strength=0., conditional_strength=0.,
                   acceptance_strength=0., hidden=(8,), dropout=0., epochs=2,
                   batch_size=32, patience=2, seed=18)
    return data, options


class CompleteNeuralTests(unittest.TestCase):
    def test_absent_training_category_is_accepted_in_validation_by_fixed_schema(self):
        data, options = fixture()
        original_encoder = legacy._category_arrays
        with self.assertRaisesRegex(ValueError, "unseen validation"):
            original_encoder(options["train_category"], options["val_category"], 64, 32)
        model, history = train_category_vertex(*data, **options, architecture="shared")
        self.assertIs(legacy._category_arrays, original_encoder)
        self.assertEqual(tuple(model.categories_), CATEGORY_NAMES)
        self.assertEqual(history["training_categories"], list(CATEGORY_NAMES))
        self.assertEqual(history["fixed_reco_schema"]["missing_training_categories"], list(CATEGORY_NAMES[1:]))
        self.assertEqual(history["category_support"]["Other"]["train"]["genuine"], 0)
        self.assertGreater(history["category_support"]["Other"]["val"]["genuine"], 0)
        self.assertEqual(model.network.trunk[0].in_features, 124)
        self.assertTrue(np.isfinite(model.predict_proba(data[4], options["val_category"])).all())
        with self.assertRaisesRegex(ValueError, "unseen validation"):
            legacy._category_arrays(options["train_category"], options["val_category"], 64, 32)

    def test_fixed_heads_pickle_reload_and_empty_inference(self):
        data, options = fixture()
        model, history = train_category_vertex(*data, **options, architecture="heads")
        self.assertEqual(model.network.heads.out_features, 5)
        self.assertEqual(model.n_features_in_, 119)
        self.assertEqual(history["epochs_completed"], 2)
        self.assertEqual(history["best_state_sha256"], history["restored_state_sha256"])
        reloaded = pickle.loads(pickle.dumps(model))
        np.testing.assert_array_equal(model.predict_proba(data[4], options["val_category"]),
                                      reloaded.predict_proba(data[4], options["val_category"]))
        self.assertEqual(model.predict_proba(np.empty((0, 119)), np.array([], dtype="U16")).shape, (0, 2))
        self.assertFalse(model.requires_gen)
        for forbidden in ("Z_train", "labels", "process", "proper_lifetime", "training_ids"):
            self.assertFalse(hasattr(model, forbidden))

    def test_original_encoder_restored_on_training_failure(self):
        data, options = fixture()
        original_encoder = legacy._category_arrays
        options["dropout"] = 2.0
        with self.assertRaisesRegex(ValueError, "dropout"):
            train_category_vertex(*data, **options)
        self.assertIs(legacy._category_arrays, original_encoder)

    def test_unknown_future_prediction_category_rejected(self):
        data, options = fixture()
        model, _ = train_category_vertex(*data, **options)
        with self.assertRaisesRegex(ValueError, "Unseen reconstructed"):
            model.predict_proba(data[4][:1], np.array(["new_category"]))
        prediction = model.predict_proba(data[4][:5], np.asarray(CATEGORY_NAMES))
        self.assertEqual(prediction.shape, (5, 2))
        self.assertTrue(np.isfinite(prediction).all())

    def test_train_validation_schema_and_features_fail_closed(self):
        ct, cv, categories = fixed_category_arrays(["Both-Both"], ["Other"], 1, 1)
        self.assertEqual(categories, list(CATEGORY_NAMES))
        self.assertEqual(ct.shape, (1, 5)); self.assertEqual(cv.shape, (1, 5))
        np.testing.assert_array_equal(ct.sum(axis=1), [1]); np.testing.assert_array_equal(cv.sum(axis=1), [1])
        with self.assertRaisesRegex(ValueError, "Unknown reconstructed"):
            fixed_category_arrays(["inclusive"], ["Other"], 1, 1)
        data, options = fixture()
        bad = list(COMPLETE_FEATURE_NAMES); bad[0] = "GenPart_ctau"
        with self.assertRaisesRegex(ValueError, "119 reconstructed"):
            train_category_vertex(*data, **{**options, "feature_names": bad})
        with self.assertRaisesRegex(ValueError, "shared or heads"):
            train_category_vertex(*data, **options, architecture="linear_experts")


if __name__ == "__main__":
    unittest.main()
