# Code style: AI-code. AI-generated complete-RECO BDT behavior tests.
from pathlib import Path
import pickle
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from vertex_complete_bdt import (
    BDT_CONFIGS, CATEGORY_NAMES, COMPLETE_FEATURE_NAMES, fit_complete_bdt,
)


def fixture():
    rng = np.random.default_rng(73)
    labels = np.tile([0, 1], 160)
    # No X information: only the reconstructed category distinguishes classes.
    values = np.zeros((320, 119))
    values[:, -1] = np.nan
    topology = np.where(labels == 1, "Both-Both", "Near-Both")
    weights = rng.uniform(.5, 1.5, len(labels))
    return values, topology, labels, weights


class CompleteBDTTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x, cls.topology, cls.y, cls.weights = fixture()
        cls.model = fit_complete_bdt(cls.x, cls.topology, cls.y, cls.weights,
                                    feature_names=COMPLETE_FEATURE_NAMES, seed=19)

    def test_category_information_is_used_without_other_predictors(self):
        score = self.model.predict_proba(self.x, self.topology)[:, 1]
        self.assertGreater(float(np.min(score[self.y == 1])), .95)
        self.assertLess(float(np.max(score[self.y == 0])), .05)
        changed = np.where(self.y == 1, "Near-Both", "Both-Both")
        swapped = self.model.predict_proba(self.x, changed)[:, 1]
        self.assertGreater(float(np.mean(score[self.y == 1] - swapped[self.y == 1])), .9)
        self.assertEqual(self.model.estimator.n_features_in_, 124)

    def test_fixed_iterations_no_internal_split_and_pickle_reload(self):
        self.assertFalse(self.model.estimator.early_stopping)
        self.assertEqual(self.model.estimator.n_iter_, 200)
        self.assertFalse(self.model.training_metadata_["validation_or_test_used"])
        restored = pickle.loads(pickle.dumps(self.model))
        np.testing.assert_array_equal(self.model.predict_proba(self.x, self.topology),
                                      restored.predict_proba(self.x, self.topology))
        self.assertFalse(self.model.requires_gen)
        for forbidden in ("labels", "nuisance", "process", "training_ids", "Z_train"):
            self.assertFalse(hasattr(self.model, forbidden))

    def test_zero_weight_rows_cannot_change_imputer_or_tree_bins(self):
        extra = np.full((1, 119), 1e90)
        augmented = fit_complete_bdt(np.vstack((self.x, extra)), np.concatenate((self.topology, ["Other"])),
                                    np.r_[self.y, 1], np.r_[self.weights, 0.],
                                    feature_names=COMPLETE_FEATURE_NAMES, seed=19)
        np.testing.assert_array_equal(self.model.imputer.statistics_, augmented.imputer.statistics_)
        np.testing.assert_array_equal(self.model.predict_proba(self.x, self.topology),
                                      augmented.predict_proba(self.x, self.topology))
        self.assertEqual(augmented.training_metadata_["positive_weight_rows"], len(self.x))

    def test_prediction_does_not_refit_missing_value_statistics(self):
        before = self.model.imputer.statistics_.copy()
        changed = self.x.copy(); changed[:, 0] = 1e20
        self.model.predict_proba(changed, self.topology)
        np.testing.assert_array_equal(before, self.model.imputer.statistics_)
        self.assertEqual(self.model.imputer.statistics_[-1], 0.)
        self.assertEqual(self.model.predict_proba(np.empty((0, 119)), np.array([], dtype="U16")).shape, (0, 2))

    def test_exact_reconstructed_feature_contract_and_categories_fail_closed(self):
        bad_names = list(COMPLETE_FEATURE_NAMES); bad_names[0] = "GenPart_properLifetime"
        with self.assertRaisesRegex(ValueError, "119 reconstructed"):
            fit_complete_bdt(self.x, self.topology, self.y, self.weights, feature_names=bad_names)
        swapped = list(COMPLETE_FEATURE_NAMES); swapped[0], swapped[1] = swapped[1], swapped[0]
        with self.assertRaisesRegex(ValueError, "119 reconstructed"):
            self.model.predict_proba(self.x, self.topology, feature_names=swapped)
        with self.assertRaisesRegex(ValueError, "Unknown reconstructed category"):
            self.model.predict_proba(self.x[:1], np.array(["inclusive"]))
        with self.assertRaisesRegex(ValueError, "119 columns"):
            self.model.predict_proba(self.x[:, :-1], self.topology)
        self.assertEqual(len(COMPLETE_FEATURE_NAMES), 119)
        self.assertEqual(len(set(COMPLETE_FEATURE_NAMES)), 119)
        self.assertEqual(len(CATEGORY_NAMES), 5)

    def test_corrupt_training_inputs_and_undeclared_configuration_rejected(self):
        for options in ({"config": "new_unlocked_config"}, {"config": {"max_iter": 5}},
                        {"seed": -1}, {"seed": True}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                fit_complete_bdt(self.x, self.topology, self.y, self.weights,
                                 feature_names=COMPLETE_FEATURE_NAMES, **options)
        for weights in (np.full(len(self.y), np.nan), -np.ones(len(self.y)), np.zeros(len(self.y))):
            with self.subTest(weights=weights[:2]), self.assertRaises(ValueError):
                fit_complete_bdt(self.x, self.topology, self.y, weights, feature_names=COMPLETE_FEATURE_NAMES)
        changed = self.y.copy(); changed[0] = -1
        with self.assertRaisesRegex(ValueError, "integer 0/1"):
            fit_complete_bdt(self.x, self.topology, changed, self.weights, feature_names=COMPLETE_FEATURE_NAMES)

    def test_settings_menu_is_predeclared_and_distinct(self):
        self.assertEqual(set(BDT_CONFIGS), {"shallow_d3", "medium_d5", "regularized_d5"})
        self.assertLess(BDT_CONFIGS["shallow_d3"]["max_depth"], BDT_CONFIGS["medium_d5"]["max_depth"])
        self.assertGreater(BDT_CONFIGS["regularized_d5"]["min_samples_leaf"], BDT_CONFIGS["medium_d5"]["min_samples_leaf"])
        self.assertGreater(BDT_CONFIGS["regularized_d5"]["l2_regularization"], BDT_CONFIGS["medium_d5"]["l2_regularization"])

    def test_last_reconstructed_feature_is_used_with_constant_category(self):
        x = np.zeros_like(self.x)
        x[:, 118] = self.y
        topology = np.full(len(self.y), "Both-Both")
        for config in ("medium_d5", "regularized_d5"):
            with self.subTest(config=config):
                model = fit_complete_bdt(x, topology, self.y, self.weights,
                                         feature_names=COMPLETE_FEATURE_NAMES, config=config, seed=19)
                score = model.predict_proba(x, topology)[:, 1]
                self.assertGreater(float(np.min(score[self.y == 1])), .9)
                self.assertLess(float(np.max(score[self.y == 0])), .1)
                self.assertEqual(model.estimator.n_iter_, 200)


if __name__ == "__main__":
    unittest.main()
