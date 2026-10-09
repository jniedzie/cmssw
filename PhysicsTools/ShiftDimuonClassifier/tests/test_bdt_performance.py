"""Scientific checks for ROC, independent events, support and held-out masking."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bdt_performance import (SCHEMA, efficiency_summary, evaluate, heldout_mask,
                             roc_summary, split_registry, validate_manifest, weighted_efficiency)


class IndependentPerformance(unittest.TestCase):
    def test_exact_physical_ids_and_process_masks(self):
        event = 2 ** 60 + 17
        registry = {("qcd", 1, 2, event): 0, ("qcd", 1, 2, event + 1): 1,
                    ("qcd", 1, 2, event + 2): 2}
        ids = np.asarray([[1, 2, event], [1, 2, event + 1], [1, 2, event + 2],
                          [1, 2, event + 3], [1, 2, event]], dtype=np.uint64)
        mask = heldout_mask(["qcd", "qcd", "qcd", "qcd", "dy"], ids, registry)
        np.testing.assert_array_equal(mask, [False, False, True, True, True])

    def test_same_event_cannot_occur_in_train_and_test(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "split.npz"
            np.savez(path, ids=np.array([[1, 2, 101], [1, 2, 101]]),
                     process=np.array(["qcd", "qcd"]), split=np.array([0, 2]))
            with self.assertRaisesRegex(ValueError, "multiple training partitions"):
                split_registry(path)

    def test_tied_and_separated_roc_have_known_auc(self):
        y, groups = np.array([0, 0, 1, 1]), [("qcd", i) for i in range(4)]
        tied = roc_summary(y, [.5, .5, .5, .5], groups, repeats=30)
        separated = roc_summary(y, [.1, .2, .8, .9], groups, repeats=30)
        reversed_order = roc_summary(y, [.8, .9, .1, .2], groups, repeats=0)
        self.assertEqual(tied["auc"], .5)
        self.assertEqual(separated["auc"], 1.)
        self.assertEqual(reversed_order["auc"], 0.)
        self.assertEqual(tied["thresholds"][0], None)
        json.dumps(tied, allow_nan=False)

    def test_auc_does_not_replace_fixed_cut_efficiency(self):
        y, groups = np.array([0, 0, 1, 1]), [("event", i) for i in range(4)]
        first, second = np.array([.1, .2, .8, .9]), np.array([.01, .02, .08, .09])
        self.assertEqual(roc_summary(y, first, groups, 0)["auc"], roc_summary(y, second, groups, 0)["auc"])
        self.assertEqual(efficiency_summary(first[y == 1] >= .5, groups[2:], 0)["efficiency"], 1.)
        self.assertEqual(efficiency_summary(second[y == 1] >= .5, groups[2:], 0)["efficiency"], 0.)

    def test_empty_bin_is_unmeasured_with_no_false_zero(self):
        eff, roc = efficiency_summary([], []), roc_summary([], [], [])
        self.assertIsNone(eff["efficiency"])
        self.assertEqual(eff["ci95_pair_wilson"], [None, None])
        self.assertIsNone(roc["auc"])
        self.assertFalse(roc["supported"])

    def test_bootstrap_samples_events_together(self):
        # Ten passing pairs in one event and one failing pair in another.
        # An event bootstrap produces exactly 0, 10/11 or 1, rather than an
        # artificially precise independent-pair interval around 10/11.
        result = efficiency_summary([True] * 10 + [False], [("event", 1)] * 10 + [("event", 2)], repeats=200)
        self.assertEqual(result["events"], 2)
        self.assertEqual(result["ci95_event_bootstrap"], [0., 1.])
        self.assertLess(result["ci95_pair_wilson"][1], 1.)
        lone = efficiency_summary([True] * 10, [("event", 1)] * 10)
        self.assertEqual(lone["bootstrap_replicates"], 0)
        self.assertFalse(lone["bootstrap_informative"])

    def test_truth_labels_are_not_invented_for_unknown_pairs(self):
        data = dict(ids=np.asarray([[1, 1, i] for i in range(5)], dtype=np.uint64),
                    process=np.array(["qcd", "jpsi", "dy", "signal", "signal"]),
                    point=np.array(["", "", "", "m15", "m15"]),
                    score=np.array([.1, .9, .8, .95, .99]), y=np.array([0, 1, 1, 1, -1]),
                    reason=np.array(["different_vertices", "common_vertex", "common_vertex", "common_vertex", "unknown_match"]),
                    signal=np.array([False, False, False, True, True]), nuisance=np.ones((5, 3)))
        points = {"m15": dict(mass_gev=15., proper_ctau_mm=10., epsilon=1.e-8),
                  "pending": dict(mass_gev=30., proper_ctau_mm=20., epsilon=1.e-8)}
        report = evaluate(data, [], points, {}, .5, repeats=0, minimum=2)
        self.assertEqual(report["signal_points"]["m15"]["common_signal_pair_efficiency"]["pairs"], 1)
        self.assertEqual(report["signal_points"]["m15"]["unknown_pair_retention"]["pairs"], 1)
        self.assertIsNone(report["signal_points"]["pending"]["common_signal_pair_efficiency"]["efficiency"])
        self.assertFalse(report["signal_points"]["m15"]["supported"])
        self.assertEqual(report["per_process"]["dy"]["genuine_retention"]["efficiency"], 1.)

    def test_simulation_and_point_metadata_fail_closed(self):
        manifest = dict(schema=SCHEMA, sample_kind="simulation_only", inputs=[
            dict(path="/tmp/signal.root", process="signal", point="m15", mass_gev=15., proper_ctau_mm=10., epsilon=1.e-8)],
            planned_signal_points=[dict(point="m15", mass_gev=15., proper_ctau_mm=20., epsilon=1.e-8)])
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            validate_manifest(manifest)
        manifest["sample_kind"] = "collision_data"
        with self.assertRaisesRegex(ValueError, "simulation-only"):
            validate_manifest(manifest)

    def test_weighted_mixture_changes_auc_with_exact_expected_value(self):
        # One easy and one hard negative; upweighting the hard process/bin
        # changes AUC from 1/2 to 1/10 without changing any predictor or score.
        y, scores = np.array([0, 0, 1]), np.array([.1, .9, .5])
        groups = [("event", i) for i in range(3)]
        empirical = roc_summary(y, scores, groups, repeats=0)
        weighted = roc_summary(y, scores, groups, repeats=0, weights=[1., 9., 2.])
        self.assertEqual(empirical["auc"], .5)
        self.assertAlmostEqual(weighted["auc"], .1)
        self.assertEqual(weighted["negative_weight_sum"], 10.)
        self.assertEqual(weighted["positive_weight_sum"], 2.)
        # Overall class normalizations cancel in a ROC.
        self.assertEqual(weighted["auc"], roc_summary(y, scores, groups, 0, weights=[100., 900., .02])["auc"])

    def test_weighted_efficiency_and_bootstrap_keep_entire_events(self):
        result = weighted_efficiency([True, False], [("event", 1), ("event", 2)], [9., 1.], repeats=100)
        self.assertEqual(result["efficiency"], .9)
        self.assertEqual(result["ci95_event_bootstrap"], [0., 1.])
        self.assertAlmostEqual(result["event_weight_effective_n"], 100 / 82)

    def test_signed_and_zero_weights_never_become_probability_roc(self):
        y, scores, groups = [0, 1], [.1, .9], [("event", 1), ("event", 2)]
        for weights in ([-1., 2.], [0., 2.], [1., float("nan")]):
            result = roc_summary(y, scores, groups, weights=weights)
            self.assertIsNone(result["auc"])
            self.assertFalse(result["available"])
            self.assertFalse(result["supported"])
            json.dumps(result, allow_nan=False)
        signed = weighted_efficiency([False, True], groups, [-1., 2.])
        self.assertIsNone(signed["efficiency"])
        self.assertEqual(signed["signed_mc_acceptance_ratio"], 2.)
        self.assertFalse(signed["available"])


if __name__ == "__main__":
    unittest.main()
