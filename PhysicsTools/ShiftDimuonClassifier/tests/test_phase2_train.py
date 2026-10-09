"""Phase-2 population weighting, acceptance ranking and holdout boundaries."""

from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import phase2_train as phase
from evaluation import event_split
from features import FEATURE_NAMES, QUALITY_NAMES


class FakeEstimator:
    """Cheap pickleable stand-in: test runner boundaries without fitting ML."""

    def __init__(self, seed=None):
        self.seed = seed

    def fit(self, X, y, **kwargs):
        self.fit_X = np.array(X, copy=True)
        self.fit_y = np.array(y, copy=True)
        self.feature_importances_ = np.full(X.shape[1], 1.0 / X.shape[1])
        return self

    def __getitem__(self, index):
        if index != -1:
            raise IndexError(index)
        return self

    def predict_proba(self, X):
        # The first quality feature stores row IDs solely in this test fixture.
        score = 0.25 + 0.5 * ((np.asarray(X)[:, 0] % 17) / 16.)
        return np.column_stack((1. - score, score))


def acceptance_fixture(failures):
    """Two processes, three supported cells each, identical aggregate scores."""
    per_cell = 60
    positive_nuisance = []
    positive_process = []
    scores = np.ones(2 * 3 * per_cell)
    fail_rows = []
    for p, process in enumerate(("jpsi", "dy")):
        for cell in range(3):
            positive_nuisance.extend([[cell + 1., 100. + cell, cell + 1.]] * per_cell)
            positive_process.extend([process] * per_cell)
            start = (3 * p + cell) * per_cell
            fail_rows.extend(range(start, start + failures[p][cell]))
    # Continuous failing scores avoid an artificial all-failing tie at the cut.
    scores[fail_rows] = np.linspace(0.01, 0.1, len(fail_rows))
    positive_nuisance = np.asarray(positive_nuisance)
    negative_nuisance = np.tile([[1., 100., 1.], [2., 101., 2.], [3., 102., 3.]], (30, 1))
    nuisance = np.vstack((positive_nuisance, negative_nuisance))
    labels = np.r_[np.ones(len(scores), dtype=int), np.zeros(90, dtype=int)]
    process = np.r_[positive_process, ["jpsi", "dy"] * 45]
    scores = np.r_[scores, np.zeros(45), np.ones(45)]
    return scores, labels, process, nuisance, np.ones(len(labels)), positive_nuisance, labels[:360]


class PhasePopulationWeights(unittest.TestCase):
    def test_single_qcd_negative_is_smoothed_and_class_balance_preserved(self):
        y = np.r_[np.zeros(201, dtype=int), np.ones(200, dtype=int)]
        process = np.asarray(["jpsi"] * 100 + ["dy"] * 100 + ["qcd"] +
                             ["jpsi"] * 100 + ["dy"] * 100)
        weights = phase.phase_weights(y, process)
        np.testing.assert_allclose([weights[y == label].sum() for label in (0, 1)],
                                   [len(y) / 2.] * 2, atol=1e-12)
        regular = weights[(y == 0) & (process == "jpsi")][0]
        rare = weights[(y == 0) & (process == "qcd")][0]
        # The floor of 50 limits this to 2x, versus 100x for naive balancing.
        self.assertAlmostEqual(rare / regular, 2.)
        self.assertLess(weights[process == "qcd"].sum() / weights[y == 0].sum(), .02)
        order = np.random.default_rng(13).permutation(len(y))
        np.testing.assert_array_equal(phase.phase_weights(y[order], process[order]), weights[order])


class PhaseFeatureGroups(unittest.TestCase):
    def test_groups_use_only_the_declared_reconstructed_columns(self):
        groups = phase.feature_groups()
        self.assertEqual(set(groups), {"bdt_quality", "bdt_geometry", "bdt_kinematics",
                                      "bdt_angle", "bdt_expanded"})
        for columns in groups.values():
            self.assertTrue(np.issubdtype(columns.dtype, np.integer))
            self.assertEqual(len(np.unique(columns)), len(columns))
            self.assertTrue(np.all((columns >= 0) & (columns < len(FEATURE_NAMES))))
            names = {FEATURE_NAMES[i] for i in columns}
            self.assertTrue(set(QUALITY_NAMES) <= names)
            for name in names:
                self.assertFalse(any(word in name.lower() for word in (
                    "genpart", "truth", "simtrack", "process", "sampling", "label",
                    "event", "nuisance", "mass")))
        np.testing.assert_array_equal(groups["bdt_expanded"], np.arange(len(FEATURE_NAMES)))
        self.assertEqual({FEATURE_NAMES[i] for i in groups["bdt_angle"]} - set(QUALITY_NAMES),
                         {"opening_angle"})
        geometry = {FEATURE_NAMES[i] for i in groups["bdt_geometry"]} - set(QUALITY_NAMES)
        kinematics = {FEATURE_NAMES[i] for i in groups["bdt_kinematics"]} - set(QUALITY_NAMES)
        self.assertTrue({"pair_vx", "pair_vy", "pair_vz", "mu_vertex_distance"} <= geometry)
        self.assertTrue({"pair_pt", "pair_pz", "mu_highPt_energy", "mu_lowPt_energy"} <= kinematics)
        self.assertFalse(geometry & kinematics)
        # Four-vectors still encode mass: this verifies column provenance only.


class ValidationAcceptanceRanking(unittest.TestCase):
    @staticmethod
    def summarize(failures):
        with patch.object(phase, "dependence", return_value={
            "common_vertex": {"joint_dcor": 0.}, "accidental": {"joint_dcor": 0.},
        }):
            # Isolate the acceptance terms from finite-sample dCor fluctuation.
            return phase.validation_summary(*acceptance_fixture(failures))

    def test_opposing_process_mass_acceptance_is_penalized_when_combined_flat(self):
        flat = self.summarize([[6, 6, 6], [6, 6, 6]])
        opposed = self.summarize([[12, 6, 0], [0, 6, 12]])
        self.assertEqual(opposed["accidental"]["efficiency"], flat["accidental"]["efficiency"])
        self.assertEqual(opposed["common_vertex"]["efficiency"], flat["common_vertex"]["efficiency"])
        for variable in ("mass", "vertex_z", "vertex_radius"):
            aggregate = [row["efficiency"] for row in opposed["flatness"][variable]["bins"]]
            self.assertLess(max(aggregate) - min(aggregate), .01)
            self.assertTrue(all(row["supported"] for row in
                                opposed["per_process_flatness"]["jpsi"][variable]["bins"]))
        self.assertGreater(opposed["maximum_marginal_efficiency_span"], .19)
        self.assertLess(opposed["selector_objective"], flat["selector_objective"] - .3)
        json.dumps(opposed, allow_nan=False)

    def test_different_process_means_are_penalized_even_with_flat_marginals(self):
        flat = self.summarize([[6, 6, 6], [6, 6, 6]])
        means = self.summarize([[12, 12, 12], [0, 0, 0]])
        self.assertEqual(means["common_vertex"]["efficiency"], flat["common_vertex"]["efficiency"])
        self.assertGreater(means["process_efficiency_span"], .19)
        self.assertLess(means["selector_objective"], flat["selector_objective"] - .3)


class PhaseRunnerBoundaries(unittest.TestCase):
    def test_invalid_mc_or_freshness_flags_stop_before_npz_reads_and_output(self):
        base = {"sample_kind": "simulation_only", "complete": True,
                "no_pilot_event_overlap": True, "no_pilot_file_overlap": True}
        invalid = [{"sample_kind": "collision_data"}, {"sample_kind": None}]
        for flag in ("complete", "no_pilot_event_overlap", "no_pilot_file_overlap"):
            invalid.extend([{flag: False}, {flag: None}, {flag: "false"}, {flag: 1}])
        with tempfile.TemporaryDirectory(prefix="shift_phase2_gate_") as directory:
            dataset = Path(directory) / "dataset"
            dataset.mkdir()
            for index, override in enumerate(invalid):
                manifest = dict(base, **override)
                (dataset / "inputs.json").write_text(json.dumps(manifest))
                output = Path(directory) / f"output_{index}"
                with self.subTest(override=override), patch.object(sys, "argv", [
                    "phase2_train.py", "--dataset", str(dataset), "--output", str(output),
                ]), patch.object(phase.np, "load", side_effect=AssertionError("must not load a table")):
                    with self.assertRaisesRegex(ValueError, "fresh MC dataset"):
                        phase.main()
                    self.assertFalse(output.exists())

    def test_all_fit_seeds_and_models_share_one_event_split_and_lock_before_test(self):
        rng = np.random.default_rng(281)
        n = 1200
        X = rng.normal(size=(n, len(FEATURE_NAMES)))
        X[:, 0] = np.arange(n)
        y = (np.arange(n) % 2).astype(np.int64)
        process = rng.choice(np.array(["dy", "jpsi", "qcd"]), n)
        ids = np.column_stack((np.ones(n), np.ones(n), np.arange(1, n + 1))).astype(np.uint64)
        nuisance = np.column_stack((rng.uniform(1, 10, n), rng.uniform(100, 200, n), rng.uniform(1, 10, n)))
        expected = event_split(ids, process, seed=71)
        expected_train = X[(expected == 0) & (y >= 0), 0]
        expected_val = X[(expected == 1) & (y >= 0), 0]
        bdt_models, other_calls, test_calls = [], [], []

        def fake_pipeline(*stages):
            model = FakeEstimator(stages[-1].random_state)
            bdt_models.append(model)
            return model

        def fake_uniform(*args, **kwargs):
            other_calls.append(("uniform", kwargs["seed"], args[0][:, 0].copy(), None))
            return FakeEstimator(kwargs["seed"]).fit(args[0], args[1]), {"rounds": []}

        def fake_neural(kind):
            def train(*args, **kwargs):
                other_calls.append((kind, kwargs["seed"], args[0][:, 0].copy(), args[4][:, 0].copy()))
                return FakeEstimator(kwargs["seed"]).fit(args[0], args[1]), {"epochs": []}
            return train

        with tempfile.TemporaryDirectory(prefix="shift_phase2_partition_") as directory:
            dataset, output = Path(directory) / "dataset", Path(directory) / "output"
            dataset.mkdir()
            np.savez_compressed(dataset / "reco.npz", X=X, ids=ids, nuisance=nuisance,
                                process=process, stratum=process, sampling_weight=np.ones(n),
                                feature_names=np.asarray(FEATURE_NAMES))
            np.savez_compressed(dataset / "labels.npz", y=y)
            manifest = {"sample_kind": "simulation_only", "complete": True,
                        "no_pilot_event_overlap": True, "no_pilot_file_overlap": True,
                        "synthetic_mechanics_fixture": True}
            for filename, key in (("reco.npz", "reco_sha256"), ("labels.npz", "labels_sha256")):
                manifest[key] = hashlib.sha256((dataset / filename).read_bytes()).hexdigest()
            (dataset / "inputs.json").write_text(json.dumps(manifest))

            def fake_evaluate(*args):
                self.assertTrue((output / "selection_locked.json").exists())
                lock = json.loads((output / "selection_locked.json").read_text())
                self.assertFalse(lock["selection_uses_test"])
                np.testing.assert_array_equal(args[3], ids[expected == 2])
                test_calls.append(1)
                return {"auc": .5, "common_vertex": {"efficiency": .9},
                        "background_rejection": .1,
                        "dependence": {"common_vertex": {"joint_dcor": .2}}}

            with patch.object(sys, "argv", ["phase2_train.py", "--dataset", str(dataset),
                                             "--output", str(output), "--epochs", "1"]), \
                    patch.object(phase, "event_split", wraps=event_split) as split_spy, \
                    patch.object(phase, "make_pipeline", side_effect=fake_pipeline), \
                    patch.object(phase, "train_uniform_bdt", side_effect=fake_uniform), \
                    patch.object(phase, "train_disco", side_effect=fake_neural("disco")), \
                    patch.object(phase, "train_adversarial", side_effect=fake_neural("adversarial")), \
                    patch.object(phase, "evaluate", side_effect=fake_evaluate), \
                    redirect_stdout(io.StringIO()):
                phase.main()
            self.assertEqual(split_spy.call_count, 1)
            self.assertEqual(split_spy.call_args.args[2], 71)
            self.assertEqual(len(bdt_models), 10)
            self.assertEqual(len(other_calls), 16)
            self.assertEqual(len(test_calls), 26)
            for model in bdt_models:
                np.testing.assert_array_equal(model.fit_X[:, 0], expected_train)
            for _, seed, train_rows, val_rows in other_calls:
                self.assertIn(seed, (71, 72))
                np.testing.assert_array_equal(train_rows, expected_train)
                if val_rows is not None:
                    np.testing.assert_array_equal(val_rows, expected_val)
            saved = np.load(output / "splits.npz", allow_pickle=False)
            np.testing.assert_array_equal(saved["split"], expected)
            result = json.loads((output / "metrics.json").read_text())
            self.assertTrue(result["selection_locked"])
            self.assertEqual(result["protocol"]["fit_seeds"], [71, 72])
            self.assertEqual(result["protocol"]["split_seed"], 71)
            self.assertEqual(len(result["models"]), 26)


if __name__ == "__main__":
    unittest.main()
