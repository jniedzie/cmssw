"""Fixed-benchmark signal audits preserve collision-only scoring and denominators."""
import copy
import json
from pathlib import Path
import pickle
import sys
import tempfile
import unittest

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from features import FEATURE_NAMES, RECO_BRANCHES
from labels import TRUTH_BRANCHES
from signal_grid_audit import (SCHEMA, digest, event_summary, locked_bundle, parent_ancestor,
                               run_audit, signal_associations, validate_manifest)
from test_contract import reconstructed_event


class FixtureEstimator:
    def predict_proba(self, X):
        score = np.where(X[:, 0] > .1, .8, .2)
        return np.column_stack((1 - score, score))


def synthetic_nano(path):
    original = reconstructed_event()
    payload = {k: np.array([original[k][0], original[k][0] + 1])
               for k in ("run", "luminosityBlock", "event")}
    muons = {key[len("ShiftMuon_"):]: ak.Array([original[key][0].tolist(), original[key][0].tolist()])
             for key in RECO_BRANCHES if key.startswith("ShiftMuon_")}
    for key, values in (("hitGenPartIdx", [0, 1]), ("hitTruthPurity", [1., 1.]),
                        ("hitTruthMatchedLayers", [6, 6]), ("hitSimTrackId", [10, 11])):
        muons[key] = ak.Array([values, values])
    payload["ShiftMuon"] = ak.zip(muons)
    payload["ShiftDimuonVertex"] = ak.zip({
        key[len("ShiftDimuonVertex_"):]: ak.Array([original[key][0].tolist(), []])
        for key in RECO_BRANCHES if key.startswith("ShiftDimuonVertex_")})
    generator = dict(pdgId=[13, -13, 32], status=[1, 1, 2], genPartIdxMother=[2, 2, -1],
                     vx=[0., 0., 0.], vy=[0., 0., 0.], vz=[14700., 14700., 14800.])
    payload["GenPart"] = ak.zip({key: ak.Array([values, values]) for key, values in generator.items()})
    with uproot.recreate(path) as root:
        root["Events"] = payload


def inputs_manifest(path, **changes):
    row = dict(path=str(path), point="m15_prompt", mass_gev=15., epsilon=1.e-5,
               proper_ctau_mm=.001, nano_sha256=digest(path))
    row.update(changes)
    return dict(schema=SCHEMA, sample_kind="simulation_only", inputs=[row])


def make_model(results, threshold=.5):
    results.mkdir()
    (results / "models").mkdir()
    name = "fixture_fixed"
    bundle = dict(model=FixtureEstimator(), columns=np.arange(len(FEATURE_NAMES)),
                  feature_names=FEATURE_NAMES, requires_gen=False, threshold=threshold, deployment_ready=False)
    with (results / "models" / (name + ".pkl")).open("wb") as source:
        pickle.dump(bundle, source)
    (results / "selection_locked.json").write_text(json.dumps(dict(
        representative_models=[name], thresholds={name: threshold}, selection_uses_test=False)))
    return name


class SignalAuditDenominators(unittest.TestCase):
    def test_zero_pair_events_and_multiple_pairs_count_once(self):
        result = event_summary(4, np.array([0, 0, 1]), np.array([1, 1, 0]),
                               np.array([False, True, True]), np.array([True, True, False]),
                               np.array([True, True, True, False]))
        all_events = result["all_nano_event_fractions"]
        self.assertEqual(all_events["any_retained_pair"]["n"], 4)
        self.assertEqual(all_events["any_retained_pair"]["passed"], 2)
        self.assertEqual(all_events["any_selected_labelled_signal_pair"]["efficiency"], .25)
        recorded = result["recorded_signal_event_fractions"]
        self.assertEqual(recorded["any_selected_labelled_signal_pair"]["n"], 3)
        self.assertEqual(recorded["any_selected_labelled_signal_pair"]["efficiency"], 1 / 3)

    def test_empty_pairs_remain_undefined_conditionally_and_zero_per_event(self):
        result = event_summary(2, np.array([], dtype=int), np.array([], dtype=int),
                               np.array([], dtype=bool), np.array([], dtype=bool), np.array([True, True]))
        self.assertEqual(result["all_nano_event_fractions"]["any_retained_pair"]["efficiency"], 0.)

    def test_invalid_pair_event_index_fails(self):
        with self.assertRaisesRegex(ValueError, "out of range"):
            event_summary(1, np.array([1]), np.array([1]), np.array([True]), np.array([True]), np.array([True]))

    def test_parent_copy_ancestry_and_corrupt_cycles(self):
        pdg = np.array([13, -13, 32, 32])
        mothers = np.array([2, 3, 3, -1])
        self.assertEqual(parent_ancestor(pdg, mothers, 0), 3)
        self.assertEqual(parent_ancestor(pdg, mothers, 1), 3)
        self.assertEqual(parent_ancestor(pdg, np.array([2, 3, 3, 2]), 0), -1)
        self.assertEqual(parent_ancestor(pdg, np.array([2, 3, 3, 99]), 0), -1)
        self.assertEqual(parent_ancestor(pdg, mothers, -1), -1)

    def test_unknown_reco_match_never_identifies_signal_pair(self):
        truth = dict(GenPart_pdgId=[np.array([13, -13, 32])], GenPart_status=[np.array([1, 1, 2])],
                     GenPart_genPartIdxMother=[np.array([2, 2, -1])], ShiftMuon_hitGenPartIdx=[np.array([0, -1])])
        reco = dict(event_index=np.array([0]), muon_indices=np.array([[0, 1]]))
        pair, event = signal_associations(truth, reco)
        self.assertFalse(pair[0])
        self.assertTrue(event[0])
        truth["ShiftMuon_hitGenPartIdx"][0][1] = 1
        truth["GenPart_pdgId"][0][1] = 13
        pair, event = signal_associations(truth, reco)
        self.assertFalse(pair[0])
        self.assertFalse(event[0])


class SignalAuditProvenance(unittest.TestCase):
    def test_manifest_refuses_missing_hash_duplicate_files_and_data(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "nano.root"
            synthetic_nano(path)
            manifest = inputs_manifest(path)
            self.assertEqual(validate_manifest(manifest)["m15_prompt"], (15., 1.e-5, .001))
            for mutate in (lambda m: m["inputs"][0].pop("nano_sha256"),
                           lambda m: m["inputs"].append(m["inputs"][0].copy()),
                           lambda m: m.update(sample_kind="data"),
                           lambda m: m["inputs"][0].update(proper_ctau_mm=-1)):
                modified = copy.deepcopy(manifest)
                mutate(modified)
                with self.assertRaises(ValueError):
                    validate_manifest(modified)

    def test_changed_threshold_and_test_selected_model_are_refused(self):
        with tempfile.TemporaryDirectory() as d:
            results = Path(d) / "results"
            name = make_model(results)
            lock_path = results / "selection_locked.json"
            lock = json.loads(lock_path.read_text())
            lock["thresholds"][name] = .6
            lock_path.write_text(json.dumps(lock))
            with self.assertRaisesRegex(ValueError, "threshold differs"):
                locked_bundle(results, name)
            lock["selection_uses_test"] = True
            lock_path.write_text(json.dumps(lock))
            with self.assertRaisesRegex(ValueError, "validation-selected"):
                locked_bundle(results, name)

    def test_end_to_end_physical_truth_removal_and_signal_event_denominator(self):
        with tempfile.TemporaryDirectory() as d:
            base = Path(d)
            path = base / "nano.root"
            synthetic_nano(path)
            manifest_path = base / "inputs.json"
            manifest_path.write_text(json.dumps(inputs_manifest(path)))
            results = base / "results"
            name = make_model(results)
            report = run_audit(manifest_path, results, name, base / "audit")
            point = report["points"]["m15_prompt"]
            self.assertTrue(report["completed"])
            self.assertEqual(report["total_nano_events"], 2)
            self.assertEqual(point["common_signal_vertex"]["n"], 1)
            self.assertEqual(point["common_signal_vertex"]["efficiency"], 1.)
            self.assertEqual(point["events"]["recorded_signal_event_fractions"]["any_selected_labelled_signal_pair"]["efficiency"], .5)
            self.assertIsNone(report["attempted_or_accepted_full_gen_efficiency"])
            self.assertFalse(report["normalization_transfer_validated"])
            with uproot.open(base / "audit" / "truth_removed" / "0000.root") as root:
                self.assertEqual(set(root["Events"].keys()), set(RECO_BRANCHES))
                self.assertFalse(set(TRUTH_BRANCHES) & set(root["Events"].keys()))
            scores = dict(np.load(base / "audit" / "m15_prompt_scores.npz"))
            self.assertNotIn("y", scores)
            self.assertNotIn("signal_pair", scores)
            with self.assertRaises(FileExistsError):
                run_audit(manifest_path, results, name, base / "audit")


if __name__ == "__main__":
    unittest.main()
