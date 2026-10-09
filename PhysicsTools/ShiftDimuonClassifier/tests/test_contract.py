"""Independent checks of the data feature and MC labeling contracts."""
import copy
from pathlib import Path
import sys
import unittest
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from features import FEATURE_NAMES, RECO_BRANCHES, extract_reco
from labels import pair_label
from evaluation import event_split
from validate_data_path import write_reco_fixture
import awkward as ak
import uproot


def reconstructed_event():
    arrays = {
        "run": np.array([1]), "luminosityBlock": np.array([2]),
        "event": np.array([101], dtype=np.uint64),
        "nShiftMuon": np.array([2]), "nShiftDimuonVertex": np.array([1]),
    }
    for branch in RECO_BRANCHES:
        if branch in arrays:
            continue
        arrays[branch] = [np.array([1.0, 2.0]) if branch.startswith("ShiftMuon_")
                          else np.array([0.2])]
    arrays.update({
        "ShiftMuon_pt": [np.array([3.0, 5.0])],
        "ShiftMuon_pz": [np.array([20.0, 25.0])],
        "ShiftMuon_phi": [np.array([0.1, -0.2])],
        "ShiftDimuonVertex_muonIdx1": [np.array([0])],
        "ShiftDimuonVertex_muonIdx2": [np.array([1])],
        "ShiftDimuonVertex_mass": [np.array([3.1])],
        "ShiftDimuonVertex_vz": [np.array([14700.0])],
    })
    return arrays


class PoisonedTruth(dict):
    def __getitem__(self, key):
        if key.startswith("GenPart_") or any(word in key for word in ("hitTruth", "hitGen", "simTrack")):
            raise AssertionError("Inference attempted to read truth: " + key)
        return super().__getitem__(key)


def truth_event():
    return {
        "ShiftMuon_hitGenPartIdx": np.array([0, 1]),
        "ShiftMuon_hitTruthPurity": np.array([1.0, 1.0]),
        "ShiftMuon_hitTruthMatchedLayers": np.array([6, 6]),
        "ShiftMuon_hitSimTrackId": np.array([20, 21]),
        "GenPart_pdgId": np.array([13, -13]),
        "GenPart_status": np.array([1, 1]),
        "GenPart_genPartIdxMother": np.array([4, 5]),
        "GenPart_vx": np.array([0.0, 0.0]),
        "GenPart_vy": np.array([0.0, 0.0]),
        "GenPart_vz": np.array([14700.0, 14700.0]),
    }


class DataFeatureContract(unittest.TestCase):
    def test_no_truth_reads_with_poisoned_extra_columns(self):
        arrays = reconstructed_event()
        expected = extract_reco(arrays)
        poisoned = PoisonedTruth(arrays)
        for key in ("GenPart_pdgId", "GenPart_vz", "ShiftMuon_hitGenPartIdx",
                    "ShiftMuon_hitTruthPurity", "ShiftMuon_simTrackP"):
            poisoned[key] = object()
        actual = extract_reco(poisoned)
        for key in expected:
            np.testing.assert_array_equal(actual[key], expected[key])
        self.assertFalse(any("GenPart" in key or "Truth" in key or "sim" in key
                             for key in RECO_BRANCHES))

    def test_pair_order_does_not_change_features(self):
        for equal_pt in (False, True):
            arrays = reconstructed_event()
            if equal_pt:
                arrays["ShiftMuon_pt"][0][:] = 4.0
            original = extract_reco(arrays)
            swapped = copy.deepcopy(arrays)
            swapped["ShiftDimuonVertex_muonIdx1"][0][0] = 1
            swapped["ShiftDimuonVertex_muonIdx2"][0][0] = 0
            result = extract_reco(swapped)
            for key in ("X", "ids", "nuisance"):
                np.testing.assert_array_equal(result[key], original[key])

    def test_mass_is_a_nuisance_not_a_direct_feature(self):
        arrays = reconstructed_event()
        original = extract_reco(arrays)
        arrays["ShiftDimuonVertex_mass"][0][0] = 15.0
        changed = extract_reco(arrays)
        np.testing.assert_array_equal(changed["X"], original["X"])
        self.assertNotEqual(changed["nuisance"][0, 0], original["nuisance"][0, 0])
        self.assertEqual(changed["X"].shape[1], len(FEATURE_NAMES))

    def test_invalid_pair_reference_fails(self):
        arrays = reconstructed_event()
        arrays["ShiftDimuonVertex_muonIdx2"][0][0] = 2
        with self.assertRaisesRegex(ValueError, "Invalid reconstructed pair"):
            extract_reco(arrays)

    def test_missing_reconstructed_column_fails(self):
        arrays = reconstructed_event()
        del arrays["ShiftMuon_pt"]
        with self.assertRaisesRegex(ValueError, "Missing reconstructed branches"):
            extract_reco(arrays)


class LabelContract(unittest.TestCase):
    def test_common_vertex_different_mothers_is_positive(self):
        self.assertEqual(pair_label(truth_event(), 0, 1), (1, "common_vertex"))

    def test_distant_vertices_same_mother_are_negative(self):
        event = truth_event()
        event["GenPart_genPartIdxMother"][:] = 4
        event["GenPart_vx"][1] = 2.0
        self.assertEqual(pair_label(event, 0, 1), (0, "different_vertices"))

    def test_duplicate_gen_muon_is_negative(self):
        event = truth_event()
        event["ShiftMuon_hitGenPartIdx"][1] = 0
        self.assertEqual(pair_label(event, 0, 1), (0, "same_muon_twice"))

    def test_duplicate_sim_track_is_negative(self):
        event = truth_event()
        event["ShiftMuon_hitSimTrackId"][1] = 20
        self.assertEqual(pair_label(event, 0, 1), (0, "same_muon_twice"))

    def test_unmatched_or_unreliable_truth_is_unknown(self):
        for branch, value in (("ShiftMuon_hitGenPartIdx", -1),
                              ("ShiftMuon_hitGenPartIdx", 99),
                              ("ShiftMuon_hitTruthPurity", 0.7),
                              ("ShiftMuon_hitTruthMatchedLayers", 2),
                              ("GenPart_pdgId", 211), ("GenPart_status", 2)):
            with self.subTest(branch=branch, value=value):
                event = truth_event()
                event[branch][1] = value
                self.assertEqual(pair_label(event, 0, 1), (-1, "unknown_match"))

    def test_position_tolerance_gap_is_unknown(self):
        event = truth_event()
        event["GenPart_vx"][1] = 0.05
        self.assertEqual(pair_label(event, 0, 1), (-1, "vertex_tolerance_gap"))

    def test_nonfinite_truth_position_is_unknown(self):
        event = truth_event()
        event["GenPart_vx"][1] = np.nan
        self.assertEqual(pair_label(event, 0, 1), (-1, "unknown_vertex"))

    def test_nonfinite_match_purity_is_unknown(self):
        event = truth_event()
        event["ShiftMuon_hitTruthPurity"][1] = np.nan
        self.assertEqual(pair_label(event, 0, 1), (-1, "unknown_match"))


class EventPartitionContract(unittest.TestCase):
    def test_all_pairs_from_same_event_stay_together(self):
        identities = np.repeat(np.array([[1, 1, x] for x in range(100)], dtype=np.uint64), 3, axis=0)
        processes = np.repeat(np.array(["jpsi"] * 50 + ["qcd"] * 50), 3)
        split = event_split(identities, processes)
        np.testing.assert_array_equal(split.reshape(-1, 3)[:, 0], split.reshape(-1, 3)[:, 1])
        np.testing.assert_array_equal(split.reshape(-1, 3)[:, 0], split.reshape(-1, 3)[:, 2])
        self.assertEqual(set(split.tolist()), {0, 1, 2})

    def test_event_partitions_are_independent_of_row_order(self):
        identities = np.array([[1, 1, x] for x in range(100)], dtype=np.uint64)
        processes = np.array(["dy"] * 100)
        order = np.random.default_rng(91).permutation(len(identities))
        original = event_split(identities, processes, seed=7)
        shuffled = event_split(identities[order], processes[order], seed=7)
        np.testing.assert_array_equal(shuffled, original[order])
        np.testing.assert_array_equal(event_split(identities, processes, seed=7), original)


class RootFixtureContract(unittest.TestCase):
    def test_exact_reco_ttree_fields_and_values(self):
        original = {k: ak.Array(v) for k, v in reconstructed_event().items()}
        before = extract_reco(original)
        with tempfile.TemporaryDirectory(prefix="shift_classifier_contract_") as directory:
            path = Path(directory) / "reco.root"
            write_reco_fixture(path, original)
            with uproot.open(path, handler=uproot.source.file.MultithreadedFileSource, num_workers=1) as root:
                self.assertEqual(root["Events"].classname, "TTree")
                self.assertEqual(set(root["Events"].keys()), set(RECO_BRANCHES))
                after = extract_reco(root["Events"].arrays(list(RECO_BRANCHES), library="ak", how=dict))
        for k in ("X", "ids", "nuisance"):
            np.testing.assert_array_equal(after[k], before[k])


if __name__ == "__main__":
    unittest.main()
