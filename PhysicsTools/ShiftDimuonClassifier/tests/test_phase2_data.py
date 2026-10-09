"""Fresh SM dataset and isolated label diagnostics, using synthetic inputs only."""

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import phase2_data
from features import FEATURE_NAMES
from label_audit import reconstruct_labels
from labels import label_diagnostics, pair_label
from plot_results import _name, _selected_models


def source_path(directory, stratum, job):
    # Synthetic paths identify the campaign convention but are never opened.
    return directory / "shift_detector_representative_20261007_v10" / stratum / f"job{job:07d}" / "nano.root"


def write_inventory(path, sources):
    path.write_text("\n".join(repr((str(source),)) for source in sources) + "\n")


def previous_dataset(directory, sources, identity=(1, 2, 3), process="jpsi"):
    directory.mkdir()
    (directory / "inputs.json").write_text(json.dumps(dict(
        sample_kind="simulation_only", inputs=[dict(path=str(path)) for path in sources])) + "\n")
    np.savez_compressed(directory / "reco.npz", ids=np.array([identity], dtype=np.uint64),
                        process=np.array([process]))


def truth_fixture():
    return {
        "ShiftMuon_hitGenPartIdx": np.array([0, 1]),
        "ShiftMuon_hitTruthPurity": np.array([1., 1.]),
        "ShiftMuon_hitTruthMatchedLayers": np.array([6, 6]),
        "ShiftMuon_hitSimTrackId": np.array([20, 21]),
        "GenPart_pdgId": np.array([13, -13, 443, 23]),
        "GenPart_status": np.array([1, 1, 2, 2]),
        "GenPart_genPartIdxMother": np.array([2, 3, -1, -1]),
        "GenPart_vx": np.array([0., 0., 0., 0.]),
        "GenPart_vy": np.array([0., 0., 0., 0.]),
        "GenPart_vz": np.array([14700., 14700., 14700., 14700.]),
    }


def diagnostics(event):
    return label_diagnostics({key: [value] for key, value in event.items()},
                             dict(event_index=np.array([0]), muon_indices=np.array([[0, 1]])))


class FreshInputContract(unittest.TestCase):
    def test_freeze_excludes_previous_sources_and_is_inventory_order_independent(self):
        with tempfile.TemporaryDirectory(prefix="shift_phase2_freeze_") as temporary:
            root = Path(temporary)
            sources = [source_path(root, stratum, job) for stratum in
                       ("jpsi_0to1", "dy_0to1", "qcd_0to1") for job in range(5)]
            excluded = [sources[0], sources[-1]]
            previous_dataset(root / "previous", excluded)
            write_inventory(root / "inventory_a.txt", sources + [sources[2]])
            write_inventory(root / "inventory_b.txt", list(reversed(sources)))
            first = phase2_data.freeze_inputs(root / "inventory_a.txt", root / "previous", 2, 1)
            second = phase2_data.freeze_inputs(root / "inventory_b.txt", root / "previous", 2, 1)
            self.assertEqual(first["inputs"], second["inputs"])
            selected = [row["path"] for row in first["inputs"]]
            self.assertEqual(len(selected), 5)
            self.assertEqual(len(set(selected)), 5)
            self.assertFalse(set(selected).intersection(map(str, excluded)))
            self.assertEqual(sum(row["process"] == "qcd" for row in first["inputs"]), 1)
            self.assertEqual(first["inventory_sha256"], hashlib.sha256((root / "inventory_a.txt").read_bytes()).hexdigest())

    def test_collision_reference_or_wrong_campaign_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix="shift_phase2_invalid_") as temporary:
            root = Path(temporary)
            previous_dataset(root / "previous", [])
            write_inventory(root / "inventory.txt", [root / "other_campaign" / "jpsi_0to1" / "job0000000" / "nano.root"])
            with self.assertRaisesRegex(ValueError, "corrected V10 MC"):
                phase2_data.freeze_inputs(root / "inventory.txt", root / "previous", 1, 1)
            write_inventory(root / "inventory.txt", [source_path(root, "jpsi_0to1", 0)])
            (root / "previous" / "inputs.json").write_text(json.dumps(dict(sample_kind="collision_data", inputs=[])))
            with self.assertRaisesRegex(ValueError, "previous MC study"):
                phase2_data.freeze_inputs(root / "inventory.txt", root / "previous", 1, 1)

    def test_different_source_and_pair_cannot_reuse_a_previous_physical_event(self):
        with tempfile.TemporaryDirectory(prefix="shift_phase2_event_") as temporary:
            root = Path(temporary)
            old_source = source_path(root, "jpsi_0to1", 0)
            new_source = source_path(root, "jpsi_0to1", 1)
            previous_dataset(root / "previous", [old_source], identity=(1, 2, 3))
            write_inventory(root / "inventory.txt", [old_source, new_source])
            reco = dict(X=np.zeros((1, len(FEATURE_NAMES))), ids=np.array([[1, 2, 3]], dtype=np.uint64),
                        event_index=np.array([0]), pair_index=np.array([7]),
                        muon_indices=np.array([[0, 1]]), nuisance=np.array([[3.1, 14700., 0.]]))
            labels = dict(y=np.array([1], dtype=np.int8), reason=np.array(["common_vertex"]))
            arguments = ["phase2_data.py", "--inventory", str(root / "inventory.txt"),
                         "--exclude-dataset", str(root / "previous"), "--output", str(root / "fresh"),
                         "--sm-files-per-bin", "1", "--qcd-files-per-bin", "1"]
            with patch.object(sys, "argv", arguments), patch.object(phase2_data, "read_source",
                    return_value=(reco, labels, np.ones(1), 1)) as reader:
                with self.assertRaisesRegex(ValueError, "physical event overlap"):
                    phase2_data.main()
                self.assertEqual(reader.call_count, 1)
                self.assertEqual(reader.call_args.args[0]["path"], str(new_source))
            self.assertFalse((root / "fresh" / "reco.npz").exists())
            self.assertFalse((root / "fresh" / "labels.npz").exists())
            self.assertFalse(json.loads((root / "fresh" / "inputs.json").read_text()).get("complete", False))


class DetailedLabelContract(unittest.TestCase):
    def test_nominal_diagnostics_reproduce_labels_and_unknown_reasons(self):
        variations = [None, ("GenPart_vx", 1, 2.), ("GenPart_vx", 1, 0.05),
                      ("GenPart_vx", 1, np.nan), ("ShiftMuon_hitGenPartIdx", 1, 0),
                      ("ShiftMuon_hitSimTrackId", 1, 20), ("ShiftMuon_hitGenPartIdx", 1, -1),
                      ("ShiftMuon_hitGenPartIdx", 1, 99), ("ShiftMuon_hitTruthPurity", 1, np.nan),
                      ("ShiftMuon_hitTruthPurity", 1, 0.7), ("ShiftMuon_hitTruthMatchedLayers", 1, 2),
                      ("GenPart_pdgId", 1, 211), ("GenPart_status", 1, 2)]
        for variation in variations:
            with self.subTest(variation=variation):
                event = truth_fixture()
                if variation:
                    event[variation[0]][variation[1]] = variation[2]
                original = copy.deepcopy(event)
                varied, reason = reconstruct_labels(diagnostics(event))
                self.assertEqual((int(varied[0]), str(reason[0])), pair_label(event, 0, 1))
                for key in original:
                    np.testing.assert_array_equal(event[key], original[key])

    def test_missing_mothers_remain_unavailable_instead_of_distinct(self):
        for mothers, expected in (([2, 3], 0), ([2, 2], 1), ([-1, 3], -1), ([2, 99], -1)):
            with self.subTest(mothers=mothers):
                event = truth_fixture()
                event["GenPart_genPartIdxMother"][:2] = mothers
                detailed = diagnostics(event)
                self.assertEqual(int(detailed["same_immediate_mother"][0]), expected)
                self.assertEqual(pair_label(event, 0, 1), (1, "common_vertex"))

    def test_looser_purity_changes_labels_only_when_the_gen_muon_is_available(self):
        event = truth_fixture()
        event["ShiftMuon_hitTruthPurity"][1] = 0.7
        detailed = diagnostics(event)
        self.assertEqual(detailed["match_state"][0, 1], "low_purity")
        self.assertEqual(int(reconstruct_labels(detailed)[0][0]), -1)
        self.assertEqual(int(reconstruct_labels(detailed, purity=0.5)[0][0]), 1)
        event["ShiftMuon_hitGenPartIdx"][1] = -1
        self.assertEqual(int(reconstruct_labels(diagnostics(event), purity=0.5)[0][0]), -1)
        event = truth_fixture()
        event["ShiftMuon_hitTruthPurity"][1] = np.nan
        self.assertEqual(int(reconstruct_labels(diagnostics(event), purity=0.)[0][0]), -1)

    def test_layer_and_vertex_sensitivity_preserves_gap_and_duplicate_conventions(self):
        event = truth_fixture()
        event["ShiftMuon_hitTruthMatchedLayers"][1] = 2
        detailed = diagnostics(event)
        self.assertEqual(detailed["match_state"][0, 1], "too_few_layers")
        self.assertEqual(int(reconstruct_labels(detailed)[0][0]), -1)
        self.assertEqual(int(reconstruct_labels(detailed, layers=2)[0][0]), 1)
        event = truth_fixture()
        event["GenPart_vx"][1] = 0.005
        detailed = diagnostics(event)
        self.assertEqual(int(reconstruct_labels(detailed)[0][0]), 1)
        varied, reason = reconstruct_labels(detailed, common_cm=0.001)
        self.assertEqual((int(varied[0]), str(reason[0])), (-1, "vertex_tolerance_gap"))
        event["ShiftMuon_hitSimTrackId"][1] = 20
        varied, reason = reconstruct_labels(diagnostics(event), common_cm=0.001)
        self.assertEqual((int(varied[0]), str(reason[0])), (0, "same_muon_twice"))


class PlotSelectionContract(unittest.TestCase):
    def test_legacy_complete_comparison_and_representative_override(self):
        metrics = dict(models={"bdt_quality": {}, "bdt_expanded": {}, "nn_disco_0p5_s71": {}})
        self.assertEqual(_selected_models(metrics), list(metrics["models"]))
        metrics["representative_models"] = ["bdt_quality", "nn_disco_0p5_s71"]
        self.assertEqual(_selected_models(metrics), metrics["representative_models"])
        self.assertEqual(_selected_models(metrics, ["bdt_expanded"]), ["bdt_expanded"])
        for invalid in ([], ["missing"], ["bdt_quality", "bdt_quality"], ["../bdt_quality"], "bdt_quality"):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                _selected_models(metrics, invalid)

    def test_names_preserve_legacy_labels_and_expose_scan_seed(self):
        self.assertEqual(_name("bdt_expanded"), "BDT expanded")
        self.assertEqual(_name("nn_lambda_0p5"), "NN lambda=0.5")
        self.assertEqual(_name("bdt_expanded_s71"), "BDT expanded (seed 71)")
        self.assertEqual(_name("s71_bdt_quality"), "BDT quality (seed 71)")
        self.assertEqual(_name("nn_disco_lambda_0p5_s71"), "NN DisCo lambda=0.5 (seed 71)")
        self.assertEqual(_name("nn_disco_0p5_s71"), "NN DisCo lambda=0.5 (seed 71)")
        self.assertEqual(_name("uniform_bdt_3p0_s71"), "Uniform BDT strength=3.0 (seed 71)")
        self.assertEqual(_name("nn_adversarial_0p5_s71"), "Adversarial NN lambda=0.5 (seed 71)")


if __name__ == "__main__":
    unittest.main()
