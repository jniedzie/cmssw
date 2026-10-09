#!/usr/bin/env python3
"""Evaluate a locked SM benchmark on an explicit signal Nano inventory.

No fitting, threshold choice or lifetime decorrelation occurs here. Scoring reads
only RECO_BRANCHES; truth is read afterwards for label/evaluation bookkeeping.
Nano denominators describe recorded events, never attempted/accepted full GEN.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import pickle

import numpy as np
import uproot

from acceptance_audit import unknown_bounds
from evaluation import fraction
from features import FEATURE_NAMES, RECO_BRANCHES, extract_reco
from inference import predict
from labels import TRUTH_BRANCHES, extract_labels
from validate_data_path import write_reco_fixture


SCHEMA = "shift-signal-classifier-inputs-v1"
POINT_FIELDS = ("mass_gev", "epsilon", "proper_ctau_mm")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as source:
        for piece in iter(lambda: source.read(1024 * 1024), b""):
            h.update(piece)
    return h.hexdigest()


def validate_manifest(manifest):
    if manifest.get("schema") != SCHEMA or manifest.get("sample_kind") != "simulation_only":
        raise ValueError("Only an explicit simulation signal Nano manifest is accepted")
    inputs = manifest.get("inputs", [])
    if not inputs:
        raise ValueError("No signal Nano inputs")
    paths, points = set(), {}
    for row in inputs:
        path = Path(row["path"])
        if not path.is_absolute() or path.suffix != ".root":
            raise ValueError("Signal inputs require explicit absolute ROOT paths")
        resolved = str(path.resolve())
        if resolved in paths:
            raise ValueError("Duplicate signal Nano file")
        paths.add(resolved)
        point = row.get("point")
        if not isinstance(point, str) or not point or Path(point).name != point:
            raise ValueError("Invalid signal point name")
        metadata = tuple(float(row[k]) for k in POINT_FIELDS)
        if not all(np.isfinite(metadata)) or metadata[0] <= 0 or metadata[1] <= 0 or metadata[2] < 0:
            raise ValueError("Invalid mass, epsilon or proper lifetime metadata")
        if point in points and points[point] != metadata:
            raise ValueError("Inconsistent signal point metadata")
        points[point] = metadata
        expected = row.get("nano_sha256")
        if not isinstance(expected, str) or len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
            raise ValueError("Every source requires a frozen Nano SHA-256 digest")
    return points


def locked_bundle(results, model_name):
    if Path(model_name).name != model_name:
        raise ValueError("Invalid model name")
    lock_path = results / "selection_locked.json"
    lock = json.loads(lock_path.read_text())
    if lock.get("selection_uses_test") is not False or model_name not in lock.get("representative_models", []):
        raise ValueError("Model must be a validation-selected representative with a locked threshold")
    model_path = results / "models" / (model_name + ".pkl")
    # Pickles are local trusted study outputs; never fetch arbitrary model files.
    with model_path.open("rb") as source:
        bundle = pickle.load(source)
    if tuple(bundle["feature_names"]) != FEATURE_NAMES or bundle.get("requires_gen") is not False:
        raise ValueError("Model must satisfy the explicit reconstructed-only contract")
    cut = float(lock["thresholds"][model_name])
    if not np.isfinite(cut) or float(bundle["threshold"]) != cut:
        raise ValueError("Saved model threshold differs from its lock receipt")
    return bundle, dict(name=model_name, threshold=cut, model_sha256=digest(model_path),
                        lock_sha256=digest(lock_path), results=str(results.resolve()),
                        deployment_ready=False)


def parent_ancestor(pdg, mothers, index, parent_pdg_id=32):
    """Audit-only ancestry, allowing muon/resonance copies and invalid links.

    A cycle is unknown even when an apparent parent was visited: corrupt
    ancestry cannot establish a signal association.
    """
    visited, found = set(), -1
    while 0 <= index < len(pdg):
        if index in visited:
            return -1
        visited.add(index)
        if abs(int(pdg[index])) == abs(parent_pdg_id):
            found = index
        index = int(mothers[index])
    return found if index == -1 else -1


def signal_associations(truth, reco, parent_pdg_id=32):
    """Label-only association to the same recorded signal parent.

    This never supplies predictor columns, adversary targets or row selection.
    Event presence means a recorded final-state opposite-sign dimuon decay,
    without a fiducial/gen-selection denominator or normalization claim.
    """
    event_present = np.zeros(len(truth["GenPart_pdgId"]), dtype=bool)
    ancestors = []
    for e in range(len(event_present)):
        pdg = np.asarray(truth["GenPart_pdgId"][e])
        status = np.asarray(truth["GenPart_status"][e])
        mothers = np.asarray(truth["GenPart_genPartIdxMother"][e])
        if len(pdg) != len(mothers) or len(pdg) != len(status):
            raise ValueError("Generator collection lengths differ")
        parent = np.asarray([parent_ancestor(pdg, mothers, i, parent_pdg_id) for i in range(len(pdg))])
        ancestors.append(parent)
        charges = defaultdict(set)
        for i in np.flatnonzero((np.abs(pdg) == 13) & (status == 1)):
            if parent[i] >= 0:
                charges[int(parent[i])].add(int(pdg[i]))
        event_present[e] = any({13, -13}.issubset(signs) for signs in charges.values())
    pair_signal = np.zeros(len(reco["event_index"]), dtype=bool)
    for row, (e, muons) in enumerate(zip(reco["event_index"], reco["muon_indices"])):
        g = [int(truth["ShiftMuon_hitGenPartIdx"][e][i]) for i in muons]
        if all(0 <= i < len(ancestors[e]) for i in g):
            p = [int(ancestors[e][i]) for i in g]
            charges = {int(truth["GenPart_pdgId"][e][i]) for i in g}
            pair_signal[row] = p[0] >= 0 and p[0] == p[1] and charges == {13, -13}
    return pair_signal, event_present


def event_summary(n_events, event_index, y, passed, signal_pair, signal_event):
    """Count each event once and keep zero-candidate events in denominators."""
    event_index, y, passed, signal_pair = map(np.asarray, (event_index, y, passed, signal_pair))
    if any(a.shape != y.shape for a in (event_index, passed, signal_pair)) or signal_event.shape != (n_events,):
        raise ValueError("Pair/event alignment differs")
    if ((event_index < 0) | (event_index >= n_events)).any():
        raise ValueError("Pair event index out of range")
    masks = {}
    for name, pairs in (("any_retained_pair", np.ones(len(y), dtype=bool)),
                        ("any_selected_pair", passed), ("any_labelled_authentic_pair", y == 1),
                        ("any_selected_labelled_authentic_pair", (y == 1) & passed),
                        ("any_labelled_signal_pair", (y == 1) & signal_pair),
                        ("any_selected_labelled_signal_pair", (y == 1) & signal_pair & passed)):
        present = np.zeros(n_events, dtype=bool)
        present[event_index[pairs].astype(int)] = True
        masks[name] = present
    all_events = np.ones(n_events, dtype=bool)
    return dict(nano_events=n_events, recorded_signal_mumu_events=int(signal_event.sum()),
                all_nano_event_fractions={name: fraction(all_events, mask) for name, mask in masks.items()},
                recorded_signal_event_fractions={name: fraction(signal_event, mask) for name, mask in masks.items()},
                interpretation="Recorded Nano-event fractions include zero-pair events; not full-GEN/trigger acceptance")


def summarize_point(events, reco, y, reasons, scores, pair_signal, signal_event, threshold):
    passed = scores >= threshold
    return dict(metadata={key: events[0][key] for key in POINT_FIELDS}, source_files=len(events),
                retained_pairs=len(y), label_reasons=dict(Counter(reasons.tolist())),
                common_vertex=fraction(y == 1, passed), accidental=fraction(y == 0, passed),
                unknown_retention=fraction(y < 0, passed),
                common_signal_vertex=fraction((y == 1) & pair_signal, passed),
                unknown_bounds=unknown_bounds(y, passed),
                events=event_summary(len(signal_event), reco["event_index"], y, passed, pair_signal, signal_event),
                normalization_transfer_validated=False, lifetime_independence_validated=False)


def run_audit(manifest_path, results, model_name, output, parent_pdg_id=32):
    manifest = json.loads(manifest_path.read_text())
    validate_manifest(manifest)
    bundle, model_receipt = locked_bundle(results, model_name)
    output.mkdir(parents=True, exist_ok=False)
    (output / "benchmark_locked.json").write_text(json.dumps(model_receipt, indent=2) + "\n")
    fixtures = output / "truth_removed"
    fixtures.mkdir()
    points = defaultdict(lambda: dict(rows=[], reco=[], labels=[], reasons=[], scores=[], signal_pair=[], signal_event=[], events=0))
    seen_events, sources = set(), []
    report = dict(schema="shift-signal-classifier-audit-v1", completed=False, sample_kind="simulation_only",
                  input_manifest=str(manifest_path.resolve()), input_manifest_sha256=digest(manifest_path), model=model_receipt,
                  scope="Fixed SM benchmark evaluation only; reconstructed-only scoring and separate MC labels",
                  collision_data_accessed=False, deployment_ready=False, normalization_transfer_validated=False,
                  attempted_or_accepted_full_gen_efficiency=None, parent_pdg_id=parent_pdg_id,
                  pair_ci_note="Wilson intervals approximate pair independence; event fractions count each event once",
                  lifetime_note="proper_ctau_mm is point metadata; no lifetime or truth value enters scoring/decorrelation",
                  sources=sources, points={})
    for index, row in enumerate(manifest["inputs"]):
        path = Path(row["path"])
        if digest(path) != row["nano_sha256"]:
            raise ValueError("Signal Nano differs from frozen source digest: " + str(path))
        with uproot.open(path, handler=uproot.source.file.MultithreadedFileSource, num_workers=1) as root:
            tree = root["Events"]
            missing = set(RECO_BRANCHES + TRUTH_BRANCHES) - set(tree.keys())
            if missing:
                raise ValueError("Signal Nano schema lacks: " + ", ".join(sorted(missing)))
            arrays = tree.arrays(list(RECO_BRANCHES), library="ak", how=dict)
            reco = extract_reco(arrays)
            score = predict(bundle, reco["X"])
            if not np.isfinite(score).all():
                raise ValueError("Nonfinite signal scores")
            # This call is deliberately after inference, and kept separate.
            truth = tree.arrays(list(TRUTH_BRANCHES), library="ak", how=dict)
            y, reasons = extract_labels(truth, reco)
            signal_pair, signal_event = signal_associations(truth, reco, parent_pdg_id)
            if len(signal_event) != len(arrays["event"]):
                raise ValueError("Truth-label event count differs from reconstructed event count")
        fixture = fixtures / f"{index:04d}.root"
        write_reco_fixture(fixture, arrays)
        with uproot.open(fixture, handler=uproot.source.file.MultithreadedFileSource, num_workers=1) as root:
            if set(root["Events"].keys()) != set(RECO_BRANCHES):
                raise AssertionError("Truth-stripped signal fixture schema differs")
            stripped = extract_reco(root["Events"].arrays(list(RECO_BRANCHES), library="ak", how=dict))
        for key in ("X", "nuisance", "ids", "event_index", "pair_index", "muon_indices"):
            np.testing.assert_array_equal(reco[key], stripped[key])
        np.testing.assert_array_equal(score, predict(bundle, stripped["X"]))
        np.testing.assert_array_equal(score >= bundle["threshold"], predict(bundle, stripped["X"]) >= bundle["threshold"])
        n_events = len(arrays["event"])
        for identity in zip(*(arrays[k].to_list() for k in ("run", "luminosityBlock", "event"))):
            identity = (row["point"], *map(int, identity))
            if identity in seen_events:
                raise ValueError("Duplicate signal point/run/lumi/event")
            seen_events.add(identity)
        point = points[row["point"]]
        reco["event_index"] = reco["event_index"] + point["events"]
        point["events"] += n_events
        point["rows"].append(row)
        point["reco"].append(reco)
        for key, values in (("labels", y), ("reasons", reasons), ("scores", score),
                            ("signal_pair", signal_pair), ("signal_event", signal_event)):
            point[key].append(values)
        sources.append(dict(path=str(path), point=row["point"], nano_sha256=row["nano_sha256"],
                            nano_events=n_events, retained_pairs=len(y), fixture=str(fixture),
                            physical_truth_removal_equal=True, score_exercised=bool(len(y))))
    for name, point in points.items():
        reco = {key: np.concatenate([r[key] for r in point["reco"]])
                for key in ("X", "nuisance", "ids", "event_index", "pair_index", "muon_indices")}
        y, reasons, scores, pair_signal, signal_event = (np.concatenate(point[key]) for key in
                                                       ("labels", "reasons", "scores", "signal_pair", "signal_event"))
        report["points"][name] = summarize_point(point["rows"], reco, y, reasons, scores,
                                                pair_signal, signal_event, bundle["threshold"])
        np.savez_compressed(output / (name + "_scores.npz"), score=scores, passed=scores >= bundle["threshold"],
                            ids=reco["ids"], pair_index=reco["pair_index"], event_index=reco["event_index"],
                            nuisance=reco["nuisance"])
        # Labels and signal ancestry are saved separately from score/predictor data.
        np.savez_compressed(output / (name + "_labels.npz"), y=y, reason=reasons,
                            signal_pair=pair_signal, recorded_signal_event=signal_event)
    report.update(completed=True, total_nano_events=sum(p["events"] for p in points.values()),
                  truth_removed_sources=len(sources), scoring_sources=sum(s["score_exercised"] for s in sources),
                  no_model_or_threshold_changes=True)
    (output / "signal_grid_audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--model", default="uniform_bdt_3p0_s71")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--parent-pdg-id", type=int, default=32)
    args = parser.parse_args()
    if args.parent_pdg_id == 0:
        parser.error("parent PDG ID must be nonzero")
    report = run_audit(args.manifest, args.results, args.model, args.output, args.parent_pdg_id)
    print(json.dumps(dict(output=str(args.output), points=list(report["points"]),
                          nano_events=report["total_nano_events"], scoring_sources=report["scoring_sources"],
                          normalization_transfer_validated=False)), flush=True)


if __name__ == "__main__":
    main()
