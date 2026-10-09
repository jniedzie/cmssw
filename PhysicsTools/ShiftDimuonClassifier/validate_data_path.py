#!/usr/bin/env python3
"""Compare MC scoring with ROOT fixtures containing only reconstructed columns."""
import argparse
from collections import Counter
import json
from pathlib import Path
import pickle

import numpy as np
import uproot
import awkward as ak

from features import FEATURE_NAMES, RECO_BRANCHES, extract_reco
from inference import predict


def choose_sources(manifest, files_per_process):
    if manifest.get("sample_kind") != "simulation_only":
        raise ValueError("This validation is restricted to explicitly identified MC inputs")
    inputs = manifest.get("inputs", [])
    selected = []
    for process in ("jpsi", "dy", "qcd"):
        candidates = [row for row in inputs if row.get("process") == process]
        # Prefer the first published files with pairs so every process exercises scoring.
        with_pairs = [row for row in candidates if row.get("pairs", 0) > 0]
        ordered = with_pairs + [row for row in candidates if row not in with_pairs]
        if not ordered:
            raise ValueError("Missing MC process in manifest: " + process)
        selected.extend(ordered[:files_per_process])
    return selected


def read_models(directory):
    paths = sorted(directory.rglob("*.pkl"))
    if not paths:
        raise ValueError("No study model bundles found in " + str(directory))
    bundles = {}
    for path in paths:
        with path.open("rb") as source:
            bundle = pickle.load(source)
        if tuple(bundle["feature_names"]) != FEATURE_NAMES:
            raise ValueError("Model feature contract differs: " + str(path))
        columns = np.asarray(bundle["columns"], dtype=int)
        if columns.ndim != 1 or len(columns) == 0 or (columns < 0).any() or (columns >= len(FEATURE_NAMES)).any():
            raise ValueError("Invalid model input columns: " + str(path))
        bundles[str(path.relative_to(directory))] = bundle
    return bundles


def write_reco_fixture(path, original):
    """Write the ordinary TTree schema with shared Nano collection counters."""
    payload = {k: original[k] for k in ("run", "luminosityBlock", "event")}
    for collection in ("ShiftMuon", "ShiftDimuonVertex"):
        prefix = collection + "_"
        payload[collection] = ak.zip({k[len(prefix):]: original[k]
                                      for k in RECO_BRANCHES if k.startswith(prefix)})
    # Each jagged record creates one n<collection> counter and the original
    # <collection>_<field> leaf names, rather than one new counter per leaf.
    with uproot.recreate(path) as root:
        root["Events"] = payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--models", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--files-per-process", type=int, default=3)
    args = parser.parse_args()
    if not 1 <= args.files_per_process <= 3:
        parser.error("Bounded validation accepts 1--3 source files per process")
    manifest = json.loads((args.dataset / "inputs.json").read_text())
    sources = choose_sources(manifest, args.files_per_process)
    models = read_models(args.models)
    args.output.mkdir(parents=True, exist_ok=False)
    fixtures = args.output / "fixtures"
    fixtures.mkdir()
    report = {
        "validated": False, "sample_kind": "simulation_only", "collision_data_accessed": False,
        "scope": "MC-only equivalence with truth columns physically removed; not approval to score collision data",
        "dataset": str(args.dataset.resolve()), "models": str(args.models.resolve()),
        "reconstructed_fields": list(RECO_BRANCHES), "files": [],
        "models_checked": list(models),
    }
    counts = Counter()
    for index, row in enumerate(sources):
        source_path = Path(row["path"])
        with uproot.open(source_path, handler=uproot.source.file.MultithreadedFileSource, num_workers=1) as root:
            tree = root["Events"]
            removed = sorted(set(tree.keys()) - set(RECO_BRANCHES))
            original = tree.arrays(list(RECO_BRANCHES), library="ak", how=dict)
        expected = extract_reco(original)
        fixture = fixtures / f"{index:02d}_{row['process']}_{source_path.parent.name}.root"
        write_reco_fixture(fixture, original)
        with uproot.open(fixture, handler=uproot.source.file.MultithreadedFileSource, num_workers=1) as root:
            actual_fields = sorted(root["Events"].keys())
            if set(actual_fields) != set(RECO_BRANCHES):
                raise AssertionError("Stripped fixture does not have the exact reconstructed-field contract")
            if any(key.startswith("GenPart_") or "hitTruth" in key or "hitGen" in key or "sim" in key
                   for key in actual_fields):
                raise AssertionError("Truth field remains in stripped ROOT fixture")
            stripped = root["Events"].arrays(list(RECO_BRANCHES), library="ak", how=dict)
        actual = extract_reco(stripped)
        for key in ("X", "nuisance", "ids", "event_index", "pair_index", "muon_indices"):
            np.testing.assert_array_equal(actual[key], expected[key], err_msg="Truth removal changed " + key)
        model_checks = {}
        for name, bundle in models.items():
            if len(expected["X"]):
                before = predict(bundle, expected["X"])
                after = predict(bundle, actual["X"])
                np.testing.assert_array_equal(before, after, err_msg="Truth removal changed model " + name)
                if not np.isfinite(after).all():
                    raise AssertionError("Nonfinite prediction from " + name)
                np.testing.assert_array_equal(before >= bundle["threshold"], after >= bundle["threshold"])
                model_checks[name] = {"prediction_equal": True, "selection_equal": True,
                                      "max_absolute_score_difference": float(np.max(np.abs(after - before)))}
            else:
                model_checks[name] = {"prediction_equal": None, "reason": "No retained pairs in source file"}
        counts[row["process"]] += len(actual["X"])
        report["files"].append({"source": str(source_path), "process": row["process"],
                                "fixture": str(fixture.resolve()), "events": len(original["event"]),
                                "pairs": len(actual["X"]), "removed_fields": removed,
                                "fixture_fields": actual_fields, "model_checks": model_checks})
    if not sum(counts.values()):
        raise ValueError("No reconstructed pairs were exercised by the stripped-file test")
    report.update(validated=True, pairs_per_process=dict(counts), files_checked=len(sources))
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"validated": True, "output": str(args.output), "models_checked": list(models),
                      "files_checked": len(sources), "pairs_per_process": dict(counts)}), flush=True)


if __name__ == "__main__":
    main()
