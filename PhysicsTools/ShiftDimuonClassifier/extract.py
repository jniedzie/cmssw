#!/usr/bin/env python3
"""Freeze a bounded, completed MC inventory; extract separate reco/label tables."""
import argparse
import ast
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import uproot
from features import RECO_BRANCHES, extract_reco
from labels import TRUTH_BRANCHES, extract_labels


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def select_inputs(inventory, files_per_bin):
    bins = defaultdict(list)
    for line in Path(inventory).read_text().splitlines():
        row = ast.literal_eval(line)
        path = Path(row[0])
        if path.name != "nano.root" or "shift_detector_representative_20261007_v10" not in path.parts:
            raise ValueError("Only explicitly identified corrected V10 MC is allowed in this pilot")
        bins[path.parent.parent.name].append(str(path))
    result = []
    for stratum, paths in sorted(bins.items()):
        for path in sorted(set(paths))[:files_per_bin]:
            result.append(dict(path=path, stratum=stratum, process=stratum.split("_", 1)[0]))
    if not result or any(row["process"] not in ("jpsi", "dy", "qcd") for row in result):
        raise ValueError("Empty or unknown SM inventory")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--files-per-bin", type=int, default=40)
    args = parser.parse_args()
    if not 1 <= args.files_per_bin <= 100:
        parser.error("Bounded pilot accepts 1--100 files/bin")
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    inputs = select_inputs(args.inventory, args.files_per_bin)
    manifest = dict(sample_kind="simulation_only", campaign="shift_detector_representative_20261007_v10",
                    selection="lexicographically first published paths per stratum; not a random production estimate",
                    inventory=str(Path(args.inventory).resolve()), inventory_sha256=digest(args.inventory),
                    files_per_bin=args.files_per_bin, inputs=inputs, physics_valid=False)
    (output / "inputs.json").write_text(json.dumps(manifest, indent=2) + "\n")
    pieces = defaultdict(list)
    labels, reasons, process, strata, weights = [], [], [], [], []
    counts = Counter()
    started = time.monotonic()
    for index, row in enumerate(inputs):
        with uproot.open(row["path"], handler=uproot.source.file.MultithreadedFileSource, num_workers=1) as root:
            tree = root["Events"]
            missing = set(RECO_BRANCHES + TRUTH_BRANCHES + ("shiftSamplingGenWeight",)) - set(tree.keys())
            if missing:
                raise ValueError(f"Missing required columns in {row['path']}: {sorted(missing)}")
            arrays = tree.arrays(list(RECO_BRANCHES + TRUTH_BRANCHES + ("shiftSamplingGenWeight",)), library="ak", how=dict)
            reco = extract_reco(arrays)
            y, why = extract_labels(arrays, reco)
            n = len(y)
            for key, values in reco.items():
                if key != "feature_names":
                    pieces[key].append(values)
            labels.append(y)
            reasons.append(why)
            process.extend([row["process"]] * n)
            strata.extend([row["stratum"]] * n)
            weights.extend(float(arrays["shiftSamplingGenWeight"][i]) for i in reco["event_index"])
            row.update(events=len(arrays["event"]), pairs=n,
                       labels=dict(Counter(why.tolist())), source_bytes=Path(row["path"]).stat().st_size)
            counts["events"] += row["events"]
            counts["pairs"] += n
            counts.update(why.tolist())
        if (index + 1) % 50 == 0:
            print(json.dumps(dict(files=index + 1, **counts)), flush=True)
    combined = {key: np.concatenate(values) for key, values in pieces.items()}
    combined["feature_names"] = np.asarray(reco["feature_names"])
    # Process identity is used for diagnostics/splitting only, never as X input.
    combined.update(process=np.asarray(process), stratum=np.asarray(strata), sampling_weight=np.asarray(weights))
    # Duplicate physical events are not allowed to re-enter this bounded study.
    identities = [(p, *map(int, ids), int(j)) for p, ids, j in
                  zip(combined["process"], combined["ids"], combined["pair_index"])]
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate process/run/lumi/event/pair in selected files")
    np.savez_compressed(output / "reco.npz", **combined)
    np.savez_compressed(output / "labels.npz", y=np.concatenate(labels), reason=np.concatenate(reasons))
    manifest.update(counts=dict(counts), extraction_seconds=time.monotonic() - started,
                    reco_sha256=digest(output / "reco.npz"), labels_sha256=digest(output / "labels.npz"))
    (output / "inputs.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(dict(output=str(output), **counts)), flush=True)


if __name__ == "__main__":
    main()
