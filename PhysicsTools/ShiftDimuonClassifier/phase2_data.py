#!/usr/bin/env python3
"""One-reader, resumable extraction of a fresh completed V10 SM subset."""
import argparse
import ast
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from features import FEATURE_NAMES, RECO_BRANCHES, extract_reco
from labels import TRUTH_BRANCHES, extract_labels, label_diagnostics


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def freeze_inputs(inventory, excluded_dataset, sm_files, qcd_files):
    previous = json.loads((excluded_dataset / "inputs.json").read_text())
    if previous.get("sample_kind") != "simulation_only":
        raise ValueError("Exclusion reference must be the previous MC study")
    excluded = {row["path"] for row in previous["inputs"]}
    bins = defaultdict(set)
    for line in inventory.read_text().splitlines():
        path = Path(ast.literal_eval(line)[0])
        if path.name != "nano.root" or "shift_detector_representative_20261007_v10" not in path.parts:
            raise ValueError("Inventory must explicitly identify corrected V10 MC")
        if str(path) not in excluded:
            bins[path.parent.parent.name].add(str(path))
    result = []
    for stratum, paths in sorted(bins.items()):
        process = stratum.split("_", 1)[0]
        if process not in ("jpsi", "dy", "qcd"):
            raise ValueError("Unknown SM process")
        # Selection depends on immutable paths only, not labels or observables.
        ordered = sorted(paths, key=lambda p: hashlib.sha256(("phase2:" + p).encode()).hexdigest())
        limit = qcd_files if process == "qcd" else sm_files
        for path in ordered[:limit]:
            result.append(dict(path=path, process=process, stratum=stratum))
    return dict(sample_kind="simulation_only", campaign="shift_detector_representative_20261007_v10",
                phase="fresh-sample-2", physics_valid=False,
                selection="path-hash ordering; no observable/label selection; bounded completed QCD and SM files",
                inventory=str(inventory.resolve()), inventory_sha256=digest(inventory),
                excluded_dataset=str(excluded_dataset.resolve()),
                excluded_manifest_sha256=digest(excluded_dataset / "inputs.json"),
                sm_files_per_bin=sm_files, qcd_files_per_bin=qcd_files,
                source_sha256={p.name: digest(p) for p in (Path(__file__), Path(__file__).with_name("features.py"),
                                                          Path(__file__).with_name("labels.py"))},
                inputs=result)


def read_source(row):
    """Read only allowed leaf data; skip full collections in zero-pair events."""
    import ROOT
    root = ROOT.TFile(row["path"], "READ")
    if not root or root.IsZombie():
        raise ValueError("Unreadable completed Nano: " + row["path"])
    try:
        tree = root.Get("Events")
        if not tree:
            raise ValueError("Events tree absent")
        keys = {b.GetName() for b in tree.GetListOfBranches()}
        wanted = tuple(dict.fromkeys(RECO_BRANCHES + TRUTH_BRANCHES + ("shiftSamplingGenWeight",)))
        missing = set(wanted) - keys
        if missing:
            raise ValueError("Required branches absent: " + repr(sorted(missing)))
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("nShiftDimuonVertex", 1)
        n_events = int(tree.GetEntries())
        selected_entries = []
        for e in range(n_events):
            if tree.GetEntry(e) <= 0:
                raise ValueError("Unreadable event counter")
            if int(tree.nShiftDimuonVertex) > 0:
                selected_entries.append(e)
        arrays = {k: [] for k in wanted}
        for k in wanted:
            tree.SetBranchStatus(k, 1)
        scalar = {"run", "luminosityBlock", "event", "nShiftMuon", "nShiftDimuonVertex", "shiftSamplingGenWeight"}
        for e in selected_entries:
            if tree.GetEntry(e) <= 0:
                raise ValueError("Unreadable reconstructed event")
            for k in wanted:
                value = getattr(tree, k)
                arrays[k].append(value if k in scalar else list(value))
        reco = extract_reco(arrays)
        y, reasons = extract_labels(arrays, reco)
        diagnostics = label_diagnostics(arrays, reco)
        weights = np.asarray([arrays["shiftSamplingGenWeight"][e] for e in reco["event_index"]], dtype=float)
        reco["event_index"] = np.asarray([selected_entries[e] for e in reco["event_index"]], dtype=int)
        reco.pop("feature_names")
        return reco, dict(y=y, reason=reasons, **diagnostics), weights, n_events
    finally:
        root.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--exclude-dataset", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sm-files-per-bin", type=int, default=80)
    parser.add_argument("--qcd-files-per-bin", type=int, default=650)
    args = parser.parse_args()
    if not 1 <= args.sm_files_per_bin <= 150 or not 1 <= args.qcd_files_per_bin <= 650:
        parser.error("Bounds: 1--150 SM files/bin, 1--650 QCD files/bin")
    manifest_path = args.output / "inputs.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("complete"):
            raise ValueError("Completed dataset is immutable; choose a new output")
        if manifest["inventory_sha256"] != digest(args.inventory) or manifest["excluded_manifest_sha256"] != digest(args.exclude_dataset / "inputs.json"):
            raise ValueError("Resume references changed")
        if (manifest["sm_files_per_bin"], manifest["qcd_files_per_bin"]) != (args.sm_files_per_bin, args.qcd_files_per_bin):
            raise ValueError("Resume limits differ from frozen manifest")
        for p in (Path(__file__), Path(__file__).with_name("features.py"), Path(__file__).with_name("labels.py")):
            if manifest["source_sha256"][p.name] != digest(p):
                raise ValueError("Extractor source changed; use a new dataset rather than mixing implementations")
    else:
        args.output.mkdir(parents=True, exist_ok=False)
        manifest = freeze_inputs(args.inventory, args.exclude_dataset, args.sm_files_per_bin, args.qcd_files_per_bin)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    chunks = args.output / "chunks"
    chunks.mkdir(exist_ok=True)
    started = time.monotonic()
    pieces, labels, process, strata, sampling_weights = defaultdict(list), defaultdict(list), [], [], []
    counts = Counter()
    for index, row in enumerate(manifest["inputs"]):
        shard = chunks / f"g{index//250:03d}"
        shard.mkdir(exist_ok=True)
        checkpoint = shard / f"source{index:05d}.npz"
        receipt = checkpoint.with_suffix(".json")
        if row.get("pairs") == 0 and isinstance(row.get("events"), int):
            empty = {k: [] for k in RECO_BRANCHES + TRUTH_BRANCHES}
            reco = extract_reco(empty)
            y, reason = extract_labels(empty, reco)
            truth_labels = dict(y=y, reason=reason, **label_diagnostics(empty, reco))
            reco.pop("feature_names")
            weights, n_events = np.empty(0), row["events"]
        elif checkpoint.exists() and receipt.exists():
            info = json.loads(receipt.read_text())
            if info.get("path") != row["path"] or info.get("sha256") != digest(checkpoint):
                raise ValueError("Checkpoint identity or content changed")
            archive = dict(np.load(checkpoint, allow_pickle=False))
            reco = {k[len("reco_"):]: v for k, v in archive.items() if k.startswith("reco_")}
            truth_labels = {k[len("label_"):]: v for k, v in archive.items() if k.startswith("label_")}
            weights, n_events = archive["sampling_weights"], info["events"]
        else:
            if checkpoint.exists() or receipt.exists():
                raise ValueError("Incomplete checkpoint needs inspection; no overwrite")
            reco, truth_labels, weights, n_events = read_source(row)
            if len(truth_labels["y"]):
                np.savez_compressed(checkpoint, sampling_weights=weights,
                                    **{"reco_"+k: v for k,v in reco.items()},
                                    **{"label_"+k: v for k,v in truth_labels.items()})
                receipt.write_text(json.dumps(dict(path=row["path"], events=n_events, sha256=digest(checkpoint))) + "\n")
        n_pairs = len(truth_labels["y"])
        for k, v in reco.items():
            pieces[k].append(v)
        for k, v in truth_labels.items():
            labels[k].append(v)
        process.extend([row["process"]] * n_pairs)
        strata.extend([row["stratum"]] * n_pairs)
        sampling_weights.append(weights)
        row.update(events=n_events, pairs=n_pairs, labels=dict(Counter(truth_labels["reason"].tolist())))
        counts["events"] += n_events
        counts["pairs"] += n_pairs
        counts.update(truth_labels["reason"].tolist())
        if (index + 1) % 100 == 0:
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
            print(json.dumps(dict(files=index+1, **counts)), flush=True)
    combined = {k: np.concatenate(v) for k,v in pieces.items()}
    combined.update(feature_names=np.asarray(FEATURE_NAMES), process=np.asarray(process), stratum=np.asarray(strata),
                    sampling_weight=np.concatenate(sampling_weights))
    all_labels = {k: np.concatenate(v) for k,v in labels.items()}
    identities = [(p, *map(int, identity), int(pair)) for p, identity, pair in
                  zip(combined["process"], combined["ids"], combined["pair_index"])]
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate physical event pair")
    previous = dict(np.load(args.exclude_dataset / "reco.npz", allow_pickle=False))
    previous_groups = {(p, *map(int, i)) for p,i in zip(previous["process"], previous["ids"])}
    fresh_groups = {(p, *map(int, i)) for p,i in zip(combined["process"], combined["ids"])}
    if previous_groups & fresh_groups:
        raise ValueError("Fresh dataset has physical event overlap with pilot")
    np.savez_compressed(args.output / "reco.npz", **combined)
    np.savez_compressed(args.output / "labels.npz", **all_labels)
    manifest.update(complete=True, counts=dict(counts), extraction_seconds=time.monotonic()-started,
                    no_pilot_file_overlap=True, no_pilot_event_overlap=True,
                    reco_sha256=digest(args.output / "reco.npz"), labels_sha256=digest(args.output / "labels.npz"))
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(dict(complete=True, files=len(manifest["inputs"]), **counts)), flush=True)


if __name__ == "__main__":
    main()
