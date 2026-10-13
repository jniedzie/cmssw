#!/usr/bin/env python3
"""ROC and conditional BDT acceptance from explicitly identified simulation Nano.

Stored scores are checked against their portable-model metadata. Truth supplies
evaluation labels and signal-parent association only. Nothing is fitted here,
and the validation-locked threshold is never changed. N_eff is an MC-to-data
normalization; this program measures classifier performance, not N_eff.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from features import PREDICTOR_BRANCHES, extract_reco
from labels import TRUTH_BRANCHES, extract_labels
from portable_inference import file_digest, load_model, predict as predict_score


SCHEMA = "shift-bdt-performance-inputs-v1"
SCORE_BRANCH = "ShiftDimuonVertex_bdtScore"
POINT_FIELDS = ("mass_gev", "proper_ctau_mm", "epsilon")


def identity(process, ids):
    """Exact integer physical identities, unaffected by source relocation/merging."""
    return (str(process), *(int(i) for i in ids))


def split_registry(path):
    with np.load(path, allow_pickle=False) as data:
        ids, process, split = data["ids"], data["process"], data["split"]
    if ids.shape != (len(split), 3) or process.shape != split.shape or not np.isin(split, (0, 1, 2)).all():
        raise ValueError("Malformed saved training event partitions")
    result = {}
    for p, event, part in zip(process, ids, split):
        key = identity(p, event)
        if key in result and result[key] != int(part):
            raise ValueError("One physical event occurs in multiple training partitions")
        result[key] = int(part)
    return result


def heldout_mask(process, ids, registry):
    """Keep saved test and genuinely unseen events; exclude ALL fit/selection events."""
    return np.asarray([registry.get(identity(p, event), 2) == 2 for p, event in zip(process, ids)], dtype=bool)


def training_reference(training_results, training_dataset):
    """Check saved split alignment against the original digest-locked dataset."""
    metrics = json.loads((training_results / "metrics.json").read_text())
    manifest_path = training_dataset / "inputs.json"
    if file_digest(manifest_path) != metrics.get("dataset_manifest_sha256"):
        raise ValueError("Original training dataset manifest differs from locked results")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("sample_kind") != "simulation_only" or manifest.get("complete") is not True:
        raise ValueError("Original training dataset is not completed simulation")
    for filename, key in (("reco.npz", "reco_sha256"), ("labels.npz", "labels_sha256")):
        if file_digest(training_dataset / filename) != manifest.get(key):
            raise ValueError("Original training dataset content differs: " + filename)
    with np.load(training_dataset / "reco.npz", allow_pickle=False) as source:
        ids, process, nuisance = source["ids"], source["process"], source["nuisance"]
    with np.load(training_dataset / "labels.npz", allow_pickle=False) as source:
        y = source["y"]
    with np.load(training_results / "splits.npz", allow_pickle=False) as source:
        split = source["split"]
        if not np.array_equal(source["ids"], ids) or not np.array_equal(source["process"], process):
            raise ValueError("Saved split identities differ from the original training dataset")
    # Check the originally declared deterministic physical-event partition too.
    from evaluation import event_split
    if not np.array_equal(split, event_split(ids, process, int(metrics["protocol"]["split_seed"]))):
        raise ValueError("Saved split differs from the originally declared event split")
    names = ("reconstructed_mass_gev", "reconstructed_vertex_z_cm", "reconstructed_vertex_radius_cm")
    train = (split == 0) & (y == 1)
    edges = {name: np.unique(np.quantile(nuisance[train, i], [1/3, 2/3])).tolist() for i, name in enumerate(names)}
    return split_registry(training_results / "splits.npz"), edges, file_digest(manifest_path)


def wilson(passed, total):
    if not total:
        return [None, None]
    z = 1.959963984540054
    p, d = passed / total, 1 + z * z / total
    center = (p + z * z / (2 * total)) / d
    delta = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / d
    return [float(max(0., center - delta)), float(min(1., center + delta))]


def _group_ids(groups):
    # Do not convert uint64 event IDs to floats or use their string truncation.
    table, inverse = {}, []
    for key in groups:
        key = tuple(key)
        if key not in table:
            table[key] = len(table)
        inverse.append(table[key])
    return np.asarray(inverse, dtype=int), len(table)


def efficiency_summary(passed, groups, repeats=200, seed=271):
    passed = np.asarray(passed, dtype=bool)
    if len(passed) != len(groups):
        raise ValueError("Efficiency event identities are not aligned")
    inverse, n_events = _group_ids(groups)
    k, n = int(passed.sum()), len(passed)
    result = dict(pairs=n, passed_pairs=k, events=n_events, efficiency=k / n if n else None,
                  ci95_pair_wilson=wilson(k, n), ci95_event_bootstrap=[None, None],
                  bootstrap_replicates=0, bootstrap_informative=False)
    # A lone event cannot establish an event-to-event sampling uncertainty.
    if n_events < 2 or repeats <= 0:
        return result
    total = np.bincount(inverse, minlength=n_events)
    accepted = np.bincount(inverse, weights=passed, minlength=n_events)
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(repeats):
        counts = np.bincount(rng.integers(n_events, size=n_events), minlength=n_events)
        values.append(float((counts @ accepted) / (counts @ total)))
    result.update(ci95_event_bootstrap=np.quantile(values, [.025, .975]).tolist(),
                  bootstrap_replicates=len(values), bootstrap_informative=len(set(values)) > 1)
    return result


def weighted_efficiency(passed, groups, weights, repeats=200, seed=271):
    """Conditional weighted acceptance; signed MC ratios are not probabilities."""
    passed, weights = np.asarray(passed, bool), np.asarray(weights, float)
    if weights.shape != passed.shape or len(groups) != len(passed):
        raise ValueError("Weighted efficiency arrays are not aligned")
    inverse, n_events = _group_ids(groups)
    finite = np.isfinite(weights).all()
    total = float(weights.sum()) if finite else None
    accepted = float(weights[passed].sum()) if finite else None
    result = dict(pairs=len(passed), events=n_events, weight_sum=total, selected_weight_sum=accepted,
                  negative_weight_pairs=int((weights < 0).sum()), zero_weight_pairs=int((weights == 0).sum()),
                  efficiency=None, ci95_event_bootstrap=[None, None], bootstrap_replicates=0,
                  bootstrap_informative=False, available=False,
                  signed_mc_acceptance_ratio=accepted / total if total is not None and total != 0 else None)
    if not finite:
        result["reason"] = "Nonfinite evaluation weights; no weighted probability is defined"
    elif (weights < 0).any():
        result["reason"] = "Signed MC weights: raw signed acceptance ratio retained separately; probability efficiency and its bootstrap are unavailable"
    elif total <= 0:
        result["reason"] = "Zero total evaluation weight; conditional efficiency is unmeasured"
    else:
        result.update(available=True, efficiency=accepted / total,
                      pair_weight_effective_n=total * total / float(weights @ weights))
        event_total = np.bincount(inverse, weights=weights, minlength=n_events)
        event_accepted = np.bincount(inverse, weights=weights * passed, minlength=n_events)
        result["event_weight_effective_n"] = total * total / float(event_total @ event_total)
        if n_events >= 2 and repeats > 0:
            rng, values = np.random.default_rng(seed), []
            for _ in range(repeats):
                counts = np.bincount(rng.integers(n_events, size=n_events), minlength=n_events)
                denominator = counts @ event_total
                if denominator > 0:
                    values.append(float((counts @ event_accepted) / denominator))
            if values:
                result.update(ci95_event_bootstrap=np.quantile(values, [.025, .975]).tolist(),
                              bootstrap_replicates=len(values), bootstrap_informative=len(set(values)) > 1)
    return result


def roc_summary(y, score, groups, repeats=200, seed=319, minimum=20, weights=None):
    """Pair ROC, with ties treated identically and entire events bootstrapped."""
    from sklearn.metrics import auc, roc_curve
    y, score = np.asarray(y), np.asarray(score, dtype=float)
    if y.shape != score.shape or y.shape != (len(groups),) or not np.isin(y, (0, 1)).all() or not np.isfinite(score).all():
        raise ValueError("ROC needs aligned finite scores and explicit binary labels")
    neg, pos = int((y == 0).sum()), int((y == 1).sum())
    positive_events = len({tuple(g) for label, g in zip(y, groups) if label == 1})
    negative_events = len({tuple(g) for label, g in zip(y, groups) if label == 0})
    result = dict(negative_pairs=neg, positive_pairs=pos, auc=None, auc_ci95_event_bootstrap=[None, None],
                  negative_events=negative_events, positive_events=positive_events,
                  supported=neg >= minimum and pos >= minimum and positive_events >= minimum and negative_events >= minimum,
                  bootstrap_replicates=0,
                  fpr=[], tpr=[], thresholds=[])
    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        if weights.shape != y.shape:
            raise ValueError("ROC weights are not aligned")
        finite = np.isfinite(weights).all()
        result.update(weighted=True, available=False, negative_weight_pairs=int((weights < 0).sum()),
                      zero_weight_pairs=int((weights == 0).sum()),
                      positive_weight_sum=float(weights[y == 1].sum()) if finite else None,
                      negative_weight_sum=float(weights[y == 0].sum()) if finite else None)
        if not finite or (weights < 0).any():
            result.update(reason="Probability ROC unavailable for nonfinite or signed weights; weights were not clipped or replaced by absolute values", supported=False)
            return result
        if weights[y == 0].sum() <= 0 or weights[y == 1].sum() <= 0:
            result.update(reason="Positive total evaluation weight is required in both classes", supported=False)
            return result
        result["available"] = True
        contributing_events = [{tuple(g) for label, g, w in zip(y, groups, weights) if label == c and w > 0} for c in (0, 1)]
        result["supported"] = result["supported"] and all(len(g) >= minimum for g in contributing_events)
        inverse_for_weights, group_count = _group_ids(groups)
        for label, name in ((0, "negative"), (1, "positive")):
            event_weights = np.bincount(inverse_for_weights, weights=weights * (y == label), minlength=group_count)
            effective_events = float(event_weights.sum() ** 2 / (event_weights @ event_weights))
            result[name + "_weight_effective_events"] = effective_events
            result["supported"] = result["supported"] and effective_events >= minimum
    if not neg or not pos:
        result["reason"] = "Both evaluated classes are required; empty cells are unmeasured"
        return result
    fpr, tpr, thresholds = roc_curve(y, score, sample_weight=weights, drop_intermediate=False)
    result.update(auc=float(auc(fpr, tpr)), fpr=fpr.tolist(), tpr=tpr.tolist(),
                  thresholds=[float(t) if np.isfinite(t) else None for t in thresholds])
    inverse, n_events = _group_ids(groups)
    rng, values = np.random.default_rng(seed), []
    for _ in range(max(0, repeats)):
        counts = np.bincount(rng.integers(n_events, size=n_events), minlength=n_events)
        replica_weights = counts[inverse] if weights is None else weights * counts[inverse]
        if (replica_weights[y == 0].sum() > 0) and (replica_weights[y == 1].sum() > 0):
            f, t, _ = roc_curve(y, score, sample_weight=replica_weights, drop_intermediate=True)
            values.append(float(auc(f, t)))
    if values:
        result.update(auc_ci95_event_bootstrap=np.quantile(values, [.025, .975]).tolist(),
                      bootstrap_replicates=len(values))
    return result


def validate_manifest(manifest):
    if manifest.get("schema") != SCHEMA or manifest.get("sample_kind") != "simulation_only":
        raise ValueError("Evaluation requires an explicit simulation-only manifest")
    inputs = manifest.get("inputs", [])
    if not isinstance(inputs, list) or not inputs:
        raise ValueError("No simulation Nano inputs")
    paths, points = set(), {}
    for point in manifest.get("planned_signal_points", []):
        _register_point(points, point)
    for row in inputs:
        path = Path(row["path"])
        if not path.is_absolute() or path.suffix != ".root" or str(path.resolve()) in paths:
            raise ValueError("Inputs must be unique absolute ROOT paths")
        paths.add(str(path.resolve()))
        if row.get("process") not in ("qcd", "dy", "jpsi", "signal"):
            raise ValueError("Every input needs a declared process")
        if row["process"] == "signal":
            _register_point(points, row)
        if "mixture_scale" in row and (not np.isfinite(float(row["mixture_scale"])) or float(row["mixture_scale"]) < 0):
            raise ValueError("Mixture scales must be finite and nonnegative")
        if "weight_branch" in row and (not isinstance(row["weight_branch"], str) or not row["weight_branch"]):
            raise ValueError("Configured weight branch must be a nonempty name")
        expected = row.get("nano_sha256")
        if expected is not None and (len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected)):
            raise ValueError("Invalid optional Nano SHA-256")
    return points


def _register_point(points, row):
    point = row.get("point")
    if not isinstance(point, str) or not point or Path(point).name != point:
        raise ValueError("Invalid signal point name")
    values = tuple(float(row[k]) for k in POINT_FIELDS)
    if not all(np.isfinite(values)) or values[0] <= 0 or values[1] < 0 or values[2] <= 0:
        raise ValueError("Invalid physical mass, proper c-tau or epsilon")
    if point in points and tuple(points[point][k] for k in POINT_FIELDS) != values:
        raise ValueError("One signal point has inconsistent physical parameters")
    points[point] = {k: value for k, value in zip(POINT_FIELDS, values)}
    if "target_mean_lab_flight_m" in row:
        points[point]["target_mean_lab_flight_m"] = row["target_mean_lab_flight_m"]


def _parent_ancestor(pdg, mothers, index, parent_pdg_id):
    found, seen = -1, set()
    while 0 <= index < len(pdg):
        if index in seen:
            return -1
        seen.add(index)
        if abs(int(pdg[index])) == abs(parent_pdg_id):
            found = index
        index = int(mothers[index])
    return found if index == -1 else -1


def signal_pairs(arrays, reco, parent_pdg_id=32):
    """MC association only, independent of stored scores and row retention."""
    result = np.zeros(len(reco["event_index"]), dtype=bool)
    for row, (e, muons) in enumerate(zip(reco["event_index"], reco["muon_indices"])):
        pdg, mothers = arrays["GenPart_pdgId"][e], arrays["GenPart_genPartIdxMother"][e]
        g = [int(arrays["ShiftMuon_hitGenPartIdx"][e][i]) for i in muons]
        if all(0 <= i < len(pdg) for i in g):
            parent = [_parent_ancestor(pdg, mothers, i, parent_pdg_id) for i in g]
            result[row] = parent[0] >= 0 and parent[0] == parent[1] and {int(pdg[i]) for i in g} == {13, -13}
    return result


def read_inputs(manifest, model, model_sha256, allow_unscored=False, chunk_events=1000):
    import awkward as ak
    import uproot
    rows, sources, seen = [], [], set()
    for source in manifest["inputs"]:
        path = Path(source["path"])
        original_stat = path.stat()
        source_signature = (original_stat.st_dev, original_stat.st_ino, original_stat.st_size, original_stat.st_mtime_ns)
        score_receipt_path = Path(source["score_receipt"]) if "score_receipt" in source else None
        if score_receipt_path is not None:
            score_receipt = json.loads(score_receipt_path.read_text())
            if (score_receipt.get("schema") != "shift-dimuon-bdt-score-receipt-v1" or
                    score_receipt.get("complete") is not True or score_receipt.get("model_sha256") != model_sha256 or
                    Path(score_receipt.get("output", "")).resolve() != path.resolve()):
                raise ValueError("Completed scoring receipt does not identify this output and model")
            source_sha = score_receipt["output_sha256"]
        else:
            source_sha = file_digest(path)
        if source.get("nano_sha256", source_sha) != source_sha:
            raise ValueError("Nano content differs from frozen input: " + str(path))
        summary = dict(path=str(path), process=source["process"], nano_sha256=source_sha, events=0, pairs=0)
        if score_receipt_path is not None:
            summary.update(score_receipt=str(score_receipt_path.resolve()), score_receipt_sha256=file_digest(score_receipt_path),
                           integrity_reference="completed preservation/scoring receipt; inode/size/mtime guarded during leaf extraction")
        else:
            summary["integrity_reference"] = "one evaluator SHA-256 read; inode/size/mtime guarded during leaf extraction"
        # Use the production scorer's bounded local-file reader; the default
        # fsspec async reader can stall on this AFS/CERN runtime combination.
        with uproot.open(path, handler=uproot.source.file.MultithreadedFileSource, num_workers=1) as root:
            tree = root["Events"]
            scored = SCORE_BRANCH in tree.keys()
            if not scored and not allow_unscored:
                raise ValueError("BDT score branch absent: " + str(path))
            if scored:
                metadata = json.loads(str(root[SCORE_BRANCH + "Metadata"]))
                if metadata.get("model_sha256") != model_sha256 or metadata.get("selection_applied") is not False:
                    raise ValueError("Stored score model/selection metadata differs: " + str(path))
            weight_branch = source.get("weight_branch", "genWeight" if source["process"] == "signal" else "shiftSamplingGenWeight")
            weight_present = weight_branch in tree.keys()
            if "weight_branch" in source and not weight_present:
                raise ValueError("Explicit evaluation weight branch absent: " + weight_branch)
            mixture_scale = source.get("mixture_scale", 1. if source["process"] == "signal" else None)
            weighted_available = weight_present and mixture_scale is not None
            wanted = tuple(dict.fromkeys(PREDICTOR_BRANCHES + TRUTH_BRANCHES +
                                        ("ShiftDimuonVertex_mass",) + ((SCORE_BRANCH,) if scored else ()) +
                                        ((weight_branch,) if weight_present else ())))
            missing = set(wanted) - set(tree.keys())
            if missing:
                raise ValueError("Evaluation branches absent: " + repr(sorted(missing)))
            summary["score_origin"] = "stored_checked_metadata" if scored else "portable_reco_inference_for_evaluation"
            summary.update(evaluation_weight_branch=weight_branch, mixture_scale=mixture_scale,
                           weighted_evaluation_available=weighted_available)
            for arrays in tree.iterate(wanted, step_size=chunk_events, library="ak", how=dict):
                # No truth is read by extract_reco or predict_score.
                reco = extract_reco({k: arrays[k] for k in PREDICTOR_BRANCHES}, include_nuisance=False)
                # Mixed signed/unsigned integer column_stack can promote to
                # Float64, losing event IDs above 2**53. Cast each field first.
                event_ids = np.column_stack([ak.to_numpy(arrays[k]).astype(np.uint64) for k in ("run", "luminosityBlock", "event")])
                for event in event_ids:
                    key = identity(source["process"] if source["process"] != "signal" else source["point"], event)
                    if key in seen:
                        raise ValueError("Duplicate physical event across input Nano files")
                    seen.add(key)
                calculated = predict_score(model, reco["X"])
                score = ak.to_numpy(ak.flatten(arrays[SCORE_BRANCH])) if scored else calculated
                if score.shape != calculated.shape or not np.isfinite(score).all() or not np.array_equal(score, calculated):
                    raise ValueError("Stored scores disagree with reconstructed-only portable inference")
                y, reason = extract_labels(arrays, reco)
                evaluation_weights = (ak.to_numpy(arrays[weight_branch])[reco["event_index"]] * float(mixture_scale)
                                      if weighted_available else np.zeros(len(y), dtype=float))
                signal = signal_pairs(arrays, reco) if source["process"] == "signal" else np.zeros(len(y), bool)
                nuisance = np.column_stack([ak.to_numpy(ak.flatten(arrays["ShiftDimuonVertex_mass"])),
                                           ak.to_numpy(ak.flatten(arrays["ShiftDimuonVertex_vz"])),
                                           np.hypot(ak.to_numpy(ak.flatten(arrays["ShiftDimuonVertex_vx"])),
                                                    ak.to_numpy(ak.flatten(arrays["ShiftDimuonVertex_vy"])))])
                rows.append(dict(ids=reco["ids"], score=score, y=y, reason=reason, signal=signal, nuisance=nuisance,
                                 evaluation_weight=evaluation_weights,
                                 evaluation_weight_available=np.full(len(y), weighted_available, dtype=bool),
                                 process=np.full(len(y), source["process"]), point=np.full(len(y), source.get("point", ""))))
                summary["events"] += len(event_ids)
                summary["pairs"] += len(y)
        final_stat = path.stat()
        if (final_stat.st_dev, final_stat.st_ino, final_stat.st_size, final_stat.st_mtime_ns) != source_signature:
            raise ValueError("Input changed during performance extraction")
        sources.append(summary)
        print(json.dumps(dict(source_complete=str(path), events=summary["events"], pairs=summary["pairs"])), flush=True)
    keys = ("ids", "score", "y", "reason", "signal", "nuisance", "process", "point", "evaluation_weight", "evaluation_weight_available")
    return {k: np.concatenate([row[k] for row in rows]) for k in keys}, sources


def _groups(data, mask):
    return [identity(p if p != "signal" else point, event) for p, point, event in
            zip(data["process"][mask], data["point"][mask], data["ids"][mask])]


def _eff(data, mask, threshold, repeats):
    result = efficiency_summary(data["score"][mask] >= threshold, _groups(data, mask), repeats)
    if "evaluation_weight_available" in data and data["evaluation_weight_available"][mask].all():
        weighted = weighted_efficiency(data["score"][mask] >= threshold, _groups(data, mask), data["evaluation_weight"][mask], repeats)
    else:
        weighted = dict(available=False, efficiency=None, ci95_event_bootstrap=[None, None],
                        reason="Evaluation weight branch or explicit SM mixture_scale unavailable")
    result["physics_mixture_weighted"] = weighted
    return result


def _roc(data, positive, negative, repeats, minimum):
    mask = positive | negative
    result = roc_summary(positive[mask].astype(np.int8), data["score"][mask], _groups(data, mask), repeats, minimum=minimum)
    if "evaluation_weight_available" in data and data["evaluation_weight_available"][mask].all():
        weighted = roc_summary(positive[mask].astype(np.int8), data["score"][mask], _groups(data, mask), repeats,
                               minimum=minimum, weights=data["evaluation_weight"][mask])
    else:
        weighted = dict(available=False, auc=None, auc_ci95_event_bootstrap=[None, None], fpr=[], tpr=[], thresholds=[],
                        reason="Evaluation weight branch or explicit SM mixture_scale unavailable")
    result["physics_mixture_weighted"] = weighted
    return result


def evaluate(data, sources, points, registry, threshold, repeats=200, minimum=20, training_edges=None):
    holdout = heldout_mask(data["process"], data["ids"], registry)
    sm = np.isin(data["process"], ("qcd", "dy", "jpsi"))
    known, genuine, accidental = data["y"] >= 0, data["y"] == 1, data["y"] == 0
    qcd_negative = holdout & (data["process"] == "qcd") & accidental
    report = dict(threshold=threshold, threshold_origin="frozen SM validation; unchanged for every process/point",
                  evaluation_weighting="Central empirical summaries use unweighted reconstructed pairs and their selected MC mixture. Explicit physics_mixture_weighted companions use the declared cross sections and full-parent-GEN sumw, without an MC-to-data fit.",
                  weighted_companion="Evaluation-only weights: SM shiftSamplingGenWeight times explicit sigma_pb / full-parent-GEN sumw; each signal point uses its native genWeight. No N_eff is fitted. Signed weights retain their raw efficiency ratio but have no probability ROC.",
                  mc_to_data_normalization_measured=False, mass_independence_validated=False,
                  lifetime_independence_validated=False,
                  support_minimum_pairs_per_class=minimum, event_bootstrap_repeats=repeats,
                  support_minimum_independent_events_per_class=minimum,
                  uncertainty_note="Pairs from the same physical event are resampled together. Wilson intervals additionally expose finite counts; saturated/small bootstrap intervals do not establish closure.",
                  heldout=dict(excluded_training_or_validation_pairs=int((sm & ~holdout).sum()),
                               saved_test_pairs=int(sum(registry.get(identity(p, event), -1) == 2 for p, event in zip(data["process"], data["ids"]))),
                               unseen_pairs=int(sum(identity(p, event) not in registry for p, event in zip(data["process"][sm], data["ids"][sm]))),
                               selection="saved test plus events absent from the entire saved fit/selection registry"),
                  label_counts=dict(Counter(data["reason"].tolist())), roc={}, per_process={}, signal_points={})
    report["roc"]["heldout_genuine_vs_all_accidental"] = _roc(data, sm & holdout & genuine, sm & holdout & accidental, repeats, minimum)
    report["roc"]["heldout_genuine_vs_qcd_accidental"] = _roc(data, sm & holdout & genuine, qcd_negative, repeats, minimum)
    saved_test = np.asarray([registry.get(identity(p, event), -1) == 2 for p, event in zip(data["process"], data["ids"])])
    report["roc"]["original_saved_test_genuine_vs_accidental"] = _roc(data, sm & saved_test & genuine, sm & saved_test & accidental, repeats, minimum)
    report["roc"]["full_sample_descriptive_genuine_vs_accidental"] = _roc(data, sm & genuine, sm & accidental, repeats, minimum)
    report["roc"]["full_sample_descriptive_genuine_vs_accidental"]["warning"] = "Includes fitted/validation events; descriptive, not held-out performance"
    for process in ("jpsi", "dy", "qcd"):
        mask = sm & holdout & (data["process"] == process)
        report["per_process"][process] = dict(genuine_retention=_eff(data, mask & genuine, threshold, repeats),
                                              accidental_retention=_eff(data, mask & accidental, threshold, repeats),
                                              unknown_retention=_eff(data, mask & ~known, threshold, repeats),
                                              unknown_labels=int((mask & ~known).sum()))
        report["roc"][process + "_genuine_vs_qcd_accidental"] = _roc(data, mask & genuine, qcd_negative, repeats, minimum)
    for point, metadata in points.items():
        mask = (data["process"] == "signal") & (data["point"] == point)
        positive = mask & genuine & data["signal"]
        entry = dict(metadata=metadata, source_files=sum(s.get("process") == "signal" and s.get("point") == point for s in sources),
                     recorded_nano_events=sum(s["events"] for s in sources if s.get("process") == "signal" and s.get("point") == point),
                     all_retained_pairs=int(mask.sum()), label_counts=dict(Counter(data["reason"][mask].tolist())),
                     common_signal_pair_efficiency=_eff(data, positive, threshold, repeats),
                     common_vertex_pair_efficiency=_eff(data, mask & genuine, threshold, repeats),
                     unknown_pair_retention=_eff(data, mask & ~known, threshold, repeats),
                     roc_vs_heldout_qcd_accidental=_roc(data, positive, qcd_negative, repeats, minimum),
                     interpretation="BDT acceptance conditional on already reconstructed, truth-associated genuine signal pairs; zero-reconstruction events are not BDT inefficiency")
        entry["supported"] = entry["roc_vs_heldout_qcd_accidental"]["supported"]
        if not mask.any():
            entry["status"] = "unmeasured: no retained pairs or no completed Nano input"
        elif not positive.any():
            entry["status"] = "unmeasured: no labelled genuine signal pairs"
        else:
            entry["status"] = "supported" if entry["supported"] else "small sample; uniformity unvalidated"
        report["signal_points"][point] = entry
    # These measure reconstructed observables, never the generated lifetime.
    profiles = {}
    for index, name in enumerate(("reconstructed_mass_gev", "reconstructed_vertex_z_cm", "reconstructed_vertex_radius_cm")):
        original_training = np.asarray([registry.get(identity(p, event), -1) == 0 for p, event in zip(data["process"], data["ids"])])
        training = sm & genuine & original_training & np.isfinite(data["nuisance"][:, index])
        values = data["nuisance"][training, index]
        if training_edges is not None and name in training_edges:
            edges = np.asarray(training_edges[name])
        elif not len(values):
            profiles[name] = dict(status="No frozen training support; bin edges unavailable")
            continue
        else:
            edges = np.unique(np.quantile(values, [1/3, 2/3]))
        bins = np.searchsorted(edges, data["nuisance"][:, index], side="right")
        entries = []
        for b in range(len(edges) + 1):
            mask = sm & holdout & genuine & (bins == b) & np.isfinite(data["nuisance"][:, index])
            entry = dict(bin=b, lower=float(edges[b-1]) if b else None,
                         upper=float(edges[b]) if b < len(edges) else None, **_eff(data, mask, threshold, repeats))
            entry["supported"] = entry["pairs"] >= minimum and entry["events"] >= minimum
            entry["per_process"] = {p: _eff(data, mask & (data["process"] == p), threshold, repeats) for p in ("jpsi", "dy", "qcd")}
            entries.append(entry)
        profiles[name] = dict(interior_edges=edges.tolist(), bins=entries, bin_definition="tertiles of original training events; never held-out scores")
    report["reconstructed_acceptance_profiles"] = profiles
    report["reconstructed_mass_radius_map"] = []
    mass_profile, radius_profile = profiles["reconstructed_mass_gev"], profiles["reconstructed_vertex_radius_cm"]
    if "interior_edges" in mass_profile and "interior_edges" in radius_profile:
        mass_bins = np.searchsorted(mass_profile["interior_edges"], data["nuisance"][:, 0], side="right")
        radius_bins = np.searchsorted(radius_profile["interior_edges"], data["nuisance"][:, 2], side="right")
        for m in range(len(mass_profile["interior_edges"]) + 1):
            for r in range(len(radius_profile["interior_edges"]) + 1):
                mask = sm & holdout & genuine & (mass_bins == m) & (radius_bins == r) & np.isfinite(data["nuisance"][:, [0, 2]]).all(axis=1)
                entry = dict(mass_bin=m, radius_bin=r, genuine_pair_efficiency=_eff(data, mask, threshold, repeats),
                             roc_vs_qcd_accidental=_roc(data, mask, qcd_negative, repeats, minimum),
                             per_process={p: _eff(data, mask & (data["process"] == p), threshold, repeats) for p in ("jpsi", "dy", "qcd")})
                entry["supported"] = entry["genuine_pair_efficiency"]["pairs"] >= minimum and entry["genuine_pair_efficiency"]["events"] >= minimum
                report["reconstructed_mass_radius_map"].append(entry)
    return report


def make_plots(report, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    def save(fig, name):
        fig.savefig(output / (name + ".pdf"), bbox_inches="tight")
        plt.close(fig)
    def support_note(roc):
        if roc.get("supported") is True:
            return ""
        pair_min = report.get("support_minimum_pairs_per_class")
        event_min = report.get("support_minimum_independent_events_per_class")
        limits = []
        for prefix, label in (("positive", "genuine"), ("negative", "QCD")):
            pairs, events = roc.get(prefix + "_pairs"), roc.get(prefix + "_events")
            if pair_min is not None and pairs is not None and pairs < pair_min:
                limits.append(f"{label} pairs {pairs} < {pair_min}")
            if event_min is not None and events is not None and events < event_min:
                limits.append(f"{label} events {events} < {event_min}")
            effective = roc.get(prefix + "_weight_effective_events")
            if event_min is not None and effective is not None and effective < event_min:
                limits.append(f"{label} effective events {effective:.2f} < {event_min}")
        return "; ".join(limits) or "insufficient statistical support"

    sm_titles = (("heldout_genuine_vs_qcd_accidental", "All SM genuine"),
                 ("jpsi_genuine_vs_qcd_accidental", "J/psi genuine"),
                 ("dy_genuine_vs_qcd_accidental", "DY genuine"))
    for weighted in (False, True):
        for signal_plot in (False, True):
            fig, axis = plt.subplots(figsize=(7, 6))
            curves = []
            if signal_plot:
                for point, entry in report["signal_points"].items():
                    metadata = entry["metadata"]
                    title = f"m = {metadata['mass_gev']:g} GeV, proper cτ = {metadata['proper_ctau_mm']:.3g} mm"
                    curves.append((title, entry["roc_vs_heldout_qcd_accidental"]))
            else:
                curves = [(title, report["roc"][name]) for name, title in sm_titles]
            unsupported, unavailable = False, []
            for title, comparison in curves:
                roc = comparison["physics_mixture_weighted"] if weighted else comparison
                if roc["auc"] is None:
                    unavailable.append(title)
                    continue
                note = support_note(roc)
                unsupported = unsupported or bool(note)
                lo, hi = roc["auc_ci95_event_bootstrap"]
                interval = f", 95% [{lo:.3f}, {hi:.3f}]" if lo is not None else ""
                label = f"{title}\nAUC {roc['auc']:.3f}{interval}"
                if note:
                    # Keep counts in the legend; styling alone must not suggest
                    # that an unsupported AUC is an established measurement.
                    import textwrap
                    label += "\n" + textwrap.fill(note, width=65)
                axis.plot(roc["fpr"], roc["tpr"], linestyle="--" if note else "-", label=label)
            axis.plot([0, 1], [0, 1], "k:", linewidth=1)
            prefix = "Physics-mixture weighted" if weighted else "Empirical unweighted"
            title = prefix + (" signal" if signal_plot else " held-out SM") + " BDT ROC"
            if unsupported:
                title += "\nDashed curves: insufficient statistical support"
            axis.set(xlabel=("Weighted " if weighted else "") + "QCD accidental-pair acceptance",
                     ylabel=("Weighted " if weighted else "") + "Genuine-pair acceptance",
                     xlim=(0, 1), ylim=(0, 1), title=title)
            if axis.get_legend_handles_labels()[0]:
                axis.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8,
                            borderaxespad=0, frameon=False, labelspacing=1.)
            else:
                axis.text(.5, .2, "Probability ROC unavailable\nSee weights and class-support audit", ha="center")
            if unavailable:
                axis.text(0, -.18, "Unmeasured: " + "; ".join(unavailable),
                          transform=axis.transAxes, fontsize=8, wrap=True)
            axis.grid(alpha=.2)
            name = "signal_roc" if signal_plot else "heldout_roc"
            save(fig, name + ("_physics_mixture_weighted" if weighted else ""))
    for metric, title, weighted in (("efficiency", "Fixed-cut genuine signal-pair efficiency", False),
                                    ("auc", "Signal vs QCD accidental-pair AUC", False),
                                    ("efficiency", "Weighted fixed-cut genuine signal-pair efficiency", True),
                                    ("auc", "Weighted signal vs QCD accidental-pair AUC", True)):
        fig, axis = plt.subplots(figsize=(8, 6))
        colored = None
        maximum_mass = max((entry["metadata"]["mass_gev"] for entry in report["signal_points"].values()), default=0)
        for point, entry in report["signal_points"].items():
            m, ctau = entry["metadata"]["mass_gev"], entry["metadata"]["proper_ctau_mm"]
            y = max(ctau, 1.e-9)
            eff = entry["common_signal_pair_efficiency"]
            eff_values = eff["physics_mixture_weighted"] if weighted else eff
            roc_values = entry["roc_vs_heldout_qcd_accidental"]
            if weighted:
                roc_values = roc_values["physics_mixture_weighted"]
            value = eff_values["efficiency"] if metric == "efficiency" else roc_values["auc"]
            if value is None:
                axis.scatter([m], [y], marker="x", c="gray", s=130)
                label = "unmeasured"
            else:
                colored = axis.scatter([m], [y], c=[value], cmap="viridis", vmin=0, vmax=1, s=330,
                                       edgecolor="black" if (entry["supported"] and (not weighted or roc_values.get("supported", False))) else "red", linewidth=1.4)
                interval = (eff_values["ci95_event_bootstrap"] if weighted else eff["ci95_pair_wilson"]) if metric == "efficiency" else roc_values["auc_ci95_event_bootstrap"]
                label = f"{value:.3f}; n={eff['pairs']}"
                if interval[0] is not None:
                    label += f"\n95% [{interval[0]:.3f}, {interval[1]:.3f}]"
            right = m == maximum_mass
            axis.annotate(f"{point}\n{label}", (m, y), xytext=(-8 if right else 8, 8),
                          ha="right" if right else "left", textcoords="offset points", fontsize=7)
        axis.set(xlabel="Generated signal mass [GeV]", ylabel="Generated proper cτ [mm]", yscale="log",
                 title=title + "\nRed borders: insufficient pairs or independent events")
        axis.margins(x=.25, y=.3)
        axis.grid(alpha=.2)
        if colored is not None:
            fig.colorbar(colored, ax=axis, label=metric)
        save(fig, "signal_mass_proper_ctau_" + metric + ("_physics_mixture_weighted" if weighted else ""))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for axis, (name, profile) in zip(axes, report["reconstructed_acceptance_profiles"].items()):
        for process, color in ((None, "black"), ("jpsi", "tab:blue"), ("dy", "tab:orange"), ("qcd", "tab:green")):
            for entry in profile.get("bins", []):
                item = entry if process is None else entry["per_process"][process]
                if item["efficiency"] is not None:
                    value, interval = item["efficiency"], item["ci95_pair_wilson"]
                    axis.errorbar(entry["bin"] + {None: -.12, "jpsi": -.04, "dy": .04, "qcd": .12}[process], value,
                                  yerr=np.array([[value-interval[0]], [interval[1]-value]]), fmt="o", color=color,
                                  label=(process or "all SM") if entry["bin"] == 0 else None, capsize=2)
        axis.set(title=name.replace("reconstructed_", "").replace("_", " "), xlabel="Frozen training tertile", ylabel="Genuine-pair retention", ylim=(0, 1.05))
        axis.grid(alpha=.2)
    if axes[0].get_legend_handles_labels()[0]:
        axes[0].legend(fontsize=8)
    fig.suptitle(f"Held-out fixed-cut acceptance (score ≥ {report['threshold']:.8f}); 95% Wilson intervals")
    save(fig, "reconstructed_acceptance_profiles")
    cells = report["reconstructed_mass_radius_map"]
    if cells:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        shape = (max(c["radius_bin"] for c in cells) + 1, max(c["mass_bin"] for c in cells) + 1)
        for axis, metric in zip(axes, ("efficiency", "auc")):
            values = np.full(shape, np.nan)
            for cell in cells:
                value = cell["genuine_pair_efficiency"]["efficiency"] if metric == "efficiency" else cell["roc_vs_qcd_accidental"]["auc"]
                if value is not None:
                    values[cell["radius_bin"], cell["mass_bin"]] = value
            rendered = axis.imshow(values, origin="lower", cmap="viridis", vmin=0, vmax=1)
            for cell in cells:
                value = values[cell["radius_bin"], cell["mass_bin"]]
                n = cell["genuine_pair_efficiency"]["pairs"]
                label = f"{value:.3f}\nn={n}" if np.isfinite(value) else "unmeasured"
                axis.text(cell["mass_bin"], cell["radius_bin"], label, ha="center", va="center", fontsize=9,
                          color="black" if value > .6 else "white")
            axis.set(xlabel="Reconstructed mass tertile", ylabel="Reconstructed radius tertile", title=metric,
                     xticks=range(shape[1]), yticks=range(shape[0]))
            fig.colorbar(rendered, ax=axis)
        fig.suptitle("Held-out SM performance; reconstructed radius is not proper lifetime")
        save(fig, "reconstructed_mass_radius_map")


def write_tables(report, output):
    columns = ["point", "mass_gev", "proper_ctau_mm", "epsilon", "pairs", "events", "passed_pairs", "efficiency", "efficiency_ci95_low", "efficiency_ci95_high", "efficiency_event_bootstrap_ci95_low", "efficiency_event_bootstrap_ci95_high", "efficiency_bootstrap_informative", "auc", "auc_ci95_low", "auc_ci95_high", "weighted_efficiency", "weighted_efficiency_ci95_low", "weighted_efficiency_ci95_high", "signed_mc_acceptance_ratio", "weighted_auc", "weighted_auc_ci95_low", "weighted_auc_ci95_high", "supported", "status"]
    with (output / "signal_mass_lifetime.csv").open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=columns)
        writer.writeheader()
        for point, entry in report["signal_points"].items():
            eff, roc = entry["common_signal_pair_efficiency"], entry["roc_vs_heldout_qcd_accidental"]
            weighted_eff, weighted_roc = eff["physics_mixture_weighted"], roc["physics_mixture_weighted"]
            writer.writerow(dict(point=point, **{k: entry["metadata"][k] for k in POINT_FIELDS},
                                 **{k: eff[k] for k in ("pairs", "events", "passed_pairs", "efficiency")},
                                 efficiency_ci95_low=eff["ci95_pair_wilson"][0], efficiency_ci95_high=eff["ci95_pair_wilson"][1],
                                 efficiency_event_bootstrap_ci95_low=eff["ci95_event_bootstrap"][0],
                                 efficiency_event_bootstrap_ci95_high=eff["ci95_event_bootstrap"][1],
                                 efficiency_bootstrap_informative=eff["bootstrap_informative"],
                                 auc=roc["auc"], auc_ci95_low=roc["auc_ci95_event_bootstrap"][0], auc_ci95_high=roc["auc_ci95_event_bootstrap"][1],
                                 weighted_efficiency=weighted_eff["efficiency"],
                                 weighted_efficiency_ci95_low=weighted_eff["ci95_event_bootstrap"][0],
                                 weighted_efficiency_ci95_high=weighted_eff["ci95_event_bootstrap"][1],
                                 signed_mc_acceptance_ratio=weighted_eff.get("signed_mc_acceptance_ratio"),
                                 weighted_auc=weighted_roc["auc"], weighted_auc_ci95_low=weighted_roc["auc_ci95_event_bootstrap"][0],
                                 weighted_auc_ci95_high=weighted_roc["auc_ci95_event_bootstrap"][1],
                                 supported=entry["supported"], status=entry["status"]))
    with (output / "roc_curves.csv").open("w", newline="") as target:
        writer = csv.writer(target)
        writer.writerow(["comparison", "false_positive_rate", "true_positive_rate", "score_threshold"])
        curves = dict(report["roc"])
        curves.update({point: entry["roc_vs_heldout_qcd_accidental"] for point, entry in report["signal_points"].items()})
        for name, roc in curves.items():
            writer.writerows((name, *values) for values in zip(roc["fpr"], roc["tpr"], roc["thresholds"]))
            weighted = roc["physics_mixture_weighted"]
            writer.writerows((name + "_physics_mixture_weighted", *values) for values in zip(weighted["fpr"], weighted["tpr"], weighted["thresholds"]))


def run(manifest_path, model_path, training_results, output, repeats=200, minimum=20, allow_unscored=False, training_dataset=None):
    manifest_sha256, model_sha256 = file_digest(manifest_path), file_digest(model_path)
    manifest = json.loads(manifest_path.read_text())
    points = validate_manifest(manifest)
    model = load_model(model_path)
    lock_path = training_results / "selection_locked.json"
    if file_digest(lock_path) != model["provenance"]["selection_lock_sha256"]:
        raise ValueError("Training selection lock differs from portable model provenance")
    if training_dataset is None:
        training_dataset = training_results.parent / "sm_v10_phase2_fresh"
    registry, training_edges, training_dataset_manifest_sha256 = training_reference(training_results, training_dataset)
    output.mkdir(parents=True, exist_ok=False)
    data, sources = read_inputs(manifest, model, model_sha256, allow_unscored)
    # Store source point identities too, including sources that reconstruct no pairs.
    for summary, source in zip(sources, manifest["inputs"]):
        if "point" in source:
            summary["point"] = source["point"]
    report = evaluate(data, sources, points, registry, float(model["threshold_metadata"]), repeats, minimum, training_edges)
    report.update(schema="shift-bdt-performance-v1", complete=False,
                  input_manifest=str(manifest_path.resolve()), input_manifest_sha256=manifest_sha256,
                  model_sha256=model_sha256, training_split_sha256=file_digest(training_results / "splits.npz"),
                  training_dataset_manifest_sha256=training_dataset_manifest_sha256,
                  source_sha256={name: file_digest(Path(__file__).with_name(name)) for name in
                                 ("bdt_performance.py", "features.py", "labels.py", "portable_inference.py")},
                  sources=sources, selection_applied=False,
                  signal_axis_note="Generated proper cτ in mm. Target mean laboratory flight, if provided, is a finite-reference boost estimate, not proper lifetime or reconstructed radius.")
    np.savez_compressed(output / "evaluation_rows.npz", **data,
                        heldout=heldout_mask(data["process"], data["ids"], registry))
    (output / "performance.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    make_plots(report, output)
    write_tables(report, output)
    if file_digest(manifest_path) != manifest_sha256 or file_digest(model_path) != model_sha256:
        raise ValueError("Manifest or model changed while evaluating")
    report["complete"] = True
    report["artifact_sha256"] = {path.name: file_digest(path) for path in output.iterdir() if path.name != "performance.json" and path.is_file()}
    (output / "performance.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(dict(complete=True, retained_pairs=len(data["score"]), output=str(output))), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--training-results", type=Path, required=True)
    parser.add_argument("--training-dataset", type=Path, help="Original digest-locked training dataset; defaults to training-results/../sm_v10_phase2_fresh")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-repeats", type=int, default=200)
    parser.add_argument("--minimum-pairs", type=int, default=20)
    parser.add_argument("--allow-unscored", action="store_true", help="Explicitly compute reco-only scores for evaluation-only sources without a stored BDT branch")
    args = parser.parse_args()
    if args.bootstrap_repeats < 0 or args.minimum_pairs < 1:
        parser.error("Nonnegative bootstrap repeats and positive support minimum required")
    run(args.manifest, args.model, args.training_results, args.output, args.bootstrap_repeats, args.minimum_pairs, args.allow_unscored, args.training_dataset)


if __name__ == "__main__":
    main()
