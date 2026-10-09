#!/usr/bin/env python3
"""Audit MC-only label completeness using reconstructed population diagnostics.

Truth is read only from labels.npz to check training-label conventions. No truth
columns are appended to X, and no classifier or adversary is fitted here.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from evaluation import wilson
from features import FEATURE_NAMES


DETAIL_SHAPES = {
    "gen_index": (2,), "sim_track_id": (2,), "match_purity": (2,),
    "matched_layers": (2,), "match_state": (2,),
    "vertex_distance_cm": (), "same_gen_muon": (), "same_sim_track": (),
    "same_immediate_mother": (),
}
VALID_GEN_STATES = ("ok", "nonfinite_purity", "low_purity", "too_few_layers")
LABEL_NAMES = {-1: "unknown", 0: "accidental", 1: "common_vertex"}
CONTINUOUS_FEATURES = (
    "mu_min_nMuonStations", "mu_min_nValidMuonHits", "mu_min_nValidHits",
    "mu_min_normalizedChi2", "mu_max_normalizedChi2", "pair_normalizedChi2",
    "pair_dca", "opening_angle", "mu_vertex_distance", "pair_vz", "pair_pt",
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_dataset(path):
    path = Path(path)
    manifest = json.loads((path / "inputs.json").read_text())
    if manifest.get("sample_kind") != "simulation_only":
        raise ValueError("Label audits may only read explicitly identified simulation")
    for filename, key in (("reco.npz", "reco_sha256"), ("labels.npz", "labels_sha256")):
        if digest(path / filename) != manifest.get(key):
            raise ValueError("Dataset digest changed: " + filename)
    with np.load(path / "reco.npz", allow_pickle=False) as archive:
        reco = {key: archive[key] for key in archive.files}
    with np.load(path / "labels.npz", allow_pickle=False) as archive:
        labels = {key: archive[key] for key in archive.files}
    if tuple(reco["feature_names"]) != FEATURE_NAMES:
        raise ValueError("Stored reconstructed feature contract differs")
    y = labels["y"]
    n = len(y)
    shapes = {"X": (n, len(FEATURE_NAMES)), "ids": (n, 3), "process": (n,),
              "stratum": (n,), "pair_index": (n,), "nuisance": (n, 3)}
    for key, shape in shapes.items():
        if reco[key].shape != shape:
            raise ValueError("Invalid reconstructed array shape: " + key)
    if y.shape != (n,) or labels["reason"].shape != (n,) or not np.isin(y, (-1, 0, 1)).all():
        raise ValueError("Invalid label or label-reason shape/value")
    present = set(DETAIL_SHAPES).intersection(labels)
    if present and present != set(DETAIL_SHAPES):
        raise ValueError("Partial label diagnostics: missing " + ", ".join(sorted(set(DETAIL_SHAPES) - present)))
    for key in present:
        if labels[key].shape != (n,) + DETAIL_SHAPES[key]:
            raise ValueError("Invalid label diagnostic shape: " + key)
    return manifest, reco, labels


def completeness(y, mask=None):
    if mask is None:
        mask = np.ones(len(y), dtype=bool)
    n = int(mask.sum())
    unknown = int((mask & (y < 0)).sum())
    return dict(n=n, common_vertex=int((mask & (y == 1)).sum()),
                accidental=int((mask & (y == 0)).sum()), unknown=unknown,
                unknown_fraction=unknown / n if n else None,
                unknown_ci95_pair_wilson=wilson(unknown, n))


def grouped(values, y):
    return {str(value): completeness(y, values == value) for value in np.unique(values)}


def category_pairs(reco, y, kind):
    names = list(reco["feature_names"])
    values = reco["X"][:, [names.index("mu_min_" + kind), names.index("mu_max_" + kind)]]
    keys = np.array(["/".join(format(float(x), ".6g") for x in pair) for pair in values])
    result = grouped(keys, y)
    for key, entry in result.items():
        entry["per_process"] = {
            str(process): completeness(y, (keys == key) & (reco["process"] == process))
            for process in np.unique(reco["process"])
        }
    return result


def distribution(values):
    finite = np.isfinite(values)
    return dict(n=int(len(values)), nonfinite=int((~finite).sum()),
                quantiles_10_50_90=np.quantile(values[finite], [0.1, 0.5, 0.9]).tolist()
                if finite.any() else [None, None, None])


def numeric_audit(values, y, process):
    """Audit-only quantiles use all reconstructed rows; no training bins change."""
    finite = np.isfinite(values)
    edges = np.unique(np.quantile(values[finite], [0.25, 0.5, 0.75])) if finite.any() else np.array([])
    bins = np.searchsorted(edges, values, side="right")
    entries = []
    for index in range(len(edges) + 1):
        mask = finite & (bins == index)
        entry = completeness(y, mask)
        entry.update(bin=index, lower=float(edges[index - 1]) if index else None,
                     upper=float(edges[index]) if index < len(edges) else None,
                     per_process={str(p): completeness(y, mask & (process == p)) for p in np.unique(process)})
        entries.append(entry)
    return dict(interior_edges=edges.tolist(), bins=entries,
                nonfinite=completeness(y, ~finite),
                known=distribution(values[y >= 0]), unknown=distribution(values[y < 0]),
                common_vertex=distribution(values[y == 1]), accidental=distribution(values[y == 0]))


def reconstruct_labels(labels, purity=0.75, layers=3, common_cm=0.01, separate_cm=0.1):
    """Re-evaluate solely the label convention from isolated truth diagnostics."""
    if not (0 <= purity <= 1 and layers >= 1 and 0 <= common_cm < separate_cm):
        raise ValueError("Invalid sensitivity convention")
    purity_values = labels["match_purity"]
    eligible_legs = (np.isin(labels["match_state"], VALID_GEN_STATES) & np.isfinite(purity_values)
                     & (purity_values >= purity) & (labels["matched_layers"] >= layers))
    matched = eligible_legs.all(axis=1)
    duplicate = ((labels["same_gen_muon"] == 1) | (labels["same_sim_track"] == 1))
    distance = labels["vertex_distance_cm"]
    y = np.full(len(matched), -1, dtype=np.int8)
    reason = np.full(len(matched), "unknown_match", dtype="U32")
    y[matched & duplicate] = 0
    reason[matched & duplicate] = "same_muon_twice"
    candidate = matched & ~duplicate
    reason[candidate & ~np.isfinite(distance)] = "unknown_vertex"
    reason[candidate & np.isfinite(distance)] = "vertex_tolerance_gap"
    common = candidate & np.isfinite(distance) & (distance <= common_cm)
    separate = candidate & np.isfinite(distance) & (distance >= separate_cm)
    y[common], reason[common] = 1, "common_vertex"
    y[separate], reason[separate] = 0, "different_vertices"
    return y, reason


def detail_audit(labels, reco):
    if "match_state" not in labels:
        return dict(available=False, reason="Initial extraction stores only coarse reasons; no diagnostic truth is invented.")
    y = labels["y"]
    nominal_y, nominal_reason = reconstruct_labels(labels)
    if not np.array_equal(nominal_y, y) or not np.array_equal(nominal_reason, labels["reason"]):
        raise ValueError("Diagnostic reconstruction does not reproduce the nominal label convention")
    unknown = y < 0
    result = dict(available=True, nominal_labels_reproduced=True,
                  leg_match_states=dict(Counter(labels["match_state"].ravel().tolist())),
                  unknown_leg_match_states=dict(Counter(labels["match_state"][unknown].ravel().tolist())),
                  unknown_pair_states=dict(Counter("/".join(sorted(pair)) for pair in labels["match_state"][unknown])),
                  vertex_distance_cm={name: distribution(labels["vertex_distance_cm"][y == label])
                                      for label, name in LABEL_NAMES.items()},
                  common_distinct_immediate_mothers=int(((y == 1) & (labels["same_immediate_mother"] == 0)).sum()),
                  common_same_immediate_mother=int(((y == 1) & (labels["same_immediate_mother"] == 1)).sum()),
                  common_mother_identity_unavailable=int(((y == 1) & ~np.isin(labels["same_immediate_mother"], (0, 1))).sum()),
                  duplicate_gen_muon=int((labels["same_gen_muon"] == 1).sum()),
                  duplicate_sim_track=int((labels["same_sim_track"] == 1).sum()),
                  sensitivity=[],
                  sensitivity_warning="Recovered labels are convention changes, not evidence that matching becomes correct.")
    conventions = []
    for purity in (0.5, 0.75, 0.9):
        for layers in (2, 3, 4):
            conventions.append(dict(purity=purity, layers=layers, common_cm=0.01, separate_cm=0.1))
    for common in (0.001, 0.01, 0.05):
        for separate in (0.1, 0.3, 1.0):
            config = dict(purity=0.75, layers=3, common_cm=common, separate_cm=separate)
            if config not in conventions:
                conventions.append(config)
    for config in conventions:
        varied, reason = reconstruct_labels(labels, **config)
        changed = varied != y
        confusion = {LABEL_NAMES[a]: {LABEL_NAMES[b]: int(((y == a) & (varied == b)).sum())
                                     for b in (-1, 0, 1)} for a in (-1, 0, 1)}
        result["sensitivity"].append(dict(convention=config, counts=completeness(varied),
                    labels_changed=int(changed.sum()),
                    known_labels_flipped=int(((y >= 0) & (varied >= 0) & changed).sum()),
                    baseline_to_variant=confusion, reasons=dict(Counter(reason.tolist())),
                    per_process=grouped(reco["process"], varied)))
    result["per_process"] = {}
    for process in np.unique(reco["process"]):
        mask = reco["process"] == process
        result["per_process"][str(process)] = dict(
            leg_match_states=dict(Counter(labels["match_state"][mask].ravel().tolist())),
            unknown_leg_match_states=dict(Counter(labels["match_state"][mask & unknown].ravel().tolist())),
            common_distinct_immediate_mothers=int((mask & (y == 1) & (labels["same_immediate_mother"] == 0)).sum()))
    return result


def freshness(dataset_manifest, reco, reference, require_disjoint):
    if reference is None:
        if require_disjoint:
            raise ValueError("--require-disjoint needs --reference-dataset")
        return dict(checked=False)
    old_manifest, old_reco, _ = load_dataset(reference)
    old_paths = {str(Path(row["path"]).resolve()) for row in old_manifest["inputs"]}
    new_paths = {str(Path(row["path"]).resolve()) for row in dataset_manifest["inputs"]}
    def event_ids(data):
        return {(str(p), *map(int, identity)) for p, identity in zip(data["process"], data["ids"])}
    overlap_paths = old_paths & new_paths
    overlap_events = event_ids(old_reco) & event_ids(reco)
    if require_disjoint and (overlap_paths or overlap_events):
        raise ValueError(f"Fresh-sample protection failed: {len(overlap_paths)} shared sources, {len(overlap_events)} shared process/run/lumi/event identities")
    return dict(checked=True, reference_dataset=str(Path(reference).resolve()),
                reference_manifest_sha256=digest(Path(reference) / "inputs.json"),
                shared_source_files=len(overlap_paths), shared_reconstructed_events=len(overlap_events),
                disjoint_required=require_disjoint, disjoint=not (overlap_paths or overlap_events))


def build_report(dataset, manifest, reco, labels, reference=None, require_disjoint=False):
    y = labels["y"]
    names = list(reco["feature_names"])
    report = dict(dataset=str(Path(dataset).resolve()),
                  dataset_manifest_sha256=digest(Path(dataset) / "inputs.json"),
                  code_sha256=digest(__file__), sample_kind="simulation_only",
                  truth_usage="Isolated label construction and convention audit only. No training X or adversary targets change.",
                  classifier_fitted=False, deployment_ready=False, lifetime_independence_validated=False,
                  denominator="Already retained reconstructed dimuon pairs, including unknown labels.",
                  uncertainty="Wilson intervals treat pairs as independent approximations; repeated pairs from events can correlate.",
                  overall=completeness(y), reasons=dict(Counter(labels["reason"].tolist())),
                  per_process=grouped(reco["process"], y), per_stratum=grouped(reco["stratum"], y),
                  topology_pairs=category_pairs(reco, y, "topology"),
                  reconstruction_algorithm_pairs=category_pairs(reco, y, "recoAlgorithm"),
                  continuous={}, details=detail_audit(labels, reco),
                  fresh_sample=freshness(manifest, reco, reference, require_disjoint))
    for name in CONTINUOUS_FEATURES:
        report["continuous"][name] = numeric_audit(reco["X"][:, names.index(name)], y, reco["process"])
    for index, name in enumerate(("reconstructed_mass", "reconstructed_vertex_z", "reconstructed_vertex_radius")):
        report["continuous"][name] = numeric_audit(reco["nuisance"][:, index], y, reco["process"])
    return report


def draw_fraction(axis, entries, labels):
    x = np.arange(len(entries))
    for index, entry in enumerate(entries):
        estimate = entry["unknown_fraction"]
        if estimate is None:
            continue
        low, high = entry["unknown_ci95_pair_wilson"]
        axis.errorbar(index, 100 * estimate, yerr=[[100 * (estimate - low)], [100 * (high - estimate)]],
                      fmt="o", color="#2166ac", capsize=3)
        axis.text(index, 100 * high + 2, f'n={entry["n"]}', ha="center", fontsize=8)
    axis.set_xticks(x, labels, rotation=25, ha="right")
    axis.set_ylim(0, 110)
    axis.set_ylabel("Unknown labels (%)")
    axis.grid(axis="y", alpha=0.25)


def draw_quantiles(axis, item, title):
    labels = []
    for entry in item["bins"]:
        low, high = entry["lower"], entry["upper"]
        labels.append(f"< {high:.3g}" if low is None and high is not None else
                      f">= {low:.3g}" if high is None and low is not None else
                      f"{low:.3g}–{high:.3g}" if low is not None else "all finite")
    draw_fraction(axis, item["bins"], labels)
    axis.set_title(title)


def make_figure(report, output):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")
    processes = sorted(report["per_process"])
    draw_fraction(axes[0, 0], [report["per_process"][p] for p in processes], processes)
    axes[0, 0].set_title("Process")
    for axis, key, title in ((axes[0, 1], "topology_pairs", "Reconstructed topology pair (codes)"),
                              (axes[0, 2], "reconstruction_algorithm_pairs", "Reconstruction algorithm pair (codes)")):
        keys = sorted(report[key], key=lambda value: -report[key][value]["n"])
        draw_fraction(axis, [report[key][value] for value in keys], keys)
        axis.set_title(title)
    draw_quantiles(axes[1, 0], report["continuous"]["mu_min_nValidMuonHits"], "Smaller muon valid-hit count")
    draw_quantiles(axes[1, 1], report["continuous"]["mu_min_normalizedChi2"], "Smaller muon normalized χ²")
    draw_quantiles(axes[1, 2], report["continuous"]["reconstructed_vertex_z"], "Reconstructed pair vertex z (cm)")
    fig.suptitle(f'Label completeness: {report["overall"]["unknown"]}/{report["overall"]["n"]} unknown reconstructed pairs\n'
                 "Simulation only; reconstruction-quality bins use all retained pairs; intervals are pair-level approximations", fontsize=13)
    fig.savefig(output / "label_audit.pdf")
    fig.savefig(output / "label_audit.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--reference-dataset", help="Fixed earlier pilot for source/event overlap audit")
    parser.add_argument("--require-disjoint", action="store_true", help="Fail closed on reference source/event reuse")
    args = parser.parse_args()
    dataset, output = Path(args.dataset), Path(args.output)
    manifest, reco, labels = load_dataset(dataset)
    report = build_report(dataset, manifest, reco, labels, args.reference_dataset, args.require_disjoint)
    output.mkdir(parents=True, exist_ok=False)
    (output / "label_audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    make_figure(report, output)
    print(json.dumps(dict(output=str(output), overall=report["overall"], per_process=report["per_process"],
                          detailed_diagnostics=report["details"]["available"], fresh_sample=report["fresh_sample"])), flush=True)


if __name__ == "__main__":
    main()
