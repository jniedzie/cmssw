#!/usr/bin/env python3
"""Plot saved predictions from the bounded, simulation-only classifier study."""

import argparse
import hashlib
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
import numpy as np
from sklearn.metrics import roc_curve

from evaluation import wilson


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _selected_models(metrics, requested=None):
    """Keep legacy comparisons complete; honor an audited scan summary."""
    selected = (requested if requested is not None else metrics["representative_models"]
                if "representative_models" in metrics else list(metrics.get("models", {})))
    if not isinstance(selected, (list, tuple)) or not selected:
        raise ValueError("Model selection must be a nonempty list")
    if not all(isinstance(name, str) and Path(name).name == name for name in selected):
        raise ValueError("Invalid model filename in selection")
    if len(set(selected)) != len(selected):
        raise ValueError("Duplicate models in selection")
    missing = set(selected) - set(metrics.get("models", {}))
    if missing:
        raise ValueError("Selected models absent from metrics: " + ", ".join(sorted(missing)))
    return list(selected)


def load_study(dataset, results, models=None):
    """Require the saved MC inventory and exact event/score alignment."""
    manifest_path = dataset / "inputs.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("sample_kind") != "simulation_only":
        raise ValueError("Plots are restricted to explicitly identified simulation")
    for filename, key in (("reco.npz", "reco_sha256"), ("labels.npz", "labels_sha256")):
        if _digest(dataset / filename) != manifest.get(key):
            raise ValueError("Dataset digest changed: " + filename)
    metrics = json.loads((results / "metrics.json").read_text())
    if _digest(manifest_path) != metrics.get("dataset_manifest_sha256"):
        raise ValueError("Results were not produced from this dataset manifest")
    with np.load(dataset / "reco.npz", allow_pickle=False) as archive:
        reco = {key: archive[key] for key in ("ids", "process", "nuisance")}
    with np.load(dataset / "labels.npz", allow_pickle=False) as archive:
        y = archive["y"]
    with np.load(results / "splits.npz", allow_pickle=False) as archive:
        split = archive["split"]
        if not np.array_equal(archive["ids"], reco["ids"]) or not np.array_equal(archive["process"], reco["process"]):
            raise ValueError("Saved split identities differ from the dataset")
    n = len(y)
    if y.shape != (n,) or split.shape != (n,) or reco["nuisance"].shape != (n, 3):
        raise ValueError("Saved label, split, or nuisance shape differs")
    if not np.isin(split, (0, 1, 2)).all() or not np.isin(y, (-1, 0, 1)).all():
        raise ValueError("Unexpected saved label or partition")
    scores = {}
    for name in _selected_models(metrics, models):
        with np.load(results / (name + "_scores.npz"), allow_pickle=False) as archive:
            score = archive["score"]
            if score.shape != (n,) or not np.isfinite(score).all():
                raise ValueError("Invalid score array for " + name)
            if not np.array_equal(archive["split"], split):
                raise ValueError("Score partition differs for " + name)
            scores[name] = score
    if not scores:
        raise ValueError("No trained models in the saved results")
    training = (split == 0) & (y == 1)
    test = (split == 2) & (y >= 0)
    if not training.any() or len(np.unique(y[test])) != 2:
        raise ValueError("Plots require genuine training pairs and both held-out classes")
    if not np.isfinite(reco["nuisance"]).all():
        raise ValueError("Nonfinite reconstructed nuisance in saved dataset")
    return reco, y, split, metrics, scores


def _name(name):
    seed = None
    prefix = re.fullmatch(r"s(\d+)_(.+)", name)
    if prefix:
        seed, name = prefix.groups()
    suffix = re.fullmatch(r"(.+)_s(\d+)", name)
    if suffix:
        name, seed = suffix.groups()
    labels = {"bdt_quality": "BDT quality", "bdt_expanded": "BDT expanded",
              "bdt_distilled": "BDT distilled"}
    if name in labels:
        label = labels[name]
    elif name.startswith("nn_lambda_"):
        label = "NN lambda=" + name[len("nn_lambda_"):].replace("p", ".")
    elif name.startswith("nn_disco_lambda_"):
        label = "NN DisCo lambda=" + name[len("nn_disco_lambda_"):].replace("p", ".")
    elif name.startswith("nn_disco_"):
        label = "NN DisCo lambda=" + name[len("nn_disco_"):].replace("p", ".")
    elif name.startswith("nn_adversarial_"):
        label = "Adversarial NN lambda=" + name[len("nn_adversarial_"):].replace("p", ".")
    elif name.startswith("uniform_bdt_"):
        label = "Uniform BDT strength=" + name[len("uniform_bdt_"):].replace("p", ".")
    elif name.startswith("bdt_"):
        label = "BDT " + re.sub(r"(?<=\d)p(?=\d)", ".", name[len("bdt_"):].replace("_", " "))
    else:
        label = name
    return label + (f" (seed {seed})" if seed is not None else "")


def _bound(value):
    return f"{value:.5g}"


def _bin_labels(edges, unit=""):
    return [f"< {_bound(edges[0])}{unit}",
            f"{_bound(edges[0])} to {_bound(edges[1])}{unit}",
            f">= {_bound(edges[1])}{unit}"]


def _efficiency(mask, passed):
    n = int(mask.sum())
    k = int((mask & passed).sum())
    return n, k, k / n if n else None, wilson(k, n)


def _errorbar(axis, x, estimate, interval, color, label=None):
    if estimate is not None:
        axis.errorbar(x, estimate,
                      yerr=np.array([[estimate - interval[0]], [interval[1] - estimate]]),
                      fmt="o", markersize=4, capsize=2, color=color, label=label)


def plot_overview(reco, y, split, metrics, scores, destination):
    test = (split == 2) & (y >= 0)
    genuine_train = (split == 0) & (y == 1)
    mass_edges = np.quantile(reco["nuisance"][genuine_train, 0], [1 / 3, 2 / 3])
    mass_bins = np.searchsorted(mass_edges, reco["nuisance"][:, 0], side="right")
    names = list(scores)
    colors = plt.get_cmap("tab10")(np.arange(len(names)) % 10)
    offsets = np.linspace(-0.2, 0.2, len(names)) if len(names) > 1 else [0.]
    figure, axes = plt.subplots(2, 2, figsize=(14, 10))
    roc_axis, mass_axis, process_axis, dependence_axis = axes.ravel()
    processes = [p for p in ("jpsi", "dy", "qcd") if ((reco["process"] == p) & test).any()]
    process_labels = {"jpsi": "J/psi", "dy": "DY", "qcd": "QCD"}
    for model_index, (name, color, offset) in enumerate(zip(names, colors, offsets)):
        score = scores[name]
        entry = metrics["models"][name]
        threshold = float(entry["threshold"])
        passed = score >= threshold
        fpr, tpr, _ = roc_curve(y[test], score[test])
        roc_axis.plot(fpr, tpr, color=color,
                      label=f"{_name(name)} (AUC {entry['auc']:.3f})")
        efficiencies = []
        for b in range(3):
            _, _, efficiency, interval = _efficiency(test & (y == 1) & (mass_bins == b), passed)
            _errorbar(mass_axis, b + offset, efficiency, interval, color)
            efficiencies.append(efficiency if efficiency is not None else np.nan)
        mass_axis.plot(np.arange(3) + offset, efficiencies, color=color, alpha=0.45, linewidth=0.8)
        for p_index, process in enumerate(processes):
            _, _, efficiency, interval = _efficiency(test & (y == 1) & (reco["process"] == process), passed)
            _errorbar(process_axis, p_index + offset, efficiency, interval, color,
                      label=_name(name) if p_index == 0 else None)
        for label_index, label in enumerate(("common_vertex", "accidental")):
            value = entry["dependence"][label].get("joint_dcor")
            if value is not None:
                dependence_axis.bar(model_index + (-0.18 if label_index == 0 else 0.18), value,
                                    width=0.34, color=color if label_index == 0 else "#b5b5b5",
                                    hatch=None if label_index == 0 else "//")
    roc_axis.plot([0, 1], [0, 1], "k:", linewidth=1)
    roc_axis.set(xlabel="Accidental-pair acceptance", ylabel="Genuine-pair efficiency",
                 title="Held-out ROC (unweighted)", xlim=(0, 1), ylim=(0, 1.02))
    roc_axis.legend(fontsize=8, loc="lower right")
    mass_axis.set(xticks=np.arange(3), xticklabels=_bin_labels(mass_edges, " GeV"),
                  xlabel="Reconstructed dimuon mass: training-derived thirds",
                  ylabel="Genuine-pair efficiency", title="Mass acceptance at saved thresholds",
                  ylim=(0, 1.05))
    for label in mass_axis.get_xticklabels():
        label.set_fontsize(8)
    process_axis.set(xticks=np.arange(len(processes)), xticklabels=[
        f"{process_labels[p]}\n(n={int((test & (y == 1) & (reco['process'] == p)).sum())})"
        for p in processes], ylabel="Genuine-pair efficiency", ylim=(0, 1.05),
        title="Preservation of authentic SM pairs (95% intervals)")
    process_axis.legend(fontsize=8, loc="lower right")
    dependence_axis.set(xticks=np.arange(len(names)), xticklabels=[_name(n) for n in names],
                       ylabel="Empirical joint distance correlation", ylim=(0, 1.05),
                       title="Score dependence on reconstructed mass, z, radius")
    dependence_axis.tick_params(axis="x", labelrotation=25, labelsize=8)
    dependence_axis.legend(handles=[Rectangle((0, 0), 1, 1, facecolor=colors[0], label="Genuine pairs"),
                                   Rectangle((0, 0), 1, 1, facecolor="#b5b5b5", hatch="//",
                                             label="Accidental pairs")], fontsize=8)
    for axis in axes.ravel():
        axis.grid(alpha=0.2)
        axis.set_axisbelow(True)
    figure.suptitle("Small held-out SM study: reliably matched pairs only\n"
                   f"{int((test & (y == 1)).sum())} genuine / {int((test & (y == 0)).sum())} accidental pairs; "
                   "mass and lifetime independence unvalidated", fontsize=13)
    figure.text(0.5, 0.025,
                "Thresholds fixed on validation. Efficiencies and ROC are unweighted. "
                "Wilson intervals approximate pair independence.\n"
                "Joint dCor uses the saved class-conditional diagnostic; finite-sample values have a nonzero bias. "
                "Reconstructed coordinates are displacement proxies.",
                ha="center", fontsize=8)
    figure.tight_layout(rect=(0, 0.075, 1, 0.94))
    for extension in ("pdf", "png"):
        figure.savefig(destination / ("overview." + extension), dpi=180)
    plt.close(figure)


def _joint_panel(axis, y, test, passed, mass_bin, z_bin, process_mask, mass_edges, z_edges, title):
    efficiency_map = np.full((3, 3), np.nan)
    annotations = []
    for m in range(3):
        for z in range(3):
            mask = test & (y == 1) & process_mask & (mass_bin == m) & (z_bin == z)
            n, k, efficiency, interval = _efficiency(mask, passed)
            if n >= 20:
                efficiency_map[z, m] = efficiency
            if n:
                annotation = f"{k}/{n}\n{100 * efficiency:.0f}% [{100 * interval[0]:.0f}, {100 * interval[1]:.0f}]"
                if n < 20:
                    annotation += "\nlow support"
            else:
                annotation = "0/0\nunvalidated"
            annotations.append((m, z, annotation, n >= 20 and efficiency > 0.55))
    color_map = plt.get_cmap("Blues").copy()
    color_map.set_bad("#dddddd")
    image = axis.imshow(np.ma.masked_invalid(efficiency_map), origin="lower", vmin=0, vmax=1,
                        cmap=color_map, interpolation="nearest", aspect="auto")
    for m, z, annotation, white in annotations:
        axis.text(m, z, annotation, ha="center", va="center", fontsize=8,
                  color="white" if white else "#222222")
    axis.set(xticks=np.arange(3), xticklabels=_bin_labels(mass_edges),
             yticks=np.arange(3), yticklabels=_bin_labels(z_edges),
             xlabel="Reconstructed dimuon mass [GeV]", ylabel="Reconstructed vertex z [cm]", title=title)
    axis.tick_params(labelsize=8)
    return image


def plot_joint(reco, y, split, metrics, scores, destination):
    training = (split == 0) & (y == 1)
    test = (split == 2) & (y >= 0)
    # Exactly three quantile intervals; tied boundaries honestly leave empty cells.
    mass_edges = np.quantile(reco["nuisance"][training, 0], [1 / 3, 2 / 3])
    z_edges = np.quantile(reco["nuisance"][training, 1], [1 / 3, 2 / 3])
    mass_bin = np.searchsorted(mass_edges, reco["nuisance"][:, 0], side="right")
    z_bin = np.searchsorted(z_edges, reco["nuisance"][:, 1], side="right")
    panels = [(None, "All SM processes"), ("jpsi", "J/psi"), ("dy", "DY"), ("qcd", "QCD")]
    with PdfPages(destination / "joint_efficiency.pdf") as pdf:
        for name, score in scores.items():
            threshold = float(metrics["models"][name]["threshold"])
            passed = score >= threshold
            figure, axes = plt.subplots(2, 2, figsize=(12, 10))
            for axis, (process, title) in zip(axes.ravel(), panels):
                process_mask = np.ones(len(y), dtype=bool) if process is None else reco["process"] == process
                image = _joint_panel(axis, y, test, passed, mass_bin, z_bin, process_mask,
                                     mass_edges, z_edges, title)
                figure.colorbar(image, ax=axis, fraction=0.045, pad=0.03, label="Genuine-pair efficiency")
            figure.suptitle(f"{_name(name)}: held-out matched-only SM efficiency\n"
                           f"Saved threshold {threshold:.5g}; all reconstructed radii combined", fontsize=13)
            figure.text(0.5, 0.025,
                        "Cells show passed/total and efficiency [95% Wilson interval]. "
                        "Grey: fewer than 20 genuine test pairs; unvalidated.\n"
                        "Bins use genuine training pairs only. Intervals approximate pair independence. "
                        "Radius is integrated out; no joint 3D or lifetime closure is claimed.",
                        ha="center", fontsize=8)
            figure.tight_layout(rect=(0, 0.07, 1, 0.94))
            pdf.savefig(figure)
            plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--models", nargs="+", help="Explicit model list; overrides representative_models in metrics")
    args = parser.parse_args()
    reco, y, split, metrics, scores = load_study(args.dataset, args.results, args.models)
    destination = args.results / "figures"
    destination.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42})
    plot_overview(reco, y, split, metrics, scores, destination)
    plot_joint(reco, y, split, metrics, scores, destination)
    print(json.dumps({"sample_kind": "simulation_only", "held_out_labelled_pairs": int(((split == 2) & (y >= 0)).sum()),
                      "models": list(scores), "figures": [str(destination / filename)
                         for filename in ("overview.pdf", "overview.png", "joint_efficiency.pdf")],
                      "lifetime_independence_validated": False}), flush=True)


if __name__ == "__main__":
    main()
