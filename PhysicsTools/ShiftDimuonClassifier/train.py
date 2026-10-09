#!/usr/bin/env python3
"""Bounded BDT / adversarial-NN comparison on reconstructed-only SM tables."""
import argparse
from collections import Counter
import hashlib
import json
import pickle
from pathlib import Path
import time
import numpy as np
import scipy
import sklearn
import torch
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
from adversarial import train_adversarial
from evaluation import balanced_weights, dependence, evaluate, event_split, weighted_quantile
from features import FEATURE_NAMES, QUALITY_NAMES, nuisance_transform
from inference import DistilledBDT, predict


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--strengths", type=float, nargs="+", default=[0., 0.1, 0.5, 2.])
    args = parser.parse_args()
    dataset, output = Path(args.dataset), Path(args.output)
    manifest = json.loads((dataset / "inputs.json").read_text())
    if manifest.get("sample_kind") != "simulation_only":
        raise ValueError("This pilot cannot open collision data")
    for filename, key in (("reco.npz", "reco_sha256"), ("labels.npz", "labels_sha256")):
        if hashlib.sha256((dataset / filename).read_bytes()).hexdigest() != manifest[key]:
            raise ValueError("Dataset digest changed")
    output.mkdir(parents=True, exist_ok=False)
    models_dir = output / "models"
    models_dir.mkdir()
    data = dict(np.load(dataset / "reco.npz", allow_pickle=False))
    if tuple(data["feature_names"]) != FEATURE_NAMES:
        raise ValueError("Stored feature column order differs from contract")
    labels = dict(np.load(dataset / "labels.npz", allow_pickle=False))
    y = labels["y"]
    split = event_split(data["ids"], data["process"], args.seed)
    np.savez_compressed(output / "splits.npz", split=split, ids=data["ids"], process=data["process"])
    labelled = y >= 0
    masks = [(split == i) & labelled for i in range(3)]
    tr, va, te = masks
    for mask in masks:
        if np.count_nonzero(mask & (y == 0)) < 10 or np.count_nonzero(mask & (y == 1)) < 30:
            raise ValueError("Insufficient labelled classes in event partition; preserve extraction")
    X = np.where(np.isfinite(data["X"]), data["X"], np.nan)
    z = nuisance_transform(data["nuisance"])
    wtr = balanced_weights(y[tr], data["process"][tr])
    wva = balanced_weights(y[va], data["process"][va])
    summary = dict(seed=args.seed, target_validation_common_efficiency=0.9,
                   runtime=dict(numpy=np.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__, torch=torch.__version__),
                   code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(Path(__file__).resolve().parent.glob("*.py"))},
                   features=list(FEATURE_NAMES), dataset_manifest_sha256=hashlib.sha256((dataset / "inputs.json").read_bytes()).hexdigest(),
                   weights="Class/process balanced training; raw and sampling-weighted diagnostics separately. No cross-section yield estimate.",
                   partitions={name: dict(rows=int(mask.sum()), labels=dict(Counter(y[mask].astype(str).tolist())))
                               for name, mask in zip(("train", "validation", "test"), masks)},
                   unknown_labels=int(np.count_nonzero(~labelled)), models={},
                   lifetime_independence_validated=False, deployment_ready=False)
    summary["per_process_labels"] = {
        str(p): {"all_pairs": int(np.count_nonzero(data["process"] == p)),
                 "common_vertex": int(np.count_nonzero((data["process"] == p) & (y == 1))),
                 "accidental": int(np.count_nonzero((data["process"] == p) & (y == 0))),
                 "unknown": int(np.count_nonzero((data["process"] == p) & (y < 0)))}
        for p in np.unique(data["process"])
    }
    summary["training_weight_effective_n"] = {}
    for label in (0, 1):
        for p in np.unique(data["process"][tr]):
            mask = (y[tr] == label) & (data["process"][tr] == p)
            if mask.any():
                w = wtr[mask]
                summary["training_weight_effective_n"][f"label{label}_{p}"] = float(w.sum() ** 2 / (w @ w))
    summary["feature_nuisance_spearman_train"] = {}
    for label in (0, 1):
        entries = {}
        for i, name in enumerate(FEATURE_NAMES):
            entry = {}
            for j, protected in enumerate(("mass", "vertex_z", "vertex_radius")):
                mask = tr & (y == label) & np.isfinite(X[:, i])
                value = np.nan
                if mask.sum() >= 5 and np.std(X[mask, i]) > 0 and np.std(z[mask, j]) > 0:
                    value = spearmanr(X[mask, i], z[mask, j]).statistic
                entry[protected] = float(value) if np.isfinite(value) else None
            entries[name] = entry
        summary["feature_nuisance_spearman_train"][str(label)] = entries
    candidate_models = {}
    def store_model(name, model, columns, details):
        bundle = dict(model=model, columns=np.asarray(columns), feature_names=FEATURE_NAMES,
                      requires_gen=False, deployment_ready=False)
        sv = predict(bundle, X[va])
        threshold = weighted_quantile(sv[y[va] == 1], wva[y[va] == 1], 0.1)
        bundle["threshold"] = threshold
        st = predict(bundle, X[te])
        result = evaluate(st, y[te], data["nuisance"][te], data["ids"][te],
                          data["process"][te], data["stratum"][te], data["sampling_weight"][te],
                          threshold, data["nuisance"][tr], y[tr])
        result.update(details)
        result["validation_auc"] = float(roc_auc_score(y[va], sv))
        result["validation_balanced_auc"] = float(roc_auc_score(y[va], sv, sample_weight=wva))
        result["validation_common_efficiency_raw"] = float(np.mean(sv[y[va] == 1] >= threshold))
        result["validation_common_efficiency_balanced"] = float(np.average(sv[y[va] == 1] >= threshold,
                                                                         weights=wva[y[va] == 1]))
        result["validation_dependence"] = dependence(sv, data["nuisance"][va], y[va])
        summary["models"][name] = result
        with (models_dir / (name + ".pkl")).open("wb") as stream:
            pickle.dump(bundle, stream)
        np.savez_compressed(output / (name + "_scores.npz"), score=predict(bundle, X), split=split)
        (output / "metrics.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        print(json.dumps(dict(model=name, auc=result["auc"], auc_ci95=result["auc_ci95"],
                              common_efficiency=result["common_vertex"]["efficiency"],
                              rejection=result["background_rejection"],
                              positive_joint_dcor=result["dependence"]["common_vertex"]["joint_dcor"])), flush=True)
        return bundle
    all_columns = np.arange(len(FEATURE_NAMES))
    started = time.monotonic()
    for name, columns in (("bdt_quality", np.arange(len(QUALITY_NAMES))), ("bdt_expanded", all_columns)):
        estimator = make_pipeline(SimpleImputer(strategy="median", keep_empty_features=True),
                                  GradientBoostingClassifier(n_estimators=120, max_depth=2, min_samples_leaf=12,
                                      learning_rate=0.05, subsample=0.8, random_state=args.seed))
        estimator.fit(X[tr][:, columns], y[tr], gradientboostingclassifier__sample_weight=wtr)
        importance = estimator[-1].feature_importances_
        store_model(name, estimator, columns, dict(kind="BDT", feature_importance={FEATURE_NAMES[i]: float(v)
                    for i, v in zip(columns, importance)}))
    for strength in args.strengths:
        name = "nn_lambda_" + str(strength).replace(".", "p")
        model, history = train_adversarial(X[tr], y[tr], wtr, z[tr], X[va], y[va], wva, z[va],
                                           seed=args.seed, strength=strength, epochs=args.epochs)
        (output / (name + "_training.json")).write_text(json.dumps(history, indent=2, allow_nan=False) + "\n")
        bundle = store_model(name, model, all_columns, dict(kind="adversarial NN" if strength else "unconstrained NN",
                                                          adversarial_strength=strength))
        if strength > 0:
            candidate_models[name] = bundle
    if candidate_models:
        def objective(name):
            value = summary["models"][name]
            deps = value["validation_dependence"]
            penalty = max(deps[k]["joint_dcor"] or 0 for k in ("common_vertex", "accidental"))
            return value["validation_auc"] - 0.5 * penalty
        teacher_name = max(candidate_models, key=objective)
        teacher = candidate_models[teacher_name]
        regression = make_pipeline(SimpleImputer(strategy="median", keep_empty_features=True),
                                   GradientBoostingRegressor(n_estimators=150, max_depth=2, min_samples_leaf=12,
                                                           learning_rate=0.05, random_state=args.seed))
        regression.fit(X[tr], predict(teacher, X[tr]), gradientboostingregressor__sample_weight=wtr)
        store_model("bdt_distilled", DistilledBDT(regression), all_columns,
                    dict(kind="BDT distillation", teacher=teacher_name,
                         teacher_selection="validation AUC - 0.5*max(class conditional joint dCor); test excluded"))
    summary["training_seconds"] = time.monotonic() - started
    (output / "metrics.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(dict(completed=str(output), seconds=summary["training_seconds"])), flush=True)


if __name__ == "__main__":
    main()
