#!/usr/bin/env python3
"""Export a trusted, locked study UniformityBDT pickle to explicit JSON.

Run in the pinned training environment. Pickle/sklearn dependencies are loaded
only during export/validation, never by portable_inference in production.
"""
import argparse
import json
from pathlib import Path
import pickle
import sys

import numpy as np

from features import FEATURE_NAMES
from portable_inference import (INPUT_CAST, SCHEMA, SCORE_FORMULA, feature_contract_digest,
                                file_digest, predict, validate_model)


def export_bundle(bundle, provenance):
    """Serialize ordered nodes, class votes, coefficients and fixed preprocessing."""
    from uniform_bdt import UniformityBDT
    model = bundle["model"]
    if not isinstance(model, UniformityBDT) or tuple(bundle["feature_names"]) != FEATURE_NAMES:
        raise ValueError("Only the declared UniformityBDT feature contract can be exported")
    if bundle.get("requires_gen") is not False:
        raise ValueError("Cannot export a model requiring truth")
    columns = np.asarray(bundle["columns"])
    if columns.ndim != 1 or columns.dtype.kind not in "iu" or len(columns) != model.n_features_in_:
        raise ValueError("Invalid stored model feature columns")
    imputer = model.imputer
    if imputer.strategy != "median" or imputer.add_indicator or not imputer.keep_empty_features:
        raise ValueError("Unsupported training imputer")
    if (list(model.classes_) != [0, 1] or len(model.estimators_) != len(model.estimator_weights_)
            or len(model.imputer.statistics_) != model.n_features_in_):
        raise ValueError("Inconsistent stored UniformityBDT model")
    trees = []
    for estimator, coefficient in zip(model.estimators_, model.estimator_weights_):
        if list(estimator.classes_) != [0, 1] or estimator.n_features_in_ != model.n_features_in_:
            raise ValueError("Stored tree class/feature contract differs")
        tree = estimator.tree_
        classes = estimator.classes_[np.argmax(tree.value[:, 0, :], axis=1)]
        leaves = tree.children_left == -1
        trees.append(dict(coefficient=float(coefficient),
                          children_left=tree.children_left.tolist(), children_right=tree.children_right.tolist(),
                          feature=tree.feature.tolist(), threshold=tree.threshold.tolist(),
                          leaf_class=[int(classes[i]) if leaves[i] else None for i in range(tree.node_count)]))
    result = dict(schema=SCHEMA, feature_names=list(FEATURE_NAMES),
                  feature_contract_sha256=feature_contract_digest(FEATURE_NAMES),
                  columns=columns.astype(int).tolist(),
                  imputer_statistics=model.imputer.statistics_.tolist(),
                  imputer_nonfinite_policy="replace-with-frozen-training-median",
                  input_cast=INPUT_CAST, score_formula=SCORE_FORMULA, classes=[0, 1],
                  sum_abs_coefficients=float(np.sum(np.abs(model.estimator_weights_))), trees=trees,
                  threshold_metadata=float(bundle["threshold"]),
                  threshold_interpretation="Locked historical SM validation cut, stored as metadata; no cut is applied",
                  score_interpretation="Uncalibrated provisional common-vertex discriminator",
                  score_storage="Float64", requires_gen=False, selection_applied=False,
                  physics_ready=False, mass_lifetime_independence_validated=False,
                  provenance=provenance)
    return validate_model(result)


def validate_scores(model, bundle, X):
    """Demand exact original margins/scores and threshold decisions on a fixture."""
    from inference import predict as original_predict
    expected, actual = original_predict(bundle, X), predict(model, X)
    np.testing.assert_array_equal(actual, expected, err_msg="Portable BDT scores differ from the frozen pickle")
    cut = model["threshold_metadata"]
    np.testing.assert_array_equal(actual >= cut, expected >= cut)
    return dict(rows=len(actual), bitwise_score_equality=True, threshold_decision_equality=True,
                maximum_absolute_score_difference=float(np.max(np.abs(actual - expected))) if len(actual) else 0.,
                score_dtype=str(actual.dtype), selection_applied=False)


def export_model(results, model_name, output, validation_table=None):
    """Freeze a JSON model and an adjacent receipt; never overwrite either."""
    results, output = Path(results), Path(output)
    receipt_path = output.with_suffix(output.suffix + ".receipt.json")
    if output.exists() or receipt_path.exists():
        raise FileExistsError("Refuse to overwrite an exported BDT or its receipt")
    if Path(model_name).name != model_name:
        raise ValueError("Invalid study model name")
    lock_path = results / "selection_locked.json"
    lock = json.loads(lock_path.read_text())
    if lock.get("selection_uses_test") is not False or model_name not in lock.get("representative_models", []):
        raise ValueError("Export requires a locked validation-selected representative")
    pickle_path = results / "models" / (model_name + ".pkl")
    # Only trusted local study files are accepted. Unpickling imports the pinned
    # original implementation here; portable production scoring loads JSON only.
    with pickle_path.open("rb") as source:
        bundle = pickle.load(source)
    if float(bundle["threshold"]) != float(lock["thresholds"][model_name]):
        raise ValueError("Pickle threshold differs from its lock receipt")
    here = Path(__file__).resolve().parent
    provenance = dict(model_name=model_name, pickle_sha256=file_digest(pickle_path),
                      selection_lock_sha256=file_digest(lock_path),
                      source_sha256={name: file_digest(here / name) for name in
                                     ("features.py", "inference.py", "uniform_bdt.py", "portable_inference.py", "export_bdt.py")},
                      export_runtime=dict(python=sys.version, numpy=np.__version__))
    model = export_bundle(bundle, provenance)
    validation = None
    if validation_table is not None:
        validation_table = Path(validation_table)
        with np.load(validation_table, allow_pickle=False) as data:
            if tuple(data["feature_names"]) != FEATURE_NAMES:
                raise ValueError("Export validation table feature contract differs")
            validation = validate_scores(model, bundle, data["X"])
        if not validation["rows"]:
            raise ValueError("Export validation table must exercise at least one feature row")
        validation.update(table=str(validation_table.resolve()), table_sha256=file_digest(validation_table))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as destination:
        destination.write(json.dumps(model, indent=2, allow_nan=False) + "\n")
    receipt = dict(schema="shift-dimuon-bdt-export-receipt-v1", model=str(output.resolve()),
                   model_sha256=file_digest(output), model_name=model_name,
                   provenance=provenance, validation=validation, selection_applied=False,
                   physics_ready=False, production_scoring_parity_validated=bool(validation))
    with receipt_path.open("x") as destination:
        destination.write(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    return model, receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--model", default="uniform_bdt_3p0_s71")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--validation-table", type=Path,
                        help="Frozen reco.npz fixture; require bitwise equality with the saved sklearn model")
    args = parser.parse_args()
    _, receipt = export_model(args.results, args.model, args.output, args.validation_table)
    print(json.dumps(dict(output=str(args.output), model=receipt["model_name"],
                          parity_validated=receipt["production_scoring_parity_validated"],
                          validation_rows=(receipt["validation"] or {}).get("rows", 0),
                          selection_applied=False, physics_ready=False)), flush=True)


if __name__ == "__main__":
    main()
