"""Explicit UniformityBDT JSON scoring with NumPy and the reco feature contract.

No pickle, sklearn, scipy, labels, nuisance targets or selection are used here.
The score is an uncalibrated study discriminator; save it as a Float64 value.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from features import FEATURE_NAMES


SCHEMA = "shift-dimuon-uniform-bdt-json-v1"
INPUT_CAST = "float32-after-float64-imputation"
SCORE_FORMULA = "1/(1+exp(-2*ordered_signed_tree_margin/sum_abs_coefficients))"


def file_digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def feature_contract_digest(names):
    return hashlib.sha256(json.dumps(list(names), separators=(",", ":")).encode()).hexdigest()


def _integer(value):
    return type(value) is int


def _finite_number(value):
    return type(value) in (int, float) and math.isfinite(value)


def _sha(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def validate_model(model):
    """Fail closed on source/feature mismatches and malformed tree topology."""
    if not isinstance(model, dict) or model.get("schema") != SCHEMA:
        raise ValueError("Unsupported portable BDT schema")
    if tuple(model.get("feature_names", ())) != FEATURE_NAMES:
        raise ValueError("Portable BDT feature names/order differ")
    if model.get("feature_contract_sha256") != feature_contract_digest(FEATURE_NAMES):
        raise ValueError("Portable BDT feature contract digest differs")
    if model.get("input_cast") != INPUT_CAST or model.get("score_formula") != SCORE_FORMULA:
        raise ValueError("Unsupported portable BDT numeric contract")
    classes = model.get("classes")
    if (classes != [0, 1] or not all(_integer(c) for c in classes)
            or model.get("requires_gen") is not False):
        raise ValueError("Portable BDT must have binary classes and reconstructed-only inputs")
    if model.get("selection_applied") is not False or model.get("physics_ready") is not False:
        raise ValueError("This model stores a provisional score without a selection")
    columns = model.get("columns")
    if (not isinstance(columns, list) or not columns or not all(_integer(i) for i in columns)
            or len(set(columns)) != len(columns) or any(i < 0 or i >= len(FEATURE_NAMES) for i in columns)):
        raise ValueError("Invalid selected feature columns")
    statistics = model.get("imputer_statistics")
    if (not isinstance(statistics, list) or len(statistics) != len(columns)
            or not all(_finite_number(x) for x in statistics)):
        raise ValueError("Invalid frozen imputer statistics")
    if model.get("imputer_nonfinite_policy") != "replace-with-frozen-training-median":
        raise ValueError("Unsupported imputer policy")
    threshold = model.get("threshold_metadata")
    if not _finite_number(threshold) or not 0 <= threshold <= 1:
        raise ValueError("Invalid locked threshold metadata")
    provenance = model.get("provenance", {})
    if not isinstance(provenance, dict):
        raise ValueError("Invalid model provenance")
    for name in ("pickle_sha256", "selection_lock_sha256"):
        if not _sha(provenance.get(name)):
            raise ValueError("Missing model/lock provenance digest")
    sources = provenance.get("source_sha256", {})
    if not isinstance(sources, dict):
        raise ValueError("Invalid source provenance")
    for name in ("features.py", "inference.py", "uniform_bdt.py", "portable_inference.py", "export_bdt.py"):
        if not _sha(sources.get(name)):
            raise ValueError("Missing source provenance digest: " + name)
    here = Path(__file__).resolve().parent
    # Historical pickle implementation hashes are provenance only. These two
    # files define the currently executed collision-data scoring contract.
    for name in ("features.py", "portable_inference.py"):
        if file_digest(here / name) != sources[name]:
            raise ValueError("Portable BDT scoring source differs: " + name)
    trees = model.get("trees")
    if not isinstance(trees, list):
        raise ValueError("Portable BDT trees must be an ordered list")
    coefficients = []
    for tree in trees:
        if not isinstance(tree, dict) or not _finite_number(tree.get("coefficient")) or tree["coefficient"] <= 0:
            raise ValueError("Invalid tree coefficient")
        coefficients.append(tree["coefficient"])
        arrays = [tree.get(key) for key in ("children_left", "children_right", "feature", "threshold", "leaf_class")]
        if not all(isinstance(a, list) for a in arrays) or not arrays[0] or len({len(a) for a in arrays}) != 1:
            raise ValueError("Tree node arrays differ")
        left, right, features, thresholds, classes = arrays
        n = len(left)
        for i in range(n):
            if not all(_integer(x) for x in (left[i], right[i], features[i])) or not _finite_number(thresholds[i]):
                raise ValueError("Invalid tree node index/threshold")
            if left[i] == right[i] == -1:
                if features[i] != -2 or type(classes[i]) is not int or classes[i] not in (0, 1):
                    raise ValueError("Invalid tree leaf class")
            elif (not 0 <= left[i] < n or not 0 <= right[i] < n or left[i] == right[i]
                  or not 0 <= features[i] < len(columns) or classes[i] is not None):
                raise ValueError("Invalid tree child/feature references")
        seen, pending = set(), [0]
        while pending:
            node = pending.pop()
            if node in seen:
                raise ValueError("Tree has a cycle or shared child")
            seen.add(node)
            if left[node] != -1:
                pending.extend((left[node], right[node]))
        if len(seen) != n:
            raise ValueError("Tree has unreachable nodes")
    normalizer = model.get("sum_abs_coefficients")
    if not _finite_number(normalizer) or normalizer < 0:
        raise ValueError("Invalid ensemble normalizer")
    expected = float(np.sum(np.abs(np.asarray(coefficients, dtype=np.float64))))
    if normalizer != expected or (bool(trees) != bool(normalizer)):
        raise ValueError("Ensemble normalizer differs from ordered coefficients")
    return model


def load_model(path):
    """Load an explicit JSON model and verify its numeric/source contract."""
    return validate_model(json.loads(Path(path).read_text()))


def predict(model, X):
    """Return one Float64 score per reconstructed feature row, without a cut."""
    validate_model(model)
    raw = np.asarray(X, dtype=np.float64)
    if raw.ndim != 2 or raw.shape[1] != len(FEATURE_NAMES):
        raise ValueError("Wrong reconstructed feature matrix shape")
    if not len(raw):
        return np.empty(0, dtype=np.float64)
    selected = raw[:, model["columns"]]
    imputed = np.where(np.isfinite(selected), selected, np.asarray(model["imputer_statistics"], dtype=np.float64))
    # sklearn trees cast to Float32 after Float64 imputation. Promote those
    # rounded values to Float64 to compare against Float64 split thresholds;
    # a Float32 array/Python-scalar comparison can round the threshold too.
    with np.errstate(over="ignore", invalid="ignore"):
        features = imputed.astype(np.float32).astype(np.float64)
    if not np.isfinite(features).all():
        raise ValueError("Features overflow sklearn-compatible Float32 tree inputs")
    margin = np.zeros(len(features), dtype=np.float64)
    for tree in model["trees"]:
        left = np.asarray(tree["children_left"], dtype=np.int64)
        right = np.asarray(tree["children_right"], dtype=np.int64)
        column = np.asarray(tree["feature"], dtype=np.int64)
        threshold = np.asarray(tree["threshold"], dtype=np.float64)
        nodes = np.zeros(len(features), dtype=np.int64)
        while True:
            active = np.flatnonzero(left[nodes] != -1)
            if not len(active):
                break
            current = nodes[active]
            take_left = features[active, column[current]] <= threshold[current]
            nodes[active] = np.where(take_left, left[current], right[current])
        classes = np.asarray([0 if c is None else c for c in tree["leaf_class"]], dtype=np.int64)
        margin += tree["coefficient"] * (2.0 * classes[nodes] - 1.0)
    normalizer = model["sum_abs_coefficients"]
    if normalizer > 0:
        margin /= normalizer
    # Preserve the original ordered arithmetic and C double expit equation.
    return np.asarray([1.0 / (1.0 + math.exp(-2.0 * value)) for value in margin], dtype=np.float64)
