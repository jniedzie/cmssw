"""Held-out, event-grouped diagnostics; no assertion of lifetime closure."""
import hashlib
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from features import nuisance_transform


def event_split(ids, process, seed=42):
    result = []
    for identity, sample in zip(ids, process):
        key = f"{seed}:{sample}:" + ":".join(str(int(x)) for x in identity)
        value = int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "big") / 2 ** 64
        result.append(0 if value < 0.6 else 1 if value < 0.8 else 2)
    return np.asarray(result, dtype=np.int8)


def balanced_weights(y, process):
    """Equal label totals; equal available processes within each label."""
    weights = np.zeros(len(y))
    for label in (0, 1):
        samples = np.unique(process[y == label])
        for sample in samples:
            mask = (y == label) & (process == sample)
            weights[mask] = 1 / (2 * len(samples) * np.count_nonzero(mask))
    return weights * len(y)


def weighted_quantile(values, weights, probability):
    order = np.argsort(values, kind="stable")
    cumulative = np.cumsum(weights[order])
    return float(values[order[np.searchsorted(cumulative, probability * cumulative[-1])]])


def wilson(k, n):
    if n == 0:
        return [None, None]
    z = 1.96
    p = k / n
    center = (p + z * z / (2 * n)) / (1 + z * z / n)
    delta = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [float(max(0, center - delta)), float(min(1, center + delta))]


def distance_correlation(score, nuisance, maximum=400):
    """Biased empirical dCor, so small nonzero values do not prove dependence."""
    score, nuisance = np.asarray(score), np.asarray(nuisance)
    if len(score) < 5:
        return None
    if len(score) > maximum:
        pick = np.random.default_rng(917).choice(len(score), maximum, replace=False)
        score, nuisance = score[pick], nuisance[pick]
    nuisance = np.atleast_2d(nuisance).reshape(len(score), -1)
    scale = np.std(nuisance, axis=0)
    nuisance = nuisance / np.where(scale > 0, scale, 1)
    a = np.abs(score[:, None] - score[None, :])
    b = np.linalg.norm(nuisance[:, None, :] - nuisance[None, :, :], axis=2)
    a = a - a.mean(0)[None, :] - a.mean(1)[:, None] + a.mean()
    b = b - b.mean(0)[None, :] - b.mean(1)[:, None] + b.mean()
    denominator = np.sqrt(np.mean(a * a) * np.mean(b * b))
    return float(np.sqrt(max(0., np.mean(a * b) / denominator))) if denominator > 0 else 0.


def dependence(score, nuisance, labels):
    z = nuisance_transform(nuisance)
    result = {}
    for label, name in ((0, "accidental"), (1, "common_vertex")):
        mask = labels == label
        entry = {"n": int(mask.sum()), "joint_dcor": distance_correlation(score[mask], z[mask])}
        for i, variable in enumerate(("mass", "vertex_z", "vertex_radius")):
            entry[variable + "_dcor"] = distance_correlation(score[mask], z[mask, i:i+1])
            value = spearmanr(score[mask], z[mask, i]).statistic if mask.sum() >= 5 else np.nan
            entry[variable + "_spearman"] = float(value) if np.isfinite(value) else None
        result[name] = entry
    return result


def fraction(mask, passed, weights=None):
    n = int(np.count_nonzero(mask))
    k = int(np.count_nonzero(mask & passed))
    result = dict(n=n, passed=k, efficiency=k / n if n else None, ci95=wilson(k, n))
    if weights is not None and n:
        w = weights[mask]
        result.update(weight_sum=float(w.sum()), effective_n=float(w.sum() ** 2 / (w @ w)) if w @ w > 0 else None,
                      sampling_weighted_efficiency=float(weights[mask & passed].sum() / w.sum()) if w.sum() else None)
    return result


def binned_efficiencies(nuisance, y, passed, training_nuisance, training_y, minimum=20):
    """Training-derived bins; explicitly flag unvalidated low-support cells."""
    result = {}
    train = training_nuisance[training_y == 1]
    for i, name in enumerate(("mass", "vertex_z", "vertex_radius")):
        edges = np.unique(np.quantile(train[:, i], [1/3, 2/3]))
        bins = np.searchsorted(edges, nuisance[:, i], side="right")
        entries = []
        for b in range(len(edges) + 1):
            item = fraction((y == 1) & (bins == b), passed)
            item.update(bin=b, supported=item["n"] >= minimum)
            entries.append(item)
        result[name] = dict(interior_edges=edges.tolist(), bins=entries)
    mass_edges = np.asarray(result["mass"]["interior_edges"])
    z_edges = np.asarray(result["vertex_z"]["interior_edges"])
    rb_edges = np.asarray(result["vertex_radius"]["interior_edges"])
    m = np.searchsorted(mass_edges, nuisance[:, 0], side="right")
    z = np.searchsorted(z_edges, nuisance[:, 1], side="right")
    r = np.searchsorted(rb_edges, nuisance[:, 2], side="right")
    cells = []
    for i in range(len(mass_edges) + 1):
        for j in range(len(z_edges) + 1):
            for k in range(len(rb_edges) + 1):
                item = fraction((y == 1) & (m == i) & (z == j) & (r == k), passed)
                item.update(mass_bin=i, z_bin=j, radius_bin=k, supported=item["n"] >= minimum)
                cells.append(item)
    result["joint_mass_z_radius"] = cells
    supported = [c["efficiency"] for c in cells if c["supported"]]
    result["supported_joint_cells"] = len(supported)
    result["supported_joint_efficiency_range"] = [min(supported), max(supported)] if supported else None
    result["all_joint_cells_supported"] = all(c["supported"] for c in cells)
    return result


def bootstrap_auc(y, score, ids, process, repeats=200):
    groups = [f"{p}:" + ":".join(str(int(x)) for x in identity) for p, identity in zip(process, ids)]
    _, inverse = np.unique(groups, return_inverse=True)
    n_groups = int(inverse.max()) + 1
    rng = np.random.default_rng(471)
    values = []
    for _ in range(repeats):
        counts = np.bincount(rng.integers(n_groups, size=n_groups), minlength=n_groups)
        w = counts[inverse]
        if len(np.unique(y[w > 0])) == 2:
            values.append(roc_auc_score(y, score, sample_weight=w))
    return np.quantile(values, [0.025, 0.975]).tolist() if values else [None, None]


def evaluate(score, y, nuisance, ids, process, strata, sampling_weights, threshold, train_nuisance, train_y):
    passed = score >= threshold
    result = dict(auc=float(roc_auc_score(y, score)), auc_ci95=bootstrap_auc(y, score, ids, process),
                  threshold=float(threshold), common_vertex=fraction(y == 1, passed, sampling_weights),
                  accidental=fraction(y == 0, passed, sampling_weights),
                  dependence=dependence(score, nuisance, y), per_process={}, per_stratum={})
    result["background_rejection"] = 1 - result["accidental"]["efficiency"]
    for p in np.unique(process):
        result["per_process"][str(p)] = {name: fraction((process == p) & (y == label), passed, sampling_weights)
                                         for label, name in ((0, "accidental"), (1, "common_vertex"))}
        result["per_process"][str(p)]["dependence"] = dependence(score[process == p], nuisance[process == p], y[process == p])
    for p in np.unique(strata):
        result["per_stratum"][str(p)] = {name: fraction((strata == p) & (y == label), passed, sampling_weights)
                                         for label, name in ((0, "accidental"), (1, "common_vertex"))}
    result["flatness"] = binned_efficiencies(nuisance, y, passed, train_nuisance, train_y)
    result["lifetime_independence_validated"] = False
    result["mass_independence_validated"] = False
    result["deployment_ready"] = False
    return result
