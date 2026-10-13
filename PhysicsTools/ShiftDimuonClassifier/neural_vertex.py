# Code style: AI-code. AI-generated reconstructed-only neural training utility.
"""Small DisCo benchmark with train-only transforms and validation checkpoints.

X contains reconstructed predictors only. Z is nuisance_transform(reco mass,
vertex z, radius), NOT generator position or lifetime. Process identifiers may
define training/validation penalty groups but never enter the predictor. There
is no test-set argument. This module does not open files or inspect truth.
"""
from __future__ import annotations

import hashlib
import math

import numpy as np
import torch
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from torch import nn

from adversarial import _labels_and_weights, _matrix
from disco import weighted_distance_correlation_squared
from features import FEATURE_NAMES, signed_log


def _transform(x, mode):
    return signed_log(x) if mode == "signedlog" else x


def _state_hash(network):
    digest = hashlib.sha256()
    for name, value in network.state_dict().items():
        digest.update(name.encode())
        digest.update(value.detach().cpu().contiguous().numpy().tobytes())
    return digest.hexdigest()


class NeuralVertexClassifier:
    """Trusted-local pickle model; prediction requires only raw reconstructed X."""

    def __init__(self, network, imputer, scaler, n_features, feature_transform):
        self.network = network.cpu().eval()
        self.imputer, self.scaler = imputer, scaler
        self.n_features_in_ = int(n_features)
        self.feature_transform = feature_transform
        self.classes_ = np.asarray([0, 1])

    def predict_proba(self, X):
        raw = _matrix(X, name="X", n_features=self.n_features_in_)
        if not len(raw):
            return np.empty((0, 2))
        x = self.scaler.transform(self.imputer.transform(_transform(raw, self.feature_transform)))
        self.network.eval()
        with torch.no_grad():
            probability = torch.sigmoid(self.network(torch.tensor(x, dtype=torch.float32)).ravel()).numpy()
        return np.column_stack((1 - probability, probability))


def _groups(labels, weights, nuisance, process=None):
    """Supported row groups; missing targets are excluded, never imputed."""
    valid = (weights > 0) & np.isfinite(nuisance).all(axis=1)
    groups = [("class_0", valid & (labels == 0)), ("class_1", valid & (labels == 1))]
    if process is not None:
        groups += [("genuine_process_" + str(p), valid & (labels == 1) & (process == p))
                   for p in np.unique(process)]
    return [(name, np.flatnonzero(mask)) for name, mask in groups]


def _penalty(scores, nuisance, weights, groups):
    """Mean class dCor² plus optional mean per-process genuine dCor²."""
    class_values, process_values, details = [], [], []
    for name, indices in groups:
        selected_weights = weights[indices]
        effective_rows = (float(selected_weights.sum().square() / selected_weights.square().sum())
                          if len(indices) else 0.)
        if len(indices) < 4:
            details.append({"group": name, "rows": len(indices), "effective_rows": effective_rows,
                            "dcor_squared": None, "skip_reason": "fewer than four complete positive-weight rows"})
            continue
        value = weighted_distance_correlation_squared(scores[indices], nuisance[indices], weights[indices])
        (class_values if name.startswith("class_") else process_values).append(value)
        details.append({"group": name, "rows": len(indices), "effective_rows": effective_rows,
                        "dcor_squared": float(value.detach()), "skip_reason": None})
    zero = scores.sum() * 0
    result = torch.stack(class_values).mean() if class_values else zero
    if process_values:
        result = result + torch.stack(process_values).mean()
    return result, details


def _bounded_groups(groups, rng, limit=256):
    return [(name, np.sort(rng.choice(idx, limit, replace=False)) if len(idx) > limit else idx)
            for name, idx in groups]


def _category_penalty(scores, category, weights, indices):
    """Genuine-only category dCor²; one populated category is undefined."""
    selected_weights = weights[indices]
    effective_rows = (float(selected_weights.sum().square() / selected_weights.square().sum())
                      if len(indices) else 0.)
    populated = int((category[indices].sum(dim=0) > 0).sum())
    reason = None
    if len(indices) < 4:
        reason = "fewer than four positive-weight genuine rows"
    elif populated < 2:
        reason = "fewer than two populated genuine categories"
    value = (scores.sum() * 0 if reason else
             weighted_distance_correlation_squared(scores[indices], category[indices], weights[indices]))
    return value, {"rows": len(indices), "effective_rows": effective_rows,
                   "populated_categories": populated,
                   "dcor_squared": None if reason else float(value.detach()), "skip_reason": reason}


def _category_arrays(train, val, n_train, n_val):
    """Fit one-hot encoding on training strings, rejecting unseen validation."""
    train, val = np.asarray(train), np.asarray(val)
    if train.shape != (n_train,) or val.shape != (n_val,):
        raise ValueError("category penalty requires aligned train_category and val_category")
    if any(not isinstance(item, (str, np.str_)) or not item for values in (train, val) for item in values):
        raise ValueError("reconstructed categories must be nonempty strings")
    categories = np.unique(train)
    unseen = sorted(set(val) - set(categories))
    if unseen:
        raise ValueError("unseen validation reconstructed categories: " + ", ".join(unseen))
    return ((train[:, None] == categories).astype(float),
            (val[:, None] == categories).astype(float), categories.tolist())


def _scale_nuisance(train, val, weights):
    valid = (weights > 0) & np.isfinite(train).all(axis=1)
    if not valid.any():
        return train, val, None
    scaler = StandardScaler().fit(train[valid], sample_weight=weights[valid])
    return ((train - scaler.mean_) / scaler.scale_, (val - scaler.mean_) / scaler.scale_,
            {"mean": scaler.mean_.tolist(), "scale": scaler.scale_.tolist(),
             "complete_training_rows": int(valid.sum())})


def train_neural_vertex(X_train, y_train, w_train, Z_train,
                        X_val, y_val, w_val, Z_val, *, train_process=None,
                        val_process=None, seed=71, strength=.3, hidden=(64, 32),
                        epochs=160, batch_size=512, patience=25,
                        feature_transform="signedlog", weight_decay=.0001,
                        dropout=.05, penalty_scope="class", energy_strength=0.,
                        learning_rate=.001, min_delta=.00001, category_strength=0.,
                        train_category=None, val_category=None):
    """Return (restored best model, history), using train/validation arrays only.

Balanced minibatches sample each class uniformly with replacement. Their BCE
weights correct the sampling probabilities, preserving the supplied empirical
training objective. dCor uses the original weights within each class/group.
Validation selection minimizes weighted BCE + fixed strength * joint dCor²
+ energy_strength * energy dCor² + category_strength * genuine category dCor²
on fixed, bounded validation subsets. The
optional process term adds equal-group mean genuine dCor² to the class mean.
Early stopping counts epochs without a min_delta improvement; the best actual
objective is always saved, even when its improvement is below min_delta.
Energy protection uses the two reconstructed energy columns of canonical X50;
the separate joint energy penalty also conditions on class. No target/metadata
or nuisance scaler is stored in the returned inference model.
Optional reconstructed topology categories use train-fitted one-hot encoding
(Euclidean distance zero or sqrt(2)); unseen validation values fail explicitly.
"""
    for name, value in (("strength", strength), ("energy_strength", energy_strength),
                        ("weight_decay", weight_decay), ("min_delta", min_delta),
                        ("category_strength", category_strength)):
        if not math.isfinite(value) or value < 0:
            raise ValueError(name + " must be finite and nonnegative")
    if not math.isfinite(learning_rate) or learning_rate <= 0:
        raise ValueError("learning_rate must be finite and positive")
    if not math.isfinite(dropout) or not 0 <= dropout < 1:
        raise ValueError("dropout must be in [0, 1)")
    for name, value, minimum in (("epochs", epochs, 1), ("patience", patience, 1), ("batch_size", batch_size, 8)):
        if not isinstance(value, (int, np.integer)) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    if feature_transform not in ("raw", "signedlog") or penalty_scope not in ("class", "class_process"):
        raise ValueError("invalid feature_transform or penalty_scope")
    if not hidden or any(not isinstance(n, (int, np.integer)) or n < 1 for n in hidden):
        raise ValueError("hidden must contain positive integer widths")
    xt = _matrix(X_train, name="X_train")
    xv = _matrix(X_val, name="X_val", n_features=xt.shape[1])
    yt, wt = _labels_and_weights(y_train, w_train, len(xt), "training")
    yv, wv = _labels_and_weights(y_val, w_val, len(xv), "validation")
    for name, y, w in (("training", yt, wt), ("validation", yv, wv)):
        if len(np.unique(y[w > 0])) != 2:
            raise ValueError(name + " requires positive-weight examples of both classes")
    zt, zv = _matrix(Z_train, name="Z_train", n_features=3), _matrix(Z_val, name="Z_val", n_features=3)
    if len(zt) != len(xt) or len(zv) != len(xv):
        raise ValueError("nuisances must align with classifier rows")
    ct, cv, categories = None, None, None
    if category_strength > 0:
        ct, cv, categories = _category_arrays(train_category, val_category, len(xt), len(xv))
    pt = pv = None
    if penalty_scope == "class_process":
        pt, pv = np.asarray(train_process), np.asarray(val_process)
        if pt.shape != yt.shape or pv.shape != yv.shape:
            raise ValueError("class_process requires aligned train_process and val_process")
        pt, pv = pt.astype(str), pv.astype(str)
    energy_scaling = None
    et, ev = np.empty((len(xt), 0)), np.empty((len(xv), 0))
    if energy_strength > 0:
        if xt.shape[1] != len(FEATURE_NAMES):
            raise ValueError("energy penalty requires canonical reconstructed 50-feature X")
        energy_columns = [FEATURE_NAMES.index("mu_" + side + "_energy") for side in ("highPt", "lowPt")]
        et, ev, energy_scaling = _scale_nuisance(signed_log(xt[:, energy_columns]), signed_log(xv[:, energy_columns]), wt)
    zt, zv, nuisance_scaling = _scale_nuisance(zt, zv, wt)
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    scaler = StandardScaler()
    tx = torch.tensor(scaler.fit_transform(imputer.fit_transform(_transform(xt, feature_transform))), dtype=torch.float32)
    vx = torch.tensor(scaler.transform(imputer.transform(_transform(xv, feature_transform))), dtype=torch.float32)
    torch.set_num_threads(1)
    torch.manual_seed(int(seed))
    rng = np.random.default_rng(seed)
    layers = []
    width = xt.shape[1]
    for next_width in hidden:
        layers.extend([nn.Linear(width, int(next_width)), nn.ReLU(), nn.Dropout(dropout)])
        width = int(next_width)
    layers.append(nn.Linear(width, 1))
    network = nn.Sequential(*layers)
    optimizer = torch.optim.AdamW(network.parameters(), lr=learning_rate, weight_decay=weight_decay)
    bce = nn.BCEWithLogitsLoss(reduction="none")
    ty, vy = torch.tensor(yt, dtype=torch.float32), torch.tensor(yv, dtype=torch.float32)
    tw, vw = torch.tensor(wt, dtype=torch.float64), torch.tensor(wv, dtype=torch.float64)
    tz, vz, te, ve = map(torch.from_numpy, (zt, zv, et, ev))
    tc, vc = (None, None) if ct is None else (torch.from_numpy(ct), torch.from_numpy(cv))
    diagnostic_rng = np.random.default_rng(int(seed) + 130363)
    diagnostic_groups, support = {}, {}
    for name, y, w, z, e, p in (("train", yt, wt, zt, et, pt), ("val", yv, wv, zv, ev, pv)):
        all_groups = _groups(y, w, z, p)
        support[name] = [{"group": group, "complete_positive_weight_rows": len(idx)} for group, idx in all_groups]
        diagnostic_groups[name] = _bounded_groups(all_groups, diagnostic_rng)
        diagnostic_groups[name + "_energy"] = _bounded_groups(_groups(y, w, e), diagnostic_rng)
    if category_strength > 0:
        for name, y, w in (("train", yt, wt), ("val", yv, wv)):
            members = np.flatnonzero((y == 1) & (w > 0))
            diagnostic_groups[name + "_category"] = _bounded_groups([("genuine", members)], diagnostic_rng)[0][1]
    pools = [np.flatnonzero((yt == label) & (wt > 0)) for label in (0, 1)]
    sizes = [batch_size // 2, batch_size - batch_size // 2]
    updates = math.ceil(sum(map(len, pools)) / batch_size)
    records, best_state, best_prediction = [], None, None
    best_objective = stopping_reference = math.inf
    stale, best_epoch = 0, 0
    for epoch in range(1, epochs + 1):
        network.train()
        for _ in range(updates):
            draws = [rng.choice(pool, size, replace=True) for pool, size in zip(pools, sizes)]
            idx = np.concatenate(draws)
            correction = np.concatenate([np.full(size, len(pool) / size) for pool, size in zip(pools, sizes)])
            logits = network(tx[idx]).ravel()
            classification = (bce(logits, ty[idx]) * tw[idx] * torch.from_numpy(correction)).sum() / tw.sum()
            objective = classification
            if strength > 0:
                groups = _bounded_groups(_groups(yt[idx], wt[idx], zt[idx], pt[idx] if pt is not None else None), rng)
                penalty, _ = _penalty(torch.sigmoid(logits), tz[idx], tw[idx], groups)
                objective = objective + strength * penalty
            if energy_strength > 0:
                groups = _bounded_groups(_groups(yt[idx], wt[idx], et[idx]), rng)
                penalty, _ = _penalty(torch.sigmoid(logits), te[idx], tw[idx], groups)
                objective = objective + energy_strength * penalty
            if category_strength > 0:
                members = np.flatnonzero((yt[idx] == 1) & (wt[idx] > 0))
                members = _bounded_groups([("genuine", members)], rng)[0][1]
                penalty, _ = _category_penalty(torch.sigmoid(logits), tc[idx], tw[idx], members)
                objective = objective + category_strength * penalty
            if not torch.isfinite(objective):
                raise RuntimeError("nonfinite neural objective")
            optimizer.zero_grad()
            objective.backward()
            if any(p.grad is not None and not torch.isfinite(p.grad).all() for p in network.parameters()):
                raise RuntimeError("nonfinite neural gradient")
            optimizer.step()
        network.eval()
        record = {"epoch": epoch}
        with torch.no_grad():
            for name, x, y, w, z, e, c in (("train", tx, ty, tw, tz, te, tc), ("val", vx, vy, vw, vz, ve, vc)):
                logits = network(x).ravel()
                score = torch.sigmoid(logits)
                loss = float((bce(logits, y) * w).sum() / w.sum())
                penalty, details = _penalty(score, z, w, diagnostic_groups[name])
                energy_penalty = score.sum() * 0
                energy_details = []
                if energy_strength > 0:
                    energy_penalty, energy_details = _penalty(score, e, w, diagnostic_groups[name + "_energy"])
                category_penalty, category_details = score.sum() * 0, None
                if category_strength > 0:
                    category_penalty, category_details = _category_penalty(score, c, w, diagnostic_groups[name + "_category"])
                record[name + "_bce"] = loss
                record[name + "_dcor_squared"] = float(penalty)
                record[name + "_energy_dcor_squared"] = float(energy_penalty)
                record[name + "_groups"] = details
                record[name + "_energy_groups"] = energy_details
                record[name + "_category_dcor_squared"] = (category_details["dcor_squared"] if category_details else None)
                record[name + "_category_diagnostic"] = category_details
                record[name + "_objective"] = (loss + strength * float(penalty) + energy_strength * float(energy_penalty)
                                               + category_strength * float(category_penalty))
                if name == "val":
                    record["val_auc"] = float(roc_auc_score(yv, score.numpy(), sample_weight=wv))
                    val_prediction = score.numpy().copy()
        value = record["val_objective"]
        if not math.isfinite(value):
            raise RuntimeError("nonfinite validation objective")
        if value < best_objective:
            best_objective, best_epoch = value, epoch
            best_state = {key: value.detach().clone() for key, value in network.state_dict().items()}
            best_prediction = val_prediction
            best_hash = _state_hash(network)
        if value < stopping_reference - min_delta:
            stopping_reference, stale = value, 0
        else:
            stale += 1
        record["state_sha256"] = _state_hash(network)
        records.append(record)
        if stale >= patience:
            break
    final_hash = _state_hash(network)
    network.load_state_dict(best_state)
    model = NeuralVertexClassifier(network, imputer, scaler, xt.shape[1], feature_transform)
    restored_hash = _state_hash(network)
    if restored_hash != best_hash or not np.array_equal(model.predict_proba(xv)[:, 1], best_prediction):
        raise RuntimeError("best validation checkpoint restoration failed")
    history = {
        "method": "class_conditional_joint_disco_with_validation_checkpoint",
        "seed": int(seed), "architecture": [xt.shape[1], *map(int, hidden), 1],
        "strength": float(strength), "energy_strength": float(energy_strength),
        "category_strength": float(category_strength), "training_categories": categories,
        "penalty_scope": penalty_scope, "feature_transform": feature_transform,
        "dropout": float(dropout), "weight_decay": float(weight_decay),
        "learning_rate": float(learning_rate), "epochs_requested": int(epochs),
        "epochs_completed": len(records), "patience": int(patience), "min_delta": float(min_delta),
        "batch_size": int(batch_size), "updates_per_epoch": updates,
        "sampling": "equal class draws with replacement; inverse-probability BCE correction",
        "checkpoint_objective": "weighted validation BCE + strength * joint dCor2 + energy_strength * energy dCor2 + category_strength * genuine category dCor2",
        "nuisance_scaling": nuisance_scaling, "energy_scaling": energy_scaling,
        "diagnostic_max_rows_per_group": 256,
        "penalty_estimator": "biased weighted empirical dCor squared; not an independence test",
        "support": support,
        "best_epoch": best_epoch, "best_val_objective": best_objective,
        "best_state_sha256": best_hash, "last_state_sha256": final_hash,
        "restored_state_sha256": restored_hash, "restored_prediction_exact": True,
        "epochs": records,
    }
    return model, history
