# Code style: AI-code. AI-generated scheduled reconstructed-only neural trainer.
"""Optimization-only penalty ramp; checkpoints always use the fixed full loss.

The original frozen neural_vertex trainer is unchanged. This function preserves
its train/validation-only API and returns the same reconstructed-X classifier.
There are no truth, test, file-loading or deployment paths. A 20-epoch ramp starts
at zero in epoch1 and reaches full strength in epoch21. Early stopping and best
checkpoint comparison always use full-strength validation, including epoch1.
"""
from __future__ import annotations

import math
import numpy as np
import torch
from torch import nn
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from adversarial import _labels_and_weights, _matrix
from features import FEATURE_NAMES, signed_log
from neural_vertex import (NeuralVertexClassifier, _bounded_groups, _category_arrays,
                           _category_penalty, _groups, _penalty, _scale_nuisance,
                           _state_hash, _transform)


def train_neural_vertex(X_train, y_train, w_train, Z_train,
                        X_val, y_val, w_val, Z_val, *, train_process=None,
                        val_process=None, seed=71, strength=.3, hidden=(64, 32),
                        epochs=160, batch_size=512, patience=25,
                        feature_transform="signedlog", weight_decay=.0001,
                        dropout=.05, penalty_scope="class", energy_strength=0.,
                        learning_rate=.001, min_delta=.00001, category_strength=0.,
                        train_category=None, val_category=None, warmup_epochs=20):
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
    if not isinstance(warmup_epochs, (int, np.integer)) or isinstance(warmup_epochs, (bool, np.bool_)) or warmup_epochs < 0:
        raise ValueError("warmup_epochs must be a nonnegative integer")
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
        # Only optimization changes with epoch; every checkpoint uses full strengths.
        ramp = min(1., (epoch - 1) / warmup_epochs) if warmup_epochs else 1.
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
                objective = objective + (ramp * strength) * penalty
            if energy_strength > 0:
                groups = _bounded_groups(_groups(yt[idx], wt[idx], et[idx]), rng)
                penalty, _ = _penalty(torch.sigmoid(logits), te[idx], tw[idx], groups)
                objective = objective + (ramp * energy_strength) * penalty
            if category_strength > 0:
                members = np.flatnonzero((yt[idx] == 1) & (wt[idx] > 0))
                members = _bounded_groups([("genuine", members)], rng)[0][1]
                penalty, _ = _category_penalty(torch.sigmoid(logits), tc[idx], tw[idx], members)
                objective = objective + (ramp * category_strength) * penalty
            if not torch.isfinite(objective):
                raise RuntimeError("nonfinite neural objective")
            optimizer.zero_grad()
            objective.backward()
            if any(p.grad is not None and not torch.isfinite(p.grad).all() for p in network.parameters()):
                raise RuntimeError("nonfinite neural gradient")
            optimizer.step()
        network.eval()
        record = {"epoch": epoch, "optimization_penalty_fraction": ramp,
                  "optimization_joint_strength": ramp * strength,
                  "optimization_energy_strength": ramp * energy_strength,
                  "optimization_category_strength": ramp * category_strength}
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
                mean_score = (score * w).sum() / w.sum()
                record[name + "_score_weighted_mean"] = float(mean_score)
                record[name + "_score_weighted_std"] = float(torch.sqrt(((score - mean_score).square() * w).sum() / w.sum()))
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
        "method": "scheduled_class_conditional_joint_disco_with_fixed_validation_checkpoint",
        "warmup_epochs": int(warmup_epochs),
        "optimization_schedule": "All penalty strengths multiplied by min(1, (epoch-1)/warmup_epochs); zero warmup means full strength from epoch1",
        "checkpoint_schedule": "Full fixed penalty strengths from epoch1, including warmup epochs",
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
