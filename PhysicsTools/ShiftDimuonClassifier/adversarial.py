"""Small CPU neural benchmark with reconstructed-only adversarial targets.

The classifier sees only ``X``. ``Z`` must contain reconstructed log1p mass,
signed-log1p vertex z, and log1p vertex radius, in that order. Two score-only
adversaries protect the joint distribution of these quantities separately in
each training-label class. Missing nuisance targets are excluded and reported;
they are never imputed into fictitious vertices. This is a neural benchmark,
not adversarial backpropagation through a conventional decision tree.
"""

from __future__ import annotations

import math

import numpy as np
import torch
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from torch import nn


def _classifier(n_features):
    return nn.Sequential(
        nn.Linear(n_features, 32), nn.ReLU(),
        nn.Linear(32, 16), nn.ReLU(), nn.Linear(16, 1),
    )


def _adversary():
    return nn.Sequential(
        nn.Linear(1, 32), nn.Tanh(),
        nn.Linear(32, 32), nn.Tanh(), nn.Linear(32, 27),
    )


def _matrix(values, *, name, n_features=None):
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 2 or (n_features is not None and result.shape[1] != n_features):
        raise ValueError(f"{name} must be a two-dimensional array with the expected columns")
    if result.shape[1] == 0:
        raise ValueError(f"{name} must contain at least one feature")
    # Infinities from an invalid reconstructed fit are missing observations.
    return np.where(np.isfinite(result), result, np.nan)


def _labels_and_weights(labels, weights, n_rows, name):
    y = np.asarray(labels)
    w = np.asarray(weights, dtype=np.float64)
    if y.shape != (n_rows,) or w.shape != (n_rows,):
        raise ValueError(f"{name} labels and weights must have one entry per row")
    if not np.isin(y, (0, 1)).all():
        raise ValueError(f"{name} labels must be binary")
    if not np.isfinite(w).all() or (w < 0).any() or w.sum() <= 0:
        raise ValueError(f"{name} weights must be finite, nonnegative, and have positive sum")
    return y.astype(np.int64), w / w.sum() * n_rows


def _weighted_quantiles(values, weights):
    order = np.argsort(values, kind="stable")
    ordered = values[order]
    cumulative = np.cumsum(weights[order]) / weights.sum()
    # Step quantiles preserve observed nuisance values and tolerate ties.
    return ordered[np.searchsorted(cumulative, (1.0 / 3.0, 2.0 / 3.0))]


def _joint_bins(z, edges):
    valid = np.isfinite(z).all(axis=1)
    result = np.full(len(z), -1, dtype=np.int64)
    if valid.any():
        indices = [np.searchsorted(edge, z[valid, j], side="right")
                   for j, edge in enumerate(edges)]
        result[valid] = 9 * indices[0] + 3 * indices[1] + indices[2]
    return result


def _entropy_and_counts(bins, weights):
    valid = (bins >= 0) & (weights > 0)
    if not valid.any():
        return None, 0
    counts = np.bincount(bins[valid], weights=weights[valid], minlength=27)
    probabilities = counts[counts > 0] / counts.sum()
    return float(-(probabilities * np.log(probabilities)).sum()), int(len(probabilities))


def _weighted_mean(loss, weights):
    return (loss * weights).sum() / weights.sum().clamp_min(1e-12)


class AdversarialClassifier:
    """Pickle-compatible local model; inference needs reconstructed X only.

    Pickle files must be loaded only from a trusted source. The module and its
    numpy/sklearn/torch dependencies must be importable when unpickling.
    """

    def __init__(self, network, imputer, scaler, n_features):
        self.network = network.cpu().eval()
        self.imputer = imputer
        self.scaler = scaler
        self.n_features_in_ = int(n_features)
        self.classes_ = np.array([0, 1], dtype=np.int64)

    def predict_proba(self, X):
        raw = _matrix(X, name="X", n_features=self.n_features_in_)
        if not len(raw):
            return np.empty((0, 2), dtype=np.float64)
        transformed = self.scaler.transform(self.imputer.transform(raw)).astype(np.float32)
        self.network.eval()
        with torch.no_grad():
            score = torch.sigmoid(self.network(torch.from_numpy(transformed))).numpy().ravel()
        return np.column_stack((1.0 - score, score))


def train_adversarial(X_train, y_train, w_train, Z_train,
                      X_val, y_val, w_val, Z_val, *, seed=42,
                      strength=0.3, epochs=100):
    """Return ``(model, history)`` for a deterministic small neural benchmark.

    ``history`` is a dict with ``epochs`` (per-epoch records), configuration,
    and ``adversaries`` (bin edges, valid counts, and skip reasons). Strength
    ramps linearly over the first 20% of epochs. A strength of zero gives the
    same neural architecture and classifier optimization without decorrelation.
    Inputs and weights are copied; no files or metadata are read.
    """
    if not math.isfinite(strength) or strength < 0:
        raise ValueError("strength must be finite and nonnegative")
    if not isinstance(epochs, (int, np.integer)) or epochs < 1:
        raise ValueError("epochs must be a positive integer")
    train_raw = _matrix(X_train, name="X_train")
    val_raw = _matrix(X_val, name="X_val", n_features=train_raw.shape[1])
    if not len(train_raw) or not len(val_raw):
        raise ValueError("training and validation sets must both be nonempty")
    yt, wt = _labels_and_weights(y_train, w_train, len(train_raw), "training")
    yv, wv = _labels_and_weights(y_val, w_val, len(val_raw), "validation")
    if len(np.unique(yt[wt > 0])) != 2:
        raise ValueError("training must have positive-weight examples in both classes")
    zt = _matrix(Z_train, name="Z_train", n_features=3)
    zv = _matrix(Z_val, name="Z_val", n_features=3)
    if len(zt) != len(yt) or len(zv) != len(yv):
        raise ValueError("Z must have one row per classifier example")

    torch.set_num_threads(1)
    torch.manual_seed(int(seed))
    rng = np.random.default_rng(seed)
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    scaler = StandardScaler()
    xt = scaler.fit_transform(imputer.fit_transform(train_raw)).astype(np.float32)
    xv = scaler.transform(imputer.transform(val_raw)).astype(np.float32)
    network = _classifier(xt.shape[1])
    optimizer = torch.optim.Adam(network.parameters(), lr=1e-3)
    bce = nn.BCEWithLogitsLoss(reduction="none")
    ce = nn.CrossEntropyLoss(reduction="none")
    tx = torch.from_numpy(xt)
    vx = torch.from_numpy(xv)
    ty = torch.from_numpy(yt.astype(np.float32))
    vy = torch.from_numpy(yv.astype(np.float32))
    tw = torch.from_numpy(wt.astype(np.float32))
    vw = torch.from_numpy(wv.astype(np.float32))
    heads = []
    metadata = []
    for label in (0, 1):
        members = (yt == label) & (wt > 0)
        valid = members & np.isfinite(zt).all(axis=1)
        n_valid = int(valid.sum())
        info = {"label": label, "train_rows": int(members.sum()),
                "train_valid_rows": n_valid,
                "train_missing_rows": int((members & ~np.isfinite(zt).all(axis=1)).sum()),
                "edges": None, "active": False, "skip_reason": None}
        if n_valid < 16:
            info["skip_reason"] = "fewer than 16 positive-weight rows with all nuisance values"
            metadata.append(info)
            continue
        edges = [_weighted_quantiles(zt[valid, j], wt[valid]) for j in range(3)]
        train_bins = _joint_bins(zt, edges)
        val_bins = _joint_bins(zv, edges)
        train_entropy, occupied = _entropy_and_counts(train_bins[members], wt[members])
        val_members = (yv == label) & (wv > 0)
        val_entropy, _ = _entropy_and_counts(val_bins[val_members], wv[val_members])
        info.update(edges=[edge.tolist() for edge in edges],
                    occupied_train_bins=occupied,
                    train_baseline_entropy=train_entropy,
                    val_baseline_entropy=val_entropy,
                    val_valid_rows=int((val_members & (val_bins >= 0)).sum()),
                    val_missing_rows=int((val_members & (val_bins < 0)).sum()))
        if occupied < 2:
            info["skip_reason"] = "nuisance variables have only one occupied joint bin"
            metadata.append(info)
            continue
        head = _adversary()
        info["active"] = True
        metadata.append(info)
        heads.append({"label": label, "network": head,
                      "optimizer": torch.optim.Adam(head.parameters(), lr=2e-3),
                      "train_bins": torch.from_numpy(train_bins),
                      "val_bins": torch.from_numpy(val_bins), "metadata": info})

    records = []
    batch_size = min(256, len(yt))
    ramp_epochs = max(1, int(math.ceil(epochs * 0.2)))
    for epoch in range(epochs):
        effective_strength = float(strength * min(1.0, epoch / ramp_epochs))
        network.train()
        for start in range(0, len(yt), batch_size):
            # Generate one permutation per epoch below, shared by both objectives.
            if start == 0:
                order = rng.permutation(len(yt))
            idx = torch.from_numpy(order[start:start + batch_size])
            if tw[idx].sum().item() <= 0:
                continue
            logits = network(tx[idx]).ravel()
            score = torch.sigmoid(logits).unsqueeze(1)
            usable = []
            for head in heads:
                mask = ((ty[idx] == head["label"]) &
                        (head["train_bins"][idx] >= 0) & (tw[idx] > 0))
                if mask.sum().item() < 4:
                    continue
                usable.append((head, mask))
                # The detached score prevents adversary updates from changing f.
                for parameter in head["network"].parameters():
                    parameter.requires_grad_(True)
                for _ in range(3):
                    head["optimizer"].zero_grad()
                    loss = _weighted_mean(
                        ce(head["network"](score[mask].detach()),
                           head["train_bins"][idx][mask]), tw[idx][mask])
                    loss.backward()
                    head["optimizer"].step()
            classification = _weighted_mean(bce(logits, ty[idx]), tw[idx])
            penalties = []
            for head, mask in usable:
                # Freeze r while retaining its derivative with respect to score.
                for parameter in head["network"].parameters():
                    parameter.requires_grad_(False)
                penalties.append(_weighted_mean(
                    ce(head["network"](score[mask]), head["train_bins"][idx][mask]),
                    tw[idx][mask]))
            penalty = torch.stack(penalties).mean() if penalties else logits.sum() * 0.0
            objective = classification - effective_strength * penalty
            optimizer.zero_grad()
            objective.backward()
            optimizer.step()

        network.eval()
        with torch.no_grad():
            train_logits = network(tx).ravel()
            val_logits = network(vx).ravel()
            train_score = torch.sigmoid(train_logits).unsqueeze(1)
            val_score = torch.sigmoid(val_logits).unsqueeze(1)
            train_bce = float(_weighted_mean(bce(train_logits, ty), tw).item())
            val_bce = float(_weighted_mean(bce(val_logits, vy), vw).item())
            train_penalties = []
            diagnostics = []
            for head in heads:
                diagnostic = {"label": head["label"]}
                for split, labels, weights, bins, scores in (
                    ("train", ty, tw, head["train_bins"], train_score),
                    ("val", vy, vw, head["val_bins"], val_score),
                ):
                    mask = (labels == head["label"]) & (bins >= 0) & (weights > 0)
                    value = None
                    if mask.any():
                        value = float(_weighted_mean(
                            ce(head["network"](scores[mask]), bins[mask]), weights[mask]).item())
                        if split == "train":
                            train_penalties.append(value)
                    baseline = head["metadata"][f"{split}_baseline_entropy"]
                    diagnostic[f"{split}_ce"] = value
                    diagnostic[f"{split}_baseline_entropy"] = baseline
                    diagnostic[f"{split}_entropy_minus_ce"] = (
                        baseline - value if baseline is not None and value is not None else None)
                diagnostics.append(diagnostic)
            valid_auc = len(np.unique(yv[wv > 0])) == 2
            auc = float(roc_auc_score(yv, val_score.numpy().ravel(), sample_weight=wv)) if valid_auc else None
            average_penalty = float(np.mean(train_penalties)) if train_penalties else 0.0
            records.append({"epoch": epoch + 1, "effective_strength": effective_strength,
                            "train_loss": train_bce - effective_strength * average_penalty,
                            "train_bce": train_bce, "val_bce": val_bce,
                            "val_auc": auc, "adversaries": diagnostics})

    model = AdversarialClassifier(network, imputer, scaler, train_raw.shape[1])
    history = {"seed": int(seed), "strength": float(strength), "epochs_requested": int(epochs),
               "architecture": [int(xt.shape[1]), 32, 16, 1],
               "nuisance_definition": ["reconstructed log1p mass", "reconstructed signedlog1p vz",
                                       "reconstructed log1p vertex radius"],
               "adversaries": metadata, "epochs": records}
    return model, history
