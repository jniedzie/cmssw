# Code style: AI-code. AI-generated reconstructed-category neural benchmark.
"""Pooled or shared-representation category models; no test or truth inputs.

X and category are reconstructed observations. Z is reconstructed mass/z/r
after nuisance_transform, used only for regularization. Labels may be truth
derived; process identifiers and weights are training/evaluation metadata only.
Inference needs exactly X and its reconstructed category, never Z/process/truth.
All model families require ONE global validation cut, supplied by the caller.
"""
from __future__ import annotations

import math
import numpy as np
import torch
from torch import nn
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

from adversarial import _matrix, _labels_and_weights
from disco import weighted_distance_correlation_squared
from vertex_acceptance_loss import (build_acceptance_plan, acceptance_memberships, training_cut,
                                    soft_acceptance_penalty, acceptance_diagnostics)
from neural_vertex import (_transform, _state_hash, _scale_nuisance, _groups, _penalty,
                           _bounded_groups, _category_arrays, _category_penalty)


def _feature_names(names, width):
    names = tuple(map(str, names))
    if len(names) != width or len(set(names)) != width:
        raise ValueError('Explicit unique feature_names must match X columns')
    forbidden = ('gen', 'truth', 'simtrack', 'simhit', 'lifetime', 'ctau', 'process', 'weight', 'sample')
    for name in names:
        lower = name.lower()
        if lower in ('run', 'lumi', 'event', 'luminosityblock', 'event_id', 'source_index', 'pair_index', 'source_id', 'event_index') or any(token in lower for token in forbidden):
            raise ValueError('Non-reconstructed predictor name: ' + name)
    return names


def _encode(categories, vocabulary, rows):
    values = np.asarray(categories)
    if values.shape != (rows,):
        raise ValueError('Reconstructed categories must align with X rows')
    if any(not isinstance(value, (str, np.str_)) or not value for value in values):
        raise ValueError('Reconstructed categories must be nonempty strings')
    unseen = sorted(set(values) - set(vocabulary))
    if unseen:
        raise ValueError('Unseen reconstructed categories: ' + ', '.join(unseen))
    return (values[:, None] == np.asarray(vocabulary)).astype(np.float32)


class CategoryNetwork(nn.Module):
    """Small heads share a trunk and a common logit; residual heads start at zero."""
    def __init__(self, width, categories, architecture, hidden, dropout):
        super().__init__()
        self.architecture = architecture
        if architecture == 'shared':
            width += categories
        layers = []
        for next_width in (() if architecture == 'linear_experts' else hidden):
            layers.extend((nn.Linear(width, int(next_width)), nn.ReLU(), nn.Dropout(dropout)))
            width = int(next_width)
        self.trunk = nn.Sequential(*layers)
        self.common = nn.Linear(width, 1)
        self.heads = None
        if architecture != 'shared':
            self.heads = nn.Linear(width, categories)
            nn.init.zeros_(self.heads.weight)
            nn.init.zeros_(self.heads.bias)

    def forward(self, x, category):
        hidden = self.trunk(torch.cat((x, category), dim=1) if self.architecture == 'shared' else x)
        score = self.common(hidden).ravel()
        if self.heads is not None:
            score = score + (self.heads(hidden) * category).sum(dim=1)
        return score


class CategoryVertexClassifier:
    """Trusted-local model; reconstructed X and category are the only inputs."""
    def __init__(self, network, imputer, scaler, names, categories, transform):
        self.network = network.cpu().eval()
        self.imputer, self.scaler = imputer, scaler
        self.feature_names_in_ = np.asarray(names)
        self.n_features_in_ = len(names)
        self.categories_ = tuple(categories)
        self.feature_transform = transform
        self.classes_ = np.asarray([0, 1])

    def predict_proba(self, X, categories):
        raw = _matrix(X, name='X', n_features=self.n_features_in_)
        onehot = _encode(categories, self.categories_, len(raw))
        if not len(raw):
            return np.empty((0, 2))
        values = self.scaler.transform(self.imputer.transform(_transform(raw, self.feature_transform)))
        self.network.eval()
        with torch.no_grad():
            score = torch.sigmoid(self.network(torch.tensor(values, dtype=torch.float32),
                                               torch.from_numpy(onehot))).numpy()
        return np.column_stack((1-score, score))


def conditional_groups(y, weights, z, category, process, specifications=None, minimum=20):
    """Group support is decided on training only; validation never invents groups."""
    valid = (weights > 0) & np.isfinite(z).all(axis=1)
    if specifications is None:
        specifications = [('category', int(label), str(cat), None)
                          for label in (0, 1) for cat in np.unique(category)]
        if process is not None:
            specifications += [('category_process', 1, str(cat), str(p))
                               for cat in np.unique(category) for p in np.unique(process)]
    groups, supported = [], []
    for kind, label, cat, p in specifications:
        mask = valid & (y == label) & (category == cat)
        if p is not None:
            mask &= process == p
        idx = np.flatnonzero(mask)
        if len(idx) >= minimum:
            specification = (kind, label, cat, p)
            supported.append(specification)
            groups.append((specification, idx))
    return groups, supported


def conditional_penalty(scores, z, weights, groups):
    """Mean within-class/category dCor² + mean genuine category/process dCor²."""
    values = {'category': [], 'category_process': []}
    details = []
    for specification, indices in groups:
        reason = 'fewer than four complete positive-weight rows' if len(indices) < 4 else None
        value = (scores.sum()*0 if reason else
                 weighted_distance_correlation_squared(scores[indices], z[indices], weights[indices]))
        if reason is None:
            values[specification[0]].append(value)
        selected = weights[indices]
        effective = float(selected.sum().square()/selected.square().sum()) if len(indices) else 0.
        details.append({'group': list(specification), 'rows': len(indices), 'effective_rows': effective,
                        'dcor_squared': None if reason else float(value.detach()), 'skip_reason': reason})
    result = scores.sum()*0
    for group_values in values.values():
        if group_values:
            result = result + torch.stack(group_values).mean()
    return result, details


def train_category_vertex(X_train, y_train, w_train, Z_train, X_val, y_val, w_val, Z_val, *,
                          feature_names, train_category, val_category, train_process=None, val_process=None,
                          architecture='shared', strength=.3, category_strength=0., conditional_strength=0.,
                          acceptance_strength=0., seed=71, hidden=(64, 32), epochs=160, batch_size=512,
                          patience=25, warmup_epochs=20, feature_transform='signedlog', weight_decay=.0005,
                          dropout=.1, learning_rate=.0005, min_delta=1e-5):
    """Train/validation-only optimization; return exactly restored best checkpoint.

Global class joint DisCo, genuine category DisCo and optional conditional DisCo
are separate coefficients. Conditional groups need20 complete positive-weight
TRAINING rows; bounded minibatches need4. Validation uses the same training-
defined groups and fixed subsets, even during warmup. Sparse categories retain
the shared representation and regularized residual heads, without large expert
networks. Uniform class sampling is corrected back to supplied BCE weights.
Optional acceptance loss uses train-defined bins and detached train cut only;
its unweighted group loss requires strictly positive phase weights. Inference
stores neither acceptance bins nor process metadata.
"""
    if architecture not in ('shared', 'heads', 'linear_experts'):
        raise ValueError('Unknown category architecture')
    for name, value in [('strength', strength), ('category_strength', category_strength),
                        ('conditional_strength', conditional_strength), ('acceptance_strength', acceptance_strength), ('weight_decay', weight_decay), ('min_delta', min_delta)]:
        if not math.isfinite(value) or value < 0:
            raise ValueError(name + ' must be finite and nonnegative')
    for name, value, minimum in [('epochs', epochs, 1), ('batch_size', batch_size, 8), ('patience', patience, 1), ('warmup_epochs', warmup_epochs, 0)]:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
            raise ValueError(name + ' has invalid integer range')
    if not math.isfinite(learning_rate) or learning_rate <= 0 or not math.isfinite(dropout) or not 0 <= dropout < 1:
        raise ValueError('Invalid learning_rate/dropout')
    if feature_transform not in ('raw', 'signedlog') or not hidden or any(int(n) != n or n < 1 for n in hidden):
        raise ValueError('Invalid feature transform/hidden widths')
    xt, xv = _matrix(X_train, name='X_train'), _matrix(X_val, name='X_val')
    names = _feature_names(feature_names, xt.shape[1])
    if xv.shape[1] != len(names):
        raise ValueError('Train/validation feature widths differ')
    yt, wt = _labels_and_weights(y_train, w_train, len(xt), 'training')
    yv, wv = _labels_and_weights(y_val, w_val, len(xv), 'validation')
    for y, w in ((yt, wt), (yv, wv)):
        if len(np.unique(y[w > 0])) != 2:
            raise ValueError('Both positive-weight classes are required')
    zt, zv = _matrix(Z_train, name='Z_train', n_features=3), _matrix(Z_val, name='Z_val', n_features=3)
    if len(zt) != len(xt) or len(zv) != len(xv):
        raise ValueError('Nuisance rows differ')
    ct, cv, vocabulary = _category_arrays(train_category, val_category, len(xt), len(xv))
    cat_t, cat_v = np.asarray(train_category), np.asarray(val_category)
    pt = pv = None
    if acceptance_strength > 0 or (conditional_strength > 0 and (train_process is not None or val_process is not None)):
        pt, pv = np.asarray(train_process), np.asarray(val_process)
        if pt.shape != yt.shape or pv.shape != yv.shape:
            raise ValueError('Conditional process metadata must align in both splits')
        pt, pv = pt.astype(str), pv.astype(str)
    acceptance_plan, acceptance_groups, acceptance_masks = None, {}, []
    if acceptance_strength > 0:
        if np.any(wt <= 0) or np.any(wv <= 0):
            raise ValueError('Unweighted acceptance loss requires strictly positive phase weights')
        acceptance_plan = build_acceptance_plan(yt, zt, cat_t, pt)
        acceptance_plan['coordinate_definition'] = 'incoming reco nuisance_transform: log1p mass, signedlog z, log1p radius'
        for split, y, z, cat, process in [('train', yt, zt, cat_t, pt), ('val', yv, zv, cat_v, pv)]:
            acceptance_groups[split] = acceptance_memberships(acceptance_plan, y, z, cat, process)
        for name, members in acceptance_groups['train']:
            mask = np.zeros(len(yt), dtype=bool)
            mask[members] = True
            acceptance_masks.append((name, mask))
    zt, zv, scaling = _scale_nuisance(zt, zv, wt)
    _, specifications = conditional_groups(yt, wt, zt, cat_t, pt)
    imputer = SimpleImputer(strategy='median', keep_empty_features=True)
    scaler = StandardScaler()
    tx = torch.tensor(scaler.fit_transform(imputer.fit_transform(_transform(xt, feature_transform))), dtype=torch.float32)
    vx = torch.tensor(scaler.transform(imputer.transform(_transform(xv, feature_transform))), dtype=torch.float32)
    tc, vc = torch.tensor(ct, dtype=torch.float32), torch.tensor(cv, dtype=torch.float32)
    ty, vy = torch.tensor(yt, dtype=torch.float32), torch.tensor(yv, dtype=torch.float32)
    tw, vw, tz, vz = map(torch.from_numpy, (wt, wv, zt, zv))
    torch.set_num_threads(1)
    torch.manual_seed(int(seed))
    rng, diagnostic_rng = np.random.default_rng(seed), np.random.default_rng(seed + 130363)
    network = CategoryNetwork(len(names), len(vocabulary), architecture, hidden, dropout)
    optimizer = torch.optim.AdamW(network.parameters(), lr=learning_rate, weight_decay=weight_decay)
    bce = nn.BCEWithLogitsLoss(reduction='none')
    diagnostic = {}
    for split, y, w, z, cat, p in [('train', yt, wt, zt, cat_t, pt), ('val', yv, wv, zv, cat_v, pv)]:
        conditional, _ = conditional_groups(y, w, z, cat, p, specifications, minimum=0)
        genuine = np.flatnonzero((y == 1) & (w > 0))
        diagnostic[split] = {'global': _bounded_groups(_groups(y, w, z), diagnostic_rng),
            'conditional': _bounded_groups(conditional, diagnostic_rng),
            'category': _bounded_groups([('genuine', genuine)], diagnostic_rng)[0][1]}
    pools = [np.flatnonzero((yt == label) & (wt > 0)) for label in (0, 1)]
    sizes = [batch_size//2, batch_size-batch_size//2]
    updates = math.ceil(sum(map(len, pools))/batch_size)
    records, best_state, best_prediction = [], None, None
    best_objective = stopping_reference = math.inf
    stale = 0
    for epoch in range(1, epochs + 1):
        ramp = min(1., (epoch-1)/warmup_epochs) if warmup_epochs else 1.
        epoch_cut = None
        if acceptance_strength > 0:
            network.eval()
            with torch.no_grad():
                epoch_cut = training_cut(torch.sigmoid(network(tx, tc)).numpy(), yt)
        network.train()
        for _ in range(updates):
            idx = np.concatenate([rng.choice(pool, size, replace=True) for pool, size in zip(pools, sizes)])
            correction = torch.tensor(np.concatenate([np.full(size, len(pool)/size) for pool, size in zip(pools, sizes)]))
            logits = network(tx[idx], tc[idx])
            score = torch.sigmoid(logits)
            objective = (bce(logits, ty[idx])*tw[idx]*correction).sum()/tw.sum()
            if strength > 0:
                penalty, _ = _penalty(score, tz[idx], tw[idx], _bounded_groups(_groups(yt[idx], wt[idx], zt[idx]), rng))
                objective = objective + ramp*strength*penalty
            if category_strength > 0:
                genuine = np.flatnonzero((yt[idx] == 1) & (wt[idx] > 0))
                genuine = _bounded_groups([('genuine', genuine)], rng)[0][1]
                penalty, _ = _category_penalty(score, tc[idx], tw[idx], genuine)
                objective = objective + ramp*category_strength*penalty
            if conditional_strength > 0:
                conditional, _ = conditional_groups(yt[idx], wt[idx], zt[idx], cat_t[idx], pt[idx] if pt is not None else None,
                                                     specifications, minimum=0)
                penalty, _ = conditional_penalty(score, tz[idx], tw[idx], _bounded_groups(conditional, rng))
                objective = objective + ramp*conditional_strength*penalty
            if acceptance_strength > 0:
                memberships = [(name, np.flatnonzero(mask[idx])) for name, mask in acceptance_masks]
                penalty, _ = soft_acceptance_penalty(score, np.flatnonzero(yt[idx] == 1), memberships, epoch_cut)
                objective = objective + ramp*acceptance_strength*penalty
            if not torch.isfinite(objective):
                raise RuntimeError('Nonfinite category-network objective')
            optimizer.zero_grad()
            objective.backward()
            if any(p.grad is not None and not torch.isfinite(p.grad).all() for p in network.parameters()):
                raise RuntimeError('Nonfinite category-network gradient')
            optimizer.step()
        network.eval()
        record = {'epoch': epoch, 'optimization_penalty_fraction': ramp, 'optimization_training_cut': epoch_cut}
        with torch.no_grad():
            for split, x, c, y, w, z in [('train', tx, tc, ty, tw, tz), ('val', vx, vc, vy, vw, vz)]:
                logits = network(x, c)
                score = torch.sigmoid(logits)
                loss = float((bce(logits, y)*w).sum()/w.sum())
                acceptance = score.sum()*0
                if acceptance_strength > 0:
                    if split == 'train':
                        checkpoint_cut = training_cut(score.numpy(), yt)
                        record['checkpoint_training_cut'] = checkpoint_cut
                    acceptance, acceptance_details = soft_acceptance_penalty(
                        score, np.flatnonzero(y.numpy() == 1), acceptance_groups[split], checkpoint_cut)
                    record[split+'_acceptance_support'] = acceptance_details
                record[split+'_soft_acceptance_penalty'] = float(acceptance)
                joint, joint_details = _penalty(score, z, w, diagnostic[split]['global'])
                category, category_details = _category_penalty(score, c, w, diagnostic[split]['category'])
                conditional, conditional_details = conditional_penalty(score, z, w, diagnostic[split]['conditional'])
                record.update({split+'_bce': loss, split+'_dcor_squared': float(joint),
                    split+'_category_dcor_squared': category_details['dcor_squared'],
                    split+'_conditional_dcor_squared': float(conditional), split+'_groups': joint_details,
                    split+'_category_diagnostic': category_details, split+'_conditional_groups': conditional_details,
                    split+'_objective': loss+strength*float(joint)+category_strength*float(category)+conditional_strength*float(conditional)+acceptance_strength*float(acceptance)})
                mean = (score*w).sum()/w.sum()
                record[split+'_score_weighted_std'] = float(torch.sqrt(((score-mean).square()*w).sum()/w.sum()))
                if split == 'val':
                    record['val_auc'] = float(roc_auc_score(yv, score.numpy(), sample_weight=wv))
                    val_prediction = score.numpy().copy()
        value = record['val_objective']
        if not math.isfinite(value):
            raise RuntimeError('Nonfinite validation objective')
        if value < best_objective:
            best_objective, best_epoch = value, epoch
            best_state = {key: tensor.detach().clone() for key, tensor in network.state_dict().items()}
            best_prediction, best_hash = val_prediction, _state_hash(network)
        if value < stopping_reference-min_delta:
            stopping_reference, stale = value, 0
        else:
            stale += 1
        record['state_sha256'] = _state_hash(network)
        records.append(record)
        if stale >= patience:
            break
    last_hash = _state_hash(network)
    network.load_state_dict(best_state)
    model = CategoryVertexClassifier(network, imputer, scaler, names, vocabulary, feature_transform)
    if _state_hash(network) != best_hash or not np.array_equal(model.predict_proba(xv, val_category)[:, 1], best_prediction):
        raise RuntimeError('Best category-network checkpoint restoration failed')
    history = {'method': 'reconstructed_category_network_with_scheduled_disco', 'architecture': architecture,
        'feature_names': list(names), 'hidden': list(hidden), 'training_categories': vocabulary,
        'seed': int(seed), 'strength': float(strength), 'category_strength': float(category_strength),
        'conditional_strength': float(conditional_strength), 'acceptance_strength': float(acceptance_strength),
        'nuisance_scaling': scaling, 'conditional_training_minimum_rows': 20,
        'conditional_training_groups': [list(spec) for spec in specifications],
        'conditional_definition': 'mean class/category joint dCor2 + mean genuine category/process joint dCor2; only train-supported groups',
        'warmup_epochs': int(warmup_epochs), 'checkpoint_schedule': 'full fixed coefficients from epoch1',
        'feature_transform': feature_transform, 'weight_decay': float(weight_decay), 'dropout': float(dropout),
        'learning_rate': float(learning_rate), 'batch_size': int(batch_size), 'updates_per_epoch': updates,
        'epochs_requested': int(epochs), 'epochs_completed': len(records), 'patience': int(patience), 'min_delta': float(min_delta),
        'best_epoch': best_epoch, 'best_val_objective': best_objective, 'best_state_sha256': best_hash,
        'last_state_sha256': last_hash, 'restored_state_sha256': _state_hash(network),
        'restored_prediction_exact': True, 'epochs': records,
        'inference_contract': 'reconstructed X plus reconstructed category; no nuisance/process/weight/truth arguments',
        'sampling': 'equal class draws with replacement and inverse-probability BCE correction',
        'penalty_estimator': 'biased weighted empirical joint dCor squared; not an independence test'}
    if acceptance_strength > 0:
        history['acceptance_plan'] = acceptance_plan
        history['acceptance_temperature'] = .05
        cut = records[best_epoch-1]['checkpoint_training_cut']
        history['best_acceptance_diagnostics'] = {
            split: acceptance_diagnostics(model.predict_proba(x, category)[:, 1], y, acceptance_groups[split], cut)
            for split, x, category, y in [('train', xt, cat_t, yt), ('val', xv, cat_v, yv)]}
    history['parameter_counts'] = {'total': sum(p.numel() for p in network.parameters()),
        'per_category_residual_head': (network.heads.in_features + 1) if network.heads is not None else 0}
    history['category_support'] = {str(cat): {split: {'genuine': int(((c == cat) & (y == 1) & (w > 0)).sum()),
        'fake': int(((c == cat) & (y == 0) & (w > 0)).sum())}
        for split, c, y, w in [('train', cat_t, yt, wt), ('val', cat_v, yv, wv)]} for cat in vocabulary}
    return model, history
