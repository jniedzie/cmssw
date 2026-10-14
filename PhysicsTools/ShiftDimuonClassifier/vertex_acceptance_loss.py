# Code style: AI-code. AI-generated training-only working-point uniformity loss.
"""Training-defined genuine groups and a smooth fixed-cut acceptance penalty.

This is an optimization surrogate. Report hard-cut acceptance separately.
No generated mass, lifetime, ancestry, identity or sample weight enters it.
Process labels select training strata only; none enter model inference.
"""
import numpy as np
import torch


def build_acceptance_plan(y, nuisance, category, process, *, bin_minimum=20,
                          category_minimum=50, process_minimum=50):
    y, nuisance = np.asarray(y), np.asarray(nuisance, dtype=float)
    category, process = np.asarray(category).astype(str), np.asarray(process).astype(str)
    if nuisance.shape != (len(y), 3) or category.shape != y.shape or process.shape != y.shape:
        raise ValueError('Acceptance-plan arrays must align')
    if not np.isin(y, [0, 1]).all() or not np.isfinite(nuisance).all():
        raise ValueError('Acceptance plan requires known labels and finite reconstructed nuisances')
    if any(not isinstance(n, int) or n < 1 for n in (bin_minimum, category_minimum, process_minimum)):
        raise ValueError('Support minima must be positive integers')
    genuine = y == 1
    if not genuine.any():
        raise ValueError('Acceptance plan requires genuine training rows')
    edges = [np.quantile(nuisance[genuine, axis], [1 / 3, 2 / 3]).tolist() for axis in range(3)]
    descriptors = []
    def add(name, scope, value, axis=None, bin_index=None):
        item = dict(name=name, scope=scope, value=value, axis=axis, bin_index=bin_index)
        mask = descriptor_mask(item, edges, y, nuisance, category, process)
        minimum = bin_minimum if axis is not None else (category_minimum if scope == 'category' else process_minimum)
        count = int(mask.sum())
        if count >= minimum:
            item['training_genuine_rows'] = count
            descriptors.append(item)
    for axis in range(3):
        for bin_index in range(3):
            add(f'pooled_axis{axis}_bin{bin_index}', 'all', '', axis, bin_index)
    for scope, values, minimum in [('category', category, category_minimum), ('process', process, process_minimum)]:
        for value in np.unique(values[genuine]):
            mask = genuine & (values == value)
            if int(mask.sum()) < minimum:
                continue
            add(scope + '_' + value, scope, str(value))
            for axis in range(3):
                for bin_index in range(3):
                    add(f'{scope}_{value}_axis{axis}_bin{bin_index}', scope, str(value), axis, bin_index)
    return dict(schema='shift-train-defined-soft-acceptance-v1', edges=edges,
                bin_minimum=bin_minimum, category_minimum=category_minimum,
                process_minimum=process_minimum, groups=descriptors,
                objective='Mean squared unweighted genuine group-minus-pooled soft acceptance')


def descriptor_mask(item, edges, y, nuisance, category, process):
    mask = np.asarray(y) == 1
    if item['scope'] != 'all':
        values = category if item['scope'] == 'category' else process
        mask = mask & (np.asarray(values).astype(str) == item['value'])
    if item['axis'] is not None:
        axis = item['axis']
        bins = np.searchsorted(edges[axis], np.asarray(nuisance)[:, axis], side='right')
        mask = mask & (bins == item['bin_index'])
    return mask


def acceptance_memberships(plan, y, nuisance, category, process):
    return [(item['name'], np.flatnonzero(descriptor_mask(item, plan['edges'], y, nuisance, category, process)))
            for item in plan['groups']]


def training_cut(scores, y):
    scores, y = np.asarray(scores), np.asarray(y)
    genuine = scores[y == 1]
    if not len(genuine) or not np.isfinite(genuine).all():
        raise ValueError('A finite genuine TRAIN population is required for the cut')
    # Same empirical inclusive operating-point convention as the study metrics.
    ordered = np.sort(genuine, kind='stable')
    index = int(np.searchsorted(np.arange(1, len(ordered) + 1), .1 * len(ordered)))
    return float(ordered[index])


def soft_acceptance_penalty(scores, genuine_indices, memberships, cut,
                            *, temperature=.05, batch_minimum=4):
    if not np.isfinite(cut) or not 0 <= cut <= 1 or not np.isfinite(temperature) or temperature <= 0:
        raise ValueError('Invalid soft-acceptance cut or temperature')
    genuine_indices = np.asarray(genuine_indices, dtype=int)
    if len(genuine_indices) < batch_minimum:
        return scores.sum() * 0, {'groups_used': 0, 'pooled_soft_acceptance': None}
    soft = torch.sigmoid((scores - float(cut)) / float(temperature))
    pooled = soft[genuine_indices].mean()
    terms = [torch.square(soft[np.asarray(indices, dtype=int)].mean() - pooled)
             for _, indices in memberships if len(indices) >= batch_minimum]
    if not terms:
        return scores.sum() * 0, {'groups_used': 0, 'pooled_soft_acceptance': float(pooled.detach())}
    return torch.stack(terms).mean(), {'groups_used': len(terms), 'pooled_soft_acceptance': float(pooled.detach())}


def acceptance_diagnostics(scores, y, memberships, cut, *, temperature=.05):
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError('Invalid diagnostic temperature')
    scores, y = np.asarray(scores, dtype=float), np.asarray(y)
    if scores.shape != y.shape or not np.isfinite(scores).all():
        raise ValueError('Invalid diagnostic score shape or values')
    soft = 1 / (1 + np.exp(np.clip((float(cut) - scores) / temperature, -700, 700)))
    groups = {}
    for name, indices in memberships:
        groups[name] = dict(n=len(indices),
                           hard_acceptance=float(np.mean(scores[indices] >= cut)) if len(indices) else None,
                           soft_acceptance=float(np.mean(soft[indices])) if len(indices) else None)
    genuine = y == 1
    return dict(cut=float(cut), temperature=temperature, score_std=float(np.std(scores)),
                genuine_score_std=float(np.std(scores[genuine])),
                pooled_hard_acceptance=float(np.mean(scores[genuine] >= cut)),
                pooled_soft_acceptance=float(np.mean(soft[genuine])), groups=groups,
                interpretation='Surrogate diagnostics at detached TRAIN cut; validation operating cut remains separately chosen')
