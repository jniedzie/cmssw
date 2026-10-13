#!/usr/bin/env python3
# Code style: AI-code. AI-generated bounded relative-geometry DNN comparison.
"""Test explicit track consistency inputs before any neural test-score review."""
import argparse
import json
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
from sklearn.metrics import roc_auc_score

from features import FEATURE_NAMES, nuisance_transform
from neural_vertex_geometry import train_geometry_vertex
from neural_vertex_study import load, prepare, reference, require, summary, write
from phase2_train import phase_weights


def recipes():
    result = []
    for mode in ('all', 'quality_relative'):
        for strength in (0., .3, 1., 3.):
            for seed in (71, 72):
                result.append({'name': f'dnn_geometry_{mode}_l{strength:g}_s{seed}', 'mode': mode,
                    'strength': strength, 'category_strength': 1., 'energy_strength': 0.,
                    'hidden': [64, 32], 'feature_transform': 'signedlog', 'penalty_scope': 'class_process', 'seed': seed})
    return result


def fit_geometry(output, spec):
    protocol = json.loads((output / 'protocol.json').read_text())
    for ref in protocol['sources'].values():
        require(reference(ref['path'])['sha256'] == ref['sha256'], 'Geometry source changed after freeze')
    for key in ('cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits'):
        ref = protocol[key]
        require(reference(ref['path'])['sha256'] == ref['sha256'], 'Geometry study input changed: ' + key)
    data, labels = load(protocol['reco']['path']), load(protocol['labels']['path'])
    split = load(protocol['splits']['path'])['expanded_split']
    y, X, Z = labels['y'], data['X'], nuisance_transform(data['nuisance'])
    tr, va = [(split == part) & (data['process'] != 'signal') & (y >= 0) for part in (0, 1)]
    weights = [phase_weights(y[mask], data['process'][mask]) for mask in (tr, va)]
    name = spec['name']
    require(not (output / (name + '_validation.json')).exists(), 'Refusing archived geometry refit')
    options = {key: value for key, value in spec.items() if key != 'name'}
    start = time.monotonic()
    model, history = train_geometry_vertex(X[tr], y[tr], weights[0], Z[tr], X[va], y[va], weights[1], Z[va],
        train_process=data['process'][tr], val_process=data['process'][va],
        train_category=data['topology'][tr], val_category=data['topology'][va], **protocol['training'], **options)
    validation = summary(model.predict_proba(X[va])[:, 1], y[va], data['process'][va], data['nuisance'][va],
        data['nuisance'][tr], y[tr], category=data['topology'][va])
    train_score = model.predict_proba(X[tr])[:, 1]
    train = {'auc': float(roc_auc_score(y[tr], train_score)),
        'genuine_retention': float(np.mean(train_score[y[tr] == 1] >= validation['threshold'])),
        'accidental_rejection': float(np.mean(train_score[y[tr] == 0] < validation['threshold']))}
    bundle = {'model': model, 'columns': np.arange(len(FEATURE_NAMES)), 'feature_names': FEATURE_NAMES,
              'threshold': validation['threshold'], 'requires_gen': False, 'deployment_ready': False}
    with (output / 'models' / (name + '.pkl')).open('xb') as stream:
        pickle.dump(bundle, stream)
    write(output / (name + '_history.json'), history)
    entry = {'recipe': spec, 'validation': validation, 'training_at_validation_cut': train,
        'training_minus_validation_auc': train['auc'] - validation['auc'], 'fit_seconds': time.monotonic() - start,
        'model': reference(output / 'models' / (name + '.pkl')), 'history': reference(output / (name + '_history.json')),
        'protocol': reference(output / 'protocol.json')}
    write(output / (name + '_validation.json'), entry)
    print(json.dumps({'model': name, 'rejection': validation['rejection'],
        'worst_span': validation['maximum_marginal_efficiency_span'], 'category_span': validation['category_efficiency_span'],
        'best_epoch': history['best_epoch'], 'seconds': entry['fit_seconds']}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--previous', type=Path, nargs='+', required=True)
    parser.add_argument('--action', choices=('prepare', 'fit'), default='prepare')
    args = parser.parse_args()
    workspace, output = args.workspace.resolve(), args.output.resolve()
    require(output.is_relative_to(workspace / 'tea_shift_cmssw/plots'), 'Final results belong in repository')
    if args.action == 'prepare':
        adaptations = []
        for previous in args.previous:
            require(not (previous / 'development_scores.npz').exists(), 'Geometry study must precede DNN test-score inspection')
            prior = json.loads((previous / 'protocol.json').read_text())
            adaptations.append({'protocol': reference(previous / 'protocol.json'),
                'validation': [reference(previous / (spec['name'] + '_validation.json')) for spec in prior['recipes']]})
        protocol, *_ = prepare(workspace, output)
        (output / 'protocol.json').rename(output / 'initial_preparation_protocol.json')
        protocol['recipes'] = recipes()
        protocol['fit_budget'] = len(protocol['recipes'])
        protocol['training'].update(epochs=240, patience=35, learning_rate=.0005, dropout=.1, weight_decay=.0005)
        protocol['adaptation'] = {'uses': 'Earlier validation only; no neural test scores', 'prior_rounds': adaptations,
            'reason': 'Penalties improve acceptance but can discard classification power. Add explicit measured track consistency; compare a quality/relative representation excluding all absolute energies, momenta and coordinates.',
            'base_preparation': reference(output / 'initial_preparation_protocol.json')}
        protocol['geometry_contract'] = 'Raw reconstructed X50 at inference. All mode adds12relative features; quality_relative uses26quality + originCompatibility +12relative features. Approximate lineDCA and scale ratios are not helixDCA or calibrated significances.'
        for path in (Path(__file__), Path(__file__).parent / 'neural_vertex_geometry.py'):
            protocol['sources'][path.name] = reference(path)
            shutil.copy2(path, output / 'source_snapshot' / path.name)
        write(output / 'protocol.json', protocol)
        print(json.dumps({'prepared': True, 'fits': len(protocol['recipes'])}), flush=True)
    else:
        protocol = json.loads((output / 'protocol.json').read_text())
        for spec in protocol['recipes']:
            if not (output / (spec['name'] + '_validation.json')).exists():
                fit_geometry(output, spec)
        print('Geometry validation fits complete; no neural test scores evaluated.', flush=True)


if __name__ == '__main__':
    main()
