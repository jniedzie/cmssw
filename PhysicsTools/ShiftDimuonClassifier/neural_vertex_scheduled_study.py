#!/usr/bin/env python3
# Code style: AI-code. AI-generated sealed validation-only scheduled DNN phase.
"""Twelve fixed warm-up recipes; never finalize or infer DNN test scores."""
import argparse
import json
from pathlib import Path
import pickle
import shutil
import time
import numpy as np
from sklearn.metrics import roc_auc_score
from features import FEATURE_NAMES, nuisance_transform
from phase2_train import phase_weights
from neural_vertex_study import prepare, reference, require, write, load, summary
from neural_vertex_scheduled import train_neural_vertex


def require_no_test_scores(study):
    require(not list(study.glob('training_round*/development_scores.npz')),
            'Scheduled validation must precede all DNN test-score extraction')


def recipes():
    return [{'name': f'dnn_r4_l{joint}_c{category}_e0_h64_s{seed}',
             'strength': float(joint), 'category_strength': float(category), 'energy_strength': 0.,
             'hidden': [64, 32], 'feature_transform': 'signedlog',
             'penalty_scope': 'class_process', 'seed': seed}
            for joint in (1, 3) for category in (1, 3, 10) for seed in (71, 72)]


def prepare_scheduled(workspace, output):
    require_no_test_scores(output.parent)
    protocol, *_ = prepare(workspace, output)
    (output / 'protocol.json').rename(output / 'initial_preparation_protocol.json')
    protocol['recipes'] = recipes()
    protocol['fit_budget'] = len(protocol['recipes'])
    protocol['training'].update(epochs=240, patience=35, learning_rate=.0005, dropout=.1,
                                weight_decay=.0005, batch_size=512, warmup_epochs=20)
    protocol['adaptation'] = {
        'uses': 'Only earlier training and validation evidence; no DNN development scores exist',
        'reason': 'Predeclared joint1/3 x category1/3/10 x seed71/72 stability comparison. Strong penalties showed seed-dependent weak validation discrimination; test an optimization-only20epoch ramp.',
        'no_previous_DNN_test_scores': True,
        'base_preparation': reference(output / 'initial_preparation_protocol.json'),
        'global_round1_round2_validation_review': reference(output.parent / 'global_validation_review.json'),
        'strong_penalty_validation_review': reference(output.parent / 'strong_penalty_validation_review.json'),
        'scheduled_trainer_review': reference(output.parent / 'scheduled_trainer_review.json')}
    protocol['checkpoint'] = 'Full-strength weighted validation BCE + joint and genuine-category penalties at EVERY epoch, including warmup; restore exact minimum checkpoint. No DNN test predictions.'
    for filename in (Path(__file__).name, 'neural_vertex_scheduled.py'):
        path = Path(__file__).parent / filename
        protocol['sources'][filename] = reference(path)
        shutil.copy2(path, output / 'source_snapshot' / filename)
    write(output / 'protocol.json', protocol)
    print(json.dumps({'prepared': True, 'fits': len(protocol['recipes']),
                      'protocol': reference(output / 'protocol.json')}), flush=True)


def fit_one(output, spec):
    require_no_test_scores(output.parent)
    protocol = json.loads((output / 'protocol.json').read_text())
    for ref in protocol['sources'].values():
        require(reference(ref['path'])['sha256'] == ref['sha256'], 'Training source changed after protocol freeze')
    for key in ('cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits'):
        ref = protocol[key]
        require(reference(ref['path'])['sha256'] == ref['sha256'], 'Fitting input changed: ' + key)
    reco, labels = load(protocol['reco']['path']), load(protocol['labels']['path'])
    split = load(output / 'splits.npz')['expanded_split']
    sm = reco['process'] != 'signal'
    tr, va = [(split == part) & sm & (labels['y'] >= 0) for part in (0, 1)]
    X, y, nuisance, process = reco['X'], labels['y'], reco['nuisance'], reco['process']
    Z = nuisance_transform(nuisance)
    wtr, wva = [phase_weights(y[mask], process[mask]) for mask in (tr, va)]
    name = spec['name']
    require(not (output / (name + '_validation.json')).exists(), 'Refusing to refit an archived candidate')
    started = time.monotonic()
    options = {key: value for key, value in spec.items() if key != 'name'}
    model, history = train_neural_vertex(X[tr], y[tr], wtr, Z[tr], X[va], y[va], wva, Z[va],
        train_process=process[tr], val_process=process[va], train_category=reco['topology'][tr],
        val_category=reco['topology'][va], **protocol['training'], **options)
    bundle = {'model': model, 'columns': np.arange(len(FEATURE_NAMES)), 'feature_names': FEATURE_NAMES,
              'requires_gen': False, 'deployment_ready': False}
    validation = summary(model.predict_proba(X[va])[:, 1], y[va], process[va], nuisance[va], nuisance[tr], y[tr],
                         category=reco['topology'][va])
    bundle['threshold'] = validation['threshold']
    train_score = model.predict_proba(X[tr])[:, 1]
    training = {'auc': float(roc_auc_score(y[tr], train_score)),
                'genuine_retention': float(np.mean(train_score[y[tr] == 1] >= bundle['threshold'])),
                'accidental_rejection': float(np.mean(train_score[y[tr] == 0] < bundle['threshold']))}
    with (output / 'models' / (name + '.pkl')).open('xb') as stream:
        pickle.dump(bundle, stream)
    write(output / (name + '_history.json'), history)
    result = {'recipe': spec, 'fit_seconds': time.monotonic() - started, 'validation': validation,
              'training_at_validation_cut': training, 'training_minus_validation_auc': training['auc'] - validation['auc'],
              'model': reference(output / 'models' / (name + '.pkl')), 'history': reference(output / (name + '_history.json')),
              'protocol': reference(output / 'protocol.json')}
    write(output / (name + '_validation.json'), result)
    print(json.dumps({'model': name, 'auc': validation['auc'], 'fake_rejection': validation['rejection'],
        'maximum_acceptance_span': validation['maximum_marginal_efficiency_span'], 'best_epoch': history.get('best_epoch'),
        'seconds': result['fit_seconds']}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--action', choices=('prepare', 'fit', 'all'), default='all')
    args = parser.parse_args()
    workspace, output = args.workspace.resolve(), args.output.resolve()
    require(output.is_relative_to(workspace / 'tea_shift_cmssw/plots'), 'Study belongs in project plots')
    if args.action in ('prepare', 'all'):
        prepare_scheduled(workspace, output)
    if args.action in ('fit', 'all'):
        protocol = json.loads((output / 'protocol.json').read_text())
        for spec in protocol['recipes']:
            if not (output / (spec['name'] + '_validation.json')).exists():
                fit_one(output, spec)
        require_no_test_scores(output.parent)
        results = []
        for spec in protocol['recipes']:
            path = output / (spec['name'] + '_validation.json')
            result = json.loads(path.read_text())
            require(result['protocol'] == reference(output / 'protocol.json'), 'Candidate protocol changed')
            for key in ('model', 'history'):
                require(reference(result[key]['path']) == result[key], 'Candidate artifact changed: ' + key)
            results.append(reference(path))
        write(output / 'validation_complete.json', {'complete': True, 'fits': len(results),
            'protocol': reference(output / 'protocol.json'), 'validation_receipts': results,
            'DNN_test_scores_evaluated': False, 'finalized': False, 'production_deployed': False})
        print(json.dumps({'validation_complete': True, 'fits': len(results), 'DNN_test_scores_evaluated': False}), flush=True)


if __name__ == '__main__':
    main()
