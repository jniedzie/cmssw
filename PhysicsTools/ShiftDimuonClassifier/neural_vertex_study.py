#!/usr/bin/env python3
# Code style: AI-code. AI-generated bounded DNN study with validation-only selection.
"""Train checkpointed vertex DNNs on the preserved expanded SM event split."""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import shutil
import sys
import time

import numpy as np
from sklearn.metrics import roc_auc_score

from evaluation import event_split, weighted_quantile
from features import FEATURE_NAMES, nuisance_transform
from phase2_train import phase_weights, validation_summary


def reference(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def load(path):
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def candidates():
    recipes = []
    for transform, hidden in [('raw', [32, 16]), ('signedlog', [64, 32])]:
        for strength in (0., .3, 1.):
            for seed in (71, 72):
                recipes.append({'name': f'dnn_{transform}_class_l{strength:g}_s{seed}',
                    'feature_transform': transform, 'hidden': hidden, 'strength': strength,
                    'penalty_scope': 'class', 'energy_strength': 0., 'seed': seed})
    for strength in (.3, 1.):
        for seed in (71, 72):
            recipes.append({'name': f'dnn_signedlog_process_l{strength:g}_s{seed}',
                'feature_transform': 'signedlog', 'hidden': [64, 32], 'strength': strength,
                'penalty_scope': 'class_process', 'energy_strength': 0., 'seed': seed})
    for seed in (71, 72):
        recipes.append({'name': f'dnn_signedlog_process_energy_l0.3_s{seed}',
            'feature_transform': 'signedlog', 'hidden': [64, 32], 'strength': .3,
            'penalty_scope': 'class_process', 'energy_strength': .3, 'seed': seed})
    for energy_strength in (0., .3):
        for seed in (71, 72):
            recipes.append({'name': f'dnn_signedlog_process_category_e{energy_strength:g}_s{seed}',
                'feature_transform': 'signedlog', 'hidden': [64, 32], 'strength': .3,
                'penalty_scope': 'class_process', 'energy_strength': energy_strength,
                'category_strength': .3, 'seed': seed})
    return recipes


def prepare(workspace, output):
    cache = workspace / 'validation/bdt_acceptance_study_20261012/comparison/extraction'
    previous = workspace / 'validation/bdt_discrimination_improvement_20261012/expanded_training'
    manifest = json.loads((cache / 'inputs.json').read_text())
    completion = json.loads((previous / 'completion.json').read_text())
    require(manifest['complete'] and manifest['sample_kind'] == 'simulation_only', 'Incomplete simulation cache')
    for name in ('reco', 'labels'):
        require(reference(cache / (name + '.npz'))['sha256'] == manifest[name + '_sha256'], 'Changed cache: ' + name)
    require(reference(previous / 'development_scores.npz')['sha256'] == completion['development']['scores']['sha256'], 'Changed BDT scores')
    reco, labels, previous_scores = load(cache / 'reco.npz'), load(cache / 'labels.npz'), load(previous / 'development_scores.npz')
    require(tuple(reco['feature_names']) == FEATURE_NAMES, 'Predictor contract changed')
    for name in ('ids', 'process', 'pair_index', 'source_index'):
        require(np.array_equal(reco[name], previous_scores[name]), 'Misaligned BDT cache: ' + name)
    split = event_split(reco['ids'], reco['process'], 71)
    require(np.array_equal(split, previous_scores['expanded_split']), 'Historical expanded split changed')
    sm = np.isin(reco['process'], ('dy', 'jpsi', 'qcd'))
    test = (split == 2) & sm
    require(np.array_equal(test, previous_scores['expanded_heldout']), 'Historical test mask changed')
    require(reco['ids'].dtype == np.uint64, 'Event identities must remain uint64')
    event_sets = [{(str(p), *map(int, identity)) for p, identity in zip(reco['process'][sm & (split == part)], reco['ids'][sm & (split == part)])} for part in range(3)]
    require(all(not event_sets[a].intersection(event_sets[b]) for a, b in [(0, 1), (0, 2), (1, 2)]), 'An event crosses partitions')
    tr, va = [(split == part) & sm & (labels['y'] >= 0) for part in (0, 1)]
    require(int(tr.sum()) == 9826 and int(va.sum()) == 3400, 'Fitting population changed')
    output.mkdir(parents=True, exist_ok=False)
    (output / 'models').mkdir()
    (output / 'source_snapshot').mkdir()
    source = {}
    for filename in ('neural_vertex_study.py', 'neural_vertex.py', 'disco.py', 'adversarial.py', 'features.py',
                     'evaluation.py', 'phase2_train.py', 'inference.py', 'bdt_performance.py'):
        path = Path(__file__).parent / filename
        source[filename] = reference(path)
        shutil.copy2(path, output / 'source_snapshot' / filename)
    protocol = {'schema': 'shift-vertex-dnn-study-v1', 'recipes': candidates(), 'fit_budget': len(candidates()),
        'cache_manifest': reference(cache / 'inputs.json'), 'reco': reference(cache / 'reco.npz'),
        'labels': reference(cache / 'labels.npz'), 'bdt_scores': reference(previous / 'development_scores.npz'),
        'sources': source, 'predictors': list(FEATURE_NAMES), 'split_seed': 71,
        'train_known_pairs': int(tr.sum()), 'validation_known_pairs': int(va.sum()),
        'excluded_signal_pairs': int((~sm).sum()),
        'training': {'epochs': 160, 'batch_size': 512, 'patience': 25, 'learning_rate': .001,
                     'dropout': .05, 'weight_decay': .0001, 'min_delta': .00001},
        'weights': 'Historical equal-class, smoothed process-balanced training weights; native and exposure weights retained for diagnostics only',
        'nuisances': 'Joint reconstructed log1p mass, signedlog1p fitted z, log1p fitted radius. Optional joint reconstructed muon-energy penalty and genuine reconstructed-category penalty. No truth lifetime in loss or inference.',
        'checkpoint': 'Lowest deterministic validation BCE plus fixed joint, optional energy and optional genuine-category penalties; patience25; restore best state. No test predictions before recipe/model lock.',
        'threshold': 'Empirical validation genuine10th percentile, inclusive >=; ties explicit, about90percent retained',
        'selection': 'Validation-only: rejection winner, old supported-flatness objective winner, and highest rejection inside predeclared acceptance envelope',
        'acceptance_envelope': {'minimum_supported_genuine_retention': .8, 'maximum_supported_marginal_span': .15,
                                'maximum_process_efficiency_span': .05, 'maximum_class_joint_dcor': .3,
                                'maximum_category_efficiency_span': .15},
        'support': 'At least20 genuine pairs/bin and50/process or large category; QCD genuine validation count1 does not establish acceptance',
        'interpretation': 'Consumed development study; all existing test data previously inspected. Unknown labels excluded from fitting and retained separately. No production deployment or physics cut.'}
    np.savez_compressed(output / 'splits.npz', ids=reco['ids'], process=reco['process'], pair_index=reco['pair_index'],
                        source_index=reco['source_index'], expanded_split=split, expanded_heldout=test)
    protocol['splits'] = reference(output / 'splits.npz')
    write(output / 'protocol.json', protocol)
    return protocol, reco, labels, split, tr, va


def summary(score, y, process, nuisance, training_nuisance, training_y, threshold=None, category=None):
    result = validation_summary(score, y, process, nuisance, np.ones(len(y)), training_nuisance, training_y)
    if threshold is not None and threshold != result['threshold']:
        # Only train/test summaries use a frozen validation cut; dependence and binning stay unchanged.
        from evaluation import binned_efficiencies, fraction
        passed = score >= threshold
        result.update(threshold=float(threshold), common_vertex=fraction(y == 1, passed), accidental=fraction(y == 0, passed),
                      flatness=binned_efficiencies(nuisance, y, passed, training_nuisance, training_y))
        result['per_process'] = {str(p): {key: fraction((process == p) & (y == label), passed)
            for label, key in [(1, 'common_vertex'), (0, 'accidental')]} for p in np.unique(process)}
        spans = []
        process_efficiencies = []
        result['per_process_flatness'] = {}
        for variable in ('mass', 'vertex_z', 'vertex_radius'):
            values = [row['efficiency'] for row in result['flatness'][variable]['bins'] if row['supported']]
            spans.append(max(values) - min(values) if len(values) >= 2 else 1.)
        for p in np.unique(process):
            mask = process == p
            if int((mask & (y == 1)).sum()) < 50:
                continue
            process_efficiencies.append(result['per_process'][str(p)]['common_vertex']['efficiency'])
            profiles = binned_efficiencies(nuisance[mask], y[mask], passed[mask], training_nuisance, training_y)
            result['per_process_flatness'][str(p)] = profiles
            for variable in ('mass', 'vertex_z', 'vertex_radius'):
                values = [row['efficiency'] for row in profiles[variable]['bins'] if row['supported']]
                if len(values) >= 2:
                    spans.append(max(values) - min(values))
        result['maximum_marginal_efficiency_span'] = max(spans)
        result['process_efficiency_span'] = max(process_efficiencies) - min(process_efficiencies) if len(process_efficiencies) >= 2 else 0.
        joint = max(result['dependence'][label]['joint_dcor'] or 0. for label in ('common_vertex', 'accidental'))
        rejection = 1 - result['accidental']['efficiency']
        result['selector_objective'] = rejection - 2 * max(spans) - 2 * result['process_efficiency_span'] - .25 * joint
        result['balanced_common_efficiency'] = result['common_vertex']['efficiency']
    result['rejection'] = 1 - result['accidental']['efficiency']
    result['auc'] = float(roc_auc_score(y, score))
    if category is not None:
        from evaluation import fraction
        passed = score >= result['threshold']
        result['categories'] = {str(c): {key: fraction((category == c) & (y == label), passed)
            for label, key in [(1, 'common_vertex'), (0, 'accidental')]} for c in np.unique(category)}
        supported = [row['common_vertex']['efficiency'] for row in result['categories'].values()
                     if row['common_vertex']['n'] >= 50]
        result['minimum_supported_category_efficiency'] = min(supported) if supported else None
        result['category_efficiency_span'] = max(supported) - min(supported) if len(supported) >= 2 else None
    return result


def envelope_eligible(entry, bounds):
    spans = entry['maximum_marginal_efficiency_span']
    efficiencies = []
    for profiles in [entry['flatness'], *entry['per_process_flatness'].values()]:
        for variable in ('mass', 'vertex_z', 'vertex_radius'):
            efficiencies.extend(row['efficiency'] for row in profiles[variable]['bins'] if row['supported'])
    joint = max(entry['dependence'][label]['joint_dcor'] or 0. for label in ('common_vertex', 'accidental'))
    category_minimum = entry.get('minimum_supported_category_efficiency')
    category_span = entry.get('category_efficiency_span')
    return (bool(efficiencies) and min(efficiencies) >= bounds['minimum_supported_genuine_retention']
            and spans <= bounds['maximum_supported_marginal_span']
            and entry['process_efficiency_span'] <= bounds['maximum_process_efficiency_span']
            and joint <= bounds['maximum_class_joint_dcor']
            and category_minimum is not None and category_span is not None
            and category_minimum >= bounds['minimum_supported_genuine_retention']
            and category_span <= bounds['maximum_category_efficiency_span'])


def fit_one(output, spec):
    from neural_vertex import train_neural_vertex
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


def finalize(output):
    protocol = json.loads((output / 'protocol.json').read_text())
    for ref in protocol['sources'].values():
        require(reference(ref['path'])['sha256'] == ref['sha256'], 'Source changed before evaluation')
    for key in ('cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits'):
        ref = protocol[key]
        require(reference(ref['path'])['sha256'] == ref['sha256'], 'Evaluation input changed: ' + key)
    results = {spec['name']: json.loads((output / (spec['name'] + '_validation.json')).read_text()) for spec in protocol['recipes']}
    for name, result in results.items():
        require(result['protocol'] == reference(output / 'protocol.json'), 'Candidate protocol differs: ' + name)
        for key in ('model', 'history'):
            ref = result[key]
            require(reference(ref['path'])['sha256'] == ref['sha256'], 'Candidate artifact changed: ' + name + '/' + key)
    rank = lambda name: (results[name]['validation']['rejection'], results[name]['validation']['auc'], name)
    winner = max(results, key=rank)
    flatness = max(results, key=lambda name: (results[name]['validation']['selector_objective'], name))
    eligible = [name for name in results if envelope_eligible(results[name]['validation'], protocol['acceptance_envelope'])]
    acceptable = max(eligible, key=rank) if eligible else None
    selected = list(dict.fromkeys([winner, flatness] + ([acceptable] if acceptable else [])))
    write(output / 'validation_candidates.json', results)
    lock = {'schema': 'shift-vertex-dnn-selection-v1', 'selection_uses_test': False,
        'rejection_winner': winner, 'flatness_winner': flatness, 'acceptance_envelope_winner': acceptable,
        'acceptance_envelope_eligible': eligible, 'representative_models': selected,
        'thresholds': {name: entry['validation']['threshold'] for name, entry in results.items()},
        'protocol': reference(output / 'protocol.json'), 'validation': reference(output / 'validation_candidates.json'),
        'interpretation': 'Validation-selected research candidates. No final independent performance or lifetime-uniformity claim.'}
    write(output / 'selection_locked.json', lock)
    reco, labels, split_data = load(protocol['reco']['path']), load(protocol['labels']['path']), load(output / 'splits.npz')
    earlier = load(protocol['bdt_scores']['path'])
    score_table = {name: earlier[name] for name in ('ids', 'process', 'pair_index', 'source_index', 'heldout', 'expanded_split', 'expanded_heldout')}
    for name in ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400'):
        score_table[name] = earlier[name]
    previous_lock = json.loads((Path(protocol['bdt_scores']['path']).parent / 'selection_locked.json').read_text())
    for name in ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400'):
        lock['thresholds'][name] = previous_lock['thresholds'][name]
    # A separate combined diagnostic lock pins all reference cuts without modifying the model-selection lock.
    write(output / 'diagnostic_selection.json', lock)
    test = split_data['expanded_heldout'] & (labels['y'] >= 0)
    train = (split_data['expanded_split'] == 0) & (reco['process'] != 'signal') & (labels['y'] >= 0)
    test_results = {}
    for name in results:
        with Path(results[name]['model']['path']).open('rb') as stream:
            bundle = pickle.load(stream)
        require(tuple(bundle['feature_names']) == FEATURE_NAMES and np.array_equal(bundle['columns'], np.arange(len(FEATURE_NAMES))),
                'Candidate predictor contract changed')
        require(bundle['threshold'] == lock['thresholds'][name], 'Candidate threshold differs from validation lock')
        predicted = bundle['model'].predict_proba(reco['X'])[:, 1]
        require(predicted.shape == (len(reco['ids']),) and np.isfinite(predicted).all(), 'Invalid DNN scores')
        score_table[name] = predicted
        if name in selected:
            test_results[name] = summary(predicted[test], labels['y'][test], reco['process'][test], reco['nuisance'][test],
                                          reco['nuisance'][train], labels['y'][train], bundle['threshold'], reco['topology'][test])
    np.savez_compressed(output / 'development_scores.npz', **score_table)
    write(output / 'selected_development_metrics.json', test_results)
    diagnostic_inputs = {'cache_manifest': protocol['cache_manifest'], 'reco': protocol['reco'], 'labels': protocol['labels'],
        'scores': reference(output / 'development_scores.npz'), 'selection': reference(output / 'diagnostic_selection.json')}
    write(output / 'diagnostic_inputs.json', diagnostic_inputs)
    write(output / 'completion.json', {'complete': True, 'fit_count': len(results), 'protocol': reference(output / 'protocol.json'),
        'selection': reference(output / 'selection_locked.json'), 'scores': reference(output / 'development_scores.npz'),
        'diagnostic_inputs': reference(output / 'diagnostic_inputs.json'), 'selected_models': selected,
        'test_interpretation': 'Previously inspected development evidence; all model and threshold choices locked first',
        'physics_ready': False, 'production_deployed': False})
    print(json.dumps({'complete': True, 'selected_models': selected, 'acceptance_eligible': eligible}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--action', choices=('prepare', 'fit', 'finalize', 'all'), default='all')
    parser.add_argument('--recipe')
    args = parser.parse_args()
    workspace, output = args.workspace.resolve(), args.output.resolve()
    require(output.is_relative_to(workspace / 'tea_shift_cmssw/plots'), 'Final study results must be inside repository plots')
    if args.action in ('prepare', 'all'):
        prepare(workspace, output)
    if args.action in ('fit', 'all'):
        protocol = json.loads((output / 'protocol.json').read_text())
        specs = [spec for spec in protocol['recipes'] if args.recipe is None or spec['name'] == args.recipe]
        require(bool(specs), 'Unknown recipe')
        for spec in specs:
            fit_one(output, spec)
    if args.action in ('finalize', 'all'):
        finalize(output)


if __name__ == '__main__':
    main()
