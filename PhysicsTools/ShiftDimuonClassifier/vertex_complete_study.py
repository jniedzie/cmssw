#!/usr/bin/env python3
# Code style: AI-code. AI-generated complete-input BDT/DNN learning study.
"""All models receive RECO119 and category; event-grouped nested sample sizes.

Current and historical samples have separate protocols, fits and validation.
There is no pooling, collision-data access, new cut in production, or claim
that an already inspected development partition is independent confirmation.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import pickle
import shutil
import sys
import time

import numpy as np
from sklearn.metrics import roc_auc_score

from evaluation import event_split, fraction
from features import nuisance_transform
from neural_vertex_study import summary, envelope_eligible, reference, write, load, require
from phase2_train import phase_weights
from vertex_complete_bdt import COMPLETE_FEATURE_NAMES, CATEGORY_NAMES, fit_complete_bdt


FRACTIONS = (.25, .5, 1.)
SEEDS = (71, 72, 73)
PROTOCOL_KIND = 'capacity_and_strong_decorrelation'
BOUNDS = {'minimum_supported_genuine_retention': .8, 'maximum_supported_marginal_span': .15,
          'maximum_process_efficiency_span': .05, 'maximum_class_joint_dcor': .3,
          'maximum_category_efficiency_span': .15}
SOURCES = ('vertex_complete_study.py', 'vertex_complete_bdt.py', 'vertex_complete_neural.py',
           'neural_vertex_categories.py', 'neural_vertex.py', 'vertex_acceptance_loss.py',
           'neural_vertex_study.py', 'features.py', 'evaluation.py', 'phase2_train.py',
           'adversarial.py', 'disco.py', 'uniform_bdt.py', 'inference.py')


def verify(ref):
    require(reference(ref['path'])['sha256'] == ref['sha256'], 'Changed artifact: ' + ref['path'])


def recipes():
    families = [dict(family='bdt_' + config, algorithm='bdt', config=config)
                for config in ('shallow_d3', 'medium_d5', 'regularized_d5')]
    for architecture, hidden, label in [('shared', [64, 32], 'shared64'),
                                         ('shared', [128, 64], 'shared128'),
                                         ('heads', [64, 32], 'heads64')]:
        for loss, joint, category, conditional in [('classification', 0., 0., 0.),
                                                    ('decorrelated', 3., 3., 1.)]:
            families.append(dict(family='dnn_' + label + '_' + loss, algorithm='dnn',
                architecture=architecture, hidden=hidden, loss=loss,
                strength=joint, category_strength=category, conditional_strength=conditional))
    return [{**family, 'fraction': size, 'seed': seed,
             'name': f"{family['family']}_f{int(size*100)}_s{seed}"}
            for family in families for size in FRACTIONS for seed in SEEDS]


def parent_keys(reco, domain):
    if domain == 'historical_proxy':
        keys = np.asarray(reco['physical_parent_key'])
        require(keys.shape == (len(reco['ids']),), 'Historical physical-parent context is missing')
        return keys
    require(domain == 'current_CMSIR5', 'Unknown physics domain')
    return np.asarray([str(process) + ':' + ':'.join(map(str, identity))
                       for process, identity in zip(reco['process'], reco['ids'])])


def hash_values(keys, seed, purpose):
    return np.asarray([int.from_bytes(hashlib.sha256(
        f'vertex-complete:{purpose}:{seed}:{key}'.encode()).digest()[:8], 'big') / 2**64
        for key in keys])


def nested_subset(keys, train, seed, size):
    require(size in FRACTIONS, 'Undeclared training fraction')
    values = hash_values(keys, seed, 'nested-training')
    return np.asarray(train, dtype=bool) & (values < size)


def counts(reco, labels, mask):
    y = labels['y'][mask]
    process, category = reco['process'][mask], reco['topology'][mask]
    weights = phase_weights(y, process) if len(y) and len(np.unique(y)) == 2 else np.zeros(len(y))
    return {'pairs': int(mask.sum()), 'genuine': int((y == 1).sum()), 'fake': int((y == 0).sum()),
        'unique_parents': len(set(map(str, reco['_parent_key'][mask]))),
        'categories': {str(cat): {'genuine': int(((category == cat) & (y == 1)).sum()),
                                  'fake': int(((category == cat) & (y == 0)).sum())} for cat in CATEGORY_NAMES},
        'processes': {str(p): {'genuine': int(((process == p) & (y == 1)).sum()),
                               'fake': int(((process == p) & (y == 0)).sum())} for p in np.unique(process)},
        'training_weight_masses': {f'{p}:{label}': float(weights[(process == p) & (y == label)].sum())
                                  for p in np.unique(process) for label in (0, 1)}}


def prepare(cache, output, domain):
    require(not output.exists(), 'Refusing study overwrite')
    manifest_path = cache / 'inputs.json'
    manifest = json.loads(manifest_path.read_text())
    require(manifest['complete'] and manifest['sample_kind'] == 'simulation_only', 'Incomplete/non-MC cache')
    refs = {name: reference(cache / (name + '.npz')) for name in ('reco', 'labels', 'events')}
    for name, ref in refs.items():
        expected = manifest.get(name + '_sha256', manifest.get('cache', {}).get(name, {}).get('sha256'))
        require(ref['sha256'] == expected, 'Cache digest mismatch: ' + name)
    reco, labels = load(refs['reco']['path']), load(refs['labels']['path'])
    require(tuple(reco['feature_names']) == COMPLETE_FEATURE_NAMES and reco['X'].shape[1] == 119,
            'Every fit requires all 119 reconstructed features')
    require(set(reco['topology']).issubset(CATEGORY_NAMES), 'Unknown reconstruction category')
    require(reco['ids'].dtype == np.uint64, 'Rounded physical event identities')
    keys = parent_keys(reco, domain)
    sm = np.isin(reco['process'], ('dy', 'jpsi', 'qcd'))
    if domain == 'current_CMSIR5':
        split = event_split(reco['ids'], reco['process'], 71)
    else:
        values = hash_values(keys, 71, 'historical-partition')
        split = np.where(values < .6, 0, np.where(values < .8, 1, 2)).astype(np.int8)
    masks = [(split == part) & sm & (labels['y'] >= 0) for part in range(3)]
    event_sets = [set(keys[mask]) for mask in masks]
    require(all(not event_sets[a] & event_sets[b] for a, b in ((0, 1), (0, 2), (1, 2))),
            'Physical parents cross train/validation/development partitions')
    if domain == 'current_CMSIR5':
        require([int(mask.sum()) for mask in masks] == [9826, 3400, 3255], 'Current full training population changed')
    for mask in masks:
        require(set(labels['y'][mask]) == {0, 1}, 'Both classes required in each known-label partition')
    reco['_parent_key'] = keys
    subsets = {}
    for seed in SEEDS:
        previous = np.zeros(len(keys), dtype=bool)
        for size in FRACTIONS:
            mask = nested_subset(keys, masks[0], seed, size)
            require(not np.any(previous & ~mask), 'Training subsets are not nested')
            require(not np.any(mask & ~masks[0]), 'Subset leaks validation/test/signal/unknown rows')
            require(set(labels['y'][mask]) == {0, 1}, 'Subset lacks a positive-weight class')
            subsets[f'f{int(size*100)}_s{seed}'] = mask
            previous = mask
    output.mkdir(parents=True, exist_ok=False)
    for directory in ('models', 'fits', 'histories', 'source_snapshot'):
        (output / directory).mkdir()
    sources = {}
    for name in SOURCES:
        original = Path(__file__).parent / name
        copied = output / 'source_snapshot' / name
        shutil.copyfile(original, copied)
        sources[name] = {**reference(copied), 'original_path': str(original.resolve())}
    np.savez_compressed(output / 'splits.npz', ids=reco['ids'], process=reco['process'],
        pair_index=reco['pair_index'], source_index=reco['source_index'], parent_key=keys, split=split, **subsets)
    protocol = {'schema': 'shift-complete-vertex-learning-study-v1', 'domain': domain,
        'setting_menu': PROTOCOL_KIND,
        'recipes': recipes(), 'fit_budget': 81, 'fractions': FRACTIONS, 'seeds': SEEDS,
        'cache_manifest': reference(manifest_path), **refs, 'splits': reference(output / 'splits.npz'),
        'sources': sources, 'input_feature_names': list(COMPLETE_FEATURE_NAMES),
        'category_names': list(CATEGORY_NAMES), 'category_required_for_every_model': True,
        'partitions': {name: counts(reco, labels, mask) for name, mask in zip(('train', 'validation', 'development'), masks)},
        'subsets': {name: counts(reco, labels, mask) for name, mask in subsets.items()},
        'training': {'epochs_at_full_size': 240, 'batch_size': 512, 'patience': 30, 'warmup_epochs': 20,
            'dropout': .1, 'learning_rate': .0005, 'weight_decay': .0005, 'min_delta': 1e-5},
        'epoch_budget': 'ceil(240/fraction) maximum; early stopping normally limits fits; optimizer updates and cap exhaustion reported',
        'training_weights': 'Subset-only equal class/process weights with historical count50 smoothing; masses recorded per subset',
        'reporting_bins': 'Fixed full TRAIN genuine reconstructed mass/z/r tertiles; same edges across all sample sizes',
        'acceptance_envelope': BOUNDS,
        'threshold': 'One empirical validation genuine10th percentile, inclusive >=; all ties/counts reported',
        'selection': 'Full-size validation maximin rejection per algorithm; acceptance ranking additionally requires every seed passing fixed gates',
        'test_policy': 'No development predictions before all81 validation receipts and locked family choices',
        'signal_training_allowed': False, 'requires_gen_for_inference': False,
        'physics_ready': False, 'interpretation': 'Development comparison; existing current validation/test already inspected. Historical proxy kept separate.'}
    write(output / 'protocol.json', protocol)
    print(json.dumps({'prepared': str(output), 'domain': domain, 'fits': 81,
                      'train': int(masks[0].sum()), 'validation': int(masks[1].sum())}), flush=True)


def fitting_data(output):
    protocol = json.loads((output / 'protocol.json').read_text())
    for key in ('reco', 'labels', 'events', 'cache_manifest', 'splits'):
        verify(protocol[key])
    for ref in protocol['sources'].values():
        verify(ref)
    require(Path(__file__).resolve() == Path(protocol['sources']['vertex_complete_study.py']['path']).resolve(),
            'Execute the frozen study driver from its source_snapshot')
    for filename, ref in protocol['sources'].items():
        module = sys.modules.get(Path(filename).stem)
        if module is not None:
            require(Path(module.__file__).resolve() == Path(ref['path']).resolve(),
                    'Imported source lies outside the sealed snapshot: ' + filename)
    reco, labels, splits = load(protocol['reco']['path']), load(protocol['labels']['path']), load(protocol['splits']['path'])
    for key in ('ids', 'process', 'pair_index', 'source_index'):
        require(np.array_equal(reco[key], splits[key]), 'Split/candidate alignment changed: ' + key)
    reco['_parent_key'] = splits['parent_key']
    return protocol, reco, labels, splits


def model_score(bundle, reco, mask):
    # The only classifier arguments are reconstructed X and topology.
    return bundle['model'].predict_proba(reco['X'][mask], reco['topology'][mask])[:, 1]


def fit(output, spec):
    require(not (output / 'selection_locked.json').exists(), 'Fitting is closed after selection lock')
    path = output / 'fits' / (spec['name'] + '.json')
    require(not path.exists(), 'Refusing archived fit overwrite')
    protocol, reco, labels, splits = fitting_data(output)
    require(spec in protocol['recipes'], 'Undeclared recipe')
    sm = np.isin(reco['process'], ('dy', 'jpsi', 'qcd')) & (labels['y'] >= 0)
    full = (splits['split'] == 0) & sm
    train = splits[f"f{int(spec['fraction']*100)}_s{spec['seed']}"]
    val = (splits['split'] == 1) & sm
    x, y, process, category, nuisance = reco['X'], labels['y'], reco['process'], reco['topology'], reco['nuisance']
    started = time.monotonic()
    weights = phase_weights(y[train], process[train])
    if spec['algorithm'] == 'bdt':
        model = fit_complete_bdt(x[train], category[train], y[train], weights,
            feature_names=COMPLETE_FEATURE_NAMES, seed=spec['seed'], config=spec['config'])
        history = model.training_metadata_
        optimizer = {'epochs_completed': None, 'best_epoch': None, 'updates_per_epoch': None,
                     'total_updates': 200, 'cap_reached': False, 'stopped_early': False,
                     'fixed_tree_iterations': 200}
    else:
        from vertex_complete_neural import train_complete_neural
        options = dict(protocol['training'])
        options['epochs'] = int(np.ceil(options.pop('epochs_at_full_size') / spec['fraction']))
        model, history = train_complete_neural(x[train], y[train], weights, nuisance_transform(nuisance[train]),
            x[val], y[val], phase_weights(y[val], process[val]), nuisance_transform(nuisance[val]),
            feature_names=COMPLETE_FEATURE_NAMES, train_category=category[train], val_category=category[val],
            train_process=process[train], val_process=process[val], seed=spec['seed'],
            architecture=spec['architecture'], hidden=spec['hidden'], strength=spec['strength'],
            category_strength=spec['category_strength'], conditional_strength=spec['conditional_strength'],
            acceptance_strength=spec.get('acceptance_strength', 0.), feature_transform='signedlog', **options)
        optimizer = {key: history[key] for key in ('epochs_completed', 'best_epoch', 'updates_per_epoch')}
        optimizer.update(total_updates=history['epochs_completed'] * history['updates_per_epoch'],
            cap_reached=history['epochs_completed'] == history['epochs_requested'],
            stopped_early=history['epochs_completed'] < history['epochs_requested'])
    bundle = {'model': model, 'feature_names': list(COMPLETE_FEATURE_NAMES), 'category_required': True,
              'recipe': spec, 'requires_gen': False, 'deployment_ready': False}
    val_score, train_score = model_score(bundle, reco, val), model_score(bundle, reco, train)
    validation = summary(val_score, y[val], process[val], nuisance[val], nuisance[full], y[full], category=category[val])
    bundle['threshold'] = validation['threshold']
    for label, label_name in ((1, 'genuine'), (0, 'fake')):
        validation['exposure_weighted_' + label_name] = fraction(y[val] == label, val_score >= bundle['threshold'],
                                                                  reco['evaluation_weight'][val])
    model_path = output / 'models' / (spec['name'] + '.pkl')
    with model_path.open('xb') as stream:
        pickle.dump(bundle, stream)
    with model_path.open('rb') as stream:
        restored = pickle.load(stream)
    require(np.array_equal(val_score, model_score(restored, reco, val)), 'Model reload parity failed')
    history_path = output / 'histories' / (spec['name'] + '.json')
    write(history_path, history)
    receipt = {'recipe': spec, 'training_counts': protocol['subsets'][f"f{int(spec['fraction']*100)}_s{spec['seed']}"],
        'validation': validation, 'training_at_validation_cut': {'auc': float(roc_auc_score(y[train], train_score)),
            'genuine_retention': float(np.mean(train_score[y[train] == 1] >= bundle['threshold'])),
            'accidental_rejection': float(np.mean(train_score[y[train] == 0] < bundle['threshold']))},
        'fit_seconds': time.monotonic() - started, 'optimizer': optimizer,
        'model': reference(model_path), 'history': reference(history_path), 'protocol': reference(output / 'protocol.json'),
        'reload_parity_exact': True}
    write(path, receipt)
    print(json.dumps({'fit': spec['name'], 'domain': protocol['domain'], 'seconds': round(receipt['fit_seconds'], 2),
        'rejection': validation['rejection'], 'auc': validation['auc'], 'genuine': validation['common_vertex']['efficiency'],
        'marginal_span': validation['maximum_marginal_efficiency_span'],
        'category_span': validation['category_efficiency_span']}), flush=True)


def complete_fits(output):
    protocol = json.loads((output / 'protocol.json').read_text())
    require(protocol['recipes'] == recipes(), 'Declared recipe grid differs from the reviewed 81 fits')
    this_protocol = reference(output / 'protocol.json')
    fits = {}
    for spec in protocol['recipes']:
        path = output / 'fits' / (spec['name'] + '.json')
        require(path.exists(), 'Missing fit before selection lock: ' + spec['name'])
        receipt = json.loads(path.read_text())
        require(receipt['recipe'] == spec, 'Fit recipe changed')
        require(receipt['protocol'] == this_protocol, 'Fit belongs to another study/domain')
        verify(receipt['protocol']); verify(receipt['model']); verify(receipt['history'])
        fits[spec['name']] = receipt
    return protocol, fits


def lock(output):
    protocol, fits = complete_fits(output)
    families = {}
    for spec in protocol['recipes']:
        if spec['fraction'] != 1.:
            continue
        entry = fits[spec['name']]['validation']
        row = families.setdefault(spec['family'], {'algorithm': spec['algorithm'], 'seeds': {}})
        row['seeds'][str(spec['seed'])] = {'name': spec['name'], 'rejection': entry['rejection'],
            'auc': entry['auc'], 'passed_acceptance':
                all(entry['dependence'][label]['joint_dcor'] is not None
                    for label in ('common_vertex', 'accidental')) and envelope_eligible(entry, BOUNDS)}
    for row in families.values():
        require(set(row['seeds']) == set(map(str, SEEDS)), 'Family missing a declared seed')
        row['minimum_rejection'] = min(seed['rejection'] for seed in row['seeds'].values())
        row['minimum_auc'] = min(seed['auc'] for seed in row['seeds'].values())
        row['all_seeds_pass_acceptance'] = all(seed['passed_acceptance'] for seed in row['seeds'].values())
    winners = {}
    for algorithm in ('bdt', 'dnn'):
        available = [name for name, row in families.items() if row['algorithm'] == algorithm]
        eligible = [name for name in available if families[name]['all_seeds_pass_acceptance']]
        rank = lambda name: (families[name]['minimum_rejection'], families[name]['minimum_auc'], name)
        winners[algorithm] = {'rejection': max(available, key=rank),
                              'acceptance': max(eligible, key=rank) if eligible else None}
    locked = {'schema': 'shift-complete-vertex-family-lock-v1', 'selection_uses_test': False,
        'domain': protocol['domain'], 'winners': winners, 'families': families,
        'protocol': reference(output / 'protocol.json'),
        'validation_fits': {name: reference(output / 'fits' / (name + '.json')) for name in fits},
        'thresholds': {name: fit['validation']['threshold'] for name, fit in fits.items()},
        'all81fits_before_development': True, 'requires_gen_for_inference': False, 'physics_ready': False}
    write(output / 'selection_locked.json', locked)
    print(json.dumps({'locked': protocol['domain'], 'winners': winners}), flush=True)


def evaluate(output):
    protocol, reco, labels, splits = fitting_data(output)
    locked = json.loads((output / 'selection_locked.json').read_text()); verify(locked['protocol'])
    require(locked['protocol'] == reference(output / 'protocol.json') and locked['domain'] == protocol['domain'],
            'Selection belongs to another study/domain')
    require(set(locked['validation_fits']) == {spec['name'] for spec in protocol['recipes']},
            'Selection did not lock exactly the declared 81 fits')
    for ref in locked['validation_fits'].values():
        verify(ref)
    sm = np.isin(reco['process'], ('dy', 'jpsi', 'qcd')) & (labels['y'] >= 0)
    test, full = (splits['split'] == 2) & sm, (splits['split'] == 0) & sm
    families = sorted({family for values in locked['winners'].values() for family in values.values() if family})
    selected = [spec for spec in protocol['recipes'] if spec['family'] in families]
    results, predicted = {}, {}
    for spec in selected:
        receipt = json.loads((output / 'fits' / (spec['name'] + '.json')).read_text()); verify(receipt['model'])
        require(receipt['protocol'] == reference(output / 'protocol.json') and receipt['recipe'] == spec,
                'Development fit belongs to another study/domain')
        require(receipt['validation']['threshold'] == locked['thresholds'][spec['name']],
                'Locked threshold differs from validation fit')
        with Path(receipt['model']['path']).open('rb') as stream:
            bundle = pickle.load(stream)
        require(bundle['requires_gen'] is False and bundle['threshold'] == locked['thresholds'][spec['name']], 'Inference cut/contract changed')
        score = model_score(bundle, reco, test)
        results[spec['name']] = summary(score, labels['y'][test], reco['process'][test], reco['nuisance'][test],
            reco['nuisance'][full], labels['y'][full], threshold=bundle['threshold'], category=reco['topology'][test])
        predicted[spec['name']] = score
    np.savez_compressed(output / 'development_scores.npz', ids=reco['ids'][test], process=reco['process'][test],
        category=reco['topology'][test], nuisance=reco['nuisance'][test], y=labels['y'][test],
        evaluation_weight=reco['evaluation_weight'][test], **predicted)
    write(output / 'development_metrics.json', {'metrics': results, 'selection': reference(output / 'selection_locked.json'),
        'domain': protocol['domain'], 'selection_uses_test': False, 'previously_inspected_development': protocol['domain'] == 'current_CMSIR5'})
    write(output / 'completion.json', {'complete': True, 'domain': protocol['domain'], 'fits': 81,
        'selection': reference(output / 'selection_locked.json'), 'protocol': reference(output / 'protocol.json'),
        'scores': reference(output / 'development_scores.npz'), 'metrics': reference(output / 'development_metrics.json'),
        'physics_ready': False, 'production_deployed': False})
    print(json.dumps({'evaluated': protocol['domain'], 'models': len(selected), 'pairs': int(test.sum())}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--action', choices=('prepare', 'fit', 'lock', 'evaluate'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--domain', choices=('current_CMSIR5', 'historical_proxy'), default='current_CMSIR5')
    parser.add_argument('--shard', type=int, default=0)
    parser.add_argument('--shards', type=int, default=1)
    args = parser.parse_args()
    require(0 <= args.shard < args.shards, 'Invalid shard partition')
    if args.action == 'prepare':
        require(args.cache is not None, 'Preparation needs an explicit complete cache')
        prepare(args.cache.resolve(), args.output.resolve(), args.domain)
    elif args.action == 'fit':
        protocol = json.loads((args.output / 'protocol.json').read_text())
        for index, spec in enumerate(protocol['recipes']):
            if index % args.shards == args.shard and not (args.output / 'fits' / (spec['name'] + '.json')).exists():
                fit(args.output, spec)
    elif args.action == 'lock':
        lock(args.output)
    else:
        evaluate(args.output)


if __name__ == '__main__':
    main()
