#!/usr/bin/env python3
# Code style: AI-code. AI-generated reconstructed-information neural study.
"""Controlled feature/category/loss ablations, selected using seed pairs only."""
import argparse
import json
from pathlib import Path
import pickle
import shutil
import sys
import time

import numpy as np
from sklearn.metrics import roc_auc_score

from evaluation import event_split
from features import FEATURE_NAMES, nuisance_transform
from neural_vertex_study import envelope_eligible, load, reference, require, summary, write
from phase2_train import phase_weights


def recipes():
    # All comparisons are declared before fitting. One global score cut/model.
    families = [
        ('global_raw50', 'global', 'raw50', 3., 3., 0., 0.),
        ('global_rich', 'global', 'rich', 3., 3., 0., 0.),
        ('shared_raw50', 'shared', 'raw50', 1., 3., 1., 0.),
        ('shared_rich', 'shared', 'rich', 1., 3., 1., 0.),
        ('heads_raw50', 'heads', 'raw50', 1., 3., 1., 0.),
        ('heads_rich', 'heads', 'rich', 1., 3., 1., 0.),
        ('heads_rich_weak', 'heads', 'rich', .3, 1., 1., 0.),
        ('heads_rich_strong', 'heads', 'rich', 3., 3., 1., 0.),
        ('shared_rich_strong', 'shared', 'rich', 3., 3., 1., 0.),
        ('global_muon', 'global', 'muon', 3., 3., 0., 0.),
        ('shared_rich_accept10', 'shared', 'rich', .3, 1., 0., 10.),
        ('shared_rich_accept30', 'shared', 'rich', .3, 1., 0., 30.),
    ]
    return [dict(name=f'info_{family}_s{seed}', family=family, architecture=architecture,
                 inputs=inputs, strength=joint, category_strength=category,
                 conditional_strength=conditional, acceptance_strength=acceptance, seed=seed)
            for family, architecture, inputs, joint, category, conditional, acceptance in families
            for seed in (71, 72)]


def verify(ref):
    require(reference(ref['path'])['sha256'] == ref['sha256'], 'Changed sealed artifact: ' + ref['path'])


def input_columns(names, mode):
    require(tuple(names[:50]) == FEATURE_NAMES, 'Canonical X50 prefix changed')
    if mode == 'raw50':
        return np.arange(50)
    if mode == 'muon':
        return np.asarray([i for i, name in enumerate(names) if not str(name).startswith('event_')])
    require(mode == 'rich', 'Unknown input ablation')
    return np.arange(len(names))


def prepare(workspace, cache, output):
    require(not output.exists(), 'Refusing to overwrite a prepared study')
    manifest = json.loads((cache / 'inputs.json').read_text())
    require(manifest['complete'] and manifest['sample_kind'] == 'simulation_only', 'Incomplete rich cache')
    require(manifest['old_X50_exact'] and manifest['truth_branches_read'] is False, 'Unqualified reconstructed extension')
    for key in ('reco', 'labels', 'events'):
        require(reference(cache / (key + '.npz'))['sha256'] == manifest[key + '_sha256'], 'Changed extracted cache')
    reco, labels = load(cache / 'reco.npz'), load(cache / 'labels.npz')
    old = workspace / 'validation/bdt_discrimination_improvement_20261012/expanded_training'
    completion = json.loads((old / 'completion.json').read_text())
    verify(completion['development']['scores'])
    verify(completion['selection'])
    earlier = load(old / 'development_scores.npz')
    for key in ('ids', 'process', 'pair_index', 'source_index'):
        require(np.array_equal(reco[key], earlier[key]), 'Old physical pair identity changed: ' + key)
    require(reco['ids'].dtype == np.uint64, 'Physical identities must remain uint64')
    names = list(map(str, reco['feature_names']))
    require(len(names) == reco['X'].shape[1] and len(set(names)) == len(names), 'Feature names disagree')
    input_columns(names, 'raw50')
    extension_source = workspace / 'tea_shift_cmssw/utils/shift_vertex_reco_information.py'
    verify(manifest['protocol'])
    verify(manifest['source_snapshot'])
    require(reference(extension_source)['sha256'] == manifest['source_snapshot']['sha256'], 'Reviewed extractor source changed')
    sys.path.insert(0, str(extension_source.parent))
    from shift_vertex_reco_information import FEATURE_NAMES as EXTENSION_FEATURE_NAMES
    require(tuple(names) == FEATURE_NAMES + EXTENSION_FEATURE_NAMES, 'Unreviewed reconstructed feature contract')
    extraction_protocol = json.loads(Path(manifest['protocol']['path']).read_text())
    require(tuple(extraction_protocol['additional_feature_names']) == EXTENSION_FEATURE_NAMES
            and extraction_protocol['truth_branches_read'] is False, 'Feature extraction contract differs')
    old_manifest = json.loads(Path(extraction_protocol['old_manifest']['path']).read_text())
    verify(extraction_protocol['old_manifest'])
    old_reco_path = Path(extraction_protocol['old_manifest']['path']).parent / 'reco.npz'
    require(reference(old_reco_path)['sha256'] == old_manifest['reco_sha256'], 'Original predictor cache changed')
    old_reco = load(old_reco_path)
    require(np.array_equal(reco['X'][:, :50], old_reco['X'], equal_nan=True), 'Original predictor values changed')
    require(manifest['labels_sha256'] == old_manifest['labels_sha256'], 'Truth labels changed during feature extension')
    split = event_split(reco['ids'], reco['process'], 71)
    require(np.array_equal(split, earlier['expanded_split']), 'Physical-event split changed')
    sm = reco['process'] != 'signal'
    test = sm & (split == 2)
    require(np.array_equal(test, earlier['expanded_heldout']), 'Historical development-test mask changed')
    masks = [(split == part) & sm & (labels['y'] >= 0) for part in range(3)]
    require([int(mask.sum()) for mask in masks] == [9826, 3400, 3255], 'Controlled fitting population changed')
    keys = [{(str(p), *map(int, identity)) for p, identity in zip(reco['process'][mask], reco['ids'][mask])}
            for mask in masks]
    require(all(not keys[a].intersection(keys[b]) for a, b in ((0, 1), (0, 2), (1, 2))), 'Event crosses partitions')
    output.mkdir(parents=True)
    (output / 'models').mkdir()
    (output / 'source_snapshot').mkdir()
    sources = {}
    filenames = ('neural_vertex_information_study.py', 'neural_vertex_categories.py',
                 'vertex_acceptance_loss.py', 'neural_vertex_scheduled.py', 'neural_vertex.py',
                 'neural_vertex_study.py', 'disco.py', 'adversarial.py', 'features.py',
                 'evaluation.py', 'phase2_train.py', 'uniform_bdt.py', 'inference.py', 'bdt_performance.py')
    for filename in filenames:
        path = Path(__file__).parent / filename
        sources[filename] = reference(path)
        shutil.copy2(path, output / 'source_snapshot' / filename)
    sources[extension_source.name] = reference(extension_source)
    shutil.copy2(extension_source, output / 'source_snapshot' / extension_source.name)
    np.savez_compressed(output / 'splits.npz', ids=reco['ids'], process=reco['process'],
                        pair_index=reco['pair_index'], source_index=reco['source_index'],
                        expanded_split=split, expanded_heldout=test)
    protocol = dict(schema='shift-rich-vertex-neural-study-v1', recipes=recipes(), fit_count=len(recipes()),
                    cache_manifest=reference(cache / 'inputs.json'), reco=reference(cache / 'reco.npz'),
                    labels=reference(cache / 'labels.npz'), events=reference(cache / 'events.npz'),
                    bdt_scores=completion['development']['scores'], previous_bdt_selection=completion['selection'],
                    sources=sources, splits=reference(output / 'splits.npz'), input_feature_names=names,
                    training=dict(epochs=200, batch_size=512, patience=25, learning_rate=.0005,
                                  dropout=.1, weight_decay=.0005, min_delta=.00001, warmup_epochs=20,
                                  hidden=[64, 32], feature_transform='signedlog'),
                    acceptance_envelope=dict(minimum_supported_genuine_retention=.8,
                                             maximum_supported_marginal_span=.15,
                                             maximum_process_efficiency_span=.05,
                                             maximum_class_joint_dcor=.3,
                                             maximum_category_efficiency_span=.15),
                    weights='Identical historical equal-class, smoothed process-balanced phase_weights; soft acceptance is empirical/unweighted.',
                    nuisances='Only reconstructed mass/z/r, with training-only transforms and bins. No generated lifetime/mass targets.',
                    loss_ablation='Global models use the previous class/process joint penalty; category models use global class joint penalty plus explicitly weighted within-class/category and genuine category/process penalties. Each raw/rich input pair has identical loss and optimizer.',
                    selection='Group recipes by family; BOTH predeclared seeds must pass every acceptance gate. Rank eligible families by minimum seed rejection, then minimum AUC. No test-based fallback.',
                    thresholds='ONE inclusive empirical validation90percent genuine cut per model; no category-specific cuts or test calibration.',
                    support='20 genuine pairs per nuisance bin,50 per process/category in hard gates; sparse regions are unvalidated.',
                    known_counts=[int(mask.sum()) for mask in masks], excluded_signal_pairs=int((~sm).sum()),
                    input_ablation_labels=dict(raw50='Original reconstructed50 predictors',
                        rich='Original50 plus associated muon quality/hit/timing and event occupancy',
                        muon='Rich reconstructed information without event occupancy; original pair predictors retained'),
                    interpretation='Previously inspected SM test and existing signal points are consumed development evidence. This is an information ablation, not a data-volume increase or independent physics validation.')
    write(output / 'protocol.json', protocol)
    print(json.dumps({'prepared': True, 'fits': len(recipes()), 'features': len(names)}), flush=True)


def fitting_data(output):
    protocol = json.loads((output / 'protocol.json').read_text())
    for ref in protocol['sources'].values():
        verify(ref)
    for key in ('cache_manifest', 'reco', 'labels', 'events', 'bdt_scores', 'splits', 'previous_bdt_selection'):
        verify(protocol[key])
    reco, labels, splits = load(protocol['reco']['path']), load(protocol['labels']['path']), load(protocol['splits']['path'])
    for key in ('ids', 'process', 'pair_index', 'source_index'):
        require(np.array_equal(reco[key], splits[key]), 'Frozen split identity changed')
    require(list(map(str, reco['feature_names'])) == protocol['input_feature_names'], 'Feature schema changed')
    sm = reco['process'] != 'signal'
    masks = [(splits['expanded_split'] == part) & sm & (labels['y'] >= 0) for part in range(3)]
    return protocol, reco, labels, splits, masks


def predict(bundle, reco):
    require(list(map(str, reco['feature_names'])) == bundle['input_feature_names'], 'Prediction feature contract changed')
    columns = bundle['columns']
    require(tuple(reco['feature_names'][columns]) == tuple(bundle['feature_names']), 'Selected feature order changed')
    x = reco['X'][:, columns]
    if bundle['category_required']:
        return bundle['model'].predict_proba(x, reco['topology'])[:, 1]
    return bundle['model'].predict_proba(x)[:, 1]


def fit(output, spec):
    require(not (output / 'selection_locked.json').exists(), 'No fitting after selection')
    protocol, reco, labels, _, masks = fitting_data(output)
    require(spec in protocol['recipes'], 'Recipe was not declared')
    name = spec['name']
    receipt = output / (name + '_validation.json')
    require(not receipt.exists(), 'Refusing to replace an archived fit')
    tr, va, _ = masks
    names = protocol['input_feature_names']
    columns = input_columns(names, spec['inputs'])
    x, y, process, category, nuisance = reco['X'][:, columns], labels['y'], reco['process'], reco['topology'], reco['nuisance']
    z = nuisance_transform(nuisance)
    arrays = (x[tr], y[tr], phase_weights(y[tr], process[tr]), z[tr],
              x[va], y[va], phase_weights(y[va], process[va]), z[va])
    options = dict(protocol['training'], seed=spec['seed'], strength=spec['strength'],
                   category_strength=spec['category_strength'], train_process=process[tr],
                   val_process=process[va], train_category=category[tr], val_category=category[va])
    started = time.monotonic()
    category_required = spec['architecture'] != 'global'
    if category_required:
        from neural_vertex_categories import train_category_vertex
        model, history = train_category_vertex(*arrays, feature_names=[names[i] for i in columns],
            architecture=spec['architecture'], conditional_strength=spec['conditional_strength'],
            acceptance_strength=spec['acceptance_strength'], **options)
    else:
        from neural_vertex_scheduled import train_neural_vertex
        model, history = train_neural_vertex(*arrays, penalty_scope='class_process', energy_strength=0., **options)
    bundle = dict(model=model, columns=columns, feature_names=[names[i] for i in columns],
                  input_feature_names=names, category_required=category_required,
                  requires_gen=False, deployment_ready=False)
    val_reco = {key: value[va] if value.ndim and len(value) == len(y) else value for key, value in reco.items()}
    train_reco = {key: value[tr] if value.ndim and len(value) == len(y) else value for key, value in reco.items()}
    validation = summary(predict(bundle, val_reco), y[va], process[va], nuisance[va], nuisance[tr], y[tr], category=category[va])
    bundle['threshold'] = validation['threshold']
    train_score = predict(bundle, train_reco)
    training = dict(auc=float(roc_auc_score(y[tr], train_score)),
                    genuine_retention=float(np.mean(train_score[y[tr] == 1] >= bundle['threshold'])),
                    accidental_rejection=float(np.mean(train_score[y[tr] == 0] < bundle['threshold'])))
    with (output / 'models' / (name + '.pkl')).open('xb') as stream:
        pickle.dump(bundle, stream)
    write(output / (name + '_history.json'), history)
    result = dict(recipe=spec, fit_seconds=time.monotonic() - started, validation=validation,
                  training_at_validation_cut=training,
                  training_minus_validation_auc=training['auc'] - validation['auc'],
                  model=reference(output / 'models' / (name + '.pkl')),
                  history=reference(output / (name + '_history.json')), protocol=reference(output / 'protocol.json'))
    write(receipt, result)
    print(json.dumps(dict(model=name, rejection=validation['rejection'], auc=validation['auc'],
                          marginal_span=validation['maximum_marginal_efficiency_span'],
                          category_span=validation['category_efficiency_span'],
                          eligible=envelope_eligible(validation, protocol['acceptance_envelope']),
                          best_epoch=history['best_epoch'], seconds=result['fit_seconds'])), flush=True)


def seal(output):
    require(not (output / 'selection_locked.json').exists(), 'Refusing to overwrite selection')
    require(not (output / 'development_scores.npz').exists(), 'Selection must precede new scores')
    protocol, _, _, _, _ = fitting_data(output)
    results = {}
    families = {}
    for spec in protocol['recipes']:
        name = spec['name']
        result = json.loads((output / (name + '_validation.json')).read_text())
        require(result['recipe'] == spec and result['protocol'] == reference(output / 'protocol.json'), 'Recipe receipt changed')
        for key in ('model', 'history'):
            verify(result[key])
        history = json.loads(Path(result['history']['path']).read_text())
        require(history['restored_prediction_exact'] and history['best_state_sha256'] == history['restored_state_sha256'], 'Checkpoint was not restored')
        require(history['best_epoch'] == min(history['epochs'], key=lambda row: row['val_objective'])['epoch'], 'Wrong checkpoint objective')
        results[name] = result
        families.setdefault(spec['family'], []).append(name)
    eligible = []
    family_report = {}
    for family, names in families.items():
        require({results[name]['recipe']['seed'] for name in names} == {71, 72}, 'Unpaired recipe family')
        passing = [envelope_eligible(results[name]['validation'], protocol['acceptance_envelope']) for name in names]
        entry = dict(models=names, all_seeds_pass=all(passing), passing=passing,
                     minimum_rejection=min(results[name]['validation']['rejection'] for name in names),
                     minimum_auc=min(results[name]['validation']['auc'] for name in names))
        family_report[family] = entry
        if entry['all_seeds_pass']:
            eligible.append(family)
    robust = max(eligible, key=lambda family: (family_report[family]['minimum_rejection'], family_report[family]['minimum_auc'], family)) if eligible else None
    rejection = max(results, key=lambda name: (results[name]['validation']['rejection'], results[name]['validation']['auc'], name))
    flatness = max(results, key=lambda name: (results[name]['validation']['selector_objective'], name))
    selected = list(dict.fromkeys(([name for name in family_report[robust]['models']] if robust else []) + [rejection, flatness]))
    paired = []
    for name in selected:
        family = results[name]['recipe']['family']
        paired.extend(other for other in family_report[family]['models'] if other not in selected and other not in paired)
    write(output / 'validation_candidates.json', results)
    write(output / 'family_validation.json', family_report)
    lock = dict(schema='shift-rich-neural-family-selection-v1', selection_uses_test=False,
                rejection_winner=rejection, flatness_winner=flatness,
                acceptance_envelope_winner=family_report[robust]['models'][0] if robust else None,
                robust_family=robust, robust_eligible_families=eligible,
                representative_models=selected, paired_seed_models=paired,
                thresholds={name: result['validation']['threshold'] for name, result in results.items()},
                protocol=reference(output / 'protocol.json'), validation=reference(output / 'validation_candidates.json'),
                family_validation=reference(output / 'family_validation.json'), physics_ready=False)
    write(output / 'selection_locked.json', lock)
    diagnostic = dict(lock, thresholds=dict(lock['thresholds']))
    old_lock = json.loads(Path(protocol['previous_bdt_selection']['path']).read_text())
    for name in ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400'):
        diagnostic['thresholds'][name] = old_lock['thresholds'][name]
    write(output / 'diagnostic_selection.json', diagnostic)
    write(output / 'pre_evaluation_receipt.json', dict(protocol=reference(output / 'protocol.json'),
          selection=reference(output / 'selection_locked.json'), diagnostic_selection=reference(output / 'diagnostic_selection.json'),
          validation=reference(output / 'validation_candidates.json'), family_validation=reference(output / 'family_validation.json'),
          new_test_scores_generated=False))
    print(json.dumps(dict(locked=True, robust_family=robust, eligible_families=eligible, selected=selected, paired=paired)), flush=True)


def evaluate(output):
    for name in ('development_scores.npz', 'diagnostic_inputs.json', 'selected_development_metrics.json', 'reload_parity.json', 'completion.json'):
        require(not (output / name).exists(), 'Refusing to overwrite evaluation: ' + name)
    before = json.loads((output / 'pre_evaluation_receipt.json').read_text())
    for key in ('protocol', 'selection', 'diagnostic_selection', 'validation', 'family_validation'):
        verify(before[key])
    protocol, reco, labels, splits, masks = fitting_data(output)
    lock = json.loads((output / 'selection_locked.json').read_text())
    diagnostic = json.loads((output / 'diagnostic_selection.json').read_text())
    require(lock['selection_uses_test'] is False, 'Unsafe selection')
    for name, cut in lock['thresholds'].items():
        require(diagnostic['thresholds'][name] == cut, 'Neural cut changed')
    old_lock = json.loads(Path(protocol['previous_bdt_selection']['path']).read_text())
    for name in ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400'):
        require(diagnostic['thresholds'][name] == old_lock['thresholds'][name], 'Reference cut changed')
    earlier = load(protocol['bdt_scores']['path'])
    scores = {key: earlier[key] for key in ('ids', 'process', 'pair_index', 'source_index', 'heldout', 'expanded_split', 'expanded_heldout')}
    for name in ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400'):
        scores[name] = earlier[name]
    results = json.loads((output / 'validation_candidates.json').read_text())
    parity = {}
    for name in lock['representative_models'] + lock['paired_seed_models']:
        result = results[name]
        for key in ('model', 'history', 'protocol'):
            verify(result[key])
        with Path(result['model']['path']).open('rb') as stream:
            bundle = pickle.load(stream)
        require(bundle['requires_gen'] is False and bundle['threshold'] == lock['thresholds'][name], 'Inference contract or cut changed')
        predicted = predict(bundle, reco)
        with Path(result['model']['path']).open('rb') as stream:
            reloaded = pickle.load(stream)
        require(np.array_equal(predicted, predict(reloaded, reco)), 'Trusted reload score parity failed')
        require(predicted.shape == labels['y'].shape and np.isfinite(predicted).all(), 'Invalid model predictions')
        scores[name] = predicted
        parity[name] = dict(rows=len(predicted), bitwise_reload_parity=True)
    np.savez_compressed(output / 'development_scores.npz', **scores)
    inputs = {key: protocol[key] for key in ('cache_manifest', 'reco', 'labels')}
    inputs.update(scores=reference(output / 'development_scores.npz'), selection=reference(output / 'diagnostic_selection.json'))
    write(output / 'diagnostic_inputs.json', inputs)
    tr, _, te = masks
    metrics = {name: summary(values[te], labels['y'][te], reco['process'][te], reco['nuisance'][te],
                            reco['nuisance'][tr], labels['y'][tr], diagnostic['thresholds'][name], reco['topology'][te])
               for name, values in scores.items() if name in diagnostic['thresholds']}
    write(output / 'selected_development_metrics.json', metrics)
    write(output / 'reload_parity.json', parity)
    write(output / 'completion.json', dict(complete=True, protocol=reference(output / 'protocol.json'),
        selection=reference(output / 'selection_locked.json'), scores=reference(output / 'development_scores.npz'),
        diagnostic_inputs=reference(output / 'diagnostic_inputs.json'), metrics=reference(output / 'selected_development_metrics.json'),
        parity=reference(output / 'reload_parity.json'), fit_count=protocol['fit_count'],
        physics_ready=False, production_deployed=False))
    print(json.dumps(dict(evaluated=True, selected=lock['representative_models'], robust_family=lock['robust_family'])), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--action', choices=('prepare', 'fit', 'lock', 'evaluate'), required=True)
    parser.add_argument('--recipe')
    parser.add_argument('--shard', type=int, default=0)
    parser.add_argument('--shards', type=int, default=1)
    args = parser.parse_args()
    workspace, output = args.workspace.resolve(), args.output.resolve()
    require(output.is_relative_to(workspace / 'tea_shift_cmssw/plots'), 'Results must be inside repository plots')
    require(0 <= args.shard < args.shards, 'Invalid fit shard')
    if args.action == 'prepare':
        require(args.cache is not None, 'Explicit rich cache required')
        prepare(workspace, args.cache.resolve(), output)
    elif args.action == 'fit':
        protocol = json.loads((output / 'protocol.json').read_text())
        chosen = [spec for index, spec in enumerate(protocol['recipes']) if index % args.shards == args.shard
                  and (args.recipe is None or spec['name'] == args.recipe)]
        require(bool(chosen), 'No matching declared recipes')
        for spec in chosen:
            if not (output / (spec['name'] + '_validation.json')).exists():
                fit(output, spec)
    elif args.action == 'lock':
        seal(output)
    else:
        evaluate(output)


if __name__ == '__main__':
    main()
