#!/usr/bin/env python3
# Code style: AI-code. AI-generated validation lock and selected-model evaluation.
"""Seal cross-round choices before evaluating any selected neural test scores."""
import argparse
import json
from pathlib import Path
import pickle
import shutil

import numpy as np

from features import FEATURE_NAMES
from neural_vertex_study import envelope_eligible, load, reference, require, summary, write


def verify_reference(ref):
    require(reference(ref['path'])['sha256'] == ref['sha256'], 'Artifact changed: ' + ref['path'])


def gather(rounds):
    results, protocols = {}, []
    for directory in rounds:
        require(not (directory / 'development_scores.npz').exists(), 'Global choices must precede neural test-score generation')
        protocol = json.loads((directory / 'protocol.json').read_text())
        protocols.append((directory, protocol))
        for ref in protocol['sources'].values():
            verify_reference(ref)
            require(reference(directory / 'source_snapshot' / Path(ref['path']).name)['sha256'] == ref['sha256'], 'Source snapshot changed')
        for key in ('cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits'):
            verify_reference(protocol[key])
        for spec in protocol['recipes']:
            name = spec['name']
            require(name not in results, 'Duplicate recipe name across rounds')
            entry = json.loads((directory / (name + '_validation.json')).read_text())
            require(entry['recipe'] == spec and entry['protocol'] == reference(directory / 'protocol.json'), 'Recipe origin changed')
            for key in ('model', 'history'):
                verify_reference(entry[key])
            history = json.loads(Path(entry['history']['path']).read_text())
            require(history['restored_prediction_exact'] and history['best_state_sha256'] == history['restored_state_sha256'], 'Checkpoint restoration failed')
            require(history['best_epoch'] == min(history['epochs'], key=lambda row: row['val_objective'])['epoch'], 'Checkpoint is not validation minimum')
            results[name] = entry
    return results, protocols


def seal(output, rounds):
    require(not output.exists(), 'Refusing to overwrite a sealed comparison')
    results, protocols = gather(rounds)
    baseline = protocols[0][1]
    for _, protocol in protocols:
        for key in ('cache_manifest', 'reco', 'labels', 'bdt_scores'):
            require(protocol[key] == baseline[key], 'Different fitting populations between rounds')
        require(protocol['acceptance_envelope'] == baseline['acceptance_envelope'], 'Acceptance target changed across rounds')
        require(reference(protocol['splits']['path'])['sha256'] == baseline['splits']['sha256'], 'Event splits changed across rounds')
    rejection = max(results, key=lambda name: (results[name]['validation']['rejection'], results[name]['validation']['auc'], name))
    flatness = max(results, key=lambda name: (results[name]['validation']['selector_objective'], name))
    eligible = [name for name in results if envelope_eligible(results[name]['validation'], baseline['acceptance_envelope'])]
    acceptable = max(eligible, key=lambda name: (results[name]['validation']['rejection'], results[name]['validation']['auc'], name)) if eligible else None
    selected = list(dict.fromkeys([rejection, flatness] + ([acceptable] if acceptable else [])))
    paired = []
    for name in selected:
        family, seed = name.rsplit('_s', 1)
        other = family + ('_s72' if seed == '71' else '_s71')
        if other in results and other not in selected:
            paired.append(other)
    output.mkdir(parents=True)
    (output / 'source_snapshot').mkdir()
    stats = rounds[0] / 'source_snapshot/bdt_performance.py'
    shutil.copy2(stats, output / 'source_snapshot/bdt_performance.py')
    shutil.copy2(Path(__file__), output / 'source_snapshot' / Path(__file__).name)
    protocol = {key: baseline[key] for key in ('cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits', 'acceptance_envelope')}
    protocol.update(schema='shift-dnn-multiround-comparison-v1', fit_count=len(results),
        rounds=[reference(directory / 'protocol.json') for directory, _ in protocols],
        sources={'bdt_performance.py': reference(stats), Path(__file__).name: reference(Path(__file__))},
        interpretation='All neural choices use validation only; no neural test scores exist before this lock. Previously inspected SM test is development evidence.',
        validation_selection='Same three predeclared roles as every round; paired seed is evaluated for stability, never reselected on test')
    write(output / 'protocol.json', protocol)
    write(output / 'validation_candidates.json', results)
    lock = {'schema': 'shift-dnn-global-validation-selection-v1', 'selection_uses_test': False,
        'rejection_winner': rejection, 'flatness_winner': flatness, 'acceptance_envelope_winner': acceptable,
        'acceptance_envelope_eligible': eligible, 'representative_models': selected, 'paired_seed_models': paired,
        'thresholds': {name: result['validation']['threshold'] for name, result in results.items()},
        'protocol': reference(output / 'protocol.json'), 'validation': reference(output / 'validation_candidates.json'),
        'fit_count': len(results), 'physics_ready': False}
    write(output / 'selection_locked.json', lock)
    # Explicitly pin the old reference lock, including its unchanged diagnostic cuts.
    previous = Path(baseline['bdt_scores']['path']).parent
    old_completion = json.loads((previous / 'completion.json').read_text())
    verify_reference(old_completion['selection'])
    old_lock = json.loads(Path(old_completion['selection']['path']).read_text())
    diagnostic = dict(lock)
    diagnostic['thresholds'] = dict(lock['thresholds'])
    for name in ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400'):
        diagnostic['thresholds'][name] = old_lock['thresholds'][name]
    diagnostic['previous_bdt_selection'] = old_completion['selection']
    write(output / 'diagnostic_selection.json', diagnostic)
    write(output / 'pre_evaluation_receipt.json', {'protocol': reference(output / 'protocol.json'),
        'selection': reference(output / 'selection_locked.json'), 'diagnostic_selection': reference(output / 'diagnostic_selection.json'),
        'validation': reference(output / 'validation_candidates.json'), 'previous_bdt_selection': old_completion['selection'],
        'neural_test_scores_generated': False})
    print(json.dumps({'locked_before_neural_test': True, 'fits': len(results), 'selected': selected,
                      'paired_seeds': paired, 'acceptance_eligible': eligible}), flush=True)


def evaluate(output):
    for filename in ('development_scores.npz', 'diagnostic_inputs.json', 'selected_development_metrics.json', 'reload_parity.json', 'completion.json'):
        require(not (output / filename).exists(), 'Refusing to overwrite evaluation artifact: ' + filename)
    before = json.loads((output / 'pre_evaluation_receipt.json').read_text())
    for key in ('protocol', 'selection', 'diagnostic_selection', 'validation', 'previous_bdt_selection'):
        verify_reference(before[key])
    protocol = json.loads((output / 'protocol.json').read_text())
    lock = json.loads((output / 'selection_locked.json').read_text())
    diagnostic = json.loads((output / 'diagnostic_selection.json').read_text())
    previous_lock = json.loads(Path(before['previous_bdt_selection']['path']).read_text())
    for name, cut in lock['thresholds'].items():
        require(diagnostic['thresholds'].get(name) == cut, 'Neural diagnostic cut changed')
    for name in ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400'):
        require(diagnostic['thresholds'].get(name) == previous_lock['thresholds'][name], 'Reference diagnostic cut changed')
    require(lock['selection_uses_test'] is False, 'Unsafe model selection')
    verify_reference(lock['protocol'])
    verify_reference(lock['validation'])
    for key in ('cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits'):
        verify_reference(protocol[key])
    for ref in protocol['sources'].values():
        verify_reference(ref)
    results = json.loads(Path(lock['validation']['path']).read_text())
    data, labels, previous = load(protocol['reco']['path']), load(protocol['labels']['path']), load(protocol['bdt_scores']['path'])
    names = ('ids', 'process', 'pair_index', 'source_index', 'heldout', 'expanded_split', 'expanded_heldout')
    scores = {name: previous[name] for name in names}
    for name in ('ids', 'process', 'pair_index', 'source_index'):
        require(np.array_equal(data[name], scores[name]), 'Evaluation identities changed')
    for name in ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400'):
        scores[name] = previous[name]
    selected = lock['representative_models'] + lock['paired_seed_models']
    tr = (scores['expanded_split'] == 0) & (scores['process'] != 'signal') & (labels['y'] >= 0)
    te = scores['expanded_heldout'] & (labels['y'] >= 0)
    metrics, parity = {}, {}
    for name in selected:
        result = results[name]
        verify_reference(result['model'])
        verify_reference(result['history'])
        verify_reference(result['protocol'])
        origin = json.loads(Path(result['protocol']['path']).read_text())
        for ref in origin['sources'].values():
            verify_reference(ref)
        with Path(result['model']['path']).open('rb') as stream:
            bundle = pickle.load(stream)
        require(tuple(bundle['feature_names']) == FEATURE_NAMES and np.array_equal(bundle['columns'], np.arange(50)), 'Inference feature contract changed')
        require(bundle['threshold'] == lock['thresholds'][name], 'Model cut differs from validation lock')
        predicted = bundle['model'].predict_proba(data['X'])[:, 1]
        with Path(result['model']['path']).open('rb') as stream:
            reloaded = pickle.load(stream)
        repeated = reloaded['model'].predict_proba(data['X'])[:, 1]
        require(np.array_equal(predicted, repeated), 'Trusted model reload score parity failed')
        require(predicted.shape == (len(data['ids']),) and np.isfinite(predicted).all(), 'Invalid neural score')
        scores[name] = predicted
        parity[name] = {'rows': len(predicted), 'bitwise_reloaded_score_parity': True, 'score_dtype': str(predicted.dtype)}
    np.savez_compressed(output / 'development_scores.npz', **scores)
    inputs = {key: protocol[key] for key in ('cache_manifest', 'reco', 'labels')}
    inputs.update(scores=reference(output / 'development_scores.npz'), selection=reference(output / 'diagnostic_selection.json'))
    write(output / 'diagnostic_inputs.json', inputs)
    for name in selected + ['reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400']:
        cut = diagnostic['thresholds'][name]
        metrics[name] = summary(scores[name][te], labels['y'][te], data['process'][te], data['nuisance'][te],
            data['nuisance'][tr], labels['y'][tr], cut, data['topology'][te])
    write(output / 'selected_development_metrics.json', metrics)
    write(output / 'reload_parity.json', parity)
    write(output / 'completion.json', {'complete': True, 'protocol': reference(output / 'protocol.json'),
        'selection': reference(output / 'selection_locked.json'), 'scores': reference(output / 'development_scores.npz'),
        'diagnostic_inputs': reference(output / 'diagnostic_inputs.json'), 'metrics': reference(output / 'selected_development_metrics.json'),
        'parity': reference(output / 'reload_parity.json'), 'selected_models': selected, 'fit_count': lock['fit_count'],
        'physics_ready': False, 'production_deployed': False})
    print(json.dumps({'evaluated': selected, 'test_known_pairs': int(te.sum()), 'score_parity': True}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--rounds', type=Path, nargs='+')
    parser.add_argument('--action', choices=('lock', 'evaluate'), required=True)
    args = parser.parse_args()
    if args.action == 'lock':
        require(args.rounds is not None, 'Explicit training rounds required')
        seal(args.output.resolve(), [path.resolve() for path in args.rounds])
    else:
        evaluate(args.output.resolve())


if __name__ == '__main__':
    main()
