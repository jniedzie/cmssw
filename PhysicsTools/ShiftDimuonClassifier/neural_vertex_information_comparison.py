#!/usr/bin/env python3
# Code style: AI-code. AI-generated cross-round paired-seed validation assembly.
"""Select families before any new study test scores; preserve every fit origin."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np

from neural_vertex_comparison import gather
from neural_vertex_information_study import evaluate as evaluate_base, verify, recipes as original_recipes
from neural_vertex_information_refinement import recipes as refinement_recipes
from neural_vertex_study import envelope_eligible, load, reference, require, write


def select_families(results, envelope):
    families = {}
    for name, result in results.items():
        families.setdefault(result['recipe']['family'], []).append(name)
    report = {}
    for family, names in families.items():
        names = sorted(names, key=lambda name: results[name]['recipe']['seed'])
        require(len(names) == 2 and {results[name]['recipe']['seed'] for name in names} == {71, 72}, 'Family requires exactly both seeds')
        base = {key: value for key, value in results[names[0]]['recipe'].items() if key not in ('name', 'seed')}
        require(all({key: value for key, value in results[name]['recipe'].items() if key not in ('name', 'seed')} == base for name in names), 'Family recipe differs across seeds')
        passing = [envelope_eligible(results[name]['validation'], envelope) for name in names]
        report[family] = dict(models=names, passing=passing, all_seeds_pass=all(passing),
            minimum_rejection=min(results[name]['validation']['rejection'] for name in names),
            minimum_auc=min(results[name]['validation']['auc'] for name in names))
    eligible = sorted(family for family, row in report.items() if row['all_seeds_pass'])
    robust = max(eligible, key=lambda f: (report[f]['minimum_rejection'], report[f]['minimum_auc'], f)) if eligible else None
    rejection = max(results, key=lambda n: (results[n]['validation']['rejection'], results[n]['validation']['auc'], n))
    flatness = max(results, key=lambda n: (results[n]['validation']['selector_objective'], n))
    selected = list(dict.fromkeys((report[robust]['models'] if robust else []) + [rejection, flatness]))
    paired = []
    for name in selected:
        for other in report[results[name]['recipe']['family']]['models']:
            if other not in selected and other not in paired:
                paired.append(other)
    return report, dict(rejection_winner=rejection, flatness_winner=flatness,
        acceptance_envelope_winner=report[robust]['models'][0] if robust else None,
        robust_family=robust, robust_eligible_families=eligible, representative_models=selected, paired_seed_models=paired)


def validate_rounds(rounds):
    require(len(rounds) == 2, 'Exactly original24 and refinement8 rounds required')
    results, protocols = gather(rounds)
    baseline = protocols[0][1]
    require(baseline['recipes'] == original_recipes() and baseline['fit_count'] == 24, 'Original24 recipe grid changed')
    refinement = protocols[1][1]
    require(refinement['recipes'] == refinement_recipes() and refinement['fit_count'] == 8, 'Refinement8 recipe grid changed')
    require(len(results) == 32, 'Expected exactly32 original plus refinement fits')
    require(refinement['adaptation']['prior_protocol'] == reference(rounds[0] / 'protocol.json'), 'Refinement adaptation points to another original round')
    splits = load(baseline['splits']['path'])
    for directory, protocol in protocols:
        for key in ('cache_manifest', 'reco', 'labels', 'events', 'bdt_scores', 'previous_bdt_selection'):
            verify(protocol[key])
            require(protocol[key] == baseline[key], 'Different fitting data: ' + key)
        for key in ('acceptance_envelope', 'input_feature_names', 'training'):
            require(protocol[key] == baseline[key], 'Different controlled protocol: ' + key)
        other = load(protocol['splits']['path'])
        require(set(splits) == set(other) and all(np.array_equal(splits[key], other[key]) for key in splits), 'Physical-event split changed')
        if 'adaptation' in protocol:
            verify(protocol['base_prepare_protocol'])
            verify(protocol['adaptation']['prior_protocol'])
            for ref in protocol['adaptation']['validation_evidence']:
                verify(ref)
    return results, protocols


def seal(output, rounds):
    require(not output.exists(), 'Refusing to overwrite combined selection')
    results, protocols = validate_rounds(rounds)
    baseline = protocols[0][1]
    families, roles = select_families(results, baseline['acceptance_envelope'])
    output.mkdir(parents=True)
    (output / 'source_snapshot').mkdir()
    sources = {}
    for directory, protocol in protocols:
        for name, ref in protocol['sources'].items():
            require(name not in sources or sources[name] == ref, 'Different source versions across controlled rounds')
            if name not in sources:
                shutil.copy2(directory / 'source_snapshot' / name, output / 'source_snapshot' / name)
                sources[name] = ref
    for filename in ('neural_vertex_information_comparison.py', 'neural_vertex_comparison.py'):
        path = Path(__file__).parent / filename
        sources[filename] = reference(path)
        shutil.copy2(path, output / 'source_snapshot' / filename)
    protocol = dict(baseline)
    protocol.update(schema='shift-rich-neural-multiround-comparison-v1', fit_count=len(results),
        recipes=[result['recipe'] for result in results.values()], sources=sources,
        rounds=[reference(directory / 'protocol.json') for directory, _ in protocols],
        validation_selection='Both predeclared seeds pass unchanged envelope; maximize minimum rejection, then minimum AUC. Preserve rejection and flatness descriptive winners. No test selection.')
    write(output / 'protocol.json', protocol)
    write(output / 'validation_candidates.json', results)
    write(output / 'family_validation.json', families)
    lock = dict(schema='shift-rich-neural-global-family-selection-v1', selection_uses_test=False, **roles,
        thresholds={name: row['validation']['threshold'] for name, row in results.items()},
        protocol=reference(output / 'protocol.json'), validation=reference(output / 'validation_candidates.json'),
        family_validation=reference(output / 'family_validation.json'), fit_count=len(results), physics_ready=False)
    write(output / 'selection_locked.json', lock)
    diagnostic = dict(lock, thresholds=dict(lock['thresholds']))
    previous = json.loads(Path(protocol['previous_bdt_selection']['path']).read_text())
    for name in ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400'):
        diagnostic['thresholds'][name] = previous['thresholds'][name]
    write(output / 'diagnostic_selection.json', diagnostic)
    write(output / 'pre_evaluation_receipt.json', dict(protocol=reference(output / 'protocol.json'),
        selection=reference(output / 'selection_locked.json'), diagnostic_selection=reference(output / 'diagnostic_selection.json'),
        validation=reference(output / 'validation_candidates.json'), family_validation=reference(output / 'family_validation.json'),
        new_test_scores_generated=False))
    print(json.dumps(dict(locked=True, fits=len(results), **roles)), flush=True)


def evaluate(output):
    before = json.loads((output / 'pre_evaluation_receipt.json').read_text())
    for ref in before.values():
        if isinstance(ref, dict):
            verify(ref)
    protocol = json.loads((output / 'protocol.json').read_text())
    for ref in protocol['rounds']:
        verify(ref)
    rounds = [Path(ref['path']).parent for ref in protocol['rounds']]
    results, _ = validate_rounds(rounds)
    require(results == json.loads((output / 'validation_candidates.json').read_text()), 'Original validation receipts changed after global lock')
    for name, ref in protocol['sources'].items():
        verify(ref)
        require(reference(output / 'source_snapshot' / name)['sha256'] == ref['sha256'], 'Combined source snapshot changed')
    evaluate_base(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--rounds', type=Path, nargs='+')
    parser.add_argument('--action', choices=('lock', 'evaluate'), required=True)
    args = parser.parse_args()
    if args.action == 'lock':
        require(args.rounds is not None and len(args.rounds) >= 2, 'Explicit original and refinement rounds required')
        seal(args.output.resolve(), [path.resolve() for path in args.rounds])
    else:
        evaluate(args.output.resolve())


if __name__ == '__main__':
    main()
