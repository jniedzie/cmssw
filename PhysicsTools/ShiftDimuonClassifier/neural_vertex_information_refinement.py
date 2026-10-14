#!/usr/bin/env python3
# Code style: AI-code. AI-generated validation-only rich-input refinement driver.
"""Eight preregistered global-rich fits; the original study remains untouched."""
import argparse
import json
from pathlib import Path
import shutil

from neural_vertex_information_study import prepare as prepare_base, fit, verify
from neural_vertex_study import reference, require, write


def recipes():
    return [dict(name=f'info_r2_global_rich_l{joint}_c{category}_s{seed}',
                 family=f'r2_global_rich_l{joint}_c{category}', architecture='global', inputs='rich',
                 strength=float(joint), category_strength=float(category), conditional_strength=0.,
                 acceptance_strength=0., seed=seed)
            for joint in (1, 2) for category in (1, 3) for seed in (71, 72)]


def prepare(workspace, cache, output, prior):
    require(not (prior / 'development_scores.npz').exists(), 'Refinement must precede new study test inference')
    prior_protocol = json.loads((prior / 'protocol.json').read_text())
    evidence = []
    for spec in prior_protocol['recipes']:
        path = prior / (spec['name'] + '_validation.json')
        if not path.exists():
            continue
        entry = json.loads(path.read_text())
        require(entry['recipe'] == spec and entry['protocol'] == reference(prior / 'protocol.json'), 'Changed prior validation origin')
        for key in ('history', 'model'):
            verify(entry[key])
        evidence.append(reference(path))
    require(bool(evidence), 'No completed validation evidence for refinement')
    prepare_base(workspace, cache, output)
    initial = output / 'protocol.json'
    initial.rename(output / 'base_prepare_protocol.json')
    protocol = json.loads((output / 'base_prepare_protocol.json').read_text())
    for key in ('cache_manifest', 'reco', 'labels', 'events', 'bdt_scores'):
        require(protocol[key] == prior_protocol[key], 'Refinement changed data: ' + key)
    path = Path(__file__).resolve()
    shutil.copy2(path, output / 'source_snapshot' / path.name)
    protocol['sources'][path.name] = reference(path)
    protocol.update(schema='shift-rich-vertex-neural-refinement-v1', recipes=recipes(), fit_count=8,
        base_prepare_protocol=reference(output / 'base_prepare_protocol.json'),
        adaptation=dict(prior_protocol=reference(prior / 'protocol.json'), validation_evidence=evidence,
            new_study_test_inference=False,
            reason='Round1 rich global class_process models pass both seeds; compare smaller joint/category coefficients on validation only.'),
        selection='Combined original24 and refinement8 families; both seeds must pass unchanged gates; rank minimum rejection, then minimum AUC. No test-based choices.')
    write(initial, protocol)
    print(json.dumps(dict(refinement_prepared=True, fits=8, validation_receipts=len(evidence))), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--prior', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--action', choices=('prepare', 'fit'), required=True)
    parser.add_argument('--recipe')
    parser.add_argument('--shard', type=int, default=0)
    parser.add_argument('--shards', type=int, default=1)
    args = parser.parse_args()
    workspace, output = args.workspace.resolve(), args.output.resolve()
    require(output.is_relative_to(workspace / 'tea_shift_cmssw/plots'), 'Results must be inside project plots')
    require(0 <= args.shard < args.shards, 'Invalid shard')
    if args.action == 'prepare':
        require(args.cache is not None and args.prior is not None, 'Explicit cache and prior training directory required')
        prepare(workspace, args.cache.resolve(), output, args.prior.resolve())
    else:
        protocol = json.loads((output / 'protocol.json').read_text())
        require(protocol['recipes'] == recipes(), 'Refinement recipe grid changed')
        verify(protocol['base_prepare_protocol'])
        verify(protocol['adaptation']['prior_protocol'])
        for ref in protocol['adaptation']['validation_evidence']:
            verify(ref)
        chosen = [spec for i, spec in enumerate(protocol['recipes']) if i % args.shards == args.shard
                  and (args.recipe is None or spec['name'] == args.recipe)]
        require(bool(chosen), 'No matching declared recipes')
        for spec in chosen:
            if not (output / (spec['name'] + '_validation.json')).exists():
                fit(output, spec)


if __name__ == '__main__':
    main()
