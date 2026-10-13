#!/usr/bin/env python3
# Code style: AI-code. AI-generated validation-guided DNN refinement protocol.
"""Refine penalties after round-one validation, before viewing any DNN test scores."""
import argparse
import json
from pathlib import Path
import shutil

from neural_vertex_study import prepare, fit_one, reference, require, write


def recipes():
    plans = []
    settings = [(0., 0., 0., [64, 32]), (1., 1., 0., [64, 32]),
                (1., 3., 0., [64, 32]), (1., 10., 0., [64, 32]),
                (3., 1., 0., [64, 32]), (3., 3., 0., [64, 32]),
                (3., 10., 0., [64, 32]), (.3, 3., 0., [64, 32]),
                (1., 3., 1., [64, 32]), (1., 3., 0., [96, 48])]
    for strength, category, energy, hidden in settings:
        for seed in (71, 72):
            plans.append({'name': f'dnn_r2_l{strength:g}_c{category:g}_e{energy:g}_h{hidden[0]}_s{seed}',
                'strength': strength, 'category_strength': category, 'energy_strength': energy,
                'hidden': hidden, 'feature_transform': 'signedlog', 'penalty_scope': 'class_process', 'seed': seed})
    return plans


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--action', choices=('prepare', 'fit'), default='prepare')
    args = parser.parse_args()
    workspace, output, previous = args.workspace.resolve(), args.output.resolve(), args.previous.resolve()
    require(output.is_relative_to(workspace / 'tea_shift_cmssw/plots'), 'Final study belongs in repository')
    if args.action == 'prepare':
        require(not (previous / 'development_scores.npz').exists(), 'Refinement must precede DNN test-score inspection')
        previous_protocol = json.loads((previous / 'protocol.json').read_text())
        previous_results = [reference(previous / (spec['name'] + '_validation.json')) for spec in previous_protocol['recipes']]
        require(len(previous_results) == 22, 'Expected complete round-one validation grid')
        protocol, *_ = prepare(workspace, output)
        # Preserve the base preparation receipt before sealing this phase's actual fit plan.
        (output / 'protocol.json').rename(output / 'initial_preparation_protocol.json')
        protocol['recipes'] = recipes()
        protocol['fit_budget'] = len(protocol['recipes'])
        protocol['training'].update(epochs=240, patience=35, learning_rate=.0005, dropout=.1, weight_decay=.0005)
        protocol['adaptation'] = {'uses': 'Only round-one training and validation diagnostics',
            'reason': 'First grid improved separation and reduced joint dependence but failed category and worst-bin acceptance gates. Test stronger genuine-category and joint penalties with slower regularized optimization.',
            'previous_protocol': reference(previous / 'protocol.json'), 'validation_receipts': previous_results,
            'no_previous_DNN_test_scores': True, 'base_preparation': reference(output / 'initial_preparation_protocol.json')}
        protocol['sources'][Path(__file__).name] = reference(Path(__file__))
        shutil.copy2(Path(__file__), output / 'source_snapshot' / Path(__file__).name)
        write(output / 'protocol.json', protocol)
        print(json.dumps({'prepared': True, 'fits': len(protocol['recipes'])}), flush=True)
    else:
        protocol = json.loads((output / 'protocol.json').read_text())
        for spec in protocol['recipes']:
            if not (output / (spec['name'] + '_validation.json')).exists():
                fit_one(output, spec)
        print('All refinement validation fits complete; no DNN test scores evaluated.', flush=True)


if __name__ == '__main__':
    main()
