#!/usr/bin/env python3
# Code style: AI-code. AI-generated complete-input loss-setting refinement.
"""Compare moderate decorrelation and direct soft acceptance with full RECO.

This is a separate frozen validation study, not a change to archived fits.
The same BDT controls, partitions, nested training subsets, bounds and ranking
are retained. No development or signal result defines this fixed menu.
"""
import vertex_complete_study as base


def recipes():
    families = [dict(family='bdt_' + config, algorithm='bdt', config=config)
                for config in ('shallow_d3', 'medium_d5', 'regularized_d5')]
    for architecture, hidden, label in [('shared', [64, 32], 'shared64'),
                                        ('shared', [128, 64], 'shared128'),
                                        ('heads', [64, 32], 'heads64')]:
        for loss, joint, category, conditional, acceptance in [
                ('moderate', .3, .3, .1, 0.),
                ('acceptance10', .3, 0., .3, 10.)]:
            families.append(dict(family='dnn_' + label + '_' + loss, algorithm='dnn',
                architecture=architecture, hidden=hidden, loss=loss,
                strength=joint, category_strength=category,
                conditional_strength=conditional, acceptance_strength=acceptance))
    return [{**family, 'fraction': size, 'seed': seed,
             'name': f"{family['family']}_f{int(size*100)}_s{seed}"}
            for family in families for size in base.FRACTIONS for seed in base.SEEDS]


def configure():
    base.recipes = recipes
    base.PROTOCOL_KIND = 'moderate_decorrelation_and_direct_acceptance'
    base.SOURCES = (*base.SOURCES, 'vertex_complete_refinement.py')


if __name__ == '__main__':
    configure()
    base.main()
