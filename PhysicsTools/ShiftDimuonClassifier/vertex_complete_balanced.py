#!/usr/bin/env python3
# Code style: AI-code. AI-generated complete-input balanced-loss study.
"""Keep explicit category inputs, test the successful pooled regularization.

Earlier explicit-category recipes added conditional penalties absent from the
successful rich pooled model. This fixed menu removes that additional loss,
with and without direct working-point regularization. All RECO inputs remain.
Choices and cuts use validation only; signals are subsequent response checks.
"""
import vertex_complete_study as base


def recipes():
    families = [dict(family='bdt_' + config, algorithm='bdt', config=config)
                for config in ('shallow_d3', 'medium_d5', 'regularized_d5')]
    for architecture, hidden, label in [('shared', [64, 32], 'shared64'),
                                        ('shared', [128, 64], 'shared128'),
                                        ('heads', [64, 32], 'heads64')]:
        for loss, acceptance in [('joint3', 0.), ('joint3_accept30', 30.)]:
            families.append(dict(family='dnn_' + label + '_' + loss, algorithm='dnn',
                architecture=architecture, hidden=hidden, loss=loss,
                strength=3., category_strength=3., conditional_strength=0.,
                acceptance_strength=acceptance))
    return [{**family, 'fraction': size, 'seed': seed,
             'name': f"{family['family']}_f{int(size*100)}_s{seed}"}
            for family in families for size in base.FRACTIONS for seed in base.SEEDS]


if __name__ == '__main__':
    base.recipes = recipes
    base.PROTOCOL_KIND = 'joint_category_penalties_without_extra_conditional_loss'
    base.SOURCES = (*base.SOURCES, 'vertex_complete_balanced.py')
    base.main()
