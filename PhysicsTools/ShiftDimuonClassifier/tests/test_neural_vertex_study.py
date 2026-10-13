# Code style: AI-code. AI-generated checks of frozen DNN evaluation cuts.
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from neural_vertex_study import candidates, envelope_eligible, summary


class StudyMetricsTests(unittest.TestCase):
    def test_supplied_validation_cut_controls_every_acceptance_profile(self):
        rng = np.random.default_rng(42)
        y = np.r_[np.ones(180, dtype=int), np.zeros(80, dtype=int)]
        process = np.array(['dy'] * 90 + ['jpsi'] * 90 + ['qcd'] * 80)
        category = np.array(['Both-Both'] * 90 + ['Near-Both'] * 90 + ['Other'] * 80)
        nuisance = rng.uniform([1, 100, 1], [10, 200, 100], (260, 3))
        score = np.linspace(.2, .8, 260)
        result = summary(score, y, process, nuisance, nuisance, y, threshold=.75, category=category)
        self.assertEqual(result['threshold'], .75)
        self.assertEqual(result['common_vertex']['passed'], 0)
        self.assertEqual(result['balanced_common_efficiency'], 0)
        for profiles in [result['flatness'], *result['per_process_flatness'].values()]:
            for dimension in ('mass', 'vertex_z', 'vertex_radius'):
                self.assertTrue(all(row['passed'] == 0 for row in profiles[dimension]['bins']))
        self.assertTrue(all(row['common_vertex']['passed'] == 0 for row in result['categories'].values()))

    def test_category_support_is_required_by_acceptance_gate(self):
        profiles = {name: {'bins': [{'n': 100, 'efficiency': .9, 'supported': True}]}
                    for name in ('mass', 'vertex_z', 'vertex_radius')}
        entry = {'maximum_marginal_efficiency_span': .05, 'process_efficiency_span': .01,
                 'flatness': profiles, 'per_process_flatness': {},
                 'dependence': {'common_vertex': {'joint_dcor': .1}, 'accidental': {'joint_dcor': .1}},
                 'minimum_supported_category_efficiency': None, 'category_efficiency_span': None}
        bounds = {'minimum_supported_genuine_retention': .8, 'maximum_supported_marginal_span': .15,
                  'maximum_process_efficiency_span': .05, 'maximum_class_joint_dcor': .3,
                  'maximum_category_efficiency_span': .15}
        self.assertFalse(envelope_eligible(entry, bounds))
        entry.update(minimum_supported_category_efficiency=.85, category_efficiency_span=.1)
        self.assertTrue(envelope_eligible(entry, bounds))
        entry['minimum_supported_category_efficiency'] = .65
        self.assertFalse(envelope_eligible(entry, bounds))

    def test_recipe_names_are_unique_and_use_paired_seeds(self):
        recipes = candidates()
        self.assertEqual(len(recipes), 22)
        self.assertEqual(len({row['name'] for row in recipes}), 22)
        families = {}
        for row in recipes:
            family = row['name'].rsplit('_s', 1)[0]
            families.setdefault(family, set()).add(row['seed'])
        self.assertTrue(all(seeds == {71, 72} for seeds in families.values()))


if __name__ == '__main__':
    unittest.main()
