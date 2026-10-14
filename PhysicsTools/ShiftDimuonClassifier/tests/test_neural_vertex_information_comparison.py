# Code style: AI-code. AI-generated paired-family selection regression tests.
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from neural_vertex_information_comparison import select_families, validate_rounds
from neural_vertex_information_refinement import recipes


def fixture():
    results = {}
    for family, reject, passing in [('onegood', [.9, .8], [True, False]), ('stable', [.5, .6], [True, True]), ('weaker', [.49, .8], [True, True])]:
        for seed, rejection, good in zip((71, 72), reject, passing):
            name = family + '_s' + str(seed)
            results[name] = dict(recipe=dict(name=name, family=family, seed=seed, architecture='global'),
                                 validation=dict(rejection=rejection, auc=.8, selector_objective=rejection, good=good))
    return results


class InformationComparison(unittest.TestCase):
    @patch('neural_vertex_information_comparison.envelope_eligible', side_effect=lambda metrics, envelope: metrics['good'])
    def test_both_seed_gate_and_minimum_rejection_rank(self, unused):
        report, lock = select_families(fixture(), {})
        self.assertFalse(report['onegood']['all_seeds_pass'])
        self.assertEqual(lock['robust_family'], 'stable')
        self.assertTrue({'stable_s71', 'stable_s72'} <= set(lock['representative_models']))
        self.assertIn('onegood_s72', lock['paired_seed_models'])
        self.assertEqual(lock['rejection_winner'], 'onegood_s71')

    @patch('neural_vertex_information_comparison.envelope_eligible', side_effect=lambda metrics, envelope: metrics['good'])
    def test_missing_seed_or_different_family_recipe_fails(self, unused):
        rows = fixture();del rows['stable_s72']
        with self.assertRaisesRegex(ValueError, 'exactly both seeds'):
            select_families(rows, {})
        rows = fixture();rows['stable_s72']['recipe']['architecture'] = 'heads'
        with self.assertRaisesRegex(ValueError, 'differs across seeds'):
            select_families(rows, {})

    @patch('neural_vertex_information_comparison.envelope_eligible', return_value=False)
    def test_no_eligible_family_has_no_test_fallback(self, unused):
        _, lock = select_families(fixture(), {})
        self.assertIsNone(lock['robust_family'])
        self.assertIsNone(lock['acceptance_envelope_winner'])
        self.assertEqual(lock['robust_eligible_families'], [])

    def test_combined_scope_requires_exact_original_and_refinement_rounds(self):
        with self.assertRaisesRegex(ValueError, 'Exactly original24'):
            validate_rounds([Path('/not/read')])

    def test_refinement_is_exact_eight_predeclared_recipes(self):
        rows = recipes()
        self.assertEqual(len(rows), 8)
        self.assertEqual(len({row['name'] for row in rows}), 8)
        self.assertEqual({(r['strength'], r['category_strength'], r['seed']) for r in rows},
                         {(j, c, s) for j in (1., 2.) for c in (1., 3.) for s in (71, 72)})
        self.assertTrue(all(r['architecture'] == 'global' and r['inputs'] == 'rich' for r in rows))


if __name__ == '__main__':
    unittest.main()
