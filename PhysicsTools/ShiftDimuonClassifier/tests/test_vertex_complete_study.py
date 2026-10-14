# Code style: AI-code. AI-generated independent split and receipt guards.
"""Synthetic protocol tests without physics fitting, model loading or ROOT."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import vertex_complete_study as study


class CompleteStudyTests(unittest.TestCase):
    def test_predeclared_complete_grid_and_zero_loss_control(self):
        recipes = study.recipes()
        self.assertEqual(len(recipes), 81)
        self.assertEqual(len({spec['name'] for spec in recipes}), 81)
        self.assertEqual(len({spec['family'] for spec in recipes}), 9)
        for family in {spec['family'] for spec in recipes}:
            members = [spec for spec in recipes if spec['family'] == family]
            self.assertEqual({(spec['fraction'], spec['seed']) for spec in members},
                             {(fraction, seed) for fraction in study.FRACTIONS for seed in study.SEEDS})
            for spec in members:
                if spec.get('loss') == 'classification':
                    self.assertEqual((spec['strength'], spec['category_strength'], spec['conditional_strength']),
                                     (0., 0., 0.))

    def test_nested_subsets_never_cross_fixed_training_mask(self):
        keys = np.asarray([f'physical:{index // 2}' for index in range(400)])
        train = np.arange(400) < 200
        for seed in study.SEEDS:
            small, middle, full = [study.nested_subset(keys, train, seed, size) for size in study.FRACTIONS]
            self.assertFalse(np.any(small & ~middle))
            self.assertFalse(np.any(middle & ~full))
            np.testing.assert_array_equal(full, train)
            for mask in (small, middle, full):
                self.assertFalse(np.any(mask & ~train))
                np.testing.assert_array_equal(mask[::2], mask[1::2])
        self.assertFalse(np.array_equal(study.nested_subset(keys, train, 71, .5),
                                       study.nested_subset(keys, train, 72, .5)))

    def test_subset_is_permutation_invariant_and_namespaced(self):
        keys = np.asarray([f'qcd:{index}' for index in range(100)])
        train = np.arange(100) % 3 == 0
        order = np.random.default_rng(14).permutation(len(keys))
        original = study.nested_subset(keys, train, 71, .5)
        shuffled = study.nested_subset(keys[order], train[order], 71, .5)
        np.testing.assert_array_equal(shuffled, original[order])
        self.assertFalse(np.array_equal(study.hash_values(keys, 71, 'nested-training'),
                                       study.hash_values(keys, 71, 'historical-partition')))

    def test_physical_keys_keep_exact_ids_and_historical_source_context(self):
        ids = np.asarray([[1, 2, 2**53], [1, 2, 2**53 + 1], [1, 2, 2**53]], dtype=np.uint64)
        reco = dict(ids=ids, process=np.asarray(['dy', 'dy', 'qcd']))
        keys = study.parent_keys(reco, 'current_CMSIR5')
        self.assertEqual(len(set(keys)), 3)
        historical = dict(reco, physical_parent_key=np.asarray(['sourceA:1', 'sourceA:2', 'sourceB:1']))
        np.testing.assert_array_equal(study.parent_keys(historical, 'historical_proxy'), historical['physical_parent_key'])
        with self.assertRaises(ValueError):
            study.parent_keys(reco, 'unknown-domain')

    def test_complete_fits_rejects_cross_domain_protocol_before_pickle(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'current'
            (output / 'fits').mkdir(parents=True)
            protocol = dict(recipes=study.recipes(), domain='current_CMSIR5')
            study.write(output / 'protocol.json', protocol)
            foreign = root / 'historical_protocol.json'
            study.write(foreign, dict(protocol, domain='historical_proxy'))
            dummy = root / 'not_a_model.pkl'
            dummy.write_text('No pickle should be loaded during validation readiness checks.')
            history = root / 'history.json'
            study.write(history, {})
            for spec in protocol['recipes']:
                study.write(output / 'fits' / (spec['name'] + '.json'), dict(recipe=spec,
                    protocol=study.reference(foreign), model=study.reference(dummy), history=study.reference(history)))
            with patch.object(study.pickle, 'load', side_effect=AssertionError('Unexpected pickle read')):
                with self.assertRaises(ValueError):
                    study.complete_fits(output)

    def test_no_fit_after_lock_even_before_input_reads(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / 'selection_locked.json').write_text('{}')
            with patch.object(study, 'fitting_data', side_effect=AssertionError('Data read after lock')):
                with self.assertRaisesRegex(ValueError, 'closed'):
                    study.fit(output, study.recipes()[0])


if __name__ == '__main__':
    unittest.main()
