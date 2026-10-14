# Code style: AI-code. AI-generated fixed-menu and base-control tests.
import importlib
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import vertex_complete_study as base
import vertex_complete_refinement as refinement


class RefinementTests(unittest.TestCase):
    def test_menu_keeps_controls_subsets_and_all_three_architectures(self):
        original = base.recipes()
        refined = refinement.recipes()
        self.assertEqual(len(refined), 81)
        self.assertEqual(len({row['name'] for row in refined}), 81)
        self.assertEqual([row for row in refined if row['algorithm'] == 'bdt'],
                         [row for row in original if row['algorithm'] == 'bdt'])
        self.assertEqual({row['architecture'] for row in refined if row['algorithm'] == 'dnn'},
                         {'shared', 'heads'})
        self.assertEqual({tuple(row['hidden']) for row in refined if row['algorithm'] == 'dnn'},
                         {(64, 32), (128, 64)})
        for family in {row['family'] for row in refined}:
            self.assertEqual({(row['fraction'], row['seed']) for row in refined if row['family'] == family},
                             {(fraction, seed) for fraction in base.FRACTIONS for seed in base.SEEDS})

    def test_import_does_not_change_original_frozen_menu(self):
        before, sources = base.recipes(), base.SOURCES
        importlib.reload(refinement)
        self.assertEqual(base.recipes(), before)
        self.assertEqual(base.SOURCES, sources)
        for row in refinement.recipes():
            if row['algorithm'] == 'dnn':
                self.assertGreater(row['strength'], 0.)
                self.assertEqual(row['acceptance_strength'], 10. if row['loss'] == 'acceptance10' else 0.)


if __name__ == '__main__':
    unittest.main()
