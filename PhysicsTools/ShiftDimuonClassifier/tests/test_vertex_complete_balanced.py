# Code style: AI-code. AI-generated fixed loss-setting and control tests.
import importlib
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import vertex_complete_study as base
import vertex_complete_balanced as balanced


class BalancedMenuTests(unittest.TestCase):
    def test_controls_and_nested_grid_are_preserved(self):
        menu = balanced.recipes()
        self.assertEqual(len(menu), 81)
        self.assertEqual(len({row['name'] for row in menu}), 81)
        self.assertEqual([row for row in menu if row['algorithm'] == 'bdt'],
                         [row for row in base.recipes() if row['algorithm'] == 'bdt'])
        for family in {row['family'] for row in menu}:
            self.assertEqual({(row['fraction'], row['seed']) for row in menu if row['family'] == family},
                             {(f, s) for f in base.FRACTIONS for s in base.SEEDS})
        for row in menu:
            if row['algorithm'] == 'dnn':
                self.assertEqual((row['strength'], row['category_strength'], row['conditional_strength']),
                                 (3., 3., 0.))
                self.assertEqual(row['acceptance_strength'], 30. if 'accept30' in row['loss'] else 0.)

    def test_import_does_not_modify_archived_driver_menu(self):
        before = base.recipes(), base.SOURCES, base.PROTOCOL_KIND
        importlib.reload(balanced)
        self.assertEqual((base.recipes(), base.SOURCES, base.PROTOCOL_KIND), before)


if __name__ == '__main__':
    unittest.main()
