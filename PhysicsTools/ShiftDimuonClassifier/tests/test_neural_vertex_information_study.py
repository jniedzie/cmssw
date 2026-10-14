# Code style: AI-code. AI-generated information-study contract tests.
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from features import FEATURE_NAMES
from neural_vertex_information_study import evaluate, input_columns, predict, recipes


class InformationStudyTest(unittest.TestCase):
    def test_recipe_families_have_exactly_paired_seeds(self):
        families = {}
        for recipe in recipes():
            families.setdefault(recipe['family'], []).append(recipe)
        self.assertEqual(len(families), 12)
        for members in families.values():
            self.assertEqual({entry['seed'] for entry in members}, {71, 72})
            first, second = ({k: v for k, v in entry.items() if k not in ('name', 'seed')} for entry in members)
            self.assertEqual(first, second)

    def test_input_ablation_preserves_pair_features_and_removes_only_event_extension(self):
        names = list(FEATURE_NAMES) + ['mu_highPt_ptErr', 'event_nShiftMuon']
        np.testing.assert_array_equal(input_columns(names, 'raw50'), np.arange(50))
        np.testing.assert_array_equal(input_columns(names, 'muon'), np.arange(51))
        np.testing.assert_array_equal(input_columns(names, 'rich'), np.arange(52))

    def test_invalid_input_prefix_rejected(self):
        names = list(FEATURE_NAMES)
        names[0] = 'GenPart_pt'
        with self.assertRaisesRegex(ValueError, 'prefix'):
            input_columns(names, 'rich')

    def test_prediction_schema_rejected_before_model_is_called(self):
        reco = {'feature_names': np.array(['wrong']), 'X': np.ones((1, 1))}
        bundle = {'input_feature_names': ['correct']}
        with self.assertRaisesRegex(ValueError, 'feature contract'):
            predict(bundle, reco)

    def test_partial_evaluation_is_preserved_before_inputs_are_read(self):
        for filename in ('development_scores.npz', 'diagnostic_inputs.json', 'selected_development_metrics.json', 'reload_parity.json', 'completion.json'):
            with TemporaryDirectory() as directory:
                output = Path(directory)
                target = output / filename
                target.write_bytes(b'previous evidence')
                with self.assertRaisesRegex(ValueError, 'overwrite'):
                    evaluate(output)
                self.assertEqual(target.read_bytes(), b'previous evidence')


if __name__ == '__main__':
    unittest.main()
