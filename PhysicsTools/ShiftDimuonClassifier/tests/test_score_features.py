"""Inference keeps all candidates and never needs the mass target."""
import copy
import unittest
import numpy as np
from features import PREDICTOR_BRANCHES, extract_reco
from test_contract import reconstructed_event


class ScoreFeatureContract(unittest.TestCase):
    def test_missing_mass_does_not_change_predictors(self):
        arrays = reconstructed_event()
        expected = extract_reco(arrays)
        del arrays['ShiftDimuonVertex_mass']
        result = extract_reco(arrays, include_nuisance=False)
        self.assertNotIn('ShiftDimuonVertex_mass', PREDICTOR_BRANCHES)
        np.testing.assert_array_equal(result['X'], expected['X'])
        np.testing.assert_array_equal(result['pair_index'], expected['pair_index'])
        self.assertIsNone(result['nuisance'])

    def test_invalid_target_cannot_reject_a_scored_row(self):
        arrays = reconstructed_event()
        expected = extract_reco(arrays)['X']
        for mass in (-100., np.nan, np.inf):
            value = copy.deepcopy(arrays)
            value['ShiftDimuonVertex_mass'][0][0] = mass
            with self.assertRaises(ValueError):
                extract_reco(value)
            np.testing.assert_array_equal(extract_reco(value, include_nuisance=False)['X'], expected)

    def test_nonfinite_geometry_or_angle_is_preserved_for_frozen_imputation(self):
        arrays = reconstructed_event()
        arrays['ShiftDimuonVertex_vz'][0][0] = np.nan
        arrays['ShiftMuon_phi'][0][0] = np.inf
        result = extract_reco(arrays, include_nuisance=False)
        self.assertEqual(len(result['X']), 1)
        self.assertTrue(np.isnan(result['X']).any())


if __name__ == '__main__':
    unittest.main()
