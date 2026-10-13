# Code style: AI-code. AI-generated deterministic reco geometry contract tests.
import pickle
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from features import FEATURE_NAMES
from neural_vertex_geometry import (GEOMETRY_FEATURE_NAMES, GeometryVertexClassifier,
    _derived, derive_geometry_features, geometry_feature_names, geometry_inputs, train_geometry_vertex)


def fixture():
    x = np.ones((3, len(FEATURE_NAMES)))
    values = {'mu_highPt_px': 3, 'mu_highPt_py': 0, 'mu_highPt_pz': 4,
              'mu_lowPt_px': 0, 'mu_lowPt_py': 4, 'mu_lowPt_pz': 3,
              'mu_highPt_energy': 6, 'mu_lowPt_energy': 5,
              'mu_highPt_vx': 1, 'mu_highPt_vy': 2, 'mu_highPt_vz': 3,
              'mu_lowPt_vx': 4, 'mu_lowPt_vy': 6, 'mu_lowPt_vz': 9,
              'pair_vx': 2, 'pair_vy': 4, 'pair_vz': 8}
    for key, value in values.items():
        x[:, FEATURE_NAMES.index(key)] = value
    return x


class InputOnlyModel:
    """Pickleable deterministic stub: wrapper tests do not train a network."""
    def predict_proba(self, x):
        score = 1 / (1 + np.exp(-x[:, -12]))
        return np.column_stack((1-score, score))


class ReconstructedGeometry(unittest.TestCase):
    def test_relative_features_and_quality_mode_are_translation_invariant(self):
        x, shifted = fixture(), fixture()
        for axis, offset in zip(('vx', 'vy', 'vz'), (1000, -40, 370)):
            for prefix in ('pair_', 'mu_highPt_', 'mu_lowPt_'):
                shifted[:, FEATURE_NAMES.index(prefix + axis)] += offset
        np.testing.assert_allclose(derive_geometry_features(x), derive_geometry_features(shifted), atol=1e-12)
        np.testing.assert_allclose(geometry_inputs(x, 'quality_relative'), geometry_inputs(shifted, 'quality_relative'), atol=1e-12)

    def test_common_momentum_energy_scale_cancels(self):
        x, scaled = fixture(), fixture()
        for side in ('highPt', 'lowPt'):
            for field in ('px', 'py', 'pz', 'energy'):
                scaled[:, FEATURE_NAMES.index('mu_' + side + '_' + field)] *= 1000
        np.testing.assert_allclose(derive_geometry_features(x), derive_geometry_features(scaled), atol=1e-12)
        np.testing.assert_allclose(geometry_inputs(x, 'quality_relative'), geometry_inputs(scaled, 'quality_relative'), atol=1e-12)
        names = geometry_feature_names('quality_relative')
        self.assertEqual(len(names), 39)
        self.assertEqual(len(geometry_feature_names('all')), 62)
        self.assertFalse(any(name.startswith(('mu_highPt_', 'mu_lowPt_')) for name in names))

    def test_direction_reversal_changes_only_signed_geometry(self):
        x, reverse = fixture(), fixture()
        for field in ('px', 'py', 'pz'):
            reverse[:, FEATURE_NAMES.index('mu_highPt_' + field)] *= -1
        before, after = derive_geometry_features(x), derive_geometry_features(reverse)
        signed = [GEOMETRY_FEATURE_NAMES.index(name) for name in ('geometry_direction_cosine', 'geometry_reference_along_highPt')]
        for i in signed:
            np.testing.assert_allclose(after[:, i], -before[:, i])
        unsigned = [i for i in range(12) if i not in signed]
        np.testing.assert_allclose(after[:, unsigned], before[:, unsigned])

    def test_exact_line_distance_and_near_parallel_limit_are_finite(self):
        x = fixture()
        for side, vector in [('highPt', (1, 0, 0)), ('lowPt', (0, 1, 0))]:
            for axis, value in zip(('px', 'py', 'pz'), vector):
                x[:, FEATURE_NAMES.index('mu_' + side + '_' + axis)] = value
        for side, vertex in [('highPt', (0, 0, 0)), ('lowPt', (0, 0, 3))]:
            for axis, value in zip(('vx', 'vy', 'vz'), vertex):
                x[:, FEATURE_NAMES.index('mu_' + side + '_' + axis)] = value
        derived = derive_geometry_features(x)
        np.testing.assert_allclose(derived[:, 8], 3)
        x[:, FEATURE_NAMES.index('mu_lowPt_px')] = 1
        x[:, FEATURE_NAMES.index('mu_lowPt_py')] = [0, 1e-12, -1e-12]
        derived, support = _derived(x)
        self.assertTrue(np.isfinite(derived).all())
        np.testing.assert_allclose(derived[:, 8], 3)
        self.assertEqual(support['near_parallel_rows'], 3)

    def test_invalid_direction_and_error_stay_missing_with_diagnostics(self):
        x = fixture()
        for field in ('px', 'py', 'pz'):
            x[0, FEATURE_NAMES.index('mu_highPt_' + field)] = 0
        for field in ('vxErr', 'vyErr', 'vzErr'):
            x[1, FEATURE_NAMES.index('pair_' + field)] = 0
        derived, support = _derived(x)
        self.assertTrue(np.isnan(derived[0, 0]))
        self.assertTrue(np.isnan(derived[1, -1]))
        self.assertEqual(support['invalid_direction_rows'], 1)
        self.assertEqual(support['invalid_error_ratio_rows'], 1)

    def test_wrapper_raw_contract_pickle_and_forbidden_energy_penalty(self):
        model = GeometryVertexClassifier(InputOnlyModel(), 'quality_relative')
        restored = pickle.loads(pickle.dumps(model))
        np.testing.assert_array_equal(model.predict_proba(fixture()), restored.predict_proba(fixture()))
        self.assertEqual(model.predict_proba(np.empty((0, 50))).shape, (0, 2))
        for bad in (np.ones(50), np.ones((2, 49)), np.ones((2, 62))):
            with self.assertRaisesRegex(ValueError, 'exactly50'):
                model.predict_proba(bad)
        with self.assertRaisesRegex(ValueError, 'energy_strength=0'):
            train_geometry_vertex(*([None] * 8), energy_strength=.3)
        self.assertEqual(set(model.__dict__), {'model', 'mode', 'n_features_in_', 'classes_'})


if __name__ == '__main__':
    unittest.main()
