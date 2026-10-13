# Code style: AI-code. AI-generated reconstructed relative-geometry benchmark.
"""Twelve deterministic reco features and a raw-X50 neural inference wrapper.

Straight-line distances are approximations to track geometry, not helix DCA.
The DCA/error scale ratio is not a calibrated significance. No per-leg pt error
is reconstructed from independently sorted min/max error columns. Absolute
coordinate translations cancel in every added feature. No truth is accessed.
"""
from __future__ import annotations

import numpy as np

from features import FEATURE_NAMES, QUALITY_NAMES


GEOMETRY_FEATURE_NAMES = (
    'geometry_direction_cosine', 'geometry_direction_sine',
    'geometry_momentum_balance', 'geometry_energy_balance',
    'geometry_reference_along_highPt', 'geometry_reference_along_lowPt',
    'geometry_reference_transverse_highPt', 'geometry_reference_transverse_lowPt',
    'geometry_straight_line_pair_dca',
    'geometry_vertex_transverse_highPt', 'geometry_vertex_transverse_lowPt',
    'geometry_pair_dca_over_vertex_error_scale',
)
PARALLEL_SINE_GUARD = 1e-8
MIN_ERROR_SCALE = 1e-12


def _raw_matrix(X):
    values = np.asarray(X, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != len(FEATURE_NAMES):
        raise ValueError('geometry wrapper requires raw reconstructed X with exactly50 columns')
    return np.where(np.isfinite(values), values, np.nan)


def _norm(vector):
    return np.hypot(np.hypot(vector[:, 0], vector[:, 1]), vector[:, 2])


def _ratio(numerator, denominator, minimum=0.):
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan),
                     where=np.isfinite(denominator) & (denominator > minimum))


def _derived(raw):
    columns = {name: raw[:, i] for i, name in enumerate(FEATURE_NAMES)}
    momentum, reference, direction, magnitude = {}, {}, {}, {}
    for side in ('highPt', 'lowPt'):
        momentum[side] = np.column_stack([columns['mu_' + side + '_' + key] for key in ('px', 'py', 'pz')])
        reference[side] = np.column_stack([columns['mu_' + side + '_' + key] for key in ('vx', 'vy', 'vz')])
        magnitude[side] = _norm(momentum[side])
        direction[side] = np.divide(momentum[side], magnitude[side][:, None],
            out=np.full_like(momentum[side], np.nan), where=(magnitude[side] > 0)[:, None])
    high, low = direction['highPt'], direction['lowPt']
    delta = reference['lowPt'] - reference['highPt']
    cross = np.cross(high, low)
    sine = np.clip(_norm(cross), 0, 1)
    cosine = np.clip(np.sum(high * low, axis=1), -1, 1)
    along_high, along_low = np.sum(delta * high, axis=1), np.sum(delta * low, axis=1)
    transverse_high, transverse_low = _norm(np.cross(delta, high)), _norm(np.cross(delta, low))
    # For exactly parallel lines their common perpendicular is the line distance.
    # Near the numerical guard, use the symmetric parallel approximation and
    # report its count. This is not the unconstrained nonparallel-line solution.
    parallel = np.isfinite(sine) & (sine <= PARALLEL_SINE_GUARD)
    straight_dca = _ratio(np.abs(np.sum(delta * cross, axis=1)), sine, PARALLEL_SINE_GUARD)
    straight_dca[parallel] = .5 * (transverse_high[parallel] + transverse_low[parallel])
    vertex = np.column_stack([columns['pair_' + key] for key in ('vx', 'vy', 'vz')])
    vertex_high = _norm(np.cross(vertex - reference['highPt'], high))
    vertex_low = _norm(np.cross(vertex - reference['lowPt'], low))
    ph, pl = magnitude['highPt'], magnitude['lowPt']
    momentum_balance = _ratio(np.minimum(ph, pl), np.maximum(ph, pl))
    eh, el = columns['mu_highPt_energy'], columns['mu_lowPt_energy']
    energy_balance = _ratio(np.minimum(eh, el), np.maximum(eh, el))
    energy_balance[(eh < 0) | (el < 0)] = np.nan
    error = np.column_stack([columns['pair_' + key + 'Err'] for key in ('vx', 'vy', 'vz')])
    error_scale = _norm(error)
    error_scale[(error < 0).any(axis=1)] = np.nan
    dca = columns['pair_dca'].copy()
    dca[dca < 0] = np.nan
    dca_scale = _ratio(dca, error_scale, MIN_ERROR_SCALE)
    derived = np.column_stack((cosine, sine, momentum_balance, energy_balance, along_high, along_low,
                               transverse_high, transverse_low, straight_dca, vertex_high, vertex_low, dca_scale))
    derived[~np.isfinite(derived)] = np.nan
    diagnostic = {
        'rows': len(raw), 'near_parallel_rows': int(parallel.sum()),
        'invalid_direction_rows': int((~np.isfinite(high).all(axis=1) | ~np.isfinite(low).all(axis=1)).sum()),
        'invalid_error_ratio_rows': int((~np.isfinite(dca_scale)).sum()),
        'missing_by_feature': {name: int((~np.isfinite(derived[:, i])).sum())
                               for i, name in enumerate(GEOMETRY_FEATURE_NAMES)},
    }
    return derived, diagnostic


def derive_geometry_features(X):
    """Return relative12 columns; invalid measurements stay missing for imputation."""
    return _derived(_raw_matrix(X))[0]


def geometry_feature_names(mode='all'):
    if mode not in ('all', 'quality_relative'):
        raise ValueError("geometry mode must be 'all' or 'quality_relative'")
    original = (FEATURE_NAMES if mode == 'all' else
                QUALITY_NAMES + ('pair_originCompatibilityNormalizedChi2',))
    return original + GEOMETRY_FEATURE_NAMES


def geometry_inputs(X, mode='all'):
    raw = _raw_matrix(X)
    names = geometry_feature_names(mode)
    original = [i for i, name in enumerate(FEATURE_NAMES) if name in names]
    return np.column_stack((raw[:, original], _derived(raw)[0]))


class GeometryVertexClassifier:
    """Raw50-only inference; reconstruction metadata/truth are never arguments."""

    def __init__(self, model, mode):
        self.model = model
        self.mode = mode
        self.n_features_in_ = len(FEATURE_NAMES)
        self.classes_ = np.asarray([0, 1])

    def predict_proba(self, X):
        return self.model.predict_proba(geometry_inputs(X, self.mode))


def train_geometry_vertex(X_train, y_train, w_train, Z_train,
                          X_val, y_val, w_val, Z_val, *, mode='all', **options):
    """Reuse unchanged train_neural_vertex and restore a raw50-only wrapper.

    The quality_relative mode uses26 quality features, precomputed origin
compatibility and12 relative features. Absolute momentum, energy and XYZ
columns are absent, but mass/displacement correlations can remain. Stored
origin compatibility is an input fit diagnostic, not recomputed after shifts.
Derived features are transformed/imputed/scaled by the train-only base trainer.
Energy regularization is forbidden because its canonical-X50 indexing would be
wrong after augmentation. Joint reco nuisances and category metadata retain the
base trainer's contracts and never enter wrapper inference.
"""
    if options.get('energy_strength', 0) != 0:
        raise ValueError('geometry augmentation requires energy_strength=0')
    from neural_vertex import train_neural_vertex
    train, val = _raw_matrix(X_train), _raw_matrix(X_val)
    model, history = train_neural_vertex(geometry_inputs(train, mode), y_train, w_train, Z_train,
                                         geometry_inputs(val, mode), y_val, w_val, Z_val, **options)
    history['geometry'] = {
        'mode': mode, 'input_features': list(FEATURE_NAMES), 'model_features': list(geometry_feature_names(mode)),
        'parallel_sine_guard': PARALLEL_SINE_GUARD, 'minimum_error_scale': MIN_ERROR_SCALE,
        'training': _derived(train)[1], 'validation': _derived(val)[1],
        'meaning': 'Approximate straight-line geometry; error scale ratios are not calibrated significances.',
        'parallel_guard': 'Symmetric mean transverse distance replaces ill-conditioned near-parallel line solution.',
        'quality_relative_caution': 'Absolute kinematic columns absent; correlations with mass/displacement can remain.',
    }
    return GeometryVertexClassifier(model, mode), history
