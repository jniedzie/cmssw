# Code style: AI-code. AI-generated working-point loss integrity tests.
import sys
import unittest
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from vertex_acceptance_loss import (acceptance_diagnostics, acceptance_memberships,
                                    build_acceptance_plan, soft_acceptance_penalty, training_cut)


class AcceptanceLossTest(unittest.TestCase):
    def test_group_construction_uses_only_training_distribution(self):
        y = np.r_[np.ones(90), np.zeros(30)].astype(int)
        nuisance = np.tile(np.arange(120)[:, None], (1, 3))
        category = np.full(120, 'Both-Both')
        process = np.full(120, 'dy')
        plan = build_acceptance_plan(y, nuisance, category, process)
        before = np.asarray(plan['edges']).copy()
        acceptance_memberships(plan, y, nuisance * 100, category, process)
        np.testing.assert_array_equal(plan['edges'], before)
        np.testing.assert_allclose(plan['edges'][0], [89 / 3, 178 / 3])

    def test_sparse_categories_are_not_targets(self):
        y = np.ones(90, dtype=int)
        z = np.tile(np.arange(90)[:, None], (1, 3))
        category = np.r_[np.full(85, 'Both-Both'), np.full(5, 'Other')]
        plan = build_acceptance_plan(y, z, category, np.full(90, 'dy'))
        self.assertFalse(any(g['scope'] == 'category' and g['value'] == 'Other' for g in plan['groups']))

    def test_uniform_group_shapes_have_zero_penalty(self):
        values = torch.tensor([.1, .4, .6, .9] * 2, requires_grad=True)
        loss, _ = soft_acceptance_penalty(values, np.arange(8), [('a', np.arange(4)), ('b', np.arange(4, 8))], .4)
        self.assertAlmostEqual(float(loss.detach()), 0., places=12)

    def test_nonuniform_groups_have_nonzero_gradient(self):
        values = torch.tensor([.2] * 4 + [.7] * 4, requires_grad=True)
        loss, details = soft_acceptance_penalty(values, np.arange(8), [('a', np.arange(4)), ('b', np.arange(4, 8))], .4)
        loss.backward()
        self.assertGreater(float(loss.detach()), .1)
        self.assertGreater(float(values.grad.abs().sum()), 0)
        self.assertEqual(details['groups_used'], 2)

    def test_accidental_scores_do_not_enter_genuine_acceptance_loss(self):
        values = torch.tensor([.2] * 4 + [.7] * 4 + [.99], requires_grad=True)
        loss, _ = soft_acceptance_penalty(values, np.arange(8), [('a', np.arange(4)), ('b', np.arange(4, 8))], .4)
        loss.backward()
        self.assertEqual(float(values.grad[-1]), 0.)

    def test_cut_ignores_fakes_and_exposes_soft_hard_difference(self):
        scores = np.array([.1, .2, .3, .4, .99])
        y = np.array([1, 1, 1, 1, 0])
        cut = training_cut(scores, y)
        self.assertAlmostEqual(cut, .1)
        detail = acceptance_diagnostics(scores, y, [('g', np.arange(4))], cut)
        self.assertEqual(detail['pooled_hard_acceptance'], 1.)
        self.assertNotEqual(detail['pooled_hard_acceptance'], detail['pooled_soft_acceptance'])

    def test_nonfinite_temperature_rejected(self):
        scores = torch.tensor([.2] * 4 + [.7] * 4)
        for temperature in (float('nan'), float('inf'), 0.):
            with self.assertRaisesRegex(ValueError, 'temperature'):
                soft_acceptance_penalty(scores, np.arange(8), [('a', np.arange(4))], .4, temperature=temperature)


if __name__ == '__main__':
    unittest.main()
