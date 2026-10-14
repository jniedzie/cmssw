# Code style: AI-code. AI-generated category-network boundary tests.
from pathlib import Path
import pickle
import sys
import unittest

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from neural_vertex_categories import CategoryNetwork, conditional_groups, conditional_penalty, train_category_vertex


def fixture():
    rng = np.random.default_rng(39)
    x = rng.normal(size=(384, 7))
    y = np.tile([0, 1], 192)
    category = np.tile(['A', 'A', 'B', 'B'], 96)
    process = np.tile(['dy']*4 + ['jpsi']*4, 48)
    x[:, 0] += 2*y-1
    x[:, -1] = np.nan
    z = x[:, :3].copy()
    w = rng.uniform(.2, 2, len(y))
    arrays = [x[:256], y[:256], w[:256], z[:256], x[256:], y[256:], w[256:], z[256:]]
    options = dict(feature_names=['reco_'+str(i) for i in range(7)], train_category=category[:256],
                   val_category=category[256:], train_process=process[:256], val_process=process[256:],
                   hidden=(8,), batch_size=64, epochs=4, seed=12)
    return arrays, options


class CategoryTraining(unittest.TestCase):
    def test_heads_route_reco_category_and_start_as_shared_logit(self):
        network = CategoryNetwork(2, 2, 'heads', (3,), 0)
        x = torch.ones((2, 2)); categories = torch.eye(2)
        np.testing.assert_array_equal(network(x, categories).detach().numpy(), network(x, categories.flip(0)).detach().numpy())
        with torch.no_grad():
            network.heads.bias[:] = torch.tensor([-2., 2.])
        self.assertAlmostEqual(float((network(x, categories)[1]-network(x, categories)[0]).detach()), 4.)

    def test_zero_penalties_ignore_nuisance_and_process_not_category(self):
        data, options = fixture()
        first, h1 = train_category_vertex(*data, **options, strength=0, category_strength=0, conditional_strength=0)
        changed = [a.copy() for a in data]; changed[3][:] = np.nan; changed[7] *= 10000
        options['train_process'] = np.full(256, 'changed'); options['val_process'] = np.full(128, 'different')
        second, h2 = train_category_vertex(*changed, **options, strength=0, category_strength=0, conditional_strength=0)
        np.testing.assert_array_equal(first.predict_proba(data[4], options['val_category']), second.predict_proba(data[4], options['val_category']))
        self.assertEqual(h1['best_state_sha256'], h2['best_state_sha256'])

    def test_checkpoint_restoration_train_only_preprocessing_and_pickle(self):
        data, options = fixture(); data[5] = 1-data[5]
        options.update(epochs=25, patience=3, learning_rate=.03, dropout=0, architecture='heads')
        model, history = train_category_vertex(*data, **options, strength=0)
        self.assertLess(history['best_epoch'], history['epochs_completed'])
        self.assertEqual(history['best_state_sha256'], history['restored_state_sha256'])
        self.assertNotEqual(history['last_state_sha256'], history['restored_state_sha256'])
        restored = pickle.loads(pickle.dumps(model))
        np.testing.assert_array_equal(model.predict_proba(data[4], options['val_category']), restored.predict_proba(data[4], options['val_category']))
        self.assertEqual(model.predict_proba(np.empty((0, 7)), np.array([], dtype=str)).shape, (0, 2))
        self.assertFalse(hasattr(model, 'train_process'))
        self.assertFalse(hasattr(model, 'nuisance_scaling'))
        changed = [a.copy() for a in data];changed[4][:, :-1] += 1000;changed[7] *= 1000
        other, h2 = train_category_vertex(*changed, **options, strength=0)
        np.testing.assert_array_equal(model.imputer.statistics_, other.imputer.statistics_)
        np.testing.assert_array_equal(model.scaler.mean_, other.scaler.mean_)
        self.assertEqual(history['nuisance_scaling'], h2['nuisance_scaling'])

    def test_conditional_support_uses_training_and_penalty_has_joint_gradient(self):
        data, options = fixture()
        groups, specs = conditional_groups(data[1], data[2], data[3], options['train_category'], options['train_process'])
        self.assertEqual(len(specs), 8)
        score = torch.linspace(.01, .99, 256, dtype=torch.float64, requires_grad=True)
        value, details = conditional_penalty(score, torch.tensor(data[3]), torch.tensor(data[2]), groups)
        value.backward()
        self.assertGreater(float(score.grad.abs().sum()), 0)
        self.assertTrue(torch.isfinite(score.grad).all())
        empty, _ = conditional_groups(data[5], data[6], data[7], np.full(128, 'untrained'), options['val_process'], specs, minimum=0)
        self.assertEqual(len(empty), 8)
        self.assertTrue(all(len(idx) == 0 for _, idx in empty))

    def test_full_validation_objective_during_warmup_and_support_metadata(self):
        data, options = fixture()
        _, history = train_category_vertex(*data, **options, architecture='heads', strength=.2,
                                           category_strength=.3, conditional_strength=.4, warmup_epochs=2)
        self.assertEqual([r['optimization_penalty_fraction'] for r in history['epochs']], [0., .5, 1., 1.])
        for row in history['epochs']:
            expected = row['val_bce'] + .2*row['val_dcor_squared'] + .3*row['val_category_dcor_squared'] + .4*row['val_conditional_dcor_squared']
            self.assertAlmostEqual(row['val_objective'], expected, places=12)
        self.assertEqual(history['parameter_counts']['per_category_residual_head'], 9)
        self.assertEqual(history['category_support']['A']['train']['genuine'], 64)

    def test_acceptance_cut_and_plan_use_training_and_full_checkpoint_loss(self):
        from vertex_acceptance_loss import training_cut, acceptance_memberships, soft_acceptance_penalty
        data, options = fixture()
        model, history = train_category_vertex(*data, **options, strength=.2, category_strength=.3,
                                               acceptance_strength=10, warmup_epochs=2)
        best = history['epochs'][history['best_epoch']-1]
        score = model.predict_proba(data[0], options['train_category'])[:, 1]
        self.assertEqual(best['checkpoint_training_cut'], training_cut(score, data[1]))
        self.assertFalse(hasattr(model, 'acceptance_plan'))
        memberships = acceptance_memberships(history['acceptance_plan'], data[5], data[7], options['val_category'], options['val_process'])
        val_score = torch.from_numpy(model.predict_proba(data[4], options['val_category'])[:, 1])
        penalty, _ = soft_acceptance_penalty(val_score, np.flatnonzero(data[5] == 1), memberships, best['checkpoint_training_cut'])
        self.assertEqual(float(penalty), best['val_soft_acceptance_penalty'])
        self.assertGreater(best['val_soft_acceptance_penalty'], 0)
        for row in history['epochs']:
            expected = row['val_bce'] + .2*row['val_dcor_squared'] + .3*row['val_category_dcor_squared'] + 10*row['val_soft_acceptance_penalty']
            self.assertAlmostEqual(row['val_objective'], expected, places=12)
        changed = [a.copy() for a in data];changed[7] += 10000
        _, other = train_category_vertex(*changed, **options, strength=.2, category_strength=.3,
                                        acceptance_strength=10, warmup_epochs=2)
        self.assertEqual(history['acceptance_plan'], other['acceptance_plan'])
        self.assertEqual([r['state_sha256'] for r in history['epochs']], [r['state_sha256'] for r in other['epochs']])
        self.assertEqual([r['optimization_training_cut'] for r in history['epochs']], [r['optimization_training_cut'] for r in other['epochs']])
        self.assertNotEqual(best['val_soft_acceptance_penalty'], other['epochs'][history['best_epoch']-1]['val_soft_acceptance_penalty'])

    def test_unseen_category_truth_feature_and_bad_weights_fail(self):
        data, options = fixture()
        model, _ = train_category_vertex(*data, **options, strength=0)
        with self.assertRaisesRegex(ValueError, 'Unseen'):
            model.predict_proba(data[4], np.full(128, 'new'))
        for name in ['ShiftMuon_genPartIdx', 'samplingWeight', 'weight', 'sample', 'source_index', 'pair_index']:
            options['feature_names'][0] = name
            with self.assertRaisesRegex(ValueError, 'Non-reconstructed'):
                train_category_vertex(*data, **options)
        data, options = fixture();data[2][0] = 0
        with self.assertRaisesRegex(ValueError, 'strictly positive'):
            train_category_vertex(*data, **options, acceptance_strength=1)
        data, options = fixture();data[2][0] = -1
        with self.assertRaisesRegex(ValueError, 'nonnegative'):
            train_category_vertex(*data, **options)


if __name__ == '__main__':
    unittest.main()
