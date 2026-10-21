# Code style: AI-code. AI-generated analytic closure tests for a reference trajectory fit.
"""These verify numerical/information contracts, not SHIFT resolution or Geant4."""
import argparse
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

from iteratedTrajectoryFit import Measurement, Prior, Transition, TrajectoryProblem, fit_trajectory

RESULTS = {}


def linear_hit(site, matrix, value, covariance):
    matrix = np.asarray(matrix, dtype=float)
    return Measurement(site, value, covariance, lambda state: matrix @ state,
                       lambda state: matrix)


def linear_transition(matrix, offset, factor):
    matrix, offset = np.asarray(matrix, dtype=float), np.asarray(offset, dtype=float)
    return Transition(lambda state: matrix @ state + offset, lambda state: matrix, factor)


def nonlinear_charge_problem(charge):
    """Synthetic bending, deterministic mean loss and fixed angular scattering.

    State is (log(abs(q/p)), direction slope, position); charge is a discrete
    hypothesis. Numbers define an analytic fixture, not a CMS material model.
    """
    distance, bending, loss = 1.0, 0.8, 0.2

    def predict(state):
        inverse_p = np.exp(state[0])
        if loss * inverse_p >= 1.:
            raise ValueError("Particle cannot traverse this fixture's material")
        middle = charge * inverse_p / (1. - 0.5 * loss * inverse_p)
        return np.array([state[0] - np.log1p(-loss * inverse_p),
                         state[1] + bending * distance * middle,
                         state[2] + distance * state[1] + 0.5 * bending * distance**2 * middle])

    def jacobian(state):
        inverse_p = np.exp(state[0])
        if loss * inverse_p >= 1.:
            raise ValueError("Particle cannot traverse this fixture's material")
        curvature_derivative = charge * inverse_p / (1. - 0.5 * loss * inverse_p)**2
        return np.array([[1. / (1. - loss * inverse_p), 0., 0.],
                         [bending * distance * curvature_derivative, 1., 0.],
                         [0.5 * bending * distance**2 * curvature_derivative, distance, 1.]])

    transition = Transition(predict, jacobian, [[0.], [0.02], [0.]])
    # Measurements are literal inputs; no true parameter is passed to the fit.
    positions = [0.02, 0.14, 0.43, 0.89, 1.52, 2.34]
    hits = [linear_hit(i, [[0., 0., 1.]], [charge * position], [[0.01**2]])
            for i, position in enumerate(positions)]
    return TrajectoryProblem(3, [transition] * (len(hits)-1), hits)


class IteratedTrajectoryTest(unittest.TestCase):
    def close(self, observed, expected, **options):
        np.testing.assert_allclose(observed, expected, rtol=2.e-9, atol=2.e-10, **options)

    def record(self, name, fit):
        RESULTS[name] = dict(valid=fit.valid, status=fit.status, parameters=fit.parameters.tolist(),
                             states=[s.tolist() for s in fit.states],
                             covariance=fit.covariance.tolist() if fit.covariance is not None else None,
                             objective=fit.objective, iterations=fit.iterations,
                             rank=fit.rank, evaluations=fit.evaluations, history=fit.history)

    def test_fixed_prior_is_not_reused_as_posterior(self):
        problem = TrajectoryProblem(1, [], [linear_hit(0, [[1.]], [2.], [[1.]])], Prior([0.], [[4.]]))
        start = [0.]
        for _ in range(4):
            fit = fit_trajectory(problem, start)
            self.assertTrue(fit.valid, fit.status)
            self.close(fit.states[0], [8./5.])
            self.close(fit.covariance, [[4./5.]])
            start = fit.states[0]
        feedback = TrajectoryProblem(1, [], problem.measurements, Prior(fit.states[0], fit.covariance))
        wrong = fit_trajectory(feedback, start)
        self.close(wrong.states[0], [16./9.])
        self.close(wrong.covariance, [[4./9.]])
        self.record("fixed_prior", fit)
        self.record("deliberate_posterior_feedback_negative_control", wrong)

    def test_affine_transport_and_process_noise(self):
        process = linear_transition([[2.]], [3.], [[3.]])
        hit = linear_hit(1, [[1.]], [11.], [[1.]])
        problem = TrajectoryProblem(1, [process], [hit], Prior([1.], [[4.]]))
        fit = fit_trajectory(problem, [5.])
        self.assertTrue(fit.valid, fit.status)
        self.close(fit.states, [[37./13.], [140./13.]])
        self.close(fit.state_covariance(0), [[20./13.]])
        self.close(fit.state_covariance(1), [[25./26.]])
        self.close(fit.state_covariance(0, 1), [[4./13.]])
        deterministic = TrajectoryProblem(1, [linear_transition([[2.]], [3.], np.empty((1, 0)))],
                                          [hit], problem.prior)
        omitted = fit_trajectory(deterministic, [5.])
        self.assertTrue(omitted.valid, omitted.status)
        self.close(omitted.state_covariance(0), [[4./17.]])
        self.record("affine_stochastic_transport", fit)
        self.record("deliberate_process_noise_omission_negative_control", omitted)

    def test_rank_one_noise_and_correlated_prior(self):
        prior = Prior([0., 1.], [[4., 1.], [1., 2.]])
        hit = linear_hit(1, [[1., 0.]], [3.], [[1.]])
        factor = np.array([[0.5], [1.]])
        fits = []
        for sign in [1., -1.]:
            process = linear_transition([[1., 1.], [0., 1.]], [0., 0.], sign * factor)
            fit = fit_trajectory(TrajectoryProblem(2, [process], [hit], prior), [-5., 6.])
            self.assertTrue(fit.valid, fit.status)
            self.close(fit.states[0], [40./37., 61./37.])
            self.close(fit.states[1], [103./37., 65./37.])
            self.close(fit.state_covariance(1), np.array([[33., 14.], [14., 62.]])/37.)
            self.close(fit.parameters[2:], [sign * 4./37.])
            # Every process realization lies in the rank-one noise support.
            defect = fit.states[1] - np.array([[1., 1.], [0., 1.]]) @ fit.states[0]
            self.close(defect[0], 0.5 * defect[1])
            fits.append(fit)
        self.close(fits[0].covariance, np.array([[48., -23., -10.], [-23., 38., -6.], [-10., -6., 36.]])/37.)
        self.close(fits[0].state_covariance(0, 1), fits[1].state_covariance(0, 1))
        self.record("singular_process_noise", fits[0])

    def test_starting_guess_is_not_an_independent_prior(self):
        hits = [linear_hit(0, [[1.]], [3.], [[4.]]), linear_hit(0, [[1.]], [5.], [[1.]])]
        fits = [fit_trajectory(TrajectoryProblem(1, [], hits), [start]) for start in [-100., 100.]]
        for fit in fits:
            self.assertTrue(fit.valid, fit.status)
            self.close(fit.states[0], [23./5.])
            self.close(fit.covariance, [[4./5.]])
        self.record("no_prior_seed_independence", fits[0])

    def test_nonlinear_affine_offset_and_covariance(self):
        hit = Measurement(0, [4.], [[1.]], lambda x: x*x, lambda x: np.array([[2.*x[0]]]))
        problem = TrajectoryProblem(1, [], [hit], Prior([1.], [[1.]]))
        first = fit_trajectory(problem, [2.], maximum_iterations=1)
        self.assertEqual(first.status, "iteration_limit")
        self.close(first.states[0], [33./17.])
        roots = np.roots([2., 0., -7., -1.])
        expected, = [x.real for x in roots if abs(x.imag) < 1.e-12 and x.real > 0.]
        fits = [fit_trajectory(problem, [start]) for start in [0.4, 3.]]
        for fit in fits:
            self.assertTrue(fit.valid, fit.status)
            self.close(fit.states[0], [expected])
            self.close(fit.covariance, [[1./(1.+4.*expected**2)]])
            costs = [h["objective"] for h in fit.history]
            self.assertTrue(all(b <= a for a,b in zip(costs, costs[1:])))
        self.record("nonlinear_fixed_objective", fits[0])

    def test_nonlinear_deterministic_transport(self):
        process = Transition(lambda x: x*x, lambda x: np.array([[2.*x[0]]]), np.empty((1, 0)))
        problem = TrajectoryProblem(1, [process], [linear_hit(1, [[1.]], [4.], [[1.]])], Prior([1.], [[1.]]))
        fit = fit_trajectory(problem, [2.])
        self.assertTrue(fit.valid, fit.status)
        self.close(fit.states[1], fit.states[0]**2)
        self.record("nonlinear_zero_noise_transport", fit)

    def test_reject_unobservable_and_stationary_maximum(self):
        unobserved = TrajectoryProblem(2, [], [linear_hit(0, [[1., 0.]], [1.], [[1.]])])
        fit = fit_trajectory(unobserved, [0., 0.])
        self.assertFalse(fit.valid)
        self.assertEqual(fit.status, "unobservable")
        self.assertIsNone(fit.covariance)
        hit = Measurement(0, [0., 1.], np.eye(2), lambda x: np.array([x[0], x[0]**2]),
                          lambda x: np.array([[1.], [2.*x[0]]]))
        saddle = fit_trajectory(TrajectoryProblem(1, [], [hit]), [0.])
        self.assertFalse(saddle.valid)
        self.assertEqual(saddle.status, "negative_curvature")
        self.assertIsNone(saddle.covariance)
        self.record("stationary_maximum_rejected", saddle)

    def test_weak_direction_does_not_falsely_converge(self):
        matrix = np.array([[1., 1.], [0., 1.e-10]])
        problem = TrajectoryProblem(2, [], [linear_hit(0, matrix, [0., 1.], np.eye(2))])
        fit = fit_trajectory(problem, [0., 0.])
        self.assertTrue(fit.valid, fit.status)
        self.assertGreater(fit.iterations, 0)
        self.close(matrix @ fit.states[0], [0., 1.])
        self.record("weak_observable_direction", fit)

    def test_failure_budget_preserves_returned_state(self):
        problem = TrajectoryProblem(1, [], [linear_hit(0, [[1.]], [2.], [[1.]])])
        for budget in [1, 2]:
            fit = fit_trajectory(problem, [0.], maximum_evaluations=budget)
            self.assertFalse(fit.valid)
            self.assertEqual(fit.status, "evaluation_limit")
            self.assertIsNone(fit.covariance)
            self.assertEqual(fit.evaluations, budget)
            if fit.states:
                self.close(fit.parameters[:1], fit.states[0])
                self.close(fit.objective, 0.5 * (2.-fit.states[0][0])**2)
            self.record("evaluation_limit_"+str(budget), fit)

    def test_unrepresentable_curvature_probe_and_overflow(self):
        large = TrajectoryProblem(1, [], [linear_hit(0, [[1.]], [1.e20], [[1.]])])
        fit = fit_trajectory(large, [1.e20])
        self.assertFalse(fit.valid)
        self.assertEqual(fit.status, "curvature_unverified")
        self.assertIsNone(fit.covariance)
        overflow = TrajectoryProblem(1, [], [linear_hit(0, [[1.]], [1.e200], [[1.]])])
        invalid = fit_trajectory(overflow, [0.])
        self.assertFalse(invalid.valid)
        self.assertEqual(invalid.status, "invalid_objective")
        self.assertIsNone(invalid.covariance)

    def test_fixed_charge_nonlinear_bending_and_loss(self):
        fits = {}
        for charge in [-1., 1.]:
            problem = nonlinear_charge_problem(charge)
            starts = [[np.log(0.1), charge*0.2, charge*0.1],
                      [np.log(0.3), charge*-0.1, charge*-0.1]]
            for index, initial in enumerate(starts):
                fit = fit_trajectory(problem, initial)
                self.assertTrue(fit.valid, fit.status)
                self.assertTrue(all(np.isfinite(s).all() for s in fit.states))
                self.assertTrue(all(0.2*np.exp(s[0]) < 1. for s in fit.states))
                self.assertTrue(all(b["objective"] <= a["objective"] for a,b in zip(fit.history, fit.history[1:])))
                fits[charge, index] = fit
                self.record("nonlinear_charge_"+str(charge)+"_start"+str(index), fit)
        for charge in [-1., 1.]:
            delta = fits[charge, 0].parameters - fits[charge, 1].parameters
            standardized = np.linalg.solve(np.linalg.cholesky(fits[charge, 0].covariance), delta)
            # The stationarity floor concerns distinguishable objective, not
            # an arbitrary absolute tolerance on differently scaled states.
            self.assertLess(np.linalg.norm(standardized), 1.e-7)
            RESULTS["nonlinear_charge_"+str(charge)+"_start0"]["alternate_start_GN_standardized_distance"] = float(np.linalg.norm(standardized))
            self.close(fits[charge, 0].covariance, fits[charge, 1].covariance)
        self.close([s[0] for s in fits[-1., 0].states], [s[0] for s in fits[1., 0].states])
        self.close([s[1:] for s in fits[-1., 0].states], -np.array([s[1:] for s in fits[1., 0].states]))
        self.close(fits[-1., 0].objective, fits[1., 0].objective)
        # Independent central differences test the actual model Jacobian.
        problem = nonlinear_charge_problem(1.)
        state = fits[1., 0].states[0]
        observed = np.column_stack([(problem.transitions[0].predict(state+np.eye(3)[i]*1.e-6)-
                                     problem.transitions[0].predict(state-np.eye(3)[i]*1.e-6))/(2.e-6)
                                    for i in range(3)])
        self.close(problem.transitions[0].jacobian(state), observed)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(IteratedTrajectoryTest)
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    if args.receipt:
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        sources = [Path(__file__), Path(__file__).with_name("iteratedTrajectoryFit.py")]
        receipt = dict(passed=outcome.wasSuccessful(), tests=outcome.testsRun,
                       failures=len(outcome.failures), errors=len(outcome.errors),
                       native_events=0, synthetic_only=True, results=RESULTS,
                       source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
        args.receipt.write_text(json.dumps(receipt, indent=2, allow_nan=False)+"\n")
    raise SystemExit(not outcome.wasSuccessful())
