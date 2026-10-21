# Code style: AI-code. AI-generated analytic reference, not a production CMSSW fitter.
"""Iterated trajectory least squares with one fixed set of information factors.

The variables are the initial state and independent standard-normal process
kicks. Each transition is x[i+1] = f[i](x[i]) + L[i] w[i]. Rectangular L, including
zero columns, represents singular process noise without inventing a variance.
L and measurement errors are fixed for one fit. This is not a solver for a
parameter-dependent material likelihood. Raw track seeds are starting points;
only an independent prior belongs in Prior. Final covariance is the local
Gauss-Newton covariance, not an exact nonlinear posterior or momentum interval.
"""
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np


def immutable_array(value):
    result = np.array(value, dtype=float, copy=True)
    if not np.isfinite(result).all():
        raise ValueError("Nonfinite fit input")
    result.flags.writeable = False
    return result


@dataclass(frozen=True)
class Prior:
    mean: np.ndarray
    covariance: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, "mean", immutable_array(self.mean))
        object.__setattr__(self, "covariance", immutable_array(self.covariance))


@dataclass(frozen=True)
class Transition:
    predict: Callable
    jacobian: Callable
    noise_factor: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, "noise_factor", immutable_array(self.noise_factor))


@dataclass(frozen=True)
class Measurement:
    site: int
    value: np.ndarray
    covariance: np.ndarray
    predict: Callable
    jacobian: Callable

    def __post_init__(self):
        object.__setattr__(self, "value", immutable_array(self.value))
        object.__setattr__(self, "covariance", immutable_array(self.covariance))


@dataclass
class FitResult:
    valid: bool
    status: str
    parameters: np.ndarray
    states: tuple
    covariance: Optional[np.ndarray]
    state_jacobians: tuple
    objective: float
    iterations: int
    rank: int
    history: tuple
    evaluations: int = 0

    def state_covariance(self, first, second=None):
        if self.covariance is None:
            raise ValueError("No covariance for an unsuccessful fit")
        if second is None:
            second = first
        return self.state_jacobians[first] @ self.covariance @ self.state_jacobians[second].T


class TrajectoryProblem:
    """A frozen-noise generative model; callbacks must be deterministic and pure."""

    def __init__(self, dimension, transitions, measurements, prior=None):
        self.dimension = dimension
        self.transitions = tuple(transitions)
        self.measurements = tuple(measurements)
        self.prior = prior
        if dimension < 1 or not self.measurements:
            raise ValueError("A state and measurements are required")
        self.kick_slices = []
        count = dimension
        for transition in self.transitions:
            factor = transition.noise_factor
            if factor.ndim != 2 or factor.shape[0] != dimension:
                raise ValueError("Process noise factor has wrong dimensions")
            width = factor.shape[1]
            if width and np.linalg.matrix_rank(factor) != width:
                raise ValueError("Process noise columns must have their declared rank")
            self.kick_slices.append(slice(count, count + width))
            count += width
        self.parameter_count = count
        self.hit_cholesky = []
        for hit in self.measurements:
            if not 0 <= hit.site <= len(self.transitions):
                raise ValueError("Measurement site outside trajectory")
            if hit.value.ndim != 1 or hit.covariance.shape != (len(hit.value), len(hit.value)):
                raise ValueError("Measurement dimensions differ")
            self.hit_cholesky.append(self._cholesky(hit.covariance))
        self.prior_cholesky = None
        if prior is not None:
            if prior.mean.shape != (dimension,) or prior.covariance.shape != (dimension, dimension):
                raise ValueError("Prior dimensions differ")
            self.prior_cholesky = self._cholesky(prior.covariance)

    @staticmethod
    def _cholesky(covariance):
        if not np.allclose(covariance, covariance.T, rtol=0., atol=1.e-14):
            raise ValueError("Covariance must be symmetric")
        return np.linalg.cholesky(covariance)

    def initial_parameters(self, initial_state):
        state = np.asarray(initial_state, dtype=float)
        if state.shape != (self.dimension,) or not np.isfinite(state).all():
            raise ValueError("Invalid initial state")
        parameters = np.zeros(self.parameter_count)
        parameters[:self.dimension] = state
        return parameters

    def evaluate(self, parameters):
        """Return y-h, dh/dtheta and full states, including all affine offsets."""
        if parameters.shape != (self.parameter_count,) or not np.isfinite(parameters).all():
            raise ValueError("Invalid trajectory parameters")
        state = parameters[:self.dimension].copy()
        derivative = np.zeros((self.dimension, self.parameter_count))
        derivative[:, :self.dimension] = np.eye(self.dimension)
        states, derivatives = [state], [derivative]
        for transition, kick in zip(self.transitions, self.kick_slices):
            matrix = np.asarray(transition.jacobian(state), dtype=float)
            predicted = np.asarray(transition.predict(state), dtype=float)
            if matrix.shape != (self.dimension, self.dimension) or predicted.shape != state.shape:
                raise ValueError("Transition dimensions differ")
            state = predicted + transition.noise_factor @ parameters[kick]
            derivative = matrix @ derivative
            derivative[:, kick] += transition.noise_factor
            if not np.isfinite(state).all() or not np.isfinite(derivative).all():
                raise ValueError("Nonfinite propagation")
            states.append(state)
            derivatives.append(derivative)
        residuals, rows = [], []
        if self.prior is not None:
            residuals.append(np.linalg.solve(self.prior_cholesky, self.prior.mean - states[0]))
            rows.append(np.linalg.solve(self.prior_cholesky, derivatives[0]))
        for kick in self.kick_slices:
            if kick.stop == kick.start:
                continue
            row = np.zeros((kick.stop-kick.start, self.parameter_count))
            row[:, kick] = np.eye(kick.stop-kick.start)
            residuals.append(-parameters[kick])
            rows.append(row)
        for hit, chol in zip(self.measurements, self.hit_cholesky):
            predicted = np.asarray(hit.predict(states[hit.site]), dtype=float)
            matrix = np.asarray(hit.jacobian(states[hit.site]), dtype=float)
            if predicted.shape != hit.value.shape or matrix.shape != (len(hit.value), self.dimension):
                raise ValueError("Measurement model dimensions differ")
            residuals.append(np.linalg.solve(chol, hit.value - predicted))
            rows.append(np.linalg.solve(chol, matrix @ derivatives[hit.site]))
        residual = np.concatenate(residuals)
        jacobian = np.vstack(rows)
        if not np.isfinite(residual).all() or not np.isfinite(jacobian).all():
            raise ValueError("Nonfinite measurement model")
        return residual, jacobian, tuple(states), tuple(derivatives)


def stationary_curvature(problem, parameters, scale, evaluate):
    """Reject clear negative curvature; never turn damping into fit information.

    This bounded reference check rejects clear negative curvature at its probe
    resolution, without certifying a minimum. It differentiates the objective gradient,
    independently of the positive Gauss-Newton matrix. It is not a production
    prescription to perform O(parameter_count) full Geant4 trajectories.
    """
    matrices, roundoff_bounds = [], []
    for step in (1.e-4, 5.e-5):
        matrix = np.empty((len(parameters), len(parameters)))
        entry_roundoff = 0.
        for column in range(len(parameters)):
            offset = np.zeros(len(parameters))
            offset[column] = step * scale[column]
            plus, minus = parameters + offset, parameters - offset
            if plus[column] == parameters[column] or minus[column] == parameters[column]:
                return "curvature_unverified", dict(lost_curvature_probe_coordinate=column)
            plus_r, plus_j, _, _ = evaluate(plus)
            minus_r, minus_j, _, _ = evaluate(minus)
            plus_g = -scale * (plus_j.T @ plus_r)
            minus_g = -scale * (minus_j.T @ minus_r)
            if not np.isfinite(plus_g).all() or not np.isfinite(minus_g).all():
                return "curvature_unverified", dict(nonfinite_curvature_gradient=True)
            matrix[:, column] = (plus_g - minus_g) / (2. * step)
            # gamma_m bounds the floating-point gradient dot products. This
            # does not certify inaccuracies inside user-supplied callbacks.
            epsilon = np.finfo(float).eps
            gamma = len(plus_r) * epsilon / (1. - len(plus_r) * epsilon)
            plus_sum = float(np.max(np.abs(plus_j * scale).T @ np.abs(plus_r)))
            minus_sum = float(np.max(np.abs(minus_j * scale).T @ np.abs(minus_r)))
            difference_roundoff = epsilon * float(np.max(np.abs(plus_g) + np.abs(minus_g)))
            bound = (gamma * (plus_sum + minus_sum) + difference_roundoff) / (2. * step)
            entry_roundoff = max(entry_roundoff, bound)
        if not np.isfinite(matrix).all() or not np.isfinite(entry_roundoff):
            return "curvature_unverified", dict(nonfinite_curvature_probe=True)
        matrices.append(matrix)
        roundoff_bounds.append(entry_roundoff)
    first, second = matrices
    difference = float(np.linalg.norm(first-second, ord=2))
    # Entry error b implies operator error <= p*b, and skew error <= 2*p*b.
    roundoff = 2. * len(parameters) * max(roundoff_bounds)
    error = 4. * difference + roundoff + 100. * np.finfo(float).eps * max(
        1., float(np.linalg.norm(second, ord=2)))
    skew = float(np.linalg.norm(second-second.T, ord=2))
    if skew > error:
        return "curvature_unverified", dict(curvature_skew=skew, curvature_probe_error=error,
                                             curvature_roundoff_estimate=roundoff,
                                             curvature_two_step_difference=difference)
    minimum = float(np.linalg.eigvalsh(0.5 * (second+second.T))[0])
    diagnostics = dict(minimum_objective_curvature=minimum, curvature_probe_error=error,
                       curvature_roundoff_estimate=roundoff, curvature_two_step_difference=difference)
    return ("negative_curvature" if minimum < -error else "nonnegative_curvature"), diagnostics


class EvaluationLimit(Exception):
    pass


def fit_trajectory(problem, initial_state, maximum_iterations=40,
                   gradient_tolerance=1.e-9, maximum_backtracks=24,
                   maximum_evaluations=1000):
    """Relinearize one fixed objective, without posterior covariance feedback.

    Accepted stationarity includes a local curvature check. It is not a global
    minimum guarantee; final covariance is the local undamped GN approximation.
    """
    if (maximum_iterations < 0 or maximum_backtracks < 1 or maximum_evaluations < 1 or
            not np.isfinite(gradient_tolerance) or gradient_tolerance <= 0.):
        raise ValueError("Invalid numerical work limits")
    parameters = problem.initial_parameters(initial_state)
    history = []
    states, derivatives = (), ()
    objective, rank, evaluations = float("inf"), 0, 0

    def evaluate(parameters):
        nonlocal evaluations
        if evaluations == maximum_evaluations:
            raise EvaluationLimit()
        evaluations += 1
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            return problem.evaluate(parameters)

    def result(valid, status, iteration, covariance=None):
        return FitResult(valid, status, parameters, states, covariance,
                         derivatives, objective, iteration, rank, tuple(history), evaluations)

    for iteration in range(maximum_iterations + 1):
        try:
            residual, jacobian, states, derivatives = evaluate(parameters)
        except EvaluationLimit:
            return result(False, "evaluation_limit", iteration)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            return result(False, "invalid_model", iteration)
        with np.errstate(over="ignore", invalid="ignore"):
            objective = 0.5 * float(residual @ residual)
            norms = np.linalg.norm(jacobian, axis=0)
        if not np.isfinite(objective):
            return result(False, "invalid_objective", iteration)
        if np.any(norms == 0.) or not np.isfinite(norms).all():
            return result(False, "unobservable", iteration)
        # Column scaling changes numerical units, never the information factors.
        with np.errstate(over="ignore"):
            scale = 1. / norms
        if not np.isfinite(scale).all():
            return result(False, "unobservable", iteration)
        matrix = jacobian * scale
        try:
            left, singular, right = np.linalg.svd(matrix, full_matrices=False)
        except np.linalg.LinAlgError:
            return result(False, "singular_linearization", iteration)
        threshold = np.finfo(float).eps * max(matrix.shape) * singular[0]
        rank = int(np.count_nonzero(singular > threshold))
        if rank < problem.parameter_count:
            return result(False, "unobservable", iteration)
        projected_residual = left.T @ residual
        gradient = float(np.max(np.abs(matrix.T @ residual)))
        decrement = float(np.linalg.norm(projected_residual))
        if not np.isfinite(gradient) or not np.isfinite(decrement):
            return result(False, "invalid_linearization", iteration)
        history.append(dict(iteration=iteration, objective=objective,
                            scaled_gradient=gradient, linearized_decrement=decrement, rank=rank))
        residual_norm = float(np.linalg.norm(residual))
        tolerance = gradient_tolerance * max(1., residual_norm)
        # Below this decrement the predicted cost reduction is smaller than
        # machine epsilon times the current cost. Label this separately; it
        # does not establish parameter accuracy or calibrated covariance.
        roundoff_tolerance = np.sqrt(np.finfo(float).eps) * residual_norm
        decrement_tolerance = max(tolerance, roundoff_tolerance)
        # A small normal-equation gradient alone can miss a large step in a
        # weak but observable direction. The SVD decrement catches that case.
        if gradient <= tolerance and decrement <= decrement_tolerance:
            history[-1]["convergence_reason"] = ("user_tolerance" if decrement <= tolerance
                                                   else "objective_roundoff")
            try:
                status, diagnostics = stationary_curvature(problem, parameters, scale, evaluate)
            except EvaluationLimit:
                return result(False, "evaluation_limit", iteration)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                return result(False, "curvature_unverified", iteration)
            history[-1].update(diagnostics)
            if status != "nonnegative_curvature":
                return result(False, status, iteration)
            inverse = (scale[:, None] * right.T) / singular
            with np.errstate(over="ignore", invalid="ignore"):
                covariance = inverse @ inverse.T
            if not np.isfinite(covariance).all():
                return result(False, "unobservable", iteration)
            return result(True, "converged", iteration, covariance)
        if iteration == maximum_iterations:
            break
        step = scale * (right.T @ (projected_residual / singular))
        descent = float(residual @ (jacobian @ step))
        accepted = False
        for backtrack in range(maximum_backtracks):
            fraction = 0.5 ** backtrack
            trial = parameters + fraction * step
            try:
                trial_evaluation = evaluate(trial)
                trial_residual = trial_evaluation[0]
                trial_objective = 0.5 * float(trial_residual @ trial_residual)
            except EvaluationLimit:
                return result(False, "evaluation_limit", iteration)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                continue
            if trial_objective <= objective - 1.e-4 * fraction * descent:
                history[-1]["step_fraction"] = fraction
                parameters = trial
                states, derivatives = trial_evaluation[2:]
                objective = trial_objective
                accepted = True
                break
        if not accepted:
            return result(False, "no_descent", iteration)
    return result(False, "iteration_limit", maximum_iterations)
