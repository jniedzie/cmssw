"""Count-only identification bounds against exhaustive unknown assignments."""

from itertools import product
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from acceptance_audit import unknown_bounds


def assignment_rates(y, passed):
    """Independent enumeration, omitting rates whose class is absent.

    A rate is undefined in an assignment with no members of that class. Such
    assignments cannot supply a numerical extremum. If every assignment lacks
    the class, both extrema are undefined. No NaN placeholder participates.
    """
    unknown = np.flatnonzero(y < 0)
    genuine, rejection = [], []
    for values in product((0, 1), repeat=len(unknown)):
        assigned = y.copy()
        assigned[unknown] = values
        good = assigned == 1
        bad = assigned == 0
        if good.any():
            genuine.append(float(passed[good].mean()))
        if bad.any():
            rejection.append(float((~passed[bad]).mean()))
    return genuine, rejection


def extrema(values):
    return [min(values), max(values)] if values else [None, None]


class UnknownIdentificationBounds(unittest.TestCase):
    def test_extrema_equal_every_feasible_assignment_for_vectors_up_to_four_rows(self):
        # 1555 observed-label/decision patterns, including all empty-class cases.
        checked = 0
        for size in range(5):
            for labels in product((-1, 0, 1), repeat=size):
                y = np.array(labels, dtype=np.int8)
                for decisions in product((False, True), repeat=size):
                    passed = np.array(decisions, dtype=bool)
                    genuine, rejection = assignment_rates(y, passed)
                    result = unknown_bounds(y, passed)
                    context = f"labels={labels}, passed={decisions}"
                    self.assertEqual(result["genuine_efficiency_extrema"], extrema(genuine), context)
                    self.assertEqual(result["accidental_rejection_extrema"], extrema(rejection), context)
                    # Endpoints come from actual assignments, rather than a
                    # hypothetical zero-denominator corner or extrapolation.
                    for key, rates in (("genuine_efficiency_extrema", genuine),
                                       ("accidental_rejection_extrema", rejection)):
                        for bound in result[key]:
                            if bound is not None:
                                self.assertIn(bound, rates, context)
                                self.assertTrue(np.isfinite(bound), context)
                    checked += 1
        self.assertEqual(checked, 1555)

    def test_no_unknowns_collapse_to_observed_class_fractions(self):
        result = unknown_bounds(np.array([1, 1, 0, 0, 0]),
                                np.array([True, False, True, False, False]))
        self.assertEqual(result["genuine_efficiency_extrema"], [0.5, 0.5])
        self.assertEqual(result["accidental_rejection_extrema"], [2 / 3, 2 / 3])
        self.assertEqual(result["unknown_passed"], 0)
        self.assertEqual(result["unknown_rejected"], 0)

    def test_class_absent_from_every_assignment_has_undefined_extrema(self):
        cases = [([], [], [None, None], [None, None]),
                 ([0, 0], [True, False], [None, None], [0.5, 0.5]),
                 ([1, 1], [True, False], [0.5, 0.5], [None, None])]
        for labels, decisions, genuine, rejection in cases:
            with self.subTest(labels=labels, decisions=decisions):
                result = unknown_bounds(np.array(labels, dtype=np.int8), np.array(decisions, dtype=bool))
                self.assertEqual(result["genuine_efficiency_extrema"], genuine)
                self.assertEqual(result["accidental_rejection_extrema"], rejection)

    def test_one_sided_unknowns_have_finite_rates_whenever_the_class_exists(self):
        for decision in (False, True):
            with self.subTest(passed=decision):
                result = unknown_bounds(np.array([-1]), np.array([decision]))
                self.assertEqual(result["genuine_efficiency_extrema"], [float(decision)] * 2)
                self.assertEqual(result["accidental_rejection_extrema"], [float(not decision)] * 2)

    def test_known_and_unknown_counts_are_retained_without_input_mutation(self):
        y = np.array([1, 1, 0, 0, -1, -1], dtype=np.int8)
        score = np.array([0.9, 0.1, 0.8, 0.2, 0.7, 0.3])
        passed = score >= 0.5
        before_y, before_score, before_passed = y.copy(), score.copy(), passed.copy()
        result = unknown_bounds(y, passed)
        self.assertEqual({key: result[key] for key in (
            "known_common", "known_common_passed", "known_accidental", "known_accidental_rejected",
            "unknown_passed", "unknown_rejected")}, dict(
            known_common=2, known_common_passed=1, known_accidental=2, known_accidental_rejected=1,
            unknown_passed=1, unknown_rejected=1))
        self.assertEqual(result["genuine_efficiency_extrema"], [1 / 3, 2 / 3])
        self.assertEqual(result["accidental_rejection_extrema"], [1 / 3, 2 / 3])
        np.testing.assert_array_equal(y, before_y)
        np.testing.assert_array_equal(score, before_score)
        np.testing.assert_array_equal(passed, before_passed)
        self.assertIn("not a physical uncertainty", result["interpretation"])


if __name__ == "__main__":
    unittest.main()
