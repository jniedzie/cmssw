"""Portable numeric/schema contracts and exact trusted-pickle export parity."""
import copy
import json
import math
from pathlib import Path
import pickle
import sys
import tempfile
import unittest

import numpy as np

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE))
from features import FEATURE_NAMES
from portable_inference import (INPUT_CAST, SCHEMA, SCORE_FORMULA, feature_contract_digest,
                                file_digest, load_model, predict, validate_model)
from export_bdt import export_bundle, export_model, validate_scores


def provenance():
    return dict(pickle_sha256="1" * 64, selection_lock_sha256="2" * 64,
                source_sha256={name: file_digest(PACKAGE / name) for name in
                               ("features.py", "inference.py", "uniform_bdt.py", "portable_inference.py", "export_bdt.py")})


def simple_model():
    return dict(schema=SCHEMA, feature_names=list(FEATURE_NAMES),
                feature_contract_sha256=feature_contract_digest(FEATURE_NAMES), columns=[0],
                imputer_statistics=[0.], imputer_nonfinite_policy="replace-with-frozen-training-median",
                input_cast=INPUT_CAST, score_formula=SCORE_FORMULA, classes=[0, 1],
                sum_abs_coefficients=1., threshold_metadata=.5, requires_gen=False,
                selection_applied=False, physics_ready=False, provenance=provenance(),
                trees=[dict(coefficient=1., children_left=[1,-1,-1], children_right=[2,-1,-1],
                            feature=[0,-2,-2], threshold=[1.00000008,-2.,-2.], leaf_class=[None,0,1])])


def matrix(values):
    result = np.zeros((len(values), len(FEATURE_NAMES)), dtype=np.float64)
    result[:, 0] = values
    return result


class PortableNumericContract(unittest.TestCase):
    def test_float32_inputs_compare_against_float64_thresholds(self):
        # The Float64 threshold rounds UP to the same Float32 value as the
        # right-hand row. A Float32-threshold comparison wrongly sends it left.
        score = predict(simple_model(), matrix([1.00000005, 1.00000009, 1.00000008]))
        np.testing.assert_array_equal(score, [1/(1+math.exp(2)), 1/(1+math.exp(-2)), 1/(1+math.exp(-2))])
        self.assertEqual(score.dtype, np.float64)

    def test_nonfinite_values_use_frozen_imputation_before_float32_cast(self):
        model = simple_model()
        model["imputer_statistics"] = [1.00000009]
        result = predict(model, matrix([np.nan, np.inf, -np.inf, 1.00000009]))
        np.testing.assert_array_equal(result, np.repeat(1/(1+math.exp(-2)), 4))
        self.assertEqual(model["imputer_statistics"], [1.00000009])

    def test_float32_overflow_and_wrong_feature_shape_fail(self):
        with self.assertRaisesRegex(ValueError, "overflow"):
            predict(simple_model(), matrix([np.finfo(np.float64).max]))
        with self.assertRaisesRegex(ValueError, "matrix shape"):
            predict(simple_model(), np.zeros((1, len(FEATURE_NAMES)-1)))

    def test_zero_rows_and_neutral_ensemble_are_defined_without_selection(self):
        model = simple_model()
        self.assertEqual(predict(model, matrix([])).shape, (0,))
        model["trees"] = []; model["sum_abs_coefficients"] = 0.
        np.testing.assert_array_equal(predict(model, matrix([0.,1.])), [.5,.5])
        self.assertFalse(model["selection_applied"])

    def test_json_roundtrip_retains_scores_and_double_precision(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/"model.json"
            path.write_text(json.dumps(simple_model()))
            restored = load_model(path)
            np.testing.assert_array_equal(predict(restored, matrix([0.,2.])), predict(simple_model(), matrix([0.,2.])))


class PortableSchemaContract(unittest.TestCase):
    def test_malformed_features_stats_coefficients_classes_and_provenance_fail(self):
        mutations = (lambda m:m.update(schema="unknown"), lambda m:m["feature_names"].reverse(),
                     lambda m:m.update(feature_contract_sha256="0"*64), lambda m:m.update(columns=[-1]),
                     lambda m:m.update(columns=[0,0]), lambda m:m.update(columns=[True]),
                     lambda m:m.update(imputer_statistics=[float("nan")]),
                     lambda m:m.update(imputer_statistics=[float("inf")]),
                     lambda m:m.update(classes=[1,0]), lambda m:m.update(classes=[False,True]),
                     lambda m:m.update(provenance=None), lambda m:m["provenance"].update(source_sha256=None),
                     lambda m:m.update(requires_gen=True),
                     lambda m:m.update(selection_applied=True), lambda m:m.update(physics_ready=True),
                     lambda m:m.update(sum_abs_coefficients=2.), lambda m:m.update(threshold_metadata=float("nan")),
                     lambda m:m["trees"][0].update(coefficient=-1.), lambda m:m["trees"][0].update(coefficient=float("inf")),
                     lambda m:m["provenance"].update(pickle_sha256="not-a-digest"),
                     lambda m:m["provenance"]["source_sha256"].update({"features.py":"0"*64}))
        for index, mutate in enumerate(mutations):
            model = simple_model(); mutate(model)
            with self.subTest(index=index), self.assertRaises(ValueError): validate_model(model)

    def test_cycles_shared_children_invalid_references_and_unreachable_nodes_fail(self):
        mutations = (lambda t:t.update(children_left=[0,-1,-1]),
                     lambda t:t.update(children_left=[3,-1,-1]),
                     lambda t:t.update(children_right=[1,-1,-1]),
                     lambda t:t.update(children_left=[-1,-1,-1]),
                     lambda t:t.update(feature=[1,-2,-2]),
                     lambda t:t.update(leaf_class=[None,2,1]),
                     lambda t:t.update(leaf_class=[1,0,1]),
                     lambda t:t.update(threshold=[float("nan"),-2.,-2.]))
        for index, mutate in enumerate(mutations):
            model = simple_model(); mutate(model["trees"][0])
            with self.subTest(index=index), self.assertRaises(ValueError): validate_model(model)
        model = simple_model(); tree = model["trees"][0]
        for key, value in (("children_left",-1),("children_right",-1),("feature",-2),("threshold",-2.),("leaf_class",0)):
            tree[key].append(value)
        with self.assertRaisesRegex(ValueError,"unreachable"):
            validate_model(model)


class TrustedExportContract(unittest.TestCase):
    def bundle(self):
        from uniform_bdt import train_uniform_bdt
        rng = np.random.default_rng(163)
        X = rng.normal(size=(100,len(FEATURE_NAMES)))
        y = (X[:,0] + .3*X[:,2] > 0).astype(np.int8)
        X[::9,0] = np.nan
        model,_ = train_uniform_bdt(X[:,[0,2]],y,np.ones(len(y)),X[:,3:6],strength=0,n_estimators=12)
        bundle = dict(model=model,columns=np.array([0,2]),feature_names=FEATURE_NAMES,threshold=.4403,
                      requires_gen=False,deployment_ready=False)
        return bundle,X

    def test_trusted_export_has_bitwise_parity_and_preserves_leaf_ties(self):
        bundle,X = self.bundle()
        exported = export_bundle(bundle,provenance())
        self.assertEqual(validate_scores(exported,bundle,X)["rows"],len(X))
        self.assertFalse(exported["requires_gen"]); self.assertFalse(exported["selection_applied"])
        self.assertEqual(exported["columns"],[0,2])
        # sklearn's argmax tie goes to its first class. Check independently on
        # every stored node after forcing a two-class tie in a synthetic tree.
        tree=bundle["model"].estimators_[0].tree_
        leaves=tree.children_left==-1
        tree.value[leaves,0,:]=.5
        tied=export_bundle(bundle,provenance())
        self.assertTrue(all(c==0 for c in tied["trees"][0]["leaf_class"] if c is not None))

    def test_export_freezes_model_and_receipt_and_refuses_overwrite_or_unlocked_cut(self):
        bundle,X=self.bundle()
        with tempfile.TemporaryDirectory() as d:
            base=Path(d); results=base/"results"; (results/"models").mkdir(parents=True)
            name="fixture_uniform"
            (results/"models"/(name+".pkl")).write_bytes(pickle.dumps(bundle))
            lock=dict(representative_models=[name],thresholds={name:bundle["threshold"]},selection_uses_test=False)
            (results/"selection_locked.json").write_text(json.dumps(lock))
            table=base/"reco.npz"; np.savez(table,X=X,feature_names=np.asarray(FEATURE_NAMES))
            path=base/"export.json"
            model,receipt=export_model(results,name,path,table)
            self.assertTrue(receipt["production_scoring_parity_validated"])
            self.assertEqual(receipt["model_sha256"],file_digest(path))
            self.assertEqual(load_model(path),model)
            with self.assertRaises(FileExistsError): export_model(results,name,path,table)
            lock["thresholds"][name]=.7; (results/"selection_locked.json").write_text(json.dumps(lock))
            with self.assertRaisesRegex(ValueError,"threshold differs"):
                export_model(results,name,base/"other.json",table)


if __name__ == "__main__":
    unittest.main()
