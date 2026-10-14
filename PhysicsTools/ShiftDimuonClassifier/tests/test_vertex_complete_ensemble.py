# Code style: AI-code. AI-generated ensemble contract and round-trip checks.
from pathlib import Path
import pickle
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from vertex_complete_ensemble import CompleteVertexEnsemble


class ConstantModel:
    n_features_in_=119
    def __init__(self, probability): self.probability=probability
    def predict_proba(self,x,*category):
        return np.tile([1-self.probability,self.probability],(len(x),1))


class EnsembleTests(unittest.TestCase):
    def test_equal_weight_mean_and_trusted_reload(self):
        model=CompleteVertexEnsemble([ConstantModel(.2),ConstantModel(.6)],[False,True])
        x=np.zeros((2,119));c=np.asarray(['Near-Both','Both-Both'])
        expected=np.asarray([[.6,.4],[.6,.4]])
        np.testing.assert_allclose(model.predict_proba(x,c),expected)
        np.testing.assert_array_equal(pickle.loads(pickle.dumps(model)).predict_proba(x,c),model.predict_proba(x,c))

    def test_incomplete_members_unknown_category_and_bad_probabilities_fail(self):
        member=ConstantModel(.5);member.n_features_in_=50
        with self.assertRaises(ValueError):CompleteVertexEnsemble([member],[False])
        model=CompleteVertexEnsemble([ConstantModel(.5)],[False])
        with self.assertRaises(ValueError):model.predict_proba(np.zeros((1,119)),['future-category'])
        with self.assertRaises(ValueError):CompleteVertexEnsemble([ConstantModel(np.nan)],[False]).predict_proba(np.zeros((1,119)),['Other'])
        self.assertEqual(model.predict_proba(np.empty((0,119)),np.asarray([],dtype=str)).shape,(0,2))


if __name__=='__main__':unittest.main()
