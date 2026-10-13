# Code style: AI-code. AI-generated synthetic staged-provenance regression tests.
"""Exercise guards without loading physics arrays, executable pickles, or fitting."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import neural_vertex_study as study


class StudyGuardTests(unittest.TestCase):
    def fixture(self, directory):
        output = Path(directory) / 'study'
        output.mkdir()
        previous = Path(directory) / 'previous'
        previous.mkdir()
        refs = {}
        for key in ('source', 'cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits', 'model', 'history'):
            path = (previous if key == 'bdt_scores' else output) / (key + '.npz')
            path.write_bytes(b'synthetic original bytes; no physics or pickle')
            refs[key] = study.reference(path)
        protocol = {key: refs[key] for key in ('cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits')}
        protocol.update(sources={'synthetic.py': refs['source']}, recipes=[{'name': 'synthetic'}], acceptance_envelope={})
        study.write(output / 'protocol.json', protocol)
        result = {'protocol': study.reference(output / 'protocol.json'), 'model': refs['model'], 'history': refs['history'],
                  'validation': {'rejection': .5, 'auc': .7, 'selector_objective': .1, 'threshold': .5}}
        study.write(output / 'synthetic_validation.json', result)
        study.write(previous / 'selection_locked.json', {'thresholds': {key: .5 for key in self.reference_names}})
        return output, refs

    reference_names = ('reference_uniform_bdt_3p0_s71', 'expanded_d3_s0_n400', 'expanded_d3_s3_n400')

    def assert_rejected_before_load(self, output, action, message):
        with patch.object(study, 'load') as arrays, patch.object(study.pickle, 'load') as unpickle:
            with self.assertRaisesRegex(ValueError, message):
                action()
            arrays.assert_not_called()
            unpickle.assert_not_called()
        self.assertFalse((output / 'selection_locked.json').exists())
        self.assertFalse((output / 'completion.json').exists())

    def test_all_input_and_source_mutations_stop_fit_before_data_read(self):
        for key in ('source', 'cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits'):
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                output, refs = self.fixture(directory)
                Path(refs[key]['path']).write_bytes(b'changed after protocol freeze')
                with patch('neural_vertex.train_neural_vertex') as train:
                    self.assert_rejected_before_load(output, lambda: study.fit_one(output, {'name': 'synthetic'}),
                        'Training source changed' if key == 'source' else 'Fitting input changed: ' + key)
                    train.assert_not_called()

    def test_all_input_and_source_mutations_stop_finalize_before_data_read(self):
        for key in ('source', 'cache_manifest', 'reco', 'labels', 'bdt_scores', 'splits'):
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                output, refs = self.fixture(directory)
                Path(refs[key]['path']).write_bytes(b'changed after candidate archive')
                self.assert_rejected_before_load(output, lambda: study.finalize(output),
                    'Source changed' if key == 'source' else 'Evaluation input changed: ' + key)

    def test_candidate_artifact_mutations_stop_before_unpickling(self):
        for key in ('model', 'history'):
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                output, refs = self.fixture(directory)
                # A threshold change serialized into the model changes these bytes too.
                Path(refs[key]['path']).write_bytes(b'changed serialized threshold or history')
                self.assert_rejected_before_load(output, lambda: study.finalize(output),
                    'Candidate artifact changed: synthetic/' + key)

    def test_candidate_protocol_mutation_stops_before_unpickling(self):
        with tempfile.TemporaryDirectory() as directory:
            output, refs = self.fixture(directory)
            protocol = json.loads((output / 'protocol.json').read_text())
            protocol['threshold'] = 'changed recipe after candidate archive'
            (output / 'protocol.json').write_text(json.dumps(protocol))
            self.assert_rejected_before_load(output, lambda: study.finalize(output), 'Candidate protocol differs')

    def test_hash_valid_bundle_threshold_disagreement_stops_before_prediction(self):
        with tempfile.TemporaryDirectory() as directory:
            output, refs = self.fixture(directory)
            # Synthetic arrays only: no input reader is used in this semantic check.
            split = {'expanded_split': np.array([0, 2]), 'expanded_heldout': np.array([False, True])}
            reco = {'ids': np.zeros((2, 3), dtype=np.uint64), 'process': np.array(['dy', 'dy']),
                    'pair_index': np.zeros(2, dtype=int), 'source_index': np.zeros(2, dtype=int),
                    'heldout': split['expanded_heldout'], 'X': np.zeros((2, len(study.FEATURE_NAMES))),
                    'nuisance': np.zeros((2, 3)), 'topology': np.array(['category', 'category']), **split}
            earlier = {**reco, **{key: np.zeros(2) for key in self.reference_names}}
            datasets = {refs['reco']['path']: reco, refs['labels']['path']: {'y': np.array([1, 0])},
                        refs['bdt_scores']['path']: earlier, str(output / 'splits.npz'): split}
            model = Mock()
            bundle = {'feature_names': study.FEATURE_NAMES, 'columns': np.arange(len(study.FEATURE_NAMES)),
                      'threshold': .4, 'model': model}
            with patch.object(study, 'load', side_effect=lambda path: datasets[str(path)]), \
                 patch.object(study, 'envelope_eligible', return_value=False), \
                 patch.object(study.pickle, 'load', return_value=bundle) as unpickle:
                with self.assertRaisesRegex(ValueError, 'Candidate threshold differs from validation lock'):
                    study.finalize(output)
                unpickle.assert_called_once()
            model.predict_proba.assert_not_called()
            self.assertFalse((output / 'development_scores.npz').exists())
            self.assertFalse((output / 'completion.json').exists())


if __name__ == '__main__':
    unittest.main()
