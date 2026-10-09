"""ROOT augmentation preserves events, arrays, weights and opaque metadata."""
import copy
from array import array
import json
from pathlib import Path
import sys
import tempfile
import unittest

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from add_bdt_score import DEFAULT_BRANCH, add_score, open_root, scores_for_arrays
from features import RECO_BRANCHES
from test_contract import reconstructed_event
from test_portable_inference import simple_model


def fixture_arrays():
    single = reconstructed_event()
    arrays = dict(run=np.array([1,1,1], dtype=np.uint32),
                  luminosityBlock=np.array([2,2,2], dtype=np.uint32),
                  event=np.array([2**53+1,2**53+2,2**53+3], dtype=np.uint64),
                  nShiftMuon=np.array([2,0,2]), nShiftDimuonVertex=np.array([2,0,1]))
    for name in RECO_BRANCHES:
        if name in arrays:
            continue
        base = single[name][0]
        arrays[name] = ak.Array([base.tolist() if name.startswith('ShiftMuon_') else [float(base[0]),float(base[0])],
                                [], base.tolist()])
    arrays['ShiftDimuonVertex_muonIdx1'] = ak.Array([[0,1],[],[0]])
    arrays['ShiftDimuonVertex_muonIdx2'] = ak.Array([[1,0],[],[1]])
    arrays['ShiftDimuonVertex_dca'] = ak.Array([[.2,2.],[],[.3]])
    return arrays


def write_fixture(path, arrays):
    # Native ROOT streamers match the CMSSW-created Nano schema. Uproot's
    # synthetic older TBranch descriptors are not the production dictionary.
    import ROOT
    root=ROOT.TFile(str(path),'RECREATE')
    tree=ROOT.TTree('Events','Independent native Nano fixture')
    scalars={}
    for name in ('run','luminosityBlock','event','nShiftMuon','nShiftDimuonVertex'):
        buf=array('Q' if name=='event' else 'I',[0]);scalars[name]=buf
        tree.Branch(name,buf,name+('/l' if name=='event' else '/i'))
    buffers={}
    for name in RECO_BRANCHES:
        if name in scalars:continue
        counter='nShiftMuon' if name.startswith('ShiftMuon_') else 'nShiftDimuonVertex'
        integer=name.endswith(('muonIdx1','muonIdx2'))
        buf=array('I' if integer else 'd',[0]*3);buffers[name]=buf
        tree.Branch(name,buf,name+'['+counter+']'+('/i' if integer else '/D'))
    gen=array('f',[0.]);sampling=array('d',[0.]);unrelated=array('d',[0.]*3);n_unrelated=array('I',[0])
    tree.Branch('genWeight',gen,'genWeight/F')
    tree.Branch('shiftSamplingGenWeight',sampling,'shiftSamplingGenWeight/D')
    tree.Branch('nUnrelated',n_unrelated,'nUnrelated/i')
    tree.Branch('unrelatedJagged',unrelated,'unrelatedJagged[nUnrelated]/D')
    for i in range(3):
        for name,buf in scalars.items():buf[0]=int(arrays[name][i])
        for name,buf in buffers.items():
            for j,value in enumerate(arrays[name][i]):buf[j]=int(value) if buf.typecode=='I' else float(value)
        gen[0]=[-1.,.5,7.][i];sampling[0]=[123.456,0.,-1.2345678912345][i]
        values=[[1.,2.,3.],[],[4.]][i];n_unrelated[0]=len(values)
        for j,value in enumerate(values):unrelated[j]=value
        tree.Fill()
    tree.Write()
    run=array('I',[1]);lumi=array('I',[2]);total=array('d',[999.25])
    runs=ROOT.TTree('Runs','Preserved run normalization');runs.Branch('run',run,'run/i')
    runs.Branch('genEventSumw',total,'genEventSumw/D');runs.Fill();runs.Write()
    lumis=ROOT.TTree('LuminosityBlocks','Preserved lumi metadata');lumis.Branch('run',run,'run/i')
    lumis.Branch('luminosityBlock',lumi,'luminosityBlock/i');lumis.Fill();lumis.Write()
    ROOT.TObjString('unchanged opaque metadata').Write('tag')
    root.Close()


def write_model(path):
    model = simple_model()
    model.update(score_interpretation='Synthetic uncalibrated discriminator', score_storage='Float64')
    path.write_text(json.dumps(model))
    return model


class PreservingScoreProduction(unittest.TestCase):
    def test_append_double_scores_with_zero_and_multiple_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory); source=base/'source.root'; output=base/'scored.root'; model_path=base/'model.json'
            arrays=fixture_arrays(); write_fixture(source,arrays); model=write_model(model_path)
            expected,_=scores_for_arrays(model,arrays)
            report=add_score(source,output,model_path,sample_kind='simulation',chunk_events=1)
            self.assertTrue(report['complete']); self.assertFalse(report['selection_applied'])
            self.assertEqual((report['events'],report['retained_pairs'],report['zero_pair_events']),(3,3,1))
            self.assertTrue(report['verification']['all_original_content_equal'])
            with open_root(output) as root:
                saved=root['Events'][DEFAULT_BRANCH].array(library='ak')
                self.assertEqual(root['Events'][DEFAULT_BRANCH].typename,'double[]')
                self.assertEqual(ak.num(saved).to_list(),[2,0,1])
                for a,b in zip(saved,expected): np.testing.assert_array_equal(np.asarray(a),b)
                np.testing.assert_array_equal(root['Events']['event'].array(library='np'),arrays['event'])
                np.testing.assert_array_equal(root['Events']['shiftSamplingGenWeight'].array(library='np'),[123.456,0.,-1.2345678912345])
                metadata=json.loads(str(root[DEFAULT_BRANCH+'Metadata']))
                self.assertNotIn('complete',metadata)
                self.assertFalse(metadata['physics_ready']); self.assertFalse(metadata['selection_applied'])

    def test_mass_nan_is_preserved_and_cannot_filter_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory); source=base/'source.root'; output=base/'scored.root'; model_path=base/'model.json'
            arrays=fixture_arrays(); arrays['ShiftDimuonVertex_mass']=ak.Array([[np.nan,-1.],[],[np.inf]])
            write_fixture(source,arrays);write_model(model_path)
            report=add_score(source,output,model_path,sample_kind='simulation')
            self.assertEqual(report['retained_pairs'],3)
            self.assertTrue(report['verification']['all_original_content_equal'])

    def test_invalid_pair_never_publishes_a_filtered_file(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory); source=base/'source.root'; output=base/'scored.root'; model_path=base/'model.json'
            arrays=fixture_arrays();arrays['ShiftDimuonVertex_muonIdx1']=ak.Array([[0,99],[],[0]])
            write_fixture(source,arrays);write_model(model_path)
            with self.assertRaisesRegex(ValueError,'pair reference'):
                add_score(source,output,model_path,sample_kind='simulation')
            self.assertFalse(output.exists())
            self.assertTrue(output.with_name(output.name+'.incomplete').exists())
            report=json.loads(output.with_suffix('.bdt.json').read_text())
            self.assertFalse(report['complete']);self.assertIn('error',report)

    def test_no_overwrite_and_simulation_scope_are_enforced(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory); source=base/'source.root'; output=base/'scored.root'; model_path=base/'model.json'
            write_fixture(source,fixture_arrays());write_model(model_path)
            output.write_bytes(b'protected existing output')
            with self.assertRaises(FileExistsError):add_score(source,output,model_path,sample_kind='simulation')
            self.assertEqual(output.read_bytes(),b'protected existing output')
            with self.assertRaises(ValueError):add_score(source,source,model_path,sample_kind='simulation')
            with self.assertRaises(ValueError):add_score(source,base/'data.root',model_path,sample_kind='collision_data')


if __name__ == '__main__':
    unittest.main()
