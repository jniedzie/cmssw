#!/usr/bin/env python3
"""Bounded BDT-capacity study; select only on the original validation partition."""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import time
import sys

import numpy as np
import sklearn
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline

from evaluation import event_split, fraction, weighted_quantile, binned_efficiencies
from export_bdt import export_model
from features import FEATURE_NAMES, nuisance_transform
from inference import predict
from phase2_train import phase_weights, validation_summary
from uniform_bdt import train_uniform_bdt


def reference(path):
    path=Path(path)
    return dict(path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def write(path,value):
    with Path(path).open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False)
        stream.write('\n')


def candidates(phase='capacity'):
    result=[]
    if phase=='expanded':
        return [dict(name='expanded_historical_d2_s3_n160',family='uniform',uniformity_mode='constant',
                     n_estimators=160,max_depth=2,min_samples_leaf=None,strength=3.,seed=71)]+[
            dict(name=f'expanded_d3_s{strength}_n400',family='uniform',uniformity_mode='constant',
                 n_estimators=400,max_depth=3,min_samples_leaf=20,strength=float(strength),seed=71)
            for strength in (0,1,3)]
    if phase=='adaptive':
        for strength in (6.,10.):
            result.append(dict(name=f'constant_d3_s{strength:g}_n400',family='uniform',
                uniformity_mode='constant',n_estimators=400,max_depth=3,min_samples_leaf=20,strength=strength,seed=71))
        for depth in (2,3):
            for strength in (.5,1.,2.):
                result.append(dict(name=f'adaptive_d{depth}_s{strength:g}_n400'.replace('.','p'),
                    family='uniform',uniformity_mode='adaptive',n_estimators=400,max_depth=depth,
                    min_samples_leaf=12 if depth==2 else 20,strength=strength,seed=71,
                    adaptive_beta_cap=6.,adaptive_error_floor=1.e-6))
        return result
    for depth,strengths in ((3,(0.,.5,1.,3.)),(4,(0.,1.))):
        for strength in strengths:
            result.append(dict(name=f'uniform_d{depth}_s{strength:g}_n400'.replace('.','p'),
                family='uniform',n_estimators=400,max_depth=depth,
                min_samples_leaf=20 if depth==3 else 30,strength=strength,seed=71))
    for trees in (300,500):
        result.append(dict(name=f'gradient_d3_n{trees}',family='gradient',n_estimators=trees,
            max_depth=3,min_samples_leaf=20 if trees==300 else 30,learning_rate=.03,seed=71))
    return result


def empirical_summary(scores,y,process,nuisance,threshold,training_nuisance,training_y):
    passed=scores>=threshold
    return dict(auc=float(roc_auc_score(y,scores)),threshold=float(threshold),
        genuine=fraction(y==1,passed),accidental=fraction(y==0,passed),
        rejection=1-fraction(y==0,passed)['efficiency'],
        flatness=binned_efficiencies(nuisance,y,passed,training_nuisance,training_y),
        per_process={str(p):{name:fraction((process==p)&(y==label),passed)
            for label,name in ((1,'genuine'),(0,'accidental'))} for p in np.unique(process)})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset',type=Path,required=True)
    parser.add_argument('--reference-results',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--development-cache',type=Path)
    parser.add_argument('--candidate-phase',choices=('capacity','adaptive','expanded'),default='capacity')
    parser.add_argument('--expanded-cache',type=Path,
                        help='Expand SM fitting population while retaining every original event hash partition.')
    parser.add_argument('--additional-reference',type=Path,
                        help='Previously validation-locked trusted research model pickle; never refit.')
    args=parser.parse_args()
    manifest=json.loads((args.dataset/'inputs.json').read_text())
    assert manifest['complete'] and manifest['sample_kind']=='simulation_only'
    for filename,key in (('reco.npz','reco_sha256'),('labels.npz','labels_sha256')):
        assert reference(args.dataset/filename)['sha256']==manifest[key]
    with np.load(args.dataset/'reco.npz',allow_pickle=False) as arrays:data=dict(arrays)
    with np.load(args.dataset/'labels.npz',allow_pickle=False) as arrays:y=arrays['y']
    assert tuple(data['feature_names'])==FEATURE_NAMES
    assert set(np.unique(data['process']))<=set(('jpsi','dy','qcd'))
    split=event_split(data['ids'],data['process'],71)
    with np.load(args.reference_results/'splits.npz',allow_pickle=False) as saved:
        for key,value in (('split',split),('ids',data['ids']),('process',data['process'])):
            np.testing.assert_array_equal(value,saved[key])
    tr=(split==0)&(y>=0);va=(split==1)&(y>=0)
    assert tr.sum()==1693 and va.sum()==573
    expansion=None
    if args.expanded_cache:
        assert args.candidate_phase=='expanded'
        expanded_manifest=json.loads((args.expanded_cache/'inputs.json').read_text())
        assert expanded_manifest['complete'] and expanded_manifest['sample_kind']=='simulation_only'
        for filename,key in (('reco.npz','reco_sha256'),('labels.npz','labels_sha256')):
            assert reference(args.expanded_cache/filename)['sha256']==expanded_manifest[key]
        with np.load(args.expanded_cache/'reco.npz',allow_pickle=False) as arrays:expanded=dict(arrays)
        with np.load(args.expanded_cache/'labels.npz',allow_pickle=False) as arrays:expanded_y=arrays['y']
        assert tuple(expanded['feature_names'])==FEATURE_NAMES
        assert expanded['ids'].dtype.kind=='u' and expanded['ids'].dtype.itemsize==8
        def keys(table):
            return [(str(p),*(int(i) for i in identity),int(pair))
                    for p,identity,pair in zip(table['process'],table['ids'],table['pair_index'])]
        original_keys,expanded_keys=keys(data),keys(expanded)
        index={key:i for i,key in enumerate(expanded_keys)}
        assert len(index)==len(expanded_keys) and len(set(original_keys))==len(original_keys)
        matched=np.asarray([index[key] for key in original_keys],dtype=int)
        expanded_split=event_split(expanded['ids'],expanded['process'],71)
        np.testing.assert_array_equal(expanded_split[matched],split)
        np.testing.assert_array_equal(expanded['X'][matched],data['X'])
        np.testing.assert_array_equal(expanded['muon_indices'][matched],data['muon_indices'])
        np.testing.assert_array_equal(expanded_y[matched],y)
        sm=np.isin(expanded['process'],('jpsi','dy','qcd'))
        data={key:expanded[key][sm] for key in ('X','ids','nuisance','process','stratum','pair_index','source_index')}
        data['feature_names']=expanded['feature_names']
        y=expanded_y[sm];split=expanded_split[sm]
        tr=(split==0)&(y>=0);va=(split==1)&(y>=0)
        events=[{(str(p),*(int(i) for i in identity)) for p,identity in zip(data['process'][split==part],data['ids'][split==part])} for part in range(3)]
        assert all(not(events[a]&events[b]) for a,b in ((0,1),(0,2),(1,2)))
        expansion=dict(manifest=reference(args.expanded_cache/'inputs.json'),
            reco=reference(args.expanded_cache/'reco.npz'),labels=reference(args.expanded_cache/'labels.npz'),
            original_pairs_preserved=len(matched),original_split_feature_label_muon_identity_parity=True,
            physical_event_partition_disjoint=True,integer_identity_dtype=str(data['ids'].dtype),
            candidate_identity='process, uint64(run,lumi,event), pair_index; all unique',
            sm_pairs=int(sm.sum()),excluded_signal_pairs=int((~sm).sum()),
            partitions={name:dict(pairs=int(np.sum(split==part)),events=len(events[part]),
                labels={str(label):int(np.sum((split==part)&(y==label))) for label in (-1,0,1)})
                for part,name in enumerate(('train','validation','development_test'))},
            interpretation='All cache rows previously inspected; expanded hash2 test is event-disjoint development evidence, never fresh final evaluation.')
    elif args.candidate_phase=='expanded':
        raise ValueError('Expanded phase requires --expanded-cache')
    X=np.where(np.isfinite(data['X']),data['X'],np.nan)
    Z=nuisance_transform(data['nuisance'])
    wtr,wva=(phase_weights(y[mask],data['process'][mask]) for mask in (tr,va))
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'models').mkdir()
    (args.output/'source_snapshot').mkdir()
    sources={}
    for path in sorted(Path(__file__).parent.glob('*.py')):
        sources[path.name]=reference(path)
        (args.output/'source_snapshot'/path.name).write_bytes(path.read_bytes())
    references={name:reference(args.reference_results/'models'/(name+'.pkl'))
                for name in ('uniform_bdt_3p0_s71','bdt_quality_s71')}
    if args.additional_reference:
        previous_lock=args.additional_reference.parent.parent/'selection_locked.json'
        previous=json.loads(previous_lock.read_text())
        assert previous['selection_uses_test'] is False and args.additional_reference.stem in previous['representative_models']
        references[args.additional_reference.stem]=reference(args.additional_reference)
    protocol=dict(schema='shift-bdt-capacity-protocol-v1',candidate_phase=args.candidate_phase,
        candidates=candidates(args.candidate_phase),fit_budget=len(candidates(args.candidate_phase)),
        dataset=reference(args.dataset/'inputs.json'),reco=reference(args.dataset/'reco.npz'),
        labels=reference(args.dataset/'labels.npz'),splits=reference(args.reference_results/'splits.npz'),
        references=references,source=sources,
        expansion=expansion,
        train_pairs=int(tr.sum()),validation_pairs=int(va.sum()),train_counts={str(k):int(np.sum(tr&(y==k))) for k in (0,1)},
        validation_counts={str(k):int(np.sum(va&(y==k))) for k in (0,1)},
        predictors=list(FEATURE_NAMES),predictor_contract='Existing reconstructed features only; no explicit mass, truth, label, process or weight columns. Momentum and energy inputs can encode mass.',
        training_weights='Unchanged phase_weights: equal classes, process smoothing count50; weights never enter predictors.',
        threshold='Primary comparison: empirical validation-positive10th quantile, inclusive >=, target90% genuine retention; ties reported. Original process-balanced90% target also reported.',
        selection='Highest primary validation accidental rejection; then validation AUC, then name. Select overall and portable-only winners independently. Original validation flatness selector supplies a separate constrained candidate.',
        acceptance_envelope='Additional descriptive gate at original balanced90% working point: worst supported marginal span, process-efficiency span and maximum class joint dCor must each be no worse than frozen uniform3 reference; no independence claim follows.',
        adaptive_method='Optional single-point variant inspired by https://arxiv.org/pdf/1305.7248 equations2.3--2.4. Post-classification weights normalized over all rows; error sums absolute global-local deviations over valid positives. beta=clip(log((1-e)/e),0,6), e clipped[1e-6,1-1e-6]; log update=clip(strength*beta*(global-local),-1,1), without extra classification learning-rate factor. Historical constant path unchanged; not full paper uBoost.',
        test_access='Original test and current full firstlook already inspected; no test scores or metrics accessed before validation lock. All later test use is consumed development diagnostics, never fresh final evaluation.',
        deploy=False,physics_selection_ready=False,signal_training=False,
        runtime=dict(python=sys.version,numpy=np.__version__,sklearn=sklearn.__version__))
    write(args.output/'protocol.json',protocol)
    np.savez_compressed(args.output/'splits.npz',split=split,ids=data['ids'],process=data['process'],pair_index=data['pair_index'])
    bundles,results={},{}
    started=time.monotonic()

    def archive(name,bundle,description,fit_seconds):
        sv=predict(bundle,X[va]);st=predict(bundle,X[tr])
        balanced=validation_summary(sv,y[va],data['process'][va],data['nuisance'][va],wva,data['nuisance'][tr],y[tr])
        cut=weighted_quantile(sv[y[va]==1],np.ones(np.sum(y[va]==1)),.1)
        primary=empirical_summary(sv,y[va],data['process'][va],data['nuisance'][va],cut,data['nuisance'][tr],y[tr])
        train=empirical_summary(st,y[tr],data['process'][tr],data['nuisance'][tr],cut,data['nuisance'][tr],y[tr])
        legacy=bundle.get('threshold')
        bundle=dict(bundle,threshold=cut)
        bundles[name]=bundle
        results[name]=dict(description=description,fit_seconds=fit_seconds,validation=primary,
            process_balanced_validation=balanced,train_at_validation_cut=train,
            training_minus_validation_auc=train['auc']-primary['auc'],legacy_threshold=legacy,
            portable=isinstance(bundle['model'],__import__('uniform_bdt').UniformityBDT))
        with (args.output/'models'/(name+'.pkl')).open('xb') as stream:pickle.dump(bundle,stream)
        write(args.output/(name+'_validation.json'),results[name])
        print(json.dumps(dict(model=name,auc=primary['auc'],rejection=primary['rejection'],
            genuine=primary['genuine']['efficiency'],train_auc=train['auc'],seconds=fit_seconds)),flush=True)

    for name,ref in references.items():
        with Path(ref['path']).open('rb') as stream:bundle=pickle.load(stream)
        archive('reference_'+name,bundle,dict(family='frozen_reference',source=ref),0.)
    for spec in protocol['candidates']:
        begin=time.monotonic()
        if spec['family']=='uniform':
            model,history=train_uniform_bdt(X[tr],y[tr],wtr,Z[tr],seed=spec['seed'],
                strength=spec['strength'],n_estimators=spec['n_estimators'],
                max_depth=spec['max_depth'],min_samples_leaf=spec['min_samples_leaf'],
                uniformity_mode=spec.get('uniformity_mode','constant'),
                adaptive_beta_cap=spec.get('adaptive_beta_cap',6.),
                adaptive_error_floor=spec.get('adaptive_error_floor',1.e-6))
            write(args.output/(spec['name']+'_history.json'),history)
        else:
            model=make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),
                GradientBoostingClassifier(n_estimators=spec['n_estimators'],max_depth=spec['max_depth'],
                    min_samples_leaf=spec['min_samples_leaf'],learning_rate=spec['learning_rate'],
                    subsample=.8,random_state=spec['seed']))
            model.fit(X[tr],y[tr],gradientboostingclassifier__sample_weight=wtr)
        archive(spec['name'],dict(model=model,columns=np.arange(len(FEATURE_NAMES)),
            feature_names=FEATURE_NAMES,requires_gen=False,deployment_ready=False),spec,time.monotonic()-begin)
    write(args.output/'validation_candidates.json',results)
    def rank(name):return (results[name]['validation']['rejection'],results[name]['validation']['auc'],name)
    overall=max(results,key=rank)
    portable=max((name for name in results if results[name]['portable']),key=rank)
    new_candidate=max((spec['name'] for spec in protocol['candidates']),key=rank)
    constrained=max(results,key=lambda name:(results[name]['process_balanced_validation']['selector_objective'],name))
    def envelope(entry):
        values=entry['process_balanced_validation']
        return (values['maximum_marginal_efficiency_span'],values['process_efficiency_span'],
                max(values['dependence'][label]['joint_dcor'] or 0. for label in ('common_vertex','accidental')))
    reference_envelope=envelope(results['reference_uniform_bdt_3p0_s71'])
    eligible=[name for name in results if all(x<=bound+1.e-12 for x,bound in zip(envelope(results[name]),reference_envelope))]
    envelope_winner=max(eligible,key=rank)
    representatives=list(dict.fromkeys((overall,portable,new_candidate,constrained,envelope_winner,
        'reference_uniform_bdt_3p0_s71','reference_bdt_quality_s71')))
    lock=dict(schema='shift-bdt-capacity-selection-v1',selection_uses_test=False,
        protocol=reference(args.output/'protocol.json'),validation=reference(args.output/'validation_candidates.json'),
        overall_winner=overall,portable_winner=portable,original_flatness_selector_winner=constrained,
        new_candidate_winner=new_candidate,
        acceptance_envelope_winner=envelope_winner,acceptance_envelope_eligible=eligible,
        acceptance_envelope_reference_values=list(reference_envelope),
        representative_models=representatives,thresholds={name:bundle['threshold'] for name,bundle in bundles.items()},
        original_balanced_thresholds={name:entry['process_balanced_validation']['threshold'] for name,entry in results.items()},
        interpretation='Validation-selected research candidates only; no production model replacement or physics cut.')
    write(args.output/'selection_locked.json',lock)
    # Selection and all thresholds are immutable before any consumed test predictions.
    score_table=dict(ids=data['ids'],process=data['process'],split=split,y=y)
    if 'pair_index' in data:score_table['pair_index']=data['pair_index']
    if 'source_index' in data:score_table['source_index']=data['source_index']
    score_table.update({name:predict(bundle,X) for name,bundle in bundles.items()})
    np.savez_compressed(args.output/'candidate_scores.npz',**score_table)
    for name in representatives:
        if results[name]['portable']:
            export_model(args.output,name,args.output/'portable'/(name+'.json'),args.dataset/'reco.npz')
    development=None
    if args.development_cache:
        cachemanifest=json.loads((args.development_cache/'inputs.json').read_text())
        assert reference(args.development_cache/'reco.npz')['sha256']==cachemanifest['reco_sha256']
        with np.load(args.development_cache/'reco.npz',allow_pickle=False) as arrays:cached=dict(arrays)
        assert tuple(cached['feature_names'])==FEATURE_NAMES
        table={key:cached[key] for key in ('ids','process','pair_index','source_index','heldout')}
        if args.expanded_cache:
            table['expanded_split']=event_split(cached['ids'],cached['process'],71)
            table['expanded_heldout']=(table['expanded_split']==2)&np.isin(cached['process'],('jpsi','dy','qcd'))
        for name in representatives:table[name]=predict(bundles[name],cached['X'])
        np.savez_compressed(args.output/'development_scores.npz',**table)
        development=dict(cache=reference(args.development_cache/'inputs.json'),scores=reference(args.output/'development_scores.npz'),
                         models=representatives,rows=len(cached['X']),interpretation='Consumed development diagnostic only, includes signal inference but never signal training or selection.')
        if args.expanded_cache:
            development.update(required_sm_test_mask='expanded_heldout; never use historical heldout for expanded-model testing',
                expanded_heldout_pairs=int(table['expanded_heldout'].sum()))
    write(args.output/'completion.json',dict(complete=True,elapsed_seconds=time.monotonic()-started,
        protocol=reference(args.output/'protocol.json'),selection=reference(args.output/'selection_locked.json'),
        validation=reference(args.output/'validation_candidates.json'),scores=reference(args.output/'candidate_scores.npz'),
        development=development,physics_ready=False))
    print(json.dumps(dict(complete=True,overall=overall,portable=portable,constrained=constrained)),flush=True)


if __name__=='__main__':main()
