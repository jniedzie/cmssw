#!/usr/bin/env python3
"""Predeclared ablations, uniform BDT and DisCo on a fresh SM holdout."""
import argparse
from collections import Counter
import hashlib
import json
import pickle
from pathlib import Path
import time
import numpy as np
import sklearn
import torch
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from adversarial import train_adversarial
from disco import train_disco
from uniform_bdt import train_uniform_bdt
from evaluation import binned_efficiencies, dependence, evaluate, event_split, fraction, weighted_quantile
from features import FEATURE_NAMES, QUALITY_NAMES, nuisance_transform
from inference import predict


def phase_weights(y, process, minimum=50):
    """Balance labels and smooth process balancing to avoid huge rare-class weights."""
    weights = np.zeros(len(y))
    for label in (0, 1):
        samples = np.unique(process[y == label])
        masses = np.array([min(1., np.count_nonzero((y == label) & (process == p)) / minimum)
                           for p in samples])
        for p, mass in zip(samples, masses):
            mask = (y == label) & (process == p)
            weights[mask] = 0.5 * mass / masses.sum() / mask.sum()
    return weights * len(y)


def feature_groups():
    quality = np.arange(len(QUALITY_NAMES))
    geometry = [i for i,n in enumerate(FEATURE_NAMES) if n in (
        "pair_vx", "pair_vy", "pair_vz", "pair_originCompatibilityNormalizedChi2", "mu_vertex_distance")
        or n.startswith(("mu_highPt_v", "mu_lowPt_v"))]
    kinematics = [i for i,n in enumerate(FEATURE_NAMES) if n in ("pair_pt", "pair_pz", "pair_eta", "pair_phi")
                  or n.startswith(("mu_highPt_p", "mu_lowPt_p", "mu_highPt_energy", "mu_lowPt_energy"))]
    return dict(bdt_quality=quality, bdt_geometry=np.unique(np.r_[quality, geometry]),
                bdt_kinematics=np.unique(np.r_[quality, kinematics]),
                bdt_angle=np.r_[quality, FEATURE_NAMES.index("opening_angle")],
                bdt_expanded=np.arange(len(FEATURE_NAMES)))


def validation_summary(score, y, process, nuisance, weights, train_nuisance, train_y):
    cut = weighted_quantile(score[y == 1], weights[y == 1], 0.1)
    passed = score >= cut
    flatness = binned_efficiencies(nuisance, y, passed, train_nuisance, train_y)
    spans = []
    for variable in ("mass", "vertex_z", "vertex_radius"):
        values = [item["efficiency"] for item in flatness[variable]["bins"] if item["supported"]]
        spans.append(max(values) - min(values) if len(values) >= 2 else 1.)
    deps = dependence(score, nuisance, y)
    per_process = {str(p): {name: fraction((process == p) & (y == label), passed)
                           for label,name in ((0,"accidental"),(1,"common_vertex"))}
                   for p in np.unique(process)}
    process_flatness = {}
    process_efficiencies = []
    for p in np.unique(process):
        pmask = process == p
        if np.count_nonzero(pmask & (y == 1)) < 50:
            continue
        process_efficiencies.append(per_process[str(p)]["common_vertex"]["efficiency"])
        pf = binned_efficiencies(nuisance[pmask],y[pmask],passed[pmask],train_nuisance,train_y)
        process_flatness[str(p)] = pf
        for variable in ("mass","vertex_z","vertex_radius"):
            values = [item["efficiency"] for item in pf[variable]["bins"] if item["supported"]]
            if len(values) >= 2:
                spans.append(max(values)-min(values))
    process_span = max(process_efficiencies)-min(process_efficiencies) if len(process_efficiencies)>=2 else 0.
    # Predeclared selector: preserve acceptance while rewarding useful rejection.
    reject = 1 - fraction(y == 0, passed)["efficiency"]
    joint = max(deps[k]["joint_dcor"] or 0. for k in ("common_vertex", "accidental"))
    objective = reject - 2 * max(spans) - 2 * process_span - 0.25 * joint
    return dict(auc=float(roc_auc_score(y, score)), threshold=cut,
                common_vertex=fraction(y == 1, passed), accidental=fraction(y == 0, passed),
                balanced_common_efficiency=float(np.average(passed[y == 1], weights=weights[y == 1])),
                dependence=deps, flatness=flatness, maximum_marginal_efficiency_span=max(spans),
                selector_objective=objective,
                process_efficiency_span=process_span, per_process=per_process,
                per_process_flatness=process_flatness)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seeds", type=int, nargs="+", default=[71, 72])
    parser.add_argument("--split-seed", type=int, default=71)
    parser.add_argument("--epochs", type=int, default=120)
    parser.add_argument("--disco-strengths", type=float, nargs="+", default=[0., 0.3, 1., 3.])
    parser.add_argument("--uniform-strengths", type=float, nargs="+", default=[0., 1., 3.])
    args = parser.parse_args()
    manifest = json.loads((args.dataset / "inputs.json").read_text())
    if (manifest.get("sample_kind") != "simulation_only" or manifest.get("complete") is not True
            or manifest.get("no_pilot_event_overlap") is not True or manifest.get("no_pilot_file_overlap") is not True):
        raise ValueError("Phase 2 requires a completed, independently checked fresh MC dataset")
    for filename,key in (("reco.npz","reco_sha256"),("labels.npz","labels_sha256")):
        if hashlib.sha256((args.dataset / filename).read_bytes()).hexdigest() != manifest[key]:
            raise ValueError("Dataset digest differs")
    args.output.mkdir(parents=True, exist_ok=False)
    models_dir = args.output / "models"
    models_dir.mkdir()
    data = dict(np.load(args.dataset / "reco.npz", allow_pickle=False))
    lab = dict(np.load(args.dataset / "labels.npz", allow_pickle=False))
    if tuple(data["feature_names"]) != FEATURE_NAMES:
        raise ValueError("Feature order differs")
    y = lab["y"]
    split = event_split(data["ids"], data["process"], args.split_seed)
    tr, va, te = [(split == i) & (y >= 0) for i in range(3)]
    for mask in (tr,va,te):
        if np.count_nonzero(mask & (y == 0)) < 20 or np.count_nonzero(mask & (y == 1)) < 50:
            raise ValueError("Insufficient labelled event partition")
    X = np.where(np.isfinite(data["X"]), data["X"], np.nan)
    Z = nuisance_transform(data["nuisance"])
    wtr, wva = (phase_weights(y[mask], data["process"][mask]) for mask in (tr,va))
    np.savez_compressed(args.output / "splits.npz", split=split, ids=data["ids"], process=data["process"])
    protocol = dict(phase="fresh-SM-2", split_seed=args.split_seed, fit_seeds=args.seeds, epochs=args.epochs,
                    target_validation_common_efficiency=0.9,
                    disco_strengths=args.disco_strengths, uniform_strengths=args.uniform_strengths,
                    adversarial_reference_strength=0.5,
                    rare_process_smoothing_count=50,
                    selection="validation rejection - 2*worst supported marginal efficiency span (combined and per process) - 2*supported process efficiency spread - 0.25*maximum class joint dCor",
                    selection_support="20 genuine pairs/bin; 50 genuine pairs/process; full joint cut-efficiency closure is checked separately and is not guaranteed by ranking",
                    test_access="Partition class counts inspected for adequacy; no test scores or performance used until all models/thresholds are locked",
                    deployment_ready=False, lifetime_independence_validated=False,
                    runtime=dict(numpy=np.__version__, sklearn=sklearn.__version__, torch=torch.__version__),
                    code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sorted(Path(__file__).resolve().parent.glob("*.py"))})
    (args.output / "protocol.json").write_text(json.dumps(protocol,indent=2)+"\n")
    result = dict(protocol=protocol, seed=args.split_seed, models={}, validation_candidates={},
                  dataset_manifest_sha256=hashlib.sha256((args.dataset/"inputs.json").read_bytes()).hexdigest(),
                  partitions={name:dict(rows=int(mask.sum()),labels=dict(Counter(y[mask].astype(str).tolist())))
                              for name,mask in zip(("train","validation","test"),(tr,va,te))},
                  unknown_labels=int(np.count_nonzero(y<0)),
                  lifetime_independence_validated=False, deployment_ready=False)
    result["per_process_labels"] = {str(p):{str(label):int(np.count_nonzero((data["process"]==p)&(y==label)))
                                            for label in (-1,0,1)} for p in np.unique(data["process"])}
    fitted = {}
    details = {}
    started = time.monotonic()
    def archive(name,model,columns,description):
        bundle = dict(model=model, columns=np.asarray(columns,dtype=int),feature_names=FEATURE_NAMES,
                      requires_gen=False,deployment_ready=False)
        sv = predict(bundle,X[va])
        validation = validation_summary(sv,y[va],data["process"][va],data["nuisance"][va],wva,
                                        data["nuisance"][tr],y[tr])
        bundle["threshold"] = validation["threshold"]
        fitted[name], details[name] = bundle,description
        result["validation_candidates"][name] = validation
        with (models_dir/(name+".pkl")).open("wb") as out:
            pickle.dump(bundle,out)
        (args.output/"validation_candidates.json").write_text(json.dumps(result["validation_candidates"],indent=2,allow_nan=False)+"\n")
        print(json.dumps(dict(model=name,validation_auc=validation["auc"],
                              validation_rejection=1-validation["accidental"]["efficiency"],
                              validation_efficiency_span=validation["maximum_marginal_efficiency_span"],
                              selector=validation["selector_objective"])),flush=True)
    all_columns=np.arange(len(FEATURE_NAMES))
    for seed in args.seeds:
        for group,columns in feature_groups().items():
            model = make_pipeline(SimpleImputer(strategy="median",keep_empty_features=True),
                                   GradientBoostingClassifier(n_estimators=120,max_depth=2,min_samples_leaf=15,
                                       learning_rate=0.05,subsample=0.8,random_state=seed))
            model.fit(X[tr][:,columns],y[tr],gradientboostingclassifier__sample_weight=wtr)
            importance = {FEATURE_NAMES[i]:float(v) for i,v in zip(columns,model[-1].feature_importances_)}
            archive(f"{group}_s{seed}",model,columns,dict(kind="BDT ablation",feature_group=group,
                                                         fit_seed=seed,feature_importance=importance))
        for strength in args.uniform_strengths:
            name=f"uniform_bdt_{str(strength).replace('.','p')}_s{seed}"
            model,history=train_uniform_bdt(X[tr],y[tr],wtr,Z[tr],seed=seed,strength=strength,n_estimators=160)
            (args.output/(name+"_training.json")).write_text(json.dumps(history,indent=2,allow_nan=False)+"\n")
            archive(name,model,all_columns,dict(kind="neighborhood-uniformity BDT variant",strength=strength,fit_seed=seed))
        for strength in args.disco_strengths:
            name=f"nn_disco_{str(strength).replace('.','p')}_s{seed}"
            model,history=train_disco(X[tr],y[tr],wtr,Z[tr],X[va],y[va],wva,Z[va],
                                       seed=seed,strength=strength,epochs=args.epochs)
            (args.output/(name+"_training.json")).write_text(json.dumps(history,indent=2,allow_nan=False)+"\n")
            archive(name,model,all_columns,dict(kind="joint DisCo NN",strength=strength,fit_seed=seed))
        name=f"nn_adversarial_0p5_s{seed}"
        model,history=train_adversarial(X[tr],y[tr],wtr,Z[tr],X[va],y[va],wva,Z[va],
                                          seed=seed,strength=0.5,epochs=args.epochs)
        (args.output/(name+"_training.json")).write_text(json.dumps(history,indent=2,allow_nan=False)+"\n")
        archive(name,model,all_columns,dict(kind="adversarial NN reference",strength=0.5,fit_seed=seed))
    def best(prefix):
        candidates=[n for n in fitted if n.startswith(prefix)]
        return max(candidates,key=lambda n:result["validation_candidates"][n]["selector_objective"])
    reference_seed=args.seeds[0]
    representatives=[f"bdt_quality_s{reference_seed}",f"bdt_expanded_s{reference_seed}",
                     f"bdt_geometry_s{reference_seed}",best("uniform_bdt_"),best("nn_disco_"),best("nn_adversarial_")]
    result["representative_models"]=representatives
    result["selection_locked"]=True
    (args.output/"selection_locked.json").write_text(json.dumps(dict(representative_models=representatives,
         thresholds={n:fitted[n]["threshold"] for n in fitted},selection_uses_test=False),indent=2)+"\n")
    # Only now evaluate the reserved test population. No fit or selection follows.
    for name,bundle in fitted.items():
        scores=predict(bundle,X)
        value=evaluate(scores[te],y[te],data["nuisance"][te],data["ids"][te],data["process"][te],
                        data["stratum"][te],data["sampling_weight"][te],bundle["threshold"],data["nuisance"][tr],y[tr])
        value.update(details[name])
        value["validation_auc"]=result["validation_candidates"][name]["auc"]
        result["models"][name]=value
        np.savez_compressed(args.output/(name+"_scores.npz"),score=scores,split=split)
        print(json.dumps(dict(model=name,test_auc=value["auc"],genuine_efficiency=value["common_vertex"]["efficiency"],
                              accidental_rejection=value["background_rejection"],
                              joint_dcor=value["dependence"]["common_vertex"]["joint_dcor"])),flush=True)
    result["training_seconds"]=time.monotonic()-started
    (args.output/"metrics.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(complete=True,models=len(fitted),seconds=result["training_seconds"],
                          representatives=representatives)),flush=True)


if __name__ == "__main__":
    main()
