#!/usr/bin/env python3
"""Cut-level process/support and unknown-label audit; never fit or tune models."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from evaluation import binned_efficiencies, fraction


def unknown_bounds(y, passed):
    """Count-only extrema over every possible assignment of unknown labels."""
    good = int(np.count_nonzero(y == 1))
    good_pass = int(np.count_nonzero((y == 1) & passed))
    bad = int(np.count_nonzero(y == 0))
    bad_reject = int(np.count_nonzero((y == 0) & ~passed))
    unknown_pass = int(np.count_nonzero((y < 0) & passed))
    unknown_reject = int(np.count_nonzero((y < 0) & ~passed))
    def ratio(n,d):
        return n/d if d else None
    def bounds(low_n,low_d,high_n,high_d):
        low,high=ratio(low_n,low_d),ratio(high_n,high_d)
        # Ignore assignments where the class is absent. If only one corner
        # defines a possible class, every defined assignment has that value.
        if low is None and high is not None:
            low=high
        if high is None and low is not None:
            high=low
        return [low,high]
    return dict(known_common=good, known_common_passed=good_pass,
                known_accidental=bad, known_accidental_rejected=bad_reject,
                unknown_passed=unknown_pass, unknown_rejected=unknown_reject,
                genuine_efficiency_extrema=bounds(good_pass,good+unknown_reject,
                                                  good_pass+unknown_pass,good+unknown_pass),
                accidental_rejection_extrema=bounds(bad_reject,bad+unknown_pass,
                                                    bad_reject+unknown_reject,bad+unknown_reject),
                interpretation="Count-only identification bounds in this sample; not a physical uncertainty or correction")


def plot_profiles(report, metrics, model, output):
    fig, axes=plt.subplots(3,3,figsize=(12,9))
    for row,process in enumerate(("all","jpsi","dy")):
        cells=report["models"][model]["per_process"][process]["flatness"]
        for col,variable in enumerate(("mass","vertex_z","vertex_radius")):
            ax=axes[row,col]
            entries=cells[variable]["bins"]
            for j,e in enumerate(entries):
                if e["n"]:
                    value=e["efficiency"]
                    ax.errorbar(j,value,yerr=[[value-e["ci95"][0]],[e["ci95"][1]-value]],fmt="o",capsize=3,
                                color="#1565c0" if e["supported"] else "#777777")
                if not e["supported"]:
                    ax.axvspan(j-.4,j+.4,color="#eeeeee",zorder=-1)
            edges=cells[variable]["interior_edges"]
            if variable=="vertex_z":
                edges=[value/100 for value in edges]
            labels=[f"<{edges[0]:.4g}",f"{edges[0]:.4g}--{edges[1]:.4g}",f">={edges[1]:.4g}"]
            ax.set_xticks(range(len(entries)),[f"{label}\n(n={e['n']})" for label,e in zip(labels,entries)])
            ax.set_ylim(0,1.04)
            process_label={"all":"All SM","jpsi":"J/psi","dy":"DY"}[process]
            variable_label={"mass":"Mass [GeV]","vertex_z":"Vertex z [m]","vertex_radius":"Vertex radius [cm]"}[variable]
            ax.set_title(f"{process_label}: {variable_label}")
            ax.set_ylabel("Genuine-pair acceptance")
            ax.grid(alpha=.2)
    fig.suptitle("Uniformity BDT: held-out acceptance by process\n"
                 "Train-defined boundaries; genuine pairs with reliable labels",fontsize=13)
    fig.text(.5,.02,"Grey bins: fewer than 20 genuine test pairs, unvalidated. 95% Wilson intervals approximate pair independence.\n"
             "Only reliably labelled SM pairs; no lifetime or normalization closure.",ha="center",fontsize=8)
    fig.tight_layout(rect=(0,.065,1,.93))
    for ext in ("pdf","png"):
        fig.savefig(output/("acceptance_profiles."+ext),dpi=160)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset",required=True,type=Path)
    parser.add_argument("--results",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    manifest_path=args.dataset/"inputs.json"
    manifest=json.loads(manifest_path.read_text())
    metrics=json.loads((args.results/"metrics.json").read_text())
    if manifest.get("sample_kind")!="simulation_only" or metrics.get("selection_locked") is not True:
        raise ValueError("Only locked SM comparisons may be audited")
    if hashlib.sha256(manifest_path.read_bytes()).hexdigest()!=metrics["dataset_manifest_sha256"]:
        raise ValueError("Dataset manifest differs")
    for filename,key in (("reco.npz","reco_sha256"),("labels.npz","labels_sha256")):
        if hashlib.sha256((args.dataset/filename).read_bytes()).hexdigest()!=manifest[key]:
            raise ValueError("Dataset digest differs")
    data=dict(np.load(args.dataset/"reco.npz",allow_pickle=False))
    y=np.load(args.dataset/"labels.npz",allow_pickle=False)["y"]
    partitions=dict(np.load(args.results/"splits.npz",allow_pickle=False))
    split=partitions["split"]
    if not np.array_equal(partitions["ids"],data["ids"]) or not np.array_equal(partitions["process"],data["process"]):
        raise ValueError("Saved event identities differ")
    test=split==2
    train=(split==0)&(y>=0)
    report=dict(dataset_manifest_sha256=metrics["dataset_manifest_sha256"],
                scope="Held-out SM cut diagnostics; no model changes; unknown bounds are not physics uncertainties",
                models={}, deployment_ready=False,lifetime_independence_validated=False)
    for name in metrics["representative_models"]:
        if Path(name).name!=name:
            raise ValueError("Invalid score filename")
        scores=dict(np.load(args.results/(name+"_scores.npz"),allow_pickle=False))
        if not np.array_equal(scores["split"],split) or scores["score"].shape!=y.shape or not np.isfinite(scores["score"]).all():
            raise ValueError("Score alignment differs")
        passed=scores["score"]>=metrics["models"][name]["threshold"]
        groups={}
        for process in ("all","jpsi","dy","qcd"):
            mask=test if process=="all" else test&(data["process"]==process)
            labelled=mask&(y>=0)
            flat=binned_efficiencies(data["nuisance"][labelled],y[labelled],passed[labelled],
                                     data["nuisance"][train],y[train])
            unsupported=sum(e["n"] for e in flat["joint_mass_z_radius"] if not e["supported"])
            groups[process]=dict(common_vertex=fraction(mask&(y==1),passed),
                                 accidental=fraction(mask&(y==0),passed),unknown_retention=fraction(mask&(y<0),passed),
                                 bounds=unknown_bounds(y[mask],passed[mask]),flatness=flat,
                                 genuine_test_pairs_in_unsupported_joint_cells=unsupported)
        report["models"][name]=dict(per_process=groups)
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/"acceptance_audit.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    uniform=next((n for n in metrics["representative_models"] if n.startswith("uniform_bdt_")),None)
    if uniform:
        plot_profiles(report,metrics,uniform,args.output)
    print(json.dumps(dict(output=str(args.output),models=list(report["models"]),
                          total_test_pairs=int(test.sum()),unknown_test_pairs=int((test&(y<0)).sum()))),flush=True)


if __name__=="__main__":
    main()
