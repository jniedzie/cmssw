#!/usr/bin/env python3
"""Summarize double-traversing precision-only momentum closure."""

import argparse
import math

import numpy as np
import ROOT


def summarize(label, values):
    values = np.asarray([value for value in values if math.isfinite(value)])
    if not len(values):
        print(f"{label}: n=0")
        return
    q16, median, q84 = np.quantile(values, (0.16, 0.5, 0.84))
    print(f"{label}: n={len(values)} median={100 * median:+.3f}% half68={50 * (q84 - q16):.3f}%")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="+")
    args = parser.parse_args()
    events = ROOT.TChain("Events")
    for input_path in args.input:
        if not events.Add(input_path):
            raise RuntimeError(f"cannot add {input_path}")
    available = {branch.GetName() for branch in events.GetListOfBranches()}
    required = {
        "nShiftMuon", "ShiftMuon_directionalRefitUsedPrecisionHits",
        "ShiftMuon_genPartIdx", "nGenPart", "GenPart_pt", "GenPart_eta", "ShiftMuon_p",
        "ShiftMuon_directionalRefitPrecisionChi2", "ShiftMuon_directionalRefitPrecisionNdof",
        "ShiftMuon_directionalRefitPrecisionHits", "ShiftMuon_nPrecisionRefitStations",
        "ShiftMuon_precisionRefitLeverArm",
        "ShiftMuon_directionalRefitPrecisionValid", "ShiftMuon_directionalRefitAllHitsValid",
        "ShiftMuon_directionalRefitPrecisionRelativeToAllQoverP",
        "ShiftMuon_directionalRefitPrecisionTargetDca",
        "ShiftMuon_directionalRefitPrecisionFirstTargetPt",
        "ShiftMuon_directionalRefitPrecisionFirstTargetPz",
        "ShiftMuon_directionalRefitPrecisionSecondTargetPt",
        "ShiftMuon_directionalRefitPrecisionSecondTargetPz",
        "ShiftMuon_directionalRefitPrecisionSecondValid",
        "ShiftMuon_directionalRefitPrecisionSecondConverged",
        "ShiftMuon_directionalRefitPrecisionRelativeQoverPChange",
    }
    optional = {
        "ShiftMuon_simTruthMatched",
        "ShiftMuon_directionalRefitVacuumTargetPt",
        "ShiftMuon_directionalRefitVacuumTargetPz",
        "ShiftMuon_directionalRefitMaterialBoundaryValid",
        "ShiftMuon_directionalRefitMaterialBoundaryPt",
        "ShiftMuon_directionalRefitMaterialBoundaryPz",
        "ShiftMuon_directionalRefitMaterialPath",
        "ShiftMuon_directionalRefitPrecisionFirstRejected",
        "ShiftMuon_directionalRefitPrecisionSecondRejected",
        "ShiftMuon_directionalRefitPrecisionSelectedIteration",
    }
    if "ShiftMuon_topology" in available:
        required.add("ShiftMuon_topology")
    elif "ShiftMuon_quality" in available:
        required.add("ShiftMuon_quality")
    elif "ShiftMuon_source" in available:
        required.add("ShiftMuon_source")
    else:
        raise RuntimeError("missing ShiftMuon topology/legacy category branch")
    missing = required - available
    if missing:
        raise RuntimeError(f"missing branches: {', '.join(sorted(missing))}")
    events.SetBranchStatus("*", False)
    for branch in required | (optional & available):
        events.SetBranchStatus(branch, True)

    canonical = []
    first_iteration = []
    second_iteration = []
    iteration_pairs = []
    quality_diagnostics = []
    truth_rows = []
    both_endcap_rows = both_endcap_precision_valid = both_endcap_all_hits_valid = 0
    both_endcap_precision_relative = []
    both_endcap_precision_dca = []
    material_delta = []
    material_path = []
    material_boundary_fraction = []
    first_rejected = []
    second_rejected = []
    selected_iterations = []
    for event in events:
        for index in range(int(event.nShiftMuon)):
            if "ShiftMuon_topology" in available:
                is_both_endcaps = int(event.ShiftMuon_topology[index]) == 2
            elif "ShiftMuon_quality" in available:
                is_both_endcaps = int(event.ShiftMuon_quality[index]) == 3
            else:
                is_both_endcaps = int(event.ShiftMuon_source[index]) == 1
            if is_both_endcaps:
                both_endcap_rows += 1
                both_endcap_precision_valid += bool(event.ShiftMuon_directionalRefitPrecisionValid[index])
                both_endcap_all_hits_valid += bool(event.ShiftMuon_directionalRefitAllHitsValid[index])
                if event.ShiftMuon_directionalRefitPrecisionValid[index]:
                    both_endcap_precision_relative.append(
                        float(event.ShiftMuon_directionalRefitPrecisionRelativeToAllQoverP[index]))
                    both_endcap_precision_dca.append(float(event.ShiftMuon_directionalRefitPrecisionTargetDca[index]))
            if (not is_both_endcaps or not event.ShiftMuon_directionalRefitUsedPrecisionHits[index]):
                continue
            gen_index = int(event.ShiftMuon_genPartIdx[index])
            if gen_index < 0:
                continue
            gen_p = float(event.GenPart_pt[gen_index]) * math.cosh(float(event.GenPart_eta[gen_index]))
            reco_p = float(event.ShiftMuon_p[index])
            residual = reco_p / gen_p - 1
            canonical.append(residual)
            if "ShiftMuon_directionalRefitPrecisionSelectedIteration" in available:
                selected_iterations.append(
                    int(event.ShiftMuon_directionalRefitPrecisionSelectedIteration[index]))
            if "ShiftMuon_directionalRefitPrecisionFirstRejected" in available:
                first_rejected.append(int(event.ShiftMuon_directionalRefitPrecisionFirstRejected[index]))
                if event.ShiftMuon_directionalRefitPrecisionSecondValid[index]:
                    second_rejected.append(int(event.ShiftMuon_directionalRefitPrecisionSecondRejected[index]))
            if "ShiftMuon_directionalRefitVacuumTargetPt" in available:
                vacuum_p = math.hypot(float(event.ShiftMuon_directionalRefitVacuumTargetPt[index]),
                                      float(event.ShiftMuon_directionalRefitVacuumTargetPz[index]))
                if vacuum_p > 0:
                    material_delta.append(reco_p / vacuum_p - 1)
                boundary_valid = bool(event.ShiftMuon_directionalRefitMaterialBoundaryValid[index])
                material_boundary_fraction.append(boundary_valid)
                if boundary_valid:
                    material_path.append(float(event.ShiftMuon_directionalRefitMaterialPath[index]))
            first_p = math.hypot(float(event.ShiftMuon_directionalRefitPrecisionFirstTargetPt[index]),
                                 float(event.ShiftMuon_directionalRefitPrecisionFirstTargetPz[index]))
            first_iteration.append(first_p / gen_p - 1)
            if event.ShiftMuon_directionalRefitPrecisionSecondValid[index]:
                second_p = math.hypot(float(event.ShiftMuon_directionalRefitPrecisionSecondTargetPt[index]),
                                      float(event.ShiftMuon_directionalRefitPrecisionSecondTargetPz[index]))
                second_residual = second_p / gen_p - 1
                second_iteration.append(second_residual)
                iteration_pairs.append((
                    first_p / gen_p - 1,
                    second_residual,
                    float(event.ShiftMuon_directionalRefitPrecisionRelativeQoverPChange[index]),
                ))
            ndof = float(event.ShiftMuon_directionalRefitPrecisionNdof[index])
            quality_diagnostics.append((
                residual,
                float(event.ShiftMuon_directionalRefitPrecisionChi2[index]) / ndof if ndof > 0 else math.inf,
                int(event.ShiftMuon_directionalRefitPrecisionHits[index]),
                int(event.ShiftMuon_nPrecisionRefitStations[index]),
                float(event.ShiftMuon_precisionRefitLeverArm[index]),
            ))
            if ("ShiftMuon_simTruthMatched" in available and
                    event.ShiftMuon_simTruthMatched[index]):
                truth_rows.append(index)

    print(f"events={events.GetEntries()} bothEndcaps={both_endcap_rows} "
          f"precisionValid={both_endcap_precision_valid} allHitsValid={both_endcap_all_hits_valid} "
          f"selected={len(canonical)} truthMatched={len(truth_rows)}")
    if both_endcap_precision_relative:
        print(f"precision/all qop median={np.median(both_endcap_precision_relative):.4f} "
              f"targetDca median={np.median(both_endcap_precision_dca):.2f} cm")
    summarize("canonical", canonical)
    summarize("first iteration", first_iteration)
    summarize("valid second iteration", second_iteration)
    if iteration_pairs:
        first_residuals = np.asarray([row[0] for row in iteration_pairs])
        second_residuals = np.asarray([row[1] for row in iteration_pairs])
        deltas = second_residuals - first_residuals
        summarize("paired second-minus-first residual", deltas)
        if len(deltas) > 1 and np.std(deltas) > 0 and np.std(first_residuals) > 0:
            print(f"correlation(second-first, first residual)="
                  f"{np.corrcoef(deltas, first_residuals)[0, 1]:+.3f}")
        for maximum_change in (0.01, 0.02, 0.05, 0.10, 0.50):
            guarded = [second if relative_change <= maximum_change else first
                       for first, second, relative_change in iteration_pairs]
            accepted = sum(relative_change <= maximum_change
                           for _, _, relative_change in iteration_pairs)
            summarize(f"q/p guard<={maximum_change:.2f} ({accepted}/{len(iteration_pairs)} second)", guarded)
    if first_rejected:
        print(f"outliers first={sum(value > 0 for value in first_rejected)}/{len(first_rejected)} "
              f"second={sum(value > 0 for value in second_rejected)}/{len(second_rejected)} "
              f"removed={sum(first_rejected)}/{sum(second_rejected)}")
    if selected_iterations:
        print(f"selected iteration 1={selected_iterations.count(1)} 2={selected_iterations.count(2)}")
    if material_boundary_fraction:
        print(f"target material boundary valid={sum(material_boundary_fraction)}/"
              f"{len(material_boundary_fraction)}")
        summarize("target material / vacuum momentum", material_delta)
        if material_path:
            print(f"target material path median={np.median(material_path):.2f} cm")
    for maximum in (1, 2, 5, 10, 20, 50):
        summarize(f"canonical chi2/ndof<{maximum}",
                  [row[0] for row in quality_diagnostics if row[1] < maximum])
    for minimum in (4, 6, 8, 10, 12):
        summarize(f"canonical hits>={minimum}",
                  [row[0] for row in quality_diagnostics if row[2] >= minimum])


if __name__ == "__main__":
    main()
