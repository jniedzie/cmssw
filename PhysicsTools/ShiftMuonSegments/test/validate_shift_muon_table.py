#!/usr/bin/env python3
"""Validate the compact dual-momentum ShiftMuon schema."""

import argparse

import ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    args = parser.parse_args()
    root_file = ROOT.TFile.Open(args.input)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open {args.input}")
    events = root_file.Get("Events")
    if not events:
        raise RuntimeError("missing Events tree")

    available = {branch.GetName() for branch in events.GetListOfBranches()}
    required = {
        "nShiftMuon", "ShiftMuon_pt", "ShiftMuon_eta", "ShiftMuon_phi",
        "ShiftMuon_vx", "ShiftMuon_vy", "ShiftMuon_vz", "ShiftMuon_chi2",
        "ShiftMuon_ndof", "ShiftMuon_constrainedPt", "ShiftMuon_constrainedEta",
        "ShiftMuon_constrainedPhi", "ShiftMuon_constrainedVx",
        "ShiftMuon_constrainedVy", "ShiftMuon_constrainedVz",
        "ShiftMuon_constrainedValid", "ShiftMuon_constrainedStatus",
        "ShiftMuon_constrainedTargetChi2", "ShiftMuon_quality", "ShiftMuon_sourceIndex",
    }
    forbidden = {
        "ShiftMuon_unconstrainedPt", "ShiftMuon_source", "ShiftMuon_isTraversing",
        "ShiftMuon_isFullLeverArm", "ShiftMuon_isDSA", "ShiftMuon_isCosmic",
    }
    missing = sorted(required - available)
    redundant = sorted(forbidden & available)
    if missing:
        raise RuntimeError(f"missing branches: {', '.join(missing)}")
    if redundant:
        raise RuntimeError(f"redundant branches remain: {', '.join(redundant)}")

    quality_counts = {index: 0 for index in range(4)}
    constrained_valid = 0
    muons = 0
    for event in events:
        for index in range(int(event.nShiftMuon)):
            muons += 1
            quality = int(event.ShiftMuon_quality[index])
            if quality not in quality_counts:
                raise RuntimeError(f"invalid quality value {quality}")
            quality_counts[quality] += 1
            constrained_valid += bool(event.ShiftMuon_constrainedValid[index])

    print(f"events={events.GetEntries()} muons={muons} constrainedValid={constrained_valid}/{muons}")
    print("quality=" + ",".join(f"{key}:{value}" for key, value in quality_counts.items()))


if __name__ == "__main__":
    main()
