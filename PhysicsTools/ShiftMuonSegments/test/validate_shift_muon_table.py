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
        "ShiftMuon_nDTHits", "ShiftMuon_nCSCHits", "ShiftMuon_nRPCHits",
        "ShiftMuon_nGEMHits", "ShiftMuon_nME0Hits", "ShiftMuon_detectorMask",
        "ShiftMuon_nHitsPlusEndcap", "ShiftMuon_nHitsBarrel", "ShiftMuon_nHitsMinusEndcap",
        "ShiftMuon_nHitsNearEndcap", "ShiftMuon_nHitsFarEndcap",
        "ShiftMuon_nStationsNearEndcap", "ShiftMuon_nStationsBarrel",
        "ShiftMuon_nStationsFarEndcap", "ShiftMuon_entryExitValid",
        "ShiftMuon_entryRegion", "ShiftMuon_exitRegion",
        "ShiftMuon_entrySubdetector", "ShiftMuon_exitSubdetector",
        "ShiftMuon_entryX", "ShiftMuon_entryY", "ShiftMuon_entryZ", "ShiftMuon_entryR",
        "ShiftMuon_exitX", "ShiftMuon_exitY", "ShiftMuon_exitZ", "ShiftMuon_exitR",
        "ShiftMuon_hitSpan",
        "nShiftDimuonVertex", "ShiftDimuonVertex_constrainedValid",
        "ShiftDimuonVertex_constrainedMass", "ShiftDimuonVertex_constrainedPt",
        "ShiftDimuonVertex_constrainedPz", "ShiftDimuonVertex_constrainedEta",
        "ShiftDimuonVertex_constrainedPhi", "ShiftDimuonVertex_constrainedVx",
        "ShiftDimuonVertex_constrainedVy", "ShiftDimuonVertex_constrainedVz",
        "ShiftDimuonVertex_isCosmicCosmic", "ShiftDimuonVertex_isCosmicDSA",
        "ShiftDimuonVertex_isCosmicTraversing", "ShiftDimuonVertex_isCosmicDoubleTraversing",
        "ShiftDimuonVertex_isDSADSA", "ShiftDimuonVertex_isDSATraversing",
        "ShiftDimuonVertex_isDSADoubleTraversing", "ShiftDimuonVertex_isTraversingTraversing",
        "ShiftDimuonVertex_isTraversingDoubleTraversing",
        "ShiftDimuonVertex_isDoubleTraversingDoubleTraversing",
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
    dimuons = 0
    constrained_dimuons = 0
    quality_flags = (
        "isCosmicCosmic", "isCosmicDSA", "isCosmicTraversing", "isCosmicDoubleTraversing",
        "isDSADSA", "isDSATraversing", "isDSADoubleTraversing", "isTraversingTraversing",
        "isTraversingDoubleTraversing", "isDoubleTraversingDoubleTraversing",
    )
    for event in events:
        for index in range(int(event.nShiftMuon)):
            muons += 1
            quality = int(event.ShiftMuon_quality[index])
            if quality not in quality_counts:
                raise RuntimeError(f"invalid quality value {quality}")
            quality_counts[quality] += 1
            constrained_valid += bool(event.ShiftMuon_constrainedValid[index])
            subdetector_counts = (
                int(event.ShiftMuon_nDTHits[index]), int(event.ShiftMuon_nCSCHits[index]),
                int(event.ShiftMuon_nRPCHits[index]), int(event.ShiftMuon_nGEMHits[index]),
                int(event.ShiftMuon_nME0Hits[index]),
            )
            region_count = (
                int(event.ShiftMuon_nHitsPlusEndcap[index]) +
                int(event.ShiftMuon_nHitsBarrel[index]) +
                int(event.ShiftMuon_nHitsMinusEndcap[index])
            )
            if region_count != sum(subdetector_counts):
                raise RuntimeError(
                    f"muon {muons - 1} has inconsistent detector/region hit counts: "
                    f"{subdetector_counts} versus {region_count}"
                )
            expected_mask = sum((1 << bit) for bit, count in enumerate(subdetector_counts) if count)
            if int(event.ShiftMuon_detectorMask[index]) != expected_mask:
                raise RuntimeError(f"muon {muons - 1} has inconsistent detectorMask")
            endpoints_valid = bool(event.ShiftMuon_entryExitValid[index])
            entry_region = int(event.ShiftMuon_entryRegion[index])
            exit_region = int(event.ShiftMuon_exitRegion[index])
            if endpoints_valid and (entry_region not in range(3) or exit_region not in range(3)):
                raise RuntimeError(
                    f"muon {muons - 1} has invalid entry/exit regions {entry_region}/{exit_region}"
                )
            if float(event.ShiftMuon_hitSpan[index]) < 0.0:
                raise RuntimeError(f"muon {muons - 1} has negative hit span")
        for index in range(int(event.nShiftDimuonVertex)):
            dimuons += 1
            constrained_dimuons += bool(event.ShiftDimuonVertex_constrainedValid[index])
            active_flags = sum(
                bool(getattr(event, f"ShiftDimuonVertex_{name}")[index]) for name in quality_flags
            )
            if active_flags != 1:
                raise RuntimeError(f"dimuon {dimuons - 1} has {active_flags} active quality-pair flags")

    print(f"events={events.GetEntries()} muons={muons} constrainedValid={constrained_valid}/{muons}")
    print("quality=" + ",".join(f"{key}:{value}" for key, value in quality_counts.items()))
    print(f"dimuons={dimuons} constrainedValid={constrained_dimuons}/{dimuons}")


if __name__ == "__main__":
    main()
