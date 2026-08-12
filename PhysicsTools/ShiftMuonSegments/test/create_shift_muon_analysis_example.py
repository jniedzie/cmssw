#!/usr/bin/env python3
"""Extract one populated event with the complete SHIFT analysis schema."""

import argparse
from pathlib import Path

import ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="NanoAOD file produced by the current ShiftMuonTableProducer")
    parser.add_argument("output", nargs="?", default="shift_muon_analysis_example.root")
    args = parser.parse_args()

    source = ROOT.TFile.Open(args.input)
    if not source or source.IsZombie():
        raise RuntimeError(f"cannot open {args.input}")
    events = source.Get("Events")
    if not events:
        raise RuntimeError("input has no Events tree")

    # Keep the complete ShiftMuon and ShiftDimuonVertex tables, along with the
    # generator information needed to prototype resolution studies.  This is
    # deliberately extracted from a real output instead of manually copying a
    # branch list, so established fields such as vx/vy/vz and chi2 cannot be
    # accidentally omitted from the example again.
    events.SetBranchStatus("*", 0)
    for pattern in (
        "run", "luminosityBlock", "event",
        "nShiftMuon", "ShiftMuon_*",
        "nShiftDimuonVertex", "ShiftDimuonVertex_*",
        "nGenPart", "GenPart_*",
    ):
        events.SetBranchStatus(pattern, 1)

    # Allow this helper to consume the transitional output used to develop the
    # final schema while ensuring that its output has only the final branches.
    for branch in (
        "ShiftMuon_unconstrainedPt", "ShiftMuon_unconstrainedEta",
        "ShiftMuon_unconstrainedPhi", "ShiftMuon_unconstrainedMass",
        "ShiftMuon_unconstrainedP", "ShiftMuon_unconstrainedPx",
        "ShiftMuon_unconstrainedPy", "ShiftMuon_unconstrainedPz",
        "ShiftMuon_source", "ShiftMuon_isTraversing", "ShiftMuon_isFullLeverArm",
        "ShiftMuon_isDSA", "ShiftMuon_isCosmic",
    ):
        events.SetBranchStatus(branch, 0)

    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output = ROOT.TFile(str(output_path), "RECREATE")
    populated_entry = next((index for index, event in enumerate(events) if int(event.nShiftMuon) > 0), None)
    if populated_entry is None:
        raise RuntimeError("input contains no event with a ShiftMuon")
    example = events.CloneTree(0)
    events.GetEntry(populated_entry)
    example.Fill()
    example.SetName("Events")
    example.Write()
    output.Close()
    source.Close()
    print(output_path)


if __name__ == "__main__":
    main()
