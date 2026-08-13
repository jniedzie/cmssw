#!/usr/bin/env python3
"""Relate ShiftMuon reconstruction quality to detector geometry diagnostics."""

import argparse
from collections import Counter, defaultdict

import ROOT


QUALITY_NAMES = {
    0: "cosmic",
    1: "DSA",
    2: "traversing",
    3: "double_traversing",
}
REGION_NAMES = {-1: "unknown", 0: "near_endcap", 1: "barrel", 2: "far_endcap"}
SUBDETECTORS = ("DT", "CSC", "RPC", "GEM", "ME0")


def topology(event, index):
    near = int(event.ShiftMuon_nHitsNearEndcap[index])
    barrel = int(event.ShiftMuon_nHitsBarrel[index])
    far = int(event.ShiftMuon_nHitsFarEndcap[index])
    if near and far:
        return "both_endcaps"
    if near and barrel:
        return "near_plus_barrel"
    if near:
        return "near_only"
    if barrel and far:
        return "barrel_plus_far"
    if barrel:
        return "barrel_only"
    if far:
        return "far_only"
    return "unclassified"


def detector_pattern(mask):
    names = [name for bit, name in enumerate(SUBDETECTORS) if mask & (1 << bit)]
    return "+".join(names) if names else "none"


def eta_bin(value):
    value = abs(value)
    if value < 2.4:
        return "<2.4"
    if value < 3.0:
        return "2.4-3"
    if value < 4.0:
        return "3-4"
    return ">=4"


def momentum_bin(value):
    if value < 20.0:
        return "<20"
    if value < 50.0:
        return "20-50"
    if value < 100.0:
        return "50-100"
    return ">=100"


def median(values):
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return 0.5 * (ordered[middle - 1] + ordered[middle])


def print_matrix(title, counts, columns):
    print(f"\n{title}")
    print(f"{'quality':20s}" + "".join(f"{column:>20s}" for column in columns))
    for quality, name in QUALITY_NAMES.items():
        print(f"{name:20s}" + "".join(f"{counts[(quality, column)]:20d}" for column in columns))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="NanoAOD containing the ShiftMuon geometry diagnostics")
    args = parser.parse_args()

    root_file = ROOT.TFile.Open(args.input)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open {args.input}")
    events = root_file.Get("Events")
    if not events:
        raise RuntimeError("missing Events tree")

    topology_counts = Counter()
    exit_counts = Counter()
    detector_counts = Counter()
    eta_counts = Counter()
    momentum_counts = Counter()
    values = defaultdict(lambda: defaultdict(list))

    for event in events:
        for index in range(int(event.nShiftMuon)):
            quality = int(event.ShiftMuon_quality[index])
            topology_name = topology(event, index)
            exit_name = REGION_NAMES.get(int(event.ShiftMuon_exitRegion[index]), "invalid")
            detectors = detector_pattern(int(event.ShiftMuon_detectorMask[index]))
            topology_counts[(quality, topology_name)] += 1
            exit_counts[(quality, exit_name)] += 1
            detector_counts[(quality, detectors)] += 1
            eta_counts[(quality, eta_bin(float(event.ShiftMuon_eta[index])))] += 1
            momentum_counts[(quality, momentum_bin(float(event.ShiftMuon_p[index])))] += 1
            values[quality]["absEta"].append(abs(float(event.ShiftMuon_eta[index])))
            values[quality]["p"].append(float(event.ShiftMuon_p[index]))
            values[quality]["hitSpan"].append(float(event.ShiftMuon_hitSpan[index]))

    topology_columns = (
        "near_only", "near_plus_barrel", "both_endcaps", "barrel_only",
        "barrel_plus_far", "far_only", "unclassified",
    )
    print(f"events={events.GetEntries()}")
    print_matrix("Recorded-hit topology", topology_counts, topology_columns)
    print_matrix("Last recorded hit region", exit_counts, ("near_endcap", "barrel", "far_endcap", "unknown"))
    detector_columns = sorted({key[1] for key in detector_counts})
    print_matrix("Muon subdetector combinations", detector_counts, detector_columns)
    print_matrix("Absolute eta bins", eta_counts, ("<2.4", "2.4-3", "3-4", ">=4"))
    print_matrix("Momentum bins [GeV]", momentum_counts, ("<20", "20-50", "50-100", ">=100"))

    print("\nMedians by reconstruction quality")
    print(f"{'quality':20s}{'|eta|':>12s}{'p [GeV]':>12s}{'hit span [cm]':>18s}")
    for quality, name in QUALITY_NAMES.items():
        print(
            f"{name:20s}{median(values[quality]['absEta']):12.3f}"
            f"{median(values[quality]['p']):12.3f}{median(values[quality]['hitSpan']):18.1f}"
        )


if __name__ == "__main__":
    main()
