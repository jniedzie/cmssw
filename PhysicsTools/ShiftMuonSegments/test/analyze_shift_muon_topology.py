#!/usr/bin/env python3
"""Relate ShiftMuon measured topology to reconstruction provenance."""

import argparse
from collections import Counter, defaultdict

import ROOT


ALGORITHM_NAMES = {
    0: "DSA",
    1: "strict_traversing",
    2: "cosmic",
}
TOPOLOGY_NAMES = {
    0: "near_endcap_only",
    1: "near_endcap_and_barrel",
    2: "both_endcaps",
    3: "far_endcap_only",
    4: "unclassified",
}
REGION_NAMES = {-1: "unknown", 0: "near_endcap", 1: "barrel", 2: "far_endcap"}
SUBDETECTORS = ("DT", "CSC", "RPC", "GEM", "ME0")


def derived_topology(event, index):
    if not bool(event.ShiftMuon_orientationValid[index]):
        return 4
    near = int(event.ShiftMuon_nHitsNearEndcap[index])
    barrel = int(event.ShiftMuon_nHitsBarrel[index])
    far = int(event.ShiftMuon_nHitsFarEndcap[index])
    if near and far:
        return 2
    if near and barrel:
        return 1
    if near:
        return 0
    if far and not barrel:
        return 3
    return 4


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
    print(f"{'reco algorithm':20s}" + "".join(f"{column:>28s}" for column in columns))
    for algorithm, name in ALGORITHM_NAMES.items():
        print(f"{name:20s}" + "".join(f"{counts[(algorithm, column)]:28d}" for column in columns))


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
            algorithm = int(event.ShiftMuon_recoAlgorithm[index])
            topology_value = int(event.ShiftMuon_topology[index])
            expected_topology = derived_topology(event, index)
            if topology_value != expected_topology:
                raise RuntimeError(
                    f"stored topology {topology_value} differs from derived value {expected_topology}"
                )
            topology_name = TOPOLOGY_NAMES[topology_value]
            exit_name = REGION_NAMES.get(int(event.ShiftMuon_exitRegion[index]), "invalid")
            detectors = detector_pattern(int(event.ShiftMuon_detectorMask[index]))
            topology_counts[(algorithm, topology_name)] += 1
            exit_counts[(algorithm, exit_name)] += 1
            detector_counts[(algorithm, detectors)] += 1
            eta_counts[(algorithm, eta_bin(float(event.ShiftMuon_eta[index])))] += 1
            momentum_counts[(algorithm, momentum_bin(float(event.ShiftMuon_p[index])))] += 1
            values[algorithm]["absEta"].append(abs(float(event.ShiftMuon_eta[index])))
            values[algorithm]["p"].append(float(event.ShiftMuon_p[index]))
            values[algorithm]["hitSpan"].append(float(event.ShiftMuon_hitSpan[index]))

    topology_columns = (
        "near_endcap_only", "near_endcap_and_barrel", "both_endcaps",
        "far_endcap_only", "unclassified",
    )
    print(f"events={events.GetEntries()}")
    print_matrix("Recorded-hit topology", topology_counts, topology_columns)
    print_matrix("Last recorded hit region", exit_counts, ("near_endcap", "barrel", "far_endcap", "unknown"))
    detector_columns = sorted({key[1] for key in detector_counts})
    print_matrix("Muon subdetector combinations", detector_counts, detector_columns)
    print_matrix("Absolute eta bins", eta_counts, ("<2.4", "2.4-3", "3-4", ">=4"))
    print_matrix("Momentum bins [GeV]", momentum_counts, ("<20", "20-50", "50-100", ">=100"))

    print("\nMedians by reconstruction algorithm")
    print(f"{'reco algorithm':20s}{'|eta|':>12s}{'p [GeV]':>12s}{'hit span [cm]':>18s}")
    for algorithm, name in ALGORITHM_NAMES.items():
        print(
            f"{name:20s}{median(values[algorithm]['absEta']):12.3f}"
            f"{median(values[algorithm]['p']):12.3f}{median(values[algorithm]['hitSpan']):18.1f}"
        )


if __name__ == "__main__":
    main()
