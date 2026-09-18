#!/usr/bin/env python3

"""Verify that enabling a simulation diagnostic did not perturb Geant4 output."""

import argparse
import hashlib
import json
import struct

from DataFormats.FWLite import Events, Handle


SIMHIT_INSTANCES = (
    "BCM1FHits", "BHMHits", "BSCHits", "CTPPSPixelHits", "CTPPSTimingHits",
    "MuonCSCHits", "MuonDTHits", "MuonGEMHits", "MuonME0Hits", "MuonRPCHits",
    "PLTHits", "TotemHitsRP", "TotemHitsT1",
    "TrackerHitsPixelBarrelHighTof", "TrackerHitsPixelBarrelLowTof",
    "TrackerHitsPixelEndcapHighTof", "TrackerHitsPixelEndcapLowTof",
    "TrackerHitsTECHighTof", "TrackerHitsTECLowTof",
    "TrackerHitsTIBHighTof", "TrackerHitsTIBLowTof",
    "TrackerHitsTIDHighTof", "TrackerHitsTIDLowTof",
    "TrackerHitsTOBHighTof", "TrackerHitsTOBLowTof",
)


def _product(event, type_name, instance=""):
    handle = Handle(type_name)
    event.getByLabel("g4SimHits", instance, handle)
    if not handle.isValid():
        return None
    return handle.product()


def _four_vector(value):
    return value.px(), value.py(), value.pz(), value.e()


def _simtrack_record(track):
    values = [
        track.eventId().rawId(), track.trackId(), track.type(), track.vertIndex(),
        track.genpartIndex(), int(track.isPrimary()), int(track.crossedBoundary()),
        int(track.isFromBackScattering()), track.charge(), *_four_vector(track.momentum()),
    ]
    if track.crossedBoundary():
        position = track.getPositionAtBoundary()
        values.extend((track.getIDAtBoundary(), position.x(), position.y(), position.z(), position.t()))
        values.extend(_four_vector(track.getMomentumAtBoundary()))
    return repr(tuple(values)).encode("ascii")


def _simvertex_record(vertex):
    position = vertex.position()
    return repr((
        vertex.eventId().rawId(), vertex.vertexId(), vertex.parentIndex(),
        vertex.processType(), position.x(), position.y(), position.z(), position.t(),
    )).encode("ascii")


def _simhit_record(hit):
    entry = hit.entryPoint()
    exit_point = hit.exitPoint()
    direction = hit.localDirection()
    momentum = hit.momentumAtEntry()
    return struct.pack(
        "<IIIiHH17f",
        int(hit.eventId().rawId()), int(hit.trackId()), int(hit.detUnitId()),
        int(hit.particleType()), int(hit.processType()), int(hit.hitProdType()),
        float(entry.x()), float(entry.y()), float(entry.z()),
        float(exit_point.x()), float(exit_point.y()), float(exit_point.z()),
        float(direction.x()), float(direction.y()), float(direction.z()),
        float(hit.pabs()), float(hit.energyLoss()), float(hit.timeOfFlight()),
        float(hit.thetaAtEntry()), float(hit.phiAtEntry()),
        float(momentum.x()), float(momentum.y()), float(momentum.z()),
    )


def _digest(records):
    digest = hashlib.sha256()
    count = 0
    for record in records:
        digest.update(record)
        count += 1
    return {"count": count, "sha256": digest.hexdigest()}


def _product_digest(product, serializer):
    if product is None:
        return {"present": False}
    result = _digest(serializer(value) for value in product)
    result["present"] = True
    return result


def _event_summary(event):
    event_id = event.eventAuxiliary().id()
    products = {
        "SimTracks": _product_digest(_product(event, "std::vector<SimTrack>"), _simtrack_record),
        "SimVertices": _product_digest(_product(event, "std::vector<SimVertex>"), _simvertex_record),
    }
    for instance in SIMHIT_INSTANCES:
        products[instance] = _product_digest(
            _product(event, "std::vector<PSimHit>", instance), _simhit_record
        )
    return {
        "run": int(event_id.run()),
        "lumi": int(event.eventAuxiliary().luminosityBlock()),
        "event": int(event_id.event()),
        "products": products,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("--output")
    args = parser.parse_args()
    if args.max_events == 0 or args.max_events < -1:
        parser.error("--max-events must be -1 or positive")

    def summaries(path):
        result = []
        for index, event in enumerate(Events(path)):
            if args.max_events > 0 and index >= args.max_events:
                break
            result.append(_event_summary(event))
        return result

    # FWLite uses process-wide ROOT state. Fully consume one file before opening
    # the other so this remains reliable for files with different branch lists.
    reference_summaries = summaries(args.reference)
    candidate_summaries = summaries(args.candidate)
    comparisons = []
    valid = True
    for reference_summary, candidate_summary in zip(reference_summaries, candidate_summaries):
        mismatches = [
            name for name in reference_summary["products"]
            if reference_summary["products"][name] != candidate_summary["products"].get(name)
        ]
        same_id = all(reference_summary[key] == candidate_summary[key] for key in ("run", "lumi", "event"))
        valid = valid and same_id and not mismatches
        comparisons.append({
            "reference_event": {key: reference_summary[key] for key in ("run", "lumi", "event")},
            "candidate_event": {key: candidate_summary[key] for key in ("run", "lumi", "event")},
            "matching_event_id": same_id,
            "mismatched_products": mismatches,
            "products": candidate_summary["products"],
        })

    report = {
        "reference": args.reference,
        "candidate": args.candidate,
        "events_compared": len(comparisons),
        "valid": (
            valid and bool(comparisons)
            and len(reference_summaries) == len(candidate_summaries)
        ),
        "comparisons": comparisons,
    }
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as output:
            output.write(text + "\n")
    print(text)
    raise SystemExit(0 if report["valid"] else 1)


if __name__ == "__main__":
    main()
