#!/usr/bin/env python3

import argparse
from collections import defaultdict
import json
import math

from DataFormats.FWLite import Events, Handle


PREFIX = "shiftMuonTruthTrail"
INT_FIELDS = ("TrackId", "ParentId", "PdgId", "Step", "Kind")
FLOAT_FIELDS = (
    "X", "Y", "Z", "Px", "Py", "Pz", "GlobalTime", "TrackLength",
    "StepLength", "Density", "RadiationLength", "NuclearInteractionLength",
    "CumulativeX0", "CumulativeInteractionLengths", "CumulativeEnergyLoss",
    "EnergyDeposit",
)
STRING_FIELDS = ("Volume", "Material", "Process")


def read_product(event, cpp_type, suffix, module="g4SimHits", prefix=PREFIX):
    handle = Handle(cpp_type)
    event.getByLabel(module, prefix + suffix, handle)
    if not handle.isValid():
        raise RuntimeError(f"missing {module}:{prefix}{suffix}")
    if cpp_type == "vector<string>":
        return [str(value) for value in handle.product()]
    return list(handle.product()) if cpp_type.startswith("vector") else handle.product()[0]


def audit(path):
    report = {"input": path, "events": [], "valid": True}
    total_checkpoints = 0
    total_tracks = 0
    for event in Events(path):
        event_id = event.eventAuxiliary().event()
        values = {}
        for field in INT_FIELDS:
            values[field] = read_product(event, "vector<int>", field)
        for field in FLOAT_FIELDS:
            values[field] = read_product(event, "vector<float>", field)
        for field in ("FieldX", "FieldY", "FieldZ"):
            values[field] = read_product(
                event,
                "vector<float>",
                field[5:],
                module="shiftMuonTruthTrailField",
                prefix="field",
            )
        for field in STRING_FIELDS:
            values[field] = read_product(event, "vector<string>", field)
        truncated = read_product(event, "int", "Truncated")
        sizes = {field: len(entries) for field, entries in values.items()}
        unique_sizes = set(sizes.values())
        if len(unique_sizes) != 1:
            raise RuntimeError(f"event {event_id} has inconsistent truth-trail vector sizes: {sizes}")
        count = unique_sizes.pop()
        if any(
            not math.isfinite(value)
            for field in FLOAT_FIELDS + ("FieldX", "FieldY", "FieldZ")
            for value in values[field]
        ):
            raise RuntimeError(f"event {event_id} has a non-finite truth-trail value")

        tracks = defaultdict(list)
        for index, track_id in enumerate(values["TrackId"]):
            tracks[track_id].append(index)
        event_report = {
            "event": event_id,
            "checkpoints": count,
            "truncated": int(truncated),
            "tracks": [],
        }
        if truncated:
            report["valid"] = False
        for track_id, indices in sorted(tracks.items()):
            lengths = [values["TrackLength"][index] for index in indices]
            x0 = [values["CumulativeX0"][index] for index in indices]
            interaction_lengths = [
                values["CumulativeInteractionLengths"][index] for index in indices
            ]
            if any(second + 1e-5 < first for first, second in zip(lengths, lengths[1:])):
                raise RuntimeError(f"event {event_id} track {track_id} has decreasing path length")
            if any(second + 1e-6 < first for first, second in zip(x0, x0[1:])):
                raise RuntimeError(f"event {event_id} track {track_id} has decreasing radiation length")
            if any(
                second + 1e-6 < first
                for first, second in zip(interaction_lengths, interaction_lengths[1:])
            ):
                raise RuntimeError(
                    f"event {event_id} track {track_id} has decreasing interaction length"
                )
            boundary_pre = sum(values["Kind"][index] == 2 for index in indices)
            boundary_post = sum(values["Kind"][index] == 3 for index in indices)
            if boundary_pre != boundary_post:
                raise RuntimeError(
                    f"event {event_id} track {track_id} has unmatched boundary checkpoints"
                )
            first = indices[0]
            last = indices[-1]
            momentum = lambda index: math.sqrt(
                values["Px"][index] ** 2
                + values["Py"][index] ** 2
                + values["Pz"][index] ** 2
            )
            materials = []
            seen_materials = set()
            for index in indices:
                material = values["Material"][index]
                if material not in seen_materials:
                    materials.append(material)
                    seen_materials.add(material)
            event_report["tracks"].append({
                "track_id": track_id,
                "pdg_id": values["PdgId"][first],
                "checkpoints": len(indices),
                "boundary_crossings": boundary_pre,
                "start": {
                    "z_cm": values["Z"][first],
                    "p_GeV": momentum(first),
                    "material": values["Material"][first],
                },
                "end": {
                    "z_cm": values["Z"][last],
                    "p_GeV": momentum(last),
                    "material": values["Material"][last],
                },
                "path_cm": lengths[-1],
                "cumulative_x0": x0[-1],
                "cumulative_interaction_lengths": interaction_lengths[-1],
                "cumulative_energy_loss_GeV": values["CumulativeEnergyLoss"][last],
                "max_field_T": max(
                    math.sqrt(
                        values["FieldX"][index] ** 2
                        + values["FieldY"][index] ** 2
                        + values["FieldZ"][index] ** 2
                    )
                    for index in indices
                ),
                "materials": materials,
            })
        total_checkpoints += count
        total_tracks += len(tracks)
        report["events"].append(event_report)
    report["event_count"] = len(report["events"])
    report["track_count"] = total_tracks
    report["checkpoint_count"] = total_checkpoints
    return report


def main():
    parser = argparse.ArgumentParser(description="Audit a SHIFT Geant4 primary-muon truth trail")
    parser.add_argument("input")
    parser.add_argument("--output")
    args = parser.parse_args()
    report = audit(args.input)
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as output:
            output.write(text + "\n")
    print(text)
    raise SystemExit(0 if report["valid"] else 1)


if __name__ == "__main__":
    main()
