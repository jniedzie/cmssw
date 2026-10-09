"""Explicit collision-data feature contract; no simulation imports or fields."""
import math
import numpy as np

MUON_MASS = 0.105658
MUON_QUALITY = (
    "normalizedChi2", "nValidHits", "nValidMuonHits", "nMuonStations",
    "nLostHits", "topology", "recoAlgorithm", "ptErr", "etaErr", "phiErr",
)
MUON_KINEMATICS = ("pt", "pz", "phi", "vx", "vy", "vz")
PAIR_QUALITY = ("dca", "normalizedChi2", "probability", "vxErr", "vyErr", "vzErr")
PAIR_GEOMETRY = ("vx", "vy", "vz", "pt", "pz", "eta", "phi", "originCompatibilityNormalizedChi2")
QUALITY_NAMES = tuple("pair_" + k for k in PAIR_QUALITY) + tuple(
    f"mu_{side}_{k}" for k in MUON_QUALITY for side in ("min", "max")
)
FEATURE_NAMES = QUALITY_NAMES + tuple("pair_" + k for k in PAIR_GEOMETRY) + (
    "opening_angle", "mu_vertex_distance",
) + tuple(f"mu_{side}_{k}" for side in ("highPt", "lowPt")
          for k in ("px", "py", "pz", "energy", "vx", "vy", "vz"))
RECO_BRANCHES = tuple(dict.fromkeys(
    ("run", "luminosityBlock", "event", "nShiftMuon", "nShiftDimuonVertex")
    + tuple("ShiftMuon_" + k for k in MUON_QUALITY + MUON_KINEMATICS)
    + tuple("ShiftDimuonVertex_" + k for k in
            ("muonIdx1", "muonIdx2", "mass") + PAIR_QUALITY + PAIR_GEOMETRY)
))
PREDICTOR_BRANCHES = tuple(k for k in RECO_BRANCHES if k != "ShiftDimuonVertex_mass")


def signed_log(values):
    values = np.asarray(values)
    return np.sign(values) * np.log1p(np.abs(values))


def required_reco_branches():
    return RECO_BRANCHES


def extract_reco(arrays, *, include_nuisance=True):
    """Only access the allowlisted branches, including when extra branches exist.

    Each output row corresponds to an already retained ShiftDimuonVertex.
    Nonfinite feature values are preserved for training-only fitted imputation.
    Invalid pair references fail closed; no truth participates in row selection.
    """
    required = RECO_BRANCHES if include_nuisance else PREDICTOR_BRANCHES
    pair_fields = (("muonIdx1", "muonIdx2", "mass") if include_nuisance
                   else ("muonIdx1", "muonIdx2")) + PAIR_QUALITY + PAIR_GEOMETRY
    missing = set(required) - set(arrays)
    if missing:
        raise ValueError("Missing reconstructed branches: " + ", ".join(sorted(missing)))
    rows, ids, event_indices, pair_indices, muon_indices, nuisances = [], [], [], [], [], []
    for event_index in range(len(arrays["event"])):
        n_mu = int(arrays["nShiftMuon"][event_index])
        n_pair = int(arrays["nShiftDimuonVertex"][event_index])
        for k in MUON_QUALITY + MUON_KINEMATICS:
            if len(arrays["ShiftMuon_" + k][event_index]) != n_mu:
                raise ValueError("Muon collection lengths differ")
        for k in pair_fields:
            if len(arrays["ShiftDimuonVertex_" + k][event_index]) != n_pair:
                raise ValueError("Pair collection lengths differ")
        for pair_index in range(n_pair):
            def pair(k):
                return float(arrays["ShiftDimuonVertex_" + k][event_index][pair_index])
            first, second = int(pair("muonIdx1")), int(pair("muonIdx2"))
            if first == second or not (0 <= first < n_mu and 0 <= second < n_mu):
                raise ValueError("Invalid reconstructed pair reference")
            def mu(i, k):
                return float(arrays["ShiftMuon_" + k][event_index][i])
            ordered = sorted((first, second), key=lambda i: (-mu(i, "pt"), i))
            vecs = [np.array([mu(i, "pt") * (math.cos(mu(i, "phi")) if math.isfinite(mu(i, "phi")) else np.nan),
                              mu(i, "pt") * (math.sin(mu(i, "phi")) if math.isfinite(mu(i, "phi")) else np.nan), mu(i, "pz")])
                    for i in ordered]
            denominator = np.linalg.norm(vecs[0]) * np.linalg.norm(vecs[1])
            angle = (math.acos(float(np.clip(np.dot(*vecs) / denominator, -1, 1)))
                     if denominator > 0 and np.isfinite(denominator) else np.nan)
            values = [pair(k) for k in PAIR_QUALITY]
            for k in MUON_QUALITY:
                a, b = mu(first, k), mu(second, k)
                values.extend((np.minimum(a, b), np.maximum(a, b)))
            values.extend(pair(k) for k in PAIR_GEOMETRY)
            distance = np.linalg.norm([mu(first, k) - mu(second, k) for k in ("vx", "vy", "vz")])
            values.extend((angle, distance))
            for i, vector in zip(ordered, vecs):
                values.extend((*vector, math.sqrt(float(np.dot(vector, vector)) + MUON_MASS ** 2),
                               mu(i, "vx"), mu(i, "vy"), mu(i, "vz")))
            rows.append(values)
            ids.append([int(arrays[k][event_index]) for k in ("run", "luminosityBlock", "event")])
            event_indices.append(event_index)
            pair_indices.append(pair_index)
            muon_indices.append((first, second))
            if include_nuisance:
                nuisances.append((pair("mass"), pair("vz"), math.hypot(pair("vx"), pair("vy"))))
    result = dict(
        X=np.asarray(rows, dtype=float).reshape(-1, len(FEATURE_NAMES)),
        ids=np.asarray(ids, dtype=np.uint64).reshape(-1, 3),
        event_index=np.asarray(event_indices, dtype=int), pair_index=np.asarray(pair_indices, dtype=int),
        muon_indices=np.asarray(muon_indices, dtype=int).reshape(-1, 2),
        nuisance=(np.asarray(nuisances, dtype=float).reshape(-1, 3) if include_nuisance else None),
        feature_names=np.asarray(FEATURE_NAMES),
    )
    if include_nuisance and result["nuisance"].size and not np.isfinite(result["nuisance"]).all():
        raise ValueError("Nonfinite reconstructed mass or vertex nuisance")
    if include_nuisance and result["nuisance"].size and (result["nuisance"][:, 0] < 0).any():
        raise ValueError("Negative reconstructed mass nuisance")
    return result


def nuisance_transform(values):
    """Mass and vertex coordinates are adversary targets, never truth values."""
    values = np.asarray(values)
    return np.column_stack((np.log1p(values[:, 0]), signed_log(values[:, 1]),
                            np.log1p(values[:, 2])))
