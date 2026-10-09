"""MC-only training labels. Never import this module from scoring code."""
import numpy as np

TRUTH_BRANCHES = (
    "ShiftMuon_hitGenPartIdx", "ShiftMuon_hitTruthPurity", "ShiftMuon_hitTruthMatchedLayers",
    "ShiftMuon_hitSimTrackId", "GenPart_pdgId", "GenPart_status", "GenPart_genPartIdxMother",
    "GenPart_vx", "GenPart_vy", "GenPart_vz",
)


def pair_label(event, first, second, common_tolerance_cm=0.01, separate_tolerance_cm=0.1):
    """Label distinct muons at a common *position*, independent of sample name.

    Unknown matches and the position-tolerance gap are not negative examples.
    Mother identity is insufficient: prompt muons from different mothers can
    share a physical vertex; muon copies in DY can have different mothers.
    """
    n_gen = len(event["GenPart_pdgId"])
    gen = [int(event["ShiftMuon_hitGenPartIdx"][i]) for i in (first, second)]
    for i, g in zip((first, second), gen):
        if not 0 <= g < n_gen or int(event["GenPart_status"][g]) != 1 or abs(int(event["GenPart_pdgId"][g])) != 13:
            return -1, "unknown_match"
        purity = float(event["ShiftMuon_hitTruthPurity"][i])
        if not np.isfinite(purity) or purity < 0.75 or int(event["ShiftMuon_hitTruthMatchedLayers"][i]) < 3:
            return -1, "unknown_match"
    same_track = (int(event["ShiftMuon_hitSimTrackId"][first]) >= 0 and
                  int(event["ShiftMuon_hitSimTrackId"][first]) == int(event["ShiftMuon_hitSimTrackId"][second]))
    if gen[0] == gen[1] or same_track:
        return 0, "same_muon_twice"
    positions = [np.array([float(event["GenPart_" + k][g]) for k in ("vx", "vy", "vz")]) for g in gen]
    if not all(np.isfinite(p).all() for p in positions):
        return -1, "unknown_vertex"
    distance = float(np.linalg.norm(positions[0] - positions[1]))
    if distance <= common_tolerance_cm:
        return 1, "common_vertex"
    if distance >= separate_tolerance_cm:
        return 0, "different_vertices"
    return -1, "vertex_tolerance_gap"


def extract_labels(arrays, reco):
    missing = set(TRUTH_BRANCHES) - set(arrays)
    if missing:
        raise ValueError("MC labels require: " + ", ".join(sorted(missing)))
    result, reasons = [], []
    for event_index, (first, second) in zip(reco["event_index"], reco["muon_indices"]):
        event = {k: arrays[k][event_index] for k in TRUTH_BRANCHES}
        label, reason = pair_label(event, int(first), int(second))
        result.append(label)
        reasons.append(reason)
    return np.asarray(result, dtype=np.int8), np.asarray(reasons, dtype="U32")


def label_diagnostics(arrays, reco):
    """Auxiliary label audits only; never predictors or adversary targets."""
    n = len(reco["event_index"])
    out = dict(gen_index=np.full((n, 2), -1, dtype=int),
               sim_track_id=np.full((n, 2), -1, dtype=int),
               match_purity=np.full((n, 2), np.nan), matched_layers=np.zeros((n, 2), dtype=int),
               match_state=np.full((n, 2), "invalid_index", dtype="U24"),
               vertex_distance_cm=np.full(n, np.nan), same_gen_muon=np.zeros(n, dtype=bool),
               same_sim_track=np.zeros(n, dtype=bool), same_immediate_mother=np.full(n, -1, dtype=np.int8))
    for row, (e, muons) in enumerate(zip(reco["event_index"], reco["muon_indices"])):
        event = {k: arrays[k][e] for k in TRUTH_BRANCHES}
        n_gen = len(event["GenPart_pdgId"])
        positions = []
        for side, i in enumerate(muons):
            g = int(event["ShiftMuon_hitGenPartIdx"][i])
            purity = float(event["ShiftMuon_hitTruthPurity"][i])
            layers = int(event["ShiftMuon_hitTruthMatchedLayers"][i])
            out["gen_index"][row, side] = g
            out["sim_track_id"][row, side] = int(event["ShiftMuon_hitSimTrackId"][i])
            out["match_purity"][row, side] = purity
            out["matched_layers"][row, side] = layers
            state = ("invalid_index" if not 0 <= g < n_gen else
                     "non_muon" if abs(int(event["GenPart_pdgId"][g])) != 13 else
                     "non_status1" if int(event["GenPart_status"][g]) != 1 else
                     "nonfinite_purity" if not np.isfinite(purity) else
                     "low_purity" if purity < 0.75 else "too_few_layers" if layers < 3 else "ok")
            out["match_state"][row, side] = state
            if 0 <= g < n_gen and abs(int(event["GenPart_pdgId"][g])) == 13 and int(event["GenPart_status"][g]) == 1:
                positions.append(np.array([float(event["GenPart_" + k][g]) for k in ("vx", "vy", "vz")]))
        g1, g2 = out["gen_index"][row]
        s1, s2 = out["sim_track_id"][row]
        out["same_gen_muon"][row] = g1 >= 0 and g1 == g2
        out["same_sim_track"][row] = s1 >= 0 and s1 == s2
        if len(positions) == 2 and all(np.isfinite(p).all() for p in positions):
            out["vertex_distance_cm"][row] = np.linalg.norm(positions[0] - positions[1])
            m1, m2 = (int(event["GenPart_genPartIdxMother"][g]) for g in (g1, g2))
            if 0 <= m1 < n_gen and 0 <= m2 < n_gen:
                out["same_immediate_mother"][row] = int(m1 == m2)
    return out
