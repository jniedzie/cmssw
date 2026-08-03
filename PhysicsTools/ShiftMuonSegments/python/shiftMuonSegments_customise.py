from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cff import addShiftMuonSegments


def customise(process):
    return addShiftMuonSegments(process)


def customiseRecoDebug(process):
    """Run the reconstruction-funnel diagnostic after the RECO sequence."""
    from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonSegmentsCounter

    process.shiftMuonRecoDebug = shiftMuonSegmentsCounter.clone(printDetails=True)
    if hasattr(process, "reconstruction_step"):
        process.reconstruction_step += process.shiftMuonRecoDebug
    else:
        raise RuntimeError("Cannot attach Shift muon debug analyzer: no reconstruction_step")
    return process
