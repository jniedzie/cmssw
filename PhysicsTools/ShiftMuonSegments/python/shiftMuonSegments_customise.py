import FWCore.ParameterSet.Config as cms

from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cff import addShiftMuonSegments


def customise(process):
    return addShiftMuonSegments(process)


def customiseRecoDebug(process):
    """Run the reconstruction-funnel diagnostic on a dedicated end path."""
    from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonSegmentsCounter

    process.shiftMuonRecoDebug = shiftMuonSegmentsCounter.clone(printDetails=True)
    process.shiftMuonRecoDebug_step = cms.EndPath(process.shiftMuonRecoDebug)
    if hasattr(process, "schedule"):
        process.schedule.append(process.shiftMuonRecoDebug_step)
    else:
        raise RuntimeError("Cannot attach Shift muon debug analyzer: no process schedule")

    # LogVerbatim categories must be explicitly admitted by MessageLogger.
    # The default INFO threshold alone does not make a custom verbatim
    # category visible.
    if hasattr(process, "MessageLogger"):
        if hasattr(process.MessageLogger, "cerr"):
            process.MessageLogger.cerr.ShiftMuonRecoDebug = cms.untracked.PSet(
                limit=cms.untracked.int32(-1)
            )
    # Keep the path/module report in these intentionally small diagnostic
    # jobs; it proves that the analyzer ran even if logger settings regress.
    process.options.wantSummary = cms.untracked.bool(True)
    return process


def customiseRecoForShiftMuons(process, numberOfSigma=5.0, maxHitChi2=100.0):
    """Tune explicitly selected DSA cuts, keeping each trial reproducible."""
    if not hasattr(process, "displacedStandAloneMuons"):
        raise RuntimeError("Cannot tune Shift DSA reconstruction: displacedStandAloneMuons is absent")
    builder = process.displacedStandAloneMuons.STATrajBuilderParameters
    builder.FilterParameters.NumberOfSigma = numberOfSigma
    builder.BWFilterParameters.NumberOfSigma = numberOfSigma
    builder.FilterParameters.MuonTrajectoryUpdatorParameters.MaxChi2 = maxHitChi2
    builder.BWFilterParameters.MuonTrajectoryUpdatorParameters.MaxChi2 = maxHitChi2
    return process
