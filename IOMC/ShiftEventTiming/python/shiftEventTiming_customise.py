import FWCore.ParameterSet.Config as cms

from IOMC.ShiftEventTiming.shiftEventTime_cfi import shiftEventTime


def customiseShiftEventTiming(
    process,
    timingMode="nominal",
    beamDirectionZ=-1,
    bxOffset=0,
    phaseNs=0.0,
    fixedOffsetNs=0.0,
    cmsReferenceZmm=0.0,
    bunchSpacingNs=25.0,
    legacyOffsetCtMm=-148000.0,
    modelVersion="fixed-target-v1",
    maxTrackTimeNs=5000.0,
    maxTrackTimeForwardNs=5000.0,
):
    """Apply one coherent time shift after vertex smearing and before Geant4."""
    if not hasattr(process, "generatorSmeared") or not hasattr(process, "pgen"):
        raise RuntimeError(
            "customiseShiftEventTiming requires the standard GEN sequence "
            "with process.pgen and process.generatorSmeared"
        )
    if maxTrackTimeNs <= 0.0 or maxTrackTimeForwardNs <= 0.0:
        raise ValueError("Geant4 maximum track times must be positive")

    process.shiftEventTime = shiftEventTime.clone(
        timingMode=timingMode,
        beamDirectionZ=beamDirectionZ,
        bxOffset=bxOffset,
        phaseNs=phaseNs,
        fixedOffsetNs=fixedOffsetNs,
        cmsReferenceZmm=cmsReferenceZmm,
        bunchSpacingNs=bunchSpacingNs,
        legacyOffsetCtMm=legacyOffsetCtMm,
        modelVersion=modelVersion,
    )
    process.pgen.replace(
        process.generatorSmeared,
        process.shiftEventTime + process.generatorSmeared,
    )
    process.generatorSmeared.currentTag = cms.untracked.InputTag("shiftEventTime")

    # The standard 500 ns central transport cutoff is comparable to the
    # 148 m source-to-CMS flight time.  Without the old negative Pythia offset,
    # an inward primary reaches CMS at roughly that cutoff and is then killed
    # by SteppingAction.  Timing must reach detector response code first; this
    # broad Geant4 transport guard must not act as an accidental readout gate.
    if hasattr(process, "g4SimHits"):
        for actionName in ("Physics", "StackingAction", "SteppingAction"):
            action = getattr(process.g4SimHits, actionName)
            action.MaxTrackTime = cms.double(maxTrackTimeNs)
            action.MaxTrackTimeForward = cms.double(maxTrackTimeForwardNs)

    # FEVTDEBUG currently keeps these products, but make the provenance contract
    # explicit so a future event-content change cannot silently drop the timing.
    for output in process.outputModules_().values():
        if hasattr(output, "outputCommands"):
            output.outputCommands.append("keep *_shiftEventTime_*_*")
    return process
