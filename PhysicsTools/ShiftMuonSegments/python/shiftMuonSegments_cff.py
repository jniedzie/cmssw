import FWCore.ParameterSet.Config as cms
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonSegments, shiftMuonSegmentsCounter

shiftMuonSegmentsSequence = cms.Sequence(shiftMuonSegmentsCounter + shiftMuonSegments)

def addShiftMuonSegments(process):
    process.shiftMuonSegmentsCounter = shiftMuonSegmentsCounter.clone()
    process.shiftMuonSegments = shiftMuonSegments.clone()
    process.shiftMuonSegmentsSequence = cms.Sequence(process.shiftMuonSegmentsCounter + process.shiftMuonSegments)
    # cmsDriver's NANO step is the path that is actually scheduled and
    # written by the NANOAODSIM output module.  Attach the producer to that
    # path, rather than only defining it in the process.  The nanoSequence
    # fallback is retained for standalone/custom NanoAOD configurations.
    if hasattr(process, "nanoAOD_step"):
        process.nanoAOD_step += process.shiftMuonSegmentsSequence
    elif hasattr(process, "nanoSequence"):
        process.nanoSequence += process.shiftMuonSegmentsSequence
    else:
        raise RuntimeError(
            "Cannot attach ShiftMuonSegments: process has neither "
            "nanoAOD_step nor nanoSequence"
        )
    return process
