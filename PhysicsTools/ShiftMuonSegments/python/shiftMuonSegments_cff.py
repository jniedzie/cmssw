import FWCore.ParameterSet.Config as cms
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonSegments, shiftMuonSegmentsCounter

shiftMuonSegmentsSequence = cms.Sequence(shiftMuonSegmentsCounter + shiftMuonSegments)

def addShiftMuonSegments(process):
    process.shiftMuonSegmentsCounter = shiftMuonSegmentsCounter.clone()
    process.shiftMuonSegments = shiftMuonSegments.clone()
    process.shiftMuonSegmentsSequence = cms.Sequence(process.shiftMuonSegmentsCounter + process.shiftMuonSegments)
    if hasattr(process, "nanoAOD_step"):
        process.nanoAOD_step += process.shiftMuonSegmentsSequence
    elif hasattr(process, "nanoSequence"):
        process.nanoSequence += process.shiftMuonSegmentsSequence
    return process
