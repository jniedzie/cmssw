import FWCore.ParameterSet.Config as cms
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonSegments, shiftMuonSegmentsCounter
process = cms.Process("SHIFTSEGMENTS")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("file:events_AOD.root"))
process.p = cms.Path(shiftMuonSegmentsCounter + shiftMuonSegments)
