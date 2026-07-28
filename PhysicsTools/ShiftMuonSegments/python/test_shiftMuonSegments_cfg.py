import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonSegments, shiftMuonSegmentsCounter

options = VarParsing("analysis")
options.register("inputFile", "file:events_AOD.root", VarParsing.multiplicity.singleton, VarParsing.varType.string, "AOD input file")
options.outputFile = "shiftMuonSegments_test.root"
options.parseArguments()

process = cms.Process("SHIFTSEGMENTS")
process.MessageLogger = cms.Service(
    "MessageLogger",
    cerr=cms.untracked.PSet(
        threshold=cms.untracked.string("INFO"),
        ShiftMuonSegments=cms.untracked.PSet(limit=cms.untracked.int32(100000)),
    ),
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFile))
process.shiftMuonSegmentsCounter = shiftMuonSegmentsCounter.clone()
process.shiftMuonSegments = shiftMuonSegments.clone()
process.p = cms.Path(process.shiftMuonSegmentsCounter + process.shiftMuonSegments)
process.out = cms.OutputModule("PoolOutputModule", fileName=cms.untracked.string(options.outputFile), SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring("p")))
process.end = cms.EndPath(process.out)
