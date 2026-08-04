import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonTable

options = VarParsing("analysis")
options.register("inputFile", "file:events_AOD.root", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "AOD input file")
options.parseArguments()

process = cms.Process("SHIFTMUON")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFile))
process.shiftMuonTable = shiftMuonTable.clone()
process.p = cms.Path(process.shiftMuonTable)
