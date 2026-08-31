import FWCore.ParameterSet.Config as cms


shiftSimHitTime = cms.EDProducer(
    "ShiftSimHitTimeProducer",
    dtSimHits=cms.InputTag("g4SimHits", "MuonDTHits"),
    cscSimHits=cms.InputTag("g4SimHits", "MuonCSCHits"),
    rpcSimHits=cms.InputTag("g4SimHits", "MuonRPCHits"),
    gemSimHits=cms.InputTag("g4SimHits", "MuonGEMHits"),
    bxOffset=cms.int32(0),
    bunchSpacingNs=cms.double(25.0),
    phaseNs=cms.double(0.0),
    modelVersion=cms.string("same-simhit-reference-v1"),
)
