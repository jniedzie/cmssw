import FWCore.ParameterSet.Config as cms

shiftMuonSegments = cms.EDProducer(
    "ShiftMuonSegmentsTableProducer",
    dtSegments=cms.InputTag("dt4DSegments"),
    cscSegments=cms.InputTag("cscSegments"),
    rpcRecHits=cms.InputTag("rpcRecHits"),
    gemSegments=cms.InputTag("gemSegments"),
)
shiftMuonSegmentsCounter = cms.EDAnalyzer(
    "ShiftMuonSegmentsCounter",
    dtRecHits=cms.InputTag("dt1DRecHits"),
    dtSegments=shiftMuonSegments.dtSegments,
    cscRecHits=cms.InputTag("csc2DRecHits"),
    cscSegments=shiftMuonSegments.cscSegments,
    rpcRecHits=shiftMuonSegments.rpcRecHits,
    gemRecHits=cms.InputTag("gemRecHits"),
    gemSegments=shiftMuonSegments.gemSegments,
    dsaSeeds=cms.InputTag("displacedMuonSeeds"),
    dsaTracks=cms.InputTag("displacedStandAloneMuons"),
    cosmicTracks=cms.InputTag("shiftCosmicMuons"),
    traversingTracks=cms.InputTag("shiftTraversingMuons"),
    cosmicTrackerTracks=cms.InputTag("generalTracks"),
    simTracks=cms.InputTag("g4SimHits"),
    simVertices=cms.InputTag("g4SimHits"),
    dtSimHits=cms.InputTag("g4SimHits", "MuonDTHits"),
    cscSimHits=cms.InputTag("g4SimHits", "MuonCSCHits"),
    rpcSimHits=cms.InputTag("g4SimHits", "MuonRPCHits"),
    gemSimHits=cms.InputTag("g4SimHits", "MuonGEMHits"),
    printDetails=cms.bool(True),
)
