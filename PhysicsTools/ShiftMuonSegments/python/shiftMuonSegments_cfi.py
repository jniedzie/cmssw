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
    printDetails=cms.bool(True),
)
