import FWCore.ParameterSet.Config as cms

shiftMuonSegments = cms.EDProducer(
    "ShiftMuonSegmentsTableProducer",
    dtSegments=cms.InputTag("dt4DSegments"),
    cscSegments=cms.InputTag("cscSegments"),
)
shiftMuonSegmentsCounter = cms.EDAnalyzer(
    "ShiftMuonSegmentsCounter",
    dtSegments=shiftMuonSegments.dtSegments,
    cscSegments=shiftMuonSegments.cscSegments,
)
