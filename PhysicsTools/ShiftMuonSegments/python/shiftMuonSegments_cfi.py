import FWCore.ParameterSet.Config as cms

shiftMuonSegments = cms.EDProducer(
    "ShiftMuonSegmentsTableProducer",
    dtSegments=cms.InputTag("dt4DSegments"),
    cscSegments=cms.InputTag("cscSegments"),
    rpcRecHits=cms.InputTag("rpcRecHits"),
    gemSegments=cms.InputTag("gemSegments"),
)
shiftMuonTable = cms.EDProducer(
    "ShiftMuonTableProducer",
    dsaTracks=cms.InputTag("displacedStandAloneMuons"),
    cosmicTracks=cms.InputTag("shiftCosmicMuons"),
    traversingTracks=cms.InputTag("shiftTraversingMuons"),
    genParticles=cms.InputTag("finalGenParticles"),
    targetLinePropagator=cms.string("SteppingHelixPropagatorAny"),
    minSharedHitFraction=cms.double(0.5),
    minSharedDetIds=cms.uint32(2),
    maxDuplicateAngle=cms.double(0.03),
    maxDuplicateLineDistance=cms.double(30.0),
    # Physical transverse pointing compatibility with the unbounded beam/
    # target line.  This deliberately imposes no requirement on origin z.
    maxTargetLineDca=cms.double(200.0),
    # SHIFT particles enter along the beam/target line.  Keep the requirement
    # symmetric in eta so the same reconstruction supports either target side.
    minAbsEta=cms.double(3.0),
    originTransverseResolution=cms.double(100.0),
    originZResolution=cms.double(2000.0),
    commonVertexLineResolution=cms.double(100.0),
    commonVertexBeamLineResolution=cms.double(100.0),
    maxPairOriginNormalizedChi2=cms.double(9.0),
    maxPairDca=cms.double(500.0),
    maxDimuonVertices=cms.uint32(1),
    requireOppositeSign=cms.bool(True),
    maxGenDeltaR=cms.double(0.5),
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
