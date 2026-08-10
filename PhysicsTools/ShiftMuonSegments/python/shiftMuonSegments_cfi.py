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
    muonRecHitBuilder=cms.string("MuonRecHitBuilder"),
    # Optional scale-improvement study: detailed Geant4e material transport,
    # propagated-path hit ordering, robust staged hit rejection, and the
    # guarded precision-only refit.  The latest validation worsened the
    # single-muon momentum scale, so retain the implementation but keep the
    # established directional refit as the production default.
    useImprovedMomentumRefit=cms.bool(False),
    # Only loosen the q/p seed. Angular and position errors retain the useful
    # original fit constraints while each nonlinear pass can re-estimate scale.
    directionalRefitSeedCurvatureErrorRescale=cms.double(100.0),
    # Backward-start uncertainty used internally by KFTrajectorySmoother.
    directionalRefitErrorRescale=cms.double(10.0),
    # Use a loose first pass to recover the unusual incoming trajectory, then
    # a substantially tighter second pass to suppress incompatible hits.
    directionalRefitInitialMaxHitChi2=cms.double(1000.0),
    directionalRefitMaxHitChi2=cms.double(100.0),
    # Reject a second pass that moves too far from the first smoothed q/p.
    directionalRefitMaxRelativeQoverPChange=cms.double(0.25),
    # Use the precision-only result canonically only when it spans at least
    # two independent DT, CSC, or GEM station layers.
    directionalRefitMinPrecisionStations=cms.uint32(2),
    # Reject a precision-only solution that is incompatible with the all-hit
    # curvature; its diagnostic branches are still retained.
    directionalRefitMaxPrecisionRelativeQoverPChange=cms.double(0.25),
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
