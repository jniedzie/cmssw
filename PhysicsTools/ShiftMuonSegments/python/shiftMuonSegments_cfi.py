import FWCore.ParameterSet.Config as cms

# Commit marker for the workflow-wide CMSSW conditions consistency test.
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
    dsaGlobalLinks=cms.InputTag("shiftGlobalDSAMuons"),
    cosmicGlobalLinks=cms.InputTag("shiftGlobalCosmicMuons"),
    traversingGlobalLinks=cms.InputTag("shiftGlobalTraversingMuons"),
    genParticles=cms.InputTag("finalGenParticles"),
    simTracks=cms.InputTag("g4SimHits"),
    simVertices=cms.InputTag("g4SimHits"),
    dtSimHits=cms.InputTag("g4SimHits", "MuonDTHits"),
    cscSimHits=cms.InputTag("g4SimHits", "MuonCSCHits"),
    rpcSimHits=cms.InputTag("g4SimHits", "MuonRPCHits"),
    gemSimHits=cms.InputTag("g4SimHits", "MuonGEMHits"),
    muonRecHitBuilder=cms.string("MuonRecHitBuilder"),
    # Optional scale-improvement study: propagated-path hit ordering and the
    # guarded precision-only refit. Keep it enabled for the current iteration.
    useImprovedMomentumRefit=cms.bool(True),
    # Use explicitly backward Geant4e transport to recover material losses
    # between the source-facing fitted state and the CMS material boundary.
    useDetailedMaterialPropagation=cms.bool(False),
    # Keep material updates in the canonical Kalman fit. The material-free
    # v9 control remains available through this switch for focused ablations.
    directionalRefitUseMaterialEffects=cms.bool(True),
    # Canonical material transport: query the Geant4 detector material at
    # every curved SteppingHelix integration point and use the same model in
    # the filter, smoother, and backward target extrapolation.
    # The exact Geant4-material provider is retained as a controlled closure
    # option.  A paired reconstruction test found it neutral relative to the
    # standard SteppingHelix material model, while costing substantially more
    # CPU, so production keeps the standard model by default.
    directionalRefitUseFirstPrinciplesMaterialEffects=cms.bool(False),
    directionalRefitFirstPrinciplesStepCm=cms.double(0.2),
    # Optional first-principles material ablation. When enabled, the Kalman
    # filter uses detailed Geant4e transport along the incoming muon momentum;
    # CMSSW's smoother constructs the matching opposite-to-momentum clone for
    # its backward pass. Keep the approximate SteppingHelix model canonical
    # until the detailed mode passes the small-sample closure test.
    directionalRefitUseDetailedMaterialEffects=cms.bool(False),
    # Geometry-sampled mean-loss experiment: retain SteppingHelix covariance
    # and scattering, disable only its approximate continuous dE/dx, and use
    # Geant4 material-specific stopping power sampled along each propagation.
    directionalRefitUseGeometryMaterialEffects=cms.bool(False),
    # Split-pass ablation switches. The umbrella option above still enables
    # geometry mean loss in both passes for backwards compatibility.
    directionalRefitUseGeometryMaterialEffectsInFitter=cms.bool(False),
    directionalRefitUseGeometryMaterialEffectsInSmoother=cms.bool(False),
    # Focused closure test: sample the DD4hep/Geant4 material only between the
    # source-facing precision state and the outer CMS material boundary. This
    # replaces the standard R-Z map where it is known to return zero at the
    # forward muon radius; it does not scale dE/dx or alter the hit fit.
    directionalRefitUseGeometryTargetMaterialEffects=cms.bool(False),
    directionalRefitGeometryMaterialStepCm=cms.double(1.0),
    directionalRefitLogGeometryMaterialComparison=cms.bool(False),
    # Isolate hit ordering from the other refit changes. Use the established
    # source-side signed-z ordering for this comparison while retaining the
    # precision refit, both iterations, and target-only Geant4e correction.
    usePropagatedPathOrdering=cms.bool(False),
    # V7 experiment: restore the established refit's full five-parameter seed
    # covariance inflation while keeping the smoother uncertainty at 10.
    directionalRefitUseFullSeedErrorRescale=cms.bool(True),
    # Scale either the full seed covariance above or only q/p in the retained
    # curvature-only implementation.
    directionalRefitSeedCurvatureErrorRescale=cms.double(100.0),
    # Pass two starts from pass one's posterior state. Vary only this inflation
    # to test whether reusing the same-hit posterior as a prior causes the
    # observed downward q/p feedback; 100 preserves the established result.
    directionalRefitSecondSeedErrorRescale=cms.double(100.0),
    # Diagnostic nonlinear-convergence probe: change only the magnitude of
    # the first-iteration momentum seed, retaining its direction and charge.
    directionalRefitSeedMomentumScale=cms.double(1.0),
    # Diagnostic only: scale the continuous SteppingHelix energy-loss term in
    # the Kalman fit while leaving scattering and target propagation intact.
    directionalRefitEnergyLossScale=cms.double(1.0),
    # Backward-start uncertainty used internally by KFTrajectorySmoother.
    directionalRefitErrorRescale=cms.double(10.0),
    # These now drive a real remove-and-refit outlier loop. The conventional
    # EstimateCut=20 probe worsened the 100-event closure, so production keeps
    # rejection effectively disabled pending a larger, detector-aware study.
    directionalRefitInitialMaxHitChi2=cms.double(100000.0),
    directionalRefitMaxHitChi2=cms.double(100000.0),
    # Reject a second pass that moves too far from the first smoothed q/p.
    directionalRefitMaxRelativeQoverPChange=cms.double(0.5),
    # Muon DT/CSC RecHits cannot improve with a track hypothesis, so a second
    # pass does not relinearize them. Reusing pass one's same-hit posterior as
    # a prior double-counts information and is retained only as a diagnostic.
    directionalRefitUseSecondIteration=cms.bool(False),
    # Use the precision-only result canonically only when it spans at least
    # two independent DT, CSC, or GEM station layers.
    directionalRefitMinPrecisionStations=cms.uint32(2),
    # Reject a precision-only solution that is incompatible with the all-hit
    # curvature; its diagnostic branches are still retained.
    directionalRefitMaxPrecisionRelativeQoverPChange=cms.double(0.5),
    # Simulation-only closure columns and source/opposite-side precision
    # refits. Disabled in production because the extra fits cost CPU.
    produceMomentumClosureDiagnostics=cms.bool(True),
    produceSplitLegRefits=cms.bool(False),
    # Backward target transport is now the canonical state contract.  This
    # compatibility parameter is retained for older configurations; it no
    # longer permits anyDirection on the source-facing-to-target leg.
    directionalRefitUseExplicitBackwardTargetPropagation=cms.bool(False),
    # Preserve the detector-only result and additionally form a prompt-target
    # hypothesis by updating the independently smoothed detector state with an
    # x/y measurement on the configured target plane.  The 50 cm z width
    # matches Beams:sigmaVertexZ in the current production and is projected
    # into x/y by the fitted track slopes.
    produceTargetConstrainedMomentum=cms.bool(True),
    targetUseInferredSide=cms.bool(True),
    targetX=cms.double(0.0),
    targetY=cms.double(0.0),
    targetZ=cms.double(14800.0),
    targetSigmaX=cms.double(0.1),
    targetSigmaY=cms.double(0.1),
    targetSigmaZ=cms.double(50.0),
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
    # MC-only, read-only transport closure probe. Keep disabled in production.
    printPropagationClosure=cms.bool(False),
)

# Persist an event-level reconstruction funnel while Step 3 still has access
# to RAW2DIGI products. A value of -1 means that a stage/collection was absent,
# which is distinct from a present but empty collection.
shiftMuonRecoDiagnostics = cms.EDProducer(
    "ShiftMuonRecoDiagnosticsProducer",
    dtDigis=cms.InputTag("muonDTDigis"),
    cscStripDigis=cms.InputTag("muonCSCDigis", "MuonCSCStripDigi"),
    cscWireDigis=cms.InputTag("muonCSCDigis", "MuonCSCWireDigi"),
    rpcDigis=cms.InputTag("muonRPCDigis"),
    gemDigis=cms.InputTag("muonGEMDigis"),
    dtRecHits=cms.InputTag("dt1DRecHits"),
    dtSegments=cms.InputTag("dt4DSegments"),
    cscRecHits=cms.InputTag("csc2DRecHits"),
    cscSegments=cms.InputTag("cscSegments"),
    rpcRecHits=cms.InputTag("rpcRecHits"),
    gemRecHits=cms.InputTag("gemRecHits"),
    gemSegments=cms.InputTag("gemSegments"),
    dtSimHits=cms.InputTag("g4SimHits", "MuonDTHits"),
    cscSimHits=cms.InputTag("g4SimHits", "MuonCSCHits"),
    rpcSimHits=cms.InputTag("g4SimHits", "MuonRPCHits"),
    gemSimHits=cms.InputTag("g4SimHits", "MuonGEMHits"),
    generalTracks=cms.InputTag("generalTracks"),
    dsaGlobalLinks=cms.InputTag("shiftGlobalDSAMuons"),
    cosmicGlobalLinks=cms.InputTag("shiftGlobalCosmicMuons"),
    traversingGlobalLinks=cms.InputTag("shiftGlobalTraversingMuons"),
)
