import FWCore.ParameterSet.Config as cms

from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cff import addShiftMuonSegments


def customise(process):
    return addShiftMuonSegments(process)


def customiseKeepShiftTruth(process):
    """Retain generator and simulation truth needed for SHIFT hit attribution."""
    truth_commands = (
        "keep edmHepMCProduct_generator_*_*",
        "keep SimTracks_g4SimHits_*_*",
        "keep SimVertexs_g4SimHits_*_*",
        "keep PSimHits_g4SimHits_*_*",
        "keep *_mix_MergedTrackTruth_*",
        "keep TrackingParticles_mix_MergedTrackTruth_*",
        "keep TrackingVertexs_mix_MergedTrackTruth_*",
    )
    for output in process.outputModules_().values():
        if not hasattr(output, "outputCommands"):
            continue
        existing = set(output.outputCommands)
        for command in truth_commands:
            if command not in existing:
                output.outputCommands.append(command)
    return process


def customiseRecoDebug(process):
    """Run the reconstruction-funnel diagnostic on a dedicated end path."""
    from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonSegmentsCounter

    process.shiftMuonRecoDebug = shiftMuonSegmentsCounter.clone(printDetails=True)
    process.shiftMuonRecoDebug_step = cms.EndPath(process.shiftMuonRecoDebug)
    if hasattr(process, "schedule"):
        process.schedule.append(process.shiftMuonRecoDebug_step)
    else:
        raise RuntimeError("Cannot attach Shift muon debug analyzer: no process schedule")

    # LogVerbatim categories must be explicitly admitted by MessageLogger.
    # The default INFO threshold alone does not make a custom verbatim
    # category visible.
    if hasattr(process, "MessageLogger"):
        if hasattr(process.MessageLogger, "cerr"):
            process.MessageLogger.cerr.ShiftMuonRecoDebug = cms.untracked.PSet(
                limit=cms.untracked.int32(-1)
            )
    # Keep the path/module report in these intentionally small diagnostic
    # jobs; it proves that the analyzer ran even if logger settings regress.
    process.options.wantSummary = cms.untracked.bool(True)
    return process


def customiseTraversingShiftMuonReco(process):
    """Add parallel cosmic-style through-going muon and tracker reconstruction."""
    from RecoMuon.MuonSeedGenerator.CosmicMuonSeedProducer_cfi import CosmicMuonSeed
    from RecoMuon.CosmicMuonProducer.cosmicMuons_cfi import cosmicMuons

    process.shiftCosmicMuonSeed = CosmicMuonSeed.clone(
        ForcePointDown=False,
        KeepAllSegments=True,
        MaxCSCChi2=1000.0,
        MaxDTChi2=1000.0,
        # A one-segment seed needs a non-zero momentum to define a valid
        # trajectory state, but this is only an initial hypothesis, not a cut.
        SingleSegmentPt=0.01,
        # Forward SHIFT muons can have tiny pT while retaining large |pz|.
        MinPairPt=0.0,
    )
    process.shiftCosmicMuons = cosmicMuons.clone(
        MuonSeedCollectionLabel="shiftCosmicMuonSeed",
        TrajectoryBuilderParameters=dict(
            BuildTraversingMuon=False,
            Strict1Leg=False,
        ),
    )
    process.shiftTraversingMuons = cosmicMuons.clone(
        MuonSeedCollectionLabel="shiftCosmicMuonSeed",
        TrajectoryBuilderParameters=dict(
            BuildTraversingMuon=True,
            Strict1Leg=False,
        ),
    )
    process.shiftTraversingMuon_step = cms.Path(
        process.shiftCosmicMuonSeed
        + process.shiftCosmicMuons
        + process.shiftTraversingMuons
    )
    if hasattr(process, "schedule"):
        process.schedule.append(process.shiftTraversingMuon_step)
    else:
        raise RuntimeError("Cannot attach traversing Shift muon reconstruction: no process schedule")

    keep_commands = (
        "keep *_shiftCosmicMuonSeed_*_*",
        "keep *_shiftCosmicMuons_*_*",
        "keep *_shiftTraversingMuons_*_*",
        "keep *_generalTracks_*_*",
    )
    for output in process.outputModules_().values():
        if hasattr(output, "outputCommands"):
            output.outputCommands.extend(keep_commands)
    return process


def customiseRecoForShiftMuons(
    process,
    numberOfSigma=5.0,
    maxHitChi2=100.0,
    seedPosition="in",
    doBackwardFilter=True,
    keepAllSeedSegments=True,
    navigationType="Standard",
    pcaPropagator="SteppingHelixPropagatorAny",
):
    """Tune explicitly selected DSA cuts, keeping each trial reproducible."""
    if not hasattr(process, "displacedStandAloneMuons"):
        raise RuntimeError("Cannot tune Shift DSA reconstruction: displacedStandAloneMuons is absent")
    builder = process.displacedStandAloneMuons.STATrajBuilderParameters
    builder.FilterParameters.NumberOfSigma = numberOfSigma
    builder.BWFilterParameters.NumberOfSigma = numberOfSigma
    builder.FilterParameters.MuonTrajectoryUpdatorParameters.MaxChi2 = maxHitChi2
    builder.BWFilterParameters.MuonTrajectoryUpdatorParameters.MaxChi2 = maxHitChi2
    builder.FilterParameters.EnableGEMMeasurement = True
    builder.BWFilterParameters.EnableGEMMeasurement = True
    builder.SeedPosition = seedPosition
    builder.DoBackwardFilter = doBackwardFilter
    builder.NavigationType = navigationType
    process.displacedMuonSeeds.KeepAllSegments = keepAllSeedSegments
    # The cosmic seed producer otherwise assigns 10 GeV to every one-segment
    # seed and rejects curvature estimates below 10 GeV.  Neither behavior is
    # appropriate for the very low-pT SHIFT signal.
    process.displacedMuonSeeds.SingleSegmentPt = 0.01
    process.displacedMuonSeeds.MinPairPt = 0.0

    # SHIFT muons originate far outside CMS and need not extrapolate to the
    # beam line.  The default standalone loader rejects a valid trajectory
    # when propagation to the PCA fails outside the tracker; retain it using
    # its geometrically innermost state, as done for cosmic-muon tracks.
    process.displacedStandAloneMuons.TrackLoaderParameters.AllowNoVertex = cms.untracked.bool(True)
    process.displacedStandAloneMuons.TrackLoaderParameters.MuonUpdatorAtVertexParameters.Propagator = pcaPropagator
    return process
