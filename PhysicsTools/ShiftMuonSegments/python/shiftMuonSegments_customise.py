import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cff import addShiftMuonSegments


def customise(
    process,
    directionalRefitUseDetailedMaterialEffects=None,
    directionalRefitUseGeometryMaterialEffects=None,
    directionalRefitUseGeometryMaterialEffectsInFitter=None,
    directionalRefitUseGeometryMaterialEffectsInSmoother=None,
    directionalRefitUseGeometryTargetMaterialEffects=None,
    enableHcalDiagnostics=False,
    enableZDCDiagnostics=False,
    augmentDTHits=False,
    augmentTrackerHits=False,
    useExtendedTiming=False,
):
    process = addShiftMuonSegments(
        process,
        directionalRefitUseDetailedMaterialEffects=directionalRefitUseDetailedMaterialEffects,
        directionalRefitUseGeometryMaterialEffects=directionalRefitUseGeometryMaterialEffects,
        directionalRefitUseGeometryMaterialEffectsInFitter=directionalRefitUseGeometryMaterialEffectsInFitter,
        directionalRefitUseGeometryMaterialEffectsInSmoother=directionalRefitUseGeometryMaterialEffectsInSmoother,
        directionalRefitUseGeometryTargetMaterialEffects=directionalRefitUseGeometryTargetMaterialEffects,
        enableHcalDiagnostics=enableHcalDiagnostics,
        enableZDCDiagnostics=enableZDCDiagnostics,
        augmentDTHits=augmentDTHits,
        augmentTrackerHits=augmentTrackerHits,
        useExtendedTiming=useExtendedTiming,
    )
    # Standard NanoAOD already provides every generator quantity used by the
    # SHIFT TEA analysis except pz and the production vertex.  Add only those
    # four columns here instead of enabling the much larger EXONanoAOD table
    # set.  Full float precision is important at the O(150 m) SHIFT baseline:
    # EXONanoAOD's eight-bit vertex precision quantises z by tens of cm.
    if hasattr(process, "genParticleTable"):
        process.genParticleTable.variables.pz = Var(
            "pz",
            float,
            doc="gen particle longitudinal momentum",
            precision=23,
        )
        for coordinate in ("vx", "vy", "vz"):
            setattr(
                process.genParticleTable.variables,
                coordinate,
                Var(
                    coordinate,
                    float,
                    doc=f"gen particle production vertex {coordinate[-1]} coordinate (cm)",
                    precision=23,
                ),
            )
    return process


def customiseKeepShiftTruth(
    process,
    keepHcalSimHits=False,
    keepZDCSimHits=False,
    zdcDigiTimePhaseOffset=0.0,
):
    """Retain generator and simulation truth needed for SHIFT hit attribution."""
    truth_commands = [
        "keep edmHepMCProduct_generator_*_*",
        "keep SimTracks_g4SimHits_*_*",
        "keep SimVertexs_g4SimHits_*_*",
        "keep PSimHits_g4SimHits_*_*",
        "keep *_mix_MergedTrackTruth_*",
        "keep TrackingParticles_mix_MergedTrackTruth_*",
        "keep TrackingVertexs_mix_MergedTrackTruth_*",
    ]
    # Direct tracker-hit association is performed in Step 4 without requiring
    # a collision-style tracker track. Preserve the rechits made in Step 3.
    truth_commands.extend((
        "keep *_siPixelRecHits_*_*",
        "keep *_siStripMatchedRecHits_*_*",
        "keep *_reducedEcalRecHitsEB_*_*",
        "keep *_reducedEcalRecHitsEE_*_*",
        "keep *_reducedHcalRecHits_*_*",
        "keep *_hbhereco_*_*",
        "keep *_hfreco_*_*",
        "keep *_horeco_*_*",
        "keep *_zdcreco_*_*",
    ))
    if keepHcalSimHits:
        truth_commands.append("keep PCaloHits_g4SimHits_HcalHits_*")
    if keepZDCSimHits:
        truth_commands.append("keep PCaloHits_g4SimHits_ZDCHITS_*")
        # Phase-I input contains legacy ZDCDataFramesSorted from
        # simHcalUnsuppressedDigis, not the QIE10 hcalDigis:ZDC product used
        # by the default Run-3 module.
        from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_zdc_cfi import zdcreco
        process.shiftZDCReco = zdcreco.clone(
            digiLabelhcal=cms.InputTag("simHcalUnsuppressedDigis"),
            digiLabelQIE10ZDC=cms.InputTag(""),
        )
        process.shiftZDCReco_step = cms.Path(process.shiftZDCReco)
        if hasattr(process, "schedule"):
            process.schedule.append(process.shiftZDCReco_step)
        else:
            raise RuntimeError("Cannot attach legacy SHIFT ZDC reconstruction: no process schedule")
        truth_commands.append("keep *_shiftZDCReco_*_*")
    # SHIFT muons cross the source-side ZDC roughly 0.93 microseconds earlier
    # than a collision particle emitted at the interaction point would reach
    # the same cell.  The ordinary ten-sample ZDC frame consequently contains
    # no signal.  A dedicated Step-2 stream can shift the digitizer phase while
    # leaving the ordinary campaign and all central-detector timing untouched.
    if zdcDigiTimePhaseOffset:
        adjusted = False
        if (
            hasattr(process, "mix")
            and hasattr(process.mix, "digitizers")
            and hasattr(process.mix.digitizers, "hcal")
            and hasattr(process.mix.digitizers.hcal, "zdc")
        ):
            hcal = process.mix.digitizers.hcal
            # Run-3 DIGI disables the legacy Phase-I ZDC digitizer by
            # default.  Changing its timing PSet alone therefore cannot make
            # a collection; explicitly instantiate it for this stream.
            if hasattr(hcal, "doZDCDigi"):
                hcal.doZDCDigi = cms.bool(True)
            zdc = hcal.zdc
            zdc.timePhase = cms.double(
                float(zdc.timePhase.value()) + float(zdcDigiTimePhaseOffset)
            )
            adjusted = True
        if not adjusted:
            for module_name in ("simHcalUnsuppressedDigis", "simHcalDigis"):
                if not hasattr(process, module_name):
                    continue
                module = getattr(process, module_name)
                if hasattr(module, "digiCfg") and hasattr(module.digiCfg, "zdc"):
                    module.digiCfg.zdc.timePhase = cms.double(
                        float(module.digiCfg.zdc.timePhase.value()) + float(zdcDigiTimePhaseOffset)
                    )
                    adjusted = True
        if not adjusted:
            raise RuntimeError("Requested SHIFT ZDC time-phase shift, but no HCAL digitizer was found")
    for output in process.outputModules_().values():
        if not hasattr(output, "outputCommands"):
            continue
        existing = set(output.outputCommands)
        for command in truth_commands:
            if command not in existing:
                output.outputCommands.append(command)
    return process


def customiseRecoDebug(
    process,
    enableDTMeasurement=True,
    enableGEMMeasurement=True,
    trackerMode="general",
    enableHcalDiagnostics=False,
    enableZDCDiagnostics=False,
    dtNavigationMode=1,
    recoVariantCode=0,
):
    """Run the reconstruction-funnel diagnostic on a dedicated end path."""
    from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import (
        shiftMuonRecoDiagnostics,
        shiftMuonSegmentsCounter,
    )

    tracker_mode_codes = {"none": 0, "general": 1, "p5": 2}
    if trackerMode not in tracker_mode_codes:
        raise ValueError(f"Unsupported SHIFT tracker mode: {trackerMode}")
    process.shiftMuonRecoDiagnostics = shiftMuonRecoDiagnostics.clone(
        enableDTMeasurement=enableDTMeasurement,
        enableGEMMeasurement=enableGEMMeasurement,
        trackerMode=tracker_mode_codes[trackerMode],
        enableHcalDiagnostics=enableHcalDiagnostics,
        enableZDCDiagnostics=enableZDCDiagnostics,
        dtNavigationMode=dtNavigationMode,
        recoVariantCode=recoVariantCode,
    )
    if enableZDCDiagnostics:
        process.shiftMuonRecoDiagnostics.zdcRecHits = cms.InputTag("shiftZDCReco")
    if trackerMode == "p5":
        process.shiftMuonRecoDiagnostics.generalTracks = cms.InputTag(
            "ctfWithMaterialTracksP5"
        )
    process.shiftMuonRecoDebug = shiftMuonSegmentsCounter.clone(printDetails=True)
    process.shiftMuonRecoDebug_step = cms.EndPath(
        process.shiftMuonRecoDiagnostics + process.shiftMuonRecoDebug
    )
    if hasattr(process, "schedule"):
        process.schedule.append(process.shiftMuonRecoDebug_step)
    else:
        raise RuntimeError("Cannot attach Shift muon debug analyzer: no process schedule")
    for output in process.outputModules_().values():
        if hasattr(output, "outputCommands"):
            output.outputCommands.append(
                "keep nanoaodFlatTable_shiftMuonRecoDiagnostics_*_*"
            )

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


def customiseTraversingShiftMuonReco(
    process,
    trackerMode="general",
    enableDTMeasurement=True,
):
    """Add parallel cosmic-style through-going muon and tracker reconstruction."""
    from RecoMuon.MuonSeedGenerator.CosmicMuonSeedProducer_cfi import CosmicMuonSeed
    from RecoMuon.CosmicMuonProducer.cosmicMuons_cfi import cosmicMuons
    from RecoMuon.CosmicMuonProducer.globalCosmicMuons_cfi import globalCosmicMuons

    process.shiftCosmicMuonSeed = CosmicMuonSeed.clone(
        ForcePointDown=False,
        TryBothDirections=True,
        KeepAllSegments=True,
        MaxCSCChi2=1000.0,
        MaxDTChi2=1000.0,
        # A one-segment seed needs a non-zero momentum to define a valid
        # trajectory state, but this is only an initial hypothesis, not a cut.
        SingleSegmentPt=0.01,
        # Forward SHIFT muons can have tiny pT while retaining large |pz|.
        MinPairPt=0.0,
        PairSegmentPt=0.01,
        UsePairPtEstimate=False,
        TryBothPairCharges=True,
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
            # Require hits in separated detector hemispheres.  With False the
            # cosmic builder accepts one-leg and mixed-muon combinations into
            # the collection that is supposed to represent traversing tracks.
            Strict1Leg=True,
        ),
    )
    for module in (process.shiftCosmicMuons, process.shiftTraversingMuons):
        module.TrajectoryBuilderParameters.EnableDTMeasurement = cms.bool(
            enableDTMeasurement
        )

    # Keep the detector-only collections canonical, and form independent
    # tracker+muon test hypotheses with CMSSW's established cosmic global fit.
    # The stock collision thresholds reject the very small transverse momenta
    # expected for SHIFT, so relax only those preselection thresholds; retain
    # the covariance-aware spatial matcher and all of its quality cuts.
    tracker_collections = {
        "general": "generalTracks",
        "p5": "ctfWithMaterialTracksP5",
    }
    if trackerMode not in ("none", *tracker_collections):
        raise ValueError(f"Unsupported SHIFT tracker mode: {trackerMode}")

    def make_global(muon_collection, tracker_collection):
        module = globalCosmicMuons.clone(
            MuonCollectionLabel=muon_collection,
            TrajectoryBuilderParameters=dict(
                TkTrackCollectionLabel=tracker_collection,
            ),
        )
        matcher = module.TrajectoryBuilderParameters.GlobalMuonTrackMatcher
        matcher.MinP = 0.01
        matcher.MinPt = 0.0
        return module

    shift_reco_sequence = (
        process.shiftCosmicMuonSeed
        + process.shiftCosmicMuons
        + process.shiftTraversingMuons
    )
    if trackerMode != "none":
        tracker_collection = tracker_collections[trackerMode]
        process.shiftGlobalDSAMuons = make_global(
            "displacedStandAloneMuons", tracker_collection
        )
        process.shiftGlobalCosmicMuons = make_global(
            "shiftCosmicMuons", tracker_collection
        )
        process.shiftGlobalTraversingMuons = make_global(
            "shiftTraversingMuons", tracker_collection
        )
        shift_reco_sequence += (
            process.shiftGlobalDSAMuons
            + process.shiftGlobalCosmicMuons
            + process.shiftGlobalTraversingMuons
        )

    if trackerMode == "p5":
        process.load("RecoTracker.Configuration.RecoTrackerP5_cff")
        # The stock P5 seeder assumes field-off cosmics in the outer TEC.
        # SHIFT is field-on and the observed signal crossing uses TID plus
        # inner TEC rings, so retain the bounded forward prototype explicitly.
        process.combinatorialcosmicseedfinderP5.requireBOFF = cms.bool(False)
        for layers in (
            process.combinatorialcosmicseedingpairsTECposP5,
            process.combinatorialcosmicseedingpairsTECnegP5,
        ):
            layers.TEC.minRing = cms.int32(1)
            layers.TEC.maxRing = cms.int32(7)
        process.shiftP5Tracker_step = cms.Path(process.ctftracksP5)
        if hasattr(process, "schedule"):
            process.schedule.append(process.shiftP5Tracker_step)
        else:
            raise RuntimeError("Cannot attach Shift P5 tracker reconstruction: no process schedule")

    process.shiftTraversingMuon_step = cms.Path(shift_reco_sequence)
    if hasattr(process, "schedule"):
        process.schedule.append(process.shiftTraversingMuon_step)
    else:
        raise RuntimeError("Cannot attach traversing Shift muon reconstruction: no process schedule")

    keep_commands = [
        "keep *_shiftCosmicMuonSeed_*_*",
        "keep *_shiftCosmicMuons_*_*",
        "keep *_shiftTraversingMuons_*_*",
    ]
    if trackerMode != "none":
        keep_commands.extend(
            (
                "keep *_shiftGlobalDSAMuons_*_*",
                "keep *_shiftGlobalCosmicMuons_*_*",
                "keep *_shiftGlobalTraversingMuons_*_*",
            )
        )
    if trackerMode == "general":
        keep_commands.append("keep *_generalTracks_*_*")
    elif trackerMode == "p5":
        keep_commands.extend(
            (
                "keep *_ctfWithMaterialTracksP5_*_*",
                "keep *_ctfWithMaterialTracksCosmics_*_*",
            )
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
    enableDTMeasurement=True,
    enableGEMMeasurement=True,
):
    """Tune explicitly selected DSA cuts, keeping each trial reproducible."""
    if not hasattr(process, "displacedStandAloneMuons"):
        raise RuntimeError("Cannot tune Shift DSA reconstruction: displacedStandAloneMuons is absent")
    builder = process.displacedStandAloneMuons.STATrajBuilderParameters
    builder.FilterParameters.NumberOfSigma = numberOfSigma
    builder.BWFilterParameters.NumberOfSigma = numberOfSigma
    builder.FilterParameters.MuonTrajectoryUpdatorParameters.MaxChi2 = maxHitChi2
    builder.BWFilterParameters.MuonTrajectoryUpdatorParameters.MaxChi2 = maxHitChi2
    builder.FilterParameters.EnableDTMeasurement = enableDTMeasurement
    builder.BWFilterParameters.EnableDTMeasurement = enableDTMeasurement
    builder.FilterParameters.EnableGEMMeasurement = enableGEMMeasurement
    builder.BWFilterParameters.EnableGEMMeasurement = enableGEMMeasurement
    builder.SeedPosition = seedPosition
    builder.DoBackwardFilter = doBackwardFilter
    builder.NavigationType = navigationType
    process.displacedMuonSeeds.KeepAllSegments = keepAllSeedSegments
    # The cosmic seed producer otherwise assigns 10 GeV to every one-segment
    # seed and rejects curvature estimates below 10 GeV.  Neither behavior is
    # appropriate for the very low-pT SHIFT signal.
    process.displacedMuonSeeds.SingleSegmentPt = 0.01
    process.displacedMuonSeeds.MinPairPt = 0.0
    process.displacedMuonSeeds.PairSegmentPt = 0.01
    process.displacedMuonSeeds.UsePairPtEstimate = False
    process.displacedMuonSeeds.TryBothPairCharges = True
    process.displacedMuonSeeds.TryBothDirections = True
    if hasattr(process, "displacedMuonReducedTrackExtras"):
        process.displacedMuonReducedTrackExtras.cut = cms.string("pt > 0")

    # SHIFT muons originate far outside CMS and need not extrapolate to the
    # beam line.  The default standalone loader rejects a valid trajectory
    # when propagation to the PCA fails outside the tracker; retain it using
    # its geometrically innermost state, as done for cosmic-muon tracks.
    process.displacedStandAloneMuons.TrackLoaderParameters.AllowNoVertex = cms.untracked.bool(True)
    process.displacedStandAloneMuons.TrackLoaderParameters.MuonUpdatorAtVertexParameters.Propagator = pcaPropagator
    return process
