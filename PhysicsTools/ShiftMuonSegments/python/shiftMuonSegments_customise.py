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
    augmentDTHits=True,
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
    keepMergedTrackTruth=True,
    keepSimMuonRPCDigis=False,
    keepPileupPlayback=False,
    zdcDigiTimePhaseOffset=0.0,
):
    """Retain generator and simulation truth needed for SHIFT hit attribution."""
    truth_commands = [
        "keep edmHepMCProduct_generator_*_*",
        "keep SimTracks_g4SimHits_*_*",
        "keep SimVertexs_g4SimHits_*_*",
        "keep PSimHits_g4SimHits_*_*",
    ]
    if keepMergedTrackTruth:
        # Retain this only for workflows that explicitly run RECOSIM. Pileup
        # makes these collections much larger than the reconstructed event;
        # canonical SHIFT reconstruction and Step 4 use g4SimHits directly.
        truth_commands.append("keep *_mix_MergedTrackTruth_*")
    else:
        # AODSIM event content keeps this collection by default. Output
        # commands are applied in order, so an explicit final drop is needed.
        truth_commands.append("drop *_mix_MergedTrackTruth_*")
    if keepSimMuonRPCDigis:
        # Run-3 RAW2DIGI configures muonRPCDigis as an RPCDigiMerger whose
        # simulation input is not packed into rawDataCollector. GENRAW keeps
        # only its DigiSimLinks, so retain the digi collection needed by the
        # downstream merger explicitly.
        truth_commands.append("keep *_simMuonRPCDigis_*_*")
    if keepPileupPlayback:
        # Preserve the exact mixed-event identities for occupancy provenance
        # without retaining the much larger mixed TrackingParticle graph.
        truth_commands.append("keep CrossingFramePlaybackInfoNew_mix_*_*")
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


def customiseRecoDiagnostics(
    process,
    enableDTMeasurement=True,
    enableGEMMeasurement=True,
    trackerMode="general",
    enableHcalDiagnostics=False,
    enableZDCDiagnostics=False,
    dtNavigationMode=1,
):
    """Persist compact reconstruction provenance on a dedicated end path."""
    from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import (
        shiftMuonRecoDiagnostics,
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
    )
    if enableZDCDiagnostics:
        process.shiftMuonRecoDiagnostics.zdcRecHits = cms.InputTag("shiftZDCReco")
    if trackerMode == "p5":
        process.shiftMuonRecoDiagnostics.generalTracks = cms.InputTag(
            "ctfWithMaterialTracksP5LHCNavigation"
        )
    process.shiftMuonRecoDiagnostics_step = cms.EndPath(process.shiftMuonRecoDiagnostics)
    if hasattr(process, "schedule"):
        process.schedule.append(process.shiftMuonRecoDiagnostics_step)
    else:
        raise RuntimeError("Cannot attach Shift muon debug analyzer: no process schedule")
    for output in process.outputModules_().values():
        if hasattr(output, "outputCommands"):
            output.outputCommands.append(
                "keep nanoaodFlatTable_shiftMuonRecoDiagnostics_*_*"
            )

    return process


# Backward-compatible alias for archived generated configurations.
customiseRecoDebug = customiseRecoDiagnostics


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
        # The P5 cosmic final selector imposes prompt/cosmic kinematic cuts
        # that remove the recovered forward tracks. Use the independently
        # fitted LHC-navigation collection as the tracker-seeded hypothesis;
        # it still enters ShiftMuon only through a successful tracker+muon
        # global match below.
        "p5": "ctfWithMaterialTracksP5LHCNavigation",
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
        # The stock two-hit P5 seed fixes the *total* momentum to 5 GeV and
        # then immediately applies a 0.5 GeV pT CKF cut.  At the observed
        # |eta|~6-7 this constructs pT~0.01 GeV and rejects every seed before
        # the first tracker measurement.  Use a forward-appropriate initial
        # momentum, retain a deliberately broad covariance, admit both
        # charges, and let hit compatibility rather than that inconsistent
        # prompt/cosmic pT threshold decide whether a candidate grows.
        process.combinatorialcosmicseedfinderP5.SeedMomentum = cms.double(500.0)
        process.combinatorialcosmicseedfinderP5.ErrorRescaling = cms.double(1000.0)
        process.combinatorialcosmicseedfinderP5.Charges = cms.vint32(-1, 1)
        process.ckfBaseTrajectoryFilterP5.minPt = cms.double(0.0)
        # The forward seed already contains two TEC measurements. Require one
        # genuine extension hit for this prototype and avoid the stock P5
        # high-|eta| five-hit rule rejecting the trajectory while it grows.
        # Candidate purity is evaluated later with independent detector and
        # truth diagnostics; two-hit seed-only trajectories are never kept.
        process.ckfBaseTrajectoryFilterP5.minimumNumberOfHits = cms.int32(3)
        process.ckfBaseTrajectoryFilterP5.minHitsAtHighEta = cms.int32(3)
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
                "keep *_combinatorialcosmicseedfinderP5_*_*",
                "keep *_ckfTrackCandidatesP5_*_*",
                "keep *_ckfTrackCandidatesP5LHCNavigation_*_*",
                "keep *_ctfWithMaterialTracksCosmics_*_*",
                "keep *_ctfWithMaterialTracksP5_*_*",
                "keep *_ctfWithMaterialTracksP5LHCNavigation_*_*",
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


def customiseShiftLssTransport(
    process,
    *,
    magneticFieldLabel="",
    baseMagneticFieldLabel="shiftLssBaseMagneticField",
    baseMagneticFieldProducer="VolumeBasedMagneticFieldESProducer",
    fieldElements=None,
    addBaseField=False,
    sumOverlaps=False,
    materialBoundaryAbsZCm=1100.0,
    geant4eMomentumLimitGeV=0.05,
    geant4eMaximumStepLengthMm=2.0,
    geant4eMaximumPathLengthCm=2500.0,
):
    """Configure simulation and reconstruction to share an optional LSS field."""
    if not hasattr(process, "shiftMuonTable"):
        raise RuntimeError("Cannot configure SHIFT LSS transport: shiftMuonTable is absent")
    positive_values = {
        "materialBoundaryAbsZCm": materialBoundaryAbsZCm,
        "geant4eMomentumLimitGeV": geant4eMomentumLimitGeV,
        "geant4eMaximumStepLengthMm": geant4eMaximumStepLengthMm,
        "geant4eMaximumPathLengthCm": geant4eMaximumPathLengthCm,
    }
    invalid = [name for name, value in positive_values.items() if not value > 0.0]
    if invalid:
        raise ValueError("SHIFT LSS transport values must be positive: " + ", ".join(invalid))
    if fieldElements is not None:
        if magneticFieldLabel:
            raise ValueError(
                "The LSS composite field must use the empty label required by g4SimHits"
            )
        if magneticFieldLabel == baseMagneticFieldLabel:
            raise ValueError("LSS composite and base magnetic-field labels must differ")
        if not hasattr(process, baseMagneticFieldProducer):
            raise RuntimeError(
                "Cannot configure SHIFT LSS field: base EventSetup producer "
                + baseMagneticFieldProducer
                + " is absent"
            )
        from PhysicsTools.ShiftMuonSegments.shiftLssMagneticField_cfi import shiftLssMagneticField

        getattr(process, baseMagneticFieldProducer).label = cms.untracked.string(baseMagneticFieldLabel)
        process.shiftLssMagneticField = shiftLssMagneticField.clone(
            baseFieldLabel=baseMagneticFieldLabel,
            outputLabel=magneticFieldLabel,
            addBaseField=addBaseField,
            sumOverlaps=sumOverlaps,
            elements=cms.VPSet(*fieldElements),
        )
        # g4SimHits consumes the empty-label product. The standard CMS field
        # has been moved to baseMagneticFieldLabel, leaving the composite as
        # the single default product consumed below by reconstruction too.
    process.shiftMuonTable.lssTransport = cms.PSet(
        magneticFieldLabel=cms.string(magneticFieldLabel),
        materialBoundaryAbsZCm=cms.double(materialBoundaryAbsZCm),
        geant4eMomentumLimitGeV=cms.double(geant4eMomentumLimitGeV),
        geant4eMaximumStepLengthMm=cms.double(geant4eMaximumStepLengthMm),
        geant4eMaximumPathLengthCm=cms.double(geant4eMaximumPathLengthCm),
    )
    return process
