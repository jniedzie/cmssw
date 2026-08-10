import FWCore.ParameterSet.Config as cms
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonSegments, shiftMuonSegmentsCounter, shiftMuonTable

shiftMuonSegmentsSequence = cms.Sequence(shiftMuonSegmentsCounter + shiftMuonTable + shiftMuonSegments)

def addShiftMuonSegments(process):
    process.load("RecoMuon.TransientTrackingRecHit.MuonTransientTrackingRecHitBuilder_cfi")
    process.shiftMuonSegmentsCounter = shiftMuonSegmentsCounter.clone()
    process.shiftMuonTable = shiftMuonTable.clone()
    # NanoAOD retains FlatTables with module labels matching ``*Table``.  Keep
    # that suffix here so the standard NanoAOD output commands write both
    # ShiftDT/CSC/GEM segment tables and ShiftRPC reconstructed-hit table.
    process.shiftMuonSegmentsTable = shiftMuonSegments.clone()
    shift_muon_modules = (
        process.shiftMuonSegmentsCounter + process.shiftMuonTable + process.shiftMuonSegmentsTable
    )
    if (
        process.shiftMuonTable.useImprovedMomentumRefit.value()
        and process.shiftMuonTable.useDetailedMaterialPropagation.value()
    ):
        # The optional detailed-material study needs the full DD4hep Geant4
        # geometry.  Do not pay its initialization or single-stream cost while
        # the study is disabled.
        process.load("Configuration.Geometry.GeometryDD4hepSimDB_cff")
        from SimG4Core.Application.g4SimHits_cfi import g4SimHits as _g4SimHits
        process.shiftMuonGeant4Geometry = cms.EDProducer(
            "GeometryProducer",
            GeoFromDD4hep=cms.bool(True),
            UseMagneticField=cms.bool(True),
            UseSensitiveDetectors=cms.bool(False),
            MagneticField=_g4SimHits.MagneticField.clone(),
        )
        process.options.numberOfThreads = cms.untracked.uint32(1)
        process.options.numberOfStreams = cms.untracked.uint32(1)
        shift_muon_modules = process.shiftMuonGeant4Geometry + shift_muon_modules
    process.shiftMuonSegmentsSequence = cms.Sequence(shift_muon_modules)

    # Be explicit as well as following the standard ``*Table`` convention:
    # this makes the tables resilient to a future restrictive NanoAOD output
    # command list.
    for output_name in ("NANOAODoutput", "NANOAODSIMoutput"):
        if hasattr(process, output_name):
            getattr(process, output_name).outputCommands.append(
                "keep nanoaodFlatTable_shiftMuonSegmentsTable_*_*"
            )
            getattr(process, output_name).outputCommands.append(
                "keep nanoaodFlatTable_shiftMuonTable_*_*"
            )
    # cmsDriver's NANO step is the path that is actually scheduled and
    # written by the NANOAODSIM output module.  Attach the producer to that
    # path, rather than only defining it in the process.  The nanoSequence
    # fallback is retained for standalone/custom NanoAOD configurations.
    if hasattr(process, "nanoAOD_step"):
        process.nanoAOD_step += process.shiftMuonSegmentsSequence
    elif hasattr(process, "nanoSequence"):
        process.nanoSequence += process.shiftMuonSegmentsSequence
    else:
        raise RuntimeError(
            "Cannot attach ShiftMuonSegments: process has neither "
            "nanoAOD_step nor nanoSequence"
        )
    # The standard NanoAOD event content currently keeps all FlatTables, but
    # retain these products explicitly as well.  This makes the contract
    # independent of future EXONanoAOD event-content changes.
    for output_name in ("NANOAODoutput", "NANOAODSIMoutput"):
        if hasattr(process, output_name):
            output = getattr(process, output_name)
            for command in (
                "keep nanoaodFlatTable_shiftMuonSegments_ShiftDT_*",
                "keep nanoaodFlatTable_shiftMuonSegments_ShiftCSC_*",
                "keep nanoaodFlatTable_shiftMuonSegments_ShiftRPC_*",
                "keep nanoaodFlatTable_shiftMuonSegments_ShiftGEM_*",
                "keep nanoaodFlatTable_shiftMuonTable__*",
            ):
                if command not in output.outputCommands:
                    output.outputCommands.append(command)
    return process
