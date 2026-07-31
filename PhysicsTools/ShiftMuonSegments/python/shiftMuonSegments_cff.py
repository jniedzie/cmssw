import FWCore.ParameterSet.Config as cms
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonSegments, shiftMuonSegmentsCounter

shiftMuonSegmentsSequence = cms.Sequence(shiftMuonSegmentsCounter + shiftMuonSegments)

def addShiftMuonSegments(process):
    process.shiftMuonSegmentsCounter = shiftMuonSegmentsCounter.clone()
    # NanoAOD retains FlatTables with module labels matching ``*Table``.  Keep
    # that suffix here so the standard NanoAOD output commands write both
    # ShiftDT and ShiftCSC tables.
    process.shiftMuonSegmentsTable = shiftMuonSegments.clone()
    process.shiftMuonSegmentsSequence = cms.Sequence(process.shiftMuonSegmentsCounter + process.shiftMuonSegmentsTable)

    # Be explicit as well as following the standard ``*Table`` convention:
    # this makes the tables resilient to a future restrictive NanoAOD output
    # command list.
    for output_name in ("NANOAODoutput", "NANOAODSIMoutput"):
        if hasattr(process, output_name):
            getattr(process, output_name).outputCommands.append(
                "keep nanoaodFlatTable_shiftMuonSegmentsTable_*_*"
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
            ):
                if command not in output.outputCommands:
                    output.outputCommands.append(command)
    return process
