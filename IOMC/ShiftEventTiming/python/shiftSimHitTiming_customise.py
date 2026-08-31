import FWCore.ParameterSet.Config as cms

from IOMC.ShiftEventTiming.shiftSimHitTime_cfi import shiftSimHitTime


_MUON_SIMHIT_INSTANCES = (
    "MuonDTHits",
    "MuonCSCHits",
    "MuonRPCHits",
    "MuonGEMHits",
)


def customiseShiftSimHitReferenceTiming(
    process,
    bxOffset=0,
    phaseNs=0.0,
    bunchSpacingNs=25.0,
    modelVersion="same-simhit-reference-v1",
):
    """Shift fixed muon PSimHit times before unchanged no-pileup digitization."""
    if not hasattr(process, "mix") or not hasattr(process, "pdigiTask"):
        raise RuntimeError(
            "customiseShiftSimHitReferenceTiming requires the standard DIGI sequence"
        )
    if hasattr(process.mix, "input"):
        raise RuntimeError(
            "same-SimHit reference timing is currently restricted to no-pileup controls"
        )

    process.shiftSimHitTime = shiftSimHitTime.clone(
        bxOffset=bxOffset,
        phaseNs=phaseNs,
        bunchSpacingNs=bunchSpacingNs,
        modelVersion=modelVersion,
    )
    process.pdigiTask.add(process.shiftSimHitTime)

    shifted_inputs = []
    replaced = set()
    for tag in process.mix.mixObjects.mixSH.input:
        instance = tag.productInstanceLabel
        if tag.moduleLabel == "g4SimHits" and instance in _MUON_SIMHIT_INSTANCES:
            shifted_inputs.append(cms.InputTag("shiftSimHitTime", instance))
            replaced.add(instance)
        else:
            shifted_inputs.append(tag)
    missing = set(_MUON_SIMHIT_INSTANCES) - replaced
    if missing:
        raise RuntimeError(
            "standard mix configuration lacks muon SimHit inputs: "
            + ", ".join(sorted(missing))
        )
    process.mix.mixObjects.mixSH.input = cms.VInputTag(*shifted_inputs)

    collection_names = {
        "simMuonDTDigis": ("InputCollection", "InputCollectionPU", "MuonDTHits"),
        "simMuonCSCDigis": ("InputCollection", "InputCollectionPU", "MuonCSCHits"),
        "simMuonRPCDigis": ("InputCollection", "InputCollectionPU", "MuonRPCHits"),
        "simMuonGEMDigis": ("inputCollection", "inputCollectionPU", "MuonGEMHits"),
    }
    for module_name, (signal_name, pileup_name, instance) in collection_names.items():
        if not hasattr(process, module_name):
            raise RuntimeError(f"standard DIGI sequence lacks {module_name}")
        module = getattr(process, module_name)
        shifted_name = f"shiftSimHitTime{instance}"
        setattr(module, signal_name, cms.string(shifted_name))
        setattr(module, pileup_name, cms.string(shifted_name))

    for output in process.outputModules_().values():
        if hasattr(output, "outputCommands"):
            output.outputCommands.append("keep *_shiftSimHitTime_*_*")
    return process
