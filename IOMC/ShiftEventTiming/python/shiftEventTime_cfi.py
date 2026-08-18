import FWCore.ParameterSet.Config as cms


shiftEventTime = cms.EDProducer(
    "ShiftEventTimeProducer",
    src=cms.InputTag("VtxSmeared"),
    timingMode=cms.string("nominal"),
    beamDirectionZ=cms.int32(-1),
    bxOffset=cms.int32(0),
    cmsReferenceZmm=cms.double(0.0),
    bunchSpacingNs=cms.double(25.0),
    phaseNs=cms.double(0.0),
    fixedOffsetNs=cms.double(0.0),
    legacyOffsetCtMm=cms.double(-148000.0),
    modelVersion=cms.string("fixed-target-v1"),
)
