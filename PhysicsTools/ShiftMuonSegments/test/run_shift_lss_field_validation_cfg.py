import FWCore.ParameterSet.Config as cms

from PhysicsTools.ShiftMuonSegments.shiftLssMagneticField_cfi import (
    shiftLssUniformFieldElement,
)
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_customise import (
    customiseShiftLssMagneticField,
)


process = cms.Process("LSSFIELDTEST")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.load("MagneticField.Engine.uniformMagneticField_cfi")
process.UniformMagneticFieldESProducer.ZFieldInTesla = cms.double(3.8)

elements = [
    shiftLssUniformFieldElement(
        "boundedDipole",
        (-10.0, -10.0, -100.0),
        (10.0, 10.0, 100.0),
        (0.0, 1.2, 0.0),
        boundsShape="cylinderZ",
        outerRadiusCm=10.0,
        excludedCylindersCm=(0.0, 0.0, 2.0),
    ),
]
process = customiseShiftLssMagneticField(
    process,
    baseMagneticFieldProducer="UniformMagneticFieldESProducer",
    fieldElements=elements,
)
process.validateShiftLssField = cms.EDAnalyzer(
    "ShiftLssMagneticFieldValidator",
    samples=cms.VPSet(
        cms.PSet(
            name=cms.string("inside element"),
            pointCm=cms.vdouble(5.0, 0.0, 0.0),
            expectedTesla=cms.vdouble(0.0, 1.2, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("excluded beam hole delegates to CMS"),
            pointCm=cms.vdouble(0.0, 0.0, 0.0),
            expectedTesla=cms.vdouble(0.0, 0.0, 3.8),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("outside element delegates to CMS"),
            pointCm=cms.vdouble(20.0, 0.0, 0.0),
            expectedTesla=cms.vdouble(0.0, 0.0, 3.8),
            toleranceTesla=cms.double(1.0e-6),
        ),
    ),
)
process.validation = cms.Path(process.validateShiftLssField)
