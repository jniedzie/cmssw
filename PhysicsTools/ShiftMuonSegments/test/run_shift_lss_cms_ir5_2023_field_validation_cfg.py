"""Validate installed CMS IR5 field maps and exact domains under a proper rotation.

This is a model-frame software test. It does not approve the absolute CMS
alignment, source-side convention, or native FLUKA field reference.
"""

import FWCore.ParameterSet.Config as cms

from PhysicsTools.ShiftMuonSegments.shiftLssCmsIr5_2023Z1100_cff import (
    _ELEMENTS,
    shiftLssCmsIr5_2023Z1100FieldElements,
)
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_customise import (
    customiseShiftLssMagneticField,
)

assert len(_ELEMENTS) == 22
process = cms.Process("CMSIR5LSSFIELDTEST")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.load("MagneticField.Engine.uniformMagneticField_cfi")
process.UniformMagneticFieldESProducer.ZFieldInTesla = cms.double(3.8)

model_origin = (100.0, 200.0, -100.0)
model_to_test = (0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
process = customiseShiftLssMagneticField(
    process,
    baseMagneticFieldProducer="UniformMagneticFieldESProducer",
    fieldElements=shiftLssCmsIr5_2023Z1100FieldElements(
        modelOriginCm=model_origin,
        modelToCms=model_to_test,
        fieldScale=1.0,
    ),
)


def global_point(local_x, local_y, model_z):
    return (model_origin[0] - local_y, model_origin[1] + local_x, model_origin[2] + model_z)


process.validateCmsIr5LssField = cms.EDAnalyzer(
    "ShiftLssMagneticFieldValidator",
    samples=cms.VPSet(
        cms.PSet(
            name=cms.string("MQXA analytic core and signed assignment"),
            pointCm=cms.vdouble(*global_point(1.0, 0.0, 2615.0)),
            expectedTesla=cms.vdouble(-1.9741973648239945, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MQXB analytic core and signed assignment"),
            pointCm=cms.vdouble(*global_point(1.0, 0.0, 3480.0)),
            expectedTesla=cms.vdouble(1.991642307086264, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MBXW analytic kick core"),
            pointCm=cms.vdouble(*global_point(1.0, 0.0, 6132.2)),
            expectedTesla=cms.vdouble(1.2553669967398966, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MQY displaced analytic origin"),
            pointCm=cms.vdouble(*global_point(10.7, 0.0, 16955.3)),
            expectedTesla=cms.vdouble(0.6439763556693612, 0.0, 0.0),
            # GlobalPoint stores centimetre coordinates at float precision; the
            # translated 210.7 cm point shifts local x by about 3e-6 cm.
            toleranceTesla=cms.double(5.0e-6),
        ),
        cms.PSet(
            name=cms.string("MCBXH constant corrector vector rotation"),
            pointCm=cms.vdouble(*global_point(1.0, 0.0, 2984.2)),
            expectedTesla=cms.vdouble(-0.3123278592383179, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MBRC constant dipole vector rotation"),
            pointCm=cms.vdouble(*global_point(1.0, 0.0, 15790.0)),
            expectedTesla=cms.vdouble(-2.7099984571548212, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("outside all LSS elements delegates to CMS"),
            pointCm=cms.vdouble(0.0, 0.0, 0.0),
            expectedTesla=cms.vdouble(0.0, 0.0, 3.8),
            toleranceTesla=cms.double(1.0e-6),
        ),
    ),
)
process.validation = cms.Path(process.validateCmsIr5LssField)
