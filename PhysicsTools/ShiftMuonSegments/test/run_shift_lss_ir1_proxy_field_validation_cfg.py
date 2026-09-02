"""Validate frozen IR1 field samples after a non-trivial proper transform.

This is a model-frame software regression.  It deliberately does not claim
that the chosen transform or the IR1 polarities describe CMS/IR5.
"""

import FWCore.ParameterSet.Config as cms

from PhysicsTools.ShiftMuonSegments.shiftLssIr1AtlasProxy_cff import (
    shiftLssIr1AtlasProxyFieldElements,
)
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_customise import (
    customiseShiftLssMagneticField,
)


process = cms.Process("LSSIR1FIELDTEST")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.load("MagneticField.Engine.uniformMagneticField_cfi")
process.UniformMagneticFieldESProducer.ZFieldInTesla = cms.double(3.8)

# A +90 degree rotation around z verifies that both sampling positions and
# field vectors obey global = modelOrigin + modelToCms * model.
model_origin = (100.0, 200.0, -100.0)
model_to_test = (0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
process = customiseShiftLssMagneticField(
    process,
    baseMagneticFieldProducer="UniformMagneticFieldESProducer",
    fieldElements=shiftLssIr1AtlasProxyFieldElements(
        modelOriginCm=model_origin,
        modelToCms=model_to_test,
        fieldScale=1.0,
    ),
)


def global_point(local_x, local_y, model_z):
    return (model_origin[0] - local_y, model_origin[1] + local_x, model_origin[2] + model_z)


process.validateShiftLssIr1Field = cms.EDAnalyzer(
    "ShiftLssMagneticFieldValidator",
    samples=cms.VPSet(
        cms.PSet(
            name=cms.string("MQXA aperture analytic quadrupole"),
            pointCm=cms.vdouble(*global_point(1.0, 0.0, 2615.0)),
            expectedTesla=cms.vdouble(-1.8928521203, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MQXB analytic quadrupole and signed assignment"),
            pointCm=cms.vdouble(*global_point(1.0, 0.0, 3480.0)),
            expectedTesla=cms.vdouble(1.8928521203, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MBXW analytic kick core"),
            pointCm=cms.vdouble(*global_point(1.0, 0.0, 6132.2)),
            expectedTesla=cms.vdouble(1.1999831604, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MQTL displaced quadrupole origin"),
            pointCm=cms.vdouble(*global_point(10.7, 0.0, 19649.0)),
            expectedTesla=cms.vdouble(0.37950611892, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-4),
        ),
        cms.PSet(
            name=cms.string("constant corrector vector rotation"),
            pointCm=cms.vdouble(*global_point(1.0, 0.0, 2984.2)),
            expectedTesla=cms.vdouble(0.0, 0.36866803721, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MQXA interpolated outer field"),
            pointCm=cms.vdouble(*global_point(5.0, 0.0, 2615.0)),
            expectedTesla=cms.vdouble(-3.7952574652511535, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-4),
        ),
        cms.PSet(
            name=cms.string("MQXB interpolated outer field"),
            pointCm=cms.vdouble(*global_point(0.0, 5.90909090909091, 3480.0)),
            expectedTesla=cms.vdouble(2.27142254436e-5, 0.8601309319855229, 0.0),
            toleranceTesla=cms.double(1.0e-4),
        ),
        cms.PSet(
            name=cms.string("MBXW interpolated field outside kick core"),
            pointCm=cms.vdouble(*global_point(5.0, 5.0, 6132.2)),
            expectedTesla=cms.vdouble(1.294253837481024, 0.1702488108649104, 0.0),
            toleranceTesla=cms.double(1.0e-4),
        ),
        cms.PSet(
            name=cms.string("MQTL interpolated displaced map"),
            pointCm=cms.vdouble(*global_point(12.191525423728812, 2.491525423728813, 19649.0)),
            expectedTesla=cms.vdouble(0.65627499786881, -0.65627499786881, 0.0),
            toleranceTesla=cms.double(1.0e-4),
        ),
        cms.PSet(
            name=cms.string("outside proxy delegates to CMS field"),
            pointCm=cms.vdouble(0.0, 0.0, 0.0),
            expectedTesla=cms.vdouble(0.0, 0.0, 3.8),
            toleranceTesla=cms.double(1.0e-6),
        ),
    ),
)
process.validation = cms.Path(process.validateShiftLssIr1Field)
