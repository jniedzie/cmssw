"""Attach the audited IR1 GDML and field in one model-frame contract.

Identity here means only that FLUKA model coordinates are retained for this
software closure.  It is not an approved IR1-to-IR5/CMS physics transform.
"""

import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from PhysicsTools.ShiftLssGeometry.shiftLssExternalGeometry_cff import (
    customiseShiftLssExternalGeometry,
)
from PhysicsTools.ShiftMuonSegments.shiftLssIr1AtlasProxy_cff import (
    shiftLssIr1AtlasProxyFieldElements,
)
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_customise import (
    customiseShiftLssMagneticField,
)


process = cms.Process("LSSIR1INTEGRATION")
process.load("Configuration.Geometry.GeometryDD4hepSimDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase1_2023_realistic", "")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

model_origin = (0.0, 0.0, 0.0)
model_to_test = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
process = customiseShiftLssExternalGeometry(
    process,
    gdmlFile=(
        "PhysicsTools/ShiftLssGeometry/data/ir1_atlas_proxy/"
        "lhc_ir1_atlas_proxy_bounded.gdml"
    ),
    artifactOriginInModelCm=(0.0, 4299.5, 14575.200000105498),
    modelOriginCm=model_origin,
    modelToCms=model_to_test,
    minimumAbsZCm=1100.0,
    checkOverlaps=True,
)
process = customiseShiftLssMagneticField(
    process,
    fieldElements=shiftLssIr1AtlasProxyFieldElements(
        modelOriginCm=model_origin,
        modelToCms=model_to_test,
        fieldScale=1.0,
    ),
)

process.verifyShiftLssGeometry = cms.EDAnalyzer("ShiftLssGeometryVerifier")
process.validateShiftLssField = cms.EDAnalyzer(
    "ShiftLssMagneticFieldValidator",
    samples=cms.VPSet(
        cms.PSet(
            name=cms.string("MQXA aperture in shared geometry-field frame"),
            pointCm=cms.vdouble(1.0, 0.0, 2615.0),
            expectedTesla=cms.vdouble(0.0, 1.8928521203, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MBXW kick in shared geometry-field frame"),
            pointCm=cms.vdouble(1.0, 0.0, 6132.2),
            expectedTesla=cms.vdouble(0.0, -1.1999831604, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("target coordinate outside assigned proxy field"),
            pointCm=cms.vdouble(-9.7, 0.0, 14800.0),
            expectedTesla=cms.vdouble(0.0, 0.0, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
    ),
)
process.validation = cms.Path(process.verifyShiftLssGeometry + process.validateShiftLssField)
