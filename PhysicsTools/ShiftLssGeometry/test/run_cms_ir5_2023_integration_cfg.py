"""Attach the ownership-partitioned CMS IR5 LSS geometry and field together.

Identity retains the delivered model frame only for software closure. This test
does not approve the absolute model-to-CMS survey transform or source side.
"""

import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from PhysicsTools.ShiftLssGeometry.shiftLssExternalGeometry_cff import (
    customiseShiftLssExternalGeometry,
)
from PhysicsTools.ShiftMuonSegments.shiftLssCmsIr5_2023Z1100_cff import (
    shiftLssCmsIr5_2023Z1100FieldElements,
)
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_customise import (
    customiseShiftLssMagneticField,
)

process = cms.Process("CMSIR5LSSINTEGRATION")
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
        "PhysicsTools/ShiftLssGeometry/data/cms_ir5_2023_z1100/"
        "lhc_ir5_2023_physical_z1100.gdml"
    ),
    artifactOriginInModelCm=(-500.0, 2550.0, 9550.5),
    modelOriginCm=model_origin,
    modelToCms=model_to_test,
    minimumAbsZCm=1100.0,
    overlapToleranceCm=0.0001,
    checkOverlaps=True,
)
process = customiseShiftLssMagneticField(
    process,
    fieldElements=shiftLssCmsIr5_2023Z1100FieldElements(
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
            name=cms.string("MQXA shared geometry-field model frame"),
            pointCm=cms.vdouble(1.0, 0.0, 2615.0),
            expectedTesla=cms.vdouble(0.0, 1.9741973648239945, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
        cms.PSet(
            name=cms.string("MBXW shared geometry-field model frame"),
            pointCm=cms.vdouble(1.0, 0.0, 6132.2),
            expectedTesla=cms.vdouble(0.0, -1.2553669967398966, 0.0),
            toleranceTesla=cms.double(1.0e-6),
        ),
    ),
)
process.validation = cms.Path(process.verifyShiftLssGeometry + process.validateShiftLssField)
