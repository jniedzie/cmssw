import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from PhysicsTools.ShiftLssGeometry.shiftLssExternalGeometry_cff import (
    customiseShiftLssExternalGeometry,
)


process = cms.Process("LSSGEOMREJECT")
process.load("Configuration.Geometry.GeometryDD4hepSimDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase1_2023_realistic", "")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

customiseShiftLssExternalGeometry(
    process,
    gdmlFile="PhysicsTools/ShiftLssGeometry/test/fixtures/multi_union_external.gdml",
    modelOriginCm=(0.0, 0.0, 3000.0),
    modelToCms=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    minimumAbsZCm=2750.0,
    checkOverlaps=False,
)

process.verifyShiftLssGeometry = cms.EDAnalyzer("ShiftLssGeometryVerifier")
process.path = cms.Path(process.verifyShiftLssGeometry)
