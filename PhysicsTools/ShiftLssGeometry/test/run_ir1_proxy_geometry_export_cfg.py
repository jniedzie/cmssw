"""Export the combined CMS plus IR1-proxy geometry for visualization.

The identity placement is an IR1 model-frame software view, not an approved
IR1-to-IR5/CMS transform. Set SHIFT_LSS_GEOMETRY_EXPORT_FILE to the desired
TGeo ROOT output path.
"""

import os

import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from PhysicsTools.ShiftLssGeometry.shiftLssExternalGeometry_cff import (
    customiseShiftLssExternalGeometry,
)


output_file = os.environ.get("SHIFT_LSS_GEOMETRY_EXPORT_FILE")
if not output_file:
    raise RuntimeError("SHIFT_LSS_GEOMETRY_EXPORT_FILE is required")

process = cms.Process("LSSGEOMETRYEXPORT")
process.load("Configuration.Geometry.GeometryDD4hepSimDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase1_2023_realistic", "")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process = customiseShiftLssExternalGeometry(
    process,
    gdmlFile=(
        "PhysicsTools/ShiftLssGeometry/data/ir1_atlas_proxy/"
        "lhc_ir1_atlas_proxy_bounded.gdml"
    ),
    artifactOriginInModelCm=(0.0, 4299.5, 14575.200000105498),
    modelOriginCm=(0.0, 0.0, 0.0),
    modelToCms=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    minimumAbsZCm=1100.0,
    checkOverlaps=False,
)

process.exportShiftLssGeometry = cms.EDAnalyzer(
    "ShiftLssGeometryVerifier",
    outputFile=cms.string(output_file),
)
process.exportGeometry = cms.Path(process.exportShiftLssGeometry)
