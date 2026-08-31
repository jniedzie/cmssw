#!/usr/bin/env python3

import unittest

import FWCore.ParameterSet.Config as cms

from PhysicsTools.ShiftLssGeometry.shiftLssExternalGeometry_cff import (
    customiseShiftLssExternalGeometry,
)


class ShiftLssExternalGeometryTest(unittest.TestCase):
    def process(self):
        process = cms.Process("TEST")
        process.DDDetectorESProducerFromDB = cms.ESSource(
            "DDDetectorESProducer",
            fromDB=cms.bool(True),
            label=cms.string("Extended"),
        )
        return process

    def test_replaces_only_the_geometry_source_and_records_contract(self):
        process = self.process()
        customiseShiftLssExternalGeometry(
            process,
            gdmlFile="PhysicsTools/ShiftLssGeometry/test/fixtures/bounded_external.gdml",
            modelOriginCm=(0.0, 0.0, 3000.0),
            modelToCms=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
            minimumAbsZCm=2750.0,
            checkOverlaps=False,
        )
        self.assertFalse(hasattr(process, "DDDetectorESProducerFromDB"))
        self.assertTrue(hasattr(process, "shiftLssGeometryESSource"))
        self.assertEqual(process.shiftLssGeometryESSource.geometryLabel.value(), "Extended")
        self.assertTrue(process.shiftLssGeometryContract.standardCmsGeometryPreserved.value())
        self.assertTrue(process.shiftLssGeometryContract.externalExtensionOnly.value())

    def test_fails_closed_without_standard_db_geometry(self):
        process = cms.Process("TEST")
        with self.assertRaisesRegex(RuntimeError, "standard DDDetector"):
            customiseShiftLssExternalGeometry(
                process,
                gdmlFile="missing.gdml",
                modelOriginCm=(0.0, 0.0, 3000.0),
                modelToCms=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
                minimumAbsZCm=2750.0,
            )

    def test_rejects_incomplete_coordinate_contract(self):
        process = self.process()
        with self.assertRaisesRegex(ValueError, "require 3 and 9"):
            customiseShiftLssExternalGeometry(
                process,
                gdmlFile="missing.gdml",
                modelOriginCm=(0.0, 0.0),
                modelToCms=(1.0,) * 9,
                minimumAbsZCm=2750.0,
            )


if __name__ == "__main__":
    unittest.main()
