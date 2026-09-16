"""Resolve real transport customisation without extending the CMS material map."""
import unittest

import FWCore.ParameterSet.Config as cms
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonTable
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_customise import customiseShiftLssTransport
from PhysicsTools.ShiftMuonSegments.shiftLssMagneticField_cfi import shiftLssUniformFieldElement


class MaterialBoundaryTest(unittest.TestCase):
    def process(self, detailed=False):
        process = cms.Process("BOUNDARY")
        process.shiftMuonTable = shiftMuonTable.clone(
            useDetailedMaterialPropagation=detailed,
            targetUseDetailedMaterialPropagation=True,
        )
        process.load("Configuration.StandardSequences.MagneticField_cff")
        return process

    def test_four_modes(self):
        for material in (False, True):
            for field in (False, True):
                with self.subTest(material=material, field=field):
                    process = self.process(material)
                    elements = [shiftLssUniformFieldElement(
                        "test_external", (-100., -100., 2000.), (100., 100., 15000.),
                        (0., 1., 0.))] if field else None
                    process = customiseShiftLssTransport(
                        process, materialBoundaryAbsZCm=14800., fieldElements=elements,
                        geant4eMaximumPathLengthCm=20000.)
                    self.assertEqual(process.shiftMuonTable.lssTransport.materialBoundaryAbsZCm.value(),
                                     14800. if material else 1100.)
                    self.assertEqual(hasattr(process, "shiftLssMagneticField"), field)
                    self.assertTrue(process.shiftMuonTable.targetUseDetailedMaterialPropagation.value())
                    if field:
                        self.assertEqual(process.shiftLssMagneticField.elements[0].maximumCm[2], 15000.)
                        self.assertEqual(process.shiftMuonTable.lssTransport.magneticFieldLabel.value(), "")

    def test_short_boundary_and_invalid_limits(self):
        process = customiseShiftLssTransport(self.process(), materialBoundaryAbsZCm=900.)
        self.assertEqual(process.shiftMuonTable.lssTransport.materialBoundaryAbsZCm.value(), 900.)
        for bad in (0., -1., float("inf"), float("nan")):
            with self.assertRaises(ValueError):
                customiseShiftLssTransport(self.process(), approximateMaterialBoundaryAbsZCm=bad)


if __name__ == "__main__":
    unittest.main()
