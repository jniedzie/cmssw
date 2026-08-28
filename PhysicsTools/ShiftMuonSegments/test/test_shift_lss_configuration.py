#!/usr/bin/env python3

import unittest

import FWCore.ParameterSet.Config as cms

from PhysicsTools.ShiftMuonSegments.shiftLssMagneticField_cfi import (
    shiftLssQuadrupoleFieldElement,
    shiftLssUniformFieldElement,
)
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonTable
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_customise import customiseShiftLssTransport


class ShiftLssConfigurationTest(unittest.TestCase):
    def process(self):
        process = cms.Process("TEST")
        process.shiftMuonTable = shiftMuonTable.clone()
        process.baseField = cms.ESProducer(
            "UniformMagneticFieldESProducer",
            ZFieldInTesla=cms.double(3.8),
            label=cms.untracked.string(""),
        )
        return process

    def test_shared_simulation_and_reconstruction_field(self):
        process = self.process()
        elements = [
            shiftLssUniformFieldElement(
                "dipole", (-10.0, -10.0, -100.0), (10.0, 10.0, 100.0), (0.0, 1.2, 0.0)
            ),
            shiftLssQuadrupoleFieldElement(
                "quadrupole", (-5.0, -5.0, -50.0), (5.0, 5.0, 50.0), 0.02
            ),
        ]
        customiseShiftLssTransport(
            process,
            baseMagneticFieldProducer="baseField",
            fieldElements=elements,
            materialBoundaryAbsZCm=14800.0,
            geant4eMaximumPathLengthCm=20000.0,
        )
        self.assertEqual(process.shiftMuonTable.lssTransport.magneticFieldLabel.value(), "")
        self.assertEqual(process.shiftMuonTable.lssTransport.materialBoundaryAbsZCm.value(), 14800.0)
        self.assertEqual(len(process.shiftLssMagneticField.elements), 2)
        self.assertEqual(process.baseField.label.value(), "shiftLssBaseMagneticField")
        self.assertEqual(process.shiftLssMagneticField.outputLabel.value(), "")

    def test_rejects_incomplete_field_configuration(self):
        process = self.process()
        with self.assertRaises(ValueError):
            customiseShiftLssTransport(
                process,
                magneticFieldLabel="namedComposite",
                baseMagneticFieldProducer="baseField",
                fieldElements=[],
            )
        with self.assertRaises(ValueError):
            customiseShiftLssTransport(process, materialBoundaryAbsZCm=0.0)


if __name__ == "__main__":
    unittest.main()
