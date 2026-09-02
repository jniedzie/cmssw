#!/usr/bin/env python3

import hashlib
from pathlib import Path
import unittest

import FWCore.ParameterSet.Config as cms

from PhysicsTools.ShiftMuonSegments.shiftLssMagneticField_cfi import (
    shiftLssFlukaMap2DFieldElement,
    shiftLssQuadrupoleFieldElement,
    shiftLssUniformFieldElement,
)
from PhysicsTools.ShiftMuonSegments.shiftLssIr1AtlasProxy_cff import (
    shiftLssIr1AtlasProxyFieldElements,
)
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cfi import shiftMuonTable
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cff import addShiftMuonSegments
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_customise import (
    customiseShiftLssMagneticField,
    customiseShiftLssTransport,
)


class ShiftLssConfigurationTest(unittest.TestCase):
    data_directory = Path(__file__).resolve().parents[1] / "data" / "lss" / "ir1_atlas_proxy"

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
        self.assertEqual(process.shiftLssFieldContract.elementCount.value(), 2)

    def test_simulation_field_does_not_require_reconstruction_table(self):
        process = self.process()
        del process.shiftMuonTable
        elements = [
            shiftLssUniformFieldElement(
                "dipole", (-10.0, -10.0, -100.0), (10.0, 10.0, 100.0), (0.0, 1.2, 0.0)
            )
        ]
        customiseShiftLssMagneticField(
            process,
            baseMagneticFieldProducer="baseField",
            fieldElements=elements,
        )
        self.assertEqual(len(process.shiftLssMagneticField.elements), 1)
        self.assertFalse(hasattr(process, "shiftMuonTable"))

    def test_detailed_target_material_is_selected_before_sequence_construction(self):
        process = cms.Process("TEST")
        process.nanoAOD_step = cms.Path()
        addShiftMuonSegments(process, useDetailedMaterialPropagation=True)
        self.assertTrue(process.shiftMuonTable.useDetailedMaterialPropagation.value())
        self.assertTrue(hasattr(process, "shiftMuonGeant4Geometry"))

    def test_fluka_map_element_keeps_all_runtime_inputs(self):
        element = shiftLssFlukaMap2DFieldElement(
            "mqxa",
            (-25.0, -25.0, -100.0),
            (25.0, 25.0, 100.0),
            "/tmp/MQXA.dat",
            1.8928521203,
            originCm=(0.0, 0.0, 2615.0),
        )
        self.assertEqual(element.type.value(), "flukaMap2D")
        self.assertEqual(element.mapFile.value(), "/tmp/MQXA.dat")
        self.assertAlmostEqual(element.fieldScaleTesla.value(), 1.8928521203)
        self.assertEqual(list(element.originCm), [0.0, 0.0, 2615.0])

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
            customiseShiftLssMagneticField(
                process,
                baseMagneticFieldProducer="baseField",
                fieldElements=[],
            )
        with self.assertRaises(ValueError):
            customiseShiftLssTransport(process, materialBoundaryAbsZCm=0.0)

    def test_ir1_proxy_requires_and_applies_coordinate_transform(self):
        rotation = (0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0)
        elements = shiftLssIr1AtlasProxyFieldElements(
            modelOriginCm=(10.0, 20.0, 30.0),
            modelToCms=rotation,
            fieldScale=-1.0,
            includeConstantFields=False,
        )
        self.assertEqual(len(elements), 17)
        self.assertEqual(
            [element.name.value() for element in elements[:4]],
            ["MQXA.1R1.outer", "MQXA.1R1.aperture", "MQXB.A2R1", "MQXB.B2R1"],
        )
        self.assertEqual(list(elements[0].originCm), [2625.0, 20.0, 30.0])
        self.assertEqual(list(elements[0].localToGlobal), list(rotation))
        self.assertAlmostEqual(elements[0].fieldScaleTesla.value(), -1.8928521203)
        self.assertEqual(elements[0].boundsShape.value(), "cylinderZ")
        self.assertEqual(elements[0].innerRadiusCm.value(), 3.45)
        self.assertEqual(elements[1].outerRadiusCm.value(), 3.3)
        self.assertEqual(
            elements[0].mapFile.value(),
            "PhysicsTools/ShiftMuonSegments/data/lss/ir1_atlas_proxy/MQXA.dat",
        )
        self.assertEqual(list(elements[-1].minimumCm), [-45.7, -53.7, -170.0])
        with self.assertRaises(ValueError):
            shiftLssIr1AtlasProxyFieldElements(
                modelOriginCm=(0.0, 0.0, 0.0),
                modelToCms=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0),
            )

        complete_elements = shiftLssIr1AtlasProxyFieldElements(
            modelOriginCm=(0.0, 0.0, 0.0),
            modelToCms=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        )
        self.assertEqual(len(complete_elements), 29)
        self.assertEqual(complete_elements[-3].name.value(), "MCBCH.5R1.outer")
        self.assertEqual(
            list(complete_elements[-3].excludedCylindersCm),
            [-9.7, 0.0, 2.8, 9.7, 0.0, 2.8],
        )
        with self.assertRaises(ValueError):
            shiftLssIr1AtlasProxyFieldElements(
                modelOriginCm=(0.0, 0.0, 0.0),
                modelToCms=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
                includeMappedFields=False,
                includeConstantFields=False,
            )

    def test_packaged_proxy_maps_match_frozen_checksums_and_shapes(self):
        expected_shapes = {
            "MBXW.dat": ("KICKINT", 199, 100),
            "MQTL.dat": ("QUADINT", 60, 60),
            "MQXA.dat": ("QUADINT", 97, 97),
            "MQXB.dat": ("QUADINT", 89, 89),
            "MQYana.dat": ("QUAD", 0, 0),
        }
        checksums = {}
        for line in (self.data_directory / "FIELD_MAP_SHA256SUMS").read_text().splitlines():
            digest, name = line.split()
            checksums[name] = digest
        self.assertEqual(set(checksums), set(expected_shapes))

        for name, expected_shape in expected_shapes.items():
            path = self.data_directory / name
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), checksums[name])
            cards = {}
            data_rows = 0
            reading_data = False
            for raw_line in path.read_text(encoding="ascii").splitlines():
                line = raw_line.split("#", 1)[0].strip()
                if not line:
                    continue
                if reading_data:
                    self.assertEqual(len(line.split()), 2)
                    data_rows += 1
                    continue
                tokens = line.split()
                cards[tokens[0]] = tokens[1:]
                reading_data = tokens[0] == "DATA"
            field_type, bins_x, bins_y = expected_shape
            self.assertEqual(cards["TYPE"], [field_type])
            if field_type == "QUAD":
                self.assertEqual(data_rows, 0)
            else:
                self.assertEqual(tuple(map(int, (cards["XGRID"][2], cards["YGRID"][2]))),
                                 (bins_x, bins_y))
                self.assertEqual(data_rows, bins_x * bins_y)


if __name__ == "__main__":
    unittest.main()
