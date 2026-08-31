#!/usr/bin/env python3

import unittest

import FWCore.ParameterSet.Config as cms

from IOMC.ShiftEventTiming.shiftSimHitTiming_customise import (
    customiseShiftSimHitReferenceTiming,
)


class ShiftSimHitTimingConfigurationTest(unittest.TestCase):
    def make_process(self, with_pileup=False):
        process = cms.Process("TEST")
        mix_parameters = {
            "mixObjects": cms.PSet(
                mixSH=cms.PSet(
                    input=cms.VInputTag(
                        cms.InputTag("g4SimHits", "MuonCSCHits"),
                        cms.InputTag("g4SimHits", "MuonDTHits"),
                        cms.InputTag("g4SimHits", "MuonRPCHits"),
                        cms.InputTag("g4SimHits", "MuonGEMHits"),
                    )
                )
            )
        }
        if with_pileup:
            mix_parameters["input"] = cms.PSet()
        process.mix = cms.EDProducer("MixingModule", **mix_parameters)
        process.simMuonDTDigis = cms.EDProducer(
            "TestDT", InputCollection=cms.string("g4SimHitsMuonDTHits"),
            InputCollectionPU=cms.string("g4SimHitsMuonDTHits")
        )
        process.simMuonCSCDigis = cms.EDProducer(
            "TestCSC", InputCollection=cms.string("g4SimHitsMuonCSCHits"),
            InputCollectionPU=cms.string("g4SimHitsMuonCSCHits")
        )
        process.simMuonRPCDigis = cms.EDProducer(
            "TestRPC", InputCollection=cms.string("g4SimHitsMuonRPCHits"),
            InputCollectionPU=cms.string("g4SimHitsMuonRPCHits")
        )
        process.simMuonGEMDigis = cms.EDProducer(
            "TestGEM", inputCollection=cms.string("g4SimHitsMuonGEMHits"),
            inputCollectionPU=cms.string("g4SimHitsMuonGEMHits")
        )
        process.pdigiTask = cms.Task(
            process.simMuonDTDigis,
            process.simMuonCSCDigis,
            process.simMuonRPCDigis,
            process.simMuonGEMDigis,
        )
        process.output = cms.OutputModule(
            "PoolOutputModule", outputCommands=cms.untracked.vstring("drop *")
        )
        return process

    def test_redirects_only_muon_simhits(self):
        process = customiseShiftSimHitReferenceTiming(
            self.make_process(), bxOffset=1, phaseNs=2.5
        )
        tags = list(process.mix.mixObjects.mixSH.input)
        self.assertEqual({tag.moduleLabel for tag in tags}, {"shiftSimHitTime"})
        self.assertEqual(process.shiftSimHitTime.bxOffset.value(), 1)
        self.assertEqual(process.shiftSimHitTime.phaseNs.value(), 2.5)
        self.assertEqual(
            process.simMuonDTDigis.InputCollection.value(),
            "shiftSimHitTimeMuonDTHits",
        )
        self.assertIn("keep *_shiftSimHitTime_*_*", process.output.outputCommands)

    def test_rejects_pileup(self):
        with self.assertRaisesRegex(RuntimeError, "restricted to no-pileup"):
            customiseShiftSimHitReferenceTiming(self.make_process(with_pileup=True))


if __name__ == "__main__":
    unittest.main()
