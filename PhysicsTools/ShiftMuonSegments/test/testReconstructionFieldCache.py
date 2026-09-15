"""Offline field derivatives must resolve nearby points without changing simulation."""
import FWCore.ParameterSet.Config as cms
from SimG4Core.Application.g4SimHits_cfi import g4SimHits
from PhysicsTools.ShiftMuonSegments.shiftMuonSegments_cff import addShiftMuonSegments

simulation_before = g4SimHits.MagneticField.dumpPython()
process = cms.Process("FIELDCACHE")
process.nanoAOD_step = cms.Path()
process = addShiftMuonSegments(process, targetUseDetailedMaterialPropagation=True)
assert process.shiftMuonGeant4Geometry.MagneticField.delta.value() == 0.
assert g4SimHits.MagneticField.dumpPython() == simulation_before

# The reconstruction clone must not alias the simulation PSet.
process.shiftMuonGeant4Geometry.MagneticField.delta = cms.double(.25)
assert g4SimHits.MagneticField.dumpPython() == simulation_before
