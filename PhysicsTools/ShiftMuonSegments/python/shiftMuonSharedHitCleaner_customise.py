# Code style: AI-code. AI-generated opt-in configuration for SHIFT offline tracks.
import FWCore.ParameterSet.Config as cms


def customiseShiftSharedHitCleaner(process):
    """Enable native shared-hit cleaning only on the two SHIFT track producers."""
    process.shiftCosmicMuons.UseSharedHitTrajectoryCleaner = cms.bool(True)
    process.shiftTraversingMuons.UseSharedHitTrajectoryCleaner = cms.bool(True)
    return process
