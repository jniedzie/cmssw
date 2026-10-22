# Code style: AI-code. AI-generated synthetic-region native stacking fixture.
import FWCore.ParameterSet.Config as cms
from SimG4Core.Application.g4SimHits_cfi import g4SimHits

# Keep the installed time, neutrino and energy guards. Only name the synthetic
# CPU dead region used by the native fixture; this is not a detector config.
stack = g4SimHits.StackingAction.clone(DeadRegions=cms.vstring('synthetic_dead'))
