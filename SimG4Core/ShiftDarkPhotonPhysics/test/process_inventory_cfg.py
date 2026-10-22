# Code style: AI-code. AI-generated explicit conversion-test point only.
import FWCore.ParameterSet.Config as cms
from SimG4Core.Application.g4SimHits_cfi import g4SimHits

physics = g4SimHits.Physics.clone(
    type="SimG4Core/Physics/FTFP_BERT_EMM",
    DarkPhoton=cms.PSet(massGeV=cms.double(0.3), properCTauMm=cms.double(100.0)),
)
