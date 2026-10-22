import FWCore.ParameterSet.Config as cms


# Input adapter for a completed Pythia PDG-32 dark-photon graph. The generator
# must provide the physical width/lifetime, polarized decay and all vertices.
shiftDarkPhotonHepMC = cms.EDProducer(
    "ShiftDarkPhotonHepMCProducer",
    src=cms.InputTag("VtxSmeared"),
    massGeV=cms.double(15.0),
    # Match these to the explicit invariant-mass support in the generator card.
    generatedMassMinGeV=cms.double(12.0),
    generatedMassMaxGeV=cms.double(18.0),
    momentumRelativeTolerance=cms.double(1.e-6),
)
