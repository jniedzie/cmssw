// Code style: AI-code. AI-generated covariance repair for nested signal decays.
#ifndef SimG4Core_ShiftDarkPhotonPhysics_PredefinedSignalDecay_h
#define SimG4Core_ShiftDarkPhotonPhysics_PredefinedSignalDecay_h

#include "G4Decay.hh"

// All interaction lengths, lifetimes, standard decays and transport are inherited.
// Only a predefined kaon chain containing PDG32 needs its nested frames updated.
class PredefinedSignalDecay : public G4Decay {
public:
  explicit PredefinedSignalDecay(const G4String& name = "Decay") : G4Decay(name) {}

protected:
  void DaughterPolarization(const G4Track&, G4DecayProducts*) override;
};

#endif
