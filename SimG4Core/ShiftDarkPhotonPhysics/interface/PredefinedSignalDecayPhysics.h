// Code style: AI-code. AI-generated registration preserving nominal process slots.
#ifndef SimG4Core_ShiftDarkPhotonPhysics_PredefinedSignalDecayPhysics_h
#define SimG4Core_ShiftDarkPhotonPhysics_PredefinedSignalDecayPhysics_h

#include "G4DecayPhysics.hh"

class PredefinedSignalDecayPhysics : public G4DecayPhysics {
public:
  PredefinedSignalDecayPhysics(const G4String& name, G4int verbosity) : G4DecayPhysics(name, verbosity) {}
  void ConstructProcess() override;
};

#endif
