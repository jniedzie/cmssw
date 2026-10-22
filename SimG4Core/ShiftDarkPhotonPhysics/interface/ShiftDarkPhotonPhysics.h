// Code style: AI-code. AI-generated particle-only dark-photon constructor.
#ifndef SimG4Core_ShiftDarkPhotonPhysics_ShiftDarkPhotonPhysics_h
#define SimG4Core_ShiftDarkPhotonPhysics_ShiftDarkPhotonPhysics_h

#include "G4VPhysicsConstructor.hh"

class G4ParticleDefinition;

class ShiftDarkPhotonPhysics : public G4VPhysicsConstructor {
public:
  ShiftDarkPhotonPhysics(double massGeV, double properCTauMm);
  void ConstructParticle() override;
  void ConstructProcess() override;

private:
  void checkDefinition(const G4ParticleDefinition&) const;
  const double mass_;
  const double lifetime_;
  const double width_;
};

#endif
