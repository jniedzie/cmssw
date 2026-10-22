// Code style: AI-code. AI-generated particle-only dark-photon constructor.
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/ShiftDarkPhotonPhysics.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "G4Decay.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

ShiftDarkPhotonPhysics::ShiftDarkPhotonPhysics(double massGeV, double properCTauMm)
    : G4VPhysicsConstructor("ShiftDarkPhoton", 0),
      mass_(massGeV * CLHEP::GeV),
      lifetime_(properCTauMm * CLHEP::mm / CLHEP::c_light),
      width_(lifetime_ > 0. ? CLHEP::hbar_Planck / lifetime_ : 0.) {
  if (!(std::isfinite(mass_) && mass_ > 0. && std::isfinite(lifetime_) && lifetime_ > 0. &&
        std::isfinite(width_) && width_ > 0.)) {
    throw cms::Exception("DarkPhotonPhysics")
        << "DarkPhoton requires an explicit finite positive massGeV and properCTauMm; got " << massGeV << " GeV and "
        << properCTauMm << " mm. LAB distances must not be used as proper c-tau.";
  }
}

void ShiftDarkPhotonPhysics::checkDefinition(const G4ParticleDefinition& particle) const {
  if (particle.GetPDGEncoding() != 32 || particle.GetParticleName() != "dark_photon" ||
      particle.GetPDGMass() != mass_ || particle.GetPDGLifeTime() != lifetime_ || particle.GetPDGWidth() != width_ ||
      particle.GetPDGCharge() != 0. || particle.GetPDGiSpin() != 2 || particle.GetPDGiParity() != -1 ||
      particle.GetPDGiConjugation() != -1 || particle.GetPDGStable() || particle.IsShortLived() ||
      particle.GetAntiPDGEncoding() != 32 || particle.GetParticleType() != "boson" ||
      particle.GetParticleSubType() != "dark_photon") {
    throw cms::Exception("DarkPhotonPhysics") << "Conflicting Geant4 PDG32 definition: " << particle.GetParticleName()
                                             << ". One physical mass and proper lifetime are required per runtime.";
  }
}

void ShiftDarkPhotonPhysics::ConstructParticle() {
  auto* table = G4ParticleTable::GetParticleTable();
  if (const auto* existing = table->FindParticle(32)) {
    checkDefinition(*existing);
    return;
  }
  if (table->FindParticle("dark_photon") != nullptr) {
    throw cms::Exception("DarkPhotonPhysics") << "The Geant4 name dark_photon already belongs to another PDG code.";
  }
  auto* particle = new G4ParticleDefinition("dark_photon",
                                            mass_, width_, 0., 2, -1, -1, 0, 0, 0,
                                            "boson", 0, 0, 32, false, lifetime_, nullptr,
                                            false, "dark_photon", 32);
  checkDefinition(*particle);
}

void ShiftDarkPhotonPhysics::ConstructProcess() {
  const auto* particle = G4ParticleTable::GetParticleTable()->FindParticle(32);
  if (particle == nullptr) {
    throw cms::Exception("DarkPhotonPhysics") << "PDG32 was not defined before process construction.";
  }
  checkDefinition(*particle);
  const auto* manager = particle->GetProcessManager();
  if (manager == nullptr) {
    throw cms::Exception("DarkPhotonPhysics") << "PDG32 has no Geant4 process manager.";
  }
  int transportation = 0;
  int decay = 0;
  for (int index = 0; index < manager->GetProcessListLength(); ++index) {
    const auto* process = (*manager->GetProcessList())[index];
    transportation += process->GetProcessType() == fTransportation;
    decay += dynamic_cast<const G4Decay*>(process) != nullptr;
  }
  // The unchanged nominal list supplies both processes. Never add a second
  // decay process, a dark interaction, or a process on a Standard Model particle.
  if (transportation != 1 || decay != 1 || manager->GetProcessListLength() != 2) {
    throw cms::Exception("DarkPhotonPhysics")
        << "PDG32 requires exactly one inherited standard G4Decay and one transportation process; got " << decay
        << " decay, " << transportation << " transportation, " << manager->GetProcessListLength() << " total.";
  }
}
