// Code style: AI-code. AI-generated opt-in preassigned-decay conversion repair.
#include "SimG4Core/Generators/interface/PreassignedDecayKinematics.h"
#include "SimG4Core/Generators/interface/SourceAncestry.h"
#include "G4DecayProducts.hh"
#include "G4DynamicParticle.hh"
#include "G4LorentzVector.hh"
#include "G4PrimaryParticle.hh"
#include "G4Track.hh"
#include <cmath>
#include <stdexcept>

namespace {
  G4LorentzVector sourceMomentum(const G4PrimaryParticle& primary) {
    const G4LorentzVector value(primary.GetMomentum(), primary.GetTotalEnergy());
    const double mass = primary.GetMass();
    if (!std::isfinite(value.x()) || !std::isfinite(value.y()) || !std::isfinite(value.z()) ||
        !std::isfinite(value.e()) || !std::isfinite(mass) || mass < 0. || value.e() < mass || value.e() < 0.)
      throw std::runtime_error("Preassigned decay repair requires finite physical source four-momenta");
    return value;
  }

  void setMomentum(G4DynamicParticle& particle, const G4LorentzVector& value, double mass) {
    // SetMass alone retains the old kinetic energy. SetMomentum then rebuilds
    // that energy using the generated mass while preserving the source vector.
    // Never infer a new mass from a boosted, nearly cancelling subtraction.
    particle.SetMass(mass);
    particle.SetMomentum(value.vect());
  }

  void repair(G4DynamicParticle& particle) {
    const auto* primary = particle.GetPrimaryParticle();
    if (!primary || primary->GetPDGcode() != particle.GetPDGcode())
      throw std::runtime_error("Preassigned decay repair lost its source particle identity");
    const auto parentLab = sourceMomentum(*primary);
    const double parentMass = primary->GetMass();
    setMomentum(particle, parentLab, parentMass);
    // The products belong to this mutable DynamicParticle. Geant4 exposes
    // their accessor as const even though its native transformer fills them.
    auto* products = const_cast<G4DecayProducts*>(particle.GetPreAssignedDecayProducts());
    if (!products)
      return;
    if (parentMass <= 0.)
      throw std::runtime_error("A preassigned decaying source particle needs positive mass");
    for (int i = 0; i < products->entries(); ++i) {
      auto* child = (*products)[i];
      repair(*child);
      const auto* childPrimary = child->GetPrimaryParticle();
      auto childRest = sourceMomentum(*childPrimary);
      childRest.boost(-parentLab.boostVector());
      setMomentum(*child, childRest, childPrimary->GetMass());
    }
    // Keep the original state as a rest-frame parent snapshot. Detach the
    // copied nested tree before SetParentParticle copies this snapshot again.
    G4DynamicParticle restParent(particle);
    const auto* copiedProducts = restParent.GetPreAssignedDecayProducts();
    restParent.SetPreAssignedDecayProducts(nullptr);
    delete copiedProducts;
    restParent.SetKineticEnergy(0.);
    restParent.SetMomentumDirection(G4ThreeVector(0., 0., 1.));
    products->SetParentParticle(restParent);
  }
}

void sim::preservePreassignedDecayKinematics(G4Track& track) {
  if (track.GetParentID() != 0 || track.GetCurrentStepNumber() != 0 || track.GetTrackLength() != 0. ||
      track.GetLocalTime() != 0. || track.GetProperTime() != 0.)
    throw std::runtime_error("Preassigned decay repair must precede primary transport");
  sim::sourceAncestry::markKinematicsRestored(&track);
  // This mutable track owns its DynamicParticle; Geant4 offers only a const
  // accessor. The function is restricted above to pre-transport primaries.
  repair(*const_cast<G4DynamicParticle*>(track.GetDynamicParticle()));
}

void sim::preservePreassignedDaughterMass(G4Track& track) {
  if (track.GetParentID() <= 0 || track.GetCurrentStepNumber() != 0 || track.GetTrackLength() != 0. ||
      track.GetLocalTime() != 0. || track.GetProperTime() != 0.)
    throw std::runtime_error("Preassigned daughter mass repair must precede secondary transport");
  auto& particle = *const_cast<G4DynamicParticle*>(track.GetDynamicParticle());
  const auto* primary = particle.GetPrimaryParticle();
  if (!primary || primary->GetPDGcode() != particle.GetPDGcode())
    throw std::runtime_error("Preassigned daughter mass repair lost its source particle identity");
  const double mass = primary->GetMass();
  const auto momentum = particle.GetMomentum();
  if (!std::isfinite(mass) || mass < 0. || !std::isfinite(momentum.x()) ||
      !std::isfinite(momentum.y()) || !std::isfinite(momentum.z()))
    throw std::runtime_error("Preassigned daughter mass repair requires finite physical kinematics");
  sim::sourceAncestry::markKinematicsRestored(&track);
  particle.SetMass(mass);
  particle.SetMomentum(momentum);
}
