// Code style: AI-code. AI-generated covariance repair for nested signal decays.
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/PredefinedSignalDecay.h"
#include "G4DecayProducts.hh"
#include "G4DynamicParticle.hh"
#include "G4Exception.hh"
#include "G4LorentzVector.hh"
#include "G4Track.hh"

#include <cmath>
#include <memory>

namespace {
  G4LorentzVector transformed(G4LorentzVector momentum, const G4ThreeVector& oldBeta, const G4ThreeVector& newBeta) {
    momentum.boost(-oldBeta);
    momentum.boost(newBeta);
    return momentum;
  }

  void transformNested(G4DynamicParticle& particle, const G4ThreeVector& oldBeta, const G4ThreeVector& newBeta) {
    const auto* original = particle.GetPreAssignedDecayProducts();
    if (original == nullptr)
      return;
    auto replacement = std::make_unique<G4DecayProducts>(*original);
    G4DynamicParticle parent(*replacement->GetParentParticle());
    parent.Set4Momentum(transformed(parent.Get4Momentum(), oldBeta, newBeta));
    replacement->SetParentParticle(parent);
    for (int index = 0; index < replacement->entries(); ++index) {
      auto* daughter = (*replacement)[index];
      daughter->Set4Momentum(transformed(daughter->Get4Momentum(), oldBeta, newBeta));
      transformNested(*daughter, oldBeta, newBeta);
    }
    // SetPreAssignedDecayProducts transfers ownership but does not delete the
    // previous graph. Work on a deep copy, then replace exactly that owned tree.
    particle.SetPreAssignedDecayProducts(replacement.release());
    delete original;
  }
}  // namespace

void PredefinedSignalDecay::DaughterPolarization(const G4Track& track, G4DecayProducts* products) {
  const auto* particle = track.GetDynamicParticle();
  const int pdg = particle->GetPDGcode();
  if (pdg != 130 && pdg != 310 && std::abs(pdg) != 321)
    return;
  const auto* original = particle->GetPreAssignedDecayProducts();
  if (original == nullptr)
    return;
  bool hasSignal = false;
  for (int index = 0; index < original->entries(); ++index)
    hasSignal = hasSignal || (*original)[index]->GetPDGcode() == 32;
  if (!hasSignal)
    return;

  const auto oldBeta = original->GetParentParticle()->Get4Momentum().boostVector();
  const auto newBeta = particle->Get4Momentum().boostVector();
  if (!std::isfinite(oldBeta.mag2()) || !std::isfinite(newBeta.mag2()) || oldBeta.mag2() >= 1. || newBeta.mag2() >= 1.) {
    G4Exception("PredefinedSignalDecay", "ShiftInvalidDecayBoost", FatalException,
                "A predefined signal kaon needs finite timelike old and current momenta");
    return;
  }
  // Base G4Decay has already boosted the immediate daughters. Its empty hook
  // leaves deeper decay frames unchanged. Apply the same Lorentz transform to
  // their parent copies and descendants, preserving spin correlations and all
  // preassigned proper times. Ordinary kaon decays never enter this branch.
  for (int index = 0; index < products->entries(); ++index)
    transformNested(*(*products)[index], oldBeta, newBeta);
}
