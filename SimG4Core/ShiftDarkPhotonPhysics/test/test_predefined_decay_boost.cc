// Code style: AI-code. AI-generated native nested-frame covariance regression.
// Uses actual G4DecayProducts and the same virtual hook called by G4Decay.
// This does not test material transport or change any Run-3 physics setting.
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/PredefinedSignalDecay.h"
#include "G4BosonConstructor.hh"
#include "G4DecayProducts.hh"
#include "G4DynamicParticle.hh"
#include "G4LeptonConstructor.hh"
#include "G4LorentzVector.hh"
#include "G4MesonConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
  class HookControl : public PredefinedSignalDecay {
  public:
    using PredefinedSignalDecay::DaughterPolarization;
  };

  void require(bool condition, const char* message) {
    if (!condition)
      throw std::runtime_error(message);
  }

  double difference(const G4LorentzVector& first, const G4LorentzVector& second) {
    double maximum = 0.;
    for (int component = 0; component < 4; ++component)
      maximum = std::max(maximum, std::abs(first[component] - second[component]));
    return maximum / GeV;
  }

  G4LorentzVector covariant(G4LorentzVector momentum, const G4LorentzVector& oldParent,
                           const G4ThreeVector& newBeta) {
    momentum.boost(-oldParent.boostVector());
    momentum.boost(newBeta);
    return momentum;
  }

  G4DynamicParticle* copyWithDecay(const G4DynamicParticle& original) {
    // G4DynamicParticle's copy constructor deliberately omits predefined
    // decay products/time. Mirror G4DecayProducts' explicit deep-copy policy.
    auto* copy = new G4DynamicParticle(original);
    copy->SetPreAssignedDecayProperTime(original.GetPreAssignedDecayProperTime());
    if (original.GetPreAssignedDecayProducts() != nullptr)
      copy->SetPreAssignedDecayProducts(new G4DecayProducts(*original.GetPreAssignedDecayProducts()));
    return copy;
  }
}  // namespace

int main() {
  try {
    G4BosonConstructor().ConstructParticle();
    G4LeptonConstructor().ConstructParticle();
    G4MesonConstructor().ConstructParticle();
    auto* table = G4ParticleTable::GetParticleTable();
    auto* darkDefinition = new G4ParticleDefinition("test_covariant_dark", .3 * GeV, 0., 0., 2, -1, -1,
        0, 0, 0, "boson", 0, 0, 32, false, 1. * ns, nullptr, false, "test", 32);
    table->SetReadiness();
    const double darkMass = darkDefinition->GetPDGMass();
    const double muonMass = table->FindParticle(13)->GetPDGMass();
    int controls = 0;
    double largestMissingRotation = 0., largestCorrectedDifference = 0.;
    for (int parentId : {130, 310, 321, -321}) {
      const double parentMass = table->FindParticle(parentId)->GetPDGMass();
      const int recoilId = std::abs(parentId) == 321 ? (parentId > 0 ? 211 : -211) : 22;
      const double recoilMass = table->FindParticle(recoilId)->GetPDGMass();
      const G4LorentzVector oldParent(0., 0., 1. * GeV, std::hypot(parentMass, 1. * GeV));
      const double darkEnergy = (parentMass * parentMass + darkMass * darkMass - recoilMass * recoilMass) / (2. * parentMass);
      const double q = std::sqrt(darkEnergy * darkEnergy - darkMass * darkMass);
      G4LorentzVector oldDark(q, 0., 0., darkEnergy), oldRecoil(-q, 0., 0., parentMass - darkEnergy);
      oldDark.boost(oldParent.boostVector());
      oldRecoil.boost(oldParent.boostVector());
      const double pMuon = std::sqrt(darkMass * darkMass / 4. - muonMass * muonMass);
      for (const auto& muonDirection : {G4ThreeVector(1., 0., 0.), G4ThreeVector(0., 1., 0.),
                                      G4ThreeVector(.3, .4, std::sqrt(.75))}) {
        G4LorentzVector oldMuon(pMuon * muonDirection, darkMass / 2.);
        oldMuon.boost(oldDark.boostVector());
        for (const auto& newBeta : {oldParent.boostVector(), G4ThreeVector(0., 0., .6),
                                   G4ThreeVector(.2, .1, .6), G4ThreeVector(0., 0., 0.)}) {
          G4DynamicParticle originalParent(table->FindParticle(parentId), oldParent);
          auto originalProducts = std::make_unique<G4DecayProducts>(originalParent);
          auto* dark = new G4DynamicParticle(darkDefinition, oldDark);
          auto* nested = new G4DecayProducts(*dark);
          nested->PushProducts(new G4DynamicParticle(table->FindParticle(13), oldMuon));
          nested->PushProducts(new G4DynamicParticle(table->FindParticle(-13), oldDark - oldMuon));
          dark->SetPreAssignedDecayProducts(nested);
          dark->SetPreAssignedDecayProperTime(7. * ns);
          originalProducts->PushProducts(dark);
          originalProducts->PushProducts(new G4DynamicParticle(table->FindParticle(recoilId), oldRecoil));
          G4DecayProducts products(*originalProducts);
          products.Boost(newBeta.x(), newBeta.y(), newBeta.z());
          G4DecayProducts uncorrected(*products[0]->GetPreAssignedDecayProducts());
          uncorrected.Boost(products[0]->GetTotalEnergy(), products[0]->GetMomentumDirection());
          const auto expectedMuon = covariant(oldMuon, oldParent, newBeta);
          largestMissingRotation = std::max(largestMissingRotation, difference(uncorrected[0]->Get4Momentum(), expectedMuon));

          auto* transportedParent = new G4DynamicParticle(table->FindParticle(parentId), covariant(oldParent, oldParent, newBeta));
          transportedParent->SetPreAssignedDecayProducts(originalProducts.release());
          transportedParent->SetPreAssignedDecayProperTime(11. * ns);
          G4Track track(transportedParent, 1. * ns, G4ThreeVector());
          HookControl process;
          process.DaughterPolarization(track, &products);
          const auto* corrected = products[0]->GetPreAssignedDecayProducts();
          require(corrected != nullptr && corrected->entries() == 2, "Nested daughters lost");
          require(products[0]->GetPreAssignedDecayProperTime() == 7. * ns &&
                  track.GetDynamicParticle()->GetPreAssignedDecayProperTime() == 11. * ns, "Proper time changed");
          require(difference(corrected->GetParentParticle()->Get4Momentum(), products[0]->Get4Momentum()) < 1.e-10,
                  "Nested parent frame differs from transported parent");
          G4DecayProducts finalMuons(*corrected);
          finalMuons.Boost(products[0]->GetTotalEnergy(), products[0]->GetMomentumDirection());
          const double residual = difference(finalMuons[0]->Get4Momentum(), expectedMuon);
          largestCorrectedDifference = std::max(largestCorrectedDifference, residual);
          require(residual < 1.e-10, "Correlated muon fails full-chain Lorentz covariance");
          require(difference(finalMuons[0]->Get4Momentum() + finalMuons[1]->Get4Momentum(), products[0]->Get4Momentum()) < 1.e-10,
                  "Corrected dimuon momentum does not close");
          require(difference(track.GetDynamicParticle()->GetPreAssignedDecayProducts()->operator[](0)->Get4Momentum(), oldDark) < 1.e-12,
                  "Original preassigned source graph mutated");

          // Exercise the inherited real decay routine, not just a direct hook
          // call. Its first decay must hand the corrected frame to ordinary
          // G4Decay for the subsequent A' -> mu mu decay.
          G4Track actualParent(copyWithDecay(*track.GetDynamicParticle()), 1. * ns, G4ThreeVector());
          // G4Decay copies the parent's touchable into each secondary. Its
          // reference-counted handle must be initialized even in this harness.
          // An empty history supplies ownership only, not physical geometry.
          G4TouchableHandle touchable(new G4TouchableHistory());
          actualParent.SetTouchableHandle(touchable);
          G4Step step;
          step.SetTrack(&actualParent);
          actualParent.SetStep(&step);
          auto* parentChange = process.PostStepDoIt(actualParent, step);
          const G4DynamicParticle* actualDark = nullptr;
          for (int index = 0; index < parentChange->GetNumberOfSecondaries(); ++index) {
            const auto* secondary = parentChange->GetSecondary(index)->GetDynamicParticle();
            if (secondary->GetPDGcode() == 32)
              actualDark = secondary;
          }
          require(actualDark != nullptr, "Actual kaon decay lost its A'");
          G4Track actualDarkTrack(copyWithDecay(*actualDark), 1. * ns, G4ThreeVector());
          actualDarkTrack.SetTouchableHandle(touchable);
          G4Step darkStep;
          darkStep.SetTrack(&actualDarkTrack);
          actualDarkTrack.SetStep(&darkStep);
          G4Decay ordinaryDarkDecay;
          auto* darkChange = ordinaryDarkDecay.PostStepDoIt(actualDarkTrack, darkStep);
          require(darkChange->GetNumberOfSecondaries() == 2, "Actual A' decay lost its muons");
          bool foundMuon = false;
          for (int index = 0; index < darkChange->GetNumberOfSecondaries(); ++index) {
            const auto* muon = darkChange->GetSecondary(index)->GetDynamicParticle();
            if (muon->GetPDGcode() == 13) {
              require(difference(muon->Get4Momentum(), expectedMuon) < 1.e-10,
                      "Actual inherited decay routines violate full-chain covariance");
              foundMuon = true;
            }
          }
          require(foundMuon, "Actual inherited decay has wrong daughter charges");
          for (int index = 0; index < darkChange->GetNumberOfSecondaries(); ++index)
            delete darkChange->GetSecondary(index);
          darkChange->Clear();
          for (int index = 0; index < parentChange->GetNumberOfSecondaries(); ++index)
            delete parentChange->GetSecondary(index);
          parentChange->Clear();

          // With no predefined signal graph the hook must leave SM products
          // byte-for-byte equal in their four-momenta and proper-time values.
          G4DynamicParticle smParent(table->FindParticle(parentId), oldParent);
          G4Track smTrack(new G4DynamicParticle(smParent), 0., G4ThreeVector());
          const auto smBefore = products[0]->Get4Momentum();
          process.DaughterPolarization(smTrack, &products);
          require(difference(smBefore, products[0]->Get4Momentum()) == 0., "No-signal SM path changed");
          ++controls;
        }
      }
    }
    require(largestMissingRotation > 1.e-3, "Regression did not expose the original nested-frame error");
    std::cout << std::setprecision(17) << "{\"controls\":" << controls
              << ",\"max_uncorrected_difference_gev\":" << largestMissingRotation
              << ",\"max_corrected_difference_gev\":" << largestCorrectedDifference
              << ",\"proper_times_preserved\":true,\"source_graphs_unchanged\":true,\"material_transport_tested\":false}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
