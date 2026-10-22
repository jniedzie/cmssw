// Code style: AI-code. AI-generated registration preserving nominal process slots.
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/PredefinedSignalDecayPhysics.h"
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/PredefinedSignalDecay.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include <memory>
#include <typeinfo>

void PredefinedSignalDecayPhysics::ConstructProcess() {
  // Retain the exact standard decay registration for every particle. Replace
  // only these four kaon instances before later hadronic constructors run.
  G4DecayPhysics::ConstructProcess();
  for (int pdg : {130, 310, 321, -321}) {
    auto* particle = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
    auto* manager = particle == nullptr ? nullptr : particle->GetProcessManager();
    if (manager == nullptr || manager->GetProcessListLength() == 0)
      throw cms::Exception("DarkPhotonPhysics") << "Missing nominal kaon processes for PDG " << pdg;
    const int slot = manager->GetProcessListLength() - 1;
    auto* original = dynamic_cast<G4Decay*>((*manager->GetProcessList())[slot]);
    int decayCount = 0;
    for (int index = 0; index <= slot; ++index)
      decayCount += dynamic_cast<G4Decay*>((*manager->GetProcessList())[index]) != nullptr;
    if (original == nullptr || typeid(*original) != typeid(G4Decay) || decayCount != 1 ||
        original->GetExtDecayer() != nullptr)
      throw cms::Exception("DarkPhotonPhysics")
          << "The unique standard decay must be the last kaon process at registration for PDG " << pdg;
    const int rest = manager->GetProcessOrdering(original, idxAtRest);
    const int along = manager->GetProcessOrdering(original, idxAlongStep);
    const int post = manager->GetProcessOrdering(original, idxPostStep);
    const bool active = manager->GetProcessActivation(original);
    auto replacement = std::make_unique<PredefinedSignalDecay>(original->GetProcessName());
    replacement->SetVerboseLevel(original->GetVerboseLevel());
    replacement->SetProcessSubType(original->GetProcessSubType());
    auto* updated = replacement.get();
    manager->RemoveProcess(original);  // Shared by other particles; do not delete.
    const int updatedSlot = manager->AddProcess(updated, rest, along, post);
    if (updatedSlot < 0)
      throw cms::Exception("DarkPhotonPhysics") << "Cannot register predefined kaon decay for PDG " << pdg;
    replacement.release();
    manager->SetProcessActivation(updated, active);
    if (updatedSlot != slot || manager->GetProcessListLength() != slot + 1 ||
        manager->GetProcessOrdering(updated, idxAtRest) != rest ||
        manager->GetProcessOrdering(updated, idxAlongStep) != along ||
        manager->GetProcessOrdering(updated, idxPostStep) != post ||
        updated->GetProcessType() != original->GetProcessType())
      throw cms::Exception("DarkPhotonPhysics") << "Predefined kaon decay changed the nominal process slot/order for PDG " << pdg;
  }
}
