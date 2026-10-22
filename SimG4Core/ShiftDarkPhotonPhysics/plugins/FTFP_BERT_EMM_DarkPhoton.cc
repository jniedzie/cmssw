// Code style: AI-code. AI-generated optional wrapper of the exact CMS factory.
#include "SimG4Core/Physics/interface/PhysicsListFactory.h"
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/ShiftDarkPhotonPhysics.h"
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/PredefinedSignalDecayPhysics.h"
#include "FWCore/Utilities/interface/Exception.h"

class FTFPBertEmmDarkPhotonMaker : public PhysicsListMakerBase {
public:
  std::unique_ptr<PhysicsList> make(const edm::ParameterSet& parameters, SimActivityRegistry& registry) const override {
    const auto& point = parameters.getParameter<edm::ParameterSet>("DarkPhoton");
    auto particle = std::make_unique<ShiftDarkPhotonPhysics>(point.getParameter<double>("massGeV"),
                                                            point.getParameter<double>("properCTauMm"));
    // Use the ordinary CMS factory and its existing observer enrollment. The
    // returned object is the actual nominal list, with its parameters unchanged.
    std::unique_ptr<PhysicsListMakerBase> nominal(PhysicsListFactory::get()->create("SimG4Core/Physics/FTFP_BERT_EMM"));
    auto physics = nominal->make(parameters, registry);
    const auto* originalDecay = dynamic_cast<const G4DecayPhysics*>(physics->GetPhysics("Decay"));
    if (originalDecay == nullptr)
      throw cms::Exception("DarkPhotonPhysics") << "The exact nominal list must provide its standard decay constructor.";
    // The nominal factory passes this parameter to G4DecayPhysics. Its private
    // verbosity is separate from the public constructor verbosity property.
    auto decay = std::make_unique<PredefinedSignalDecayPhysics>(originalDecay->GetPhysicsName(),
                                                               parameters.getUntrackedParameter<int>("Verbosity", 0));
    decay->SetVerboseLevel(originalDecay->GetVerboseLevel());
    physics->ReplacePhysics(decay.get());
    if (physics->GetPhysics("Decay") != decay.get())
      throw cms::Exception("DarkPhotonPhysics") << "Cannot preserve the nominal decay-constructor slot.";
    decay.release();
    physics->RegisterPhysics(particle.get());
    if (physics->GetPhysics("ShiftDarkPhoton") != particle.get()) {
      throw cms::Exception("DarkPhotonPhysics") << "Cannot append the dark-photon constructor in the current G4 state.";
    }
    particle.release();  // The modular physics list owns registered constructors.
    return physics;
  }
};

DEFINE_EDM_PLUGIN(PhysicsListFactory, FTFPBertEmmDarkPhotonMaker, "SimG4Core/Physics/FTFP_BERT_EMM_DarkPhoton");
