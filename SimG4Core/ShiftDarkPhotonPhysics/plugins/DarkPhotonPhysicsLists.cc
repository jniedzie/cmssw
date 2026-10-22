// Code style: AI-code. AI-generated optional wrapper of the exact CMS factory.
#include "SimG4Core/Physics/interface/PhysicsListFactory.h"
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/ShiftDarkPhotonPhysics.h"
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/PredefinedSignalDecayPhysics.h"
#include "FWCore/Utilities/interface/Exception.h"

class DarkPhotonPhysicsListMaker : public PhysicsListMakerBase {
public:
  explicit DarkPhotonPhysicsListMaker(std::string nominalFactory) : nominalFactory_(std::move(nominalFactory)) {}

  std::unique_ptr<PhysicsList> make(const edm::ParameterSet& parameters, SimActivityRegistry& registry) const override {
    const auto& point = parameters.getParameter<edm::ParameterSet>("DarkPhoton");
    auto particle = std::make_unique<ShiftDarkPhotonPhysics>(point.getParameter<double>("massGeV"),
                                                            point.getParameter<double>("properCTauMm"));
    // Use the ordinary CMS factory and its existing observer enrollment. The
    // returned object is the actual nominal list, with its parameters unchanged.
    // EMH uses this exact type string to enable G4HepEm. Forwarding the wrapper
    // name would silently select different electromagnetic physics.
    auto nominalParameters = parameters;
    nominalParameters.addParameter<std::string>("type", nominalFactory_);
    std::unique_ptr<PhysicsListMakerBase> nominal(PhysicsListFactory::get()->create(nominalFactory_));
    auto physics = nominal->make(nominalParameters, registry);
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

private:
  const std::string nominalFactory_;
};

class FTFPBertEmmDarkPhotonMaker : public DarkPhotonPhysicsListMaker {
public:
  FTFPBertEmmDarkPhotonMaker() : DarkPhotonPhysicsListMaker("SimG4Core/Physics/FTFP_BERT_EMM") {}
};

class FTFPBertEmhDarkPhotonMaker : public DarkPhotonPhysicsListMaker {
public:
  FTFPBertEmhDarkPhotonMaker() : DarkPhotonPhysicsListMaker("SimG4Core/Physics/FTFP_BERT_EMH") {}
};

DEFINE_EDM_PLUGIN(PhysicsListFactory, FTFPBertEmmDarkPhotonMaker, "SimG4Core/Physics/FTFP_BERT_EMM_DarkPhoton");
DEFINE_EDM_PLUGIN(PhysicsListFactory, FTFPBertEmhDarkPhotonMaker, "SimG4Core/Physics/FTFP_BERT_EMH_DarkPhoton");
