// Code style: AI-code. AI-generated geometry-free native process inventory.
#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/standard.h"
#include "FWCore/PythonParameterSet/interface/MakePyBind11ParameterSets.h"
#include "SimG4Core/Notification/interface/SimActivityRegistry.h"
#include "SimG4Core/Physics/interface/PhysicsListFactory.h"
#include "SimG4Core/ShiftDarkPhotonPhysics/interface/PredefinedSignalDecay.h"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4Version.hh"
#include "G4HepEmTrackingManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"

#include <fstream>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <link.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <typeinfo>

class NoopExactNominal : public PhysicsListMakerBase {
public:
  std::unique_ptr<PhysicsList> make(const edm::ParameterSet& parameters, SimActivityRegistry& registry) const override {
    std::unique_ptr<PhysicsListMakerBase> nominal(PhysicsListFactory::get()->create(parameters.getParameter<std::string>("type")));
    return nominal->make(parameters,registry);
  }
};

void dumpParticle(std::ostream& out, const G4ParticleDefinition& particle) {
  out << "{\"name\":" << std::quoted(std::string(particle.GetParticleName()))
      << ",\"pdg_id\":" << particle.GetPDGEncoding() << ",\"mass_gev\":" << particle.GetPDGMass()/CLHEP::GeV
      << ",\"proper_lifetime_ns\":" << particle.GetPDGLifeTime()/CLHEP::ns
      << ",\"twice_spin\":" << particle.GetPDGiSpin() << ",\"charge_e\":" << particle.GetPDGCharge()/CLHEP::eplus
      << ",\"stable\":" << (particle.GetPDGStable()?"true":"false")
      << ",\"short_lived\":" << (particle.IsShortLived()?"true":"false")
      << ",\"g4hepem_tracking_manager\":" << (dynamic_cast<const G4HepEmTrackingManager*>(particle.GetTrackingManager())?"true":"false")
      << ",\"tracking_manager_class\":";
  const auto* tracking = particle.GetTrackingManager();
  if (tracking) out << std::quoted(std::string(typeid(*tracking).name())); else out << "null";
  out
      << ",\"processes\":";
  auto* manager = particle.GetProcessManager();
  if (manager == nullptr) {
    out << "null}";
    return;
  }
  out << '[';
  for (int index=0; index<manager->GetProcessListLength(); ++index) {
    if (index) out << ',';
    auto* process = (*manager->GetProcessList())[index];
    out << "{\"name\":" << std::quoted(std::string(process->GetProcessName()))
        << ",\"type\":" << process->GetProcessType() << ",\"subtype\":" << process->GetProcessSubType()
        << ",\"at_rest_order\":" << manager->GetProcessOrdering(process,idxAtRest)
        << ",\"along_step_order\":" << manager->GetProcessOrdering(process,idxAlongStep)
        << ",\"post_step_order\":" << manager->GetProcessOrdering(process,idxPostStep) << '}';
  }
  out << "]}";
}

int main(int argc, char** argv) {
  try {
    if (argc!=4) throw std::runtime_error("Usage: inventory config.py nominal|noop|dark-photon output.json");
    const std::string mode = argv[2];
    if (mode!="nominal" && mode!="noop" && mode!="dark-photon" && mode!="conflict")
      throw std::runtime_error("Unknown control mode");
    edmplugin::PluginManager::configure(edmplugin::standard::config());
    auto parameters = edm::cmspybind11::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("physics");
    const auto nominalFactory = parameters.getParameter<std::string>("type");
    if (nominalFactory != "SimG4Core/Physics/FTFP_BERT_EMM" && nominalFactory != "SimG4Core/Physics/FTFP_BERT_EMH")
      throw std::runtime_error("The control requires the explicit nominal EMM or EMH physics type");
    const std::string factory = nominalFactory + (mode=="dark-photon" || mode=="conflict" ? "_DarkPhoton" : "");
    parameters.addParameter<std::string>("type",factory);
    SimActivityRegistry registry;
    // Provide the normal Geant4 kernel and ownership, without initializing
    // geometry, building physics tables, or transporting any event.
    G4RunManager runManager;
    std::unique_ptr<PhysicsList> physicsOwner;
    if (mode=="noop") {
      physicsOwner = NoopExactNominal().make(parameters,registry);
    } else {
      std::unique_ptr<PhysicsListMakerBase> maker(PhysicsListFactory::get()->create(factory));
      physicsOwner = maker->make(parameters,registry);
    }
    auto* physics = physicsOwner.release();
    if (mode=="conflict") {
      const auto& point = parameters.getParameter<edm::ParameterSet>("DarkPhoton");
      const double lifetime = point.getParameter<double>("properCTauMm")*CLHEP::mm/CLHEP::c_light;
      new G4ParticleDefinition("dark_photon", .4*CLHEP::GeV, CLHEP::hbar_Planck/lifetime,
                               0.,2,-1,-1,0,0,0,"boson",0,0,32,false,lifetime,nullptr,false,"dark_photon",32);
    }
    runManager.SetUserInitialization(physics);
    physics->ConstructParticle();
    const auto initialCount = G4ParticleTable::GetParticleTable()->entries();
    physics->ConstructParticle();
    if (G4ParticleTable::GetParticleTable()->entries()!=initialCount)
      throw std::runtime_error("Repeated particle construction duplicated definitions");
    G4ParticleTable::GetParticleTable()->SetReadiness(true);
    // G4HepEm initializes its safety navigator during process construction.
    // Supply a synthetic empty world; never load CMS/LSS geometry or run events.
    auto* material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    auto* solid = new G4Box("ProcessInventoryWorld",CLHEP::m,CLHEP::m,CLHEP::m);
    auto* logical = new G4LogicalVolume(solid,material,"ProcessInventoryWorld");
    auto* world = new G4PVPlacement(nullptr,G4ThreeVector(),logical,"ProcessInventoryWorld",nullptr,false,0);
    runManager.DefineWorldVolume(world);
    // Construct initializes process managers and calls the nominal constructors;
    // no geometry, physics tables, tracking, detector or event processing is run.
    physics->Construct();
    std::ofstream out(argv[3]);
    if (!out) throw std::runtime_error("Cannot write process inventory");
    out << std::setprecision(17) << "{\"schema\":\"shift-dp-native-process-inventory-v1\",\"mode\":"
        << std::quoted(mode) << ",\"nominal_factory\":" << std::quoted(nominalFactory)
        << ",\"geant4_version_number\":" << G4VERSION_NUMBER << ",\"constructors\":[";
    for (int index=0; const auto* constructor=physics->GetPhysics(index); ++index) {
      if (index) out << ',';
      out << "{\"name\":" << std::quoted(std::string(constructor->GetPhysicsName()))
          << ",\"type\":" << constructor->GetPhysicsType() << '}';
    }
    out << "],\"sm_particles\":[";
    std::vector<const G4ParticleDefinition*> ordinaryParticles;
    auto* iterator = G4ParticleTable::GetParticleTable()->GetIterator();
    iterator->reset();
    while ((*iterator)()) {
      const auto* particle = iterator->value();
      if (particle->GetPDGEncoding() != 32) ordinaryParticles.push_back(particle);
    }
    std::sort(ordinaryParticles.begin(),ordinaryParticles.end(),[](const auto* left,const auto* right) {
      return left->GetParticleName() < right->GetParticleName();
    });
    bool first = true;
    for (const auto* particle : ordinaryParticles) {
      if (!first) out << ',';
      first=false;
      dumpParticle(out,*particle);
    }
    out << "],\"signal_decay_process_pdgs\":[";
    first = true;
    for (const auto* particle : ordinaryParticles) {
      const auto* manager = particle->GetProcessManager();
      if (manager == nullptr) continue;
      for (int index=0; index<manager->GetProcessListLength(); ++index) {
        if (dynamic_cast<const PredefinedSignalDecay*>((*manager->GetProcessList())[index]) == nullptr) continue;
        if (!first) out << ',';
        first = false;
        out << particle->GetPDGEncoding();
      }
    }
    out << "],\"dark_photon\":";
    const auto* dark = G4ParticleTable::GetParticleTable()->FindParticle(32);
    if (dark) dumpParticle(out,*dark); else out << "null";
    std::vector<std::string> libraries;
    dl_iterate_phdr([](dl_phdr_info* information, std::size_t, void* data) {
      if (information->dlpi_name[0] != '\0')
        static_cast<std::vector<std::string>*>(data)->emplace_back(information->dlpi_name);
      return 0;
    }, &libraries);
    std::sort(libraries.begin(),libraries.end());
    out << ",\"runtime_library_paths\":[";
    first = true;
    for (const auto& library : libraries) {
      if (!first) out << ',';
      first = false;
      out << std::quoted(library);
    }
    out << "],\"synthetic_vacuum_world\":true,\"detector_geometry_loaded\":false,\"events_transported\":0}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
