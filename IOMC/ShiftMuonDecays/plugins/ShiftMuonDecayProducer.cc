#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/AbstractServices/interface/RandomNumberGenerator.h"
#include "GeneratorInterface/Pythia8Interface/interface/P8RndmEngine.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

#include <cmath>
#include <cstdlib>
#include <memory>
#include <set>
#include <sstream>
#include <vector>

// Complete undecayed muon-producing hadron chains in a persisted GEN event.
// The hard event, existing decays, event identity and native weight stay intact.
// Channels retain their physical branching ratios: no muon-enrichment filter.
class ShiftMuonDecayProducer : public edm::stream::EDProducer<> {
public:
  explicit ShiftMuonDecayProducer(edm::ParameterSet const& config)
      : token_(consumes<edm::HepMCProduct>(config.getParameter<edm::InputTag>("src"))),
        pythia_(std::getenv("PYTHIA8DATA") ? std::getenv("PYTHIA8DATA") : "", false),
        random_(std::make_shared<gen::P8RndmEngine>()) {
    auto xy = config.getParameter<double>("decayCylinderRadiusMm");
    auto z = config.getParameter<double>("decayCylinderHalfLengthMm");
    if (!(std::isfinite(xy) && xy > 0 && std::isfinite(z) && z > 0))
      throw cms::Exception("Configuration") << "Finite positive decay-cylinder dimensions required";
    pythia_.readString("ProcessLevel:all = off");
    pythia_.readString("PartonLevel:all = off");
    pythia_.readString("HadronLevel:Hadronize = off");
    pythia_.readString("ParticleDecays:limitTau0 = off");
    pythia_.readString("ParticleDecays:limitCylinder = on");
    pythia_.settings.parm("ParticleDecays:xyMax", xy);
    pythia_.settings.parm("ParticleDecays:zMax", z);
    pythia_.readString("Init:showProcesses = off");
    pythia_.readString("Init:showChangedSettings = off");
    pythia_.readString("Init:showChangedParticleData = off");
    pythia_.setRndmEnginePtr(random_);

    // Transitive closure over all enabled, positive-BR Pythia decay channels.
    // This includes indirect chains (strange baryons -> pions -> muons),
    // rather than assuming that three named species exhaust the parents.
    std::set<int> reachesMuon{13};
    bool changed = true;
    while (changed) {
      changed = false;
      for (auto const& item : pythia_.particleData) {
        int id = item.first;
        auto const& entry = item.second;
        // canDecay is initialized by Pythia::init(), after configuration.
        // The XML channel table is already available at this point.
        if (reachesMuon.count(id) || entry->sizeChannels() == 0)
          continue;
        for (int c = 0; c < entry->sizeChannels(); ++c) {
          auto const& channel = entry->channel(c);
          if (!channel.onMode() || !(channel.bRatio() > 0))
            continue;
          for (int d = 0; d < channel.multiplicity(); ++d) {
            if (reachesMuon.count(std::abs(channel.product(d)))) {
              reachesMuon.insert(id);
              changed = true;
              break;
            }
          }
          if (reachesMuon.count(id))
            break;
        }
      }
    }
    for (int id : reachesMuon) {
      if (id != 13) {
        enabled_.insert(id);
        pythia_.particleData.mayDecay(id, true);
      }
    }
    // Muons remain available to detector transport; do not decay the observable.
    pythia_.particleData.mayDecay(13, false);
    converter_.set_store_pdf(false);
    converter_.set_store_proc(false);
    converter_.set_store_xsec(false);
    converter_.set_store_weights(false);
    produces<edm::HepMCProduct>();
    produces<unsigned int>("decayedParents");
    produces<unsigned int>("addedMuons");
    produces<std::vector<int>>("enabledPdgIds");
    std::ostringstream ids;
    for (int id : enabled_)
      if (pythia_.particleData.tau0(id) >= 10.)
        ids << ' ' << id;
    edm::LogVerbatim("ShiftMuonDecays") << "Muon-producing parents with ctau >= 10 mm:" << ids.str();
  }

  void produce(edm::Event& event, edm::EventSetup const&) override {
    edm::Service<edm::RandomNumberGenerator> rng;
    if (!rng.isAvailable())
      throw cms::Exception("Configuration") << "ShiftMuonDecayProducer needs RandomNumberGeneratorService";
    random_->setRandomEngine(&rng->getEngine(event.streamID()));
    if (!initialized_) {
      if (!pythia_.init())
        throw cms::Exception("PythiaError") << "Could not initialize residual decays";
      initialized_ = true;
    }
    auto const& input = event.get(token_);
    auto output = std::make_unique<HepMC::GenEvent>(*input.GetEvent());
    if (output->momentum_unit() != HepMC::Units::GEV || output->length_unit() != HepMC::Units::MM)
      throw cms::Exception("EventCorruption") << "Residual decays require HepMC GeV/mm units";
    std::vector<HepMC::GenParticle*> parents;
    int maximumBarcode = 0;
    for (auto it = output->particles_begin(); it != output->particles_end(); ++it) {
      auto* particle = *it;
      maximumBarcode = std::max(maximumBarcode, particle->barcode());
      if (particle->status() == 1 && !particle->end_vertex() && enabled_.count(std::abs(particle->pdg_id())))
        parents.push_back(particle);
    }
    unsigned int decayed = 0, muons = 0;
    for (auto* parent : parents) {
      if (!parent->production_vertex())
        throw cms::Exception("EventCorruption") << "Missing parent production vertex";
      auto const& p = parent->momentum();
      auto const& v = parent->production_vertex()->position();
      pythia_.event.reset();
      int index = pythia_.event.append(parent->pdg_id(), 93, 0, 0, p.px(), p.py(), p.pz(), p.e(),
                                     parent->generated_mass());
      pythia_.event[index].vProd(v.x(), v.y(), v.z(), v.t());
      pythia_.event[index].tau(pythia_.rndm.exp() * pythia_.particleData.tau0(parent->pdg_id()));
      if (!pythia_.moreDecays())
        throw cms::Exception("PythiaError") << "Residual decay failed for PDG " << parent->pdg_id();
      if (pythia_.event[index].isFinal())
        continue;  // Sampled decay lies outside the established flight cylinder.
      parent->set_status(2);
      for (int i = 2; i < pythia_.event.size(); ++i)
        if (pythia_.event[i].isFinal() && std::abs(pythia_.event[i].id()) == 13)
          ++muons;
      if (!converter_.fill_next_event(pythia_.event, output.get(), -1, nullptr, nullptr, true, parent,
                                      maximumBarcode))
        throw cms::Exception("EventCorruption") << "Could not append residual decay to HepMC";
      maximumBarcode += pythia_.event.size() - 2;
      ++decayed;
    }
    event.put(std::make_unique<edm::HepMCProduct>(output.release()));
    event.put(std::make_unique<unsigned int>(decayed), "decayedParents");
    event.put(std::make_unique<unsigned int>(muons), "addedMuons");
    event.put(std::make_unique<std::vector<int>>(enabled_.begin(), enabled_.end()), "enabledPdgIds");
  }

private:
  edm::EDGetTokenT<edm::HepMCProduct> token_;
  Pythia8::Pythia pythia_;
  gen::P8RndmEnginePtr random_;
  HepMC::Pythia8ToHepMC converter_;
  std::set<int> enabled_;
  bool initialized_ = false;
};

DEFINE_FWK_MODULE(ShiftMuonDecayProducer);
