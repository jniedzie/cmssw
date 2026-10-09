#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Units.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

namespace {
  bool finite(HepMC::FourVector const& vector) {
    return std::isfinite(vector.x()) && std::isfinite(vector.y()) && std::isfinite(vector.z()) &&
           std::isfinite(vector.t());
  }

  std::array<long double, 4> components(HepMC::FourVector const& vector) {
    return {vector.px(), vector.py(), vector.pz(), vector.e()};
  }
}  // namespace

// Pythia owns the physical production, polarization and complete decay graph.
// Mark its completed neutral signal parents as generator-only (HepMC status 3).
// The standard CMSSW generator otherwise can try to transport a status-2 PDG 32
// with a predefined decay, although no Geant4 dark-photon definition is present.
// The CMS Pythia interface also retains native hard/recoil-copy statuses.
// Daughters retain their original displaced vertices; no decay is resampled.
class ShiftDarkPhotonHepMCProducer : public edm::global::EDProducer<> {
public:
  explicit ShiftDarkPhotonHepMCProducer(edm::ParameterSet const& config)
      : token_(consumes<edm::HepMCProduct>(config.getParameter<edm::InputTag>("src"))),
        massGeV_(config.getParameter<double>("massGeV")),
        generatedMassMinGeV_(config.getParameter<double>("generatedMassMinGeV")),
        generatedMassMaxGeV_(config.getParameter<double>("generatedMassMaxGeV")),
        momentumTolerance_(config.getParameter<double>("momentumRelativeTolerance")) {
    if (!(std::isfinite(massGeV_) && massGeV_ > 0. && std::isfinite(generatedMassMinGeV_) &&
          generatedMassMinGeV_ > 0. && std::isfinite(generatedMassMaxGeV_) &&
          generatedMassMaxGeV_ > generatedMassMinGeV_ && massGeV_ >= generatedMassMinGeV_ &&
          massGeV_ <= generatedMassMaxGeV_ && std::isfinite(momentumTolerance_) && momentumTolerance_ > 0. &&
          momentumTolerance_ < 1.))
      throw cms::Exception("Configuration")
          << "ShiftDarkPhotonHepMCProducer requires finite positive pole mass inside the explicit generated-mass "
             "interval and momentum tolerance in (0,1)";
    produces<edm::HepMCProduct>();
    produces<unsigned int>("adaptedParents");
  }

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const& input = event.get(token_);
    auto const* original = input.GetEvent();
    if (original == nullptr || original->vertices_empty())
      throw cms::Exception("EventCorruption") << "Dark-photon input HepMC has no vertices";
    // The common timing and simulation contracts use GeV and mm (including ct).
    // Refuse other units instead of silently changing any original coordinates.
    if (original->momentum_unit() != HepMC::Units::GEV || original->length_unit() != HepMC::Units::MM)
      throw cms::Exception("EventCorruption") << "Dark-photon input HepMC must use GeV/mm units";

    auto output = std::make_unique<HepMC::GenEvent>(*original);
    unsigned int signalParents = 0, physicalDecays = 0, adaptedParents = 0;
    for (auto particle = output->particles_begin(); particle != output->particles_end(); ++particle) {
      auto* parent = *particle;
      if (std::abs(parent->pdg_id()) != 32)
        continue;
      ++signalParents;
      auto const* production = parent->production_vertex();
      auto const* decay = parent->end_vertex();
      auto const& momentum = parent->momentum();
      double mass = parent->generated_mass();
      bool hasSignalContinuation = false;
      if (decay != nullptr) {
        for (auto child = decay->particles_out_const_begin(); child != decay->particles_out_const_end(); ++child)
          hasSignalContinuation = hasSignalContinuation || std::abs((*child)->pdg_id()) == 32;
      }
      // Pythia status 44 is an outgoing particle shifted by an ISR branching.
      // Accept it only as completed recoil-copy history, not as a new decay
      // prescription. Unknown statuses and unfinished parents remain rejected.
      bool knownStatus = parent->status() == 2 || parent->status() == 3 || parent->status() == 22 ||
                         parent->status() == 62 || (parent->status() == 44 && hasSignalContinuation);
      if (!knownStatus ||
          production == nullptr || decay == nullptr ||
          decay == production || decay->particles_in_size() != 1 || decay->particles_out_size() == 0)
        throw cms::Exception("EventCorruption")
            << "PDG 32 barcode " << parent->barcode() << " status " << parent->status()
            << " must have a complete single-parent decay graph with a recognized status";
      // A finite-width generator samples a physical invariant mass inside its
      // declared hard-process support. The pole mass is not a tolerance cut.
      // The caller must supply the same frozen support as the generator card.
      if (!finite(momentum) || momentum.e() <= 0. || !std::isfinite(mass) || mass < generatedMassMinGeV_ ||
          mass > generatedMassMaxGeV_ || !finite(production->position()) ||
          !finite(decay->position()))
        throw cms::Exception("EventCorruption")
            << "PDG 32 barcode " << parent->barcode() << " status " << parent->status() << " mass=" << mass
            << " GeV has invalid kinematics, mass or spacetime; pole=" << massGeV_ << " GeV, generated support=["
            << generatedMassMinGeV_ << "," << generatedMassMaxGeV_ << "] GeV; momentum=(" << momentum.px() << ","
            << momentum.py() << "," << momentum.pz() << "," << momentum.e() << "); production=("
            << production->position().x() << "," << production->position().y() << "," << production->position().z()
            << "," << production->position().t() << "); decay=(" << decay->position().x() << ","
            << decay->position().y() << "," << decay->position().z() << "," << decay->position().t() << ") mm";

      std::array<long double, 4> incoming{}, outgoing{};
      for (auto child = decay->particles_in_const_begin(); child != decay->particles_in_const_end(); ++child) {
        if (!finite((*child)->momentum()))
          throw cms::Exception("EventCorruption") << "Nonfinite momentum at a dark-photon decay";
        auto values = components((*child)->momentum());
        for (unsigned int component = 0; component < values.size(); ++component)
          incoming[component] += values[component];
      }
      for (auto child = decay->particles_out_const_begin(); child != decay->particles_out_const_end(); ++child) {
        if (!finite((*child)->momentum()))
          throw cms::Exception("EventCorruption") << "Nonfinite momentum at a dark-photon decay";
        auto values = components((*child)->momentum());
        for (unsigned int component = 0; component < values.size(); ++component)
          outgoing[component] += values[component];
      }
      // A Pythia hard-to-recoil copy records shower recoil, not an isolated
      // physical decay. Preserve this history even when its single copy edge
      // does not conserve momentum. Only the terminal decay must conserve it.
      if (!hasSignalContinuation) {
        ++physicalDecays;
        if (decay->particles_out_size() < 2)
          throw cms::Exception("EventCorruption") << "Dark photon must have at least two physical decay daughters";
        for (unsigned int component = 0; component < incoming.size(); ++component) {
          long double scale = std::max({1.L, std::abs(incoming[component]), std::abs(outgoing[component])});
          if (std::abs(incoming[component] - outgoing[component]) > momentumTolerance_ * scale)
            throw cms::Exception("EventCorruption")
                << "PDG 32 barcode " << parent->barcode() << " decay violates four-momentum conservation";
        }
      }
      if (parent->status() != 3) {
        parent->set_status(3);
        ++adaptedParents;
      }
    }
    if (signalParents == 0 || physicalDecays == 0)
      throw cms::Exception("EventCorruption") << "Dark-photon input contains no completed PDG 32 parent";
    event.put(std::make_unique<edm::HepMCProduct>(output.release()));
    event.put(std::make_unique<unsigned int>(adaptedParents), "adaptedParents");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription description;
    description.add<edm::InputTag>("src", edm::InputTag("VtxSmeared"));
    description.add<double>("massGeV");
    description.add<double>("generatedMassMinGeV");
    description.add<double>("generatedMassMaxGeV");
    description.add<double>("momentumRelativeTolerance", 1.e-6);
    descriptions.add("shiftDarkPhotonHepMC", description);
  }

private:
  edm::EDGetTokenT<edm::HepMCProduct> token_;
  double const massGeV_, generatedMassMinGeV_, generatedMassMaxGeV_, momentumTolerance_;
};

DEFINE_FWK_MODULE(ShiftDarkPhotonHepMCProducer);
