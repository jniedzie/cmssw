#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>

namespace {
  constexpr std::array<char const*, 4> kParameterNames = {"dtSimHits", "cscSimHits", "rpcSimHits", "gemSimHits"};
  constexpr std::array<char const*, 4> kOutputInstances = {"MuonDTHits", "MuonCSCHits", "MuonRPCHits", "MuonGEMHits"};
}  // namespace

class ShiftSimHitTimeProducer : public edm::stream::EDProducer<> {
public:
  explicit ShiftSimHitTimeProducer(edm::ParameterSet const&);
  void produce(edm::Event&, edm::EventSetup const&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  std::array<edm::EDGetTokenT<edm::PSimHitContainer>, 4> tokens_;
  std::int32_t bxOffset_;
  double bunchSpacingNs_;
  double phaseNs_;
  double appliedShiftNs_;
  std::string modelVersion_;
};

ShiftSimHitTimeProducer::ShiftSimHitTimeProducer(edm::ParameterSet const& config)
    : bxOffset_(config.getParameter<std::int32_t>("bxOffset")),
      bunchSpacingNs_(config.getParameter<double>("bunchSpacingNs")),
      phaseNs_(config.getParameter<double>("phaseNs")),
      appliedShiftNs_(static_cast<double>(bxOffset_) * bunchSpacingNs_ + phaseNs_),
      modelVersion_(config.getParameter<std::string>("modelVersion")) {
  if (!(bunchSpacingNs_ > 0.0) || !std::isfinite(bunchSpacingNs_))
    throw cms::Exception("Configuration")
        << "ShiftSimHitTimeProducer: bunchSpacingNs must be finite and positive, got " << bunchSpacingNs_;
  if (!std::isfinite(phaseNs_) || !std::isfinite(appliedShiftNs_))
    throw cms::Exception("Configuration") << "ShiftSimHitTimeProducer: configured time shift is not finite";

  for (std::size_t index = 0; index < tokens_.size(); ++index) {
    tokens_[index] = consumes<edm::PSimHitContainer>(config.getParameter<edm::InputTag>(kParameterNames[index]));
    produces<edm::PSimHitContainer>(kOutputInstances[index]);
  }
  produces<double>("appliedShiftNs");
  produces<double>("phaseNs");
  produces<std::int32_t>("bxOffset");
  produces<std::string>("modelVersion");
}

void ShiftSimHitTimeProducer::produce(edm::Event& event, edm::EventSetup const&) {
  for (std::size_t index = 0; index < tokens_.size(); ++index) {
    auto const& input = event.get(tokens_[index]);
    auto output = std::make_unique<edm::PSimHitContainer>(input);
    for (auto& hit : *output)
      hit.setTof(hit.timeOfFlight() + appliedShiftNs_);
    event.put(std::move(output), kOutputInstances[index]);
  }
  event.put(std::make_unique<double>(appliedShiftNs_), "appliedShiftNs");
  event.put(std::make_unique<double>(phaseNs_), "phaseNs");
  event.put(std::make_unique<std::int32_t>(bxOffset_), "bxOffset");
  event.put(std::make_unique<std::string>(modelVersion_), "modelVersion");
}

void ShiftSimHitTimeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("dtSimHits", edm::InputTag("g4SimHits", "MuonDTHits"));
  description.add<edm::InputTag>("cscSimHits", edm::InputTag("g4SimHits", "MuonCSCHits"));
  description.add<edm::InputTag>("rpcSimHits", edm::InputTag("g4SimHits", "MuonRPCHits"));
  description.add<edm::InputTag>("gemSimHits", edm::InputTag("g4SimHits", "MuonGEMHits"));
  description.add<std::int32_t>("bxOffset", 0);
  description.add<double>("bunchSpacingNs", 25.0);
  description.add<double>("phaseNs", 0.0);
  description.add<std::string>("modelVersion", "same-simhit-reference-v1");
  descriptions.add("shiftSimHitTime", description);
}

DEFINE_FWK_MODULE(ShiftSimHitTimeProducer);
