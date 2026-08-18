#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMC3Product.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>

namespace {
  constexpr double kSpeedOfLightMmPerNs = 299.792458;

  enum class TimingMode : std::int32_t { nominal = 0, legacy = 1, fixed = 2 };

  TimingMode parseMode(std::string const& mode) {
    if (mode == "nominal")
      return TimingMode::nominal;
    if (mode == "legacy")
      return TimingMode::legacy;
    if (mode == "fixed")
      return TimingMode::fixed;
    throw cms::Exception("Configuration")
        << "ShiftEventTimeProducer: unsupported timingMode '" << mode
        << "'; expected nominal, legacy, or fixed";
  }

  struct TimingResult {
    double sourceXmm;
    double sourceYmm;
    double sourceZmm;
    double sourceCtBeforeMm;
    double appliedShiftCtMm;
  };
}  // namespace

class ShiftEventTimeProducer : public edm::stream::EDProducer<> {
public:
  explicit ShiftEventTimeProducer(edm::ParameterSet const&);
  void produce(edm::Event&, edm::EventSetup const&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  double commonShiftCtMm(double sourceZmm) const;
  void putMetadata(edm::Event&, TimingResult const&) const;

  edm::EDGetTokenT<edm::HepMCProduct> srcToken_;
  edm::EDGetTokenT<edm::HepMC3Product> srcToken3_;
  std::string modeName_;
  TimingMode mode_;
  std::int32_t beamDirectionZ_;
  std::int32_t bxOffset_;
  double cmsReferenceZmm_;
  double bunchSpacingNs_;
  double phaseNs_;
  double fixedOffsetNs_;
  double legacyOffsetCtMm_;
  std::string modelVersion_;
};

ShiftEventTimeProducer::ShiftEventTimeProducer(edm::ParameterSet const& config)
    : srcToken_(consumes<edm::HepMCProduct>(config.getParameter<edm::InputTag>("src"))),
      srcToken3_(consumes<edm::HepMC3Product>(config.getParameter<edm::InputTag>("src"))),
      modeName_(config.getParameter<std::string>("timingMode")),
      mode_(parseMode(modeName_)),
      beamDirectionZ_(config.getParameter<std::int32_t>("beamDirectionZ")),
      bxOffset_(config.getParameter<std::int32_t>("bxOffset")),
      cmsReferenceZmm_(config.getParameter<double>("cmsReferenceZmm")),
      bunchSpacingNs_(config.getParameter<double>("bunchSpacingNs")),
      phaseNs_(config.getParameter<double>("phaseNs")),
      fixedOffsetNs_(config.getParameter<double>("fixedOffsetNs")),
      legacyOffsetCtMm_(config.getParameter<double>("legacyOffsetCtMm")),
      modelVersion_(config.getParameter<std::string>("modelVersion")) {
  if (beamDirectionZ_ != -1 && beamDirectionZ_ != 1) {
    throw cms::Exception("Configuration")
        << "ShiftEventTimeProducer: beamDirectionZ must be -1 or +1, got " << beamDirectionZ_;
  }
  if (!(bunchSpacingNs_ > 0.0)) {
    throw cms::Exception("Configuration")
        << "ShiftEventTimeProducer: bunchSpacingNs must be positive, got " << bunchSpacingNs_;
  }

  produces<edm::HepMCProduct>();
  produces<edm::HepMC3Product>();
  produces<double>("sourceXmm");
  produces<double>("sourceYmm");
  produces<double>("sourceZmm");
  produces<double>("sourceCtBeforeMm");
  produces<double>("appliedShiftCtMm");
  produces<double>("appliedShiftNs");
  produces<double>("phaseNs");
  produces<std::int32_t>("beamDirectionZ");
  produces<std::int32_t>("bxOffset");
  produces<std::int32_t>("timingModeCode");
  produces<std::string>("timingMode");
  produces<std::string>("modelVersion");
}

double ShiftEventTimeProducer::commonShiftCtMm(double sourceZmm) const {
  double shiftCtMm = 0.0;
  switch (mode_) {
    case TimingMode::nominal:
      // The beam travels in beamDirectionZ.  Relative to a CMS crossing at
      // cmsReferenceZmm, a source upstream along that direction occurred at
      // ct = beamDirectionZ * (sourceZ - cmsReferenceZ).
      shiftCtMm = beamDirectionZ_ * (sourceZmm - cmsReferenceZmm_);
      break;
    case TimingMode::legacy:
      // Regression mode for the former Pythia Beams:offsetTime value.  It is
      // an additive common shift, exactly as that Pythia setting was.
      shiftCtMm = legacyOffsetCtMm_;
      break;
    case TimingMode::fixed:
      shiftCtMm = fixedOffsetNs_ * kSpeedOfLightMmPerNs;
      break;
  }
  return shiftCtMm + (static_cast<double>(bxOffset_) * bunchSpacingNs_ + phaseNs_) * kSpeedOfLightMmPerNs;
}

void ShiftEventTimeProducer::putMetadata(edm::Event& event, TimingResult const& result) const {
  event.put(std::make_unique<double>(result.sourceXmm), "sourceXmm");
  event.put(std::make_unique<double>(result.sourceYmm), "sourceYmm");
  event.put(std::make_unique<double>(result.sourceZmm), "sourceZmm");
  event.put(std::make_unique<double>(result.sourceCtBeforeMm), "sourceCtBeforeMm");
  event.put(std::make_unique<double>(result.appliedShiftCtMm), "appliedShiftCtMm");
  event.put(std::make_unique<double>(result.appliedShiftCtMm / kSpeedOfLightMmPerNs), "appliedShiftNs");
  event.put(std::make_unique<double>(phaseNs_), "phaseNs");
  event.put(std::make_unique<std::int32_t>(beamDirectionZ_), "beamDirectionZ");
  event.put(std::make_unique<std::int32_t>(bxOffset_), "bxOffset");
  event.put(std::make_unique<std::int32_t>(static_cast<std::int32_t>(mode_)), "timingModeCode");
  event.put(std::make_unique<std::string>(modeName_), "timingMode");
  event.put(std::make_unique<std::string>(modelVersion_), "modelVersion");
}

void ShiftEventTimeProducer::produce(edm::Event& event, edm::EventSetup const&) {
  edm::Handle<edm::HepMCProduct> input;
  if (event.getByToken(srcToken_, input)) {
    HepMC::GenEvent const* inputEvent = input->GetEvent();
    if (inputEvent == nullptr || inputEvent->vertices_empty())
      throw cms::Exception("EventCorruption") << "ShiftEventTimeProducer: input HepMC event has no vertices";

    HepMC::GenVertex const* sourceVertex = inputEvent->signal_process_vertex();
    if (sourceVertex == nullptr)
      sourceVertex = *inputEvent->vertices_begin();
    auto const& source = sourceVertex->position();
    TimingResult result{source.x(), source.y(), source.z(), source.t(), commonShiftCtMm(source.z())};

    auto shiftedEvent = std::make_unique<HepMC::GenEvent>(*inputEvent);
    for (auto vertex = shiftedEvent->vertices_begin(); vertex != shiftedEvent->vertices_end(); ++vertex) {
      auto const& position = (*vertex)->position();
      (*vertex)->set_position(
          HepMC::FourVector(position.x(), position.y(), position.z(), position.t() + result.appliedShiftCtMm));
    }
    event.put(std::make_unique<edm::HepMCProduct>(shiftedEvent.release()));
    putMetadata(event, result);
    edm::LogVerbatim("ShiftEventTiming")
        << "mode=" << modeName_ << " modelVersion=" << modelVersion_ << " source=(" << result.sourceXmm << ","
        << result.sourceYmm << "," << result.sourceZmm << "," << result.sourceCtBeforeMm
        << ") mm appliedShiftCt=" << result.appliedShiftCtMm
        << " mm appliedShift=" << result.appliedShiftCtMm / kSpeedOfLightMmPerNs << " ns bxOffset=" << bxOffset_
        << " phase=" << phaseNs_ << " ns";
    return;
  }

  edm::Handle<edm::HepMC3Product> input3;
  if (!event.getByToken(srcToken3_, input3))
    throw cms::Exception("ProductNotFound") << "ShiftEventTimeProducer: neither HepMCProduct nor HepMC3Product found";

  HepMC3::GenEvent shiftedEvent;
  shiftedEvent.read_data(*input3->GetEvent());
  if (shiftedEvent.vertices().empty())
    throw cms::Exception("EventCorruption") << "ShiftEventTimeProducer: input HepMC3 event has no vertices";
  auto const sourceVertex = shiftedEvent.vertices().front();
  auto const& source = sourceVertex->position();
  TimingResult result{source.x(), source.y(), source.z(), source.t(), commonShiftCtMm(source.z())};
  shiftedEvent.shift_position_by(HepMC3::FourVector(0.0, 0.0, 0.0, result.appliedShiftCtMm));
  event.put(std::make_unique<edm::HepMC3Product>(shiftedEvent));
  putMetadata(event, result);
}

void ShiftEventTimeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("src", edm::InputTag("VtxSmeared"));
  description.add<std::string>("timingMode", "nominal");
  description.add<std::int32_t>("beamDirectionZ", -1);
  description.add<std::int32_t>("bxOffset", 0);
  description.add<double>("cmsReferenceZmm", 0.0);
  description.add<double>("bunchSpacingNs", 25.0);
  description.add<double>("phaseNs", 0.0);
  description.add<double>("fixedOffsetNs", 0.0);
  description.add<double>("legacyOffsetCtMm", -148000.0);
  description.add<std::string>("modelVersion", "fixed-target-v1");
  descriptions.add("shiftEventTime", description);
}

DEFINE_FWK_MODULE(ShiftEventTimeProducer);
