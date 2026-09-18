#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <memory>
#include <vector>

class ShiftMuonTruthTrailFieldProducer : public edm::stream::EDProducer<> {
public:
  explicit ShiftMuonTruthTrailFieldProducer(edm::ParameterSet const& parameters)
      : xToken_(consumes<std::vector<float>>(parameters.getParameter<edm::InputTag>("x"))),
        yToken_(consumes<std::vector<float>>(parameters.getParameter<edm::InputTag>("y"))),
        zToken_(consumes<std::vector<float>>(parameters.getParameter<edm::InputTag>("z"))),
        fieldToken_(esConsumes(parameters.getParameter<edm::ESInputTag>("magneticField"))) {
    produces<std::vector<float>>("fieldX");
    produces<std::vector<float>>("fieldY");
    produces<std::vector<float>>("fieldZ");
  }

  void produce(edm::Event& event, edm::EventSetup const& setup) override {
    auto const& x = event.get(xToken_);
    auto const& y = event.get(yToken_);
    auto const& z = event.get(zToken_);
    if (x.size() != y.size() || x.size() != z.size())
      throw cms::Exception("LogicError") << "SHIFT muon truth-trail position vectors have different sizes";
    auto const& field = setup.getData(fieldToken_);
    auto fieldX = std::make_unique<std::vector<float>>();
    auto fieldY = std::make_unique<std::vector<float>>();
    auto fieldZ = std::make_unique<std::vector<float>>();
    fieldX->reserve(x.size());
    fieldY->reserve(x.size());
    fieldZ->reserve(x.size());
    for (std::size_t index = 0; index < x.size(); ++index) {
      auto const value = field.inTesla(GlobalPoint(x[index], y[index], z[index]));
      fieldX->push_back(value.x());
      fieldY->push_back(value.y());
      fieldZ->push_back(value.z());
    }
    event.put(std::move(fieldX), "fieldX");
    event.put(std::move(fieldY), "fieldY");
    event.put(std::move(fieldZ), "fieldZ");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription description;
    description.add<edm::InputTag>("x", edm::InputTag("g4SimHits", "shiftMuonTruthTrailX"));
    description.add<edm::InputTag>("y", edm::InputTag("g4SimHits", "shiftMuonTruthTrailY"));
    description.add<edm::InputTag>("z", edm::InputTag("g4SimHits", "shiftMuonTruthTrailZ"));
    description.add<edm::ESInputTag>("magneticField", edm::ESInputTag("", ""));
    descriptions.add("shiftMuonTruthTrailField", description);
  }

private:
  edm::EDGetTokenT<std::vector<float>> xToken_, yToken_, zToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;
};

DEFINE_FWK_MODULE(ShiftMuonTruthTrailFieldProducer);
