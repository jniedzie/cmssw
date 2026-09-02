#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <array>
#include <cmath>
#include <string>
#include <vector>

class ShiftLssMagneticFieldValidator : public edm::one::EDAnalyzer<> {
public:
  explicit ShiftLssMagneticFieldValidator(edm::ParameterSet const& parameters)
      : fieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>(
            edm::ESInputTag("", parameters.getParameter<std::string>("fieldLabel")))) {
    for (auto const& parameters : parameters.getParameter<std::vector<edm::ParameterSet>>("samples")) {
      auto const point = parameters.getParameter<std::vector<double>>("pointCm");
      auto const expected = parameters.getParameter<std::vector<double>>("expectedTesla");
      double const tolerance = parameters.getParameter<double>("toleranceTesla");
      if (point.size() != 3 || expected.size() != 3)
        throw cms::Exception("Configuration") << "Each field sample requires three point and field components";
      if (!(tolerance > 0.) || !std::isfinite(tolerance))
        throw cms::Exception("Configuration") << "Field-sample tolerance must be finite and positive";
      for (double value : point)
        if (!std::isfinite(value))
          throw cms::Exception("Configuration") << "Field-sample coordinates must be finite";
      for (double value : expected)
        if (!std::isfinite(value))
          throw cms::Exception("Configuration") << "Expected field components must be finite";
      samples_.push_back({parameters.getParameter<std::string>("name"),
                          {point[0], point[1], point[2]},
                          {expected[0], expected[1], expected[2]},
                          tolerance});
    }
    if (samples_.empty())
      throw cms::Exception("Configuration") << "At least one field sample is required";
  }

  void analyze(edm::Event const&, edm::EventSetup const& setup) override {
    auto const& field = setup.getData(fieldToken_);
    for (auto const& sample : samples_) {
      auto const actual = field.inTesla(GlobalPoint(sample.point[0], sample.point[1], sample.point[2]));
      std::array<double, 3> const components{actual.x(), actual.y(), actual.z()};
      for (unsigned int component = 0; component < components.size(); ++component) {
        if (!std::isfinite(components[component]) ||
            std::abs(components[component] - sample.expected[component]) > sample.tolerance) {
          throw cms::Exception("FieldValidation")
              << "SHIFT LSS field sample '" << sample.name << "' failed at (" << sample.point[0] << ", "
              << sample.point[1] << ", " << sample.point[2] << ") cm: expected (" << sample.expected[0] << ", "
              << sample.expected[1] << ", " << sample.expected[2] << ") T, got (" << components[0] << ", "
              << components[1] << ", " << components[2] << ") T with tolerance " << sample.tolerance << " T";
        }
      }
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription sample;
    sample.add<std::string>("name");
    sample.add<std::vector<double>>("pointCm");
    sample.add<std::vector<double>>("expectedTesla");
    sample.add<double>("toleranceTesla", 1.e-9);
    edm::ParameterSetDescription description;
    description.add<std::string>("fieldLabel", "");
    description.addVPSet("samples", sample, {});
    descriptions.addWithDefaultLabel(description);
  }

private:
  struct Sample {
    std::string name;
    std::array<double, 3> point;
    std::array<double, 3> expected;
    double tolerance;
  };

  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;
  std::vector<Sample> samples_;
};

DEFINE_FWK_MODULE(ShiftLssMagneticFieldValidator);
