#include <string>

#include "DD4hep/Detector.h"
#include "DetectorDescription/DDCMS/interface/DDDetector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "TGeoVolume.h"

class ShiftLssGeometryVerifier : public edm::one::EDAnalyzer<> {
public:
  explicit ShiftLssGeometryVerifier(edm::ParameterSet const& parameters)
      : minimumWorldDaughters_(parameters.getParameter<unsigned int>("minimumWorldDaughters")),
        externalElementPath_(parameters.getParameter<std::string>("externalElementPath")),
        detectorToken_(esConsumes<cms::DDDetector, IdealGeometryRecord>(
            edm::ESInputTag("", parameters.getParameter<std::string>("productLabel")))) {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription description;
    description.add<std::string>("productLabel", "");
    description.add<unsigned int>("minimumWorldDaughters", 2);
    description.add<std::string>("externalElementPath", "/world/shiftLssExternal");
    descriptions.add("shiftLssGeometryVerifier", description);
  }

  void analyze(edm::Event const&, edm::EventSetup const& setup) override {
    cms::DDDetector const& detector = setup.getData(detectorToken_);
    int const worldDaughters = detector.worldVolume().ptr()->GetNdaughters();
    if (worldDaughters < static_cast<int>(minimumWorldDaughters_)) {
      throw cms::Exception("GeometryVerification")
          << "Extended geometry has only " << worldDaughters << " world daughters; expected at least "
          << minimumWorldDaughters_ << " (standard CMS geometry plus external LSS)";
    }
    dd4hep::DetElement const external = detector.findElement(externalElementPath_);
    if (!external.isValid() || !external.placement().isValid()) {
      throw cms::Exception("GeometryVerification")
          << "External LSS detector element is absent: " << externalElementPath_;
    }
    edm::LogInfo("ShiftLssGeometry") << "Verified " << worldDaughters << " world daughters and external element "
                                     << externalElementPath_ << " in the Extended geometry";
  }

private:
  unsigned int minimumWorldDaughters_;
  std::string externalElementPath_;
  edm::ESGetToken<cms::DDDetector, IdealGeometryRecord> detectorToken_;
};

DEFINE_FWK_MODULE(ShiftLssGeometryVerifier);
