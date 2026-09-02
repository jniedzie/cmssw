#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "CondFormats/Common/interface/FileBlob.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"
#include "DetectorDescription/DDCMS/interface/DDDetector.h"
#include "FWCore/Concurrency/interface/SharedResourceNames.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/FileInPath.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "TGeoManager.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TObjArray.h"

namespace {
  struct GlobalGeoManagerGuard {
    explicit GlobalGeoManagerGuard(TGeoManager* manager) : previous_(gGeoManager) { gGeoManager = manager; }
    ~GlobalGeoManagerGuard() { gGeoManager = previous_; }
    TGeoManager* previous_;
  };

  std::array<double, 9> checkedRotation(std::vector<double> const& values) {
    if (values.size() != 9) {
      throw cms::Exception("Configuration") << "modelToCms must contain nine row-major values";
    }
    std::array<double, 9> rotation;
    std::copy(values.begin(), values.end(), rotation.begin());
    for (unsigned int row = 0; row < 3; ++row) {
      for (unsigned int other = 0; other < 3; ++other) {
        double product = 0.0;
        for (unsigned int column = 0; column < 3; ++column) {
          product += rotation[3 * row + column] * rotation[3 * other + column];
        }
        double const expected = row == other ? 1.0 : 0.0;
        if (std::abs(product - expected) > 1.e-9) {
          throw cms::Exception("Configuration") << "modelToCms must be orthonormal";
        }
      }
    }
    double const determinant = rotation[0] * (rotation[4] * rotation[8] - rotation[5] * rotation[7]) -
                               rotation[1] * (rotation[3] * rotation[8] - rotation[5] * rotation[6]) +
                               rotation[2] * (rotation[3] * rotation[7] - rotation[4] * rotation[6]);
    if (std::abs(determinant - 1.0) > 1.e-9) {
      throw cms::Exception("Configuration") << "modelToCms must be a proper rotation with determinant +1";
    }
    return rotation;
  }

  void preflightGdml(std::string const& path) {
    std::ifstream input(path);
    if (!input) {
      throw cms::Exception("FileOpenError") << "Cannot open external LSS GDML " << path;
    }
    std::string line;
    while (std::getline(input, line)) {
      if (line.find("<multiUnion") != std::string::npos) {
        throw cms::Exception("UnsupportedGeometry")
            << "External LSS GDML contains multiUnion, unsupported by CMSSW ROOT 6.36 TGDMLParse; "
               "provide a validated lossless ROOT-compatible conversion";
      }
    }
  }

  std::array<double, 6> transformedBounds(TGeoVolume const& volume,
                                          std::array<double, 9> const& rotation,
                                          std::array<double, 3> const& translation) {
    std::array<double, 3> low;
    std::array<double, 3> high;
    for (int axis = 0; axis < 3; ++axis) {
      volume.GetShape()->GetAxisRange(axis + 1, low[axis], high[axis]);
      if (!(low[axis] < high[axis])) {
        throw cms::Exception("UnsupportedGeometry") << "External LSS top volume has invalid local bounds";
      }
    }
    std::array<double, 6> bounds = {
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
    };
    for (unsigned int corner = 0; corner < 8; ++corner) {
      std::array<double, 3> local = {
          corner & 1 ? high[0] : low[0],
          corner & 2 ? high[1] : low[1],
          corner & 4 ? high[2] : low[2],
      };
      for (unsigned int row = 0; row < 3; ++row) {
        double global = translation[row];
        for (unsigned int column = 0; column < 3; ++column) {
          global += rotation[3 * row + column] * local[column];
        }
        bounds[2 * row] = std::min(bounds[2 * row], global);
        bounds[2 * row + 1] = std::max(bounds[2 * row + 1], global);
      }
    }
    return bounds;
  }
}  // namespace

class ShiftLssGeometryESSource : public edm::ESProducer, public edm::EventSetupRecordIntervalFinder {
public:
  explicit ShiftLssGeometryESSource(edm::ParameterSet const& parameters)
      : gdmlFile_(parameters.getParameter<edm::FileInPath>("gdmlFile").fullPath()),
        geometryLabel_(parameters.getParameter<std::string>("geometryLabel")),
        detectorElementName_(parameters.getParameter<std::string>("detectorElementName")),
        artifactOriginInModelCm_(parameters.getParameter<std::vector<double>>("artifactOriginInModelCm")),
        modelOriginCm_(parameters.getParameter<std::vector<double>>("modelOriginCm")),
        rotation_(checkedRotation(parameters.getParameter<std::vector<double>>("modelToCms"))),
        minimumAbsZCm_(parameters.getParameter<double>("minimumAbsZCm")),
        overlapToleranceCm_(parameters.getParameter<double>("overlapToleranceCm")),
        checkOverlaps_(parameters.getParameter<bool>("checkOverlaps")) {
    if (artifactOriginInModelCm_.size() != 3 || modelOriginCm_.size() != 3) {
      throw cms::Exception("Configuration") << "artifactOriginInModelCm and modelOriginCm must contain three values";
    }
    if (!(minimumAbsZCm_ > 0.0) || !(overlapToleranceCm_ > 0.0)) {
      throw cms::Exception("Configuration") << "minimumAbsZCm and overlapToleranceCm must be positive";
    }
    if (detectorElementName_.empty() || detectorElementName_ == "world") {
      throw cms::Exception("Configuration") << "detectorElementName must be a non-world child name";
    }
    usesResources({edm::ESSharedResourceNames::kDD4hep});
    auto collector = setWhatProduced(this, &ShiftLssGeometryESSource::produce);
    geometryToken_ = collector.consumes(edm::ESInputTag("", geometryLabel_));
    findingRecord<IdealGeometryRecord>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription description;
    description.add<edm::FileInPath>("gdmlFile");
    description.add<std::string>("geometryLabel", "Extended");
    description.add<std::string>("detectorElementName", "shiftLssExternal");
    description.add<std::vector<double>>("artifactOriginInModelCm");
    description.add<std::vector<double>>("modelOriginCm");
    description.add<std::vector<double>>("modelToCms");
    description.add<double>("minimumAbsZCm");
    description.add<double>("overlapToleranceCm", 0.001);
    description.add<bool>("checkOverlaps", true);
    descriptions.add("shiftLssGeometryESSource", description);
  }

private:
  std::unique_ptr<cms::DDDetector> produce(IdealGeometryRecord const& record) {
    preflightGdml(gdmlFile_);
    edm::ESTransientHandle<FileBlob> geometryBlob = record.getTransientHandle(geometryToken_);
    std::unique_ptr<std::vector<unsigned char>> payload = geometryBlob->getUncompressedBlob();
    auto detector =
        std::make_unique<cms::DDDetector>(geometryLabel_, std::string(payload->begin(), payload->end()), true);

    dd4hep::Detector const* description = detector->description();
    TGeoManager& manager = const_cast<TGeoManager&>(detector->manager());
    GlobalGeoManagerGuard managerGuard(&manager);
    TGeoVolume* cmsWorld = detector->worldVolume().ptr();
    int const baselineWorldDaughters = cmsWorld->GetNdaughters();
    if (baselineWorldDaughters < 1) {
      throw cms::Exception("GeometryVerification") << "Standard CMSSW geometry payload has no world daughters";
    }
    std::vector<TGeoNode*> baselineWorldNodes;
    baselineWorldNodes.reserve(baselineWorldDaughters);
    for (int index = 0; index < baselineWorldDaughters; ++index) {
      baselineWorldNodes.push_back(cmsWorld->GetNode(index));
    }
    int baselineOverlaps = 0;
    if (checkOverlaps_) {
      manager.CheckOverlaps(overlapToleranceCm_ * dd4hep::cm, "s");
      baselineOverlaps = manager.GetListOfOverlaps()->GetEntries();
    }
    std::string const path = "/world/" + detectorElementName_;
    char const* arguments[] = {"-input", gdmlFile_.c_str(), "-path", path.c_str(), nullptr};
    description->apply("DD4hep_ROOTGDMLParse", 4, const_cast<char**>(arguments));

    dd4hep::DetElement child = description->world().child(detectorElementName_, false);
    if (!child.isValid() || !child.placement().isValid()) {
      throw cms::Exception("UnsupportedGeometry") << "GDML importer did not create " << path;
    }
    dd4hep::PlacedVolume oldPlacement = child.placement();
    dd4hep::Volume importedVolume = oldPlacement.volume();
    // The bounded converter recentres the source model around an artifact
    // origin.  Place that artifact origin at the transformed source-model
    // coordinate; fields use the same modelOrigin + R * modelPoint contract.
    std::array<double, 3> translation;
    for (unsigned int row = 0; row < 3; ++row) {
      translation[row] = modelOriginCm_[row];
      for (unsigned int column = 0; column < 3; ++column) {
        translation[row] += rotation_[3 * row + column] * artifactOriginInModelCm_[column];
      }
      translation[row] *= dd4hep::cm;
    }
    auto const bounds = transformedBounds(importedVolume, rotation_, translation);
    double const boundary = minimumAbsZCm_ * dd4hep::cm;
    if (!(bounds[4] >= boundary || bounds[5] <= -boundary)) {
      throw cms::Exception("UnsupportedGeometry")
          << "External LSS bounds [" << bounds[4] / dd4hep::cm << ", " << bounds[5] / dd4hep::cm
          << "] cm cross the protected |z| < " << minimumAbsZCm_ << " cm CMS region";
    }

    dd4hep::Volume mother = description->worldVolume();
    mother.ptr()->RemoveNode(oldPlacement.ptr());
    dd4hep::Rotation3D modelRotation(rotation_.begin(), rotation_.end());
    dd4hep::Transform3D modelTransform(
        modelRotation, dd4hep::Position(translation[0], translation[1], translation[2]));
    dd4hep::PlacedVolume placed = mother.placeVolume(importedVolume, 1, modelTransform);
    child.setPlacement(placed);

    if (cmsWorld->GetNdaughters() != baselineWorldDaughters + 1) {
      throw cms::Exception("GeometryVerification")
          << "External attachment changed the CMS world daughter count unexpectedly: " << baselineWorldDaughters
          << " before, " << cmsWorld->GetNdaughters() << " after";
    }
    for (int index = 0; index < baselineWorldDaughters; ++index) {
      if (cmsWorld->GetNode(index) != baselineWorldNodes[index]) {
        throw cms::Exception("GeometryVerification")
            << "External attachment replaced or reordered standard CMS world daughter " << index;
      }
    }

    if (checkOverlaps_) {
      manager.CheckOverlaps(overlapToleranceCm_ * dd4hep::cm, "s");
      int const finalOverlaps = manager.GetListOfOverlaps()->GetEntries();
      if (finalOverlaps > baselineOverlaps) {
        throw cms::Exception("GeometryOverlap")
            << "External LSS geometry introduced " << finalOverlaps - baselineOverlaps << " overlap(s) at tolerance "
            << overlapToleranceCm_ << " cm";
      }
    }
    edm::LogInfo("ShiftLssGeometry") << "Preserved the standard CMSSW Extended geometry and attached external element "
                                     << detectorElementName_ << " with transformed z bounds [" << bounds[4] / dd4hep::cm
                                     << ", " << bounds[5] / dd4hep::cm << "] cm; all " << baselineWorldDaughters
                                     << " pre-existing CMS world daughter(s) remain unchanged";
    return detector;
  }

  void setIntervalFor(edm::eventsetup::EventSetupRecordKey const&,
                      edm::IOVSyncValue const&,
                      edm::ValidityInterval& interval) override {
    interval = edm::ValidityInterval(edm::IOVSyncValue::beginOfTime(), edm::IOVSyncValue::endOfTime());
  }

  std::string gdmlFile_;
  std::string geometryLabel_;
  std::string detectorElementName_;
  std::vector<double> artifactOriginInModelCm_;
  std::vector<double> modelOriginCm_;
  std::array<double, 9> rotation_;
  double minimumAbsZCm_;
  double overlapToleranceCm_;
  bool checkOverlaps_;
  edm::ESGetToken<FileBlob, GeometryFileRcd> geometryToken_;
};

DEFINE_FWK_EVENTSETUP_SOURCE(ShiftLssGeometryESSource);
