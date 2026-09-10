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
#include "PhysicsTools/ShiftLssGeometry/interface/BooleanFrameNormalizer.h"
#include "PhysicsTools/ShiftLssGeometry/interface/ExternalVolumeClipper.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
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

  std::array<double, 6> transformedNodeBounds(TGeoNode const& node,
                                              std::array<double, 9> const& rotation,
                                              std::array<double, 3> const& translation) {
    TGeoVolume const* volume = node.GetVolume();
    if (!volume || !volume->GetShape() || !node.GetMatrix()) {
      throw cms::Exception("UnsupportedGeometry") << "External LSS continuation node is malformed";
    }
    std::array<double, 3> low;
    std::array<double, 3> high;
    for (int axis = 0; axis < 3; ++axis) {
      volume->GetShape()->GetAxisRange(axis + 1, low[axis], high[axis]);
      if (!(low[axis] < high[axis])) {
        throw cms::Exception("UnsupportedGeometry") << "External LSS continuation has invalid local bounds";
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
      double local[3] = {corner & 1 ? high[0] : low[0], corner & 2 ? high[1] : low[1], corner & 4 ? high[2] : low[2]};
      double artifact[3];
      node.GetMatrix()->LocalToMaster(local, artifact);
      for (unsigned int row = 0; row < 3; ++row) {
        double global = translation[row];
        for (unsigned int column = 0; column < 3; ++column) {
          global += rotation[3 * row + column] * artifact[column];
        }
        bounds[2 * row] = std::min(bounds[2 * row], global);
        bounds[2 * row + 1] = std::max(bounds[2 * row + 1], global);
      }
    }
    return bounds;
  }

  std::array<double, 6> transformedDaughterBounds(TGeoVolume const& volume,
                                                  std::array<double, 9> const& rotation,
                                                  std::array<double, 3> const& translation) {
    if (volume.GetNdaughters() < 1) {
      throw cms::Exception("UnsupportedGeometry") << "External LSS world has no placed daughter volumes";
    }
    std::array<double, 6> bounds = {
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
    };
    for (int index = 0; index < volume.GetNdaughters(); ++index) {
      TGeoNode const* node = volume.GetNode(index);
      if (!node || !node->GetVolume() || !node->GetVolume()->GetShape() || !node->GetMatrix()) {
        throw cms::Exception("UnsupportedGeometry") << "External LSS world has an invalid placed daughter";
      }
      std::array<double, 3> low;
      std::array<double, 3> high;
      for (int axis = 0; axis < 3; ++axis) {
        node->GetVolume()->GetShape()->GetAxisRange(axis + 1, low[axis], high[axis]);
        if (!(low[axis] < high[axis])) {
          throw cms::Exception("UnsupportedGeometry")
              << "External LSS daughter " << node->GetName() << " has invalid local bounds";
        }
      }
      for (unsigned int corner = 0; corner < 8; ++corner) {
        double local[3] = {
            corner & 1 ? high[0] : low[0],
            corner & 2 ? high[1] : low[1],
            corner & 4 ? high[2] : low[2],
        };
        double artifact[3];
        node->GetMatrix()->LocalToMaster(local, artifact);
        for (unsigned int row = 0; row < 3; ++row) {
          double global = translation[row];
          for (unsigned int column = 0; column < 3; ++column) {
            global += rotation[3 * row + column] * artifact[column];
          }
          bounds[2 * row] = std::min(bounds[2 * row], global);
          bounds[2 * row + 1] = std::max(bounds[2 * row + 1], global);
        }
      }
    }
    return bounds;
  }

  TGeoVolume* findUniqueVolume(TGeoManager const& manager, std::string const& name) {
    TGeoVolume* result = nullptr;
    TObjArray const* volumes = manager.GetListOfVolumes();
    for (int index = 0; index < volumes->GetEntriesFast(); ++index) {
      auto* volume = dynamic_cast<TGeoVolume*>(volumes->UncheckedAt(index));
      if (!volume || name != volume->GetName()) {
        continue;
      }
      if (result && result != volume) {
        throw cms::Exception("GeometryVerification") << "CMS geometry contains multiple volumes named " << name;
      }
      result = volume;
    }
    if (!result) {
      throw cms::Exception("GeometryVerification") << "CMS geometry has no volume named " << name;
    }
    return result;
  }
}  // namespace

class ShiftLssGeometryESSource : public edm::ESProducer, public edm::EventSetupRecordIntervalFinder {
public:
  explicit ShiftLssGeometryESSource(edm::ParameterSet const& parameters)
      : gdmlFile_(parameters.getParameter<edm::FileInPath>("gdmlFile").fullPath()),
        geometryLabel_(parameters.getParameter<std::string>("geometryLabel")),
        detectorElementName_(parameters.getParameter<std::string>("detectorElementName")),
        externalMotherVolumeName_(parameters.getParameter<std::string>("externalMotherVolumeName")),
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
    if (externalMotherVolumeName_.empty()) {
      throw cms::Exception("Configuration") << "externalMotherVolumeName must not be empty";
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
    description.add<std::string>("externalMotherVolumeName", "cms:CMSE");
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
    TGeoVolume* externalMother = findUniqueVolume(manager, externalMotherVolumeName_);
    int const baselineWorldDaughters = cmsWorld->GetNdaughters();
    int const baselineMotherDaughters = externalMother->GetNdaughters();
    if (baselineWorldDaughters < 1) {
      throw cms::Exception("GeometryVerification") << "Standard CMSSW geometry payload has no world daughters";
    }
    std::vector<TGeoNode*> baselineWorldNodes;
    baselineWorldNodes.reserve(baselineWorldDaughters);
    for (int index = 0; index < baselineWorldDaughters; ++index) {
      baselineWorldNodes.push_back(cmsWorld->GetNode(index));
    }
    std::vector<TGeoNode*> baselineMotherNodes;
    baselineMotherNodes.reserve(baselineMotherDaughters);
    for (int index = 0; index < baselineMotherDaughters; ++index) {
      baselineMotherNodes.push_back(externalMother->GetNode(index));
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
    shift::BooleanFrameNormalizer normalizer;
    auto const importedFrame = normalizer.volume(importedVolume.ptr());
    if (!importedFrame.IsIdentity())
      throw cms::Exception("UnsupportedGeometry") << "External GDML world must have an untransformed solid frame";
    edm::LogInfo("ShiftLssGeometry") << "Normalized " << normalizer.changed()
                                     << " Boolean frames for faithful DDG4 conversion";
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
    auto const artifactBounds = transformedBounds(importedVolume, rotation_, translation);
    // The GDML world is a bookkeeping container and may be larger than its
    // physical daughters when it carries an exterior continuation shell. Raw
    // AABBs of legacy Boolean daughters are not safe for an aggregate gate.
    // The converter constructs this shell around the original bounded world;
    // check that physical shell and the original model box separately.
    TGeoNode const* rockNode = nullptr;
    for (int index = 0; index < importedVolume.ptr()->GetNdaughters(); ++index) {
      TGeoNode const* node = importedVolume.ptr()->GetNode(index);
      if (node && std::string(node->GetName()) == "shift_rock_continuation_pv") {
        if (rockNode) {
          throw cms::Exception("UnsupportedGeometry") << "External LSS has multiple rock continuation nodes";
        }
        rockNode = node;
      }
    }
    auto bounds = artifactBounds;
    if (rockNode) {
      auto const rockBounds = transformedNodeBounds(*rockNode, rotation_, translation);
      // The inner world is centered at the artifact origin. Its extent is
      // recovered from the shell subtraction's second operand by the
      // converter, so the continuation itself is the only new protected-zone
      // risk. The original artifact passed this gate before augmentation.
      bounds = rockBounds;
    }
    double const boundary = minimumAbsZCm_ * dd4hep::cm;
    if (!(bounds[4] >= boundary || bounds[5] <= -boundary)) {
      throw cms::Exception("UnsupportedGeometry")
          << "External LSS bounds [" << bounds[4] / dd4hep::cm << ", " << bounds[5] / dd4hep::cm
          << "] cm cross the protected |z| < " << minimumAbsZCm_ << " cm CMS region";
    }

    cmsWorld->RemoveNode(oldPlacement.ptr());
    shift::ExternalVolumeClipper clipper;
    for (auto const* node : baselineMotherNodes)
      clipper.protect(*node->GetVolume(), TGeoHMatrix(*node->GetMatrix()));
    TGeoHMatrix assemblyPlacement;
    assemblyPlacement.SetRotation(rotation_.data());
    assemblyPlacement.SetTranslation(translation.data());
    dd4hep::Assembly importedAssembly(detectorElementName_ + "_assembly");
    for (int index = 0; index < importedVolume.ptr()->GetNdaughters(); ++index) {
      TGeoNode* sourceNode = importedVolume.ptr()->GetNode(index);
      auto placement = assemblyPlacement;
      placement.Multiply(sourceNode->GetMatrix());
      importedAssembly.ptr()->AddNode(clipper.clip(sourceNode->GetVolume(), placement),
                                      sourceNode->GetNumber(),
                                      new TGeoHMatrix(*sourceNode->GetMatrix()));
    }
    // The protected CMS solids can themselves contain transformed Booleans.
    // Normalize the new difference expressions without mutating CMS volumes.
    shift::BooleanFrameNormalizer clippedNormalizer;
    clippedNormalizer.volume(importedAssembly.ptr());
    edm::LogInfo("ShiftLssGeometry") << "Applied " << clipper.subtractions()
                                     << " exact solid exclusions to preserve existing CMS volume ownership";
    dd4hep::Rotation3D modelRotation(rotation_.begin(), rotation_.end());
    dd4hep::Transform3D modelTransform(modelRotation, dd4hep::Position(translation[0], translation[1], translation[2]));
    dd4hep::PlacedVolume placed = dd4hep::Volume(externalMother).placeVolume(importedAssembly, 1, modelTransform);
    child.setPlacement(placed);

    if (cmsWorld->GetNdaughters() != baselineWorldDaughters) {
      throw cms::Exception("GeometryVerification")
          << "External attachment changed the CMS world daughter count: " << baselineWorldDaughters << " before, "
          << cmsWorld->GetNdaughters() << " after";
    }
    if (externalMother->GetNdaughters() != baselineMotherDaughters + 1) {
      throw cms::Exception("GeometryVerification")
          << "External attachment did not add exactly one assembly below " << externalMotherVolumeName_ << ": "
          << baselineMotherDaughters << " daughters before, " << externalMother->GetNdaughters() << " after";
    }
    for (int index = 0; index < baselineMotherDaughters; ++index) {
      if (externalMother->GetNode(index) != baselineMotherNodes[index]) {
        throw cms::Exception("GeometryVerification") << "External attachment replaced or reordered existing "
                                                     << externalMotherVolumeName_ << " daughter " << index;
      }
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
        cms::Exception error("GeometryOverlap");
        error << "External LSS geometry introduced " << finalOverlaps - baselineOverlaps << " overlap(s) at tolerance "
              << overlapToleranceCm_ << " cm";
        for (auto const* overlap : *manager.GetListOfOverlaps())
          error << "\n" << overlap->GetTitle();
        throw error;
      }
    }
    edm::LogInfo("ShiftLssGeometry")
        << "Preserved the standard CMSSW Extended geometry and attached the unwrapped external assembly below "
        << externalMotherVolumeName_ << " with transformed z bounds [" << bounds[4] / dd4hep::cm << ", "
        << bounds[5] / dd4hep::cm << "] cm (GDML container bounds [" << artifactBounds[4] / dd4hep::cm << ", "
        << artifactBounds[5] / dd4hep::cm << "] cm); all " << baselineWorldDaughters
        << " pre-existing CMS world daughter(s) and " << baselineMotherDaughters
        << " pre-existing mother-volume daughter(s) remain unchanged";
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
  std::string externalMotherVolumeName_;
  std::vector<double> artifactOriginInModelCm_;
  std::vector<double> modelOriginCm_;
  std::array<double, 9> rotation_;
  double minimumAbsZCm_;
  double overlapToleranceCm_;
  bool checkOverlaps_;
  edm::ESGetToken<FileBlob, GeometryFileRcd> geometryToken_;
};

DEFINE_FWK_EVENTSETUP_SOURCE(ShiftLssGeometryESSource);
