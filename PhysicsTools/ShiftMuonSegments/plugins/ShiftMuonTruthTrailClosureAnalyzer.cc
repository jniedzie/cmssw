#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {
  struct Checkpoint {
    int kind;
    double x, y, z, px, py, pz, length, x0, energyLoss;
    std::string material, volume;
  };

  Surface::RotationType normalTo(GlobalVector const& direction) {
    GlobalVector const zAxis = direction.unit();
    GlobalVector yAxis(zAxis.y(), -zAxis.x(), 0.);
    if (!(yAxis.mag2() > 1.e-12))
      yAxis = GlobalVector(0., 1., 0.);
    yAxis = yAxis.unit();
    GlobalVector const xAxis = yAxis.cross(zAxis).unit();
    return Surface::RotationType(xAxis, yAxis, zAxis);
  }

  double angle(GlobalVector const& first, GlobalVector const& second) {
    double const cosine = std::clamp(static_cast<double>(first.unit().dot(second.unit())), -1., 1.);
    return std::acos(cosine);
  }
}  // namespace

class ShiftMuonTruthTrailClosureAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit ShiftMuonTruthTrailClosureAnalyzer(edm::ParameterSet const& parameters)
      : simTracks_(consumes<edm::SimTrackContainer>(parameters.getParameter<edm::InputTag>("simTracks"))),
        simVertices_(consumes<edm::SimVertexContainer>(parameters.getParameter<edm::InputTag>("simVertices"))),
        cscSimHits_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("cscSimHits"))),
        trailTrackId_(consumes<std::vector<int>>(parameters.getParameter<edm::InputTag>("trailTrackId"))),
        trailKind_(consumes<std::vector<int>>(parameters.getParameter<edm::InputTag>("trailKind"))),
        trailX_(consumes<std::vector<double>>(parameters.getParameter<edm::InputTag>("trailX"))),
        trailY_(consumes<std::vector<double>>(parameters.getParameter<edm::InputTag>("trailY"))),
        trailZ_(consumes<std::vector<double>>(parameters.getParameter<edm::InputTag>("trailZ"))),
        trailPx_(consumes<std::vector<double>>(parameters.getParameter<edm::InputTag>("trailPx"))),
        trailPy_(consumes<std::vector<double>>(parameters.getParameter<edm::InputTag>("trailPy"))),
        trailPz_(consumes<std::vector<double>>(parameters.getParameter<edm::InputTag>("trailPz"))),
        trailLength_(consumes<std::vector<double>>(parameters.getParameter<edm::InputTag>("trailLength"))),
        trailX0_(consumes<std::vector<double>>(parameters.getParameter<edm::InputTag>("trailX0"))),
        trailEnergyLoss_(consumes<std::vector<double>>(parameters.getParameter<edm::InputTag>("trailEnergyLoss"))),
        trailMaterial_(consumes<std::vector<std::string>>(parameters.getParameter<edm::InputTag>("trailMaterial"))),
        trailVolume_(consumes<std::vector<std::string>>(parameters.getParameter<edm::InputTag>("trailVolume"))),
        field_(esConsumes(parameters.getParameter<edm::ESInputTag>("magneticField"))),
        geometry_(esConsumes()),
        output_(parameters.getParameter<std::string>("output")),
        maximumStepLengthMm_(parameters.getParameter<double>("maximumStepLengthMm")),
        maximumPathLengthCm_(parameters.getParameter<double>("maximumPathLengthCm")),
        cumulativeStrideCm_(parameters.getParameter<double>("cumulativeStrideCm")),
        externalBoundaryAbsZCm_(parameters.getParameter<double>("externalBoundaryAbsZCm")),
        useMeanEnergyLossJacobian_(parameters.getParameter<bool>("useMeanEnergyLossJacobian")),
        useFieldGradientJacobian_(parameters.getParameter<bool>("useFieldGradientJacobian")),
        useUnquenchedIonizationVariance_(parameters.getParameter<bool>("useUnquenchedIonizationVariance")) {
    if (!(maximumStepLengthMm_ > 0.) || !(maximumPathLengthCm_ > 0.) || !(cumulativeStrideCm_ > 0.) ||
        !(externalBoundaryAbsZCm_ > 0.) ||
        !std::isfinite(maximumStepLengthMm_) || !std::isfinite(maximumPathLengthCm_) ||
        !std::isfinite(cumulativeStrideCm_) || !std::isfinite(externalBoundaryAbsZCm_))
      throw cms::Exception("Configuration") << "Truth-trail closure lengths must be finite and positive";
  }

  void analyze(edm::Event const& event, edm::EventSetup const& setup) override {
    auto const& simTracks = event.get(simTracks_);
    auto const& simVertices = event.get(simVertices_);
    auto const& cscHits = event.get(cscSimHits_);
    auto const& ids = event.get(trailTrackId_);
    auto const& kinds = event.get(trailKind_);
    auto const& xs = event.get(trailX_);
    auto const& ys = event.get(trailY_);
    auto const& zs = event.get(trailZ_);
    auto const& pxs = event.get(trailPx_);
    auto const& pys = event.get(trailPy_);
    auto const& pzs = event.get(trailPz_);
    auto const& lengths = event.get(trailLength_);
    auto const& x0s = event.get(trailX0_);
    auto const& losses = event.get(trailEnergyLoss_);
    auto const& materials = event.get(trailMaterial_);
    auto const& volumes = event.get(trailVolume_);
    std::size_t const size = ids.size();
    for (auto const actual : {kinds.size(), xs.size(), ys.size(), zs.size(), pxs.size(), pys.size(), pzs.size(),
                              lengths.size(), x0s.size(), losses.size(), materials.size(), volumes.size()})
      if (actual != size)
        throw cms::Exception("LogicError") << "Truth-trail vectors have different sizes";

    auto const& geometry = setup.getData(geometry_);
    auto const& field = setup.getData(field_);
    Geant4ePropagator propagator(
        &field, "mu", oppositeToMomentum, 0.05, maximumStepLengthMm_, maximumPathLengthCm_);
    propagator.setUseConsistentBackwardCovariance(true);
    propagator.setUseMeanEnergyLossJacobian(useMeanEnergyLossJacobian_);
    propagator.setUseFieldGradientJacobian(useFieldGradientJacobian_);
    propagator.setUseUnquenchedIonizationVariance(useUnquenchedIonizationVariance_);
    Geant4ePropagator forwardPropagator(propagator);
    forwardPropagator.setPropagationDirection(alongMomentum);

    std::ostringstream json;
    json << std::setprecision(9) << "{\"run\":" << event.id().run() << ",\"lumi\":"
         << event.luminosityBlock() << ",\"event\":" << event.id().event() << ",\"tracks\":[";
    bool firstTrack = true;
    for (auto const& simTrack : simTracks) {
      if (simTrack.eventId().rawId() != 0 || simTrack.genpartIndex() < 0 || std::abs(simTrack.type()) != 13 ||
          simTrack.vertIndex() < 0 || static_cast<std::size_t>(simTrack.vertIndex()) >= simVertices.size())
        continue;
      std::vector<Checkpoint> checkpoints;
      for (std::size_t index = 0; index < size; ++index)
        if (ids[index] == static_cast<int>(simTrack.trackId()))
          checkpoints.push_back({kinds[index], xs[index], ys[index], zs[index], pxs[index], pys[index], pzs[index],
                                 lengths[index], x0s[index], losses[index], materials[index], volumes[index]});
      if (checkpoints.empty())
        continue;
      std::sort(checkpoints.begin(), checkpoints.end(), [](auto const& first, auto const& second) {
        return first.length < second.length;
      });

      auto const& vertex = simVertices[simTrack.vertIndex()].position();
      GlobalPoint const vertexPoint(vertex.x(), vertex.y(), vertex.z());
      GlobalVector const productionMomentum(
          simTrack.momentum().px(), simTrack.momentum().py(), simTrack.momentum().pz());
      PSimHit const* firstHit = nullptr;
      GlobalPoint hitPoint;
      GlobalVector hitMomentum;
      double firstPath = std::numeric_limits<double>::infinity();
      for (auto const& hit : cscHits) {
        if (hit.eventId() != simTrack.eventId() || hit.trackId() != simTrack.trackId())
          continue;
        auto const* det = geometry.idToDetUnit(DetId(hit.detUnitId()));
        if (!det)
          continue;
        auto const point = det->surface().toGlobal(hit.entryPoint());
        double const path = (point - vertexPoint).dot(productionMomentum.unit());
        if (path >= 0. && path < firstPath) {
          firstPath = path;
          firstHit = &hit;
          hitPoint = point;
          hitMomentum = det->surface().toGlobal(hit.momentumAtEntry());
        }
      }
      if (!firstHit)
        continue;

      std::size_t nearest = 0;
      double nearestDistance = std::numeric_limits<double>::infinity();
      for (std::size_t index = 0; index < checkpoints.size(); ++index) {
        auto const& point = checkpoints[index];
        double const distance = (GlobalPoint(point.x, point.y, point.z) - hitPoint).mag();
        if (distance < nearestDistance) {
          nearestDistance = distance;
          nearest = index;
        }
      }
      double const hitLength = checkpoints[nearest].length;
      std::vector<std::size_t> periodic;
      for (std::size_t index = 0; index <= nearest; ++index)
        if (checkpoints[index].kind == 0 || checkpoints[index].kind == 1)
          periodic.push_back(index);
      if (periodic.empty() || periodic.back() != nearest)
        periodic.push_back(nearest);

      int const charge = firstHit->particleType() > 0 ? -1 : 1;
      auto propagate = [&](GlobalPoint const& startPoint,
                           GlobalVector const& startMomentum,
                           Checkpoint const& target) {
        GlobalPoint const truthPoint(target.x, target.y, target.z);
        GlobalVector const truthMomentum(target.px, target.py, target.pz);
        auto const targetPlane = Plane::build(truthPoint, normalTo(truthMomentum));
        AlgebraicSymMatrix55 covariance;
        for (unsigned int coordinate = 0; coordinate < 5; ++coordinate)
          covariance(coordinate, coordinate) = 1.e-12;
        FreeTrajectoryState const start(GlobalTrajectoryParameters(startPoint, startMomentum, charge, &field),
                                        CurvilinearTrajectoryError(covariance));
        auto const result = propagator.propagateWithPath(start, *targetPlane);
        std::ostringstream record;
        record << std::setprecision(9) << "{\"valid\":" << (result.first.isValid() ? "true" : "false")
               << ",\"targetLengthCm\":" << target.length << ",\"targetZCm\":" << target.z
               << ",\"targetX0\":" << target.x0 << ",\"targetEnergyLossGeV\":" << target.energyLoss
               << ",\"material\":" << std::quoted(target.material)
               << ",\"volume\":" << std::quoted(target.volume);
        if (result.first.isValid()) {
          auto const position = result.first.globalPosition();
          auto const momentum = result.first.globalMomentum();
          record << ",\"pathCm\":" << result.second
                 << ",\"positionResidualCm\":" << (position - truthPoint).mag()
                 << ",\"angleResidualRad\":" << angle(momentum, truthMomentum)
                 << ",\"relativeMomentumResidual\":" << momentum.mag() / truthMomentum.mag() - 1.;
          if (result.first.hasError()) {
            FreeTrajectoryState const truth(
                GlobalTrajectoryParameters(truthPoint, truthMomentum, charge, &field),
                CurvilinearTrajectoryError(covariance));
            TrajectoryStateOnSurface const actual(truth, result.first.surface());
            Eigen::Matrix<double, 5, 1> difference;
            Eigen::Matrix<double, 5, 5> noise;
            for (unsigned int row = 0; row < 5; ++row) {
              difference[row] =
                  actual.localParameters().vector()[row] - result.first.localParameters().vector()[row];
              for (unsigned int column = 0; column < 5; ++column)
                noise(row, column) = result.first.localError().matrix()(row, column);
            }
            Eigen::LLT<Eigen::Matrix<double, 5, 5>> solve(noise);
            if (solve.info() == Eigen::Success && difference.allFinite()) {
              record << ",\"noiseChi2\":" << difference.dot(solve.solve(difference)) << ",\"pulls\":[";
              for (unsigned int coordinate = 0; coordinate < 5; ++coordinate) {
                if (coordinate)
                  record << ',';
                record << difference[coordinate] / std::sqrt(noise(coordinate, coordinate));
              }
              record << ']';
            }
          }
        }
        record << '}';
        return record.str();
      };

      if (!firstTrack)
        json << ',';
      firstTrack = false;
      json << "{\"trackId\":" << simTrack.trackId() << ",\"pdgId\":" << simTrack.type()
           << ",\"firstCscDetId\":" << firstHit->detUnitId() << ",\"firstCscPathCm\":" << firstPath
           << ",\"firstCscPositionCm\":[" << hitPoint.x() << ',' << hitPoint.y() << ',' << hitPoint.z() << ']'
           << ",\"firstCscMomentumGeV\":[" << hitMomentum.x() << ',' << hitMomentum.y() << ','
           << hitMomentum.z() << "]"
           << ",\"nearestCheckpointDistanceCm\":" << nearestDistance
           << ",\"nearestCheckpointLengthCm\":" << hitLength << ",\"localIntervals\":[";
      bool first = true;
      for (std::size_t item = 1; item < periodic.size(); ++item) {
        auto const& earlier = checkpoints[periodic[item - 1]];
        auto const& later = checkpoints[periodic[item]];
        if (!(later.length - earlier.length > 1.e-3))
          continue;
        if (!first)
          json << ',';
        first = false;
        json << "{\"startLengthCm\":" << later.length << ",\"deltaLengthCm\":"
             << later.length - earlier.length << ",\"deltaX0\":" << later.x0 - earlier.x0
             << ",\"deltaEnergyLossGeV\":" << later.energyLoss - earlier.energyLoss
             << ",\"startMaterial\":" << std::quoted(later.material) << ",\"result\":"
             << propagate(GlobalPoint(later.x, later.y, later.z), GlobalVector(later.px, later.py, later.pz), earlier)
             << '}';
      }
      json << "],\"cumulative\":[";
      first = true;
      double lastDistance = -std::numeric_limits<double>::infinity();
      for (auto item = periodic.rbegin(); item != periodic.rend(); ++item) {
        auto const& target = checkpoints[*item];
        double const distance = hitLength - target.length;
        if (distance < 100.)
          continue;
        bool const endpoint = *item == periodic.front() || *item == periodic.back();
        if (!endpoint && distance < lastDistance + cumulativeStrideCm_)
          continue;
        if (!first)
          json << ',';
        first = false;
        json << propagate(hitPoint, hitMomentum, target);
        lastDistance = distance;
      }
      json << ']';

      // The CSC-seeded cumulative propagation necessarily includes stochastic
      // detector material.  Seed separately from the last truth checkpoint on
      // the source side of the requested |z| boundary so that the long LSS leg
      // can be tested without conflating it with CMS material.
      auto externalSeed = periodic.end();
      double const sourceSide = checkpoints.front().z < 0. ? -1. : 1.;
      for (auto item = periodic.begin(); item != periodic.end(); ++item) {
        auto const& candidate = checkpoints[*item];
        if (sourceSide * candidate.z >= externalBoundaryAbsZCm_)
          externalSeed = item;
      }
      json << ",\"externalBoundaryAbsZCm\":" << externalBoundaryAbsZCm_;
      if (externalSeed == periodic.end()) {
        json << ",\"externalSeedFound\":false,\"externalCumulative\":[]}";
        continue;
      }

      auto const& seed = checkpoints[*externalSeed];
      json << ",\"externalSeedFound\":true,\"externalSeedLengthCm\":" << seed.length
           << ",\"externalSeedPositionCm\":[" << seed.x << ',' << seed.y << ',' << seed.z << ']'
           << ",\"externalSeedMomentumGeV\":[" << seed.px << ',' << seed.py << ',' << seed.pz << ']'
           << ",\"externalSeedX0\":" << seed.x0
           << ",\"externalSeedEnergyLossGeV\":" << seed.energyLoss
           << ",\"externalSeedMaterial\":" << std::quoted(seed.material)
           << ",\"externalSeedVolume\":" << std::quoted(seed.volume)
           << ",\"externalCumulative\":[";
      first = true;
      lastDistance = -std::numeric_limits<double>::infinity();
      for (auto item = std::make_reverse_iterator(externalSeed + 1); item != periodic.rend(); ++item) {
        auto const& target = checkpoints[*item];
        double const distance = seed.length - target.length;
        if (distance < 100.)
          continue;
        bool const endpoint = *item == periodic.front();
        if (!endpoint && distance < lastDistance + cumulativeStrideCm_)
          continue;
        if (!first)
          json << ',';
        first = false;
        json << propagate(GlobalPoint(seed.x, seed.y, seed.z), GlobalVector(seed.px, seed.py, seed.pz), target);
        lastDistance = distance;
      }
      json << "]";

      // This is the genuinely deterministic field check.  Use Geant4e for
      // both directions, so the return residual cannot contain the random
      // Geant4 scattering encoded in the truth checkpoints.
      auto const sourcePlane = Plane::build(vertexPoint, normalTo(productionMomentum));
      GlobalPoint const seedPoint(seed.x, seed.y, seed.z);
      GlobalVector const seedMomentum(seed.px, seed.py, seed.pz);
      auto const seedPlane = Plane::build(seedPoint, normalTo(seedMomentum));
      AlgebraicSymMatrix55 roundTripCovariance;
      for (unsigned int coordinate = 0; coordinate < 5; ++coordinate)
        roundTripCovariance(coordinate, coordinate) = 1.e-12;
      FreeTrajectoryState const sourceState(
          GlobalTrajectoryParameters(vertexPoint, productionMomentum, charge, &field),
          CurvilinearTrajectoryError(roundTripCovariance));
      auto const forwardResult = forwardPropagator.propagateWithPath(sourceState, *seedPlane);
      json << ",\"deterministicRoundTrip\":{\"forwardValid\":"
           << (forwardResult.first.isValid() ? "true" : "false");
      if (forwardResult.first.isValid()) {
        auto const forwardPosition = forwardResult.first.globalPosition();
        auto const forwardMomentum = forwardResult.first.globalMomentum();
        json << ",\"forwardPathCm\":" << forwardResult.second
             << ",\"forwardTruthPositionResidualCm\":" << (forwardPosition - seedPoint).mag()
             << ",\"forwardTruthAngleResidualRad\":" << angle(forwardMomentum, seedMomentum)
             << ",\"forwardTruthRelativeMomentumResidual\":"
             << forwardMomentum.mag() / seedMomentum.mag() - 1.;
        auto const backwardResult = propagator.propagateWithPath(forwardResult.first, *sourcePlane);
        json << ",\"backwardValid\":" << (backwardResult.first.isValid() ? "true" : "false");
        if (backwardResult.first.isValid()) {
          auto const backwardPosition = backwardResult.first.globalPosition();
          auto const backwardMomentum = backwardResult.first.globalMomentum();
          json << ",\"backwardPathCm\":" << backwardResult.second
               << ",\"returnPositionResidualCm\":" << (backwardPosition - vertexPoint).mag()
               << ",\"returnAngleResidualRad\":" << angle(backwardMomentum, productionMomentum)
               << ",\"returnRelativeMomentumResidual\":"
               << backwardMomentum.mag() / productionMomentum.mag() - 1.;
        }
      }
      json << "}}";
    }
    json << "]}";
    events_.push_back(json.str());
  }

  void endJob() override {
    std::ofstream output(output_);
    if (!output)
      throw cms::Exception("FileWriteError") << "Cannot write " << output_;
    output << std::setprecision(17)
           << "{\"schema_version\":2,\"maximum_step_length_mm\":" << maximumStepLengthMm_
           << ",\"maximum_path_length_cm\":" << maximumPathLengthCm_
           << ",\"cumulative_stride_cm\":" << cumulativeStrideCm_
           << ",\"external_boundary_abs_z_cm\":" << externalBoundaryAbsZCm_
           << ",\"use_mean_energy_loss_jacobian\":"
           << (useMeanEnergyLossJacobian_ ? "true" : "false")
           << ",\"use_field_gradient_jacobian\":" << (useFieldGradientJacobian_ ? "true" : "false")
           << ",\"use_unquenched_ionization_variance\":"
           << (useUnquenchedIonizationVariance_ ? "true" : "false") << ",\"events\":[\n";
    for (std::size_t index = 0; index < events_.size(); ++index) {
      if (index)
        output << ",\n";
      output << events_[index];
    }
    output << "\n]}\n";
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription description;
    description.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
    description.add<edm::InputTag>("simVertices", edm::InputTag("g4SimHits"));
    description.add<edm::InputTag>("cscSimHits", edm::InputTag("g4SimHits", "MuonCSCHits"));
    auto trail = [&description](char const* name, char const* suffix) {
      description.add<edm::InputTag>(name, edm::InputTag("g4SimHits", std::string("shiftMuonTruthTrail") + suffix));
    };
    trail("trailTrackId", "TrackId");
    trail("trailKind", "Kind");
    trail("trailX", "X");
    trail("trailY", "Y");
    trail("trailZ", "Z");
    trail("trailPx", "Px");
    trail("trailPy", "Py");
    trail("trailPz", "Pz");
    trail("trailLength", "TrackLength");
    trail("trailX0", "CumulativeX0");
    trail("trailEnergyLoss", "CumulativeEnergyLoss");
    trail("trailMaterial", "Material");
    trail("trailVolume", "Volume");
    description.add<edm::ESInputTag>("magneticField", edm::ESInputTag("", ""));
    description.add<std::string>("output");
    description.add<double>("maximumStepLengthMm", 10.);
    description.add<double>("maximumPathLengthCm", 20000.);
    description.add<double>("cumulativeStrideCm", 1000.);
    description.add<double>("externalBoundaryAbsZCm", 3000.);
    description.add<bool>("useMeanEnergyLossJacobian", true);
    description.add<bool>("useFieldGradientJacobian", true);
    description.add<bool>("useUnquenchedIonizationVariance", true);
    descriptions.add("shiftMuonTruthTrailClosureAnalyzer", description);
  }

private:
  edm::EDGetTokenT<edm::SimTrackContainer> simTracks_;
  edm::EDGetTokenT<edm::SimVertexContainer> simVertices_;
  edm::EDGetTokenT<edm::PSimHitContainer> cscSimHits_;
  edm::EDGetTokenT<std::vector<int>> trailTrackId_, trailKind_;
  edm::EDGetTokenT<std::vector<double>> trailX_, trailY_, trailZ_, trailPx_, trailPy_, trailPz_, trailLength_,
      trailX0_, trailEnergyLoss_;
  edm::EDGetTokenT<std::vector<std::string>> trailMaterial_, trailVolume_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> field_;
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geometry_;
  std::string output_;
  double maximumStepLengthMm_, maximumPathLengthCm_, cumulativeStrideCm_, externalBoundaryAbsZCm_;
  bool useMeanEnergyLossJacobian_, useFieldGradientJacobian_, useUnquenchedIonizationVariance_;
  std::vector<std::string> events_;
};

DEFINE_FWK_MODULE(ShiftMuonTruthTrailClosureAnalyzer);
