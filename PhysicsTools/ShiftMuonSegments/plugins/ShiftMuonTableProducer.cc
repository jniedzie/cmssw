#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkClonerImpl.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace {
  struct PropagatedState {
    bool valid = false;
    GlobalPoint position;
    GlobalVector momentum;
    double path = 0.;
  };

  struct HitFingerprint {
    uint32_t detId;
    float localX;
    float localY;
  };

  struct Candidate {
    reco::Track const* track;
    int source;
    unsigned int sourceIndex;
    std::vector<HitFingerprint> hitFingerprints;
    PropagatedState targetLineState;
    int timingDirectionSign = 0;
    unsigned int timingMeasurements = 0;
    double timingChi2 = 0.;
    double timingDeltaChi2 = 0.;
    int physicalDirectionSign = 0;
    bool directionalRefitAttempted = false;
    bool directionalRefitValid = false;
    unsigned int directionalRefitHits = 0;
    double directionalRefitChi2 = 0.;
    double directionalRefitNdof = 0.;
    double directionalRefitUpstreamPt = 0.;
    double preRefitPt = 0.;
  };

  PropagatedState propagateStateToTargetLine(GlobalPoint const& position,
                                             GlobalVector const& momentum,
                                             int charge,
                                             Propagator const& vacuumPropagator,
                                             Propagator const* materialPropagator = nullptr,
                                             int sourceSide = 0) {
    FreeTrajectoryState const start(position, momentum, charge, vacuumPropagator.magneticField());
    // SteppingHelixPropagator uses an internal approximate R-Z material map,
    // not the detailed detector geometry.  Its last modeled CMS endcap
    // structures end near |z|=1073 cm, while more distant regions are not a
    // reliable description of the evacuated SHIFT-to-CMS flight path.  Stop
    // the material-aware transport just outside that envelope, then continue
    // in vacuum to the external target line.
    constexpr double materialBoundaryZ = 1100.;
    FreeTrajectoryState vacuumStart = start;
    double path = 0.;
    if (materialPropagator && sourceSide != 0 && sourceSide * position.z() < materialBoundaryZ) {
      auto const boundary = Plane::build(GlobalPoint(0., 0., sourceSide * materialBoundaryZ),
                                         Surface::RotationType());
      auto const toBoundary = materialPropagator->propagateWithPath(start, *boundary);
      if (!toBoundary.first.isValid() || !toBoundary.first.freeState())
        return {};
      vacuumStart = *toBoundary.first.freeState();
      path += toBoundary.second;
    }

    auto const propagated = vacuumPropagator.propagateWithPath(
        vacuumStart, GlobalPoint(0., 0., -1.), GlobalPoint(0., 0., 1.));
    // A failed FreeTrajectoryState propagation is returned default-constructed
    // (zero charge and a null field pointer), unlike a TSOS it has no
    // isValid() accessor.  Check its sentinel before accessing position.
    if (propagated.first.charge() == 0)
      return {};
    auto const resultPosition = propagated.first.position();
    auto const resultMomentum = propagated.first.momentum();
    bool const valid = resultMomentum.mag2() > 0. && std::isfinite(resultPosition.x()) &&
                       std::isfinite(resultPosition.y()) && std::isfinite(resultPosition.z()) &&
                       std::isfinite(path + propagated.second);
    return {valid, resultPosition, resultMomentum, path + propagated.second};
  }

  PropagatedState propagateToTargetLine(reco::Track const& track,
                                        Propagator const& vacuumPropagator,
                                        int travelSign = 0,
                                        Propagator const* materialPropagator = nullptr,
                                        int sourceSide = 0) {
    // Once the event side is known, start from the fitted endpoint closest in
    // z to the external source.  This avoids beginning the backward transport
    // from a state that has already crossed additional detector material.
    auto const endpointDelta = track.outerPosition() - track.innerPosition();
    bool const useOuter = travelSign != 0 && travelSign * track.innerMomentum().Dot(endpointDelta) < 0.;
    auto const& endpoint = useOuter ? track.outerPosition() : track.innerPosition();
    auto const& endpointMomentum = useOuter ? track.outerMomentum() : track.innerMomentum();
    GlobalPoint const position(endpoint.x(), endpoint.y(), endpoint.z());
    // Momentum and charge must be reversed together to describe the same
    // fitted helix with the physical time orientation.  Since the external
    // target lies behind the incoming particle, the any-direction propagator
    // then selects the opposite-to-momentum geometrical solution.
    double const sign = travelSign == 0 ? 1. : travelSign;
    GlobalVector const momentum(
        sign * endpointMomentum.x(), sign * endpointMomentum.y(), sign * endpointMomentum.z());
    return propagateStateToTargetLine(position,
                                      momentum,
                                      sign * track.charge(),
                                      vacuumPropagator,
                                      materialPropagator,
                                      sourceSide);
  }

  struct DirectionalRefitResult {
    bool valid = false;
    unsigned int hits = 0;
    double chi2 = 0.;
    double ndof = 0.;
    double upstreamPt = 0.;
    PropagatedState targetLineState;
  };

  DirectionalRefitResult directionalRefit(reco::Track const& track,
                                          int directionSign,
                                          int sourceSide,
                                          MagneticField const& magneticField,
                                          TransientTrackingRecHitBuilder const& hitBuilder,
                                          TkCloner const& hitCloner,
                                          Propagator const& vacuumPropagator,
                                          Propagator const& materialPropagator,
                                          double errorRescale,
                                          double maxHitChi2) {
    struct OrderedHit {
      TransientTrackingRecHit::RecHitPointer hit;
      double sourceCoordinate;
    };
    std::vector<OrderedHit> orderedHits;
    for (auto hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
      if (!(*hit)->isValid() || (*hit)->geographicalId().det() != DetId::Muon)
        continue;
      auto transientHit = hitBuilder.build(&**hit);
      if (!transientHit || !transientHit->isValid() || !transientHit->det())
        continue;
      orderedHits.push_back({transientHit, sourceSide * transientHit->globalPosition().z()});
    }
    if (orderedHits.size() < 3)
      return {};
    std::stable_sort(orderedHits.begin(), orderedHits.end(), [](OrderedHit const& first, OrderedHit const& second) {
      return first.sourceCoordinate > second.sourceCoordinate;
    });

    bool const useOuter = sourceSide * track.outerPosition().z() > sourceSide * track.innerPosition().z();
    auto const original = useOuter ? trajectoryStateTransform::outerFreeState(track, &magneticField)
                                   : trajectoryStateTransform::innerFreeState(track, &magneticField);
    if (!original.hasError())
      return {};
    double const sign = directionSign == 0 ? 1. : directionSign;
    auto const rawMomentum = original.momentum();
    GlobalVector const momentum(sign * rawMomentum.x(), sign * rawMomentum.y(), sign * rawMomentum.z());
    GlobalTrajectoryParameters const parameters(
        original.position(), momentum, sign * original.charge(), &magneticField);
    FreeTrajectoryState start(parameters, original.curvilinearError());
    start.rescaleError(errorRescale);

    Trajectory::RecHitContainer fitHits;
    fitHits.reserve(orderedHits.size());
    for (auto const& orderedHit : orderedHits)
      fitHits.push_back(orderedHit.hit);

    auto const firstPredicted = materialPropagator.propagate(start, orderedHits.front().hit->det()->surface());
    if (!firstPredicted.isValid())
      return {};

    KFUpdator updator;
    Chi2MeasurementEstimator estimator(maxHitChi2);
    KFTrajectoryFitter fitter(materialPropagator, updator, estimator, 3, nullptr, &hitCloner);
    TrajectorySeed const seed(
        PTrajectoryStateOnDet(), TrajectorySeed::RecHitContainer(), alongMomentum);
    auto const filtered = fitter.fitOne(seed, fitHits, firstPredicted, TrajectoryFitter::standard);
    if (!filtered.isValid() || filtered.foundHits() < 3)
      return {};

    KFTrajectorySmoother smoother(
        materialPropagator, updator, estimator, static_cast<float>(errorRescale), 3);
    smoother.setHitCloner(&hitCloner);
    auto const smoothed = smoother.trajectory(filtered);
    if (!smoothed.isValid() || smoothed.foundHits() < 3 || smoothed.empty())
      return {};

    // The smoother returns measurements in the reverse order of the forward
    // source-to-CMS fit.  Its last measurement is therefore the source-facing
    // state, consistently combining the forward and backward information
    // without applying the same measurement twice.
    auto const& upstream = smoothed.lastMeasurement().updatedState();
    if (!upstream.isValid())
      return {};
    auto const upstreamMomentum = upstream.globalMomentum();
    auto const targetLineState = propagateStateToTargetLine(upstream.globalPosition(),
                                                            upstreamMomentum,
                                                            upstream.charge(),
                                                            vacuumPropagator,
                                                            &materialPropagator,
                                                            sourceSide);
    if (!targetLineState.valid)
      return {};
    auto const refitHits = static_cast<unsigned int>(smoothed.foundHits());
    return {true,
            refitHits,
            smoothed.chiSquared(),
            std::max(1., 2. * static_cast<double>(refitHits) - 5.),
            upstreamMomentum.perp(),
            targetLineState};
  }

  struct TimingResult {
    int directionSign = 0;
    unsigned int measurements = 0;
    double chi2 = 0.;
    double deltaChi2 = 0.;
  };

  TimingResult timingDirection(reco::Track const& track,
                               GlobalTrackingGeometry const& geometry,
                               unsigned int minMeasurements,
                               double minDeltaChi2) {
    struct TimedPoint { GlobalPoint position; double time; double sigma; };
    std::vector<TimedPoint> points;
    for (auto hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
      if (!(*hit)->isValid())
        continue;
      double time = 0., sigma = 0.;
      if (auto const* segment = dynamic_cast<CSCSegment const*>(&**hit)) {
        time = segment->time();
        sigma = 7.;
      } else if (auto const* segment = dynamic_cast<DTRecSegment4D const*>(&**hit)) {
        if (!segment->hasPhi() || !segment->phiSegment()->ist0Valid())
          continue;
        time = segment->phiSegment()->t0();
        sigma = 3.;
      } else {
        continue;
      }
      auto const* detector = geometry.idToDet((*hit)->geographicalId());
      if (!detector || !std::isfinite(time) || std::abs(time) > 1.e4)
        continue;
      points.push_back({detector->surface().toGlobal((*hit)->localPosition()), time, sigma});
    }
    if (points.size() < minMeasurements)
      return {0, static_cast<unsigned int>(points.size()), 0., 0.};

    constexpr double inverseSpeedOfLight = 1. / 29.9792458;  // ns/cm
    auto const axis = track.momentum().unit();
    auto score = [&](double sign) {
      auto coordinate = [&axis](GlobalPoint const& position) {
        return axis.x() * position.x() + axis.y() * position.y() + axis.z() * position.z();
      };
      std::vector<double> offsets;
      offsets.reserve(points.size());
      for (auto const& point : points)
        offsets.push_back(point.time - sign * coordinate(point.position) * inverseSpeedOfLight);
      std::sort(offsets.begin(), offsets.end());
      double const offset = offsets.size() % 2 != 0
                                ? offsets[offsets.size() / 2]
                                : 0.5 * (offsets[offsets.size() / 2 - 1] + offsets[offsets.size() / 2]);
      double robustChi2 = 0.;
      for (auto const& point : points) {
        double const residual = point.time - offset - sign * coordinate(point.position) * inverseSpeedOfLight;
        double const pull = std::abs(residual / point.sigma);
        // Huber loss keeps one mistimed chamber from deciding the direction.
        robustChi2 += pull <= 3. ? pull * pull : 6. * pull - 9.;
      }
      return robustChi2;
    };
    double const alongChi2 = score(1.);
    double const oppositeChi2 = score(-1.);
    double const delta = std::abs(alongChi2 - oppositeChi2);
    int const sign = delta >= minDeltaChi2 ? (alongChi2 < oppositeChi2 ? 1 : -1) : 0;
    return {sign,
            static_cast<unsigned int>(points.size()),
            std::min(alongChi2, oppositeChi2),
            delta};
  }

  void appendHitFingerprints(TrackingRecHit const& hit, std::vector<HitFingerprint>& result) {
    auto const components = hit.recHits();
    if (!components.empty()) {
      for (auto const* component : components)
        if (component && component->isValid())
          appendHitFingerprints(*component, result);
      return;
    }
    auto const position = hit.localPosition();
    result.push_back({hit.rawId(), position.x(), position.y()});
  }

  std::vector<HitFingerprint> hitFingerprints(reco::Track const& track) {
    std::vector<HitFingerprint> result;
    for (auto hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit)
      if ((*hit)->isValid())
        appendHitFingerprints(**hit, result);
    return result;
  }

  double directionAngle(reco::Track const& first, reco::Track const& second) {
    auto const firstDirection = first.momentum().unit();
    auto const secondDirection = second.momentum().unit();
    double const cosine = std::clamp(std::abs(firstDirection.Dot(secondDirection)), 0., 1.);
    return std::acos(cosine);
  }

  double directionDeltaR(GlobalVector const& momentum, reco::GenParticle const& particle) {
    auto deltaPhi = [](double first, double second) {
      double value = std::remainder(first - second, 2. * M_PI);
      return std::abs(value);
    };
    double const pt = momentum.perp();
    double const eta = pt > 0. ? std::asinh(momentum.z() / pt)
                               : std::copysign(std::numeric_limits<double>::infinity(), momentum.z());
    double const phi = std::atan2(momentum.y(), momentum.x());
    double const direct = std::hypot(eta - particle.eta(), deltaPhi(phi, particle.phi()));
    double const reverse = std::hypot(eta + particle.eta(), deltaPhi(phi + M_PI, particle.phi()));
    return std::min(direct, reverse);
  }

  double pointLineDistance(reco::Track::Point const& point, reco::Track const& line) {
    auto const delta = point - line.referencePoint();
    return delta.Cross(line.momentum().unit()).R();
  }

  double symmetricLineDistance(reco::Track const& first, reco::Track const& second) {
    // Comparing detector states avoids merging two distinct tracks which
    // merely intersect at the common SHIFT production region.
    return std::max({pointLineDistance(first.innerPosition(), second),
                     pointLineDistance(first.outerPosition(), second),
                     pointLineDistance(second.innerPosition(), first),
                     pointLineDistance(second.outerPosition(), first)});
  }

  std::pair<double, double> transverseLinePca(reco::Track const& track) {
    auto const& position = track.innerPosition();
    auto const& momentum = track.innerMomentum();
    double const pt2 = momentum.x() * momentum.x() + momentum.y() * momentum.y();
    if (pt2 == 0.)
      return {position.rho(), position.z()};
    double const scale = -(position.x() * momentum.x() + position.y() * momentum.y()) / pt2;
    return {std::hypot(position.x() + scale * momentum.x(), position.y() + scale * momentum.y()),
            position.z() + scale * momentum.z()};
  }

  std::pair<double, double> chordLinePca(reco::Track const& track) {
    // Independent direction diagnostic made only from the fitted endpoint
    // positions.  It is insensitive to the momentum parameterisation used by
    // the standalone Kalman fit and exposes polar-slope biases directly.
    auto const& inner = track.innerPosition();
    auto const delta = track.outerPosition() - inner;
    double const transverse2 = delta.x() * delta.x() + delta.y() * delta.y();
    if (transverse2 == 0.)
      return {inner.rho(), inner.z()};
    double const scale = -(inner.x() * delta.x() + inner.y() * delta.y()) / transverse2;
    return {std::hypot(inner.x() + scale * delta.x(), inner.y() + scale * delta.y()),
            inner.z() + scale * delta.z()};
  }

  double shiftDirectionSign(GlobalVector const& momentum, int sourceSide) {
    return sourceSide * momentum.z() > 0. ? -1. : 1.;
  }

  struct StraightLineApproach {
    bool valid;
    GlobalPoint midpoint;
    double distance;
  };

  struct CommonLineVertex {
    bool valid = false;
    GlobalPoint position;
    std::array<double, 3> error{{0., 0., 0.}};
    double chi2 = 0.;
    double ndof = 3.;
  };

  bool invertSymmetric3x3(std::array<std::array<double, 3>, 3> const& matrix,
                          std::array<std::array<double, 3>, 3>& inverse) {
    double const a = matrix[0][0], b = matrix[0][1], c = matrix[0][2];
    double const d = matrix[1][1], e = matrix[1][2], f = matrix[2][2];
    double const determinant = a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
    if (!std::isfinite(determinant) || std::abs(determinant) < 1.e-18)
      return false;
    inverse = {{{(d * f - e * e) / determinant, (c * e - b * f) / determinant, (b * e - c * d) / determinant},
                {(c * e - b * f) / determinant, (a * f - c * c) / determinant, (b * c - a * e) / determinant},
                {(b * e - c * d) / determinant, (b * c - a * e) / determinant, (a * d - b * b) / determinant}}};
    return true;
  }

  CommonLineVertex commonLineVertex(PropagatedState const& first,
                                    PropagatedState const& second,
                                    double lineResolution,
                                    double beamLineResolution) {
    // Minimise the distances to the two unbounded track lines together with
    // a loose x=y=0 target-line prior.  Unlike TransientTrack propagation,
    // this analytic fit has no detector-volume boundary at |z|~10 m.
    std::array<std::array<double, 3>, 3> normal{};
    std::array<double, 3> rhs{};
    double const lineWeight = 1. / (lineResolution * lineResolution);
    for (auto const* state : {&first, &second}) {
      auto const point = state->position;
      auto const direction = state->momentum.unit();
      std::array<double, 3> const p{{point.x(), point.y(), point.z()}};
      std::array<double, 3> const u{{direction.x(), direction.y(), direction.z()}};
      for (unsigned int row = 0; row < 3; ++row)
        for (unsigned int column = 0; column < 3; ++column) {
          double const projector = (row == column ? 1. : 0.) - u[row] * u[column];
          normal[row][column] += lineWeight * projector;
          rhs[row] += lineWeight * projector * p[column];
        }
    }
    double const beamWeight = 1. / (beamLineResolution * beamLineResolution);
    normal[0][0] += beamWeight;
    normal[1][1] += beamWeight;
    std::array<std::array<double, 3>, 3> covariance{};
    if (!invertSymmetric3x3(normal, covariance))
      return {};
    std::array<double, 3> fitted{};
    for (unsigned int row = 0; row < 3; ++row)
      for (unsigned int column = 0; column < 3; ++column)
        fitted[row] += covariance[row][column] * rhs[column];
    if (!std::isfinite(fitted[0]) || !std::isfinite(fitted[1]) || !std::isfinite(fitted[2]))
      return {};
    GlobalPoint const position(fitted[0], fitted[1], fitted[2]);
    double chi2 = (fitted[0] * fitted[0] + fitted[1] * fitted[1]) * beamWeight;
    for (auto const* state : {&first, &second}) {
      GlobalVector const delta(fitted[0] - state->position.x(),
                               fitted[1] - state->position.y(),
                               fitted[2] - state->position.z());
      chi2 += std::pow(delta.cross(state->momentum.unit()).mag(), 2) * lineWeight;
    }
    return {true,
            position,
            {{std::sqrt(std::max(0., covariance[0][0])),
              std::sqrt(std::max(0., covariance[1][1])),
              std::sqrt(std::max(0., covariance[2][2]))}},
            chi2,
            3.};
  }

  StraightLineApproach straightLineApproach(PropagatedState const& first, PropagatedState const& second) {
    auto const firstPoint = first.position;
    auto const secondPoint = second.position;
    auto const firstDirection = first.momentum.unit();
    auto const secondDirection = second.momentum.unit();
    auto const separation = firstPoint - secondPoint;
    double const dot = firstDirection.dot(secondDirection);
    double const denominator = 1. - dot * dot;
    if (denominator < 1.e-10) {
      GlobalPoint const midpoint(0.5 * (firstPoint.x() + secondPoint.x()),
                                 0.5 * (firstPoint.y() + secondPoint.y()),
                                 0.5 * (firstPoint.z() + secondPoint.z()));
      return {true, midpoint, (firstPoint - secondPoint).mag()};
    }
    double const firstScale = (dot * secondDirection.dot(separation) - firstDirection.dot(separation)) / denominator;
    double const secondScale = (secondDirection.dot(separation) - dot * firstDirection.dot(separation)) / denominator;
    auto const firstClosest = firstPoint + firstScale * firstDirection;
    auto const secondClosest = secondPoint + secondScale * secondDirection;
    GlobalPoint const midpoint(0.5 * (firstClosest.x() + secondClosest.x()),
                               0.5 * (firstClosest.y() + secondClosest.y()),
                               0.5 * (firstClosest.z() + secondClosest.z()));
    double const distance = (firstClosest - secondClosest).mag();
    return {std::isfinite(distance) && std::isfinite(midpoint.z()), midpoint, distance};
  }
}  // namespace

class ShiftMuonTableProducer : public edm::stream::EDProducer<> {
public:
  explicit ShiftMuonTableProducer(edm::ParameterSet const& parameters)
      : dsaToken_(consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("dsaTracks"))),
        cosmicToken_(consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("cosmicTracks"))),
        traversingToken_(consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("traversingTracks"))),
        genParticlesToken_(
            consumes<reco::GenParticleCollection>(parameters.getParameter<edm::InputTag>("genParticles"))),
        magneticFieldToken_(esConsumes()),
        trackingGeometryToken_(esConsumes()),
        muonRecHitBuilderToken_(esConsumes(edm::ESInputTag(
            "", parameters.getParameter<std::string>("muonRecHitBuilder")))),
        directionalRefitErrorRescale_(parameters.getParameter<double>("directionalRefitErrorRescale")),
        directionalRefitMaxHitChi2_(parameters.getParameter<double>("directionalRefitMaxHitChi2")),
        minSharedHitFraction_(parameters.getParameter<double>("minSharedHitFraction")),
        minSharedDetIds_(parameters.getParameter<unsigned int>("minSharedDetIds")),
        maxDuplicateAngle_(parameters.getParameter<double>("maxDuplicateAngle")),
        maxDuplicateLineDistance_(parameters.getParameter<double>("maxDuplicateLineDistance")),
        minTimingMeasurements_(parameters.getParameter<unsigned int>("minTimingMeasurements")),
        minTimingDeltaChi2_(parameters.getParameter<double>("minTimingDeltaChi2")),
        maxTargetLineDca_(parameters.getParameter<double>("maxTargetLineDca")),
        minAbsEta_(parameters.getParameter<double>("minAbsEta")),
        originTransverseResolution_(parameters.getParameter<double>("originTransverseResolution")),
        originZResolution_(parameters.getParameter<double>("originZResolution")),
        commonVertexLineResolution_(parameters.getParameter<double>("commonVertexLineResolution")),
        commonVertexBeamLineResolution_(parameters.getParameter<double>("commonVertexBeamLineResolution")),
        maxPairOriginNormalizedChi2_(parameters.getParameter<double>("maxPairOriginNormalizedChi2")),
        maxPairDca_(parameters.getParameter<double>("maxPairDca")),
        maxDimuonVertices_(parameters.getParameter<unsigned int>("maxDimuonVertices")),
        requireOppositeSign_(parameters.getParameter<bool>("requireOppositeSign")),
        maxGenDeltaR_(parameters.getParameter<double>("maxGenDeltaR")) {
    produces<nanoaod::FlatTable>();
    produces<nanoaod::FlatTable>("ShiftDimuonVertex");
  }

  void produce(edm::Event& event, edm::EventSetup const& setup) override {
    auto const dsa = event.getHandle(dsaToken_);
    auto const cosmic = event.getHandle(cosmicToken_);
    auto const traversing = event.getHandle(traversingToken_);
    auto const genParticles = event.getHandle(genParticlesToken_);
    auto const& magneticField = setup.getData(magneticFieldToken_);
    // Use material effects only within the modeled CMS envelope.  Outside it,
    // propagate through the evacuated SHIFT-to-CMS flight path without dE/dx.
    SteppingHelixPropagator materialPropagator(&magneticField, anyDirection);
    // In this API the argument is `noMaterial`, so false enables dE/dx.
    materialPropagator.setMaterialMode(false);
    materialPropagator.setUseMagVolumes(true);
    materialPropagator.setUseMatVolumes(true);
    materialPropagator.applyRadX0Correction(true);
    SteppingHelixPropagator vacuumPropagator(&magneticField, anyDirection);
    vacuumPropagator.setMaterialMode(true);
    vacuumPropagator.setUseMagVolumes(true);
    auto const& trackingGeometry = setup.getData(trackingGeometryToken_);
    auto const& muonRecHitBuilder = setup.getData(muonRecHitBuilderToken_);
    // KFTrajectoryFitter/Smoother use TkCloner as a generic rechit clone-or-
    // reuse helper. Muon transient rechits cannot improve with a tracker CPE,
    // so a default cloner simply returns their existing shared pointers and
    // never accesses tracker calibration objects.
    TkClonerImpl const hitCloner;

    std::vector<Candidate> candidates;
    auto append = [&candidates, &vacuumPropagator](auto const& handle, int source) {
      if (!handle.isValid())
        return;
      unsigned int index = 0;
      for (auto const& track : *handle) {
        candidates.push_back({&track,
                              source,
                              index++,
                              hitFingerprints(track),
                              propagateToTargetLine(track, vacuumPropagator)});
      }
    };
    // Preserve the public source numbering; representative precedence is
    // applied explicitly below and need not follow this append order.
    append(dsa, 0);
    append(traversing, 1);
    append(cosmic, 2);

    for (auto& candidate : candidates) {
      auto const timing = timingDirection(*candidate.track,
                                          trackingGeometry,
                                          minTimingMeasurements_,
                                          minTimingDeltaChi2_);
      candidate.timingDirectionSign = timing.directionSign;
      candidate.timingMeasurements = timing.measurements;
      candidate.timingChi2 = timing.chi2;
      candidate.timingDeltaChi2 = timing.deltaChi2;
    }

    // Determine the common source side from all reconstructed hypotheses.
    // Its magnitude is never used as a selection, while the event-level sign
    // resolves the no-timing ambiguity consistently for both detector sides.
    int eventSourceSide = 1;
    double largestAbsOriginZ = 0.;
    unsigned int positiveOrigins = 0, negativeOrigins = 0;
    for (auto const& candidate : candidates) {
      if (!candidate.targetLineState.valid)
        continue;
      double const originZ = candidate.targetLineState.position.z();
      positiveOrigins += originZ > 0.;
      negativeOrigins += originZ < 0.;
      if (std::abs(originZ) > largestAbsOriginZ) {
        largestAbsOriginZ = std::abs(originZ);
        eventSourceSide = originZ < 0. ? -1 : 1;
      }
    }
    if (positiveOrigins != negativeOrigins)
      eventSourceSide = positiveOrigins > negativeOrigins ? 1 : -1;

    // Repeat the field-aware transport from the upstream fitted endpoint.
    // The preliminary inner-state result above is used only to infer the
    // common +/-z side and is never stored.
    for (auto& candidate : candidates) {
      int const directionSign = candidate.timingDirectionSign != 0
                                    ? candidate.timingDirectionSign
                                    : shiftDirectionSign(candidate.targetLineState.momentum, eventSourceSide);
      candidate.physicalDirectionSign = directionSign;
      auto const preRefitState = propagateToTargetLine(
          *candidate.track, vacuumPropagator, directionSign, &materialPropagator, eventSourceSide);
      candidate.preRefitPt = preRefitState.valid ? preRefitState.momentum.perp() : 0.;
      candidate.targetLineState = preRefitState;

      // Refit the same reconstructed hits in their measured time-of-flight
      // direction.  This changes the direction in which the Kalman filter
      // applies material updates; flipping a completed fit cannot do that.
      candidate.directionalRefitAttempted = true;
      auto const refit = directionalRefit(*candidate.track,
                                          directionSign,
                                          eventSourceSide,
                                          magneticField,
                                          muonRecHitBuilder,
                                          hitCloner,
                                          vacuumPropagator,
                                          materialPropagator,
                                          directionalRefitErrorRescale_,
                                          directionalRefitMaxHitChi2_);
      candidate.directionalRefitValid = refit.valid;
      candidate.directionalRefitHits = refit.hits;
      candidate.directionalRefitChi2 = refit.chi2;
      candidate.directionalRefitNdof = refit.ndof;
      candidate.directionalRefitUpstreamPt = refit.upstreamPt;
      if (refit.valid)
        candidate.targetLineState = refit.targetLineState;
    }

    // Do not select on reconstructed z.  Reject only invalid/empty fits and
    // trajectories which do not point back to the unbounded target line in
    // the transverse plane, or which lack the forward topology of particles
    // arriving along that line.  These remove mixed-hit cosmic candidates
    // without preferring either source side or any longitudinal location.
    candidates.erase(std::remove_if(candidates.begin(),
                                    candidates.end(),
                                    [this](Candidate const& candidate) {
                                      double const pt = candidate.targetLineState.momentum.perp();
                                      double const absEta =
                                          pt > 0. ? std::abs(std::asinh(candidate.targetLineState.momentum.z() / pt))
                                                  : std::numeric_limits<double>::infinity();
                                      return !candidate.targetLineState.valid ||
                                             candidate.track->numberOfValidHits() == 0 ||
                                             candidate.targetLineState.position.perp() > maxTargetLineDca_ ||
                                             absEta < minAbsEta_;
                                    }),
                     candidates.end());

    std::stable_sort(candidates.begin(), candidates.end(), [](Candidate const& first, Candidate const& second) {
      if (first.source != second.source)
        return first.source < second.source;
      if (first.track->numberOfValidHits() != second.track->numberOfValidHits())
        return first.track->numberOfValidHits() > second.track->numberOfValidHits();
      return first.track->normalizedChi2() < second.track->normalizedChi2();
    });

    auto duplicates = [this](Candidate const& first, Candidate const& second) {
      // A chamber/detId is much coarser than a hit.  Two nearby J/psi muons,
      // or the two legs of a through-going track, often cross the same
      // chambers at different local positions and must not be merged for
      // that reason alone.
      unsigned int sharedInputs = 0;
      for (auto const& firstHit : first.hitFingerprints) {
        bool sharesThisHit = false;
        for (auto const& secondHit : second.hitFingerprints)
          if (firstHit.detId == secondHit.detId &&
              std::hypot(firstHit.localX - secondHit.localX, firstHit.localY - secondHit.localY) < 0.1) {
            sharesThisHit = true;
            break;
          }
        sharedInputs += sharesThisHit;
      }
      auto const smallerRecHitSet = std::min(first.hitFingerprints.size(), second.hitFingerprints.size());
      bool const exactSharedInputs = smallerRecHitSet > 0 && sharedInputs >= minSharedDetIds_ &&
                                     static_cast<double>(sharedInputs) / smallerRecHitSet >= minSharedHitFraction_;
      bool const sameLine = directionAngle(*first.track, *second.track) < maxDuplicateAngle_ &&
                            symmetricLineDistance(*first.track, *second.track) < maxDuplicateLineDistance_;
      return exactSharedInputs || sameLine;
    };

    // Collapse complete connected components rather than applying an
    // order-dependent greedy veto.  This handles A<->B<->C duplicate chains
    // even when the two endpoint fits narrowly fail the direct comparison.
    std::vector<unsigned int> parent(candidates.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto findRoot = [&parent](unsigned int index) {
      while (parent[index] != index) {
        parent[index] = parent[parent[index]];
        index = parent[index];
      }
      return index;
    };
    for (unsigned int first = 0; first < candidates.size(); ++first)
      for (unsigned int second = first + 1; second < candidates.size(); ++second)
        if (duplicates(candidates[first], candidates[second])) {
          unsigned int const firstRoot = findRoot(first);
          unsigned int const secondRoot = findRoot(second);
          if (firstRoot != secondRoot)
            parent[secondRoot] = firstRoot;
        }

    std::vector<std::vector<unsigned int>> groups(candidates.size());
    for (unsigned int index = 0; index < candidates.size(); ++index)
      groups[findRoot(index)].push_back(index);
    auto betterRepresentative = [&candidates, eventSourceSide](unsigned int first, unsigned int second) {
      auto const& a = candidates[first];
      auto const& b = candidates[second];
      bool const aHasTiming = a.timingDirectionSign != 0;
      bool const bHasTiming = b.timingDirectionSign != 0;
      if (aHasTiming != bHasTiming)
        return aHasTiming;
      // Both seed orientations were fitted.  Prefer the fit whose native
      // propagation direction agrees with time of flight; merely flipping a
      // completed wrong-way fit would retain its biased material corrections
      // and is the source of the near-zero-pT solution.
      if (aHasTiming && (a.timingDirectionSign > 0) != (b.timingDirectionSign > 0))
        return a.timingDirectionSign > 0;
      if (aHasTiming && bHasTiming && a.timingChi2 != b.timingChi2)
        return a.timingChi2 < b.timingChi2;
      auto geometryPriority = [](int source) {
        // Traversing and cosmic fits use both detector legs and give much
        // better direction/origin resolution.  DSA remains available when no
        // such fit exists, but must not replace a better geometrical fit just
        // because its momentum magnitude is less biased.
        return source == 1 ? 0 : (source == 2 ? 1 : 2);
      };
      if (geometryPriority(a.source) != geometryPriority(b.source))
        return geometryPriority(a.source) < geometryPriority(b.source);
      // For a far source, a standalone leg in the endcap on the same side as
      // its propagated target-line PCA is the upstream leg.  It has crossed
      // less field and material than a downstream leg and gives a markedly
      // less ambiguous backward transport.  This uses reconstructed
      // geometry only and works symmetrically for +/-z sources.
      bool const aUpstream = eventSourceSide * a.track->innerPosition().z() > 0. &&
                             eventSourceSide * a.track->outerPosition().z() > 0.;
      bool const bUpstream = eventSourceSide * b.track->innerPosition().z() > 0. &&
                             eventSourceSide * b.track->outerPosition().z() > 0.;
      if (aUpstream != bUpstream)
        return aUpstream;
      if (a.track->hitPattern().muonStationsWithValidHits() != b.track->hitPattern().muonStationsWithValidHits())
        return a.track->hitPattern().muonStationsWithValidHits() > b.track->hitPattern().muonStationsWithValidHits();
      if (a.track->numberOfValidHits() != b.track->numberOfValidHits())
        return a.track->numberOfValidHits() > b.track->numberOfValidHits();
      return a.track->normalizedChi2() < b.track->normalizedChi2();
    };
    std::vector<Candidate const*> selected;
    std::vector<unsigned int> duplicateGroupSize;
    for (auto const& group : groups) {
      if (group.empty())
        continue;
      unsigned int representative = group.front();
      for (auto const member : group)
        if (betterRepresentative(member, representative))
          representative = member;
      selected.push_back(&candidates[representative]);
      duplicateGroupSize.push_back(group.size());
    }
    auto candidateDirectionSign = [eventSourceSide](Candidate const& candidate) {
      return candidate.physicalDirectionSign != 0
                 ? candidate.physicalDirectionSign
                 : static_cast<int>(shiftDirectionSign(candidate.targetLineState.momentum, eventSourceSide));
    };

    // Build an optional one-to-one MC association only after cleaning.  This
    // block cannot affect reconstruction or the duplicate decision, and an
    // absent collection (as in collision data) leaves every index at -1.
    std::vector<int> genPartIdx(selected.size(), -1);
    std::vector<float> genPartDeltaR(selected.size(), -1.f);
    if (genParticles.isValid()) {
      struct Match {
        double deltaR;
        unsigned int selectedIndex;
        unsigned int genIndex;
      };
      std::vector<Match> matches;
      for (unsigned int selectedIndex = 0; selectedIndex < selected.size(); ++selectedIndex)
        for (unsigned int genIndex = 0; genIndex < genParticles->size(); ++genIndex) {
          auto const& particle = (*genParticles)[genIndex];
          if (std::abs(particle.pdgId()) != 13 || particle.status() != 1)
            continue;
          double const deltaR = directionDeltaR(selected[selectedIndex]->targetLineState.momentum, particle);
          if (deltaR < maxGenDeltaR_)
            matches.push_back({deltaR, selectedIndex, genIndex});
        }
      // This is a diagnostic association, not part of reconstruction.  Allow
      // several reconstructed candidates to reference the same primary muon;
      // one-to-one assignment hid overlaps by leaving all but one at -1.
      for (auto const& match : matches)
        if (genPartIdx[match.selectedIndex] < 0 || match.deltaR < genPartDeltaR[match.selectedIndex]) {
          genPartIdx[match.selectedIndex] = match.genIndex;
          genPartDeltaR[match.selectedIndex] = match.deltaR;
        }
    }

    std::vector<float> pt, eta, phi, mass, p, px, py, pz, trackPt, innerPt, outerPt, upstreamPt, preRefitPt,
        directionalRefitUpstreamPt, directionalRefitChi2, directionalRefitNdof, ptError, etaError,
        phiError, vx, vy, vz, trackVx, trackVy, trackVz, dxy, dz, innerR, innerZ, outerR, outerZ, chi2, ndof,
        normalizedChi2, linePcaR, linePcaZ,
        chordLinePcaR, chordLinePcaZ, targetLinePath, timingChi2, timingDeltaChi2;
    std::vector<int> charge, source, sourceIndex, validHits, validMuonHits, muonStations, lostHits, directionFlipped,
        inferredSourceSide, chargeMatchesGen, timingDirectionSign, nTimingMeasurements, directionalRefitAttempted,
        directionalRefitValid, directionalRefitHits;
    for (unsigned int selectedIndex = 0; selectedIndex < selected.size(); ++selectedIndex) {
      auto const* candidate = selected[selectedIndex];
      auto const& track = *candidate->track;
      auto const& propagated = candidate->targetLineState;
      auto const [pcaR, pcaZ] = transverseLinePca(track);
      auto const [chordPcaR, chordPcaZ] = chordLinePca(track);
      // Cosmic-style fits do not determine which way along the fitted helix
      // the particle travelled.  Orient it from the inferred source side
      // toward CMS, supporting both +z and -z SHIFT locations without truth.
      double const sign = candidateDirectionSign(*candidate);
      bool const flip = sign < 0.;
      // The propagated state was already canonicalised before transport.
      double const storedPx = propagated.momentum.x();
      double const storedPy = propagated.momentum.y();
      double const storedPz = propagated.momentum.z();
      double const storedPt = std::hypot(storedPx, storedPy);
      double const storedP = propagated.momentum.mag();
      pt.push_back(storedPt);
      eta.push_back(storedPt > 0. ? std::asinh(storedPz / storedPt)
                                  : std::copysign(std::numeric_limits<float>::infinity(), storedPz));
      phi.push_back(std::atan2(storedPy, storedPx));
      mass.push_back(0.105658f);
      p.push_back(storedP);
      px.push_back(storedPx);
      py.push_back(storedPy);
      pz.push_back(storedPz);
      trackPt.push_back(track.pt());
      innerPt.push_back(track.innerMomentum().rho());
      outerPt.push_back(track.outerMomentum().rho());
      auto const endpointDelta = track.outerPosition() - track.innerPosition();
      bool const upstreamIsOuter = sign * track.innerMomentum().Dot(endpointDelta) < 0.;
      upstreamPt.push_back(upstreamIsOuter ? track.outerMomentum().rho() : track.innerMomentum().rho());
      preRefitPt.push_back(candidate->preRefitPt);
      directionalRefitUpstreamPt.push_back(candidate->directionalRefitUpstreamPt);
      directionalRefitAttempted.push_back(candidate->directionalRefitAttempted);
      directionalRefitValid.push_back(candidate->directionalRefitValid);
      directionalRefitHits.push_back(candidate->directionalRefitHits);
      directionalRefitChi2.push_back(candidate->directionalRefitChi2);
      directionalRefitNdof.push_back(candidate->directionalRefitNdof);
      ptError.push_back(track.ptError());
      etaError.push_back(track.etaError());
      phiError.push_back(track.phiError());
      // Momentum and charge are the simultaneous two-fold ambiguity of a
      // no-timing cosmic-style helix fit.  Reverse both to preserve curvature.
      int const storedCharge = sign * track.charge();
      charge.push_back(storedCharge);
      int chargeMatch = -1;
      if (genParticles.isValid() && genPartIdx[selectedIndex] >= 0) {
        int const genCharge = (*genParticles)[genPartIdx[selectedIndex]].pdgId() == 13 ? -1 : 1;
        chargeMatch = storedCharge == genCharge;
      }
      chargeMatchesGen.push_back(chargeMatch);
      directionFlipped.push_back(flip);
      inferredSourceSide.push_back(eventSourceSide);
      timingDirectionSign.push_back(candidate->timingDirectionSign);
      nTimingMeasurements.push_back(candidate->timingMeasurements);
      timingChi2.push_back(candidate->timingChi2);
      timingDeltaChi2.push_back(candidate->timingDeltaChi2);
      vx.push_back(propagated.position.x());
      vy.push_back(propagated.position.y());
      vz.push_back(propagated.position.z());
      trackVx.push_back(track.vx());
      trackVy.push_back(track.vy());
      trackVz.push_back(track.vz());
      dxy.push_back(track.dxy());
      dz.push_back(track.dz());
      innerR.push_back(track.innerPosition().rho());
      innerZ.push_back(track.innerPosition().z());
      outerR.push_back(track.outerPosition().rho());
      outerZ.push_back(track.outerPosition().z());
      linePcaR.push_back(pcaR);
      linePcaZ.push_back(pcaZ);
      chordLinePcaR.push_back(chordPcaR);
      chordLinePcaZ.push_back(chordPcaZ);
      targetLinePath.push_back(propagated.path);
      chi2.push_back(track.chi2());
      ndof.push_back(track.ndof());
      normalizedChi2.push_back(track.normalizedChi2());
      validHits.push_back(track.numberOfValidHits());
      validMuonHits.push_back(track.hitPattern().numberOfValidMuonHits());
      muonStations.push_back(track.hitPattern().muonStationsWithValidHits());
      lostHits.push_back(track.numberOfLostHits());
      source.push_back(candidate->source);
      sourceIndex.push_back(candidate->sourceIndex);
    }

    auto table = std::make_unique<nanoaod::FlatTable>(selected.size(), "ShiftMuon", false, false);
    table->addColumn<float>("pt", pt, "transverse momentum");
    table->addColumn<float>("eta", eta, "pseudorapidity");
    table->addColumn<float>("phi", phi, "azimuthal angle");
    table->addColumn<float>("mass", mass, "muon mass");
    table->addColumn<float>("p", p, "momentum magnitude");
    table->addColumn<float>("px", px, "momentum x component");
    table->addColumn<float>("py", py, "momentum y component");
    table->addColumn<float>("pz", pz, "momentum z component");
    table->addColumn<float>("trackPt", trackPt, "pT at the original CMSSW track reference state");
    table->addColumn<float>("innerPt", innerPt, "pT at the geometrically inner detector state");
    table->addColumn<float>("outerPt", outerPt, "pT at the geometrically outer detector state");
    table->addColumn<float>("upstreamPt", upstreamPt, "pT at the fitted endpoint nearest the inferred source side");
    table->addColumn<float>("preRefitPt", preRefitPt, "target-line pT before the timing-directed Kalman refit");
    table->addColumn<float>(
        "directionalRefitUpstreamPt", directionalRefitUpstreamPt, "pT at the source-facing refitted hit state");
    table->addColumn<int>("directionalRefitAttempted",
                          directionalRefitAttempted,
                          "1 when the timing/inferred-direction Kalman refit was attempted");
    table->addColumn<int>("directionalRefitValid",
                          directionalRefitValid,
                          "1 when final kinematics come from the timing/inferred-direction Kalman refit");
    table->addColumn<int>("directionalRefitHits", directionalRefitHits, "valid hits retained by the directional refit");
    table->addColumn<float>("directionalRefitChi2", directionalRefitChi2, "directional-refit trajectory chi2");
    table->addColumn<float>("directionalRefitNdof", directionalRefitNdof, "directional-refit trajectory ndof");
    table->addColumn<float>("ptErr", ptError, "transverse-momentum uncertainty");
    table->addColumn<float>("etaErr", etaError, "pseudorapidity uncertainty");
    table->addColumn<float>("phiErr", phiError, "azimuthal-angle uncertainty");
    table->addColumn<int>("charge", charge, "electric charge");
    table->addColumn<int>("chargeMatchesGen",
                          chargeMatchesGen,
                          "MC diagnostic: 1/0 if charge agrees/disagrees with matched GenPart, -1 on data/unmatched");
    table->addColumn<int>("directionFlipped",
                          directionFlipped,
                          "1 when momentum and charge were reversed to point from inferred source toward CMS");
    table->addColumn<int>(
        "inferredSourceSide", inferredSourceSide, "sign of reconstructed linePcaZ: -1=-z source, +1=+z source");
    table->addColumn<int>("timingDirectionSign",
                          timingDirectionSign,
                          "chosen sign relative to fitted momentum; 0 when timing is inconclusive");
    table->addColumn<int>("nTimingMeasurements", nTimingMeasurements, "number of timed DT/CSC segments");
    table->addColumn<float>("timingChi2", timingChi2, "chi2 of the preferred time-of-flight direction");
    table->addColumn<float>(
        "timingDeltaChi2", timingDeltaChi2, "chi2 separation between the two time-of-flight directions");
    table->addColumn<float>("vx", vx, "x at field-aware stepping-helix PCA to the target line");
    table->addColumn<float>("vy", vy, "y at field-aware stepping-helix PCA to the target line");
    table->addColumn<float>("vz", vz, "z at field-aware stepping-helix PCA to the target line");
    table->addColumn<float>("trackVx", trackVx, "original CMSSW track reference-point x");
    table->addColumn<float>("trackVy", trackVy, "original CMSSW track reference-point y");
    table->addColumn<float>("trackVz", trackVz, "original CMSSW track reference-point z");
    table->addColumn<float>("dxy", dxy, "transverse impact parameter relative to the origin");
    table->addColumn<float>("dz", dz, "longitudinal impact parameter relative to the origin");
    table->addColumn<float>("innerR", innerR, "inner-state cylindrical radius");
    table->addColumn<float>("innerZ", innerZ, "inner-state z");
    table->addColumn<float>("outerR", outerR, "outer-state cylindrical radius");
    table->addColumn<float>("outerZ", outerZ, "outer-state z");
    table->addColumn<float>("linePcaR", linePcaR, "straight-line transverse PCA radius");
    table->addColumn<float>("linePcaZ", linePcaZ, "z at straight-line transverse PCA");
    table->addColumn<float>("chordLinePcaR", chordLinePcaR, "transverse PCA radius from the inner-to-outer chord");
    table->addColumn<float>("chordLinePcaZ", chordLinePcaZ, "z at transverse PCA from the inner-to-outer chord");
    table->addColumn<float>(
        "targetLinePath", targetLinePath, "signed stepping-helix path length from the inner state to target-line PCA");
    table->addColumn<float>("chi2", chi2, "track chi2");
    table->addColumn<float>("ndof", ndof, "track fit degrees of freedom");
    table->addColumn<float>("normalizedChi2", normalizedChi2, "track chi2 divided by ndof");
    table->addColumn<int>("nValidHits", validHits, "number of valid track hits");
    table->addColumn<int>("nValidMuonHits", validMuonHits, "number of valid muon-system hits");
    table->addColumn<int>("nMuonStations", muonStations, "muon stations with valid hits");
    table->addColumn<int>("nLostHits", lostHits, "number of lost track hits");
    table->addColumn<int>("source", source, "0=DSA, 1=traversing, 2=cosmic");
    table->addColumn<int>("sourceIndex", sourceIndex, "index in the source track collection");
    table->addColumn<unsigned int>("duplicateGroupSize",
                                   duplicateGroupSize,
                                   "number of transitive input-track duplicates represented by this row");
    table->addColumn<int>("genPartIdx", genPartIdx, "index in GenPart, or -1 when unmatched or on data");
    table->addColumn<float>(
        "genPartDeltaR", genPartDeltaR, "direction-ambiguous deltaR to matched GenPart, or -1 when unmatched/on data");
    event.put(std::move(table));

    // Fit every cleaned pair directly from its retained source tracks.  The
    // resulting indices always refer to ShiftMuon rows and therefore do not
    // depend on keeping any of the input collections in NanoAOD.
    std::vector<int> vertexMuonIdx1, vertexMuonIdx2, vertexIsOS;
    std::vector<int> vertexDcaStatus, vertexKalmanAttempted, vertexKalmanValid, vertexUsesLineFallback,
        vertexSameGenMuon, vertexGenIsOS;
    std::vector<float> vertexVx, vertexVy, vertexVz, vertexVxError, vertexVyError, vertexVzError, vertexChi2,
        vertexNdof, vertexNormalizedChi2, vertexProbability, vertexMass, vertexPt, vertexEta, vertexPhi, vertexPz,
        vertexDca, vertexDcaX, vertexDcaY, vertexDcaZ, vertexOriginCompatibilityChi2,
        vertexOriginCompatibilityNormalizedChi2;
    constexpr double muonMass = 0.105658;

    struct PairChoice {
      unsigned int first;
      unsigned int second;
      StraightLineApproach approach;
      CommonLineVertex fit;
      double originChi2;
      double score;
    };
    std::vector<PairChoice> pairChoices;
    for (unsigned int first = 0; first < selected.size(); ++first) {
      for (unsigned int second = first + 1; second < selected.size(); ++second) {
        auto const& firstState = selected[first]->targetLineState;
        auto const& secondState = selected[second]->targetLineState;
        auto const lineApproach = straightLineApproach(firstState, secondState);
        if (!lineApproach.valid || lineApproach.distance > maxPairDca_)
          continue;

        int const firstCharge = candidateDirectionSign(*selected[first]) * selected[first]->track->charge();
        int const secondCharge = candidateDirectionSign(*selected[second]) * selected[second]->track->charge();
        if (requireOppositeSign_ && firstCharge == secondCharge)
          continue;

        auto const firstOrigin = firstState.position;
        auto const secondOrigin = secondState.position;
        // A physical SHIFT pair must point back to the same side.  This is a
        // reconstruction-only condition and prevents +z/-z averages near 0.
        if (firstOrigin.z() * secondOrigin.z() <= 0.)
          continue;
        double const deltaX = firstOrigin.x() - secondOrigin.x();
        double const deltaY = firstOrigin.y() - secondOrigin.y();
        double const deltaZ = firstOrigin.z() - secondOrigin.z();
        double const originChi2 =
            (deltaX * deltaX + deltaY * deltaY) / (2. * originTransverseResolution_ * originTransverseResolution_) +
            deltaZ * deltaZ / (2. * originZResolution_ * originZResolution_);
        constexpr double originNdof = 3.;
        if (originChi2 / originNdof > maxPairOriginNormalizedChi2_)
          continue;
        auto const fit = commonLineVertex(firstState,
                                          secondState,
                                          commonVertexLineResolution_,
                                          commonVertexBeamLineResolution_);
        if (!fit.valid)
          continue;
        double const score = originChi2 / originNdof +
                             std::pow(lineApproach.distance / commonVertexLineResolution_, 2) + fit.chi2 / fit.ndof;
        pairChoices.push_back({first, second, lineApproach, fit, originChi2, score});
      }
    }
    std::stable_sort(pairChoices.begin(), pairChoices.end(), [](PairChoice const& first, PairChoice const& second) {
      return first.score < second.score;
    });

    std::vector<bool> muonAlreadyUsed(selected.size(), false);
    unsigned int retainedVertices = 0;
    for (auto const& choice : pairChoices) {
      unsigned int const first = choice.first;
      unsigned int const second = choice.second;
      if (muonAlreadyUsed[first] || muonAlreadyUsed[second])
        continue;
      if (maxDimuonVertices_ > 0 && retainedVertices >= maxDimuonVertices_)
        break;
      muonAlreadyUsed[first] = true;
      muonAlreadyUsed[second] = true;
      ++retainedVertices;

      auto const& firstState = selected[first]->targetLineState;
      auto const& secondState = selected[second]->targetLineState;
      int const firstCharge = candidateDirectionSign(*selected[first]) * selected[first]->track->charge();
      int const secondCharge = candidateDirectionSign(*selected[second]) * selected[second]->track->charge();

      auto canonicalMomentum = [](Candidate const& candidate) {
        auto const& state = candidate.targetLineState;
        return std::array<double, 3>{state.momentum.x(), state.momentum.y(), state.momentum.z()};
      };
      auto const firstP = canonicalMomentum(*selected[first]);
      auto const secondP = canonicalMomentum(*selected[second]);
      double const pairPx = firstP[0] + secondP[0];
      double const pairPy = firstP[1] + secondP[1];
      double const pairPz = firstP[2] + secondP[2];
      double const firstEnergy = std::sqrt(firstState.momentum.mag2() + muonMass * muonMass);
      double const secondEnergy = std::sqrt(secondState.momentum.mag2() + muonMass * muonMass);
      double const mass2 =
          std::pow(firstEnergy + secondEnergy, 2) - pairPx * pairPx - pairPy * pairPy - pairPz * pairPz;
      double const pairPt = std::hypot(pairPx, pairPy);

      vertexMuonIdx1.push_back(first);
      vertexMuonIdx2.push_back(second);
      vertexIsOS.push_back(firstCharge != secondCharge);
      int sameGenMuon = -1;
      int genIsOS = -1;
      if (genParticles.isValid() && genPartIdx[first] >= 0 && genPartIdx[second] >= 0) {
        sameGenMuon = genPartIdx[first] == genPartIdx[second];
        int const firstGenPdgId = (*genParticles)[genPartIdx[first]].pdgId();
        int const secondGenPdgId = (*genParticles)[genPartIdx[second]].pdgId();
        genIsOS = firstGenPdgId * secondGenPdgId < 0;
      }
      vertexSameGenMuon.push_back(sameGenMuon);
      vertexGenIsOS.push_back(genIsOS);
      vertexKalmanAttempted.push_back(0);
      vertexKalmanValid.push_back(0);
      vertexUsesLineFallback.push_back(1);
      vertexVx.push_back(choice.fit.position.x());
      vertexVy.push_back(choice.fit.position.y());
      vertexVz.push_back(choice.fit.position.z());
      vertexVxError.push_back(choice.fit.error[0]);
      vertexVyError.push_back(choice.fit.error[1]);
      vertexVzError.push_back(choice.fit.error[2]);
      vertexChi2.push_back(choice.fit.chi2);
      vertexNdof.push_back(choice.fit.ndof);
      vertexNormalizedChi2.push_back(choice.fit.chi2 / choice.fit.ndof);
      vertexProbability.push_back(ChiSquaredProbability(choice.fit.chi2, choice.fit.ndof));
      vertexOriginCompatibilityChi2.push_back(choice.originChi2);
      vertexOriginCompatibilityNormalizedChi2.push_back(choice.originChi2 / 3.);
      vertexDcaStatus.push_back(choice.approach.valid);
      vertexDca.push_back(choice.approach.distance);
      vertexDcaX.push_back(choice.approach.midpoint.x());
      vertexDcaY.push_back(choice.approach.midpoint.y());
      vertexDcaZ.push_back(choice.approach.midpoint.z());
      vertexMass.push_back(std::sqrt(std::max(0., mass2)));
      vertexPt.push_back(pairPt);
      vertexPz.push_back(pairPz);
      vertexEta.push_back(pairPt > 0. ? std::asinh(pairPz / pairPt)
                                      : std::copysign(std::numeric_limits<float>::infinity(), pairPz));
      vertexPhi.push_back(std::atan2(pairPy, pairPx));
    }

    auto vertexTable = std::make_unique<nanoaod::FlatTable>(vertexMuonIdx1.size(), "ShiftDimuonVertex", false, false);
    vertexTable->addColumn<int>("muonIdx1", vertexMuonIdx1, "index of first muon in ShiftMuon");
    vertexTable->addColumn<int>("muonIdx2", vertexMuonIdx2, "index of second muon in ShiftMuon");
    vertexTable->addColumn<int>("isOS", vertexIsOS, "1 for an opposite-sign pair");
    vertexTable->addColumn<int>(
        "kalmanAttempted", vertexKalmanAttempted, "1 when the pair lies inside safe transient-track propagation range");
    vertexTable->addColumn<int>("sameGenMuon",
                                vertexSameGenMuon,
                                "MC diagnostic: 1 when both legs match the same GenPart, -1 on data/unmatched");
    vertexTable->addColumn<int>(
        "genIsOS", vertexGenIsOS, "MC diagnostic: generator pair is opposite-sign, -1 on data/unmatched");
    vertexTable->addColumn<int>("kalmanValid",
                                vertexKalmanValid,
                                "1 when the unconstrained far-vertex Kalman fit converged near its line seed");
    vertexTable->addColumn<int>("usesLineFallback",
                                vertexUsesLineFallback,
                                "1 when position comes from straight-line closest approach after Kalman failure");
    vertexTable->addColumn<float>("vx", vertexVx, "unbounded common-line fit vertex x");
    vertexTable->addColumn<float>("vy", vertexVy, "unbounded common-line fit vertex y");
    vertexTable->addColumn<float>("vz", vertexVz, "unbounded common-line fit vertex z");
    vertexTable->addColumn<float>("vxErr", vertexVxError, "common-line fit vertex x uncertainty");
    vertexTable->addColumn<float>("vyErr", vertexVyError, "common-line fit vertex y uncertainty");
    vertexTable->addColumn<float>("vzErr", vertexVzError, "common-line fit vertex z uncertainty");
    vertexTable->addColumn<float>("chi2", vertexChi2, "common-line vertex fit chi2");
    vertexTable->addColumn<float>("ndof", vertexNdof, "common-line vertex fit degrees of freedom");
    vertexTable->addColumn<float>("normalizedChi2", vertexNormalizedChi2, "common-line vertex chi2 divided by ndof");
    vertexTable->addColumn<float>("probability", vertexProbability, "common-line vertex fit probability");
    vertexTable->addColumn<float>("originCompatibilityChi2",
                                  vertexOriginCompatibilityChi2,
                                  "compatibility chi2 of the two independent ShiftMuon origins");
    vertexTable->addColumn<float>("originCompatibilityNormalizedChi2",
                                  vertexOriginCompatibilityNormalizedChi2,
                                  "independent-origin compatibility chi2 divided by three");
    vertexTable->addColumn<int>("dcaValid", vertexDcaStatus, "1 when the two-track DCA calculation succeeded");
    vertexTable->addColumn<float>("dca", vertexDca, "three-dimensional distance of closest approach");
    vertexTable->addColumn<float>("dcaX", vertexDcaX, "x of the two-track closest-approach crossing point");
    vertexTable->addColumn<float>("dcaY", vertexDcaY, "y of the two-track closest-approach crossing point");
    vertexTable->addColumn<float>("dcaZ", vertexDcaZ, "z of the two-track closest-approach crossing point");
    vertexTable->addColumn<float>("mass", vertexMass, "dimuon invariant mass using canonical SHIFT directions");
    vertexTable->addColumn<float>("pt", vertexPt, "dimuon transverse momentum");
    vertexTable->addColumn<float>("pz", vertexPz, "dimuon longitudinal momentum");
    vertexTable->addColumn<float>("eta", vertexEta, "dimuon pseudorapidity");
    vertexTable->addColumn<float>("phi", vertexPhi, "dimuon azimuthal angle");
    event.put(std::move(vertexTable), "ShiftDimuonVertex");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription description;
    description.add<edm::InputTag>("dsaTracks", edm::InputTag("displacedStandAloneMuons"));
    description.add<edm::InputTag>("cosmicTracks", edm::InputTag("shiftCosmicMuons"));
    description.add<edm::InputTag>("traversingTracks", edm::InputTag("shiftTraversingMuons"));
    description.add<edm::InputTag>("genParticles", edm::InputTag("finalGenParticles"));
    description.add<std::string>("muonRecHitBuilder", "MuonRecHitBuilder");
    description.add<double>("directionalRefitErrorRescale", 100.0);
    description.add<double>("directionalRefitMaxHitChi2", 100000.0);
    description.add<double>("minSharedHitFraction", 0.5);
    description.add<unsigned int>("minSharedDetIds", 2);
    description.add<double>("maxDuplicateAngle", 0.03);
    description.add<double>("maxDuplicateLineDistance", 30.0);
    description.add<unsigned int>("minTimingMeasurements", 2);
    description.add<double>("minTimingDeltaChi2", 4.0);
    description.add<double>("maxTargetLineDca", 200.0);
    description.add<double>("minAbsEta", 3.0);
    description.add<double>("originTransverseResolution", 100.0);
    description.add<double>("originZResolution", 2000.0);
    description.add<double>("commonVertexLineResolution", 100.0);
    description.add<double>("commonVertexBeamLineResolution", 100.0);
    description.add<double>("maxPairOriginNormalizedChi2", 9.0);
    description.add<double>("maxPairDca", 500.0);
    description.add<unsigned int>("maxDimuonVertices", 1);
    description.add<bool>("requireOppositeSign", true);
    description.add<double>("maxGenDeltaR", 0.5);
    descriptions.add("shiftMuonTable", description);
  }

private:
  edm::EDGetTokenT<reco::TrackCollection> dsaToken_;
  edm::EDGetTokenT<reco::TrackCollection> cosmicToken_;
  edm::EDGetTokenT<reco::TrackCollection> traversingToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> trackingGeometryToken_;
  edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> muonRecHitBuilderToken_;
  double directionalRefitErrorRescale_;
  double directionalRefitMaxHitChi2_;
  double minSharedHitFraction_;
  unsigned int minSharedDetIds_;
  double maxDuplicateAngle_;
  double maxDuplicateLineDistance_;
  unsigned int minTimingMeasurements_;
  double minTimingDeltaChi2_;
  double maxTargetLineDca_;
  double minAbsEta_;
  double originTransverseResolution_;
  double originZResolution_;
  double commonVertexLineResolution_;
  double commonVertexBeamLineResolution_;
  double maxPairOriginNormalizedChi2_;
  double maxPairDca_;
  unsigned int maxDimuonVertices_;
  bool requireOppositeSign_;
  double maxGenDeltaR_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonTableProducer);
