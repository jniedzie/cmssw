#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanSmoothedVertexChi2Estimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanTrackToTrackCovCalculator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/LinearizationPointFinders/interface/FsmwLinearizationPointFinder.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/LinearizedTrackStateFactory.h"
#include "RecoVertex/VertexTools/interface/SequentialVertexFitter.h"
#include "RecoVertex/VertexTools/interface/SequentialVertexSmoother.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

namespace {
  struct Candidate {
    reco::Track const* track;
    int source;
    unsigned int sourceIndex;
    std::unordered_set<uint32_t> hitDetIds;
  };

  std::unordered_set<uint32_t> validHitDetIds(reco::Track const& track) {
    std::unordered_set<uint32_t> result;
    for (auto hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit)
      if ((*hit)->isValid())
        result.insert((*hit)->geographicalId().rawId());
    return result;
  }

  double directionAngle(reco::Track const& first, reco::Track const& second) {
    auto const firstDirection = first.momentum().unit();
    auto const secondDirection = second.momentum().unit();
    double const cosine = std::clamp(std::abs(firstDirection.Dot(secondDirection)), 0., 1.);
    return std::acos(cosine);
  }

  double directionDeltaR(reco::Track const& track, reco::GenParticle const& particle) {
    auto deltaPhi = [](double first, double second) {
      double value = std::remainder(first - second, 2. * M_PI);
      return std::abs(value);
    };
    double const direct = std::hypot(track.eta() - particle.eta(), deltaPhi(track.phi(), particle.phi()));
    double const reverse =
        std::hypot(track.eta() + particle.eta(), deltaPhi(track.phi() + M_PI, particle.phi()));
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

  double shiftDirectionSign(reco::Track const& track) {
    // Orient the track from its reconstructed beam-axis crossing toward CMS.
    // A source at +z must have pz<0, while a source at -z must have pz>0.
    // The line PCA is reconstructed information and is available on data.
    double const sourceZ = transverseLinePca(track).second;
    return sourceZ * track.pz() > 0. ? -1. : 1.;
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
        transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
        magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
        minSharedHitFraction_(parameters.getParameter<double>("minSharedHitFraction")),
        minSharedDetIds_(parameters.getParameter<unsigned int>("minSharedDetIds")),
        maxDuplicateAngle_(parameters.getParameter<double>("maxDuplicateAngle")),
        maxDuplicateLineDistance_(parameters.getParameter<double>("maxDuplicateLineDistance")),
        maxGenDeltaR_(parameters.getParameter<double>("maxGenDeltaR")) {
    produces<nanoaod::FlatTable>();
    produces<nanoaod::FlatTable>("ShiftDimuonVertex");
  }

  void produce(edm::Event& event, edm::EventSetup const& setup) override {
    auto const dsa = event.getHandle(dsaToken_);
    auto const cosmic = event.getHandle(cosmicToken_);
    auto const traversing = event.getHandle(traversingToken_);
    auto const genParticles = event.getHandle(genParticlesToken_);

    std::vector<Candidate> candidates;
    auto append = [&candidates](auto const& handle, int source) {
      if (!handle.isValid())
        return;
      unsigned int index = 0;
      for (auto const& track : *handle) {
        candidates.push_back({&track, source, index++, validHitDetIds(track)});
      }
    };
    // Lower source values have precedence.  MC studies show that DSA gives
    // the best momentum estimate, while traversing tracks outperform the
    // ordinary cosmic fit in direction and production-line precision.
    append(dsa, 0);
    append(traversing, 1);
    append(cosmic, 2);

    std::stable_sort(candidates.begin(), candidates.end(), [](Candidate const& first, Candidate const& second) {
      if (first.source != second.source)
        return first.source < second.source;
      if (first.track->numberOfValidHits() != second.track->numberOfValidHits())
        return first.track->numberOfValidHits() > second.track->numberOfValidHits();
      return first.track->normalizedChi2() < second.track->normalizedChi2();
    });

    auto duplicates = [this](Candidate const& first, Candidate const& second) {
      unsigned int shared = 0;
      for (auto const detId : first.hitDetIds)
        shared += second.hitDetIds.count(detId);
      auto const smallerHitSet = std::min(first.hitDetIds.size(), second.hitDetIds.size());
      bool const sharedHits = smallerHitSet > 0 && shared >= minSharedDetIds_ &&
                              static_cast<double>(shared) / smallerHitSet >= minSharedHitFraction_;
      bool const sameLine = directionAngle(*first.track, *second.track) < maxDuplicateAngle_ &&
                            symmetricLineDistance(*first.track, *second.track) < maxDuplicateLineDistance_;
      return sharedHits || sameLine;
    };

    std::vector<Candidate const*> selected;
    for (auto const& candidate : candidates) {
      bool overlaps = false;
      for (auto const* retained : selected)
        if (duplicates(candidate, *retained)) {
          overlaps = true;
          break;
        }
      if (!overlaps)
        selected.push_back(&candidate);
    }

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
          double const deltaR = directionDeltaR(*selected[selectedIndex]->track, particle);
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

    std::vector<float> pt, eta, phi, mass, p, px, py, pz, ptError, etaError, phiError, vx, vy, vz, dxy,
        dz, innerR, innerZ, outerR, outerZ, chi2, ndof, normalizedChi2, linePcaR, linePcaZ;
    std::vector<int> charge, source, sourceIndex, validHits, validMuonHits, muonStations, lostHits,
        directionFlipped, inferredSourceSide;
    for (auto const* candidate : selected) {
      auto const& track = *candidate->track;
      auto const [pcaR, pcaZ] = transverseLinePca(track);
      // Cosmic-style fits do not determine which way along the fitted helix
      // the particle travelled.  Orient it from the inferred source side
      // toward CMS, supporting both +z and -z SHIFT locations without truth.
      double const sign = shiftDirectionSign(track);
      bool const flip = sign < 0.;
      double const storedPx = sign * track.px();
      double const storedPy = sign * track.py();
      double const storedPz = sign * track.pz();
      pt.push_back(track.pt());
      eta.push_back(sign * track.eta());
      phi.push_back(std::atan2(storedPy, storedPx));
      mass.push_back(0.105658f);
      p.push_back(track.p());
      px.push_back(storedPx);
      py.push_back(storedPy);
      pz.push_back(storedPz);
      ptError.push_back(track.ptError());
      etaError.push_back(track.etaError());
      phiError.push_back(track.phiError());
      // Momentum and charge are the simultaneous two-fold ambiguity of a
      // no-timing cosmic-style helix fit.  Reverse both to preserve curvature.
      charge.push_back(sign * track.charge());
      directionFlipped.push_back(flip);
      inferredSourceSide.push_back((pcaZ > 0.) - (pcaZ < 0.));
      vx.push_back(track.vx());
      vy.push_back(track.vy());
      vz.push_back(track.vz());
      dxy.push_back(track.dxy());
      dz.push_back(track.dz());
      innerR.push_back(track.innerPosition().rho());
      innerZ.push_back(track.innerPosition().z());
      outerR.push_back(track.outerPosition().rho());
      outerZ.push_back(track.outerPosition().z());
      linePcaR.push_back(pcaR);
      linePcaZ.push_back(pcaZ);
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
    table->addColumn<float>("ptErr", ptError, "transverse-momentum uncertainty");
    table->addColumn<float>("etaErr", etaError, "pseudorapidity uncertainty");
    table->addColumn<float>("phiErr", phiError, "azimuthal-angle uncertainty");
    table->addColumn<int>("charge", charge, "electric charge");
    table->addColumn<int>("directionFlipped", directionFlipped,
                          "1 when momentum and charge were reversed to point from inferred source toward CMS");
    table->addColumn<int>("inferredSourceSide", inferredSourceSide,
                          "sign of reconstructed linePcaZ: -1=-z source, +1=+z source");
    table->addColumn<float>("vx", vx, "track reference-point x");
    table->addColumn<float>("vy", vy, "track reference-point y");
    table->addColumn<float>("vz", vz, "track reference-point z");
    table->addColumn<float>("dxy", dxy, "transverse impact parameter relative to the origin");
    table->addColumn<float>("dz", dz, "longitudinal impact parameter relative to the origin");
    table->addColumn<float>("innerR", innerR, "inner-state cylindrical radius");
    table->addColumn<float>("innerZ", innerZ, "inner-state z");
    table->addColumn<float>("outerR", outerR, "outer-state cylindrical radius");
    table->addColumn<float>("outerZ", outerZ, "outer-state z");
    table->addColumn<float>("linePcaR", linePcaR, "straight-line transverse PCA radius");
    table->addColumn<float>("linePcaZ", linePcaZ, "z at straight-line transverse PCA");
    table->addColumn<float>("chi2", chi2, "track chi2");
    table->addColumn<float>("ndof", ndof, "track fit degrees of freedom");
    table->addColumn<float>("normalizedChi2", normalizedChi2, "track chi2 divided by ndof");
    table->addColumn<int>("nValidHits", validHits, "number of valid track hits");
    table->addColumn<int>("nValidMuonHits", validMuonHits, "number of valid muon-system hits");
    table->addColumn<int>("nMuonStations", muonStations, "muon stations with valid hits");
    table->addColumn<int>("nLostHits", lostHits, "number of lost track hits");
    table->addColumn<int>("source", source, "0=DSA, 1=traversing, 2=cosmic");
    table->addColumn<int>("sourceIndex", sourceIndex, "index in the source track collection");
    table->addColumn<int>("genPartIdx", genPartIdx, "index in GenPart, or -1 when unmatched or on data");
    table->addColumn<float>("genPartDeltaR", genPartDeltaR,
                            "direction-ambiguous deltaR to matched GenPart, or -1 when unmatched/on data");
    event.put(std::move(table));

    // Fit every cleaned pair directly from its retained source tracks.  The
    // resulting indices always refer to ShiftMuon rows and therefore do not
    // depend on keeping any of the input collections in NanoAOD.
    std::vector<int> vertexMuonIdx1, vertexMuonIdx2, vertexIsOS;
    std::vector<int> vertexDcaStatus;
    std::vector<float> vertexX, vertexY, vertexZ, vertexXError, vertexYError, vertexZError,
        vertexChi2, vertexNdof, vertexNormalizedChi2, vertexProbability, vertexMass, vertexPt, vertexEta, vertexPhi,
        vertexDca, vertexDcaX, vertexDcaY, vertexDcaZ;
    auto const& builder = setup.getData(transientTrackBuilderToken_);
    auto const& magneticField = setup.getData(magneticFieldToken_);
    edm::ParameterSet vertexParameters;
    vertexParameters.addParameter<double>("maxDistance", 0.01);
    vertexParameters.addParameter<int>("maxNbrOfIterations", 100);
    KalmanVertexTrackUpdator<5> vertexTrackUpdator;
    KalmanSmoothedVertexChi2Estimator<5> vertexChi2Estimator;
    KalmanTrackToTrackCovCalculator<5> covarianceCalculator;
    SequentialVertexSmoother<5> vertexSmoother(vertexTrackUpdator, vertexChi2Estimator, covarianceCalculator);
    SequentialVertexFitter<5> vertexFitter(vertexParameters,
                                           FsmwLinearizationPointFinder(20, -2., 0.4, 10.),
                                           KalmanVertexUpdator<5>(),
                                           vertexSmoother,
                                           LinearizedTrackStateFactory());
    // The standard fitter rejects vertices outside either the tracker or,
    // with its DSA option, |z|<960 cm muon-system bounds.  SHIFT production
    // is around |z|=14800 cm, so retain the Kalman algorithm but remove this
    // detector-size acceptance restriction.
    vertexFitter.setTrackerBounds(1.e6, 1.e6);
    constexpr double muonMass = 0.105658;
    for (unsigned int first = 0; first < selected.size(); ++first) {
      for (unsigned int second = first + 1; second < selected.size(); ++second) {
        std::vector<reco::TransientTrack> transientTracks{
            builder.build(selected[first]->track), builder.build(selected[second]->track)};
        TwoTrackMinimumDistance minimumDistance;
        FreeTrajectoryState firstState(GlobalPoint(selected[first]->track->vx(),
                                                   selected[first]->track->vy(),
                                                   selected[first]->track->vz()),
                                       GlobalVector(selected[first]->track->px(),
                                                    selected[first]->track->py(),
                                                    selected[first]->track->pz()),
                                       selected[first]->track->charge(),
                                       &magneticField);
        FreeTrajectoryState secondState(GlobalPoint(selected[second]->track->vx(),
                                                    selected[second]->track->vy(),
                                                    selected[second]->track->vz()),
                                        GlobalVector(selected[second]->track->px(),
                                                     selected[second]->track->py(),
                                                     selected[second]->track->pz()),
                                        selected[second]->track->charge(),
                                        &magneticField);
        bool const dcaIsValid = minimumDistance.calculate(firstState, secondState);
        auto const dcaCrossingPoint = dcaIsValid ? minimumDistance.crossingPoint() : GlobalPoint();
        double const fallbackZ = 0.5 * (transverseLinePca(*selected[first]->track).second +
                                        transverseLinePca(*selected[second]->track).second);
        GlobalPoint const linearizationPoint = dcaIsValid ? dcaCrossingPoint : GlobalPoint(0., 0., fallbackZ);
        // This overload uses the supplied point only for linearization; it is
        // not a positional prior or beam-spot constraint.
        TransientVertex const vertex = vertexFitter.vertex(transientTracks, linearizationPoint);
        if (!vertex.isValid())
          continue;
        reco::Vertex const persistentVertex(vertex);

        auto canonicalMomentum = [](reco::Track const& track) {
          double const sign = shiftDirectionSign(track);
          return std::array<double, 3>{sign * track.px(), sign * track.py(), sign * track.pz()};
        };
        auto const firstP = canonicalMomentum(*selected[first]->track);
        auto const secondP = canonicalMomentum(*selected[second]->track);
        double const pairPx = firstP[0] + secondP[0];
        double const pairPy = firstP[1] + secondP[1];
        double const pairPz = firstP[2] + secondP[2];
        double const firstEnergy = std::sqrt(selected[first]->track->p() * selected[first]->track->p() +
                                             muonMass * muonMass);
        double const secondEnergy = std::sqrt(selected[second]->track->p() * selected[second]->track->p() +
                                              muonMass * muonMass);
        double const mass2 = std::pow(firstEnergy + secondEnergy, 2) -
                             pairPx * pairPx - pairPy * pairPy - pairPz * pairPz;
        double const pairPt = std::hypot(pairPx, pairPy);

        vertexMuonIdx1.push_back(first);
        vertexMuonIdx2.push_back(second);
        int const firstCharge = shiftDirectionSign(*selected[first]->track) * selected[first]->track->charge();
        int const secondCharge = shiftDirectionSign(*selected[second]->track) * selected[second]->track->charge();
        vertexIsOS.push_back(firstCharge != secondCharge);
        vertexX.push_back(vertex.position().x());
        vertexY.push_back(vertex.position().y());
        vertexZ.push_back(vertex.position().z());
        vertexXError.push_back(persistentVertex.xError());
        vertexYError.push_back(persistentVertex.yError());
        vertexZError.push_back(persistentVertex.zError());
        vertexChi2.push_back(vertex.totalChiSquared());
        vertexNdof.push_back(vertex.degreesOfFreedom());
        vertexNormalizedChi2.push_back(vertex.degreesOfFreedom() > 0.
                                           ? vertex.totalChiSquared() / vertex.degreesOfFreedom()
                                           : std::numeric_limits<float>::infinity());
        vertexProbability.push_back(vertex.degreesOfFreedom() > 0.
                                        ? ChiSquaredProbability(vertex.totalChiSquared(), vertex.degreesOfFreedom())
                                        : 0.f);
        vertexDcaStatus.push_back(dcaIsValid);
        vertexDca.push_back(dcaIsValid ? minimumDistance.distance() : -1.f);
        vertexDcaX.push_back(dcaIsValid ? dcaCrossingPoint.x() : 0.f);
        vertexDcaY.push_back(dcaIsValid ? dcaCrossingPoint.y() : 0.f);
        vertexDcaZ.push_back(dcaIsValid ? dcaCrossingPoint.z() : 0.f);
        vertexMass.push_back(std::sqrt(std::max(0., mass2)));
        vertexPt.push_back(pairPt);
        vertexEta.push_back(pairPt > 0. ? std::asinh(pairPz / pairPt)
                                       : std::copysign(std::numeric_limits<float>::infinity(), pairPz));
        vertexPhi.push_back(std::atan2(pairPy, pairPx));
      }
    }

    auto vertexTable =
        std::make_unique<nanoaod::FlatTable>(vertexMuonIdx1.size(), "ShiftDimuonVertex", false, false);
    vertexTable->addColumn<int>("muonIdx1", vertexMuonIdx1, "index of first muon in ShiftMuon");
    vertexTable->addColumn<int>("muonIdx2", vertexMuonIdx2, "index of second muon in ShiftMuon");
    vertexTable->addColumn<int>("isOS", vertexIsOS, "1 for an opposite-sign pair");
    vertexTable->addColumn<float>("x", vertexX, "Kalman vertex x");
    vertexTable->addColumn<float>("y", vertexY, "Kalman vertex y");
    vertexTable->addColumn<float>("z", vertexZ, "Kalman vertex z");
    vertexTable->addColumn<float>("xErr", vertexXError, "Kalman vertex x uncertainty");
    vertexTable->addColumn<float>("yErr", vertexYError, "Kalman vertex y uncertainty");
    vertexTable->addColumn<float>("zErr", vertexZError, "Kalman vertex z uncertainty");
    vertexTable->addColumn<float>("chi2", vertexChi2, "Kalman vertex chi2");
    vertexTable->addColumn<float>("ndof", vertexNdof, "Kalman vertex degrees of freedom");
    vertexTable->addColumn<float>("normalizedChi2", vertexNormalizedChi2, "Kalman vertex chi2 divided by ndof");
    vertexTable->addColumn<float>("probability", vertexProbability, "Kalman vertex chi2 probability");
    vertexTable->addColumn<int>("dcaValid", vertexDcaStatus, "1 when the two-track DCA calculation succeeded");
    vertexTable->addColumn<float>("dca", vertexDca, "three-dimensional distance of closest approach");
    vertexTable->addColumn<float>("dcaX", vertexDcaX, "x of the two-track closest-approach crossing point");
    vertexTable->addColumn<float>("dcaY", vertexDcaY, "y of the two-track closest-approach crossing point");
    vertexTable->addColumn<float>("dcaZ", vertexDcaZ, "z of the two-track closest-approach crossing point");
    vertexTable->addColumn<float>("mass", vertexMass, "dimuon invariant mass using canonical SHIFT directions");
    vertexTable->addColumn<float>("pt", vertexPt, "dimuon transverse momentum");
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
    description.add<double>("minSharedHitFraction", 0.5);
    description.add<unsigned int>("minSharedDetIds", 2);
    description.add<double>("maxDuplicateAngle", 0.03);
    description.add<double>("maxDuplicateLineDistance", 30.0);
    description.add<double>("maxGenDeltaR", 0.5);
    descriptions.add("shiftMuonTable", description);
  }

private:
  edm::EDGetTokenT<reco::TrackCollection> dsaToken_;
  edm::EDGetTokenT<reco::TrackCollection> cosmicToken_;
  edm::EDGetTokenT<reco::TrackCollection> traversingToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  double minSharedHitFraction_;
  unsigned int minSharedDetIds_;
  double maxDuplicateAngle_;
  double maxDuplicateLineDistance_;
  double maxGenDeltaR_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonTableProducer);
