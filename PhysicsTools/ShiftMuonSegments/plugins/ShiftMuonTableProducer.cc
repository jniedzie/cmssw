#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

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

  reco::Track::Point transverseLinePcaPoint(reco::Track const& track) {
    auto const& position = track.innerPosition();
    auto const& momentum = track.innerMomentum();
    double const pt2 = momentum.x() * momentum.x() + momentum.y() * momentum.y();
    if (pt2 == 0.)
      return position;
    double const scale = -(position.x() * momentum.x() + position.y() * momentum.y()) / pt2;
    return position + scale * momentum;
  }

  double shiftDirectionSign(reco::Track const& track, int sourceSide) {
    return sourceSide * track.pz() > 0. ? -1. : 1.;
  }

  struct StraightLineApproach {
    bool valid;
    GlobalPoint midpoint;
    double distance;
  };

  StraightLineApproach straightLineApproach(reco::Track const& first, reco::Track const& second) {
    auto const firstPoint = first.innerPosition();
    auto const secondPoint = second.innerPosition();
    auto const firstDirection = first.innerMomentum().unit();
    auto const secondDirection = second.innerMomentum().unit();
    auto const separation = firstPoint - secondPoint;
    double const dot = firstDirection.Dot(secondDirection);
    double const denominator = 1. - dot * dot;
    if (denominator < 1.e-10) {
      auto const firstPca = transverseLinePcaPoint(first);
      auto const secondPca = transverseLinePcaPoint(second);
      GlobalPoint const midpoint(0.5 * (firstPca.x() + secondPca.x()),
                                 0.5 * (firstPca.y() + secondPca.y()),
                                 0.5 * (firstPca.z() + secondPca.z()));
      return {true, midpoint, (firstPca - secondPca).R()};
    }
    double const firstScale = (dot * secondDirection.Dot(separation) - firstDirection.Dot(separation)) /
                              denominator;
    double const secondScale = (secondDirection.Dot(separation) - dot * firstDirection.Dot(separation)) /
                               denominator;
    auto const firstClosest = firstPoint + firstScale * firstDirection;
    auto const secondClosest = secondPoint + secondScale * secondDirection;
    GlobalPoint const midpoint(0.5 * (firstClosest.x() + secondClosest.x()),
                               0.5 * (firstClosest.y() + secondClosest.y()),
                               0.5 * (firstClosest.z() + secondClosest.z()));
    double const distance = (firstClosest - secondClosest).R();
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
        minSharedHitFraction_(parameters.getParameter<double>("minSharedHitFraction")),
        minSharedDetIds_(parameters.getParameter<unsigned int>("minSharedDetIds")),
        maxDuplicateAngle_(parameters.getParameter<double>("maxDuplicateAngle")),
        maxDuplicateLineDistance_(parameters.getParameter<double>("maxDuplicateLineDistance")),
        minAbsOriginZ_(parameters.getParameter<double>("minAbsOriginZ")),
        originTransverseResolution_(parameters.getParameter<double>("originTransverseResolution")),
        originZResolution_(parameters.getParameter<double>("originZResolution")),
        requireOppositeSign_(parameters.getParameter<bool>("requireOppositeSign")),
        maxGenDeltaR_(parameters.getParameter<double>("maxGenDeltaR")) {
    produces<nanoaod::FlatTable>();
    produces<nanoaod::FlatTable>("ShiftDimuonVertex");
  }

  void produce(edm::Event& event, edm::EventSetup const&) override {
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

    // This is a SHIFT-specific topology requirement, not an MC requirement:
    // reject detector-edge PCA solutions before overlap removal so they do
    // not suppress a far-origin reconstruction of the same hits.
    candidates.erase(std::remove_if(candidates.begin(),
                                    candidates.end(),
                                    [this](Candidate const& candidate) {
                                      return std::abs(transverseLinePca(*candidate.track).second) < minAbsOriginZ_;
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

    // All primary muons in an event share the target side.  Use the most
    // displaced reconstructed origin to define that side once per event,
    // preventing a bad low-z leg from independently reversing eta and charge.
    int eventSourceSide = 1;
    double largestAbsOriginZ = 0.;
    unsigned int positiveOrigins = 0, negativeOrigins = 0;
    for (auto const* candidate : selected) {
      double const originZ = transverseLinePca(*candidate->track).second;
      positiveOrigins += originZ > 0.;
      negativeOrigins += originZ < 0.;
      if (std::abs(originZ) > largestAbsOriginZ) {
        largestAbsOriginZ = std::abs(originZ);
        eventSourceSide = originZ < 0. ? -1 : 1;
      }
    }
    if (positiveOrigins != negativeOrigins)
      eventSourceSide = positiveOrigins > negativeOrigins ? 1 : -1;

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

    std::vector<float> pt, eta, phi, mass, p, px, py, pz, ptError, etaError, phiError, vx, vy, vz,
        trackVx, trackVy, trackVz, dxy, dz, innerR, innerZ, outerR, outerZ, chi2, ndof, normalizedChi2,
        linePcaR, linePcaZ;
    std::vector<int> charge, source, sourceIndex, validHits, validMuonHits, muonStations, lostHits,
        directionFlipped, inferredSourceSide, chargeMatchesGen;
    for (unsigned int selectedIndex = 0; selectedIndex < selected.size(); ++selectedIndex) {
      auto const* candidate = selected[selectedIndex];
      auto const& track = *candidate->track;
      auto const [pcaR, pcaZ] = transverseLinePca(track);
      auto const pcaPoint = transverseLinePcaPoint(track);
      // Cosmic-style fits do not determine which way along the fitted helix
      // the particle travelled.  Orient it from the inferred source side
      // toward CMS, supporting both +z and -z SHIFT locations without truth.
      double const sign = shiftDirectionSign(track, eventSourceSide);
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
      vx.push_back(pcaPoint.x());
      vy.push_back(pcaPoint.y());
      vz.push_back(pcaPoint.z());
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
    table->addColumn<int>("chargeMatchesGen", chargeMatchesGen,
                          "MC diagnostic: 1/0 if charge agrees/disagrees with matched GenPart, -1 on data/unmatched");
    table->addColumn<int>("directionFlipped", directionFlipped,
                          "1 when momentum and charge were reversed to point from inferred source toward CMS");
    table->addColumn<int>("inferredSourceSide", inferredSourceSide,
                          "sign of reconstructed linePcaZ: -1=-z source, +1=+z source");
    table->addColumn<float>("vx", vx, "estimated origin x at straight-line transverse PCA");
    table->addColumn<float>("vy", vy, "estimated origin y at straight-line transverse PCA");
    table->addColumn<float>("vz", vz, "estimated origin z at straight-line transverse PCA");
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
    std::vector<int> vertexDcaStatus, vertexKalmanAttempted, vertexKalmanValid, vertexUsesLineFallback,
        vertexSameGenMuon, vertexGenIsOS;
    std::vector<float> vertexX, vertexY, vertexZ, vertexXError, vertexYError, vertexZError,
        vertexChi2, vertexNdof, vertexNormalizedChi2, vertexProbability, vertexMass, vertexPt, vertexEta, vertexPhi,
        vertexDca, vertexDcaX, vertexDcaY, vertexDcaZ;
    constexpr double muonMass = 0.105658;
    for (unsigned int first = 0; first < selected.size(); ++first) {
      for (unsigned int second = first + 1; second < selected.size(); ++second) {
        auto const lineApproach = straightLineApproach(*selected[first]->track, *selected[second]->track);
        if (!lineApproach.valid)
          continue;

        int const firstCharge = shiftDirectionSign(*selected[first]->track, eventSourceSide) *
                                selected[first]->track->charge();
        int const secondCharge = shiftDirectionSign(*selected[second]->track, eventSourceSide) *
                                 selected[second]->track->charge();
        if (requireOppositeSign_ && firstCharge == secondCharge)
          continue;

        auto canonicalMomentum = [eventSourceSide](reco::Track const& track) {
          double const sign = shiftDirectionSign(track, eventSourceSide);
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
        auto const firstOrigin = transverseLinePcaPoint(*selected[first]->track);
        auto const secondOrigin = transverseLinePcaPoint(*selected[second]->track);
        GlobalPoint const storedPosition(0.5 * (firstOrigin.x() + secondOrigin.x()),
                                         0.5 * (firstOrigin.y() + secondOrigin.y()),
                                         0.5 * (firstOrigin.z() + secondOrigin.z()));
        double const deltaX = firstOrigin.x() - secondOrigin.x();
        double const deltaY = firstOrigin.y() - secondOrigin.y();
        double const deltaZ = firstOrigin.z() - secondOrigin.z();
        double const geometricChi2 =
            (deltaX * deltaX + deltaY * deltaY) /
                (2. * originTransverseResolution_ * originTransverseResolution_) +
            deltaZ * deltaZ / (2. * originZResolution_ * originZResolution_);
        constexpr double geometricNdof = 3.;
        vertexKalmanAttempted.push_back(0);
        vertexKalmanValid.push_back(0);
        vertexUsesLineFallback.push_back(1);
        vertexX.push_back(storedPosition.x());
        vertexY.push_back(storedPosition.y());
        vertexZ.push_back(storedPosition.z());
        vertexXError.push_back(originTransverseResolution_ / std::sqrt(2.));
        vertexYError.push_back(originTransverseResolution_ / std::sqrt(2.));
        vertexZError.push_back(originZResolution_ / std::sqrt(2.));
        vertexChi2.push_back(geometricChi2);
        vertexNdof.push_back(geometricNdof);
        vertexNormalizedChi2.push_back(geometricChi2 / geometricNdof);
        vertexProbability.push_back(ChiSquaredProbability(geometricChi2, geometricNdof));
        vertexDcaStatus.push_back(lineApproach.valid);
        vertexDca.push_back(lineApproach.valid ? lineApproach.distance : -1.f);
        vertexDcaX.push_back(lineApproach.valid ? lineApproach.midpoint.x() : 0.f);
        vertexDcaY.push_back(lineApproach.valid ? lineApproach.midpoint.y() : 0.f);
        vertexDcaZ.push_back(lineApproach.valid ? lineApproach.midpoint.z() : 0.f);
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
    vertexTable->addColumn<int>("kalmanAttempted", vertexKalmanAttempted,
                                "1 when the pair lies inside safe transient-track propagation range");
    vertexTable->addColumn<int>("sameGenMuon", vertexSameGenMuon,
                                "MC diagnostic: 1 when both legs match the same GenPart, -1 on data/unmatched");
    vertexTable->addColumn<int>("genIsOS", vertexGenIsOS,
                                "MC diagnostic: generator pair is opposite-sign, -1 on data/unmatched");
    vertexTable->addColumn<int>("kalmanValid", vertexKalmanValid,
                                "1 when the unconstrained far-vertex Kalman fit converged near its line seed");
    vertexTable->addColumn<int>("usesLineFallback", vertexUsesLineFallback,
                                "1 when position comes from straight-line closest approach after Kalman failure");
    vertexTable->addColumn<float>("x", vertexX, "average of the two ShiftMuon transverse-PCA origins in x");
    vertexTable->addColumn<float>("y", vertexY, "average of the two ShiftMuon transverse-PCA origins in y");
    vertexTable->addColumn<float>("z", vertexZ, "average of the two ShiftMuon transverse-PCA origins in z");
    vertexTable->addColumn<float>("xErr", vertexXError, "configured geometric origin x uncertainty");
    vertexTable->addColumn<float>("yErr", vertexYError, "configured geometric origin y uncertainty");
    vertexTable->addColumn<float>("zErr", vertexZError, "configured geometric origin z uncertainty");
    vertexTable->addColumn<float>("chi2", vertexChi2, "compatibility chi2 of the two independent origin estimates");
    vertexTable->addColumn<float>("ndof", vertexNdof, "degrees of freedom of geometric origin compatibility");
    vertexTable->addColumn<float>("normalizedChi2", vertexNormalizedChi2, "geometric compatibility chi2 divided by ndof");
    vertexTable->addColumn<float>("probability", vertexProbability, "geometric origin compatibility probability");
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
    description.add<double>("minAbsOriginZ", 2000.0);
    description.add<double>("originTransverseResolution", 100.0);
    description.add<double>("originZResolution", 2000.0);
    description.add<bool>("requireOppositeSign", true);
    description.add<double>("maxGenDeltaR", 0.5);
    descriptions.add("shiftMuonTable", description);
  }

private:
  edm::EDGetTokenT<reco::TrackCollection> dsaToken_;
  edm::EDGetTokenT<reco::TrackCollection> cosmicToken_;
  edm::EDGetTokenT<reco::TrackCollection> traversingToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  double minSharedHitFraction_;
  unsigned int minSharedDetIds_;
  double maxDuplicateAngle_;
  double maxDuplicateLineDistance_;
  double minAbsOriginZ_;
  double originTransverseResolution_;
  double originZResolution_;
  bool requireOppositeSign_;
  double maxGenDeltaR_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonTableProducer);
