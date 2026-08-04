#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
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
        maxGenDeltaR_(parameters.getParameter<double>("maxGenDeltaR")) {
    produces<nanoaod::FlatTable>();
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
      std::sort(matches.begin(), matches.end(), [](Match const& first, Match const& second) {
        return first.deltaR < second.deltaR;
      });
      std::vector<bool> usedGen(genParticles->size(), false);
      for (auto const& match : matches)
        if (genPartIdx[match.selectedIndex] < 0 && !usedGen[match.genIndex]) {
          genPartIdx[match.selectedIndex] = match.genIndex;
          usedGen[match.genIndex] = true;
        }
    }

    std::vector<float> pt, eta, phi, mass, p, px, py, pz, ptError, etaError, phiError, vx, vy, vz, dxy,
        dz, innerR, innerZ, outerR, outerZ, chi2, ndof, normalizedChi2, linePcaR, linePcaZ;
    std::vector<int> charge, source, sourceIndex, validHits, validMuonHits, muonStations, lostHits;
    for (auto const* candidate : selected) {
      auto const& track = *candidate->track;
      auto const [pcaR, pcaZ] = transverseLinePca(track);
      pt.push_back(track.pt());
      eta.push_back(track.eta());
      phi.push_back(track.phi());
      mass.push_back(0.105658f);
      p.push_back(track.p());
      px.push_back(track.px());
      py.push_back(track.py());
      pz.push_back(track.pz());
      ptError.push_back(track.ptError());
      etaError.push_back(track.etaError());
      phiError.push_back(track.phiError());
      charge.push_back(track.charge());
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
    event.put(std::move(table));
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
  double minSharedHitFraction_;
  unsigned int minSharedDetIds_;
  double maxDuplicateAngle_;
  double maxDuplicateLineDistance_;
  double maxGenDeltaR_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonTableProducer);
