#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
  template <typename T>
  std::string sizeOrMissing(edm::Handle<T> const& handle) {
    return handle.isValid() ? std::to_string(handle->size()) : "MISSING";
  }
}

class ShiftMuonSegmentsCounter : public edm::one::EDAnalyzer<> {
public:
  explicit ShiftMuonSegmentsCounter(edm::ParameterSet const& p)
      : dtHits_(consumes<DTRecHitCollection>(p.getParameter<edm::InputTag>("dtRecHits"))),
        dtSegments_(consumes<DTRecSegment4DCollection>(p.getParameter<edm::InputTag>("dtSegments"))),
        cscHits_(consumes<CSCRecHit2DCollection>(p.getParameter<edm::InputTag>("cscRecHits"))),
        cscSegments_(consumes<CSCSegmentCollection>(p.getParameter<edm::InputTag>("cscSegments"))),
        rpcHits_(consumes<RPCRecHitCollection>(p.getParameter<edm::InputTag>("rpcRecHits"))),
        gemHits_(consumes<GEMRecHitCollection>(p.getParameter<edm::InputTag>("gemRecHits"))),
        gemSegments_(consumes<GEMSegmentCollection>(p.getParameter<edm::InputTag>("gemSegments"))),
        seeds_(consumes<TrajectorySeedCollection>(p.getParameter<edm::InputTag>("dsaSeeds"))),
        tracks_(consumes<reco::TrackCollection>(p.getParameter<edm::InputTag>("dsaTracks"))),
        cosmicTracks_(consumes<reco::TrackCollection>(p.getParameter<edm::InputTag>("cosmicTracks"))),
        traversingTracks_(consumes<reco::TrackCollection>(p.getParameter<edm::InputTag>("traversingTracks"))),
        cosmicTrackerTracks_(
            consumes<reco::TrackCollection>(p.getParameter<edm::InputTag>("cosmicTrackerTracks"))),
        simTracks_(consumes<edm::SimTrackContainer>(p.getParameter<edm::InputTag>("simTracks"))),
        simVertices_(consumes<edm::SimVertexContainer>(p.getParameter<edm::InputTag>("simVertices"))),
        dtSimHits_(consumes<edm::PSimHitContainer>(p.getParameter<edm::InputTag>("dtSimHits"))),
        cscSimHits_(consumes<edm::PSimHitContainer>(p.getParameter<edm::InputTag>("cscSimHits"))),
        rpcSimHits_(consumes<edm::PSimHitContainer>(p.getParameter<edm::InputTag>("rpcSimHits"))),
        gemSimHits_(consumes<edm::PSimHitContainer>(p.getParameter<edm::InputTag>("gemSimHits"))),
        magneticField_(esConsumes()),
        trackingGeometry_(esConsumes()),
        printDetails_(p.getParameter<bool>("printDetails")),
        printPropagationClosure_(p.getParameter<bool>("printPropagationClosure")) {}

  void analyze(edm::Event const& event, edm::EventSetup const& setup) override {
    auto dtHits = event.getHandle(dtHits_);
    auto dtSegments = event.getHandle(dtSegments_);
    auto cscHits = event.getHandle(cscHits_);
    auto cscSegments = event.getHandle(cscSegments_);
    auto rpcHits = event.getHandle(rpcHits_);
    auto gemHits = event.getHandle(gemHits_);
    auto gemSegments = event.getHandle(gemSegments_);
    auto seeds = event.getHandle(seeds_);
    auto tracks = event.getHandle(tracks_);
    auto cosmicTracks = event.getHandle(cosmicTracks_);
    auto traversingTracks = event.getHandle(traversingTracks_);
    auto cosmicTrackerTracks = event.getHandle(cosmicTrackerTracks_);
    auto simTracks = event.getHandle(simTracks_);
    auto simVertices = event.getHandle(simVertices_);
    auto dtSimHits = event.getHandle(dtSimHits_);
    auto cscSimHits = event.getHandle(cscSimHits_);
    auto rpcSimHits = event.getHandle(rpcSimHits_);
    auto gemSimHits = event.getHandle(gemSimHits_);

    if (printPropagationClosure_) {
      auto const& field = setup.getData(magneticField_);
      auto const& geometry = setup.getData(trackingGeometry_);
      SteppingHelixPropagator material(&field, alongMomentum);
      material.setMaterialMode(false);
      material.setUseMagVolumes(true);
      material.setUseMatVolumes(true);
      material.applyRadX0Correction(true);
      SteppingHelixPropagator vacuum(&field, alongMomentum);
      vacuum.setMaterialMode(true);
      vacuum.setUseMagVolumes(true);
      SteppingHelixPropagator backwardMaterial(&field, oppositeToMomentum);
      backwardMaterial.setMaterialMode(false);
      backwardMaterial.setUseMagVolumes(true);
      backwardMaterial.setUseMatVolumes(true);
      backwardMaterial.applyRadX0Correction(true);

      std::unordered_map<unsigned int, std::vector<PSimHit const*>> hitsByTrack;
      auto collectHits = [&hitsByTrack](auto const& handle) {
        if (!handle.isValid())
          return;
        for (auto const& hit : *handle)
          if (std::abs(hit.particleType()) == 13)
            hitsByTrack[hit.trackId()].push_back(&hit);
      };
      collectHits(dtSimHits);
      collectHits(cscSimHits);
      collectHits(rpcSimHits);
      collectHits(gemSimHits);

      unsigned int pairs = 0, materialValid = 0, vacuumValid = 0;
      double materialResidualSum = 0., materialResidual2Sum = 0.;
      double vacuumResidualSum = 0., truthChangeSum = 0.;
      unsigned int targetStates = 0, targetValid = 0;
      double targetResidualSum = 0., targetResidual2Sum = 0.;
      for (auto& [trackId, hits] : hitsByTrack) {
        if (!simTracks.isValid() || !simVertices.isValid())
          continue;
        auto const simTrack = std::find_if(simTracks->begin(), simTracks->end(), [trackId](auto const& track) {
          return track.trackId() == trackId && std::abs(track.type()) == 13;
        });
        // Closure is defined for the primary SHIFT muons. Secondary muons can
        // start inside CMS material and test a different transport problem.
        if (simTrack == simTracks->end() || simTrack->vertIndex() != 0 || simVertices->empty() ||
            !(simTrack->momentum().P() > 0.))
          continue;
        auto const& vertex = (*simVertices)[0].position();
        GlobalPoint const vertexPosition(vertex.x(), vertex.y(), vertex.z());
        GlobalVector const truthDirection(
            simTrack->momentum().px(), simTrack->momentum().py(), simTrack->momentum().pz());
        auto const pathFromVertex = [&geometry, &vertexPosition, &truthDirection](PSimHit const* hit) {
          auto const* det = geometry.idToDetUnit(DetId(hit->detUnitId()));
          if (!det)
            return std::numeric_limits<double>::infinity();
          return static_cast<double>(
              (det->surface().toGlobal(hit->entryPoint()) - vertexPosition).dot(truthDirection.unit()));
        };
        std::sort(hits.begin(), hits.end(), [&pathFromVertex](PSimHit const* first, PSimHit const* second) {
          return pathFromVertex(first) < pathFromVertex(second);
        });
        for (unsigned int index = 0; index + 1 < hits.size(); ++index) {
          auto const& first = *hits[index];
          auto const& second = *hits[index + 1];
          auto const* firstDet = geometry.idToDetUnit(DetId(first.detUnitId()));
          auto const* secondDet = geometry.idToDetUnit(DetId(second.detUnitId()));
          if (!firstDet || !secondDet || !(second.pabs() > 0.))
            continue;
          GlobalPoint const position = firstDet->surface().toGlobal(first.entryPoint());
          GlobalVector const momentum = firstDet->surface().toGlobal(first.momentumAtEntry());
          int const charge = first.particleType() > 0 ? -1 : 1;
          FreeTrajectoryState const state(GlobalTrajectoryParameters(position, momentum, charge, &field));
          ++pairs;
          truthChangeSum += (first.pabs() - second.pabs()) / second.pabs();
          auto const materialState = material.propagate(state, secondDet->surface());
          if (materialState.isValid()) {
            double const residual = (materialState.globalMomentum().mag() - second.pabs()) / second.pabs();
            ++materialValid;
            materialResidualSum += residual;
            materialResidual2Sum += residual * residual;
          }
          auto const vacuumState = vacuum.propagate(state, secondDet->surface());
          if (vacuumState.isValid()) {
            ++vacuumValid;
            vacuumResidualSum += (vacuumState.globalMomentum().mag() - second.pabs()) / second.pabs();
          }
        }

        if (hits.empty())
          continue;
        auto const& first = *hits.front();
        auto const* firstDet = geometry.idToDetUnit(DetId(first.detUnitId()));
        if (!firstDet)
          continue;
        GlobalPoint const position = firstDet->surface().toGlobal(first.entryPoint());
        GlobalVector const momentum = firstDet->surface().toGlobal(first.momentumAtEntry());
        int const charge = simTrack->type() > 0 ? -1 : 1;
        FreeTrajectoryState const state(GlobalTrajectoryParameters(position, momentum, charge, &field));
        int const sourceSide = (*simVertices)[0].position().z() < 0. ? -1 : 1;
        auto const boundary = Plane::build(GlobalPoint(0., 0., sourceSide * 1100.), Surface::RotationType());
        ++targetStates;
        auto const propagated = backwardMaterial.propagate(state, *boundary);
        if (propagated.isValid()) {
          double const residual =
              (propagated.globalMomentum().mag() - simTrack->momentum().P()) / simTrack->momentum().P();
          ++targetValid;
          targetResidualSum += residual;
          targetResidual2Sum += residual * residual;
        }
      }
      edm::LogPrint("ShiftMuonPropagationClosure")
          << "[ShiftMuonPropagationClosure] run=" << event.id().run() << " event=" << event.id().event()
          << " pairs=" << pairs << " materialValid=" << materialValid
          << " materialResidualSum=" << materialResidualSum
          << " materialResidual2Sum=" << materialResidual2Sum << " vacuumValid=" << vacuumValid
          << " vacuumResidualSum=" << vacuumResidualSum << " truthChangeSum=" << truthChangeSum
          << " targetStates=" << targetStates << " targetValid=" << targetValid
          << " targetResidualSum=" << targetResidualSum << " targetResidual2Sum=" << targetResidual2Sum;
    }

    edm::LogPrint("ShiftMuonRecoDebug")
        << "[ShiftMuonRecoDebug][summary] run=" << event.id().run() << " lumi=" << event.luminosityBlock()
        << " event=" << event.id().event() << " dtRecHits=" << sizeOrMissing(dtHits)
        << " dtSegments=" << sizeOrMissing(dtSegments) << " cscRecHits=" << sizeOrMissing(cscHits)
        << " cscSegments=" << sizeOrMissing(cscSegments) << " rpcRecHits=" << sizeOrMissing(rpcHits)
        << " gemRecHits=" << sizeOrMissing(gemHits) << " gemSegments=" << sizeOrMissing(gemSegments)
        << " dsaSeeds=" << sizeOrMissing(seeds) << " dsaTracks=" << sizeOrMissing(tracks)
        << " cosmicTracks=" << sizeOrMissing(cosmicTracks)
        << " traversingTracks=" << sizeOrMissing(traversingTracks)
        << " cosmicTrackerTracks=" << sizeOrMissing(cosmicTrackerTracks)
        << " simTracks=" << sizeOrMissing(simTracks) << " simVertices=" << sizeOrMissing(simVertices)
        << " dtSimHits=" << sizeOrMissing(dtSimHits) << " cscSimHits=" << sizeOrMissing(cscSimHits)
        << " rpcSimHits=" << sizeOrMissing(rpcSimHits) << " gemSimHits=" << sizeOrMissing(gemSimHits);

    if (!printDetails_)
      return;

    struct SimHitRange {
      PSimHit const* first = nullptr;
      PSimHit const* last = nullptr;
      unsigned int firstDetector = 4;
      unsigned int lastDetector = 4;
    };
    constexpr std::array<char const*, 5> detectorNames{{"DT", "CSC", "RPC", "GEM", "none"}};
    std::unordered_map<unsigned int, std::array<unsigned int, 4>> muonHitCounts;
    std::unordered_map<unsigned int, SimHitRange> muonHitRanges;
    std::unordered_map<unsigned int, SimHitRange> precisionHitRanges;
    auto updateRange = [](SimHitRange& range, PSimHit const& hit, unsigned int detectorIndex) {
      if (!range.first || hit.timeOfFlight() < range.first->timeOfFlight()) {
        range.first = &hit;
        range.firstDetector = detectorIndex;
      }
      if (!range.last || hit.timeOfFlight() > range.last->timeOfFlight()) {
        range.last = &hit;
        range.lastDetector = detectorIndex;
      }
    };
    auto countSimHits = [&muonHitCounts, &muonHitRanges, &precisionHitRanges, &updateRange](
                            auto const& handle, unsigned int detectorIndex) {
      if (!handle.isValid())
        return;
      for (auto const& hit : *handle) {
        ++muonHitCounts[hit.trackId()][detectorIndex];
        updateRange(muonHitRanges[hit.trackId()], hit, detectorIndex);
        if (detectorIndex < 2)
          updateRange(precisionHitRanges[hit.trackId()], hit, detectorIndex);
      }
    };
    countSimHits(dtSimHits, 0);
    countSimHits(cscSimHits, 1);
    countSimHits(rpcSimHits, 2);
    countSimHits(gemSimHits, 3);
    if (simTracks.isValid()) {
      for (auto const& simTrack : *simTracks) {
        if (std::abs(simTrack.type()) != 13)
          continue;
        auto const rangeIt = muonHitRanges.find(simTrack.trackId());
        SimHitRange const emptyRange;
        auto const& range = rangeIt != muonHitRanges.end() ? rangeIt->second : emptyRange;
        PSimHit const* firstMuonHit = range.first;
        PSimHit const* lastMuonHit = range.last;
        char const* firstMuonDetector = detectorNames[range.firstDetector];
        char const* lastMuonDetector = detectorNames[range.lastDetector];
        auto const precisionRangeIt = precisionHitRanges.find(simTrack.trackId());
        auto const& precisionRange =
            precisionRangeIt != precisionHitRanges.end() ? precisionRangeIt->second : emptyRange;
        PSimHit const* firstPrecisionHit = precisionRange.first;
        PSimHit const* lastPrecisionHit = precisionRange.last;
        double vertexZ = -999999.;
        if (simVertices.isValid() && simTrack.vertIndex() >= 0 &&
            static_cast<std::size_t>(simTrack.vertIndex()) < simVertices->size())
          vertexZ = (*simVertices)[simTrack.vertIndex()].position().z();
        auto const counts = muonHitCounts[simTrack.trackId()];
        double const simMomentum = simTrack.momentum().P();
        double const firstHitMomentum = firstMuonHit ? firstMuonHit->pabs() : -1.;
        double const lastHitMomentum = lastMuonHit ? lastMuonHit->pabs() : -1.;
        double const firstPrecisionHitMomentum = firstPrecisionHit ? firstPrecisionHit->pabs() : -1.;
        double const lastPrecisionHitMomentum = lastPrecisionHit ? lastPrecisionHit->pabs() : -1.;
        double const boundaryMomentum = simTrack.crossedBoundary() ? simTrack.getMomentumAtBoundary().P() : -1.;
        edm::LogPrint("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][SimMuon] event=" << event.id().event()
            << " trackId=" << simTrack.trackId() << " type=" << simTrack.type()
            << " vertIndex=" << simTrack.vertIndex() << " vertexZ=" << vertexZ
            << " pt=" << simTrack.momentum().pt() << " eta=" << simTrack.momentum().eta()
            << " phi=" << simTrack.momentum().phi() << " dtHits=" << counts[0]
            << " cscHits=" << counts[1] << " rpcHits=" << counts[2] << " gemHits=" << counts[3]
            << " simP=" << simMomentum << " crossedBoundary=" << simTrack.crossedBoundary()
            << " boundaryP=" << boundaryMomentum << " firstHitDetector=" << firstMuonDetector
            << " firstHitTof=" << (firstMuonHit ? firstMuonHit->timeOfFlight() : -1.)
            << " firstHitP=" << firstHitMomentum << " lastHitDetector=" << lastMuonDetector
            << " lastHitTof=" << (lastMuonHit ? lastMuonHit->timeOfFlight() : -1.)
            << " lastHitP=" << lastHitMomentum
            << " firstPrecisionHitDetector=" << detectorNames[precisionRange.firstDetector]
            << " firstPrecisionHitP=" << firstPrecisionHitMomentum
            << " lastPrecisionHitDetector=" << detectorNames[precisionRange.lastDetector]
            << " lastPrecisionHitP=" << lastPrecisionHitMomentum
            << " fractionalLossToFirstHit="
            << (simMomentum > 0. && firstHitMomentum >= 0. ? (simMomentum - firstHitMomentum) / simMomentum : -1.)
            << " fractionalLossAcrossMuonHits="
            << (firstHitMomentum > 0. && lastHitMomentum >= 0.
                    ? (firstHitMomentum - lastHitMomentum) / firstHitMomentum
                    : -1.)
            << " fractionalLossAcrossPrecisionHits="
            << (firstPrecisionHitMomentum > 0. && lastPrecisionHitMomentum >= 0.
                    ? (firstPrecisionHitMomentum - lastPrecisionHitMomentum) / firstPrecisionHitMomentum
                    : -1.);
      }
    }
    if (dtSegments.isValid())
      for (auto const& segment : *dtSegments) {
        auto id = segment.chamberId();
        auto direction = segment.localDirection();
        edm::LogPrint("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][DTsegment] event=" << event.id().event() << " rawId=" << id.rawId()
            << " wheel=" << id.wheel() << " station=" << id.station() << " sector=" << id.sector()
            << " nHits=" << segment.recHits().size() << " chi2=" << segment.chi2()
            << " direction=(" << direction.x() << ',' << direction.y() << ',' << direction.z() << ')';
      }
    if (cscSegments.isValid())
      for (auto const& segment : *cscSegments) {
        auto id = segment.cscDetId();
        auto direction = segment.localDirection();
        edm::LogPrint("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][CSCsegment] event=" << event.id().event() << " rawId=" << id.rawId()
            << " endcap=" << id.endcap() << " station=" << id.station() << " ring=" << id.ring()
            << " chamber=" << id.chamber() << " nHits=" << segment.nRecHits() << " chi2=" << segment.chi2()
            << " direction=(" << direction.x() << ',' << direction.y() << ',' << direction.z() << ')';
      }
    if (rpcHits.isValid())
      for (auto const& hit : *rpcHits) {
        auto id = hit.rpcId();
        edm::LogPrint("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][RPChit] event=" << event.id().event() << " rawId=" << id.rawId()
            << " region=" << id.region() << " station=" << id.station() << " sector=" << id.sector()
            << " layer=" << id.layer() << " bx=" << hit.BunchX() << " time=" << hit.time();
      }
    if (gemSegments.isValid())
      for (auto const& segment : *gemSegments) {
        auto id = segment.gemDetId();
        auto direction = segment.localDirection();
        edm::LogPrint("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][GEMsegment] event=" << event.id().event() << " rawId=" << id.rawId()
            << " region=" << id.region() << " station=" << id.station() << " ring=" << id.ring()
            << " chamber=" << id.chamber() << " nHits=" << segment.nRecHits() << " chi2=" << segment.chi2()
            << " direction=(" << direction.x() << ',' << direction.y() << ',' << direction.z() << ')';
      }
    if (seeds.isValid()) {
      unsigned index = 0;
      for (auto const& seed : *seeds) {
        edm::LogPrint("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][DSAseed] event=" << event.id().event() << " index=" << index++
            << " nHits=" << seed.nHits() << " detId=" << seed.startingState().detId()
            << " direction=" << static_cast<int>(seed.direction());
      }
    }
    if (tracks.isValid()) {
      unsigned index = 0;
      for (auto const& track : *tracks) {
        auto const& innerPosition = track.innerPosition();
        auto const& innerMomentum = track.innerMomentum();
        double const transverseMomentum2 = innerMomentum.x() * innerMomentum.x() +
                                           innerMomentum.y() * innerMomentum.y();
        double linePcaZ = innerPosition.z();
        double linePcaR = innerPosition.rho();
        if (transverseMomentum2 > 0.) {
          double const scale = -(innerPosition.x() * innerMomentum.x() +
                                 innerPosition.y() * innerMomentum.y()) /
                               transverseMomentum2;
          double const linePcaX = innerPosition.x() + scale * innerMomentum.x();
          double const linePcaY = innerPosition.y() + scale * innerMomentum.y();
          linePcaZ += scale * innerMomentum.z();
          linePcaR = std::hypot(linePcaX, linePcaY);
        }
        edm::LogPrint("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][DSAtrack] event=" << event.id().event() << " index=" << index++
            << " pt=" << track.pt() << " eta=" << track.eta() << " phi=" << track.phi()
            << " vx=" << track.vx() << " vy=" << track.vy() << " vz=" << track.vz()
            << " linePcaR=" << linePcaR << " linePcaZ=" << linePcaZ
            << " innerR=" << track.innerPosition().rho() << " innerZ=" << track.innerPosition().z()
            << " outerR=" << track.outerPosition().rho() << " outerZ=" << track.outerPosition().z()
            << " validHits=" << track.numberOfValidHits() << " lostHits=" << track.numberOfLostHits()
            << " chi2=" << track.chi2() << " ndof=" << track.ndof();
      }
    }
    auto printTrackCollection = [&event](auto const& handle, char const* collection) {
      if (!handle.isValid())
        return;
      unsigned int index = 0;
      for (auto const& track : *handle) {
        edm::LogPrint("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][" << collection << "] event=" << event.id().event()
            << " index=" << index++ << " pt=" << track.pt() << " eta=" << track.eta()
            << " phi=" << track.phi() << " vx=" << track.vx() << " vy=" << track.vy()
            << " vz=" << track.vz() << " innerR=" << track.innerPosition().rho()
            << " innerZ=" << track.innerPosition().z() << " outerR=" << track.outerPosition().rho()
            << " outerZ=" << track.outerPosition().z() << " validHits=" << track.numberOfValidHits()
            << " lostHits=" << track.numberOfLostHits() << " chi2=" << track.chi2()
            << " ndof=" << track.ndof();
      }
    };
    printTrackCollection(cosmicTracks, "CosmicTrack");
    printTrackCollection(traversingTracks, "TraversingTrack");
    printTrackCollection(cosmicTrackerTracks, "TrackerTrack");
  }

private:
  edm::EDGetTokenT<DTRecHitCollection> dtHits_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegments_;
  edm::EDGetTokenT<CSCRecHit2DCollection> cscHits_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcHits_;
  edm::EDGetTokenT<GEMRecHitCollection> gemHits_;
  edm::EDGetTokenT<GEMSegmentCollection> gemSegments_;
  edm::EDGetTokenT<TrajectorySeedCollection> seeds_;
  edm::EDGetTokenT<reco::TrackCollection> tracks_;
  edm::EDGetTokenT<reco::TrackCollection> cosmicTracks_;
  edm::EDGetTokenT<reco::TrackCollection> traversingTracks_;
  edm::EDGetTokenT<reco::TrackCollection> cosmicTrackerTracks_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTracks_;
  edm::EDGetTokenT<edm::SimVertexContainer> simVertices_;
  edm::EDGetTokenT<edm::PSimHitContainer> dtSimHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> cscSimHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> rpcSimHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> gemSimHits_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticField_;
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> trackingGeometry_;
  bool printDetails_;
  bool printPropagationClosure_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonSegmentsCounter);
