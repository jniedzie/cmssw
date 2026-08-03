#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <sstream>
#include <string>

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
        seeds_(consumes<L2MuonTrajectorySeedCollection>(p.getParameter<edm::InputTag>("dsaSeeds"))),
        tracks_(consumes<reco::TrackCollection>(p.getParameter<edm::InputTag>("dsaTracks"))),
        printDetails_(p.getParameter<bool>("printDetails")) {}

  void analyze(edm::Event const& event, edm::EventSetup const&) override {
    auto dtHits = event.getHandle(dtHits_);
    auto dtSegments = event.getHandle(dtSegments_);
    auto cscHits = event.getHandle(cscHits_);
    auto cscSegments = event.getHandle(cscSegments_);
    auto rpcHits = event.getHandle(rpcHits_);
    auto gemHits = event.getHandle(gemHits_);
    auto gemSegments = event.getHandle(gemSegments_);
    auto seeds = event.getHandle(seeds_);
    auto tracks = event.getHandle(tracks_);

    edm::LogVerbatim("ShiftMuonRecoDebug")
        << "[ShiftMuonRecoDebug][summary] run=" << event.id().run() << " lumi=" << event.luminosityBlock()
        << " event=" << event.id().event() << " dtRecHits=" << sizeOrMissing(dtHits)
        << " dtSegments=" << sizeOrMissing(dtSegments) << " cscRecHits=" << sizeOrMissing(cscHits)
        << " cscSegments=" << sizeOrMissing(cscSegments) << " rpcRecHits=" << sizeOrMissing(rpcHits)
        << " gemRecHits=" << sizeOrMissing(gemHits) << " gemSegments=" << sizeOrMissing(gemSegments)
        << " dsaSeeds=" << sizeOrMissing(seeds) << " dsaTracks=" << sizeOrMissing(tracks);

    if (!printDetails_)
      return;
    if (dtSegments.isValid())
      for (auto const& segment : *dtSegments) {
        auto id = segment.chamberId();
        auto direction = segment.localDirection();
        edm::LogVerbatim("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][DTsegment] event=" << event.id().event() << " rawId=" << id.rawId()
            << " wheel=" << id.wheel() << " station=" << id.station() << " sector=" << id.sector()
            << " nHits=" << segment.recHits().size() << " chi2=" << segment.chi2()
            << " direction=(" << direction.x() << ',' << direction.y() << ',' << direction.z() << ')';
      }
    if (cscSegments.isValid())
      for (auto const& segment : *cscSegments) {
        auto id = segment.cscDetId();
        auto direction = segment.localDirection();
        edm::LogVerbatim("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][CSCsegment] event=" << event.id().event() << " rawId=" << id.rawId()
            << " endcap=" << id.endcap() << " station=" << id.station() << " ring=" << id.ring()
            << " chamber=" << id.chamber() << " nHits=" << segment.nRecHits() << " chi2=" << segment.chi2()
            << " direction=(" << direction.x() << ',' << direction.y() << ',' << direction.z() << ')';
      }
    if (rpcHits.isValid())
      for (auto const& hit : *rpcHits) {
        auto id = hit.rpcId();
        edm::LogVerbatim("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][RPChit] event=" << event.id().event() << " rawId=" << id.rawId()
            << " region=" << id.region() << " station=" << id.station() << " sector=" << id.sector()
            << " layer=" << id.layer() << " bx=" << hit.BunchX() << " time=" << hit.time();
      }
    if (gemSegments.isValid())
      for (auto const& segment : *gemSegments) {
        auto id = segment.gemDetId();
        auto direction = segment.localDirection();
        edm::LogVerbatim("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][GEMsegment] event=" << event.id().event() << " rawId=" << id.rawId()
            << " region=" << id.region() << " station=" << id.station() << " ring=" << id.ring()
            << " chamber=" << id.chamber() << " nHits=" << segment.nRecHits() << " chi2=" << segment.chi2()
            << " direction=(" << direction.x() << ',' << direction.y() << ',' << direction.z() << ')';
      }
    if (tracks.isValid()) {
      unsigned index = 0;
      for (auto const& track : *tracks) {
        edm::LogVerbatim("ShiftMuonRecoDebug")
            << "[ShiftMuonRecoDebug][DSAtrack] event=" << event.id().event() << " index=" << index++
            << " pt=" << track.pt() << " eta=" << track.eta() << " phi=" << track.phi()
            << " innerR=" << track.innerPosition().rho() << " innerZ=" << track.innerPosition().z()
            << " outerR=" << track.outerPosition().rho() << " outerZ=" << track.outerPosition().z()
            << " validHits=" << track.numberOfValidHits() << " lostHits=" << track.numberOfLostHits()
            << " chi2=" << track.chi2() << " ndof=" << track.ndof();
      }
    }
  }

private:
  edm::EDGetTokenT<DTRecHitCollection> dtHits_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegments_;
  edm::EDGetTokenT<CSCRecHit2DCollection> cscHits_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcHits_;
  edm::EDGetTokenT<GEMRecHitCollection> gemHits_;
  edm::EDGetTokenT<GEMSegmentCollection> gemSegments_;
  edm::EDGetTokenT<L2MuonTrajectorySeedCollection> seeds_;
  edm::EDGetTokenT<reco::TrackCollection> tracks_;
  bool printDetails_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonSegmentsCounter);
