#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include <cmath>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

namespace {
  template <typename Collection>
  int collectionSize(edm::Handle<Collection> const& handle) {
    return handle.isValid() ? static_cast<int>(handle->size()) : -1;
  }

  template <typename Collection>
  int digiCount(edm::Handle<Collection> const& handle) {
    if (!handle.isValid())
      return -1;
    int result = 0;
    for (auto const& item : *handle)
      result += std::distance(item.second.first, item.second.second);
    return result;
  }

  int muonSimHitCount(edm::Handle<edm::PSimHitContainer> const& handle) {
    if (!handle.isValid())
      return -1;
    int result = 0;
    for (auto const& hit : *handle)
      result += std::abs(hit.particleType()) == 13;
    return result;
  }
}  // namespace

class ShiftMuonRecoDiagnosticsProducer : public edm::stream::EDProducer<> {
public:
  explicit ShiftMuonRecoDiagnosticsProducer(edm::ParameterSet const& parameters)
      : dtDigis_(consumes<DTDigiCollection>(parameters.getParameter<edm::InputTag>("dtDigis"))),
        cscStripDigis_(
            consumes<CSCStripDigiCollection>(parameters.getParameter<edm::InputTag>("cscStripDigis"))),
        cscWireDigis_(consumes<CSCWireDigiCollection>(parameters.getParameter<edm::InputTag>("cscWireDigis"))),
        rpcDigis_(consumes<RPCDigiCollection>(parameters.getParameter<edm::InputTag>("rpcDigis"))),
        gemDigis_(consumes<GEMDigiCollection>(parameters.getParameter<edm::InputTag>("gemDigis"))),
        dtRecHits_(consumes<DTRecHitCollection>(parameters.getParameter<edm::InputTag>("dtRecHits"))),
        dtSegments_(
            consumes<DTRecSegment4DCollection>(parameters.getParameter<edm::InputTag>("dtSegments"))),
        cscRecHits_(
            consumes<CSCRecHit2DCollection>(parameters.getParameter<edm::InputTag>("cscRecHits"))),
        cscSegments_(consumes<CSCSegmentCollection>(parameters.getParameter<edm::InputTag>("cscSegments"))),
        rpcRecHits_(consumes<RPCRecHitCollection>(parameters.getParameter<edm::InputTag>("rpcRecHits"))),
        gemRecHits_(consumes<GEMRecHitCollection>(parameters.getParameter<edm::InputTag>("gemRecHits"))),
        gemSegments_(consumes<GEMSegmentCollection>(parameters.getParameter<edm::InputTag>("gemSegments"))),
        dtSimHits_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("dtSimHits"))),
        cscSimHits_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("cscSimHits"))),
        rpcSimHits_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("rpcSimHits"))),
        gemSimHits_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("gemSimHits"))),
        generalTracks_(consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("generalTracks"))),
        dsaGlobalLinks_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("dsaGlobalLinks"))),
        cosmicGlobalLinks_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("cosmicGlobalLinks"))),
        traversingGlobalLinks_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("traversingGlobalLinks"))) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& event, edm::EventSetup const&) override {
    auto table = std::make_unique<nanoaod::FlatTable>(1, "ShiftRecoDiag", true, false);
    auto add = [&table](std::string const& name, int value, std::string const& documentation) {
      table->addColumn<int>(name, std::vector<int>{value}, documentation);
    };

    add("nDTSimHits", muonSimHitCount(event.getHandle(dtSimHits_)), "muon SimHits in DT sensitive volumes");
    add("nCSCSimHits", muonSimHitCount(event.getHandle(cscSimHits_)), "muon SimHits in CSC sensitive volumes");
    add("nRPCSimHits", muonSimHitCount(event.getHandle(rpcSimHits_)), "muon SimHits in RPC sensitive volumes");
    add("nGEMSimHits", muonSimHitCount(event.getHandle(gemSimHits_)), "muon SimHits in GEM sensitive volumes");
    add("nDTDigis", digiCount(event.getHandle(dtDigis_)), "DT digis, or -1 when unavailable");
    add("nCSCStripDigis", digiCount(event.getHandle(cscStripDigis_)), "CSC strip digis, or -1 when unavailable");
    add("nCSCWireDigis", digiCount(event.getHandle(cscWireDigis_)), "CSC wire digis, or -1 when unavailable");
    add("nRPCDigis", digiCount(event.getHandle(rpcDigis_)), "RPC digis, or -1 when unavailable");
    add("nGEMDigis", digiCount(event.getHandle(gemDigis_)), "GEM digis, or -1 when unavailable");
    add("nDTRecHits", collectionSize(event.getHandle(dtRecHits_)), "DT reconstructed hits");
    add("nDTSegments", collectionSize(event.getHandle(dtSegments_)), "DT 4D segments");
    add("nCSCRecHits", collectionSize(event.getHandle(cscRecHits_)), "CSC reconstructed 2D hits");
    add("nCSCSegments", collectionSize(event.getHandle(cscSegments_)), "CSC segments");
    add("nRPCRecHits", collectionSize(event.getHandle(rpcRecHits_)), "RPC reconstructed hits");
    add("nGEMRecHits", collectionSize(event.getHandle(gemRecHits_)), "GEM reconstructed hits");
    add("nGEMSegments", collectionSize(event.getHandle(gemSegments_)), "GEM segments");
    add("nGeneralTracks", collectionSize(event.getHandle(generalTracks_)), "generalTracks available to global matching");
    add("nDSATrackerMatches", collectionSize(event.getHandle(dsaGlobalLinks_)), "DSA tracker+muon links");
    add("nCosmicTrackerMatches", collectionSize(event.getHandle(cosmicGlobalLinks_)), "cosmic tracker+muon links");
    add("nTraversingTrackerMatches",
        collectionSize(event.getHandle(traversingGlobalLinks_)),
        "strict-traversing tracker+muon links");
    event.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription description;
    description.add<edm::InputTag>("dtDigis", edm::InputTag("muonDTDigis"));
    description.add<edm::InputTag>("cscStripDigis", edm::InputTag("muonCSCDigis", "MuonCSCStripDigi"));
    description.add<edm::InputTag>("cscWireDigis", edm::InputTag("muonCSCDigis", "MuonCSCWireDigi"));
    description.add<edm::InputTag>("rpcDigis", edm::InputTag("muonRPCDigis"));
    description.add<edm::InputTag>("gemDigis", edm::InputTag("muonGEMDigis"));
    description.add<edm::InputTag>("dtRecHits", edm::InputTag("dt1DRecHits"));
    description.add<edm::InputTag>("dtSegments", edm::InputTag("dt4DSegments"));
    description.add<edm::InputTag>("cscRecHits", edm::InputTag("csc2DRecHits"));
    description.add<edm::InputTag>("cscSegments", edm::InputTag("cscSegments"));
    description.add<edm::InputTag>("rpcRecHits", edm::InputTag("rpcRecHits"));
    description.add<edm::InputTag>("gemRecHits", edm::InputTag("gemRecHits"));
    description.add<edm::InputTag>("gemSegments", edm::InputTag("gemSegments"));
    description.add<edm::InputTag>("dtSimHits", edm::InputTag("g4SimHits", "MuonDTHits"));
    description.add<edm::InputTag>("cscSimHits", edm::InputTag("g4SimHits", "MuonCSCHits"));
    description.add<edm::InputTag>("rpcSimHits", edm::InputTag("g4SimHits", "MuonRPCHits"));
    description.add<edm::InputTag>("gemSimHits", edm::InputTag("g4SimHits", "MuonGEMHits"));
    description.add<edm::InputTag>("generalTracks", edm::InputTag("generalTracks"));
    description.add<edm::InputTag>("dsaGlobalLinks", edm::InputTag("shiftGlobalDSAMuons"));
    description.add<edm::InputTag>("cosmicGlobalLinks", edm::InputTag("shiftGlobalCosmicMuons"));
    description.add<edm::InputTag>("traversingGlobalLinks", edm::InputTag("shiftGlobalTraversingMuons"));
    descriptions.add("shiftMuonRecoDiagnostics", description);
  }

private:
  edm::EDGetTokenT<DTDigiCollection> dtDigis_;
  edm::EDGetTokenT<CSCStripDigiCollection> cscStripDigis_;
  edm::EDGetTokenT<CSCWireDigiCollection> cscWireDigis_;
  edm::EDGetTokenT<RPCDigiCollection> rpcDigis_;
  edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  edm::EDGetTokenT<DTRecHitCollection> dtRecHits_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegments_;
  edm::EDGetTokenT<CSCRecHit2DCollection> cscRecHits_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHits_;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<GEMSegmentCollection> gemSegments_;
  edm::EDGetTokenT<edm::PSimHitContainer> dtSimHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> cscSimHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> rpcSimHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> gemSimHits_;
  edm::EDGetTokenT<reco::TrackCollection> generalTracks_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> dsaGlobalLinks_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> cosmicGlobalLinks_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> traversingGlobalLinks_;
};

DEFINE_FWK_MODULE(ShiftMuonRecoDiagnosticsProducer);
