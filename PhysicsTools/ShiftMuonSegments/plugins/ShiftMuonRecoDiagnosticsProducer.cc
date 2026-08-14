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
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>
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

  struct CaloSummary {
    int allHits = -1;
    int signalMuonHits = -1;
    float signalMuonEnergy = -1.f;
    float signalMuonFirstTime = 0.f;
    float signalMuonLastTime = 0.f;
  };

  CaloSummary caloSummary(edm::Handle<edm::PCaloHitContainer> const& hits,
                          std::unordered_set<unsigned int> const& signalMuonTrackIds,
                          bool enabled) {
    if (!enabled)
      return {-2, -2, -2.f, 0.f, 0.f};
    if (!hits.isValid())
      return {};
    CaloSummary result;
    result.allHits = hits->size();
    result.signalMuonHits = 0;
    result.signalMuonEnergy = 0.f;
    result.signalMuonFirstTime = std::numeric_limits<float>::infinity();
    result.signalMuonLastTime = -std::numeric_limits<float>::infinity();
    for (auto const& hit : *hits) {
      if (!signalMuonTrackIds.count(hit.geantTrackId()))
        continue;
      ++result.signalMuonHits;
      result.signalMuonEnergy += hit.energy();
      result.signalMuonFirstTime = std::min(result.signalMuonFirstTime, static_cast<float>(hit.time()));
      result.signalMuonLastTime = std::max(result.signalMuonLastTime, static_cast<float>(hit.time()));
    }
    if (result.signalMuonHits == 0)
      result.signalMuonFirstTime = result.signalMuonLastTime = 0.f;
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
        hcalSimHits_(consumes<edm::PCaloHitContainer>(parameters.getParameter<edm::InputTag>("hcalSimHits"))),
        zdcSimHits_(consumes<edm::PCaloHitContainer>(parameters.getParameter<edm::InputTag>("zdcSimHits"))),
        simTracks_(consumes<edm::SimTrackContainer>(parameters.getParameter<edm::InputTag>("simTracks"))),
        generalTracks_(consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("generalTracks"))),
        dsaGlobalLinks_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("dsaGlobalLinks"))),
        cosmicGlobalLinks_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("cosmicGlobalLinks"))),
        traversingGlobalLinks_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("traversingGlobalLinks"))),
        enableDTMeasurement_(parameters.getParameter<bool>("enableDTMeasurement")),
        enableGEMMeasurement_(parameters.getParameter<bool>("enableGEMMeasurement")),
        trackerMode_(parameters.getParameter<int>("trackerMode")),
        enableHcalDiagnostics_(parameters.getParameter<bool>("enableHcalDiagnostics")),
        enableZDCDiagnostics_(parameters.getParameter<bool>("enableZDCDiagnostics")),
        dtNavigationMode_(parameters.getParameter<int>("dtNavigationMode")),
        recoVariantCode_(parameters.getParameter<int>("recoVariantCode")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& event, edm::EventSetup const&) override {
    auto table = std::make_unique<nanoaod::FlatTable>(1, "ShiftRecoDiag", true, false);
    auto add = [&table](std::string const& name, int value, std::string const& documentation) {
      table->addColumn<int>(name, std::vector<int>{value}, documentation);
    };
    auto addFloat = [&table](std::string const& name, float value, std::string const& documentation) {
      table->addColumn<float>(name, std::vector<float>{value}, documentation);
    };

    std::unordered_set<unsigned int> signalMuonTrackIds;
    auto const simTracks = event.getHandle(simTracks_);
    if (simTracks.isValid())
      for (auto const& track : *simTracks)
        if (std::abs(track.type()) == 13 && track.genpartIndex() >= 0)
          signalMuonTrackIds.insert(track.trackId());
    auto const hcal = caloSummary(event.getHandle(hcalSimHits_), signalMuonTrackIds, enableHcalDiagnostics_);
    auto const zdc = caloSummary(event.getHandle(zdcSimHits_), signalMuonTrackIds, enableZDCDiagnostics_);

    add("recoVariantCode", recoVariantCode_, "workflow reconstruction-variant code");
    add("enableDTMeasurement", enableDTMeasurement_, "DT measurements enabled in Shift muon reconstruction");
    add("dtNavigationMode", dtNavigationMode_, "DT navigation: 0=off, 1=Standard, 2=Direct");
    add("enableGEMMeasurement", enableGEMMeasurement_, "GEM measurements enabled in Shift DSA reconstruction");
    add("trackerMode", trackerMode_, "tracker combination: 0=off, 1=generalTracks, 2=forward P5 prototype");
    add("enableHcalDiagnostics", enableHcalDiagnostics_, "HCAL SimHit association enabled");
    add("enableZDCDiagnostics", enableZDCDiagnostics_, "ZDC SimHit association enabled");

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
    add("nHcalSimHits", hcal.allHits, "all HCAL SimHits; -2 when this diagnostic is disabled");
    add("nSignalMuonHcalSimHits", hcal.signalMuonHits, "primary signal-muon HCAL SimHits");
    addFloat("signalMuonHcalSimEnergy", hcal.signalMuonEnergy, "summed primary signal-muon HCAL SimHit energy");
    addFloat("signalMuonHcalFirstTime", hcal.signalMuonFirstTime, "earliest primary signal-muon HCAL SimHit time");
    addFloat("signalMuonHcalLastTime", hcal.signalMuonLastTime, "latest primary signal-muon HCAL SimHit time");
    add("nZDCSimHits", zdc.allHits, "all ZDC SimHits; -2 when this diagnostic is disabled");
    add("nSignalMuonZDCSimHits", zdc.signalMuonHits, "primary signal-muon ZDC SimHits");
    addFloat("signalMuonZDCSimEnergy", zdc.signalMuonEnergy, "summed primary signal-muon ZDC SimHit energy");
    addFloat("signalMuonZDCFirstTime", zdc.signalMuonFirstTime, "earliest primary signal-muon ZDC SimHit time");
    addFloat("signalMuonZDCLastTime", zdc.signalMuonLastTime, "latest primary signal-muon ZDC SimHit time");
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
    description.add<edm::InputTag>("hcalSimHits", edm::InputTag("g4SimHits", "HcalHits"));
    description.add<edm::InputTag>("zdcSimHits", edm::InputTag("g4SimHits", "ZDCHITS"));
    description.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
    description.add<edm::InputTag>("generalTracks", edm::InputTag("generalTracks"));
    description.add<edm::InputTag>("dsaGlobalLinks", edm::InputTag("shiftGlobalDSAMuons"));
    description.add<edm::InputTag>("cosmicGlobalLinks", edm::InputTag("shiftGlobalCosmicMuons"));
    description.add<edm::InputTag>("traversingGlobalLinks", edm::InputTag("shiftGlobalTraversingMuons"));
    description.add<bool>("enableDTMeasurement", true);
    description.add<bool>("enableGEMMeasurement", true);
    description.add<int>("trackerMode", 1);
    description.add<bool>("enableHcalDiagnostics", false);
    description.add<bool>("enableZDCDiagnostics", false);
    description.add<int>("dtNavigationMode", 1);
    description.add<int>("recoVariantCode", 0);
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
  edm::EDGetTokenT<edm::PCaloHitContainer> hcalSimHits_;
  edm::EDGetTokenT<edm::PCaloHitContainer> zdcSimHits_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTracks_;
  edm::EDGetTokenT<reco::TrackCollection> generalTracks_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> dsaGlobalLinks_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> cosmicGlobalLinks_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> traversingGlobalLinks_;
  bool enableDTMeasurement_;
  bool enableGEMMeasurement_;
  int trackerMode_;
  bool enableHcalDiagnostics_;
  bool enableZDCDiagnostics_;
  int dtNavigationMode_;
  int recoVariantCode_;
};

DEFINE_FWK_MODULE(ShiftMuonRecoDiagnosticsProducer);
