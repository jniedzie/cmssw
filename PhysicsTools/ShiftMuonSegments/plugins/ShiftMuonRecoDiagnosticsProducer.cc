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
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
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
#include <cctype>
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

  struct RecHitEnergySummary {
    int aboveThreshold = -1;
    float positiveEnergy = -1.f;
    float maximumEnergy = -1.f;
    float maximumTime = 0.f;
  };

  template <typename Collection>
  RecHitEnergySummary recHitEnergySummary(edm::Handle<Collection> const& hits, double threshold) {
    if (!hits.isValid())
      return {};
    RecHitEnergySummary result{0, 0.f, -std::numeric_limits<float>::infinity(), 0.f};
    for (auto const& hit : *hits) {
      float const energy = hit.energy();
      if (energy > 0.f)
        result.positiveEnergy += energy;
      result.aboveThreshold += energy >= threshold;
      if (energy > result.maximumEnergy) {
        result.maximumEnergy = energy;
        result.maximumTime = hit.time();
      }
    }
    if (hits->empty())
      result.maximumEnergy = -1.f;
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
        trackerSeeds_(consumes<TrajectorySeedCollection>(parameters.getParameter<edm::InputTag>("trackerSeeds"))),
        trackerTrackCandidates_(consumes<TrackCandidateCollection>(
            parameters.getParameter<edm::InputTag>("trackerTrackCandidates"))),
        trackerRawTracks_(
            consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("trackerRawTracks"))),
        trackerSelectedTracks_(
            consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("trackerSelectedTracks"))),
        trackerLHCTrackCandidates_(consumes<TrackCandidateCollection>(
            parameters.getParameter<edm::InputTag>("trackerLHCTrackCandidates"))),
        trackerLHCTracks_(
            consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("trackerLHCTracks"))),
        dsaGlobalLinks_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("dsaGlobalLinks"))),
        cosmicGlobalLinks_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("cosmicGlobalLinks"))),
        traversingGlobalLinks_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("traversingGlobalLinks"))),
        pixelRecHits_(consumes<SiPixelRecHitCollection>(parameters.getParameter<edm::InputTag>("pixelRecHits"))),
        stripMatchedRecHits_(consumes<SiStripMatchedRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripMatchedRecHits"))),
        stripRphiRecHits_(consumes<SiStripRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripRphiRecHits"))),
        stripRphiUnmatchedRecHits_(consumes<SiStripRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripRphiUnmatchedRecHits"))),
        stripStereoRecHits_(consumes<SiStripRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripStereoRecHits"))),
        stripStereoUnmatchedRecHits_(consumes<SiStripRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripStereoUnmatchedRecHits"))),
        hbheDigis_(consumes<HBHEDigiCollection>(parameters.getParameter<edm::InputTag>("hbheDigis"))),
        hfDigis_(consumes<HFDigiCollection>(parameters.getParameter<edm::InputTag>("hfDigis"))),
        hoDigis_(consumes<HODigiCollection>(parameters.getParameter<edm::InputTag>("hoDigis"))),
        zdcDigis_(consumes<ZDCDigiCollection>(parameters.getParameter<edm::InputTag>("zdcDigis"))),
        hbheQIE11Digis_(consumes<QIE11DigiCollection>(parameters.getParameter<edm::InputTag>("hbheQIE11Digis"))),
        hfQIE10Digis_(consumes<QIE10DigiCollection>(parameters.getParameter<edm::InputTag>("hfQIE10Digis"))),
        ecalBarrelRecHits_(consumes<EcalRecHitCollection>(parameters.getParameter<edm::InputTag>("ecalBarrelRecHits"))),
        ecalEndcapRecHits_(consumes<EcalRecHitCollection>(parameters.getParameter<edm::InputTag>("ecalEndcapRecHits"))),
        hbheRecHits_(consumes<HBHERecHitCollection>(parameters.getParameter<edm::InputTag>("hbheRecHits"))),
        hfRecHits_(consumes<HFRecHitCollection>(parameters.getParameter<edm::InputTag>("hfRecHits"))),
        hoRecHits_(consumes<HORecHitCollection>(parameters.getParameter<edm::InputTag>("hoRecHits"))),
        zdcRecHits_(consumes<ZDCRecHitCollection>(parameters.getParameter<edm::InputTag>("zdcRecHits"))),
        bcm1fSimHits_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("bcm1fSimHits"))),
        bhmSimHits_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("bhmSimHits"))),
        pltSimHits_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("pltSimHits"))),
        enableDTMeasurement_(parameters.getParameter<bool>("enableDTMeasurement")),
        enableGEMMeasurement_(parameters.getParameter<bool>("enableGEMMeasurement")),
        trackerMode_(parameters.getParameter<int>("trackerMode")),
        enableHcalDiagnostics_(parameters.getParameter<bool>("enableHcalDiagnostics")),
        enableZDCDiagnostics_(parameters.getParameter<bool>("enableZDCDiagnostics")),
        dtNavigationMode_(parameters.getParameter<int>("dtNavigationMode")) {
    for (auto const& tag : parameters.getParameter<std::vector<edm::InputTag>>("trackerSimHits"))
      trackerSimHits_.push_back(consumes<edm::PSimHitContainer>(tag));
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& event, edm::EventSetup const&) override {
    // Preserve the established NanoAOD schema: singleton table name provides
    // the ``ShiftRecoDiag_`` prefix and columns retain their lower-camel names.
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
    std::unordered_map<unsigned int, unsigned int> signalTrackerHitCounts;
    std::unordered_map<unsigned int, unsigned int> signalMuonSystemHitCounts;
    for (auto const& token : trackerSimHits_) {
      auto const hits = event.getHandle(token);
      if (!hits.isValid())
        continue;
      for (auto const& hit : *hits)
        if (signalMuonTrackIds.count(hit.trackId()))
          ++signalTrackerHitCounts[hit.trackId()];
    }
    auto countMuonSystemHits = [&](auto const& hits) {
      if (!hits.isValid())
        return;
      for (auto const& hit : *hits)
        if (signalMuonTrackIds.count(hit.trackId()))
          ++signalMuonSystemHitCounts[hit.trackId()];
    };
    countMuonSystemHits(event.getHandle(dtSimHits_));
    countMuonSystemHits(event.getHandle(cscSimHits_));
    countMuonSystemHits(event.getHandle(rpcSimHits_));
    countMuonSystemHits(event.getHandle(gemSimHits_));
    int signalWithTrackerHits = 0, signalWithMuonSystemHits = 0, signalWithBoth = 0;
    for (auto const trackId : signalMuonTrackIds) {
      bool const hasTracker = signalTrackerHitCounts[trackId] > 0;
      bool const hasMuonSystem = signalMuonSystemHitCounts[trackId] > 0;
      signalWithTrackerHits += hasTracker;
      signalWithMuonSystemHits += hasMuonSystem;
      signalWithBoth += hasTracker && hasMuonSystem;
    }
    auto const hcal = caloSummary(event.getHandle(hcalSimHits_), signalMuonTrackIds, enableHcalDiagnostics_);
    auto const zdc = caloSummary(event.getHandle(zdcSimHits_), signalMuonTrackIds, enableZDCDiagnostics_);
    auto const ecalBarrelRecHits = event.getHandle(ecalBarrelRecHits_);
    auto const ecalEndcapRecHits = event.getHandle(ecalEndcapRecHits_);
    auto const hbheRecHits = event.getHandle(hbheRecHits_);
    auto const hfRecHits = event.getHandle(hfRecHits_);
    auto const hoRecHits = event.getHandle(hoRecHits_);
    auto const zdcRecHits = event.getHandle(zdcRecHits_);
    auto const ebReco = recHitEnergySummary(ecalBarrelRecHits, 0.05);
    auto const eeReco = recHitEnergySummary(ecalEndcapRecHits, 0.05);
    auto const hbheReco = recHitEnergySummary(hbheRecHits, 0.1);
    auto const hfReco = recHitEnergySummary(hfRecHits, 0.1);
    auto const hoReco = recHitEnergySummary(hoRecHits, 0.1);
    auto const zdcReco = recHitEnergySummary(zdcRecHits, 0.01);

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
    add("nTrackerSeeds", collectionSize(event.getHandle(trackerSeeds_)), "dedicated cosmic tracker seeds");
    add("nTrackerTrackCandidates",
        collectionSize(event.getHandle(trackerTrackCandidates_)),
        "dedicated cosmic tracker CKF candidates");
    add("nTrackerRawTracks", collectionSize(event.getHandle(trackerRawTracks_)), "tracks before the cosmic selector");
    add("nTrackerSelectedTracks",
        collectionSize(event.getHandle(trackerSelectedTracks_)),
        "selected dedicated cosmic tracker tracks");
    add("nTrackerLHCTrackCandidates",
        collectionSize(event.getHandle(trackerLHCTrackCandidates_)),
        "dedicated P5 CKF candidates using the LHC navigation school");
    add("nTrackerLHCTracks",
        collectionSize(event.getHandle(trackerLHCTracks_)),
        "dedicated P5 tracks using the LHC navigation school");
    add("nSignalMuonSimTracks", signalMuonTrackIds.size(), "primary signal-muon SimTracks");
    add("nSignalMuonWithTrackerSimHits", signalWithTrackerHits, "primary signal muons with tracker SimHits");
    add("nSignalMuonWithMuonSystemSimHits",
        signalWithMuonSystemHits,
        "primary signal muons with DT, CSC, RPC, or GEM SimHits");
    add("nSignalMuonWithTrackerAndMuonSystemSimHits",
        signalWithBoth,
        "primary signal muons with both tracker and muon-system SimHits on the same SimTrack");
    add("nDSATrackerMatches", collectionSize(event.getHandle(dsaGlobalLinks_)), "DSA tracker+muon links");
    add("nCosmicTrackerMatches", collectionSize(event.getHandle(cosmicGlobalLinks_)), "cosmic tracker+muon links");
    add("nTraversingTrackerMatches",
        collectionSize(event.getHandle(traversingGlobalLinks_)),
        "strict-traversing tracker+muon links");
    add("nPixelRecHits", collectionSize(event.getHandle(pixelRecHits_)), "pixel rechits available for direct association");
    add("nStripMatchedRecHits", collectionSize(event.getHandle(stripMatchedRecHits_)), "matched strip rechits available for direct association");
    add("nStripRphiRecHits", collectionSize(event.getHandle(stripRphiRecHits_)), "rphi strip rechits available for direct association");
    add("nStripRphiUnmatchedRecHits", collectionSize(event.getHandle(stripRphiUnmatchedRecHits_)), "unmatched rphi strip rechits available for direct association");
    add("nStripStereoRecHits", collectionSize(event.getHandle(stripStereoRecHits_)), "stereo strip rechits available for direct association");
    add("nStripStereoUnmatchedRecHits", collectionSize(event.getHandle(stripStereoUnmatchedRecHits_)), "unmatched stereo strip rechits available for direct association");
    add("nHBHEDigis", collectionSize(event.getHandle(hbheDigis_)), "legacy HBHE digi frames");
    add("nHFDigis", collectionSize(event.getHandle(hfDigis_)), "legacy HF digi frames");
    add("nHODigis", collectionSize(event.getHandle(hoDigis_)), "HO digi frames");
    add("nZDCDigis", collectionSize(event.getHandle(zdcDigis_)), "legacy ZDC digi frames");
    add("nHBHEQIE11Digis", collectionSize(event.getHandle(hbheQIE11Digis_)), "Run-3 HBHE QIE11 digi frames");
    add("nHFQIE10Digis", collectionSize(event.getHandle(hfQIE10Digis_)), "Run-3 HF QIE10 digi frames");
    add("nEcalBarrelRecHits", collectionSize(ecalBarrelRecHits), "ECAL barrel rechits");
    add("nEcalEndcapRecHits", collectionSize(ecalEndcapRecHits), "ECAL endcap rechits");
    add("nHBHERecHits", collectionSize(hbheRecHits), "full HBHE rechits");
    add("nHFRecHits", collectionSize(hfRecHits), "full HF rechits");
    add("nHORecHits", collectionSize(hoRecHits), "full HO rechits");
    add("nZDCRecHits", collectionSize(zdcRecHits), "ZDC rechits, or -1 when unavailable");
    auto addRecoSummary = [&add, &addFloat](std::string const& prefix,
                                           RecHitEnergySummary const& summary,
                                           std::string const& threshold) {
      add("n" + prefix + "RecHitsAboveThreshold",
          summary.aboveThreshold,
          prefix + " rechits above " + threshold + " GeV");
      addFloat(prefix + "PositiveRecEnergy", summary.positiveEnergy, "summed positive " + prefix + " rechit energy");
      addFloat(prefix + "MaxRecEnergy", summary.maximumEnergy, "maximum " + prefix + " rechit energy");
      addFloat(prefix + "MaxRecTime", summary.maximumTime, "time of the maximum-energy " + prefix + " rechit");
    };
    addRecoSummary("EcalBarrel", ebReco, "0.05");
    addRecoSummary("EcalEndcap", eeReco, "0.05");
    addRecoSummary("HBHE", hbheReco, "0.1");
    addRecoSummary("HF", hfReco, "0.1");
    addRecoSummary("HO", hoReco, "0.1");
    addRecoSummary("ZDC", zdcReco, "0.01");
    add("nBCM1FSimHits", muonSimHitCount(event.getHandle(bcm1fSimHits_)), "muon SimHits in BCM1F");
    add("nBHMSimHits", muonSimHitCount(event.getHandle(bhmSimHits_)), "muon SimHits in BHM");
    add("nPLTSimHits", muonSimHitCount(event.getHandle(pltSimHits_)), "muon SimHits in PLT");
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
    description.add<edm::InputTag>("trackerSeeds", edm::InputTag("combinatorialcosmicseedfinderP5"));
    description.add<edm::InputTag>("trackerTrackCandidates", edm::InputTag("ckfTrackCandidatesP5"));
    description.add<edm::InputTag>("trackerRawTracks", edm::InputTag("ctfWithMaterialTracksCosmics"));
    description.add<edm::InputTag>("trackerSelectedTracks", edm::InputTag("ctfWithMaterialTracksP5"));
    description.add<edm::InputTag>("trackerLHCTrackCandidates", edm::InputTag("ckfTrackCandidatesP5LHCNavigation"));
    description.add<edm::InputTag>("trackerLHCTracks", edm::InputTag("ctfWithMaterialTracksP5LHCNavigation"));
    description.add<std::vector<edm::InputTag>>(
        "trackerSimHits",
        {edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
         edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof"),
         edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"),
         edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof"),
         edm::InputTag("g4SimHits", "TrackerHitsTIBLowTof"),
         edm::InputTag("g4SimHits", "TrackerHitsTIBHighTof"),
         edm::InputTag("g4SimHits", "TrackerHitsTIDLowTof"),
         edm::InputTag("g4SimHits", "TrackerHitsTIDHighTof"),
         edm::InputTag("g4SimHits", "TrackerHitsTOBLowTof"),
         edm::InputTag("g4SimHits", "TrackerHitsTOBHighTof"),
         edm::InputTag("g4SimHits", "TrackerHitsTECLowTof"),
         edm::InputTag("g4SimHits", "TrackerHitsTECHighTof")});
    description.add<edm::InputTag>("dsaGlobalLinks", edm::InputTag("shiftGlobalDSAMuons"));
    description.add<edm::InputTag>("cosmicGlobalLinks", edm::InputTag("shiftGlobalCosmicMuons"));
    description.add<edm::InputTag>("traversingGlobalLinks", edm::InputTag("shiftGlobalTraversingMuons"));
    description.add<edm::InputTag>("pixelRecHits", edm::InputTag("siPixelRecHits"));
    description.add<edm::InputTag>("stripMatchedRecHits", edm::InputTag("siStripMatchedRecHits", "matchedRecHit"));
    description.add<edm::InputTag>("stripRphiRecHits", edm::InputTag("siStripMatchedRecHits", "rphiRecHit"));
    description.add<edm::InputTag>("stripRphiUnmatchedRecHits",
                                   edm::InputTag("siStripMatchedRecHits", "rphiRecHitUnmatched"));
    description.add<edm::InputTag>("stripStereoRecHits", edm::InputTag("siStripMatchedRecHits", "stereoRecHit"));
    description.add<edm::InputTag>("stripStereoUnmatchedRecHits",
                                   edm::InputTag("siStripMatchedRecHits", "stereoRecHitUnmatched"));
    description.add<edm::InputTag>("hbheDigis", edm::InputTag("simHcalUnsuppressedDigis"));
    description.add<edm::InputTag>("hfDigis", edm::InputTag("simHcalUnsuppressedDigis"));
    description.add<edm::InputTag>("hoDigis", edm::InputTag("simHcalUnsuppressedDigis"));
    description.add<edm::InputTag>("zdcDigis", edm::InputTag("simHcalUnsuppressedDigis"));
    description.add<edm::InputTag>("hbheQIE11Digis",
                                   edm::InputTag("simHcalUnsuppressedDigis", "HBHEQIE11DigiCollection"));
    description.add<edm::InputTag>("hfQIE10Digis",
                                   edm::InputTag("simHcalUnsuppressedDigis", "HFQIE10DigiCollection"));
    description.add<edm::InputTag>("ecalBarrelRecHits", edm::InputTag("reducedEcalRecHitsEB"));
    description.add<edm::InputTag>("ecalEndcapRecHits", edm::InputTag("reducedEcalRecHitsEE"));
    description.add<edm::InputTag>("hbheRecHits", edm::InputTag("reducedHcalRecHits", "hbhereco"));
    description.add<edm::InputTag>("hfRecHits", edm::InputTag("reducedHcalRecHits", "hfreco"));
    description.add<edm::InputTag>("hoRecHits", edm::InputTag("reducedHcalRecHits", "horeco"));
    description.add<edm::InputTag>("zdcRecHits", edm::InputTag("zdcreco"));
    description.add<edm::InputTag>("bcm1fSimHits", edm::InputTag("g4SimHits", "BCM1FHits"));
    description.add<edm::InputTag>("bhmSimHits", edm::InputTag("g4SimHits", "BHMHits"));
    description.add<edm::InputTag>("pltSimHits", edm::InputTag("g4SimHits", "PLTHits"));
    description.add<bool>("enableDTMeasurement", true);
    description.add<bool>("enableGEMMeasurement", true);
    description.add<int>("trackerMode", 1);
    description.add<bool>("enableHcalDiagnostics", false);
    description.add<bool>("enableZDCDiagnostics", false);
    description.add<int>("dtNavigationMode", 1);
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
  edm::EDGetTokenT<TrajectorySeedCollection> trackerSeeds_;
  edm::EDGetTokenT<TrackCandidateCollection> trackerTrackCandidates_;
  edm::EDGetTokenT<reco::TrackCollection> trackerRawTracks_;
  edm::EDGetTokenT<reco::TrackCollection> trackerSelectedTracks_;
  edm::EDGetTokenT<TrackCandidateCollection> trackerLHCTrackCandidates_;
  edm::EDGetTokenT<reco::TrackCollection> trackerLHCTracks_;
  std::vector<edm::EDGetTokenT<edm::PSimHitContainer>> trackerSimHits_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> dsaGlobalLinks_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> cosmicGlobalLinks_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> traversingGlobalLinks_;
  edm::EDGetTokenT<SiPixelRecHitCollection> pixelRecHits_;
  edm::EDGetTokenT<SiStripMatchedRecHit2DCollection> stripMatchedRecHits_;
  edm::EDGetTokenT<SiStripRecHit2DCollection> stripRphiRecHits_;
  edm::EDGetTokenT<SiStripRecHit2DCollection> stripRphiUnmatchedRecHits_;
  edm::EDGetTokenT<SiStripRecHit2DCollection> stripStereoRecHits_;
  edm::EDGetTokenT<SiStripRecHit2DCollection> stripStereoUnmatchedRecHits_;
  edm::EDGetTokenT<HBHEDigiCollection> hbheDigis_;
  edm::EDGetTokenT<HFDigiCollection> hfDigis_;
  edm::EDGetTokenT<HODigiCollection> hoDigis_;
  edm::EDGetTokenT<ZDCDigiCollection> zdcDigis_;
  edm::EDGetTokenT<QIE11DigiCollection> hbheQIE11Digis_;
  edm::EDGetTokenT<QIE10DigiCollection> hfQIE10Digis_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalBarrelRecHits_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalEndcapRecHits_;
  edm::EDGetTokenT<HBHERecHitCollection> hbheRecHits_;
  edm::EDGetTokenT<HFRecHitCollection> hfRecHits_;
  edm::EDGetTokenT<HORecHitCollection> hoRecHits_;
  edm::EDGetTokenT<ZDCRecHitCollection> zdcRecHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> bcm1fSimHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> bhmSimHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> pltSimHits_;
  bool enableDTMeasurement_;
  bool enableGEMMeasurement_;
  int trackerMode_;
  bool enableHcalDiagnostics_;
  bool enableZDCDiagnostics_;
  int dtNavigationMode_;
};

DEFINE_FWK_MODULE(ShiftMuonRecoDiagnosticsProducer);
