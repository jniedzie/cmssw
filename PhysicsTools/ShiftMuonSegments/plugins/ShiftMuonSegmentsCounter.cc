#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

class ShiftMuonSegmentsCounter : public edm::one::EDAnalyzer<> {
public:
  explicit ShiftMuonSegmentsCounter(edm::ParameterSet const& p)
      : dt_(consumes<DTRecSegment4DCollection>(p.getParameter<edm::InputTag>("dtSegments"))),
        csc_(consumes<CSCSegmentCollection>(p.getParameter<edm::InputTag>("cscSegments"))) {}
  void analyze(edm::Event const& e, edm::EventSetup const&) override {
    auto dh = e.getHandle(dt_);
    auto ch = e.getHandle(csc_);
    edm::LogInfo("ShiftMuonSegments") << "run " << e.id().run() << " lumi " << e.luminosityBlock() << " event "
                                      << e.id().event()
                                      << ": DT=" << (dh.isValid() ? std::to_string(dh->size()) : "MISSING")
                                      << " CSC=" << (ch.isValid() ? std::to_string(ch->size()) : "MISSING");
  }

private:
  edm::EDGetTokenT<DTRecSegment4DCollection> dt_;
  edm::EDGetTokenT<CSCSegmentCollection> csc_;
};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonSegmentsCounter);
