#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include <vector>

class ShiftMuonSegmentsTableProducer : public edm::stream::EDProducer<> {
public:
  explicit ShiftMuonSegmentsTableProducer(edm::ParameterSet const& p)
      : dt_(consumes<DTRecSegment4DCollection>(p.getParameter<edm::InputTag>("dtSegments"))),
        csc_(consumes<CSCSegmentCollection>(p.getParameter<edm::InputTag>("cscSegments"))) {
    produces<nanoaod::FlatTable>("ShiftDT");
    produces<nanoaod::FlatTable>("ShiftCSC");
  }
  void produce(edm::Event& e, edm::EventSetup const&) override {
    std::vector<unsigned int> dr;
    std::vector<int> dw, ds, dse, dnh;
    std::vector<float> dx, dy, dvx, dvy, dc;
    std::vector<unsigned int> cr;
    std::vector<int> ce, cs, crg, cch, cnh;
    std::vector<float> cx, cy, cvx, cvy, cc;
    auto dh = e.getHandle(dt_);
    auto ch = e.getHandle(csc_);
    if (dh.isValid()) for (auto const& s : *dh) { auto id = s.chamberId(); auto pos = s.localPosition(); auto dir = s.localDirection(); dr.push_back(id.rawId()); dw.push_back(id.wheel()); ds.push_back(id.station()); dse.push_back(id.sector()); dx.push_back(pos.x()); dy.push_back(pos.y()); dvx.push_back(dir.x()); dvy.push_back(dir.y()); dc.push_back(s.chi2()); dnh.push_back(s.recHits().size()); }
    if (ch.isValid()) for (auto const& s : *ch) { auto id = s.cscDetId(); auto pos = s.localPosition(); auto dir = s.localDirection(); cr.push_back(id.rawId()); ce.push_back(id.endcap()); cs.push_back(id.station()); crg.push_back(id.ring()); cch.push_back(id.chamber()); cx.push_back(pos.x()); cy.push_back(pos.y()); cvx.push_back(dir.x()); cvy.push_back(dir.y()); cc.push_back(s.chi2()); cnh.push_back(s.nRecHits()); }
    // These are standalone tables, not extensions of tables produced by
    // another module.  Marking them as extensions leaves NanoAOD without the
    // owning nShiftDT/nShiftCSC counter branch and makes the output module
    // dereference a null counter-branch address while filling the event.
    auto dt = std::make_unique<nanoaod::FlatTable>(dr.size(), "ShiftDT", false, false);
    auto ct = std::make_unique<nanoaod::FlatTable>(cr.size(), "ShiftCSC", false, false);
    dt->addColumn<unsigned int>("rawId", dr, "detector raw id"); dt->addColumn<int>("wheel", dw, "DT wheel"); dt->addColumn<int>("station", ds, "DT station"); dt->addColumn<int>("sector", dse, "DT sector"); dt->addColumn<float>("x", dx, "local x"); dt->addColumn<float>("y", dy, "local y"); dt->addColumn<float>("directionX", dvx, "local direction x"); dt->addColumn<float>("directionY", dvy, "local direction y"); dt->addColumn<float>("chi2", dc, "segment chi2"); dt->addColumn<int>("nHits", dnh, "number of hits");
    ct->addColumn<unsigned int>("rawId", cr, "detector raw id"); ct->addColumn<int>("endcap", ce, "CSC endcap"); ct->addColumn<int>("station", cs, "CSC station"); ct->addColumn<int>("ring", crg, "CSC ring"); ct->addColumn<int>("chamber", cch, "CSC chamber"); ct->addColumn<float>("x", cx, "local x"); ct->addColumn<float>("y", cy, "local y"); ct->addColumn<float>("directionX", cvx, "local direction x"); ct->addColumn<float>("directionY", dvy, "local direction y"); ct->addColumn<float>("chi2", cc, "segment chi2"); ct->addColumn<int>("nHits", cnh, "number of hits");
    e.put(std::move(dt), "ShiftDT"); e.put(std::move(ct), "ShiftCSC");
  }
private:
  edm::EDGetTokenT<DTRecSegment4DCollection> dt_;
  edm::EDGetTokenT<CSCSegmentCollection> csc_;
};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonSegmentsTableProducer);
