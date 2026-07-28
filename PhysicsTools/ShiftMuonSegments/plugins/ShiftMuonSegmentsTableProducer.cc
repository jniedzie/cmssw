#include "FWCore/Framework/interface/Frameworkfwd.h"
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
    if (dh.isValid())
      for (auto const& p : *dh) {
        auto const& s = p.second;
        auto id = p.first.rawId();
        auto pos = s.localPosition();
        auto dir = s.localDirection();
        dr.push_back(id);
        dw.push_back(p.first.wheel());
        ds.push_back(p.first.station());
        dse.push_back(p.first.sector());
        dx.push_back(pos.x());
        dy.push_back(pos.y());
        dvx.push_back(dir.x());
        dvy.push_back(dir.y());
        dc.push_back(s.chi2());
        dnh.push_back(s.specificRecHits().size());
      }
    if (ch.isValid())
      for (auto const& p : *ch) {
        auto const& s = p.second;
        auto id = p.first.rawId();
        auto pos = s.localPosition();
        auto dir = s.localDirection();
        cr.push_back(id);
        ce.push_back(p.first.endcap());
        cs.push_back(p.first.station());
        crg.push_back(p.first.ring());
        cch.push_back(p.first.chamber());
        cx.push_back(pos.x());
        cy.push_back(pos.y());
        cvx.push_back(dir.x());
        cvy.push_back(dir.y());
        cc.push_back(s.chi2());
        cnh.push_back(s.specificRecHits().size());
      }
    auto dt = std::make_unique<nanoaod::FlatTable>(dr.size(), "ShiftDT", false, true);
    auto ct = std::make_unique<nanoaod::FlatTable>(cr.size(), "ShiftCSC", false, true);
    dt->addColumn("rawId", dr, "detector raw id");
    dt->addColumn("wheel", dw, "DT wheel");
    dt->addColumn("station", ds, "DT station");
    dt->addColumn("sector", dse, "DT sector");
    dt->addColumn("x", dx, "local x");
    dt->addColumn("y", dy, "local y");
    dt->addColumn("directionX", dvx, "local direction x");
    dt->addColumn("directionY", dvy, "local direction y");
    dt->addColumn("chi2", dc, "segment chi2");
    dt->addColumn("nHits", dnh, "number of hits");
    ct->addColumn("rawId", cr, "detector raw id");
    ct->addColumn("endcap", ce, "CSC endcap");
    ct->addColumn("station", cs, "CSC station");
    ct->addColumn("ring", crg, "CSC ring");
    ct->addColumn("chamber", cch, "CSC chamber");
    ct->addColumn("x", cx, "local x");
    ct->addColumn("y", cy, "local y");
    ct->addColumn("directionX", cvx, "local direction x");
    ct->addColumn("directionY", cvy, "local direction y");
    ct->addColumn("chi2", cc, "segment chi2");
    ct->addColumn("nHits", cnh, "number of hits");
    e.put(std::move(dt), "ShiftDT");
    e.put(std::move(ct), "ShiftCSC");
  }

private:
  edm::EDGetTokenT<DTRecSegment4DCollection> dt_;
  edm::EDGetTokenT<CSCSegmentCollection> csc_;
};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonSegmentsTableProducer);
