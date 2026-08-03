#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include <vector>

class ShiftMuonSegmentsTableProducer : public edm::stream::EDProducer<> {
public:
  explicit ShiftMuonSegmentsTableProducer(edm::ParameterSet const& p)
      : dt_(consumes<DTRecSegment4DCollection>(p.getParameter<edm::InputTag>("dtSegments"))),
        csc_(consumes<CSCSegmentCollection>(p.getParameter<edm::InputTag>("cscSegments"))),
        rpc_(consumes<RPCRecHitCollection>(p.getParameter<edm::InputTag>("rpcRecHits"))),
        gem_(consumes<GEMSegmentCollection>(p.getParameter<edm::InputTag>("gemSegments"))) {
    produces<nanoaod::FlatTable>("ShiftDT");
    produces<nanoaod::FlatTable>("ShiftCSC");
    produces<nanoaod::FlatTable>("ShiftRPC");
    produces<nanoaod::FlatTable>("ShiftGEM");
  }
  void produce(edm::Event& e, edm::EventSetup const&) override {
    std::vector<unsigned int> dr;
    std::vector<int> dw, ds, dse, dnh;
    std::vector<float> dx, dy, dvx, dvy, dc;
    std::vector<unsigned int> cr;
    std::vector<int> ce, cs, crg, cch, cnh;
    std::vector<float> cx, cy, cvx, cvy, cc;
    std::vector<unsigned int> rr, gr;
    std::vector<int> rrg, rst, rse, rly, rbx, ge, gs, grg, gch, gly, gnh;
    std::vector<float> rx, ry, rtime, gx, gy, gvx, gvy, gc;
    auto dh = e.getHandle(dt_);
    auto ch = e.getHandle(csc_);
    auto rh = e.getHandle(rpc_);
    auto gh = e.getHandle(gem_);
    if (dh.isValid()) for (auto const& s : *dh) { auto id = s.chamberId(); auto pos = s.localPosition(); auto dir = s.localDirection(); dr.push_back(id.rawId()); dw.push_back(id.wheel()); ds.push_back(id.station()); dse.push_back(id.sector()); dx.push_back(pos.x()); dy.push_back(pos.y()); dvx.push_back(dir.x()); dvy.push_back(dir.y()); dc.push_back(s.chi2()); dnh.push_back(s.recHits().size()); }
    if (ch.isValid()) for (auto const& s : *ch) { auto id = s.cscDetId(); auto pos = s.localPosition(); auto dir = s.localDirection(); cr.push_back(id.rawId()); ce.push_back(id.endcap()); cs.push_back(id.station()); crg.push_back(id.ring()); cch.push_back(id.chamber()); cx.push_back(pos.x()); cy.push_back(pos.y()); cvx.push_back(dir.x()); cvy.push_back(dir.y()); cc.push_back(s.chi2()); cnh.push_back(s.nRecHits()); }
    if (rh.isValid()) for (auto const& h : *rh) { auto id = h.rpcId(); auto pos = h.localPosition(); rr.push_back(id.rawId()); rrg.push_back(id.region()); rst.push_back(id.station()); rse.push_back(id.sector()); rly.push_back(id.layer()); rbx.push_back(h.BunchX()); rx.push_back(pos.x()); ry.push_back(pos.y()); rtime.push_back(h.time()); }
    if (gh.isValid()) for (auto const& s : *gh) { auto id = s.gemDetId(); auto pos = s.localPosition(); auto dir = s.localDirection(); gr.push_back(id.rawId()); ge.push_back(id.region()); gs.push_back(id.station()); grg.push_back(id.ring()); gch.push_back(id.chamber()); gly.push_back(id.layer()); gx.push_back(pos.x()); gy.push_back(pos.y()); gvx.push_back(dir.x()); gvy.push_back(dir.y()); gc.push_back(s.chi2()); gnh.push_back(s.nRecHits()); }
    // These are standalone tables, not extensions of tables produced by
    // another module.  Marking them as extensions leaves NanoAOD without the
    // owning nShiftDT/nShiftCSC counter branch and makes the output module
    // dereference a null counter-branch address while filling the event.
    auto dt = std::make_unique<nanoaod::FlatTable>(dr.size(), "ShiftDT", false, false);
    auto ct = std::make_unique<nanoaod::FlatTable>(cr.size(), "ShiftCSC", false, false);
    auto rpcTable = std::make_unique<nanoaod::FlatTable>(rr.size(), "ShiftRPC", false, false);
    auto gt = std::make_unique<nanoaod::FlatTable>(gr.size(), "ShiftGEM", false, false);
    dt->addColumn<unsigned int>("rawId", dr, "detector raw id"); dt->addColumn<int>("wheel", dw, "DT wheel"); dt->addColumn<int>("station", ds, "DT station"); dt->addColumn<int>("sector", dse, "DT sector"); dt->addColumn<float>("x", dx, "local x"); dt->addColumn<float>("y", dy, "local y"); dt->addColumn<float>("directionX", dvx, "local direction x"); dt->addColumn<float>("directionY", dvy, "local direction y"); dt->addColumn<float>("chi2", dc, "segment chi2"); dt->addColumn<int>("nHits", dnh, "number of hits");
    ct->addColumn<unsigned int>("rawId", cr, "detector raw id"); ct->addColumn<int>("endcap", ce, "CSC endcap"); ct->addColumn<int>("station", cs, "CSC station"); ct->addColumn<int>("ring", crg, "CSC ring"); ct->addColumn<int>("chamber", cch, "CSC chamber"); ct->addColumn<float>("x", cx, "local x"); ct->addColumn<float>("y", cy, "local y"); ct->addColumn<float>("directionX", cvx, "local direction x"); ct->addColumn<float>("directionY", cvy, "local direction y"); ct->addColumn<float>("chi2", cc, "segment chi2"); ct->addColumn<int>("nHits", cnh, "number of hits");
    rpcTable->addColumn<unsigned int>("rawId", rr, "detector raw id"); rpcTable->addColumn<int>("region", rrg, "RPC region"); rpcTable->addColumn<int>("station", rst, "RPC station"); rpcTable->addColumn<int>("sector", rse, "RPC sector"); rpcTable->addColumn<int>("layer", rly, "RPC layer"); rpcTable->addColumn<int>("bx", rbx, "reconstructed bunch crossing"); rpcTable->addColumn<float>("x", rx, "local x"); rpcTable->addColumn<float>("y", ry, "local y"); rpcTable->addColumn<float>("time", rtime, "reconstructed time");
    gt->addColumn<unsigned int>("rawId", gr, "detector raw id"); gt->addColumn<int>("region", ge, "GEM region"); gt->addColumn<int>("station", gs, "GEM station"); gt->addColumn<int>("ring", grg, "GEM ring"); gt->addColumn<int>("chamber", gch, "GEM chamber"); gt->addColumn<int>("layer", gly, "GEM layer"); gt->addColumn<float>("x", gx, "local x"); gt->addColumn<float>("y", gy, "local y"); gt->addColumn<float>("directionX", gvx, "local direction x"); gt->addColumn<float>("directionY", gvy, "local direction y"); gt->addColumn<float>("chi2", gc, "segment chi2"); gt->addColumn<int>("nHits", gnh, "number of hits");
    e.put(std::move(dt), "ShiftDT"); e.put(std::move(ct), "ShiftCSC"); e.put(std::move(rpcTable), "ShiftRPC"); e.put(std::move(gt), "ShiftGEM");
  }
private:
  edm::EDGetTokenT<DTRecSegment4DCollection> dt_;
  edm::EDGetTokenT<CSCSegmentCollection> csc_;
  edm::EDGetTokenT<RPCRecHitCollection> rpc_;
  edm::EDGetTokenT<GEMSegmentCollection> gem_;
};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonSegmentsTableProducer);
