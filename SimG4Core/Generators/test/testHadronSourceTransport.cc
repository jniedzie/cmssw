// Code style: AI-code. AI-generated native finite-hadron and source-lineage controls.
// Primary-tree tests do not establish Geant4 material or detector acceptance.
#include "SimG4Core/Generators/interface/Generator.h"
#include "SimG4Core/Generators/interface/Generator3.h"
#include "SimG4Core/Generators/interface/SourceAncestry.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HepMC/GenVertex.h"
#include "HepMC3/GenVertex.h"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
  void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
  }
  void close(double value, double expected) {
    require(std::isfinite(value) && std::abs(value-expected)<1.e-10*std::max(1.,std::abs(expected)),
            "Wrong native preassigned lifetime or momentum");
  }
  edm::ParameterSet parameters(bool enabled, bool explicitOption=true) {
    edm::ParameterSet p;
    for (const char* name : {"ApplyPCuts", "ApplyEtaCuts"}) p.addParameter<bool>(name,true);
    for (const char* name : {"ApplyPtransCut", "ApplyPhiCuts", "ApplyLumiMonitorCuts", "IsSlepton"})
      p.addParameter<bool>(name,false);
    p.addParameter<double>("MinPhiCut",-3.141592653589793); p.addParameter<double>("MaxPhiCut",3.141592653589793);
    p.addParameter<double>("MinEtaCut",-5.5); p.addParameter<double>("MaxEtaCut",5.5);
    p.addParameter<double>("MinPCut",.04); p.addParameter<double>("MaxPCut",99999.);
    p.addParameter<double>("EtaCutForHector",8.2); p.addParameter<double>("RDecLenCut",2.9);
    p.addParameter<double>("LDecLenCut",10.); p.addParameter<double>("MaxZCentralCMS",50.);
    if (explicitOption) p.addUntrackedParameter<bool>("TrackAllHadronDecays",enabled);
    return p;
  }
  struct Particle { int pdg,status,birth,end; G4LorentzVector p4; };
  struct Graph { std::vector<G4LorentzVector> vertices; std::vector<Particle> particles; };
  Graph cascade(double z, double beta, bool prompt=false) {
    const auto mass=[](int pdg){ return G4ParticleTable::GetParticleTable()->FindParticle(pdg)->GetPDGMass()/GeV; };
    const double mb=mass(511),md=mass(411),mm=mass(13);
    const double first=(std::pow(mb-md,2)-mm*mm)/(2*(mb-md));
    const double second=(md*md-mm*mm)/(2*md);
    Graph g;
    g.particles={{511,2,0,1,{0.,0.,0.,mb}},{411,2,1,2,{0.,0.,0.,md}},
                 {13,1,1,-1,{first,0.,0.,std::hypot(first,mm)}},{-14,1,1,-1,{-first,0.,0.,first}},
                 {-13,1,2,-1,{0.,second,0.,std::hypot(second,mm)}},{14,1,2,-1,{0.,-second,0.,second}}};
    for (auto& p:g.particles) p.p4.boost(0.,0.,beta);
    const G4LorentzVector origin(0.,0.,z,0.);
    const auto firstDecay=prompt?origin:origin+g.particles[0].p4*(30./mb);
    const auto secondDecay=prompt?origin:firstDecay+g.particles[1].p4*(40./md);
    g.vertices={origin,firstDecay,secondDecay};
    return g;
  }
  std::unique_ptr<HepMC::GenEvent> event2(const Graph& g, bool reverse) {
    auto e=std::make_unique<HepMC::GenEvent>(HepMC::Units::GEV,HepMC::Units::MM);
    std::vector<HepMC::GenVertex*> v;
    for (const auto& x:g.vertices) v.push_back(new HepMC::GenVertex({x.x(),x.y(),x.z(),x.t()}));
    for (const auto& p:g.particles) {
      auto* q=new HepMC::GenParticle({p.p4.x(),p.p4.y(),p.p4.z(),p.p4.t()},p.pdg,p.status);
      v[p.birth]->add_particle_out(q); if(p.end>=0)v[p.end]->add_particle_in(q);
    }
    if(reverse)std::reverse(v.begin(),v.end());
    for(auto* x:v)e->add_vertex(x);
    e->weights().push_back(-.03125); return e;
  }
  std::unique_ptr<HepMC3::GenEvent> event3(const Graph& g, bool reverse) {
    auto e=std::make_unique<HepMC3::GenEvent>(HepMC3::Units::GEV,HepMC3::Units::MM);
    std::vector<HepMC3::GenVertexPtr> v;
    for (const auto& x:g.vertices) v.push_back(std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(x.x(),x.y(),x.z(),x.t())));
    for (const auto& p:g.particles) {
      auto q=std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(p.p4.x(),p.p4.y(),p.p4.z(),p.p4.t()),p.pdg,p.status);
      v[p.birth]->add_particle_out(q); if(p.end>=0)v[p.end]->add_particle_in(q);
    }
    if(reverse)std::reverse(v.begin(),v.end());
    for(const auto& x:v)e->add_vertex(x);
    e->weights().push_back(-.03125); return e;
  }
  std::unique_ptr<G4Event> convert(const Graph& g, bool hep3, bool enabled, bool reverse=false, bool explicitOption=true, double tinyFlightMm=0.) {
    auto config=parameters(enabled,explicitOption);
    if(tinyFlightMm!=0.) config.addUntrackedParameter<double>("MaxUntrackedHadronFlightMm",tinyFlightMm);
    auto out=std::make_unique<G4Event>(1);
    if(hep3) {
      auto e=event3(g,reverse); Generator3 gen(config); gen.HepMC2G4(e.get(),out.get());
      close(gen.eventWeight(),-.03125); require(e->weights()[0]==-.03125,"HepMC3 native weight changed");
      for(const auto& p:e->particles())require(p->status()<1000,"HepMC3 source graph mutated");
    } else {
      auto e=event2(g,reverse); Generator gen(config); gen.HepMC2G4(e.get(),out.get());
      close(gen.eventWeight(),-.03125); require(e->weights()[0]==-.03125,"HepMC2 native weight changed");
      for(auto p=e->particles_begin();p!=e->particles_end();++p)require((*p)->status()<1000,"HepMC2 source graph mutated");
    }
    return out;
  }
  std::vector<const G4PrimaryParticle*> roots(const G4Event& e) {
    std::vector<const G4PrimaryParticle*> result;
    for(int i=0;i<e.GetNumberOfPrimaryVertex();++i)
      for(auto* p=e.GetPrimaryVertex(i)->GetPrimary();p;p=p->GetNext())result.push_back(p);
    return result;
  }
  std::vector<double> signature(const G4Event& e) {
    std::vector<double> result;
    const auto append=[&](const auto& self,const G4PrimaryParticle* p)->void {
      for(;p;p=p->GetNext()) {
        const auto momentum=p->GetMomentum();
        result.insert(result.end(),{double(p->GetPDGcode()),momentum.x(),momentum.y(),momentum.z(),p->GetTotalEnergy(),p->GetProperTime()});
        self(self,p->GetDaughter());result.push_back(-999999.);
      }
    };
    for(int i=0;i<e.GetNumberOfPrimaryVertex();++i) {
      const auto* v=e.GetPrimaryVertex(i);
      result.insert(result.end(),{v->GetX0(),v->GetY0(),v->GetZ0(),v->GetT0()});append(append,v->GetPrimary());
    }
    return result;
  }
  void checkMomenta(const G4Event& e,const Graph& g) {
    const auto check=[&](const auto& self,const G4PrimaryParticle* p)->void {
      for(;p;p=p->GetNext()) {
        const auto expected=std::find_if(g.particles.begin(),g.particles.end(),[&](const auto& q){return q.pdg==p->GetPDGcode();});
        require(expected!=g.particles.end(),"Unexpected primary-tree particle");
        const auto momentum=p->GetMomentum();
        close(momentum.x()/GeV,expected->p4.x());close(momentum.y()/GeV,expected->p4.y());
        close(momentum.z()/GeV,expected->p4.z());close(p->GetTotalEnergy()/GeV,expected->p4.t());
        self(self,p->GetDaughter());
      }
    };
    for(const auto* p:roots(e))check(check,p);
  }
  void checkCascade(const G4Event& e) {
    auto primaries=roots(e); require(primaries.size()==1 && primaries[0]->GetPDGcode()==511,"Missing/duplicate finite B root");
    close(primaries[0]->GetProperTime(),30.*mm/CLHEP::c_light);
    const G4PrimaryParticle* d=nullptr; int children=0;
    for(auto* q=primaries[0]->GetDaughter();q;q=q->GetNext()){++children;if(q->GetPDGcode()==411)d=q;}
    require(children==3 && d,"Missing B daughters or nested D"); close(d->GetProperTime(),40.*mm/CLHEP::c_light);
    children=0; for(auto* q=d->GetDaughter();q;q=q->GetNext())++children;
    require(children==2,"Missing/duplicate D decay daughters");
  }
  void lineage() {
    const auto track=[](int id,int parent,double z) {
      auto t=std::make_unique<G4Track>(new G4DynamicParticle(G4ParticleTable::GetParticleTable()->FindParticle(13),{0.,0.,1.},40.*GeV),0.,G4ThreeVector(0.,0.,z));
      t->SetTrackID(id);t->SetParentID(parent);return t;
    };
    for(double z:{-148000.,148000.}) {
      sim::sourceAncestry::reset(true,50000.);
      auto root=track(1,0,z);sim::sourceAncestry::remember(root.get());require(sim::sourceAncestry::isExternal(root.get()),"External root not recognized");
      root->SetPosition({0.,0.,0.});sim::sourceAncestry::remember(root.get());
      auto child=track(2,1,8000.);sim::sourceAncestry::remember(child.get());require(sim::sourceAncestry::isExternal(child.get()),"Central child lost its external source");
    }
    sim::sourceAncestry::reset(true,50000.);auto ip=track(1,0,0.);sim::sourceAncestry::remember(ip.get());ip->SetPosition({0.,0.,148000.});sim::sourceAncestry::remember(ip.get());require(!sim::sourceAncestry::isExternal(ip.get()),"IP root changed origin during propagation");
    sim::sourceAncestry::reset(true,50000.);bool missing=false;try{sim::sourceAncestry::isExternal(ip.get());}catch(...){missing=true;}require(missing,"Old event lineage survived reset");
    sim::sourceAncestry::reset(false,50000.);require(!sim::sourceAncestry::isExternal(ip.get()),"Disabled mode retained external source behavior");
  }
}
int main() {
  try {
    G4BosonConstructor().ConstructParticle();G4LeptonConstructor().ConstructParticle();G4MesonConstructor().ConstructParticle();G4ParticleTable::GetParticleTable()->SetReadiness();
    int checked=0;
    for(bool hep3:{false,true}) {
      for(double z:{-148000.,0.,148000.})for(double beta:{-.8,.8})for(bool reverse:{false,true}) {
        auto g=cascade(z,beta);auto converted=convert(g,hep3,true,reverse);checkCascade(*converted);checkMomenta(*converted,g);++checked;
      }
      checkCascade(*convert(cascade(148000.,0.),hep3,true));++checked;
      {
        auto tiny=cascade(148000.,-.8);tiny.particles[1].pdg=413;
        tiny.vertices[2]=tiny.vertices[1]+G4LorentzVector(0.,0.,-1.2e-6,1.2e-6);
        bool failed=false;try{convert(tiny,hep3,true);}catch(...){failed=true;}
        require(failed,"Default accepted finite unknown resonance");
        auto projected=convert(tiny,hep3,true,false,true,1.e-5);
        require(roots(*projected).size()==1,"Bounded unknown became an independent primary");
        checkMomenta(*projected,tiny);++checked;
        tiny.vertices[2]=tiny.vertices[1]+G4LorentzVector(0.,0.,-1.e-3,1.e-3);
        failed=false;try{convert(tiny,hep3,true,false,true,1.e-5);}catch(...){failed=true;}
        require(failed,"Over-budget finite unknown resonance accepted");++checked;
      }
      auto prompt=cascade(148000.,-.8,true);auto omitted=convert(prompt,hep3,false,false,false);auto disabled=convert(prompt,hep3,false);auto enabled=convert(prompt,hep3,true);
      require(signature(*omitted)==signature(*disabled),"Default zero4 prompt history changed");
      // Opt-in conversion now preserves GEN p4/mass explicitly. Compare the
      // enabled physical values to the source, retaining exact default parity.
      require(roots(*enabled).size()==roots(*disabled).size(),"Opt-in prompt primary count changed");
      checkMomenta(*enabled,prompt);
      for(int pdg:{13,211})for(double direction:{-1.,1.}) {
        Graph stable{{{0.,0.,148000.,0.}},{{pdg,1,0,-1,{0.,0.,10.*direction,std::hypot(10.,G4ParticleTable::GetParticleTable()->FindParticle(pdg)->GetPDGMass()/GeV)}}}};
        auto result=convert(stable,hep3,true);require(roots(*result).size()==1 && roots(*result)[0]->GetPDGcode()==pdg,"External stable muon/hadron lost");++checked;
      }
      for(int pdg:{211,321}) {
        const double parentMass=G4ParticleTable::GetParticleTable()->FindParticle(pdg)->GetPDGMass()/GeV;
        const double muMass=G4ParticleTable::GetParticleTable()->FindParticle(13)->GetPDGMass()/GeV;
        const double q=(parentMass*parentMass-muMass*muMass)/(2.*parentMass);
        Graph rest{{{0.,0.,148000.,0.},{0.,0.,148000.,25.}},
          {{pdg,2,0,1,{0.,0.,0.,parentMass}},{-13,1,1,-1,{q,0.,0.,std::hypot(q,muMass)}},{14,1,1,-1,{-q,0.,0.,q}}}};
        auto result=convert(rest,hep3,true);
        require(roots(*result).size()==1 && roots(*result)[0]->GetPDGcode()==pdg,"At-rest pion/kaon root lost");
        close(roots(*result)[0]->GetProperTime(),25.*mm/CLHEP::c_light);checkMomenta(*result,rest);++checked;
      }
      for(int bad:{0,1,2,3,4}) {
        auto g=cascade(148000.,0.);
        if(bad==0)g.vertices[1].setT(-1.);
        if(bad==1)g.particles[0].pdg=100113;
        if(bad==2)g.vertices[1].setT(std::numeric_limits<double>::infinity());
        if(bad==3)g.particles[2].pdg=100113;
        if(bad==4)g.particles[1].status=3;
        bool rejected=false;try{convert(g,hep3,true);}catch(...){rejected=true;}require(rejected,"Invalid clock, undefined particle or unsupported subtree accepted");
      }
    }
    lineage();std::cout<<"Finite hadron source: "<<checked<<" HepMC2/3 tree controls and native source-lineage controls passed; no material claim.\n";return 0;
  } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
