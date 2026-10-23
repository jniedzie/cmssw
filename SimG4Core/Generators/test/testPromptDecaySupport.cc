// Code style: AI-code. AI-generated exact-zero prompt-support native controls.
#include "SimG4Core/Generators/interface/PromptDecaySupport.h"
#include "HepMC/GenEvent.h"
#include "HepMC3/GenEvent.h"
#include <iostream>
#include <limits>
#include <memory>

void require(bool x, const char* why) { if (!x) throw std::runtime_error(why); }
template<typename Call> void rejects(Call call, const char* reason) {
  try { call(); } catch (const std::runtime_error& e) {
    require(std::string(e.what()).find(reason) != std::string::npos,"Wrong rejection reason"); return;
  }
  throw std::runtime_error("An unsupported prompt tree was accepted");
}
int main() {
  const auto known=[](int id){return id!=413 && id!=-423 && id!=999999;};
  int checks=0;
  {
    HepMC::GenEvent event(HepMC::Units::GEV,HepMC::Units::MM);
    event.weights().push_back(-.03125);
    auto* birth=new HepMC::GenVertex({0.,0.,148000.,-148000.});
    auto* end=new HepMC::GenVertex({0.,0.,148000.,-148000.});
    auto* resonance=new HepMC::GenParticle({0.,0.,10.,11.},-423,2);
    auto* daughter=new HepMC::GenParticle({0.,0.,8.,9.},-421,2);
    auto* photon=new HepMC::GenParticle({0.,0.,2.,2.},22,1);
    birth->add_particle_out(resonance);end->add_particle_in(resonance);
    end->add_particle_out(daughter);end->add_particle_out(photon);event.add_vertex(birth);event.add_vertex(end);
    const auto use=[&]{return sim::promptDecay::frontier(*resonance,known);};
    require(use()==std::vector<const HepMC::GenParticle*>({daughter,photon}),"Wrong HepMC2 frontier");++checks;
    require(resonance->status()==2 && daughter->status()==2 && event.weights()[0]==-.03125,"Source mutated");++checks;
    require(use()[0]->momentum()==daughter->momentum(),"Source p4 changed");++checks;
    end->set_position({0.,0.,std::nextafter(148000.,148001.),-148000.});rejects(use,"nonzero");++checks;
    end->set_position({0.,0.,148000.,std::nextafter(-148000.,-147999.)});rejects(use,"nonzero");++checks;
    end->set_position({0.,0.,148000.,-148000.});resonance->set_momentum({0.,0.,0.,2.});
    require(sim::promptDecay::exactlyPrompt(*resonance),"At-rest prompt misclassified");++checks;
    resonance->set_status(1);rejects(use,"terminal");++checks;resonance->set_status(2);
    auto* other=new HepMC::GenParticle({0.,0.,0.,1.},22,1);end->add_particle_in(other);rejects(use,"incoming");++checks;end->remove_particle(other);delete other;
    end->set_position({std::numeric_limits<double>::quiet_NaN(),0.,148000.,-148000.});rejects(use,"nonfinite");++checks;
  }
  {
    auto a=std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0.,0.,10.,11.),413,2);
    auto b=std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0.,0.,8.,9.),-423,2);
    auto d=std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0.,0.,8.,9.),421,2);
    auto v0=std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(0.,0.,148000.,-148000.));
    auto v1=std::make_shared<HepMC3::GenVertex>(v0->position());auto v2=std::make_shared<HepMC3::GenVertex>(v0->position());
    v0->add_particle_out(a);v1->add_particle_in(a);v1->add_particle_out(b);v2->add_particle_in(b);v2->add_particle_out(d);
    HepMC3::GenEvent event(HepMC3::Units::GEV,HepMC3::Units::MM);event.add_vertex(v0);event.add_vertex(v1);event.add_vertex(v2);event.weights().push_back(-.03125);
    const auto use=[&]{return sim::promptDecay::frontier(*a,known);};
    require(use()==std::vector<const HepMC3::GenParticle*>({d.get()}),"Wrong nested HepMC3 frontier");++checks;
    require(event.weights()[0]==-.03125 && a->status()==2 && b->status()==2,"HepMC3 source mutated");++checks;
    v2->set_position({0.,0.,148000.,-147999.});rejects(use,"nonzero");++checks;v2->set_position(v1->position());
    d->set_pid(999999);rejects(use,"terminal");++checks;d->set_pid(421);
  }
  {
    // No owning GenEvent: remove the cycle before deleting these two vertices.
    auto* a=new HepMC::GenParticle({0.,0.,0.,2.},413,2);auto* b=new HepMC::GenParticle({0.,0.,0.,2.},-423,2);
    auto* v1=new HepMC::GenVertex({0.,0.,0.,0.});auto* v2=new HepMC::GenVertex({0.,0.,0.,0.});
    v1->add_particle_in(a);v1->add_particle_out(b);v2->add_particle_in(b);v2->add_particle_out(a);
    rejects([&]{return sim::promptDecay::frontier(*a,known);},"cycle");++checks;
    v1->remove_particle(b);v2->remove_particle(b);v1->remove_particle(a);v2->remove_particle(a);delete a;delete b;delete v1;delete v2;
  }
  {
    HepMC::GenEvent event(HepMC::Units::GEV,HepMC::Units::MM);
    auto* a=new HepMC::GenVertex({0.,0.,0.,0.});auto* z=new HepMC::GenVertex({0.,0.,1.2e-6,1.2e-6});
    auto* p=new HepMC::GenParticle({0.,0.,10.,11.},413,2);
    auto* d=new HepMC::GenParticle({0.,0.,10.,11.},421,1);
    a->add_particle_out(p);z->add_particle_in(p);z->add_particle_out(d);event.add_vertex(a);event.add_vertex(z);
    const auto bounded=[&]{return sim::promptDecay::frontier(*p,known,1.e-5);};
    rejects([&]{sim::promptDecay::frontier(*p,known);},"nonzero");++checks;
    require(bounded()==std::vector<const HepMC::GenParticle*>({d}),"Tiny H2 flight rejected");++checks;
    require(z->position().z()==1.2e-6 && p->status()==2,"Tiny H2 source mutated");++checks;
    z->set_position({0.,0.,1.e-5,1.e-5});require(bounded().size()==1,"H2 budget edge rejected");++checks;
    z->set_position({0.,0.,std::nextafter(1.e-5,1.),1.e-5});rejects(bounded,"nonzero");++checks;
    z->set_position({0.,0.,0.,std::nextafter(1.e-5,1.)});rejects(bounded,"nonzero");++checks;
    z->set_position({0.,0.,0.,-1.e-6});rejects(bounded,"nonzero");++checks;
    rejects([&]{sim::promptDecay::frontier(*p,known,-1.);},"Invalid");++checks;
  }
  {
    auto a=std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(0.,0.,0.,0.));
    auto z=std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(0.,0.,1.2e-6,1.2e-6));
    auto p=std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0.,0.,10.,11.),413,2);
    auto d=std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0.,0.,10.,11.),421,1);
    a->add_particle_out(p);z->add_particle_in(p);z->add_particle_out(d);
    HepMC3::GenEvent event(HepMC3::Units::GEV,HepMC3::Units::MM);event.add_vertex(a);event.add_vertex(z);
    const auto bounded=[&]{return sim::promptDecay::frontier(*p,known,1.e-5);};
    rejects([&]{sim::promptDecay::frontier(*p,known);},"nonzero");++checks;
    require(bounded()==std::vector<const HepMC3::GenParticle*>({d.get()}),"Tiny H3 flight rejected");++checks;
    require(z->position().z()==1.2e-6 && p->status()==2,"Tiny H3 source mutated");++checks;
    z->set_position({0.,0.,1.e-5,1.e-5});require(bounded().size()==1,"H3 budget edge rejected");++checks;
    z->set_position({0.,0.,std::nextafter(1.e-5,1.),1.e-5});rejects(bounded,"nonzero");++checks;
    z->set_position({0.,0.,0.,std::nextafter(1.e-5,1.)});rejects(bounded,"nonzero");++checks;
    z->set_position({0.,0.,0.,-1.e-6});rejects(bounded,"nonzero");++checks;
    rejects([&]{sim::promptDecay::frontier(*p,known,std::numeric_limits<double>::infinity());},"Invalid");++checks;
  }
  std::cout<<"{\"checks\":"<<checks<<",\"passed\":true,\"default_nonzero_flight_never_collapsed\":true}"<<std::endl;
}
