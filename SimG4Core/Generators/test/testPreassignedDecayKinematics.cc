// Code style: AI-code. AI-generated native source-to-decay conversion regression.
// Changes to parent energy/direction below are controlled transport states,
// not measurements of material loss.
#include "SimG4Core/Generators/interface/PreassignedDecayKinematics.h"
#include "SimG4Core/Generators/interface/RecordedPrimaryMomentum.h"
#include "SimG4Core/Generators/interface/SourceAncestry.h"
#include "G4Decay.hh"
#include "G4DecayProducts.hh"
#include "G4Event.hh"
#include "G4MesonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryTransformer.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleTable.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
  void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
  }
  G4PrimaryParticle* primary(int pdg, const G4LorentzVector& p, double mass) {
    auto* result = new G4PrimaryParticle(pdg);
    result->Set4Momentum(p.x(), p.y(), p.z(), p.e());
    result->SetMass(mass);
    result->SetTotalEnergy(p.e());
    result->SetProperTime(.123 * ns);
    return result;
  }
  std::pair<G4LorentzVector, G4LorentzVector> twoBody(
      const G4LorentzVector& parent, double mass, double a, double b) {
    const double q = std::sqrt((mass*mass-(a+b)*(a+b))*(mass*mass-(a-b)*(a-b)))/(2*mass);
    G4LorentzVector first(.6*q, 0., .8*q, std::hypot(q,a));
    G4LorentzVector second(-.6*q, 0., -.8*q, std::hypot(q,b));
    first.boost(parent.boostVector()); second.boost(parent.boostVector());
    return {first, second};
  }
  std::vector<G4Track*> decay(G4Track& parent, int state) {
    if (state) {
      const double energy = state == 2 ? parent.GetDynamicParticle()->GetMass() : 80.*GeV;
      parent.SetKineticEnergy(energy-parent.GetDynamicParticle()->GetMass());
      parent.SetMomentumDirection({.3,.4,-std::sqrt(.75)});
    }
    parent.SetPosition({2.*cm,3.*cm,4.*cm});
    parent.SetTouchableHandle(G4TouchableHandle(new G4TouchableHistory()));
    parent.SetTrackStatus(state==2 ? fStopButAlive : fAlive);
    G4Decay process;
    process.StartTracking(&parent);
    G4Step step; step.SetTrack(&parent); parent.SetStep(&step);
    std::cerr << std::setprecision(17) << "decay entry actualMassGeV=" << parent.GetDynamicParticle()->GetMass()/GeV << " productsSnapshotMassGeV=" << parent.GetDynamicParticle()->GetPreAssignedDecayProducts()->GetParentParticle()->GetMass()/GeV << "\n";
    G4VParticleChange* result;
    if (state==2) {
      G4ForceCondition condition;
      process.AtRestGetPhysicalInteractionLength(parent,&condition);
      result=process.AtRestDoIt(parent,step);
    } else result=process.PostStepDoIt(parent,step);
    std::vector<G4Track*> children;
    for(int i=0;i<result->GetNumberOfSecondaries();++i) {
      auto* child=result->GetSecondary(i);
      require(child->GetPosition()==parent.GetPosition(),"Decay birth vertex changed");
      children.push_back(child);
    }
    result->Clear(); parent.SetStep(nullptr);
    static int nextId = 10;
    for (auto* child : children) {
      child->SetParentID(parent.GetTrackID());
      child->SetTrackID(nextId++);
      sim::sourceAncestry::remember(child);
      const auto momentum = child->GetMomentum();
      const auto weight = child->GetWeight();
      const auto properTime = child->GetDynamicParticle()->GetPreAssignedDecayProperTime();
      const auto* products = child->GetDynamicParticle()->GetPreAssignedDecayProducts();
      const auto* source = child->GetDynamicParticle()->GetPrimaryParticle();
      const auto sourceMomentum = source->GetMomentum();
      const auto sourceEnergy = source->GetTotalEnergy();
      sim::preservePreassignedDaughterMass(*child);
      const double roundtrip = (child->GetMomentum()-momentum).mag();
      require(roundtrip <= 32*std::numeric_limits<double>::epsilon()*std::max(1.*GeV,momentum.mag()),
              "Daughter repair changed the current momentum beyond native floating roundtrip");
      require(child->GetWeight()==weight, "Daughter repair changed its weight");
      require(child->GetDynamicParticle()->GetPreAssignedDecayProperTime()==properTime,
              "Daughter repair changed its assigned clock");
      require(child->GetDynamicParticle()->GetPreAssignedDecayProducts()==products,
              "Daughter repair replaced its decay tree");
      require(source->GetMomentum()==sourceMomentum && source->GetTotalEnergy()==sourceEnergy,
              "Daughter repair changed its source record");
      if (child->GetDynamicParticle()->GetPreAssignedDecayProducts()) {
        const auto kinetic = child->GetKineticEnergy();
        const auto mass = child->GetDynamicParticle()->GetMass();
        const auto direction = child->GetMomentumDirection();
        bool repeated = false;
        try { sim::preservePreassignedDaughterMass(*child); }
        catch (const std::runtime_error&) { repeated = true; }
        require(repeated, "Daughter mass repair was applied twice");
        require(child->GetKineticEnergy()==kinetic && child->GetDynamicParticle()->GetMass()==mass &&
                child->GetMomentumDirection()==direction, "Rejected repeated repair changed the stored state");
      }
    }
    return children;
  }
  double residual(const G4Track& parent, const std::vector<G4Track*>& children) {
    G4LorentzVector sum;
    for(const auto* child:children) sum+=child->GetDynamicParticle()->Get4Momentum();
    const auto delta=sum-parent.GetDynamicParticle()->Get4Momentum();
    return std::max({std::abs(delta.x()),std::abs(delta.y()),std::abs(delta.z()),std::abs(delta.e())})/GeV;
  }
  void closes(const G4Track& parent,const std::vector<G4Track*>& children) {
    require(children.size()==2,"Native transformer lost the two-body frontier");
    std::cerr << "closure PDG=" << parent.GetDefinition()->GetPDGEncoding() << " E=" << parent.GetTotalEnergy()/GeV << " residualGeV=" << residual(parent,children) << " parent=" << parent.GetDynamicParticle()->Get4Momentum()/GeV << "\n";
    require(residual(parent,children)<=1.e-7*std::max(1.,parent.GetTotalEnergy()/GeV),"Native post-loss decay four-momentum mismatch");
  }
  template<class Action> void rejected(Action action) {
    bool failed=false;
    try { action(); } catch(const std::runtime_error&) { failed=true; }
    require(failed,"Invalid or repeated conversion was accepted");
  }
}
int main() {
  try {
    G4MesonConstructor().ConstructParticle();
    G4ShortLivedConstructor().ConstructParticle();
    auto* table=G4ParticleTable::GetParticleTable(); table->SetReadiness();
    // Exact real-source inputs that exposed prior-boost rounding in SIM.
    for(const auto& p : std::vector<G4LorentzVector>{
          {0.002691589156319015,0.0021077224618635975,-0.10421167580075427,0.10426773479540079},
          {-0.012968645193069017,0.00913319433658286,-0.47135335726589983,0.471620173871694},
          {-0.02518569704765688,0.006841913403374407,-1.245545757584181,1.245819154352164},
          {0.01154075301341026,0.0036364985561670674,-0.15562588708756547,0.1560955791611571}}) {
      G4PrimaryParticle particle(22);
      sim::setRecordedPrimaryMomentum(particle,p,0.);
      require(particle.GetMass()==0.,"Rounded massless source gained mass");
      require((particle.GetMomentum()/GeV-p.vect()).mag()<1.e-14,"Source input momentum changed");
    }
    { G4PrimaryParticle particle(22);
      rejected([&]{sim::setRecordedPrimaryMomentum(particle,G4LorentzVector(1.,0.,0.,.5),0.);}); }
    {
      // A real boosted pion whose local energy hides accumulated input error.
      const G4LorentzVector source(.5106473297550549,.24222631765634767,
                                   -23.10845736398936,23.115789318382205);
      G4PrimaryParticle particle(211);
      rejected([&]{sim::setRecordedPrimaryMomentum(particle,source,.13957);});
      sim::setRecordedPrimaryMomentum(particle,source,.13957,6800.);
      require(particle.GetMass()==.13957*GeV,"Event-scale check changed the source mass");
      rejected([&]{sim::setRecordedPrimaryMomentum(particle,G4LorentzVector(1.,0.,0.,.5),0.,6800.);});
    }
    unsigned cases=0;
    double maximumResidual=0.;
    for(bool matched:{false,true}) for(int shape=0;shape<3;++shape) for(int state=0;state<3;++state) {
      std::cerr << "CASE matched=" << matched << " shape=" << shape << " state=" << state << "\n";
      const double dMass=(matched?table->FindParticle(421)->GetPDGMass():1.86486*GeV);
      const double aMass=(matched?table->FindParticle(313)->GetPDGMass():.8925589953787597*GeV);
      const double bMass=(matched?table->FindParticle(223)->GetPDGMass():.7819656374787036*GeV);
      const double b0Mass=table->FindParticle(511)->GetPDGMass();
      const double piMass=table->FindParticle(211)->GetPDGMass();
      G4LorentzVector rootP4(0.,0.,-std::sqrt(std::pow(120.44811847523788*GeV,2)-dMass*dMass),120.44811847523788*GeV);
      G4LorentzVector dP4=rootP4, companion;
      if(shape==1) {
        rootP4=G4LorentzVector(0.,0.,-std::sqrt(std::pow(200.*GeV,2)-b0Mass*b0Mass),200.*GeV);
        auto pair=twoBody(rootP4,b0Mass,dMass,piMass); dP4=pair.first; companion=pair.second;
      }
      auto children=twoBody(dP4,dMass,aMass,bMass);
      auto* d=primary(421,dP4,dMass);
      auto* first=primary(313,children.first,aMass);
      auto* second=primary(223,children.second,bMass);
      if(shape==2) {
        // Native Geant4 promotes this undefined prompt node's known children.
        auto* prompt=primary(90000001,dP4,dMass); prompt->SetProperTime(0.);
        prompt->SetDaughter(first);prompt->SetDaughter(second);d->SetDaughter(prompt);
      } else { d->SetDaughter(first);d->SetDaughter(second); }
      auto* root=d;
      if(shape==1) { root=primary(511,rootP4,b0Mass);root->SetDaughter(d);root->SetDaughter(primary(211,companion,piMass)); }
      G4Event event(cases); auto* vertex=new G4PrimaryVertex(G4ThreeVector(),0.);
      vertex->SetPrimary(root);event.AddPrimaryVertex(vertex);
      G4PrimaryTransformer transformer;auto* tracks=transformer.GimmePrimaries(&event);
      require(tracks->size()==1,"Unexpected native primary count");auto* track=tracks->front();
      track->SetTrackID(1);track->SetWeight(3.75);
      sim::sourceAncestry::reset(true,50000.);
      sim::sourceAncestry::remember(track);
      track->SetParentID(1);rejected([&]{sim::preservePreassignedDecayKinematics(*track);});track->SetParentID(0);
      { G4Track stepped(*track); stepped.IncrementCurrentStepNumber();
        rejected([&]{sim::preservePreassignedDecayKinematics(stepped);}); }
      rejected([&]{sim::preservePreassignedDaughterMass(*track);});
      sim::preservePreassignedDecayKinematics(*track);
      rejected([&]{sim::preservePreassignedDecayKinematics(*track);});
      require(track->GetWeight()==3.75,"Conversion changed the native weight");
      require(track->GetDynamicParticle()->GetPreAssignedDecayProperTime()==.123*ns,"Conversion changed assigned proper time");
      auto result=decay(*track,state);closes(*track,result);maximumResidual=std::max(maximumResidual,residual(*track,result));
      if(shape==1) {
        auto found=std::find_if(result.begin(),result.end(),[](auto* t){return t->GetDynamicParticle()->GetPDGcode()==421;});
        require(found!=result.end(),"Nested D0 missing");
        const auto afterParentDecay = (*found)->GetDynamicParticle()->Get4Momentum();
        { G4Track stepped(**found); stepped.IncrementCurrentStepNumber();
          rejected([&]{sim::preservePreassignedDaughterMass(stepped);}); }
        { G4Track flown(**found); flown.AddTrackLength(1.*mm);
          rejected([&]{sim::preservePreassignedDaughterMass(flown);}); }
        { G4Track aged(**found); aged.SetLocalTime(1.*ns);
          rejected([&]{sim::preservePreassignedDaughterMass(aged);}); }
        require((*found)->GetDynamicParticle()->Get4Momentum()==afterParentDecay,
                "Rejected late repair changed the physical daughter");
        auto grandchildren=decay(**found,(state+1)%3);closes(**found,grandchildren);
        maximumResidual=std::max(maximumResidual,residual(**found,grandchildren));
        for(auto* child:grandchildren)delete child;
      }
      for(auto* child:result)delete child;
      delete track;
      ++cases;
    }
    std::cout<<"PreassignedDecayKinematics native cases="<<cases<<" maximum_four_residual_GeV="<<maximumResidual<<"\n";
    return 0;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
