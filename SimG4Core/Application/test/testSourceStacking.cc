// Code style: AI-code. AI-generated actual ClassifyNewTrack corridor controls.
// IDs are assigned by this fixture; real event-manager ID ordering remains a SIM gate.
#include "SimG4Core/Application/interface/StackingAction.h"
#include "SimG4Core/Generators/interface/SourceAncestry.h"
#include "SimG4Core/Physics/interface/CMSG4TrackInterface.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "G4Box.hh"
#include "G4Decay.hh"
#include "G4DynamicParticle.hh"
#include "G4GeometryManager.hh"
#include "G4LeptonConstructor.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleTable.hh"
#include "G4Region.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
  void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
}
int main(int argc,char** argv) {
  try {
    require(argc==2,"Need the synthetic stacking parameter fixture");
    const auto config=edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("stack");
    require(!config.getParameter<bool>("TrackNeutrino"),"Fixture must retain the installed neutrino guard");
    G4LeptonConstructor().ConstructParticle();
    G4ParticleTable::GetParticleTable()->SetReadiness();
    auto* vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    auto* worldLogical=new G4LogicalVolume(new G4Box("synthetic_world",1.*km,1.*km,1.*km),vacuum,"synthetic_world");
    auto* world=new G4PVPlacement(nullptr,{},worldLogical,"synthetic_world",nullptr,false,0);
    auto* normal=new G4Region("synthetic_normal");normal->AddRootLogicalVolume(worldLogical);
    auto* deadLogical=new G4LogicalVolume(new G4Box("synthetic_dead",20.*m,20.*m,20.*m),vacuum,"synthetic_dead");
    new G4PVPlacement(nullptr,{},deadLogical,"synthetic_dead",worldLogical,false,0);
    auto* dead=new G4Region("synthetic_dead");dead->AddRootLogicalVolume(deadLogical);
    auto* nav=G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    nav->SetWorldVolume(world);G4GeometryManager::GetInstance()->CloseGeometry();
    StackingAction on(config,nullptr,true),off(config,nullptr,false);G4Decay decay;
    const double boundary=config.getParameter<double>("MaxZCentralCMS")*m;
    auto makeTrack=[&](int id,int parent,int pdg,double z,double time,double pz) {
      auto result=std::make_unique<G4Track>(new G4DynamicParticle(G4ParticleTable::GetParticleTable()->FindParticle(pdg),{0.,0.,pz},40.*GeV),time,G4ThreeVector(0.,0.,z));
      result->SetTrackID(id);result->SetParentID(parent);
      if(parent)result->SetCreatorProcess(&decay);
      nav->LocateGlobalPointAndSetup(result->GetPosition());
      result->SetTouchableHandle(nav->CreateTouchableHistoryHandle());return result;
    };
    int checked=0;
    for(bool active:{false,true})for(double rootZ:{-148.*m,0.,148.*m})for(double sign:{-1.,1.}) {
      for(bool corridor:{false,true}) {
        sim::sourceAncestry::reset(active,boundary);
        auto root=makeTrack(1,0,13,rootZ,0.,sign);auto& action=corridor?on:off;
        require(action.ClassifyNewTrack(root.get())==fUrgent,"Root unexpectedly rejected");
        CMSG4TrackInterface::instance()->setCurrentTrack(root.get());
        auto child=makeTrack(2,1,13,10.*m,0.,1.);
        const bool protectedSource=active && corridor && std::abs(rootZ)>=boundary;
        require(action.ClassifyNewTrack(child.get())==(protectedSource?fUrgent:fKill),"Central dead-region source corridor decision differs");
        auto normalChild=makeTrack(3,1,13,30.*m,0.,1.);
        require(action.ClassifyNewTrack(normalChild.get())==fUrgent,"Normal region child changed");
        auto late=makeTrack(4,1,13,10.*m,1.e12*ns,1.);
        require(action.ClassifyNewTrack(late.get())==fKill,"Central time guard changed");
        auto lateForward=makeTrack(5,1,13,100.*m,1.e12*ns,1.);
        require(action.ClassifyNewTrack(lateForward.get())==fKill,"Forward time guard changed");
        auto neutrino=makeTrack(6,1,14,10.*m,0.,1.);
        require(action.ClassifyNewTrack(neutrino.get())==fKill,"Neutrino guard changed");++checked;
      }
    }
    sim::sourceAncestry::reset(false,boundary);
    std::cout<<"Actual stacking: "<<checked<<" corridor/source/direction cases plus normal-region, time and neutrino controls passed. IDs were fixture-assigned.\n";return 0;
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
