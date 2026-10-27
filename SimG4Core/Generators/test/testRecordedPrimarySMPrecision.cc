// Code style: AI-code. Exact SM input failures and strict malformed-input controls.
#include "SimG4Core/Generators/interface/RecordedPrimaryMomentum.h"
#include "G4BosonConstructor.hh"
#include "G4ParticleTable.hh"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
struct P4 { double x,y,z,t;double px()const{return x;}double py()const{return y;}double pz()const{return z;}double e()const{return t;} };
struct Case {P4 p;double mass;};
void require(bool condition,const char* reason){if(!condition)throw std::runtime_error(reason);}
bool accepted(P4 p,double m,double scale){G4PrimaryParticle g(22);try{sim::setRecordedPrimaryMomentum(g,p,m,scale);return true;}catch(const std::runtime_error&){return false;}}
int main(int argc,char**argv){
 const bool legacy=argc==2&&std::string(argv[1])=="--expect-legacy";
 G4BosonConstructor().ConstructParticle();G4ParticleTable::GetParticleTable()->SetReadiness();
 const std::vector<Case> cases={
  {{0x1.802874db30438p-3,-0x1.16875398621acp-3,-0x1.1390bf2d49258p+4,0x1.1399722ae1e64p+4},0x1.1dd6e04c05921p-3},
  {{0x1.5d17202741207p-2,-0x1.5d86f1d61b114p-5,-0x1.49237acdba415p+3,0x1.495a39c83c8d6p+3},0x1.1dd6e04c05921p-3},
  {{0x1.57fcf7161a0f8p-2,-0x1.7b8c330f4c78ap-5,-0x1.47fe548dae730p+3,0x1.4830bacd2fc4dp+3},0x1.b0c88a47ecfeap-4},
  {{0x1.02e78e4ce2fd0p+0,0x1.886666b80f3d2p-3,-0x1.6eca7d9f1f740p+5,0x1.6ef5ce78263b2p+5},0x1.e064ece9a2c67p-1},
  {{-0x1.72e9696640ca8p-1,0x1.679a3d6bf3ff4p-2,-0x1.05edc292f0954p+4,0x1.06acb4a0722cfp+4},0x1.e064ece9a2c67p-1},
  {{0x1.8c78ba5e8e435p-4,-0x1.7f058b50881fdp-2,-0x1.7eedac69aeb91p+4,0x1.7efc01703f697p+4},0x1.147064ece9a2cp-3},
  {{-0x1.60ac6bd39d560p-6,-0x1.291daf40b4508p-4,-0x1.6902b3c2a02f7p+1,0x1.6923f43ee2713p+1},0x0.0p+0},
  {{0x1.e4a3d5537598ep-4,-0x1.34be1f805b0bbp-2,-0x1.51cd55f15ab32p+4,0x1.51d782e8631b5p+4},0x0.0p+0},
  {{-0x1.1fc1e2f0d8692p-3,-0x1.83b5b975c07acp-4,-0x1.0f9ae386754b6p+2,0x1.0ff5675b64efap+2},0x1.147064ece9a2cp-3},
 };
 for(const auto& row:cases){
  require(accepted(row.p,row.mass,6800.)!=legacy,"Exact SM failure acceptance differs");
 }
 // These good vectors must keep the exact same setter outcome in old and fixed code.
 const std::vector<Case> good={{{0.,0.,2.,2.},0.},{{.3,.4,-100.,std::hypot(std::hypot(.3,.4,100.),.13957)},.13957},{{0.,0.,0.,.105658},.105658}};
 std::cout<<std::setprecision(17)<<"{\"saved_failures\":"<<cases.size()<<",\"legacy\":"<<(legacy?"true":"false")<<",\"good_converted_p4\":[";
 bool first=true;
 for(const auto& row:good){
  require(accepted(row.p,row.mass,6800.),"Good source rejected");
  G4PrimaryParticle g(22);sim::setRecordedPrimaryMomentum(g,row.p,row.mass,6800.);
  if(!first)std::cout<<",";first=false;
  std::cout<<"["<<g.GetPx()/GeV<<","<<g.GetPy()/GeV<<","<<g.GetPz()/GeV<<","<<g.GetTotalEnergy()/GeV<<","<<g.GetMass()/GeV<<"]";
 }
 const double nan=std::numeric_limits<double>::quiet_NaN();
 require(!accepted({1.,0.,0.,.5},0.,6800.),"Genuine spacelike source accepted");
 require(!accepted({1.,0.,0.,1.+2.e-7},0.,6800.),"Conservation cap exceeded");
 require(!accepted({0.,0.,0.,1.},2.,6800.),"Energy below mass accepted");
 require(!accepted({0.,0.,0.,1.},-1.,6800.),"Negative mass accepted");
 require(!accepted({nan,0.,0.,1.},0.,6800.),"Nonfinite momentum accepted");
 std::cout<<"],\"negative_controls\":5,\"passed\":true}\n";
}
