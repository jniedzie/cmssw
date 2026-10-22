// Code style: AI-code. AI-generated private generative iteration contract tests.
// Experimental mathematics reference for tests; not a production fitter.
#include "IteratedDetectorFitIteration.h"
#include <iostream>
#include <stdexcept>
using namespace shift::iterated;
int checks=0,groups=0;
void require(bool value,char const* message) {++checks;if(!value)throw std::runtime_error(message);}
Hit full(std::size_t site,State const& z) {Hit h;h.site=site;h.z=z;h.H=Matrix5::Identity();h.R=Matrix5::Identity();return h;}
void deterministicLinear() {
  State offset=State::Constant(.3),y0=State::Constant(1.),y1=State::Constant(3.);
  Matrix5 F=Matrix5::Identity()*2.;
  Mean mean=[&](std::size_t,State const& x,State& result,std::string&){result=F*x+offset;return true;};
  Linearize model=[&](States const& refs,std::vector<Edge>& edges,std::vector<Hit>& hits,std::string&){Edge e;e.mean=F*refs[0]+offset;e.F=F;edges={e};hits={full(0,y0),full(1,y1)};return true;};
  State initial=State::Constant(-5.);State next=F*initial+offset;
  auto fit=iterate({initial,next},model,mean);
  require(fit.valid,"linear convergence");require(fit.iterations==2,"linear one update then fixed point");
  State expected=(y0+2.*(y1-offset))/5.;
  require((fit.states[0]-expected).norm()<1.e-10,"linear same-data mean");
  require((fit.states[1]-F*fit.states[0]-offset).norm()<1.e-12,"exact deterministic closure");
  require((fit.fit.covariances[0]-Matrix5::Identity()/5.).norm()<1.e-10,"no-prior covariance");
  auto again=iterate(fit.states,model,mean);
  require(again.valid,"repeat converges");
  require((fit.fit.covariances[0]-again.fit.covariances[0]).norm()<1.e-12,"no posterior accumulation");
  require(again.meanEvaluations==1,"one fixed generative check");++groups;
}
void nonlinearSingular(bool noisy) {
  State y0=State::Constant(.4),y1=State::Constant(.8);y0[0]=1.8;y1[0]=4.;
  Matrix5 Q=Matrix5::Zero();if(noisy)Q(1,1)=.1;
  Mean mean=[](std::size_t,State const& x,State& result,std::string&error){
    if(x[0]<=0.) {error="fixed positive branch";return false;}
    result=x;result[0]=x[0]*x[0];return true;};
  Linearize model=[&](States const& refs,std::vector<Edge>& edges,std::vector<Hit>& hits,std::string&error){
    State prediction;if(!mean(0,refs[0],prediction,error))return false;
    Edge e;e.mean=prediction;e.F=Matrix5::Identity();e.F(0,0)=2.*refs[0][0];e.Q=Q;
    edges={e};hits={full(0,y0),full(1,y1)};return true;};
  State first=State::Constant(-.2);first[0]=2.;State next;std::string error;mean(0,first,next,error);if(noisy)next[1]+=.3;
  auto fit=iterate({first,next},model,mean);
  require(fit.valid,"nonlinear convergence");
  require(fit.stateStep<=1.e-3 && fit.momentChange<=1.e-3,"full-step and moment convergence");
  State predicted;mean(0,fit.states[0],predicted,error);
  require(std::abs(fit.states[1][0]-predicted[0])<1.e-12,"nonlinear deterministic nullspace");
  for(int i=2;i<5;++i)require(std::abs(fit.states[1][i]-predicted[i])<1.e-12,"other deterministic directions");
  require(fit.kicks[0].size()==(noisy?1:0),"rectangular process kick");
  require(fit.iterations<=8,"bounded outer work");
  for(std::size_t i=1;i<fit.history.size();++i)
    require(fit.history[i].objective<=fit.history[i-1].objective+1.e-12,"fixed factors descend");
  State other=State::Constant(1.5);other[0]=1.5;State otherNext;mean(0,other,otherNext,error);
  auto second=iterate({other,otherNext},model,mean);
  require(second.valid,"second initialization converges");
  require((fit.states[0]-second.states[0]).norm()<1.e-3,"initialization only");++groups;
}
void failureAccounting() {
  Mean mean=[](std::size_t,State const& x,State& y,std::string&){y=x;return true;};
  Linearize model=[](States const& refs,std::vector<Edge>& edges,std::vector<Hit>& hits,std::string&){Edge e;e.mean=refs[0];edges={e};hits={full(0,State::Zero()),full(1,State::Zero())};return true;};
  auto unsupported=iterate({State::Zero(),State::Ones()},model,mean);
  require(!unsupported.valid && unsupported.status==-11,"unsupported null residual fails");
  IterationOptions options;options.maximumOuter=1;
  auto capped=iterate({State::Ones(),State::Ones()},model,mean,options);
  require(!capped.valid && capped.status==-15,"nonconvergence explicit");
  Linearize failed=[](States const&,std::vector<Edge>&,std::vector<Hit>&,std::string&error){error="missing required edge";return false;};
  auto missing=iterate({State::Zero(),State::Zero()},failed,mean);
  require(!missing.valid && missing.reason.find("missing required edge")!=std::string::npos,"failed callback retained");++groups;
}
int main() {
  try{deterministicLinear();nonlinearSingular(false);nonlinearSingular(true);failureAccounting();}
  catch(std::exception const&error){std::cerr<<"FAIL "<<error.what()<<"\n";return 1;}
  std::cout<<"PASS groups="<<groups<<" checks="<<checks<<" active_under_NDEBUG=true\n";
}
