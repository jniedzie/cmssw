// Code style: AI-code. AI-generated standalone fixed-linear contract tests.
// Experimental mathematics reference for tests; not a production fitter.
#include "IteratedDetectorFit.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
using namespace shift::iterated;

int checks=0, groups=0;
void require(bool condition, char const* message) {
  ++checks;
  if (!condition) throw std::runtime_error(message); // active under -DNDEBUG
}
template<class A,class B> void close(A const& a,B const& b,char const* message) {
  require(a.rows()==b.rows() && a.cols()==b.cols(),message);
  double const error=(a-b).norm();
  require(error<2.e-10*std::max(1.,b.norm()),message);
}
Hit fullHit(std::size_t site,State const& z,Matrix5 const& R) {
  Hit h;h.site=site;h.z=z;h.H=Matrix5::Identity();h.R=R;return h;
}
void singleNode() {
  State z;z<<1.,-2.,3.,-.4,.05;
  Matrix5 R=Matrix5::Identity();R(0,1)=R(1,0)=.2;
  std::vector<Hit> hits{fullHit(0,z,R)};
  auto a=solve({State::Constant(50.)},{},hits);
  auto b=solve({State::Constant(-300.)},{},hits);
  require(a.valid && b.valid,"single node solve");
  close(a.states[0],z,"single mean");close(a.covariances[0],R,"correlated R covariance");
  close(a.states[0],b.states[0],"reference not prior");
  close(a.covariances[0],b.covariances[0],"reference not covariance");
  require(a.objective<1.e-20,"single exact residual");++groups;
}
void affineNoPrior() {
  State ref;ref<<10.,20.,30.,40.,50.;
  Matrix5 F=Matrix5::Identity();F(0,1)=.7;F(2,0)=-.3;
  State offset;offset<<3.,-2.,.5,1.,-3.;
  Edge e;e.F=F;e.mean=F*ref+offset;
  State y0; y0<<.2,-.8,.9,-1.,2.;State y1; y1<<1.,.5,-.4,2.,3.;
  Matrix5 R0=Matrix5::Identity(),R1=Matrix5::Identity()*2.;
  auto fit=solve({ref,State::Constant(-900.)},{e},{fullHit(0,y0,R0),fullHit(1,y1,R1)});
  require(fit.valid && fit.processRanks[0]==0,"zero Q solve");
  Matrix5 information=R0.inverse()+F.transpose()*R1.inverse()*F;
  Matrix5 expectedCov=information.inverse();
  State expected=expectedCov*(R0.inverse()*y0+F.transpose()*R1.inverse()*(y1-offset));
  close(fit.states[0],expected,"affine mean batch oracle");
  close(fit.states[1],F*expected+offset,"exact deterministic edge");
  close(fit.covariances[0],expectedCov,"batch covariance");
  close(fit.crossCovariance(0,1),expectedCov*F.transpose(),"cross covariance");++groups;
}
void rtsRankOne() {
  State ref;ref<<4.,-3.,2.,1.,.5;
  State offset;offset<<.2,.7,-.5,1.,-.3;
  Matrix5 F=Matrix5::Identity();F(0,1)=.4;F(2,0)=-.2;
  State kick;kick<<.2,-.1,.05,0.,0.;Matrix5 Q=kick*kick.transpose();
  Edge e;e.mean=F*ref+offset;e.F=F;e.Q=Q;
  Prior prior;prior.enabled=true;prior.mean=State::Constant(.3);
  prior.covariance=Matrix5::Identity()*2.;prior.covariance(0,1)=prior.covariance(1,0)=.25;
  State y0; y0<<.7,-.1,.4,.8,-1.;State y1;y1<<1.,.2,-.2,-.5,.9;
  Matrix5 R0=Matrix5::Identity();R0(2,3)=R0(3,2)=.1;
  Matrix5 R1=Matrix5::Identity()*1.5;
  auto fit=solve({ref,State::Constant(-500.)},{e},{fullHit(0,y0,R0),fullHit(1,y1,R1)},prior);
  require(fit.valid && fit.processRanks[0]==1,"rank1 process solve");
  Matrix5 K0=prior.covariance*(prior.covariance+R0).inverse();
  State m0=prior.mean+K0*(y0-prior.mean);
  Matrix5 C0=(Matrix5::Identity()-K0)*prior.covariance;
  State pred=F*m0+offset;Matrix5 P=F*C0*F.transpose()+Q;
  Matrix5 K1=P*(P+R1).inverse();State m1=pred+K1*(y1-pred);
  Matrix5 C1=(Matrix5::Identity()-K1)*P;
  Matrix5 G=C0*F.transpose()*P.inverse();
  close(fit.states[1],m1,"filter endpoint mean");
  close(fit.states[0],m0+G*(m1-pred),"RTS source mean");
  close(fit.covariances[1],C1,"RTS endpoint covariance");
  close(fit.covariances[0],C0+G*(C1-P)*G.transpose(),"RTS source covariance");
  close(fit.crossCovariance(0,1),G*C1,"RTS cross covariance");
  auto repeated=solve({ref,State::Zero()},{e},{fullHit(0,y0,R0),fullHit(1,y1,R1)},prior);
  close(fit.states[0],repeated.states[0],"fixed prior repeated mean");
  close(fit.covariances[0],repeated.covariances[0],"fixed prior repeated covariance");++groups;
}
void smallNoiseUnits() {
  Edge e;e.Q(0,0)=1.e-30;e.Q(1,1)=1.e-4;
  auto fit=solve({State::Zero(),State::Zero()},{e},
      {fullHit(0,State::Zero(),Matrix5::Identity()),fullHit(1,State::Zero(),Matrix5::Identity())});
  require(fit.valid && fit.processRanks[0]==2,"tiny distinct physical noise not erased by units");++groups;
}
void failures() {
  Hit h;h.site=0;h.z=Eigen::VectorXd::Zero(1);h.H=Eigen::MatrixXd::Zero(1,5);h.H(0,0)=1.;h.R=Eigen::MatrixXd::Identity(1,1);
  auto rank=solve({State::Zero()},{},{h});
  require(!rank.valid && rank.status==unobservable && rank.states.empty(),"unobservable no partial fit");
  Edge e;e.Q(0,0)=-1.e-30;
  auto neg=solve({State::Zero(),State::Zero()},{e},{fullHit(0,State::Zero(),Matrix5::Identity())});
  require(!neg.valid && neg.status==invalidNoise,"negative process diagonal");
  e.Q.setZero();e.Q(0,0)=e.Q(1,1)=1.;e.Q(0,1)=e.Q(1,0)=2.;
  require(solve({State::Zero(),State::Zero()},{e},{fullHit(0,State::Zero(),Matrix5::Identity())}).status==invalidNoise,"indefinite process covariance");
  e.Q.setZero();e.Q(0,1)=e.Q(1,0)=1.e-20;
  require(solve({State::Zero(),State::Zero()},{e},{fullHit(0,State::Zero(),Matrix5::Identity())}).status==invalidNoise,"nullspace cross noise");
  h=fullHit(0,State::Zero(),Matrix5::Identity());h.valid=false;
  require(solve({State::Zero()},{},{h}).status==invalidMeasurement,"invalid hit not skipped");
  Options limits;limits.maximumNodes=1;e=Edge{};
  require(solve({State::Zero(),State::Zero()},{e},{fullHit(0,State::Zero(),Matrix5::Identity())},Prior{},limits).status==budgetExceeded,"bounded core");++groups;
}
int main() {
  try { singleNode();affineNoPrior();rtsRankOne();smallNoiseUnits();failures(); }
  catch(std::exception const& e) {std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}
  std::cout<<"PASS groups="<<groups<<" checks="<<checks<<" active_under_NDEBUG=true\n";
}
