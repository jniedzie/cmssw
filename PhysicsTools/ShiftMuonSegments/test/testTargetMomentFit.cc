#include "PhysicsTools/ShiftMuonSegments/interface/TargetMomentFit.h"
#include <Eigen/LU>
#include <stdexcept>

int main() {
  using namespace shift;
  auto check=[](bool x){if(!x)throw std::runtime_error("Target moment fit regression");};
  TargetCovariance j=TargetCovariance::Identity();j(3,1)=-14000.;j(4,2)=-14000.;
  TargetCovariance measurement=TargetCovariance::Identity();measurement(0,0)=1.e-6;
  measurement(1,1)=measurement(2,2)=1.e-8;
  TargetCovariance noise=.3*measurement;
  TargetParameters detector;detector<<.01,.01,-.02,-145.,282.;
  auto transport=[&](auto const& x){return ForwardTargetPrediction{true,j*x,noise};};
  TargetCovariance const weight=(measurement+noise).inverse();
  TargetCovariance normal=j.transpose()*weight*j;normal(3,3)+=100.;normal(4,4)+=100.;
  TargetCovariance const expectedCovariance=normal.inverse();
  TargetParameters const expected=expectedCovariance*j.transpose()*weight*detector;
  TargetParameters seed=j.inverse()*detector;
  for (double scale:{1.,2.}) {
    auto changed=seed;changed[0]*=scale;
    auto fit=fitTargetMoments(changed,detector,measurement,.1,.1,0.,transport);
    check(fit.valid && (fit.parameters-expected).norm()<1.e-6 && (fit.covariance-expectedCovariance).norm()<1.e-8);
  }
  check(!fitTargetMoments(seed,detector,measurement,.1,.1,0.,[](auto const&){return ForwardTargetPrediction{};}).valid);
  unsigned int calls=0;
  auto analytic=[&](TargetParameters const& x) {
    ++calls;ForwardTargetPrediction p{true,j*x,noise};p.jacobianValid=true;p.jacobian=j;return p;
  };
  auto analyticFit=fitTargetMoments(expected,detector,measurement,.1,.1,0.,analytic);
  check(analyticFit.valid && (analyticFit.parameters-expected).norm()<1.e-6 &&
        (analyticFit.covariance-expectedCovariance).norm()<1.e-8 && calls==1);
  // A linear problem must converge from a distant seed; a step cap in
  // posterior sigmas incorrectly makes convergence depend on the guess.
  auto distant=seed;distant[1]+=.1;distant[2]-=.1;distant[3]=30.;distant[4]=100.;
  auto distantFit=fitTargetMoments(distant,detector,measurement,.1,.1,0.,analytic);
  check(distantFit.valid && distantFit.iterations<=3 &&
        (distantFit.parameters-expected).norm()<1.e-6 &&
        (distantFit.covariance-expectedCovariance).norm()<1.e-8);
  // A known mean with state-dependent variance: the quasi-score solution
  // is the observed mean. It deliberately differs from a Gaussian-density
  // maximum which uses the variance itself to infer the parameter.
  detector.setZero();detector[0]=.02;measurement.setIdentity();measurement(0,0)=1.e-6;
  auto changing=[](TargetParameters const& x) {
    ForwardTargetPrediction p{true,x,TargetCovariance::Zero()};p.noise(0,0)=.2*x[0]*x[0];return p;
  };
  seed=detector;seed[0]=.03;
  auto fit=fitTargetMoments(seed,detector,measurement,.1,.1,0.,changing);
  check(fit.valid && std::abs(fit.parameters[0]-.02)<1.e-8);
  auto nonlinear=[](TargetParameters const& x) {
    ForwardTargetPrediction p{true,x,TargetCovariance::Zero()};p.parameters[0]=x[0]*x[0];
    p.noise(0,0)=.01*std::pow(x[0],4);return p;
  };
  detector[0]=.0004;measurement(0,0)=1.e-10;seed[0]=.035;
  fit=fitTargetMoments(seed,detector,measurement,.1,.1,0.,nonlinear);
  check(fit.valid && std::abs(fit.parameters[0]-.02)<2.e-5);
}
