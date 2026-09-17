#include "PhysicsTools/ShiftMuonSegments/interface/MultiscaleDerivativeAudit.h"
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
  using namespace shift;
  void check(bool condition,char const* message) {
    if (!condition) throw std::runtime_error(message);
  }
}
int main() {
  DerivativeAuditVector derivative;derivative<<1.,2.,-3.,.5,4.;
  DerivativeAuditMatrix factor=DerivativeAuditMatrix::Identity();
  factor.diagonal()<<2.,.2,10.,4.,.7;factor(3,1)=.2;
  DerivativeAuditMatrix const weight=factor*factor.transpose();
  unsigned int calls=0;
  auto linear=[&](double offset) {
    ++calls;return DerivativeAuditProbe{true,derivative*offset,offset};
  };
  auto result=auditDerivativeMultiscale(derivative,weight,linear,.1);
  check(result.valid && result.resolved && result.status==DerivativeAuditResult::passed &&
        result.calls==4 && calls==4 && result.attempts.size()==1 && result.chosenStep==.1 &&
        (result.reference-derivative).norm()<1.e-12,"exact derivative did not pass first scale");
  std::cout<<"Exact derivative and fixed four-call first attempt: PASS\n";

  // A cubic derivative resolves first at scale .125. Earlier scales are
  // explicitly attempted in the declared order, regardless of their size.
  calls=0;
  auto nonlinear=[&](double offset) {
    ++calls;return DerivativeAuditProbe{true,derivative*(offset+.1*offset*offset*offset),offset};
  };
  result=auditDerivativeMultiscale(derivative,weight,nonlinear,1.);
  check(result.valid && result.chosenStep==.125 && result.attempts.size()==7 && result.calls==28 && calls==28,
        "nonlinear reference did not select first stable prescribed scale");
  for (unsigned int i=0;i<6;++i) check(!result.attempts[i].stable && result.attempts[i].disagreement<0.,
                                    "analytic comparison influenced unstable scale selection");
  std::cout<<"Nonlinear reference follows prescribed scales and independent stability: PASS\n";

  // A float-backed source has ULP .0625 at 1e6: small input steps collapse,
  // and some coarse/half probes land on the SAME representable coordinate.
  // Only the last requested scale has four distinct usable coordinates.
  // The float-backed means are evaluated from those actual coordinates.
  calls=0;
  auto quantized=[&](double offset) {
    ++calls;
    float const source=1000000.f;
    float const actualSource=float(double(source)+offset);
    double const actual=double(actualSource)-source;
    DerivativeAuditVector mean;
    for(unsigned int i=0;i<5;++i) mean[i]=float(1024.+derivative[i]*actual);
    return DerivativeAuditProbe{true,mean,actual};
  };
  result=auditDerivativeMultiscale(derivative,weight,quantized,.01);
  check(result.valid && result.chosenStep==.16 && result.attempts.size()==8 && result.calls<=32 && calls==result.calls &&
        (result.reference-derivative).norm()<1.e-12,"float quantization/actual-offset handling failed");
  check(result.attempts.back().actualOffsets[0]!=result.chosenStep,"quantized fixture did not exercise actual denominator");
  std::cout<<"Float quantization, collapsed probes and actual input separations: PASS\n";

  calls=0;
  auto biased=auditDerivativeMultiscale(1.1*derivative,weight,linear,.1);
  check(!biased.valid && biased.resolved && biased.status==DerivativeAuditResult::stableMismatch &&
        biased.calls==4 && calls==4 && biased.attempts.size()==1 && std::abs(biased.disagreement-.1)<1.e-12,
        "stable but biased analytic derivative triggered a forbidden later-scale search");
  std::cout<<"Stable numerical reference with biased analytic derivative stops immediately: PASS\n";

  calls=0;
  auto unresolved=[&](double offset) {
    ++calls;return DerivativeAuditProbe{true,derivative*(offset+1.e6*offset*offset*offset),offset};
  };
  result=auditDerivativeMultiscale(derivative,weight,unresolved,.1);
  check(!result.valid && !result.resolved && result.status==DerivativeAuditResult::unresolvedReference &&
        result.attempts.size()==8 && result.calls==32 && calls==32,"unresolved reference exceeded deterministic budget");
  for(auto const& attempt:result.attempts) check(!attempt.stable && attempt.disagreement<0.,"unresolved reference compared to analytic J");
  std::cout<<"No stable reference stays unresolved at exact 32-call bound: PASS\n";

  auto wrongSign=[&](double offset){return DerivativeAuditProbe{true,derivative*offset,-offset};};
  result=auditDerivativeMultiscale(derivative,weight,wrongSign,.1);
  check(!result.resolved && result.calls==8,"wrong-sign offsets accepted");
  auto collapsed=[&](double offset){return DerivativeAuditProbe{true,derivative*offset,std::copysign(.1,offset)};};
  result=auditDerivativeMultiscale(derivative,weight,collapsed,.1);
  check(!result.resolved && result.calls==32,"collapsed coarse/half offsets accepted");
  auto nonfinite=[&](double offset){auto mean=derivative;mean[0]=std::numeric_limits<double>::quiet_NaN();return DerivativeAuditProbe{true,mean,offset};};
  result=auditDerivativeMultiscale(derivative,weight,nonfinite,.1);
  check(!result.resolved && result.calls==8,"nonfinite callback accepted");
  auto constant=[](double offset){return DerivativeAuditProbe{true,DerivativeAuditVector::Ones(),offset};};
  result=auditDerivativeMultiscale(derivative,weight,constant,.1);
  check(!result.resolved && result.calls==32,"zero-information relative reference accepted");
  auto badWeight=weight;badWeight(0,0)=-1.;
  check(auditDerivativeMultiscale(derivative,badWeight,linear,.1).calls==0,"indefinite metric invoked callback");
  check(auditDerivativeMultiscale(derivative,weight,linear,0.).calls==0,"zero initial step invoked callback");
  auto badAnalytic=derivative;badAnalytic[0]=std::numeric_limits<double>::infinity();
  check(auditDerivativeMultiscale(badAnalytic,weight,linear,.1).calls==0,"nonfinite analytic derivative accepted");
  std::cout<<"Invalid, wrong-sign, collapsed, nonfinite and zero-information probes: PASS\n";
  std::cout<<"Multiscale derivative audit regression suite passed\n";
}
