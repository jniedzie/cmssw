#include "PhysicsTools/ShiftMuonSegments/interface/ShortFlowDerivative.h"
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace shift;
void require(bool value,char const* text) {if (!value) throw std::runtime_error(text);}

int main() {
  // An exactly integrated smooth local flow: q/p'=loss, tx'=field*q/p,
  // x'=tx, y'=ty. Nonzero loss tests the curvature part of -J*f, not only
  // straight-line position derivatives. This is a mathematical fixture,
  // not a physical loss model or an empirical reconstruction correction.
  for(double field:{-.0003,0.,.0003}) for(double loss:{-1.e-7,0.,1.e-7})
    for(double sourceSign:{-1.,1.}) for(double charge:{-1.,1.}) {
      double const z0=sourceSign*14800.,zend=-sourceSign*800.;
      FlowVector initial;initial<<charge*.012,.015,-.009,2.,-3.;
      auto flow=[&](FlowVector const& u,double length) {
        FlowVector out=u;
        out[0]+=loss*length;
        out[1]+=field*(u[0]*length+.5*loss*length*length);
        out[3]+=u[1]*length+.5*field*u[0]*length*length+
                field*loss*length*length*length/6.;
        out[4]+=u[2]*length;
        return out;
      };
      double const length=zend-z0;
      FlowMatrix jacobian=FlowMatrix::Identity();
      jacobian(1,0)=field*length;
      jacobian(3,0)=.5*field*length*length;
      jacobian(3,1)=jacobian(4,2)=length;
      FlowMatrix weight=FlowMatrix::Identity();
      weight(0,0)=1.e10;weight(1,1)=weight(2,2)=1.e6;
      auto probe=[&](double requested) {
        // Reproduce representable detector-library plane coordinates.
        double const actual=double(float(z0+requested))-z0;
        return FlowProbe{true,actual,flow(initial,actual)};
      };
      auto const estimate=shortFlowDerivative(jacobian,weight,1.,.01,probe);
      require(estimate.valid && estimate.calls==4,"smooth flow rejected");
      FlowVector const full=(flow(initial,zend-(z0+1.))-flow(initial,zend-(z0-1.)))/2.;
      FlowVector const half=(flow(initial,zend-(z0+.5))-flow(initial,zend-(z0-.5)));
      require((estimate.derivative-half).norm()<1.e-9,"short-flow/full-leg derivative mismatch");
      require((full-half).norm()<1.e-9,"independent full-leg half-step mismatch");
  }
  // Starting exactly at a material interface: the central average is stable
  // under halving, but the two one-sided derivatives disagree. Reject it.
  FlowMatrix const unit=FlowMatrix::Identity();
  auto boundary=[](double h) {
    FlowVector p=FlowVector::Zero();p[0]=.01+(h>0. ? 1.e-4 : 3.e-4)*h;
    return FlowProbe{true,h,p};
  };
  auto const edge=shortFlowDerivative(unit,unit,1.,.01,boundary);
  require(!edge.valid && edge.centralRelativeError<1.e-9 && edge.oneSidedRelativeError>.4,
          "material-edge central-difference false convergence not rejected");
  auto collapsed=[](double h) {return FlowProbe{true,std::copysign(1.,h),FlowVector::Zero()};};
  require(!shortFlowDerivative(unit,unit,1.,.01,collapsed).valid,"collapsed float offsets accepted");
  auto failed=[](double) {return FlowProbe{};};
  auto const failure=shortFlowDerivative(unit,unit,1.,.01,failed);
  require(!failure.valid && failure.calls==1,"failed short transport ignored");
  std::cout<<"Short-flow identity, full-leg half-step closure, both signs/charges, material-edge and precision gates: PASS\n";
}
