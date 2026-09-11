#include "PhysicsTools/ShiftMuonSegments/interface/TargetAngularError.h"
#include <cassert>
#include <cmath>
int main() {
  AlgebraicSymMatrix55 c;
  c(1,1)=4.e-6; c(2,2)=9.e-6; c(1,2)=1.e-6;
  for (double sign : {-1., 1.}) {
    double const x=.03, y=-.02, h=1.e-7;
    auto eta=[sign](double a,double b) {return std::asinh(sign/std::hypot(a,b));};
    auto phi=[sign](double a,double b) {return std::atan2(sign*b,sign*a);};
    auto out=shift::targetAngularError(x,y,sign,c);
    double ex=(eta(x+h,y)-eta(x-h,y))/(2*h), ey=(eta(x,y+h)-eta(x,y-h))/(2*h);
    double px=(phi(x+h,y)-phi(x-h,y))/(2*h), py=(phi(x,y+h)-phi(x,y-h))/(2*h);
    double ve=ex*ex*c(1,1)+2*ex*ey*c(1,2)+ey*ey*c(2,2);
    double vp=px*px*c(1,1)+2*px*py*c(1,2)+py*py*c(2,2);
    assert(std::abs(out.eta* out.eta/ve-1.) < 1.e-7);
    assert(std::abs(out.phi* out.phi/vp-1.) < 1.e-7);
  }
  assert(shift::targetAngularError(0.,0.,1.,c).eta < 0.);
}
