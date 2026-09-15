#include "TrackPropagation/Geant4e/interface/FieldGradientJacobian.h"
#include <cmath>
#include <stdexcept>

int main() {
  // A straight on-axis quadrupole has zero nominal field. Its transverse
  // variational equation nevertheless focuses one plane and defocuses the
  // other. Compare the integrated Jacobian with the closed-form cos/sin or
  // cosh/sinh transfer matrix, for both charges and propagation directions.
  constexpr double slope = 1.e-5, gradient = .1, momentum = 50., length = 500.;
  constexpr unsigned int steps = 20000;
  double const ds = length/steps;
  for (double sign : {-1., 1.}) for (double charge : {-1., 1.}) {
    G4ThreeVector const direction(slope, 0., sign*std::sqrt(1.-slope*slope));
    G4ThreeVector const u(0.,1.,0.), v = direction.cross(u);
    std::array<G4ThreeVector,2> const derivatives{
      G4ThreeVector(gradient,0.,0.),G4ThreeVector(0.,gradient*v.x(),0.)};
    auto correction = TrackPropagation::fieldGradientJacobianCorrection(direction,derivatives,ds,charge/momentum);
    G4ErrorMatrix step(5,5,1), total(5,5,1);
    step(4,3)=ds*slope; step(5,2)=ds;
    step+=correction;
    for (unsigned int j=0;j<steps;++j) total=step*total;
    double const strength=.00299792458*charge*sign*gradient/momentum;
    double const w=std::sqrt(std::abs(strength));
    double const c=strength>0.?std::cos(w*length):std::cosh(w*length);
    double const s=(strength>0.?std::sin(w*length):std::sinh(w*length))/w;
    if (std::abs(total(5,5)-c)>5.e-4 || std::abs(total(5,2)-s)/length>5.e-4 ||
        std::abs(total(2,5)+strength*s)*length>5.e-4 || std::abs(total(2,2)-c)>5.e-4)
      throw std::runtime_error("Quadrupole gradient transport does not match analytic focusing");
    std::array<G4ThreeVector,2> zero{};
    auto const uniform=TrackPropagation::fieldGradientJacobianCorrection(direction,zero,ds,charge/momentum);
    for (int i=1;i<=5;++i) for (int j=1;j<=5;++j)
      if (uniform(i,j)!=0.) throw std::runtime_error("Uniform field must have no gradient correction");
  }
}
