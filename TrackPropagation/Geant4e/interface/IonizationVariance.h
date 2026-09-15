#ifndef TrackPropagation_Geant4e_IonizationVariance_h
#define TrackPropagation_Geant4e_IonizationVariance_h
#include <algorithm>
#include <cmath>

namespace TrackPropagation {
  // Unrestricted second moment of ionization energy transfer, in GeV^2.
  // This is the same expression as G4ErrorFreeTrajState::PropagateErrorIoni
  // before its per-step kappa suppression. It describes the second moment,
  // including rare delta rays, rather than the width of a Gaussian core.
  struct IonizationVariance { double physical = 0.; double legacy = 0.; };
  inline IonizationVariance ionizationVariance(double totalEnergyGeV, double massGeV,
                                               double zOverA, double densityGPerCm3, double lengthCm) {
    if (!(totalEnergyGeV > massGeV && massGeV > 0. && lengthCm > 0.)) return {};
    double const gamma = totalEnergyGeV/massGeV;
    double const beta2 = 1.-1./(gamma*gamma);
    constexpr double electronMassGeV = .00051099906;
    double const ratio = electronMassGeV/massGeV;
    double const maximumTransfer = 2.*electronMassGeV*(gamma*gamma-1.) /
        (1.+2.*ratio*gamma+ratio*ratio);
    double const xi = .0001535*zOverA*densityGPerCm3*lengthCm/beta2;
    double const variance = xi*maximumTransfer*(1.-.5*beta2);
    return {variance, variance*std::min(1.,100.*xi/maximumTransfer)};
  }
}
#endif
