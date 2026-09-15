#ifndef TrackPropagation_Geant4e_FieldGradientJacobian_h
#define TrackPropagation_Geant4e_FieldGradientJacobian_h

#include "G4ErrorMatrix.hh"
#include "G4ThreeVector.hh"
#include <array>

namespace TrackPropagation {
  // Missing first-order term in the locally uniform-field helix Jacobian:
  // d(delta u)/ds = k q/p u x ((grad B) delta r). Coordinates are
  // (1/p, lambda, phi, xT, yT), B in tesla, distance in cm, p in GeV.
  // qOverP includes the effective reversed charge for backward Geant4 motion.
  inline G4ErrorMatrix fieldGradientJacobianCorrection(
      G4ThreeVector const& direction,
      std::array<G4ThreeVector, 2> const& fieldDerivatives,
      double lengthCm,
      double qOverP) {
    G4ErrorMatrix correction(5, 5, 0);
    double const transverse = direction.perp();
    if (!(transverse > 0.)) return correction;
    G4ThreeVector const u(-direction.y()/transverse, direction.x()/transverse, 0.);
    G4ThreeVector const v = direction.cross(u);
    double const scale = .00299792458 * qOverP * lengthCm;
    for (unsigned int j = 0; j < 2; ++j) {
      G4ThreeVector const change = scale * direction.cross(fieldDerivatives[j]);
      correction(2, j+4) = v.dot(change);
      correction(3, j+4) = u.dot(change)/transverse;
    }
    return correction;
  }
}
#endif
