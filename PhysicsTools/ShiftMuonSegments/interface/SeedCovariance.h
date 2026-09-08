#ifndef PhysicsTools_ShiftMuonSegments_SeedCovariance_h
#define PhysicsTools_ShiftMuonSegments_SeedCovariance_h

#include <cmath>
#include <stdexcept>
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"

namespace shift {
  // Parameters are (q/p, lambda, phi, xT, yT). Reversing both p and q
  // preserves the geometrical helix but changes the curvilinear frame.
  inline AlgebraicSymMatrix55 seedCovariance(AlgebraicSymMatrix55 const& covariance,
                                             double directionSign,
                                             double momentumScale = 1.) {
    if ((directionSign != 1. && directionSign != -1.) || !(momentumScale > 0.) || !std::isfinite(momentumScale))
      throw std::invalid_argument("Invalid seed frame transformation");
    double const jacobian[5] = {directionSign / momentumScale, directionSign, 1., directionSign, 1.};
    auto result = covariance;
    for (unsigned int row = 0; row < 5; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        result(row, column) *= jacobian[row] * jacobian[column];
    return result;
  }
}  // namespace shift
#endif
