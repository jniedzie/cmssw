#include "TrackPropagation/Geant4e/interface/ConvertFromToCLHEP.h"
#include <stdexcept>

void check(bool value) {
  if (!value) throw std::runtime_error("Backward covariance frame regression");
}
#include <cmath>

int main() {
  // A correlated positive-definite covariance makes a missed frame reversal
  // visible; testing only variances cannot detect it.
  AlgebraicSymMatrix55 covariance;
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j <= i; ++j)
      covariance(i, j) = (i == j ? 1. : 0.) + (i + 1.) * (j + 1.) * 0.01;
  auto reversed = TrackPropagation::reverseMomentumCovariance(covariance);
  auto restored = TrackPropagation::reverseMomentumCovariance(reversed);
  // Derive the signs independently from the transverse bases for p and -p.
  GlobalVector p(3., 4., 12.);
  auto x = GlobalVector(-p.y(), p.x(), 0.).unit();
  auto y = p.unit().cross(x);
  auto reverseX = GlobalVector(p.y(), -p.x(), 0.).unit();
  auto reverseY = (-p).unit().cross(reverseX);
  double jacobian[5] = {1., -1., 1., x.dot(reverseX), y.dot(reverseY)};
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j <= i; ++j) {
      check(std::abs(reversed(i, j) - jacobian[i] * jacobian[j] * covariance(i, j)) < 1.e-6);
      check(std::abs(restored(i, j) - covariance(i, j)) < 1.e-12);
    }
  // A positive spatial step in the reversed frame must map to a negative
  // physical path. A deflated (negative) step here would double-reverse it.
  double length = 100., cosLambda = p.perp() / p.mag();
  check(std::abs(jacobian[3] * length * cosLambda * jacobian[2] + length * cosLambda) < 1.e-4);
  check(std::abs(jacobian[4] * length * jacobian[1] + length) < 1.e-4);
}
