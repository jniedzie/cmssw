#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "PhysicsTools/ShiftMuonSegments/interface/SeedCovariance.h"

using Parameters = std::array<double, 5>;
using Vector = std::array<double, 3>;
double dot(Vector const& a, Vector const& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
std::array<Vector, 2> basis(double lambda, double phi) {
  return {Vector{-std::sin(phi), std::cos(phi), 0.},
          Vector{-std::sin(lambda) * std::cos(phi), -std::sin(lambda) * std::sin(phi), std::cos(lambda)}};
}
Parameters physicalTransform(Parameters const& a, double charge, double sign, double scale) {
  double const p = charge / a[0];
  Vector momentum{sign * scale * p * std::cos(a[1]) * std::cos(a[2]),
                  sign * scale * p * std::cos(a[1]) * std::sin(a[2]),
                  sign * scale * p * std::sin(a[1])};
  double const norm = std::sqrt(dot(momentum, momentum));
  double const lambda = std::atan2(momentum[2], std::hypot(momentum[0], momentum[1]));
  double const phi = std::atan2(momentum[1], momentum[0]);
  auto oldBasis = basis(a[1], a[2]), newBasis = basis(lambda, phi);
  Vector position;
  for (unsigned int k = 0; k < 3; ++k)
    position[k] = a[3] * oldBasis[0][k] + a[4] * oldBasis[1][k];
  return {charge * sign / norm, lambda, phi, dot(position, newBasis[0]), dot(position, newBasis[1])};
}
int main() {
  AlgebraicSymMatrix55 covariance;
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j <= i; ++j)
      covariance(i, j) = i == j ? 2. + i : 0.03 * (i + 1) * (j + 1);
  for (double charge : {-1., 1.})
    for (double sign : {-1., 1.})
      for (double scale : {0.5, 1., 2.}) {
        Parameters const nominal{charge / 50., 0.6, 2.4, 0.3, -0.2};
        double jacobian[5][5];
        for (unsigned int column = 0; column < 5; ++column) {
          auto plus = nominal, minus = nominal;
          double const step = 1.e-6;
          plus[column] += step;
          minus[column] -= step;
          auto high = physicalTransform(plus, charge, sign, scale), low = physicalTransform(minus, charge, sign, scale);
          for (unsigned int row = 0; row < 5; ++row) {
            double difference = high[row] - low[row];
            if (row == 2)
              difference = std::remainder(difference, 2. * std::acos(-1.));
            jacobian[row][column] = difference / (2. * step);
          }
        }
        auto actual = shift::seedCovariance(covariance, sign, scale);
        for (unsigned int row = 0; row < 5; ++row)
          for (unsigned int column = 0; column <= row; ++column) {
            double expected = 0.;
            for (unsigned int i = 0; i < 5; ++i)
              for (unsigned int j = 0; j < 5; ++j)
                expected += jacobian[row][i] * covariance(i, j) * jacobian[column][j];
            if (std::abs(actual(row, column) - expected) > 1.e-7)
              throw std::runtime_error("Seed covariance disagrees with Cartesian finite-difference transport");
          }
      }
  std::cout << "Seed covariance matches Cartesian momentum/charge reversal and rescaling\n";
}
