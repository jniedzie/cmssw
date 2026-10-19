// Code style: AI-code. AI-generated numerical regression of measured dimuon kinematics.
#include "PhysicsTools/ShiftMuonSegments/interface/DimuonInvariantMass.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <type_traits>

namespace {
  constexpr double muonMass = 0.105658;

  void check(bool passed, char const* message) {
    if (!passed)
      throw std::runtime_error(message);
  }

  double legacyMass(GlobalVector const& first, GlobalVector const& second) {
    double const x = double(first.x()) + second.x();
    double const y = double(first.y()) + second.y();
    double const z = double(first.z()) + second.z();
    double const e = std::sqrt(first.mag2() + muonMass * muonMass) +
                     std::sqrt(second.mag2() + muonMass * muonMass);
    return std::sqrt(std::max(0., e * e - x * x - y * y - z * z));
  }

  // Independent extended-precision four-vector subtraction, evaluated on
  // the same float-rounded measured inputs as the reconstruction helper.
  long double referenceMass(GlobalVector const& first, GlobalVector const& second) {
    long double energy = 0.;
    std::array<long double, 3> total{};
    for (auto const& momentum : {first, second}) {
      long double const x = momentum.x(), y = momentum.y(), z = momentum.z();
      long double const mass = muonMass;
      energy += std::sqrt(x * x + y * y + z * z + mass * mass);
      total[0] += x;
      total[1] += y;
      total[2] += z;
    }
    return std::sqrt(energy * energy - total[0] * total[0] - total[1] * total[1] - total[2] * total[2]);
  }
}

int main() {
  static_assert(std::is_same_v<decltype(GlobalVector().mag2()), float>);
  unsigned int assertions = 0;
  unsigned int legacyFailures = 0;
  for (float longitudinal : {0.f, 1.f, 100.f, 1000.f, 10000.f, 1000000.f, 1000000000.f}) {
    for (float transverse : {0.f, 0.01f, 0.1f, 1.f, 10.f, 100.f}) {
      GlobalVector const first(transverse, 0.f, longitudinal);
      GlobalVector const second(-transverse, 0.f, longitudinal);
      double const expected = 2. * std::hypot(double(transverse), muonMass);
      double const mass = shift::dimuonInvariantMass(first, second, muonMass);
      check(std::abs(mass - expected) < 2.e-12 * expected, "Boosted equal-momentum pair mass changed");
      check(shift::dimuonInvariantMass(second, first, muonMass) == mass, "Muon exchange changed invariant mass");
      legacyFailures += std::abs(legacyMass(first, second) - expected) > 0.01 * expected;
      assertions += 2;
    }
  }
  check(legacyFailures > 0, "Regression did not expose float-magnitude cancellation");
  ++assertions;

  for (float ratio : {0.01f, 0.1f, 0.7f, 1.f, 2.f, 10.f}) {
    GlobalVector const first(0.f, 0.f, 1000.f);
    GlobalVector const second(0.f, 0.f, 1000.f * ratio);
    double const mass = shift::dimuonInvariantMass(first, second, muonMass);
    long double const expected = referenceMass(first, second);
    check(std::abs(mass - expected) < 1.e-8 * expected, "Unequal collinear momenta failed");
    ++assertions;
  }

  std::mt19937 generator(19471019);
  std::uniform_real_distribution<float> component(-10000.f, 10000.f);
  for (unsigned int i = 0; i < 2000; ++i) {
    GlobalVector const first(component(generator), component(generator), component(generator));
    GlobalVector const second(component(generator), component(generator), component(generator));
    double const mass = shift::dimuonInvariantMass(first, second, muonMass);
    long double const expected = referenceMass(first, second);
    check(std::abs(mass - expected) < 1.e-10 * expected, "General pair differs from independent reference");
    check(mass >= 2. * muonMass, "Physical dimuon fell below two-muon threshold");
    assertions += 2;
  }
  for (GlobalVector const& moving : {GlobalVector(1.f, 2.f, 3.f), GlobalVector(0.f, 0.f, 10000.f)}) {
    GlobalVector const resting;
    double const momentum = std::hypot(double(moving.x()), double(moving.y()), double(moving.z()));
    double const expected = std::sqrt(2. * muonMass * (muonMass + std::hypot(momentum, muonMass)));
    check(std::abs(shift::dimuonInvariantMass(resting, moving, muonMass) - expected) < 1.e-12 * expected,
          "Resting-leg limit failed");
    ++assertions;
  }
  std::cout << assertions << " mass checks passed; " << legacyFailures
            << " boosted fixtures expose legacy errors above 1%.\n";
  GlobalVector const first(1.f, 0.f, 10000.f), second(-1.f, 0.f, 10000.f);
  std::cout << "Measured example [(1,0,10000),(-1,0,10000)] GeV: legacy="
            << legacyMass(first, second) << " stable=" << shift::dimuonInvariantMass(first, second, muonMass)
            << " GeV.\n";
}
