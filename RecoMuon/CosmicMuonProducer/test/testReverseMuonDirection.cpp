// Code style: human-code.
#include "RecoMuon/CosmicMuonProducer/interface/ReverseMuonDirection.h"
#include "MagneticField/UniformEngine/interface/UniformMagneticField.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryError.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace {
  void requireAt(bool condition, int line) {
    if (!condition) {
      std::cerr << "Requirement failed at line " << line << '\n';
      std::abort();
    }
  }
}
#define require(condition) requireAt(condition, __LINE__)

int main() {
  CosmicMuonUtilities utilities;
  unsigned int cases = 0;
  double largestCorrectedError = 0.;
  double largestInvolutionError = 0.;
  for (float fieldValue : {0.f, 3.8f}) {
    UniformMagneticField field(fieldValue);
    auto surface = Plane::build(GlobalPoint(10., 20., 1030.), Surface::RotationType());
    for (float charge : {-1.f, 1.f}) {
      for (float sign : {-1.f, 1.f}) {
        AlgebraicSymMatrix55 covariance;
        double const scale[5] = {1.e-4, 3.e-4, 4.e-4, .1, .2};
        for (unsigned int row = 0; row < 5; ++row)
          for (unsigned int column = 0; column <= row; ++column)
            covariance(row, column) = row == column ? scale[row] * scale[row] : 0.;
        covariance(3, 0) = .1 * scale[3] * scale[0];
        covariance(4, 1) = -.1 * scale[4] * scale[1];
        LocalTrajectoryParameters const parameters(charge / 200.f, .01f, -.02f, .3f, -.4f, sign);
        TrajectoryStateOnSurface const original(parameters, LocalTrajectoryError(covariance), *surface, &field);
        auto native = original;
        utilities.reverseDirection(native, &field);
        auto off = original;
        cosmic_muon::reverseMuonDirection(off, &field, utilities, false);
        require((off.globalPosition() - native.globalPosition()).mag() == 0.);
        require((off.globalMomentum() - native.globalMomentum()).mag() == 0.);
        require(off.charge() == native.charge());
        auto corrected = original;
        cosmic_muon::reverseMuonDirection(corrected, &field, utilities, true);
        require(corrected.isValid() && corrected.hasError() && corrected.localError().posDef());
        require((corrected.globalPosition() - original.globalPosition()).mag() == 0.);
        require((corrected.globalMomentum() + original.globalMomentum()).mag() == 0.);
        require(corrected.charge() == -original.charge());
        auto const input = original.cartesianError().matrix();
        auto const output = corrected.cartesianError().matrix();
        Eigen::Matrix<double, 6, 6> dense;
        double largestEntry = 1.;
        double correctedError = 0.;
        double stockError = 0.;
        AlgebraicSymMatrix66 expected;
        double const diagonal[6] = {1., 1., 1., -1., -1., -1.};
        for (unsigned int row = 0; row < 6; ++row)
          for (unsigned int column = 0; column <= row; ++column)
            expected(row, column) = diagonal[row] * input(row, column) * diagonal[column];
        GlobalTrajectoryParameters const reversedMean(
            original.globalPosition(), -original.globalMomentum(), -original.charge(), &field);
        TrajectoryStateOnSurface const reference(reversedMean, CartesianTrajectoryError(expected), *surface);
        for (unsigned int row = 0; row < 6; ++row) {
          for (unsigned int column = 0; column < 6; ++column)
            dense(row, column) = output(row, column);
          for (unsigned int column = 0; column <= row; ++column) {
            require(off.cartesianError().matrix()(row, column) == native.cartesianError().matrix()(row, column));
            require(output(row, column) == reference.cartesianError().matrix()(row, column));
            correctedError = std::max(correctedError, std::abs(output(row, column) - expected(row, column)));
            stockError = std::max(stockError, std::abs(native.cartesianError().matrix()(row, column) - expected(row, column)));
            largestEntry = std::max(largestEntry, std::abs(expected(row, column)));
          }
        }
        // Native local/global chart conversions use float momentum components.
        double const projectionBound = 4. * std::numeric_limits<float>::epsilon() * largestEntry;
        require(correctedError <= projectionBound);
        require(stockError > projectionBound);
        largestCorrectedError = std::max(largestCorrectedError, correctedError);
        // The fixed-surface five-parameter Cartesian covariance has one null direction.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> eigenvalues(dense);
        require(eigenvalues.info() == Eigen::Success);
        require(eigenvalues.eigenvalues().minCoeff() >= -1.e-10 * largestEntry);
        cosmic_muon::reverseMuonDirection(corrected, &field, utilities, true);
        require((corrected.globalPosition() - original.globalPosition()).mag() == 0.);
        require((corrected.globalMomentum() - original.globalMomentum()).mag() == 0.);
        require(corrected.charge() == original.charge());
        GlobalTrajectoryParameters const originalMean(
            original.globalPosition(), original.globalMomentum(), original.charge(), &field);
        TrajectoryStateOnSurface const roundtrip(originalMean, CartesianTrajectoryError(input), *surface);
        TrajectoryStateOnSurface const twiceRoundtrip(
            originalMean, roundtrip.cartesianError(), *surface);
        double involutionError = 0.;
        for (unsigned int row = 0; row < 6; ++row)
          for (unsigned int column = 0; column <= row; ++column) {
            require(corrected.cartesianError().matrix()(row, column) ==
                    twiceRoundtrip.cartesianError().matrix()(row, column));
            involutionError = std::max(
                involutionError, std::abs(corrected.cartesianError().matrix()(row, column) - input(row, column)));
          }
        require(involutionError <= 8. * std::numeric_limits<float>::epsilon() * largestEntry);
        largestInvolutionError = std::max(largestInvolutionError, involutionError);
        ++cases;
      }
    }
  }
  require(cases == 8);
  std::cout << "{\"cases\":" << cases << ",\"max_corrected_error\":" << largestCorrectedError
            << ",\"max_involution_error\":" << largestInvolutionError << "}\n";
}
