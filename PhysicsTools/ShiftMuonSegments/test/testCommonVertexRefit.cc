#include "PhysicsTools/ShiftMuonSegments/interface/CommonVertexRefit.h"
#include <Eigen/LU>
#include <iostream>
#include <stdexcept>

void require(bool value, char const* message) {
  if (!value)
    throw std::runtime_error(message);
}

int main() {
  // Independent full constrained least-squares solution: 10 track parameters,
  // 3 vertex coordinates and 4 Lagrange multipliers enforcing a common vertex.
  // Compare both means and ALL marginal/cross covariance blocks against the
  // reduced fitter. Correlated curvature-position errors exercise p updates.
  for (double sign : {-1., 1.}) {
    std::array<shift::VertexTrack, 2> tracks;
    for (unsigned int i = 0; i < 2; ++i) {
      tracks[i].parameters << sign * (i ? -.02 : .01), .02 * (i ? -1 : 1), .01, i ? -2. : 3., 1.;
      shift::VertexTrackCovariance factor = shift::VertexTrackCovariance::Identity();
      factor.diagonal() << .001, .002, .003, 2., 3.;
      factor(3, 0) = .5;
      factor(4, 1) = -.8;
      tracks[i].covariance = factor * factor.transpose();
    }
    double const beamSigma = 2.;
    auto const result = shift::refitCommonVertex(tracks, beamSigma);
    require(result.valid, "valid correlated input rejected");
    Eigen::MatrixXd system = Eigen::MatrixXd::Zero(17, 17);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(17);
    for (unsigned int i = 0; i < 2; ++i) {
      auto precision = tracks[i].covariance.inverse().eval();
      system.block<5, 5>(5 * i, 5 * i) = precision;
      rhs.segment<5>(5 * i) = precision * tracks[i].parameters;
      system(13 + 2 * i, 5 * i + 3) = 1.;
      system(14 + 2 * i, 5 * i + 4) = 1.;
      system(13 + 2 * i, 10) = -1.;
      system(14 + 2 * i, 11) = -1.;
      system(13 + 2 * i, 12) = tracks[i].parameters[1];
      system(14 + 2 * i, 12) = tracks[i].parameters[2];
    }
    system.block(0, 13, 13, 4) = system.block(13, 0, 4, 13).transpose().eval();
    system(10, 10) = system(11, 11) = 1. / (beamSigma * beamSigma);
    auto const inverse = system.inverse().eval();
    auto const fitted = (inverse * rhs).eval();
    require((fitted.segment<3>(10) - result.displacement).norm() < 1.e-8, "vertex mean mismatch");
    require((inverse.block<3, 3>(10, 10) - result.covariance).norm() < 1.e-7, "vertex covariance mismatch");
    for (unsigned int i = 0; i < 2; ++i) {
      require((fitted.segment<5>(5 * i) - result.tracks[i].parameters).norm() < 1.e-8, "track mean mismatch");
      require((inverse.block<5, 5>(5 * i, 5 * i) - result.tracks[i].covariance).norm() < 1.e-8,
              "track covariance mismatch");
    }
    require((inverse.block<5, 5>(0, 5) - result.crossCovariance).norm() < 1.e-8, "cross covariance mismatch");
    std::swap(tracks[0], tracks[1]);
    auto swapped = shift::refitCommonVertex(tracks, beamSigma);
    require(swapped.valid && (swapped.displacement - result.displacement).norm() < 1.e-9, "track-order dependence");
    require(!shift::refitCommonVertex(tracks, 0.).valid, "invalid prior accepted");
    tracks[0].covariance(0, 0) = -1.;
    require(!shift::refitCommonVertex(tracks, beamSigma).valid, "indefinite covariance accepted");
    tracks[0] = tracks[1];
    tracks[0].parameters[1] = tracks[1].parameters[1] = 0.;
    tracks[0].parameters[2] = tracks[1].parameters[2] = 0.;
    require(!shift::refitCommonVertex(tracks, beamSigma).valid, "unidentifiable vertex z accepted");
  }
  std::cout << "Joint vertex means, marginal/cross covariances and failure checks passed\n";
}
