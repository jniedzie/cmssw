// Code style: human-code.
#include "PhysicsTools/ShiftMuonSegments/interface/PositionOnlyVertexRefit.h"
#include <iostream>
#include <stdexcept>

namespace {
  void require(bool value, char const* message) {
    if (!value) throw std::runtime_error(message);
  }

  std::array<shift::VertexTrack, 2> crossing(double z) {
    std::array<shift::VertexTrack, 2> tracks;
    for (unsigned int i = 0; i < 2; ++i) {
      double const tx = i ? -.015 : .02;
      double const ty = i ? .01 : -.005;
      tracks[i].parameters << (i ? -.002 : .001), tx, ty,
          3. + tx * (z - 7000.), -4. + ty * (z - 7000.);
      Eigen::Matrix<double, 5, 5> factor = Eigen::Matrix<double, 5, 5>::Identity();
      factor.diagonal() << .0002, .0003, .0004, (i ? 4. : 2.), (i ? 1. : 3.);
      factor(3, 0) = .5;  // Curvature/position correlation must not update q/p.
      factor(4, 3) = .4;
      tracks[i].covariance = factor * factor.transpose();
    }
    return tracks;
  }

  void straightGeometry() {
    for (double seed : {-3000., 6999., 17000.}) {
      auto fit = shift::fitPositionOnlyVertex(seed, [](double z, std::array<shift::VertexTrack, 2>& tracks) {
        tracks = crossing(z); return true;
      });
      require(fit.valid, "Straight common-plane fit failed");
      require((fit.position - Eigen::Vector3d(3., -4., 7000.)).norm() < 1.e-7,
              "Unconstrained displaced intersection was not recovered");
      auto const expected = crossing(fit.position[2]);
      for (unsigned int i = 0; i < 2; ++i) {
        require((fit.tracks[i].parameters - expected[i].parameters).norm() == 0.,
                "Position fit changed transported track parameters");
        require((fit.tracks[i].covariance - expected[i].covariance).norm() == 0.,
                "Position fit changed transported track covariance");
      }
      require(fit.covariance.llt().info() == Eigen::Success, "Vertex covariance is indefinite");
      require(fit.maximumUpdateCm < .1 && fit.normalizedUpdate < .02 && fit.jointUpdate < .02 * std::sqrt(3.),
              "Success does not meet separate geometric and normalized convergence");
    }
    std::cout << "Displaced intersection, distant seeds, covariance and unchanged momentum states: PASS\n";
  }

  void weightedReference() {
    auto callback = [](double z, std::array<shift::VertexTrack, 2>& tracks) {
      tracks = crossing(z);
      tracks[1].parameters[4] += 2.;
      return true;
    };
    auto fit = shift::fitPositionOnlyVertex(6500., callback);
    require(fit.valid, "Skew-line weighted fit failed");
    std::array<shift::VertexTrack, 2> atZero;
    callback(0., atZero);
    // Independent 4x3 weighted least-squares reference from all observations.
    Eigen::Matrix<double, 4, 3> design;
    Eigen::Matrix4d covariance = Eigen::Matrix4d::Zero();
    Eigen::Vector4d observations;
    for (unsigned int i = 0; i < 2; ++i) {
      design.block<2, 3>(2 * i, 0) << 1., 0., -atZero[i].parameters[1],
                                               0., 1., -atZero[i].parameters[2];
      covariance.block<2, 2>(2 * i, 2 * i) = atZero[i].covariance.bottomRightCorner<2, 2>();
      observations.segment<2>(2 * i) = atZero[i].parameters.tail<2>();
    }
    Eigen::Matrix4d const weight = covariance.llt().solve(Eigen::Matrix4d::Identity());
    Eigen::Matrix3d const normal = design.transpose() * weight * design;
    Eigen::Vector3d const expected = normal.llt().solve(design.transpose() * weight * observations);
    require((fit.position - expected).norm() < 1.e-6, "Position differs from independent GLS reference");
    require((fit.covariance - normal.llt().solve(Eigen::Matrix3d::Identity())).norm() < 1.e-7,
            "Position covariance differs from independent GLS reference");
    auto swapped = shift::fitPositionOnlyVertex(6500., [&](double z, std::array<shift::VertexTrack, 2>& tracks) {
      callback(z, tracks); std::swap(tracks[0], tracks[1]); return true;
    });
    require(swapped.valid && (swapped.position - fit.position).norm() < 1.e-8,
            "Track order changes fitted position");
    auto translated = shift::fitPositionOnlyVertex(6500., [&](double z, std::array<shift::VertexTrack, 2>& tracks) {
      callback(z, tracks);
      for (auto& track : tracks) track.parameters.tail<2>() += Eigen::Vector2d(150., -230.);
      return true;
    });
    require(translated.valid && (translated.position - fit.position - Eigen::Vector3d(150., -230., 0.)).norm() < 1.e-6,
            "Prior-free estimator is not translation invariant");
    std::cout << "Independent correlated GLS reference, track swap and transverse translation: PASS\n";
  }

  void curvedAndFrozenWeights() {
    unsigned int probes = 0;
    auto curved = [&](double z, std::array<shift::VertexTrack, 2>& tracks) {
      ++probes;
      tracks = crossing(z);
      double const dz = z - 7000.;
      for (unsigned int i = 0; i < 2; ++i) {
        double const bend = i ? -2.e-6 : 1.e-6;
        tracks[i].parameters[3] += bend * dz * dz;
        tracks[i].parameters[1] += 2. * bend * dz;
      }
      return true;
    };
    auto fit = shift::fitPositionOnlyVertex(8000., curved);
    require(fit.valid && (fit.position - Eigen::Vector3d(3., -4., 7000.)).norm() < .11,
            "Curved mean geometry did not close");
    require(fit.planeEvaluations == probes, "Plane-evaluation accounting differs from actual calls");
    auto original = crossing(7000.);
    std::array<Eigen::Matrix2d, 2> weights;
    require(shift::positionOnlyWeights(original, weights), "Valid weights rejected");
    auto trial = original; trial[1].parameters[3] += 100.;
    Eigen::Vector2d vertex;
    double const first = shift::positionOnlyObjective(trial, weights, 0., vertex);
    for (auto& track : trial) track.covariance *= 1.e12;
    double const inflated = shift::positionOnlyObjective(trial, weights, 0., vertex);
    require(first == inflated, "Trial can lower frozen objective by inflating covariance");
    std::cout << "Curved transport closure, exact call bounds and covariance-inflation immunity: PASS\n";
  }

  void failuresAndPrior() {
    auto parallel = [](double z, std::array<shift::VertexTrack, 2>& tracks) {
      tracks = crossing(z);
      tracks[1].parameters.segment<2>(1) = tracks[0].parameters.segment<2>(1);
      return true;
    };
    auto fit = shift::fitPositionOnlyVertex(7000., parallel);
    require(!fit.valid && fit.status == shift::PositionOnlyVertexFit::singularGeometry,
            "Unidentified parallel-track z accepted without a prior");
    fit = shift::fitPositionOnlyVertex(7000., [](double, auto&) { return false; });
    require(!fit.valid && fit.status == shift::PositionOnlyVertexFit::invalidTransport,
            "Failed transport accepted");
    fit = shift::fitPositionOnlyVertex(7000., [](double z, auto& tracks) {
      tracks = crossing(z); tracks[0].covariance(0, 0) = -1.; return true;
    });
    require(!fit.valid && fit.status == shift::PositionOnlyVertexFit::invalidCovariance,
            "Indefinite original covariance accepted");
    shift::PositionOnlyVertexOptions bounded; bounded.maximumPlaneEvaluations = 1;
    fit = shift::fitPositionOnlyVertex(-3000., [](double z, auto& tracks) { tracks = crossing(z); return true; }, bounded);
    require(!fit.valid && fit.status == shift::PositionOnlyVertexFit::transportLimit && fit.planeEvaluations == 1,
            "Transport budget exceeded or relabelled success");
    shift::PositionOnlyVertexOptions prior; prior.beamLineSigma = 100.;
    fit = shift::fitPositionOnlyVertex(7000., [](double z, auto& tracks) { tracks = crossing(z); return true; }, prior);
    require(fit.valid && fit.chi2 > 0., "Explicit optional beamline prior did not enter the objective");
    require(!shift::fitPositionOnlyVertex(std::numeric_limits<double>::quiet_NaN(), parallel).valid,
            "Nonfinite seed accepted");
    std::cout << "Explicit prior, singular geometry, covariance, transport and budget failures: PASS\n";
  }
}

int main() {
  straightGeometry();
  weightedReference();
  curvedAndFrozenWeights();
  failuresAndPrior();
}
