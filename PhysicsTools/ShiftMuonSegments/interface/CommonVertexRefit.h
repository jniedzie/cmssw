#ifndef PhysicsTools_ShiftMuonSegments_CommonVertexRefit_h
#define PhysicsTools_ShiftMuonSegments_CommonVertexRefit_h

#include <Eigen/Cholesky>
#include <array>
#include <cmath>

namespace shift {
  using VertexTrackParameters = Eigen::Matrix<double, 5, 1>;
  using VertexTrackCovariance = Eigen::Matrix<double, 5, 5>;
  struct VertexTrack {
    // Local (q/p, dx/dz, dy/dz, x, y) on a global-z plane.
    VertexTrackParameters parameters;
    VertexTrackCovariance covariance;
  };
  struct CommonVertexRefit {
    bool valid = false;
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero();
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    std::array<VertexTrack, 2> tracks;
    VertexTrackCovariance crossCovariance = VertexTrackCovariance::Zero();
    double chi2 = 0.;
  };

  // Joint linearized fit: two independent track priors, a common vertex and
  // the existing transverse beam-line prior. There is no target-z prior.
  // Eliminate each track's parameters, solve the 3D vertex, then smooth the
  // tracks. The vertex is NOT reused as an independent measurement of its
  // input tracks. Its uncertainty and the induced inter-track covariance are
  // retained explicitly. Caller relinearizes from original detector priors.
  inline CommonVertexRefit refitCommonVertex(std::array<VertexTrack, 2> const& inputs, double beamLineSigma) {
    CommonVertexRefit result;
    if (!(beamLineSigma > 0.) || !std::isfinite(beamLineSigma))
      return result;
    std::array<Eigen::Matrix<double, 2, 3>, 2> design;
    std::array<Eigen::Matrix2d, 2> weight;
    std::array<Eigen::Matrix<double, 5, 2>, 2> gain;
    Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
    normal(0, 0) = normal(1, 1) = 1. / (beamLineSigma * beamLineSigma);
    Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
    for (unsigned int i = 0; i < 2; ++i) {
      auto const& track = inputs[i];
      if (!track.parameters.allFinite() || !track.covariance.allFinite() ||
          track.covariance.llt().info() != Eigen::Success)
        return result;
      design[i] << 1., 0., -track.parameters[1], 0., 1., -track.parameters[2];
      weight[i] = track.covariance.bottomRightCorner<2, 2>().llt().solve(Eigen::Matrix2d::Identity());
      gain[i] = track.covariance.rightCols<2>() * weight[i];
      normal += design[i].transpose() * weight[i] * design[i];
      rhs += design[i].transpose() * weight[i] * track.parameters.tail<2>();
    }
    Eigen::LLT<Eigen::Matrix3d> solve(normal);
    if (solve.info() != Eigen::Success)
      return result;
    result.displacement = solve.solve(rhs);
    result.covariance = solve.solve(Eigen::Matrix3d::Identity());
    if (!result.displacement.allFinite() || !result.covariance.allFinite())
      return result;
    result.chi2 = result.displacement.head<2>().squaredNorm() / (beamLineSigma * beamLineSigma);
    for (unsigned int i = 0; i < 2; ++i) {
      auto const residual = (design[i] * result.displacement - inputs[i].parameters.tail<2>()).eval();
      result.tracks[i].parameters = inputs[i].parameters + gain[i] * residual;
      result.tracks[i].covariance =
          inputs[i].covariance - gain[i] * inputs[i].covariance.bottomRows<2>() +
          gain[i] * design[i] * result.covariance * design[i].transpose() * gain[i].transpose();
      result.chi2 += residual.dot(weight[i] * residual);
      if (!result.tracks[i].parameters.allFinite() || !result.tracks[i].covariance.allFinite() ||
          result.tracks[i].parameters[0] * inputs[i].parameters[0] <= 0.)
        return result;
    }
    result.crossCovariance = gain[0] * design[0] * result.covariance * design[1].transpose() * gain[1].transpose();
    result.valid = std::isfinite(result.chi2) && result.chi2 >= 0.;
    return result;
  }
}  // namespace shift
#endif
