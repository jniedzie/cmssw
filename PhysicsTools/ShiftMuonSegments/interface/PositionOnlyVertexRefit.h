// Code style: human-code.
#ifndef PhysicsTools_ShiftMuonSegments_PositionOnlyVertexRefit_h
#define PhysicsTools_ShiftMuonSegments_PositionOnlyVertexRefit_h

#include "PhysicsTools/ShiftMuonSegments/interface/CommonVertexRefit.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <limits>

namespace shift {
  struct PositionOnlyVertexOptions {
    unsigned int maximumIterations = 16;
    unsigned int maximumPlaneEvaluations = 64;
    unsigned int maximumBacktracks = 8;
    double beamLineSigma = 0.;  // Zero explicitly disables the transverse prior.
    double absoluteToleranceCm = .1;
    double normalizedTolerance = .02;
    double maximumZStepCm = 2000.;  // Numerical trust step, not a vertex window.
  };

  struct PositionOnlyVertexFit {
    enum Status { success = 1, invalidInput = -1, invalidTransport = -2,
                  invalidCovariance = -3, singularGeometry = -4,
                  noDescent = -5, iterationLimit = -6, transportLimit = -7 };
    bool valid = false;
    int status = invalidInput;
    Eigen::Vector3d position = Eigen::Vector3d::Zero();
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    // These are the ORIGINAL transported states at position.z(). No q/p,
    // slope, or covariance update is applied to either track.
    std::array<VertexTrack, 2> tracks;
    double chi2 = -1.;
    double maximumUpdateCm = -1.;
    double normalizedUpdate = -1.;
    double jointUpdate = -1.;
    unsigned int iterations = 0;
    unsigned int planeEvaluations = 0;
  };

  inline bool positionOnlyWeights(std::array<VertexTrack, 2> const& tracks,
                                  std::array<Eigen::Matrix2d, 2>& weights) {
    for (unsigned int i = 0; i < 2; ++i) {
      auto const& covariance = tracks[i].covariance;
      if (!tracks[i].parameters.allFinite() || !covariance.allFinite() ||
          (covariance - covariance.transpose()).norm() > 1.e-10 * std::max(covariance.norm(), 1.e-30) ||
          covariance.llt().info() != Eigen::Success)
        return false;
      Eigen::LLT<Eigen::Matrix2d> solve(covariance.bottomRightCorner<2, 2>());
      if (solve.info() != Eigen::Success) return false;
      weights[i] = solve.solve(Eigen::Matrix2d::Identity());
    }
    return true;
  }

  // Profile x/y at one z plane using FIXED weights. This is also the actual
  // transported-geometry objective used during the nonlinear line search.
  inline double positionOnlyObjective(std::array<VertexTrack, 2> const& tracks,
                                      std::array<Eigen::Matrix2d, 2> const& weights,
                                      double priorWeight, Eigen::Vector2d& vertex) {
    Eigen::Matrix2d normal = weights[0] + weights[1];
    normal.diagonal().array() += priorWeight;
    Eigen::Vector2d rhs = Eigen::Vector2d::Zero();
    for (unsigned int i = 0; i < 2; ++i) {
      if (!tracks[i].parameters.allFinite() || !weights[i].allFinite())
        return std::numeric_limits<double>::quiet_NaN();
      rhs += weights[i] * tracks[i].parameters.tail<2>();
    }
    Eigen::LLT<Eigen::Matrix2d> solve(normal);
    if (solve.info() != Eigen::Success) return std::numeric_limits<double>::quiet_NaN();
    vertex = solve.solve(rhs);
    double chi2 = priorWeight * vertex.squaredNorm();
    for (unsigned int i = 0; i < 2; ++i) {
      Eigen::Vector2d const residual = vertex - tracks[i].parameters.tail<2>();
      chi2 += residual.dot(weights[i] * residual);
    }
    return vertex.allFinite() && std::isfinite(chi2) ? chi2 : std::numeric_limits<double>::quiet_NaN();
  }

  // Callback transport(z, tracks) must start afresh from the two original,
  // independent detector posteriors, returning full states on the SAME z
  // plane. It must never transport a previously vertex-updated state.
  //
  // This is an iterated position-moment fit, not a likelihood containing
  // log(det covariance). Covariance is frozen throughout each line search;
  // changing process noise cannot make a rejected trial cheaper. A local
  // covariance is reported, but its coverage must be validated separately.
  template <class Transport>
  PositionOnlyVertexFit fitPositionOnlyVertex(double seedZ, Transport const& transport,
                                              PositionOnlyVertexOptions const& options = {}) {
    PositionOnlyVertexFit result;
    if (!std::isfinite(seedZ) || !(options.beamLineSigma >= 0.) || !std::isfinite(options.beamLineSigma) ||
        !(options.absoluteToleranceCm > 0.) || !std::isfinite(options.absoluteToleranceCm) ||
        !(options.normalizedTolerance > 0.) || !std::isfinite(options.normalizedTolerance) ||
        !(options.maximumZStepCm > 0.) || !std::isfinite(options.maximumZStepCm) ||
        !options.maximumIterations || !options.maximumPlaneEvaluations || !options.maximumBacktracks)
      return result;
    double const priorWeight = options.beamLineSigma > 0. ? 1. / (options.beamLineSigma * options.beamLineSigma) : 0.;
    if (!std::isfinite(priorWeight)) return result;
    auto evaluate = [&](double z, std::array<VertexTrack, 2>& tracks) {
      if (result.planeEvaluations >= options.maximumPlaneEvaluations) {
        result.status = PositionOnlyVertexFit::transportLimit;
        return false;
      }
      ++result.planeEvaluations;
      return transport(z, tracks);
    };
    double z = seedZ;
    for (unsigned int iteration = 0; iteration < options.maximumIterations; ++iteration) {
      result.iterations = iteration + 1;
      std::array<VertexTrack, 2> tracks;
      if (!evaluate(z, tracks)) {
        if (result.status != PositionOnlyVertexFit::transportLimit) result.status = PositionOnlyVertexFit::invalidTransport;
        return result;
      }
      std::array<Eigen::Matrix2d, 2> weights;
      if (!positionOnlyWeights(tracks, weights)) {
        result.status = PositionOnlyVertexFit::invalidCovariance;
        return result;
      }
      Eigen::Vector2d transverse;
      double const chi2 = positionOnlyObjective(tracks, weights, priorWeight, transverse);
      Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
      normal(0, 0) = normal(1, 1) = priorWeight;
      Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
      for (unsigned int i = 0; i < 2; ++i) {
        Eigen::Matrix<double, 2, 3> design;
        design << 1., 0., -tracks[i].parameters[1], 0., 1., -tracks[i].parameters[2];
        normal += design.transpose() * weights[i] * design;
        rhs += design.transpose() * weights[i] * tracks[i].parameters.tail<2>();
      }
      Eigen::Vector3d const scales = normal.diagonal().array().sqrt().inverse();
      result.status = PositionOnlyVertexFit::singularGeometry;
      if (!scales.allFinite() || !rhs.allFinite() || !std::isfinite(chi2)) return result;
      Eigen::Matrix3d const scaledNormal = scales.asDiagonal() * normal * scales.asDiagonal();
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenvalues(scaledNormal, Eigen::EigenvaluesOnly);
      Eigen::LLT<Eigen::Matrix3d> solve(scaledNormal);
      if (eigenvalues.info() != Eigen::Success || solve.info() != Eigen::Success ||
          eigenvalues.eigenvalues().minCoeff() <= 1.e-12 * eigenvalues.eigenvalues().maxCoeff())
        return result;
      Eigen::Matrix3d const covariance = scales.asDiagonal() * solve.solve(Eigen::Matrix3d::Identity()) * scales.asDiagonal();
      Eigen::Vector3d const target = scales.asDiagonal() * solve.solve((scales.array() * rhs.array()).matrix());
      Eigen::Vector3d current(transverse[0], transverse[1], 0.);
      Eigen::Vector3d const delta = target - current;
      double const normalized = (delta.array().abs() / covariance.diagonal().array().sqrt()).maxCoeff();
      double const reduction = delta.dot(normal * delta);
      if (!covariance.allFinite() || !delta.allFinite() || !std::isfinite(normalized) ||
          !(reduction >= 0.) || !std::isfinite(reduction)) return result;
      result.position << transverse[0], transverse[1], z;
      result.covariance = covariance;
      result.tracks = tracks;
      result.chi2 = chi2;
      result.maximumUpdateCm = delta.cwiseAbs().maxCoeff();
      result.normalizedUpdate = normalized;
      result.jointUpdate = std::sqrt(reduction);
      if (result.maximumUpdateCm < options.absoluteToleranceCm && normalized < options.normalizedTolerance &&
          result.jointUpdate < options.normalizedTolerance * std::sqrt(3.)) {
        result.valid = true;
        result.status = PositionOnlyVertexFit::success;
        return result;
      }
      double fraction = std::abs(delta[2]) > options.maximumZStepCm ? options.maximumZStepCm / std::abs(delta[2]) : 1.;
      bool accepted = false;
      for (unsigned int backtrack = 0; backtrack < options.maximumBacktracks; ++backtrack, fraction *= .5) {
        double const trialZ = z + fraction * delta[2];
        if (!std::isfinite(trialZ) || trialZ == z) continue;
        std::array<VertexTrack, 2> trial;
        if (!evaluate(trialZ, trial)) {
          if (result.status == PositionOnlyVertexFit::transportLimit) return result;
          continue;
        }
        Eigen::Vector2d trialTransverse;
        double const trialChi2 = positionOnlyObjective(trial, weights, priorWeight, trialTransverse);
        if (std::isfinite(trialChi2) && trialChi2 < chi2 - 2.e-4 * fraction * reduction) {
          z = trialZ;
          accepted = true;
          break;
        }
      }
      if (!accepted) {
        result.status = PositionOnlyVertexFit::noDescent;
        return result;
      }
    }
    result.status = PositionOnlyVertexFit::iterationLimit;
    return result;
  }
}
#endif
