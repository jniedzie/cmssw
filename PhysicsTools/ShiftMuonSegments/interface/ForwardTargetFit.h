#ifndef PhysicsTools_ShiftMuonSegments_ForwardTargetFit_h
#define PhysicsTools_ShiftMuonSegments_ForwardTargetFit_h

#include <Eigen/Cholesky>
#include <array>
#include <cmath>
#include <vector>

namespace shift {
  using TargetParameters = Eigen::Matrix<double, 5, 1>;
  using TargetCovariance = Eigen::Matrix<double, 5, 5>;
  struct ForwardTargetPrediction {
    bool valid = false;
    TargetParameters parameters;
    TargetCovariance noise;
  };
  struct ForwardTargetFit {
    bool valid = false;
    int status = -1;
    TargetParameters parameters;
    TargetCovariance covariance;
    double chi2 = -1.;
    unsigned int iterations = 0;
    struct Step { TargetParameters parameters; double chi2; double maxUpdate; double objective; };
    std::vector<Step> steps;
  };

  // Fit one detector posterior and one target-position measurement. The
  // transported detector posterior is used only to initialize the search;
  // it is NOT reused as an additional prior on target momentum. The callback
  // transports zero input covariance from target to detector, providing the
  // mean and process noise for the current trajectory through the geometry.
  template <typename Transport>
  ForwardTargetFit fitForwardTarget(TargetParameters const& seed,
                                    TargetParameters const& detector,
                                    TargetCovariance const& detectorCovariance,
                                    double sigmaX, double sigmaY, double sigmaZ,
                                    Transport const& transport,
                                    unsigned int maximumIterations = 8) {
    ForwardTargetFit result;
    if (!(sigmaX > 0.) || !(sigmaY > 0.) || !(sigmaZ >= 0.) ||
        !std::isfinite(sigmaX) || !std::isfinite(sigmaY) || !std::isfinite(sigmaZ) || !seed.allFinite() ||
        !detector.allFinite() || !detectorCovariance.allFinite() || detectorCovariance.llt().info() != Eigen::Success)
      return result;
    TargetParameters state = seed;
    auto hitCovariance = [&](TargetParameters const& parameters) {
      Eigen::Matrix2d hit;
      hit << sigmaX*sigmaX + parameters[1]*parameters[1]*sigmaZ*sigmaZ, parameters[1]*parameters[2]*sigmaZ*sigmaZ,
             parameters[1]*parameters[2]*sigmaZ*sigmaZ, sigmaY*sigmaY + parameters[2]*parameters[2]*sigmaZ*sigmaZ;
      return hit;
    };
    for (unsigned int iteration = 0; iteration < maximumIterations; ++iteration) {
      result.iterations = iteration + 1;
      result.status = -2;
      auto const prediction = transport(state);
      if (!prediction.valid) return result;
      result.status = -3;
      TargetCovariance const measurementCovariance = detectorCovariance + prediction.noise;
      Eigen::LLT<TargetCovariance> weightSolve(measurementCovariance);
      if (weightSolve.info() != Eigen::Success) return result;
      TargetCovariance const weight = weightSolve.solve(TargetCovariance::Identity());
      Eigen::LLT<Eigen::Matrix2d> hitSolve(hitCovariance(state));
      if (hitSolve.info() != Eigen::Success) return result;
      Eigen::Matrix2d const hitWeight = hitSolve.solve(Eigen::Matrix2d::Identity());
      TargetCovariance jacobian;
      std::array<TargetCovariance, 5> noiseDerivatives;
      std::array<Eigen::Matrix2d, 5> hitDerivatives;
      result.status = -4;
      for (unsigned int j = 0; j < 5; ++j) {
        double const h = j == 0 ? std::abs(state[0])*1.e-3 : (j < 3 ? 1.e-5 : .01);
        if (!(h > 0.)) return result;
        auto plus = state, minus = state;
        plus[j] += h; minus[j] -= h;
        auto const a = transport(plus), b = transport(minus);
        if (!a.valid || !b.valid) return result;
        jacobian.col(j) = (a.parameters - b.parameters) / (2.*h);
        noiseDerivatives[j] = (a.noise - b.noise) / (2.*h);
        hitDerivatives[j] = (hitCovariance(plus) - hitCovariance(minus)) / (2.*h);
      }
      TargetParameters const residual = detector - prediction.parameters;
      result.status = -5;
      TargetCovariance normal = jacobian.transpose()*weight*jacobian;
      normal.template bottomRightCorner<2, 2>() += hitWeight;
      TargetParameters rhs = jacobian.transpose()*weight*residual;
      rhs.template tail<2>() -= hitWeight*state.template tail<2>();
      // The material path and scattering noise depend on the fitted state.
      // Use the full Gaussian likelihood, including log determinants. Merely
      // reweighting least squares can improve its old weights while making
      // the new trajectory likelihood dramatically worse at a boundary.
      // The covariance terms below are its score and Fisher information.
      TargetParameters const weightedResidual = weight*residual;
      Eigen::Vector2d const weightedHit = hitWeight*state.template tail<2>();
      for (unsigned int i = 0; i < 5; ++i) {
        rhs[i] += .5*(weightedResidual.dot(noiseDerivatives[i]*weightedResidual) -
                      (weight*noiseDerivatives[i]).trace() +
                      weightedHit.dot(hitDerivatives[i]*weightedHit) - (hitWeight*hitDerivatives[i]).trace());
        for (unsigned int j = 0; j < 5; ++j)
          normal(i,j) += .5*((weight*noiseDerivatives[i]*weight*noiseDerivatives[j]).trace() +
                             (hitWeight*hitDerivatives[i]*hitWeight*hitDerivatives[j]).trace());
      }
      TargetParameters const scales = normal.diagonal().array().sqrt().inverse();
      if (!scales.allFinite()) return result;
      TargetCovariance const scaledNormal = scales.asDiagonal()*normal*scales.asDiagonal();
      Eigen::LLT<TargetCovariance> solve(scaledNormal);
      if (solve.info() != Eigen::Success) return result;
      TargetCovariance const covariance =
          scales.asDiagonal()*solve.solve(TargetCovariance::Identity())*scales.asDiagonal();
      TargetParameters const delta = covariance*rhs;
      double const chi2 = residual.dot(weight*residual) +
          state.template tail<2>().dot(hitWeight*state.template tail<2>());
      double const objective = chi2 +
          2.*weightSolve.matrixL().toDenseMatrix().diagonal().array().log().sum() +
          2.*hitSolve.matrixL().toDenseMatrix().diagonal().array().log().sum();
      double const change = (delta.array().abs()/covariance.diagonal().array().sqrt()).maxCoeff();
      if (!delta.allFinite() || !std::isfinite(chi2) || !std::isfinite(change) || !std::isfinite(objective)) return result;
      result.steps.push_back({state, chi2, change, objective});
      result.iterations = iteration + 1;
      if (change < .02) {
        result.valid = true;
        result.status = 1;
        result.parameters = state;
        result.covariance = covariance;
        result.chi2 = chi2;
        return result;
      }
      bool accepted = false;
      result.status = -6;
      for (double fraction = 1.; fraction >= 1./64.; fraction *= .5) {
        TargetParameters const trial = state + fraction*delta;
        if (trial[0]*seed[0] <= 0.) continue;
        auto const next = transport(trial);
        if (!next.valid) continue;
        TargetParameters const nextResidual = detector - next.parameters;
        Eigen::LLT<TargetCovariance> nextSolve(detectorCovariance + next.noise);
        Eigen::LLT<Eigen::Matrix2d> nextHitSolve(hitCovariance(trial));
        if (nextSolve.info() != Eigen::Success || nextHitSolve.info() != Eigen::Success) continue;
        double const nextObjective = nextResidual.dot(nextSolve.solve(nextResidual)) +
            trial.template tail<2>().dot(nextHitSolve.solve(trial.template tail<2>())) +
            2.*nextSolve.matrixL().toDenseMatrix().diagonal().array().log().sum() +
            2.*nextHitSolve.matrixL().toDenseMatrix().diagonal().array().log().sum();
        if (std::isfinite(nextObjective) && nextObjective < objective) {
          state = trial;
          accepted = true;
          break;
        }
      }
      if (!accepted) return result;
    }
    result.status = -7;
    return result;
  }
}  // namespace shift
#endif
