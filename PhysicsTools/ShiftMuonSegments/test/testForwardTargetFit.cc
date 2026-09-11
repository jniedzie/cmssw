#include "PhysicsTools/ShiftMuonSegments/interface/ForwardTargetFit.h"
#include <Eigen/LU>
#include <stdexcept>

int main() {
  using namespace shift;
  auto check = [](bool value) { if (!value) throw std::runtime_error("Forward target fit regression"); };
  TargetCovariance jacobian = TargetCovariance::Identity();
  jacobian(3,1) = -14000.; jacobian(4,2) = -14000.;
  TargetCovariance detectorCovariance = TargetCovariance::Identity();
  detectorCovariance(0,0) = 1.e-6;
  detectorCovariance(1,1) = detectorCovariance(2,2) = 1.e-8;
  TargetCovariance noise = detectorCovariance*.3;
  TargetParameters detector;
  detector << .01, .01, -.02, -145., 282.;
  auto transport = [&](TargetParameters const& state) {
    return ForwardTargetPrediction{true, jacobian*state, noise};
  };
  TargetParameters seed = jacobian.inverse()*detector;
  auto fitted = fitForwardTarget(seed, detector, detectorCovariance, .1, .1, 0., transport);
  check(fitted.valid);
  // Independent closed-form GLS solution with each observation counted once.
  TargetCovariance weight = (detectorCovariance + noise).inverse();
  TargetCovariance normal = jacobian.transpose()*weight*jacobian;
  normal(3,3) += 100.; normal(4,4) += 100.;
  TargetCovariance expectedCovariance = normal.inverse();
  TargetParameters expected = expectedCovariance*jacobian.transpose()*weight*detector;
  check((fitted.parameters-expected).norm() < 1.e-6);
  check((fitted.covariance-expectedCovariance).norm() < 1.e-8);
  // A deliberately different seed must not become another momentum prior.
  seed[0] *= 2.; seed[1] += .005;
  auto repeated = fitForwardTarget(seed, detector, detectorCovariance, .1, .1, 0., transport);
  check(repeated.valid && (repeated.parameters-fitted.parameters).norm() < 1.e-6);
  check(!fitForwardTarget(seed, detector, detectorCovariance, .1, .1, 0.,
      [](TargetParameters const&) { return ForwardTargetPrediction{}; }).valid);
  // Independent scalar likelihood minimum when process variance depends on
  // momentum. A fit that ignores the determinant would incorrectly return y.
  TargetParameters observation = TargetParameters::Zero();
  observation[0] = .02;
  TargetCovariance measurement = TargetCovariance::Identity();
  measurement(0,0) = 1.e-6;
  auto varyingNoise = [](TargetParameters const& state) {
    ForwardTargetPrediction prediction{true, state, TargetCovariance::Zero()};
    prediction.noise(0,0) = .2*state[0]*state[0];
    return prediction;
  };
  auto varying = fitForwardTarget(observation, observation, measurement, .1, .1, 0., varyingNoise, 64);
  check(varying.valid);
  auto objective = [](double x) {
    double const variance = 1.e-6 + .2*x*x;
    return (x-.02)*(x-.02)/variance + std::log(variance);
  };
  double lower = .005, upper = .04;
  for (unsigned int i = 0; i < 100; ++i) {
    double const a = (2.*lower+upper)/3., b = (lower+2.*upper)/3.;
    if (objective(a) < objective(b)) upper = b; else lower = a;
  }
  check(std::abs(varying.parameters[0] - (lower+upper)/2.) < .02*std::sqrt(varying.covariance(0,0)));
  check(std::abs(varying.parameters[0] - observation[0]) > .0005);
  for (unsigned int i = 1; i < varying.steps.size(); ++i)
    check(varying.steps[i].objective < varying.steps[i-1].objective);
}
