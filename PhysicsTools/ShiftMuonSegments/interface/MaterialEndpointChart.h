#ifndef PhysicsTools_ShiftMuonSegments_MaterialEndpointChart_h
#define PhysicsTools_ShiftMuonSegments_MaterialEndpointChart_h

#include <Eigen/Core>
#include <cmath>

namespace shift {
  using MaterialEndpointMatrix = Eigen::Matrix<double, 5, 5>;

  // Endpoint quantities refer to physical momentum, NOT the propagation
  // stepping direction. State order is (q/p, tx, ty, x, y).
  struct MaterialEndpoint {
    double inverseMomentumFlow = 0.;  // d(q/p)/ds, in (1/GeV)/cm
    Eigen::Vector3d localPhysicalDirection = Eigen::Vector3d::UnitZ();
  };

  struct MaterialEndpointChartResult {
    enum Status { success = 1, invalidInput = -1, nonfiniteResult = -2 };
    bool valid = false;
    int status = invalidInput;
    MaterialEndpointMatrix inputChart = MaterialEndpointMatrix::Identity();
    MaterialEndpointMatrix outputChart = MaterialEndpointMatrix::Identity();
    MaterialEndpointMatrix jacobian = MaterialEndpointMatrix::Zero();
    MaterialEndpointMatrix noise = MaterialEndpointMatrix::Zero();
  };

  // Positive stopping-power convention:
  // a = q E/p^3 * (dE_loss/ds) = d(q/p)/ds along physical momentum.
  // This helper does not infer material, average stopping power over a step,
  // or change an existing free-propagation energy-loss derivative.
  inline double materialInverseMomentumFlow(double charge,
                                            double momentum,
                                            double mass,
                                            double positiveStoppingPower) {
    return charge * std::hypot(momentum, mass) / (momentum * momentum * momentum) * positiveStoppingPower;
  }

  // Add the missing endpoint material terms to
  // otherwise complete free/local Jacobians. Jraw already includes the
  // changing-mean q/p derivative. Qraw is zero-input process noise in the
  // uncorrected destination chart, so it receives ONLY the output congruence.
  // No clipping, symmetrisation or PSD repair is performed; covariance
  // validity remains the caller's responsibility.
  inline MaterialEndpointChartResult correctMaterialEndpointCharts(
      MaterialEndpointMatrix const& rawJacobian,
      MaterialEndpointMatrix const& rawNoise,
      MaterialEndpoint const& start,
      MaterialEndpoint const& end) {
    MaterialEndpointChartResult result;
    auto const validEndpoint = [](MaterialEndpoint const& endpoint) {
      return std::isfinite(endpoint.inverseMomentumFlow) && endpoint.localPhysicalDirection.allFinite() &&
             std::abs(endpoint.localPhysicalDirection.squaredNorm() - 1.) <= 1.e-10;
    };
    if (!rawJacobian.allFinite() || !rawNoise.allFinite() || !validEndpoint(start) || !validEndpoint(end))
      return result;
    for (unsigned int coordinate = 0; coordinate < 2; ++coordinate) {
      result.inputChart(0, 3 + coordinate) =
          -start.inverseMomentumFlow * start.localPhysicalDirection[coordinate];
      result.outputChart(0, 3 + coordinate) =
          end.inverseMomentumFlow * end.localPhysicalDirection[coordinate];
    }
    result.jacobian = result.outputChart * rawJacobian * result.inputChart;
    result.noise = result.outputChart * rawNoise * result.outputChart.transpose();
    result.valid = result.jacobian.allFinite() && result.noise.allFinite();
    result.status = result.valid ? MaterialEndpointChartResult::success : MaterialEndpointChartResult::nonfiniteResult;
    return result;
  }
}  // namespace shift

#endif
