#include "PhysicsTools/ShiftMuonSegments/interface/MaterialEndpointChart.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>

namespace {
  using Vector = Eigen::Matrix<double, 5, 1>;
  using Matrix = shift::MaterialEndpointMatrix;
  using Derivatives = Eigen::Matrix<double, 3, 5>;
  struct Plane {
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    Eigen::Matrix3d axes = Eigen::Matrix3d::Identity();
  };
  struct StraightMaterial {
    Plane source;
    Plane target;
    double localPzSign = 1.;
    double momentumLossPerCm = .007;  // dp/ds = -b, physical-direction convention

    Eigen::Vector3d localDirection(Vector const& state) const {
      return localPzSign * Eigen::Vector3d(state[1], state[2], 1.).normalized();
    }
    Eigen::Vector3d position(Vector const& state) const {
      return source.origin + source.axes.col(0) * state[3] + source.axes.col(1) * state[4];
    }
    double path(Vector const& state) const {
      Eigen::Vector3d const normal = target.axes.col(2);
      return normal.dot(target.origin - position(state)) / normal.dot(source.axes * localDirection(state));
    }
    Vector mean(Vector const& state) const {
      Eigen::Vector3d const direction = source.axes * localDirection(state);
      double const distance = path(state);
      double const finalMomentum = std::abs(1. / state[0]) - momentumLossPerCm * distance;
      assert(finalMomentum > 0.);
      Eigen::Vector3d const endPosition = target.axes.transpose() * (position(state) + distance * direction - target.origin);
      Eigen::Vector3d const endDirection = target.axes.transpose() * direction;
      Vector output;
      output << std::copysign(1. / finalMomentum, state[0]), endDirection.x() / endDirection.z(),
          endDirection.y() / endDirection.z(), endPosition.x(), endPosition.y();
      return output;
    }
    // Independent differentiation of the exact plane-intersection mean,
    // including path-dependent energy loss. No chart correction is used.
    Matrix exactJacobian(Vector const& state) const {
      Eigen::Vector3d const local = localDirection(state);
      Eigen::Vector3d const direction = source.axes * local;
      Eigen::Vector3d const normal = target.axes.col(2);
      double const distance = path(state);
      double const p = std::abs(1. / state[0]);
      double const finalP = p - momentumLossPerCm * distance;
      double const charge = std::copysign(1., state[0]);
      Derivatives dr = Derivatives::Zero(), du = Derivatives::Zero();
      dr.col(3) = source.axes.col(0);
      dr.col(4) = source.axes.col(1);
      Eigen::Vector3d const v(state[1], state[2], 1.);
      double const norm = v.norm();
      du.col(1) = source.axes * (localPzSign / norm * (Eigen::Vector3d::UnitX() - v * state[1] / v.squaredNorm()));
      du.col(2) = source.axes * (localPzSign / norm * (Eigen::Vector3d::UnitY() - v * state[2] / v.squaredNorm()));
      Eigen::Matrix<double, 1, 5> const dLength = -normal.transpose() * (dr + distance * du) / normal.dot(direction);
      Matrix jacobian = Matrix::Zero();
      jacobian.row(0) = charge * momentumLossPerCm / (finalP * finalP) * dLength;
      jacobian(0, 0) += p * p / (finalP * finalP);
      Derivatives const targetDu = target.axes.transpose() * du;
      Eigen::Vector3d const targetDirection = target.axes.transpose() * direction;
      jacobian.row(1) = (targetDu.row(0) - targetDirection.x() / targetDirection.z() * targetDu.row(2)) / targetDirection.z();
      jacobian.row(2) = (targetDu.row(1) - targetDirection.y() / targetDirection.z() * targetDu.row(2)) / targetDirection.z();
      jacobian.bottomRows<2>() = (target.axes.transpose() * (dr + direction * dLength + distance * du)).topRows<2>();
      return jacobian;
    }
    // Raw free transport between transverse planes: the first-order path
    // length is fixed. The incoming/outgoing geometric chart changes project
    // positions but OMIT the material q/p change caused by those projections.
    // This construction is independent of both the corrected and exact J.
    Matrix rawJacobian(Vector const& state) const {
      Eigen::Vector3d const direction = source.axes * localDirection(state);
      Eigen::Vector3d const normal = target.axes.col(2);
      Derivatives dr = Derivatives::Zero(), du = Derivatives::Zero();
      dr.col(3) = source.axes.col(0);
      dr.col(4) = source.axes.col(1);
      Eigen::Vector3d const v(state[1], state[2], 1.);
      du.col(1) = source.axes * (localPzSign / v.norm() * (Eigen::Vector3d::UnitX() - v * state[1] / v.squaredNorm()));
      du.col(2) = source.axes * (localPzSign / v.norm() * (Eigen::Vector3d::UnitY() - v * state[2] / v.squaredNorm()));
      Derivatives const transverseStart = (Eigen::Matrix3d::Identity() - direction * direction.transpose()) * dr;
      Derivatives const transverseEnd = transverseStart + path(state) * du;
      Derivatives const surfaceEnd = (Eigen::Matrix3d::Identity() - direction * normal.transpose() / normal.dot(direction)) * transverseEnd;
      Eigen::Vector3d const targetDirection = target.axes.transpose() * direction;
      Derivatives const targetDu = target.axes.transpose() * du;
      double const p = std::abs(1. / state[0]);
      double const finalP = p - momentumLossPerCm * path(state);
      Matrix jacobian = Matrix::Zero();
      jacobian(0, 0) = p * p / (finalP * finalP);  // Already present in free transport!
      jacobian.row(1) = (targetDu.row(0) - targetDirection.x() / targetDirection.z() * targetDu.row(2)) / targetDirection.z();
      jacobian.row(2) = (targetDu.row(1) - targetDirection.y() / targetDirection.z() * targetDu.row(2)) / targetDirection.z();
      jacobian.bottomRows<2>() = (target.axes.transpose() * surfaceEnd).topRows<2>();
      return jacobian;
    }
    shift::MaterialEndpoint endpoint(Vector const& state, bool end) const {
      double const charge = std::copysign(1., state[0]);
      double const p = std::abs(1. / state[0]) - (end ? momentumLossPerCm * path(state) : 0.);
      double const mass = .105658;
      double const stoppingPower = p / std::hypot(p, mass) * momentumLossPerCm;
      Eigen::Vector3d const local = end ? Eigen::Vector3d(target.axes.transpose() * source.axes * localDirection(state)) : localDirection(state);
      return {shift::materialInverseMomentumFlow(charge, p, mass, stoppingPower), local};
    }
  };

  void compare(Matrix const& a, Matrix const& b, double relativeTolerance, char const* label) {
    double const relative = (a - b).norm() / std::max(1., b.norm());
    if (!(relative <= relativeTolerance)) {
      std::cerr << label << " relative difference " << relative << "\nactual:\n" << a << "\nexpected:\n" << b << '\n';
      std::abort();
    }
  }

  void analyticPlanes() {
    unsigned int cases = 0;
    double worstExact = 0., worstNumerical = 0.;
    for (double charge : {-1., 1.}) for (double sign : {-1., 1.})
      for (double travelSign : {-1., 1.}) for (bool rotated : {false, true}) {
      StraightMaterial fixture;
      fixture.localPzSign = sign;
      Vector state;
      state << charge / 7., .35, -.22, 1.2, -2.1;
      if (rotated) {
        fixture.source.origin << 5., -3., 8.;
        fixture.source.axes = (Eigen::AngleAxisd(.3, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(-.2, Eigen::Vector3d::UnitX())).toRotationMatrix();
        fixture.target.axes = (Eigen::AngleAxisd(-.4, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(.1, Eigen::Vector3d::UnitZ())).toRotationMatrix();
        fixture.target.origin = fixture.position(state) + travelSign * 320. * fixture.source.axes * fixture.localDirection(state);
      } else {
        fixture.target.origin.z() = travelSign * sign * 300.;
      }
      auto const corrected = shift::correctMaterialEndpointCharts(fixture.rawJacobian(state), Matrix::Zero(), fixture.endpoint(state, false), fixture.endpoint(state, true));
      assert(corrected.valid && corrected.status == shift::MaterialEndpointChartResult::success);
      Matrix const exact = fixture.exactJacobian(state);
      compare(corrected.jacobian, exact, 3.e-15, "Exact plane/material derivative");
      assert((corrected.jacobian.row(0) - exact.row(0)).norm() < 1.e-13);
      worstExact = std::max(worstExact, (corrected.jacobian - exact).norm() / exact.norm());
      Matrix numerical;
      for (unsigned int coordinate = 0; coordinate < 5; ++coordinate) {
        double const h = coordinate == 0 ? 1.e-6 : 1.e-5;
        Vector plus = state, minus = state;
        plus[coordinate] += h;
        minus[coordinate] -= h;
        numerical.col(coordinate) = (fixture.mean(plus) - fixture.mean(minus)) / (2. * h);
      }
      compare(corrected.jacobian, numerical, 1.e-9, "Numerical plane/material derivative");
      assert((corrected.jacobian.row(0) - numerical.row(0)).norm() < 1.e-8);
      worstNumerical = std::max(worstNumerical, (corrected.jacobian - numerical).norm() / numerical.norm());
      assert((fixture.rawJacobian(state).row(0) - exact.row(0)).norm() > .005);
      assert(corrected.jacobian(0, 0) == fixture.rawJacobian(state)(0, 0));
      if (!rotated) {
        assert(std::abs(corrected.jacobian(0, 3)) < 1.e-17);
        assert(std::abs(corrected.jacobian(0, 4)) < 1.e-17);
        double const a = fixture.endpoint(state, true).inverseMomentumFlow;
        assert(std::abs(corrected.jacobian(0, 1) - a * travelSign * 300. * state[1] / std::sqrt(1. + state[1]*state[1] + state[2]*state[2])) < 1.e-15);
      }
      ++cases;
    }
    std::cout << "Both charges, both local-pz signs, forward/backward signed paths, global-z/rotated planes: " << cases
              << " cases PASS; exact relative error=" << worstExact << ", numerical=" << worstNumerical << '\n';
  }

  void identitiesAndNoise() {
    Matrix raw;
    raw << 2., .1, -.2, 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0.,
        0., 100., 2., 1., 0., 0., -3., 120., 0., 1.;
    Matrix factor;
    factor << .01, 0., 0., 0., 0., .002, .03, 0., 0., 0., -.001, .004, .02, 0., 0.,
        .1, 2., .3, 1., 0., -.2, -.7, 1.1, .4, .8;
    for (bool rankDeficient : {false, true}) {
      if (rankDeficient) factor.col(4).setZero();
      Matrix const noise = factor * factor.transpose();
      shift::MaterialEndpoint start{.003, Eigen::Vector3d(.3, -.2, 1.).normalized()};
      shift::MaterialEndpoint end{-.008, Eigen::Vector3d(-.4, .25, -1.).normalized()};
      auto const result = shift::correctMaterialEndpointCharts(raw, noise, start, end);
      assert(result.valid);
      Matrix const transformedFactor = result.outputChart * factor;
      compare(result.noise, transformedFactor * transformedFactor.transpose(), 2.e-15, "Output-only noise congruence");
      assert((result.noise - result.noise.transpose()).norm() < 1.e-13);
      assert(result.noise.selfadjointView<Eigen::Lower>().eigenvalues().minCoeff() >= -1.e-13);
      assert((result.noise - result.jacobian * noise * result.jacobian.transpose()).norm() > 1.);
      auto const samePlane = shift::correctMaterialEndpointCharts(Matrix::Identity(), noise, start, start);
      compare(samePlane.jacobian, Matrix::Identity(), 0., "Same-plane chart cancellation");
      start.inverseMomentumFlow = 0.;
      end.inverseMomentumFlow = 0.;
      auto const vacuum = shift::correctMaterialEndpointCharts(raw, noise, start, end);
      compare(vacuum.jacobian, raw, 0., "Vacuum Jacobian identity");
      compare(vacuum.noise, noise, 0., "Vacuum noise identity");
    }
    shift::MaterialEndpoint const normal{.3, Eigen::Vector3d::UnitZ()};
    auto const normalIncidence = shift::correctMaterialEndpointCharts(raw, Matrix::Zero(), normal, normal);
    compare(normalIncidence.jacobian, raw, 0., "Normal-incidence identity");
    assert(normalIncidence.noise.isZero());
    std::cout << "Vacuum/normal/same-plane identities, PSD and rank-deficient noise, no double mean correction: PASS\n";
  }

  void invalidInputs() {
    shift::MaterialEndpoint endpoint;
    Matrix matrix = Matrix::Identity();
    endpoint.inverseMomentumFlow = std::numeric_limits<double>::quiet_NaN();
    assert(!shift::correctMaterialEndpointCharts(matrix, matrix, endpoint, {}).valid);
    endpoint = {};
    endpoint.localPhysicalDirection *= 2.;
    assert(!shift::correctMaterialEndpointCharts(matrix, matrix, {}, endpoint).valid);
    matrix(0, 0) = std::numeric_limits<double>::infinity();
    assert(!shift::correctMaterialEndpointCharts(matrix, Matrix::Zero(), {}, {}).valid);
    assert(!shift::correctMaterialEndpointCharts(Matrix::Identity(), matrix, {}, {}).valid);
    std::cout << "Nonfinite matrices/flow and nonunit endpoint directions rejected: PASS\n";
  }
}  // namespace

int main() {
  analyticPlanes();
  identitiesAndNoise();
  invalidInputs();
  std::cout << "Material endpoint chart regression suite passed\n";
}
