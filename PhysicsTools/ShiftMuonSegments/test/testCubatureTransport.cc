#include "PhysicsTools/ShiftMuonSegments/interface/CubatureTransport.h"
#include <stdexcept>
int main() {
  using namespace shift;
  TransportMoments input;
  input.valid=true; input.mean<<.02,.01,-.02,1.,2.;
  input.covariance.setIdentity();
  Eigen::Matrix<double,5,5> jacobian=Eigen::Matrix<double,5,5>::Identity();
  jacobian(3,1)=100.; jacobian(4,2)=100.; jacobian(1,0)=2.;
  Eigen::Matrix<double,5,5> noise=.2*Eigen::Matrix<double,5,5>::Identity();
  auto linear=cubatureTransport(input,[&](auto const& x){return TransportMoments{true,jacobian*x,noise};});
  if (!linear.valid || (linear.mean-jacobian*input.mean).norm()>1.e-12 ||
      (linear.covariance-jacobian*input.covariance*jacobian.transpose()-noise).norm()>1.e-10)
    throw std::runtime_error("Cubature must reproduce linear Gaussian transport");
  auto nonlinear=cubatureTransport(input,[&](auto const& x){
    TransportMoments r{true,x,noise}; r.mean[3]+=x[1]*x[1]; return r;
  });
  if (!nonlinear.valid || std::abs(nonlinear.mean[3]-input.mean[3]-input.mean[1]*input.mean[1]-1.)>1.e-12)
    throw std::runtime_error("Cubature must include nonlinear mean shift");
  auto failed=cubatureTransport(input,[](auto const&){return TransportMoments{};});
  if (failed.valid) throw std::runtime_error("Invalid quadrature point must reject entire transport");
}
