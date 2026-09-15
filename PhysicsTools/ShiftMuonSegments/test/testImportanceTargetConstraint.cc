#include "PhysicsTools/ShiftMuonSegments/interface/ImportanceTargetConstraint.h"
#include <stdexcept>

int main() {
  using namespace shift;
  using Matrix=Eigen::Matrix<double,5,5>;
  using Vector=Eigen::Matrix<double,5,1>;
  TransportMoments input;input.valid=true;input.mean<<.02,.01,-.02,.8,-.4;
  Vector scales,axis;scales<<.001,.0002,.0003,.5,.4;axis<<1.,.3,-.2,.1,-.4;
  input.covariance=scales.asDiagonal()*(Matrix::Identity()+axis*axis.transpose())*scales.asDiagonal();
  Matrix jacobian=Matrix::Identity();jacobian(0,0)=1.4;jacobian(1,0)=.4;
  jacobian(3,1)=-14000.;jacobian(4,2)=-14000.;
  Matrix noise=Matrix::Zero();noise.diagonal()<<1.e-6,1.e-6,1.e-6,100.,100.;
  Vector correlation;correlation<<.001,.001,.002,1.,1.;noise+=correlation*correlation.transpose();
  auto transport=[&](Vector const& x) {return TransportMoments{true,jacobian*x,noise};};
  Vector const priorMean=jacobian*input.mean;
  Matrix const priorCovariance=jacobian*input.covariance*jacobian.transpose()+noise;
  Eigen::Matrix2d hit=.01*Eigen::Matrix2d::Identity();
  Eigen::LLT<Eigen::Matrix2d> solve(priorCovariance.bottomRightCorner<2,2>()+hit);
  Eigen::Matrix<double,5,2> gain=solve.solve(priorCovariance.bottomRows<2>()).transpose();
  Vector const expectedMean=priorMean-gain*priorMean.tail<2>();
  Matrix const expectedCovariance=priorCovariance-gain*priorCovariance.bottomRows<2>();
  for(unsigned int clouds:{1u,2u,4u}) {
    auto result=importanceTargetConstraint(input,transport,.1,.1,0.,clouds);
    if(!result.valid || (result.posterior.mean-expectedMean).norm()>1.e-9 ||
       (result.posterior.covariance-expectedCovariance).norm()>1.e-9 ||
       std::abs(result.effectivePoints-10.*clouds)>1.e-9)
      throw std::runtime_error("Importance proposal must reproduce Gaussian conditioning without double counting the hit");
  }
  auto failed=importanceTargetConstraint(input,[](Vector const&){return TransportMoments{};},.1,.1,0.);
  if(failed.valid) throw std::runtime_error("A failed quadrature propagation cannot be silently omitted");
}
