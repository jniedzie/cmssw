#ifndef PhysicsTools_ShiftMuonSegments_CubatureTransport_h
#define PhysicsTools_ShiftMuonSegments_CubatureTransport_h
#include <Eigen/Cholesky>
#include <array>
#include <cmath>

namespace shift {
  struct TransportMoments {
    bool valid = false;
    Eigen::Matrix<double,5,1> mean = Eigen::Matrix<double,5,1>::Zero();
    Eigen::Matrix<double,5,5> covariance = Eigen::Matrix<double,5,5>::Zero();
  };
  // Third-degree spherical cubature. The callback supplies conditional mean
  // and process noise. The total covariance is E[Q]+Cov[E[x|input]], so noise
  // accumulated in earlier segments can sample different material paths.
  // All points have positive equal weight; none may be silently discarded.
  template<class Transport>
  TransportMoments cubatureTransport(TransportMoments const& input, Transport const& transport) {
    TransportMoments output;
    if (!input.valid || !input.mean.allFinite() || !input.covariance.allFinite()) return output;
    Eigen::LLT<Eigen::Matrix<double,5,5>> factor(input.covariance);
    if (factor.info()!=Eigen::Success) return output;
    Eigen::Matrix<double,5,5> const spread=std::sqrt(5.)*factor.matrixL().toDenseMatrix();
    std::array<TransportMoments,10> points;
    for (unsigned int j=0;j<10;++j) {
      Eigen::Matrix<double,5,1> const state=input.mean+(j%2 ? -1. : 1.)*spread.col(j/2);
      points[j]=transport(state);
      if (!points[j].valid || !points[j].mean.allFinite() || !points[j].covariance.allFinite()) return output;
      output.mean+=points[j].mean/10.;
    }
    for (auto const& point:points) {
      Eigen::Matrix<double,5,1> const residual=point.mean-output.mean;
      output.covariance+=(point.covariance+residual*residual.transpose())/10.;
    }
    output.valid=output.mean.allFinite() && output.covariance.allFinite() &&
        output.covariance.llt().info()==Eigen::Success;
    return output;
  }
}
#endif
