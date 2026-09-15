#include "PhysicsTools/ShiftMuonSegments/interface/DimuonKinematicCovariance.h"
#include <stdexcept>

int main() {
  using namespace shift;
  std::array<VertexTrack,2> tracks;
  tracks[0].parameters << .01, .02, -.03, 0., 0.;
  tracks[1].parameters << -.02, -.01, .04, 0., 0.;
  VertexTrackParameters scales; scales << 1.e-5, 1.e-4, 1.e-4, .01, .01;
  for (auto& track : tracks) track.covariance=scales.array().square().matrix().asDiagonal();
  VertexTrackCovariance cross=.3*tracks[0].covariance;
  auto check=[](bool value) { if (!value) throw std::runtime_error("Dimuon covariance regression"); };
  for (std::array<double,2> signs : {std::array<double,2>{-1.,-1.},std::array<double,2>{1.,-1.}}) {
    auto const result=dimuonKinematicCovariance(tracks,cross,signs);
    check(result.valid);
    Eigen::Matrix<double,10,10> covariance=Eigen::Matrix<double,10,10>::Zero();
    covariance.topLeftCorner<5,5>()=tracks[0].covariance;
    covariance.bottomRightCorner<5,5>()=tracks[1].covariance;
    covariance.topRightCorner<5,5>()=cross; covariance.bottomLeftCorner<5,5>()=cross.transpose();
    // Differentiate the physical two-four-vector observables independently
    // of the analytic uncertainty Jacobian, including cross-track covariance.
    auto observe=[&](std::array<VertexTrack,2> const& values) {
      Eigen::Vector3d sum=Eigen::Vector3d::Zero(); double energy=0.;
      for (unsigned int i=0;i<2;++i) {
        auto const& a=values[i].parameters;
        Eigen::Vector3d p(a[1],a[2],1.);p*=signs[i]/(std::abs(a[0])*p.norm());
        sum+=p;energy+=std::sqrt(p.squaredNorm()+.1056583745*.1056583745);
      }
      VertexTrackParameters value;double const pt=sum.head<2>().norm();
      value << std::sqrt(energy*energy-sum.squaredNorm()),pt,sum[2],std::asinh(sum[2]/pt),std::atan2(sum[1],sum[0]);
      return value;
    };
    Eigen::Matrix<double,5,10> jacobian;
    for (unsigned int j=0;j<10;++j) {
      auto plus=tracks,minus=tracks;double const h=scales[j%5]*.01;
      plus[j/5].parameters[j%5]+=h;minus[j/5].parameters[j%5]-=h;
      jacobian.col(j)=(observe(plus)-observe(minus))/(2.*h);
    }
    auto const expected=(jacobian*covariance*jacobian.transpose()).eval();
    for (unsigned int i=0;i<5;++i) for (unsigned int j=0;j<5;++j)
      check(std::abs(expected(i,j)-result.covariance(i,j)) < 1.e-5*std::sqrt(expected(i,i)*expected(j,j)));
    check(std::abs(result.minQoverPSignificance-1000.)<1.e-8);
    auto const independent=dimuonKinematicCovariance(tracks,VertexTrackCovariance::Zero(),signs);
    check(std::abs(independent.covariance(0,0)-result.covariance(0,0))>1.e-6);
  }
  tracks[0].parameters[0]=0.;
  check(!dimuonKinematicCovariance(tracks,cross,{-1.,-1.}).valid);
}
