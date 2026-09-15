#ifndef PhysicsTools_ShiftMuonSegments_DimuonKinematicCovariance_h
#define PhysicsTools_ShiftMuonSegments_DimuonKinematicCovariance_h

#include "PhysicsTools/ShiftMuonSegments/interface/CommonVertexRefit.h"
#include <limits>

namespace shift {
  struct DimuonKinematicCovariance {
    // Observable order: mass, pT, pZ, eta, phi. Momentum units are GeV.
    bool valid = false;
    VertexTrackParameters observables = VertexTrackParameters::Zero();
    VertexTrackCovariance covariance = VertexTrackCovariance::Zero();
    double minQoverPSignificance = -1.;
  };

  inline DimuonKinematicCovariance dimuonKinematicCovariance(
      std::array<VertexTrack, 2> const& tracks,
      VertexTrackCovariance const& crossCovariance,
      std::array<double, 2> const& pzSigns,
      double mass = .1056583745) {
    DimuonKinematicCovariance result;
    std::array<Eigen::Vector3d, 2> momenta;
    std::array<double, 2> energies;
    Eigen::Vector3d total = Eigen::Vector3d::Zero();
    double energy = 0.;
    double minSignificance = std::numeric_limits<double>::infinity();
    for (unsigned int i=0; i<2; ++i) {
      auto const& a = tracks[i].parameters;
      if (!a.allFinite() || !tracks[i].covariance.allFinite() || a[0] == 0. ||
          !(tracks[i].covariance(0,0)>0.) || pzSigns[i] == 0.) return result;
      double const pz = std::copysign(1.,pzSigns[i]) /
          (std::abs(a[0])*std::sqrt(1.+a[1]*a[1]+a[2]*a[2]));
      momenta[i] = Eigen::Vector3d(a[1]*pz,a[2]*pz,pz);
      energies[i] = std::sqrt(momenta[i].squaredNorm()+mass*mass);
      total += momenta[i]; energy += energies[i];
      minSignificance = std::min(minSignificance,
          std::abs(a[0])/std::sqrt(tracks[i].covariance(0,0)));
    }
    result.minQoverPSignificance = minSignificance;
    double const pt2=total.head<2>().squaredNorm(), pt=std::sqrt(pt2), p=total.norm();
    double const m2=energy*energy-total.squaredNorm();
    if (!(pt>0.) || !(m2>0.) || !crossCovariance.allFinite()) return result;
    double const m=std::sqrt(m2);
    result.observables << m, pt, total[2], std::asinh(total[2]/pt), std::atan2(total[1],total[0]);
    std::array<VertexTrackCovariance,2> jacobian;
    for (unsigned int i=0;i<2;++i) {
      jacobian[i].setZero();
      auto const& a=tracks[i].parameters;
      auto const& v=momenta[i];
      double const d=1.+a[1]*a[1]+a[2]*a[2];
      std::array<Eigen::Vector3d,3> const derivative{
          -v/a[0], Eigen::Vector3d(v[2],0.,0.)-v*a[1]/d,
          Eigen::Vector3d(0.,v[2],0.)-v*a[2]/d};
      for (unsigned int j=0;j<3;++j) {
        auto const& dv=derivative[j];
        double const dptNumerator=total[0]*dv[0]+total[1]*dv[1];
        jacobian[i](0,j)=(energy*v.dot(dv)/energies[i]-total.dot(dv))/m;
        jacobian[i](1,j)=dptNumerator/pt;
        jacobian[i](2,j)=dv[2];
        jacobian[i](3,j)=dv[2]/p-total[2]*dptNumerator/(p*pt2);
        jacobian[i](4,j)=(total[0]*dv[1]-total[1]*dv[0])/pt2;
      }
    }
    result.covariance=jacobian[0]*tracks[0].covariance*jacobian[0].transpose()+
        jacobian[1]*tracks[1].covariance*jacobian[1].transpose();
    VertexTrackCovariance const cross=jacobian[0]*crossCovariance*jacobian[1].transpose();
    result.covariance+=cross+cross.transpose();
    result.valid=result.observables.allFinite() && result.covariance.allFinite() &&
        (result.covariance.diagonal().array()>=0.).all();
    return result;
  }
}
#endif
