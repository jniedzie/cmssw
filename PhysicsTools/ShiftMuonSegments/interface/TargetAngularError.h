#ifndef PhysicsTools_ShiftMuonSegments_TargetAngularError_h
#define PhysicsTools_ShiftMuonSegments_TargetAngularError_h
#include <cmath>
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
namespace shift {
  struct TargetAngularError { double eta = -1.; double phi = -1.; };
  // Local parameters on a global-z plane are (q/p, px/pz, py/pz, x, y).
  // Evaluate the angular Jacobian at the SAME state as its covariance.
  inline TargetAngularError targetAngularError(double tx, double ty, double pz,
                                                AlgebraicSymMatrix55 const& c) {
    TargetAngularError result;
    double const r2 = tx * tx + ty * ty;
    if (!(r2 > 0.) || pz == 0.) return result;
    double const f = -std::copysign(1., pz) / (r2 * std::sqrt(1. + r2));
    double const ex = f * tx, ey = f * ty, px = -ty / r2, py = tx / r2;
    double const ve = ex*ex*c(1,1) + 2.*ex*ey*c(1,2) + ey*ey*c(2,2);
    double const vp = px*px*c(1,1) + 2.*px*py*c(1,2) + py*py*c(2,2);
    if (std::isfinite(ve) && ve >= 0.) result.eta = std::sqrt(ve);
    if (std::isfinite(vp) && vp >= 0.) result.phi = std::sqrt(vp);
    return result;
  }
}
#endif
