// Code style: human-code.
#ifndef PhysicsTools_ShiftMuonSegments_DimuonInvariantMass_h
#define PhysicsTools_ShiftMuonSegments_DimuonInvariantMass_h

#include <cmath>

namespace shift {
  template <typename Vector>
  inline double dimuonInvariantMass(Vector const& first, Vector const& second, double muonMass = 0.105658) {
    // GlobalVector stores floats. Promote before every product and sum.
    double const x1 = first.x(), y1 = first.y(), z1 = first.z();
    double const x2 = second.x(), y2 = second.y(), z2 = second.z();
    double const p1 = std::hypot(x1, y1, z1);
    double const p2 = std::hypot(x2, y2, z2);
    double const e1 = std::hypot(p1, muonMass);
    double const e2 = std::hypot(p2, muonMass);
    double const mass2 = muonMass * muonMass;
    if (p1 == 0. || p2 == 0.)
      return std::sqrt(2. * (mass2 + e1 * e2));

    // Separate E1 E2 - p1.p2 into positive radial and angular terms. This
    // avoids subtracting TeV-scale squared energies to obtain a small mass.
    double const radial = mass2 * (p1 * p1 + p2 * p2 + mass2) / (e1 * e2 + p1 * p2);
    double const dx = x1 / p1 - x2 / p2;
    double const dy = y1 / p1 - y2 / p2;
    double const dz = z1 / p1 - z2 / p2;
    double const angular = 0.5 * p1 * p2 * (dx * dx + dy * dy + dz * dz);
    return std::sqrt(2. * (mass2 + radial + angular));
  }
}
#endif
