// Code style: human-code.
#ifndef PhysicsTools_ShiftMuonSegments_DetectorFitAudit_h
#define PhysicsTools_ShiftMuonSegments_DetectorFitAudit_h

#include <array>
#include <cmath>
#include <limits>

namespace shift {
  // Evaluation only. The caller supplies independently associated SimHit
  // endpoints in one detector-surface chart. No extrapolation is allowed.
  struct DetectorTruthEndpoint {
    std::array<double, 3> position{};
    std::array<double, 3> momentum{};
    double localZ = 0.;
    int kind = 0;  // 0 entry, 1 exit; exit uses the SimHit entry momentum proxy.
  };

  struct DetectorTruthBracket {
    bool hasLower = false;
    bool hasUpper = false;
    DetectorTruthEndpoint lower;
    DetectorTruthEndpoint upper;

    void add(DetectorTruthEndpoint const& point) {
      if (!std::isfinite(point.localZ))
        return;
      for (unsigned int i = 0; i < 3; ++i)
        if (!std::isfinite(point.position[i]) || !std::isfinite(point.momentum[i]))
          return;
      if (point.localZ <= 0. && (!hasLower || point.localZ > lower.localZ)) {
        lower = point;
        hasLower = true;
      }
      if (point.localZ >= 0. && (!hasUpper || point.localZ < upper.localZ)) {
        upper = point;
        hasUpper = true;
      }
    }

    bool interpolate(std::array<double, 3>& position,
                     std::array<double, 3>& momentum,
                     double& weight) const {
      if (!hasLower || !hasUpper)
        return false;
      double const span = upper.localZ - lower.localZ;
      if (!(span >= 0.) || !std::isfinite(span))
        return false;
      weight = span > 0. ? -lower.localZ / span : 0.;
      if (!(weight >= 0. && weight <= 1.) || !std::isfinite(weight))
        return false;
      double p2 = 0.;
      for (unsigned int i = 0; i < 3; ++i) {
        position[i] = lower.position[i] + weight * (upper.position[i] - lower.position[i]);
        momentum[i] = lower.momentum[i] + weight * (upper.momentum[i] - lower.momentum[i]);
        if (!std::isfinite(position[i]) || !std::isfinite(momentum[i]))
          return false;
        p2 += momentum[i] * momentum[i];
      }
      return p2 > 0. && std::isfinite(p2);
    }
  };
}  // namespace shift
#endif
