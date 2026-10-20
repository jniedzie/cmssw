// Code style: human-code.
#ifndef PhysicsTools_ShiftMuonSegments_ExactMuonFitInput_h
#define PhysicsTools_ShiftMuonSegments_ExactMuonFitInput_h

#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"

#include <cmath>
#include <typeinfo>

namespace shift {
  // Only remove copies of a supported leaf measurement. Channel sharing or
  // spatial proximity alone does not identify the same detector pulse.
  inline bool isExactCSCMeasurementCopy(TrackingRecHit const& first, TrackingRecHit const& second) {
    if (!first.isValid() || !second.isValid() || first.geographicalId().rawId() == 0 ||
        first.geographicalId() != second.geographicalId() || first.dimension() != 2 || second.dimension() != 2 ||
        typeid(first) != typeid(CSCRecHit2D) || typeid(second) != typeid(CSCRecHit2D))
      return false;

    auto const firstParameters = first.parameters();
    auto const secondParameters = second.parameters();
    auto const firstError = first.parametersError();
    auto const secondError = second.parametersError();
    for (int row = 0; row < 2; ++row) {
      if (!std::isfinite(firstParameters[row]) || !std::isfinite(secondParameters[row]) ||
          firstParameters[row] != secondParameters[row])
        return false;
      for (int column = 0; column < 2; ++column)
        if (!std::isfinite(firstError[row][column]) || !std::isfinite(secondError[row][column]) ||
            firstError[row][column] != secondError[row][column])
          return false;
    }

    auto const& a = static_cast<CSCRecHit2D const&>(first);
    auto const& b = static_cast<CSCRecHit2D const&>(second);
    auto equalFinite = [](double x, double y) { return std::isfinite(x) && std::isfinite(y) && x == y; };
    if (!equalFinite(a.localPosition().z(), b.localPosition().z()) || !equalFinite(a.tpeak(), b.tpeak()) ||
        !equalFinite(a.positionWithinStrip(), b.positionWithinStrip()) ||
        !equalFinite(a.errorWithinStrip(), b.errorWithinStrip()) || !equalFinite(a.wireTime(), b.wireTime()) ||
        !equalFinite(a.energyDepositedInLayer(), b.energyDepositedInLayer()) || a.quality() != b.quality() ||
        a.badStrip() != b.badStrip() || a.badWireGroup() != b.badWireGroup() || a.hitWire() != b.hitWire() ||
        a.wgroupsBX() != b.wgroupsBX() || a.nWireGroups() != b.nWireGroups() || a.nStrips() != b.nStrips() ||
        a.nTimeBins() != b.nTimeBins() || a.nStrips() == 0 || a.nStrips() > CSCRecHit2D::MAXSTRIPS ||
        a.nTimeBins() == 0 || a.nTimeBins() > CSCRecHit2D::MAXTIMEBINS)
      return false;
    for (unsigned int strip = 0; strip < a.nStrips(); ++strip) {
      if (a.channels(strip) != b.channels(strip) || a.channelsl1a(strip) != b.channelsl1a(strip))
        return false;
      for (unsigned int bin = 0; bin < a.nTimeBins(); ++bin)
        if (!equalFinite(a.adcs(strip, bin), b.adcs(strip, bin)))
          return false;
    }
    // CSC sharesInput assumes the other hit is a CSCRecHit2D. The exact type
    // checks above must precede these calls; require symmetry as well.
    return a.sharesInput(&b, TrackingRecHit::all) && b.sharesInput(&a, TrackingRecHit::all);
  }
}  // namespace shift

#endif
