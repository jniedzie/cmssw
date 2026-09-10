#ifndef PhysicsTools_ShiftMuonSegments_HitTruthAssociation_h
#define PhysicsTools_ShiftMuonSegments_HitTruthAssociation_h

#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <vector>

namespace shift {
  struct TruthMeasurement {
    uint32_t layer;
    double x, y, z, xx, xy, yy;
  };
  struct TruthCrossing {
    uint32_t layer, track;
    double x, y, z, px, py, pz;
  };
  struct HitTruthAssociation {
    int track = -1;
    unsigned int total = 0, matched = 0, ambiguous = 0;
    double purity = 0.;
  };

  // Diagnostic policy: at least three distinct layers, >=75% uniquely
  // compatible layers and a unique winner. Chi2=25 is a compatibility gate,
  // not a nearest-truth choice: two compatible tracks cast no vote.
  inline HitTruthAssociation associateHits(std::vector<TruthMeasurement> const& measurements,
                                          std::vector<TruthCrossing> const& crossings) {
    HitTruthAssociation result;
    std::map<uint32_t, std::vector<TruthMeasurement>> layers;
    for (auto const& hit : measurements) {
      auto& existing = layers[hit.layer];
      bool duplicate = false;
      for (auto const& other : existing)
        duplicate |= hit.x == other.x && hit.y == other.y && hit.z == other.z;
      if (!duplicate)
        existing.push_back(hit);
    }
    std::map<uint32_t, unsigned int> votes;
    for (auto const& [layer, hits] : layers) {
      ++result.total;
      if (hits.size() != 1) {
        ++result.ambiguous;
        continue;
      }
      auto const& h = hits.front();
      double const determinant = h.xx * h.yy - h.xy * h.xy;
      if (!(h.xx > 0. && h.yy > 0. && determinant > 0.) || !std::isfinite(determinant))
        continue;
      std::set<uint32_t> compatible;
      for (auto const& s : crossings) {
        if (s.layer != layer || !std::isfinite(s.pz) || std::abs(s.pz) < 1.e-12)
          continue;
        double const dz = h.z - s.z;
        double const dx = h.x - s.x - dz * s.px / s.pz;
        double const dy = h.y - s.y - dz * s.py / s.pz;
        double const chi2 = (h.yy * dx * dx - 2. * h.xy * dx * dy + h.xx * dy * dy) / determinant;
        if (std::isfinite(chi2) && chi2 >= 0. && chi2 <= 25.)
          compatible.insert(s.track);
      }
      if (compatible.size() == 1)
        ++votes[*compatible.begin()];
      else if (compatible.size() > 1)
        ++result.ambiguous;
    }
    int winner = -1;
    bool tied = false;
    for (auto const& [track, count] : votes) {
      if (count > result.matched) {
        result.matched = count;
        winner = track;
        tied = false;
      } else if (count == result.matched)
        tied = true;
    }
    result.purity = result.total ? double(result.matched) / result.total : 0.;
    if (!tied && result.matched >= 3 && result.purity >= .75)
      result.track = winner;
    return result;
  }
}
#endif
