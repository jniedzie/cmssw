// Code style: AI-code. AI-generated helpers for opt-in measured-track validation.
#ifndef PhysicsTools_ShiftMuonSegments_MeasuredMuonClonePolicy_h
#define PhysicsTools_ShiftMuonSegments_MeasuredMuonClonePolicy_h

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <tuple>
#include <vector>

namespace shift {
  struct MuonHitFingerprint {
    uint32_t detId;
    float localX;
    float localY;
  };

  inline std::vector<MuonHitFingerprint> distinctMuonHitFingerprints(std::vector<MuonHitFingerprint> hits) {
    hits.erase(std::remove_if(hits.begin(), hits.end(), [](auto const& hit) {
      return !std::isfinite(hit.localX) || !std::isfinite(hit.localY);
    }), hits.end());
    auto key = [](auto const& hit) { return std::tie(hit.detId, hit.localX, hit.localY); };
    std::sort(hits.begin(), hits.end(), [&](auto const& a, auto const& b) { return key(a) < key(b); });
    hits.erase(std::unique(hits.begin(), hits.end(), [&](auto const& a, auto const& b) {
      return key(a) == key(b);
    }), hits.end());
    return hits;
  }

  // Maximum one-to-one matching makes the fraction symmetric and <=1. Exact
  // repeated fingerprints cannot add votes. The position tolerance is the
  // existing Nano criterion; this is not full detector pulse/channel identity.
  inline bool shareMuonHitFingerprints(std::vector<MuonHitFingerprint> const& first,
                                      std::vector<MuonHitFingerprint> const& second,
                                      unsigned int minimumShared,
                                      double minimumFraction) {
    auto const a = distinctMuonHitFingerprints(first);
    auto const b = distinctMuonHitFingerprints(second);
    auto const smaller = std::min(a.size(), b.size());
    if (smaller == 0)
      return false;
    std::vector<int> matched(b.size(), -1);
    std::function<bool(unsigned int, std::vector<bool>&)> augment = [&](unsigned int i, std::vector<bool>& visited) {
      for (unsigned int j = 0; j < b.size(); ++j) {
        if (visited[j] || a[i].detId != b[j].detId ||
            std::hypot(a[i].localX - b[j].localX, a[i].localY - b[j].localY) >= 0.1)
          continue;
        visited[j] = true;
        if (matched[j] < 0 || augment(matched[j], visited)) {
          matched[j] = i;
          return true;
        }
      }
      return false;
    };
    unsigned int shared = 0;
    for (unsigned int i = 0; i < a.size(); ++i) {
      std::vector<bool> visited(b.size(), false);
      shared += augment(i, visited);
    }
    return shared >= minimumShared && static_cast<double>(shared)/smaller >= minimumFraction;
  }

  struct BothMuonCloneQuality {
    bool directionalRefitValid;
    int timingDirectionSign;
    double normalizedChi2;
    double ndof;
    unsigned int stations;
    unsigned int validHits;
    int source;
    unsigned int sourceIndex;
  };

  // Used only after a whole duplicate group has been restricted to measured
  // Both-endcap tracks. A single total order avoids pair-dependent preferences
  // cycling between Both/near tracks and strict/cosmic producers.
  inline bool betterBothMuonClone(BothMuonCloneQuality const& a, BothMuonCloneQuality const& b) {
    auto key = [](BothMuonCloneQuality const& candidate) {
      int const timingRank = candidate.timingDirectionSign > 0 ? 0 : candidate.timingDirectionSign == 0 ? 1 : 2;
      bool const finiteFit = candidate.ndof > 0. && std::isfinite(candidate.normalizedChi2) && candidate.normalizedChi2 >= 0.;
      double const chi2 = finiteFit ? candidate.normalizedChi2 : std::numeric_limits<double>::infinity();
      int const sourceRank = candidate.source == 1 ? 0 : candidate.source == 0 ? 1 : 2;
      return std::make_tuple(!candidate.directionalRefitValid, timingRank, !finiteFit, chi2,
                             -static_cast<int>(candidate.stations), -static_cast<int>(candidate.validHits),
                             sourceRank, candidate.sourceIndex);
    };
    return key(a) < key(b);
  }
}
#endif
