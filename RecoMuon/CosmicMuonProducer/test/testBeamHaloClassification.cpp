// Code style: human-code.
#include "RecoMuon/CosmicMuonProducer/interface/BeamHaloClassification.h"

#include <cstdlib>
#include <limits>

namespace {
  void require(bool condition) {
    if (!condition)
      std::abort();
  }
}

int main() {
  using cosmic_muon::isBeamHaloForSmoothing;
  double const values[] = {-8., -4.01, -4., -1., 0., 1., 4., 4.01, 8.};
  for (double first : values) {
    for (double last : values) {
      require(isBeamHaloForSmoothing(first, last, false) == (first > 4. || last > 4.));
      require(isBeamHaloForSmoothing(first, last, true) ==
              isBeamHaloForSmoothing(-first, -last, true));
      if (first >= 0. && last >= 0.)
        require(isBeamHaloForSmoothing(first, last, false) == isBeamHaloForSmoothing(first, last, true));
      if (std::abs(first) <= 4. && std::abs(last) <= 4.)
        require(!isBeamHaloForSmoothing(first, last, true));
    }
  }
  require(!isBeamHaloForSmoothing(-5., -6., false));
  require(isBeamHaloForSmoothing(-5., -6., true));
  require(isBeamHaloForSmoothing(-1., -5., true));
  require(!isBeamHaloForSmoothing(-4., 4., true));
  double const nan = std::numeric_limits<double>::quiet_NaN();
  require(!isBeamHaloForSmoothing(nan, nan, false));
  require(!isBeamHaloForSmoothing(nan, nan, true));
}
