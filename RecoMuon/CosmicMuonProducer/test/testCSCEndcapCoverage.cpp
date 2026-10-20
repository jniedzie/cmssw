// Code style: human-code.
#include "RecoMuon/CosmicMuonProducer/interface/CSCEndcapCoverage.h"

#include <cstdlib>
#include <initializer_list>

namespace {
  void require(bool condition) {
    if (!condition)
      std::abort();
  }
}  // namespace

int main() {
  CSCDetId const positive(1, 1, 1, 1, 1);
  CSCDetId const negative(2, 1, 1, 1, 1);
  // Disabled always preserves the prior second-half build, including when
  // both CSC endcaps are already present. Enabled reproduces the prototype.
  for (unsigned int endcaps = 0; endcaps < 4; ++endcaps) {
    cosmic_muon::CSCEndcapCoverage selected;
    if (endcaps & 1)
      selected.addMeasurement(true, positive);
    if (endcaps & 2)
      selected.addMeasurement(true, negative);
    for (bool traversing : {false, true})
      for (bool strictOneLeg : {false, true}) {
        require(!selected.skipSecondHalf(false, traversing, strictOneLeg));
        require(selected.skipSecondHalf(true, traversing, strictOneLeg) ==
                (traversing && strictOneLeg && endcaps == 3));
      }
  }

  cosmic_muon::CSCEndcapCoverage coverage;
  require(!coverage.skipSecondHalf(true, true, true));
  coverage.addMeasurement(true, positive);
  coverage.addMeasurement(true, positive);
  require(!coverage.skipSecondHalf(true, true, true));
  coverage.addMeasurement(false, negative);
  coverage.addMeasurement(true, DetId(DetId::Tracker, 1));
  coverage.addMeasurement(true, DetId(DetId::Muon, MuonSubdetId::DT));
  coverage.addMeasurement(true, DetId(DetId::Muon, MuonSubdetId::RPC));
  require(!coverage.skipSecondHalf(true, true, true));
  coverage.addMeasurement(true, negative);
  require(coverage.skipSecondHalf(true, true, true));
  require(!coverage.skipSecondHalf(true, false, true));
  require(!coverage.skipSecondHalf(true, true, false));
  require(!coverage.skipSecondHalf(true, false, false));

  cosmic_muon::CSCEndcapCoverage reversed;
  reversed.addMeasurement(true, negative);
  reversed.addMeasurement(true, negative);
  require(!reversed.skipSecondHalf(true, true, true));
  reversed.addMeasurement(true, positive);
  require(reversed.skipSecondHalf(true, true, true));
  return 0;
}
