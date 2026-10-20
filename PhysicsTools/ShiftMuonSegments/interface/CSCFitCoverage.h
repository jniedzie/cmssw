// Code style: human-code.
#ifndef PhysicsTools_ShiftMuonSegments_CSCFitCoverage_h
#define PhysicsTools_ShiftMuonSegments_CSCFitCoverage_h

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

namespace shift {
  struct CSCFitCoverage {
    unsigned int hits = 0;
    unsigned int endcapMask = 0;
    unsigned int stationMask = 0;

    void addMeasurement(bool valid, DetId id) {
      if (!valid || id.det() != DetId::Muon || id.subdetId() != MuonSubdetId::CSC)
        return;
      CSCDetId const csc(id);
      int const endcap = csc.endcap();
      int const station = csc.station();
      if (endcap < 1 || endcap > 2 || station < 1 || station > 4)
        return;
      ++hits;
      endcapMask |= 1U << (endcap - 1);
      // Bits 0..3 are +z stations 1..4; bits 4..7 are -z stations 1..4.
      stationMask |= 1U << (4 * (endcap - 1) + station - 1);
    }
  };

  inline bool passesCSCFitCoverage(bool enabled, unsigned int inputMask, unsigned int fittedMask) {
    return !enabled || (inputMask & 3U) != 3U || (fittedMask & 3U) == 3U;
  }
}  // namespace shift

#endif
