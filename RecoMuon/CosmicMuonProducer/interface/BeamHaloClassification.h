// Code style: human-code.
#ifndef RecoMuon_CosmicMuonProducer_BeamHaloClassification_h
#define RecoMuon_CosmicMuonProducer_BeamHaloClassification_h

#include <cmath>

namespace cosmic_muon {
  inline bool isBeamHaloForSmoothing(double firstEta, double lastEta, bool useAbsoluteEta) {
    // Keep the historical branch exact unless the symmetric study is enabled.
    if (!useAbsoluteEta)
      return firstEta > 4.0 || lastEta > 4.0;
    return std::abs(firstEta) > 4.0 || std::abs(lastEta) > 4.0;
  }
}  // namespace cosmic_muon
#endif
