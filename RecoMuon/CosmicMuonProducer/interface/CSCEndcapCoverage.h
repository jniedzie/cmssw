// Code style: human-code.
#ifndef RecoMuon_CosmicMuonProducer_CSCEndcapCoverage_h
#define RecoMuon_CosmicMuonProducer_CSCEndcapCoverage_h

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

namespace cosmic_muon {
  class CSCEndcapCoverage {
  public:
    void addMeasurement(bool valid, DetId id) {
      if (!valid || id.det() != DetId::Muon || id.subdetId() != MuonSubdetId::CSC)
        return;
      int const endcap = CSCDetId(id.rawId()).endcap();
      firstEndcap_ = firstEndcap_ || endcap == 1;
      secondEndcap_ = secondEndcap_ || endcap == 2;
    }

    bool skipSecondHalf(bool enabled, bool traversing, bool strictOneLeg) const {
      return enabled && traversing && strictOneLeg && firstEndcap_ && secondEndcap_;
    }

  private:
    bool firstEndcap_ = false;
    bool secondEndcap_ = false;
  };
}  // namespace cosmic_muon

#endif
