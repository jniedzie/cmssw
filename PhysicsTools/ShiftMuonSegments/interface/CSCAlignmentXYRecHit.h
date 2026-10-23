// Code style: human-code.
// CSC point measurements use the positional alignment-error block.
#ifndef PhysicsTools_ShiftMuonSegments_CSCAlignmentXYRecHit_h
#define PhysicsTools_ShiftMuonSegments_CSCAlignmentXYRecHit_h

#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/AlignmentPositionError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "Geometry/CommonTopologies/interface/GeomDet.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/GenericTransientTrackingRecHit.h"
#include <typeinfo>

namespace shift::alignment_xy {
  class CSCAlignmentXYRecHit final : public GenericTransientTrackingRecHit {
  public:
    CSCAlignmentXYRecHit(GeomDet const& geometry, TrackingRecHit const& persistent)
        : GenericTransientTrackingRecHit(geometry, persistent) {}

    AlgebraicSymMatrix parametersError() const override {
      AlgebraicSymMatrix error = GenericTransientTrackingRecHit::parametersError();
      auto const* alignment = det()->alignmentPositionError();
      if (!alignment)
        return error;

      // Keep the original MuonTransientTrackingRecHit pose convention exactly.
      // Only select its positional (x,y) block for the CSC point measurement.
      AlgebraicVector positions(2, 0), directions(2, 0);
      positions[0] = localPosition().x();
      auto const projected = ErrorFrameTransformer::transform46(alignment->globalError(), positions, directions);
      error[0][0] += projected.cxx();
      error[0][1] += projected.cyx();
      error[1][1] += projected.cyy();
      return error;
    }

    ConstRecHitContainer transientHits() const override {
      return {ConstRecHitPointer(clone())};
    }

  private:
    CSCAlignmentXYRecHit* clone() const override {
      return new CSCAlignmentXYRecHit(*this);
    }
  };

  inline TransientTrackingRecHit::RecHitPointer correctedHit(
      TransientTrackingRecHit::RecHitPointer const& original, bool enabled) {
    // The disabled path keeps the exact original object and all dispatch.
    if (!enabled || !original || !original->isValid() || !original->det())
      return original;
    auto const* persistent = original->hit();
    if (!dynamic_cast<MuonTransientTrackingRecHit const*>(original.get()) || !persistent ||
        typeid(*persistent) != typeid(CSCRecHit2D) || persistent->dimension() != 2 ||
        !persistent->recHits().empty())
      return original;
    return TransientTrackingRecHit::RecHitPointer(
        new CSCAlignmentXYRecHit(*original->det(), *persistent));
  }
}
#endif
