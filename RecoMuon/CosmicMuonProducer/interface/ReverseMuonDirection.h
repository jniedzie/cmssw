// Code style: human-code.
#ifndef RecoMuon_CosmicMuonProducer_ReverseMuonDirection_h
#define RecoMuon_CosmicMuonProducer_ReverseMuonDirection_h

#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonUtilities.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"

namespace cosmic_muon {
  inline void reverseMuonDirection(TrajectoryStateOnSurface& state,
                                   MagneticField const* field,
                                   CosmicMuonUtilities const& utilities,
                                   bool useConsistentCovariance) {
    if (!useConsistentCovariance) {
      utilities.reverseDirection(state, field);
      return;
    }
    auto covariance = state.cartesianError().matrix();
    // x is unchanged and p changes sign, so the x-p covariance blocks do too.
    for (unsigned int row = 0; row < 6; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        if ((row < 3) != (column < 3))
          covariance(row, column) = -covariance(row, column);
    GlobalTrajectoryParameters const parameters(state.globalPosition(), -state.globalMomentum(), -state.charge(), field);
    state = TrajectoryStateOnSurface(parameters, CartesianTrajectoryError(covariance), state.surface());
  }
}  // namespace cosmic_muon
#endif
