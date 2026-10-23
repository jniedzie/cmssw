// Code style: AI-code. AI-generated opt-in preservation of native GEN kinematics.
#ifndef SimG4Core_Generators_RecordedPrimaryMomentum_h
#define SimG4Core_Generators_RecordedPrimaryMomentum_h
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>

namespace sim {
  template <typename Momentum>
  void setRecordedPrimaryMomentum(G4PrimaryParticle& target, const Momentum& p, double generatedMass,
                                  double sourcePrecisionScaleGeV=1.) {
    const double energy=p.e(),px=p.px(),py=p.py(),pz=p.pz();
    if (!std::isfinite(energy) || !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(generatedMass) || generatedMass<0. || energy<generatedMass ||
        !std::isfinite(sourcePrecisionScaleGeV) || sourcePrecisionScaleGeV<0.)
      throw std::runtime_error("Preassigned source requires finite physical momentum and generated mass");
    const double physicalEnergy=std::hypot(std::hypot(px,py,pz),generatedMass);
    // Saved soft daughters inherit absolute error from earlier boosts of much
    // more energetic parents. The caller supplies the largest original event
    // energy, not the soft daughter's own energy. Compare energy directly to
    // avoid an additional nearly cancelling invariant-mass subtraction.
    const double rounding=32.*std::numeric_limits<double>::epsilon()*
        std::max({1.,sourcePrecisionScaleGeV,energy,physicalEnergy});
    if (!std::isfinite(physicalEnergy) || std::abs(physicalEnergy-energy)>rounding)
      throw std::runtime_error("Preassigned source is off shell beyond its event-scale numerical precision");
    target.Set4Momentum(px*GeV,py*GeV,pz*GeV,energy*GeV);
    target.SetMass(generatedMass*GeV);
    target.SetTotalEnergy(energy*GeV);
  }
}
#endif
