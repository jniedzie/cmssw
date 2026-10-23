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
  void setRecordedPrimaryMomentum(G4PrimaryParticle& target, const Momentum& p, double generatedMass) {
    const double energy = p.e(), px = p.px(), py = p.py(), pz = p.pz();
    if (!std::isfinite(energy) || !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(generatedMass) || generatedMass < 0. || energy < generatedMass)
      throw std::runtime_error("Preassigned source requires finite physical momentum and generated mass");
    const double p2 = px*px + py*py + pz*pz;
    const double e2 = energy*energy;
    // These saved components already contain Lorentz-boost roundoff, not
    // only the operations below. Soft photons/neutrinos can have a tiny
    // negative norm. Use a double-precision input tolerance with a 1 GeV^2
    // absolute scale; the native decay-conservation test remains separate.
    // Keep the recorded mass rather than deriving it by cancellation.
    const double rounding = 32.*std::numeric_limits<double>::epsilon()*std::max(1.,e2+p2);
    if (!std::isfinite(p2) || !std::isfinite(e2) || e2-p2 < -rounding)
      throw std::runtime_error("Preassigned source has genuinely spacelike momentum");
    target.Set4Momentum(px*GeV, py*GeV, pz*GeV, energy*GeV);
    // Set4Momentum recomputes/clamps mass. Restore the native source authority
    // and kinetic energy explicitly before the dynamic-particle conversion.
    target.SetMass(generatedMass*GeV);
    target.SetTotalEnergy(energy*GeV);
  }
}
#endif
