// Code style: AI-code. AI-generated opt-in preservation of native GEN kinematics.
#ifndef SimG4Core_Generators_RecordedPrimaryMomentum_h
#define SimG4Core_Generators_RecordedPrimaryMomentum_h
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace sim {
  inline double recordedPrimaryEnergyTolerance(double energy, double physicalEnergy, double sourceScaleGeV) {
    const double scale=std::max({1.,sourceScaleGeV,energy,physicalEnergy});
    const double energyRoundoff=32.*std::numeric_limits<double>::epsilon()*scale;
    // Native generator boosts/fragmentation carry an event-scale invariant
    // uncertainty in GeV^2. Propagate it through E_on^2-E^2, rather than
    // demanding that each soft descendant has only local energy roundoff.
    // E^2-px^2-py^2-pz^2 has four signed square terms. Each
    // component is bounded by the original event scale; retain all four
    // contributions rather than using a single-component invariant bound.
    const double invariantRoundoff=4.*energyRoundoff*scale;
    if (!std::isfinite(invariantRoundoff))
      throw std::runtime_error("Preassigned source numerical envelope overflow");
    const double denominator=energy+physicalEnergy;
    const double propagated=denominator>0. ? invariantRoundoff/denominator : energyRoundoff;
    // Never enlarge the existing physical four-momentum conservation scale.
    const double conservation=1.e-7*std::max({1.,energy,physicalEnergy});
    return std::min(conservation,std::max(energyRoundoff,propagated));
  }

  template <typename Momentum>
  void setRecordedPrimaryMomentum(G4PrimaryParticle& target, const Momentum& p, double generatedMass,
                                  double sourcePrecisionScaleGeV=1.) {
    const double energy=p.e(),px=p.px(),py=p.py(),pz=p.pz();
    if (!std::isfinite(energy) || !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(generatedMass) || generatedMass<0. || energy<generatedMass ||
        !std::isfinite(sourcePrecisionScaleGeV) || sourcePrecisionScaleGeV<0.)
      throw std::runtime_error("Preassigned source requires finite physical momentum and generated mass");
    const double physicalEnergy=std::hypot(std::hypot(px,py,pz),generatedMass);
    const double rounding=recordedPrimaryEnergyTolerance(energy,physicalEnergy,sourcePrecisionScaleGeV);
    if (!std::isfinite(physicalEnergy) || std::abs(physicalEnergy-energy)>rounding) {
      std::ostringstream error;
      error << std::setprecision(17) << "Preassigned source is off shell beyond its invariant precision: pdg="
            << target.GetPDGcode() << " p4_GeV=(" << px << "," << py << "," << pz << "," << energy
            << ") generated_mass_GeV=" << generatedMass << " energy_residual_GeV=" << physicalEnergy-energy
            << " allowed_GeV=" << rounding << " original_event_scale_GeV=" << sourcePrecisionScaleGeV;
      throw std::runtime_error(error.str());
    }
    target.Set4Momentum(px*GeV,py*GeV,pz*GeV,energy*GeV);
    target.SetMass(generatedMass*GeV);
    target.SetTotalEnergy(energy*GeV);
  }
}
#endif
