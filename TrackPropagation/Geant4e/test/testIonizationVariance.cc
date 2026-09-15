#include "TrackPropagation/Geant4e/interface/IonizationVariance.h"
#include <cmath>
#include <stdexcept>
int main() {
  for (double energy : {5.,50.,500.}) for (double density : {.0012,7.87}) {
    auto const whole = TrackPropagation::ionizationVariance(energy,.105658,26./55.845,density,100.);
    for (double count : {10.,100.,10000.}) {
      auto const step = TrackPropagation::ionizationVariance(energy,.105658,26./55.845,density,100./count);
      if (std::abs(count*step.physical/whole.physical-1.)>1.e-12 ||
          !(step.physical >= step.legacy && step.legacy >= 0.))
        throw std::runtime_error("Ionization second moment must be additive across steps");
    }
  }
}
