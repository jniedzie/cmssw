// Code style: human-code.
#include "TrackPropagation/Geant4e/interface/Geant4eFieldContext.h"
#include "G4ChargeState.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4UniformMagField.hh"
#include <array>
#include <cmath>
#include <stdexcept>

namespace {
  void check(bool value, char const* message) {
    if (!value) throw std::runtime_error(message);
  }

  G4EquationOfMotion* equation(G4FieldManager& manager) {
    return manager.GetChordFinder()->GetIntegrationDriver()->GetEquationOfMotion();
  }

  std::array<double, 6> derivatives(G4EquationOfMotion& equation, G4ErrorMode mode) {
    G4ErrorPropagatorData::GetErrorPropagatorData()->SetMode(mode);
    equation.SetChargeMomentumMass(G4ChargeState(1., 0., 0., 0., 0.), 1000., 105.658);
    double state[8] = {0., 0., 0., 1000., 0., 0., 0., 0.};
    std::array<double, 6> result{};
    equation.RightHandSide(state, result.data());
    return result;
  }

  std::array<double, 4> tolerances(G4FieldManager const& manager) {
    return {manager.GetDeltaOneStep(), manager.GetDeltaIntersection(),
            manager.GetMinimumEpsilonStep(), manager.GetMaximumEpsilonStep()};
  }
}

int main() {
  auto* transport = G4TransportationManager::GetTransportationManager();
  auto* propagation = G4ErrorPropagatorManager::GetErrorPropagatorManager();
  G4UniformMagField firstField(G4ThreeVector(0., 0., 1.));
  G4FieldManager firstManager(&firstField);
  firstManager.GetChordFinder()->SetDeltaChord(.007);
  transport->SetFieldManager(&firstManager);
  check(!dynamic_cast<G4ErrorMag_UsualEqRhs*>(equation(firstManager)), "Initial equation unexpectedly backward aware");
  check(TrackPropagation::ensureGeant4eFieldContext(*propagation), "First field initialization failed");
  check(dynamic_cast<G4ErrorMag_UsualEqRhs*>(equation(firstManager)), "First equation is not backward aware");
  double const firstBackwardDeltaChord = firstManager.GetChordFinder()->GetDeltaChord();
  auto const forward = derivatives(*equation(firstManager), G4ErrorMode_PropForwards);
  auto const backward = derivatives(*equation(firstManager), G4ErrorMode_PropBackwards);
  check(std::abs(forward[4]) > 0., "Test field produced no curvature");
  check(std::abs(backward[4] + forward[4]) < 1.e-12, "Backward curvature did not reverse");

  // Simulate a new Run with a new manager and a different field provider.
  G4UniformMagField secondField(G4ThreeVector(0., 0., 2.));
  G4FieldManager secondManager(&secondField);
  secondManager.GetChordFinder()->SetDeltaChord(.041);
  secondManager.SetDeltaOneStep(.0031);
  secondManager.SetDeltaIntersection(.0007);
  secondManager.SetMinimumEpsilonStep(1.e-6);
  secondManager.SetMaximumEpsilonStep(2.e-5);
  auto const settings = tolerances(secondManager);
  transport->SetFieldManager(&secondManager);
  check(!dynamic_cast<G4ErrorMag_UsualEqRhs*>(equation(secondManager)), "Run transition did not replace equation");
  check(TrackPropagation::ensureGeant4eFieldContext(*propagation), "Run transition repair failed");
  check(dynamic_cast<G4ErrorMag_UsualEqRhs*>(equation(secondManager)), "Repaired equation is not backward aware");
  check(equation(secondManager)->GetFieldObj() == &secondField, "Repaired equation retained old field provider");
  check(tolerances(secondManager) == settings, "Field-manager tolerances changed");
  check(secondManager.GetChordFinder()->GetDeltaChord() == firstBackwardDeltaChord,
        "Run transition did not restore first-Run backward chord contract");
  auto const secondBackward = derivatives(*equation(secondManager), G4ErrorMode_PropBackwards);
  check(std::abs(secondBackward[4] - 2. * backward[4]) < 1.e-12, "Propagation did not use changed field value");
  auto* const repairedChord = secondManager.GetChordFinder();
  for (unsigned int i = 0; i < 100; ++i)
    check(TrackPropagation::ensureGeant4eFieldContext(*propagation), "Repeated field-context check failed");
  check(secondManager.GetChordFinder() == repairedChord, "Repeated checks rebuilt chord finder");
  check(derivatives(*equation(secondManager), G4ErrorMode_PropBackwards) == secondBackward,
        "Repeated checks changed propagated derivatives");
}
