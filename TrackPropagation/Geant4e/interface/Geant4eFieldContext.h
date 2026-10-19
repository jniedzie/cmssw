// Code style: human-code.
#ifndef TrackPropagation_Geant4e_Geant4eFieldContext_h
#define TrackPropagation_Geant4e_Geant4eFieldContext_h

#include "G4ChordFinder.hh"
#include "G4EquationOfMotion.hh"
#include "G4ErrorMag_UsualEqRhs.hh"
#include "G4ErrorPropagatorManager.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4VIntegrationDriver.hh"

namespace TrackPropagation {
  inline bool ensureGeant4eFieldContext(G4ErrorPropagatorManager& manager) {
    auto* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    auto const* field = fieldManager ? fieldManager->GetDetectorField() : nullptr;
    if (!dynamic_cast<G4MagneticField const*>(field)) return false;

    auto* chord = fieldManager->GetChordFinder();
    auto* driver = chord ? chord->GetIntegrationDriver() : nullptr;
    if (!driver) return false;
    auto* equation = driver->GetEquationOfMotion();

    // GeometryProducer can replace the field manager at a Run boundary after
    // Geant4e has initialized. Restore its normal backward-transport setup.
    if (!dynamic_cast<G4ErrorMag_UsualEqRhs*>(equation)) {
      if (!manager.InitFieldForBackwards()) return false;
      chord = fieldManager->GetChordFinder();
      driver = chord ? chord->GetIntegrationDriver() : nullptr;
      equation = driver ? driver->GetEquationOfMotion() : nullptr;
    }
    if (!dynamic_cast<G4ErrorMag_UsualEqRhs*>(equation)) return false;

    // InitFieldForBackwards reuses its equation across Runs. Its field must
    // follow the active manager, including changes of the EventSetup provider.
    if (equation->GetFieldObj() != field)
      equation->SetFieldObj(const_cast<G4Field*>(field));
    return true;
  }
}
#endif
