// Code style: AI-code. AI-generated event-local source-lineage interface.
#ifndef SimG4Core_Generators_SourceAncestry_h
#define SimG4Core_Generators_SourceAncestry_h

class G4Track;

namespace sim::sourceAncestry {
  void reset(bool enabled, double externalBoundary);
  void remember(const G4Track* track);
  // Fail if a primary conversion is repeated within this GEN event.
  void markKinematicsRestored(const G4Track* track);
  bool enabled();
  bool isExternal(const G4Track* track);
}

#endif
