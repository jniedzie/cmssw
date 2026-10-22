// Code style: AI-code. AI-generated event-local source-lineage registry.
#include "SimG4Core/Generators/interface/SourceAncestry.h"
#include "G4Track.hh"

#include <cmath>
#include <map>
#include <stdexcept>

namespace {
  thread_local bool active = false;
  thread_local double boundary = 0.;
  thread_local std::map<int, double> rootBirthZ;
}

void sim::sourceAncestry::reset(bool enabled, double externalBoundary) {
  rootBirthZ.clear();
  active = enabled;
  boundary = externalBoundary;
  if (active && (!std::isfinite(boundary) || boundary <= 0.))
    throw std::runtime_error("Source ancestry requires a finite positive external boundary");
}

bool sim::sourceAncestry::enabled() { return active; }

void sim::sourceAncestry::remember(const G4Track* track) {
  if (!active)
    return;
  const int id = track->GetTrackID();
  if (id <= 0)
    throw std::runtime_error("Source ancestry requires a positive classified Geant4 track ID");
  if (rootBirthZ.count(id))
    return;
  if (track->GetParentID() == 0) {
    // At stacking time vertex metadata need not yet be initialized. Current
    // position is the birth position; retain it before this track propagates.
    const double z = track->GetPosition().z();
    if (!std::isfinite(z))
      throw std::runtime_error("Source ancestry requires a finite primary birth position");
    rootBirthZ.emplace(id, z);
  } else {
    const auto parent = rootBirthZ.find(track->GetParentID());
    if (parent == rootBirthZ.end())
      throw std::runtime_error("Source ancestry parent was not registered before its secondary");
    rootBirthZ.emplace(id, parent->second);
  }
}

bool sim::sourceAncestry::isExternal(const G4Track* track) {
  if (!active)
    return false;
  const auto found = rootBirthZ.find(track->GetTrackID());
  if (found == rootBirthZ.end())
    throw std::runtime_error("Source ancestry track was not registered before classification");
  return std::abs(found->second) >= boundary;
}
