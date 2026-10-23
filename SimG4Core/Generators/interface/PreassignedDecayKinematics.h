// Code style: AI-code. AI-generated opt-in preassigned-decay conversion repair.
#ifndef SimG4Core_Generators_PreassignedDecayKinematics_h
#define SimG4Core_Generators_PreassignedDecayKinematics_h
class G4Track;
namespace sim {
  // Call on an original primary before its first transport step only. Every
  // retained PrimaryParticle must carry the complete source four-momentum.
  // Preassigned products are rebuilt in each mother's own rest frame, so
  // native G4Decay can boost them using the parent's eventual transported state.
  void preservePreassignedDecayKinematics(G4Track& track);
  // For a native decay daughter before its first transport step, retain its
  // current momentum and restore only the recorded mass. Never restore its
  // original GEN momentum: the parent may already have slowed or scattered.
  void preservePreassignedDaughterMass(G4Track& track);
}
#endif
