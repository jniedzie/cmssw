// Code style: AI-code. AI-generated exact-zero-flight preassigned-decay support.
#ifndef SimG4Core_Generators_PromptDecaySupport_h
#define SimG4Core_Generators_PromptDecaySupport_h
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

namespace sim::promptDecay {
  inline int pdg(const HepMC::GenParticle& p) { return p.pdg_id(); }
  inline int pdg(const HepMC3::GenParticle& p) { return p.pid(); }
  inline int barcode(const HepMC::GenParticle& p) { return p.barcode(); }
  inline int barcode(const HepMC3::GenParticle& p) { return p.id(); }
  inline unsigned incoming(const HepMC::GenParticle& p) { return p.end_vertex()->particles_in_size(); }
  inline unsigned incoming(const HepMC3::GenParticle& p) { return p.end_vertex()->particles_in().size(); }
  inline std::vector<const HepMC::GenParticle*> children(const HepMC::GenParticle& p) {
    return {p.end_vertex()->particles_out_const_begin(), p.end_vertex()->particles_out_const_end()};
  }
  inline std::vector<const HepMC3::GenParticle*> children(const HepMC3::GenParticle& p) {
    std::vector<const HepMC3::GenParticle*> result;
    for (const auto& daughter : p.end_vertex()->particles_out()) result.push_back(daughter.get());
    return result;
  }

  template <typename Particle>
  bool exactlyPrompt(const Particle& p) {
    if (!p.production_vertex() || !p.end_vertex()) return false;
    const auto& a = p.production_vertex()->position();
    const auto& b = p.end_vertex()->position();
    // There is no distance/time tolerance. Even one nonzero coordinate is flight.
    for (const auto values : {std::pair{a.x(), b.x()}, std::pair{a.y(), b.y()},
                              std::pair{a.z(), b.z()}, std::pair{a.t(), b.t()}})
      if (!std::isfinite(values.first) || !std::isfinite(values.second) || values.first != values.second)
        return false;
    return true;
  }

  template <typename Particle>
  [[noreturn]] void unsupported(const Particle& p, const char* reason) {
    throw std::runtime_error("Undefined preassigned particle PDG=" + std::to_string(pdg(p)) +
                             " barcode=" + std::to_string(barcode(p)) + ": " + reason);
  }

  template <typename Particle, typename Known>
  void appendFrontier(const Particle& p, Known known, std::set<const Particle*>& path,
                      std::vector<const Particle*>& result) {
    if (!path.insert(&p).second) unsupported(p, "cycle in prompt decay tree");
    if (known(pdg(p))) {
      result.push_back(&p);
    } else {
      if (p.status() != 2 || !p.production_vertex() || !p.end_vertex())
        unsupported(p, "unknown terminal particle or missing decay vertices");
      if (!exactlyPrompt(p)) unsupported(p, "nonzero or nonfinite flight cannot be collapsed");
      if (incoming(p) != 1) unsupported(p, "decay vertex does not have exactly one incoming particle");
      const auto daughters = children(p);
      if (daughters.empty()) unsupported(p, "unknown decay has no daughters");
      for (const auto* daughter : daughters) appendFrontier(*daughter, known, path, result);
    }
    path.erase(&p);
  }

  // Return original pointers, without changing the GEN graph, statuses or weights.
  // Geant4 already promotes undefined nodes; this validates that promotion loses
  // exactly zero flight/time. Supported nodes keep their own preassigned decay.
  template <typename Particle, typename Known>
  std::vector<const Particle*> frontier(const Particle& p, Known known) {
    std::set<const Particle*> path;
    std::vector<const Particle*> result;
    appendFrontier(p, known, path, result);
    return result;
  }
}
#endif
