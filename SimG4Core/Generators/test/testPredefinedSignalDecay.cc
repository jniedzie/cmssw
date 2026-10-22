// Code style: AI-code. AI-generated regression using actual HepMC2/3 converters.
// This checks primary admission and predefined lifetimes, not material transport.
#include "SimG4Core/Generators/interface/Generator.h"
#include "SimG4Core/Generators/interface/Generator3.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HepMC/GenVertex.h"
#include "HepMC3/GenVertex.h"
#include "G4BosonConstructor.hh"
#include "G4Event.hh"
#include "G4LeptonConstructor.hh"
#include "G4LorentzVector.hh"
#include "G4MesonConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  constexpr double darkMass = .3;  // GeV
  constexpr double kaonCTau = 30.;  // mm, a fixed sampled diagnostic lifetime
  constexpr double darkCTau = 300.;  // mm, a separate sampled lifetime

  void require(bool condition, const std::string& message) {
    if (!condition)
      throw std::runtime_error(message);
  }

  void close(double actual, double expected, const std::string& message) {
    require(std::isfinite(actual) && std::abs(actual - expected) < 1.e-9 * std::max(1., std::abs(expected)),
            message);
  }

  edm::ParameterSet parameters(const std::vector<int>& decayIds, bool explicitOption = true) {
    edm::ParameterSet config;
    config.addParameter<bool>("ApplyPCuts", true);
    config.addParameter<bool>("ApplyPtransCut", false);
    config.addParameter<bool>("ApplyEtaCuts", true);
    config.addParameter<bool>("ApplyPhiCuts", false);
    config.addParameter<bool>("ApplyLumiMonitorCuts", false);
    config.addParameter<bool>("IsSlepton", false);
    config.addParameter<double>("MinPhiCut", -3.141592653589793);
    config.addParameter<double>("MaxPhiCut", 3.141592653589793);
    config.addParameter<double>("MinEtaCut", -5.5);
    config.addParameter<double>("MaxEtaCut", 5.5);
    config.addParameter<double>("MinPCut", .04);
    config.addParameter<double>("MaxPCut", 99999.);
    config.addParameter<double>("EtaCutForHector", 8.2);
    config.addParameter<double>("RDecLenCut", 2.9);
    config.addParameter<double>("LDecLenCut", 10.);
    config.addParameter<double>("MaxZCentralCMS", 12.);
    if (explicitOption)
      config.addUntrackedParameter<std::vector<int>>("AlwaysTrackDecaysTo", decayIds);
    return config;
  }

  struct Chain {
    int parentId, recoilId;
    G4LorentzVector parent, recoil, dark, positiveMuon, negativeMuon;
    G4LorentzVector production, kaonDecay, darkDecay;
  };

  Chain chain(int parentId, double startZ, bool signal) {
    Chain result;
    result.parentId = parentId;
    result.recoilId = std::abs(parentId) == 321 ? (parentId > 0 ? 211 : -211) : 22;
    double mass = G4ParticleTable::GetParticleTable()->FindParticle(parentId)->GetPDGMass() / GeV;
    double recoilMass = G4ParticleTable::GetParticleTable()->FindParticle(result.recoilId)->GetPDGMass() / GeV;
    result.parent = G4LorentzVector(0., 0., 1., std::hypot(mass, 1.));
    double finalMass = signal ? darkMass : .25;  // SM-only virtual-photon control
    double eDark = (mass * mass + finalMass * finalMass - recoilMass * recoilMass) / (2. * mass);
    double momentum = std::sqrt(eDark * eDark - finalMass * finalMass);
    result.dark = G4LorentzVector(momentum, 0., 0., eDark);
    result.recoil = G4LorentzVector(-momentum, 0., 0., mass - eDark);
    result.dark.boost(result.parent.boostVector());
    result.recoil.boost(result.parent.boostVector());
    double muonMass = G4ParticleTable::GetParticleTable()->FindParticle(13)->GetPDGMass() / GeV;
    double muonMomentum = std::sqrt(finalMass * finalMass / 4. - muonMass * muonMass);
    result.positiveMuon = G4LorentzVector(0., muonMomentum, 0., finalMass / 2.);
    result.negativeMuon = G4LorentzVector(0., -muonMomentum, 0., finalMass / 2.);
    result.positiveMuon.boost(result.dark.boostVector());
    result.negativeMuon.boost(result.dark.boostVector());
    result.production = G4LorentzVector(0., 0., startZ, 0.);
    result.kaonDecay = result.production + result.parent * (kaonCTau / mass);
    result.darkDecay = result.kaonDecay + result.dark * (darkCTau / finalMass);
    return result;
  }

  HepMC::FourVector hep2(const G4LorentzVector& vector) {
    return HepMC::FourVector(vector.x(), vector.y(), vector.z(), vector.t());
  }

  HepMC3::FourVector hep3(const G4LorentzVector& vector) {
    return HepMC3::FourVector(vector.x(), vector.y(), vector.z(), vector.t());
  }

  std::unique_ptr<HepMC::GenEvent> event2(const Chain& source, bool signal) {
    auto event = std::make_unique<HepMC::GenEvent>(HepMC::Units::GEV, HepMC::Units::MM);
    auto* origin = new HepMC::GenVertex(hep2(source.production));
    auto* kaonDecay = new HepMC::GenVertex(hep2(source.kaonDecay));
    auto* finalDecay = new HepMC::GenVertex(hep2(source.darkDecay));
    auto* kaon = new HepMC::GenParticle(hep2(source.parent), source.parentId, 2);
    auto* recoil = new HepMC::GenParticle(hep2(source.recoil), source.recoilId, 1);
    auto* intermediate = new HepMC::GenParticle(hep2(source.dark), signal ? 32 : 22, signal ? 2 : 3);
    origin->add_particle_out(kaon);
    kaonDecay->add_particle_in(kaon);
    kaonDecay->add_particle_out(recoil);
    kaonDecay->add_particle_out(intermediate);
    finalDecay->add_particle_in(intermediate);
    finalDecay->add_particle_out(new HepMC::GenParticle(hep2(source.positiveMuon), -13, 1));
    finalDecay->add_particle_out(new HepMC::GenParticle(hep2(source.negativeMuon), 13, 1));
    event->add_vertex(origin);
    event->add_vertex(kaonDecay);
    event->add_vertex(finalDecay);
    event->weights().push_back(1.);
    return event;
  }

  std::unique_ptr<HepMC3::GenEvent> event3(const Chain& source, bool signal) {
    auto event = std::make_unique<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);
    auto origin = std::make_shared<HepMC3::GenVertex>(hep3(source.production));
    auto kaonDecay = std::make_shared<HepMC3::GenVertex>(hep3(source.kaonDecay));
    auto finalDecay = std::make_shared<HepMC3::GenVertex>(hep3(source.darkDecay));
    auto kaon = std::make_shared<HepMC3::GenParticle>(hep3(source.parent), source.parentId, 2);
    auto recoil = std::make_shared<HepMC3::GenParticle>(hep3(source.recoil), source.recoilId, 1);
    auto intermediate = std::make_shared<HepMC3::GenParticle>(hep3(source.dark), signal ? 32 : 22, signal ? 2 : 3);
    origin->add_particle_out(kaon);
    kaonDecay->add_particle_in(kaon);
    kaonDecay->add_particle_out(recoil);
    kaonDecay->add_particle_out(intermediate);
    finalDecay->add_particle_in(intermediate);
    finalDecay->add_particle_out(std::make_shared<HepMC3::GenParticle>(hep3(source.positiveMuon), -13, 1));
    finalDecay->add_particle_out(std::make_shared<HepMC3::GenParticle>(hep3(source.negativeMuon), 13, 1));
    event->add_vertex(origin);
    event->add_vertex(kaonDecay);
    event->add_vertex(finalDecay);
    event->weights().push_back(1.);
    return event;
  }

  void particleSignature(std::ostream& output, const G4PrimaryParticle* particle) {
    output << particle->GetPDGcode() << ':' << particle->GetPx() << ',' << particle->GetPy() << ','
           << particle->GetPz() << ',' << particle->GetMass() << ',' << particle->GetProperTime() << '[';
    for (auto* child = particle->GetDaughter(); child != nullptr; child = child->GetNext())
      particleSignature(output, child);
    output << ']';
  }

  std::vector<std::string> signature(const G4Event& event) {
    std::vector<std::string> rows;
    for (int index = 0; index < event.GetNumberOfPrimaryVertex(); ++index) {
      auto* vertex = event.GetPrimaryVertex(index);
      for (auto* particle = vertex->GetPrimary(); particle != nullptr; particle = particle->GetNext()) {
        std::ostringstream row;
        row << std::setprecision(14) << vertex->GetX0() << ',' << vertex->GetY0() << ',' << vertex->GetZ0() << ','
            << vertex->GetT0() << ':';
        particleSignature(row, particle);
        rows.push_back(row.str());
      }
    }
    std::sort(rows.begin(), rows.end());
    return rows;
  }

  std::unique_ptr<G4Event> convert(const Chain& source, bool signal, bool hepMC3, const edm::ParameterSet& config) {
    auto result = std::make_unique<G4Event>(1);
    if (hepMC3) {
      auto input = event3(source, signal);
      Generator3 converter(config);
      converter.HepMC2G4(input.get(), result.get());
      require(input->particles().front()->status() == 2 && input->weights().at(0) == 1., "HepMC3 input mutated");
    } else {
      auto input = event2(source, signal);
      Generator converter(config);
      converter.HepMC2G4(input.get(), result.get());
      require(input->weights()[0] == 1., "HepMC2 input weight mutated");
      for (auto particle = input->particles_begin(); particle != input->particles_end(); ++particle)
        require((*particle)->status() < 1000, "HepMC2 input status mutated");
    }
    return result;
  }

  void checkNested(const G4Event& event, int parentId, double startZ) {
    const G4PrimaryParticle* parent = nullptr;
    int roots = 0;
    for (int index = 0; index < event.GetNumberOfPrimaryVertex(); ++index) {
      auto* vertex = event.GetPrimaryVertex(index);
      for (auto* particle = vertex->GetPrimary(); particle != nullptr; particle = particle->GetNext()) {
        ++roots;
        parent = particle;
        close(vertex->GetZ0() / mm, startZ, "Parent not injected at production vertex");
      }
    }
    require(roots == 1 && parent != nullptr && parent->GetPDGcode() == parentId, "Missing/duplicate parent roots");
    close(parent->GetProperTime() / ns, kaonCTau / (c_light / (mm / ns)), "Kaon lifetime lost");
    int children = 0, darkCount = 0;
    for (auto* child = parent->GetDaughter(); child != nullptr; child = child->GetNext()) {
      ++children;
      if (child->GetPDGcode() != 32)
        continue;
      ++darkCount;
      close(child->GetProperTime() / ns, darkCTau / (c_light / (mm / ns)), "Dark-photon lifetime lost");
      int muons = 0;
      for (auto* muon = child->GetDaughter(); muon != nullptr; muon = muon->GetNext()) {
        require(std::abs(muon->GetPDGcode()) == 13, "Wrong dark-photon daughter");
        ++muons;
      }
      require(muons == 2, "Missing/duplicate dimuon daughters");
    }
    require(children == 2 && darkCount == 1, "Missing/duplicate nested dark-photon chain");
  }
}  // namespace

int main() {
  try {
    G4BosonConstructor().ConstructParticle();
    G4LeptonConstructor().ConstructParticle();
    G4MesonConstructor().ConstructParticle();
    G4ParticleTable::GetParticleTable()->SetReadiness();
    new G4ParticleDefinition("test_dark_photon", darkMass * GeV, 0., 0., 2, -1, -1, 0, 0, 0,
                             "boson", 0, 0, 32, false, 1. * ns, nullptr, false, "test", 32);
    int cases = 0;
    for (int parentId : {130, 310, 321, -321}) {
      for (double startZ : {-148000., 148000.}) {
        for (bool useHepMC3 : {false, true}) {
          auto source = chain(parentId, startZ, true);
          auto omitted = convert(source, true, useHepMC3, parameters({}, false));
          auto empty = convert(source, true, useHepMC3, parameters({}));
          auto unrelated = convert(source, true, useHepMC3, parameters({999}));
          require(signature(*omitted) == signature(*empty), "Empty option changed default primary tree");
          require(signature(*omitted) == signature(*unrelated), "Unrelated tag changed primary tree");
          auto tracked = convert(source, true, useHepMC3, parameters({32}));
          checkNested(*tracked, parentId, startZ);
          auto smSource = chain(parentId, startZ, false);
          auto smDefault = convert(smSource, false, useHepMC3, parameters({}));
          auto smSelected = convert(smSource, false, useHepMC3, parameters({32}));
          require(signature(*smDefault) == signature(*smSelected), "Unrelated SM control changed");
          ++cases;
        }
      }
    }
    for (int invalidId : {0, -32}) {
      bool rejected2 = false, rejected3 = false;
      try { Generator converter(parameters({invalidId})); } catch (...) { rejected2 = true; }
      try { Generator3 converter(parameters({invalidId})); } catch (...) { rejected3 = true; }
      require(rejected2 && rejected3, "Nonpositive absolute-PDG selector accepted");
    }
    std::cout << "Predefined signal decay: " << cases << " HepMC2/3 graph controls passed; no material transport claim.\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Predefined signal decay regression failed: " << error.what() << '\n';
    return 1;
  }
}
