// Code style: human-code.
#include "SimG4Core/Generators/interface/Generator.h"
#include "SimG4Core/Generators/interface/HepMCParticle.h"
#include "SimG4Core/Generators/interface/LumiMonitorFilter.h"
#include "SimG4Core/Generators/interface/SourceAncestry.h"
#include "SimG4Core/Generators/interface/PromptDecaySupport.h"
#include "SimG4Core/Generators/interface/RecordedPrimaryMomentum.h"

#include "SimG4Core/Notification/interface/SimG4Exception.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "HepPDT/ParticleID.hh"

#include "G4Event.hh"
#include "G4HEPEvtParticle.hh"
#include "G4Log.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4UnitsTable.hh"

#include <algorithm>
#include <cmath>
#include <set>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace edm;

Generator::Generator(const ParameterSet &p)
    : fPCuts(p.getParameter<bool>("ApplyPCuts")),
      fPtransCut(p.getParameter<bool>("ApplyPtransCut")),
      fEtaCuts(p.getParameter<bool>("ApplyEtaCuts")),
      fPhiCuts(p.getParameter<bool>("ApplyPhiCuts")),
      fSlepton(p.getParameter<bool>("IsSlepton")),
      theMinPhiCut(p.getParameter<double>("MinPhiCut")),  // in radians (CMS standard)
      theMaxPhiCut(p.getParameter<double>("MaxPhiCut")),
      theMinEtaCut(p.getParameter<double>("MinEtaCut")),
      theMaxEtaCut(p.getParameter<double>("MaxEtaCut")),
      theMinPCut(p.getParameter<double>("MinPCut")),  // in GeV (CMS standard)
      theMaxPCut(p.getParameter<double>("MaxPCut")),
      theEtaCutForHector(p.getParameter<double>("EtaCutForHector")),
      verbose(p.getUntrackedParameter<int>("Verbosity", 0)),
      fLumiFilter(nullptr),
      evt_(nullptr),
      vtx_(nullptr),
      weight_(0),
      Z_lmin(0),
      Z_lmax(0),
      Z_hector(0),
      pdgFilterSel(false),
      fPDGFilter(false),
      debugMuonPrimaries_(p.getUntrackedParameter<bool>("DebugMuonPrimaries", false)),
      alwaysTrackDecaysTo_(p.getUntrackedParameter<std::vector<int>>("AlwaysTrackDecaysTo", {})),
      trackAllHadronDecays_(p.getUntrackedParameter<bool>("TrackAllHadronDecays", false)),
      maxUntrackedHadronFlightMm_(p.getUntrackedParameter<double>("MaxUntrackedHadronFlightMm", 0.)) {
  if (!std::isfinite(maxUntrackedHadronFlightMm_) || maxUntrackedHadronFlightMm_ < 0. ||
      maxUntrackedHadronFlightMm_ > 1.e-5 || (!trackAllHadronDecays_ && maxUntrackedHadronFlightMm_ != 0.))
    throw SimG4Exception("Tiny-flight approximation requires complete transport and a budget between zero and 10 nm");
  for (int pdg : alwaysTrackDecaysTo_) {
    if (pdg <= 0)
      throw SimG4Exception("AlwaysTrackDecaysTo requires positive absolute PDG IDs");
  }
  bool lumi = p.getParameter<bool>("ApplyLumiMonitorCuts");
  if (lumi) {
    fLumiFilter = new LumiMonitorFilter();
  }

  double theRDecLenCut = p.getParameter<double>("RDecLenCut") * CLHEP::cm;
  theDecRCut2 = theRDecLenCut * theRDecLenCut;

  theMinPtCut2 = theMinPCut * theMinPCut;

  double theDecLenCut = p.getParameter<double>("LDecLenCut") * CLHEP::cm;

  maxZCentralCMS = p.getParameter<double>("MaxZCentralCMS") * CLHEP::m;

  fFiductialCuts = (fPCuts || fPtransCut || fEtaCuts || fPhiCuts);

  pdgFilter.resize(0);
  if (p.exists("PDGselection")) {
    pdgFilterSel = (p.getParameter<edm::ParameterSet>("PDGselection")).getParameter<bool>("PDGfilterSel");
    pdgFilter = (p.getParameter<edm::ParameterSet>("PDGselection")).getParameter<std::vector<int>>("PDGfilter");
    if (!pdgFilter.empty()) {
      fPDGFilter = true;
      std::stringstream ss;
      ss << "SimG4Core/Generator: ";
      if (pdgFilterSel) {
        ss << " Selecting only PDG ID = ";
      } else {
        ss << " Filtering out PDG ID = ";
      }
      for (unsigned int ii = 0; ii < pdgFilter.size(); ++ii) {
        ss << pdgFilter[ii] << "  ";
      }
      edm::LogVerbatim("SimG4CoreGenerator") << ss.str();
    }
  }

  if (fEtaCuts) {
    Z_lmax = theRDecLenCut * ((1 - exp(-2 * theMaxEtaCut)) / (2 * exp(-theMaxEtaCut)));
    Z_lmin = theRDecLenCut * ((1 - exp(-2 * theMinEtaCut)) / (2 * exp(-theMinEtaCut)));
  }

  Z_hector = theRDecLenCut * ((1 - exp(-2 * theEtaCutForHector)) / (2 * exp(-theEtaCutForHector)));

  edm::LogVerbatim("SimG4CoreGenerator") << "SimG4Core/Generator: Rdecaycut= " << theRDecLenCut / CLHEP::cm
                                         << " cm;  Zdecaycut= " << theDecLenCut / CLHEP::cm
                                         << "Z_min= " << Z_lmin / CLHEP::cm << " cm; Z_max= " << Z_lmax / CLHEP::cm
                                         << " cm;\n"
                                         << "                     MaxZCentralCMS = " << maxZCentralCMS / CLHEP::m
                                         << " m;"
                                         << " Z_hector = " << Z_hector / CLHEP::cm << " cm\n"
                                         << "                     ApplyCuts: " << fFiductialCuts
                                         << "  PCuts: " << fPCuts << "  PtransCut: " << fPtransCut
                                         << "  EtaCut: " << fEtaCuts << "  PhiCut: " << fPhiCuts
                                         << "  LumiMonitorCut: " << lumi;
  if (fFiductialCuts) {
    edm::LogVerbatim("SimG4CoreGenerator")
        << "SimG4Core/Generator: Pmin(GeV)= " << theMinPCut << "; Pmax(GeV)= " << theMaxPCut
        << "; EtaMin= " << theMinEtaCut << "; EtaMax= " << theMaxEtaCut << "; PhiMin(rad)= " << theMinPhiCut
        << "; PhiMax(rad)= " << theMaxPhiCut;
  }
  if (lumi) {
    fLumiFilter->Describe();
  }
}

Generator::~Generator() { delete fLumiFilter; }

bool Generator::hasTrackedDecayDaughter(const HepMC::GenParticle *particle) const {
  if (particle->end_vertex() == nullptr || alwaysTrackDecaysTo_.empty())
    return false;
  for (auto daughter = particle->end_vertex()->particles_out_const_begin();
       daughter != particle->end_vertex()->particles_out_const_end(); ++daughter) {
    if (std::find(alwaysTrackDecaysTo_.begin(), alwaysTrackDecaysTo_.end(), std::abs((*daughter)->pdg_id())) !=
        alwaysTrackDecaysTo_.end())
      return true;
  }
  return false;
}

bool Generator::hasFiniteHadronFlight(const HepMC::GenParticle *particle) const {
  if (!trackAllHadronDecays_ || particle->status() != 2 || !particle->production_vertex() ||
      !particle->end_vertex() || !HepPDT::ParticleID(particle->pdg_id()).isHadron())
    return false;
  if (fPDGFilter) {
    const bool listed = IsInTheFilterList(particle->pdg_id());
    if ((!pdgFilterSel && listed) || (pdgFilterSel && !listed))
      return false;
  }
  const auto &birth = particle->production_vertex()->position();
  const auto &decay = particle->end_vertex()->position();
  const double length = std::hypot(decay.x() - birth.x(), decay.y() - birth.y(), decay.z() - birth.z());
  const double deltaCt = decay.t() - birth.t();
  const auto &p4 = particle->momentum();
  const double momentum = std::hypot(p4.px(), p4.py(), p4.pz());
  if (!std::isfinite(length) || !std::isfinite(deltaCt) || deltaCt < 0. ||
      !std::isfinite(momentum) || !std::isfinite(p4.e()) || p4.e() <= 0.)
    throw SimG4Exception("TrackAllHadronDecays requires finite forward-time hadron flight and momentum");
  // At rest a finite positive clock is a real lifetime, not prompt history.
  if (length == 0. && deltaCt == 0.)
    return false;
  if ((length == 0.) != (momentum == 0.))
    throw SimG4Exception("TrackAllHadronDecays found inconsistent spatial flight and parent momentum");
  if (particle->end_vertex()->particles_in_size() != 1)
    throw SimG4Exception("TrackAllHadronDecays requires a single parent at a finite decay vertex");
  if (G4ParticleTable::GetParticleTable()->FindParticle(particle->pdg_id()) == nullptr) {
    if (sim::promptDecay::withinPromptBudget(*particle,maxUntrackedHadronFlightMm_)) return false;
    throw SimG4Exception("TrackAllHadronDecays cannot transport an undefined Geant4 hadron beyond its prompt budget");
  }
  return true;
}

void Generator::HepMC2G4(const HepMC::GenEvent *evt_orig, G4Event *g4evt) {
  sim::sourceAncestry::reset(trackAllHadronDecays_, maxZCentralCMS);
  HepMC::GenEvent *evt = new HepMC::GenEvent(*evt_orig);

  if (*(evt->vertices_begin()) == nullptr) {
    std::stringstream ss;
    ss << "SimG4Core/Generator: in event " << g4evt->GetEventID() << " Corrupted Event - GenEvent with no vertex \n";
    throw SimG4Exception(ss.str());
  }

  if (!evt->weights().empty()) {
    weight_ = evt->weights()[0];
    for (unsigned int iw = 1; iw < evt->weights().size(); ++iw) {
      // terminate if the vector of weights contains a zero-weight
      if (evt->weights()[iw] <= 0)
        break;
      weight_ *= evt->weights()[iw];
    }
  }

  if (vtx_ != nullptr) {
    delete vtx_;
  }
  vtx_ = new math::XYZTLorentzVector((*(evt->vertices_begin()))->position().x(),
                                     (*(evt->vertices_begin()))->position().y(),
                                     (*(evt->vertices_begin()))->position().z(),
                                     (*(evt->vertices_begin()))->position().t());

  edm::LogVerbatim("SimG4CoreGenerator") << "Generator: primary Vertex = (" << vtx_->x() << ", " << vtx_->y() << ", "
                                         << vtx_->z() << ")";

  unsigned int ng4vtx = 0;
  unsigned int ng4par = 0;
  unsigned int nHepMCMuons = 0;
  unsigned int nG4MuonPrimaries = 0;
  std::map<std::string, unsigned int> muonFirstFailure;

  // Attach each finite hadron tree once, independent of HepMC vertex order.
  std::set<const HepMC::GenParticle *> transportedDescendants;
  for (auto parent = evt->particles_begin(); parent != evt->particles_end(); ++parent) {
    if (trackAllHadronDecays_ && maxUntrackedHadronFlightMm_ > 0. && (*parent)->status()==2 &&
        G4ParticleTable::GetParticleTable()->FindParticle((*parent)->pdg_id()) == nullptr &&
        !sim::promptDecay::exactlyPrompt(*(*parent)) &&
        sim::promptDecay::withinPromptBudget(*(*parent),maxUntrackedHadronFlightMm_)) {
      const auto& a=(*parent)->production_vertex()->position();
      const auto& z=(*parent)->end_vertex()->position();
      std::cout << "[SHIFT tiny-flight approximation] event=" << g4evt->GetEventID()
                << " barcode=" << (*parent)->barcode() << " pdg=" << (*parent)->pdg_id()
                << " flight_mm=" << std::hypot(z.x()-a.x(),z.y()-a.y(),z.z()-a.z())
                << " delta_ct_mm=" << z.t()-a.t() << " maximum_mm=" << maxUntrackedHadronFlightMm_ << std::endl;
    }
    if (!hasFiniteHadronFlight(*parent))
      continue;
    std::vector<HepMC::GenParticle *> pending;
    for (auto child = (*parent)->end_vertex()->particles_begin(HepMC::children);
         child != (*parent)->end_vertex()->particles_end(HepMC::children); ++child)
      pending.push_back(*child);
    while (!pending.empty()) {
      auto *child = pending.back();
      pending.pop_back();
      if (child == *parent)
        throw SimG4Exception("TrackAllHadronDecays found a cyclic decay graph");
      if (!transportedDescendants.insert(child).second)
        continue;
      if (G4ParticleTable::GetParticleTable()->FindParticle(child->pdg_id()) == nullptr) {
        // Native Geant4 promotes undefined prompt nodes. Validate that this
        // loses exactly zero source flight/time; finite unknowns still fail.
        sim::promptDecay::frontier(*child, [](int id) {
          return G4ParticleTable::GetParticleTable()->FindParticle(id) != nullptr;
        }, maxUntrackedHadronFlightMm_);
      }
      if (auto *end = child->end_vertex()) {
        if (child->status() != 2 && !(child->status() > 50 && child->status() < 100))
          throw SimG4Exception("TrackAllHadronDecays found an unsupported preassigned daughter status");
        for (auto next = end->particles_begin(HepMC::children); next != end->particles_end(HepMC::children); ++next)
          pending.push_back(*next);
      }
    }
  }


  for (HepMC::GenEvent::vertex_const_iterator vitr = evt->vertices_begin(); vitr != evt->vertices_end(); ++vitr) {
    // loop for vertex, is it a real vertex?
    // Set qvtx to true for any particles that should be propagated by GEANT,
    // i.e., status 1 particles or status 2 particles that decay outside the
    // beampipe.
    G4bool qvtx = false;
    HepMC::GenVertex::particle_iterator pitr;
    for (pitr = (*vitr)->particles_begin(HepMC::children); pitr != (*vitr)->particles_end(HepMC::children); ++pitr) {
      if (transportedDescendants.count(*pitr))
        continue;
      // For purposes of this function, the status is defined as follows:
      // 1:  particles are not decayed by generator
      // 2:  particles are decayed by generator but need to be propagated by GEANT
      // 3:  particles are decayed by generator and do not need to be propagated by GEANT
      int status = (*pitr)->status();
      int pdg = (*pitr)->pdg_id();
      if (debugMuonPrimaries_ && std::abs(pdg) == 13 && status == 1) {
        ++nHepMCMuons;
        const auto &p4 = (*pitr)->momentum();
        const auto &xyz = (*vitr)->position();
        std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                  << " stage=inventory decision=seen barcode=" << (*pitr)->barcode() << " pdgId=" << pdg
                  << " status=" << status << " px=" << p4.px() << " py=" << p4.py() << " pz=" << p4.pz()
                  << " E=" << p4.e() << " vertex=(" << xyz.x() << "," << xyz.y() << "," << xyz.z() << ")"
                  << std::endl;
      }
      if (status > 3 && isExotic(pdg) && (!(isExoticNonDetectable(pdg)))) {
        // In Pythia 8, there are many status codes besides 1, 2, 3.
        // By setting the status to 2 for exotic particles, they will be
        // checked: if its decay vertex is outside the beampipe, it will be
        // propagated by GEANT. Some Standard Model particles, e.g., K0, cannot
        // be propagated by GEANT, so do not change their status code.
        status = 2;
      }
      if (status == 2 && abs(pdg) == 9900015) {
        status = 3;
      }

      // Particles which are not decayed by generator
      if (status == 1) {
        // filter out unwanted particles and vertices
        if (fPDGFilter && !pdgFilterSel && IsInTheFilterList(pdg)) {
          continue;
        }

        qvtx = true;
        if (verbose > 2)
          LogDebug("SimG4CoreGenerator") << "GenVertex barcode = " << (*vitr)->barcode()
                                         << " selected for GenParticle barcode = " << (*pitr)->barcode();
        break;
      }
      // The selection is made considering if the partcile with status = 2
      // have the end_vertex with a radius greater than the radius of beampipe
      // cylinder (no requirement on the Z of the vertex is applyed).
      else if (status == 2) {
        // A selected signal's parent must traverse material before decaying,
        // even when its generator decay vertex lies outside the CMS cylinder.
        if (hasTrackedDecayDaughter(*pitr) || hasFiniteHadronFlight(*pitr)) {
          qvtx = true;
          break;
        }
        if ((*pitr)->end_vertex() != nullptr) {
          double xx = (*pitr)->end_vertex()->position().x();
          double yy = (*pitr)->end_vertex()->position().y();
          double r_dd = xx * xx + yy * yy;
          if (r_dd > theDecRCut2) {
            qvtx = true;
            if (verbose > 2)
              LogDebug("SimG4CoreGenerator")
                  << "GenVertex barcode = " << (*vitr)->barcode()
                  << " selected for GenParticle barcode = " << (*pitr)->barcode() << " radius = " << std::sqrt(r_dd);
            break;
          }
        } else {
          // particles with status 2 without end_vertex are
          // equivalent to stable
          qvtx = true;
          break;
        }
      }
    }

    // if this vertex is inside fiductial volume inside the beam pipe
    // and has no long-lived secondary the vertex is not saved
    if (!qvtx) {
      if (debugMuonPrimaries_) {
        for (pitr = (*vitr)->particles_begin(HepMC::children); pitr != (*vitr)->particles_end(HepMC::children); ++pitr) {
          if (std::abs((*pitr)->pdg_id()) == 13 && (*pitr)->status() == 1)
            ++muonFirstFailure["vertex_not_selected"];
          if (std::abs((*pitr)->pdg_id()) == 13 && (*pitr)->status() == 1)
            std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                      << " stage=vertex-selection decision=rejected barcode=" << (*pitr)->barcode()
                      << " reason=no_status1_or_long_lived_status2_daughter_at_vertex" << std::endl;
        }
      }
      continue;
    }

    double x1 = (*vitr)->position().x() * CLHEP::mm;
    double y1 = (*vitr)->position().y() * CLHEP::mm;
    double z1 = (*vitr)->position().z() * CLHEP::mm;
    double t1 = (*vitr)->position().t() * CLHEP::mm / CLHEP::c_light;

    G4PrimaryVertex *g4vtx = new G4PrimaryVertex(x1, y1, z1, t1);

    for (pitr = (*vitr)->particles_begin(HepMC::children); pitr != (*vitr)->particles_end(HepMC::children); ++pitr) {
      if (transportedDescendants.count(*pitr))
        continue;
      int status = (*pitr)->status();
      int pdg = (*pitr)->pdg_id();
      bool hasDecayVertex = (nullptr != (*pitr)->end_vertex());
      bool isMuon = debugMuonPrimaries_ && std::abs(pdg) == 13 && status == 1;

      // Filter on allowed particle species if required
      if (fPDGFilter) {
        bool isInTheList = IsInTheFilterList(pdg);
        if ((!pdgFilterSel && isInTheList) || (pdgFilterSel && !isInTheList)) {
          if (isMuon)
            ++muonFirstFailure["pdg_filter"];
          if (isMuon)
            std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                      << " stage=pdg-filter decision=rejected barcode=" << (*pitr)->barcode() << " pdgId=" << pdg
                      << " status=" << status << " inList=" << isInTheList << " selectList=" << pdgFilterSel
                      << " reason=pdg_filter" << std::endl;
          if (0 < verbose)
            edm::LogVerbatim("SimG4CoreGenerator")
                << " Skiped GenParticle barcode= " << (*pitr)->barcode() << " PDGid= " << pdg << " status= " << status
                << " isExotic: " << isExotic(pdg) << " isExoticNotDet: " << isExoticNonDetectable(pdg)
                << " isInTheList: " << isInTheList << " hasDecayVertex: " << hasDecayVertex;
          continue;
        }
      }

      if (0 < verbose) {
        edm::LogVerbatim("SimG4CoreGenerator")
            << "Generator: pdg= " << pdg << " status= " << status << " hasPreDefinedDecay: " << hasDecayVertex
            << "\n           isExotic: " << isExotic(pdg) << " isExoticNotDet: " << isExoticNonDetectable(pdg)
            << " isInTheList: " << IsInTheFilterList(pdg) << "\n"
            << " MaxZCentralCMS = " << maxZCentralCMS / CLHEP::m << " m;  (x,y,z,t): (" << x1 << "," << y1 << "," << z1
            << "," << t1 << ")";
      }
      if (status > 3 && isExotic(pdg) && (!(isExoticNonDetectable(pdg)))) {
        status = hasDecayVertex ? 2 : 1;
      }
      if (status == 2 && abs(pdg) == 9900015) {
        status = 3;
      }

      // this particle has predefined decay but has no vertex
      if (2 == status && !hasDecayVertex) {
        edm::LogWarning("SimG4CoreGenerator: in event ")
            << g4evt->GetEventID() << " a particle "
            << " pdgid= " << pdg << " has status=2 but has no decay vertex, so will be fully tracked by Geant4";
        status = 1;
      }

      double x2 = x1;
      double y2 = y1;
      double z2 = z1;
      double decay_length = 0.0;
      if (2 == status) {
        x2 = (*pitr)->end_vertex()->position().x();
        y2 = (*pitr)->end_vertex()->position().y();
        z2 = (*pitr)->end_vertex()->position().z();
        decay_length = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
      }

      bool toBeAdded = !fFiductialCuts;

      double px = (*pitr)->momentum().px();
      double py = (*pitr)->momentum().py();
      double pz = (*pitr)->momentum().pz();
      double ptot = std::sqrt(px * px + py * py + pz * pz);
      math::XYZTLorentzVector p(px, py, pz, (*pitr)->momentum().e());

      double ximpact = x1;
      double yimpact = y1;
      double zimpact = z1;

      // protection against numerical problems for extremely low momenta
      // compute impact point at transition to Hector
      const double minTan = 1.e-20;
      if (std::abs(z1) < Z_hector && std::abs(pz) >= minTan * ptot) {
        zimpact = (pz > 0.0) ? Z_hector : -Z_hector;
        double del = (zimpact - z1) / pz;
        ximpact += del * px;
        yimpact += del * py;
      }
      double rimpact2 = ximpact * ximpact + yimpact * yimpact;

      if (verbose > 2)
        LogDebug("SimG4CoreGenerator") << "Processing GenParticle barcode= " << (*pitr)->barcode() << " pdg= " << pdg
                                       << " status= " << (*pitr)->status() << " st= " << status
                                       << " rimpact(cm)= " << std::sqrt(rimpact2) / CLHEP::cm
                                       << " zimpact(cm)= " << zimpact / CLHEP::cm << " ptot(GeV)= " << ptot
                                       << " pz(GeV)= " << pz;

      // Particles of status 1 trasnported along the beam pipe
      // HECTOR transport of protons are done in corresponding PPS producer
      if (1 == status && std::abs(zimpact) >= Z_hector && rimpact2 <= theDecRCut2) {
        // Very forward n, nbar, and gamma are allowed. The SHIFT extension is
        // limited to stable muons produced beyond the HECTOR boundary and
        // directed from either SHIFT side back toward CMS.
        const bool inwardSHIFTMuon = (13 == std::abs(pdg) && std::abs(z1) >= Z_hector && z1 * pz < 0.0);
        const bool externalSourceParticle = trackAllHadronDecays_ && std::abs(z1) >= maxZCentralCMS &&
                                            (std::abs(pdg) == 13 || HepPDT::ParticleID(pdg).isHadron());
        toBeAdded = (2112 == std::abs(pdg) || 22 == pdg || inwardSHIFTMuon || externalSourceParticle);
        if (externalSourceParticle && G4ParticleTable::GetParticleTable()->FindParticle(pdg) == nullptr)
          throw SimG4Exception("TrackAllHadronDecays cannot transport an undefined external Geant4 particle");
        if (verbose > 1) {
          edm::LogVerbatim("SimG4CoreGenerator")
              << "GenParticle barcode = " << (*pitr)->barcode() << " very forward; to be added: " << toBeAdded;
        }
      } else {
        // Standard case: particles not decayed by the generator and not forward
        if (1 == status && (std::abs(zimpact) < Z_hector || rimpact2 > theDecRCut2)) {
          // Ptot cut for all particles
          if (fPCuts && (ptot < theMinPCut || ptot > theMaxPCut)) {
            if (isMuon)
              ++muonFirstFailure["p_cut"];
            if (isMuon)
              std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                        << " stage=momentum-cut decision=rejected barcode=" << (*pitr)->barcode() << " p=" << ptot
                        << " GeV allowed=[" << theMinPCut << "," << theMaxPCut << "] GeV reason=p_cut" << std::endl;
            continue;
          }
          // phi cut is applied mainly for particle gun
          if (fPhiCuts) {
            double phi = p.phi();
            if (phi < theMinPhiCut || phi > theMaxPhiCut) {
              if (isMuon)
                ++muonFirstFailure["phi_cut"];
              if (isMuon)
                std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                          << " stage=phi-cut decision=rejected barcode=" << (*pitr)->barcode() << " phi=" << phi
                          << " allowed=[" << theMinPhiCut << "," << theMaxPhiCut << "] reason=phi_cut" << std::endl;
              continue;
            }
          }
          // eta cut is applied if position of the decay
          // is within vacuum chamber and limited in Z
          if (fEtaCuts) {
            // eta cut
            double xi = x1;
            double yi = y1;
            double zi = z1;

            // can be propagated along Z
            if (std::abs(pz) >= minTan * ptot) {
              if ((zi >= Z_lmax) & (pz < 0.0)) {
                zi = Z_lmax;
              } else if ((zi <= Z_lmin) & (pz > 0.0)) {
                zi = Z_lmin;
              } else {
                if (pz > 0) {
                  zi = Z_lmax;
                } else {
                  zi = Z_lmin;
                }
              }
              double del = (zi - z1) / pz;
              xi += del * px;
              yi += del * py;
            }
            // check eta cut
            if ((zi >= Z_lmin) & (zi <= Z_lmax) & (xi * xi + yi * yi < theDecRCut2)) {
              if (isMuon)
                ++muonFirstFailure["eta_fiducial_cut"];
              if (isMuon)
                std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                          << " stage=eta-cut decision=rejected barcode=" << (*pitr)->barcode() << " projected=("
                          << xi / CLHEP::mm << "," << yi / CLHEP::mm << "," << zi / CLHEP::mm
                          << ") mm reason=eta_fiducial_cut" << std::endl;
              continue;
            }
          }
          if (fLumiFilter && !fLumiFilter->isGoodForLumiMonitor(*pitr)) {
            if (isMuon)
              ++muonFirstFailure["lumi_monitor_filter"];
            if (isMuon)
              std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                        << " stage=lumi-filter decision=rejected barcode=" << (*pitr)->barcode()
                        << " reason=lumi_monitor_filter" << std::endl;
            continue;
          }
          toBeAdded = true;
          if (verbose > 1)
            edm::LogVerbatim("SimG4CoreGenerator")
                << "GenParticle barcode = " << (*pitr)->barcode() << " passed case 1";

          // Decay chain outside the fiducial cylinder defined by theRDecLenCut
          // are used for Geant4 tracking with predefined decay channel
          // In the case of decay in vacuum particle is not tracked by Geant4
        } else if (2 == status &&
                   (hasTrackedDecayDaughter(*pitr) || hasFiniteHadronFlight(*pitr) ||
                    (x2 * x2 + y2 * y2 >= theDecRCut2 && (fSlepton || std::abs(z2) < Z_hector)))) {
          toBeAdded = true;
          if (verbose > 1)
            edm::LogVerbatim("SimG4CoreGenerator") << "GenParticle barcode = " << (*pitr)->barcode() << " passed case 2"
                                                   << " decay_length(cm)= " << decay_length / CLHEP::cm;
        }
      }
      if (isMuon && !toBeAdded) {
        const char *reason = "forward_transport_not_allowed_for_muon";
        ++muonFirstFailure[reason];
        std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                  << " stage=primary-selection decision=rejected barcode="
                  << (*pitr)->barcode() << " status=" << status << " vertex=(" << x1 / CLHEP::mm << ","
                  << y1 / CLHEP::mm << "," << z1 / CLHEP::mm << ") mm impact=(" << ximpact / CLHEP::mm << ","
                  << yimpact / CLHEP::mm << "," << zimpact / CLHEP::mm << ") mm rImpact="
                  << std::sqrt(rimpact2) / CLHEP::mm << " mm Z_hector=" << Z_hector / CLHEP::mm
                  << " mm reason=" << reason << std::endl;
      }
      if (toBeAdded) {
        G4PrimaryParticle *g4prim = new G4PrimaryParticle(pdg, px * CLHEP::GeV, py * CLHEP::GeV, pz * CLHEP::GeV);

        if (g4prim->GetG4code() != nullptr) {
          g4prim->SetMass(g4prim->GetG4code()->GetPDGMass());
          double charge = g4prim->GetG4code()->GetPDGCharge();

          // apply Pt cut
          if (fPtransCut && 1 == status && 0.0 != charge && px * px + py * py < theMinPtCut2) {
            if (isMuon)
              ++muonFirstFailure["pt_cut"];
            if (isMuon)
              std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                        << " stage=pt-cut decision=rejected barcode=" << (*pitr)->barcode()
                        << " pt=" << std::sqrt(px * px + py * py) << " GeV minPt=" << std::sqrt(theMinPtCut2)
                        << " GeV charge=" << charge / CLHEP::eplus << " reason=pt_cut" << std::endl;
            delete g4prim;
            continue;
          }
          g4prim->SetCharge(charge);
        }

        if (trackAllHadronDecays_)
          sim::setRecordedPrimaryMomentum(*g4prim, (*pitr)->momentum(), (*pitr)->generated_mass());

        // V.I. do not use SetWeight but the same code
        // value of the code compute inside TrackWithHistory
        // g4prim->SetWeight( 10000*(*vpitr)->barcode() ) ;
        setGenId(g4prim, (*pitr)->barcode());

        if (2 == status) {
          particleAssignDaughters(g4prim, (HepMC::GenParticle *)*pitr, decay_length);
        }
        if (verbose > 1)
          g4prim->Print();

        ++ng4par;
        g4vtx->SetPrimary(g4prim);
        if (isMuon) {
          ++nG4MuonPrimaries;
          const auto *definition = g4prim->GetG4code();
          std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                    << " stage=g4-primary decision=accepted barcode=" << (*pitr)->barcode()
                    << " pdgId=" << g4prim->GetPDGcode()
                    << " name=" << (definition != nullptr ? definition->GetParticleName() : "unknown")
                    << " px=" << g4prim->GetPx() / CLHEP::GeV << " py=" << g4prim->GetPy() / CLHEP::GeV
                    << " pz=" << g4prim->GetPz() / CLHEP::GeV << " mass=" << g4prim->GetMass() / CLHEP::GeV
                    << " charge=" << g4prim->GetCharge() / CLHEP::eplus << " vertex=("
                    << g4vtx->GetX0() / CLHEP::mm << "," << g4vtx->GetY0() / CLHEP::mm << ","
                    << g4vtx->GetZ0() / CLHEP::mm << ") mm time=" << g4vtx->GetT0() / CLHEP::ns << " ns"
                    << std::endl;
        }
        edm::LogVerbatim("SimG4CoreGenerator") << "   " << ng4par << ". new Geant4 particle pdg= " << pdg
                                               << " Ptot(GeV/c)= " << ptot << " Pt= " << std::sqrt(px * px + py * py)
                                               << " status= " << status << "; dir= " << g4prim->GetMomentumDirection();
      }
    }

    if (verbose > 1)
      g4vtx->Print();
    g4evt->AddPrimaryVertex(g4vtx);
    ++ng4vtx;
  }

  // Add a protection for completely empty events (produced by LHCTransport):
  // add a dummy vertex with no particle attached to it
  if (ng4vtx == 0) {
    G4PrimaryVertex *g4vtx = new G4PrimaryVertex(0.0, 0.0, 0.0, 0.0);
    if (verbose > 1)
      g4vtx->Print();

    g4evt->AddPrimaryVertex(g4vtx);
  }

  if (debugMuonPrimaries_)
    std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
              << " stage=event-summary decision=complete status1_muons=" << nHepMCMuons
              << " GEANT4_muon_primaries=" << nG4MuonPrimaries << " GEANT4_primaries=" << ng4par
              << " GEANT4_vertices=" << ng4vtx << std::endl;
  if (debugMuonPrimaries_) {
    for (const auto &failure : muonFirstFailure)
      std::cout << "[FixedTargetMuonDebug][HepMC] event=" << g4evt->GetEventID()
                << " stage=first-failure-summary cut=" << failure.first << " muons=" << failure.second << std::endl;
  }

  edm::LogVerbatim("SimG4CoreGenerator") << "The list of Geant4 primaries includes " << ng4par << " particles in "
                                         << ng4vtx << " vertex";
  delete evt;
}

void Generator::particleAssignDaughters(G4PrimaryParticle *g4p, HepMC::GenParticle *vp, double decaylength) {
  if (verbose > 1) {
    LogDebug("SimG4CoreGenerator") << "Special case of long decay length \n"
                                   << "Assign daughters with to mother with decaylength=" << decaylength / CLHEP::cm
                                   << " cm";
  }
  math::XYZTLorentzVector p(vp->momentum().px(), vp->momentum().py(), vp->momentum().pz(), vp->momentum().e());

  // defined preassigned decay time
  double proper_time;
  if (trackAllHadronDecays_ && (sim::promptDecay::exactlyPrompt(*vp) ||
      (g4p->GetG4code()==nullptr && sim::promptDecay::withinPromptBudget(*vp,maxUntrackedHadronFlightMm_)))) {
    // A recorded prompt particle at rest has zero lifetime, not 0/0.
    proper_time = 0.;
  } else if (hasFiniteHadronFlight(vp) && p.P2() == 0.) {
    // The original at-rest lifetime is encoded in ct; the moving formula is
    // undefined at beta=0. Do not change the old/default moving-parent clock.
    proper_time = (vp->end_vertex()->position().t() - vp->production_vertex()->position().t()) * CLHEP::mm / c_light;
  } else {
    proper_time = decaylength / (p.Beta() * p.Gamma() * c_light);
  }
  g4p->SetProperTime(proper_time);

  if (verbose > 2) {
    LogDebug("SimG4CoreGenerator") << " px= " << p.px() << " py= " << p.py() << " pz= " << p.pz() << " e= " << p.e()
                                   << " beta= " << p.Beta() << " gamma= " << p.Gamma()
                                   << " Proper time= " << proper_time / CLHEP::ns << " ns";
  }

  // the particle will decay after the same length if it
  // has not interacted before
  double x1 = vp->end_vertex()->position().x();
  double y1 = vp->end_vertex()->position().y();
  double z1 = vp->end_vertex()->position().z();

  for (HepMC::GenVertex::particle_iterator vpdec = vp->end_vertex()->particles_begin(HepMC::children);
       vpdec != vp->end_vertex()->particles_end(HepMC::children);
       ++vpdec) {
    // transform decay products such that in the rest frame of mother
    math::XYZTLorentzVector pdec(
        (*vpdec)->momentum().px(), (*vpdec)->momentum().py(), (*vpdec)->momentum().pz(), (*vpdec)->momentum().e());

    // children should only be taken into account once
    G4PrimaryParticle *g4daught =
        new G4PrimaryParticle((*vpdec)->pdg_id(), pdec.x() * CLHEP::GeV, pdec.y() * CLHEP::GeV, pdec.z() * CLHEP::GeV);

    if (g4daught->GetG4code() != nullptr) {
      g4daught->SetMass(g4daught->GetG4code()->GetPDGMass());
      g4daught->SetCharge(g4daught->GetG4code()->GetPDGCharge());
    }

    if (trackAllHadronDecays_)
      sim::setRecordedPrimaryMomentum(*g4daught, pdec, (*vpdec)->generated_mass());

    // V.I. do not use SetWeight but the same code
    // value of the code compute inside TrackWithHistory
    setGenId(g4daught, (*vpdec)->barcode());

    int status = (*vpdec)->status();
    if (verbose > 1)
      LogDebug("SimG4CoreGenerator::::particleAssignDaughters")
          << "Assigning a " << (*vpdec)->pdg_id() << " as daughter of a " << vp->pdg_id() << " status=" << status;

    bool isInList;
    if (fSlepton) {
      std::vector<int> fParticleList = {1000011, 1000013, 1000015, 2000011, 2000013, 2000015};
      isInList = (std::find(fParticleList.begin(), fParticleList.end(), std::abs(vp->pdg_id())) != fParticleList.end());
    } else {
      isInList = std::abs(vp->pdg_id()) == 1000015;
    }
    if ((status == 2 || (status == 23 && isInList) || (status > 50 && status < 100)) &&
        (*vpdec)->end_vertex() != nullptr) {
      double x2 = (*vpdec)->end_vertex()->position().x();
      double y2 = (*vpdec)->end_vertex()->position().y();
      double z2 = (*vpdec)->end_vertex()->position().z();
      double dd = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
      particleAssignDaughters(g4daught, *vpdec, dd);
    }
    (*vpdec)->set_status(1000 + status);
    g4p->SetDaughter(g4daught);

    if (verbose > 1)
      g4daught->Print();
  }
}

// Used for non-beam particles
bool Generator::particlePassesPrimaryCuts(const G4ThreeVector &p) const {
  bool flag = true;
  double ptot = p.mag();
  if (fPCuts && (ptot < theMinPCut * CLHEP::GeV || ptot > theMaxPCut * CLHEP::GeV)) {
    flag = false;
  }
  if (fEtaCuts && flag) {
    double pz = p.z();
    if (ptot < pz + 1.e-10) {
      flag = false;

    } else {
      double eta = 0.5 * G4Log((ptot + pz) / (ptot - pz));
      if (eta < theMinEtaCut || eta > theMaxEtaCut) {
        flag = false;
      }
    }
  }
  if (fPhiCuts && flag) {
    double phi = p.phi();
    if (phi < theMinPhiCut || phi > theMaxPhiCut) {
      flag = false;
    }
  }

  if (verbose > 2)
    LogDebug("SimG4CoreGenerator") << "Generator ptot(GeV)= " << ptot / CLHEP::GeV << " eta= " << p.eta()
                                   << "  phi= " << p.phi() << " Flag= " << flag;

  return flag;
}

bool Generator::isExotic(int pdgcode) const {
  int pdgid = std::abs(pdgcode);
  return ((pdgid >= 1000000 && pdgid < 4000000 && pdgid != 3000022) ||  // SUSY, R-hadron, and technicolor particles
          pdgid == 17 ||                                                // 4th generation lepton
          pdgid == 34 ||                                                // W-prime
          pdgid == 37);                                                 // charged Higgs
}

bool Generator::isExoticNonDetectable(int pdgcode) const {
  int pdgid = std::abs(pdgcode);
  HepPDT::ParticleID pid(pdgcode);
  int charge = pid.threeCharge();
  return ((charge == 0) && (pdgid >= 1000000 && pdgid < 1000040));  // SUSY
}

bool Generator::IsInTheFilterList(int pdgcode) const {
  int pdgid = std::abs(pdgcode);
  for (auto &pdg : pdgFilter) {
    if (std::abs(pdg) == pdgid) {
      return true;
    }
  }
  return false;
}

void Generator::nonCentralEvent2G4(const HepMC::GenEvent *evt, G4Event *g4evt) {
  int i = g4evt->GetNumberOfPrimaryVertex();
  for (HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
    ++i;
    HepMC::GenParticle *gp = (*it);

    // storing only particle with status == 1
    if (gp->status() != 1)
      continue;

    int pdg = gp->pdg_id();
    G4PrimaryParticle *g4p = new G4PrimaryParticle(
        pdg, gp->momentum().px() * CLHEP::GeV, gp->momentum().py() * CLHEP::GeV, gp->momentum().pz() * CLHEP::GeV);
    if (g4p->GetG4code() != nullptr) {
      g4p->SetMass(g4p->GetG4code()->GetPDGMass());
      g4p->SetCharge(g4p->GetG4code()->GetPDGCharge());
    }
    setGenId(g4p, i);
    G4PrimaryVertex *v = new G4PrimaryVertex(gp->production_vertex()->position().x() * CLHEP::mm,
                                             gp->production_vertex()->position().y() * CLHEP::mm,
                                             gp->production_vertex()->position().z() * CLHEP::mm,
                                             gp->production_vertex()->position().t() * CLHEP::mm / CLHEP::c_light);
    v->SetPrimary(g4p);
    g4evt->AddPrimaryVertex(v);
    if (verbose > 0)
      v->Print();
  }  // end loop on HepMC particles
}
