#include "SimG4Core/Application/interface/TrackingAction.h"
#include "SimG4Core/Physics/interface/CMSG4TrackInterface.h"

#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfTrack.h"
#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "SimG4Core/Notification/interface/TrackWithHistory.h"
#include "SimG4Core/Notification/interface/SimTrackManager.h"
#include "SimG4Core/Notification/interface/CMSSteppingVerbose.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4UImanager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4TrackingManager.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#include <iostream>

//#define EDM_ML_DEBUG

TrackingAction::TrackingAction(SimTrackManager* stm, CMSSteppingVerbose* sv, const edm::ParameterSet& p)
    : trackManager_(stm),
      steppingVerbose_(sv),
      endPrintTrackID_(p.getParameter<int>("EndPrintTrackID")),
      checkTrack_(p.getUntrackedParameter<bool>("CheckTrack", false)),
      debugMuonTracking_(p.getUntrackedParameter<bool>("DebugMuonTracking", false)),
      debugMuonPrimaryFates_(p.getUntrackedParameter<bool>("DebugMuonPrimaryFates", false)),
      doFineCalo_(p.getParameter<bool>("DoFineCalo")),
      saveCaloBoundaryInformation_(p.getParameter<bool>("SaveCaloBoundaryInformation")),
      ekinMin_(p.getParameter<double>("PersistencyEmin") * CLHEP::GeV),
      ekinMinRegion_(p.getParameter<std::vector<double>>("RegionEmin")) {
  interface_ = CMSG4TrackInterface::instance();
  double eth = p.getParameter<double>("EminFineTrack") * CLHEP::MeV;
  if (doFineCalo_ && eth < ekinMin_) {
    ekinMin_ = eth;
  }
  edm::LogVerbatim("SimG4CoreApplication") << "TrackingAction: boundary: " << saveCaloBoundaryInformation_
                                           << "; DoFineCalo: " << doFineCalo_ << "; ekinMin(MeV)=" << ekinMin_;
  if (!ekinMinRegion_.empty()) {
    ptrRegion_.resize(ekinMinRegion_.size(), nullptr);
  }
}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
  g4Track_ = aTrack;
  currentHistory_ = new TrackWithHistory(aTrack, aTrack->GetParentID());
  interface_->setCurrentTrack(aTrack);

  BeginOfTrack bt(aTrack);
  m_beginOfTrackSignal(&bt);

  trkInfo_ = dynamic_cast<TrackInformation*>(aTrack->GetUserInformation());

  // Always save primaries
  if (nullptr != trkInfo_ && trkInfo_->isPrimary()) {
    trackManager_->cleanTracksWithHistory();
    currentHistory_->setToBeSaved();
  }

  if (nullptr != steppingVerbose_) {
    steppingVerbose_->trackStarted(aTrack, false);
    if (aTrack->GetTrackID() == endPrintTrackID_) {
      steppingVerbose_->stopEventPrint();
    }
  }
  double ekin = aTrack->GetKineticEnergy();

  bool const primaryMuon = trkInfo_ && trkInfo_->isPrimary() &&
                           std::abs(aTrack->GetDefinition()->GetPDGEncoding()) == 13;
  if ((debugMuonTracking_ && std::abs(aTrack->GetDefinition()->GetPDGEncoding()) == 13) ||
      (debugMuonPrimaryFates_ && primaryMuon)) {
    const G4Event* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    const G4VPhysicalVolume* volume = aTrack->GetVolume();
    const G4ThreeVector& position = aTrack->GetPosition();
    const G4ThreeVector& momentum = aTrack->GetMomentum();
    std::cout << "[FixedTargetMuonDebug][" << (primaryMuon ? "PrimaryFate" : "G4Track")
              << "] event=" << (event ? event->GetEventID() : -1)
              << " stage=track-start track_id=" << aTrack->GetTrackID() << " parent_id=" << aTrack->GetParentID()
              << " pdg_id=" << aTrack->GetDefinition()->GetPDGEncoding()
              << " primary=" << (trkInfo_ && trkInfo_->isPrimary())
              << " volume=" << (volume ? volume->GetName() : "outside-world")
              << " position_mm=(" << position.x() / CLHEP::mm << "," << position.y() / CLHEP::mm << ","
              << position.z() / CLHEP::mm << ")"
              << " momentum_GeV=(" << momentum.x() / CLHEP::GeV << "," << momentum.y() / CLHEP::GeV << ","
              << momentum.z() / CLHEP::GeV << ") kinetic_energy_GeV=" << ekin / CLHEP::GeV << std::endl;
  }

#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("DoFineCalo") << "PreUserTrackingAction: Start processing track " << aTrack->GetTrackID()
                                 << " pdgid=" << aTrack->GetDefinition()->GetPDGEncoding()
                                 << " ekin[GeV]=" << ekin / CLHEP::GeV << " vertex[cm]=("
                                 << aTrack->GetVertexPosition().x() / CLHEP::cm << ","
                                 << aTrack->GetVertexPosition().y() / CLHEP::cm << ","
                                 << aTrack->GetVertexPosition().z() / CLHEP::cm << ")"
                                 << " parentid=" << aTrack->GetParentID();
#endif
  if (nullptr != trkInfo_ && ekin > ekinMin_) {
    // Each track with energy above the threshold should be saved
    trkInfo_->putInHistory();
  }
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack) {
  // Tracks in history may be upgraded to stored secondary tracks,
  // which cross the boundary between Tracker and Calo
  int id = aTrack->GetTrackID();
  bool const primaryMuon = trkInfo_ && trkInfo_->isPrimary() &&
                           std::abs(aTrack->GetDefinition()->GetPDGEncoding()) == 13;
  if ((debugMuonTracking_ && std::abs(aTrack->GetDefinition()->GetPDGEncoding()) == 13) ||
      (debugMuonPrimaryFates_ && primaryMuon)) {
    const G4Event* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    const G4VPhysicalVolume* volume = aTrack->GetVolume();
    const G4Step* step = aTrack->GetStep();
    const G4VProcess* process = step ? step->GetPostStepPoint()->GetProcessDefinedStep() : nullptr;
    const G4ThreeVector& position = aTrack->GetPosition();
    std::cout << "[FixedTargetMuonDebug][" << (primaryMuon ? "PrimaryFate" : "G4Track")
              << "] event=" << (event ? event->GetEventID() : -1)
              << " stage=track-end track_id=" << id << " parent_id=" << aTrack->GetParentID()
              << " pdg_id=" << aTrack->GetDefinition()->GetPDGEncoding() << " g4_status=" << aTrack->GetTrackStatus()
              << " steps=" << aTrack->GetCurrentStepNumber()
              << " volume=" << (volume ? volume->GetName() : "outside-world")
              << " position_mm=(" << position.x() / CLHEP::mm << "," << position.y() / CLHEP::mm << ","
              << position.z() / CLHEP::mm << ") kinetic_energy_GeV=" << aTrack->GetKineticEnergy() / CLHEP::GeV
              << " global_time_ns=" << aTrack->GetGlobalTime() / CLHEP::ns
              << " last_process=" << (process ? process->GetProcessName() : "none") << std::endl;
  }
  bool ok = currentHistory_->saved();
  if (nullptr != trkInfo_) {
    ok = (ok || trkInfo_->storeTrack());
    if (trkInfo_->crossedBoundary()) {
      currentHistory_->setCrossedBoundaryPosMom(
          id, trkInfo_->getPositionAtBoundary(), trkInfo_->getMomentumAtBoundary());
      ok = (ok || saveCaloBoundaryInformation_ || doFineCalo_);
    }
  }
  if (ok) {
    currentHistory_->setToBeSaved();
  }

  bool withAncestor = false;
  bool isInHistory = false;
  if (nullptr != trkInfo_) {
    withAncestor = (trkInfo_->getIDonCaloSurface() == id || trkInfo_->isAncestor());
    isInHistory = trkInfo_->isInHistory();
  }

  trackManager_->addTrack(currentHistory_, aTrack, isInHistory, withAncestor);

#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("TrackingAction") << "TrackingAction end track=" << id << "  "
                                     << aTrack->GetDefinition()->GetParticleName() << " proposed to be saved= " << ok
                                     << " end point " << aTrack->GetPosition();
#endif

  if (!isInHistory) {
    delete currentHistory_;
  }

  EndOfTrack et(aTrack);
  m_endOfTrackSignal(&et);
}
