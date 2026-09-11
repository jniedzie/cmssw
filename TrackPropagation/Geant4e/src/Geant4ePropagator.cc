#include <sstream>
#include <array>
#include <map>
#include <iomanip>

// Geant4e
#include "TrackPropagation/Geant4e/interface/ConvertFromToCLHEP.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"

// CMSSW
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TrajectoryState/interface/SurfaceSideDefinition.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/AnalyticalJacobians/interface/AnalyticalCurvilinearJacobian.h"

// Geant4
#include "G4Box.hh"
#include "G4ErrorCylSurfaceTarget.hh"
#include "G4ErrorFreeTrajState.hh"
#include "G4ErrorPlaneSurfaceTarget.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorRunManagerHelper.hh"
#include "G4EventManager.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4SteppingControl.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4ErrorPropagationNavigator.hh"
#include "G4RunManagerKernel.hh"
#include "G4StateManager.hh"
#include "G4Step.hh"
#include "G4EnergyLossForExtrapolator.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4VSolid.hh"

// CLHEP
#include <CLHEP/Units/SystemOfUnits.h>

/** Constructor.
 */
Geant4ePropagator::Geant4ePropagator(const MagneticField *field,
                                     std::string particleName,
                                     PropagationDirection dir,
                                     double plimit,
                                     double maximumStepLengthMm,
                                     double maximumPathLengthCm)
    : Propagator(dir),
      theField(field),
      theParticleName(particleName),
      theG4eManager(G4ErrorPropagatorManager::GetErrorPropagatorManager()),
      theG4eData(G4ErrorPropagatorData::GetErrorPropagatorData()),
      plimit_(plimit),
      maximumStepLengthMm_(maximumStepLengthMm),
      maximumPathLengthCm_(maximumPathLengthCm) {
  LogDebug("Geant4e") << "Geant4e Propagator initialized";

  // has to be called here, doing it later will not load the G4 physics list
  // properly when using the G4 ES Producer. Reason: unclear
  ensureGeant4eIsInitilized(true);
}

/** Destructor.
 */
Geant4ePropagator::~Geant4ePropagator() {
  LogDebug("Geant4e") << "Geant4ePropagator::~Geant4ePropagator()" << std::endl;

  // don't close the g4 Geometry here, because the propagator might have been
  // cloned
  // but there is only one, globally-shared Geometry
}

//
////////////////////////////////////////////////////////////////////////////
//

/** Propagate from a free state (e.g. position and momentum in
 *  in global cartesian coordinates) to a plane.
 */

void Geant4ePropagator::ensureGeant4eIsInitilized(bool) const {
  LogDebug("Geant4ePropagator") << "G4 propagator starts isInitialized, theField: " << theField;

  auto man = G4RunManagerKernel::GetRunManagerKernel();
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_PreInit) {
    man->SetVerboseLevel(0);
    theG4eManager->InitGeant4e();

    // A smaller step can be required when crossing sharp magnetic-field
    // transitions; expose the limit instead of fixing it globally at 10 mm.
    std::ostringstream command;
    command << "/geant4e/limits/stepLength " << maximumStepLengthMm_ << " mm";
    G4UImanager::GetUIpointer()->ApplyCommand(command.str());
  }
  const G4Field *field = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
  if (field == nullptr) {
    edm::LogError("Geant4e") << "No G4 magnetic field defined";
  }
  LogDebug("Geant4ePropagator") << "G4 propagator initialized; field: " << field;
}

template <>
Geant4ePropagator::ErrorTargetPair Geant4ePropagator::transformToG4SurfaceTarget(const Plane &pDest,
                                                                                 bool moveTargetToEndOfSurface) const {
  //* Get position and normal (orientation) of the destination plane
  GlobalPoint posPlane = pDest.toGlobal(LocalPoint(0, 0, 0));
  GlobalVector normalPlane = pDest.toGlobal(LocalVector(0, 0, 1.));
  normalPlane = normalPlane.unit();

  //* Transform this into HepGeom::Point3D<double>  and
  // HepGeom::Normal3D<double>  that define a plane for
  //  Geant4e.
  //  CMS uses cm and GeV while Geant4 uses mm and MeV
  HepGeom::Point3D<double> surfPos = TrackPropagation::globalPointToHepPoint3D(posPlane);
  HepGeom::Normal3D<double> surfNorm = TrackPropagation::globalVectorToHepNormal3D(normalPlane);

  //* Set the target surface
  return ErrorTargetPair(false, std::make_shared<G4ErrorPlaneSurfaceTarget>(surfNorm, surfPos));
}

template <>
Geant4ePropagator::ErrorTargetPair Geant4ePropagator::transformToG4SurfaceTarget(const Cylinder &pDest,
                                                                                 bool moveTargetToEndOfSurface) const {
  // Get Cylinder parameters.
  // CMS uses cm and GeV while Geant4 uses mm and MeV.
  // - Radius
  G4float radCyl = pDest.radius() * CLHEP::cm;
  // - Position: PositionType & GlobalPoint are Basic3DPoint<float,GlobalTag>
  G4ThreeVector posCyl = TrackPropagation::globalPointToHep3Vector(pDest.position());
  // - Rotation: Type in CMSSW is RotationType == TkRotation<T>, T=float
  G4RotationMatrix rotCyl = TrackPropagation::tkRotationFToHepRotation(pDest.rotation());

  // DEBUG
  TkRotation<float> rotation = pDest.rotation();
  LogDebug("Geant4e") << "G4e -  TkRotation" << rotation;
  LogDebug("Geant4e") << "G4e -  G4Rotation" << rotCyl << "mm";

  return ErrorTargetPair(!moveTargetToEndOfSurface, std::make_shared<G4ErrorCylSurfaceTarget>(radCyl, posCyl, rotCyl));
}

template <>
std::string Geant4ePropagator::getSurfaceType(Cylinder const &c) const {
  return "Cylinder";
}

template <>
std::string Geant4ePropagator::getSurfaceType(Plane const &c) const {
  return "Plane";
}

std::string Geant4ePropagator::generateParticleName(int charge) const {
  std::string particleName = theParticleName;

  if (charge > 0) {
    particleName += "+";
  }
  if (charge < 0) {
    particleName += "-";
  }

  LogDebug("Geant4e") << "G4e -  Particle name: " << particleName;

  return particleName;
}

template <>
bool Geant4ePropagator::configureAnyPropagation(G4ErrorMode &mode,
                                                Plane const &pDest,
                                                GlobalPoint const &cmsInitPos,
                                                GlobalVector const &cmsInitMom) const {
  if (cmsInitMom.mag() < plimit_)
    return false;
  if (pDest.localZ(cmsInitPos) * pDest.localZ(cmsInitMom) < 0) {
    mode = G4ErrorMode_PropForwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'forwards\' indirect "
                           "via the Any direction"
                        << std::endl;
  } else {
    mode = G4ErrorMode_PropBackwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'backwards\' indirect "
                           "via the Any direction"
                        << std::endl;
  }

  return true;
}

template <>
bool Geant4ePropagator::configureAnyPropagation(G4ErrorMode &mode,
                                                Cylinder const &pDest,
                                                GlobalPoint const &cmsInitPos,
                                                GlobalVector const &cmsInitMom) const {
  if (cmsInitMom.mag() < plimit_)
    return false;
  //------------------------------------
  // For cylinder assume outside is backwards, inside is along
  // General use for particles from collisions
  LocalPoint lpos = pDest.toLocal(cmsInitPos);
  Surface::Side theSide = pDest.side(lpos, 0);
  if (theSide == SurfaceOrientation::positiveSide) {  // outside cylinder
    mode = G4ErrorMode_PropBackwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'backwards\' indirect "
                           "via the Any direction";
  } else {  // inside cylinder
    mode = G4ErrorMode_PropForwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'forwards\' indirect "
                           "via the Any direction";
  }

  return true;
}

template <class SurfaceType>
bool Geant4ePropagator::configurePropagation(G4ErrorMode &mode,
                                             SurfaceType const &pDest,
                                             GlobalPoint const &cmsInitPos,
                                             GlobalVector const &cmsInitMom) const {
  if (cmsInitMom.mag() < plimit_)
    return false;
  if (propagationDirection() == oppositeToMomentum) {
    mode = G4ErrorMode_PropBackwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'backwards\' " << std::endl;
  } else if (propagationDirection() == alongMomentum) {
    mode = G4ErrorMode_PropForwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'forwards\'" << std::endl;
  } else if (propagationDirection() == anyDirection) {
    if (configureAnyPropagation(mode, pDest, cmsInitPos, cmsInitMom) == false)
      return false;
  } else {
    edm::LogError("Geant4e") << "G4e - Unsupported propagation mode";
    return false;
  }
  return true;
}

template <class SurfaceType>
std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateGeneric(const FreeTrajectoryState &ftsStart,
                                                                                const SurfaceType &pDest) const {
  ///////////////////////////////
  // Construct the target surface
  //
  //* Set the target surface

  ErrorTargetPair g4eTarget_center = transformToG4SurfaceTarget(pDest, false);

  // * Get the starting point and direction and convert them to
  // CLHEP::Hep3Vector
  //   for G4. CMS uses cm and GeV while Geant4 uses mm and MeV
  GlobalPoint cmsInitPos = ftsStart.position();
  GlobalVector cmsInitMom = ftsStart.momentum();
  bool flipped = false;
  if (propagationDirection() == oppositeToMomentum) {
    // flip the momentum vector as Geant4 will not do this
    // on it's own in a backward propagation
    cmsInitMom = -cmsInitMom;
    flipped = true;
  }

  // Set the mode of propagation according to the propagation direction
  G4ErrorMode mode = G4ErrorMode_PropForwards;
  if (!configurePropagation(mode, pDest, cmsInitPos, cmsInitMom))
    return TsosPP(TrajectoryStateOnSurface(), 0.0f);

  // re-check propagation direction chosen in case of AnyDirection
  if (mode == G4ErrorMode_PropBackwards && !flipped)
    cmsInitMom = -cmsInitMom;

  CLHEP::Hep3Vector g4InitPos = TrackPropagation::globalPointToHep3Vector(cmsInitPos);
  CLHEP::Hep3Vector g4InitMom = TrackPropagation::globalVectorToHep3Vector(cmsInitMom * CLHEP::GeV);

  debugReportTrackState("intitial", cmsInitPos, g4InitPos, cmsInitMom, g4InitMom, pDest);

  // Set the mode of propagation according to the propagation direction
  // G4ErrorMode mode = G4ErrorMode_PropForwards;

  // if (!configurePropagation(mode, pDest, cmsInitPos, cmsInitMom))
  //	return TsosPP(TrajectoryStateOnSurface(), 0.0f);

  ///////////////////////////////
  // Set the error and trajectories, and finally propagate
  //
  G4ErrorTrajErr g4error(5, 1);
  if (ftsStart.hasError()) {
    CurvilinearTrajectoryError initErr;
    initErr = ftsStart.curvilinearError();
    if (consistentBackwardCovariance_ && mode == G4ErrorMode_PropBackwards)
      initErr = CurvilinearTrajectoryError(TrackPropagation::reverseMomentumCovariance(initErr.matrix()));
    g4error = TrackPropagation::algebraicSymMatrix55ToG4ErrorTrajErr(initErr, ftsStart.charge());
    LogDebug("Geant4e") << "CMS -  Error matrix: " << std::endl << initErr.matrix();
  } else {
    LogDebug("Geant4e") << "No error matrix available" << std::endl;
    return TsosPP(TrajectoryStateOnSurface(), 0.0f);
  }

  LogDebug("Geant4e") << "G4e -  Error matrix: " << std::endl << g4error;

  // in CMSSW, the state errors are deflated when performing the backward
  // propagation
  if (mode == G4ErrorMode_PropForwards || consistentBackwardCovariance_) {
    // With an explicitly reversed frame, Geant4 advances along its own
    // momentum. Deflation would reverse the Jacobian a second time.
    G4ErrorPropagatorData::GetErrorPropagatorData()->SetStage(G4ErrorStage_Inflation);
  } else if (mode == G4ErrorMode_PropBackwards) {
    G4ErrorPropagatorData::GetErrorPropagatorData()->SetStage(G4ErrorStage_Deflation);
  }

  G4ErrorFreeTrajState g4eTrajState(generateParticleName(ftsStart.charge()), g4InitPos, g4InitMom, g4error);
  LogDebug("Geant4e") << "G4e -  Traj. State: " << (g4eTrajState);

  //////////////////////////////
  // Propagate
  int iterations = 0;
  double finalPathLength = 0;
  std::map<std::string, std::array<double, 3>> materialAudit;
  if (transportAudit_)
    edm::LogVerbatim("Geant4eTransportAudit") << std::setprecision(12)
        << "BEGIN mode=" << int(mode) << " charge=" << ftsStart.charge()
        << " position=" << cmsInitPos << " momentum=" << cmsInitMom;

  HepGeom::Point3D<double> finalRecoPos;

  G4ErrorPropagatorData::GetErrorPropagatorData()->SetMode(mode);

  theG4eData->SetTarget(g4eTarget_center.second.get());
  LogDebug("Geant4e") << "Running Propagation to the RECO surface" << std::endl;

  theG4eManager->InitTrackPropagation();

  // re-initialize navigator to avoid mismatches and/or segfaults
  theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointAndSetup(
      g4InitPos, &g4InitMom, /*pRelativeSearch = */ false, /*ignoreDirection = */ false);

  std::unique_ptr<G4EnergyLossForExtrapolator> energyLoss;
  if (meanEnergyLossJacobian_)
    energyLoss = std::make_unique<G4EnergyLossForExtrapolator>(0);
  G4ErrorTrajErr previousError(5, 0);
  bool continuePropagation = true;
  while (continuePropagation) {
    iterations++;
    LogDebug("Geant4e") << std::endl << "step count " << iterations << " step length " << finalPathLength;

    // re-initialize navigator to avoid mismatches and/or segfaults
    theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointWithinVolume(g4eTrajState.GetPosition());

    if (meanEnergyLossJacobian_) previousError = g4eTrajState.GetError();
    const int ierr = theG4eManager->PropagateOneStep(&g4eTrajState, mode);
    if (ierr == 0 && energyLoss && g4eTrajState.GetG4Track()) {
      auto const* track = g4eTrajState.GetG4Track();
      auto const* step = track->GetStep();
      auto const* before = step->GetPreStepPoint();
      auto const* after = step->GetPostStepPoint();
      double const e0 = before->GetKineticEnergy(), e1 = after->GetKineticEnergy();
      double const length = step->GetStepLength();
      if (length > 0. && e0 > 0. && e1 > 0. && std::abs(e1-e0) > 1.e-12 * e0) {
        // G4ErrorFreeTrajState's helix Jacobian has unit inverse-momentum
        // derivative in field-free material, omitting the changing mean
        // energy. Differentiate G4ErrorEnergyLoss's SAME midpoint rule.
        // Reference: Geant4 v11.4.1, source/processes/electromagnetic/muons/
        // src/G4ErrorEnergyLoss.cc. Mean transport and process noise stay fixed.
        auto meanEnergy = [&](double e) {
          auto advance = [&](double v) {
            return mode == G4ErrorMode_PropBackwards
                ? energyLoss->EnergyBeforeStep(v, length, before->GetMaterial(), track->GetParticleDefinition())
                : energyLoss->EnergyAfterStep(v, length, before->GetMaterial(), track->GetParticleDefinition());
          };
          double const half = .5 * (e + advance(e));
          return e - (half - advance(half));
        };
        double const h = e0 * 1.e-4;
        double const derivative = (meanEnergy(e0+h) - meanEnergy(e0-h)) / (2.*h);
        double const mass = track->GetParticleDefinition()->GetPDGMass();
        double const p0 = std::sqrt(e0*(e0+2.*mass)), p1 = std::sqrt(e1*(e1+2.*mass));
        double const jacobian = std::pow(p0/p1,3) * (e1+mass)/(e0+mass) * derivative;
        auto const transport = g4eTrajState.GetTransfMat();
        double const delta = jacobian - transport(1,1);
        if (!std::isfinite(jacobian) || !(jacobian > 0.))
          return TsosPP(TrajectoryStateOnSurface(), 0.0f);
        auto corrected = g4eTrajState.GetError();
        // T' = T + delta e0 e0^T. Keep Q exactly, instead of rescaling it
        // together with the transported covariance.
        for (int j = 1; j <= 5; ++j) {
          double cross = 0.;
          for (int k = 1; k <= 5; ++k) cross += previousError(1,k) * transport(j,k);
          corrected(1,j) += delta * cross * (j == 1 ? 2. : 1.);
        }
        corrected(1,1) += delta*delta*previousError(1,1);
        g4eTrajState.SetError(corrected);
      }
      if (after->GetStepStatus() == fGeomBoundary && before->GetMaterial() && after->GetMaterial() &&
          before->GetMaterial() != after->GetMaterial()) {
        // At a material interface the crossing distance itself varies with
        // the incoming state. The fixed-length helix Jacobian omits this
        // term. At fixed path coordinate the interface map is
        // S = I + (f_after - f_before) n^T / (n.u).
        // Only d(1/p)/ds jumps at a material-only interface. Apply S to the
        // full covariance, including noise acquired before the boundary.
        // This is the same derivative for either orientation of the normal.
        G4bool foundNormal = false;
        G4ThreeVector normal = theG4eManager->GetErrorPropagationNavigator()->GetGlobalExitNormal(
            after->GetPosition(), &foundNormal);
        for (auto const* point : {before, after}) {
          if (foundNormal) break;
          auto const* touchable = dynamic_cast<G4TouchableHistory const*>(point->GetTouchable());
          if (!touchable || !point->GetPhysicalVolume()) continue;
          auto const& transform = touchable->GetHistory()->GetTopTransform();
          auto const local = transform.TransformPoint(after->GetPosition());
          auto const* solid = point->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
          if (solid->Inside(local) != kSurface) continue;
          normal = transform.Inverse().TransformAxis(solid->SurfaceNormal(local));
          foundNormal = normal.mag2() > 0.;
          if (foundNormal) break;
        }
        auto const direction = after->GetMomentumDirection();
        double const incidence = normal.dot(direction);
        if (!foundNormal || std::abs(incidence) < 1.e-9) {
          edm::LogWarning("Geant4eBoundaryCovariance") << "Undefined material-interface derivative at "
              << after->GetPosition() << "; rejecting propagation";
          return TsosPP(TrajectoryStateOnSurface(), 0.0f);
        }
        double const mass = track->GetParticleDefinition()->GetPDGMass();
        double const momentum = std::sqrt(e1 * (e1 + 2. * mass));
        double const stoppingBefore = energyLoss->ComputeDEDX(e1, track->GetParticleDefinition(), before->GetMaterial());
        double const stoppingAfter = energyLoss->ComputeDEDX(e1, track->GetParticleDefinition(), after->GetMaterial());
        // Geant4 error coordinates use 1/GeV and cm; stopping powers use
        // MeV/mm. Backward transport gains energy, hence the opposite sign.
        double const jump = (mode == G4ErrorMode_PropBackwards ? -1. : 1.) *
            (e1 + mass) / std::pow(momentum, 3) * (stoppingAfter - stoppingBefore) * CLHEP::GeV * CLHEP::cm;
        G4ThreeVector const transverse = G4ThreeVector(-direction.y(), direction.x(), 0.).unit();
        G4ThreeVector const vertical = direction.cross(transverse);
        G4ErrorMatrix interfaceMap(5, 5, 1);
        interfaceMap(1, 4) = jump * normal.dot(transverse) / incidence;
        interfaceMap(1, 5) = jump * normal.dot(vertical) / incidence;
        g4eTrajState.SetError(g4eTrajState.GetError().similarity(interfaceMap));
      }
    }

    if (transportAudit_ && g4eTrajState.GetG4Track()) {
      auto const* step = g4eTrajState.GetG4Track()->GetStep();
      if (step && step->GetPreStepPoint()->GetMaterial()) {
        if (step->GetPreStepPoint()->GetMaterial() != step->GetPostStepPoint()->GetMaterial()) {
          auto const* before = step->GetPreStepPoint();
          auto const* after = step->GetPostStepPoint();
          edm::LogVerbatim("Geant4eTransportAudit") << std::setprecision(12)
              << "BOUNDARY from=" << before->GetMaterial()->GetName()
              << " volume=" << (before->GetPhysicalVolume() ? before->GetPhysicalVolume()->GetName() : "none")
              << " to=" << (after->GetMaterial() ? after->GetMaterial()->GetName() : "none")
              << " volume=" << (after->GetPhysicalVolume() ? after->GetPhysicalVolume()->GetName() : "none")
              << " preMm=" << before->GetPosition() << " postMm=" << after->GetPosition();
        }
        auto& audit = materialAudit[step->GetPreStepPoint()->GetMaterial()->GetName()];
        audit[0] += step->GetStepLength() / CLHEP::cm;
        audit[1] += (step->GetPreStepPoint()->GetKineticEnergy() -
                     step->GetPostStepPoint()->GetKineticEnergy()) / CLHEP::GeV;
        audit[2] += 1.;
      }
    }

    if (ierr != 0) {
      if (transportAudit_)
        edm::LogVerbatim("Geant4eTransportAudit") << "FAILED code=" << ierr;
      // propagation failed, return invalid track state
      return TsosPP(TrajectoryStateOnSurface(), 0.0f);
    }

    const float thisPathLength = TrackPropagation::g4doubleToCmsDouble(g4eTrajState.GetG4Track()->GetStepLength());

    LogDebug("Geant4e") << "step Length was " << thisPathLength << " cm, current global position: "
                        << TrackPropagation::hepPoint3DToGlobalPoint(g4eTrajState.GetPosition()) << std::endl;

    finalPathLength += thisPathLength;

    if (std::fabs(finalPathLength) > maximumPathLengthCm_) {
      LogDebug("Geant4e") << "ERROR: Quitting propagation: path length mega large" << std::endl;
      theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
      continuePropagation = false;
      LogDebug("Geant4e") << "WARNING: Quitting propagation: max path length "
                             "exceeded, returning invalid state"
                          << std::endl;

      // reached maximum path length, bail out
      return TsosPP(TrajectoryStateOnSurface(), 0.0f);
    }

    if (theG4eManager->GetPropagator()->CheckIfLastStep(g4eTrajState.GetG4Track())) {
      theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
      continuePropagation = false;
    }
  }

  // CMSSW Tracking convention, backward propagations have negative path length
  if (propagationDirection() == oppositeToMomentum)
    finalPathLength = -finalPathLength;

  // store the correct location for the hit on the RECO surface
  LogDebug("Geant4e") << "Position on the RECO surface" << g4eTrajState.GetPosition() << std::endl;
  finalRecoPos = g4eTrajState.GetPosition();

  theG4eManager->EventTermination();

  LogDebug("Geant4e") << "Final position of the Track :" << g4eTrajState.GetPosition() << std::endl;

  //////////////////////////////
  // Retrieve the state in the end from Geant4e, convert them to CMS vectors
  // and points, and build global trajectory parameters.
  // CMS uses cm and GeV while Geant4 uses mm and MeV
  //
  const HepGeom::Vector3D<double> momEnd = g4eTrajState.GetMomentum();
  if (transportAudit_) {
    for (auto const& [material, values] : materialAudit)
      edm::LogVerbatim("Geant4eTransportAudit") << std::setprecision(12)
          << "MATERIAL name=" << material << " lengthCm=" << values[0]
          << " lossGeV=" << values[1] << " steps=" << values[2];
    edm::LogVerbatim("Geant4eTransportAudit") << std::setprecision(12)
        << "END positionMm=" << finalRecoPos << " momentumMeV=" << momEnd
        << " pathCm=" << finalPathLength;
  }

  // use the hit on the the RECO plane as the final position to be d'accor with
  // the RecHit measurements
  const GlobalPoint posEndGV = TrackPropagation::hepPoint3DToGlobalPoint(finalRecoPos);
  GlobalVector momEndGV = TrackPropagation::hep3VectorToGlobalVector(momEnd) / CLHEP::GeV;

  debugReportTrackState("final", posEndGV, finalRecoPos, momEndGV, momEnd, pDest);

  // Get the error covariance matrix from Geant4e. It comes in curvilinear
  // coordinates so use the appropiate CMS class
  G4ErrorTrajErr g4errorEnd = g4eTrajState.GetError();

  CurvilinearTrajectoryError curvError(
      TrackPropagation::g4ErrorTrajErrToAlgebraicSymMatrix55(g4errorEnd, ftsStart.charge()));

  if (mode == G4ErrorMode_PropBackwards) {
    if (consistentBackwardCovariance_)
      curvError = CurvilinearTrajectoryError(TrackPropagation::reverseMomentumCovariance(curvError.matrix()));
    GlobalTrajectoryParameters endParm(
        posEndGV, momEndGV, ftsStart.parameters().charge(), &ftsStart.parameters().magneticField());

    // flip the momentum direction because it has been flipped before running
    // G4's backwards prop
    momEndGV = -momEndGV;
  }

  LogDebug("Geant4e") << "G4e -  Error matrix after propagation: " << std::endl << g4errorEnd;

  LogDebug("Geant4e") << "CMS -  Error matrix after propagation: " << std::endl << curvError.matrix();

  GlobalTrajectoryParameters tParsDest(posEndGV, momEndGV, ftsStart.charge(), theField);

  SurfaceSideDefinition::SurfaceSide side;

  side = propagationDirection() == alongMomentum ? SurfaceSideDefinition::afterSurface
                                                 : SurfaceSideDefinition::beforeSurface;

  return TsosPP(TrajectoryStateOnSurface(tParsDest, curvError, pDest, side), finalPathLength);
}

//
////////////////////////////////////////////////////////////////////////////
//

/** The methods propagateWithPath() are identical to the corresponding
 *  methods propagate() in what concerns the resulting
 *  TrajectoryStateOnSurface, but they provide in addition the
 *  exact path length along the trajectory.
 */

std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateWithPath(const FreeTrajectoryState &ftsStart,
                                                                                 const Plane &pDest) const {
  // Finally build the pair<...> that needs to be returned where the second
  // parameter is the exact path length. Currently calculated with a stepping
  // action that adds up the length of every step
  return propagateGeneric(ftsStart, pDest);
}

std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateWithPath(const FreeTrajectoryState &ftsStart,
                                                                                 const Cylinder &cDest) const {
  // Finally build the pair<...> that needs to be returned where the second
  // parameter is the exact path length.
  return propagateGeneric(ftsStart, cDest);
}

std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateWithPath(
    const TrajectoryStateOnSurface &tsosStart, const Plane &pDest) const {
  // Finally build the pair<...> that needs to be returned where the second
  // parameter is the exact path length.
  const FreeTrajectoryState ftsStart = *tsosStart.freeState();
  return propagateGeneric(ftsStart, pDest);
}

std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateWithPath(
    const TrajectoryStateOnSurface &tsosStart, const Cylinder &cDest) const {
  const FreeTrajectoryState ftsStart = *tsosStart.freeState();
  // Finally build the pair<...> that needs to be returned where the second
  // parameter is the exact path length.
  return propagateGeneric(ftsStart, cDest);
}

void Geant4ePropagator::debugReportPlaneSetup(GlobalPoint const &posPlane,
                                              HepGeom::Point3D<double> const &surfPos,
                                              GlobalVector const &normalPlane,
                                              HepGeom::Normal3D<double> const &surfNorm,
                                              const Plane &pDest) const {
  LogDebug("Geant4e") << "G4e -  Destination CMS plane position:" << posPlane << "cm\n"
                      << "G4e -                  (Ro, eta, phi): (" << posPlane.perp() << " cm, " << posPlane.eta()
                      << ", " << posPlane.phi().degrees() << " deg)\n"
                      << "G4e -  Destination G4  plane position: " << surfPos << " mm, Ro = " << surfPos.perp()
                      << " mm";
  LogDebug("Geant4e") << "G4e -  Destination CMS plane normal  : " << normalPlane << "\n"
                      << "G4e -  Destination G4  plane normal  : " << normalPlane;
  LogDebug("Geant4e") << "G4e -  Distance from plane position to plane: " << pDest.localZ(posPlane) << " cm";
}

template <class SurfaceType>
void Geant4ePropagator::debugReportTrackState(std::string const &currentContext,
                                              GlobalPoint const &cmsInitPos,
                                              CLHEP::Hep3Vector const &g4InitPos,
                                              GlobalVector const &cmsInitMom,
                                              CLHEP::Hep3Vector const &g4InitMom,
                                              const SurfaceType &pDest) const {
  LogDebug("Geant4e") << "G4e - Current Context: " << currentContext;
  LogDebug("Geant4e") << "G4e -  CMS point position:" << cmsInitPos << "cm\n"
                      << "G4e -              (Ro, eta, phi): (" << cmsInitPos.perp() << " cm, " << cmsInitPos.eta()
                      << ", " << cmsInitPos.phi().degrees() << " deg)\n"
                      << "G4e -   G4  point position: " << g4InitPos << " mm, Ro = " << g4InitPos.perp() << " mm";
  LogDebug("Geant4e") << "G4e -   CMS momentum      :" << cmsInitMom << "GeV\n"
                      << " pt: " << cmsInitMom.perp() << "G4e -  G4  momentum      : " << g4InitMom << " MeV";
}
