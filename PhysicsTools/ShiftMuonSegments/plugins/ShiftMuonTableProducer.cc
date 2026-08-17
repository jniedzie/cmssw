#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalTestNumbering.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkClonerImpl.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TRecHit2DPosConstraint.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "G4EnergyLossForExtrapolator.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Navigator.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
  class Geant4MaterialEffectsProvider final : public SteppingHelixMaterialEffects {
  public:
    bool compute(GlobalPoint const& position,
                 GlobalVector const& momentum,
                 int charge,
                 double& momentumLossPerCm,
                 double& momentumLossDerivative,
                 double& radiationLengthCm) const override {
      auto* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
      if (!navigator || !navigator->GetWorldVolume() || !(momentum.mag() > 0.))
        return false;
      G4ThreeVector const point(position.x() * cm, position.y() * cm, position.z() * cm);
      G4ThreeVector const direction(momentum.x() / momentum.mag(),
                                    momentum.y() / momentum.mag(),
                                    momentum.z() / momentum.mag());
      auto* volume = navigator->LocateGlobalPointAndSetup(point, &direction, false, false);
      if (!volume || !volume->GetLogicalVolume() || !volume->GetLogicalVolume()->GetMaterial())
        return false;
      auto const* material = volume->GetLogicalVolume()->GetMaterial();
      G4ParticleDefinition const* particle =
          charge > 0 ? static_cast<G4ParticleDefinition const*>(G4MuonPlus::MuonPlusDefinition())
                     : static_cast<G4ParticleDefinition const*>(G4MuonMinus::MuonMinusDefinition());
      constexpr double muonMassGeV = 0.1056583755;
      thread_local G4EnergyLossForExtrapolator energyLoss(0);
      auto momentumRate = [&](double momentumGeV) {
        double const energyGeV = std::sqrt(momentumGeV * momentumGeV + muonMassGeV * muonMassGeV);
        double const kineticEnergy = (energyGeV - muonMassGeV) * GeV;
        double const energyLossGeVPerCm = energyLoss.ComputeDEDX(kineticEnergy, particle, material) / (GeV / cm);
        return -energyLossGeVPerCm * energyGeV / momentumGeV;
      };
      double const momentumGeV = momentum.mag();
      momentumLossPerCm = momentumRate(momentumGeV);
      double const delta = std::max(1.e-4, 1.e-3 * momentumGeV);
      double const lower = std::max(1.e-5, momentumGeV - delta);
      double const upper = momentumGeV + delta;
      momentumLossDerivative = (momentumRate(upper) - momentumRate(lower)) / (upper - lower);
      radiationLengthCm = material->GetRadlen() / cm;
      return std::isfinite(momentumLossPerCm) && std::isfinite(momentumLossDerivative) &&
             radiationLengthCm > 0. && std::isfinite(radiationLengthCm);
    }
  };

  class GeometryMaterialPropagator final : public Propagator {
  public:
    GeometryMaterialPropagator(SteppingHelixPropagator const& covariancePropagator,
                               SteppingHelixPropagator const& comparisonPropagator,
                               double integrationStepCm,
                               bool logComparison)
        : Propagator(covariancePropagator.propagationDirection()),
          covariancePropagator_(covariancePropagator),
          comparisonPropagator_(comparisonPropagator),
          integrationStepCm_(integrationStepCm),
          logComparison_(logComparison) {}

    GeometryMaterialPropagator* clone() const override { return new GeometryMaterialPropagator(*this); }

    MagneticField const* magneticField() const override { return covariancePropagator_.magneticField(); }

    void setPropagationDirection(PropagationDirection direction) override {
      Propagator::setPropagationDirection(direction);
      covariancePropagator_.setPropagationDirection(direction);
      comparisonPropagator_.setPropagationDirection(direction);
    }

    std::pair<TrajectoryStateOnSurface, double> propagateWithPath(FreeTrajectoryState const& state,
                                                                  Plane const& surface) const override {
      return applyGeometryEnergyLoss(state,
                                     covariancePropagator_.propagateWithPath(state, surface),
                                     comparisonPropagator_.propagateWithPath(state, surface).first);
    }

    std::pair<TrajectoryStateOnSurface, double> propagateWithPath(FreeTrajectoryState const& state,
                                                                  Cylinder const& surface) const override {
      return applyGeometryEnergyLoss(state,
                                     covariancePropagator_.propagateWithPath(state, surface),
                                     comparisonPropagator_.propagateWithPath(state, surface).first);
    }

    std::pair<TrajectoryStateOnSurface, double> propagateWithPath(TrajectoryStateOnSurface const& state,
                                                                  Plane const& surface) const override {
      if (!state.isValid() || !state.freeState())
        return {TrajectoryStateOnSurface(), 0.};
      return applyGeometryEnergyLoss(*state.freeState(),
                                     covariancePropagator_.propagateWithPath(state, surface),
                                     comparisonPropagator_.propagateWithPath(state, surface).first);
    }

    std::pair<TrajectoryStateOnSurface, double> propagateWithPath(TrajectoryStateOnSurface const& state,
                                                                  Cylinder const& surface) const override {
      if (!state.isValid() || !state.freeState())
        return {TrajectoryStateOnSurface(), 0.};
      return applyGeometryEnergyLoss(*state.freeState(),
                                     covariancePropagator_.propagateWithPath(state, surface),
                                     comparisonPropagator_.propagateWithPath(state, surface).first);
    }

  private:
    std::pair<TrajectoryStateOnSurface, double> applyGeometryEnergyLoss(
        FreeTrajectoryState const& start,
        std::pair<TrajectoryStateOnSurface, double> propagated,
        TrajectoryStateOnSurface const& comparisonState) const {
      auto const& state = propagated.first;
      double const path = propagated.second;
      if (!state.isValid() || !(integrationStepCm_ > 0.) || !std::isfinite(path) || path == 0.)
        return propagated;

      auto* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
      if (!navigator || !navigator->GetWorldVolume())
        return propagated;

      GlobalPoint const startPosition = start.position();
      GlobalPoint const endPosition = state.globalPosition();
      GlobalVector const chord = endPosition - startPosition;
      double const chordLength = chord.mag();
      if (!(chordLength > 0.) || !std::isfinite(chordLength))
        return propagated;

      double const startMomentum = start.momentum().mag();
      constexpr double muonMassGeV = 0.1056583755;
      if (!(startMomentum > 0.) || !std::isfinite(startMomentum))
        return propagated;
      double energyGeV = std::sqrt(startMomentum * startMomentum + muonMassGeV * muonMassGeV);
      // The sign returned by a Propagator denotes its configured solution,
      // not necessarily the physical flight direction of a reversed SHIFT
      // trajectory.  Infer energy loss/gain directly from whether the spatial
      // displacement follows or opposes the state's momentum.
      bool const along = chord.dot(start.momentum()) > 0.;
      // Material is sampled on this chord, so its integration measure must be
      // the chord length as well.  Weighting these samples by the propagated
      // helix length can count the same material repeatedly for a curling or
      // poorly constrained trial state.
      unsigned int const steps =
          std::max(1U, static_cast<unsigned int>(std::ceil(chordLength / integrationStepCm_)));
      double const pathStepCm = chordLength / static_cast<double>(steps);
      G4ThreeVector const g4Direction(chord.x() / chordLength, chord.y() / chordLength, chord.z() / chordLength);
      G4ParticleDefinition const* particle = nullptr;
      if (start.charge() > 0)
        particle = G4MuonPlus::MuonPlusDefinition();
      else
        particle = G4MuonMinus::MuonMinusDefinition();
      // Geant4e provides its own energy-loss tables for extrapolation.  The
      // generic G4EmCalculator depends on the processes registered in the
      // active physics list and returns zero dE/dx in the Geant4e setup used
      // here, effectively turning this branch into vacuum propagation.
      thread_local G4EnergyLossForExtrapolator energyLoss(0);

      for (unsigned int index = 0; index < steps; ++index) {
        double const fraction = (static_cast<double>(index) + 0.5) / static_cast<double>(steps);
        GlobalPoint const point = startPosition + fraction * chord;
        G4ThreeVector const g4Point(point.x() * cm, point.y() * cm, point.z() * cm);
        auto* volume = navigator->LocateGlobalPointAndSetup(g4Point, &g4Direction, false, false);
        if (!volume || !volume->GetLogicalVolume() || !volume->GetLogicalVolume()->GetMaterial())
          continue;
        double const kineticEnergy = (energyGeV - muonMassGeV) * GeV;
        if (!(kineticEnergy > 0.))
          return {TrajectoryStateOnSurface(), path};
        auto const* material = volume->GetLogicalVolume()->GetMaterial();
        double const updatedKineticEnergy = along
                                                ? energyLoss.EnergyAfterStep(
                                                      kineticEnergy, pathStepCm * cm, material, particle)
                                                : energyLoss.EnergyBeforeStep(
                                                      kineticEnergy, pathStepCm * cm, material, particle);
        if (!(updatedKineticEnergy > 0.) || !std::isfinite(updatedKineticEnergy))
          return {TrajectoryStateOnSurface(), path};
        energyGeV = updatedKineticEnergy / GeV + muonMassGeV;
        if (!(energyGeV > muonMassGeV) || !std::isfinite(energyGeV))
          return {TrajectoryStateOnSurface(), path};
      }

      double const momentum = std::sqrt(energyGeV * energyGeV - muonMassGeV * muonMassGeV);
      GlobalVector const direction = state.globalMomentum().unit();
      GlobalTrajectoryParameters const parameters(
          state.globalPosition(), momentum * direction, state.charge(), state.magneticField());
      TrajectoryStateOnSurface corrected = state.hasError()
                                               ? TrajectoryStateOnSurface(
                                                     parameters, state.curvilinearError(), state.surface(), state.surfaceSide())
                                               : TrajectoryStateOnSurface(parameters, state.surface(), state.surfaceSide());
      if (logComparison_ && comparisonState.isValid()) {
        edm::LogWarning("ShiftMuonGeometryMaterialComparison")
            << "direction=" << (propagationDirection() == oppositeToMomentum ? "opposite" : "along")
            << " startP=" << startMomentum << " approximateP=" << comparisonState.globalMomentum().mag()
            << " geometryP=" << momentum << " approximateDeltaP="
            << comparisonState.globalMomentum().mag() - startMomentum << " geometryDeltaP=" << momentum - startMomentum
            << " path=" << path << " chord=" << chordLength;
      }
      return {corrected, path};
    }

    SteppingHelixPropagator covariancePropagator_;
    SteppingHelixPropagator comparisonPropagator_;
    double integrationStepCm_;
    bool logComparison_;
  };

  struct PropagatedState {
    bool valid = false;
    GlobalPoint position;
    GlobalVector momentum;
    double path = 0.;
    bool materialBoundaryValid = false;
    GlobalVector materialBoundaryMomentum;
    double materialPath = 0.;
  };

  struct HitFingerprint {
    uint32_t detId;
    float localX;
    float localY;
  };

  struct Candidate {
    reco::Track const* track;
    int source;
    unsigned int sourceIndex;
    std::vector<HitFingerprint> hitFingerprints;
    PropagatedState targetLineState;
    int timingDirectionSign = 0;
    unsigned int timingMeasurements = 0;
    double timingChi2 = 0.;
    double timingDeltaChi2 = 0.;
    int physicalDirectionSign = 0;
    bool directionalRefitAttempted = false;
    bool directionalRefitValid = false;
    unsigned int directionalRefitHits = 0;
    double directionalRefitChi2 = 0.;
    double directionalRefitNdof = 0.;
    double directionalRefitUpstreamPt = 0.;
    double directionalRefitUpstreamP = 0.;
    double directionalRefitDownstreamP = 0.;
    double directionalRefitFractionalLossAcrossHits = 0.;
    double directionalRefitSourceFacingSurfaceDistance = -1.;
    double directionalRefitFilterDirectionCosine = -2.;
    double directionalRefitTargetDirectionCosine = -2.;
    double preRefitPt = 0.;
    double preRefitPz = 0.;
    bool directionalRefitFirstValid = false;
    unsigned int directionalRefitFirstHits = 0;
    double directionalRefitFirstChi2 = 0.;
    double directionalRefitFirstNdof = 0.;
    double directionalRefitFirstUpstreamPt = 0.;
    double directionalRefitFirstQoverP = 0.;
    double directionalRefitFirstTargetPt = 0.;
    double directionalRefitFirstTargetPz = 0.;
    bool directionalRefitSecondValid = false;
    unsigned int directionalRefitSecondHits = 0;
    double directionalRefitSecondChi2 = 0.;
    double directionalRefitSecondNdof = 0.;
    double directionalRefitSecondUpstreamPt = 0.;
    double directionalRefitSecondQoverP = 0.;
    double directionalRefitSecondTargetPt = 0.;
    double directionalRefitSecondTargetPz = 0.;
    bool directionalRefitSecondConverged = false;
    double directionalRefitRelativeQoverPChange = -1.;
    int directionalRefitSelectedIteration = 0;
    double directionalRefitVacuumTargetPt = 0.;
    double directionalRefitVacuumTargetPz = 0.;
    double directionalRefitMaterialDeltaPt = 0.;
    bool directionalRefitMaterialBoundaryValid = false;
    double directionalRefitMaterialBoundaryPt = 0.;
    double directionalRefitMaterialBoundaryPz = 0.;
    double directionalRefitMaterialPath = 0.;
    unsigned int directionalRefitAllHitsInput = 0;
    unsigned int directionalRefitAllHitsOrderingFallback = 0;
    double directionalRefitAllHitsOrderingSpan = 0.;
    unsigned int directionalRefitAllHitsFirstRejected = 0;
    unsigned int directionalRefitAllHitsSecondRejected = 0;
    unsigned int directionalRefitPrecisionInput = 0;
    unsigned int directionalRefitPrecisionOrderingFallback = 0;
    double directionalRefitPrecisionOrderingSpan = 0.;
    unsigned int directionalRefitPrecisionFirstRejected = 0;
    unsigned int directionalRefitPrecisionSecondRejected = 0;
    unsigned int nDTRefitHits = 0;
    unsigned int nCSCRefitHits = 0;
    unsigned int nRPCRefitHits = 0;
    unsigned int nGEMRefitHits = 0;
    unsigned int nPrecisionRefitStations = 0;
    unsigned int nCompatibleDTSegments = 0;
    unsigned int nCompatiblePixelHits = 0;
    unsigned int nCompatibleStripHits = 0;
    unsigned int nPropagatedDTSegments = 0;
    unsigned int nPropagatedPixelHits = 0;
    unsigned int nPropagatedStripHits = 0;
    double minDTResidual = -1.;
    double minDTEstimatorChi2 = -1.;
    double minTrackerResidual = -1.;
    double minTrackerEstimatorChi2 = -1.;
    unsigned int nMatchedEcalRecHits = 0;
    unsigned int nMatchedHBHERecHits = 0;
    unsigned int nMatchedHFRecHits = 0;
    unsigned int nMatchedHORecHits = 0;
    unsigned int nMatchedZDCRecHits = 0;
    double matchedEcalEnergy = 0.;
    double matchedHBHEEnergy = 0.;
    double matchedHFEnergy = 0.;
    double matchedHOEnergy = 0.;
    double matchedZDCEnergy = 0.;
    double matchedEcalTime = 0.;
    double matchedHBHETime = 0.;
    double matchedHFTime = 0.;
    double matchedHOTime = 0.;
    double matchedZDCTime = 0.;
    double minEcalLineDistance = -1.;
    double minHBHELineDistance = -1.;
    double minHFLineDistance = -1.;
    double minHOLineDistance = -1.;
    double minZDCLineDistance = -1.;
    int caloTimingDirectionSign = 0;
    unsigned int nCaloTimingMeasurements = 0;
    double caloTimingDeltaChi2 = 0.;
    unsigned int nAddedDTRefitHits = 0;
    unsigned int nAddedTrackerRefitHits = 0;
    std::set<uint32_t> addedDTChamberIds;
    double precisionRefitLeverArm = 0.;
    bool directionalRefitUsedPrecisionHits = false;
    bool directionalRefitAllHitsValid = false;
    unsigned int directionalRefitAllHits = 0;
    double directionalRefitAllHitsChi2 = 0.;
    double directionalRefitAllHitsNdof = 0.;
    double directionalRefitAllHitsTargetPt = 0.;
    double directionalRefitAllHitsTargetPz = 0.;
    int directionalRefitAllHitsSelectedIteration = 0;
    bool directionalRefitPrecisionValid = false;
    unsigned int directionalRefitPrecisionHits = 0;
    double directionalRefitPrecisionChi2 = 0.;
    double directionalRefitPrecisionNdof = 0.;
    double directionalRefitPrecisionUpstreamPt = 0.;
    double directionalRefitPrecisionQoverP = 0.;
    double directionalRefitPrecisionTargetPt = 0.;
    double directionalRefitPrecisionTargetPz = 0.;
    bool directionalRefitPrecisionSecondValid = false;
    bool directionalRefitPrecisionSecondConverged = false;
    double directionalRefitPrecisionRelativeQoverPChange = -1.;
    int directionalRefitPrecisionSelectedIteration = 0;
    double directionalRefitPrecisionFirstTargetPt = 0.;
    double directionalRefitPrecisionFirstTargetPz = 0.;
    double directionalRefitPrecisionSecondTargetPt = 0.;
    double directionalRefitPrecisionSecondTargetPz = 0.;
    double directionalRefitPrecisionRelativeToAllQoverP = -1.;
    double directionalRefitPrecisionTargetDca = -1.;
    bool directionalRefitSourceLegPrecisionValid = false;
    unsigned int directionalRefitSourceLegPrecisionHits = 0;
    double directionalRefitSourceLegPrecisionQoverP = 0.;
    double directionalRefitSourceLegPrecisionTargetPt = 0.;
    double directionalRefitSourceLegPrecisionTargetPz = 0.;
    bool directionalRefitOppositeLegPrecisionValid = false;
    unsigned int directionalRefitOppositeLegPrecisionHits = 0;
    double directionalRefitOppositeLegPrecisionQoverP = 0.;
    double directionalRefitOppositeLegPrecisionTargetPt = 0.;
    double directionalRefitOppositeLegPrecisionTargetPz = 0.;
    bool constrainedValid = false;
    unsigned int constrainedHits = 0;
    double constrainedChi2 = 0.;
    double constrainedNdof = 0.;
    double constrainedTargetChi2 = -1.;
    int constrainedStatus = 0;
    PropagatedState constrainedState;
  };

  struct MuonHitTopology {
    unsigned int dt = 0;
    unsigned int csc = 0;
    unsigned int rpc = 0;
    unsigned int gem = 0;
    unsigned int precisionStations = 0;
  };

  // Geometry-only description of the valid muon-system measurements on the
  // selected track. Region codes are oriented with the reconstructed source
  // side: 0=near endcap, 1=barrel, 2=far endcap, -1=unknown. This deliberately
  // does not use generator information or decide which reconstruction is best.
  struct MuonGeometryDiagnostics {
    unsigned int dt = 0;
    unsigned int csc = 0;
    unsigned int rpc = 0;
    unsigned int gem = 0;
    unsigned int me0 = 0;
    unsigned int plusEndcap = 0;
    unsigned int barrel = 0;
    unsigned int minusEndcap = 0;
    unsigned int nearEndcap = 0;
    unsigned int farEndcap = 0;
    unsigned int nearEndcapStations = 0;
    unsigned int barrelStations = 0;
    unsigned int farEndcapStations = 0;
    unsigned int detectorMask = 0;
    bool endpointsValid = false;
    int entryRegion = -1;
    int exitRegion = -1;
    int entrySubdetector = 0;
    int exitSubdetector = 0;
    GlobalPoint entryPosition;
    GlobalPoint exitPosition;
    double hitSpan = 0.;
  };

  // Public topology values describe recorded muon-system coverage only.  They
  // deliberately do not encode which reconstruction produced the track.
  enum class MuonTopology : int {
    nearEndcapOnly = 0,
    nearEndcapAndBarrel = 1,
    bothEndcaps = 2,
    farEndcapOnly = 3,
    unclassified = 4,
  };

  int measuredMuonTopology(MuonGeometryDiagnostics const& diagnostics, bool orientationValid) {
    if (!orientationValid)
      return static_cast<int>(MuonTopology::unclassified);
    bool const near = diagnostics.nearEndcap > 0;
    bool const barrel = diagnostics.barrel > 0;
    bool const far = diagnostics.farEndcap > 0;
    if (near && far)
      return static_cast<int>(MuonTopology::bothEndcaps);
    if (near && barrel)
      return static_cast<int>(MuonTopology::nearEndcapAndBarrel);
    if (near)
      return static_cast<int>(MuonTopology::nearEndcapOnly);
    if (far && !barrel)
      return static_cast<int>(MuonTopology::farEndcapOnly);
    return static_cast<int>(MuonTopology::unclassified);
  }

  int muonHitSide(DetId const& id) {
    switch (id.subdetId()) {
      case MuonSubdetId::DT:
        return 0;
      case MuonSubdetId::CSC:
        return CSCDetId(id).endcap() == 1 ? 1 : -1;
      case MuonSubdetId::RPC:
        return RPCDetId(id).region();
      case MuonSubdetId::GEM:
        return GEMDetId(id).region();
      case MuonSubdetId::ME0:
        return ME0DetId(id).region();
      default:
        return 2;
    }
  }

  int muonHitStation(DetId const& id) {
    switch (id.subdetId()) {
      case MuonSubdetId::DT:
        return DTChamberId(id).station();
      case MuonSubdetId::CSC:
        return CSCDetId(id).station();
      case MuonSubdetId::RPC:
        return RPCDetId(id).station();
      case MuonSubdetId::GEM:
        return GEMDetId(id).station();
      case MuonSubdetId::ME0:
        return ME0DetId(id).station();
      default:
        return -1;
    }
  }

  int orientedMuonRegion(DetId const& id, int sourceSide) {
    int const side = muonHitSide(id);
    if (side == 0)
      return 1;
    if (side == sourceSide)
      return 0;
    if (side == -sourceSide)
      return 2;
    return -1;
  }

  MuonGeometryDiagnostics muonGeometryDiagnostics(reco::Track const& track,
                                                  GlobalTrackingGeometry const& geometry,
                                                  int sourceSide,
                                                  double directionSign) {
    MuonGeometryDiagnostics result;
    std::array<std::set<unsigned int>, 3> stations;
    auto const endpointDelta = track.outerPosition() - track.innerPosition();
    bool const entryIsOuter = directionSign * track.innerMomentum().Dot(endpointDelta) < 0.;
    auto const entryTrackPosition = entryIsOuter ? track.outerPosition() : track.innerPosition();
    auto const exitTrackPosition = entryIsOuter ? track.innerPosition() : track.outerPosition();
    GlobalPoint const entryReference(entryTrackPosition.x(), entryTrackPosition.y(), entryTrackPosition.z());
    GlobalPoint const exitReference(exitTrackPosition.x(), exitTrackPosition.y(), exitTrackPosition.z());
    GlobalVector measurementDirection = exitReference - entryReference;
    if (!(measurementDirection.mag2() > 0.))
      measurementDirection = GlobalVector(0., 0., -sourceSide);
    measurementDirection = measurementDirection.unit();
    double minimumProjection = std::numeric_limits<double>::infinity();
    double maximumProjection = -std::numeric_limits<double>::infinity();

    for (auto hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
      if (!(*hit)->isValid() || (*hit)->geographicalId().det() != DetId::Muon)
        continue;
      DetId const id = (*hit)->geographicalId();
      switch (id.subdetId()) {
        case MuonSubdetId::DT:
          ++result.dt;
          break;
        case MuonSubdetId::CSC:
          ++result.csc;
          break;
        case MuonSubdetId::RPC:
          ++result.rpc;
          break;
        case MuonSubdetId::GEM:
          ++result.gem;
          break;
        case MuonSubdetId::ME0:
          ++result.me0;
          break;
        default:
          continue;
      }
      result.detectorMask |= 1u << (id.subdetId() - 1);

      int const side = muonHitSide(id);
      result.plusEndcap += side == 1;
      result.barrel += side == 0;
      result.minusEndcap += side == -1;
      result.nearEndcap += side == sourceSide;
      result.farEndcap += side == -sourceSide;
      int const region = orientedMuonRegion(id, sourceSide);
      int const station = muonHitStation(id);
      if (region >= 0 && station >= 0)
        stations[region].insert(1000u * id.subdetId() + 100u * (side + 1) + station);

      auto const* detector = geometry.idToDet(id);
      if (!detector)
        continue;
      GlobalPoint const position = detector->surface().toGlobal((*hit)->localPosition());
      double const projection = (position - entryReference).dot(measurementDirection);
      if (projection < minimumProjection) {
        minimumProjection = projection;
        result.entryPosition = position;
        result.entryRegion = region;
        result.entrySubdetector = id.subdetId();
      }
      if (projection > maximumProjection) {
        maximumProjection = projection;
        result.exitPosition = position;
        result.exitRegion = region;
        result.exitSubdetector = id.subdetId();
      }
    }
    result.nearEndcapStations = stations[0].size();
    result.barrelStations = stations[1].size();
    result.farEndcapStations = stations[2].size();
    result.endpointsValid = std::isfinite(minimumProjection) && std::isfinite(maximumProjection);
    if (result.endpointsValid)
      result.hitSpan = (result.exitPosition - result.entryPosition).mag();
    return result;
  }

  bool isPrecisionMuonSubdetector(int subdetector) {
    return subdetector == MuonSubdetId::DT || subdetector == MuonSubdetId::CSC ||
           subdetector == MuonSubdetId::GEM;
  }

  unsigned int precisionStationKey(DetId const& id) {
    switch (id.subdetId()) {
      case MuonSubdetId::DT:
        return 100 + DTChamberId(id).station();
      case MuonSubdetId::CSC: {
        CSCDetId const csc(id);
        return 200 + 10 * csc.endcap() + csc.station();
      }
      case MuonSubdetId::GEM: {
        GEMDetId const gem(id);
        return 400 + 10 * (gem.region() + 1) + gem.station();
      }
      default:
        return 0;
    }
  }

  MuonHitTopology muonHitTopology(reco::Track const& track) {
    MuonHitTopology result;
    std::set<unsigned int> stations;
    for (auto hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
      if (!(*hit)->isValid() || (*hit)->geographicalId().det() != DetId::Muon)
        continue;
      auto const id = (*hit)->geographicalId();
      switch (id.subdetId()) {
        case MuonSubdetId::DT:
          ++result.dt;
          break;
        case MuonSubdetId::CSC:
          ++result.csc;
          break;
        case MuonSubdetId::RPC:
          ++result.rpc;
          break;
        case MuonSubdetId::GEM:
          ++result.gem;
          break;
      }
      auto const station = precisionStationKey(id);
      if (station != 0)
        stations.insert(station);
    }
    result.precisionStations = stations.size();
    return result;
  }

  PropagatedState propagateStateToTargetLine(FreeTrajectoryState const& start,
                                             Propagator const& vacuumPropagator,
                                             Propagator const* materialPropagator = nullptr,
                                             int sourceSide = 0) {
    auto const position = start.position();
    // Geant4e uses the detailed CMS detector geometry and material.  Stop its
    // transport just beyond the last endcap structures, while still inside
    // the Geant4 world, then continue through the evacuated SHIFT-to-CMS
    // flight path with the no-material propagator.
    constexpr double materialBoundaryZ = 1100.;
    FreeTrajectoryState vacuumStart = start;
    double path = 0.;
    bool materialBoundaryValid = false;
    GlobalVector materialBoundaryMomentum;
    double materialPath = 0.;
    if (materialPropagator && sourceSide != 0 && sourceSide * position.z() < materialBoundaryZ) {
      auto const boundary = Plane::build(GlobalPoint(0., 0., sourceSide * materialBoundaryZ),
                                         Surface::RotationType());
      auto const toBoundary = materialPropagator->propagateWithPath(start, *boundary);
      if (!toBoundary.first.isValid() || !toBoundary.first.freeState()) {
        edm::LogWarning("ShiftMuonMaterialPropagation")
            << "Detailed-material propagation failed before the CMS boundary: position=" << position
            << " momentum=" << start.momentum() << " charge=" << start.charge() << " sourceSide=" << sourceSide
            << " boundaryZ=" << sourceSide * materialBoundaryZ;
        return {};
      }
      vacuumStart = *toBoundary.first.freeState();
      path += toBoundary.second;
      materialBoundaryValid = true;
      materialBoundaryMomentum = vacuumStart.momentum();
      materialPath = toBoundary.second;
    }

    auto const propagated = vacuumPropagator.propagateWithPath(
        vacuumStart, GlobalPoint(0., 0., -1.), GlobalPoint(0., 0., 1.));
    // A failed FreeTrajectoryState propagation is returned default-constructed
    // (zero charge and a null field pointer), unlike a TSOS it has no
    // isValid() accessor.  Check its sentinel before accessing position.
    if (propagated.first.charge() == 0)
      return {};
    auto const resultPosition = propagated.first.position();
    auto const resultMomentum = propagated.first.momentum();
    bool const valid = resultMomentum.mag2() > 0. && std::isfinite(resultPosition.x()) &&
                       std::isfinite(resultPosition.y()) && std::isfinite(resultPosition.z()) &&
                       std::isfinite(path + propagated.second);
    PropagatedState result;
    result.valid = valid;
    result.position = resultPosition;
    result.momentum = resultMomentum;
    result.path = path + propagated.second;
    result.materialBoundaryValid = materialBoundaryValid;
    result.materialBoundaryMomentum = materialBoundaryMomentum;
    result.materialPath = materialPath;
    return result;
  }

  PropagatedState propagateToTargetLine(reco::Track const& track,
                                        Propagator const& vacuumPropagator,
                                        int travelSign = 0,
                                        Propagator const* materialPropagator = nullptr,
                                        int sourceSide = 0) {
    // Once the event side is known, start from the fitted endpoint closest in
    // z to the external source.  This avoids beginning the backward transport
    // from a state that has already crossed additional detector material.
    auto const endpointDelta = track.outerPosition() - track.innerPosition();
    bool const useOuter = travelSign != 0 && travelSign * track.innerMomentum().Dot(endpointDelta) < 0.;
    auto const original = useOuter ? trajectoryStateTransform::outerFreeState(track, vacuumPropagator.magneticField())
                                   : trajectoryStateTransform::innerFreeState(track, vacuumPropagator.magneticField());
    if (!original.hasError())
      return {};
    // Momentum and charge must be reversed together to describe the same
    // fitted helix with the physical time orientation.  Since the external
    // target lies behind the incoming particle, the any-direction propagator
    // then selects the opposite-to-momentum geometrical solution.
    double const sign = travelSign == 0 ? 1. : travelSign;
    auto const endpointMomentum = original.momentum();
    GlobalVector const momentum(sign * endpointMomentum.x(), sign * endpointMomentum.y(), sign * endpointMomentum.z());
    GlobalTrajectoryParameters const parameters(
        original.position(), momentum, sign * original.charge(), vacuumPropagator.magneticField());
    FreeTrajectoryState const physicalState(parameters, original.curvilinearError());
    return propagateStateToTargetLine(physicalState, vacuumPropagator, materialPropagator, sourceSide);
  }

  struct RefitIterationResult {
    bool valid = false;
    int status = 0;
    unsigned int hits = 0;
    double chi2 = 0.;
    double ndof = 0.;
    double upstreamPt = 0.;
    double upstreamP = 0.;
    double downstreamP = 0.;
    double signedInverseMomentum = 0.;
    // The KF filter runs from the external source toward CMS.  Keep the
    // smoother endpoint convention explicit instead of relying on the
    // historical "last measurement" assumption when interpreting its state.
    double sourceFacingSurfaceDistance = -1.;
    double filterDirectionCosine = -2.;
    double targetDirectionCosine = -2.;
    TrajectoryStateOnSurface upstreamState;
    PropagatedState materialTargetState;
    PropagatedState vacuumTargetState;
  };

  struct TargetConstraint {
    GlobalPoint position;
    double sigmaX = 0.;
    double sigmaY = 0.;
    double sigmaZ = 0.;
  };

  struct TargetConstraintResult {
    bool valid = false;
    int status = 0;
    double chi2 = -1.;
    PropagatedState state;
  };

  TargetConstraintResult applyTargetConstraint(TrajectoryStateOnSurface const& upstream,
                                               int sourceSide,
                                               Propagator const& vacuumPropagator,
                                               Propagator const& materialPropagator,
                                               TargetConstraint const& constraint) {
    TargetConstraintResult result;
    if (!upstream.isValid() || !upstream.freeState()) {
      result.status = -1;
      return result;
    }

    FreeTrajectoryState transported = *upstream.freeState();
    constexpr double materialBoundaryZ = 1100.;
    if (sourceSide != 0 && sourceSide * transported.position().z() < materialBoundaryZ) {
      auto const boundary = Plane::build(GlobalPoint(0., 0., sourceSide * materialBoundaryZ), Surface::RotationType());
      auto const toBoundary = materialPropagator.propagate(transported, *boundary);
      if (!toBoundary.isValid() || !toBoundary.freeState()) {
        result.status = -2;
        return result;
      }
      transported = *toBoundary.freeState();
    }

    auto const targetPlane = Plane::build(constraint.position, Surface::RotationType());
    auto const predicted = vacuumPropagator.propagate(transported, *targetPlane);
    if (!predicted.isValid() || !predicted.freeState()) {
      result.status = -3;
      return result;
    }
    auto const momentum = predicted.globalMomentum();
    if (!(constraint.sigmaX > 0.) || !(constraint.sigmaY > 0.) || !(constraint.sigmaZ >= 0.) ||
        std::abs(momentum.z()) < 1.e-9) {
      result.status = -4;
      return result;
    }
    double const slopeX = momentum.x() / momentum.z();
    double const slopeY = momentum.y() / momentum.z();
    double const sigmaZ2 = constraint.sigmaZ * constraint.sigmaZ;
    double const hitXX = constraint.sigmaX * constraint.sigmaX + slopeX * slopeX * sigmaZ2;
    double const hitYY = constraint.sigmaY * constraint.sigmaY + slopeY * slopeY * sigmaZ2;
    double const hitXY = slopeX * slopeY * sigmaZ2;
    if (!(hitXX > 0.) || !(hitYY > 0.) || !std::isfinite(hitXX) || !std::isfinite(hitYY) ||
        !std::isfinite(hitXY)) {
      result.status = -4;
      return result;
    }

    auto const targetHit =
        TRecHit2DPosConstraint::build(LocalPoint(0., 0., 0.), LocalError(hitXX, hitXY, hitYY), targetPlane.get());
    auto const residual = predicted.localPosition() - targetHit->localPosition();
    auto const predictedError = predicted.localError().positionError();
    double const xx = predictedError.xx() + hitXX;
    double const xy = predictedError.xy() + hitXY;
    double const yy = predictedError.yy() + hitYY;
    double const determinant = xx * yy - xy * xy;
    if (!(determinant > 0.) || !std::isfinite(determinant)) {
      result.status = -5;
      return result;
    }
    result.chi2 = (yy * residual.x() * residual.x() - 2. * xy * residual.x() * residual.y() +
                   xx * residual.y() * residual.y()) /
                  determinant;
    KFUpdator updator;
    auto const constrained = updator.update(predicted, *targetHit);
    if (!constrained.isValid() || !constrained.freeState() || !std::isfinite(result.chi2)) {
      result.status = -6;
      return result;
    }
    result.valid = true;
    result.status = 1;
    result.state.valid = true;
    result.state.position = constrained.globalPosition();
    result.state.momentum = constrained.globalMomentum();
    return result;
  }

  struct DirectionalRefitResult {
    bool valid = false;
    int status = 0;
    RefitIterationResult first;
    RefitIterationResult second;
    bool secondConverged = false;
    double relativeQoverPChange = -1.;
    int selectedIteration = 0;
    RefitIterationResult selected;
    unsigned int inputHits = 0;
    unsigned int inputStations = 0;
    double inputLeverArm = 0.;
    unsigned int pathOrderingFallbackHits = 0;
    double pathOrderingSpan = 0.;
  };

  FreeTrajectoryState inflateCurvatureError(FreeTrajectoryState const& state, double scale) {
    auto covariance = state.curvilinearError().matrix();
    // Curvilinear parameter zero is signed inverse momentum.  Inflate its
    // uncertainty, and its covariances so correlation coefficients remain
    // unchanged, without weakening the already useful angular/position seed.
    covariance(0, 0) *= scale * scale;
    for (unsigned int index = 1; index < 5; ++index)
      covariance(index, 0) *= scale;
    return FreeTrajectoryState(state.parameters(), CurvilinearTrajectoryError(covariance));
  }

  FreeTrajectoryState inflateSeedError(FreeTrajectoryState const& state, double scale, bool rescaleFullError) {
    if (!rescaleFullError)
      return inflateCurvatureError(state, scale);
    auto inflated = state;
    inflated.rescaleError(scale);
    return inflated;
  }

  bool finiteTrajectoryState(TrajectoryStateOnSurface const& state) {
    if (!state.isValid())
      return false;
    auto const& parameters = state.localParameters().vector();
    for (int index = 0; index < 5; ++index)
      if (!std::isfinite(parameters[index]))
        return false;
    auto const& error = state.curvilinearError();
    if (!error.posDef())
      return false;
    auto const& covariance = error.matrix();
    for (int row = 0; row < 5; ++row)
      for (int column = 0; column <= row; ++column)
        if (!std::isfinite(covariance(row, column)))
          return false;
    return true;
  }

  bool finiteTrajectory(Trajectory const& trajectory) {
    if (!trajectory.isValid() || !std::isfinite(trajectory.chiSquared()))
      return false;
    for (auto const& measurement : trajectory.measurements()) {
      if (!std::isfinite(measurement.estimate()) ||
          !finiteTrajectoryState(measurement.forwardPredictedState()) ||
          !finiteTrajectoryState(measurement.updatedState()))
        return false;
    }
    return true;
  }

  DirectionalRefitResult directionalRefit(reco::Track const& track,
                                          int directionSign,
                                          int sourceSide,
                                          MagneticField const& magneticField,
                                          TransientTrackingRecHitBuilder const& hitBuilder,
                                          TkCloner const& hitCloner,
                                          Propagator const& vacuumPropagator,
                                          Propagator const& fitMaterialPropagator,
                                          Propagator const& smootherMaterialPropagator,
                                          Propagator const& targetMaterialPropagator,
                                          double seedCurvatureErrorRescale,
                                          double secondSeedErrorRescale,
                                          double seedMomentumScale,
                                          bool useFullSeedErrorRescale,
                                          double smootherErrorRescale,
                                          double initialMaxHitChi2,
                                          double maxHitChi2,
                                          double maxRelativeQoverPChange,
                                          bool useSecondIteration,
                                          bool precisionHitsOnly = false,
                                          bool usePathOrdering = true,
                                          int hitSideSelection = 0,
                                          std::vector<std::pair<TransientTrackingRecHit::RecHitPointer, bool>> const* extraHits = nullptr) {
    struct OrderedHit {
      TransientTrackingRecHit::RecHitPointer hit;
      double path;
    };
    std::vector<OrderedHit> orderedHits;
    for (auto hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
      if (!(*hit)->isValid() || (*hit)->geographicalId().det() != DetId::Muon)
        continue;
      if (precisionHitsOnly && !isPrecisionMuonSubdetector((*hit)->geographicalId().subdetId()))
        continue;
      auto transientHit = hitBuilder.build(&**hit);
      if (!transientHit || !transientHit->isValid() || !transientHit->det())
        continue;
      // Optional diagnostic split for full-lever-arm tracks. +1 keeps the
      // source-facing detector half and -1 the opposite half. The canonical
      // fit always uses zero and therefore remains unchanged.
      if (hitSideSelection != 0 &&
          hitSideSelection * sourceSide * transientHit->globalPosition().z() <= 0.)
        continue;
      orderedHits.push_back({transientHit, usePathOrdering ? 0. : sourceSide * transientHit->globalPosition().z()});
    }
    if (extraHits)
      for (auto const& [transientHit, precision] : *extraHits) {
        if (!transientHit || !transientHit->isValid() || !transientHit->det() ||
            (precisionHitsOnly && !precision))
          continue;
        orderedHits.push_back({transientHit, usePathOrdering ? 0. : sourceSide * transientHit->globalPosition().z()});
      }

    DirectionalRefitResult result;
    result.inputHits = orderedHits.size();
    std::set<unsigned int> inputStations;
    for (auto const& orderedHit : orderedHits) {
      auto const station = precisionStationKey(orderedHit.hit->geographicalId());
      if (station != 0)
        inputStations.insert(station);
    }
    result.inputStations = inputStations.size();
    if (orderedHits.size() < 3 || !(seedMomentumScale > 0.) || !std::isfinite(seedMomentumScale)) {
      result.status = -1;
      return result;
    }

    bool const useOuter = sourceSide * track.outerPosition().z() > sourceSide * track.innerPosition().z();
    auto const original = useOuter ? trajectoryStateTransform::outerFreeState(track, &magneticField)
                                   : trajectoryStateTransform::innerFreeState(track, &magneticField);
    if (!original.hasError()) {
      result.status = -2;
      return result;
    }
    double const sign = directionSign == 0 ? 1. : directionSign;
    auto const rawMomentum = original.momentum();
    GlobalVector const momentum(sign * seedMomentumScale * rawMomentum.x(),
                                sign * seedMomentumScale * rawMomentum.y(),
                                sign * seedMomentumScale * rawMomentum.z());
    GlobalTrajectoryParameters const parameters(
        original.position(), momentum, sign * original.charge(), &magneticField);
    FreeTrajectoryState const originalSeed(parameters, original.curvilinearError());
    FreeTrajectoryState iterationSeed = originalSeed;
    if (usePathOrdering) {
      // Order measurements by actual propagated path from the source-facing
      // endpoint.  Global z is not a valid trajectory coordinate for barrel,
      // bending, or traversing tracks.  If the propagation to a measurement
      // surface fails, retain that hit using its signed projection along the
      // seed direction and expose the fallback count as a diagnostic.
      auto const seedDirection = momentum.unit();
      for (auto& orderedHit : orderedHits) {
        auto const propagated = vacuumPropagator.propagateWithPath(originalSeed, *orderedHit.hit->surface());
        if (propagated.first.isValid() && std::isfinite(propagated.second)) {
          orderedHit.path = propagated.second;
        } else {
          auto const displacement = orderedHit.hit->globalPosition() - original.position();
          orderedHit.path = displacement.dot(seedDirection);
          ++result.pathOrderingFallbackHits;
        }
      }
    }
    std::stable_sort(orderedHits.begin(), orderedHits.end(), [usePathOrdering](OrderedHit const& first,
                                                                              OrderedHit const& second) {
      return usePathOrdering ? first.path < second.path : first.path > second.path;
    });
    result.pathOrderingSpan = orderedHits.back().path - orderedHits.front().path;
    result.inputLeverArm =
        (orderedHits.front().hit->globalPosition() - orderedHits.back().hit->globalPosition()).mag();

    Trajectory::RecHitContainer fitHits;
    fitHits.reserve(orderedHits.size());
    for (auto const& orderedHit : orderedHits)
      fitHits.push_back(orderedHit.hit);

    KFUpdator updator;
    auto runIteration = [&](FreeTrajectoryState const& uninflatedSeed,
                            double iterationSeedErrorRescale,
                            double iterationMaxHitChi2) {
      RefitIterationResult result;
      Chi2MeasurementEstimator estimator(iterationMaxHitChi2);
      KFTrajectoryFitter fitter(fitMaterialPropagator, updator, estimator, 3, nullptr, &hitCloner);
      KFTrajectorySmoother smoother(
          smootherMaterialPropagator, updator, estimator, static_cast<float>(smootherErrorRescale), 3);
      smoother.setHitCloner(&hitCloner);
      auto const start = inflateSeedError(uninflatedSeed, iterationSeedErrorRescale, useFullSeedErrorRescale);
      TrajectorySeed const seed(
          PTrajectoryStateOnDet(), TrajectorySeed::RecHitContainer(), alongMomentum);
      int const minimumFoundHits = 3;
      auto fitAndSmooth = [&](Trajectory::RecHitContainer const& iterationHits, int& status) {
        Trajectory invalid;
        auto const firstPredicted = fitMaterialPropagator.propagate(start, *iterationHits.front()->surface());
        if (!firstPredicted.isValid()) {
          status = 10;
          return invalid;
        }
        auto const filtered = fitter.fitOne(seed, iterationHits, firstPredicted, TrajectoryFitter::standard);
        // KFTrajectoryFitter deliberately has a no-fail policy: after a failed
        // propagation/update it can return a trajectory with a poisoned final
        // state. Reject it before the smoother can enter Geant4e.
        if (filtered.foundHits() < minimumFoundHits) {
          status = 11;
          return invalid;
        }
        if (!finiteTrajectory(filtered)) {
          status = 12;
          return invalid;
        }
        auto smoothedCandidate = smoother.trajectory(filtered);
        if (smoothedCandidate.foundHits() < minimumFoundHits || smoothedCandidate.empty()) {
          status = 13;
          return invalid;
        }
        if (!finiteTrajectory(smoothedCandidate)) {
          status = 14;
          return invalid;
        }
        status = 1;
        return smoothedCandidate;
      };

      // KFTrajectoryFitter's estimator records hit compatibility but does not
      // reject incompatible measurements. Implement the same remove-and-refit
      // principle as KFFittingSmoother: remove at most two hits and never more
      // than 20% of the input, choosing the candidate with the smallest total
      // smoothed chi2. The canonical EstimateCut=20 used by CMSSW material
      // refitters is supplied through iterationMaxHitChi2.
      Trajectory::RecHitContainer iterationHits = fitHits;
      int fitStatus = 0;
      auto smoothed = fitAndSmooth(iterationHits, fitStatus);
      if (!smoothed.isValid()) {
        result.status = fitStatus;
        return result;
      }
      unsigned int const maximumOutliers =
          std::min<unsigned int>(2, static_cast<unsigned int>(fitHits.size() / 5));
      unsigned int removedOutliers = 0;
      while (iterationMaxHitChi2 > 0. && removedOutliers < maximumOutliers &&
             smoothed.measurements().size() == iterationHits.size()) {
        std::vector<unsigned int> badMeasurements;
        for (unsigned int index = 0; index < smoothed.measurements().size(); ++index) {
          auto const& measurement = smoothed.measurements()[index];
          if (measurement.recHitR().isValid() && measurement.recHitR().det() &&
              measurement.estimate() > iterationMaxHitChi2)
            badMeasurements.push_back(index);
        }
        if (badMeasurements.empty())
          break;

        bool foundCandidate = false;
        double bestChi2 = std::numeric_limits<double>::infinity();
        unsigned int bestHitIndex = 0;
        Trajectory bestTrajectory;
        for (auto const measurementIndex : badMeasurements) {
          unsigned int const hitIndex = iterationHits.size() - measurementIndex - 1;
          auto candidateHits = iterationHits;
          candidateHits.erase(candidateHits.begin() + hitIndex);
          if (candidateHits.size() < static_cast<unsigned int>(minimumFoundHits))
            continue;
          int candidateStatus = 0;
          auto candidate = fitAndSmooth(candidateHits, candidateStatus);
          if (candidate.isValid() && candidate.chiSquared() < bestChi2) {
            foundCandidate = true;
            bestChi2 = candidate.chiSquared();
            bestHitIndex = hitIndex;
            bestTrajectory = std::move(candidate);
          }
        }
        if (!foundCandidate || !(bestChi2 < smoothed.chiSquared()))
          break;
        iterationHits.erase(iterationHits.begin() + bestHitIndex);
        smoothed = std::move(bestTrajectory);
        ++removedOutliers;
      }

      // KFTrajectorySmoother returns measurements in reverse filter order.
      // Since this fit was seeded source -> CMS, lastMeasurement is therefore
      // the source-facing state.  Record both surface and direction closure
      // quantities here: they make a convention mismatch visible without
      // hiding it with a momentum calibration.
      auto const* upstream = &smoothed.lastMeasurement().updatedState();
      if (!upstream->isValid() || !upstream->freeState()) {
        result.status = 15;
        return result;
      }
      auto const& downstream = smoothed.firstMeasurement().updatedState();
      if (!downstream.isValid()) {
        result.status = 16;
        return result;
      }
      auto const upstreamMomentum = upstream->globalMomentum();
      auto const downstreamMomentum = downstream.globalMomentum();
      // A formally valid Kalman state with effectively zero momentum is not
      // transportable over the CMS-to-target geometry. Reject it before the
      // propagator; this is particularly important for diagnostic half-track
      // fits with too little curvature lever arm.
      if (!(upstreamMomentum.mag() > 0.1) || !(downstreamMomentum.mag() > 0.1) ||
          !std::isfinite(upstreamMomentum.mag()) || !std::isfinite(downstreamMomentum.mag())) {
        result.status = 17;
        return result;
      }
      auto const sourceFacingPosition = upstream->globalPosition();
      auto const downstreamPosition = downstream.globalPosition();
      result.sourceFacingSurfaceDistance =
          (sourceFacingPosition - iterationHits.front()->globalPosition()).mag();
      auto const filterDisplacement = downstreamPosition - sourceFacingPosition;
      if (filterDisplacement.mag2() > 0.)
        result.filterDirectionCosine = filterDisplacement.dot(upstreamMomentum) /
                                       (filterDisplacement.mag() * upstreamMomentum.mag());
      // The source-facing state must point away from the external target.  A
      // transport back to that target is consequently opposite-to-momentum.
      // A non-negative value indicates that the state was selected with the
      // wrong physical orientation, so do not publish it as a valid refit.
      result.targetDirectionCosine = sourceSide == 0
                                         ? -2.
                                         : sourceSide * upstreamMomentum.z() / upstreamMomentum.mag();
      if (sourceSide != 0 && (!std::isfinite(result.targetDirectionCosine) ||
                              result.targetDirectionCosine >= 0.)) {
        result.status = 18;
        return result;
      }
      auto const upstreamFreeState = *upstream->freeState();
      PropagatedState materialTargetState;
      PropagatedState vacuumTargetState;
      materialTargetState =
          propagateStateToTargetLine(upstreamFreeState, vacuumPropagator, &targetMaterialPropagator, sourceSide);
      vacuumTargetState = propagateStateToTargetLine(upstreamFreeState, vacuumPropagator, nullptr, sourceSide);
      if (!materialTargetState.valid) {
        result.status = 17;
        return result;
      }
      auto const refitHits = static_cast<unsigned int>(smoothed.foundHits());
      result.valid = true;
      result.status = 1;
      result.hits = refitHits;
      result.chi2 = smoothed.chiSquared();
      result.ndof = std::max(1., 2. * static_cast<double>(refitHits) - 5.);
      result.upstreamPt = upstreamMomentum.perp();
      result.upstreamP = upstreamMomentum.mag();
      result.downstreamP = downstreamMomentum.mag();
      result.signedInverseMomentum = upstream->signedInverseMomentum();
      result.upstreamState = *upstream;
      result.materialTargetState = materialTargetState;
      result.vacuumTargetState = vacuumTargetState;
      return result;
    };

    result.first = runIteration(iterationSeed, seedCurvatureErrorRescale, initialMaxHitChi2);
    if (!result.first.valid) {
      result.status = 100 + result.first.status;
      return result;
    }

    if (!useSecondIteration) {
      result.selectedIteration = 1;
      result.selected = result.first;
      result.valid = true;
      result.status = 1;
      return result;
    }

    // Pass two reuses pass one's posterior state as a nonlinear seed, but it
    // must not silently reuse the information from the same hits as a strong
    // prior. Keep its covariance inflation independently controllable so the
    // repeated-measurement feedback can be tested without changing pass one.
    result.second =
        runIteration(*result.first.upstreamState.freeState(), secondSeedErrorRescale, maxHitChi2);
    if (result.second.valid) {
      double const denominator = std::max(std::abs(result.first.signedInverseMomentum), 1.e-12);
      result.relativeQoverPChange =
          std::abs(result.second.signedInverseMomentum - result.first.signedInverseMomentum) / denominator;
      result.secondConverged = std::isfinite(result.relativeQoverPChange) &&
                               result.relativeQoverPChange <= maxRelativeQoverPChange;
    }
    result.selectedIteration = result.secondConverged ? 2 : 1;
    result.selected = result.secondConverged ? result.second : result.first;
    result.valid = result.selected.valid;
    result.status = result.valid ? 1 : 200 + result.second.status;
    return result;
  }

  struct TimingResult {
    int directionSign = 0;
    unsigned int measurements = 0;
    double chi2 = 0.;
    double deltaChi2 = 0.;
  };

  TimingResult timingDirection(reco::Track const& track,
                               GlobalTrackingGeometry const& geometry,
                               unsigned int minMeasurements,
                               double minDeltaChi2,
                               bool useExtendedTiming) {
    struct TimedPoint { GlobalPoint position; double time; double sigma; };
    std::vector<TimedPoint> points;
    for (auto hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
      if (!(*hit)->isValid())
        continue;
      double time = 0., sigma = 0.;
      if (auto const* segment = dynamic_cast<CSCSegment const*>(&**hit)) {
        time = segment->time();
        sigma = 7.;
      } else if (auto const* segment = dynamic_cast<DTRecSegment4D const*>(&**hit)) {
        if (!segment->hasPhi() || !segment->phiSegment()->ist0Valid())
          continue;
        time = segment->phiSegment()->t0();
        sigma = 3.;
      } else if (useExtendedTiming) {
        if (auto const* rpcHit = dynamic_cast<RPCRecHit const*>(&**hit)) {
          time = 25. * rpcHit->BunchX();
          sigma = 12.5;
        } else if (auto const* gemHit = dynamic_cast<GEMRecHit const*>(&**hit)) {
          time = 25. * gemHit->BunchX();
          sigma = 12.5;
        } else {
          continue;
        }
      } else {
        continue;
      }
      auto const* detector = geometry.idToDet((*hit)->geographicalId());
      if (!detector || !std::isfinite(time) || std::abs(time) > 1.e4)
        continue;
      points.push_back({detector->surface().toGlobal((*hit)->localPosition()), time, sigma});
    }
    if (points.size() < minMeasurements)
      return {0, static_cast<unsigned int>(points.size()), 0., 0.};

    constexpr double inverseSpeedOfLight = 1. / 29.9792458;  // ns/cm
    auto const axis = track.momentum().unit();
    auto score = [&](double sign) {
      auto coordinate = [&axis](GlobalPoint const& position) {
        return axis.x() * position.x() + axis.y() * position.y() + axis.z() * position.z();
      };
      std::vector<double> offsets;
      offsets.reserve(points.size());
      for (auto const& point : points)
        offsets.push_back(point.time - sign * coordinate(point.position) * inverseSpeedOfLight);
      std::sort(offsets.begin(), offsets.end());
      double const offset = offsets.size() % 2 != 0
                                ? offsets[offsets.size() / 2]
                                : 0.5 * (offsets[offsets.size() / 2 - 1] + offsets[offsets.size() / 2]);
      double robustChi2 = 0.;
      for (auto const& point : points) {
        double const residual = point.time - offset - sign * coordinate(point.position) * inverseSpeedOfLight;
        double const pull = std::abs(residual / point.sigma);
        // Huber loss keeps one mistimed chamber from deciding the direction.
        robustChi2 += pull <= 3. ? pull * pull : 6. * pull - 9.;
      }
      return robustChi2;
    };
    double const alongChi2 = score(1.);
    double const oppositeChi2 = score(-1.);
    double const delta = std::abs(alongChi2 - oppositeChi2);
    int const sign = delta >= minDeltaChi2 ? (alongChi2 < oppositeChi2 ? 1 : -1) : 0;
    return {sign,
            static_cast<unsigned int>(points.size()),
            std::min(alongChi2, oppositeChi2),
            delta};
  }

  void appendHitFingerprints(TrackingRecHit const& hit, std::vector<HitFingerprint>& result) {
    auto const components = hit.recHits();
    if (!components.empty()) {
      for (auto const* component : components)
        if (component && component->isValid())
          appendHitFingerprints(*component, result);
      return;
    }
    auto const position = hit.localPosition();
    result.push_back({hit.rawId(), position.x(), position.y()});
  }

  std::vector<HitFingerprint> hitFingerprints(reco::Track const& track) {
    std::vector<HitFingerprint> result;
    for (auto hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit)
      if ((*hit)->isValid())
        appendHitFingerprints(**hit, result);
    return result;
  }

  double directionAngle(reco::Track const& first, reco::Track const& second) {
    auto const firstDirection = first.momentum().unit();
    auto const secondDirection = second.momentum().unit();
    double const cosine = std::clamp(std::abs(firstDirection.Dot(secondDirection)), 0., 1.);
    return std::acos(cosine);
  }

  double directionDeltaR(GlobalVector const& momentum, reco::GenParticle const& particle) {
    auto deltaPhi = [](double first, double second) {
      double value = std::remainder(first - second, 2. * M_PI);
      return std::abs(value);
    };
    double const pt = momentum.perp();
    double const eta = pt > 0. ? std::asinh(momentum.z() / pt)
                               : std::copysign(std::numeric_limits<double>::infinity(), momentum.z());
    double const phi = std::atan2(momentum.y(), momentum.x());
    double const direct = std::hypot(eta - particle.eta(), deltaPhi(phi, particle.phi()));
    double const reverse = std::hypot(eta + particle.eta(), deltaPhi(phi + M_PI, particle.phi()));
    return std::min(direct, reverse);
  }

  double pointLineDistance(reco::Track::Point const& point, reco::Track const& line) {
    auto const delta = point - line.referencePoint();
    return delta.Cross(line.momentum().unit()).R();
  }

  double symmetricLineDistance(reco::Track const& first, reco::Track const& second) {
    // Comparing detector states avoids merging two distinct tracks which
    // merely intersect at the common SHIFT production region.
    return std::max({pointLineDistance(first.innerPosition(), second),
                     pointLineDistance(first.outerPosition(), second),
                     pointLineDistance(second.innerPosition(), first),
                     pointLineDistance(second.outerPosition(), first)});
  }

  std::pair<double, double> transverseLinePca(reco::Track const& track) {
    auto const& position = track.innerPosition();
    auto const& momentum = track.innerMomentum();
    double const pt2 = momentum.x() * momentum.x() + momentum.y() * momentum.y();
    if (pt2 == 0.)
      return {position.rho(), position.z()};
    double const scale = -(position.x() * momentum.x() + position.y() * momentum.y()) / pt2;
    return {std::hypot(position.x() + scale * momentum.x(), position.y() + scale * momentum.y()),
            position.z() + scale * momentum.z()};
  }

  std::pair<double, double> chordLinePca(reco::Track const& track) {
    // Independent direction diagnostic made only from the fitted endpoint
    // positions.  It is insensitive to the momentum parameterisation used by
    // the standalone Kalman fit and exposes polar-slope biases directly.
    auto const& inner = track.innerPosition();
    auto const delta = track.outerPosition() - inner;
    double const transverse2 = delta.x() * delta.x() + delta.y() * delta.y();
    if (transverse2 == 0.)
      return {inner.rho(), inner.z()};
    double const scale = -(inner.x() * delta.x() + inner.y() * delta.y()) / transverse2;
    return {std::hypot(inner.x() + scale * delta.x(), inner.y() + scale * delta.y()),
            inner.z() + scale * delta.z()};
  }

  double shiftDirectionSign(GlobalVector const& momentum, int sourceSide) {
    return sourceSide * momentum.z() > 0. ? -1. : 1.;
  }

  struct StraightLineApproach {
    bool valid;
    GlobalPoint midpoint;
    double distance;
  };

  struct CommonLineVertex {
    bool valid = false;
    GlobalPoint position;
    std::array<double, 3> error{{0., 0., 0.}};
    double chi2 = 0.;
    double ndof = 3.;
  };

  bool invertSymmetric3x3(std::array<std::array<double, 3>, 3> const& matrix,
                          std::array<std::array<double, 3>, 3>& inverse) {
    double const a = matrix[0][0], b = matrix[0][1], c = matrix[0][2];
    double const d = matrix[1][1], e = matrix[1][2], f = matrix[2][2];
    double const determinant = a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
    if (!std::isfinite(determinant) || std::abs(determinant) < 1.e-18)
      return false;
    inverse = {{{(d * f - e * e) / determinant, (c * e - b * f) / determinant, (b * e - c * d) / determinant},
                {(c * e - b * f) / determinant, (a * f - c * c) / determinant, (b * c - a * e) / determinant},
                {(b * e - c * d) / determinant, (b * c - a * e) / determinant, (a * d - b * b) / determinant}}};
    return true;
  }

  CommonLineVertex commonLineVertex(PropagatedState const& first,
                                    PropagatedState const& second,
                                    double lineResolution,
                                    double beamLineResolution) {
    // Minimise the distances to the two unbounded track lines together with
    // a loose x=y=0 target-line prior.  Unlike TransientTrack propagation,
    // this analytic fit has no detector-volume boundary at |z|~10 m.
    std::array<std::array<double, 3>, 3> normal{};
    std::array<double, 3> rhs{};
    double const lineWeight = 1. / (lineResolution * lineResolution);
    for (auto const* state : {&first, &second}) {
      auto const point = state->position;
      auto const direction = state->momentum.unit();
      std::array<double, 3> const p{{point.x(), point.y(), point.z()}};
      std::array<double, 3> const u{{direction.x(), direction.y(), direction.z()}};
      for (unsigned int row = 0; row < 3; ++row)
        for (unsigned int column = 0; column < 3; ++column) {
          double const projector = (row == column ? 1. : 0.) - u[row] * u[column];
          normal[row][column] += lineWeight * projector;
          rhs[row] += lineWeight * projector * p[column];
        }
    }
    double const beamWeight = 1. / (beamLineResolution * beamLineResolution);
    normal[0][0] += beamWeight;
    normal[1][1] += beamWeight;
    std::array<std::array<double, 3>, 3> covariance{};
    if (!invertSymmetric3x3(normal, covariance))
      return {};
    std::array<double, 3> fitted{};
    for (unsigned int row = 0; row < 3; ++row)
      for (unsigned int column = 0; column < 3; ++column)
        fitted[row] += covariance[row][column] * rhs[column];
    if (!std::isfinite(fitted[0]) || !std::isfinite(fitted[1]) || !std::isfinite(fitted[2]))
      return {};
    GlobalPoint const position(fitted[0], fitted[1], fitted[2]);
    double chi2 = (fitted[0] * fitted[0] + fitted[1] * fitted[1]) * beamWeight;
    for (auto const* state : {&first, &second}) {
      GlobalVector const delta(fitted[0] - state->position.x(),
                               fitted[1] - state->position.y(),
                               fitted[2] - state->position.z());
      chi2 += std::pow(delta.cross(state->momentum.unit()).mag(), 2) * lineWeight;
    }
    return {true,
            position,
            {{std::sqrt(std::max(0., covariance[0][0])),
              std::sqrt(std::max(0., covariance[1][1])),
              std::sqrt(std::max(0., covariance[2][2]))}},
            chi2,
            3.};
  }

  StraightLineApproach straightLineApproach(PropagatedState const& first, PropagatedState const& second) {
    auto const firstPoint = first.position;
    auto const secondPoint = second.position;
    auto const firstDirection = first.momentum.unit();
    auto const secondDirection = second.momentum.unit();
    auto const separation = firstPoint - secondPoint;
    double const dot = firstDirection.dot(secondDirection);
    double const denominator = 1. - dot * dot;
    if (denominator < 1.e-10) {
      GlobalPoint const midpoint(0.5 * (firstPoint.x() + secondPoint.x()),
                                 0.5 * (firstPoint.y() + secondPoint.y()),
                                 0.5 * (firstPoint.z() + secondPoint.z()));
      return {true, midpoint, (firstPoint - secondPoint).mag()};
    }
    double const firstScale = (dot * secondDirection.dot(separation) - firstDirection.dot(separation)) / denominator;
    double const secondScale = (secondDirection.dot(separation) - dot * firstDirection.dot(separation)) / denominator;
    auto const firstClosest = firstPoint + firstScale * firstDirection;
    auto const secondClosest = secondPoint + secondScale * secondDirection;
    GlobalPoint const midpoint(0.5 * (firstClosest.x() + secondClosest.x()),
                               0.5 * (firstClosest.y() + secondClosest.y()),
                               0.5 * (firstClosest.z() + secondClosest.z()));
    double const distance = (firstClosest - secondClosest).mag();
    return {std::isfinite(distance) && std::isfinite(midpoint.z()), midpoint, distance};
  }
}  // namespace

class ShiftMuonTableProducer : public edm::stream::EDProducer<> {
public:
  explicit ShiftMuonTableProducer(edm::ParameterSet const& parameters)
      : dsaToken_(consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("dsaTracks"))),
        cosmicToken_(consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("cosmicTracks"))),
        traversingToken_(consumes<reco::TrackCollection>(parameters.getParameter<edm::InputTag>("traversingTracks"))),
        dsaGlobalLinksToken_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("dsaGlobalLinks"))),
        cosmicGlobalLinksToken_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("cosmicGlobalLinks"))),
        traversingGlobalLinksToken_(consumes<reco::MuonTrackLinksCollection>(
            parameters.getParameter<edm::InputTag>("traversingGlobalLinks"))),
        genParticlesToken_(
            consumes<reco::GenParticleCollection>(parameters.getParameter<edm::InputTag>("genParticles"))),
        simTracksToken_(consumes<edm::SimTrackContainer>(parameters.getParameter<edm::InputTag>("simTracks"))),
        simVerticesToken_(consumes<edm::SimVertexContainer>(parameters.getParameter<edm::InputTag>("simVertices"))),
        dtSimHitsToken_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("dtSimHits"))),
        cscSimHitsToken_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("cscSimHits"))),
        rpcSimHitsToken_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("rpcSimHits"))),
        gemSimHitsToken_(consumes<edm::PSimHitContainer>(parameters.getParameter<edm::InputTag>("gemSimHits"))),
        hcalSimHitsToken_(consumes<edm::PCaloHitContainer>(parameters.getParameter<edm::InputTag>("hcalSimHits"))),
        zdcSimHitsToken_(consumes<edm::PCaloHitContainer>(parameters.getParameter<edm::InputTag>("zdcSimHits"))),
        ecalBarrelRecHitsToken_(consumes<EcalRecHitCollection>(
            parameters.getParameter<edm::InputTag>("ecalBarrelRecHits"))),
        ecalEndcapRecHitsToken_(consumes<EcalRecHitCollection>(
            parameters.getParameter<edm::InputTag>("ecalEndcapRecHits"))),
        hbheRecHitsToken_(consumes<HBHERecHitCollection>(parameters.getParameter<edm::InputTag>("hbheRecHits"))),
        hfRecHitsToken_(consumes<HFRecHitCollection>(parameters.getParameter<edm::InputTag>("hfRecHits"))),
        hoRecHitsToken_(consumes<HORecHitCollection>(parameters.getParameter<edm::InputTag>("hoRecHits"))),
        zdcRecHitsToken_(consumes<ZDCRecHitCollection>(parameters.getParameter<edm::InputTag>("zdcRecHits"))),
        dtSegmentsToken_(consumes<DTRecSegment4DCollection>(parameters.getParameter<edm::InputTag>("dtSegments"))),
        pixelRecHitsToken_(consumes<SiPixelRecHitCollection>(parameters.getParameter<edm::InputTag>("pixelRecHits"))),
        stripMatchedRecHitsToken_(consumes<SiStripMatchedRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripMatchedRecHits"))),
        stripRphiRecHitsToken_(consumes<SiStripRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripRphiRecHits"))),
        stripRphiUnmatchedRecHitsToken_(consumes<SiStripRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripRphiUnmatchedRecHits"))),
        stripStereoRecHitsToken_(consumes<SiStripRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripStereoRecHits"))),
        stripStereoUnmatchedRecHitsToken_(consumes<SiStripRecHit2DCollection>(
            parameters.getParameter<edm::InputTag>("stripStereoUnmatchedRecHits"))),
        magneticFieldToken_(esConsumes()),
        trackingGeometryToken_(esConsumes()),
        muonRecHitBuilderToken_(esConsumes(edm::ESInputTag(
            "", parameters.getParameter<std::string>("muonRecHitBuilder")))),
        trackerRecHitBuilderToken_(esConsumes(edm::ESInputTag(
            "", parameters.getParameter<std::string>("trackerRecHitBuilder")))),
        caloGeometryToken_(esConsumes()),
        augmentDTHits_(parameters.getParameter<bool>("augmentDTHits")),
        augmentTrackerHits_(parameters.getParameter<bool>("augmentTrackerHits")),
        useExtendedTiming_(parameters.getParameter<bool>("useExtendedTiming")),
        additionalDTMaxDistance_(parameters.getParameter<double>("additionalDTMaxDistance")),
        additionalTrackerMaxDistance_(parameters.getParameter<double>("additionalTrackerMaxDistance")),
        additionalHitMaxChi2_(parameters.getParameter<double>("additionalHitMaxChi2")),
        ecalMatchMaxDistance_(parameters.getParameter<double>("ecalMatchMaxDistance")),
        hcalMatchMaxDistance_(parameters.getParameter<double>("hcalMatchMaxDistance")),
        zdcMatchMaxDistance_(parameters.getParameter<double>("zdcMatchMaxDistance")),
        caloMatchMinEnergy_(parameters.getParameter<double>("caloMatchMinEnergy")),
        useImprovedMomentumRefit_(parameters.getParameter<bool>("useImprovedMomentumRefit")),
        useDetailedMaterialPropagation_(parameters.getParameter<bool>("useDetailedMaterialPropagation")),
        directionalRefitUseMaterialEffects_(
            parameters.getParameter<bool>("directionalRefitUseMaterialEffects")),
        directionalRefitUseFirstPrinciplesMaterialEffects_(
            parameters.getParameter<bool>("directionalRefitUseFirstPrinciplesMaterialEffects")),
        directionalRefitFirstPrinciplesStepCm_(
            parameters.getParameter<double>("directionalRefitFirstPrinciplesStepCm")),
        directionalRefitUseDetailedMaterialEffects_(
            parameters.getParameter<bool>("directionalRefitUseDetailedMaterialEffects")),
        directionalRefitUseGeometryMaterialEffects_(
            parameters.getParameter<bool>("directionalRefitUseGeometryMaterialEffects")),
        directionalRefitUseGeometryMaterialEffectsInFitter_(
            parameters.getParameter<bool>("directionalRefitUseGeometryMaterialEffectsInFitter")),
        directionalRefitUseGeometryMaterialEffectsInSmoother_(
            parameters.getParameter<bool>("directionalRefitUseGeometryMaterialEffectsInSmoother")),
        directionalRefitUseGeometryTargetMaterialEffects_(
            parameters.getParameter<bool>("directionalRefitUseGeometryTargetMaterialEffects")),
        directionalRefitGeometryMaterialStepCm_(
            parameters.getParameter<double>("directionalRefitGeometryMaterialStepCm")),
        directionalRefitLogGeometryMaterialComparison_(
            parameters.getParameter<bool>("directionalRefitLogGeometryMaterialComparison")),
        usePropagatedPathOrdering_(parameters.getParameter<bool>("usePropagatedPathOrdering")),
        directionalRefitUseFullSeedErrorRescale_(
            parameters.getParameter<bool>("directionalRefitUseFullSeedErrorRescale")),
        directionalRefitSeedCurvatureErrorRescale_(
            parameters.getParameter<double>("directionalRefitSeedCurvatureErrorRescale")),
        directionalRefitSecondSeedErrorRescale_(
            parameters.getParameter<double>("directionalRefitSecondSeedErrorRescale")),
        directionalRefitSeedMomentumScale_(parameters.getParameter<double>("directionalRefitSeedMomentumScale")),
        directionalRefitEnergyLossScale_(parameters.getParameter<double>("directionalRefitEnergyLossScale")),
        directionalRefitErrorRescale_(parameters.getParameter<double>("directionalRefitErrorRescale")),
        directionalRefitInitialMaxHitChi2_(
            parameters.getParameter<double>("directionalRefitInitialMaxHitChi2")),
        directionalRefitMaxHitChi2_(parameters.getParameter<double>("directionalRefitMaxHitChi2")),
        directionalRefitMaxRelativeQoverPChange_(
            parameters.getParameter<double>("directionalRefitMaxRelativeQoverPChange")),
        directionalRefitUseSecondIteration_(parameters.getParameter<bool>("directionalRefitUseSecondIteration")),
        directionalRefitMinPrecisionStations_(
            parameters.getParameter<unsigned int>("directionalRefitMinPrecisionStations")),
        directionalRefitMaxPrecisionRelativeQoverPChange_(
            parameters.getParameter<double>("directionalRefitMaxPrecisionRelativeQoverPChange")),
        produceMomentumClosureDiagnostics_(parameters.getParameter<bool>("produceMomentumClosureDiagnostics")),
        enableHcalDiagnostics_(parameters.getParameter<bool>("enableHcalDiagnostics")),
        enableZDCDiagnostics_(parameters.getParameter<bool>("enableZDCDiagnostics")),
        produceSplitLegRefits_(parameters.getParameter<bool>("produceSplitLegRefits")),
        directionalRefitUseExplicitBackwardTargetPropagation_(
            parameters.getParameter<bool>("directionalRefitUseExplicitBackwardTargetPropagation")),
        produceTargetConstrainedMomentum_(parameters.getParameter<bool>("produceTargetConstrainedMomentum")),
        targetUseInferredSide_(parameters.getParameter<bool>("targetUseInferredSide")),
        targetX_(parameters.getParameter<double>("targetX")),
        targetY_(parameters.getParameter<double>("targetY")),
        targetZ_(parameters.getParameter<double>("targetZ")),
        targetSigmaX_(parameters.getParameter<double>("targetSigmaX")),
        targetSigmaY_(parameters.getParameter<double>("targetSigmaY")),
        targetSigmaZ_(parameters.getParameter<double>("targetSigmaZ")),
        minSharedHitFraction_(parameters.getParameter<double>("minSharedHitFraction")),
        minSharedDetIds_(parameters.getParameter<unsigned int>("minSharedDetIds")),
        maxDuplicateAngle_(parameters.getParameter<double>("maxDuplicateAngle")),
        maxDuplicateLineDistance_(parameters.getParameter<double>("maxDuplicateLineDistance")),
        minTimingMeasurements_(parameters.getParameter<unsigned int>("minTimingMeasurements")),
        minTimingDeltaChi2_(parameters.getParameter<double>("minTimingDeltaChi2")),
        maxTargetLineDca_(parameters.getParameter<double>("maxTargetLineDca")),
        minAbsEta_(parameters.getParameter<double>("minAbsEta")),
        originTransverseResolution_(parameters.getParameter<double>("originTransverseResolution")),
        originZResolution_(parameters.getParameter<double>("originZResolution")),
        commonVertexLineResolution_(parameters.getParameter<double>("commonVertexLineResolution")),
        commonVertexBeamLineResolution_(parameters.getParameter<double>("commonVertexBeamLineResolution")),
        maxPairOriginNormalizedChi2_(parameters.getParameter<double>("maxPairOriginNormalizedChi2")),
        maxPairDca_(parameters.getParameter<double>("maxPairDca")),
        maxDimuonVertices_(parameters.getParameter<unsigned int>("maxDimuonVertices")),
        requireOppositeSign_(parameters.getParameter<bool>("requireOppositeSign")),
        maxGenDeltaR_(parameters.getParameter<double>("maxGenDeltaR")) {
    bool const anyGeometryFitMaterialEffects = directionalRefitUseGeometryMaterialEffects_ ||
                                               directionalRefitUseGeometryMaterialEffectsInFitter_ ||
                                               directionalRefitUseGeometryMaterialEffectsInSmoother_;
    bool const anyGeometryMaterialEffects =
        anyGeometryFitMaterialEffects || directionalRefitUseGeometryTargetMaterialEffects_;
    if (directionalRefitUseFirstPrinciplesMaterialEffects_ &&
        (useDetailedMaterialPropagation_ || directionalRefitUseDetailedMaterialEffects_ ||
         anyGeometryMaterialEffects))
      throw cms::Exception("Configuration")
          << "directionalRefitUseFirstPrinciplesMaterialEffects is mutually exclusive with the legacy "
             "Geant4e and geometry-material ablations";
    if (directionalRefitUseDetailedMaterialEffects_ && anyGeometryFitMaterialEffects)
      throw cms::Exception("Configuration")
          << "directionalRefitUseDetailedMaterialEffects and directionalRefitUseGeometryMaterialEffects "
             "are mutually exclusive";
    if (anyGeometryFitMaterialEffects && !directionalRefitUseMaterialEffects_)
      throw cms::Exception("Configuration")
          << "directionalRefitUseGeometryMaterialEffects requires directionalRefitUseMaterialEffects";
    if (directionalRefitUseGeometryTargetMaterialEffects_ && useDetailedMaterialPropagation_)
      throw cms::Exception("Configuration")
          << "directionalRefitUseGeometryTargetMaterialEffects and useDetailedMaterialPropagation are mutually "
             "exclusive target-leg material models";
    if (anyGeometryMaterialEffects &&
        (!(directionalRefitGeometryMaterialStepCm_ > 0.) || !std::isfinite(directionalRefitGeometryMaterialStepCm_)))
      throw cms::Exception("Configuration") << "directionalRefitGeometryMaterialStepCm must be finite and positive";
    if (directionalRefitUseFirstPrinciplesMaterialEffects_ &&
        (!(directionalRefitFirstPrinciplesStepCm_ > 0.) || !std::isfinite(directionalRefitFirstPrinciplesStepCm_)))
      throw cms::Exception("Configuration") << "directionalRefitFirstPrinciplesStepCm must be finite and positive";
    if (!(directionalRefitSecondSeedErrorRescale_ > 0.) ||
        !std::isfinite(directionalRefitSecondSeedErrorRescale_))
      throw cms::Exception("Configuration")
          << "directionalRefitSecondSeedErrorRescale must be finite and positive";
    if (produceTargetConstrainedMomentum_ &&
        (!(targetSigmaX_ > 0.) || !(targetSigmaY_ > 0.) || !(targetSigmaZ_ >= 0.) ||
         !std::isfinite(targetX_) || !std::isfinite(targetY_) || !std::isfinite(targetZ_) ||
         !std::isfinite(targetSigmaX_) || !std::isfinite(targetSigmaY_) || !std::isfinite(targetSigmaZ_)))
      throw cms::Exception("Configuration")
          << "target position must be finite, targetSigmaX/Y positive, and targetSigmaZ non-negative";
    for (auto const& tag : parameters.getParameter<std::vector<edm::InputTag>>("trackerSimHits"))
      trackerSimHitsTokens_.push_back(consumes<edm::PSimHitContainer>(tag));
    produces<nanoaod::FlatTable>();
    produces<nanoaod::FlatTable>("ShiftDimuonVertex");
  }

  void produce(edm::Event& event, edm::EventSetup const& setup) override {
    auto const dsa = event.getHandle(dsaToken_);
    auto const cosmic = event.getHandle(cosmicToken_);
    auto const traversing = event.getHandle(traversingToken_);
    auto const dsaGlobalLinks = event.getHandle(dsaGlobalLinksToken_);
    auto const cosmicGlobalLinks = event.getHandle(cosmicGlobalLinksToken_);
    auto const traversingGlobalLinks = event.getHandle(traversingGlobalLinksToken_);
    auto const genParticles = event.getHandle(genParticlesToken_);
    auto const simTracks = event.getHandle(simTracksToken_);
    auto const simVertices = event.getHandle(simVerticesToken_);
    auto const dtSimHits = event.getHandle(dtSimHitsToken_);
    auto const cscSimHits = event.getHandle(cscSimHitsToken_);
    auto const rpcSimHits = event.getHandle(rpcSimHitsToken_);
    auto const gemSimHits = event.getHandle(gemSimHitsToken_);
    std::vector<edm::Handle<edm::PSimHitContainer>> trackerSimHits;
    trackerSimHits.reserve(trackerSimHitsTokens_.size());
    for (auto const& token : trackerSimHitsTokens_)
      trackerSimHits.push_back(event.getHandle(token));
    auto const hcalSimHits = event.getHandle(hcalSimHitsToken_);
    auto const zdcSimHits = event.getHandle(zdcSimHitsToken_);
    auto const ecalBarrelRecHits = event.getHandle(ecalBarrelRecHitsToken_);
    auto const ecalEndcapRecHits = event.getHandle(ecalEndcapRecHitsToken_);
    auto const hbheRecHits = event.getHandle(hbheRecHitsToken_);
    auto const hfRecHits = event.getHandle(hfRecHitsToken_);
    auto const hoRecHits = event.getHandle(hoRecHitsToken_);
    auto const zdcRecHits = event.getHandle(zdcRecHitsToken_);
    auto const dtSegments = event.getHandle(dtSegmentsToken_);
    auto const pixelRecHits = event.getHandle(pixelRecHitsToken_);
    auto const stripMatchedRecHits = event.getHandle(stripMatchedRecHitsToken_);
    auto const stripRphiRecHits = event.getHandle(stripRphiRecHitsToken_);
    auto const stripRphiUnmatchedRecHits = event.getHandle(stripRphiUnmatchedRecHitsToken_);
    auto const stripStereoRecHits = event.getHandle(stripStereoRecHitsToken_);
    auto const stripStereoUnmatchedRecHits = event.getHandle(stripStereoUnmatchedRecHitsToken_);
    auto const& magneticField = setup.getData(magneticFieldToken_);
    bool const useGeometryMaterialInFitter = directionalRefitUseGeometryMaterialEffects_ ||
                                             directionalRefitUseGeometryMaterialEffectsInFitter_;
    bool const useGeometryMaterialInSmoother = directionalRefitUseGeometryMaterialEffects_ ||
                                               directionalRefitUseGeometryMaterialEffectsInSmoother_;
    // Keep the established R-Z material propagator as the default Kalman
    // transport.  A focused ablation can replace it with detailed Geant4e
    // transport without changing hit selection, ordering, or fit thresholds.
    SteppingHelixPropagator approximateMaterialPropagator(&magneticField, anyDirection);
    approximateMaterialPropagator.setMaterialMode(false);
    approximateMaterialPropagator.setUseMagVolumes(true);
    approximateMaterialPropagator.setUseMatVolumes(true);
    approximateMaterialPropagator.setEnergyLossScale(directionalRefitEnergyLossScale_);
    approximateMaterialPropagator.applyRadX0Correction(true);
    SteppingHelixPropagator explicitBackwardMaterialPropagator(&magneticField, oppositeToMomentum);
    explicitBackwardMaterialPropagator.setMaterialMode(false);
    explicitBackwardMaterialPropagator.setUseMagVolumes(true);
    explicitBackwardMaterialPropagator.setUseMatVolumes(true);
    explicitBackwardMaterialPropagator.setEnergyLossScale(directionalRefitEnergyLossScale_);
    explicitBackwardMaterialPropagator.applyRadX0Correction(true);
    auto firstPrinciplesProvider = std::make_shared<Geant4MaterialEffectsProvider>();
    SteppingHelixPropagator firstPrinciplesMaterialPropagator(&magneticField, anyDirection);
    firstPrinciplesMaterialPropagator.setMaterialMode(false);
    firstPrinciplesMaterialPropagator.setUseMagVolumes(true);
    firstPrinciplesMaterialPropagator.setUseMatVolumes(false);
    firstPrinciplesMaterialPropagator.setMaterialEffectsProvider(
        firstPrinciplesProvider, directionalRefitFirstPrinciplesStepCm_);
    firstPrinciplesMaterialPropagator.applyRadX0Correction(true);
    SteppingHelixPropagator firstPrinciplesBackwardMaterialPropagator(&magneticField, oppositeToMomentum);
    firstPrinciplesBackwardMaterialPropagator.setMaterialMode(false);
    firstPrinciplesBackwardMaterialPropagator.setUseMagVolumes(true);
    firstPrinciplesBackwardMaterialPropagator.setUseMatVolumes(false);
    firstPrinciplesBackwardMaterialPropagator.setMaterialEffectsProvider(
        firstPrinciplesProvider, directionalRefitFirstPrinciplesStepCm_);
    firstPrinciplesBackwardMaterialPropagator.applyRadX0Correction(true);
    std::unique_ptr<Geant4ePropagator> detailedMaterialPropagator;
    std::unique_ptr<Geant4ePropagator> detailedRefitMaterialPropagator;
    // Candidate preselection starts from a track endpoint whose direction is
    // canonicalised by propagateToTargetLine().  Preserve its established
    // anyDirection transport.  In contrast, the source-facing smoothed refit
    // state has a known physical orientation and its target leg must be
    // opposite-to-momentum.  These are different state contracts; sharing a
    // propagator silently applied the latter twice and removed all candidates.
    Propagator const* preRefitMaterialPropagator = &approximateMaterialPropagator;
    Propagator const* sourceFacingTargetMaterialPropagator = &explicitBackwardMaterialPropagator;
    if (useImprovedMomentumRefit_ && directionalRefitUseFirstPrinciplesMaterialEffects_) {
      preRefitMaterialPropagator = &firstPrinciplesMaterialPropagator;
      sourceFacingTargetMaterialPropagator = &firstPrinciplesBackwardMaterialPropagator;
    }
    if (useImprovedMomentumRefit_ &&
        (useDetailedMaterialPropagation_ || useGeometryMaterialInFitter || useGeometryMaterialInSmoother)) {
      // The stock Geant4e limits (10 mm steps and 200 cm total path) are too
      // coarse/short for the source-facing state to material-boundary leg.
      // This leg is geometrically behind the incoming muon's momentum, so use
      // an explicit direction instead of Geant4e's ambiguous anyDirection.
      detailedMaterialPropagator =
          std::make_unique<Geant4ePropagator>(&magneticField, "mu", oppositeToMomentum, 0.05, 2.0, 2500.0);
    }
    if (useImprovedMomentumRefit_ && useDetailedMaterialPropagation_) {
      preRefitMaterialPropagator = detailedMaterialPropagator.get();
      sourceFacingTargetMaterialPropagator = detailedMaterialPropagator.get();
    }
    SteppingHelixPropagator vacuumPropagator(&magneticField, anyDirection);
    vacuumPropagator.setMaterialMode(true);
    vacuumPropagator.setUseMagVolumes(true);
    SteppingHelixPropagator geometryCovariancePropagator(&magneticField, anyDirection);
    geometryCovariancePropagator.setMaterialMode(false);
    geometryCovariancePropagator.setUseMagVolumes(true);
    geometryCovariancePropagator.setUseMatVolumes(true);
    geometryCovariancePropagator.setEnergyLossScale(0.);
    geometryCovariancePropagator.applyRadX0Correction(true);
    GeometryMaterialPropagator geometryMaterialPropagator(
        geometryCovariancePropagator,
        approximateMaterialPropagator,
        directionalRefitGeometryMaterialStepCm_,
        directionalRefitLogGeometryMaterialComparison_);
    std::unique_ptr<SteppingHelixPropagator> geometryBackwardCovariancePropagator;
    std::unique_ptr<GeometryMaterialPropagator> geometryTargetMaterialPropagator;
    if (useImprovedMomentumRefit_ && directionalRefitUseGeometryTargetMaterialEffects_) {
      geometryBackwardCovariancePropagator =
          std::make_unique<SteppingHelixPropagator>(&magneticField, oppositeToMomentum);
      geometryBackwardCovariancePropagator->setMaterialMode(false);
      geometryBackwardCovariancePropagator->setUseMagVolumes(true);
      geometryBackwardCovariancePropagator->setUseMatVolumes(true);
      geometryBackwardCovariancePropagator->setEnergyLossScale(0.);
      geometryBackwardCovariancePropagator->applyRadX0Correction(true);
      geometryTargetMaterialPropagator = std::make_unique<GeometryMaterialPropagator>(
          *geometryBackwardCovariancePropagator,
          explicitBackwardMaterialPropagator,
          directionalRefitGeometryMaterialStepCm_,
          directionalRefitLogGeometryMaterialComparison_);
      sourceFacingTargetMaterialPropagator = geometryTargetMaterialPropagator.get();
    }
    Propagator const* directionalRefitFitterPropagator = directionalRefitUseMaterialEffects_
                                                             ? static_cast<Propagator const*>(&approximateMaterialPropagator)
                                                             : static_cast<Propagator const*>(&vacuumPropagator);
    Propagator const* directionalRefitSmootherPropagator = directionalRefitFitterPropagator;
    if (useImprovedMomentumRefit_ && directionalRefitUseMaterialEffects_ &&
        directionalRefitUseFirstPrinciplesMaterialEffects_) {
      directionalRefitFitterPropagator = &firstPrinciplesMaterialPropagator;
      directionalRefitSmootherPropagator = &firstPrinciplesMaterialPropagator;
    } else if (useImprovedMomentumRefit_ && directionalRefitUseMaterialEffects_ &&
        directionalRefitUseDetailedMaterialEffects_) {
      // Hits are ordered from the source-facing side toward CMS.  The fitter
      // therefore transports along the physical momentum.  KFTrajectoryFitter
      // clones this prototype with the seed direction, while
      // KFTrajectorySmoother creates its own alongMomentum and
      // oppositeToMomentum clones and selects the latter for the backward
      // smoothing pass.  Supplying an explicitly forward Geant4e prototype
      // thus gives the two Kalman roles physically correct material signs.
      detailedRefitMaterialPropagator =
          std::make_unique<Geant4ePropagator>(&magneticField, "mu", alongMomentum, 0.05, 2.0, 2500.0);
      directionalRefitFitterPropagator = detailedRefitMaterialPropagator.get();
      directionalRefitSmootherPropagator = detailedRefitMaterialPropagator.get();
    } else if (useImprovedMomentumRefit_ && directionalRefitUseMaterialEffects_) {
      if (useGeometryMaterialInFitter)
        directionalRefitFitterPropagator = &geometryMaterialPropagator;
      if (useGeometryMaterialInSmoother)
        directionalRefitSmootherPropagator = &geometryMaterialPropagator;
    }
    auto const& trackingGeometry = setup.getData(trackingGeometryToken_);
    auto const& caloGeometry = setup.getData(caloGeometryToken_);
    auto const& muonRecHitBuilder = setup.getData(muonRecHitBuilderToken_);
    auto const& trackerRecHitBuilder = setup.getData(trackerRecHitBuilderToken_);
    // KFTrajectoryFitter/Smoother use TkCloner as a generic rechit clone-or-
    // reuse helper. Muon transient rechits cannot improve with a tracker CPE,
    // so a default cloner simply returns their existing shared pointers and
    // never accesses tracker calibration objects.
    TkClonerImpl const hitCloner;

    std::vector<Candidate> candidates;
    auto append = [&candidates, &vacuumPropagator](auto const& handle, int source) {
      if (!handle.isValid())
        return;
      unsigned int index = 0;
      for (auto const& track : *handle) {
        Candidate candidate{};
        candidate.track = &track;
        candidate.source = source;
        candidate.sourceIndex = index++;
        candidate.hitFingerprints = hitFingerprints(track);
        candidate.targetLineState = propagateToTargetLine(track, vacuumPropagator);
        candidates.push_back(std::move(candidate));
      }
    };
    // Preserve the public source numbering; representative precedence is
    // applied explicitly below and need not follow this append order.
    append(dsa, 0);
    append(traversing, 1);
    append(cosmic, 2);

    for (auto& candidate : candidates) {
      auto const timing = timingDirection(*candidate.track,
                                          trackingGeometry,
                                          minTimingMeasurements_,
                                          minTimingDeltaChi2_,
                                          useExtendedTiming_);
      candidate.timingDirectionSign = timing.directionSign;
      candidate.timingMeasurements = timing.measurements;
      candidate.timingChi2 = timing.chi2;
      candidate.timingDeltaChi2 = timing.deltaChi2;
    }

    // Determine the common source side from all reconstructed hypotheses.
    // Its magnitude is never used as a selection, while the event-level sign
    // resolves the no-timing ambiguity consistently for both detector sides.
    int eventSourceSide = 1;
    bool eventSourceSideValid = false;
    double largestAbsOriginZ = 0.;
    unsigned int positiveOrigins = 0, negativeOrigins = 0;
    for (auto const& candidate : candidates) {
      if (!candidate.targetLineState.valid)
        continue;
      eventSourceSideValid = true;
      double const originZ = candidate.targetLineState.position.z();
      positiveOrigins += originZ > 0.;
      negativeOrigins += originZ < 0.;
      if (std::abs(originZ) > largestAbsOriginZ) {
        largestAbsOriginZ = std::abs(originZ);
        eventSourceSide = originZ < 0. ? -1 : 1;
      }
    }
    if (positiveOrigins != negativeOrigins)
      eventSourceSide = positiveOrigins > negativeOrigins ? 1 : -1;

    // Repeat the field-aware transport from the upstream fitted endpoint.
    // The preliminary inner-state result above is used only to infer the
    // common +/-z side and is never stored.
    for (auto& candidate : candidates) {
      int const directionSign = candidate.timingDirectionSign != 0
                                    ? candidate.timingDirectionSign
                                    : shiftDirectionSign(candidate.targetLineState.momentum, eventSourceSide);
      candidate.physicalDirectionSign = directionSign;
      auto const preRefitState = propagateToTargetLine(
          *candidate.track, vacuumPropagator, directionSign, preRefitMaterialPropagator, eventSourceSide);
      candidate.preRefitPt = preRefitState.valid ? preRefitState.momentum.perp() : 0.;
      candidate.preRefitPz = preRefitState.valid ? preRefitState.momentum.z() : 0.;
      candidate.targetLineState = preRefitState;
      auto const topology = muonHitTopology(*candidate.track);
      candidate.nDTRefitHits = topology.dt;
      candidate.nCSCRefitHits = topology.csc;
      candidate.nRPCRefitHits = topology.rpc;
      candidate.nGEMRefitHits = topology.gem;
      candidate.nPrecisionRefitStations = topology.precisionStations;

      // Start from the already reconstructed ShiftMuon and independently ask
      // which unassigned DT segments and tracker rechits its field-aware
      // trajectory actually crosses. This avoids any beam-spot/vertex or
      // collision-track prerequisite. The accepted measurements are then
      // included in the same two-pass directional Kalman refit below.
      std::vector<std::pair<TransientTrackingRecHit::RecHitPointer, bool>> additionalHits;
      std::set<uint32_t> existingDetIds;
      for (auto hit = candidate.track->recHitsBegin(); hit != candidate.track->recHitsEnd(); ++hit)
        if ((*hit)->isValid())
          existingDetIds.insert((*hit)->geographicalId().rawId());
      std::unordered_map<uint32_t, std::pair<double, TransientTrackingRecHit::RecHitPointer>> bestAdditionalHits;
      Chi2MeasurementEstimator const additionalHitEstimator(additionalHitMaxChi2_);
      bool const useOuterAssociationState = eventSourceSide * candidate.track->outerPosition().z() >
                                            eventSourceSide * candidate.track->innerPosition().z();
      auto associationOriginal =
          useOuterAssociationState ? trajectoryStateTransform::outerFreeState(*candidate.track, &magneticField)
                                   : trajectoryStateTransform::innerFreeState(*candidate.track, &magneticField);
      std::unique_ptr<FreeTrajectoryState> associationSeed;
      if (associationOriginal.hasError()) {
        double const associationSign = directionSign == 0 ? 1. : directionSign;
        GlobalTrajectoryParameters const parameters(associationOriginal.position(),
                                                    associationSign * associationOriginal.momentum(),
                                                    associationSign * associationOriginal.charge(),
                                                    &magneticField);
        associationSeed =
            std::make_unique<FreeTrajectoryState>(parameters, associationOriginal.curvilinearError());
      }
      auto considerHit = [&](TransientTrackingRecHit::RecHitPointer const& transientHit,
                             double maxDistance,
                             unsigned int& propagatedCount,
                             unsigned int& compatibleCount,
                             double& minimumResidual,
                             double& minimumEstimatorChi2) {
        if (!transientHit || !transientHit->isValid() || !transientHit->det() ||
            existingDetIds.count(transientHit->geographicalId().rawId()) || !associationSeed)
          return;
        auto const propagated = vacuumPropagator.propagate(*associationSeed, *transientHit->surface());
        if (!propagated.isValid())
          return;
        ++propagatedCount;
        double const residual = (propagated.globalPosition() - transientHit->globalPosition()).mag();
        if (std::isfinite(residual) && (minimumResidual < 0. || residual < minimumResidual))
          minimumResidual = residual;
        auto const estimate = additionalHitEstimator.estimate(propagated, *transientHit);
        if (std::isfinite(estimate.second) &&
            (minimumEstimatorChi2 < 0. || estimate.second < minimumEstimatorChi2))
          minimumEstimatorChi2 = estimate.second;
        if (!std::isfinite(residual) || (maxDistance >= 0. && residual > maxDistance) || !estimate.first)
          return;
        ++compatibleCount;
        auto& best = bestAdditionalHits[transientHit->geographicalId().rawId()];
        // Do not reserve a detector after the first passing object.  DT often
        // has several reconstructed segment hypotheses in one chamber; keep
        // scanning and retain the covariance-compatible segment with the
        // smallest estimator chi2.
        if (!best.second || estimate.second < best.first)
          best = {estimate.second, transientHit};
      };
      if (augmentDTHits_ && dtSegments.isValid())
        for (auto segment = dtSegments->begin(); segment != dtSegments->end(); ++segment)
          considerHit(muonRecHitBuilder.build(&*segment),
                      additionalDTMaxDistance_,
                      candidate.nPropagatedDTSegments,
                      candidate.nCompatibleDTSegments,
                      candidate.minDTResidual,
                      candidate.minDTEstimatorChi2);
      auto considerTrackerCollection = [&](auto const& handle) {
        if (!handle.isValid())
          return;
        for (auto const& detSet : *handle)
          for (auto const& hit : detSet)
            considerHit(trackerRecHitBuilder.build(&hit),
                        additionalTrackerMaxDistance_,
                        candidate.nPropagatedStripHits,
                        candidate.nCompatibleStripHits,
                        candidate.minTrackerResidual,
                        candidate.minTrackerEstimatorChi2);
      };
      if (augmentTrackerHits_) {
        if (pixelRecHits.isValid())
          for (auto const& detSet : *pixelRecHits)
            for (auto const& hit : detSet)
              considerHit(trackerRecHitBuilder.build(&hit),
                          additionalTrackerMaxDistance_,
                          candidate.nPropagatedPixelHits,
                          candidate.nCompatiblePixelHits,
                          candidate.minTrackerResidual,
                          candidate.minTrackerEstimatorChi2);
        considerTrackerCollection(stripMatchedRecHits);
        considerTrackerCollection(stripRphiRecHits);
        considerTrackerCollection(stripRphiUnmatchedRecHits);
        considerTrackerCollection(stripStereoRecHits);
        considerTrackerCollection(stripStereoUnmatchedRecHits);
      }
      for (auto const& [rawId, best] : bestAdditionalHits) {
        bool const precision = DetId(rawId).det() == DetId::Muon;
        additionalHits.emplace_back(best.second, precision);
        candidate.nAddedDTRefitHits += precision;
        candidate.nAddedTrackerRefitHits += !precision;
        if (precision)
          candidate.addedDTChamberIds.insert(DTChamberId(rawId).rawId());
      }

      // Calorimeters and ZDC are independent timing/direction tags, not
      // precision trajectory measurements.  Associate their reconstructed
      // cells to the already fitted detector chord and retain enough
      // information to validate MIP efficiency and time ordering before any
      // covariance-bearing constraint is considered.
      struct TimedCaloPoint {
        GlobalPoint position;
        double time = 0.;
      };
      std::vector<TimedCaloPoint> timedCaloPoints;
      auto associateCalo = [&](auto const& handle,
                               int surfaceMode,
                               bool useForDirectionTiming,
                               double maxDistance,
                               unsigned int& matched,
                               double& energy,
                               double& meanTime,
                               double& minimumDistance) {
        if (!handle.isValid() || !associationSeed)
          return;
        // Rechits in a calorimeter section share only a small number of
        // cylindrical radii or endcap z planes. Cache the propagated helix
        // intersection per centimetre instead of propagating once per cell.
        std::unordered_map<int, std::pair<bool, GlobalPoint>> intersections;
        double const previousTimeWeight = energy;
        double timeWeight = 0.;
        double weightedTime = 0.;
        double weightedX = 0., weightedY = 0., weightedZ = 0.;
        for (auto const& hit : *handle) {
          if (hit.energy() < caloMatchMinEnergy_)
            continue;
          auto const cell = caloGeometry.getGeometry(hit.id());
          if (!cell)
            continue;
          auto const position = cell->getPosition();
          bool const usePlane = surfaceMode > 0 || (surfaceMode < 0 && std::abs(position.z()) > position.perp());
          double const coordinate = usePlane ? position.z() : position.perp();
          int const key = (usePlane ? 1000000 : 0) + static_cast<int>(std::lround(coordinate));
          auto intersection = intersections.find(key);
          if (intersection == intersections.end()) {
            TrajectoryStateOnSurface propagated;
            if (usePlane) {
              auto const surface = Plane::build(GlobalPoint(0., 0., coordinate), Surface::RotationType());
              propagated = vacuumPropagator.propagate(*associationSeed, *surface);
            } else {
              auto const surface = Cylinder::build(GlobalPoint(0., 0., 0.), Surface::RotationType(), coordinate);
              propagated = vacuumPropagator.propagate(*associationSeed, *surface);
            }
            intersection = intersections
                               .emplace(key,
                                        std::make_pair(propagated.isValid(),
                                                       propagated.isValid() ? propagated.globalPosition()
                                                                            : GlobalPoint()))
                               .first;
          }
          if (!intersection->second.first)
            continue;
          double const distance = (intersection->second.second - position).mag();
          if (std::isfinite(distance) && (minimumDistance < 0. || distance < minimumDistance))
            minimumDistance = distance;
          if (!std::isfinite(distance) || distance > maxDistance)
            continue;
          ++matched;
          energy += hit.energy();
          double const weight = std::max(std::abs(static_cast<double>(hit.energy())), 1.e-6);
          timeWeight += weight;
          weightedTime += weight * hit.time();
          weightedX += weight * position.x();
          weightedY += weight * position.y();
          weightedZ += weight * position.z();
        }
        if (timeWeight > 0.) {
          meanTime = (previousTimeWeight * meanTime + weightedTime) / (previousTimeWeight + timeWeight);
          if (useForDirectionTiming)
            timedCaloPoints.push_back(
                {GlobalPoint(weightedX / timeWeight, weightedY / timeWeight, weightedZ / timeWeight), meanTime});
        }
      };
      associateCalo(ecalBarrelRecHits,
                    0,
                    true,
                    ecalMatchMaxDistance_,
                    candidate.nMatchedEcalRecHits,
                    candidate.matchedEcalEnergy,
                    candidate.matchedEcalTime,
                    candidate.minEcalLineDistance);
      associateCalo(ecalEndcapRecHits,
                    1,
                    true,
                    ecalMatchMaxDistance_,
                    candidate.nMatchedEcalRecHits,
                    candidate.matchedEcalEnergy,
                    candidate.matchedEcalTime,
                    candidate.minEcalLineDistance);
      associateCalo(hbheRecHits,
                    -1,
                    true,
                    hcalMatchMaxDistance_,
                    candidate.nMatchedHBHERecHits,
                    candidate.matchedHBHEEnergy,
                    candidate.matchedHBHETime,
                    candidate.minHBHELineDistance);
      associateCalo(hfRecHits,
                    1,
                    true,
                    hcalMatchMaxDistance_,
                    candidate.nMatchedHFRecHits,
                    candidate.matchedHFEnergy,
                    candidate.matchedHFTime,
                    candidate.minHFLineDistance);
      associateCalo(hoRecHits,
                    0,
                    true,
                    hcalMatchMaxDistance_,
                    candidate.nMatchedHORecHits,
                    candidate.matchedHOEnergy,
                    candidate.matchedHOTime,
                    candidate.minHOLineDistance);
      associateCalo(zdcRecHits,
                    1,
                    false,
                    zdcMatchMaxDistance_,
                    candidate.nMatchedZDCRecHits,
                    candidate.matchedZDCEnergy,
                    candidate.matchedZDCTime,
                    candidate.minZDCLineDistance);

      if (timedCaloPoints.size() >= 2) {
        constexpr double inverseSpeedOfLightNsPerCm = 1. / 29.9792458;
        constexpr double timingSigmaNs = 5.;
        GlobalPoint const reference(candidate.track->innerPosition().x(),
                                    candidate.track->innerPosition().y(),
                                    candidate.track->innerPosition().z());
        GlobalVector const axis(candidate.track->innerMomentum().x(),
                               candidate.track->innerMomentum().y(),
                               candidate.track->innerMomentum().z());
        if (axis.mag2() > 0.) {
          auto directionChi2 = [&](double sign) {
            double intercept = 0.;
            for (auto const& point : timedCaloPoints)
              intercept += point.time -
                           sign * (point.position - reference).dot(axis.unit()) * inverseSpeedOfLightNsPerCm;
            intercept /= timedCaloPoints.size();
            double chi2 = 0.;
            for (auto const& point : timedCaloPoints) {
              double const residual = point.time - intercept -
                                      sign * (point.position - reference).dot(axis.unit()) *
                                          inverseSpeedOfLightNsPerCm;
              chi2 += residual * residual / (timingSigmaNs * timingSigmaNs);
            }
            return chi2;
          };
          double const alongChi2 = directionChi2(1.);
          double const oppositeChi2 = directionChi2(-1.);
          candidate.nCaloTimingMeasurements = timedCaloPoints.size();
          candidate.caloTimingDeltaChi2 = std::abs(alongChi2 - oppositeChi2);
          if (candidate.caloTimingDeltaChi2 >= minTimingDeltaChi2_)
            candidate.caloTimingDirectionSign = alongChi2 < oppositeChi2 ? 1 : -1;
        }
      }

      // Refit the same reconstructed hits in their measured time-of-flight
      // direction.  This changes the direction in which the Kalman filter
      // applies material updates; flipping a completed fit cannot do that.
      candidate.directionalRefitAttempted = true;
      auto runDirectionalRefit = [&](bool precisionHitsOnly, int hitSideSelection = 0) {
        return directionalRefit(*candidate.track,
                                directionSign,
                                eventSourceSide,
                                magneticField,
                                muonRecHitBuilder,
                                hitCloner,
                                vacuumPropagator,
                                *directionalRefitFitterPropagator,
                                *directionalRefitSmootherPropagator,
                                *sourceFacingTargetMaterialPropagator,
                                directionalRefitSeedCurvatureErrorRescale_,
                                directionalRefitSecondSeedErrorRescale_,
                                directionalRefitSeedMomentumScale_,
                                directionalRefitUseFullSeedErrorRescale_,
                                useImprovedMomentumRefit_ ? directionalRefitErrorRescale_ : 100.,
                                useImprovedMomentumRefit_ ? directionalRefitInitialMaxHitChi2_ : 100000.,
                                useImprovedMomentumRefit_ ? directionalRefitMaxHitChi2_ : 100000.,
                                useImprovedMomentumRefit_ ? directionalRefitMaxRelativeQoverPChange_ : 0.5,
                                directionalRefitUseSecondIteration_,
                                precisionHitsOnly,
                                useImprovedMomentumRefit_ && usePropagatedPathOrdering_,
                                hitSideSelection,
                                &additionalHits);
      };
      auto const allHitsRefit = runDirectionalRefit(false);
      auto const precisionRefit = useImprovedMomentumRefit_ ? runDirectionalRefit(true) : DirectionalRefitResult{};
      auto const sourceLegPrecisionRefit = produceSplitLegRefits_ && useImprovedMomentumRefit_
                                               ? runDirectionalRefit(true, 1)
                                               : DirectionalRefitResult{};
      auto const oppositeLegPrecisionRefit = produceSplitLegRefits_ && useImprovedMomentumRefit_
                                                 ? runDirectionalRefit(true, -1)
                                                 : DirectionalRefitResult{};
      double precisionRelativeToAllQoverP = -1.;
      bool precisionAgreesWithAllHits = true;
      if (precisionRefit.valid && allHitsRefit.valid) {
        double const denominator = std::max(std::abs(allHitsRefit.selected.signedInverseMomentum), 1.e-12);
        precisionRelativeToAllQoverP =
            std::abs(precisionRefit.selected.signedInverseMomentum - allHitsRefit.selected.signedInverseMomentum) /
            denominator;
        precisionAgreesWithAllHits = std::isfinite(precisionRelativeToAllQoverP) &&
                                     precisionRelativeToAllQoverP <=
                                         directionalRefitMaxPrecisionRelativeQoverPChange_;
      }
      double precisionTargetDca = -1.;
      bool precisionHasShiftTopology = false;
      if (precisionRefit.valid) {
        auto const& precisionTarget = precisionRefit.selected.materialTargetState;
        precisionTargetDca = precisionTarget.position.perp();
        double const precisionPt = precisionTarget.momentum.perp();
        double const precisionAbsEta =
            precisionPt > 0. ? std::abs(std::asinh(precisionTarget.momentum.z() / precisionPt))
                             : std::numeric_limits<double>::infinity();
        precisionHasShiftTopology = precisionTargetDca <= maxTargetLineDca_ && precisionAbsEta >= minAbsEta_;
      }
      // Enabling the tracker study must be a strict no-op when no tracker
      // measurement was accepted.  Only a genuinely augmented fit switches
      // away from the established precision-muon result.
      bool const usePrecisionRefit = useImprovedMomentumRefit_ && candidate.nAddedTrackerRefitHits == 0 &&
                                     precisionRefit.valid &&
                                     precisionRefit.inputStations >= directionalRefitMinPrecisionStations_ &&
                                     precisionHasShiftTopology && precisionAgreesWithAllHits;
      auto const& refit = usePrecisionRefit ? precisionRefit : allHitsRefit;

      // Produce a second, explicitly prompt-target hypothesis without ever
      // replacing the unconstrained result above.  A smoothed state is the
      // detector-hit likelihood expressed at one surface; transport it to the
      // target and make one Kalman update with the production measurement.
      // This is the linear-Gaussian equivalent of adding the constraint to the
      // full fit, without asking KFTrajectoryFitter to treat a detId=0 prior as
      // an ordinary detector hit.
      if (produceTargetConstrainedMomentum_ && refit.valid) {
        double const configuredTargetZ = targetUseInferredSide_ ? eventSourceSide * std::abs(targetZ_) : targetZ_;
        TargetConstraint const constraint{
            GlobalPoint(targetX_, targetY_, configuredTargetZ), targetSigmaX_, targetSigmaY_, targetSigmaZ_};
        auto const constrained = applyTargetConstraint(refit.selected.upstreamState,
                                                       eventSourceSide,
                                                       vacuumPropagator,
                                                       *sourceFacingTargetMaterialPropagator,
                                                       constraint);
        candidate.constrainedStatus = constrained.status;
        if (constrained.valid) {
          candidate.constrainedValid = true;
          candidate.constrainedHits = refit.selected.hits;
          candidate.constrainedChi2 = refit.selected.chi2 + constrained.chi2;
          candidate.constrainedNdof = refit.selected.ndof + 2.;
          candidate.constrainedTargetChi2 = constrained.chi2;
          candidate.constrainedState = constrained.state;
        }
      } else if (produceTargetConstrainedMomentum_) {
        candidate.constrainedStatus = -10;
      }

      candidate.precisionRefitLeverArm = precisionRefit.inputLeverArm;
      candidate.directionalRefitUsedPrecisionHits = usePrecisionRefit;
      candidate.directionalRefitAllHitsValid = allHitsRefit.valid;
      candidate.directionalRefitAllHits = allHitsRefit.selected.hits;
      candidate.directionalRefitAllHitsChi2 = allHitsRefit.selected.chi2;
      candidate.directionalRefitAllHitsNdof = allHitsRefit.selected.ndof;
      candidate.directionalRefitAllHitsTargetPt =
          allHitsRefit.valid ? allHitsRefit.selected.materialTargetState.momentum.perp() : 0.;
      candidate.directionalRefitAllHitsTargetPz =
          allHitsRefit.valid ? allHitsRefit.selected.materialTargetState.momentum.z() : 0.;
      candidate.directionalRefitAllHitsSelectedIteration = allHitsRefit.selectedIteration;
      candidate.directionalRefitAllHitsInput = allHitsRefit.inputHits;
      candidate.directionalRefitAllHitsOrderingFallback = allHitsRefit.pathOrderingFallbackHits;
      candidate.directionalRefitAllHitsOrderingSpan = allHitsRefit.pathOrderingSpan;
      candidate.directionalRefitAllHitsFirstRejected =
          allHitsRefit.inputHits > allHitsRefit.first.hits ? allHitsRefit.inputHits - allHitsRefit.first.hits : 0;
      candidate.directionalRefitAllHitsSecondRejected =
          allHitsRefit.second.valid && allHitsRefit.inputHits > allHitsRefit.second.hits
              ? allHitsRefit.inputHits - allHitsRefit.second.hits
              : 0;
      candidate.directionalRefitPrecisionValid = precisionRefit.valid;
      candidate.directionalRefitPrecisionHits = precisionRefit.selected.hits;
      candidate.directionalRefitPrecisionChi2 = precisionRefit.selected.chi2;
      candidate.directionalRefitPrecisionNdof = precisionRefit.selected.ndof;
      candidate.directionalRefitPrecisionUpstreamPt = precisionRefit.selected.upstreamPt;
      candidate.directionalRefitPrecisionQoverP = precisionRefit.selected.signedInverseMomentum;
      candidate.directionalRefitPrecisionTargetPt =
          precisionRefit.valid ? precisionRefit.selected.materialTargetState.momentum.perp() : 0.;
      candidate.directionalRefitPrecisionTargetPz =
          precisionRefit.valid ? precisionRefit.selected.materialTargetState.momentum.z() : 0.;
      candidate.directionalRefitPrecisionSecondValid = precisionRefit.second.valid;
      candidate.directionalRefitPrecisionSecondConverged = precisionRefit.secondConverged;
      candidate.directionalRefitPrecisionRelativeQoverPChange = precisionRefit.relativeQoverPChange;
      candidate.directionalRefitPrecisionSelectedIteration = precisionRefit.selectedIteration;
      candidate.directionalRefitPrecisionFirstTargetPt =
          precisionRefit.first.valid ? precisionRefit.first.materialTargetState.momentum.perp() : 0.;
      candidate.directionalRefitPrecisionFirstTargetPz =
          precisionRefit.first.valid ? precisionRefit.first.materialTargetState.momentum.z() : 0.;
      candidate.directionalRefitPrecisionSecondTargetPt =
          precisionRefit.second.valid ? precisionRefit.second.materialTargetState.momentum.perp() : 0.;
      candidate.directionalRefitPrecisionSecondTargetPz =
          precisionRefit.second.valid ? precisionRefit.second.materialTargetState.momentum.z() : 0.;
      candidate.directionalRefitPrecisionRelativeToAllQoverP = precisionRelativeToAllQoverP;
      candidate.directionalRefitPrecisionTargetDca = precisionTargetDca;
      candidate.directionalRefitSourceLegPrecisionValid = sourceLegPrecisionRefit.valid;
      candidate.directionalRefitSourceLegPrecisionHits = sourceLegPrecisionRefit.selected.hits;
      candidate.directionalRefitSourceLegPrecisionQoverP = sourceLegPrecisionRefit.selected.signedInverseMomentum;
      candidate.directionalRefitSourceLegPrecisionTargetPt =
          sourceLegPrecisionRefit.valid ? sourceLegPrecisionRefit.selected.materialTargetState.momentum.perp() : 0.;
      candidate.directionalRefitSourceLegPrecisionTargetPz =
          sourceLegPrecisionRefit.valid ? sourceLegPrecisionRefit.selected.materialTargetState.momentum.z() : 0.;
      candidate.directionalRefitOppositeLegPrecisionValid = oppositeLegPrecisionRefit.valid;
      candidate.directionalRefitOppositeLegPrecisionHits = oppositeLegPrecisionRefit.selected.hits;
      candidate.directionalRefitOppositeLegPrecisionQoverP = oppositeLegPrecisionRefit.selected.signedInverseMomentum;
      candidate.directionalRefitOppositeLegPrecisionTargetPt =
          oppositeLegPrecisionRefit.valid ? oppositeLegPrecisionRefit.selected.materialTargetState.momentum.perp() : 0.;
      candidate.directionalRefitOppositeLegPrecisionTargetPz =
          oppositeLegPrecisionRefit.valid ? oppositeLegPrecisionRefit.selected.materialTargetState.momentum.z() : 0.;
      candidate.directionalRefitPrecisionInput = precisionRefit.inputHits;
      candidate.directionalRefitPrecisionOrderingFallback = precisionRefit.pathOrderingFallbackHits;
      candidate.directionalRefitPrecisionOrderingSpan = precisionRefit.pathOrderingSpan;
      candidate.directionalRefitPrecisionFirstRejected =
          precisionRefit.inputHits > precisionRefit.first.hits ? precisionRefit.inputHits - precisionRefit.first.hits : 0;
      candidate.directionalRefitPrecisionSecondRejected =
          precisionRefit.second.valid && precisionRefit.inputHits > precisionRefit.second.hits
              ? precisionRefit.inputHits - precisionRefit.second.hits
              : 0;

      // Keep the established first/second diagnostic branches tied to the
      // all-hit control fit so this production can compare it directly with
      // the new precision-only variant.
      candidate.directionalRefitValid = refit.valid;
      candidate.directionalRefitFirstValid = allHitsRefit.first.valid;
      candidate.directionalRefitFirstHits = allHitsRefit.first.hits;
      candidate.directionalRefitFirstChi2 = allHitsRefit.first.chi2;
      candidate.directionalRefitFirstNdof = allHitsRefit.first.ndof;
      candidate.directionalRefitFirstUpstreamPt = allHitsRefit.first.upstreamPt;
      candidate.directionalRefitFirstQoverP = allHitsRefit.first.signedInverseMomentum;
      candidate.directionalRefitFirstTargetPt =
          allHitsRefit.first.materialTargetState.valid ? allHitsRefit.first.materialTargetState.momentum.perp() : 0.;
      candidate.directionalRefitFirstTargetPz =
          allHitsRefit.first.materialTargetState.valid ? allHitsRefit.first.materialTargetState.momentum.z() : 0.;
      candidate.directionalRefitSecondValid = allHitsRefit.second.valid;
      candidate.directionalRefitSecondHits = allHitsRefit.second.hits;
      candidate.directionalRefitSecondChi2 = allHitsRefit.second.chi2;
      candidate.directionalRefitSecondNdof = allHitsRefit.second.ndof;
      candidate.directionalRefitSecondUpstreamPt = allHitsRefit.second.upstreamPt;
      candidate.directionalRefitSecondQoverP = allHitsRefit.second.signedInverseMomentum;
      candidate.directionalRefitSecondTargetPt =
          allHitsRefit.second.materialTargetState.valid ? allHitsRefit.second.materialTargetState.momentum.perp() : 0.;
      candidate.directionalRefitSecondTargetPz =
          allHitsRefit.second.materialTargetState.valid ? allHitsRefit.second.materialTargetState.momentum.z() : 0.;
      candidate.directionalRefitSecondConverged = allHitsRefit.secondConverged;
      candidate.directionalRefitRelativeQoverPChange = allHitsRefit.relativeQoverPChange;
      candidate.directionalRefitSelectedIteration = allHitsRefit.selectedIteration;
      candidate.directionalRefitHits = refit.selected.hits;
      candidate.directionalRefitChi2 = refit.selected.chi2;
      candidate.directionalRefitNdof = refit.selected.ndof;
      candidate.directionalRefitUpstreamPt = refit.selected.upstreamPt;
      candidate.directionalRefitUpstreamP = refit.selected.upstreamP;
      candidate.directionalRefitDownstreamP = refit.selected.downstreamP;
      candidate.directionalRefitSourceFacingSurfaceDistance = refit.selected.sourceFacingSurfaceDistance;
      candidate.directionalRefitFilterDirectionCosine = refit.selected.filterDirectionCosine;
      candidate.directionalRefitTargetDirectionCosine = refit.selected.targetDirectionCosine;
      if (refit.selected.upstreamP > 0.)
        candidate.directionalRefitFractionalLossAcrossHits =
            (refit.selected.upstreamP - refit.selected.downstreamP) / refit.selected.upstreamP;
      if (refit.valid) {
        candidate.targetLineState = refit.selected.materialTargetState;
        if (refit.selected.vacuumTargetState.valid) {
          candidate.directionalRefitVacuumTargetPt = refit.selected.vacuumTargetState.momentum.perp();
          candidate.directionalRefitVacuumTargetPz = refit.selected.vacuumTargetState.momentum.z();
          candidate.directionalRefitMaterialDeltaPt =
              refit.selected.materialTargetState.momentum.perp() -
              refit.selected.vacuumTargetState.momentum.perp();
        }
        auto const& boundaryState = refit.selected.materialTargetState;
        candidate.directionalRefitMaterialBoundaryValid = boundaryState.materialBoundaryValid;
        if (boundaryState.materialBoundaryValid) {
          candidate.directionalRefitMaterialBoundaryPt = boundaryState.materialBoundaryMomentum.perp();
          candidate.directionalRefitMaterialBoundaryPz = boundaryState.materialBoundaryMomentum.z();
          candidate.directionalRefitMaterialPath = boundaryState.materialPath;
        }
      }
    }

    // Do not select on reconstructed z.  Reject only invalid/empty fits and
    // trajectories which do not point back to the unbounded target line in
    // the transverse plane, or which lack the forward topology of particles
    // arriving along that line.  These remove mixed-hit cosmic candidates
    // without preferring either source side or any longitudinal location.
    candidates.erase(std::remove_if(candidates.begin(),
                                    candidates.end(),
                                    [this](Candidate const& candidate) {
                                      double const pt = candidate.targetLineState.momentum.perp();
                                      double const absEta =
                                          pt > 0. ? std::abs(std::asinh(candidate.targetLineState.momentum.z() / pt))
                                                  : std::numeric_limits<double>::infinity();
                                      return !candidate.targetLineState.valid ||
                                             candidate.track->numberOfValidHits() == 0 ||
                                             candidate.targetLineState.position.perp() > maxTargetLineDca_ ||
                                             absEta < minAbsEta_;
                                    }),
                     candidates.end());

    std::stable_sort(candidates.begin(), candidates.end(), [](Candidate const& first, Candidate const& second) {
      if (first.source != second.source)
        return first.source < second.source;
      if (first.track->numberOfValidHits() != second.track->numberOfValidHits())
        return first.track->numberOfValidHits() > second.track->numberOfValidHits();
      return first.track->normalizedChi2() < second.track->normalizedChi2();
    });

    auto duplicates = [this](Candidate const& first, Candidate const& second) {
      // A chamber/detId is much coarser than a hit.  Two nearby J/psi muons,
      // or the two legs of a through-going track, often cross the same
      // chambers at different local positions and must not be merged for
      // that reason alone.
      unsigned int sharedInputs = 0;
      for (auto const& firstHit : first.hitFingerprints) {
        bool sharesThisHit = false;
        for (auto const& secondHit : second.hitFingerprints)
          if (firstHit.detId == secondHit.detId &&
              std::hypot(firstHit.localX - secondHit.localX, firstHit.localY - secondHit.localY) < 0.1) {
            sharesThisHit = true;
            break;
          }
        sharedInputs += sharesThisHit;
      }
      auto const smallerRecHitSet = std::min(first.hitFingerprints.size(), second.hitFingerprints.size());
      bool const exactSharedInputs = smallerRecHitSet > 0 && sharedInputs >= minSharedDetIds_ &&
                                     static_cast<double>(sharedInputs) / smallerRecHitSet >= minSharedHitFraction_;
      bool const sameLine = directionAngle(*first.track, *second.track) < maxDuplicateAngle_ &&
                            symmetricLineDistance(*first.track, *second.track) < maxDuplicateLineDistance_;
      return exactSharedInputs || sameLine;
    };

    // Collapse complete connected components rather than applying an
    // order-dependent greedy veto.  This handles A<->B<->C duplicate chains
    // even when the two endpoint fits narrowly fail the direct comparison.
    std::vector<unsigned int> parent(candidates.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto findRoot = [&parent](unsigned int index) {
      while (parent[index] != index) {
        parent[index] = parent[parent[index]];
        index = parent[index];
      }
      return index;
    };
    for (unsigned int first = 0; first < candidates.size(); ++first)
      for (unsigned int second = first + 1; second < candidates.size(); ++second)
        if (duplicates(candidates[first], candidates[second])) {
          unsigned int const firstRoot = findRoot(first);
          unsigned int const secondRoot = findRoot(second);
          if (firstRoot != secondRoot)
            parent[secondRoot] = firstRoot;
        }

    std::vector<std::vector<unsigned int>> groups(candidates.size());
    for (unsigned int index = 0; index < candidates.size(); ++index)
      groups[findRoot(index)].push_back(index);
    auto betterRepresentative = [&candidates, eventSourceSide](unsigned int first, unsigned int second) {
      auto const& a = candidates[first];
      auto const& b = candidates[second];
      auto geometryPriority = [](int source) {
        // A traversing fit constrains one curvature with measurements on both
        // sides of CMS.  That information cannot be replaced by timing or by
        // a successful refit of only one leg.  Prefer the longest physical
        // hypothesis first, then use timing and fit quality to choose among
        // hypotheses with the same topology.
        // The non-traversing cosmic fit is retained only as a last-resort
        // seed: in this forward topology it has no longer curvature lever arm
        // than DSA and is less constrained by the SHIFT-specific seeding.
        return source == 1 ? 0 : (source == 0 ? 1 : 2);
      };
      if (geometryPriority(a.source) != geometryPriority(b.source))
        return geometryPriority(a.source) < geometryPriority(b.source);
      bool const aHasTiming = a.timingDirectionSign != 0;
      bool const bHasTiming = b.timingDirectionSign != 0;
      if (aHasTiming != bHasTiming)
        return aHasTiming;
      // Both seed orientations were fitted.  Prefer the fit whose native
      // propagation direction agrees with time of flight; merely flipping a
      // completed wrong-way fit would retain its biased material corrections
      // and is the source of the near-zero-pT solution.
      if (aHasTiming && (a.timingDirectionSign > 0) != (b.timingDirectionSign > 0))
        return a.timingDirectionSign > 0;
      if (aHasTiming && bHasTiming && a.timingChi2 != b.timingChi2)
        return a.timingChi2 < b.timingChi2;
      // A successfully refitted duplicate is preferable to a candidate that
      // would fall back to the original biased endpoint state. This remains a
      // reconstruction-only decision and is evaluated before source labels.
      if (a.directionalRefitValid != b.directionalRefitValid)
        return a.directionalRefitValid;
      // For a far source, a standalone leg in the endcap on the same side as
      // its propagated target-line PCA is the upstream leg.  It has crossed
      // less field and material than a downstream leg and gives a markedly
      // less ambiguous backward transport.  This uses reconstructed
      // geometry only and works symmetrically for +/-z sources.
      bool const aUpstream = eventSourceSide * a.track->innerPosition().z() > 0. &&
                             eventSourceSide * a.track->outerPosition().z() > 0.;
      bool const bUpstream = eventSourceSide * b.track->innerPosition().z() > 0. &&
                             eventSourceSide * b.track->outerPosition().z() > 0.;
      if (aUpstream != bUpstream)
        return aUpstream;
      if (a.track->hitPattern().muonStationsWithValidHits() != b.track->hitPattern().muonStationsWithValidHits())
        return a.track->hitPattern().muonStationsWithValidHits() > b.track->hitPattern().muonStationsWithValidHits();
      if (a.track->numberOfValidHits() != b.track->numberOfValidHits())
        return a.track->numberOfValidHits() > b.track->numberOfValidHits();
      return a.track->normalizedChi2() < b.track->normalizedChi2();
    };
    std::vector<Candidate const*> selected;
    std::vector<unsigned int> duplicateGroupSize, duplicateSourceMask, duplicateDSACount,
        duplicateTraversingCount, duplicateCosmicCount;
    for (auto const& group : groups) {
      if (group.empty())
        continue;
      unsigned int representative = group.front();
      for (auto const member : group)
        if (betterRepresentative(member, representative))
          representative = member;
      selected.push_back(&candidates[representative]);
      duplicateGroupSize.push_back(group.size());
      unsigned int sourceMask = 0, dsaCount = 0, traversingCount = 0, cosmicCount = 0;
      for (auto const member : group) {
        int const source = candidates[member].source;
        if (source >= 0 && source < 3)
          sourceMask |= 1u << source;
        dsaCount += source == 0;
        traversingCount += source == 1;
        cosmicCount += source == 2;
      }
      duplicateSourceMask.push_back(sourceMask);
      duplicateDSACount.push_back(dsaCount);
      duplicateTraversingCount.push_back(traversingCount);
      duplicateCosmicCount.push_back(cosmicCount);
    }
    auto candidateDirectionSign = [eventSourceSide](Candidate const& candidate) {
      return candidate.physicalDirectionSign != 0
                 ? candidate.physicalDirectionSign
                 : static_cast<int>(shiftDirectionSign(candidate.targetLineState.momentum, eventSourceSide));
    };

    // Build an optional one-to-one MC association only after cleaning.  This
    // block cannot affect reconstruction or the duplicate decision, and an
    // absent collection (as in collision data) leaves every index at -1.
    std::vector<int> genPartIdx(selected.size(), -1);
    std::vector<float> genPartDeltaR(selected.size(), -1.f);
    if (genParticles.isValid()) {
      struct Match {
        double deltaR;
        unsigned int selectedIndex;
        unsigned int genIndex;
      };
      std::vector<Match> matches;
      for (unsigned int selectedIndex = 0; selectedIndex < selected.size(); ++selectedIndex)
        for (unsigned int genIndex = 0; genIndex < genParticles->size(); ++genIndex) {
          auto const& particle = (*genParticles)[genIndex];
          if (std::abs(particle.pdgId()) != 13 || particle.status() != 1)
            continue;
          double const deltaR = directionDeltaR(selected[selectedIndex]->targetLineState.momentum, particle);
          if (deltaR < maxGenDeltaR_)
            matches.push_back({deltaR, selectedIndex, genIndex});
        }
      // This is a diagnostic association, not part of reconstruction.  Allow
      // several reconstructed candidates to reference the same primary muon;
      // one-to-one assignment hid overlaps by leaving all but one at -1.
      for (auto const& match : matches)
        if (genPartIdx[match.selectedIndex] < 0 || match.deltaR < genPartDeltaR[match.selectedIndex]) {
          genPartIdx[match.selectedIndex] = match.genIndex;
          genPartDeltaR[match.selectedIndex] = match.deltaR;
        }
    }

    // MC-only, reconstruction-neutral truth closure diagnostics. Associate a
    // selected reconstructed row to the primary SimTrack through its already
    // established GenPart match, then retain the first/last precision SimHit
    // momenta and propagate the true first-hit state with both target
    // direction contracts. Missing simulation products leave sentinels and
    // therefore keep this producer usable on data.
    std::vector<int> simTruthMatched(selected.size(), 0), simTrackId(selected.size(), -1);
    std::vector<int> simPixelHits(selected.size(), 0), simStripHits(selected.size(), 0),
        simDTHits(selected.size(), 0), simCSCHits(selected.size(), 0),
        simRPCHits(selected.size(), 0), simGEMHits(selected.size(), 0), simMuonDetectorMask(selected.size(), 0);
    std::vector<int> simHBHEHits(selected.size(), 0), simHFHits(selected.size(), 0), simHOHits(selected.size(), 0),
        nAddedDTTruthChamberMatches(selected.size(), -1);
    std::vector<int> simHcalHits(selected.size(), enableHcalDiagnostics_ ? -1 : -2),
        simZDCHits(selected.size(), enableZDCDiagnostics_ ? -1 : -2);
    std::vector<float> simTrackP(selected.size(), -1.f), simFirstPrecisionHitP(selected.size(), -1.f),
        simLastPrecisionHitP(selected.size(), -1.f), simLossToFirstPrecisionHit(selected.size(), -1.f),
        simLossAcrossPrecisionHits(selected.size(), -1.f), simFirstPrecisionPath(selected.size(), -1.f);
    std::vector<float> simHcalEnergy(selected.size(), enableHcalDiagnostics_ ? -1.f : -2.f),
        simHcalFirstTime(selected.size(), 0.f), simHcalLastTime(selected.size(), 0.f),
        simZDCEnergy(selected.size(), enableZDCDiagnostics_ ? -1.f : -2.f),
        simZDCFirstTime(selected.size(), 0.f), simZDCLastTime(selected.size(), 0.f);
    if (produceMomentumClosureDiagnostics_ && genParticles.isValid() && simTracks.isValid() &&
        simVertices.isValid()) {
      std::unordered_map<unsigned int, std::vector<PSimHit const*>> precisionHitsByTrack;
      std::unordered_map<unsigned int, std::array<unsigned int, 4>> detectorHitCounts;
      std::unordered_map<unsigned int, std::array<unsigned int, 2>> trackerHitCounts;
      std::unordered_map<unsigned int, std::set<uint32_t>> dtChambersByTrack;
      auto collectHits = [&precisionHitsByTrack, &detectorHitCounts](auto const& handle,
                                                                    unsigned int detectorIndex,
                                                                    bool precision) {
        if (!handle.isValid())
          return;
        for (auto const& hit : *handle)
          if (std::abs(hit.particleType()) == 13) {
            ++detectorHitCounts[hit.trackId()][detectorIndex];
            if (precision)
              precisionHitsByTrack[hit.trackId()].push_back(&hit);
          }
      };
      collectHits(dtSimHits, 0, true);
      collectHits(cscSimHits, 1, true);
      collectHits(rpcSimHits, 2, false);
      collectHits(gemSimHits, 3, true);
      if (dtSimHits.isValid())
        for (auto const& hit : *dtSimHits)
          if (std::abs(hit.particleType()) == 13)
            dtChambersByTrack[hit.trackId()].insert(DTWireId(hit.detUnitId()).chamberId().rawId());
      for (auto const& handle : trackerSimHits) {
        if (!handle.isValid())
          continue;
        for (auto const& hit : *handle) {
          if (std::abs(hit.particleType()) != 13)
            continue;
          int const subdetector = DetId(hit.detUnitId()).subdetId();
          unsigned int const index = subdetector <= 2 ? 0 : 1;
          ++trackerHitCounts[hit.trackId()][index];
        }
      }

      struct CaloHitSummary {
        int hits = 0;
        std::array<int, 3> hcalSubdetectorHits{{0, 0, 0}};
        float energy = 0.f;
        float firstTime = std::numeric_limits<float>::infinity();
        float lastTime = -std::numeric_limits<float>::infinity();
      };
      std::unordered_map<unsigned int, CaloHitSummary> hcalHitsByTrack;
      std::unordered_map<unsigned int, CaloHitSummary> zdcHitsByTrack;
      auto collectCaloHits = [](auto const& handle, auto& summaries, bool classifyHcalSubdetector) {
        if (!handle.isValid())
          return;
        for (auto const& hit : *handle) {
          auto& summary = summaries[hit.geantTrackId()];
          ++summary.hits;
          if (classifyHcalSubdetector) {
            int subdetector = 0, zside = 0, depth = 0, eta = 0, phi = 0, layer = 0;
            HcalTestNumbering::unpackHcalIndex(
                hit.id(), subdetector, zside, depth, eta, phi, layer);
            if (subdetector == HcalBarrel || subdetector == HcalEndcap)
              ++summary.hcalSubdetectorHits[0];
            else if (subdetector == HcalForward)
              ++summary.hcalSubdetectorHits[1];
            else if (subdetector == HcalOuter)
              ++summary.hcalSubdetectorHits[2];
          }
          summary.energy += hit.energy();
          summary.firstTime = std::min(summary.firstTime, static_cast<float>(hit.time()));
          summary.lastTime = std::max(summary.lastTime, static_cast<float>(hit.time()));
        }
      };
      if (enableHcalDiagnostics_)
        collectCaloHits(hcalSimHits, hcalHitsByTrack, true);
      if (enableZDCDiagnostics_)
        collectCaloHits(zdcSimHits, zdcHitsByTrack, false);

      for (unsigned int selectedIndex = 0; selectedIndex < selected.size(); ++selectedIndex) {
        int const genIndex = genPartIdx[selectedIndex];
        if (genIndex < 0)
          continue;
        auto const& gen = (*genParticles)[genIndex];
        SimTrack const* matchedSimTrack = nullptr;
        double bestScore = std::numeric_limits<double>::infinity();
        for (auto const& simTrack : *simTracks) {
          if (simTrack.type() != gen.pdgId() || simTrack.vertIndex() < 0 ||
              static_cast<std::size_t>(simTrack.vertIndex()) >= simVertices->size() ||
              !(simTrack.momentum().P() > 0.) || !(gen.p() > 0.))
            continue;
          GlobalVector const simDirection(
              simTrack.momentum().px(), simTrack.momentum().py(), simTrack.momentum().pz());
          GlobalVector const genDirection(gen.px(), gen.py(), gen.pz());
          double const directionDot = simDirection.unit().dot(genDirection.unit());
          double const angle = std::acos(std::clamp(directionDot, -1., 1.));
          double const relativeMomentum = std::abs(std::log(simTrack.momentum().P() / gen.p()));
          double const score = angle + relativeMomentum;
          if (score < bestScore) {
            bestScore = score;
            matchedSimTrack = &simTrack;
          }
        }
        if (!matchedSimTrack || bestScore > 0.1)
          continue;
        simTruthMatched[selectedIndex] = 1;
        simTrackId[selectedIndex] = matchedSimTrack->trackId();
        simTrackP[selectedIndex] = matchedSimTrack->momentum().P();
        auto const trackerCounts = trackerHitCounts[matchedSimTrack->trackId()];
        simPixelHits[selectedIndex] = trackerCounts[0];
        simStripHits[selectedIndex] = trackerCounts[1];
        auto const detectorCounts = detectorHitCounts[matchedSimTrack->trackId()];
        simDTHits[selectedIndex] = detectorCounts[0];
        simCSCHits[selectedIndex] = detectorCounts[1];
        simRPCHits[selectedIndex] = detectorCounts[2];
        simGEMHits[selectedIndex] = detectorCounts[3];
        nAddedDTTruthChamberMatches[selectedIndex] = 0;
        auto const truthDTChambers = dtChambersByTrack.find(matchedSimTrack->trackId());
        if (truthDTChambers != dtChambersByTrack.end())
          for (auto const chamber : selected[selectedIndex]->addedDTChamberIds)
            nAddedDTTruthChamberMatches[selectedIndex] += truthDTChambers->second.count(chamber);
        auto assignCalo = [selectedIndex](bool enabled,
                                          bool available,
                                          auto const& summaries,
                                          unsigned int trackId,
                                          auto& counts,
                                          auto& energies,
                                          auto& firstTimes,
                                          auto& lastTimes) {
          if (!enabled || !available)
            return;
          counts[selectedIndex] = 0;
          energies[selectedIndex] = 0.f;
          auto const found = summaries.find(trackId);
          if (found == summaries.end())
            return;
          counts[selectedIndex] = found->second.hits;
          energies[selectedIndex] = found->second.energy;
          firstTimes[selectedIndex] = found->second.firstTime;
          lastTimes[selectedIndex] = found->second.lastTime;
        };
        assignCalo(enableHcalDiagnostics_,
                   hcalSimHits.isValid(),
                   hcalHitsByTrack,
                   matchedSimTrack->trackId(),
                   simHcalHits,
                   simHcalEnergy,
                   simHcalFirstTime,
                   simHcalLastTime);
        auto const hcalSummary = hcalHitsByTrack.find(matchedSimTrack->trackId());
        if (hcalSummary != hcalHitsByTrack.end()) {
          simHBHEHits[selectedIndex] = hcalSummary->second.hcalSubdetectorHits[0];
          simHFHits[selectedIndex] = hcalSummary->second.hcalSubdetectorHits[1];
          simHOHits[selectedIndex] = hcalSummary->second.hcalSubdetectorHits[2];
        }
        assignCalo(enableZDCDiagnostics_,
                   zdcSimHits.isValid(),
                   zdcHitsByTrack,
                   matchedSimTrack->trackId(),
                   simZDCHits,
                   simZDCEnergy,
                   simZDCFirstTime,
                   simZDCLastTime);
        for (unsigned int detectorIndex = 0; detectorIndex < detectorCounts.size(); ++detectorIndex)
          if (detectorCounts[detectorIndex] > 0)
            simMuonDetectorMask[selectedIndex] |= 1u << detectorIndex;
        auto hitIt = precisionHitsByTrack.find(matchedSimTrack->trackId());
        if (hitIt == precisionHitsByTrack.end() || hitIt->second.empty())
          continue;
        auto const& simVertex = (*simVertices)[matchedSimTrack->vertIndex()];
        GlobalPoint const vertexPosition(simVertex.position().x(), simVertex.position().y(), simVertex.position().z());
        GlobalVector const truthMomentum(matchedSimTrack->momentum().px(),
                                         matchedSimTrack->momentum().py(),
                                         matchedSimTrack->momentum().pz());
        auto pathFromVertex = [&trackingGeometry, &vertexPosition, &truthMomentum](PSimHit const* hit) {
          auto const* det = trackingGeometry.idToDetUnit(DetId(hit->detUnitId()));
          if (!det)
            return std::numeric_limits<double>::infinity();
          return static_cast<double>(
              (det->surface().toGlobal(hit->entryPoint()) - vertexPosition).dot(truthMomentum.unit()));
        };
        auto& hits = hitIt->second;
        std::sort(hits.begin(), hits.end(), [&pathFromVertex](PSimHit const* first, PSimHit const* second) {
          return pathFromVertex(first) < pathFromVertex(second);
        });
        auto const* firstHit = hits.front();
        auto const* lastHit = hits.back();
        auto const* firstDet = trackingGeometry.idToDetUnit(DetId(firstHit->detUnitId()));
        if (!firstDet || !(firstHit->pabs() > 0.) || !(lastHit->pabs() > 0.))
          continue;
        simFirstPrecisionHitP[selectedIndex] = firstHit->pabs();
        simLastPrecisionHitP[selectedIndex] = lastHit->pabs();
        simLossToFirstPrecisionHit[selectedIndex] =
            (matchedSimTrack->momentum().P() - firstHit->pabs()) / matchedSimTrack->momentum().P();
        simLossAcrossPrecisionHits[selectedIndex] = (firstHit->pabs() - lastHit->pabs()) / firstHit->pabs();
        simFirstPrecisionPath[selectedIndex] = pathFromVertex(firstHit);
      }
    }

    std::vector<float> pt, eta, phi, mass, p, px, py, pz, trackPt, innerPt, outerPt, upstreamPt, preRefitPt,
        trackerPt, trackerEta, trackerPhi, trackerP, trackerPz, trackerChi2, trackerNdof,
        combinedPt, combinedEta, combinedPhi, combinedP, combinedPz, combinedChi2, combinedNdof,
        combinedTargetPt, combinedTargetPz, combinedTargetDca,
        constrainedPt, constrainedEta, constrainedPhi, constrainedMass, constrainedP, constrainedPx,
        constrainedPy, constrainedPz, constrainedVx, constrainedVy, constrainedVz, constrainedChi2,
        constrainedNdof, constrainedTargetChi2,
        preRefitPz, directionalRefitUpstreamPt, directionalRefitUpstreamP, directionalRefitDownstreamP,
        directionalRefitFractionalLossAcrossHits, directionalRefitSourceFacingSurfaceDistance,
        directionalRefitFilterDirectionCosine, directionalRefitTargetDirectionCosine,
        directionalRefitChi2, directionalRefitNdof,
        directionalRefitFirstChi2, directionalRefitFirstNdof, directionalRefitFirstUpstreamPt,
        directionalRefitFirstQoverP,
        directionalRefitFirstTargetPt, directionalRefitFirstTargetPz, directionalRefitSecondChi2,
        directionalRefitSecondNdof, directionalRefitSecondUpstreamPt, directionalRefitSecondQoverP,
        directionalRefitSecondTargetPt,
        directionalRefitSecondTargetPz, directionalRefitRelativeQoverPChange, directionalRefitVacuumTargetPt,
        directionalRefitVacuumTargetPz, directionalRefitMaterialDeltaPt,
        directionalRefitMaterialBoundaryPt, directionalRefitMaterialBoundaryPz, directionalRefitMaterialPath,
        directionalRefitAllHitsOrderingSpan, directionalRefitPrecisionOrderingSpan, precisionRefitLeverArm,
        directionalRefitAllHitsChi2, directionalRefitAllHitsNdof, directionalRefitAllHitsTargetPt,
        directionalRefitAllHitsTargetPz, directionalRefitPrecisionChi2, directionalRefitPrecisionNdof,
        directionalRefitPrecisionUpstreamPt, directionalRefitPrecisionQoverP, directionalRefitPrecisionTargetPt,
        directionalRefitPrecisionTargetPz, directionalRefitPrecisionRelativeQoverPChange,
        directionalRefitPrecisionFirstTargetPt, directionalRefitPrecisionFirstTargetPz,
        directionalRefitPrecisionSecondTargetPt, directionalRefitPrecisionSecondTargetPz,
        directionalRefitPrecisionRelativeToAllQoverP, directionalRefitPrecisionTargetDca,
        directionalRefitSourceLegPrecisionQoverP, directionalRefitSourceLegPrecisionTargetPt,
        directionalRefitSourceLegPrecisionTargetPz, directionalRefitOppositeLegPrecisionQoverP,
        directionalRefitOppositeLegPrecisionTargetPt, directionalRefitOppositeLegPrecisionTargetPz,
        ptError, etaError, phiError, vx, vy, vz,
        trackVx, trackVy, trackVz, dxy, dz, innerR, innerZ, outerR, outerZ, chi2, ndof, normalizedChi2, linePcaR, linePcaZ,
        chordLinePcaR, chordLinePcaZ, targetLinePath, timingChi2, timingDeltaChi2,
        entryX, entryY, entryZ, entryR, exitX, exitY, exitZ, exitR, hitSpan, hitDeltaZ,
        minDTResidual, minDTEstimatorChi2, minTrackerResidual, minTrackerEstimatorChi2,
        matchedEcalEnergy, matchedHBHEEnergy, matchedHFEnergy, matchedHOEnergy, matchedZDCEnergy,
        matchedEcalTime, matchedHBHETime, matchedHFTime, matchedHOTime, matchedZDCTime,
        minEcalLineDistance, minHBHELineDistance, minHFLineDistance, minHOLineDistance, minZDCLineDistance,
        caloTimingDeltaChi2;
    std::vector<int> charge, sourceIndex, recoAlgorithm, topology, orientationValid, validHits, validMuonHits,
        trackerMatchValid, trackerTrackIndex, trackerValidHits, trackerPixelHits, trackerStripHits,
        combinedTrackValid, combinedValidHits, combinedTrackerHits, combinedMuonHits,
        muonStations, lostHits, directionFlipped, quality, constrainedValid, constrainedHits, constrainedStatus,
        inferredSourceSide, chargeMatchesGen, timingDirectionSign, nTimingMeasurements, directionalRefitAttempted,
        directionalRefitValid, directionalRefitHits, directionalRefitFirstValid, directionalRefitFirstHits,
        directionalRefitSecondValid, directionalRefitSecondHits, directionalRefitSecondConverged,
        directionalRefitSelectedIteration, nDTRefitHits, nCSCRefitHits, nRPCRefitHits, nGEMRefitHits,
        nCompatibleDTSegments, nCompatiblePixelHits, nCompatibleStripHits,
        nPropagatedDTSegments, nPropagatedPixelHits, nPropagatedStripHits,
        nMatchedEcalRecHits, nMatchedHBHERecHits, nMatchedHFRecHits, nMatchedHORecHits, nMatchedZDCRecHits,
        caloTimingDirectionSign, nCaloTimingMeasurements,
        nAddedDTRefitHits, nAddedTrackerRefitHits,
        nPrecisionRefitStations, directionalRefitUsedPrecisionHits, directionalRefitAllHitsValid,
        directionalRefitAllHits, directionalRefitAllHitsSelectedIteration, directionalRefitPrecisionValid,
        directionalRefitPrecisionHits, directionalRefitPrecisionSecondValid,
        directionalRefitPrecisionSecondConverged, directionalRefitPrecisionSelectedIteration,
        directionalRefitMaterialBoundaryValid, directionalRefitAllHitsInput,
        directionalRefitAllHitsOrderingFallback, directionalRefitAllHitsFirstRejected,
        directionalRefitAllHitsSecondRejected, directionalRefitPrecisionInput,
        directionalRefitPrecisionOrderingFallback, directionalRefitPrecisionFirstRejected,
        directionalRefitPrecisionSecondRejected, directionalRefitSourceLegPrecisionValid,
        directionalRefitSourceLegPrecisionHits, directionalRefitOppositeLegPrecisionValid,
        directionalRefitOppositeLegPrecisionHits, nDTHits, nCSCHits, nRPCHits, nGEMHits, nME0Hits,
        nHitsPlusEndcap, nHitsBarrel, nHitsMinusEndcap, nHitsNearEndcap, nHitsFarEndcap,
        nStationsNearEndcap, nStationsBarrel, nStationsFarEndcap, detectorMask,
        entryExitValid, entryRegion, exitRegion, entrySubdetector, exitSubdetector;
    for (unsigned int selectedIndex = 0; selectedIndex < selected.size(); ++selectedIndex) {
      auto const* candidate = selected[selectedIndex];
      auto const& track = *candidate->track;
      auto const& propagated = candidate->targetLineState;
      auto const [pcaR, pcaZ] = transverseLinePca(track);
      auto const [chordPcaR, chordPcaZ] = chordLinePca(track);
      // Cosmic-style fits do not determine which way along the fitted helix
      // the particle travelled.  Orient it from the inferred source side
      // toward CMS, supporting both +z and -z SHIFT locations without truth.
      double const sign = candidateDirectionSign(*candidate);
      bool const flip = sign < 0.;
      auto const geometryDiagnostics =
          muonGeometryDiagnostics(track, trackingGeometry, eventSourceSide, sign);

      // Associate the selected detector-only row to the independent CMSSW
      // cosmic-global fit made from the same source collection.  These values
      // are diagnostics only: they never enter cleaning, topology, momentum,
      // or dimuon selection in this pass.
      reco::MuonTrackLinksCollection const* sourceLinks = nullptr;
      if (candidate->source == 0 && dsaGlobalLinks.isValid())
        sourceLinks = dsaGlobalLinks.product();
      else if (candidate->source == 1 && traversingGlobalLinks.isValid())
        sourceLinks = traversingGlobalLinks.product();
      else if (candidate->source == 2 && cosmicGlobalLinks.isValid())
        sourceLinks = cosmicGlobalLinks.product();
      reco::MuonTrackLinks const* matchedLink = nullptr;
      if (sourceLinks)
        for (auto const& link : *sourceLinks)
          if (link.standAloneTrack().isNonnull() && link.standAloneTrack().key() == candidate->sourceIndex) {
            matchedLink = &link;
            break;
          }
      bool const hasTracker = matchedLink && matchedLink->trackerTrack().isNonnull();
      bool const hasCombined = matchedLink && matchedLink->globalTrack().isNonnull();
      trackerMatchValid.push_back(hasTracker);
      trackerTrackIndex.push_back(hasTracker ? static_cast<int>(matchedLink->trackerTrack().key()) : -1);
      trackerValidHits.push_back(hasTracker ? matchedLink->trackerTrack()->hitPattern().numberOfValidTrackerHits() : 0);
      trackerPixelHits.push_back(hasTracker ? matchedLink->trackerTrack()->hitPattern().numberOfValidPixelHits() : 0);
      trackerStripHits.push_back(hasTracker ? matchedLink->trackerTrack()->hitPattern().numberOfValidStripHits() : 0);
      trackerPt.push_back(hasTracker ? matchedLink->trackerTrack()->pt() : 0.);
      trackerEta.push_back(hasTracker ? matchedLink->trackerTrack()->eta() : 0.);
      trackerPhi.push_back(hasTracker ? matchedLink->trackerTrack()->phi() : 0.);
      trackerP.push_back(hasTracker ? matchedLink->trackerTrack()->p() : 0.);
      trackerPz.push_back(hasTracker ? matchedLink->trackerTrack()->pz() : 0.);
      trackerChi2.push_back(hasTracker ? matchedLink->trackerTrack()->chi2() : 0.);
      trackerNdof.push_back(hasTracker ? matchedLink->trackerTrack()->ndof() : 0.);
      combinedTrackValid.push_back(hasCombined);
      combinedValidHits.push_back(hasCombined ? matchedLink->globalTrack()->numberOfValidHits() : 0);
      combinedTrackerHits.push_back(
          hasCombined ? matchedLink->globalTrack()->hitPattern().numberOfValidTrackerHits() : 0);
      combinedMuonHits.push_back(hasCombined ? matchedLink->globalTrack()->hitPattern().numberOfValidMuonHits() : 0);
      combinedPt.push_back(hasCombined ? matchedLink->globalTrack()->pt() : 0.);
      combinedEta.push_back(hasCombined ? matchedLink->globalTrack()->eta() : 0.);
      combinedPhi.push_back(hasCombined ? matchedLink->globalTrack()->phi() : 0.);
      combinedP.push_back(hasCombined ? matchedLink->globalTrack()->p() : 0.);
      combinedPz.push_back(hasCombined ? matchedLink->globalTrack()->pz() : 0.);
      combinedChi2.push_back(hasCombined ? matchedLink->globalTrack()->chi2() : 0.);
      combinedNdof.push_back(hasCombined ? matchedLink->globalTrack()->ndof() : 0.);
      PropagatedState combinedTarget;
      if (hasCombined)
        combinedTarget = propagateToTargetLine(
            *matchedLink->globalTrack(), vacuumPropagator, sign, preRefitMaterialPropagator, eventSourceSide);
      combinedTargetPt.push_back(combinedTarget.valid ? combinedTarget.momentum.perp() : 0.);
      combinedTargetPz.push_back(combinedTarget.valid ? combinedTarget.momentum.z() : 0.);
      combinedTargetDca.push_back(combinedTarget.valid ? combinedTarget.position.perp() : -1.);
      // The propagated state was already canonicalised before transport.
      double const storedPx = propagated.momentum.x();
      double const storedPy = propagated.momentum.y();
      double const storedPz = propagated.momentum.z();
      double const storedPt = std::hypot(storedPx, storedPy);
      double const storedP = propagated.momentum.mag();
      pt.push_back(storedPt);
      eta.push_back(storedPt > 0. ? std::asinh(storedPz / storedPt)
                                  : std::copysign(std::numeric_limits<float>::infinity(), storedPz));
      phi.push_back(std::atan2(storedPy, storedPx));
      mass.push_back(0.105658f);
      p.push_back(storedP);
      px.push_back(storedPx);
      py.push_back(storedPy);
      pz.push_back(storedPz);
      auto const& constrained = candidate->constrainedState;
      double const constrainedPxValue = candidate->constrainedValid ? constrained.momentum.x() : 0.;
      double const constrainedPyValue = candidate->constrainedValid ? constrained.momentum.y() : 0.;
      double const constrainedPzValue = candidate->constrainedValid ? constrained.momentum.z() : 0.;
      double const constrainedPtValue = std::hypot(constrainedPxValue, constrainedPyValue);
      constrainedValid.push_back(candidate->constrainedValid);
      constrainedHits.push_back(candidate->constrainedHits);
      constrainedStatus.push_back(candidate->constrainedStatus);
      constrainedPt.push_back(constrainedPtValue);
      constrainedEta.push_back(candidate->constrainedValid && constrainedPtValue > 0.
                                   ? std::asinh(constrainedPzValue / constrainedPtValue)
                                   : 0.);
      constrainedPhi.push_back(candidate->constrainedValid
                                   ? std::atan2(constrainedPyValue, constrainedPxValue)
                                   : 0.);
      constrainedMass.push_back(0.105658f);
      constrainedP.push_back(candidate->constrainedValid ? constrained.momentum.mag() : 0.);
      constrainedPx.push_back(constrainedPxValue);
      constrainedPy.push_back(constrainedPyValue);
      constrainedPz.push_back(constrainedPzValue);
      constrainedVx.push_back(candidate->constrainedValid ? constrained.position.x() : 0.);
      constrainedVy.push_back(candidate->constrainedValid ? constrained.position.y() : 0.);
      constrainedVz.push_back(candidate->constrainedValid ? constrained.position.z() : 0.);
      constrainedChi2.push_back(candidate->constrainedChi2);
      constrainedNdof.push_back(candidate->constrainedNdof);
      constrainedTargetChi2.push_back(candidate->constrainedTargetChi2);
      trackPt.push_back(track.pt());
      innerPt.push_back(track.innerMomentum().rho());
      outerPt.push_back(track.outerMomentum().rho());
      auto const endpointDelta = track.outerPosition() - track.innerPosition();
      bool const upstreamIsOuter = sign * track.innerMomentum().Dot(endpointDelta) < 0.;
      upstreamPt.push_back(upstreamIsOuter ? track.outerMomentum().rho() : track.innerMomentum().rho());
      preRefitPt.push_back(candidate->preRefitPt);
      preRefitPz.push_back(candidate->preRefitPz);
      directionalRefitUpstreamPt.push_back(candidate->directionalRefitUpstreamPt);
      directionalRefitUpstreamP.push_back(candidate->directionalRefitUpstreamP);
      directionalRefitDownstreamP.push_back(candidate->directionalRefitDownstreamP);
      directionalRefitFractionalLossAcrossHits.push_back(candidate->directionalRefitFractionalLossAcrossHits);
      directionalRefitSourceFacingSurfaceDistance.push_back(candidate->directionalRefitSourceFacingSurfaceDistance);
      directionalRefitFilterDirectionCosine.push_back(candidate->directionalRefitFilterDirectionCosine);
      directionalRefitTargetDirectionCosine.push_back(candidate->directionalRefitTargetDirectionCosine);
      directionalRefitAttempted.push_back(candidate->directionalRefitAttempted);
      directionalRefitValid.push_back(candidate->directionalRefitValid);
      directionalRefitHits.push_back(candidate->directionalRefitHits);
      directionalRefitChi2.push_back(candidate->directionalRefitChi2);
      directionalRefitNdof.push_back(candidate->directionalRefitNdof);
      directionalRefitFirstValid.push_back(candidate->directionalRefitFirstValid);
      directionalRefitFirstHits.push_back(candidate->directionalRefitFirstHits);
      directionalRefitFirstChi2.push_back(candidate->directionalRefitFirstChi2);
      directionalRefitFirstNdof.push_back(candidate->directionalRefitFirstNdof);
      directionalRefitFirstUpstreamPt.push_back(candidate->directionalRefitFirstUpstreamPt);
      directionalRefitFirstQoverP.push_back(candidate->directionalRefitFirstQoverP);
      directionalRefitFirstTargetPt.push_back(candidate->directionalRefitFirstTargetPt);
      directionalRefitFirstTargetPz.push_back(candidate->directionalRefitFirstTargetPz);
      directionalRefitSecondValid.push_back(candidate->directionalRefitSecondValid);
      directionalRefitSecondHits.push_back(candidate->directionalRefitSecondHits);
      directionalRefitSecondChi2.push_back(candidate->directionalRefitSecondChi2);
      directionalRefitSecondNdof.push_back(candidate->directionalRefitSecondNdof);
      directionalRefitSecondUpstreamPt.push_back(candidate->directionalRefitSecondUpstreamPt);
      directionalRefitSecondQoverP.push_back(candidate->directionalRefitSecondQoverP);
      directionalRefitSecondTargetPt.push_back(candidate->directionalRefitSecondTargetPt);
      directionalRefitSecondTargetPz.push_back(candidate->directionalRefitSecondTargetPz);
      directionalRefitSecondConverged.push_back(candidate->directionalRefitSecondConverged);
      directionalRefitRelativeQoverPChange.push_back(candidate->directionalRefitRelativeQoverPChange);
      directionalRefitSelectedIteration.push_back(candidate->directionalRefitSelectedIteration);
      directionalRefitVacuumTargetPt.push_back(candidate->directionalRefitVacuumTargetPt);
      directionalRefitVacuumTargetPz.push_back(candidate->directionalRefitVacuumTargetPz);
      directionalRefitMaterialDeltaPt.push_back(candidate->directionalRefitMaterialDeltaPt);
      directionalRefitMaterialBoundaryValid.push_back(candidate->directionalRefitMaterialBoundaryValid);
      directionalRefitMaterialBoundaryPt.push_back(candidate->directionalRefitMaterialBoundaryPt);
      directionalRefitMaterialBoundaryPz.push_back(candidate->directionalRefitMaterialBoundaryPz);
      directionalRefitMaterialPath.push_back(candidate->directionalRefitMaterialPath);
      directionalRefitAllHitsInput.push_back(candidate->directionalRefitAllHitsInput);
      directionalRefitAllHitsOrderingFallback.push_back(candidate->directionalRefitAllHitsOrderingFallback);
      directionalRefitAllHitsOrderingSpan.push_back(candidate->directionalRefitAllHitsOrderingSpan);
      directionalRefitAllHitsFirstRejected.push_back(candidate->directionalRefitAllHitsFirstRejected);
      directionalRefitAllHitsSecondRejected.push_back(candidate->directionalRefitAllHitsSecondRejected);
      directionalRefitPrecisionInput.push_back(candidate->directionalRefitPrecisionInput);
      directionalRefitPrecisionOrderingFallback.push_back(candidate->directionalRefitPrecisionOrderingFallback);
      directionalRefitPrecisionOrderingSpan.push_back(candidate->directionalRefitPrecisionOrderingSpan);
      directionalRefitPrecisionFirstRejected.push_back(candidate->directionalRefitPrecisionFirstRejected);
      directionalRefitPrecisionSecondRejected.push_back(candidate->directionalRefitPrecisionSecondRejected);
      nDTRefitHits.push_back(candidate->nDTRefitHits);
      nCSCRefitHits.push_back(candidate->nCSCRefitHits);
      nRPCRefitHits.push_back(candidate->nRPCRefitHits);
      nGEMRefitHits.push_back(candidate->nGEMRefitHits);
      nCompatibleDTSegments.push_back(candidate->nCompatibleDTSegments);
      nCompatiblePixelHits.push_back(candidate->nCompatiblePixelHits);
      nCompatibleStripHits.push_back(candidate->nCompatibleStripHits);
      nPropagatedDTSegments.push_back(candidate->nPropagatedDTSegments);
      nPropagatedPixelHits.push_back(candidate->nPropagatedPixelHits);
      nPropagatedStripHits.push_back(candidate->nPropagatedStripHits);
      minDTResidual.push_back(candidate->minDTResidual);
      minDTEstimatorChi2.push_back(candidate->minDTEstimatorChi2);
      minTrackerResidual.push_back(candidate->minTrackerResidual);
      minTrackerEstimatorChi2.push_back(candidate->minTrackerEstimatorChi2);
      nMatchedEcalRecHits.push_back(candidate->nMatchedEcalRecHits);
      nMatchedHBHERecHits.push_back(candidate->nMatchedHBHERecHits);
      nMatchedHFRecHits.push_back(candidate->nMatchedHFRecHits);
      nMatchedHORecHits.push_back(candidate->nMatchedHORecHits);
      nMatchedZDCRecHits.push_back(candidate->nMatchedZDCRecHits);
      matchedEcalEnergy.push_back(candidate->matchedEcalEnergy);
      matchedHBHEEnergy.push_back(candidate->matchedHBHEEnergy);
      matchedHFEnergy.push_back(candidate->matchedHFEnergy);
      matchedHOEnergy.push_back(candidate->matchedHOEnergy);
      matchedZDCEnergy.push_back(candidate->matchedZDCEnergy);
      matchedEcalTime.push_back(candidate->matchedEcalTime);
      matchedHBHETime.push_back(candidate->matchedHBHETime);
      matchedHFTime.push_back(candidate->matchedHFTime);
      matchedHOTime.push_back(candidate->matchedHOTime);
      matchedZDCTime.push_back(candidate->matchedZDCTime);
      minEcalLineDistance.push_back(candidate->minEcalLineDistance);
      minHBHELineDistance.push_back(candidate->minHBHELineDistance);
      minHFLineDistance.push_back(candidate->minHFLineDistance);
      minHOLineDistance.push_back(candidate->minHOLineDistance);
      minZDCLineDistance.push_back(candidate->minZDCLineDistance);
      caloTimingDirectionSign.push_back(candidate->caloTimingDirectionSign);
      nCaloTimingMeasurements.push_back(candidate->nCaloTimingMeasurements);
      caloTimingDeltaChi2.push_back(candidate->caloTimingDeltaChi2);
      nAddedDTRefitHits.push_back(candidate->nAddedDTRefitHits);
      nAddedTrackerRefitHits.push_back(candidate->nAddedTrackerRefitHits);
      nPrecisionRefitStations.push_back(candidate->nPrecisionRefitStations);
      precisionRefitLeverArm.push_back(candidate->precisionRefitLeverArm);
      directionalRefitUsedPrecisionHits.push_back(candidate->directionalRefitUsedPrecisionHits);
      directionalRefitAllHitsValid.push_back(candidate->directionalRefitAllHitsValid);
      directionalRefitAllHits.push_back(candidate->directionalRefitAllHits);
      directionalRefitAllHitsChi2.push_back(candidate->directionalRefitAllHitsChi2);
      directionalRefitAllHitsNdof.push_back(candidate->directionalRefitAllHitsNdof);
      directionalRefitAllHitsTargetPt.push_back(candidate->directionalRefitAllHitsTargetPt);
      directionalRefitAllHitsTargetPz.push_back(candidate->directionalRefitAllHitsTargetPz);
      directionalRefitAllHitsSelectedIteration.push_back(candidate->directionalRefitAllHitsSelectedIteration);
      directionalRefitPrecisionValid.push_back(candidate->directionalRefitPrecisionValid);
      directionalRefitPrecisionHits.push_back(candidate->directionalRefitPrecisionHits);
      directionalRefitPrecisionChi2.push_back(candidate->directionalRefitPrecisionChi2);
      directionalRefitPrecisionNdof.push_back(candidate->directionalRefitPrecisionNdof);
      directionalRefitPrecisionUpstreamPt.push_back(candidate->directionalRefitPrecisionUpstreamPt);
      directionalRefitPrecisionQoverP.push_back(candidate->directionalRefitPrecisionQoverP);
      directionalRefitPrecisionTargetPt.push_back(candidate->directionalRefitPrecisionTargetPt);
      directionalRefitPrecisionTargetPz.push_back(candidate->directionalRefitPrecisionTargetPz);
      directionalRefitPrecisionSecondValid.push_back(candidate->directionalRefitPrecisionSecondValid);
      directionalRefitPrecisionSecondConverged.push_back(candidate->directionalRefitPrecisionSecondConverged);
      directionalRefitPrecisionRelativeQoverPChange.push_back(
          candidate->directionalRefitPrecisionRelativeQoverPChange);
      directionalRefitPrecisionSelectedIteration.push_back(candidate->directionalRefitPrecisionSelectedIteration);
      directionalRefitPrecisionFirstTargetPt.push_back(candidate->directionalRefitPrecisionFirstTargetPt);
      directionalRefitPrecisionFirstTargetPz.push_back(candidate->directionalRefitPrecisionFirstTargetPz);
      directionalRefitPrecisionSecondTargetPt.push_back(candidate->directionalRefitPrecisionSecondTargetPt);
      directionalRefitPrecisionSecondTargetPz.push_back(candidate->directionalRefitPrecisionSecondTargetPz);
      directionalRefitPrecisionRelativeToAllQoverP.push_back(
          candidate->directionalRefitPrecisionRelativeToAllQoverP);
      directionalRefitPrecisionTargetDca.push_back(candidate->directionalRefitPrecisionTargetDca);
      directionalRefitSourceLegPrecisionValid.push_back(candidate->directionalRefitSourceLegPrecisionValid);
      directionalRefitSourceLegPrecisionHits.push_back(candidate->directionalRefitSourceLegPrecisionHits);
      directionalRefitSourceLegPrecisionQoverP.push_back(candidate->directionalRefitSourceLegPrecisionQoverP);
      directionalRefitSourceLegPrecisionTargetPt.push_back(candidate->directionalRefitSourceLegPrecisionTargetPt);
      directionalRefitSourceLegPrecisionTargetPz.push_back(candidate->directionalRefitSourceLegPrecisionTargetPz);
      directionalRefitOppositeLegPrecisionValid.push_back(candidate->directionalRefitOppositeLegPrecisionValid);
      directionalRefitOppositeLegPrecisionHits.push_back(candidate->directionalRefitOppositeLegPrecisionHits);
      directionalRefitOppositeLegPrecisionQoverP.push_back(candidate->directionalRefitOppositeLegPrecisionQoverP);
      directionalRefitOppositeLegPrecisionTargetPt.push_back(
          candidate->directionalRefitOppositeLegPrecisionTargetPt);
      directionalRefitOppositeLegPrecisionTargetPz.push_back(
          candidate->directionalRefitOppositeLegPrecisionTargetPz);
      ptError.push_back(track.ptError());
      etaError.push_back(track.etaError());
      phiError.push_back(track.phiError());
      // Momentum and charge are the simultaneous two-fold ambiguity of a
      // no-timing cosmic-style helix fit.  Reverse both to preserve curvature.
      int const storedCharge = sign * track.charge();
      charge.push_back(storedCharge);
      int chargeMatch = -1;
      if (genParticles.isValid() && genPartIdx[selectedIndex] >= 0) {
        int const genCharge = (*genParticles)[genPartIdx[selectedIndex]].pdgId() == 13 ? -1 : 1;
        chargeMatch = storedCharge == genCharge;
      }
      chargeMatchesGen.push_back(chargeMatch);
      directionFlipped.push_back(flip);
      inferredSourceSide.push_back(eventSourceSide);
      timingDirectionSign.push_back(candidate->timingDirectionSign);
      nTimingMeasurements.push_back(candidate->timingMeasurements);
      timingChi2.push_back(candidate->timingChi2);
      timingDeltaChi2.push_back(candidate->timingDeltaChi2);
      vx.push_back(propagated.position.x());
      vy.push_back(propagated.position.y());
      vz.push_back(propagated.position.z());
      trackVx.push_back(track.vx());
      trackVy.push_back(track.vy());
      trackVz.push_back(track.vz());
      dxy.push_back(track.dxy());
      dz.push_back(track.dz());
      innerR.push_back(track.innerPosition().rho());
      innerZ.push_back(track.innerPosition().z());
      outerR.push_back(track.outerPosition().rho());
      outerZ.push_back(track.outerPosition().z());
      nDTHits.push_back(geometryDiagnostics.dt);
      nCSCHits.push_back(geometryDiagnostics.csc);
      nRPCHits.push_back(geometryDiagnostics.rpc);
      nGEMHits.push_back(geometryDiagnostics.gem);
      nME0Hits.push_back(geometryDiagnostics.me0);
      nHitsPlusEndcap.push_back(geometryDiagnostics.plusEndcap);
      nHitsBarrel.push_back(geometryDiagnostics.barrel);
      nHitsMinusEndcap.push_back(geometryDiagnostics.minusEndcap);
      nHitsNearEndcap.push_back(geometryDiagnostics.nearEndcap);
      nHitsFarEndcap.push_back(geometryDiagnostics.farEndcap);
      nStationsNearEndcap.push_back(geometryDiagnostics.nearEndcapStations);
      nStationsBarrel.push_back(geometryDiagnostics.barrelStations);
      nStationsFarEndcap.push_back(geometryDiagnostics.farEndcapStations);
      detectorMask.push_back(geometryDiagnostics.detectorMask);
      entryExitValid.push_back(geometryDiagnostics.endpointsValid);
      entryRegion.push_back(geometryDiagnostics.entryRegion);
      exitRegion.push_back(geometryDiagnostics.exitRegion);
      entrySubdetector.push_back(geometryDiagnostics.entrySubdetector);
      exitSubdetector.push_back(geometryDiagnostics.exitSubdetector);
      entryX.push_back(geometryDiagnostics.endpointsValid ? geometryDiagnostics.entryPosition.x() : 0.);
      entryY.push_back(geometryDiagnostics.endpointsValid ? geometryDiagnostics.entryPosition.y() : 0.);
      entryZ.push_back(geometryDiagnostics.endpointsValid ? geometryDiagnostics.entryPosition.z() : 0.);
      entryR.push_back(geometryDiagnostics.endpointsValid ? geometryDiagnostics.entryPosition.perp() : 0.);
      exitX.push_back(geometryDiagnostics.endpointsValid ? geometryDiagnostics.exitPosition.x() : 0.);
      exitY.push_back(geometryDiagnostics.endpointsValid ? geometryDiagnostics.exitPosition.y() : 0.);
      exitZ.push_back(geometryDiagnostics.endpointsValid ? geometryDiagnostics.exitPosition.z() : 0.);
      exitR.push_back(geometryDiagnostics.endpointsValid ? geometryDiagnostics.exitPosition.perp() : 0.);
      hitSpan.push_back(geometryDiagnostics.hitSpan);
      hitDeltaZ.push_back(geometryDiagnostics.endpointsValid
                              ? std::abs(geometryDiagnostics.exitPosition.z() - geometryDiagnostics.entryPosition.z())
                              : 0.);
      linePcaR.push_back(pcaR);
      linePcaZ.push_back(pcaZ);
      chordLinePcaR.push_back(chordPcaR);
      chordLinePcaZ.push_back(chordPcaZ);
      targetLinePath.push_back(propagated.path);
      chi2.push_back(track.chi2());
      ndof.push_back(track.ndof());
      normalizedChi2.push_back(track.normalizedChi2());
      validHits.push_back(track.numberOfValidHits());
      validMuonHits.push_back(track.hitPattern().numberOfValidMuonHits());
      muonStations.push_back(track.hitPattern().muonStationsWithValidHits());
      lostHits.push_back(track.numberOfLostHits());
      sourceIndex.push_back(candidate->sourceIndex);
      recoAlgorithm.push_back(candidate->source);
      orientationValid.push_back(eventSourceSideValid);
      topology.push_back(measuredMuonTopology(geometryDiagnostics, eventSourceSideValid));
      bool const traversingCategory = candidate->source == 1;
      bool const dsaCategory = candidate->source == 0;
      bool const fullLeverArm = traversingCategory && std::abs(track.outerPosition().z() - track.innerPosition().z()) > 500.;
      quality.push_back(fullLeverArm ? 3 : (traversingCategory ? 2 : (dsaCategory ? 1 : 0)));
    }

    auto table = std::make_unique<nanoaod::FlatTable>(selected.size(), "ShiftMuon", false, false);
    table->addColumn<float>("pt", pt, "transverse momentum");
    table->addColumn<float>("eta", eta, "pseudorapidity");
    table->addColumn<float>("phi", phi, "azimuthal angle");
    table->addColumn<float>("mass", mass, "muon mass");
    table->addColumn<float>("p", p, "momentum magnitude");
    table->addColumn<float>("px", px, "momentum x component");
    table->addColumn<float>("py", py, "momentum y component");
    table->addColumn<float>("pz", pz, "momentum z component");
    // The historical unqualified branches are the unconstrained hypothesis.
    // Only the target-restricted alternative receives a prefix.
    table->addColumn<int>("constrainedValid", constrainedValid, "1 when the prompt-target Kalman fit is valid");
    table->addColumn<int>("constrainedHits", constrainedHits, "detector hits retained by the prompt-target fit");
    table->addColumn<int>("constrainedStatus",
                          constrainedStatus,
                          "prompt-target status: 1=valid, -10=no valid unconstrained refit, other negative=constraint failure");
    table->addColumn<float>("constrainedPt", constrainedPt, "prompt-target-constrained transverse momentum");
    table->addColumn<float>("constrainedEta", constrainedEta, "prompt-target-constrained pseudorapidity");
    table->addColumn<float>("constrainedPhi", constrainedPhi, "prompt-target-constrained azimuth");
    table->addColumn<float>("constrainedMass", constrainedMass, "muon mass for the constrained hypothesis");
    table->addColumn<float>("constrainedP", constrainedP, "prompt-target-constrained momentum magnitude");
    table->addColumn<float>("constrainedPx", constrainedPx, "prompt-target-constrained momentum x component");
    table->addColumn<float>("constrainedPy", constrainedPy, "prompt-target-constrained momentum y component");
    table->addColumn<float>("constrainedPz", constrainedPz, "prompt-target-constrained momentum z component");
    table->addColumn<float>("constrainedVx", constrainedVx, "fitted x on the configured target plane");
    table->addColumn<float>("constrainedVy", constrainedVy, "fitted y on the configured target plane");
    table->addColumn<float>("constrainedVz", constrainedVz, "configured target-plane z");
    table->addColumn<float>("constrainedChi2", constrainedChi2, "total prompt-target trajectory chi2");
    table->addColumn<float>("constrainedNdof", constrainedNdof, "prompt-target trajectory degrees of freedom");
    table->addColumn<float>("constrainedTargetChi2", constrainedTargetChi2, "chi2 contribution of the target hit");
    table->addColumn<float>("trackPt", trackPt, "pT at the original CMSSW track reference state");
    table->addColumn<int>("trackerMatchValid",
                          trackerMatchValid,
                          "1 when CMSSW cosmic-global matching found a generalTracks partner");
    table->addColumn<int>("trackerTrackIndex", trackerTrackIndex, "index in generalTracks, or -1");
    table->addColumn<int>("trackerValidHits", trackerValidHits, "valid tracker hits on the matched tracker track");
    table->addColumn<int>("trackerPixelHits", trackerPixelHits, "valid pixel hits on the matched tracker track");
    table->addColumn<int>("trackerStripHits", trackerStripHits, "valid strip hits on the matched tracker track");
    table->addColumn<float>("trackerPt", trackerPt, "matched tracker-track transverse momentum");
    table->addColumn<float>("trackerEta", trackerEta, "matched tracker-track pseudorapidity");
    table->addColumn<float>("trackerPhi", trackerPhi, "matched tracker-track azimuth");
    table->addColumn<float>("trackerP", trackerP, "matched tracker-track momentum magnitude");
    table->addColumn<float>("trackerPz", trackerPz, "matched tracker-track longitudinal momentum");
    table->addColumn<float>("trackerChi2", trackerChi2, "matched tracker-track chi2");
    table->addColumn<float>("trackerNdof", trackerNdof, "matched tracker-track degrees of freedom");
    table->addColumn<int>("combinedTrackValid",
                          combinedTrackValid,
                          "1 when the independent tracker+muon cosmic-global fit is available");
    table->addColumn<int>("combinedValidHits", combinedValidHits, "valid hits on the tracker+muon fit");
    table->addColumn<int>("combinedTrackerHits", combinedTrackerHits, "valid tracker hits on the combined fit");
    table->addColumn<int>("combinedMuonHits", combinedMuonHits, "valid muon-system hits on the combined fit");
    table->addColumn<float>("combinedPt", combinedPt, "tracker+muon fit transverse momentum at its reference state");
    table->addColumn<float>("combinedEta", combinedEta, "tracker+muon fit pseudorapidity at its reference state");
    table->addColumn<float>("combinedPhi", combinedPhi, "tracker+muon fit azimuth at its reference state");
    table->addColumn<float>("combinedP", combinedP, "tracker+muon fit momentum magnitude at its reference state");
    table->addColumn<float>("combinedPz", combinedPz, "tracker+muon fit longitudinal momentum at its reference state");
    table->addColumn<float>("combinedChi2", combinedChi2, "tracker+muon fit chi2");
    table->addColumn<float>("combinedNdof", combinedNdof, "tracker+muon fit degrees of freedom");
    table->addColumn<float>("combinedTargetPt",
                            combinedTargetPt,
                            "tracker+muon test hypothesis transverse momentum at the target line");
    table->addColumn<float>("combinedTargetPz",
                            combinedTargetPz,
                            "tracker+muon test hypothesis longitudinal momentum at the target line");
    table->addColumn<float>("combinedTargetDca",
                            combinedTargetDca,
                            "tracker+muon test hypothesis transverse DCA to the target line, or -1");
    table->addColumn<float>("innerPt", innerPt, "pT at the geometrically inner detector state");
    table->addColumn<float>("outerPt", outerPt, "pT at the geometrically outer detector state");
    table->addColumn<float>("upstreamPt", upstreamPt, "pT at the fitted endpoint nearest the inferred source side");
    table->addColumn<float>("preRefitPt", preRefitPt, "target-line pT before the timing-directed Kalman refit");
    table->addColumn<float>("preRefitPz", preRefitPz, "target-line pz before the timing-directed Kalman refit");
    table->addColumn<float>(
        "directionalRefitUpstreamPt", directionalRefitUpstreamPt, "pT at the source-facing refitted hit state");
    table->addColumn<float>(
        "directionalRefitUpstreamP", directionalRefitUpstreamP, "momentum at the source-facing refitted hit state");
    table->addColumn<float>("directionalRefitDownstreamP",
                            directionalRefitDownstreamP,
                            "momentum at the downstream refitted hit state");
    table->addColumn<float>("directionalRefitFractionalLossAcrossHits",
                            directionalRefitFractionalLossAcrossHits,
                            "fitted fractional momentum loss from the source-facing to downstream hit state");
    table->addColumn<float>("directionalRefitSourceFacingSurfaceDistance",
                            directionalRefitSourceFacingSurfaceDistance,
                            "distance in cm between the selected source-facing smoother state and first fit-hit surface");
    table->addColumn<float>("directionalRefitFilterDirectionCosine",
                            directionalRefitFilterDirectionCosine,
                            "cosine between source-to-downstream state displacement and source-facing momentum");
    table->addColumn<float>("directionalRefitTargetDirectionCosine",
                            directionalRefitTargetDirectionCosine,
                            "cosine of source-facing momentum along the target-side z direction; valid refits are negative");
    table->addColumn<int>("directionalRefitAttempted",
                          directionalRefitAttempted,
                          "1 when the timing/inferred-direction Kalman refit was attempted");
    table->addColumn<int>("directionalRefitValid",
                          directionalRefitValid,
                          "1 when final kinematics come from the timing/inferred-direction Kalman refit");
    table->addColumn<int>("directionalRefitHits", directionalRefitHits, "valid hits retained by the directional refit");
    table->addColumn<float>("directionalRefitChi2", directionalRefitChi2, "directional-refit trajectory chi2");
    table->addColumn<float>("directionalRefitNdof", directionalRefitNdof, "directional-refit trajectory ndof");
    table->addColumn<int>("directionalRefitFirstValid",
                          directionalRefitFirstValid,
                          "1 when the first all-hit control fitter-smoother iteration is valid");
    table->addColumn<int>("directionalRefitFirstHits",
                          directionalRefitFirstHits,
                          "valid hits in the first all-hit control fitter-smoother iteration");
    table->addColumn<float>("directionalRefitFirstChi2",
                            directionalRefitFirstChi2,
                            "trajectory chi2 from the first all-hit control fitter-smoother iteration");
    table->addColumn<float>("directionalRefitFirstNdof",
                            directionalRefitFirstNdof,
                            "trajectory ndof from the first all-hit control fitter-smoother iteration");
    table->addColumn<float>("directionalRefitFirstUpstreamPt",
                            directionalRefitFirstUpstreamPt,
                            "source-facing pT from the first all-hit control fitter-smoother iteration");
    table->addColumn<float>("directionalRefitFirstQoverP",
                            directionalRefitFirstQoverP,
                            "source-facing signed inverse momentum from all-hit control iteration one");
    table->addColumn<float>("directionalRefitFirstTargetPt",
                            directionalRefitFirstTargetPt,
                            "material-aware target-line pT from all-hit control iteration one");
    table->addColumn<float>("directionalRefitFirstTargetPz",
                            directionalRefitFirstTargetPz,
                            "material-aware target-line pz from all-hit control iteration one");
    table->addColumn<int>("directionalRefitSecondValid",
                          directionalRefitSecondValid,
                          "1 when the second all-hit control nonlinear iteration is valid");
    table->addColumn<int>("directionalRefitSecondHits",
                          directionalRefitSecondHits,
                          "valid hits in the second all-hit control nonlinear iteration");
    table->addColumn<float>("directionalRefitSecondChi2",
                            directionalRefitSecondChi2,
                            "trajectory chi2 from the second all-hit control nonlinear iteration");
    table->addColumn<float>("directionalRefitSecondNdof",
                            directionalRefitSecondNdof,
                            "trajectory ndof from the second all-hit control nonlinear iteration");
    table->addColumn<float>("directionalRefitSecondUpstreamPt",
                            directionalRefitSecondUpstreamPt,
                            "source-facing pT from the second all-hit control nonlinear iteration");
    table->addColumn<float>("directionalRefitSecondQoverP",
                            directionalRefitSecondQoverP,
                            "source-facing signed inverse momentum from all-hit control iteration two");
    table->addColumn<float>("directionalRefitSecondTargetPt",
                            directionalRefitSecondTargetPt,
                            "material-aware target-line pT from all-hit control iteration two");
    table->addColumn<float>("directionalRefitSecondTargetPz",
                            directionalRefitSecondTargetPz,
                            "material-aware target-line pz from all-hit control iteration two");
    table->addColumn<int>("directionalRefitSecondConverged",
                          directionalRefitSecondConverged,
                          "1 when all-hit control iteration two passes the relative q-over-p guard");
    table->addColumn<float>("directionalRefitRelativeQoverPChange",
                            directionalRefitRelativeQoverPChange,
                            "absolute all-hit control q-over-p change from iteration one to two, or -1");
    table->addColumn<int>("directionalRefitSelectedIteration",
                          directionalRefitSelectedIteration,
                          "selected all-hit control iteration: 0=none, 1=first, 2=second");
    table->addColumn<float>("directionalRefitVacuumTargetPt",
                            directionalRefitVacuumTargetPt,
                            "selected-iteration target-line pT with vacuum-only propagation");
    table->addColumn<float>("directionalRefitVacuumTargetPz",
                            directionalRefitVacuumTargetPz,
                            "selected-iteration target-line pz with vacuum-only propagation");
    table->addColumn<float>("directionalRefitMaterialDeltaPt",
                            directionalRefitMaterialDeltaPt,
                            "selected material-aware target pT minus vacuum-only target pT");
    table->addColumn<int>("directionalRefitMaterialBoundaryValid",
                          directionalRefitMaterialBoundaryValid,
                          "1 when detailed Geant4e propagation reached the outer CMS boundary");
    table->addColumn<float>("directionalRefitMaterialBoundaryPt",
                            directionalRefitMaterialBoundaryPt,
                            "pT at the outer CMS boundary after detailed Geant4e back-propagation");
    table->addColumn<float>("directionalRefitMaterialBoundaryPz",
                            directionalRefitMaterialBoundaryPz,
                            "pz at the outer CMS boundary after detailed Geant4e back-propagation");
    table->addColumn<float>("directionalRefitMaterialPath",
                            directionalRefitMaterialPath,
                            "signed Geant4e path length from the source-facing refit state to the CMS boundary");
    table->addColumn<int>("nDTRefitHits", nDTRefitHits, "valid DT inputs available to the directional refit");
    table->addColumn<int>("nCompatibleDTSegments", nCompatibleDTSegments, "unassigned DT segments geometrically compatible with this ShiftMuon");
    table->addColumn<int>("nCompatiblePixelHits", nCompatiblePixelHits, "pixel rechits geometrically compatible with this ShiftMuon");
    table->addColumn<int>("nCompatibleStripHits", nCompatibleStripHits, "matched strip rechits geometrically compatible with this ShiftMuon");
    table->addColumn<int>("nPropagatedDTSegments", nPropagatedDTSegments, "unassigned DT segments reached by propagation before compatibility cuts");
    table->addColumn<int>("nPropagatedPixelHits", nPropagatedPixelHits, "pixel rechits reached by propagation before compatibility cuts");
    table->addColumn<int>("nPropagatedStripHits", nPropagatedStripHits, "all matched and unmatched strip rechits reached by propagation before compatibility cuts");
    table->addColumn<float>("minDTResidual", minDTResidual, "minimum propagated-to-DT-segment global residual in cm; -1 when none propagated");
    table->addColumn<float>("minDTEstimatorChi2", minDTEstimatorChi2, "minimum covariance-aware DT compatibility chi2; -1 when none propagated");
    table->addColumn<float>("minTrackerResidual", minTrackerResidual, "minimum propagated-to-tracker-rechit global residual in cm; -1 when none propagated");
    table->addColumn<float>("minTrackerEstimatorChi2", minTrackerEstimatorChi2, "minimum covariance-aware tracker compatibility chi2; -1 when none propagated");
    table->addColumn<int>("nMatchedEcalRecHits", nMatchedEcalRecHits, "ECAL rechits associated to the fitted detector chord");
    table->addColumn<int>("nMatchedHBHERecHits", nMatchedHBHERecHits, "HBHE rechits associated to the fitted detector chord");
    table->addColumn<int>("nMatchedHFRecHits", nMatchedHFRecHits, "HF rechits associated to the fitted detector chord");
    table->addColumn<int>("nMatchedHORecHits", nMatchedHORecHits, "HO rechits associated to the fitted detector chord");
    table->addColumn<int>("nMatchedZDCRecHits", nMatchedZDCRecHits, "ZDC rechits associated to the fitted detector chord");
    table->addColumn<float>("matchedEcalEnergy", matchedEcalEnergy, "summed associated ECAL rechit energy");
    table->addColumn<float>("matchedHBHEEnergy", matchedHBHEEnergy, "summed associated HBHE rechit energy");
    table->addColumn<float>("matchedHFEnergy", matchedHFEnergy, "summed associated HF rechit energy");
    table->addColumn<float>("matchedHOEnergy", matchedHOEnergy, "summed associated HO rechit energy");
    table->addColumn<float>("matchedZDCEnergy", matchedZDCEnergy, "summed associated ZDC rechit energy");
    table->addColumn<float>("matchedEcalTime", matchedEcalTime, "energy-weighted associated ECAL rechit time");
    table->addColumn<float>("matchedHBHETime", matchedHBHETime, "energy-weighted associated HBHE rechit time");
    table->addColumn<float>("matchedHFTime", matchedHFTime, "energy-weighted associated HF rechit time");
    table->addColumn<float>("matchedHOTime", matchedHOTime, "energy-weighted associated HO rechit time");
    table->addColumn<float>("matchedZDCTime", matchedZDCTime, "energy-weighted associated ZDC rechit time in the shifted readout frame");
    table->addColumn<float>("minEcalLineDistance", minEcalLineDistance, "minimum ECAL-cell distance to the fitted detector chord in cm");
    table->addColumn<float>("minHBHELineDistance", minHBHELineDistance, "minimum HBHE-cell distance to the fitted detector chord in cm");
    table->addColumn<float>("minHFLineDistance", minHFLineDistance, "minimum HF-cell distance to the fitted detector chord in cm");
    table->addColumn<float>("minHOLineDistance", minHOLineDistance, "minimum HO-cell distance to the fitted detector chord in cm");
    table->addColumn<float>("minZDCLineDistance", minZDCLineDistance, "minimum ZDC-cell distance to the fitted detector chord in cm");
    table->addColumn<int>("caloTimingDirectionSign", caloTimingDirectionSign, "calorimeter-only time-order sign relative to fitted momentum; zero when inconclusive");
    table->addColumn<int>("nCaloTimingMeasurements", nCaloTimingMeasurements, "independent calorimeter subsystem timing centroids");
    table->addColumn<float>("caloTimingDeltaChi2", caloTimingDeltaChi2, "calorimeter time-order chi2 separation between directions");
    table->addColumn<int>("nAddedDTRefitHits", nAddedDTRefitHits, "compatible DT segments supplied to the augmented refit");
    table->addColumn<int>("nAddedTrackerRefitHits", nAddedTrackerRefitHits, "compatible tracker rechits supplied to the augmented refit");
    table->addColumn<int>("nCSCRefitHits", nCSCRefitHits, "valid CSC inputs available to the directional refit");
    table->addColumn<int>("nRPCRefitHits", nRPCRefitHits, "valid RPC inputs available to the directional refit");
    table->addColumn<int>("nGEMRefitHits", nGEMRefitHits, "valid GEM inputs available to the directional refit");
    table->addColumn<int>("nPrecisionRefitStations",
                          nPrecisionRefitStations,
                          "distinct DT, CSC, or GEM station layers available to the precision-only refit");
    table->addColumn<float>("precisionRefitLeverArm",
                            precisionRefitLeverArm,
                            "three-dimensional separation of the endpoint precision-hit positions in cm");
    table->addColumn<int>("directionalRefitUsedPrecisionHits",
                          directionalRefitUsedPrecisionHits,
                          "1 when canonical momentum uses the guarded precision-only refit");
    table->addColumn<int>("directionalRefitAllHitsValid",
                          directionalRefitAllHitsValid,
                          "1 when the all-muon-hit control refit is valid");
    table->addColumn<int>("directionalRefitAllHits",
                          directionalRefitAllHits,
                          "hits retained by the selected all-hit control iteration");
    table->addColumn<float>("directionalRefitAllHitsChi2",
                            directionalRefitAllHitsChi2,
                            "trajectory chi2 of the selected all-hit control iteration");
    table->addColumn<float>("directionalRefitAllHitsNdof",
                            directionalRefitAllHitsNdof,
                            "trajectory ndof of the selected all-hit control iteration");
    table->addColumn<float>("directionalRefitAllHitsTargetPt",
                            directionalRefitAllHitsTargetPt,
                            "target-line pT from the selected all-hit control iteration");
    table->addColumn<float>("directionalRefitAllHitsTargetPz",
                            directionalRefitAllHitsTargetPz,
                            "target-line pz from the selected all-hit control iteration");
    table->addColumn<int>("directionalRefitAllHitsSelectedIteration",
                          directionalRefitAllHitsSelectedIteration,
                          "selected all-hit control iteration: 0=none, 1=first, 2=second");
    table->addColumn<int>("directionalRefitAllHitsInput",
                          directionalRefitAllHitsInput,
                          "valid all-muon-system hits supplied to the refit");
    table->addColumn<int>("directionalRefitAllHitsOrderingFallback",
                          directionalRefitAllHitsOrderingFallback,
                          "all-hit measurements ordered by projection after path propagation failed");
    table->addColumn<float>("directionalRefitAllHitsOrderingSpan",
                            directionalRefitAllHitsOrderingSpan,
                            "ordered all-hit path span in cm");
    table->addColumn<int>("directionalRefitAllHitsFirstRejected",
                          directionalRefitAllHitsFirstRejected,
                          "all-hit inputs rejected by the loose first fit");
    table->addColumn<int>("directionalRefitAllHitsSecondRejected",
                          directionalRefitAllHitsSecondRejected,
                          "all-hit inputs rejected by the tight second fit");
    table->addColumn<int>("directionalRefitPrecisionValid",
                          directionalRefitPrecisionValid,
                          "1 when the DT, CSC, and GEM-only refit is valid");
    table->addColumn<int>("directionalRefitPrecisionHits",
                          directionalRefitPrecisionHits,
                          "precision hits retained by its selected nonlinear iteration");
    table->addColumn<float>("directionalRefitPrecisionChi2",
                            directionalRefitPrecisionChi2,
                            "trajectory chi2 of the selected precision-only iteration");
    table->addColumn<float>("directionalRefitPrecisionNdof",
                            directionalRefitPrecisionNdof,
                            "trajectory ndof of the selected precision-only iteration");
    table->addColumn<float>("directionalRefitPrecisionUpstreamPt",
                            directionalRefitPrecisionUpstreamPt,
                            "source-facing pT from the selected precision-only iteration");
    table->addColumn<float>("directionalRefitPrecisionQoverP",
                            directionalRefitPrecisionQoverP,
                            "source-facing signed inverse momentum from the precision-only refit");
    table->addColumn<float>("directionalRefitPrecisionTargetPt",
                            directionalRefitPrecisionTargetPt,
                            "target-line pT from the selected precision-only iteration");
    table->addColumn<float>("directionalRefitPrecisionTargetPz",
                            directionalRefitPrecisionTargetPz,
                            "target-line pz from the selected precision-only iteration");
    table->addColumn<int>("directionalRefitPrecisionSecondValid",
                          directionalRefitPrecisionSecondValid,
                          "1 when the second precision-only nonlinear iteration is valid");
    table->addColumn<int>("directionalRefitPrecisionSecondConverged",
                          directionalRefitPrecisionSecondConverged,
                          "1 when the second precision-only iteration passes the q-over-p guard");
    table->addColumn<float>("directionalRefitPrecisionRelativeQoverPChange",
                            directionalRefitPrecisionRelativeQoverPChange,
                            "absolute relative q-over-p change between precision-only iterations");
    table->addColumn<int>("directionalRefitPrecisionSelectedIteration",
                          directionalRefitPrecisionSelectedIteration,
                          "selected precision-only iteration: 0=none, 1=first, 2=second");
    table->addColumn<float>("directionalRefitPrecisionFirstTargetPt",
                            directionalRefitPrecisionFirstTargetPt,
                            "target-line pT from precision-only iteration one");
    table->addColumn<float>("directionalRefitPrecisionFirstTargetPz",
                            directionalRefitPrecisionFirstTargetPz,
                            "target-line pz from precision-only iteration one");
    table->addColumn<float>("directionalRefitPrecisionSecondTargetPt",
                            directionalRefitPrecisionSecondTargetPt,
                            "target-line pT from precision-only iteration two");
    table->addColumn<float>("directionalRefitPrecisionSecondTargetPz",
                            directionalRefitPrecisionSecondTargetPz,
                            "target-line pz from precision-only iteration two");
    table->addColumn<float>("directionalRefitPrecisionRelativeToAllQoverP",
                            directionalRefitPrecisionRelativeToAllQoverP,
                            "relative selected q-over-p difference between precision-only and all-hit refits");
    table->addColumn<float>("directionalRefitPrecisionTargetDca",
                            directionalRefitPrecisionTargetDca,
                            "precision-only selected-state transverse DCA to the target line in cm");
    table->addColumn<int>("directionalRefitSourceLegPrecisionValid",
                          directionalRefitSourceLegPrecisionValid,
                          "diagnostic source-side-only precision refit validity");
    table->addColumn<int>("directionalRefitSourceLegPrecisionHits",
                          directionalRefitSourceLegPrecisionHits,
                          "hits in the diagnostic source-side-only precision refit");
    table->addColumn<float>("directionalRefitSourceLegPrecisionQoverP",
                            directionalRefitSourceLegPrecisionQoverP,
                            "source-facing signed inverse momentum from the source-side-only precision refit");
    table->addColumn<float>("directionalRefitSourceLegPrecisionTargetPt",
                            directionalRefitSourceLegPrecisionTargetPt,
                            "target-line pT from the source-side-only precision refit");
    table->addColumn<float>("directionalRefitSourceLegPrecisionTargetPz",
                            directionalRefitSourceLegPrecisionTargetPz,
                            "target-line pz from the source-side-only precision refit");
    table->addColumn<int>("directionalRefitOppositeLegPrecisionValid",
                          directionalRefitOppositeLegPrecisionValid,
                          "diagnostic opposite-side-only precision refit validity");
    table->addColumn<int>("directionalRefitOppositeLegPrecisionHits",
                          directionalRefitOppositeLegPrecisionHits,
                          "hits in the diagnostic opposite-side-only precision refit");
    table->addColumn<float>("directionalRefitOppositeLegPrecisionQoverP",
                            directionalRefitOppositeLegPrecisionQoverP,
                            "source-facing signed inverse momentum from the opposite-side-only precision refit");
    table->addColumn<float>("directionalRefitOppositeLegPrecisionTargetPt",
                            directionalRefitOppositeLegPrecisionTargetPt,
                            "target-line pT from the opposite-side-only precision refit");
    table->addColumn<float>("directionalRefitOppositeLegPrecisionTargetPz",
                            directionalRefitOppositeLegPrecisionTargetPz,
                            "target-line pz from the opposite-side-only precision refit");
    table->addColumn<int>("directionalRefitPrecisionInput",
                          directionalRefitPrecisionInput,
                          "valid DT, CSC, and GEM hits supplied to the precision refit");
    table->addColumn<int>("directionalRefitPrecisionOrderingFallback",
                          directionalRefitPrecisionOrderingFallback,
                          "precision measurements ordered by projection after path propagation failed");
    table->addColumn<float>("directionalRefitPrecisionOrderingSpan",
                            directionalRefitPrecisionOrderingSpan,
                            "ordered precision-hit path span in cm");
    table->addColumn<int>("directionalRefitPrecisionFirstRejected",
                          directionalRefitPrecisionFirstRejected,
                          "precision inputs rejected by the loose first fit");
    table->addColumn<int>("directionalRefitPrecisionSecondRejected",
                          directionalRefitPrecisionSecondRejected,
                          "precision inputs rejected by the tight second fit");
    table->addColumn<float>("ptErr", ptError, "transverse-momentum uncertainty");
    table->addColumn<float>("etaErr", etaError, "pseudorapidity uncertainty");
    table->addColumn<float>("phiErr", phiError, "azimuthal-angle uncertainty");
    table->addColumn<int>("charge", charge, "electric charge");
    table->addColumn<int>("chargeMatchesGen",
                          chargeMatchesGen,
                          "MC diagnostic: 1/0 if charge agrees/disagrees with matched GenPart, -1 on data/unmatched");
    table->addColumn<int>("directionFlipped",
                          directionFlipped,
                          "1 when momentum and charge were reversed to point from inferred source toward CMS");
    table->addColumn<int>(
        "inferredSourceSide", inferredSourceSide, "sign of reconstructed linePcaZ: -1=-z source, +1=+z source");
    table->addColumn<int>("timingDirectionSign",
                          timingDirectionSign,
                          "chosen sign relative to fitted momentum; 0 when timing is inconclusive");
    table->addColumn<int>("nTimingMeasurements", nTimingMeasurements, "number of timed DT/CSC segments");
    table->addColumn<float>("timingChi2", timingChi2, "chi2 of the preferred time-of-flight direction");
    table->addColumn<float>(
        "timingDeltaChi2", timingDeltaChi2, "chi2 separation between the two time-of-flight directions");
    table->addColumn<float>("vx", vx, "x at field-aware stepping-helix PCA to the target line");
    table->addColumn<float>("vy", vy, "y at field-aware stepping-helix PCA to the target line");
    table->addColumn<float>("vz", vz, "z at field-aware stepping-helix PCA to the target line");
    table->addColumn<float>("trackVx", trackVx, "original CMSSW track reference-point x");
    table->addColumn<float>("trackVy", trackVy, "original CMSSW track reference-point y");
    table->addColumn<float>("trackVz", trackVz, "original CMSSW track reference-point z");
    table->addColumn<float>("dxy", dxy, "transverse impact parameter relative to the origin");
    table->addColumn<float>("dz", dz, "longitudinal impact parameter relative to the origin");
    table->addColumn<float>("innerR", innerR, "inner-state cylindrical radius");
    table->addColumn<float>("innerZ", innerZ, "inner-state z");
    table->addColumn<float>("outerR", outerR, "outer-state cylindrical radius");
    table->addColumn<float>("outerZ", outerZ, "outer-state z");
    table->addColumn<int>("nDTHits", nDTHits, "valid DT track measurements");
    table->addColumn<int>("nCSCHits", nCSCHits, "valid CSC track measurements");
    table->addColumn<int>("nRPCHits", nRPCHits, "valid RPC track measurements");
    table->addColumn<int>("nGEMHits", nGEMHits, "valid GEM track measurements");
    table->addColumn<int>("nME0Hits", nME0Hits, "valid ME0 track measurements");
    table->addColumn<int>("nHitsPlusEndcap", nHitsPlusEndcap, "valid track measurements in the +z muon endcap");
    table->addColumn<int>("nHitsBarrel", nHitsBarrel, "valid track measurements in the muon barrel");
    table->addColumn<int>("nHitsMinusEndcap", nHitsMinusEndcap, "valid track measurements in the -z muon endcap");
    table->addColumn<int>("nHitsNearEndcap",
                          nHitsNearEndcap,
                          "valid endcap measurements on the inferred source side");
    table->addColumn<int>("nHitsFarEndcap",
                          nHitsFarEndcap,
                          "valid endcap measurements opposite the inferred source side");
    table->addColumn<int>("nStationsNearEndcap",
                          nStationsNearEndcap,
                          "distinct detector-system stations with valid measurements in the source-side endcap");
    table->addColumn<int>("nStationsBarrel",
                          nStationsBarrel,
                          "distinct detector-system stations with valid measurements in the barrel");
    table->addColumn<int>("nStationsFarEndcap",
                          nStationsFarEndcap,
                          "distinct detector-system stations with valid measurements in the opposite endcap");
    table->addColumn<int>("detectorMask",
                          detectorMask,
                          "muon subdetector bit mask: bit0=DT, bit1=CSC, bit2=RPC, bit3=GEM, bit4=ME0");
    table->addColumn<int>("entryExitValid", entryExitValid, "1 when valid muon-hit entry and exit positions exist");
    table->addColumn<int>("entryRegion",
                          entryRegion,
                          "first recorded muon-hit region along inferred flight: -1=unknown, 0=near endcap, 1=barrel, 2=far endcap");
    table->addColumn<int>("exitRegion",
                          exitRegion,
                          "last recorded muon-hit region along inferred flight: -1=unknown, 0=near endcap, 1=barrel, 2=far endcap; not a stopping-point claim");
    table->addColumn<int>("entrySubdetector",
                          entrySubdetector,
                          "subdetector of first recorded muon hit: 1=DT, 2=CSC, 3=RPC, 4=GEM, 5=ME0");
    table->addColumn<int>("exitSubdetector",
                          exitSubdetector,
                          "subdetector of last recorded muon hit: 1=DT, 2=CSC, 3=RPC, 4=GEM, 5=ME0");
    table->addColumn<float>("entryX", entryX, "x of first recorded valid muon hit along inferred flight in cm");
    table->addColumn<float>("entryY", entryY, "y of first recorded valid muon hit along inferred flight in cm");
    table->addColumn<float>("entryZ", entryZ, "z of first recorded valid muon hit along inferred flight in cm");
    table->addColumn<float>("entryR", entryR, "radius of first recorded valid muon hit along inferred flight in cm");
    table->addColumn<float>("exitX", exitX, "x of last recorded valid muon hit along inferred flight in cm");
    table->addColumn<float>("exitY", exitY, "y of last recorded valid muon hit along inferred flight in cm");
    table->addColumn<float>("exitZ", exitZ, "z of last recorded valid muon hit along inferred flight in cm");
    table->addColumn<float>("exitR", exitR, "radius of last recorded valid muon hit along inferred flight in cm");
    table->addColumn<float>("hitSpan", hitSpan, "straight-line distance between first and last recorded muon hits in cm");
    table->addColumn<float>("hitDeltaZ", hitDeltaZ, "absolute z separation between first and last recorded muon hits in cm");
    table->addColumn<float>("linePcaR", linePcaR, "straight-line transverse PCA radius");
    table->addColumn<float>("linePcaZ", linePcaZ, "z at straight-line transverse PCA");
    table->addColumn<float>("chordLinePcaR", chordLinePcaR, "transverse PCA radius from the inner-to-outer chord");
    table->addColumn<float>("chordLinePcaZ", chordLinePcaZ, "z at transverse PCA from the inner-to-outer chord");
    table->addColumn<float>(
        "targetLinePath", targetLinePath, "signed stepping-helix path length from the inner state to target-line PCA");
    table->addColumn<float>("chi2", chi2, "track chi2");
    table->addColumn<float>("ndof", ndof, "track fit degrees of freedom");
    table->addColumn<float>("normalizedChi2", normalizedChi2, "track chi2 divided by ndof");
    table->addColumn<int>("nValidHits", validHits, "number of valid track hits");
    table->addColumn<int>("nValidMuonHits", validMuonHits, "number of valid muon-system hits");
    table->addColumn<int>("nMuonStations", muonStations, "muon stations with valid hits");
    table->addColumn<int>("nLostHits", lostHits, "number of lost track hits");
    table->addColumn<int>("topology",
                          topology,
                          "recorded-hit topology relative to inferred source side: 0=near-endcap-only, "
                          "1=near-endcap-and-barrel, 2=both-endcaps, 3=far-endcap-only, 4=unclassified");
    table->addColumn<int>("orientationValid",
                          orientationValid,
                          "1 when the event source side was inferred from at least one valid target-line state");
    table->addColumn<int>("recoAlgorithm",
                          recoAlgorithm,
                          "reconstruction provenance: 0=DSA, 1=strict traversing, 2=ordinary cosmic");
    table->addColumn<int>("quality",
                          quality,
                          "deprecated compatibility category: 0=cosmic, 1=DSA, 2=traversing, "
                          "3=traversing with abs(outerZ-innerZ)>500 cm");
    table->addColumn<int>("sourceIndex", sourceIndex, "index in the source collection identified by recoAlgorithm");
    table->addColumn<unsigned int>("duplicateGroupSize",
                                   duplicateGroupSize,
                                   "number of transitive input-track duplicates represented by this row");
    table->addColumn<unsigned int>("duplicateSourceMask",
                                   duplicateSourceMask,
                                   "algorithms present in the duplicate component: bit0=DSA, bit1=strict traversing, bit2=cosmic");
    table->addColumn<unsigned int>("duplicateDSACount",
                                   duplicateDSACount,
                                   "number of DSA candidates in the duplicate component");
    table->addColumn<unsigned int>("duplicateTraversingCount",
                                   duplicateTraversingCount,
                                   "number of strict-traversing candidates in the duplicate component");
    table->addColumn<unsigned int>("duplicateCosmicCount",
                                   duplicateCosmicCount,
                                   "number of ordinary-cosmic candidates in the duplicate component");
    table->addColumn<int>("genPartIdx", genPartIdx, "index in GenPart, or -1 when unmatched or on data");
    table->addColumn<float>(
        "genPartDeltaR", genPartDeltaR, "direction-ambiguous deltaR to matched GenPart, or -1 when unmatched/on data");
    table->addColumn<int>("simTruthMatched",
                          simTruthMatched,
                          "MC closure diagnostic: selected row matched to a primary SimTrack with precision SimHits");
    table->addColumn<int>("simTrackId", simTrackId, "matched Geant4 SimTrack id, or -1");
    table->addColumn<int>("simPixelHits", simPixelHits, "matched muon Geant4 hits in pixel sensitive volumes");
    table->addColumn<int>("simStripHits", simStripHits, "matched muon Geant4 hits in strip sensitive volumes");
    table->addColumn<float>("simTrackP", simTrackP, "matched SimTrack momentum at the production vertex");
    table->addColumn<int>("simDTHits", simDTHits, "matched muon Geant4 hits in DT sensitive volumes");
    table->addColumn<int>("simCSCHits", simCSCHits, "matched muon Geant4 hits in CSC sensitive volumes");
    table->addColumn<int>("simRPCHits", simRPCHits, "matched muon Geant4 hits in RPC sensitive volumes");
    table->addColumn<int>("simGEMHits", simGEMHits, "matched muon Geant4 hits in GEM sensitive volumes");
    table->addColumn<int>("nAddedDTTruthChamberMatches",
                          nAddedDTTruthChamberMatches,
                          "MC-only number of added DT segments in a chamber crossed by the matched signal muon");
    table->addColumn<int>("simHBHEHits", simHBHEHits, "matched muon Geant4 hits in HB or HE");
    table->addColumn<int>("simHFHits", simHFHits, "matched muon Geant4 hits in HF");
    table->addColumn<int>("simHOHits", simHOHits, "matched muon Geant4 hits in HO");
    table->addColumn<int>("simHcalHits",
                          simHcalHits,
                          "matched primary-muon HCAL SimHits; -1 unavailable, -2 diagnostic disabled");
    table->addColumn<float>("simHcalEnergy", simHcalEnergy, "summed matched primary-muon HCAL SimHit energy");
    table->addColumn<float>("simHcalFirstTime", simHcalFirstTime, "earliest matched primary-muon HCAL SimHit time");
    table->addColumn<float>("simHcalLastTime", simHcalLastTime, "latest matched primary-muon HCAL SimHit time");
    table->addColumn<int>("simZDCHits",
                          simZDCHits,
                          "matched primary-muon ZDC SimHits; -1 unavailable, -2 diagnostic disabled");
    table->addColumn<float>("simZDCEnergy", simZDCEnergy, "summed matched primary-muon ZDC SimHit energy");
    table->addColumn<float>("simZDCFirstTime", simZDCFirstTime, "earliest matched primary-muon ZDC SimHit time");
    table->addColumn<float>("simZDCLastTime", simZDCLastTime, "latest matched primary-muon ZDC SimHit time");
    table->addColumn<int>("simMuonDetectorMask",
                          simMuonDetectorMask,
                          "sensitive-volume crossing mask from SimHits: bit0=DT, bit1=CSC, bit2=RPC, bit3=GEM");
    table->addColumn<float>(
        "simFirstPrecisionHitP", simFirstPrecisionHitP, "true momentum at the first source-facing precision SimHit");
    table->addColumn<float>(
        "simLastPrecisionHitP", simLastPrecisionHitP, "true momentum at the last precision SimHit");
    table->addColumn<float>("simLossToFirstPrecisionHit",
                            simLossToFirstPrecisionHit,
                            "true fractional momentum loss from production to the first precision SimHit");
    table->addColumn<float>("simLossAcrossPrecisionHits",
                            simLossAcrossPrecisionHits,
                            "true fractional momentum loss between first and last precision SimHits");
    table->addColumn<float>(
        "simFirstPrecisionPath", simFirstPrecisionPath, "projected path from SimVertex to first precision SimHit in cm");
    event.put(std::move(table));

    // Fit every cleaned pair directly from its retained source tracks.  The
    // resulting indices always refer to ShiftMuon rows and therefore do not
    // depend on keeping any of the input collections in NanoAOD.
    std::vector<int> vertexMuonIdx1, vertexMuonIdx2, vertexIsOS, vertexTopologyMin, vertexTopologyMax;
    std::vector<int> vertexDcaStatus, vertexKalmanAttempted, vertexKalmanValid, vertexUsesLineFallback,
        vertexSameGenMuon, vertexGenIsOS, vertexConstrainedValid, vertexConstrainedDcaStatus,
        vertexConstrainedUsesLineFallback, vertexIsCosmicCosmic, vertexIsCosmicDSA, vertexIsCosmicTraversing,
        vertexIsCosmicDoubleTraversing, vertexIsDSADSA, vertexIsDSATraversing, vertexIsDSADoubleTraversing,
        vertexIsTraversingTraversing, vertexIsTraversingDoubleTraversing, vertexIsDoubleTraversingDoubleTraversing;
    std::vector<float> vertexVx, vertexVy, vertexVz, vertexVxError, vertexVyError, vertexVzError, vertexChi2,
        vertexNdof, vertexNormalizedChi2, vertexProbability, vertexMass, vertexPt, vertexEta, vertexPhi, vertexPz,
        vertexDca, vertexDcaX, vertexDcaY, vertexDcaZ, vertexOriginCompatibilityChi2,
        vertexOriginCompatibilityNormalizedChi2, vertexConstrainedVx, vertexConstrainedVy, vertexConstrainedVz,
        vertexConstrainedVxError, vertexConstrainedVyError, vertexConstrainedVzError, vertexConstrainedChi2,
        vertexConstrainedNdof, vertexConstrainedNormalizedChi2, vertexConstrainedProbability,
        vertexConstrainedMass, vertexConstrainedPt, vertexConstrainedEta, vertexConstrainedPhi, vertexConstrainedPz,
        vertexConstrainedDca, vertexConstrainedDcaX, vertexConstrainedDcaY, vertexConstrainedDcaZ,
        vertexConstrainedOriginCompatibilityChi2, vertexConstrainedOriginCompatibilityNormalizedChi2;
    constexpr double muonMass = 0.105658;

    struct PairChoice {
      unsigned int first;
      unsigned int second;
      StraightLineApproach approach;
      CommonLineVertex fit;
      double originChi2;
      double score;
    };
    std::vector<PairChoice> pairChoices;
    for (unsigned int first = 0; first < selected.size(); ++first) {
      for (unsigned int second = first + 1; second < selected.size(); ++second) {
        auto const& firstState = selected[first]->targetLineState;
        auto const& secondState = selected[second]->targetLineState;
        auto const lineApproach = straightLineApproach(firstState, secondState);
        if (!lineApproach.valid || lineApproach.distance > maxPairDca_)
          continue;

        int const firstCharge = candidateDirectionSign(*selected[first]) * selected[first]->track->charge();
        int const secondCharge = candidateDirectionSign(*selected[second]) * selected[second]->track->charge();
        if (requireOppositeSign_ && firstCharge == secondCharge)
          continue;

        auto const firstOrigin = firstState.position;
        auto const secondOrigin = secondState.position;
        // A physical SHIFT pair must point back to the same side.  This is a
        // reconstruction-only condition and prevents +z/-z averages near 0.
        if (firstOrigin.z() * secondOrigin.z() <= 0.)
          continue;
        double const deltaX = firstOrigin.x() - secondOrigin.x();
        double const deltaY = firstOrigin.y() - secondOrigin.y();
        double const deltaZ = firstOrigin.z() - secondOrigin.z();
        double const originChi2 =
            (deltaX * deltaX + deltaY * deltaY) / (2. * originTransverseResolution_ * originTransverseResolution_) +
            deltaZ * deltaZ / (2. * originZResolution_ * originZResolution_);
        constexpr double originNdof = 3.;
        if (originChi2 / originNdof > maxPairOriginNormalizedChi2_)
          continue;
        auto const fit = commonLineVertex(firstState,
                                          secondState,
                                          commonVertexLineResolution_,
                                          commonVertexBeamLineResolution_);
        if (!fit.valid)
          continue;
        double const score = originChi2 / originNdof +
                             std::pow(lineApproach.distance / commonVertexLineResolution_, 2) + fit.chi2 / fit.ndof;
        pairChoices.push_back({first, second, lineApproach, fit, originChi2, score});
      }
    }
    std::stable_sort(pairChoices.begin(), pairChoices.end(), [](PairChoice const& first, PairChoice const& second) {
      return first.score < second.score;
    });

    std::vector<bool> muonAlreadyUsed(selected.size(), false);
    unsigned int retainedVertices = 0;
    for (auto const& choice : pairChoices) {
      unsigned int const first = choice.first;
      unsigned int const second = choice.second;
      if (muonAlreadyUsed[first] || muonAlreadyUsed[second])
        continue;
      if (maxDimuonVertices_ > 0 && retainedVertices >= maxDimuonVertices_)
        break;
      muonAlreadyUsed[first] = true;
      muonAlreadyUsed[second] = true;
      ++retainedVertices;

      auto const& firstState = selected[first]->targetLineState;
      auto const& secondState = selected[second]->targetLineState;
      int const firstCharge = candidateDirectionSign(*selected[first]) * selected[first]->track->charge();
      int const secondCharge = candidateDirectionSign(*selected[second]) * selected[second]->track->charge();

      auto canonicalMomentum = [](Candidate const& candidate) {
        auto const& state = candidate.targetLineState;
        return std::array<double, 3>{state.momentum.x(), state.momentum.y(), state.momentum.z()};
      };
      auto const firstP = canonicalMomentum(*selected[first]);
      auto const secondP = canonicalMomentum(*selected[second]);
      double const pairPx = firstP[0] + secondP[0];
      double const pairPy = firstP[1] + secondP[1];
      double const pairPz = firstP[2] + secondP[2];
      double const firstEnergy = std::sqrt(firstState.momentum.mag2() + muonMass * muonMass);
      double const secondEnergy = std::sqrt(secondState.momentum.mag2() + muonMass * muonMass);
      double const mass2 =
          std::pow(firstEnergy + secondEnergy, 2) - pairPx * pairPx - pairPy * pairPy - pairPz * pairPz;
      double const pairPt = std::hypot(pairPx, pairPy);

      vertexMuonIdx1.push_back(first);
      vertexMuonIdx2.push_back(second);
      vertexIsOS.push_back(firstCharge != secondCharge);
      int sameGenMuon = -1;
      int genIsOS = -1;
      if (genParticles.isValid() && genPartIdx[first] >= 0 && genPartIdx[second] >= 0) {
        sameGenMuon = genPartIdx[first] == genPartIdx[second];
        int const firstGenPdgId = (*genParticles)[genPartIdx[first]].pdgId();
        int const secondGenPdgId = (*genParticles)[genPartIdx[second]].pdgId();
        genIsOS = firstGenPdgId * secondGenPdgId < 0;
      }
      vertexSameGenMuon.push_back(sameGenMuon);
      vertexGenIsOS.push_back(genIsOS);
      vertexKalmanAttempted.push_back(0);
      vertexKalmanValid.push_back(0);
      vertexUsesLineFallback.push_back(1);
      vertexVx.push_back(choice.fit.position.x());
      vertexVy.push_back(choice.fit.position.y());
      vertexVz.push_back(choice.fit.position.z());
      vertexVxError.push_back(choice.fit.error[0]);
      vertexVyError.push_back(choice.fit.error[1]);
      vertexVzError.push_back(choice.fit.error[2]);
      vertexChi2.push_back(choice.fit.chi2);
      vertexNdof.push_back(choice.fit.ndof);
      vertexNormalizedChi2.push_back(choice.fit.chi2 / choice.fit.ndof);
      vertexProbability.push_back(ChiSquaredProbability(choice.fit.chi2, choice.fit.ndof));
      vertexOriginCompatibilityChi2.push_back(choice.originChi2);
      vertexOriginCompatibilityNormalizedChi2.push_back(choice.originChi2 / 3.);
      vertexDcaStatus.push_back(choice.approach.valid);
      vertexDca.push_back(choice.approach.distance);
      vertexDcaX.push_back(choice.approach.midpoint.x());
      vertexDcaY.push_back(choice.approach.midpoint.y());
      vertexDcaZ.push_back(choice.approach.midpoint.z());
      vertexMass.push_back(std::sqrt(std::max(0., mass2)));
      vertexPt.push_back(pairPt);
      vertexPz.push_back(pairPz);
      vertexEta.push_back(pairPt > 0. ? std::asinh(pairPz / pairPt)
                                      : std::copysign(std::numeric_limits<float>::infinity(), pairPz));
      vertexPhi.push_back(std::atan2(pairPy, pairPx));

      // The pair identity is fixed by the unconstrained reconstruction above.
      // Build exactly one alternative from the two constrained muon states;
      // deliberately do not create constrained-unconstrained combinations.
      bool const bothConstrained = selected[first]->constrainedValid && selected[second]->constrainedValid;
      StraightLineApproach constrainedApproach{};
      CommonLineVertex constrainedFit{};
      if (bothConstrained) {
        constrainedApproach =
            straightLineApproach(selected[first]->constrainedState, selected[second]->constrainedState);
        constrainedFit = commonLineVertex(selected[first]->constrainedState,
                                          selected[second]->constrainedState,
                                          commonVertexLineResolution_,
                                          commonVertexBeamLineResolution_);
      }
      bool const constrainedValid = bothConstrained && constrainedApproach.valid && constrainedFit.valid;
      vertexConstrainedValid.push_back(constrainedValid);
      vertexConstrainedUsesLineFallback.push_back(constrainedValid);
      vertexConstrainedDcaStatus.push_back(constrainedApproach.valid);
      vertexConstrainedVx.push_back(constrainedValid ? constrainedFit.position.x() : 0.);
      vertexConstrainedVy.push_back(constrainedValid ? constrainedFit.position.y() : 0.);
      vertexConstrainedVz.push_back(constrainedValid ? constrainedFit.position.z() : 0.);
      vertexConstrainedVxError.push_back(constrainedValid ? constrainedFit.error[0] : 0.);
      vertexConstrainedVyError.push_back(constrainedValid ? constrainedFit.error[1] : 0.);
      vertexConstrainedVzError.push_back(constrainedValid ? constrainedFit.error[2] : 0.);
      vertexConstrainedChi2.push_back(constrainedValid ? constrainedFit.chi2 : 0.);
      vertexConstrainedNdof.push_back(constrainedValid ? constrainedFit.ndof : 0.);
      vertexConstrainedNormalizedChi2.push_back(constrainedValid ? constrainedFit.chi2 / constrainedFit.ndof : 0.);
      vertexConstrainedProbability.push_back(
          constrainedValid ? ChiSquaredProbability(constrainedFit.chi2, constrainedFit.ndof) : 0.);
      vertexConstrainedDca.push_back(constrainedApproach.valid ? constrainedApproach.distance : 0.);
      vertexConstrainedDcaX.push_back(constrainedApproach.valid ? constrainedApproach.midpoint.x() : 0.);
      vertexConstrainedDcaY.push_back(constrainedApproach.valid ? constrainedApproach.midpoint.y() : 0.);
      vertexConstrainedDcaZ.push_back(constrainedApproach.valid ? constrainedApproach.midpoint.z() : 0.);

      double constrainedOriginChi2 = 0.;
      if (bothConstrained) {
        auto const& constrainedFirstOrigin = selected[first]->constrainedState.position;
        auto const& constrainedSecondOrigin = selected[second]->constrainedState.position;
        double const constrainedDeltaX = constrainedFirstOrigin.x() - constrainedSecondOrigin.x();
        double const constrainedDeltaY = constrainedFirstOrigin.y() - constrainedSecondOrigin.y();
        double const constrainedDeltaZ = constrainedFirstOrigin.z() - constrainedSecondOrigin.z();
        constrainedOriginChi2 =
            (constrainedDeltaX * constrainedDeltaX + constrainedDeltaY * constrainedDeltaY) /
                (2. * originTransverseResolution_ * originTransverseResolution_) +
            constrainedDeltaZ * constrainedDeltaZ / (2. * originZResolution_ * originZResolution_);
      }
      vertexConstrainedOriginCompatibilityChi2.push_back(constrainedOriginChi2);
      vertexConstrainedOriginCompatibilityNormalizedChi2.push_back(constrainedOriginChi2 / 3.);

      if (bothConstrained) {
        auto const& constrainedFirstMomentum = selected[first]->constrainedState.momentum;
        auto const& constrainedSecondMomentum = selected[second]->constrainedState.momentum;
        double const constrainedPairPx = constrainedFirstMomentum.x() + constrainedSecondMomentum.x();
        double const constrainedPairPy = constrainedFirstMomentum.y() + constrainedSecondMomentum.y();
        double const constrainedPairPz = constrainedFirstMomentum.z() + constrainedSecondMomentum.z();
        double const constrainedPairPt = std::hypot(constrainedPairPx, constrainedPairPy);
        double const constrainedFirstEnergy = std::sqrt(constrainedFirstMomentum.mag2() + muonMass * muonMass);
        double const constrainedSecondEnergy = std::sqrt(constrainedSecondMomentum.mag2() + muonMass * muonMass);
        double const constrainedMass2 = std::pow(constrainedFirstEnergy + constrainedSecondEnergy, 2) -
                                        constrainedPairPx * constrainedPairPx -
                                        constrainedPairPy * constrainedPairPy -
                                        constrainedPairPz * constrainedPairPz;
        vertexConstrainedMass.push_back(std::sqrt(std::max(0., constrainedMass2)));
        vertexConstrainedPt.push_back(constrainedPairPt);
        vertexConstrainedPz.push_back(constrainedPairPz);
        vertexConstrainedEta.push_back(
            constrainedPairPt > 0. ? std::asinh(constrainedPairPz / constrainedPairPt)
                                   : std::copysign(std::numeric_limits<float>::infinity(), constrainedPairPz));
        vertexConstrainedPhi.push_back(std::atan2(constrainedPairPy, constrainedPairPx));
      } else {
        vertexConstrainedMass.push_back(0.);
        vertexConstrainedPt.push_back(0.);
        vertexConstrainedPz.push_back(0.);
        vertexConstrainedEta.push_back(0.);
        vertexConstrainedPhi.push_back(0.);
      }

      int const lowQuality = std::min(quality[first], quality[second]);
      int const highQuality = std::max(quality[first], quality[second]);
      vertexTopologyMin.push_back(std::min(topology[first], topology[second]));
      vertexTopologyMax.push_back(std::max(topology[first], topology[second]));
      vertexIsCosmicCosmic.push_back(lowQuality == 0 && highQuality == 0);
      vertexIsCosmicDSA.push_back(lowQuality == 0 && highQuality == 1);
      vertexIsCosmicTraversing.push_back(lowQuality == 0 && highQuality == 2);
      vertexIsCosmicDoubleTraversing.push_back(lowQuality == 0 && highQuality == 3);
      vertexIsDSADSA.push_back(lowQuality == 1 && highQuality == 1);
      vertexIsDSATraversing.push_back(lowQuality == 1 && highQuality == 2);
      vertexIsDSADoubleTraversing.push_back(lowQuality == 1 && highQuality == 3);
      vertexIsTraversingTraversing.push_back(lowQuality == 2 && highQuality == 2);
      vertexIsTraversingDoubleTraversing.push_back(lowQuality == 2 && highQuality == 3);
      vertexIsDoubleTraversingDoubleTraversing.push_back(lowQuality == 3 && highQuality == 3);
    }

    auto vertexTable = std::make_unique<nanoaod::FlatTable>(vertexMuonIdx1.size(), "ShiftDimuonVertex", false, false);
    vertexTable->addColumn<int>("muonIdx1", vertexMuonIdx1, "index of first muon in ShiftMuon");
    vertexTable->addColumn<int>("muonIdx2", vertexMuonIdx2, "index of second muon in ShiftMuon");
    vertexTable->addColumn<int>("isOS", vertexIsOS, "1 for an opposite-sign pair");
    vertexTable->addColumn<int>("topologyMin",
                                vertexTopologyMin,
                                "lower of the two unordered ShiftMuon topology values");
    vertexTable->addColumn<int>("topologyMax",
                                vertexTopologyMax,
                                "higher of the two unordered ShiftMuon topology values");
    vertexTable->addColumn<int>("isCosmicCosmic", vertexIsCosmicCosmic, "both muons have quality 0 (cosmic)");
    vertexTable->addColumn<int>("isCosmicDSA", vertexIsCosmicDSA, "one cosmic and one DSA muon");
    vertexTable->addColumn<int>("isCosmicTraversing", vertexIsCosmicTraversing, "one cosmic and one traversing muon");
    vertexTable->addColumn<int>("isCosmicDoubleTraversing",
                                vertexIsCosmicDoubleTraversing,
                                "one cosmic and one full-lever-arm traversing muon");
    vertexTable->addColumn<int>("isDSADSA", vertexIsDSADSA, "both muons have quality 1 (DSA)");
    vertexTable->addColumn<int>("isDSATraversing", vertexIsDSATraversing, "one DSA and one traversing muon");
    vertexTable->addColumn<int>("isDSADoubleTraversing",
                                vertexIsDSADoubleTraversing,
                                "one DSA and one full-lever-arm traversing muon");
    vertexTable->addColumn<int>("isTraversingTraversing",
                                vertexIsTraversingTraversing,
                                "both muons have quality 2 (traversing)");
    vertexTable->addColumn<int>("isTraversingDoubleTraversing",
                                vertexIsTraversingDoubleTraversing,
                                "one traversing and one full-lever-arm traversing muon");
    vertexTable->addColumn<int>("isDoubleTraversingDoubleTraversing",
                                vertexIsDoubleTraversingDoubleTraversing,
                                "both muons have quality 3 (full-lever-arm traversing)");
    vertexTable->addColumn<int>(
        "kalmanAttempted", vertexKalmanAttempted, "1 when the pair lies inside safe transient-track propagation range");
    vertexTable->addColumn<int>("sameGenMuon",
                                vertexSameGenMuon,
                                "MC diagnostic: 1 when both legs match the same GenPart, -1 on data/unmatched");
    vertexTable->addColumn<int>(
        "genIsOS", vertexGenIsOS, "MC diagnostic: generator pair is opposite-sign, -1 on data/unmatched");
    vertexTable->addColumn<int>("kalmanValid",
                                vertexKalmanValid,
                                "1 when the unconstrained far-vertex Kalman fit converged near its line seed");
    vertexTable->addColumn<int>("usesLineFallback",
                                vertexUsesLineFallback,
                                "1 when position comes from straight-line closest approach after Kalman failure");
    vertexTable->addColumn<float>("vx", vertexVx, "unbounded common-line fit vertex x");
    vertexTable->addColumn<float>("vy", vertexVy, "unbounded common-line fit vertex y");
    vertexTable->addColumn<float>("vz", vertexVz, "unbounded common-line fit vertex z");
    vertexTable->addColumn<float>("vxErr", vertexVxError, "common-line fit vertex x uncertainty");
    vertexTable->addColumn<float>("vyErr", vertexVyError, "common-line fit vertex y uncertainty");
    vertexTable->addColumn<float>("vzErr", vertexVzError, "common-line fit vertex z uncertainty");
    vertexTable->addColumn<float>("chi2", vertexChi2, "common-line vertex fit chi2");
    vertexTable->addColumn<float>("ndof", vertexNdof, "common-line vertex fit degrees of freedom");
    vertexTable->addColumn<float>("normalizedChi2", vertexNormalizedChi2, "common-line vertex chi2 divided by ndof");
    vertexTable->addColumn<float>("probability", vertexProbability, "common-line vertex fit probability");
    vertexTable->addColumn<float>("originCompatibilityChi2",
                                  vertexOriginCompatibilityChi2,
                                  "compatibility chi2 of the two independent ShiftMuon origins");
    vertexTable->addColumn<float>("originCompatibilityNormalizedChi2",
                                  vertexOriginCompatibilityNormalizedChi2,
                                  "independent-origin compatibility chi2 divided by three");
    vertexTable->addColumn<int>("dcaValid", vertexDcaStatus, "1 when the two-track DCA calculation succeeded");
    vertexTable->addColumn<float>("dca", vertexDca, "three-dimensional distance of closest approach");
    vertexTable->addColumn<float>("dcaX", vertexDcaX, "x of the two-track closest-approach crossing point");
    vertexTable->addColumn<float>("dcaY", vertexDcaY, "y of the two-track closest-approach crossing point");
    vertexTable->addColumn<float>("dcaZ", vertexDcaZ, "z of the two-track closest-approach crossing point");
    vertexTable->addColumn<float>("mass", vertexMass, "dimuon invariant mass using canonical SHIFT directions");
    vertexTable->addColumn<float>("pt", vertexPt, "dimuon transverse momentum");
    vertexTable->addColumn<float>("pz", vertexPz, "dimuon longitudinal momentum");
    vertexTable->addColumn<float>("eta", vertexEta, "dimuon pseudorapidity");
    vertexTable->addColumn<float>("phi", vertexPhi, "dimuon azimuthal angle");
    vertexTable->addColumn<int>("constrainedValid",
                                vertexConstrainedValid,
                                "1 when both constrained muon fits and the constrained common-line fit are valid");
    vertexTable->addColumn<int>("constrainedUsesLineFallback",
                                vertexConstrainedUsesLineFallback,
                                "1 when constrained position comes from the common-line fit");
    vertexTable->addColumn<float>("constrainedVx", vertexConstrainedVx, "constrained common-line vertex x");
    vertexTable->addColumn<float>("constrainedVy", vertexConstrainedVy, "constrained common-line vertex y");
    vertexTable->addColumn<float>("constrainedVz", vertexConstrainedVz, "constrained common-line vertex z");
    vertexTable->addColumn<float>("constrainedVxErr", vertexConstrainedVxError, "constrained vertex x uncertainty");
    vertexTable->addColumn<float>("constrainedVyErr", vertexConstrainedVyError, "constrained vertex y uncertainty");
    vertexTable->addColumn<float>("constrainedVzErr", vertexConstrainedVzError, "constrained vertex z uncertainty");
    vertexTable->addColumn<float>("constrainedChi2", vertexConstrainedChi2, "constrained common-line fit chi2");
    vertexTable->addColumn<float>("constrainedNdof", vertexConstrainedNdof, "constrained common-line fit degrees of freedom");
    vertexTable->addColumn<float>("constrainedNormalizedChi2",
                                  vertexConstrainedNormalizedChi2,
                                  "constrained common-line chi2 divided by ndof");
    vertexTable->addColumn<float>("constrainedProbability",
                                  vertexConstrainedProbability,
                                  "constrained common-line fit probability");
    vertexTable->addColumn<float>("constrainedOriginCompatibilityChi2",
                                  vertexConstrainedOriginCompatibilityChi2,
                                  "compatibility chi2 of the two constrained ShiftMuon origins");
    vertexTable->addColumn<float>("constrainedOriginCompatibilityNormalizedChi2",
                                  vertexConstrainedOriginCompatibilityNormalizedChi2,
                                  "constrained-origin compatibility chi2 divided by three");
    vertexTable->addColumn<int>("constrainedDcaValid",
                                vertexConstrainedDcaStatus,
                                "1 when the constrained two-track DCA calculation succeeded");
    vertexTable->addColumn<float>("constrainedDca", vertexConstrainedDca, "constrained three-dimensional DCA");
    vertexTable->addColumn<float>("constrainedDcaX", vertexConstrainedDcaX, "constrained DCA midpoint x");
    vertexTable->addColumn<float>("constrainedDcaY", vertexConstrainedDcaY, "constrained DCA midpoint y");
    vertexTable->addColumn<float>("constrainedDcaZ", vertexConstrainedDcaZ, "constrained DCA midpoint z");
    vertexTable->addColumn<float>("constrainedMass", vertexConstrainedMass, "dimuon mass from two constrained muons");
    vertexTable->addColumn<float>("constrainedPt", vertexConstrainedPt, "dimuon pT from two constrained muons");
    vertexTable->addColumn<float>("constrainedPz", vertexConstrainedPz, "dimuon pz from two constrained muons");
    vertexTable->addColumn<float>("constrainedEta", vertexConstrainedEta, "dimuon eta from two constrained muons");
    vertexTable->addColumn<float>("constrainedPhi", vertexConstrainedPhi, "dimuon phi from two constrained muons");
    event.put(std::move(vertexTable), "ShiftDimuonVertex");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription description;
    description.add<edm::InputTag>("dsaTracks", edm::InputTag("displacedStandAloneMuons"));
    description.add<edm::InputTag>("cosmicTracks", edm::InputTag("shiftCosmicMuons"));
    description.add<edm::InputTag>("traversingTracks", edm::InputTag("shiftTraversingMuons"));
    description.add<edm::InputTag>("dsaGlobalLinks", edm::InputTag("shiftGlobalDSAMuons"));
    description.add<edm::InputTag>("cosmicGlobalLinks", edm::InputTag("shiftGlobalCosmicMuons"));
    description.add<edm::InputTag>("traversingGlobalLinks", edm::InputTag("shiftGlobalTraversingMuons"));
    description.add<edm::InputTag>("genParticles", edm::InputTag("finalGenParticles"));
    description.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
    description.add<edm::InputTag>("simVertices", edm::InputTag("g4SimHits"));
    description.add<edm::InputTag>("dtSimHits", edm::InputTag("g4SimHits", "MuonDTHits"));
    description.add<edm::InputTag>("cscSimHits", edm::InputTag("g4SimHits", "MuonCSCHits"));
    description.add<edm::InputTag>("rpcSimHits", edm::InputTag("g4SimHits", "MuonRPCHits"));
    description.add<edm::InputTag>("gemSimHits", edm::InputTag("g4SimHits", "MuonGEMHits"));
    description.add<std::vector<edm::InputTag>>("trackerSimHits", {});
    description.add<edm::InputTag>("hcalSimHits", edm::InputTag("g4SimHits", "HcalHits"));
    description.add<edm::InputTag>("zdcSimHits", edm::InputTag("g4SimHits", "ZDCHITS"));
    description.add<edm::InputTag>("ecalBarrelRecHits", edm::InputTag("reducedEcalRecHitsEB"));
    description.add<edm::InputTag>("ecalEndcapRecHits", edm::InputTag("reducedEcalRecHitsEE"));
    description.add<edm::InputTag>("hbheRecHits", edm::InputTag("hbhereco"));
    description.add<edm::InputTag>("hfRecHits", edm::InputTag("hfreco"));
    description.add<edm::InputTag>("hoRecHits", edm::InputTag("horeco"));
    description.add<edm::InputTag>("zdcRecHits", edm::InputTag("shiftZDCReco"));
    description.add<edm::InputTag>("dtSegments", edm::InputTag("dt4DSegments"));
    description.add<edm::InputTag>("pixelRecHits", edm::InputTag("siPixelRecHits"));
    description.add<edm::InputTag>("stripMatchedRecHits", edm::InputTag("siStripMatchedRecHits", "matchedRecHit"));
    description.add<edm::InputTag>("stripRphiRecHits", edm::InputTag("siStripMatchedRecHits", "rphiRecHit"));
    description.add<edm::InputTag>("stripRphiUnmatchedRecHits",
                                   edm::InputTag("siStripMatchedRecHits", "rphiRecHitUnmatched"));
    description.add<edm::InputTag>("stripStereoRecHits", edm::InputTag("siStripMatchedRecHits", "stereoRecHit"));
    description.add<edm::InputTag>("stripStereoUnmatchedRecHits",
                                   edm::InputTag("siStripMatchedRecHits", "stereoRecHitUnmatched"));
    description.add<std::string>("muonRecHitBuilder", "MuonRecHitBuilder");
    description.add<std::string>("trackerRecHitBuilder", "WithTrackAngle");
    description.add<bool>("augmentDTHits", false);
    description.add<bool>("augmentTrackerHits", false);
    description.add<bool>("useExtendedTiming", false);
    description.add<double>("additionalDTMaxDistance", 200.0);
    description.add<double>("additionalTrackerMaxDistance", 30.0);
    description.add<double>("additionalHitMaxChi2", 100.0);
    description.add<double>("ecalMatchMaxDistance", 15.0);
    description.add<double>("hcalMatchMaxDistance", 40.0);
    description.add<double>("zdcMatchMaxDistance", 200.0);
    description.add<double>("caloMatchMinEnergy", 0.01);
    description.add<bool>("useImprovedMomentumRefit", false);
    description.add<bool>("useDetailedMaterialPropagation", false);
    description.add<bool>("directionalRefitUseMaterialEffects", false);
    description.add<bool>("directionalRefitUseFirstPrinciplesMaterialEffects", false);
    description.add<double>("directionalRefitFirstPrinciplesStepCm", 0.2);
    description.add<bool>("directionalRefitUseDetailedMaterialEffects", false);
    description.add<bool>("directionalRefitUseGeometryMaterialEffects", false);
    description.add<bool>("directionalRefitUseGeometryMaterialEffectsInFitter", false);
    description.add<bool>("directionalRefitUseGeometryMaterialEffectsInSmoother", false);
    description.add<bool>("directionalRefitUseGeometryTargetMaterialEffects", false);
    description.add<double>("directionalRefitGeometryMaterialStepCm", 1.0);
    description.add<bool>("directionalRefitLogGeometryMaterialComparison", false);
    description.add<bool>("usePropagatedPathOrdering", false);
    description.add<bool>("directionalRefitUseFullSeedErrorRescale", false);
    description.add<double>("directionalRefitSeedCurvatureErrorRescale", 100.0);
    description.add<double>("directionalRefitSecondSeedErrorRescale", 100.0);
    description.add<double>("directionalRefitSeedMomentumScale", 1.0);
    description.add<double>("directionalRefitEnergyLossScale", 1.0);
    description.add<double>("directionalRefitErrorRescale", 10.0);
    description.add<double>("directionalRefitInitialMaxHitChi2", 100000.0);
    description.add<double>("directionalRefitMaxHitChi2", 100000.0);
    description.add<double>("directionalRefitMaxRelativeQoverPChange", 0.5);
    description.add<bool>("directionalRefitUseSecondIteration", false);
    description.add<unsigned int>("directionalRefitMinPrecisionStations", 2);
    description.add<double>("directionalRefitMaxPrecisionRelativeQoverPChange", 0.5);
    description.add<bool>("produceMomentumClosureDiagnostics", false);
    description.add<bool>("enableHcalDiagnostics", false);
    description.add<bool>("enableZDCDiagnostics", false);
    description.add<bool>("produceSplitLegRefits", false);
    description.add<bool>("directionalRefitUseExplicitBackwardTargetPropagation", false);
    description.add<bool>("produceTargetConstrainedMomentum", true);
    description.add<bool>("targetUseInferredSide", true);
    description.add<double>("targetX", 0.0);
    description.add<double>("targetY", 0.0);
    description.add<double>("targetZ", 14800.0);
    description.add<double>("targetSigmaX", 0.1);
    description.add<double>("targetSigmaY", 0.1);
    description.add<double>("targetSigmaZ", 50.0);
    description.add<double>("minSharedHitFraction", 0.5);
    description.add<unsigned int>("minSharedDetIds", 2);
    description.add<double>("maxDuplicateAngle", 0.03);
    description.add<double>("maxDuplicateLineDistance", 30.0);
    description.add<unsigned int>("minTimingMeasurements", 2);
    description.add<double>("minTimingDeltaChi2", 4.0);
    description.add<double>("maxTargetLineDca", 200.0);
    description.add<double>("minAbsEta", 3.0);
    description.add<double>("originTransverseResolution", 100.0);
    description.add<double>("originZResolution", 2000.0);
    description.add<double>("commonVertexLineResolution", 100.0);
    description.add<double>("commonVertexBeamLineResolution", 100.0);
    description.add<double>("maxPairOriginNormalizedChi2", 9.0);
    description.add<double>("maxPairDca", 500.0);
    description.add<unsigned int>("maxDimuonVertices", 1);
    description.add<bool>("requireOppositeSign", true);
    description.add<double>("maxGenDeltaR", 0.5);
    descriptions.add("shiftMuonTable", description);
  }

private:
  edm::EDGetTokenT<reco::TrackCollection> dsaToken_;
  edm::EDGetTokenT<reco::TrackCollection> cosmicToken_;
  edm::EDGetTokenT<reco::TrackCollection> traversingToken_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> dsaGlobalLinksToken_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> cosmicGlobalLinksToken_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> traversingGlobalLinksToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTracksToken_;
  edm::EDGetTokenT<edm::SimVertexContainer> simVerticesToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> dtSimHitsToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> cscSimHitsToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> rpcSimHitsToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> gemSimHitsToken_;
  std::vector<edm::EDGetTokenT<edm::PSimHitContainer>> trackerSimHitsTokens_;
  edm::EDGetTokenT<edm::PCaloHitContainer> hcalSimHitsToken_;
  edm::EDGetTokenT<edm::PCaloHitContainer> zdcSimHitsToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalBarrelRecHitsToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalEndcapRecHitsToken_;
  edm::EDGetTokenT<HBHERecHitCollection> hbheRecHitsToken_;
  edm::EDGetTokenT<HFRecHitCollection> hfRecHitsToken_;
  edm::EDGetTokenT<HORecHitCollection> hoRecHitsToken_;
  edm::EDGetTokenT<ZDCRecHitCollection> zdcRecHitsToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsToken_;
  edm::EDGetTokenT<SiPixelRecHitCollection> pixelRecHitsToken_;
  edm::EDGetTokenT<SiStripMatchedRecHit2DCollection> stripMatchedRecHitsToken_;
  edm::EDGetTokenT<SiStripRecHit2DCollection> stripRphiRecHitsToken_;
  edm::EDGetTokenT<SiStripRecHit2DCollection> stripRphiUnmatchedRecHitsToken_;
  edm::EDGetTokenT<SiStripRecHit2DCollection> stripStereoRecHitsToken_;
  edm::EDGetTokenT<SiStripRecHit2DCollection> stripStereoUnmatchedRecHitsToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> trackingGeometryToken_;
  edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> muonRecHitBuilderToken_;
  edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> trackerRecHitBuilderToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  bool augmentDTHits_;
  bool augmentTrackerHits_;
  bool useExtendedTiming_;
  double additionalDTMaxDistance_;
  double additionalTrackerMaxDistance_;
  double additionalHitMaxChi2_;
  double ecalMatchMaxDistance_;
  double hcalMatchMaxDistance_;
  double zdcMatchMaxDistance_;
  double caloMatchMinEnergy_;
  bool useImprovedMomentumRefit_;
  bool useDetailedMaterialPropagation_;
  bool directionalRefitUseMaterialEffects_;
  bool directionalRefitUseFirstPrinciplesMaterialEffects_;
  double directionalRefitFirstPrinciplesStepCm_;
  bool directionalRefitUseDetailedMaterialEffects_;
  bool directionalRefitUseGeometryMaterialEffects_;
  bool directionalRefitUseGeometryMaterialEffectsInFitter_;
  bool directionalRefitUseGeometryMaterialEffectsInSmoother_;
  bool directionalRefitUseGeometryTargetMaterialEffects_;
  double directionalRefitGeometryMaterialStepCm_;
  bool directionalRefitLogGeometryMaterialComparison_;
  bool usePropagatedPathOrdering_;
  bool directionalRefitUseFullSeedErrorRescale_;
  double directionalRefitSeedCurvatureErrorRescale_;
  double directionalRefitSecondSeedErrorRescale_;
  double directionalRefitSeedMomentumScale_;
  double directionalRefitEnergyLossScale_;
  double directionalRefitErrorRescale_;
  double directionalRefitInitialMaxHitChi2_;
  double directionalRefitMaxHitChi2_;
  double directionalRefitMaxRelativeQoverPChange_;
  bool directionalRefitUseSecondIteration_;
  unsigned int directionalRefitMinPrecisionStations_;
  double directionalRefitMaxPrecisionRelativeQoverPChange_;
  bool produceMomentumClosureDiagnostics_;
  bool enableHcalDiagnostics_;
  bool enableZDCDiagnostics_;
  bool produceSplitLegRefits_;
  bool directionalRefitUseExplicitBackwardTargetPropagation_;
  bool produceTargetConstrainedMomentum_;
  bool targetUseInferredSide_;
  double targetX_;
  double targetY_;
  double targetZ_;
  double targetSigmaX_;
  double targetSigmaY_;
  double targetSigmaZ_;
  double minSharedHitFraction_;
  unsigned int minSharedDetIds_;
  double maxDuplicateAngle_;
  double maxDuplicateLineDistance_;
  unsigned int minTimingMeasurements_;
  double minTimingDeltaChi2_;
  double maxTargetLineDca_;
  double minAbsEta_;
  double originTransverseResolution_;
  double originZResolution_;
  double commonVertexLineResolution_;
  double commonVertexBeamLineResolution_;
  double maxPairOriginNormalizedChi2_;
  double maxPairDca_;
  unsigned int maxDimuonVertices_;
  bool requireOppositeSign_;
  double maxGenDeltaR_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftMuonTableProducer);
