#ifndef TrackPropagation_Geant4ePropagator_h
#define TrackPropagation_Geant4ePropagator_h

#include <memory>
#include <limits>
#include <string>

// CMS includes
// - Propagator
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

// - Geant4e
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorPropagatorManager.hh"
#include "G4ErrorSurfaceTarget.hh"

/** Propagator based on the Geant4e package. Uses the Propagator class
 *  in the TrackingTools/GeomPropagators package to define the interface.
 *  See that class for more details.
 */

class Geant4ePropagator : public Propagator {
public:
  /** Constructor. Takes as arguments:
   *  * The magnetic field
   *  * The particle name whose properties will be used in the propagation.
   * Without the charge, i.e. "mu", "pi", ...
   *  * The propagation direction. It may be: alongMomentum, oppositeToMomentum
   */
  Geant4ePropagator(const MagneticField *field = nullptr,
                    std::string particleName = "mu",
                    PropagationDirection dir = alongMomentum,
                    double plimit = 1.0,
                    double maximumStepLengthMm = 10.0,
                    double maximumPathLengthCm = 200.0);

  ~Geant4ePropagator() override;

  /** The methods propagateWithPath() are identical to the corresponding
   *  methods propagate() in what concerns the resulting
   *  TrajectoryStateOnSurface, but they provide in addition the
   *  exact path length along the trajectory.
   *  All of these method calls are internally mapped to
   */

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const FreeTrajectoryState &,
                                                                const Plane &) const override;

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const FreeTrajectoryState &,
                                                                const Cylinder &) const override;

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const TrajectoryStateOnSurface &,
                                                                const Plane &) const override;

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const TrajectoryStateOnSurface &,
                                                                const Cylinder &) const override;

  Geant4ePropagator *clone() const override { return new Geant4ePropagator(*this); }

  const MagneticField *magneticField() const override { return theField; }

  // Opt-in physical backward covariance transport. The momentum reversed for
  // Geant4 must carry its own curvilinear frame; propagation follows that
  // frame with a positive step and the returned covariance is reversed back.
  void setUseConsistentBackwardCovariance(bool value) { consistentBackwardCovariance_ = value; }
  void setUseMeanEnergyLossJacobian(bool value) { meanEnergyLossJacobian_ = value; }
  void setUseFieldGradientJacobian(bool value) { fieldGradientJacobian_ = value; }
  void setUseUnquenchedIonizationVariance(bool value) { unquenchedIonizationVariance_ = value; }
  void setUseMaterialInterfaceJacobian(bool value) { materialInterfaceJacobian_ = value; }
  double maximumInterfaceRelativeCurvatureSigma() const { return maximumInterfaceRelativeCurvatureSigma_; }
  double maximumInterfaceZCm() const { return maximumInterfaceZCm_; }
  void setRecordTransportJacobian(bool value) { recordTransportJacobian_ = value; }
  bool transportJacobianValid() const { return transportJacobianValid_; }
  AlgebraicMatrix55 const& curvilinearTransportJacobian() const { return curvilinearTransportJacobian_; }
  // Diagnostic endpoint flows only: d(q/p)/ds along PHYSICAL momentum,
  // in GeV^-1 cm^-1. No change to the propagated mean, covariance or J.
  // Source uses the first pre-step material; destination uses the final
  // approached (pre-step) material. A geometric endpoint boundary is
  // conservatively ambiguous even if both named materials happen to agree.
  bool meanCurvatureFlowValid() const { return meanCurvatureFlowValid_; }
  double startMeanCurvatureFlow() const { return startMeanCurvatureFlow_; }
  double endMeanCurvatureFlow() const { return endMeanCurvatureFlow_; }
  bool materialEndpointBoundaryAmbiguous() const { return materialEndpointBoundaryAmbiguous_; }
  std::string const& startMeanCurvatureMaterial() const { return startMeanCurvatureMaterial_; }
  std::string const& endMeanCurvatureMaterial() const { return endMeanCurvatureMaterial_; }
  void setTransportAudit(bool value) { transportAudit_ = value; }

private:
  typedef std::pair<TrajectoryStateOnSurface, double> TsosPP;
  typedef std::pair<bool, std::shared_ptr<G4ErrorTarget>> ErrorTargetPair;

  // Magnetic field
  const MagneticField *theField;

  // Name of the particle whose properties will be used in the propagation
  std::string theParticleName;

  // The Geant4e manager. Does the real propagation
  G4ErrorPropagatorManager *theG4eManager;
  G4ErrorPropagatorData *theG4eData;
  bool consistentBackwardCovariance_ = false;
  bool transportAudit_ = false;
  bool meanEnergyLossJacobian_ = false;
  bool fieldGradientJacobian_ = false;
  bool unquenchedIonizationVariance_ = false;
  bool materialInterfaceJacobian_ = true;
  mutable double maximumInterfaceRelativeCurvatureSigma_ = 0.;
  mutable double maximumInterfaceZCm_ = 0.;
  bool recordTransportJacobian_ = false;
  mutable bool transportJacobianValid_ = false;
  mutable AlgebraicMatrix55 curvilinearTransportJacobian_;
  mutable bool meanCurvatureFlowValid_ = false;
  mutable double startMeanCurvatureFlow_ = std::numeric_limits<double>::quiet_NaN();
  mutable double endMeanCurvatureFlow_ = std::numeric_limits<double>::quiet_NaN();
  mutable bool materialEndpointBoundaryAmbiguous_ = false;
  mutable std::string startMeanCurvatureMaterial_, endMeanCurvatureMaterial_;
  double plimit_;
  double maximumStepLengthMm_;
  double maximumPathLengthCm_;

  // Transform a CMS Reco detector surface into a Geant4 Target for the error
  // propagation
  template <class SurfaceType>
  ErrorTargetPair transformToG4SurfaceTarget(const SurfaceType &pDest, bool moveTargetToEndOfSurface) const;

  // generates the Geant4 name for a particle from the
  // string stored in theParticleName ( set via constructor )
  // and the particle charge.
  // 'mu' as a basis for muon becomes 'mu+' or 'mu-', depening on the charge
  // This method only supports neutral and +/- 1e charges so far
  //
  // returns the generated string
  std::string generateParticleName(int charge) const;

  // flexible method which performs the actual propagation either for a plane or
  // cylinder surface type
  //
  // returns TSOS after the propagation and the path length
  template <class SurfaceType>
  std::pair<TrajectoryStateOnSurface, double> propagateGeneric(const FreeTrajectoryState &ftsStart,
                                                               const SurfaceType &pDest) const;

  // saves the Geant4 propagation direction (Forward or Backward) in the
  // provided variable reference mode and returns true if the propagation
  // direction could be set
  template <class SurfaceType>
  bool configurePropagation(G4ErrorMode &mode,
                            SurfaceType const &pDest,
                            GlobalPoint const &cmsInitPos,
                            GlobalVector const &cmsInitMom) const;

  // special case to determine the propagation direction if the CMS propagation
  // direction 'anyDirection' was set. This method is called by
  // configurePropagation and provides specific implementations for Plane and
  // Cylinder classes
  template <class SurfaceType>
  bool configureAnyPropagation(G4ErrorMode &mode,
                               SurfaceType const &pDest,
                               GlobalPoint const &cmsInitPos,
                               GlobalVector const &cmsInitMom) const;

  // Ensure Geant4 Error propagation is initialized, if not done so, yet
  // if the forceInit parameter is set to true, the initialization is performed,
  // even if already done before.
  // This can be necessary, when Geant4 needs to read in a new MagneticField
  // object, which changed during lumi section crossing
  void ensureGeant4eIsInitilized(bool forceInit) const;

  // returns the name of the SurfaceType. Mostly for debug outputs
  template <class SurfaceType>
  std::string getSurfaceType(SurfaceType const &surface) const;

  void debugReportPlaneSetup(GlobalPoint const &posPlane,
                             HepGeom::Point3D<double> const &surfPos,
                             GlobalVector const &normalPlane,
                             HepGeom::Normal3D<double> const &surfNorm,
                             const Plane &pDest) const;

  template <class SurfaceType>
  void debugReportTrackState(std::string const &currentContext,
                             GlobalPoint const &cmsInitPos,
                             CLHEP::Hep3Vector const &g4InitPos,
                             GlobalVector const &cmsInitMom,
                             CLHEP::Hep3Vector const &g4InitMom,
                             const SurfaceType &pDest) const;
};

#endif
