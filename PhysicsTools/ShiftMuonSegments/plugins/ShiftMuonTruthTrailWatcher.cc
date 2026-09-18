#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Watcher/interface/SimProducer.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace {
  constexpr char const* prefix = "shiftMuonTruthTrail";

  std::string instance(std::string const& suffix) { return std::string(prefix) + suffix; }
}  // namespace

// Simulation-only diagnostic. It stores a sparse Geant4 state trail for
// primary inward muons. The products are deliberately primitive parallel
// vectors so a bounded diagnostic can be read without a new persistent data
// format. All coordinates use cm, momenta GeV and time ns.
class ShiftMuonTruthTrailWatcher : public SimProducer,
                                   public Observer<const BeginOfEvent*>,
                                   public Observer<const G4Step*> {
public:
  explicit ShiftMuonTruthTrailWatcher(edm::ParameterSet const& parameters) {
    auto const config = parameters.getParameter<edm::ParameterSet>("ShiftMuonTruthTrailWatcher");
    interval_ = config.getUntrackedParameter<double>("intervalCm", 100.0) * cm;
    maxCheckpoints_ = config.getUntrackedParameter<unsigned int>("maxCheckpoints", 10000);
    inwardOnly_ = config.getUntrackedParameter<bool>("inwardOnly", true);
    if (!(interval_ > 0.) || !std::isfinite(interval_))
      throw cms::Exception("Configuration") << "ShiftMuonTruthTrailWatcher intervalCm must be finite and positive";
    if (maxCheckpoints_ == 0)
      throw cms::Exception("Configuration") << "ShiftMuonTruthTrailWatcher maxCheckpoints must be positive";

    setMT(true);
    for (auto const& suffix : {"TrackId", "ParentId", "PdgId", "Step", "Kind"})
      produces<std::vector<int>>(instance(suffix));
    for (auto const& suffix : {"X", "Y", "Z", "Px", "Py", "Pz", "GlobalTime", "TrackLength",
                               "StepLength", "Density", "RadiationLength", "NuclearInteractionLength",
                               "CumulativeX0", "CumulativeInteractionLengths", "CumulativeEnergyLoss",
                               "EnergyDeposit"})
      produces<std::vector<float>>(instance(suffix));
    for (auto const& suffix : {"Volume", "Material", "Process"})
      produces<std::vector<std::string>>(instance(suffix));
    produces<int>(instance("Truncated"));
  }

  void produce(edm::Event& event, edm::EventSetup const&) override {
    put(event, "TrackId", std::move(trackId_));
    put(event, "ParentId", std::move(parentId_));
    put(event, "PdgId", std::move(pdgId_));
    put(event, "Step", std::move(step_));
    put(event, "Kind", std::move(kind_));
    put(event, "X", std::move(x_));
    put(event, "Y", std::move(y_));
    put(event, "Z", std::move(z_));
    put(event, "Px", std::move(px_));
    put(event, "Py", std::move(py_));
    put(event, "Pz", std::move(pz_));
    put(event, "GlobalTime", std::move(globalTime_));
    put(event, "TrackLength", std::move(trackLength_));
    put(event, "StepLength", std::move(stepLength_));
    put(event, "Density", std::move(density_));
    put(event, "RadiationLength", std::move(radiationLength_));
    put(event, "NuclearInteractionLength", std::move(nuclearInteractionLength_));
    put(event, "CumulativeX0", std::move(cumulativeX0Values_));
    put(event, "CumulativeInteractionLengths", std::move(cumulativeInteractionLengthValues_));
    put(event, "CumulativeEnergyLoss", std::move(cumulativeEnergyLoss_));
    put(event, "EnergyDeposit", std::move(energyDeposit_));
    put(event, "Volume", std::move(volume_));
    put(event, "Material", std::move(material_));
    put(event, "Process", std::move(process_));
    event.put(std::make_unique<int>(truncated_ ? 1 : 0), instance("Truncated"));
  }

private:
  template <typename T>
  void put(edm::Event& event, char const* suffix, std::vector<T>&& values) {
    event.put(std::make_unique<std::vector<T>>(std::move(values)), instance(suffix));
  }

  void update(BeginOfEvent const*) override {
    trackId_.clear();
    parentId_.clear();
    pdgId_.clear();
    step_.clear();
    kind_.clear();
    for (auto* values : {&x_, &y_, &z_, &px_, &py_, &pz_, &globalTime_, &trackLength_, &stepLength_, &density_,
                         &radiationLength_, &nuclearInteractionLength_, &cumulativeX0Values_,
                         &cumulativeInteractionLengthValues_, &cumulativeEnergyLoss_, &energyDeposit_})
      values->clear();
    volume_.clear();
    material_.clear();
    process_.clear();
    acceptedTracks_.clear();
    rejectedTracks_.clear();
    initialKineticEnergy_.clear();
    cumulativeX0_.clear();
    cumulativeInteractionLengths_.clear();
    nextPeriodicLength_.clear();
    truncated_ = false;
  }

  bool accept(G4Track const& track) {
    int const id = track.GetTrackID();
    if (acceptedTracks_.count(id))
      return true;
    if (rejectedTracks_.count(id))
      return false;
    bool accepted = track.GetParentID() == 0 && std::abs(track.GetDefinition()->GetPDGEncoding()) == 13;
    if (accepted && inwardOnly_)
      accepted = track.GetVertexPosition().z() * track.GetVertexMomentumDirection().z() < 0.;
    (accepted ? acceptedTracks_ : rejectedTracks_).insert(id);
    if (accepted) {
      initialKineticEnergy_[id] = track.GetKineticEnergy();
      cumulativeX0_[id] = 0.;
      cumulativeInteractionLengths_[id] = 0.;
      nextPeriodicLength_[id] = interval_;
    }
    return accepted;
  }

  void update(G4Step const* step) override {
    if (!step || truncated_)
      return;
    G4Track const& track = *step->GetTrack();
    if (!accept(track))
      return;
    G4StepPoint const* pre = step->GetPreStepPoint();
    G4StepPoint const* post = step->GetPostStepPoint();
    if (!pre || !post)
      return;

    int const id = track.GetTrackID();
    G4Material const* preMaterial = pre->GetMaterial();
    double const length = step->GetStepLength();
    if (track.GetCurrentStepNumber() == 1)
      initialKineticEnergy_[id] = pre->GetKineticEnergy();
    double const preX0 = cumulativeX0_[id];
    double const preLambda = cumulativeInteractionLengths_[id];
    if (preMaterial && preMaterial->GetRadlen() > 0.)
      cumulativeX0_[id] += length / preMaterial->GetRadlen();
    if (preMaterial && preMaterial->GetNuclearInterLength() > 0.)
      cumulativeInteractionLengths_[id] += length / preMaterial->GetNuclearInterLength();

    double const preTrackLength = std::max(0., track.GetTrackLength() - length);
    if (track.GetCurrentStepNumber() == 1)
      append(track, *pre, 0, preTrackLength, 0., preX0, preLambda);

    bool const boundary =
        pre->GetPhysicalVolume() != post->GetPhysicalVolume() || pre->GetMaterial() != post->GetMaterial();
    if (boundary) {
      append(track, *pre, 2, preTrackLength, 0., preX0, preLambda);
      append(track, *post, 3, track.GetTrackLength(), length, cumulativeX0_[id], cumulativeInteractionLengths_[id]);
    }

    if (track.GetTrackLength() >= nextPeriodicLength_[id]) {
      append(track, *post, 1, track.GetTrackLength(), length, cumulativeX0_[id], cumulativeInteractionLengths_[id]);
      while (nextPeriodicLength_[id] <= track.GetTrackLength())
        nextPeriodicLength_[id] += interval_;
    }

    if (track.GetTrackStatus() != fAlive || post->GetPhysicalVolume() == nullptr)
      append(track, *post, 4, track.GetTrackLength(), length, cumulativeX0_[id], cumulativeInteractionLengths_[id]);
  }

  void append(G4Track const& track,
              G4StepPoint const& point,
              int checkpointKind,
              double trackLength,
              double stepLength,
              double cumulativeX0,
              double cumulativeInteractionLengths) {
    if (trackId_.size() >= maxCheckpoints_) {
      truncated_ = true;
      return;
    }
    auto const& position = point.GetPosition();
    auto const& momentum = point.GetMomentum();
    G4Material const* material = point.GetMaterial();
    G4VPhysicalVolume const* volume = point.GetPhysicalVolume();
    G4VProcess const* process = point.GetProcessDefinedStep();
    trackId_.push_back(track.GetTrackID());
    parentId_.push_back(track.GetParentID());
    pdgId_.push_back(track.GetDefinition()->GetPDGEncoding());
    step_.push_back(track.GetCurrentStepNumber());
    kind_.push_back(checkpointKind);
    x_.push_back(position.x() / cm);
    y_.push_back(position.y() / cm);
    z_.push_back(position.z() / cm);
    px_.push_back(momentum.x() / GeV);
    py_.push_back(momentum.y() / GeV);
    pz_.push_back(momentum.z() / GeV);
    globalTime_.push_back(point.GetGlobalTime() / ns);
    trackLength_.push_back(trackLength / cm);
    stepLength_.push_back(stepLength / cm);
    density_.push_back(material ? material->GetDensity() / (g / cm3) : 0.f);
    radiationLength_.push_back(material ? material->GetRadlen() / cm : 0.f);
    nuclearInteractionLength_.push_back(material ? material->GetNuclearInterLength() / cm : 0.f);
    cumulativeX0Values_.push_back(cumulativeX0);
    cumulativeInteractionLengthValues_.push_back(cumulativeInteractionLengths);
    cumulativeEnergyLoss_.push_back((initialKineticEnergy_[track.GetTrackID()] - point.GetKineticEnergy()) / GeV);
    energyDeposit_.push_back(stepLength > 0. && track.GetStep() ? track.GetStep()->GetTotalEnergyDeposit() / GeV : 0.f);
    volume_.push_back(volume ? volume->GetName() : "outside-world");
    material_.push_back(material ? material->GetName() : "none");
    process_.push_back(process ? process->GetProcessName() : "none");
  }

  double interval_;
  std::size_t maxCheckpoints_;
  bool inwardOnly_;
  bool truncated_ = false;
  std::set<int> acceptedTracks_;
  std::set<int> rejectedTracks_;
  std::map<int, double> initialKineticEnergy_;
  std::map<int, double> cumulativeX0_;
  std::map<int, double> cumulativeInteractionLengths_;
  std::map<int, double> nextPeriodicLength_;
  std::vector<int> trackId_, parentId_, pdgId_, step_, kind_;
  std::vector<float> x_, y_, z_, px_, py_, pz_, globalTime_, trackLength_, stepLength_, density_, radiationLength_,
      nuclearInteractionLength_, cumulativeX0Values_, cumulativeInteractionLengthValues_, cumulativeEnergyLoss_,
      energyDeposit_;
  std::vector<std::string> volume_, material_, process_;
};

DEFINE_SIMWATCHER(ShiftMuonTruthTrailWatcher);
