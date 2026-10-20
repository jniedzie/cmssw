// Code style: human-code.
#include "PhysicsTools/ShiftMuonSegments/interface/ExactMuonFitInput.h"

#include <stdexcept>
#include <limits>
#include <vector>

namespace {
  void require(bool condition) {
    if (!condition)
      throw std::runtime_error("exact CSC copy regression failed");
  }
  struct HitValues {
    CSCDetId id{1, 2, 1, 3, 4};
    float x = 1.f;
    float z = 0.f;
    float xy = .01f;
    float adc = 10.f;
    float time = 5.f;
    float stripPosition = .2f;
    float stripError = .1f;
    float energy = 3.f;
    int strip = 10;
    int wire = 12;
    int quality = 1;
    short badStrip = 0;
    short badWire = 0;
    int wireTime = 100;
  };

  CSCRecHit2D makeHit(HitValues const& v = {}) {
    CSCRecHit2D::ADCContainer adcs;
    std::vector<float> samples{v.adc, 11.f, 12.f, 13.f};
    adcs.put(v.strip, samples.begin(), samples.end());
    return CSCRecHit2D(v.id, LocalPoint(v.x, 2.f, v.z), LocalError(.1f, v.xy, .2f), {v.strip}, adcs,
                       {v.wire}, v.time, v.stripPosition, v.stripError, v.quality, v.badStrip,
                       v.badWire, v.wireTime, v.energy);
  }

  CSCRecHit2D makeCluster(std::vector<int> const& strips) {
    CSCRecHit2D::ADCContainer adcs;
    std::vector<float> samples{10.f, 11.f, 12.f, 13.f};
    for (auto strip : strips)
      adcs.put(strip, samples.begin(), samples.end());
    return CSCRecHit2D(CSCDetId(1, 2, 1, 3, 4), LocalPoint(1.f, 2.f, 0.f), LocalError(.1f, .01f, .2f),
                       strips, adcs, {12}, 5.f, .2f, .1f, 1, 0, 0, 100, 3.f);
  }

  // A different dynamic class with a CSC DetId must never reach CSC's cast.
  class OtherHit final : public RecHit2DLocalPos {
  public:
    OtherHit() : RecHit2DLocalPos(CSCDetId(1, 2, 1, 3, 4)) {}
    OtherHit* clone() const override { return new OtherHit(*this); }
    LocalPoint localPosition() const override { return LocalPoint(1.f, 2.f, 0.f); }
    LocalError localPositionError() const override { return LocalError(.1f, .01f, .2f); }
  };
}  // namespace

int main() {
  auto const original = makeHit();
  auto const copy = original;
  require(shift::isExactCSCMeasurementCopy(original, copy));
  require(shift::isExactCSCMeasurementCopy(copy, original));
  auto different = [&](HitValues const& values) { require(!shift::isExactCSCMeasurementCopy(original, makeHit(values))); };
  HitValues v;
  v.x = std::nextafter(v.x, 2.f); different(v);
  v = {}; v.xy = std::nextafter(v.xy, 1.f); different(v);
  v = {}; v.id = CSCDetId(1, 2, 1, 3, 5); different(v);
  v = {}; v.z = 1.f; different(v);
  v = {}; v.strip = 11; different(v);
  v = {}; v.wire = 13; different(v);
  v = {}; v.wire |= (1 << 16); different(v);
  v = {}; v.adc = std::nextafter(v.adc, 20.f); different(v);
  v = {}; v.time = std::nextafter(v.time, 10.f); different(v);
  v = {}; v.wireTime = 101; different(v);
  v = {}; v.stripPosition = .3f; different(v);
  v = {}; v.stripError = .2f; different(v);
  v = {}; v.energy = 4.f; different(v);
  v = {}; v.quality = 2; different(v);
  v = {}; v.badStrip = 1; different(v);
  v = {}; v.badWire = 1; different(v);
  v = {}; v.x = std::numeric_limits<float>::quiet_NaN(); different(v);
  v = {}; v.xy = std::numeric_limits<float>::infinity(); different(v);
  v = {}; v.adc = std::numeric_limits<float>::quiet_NaN(); different(v);
  auto const nonfinite = makeHit(v);
  require(!shift::isExactCSCMeasurementCopy(nonfinite, nonfinite));
  v = {}; v.xy = std::numeric_limits<float>::infinity();
  auto const infiniteCovariance = makeHit(v);
  require(!shift::isExactCSCMeasurementCopy(infiniteCovariance, infiniteCovariance));
  require(!shift::isExactCSCMeasurementCopy(original, CSCRecHit2D()));
  OtherHit other;
  require(!shift::isExactCSCMeasurementCopy(original, other));
  require(!shift::isExactCSCMeasurementCopy(other, original));

  auto const partial = makeCluster({10, 11});
  require(original.sharesInput(&partial, TrackingRecHit::some));
  require(!original.sharesInput(&partial, TrackingRecHit::all));
  require(!shift::isExactCSCMeasurementCopy(original, partial));
  require(!shift::isExactCSCMeasurementCopy(partial, original));

  // Distinct pulses in the same layer remain; keep the original copy first.
  v = {}; v.time += 1.f;
  auto const secondPulse = makeHit(v);
  std::vector<TrackingRecHit const*> inputs{&original, &copy, &secondPulse, &copy};
  std::vector<TrackingRecHit const*> kept;
  for (auto const* hit : inputs) {
    bool duplicate = false;
    for (auto const* previous : kept)
      duplicate = duplicate || shift::isExactCSCMeasurementCopy(*hit, *previous);
    if (!duplicate)
      kept.push_back(hit);
  }
  require(kept.size() == 2 && kept[0] == &original && kept[1] == &secondPulse);
  require(!shift::isExactCSCMeasurementCopy(*kept[0], *kept[1]));
}
