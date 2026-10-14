// Code style: AI-code. AI-generated regression test of RPC readout input identity.
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <cassert>
#include <memory>

namespace {
  RPCRecHit cluster(unsigned roll, int bx, int first, int size) {
    return RPCRecHit(RPCDetId(roll), bx, first, size, LocalPoint(2, 3, 0), LocalError(1, 0, 1));
  }

  bool shares(const TrackingRecHit& first, const TrackingRecHit* second, TrackingRecHit::SharedInputType policy) {
    return first.sharesInput(second, policy);
  }
}

int main() {
  auto first = cluster(637567033, 0, 10, 3);
  auto identical = cluster(637567033, 0, 10, 3);
  auto contained = cluster(637567033, 0, 11, 1);
  auto overlap = cluster(637567033, 0, 12, 3);
  auto adjacent = cluster(637567033, 0, 13, 3);
  auto otherBX = cluster(637567033, 1, 10, 3);
  auto otherRoll = cluster(637567034, 0, 10, 3);

  assert(shares(first, &identical, TrackingRecHit::all));
  assert(shares(first, &first, TrackingRecHit::all));
  assert(shares(first, &contained, TrackingRecHit::some));
  assert(shares(contained, &first, TrackingRecHit::some));
  assert(!shares(first, &contained, TrackingRecHit::all));
  assert(shares(first, &overlap, TrackingRecHit::some));
  assert(!shares(first, &overlap, TrackingRecHit::all));
  assert(!shares(first, &adjacent, TrackingRecHit::some));
  assert(!shares(first, &otherBX, TrackingRecHit::some));
  assert(!shares(first, &otherRoll, TrackingRecHit::some));
  assert(!shares(first, nullptr, TrackingRecHit::some));
  std::unique_ptr<RPCRecHit> copy(first.clone());
  assert(shares(first, copy.get(), TrackingRecHit::all));

  // Missing strip metadata must remain visible rather than merge arbitrary hits.
  RPCRecHit missing(RPCDetId(637567033), 0);
  bool rejected = false;
  try {
    shares(first, &missing, TrackingRecHit::some);
  } catch (cms::Exception const& exception) {
    rejected = exception.category() == "UnsupportedRPCInputIdentity";
  }
  assert(rejected);
}
