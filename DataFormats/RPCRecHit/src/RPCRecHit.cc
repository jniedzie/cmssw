/*
 *  See header file for a description of this class.
 *
 *  \author M. Maggi -- INFN Bari
 */

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <algorithm>
#include <cstdint>

RPCRecHit::RPCRecHit(const RPCDetId& rpcId, int bx)
    : RecHit2DLocalPos(rpcId),
      theRPCId(rpcId),
      theBx(bx),
      theFirstStrip(99),
      theClusterSize(99),
      theLocalPosition(),
      theLocalError(),
      theTime(0),
      theTimeError(-1) {}

RPCRecHit::RPCRecHit()
    : RecHit2DLocalPos(),
      theRPCId(),
      theBx(99),
      theFirstStrip(99),
      theClusterSize(99),
      theLocalPosition(),
      theLocalError(),
      theTime(0),
      theTimeError(-1) {}

RPCRecHit::RPCRecHit(const RPCDetId& rpcId, int bx, const LocalPoint& pos)
    : RecHit2DLocalPos(rpcId),
      theRPCId(rpcId),
      theBx(bx),
      theFirstStrip(99),
      theClusterSize(99),
      theLocalPosition(pos),
      theTime(0),
      theTimeError(-1) {
  float stripResolution = 3.0;  //cm  this sould be taken from trimmed cluster size times strip size
                                //    taken out from geometry service i.e. topology
  theLocalError = LocalError(stripResolution * stripResolution, 0., 0.);  //FIXME: is it really needed?
}

// Constructor from a local position and error, wireId and digi time.
RPCRecHit::RPCRecHit(const RPCDetId& rpcId, int bx, const LocalPoint& pos, const LocalError& err)
    : RecHit2DLocalPos(rpcId),
      theRPCId(rpcId),
      theBx(bx),
      theFirstStrip(99),
      theClusterSize(99),
      theLocalPosition(pos),
      theLocalError(err),
      theTime(0),
      theTimeError(-1) {}

// Constructor from a local position and error, wireId, bx and cluster size.
RPCRecHit::RPCRecHit(
    const RPCDetId& rpcId, int bx, int firstStrip, int clustSize, const LocalPoint& pos, const LocalError& err)
    : RecHit2DLocalPos(rpcId),
      theRPCId(rpcId),
      theBx(bx),
      theFirstStrip(firstStrip),
      theClusterSize(clustSize),
      theLocalPosition(pos),
      theLocalError(err),
      theTime(0),
      theTimeError(-1) {}

// Destructor
RPCRecHit::~RPCRecHit() {}

RPCRecHit* RPCRecHit::clone() const { return new RPCRecHit(*this); }

bool RPCRecHit::sharesInput(const TrackingRecHit* other, SharedInputType what) const {
  auto const* rpc = dynamic_cast<const RPCRecHit*>(other);
  if (!rpc || geographicalId() != rpc->geographicalId() || BunchX() != rpc->BunchX())
    return false;

  std::int64_t first = firstClusterStrip(), otherFirst = rpc->firstClusterStrip();
  std::int64_t size = clusterSize(), otherSize = rpc->clusterSize();
  // Constructors without cluster information use (99,99). Do not invent an identity.
  if (first < 1 || otherFirst < 1 || size < 1 || otherSize < 1 || (first == 99 && size == 99) ||
      (otherFirst == 99 && otherSize == 99))
    throw cms::Exception("UnsupportedRPCInputIdentity") << "RPC cluster has missing or invalid strip metadata";

  if (what == TrackingRecHit::all)
    return first == otherFirst && size == otherSize;
  if (what == TrackingRecHit::some)
    return std::max(first, otherFirst) <= std::min(first + size - 1, otherFirst + otherSize - 1);
  throw cms::Exception("UnsupportedRPCInputIdentity") << "Unknown shared-input policy";
}

// Access to component RecHits.
// No components rechits: it returns a null vector
std::vector<const TrackingRecHit*> RPCRecHit::recHits() const {
  std::vector<const TrackingRecHit*> nullvector;
  return nullvector;
}

// Non-const access to component RecHits.
// No components rechits: it returns a null vector
std::vector<TrackingRecHit*> RPCRecHit::recHits() {
  std::vector<TrackingRecHit*> nullvector;
  return nullvector;
}

// Comparison operator, based on the wireId and the digi time
bool RPCRecHit::operator==(const RPCRecHit& hit) const { return this->geographicalId() == hit.geographicalId(); }

// The ostream operator
std::ostream& operator<<(std::ostream& os, const RPCRecHit& hit) {
  os << "pos: " << hit.localPosition().x();
  os << " +/- " << sqrt(hit.localPositionError().xx());
  return os;
}
