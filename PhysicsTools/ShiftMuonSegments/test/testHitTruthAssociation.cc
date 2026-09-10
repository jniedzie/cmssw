#include "PhysicsTools/ShiftMuonSegments/interface/HitTruthAssociation.h"
#include <cassert>
#include <iostream>

int main() {
  std::vector<shift::TruthMeasurement> hits;
  std::vector<shift::TruthCrossing> truth;
  for (unsigned int layer = 1; layer <= 4; ++layer) {
    hits.push_back({layer, 1., 2., 0., .01, 0., .01});
    truth.push_back({layer, 7, 0., 2., -1., 1., 0., 1.});
    truth.push_back({layer, 9, 10., 20., 0., 0., 0., 1.});
  }
  auto result = shift::associateHits(hits, truth);
  assert(result.track == 7 && result.total == 4 && result.matched == 4);
  hits.push_back(hits.front());
  assert(shift::associateHits(hits, truth).total == 4);
  truth.push_back({1, 8, 1., 2., 0., 0., 0., 1.});
  result = shift::associateHits(hits, truth);
  assert(result.track == 7 && result.ambiguous == 1 && result.purity == .75);
  truth.push_back({2, 8, 1., 2., 0., 0., 0., 1.});
  assert(shift::associateHits(hits, truth).track == -1);
  assert(shift::associateHits(hits, {}).track == -1);
  hits.resize(2);
  assert(shift::associateHits(hits, truth).track == -1);
  for (auto& h : hits) h.xx = 0.;
  assert(shift::associateHits(hits, truth).track == -1);
  std::cout << "Layer projection, ambiguity, purity, duplicate and missing-hit checks passed\n";
}
