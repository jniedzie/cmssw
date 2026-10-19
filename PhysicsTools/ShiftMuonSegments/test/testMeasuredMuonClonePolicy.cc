// Code style: AI-code. AI-generated regression for measured-only clone policies.
#include "PhysicsTools/ShiftMuonSegments/interface/MeasuredMuonClonePolicy.h"
#include <stdexcept>
#include <iostream>
#include <numeric>

void require(bool condition) {
  if (!condition)
    throw std::runtime_error("Measured clone policy regression failed");
}

int main() {
  using shift::MuonHitFingerprint;
  std::vector<MuonHitFingerprint> a{{11, 0.f, 0.f}, {12, 1.f, 0.f}, {13, 2.f, 0.f}};
  std::vector<MuonHitFingerprint> disjoint{{11, 1.f, 0.f}, {12, 2.f, 0.f}, {13, 3.f, 0.f}};
  require(!shift::shareMuonHitFingerprints(a, disjoint, 2, .5));
  require(shift::shareMuonHitFingerprints(a, a, 2, 1.));
  auto repeated = std::vector<MuonHitFingerprint>(40, a.front());
  require(!shift::shareMuonHitFingerprints(repeated, a, 2, .5));
  require(!shift::shareMuonHitFingerprints(a, repeated, 2, .5));
  // One-to-one matching must reassign an ambiguous first hit. A greedy match
  // would find one shared record; the maximum matching contains two.
  std::vector<MuonHitFingerprint> ambiguousA{{11, 0.f, 0.f}, {11, -.09f, 0.f}};
  std::vector<MuonHitFingerprint> ambiguousB{{11, -.01f, 0.f}, {11, .09f, 0.f}};
  require(shift::shareMuonHitFingerprints(ambiguousA, ambiguousB, 2, 1.));
  require(shift::shareMuonHitFingerprints(ambiguousB, ambiguousA, 2, 1.));
  require(!shift::shareMuonHitFingerprints({}, a, 1, .5));
  auto twice = a; twice.insert(twice.end(), a.begin(), a.end());
  require(shift::distinctMuonHitFingerprints(twice).size() == a.size());
  require(shift::shareMuonHitFingerprints(twice, a, 2, 1.));

  using shift::BothMuonCloneQuality;
  std::vector<BothMuonCloneQuality> quality{
      {true, 1, 2.25, 143., 8, 74, 1, 0},  // strict, repeated inputs
      {true, 1, .57, 63., 8, 34, 2, 0},   // same measured Both track, better fit
      {false, 1, .01, 63., 8, 34, 2, 1}, // invalid directional refit
      {true, -1, .01, 63., 8, 34, 2, 2}, // wrong native timing direction
      {true, 1, std::numeric_limits<double>::quiet_NaN(), 63., 8, 34, 2, 3}};
  std::vector<unsigned int> order(quality.size());
  std::iota(order.begin(), order.end(), 0);
  unsigned int permutations = 0;
  do {
    auto winner = order.front();
    for (auto member : order)
      if (shift::betterBothMuonClone(quality[member], quality[winner]))
        winner = member;
    require(winner == 1);
    for (auto i : order) {
      require(!shift::betterBothMuonClone(quality[i], quality[i]));
      for (auto j : order)
        for (auto k : order)
          if (shift::betterBothMuonClone(quality[i], quality[j]) && shift::betterBothMuonClone(quality[j], quality[k]))
            require(shift::betterBothMuonClone(quality[i], quality[k]));
    }
    ++permutations;
  } while (std::next_permutation(order.begin(), order.end()));
  auto tied = quality[1]; tied.sourceIndex = 5;
  require(shift::betterBothMuonClone(quality[1], tied));
  std::cout << "Measured clone policy: shared-hit invariants and " << permutations
            << " order permutations pass\n";
}
