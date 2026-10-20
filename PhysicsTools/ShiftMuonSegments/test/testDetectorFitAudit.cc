// Code style: human-code.
#include "PhysicsTools/ShiftMuonSegments/interface/DetectorFitAudit.h"
#include <cassert>
#include <limits>
int main() {
  shift::DetectorTruthBracket bracket;
  std::array<double, 3> x{}, p{};
  double weight = -1.;
  assert(!bracket.interpolate(x, p, weight));
  bracket.add({{1., 2., -2.}, {3., 0., 4.}, -2., 0});
  assert(!bracket.interpolate(x, p, weight));
  bracket.add({{5., 4., 2.}, {5., 0., 4.}, 2., 1});
  assert(bracket.interpolate(x, p, weight));
  assert(weight == .5 && x[0] == 3. && x[1] == 3. && x[2] == 0. && p[0] == 4.);
  bracket.add({{0., 0., 0.}, {0., 0., 8.}, 0., 0});
  assert(bracket.interpolate(x, p, weight) && weight == 0. && p[2] == 8.);
  bracket.add({{}, {}, std::numeric_limits<double>::quiet_NaN(), 0});
  assert(bracket.interpolate(x, p, weight) && p[2] == 8.);
  shift::DetectorTruthBracket degenerate;
  degenerate.add({{}, {}, 0., 0});
  assert(!degenerate.interpolate(x, p, weight));
  shift::DetectorTruthBracket sameSide;
  sameSide.add({{}, {1., 1., 1.}, 1., 0});
  sameSide.add({{}, {1., 1., 1.}, 2., 1});
  assert(!sameSide.interpolate(x, p, weight));
}
