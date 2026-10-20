// Code style: human-code.
#include "PhysicsTools/ShiftMuonSegments/interface/CSCFitCoverage.h"

#include <cstdlib>
#include <iostream>

namespace {
  void require(bool condition, char const* message) {
    if (!condition) {
      std::cerr << message << '\n';
      std::exit(1);
    }
  }
}

int main() {
  shift::CSCFitCoverage coverage;
  coverage.addMeasurement(false, CSCDetId(1, 2, 1, 7, 3));
  coverage.addMeasurement(true, DetId(DetId::Tracker, 1));
  coverage.addMeasurement(true, DetId(DetId::Muon, MuonSubdetId::DT));
  coverage.addMeasurement(true, DetId(DetId::Muon, MuonSubdetId::RPC));
  require(coverage.hits == 0 && coverage.endcapMask == 0 && coverage.stationMask == 0,
          "invalid and non-CSC measurements must not supply CSC coverage");

  coverage.addMeasurement(true, CSCDetId(1, 2, 1, 7, 3));
  coverage.addMeasurement(true, CSCDetId(1, 2, 1, 7, 4));
  require(coverage.hits == 2 && coverage.endcapMask == 1 && coverage.stationMask == 2,
          "two layers in one +z station supply one station bit");
  coverage.addMeasurement(true, CSCDetId(2, 3, 1, 6, 5));
  require(coverage.hits == 3 && coverage.endcapMask == 3 && coverage.stationMask == 66,
          "opposite endcap station bits must remain distinct");

  shift::CSCFitCoverage allStations;
  for (int endcap = 1; endcap <= 2; ++endcap)
    for (int station = 1; station <= 4; ++station)
      allStations.addMeasurement(true, CSCDetId(endcap, station, 1, 1, 1));
  require(allStations.hits == 8 && allStations.endcapMask == 3 && allStations.stationMask == 255,
          "all eight CSC station/endcap combinations must be represented");

  for (unsigned int input = 0; input != 4; ++input)
    for (unsigned int fitted = 0; fitted != 4; ++fitted) {
      require(shift::passesCSCFitCoverage(false, input, fitted),
              "disabled guard must accept every native coverage combination");
      bool const expected = input != 3 || fitted == 3;
      require(shift::passesCSCFitCoverage(true, input, fitted) == expected,
              "enabled guard must reject only lost two-endcap CSC input coverage");
    }
  require(!shift::passesCSCFitCoverage(true, 3, 1), "34 surviving +z hits cannot replace a both-endcap input");
  require(!shift::passesCSCFitCoverage(true, 3, 2), "the corresponding -z-only failure must also be rejected");
  std::cout << "CSC fitted coverage checks passed\n";
}
