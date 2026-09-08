#include <iostream>
#include <array>
#include <stdexcept>
#include "PhysicsTools/ShiftLssGeometry/interface/BooleanFrameNormalizer.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDG4/Geant4Converter.h"
#include "G4VSolid.hh"
#include "G4SystemOfUnits.hh"
#include "TGeoBBox.h"
#include "TGeoManager.h"

int main() {
  auto& detector = dd4hep::Detector::getInstance();
  dd4hep::sim::Geant4Converter converter(detector);
  auto box = new TGeoBBox("frame_box", 2., 3., 4.);
  auto hole = new TGeoBBox("frame_hole", 1., 1., 7.);
  auto rotation = new TGeoRotation();
  rotation->RotateY(31.);
  rotation->RotateZ(17.);
  auto displaced = new TGeoCompositeShape(
      "displaced_difference",
      new TGeoSubtraction(box, hole, new TGeoCombiTrans(12., -7., 4., rotation), new TGeoTranslation(12., -7., 4.)));
  auto combined = new TGeoCompositeShape(
      "nested_union",
      new TGeoUnion(displaced, box, new TGeoTranslation(-3., 8., 2.), new TGeoCombiTrans(-5., 6., 1., rotation)));
  auto clipped = new TGeoCompositeShape(
      "nested_intersection",
      new TGeoIntersection(
          combined, new TGeoBBox("clip_box", 20., 20., 20.), new TGeoCombiTrans(2., 4., -3., rotation), nullptr));
  shift::BooleanFrameNormalizer normalizer;
  for (TGeoShape* original : std::array<TGeoShape*, 3>{displaced, combined, clipped}) {
    auto normalized = normalizer.solid(original);
    auto converted = static_cast<G4VSolid*>(converter.handleSolid(normalized.shape->GetName(), normalized.shape));
    unsigned int inside = 0;
    for (double x = -24.13; x < 25.; x += 0.77)
      for (double y = -24.27; y < 25.; y += 0.81)
        for (double z = -24.39; z < 25.; z += 0.83) {
          double point[3] = {x, y, z}, local[3];
          normalized.frame.MasterToLocal(point, local);
          bool expected = original->Contains(point);
          inside += expected;
          bool actual =
              converted->Inside(G4ThreeVector(local[0], local[1], local[2]) * (CLHEP::cm / dd4hep::cm)) != kOutside;
          if (expected != normalized.shape->Contains(local) || expected != actual)
            throw std::runtime_error(std::string("Boolean frame occupancy mismatch: ") + original->GetName());
        }
    if (!inside)
      throw std::runtime_error("Empty Boolean test fixture");
  }
  std::cout << "Nested union, subtraction and intersection match ROOT and Geant4\n";
}
