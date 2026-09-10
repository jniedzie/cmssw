#include "PhysicsTools/ShiftLssGeometry/interface/ExternalVolumeClipper.h"
#include "PhysicsTools/ShiftLssGeometry/interface/BooleanFrameNormalizer.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDG4/Geant4Converter.h"
#include "G4VSolid.hh"
#include "G4SystemOfUnits.hh"
#include "TGeoBBox.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include <iostream>
#include <stdexcept>

void require(bool value, char const* message) {
  if (!value)
    throw std::runtime_error(message);
}

int main() {
  auto& detector = dd4hep::Detector::getInstance();
  dd4hep::sim::Geant4Converter converter(detector);
  auto medium = new TGeoMedium("testMedium", 1, new TGeoMaterial("testMaterial", 0., 0., 0.));
  auto box = new TGeoBBox("ownerBox", 5., 4., 3.);
  auto hole = new TGeoBBox("ownerHole", 1., 1., 10.);
  auto ownerShape = new TGeoCompositeShape("ownerWithHole", new TGeoSubtraction(box, hole));
  auto& owner = *new TGeoVolume("owner", ownerShape, medium);
  TGeoRotation rotation;
  rotation.RotateZ(23.);
  rotation.RotateY(11.);
  TGeoCombiTrans ownerPlacement(5., -2., 1., &rotation);
  shift::ExternalVolumeClipper clipper;
  auto& ownerAssembly = *new TGeoVolumeAssembly("ownerAssembly");
  ownerAssembly.AddNode(&owner, 1, new TGeoHMatrix(ownerPlacement));
  clipper.protect(ownerAssembly, TGeoHMatrix());

  auto externalShape = new TGeoBBox("externalBox", 12., 10., 8.);
  auto& external = *new TGeoVolume("external", externalShape, medium);
  auto& daughter = *new TGeoVolume("daughter", new TGeoBBox("daughterBox", .5, .5, .5), medium);
  external.AddNode(&daughter, 1, new TGeoTranslation(3., 0., 0.));
  unsigned int subtracted = 0, preserved = 0;
  for (double x : {0., 7., 50.}) {
    TGeoCombiTrans externalPlacement(x, 1., -2., &rotation);
    auto* clipped = clipper.clip(&external, externalPlacement);
    require(external.GetNdaughters() == 1 && external.GetShape() == externalShape, "clipping mutated source geometry");
    require(clipped->GetNdaughters() == 1 && clipped->GetNode(0)->GetVolume() != &daughter,
            "nested shared volume was not cloned");
    shift::BooleanFrameNormalizer normalizer;
    auto frame = normalizer.volume(clipped);
    auto* g4 = static_cast<G4VSolid*>(converter.handleSolid(clipped->GetShape()->GetName(), clipped->GetShape()));
    for (double a = -12.13; a < 12.5; a += .61)
      for (double b = -10.17; b < 10.5; b += .67)
        for (double c = -8.19; c < 8.5; c += .71) {
          double local[3] = {a, b, c}, global[3], protectedLocal[3], normalized[3];
          externalPlacement.LocalToMaster(local, global);
          ownerPlacement.MasterToLocal(global, protectedLocal);
          frame.MasterToLocal(local, normalized);
          bool occupied = ownerShape->Contains(protectedLocal);
          bool expected = externalShape->Contains(local) && !occupied;
          bool actual = g4->Inside(G4ThreeVector(normalized[0], normalized[1], normalized[2]) *
                                   (CLHEP::cm / dd4hep::cm)) != kOutside;
          require(clipped->GetShape()->Contains(normalized) == expected && actual == expected,
                  "exact ownership differs between source expectation, ROOT and Geant4");
          subtracted += externalShape->Contains(local) && occupied;
          preserved += expected;
        }
  }
  require(subtracted && preserved, "fixture did not exercise both sides of the ownership boundary");
  std::cout << "Transformed ownership, holes, shared/nested volumes and ROOT/Geant4 agreement passed\n";
}
