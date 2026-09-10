#ifndef PhysicsTools_ShiftLssGeometry_ExternalVolumeClipper_h
#define PhysicsTools_ShiftLssGeometry_ExternalVolumeClipper_h

#include <array>
#include <algorithm>
#include <string>
#include <vector>
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

namespace shift {
  // Merge policy: existing daughter envelopes own their space, including
  // their internal air and daughters. External material fills the complement.
  // Work from actual solids/placements; bounding boxes only reject disjoint
  // candidates and never define the physical subtraction.
  class ExternalVolumeClipper {
  public:
    struct PlacedSolid {
      TGeoShape* shape;
      TGeoHMatrix toMother;
    };

    void protect(TGeoVolume const& volume, TGeoHMatrix const& toMother) {
      if (!volume.IsAssembly()) {
        protected_.push_back({volume.GetShape(), toMother});
        return;
      }
      for (int i = 0; i < volume.GetNdaughters(); ++i) {
        auto const* node = volume.GetNode(i);
        auto placement = toMother;
        placement.Multiply(node->GetMatrix());
        protect(*node->GetVolume(), placement);
      }
    }

    TGeoVolume* clip(TGeoVolume* original, TGeoHMatrix const& toMother) {
      TGeoShape* shape = original->GetShape();
      auto const originalBounds = bounds(shape, toMother);
      if (!original->IsAssembly()) {
        for (auto const& item : protected_) {
          auto const otherBounds = bounds(item.shape, item.toMother);
          bool disjoint = false;
          for (unsigned int axis = 0; axis < 3; ++axis)
            disjoint |= originalBounds[2 * axis + 1] < otherBounds[2 * axis] ||
                        otherBounds[2 * axis + 1] < originalBounds[2 * axis];
          if (disjoint)
            continue;
          auto* relative = new TGeoHMatrix(toMother.Inverse());
          relative->Multiply(&item.toMother);
          auto name = std::string(original->GetName()) + "_cms_exclusion_" + std::to_string(clipped_++);
          shape = new TGeoCompositeShape(name.c_str(), new TGeoSubtraction(shape, item.shape, nullptr, relative));
        }
      }
      // Cloning per placement is necessary when a shared logical volume is
      // placed at different positions relative to the protected geometry.
      auto* result = original->CloneVolume();
      auto name = std::string(original->GetName()) + "_cms_owned_" + std::to_string(volumes_++);
      result->SetName(name.c_str());
      if (!original->IsAssembly())
        result->SetShape(shape);
      while (result->GetNdaughters())
        result->RemoveNode(result->GetNode(0));
      for (int i = 0; i < original->GetNdaughters(); ++i) {
        auto const* node = original->GetNode(i);
        auto placement = toMother;
        placement.Multiply(node->GetMatrix());
        result->AddNode(clip(node->GetVolume(), placement), node->GetNumber(), new TGeoHMatrix(*node->GetMatrix()));
      }
      if (result->IsAssembly())
        result->GetShape()->ComputeBBox();
      if (result->GetNdaughters())
        result->Voxelize("");
      return result;
    }

    unsigned int subtractions() const { return clipped_; }

  private:
    static std::array<double, 6> bounds(TGeoShape const* shape, TGeoHMatrix const& placement) {
      std::array<double, 3> low{}, high{};
      for (int axis = 0; axis < 3; ++axis)
        shape->GetAxisRange(axis + 1, low[axis], high[axis]);
      std::array<double, 6> result{};
      for (unsigned int corner = 0; corner < 8; ++corner) {
        double local[3], global[3];
        for (unsigned int axis = 0; axis < 3; ++axis)
          local[axis] = corner & (1 << axis) ? high[axis] : low[axis];
        placement.LocalToMaster(local, global);
        for (unsigned int axis = 0; axis < 3; ++axis) {
          result[2 * axis] = corner ? std::min(result[2 * axis], global[axis]) : global[axis];
          result[2 * axis + 1] = corner ? std::max(result[2 * axis + 1], global[axis]) : global[axis];
        }
      }
      return result;
    }
    std::vector<PlacedSolid> protected_;
    unsigned int clipped_ = 0;
    unsigned int volumes_ = 0;
  };
}  // namespace shift
#endif
