#ifndef PhysicsTools_ShiftLssGeometry_BooleanFrameNormalizer_h
#define PhysicsTools_ShiftLssGeometry_BooleanFrameNormalizer_h

#include <map>
#include <stdexcept>
#include <string>
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

namespace shift {
  // DDG4 converts only the right operand's transform. Re-express each Boolean
  // in its left operand's frame, carrying that frame through enclosing solids
  // and physical placements. No clipping or material changes are involved.
  class BooleanFrameNormalizer {
  public:
    struct Solid {
      TGeoShape* shape;
      TGeoHMatrix frame;  // normalized coordinates -> original coordinates
    };

    Solid solid(TGeoShape* shape) {
      auto found = solids_.find(shape);
      if (found != solids_.end())
        return found->second;
      Solid result{shape, TGeoHMatrix()};
      if (auto composite = dynamic_cast<TGeoCompositeShape*>(shape)) {
        auto node = composite->GetBoolNode();
        auto left = solid(node->GetLeftShape());
        auto right = solid(node->GetRightShape());
        TGeoHMatrix leftFrame(*node->GetLeftMatrix());
        leftFrame.Multiply(&left.frame);
        TGeoHMatrix rightFrame(*node->GetRightMatrix());
        rightFrame.Multiply(&right.frame);
        TGeoHMatrix relative(leftFrame.Inverse());
        relative.Multiply(&rightFrame);
        if (!leftFrame.IsIdentity() || left.shape != node->GetLeftShape() || right.shape != node->GetRightShape() ||
            !right.frame.IsIdentity()) {
          auto matrix = new TGeoHMatrix(relative);
          TGeoBoolNode* normalized = nullptr;
          switch (node->GetBooleanOperator()) {
            case TGeoBoolNode::kGeoUnion:
              normalized = new TGeoUnion(left.shape, right.shape, nullptr, matrix);
              break;
            case TGeoBoolNode::kGeoIntersection:
              normalized = new TGeoIntersection(left.shape, right.shape, nullptr, matrix);
              break;
            case TGeoBoolNode::kGeoSubtraction:
              normalized = new TGeoSubtraction(left.shape, right.shape, nullptr, matrix);
              break;
            default:
              throw std::runtime_error("Unsupported external Boolean operation");
          }
          auto name = std::string(shape->GetName()) + "_ddg4_frame";
          result = {new TGeoCompositeShape(name.c_str(), normalized), leftFrame};
          ++changed_;
        }
      }
      solids_.emplace(shape, result);
      return result;
    }

    TGeoHMatrix volume(TGeoVolume* value) {
      auto found = volumes_.find(value);
      if (found != volumes_.end())
        return found->second;
      Solid normalized{value->GetShape(), TGeoHMatrix()};
      if (!value->IsAssembly())
        normalized = solid(value->GetShape());
      volumes_.emplace(value, normalized.frame);
      value->SetShape(normalized.shape);
      for (int index = 0; index < value->GetNdaughters(); ++index) {
        auto node = dynamic_cast<TGeoNodeMatrix*>(value->GetNode(index));
        if (!node)
          throw std::runtime_error("External geometry has a non-matrix placement");
        auto childFrame = volume(node->GetVolume());
        TGeoHMatrix placement(normalized.frame.Inverse());
        placement.Multiply(node->GetMatrix());
        placement.Multiply(&childFrame);
        node->SetMatrix(new TGeoHMatrix(placement));
      }
      if (value->IsAssembly())
        value->GetShape()->ComputeBBox();
      if (value->GetNdaughters())
        value->Voxelize("");
      return normalized.frame;
    }

    unsigned int changed() const { return changed_; }

  private:
    std::map<TGeoShape*, Solid> solids_;
    std::map<TGeoVolume*, TGeoHMatrix> volumes_;
    unsigned int changed_ = 0;
  };
}  // namespace shift
#endif
