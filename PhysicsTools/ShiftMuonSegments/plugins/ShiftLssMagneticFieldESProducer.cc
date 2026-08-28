#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {
  using Vector3 = std::array<double, 3>;
  using Matrix3 = std::array<double, 9>;

  std::string trim(std::string value) {
    auto const first =
        std::find_if_not(value.begin(), value.end(), [](unsigned char character) { return std::isspace(character); });
    auto const last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char character) {
                        return std::isspace(character);
                      }).base();
    return first < last ? std::string(first, last) : std::string();
  }

  Vector3 vector3(std::vector<double> const& values, std::string const& parameter) {
    if (values.size() != 3)
      throw cms::Exception("Configuration") << parameter << " must contain exactly 3 values";
    return {values[0], values[1], values[2]};
  }

  Matrix3 matrix3(std::vector<double> const& values, std::string const& parameter) {
    if (values.size() != 9)
      throw cms::Exception("Configuration") << parameter << " must contain exactly 9 row-major values";
    Matrix3 result;
    std::copy(values.begin(), values.end(), result.begin());
    for (double value : result)
      if (!std::isfinite(value))
        throw cms::Exception("Configuration") << parameter << " must contain only finite values";
    for (unsigned int row = 0; row < 3; ++row) {
      for (unsigned int other = 0; other < 3; ++other) {
        double dot = 0.;
        for (unsigned int column = 0; column < 3; ++column)
          dot += result[3 * row + column] * result[3 * other + column];
        double const expected = row == other ? 1. : 0.;
        if (std::abs(dot - expected) > 1.e-9)
          throw cms::Exception("Configuration") << parameter << " must be an orthonormal rotation matrix";
      }
    }
    double const determinant = result[0] * (result[4] * result[8] - result[5] * result[7]) -
                               result[1] * (result[3] * result[8] - result[5] * result[6]) +
                               result[2] * (result[3] * result[7] - result[4] * result[6]);
    if (std::abs(determinant - 1.) > 1.e-9)
      throw cms::Exception("Configuration") << parameter << " must be a proper rotation with determinant +1";
    return result;
  }

  class FlukaFieldMap2D {
  public:
    enum class Type { quadrupole, quadrupoleInterpolated, interpolated, kickInterpolated };

    explicit FlukaFieldMap2D(std::string const& configuredPath) {
      if (configuredPath.empty())
        throw cms::Exception("Configuration") << "flukaMap2D requires a non-empty mapFile";
      std::string const path =
          configuredPath.front() == '/' ? configuredPath : edm::FileInPath(configuredPath).fullPath();
      std::ifstream input(path);
      if (!input)
        throw cms::Exception("FileOpenError") << "Cannot open FLUKA field map " << path;

      bool readingData = false;
      bool haveType = false;
      bool haveXGrid = false;
      bool haveYGrid = false;
      std::string line;
      unsigned int lineNumber = 0;
      while (std::getline(input, line)) {
        ++lineNumber;
        auto const comment = line.find('#');
        if (comment != std::string::npos)
          line.erase(comment);
        line = trim(line);
        if (line.empty())
          continue;
        std::istringstream values(line);
        if (readingData) {
          double bx = 0.;
          double by = 0.;
          std::string extra;
          if (!(values >> bx >> by) || values >> extra || !std::isfinite(bx) || !std::isfinite(by))
            throw cms::Exception("Configuration")
                << path << ':' << lineNumber << ": expected exactly two finite field components";
          bx_.push_back(bx);
          by_.push_back(by);
          continue;
        }

        std::string card;
        values >> card;
        if (card == "TYPE") {
          std::string type;
          values >> type;
          if (type == "QUAD")
            type_ = Type::quadrupole;
          else if (type == "QUADINT")
            type_ = Type::quadrupoleInterpolated;
          else if (type == "INTER2D")
            type_ = Type::interpolated;
          else if (type == "KICKINT")
            type_ = Type::kickInterpolated;
          else
            throw cms::Exception("Configuration") << path << ':' << lineNumber << ": unsupported TYPE " << type;
          haveType = true;
        } else if (card == "SYMMETRY") {
          values >> symmetry_;
          if (symmetry_ != "NONE" && symmetry_.find_first_not_of("XYZ") != std::string::npos)
            throw cms::Exception("Configuration") << path << ':' << lineNumber << ": invalid SYMMETRY";
        } else if (card == "QORIGIN") {
          if (!(values >> quadrupoleOrigin_[0] >> quadrupoleOrigin_[1] >> quadrupoleOrigin_[2]))
            throw cms::Exception("Configuration") << path << ':' << lineNumber << ": invalid QORIGIN";
        } else if (card == "QRADIUS") {
          if (!(values >> quadrupoleRadiusCm_))
            throw cms::Exception("Configuration") << path << ':' << lineNumber << ": invalid QRADIUS";
        } else if (card == "XGRID") {
          if (!(values >> minimumXcm_ >> maximumXcm_ >> binsX_))
            throw cms::Exception("Configuration") << path << ':' << lineNumber << ": invalid XGRID";
          haveXGrid = true;
        } else if (card == "YGRID") {
          if (!(values >> minimumYcm_ >> maximumYcm_ >> binsY_))
            throw cms::Exception("Configuration") << path << ':' << lineNumber << ": invalid YGRID";
          haveYGrid = true;
        } else if (card == "DATA") {
          readingData = true;
        } else {
          throw cms::Exception("Configuration") << path << ':' << lineNumber << ": unknown card " << card;
        }
      }

      if (!haveType)
        throw cms::Exception("Configuration") << path << ": missing TYPE";
      if (!(quadrupoleRadiusCm_ >= 0.) || !std::isfinite(quadrupoleRadiusCm_))
        throw cms::Exception("Configuration") << path << ": invalid QRADIUS";
      bool const needsGrid = type_ != Type::quadrupole;
      if (needsGrid && (!haveXGrid || !haveYGrid || binsX_ < 2 || binsY_ < 2 || !(minimumXcm_ < maximumXcm_) ||
                        !(minimumYcm_ < maximumYcm_)))
        throw cms::Exception("Configuration") << path << ": incomplete or invalid interpolation grid";
      std::size_t const expectedValues = needsGrid ? binsX_ * binsY_ : 0;
      if (bx_.size() != expectedValues || by_.size() != expectedValues)
        throw cms::Exception("Configuration")
            << path << ": expected " << expectedValues << " DATA rows, found " << bx_.size();
      stepXcm_ = needsGrid ? (maximumXcm_ - minimumXcm_) / static_cast<double>(binsX_ - 1) : 0.;
      stepYcm_ = needsGrid ? (maximumYcm_ - minimumYcm_) / static_cast<double>(binsY_ - 1) : 0.;
    }

    Vector3 field(Vector3 point, double scaleTesla) const {
      bool const mirroredX = point[0] < 0. && symmetry_.find('X') != std::string::npos;
      bool const mirroredY = point[1] < 0. && symmetry_.find('Y') != std::string::npos;
      if (mirroredX)
        point[0] = -point[0];
      if (mirroredY)
        point[1] = -point[1];

      double bx = 0.;
      double by = 0.;
      double const qx = point[0] - quadrupoleOrigin_[0];
      double const qy = point[1] - quadrupoleOrigin_[1];
      double const radiusSquared = qx * qx + qy * qy;
      if (type_ == Type::quadrupole ||
          (type_ == Type::quadrupoleInterpolated && radiusSquared <= quadrupoleRadiusCm_ * quadrupoleRadiusCm_)) {
        if (radiusSquared <= quadrupoleRadiusCm_ * quadrupoleRadiusCm_) {
          bx = qy * scaleTesla;
          by = qx * scaleTesla;
        }
      } else if (type_ == Type::kickInterpolated &&
                 point[0] * point[0] + point[1] * point[1] <= quadrupoleRadiusCm_ * quadrupoleRadiusCm_) {
        by = scaleTesla;
      } else {
        auto const interpolated = interpolate(point[0], point[1]);
        bx = interpolated[0] * scaleTesla;
        by = interpolated[1] * scaleTesla;
      }

      if (std::hypot(bx, by) < 1.e-8)
        return {};
      if (mirroredX)
        by = -by;
      if (mirroredY)
        bx = -bx;
      return {bx, by, 0.};
    }

  private:
    Vector3 interpolate(double x, double y) const {
      double const gridX = (x - minimumXcm_) / stepXcm_;
      double const gridY = (y - minimumYcm_) / stepYcm_;
      if (gridX < 0. || gridY < 0.)
        return {};
      auto const indexX = static_cast<std::size_t>(gridX);
      auto const indexY = static_cast<std::size_t>(gridY);
      if (indexX >= binsX_ - 1 || indexY >= binsY_ - 1)
        return {};
      double const fractionX = gridX - static_cast<double>(indexX);
      double const fractionY = gridY - static_cast<double>(indexY);
      auto bilinear = [&](std::vector<double> const& values) {
        auto const lower = indexY * binsX_ + indexX;
        auto const upper = lower + binsX_;
        double const atLowerY = values[lower] + fractionX * (values[lower + 1] - values[lower]);
        double const atUpperY = values[upper] + fractionX * (values[upper + 1] - values[upper]);
        return atLowerY + fractionY * (atUpperY - atLowerY);
      };
      return {bilinear(bx_), bilinear(by_), 0.};
    }

    Type type_ = Type::interpolated;
    std::string symmetry_ = "NONE";
    Vector3 quadrupoleOrigin_ = {0., 0., 0.};
    double quadrupoleRadiusCm_ = 0.;
    double minimumXcm_ = 0.;
    double maximumXcm_ = 0.;
    double minimumYcm_ = 0.;
    double maximumYcm_ = 0.;
    std::size_t binsX_ = 0;
    std::size_t binsY_ = 0;
    double stepXcm_ = 0.;
    double stepYcm_ = 0.;
    std::vector<double> bx_;
    std::vector<double> by_;
  };

  struct FieldElement {
    enum class Type { uniform, quadrupole, flukaMap2D };

    explicit FieldElement(edm::ParameterSet const& parameters)
        : name(parameters.getParameter<std::string>("name")),
          minimum(vector3(parameters.getParameter<std::vector<double>>("minimumCm"), name + ".minimumCm")),
          maximum(vector3(parameters.getParameter<std::vector<double>>("maximumCm"), name + ".maximumCm")),
          origin(vector3(parameters.getParameter<std::vector<double>>("originCm"), name + ".originCm")),
          boundsCenter(
              vector3(parameters.getParameter<std::vector<double>>("boundsCenterCm"), name + ".boundsCenterCm")),
          localToGlobal(
              matrix3(parameters.getParameter<std::vector<double>>("localToGlobal"), name + ".localToGlobal")),
          uniformField(vector3(parameters.getParameter<std::vector<double>>("fieldTesla"), name + ".fieldTesla")),
          gradientTeslaPerCm(parameters.getParameter<double>("gradientTeslaPerCm")),
          fieldScaleTesla(parameters.getParameter<double>("fieldScaleTesla")),
          innerRadiusCm(parameters.getParameter<double>("innerRadiusCm")),
          outerRadiusCm(parameters.getParameter<double>("outerRadiusCm")),
          excludedCylindersCm(parameters.getParameter<std::vector<double>>("excludedCylindersCm")),
          cylindricalBounds(parameters.getParameter<std::string>("boundsShape") == "cylinderZ") {
      auto const boundsShape = parameters.getParameter<std::string>("boundsShape");
      if (boundsShape != "box" && boundsShape != "cylinderZ")
        throw cms::Exception("Configuration") << name << ".boundsShape must be 'box' or 'cylinderZ'";
      auto const typeName = parameters.getParameter<std::string>("type");
      if (typeName == "uniform")
        type = Type::uniform;
      else if (typeName == "quadrupole")
        type = Type::quadrupole;
      else if (typeName == "flukaMap2D") {
        type = Type::flukaMap2D;
        flukaMap = std::make_shared<FlukaFieldMap2D>(parameters.getParameter<std::string>("mapFile"));
      } else
        throw cms::Exception("Configuration") << name << ".type must be 'uniform', 'quadrupole', or 'flukaMap2D'";
      for (unsigned int axis = 0; axis < 3; ++axis) {
        if (!std::isfinite(minimum[axis]) || !std::isfinite(maximum[axis]) || !std::isfinite(origin[axis]) ||
            !std::isfinite(boundsCenter[axis]) || !(minimum[axis] < maximum[axis]))
          throw cms::Exception("Configuration") << name << " has invalid bounds or origin";
      }
      if (cylindricalBounds &&
          (!(innerRadiusCm >= 0.) || !(outerRadiusCm > innerRadiusCm) || !std::isfinite(innerRadiusCm) ||
           !std::isfinite(outerRadiusCm)))
        throw cms::Exception("Configuration") << name << " has invalid cylindrical radii";
      if (excludedCylindersCm.size() % 3 != 0)
        throw cms::Exception("Configuration") << name << ".excludedCylindersCm must contain x, y, radius triples";
      for (std::size_t index = 0; index < excludedCylindersCm.size(); index += 3)
        if (!std::isfinite(excludedCylindersCm[index]) || !std::isfinite(excludedCylindersCm[index + 1]) ||
            !(excludedCylindersCm[index + 2] > 0.) || !std::isfinite(excludedCylindersCm[index + 2]))
          throw cms::Exception("Configuration") << name << " has an invalid excluded cylinder";
      for (double value : uniformField)
        if (!std::isfinite(value))
          throw cms::Exception("Configuration") << name << ".fieldTesla must be finite";
      if (!std::isfinite(gradientTeslaPerCm))
        throw cms::Exception("Configuration") << name << ".gradientTeslaPerCm must be finite";
      if (!std::isfinite(fieldScaleTesla))
        throw cms::Exception("Configuration") << name << ".fieldScaleTesla must be finite";
    }

    Vector3 localPosition(GlobalPoint const& point) const {
      Vector3 const displacement = {point.x() - origin[0], point.y() - origin[1], point.z() - origin[2]};
      Vector3 result{};
      for (unsigned int local = 0; local < 3; ++local)
        for (unsigned int global = 0; global < 3; ++global)
          result[local] += localToGlobal[3 * global + local] * displacement[global];
      return result;
    }

    bool contains(Vector3 const& point) const {
      if (!cylindricalBounds) {
        for (unsigned int axis = 0; axis < 3; ++axis)
          if (point[axis] < minimum[axis] || point[axis] > maximum[axis])
            return false;
        return true;
      }
      if (point[2] < minimum[2] || point[2] > maximum[2])
        return false;
      double const dx = point[0] - boundsCenter[0];
      double const dy = point[1] - boundsCenter[1];
      double const radiusSquared = dx * dx + dy * dy;
      if (radiusSquared < innerRadiusCm * innerRadiusCm || radiusSquared > outerRadiusCm * outerRadiusCm)
        return false;
      for (std::size_t index = 0; index < excludedCylindersCm.size(); index += 3) {
        double const excludedDx = point[0] - excludedCylindersCm[index];
        double const excludedDy = point[1] - excludedCylindersCm[index + 1];
        if (excludedDx * excludedDx + excludedDy * excludedDy <
            excludedCylindersCm[index + 2] * excludedCylindersCm[index + 2])
          return false;
      }
      return true;
    }

    GlobalVector field(Vector3 const& point) const {
      Vector3 localField = uniformField;
      if (type == Type::quadrupole)
        localField = {gradientTeslaPerCm * point[1], gradientTeslaPerCm * point[0], 0.};
      else if (type == Type::flukaMap2D)
        localField = flukaMap->field(point, fieldScaleTesla);
      Vector3 globalField{};
      for (unsigned int global = 0; global < 3; ++global)
        for (unsigned int local = 0; local < 3; ++local)
          globalField[global] += localToGlobal[3 * global + local] * localField[local];
      return GlobalVector(globalField[0], globalField[1], globalField[2]);
    }

    std::string name;
    Type type;
    Vector3 minimum;
    Vector3 maximum;
    Vector3 origin;
    Vector3 boundsCenter;
    Matrix3 localToGlobal;
    Vector3 uniformField;
    double gradientTeslaPerCm;
    double fieldScaleTesla;
    double innerRadiusCm;
    double outerRadiusCm;
    std::vector<double> excludedCylindersCm;
    bool cylindricalBounds;
    std::shared_ptr<FlukaFieldMap2D const> flukaMap;
  };

  class ShiftLssMagneticField final : public MagneticField {
  public:
    ShiftLssMagneticField(MagneticField const& baseField,
                          std::vector<FieldElement> elements,
                          bool addBaseField,
                          bool sumOverlaps)
        : baseField_(&baseField),
          elements_(std::move(elements)),
          addBaseField_(addBaseField),
          sumOverlaps_(sumOverlaps) {
      setNominalValue();
    }

    ShiftLssMagneticField* clone() const override { return new ShiftLssMagneticField(*this); }

    GlobalVector inTesla(GlobalPoint const& point) const override {
      GlobalVector result = addBaseField_ ? baseField_->inTesla(point) : GlobalVector();
      bool matched = false;
      for (auto const& element : elements_) {
        auto const local = element.localPosition(point);
        if (!element.contains(local))
          continue;
        auto const elementField = element.field(local);
        if (sumOverlaps_ || !matched)
          result += elementField;
        matched = true;
        if (!sumOverlaps_)
          break;
      }
      if (!addBaseField_ && !matched)
        return baseField_->inTesla(point);
      return result;
    }

    bool isDefined(GlobalPoint const& point) const override {
      for (auto const& element : elements_)
        if (element.contains(element.localPosition(point)))
          return true;
      return baseField_->isDefined(point);
    }

  private:
    MagneticField const* baseField_;
    std::vector<FieldElement> elements_;
    bool addBaseField_;
    bool sumOverlaps_;
  };
}  // namespace

class ShiftLssMagneticFieldESProducer : public edm::ESProducer {
public:
  explicit ShiftLssMagneticFieldESProducer(edm::ParameterSet const& parameters)
      : addBaseField_(parameters.getParameter<bool>("addBaseField")),
        sumOverlaps_(parameters.getParameter<bool>("sumOverlaps")) {
    for (auto const& element : parameters.getParameter<std::vector<edm::ParameterSet>>("elements"))
      elements_.emplace_back(element);
    auto const baseFieldLabel = parameters.getParameter<std::string>("baseFieldLabel");
    auto const outputLabel = parameters.getParameter<std::string>("outputLabel");
    if (baseFieldLabel == outputLabel)
      throw cms::Exception("Configuration") << "baseFieldLabel and outputLabel must differ";
    auto collector = setWhatProduced(this, outputLabel);
    baseFieldToken_ = collector.consumes(edm::ESInputTag("", baseFieldLabel));
  }

  std::unique_ptr<MagneticField> produce(IdealMagneticFieldRecord const& record) {
    return std::make_unique<ShiftLssMagneticField>(record.get(baseFieldToken_), elements_, addBaseField_, sumOverlaps_);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription element;
    element.add<std::string>("name");
    element.add<std::string>("type");
    element.add<std::vector<double>>("minimumCm");
    element.add<std::vector<double>>("maximumCm");
    element.add<std::vector<double>>("originCm", {0., 0., 0.});
    element.add<std::string>("boundsShape", "box");
    element.add<std::vector<double>>("boundsCenterCm", {0., 0., 0.});
    element.add<double>("innerRadiusCm", 0.);
    element.add<double>("outerRadiusCm", 0.);
    element.add<std::vector<double>>("excludedCylindersCm", {});
    element.add<std::vector<double>>("localToGlobal", {1., 0., 0., 0., 1., 0., 0., 0., 1.});
    element.add<std::vector<double>>("fieldTesla", {0., 0., 0.});
    element.add<double>("gradientTeslaPerCm", 0.);
    element.add<std::string>("mapFile", "");
    element.add<double>("fieldScaleTesla", 1.);

    edm::ParameterSetDescription description;
    description.add<std::string>("baseFieldLabel", "");
    description.add<std::string>("outputLabel", "shiftLssMagneticField");
    description.add<bool>("addBaseField", false);
    description.add<bool>("sumOverlaps", false);
    description.addVPSet("elements", element, {});
    descriptions.addWithDefaultLabel(description);
  }

private:
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> baseFieldToken_;
  std::vector<FieldElement> elements_;
  bool addBaseField_;
  bool sumOverlaps_;
};

DEFINE_FWK_EVENTSETUP_MODULE(ShiftLssMagneticFieldESProducer);
