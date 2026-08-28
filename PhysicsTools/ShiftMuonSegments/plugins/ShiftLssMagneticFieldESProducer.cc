#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {
  using Vector3 = std::array<double, 3>;
  using Matrix3 = std::array<double, 9>;

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

  struct FieldElement {
    enum class Type { uniform, quadrupole };

    explicit FieldElement(edm::ParameterSet const& parameters)
        : name(parameters.getParameter<std::string>("name")),
          minimum(vector3(parameters.getParameter<std::vector<double>>("minimumCm"), name + ".minimumCm")),
          maximum(vector3(parameters.getParameter<std::vector<double>>("maximumCm"), name + ".maximumCm")),
          origin(vector3(parameters.getParameter<std::vector<double>>("originCm"), name + ".originCm")),
          localToGlobal(
              matrix3(parameters.getParameter<std::vector<double>>("localToGlobal"), name + ".localToGlobal")),
          uniformField(vector3(parameters.getParameter<std::vector<double>>("fieldTesla"), name + ".fieldTesla")),
          gradientTeslaPerCm(parameters.getParameter<double>("gradientTeslaPerCm")) {
      auto const typeName = parameters.getParameter<std::string>("type");
      if (typeName == "uniform")
        type = Type::uniform;
      else if (typeName == "quadrupole")
        type = Type::quadrupole;
      else
        throw cms::Exception("Configuration") << name << ".type must be 'uniform' or 'quadrupole'";
      for (unsigned int axis = 0; axis < 3; ++axis) {
        if (!std::isfinite(minimum[axis]) || !std::isfinite(maximum[axis]) || !std::isfinite(origin[axis]) ||
            !(minimum[axis] < maximum[axis]))
          throw cms::Exception("Configuration") << name << " has invalid bounds or origin";
      }
      for (double value : uniformField)
        if (!std::isfinite(value))
          throw cms::Exception("Configuration") << name << ".fieldTesla must be finite";
      if (!std::isfinite(gradientTeslaPerCm))
        throw cms::Exception("Configuration") << name << ".gradientTeslaPerCm must be finite";
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
      for (unsigned int axis = 0; axis < 3; ++axis)
        if (point[axis] < minimum[axis] || point[axis] > maximum[axis])
          return false;
      return true;
    }

    GlobalVector field(Vector3 const& point) const {
      Vector3 localField = uniformField;
      if (type == Type::quadrupole)
        localField = {gradientTeslaPerCm * point[1], gradientTeslaPerCm * point[0], 0.};
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
    Matrix3 localToGlobal;
    Vector3 uniformField;
    double gradientTeslaPerCm;
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
    element.add<std::vector<double>>("localToGlobal", {1., 0., 0., 0., 1., 0., 0., 0., 1.});
    element.add<std::vector<double>>("fieldTesla", {0., 0., 0.});
    element.add<double>("gradientTeslaPerCm", 0.);

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
