#ifndef PhysicsTools_ShiftMuonSegments_MultiscaleDerivativeAudit_h
#define PhysicsTools_ShiftMuonSegments_MultiscaleDerivativeAudit_h

#include <Eigen/Cholesky>
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace shift {
  using DerivativeAuditVector = Eigen::Matrix<double, 5, 1>;
  using DerivativeAuditMatrix = Eigen::Matrix<double, 5, 5>;
  struct DerivativeAuditProbe {
    bool valid = false;
    DerivativeAuditVector parameters = DerivativeAuditVector::Zero();
    double offset = 0.;  // ACTUAL input coordinate minus the common base coordinate
  };
  struct DerivativeAuditResult {
    enum Status {
      passed = 1, invalidInput = -1, stableMismatch = -2,
      unresolvedReference = -3, invalidComparison = -4
    };
    enum ProbeStatus {
      usable = 1, invalidProbe = -1, wrongSign = -2,
      collapsedOffsets = -3, invalidMetric = -4, noReferenceInformation = -5
    };
    bool valid = false;
    bool resolved = false;
    int status = invalidInput;
    unsigned int calls = 0;
    double chosenStep = 0.;
    double stability = -1.;
    double disagreement = -1.;
    DerivativeAuditVector reference = DerivativeAuditVector::Zero();
    struct Attempt {
      double scale;
      double step;
      std::array<double, 4> actualOffsets{{0., 0., 0., 0.}};
      int probeStatus = invalidProbe;
      double stability = -1.;
      double disagreement = -1.;
      bool stable = false;
    };
    std::vector<Attempt> attempts;
  };

  // Diagnostic only. Select the FIRST independently stable numerical
  // reference while scanning from the largest to the smallest scale; never
  // search for a step which agrees with the analytic J. Starting coarse
  // avoids a false plateau made by float-quantized output coordinates. The
  // fixed order is reproducible and bounds calls to 8*4=32. References use
  // ACTUAL representable coordinate separations. Instability, collapse and
  // zero derivative information remain explicit.
  template <class Probe>
  DerivativeAuditResult auditDerivativeMultiscale(
      DerivativeAuditVector const& analytic,
      DerivativeAuditMatrix const& weight,
      Probe const& probe,
      double initialStep) {
    constexpr std::array<double, 8> scales{{16., 8., 4., 2., 1., .5, .25, .125}};
    constexpr double stabilityTolerance = .002;
    constexpr double agreementTolerance = .01;
    DerivativeAuditResult result;
    if (!analytic.allFinite() || !weight.allFinite() || !(initialStep>0.) || !std::isfinite(initialStep) ||
        (weight-weight.transpose()).norm() > 1.e-10*std::max(weight.norm(),1.e-30) ||
        weight.llt().info()!=Eigen::Success)
      return result;
    for (double scale : scales) {
      double const h=initialStep*scale;
      result.attempts.push_back({scale,h});
      auto& attempt=result.attempts.back();
      if (!(h>0.) || !std::isfinite(h)) continue;
      std::array<double,4> const requested{{h,-h,.5*h,-.5*h}};
      std::array<DerivativeAuditProbe,4> p;
      bool usable=true;
      for (unsigned int i=0;i<4;++i) {
        ++result.calls;
        p[i]=probe(requested[i]);
        attempt.actualOffsets[i]=p[i].offset;
        if (!p[i].valid || !p[i].parameters.allFinite() || !std::isfinite(p[i].offset)) {
          attempt.probeStatus=DerivativeAuditResult::invalidProbe;usable=false;break;
        }
        if (!(p[i].offset*requested[i]>0.)) {
          attempt.probeStatus=p[i].offset==0. ? DerivativeAuditResult::collapsedOffsets : DerivativeAuditResult::wrongSign;
          usable=false;break;
        }
      }
      if (!usable) continue;
      if (!(p[0].offset>p[2].offset && p[2].offset>0. &&
            p[1].offset<p[3].offset && p[3].offset<0.)) {
        attempt.probeStatus=DerivativeAuditResult::collapsedOffsets;continue;
      }
      DerivativeAuditVector const coarse=(p[0].parameters-p[1].parameters)/(p[0].offset-p[1].offset);
      DerivativeAuditVector const fine=(p[2].parameters-p[3].parameters)/(p[2].offset-p[3].offset);
      DerivativeAuditVector const difference=coarse-fine;
      double const information=fine.dot(weight*fine);
      double const error=difference.dot(weight*difference);
      if (!coarse.allFinite() || !fine.allFinite() || !std::isfinite(information) ||
          !std::isfinite(error) || information<0. || error<0.) {
        attempt.probeStatus=DerivativeAuditResult::invalidMetric;continue;
      }
      if (information==0.) {
        attempt.probeStatus=DerivativeAuditResult::noReferenceInformation;continue;
      }
      attempt.probeStatus=DerivativeAuditResult::usable;
      attempt.stability=std::sqrt(error/information);
      attempt.stable=std::isfinite(attempt.stability) && attempt.stability<=stabilityTolerance;
      if (!attempt.stable) continue;

      // From this point there is NO retry, regardless of analytic agreement.
      result.resolved=true;
      result.chosenStep=h;
      result.stability=attempt.stability;
      result.reference=fine;
      DerivativeAuditVector const disagreement=analytic-fine;
      double const mismatch=disagreement.dot(weight*disagreement);
      if (!std::isfinite(mismatch) || mismatch<0.) {
        result.status=DerivativeAuditResult::invalidComparison;return result;
      }
      attempt.disagreement=std::sqrt(mismatch/information);
      result.disagreement=attempt.disagreement;
      if (!std::isfinite(result.disagreement)) {
        result.status=DerivativeAuditResult::invalidComparison;return result;
      }
      result.valid=result.disagreement<=agreementTolerance;
      result.status=result.valid ? DerivativeAuditResult::passed : DerivativeAuditResult::stableMismatch;
      return result;
    }
    result.status=DerivativeAuditResult::unresolvedReference;
    return result;
  }
}  // namespace shift
#endif
