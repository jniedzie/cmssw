#ifndef PhysicsTools_ShiftMuonSegments_ShortFlowDerivative_h
#define PhysicsTools_ShiftMuonSegments_ShortFlowDerivative_h

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>

namespace shift {
  using FlowVector = Eigen::Matrix<double, 5, 1>;
  using FlowMatrix = Eigen::Matrix<double, 5, 5>;
  struct FlowProbe {
    bool valid = false;
    double offset = 0.;  // ACTUAL destination-plane z minus actual starting-plane z
    FlowVector parameters = FlowVector::Zero();
  };
  struct FlowDerivative {
    bool valid = false;
    FlowVector derivative = FlowVector::Zero();
    double centralRelativeError = -1.;
    double oneSidedRelativeError = -1.;
    unsigned int calls = 0;
  };

  // J is the copied full source-to-detector mean Jacobian. W is the weight
  // of that SAME prediction: (C_detector+Q_long)^-1. The callback propagates
  // the SAME original source state to the requested SHORT signed delta-z;
  // it must return the actual representable plane offset. No process-noise
  // covariance is differentiated. The caller may refine h, but must enforce
  // its representability/precision floor and record unresolved checks.
  template<class Probe>
  FlowDerivative shortFlowDerivative(FlowMatrix const& jacobian,
                                      FlowMatrix const& weight,
                                      double h, double relativeTolerance,
                                      Probe const& probe) {
    FlowDerivative result;
    if (!jacobian.allFinite() || !weight.allFinite() || !(h>0.) || !std::isfinite(h) ||
        !(relativeTolerance>0.) || !std::isfinite(relativeTolerance)) return result;
    std::array<FlowProbe,4> p;
    std::array<double,4> const requested{{h,-h,.5*h,-.5*h}};
    for (unsigned int i=0;i<4;++i) {
      ++result.calls;
      p[i]=probe(requested[i]);
      if (!p[i].valid || !p[i].parameters.allFinite() || !std::isfinite(p[i].offset) ||
          !(p[i].offset*requested[i]>0.)) return result;
    }
    // Besides the two central differences, compare derivatives on either
    // side of the origin. An exact material step has identical coarse/fine
    // CENTRAL averages even though its derivative at the origin is undefined.
    if (!(p[0].offset>p[2].offset && p[2].offset>0. &&
          p[1].offset<p[3].offset && p[3].offset<0.)) return result;
    FlowVector const coarse=(p[0].parameters-p[1].parameters)/(p[0].offset-p[1].offset);
    FlowVector const fine=(p[2].parameters-p[3].parameters)/(p[2].offset-p[3].offset);
    FlowVector const plus=(p[0].parameters-p[2].parameters)/(p[0].offset-p[2].offset);
    FlowVector const minus=(p[3].parameters-p[1].parameters)/(p[3].offset-p[1].offset);
    result.derivative=-jacobian*fine;
    double const information=result.derivative.dot(weight*result.derivative);
    if (!result.derivative.allFinite() || !std::isfinite(information) || information<0.) return result;
    auto relativeError=[&](FlowVector const& difference) {
      FlowVector const transported=jacobian*difference;
      double const error=transported.dot(weight*transported);
      if (!std::isfinite(error) || error<0.) return -1.;
      if (information==0.) return error==0. ? 0. : -1.;
      return std::sqrt(error/information);
    };
    result.centralRelativeError=relativeError(coarse-fine);
    double const plusError=relativeError(plus-fine),minusError=relativeError(minus-fine);
    if (plusError<0. || minusError<0.) return result;
    result.oneSidedRelativeError=std::max(plusError,minusError);
    result.valid=result.centralRelativeError>=0. && result.centralRelativeError<=relativeTolerance &&
                 result.oneSidedRelativeError<=relativeTolerance;
    return result;
  }
}
#endif
