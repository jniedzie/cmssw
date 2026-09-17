#ifndef PhysicsTools_ShiftMuonSegments_JointForwardVertexFit_h
#define PhysicsTools_ShiftMuonSegments_JointForwardVertexFit_h

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

namespace shift {
  using JointVertexTrackParameters = Eigen::Matrix<double, 5, 1>;
  using JointVertexTrackCovariance = Eigen::Matrix<double, 5, 5>;
  using JointVertexParameters = Eigen::Matrix<double, 9, 1>;
  using JointVertexCovariance = Eigen::Matrix<double, 9, 9>;
  struct JointVertexTrack {
    JointVertexTrackParameters parameters = JointVertexTrackParameters::Zero();
    JointVertexTrackCovariance covariance = JointVertexTrackCovariance::Zero();
  };
  struct JointVertexPrediction {
    bool valid = false;
    JointVertexTrackParameters parameters = JointVertexTrackParameters::Zero();
    JointVertexTrackCovariance noise = JointVertexTrackCovariance::Zero();
    bool jacobianValid = false;
    JointVertexTrackCovariance jacobian = JointVertexTrackCovariance::Zero();
    bool zDerivativeValid = false;
    JointVertexTrackParameters zDerivative = JointVertexTrackParameters::Zero();
  };
  struct JointVertexOptions {
    unsigned int maximumIterations = 32;
    unsigned int maximumTransportCalls = 600;
    unsigned int maximumDerivativeRefinements = 8;
    unsigned int maximumBacktracks = 13;
    double maximumSeconds = 0.;  // optional diagnostic guard; disabled for reproducibility
    double convergence = .02;
    double derivativeRelativeTolerance = .01;
    double zDerivativeStep = .1;  // cm, checked against a half-size step
  };
  struct JointForwardVertexFit {
    enum Status {
      success = 1, invalidInput = -1, invalidTransport = -2, invalidCovariance = -3,
      singularNormal = -5, noDescent = -6, iterationLimit = -7,
      invalidDerivative = -10, transportLimit = -11, invalidCharge = -12
    };
    bool valid = false;
    int status = invalidInput;
    // (x,y,z,q/p_1,dx/dz_1,dy/dz_1,q/p_2,dx/dz_2,dy/dz_2).
    JointVertexParameters parameters = JointVertexParameters::Zero();
    JointVertexCovariance covariance = JointVertexCovariance::Zero();
    double chi2 = -1.;
    unsigned int iterations = 0;
    unsigned int transportCalls = 0;
    struct Trial {
      double fraction;
      double chi2;  // -1 when charge or incomplete transport prevents evaluation
      bool valid;
    };
    struct Step {
      JointVertexParameters parameters;
      double chi2;
      double maxUpdate;
      double acceptedFraction;
      double acceptedFrozenChi2;
      JointVertexParameters direction = JointVertexParameters::Zero();
      double expectedReduction = 0.;  // rhs.dot(direction), half the predicted negative slope
      double jointUpdate = 0.;  // sqrt(direction^T covariance^-1 direction)
      std::array<JointVertexTrackCovariance, 2> weights;
      std::array<Eigen::Matrix<double, 5, 9>, 2> designs;
      std::array<JointVertexTrackParameters, 2> predictedParameters;
      std::vector<Trial> trials;
    };
    std::vector<Step> steps;
    struct DerivativeCheck {
      unsigned int track;
      unsigned int coordinate;  // 0..4 state, 5 vertex z
      double step;
      double relativeError;
      bool valid;
    };
    std::vector<DerivativeCheck> derivativeChecks;
  };

  // Fit ORIGINAL, independent detector measurements by forward transport
  // from a common fitted vertex. Q is process noise from ZERO input error.
  // The sole external constraint is a transverse beam-line prior. There is
  // no z observation, no backward material correction and no repeated use of
  // a fitted track posterior. At each iteration solve the moment quasi-score
  // J^T (C_detector+Q)^-1 (measurement-prediction)=0 with Q frozen in its line
  // search. Covariance derivatives and log(det V) are not observations.
  // The resulting local covariance does not by itself establish coverage.
  // Callback: transport(track index, (q/p,tx,ty,x,y), vertex z).
  template <class Transport>
  JointForwardVertexFit fitJointForwardVertex(
      JointVertexParameters const& seed,
      std::array<JointVertexTrack, 2> const& detector,
      double beamLineSigma,
      Transport const& transport,
      JointVertexOptions const& options = JointVertexOptions{}) {
    JointForwardVertexFit result;
    result.parameters = seed;
    if (!seed.allFinite() || !(beamLineSigma > 0.) || !std::isfinite(beamLineSigma) ||
        !(options.maximumSeconds >= 0.) || !std::isfinite(options.maximumSeconds) ||
        !(options.convergence > 0.) || !std::isfinite(options.convergence) ||
        !(options.derivativeRelativeTolerance > 0.) || !std::isfinite(options.derivativeRelativeTolerance) ||
        !(options.zDerivativeStep > 0.) || !std::isfinite(options.zDerivativeStep) ||
        options.maximumIterations == 0 || options.maximumTransportCalls == 0 ||
        options.maximumDerivativeRefinements == 0 || options.maximumBacktracks == 0)
      return result;
    auto symmetric = [](JointVertexTrackCovariance const& c) {
      return c.allFinite() && (c-c.transpose()).norm() <= 1.e-10*std::max(c.norm(), 1.e-30);
    };
    for (unsigned int i=0; i<2; ++i) {
      auto const& d=detector[i];
      if (!d.parameters.allFinite() || !symmetric(d.covariance) ||
          d.covariance.llt().info()!=Eigen::Success)
        return result;
      if (!(seed[3+3*i]*d.parameters[0] > 0.)) {
        result.status=JointForwardVertexFit::invalidCharge;return result;
      }
    }
    auto const started=std::chrono::steady_clock::now();
    auto invoke = [&](unsigned int track, JointVertexTrackParameters const& state, double z,
                      JointVertexPrediction& prediction) {
      if (result.transportCalls >= options.maximumTransportCalls ||
          (options.maximumSeconds>0. &&
           std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count() >= options.maximumSeconds)) {
        result.status=JointForwardVertexFit::transportLimit;return false;
      }
      ++result.transportCalls;
      prediction=transport(track,state,z);
      return true;
    };
    auto trackState = [](JointVertexParameters const& state, unsigned int i) {
      JointVertexTrackParameters value;
      value.template head<3>()=state.template segment<3>(3+3*i);
      value.template tail<2>()=state.template head<2>();
      return value;
    };
    auto validMean = [&](JointVertexPrediction const& p, unsigned int i) {
      return p.valid && p.parameters.allFinite() && p.parameters[0]*seed[3+3*i]>0.;
    };
    double const priorWeight=1./(beamLineSigma*beamLineSigma);
    if (!(priorWeight>0.) || !std::isfinite(priorWeight)) return result;
    JointVertexParameters state=seed;
    for (unsigned int iteration=0; iteration<options.maximumIterations; ++iteration) {
      result.iterations=iteration+1;
      result.parameters=state;
      std::array<JointVertexPrediction, 2> predictions;
      std::array<JointVertexTrackCovariance, 2> weights;
      std::array<Eigen::Matrix<double, 5, 9>, 2> designs;
      JointVertexCovariance normal=JointVertexCovariance::Zero();
      normal(0,0)=normal(1,1)=priorWeight;
      JointVertexParameters rhs=JointVertexParameters::Zero();
      rhs.template head<2>()=-priorWeight*state.template head<2>();
      double chi2=priorWeight*state.template head<2>().squaredNorm();
      for (unsigned int i=0; i<2; ++i) {
        JointVertexTrackParameters const local=trackState(state,i);
        auto& p=predictions[i];
        if (!invoke(i,local,state[2],p)) return result;
        if (!validMean(p,i) || !symmetric(p.noise)) {
          result.status=JointForwardVertexFit::invalidTransport;return result;
        }
        // Q must be positive semidefinite, not merely rescued by C_detector.
        // Whiten before testing to avoid mixed position/curvature units.
        Eigen::LLT<JointVertexTrackCovariance> detectorSolve(detector[i].covariance);
        JointVertexTrackCovariance const left=detectorSolve.matrixL().solve(p.noise);
        JointVertexTrackCovariance const whitened=detectorSolve.matrixL().solve(left.transpose()).transpose();
        Eigen::SelfAdjointEigenSolver<JointVertexTrackCovariance> noiseSolve(whitened, Eigen::EigenvaluesOnly);
        Eigen::LLT<JointVertexTrackCovariance> solve(detector[i].covariance+p.noise);
        if (noiseSolve.info()!=Eigen::Success ||
            noiseSolve.eigenvalues().minCoeff() < -1.e-10*std::max(1.,noiseSolve.eigenvalues().cwiseAbs().maxCoeff()) ||
            solve.info()!=Eigen::Success) {
          result.status=JointForwardVertexFit::invalidCovariance;return result;
        }
        weights[i]=solve.solve(JointVertexTrackCovariance::Identity());
        JointVertexTrackCovariance jacobian=p.jacobian;
        JointVertexTrackParameters zDerivative=p.zDerivative;
        if ((p.jacobianValid && !jacobian.allFinite()) || (p.zDerivativeValid && !zDerivative.allFinite())) {
          result.status=JointForwardVertexFit::invalidDerivative;return result;
        }
        for (unsigned int j=0; j<6; ++j) {
          if ((j<5 && p.jacobianValid) || (j==5 && p.zDerivativeValid)) continue;
          double h=j==5 ? options.zDerivativeStep :
              (j==0 ? std::abs(local[0])*1.e-3 : (j<3 ? 1.e-5 : .01));
          bool resolved=false;
          for (unsigned int refinement=0; refinement<options.maximumDerivativeRefinements; ++refinement) {
            std::array<JointVertexPrediction, 4> perturbed;
            std::array<double, 4> const offsets{{h,-h,.5*h,-.5*h}};
            bool valid=true;
            for (unsigned int k=0; k<4; ++k) {
              auto trial=local;double z=state[2];
              if (j<5) trial[j]+=offsets[k];else z+=offsets[k];
              if (!trial.allFinite() || trial[0]*local[0]<=0. || !std::isfinite(z) ||
                  (j<5 ? trial[j]==local[j] : z==state[2])) {valid=false;break;}
              if (!invoke(i,trial,z,perturbed[k])) return result;
              if (!validMean(perturbed[k],i)) {valid=false;break;}
            }
            double relativeError=-1.;
            if (valid) {
              JointVertexTrackParameters const coarse=(perturbed[0].parameters-perturbed[1].parameters)/(2.*h);
              JointVertexTrackParameters const fine=(perturbed[2].parameters-perturbed[3].parameters)/h;
              JointVertexTrackParameters const difference=fine-coarse;
              double const information=fine.dot(weights[i]*fine)+(j==3 || j==4 ? priorWeight : 0.);
              double const error=difference.dot(weights[i]*difference);
              if (fine.allFinite() && std::isfinite(information) && std::isfinite(error) && information>=0. && error>=0.) {
                relativeError=information>0. ? std::sqrt(error/information) : (error==0. ? 0. : -1.);
                if ((information>0. && error <= options.derivativeRelativeTolerance*options.derivativeRelativeTolerance*information) ||
                    (information==0. && error==0.)) {
                  if (j<5) jacobian.col(j)=fine;else zDerivative=fine;
                  resolved=true;
                }
              }
            }
            result.derivativeChecks.push_back({i,j,h,relativeError,valid});
            if (resolved) break;
            h*=.5;
          }
          if (!resolved) {result.status=JointForwardVertexFit::invalidDerivative;return result;}
        }
        auto& design=designs[i];design.setZero();
        design.col(0)=jacobian.col(3);design.col(1)=jacobian.col(4);design.col(2)=zDerivative;
        design.template block<5,3>(0,3+3*i)=jacobian.template leftCols<3>();
        JointVertexTrackParameters const residual=detector[i].parameters-p.parameters;
        normal+=design.transpose()*weights[i]*design;
        rhs+=design.transpose()*weights[i]*residual;
        chi2+=residual.dot(weights[i]*residual);
      }
      JointVertexParameters const scales=normal.diagonal().array().sqrt().inverse();
      result.status=JointForwardVertexFit::singularNormal;
      if (!scales.allFinite() || !rhs.allFinite() || !std::isfinite(chi2)) return result;
      Eigen::LLT<JointVertexCovariance> normalSolve(scales.asDiagonal()*normal*scales.asDiagonal());
      if (normalSolve.info()!=Eigen::Success || normalSolve.matrixL().toDenseMatrix().diagonal().minCoeff()<1.e-6)
        return result;
      JointVertexCovariance const covariance=scales.asDiagonal()*normalSolve.solve(JointVertexCovariance::Identity())*scales.asDiagonal();
      JointVertexParameters const delta=scales.asDiagonal()*normalSolve.solve((scales.array()*rhs.array()).matrix());
      double const change=(delta.array().abs()/covariance.diagonal().array().sqrt()).maxCoeff();
      if (!covariance.allFinite() || !delta.allFinite() || !std::isfinite(change)) return result;
      result.covariance=covariance;result.chi2=chi2;
      double const expectedReduction=rhs.dot(delta);
      if (!(expectedReduction>=0.) || !std::isfinite(expectedReduction)) return result;
      double const jointUpdate=std::sqrt(expectedReduction);
      result.steps.push_back({state,chi2,change,0.,chi2,delta,expectedReduction,jointUpdate,weights,designs,
                              {{predictions[0].parameters,predictions[1].parameters}}, {}});
      // Marginal errors alone can hide motion along a tightly constrained
      // correlated direction. Retain the existing coordinate-wise condition
      // and bound the full Mahalanobis motion at the same per-coordinate
      // tolerance in this nine-dimensional fit.
      if (change < options.convergence &&
          jointUpdate < options.convergence*std::sqrt(double(JointVertexParameters::RowsAtCompileTime))) {
        result.valid=true;result.status=JointForwardVertexFit::success;return result;
      }
      if (!(expectedReduction>0.) || !std::isfinite(expectedReduction)) return result;
      bool accepted=false;double fraction=1.;
      for (unsigned int backtrack=0; backtrack<options.maximumBacktracks; ++backtrack,fraction*=.5) {
        result.steps.back().trials.push_back({fraction,-1.,false});
        JointVertexParameters const trial=state+fraction*delta;
        if (!trial.allFinite() || trial[3]*seed[3]<=0. || trial[6]*seed[6]<=0.) continue;
        double trialChi2=priorWeight*trial.template head<2>().squaredNorm();
        bool valid=true;
        for (unsigned int i=0; i<2; ++i) {
          JointVertexPrediction p;
          if (!invoke(i,trackState(trial,i),trial[2],p)) return result;
          if (!validMean(p,i) || !symmetric(p.noise)) {valid=false;break;}
          JointVertexTrackParameters const residual=detector[i].parameters-p.parameters;
          trialChi2+=residual.dot(weights[i]*residual);
        }
        if (valid) {
          result.steps.back().trials.back().chi2=trialChi2;
          result.steps.back().trials.back().valid=std::isfinite(trialChi2);
        }
        if (valid && std::isfinite(trialChi2) && trialChi2 < chi2-2.e-4*fraction*expectedReduction) {
          state=trial;accepted=true;
          result.steps.back().acceptedFraction=fraction;
          result.steps.back().acceptedFrozenChi2=trialChi2;
          break;
        }
      }
      if (!accepted) {result.status=JointForwardVertexFit::noDescent;return result;}
    }
    result.parameters=state;
    result.status=JointForwardVertexFit::iterationLimit;return result;
  }
}  // namespace shift
#endif
