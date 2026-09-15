#ifndef PhysicsTools_ShiftMuonSegments_TargetMomentFit_h
#define PhysicsTools_ShiftMuonSegments_TargetMomentFit_h
#include "PhysicsTools/ShiftMuonSegments/interface/ForwardTargetFit.h"

namespace shift {
  // Quasi-score estimator: solve J^T V^-1(y-f)=0, with one independent
  // target measurement. V is re-evaluated on each new material path, but
  // neither its determinant nor its derivatives are observations. This
  // uses a conditional mean and covariance, not a claim of Gaussian tails.
  // Its covariance is local; nonlinear coverage still requires validation.
  // The transported detector posterior initializes the search only.
  template<class Transport>
  ForwardTargetFit fitTargetMoments(TargetParameters const& seed, TargetParameters const& detector,
      TargetCovariance const& detectorCovariance, double sigmaX, double sigmaY, double sigmaZ,
      Transport const& transport, unsigned int maximumIterations=32) {
    ForwardTargetFit result;
    if (!seed.allFinite() || !detector.allFinite() || !detectorCovariance.allFinite() ||
        !(sigmaX>0.) || !(sigmaY>0.) || !(sigmaZ>=0.) || !std::isfinite(sigmaX+sigmaY+sigmaZ) ||
        detectorCovariance.llt().info()!=Eigen::Success) return result;
    TargetParameters state=seed;
    for (unsigned int iteration=0;iteration<maximumIterations;++iteration) {
      result.iterations=iteration+1;result.status=-2;
      auto const prediction=transport(state);
      if (!prediction.valid || !prediction.parameters.allFinite() || !prediction.noise.allFinite()) return result;
      Eigen::LLT<TargetCovariance> solve(detectorCovariance+prediction.noise);
      Eigen::Matrix2d hit=Eigen::Matrix2d::Zero();hit(0,0)=sigmaX*sigmaX;hit(1,1)=sigmaY*sigmaY;
      hit+=sigmaZ*sigmaZ*state.template segment<2>(1)*state.template segment<2>(1).transpose();
      Eigen::LLT<Eigen::Matrix2d> hitSolve(hit);
      result.status=-3;
      if (solve.info()!=Eigen::Success || hitSolve.info()!=Eigen::Success) return result;
      TargetCovariance const weight=solve.solve(TargetCovariance::Identity());
      Eigen::Matrix2d const hitWeight=hitSolve.solve(Eigen::Matrix2d::Identity());
      TargetCovariance jacobian=prediction.jacobian;
      if (prediction.jacobianValid && !jacobian.allFinite()) {result.status=-10;return result;}
      if (!prediction.jacobianValid) for (unsigned int j=0;j<5;++j) {
        double h=j==0 ? std::abs(state[0])*1.e-3 : (j<3 ? 1.e-5 : .01);
        bool resolved=false;
        for (unsigned int refinement=0;refinement<8;++refinement) {
          auto plus=state,minus=state;plus[j]+=h;minus[j]-=h;
          auto const a=transport(plus),b=transport(minus);
          plus=state;minus=state;plus[j]+=.5*h;minus[j]-=.5*h;
          auto const c=transport(plus),d=transport(minus);
          double relativeError=-1.;
          bool const valid=a.valid && b.valid && c.valid && d.valid;
          if (a.valid && b.valid && c.valid && d.valid && h>0.) {
            TargetParameters const coarse=(a.parameters-b.parameters)/(2.*h);
            TargetParameters const fine=(c.parameters-d.parameters)/h;
            double const information=fine.dot(weight*fine)+(j>=3 ? hitWeight(j-3,j-3) : 0.);
            TargetParameters const difference=fine-coarse;
            if (information>0.) relativeError=std::sqrt(std::max(0.,difference.dot(weight*difference)/information));
            result.derivativeChecks.push_back({state,j,h,relativeError,valid});
            if (fine.allFinite() && information>0. && difference.dot(weight*difference)<1.e-4*information) {
              jacobian.col(j)=fine;resolved=true;break;
            }
          } else result.derivativeChecks.push_back({state,j,h,relativeError,valid});
          h*=.5;
        }
        if (!resolved) {result.status=-10;return result;}
      }
      TargetCovariance normal=jacobian.transpose()*weight*jacobian;
      normal.template bottomRightCorner<2,2>()+=hitWeight;
      TargetParameters const scales=normal.diagonal().array().sqrt().inverse();
      result.status=-5;
      if (!scales.allFinite()) return result;
      Eigen::LLT<TargetCovariance> normalSolve(scales.asDiagonal()*normal*scales.asDiagonal());
      if (normalSolve.info()!=Eigen::Success) return result;
      TargetCovariance const covariance=scales.asDiagonal()*normalSolve.solve(TargetCovariance::Identity())*scales.asDiagonal();
      TargetParameters const residual=detector-prediction.parameters;
      TargetParameters rhs=jacobian.transpose()*weight*residual;
      rhs.template tail<2>()-=hitWeight*state.template tail<2>();
      TargetParameters const delta=covariance*rhs;
      double const change=(delta.array().abs()/covariance.diagonal().array().sqrt()).maxCoeff();
      auto frozenChi2=[&](TargetParameters const& x, ForwardTargetPrediction const& p) {
        TargetParameters const r=detector-p.parameters;
        return r.dot(weight*r)+x.template tail<2>().dot(hitWeight*x.template tail<2>());
      };
      double const chi2=frozenChi2(state,prediction);
      if (!delta.allFinite() || !std::isfinite(change) || !std::isfinite(chi2)) return result;
      result.steps.push_back({state,chi2,change,chi2,0.});
      if (change<.02) {
        result.valid=true;result.status=1;result.parameters=state;result.covariance=covariance;result.chi2=chi2;return result;
      }
      bool accepted=false;
      // Start with the full Newton step. A cap measured in posterior sigmas
      // can prevent even an exactly linear fit from reaching its solution
      // from a distant seed. Backtracking handles nonlinearity and charge.
      for (double fraction=1.;fraction>=1./4096.;fraction*=.5) {
        TargetParameters const trial=state+fraction*delta;
        if (trial[0]*seed[0]<=0.) continue;
        auto const next=transport(trial);
        if (next.valid && next.parameters.allFinite() &&
            frozenChi2(trial,next)<chi2-2.e-4*fraction*rhs.dot(delta)) {
          state=trial;accepted=true;break;
        }
      }
      if (!accepted) {result.status=-6;return result;}
    }
    result.status=-7;return result;
  }
}
#endif
