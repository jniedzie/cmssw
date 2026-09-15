#ifndef PhysicsTools_ShiftMuonSegments_ImportanceTargetConstraint_h
#define PhysicsTools_ShiftMuonSegments_ImportanceTargetConstraint_h

#include "PhysicsTools/ShiftMuonSegments/interface/CubatureTransport.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <vector>

namespace shift {
  struct ImportanceTargetResult {
    bool valid = false;
    TransportMoments prior, posterior;
    double effectivePoints = 0., maximumWeight = 0.;
  };

  // Integrate over an intermediate track state. A Gaussian approximation to
  // its target-conditioned distribution is used ONLY as an importance
  // proposal. Dividing by the proposal's target likelihood counts the target
  // measurement once. Each actual material path retains its own conditional
  // covariance until after the target update. Affine Gaussian transport is
  // reproduced exactly, including correlations and process noise.
  template<class Transport>
  ImportanceTargetResult importanceTargetConstraint(TransportMoments const& input, Transport const& transport,
                                                     double sigmaX, double sigmaY, double sigmaZ,
                                                     unsigned int clouds = 1) {
    using Vector=Eigen::Matrix<double,5,1>;
    using Matrix=Eigen::Matrix<double,5,5>;
    ImportanceTargetResult result;
    if (!input.valid || !input.mean.allFinite() || !input.covariance.allFinite() ||
        !(sigmaX>0. && sigmaY>0. && sigmaZ>=0.) || clouds<1 || clouds>4) return result;
    auto hitCovariance=[&](Vector const& mean) {
      Eigen::Matrix2d hit=Eigen::Matrix2d::Zero();
      hit(0,0)=sigmaX*sigmaX;hit(1,1)=sigmaY*sigmaY;
      hit+=sigmaZ*sigmaZ*mean.template segment<2>(1)*mean.template segment<2>(1).transpose();
      return hit;
    };
    auto points=[&](TransportMoments const& state) {
      std::vector<Vector> values;
      Eigen::LLT<Matrix> factor(state.covariance);
      if (factor.info()!=Eigen::Success) return values;
      Matrix const root=std::sqrt(5.)*factor.matrixL().toDenseMatrix();
      for (unsigned int cloud=0;cloud<clouds;++cloud) {
        Matrix rotation=Matrix::Identity();
        if (cloud) {
          Vector axis;
          for(unsigned int j=0;j<5;++j) axis[j]=1.+((j+cloud-1)%5);
          axis.normalize();rotation-=2.*axis*axis.transpose();
        }
        Matrix const spread=root*rotation;
        for(unsigned int j=0;j<5;++j) for(double sign:{-1.,1.}) values.push_back(state.mean+sign*spread.col(j));
      }
      return values;
    };
    auto const pilotPoints=points(input);
    if (pilotPoints.empty()) return result;
    std::vector<TransportMoments> pilot;
    Vector mean=Vector::Zero();
    for(auto const& point:pilotPoints) {
      auto prediction=transport(point);
      if(!prediction.valid || !prediction.mean.allFinite() || !prediction.covariance.allFinite()) return result;
      mean+=prediction.mean/pilotPoints.size();pilot.push_back(prediction);
    }
    Matrix covariance=Matrix::Zero(),cross=Matrix::Zero();
    for(unsigned int i=0;i<pilot.size();++i) {
      Vector const delta=pilot[i].mean-mean;
      covariance+=(pilot[i].covariance+delta*delta.transpose())/pilot.size();
      cross+=delta*(pilotPoints[i]-input.mean).transpose()/pilot.size();
    }
    Eigen::LLT<Matrix> inputSolve(input.covariance);
    Matrix const regression=inputSolve.solve(cross.transpose()).transpose();
    Matrix conditionalNoise=covariance-regression*input.covariance*regression.transpose();
    conditionalNoise=.5*(conditionalNoise+conditionalNoise.transpose()).eval();
    Eigen::Matrix2d const referenceHit=hitCovariance(mean);
    Eigen::Matrix2d const innovation=covariance.template bottomRightCorner<2,2>()+referenceHit;
    Eigen::LLT<Eigen::Matrix2d> innovationSolve(innovation);
    if(innovationSolve.info()!=Eigen::Success) return result;
    Eigen::Matrix<double,5,2> const gain=innovationSolve.solve(cross.template bottomRows<2>()).transpose();
    TransportMoments proposal{true,input.mean-gain*mean.template tail<2>(),
                             input.covariance-gain*cross.template bottomRows<2>()};
    proposal.covariance=.5*(proposal.covariance+proposal.covariance.transpose()).eval();
    auto const posteriorPoints=points(proposal);
    if(posteriorPoints.empty()) return result;
    Eigen::LLT<Eigen::Matrix2d> referenceSolve(conditionalNoise.template bottomRightCorner<2,2>()+referenceHit);
    if(referenceSolve.info()!=Eigen::Success) return result;
    auto logLikelihood=[](Eigen::Vector2d const& residual,Eigen::LLT<Eigen::Matrix2d> const& solve) {
      return -.5*solve.matrixL().solve(residual).squaredNorm()-solve.matrixL().toDenseMatrix().diagonal().array().log().sum();
    };
    std::vector<TransportMoments> conditioned;
    std::vector<double> logs;
    for(auto const& point:posteriorPoints) {
      auto prediction=transport(point);
      if(!prediction.valid || !prediction.mean.allFinite() || !prediction.covariance.allFinite()) return result;
      Eigen::Matrix2d const hit=hitCovariance(prediction.mean);
      Eigen::LLT<Eigen::Matrix2d> solve(prediction.covariance.template bottomRightCorner<2,2>()+hit);
      if(solve.info()!=Eigen::Success) return result;
      Vector const reference=mean+regression*(point-input.mean);
      double const log=logLikelihood(prediction.mean.template tail<2>(),solve)-
                       logLikelihood(reference.template tail<2>(),referenceSolve);
      if(!std::isfinite(log)) return result;
      logs.push_back(log);
      Eigen::Matrix<double,5,2> const update=solve.solve(prediction.covariance.template bottomRows<2>()).transpose();
      Matrix transform=Matrix::Identity();transform.template rightCols<2>()-=update;
      prediction.mean-=update*prediction.mean.template tail<2>().eval();
      if(prediction.mean[0]*input.mean[0]<=0.) return result;
      prediction.covariance=(transform*prediction.covariance*transform.transpose()+update*hit*update.transpose()).eval();
      conditioned.push_back(prediction);
    }
    double const maximum=*std::max_element(logs.begin(),logs.end());
    double total=0.;for(double& log:logs) {log=std::exp(log-maximum);total+=log;}
    Vector posteriorMean=Vector::Zero();
    for(unsigned int i=0;i<logs.size();++i) {logs[i]/=total;posteriorMean+=logs[i]*conditioned[i].mean;}
    Matrix posteriorCovariance=Matrix::Zero();double squareWeights=0.;
    for(unsigned int i=0;i<logs.size();++i) {
      Vector const delta=conditioned[i].mean-posteriorMean;
      posteriorCovariance+=logs[i]*(conditioned[i].covariance+delta*delta.transpose());
      squareWeights+=logs[i]*logs[i];
    }
    result.prior={true,mean,covariance};result.posterior={true,posteriorMean,posteriorCovariance};
    result.effectivePoints=1./squareWeights;result.maximumWeight=*std::max_element(logs.begin(),logs.end());
    result.valid=posteriorMean.allFinite() && posteriorCovariance.allFinite() && posteriorCovariance.llt().info()==Eigen::Success;
    return result;
  }
}
#endif
