#include "PhysicsTools/ShiftMuonSegments/interface/JointForwardVertexFit.h"
#include "PhysicsTools/ShiftMuonSegments/interface/CommonVertexRefit.h"
#include <Eigen/LU>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
  using namespace shift;
  void require(bool value, std::string const& message) {
    if (!value) throw std::runtime_error(message);
  }
  JointVertexTrackParameters localState(JointVertexParameters const& p, unsigned int i) {
    JointVertexTrackParameters local;
    local << p[3+3*i],p[4+3*i],p[5+3*i],p[0],p[1];return local;
  }
  JointVertexTrackCovariance measurement() {
    JointVertexTrackCovariance factor=JointVertexTrackCovariance::Zero();
    factor.diagonal()<<.0002,.00003,.00005,.3,.4;
    factor(3,0)=.08;factor(4,1)=-.05;
    return factor*factor.transpose();
  }
  struct ParaxialField {
    double detectorZ=15000.;
    double curvature=.0003;
    bool derivatives=true;
    JointVertexPrediction operator()(unsigned int, JointVertexTrackParameters const& p, double z) const {
      double const length=detectorZ-z;
      JointVertexPrediction prediction;
      prediction.valid=true;
      prediction.parameters=p;
      prediction.parameters[1]+=curvature*p[0]*length;
      prediction.parameters[3]+=p[1]*length+.5*curvature*p[0]*length*length;
      prediction.parameters[4]+=p[2]*length;
      // Correlated process noise, separate from the measured detector error.
      JointVertexTrackCovariance factor=JointVertexTrackCovariance::Zero();
      factor.diagonal()<<.00002,.00002,.00003,.1,.2;
      factor(3,1)=.1;factor(4,2)=.2;
      prediction.noise=factor*factor.transpose();
      prediction.jacobianValid=prediction.zDerivativeValid=derivatives;
      prediction.jacobian=JointVertexTrackCovariance::Identity();
      prediction.jacobian(1,0)=curvature*length;
      prediction.jacobian(3,0)=.5*curvature*length*length;
      prediction.jacobian(3,1)=prediction.jacobian(4,2)=length;
      prediction.zDerivative<<0.,-curvature*p[0],0.,-p[1]-curvature*p[0]*length,-p[2];
      return prediction;
    }
  };
  std::array<JointVertexTrack,2> measure(JointVertexParameters const& p, ParaxialField const& transport) {
    std::array<JointVertexTrack,2> data;
    for (unsigned int i=0;i<2;++i) {
      data[i].parameters=transport(i,localState(p,i),p[2]).parameters;
      data[i].covariance=measurement();
    }
    return data;
  }
  void linearReference() {
    std::array<VertexTrack,2> reference;
    std::array<JointVertexTrack,2> detector;
    for (unsigned int i=0;i<2;++i) {
      reference[i].parameters<< (i ? -.02 : .01), (i ? -.02 : .02), .01, (i ? -2. : 3.),1.;
      JointVertexTrackCovariance factor=JointVertexTrackCovariance::Identity();
      factor.diagonal()<<.001,.002,.003,2.,3.;
      factor(3,0)=.5;factor(4,1)=-.8;
      reference[i].covariance=factor*factor.transpose();
      detector[i]={reference[i].parameters,reference[i].covariance};
    }
    auto const expected=refitCommonVertex(reference,2.);
    require(expected.valid,"CommonVertexRefit reference rejected");
    auto transport=[&](unsigned int i, JointVertexTrackParameters const& p, double z) {
      JointVertexPrediction prediction;
      prediction.valid=true;prediction.parameters=p;
      prediction.parameters.tail<2>()-=reference[i].parameters.segment<2>(1)*z;
      prediction.jacobianValid=prediction.zDerivativeValid=true;
      prediction.jacobian.setIdentity();
      prediction.zDerivative.tail<2>()=-reference[i].parameters.segment<2>(1);
      return prediction;
    };
    JointVertexParameters seed=JointVertexParameters::Zero();
    for (unsigned int i=0;i<2;++i) seed.segment<3>(3+3*i)=reference[i].parameters.head<3>();
    for (bool distant:{false,true}) {
      if(distant) {seed[0]=500.;seed[1]=-400.;seed[2]=5000.;seed[4]+=.2;seed[8]-=.3;}
      auto const fit=fitJointForwardVertex(seed,detector,2.,transport);
      require(fit.valid,"linear joint fit rejected: "+std::to_string(fit.status));
      require(fit.iterations<=3,"linear distant seed needs unnecessary iterations");
      require((fit.parameters.head<3>()-expected.displacement).norm()<1.e-7,"linear reference vertex mean mismatch");
      require((fit.covariance.topLeftCorner<3,3>()-expected.covariance).norm()<1.e-7,"linear reference vertex covariance mismatch");
      require(std::abs(fit.chi2-expected.chi2)<1.e-9,"linear reference chi2 mismatch");
      std::array<Eigen::Matrix<double,5,9>,2> designs;
      for(unsigned int i=0;i<2;++i) {
        auto& design=designs[i];design.setZero();
        design.block<3,3>(0,3+3*i).setIdentity();
        design(3,0)=design(4,1)=1.;design.block<2,1>(3,2)=-reference[i].parameters.segment<2>(1);
        require((design*fit.parameters-expected.tracks[i].parameters).norm()<1.e-8,"linear reference track mean mismatch");
        require((design*fit.covariance*design.transpose()-expected.tracks[i].covariance).norm()<1.e-7,"linear reference track covariance mismatch");
      }
      require((designs[0]*fit.covariance*designs[1].transpose()-expected.crossCovariance).norm()<1.e-7,"linear reference cross covariance mismatch");
      for (auto const& step : fit.steps) {
        require(step.direction.allFinite() && std::isfinite(step.expectedReduction),"missing finite direction telemetry");
        require(std::isfinite(step.jointUpdate) && step.jointUpdate>=0. &&
                std::abs(step.jointUpdate*step.jointUpdate-step.expectedReduction)<1.e-10*std::max(1.,step.expectedReduction),
                "joint motion telemetry disagrees with normal-equation decrement");
        auto frozenChi2=[&](JointVertexParameters const& parameters) {
          double chi2=.25*parameters.head<2>().squaredNorm();
          for(unsigned int i=0;i<2;++i) {
            JointVertexTrackParameters const residual=detector[i].parameters-designs[i]*parameters;
            chi2+=residual.dot(step.weights[i]*residual);
          }
          return chi2;
        };
        double analyticSlope=.5*step.parameters.head<2>().dot(step.direction.head<2>());
        for(unsigned int i=0;i<2;++i) {
          require((step.weights[i]-detector[i].covariance.inverse()).norm()<1.e-9*step.weights[i].norm(),
                  "recorded frozen weight differs from detector precision");
          require((step.designs[i]-designs[i]).norm()<1.e-12,"recorded design differs from actual linear mean derivative");
          require((step.predictedParameters[i]-designs[i]*step.parameters).norm()<1.e-9,
                  "recorded predicted mean differs from actual mean at the step");
          JointVertexTrackParameters const residual=detector[i].parameters-step.predictedParameters[i];
          JointVertexTrackParameters const detectorDirection=step.designs[i]*step.direction;
          analyticSlope-=2.*residual.dot(step.weights[i]*detectorDirection);
        }
        require(std::abs(analyticSlope+2.*step.expectedReduction)<1.e-8*std::max(1.,2.*step.expectedReduction),
                "recorded per-track means/designs/residuals disagree with analytic directional score");
        require(std::abs(frozenChi2(step.parameters)-step.chi2)<1.e-8*std::max(1.,step.chi2),
                "recorded state and frozen weights do not reproduce chi2");
        if(step.acceptedFraction>0.) {
          require(!step.trials.empty() && step.trials.back().valid &&
                  step.trials.back().fraction==step.acceptedFraction &&
                  step.trials.back().chi2==step.acceptedFrozenChi2,"accepted trial telemetry differs from accepted step");
          for(double h:{1./8192.,1./65536.}) {
            double const slope=(frozenChi2(step.parameters+h*step.direction)-
                                frozenChi2(step.parameters-h*step.direction))/(2.*h);
            require(std::abs(slope+2.*step.expectedReduction)<1.e-6*std::max(1.,2.*step.expectedReduction),
                    "recorded direction/reduction disagrees with frozen directional derivative");
          }
        } else require(step.trials.empty(),"converged step unexpectedly has line-search trials");
        for(auto const& trial:step.trials) {
          require(trial.valid && std::abs(frozenChi2(step.parameters+trial.fraction*step.direction)-trial.chi2)<
                  1.e-8*std::max(1.,trial.chi2),"recorded trial does not reproduce frozen chi2");
        }
      }
    }
    std::cout<<"Linear CommonVertexRefit means, full covariance, distant seeds and frozen-direction telemetry: PASS\n";
  }
  void correlatedConvergence() {
    auto transport=[](unsigned int i,JointVertexTrackParameters const& x,double z) {
      JointVertexPrediction p;p.valid=true;p.parameters=x;
      double const tx=i ? -.02 : .02;
      p.parameters[3]-=tx*z;p.parameters[4]-=.01*z;
      p.jacobianValid=p.zDerivativeValid=true;p.jacobian.setIdentity();
      p.zDerivative.tail<2>()<<-tx,-.01;return p;
    };
    JointVertexParameters truth;truth<<0.,0.,-1000.,.01,.02,.01,-.02,-.02,.01;
    std::array<JointVertexTrack,2> data;
    for(unsigned int i=0;i<2;++i) {
      data[i].parameters=transport(i,localState(truth,i),truth[2]).parameters;
      data[i].covariance=JointVertexTrackCovariance::Identity();
      data[i].covariance.diagonal().head<3>().setConstant(1.e-6);
    }
    double const rho=1.-1.e-6;
    data[0].covariance(0,1)=data[0].covariance(1,0)=rho*1.e-6;
    auto seed=truth;seed[3]+=.019*.001;seed[4]-=.019*.001;
    auto const fit=fitJointForwardVertex(seed,data,100.,transport);
    require(fit.valid && fit.iterations==2,"correlated fit did not perform its required correcting iteration");
    require(fit.steps.front().maxUpdate<.02 && fit.steps.front().jointUpdate>20.,
            "regression does not expose marginal-only false convergence");
    require(std::abs(fit.steps.front().jointUpdate-std::sqrt(2.*.019*.019/(1.-rho)))<1.e-5,
            "correlated joint motion differs from independent analytic result");
    require(fit.steps.front().acceptedFraction==1. && fit.chi2<1.e-8 && (fit.parameters-truth).norm()<1.e-8,
            "correlated correction did not reach the exact linear solution");
    require(fit.steps.back().maxUpdate<.02 && fit.steps.back().jointUpdate<.06,
            "reported convergence exceeds the combined marginal/joint tolerance");
    // Independent coordinates already inside the old tolerance must retain
    // their ordinary one-iteration convergence behavior.
    data[0].covariance(0,1)=data[0].covariance(1,0)=0.;
    auto const independent=fitJointForwardVertex(seed,data,100.,transport);
    require(independent.valid && independent.iterations==1 &&
            (independent.parameters-seed).norm()==0. && independent.steps.front().jointUpdate<.06,
            "joint guard changed independent-coordinate convergence inside tolerance");
    std::cout<<"Correlated marginal-only false convergence blocked; independent tolerance preserved: PASS\n";
  }
  void displacedCurved() {
    for(double curvature:{0.,.0003,-.0003}) for(double z:{-3500.,1700.}) {
      ParaxialField transport;transport.curvature=curvature;
      JointVertexParameters truth;truth<<.7,-.4,z,.01,.008,.012,-.015,-.006,-.009;
      auto const data=measure(truth,transport);
      auto seed=truth;seed[0]+=20.;seed[1]-=30.;seed[2]+=2500.;
      seed[3]*=1.15;seed[6]*=.8;seed[4]+=.005;seed[8]-=.007;
      JointVertexOptions precise;precise.convergence=1.e-6;
      auto fit=fitJointForwardVertex(seed,data,10000.,transport,precise);
      require(fit.valid,"curved displaced fit rejected: "+std::to_string(fit.status));
      require((fit.parameters-truth).norm()<.005,"curved displaced solution differs by "+std::to_string((fit.parameters-truth).norm())+" cm at curvature "+std::to_string(curvature)+" z "+std::to_string(z)+" fitted z "+std::to_string(fit.parameters[2]));
      require(fit.covariance.llt().info()==Eigen::Success,"joint covariance is not positive definite");
      for(auto const& step:fit.steps) require(step.acceptedFraction==0. || step.acceptedFrozenChi2<step.chi2,"line search did not descend");
      transport.derivatives=false;
      JointVertexOptions numericalOptions=precise;numericalOptions.maximumTransportCalls=4000;
      auto const numeric=fitJointForwardVertex(seed,data,10000.,transport,numericalOptions);
      require(numeric.valid,"numerical derivative fit rejected: "+std::to_string(numeric.status));
      require((numeric.parameters-fit.parameters).norm()<1.e-4,"analytic/numeric mean mismatch");
      require((numeric.covariance-fit.covariance).norm()<1.e-6*fit.covariance.norm(),"analytic/numeric covariance mismatch");
      require(!numeric.derivativeChecks.empty(),"finite difference checks missing");
      // Track order must not change the physical fit or vertex covariance.
      auto swappedData=data;std::swap(swappedData[0],swappedData[1]);
      auto swappedSeed=seed;swappedSeed.segment<3>(3)=seed.segment<3>(6);swappedSeed.segment<3>(6)=seed.segment<3>(3);
      transport.derivatives=true;
      auto const swapped=fitJointForwardVertex(swappedSeed,swappedData,10000.,transport,precise);
      require(swapped.valid && (swapped.parameters.head<3>()-fit.parameters.head<3>()).norm()<1.e-7,"track order changes vertex");
    }
    std::cout<<"Displaced straight/curved vertices, process noise, analytic/numeric derivatives, track swap: PASS\n";
  }
  void backtracking() {
    auto transport=[](unsigned int i, JointVertexTrackParameters const& p,double z) {
      JointVertexPrediction out;out.valid=true;out.parameters=p;
      double const tx=i ? -.02 : .01;
      out.parameters[0]=p[0]*p[0]*p[0];
      out.parameters[3]-=tx*z;out.parameters[4]-=.01*z;
      out.jacobianValid=out.zDerivativeValid=true;out.jacobian.setIdentity();
      out.jacobian(0,0)=3.*p[0]*p[0];out.zDerivative.tail<2>()<<-tx,-.01;return out;
    };
    JointVertexParameters truth;truth<<0.,0.,-1200.,.02,.01,.01,-.03,-.02,.01;
    std::array<JointVertexTrack,2> data;
    for(unsigned int i=0;i<2;++i) {data[i].parameters=transport(i,localState(truth,i),truth[2]).parameters;data[i].covariance=measurement();data[i].covariance(0,0)=1.e-14;data[i].covariance(0,3)=data[i].covariance(3,0)=0.;}
    auto seed=truth;seed[0]=100.;seed[1]=-80.;seed[2]=10000.;seed[3]=.001;seed[6]=-.001;
    auto const fit=fitJointForwardVertex(seed,data,100.,transport);
    require(fit.valid,"damped nonlinear fit failed: "+std::to_string(fit.status));
    require((fit.parameters-truth).norm()<1.e-5,"damped fit incorrect");
    require(std::any_of(fit.steps.begin(),fit.steps.end(),[](auto const& x){return x.acceptedFraction>0. && x.acceptedFraction<1.;}),"backtracking regression does not exercise damping");
    for(auto const& step:fit.steps) if(step.acceptedFraction>0.) {
      require(!step.trials.empty() && step.trials.back().valid &&
              step.trials.back().fraction==step.acceptedFraction &&
              step.trials.back().chi2==step.acceptedFrozenChi2,"nonlinear accepted trial telemetry mismatch");
      require(step.trials.size()>1 || step.acceptedFraction==1.,"rejected backtracking trials were lost");
    }
    std::cout<<"Distant nonlinear seed requires and passes backtracking: PASS\n";
  }
  void failures() {
    ParaxialField transport;
    JointVertexParameters seed;seed<<0.,0.,-1000.,.01,.01,.01,-.02,-.01,-.01;
    auto const data=measure(seed,transport);
    auto invalid=fitJointForwardVertex(seed,data,100.,[](unsigned int,auto const&,double){return JointVertexPrediction{};});
    require(!invalid.valid && invalid.status==JointForwardVertexFit::invalidTransport,"invalid transport accepted");
    auto opposite=seed;opposite[3]*=-1.;
    invalid=fitJointForwardVertex(opposite,data,100.,transport);
    require(!invalid.valid && invalid.status==JointForwardVertexFit::invalidCharge && invalid.transportCalls==0,"seed charge change accepted");
    auto zero=seed;zero[3]=0.;
    require(!fitJointForwardVertex(zero,data,100.,transport).valid,"zero q/p accepted");
    auto badData=data;badData[0].covariance(0,0)=-1.;
    require(!fitJointForwardVertex(seed,badData,100.,transport).valid,"indefinite detector covariance accepted");
    invalid=fitJointForwardVertex(seed,data,100.,[&](unsigned int i,auto const& x,double z){auto p=transport(i,x,z);p.noise=-.1*data[i].covariance;return p;});
    require(!invalid.valid && invalid.status==JointForwardVertexFit::invalidCovariance,"negative process noise rescued by measurement covariance");
    invalid=fitJointForwardVertex(seed,data,100.,[&](unsigned int i,auto const& x,double z){auto p=transport(i,x,z);p.jacobian(0,0)=std::numeric_limits<double>::quiet_NaN();return p;});
    require(!invalid.valid && invalid.status==JointForwardVertexFit::invalidDerivative,"nonfinite analytic derivative accepted");
    invalid=fitJointForwardVertex(seed,data,100.,[&](unsigned int i,auto const& x,double z){auto p=transport(i,x,z);p.zDerivativeValid=false;p.valid=z==seed[2];return p;});
    require(!invalid.valid && invalid.status==JointForwardVertexFit::invalidDerivative,"unresolved numerical derivative accepted");
    auto singular=[](unsigned int,auto const& x,double){JointVertexPrediction p;p.valid=true;p.parameters=x;p.jacobianValid=p.zDerivativeValid=true;p.jacobian.setIdentity();return p;};
    invalid=fitJointForwardVertex(seed,data,100.,singular);
    require(!invalid.valid && invalid.status==JointForwardVertexFit::singularNormal,"unidentifiable z accepted");
    JointVertexOptions bounded;bounded.maximumTransportCalls=1;
    invalid=fitJointForwardVertex(seed,data,100.,transport,bounded);
    require(!invalid.valid && invalid.status==JointForwardVertexFit::transportLimit && invalid.transportCalls==1,"transport call budget exceeded");
    bounded.maximumTransportCalls=10;
    transport.derivatives=false;
    invalid=fitJointForwardVertex(seed,data,100.,transport,bounded);
    require(!invalid.valid && invalid.status==JointForwardVertexFit::transportLimit && invalid.transportCalls==10,"numerical derivative budget exceeded");
    std::cout<<"Invalid transport, charge, covariance, derivative, singularity and exact call limits: PASS\n";
  }
  void changingVariance() {
    auto transport=[](unsigned int i,JointVertexTrackParameters const& x,double z) {
      JointVertexPrediction p;p.valid=true;p.parameters=x;
      double const tx=i ? -.02 : .01;
      p.parameters[3]-=tx*z;p.parameters[4]-=.01*z;
      p.jacobianValid=p.zDerivativeValid=true;p.jacobian.setIdentity();
      p.zDerivative.tail<2>()<<-tx,-.01;
      p.noise(0,0)=.2*x[0]*x[0];return p;
    };
    JointVertexParameters truth;truth<<0.,0.,-1500.,.02,.01,.01,-.03,-.02,.01;
    std::array<JointVertexTrack,2> data;
    for(unsigned int i=0;i<2;++i) {data[i].parameters=transport(i,localState(truth,i),truth[2]).parameters;data[i].covariance=JointVertexTrackCovariance::Identity();data[i].covariance(0,0)=1.e-6;}
    auto seed=truth;seed[3]*=2.;seed[6]*=2.;seed[2]+=500.;
    auto const fit=fitJointForwardVertex(seed,data,100.,transport);
    require(fit.valid && (fit.parameters-truth).norm()<1.e-8,"state-dependent variance moved the moment estimate");
    std::cout<<"State-dependent Q does not become an additional measurement: PASS\n";
  }
}
int main() {
  linearReference();correlatedConvergence();displacedCurved();backtracking();failures();changingVariance();
  std::cout<<"Joint forward vertex fit regression suite passed\n";
}
