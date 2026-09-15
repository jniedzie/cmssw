// Diagnostic-only finite-difference audit. No truth products are consumed.
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryError.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCurvilinear.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToLocal.h"
#include <Eigen/Core>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

class ShiftTargetCovarianceValidator : public edm::one::EDAnalyzer<> {
public:
  explicit ShiftTargetCovarianceValidator(edm::ParameterSet const& p)
      : tracks_(consumes<reco::TrackCollection>(p.getParameter<edm::InputTag>("tracks"))),
        field_(esConsumes()), output_(p.getParameter<std::string>("output")),
        targetZ_(p.getParameter<double>("targetZ")), done_(false),
        sourceParameters_(p.existsAs<std::vector<double>>("sourceParameters") ? p.getParameter<std::vector<double>>("sourceParameters") : std::vector<double>{}),
        sourceZ_(p.existsAs<double>("sourceZ") ? p.getParameter<double>("sourceZ") : 0.),
        roundTripNoise_(p.existsAs<bool>("roundTripNoise") && p.getParameter<bool>("roundTripNoise")),
        useFieldGradient_(p.existsAs<bool>("useFieldGradient") && p.getParameter<bool>("useFieldGradient")) {}
  void analyze(edm::Event const& event, edm::EventSetup const& setup) override {
    if (done_) return;
    if (!sourceParameters_.empty()) {
      if (sourceParameters_.size()!=5) throw cms::Exception("Configuration")<<"Five source parameters required";
      auto const& field=setup.getData(field_);
      auto const source=Plane::build(GlobalPoint(0.,0.,sourceZ_),Surface::RotationType());
      auto const target=Plane::build(GlobalPoint(0.,0.,targetZ_),Surface::RotationType());
      std::ofstream out(output_);
      if (!out) throw cms::Exception("FileWriteError")<<output_;
      out<<std::setprecision(17)<<"[";bool first=true;
      if (roundTripNoise_) {
        // On one reversible mean path, Q_back = J_back Q_forward J_back^T.
        // Propagating the noisy forward endpoint backwards therefore gives
        // twice the covariance of propagating that same endpoint with C=0.
        // This isolates forward/backward covariance consistency without any
        // generated tracks, residual cuts, or finite-difference step choice.
        for (bool interface : {false,true}) {
          Geant4ePropagator forward(&field,"mu",alongMomentum,1.,.2,20000.);
          forward.setUseConsistentBackwardCovariance(true);
          forward.setUseMeanEnergyLossJacobian(true);
          forward.setUseUnquenchedIonizationVariance(true);
          forward.setUseFieldGradientJacobian(useFieldGradient_);
          forward.setUseMaterialInterfaceJacobian(interface);
          Geant4ePropagator backward(forward);
          backward.setPropagationDirection(oppositeToMomentum);
          AlgebraicVector5 parameters;AlgebraicSymMatrix55 zero;
          for (unsigned int i=0;i<5;++i) parameters[i]=sourceParameters_[i];
          TrajectoryStateOnSurface start(LocalTrajectoryParameters(parameters,targetZ_>sourceZ_?1.:-1.),
              LocalTrajectoryError(zero),*source,&field);
          auto const end=forward.propagate(start,*target);
          if (!end.isValid() || !end.hasError()) throw cms::Exception("CovarianceAudit")<<"Forward path failed";
          auto const back=backward.propagate(end,*source);
          TrajectoryStateOnSurface clean(end.localParameters(),LocalTrajectoryError(zero),*target,&field);
          auto const noise=backward.propagate(clean,*source);
          if (!back.isValid() || !noise.isValid() || !back.hasError() || !noise.hasError())
            throw cms::Exception("CovarianceAudit")<<"Backward path failed";
          if (!first) out<<',';first=false;
          out<<"{\"interface\":"<<interface;
          auto record=[&](char const* name,auto const& state) {
            out<<",\""<<name<<"\":{\"mean\":[";
            for (unsigned int i=0;i<5;++i) {if(i)out<<',';out<<state.localParameters().vector()[i];}
            out<<"],\"covariance\":[";
            for (unsigned int i=0;i<5;++i) {if(i)out<<',';out<<'[';
              for (unsigned int j=0;j<5;++j) {if(j)out<<',';out<<state.localError().matrix()(i,j);}out<<']';}
            out<<"]}";
          };
          record("start",start);record("forward",end);record("roundTrip",back);record("backwardNoise",noise);
          out<<'}';out.flush();
        }
        out<<"]\n";done_=true;return;
      }
      for (bool interface : {false,true}) {
        Geant4ePropagator prop(&field,"mu",alongMomentum,1.,.2,20000.);
        prop.setUseConsistentBackwardCovariance(true);prop.setUseMeanEnergyLossJacobian(true);
        prop.setUseUnquenchedIonizationVariance(true);prop.setUseMaterialInterfaceJacobian(interface);
        prop.setUseFieldGradientJacobian(useFieldGradient_);
        prop.setRecordTransportJacobian(true);
        for (int coordinate=-1;coordinate<5;++coordinate) for (double offset : (coordinate<0 ? std::vector<double>{0.} : std::vector<double>{-1.,-.1,.1,1.})) {
          AlgebraicVector5 parameters;AlgebraicSymMatrix55 zero;
          for (unsigned int i=0;i<5;++i) parameters[i]=sourceParameters_[i];
          if (coordinate>=0) parameters[coordinate]+=offset*(coordinate==0?std::abs(parameters[0])*1.e-3:(coordinate<3?1.e-5:.01));
          TrajectoryStateOnSurface state(LocalTrajectoryParameters(parameters,targetZ_>sourceZ_?1.:-1.),LocalTrajectoryError(zero),*source,&field);
          auto const result=prop.propagate(state,*target);
          if (!first) out<<',';first=false;
          out<<"{\"interface\":"<<interface<<",\"coordinate\":"<<coordinate<<",\"offset\":"<<offset<<",\"valid\":"<<result.isValid();
          if (result.isValid() && result.hasError()) {
            if (prop.transportJacobianValid()) {
              JacobianLocalToCurvilinear const input(state.surface(),state.localParameters(),field);
              JacobianCurvilinearToLocal const output(result.surface(),result.localParameters(),field);
              AlgebraicMatrix55 const jacobian=output.jacobian()*prop.curvilinearTransportJacobian()*input.jacobian();
              out<<",\"jacobian\":[";
              for (unsigned int i=0;i<5;++i) {if(i)out<<',';out<<'[';
                for(unsigned int j=0;j<5;++j) {if(j)out<<',';out<<jacobian(i,j);}out<<']';}out<<']';
            }
            out<<",\"maximumInterfaceRelativeCurvatureSigma\":"<<prop.maximumInterfaceRelativeCurvatureSigma();
            out<<",\"mean\":[";for (unsigned int i=0;i<5;++i) {if(i)out<<',';out<<result.localParameters().vector()[i];}out<<"],\"covariance\":[";
            for (unsigned int i=0;i<5;++i) {if(i)out<<',';out<<'[';for(unsigned int j=0;j<5;++j) {if(j)out<<',';out<<result.localError().matrix()(i,j);}out<<']';}out<<']';
          }
          out<<'}';out.flush();
        }
      }
      out<<"]\n";done_=true;return;
    }
    auto const tracks=event.getHandle(tracks_);
    if (!tracks.isValid()) return;
    auto const& field=setup.getData(field_);
    for (auto const& track : *tracks) {
      bool const outer=targetZ_*(track.outerPosition().z()-track.innerPosition().z())>0.;
      auto const pos=outer?track.outerPosition():track.innerPosition();
      auto const mom=outer?track.outerMomentum():track.innerMomentum();
      if (!(std::sqrt(mom.mag2())>20.) || std::abs(mom.z())<1.) continue;
      double const sign=targetZ_*mom.z()>0.?-1.:1.;
      auto const plane=Plane::build(GlobalPoint(pos.x(),pos.y(),pos.z()),Surface::RotationType());
      auto const target=Plane::build(GlobalPoint(0.,0.,targetZ_),Surface::RotationType());
      AlgebraicVector5 par;
      par[0]=sign*track.charge()/std::sqrt(mom.mag2()); par[1]=mom.x()/mom.z(); par[2]=mom.y()/mom.z(); par[3]=0.; par[4]=0.;
      double const pzSign=sign*mom.z()>0.?1.:-1.;
      // Controlled, small uncertainties; subtract common process noise to
      // isolate transport of the input covariance. Mean state is identical.
      AlgebraicSymMatrix55 covariance, tiny;
      double const sigma[5]={std::abs(par[0])*.01, 1.e-4, 1.e-4, .01, .01};
      for(unsigned i=0;i<5;++i) { covariance(i,i)=sigma[i]*sigma[i]; tiny(i,i)=covariance(i,i)*1.e-6; }
      std::ofstream out(output_);
      if (!out) throw cms::Exception("FileWriteError")<<output_;
      out<<std::setprecision(17)<<"{\"event\":"<<event.id().event()<<",\"pzSign\":"<<pzSign<<",\"momentum\":"<<std::sqrt(mom.mag2())<<",\"variants\":[";
      bool first=true;
      for(int variant : {0,1,2,3}) {
        bool const consistent=variant>0;
        Geant4ePropagator prop(&field,"mu",oppositeToMomentum,1.,.2,50000.);
        prop.setUseConsistentBackwardCovariance(consistent);
        prop.setUseMeanEnergyLossJacobian(variant>=2);
        prop.setUseFieldGradientJacobian(variant==3);
        prop.setRecordTransportJacobian(true);
        auto propagate=[&](AlgebraicVector5 const& v, AlgebraicSymMatrix55 const& c) {
          TrajectoryStateOnSurface state(LocalTrajectoryParameters(v,pzSign),LocalTrajectoryError(c),*plane,&field);
          auto result=prop.propagate(state,*target);
          if(!result.isValid() || !result.hasError()) throw cms::Exception("CovarianceAudit")<<"Propagation failed";
          return result;
        };
        auto baseline=propagate(par,covariance);
        Eigen::Matrix<double,5,5> accumulated=Eigen::Matrix<double,5,5>::Zero();
        bool const accumulatedValid=prop.transportJacobianValid();
        if (accumulatedValid) {
          JacobianLocalToCurvilinear const input(*plane,LocalTrajectoryParameters(par,pzSign),field);
          JacobianCurvilinearToLocal const output(baseline.surface(),baseline.localParameters(),field);
          AlgebraicMatrix55 const matrix=output.jacobian()*prop.curvilinearTransportJacobian()*input.jacobian();
          for (unsigned int i=0;i<5;++i) for (unsigned int j=0;j<5;++j) accumulated(i,j)=matrix(i,j);
        }
        auto zero=propagate(par,tiny);
        Eigen::Matrix<double,5,5> transported, input;
        input.setZero();
        for(unsigned i=0;i<5;++i) {
          input(i,i)=covariance(i,i)-tiny(i,i);
          for(unsigned j=0;j<5;++j) transported(i,j)=baseline.localError().matrix()(i,j)-zero.localError().matrix()(i,j);
        }
        for(double stepScale : {1.,.5}) {
          Eigen::Matrix<double,5,5> jacobian;
          for(unsigned j=0;j<5;++j) {
            double const h=stepScale*(j==0?std::abs(par[0])*1.e-3:(j<3?1.e-5:.01));
            auto plus=par,minus=par; plus[j]+=h;minus[j]-=h;
            auto a=propagate(plus,tiny),b=propagate(minus,tiny);
            for(unsigned i=0;i<5;++i) jacobian(i,j)=(a.localParameters().vector()[i]-b.localParameters().vector()[i])/(2.*h);
          }
          Eigen::Matrix<double,5,5> numerical=jacobian*input*jacobian.transpose();
          double error=0.;
          for(unsigned i=0;i<5;++i) for(unsigned j=0;j<5;++j)
            error=std::max(error,std::abs(transported(i,j)-numerical(i,j))/std::sqrt(numerical(i,i)*numerical(j,j)));
          if(!first) out<<',';
          first=false;
          out<<"{\"consistent\":"<<(consistent?"true":"false")<<",\"energyJacobian\":"<<(variant>=2?"true":"false")<<",\"fieldGradient\":"<<(variant==3?"true":"false")<<",\"stepScale\":"<<stepScale<<",\"maximumNormalizedError\":"<<error;
          auto matrix=[&](char const* name, auto const& m){out<<",\""<<name<<"\":[";for(unsigned i=0;i<5;++i){if(i)out<<',';out<<'[';for(unsigned j=0;j<5;++j){if(j)out<<',';out<<m(i,j);}out<<']';}out<<']';};
          matrix("numerical",numerical);matrix("transported",transported);matrix("jacobian",jacobian);
          if (accumulatedValid) matrix("accumulatedJacobian",accumulated);
          out<<'}';
        }
      }
      out<<"]}\n"; done_=true; return;
    }
  }
  void endJob() override { if(!done_) throw cms::Exception("CovarianceAudit")<<"No eligible detector track"; }
private:
  edm::EDGetTokenT<reco::TrackCollection> tracks_;
  edm::ESGetToken<MagneticField,IdealMagneticFieldRecord> field_;
  std::string output_; double targetZ_; bool done_;
  std::vector<double> sourceParameters_;
  double sourceZ_;
  bool roundTripNoise_;
  bool useFieldGradient_;
};
DEFINE_FWK_MODULE(ShiftTargetCovarianceValidator);
