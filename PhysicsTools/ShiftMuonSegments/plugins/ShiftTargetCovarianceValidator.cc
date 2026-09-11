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
#include <Eigen/Core>
#include <fstream>
#include <iomanip>
#include <cmath>

class ShiftTargetCovarianceValidator : public edm::one::EDAnalyzer<> {
public:
  explicit ShiftTargetCovarianceValidator(edm::ParameterSet const& p)
      : tracks_(consumes<reco::TrackCollection>(p.getParameter<edm::InputTag>("tracks"))),
        field_(esConsumes()), output_(p.getParameter<std::string>("output")),
        targetZ_(p.getParameter<double>("targetZ")), done_(false) {}
  void analyze(edm::Event const& event, edm::EventSetup const& setup) override {
    if (done_) return;
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
      for(int variant : {0,1,2}) {
        bool const consistent=variant>0;
        Geant4ePropagator prop(&field,"mu",oppositeToMomentum,1.,.2,50000.);
        prop.setUseConsistentBackwardCovariance(consistent);
        prop.setUseMeanEnergyLossJacobian(variant==2);
        auto propagate=[&](AlgebraicVector5 const& v, AlgebraicSymMatrix55 const& c) {
          TrajectoryStateOnSurface state(LocalTrajectoryParameters(v,pzSign),LocalTrajectoryError(c),*plane,&field);
          auto result=prop.propagate(state,*target);
          if(!result.isValid() || !result.hasError()) throw cms::Exception("CovarianceAudit")<<"Propagation failed";
          return result;
        };
        auto baseline=propagate(par,covariance), zero=propagate(par,tiny);
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
          out<<"{\"consistent\":"<<(consistent?"true":"false")<<",\"energyJacobian\":"<<(variant==2?"true":"false")<<",\"stepScale\":"<<stepScale<<",\"maximumNormalizedError\":"<<error;
          auto matrix=[&](char const* name, auto const& m){out<<",\""<<name<<"\":[";for(unsigned i=0;i<5;++i){if(i)out<<',';out<<'[';for(unsigned j=0;j<5;++j){if(j)out<<',';out<<m(i,j);}out<<']';}out<<']';};
          matrix("numerical",numerical);matrix("transported",transported);matrix("jacobian",jacobian);out<<'}';
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
};
DEFINE_FWK_MODULE(ShiftTargetCovarianceValidator);
