// Code style: AI-code. AI-generated private generative moment-iteration engine.
// Experimental mathematics reference for tests; not a production fitter.
// Frozen inner factors; this is NOT a parameter-dependent likelihood optimizer.
#ifndef SHIFT_PRIVATE_ITERATED_DETECTOR_FIT_ITERATION_H
#define SHIFT_PRIVATE_ITERATED_DETECTOR_FIT_ITERATION_H
#include "IteratedDetectorFit.h"
#include <functional>

namespace shift { namespace iterated {
using States = std::vector<State>;
using Linearize = std::function<bool(States const&,std::vector<Edge>&,std::vector<Hit>&,std::string&)>;
using Mean = std::function<bool(std::size_t,State const&,State&,std::string&)>;
struct IterationOptions {
  int maximumOuter = 8;
  int maximumHalvings = 8;
  double stateTolerance = 1.e-3;
  double momentTolerance = 1.e-3;
  Options linear;
};
struct IterationRecord {
  int iteration = 0;
  double objective = 0., stateStep = 0., momentChange = 0., fraction = 0.;
};
struct IterationResult {
  bool valid = false;
  int status = -10;
  std::string reason;
  // Physical generative states are here, not the slightly different linearized
  // proposal in fit.states. fit supplies the local fixed-model covariance.
  States states;
  std::vector<Eigen::VectorXd> kicks;
  Result fit;
  int iterations=0,linearizations=0,meanEvaluations=0,backtracks=0;
  double objective=std::numeric_limits<double>::infinity();
  double stateStep=std::numeric_limits<double>::infinity();
  double momentChange=std::numeric_limits<double>::infinity();
  std::vector<IterationRecord> history;
};

namespace iteration_detail {
inline bool project(States const& refs,std::vector<Edge> const& edges,
                    std::vector<Eigen::MatrixXd>& factors,std::vector<Eigen::VectorXd>& kicks) {
  factors.clear(); kicks.clear();
  if (edges.size()+1!=refs.size()) return false;
  double const eps=std::numeric_limits<double>::epsilon();
  for (std::size_t i=0;i<edges.size();++i) {
    Eigen::MatrixXd L;
    if (!edges[i].valid || !detail::noiseFactor(edges[i].Q,L)) return false;
    State const defect=refs[i+1]-edges[i].mean;
    Eigen::VectorXd kick=Eigen::VectorXd::Zero(L.cols());
    if (L.cols()) {
      Eigen::MatrixXd normalized=L;
      State scaled=defect;
      for(int row=0;row<5;++row) {
        double const s=std::sqrt(std::max(0.,edges[i].Q(row,row)));
        if(s>0.) { normalized.row(row)/=s; scaled[row]/=s; }
        else { normalized.row(row).setZero(); scaled[row]=0.; }
      }
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(normalized,Eigen::ComputeThinU|Eigen::ComputeThinV);
      if(svd.info()!=Eigen::Success || !normalized.allFinite() || !scaled.allFinite()) return false;
      double const limit=eps*std::max(normalized.rows(),normalized.cols())*svd.singularValues()[0];
      if((svd.singularValues().array()>limit).count()!=L.cols()) return false;
      kick=svd.matrixV()*svd.singularValues().cwiseInverse().asDiagonal()*svd.matrixU().transpose()*scaled;
    }
    State const represented=L*kick;
    for(int row=0;row<5;++row) {
      // Only roundoff in constructing/subtracting the same native mean is
      // tolerated in a deterministic direction; no covariance jitter/cut.
      double const roundoff=128.*eps*std::max({1.,std::abs(refs[i+1][row]),
                          std::abs(edges[i].mean[row]),std::abs(represented[row])});
      if(!std::isfinite(represented[row]) || std::abs(defect[row]-represented[row])>roundoff) return false;
    }
    factors.push_back(std::move(L)); kicks.push_back(std::move(kick));
  }
  return true;
}
inline bool objective(States const& states,std::vector<Eigen::VectorXd> const& kicks,
                      std::vector<Hit> const& hits,double& value) {
  value=0.;
  for(auto const& w:kicks) value+=0.5*w.squaredNorm();
  for(auto const& h:hits) {
    Eigen::MatrixXd W;
    if(!h.valid || h.site>=states.size() || !detail::whitener(h.R,W)) return false;
    Eigen::VectorXd const residual=W*(h.z-h.H*states[h.site]);
    value+=0.5*residual.squaredNorm();
  }
  return std::isfinite(value);
}
inline double relativeCovariance(Eigen::MatrixXd const& first,Eigen::MatrixXd const& second) {
  if(first.rows()!=second.rows() || first.cols()!=second.cols() || !first.allFinite() || !second.allFinite())
    return std::numeric_limits<double>::infinity();
  Eigen::VectorXd scale(first.rows());
  for(Eigen::Index i=0;i<scale.size();++i) {
    double const v=std::max(std::abs(first(i,i)),std::abs(second(i,i)));
    scale[i]=v>0. ? std::sqrt(v) : 1.;
  }
  Eigen::MatrixXd a=scale.cwiseInverse().asDiagonal()*first*scale.cwiseInverse().asDiagonal();
  Eigen::MatrixXd b=scale.cwiseInverse().asDiagonal()*second*scale.cwiseInverse().asDiagonal();
  return (a-b).norm()/std::max({1.,a.norm(),b.norm()});
}
inline double moments(States const& previousRefs,std::vector<Edge> const& previousEdges,
                      std::vector<Hit> const& previousHits,States const& refs,
                      std::vector<Edge> const& edges,std::vector<Hit> const& hits,Result const& fit) {
  if(edges.size()!=previousEdges.size() || hits.size()!=previousHits.size())
    return std::numeric_limits<double>::infinity();
  std::vector<Matrix5> scale(refs.size()),inverse(refs.size());
  for(std::size_t i=0;i<refs.size();++i) {
    for(int j=0;j<5;++j) if(!(fit.covariances[i](j,j)>0.)) return std::numeric_limits<double>::infinity();
    scale[i]=fit.covariances[i].diagonal().cwiseSqrt().asDiagonal();
    inverse[i]=fit.covariances[i].diagonal().cwiseSqrt().cwiseInverse().asDiagonal();
  }
  double change=0.;
  for(std::size_t i=0;i<edges.size();++i) {
    change=std::max(change,relativeCovariance(previousEdges[i].Q,edges[i].Q));
    Matrix5 const a=inverse[i+1]*previousEdges[i].F*scale[i];
    Matrix5 const b=inverse[i+1]*edges[i].F*scale[i];
    change=std::max(change,(a-b).norm()/std::max({1.,a.norm(),b.norm()}));
    State const affine=previousEdges[i].mean+previousEdges[i].F*(refs[i]-previousRefs[i]);
    change=std::max(change,(inverse[i+1]*(edges[i].mean-affine)).norm());
  }
  for(std::size_t i=0;i<hits.size();++i) {
    auto const& a=previousHits[i]; auto const& b=hits[i];
    if(a.site!=b.site || a.z.size()!=b.z.size()) return std::numeric_limits<double>::infinity();
    change=std::max(change,relativeCovariance(a.R,b.R));
    Eigen::MatrixXd W;
    if(!detail::whitener(b.R,W)) return std::numeric_limits<double>::infinity();
    Eigen::MatrixXd const ha=W*a.H*scale[b.site],hb=W*b.H*scale[b.site];
    change=std::max(change,(ha-hb).norm()/std::max({1.,ha.norm(),hb.norm()}));
    change=std::max(change,(W*((a.z-a.H*refs[a.site])-(b.z-b.H*refs[b.site]))).norm());
  }
  return change;
}
} // namespace iteration_detail

inline IterationResult iterate(States const& initial,Linearize const& linearize,Mean const& mean,
                               IterationOptions const& options=IterationOptions()) {
  IterationResult out;out.states=initial;
  auto fail=[&](int status,std::string const& reason) {out.valid=false;out.status=status;out.reason=reason;return out;};
  if(initial.empty() || options.maximumOuter<1 || options.maximumOuter>8 ||
      options.maximumHalvings<0 || options.maximumHalvings>8 || !(options.stateTolerance>0.) ||
      !(options.momentTolerance>0.)) return fail(-10,"invalid iteration inputs/limits");
  States previousRefs;
  std::vector<Edge> previousEdges;
  std::vector<Hit> previousHits;
  for(int iteration=0;iteration<options.maximumOuter;++iteration) {
    out.iterations=iteration+1;
    std::vector<Edge> edges;std::vector<Hit> hits;std::string error;
    ++out.linearizations;
    if(!linearize(out.states,edges,hits,error)) return fail(-10,"linearization: "+error);
    out.fit=solve(out.states,edges,hits,Prior{},options.linear);
    if(!out.fit.valid) return fail(-12,"linear solve: "+out.fit.reason);
    std::vector<Eigen::MatrixXd> L;
    if(!iteration_detail::project(out.states,edges,L,out.kicks)) return fail(-11,"reference outside process support");
    if(!iteration_detail::objective(out.states,out.kicks,hits,out.objective)) return fail(-10,"invalid frozen objective");
    out.stateStep=0.;
    for(std::size_t i=0;i<out.states.size();++i) {
      Eigen::MatrixXd W;
      if(!detail::whitener(out.fit.covariances[i],W)) return fail(-13,"state covariance is not SPD");
      out.stateStep=std::max(out.stateStep,(W*(out.fit.states[i]-out.states[i])).norm());
    }
    out.momentChange=iteration ? iteration_detail::moments(previousRefs,previousEdges,previousHits,
                    out.states,edges,hits,out.fit) : std::numeric_limits<double>::infinity();
    out.history.push_back({iteration,out.objective,out.stateStep,out.momentChange,0.});
    if(iteration && out.stateStep<=options.stateTolerance && out.momentChange<=options.momentTolerance) {
      out.valid=true;out.status=1;out.reason="generative moment fixed point";return out;
    }
    if(iteration+1==options.maximumOuter) return fail(-15,"outer iteration limit");
    std::vector<Eigen::VectorXd> targetKicks;
    Eigen::Index cursor=5;
    for(auto const& factor:L) {targetKicks.push_back(out.fit.parameters.segment(cursor,factor.cols()));cursor+=factor.cols();}
    double const prediction=out.objective-out.fit.objective;
    double const costRoundoff=64.*std::numeric_limits<double>::epsilon()*std::max(1.,out.objective);
    if(prediction < -costRoundoff) return fail(-10,"linear solve increases frozen objective");
    bool accepted=false;States candidate;
    for(int half=0;half<=options.maximumHalvings;++half) {
      double const fraction=std::ldexp(1.,-half);
      candidate.clear();candidate.push_back(out.states[0]+fraction*(out.fit.states[0]-out.states[0]));
      std::vector<Eigen::VectorXd> trialKicks;
      bool physical=candidate[0].allFinite();
      for(std::size_t i=0;physical && i<edges.size();++i) {
        State predicted;
        ++out.meanEvaluations;
        if(!mean(i,candidate[i],predicted,error) || !predicted.allFinite()) {physical=false;break;}
        trialKicks.push_back(out.kicks[i]+fraction*(targetKicks[i]-out.kicks[i]));
        candidate.push_back(predicted+L[i]*trialKicks.back());
        physical=candidate.back().allFinite();
      }
      double cost=0.;
      if(physical && iteration_detail::objective(candidate,trialKicks,hits,cost) &&
          cost<=out.objective-1.e-4*fraction*std::max(0.,prediction)+costRoundoff) {
        previousRefs=out.states;previousEdges=edges;previousHits=hits;
        out.states=std::move(candidate);out.kicks=std::move(trialKicks);out.objective=cost;
        out.history.back().fraction=fraction;accepted=true;break;
      }
      ++out.backtracks;
    }
    if(!accepted) return fail(-14,"no admissible generative descent");
  }
  return fail(-15,"outer iteration limit");
}
}} // namespace shift::iterated
#endif
