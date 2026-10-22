// Code style: AI-code. AI-generated private fixed-linear trajectory reference.
// Experimental mathematics reference for tests; not a production fitter.
// No transport, hit selection, charge choice, or outer/noise iteration lives here.
#ifndef SHIFT_PRIVATE_ITERATED_DETECTOR_FIT_H
#define SHIFT_PRIVATE_ITERATED_DETECTOR_FIT_H

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace shift {
namespace iterated {
using State = Eigen::Matrix<double, 5, 1>;
using Matrix5 = Eigen::Matrix<double, 5, 5>;

struct Edge {
  State mean = State::Zero();  // f_i(refs[i]) at the next measured surface
  Matrix5 F = Matrix5::Identity();
  Matrix5 Q = Matrix5::Zero(); // process covariance in NEXT surface's chart
  bool valid = true;
};
struct Hit {
  std::size_t site = 0;
  Eigen::VectorXd z;  // affine pseudo-observation z = y-h(ref)+H*ref
  Eigen::MatrixXd H;  // z approximately H*x[site], in that site's LOCAL chart
  Eigen::MatrixXd R;
  bool valid = true;
};
struct Prior {
  State mean = State::Zero();
  Matrix5 covariance = Matrix5::Identity();
  bool enabled = false; // only a genuinely independent prior; default NONE
};
struct Options {
  std::size_t maximumNodes = 256;
  Eigen::Index maximumParameters = 1024;
  Eigen::Index maximumRows = 4096;
};
enum Status {
  success = 1, invalidInput = -1, invalidNoise = -2,
  invalidMeasurement = -3, unobservable = -4, numericalFailure = -5,
  budgetExceeded = -6
};
struct Result {
  bool valid = false;
  int status = invalidInput;
  std::string reason;
  std::vector<State> states;
  std::vector<Matrix5> covariances;
  std::vector<Eigen::MatrixXd> stateJacobians;
  Eigen::MatrixXd latentCovariance;
  Eigen::VectorXd parameters; // delta-x0 followed by standard-normal kicks
  double objective = std::numeric_limits<double>::infinity();
  int rank = 0;
  Eigen::Index parameterCount = 0;
  Eigen::Index rowCount = 0;
  std::vector<int> processRanks;
  Matrix5 crossCovariance(std::size_t first, std::size_t second) const {
    return stateJacobians.at(first) * latentCovariance * stateJacobians.at(second).transpose();
  }
};

// Solve ONE immutable affine Gaussian model. No original fitted-track posterior
// is silently used as a prior. States satisfy exactly:
// x[i+1] = edge.mean + edge.F*(x[i]-refs[i]) + L[i]*w[i].
// The only penalties are each supplied hit, each w_i, and optional Prior ONCE.
inline Result solve(std::vector<State> const& refs,
                    std::vector<Edge> const& edges,
                    std::vector<Hit> const& hits,
                    Prior const& prior = Prior(), Options const& options = Options());

namespace detail {
// Dimensionless correlation charts make PSD/rank decisions independent of the
// different physical units of q/p, slopes, and positions. The tolerance only
// resolves roundoff-sized eigenvalues; it is not added to any covariance.
inline double tolerance(Eigen::MatrixXd const& matrix) {
  return 64. * std::numeric_limits<double>::epsilon() * matrix.rows() *
         std::max(1., matrix.cwiseAbs().maxCoeff());
}
inline bool symmetric(Eigen::MatrixXd const& matrix) {
  return matrix.allFinite() &&
         (matrix-matrix.transpose()).cwiseAbs().maxCoeff() <= tolerance(matrix);
}

inline bool noiseFactor(Matrix5 const& covariance, Eigen::MatrixXd& factor) {
  if (!covariance.allFinite()) return false;
  std::vector<int> active;
  State scales = State::Zero();
  for (int i=0; i<5; ++i) {
    if (covariance(i,i) < 0.) return false;
    if (covariance(i,i) > 0.) {
      active.push_back(i);
      scales[i] = std::sqrt(covariance(i,i));
    } else {
      // A PSD covariance with zero diagonal has an exactly deterministic row.
      if (covariance.row(i).cwiseAbs().maxCoeff() != 0. ||
          covariance.col(i).cwiseAbs().maxCoeff() != 0.) return false;
    }
  }
  if (active.empty()) { factor.resize(5,0); return true; }
  Eigen::Index const count = active.size();
  Eigen::MatrixXd correlation(count,count);
  for (Eigen::Index i=0; i<count; ++i)
    for (Eigen::Index j=0; j<count; ++j)
      correlation(i,j)=covariance(active[i],active[j])/scales[active[i]]/scales[active[j]];
  if (!symmetric(correlation)) return false;
  correlation=(0.5*(correlation+correlation.transpose())).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(correlation);
  if (eigen.info()!=Eigen::Success || !eigen.eigenvalues().allFinite()) return false;
  double const limit=tolerance(correlation);
  if (eigen.eigenvalues().minCoeff() < -limit) return false;
  int const rank=(eigen.eigenvalues().array()>limit).count();
  factor=Eigen::MatrixXd::Zero(5,rank);
  int column=0;
  for (Eigen::Index j=0; j<count; ++j) {
    if (eigen.eigenvalues()[j] <= limit) continue;
    for (Eigen::Index i=0; i<count; ++i)
      factor(active[i],column)=scales[active[i]]*eigen.eigenvectors()(i,j)*std::sqrt(eigen.eigenvalues()[j]);
    ++column;
  }
  return factor.allFinite();
}

inline bool whitener(Eigen::MatrixXd const& covariance, Eigen::MatrixXd& whitening) {
  if (!covariance.rows() || covariance.rows()!=covariance.cols() || !covariance.allFinite()) return false;
  Eigen::VectorXd scales(covariance.rows());
  for (Eigen::Index i=0; i<covariance.rows(); ++i) {
    if (!(covariance(i,i)>0.)) return false;
    scales[i]=std::sqrt(covariance(i,i));
  }
  Eigen::MatrixXd correlation(covariance.rows(),covariance.cols());
  for (Eigen::Index i=0; i<covariance.rows(); ++i)
    for (Eigen::Index j=0; j<covariance.cols(); ++j)
      correlation(i,j)=covariance(i,j)/scales[i]/scales[j];
  if (!symmetric(correlation)) return false;
  correlation=(0.5*(correlation+correlation.transpose())).eval();
  Eigen::LLT<Eigen::MatrixXd> llt(correlation);
  if (llt.info()!=Eigen::Success) return false;
  Eigen::MatrixXd inverseScale=scales.cwiseInverse().asDiagonal();
  whitening=llt.matrixL().solve(inverseScale);
  return whitening.allFinite();
}
} // namespace detail

inline Result solve(std::vector<State> const& refs, std::vector<Edge> const& edges,
                    std::vector<Hit> const& hits, Prior const& prior, Options const& options) {
  Result out;
  auto fail=[&](int status, char const* reason) {
    out.status=status; out.reason=reason; out.valid=false;
    // A failed model never publishes a partial path/covariance as a fit.
    out.states.clear(); out.covariances.clear(); out.stateJacobians.clear();
    out.latentCovariance.resize(0,0); out.parameters.resize(0);
    return out;
  };
  if (refs.empty() || edges.size()+1!=refs.size() || hits.empty())
    return fail(invalidInput,"node/edge/hit counts");
  if (refs.size()>options.maximumNodes)
    return fail(budgetExceeded,"node budget");
  for (auto const& ref:refs)
    if (!ref.allFinite()) return fail(invalidInput,"nonfinite reference");
  std::vector<Eigen::MatrixXd> factors;
  std::vector<Eigen::Index> starts;
  Eigen::Index parameters=5, rows=prior.enabled ? 5 : 0;
  for (auto const& edge:edges) {
    if (!edge.valid || !edge.mean.allFinite() || !edge.F.allFinite())
      return fail(invalidInput,"invalid edge");
    Eigen::MatrixXd factor;
    if (!detail::noiseFactor(edge.Q,factor)) return fail(invalidNoise,"non-PSD process covariance");
    starts.push_back(parameters); parameters+=factor.cols(); rows+=factor.cols();
    out.processRanks.push_back(factor.cols()); factors.push_back(std::move(factor));
  }
  std::vector<Eigen::MatrixXd> whitening;
  for (auto const& hit:hits) {
    if (!hit.valid || hit.site>=refs.size() || hit.z.size()==0 || hit.H.cols()!=5 ||
        hit.H.rows()!=hit.z.size() || hit.R.rows()!=hit.z.size() || hit.R.cols()!=hit.z.size() ||
        !hit.z.allFinite() || !hit.H.allFinite())
      return fail(invalidMeasurement,"hit dimensions/validity");
    Eigen::MatrixXd weight;
    if (!detail::whitener(hit.R,weight)) return fail(invalidMeasurement,"non-SPD measurement covariance");
    rows+=hit.z.size(); whitening.push_back(std::move(weight));
  }
  out.parameterCount=parameters; out.rowCount=rows;
  if (parameters>options.maximumParameters || rows>options.maximumRows)
    return fail(budgetExceeded,"dense solve budget");
  Eigen::MatrixXd priorWhitening;
  if (prior.enabled && (!prior.mean.allFinite() || !detail::whitener(prior.covariance,priorWhitening)))
    return fail(invalidInput,"non-SPD independent prior");

  std::vector<State> constants(refs.size());
  std::vector<Eigen::MatrixXd> derivatives(refs.size());
  constants[0]=refs[0];
  derivatives[0]=Eigen::MatrixXd::Zero(5,parameters);
  derivatives[0].leftCols(5).setIdentity();
  for (std::size_t i=0; i<edges.size(); ++i) {
    constants[i+1]=edges[i].mean+edges[i].F*(constants[i]-refs[i]);
    derivatives[i+1]=edges[i].F*derivatives[i];
    if (factors[i].cols()) derivatives[i+1].middleCols(starts[i],factors[i].cols())+=factors[i];
    if (!constants[i+1].allFinite() || !derivatives[i+1].allFinite())
      return fail(numericalFailure,"nonfinite affine trajectory");
  }
  Eigen::MatrixXd design=Eigen::MatrixXd::Zero(rows,parameters);
  Eigen::VectorXd rhs=Eigen::VectorXd::Zero(rows);
  Eigen::Index cursor=0;
  if (prior.enabled) {
    design.middleRows(cursor,5)=priorWhitening*derivatives[0];
    rhs.segment(cursor,5)=priorWhitening*(prior.mean-constants[0]); cursor+=5;
  }
  for (std::size_t i=0; i<factors.size(); ++i) {
    Eigen::Index const rank=factors[i].cols();
    if (rank) design.block(cursor,starts[i],rank,rank).setIdentity();
    cursor+=rank;
  }
  for (std::size_t i=0; i<hits.size(); ++i) {
    Eigen::Index const count=hits[i].z.size();
    design.middleRows(cursor,count)=whitening[i]*hits[i].H*derivatives[hits[i].site];
    rhs.segment(cursor,count)=whitening[i]*(hits[i].z-hits[i].H*constants[hits[i].site]);
    cursor+=count;
  }
  if (!design.allFinite() || !rhs.allFinite()) return fail(numericalFailure,"nonfinite weighted factors");
  Eigen::VectorXd scale(parameters);
  for (Eigen::Index i=0; i<parameters; ++i) {
    double const norm=design.col(i).stableNorm();
    if (!(norm>0.) || !std::isfinite(norm)) return fail(unobservable,"unobserved parameter");
    scale[i]=1./norm;
  }
  if (!scale.allFinite()) return fail(unobservable,"unresolved parameter scale");
  Eigen::MatrixXd scaled=design*scale.asDiagonal();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(scaled,Eigen::ComputeThinU|Eigen::ComputeThinV);
  if (svd.info()!=Eigen::Success || !svd.singularValues().allFinite())
    return fail(numericalFailure,"SVD failure");
  double const threshold=std::numeric_limits<double>::epsilon()*std::max(rows,parameters)*svd.singularValues()[0];
  out.rank=(svd.singularValues().array()>threshold).count();
  if (out.rank!=parameters) return fail(unobservable,"rank-deficient measured trajectory");
  Eigen::VectorXd singularInverse=svd.singularValues().cwiseInverse();
  out.parameters=scale.asDiagonal()*svd.matrixV()*singularInverse.asDiagonal()*svd.matrixU().transpose()*rhs;
  Eigen::MatrixXd covarianceRoot=scale.asDiagonal()*svd.matrixV()*singularInverse.asDiagonal();
  out.latentCovariance=covarianceRoot*covarianceRoot.transpose();
  Eigen::VectorXd residual=design*out.parameters-rhs;
  out.objective=0.5*residual.squaredNorm();
  if (!out.parameters.allFinite() || !out.latentCovariance.allFinite() || !std::isfinite(out.objective))
    return fail(numericalFailure,"nonfinite fit/covariance");
  out.stateJacobians=derivatives;
  for (std::size_t i=0; i<refs.size(); ++i) {
    State const state=constants[i]+derivatives[i]*out.parameters;
    Matrix5 const covariance=derivatives[i]*out.latentCovariance*derivatives[i].transpose();
    if (!state.allFinite() || !covariance.allFinite()) return fail(numericalFailure,"nonfinite fitted state");
    out.states.push_back(state);
    out.covariances.push_back((0.5*(covariance+covariance.transpose())).eval());
  }
  out.valid=true; out.status=success; out.reason="fixed linear factors solved";
  return out;
}

} // namespace iterated
} // namespace shift
#endif
