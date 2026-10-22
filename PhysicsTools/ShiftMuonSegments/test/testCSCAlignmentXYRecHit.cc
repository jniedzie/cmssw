// Code style: human-code.
// Event-free private covariance/dispatch oracle; no fitting or propagation.
#include "CSCAlignmentXYRecHit.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"
#include "PhysicsTools/ShiftMuonSegments/interface/ExactMuonFitInput.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace {
  unsigned int checks = 0;
  void require(bool value, char const* message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
  }
  void close(double value, double expected, char const* message) {
    double const bound = 64. * std::numeric_limits<double>::epsilon() *
                         std::max(1., std::abs(expected));
    require(std::isfinite(value) && std::abs(value - expected) <= bound, message);
  }
  class TestGeometry final : public GeomDet {
  public:
    TestGeometry() : GeomDet(Plane::build(Surface::PositionType(0., 0., 0.), Surface::RotationType())) {
      setDetId(CSCDetId(1, 1, 1, 1, 1));
    }
    void alignment(AlgebraicSymMatrix66 const& covariance) {
      AlignmentPositionError const alignment{GlobalErrorExtended(covariance)};
      // Native GeomDet allocates only for initially valid APE. Install a seed
      // before zeroing so zero-present and absent-pointer cases are distinct.
      if (!alignment.valid())
        setAlignmentPositionError(AlignmentPositionError(.1f, .1f, .1f));
      require(setAlignmentPositionError(alignment) == alignment.valid(), "test APE installation");
    }
  };
  class OtherPoint final : public RecHit2DLocalPos {
  public:
    OtherPoint() : RecHit2DLocalPos(CSCDetId(1, 1, 1, 1, 1)) {}
    OtherPoint* clone() const override { return new OtherPoint(*this); }
    LocalPoint localPosition() const override { return LocalPoint(1.f, 2.f, 0.f); }
    LocalError localPositionError() const override { return LocalError(.04f, .007f, .09f); }
  };
  CSCRecHit2D makeHit(float x) {
    CSCRecHit2D::ADCContainer adc;
    std::vector<float> samples{10.f, 11.f, 12.f, 13.f};
    adc.put(10, samples.begin(), samples.end());
    return CSCRecHit2D(CSCDetId(1, 1, 1, 1, 1), LocalPoint(x, 8.f, 0.f),
                      LocalError(.04f, .007f, .09f), {10}, adc, {11},
                      5.f, .2f, .01f, 7);
  }
  void matrixEqual(AlgebraicSymMatrix const& a, AlgebraicSymMatrix const& b, char const* message) {
    require(a.num_row() == b.num_row(), message);
    for (int i = 0; i < a.num_row(); ++i)
      for (int j = 0; j <= i; ++j) require(a[i][j] == b[i][j], message);
  }
  void identity(TransientTrackingRecHit const& a, TransientTrackingRecHit const& b) {
    require(a.det() == b.det() && a.geographicalId() == b.geographicalId(), "geometry identity");
    require(a.isValid() == b.isValid() && a.dimension() == b.dimension(), "validity/dimension");
    auto const pa = a.parameters(), pb = b.parameters();
    for (int i = 0; i < pa.num_row(); ++i) require(pa[i] == pb[i], "measurement unchanged");
    auto const ha = a.projectionMatrix(), hb = b.projectionMatrix();
    for (int i = 0; i < ha.num_row(); ++i)
      for (int j = 0; j < ha.num_col(); ++j) require(ha[i][j] == hb[i][j], "H unchanged");
    require(typeid(*a.hit()) == typeid(CSCRecHit2D) && typeid(*b.hit()) == typeid(CSCRecHit2D),
            "persistent type");
    require(shift::isExactCSCMeasurementCopy(*a.hit(), *b.hit()), "complete raw CSC identity");
  }
  void dispatch(TransientTrackingRecHit const& hit, AlgebraicSymMatrix const& covariance) {
    AlgebraicVector2 z, projected;
    AlgebraicSymMatrix22 r, projectedError;
    ProjectMatrix<double, 5, 2> projection;
    AlgebraicVector5 state;
    AlgebraicSymMatrix55 stateError;
    for (int i = 0; i < 5; ++i) {
      state[i] = 0.;
      for (int j = 0; j <= i; ++j) stateError(i,j) = 0.;
    }
    KfComponentsHolder holder;
    holder.setup<2>(&z, &r, &projection, &projected, &projectedError, state, stateError);
    hit.getKfComponents(holder);
    for (int i = 0; i < 2; ++i) {
      require(z[i] == hit.parameters()[i], "KF measurement dispatch");
      for (int j = 0; j < 2; ++j) require(r(i,j) == covariance[i][j], "KF/batch R dispatch");
      for (int j = 0; j < 5; ++j)
        require(projection.matrix()(i,j) == hit.projectionMatrix()[i][j], "KF H dispatch");
    }
  }
}

int main() {
  try {
    // Full SPD pose covariance, including position/rotation cross terms.
    double lower[6][6] = {{.10,0,0,0,0,0}, {.02,.12,0,0,0,0},
                         {.03,-.01,.14,0,0,0}, {.001,.002,0,.003,0,0},
                         {-.002,.001,.003,.001,.004,0}, {.002,-.003,.001,.002,-.001,.005}};
    AlgebraicSymMatrix66 full;
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j <= i; ++j) {
        full(i,j) = 0.;
        for (int k = 0; k < 6; ++k) full(i,j) += lower[i][k] * lower[j][k];
      }
    AlgebraicSymMatrix66 translation, zero;
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j <= i; ++j) translation(i,j) = zero(i,j) = 0.;
    translation(0,0)=.03; translation(1,1)=.08; translation(0,1)=.007;
    for (auto const& c : {full, translation, zero}) {
      TestGeometry geometry;
      geometry.alignment(c);
      for (float x : {-17.f, 0.f, 17.f}) {
        auto const raw = makeHit(x);
        auto const original = MuonTransientTrackingRecHit::build(&geometry, &raw);
        require(shift::alignment_xy::correctedHit(original, false).get() == original.get(),
                "disabled preserves original pointer");
        auto const corrected = shift::alignment_xy::correctedHit(original, true);
        require(corrected.get() != original.get() && !corrected->canImproveWithTrack(), "enabled wrapper");
        identity(*original, *corrected);
        auto const r = corrected->parametersError(), rawR = raw.parametersError();
        // Independent rigid-motion relation at frozen pose (x,0), slopes(0,0):
        // delta u=delta X, delta v=delta Y+x*delta phi_z.
        double const xx = c(0,0);
        double const xy = c(0,1) + double(x) * c(0,5);
        double const yy = c(1,1) + 2. * double(x) * c(1,5) + double(x)*x*c(5,5);
        close(r[0][0], rawR[0][0]+xx, "independent xx oracle");
        close(r[0][1], rawR[0][1]+xy, "independent xy cross-term oracle");
        close(r[1][1], rawR[1][1]+yy, "independent yy oracle");
        require(r[0][0]>0. && r[1][1]>0. && r[0][0]*r[1][1]>r[0][1]*r[0][1], "effective R SPD");
        dispatch(*corrected, r);
        std::unique_ptr<TrackingRecHit> cloned(static_cast<TrackingRecHit const&>(*corrected).clone());
        require(dynamic_cast<shift::alignment_xy::CSCAlignmentXYRecHit20261022 const*>(cloned.get()),
                "virtual clone keeps wrapper");
        auto const* transient = dynamic_cast<TransientTrackingRecHit const*>(cloned.get());
        require(transient, "clone transient type");
        identity(*corrected, *transient);
        matrixEqual(transient->parametersError(), r, "clone covariance");
        dispatch(*transient, r);
        if (c(0,0)==0. && c(1,1)==0.) matrixEqual(r, rawR, "zero APE unchanged");
      }
    }
    TestGeometry noAlignment;
    auto const raw = makeHit(11.f);
    auto const original = MuonTransientTrackingRecHit::build(&noAlignment, &raw);
    auto const corrected = shift::alignment_xy::correctedHit(original, true);
    matrixEqual(corrected->parametersError(), original->parametersError(), "absent APE unchanged");
    auto const generic = GenericTransientTrackingRecHit::build(&noAlignment, &raw);
    require(shift::alignment_xy::correctedHit(generic,true).get()==generic.get(), "other transient class unchanged");
    OtherPoint const otherRaw;
    auto const other = MuonTransientTrackingRecHit::build(&noAlignment, &otherRaw);
    require(shift::alignment_xy::correctedHit(other,true).get()==other.get(), "other persistent class unchanged");
    auto const invalid = MuonTransientTrackingRecHit::specificBuild(&noAlignment, &raw);
    invalid->invalidateHit();
    require(shift::alignment_xy::correctedHit(invalid,true).get()==invalid.get(), "invalid hit unchanged");
    require(!shift::alignment_xy::correctedHit({},true), "null unchanged");
    std::cout << "{\"passed\":true,\"checks\":" << checks
              << ",\"fixture_pose\":\"native x,0; slopes0,0\",\"native_events\":0,\"propagations\":0}\n";
    return 0;
  } catch (std::exception const& error) {
    std::cerr << "oracle failed: " << error.what() << '\n';
    return 1;
  }
}
