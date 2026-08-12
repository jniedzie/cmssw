# SHIFT muon momentum reconstruction: next steps

## Stopping point and restart checklist (2026-08-12)

The current canonical momentum is the unconstrained, precision-only,
full-lever-arm result from the **first** directional Kalman pass. There is no
empirical momentum calibration: the seed-momentum and energy-loss scales are
both `1.0`. Target transport has an explicit backward direction. Geometry
material on the target leg remains an optional diagnostic and is off by
default.

The last grouped 100-event test contains 24 generator-matched
double-traversing precision rows. Its total-momentum residual is `-4.057%`
with a central-68% half-width of `9.652%`. The readable output is:

`/tmp/shift_muon_iteration_ab/single_pass/step4/events_NanoAOD_part_0000.root`

The package builds successfully. Ordinary PTY-launched `cmsRun` still exits
with a startup ROOT segmentation signal in this environment, while the same
configuration completes normally under `gdb`; keep this execution issue
separate from the physics result.

Resume with the first Kalman pass itself. The most useful next test is to
publish same-surface q/p and covariance at four points: input seed, first-hit
prediction, forward-filter update, and source-facing smoothed state. Compare
those states to the true first precision SimHit, then split the paired change
by detector/station and integrated magnetic bending. This should tell us
whether the remaining bias enters through the seed prior, propagation between
measurements, or the smoother. Do not tune an energy-loss factor, q/p guard,
or truth-derived scale.

Before a larger production:

1. Run the same grouped 100-event closure and require the same 24-row
   population unless a reconstruction-validity change is understood.
2. Report pT, pz, total-p, and half-68% width, plus paired state-to-state
   changes; independent medians alone hid the pass-two tail migration.
3. Recheck constrained and unconstrained dimuon mass after any genuine
   first-pass improvement.
4. Run an LLP-like displaced control before changing the unconstrained
   analysis object.

## Status after v10 (2026-08-12)

The target-constrained muon hypothesis and the four exclusive muon quality
categories are now implemented. `ShiftDimuonVertex` also has an unqualified
unconstrained hypothesis, a `constrained*` hypothesis built only when both
muons are constrained, and ten unordered quality-pair flags. Mixed
constrained-unconstrained pairs are deliberately absent.

The v10 merged sample contains 5,770 events from commit `4df01fbe263`. The
missing events were not lost to a reconstruction crash: all jobs in Condor
cluster 17288419 failed before event processing because the shared
`pluginPhysicsToolsNanoAODPlugins.so` was being relinked and was observed as
`file too short`. Cluster 17288425 later recovered all 1,000 files and 10,000
events. The workflow now fingerprints the local CMSSW runtime at submission
and checks it before and after Step 4, so an overlapping `scram b` fails with
an explicit diagnostic instead of being mistaken for a physics failure.

Validity-aware v10 medians and central-68% half-widths for total momentum are:

| quality | unconstrained median / half-width | constrained median / half-width | matched constrained |
|---|---:|---:|---:|
| double traversing | -5.3% / 15.8% | -5.0% / 15.7% | 1,479 |
| DSA | -5.1% / 41.2% | -6.4% / 42.2% | 524 |
| cosmic | -12.9% / 43.4% | -9.2% / 39.2% | 124 |
| traversing | -7.1% / 24.7% | -6.5% / 14.2% | 4 (not interpretable) |

The unconstrained values in this table use exactly the rows for which the
constrained fit is valid. For double-traversing pairs, constrained dimuon mass
changes from -6.1% / 15.8% to -5.7% / 14.2% (140 matched pairs). Thus the
target constraint improves angular and some tail resolution, but does not
remove the roughly 5% core momentum-scale bias.

The most useful internal split is whether the precision-only refit was
selected. For double-traversing muons it gives -5.2% / 15.2% in total momentum
(1,410 rows), while the all-hit fallback gives -10.3% / 65.4% (98 rows). DSA
shows the same pattern: precision-only -4.7% / 34.0% versus all-hit fallback
-30.4% / 76.5%. Do not calibrate these populations together. The next focused
work should first diagnose the residual low-momentum bias in the
double-traversing precision population, then decide whether all-hit fallback
and low-station DSA belong in lower quality tiers rather than attempting one
global correction.

## Double-traversing precision closure (2026-08-12)

Same-surface diagnostics now separate production-to-first-hit loss from the
Kalman curvature estimate and the fitted loss across precision hits. On the
100-event probe, 24 truth-matched double-traversing precision-only rows give:

- canonical target momentum: -5.639% median, 9.483% central-68 half-width;
- true first-precision-hit / production momentum: -1.774% median;
- fitted upstream / true first-precision-hit momentum: -1.937% median;
- true loss across the precision hits: 39.084% median;
- fitted loss across the precision hits: 34.558% median.

These observations localize where the discrepancy appears, but they do not
justify an empirical momentum correction. A trial additive/fractional
calibration was explicitly rejected and removed from the producer: there are
no calibrated momentum branches or calibration constants in the data model.
The temporary filter-only/smoother-only energy-loss scale controls were also
removed; the filter and smoother again use the same unmodified propagator.

The procedural audit found two concrete limitations which should be fixed
before any further material-model scan:

- The configured `directionalRefitInitialMaxHitChi2` and
  `directionalRefitMaxHitChi2` values do not reject hits when calling
  `KFTrajectoryFitter` directly. The estimator only records each measurement
  estimate. Real outlier removal requires the `KFFittingSmoother` procedure
  (or an equivalent explicit remove-and-refit loop), with a minimum-hit and
  maximum-outlier policy.
- The second nonlinear pass is not an improvement for this sample. In v10 the
  first-pass median is -5.043%, versus -5.516% for valid second passes; on the
  100-event closure probe it is -4.057%, versus -5.639%. Iteration should be
  selected by an objective fit-quality/convergence criterion, not accepted
  merely because its q/p change is below 50%.

The detector-to-target material leg is valid for all 24 selected rows in the
100-event probe, but changes total momentum by only +0.001% (median), with a
median material-propagation path of 80.85 cm before the explicit CMS boundary.
Thus the remaining target extrapolation is not recovering the truth-level
production-to-first-hit loss; this is a transport-closure issue, not a reason
to add that loss back as a fitted constant.

The implemented state-contract check makes the convention explicit: the source-facing
smoothed state is checked against the first ordered fit-hit surface; its
filter-direction and target-direction cosines are published; and the
source-facing-to-target material leg is always propagated
opposite-to-momentum.  A state pointing toward the target is rejected rather
than silently transported with an ambiguous direction.  This deliberately
does not change an energy-loss scale or add a momentum correction.  Its effect
must be assessed on the 20-event smoke test and the grouped 100-event closure
before it is claimed as a scale improvement.

The local smoke test on Step-3 part 0000 completed with its ten available
events and a readable 1.4 MB EXONanoAOD output.  It contains four populated,
precision-refit ShiftMuon rows: all remain valid, their target-direction
cosines range from -0.9998 to -0.9968, and their filter-direction cosines from
0.9982 to 0.9999.  This validates the state-contract implementation and new
schema only; it is not a momentum-scale closure result.

### Robust-fit implementation test

An explicit `KFFittingSmoother`-style remove-and-refit loop is now implemented:
it tests every measurement above the configured estimate cut, removes the
candidate giving the lowest refitted chi2, and permits at most two removed
hits and at most 20% of the input. This fixes the previous no-op semantics of
the two max-hit-chi2 parameters without introducing a momentum correction.

The conventional CMSSW material-refit cut of 20 is not suitable as a
production default for these muon measurements. On the same 100-event probe
it removed hits in 8/24 selected double-traversing precision tracks and changed
the result from -5.639% / 9.483% to -7.270% / 11.477%. The loop is retained,
but the production cut remains 100000 (effectively disabled) until hit
compatibility is studied separately by detector and measurement dimension.

Requiring iteration two to improve reduced chi2 selected iteration one for
14/24 rows, but was neutral: -5.655% / 9.660%. That acceptance change was also
not made canonical. The next scale work should therefore target the transport
and measurement model, not tighter hit or iteration selection.

Changing target propagation to an explicit backward contract had no
measurable effect. Increasing the common Kalman energy-loss scale to 1.2 moved
the core toward zero but worsened the width and fallback tails; filter-only
and smoother-only tests did not close the bias. Do not make the 1.2 refit
scale canonical. The next step must instead repair the reconstruction
procedure: establish a consistent state direction and surface at every fit,
smoothing, and propagation boundary; ensure material is applied exactly once
with the correct propagation direction; and validate hit inclusion and fit
quality before changing any physics constants.

A 20-event split-leg stress test completed without a crash after rejecting
formally valid states below 0.1 GeV. Split legs remain optional diagnostics
because their curvature lever arm is poor; they are not a candidate momentum
estimator.

## Current evidence

- Truth SimHit-to-SimHit propagation closes locally to about -0.4% on the
  100-event probe. Backward propagation from the first hit to the target is
  low by about 3.4%, but that is too small to explain the poor one-leg tracks.
- Replacing the approximate SteppingHelix material model with material queried
  from the Geant4 geometry was stable but had essentially no paired effect on
  reconstructed momentum. Keep this implementation optional and off by
  default.
- A target-leg-only geometry-material ablation was then run on the grouped
  100-event closure sample. It left the hit fit unchanged and replaced only
  the approximate material correction from the source-facing precision state
  to the outer boundary. All 24 selected truth-matched double-traversing rows
  survived. The total-p residual changed from `-5.639%` (half-68% `9.483%`)
  to `-5.281%` (half-68% `9.356%`); the paired residual shift was only
  `+0.432%` (half-68% `0.153%`). Actual geometry raised target momentum by
  `0.476%`, versus `0.001%` for the approximate model, but this is much too
  small to close the scale. Keep `directionalRefitUseGeometryTargetMaterialEffects`
  optional and off by default; it repairs a real omitted material contribution
  but is not the scale solution.
- The two-pass Kalman logic was isolated next. Both passes use the same hits
  and outlier rejection is effectively disabled. CMSSW's DT/CSC muon RecHits
  retain the default `canImproveWithTrack() = false`, so pass two cannot
  relinearize them; it was feeding pass one's same-hit posterior back as a
  prior and processing the same measurements again. The established pass-two
  covariance inflation (`100`) gave `-5.639%`, a stronger prior (`10`) gave
  `-5.849%`, and a weak prior (`10000`) became unstable (`-53.969%` for valid
  pass-two rows). Thus pass two has no unique hit-only fixed point and depends
  on an arbitrary recycled covariance. In the established sample its paired
  median change was only `-0.032%`, but it was strongly anti-correlated with
  the pass-one residual (`rho=-0.897`) and redistributed a few tails enough to
  move the independent median. Canonical reconstruction now uses pass one
  (`directionalRefitUseSecondIteration=False`) while retaining pass two as an
  opt-in diagnostic. The resulting 24-row total-p closure is `-4.057%` with
  half-68% `9.652%`; this removes an invalid repeated-conditioning effect but
  does not yet close the remaining scale.
- Preferring a genuine two-leg traversing track increased the traversing
  population from 5 to 25 truth-matched muons in the same 100-event sample.
  Their median residuals were approximately -0.7% in pT and -7.2% in pz/total
  momentum. The remaining inclusive pull came from DSA and cosmic fallbacks.
- The retained representative order is traversing, then DSA, then cosmic.
  The `quality` labels are `0=cosmic`, `1=DSA`, `2=one-leg traversing`, and
  `3=full double-traversing`; `sourceIndex` retains the index in the source
  collection. Do not remove fallback rows or obscure this provenance.

Validation outputs are under:

- `/eos/home-j/jniedzie/shift_cmssw/experiments/topology_first_100/`
- `/eos/home-j/jniedzie/shift_cmssw/experiments/topology_first_v2_100/`
- `/eos/home-j/jniedzie/shift_cmssw/experiments/first_principles_100/`

## Data model requirements (implemented)

The quality categories and separate unconstrained/constrained hypotheses
described below are implemented. Preserve them during the remaining scale
work.

Expose reconstruction quality explicitly so an analyzer can choose its own
efficiency-versus-resolution working point. Keep these distinctions:

1. Cosmic fallback (`quality=0`).
2. DSA fallback (`quality=1`).
3. One-leg traversing track (`quality=2`).
4. Full double-traversing track (`quality=3`).

The `quality` column is the authoritative category and is derived only from
reconstructed topology, never generator matching. `sourceIndex` identifies the
row in the original track collection. Do not drop refit-validity,
station-count, hit-count, or lever-arm diagnostics.

For each accepted source hypothesis, preserve two distinct momentum variants:

- **Unconstrained:** detector hits, magnetic field, geometry, and material only.
  This is the LLP-safe result and must remain available for every source.
- **Target constrained:** include the known SHIFT production point/plane as a
  measurement in the fit. This is appropriate only for particles assumed to
  originate promptly at the target, for example prompt J/psi decays. It must
  never silently overwrite the unconstrained momentum.

Use explicit NanoAOD names (for example, the current unqualified momentum plus
`constrainedPt`, `constrainedEta`, `constrainedPhi`, `constrainedPz`,
`constrainedValid`, and `constrainedChi2`). If a separate table is cleaner,
keep a stable one-to-one index back to the unconstrained `ShiftMuon` row.

## Target-constrained fit status

The separate prompt-target hypothesis is implemented for both `ShiftMuon` and
`ShiftDimuonVertex`. It improves some angular and tail resolutions but leaves
the double-traversing core momentum scale near -5%, so it is not the remaining
scale solution. The design rules below remain requirements for future edits.

The current implementation applies the target measurement to the independently
smoothed detector state. Any future replacement with a fully simultaneous fit
should be parameterized at the known target position (nominally x=y=0 and
z=+14800 cm for the current beam-B campaign, with the opposite side supported
through configuration), then forward-propagate through the configured field,
geometry, and material to the retained muon measurements.

Do not implement the constraint by rescaling the final momentum, resetting
eta from a straight line, or applying a post-fit empirical correction. A
standard `SingleTrackVertexConstraint` is not automatically suitable here: its
implementation assumes tracker-scale vertices and can reject a state in the
near-zero field at the remote target. Verify the propagation and linearization
range before reusing it.

Make the target position and covariance explicit configuration parameters.
The constrained solution should report validity and constraint chi2, and it
should fall back to no constrained result rather than replacing the
unconstrained state after a failed or incompatible constraint.

## Validation order for future changes

1. Add diagnostics or an opt-in variant without changing the established
   unconstrained output.
2. Run a 20-event smoke test for validity, schema, and crashes.
3. Run the grouped 100-event sample used above. Report matched counts and
   pT, pz, total-p, eta, and dimuon-mass residuals separately for traversing,
   DSA, and cosmic sources, and separately for constrained/unconstrained fits.
4. Check an LLP-like displaced control: the constrained fit should be visibly
   incompatible or disabled, while the unconstrained result remains unchanged.
5. Request a full production only if a mechanism-based change materially
   improves the paired 100-event result without sculpting efficiency.

The final analysis policy should be configurable: a precision prompt selection
may require full-lever-arm tracks and/or valid target-constrained fits, whereas
an inclusive or LLP analysis may retain the explicitly labeled fallback and
unconstrained tracks.
