# SHIFT muon momentum reconstruction: next steps

This note is the handoff for the next development session. Do not start by
retuning energy loss or applying a momentum calibration. The tests through
2026-08-11 localized the dominant pull to track topology and curvature lever
arm, not to the mean material model.

## Current evidence

- Truth SimHit-to-SimHit propagation closes locally to about -0.4% on the
  100-event probe. Backward propagation from the first hit to the target is
  low by about 3.4%, but that is too small to explain the poor one-leg tracks.
- Replacing the approximate SteppingHelix material model with material queried
  from the Geant4 geometry was stable but had essentially no paired effect on
  reconstructed momentum. Keep this implementation optional and off by
  default.
- Preferring a genuine two-leg traversing track increased the traversing
  population from 5 to 25 truth-matched muons in the same 100-event sample.
  Their median residuals were approximately -0.7% in pT and -7.2% in pz/total
  momentum. The remaining inclusive pull came from DSA and cosmic fallbacks.
- The retained representative order is traversing, then DSA, then cosmic.
  Source labels in the NanoAOD table are already `0=DSA`, `1=traversing`, and
  `2=cosmic`; do not remove fallback rows or obscure this provenance.

Validation outputs are under:

- `/eos/home-j/jniedzie/shift_cmssw/experiments/topology_first_100/`
- `/eos/home-j/jniedzie/shift_cmssw/experiments/topology_first_v2_100/`
- `/eos/home-j/jniedzie/shift_cmssw/experiments/first_principles_100/`

## Required data model

Expose reconstruction quality explicitly so an analyzer can choose its own
efficiency-versus-resolution working point. Keep at least these distinctions:

1. A high-quality, full-lever-arm traversing track.
2. A lower-quality DSA fallback.
3. A last-resort cosmic fallback.

The existing `source` column is the authoritative provenance. Add a stable,
self-documenting quality column such as `isFullLeverArm` or a small enumerated
`quality` value only if it makes analyzer selection clearer; it must be derived
from reconstructed topology and never from generator matching. Do not drop the
existing source, refit-validity, station-count, hit-count, and lever-arm
diagnostics.

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

## Target-constrained fit design

Implement a single coherent fit parameterized at the known target position
(nominally x=y=0 and z=+14800 cm for the current beam-B campaign; support the
opposite target side through configuration). Forward-propagate that state
through the configured magnetic field, geometry, and material to all retained
muon measurements and minimize their residuals together with the target
constraint.

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

## Validation order

1. Add columns/variants without changing the existing unconstrained values.
2. Run a 20-event smoke test for validity, schema, and crashes.
3. Run the grouped 100-event sample used above. Report matched counts and
   pT, pz, total-p, eta, and dimuon-mass residuals separately for traversing,
   DSA, and cosmic sources, and separately for constrained/unconstrained fits.
4. Check an LLP-like displaced control: the constrained fit should be visibly
   incompatible or disabled, while the unconstrained result remains unchanged.
5. Request a full production only if the constrained prompt result materially
   improves the 100-event paired comparison without sculpting efficiency or
   changing the LLP-safe branches.

The final analysis policy should be configurable: a precision prompt selection
may require full-lever-arm tracks and/or valid target-constrained fits, whereas
an inclusive or LLP analysis may retain the explicitly labeled fallback and
unconstrained tracks.
