# Topology-based ShiftMuon reconstruction plan

## Purpose and stopping point (2026-08-13)

Implementation status (2026-08-14): Phase 1 is implemented with the old
`quality` and quality-pair branches retained temporarily for compatibility.
Reconstruction inputs, duplicate cleaning, representative precedence, and
acceptance settings are unchanged. Compact duplicate-source availability
diagnostics from Phase 2 are included; the algorithm-removal study is not.

The goal is to replace the current mixture of reconstruction provenance and
geometry in `ShiftMuon_quality` with a small, stable set of measured detector
topologies. In parallel, we want to determine whether the three current input
track collections (DSA, strict traversing, and ordinary cosmic) can safely be
reduced to two reconstruction strategies.

This document is a plan only. No producer, cleaning, precedence, or acceptance
setting is changed as part of writing it.

The central distinction is:

- a **topology** describes where the selected track has recorded measurements;
- an **algorithm/provenance** records which CMSSW reconstruction produced it.

These should remain separate until controlled tests show that one algorithm is
reliably optimal for each topology.

## Current evidence

The first complete topology summary contains 3,912 selected ShiftMuons:

| Current quality | Total | near only | near + barrel | both endcaps | far only | other |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cosmic | 276 | 156 | 18 | 93 | 5 | 4 |
| DSA | 1,024 | 834 | 0 | 0 | 190 | 0 |
| traversing | 7 | 1 | 0 | 6 | 0 | 0 |
| double traversing | 2,605 | 0 | 2 | 2,603 | 0 | 0 |

The geometry correlations are strong:

- 2,603/2,605 current double-traversing rows have measurements in both endcaps;
- all 1,024 DSA rows have measurements in only one endcap;
- six of the seven current single-traversing rows already have both-endcap
  measurements and a median measured hit span of about 18 m;
- ordinary cosmic reconstruction is not a unique topology: it supplies
  near-only, near-plus-barrel, and 93 both-endcap rows.

The `abs(outerZ-innerZ) > 500 cm` condition only divides source-1 tracks into
quality 2 and quality 3. Removing it would merge the two labels, but it would
not alter reconstruction, remove a producer, or recover a failed traversing
candidate. Measured hit coverage (`entry/exit`, region counts, and hit span) is
a better geometry definition than the track inner/outer reference states.

Eta and momentum correlate with reach but do not define clean categories. The
four current quality populations overlap strongly in both variables, so eta or
momentum should be diagnostics and performance coordinates, not category
definitions.

## Proposed measured topology

Add one authoritative, mutually exclusive `topology` value derived only from
valid recorded muon-system measurements and `inferredSourceSide`:

| Value | Proposed name | Definition |
| ---: | --- | --- |
| 0 | `nearEndcapOnly` | near-endcap hits, no barrel hits, no far-endcap hits |
| 1 | `nearEndcapAndBarrel` | near-endcap and barrel hits, no far-endcap hits |
| 2 | `bothEndcaps` | hits in both endcaps, regardless of barrel hits |
| 3 | `farEndcapOnly` | far-endcap hits, no near-endcap or barrel hits |
| 4 | `unclassified` | barrel-only, barrel-plus-far, no classified hits, or an invalid orientation |

The requested four physics tags remain the main analysis categories. `unclassified`
is an explicit safety category rather than silently dropping the nine rare
cosmic rows currently outside the four definitions. Its population should be
reported in every validation. It can be reconsidered only after understanding
those events.

Keep the existing underlying diagnostics:

- raw `+z`, barrel, and `-z` hit counts;
- oriented near/barrel/far hit and station counts;
- entry/exit regions, positions, subdetectors, and measured hit span;
- direction/timing status;
- eta, momentum, refit validity, fit quality, and target-line diagnostics.

Do not interpret `exitRegion` as a physical stopping detector. Missing or
inefficient later measurements can make a through-going particle appear to
end early.

## Proposed reconstruction vocabulary

If the two-algorithm hypothesis is validated, use names describing fit intent
rather than the historical source:

- `throughGoingFit`: the current strict traversing reconstruction, intended
  for `bothEndcaps` tracks and fitting both detector legs as one trajectory;
- `singleLegFit`: the current DSA reconstruction, intended for `nearEndcapOnly`,
  `nearEndcapAndBarrel`, and provisionally `farEndcapOnly` tracks;
- keep `cosmicReco` as the ordinary-cosmic provenance name until its unique
  recovery contribution has been measured. Calling it `recoveryFit` now would
  assume the conclusion of that study.

Do not call a one-track both-endcap fit "double muon." If a user-facing object
name is needed, `throughGoing` or `bothEndcaps` is clearer than
`doubleTraversing`.

Keep provenance independently of topology during the study, for example:

- `recoAlgorithm`: DSA, traversing, or cosmic;
- `sourceIndex`: original source-collection index;
- `duplicateSourceMask`: algorithms represented in the duplicate component.

The final public data model may reduce these after validation, but the study
cannot be performed if provenance is removed first.

## Work plan

### Phase 1: simplify the labels without changing reconstruction

1. Derive and publish the five-value topology above.
2. Remove the 500 cm `outerZ-innerZ` split from the public category logic.
   Preserve `hitSpan` and add/use `abs(exitZ-entryZ)` as continuous diagnostics.
3. Keep all three source collections and the current representative selection
   unchanged for this phase.
4. Replace quality-pair flags in `ShiftDimuonVertex` with topology-pair
   information, or derive those pairs in analysis. Retain the old fields during
   a short compatibility period if existing analysis depends on them.
5. Validate that the topology counts are exhaustive, mutually exclusive, and
   unchanged by the label-only rewrite.

This phase answers whether the data model can be made simpler independently of
the more consequential algorithm-removal question.

### Phase 2: expose which algorithms were available before cleaning

For every duplicate component, record enough information to distinguish
"algorithm absent" from "algorithm existed but lost":

- source mask for DSA/traversing/cosmic candidates;
- candidate count from each source;
- best candidate from each source: topology, valid hits, stations, hit span,
  normalized chi2, directional-refit validity, timing status, and target DCA;
- reason a candidate was rejected before duplicate cleaning, using compact
  counters/status codes rather than verbose per-event logs.

This is required to understand:

- the 93 both-endcap rows currently retained as cosmic;
- the six both-endcap rows currently called single traversing;
- the 190 far-only DSA rows;
- the near-only and near-plus-barrel rows recovered only by cosmic.

For the fixed `+z` target sample, explicitly cross-check `inferredSourceSide`
against raw `+z/-z` hit counts and timing. This will determine whether
`farEndcapOnly` means a genuinely downstream-only reconstruction or a direction-side
mistake.

### Phase 3: test the two-algorithm hypothesis

Run the same Step-3 AOD through controlled Step-4 variants. Do not regenerate
simulation or Step 3.

1. **Baseline:** DSA + strict traversing + ordinary cosmic.
2. **No-cosmic:** DSA + strict traversing, with no other setting changed.
3. **Traversing recovery scan:** start from no-cosmic and vary only settings
   that could recover the 93 both-endcap cosmic rows. Keep `BuildTraversingMuon`
   enabled and validate that accepted tracks truly have both-endcap topology.
4. **DSA recovery scan:** test whether near-only and near-plus-barrel cosmic
   rows can be reconstructed by DSA without degrading the established DSA
   population.
5. If necessary, retain a narrowly defined recovery path rather than the
   unrestricted ordinary cosmic collection.

Do not simply route an already reconstructed cosmic track to a DSA label. The
question is whether the DSA fit actually exists and has equal or better
kinematics and fit quality.

### Phase 4: select algorithms by topology only if the evidence supports it

The working hypothesis is:

- `bothEndcaps` -> prefer `throughGoingFit`;
- `nearEndcapOnly` -> prefer `singleLegFit`;
- `nearEndcapAndBarrel` -> tentatively prefer `singleLegFit`, but compare directly
  with recovery/cosmic fits first;
- `farEndcapOnly` -> do not assign automatically until source-side and direction
  inference are validated;
- `unclassified` -> retain as diagnostic/recovery, not a hidden rejection.

Within the same intended topology, continue to use reconstruction-only quality
information: valid directional refit, timing compatibility, station/hit count,
lever arm, and normalized chi2. Topology should determine which fit family is
appropriate, not force acceptance of a failed or pathological fit.

## Validation and decision gates

Use a small populated test first, then the same complete production sample for
all variants. Report denominators as well as selected counts.

For each topology and algorithm, compare:

- total selected rows and events;
- generator-matched yield and unmatched fraction;
- duplicate-component size and source mask;
- valid hit and station distributions;
- timing-direction and directional-refit validity;
- pT, pz, total-p, eta, and charge residuals;
- central-68% width and outlier fraction, not only the median;
- constrained and unconstrained dimuon mass/vertex performance;
- CPU time and failure rate.

Algorithm-removal gates:

1. Do not remove ordinary cosmic until its unique selected rows are either
   recovered by DSA/traversing or explicitly accepted as a quantified loss.
2. Require the recovered rows to have comparable or better resolution and
   outlier rates; equal raw yield is insufficient.
3. Require no unexplained migration in the established DSA and both-endcap
   populations.
4. Validate the result on an LLP-like displaced control. The unconstrained
   `ShiftMuon` must remain available, and fixed-target pointing or target-plane
   compatibility must not become an implicit topology requirement.
5. Keep `unclassified` and invalid-direction rates visible. A simpler label set must
   not hide reconstruction failures.

## Expected outcomes

The likely safe first result is a simpler topology-based public category while
all three reconstruction sources remain internally available.

The stronger two-algorithm outcome is plausible but not established. The
current data already support `bothEndcaps -> throughGoingFit` and
`nearEndcapOnly -> singleLegFit` as the main architecture. The unresolved cosmic
recovery rows and far-only DSA rows determine whether ordinary cosmic can be
removed completely or must survive as a narrowly scoped recovery algorithm.

## Tomorrow's restart checklist

1. Confirm that the topology-diagnostic changes are committed or identify the
   exact dirty-tree diff before further edits.
2. Add the mutually exclusive topology value and `abs(exitZ-entryZ)` without
   changing reconstruction or cleaning.
3. Add duplicate-component source availability diagnostics.
4. Produce focused event lists for:
   - 93 both-endcap cosmic rows;
   - six both-endcap single-traversing rows;
   - 190 far-only DSA rows;
   - 18 near-plus-barrel cosmic rows;
   - rare `unclassified` rows.
5. Inspect several examples from each list, then run the baseline/no-cosmic A/B
   on the same small input.
6. Only after the small A/B is understood, run the complete-sample comparison
   and decide whether to keep, narrow, or remove the ordinary cosmic recovery.
