# Future SHIFT muon and dimuon detector studies

## Scope and coordinate warning

This is the backlog after the Phase-I muon-diagnostics and tracker-combination
pass. The reconstructed target-state pseudorapidity is not a detector-coverage
label. An incoming SHIFT trajectory can have large `abs(eta)` at its target
state and still cross a near endcap followed by DT/RPC barrel layers. Detector
reach must be determined from propagated crossings and recorded measurements,
not from a single eta cut.

Keep detector-only, tracker-only, combined, unconstrained, and target-constrained
hypotheses separate until physics validation shows which one should be used.
Simulation truth is allowed only for diagnostics and performance denominators.

## Bounded probes on 2026-08-14

These tests used v15 plus small ten-event Step-2 inputs.  The 100-event detector
sample was deliberately enriched with files in which v15 contained at least one
`generalTracks` object, so its absolute acceptance fractions must not be
extrapolated to the full sample.

- No tracker information entered the v15 ShiftMuon fit: all 3,896 selected
  muons have invalid tracker and combined-track links and zero attached tracker
  hits.  Only 9 of 9,960 events contain any `generalTracks`, and only 6 of those
  events also contain a ShiftMuon.
- In the enriched 100-event Step-2 sample, primary signal muons made 53 tracker
  SimHits in 2 events.  One inspected leg crossed TID and TEC with 25 SimHits,
  but collision tracking, the stock P5 cosmic CTF pass, and a field-on/wider-TEC
  P5 variant all reconstructed zero tracks for it.  The next tracker prototype
  should therefore be a dedicated field-on TID/TEC or outside-in muon-seeded
  pass, not another global-matcher tuning.
- Of 2,123 v15 ShiftMuons with matched signal DT SimHits, only 13 have attached
  DT hits.  A ten-event `navigationType=Direct` probe recovered four DT hits in
  one DSA candidate that had none with Standard navigation, but also changed
  its momentum and removed the other DSA candidate.  Test Direct navigation as
  a parallel diagnostic collection with truth matching and source-selection
  checks before changing production.
- Primary signal muons in the enriched 100-event sample made ZDC SimHits in 39
  events, HCAL SimHits in 8, and no ECAL, PLT, BHM, or BCM1F SimHits.  ZDC times
  near -468 to -463 ns look potentially useful for direction diagnostics.
  Neither Run-3 QIE10 nor legacy ZDC reconstruction produced RecHits because
  both corresponding ZDC digi collections are empty; ZDC digitization must be
  enabled or repaired before testing reconstructed observables.
- GEM remains an acceptance rather than a fitting problem in this campaign:
  only 5 selected v15 muons have matched signal GEM SimHits and none has an
  attached GEM hit.

## Muon-system reconstruction

- Complete a per-muon funnel for DT, CSC, RPC, and GEM: sensitive-volume
  crossing/SimHit, digi, RecHit, segment, seed, track attachment, refit input,
  accepted measurement, and smoothed-state update.
- For the near-endcap-plus-barrel population, measure the DT loss at each stage.
  Check non-collision timing windows, DT 1D hit pairing, 2D/4D segment building,
  navigation from the endcap into the barrel, and fit compatibility separately.
- Add GEM navigation and measurements to ordinary-cosmic and strict-traversing
  reconstruction. The current DSA builder enables GEM, but the stock cosmic
  builder is configured only for DT, CSC, and RPC. Study GEM-assisted seeds for
  muons with too few CSC segments.
- Compare RPC in three roles: spatial fit measurement, association-only hit,
  and timing/direction measurement. Do not assume that its coarser position
  improves curvature merely because it improves hit efficiency.
- Audit whether CSC/DT track RecHits are segments, component hits, or mixtures,
  and whether correlated measurements are being double counted.
- Validate detector alignment, local-error scaling, timing calibration, and
  chamber-by-chamber residual pulls before changing global fit uncertainties.
- Develop detector-aware outlier handling. The generic hit-chi2 removal test
  worsened closure and must not be re-enabled as a universal cut.

## Tracker and combined fits

- Measure tracker crossing -> tracker SimHit -> cluster -> RecHit -> seed ->
  tracker-track -> muon match -> combined-fit efficiency. Split pixel/strip,
  barrel/endcap, layer count, and source-side/opposite-side legs.
- Validate the initial `generalTracks` plus `GlobalCosmicMuonProducer` test for
  match purity, direction, charge, covariance, and material ordering. Keep the
  combined result diagnostic until it improves q/p and dimuon observables.
- If collision `generalTracks` has poor crossing-conditioned efficiency, test
  cosmic tracker reconstruction and a dedicated outside-in pass seeded from
  the source-facing muon state. Identify the exact prompt-origin or seeding cut
  responsible before widening it.
- Consider separate source-leg and through-going tracker combinations; one bad
  tracker leg must not invalidate a well-measured muon-system trajectory.
- Propagate the combined covariance, not only its central momentum, into the
  target and dimuon fits.

## Momentum, direction, and dimuons

- Continue truth-assisted closure at every measurement surface: predicted,
  updated, filtered, smoothed, and true states, split by detector and transition.
- Verify material is applied exactly once with the correct physical direction.
  Retain detailed/geometry material models as controlled ablations.
- Use DT/CSC timing, and test RPC/GEM/calo timing, to improve direction
  efficiency. Report the ambiguous population rather than forcing a direction.
- Compare dimuon efficiency, pair correctness, mass scale/resolution, vertex
  resolution, and tails for muon-only and tracker-combined inputs.
- Preserve unconstrained and target-constrained dimuon hypotheses. A J/psi mass
  constraint may be a diagnostic but must not calibrate the production momentum.
- Revisit pair ranking only after single-muon covariance and tracker combination
  are validated; otherwise ranking changes can hide a single-track defect.

## ECAL, HCAL, HO, and HF

- Propagate each ShiftMuon through calorimeter cells and save matched RecHit
  energy, time, depth, transverse distance, and neighboring-cell sums.
- Determine crossing-conditioned MIP detection efficiency and timing resolution
  in ECAL and HCAL. Separate normal ionization from bremsstrahlung or showers.
- Test calorimeter information first as identification, direction/timing, noise
  rejection, and material-validation observables. Its granularity and Landau
  fluctuations make it a poor default precision tracking measurement.
- Give HF a dedicated study because many SHIFT paths are forward. Check active
  fiber intersection, long/short-fiber response, timing, and whether an HF MIP
  tag provides an efficiency denominator independent of CSC/RPC reconstruction.
- HO is relevant only for trajectories reaching the outer barrel. Treat it as
  a coarse additional tag rather than a curvature measurement.

## ZDC and other very-forward systems

- Establish whether the campaign geometry, digitization, event content, and
  target-to-CMS trajectory include the ZDC sensitive volumes. Map trajectories
  through the beam aperture and ZDC transverse segmentation before interpreting
  a missing signal.
- If available, test ZDC energy and timing as an upstream tag or time reference.
  Its coarse spatial information should not enter the main track fit without a
  demonstrated covariance-aware benefit.
- Check PPS only geometrically. Its location and beam optics make it unlikely to
  help source-side muons in this setup; do not add reconstruction infrastructure
  until an actual active-volume intersection is demonstrated.

## BRIL detectors

- Inventory which BRIL systems are present in Phase-I geometry/simulation and
  which provide event-by-event data in the CMS stream.
- PLT is the most interesting spatial candidate because it is a silicon
  telescope. Study acceptance, hit availability, alignment, timestamps, and
  whether tracks can be associated to CMS events before considering fit input.
- BHM is conceptually relevant for beam-halo-like trajectories. Test it as a
  tag/timing detector if event-level information is accessible; its coarse
  geometry is unlikely to improve momentum resolution.
- Treat BCM1F and related luminometry/radiation monitors as rate/timing systems
  unless their readout and geometry support unambiguous per-muon association.

## Validation gates for every future test

Use a populated bounded sample first, followed by an unchanged full-production
baseline. Report physical-crossing denominators, reconstruction efficiency,
fake/duplicate rate, q/p and pT/pz scale, central-68 and tail resolution,
direction/charge correctness, target DCA, dimuon mass/vertex performance, fit
failures, and CPU. Split results by topology, reconstruction provenance, source
side, detector composition, and tracker-match status. Never promote a new
hypothesis because it builds or increases raw yield alone.
