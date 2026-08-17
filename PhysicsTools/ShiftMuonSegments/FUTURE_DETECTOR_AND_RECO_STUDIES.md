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
  The later product audit found legacy `ZDCDataFramesSorted` under
  `simHcalUnsuppressedDigis`. The new `shiftZDCReco` path consumes that actual
  Phase-I product; its crossing-conditioned occupancy must now be measured on
  the full sample before attempting a spatial constraint.
- GEM remains an acceptance rather than a fitting problem in this campaign:
  only 5 selected v15 muons have matched signal GEM SimHits and none has an
  attached GEM hit.

## Integrated implementation and probes on 2026-08-17

The `detector_integration_v3` preset is the first production candidate that
contains the whole diagnostic chain in one schema. It adds all pixel and strip
RecHit collections after a detector-only ShiftMuon exists, covariance-selects
unassigned DT segments, retains full ECAL/HCAL/HF/HO RecHits, explicitly
enables the Phase-I ZDC digitizer in a separately named Step-2 stream, and
records propagated calorimeter/ZDC associations and timing hypotheses. Calo
and ZDC measurements remain diagnostic and do not alter the precision fit.

The targeted 30-event probe established several mechanism-level results:

- A tracker-dense event contained 441 pixel and 772 matched-strip RecHits, but
  the matched signal muon had zero pixel and zero strip SimHits. Zero attached
  tracker hits was therefore correct, not an association failure. The new
  `simPixelHits` and `simStripHits` columns make this denominator explicit in
  the large test.
- Removing the redundant global-distance veto admitted three DT segments that
  all pass the covariance estimator. All three are in DT chambers crossed by
  the matched signal muon. In that event the fit input grew from four to seven
  measurements and the relative pT error improved from 30.1% to 25.4%. This is
  encouraging but only one truth-matched case; use
  `nAddedDTTruthChamberMatches` and the resolution histograms to measure purity
  and net performance in the larger test.
- Run-3 DIGI had `doZDCDigi=False`; changing the time phase alone could never
  produce a signal. Enabling the legacy Phase-I digitizer plus a -930 ns phase
  produces 22 ZDC digi frames and 22 reconstructed channels per event. Empty
  frames and pedestal/noise make raw occupancy meaningless, so the new
  threshold, energy, maximum-energy, and timing summaries must be used.
- Propagating the fitted helix to detector cylinders/endcap planes is materially
  better defined than the former infinite straight-line match. The small probe
  still shows substantial accidental HO/ZDC associations. Treat calo direction
  as a test observable, not a production direction decision, until its
  truth-crossing efficiency and false-tag rate are measured.

The completed 10k `detector_integration_v3` study then changed the priorities:

- DT augmentation is populated for about 20% of selected ShiftMuons with high
  truth-chamber purity and is retained in the standard production candidate.
- Tracker attachment is empty because the selected muon-system population and
  tracker-crossing population do not overlap in this sample. None of 3,732
  truth-matched ShiftMuons enters the approximate Phase-I tracker envelope,
  although 5,798 status-1 generated muons do. A tracker-only or dedicated
  outside-in candidate path is needed to recover that population; loosening
  rechit association would only attach unrelated strip occupancy.
- ZDC, GEM, ECAL, HF, PLT, BHM, and BCM1F are removed from the active
  optimization scope. HBHE and HO remain diagnostic timing/tag inputs, with
  HO providing the only clearly truth-correlated reconstructed calo subset.
- The next `dt_hcal_tracker_v1` candidate keeps DT and guarded tracker
  augmentation, but actively associates only HBHE/HO. It records a combined
  muon-system plus HBHE/HO timing hypothesis without allowing it to change the
  canonical direction or momentum before a closure test.

A local two-file Step-4 validation of `dt_hcal_tracker_v1` completed on 20
events. The HO-positive file contained two associated ShiftMuons: adding one HO
centroid increased each timing fit from six to seven measurements, preserved
the muon-only direction in both cases, and changed the direction-separation
Delta chi2 from 78.4 to 80.6 and from 15.5 to 117.4. Canonical momentum,
topology, DT augmentation, and tracker augmentation were event-by-event
identical to the prior reconstruction. The TEA histogrammer also read the new
NanoAOD schema and produced the combined-timing diagnostics successfully.

The subsequent tracker-seeded and detector-aware HCAL prototype established
three additional constraints before a larger production:

- The forward tracker failure was localised to seed-to-CKF growth. The stock
  P5 seed assigns 5 GeV total momentum, which is only about 0.01 GeV in pT at
  the observed extreme pseudorapidity and is rejected by the 0.5 GeV CKF cut.
  A 500 GeV seed with broad covariance, both charge hypotheses, a zero CKF pT
  threshold, and a three-hit minimum recovers one eight-hit track on the known
  25-SimHit signal leg. All eight fitted detector units overlap signal-muon
  SimHits. The fitted momentum remains badly underconstrained (about 2 GeV
  rather than hundreds of GeV), so this track is only a spatial hypothesis.
- A nominal cosmic-global DSA match was nevertheless produced in that event,
  even though no primary signal muon had any DT/CSC/RPC/GEM SimHits. The
  tracker and standalone detector lines were about 480 cm apart. Raw links
  are now recorded separately and must pass a 100 cm symmetric line-distance
  and 0.15 rad unsigned axis-angle gate before tracker/combined diagnostics are
  accepted. Tracker-only tracks are not promoted to ShiftMuons.
- The legacy 40 cm HCAL/HO cell-centre matcher accepted finite HBHE time
  sentinels of -999 ns. Invalid times are now rejected. A detector-aware
  TrackDetectorAssociator path records exact crossed HBHE/HO DetIds, rechits,
  crossed and 3x3 energies, maximum-cell energy/time, and valid-time counts.
  It returned no exact crossings for the two previously loose HO-positive
  tracks, for conventional and explicit source-facing propagation. Therefore
  the old HO associations remain A/B spatial diagnostics and no longer enter
  combined timing until a truth-efficiency and accidental-rate gate supports
  a looser association.

On 3,738 truth-matched ShiftMuons the legacy HO threshold scan confirms that
distance widening is useful only as an efficiency/purity diagnostic: 40 cm
gives 9.7% efficiency and 69.4% purity, 200 cm gives 33.2% and 70.4%, and
400 cm gives 46.6% and 70.3%. HO energy does not improve separation. A loose
time requirement such as greater than -20 ns raises the 40 cm matched-sample
purity to 82.4%, but must not be used before subsystem clock and time-of-flight
calibration. HBHE cannot be rescued this way: 120 cm gives 100% efficiency but
only 3.9% purity, and all four existing 40 cm matches are false with the -999
ns time sentinel. Keep HBHE timing disabled and do not widen its raw-distance
gate.

The `dt_hcal_tracker_seeded_v1` variant runs the retained DT augmentation, the
forward P5 tracker prototype, and HCAL/HO diagnostics together. Its intended
next test is bounded: measure tracker seed-to-fit truth purity, raw-to-gated
tracker/muon link purity, exact/loose HO efficiency and fakes, valid timing
fractions, unchanged canonical momentum/dimuon output, and CPU before any
tracker or calorimeter information becomes canonical.

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
- PPS is upstream of the SHIFT target and is outside the physical flight path.
  It is excluded from reconstruction and from the detector-study backlog.

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
