# SHIFT common-vertex classifier study

Started 2026-10-09. Owner: this task log for standalone classifier experiments.
The canonical analysis contract remains in `../../../../SHIFT_ANALYSIS.md`;
reconstruction findings remain in `SHIFT_RECONSTRUCTION.md`.

## Aim and current status

The first seven-model pilot and fresh 26-fit SM comparison are complete.
The uniformity BDT is the strongest current candidate for further study, but
no model passes the required independence gates. Quantitative results follow.

Select two distinct genuine muons from a common physical vertex. Preserve genuine
DY and genuine QCD vertices alongside J/psi. Reject accidental different-vertex
pairs and duplicated muon legs. This is an exploratory SM-only study; no final
selection, normalization, collision-data access or physics readiness is claimed.

Only reconstructed collision-data quantities may be classifier inputs or
adversary targets. Truth is permitted only in the separate training-label
construction. Process, source-bin, event identity, sampling/generator weights,
truth matches, simulated momenta, true vertices and ancestry never enter X.
MC sampling weights are retained separately for diagnostics/bookkeeping.

## Design

- Compare shallow quality-only and expanded BDTs. The expanded inputs include
  reconstructed single-muon reference points, dimuon vertex, four-vectors and
  opening angle, as explicitly requested. Their correlations are measured.
- The explicit invariant-mass column is an adversary/validation target, not a
  predictor. Four-vectors nevertheless encode mass; hiding the mass column
  does not establish independence.
- A literal neural adversary needs a differentiable neural classifier. Train
  it separately from hard decision trees. Compare an identical unconstrained
  neural baseline, several adversarial strengths, and (if statistics permit)
  a BDT distilled from the adversarial model. Distillation must be validated
  again; it need not preserve decorrelation.
- Protect the joint reconstructed mass, vertex z and transverse vertex radius
  distribution separately within each truth label. The adversary only sees
  the scalar classifier score. These are displacement proxies, not lifetime.
- Use event-grouped training/validation/test partitions. Fit preprocessing and
  nuisance bins on training only; choose score thresholds on validation only.
  Report test counts, intervals, per-process efficiency, nonlinear dependence,
  joint efficiency maps and regions with insufficient support.

## Labels and sample restrictions

Use corrected V10 SM Nano products from a frozen explicit published-file
inventory. Never combine V8/V10 versions. V10 adds the previously missing
muon-producing hadron decays. A published subset is adequate for machinery
tests but is not a complete or unbiased production-yield estimate.

Use `ShiftMuon_hitGenPartIdx`, purity >=0.75 and >=3 matched CSC layers. Require
two status-1 generator muons. Independent hit matches avoid defining labels
through fitted angles. Distinct generator muons within 0.01 cm in production
position are positive; separated by >=0.1 cm are negative; the gap, unmatched
or ambiguous states are unknown and excluded from supervised training. A
reused generator/SimTrack muon is negative. The position tolerance reflects
the stored-coordinate precision and is a provisional label convention.
Different immediate mothers are insufficient: DY muon-copy chains can differ
while referring to the same physical vertex. Authentic common-vertex QCD pairs
must be positive. Match completeness and label-tolerance sensitivity need audit.

## Physical limits

SM samples cannot establish efficiency independence versus hypothetical BSM
mass or lifetime. Prompt SM support can test only reconstructed displacement
dependence. Decorrelation is a training objective, never an automatic guarantee.
Flatness of the added classifier cannot remove mass/lifetime dependence of
fixed transport, electronics, trigger or reconstruction. Transfer of a J/psi
normalization requires the relevant acceptance/efficiency ratio or proof that
its factors are shared. Unsupported regions are unvalidated.

Existing Nano rows have already passed opposite-sign/common-side/DCA/origin
compatibility conditions and a common-line fit with a transverse beam-line
prior, followed by greedy pair ranking. This pilot measures incremental
classifier performance conditional on those choices. Single-muon vx/vy/vz are
track PCA/reference coordinates, not measurements of a production vertex.
The current common-line pair vertex is not a fully prior-free fit.

## Activity log

### 2026-10-09: initial inspection and preparation

- Verified CMSSW branch `shift-muon-segments-pre4`. Existing untracked
  `IOMC/ShiftDarkPhoton` and `IOMC/ShiftMuonDecays` are preserved.
- Created an independent Python directory under PhysicsTools. No BuildFile,
  shared CMSSW rebuild, campaign edits, production submissions or data reads.
- Runtime available read-only from
  `/cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh`:
  Python 3.12.11, numpy 2.1.3, sklearn 1.5.2, torch 2.7.1, uproot 5.6.0.
- Initial read-only 10-file-per-bin inspection: 190 files, 3635 events, 390
  retained pairs. Hit-label candidates: 281 common vertices, 47 different
  vertices, one duplicated leg, 61 unmatched. J/psi and DY supply accidental
  negatives too; QCD alone has only eight retained pairs at this size.
- Started a bounded 40-file-per-bin extraction, using the first explicit
  published files from `validation/histogram_rerun_20261009/histogram_inputs.txt`.
  Record exact paths, source counts, feature contract and output digests in
  the artifact manifest. No experiment geometry payload is copied here.
- The default Uproot reader stalled while reading the mounted EOS files.
  Stopped only those two study-owned read processes. An explicit local
  `MultithreadedFileSource` with one reader opens the file in 0.002 seconds
  and reads all requested columns in 0.46 seconds. The initial manifest-only
  attempt remains at `artifacts/sm_v10_40files/`; the actual extraction uses
  `artifacts/sm_v10_40files_local/`.
- Fourteen contract tests pass. A separate synthetic neural gradient check
  confirmed that the adversarial objective changes scores and trades some
  separation for reduced adversary information. Synthetic performance is not
  reported as SM performance.
- Training weights equalize the two labels and available processes within
  each label. The illustrative threshold targets 90% balanced validation
  positive efficiency. Raw held-out efficiencies are reported separately.
  Sampling-weight diagnostics retain W/p but omit source cross sections and
  full-GEN exposure: integrated values are not physical production predictions.
  Efficiency Wilson intervals are pair-level approximations; AUC intervals
  resample complete events. No rate or normalization is inferred here.

### 2026-10-09: bounded SM training results

Completed extraction from exactly 760 corrected V10 Nano files (40/bin across
19 bins): 14,327 events and 1,551 retained dimuons. Labels are 1,029 common
vertices, 220 different vertices, three duplicate legs and 299 unknown matches
(19.3%). There are no rows in the label-position tolerance gap in this sample.
J/psi contributes 641 positives/129 negatives/194 unknown; DY 388/69/97; QCD
0/25/8. This supports a pair-authenticity comparison, not process classification.

The fixed event split has 620 positives/136 negatives for training, 203/49 for
validation and 206/38 for test. Test negatives comprise 22 from J/psi, 12 from
DY and only four from QCD. Unknown labels are not trained or used to claim
efficiency; all 1,551 reconstructed rows receive saved data-only predictions.

The following raw held-out results use each model's threshold chosen solely
to retain approximately 90% of process-balanced validation positives. Their
test efficiencies differ; these are not rejection estimates at identical test
efficiency. Distance correlation below is the empirical nonlinear dependence
of score on the joint reconstructed mass, vertex z and transverse radius for
genuine pairs. It has finite-sample bias and is a diagnostic, not a gate alone.

| Model | Test AUC (95% event bootstrap interval) | Genuine retained | Accidental rejected | Genuine joint dCor |
| --- | --- | --- | --- | --- |
| Quality BDT | 0.713 [0.619, 0.789] | 83.0% | 36.8% | 0.638 |
| Expanded BDT | 0.837 [0.750, 0.895] | 86.4% | 65.8% | 0.661 |
| Expanded NN, lambda=0 | 0.785 [0.703, 0.858] | 92.2% | 39.5% | 0.551 |
| Expanded NN, lambda=0.1 | 0.758 [0.678, 0.841] | 91.3% | 39.5% | 0.481 |
| Expanded NN, lambda=0.5 | 0.705 [0.617, 0.790] | 91.7% | 31.6% | 0.334 |
| Expanded NN, lambda=2 | 0.595 [0.487, 0.676] | 86.9% | 23.7% | 0.216 |
| BDT distilled from lambda=0.5 NN | 0.825 [0.748, 0.886] | 85.0% | 63.2% | 0.626 |

The expanded BDT's largest split importances are pair vertex z (33.6%), pair
pz (16.4%), the lower-pT muon's pz (12.9%) and energy (12.5%). This suggests
that its improved separation exploits geometry and momentum information;
split importance alone is not a causal attribution.

Mass acceptance shows the problem directly. In training-defined mass thirds
(below 3.069 GeV, 3.069--5.459 GeV, above 5.459 GeV), the expanded BDT retains
71.8%, 92.1%, 95.8% of genuine test pairs. Lambda=0.5 improves these to 90.1%,
92.1%, 93.1%, while its vertex-z efficiencies remain 82.9%, 95.2%, 97.3% and
radius efficiencies 98.6%, 88.9%, 87.8%. Lambda=2 still retains 97.1%, 84.1%,
79.7% across radius thirds. A lower correlation statistic does not establish
flat acceptance at the chosen cut.

Authentic DY/J/psi efficiencies are 82.1%/89.3% for the expanded BDT and
89.3%/93.4% for lambda=0.5. The latter retains 75/84 DY and 114/122 J/psi
positives but rejects only 12/38 accidental pairs. The expanded BDT's zero
accepted QCD pairs means 0/4, with Wilson 95% acceptance interval [0%,49%];
it cannot be reported as demonstrated strong QCD rejection. No genuine QCD
pair was available to validate its preservation.

Every one of the 27 joint mass/z/radius cells contains fewer than 20 genuine
test pairs. Joint acceptance is unvalidated. At lambda=0.5 the negative
adversary's validation cross-entropy is 3.864 versus marginal entropy 2.893;
the positive values are 3.198 versus 3.116. This reveals weak adversary
generalization, rather than proving independence from a fooled adversary.
The distilled BDT regains unwanted dependence (0.334 to 0.626); do not use
distillation as a shortcut to a validated decorrelated BDT.

Artifacts are under `artifacts/pilot_seed42/`: seven model bundles, all-row
predictions, fixed event identities/splits, full training histories, feature
correlations, per-process/per-bin efficiencies, bootstrap intervals and
`metrics.json`. Training took 10.7 seconds on one CPU thread using the pinned
LCG runtime. Extraction took approximately seven minutes with one EOS reader.
The exact measured extraction time is in the input manifest.

Figures (visually checked):

- [Overview PDF](artifacts/pilot_seed42/figures/overview.pdf) and
  [PNG](artifacts/pilot_seed42/figures/overview.png).
- [Joint efficiency maps](artifacts/pilot_seed42/figures/joint_efficiency.pdf),
  with combined and per-process mass/z maps; radius is integrated out in these
  display maps. Full 3D cell counts remain in `metrics.json`.

### 2026-10-09: scoring contract and final checks

All seven models passed physical truth removal on nine SM ROOT files, covering
11 retained pairs (four J/psi, four DY, three QCD). The stripped ordinary TTrees
contain exactly the reconstructed allowlist: no GenPart, hit-truth, simulation
diagnostics or generator weights. Feature arrays, auxiliary reconstructed
values, identities, prediction scores and threshold decisions are identical
before/after stripping; the maximum absolute score difference is zero for
every model. The valid receipt is
[`truth_removed_ttree/report.json`](artifacts/pilot_seed42/truth_removed_ttree/report.json).
This proves standalone inference independence on this MC fixture, not complete
collision-data readiness or upstream reconstruction independence.

The initial writer generated extra per-column jagged-array counters and the
exact-schema audit correctly stopped. Two attempted RNTuple workarounds exposed
the older Uproot record-input and read API differences. The final writer groups
each reconstructed collection as a jagged record, producing the ordinary Nano
shared counters and exact original leaf names. The failed fixture attempts
remain under `truth_removed/`, `truth_removed_rntuple/` and
`truth_removed_validated/` for provenance; none has a valid report. A synthetic
TTree round-trip regression test now covers this exact serialization boundary.

Sixteen contract tests pass, including poisoned-truth access, pair ordering,
labels, event grouping, nonfinite match purity and exact ROOT field/value
round-trip. All Python sources compile and the runner's shell syntax passes.
The study stays on the required CMSSW branch. Only this standalone directory
and the classifier cross-reference in `SHIFT_ANALYSIS.md` changed; existing
untracked generator packages and all campaign inputs/outputs are preserved.
No shared CMSSW build, production configuration mutation or collision-data
inspection occurred. Model/table/plot artifacts are ignored by Git; nothing
was staged, committed or pushed.

### 2026-10-09: second phase authorized and prepared

The user requested continuation. The frozen second-phase inventory contains
4,135 completed V10 files: 3,095 QCD files and 80 additional files in each of
the 13 J/psi/DY bins (1,040 SM control files). All 760 pilot source files are
excluded. Path-hash ordering spreads the control subset over the published
inventory without inspecting candidate observables or training labels.
The dataset also refuses any physical event overlap with the pilot before
marking completion. Campaign production files and configurations are read-only.

`phase2_data.py` uses one ROOT reader and reads full allowed collections only
for events with reconstructed dimuon rows. It preserves original event indices
and total source event counts. ROOT and Uproot feature/label equality passed on
19 representative files (50 pairs), with an independent three-file/four-pair
check including diagnostics and weights. Nonempty source checkpoints are
checksum-verified and sharded; zero-pair counts are checkpointed in the manifest.
The active dataset is `artifacts/sm_v10_phase2_fresh/`.

Added two independently reviewed reconstructed-only methods:

- `uniform_bdt.py`: a single-operating-point neighborhood-uniformity AdaBoost
  variant inspired by uBoost, rather than the publication's full ensemble.
  Underaccepted genuine-vertex neighborhoods in joint mass/z/radius get more
  training weight. Strength zero is the paired equal-class AdaBoost reference.
  Normalized margins preserve ranking without sigmoid saturation. No nuisance
  values or labels are retained in the inference object.
- `disco.py`: the same small neural architecture with class-conditional,
  weighted joint distance-correlation-squared regularization. Training-only
  nuisance scaling and multivariate distances protect joint dependence.
  Tests verify weighted empirical identities, finite nonzero gradients,
  invariance properties and an XOR example where marginal checks miss the
  joint dependence. Constant scores are flagged, not called independence.

Forty-seven contract/mechanism tests pass. Independent reviews and numerical
gradient checks found no blocking implementation errors. All methods still
require held-out efficiency validation; no automatic flatness claim is made.

`phase2_train.py` declares 26 candidates before test performance is computed:
five BDT feature-group ablations, uniform strengths 0/1/3, DisCo strengths
0/0.3/1/3 and an adversarial lambda=0.5 reference, each at fit seeds 71/72.
Both seeds share the fixed event split seeded 71. Thresholds target 90%
smoothed process-balanced validation-positive acceptance. Smoothing avoids
giving one rare genuine QCD pair a third of all positive training weight.
Representative choices use only validation rejection, the worst supported
marginal efficiency variation including per-process checks, process efficiency
spread and joint score distance correlation. Joint cut maps remain a separate
validation requirement. All thresholds/choices are written to a lock receipt
before test scores are computed. Test class counts are inspected for adequacy;
no test score/performance is used for training or selection.

The first label audit is in `artifacts/label_audit_phase1/`: unknown labels
vary with reconstructed coverage, quality and geometry. Below 18 hits on the
less-measured muon, 27.7% of pairs are unknown versus 15.0% at >=58 hits; vertex-z
quartiles span 26.8% to 13.1% unknown. This population bias is unresolved and
must accompany efficiency claims. `label_audit.py` now also audits match
purity/layers and position-tolerance sensitivity using separate label-only
diagnostics. Missing mother IDs are recorded as unknown, never distinct mothers.

### 2026-10-09: fresh comparison results

Completed 4,135 source files / 75,957 events / 3,494 reconstructed pairs in
747.8 seconds with one reader. The independent freshness audit verifies zero
shared files and reconstructed event identities with the pilot. Labels are
2,047 common vertices, 780 different vertices, ten duplicate legs and 657
unknown matches. QCD contributes 339 negatives, one positive and 76 unknown;
DY 715/140/179; J/psi 1,331/311/402 (positive/negative/unknown).

Fixed partitions: train 1,231 positives / 462 negatives, validation 420/153,
test 396/175. Test includes 122 genuine DY and 274 J/psi, 82 QCD negatives,
93 other SM negatives and 133 additional unknown pairs. The sole genuine QCD
pair is outside test; preservation of genuine QCD remains unvalidated.

All 26 fits completed in 49.9 seconds. Thresholds and representative model
names were locked before test scores/performance were computed. The selected
representatives below use validation-only choices. Raw test counts are shown,
with no cross-section or physical-rate claim.

| Representative | Test AUC (95% event bootstrap interval) | Genuine kept | Accidentals rejected | Genuine joint dCor |
| --- | --- | --- | --- | --- |
| Quality BDT, seed 71 | 0.819 [0.782, 0.855] | 86.9% | 60.6% | 0.619 |
| Expanded BDT, seed 71 | 0.926 [0.903, 0.951] | 90.2% | 77.1% | 0.704 |
| Quality + geometry BDT, seed 71 | 0.879 [0.848, 0.908] | 89.9% | 68.0% | 0.717 |
| Uniformity BDT, strength 3, seed 71 | 0.856 [0.821, 0.887] | 90.2% | 58.9% | 0.271 |
| Joint DisCo NN, strength 3, seed 71 | 0.786 [0.748, 0.828] | 90.4% | 47.4% | 0.250 |
| Adversarial NN, lambda 0.5, seed 71 | 0.851 [0.814, 0.879] | 86.6% | 65.7% | 0.342 |

Uniformity and expanded BDT retain exactly the same 357/396 common vertices,
allowing a direct operating-point comparison: 103/175 versus 135/175
accidentals rejected. Uniformity BDT retains 113/122 DY (92.6%, Wilson 95%
interval 86.6--96.1%) and 244/274 J/psi (89.1%, 84.8--92.2%). It rejects 56/82
QCD negatives (68.3%, 57.6--77.4%), compared with 73/82 (89.0%, 80.4--94.1%)
for expanded BDT. QCD counts now support a useful conditional rejection
estimate, unlike the original four-pair test.

Both geometry and momentum add separation: geometry raises AUC 0.819 to
0.879, momenta to 0.896, while opening angle alone gives 0.816. Expanded input
gives 0.926. These comparisons do not demonstrate independence; quality-only
inputs remain correlated too. Uniformity seeds 71/72 give identical results.
Other BDT seeds have similar AUCs; neural efficiencies/rejection vary. This
tests optimization stability on one partition, not independent-sample stability.

The uniformity BDT's combined acceptance profile is:

| Reconstructed quantity | Training-defined boundaries | Genuine test efficiency |
| --- | --- | --- |
| Mass | 2.896, 5.359 GeV | 88.8%, 90.0%, 91.7% |
| Vertex z | 12,827, 14,783 cm | 83.5%, 93.3%, 92.7% |
| Vertex radius | 13.449, 32.196 cm | 94.6%, 93.2%, 81.7% |

Mass acceptance is much flatter than expanded BDT's 82.4%, 90.7%, 97.5%.
Position acceptance still fails: lowest/highest radius-bin confidence intervals
do not overlap. Only eight of 27 three-dimensional cells reach 20 genuine
test pairs; their efficiency spans 80--100%. Another 191/396 genuine pairs
occupy unsupported cells. Low-mass DY lacks adequate test support: its thirds
retain 10/13, 23/24 and 80/85, while validation efficiencies were 69.6%, 89.7%,
95.7%. Combined SM mass flatness hides this process dependence. Joint cut
efficiency and process transfer remain unvalidated for every model.

All 26 models pass physical truth removal on nine files (21 reconstructed
pairs): scores and decisions are exactly identical without generator,
hit-truth, simulation diagnostic or generator-weight branches. No new truth
enters X or decorrelation targets. The valid receipt is
`artifacts/phase2_seed71_72/truth_removed/report.json`.

The fresh label audit reproduces every nominal label. All 657 unknown pairs
have an unavailable hit generator index: 588 one-leg failures, 69 two-leg
failures. Looser purity/layer thresholds recover none. Raising purity to 0.9
discards another 1,595 labelled pairs, including 1,189 positives, and raises
the unknown fraction to 64.5%. No known label flips between positive and
negative. Tested vertex
tolerances (common 0.001--0.05 cm; separated 0.1--1 cm) change zero labels.
This checks convention stability among available matches, not match accuracy.
Of 2,047 positives, 824 have distinct immediate mothers: mother equality
would wrongly remove 40.3% of positives.

Matching bias persists: low/high hit-count groups have 24.6%/13.1% unknown,
low/high track-chi2 groups 13.2%/32.2%, z groups 24.9%/13.4%, radius groups
14.8%/24.9%. Boundaries and process mixture differ from the pilot; compare
trends, not identical cuts. Measured efficiencies describe reliably matched
candidates, not every reconstructed candidate. The 18.8% missing fraction is
not itself a physics uncertainty.

Added `acceptance_audit.py` for per-process cut maps, unsupported populations
and unknown-score retention. Uniformity BDT keeps 100/133 unknown test pairs.
Allowing every possible unknown-label assignment gives count-only genuine
efficiency bounds 83.2--92.1% and accidental-rejection bounds 37.5--65.4% in
this unweighted sample. These are identification limits conditional on correct
known labels, not a physical uncertainty or correction. Independent exhaustive
tests cover 1,555 synthetic label/decision patterns. They exposed and verified
the fix for absent-class corner cases; the original diagnostic receipt remains
for provenance and `acceptance_audit_v2/` contains the final report.

Reviewable artifacts:

- [Full metrics](artifacts/phase2_seed71_72/metrics.json), protocol,
  validation candidates and `selection_locked.json`, alongside 26 model bundles
  and all-row predictions.
- [Overview](artifacts/phase2_seed71_72/figures/overview.pdf) and
  [joint maps](artifacts/phase2_seed71_72/figures/joint_efficiency.pdf).
- [Uniformity acceptance profiles](artifacts/phase2_seed71_72/acceptance_audit_v2/acceptance_profiles.pdf)
  and [per-process/unknown audit](artifacts/phase2_seed71_72/acceptance_audit_v2/acceptance_audit.json).
- [Label audit](artifacts/label_audit_phase2/label_audit.json).

Keep uniformity strength 3 as the next benchmark, not a frozen physics cut.
It improves flatness with useful rejection, but radius, process and support
failures prevent universal N_eff. SM-only tests still cannot establish lifetime
independence; detector/trigger/reconstruction acceptance remain upstream gates.

Final checks: 58 tests pass in the pinned runtime. All Python sources compile;
both reproduction runners pass shell syntax and environment initialization.
The latter found that the LCG setup expects some optional variables to be
unset; runners now enable strict unset-variable checks after sourcing it.
Plots were visually inspected, with explicit units and grey unsupported bins.
The current CMSSW branch remains `shift-muon-segments-pre4`. Model/data artifacts
are ignored, every model's deployment flag remains false, and no shared build,
campaign configuration change, collision-data access, commit or push occurred.

## References

- [Learning to Pivot with Adversarial Networks](https://arxiv.org/abs/1611.01046):
  adversarial optimization of predictive independence.
- [uBoost](https://arxiv.org/abs/1305.7248): uniform-efficiency boosted trees;
  a possible direct BDT comparison for the next phase.
- [DisCo](https://arxiv.org/abs/2001.05310): distance correlation for nonlinear
  dependence checks and an alternative differentiable regularizer.

## Remaining work

### 2026-10-09: committing and preparing the first signal grid

The user authorized commits and requested readiness for a 3x3 mass/epsilon
grid spanning prompt, metre-scale and tens-of-metres decays. Source, tests,
reproduction runners and this log are durable and belong in Git. Generated
tables/models, checksummed receipts and validation plots remain outside Git;
retain them until independent displaced-grid closure. Failed fixture attempts,
manifest-only extraction attempts and canaries are future cleanup candidates,
not production inputs. No deletion occurs while these replacement/closure
conditions are unmet. The exact cleanup/retention plan will accompany the grid.

The frozen uniformity-strength-3 model is the first grid benchmark. Its current
radius/process failures remain explicit; do not retune it on the consumed SM
test population. Signal scoring needs an explicit simulation manifest rather
than weakening the SM reader's V10/process allowlist. Mass/epsilon, proper
decay length and expected lab flight belong to point metadata and acceptance
audits, never classifier inputs. The physical signal prescription and grid
configuration are owned by `SHIFT_BSM.md` and the workflow repository.

1. Freeze both completed comparisons. Further optimization needs fresh reserved
   tests with additional low-mass genuine DY and extreme reconstructed z/radius
   support. Do not tune neighborhoods or thresholds against the phase-two test
   failures and report the same test as independent.
2. Match completeness remains unresolved despite stable label tolerances.
   Inspect unavailable hit-index mechanisms with independent MC association
   evidence. Extending truth-based labels must never add truth to predictors.
   Label repair may need preserved intermediate truth products and bounded
   validation, rather than weaker purity or guessing labels from fitted angles.
3. Compare prespecified genuine-neighborhood sizes, multiple efficiency targets
   or stronger joint constraints. Require per-process and joint cut acceptance
   within the intended normalization precision; lower score correlation alone
   is insufficient. Retain a transparent compatibility-statistic baseline.
4. Validate displaced signals and held-out mass/lifetime grid points, together
   with upstream pair retention and fixed detector/trigger/reco acceptance.
   This is mandatory before interpreting the selection as common to J/psi and
   BSM or transferring N_eff.
5. Integrate canonical blinding and data-path gates before any collision-data
   run. No final physics threshold is frozen; all model deployment flags stay
   false. Model fitting uses SM simulation; signal MC supplies held-out
   transfer checks, without retuning the model or threshold.

## Reproduction

### 2026-10-09: real-signal scoring canary and evidence retention

The reusable classifier study was committed as `f7a521c6bd4`. A separate
`signal_grid_audit.py` evaluates the locked SM benchmark on explicit, hashed
signal Nano inputs without retraining or choosing a threshold on signal.
It reads only reconstructed columns for prediction, then reads truth for
labels and parent-ancestry auditing. Removing all truth branches physically
preserves exact scores and decisions on both tested signal files.

The completed two-event prompt 15 GeV replay contains two authentic A' dimuon
pairs, both passing `uniform_bdt_3p0_s71` at the locked threshold
`0.4403188643890805`. The two-event displaced replay has no retained dimuon
pairs upstream (zero or one ShiftMuon per event); conditional classifier
efficiency is undefined. These tiny files validate scoring and denominator
bookkeeping only. Neither sample demonstrates mass/lifetime closure.

Receipts are in `artifacts/signal_canary_v1/`, with the frozen source inventory
in `artifacts/signal_canary_v1_inputs.json`. The complete suite has 66 passing
tests. Nano-event fractions include zero-pair events; full-GEN and recording
efficiencies remain undefined without a reconciled production ledger.

Keep source, tests, runners and this log in Git. Keep frozen tables, split
identities, model/threshold locks, predictions, source inventories and final
audits outside Git as reproducibility evidence. Preserve existing failed-run
logs until their replacement and failure explanation are recorded. Only
truth-stripped fixture ROOT files, caches and superseded scratch plots may be
removed after successful grid validation and after confirming they are not
needed by a reproduction command. No study or campaign evidence was deleted.

The first proposed signal grid uses 15, 30 and 50 GeV, each with prompt,
approximately 5 m and approximately 30 m **mean lab flight**. Its physical
epsilon values and proper decay lengths are derived consistently from the
native width authority. The signal model and grid execution plan are owned by
the workspace `SHIFT_BSM.md` and the workflow generation instructions.

```bash
python3 signal_grid_audit.py --manifest artifacts/signal_canary_v1_inputs.json \
  --results artifacts/phase2_seed71_72 --model uniform_bdt_3p0_s71 \
  --output artifacts/signal_canary_reproduction
```

### 2026-10-09: nine-point GEN qualification and second signal canary

The first 3 by 3 grid is frozen in the workspace
`validation/bsm_grid_20261009/grid_pilot/plan.json`. Columns refer to prompt,
approximately 5 m and approximately 30 m **mean lab flight**, not proper decay
lengths. This assumption was stated while the optional parameterization
question remained unanswered. A strict Cartesian epsilon grid is also
supported by the preparer's `--epsilons` option.

| Mass (GeV) | Prompt epsilon | 5 m epsilon | 30 m epsilon | Proper mean length for 5 m / 30 m (mm) |
| --- | --- | --- | --- | --- |
| 15 | 1e-3 | 1.0955803944e-7 | 4.4726882308e-8 | 62.178 / 373.071 |
| 30 | 1e-3 | 6.9521347018e-8 | 2.8381971071e-8 | 73.529 / 441.174 |
| 50 | 1e-3 | 4.7888347820e-8 | 1.9550336131e-8 | 78.087 / 468.521 |

Measured unfiltered reference boosts are 80.4137 +/- 3.4871, 68.0003 +/- 3.2296
and 64.0313 +/- 2.2056 (standard errors). Native physical widths, independently
checked partial widths, physical branching fractions and resolvable proposal
support passed for all nine points. The physical exponential is retained;
the approximate mean flights inherit the finite reference uncertainty.

All nine 100-event GEN pilots passed in the verified native CMSSW context:
900 distinct persisted event identities, unit weights and matching complete
contracts; 5,402 technical Pythia trials, 900 selected and 900 accepted events.
The qualification took 153.9 seconds and the whole preparation/evidence
directory uses about 30 MB. No outside-cylinder decay was observed in these
900 draws; the full unfiltered denominator and escaped-outcome retention
policy are unchanged. Historical long-lived references retain such outcomes.
The native attempt followed an independently preserved sandbox startup crash
before any events. Neither attempt rebuilt the release or altered live jobs.

The authoritative receipt is `grid_pilot/gen_qualification.json`, with all
per-point payload/audit hashes, source equality and timing. Native seeds/run
numbers are 24694357 through 24694365. `grid/` is the initial preflight and
`grid_ready/` preserves the failed sandbox qualification; keep these receipts.
The added finite-sample exponential mean gate also passed on all 900 stored
proper lengths in `grid_pilot/independent_lifetime_check.json`, preserving the
original validation reports. This GEN-only qualification launched no detector
jobs; it is separate from the workspace's bounded detector production plan.

The additional `artifacts/signal_canary_v2/` audit checked six existing, fully
audited Nano events. Four 15 GeV, epsilon=1e-8 events have no ShiftMuon rows or
pairs, so conditional classifier efficiency is undefined. Two 50 GeV,
epsilon=3e-8 events contain one authentic direct A' pair. It passes with score
0.51874605 at the fixed 0.44031886 threshold (1/1 pairs, Wilson 95% interval
20.7-100%; 1/2 recorded events). Exact scores and decisions survive physical
truth removal for both files. These counts validate machinery only.

Durable code is committed on the required branches: CMSSW classifier
`f7a521c6bd4` and signal audit `31e6942eeb0`, signal adapter `5dc3b2cbc50`,
natural-hadron decayer `82134fbfe32`; workflow tools `2c8bec1` and finite-sample
GEN gate `5655e88`; TEA histogram worker `2eacfaa`. Focused validation passed:
66 classifier, 63 dark-photon/model/grid/GEN-gate, 13 proposal-ledger, 13 common
worker and seven histogram-worker tests (162 total). No push was performed.
Actively changing plotting, merging and retirement work is preserved for its
own review. Production cards, ROOT files, plots, models and runtime records
remain outside Git.

AFS briefly filled and blocked a commit. Removed only 317 untracked,
regenerable bytecode cache files (3,460,053 bytes); the receipt is retained in
`validation/bsm_grid_20261009/bytecode_cleanup.json`. An attempted test-fixture
archive stopped on an unexpected nested directory before any deletion; all
original fixtures and other study evidence remain. One-time reference and
qualification probes are retained as hashed evidence copies under validation;
their scratch copies are retired only after byte-identical verification.

Mass/lifetime transfer, full detector/recording efficiency and a common N_eff
remain unvalidated. The next step is to score grid Nano outputs with this same
locked benchmark, including zero-pair events and reconciled full-GEN exposure,
then assess supported conditional acceptance at each mass/lifetime separately.

Run from this directory. The environment is read-only; packages are not
installed into or built with the shared CMSSW release. Generated tables,
models and plots remain under ignored `artifacts/`.

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 extract.py \
  --inventory /afs/cern.ch/work/j/jniedzie/private/shift_cmssw/validation/histogram_rerun_20261009/histogram_inputs.txt \
  --output artifacts/sm_v10_40files_local --files-per-bin 40
bash run_smoke.sh artifacts/sm_v10_40files_local artifacts/pilot_seed42
```

Existing artifact directories are refused to avoid overwriting evidence; use
a new output name for another run. `inputs.json` records the frozen input
inventory and table digests. `metrics.json`, saved split identities, model
bundles, training histories and all-row scores make the comparison reviewable.
`inference.predict` takes only X and the exact feature-name contract, with no
labels or nuisance targets. Model pickle files are local trusted study outputs.

The fresh comparison is reproduced using the completed second dataset:

```bash
python3 phase2_data.py \
  --inventory /afs/cern.ch/work/j/jniedzie/private/shift_cmssw/validation/histogram_rerun_20261009/histogram_inputs.txt \
  --exclude-dataset artifacts/sm_v10_40files_local \
  --output artifacts/sm_v10_phase2_fresh --sm-files-per-bin 80 --qcd-files-per-bin 650
bash run_phase2.sh artifacts/sm_v10_phase2_fresh artifacts/phase2_seed71_72
python3 label_audit.py --dataset artifacts/sm_v10_phase2_fresh \
  --output artifacts/label_audit_phase2 \
  --reference-dataset artifacts/sm_v10_40files_local --require-disjoint
```

These output names already exist and are protected. Use new output names for
reproduction with the frozen inventory and implementation; reuse the saved
complete dataset for model reproduction. Existing test predictions remain
historical evidence, not fresh data for further tuning.
