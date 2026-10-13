# SHIFT dimuon vertex classifier

Started 2026-10-09. Owner: this task log for standalone classifier experiments.
The canonical analysis contract remains in `../../../../SHIFT_ANALYSIS.md`;
reconstruction findings remain in `SHIFT_RECONSTRUCTION.md`.

## Aim and current status

The first seven-model pilot and fresh 26-fit SM comparison are complete.
The original uniformity BDT remains the frozen ntuple benchmark. A subsequent
25-fit capacity/uniformity/data-size study improves rejection substantially,
but no candidate passes all acceptance gates. Quantitative results follow.
Provisional scores, histograms and held-out evaluation are complete for 510
simulation Nano files spanning all 19 SM bins and nine signal points (October
12 first look below). No score cut was applied. A mass/lifetime-independent final classifier selection is not
validated. N_eff is the effective MC normalization fitted to collision data;
it is separate from the classifier efficiency/AUC reported in this study.

The October 13 follow-up completes 70 neural fits, a BDT overtraining audit
and category-resolved diagnostics. A penalized DNN improves genuine Near-Both
retention at approximately the same rejection as the flatter BDT. It does not
improve pooled rejection substantially, and category/seed dependence remains.
All new results and PDF figures are inside
`tea_shift_cmssw/plots/DNN_vertex_20261013/`. The original deployed BDT remains
unchanged. See the new activity entry and the study README for exact results.

Select two distinct genuine muons from a common physical vertex. Preserve genuine
DY and genuine QCD vertices alongside J/psi. Reject accidental different-vertex
pairs and duplicated muon legs. Training uses SM samples; the first evaluation
also includes the small signal grid. No final
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

### 2026-10-13: overtraining audit, vertex categories and penalized DNNs

The strong BDT's repeated energy boundary is already present in early trees.
At 50 trees the training/validation AUC is 0.8925/0.8853; at 400 it is
0.9083/0.8894, while the development-test AUC improves from 0.8979 to 0.9029.
This is a modest generalization gap, not evidence of a numerical bug or
runaway overtraining. The validation-selected 150-tree prefix gains only
1.32 percentage points in fake rejection over 400 trees, with paired 95%
interval [-0.33, 2.98]. It retains the energy-related acceptance bias.
Trimming late trees does not resolve the physics problem. The two genuine
score peaks also persist within Both-Both; category mixing alone cannot
explain them. Audit receipts and four PDFs are in `bdt_overtraining/`, and
the category-resolved BDT PDFs are in `bdt_categories/` under the new study.

The neural study preserves the existing labels, raw 50-input inference
contract and physical-event split; geometry variants transform or subset
these inputs internally. Its four validation-only rounds contain 22 initial,
20 stronger-penalty, 16 relative-geometry and 12 scheduled-penalty fits.
Signals are excluded from fitting, checkpoint choice and model selection.
The representative DNN has two hidden layers (64 and 32 nodes), a fixed
signed-log transform, imputation/scaling fitted on training only, dropout
and weight decay. Its loss
penalizes joint reconstructed mass/z/r dependence within labels, additional
genuine dependence within each supported SM process, and genuine score
dependence on vertex category. Process/category information is used for
training and validation/checkpointing, never passed to inference, which
receives reconstructed inputs only. These targets do not directly
measure lifetime or establish signal independence.

Balanced training draws use replacement and correct inverse-probability
weights. Best checkpoints are restored exactly. Introducing penalties over
20 optimization epochs improves some seed pairs, while checkpoint choice
uses the full intended validation loss from the first epoch. Stronger losses
can still collapse classification, so neither training software correctness
nor a low dependence statistic establishes useful separation. Relative
direction/DCA/fit-residual features did not recover the lost rejection in
these trials. They are straight-line proxies, not calibrated uncertainties.

The authoritative global lock is `comparison/selection_locked.json`, sealed
before any neural development-test score generation. It selects three roles
over all 70 validation fits and then evaluates those models and their paired
seeds only. The original expanded SM test was inspected during prior BDT
work; these are development results, not a new independent background test.
All cuts retain exactly 2,516/2,795 validation genuine pairs and are frozen.

| Expanded SM development test | Genuine retained | Fake rejected | Near-Both genuine retained | Lowest DY mass bin retained |
| --- | ---: | ---: | ---: | ---: |
| Current uniform BDT | 90.73% | 49.10% | 66.03% | 72.79% |
| Strong BDT | 90.61% | 70.80% | 64.38% | 47.79% |
| Flatter BDT | 91.52% | 40.62% | 78.36% | 85.29% |
| Rejection-selected DNN | 90.35% | 72.10% | 64.93% | 54.41% |
| Flatness-objective DNN | 90.88% | 45.19% | 84.38% | 84.56% |
| Full-envelope validation-selected DNN | 90.95% | 40.95% | 88.77% | 86.76% |

Names of the three selected neural models are respectively
`dnn_signedlog_class_l0_s71`, `dnn_r2_l3_c1_e0_h64_s72` and
`dnn_r4_l3_c3_e0_h64_s71`. The last model improves Near-Both by 38/365 pairs
relative to the flatter BDT: +10.41 percentage points, paired event-bootstrap
95% interval [6.03, 14.79]. Its pooled fake-rejection difference is only
+0.33 points, interval [-4.08, 4.57]. There is no demonstrated pooled-power
gain. Both-Both genuine retention is 90.97%, Both-Far 96.62%, Near-Far
94.12% (34 genuine pairs), and Other 81.13% (53 genuine pairs).

The last model reduces the worst supported one-dimensional acceptance span
to 6.59 points and maximum class-conditional joint distance correlation to
0.104. However, its test category span is 15.49 points, just beyond the
predeclared 15-point validation limit. Its paired seed gives 48.45% rejection
but retains only 73.58% of Other genuine pairs. A supported J/psi joint
mass/z/r cell retains 20/28 pairs (71.43%, Wilson interval 52.94--84.75%);
not every cell is supported. Therefore the
model is promising for category preservation, not certified acceptance-safe.
No test or signal result triggers further model selection in this study.

QCD contributes 89 fake and 14 unknown test pairs, with zero genuine pairs.
At the frozen diagnostic cuts the rejection-selected DNN rejects 91.01% of
QCD fakes; the full-envelope candidate rejects 44.94%. No genuine-QCD
efficiency, ROC, signal/background rate or normalization closure is claimed.
The empirical fractions above are unweighted; exposure-weighted comparisons
and finite-support flags are shown separately in the PDF diagnostics.

The new signal-transfer cache contains six complete production points,
1,200 unique GEN parents, 960 zero-pair parents, and 240 reconstructed pairs:
173 genuine signal-associated, 11 accidental and 56 unknown. It has no event
overlap with the original cache. The six-point inference is performed only
after the global lock, with all validation cuts unchanged. The unfinished
nine-point production cannot establish full mass/lifetime uniformity.
At unchanged cuts the full-envelope DNN keeps 160/173 genuine signal pairs,
compared with 147/173 for the flatter BDT and 173/173 for the strong BDT.
For the 50 GeV prompt point it keeps 26/29 versus the flatter BDT's 17/29.
Observed point retention spans 87.88--96.97%, compared with 58.62--100% for
the flatter BDT. Each point has only 24--33 genuine pairs; uncertainties
are wide. This is encouraging relative to the flatter BDT, not a general
advantage over the strong BDT. Signal results, full denominators, source
hashes and plots are linked by the study README. Pair retention is conditional
on existing reconstruction; it is not efficiency from the full GEN denominator.

All 132 standalone classifier tests and 18 targeted TEA tests pass, including
checkpoint restoration,
geometry contracts, scheduled-loss parity, provenance mutation guards and
existing ROOT augmentation tests. No shared CMSSW rebuild, deployed model
replacement, new ntuple selection or electronics change was made.

### 2026-10-12: stronger discrimination and the acceptance tradeoff

The user requested harder work after the broad process-score overlay. The
follow-up separates genuine, accidental and unknown truth classes explicitly,
then freezes each candidate protocol before fitting. It completes 25 new fits:
eight capacity/ordinary-boosting candidates, eight constant/adaptive uniformity
variants, five quality-only/quality-plus-angle capacity controls, and four
fixed recipes trained with the expanded available SM cache. Original training
and validation partitions alone select the first three phases. Previously
inspected test populations are consumed development evidence throughout.

The first portable winner, `uniform_d3_s0_n400`, rejects 71.24% of validation
accidentals versus 52.29% for the benchmark, both retaining exactly 379/420
genuine pairs (90.24%). Its AUC is 0.903 versus 0.798. On the larger original
held-out development population, the corresponding rejection is 65.76% versus
49.31%, with genuine retention 88.90% versus 89.01%. These cuts are the
validation empirical 90% cuts, distinct from the immutable published benchmark
cut 0.4403188643890805. Quality-only capacity controls do not improve the
acceptance-aware selection objective.

Adaptive local-efficiency reweighting tests a stabilized single-operating-point
variant inspired by [uBoost equations 2.3--2.4](https://arxiv.org/pdf/1305.7248).
It is not the paper's full ensemble. The better-flatness adaptive candidate
reduces the validation worst supported marginal acceptance span from 26.1 to
13.3 percentage points, while accidental rejection falls from 52.3% to 38.6%.
The historical constant-update defaults and scores are bitwise unchanged.

Expanded training uses exactly 20,173 SM retained pairs from the frozen cache;
39 signal pairs are excluded from fitting and selection. Deterministic
physical-event split seed 71 preserves every original pair's partition and
integer identity. Training has 9,826 known labels (8,046 genuine, 1,780
accidental), validation 3,400 (2,795/605), and development test 3,255
(2,642/613). Unknown labels number 2,245/714/733 in those partitions and are
excluded from supervised fitting. The expanded test excludes all newly used
training/validation events; 13,395 rows previously called held out are now
explicitly excluded. Signals remain inference/transfer checks only.

| Expanded development test | Genuine retained | Accidental rejected | AUC |
| --- | ---: | ---: | ---: |
| Frozen benchmark, expanded-validation diagnostic cut | 90.73% | 49.10% | 0.8294 |
| New `expanded_d3_s0_n400` | 90.61% | 70.80% | 0.9029 |
| New `expanded_d3_s3_n400` | 91.52% | 40.62% | 0.7548 |

All three diagnostic cuts retain exactly 2,516/2,795 validation genuine pairs
(90.018%). Thresholds are respectively 0.43534065598545313,
0.3529367670215394 and 0.2941898337347876, with inclusive comparison. None
changes the benchmark ntuple score or published threshold. The previous
first-phase winner remains the overall expanded-validation rejection winner
(66.28% versus 65.95% for the best newly trained recipe). The figures compare
the benchmark and the two newly fitted validation-selected candidates.

The stronger rejection is not acceptance-safe. Its pooled z-third retention
is 74.15/98.00/98.27%, versus benchmark 83.33/95.99/92.03%. Its lowest DY mass
third retains only about half of genuine pairs. The flatter candidate gives
pooled mass 92.17/92.05/90.38%, z 89.98/92.30/92.15% and radius
92.55/92.04/89.94%, but rejects fewer accidentals. Process dependence and
unsupported weighted populations remain; no new model passes the full
benchmark acceptance envelope. This is useful separation improvement with a
measured tradeoff, not a final physics selection.

The existing sigmoid score uses normalized tree votes and is uncalibrated;
its mathematical range is approximately 0.1192--0.8808. All 1,010,600 cached
predictor values are finite, so missing-value imputation is inactive. Score
stretching or endpoint peaks are not used as improvement criteria. The
feature audit finds lost-hit counts constant and valid-hit/muon-hit counts
identical. Independent min/max ordering loses association between per-leg
pt errors and momenta; exact relative per-leg errors require the original
associated Nano fields. Tested symmetric relative uncertainty/quality scales
do not individually outperform the existing fit-quality observables and are
not calibrated significances. No feature, label or reconstruction convention
was changed in this round.

At the already locked expanded-validation cuts, the original nine-point
signal cache retains 26/26 signal-associated genuine pairs for both the
benchmark and stronger rejection candidate, and 24/26 for the flatter
candidate. The latter loses one pair each at 15 GeV prompt and 30 GeV prompt.
Each point has only 1--5 genuine signal pairs, so this does not establish
mass/lifetime efficiency. The 39-pair population also contains one genuine
non-signal pair, 11 unknowns and one accidental; all remain explicit rather
than changing the signal denominator. The separate PDF/receipt is under
`signal_transfer/` and the project's `signal_transfer_genuine_retention.pdf`.

Read-only scheduler/publication checks verify all 60 completed QCD V11 pilot
jobs, 1,200 events, Nano hashes, identities and native/sampling weights.
Before BDT label/score inspection, 1,085 parents were reserved after excluding
115 overlaps with the old cache. They contain only nine retained pairs and
had prior non-BDT physics inspection. Four models and expanded-validation
cuts are locked before unblinding; this is a small independent-parent BDT
pipeline check, not a complete or pristine SM performance test. It contains
zero genuine, seven accidental and two unknown pairs. Accidental passage is
4/7 for the current model, 2/7 for the first stronger winner, 3/7 for the new
rejection candidate and 5/7 for the new flatter candidate. Genuine efficiency
and AUC are unavailable. The test preserves all 1,076 zero-pair parents and
has only about six effective accidental observations. Its exact receipt is
under `data_readiness/` in the study directory; no fitting follows unblinding.

The active larger signal stage is replacement DAG 12903259, with retained
workers from intentionally held original manager 12902110. Its current
ownership/collection authority is
[`validation/bsm_repartition_20261012/HANDOFF.md`](../../../../validation/bsm_repartition_20261012/HANDOFF.md).
Do not release the original manager or combine duplicate parents. This study
does not alter either campaign or the held QCD production supervisors.
The 14:09 snapshot has 1,770/1,800 successful detector events and five completed
points; three failed Nano-processing chunks require qualified recovery before
full collection. The mechanism finding and next diagnostic step are recorded
only in [SHIFT_RECONSTRUCTION.md](../../../../SHIFT_RECONSTRUCTION.md).

Protocols, source snapshots, selection locks, native/portable models, scores
and reports are in workspace
`validation/bdt_discrimination_improvement_20261012/`. Four expanded exports
reproduce native scores and threshold decisions bitwise on all 20,212 cached
pairs. The 103 classifier tests, 52 analysis tests and eight reserved-test
contracts pass; ten existing mixing tests are skipped in this environment.
Their logs are in the same directory. Reviewed PDF-only plots
are directly available through
[`tea_shift_cmssw/plots/BDT_improvement_20261012/README.md`](../../../../tea_shift_cmssw/plots/BDT_improvement_20261012/README.md).

Next: preserve these two ends of the rejection/acceptance tradeoff, investigate
reconstructed quality information with correct leg associations and
covariances, and obtain sufficient independent labelled SM support before
calling any replacement validated. Collect the new signal stage only through
its current handoff's publication/identity gates; do not tune on it or count
unknown labels as fake.

### 2026-10-12: fixed-model acceptance study and independent signal expansion

The one-hour follow-up study compares the saved quality-only BDT with the
uniformity-strength-3 benchmark on exactly the same physical events. The
selective scan reproduces every first-look identity, label, stored score,
mass/position nuisance and evaluation weight: 510 files, 187,394 events and
20,212 retained candidates. Audit-only truth positions and association details
are stored separately from the reconstructed predictors. The quality-only
bundle reproduces all 3,494 original saved predictions exactly. Neither model
is fitted or retuned; their original validation thresholds remain 0.3562883480311365
and 0.4403188643890805 respectively.

The held-out SM comparison contains 11,832 genuine, 2,383 accidental and 3,168
unknown pairs. It excludes all original training/validation identities and uses
200 physical-event bootstrap replicates, including paired model differences.
This is the previously inspected diagnostic population, not a fresh final test.

| Held-out empirical result | Uniformity BDT | Quality-only BDT |
| --- | ---: | ---: |
| Genuine retention | 88.83% | 86.65% |
| All-accidental rejection | 49.85% | 53.04% |
| Genuine SM radius thirds | 93.11%, 91.32%, 81.62% | 97.44%, 91.49%, 70.05% |
| Genuine DY mass thirds | 67.96%, 86.94%, 96.10% | 67.79%, 84.65%, 94.32% |

Uniformity gains 2.18 percentage points of genuine retention (paired 95%
interval 1.51--2.93), while retaining 3.19 points more accidentals
(0.88--5.75). The simpler model does not solve mass dependence and has stronger
radius dependence. Keep the uniformity model as the diagnostic benchmark;
neither cut passes the physics acceptance gates. With the declared MC mixture,
genuine DY retention is 63.56% for uniformity and 65.57% for quality-only,
demonstrating that the selected unweighted mixture is insufficient for a
physics-efficiency claim. Weighted accidentals have only 12.84 effective
independent events and fail the minimum-support requirement.

The dominant genuine QCD pair is traced to two distinct hit-associated muons
from an omega decay, with identical production positions. Original Nano,
generator ancestry, source descriptor, sampling ledger and normalization
records agree. Its weight arithmetic and 96.50% contribution to full-sample
genuine-SM weight are correct; it belongs to the original validation split
and is absent from held-out performance. The source still has provisional
physics/normalization gates: checking arithmetic does not close them. A
common-origin pair label does not certify the fitted mass or vertex position;
the new fit-precision finding is recorded only in
[SHIFT_RECONSTRUCTION.md](../../../../SHIFT_RECONSTRUCTION.md).

The full retained population has 3,703 unknown labels (18.32%): 17.95% DY,
18.46% J/psi, 18.71% QCD and 28.21% signal. Unknowns become more frequent at
larger radius and with poorer detector coverage. Most unavailable legs fail
the producer's hit-purity gate; relaxing downstream conventions cannot restore
discarded identities. Tighter purity 0.90 removes 9,378 labels, while tested
common-position tolerances change no known class. Angular matching remains a
diagnostic and does not replace hit-based labels. A separate identification
limits report reuses the exhaustively tested unknown-assignment calculation;
its bounds are conditional on correct known labels, not physics uncertainties
or efficiency corrections.
For the held-out SM population, uniformity's count-only genuine-retention
extrema are 84.07--90.78% and accidental-rejection extrema 24.33--60.85%.
These marginal endpoints need not share one unknown-label assignment and
must not be used as a paired model comparison. `paired_support.json` separately
flags all 180 empirical/weighted model differences, preserving the original
report; weighted interval availability is distinct from sufficient support.

Evidence is in workspace `validation/bdt_acceptance_study_20261012/`:
`comparison/results/comparison.json`, acceptance CSVs, paired differences,
joint/topology maps and before/after shapes; the inspected variants are in
`comparison/plots_reviewed/`. `comparison/validation.json` checks population,
partition, weight/flow and source-snapshot closure. `label_audit/` contains the
QCD source trace, match-completeness and convention audits, and full/held-out
fit-residual diagnostics. All original first-look outputs are preserved.

Signal DAG `12902110` was submitted at 10:18:33 UTC after nine new GEN audits,
resolved-config/two-event detector checks, empty/nonempty score/histogram probes
and seven publication/union contract tests. It contains 1,800 fresh GEN parents
across nine points, 90 twenty-event detector chunks and nine automatic complete-
point merging/scoring/histogram finalizers. Runs 25101201--25101209 are disjoint
from the original grid. All events and native weights remain; no BDT cut is
applied. Signals are validation-only, with roughly 10--50 genuine pairs per
point forecast from the small pilot, not the eventual 400-pair precision target.
At 10:20 UTC all 20 initial detector workers were running without holds or
failures. See `production/HANDOFF.md` and `submission_receipt.json` for current
scheduler reconciliation and post-completion collection. The existing QCD v11
pilot and bulk hold remain unchanged; overlapping V10/V11 parents prevent naive
combination or an independence claim. No shared release was rebuilt.

Fourteen focused comparison/statistics tests, seven label tests, seven
production contract tests, five identification-bound tests (including 1,555
exhaustive patterns) and nine post-completion reconciliation tests pass.
Independent production review verifies all 204 sealed files, 21 runtime
fingerprints, exact chunk unions and current shared/private dependencies.
Worker checks do not enforce every shared Python/private geometry digest;
current agreement is verified separately and the standing no-edit rule remains
essential. Before collecting new inputs, use the archived reconciliation tool
in `label_audit/postcompletion_tools/reconcile_signal.py`: it checks current
dependencies, all EOS publications and point receipts, then permits the sealed
collector with `--collect`. Its validation receipt gives the exact command.
The current preflight stops as incomplete without writing manifests. No live
or sealed file was modified to add this gate.

Next: inspect this independently generated signal stage, reconcile exact
publication receipts/hashes before collecting its evaluator inputs, and reserve
fresh SM final-test parents before any optimization of joint mass/position
acceptance. Preserve unknown regions and keep the score diagnostic throughout.

### 2026-10-12: completed score recovery and full signal-grid extension

The October 12 review found the original grouped production finished, but the
expired continuation had never submitted its repaired deployment. EOS marker
and receipt checks verified 494/503 publications. The old cluster's event log
records 97 successful batches and four failures, followed by removal of its
remaining jobs and factory on October 10.

After immutable-plan/dependency checks and scheduler reconciliation, recovery
`12901333` ran all 101 sealed v3 batches. All exited zero, independently
confirmed in Condor history. They verified/reused 494 outputs and completed
the nine missing groups (17--19, 128--129, 196--199). All 503 publications
preserve 187,254 events and 20,184 candidate pairs. The locked model, labels,
original branches, pair ordering, weights and no-cut policy are unchanged.

All nine physical signal points now have completed detector inputs. A separate
seven-point extension adds 140 events and 28 retained candidates, with 18
common-vertex genuine labels and 10 unknown labels. The upstream direct-parent
diagnostic identifies 17 authentic pairs in these seven files; that narrower
definition is kept distinct. Frozen scoring and histogram tools preserve all
original branches and streamer schemas, reproduce every score exactly, and
apply native `genWeight` once. Evidence and immutable manifests are under
`validation/bdt_signal_extension_20261012/`. One sandbox histogram attempt
timed out; its evidence is retained and the unchanged native rerun passed.

The combined first-look manifest covers 510 files, 187,394 events and 20,212
retained candidates across all 19 SM bins and all nine signal points. The SM
population remains the original 9,911-job snapshot, excluding 149 jobs from
the detector target; it is not a completed production-yield measurement.
The evaluator excludes final-model training/validation identities, uses the
unchanged threshold, and reports empirical and correctly weighted companions.
Results and figures are under `validation/bdt_first_look_20261012/`.

The full evaluation completed with 200 event-group bootstrap replicas.
It excludes 2,790 training/validation pairs; the remaining SM population
contains 11,832 genuine and 2,383 accidental pairs, including 107 accidental
QCD pairs. Empirical genuine-versus-all-accidental AUC is 0.803116
(95% interval 0.791711--0.812122); versus accidental QCD it is 0.876269
(0.842897--0.903102). At the frozen score threshold 0.4403188643890805:

| Held-out category | Passing / total | Retention |
| --- | ---: | ---: |
| Genuine J/psi | 6,813 / 7,745 | 87.97% |
| Genuine DY | 3,697 / 4,087 | 90.46% |
| All genuine SM | 10,510 / 11,832 | 88.83% |
| Accidental QCD | 37 / 107 | 34.58% |
| All accidental SM | 1,195 / 2,383 | 50.15% |

These empirical results describe the selected reconstructed sample. The
physics-mixture-weighted genuine-versus-QCD AUC is 0.895516
(0.849311--0.941967), but its QCD side has only 12.819 effective independent
events, below the required 20: `supported` remains false. This statistical
weight measure is distinct from the MC-to-data normalization N_eff.
The all-retained weighted genuine-SM score distribution is dominated by
one QCD genuine pair (96.5% of the class weight). It is excluded by the
training/validation veto; held-out genuine QCD performance is unmeasured.
The annotated label-comparison plots state this dominance explicitly.

All nine signal points are populated but have only 1--5 genuine direct-parent
pairs each: m15 prompt retains 4/5 and the other points retain 21/21 in total.
The aggregate 25/26 is descriptive only. The broad finite-count intervals
and unsupported grid cells do not establish mass/lifetime independence.
The completed report, input/source hashes, evaluation rows, CSV tables and
original figures are preserved in `performance/`; direct class overlays
are in `label_separation_annotated/` and all-sample shapes in `plots_all/`.
`performance_annotated/` rerenders that same frozen report with explicit
insufficient-support annotations and readable ROC labels. It preserves the
original figures and records a separate plotter/rendering receipt; no scores,
metrics, weights or evaluation populations are changed.

The histogram-merger/plot-preparer interface now accepts the production
`rows` manifest, keeps EOS path provenance, and labels combined SM processes
correctly. Full-grid and SM-only plot scopes have readable legends. Direct
genuine/accidental/unknown overlays use all retained pairs and are explicitly
separate from held-out performance. The 55 workflow and 12 ROOT plot-preparer
tests passed, together with actual signal scoring/histogram pilots and a signed
histogram display check. No shared CMSSW or TEA executable was rebuilt.

### 2026-10-09: merged-group first test and performance terminology

The user clarified that N_eff is the effective MC normalization fitted to
data. Classifier AUC and fixed-cut acceptance are separate quantities. The
new evaluation code and plots use those names. AUC constancy alone does not
establish uniform efficiency at a fixed cut, so both quantities are reported.
Historical frozen metadata flags are retained for provenance, not interpreted
as measurements of MC-to-data normalization.

The first test now uses the existing small merged groups, not whole-bin event
trees. Frozen source inventory is
`validation/dimuon_merged_score_test_20261009/source_inventory.json`: 501 SM
group files, 9,911 original Nano jobs, 187,214 events, 12,019,841,124 bytes and
all 19 SM bins. Each group contains up to 20 original files; every group keeps
its own scored output and histogram output. The inventory pins EOS identity,
Adler checksum, size and original receipt event counts. Two validated prompt
signal files add 40 events and 11 candidate pairs. Seven signal grid points
have no completed input yet and are explicitly unmeasured. Full V10 remains
incomplete by 149 detector jobs; original full GEN exposure remains the
normalization authority.

Added standalone BDT histogram and plotting tools without modifying the
live/shared TEA histogram configuration or C++ implementation. They retain
all candidates, require a Double score even for zero-pair files, and use
sampling weights once. Genuine/accidental/unknown score spectra and score
versus reconstructed mass/z/r are booked for all seven existing topologies.
SM bins are combined with sigma/full-parent-GEN-sumw before process shape
normalization. Each signal point stays separate and carries no rate claim.
Fifteen histogram/plot tests and real J/psi, DY and QCD pilots passed, along
with actual TEA PDF/PNG rendering.

New performance evaluation validates the frozen model and saved training
partitions, excludes fitted/validation events, reports empirical and
physics-mixture weighted ROC separately, and resamples whole events. It uses
proper c-tau from signal contracts, not reconstructed radius or estimated
laboratory mean flight. Unknown labels are excluded with counts retained;
empty/unsupported cells are not assigned a fake efficiency. Eleven tests and
real ROOT canaries passed. The first four-group canary has only three QCD
negatives, so its AUC estimates are machinery checks, not precision results.
Its 15 GeV prompt signal retains four of five genuine pairs at the locked cut.

The merged scorer exposed ROOT omission of an unused original STL descriptor
(`vector<edm::ProcessHistory>`, version 6). Restoring the original descriptor
kept all 32 source schemas, all 1,898 event branches, original key payloads
and exact recomputed scores. A second source layout has no writer-readable
directory TKey at fBEGIN, already present in the input before scoring.
Metadata-only repair now bootstraps the file header/free-space/streamer records
without parsing that obsolete directory key. Both actual failed-file probes
and six synthetic streamer tests plus four original scorer tests passed.
Existing schema changes still fail before publication; no guard was weakened.

The combined worker uses the frozen standard TEA runtime/application plus
the new score histograms. It verifies the legacy ordered Float32 cutflow sum
separately from exact Double BDT weight sums. Twenty-six worker tests passed,
including actual ROOT weight branches and native C++ arithmetic. The local
four-group canary completed 1,209 events and 58 pairs, with every published
object read back before its completion marker.

Initial grouped production was submitted as `12879497` from frozen
`deployment_v2`, 101 batches of five groups, capped at 20 materialized jobs.
After the user's small-group clarification, the earlier individual-file
campaign `12879319` was reversibly held: all 20 materialized jobs stopped;
completed outputs and logs remain. The small-group campaign never merges
whole-bin ntuples; only final histogram outputs will be merged by bin.

The two unsupported-directory cases failed in batches 2 and 3 before
publication. Good completed groups remain intact. `deployment_v3` freezes the
qualified metadata-only repair and proves identical 503 input/output rows,
101 batches, model, predictors, labels, weights, standard app and histogram
definitions before accepting old completed publications. Only writer/resume
code can differ. A bounded recovery canary is processing the five groups in
the first failed batch; no live/frozen configuration was edited.

### 2026-10-09: wrap-up and continuation checkpoint

The five-group recovery canary completed all files 10--14: 2,000 events,
65 pairs and no pending input. Standard histogram counts and legacy Float32
cutflow arithmetic, Double BDT histograms, original source contents and all
serialized schemas passed. `deployment_v3/recovery_canary_batch_status.json`
is the immutable gate for continuation. The repaired scorer uses only ROOT
header/free-space/streamer records; it never edits an event value or changes
the model, threshold or histogram selections.

Leave the valid grouped production running while it drains. The bounded
continuation watcher waits for cluster `12879497` to have no active/idle/
transferring/suspended jobs, then holds only that obsolete cluster and submits
the qualified v3 recovery once. It preserves all completed publications and
logs. Query failures are unknown, and ambiguous submission results stop
without a duplicate. Watcher state/log/submission receipts are under
`validation/dimuon_merged_score_test_20261009/continuation/`. The older
individual-file cluster `12879319` stays held; it must not be released for
this merged-group test. No jobs were cancelled and no data/logs were deleted.

The watcher started on `lxplus988` at 14:48:46 UTC on 2026-10-09, PID
`1037514`. Its first verified state is `waiting_for_old_jobs`, with nine
running, one idle and four held jobs in cluster `12879497`. It will wait up
to eight hours, then stop if the old jobs have not drained. The latest
validated publication snapshot is 240/503 groups at 14:41:35 UTC; this is a
checkpoint count, not a final completion count. `continuation/launch.json`
records the controller hash and PID; inspect `state.json` and any submission
receipt before restarting it. The watcher does not run the final histogram
merge or performance analysis. The detailed evaluation handoff is
`validation/dimuon_merged_score_test_20261009/PERFORMANCE_HANDOFF.md`.

On restart:

1. Read this checkpoint and workflow `docs/dimuon_classifier_scores.md`.
   Inspect continuation state/submission receipts before any submit/release.
   Reconcile queue/history, all 101 v3 returned statuses and exact EOS markers.
   Old publication hashes are accepted only through the frozen compatibility
   proof; every recovery batch status must have v3 plan SHA
   `9fb387d477c10f5c3ba99205a866cd4fcb475b65deca0d741f33c86e22fad6b5`.
2. If every one of the 503 groups is verified, run
   `merge_merged_dimuon_bdt_histograms.py` with that plan and a new work
   directory. It validates all source histograms and merges only histograms
   into 19 SM-bin and two signal-point outputs. Preserve missing/corrupt
   inputs as explicit failures; never silently skip a group.
3. Prepare process shape inputs using `prepare_shift_bdt_plots.py` and the
   merger's `bin_histogram_inputs.json`. Its sigma/full-parent-GEN exposure
   mixture preserves the original normalization ledger. Use
   `SHIFT_BDT_PLOT_INPUTS` with the dedicated `shift_bdt_plotter_config.py`.
   Tell the user the exact ready manifest and plot command only after merging
   and checks pass. Standard plus BDT histograms are present in the merged
   outputs; original small scored ntuples remain split by group.
4. Refresh `performance_inputs.json` using the preserved manifest-preparation
   tool and v3/compatible-old markers. Run `bdt_performance.py` with the frozen
   portable model and phase2 seed-71/72 training artifacts. Produce empirical
   and MC-mixture weighted ROC, fixed-cut efficiency and AUC maps versus
   proper c-tau/mass, reconstructed displacement maps and counts/intervals.
   Keep train/validation events excluded; weights and truth enter evaluation
   only. Empty/pending/unsupported cells remain unmeasured, with no claim of
   uniformity from a tiny or weight-dominated sample.
5. Re-query the remaining seven signal detector points. Add only completed,
   checksum-validated points via a new frozen inventory. Then increase signal
   statistics enough to compare common-threshold acceptance and AUC across
   mass/lifetime. Keep this classifier work separate from N_eff fitted to
   collision data; no collision-data access or N_eff fit is part of this test.

No final readiness announcement has been made: group production/recovery,
21 histogram merges and the full ROC/maps remain to be completed. Code is
committed on the required branches; runtime datasets, models, frozen plans,
ROOT files, receipts and logs remain outside Git. Shared CMSSW and live
campaign configurations were not rebuilt or changed.

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

### 2026-10-09: descriptive log name and scored-Nano preparation

Renamed this task log from the generic `STUDY.md` to
`DIMUON_VERTEX_CLASSIFIER.md` at the user's request. Updated active analysis
and BSM references; historical frozen receipts are unchanged.

The locked benchmark is technically ready to add a score column. It retains
357/396 genuine held-out SM pairs (90.2%), 113/122 genuine DY (92.6%) and
244/274 genuine J/psi (89.1%), and rejects 56/82 accidental QCD pairs (68.3%).
These are conditional pair efficiencies after the existing reconstruction and
pair selection. Genuine radius-bin efficiency still spans 81.7-94.6%, z-bin
efficiency 83.5-93.3%, with inadequate joint/per-process support. A final cut
with validated mass/lifetime-independent classifier performance is therefore
not ready. These metrics do not measure the MC-to-data normalization N_eff.

The v10 workers publish only Nano, their completed receipt and evidence.
No retained SM v10 AOD archive was found; four local BSM AOD pilots exist.
The existing Nano already supplies every predictor, so direct augmentation
preserves the exact sample and its sampling weights. Re-running plain Step 4
would also require recovering the sampling metadata by event identity.
An active grid worker uses the shared CMSSW release; no rebuild is performed.

Added an explicit JSON export of the 160-tree uniform BDT. `export_bdt.py`
and `portable_inference.py` preserve the frozen medians, Float32 tree-input
rounding, Float64 thresholds, ordered coefficients and sigmoid exactly.
All 3,494 comparison scores and historical threshold decisions agree bitwise
with the trusted pickle (maximum difference zero). Production imports neither
sklearn nor SciPy and does not load pickle models. Save scores as Float64:
casting to Float32 changes one decision at an observed threshold tie.
The model and parity receipt are under `artifacts/portable_bdt_v1/`.

`add_bdt_score.py` appends `ShiftDimuonVertex_bdtScore` as a Double array in a
distinct copy, with one value per original pair. No threshold, filter, pairing
or ordering change is applied. All original event values/types and non-Events
payloads are checked, along with ROOT type/checksum/streamer schemas (excluding
only runtime TObject `fBits`). Native ROOT fixtures exercise zero/multiple
pairs, large event IDs, nonfinite mass targets, weight preservation, failure
without filtered publication and no overwrite. Prediction no longer requires
the unused mass target; invalid geometry features receive frozen imputation.

The complete classifier suite passed 82 tests. Five actual-file pilots
preserved 82 events and ten pairs across v10 J/psi, DY, QCD, a QCD file with
zero pairs and an existing 50 GeV signal. Every original branch value/type,
key payload and all 32 original streamer schemas per file are preserved.
Receipts are in `artifacts/scored_nano_pilot_v2/`, with a separate dictionary
check. The earlier prototype is retained in `scored_nano_pilot_v1/`.

The workflow preparer passed 15 tests and froze
`validation/dimuon_score_production_20261009/plan.json`: 9,911 existing files
across all 19 SM bins, 199 batches of up to 50 files, plus a three-file pilot.
The 149 missing files relative to the 10,060 expected inventory are explicit.
The deployment pins the model, source, original inventory and LCG_108 runtime.
Its worker uses XRootD staging, verifies the original published source checksum,
and publishes new Nano, score receipt and scorer log, with readbacks before a
completion marker. Complete exact-plan outputs can be resumed; orphans and
unknown state stop without overwrite. Failed-file evidence is retained.

The three-file production pilot, Condor `12878715.0`, passed: exit code zero,
normal termination, 175 seconds on the worker, 60 events and one retained pair.
The returned status has no pending files. All three published Nano, receipt
and scorer-log checksums were independently verified; their source weights,
original content and exact recomputed scores agree. Two zero-pair files are
retained. Sandbox Condor DNS resolution failed; the verified native context
successfully submitted and queried the jobs.

After that check, started the complete 199-batch score production as Condor
`12879017`, with at most 20 materialized jobs. Initial status was five idle
jobs, no holds. The new EOS tree is the existing v10 campaign name with
`_bdt_v1` appended; no source ntuple is modified. Exact pilot publications
are reused by the same-plan campaign. The frozen plan remains immutable;
launch/status information is recorded separately in `pilot_submission.json`,
`production_submission.json`, returned batch statuses and EOS completion
markers under `validation/dimuon_score_production_20261009/`.
Instructions are in the workflow `docs/dimuon_classifier_scores.md`.

The bulk launch then encountered a concurrent EOS metadata-only storage
migration. Fifteen materialized batches became held after original source
Nano paths returned XRootD [3011], while source completion receipts were
still reachable. Eleven returned file records had passed before the failure,
including reuse of one pilot output. Together with the pilot, these represent
13 distinct scored files. This is an input-location interruption, not a ROOT
writer, score or memory failure; no original data was deleted by this task.

The storage migration has an explicit incomplete `path_map.json` guard.
Its audit reports no missing originals or metadata mismatches, and 52 new
files from this score campaign. A separate exact additions manifest is being
retained for that audit. Existing frozen jobs/configuration, successful files,
logs and failure archives remain unchanged. No release, retry or cancellation
was performed. The reusable preparer is being extended to freeze a completed
map and the canonical storage helper and to use the dated process-bin output
layout. Production recovery must wait for the completed migration audit;
the BDT remains technically ready for scores and unready for a physics cut.

The migration-aware preparer is now complete and committed as workflow
`255d359`, following initial scoring deployment commit `e31d187`. Its 19 tests
passed. Completed maps and the canonical storage helper are frozen by digest;
logical source/receipt identities stay unchanged, while transport paths use
the map. Outputs use a dated process-bin campaign with direct `nanoAOD` and
separate metadata tiers. The worker checks the map at startup and before each
file. An actual preparation against the workspace's incomplete guard exited
nonzero and created no deployment, recorded in `migration_gate_check.json`.
No recovery jobs were submitted. The portable scorer, model export and log
rename are committed as CMSSW `59c392b255e`.

The storage owner then completed the migration. Its final audit passes with
649,598 files, 447,355,565,782 bytes and 117,097 ROOT files on both sides, with
no missing/unexpected files, metadata mismatches, duplicated IDs or layout
errors. The completed map's audit digest and the audit's planned-map digest
match; all four mapping sections are identical to the audited plan.

Independent reconciliation in `concurrent_additions.json` verifies all 13
scored publications: 52 artifacts, 250 events, two pairs and 18,681,525 bytes.
Every relocated Nano, receipt and log SHA-256 and every completion-marker byte
matches its original marker. Sizes, Adler checksums and mtimes agree with the
migration snapshot. The historical unexpected-52 set is explained exactly;
the completed audit has zero unexpected files.

Prepared the fresh immutable
`validation/dimuon_score_production_20261009_v2/` deployment with 9,911 files,
199 batches and the completed 221 MiB map. Output campaign
`shift_detector_20261009_v10_bdt_v2` uses the canonical dated process-bin
layout. New pilot `12879278.0` was submitted at 12:27:32 UTC. Bulk recovery
waits for this pilot's complete output validation and measured peak memory;
the large map makes measured memory use a required runtime check.

The new pilot passed, normal exit zero at 12:38:18 UTC after 473 seconds.
It completed exactly three publications, 60 events and one pair, with no
pending input. Independent native XRootD readback checked all 12 published
objects. Nano/receipt/log hashes and exact completion-marker bytes match.
The actual score arrays are finite `double[]`/Float64, aligned with existing
vertex counts; the J/psi pair score is `0.5272693221696033`. Original-content
and exact re-scoring verification passed in the worker. The independent
receipt is `pilot_readback.json` in the v2 deployment.

The site's actual allocation/request was 3,000 MB despite the submit file's
2,000 MB request. Terminal RSS/image size was 1,500,000 KB, corresponding to
about 1,465 MB reported peak usage and 1,535 MB observed allocation margin.
PSS was unavailable. The 2 GB value describes the submission request, not a
measured hard limit. These runtime results allow bulk launch with the existing
sealed plan and resource submission. `selection_applied`, `physics_ready` and
`normalization_transfer_validated` remain false.

After the independent pilot checks, submitted the complete v2 production as
Condor `12879319` at 12:46:04 UTC: 199 logical batches covering 9,911 existing
v10 Nano files. At 12:47:36 UTC, 12 jobs were running and five idle, with no
holds. Materialization is capped at 20 jobs and idle jobs at five. The exact
pilot publications are reusable by their same-plan bulk batches. Submission
and initial queue receipts are `production_submission.json` and
`production_queue_initial.json` in the v2 deployment. This is a running
campaign, not evidence that all 9,911 outputs are complete. The earlier held
campaign `12879017` and its outputs/evidence remain untouched. The source
files and original sampling weights remain the normalization authority.

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
