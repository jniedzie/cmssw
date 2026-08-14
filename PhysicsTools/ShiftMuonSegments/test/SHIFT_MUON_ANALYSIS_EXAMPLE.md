# ShiftMuon analysis example

The unqualified `ShiftMuon_pt`, `eta`, `phi`, `p`, `px`, `py`, `pz`,
`vx`, `vy`, `vz`, `chi2`, `ndof`, and all other established branches are the
default unconstrained reconstruction. Only the target-restricted alternative
uses the `ShiftMuon_constrained*` prefix. Use it only when
`ShiftMuon_constrainedValid == 1`; `constrainedTargetChi2` records compatibility
with the configured production target.

Measured topology and reconstruction provenance are separate:

- `topology`: `0=near-endcap-only`, `1=near-endcap-and-barrel`,
  `2=both-endcaps`, `3=far-endcap-only`, `4=unclassified`, relative to the
  inferred source side;
- `orientationValid`: whether that source-side orientation could be inferred;
- `recoAlgorithm`: `0=DSA`, `1=strict traversing`, `2=ordinary cosmic`;
- `sourceIndex`: index in the corresponding reconstruction collection;
- `duplicateSourceMask` and the three `duplicate*Count` fields describe which
  reconstruction algorithms were present before duplicate cleaning.

`quality` remains temporarily for compatibility. It is deprecated because it
mixes reconstruction provenance with the historical 500 cm `outerZ-innerZ`
split and must not be used for new categories.

The category can be compared with reconstruction-neutral detector geometry:

- `nDTHits`, `nCSCHits`, `nRPCHits`, `nGEMHits`, and `nME0Hits` count valid
  track measurements, while `detectorMask` records which systems occur;
- `nHitsPlusEndcap`, `nHitsBarrel`, and `nHitsMinusEndcap` retain the raw
  detector-side counts;
- `nHitsNearEndcap`, `nHitsFarEndcap`, and the corresponding `nStations*`
  fields orient those counts using `inferredSourceSide`;
- `entry*` and `exit*` describe the first and last recorded valid muon hit in
  the inferred flight direction. Region codes are `0=near endcap`, `1=barrel`,
  `2=far endcap`, and `-1=unknown`. These are measurement endpoints, not proof
  that a particle physically entered or stopped at those positions;
- `hitSpan` is the straight-line separation between those recorded endpoints;
  `hitDeltaZ` is their absolute z separation.

To print reconstruction-algorithm-versus-topology, detector, eta, momentum, and hit-span
summaries from a current NanoAOD, run:

```bash
python3 analyze_shift_muon_topology.py input.root
```

`ShiftDimuonVertex` follows the same fit naming convention. Its unqualified
kinematics and vertex quantities use the two unconstrained muons, while the
`constrained*` alternatives exist only for the same pair with both constrained
muon fits; `constrainedValid` must be required. Mixed constrained-unconstrained
pairs are not produced.

`topologyMin` and `topologyMax` give the compact unordered topology pair for
each dimuon. During the compatibility period, exactly one legacy unordered
quality-pair flag is also true:
`isCosmicCosmic`, `isCosmicDSA`, `isCosmicTraversing`,
`isCosmicDoubleTraversing`, `isDSADSA`, `isDSATraversing`,
`isDSADoubleTraversing`, `isTraversingTraversing`,
`isTraversingDoubleTraversing`, or
`isDoubleTraversingDoubleTraversing`. The order of `muonIdx1` and `muonIdx2`
does not affect these flags.

The checked-in example contains one real populated event and the complete
`ShiftMuon`, `ShiftDimuonVertex`, and `GenPart` schemas. Regenerate it from any
current NanoAOD with:

```bash
python3 create_shift_muon_analysis_example.py input.root shift_muon_analysis_example.root
```
