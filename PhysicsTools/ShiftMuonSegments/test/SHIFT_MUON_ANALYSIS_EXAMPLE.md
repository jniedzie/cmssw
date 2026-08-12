# ShiftMuon analysis example

The unqualified `ShiftMuon_pt`, `eta`, `phi`, `p`, `px`, `py`, `pz`,
`vx`, `vy`, `vz`, `chi2`, `ndof`, and all other established branches are the
default unconstrained reconstruction. Only the target-restricted alternative
uses the `ShiftMuon_constrained*` prefix. Use it only when
`ShiftMuon_constrainedValid == 1`; `constrainedTargetChi2` records compatibility
with the configured production target.

There are two reconstruction-category fields:

- `quality`: exclusive category `0=cosmic`, `1=DSA`, `2=traversing`,
  `3=full-lever-arm traversing`;
- `sourceIndex`: index in the corresponding CMSSW source collection.

The checked-in example contains one real populated event and the complete
`ShiftMuon`, `ShiftDimuonVertex`, and `GenPart` schemas. Regenerate it from any
current NanoAOD with:

```bash
python3 create_shift_muon_analysis_example.py input.root shift_muon_analysis_example.root
```
